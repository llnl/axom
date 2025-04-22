// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file convert_sidre_protocol.cpp
 * \brief This file contains a utility to convert a Sidre datastore
 *        from the sidre_hdf5 protocol to another supported protocol.
 *
 * Users must supply a path to a sidre_hdf5 rootfile and base name for
 * the output datastores.  Optional command line arguments include
 * a '--protocol' option (the default is 'json')
 * and a '--strip' option to truncate the array data to at most N elements.
 * The strip option also prepends each array with its original size, the new
 * size and a filler entry of 0 for integer arrays or nan for floating point
 * arrays. E.g. if the array had 6 entries [1.01. 2.02, 3.03, 4.04, 5.05, 6.06]
 * and the user passed in --strip 3, the array would be converted to
 * [6, 3, nan, 1.01, 2.02, 3.03].
 *
 * \note The strip option is intended as a temporary solution to truncating
 * a dataset to allow easier debugging.  In the future, we intend to separate
 * the conversion and truncation/display functionality into separate utilities.
 */

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"
#include "axom/slam.hpp"

#include "mpi.h"
#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <limits>   // for numeric_limits<int>
#include <cstdlib>  // for atoi
#include <sstream>  // for stringstream

namespace sidre = axom::sidre;
namespace slam = axom::slam;
namespace slic = axom::slic;

using PosType = slam::DefaultPositionType;
using ElemType = sidre::IndexType;

using SzPol = slam::policies::RuntimeSize<PosType>;
using OffPol = slam::policies::ZeroOffset<PosType>;
using StrPol = slam::policies::RuntimeStride<PosType>;
using ViewSet = slam::OrderedSet<PosType, ElemType, SzPol, OffPol, StrPol>;

void setupLogging();
void teardownLogging();

/** Simple structure to hold the parsed command line arguments */
struct CommandLineArguments
{
  static const std::set<std::string> s_validProtocols;

  std::string m_inputName;
  std::string m_outputName;
  std::string m_protocol;
  int m_numStripElts;

  CommandLineArguments() : m_inputName(""), m_outputName(""), m_protocol("json"), m_numStripElts(-1)
  { }

  void parse(int argc, char** argv, axom::CLI::App& app);

  bool shouldStripData() const { return m_numStripElts >= 0; }

  /**  Returns the maximum allowed elements in a view of the output datastore */
  int maxEltsPerView() const
  {
    return shouldStripData() ? m_numStripElts : std::numeric_limits<int>::max();
  }
};

const std::set<std::string> CommandLineArguments::s_validProtocols({"json",
                                                                    "sidre_hdf5",
                                                                    "sidre_conduit_json",
                                                                    "sidre_json",
                                                                    "conduit_hdf5",
                                                                    "conduit_bin",
                                                                    "conduit_json"});

/** Terminates execution */
void quitProgram(int exitCode = 0)
{
  teardownLogging();
  MPI_Finalize();
  exit(exitCode);
}

/** Parse the command line arguments */
void CommandLineArguments::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i,--input", m_inputName, "Filename of input sidre-hdf5 datastore")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_option("-o,--output", m_outputName, "Filename of output datastore (without extension)")
    ->required();

  app.add_option("-p,--protocol", m_protocol, "Desired protocol for output datastore")
    ->capture_default_str()
    ->check(axom::CLI::IsMember {CommandLineArguments::s_validProtocols});

  app
    .add_option("-s,--strip",
                m_numStripElts,
                "If provided, output arrays will be stripped to first N entries")
    ->check(axom::CLI::PositiveNumber);

  bool verboseOutput = false;
  app.add_flag("-v,--verbose", verboseOutput, "Sets output to verbose")->capture_default_str();

  app.get_formatter()->column_width(35);

  // Could throw an exception
  app.parse(argc, argv);

  if(verboseOutput)
  {
    slic::setLoggingMsgLevel(slic::message::Debug);
  }
}

/**
 * \brief Allocate storage for external data of the input datastore
 *
 * Iterates recursively through the views and groups of the provided group to
 * find the external data views and allocates the required storage within the
 * extPtrs vector
 *
 * \param grp  The group to traverse
 * \param extPtrs [out] A vector to hold pointers to the allocated data
 *
 * \note Also initializes the data in each allocated array to zeros
 */
void allocateExternalData(sidre::Group* grp, std::vector<void*>& extPtrs)
{
  // for each view
  for(auto idx = grp->getFirstValidViewIndex(); sidre::indexIsValid(idx);
      idx = grp->getNextValidViewIndex(idx))
  {
    sidre::View* view = grp->getView(idx);
    if(view->isExternal())
    {
      SLIC_DEBUG("External view " << view->getPathName() << " has " << view->getNumElements()
                                  << " elements "
                                  << "(" << view->getTotalBytes() << " bytes).");

      const int idx = extPtrs.size();
      const int sz = view->getTotalBytes();
      extPtrs.push_back(new char[sz]);
      std::memset(extPtrs[idx], 0, sz);
      view->setExternalDataPtr(extPtrs[idx]);
    }
  }

  // for each group
  for(auto idx = grp->getFirstValidGroupIndex(); sidre::indexIsValid(idx);
      idx = grp->getNextValidGroupIndex(idx))
  {
    allocateExternalData(grp->getGroup(idx), extPtrs);
  }
}

/**
 * \brief Shifts the data to the right by three elements
 *
 * The new first value will be the size of the original array.
 * The next value will be the number of retained elements and the
 * third value will be 0 for integer data and Nan for float data.
 * This is followed by the values in the truncated original dataset.
 *
 * \note This function creates a copy of the data since there could be
 * several views in the original dataset pointing to the same memory
 *
 * \param view The array view on which we are operating
 * \param origSize The size of the original array
 */
template <typename sidre_type>
void modifyFinalValuesImpl(sidre::View* view, int origSize)
{
  sidre_type* arr = view->getData();
  int sz = view->getNumElements();

  // Uses a Slam set to help manage the indirection to the view data
  // Note: offset is zero since getData() already accounts for the offset
  ViewSet viewInds = ViewSet::SetBuilder().size(sz).stride(view->getStride());

#ifdef AXOM_DEBUG
  {
    axom::fmt::memory_buffer out;
    for(auto i : viewInds.positions())
    {
      axom::fmt::format_to(std::back_inserter(out),
                           "\n\ti: {0}; index: {1}; arr[{1}] = {2}",
                           i,
                           viewInds[i],
                           arr[viewInds[i]]);
    }
    SLIC_DEBUG("Before truncation" << axom::fmt::to_string(out));
  }
#endif

  // Create a new buffer for copied data
  auto* ds = view->getOwningGroup()->getDataStore();
  auto type = view->getTypeID();
  auto newSz = sz + 3;
  auto* buff = ds->createBuffer(type, newSz)->allocate();

  // Explicitly set the first two elements and copy elements over
  sidre_type* newArr = buff->getData();

  newArr[0] = static_cast<sidre_type>(origSize);
  newArr[1] = static_cast<sidre_type>(sz);
  newArr[2] = std::numeric_limits<sidre_type>::quiet_NaN();
  for(auto i : viewInds.positions())
  {
    newArr[i + 3] = arr[viewInds[i]];
  }

  // Update view's buffer to the new data
  view->detachBuffer();
  view->attachBuffer(type, newSz, buff);

#ifdef AXOM_DEBUG
  {
    axom::fmt::memory_buffer out;
    for(auto i : ViewSet(newSz).positions())
    {
      axom::fmt::format_to(std::back_inserter(out), "\n\ti: {0}; arr[{0}] = {1}", i, newArr[i]);
    }
    SLIC_DEBUG("After truncation" << axom::fmt::to_string(out));
  }
#endif
}

void modifyFinalValues(sidre::View* view, int origSize)
{
  SLIC_DEBUG("Truncating view " << view->getPathName());

  switch(view->getTypeID())
  {
  case sidre::INT8_ID:
    modifyFinalValuesImpl<std::int8_t>(view, origSize);
    break;
  case sidre::INT16_ID:
    modifyFinalValuesImpl<std::int16_t>(view, origSize);
    break;
  case sidre::INT32_ID:
    modifyFinalValuesImpl<std::int32_t>(view, origSize);
    break;
  case sidre::INT64_ID:
    modifyFinalValuesImpl<std::int64_t>(view, origSize);
    break;
  case sidre::UINT8_ID:
    modifyFinalValuesImpl<std::uint8_t>(view, origSize);
    break;
  case sidre::UINT16_ID:
    modifyFinalValuesImpl<std::uint16_t>(view, origSize);
    break;
  case sidre::UINT32_ID:
    modifyFinalValuesImpl<std::uint32_t>(view, origSize);
    break;
  case sidre::UINT64_ID:
    modifyFinalValuesImpl<std::uint64_t>(view, origSize);
    break;
  case sidre::FLOAT32_ID:
    modifyFinalValuesImpl<axom::float32>(view, origSize);
    break;
  case sidre::FLOAT64_ID:
    modifyFinalValuesImpl<axom::float64>(view, origSize);
    break;
  default:
    break;
  }
}

/**
 * \brief Recursively traverse views and groups in grp and truncate views to
 * have at most maxSize+2 elements.
 *
 * Within the truncated arrays, the first element will be the size of the
 * original array and the second will
 * be 0 for integers and nan for floating points.
 * This will be followed by (at most) the first maxSize elements of the original
 * array
 */
void truncateBulkData(sidre::Group* grp, int maxSize)
{
  // Add two to maxSize
  for(auto idx = grp->getFirstValidViewIndex(); sidre::indexIsValid(idx);
      idx = grp->getNextValidViewIndex(idx))
  {
    sidre::View* view = grp->getView(idx);
    bool isArray = view->hasBuffer() || view->isExternal();

    if(isArray)
    {
      const int numOrigElts = view->getNumElements();
      const int newSize = std::min(maxSize, numOrigElts);

      if(view->hasBuffer() && numOrigElts > newSize)
      {
        const int viewStride = view->getStride();
        const int viewOffset = view->getOffset();
        view->apply(newSize, viewOffset, viewStride);
      }
      // external
      else if(view->isExternal() && numOrigElts > newSize)
      {
        void* dataPtr = view->getVoidPtr();
        view->setExternalDataPtr(view->getTypeID(), newSize, dataPtr);
      }

      modifyFinalValues(view, numOrigElts);
    }
  }

  // for each group
  for(auto idx = grp->getFirstValidGroupIndex(); sidre::indexIsValid(idx);
      idx = grp->getNextValidGroupIndex(idx))
  {
    truncateBulkData(grp->getGroup(idx), maxSize);
  }
}

/** Sets up the logging using lumberjack */
void setupLogging()
{
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

#ifdef AXOM_USE_LUMBERJACK
  std::string rankStr = "[<RANK>]";
#else
  std::string rankStr = "";
#endif

  // Formatting for warning, errors and fatal message
  std::stringstream wefFmt;
  wefFmt << "\n***********************************\n"
         << rankStr << "[<LEVEL> in line <LINE> of file <FILE>]\n"
         << "MESSAGE=<MESSAGE>\n"
         << "***********************************\n";
  std::string wefFormatStr = wefFmt.str();

  // Simple formatting for debug and info messages
  std::string diFormatStr = rankStr + "[<LEVEL>]: <MESSAGE>\n";

  slic::LogStream* wefStream;
  slic::LogStream* diStream;

#ifdef AXOM_USE_LUMBERJACK
  const int ranksLimit = 16;
  wefStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, ranksLimit, wefFormatStr);
  diStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, ranksLimit, diFormatStr);
#else
  wefStream = new slic::GenericOutputStream(&std::cout, wefFormatStr);
  diStream = new slic::GenericOutputStream(&std::cout, diFormatStr);
#endif

  slic::addStreamToMsgLevel(wefStream, slic::message::Error);
  slic::addStreamToMsgLevel(wefStream, slic::message::Warning);
  slic::addStreamToMsgLevel(diStream, slic::message::Info);
  slic::addStreamToMsgLevel(diStream, slic::message::Debug);

  // the following is helpful for debugging
  // slic::debug::checksAreErrors = true;
}

/** Finalizes logging and flushes streams */
void teardownLogging() { slic::finalize(); }

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  setupLogging();

  // parse the command line arguments
  CommandLineArguments args;
  axom::CLI::App app {"Sidre protocol converter"};

  try
  {
    args.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(my_rank == 0)
    {
      retval = app.exit(e);
    }
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    quitProgram(retval);
  }

  // Load the original datastore
  SLIC_INFO("Loading datastore from " << args.m_inputName);
  sidre::DataStore ds;
  sidre::IOManager manager(MPI_COMM_WORLD);
  manager.read(ds.getRoot(), args.m_inputName);
  int num_files = manager.getNumFilesFromRoot(args.m_inputName);

  // Restore any external data pointers
  SLIC_INFO("Loading external data from datastore");
  std::vector<void*> externalDataPointers;
  allocateExternalData(ds.getRoot(), externalDataPointers);
  manager.loadExternalData(ds.getRoot(), args.m_inputName);

  slic::flushStreams();

  // Internal processing
  if(args.shouldStripData())
  {
    const int numElts = args.maxEltsPerView();
    SLIC_INFO("Truncating views to at most " << numElts << " elements.");

    truncateBulkData(ds.getRoot(), numElts);

    // Add a string view to the datastore to indicate that we modified the data
    std::stringstream sstr;
    sstr << "This datastore was created by axom's 'convert_sidre_protocol' "
         << "utility with option '--strip " << numElts << "'. "
         << "To simplify debugging, the bulk data in this datastore has been "
         << "truncated to have at most " << numElts << " original values "
         << "per array. Three values have been prepended to each array: "
         << "the size of the original array, the number of retained elements "
         << "and a zero/Nan.";

    ds.getRoot()->createViewString("Note", sstr.str());
  }

  // Write out datastore to the output file in the specified protocol
  SLIC_INFO("Writing out datastore in '"
            << args.m_protocol << "' protocol to file(s) with base name " << args.m_outputName);
  manager.write(ds.getRoot(), num_files, args.m_outputName, args.m_protocol);

  // Free up memory associated with external data
  for(std::vector<void*>::iterator it = externalDataPointers.begin();
      it != externalDataPointers.end();
      ++it)
  {
    delete[] static_cast<char*>(*it);
    *it = nullptr;
  }

  teardownLogging();
  MPI_Finalize();
  return 0;
}
