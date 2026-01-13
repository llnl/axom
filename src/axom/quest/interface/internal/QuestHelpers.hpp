// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_HELPERS_HPP_
#define QUEST_HELPERS_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/primal.hpp"
#include "axom/slic.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/io/STLReader.hpp"
#include "axom/quest/io/ProEReader.hpp"

#if defined(AXOM_USE_MFEM)
  #include <mfem.hpp>
#endif

// C/C++ includes
#include <string>

/*!
 * \file QuestHelpers.hpp
 *
 * \brief Helper methods that can be used across the different Quest queries.
 */
namespace axom
{
namespace quest
{
namespace internal
{
constexpr int READ_FAILED = -1;
constexpr int READ_SUCCESS = 0;

/*!
 * \brief Simple RAII-based utility class to update the slic logging level within a scope
 *
 * The original logging level will be restored when this instance goes out of scope
 */
class ScopedLogLevelChanger
{
public:
  ScopedLogLevelChanger(slic::message::Level newLevel)
  {
    if(slic::isInitialized())
    {
      m_previousLevel = slic::getLoggingMsgLevel();
      slic::setLoggingMsgLevel(newLevel);
    }
  }

  ~ScopedLogLevelChanger()
  {
    if(slic::isInitialized())
    {
      slic::setLoggingMsgLevel(m_previousLevel);
    }
  }

private:
  slic::message::Level m_previousLevel {slic::message::Level::Debug};
};

/// \name Mesh I/O methods
/// @{

#if defined(AXOM_USE_UMPIRE_SHARED_MEMORY)

/*!
 * \brief Reads in the surface mesh from the specified file into a shared
 *  memory buffer.
 *
 * \param [in] file the file consisting of the surface mesh
 * \param [in] global_comm handle to the global MPI communicator
 * \param [out] mesh_buffer pointer to the raw mesh buffer
 * \param [out] m pointer to the mesh object
 *
 * \return status set to READ_SUCCESS, or READ_FAILED on error.
 *
 * \note Each rank has a unique mint::Mesh object instance, however, the
 *  mint::Mesh object is constructed using external pointers that point into
 *  the supplied mesh_buffer, an on-node data-structure shared across all
 *  MPI ranks within the same compute node.
 *
 * \pre global_comm != MPI_COMM_NULL
 * \pre mesh_buffer == nullptr
 * \pre m == nullptr
 *
 * \post m != nullptr
 * \post m->isExternal() == true
 * \post mesh_buffer != nullptr
 */
int read_stl_mesh_shared(const std::string& file,
                         MPI_Comm global_comm,
                         unsigned char*& mesh_buffer,
                         mint::Mesh*& m);
#endif

/*!
 * \brief Reads in the surface mesh from the specified file.
 *
 * \param [in] file the file consisting of the surface
 * \param [out] m user-supplied pointer to point to the mesh object.
 * \param [in] comm the MPI communicator, only applicable when MPI is available.
 *
 * \note This method currently expects the surface mesh to be given in STL format.
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function.
 *
 * \return status set to zero on success, or to a non-zero value otherwise.
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::TRIANGLE
 *
 * \see STLReader
 * \see PSTLReader
 */
int read_stl_mesh(const std::string& file, mint::Mesh*& m, MPI_Comm comm = MPI_COMM_SELF);

#ifdef AXOM_USE_C2C
/*!
 * \brief Reads in the contour mesh from the specified file and linearize it.
 *
 * \param [in] file the file consisting of a C2C contour defined by one or more c2c::Piece
 * \param [in] uniform If true, the curves will be linearized uniformly, according to segmentsPerPiece.
 *                     Otherwise, the linearization will be non-uniform based on \a percentError.
 * \param [in] transform A 4x4 matrix that contains a transform to be applied to points.
 * \param [in] segmentsPerPiece number of segments to sample per contour Piece
 * \param [in] vertexWeldThreshold threshold for welding vertices of adjacent curves
 * \param [in] percentError An error tolerance (percent of lgneth) used in non-uniform curve linearization.
 * \param [out] m user-supplied pointer to point to the mesh object
 * \param [out] revolvedVolume An approximation of the revolved volume of the contour
 *                             or 0 if it could not be computed.
 * \param [in] comm the MPI communicator, only applicable when MPI is available
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function
 *
 * \return status set to zero on success, or to a non-zero value otherwise
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::SEGMENT
 * \post revolvedVolume > 0 if it could be computed.
 *
 * \see C2CReader
 * \see PC2CReader
 * \see LinearizeCurves
 */
int read_c2c_mesh(const std::string& file,
                  bool uniform,
                  const numerics::Matrix<double>& transform,
                  int segmentsPerPiece,
                  double vertexWeldThreshold,
                  double percentError,
                  mint::Mesh*& m,
                  double& revolvedVolume,
                  MPI_Comm comm = MPI_COMM_SELF);

#endif  // AXOM_USE_C2C

#ifdef AXOM_USE_MFEM
/*!
 * \brief Reads in the contour mesh from the specified file and linearize it.
 *
 * \param [in] file the file consisting of a contour defined by one or more bezier or NURBS zones.
 * \param [in] uniform If true, the curves will be linearized uniformly, according to segmentsPerPiece.
 *                     Otherwise, the linearization will be non-uniform based on \a percentError.
 * \param [in] transform A 4x4 matrix that contains a transform to be applied to points.
 * \param [in] segmentsPerPiece number of segments to sample per contour Piece
 * \param [in] vertexWeldThreshold threshold for welding vertices of adjacent curves
 * \param [in] percentError An error tolerance (percent of lgneth) used in non-uniform curve linearization.
 * \param [out] m user-supplied pointer to point to the mesh object
 * \param [out] revolvedVolume An approximation of the revolved volume of the contour
 *                             or 0 if it could not be computed.
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function
 *
 * \return status set to zero on success, or to a non-zero value otherwise
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::SEGMENT
 * \post revolvedVolume > 0 if it could be computed.
 *
 * \see MFEMReader
 * \see LinearizeCurves
 */
int read_mfem_mesh(const std::string& file,
                   bool uniform,
                   const numerics::Matrix<double>& transform,
                   int segmentsPerPiece,
                   double vertexWeldThreshold,
                   double percentError,
                   mint::Mesh*& m,
                   double& revolvedVolume);
#endif  // AXOM_USE_MFEM

/*!
 * \brief Reads in the Pro/E tetrahedral mesh from the specified file.
 *
 * \param [in] file the file consisting of the Pro/E mesh
 * \param [out] m user-supplied pointer to point to the mesh object.
 * \param [in] comm the MPI communicator, only applicable when MPI is available.
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function.
 *
 * \return zero on success, or a non-zero value otherwise.
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::TET
 *
 * \see ProEReader
 * \see PProEReader
 */
int read_pro_e_mesh(const std::string& file, mint::Mesh*& m, MPI_Comm comm = MPI_COMM_SELF);

/// @}

/// \name Mesh Helper Methods
/// @{

/*!
 * \brief Computes the bounds of the given mesh.
 *
 * \param [in] mesh pointer to the mesh whose bounds will be computed.
 * \param [out] lo buffer to store the lower bound mesh coordinates
 * \param [out] hi buffer to store the upper bound mesh coordinates
 *
 * \pre mesh != nullptr
 * \pre lo != nullptr
 * \pre hi != nullptr
 * \pre hi & lo must point to buffers that are at least N long, where N
 *  corresponds to the mesh dimension.
 */
void compute_mesh_bounds(const mint::Mesh* mesh, double* lo, double* hi);
/// @}

/// \name Logger Initialize/Finalize Methods
/// @{

/*!
 * \brief Helper method to initialize the Slic logger if needed.
 *
 * \param [in,out] isInitialized indicates if Slic is already initialized.
 * \param [out] mustFinalize inidicates if the caller would be responsible
 *  for finalizing the Slic logger.
 * \param [in] verbose flag to control the verbosity
 * \param [in] comm the MPI communicator (applicable when compiled with MPI)
 *
 * \note If Slic is not already initialized, this method will initialize the
 *  Slic Logging environment and set the `isInitialized` flag to true.
 *
 * \note The 'verbose' flag is only applicable when the Slic logging environment
 *  is not already initialized by the calling application. In that case, when
 *  'verbose' is true, all messages will get logged to the console, including,
 *  Info and debug messages. Otherwise, if 'false', only errors will be printed
 *  out.
 *
 *  \see logger_finalize
 */
void logger_init(bool& isInitialized, bool& mustFinalize, bool verbose, MPI_Comm comm);

/*!
 * \brief Finalizes the Slic logger (if needed)
 *
 * \param [in] mustFinalize flag that indicates whether the query is responsible
 *  for finalizing the Slic logger.
 *
 * \see logger_init
 */
void logger_finalize(bool mustFinalize);
/// @}

} /* end namespace internal */
} /* end namespace quest    */
} /* end namespace axom     */

#endif /* QUEST_HELPERS_HPP_ */
