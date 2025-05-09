// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>

#include <map>

#include "lulesh.hpp"

#include "axom/core/Macros.hpp"

#ifdef AXOM_USE_MPI
#include <mpi.h>
#endif


namespace slamLulesh {

/* Helper function for converting strings to ints, with error checking */
  int StrToInt(const char *token, Int_t *retVal)
  {
    const char *c;
    char *endptr;
    const int decimal_base = 10;

    if (token == NULL)
      return 0;

    c = token;
    *retVal = (int)strtol(c, &endptr, decimal_base);
    if((endptr != c) && ((*endptr == ' ') || (*endptr == '\0')))
      return 1;
    else
      return 0;
  }

  static void PrintCommandLineOptions(char *execname, int myRank)
  {
    if (myRank == 0)
    {

      SLIC_INFO( "Usage: "
          << execname << " [opts]\n"
          << " where [opts] is one or more of:\n"
          << " -q              : quiet mode - suppress all stdout\n"
          << " -i <iterations> : number of cycles to run\n"
          << " -s <size>       : length of cube mesh along side\n"
          << " -r <numregions> : Number of distinct regions (def: 11)\n"
          << " -b <balance>    : Load balance between regions of a domain (def: 1)\n"
          << " -c <cost>       : Extra cost of more expensive regions (def: 1)\n"
          << " -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n"
          << " -p              : Print out progress\n"
          << " -v              : Output viz file (requires compiling with -DVIZ_MESH\n"
          << " -h              : This message\n"
          << "\n\n");
    }
    axom::slic::flushStreams();
  }

  static void ParseError(const char *message, int myRank)
  {
    if (myRank == 0)
    {
      SLIC_WARNING(message);
    }
    axom::slic::flushStreams();
#ifdef AXOM_USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#else
    exit(-1);
#endif
  }

  void ParseCommandLineOptions(int argc, char *argv[],
      int myRank, struct cmdLineOpts *opts)
  {
    if(argc > 1)
    {
      int i = 1;

      while(i < argc)
      {
        int ok;
        /* -i <iterations> */
        if(strcmp(argv[i], "-i") == 0)
        {
          if (i + 1 >= argc)
          {
            ParseError("Missing integer argument to -i", myRank);
          }
          ok = StrToInt(argv[i + 1], &(opts->its));
          if(!ok)
          {
            ParseError("Parse Error on option -i integer value required after argument\n", myRank);
          }
          i += 2;
        }
        /* -s <size, sidelength> */
        else if(strcmp(argv[i], "-s") == 0)
        {
          if (i + 1 >= argc)
          {
            ParseError("Missing integer argument to -s\n", myRank);
          }
          ok = StrToInt(argv[i + 1], &(opts->nx));
          if(!ok)
          {
            ParseError("Parse Error on option -s integer value required after argument\n", myRank);
          }
          i += 2;
        }
        /* -r <numregions> */
        else if (strcmp(argv[i], "-r") == 0)
        {
          if (i + 1 >= argc)
          {
            ParseError("Missing integer argument to -r\n", myRank);
          }
          ok = StrToInt(argv[i + 1], &(opts->numReg));
          if (!ok)
          {
            ParseError("Parse Error on option -r integer value required after argument\n", myRank);
          }
          i += 2;
        }
        /* -f <numfilepieces> */
        else if (strcmp(argv[i], "-f") == 0)
        {
          if (i + 1 >= argc)
          {
            ParseError("Missing integer argument to -f\n", myRank);
          }
          ok = StrToInt(argv[i + 1], &(opts->numFiles));
          if (!ok)
          {
            ParseError("Parse Error on option -f integer value required after argument\n", myRank);
          }
          i += 2;
        }
        /* -p */
        else if (strcmp(argv[i], "-p") == 0)
        {
          opts->showProg = 1;
          i++;
        }
        /* -q */
        else if (strcmp(argv[i], "-q") == 0)
        {
          opts->quiet = 1;
          i++;
        }
        else if (strcmp(argv[i], "-b") == 0)
        {
          if (i + 1 >= argc)
          {
            ParseError("Missing integer argument to -b\n", myRank);
          }
          ok = StrToInt(argv[i + 1], &(opts->balance));
          if (!ok)
          {
            ParseError("Parse Error on option -b integer value required after argument\n", myRank);
          }
          i += 2;
        }
        else if (strcmp(argv[i], "-c") == 0)
        {
          if (i + 1 >= argc)
          {
            ParseError("Missing integer argument to -c\n", myRank);
          }
          ok = StrToInt(argv[i + 1], &(opts->cost));
          if (!ok)
          {
            ParseError("Parse Error on option -c integer value required after argument\n", myRank);
          }
          i += 2;
        }
        /* -v */
        else if (strcmp(argv[i], "-v") == 0)
        {
#if VIZ_MESH
          opts->viz = 1;
#else
          ParseError("Use of -v requires compiling with -DVIZ_MESH\n", myRank);
#endif
          i++;
        }
        /* -h */
        else if (strcmp(argv[i], "-h") == 0)
        {
          PrintCommandLineOptions(argv[0], myRank);
#ifdef AXOM_USE_MPI
          MPI_Abort(MPI_COMM_WORLD, 0);
#else
          exit(0);
#endif
        }
        else {
          char msg[80];
          PrintCommandLineOptions(argv[0], myRank);
          sprintf(msg, "ERROR: Unknown command line argument: %s\n", argv[i]);
          ParseError(msg, myRank);
        }
      }
    }
  }

/////////////////////////////////////////////////////////////////////

  void VerifyAndWriteFinalOutput(Real_t elapsed_time,
      Domain& locDom,
      Int_t nx,
      Int_t numRanks)
  {
    // GrindTime1 only takes a single domain into account, and is thus a good way to measure
    // processor speed independent of MPI parallelism.
    // GrindTime2 takes into account speedups from MPI parallelism
    Real_t grindTime1 = ((elapsed_time * 1e6) / locDom.cycle()) / (nx * nx * nx);
    Real_t grindTime2 = ((elapsed_time * 1e6) / locDom.cycle()) / (nx * nx * nx * numRanks);

    Index_t ElemId = 0;

    SLIC_INFO( "Run completed:"
        << "\n\tProblem size        =  " << nx
        << "\n\tMPI tasks           =  " << numRanks
        << "\n\tIteration count     =  " << locDom.cycle()
        << "\n\tFinal Origin Energy =  " << locDom.e(ElemId)
        << "\n");

    Real_t MaxAbsDiff = Real_t(0.0);
    Real_t TotalAbsDiff = Real_t(0.0);
    Real_t MaxRelDiff = Real_t(0.0);

    for (Index_t j = 0; j<nx; ++j)
    {
      for (Index_t k = j + 1; k<nx; ++k)
      {
        Real_t AbsDiff = FABS(locDom.e(j * nx + k) - locDom.e(k * nx + j));
        TotalAbsDiff  += AbsDiff;

        if (MaxAbsDiff <AbsDiff)
          MaxAbsDiff = AbsDiff;

        Real_t RelDiff = AbsDiff / locDom.e(k * nx + j);

        if (MaxRelDiff <RelDiff)
          MaxRelDiff = RelDiff;
      }
    }

    // Quick symmetry check
    SLIC_INFO( "Testing plane 0 of energy array on rank 0:"
        << "\n\tMaxAbsDiff   = " << MaxAbsDiff
        << "\n\tTotalAbsDiff = " << TotalAbsDiff
        << "\n\tMaxRelDiff   = " << MaxRelDiff
        << "\n");

    // Timing information
    SLIC_INFO( "Elapsed time         = "
        << elapsed_time << "(s)"
        << "\n\tGrind time (us/z/c)  = " << grindTime1 << " (per dom)  (" << grindTime2 << " overall)"
        << "\n\tFOM                  = " << 1000.0 / grindTime2 << " (z/s)"           // zones per second
        << "\n");

    /// Compare to given values for cycles and origin energy from lulesh 2.0 paper:
    std::map<Int_t, std::pair<Index_t, Real_t> > resultCheckMap;

    resultCheckMap[ 5] =  std::make_pair(   72, 7.853665e+03);
    resultCheckMap[10] =  std::make_pair(  231, 2.720531e+04);
    resultCheckMap[30] =  std::make_pair(  932, 2.025075e+05);
    resultCheckMap[45] =  std::make_pair( 1477, 4.234875e+05);
    resultCheckMap[50] =  std::make_pair( 1662, 5.124778e+05);
    resultCheckMap[70] =  std::make_pair( 2402, 9.417145e+05);
    resultCheckMap[90] =  std::make_pair( 3145, 1.482403e+06);


    // find the overall number of edges on the problem domain (domains per edge * num elems per domain
    Int_t domainsPerSide = Int_t(cbrt(Real_t(numRanks)) + 0.5);
    Int_t gEdge = nx * domainsPerSide;
    if(resultCheckMap.find(gEdge) != resultCheckMap.end() )
    {
      SLIC_CHECK_MSG( resultCheckMap[gEdge].first == locDom.cycle(),
          "Specs state that num cycles should be "
          << resultCheckMap[gEdge].first
          << " actual number of cycles was " << locDom.cycle() << "." );

      SLIC_CHECK_MSG(
        axom::utilities::isNearlyEqualRelative( resultCheckMap[gEdge].second, locDom.e(ElemId)),
        "Specs state that final energy at origin must be "
        << resultCheckMap[gEdge].second
        << " actual energy at origin was " << locDom.e(ElemId)
        << ". Difference was " << std::fabs(resultCheckMap[gEdge].second - locDom.e(ElemId) ) );

      double diff = std::fabs(resultCheckMap[gEdge].second - locDom.e(ElemId) );
      double maxFabs = std::max( std::fabs(resultCheckMap[gEdge].second), std::fabs(locDom.e(ElemId) ) );
      double relMaxFabs = 1.0e-6 * maxFabs;
      double relMaxFabsWithAbsolute = relMaxFabs + 1.0e-8;

      AXOM_UNUSED_VAR( diff);
      AXOM_UNUSED_VAR( maxFabs);
      AXOM_UNUSED_VAR( relMaxFabs);
      AXOM_UNUSED_VAR( relMaxFabsWithAbsolute);
      SLIC_DEBUG("**  comparing "
          << resultCheckMap[gEdge].second << " with " << locDom.e(ElemId)
          << "\n\tfabs difference: " << diff
          << "\n\tmaxFabs * rel (1e-6): " <<  relMaxFabs
          << "\n\t<above> with abs 1e-8: " << relMaxFabsWithAbsolute
          << "\n\tdiff (again): " << diff
          << "\n\tdiff of last two: " << relMaxFabsWithAbsolute - diff
          << "\n\tNearly equal: " << ( diff <= relMaxFabsWithAbsolute ? "TRUE" : "FALSE" )
      );
    }

    return;
  }

} // end namespace slamLulesh
