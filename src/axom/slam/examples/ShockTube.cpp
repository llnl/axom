// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 *
 * \brief 1D shock tube, split flux Euler equations
 *
 * \author J. Keasler (original)
 * \author K. Weiss (modified to use axom's Slam component)
 *
 * \details  Developing example to use and demo features of Slam on shock tube
 *  example over structured 1D mesh.
 *
 *  Tests:
 *      * Sets and subsets.
 *      * Relations over regular grid (should be implicit; currently
 *        implemented as explicit static constant relations)
 *      * Fields/maps over the data -- and access to a local datastore
 *        (using Slam::FieldRegistry, can easily convert to Sidre).
 *
 * \verbatim
 *         | m  |            |    mv    |
 *     Q = | mv |        F = | mv^2 + P |
 *         | E  |            |  v(E+P)  |
 *
 *     P = (gamma - 1.0)[E - 0.5 mv^2 ]
 *
 *             Cp
 *     gamma = --     m = mass/volume   v = velocity
 *             Cv
 *
 *     All quantities are non-dimensionalized.
 *
 *     @Q   @F    @Q   @F @Q
 *     -- + -- =  -- + -- -- = 0
 *     @t   @x    @t   @Q @x
 * \endverbatim
 */

#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>

namespace slamShocktube
{
namespace slam = axom::slam;

using PositionType = slam::DefaultPositionType;
using ElementType = slam::DefaultElementType;

using BaseSet = slam::Set<PositionType, ElementType>;

PositionType const UPWIND = 0;
PositionType const DOWNWIND = 1;

const double gammaa = M_SQRT2;
const double gammaaInverse = M_SQRT1_2;

const int INIT_NUM_ELEMS = 100;
const int INIT_NUM_OUTPUT_DUMPS = 5;
const int INIT_NUM_CYCLES_PER_DUMP = 200;

const double INIT_P_RATIO = 0.5;
const double INIT_D_RATIO = 0.5;

#ifdef AXOM_DEBUG
const bool verboseOutput = false;
#endif

/**
 * \brief Simple representation of the mesh for this 1D example
 *
 * \details Mesh contains a set of elements and a set of faces between elements.
 *  It also contains three subsets: a single inflow; a single outflow
 *  element; all internal 'tube' elements. (added 5/2015)
 *  The mesh contains the relations from faces to elements and from tube
 *  elements to faces.
 *
 * \note (5/2015) We are currently missing an implicit constant grid relation.
 *
 * \note For current implementation with explicit static (constant) relations.
 * \note We are missing a nice way to set the relation elements.
 *  It should not have to be done explicitly in each user code --especially
    in common use cases
 *
 * Idea: We could have a relationInverter function that takes a relation
 *       from sets A to B and generates a relation from set B to set A with all
         the arrows reversed.
 */
class ShockTubeMesh
{
public:
  /// types for Element and Face sets
  using PositionType = slamShocktube::PositionType;

  using ElemSet = slam::PositionSet<PositionType, ElementType>;
  using FaceSet = slam::PositionSet<PositionType, ElementType>;

  using IndexType = PositionType;

  /// types for Tube and {In,Out}Flow subsets
  using RuntimeOffsetPolicy = slam::policies::RuntimeOffset<PositionType>;
  using StrideOnePolicy = slam::policies::StrideOne<PositionType>;
  using NoIndirectionPolicy = slam::policies::NoIndirection<PositionType, ElementType>;
  using TubeSubsetPolicy = slam::policies::ConcreteParentSubset<ElemSet>;
  using ElemSubset =
    slam::GenericRangeSet<PositionType, ElementType, RuntimeOffsetPolicy, StrideOnePolicy, NoIndirectionPolicy, TubeSubsetPolicy>;
  using RangeSet = slam::RangeSet<PositionType, ElementType>;

  /// types for relations
  enum
  {
    ELEMS_PER_FACE = 2,
    FACES_PER_ELEM = 2
  };

  using EFStride = slam::policies::CompileTimeStride<PositionType, FACES_PER_ELEM>;
  using FEStride = slam::policies::CompileTimeStride<PositionType, ELEMS_PER_FACE>;

  using EFCard = slam::policies::ConstantCardinality<PositionType, EFStride>;
  using FECard = slam::policies::ConstantCardinality<PositionType, FEStride>;
  using STLIndirection = slam::policies::STLVectorIndirection<PositionType, ElementType>;
  using IndexVec = typename STLIndirection::IndirectionBufferType;

  using TubeElemToFaceRelation =
    slam::StaticRelation<PositionType, ElementType, EFCard, STLIndirection, ElemSubset, FaceSet>;
  using FaceToElemRelation =
    slam::StaticRelation<PositionType, ElementType, FECard, STLIndirection, FaceSet, ElemSet>;

public:
  ElemSet elems;            // The entire set of elements
  ElemSubset tubeElems;     // Subset of internal elements
  ElemSubset inFlowElems;   // Subset of inflow elements
  ElemSubset outFlowElems;  // Subset of outflow elements
  FaceSet faces;            // Faces between adjacent pairs of elements

  FaceToElemRelation relationFaceElem;      // co-boundary relation of faces
                                            // to their elements
  TubeElemToFaceRelation relationTubeFace;  // boundary relation of
                                            // internal 'tube' elements
};

// Define explicit instances of local (key/value) datastore for int and double
using IntsRegistry = slam::FieldRegistry<BaseSet, BaseSet::ElementType>;
using RealsRegistry = slam::FieldRegistry<BaseSet, double>;
using IntField = IntsRegistry::MapType;
using RealField = RealsRegistry::MapType;

IntsRegistry intsRegistry;
RealsRegistry realsRegistry;

/**************************************************************************
 * Subroutine:  GetUserInput
 * Purpose   :  Ask for control and output information
 *************************************************************************/

void GetUserInput()
{
  /**********************************/
  /* Get mesh info, and create mesh */
  /**********************************/
  {
    int numElems;
    int numFaces;

    SLIC_INFO("How many zones for the 1D shock tube? ");
    numElems = INIT_NUM_ELEMS;
    SLIC_INFO("\t\t" << numElems);

    // add an inflow and outflow zone
    numElems += 2;
    numFaces = numElems - 1;

    intsRegistry.addScalar("numElems", numElems);
    intsRegistry.addScalar("numFaces", numFaces);
  }

  /********************/
  /* Get physics info */
  /********************/
  {
    double pratio = -1.0;
    double dratio = -1.0;

    while(pratio < 0.0 || pratio > 1.0)
    {
      SLIC_INFO("What pressure ratio would you like (0 <= x <= 1)? ");
      pratio = INIT_P_RATIO;
      SLIC_INFO("\t\t" << pratio);
    }

    while(dratio < 0.0 || dratio > 1.0)
    {
      SLIC_INFO("What density ratio would you like (0 <= x <= 1)? ");
      dratio = INIT_D_RATIO;
      SLIC_INFO("\t\t" << dratio);
    }

    realsRegistry.addScalar("pressureRatio", pratio);
    realsRegistry.addScalar("densityRatio", dratio);
  }

  /********************/
  /* Get output  info */
  /********************/
  {
    int numOutputDumps;
    int numCyclesPerDump;

    SLIC_INFO("How many dumps would you like? ");
    numOutputDumps = INIT_NUM_OUTPUT_DUMPS;
    SLIC_INFO("\t\t" << numOutputDumps);

    SLIC_INFO("How many cycles between dumps would you like? ");
    numCyclesPerDump = INIT_NUM_CYCLES_PER_DUMP;
    SLIC_INFO("\t\t" << numCyclesPerDump);

    int numTotalCycles = numOutputDumps * numCyclesPerDump;
    SLIC_INFO("Simulation will run for " << numTotalCycles << " cycles.\n");

    intsRegistry.addScalar("numOutputDumps", numOutputDumps);
    intsRegistry.addScalar("numCyclesPerDump", numCyclesPerDump);
    intsRegistry.addScalar("numTotalCycles", numTotalCycles);
  }

  return;
}

/**
 * \brief Build an empty mesh for the shock tube
 * \details Shocktube mesh layout
 * \verbatim
 *      Gaps between elements are faces
 *                     |
 *      -------------------------------
 *      |   |   |             |   |   |
 *
 *   ### ### ### ###       ### ### ### ###
 *   ### ### ### ###  ...  ### ### ### ###  <--- 1D Shock tube model
 *   ### ### ### ###       ### ### ### ###
 *
 *    |  |                           |  |
 *    |  -----------------------------  |
 *   Inflow           |               Outflow
 *   Element      Tube Elements       Element
 *
 * \endverbatim
 */
void CreateShockTubeMesh(ShockTubeMesh* mesh)
{
  // ------------ Generate the Sets and Subsets

  // create element and face sets
  mesh->elems = ShockTubeMesh::ElemSet(intsRegistry.getScalar("numElems"));
  mesh->faces = ShockTubeMesh::FaceSet(intsRegistry.getScalar("numFaces"));

  // construct the element subsets
  ShockTubeMesh::PositionType numElems = mesh->elems.size();
  using ElemSubsetBuilder = ShockTubeMesh::ElemSubset::SetBuilder;
  mesh->inFlowElems = ElemSubsetBuilder()  //
                        .range(0, 1)       //
                        .parent(&mesh->elems);
  mesh->tubeElems = ElemSubsetBuilder()        //
                      .range(1, numElems - 1)  //
                      .parent(&mesh->elems);
  mesh->outFlowElems = ElemSubsetBuilder()               //
                         .range(numElems - 1, numElems)  //
                         .parent(&mesh->elems);

  // ------------ Generate the Relations

  // TODO: Need to define DynamicConstantRelation
  //       which will allow modifying the relation elements
  // TODO: Need to define ImplicitConstantRelation
  //       since the relations are actually implicit
  //       no storage should be necessary for regular grid neighbors
  //       For now, we use explicitly storage for the relation data.

  using IndexVec = ShockTubeMesh::IndexVec;

  /// Setup the FaceToElem relation
  IndexVec& feRelVec =
    intsRegistry.addBuffer("feRel", ShockTubeMesh::FACES_PER_ELEM * mesh->faces.size());
  IndexVec::iterator relIt = feRelVec.begin();
  for(ShockTubeMesh::IndexType idx = 0;
      idx < static_cast<ShockTubeMesh::IndexType>(mesh->faces.size());
      ++idx)
  {
    *relIt++ = mesh->faces[idx];
    *relIt++ = mesh->faces[idx] + 1;
  }

  mesh->relationFaceElem = ShockTubeMesh::FaceToElemRelation(&mesh->faces, &mesh->elems);
  mesh->relationFaceElem.bindIndices(static_cast<int>(feRelVec.size()), &feRelVec);
  SLIC_ASSERT(mesh->relationFaceElem.isValid(verboseOutput));

  /// Setup the TubeElementToFace relation:
  /// A relation from the tubes subset of the elements to their incident faces
  ShockTubeMesh::PositionType numTubeElems = mesh->tubeElems.size();
  IndexVec& efRelVec = intsRegistry.addBuffer("efRel", ShockTubeMesh::ELEMS_PER_FACE * numTubeElems);
  relIt = efRelVec.begin();
  for(ShockTubeMesh::IndexType idx = 0; idx < static_cast<ShockTubeMesh::IndexType>(numTubeElems);
      ++idx)
  {
    *relIt++ = mesh->tubeElems[idx] - 1;
    *relIt++ = mesh->tubeElems[idx];
  }

  mesh->relationTubeFace = ShockTubeMesh::TubeElemToFaceRelation(&mesh->tubeElems, &mesh->faces);
  mesh->relationTubeFace.bindIndices(static_cast<int>(efRelVec.size()), &efRelVec);
  SLIC_ASSERT(mesh->relationTubeFace.isValid(verboseOutput));
}

/**************************************************************************
 * Subroutine:  InitializeShockTube
 * Purpose   :  Populate the mesh with values
 *************************************************************************/
void InitializeShockTube(ShockTubeMesh const& mesh)
{
  using IndexType = ShockTubeMesh::IndexType;
  using RangeSet = ShockTubeMesh::RangeSet;

  // Create element centered fields
  RealField& mass = realsRegistry.addField("mass", &mesh.elems);
  RealField& momentum = realsRegistry.addField("momentum", &mesh.elems);
  RealField& energy = realsRegistry.addField("energy", &mesh.elems);
  RealField& pressure = realsRegistry.addField("pressure", &mesh.elems);

  // Create face centered fields
  realsRegistry.addField("F0", &mesh.faces);  // mv
  realsRegistry.addField("F1", &mesh.faces);  // mv^2+P
  realsRegistry.addField("F2", &mesh.faces);  // v(E+P)

  // Fill left half with high pressure, right half with low pressure
  IndexType endTube = mesh.elems.size();
  IndexType midTube = endTube / 2;

  // Non-dimensionalized reference values
  double massInitial = 1.0;
  double momentumInitial = 0.0;
  double pressureInitial = gammaaInverse;
  double energyInitial = pressureInitial / (gammaa - 1.0);

  // Initialize zonal quantities
  RangeSet lowerTube(0, midTube);
  for(IndexType i = 0; i < lowerTube.size(); ++i)
  {
    IndexType ind = lowerTube[i];
    mass[ind] = massInitial;
    momentum[ind] = momentumInitial;
    pressure[ind] = pressureInitial;
    energy[ind] = energyInitial;
  }

  // adjust parameters for low pressure portion of tube
  massInitial *= realsRegistry.getScalar("densityRatio");
  pressureInitial *= realsRegistry.getScalar("pressureRatio");
  energyInitial = pressureInitial / (gammaa - 1.0);

  RangeSet upperTube(midTube, mesh.elems.size());
  for(IndexType i = 0; i < upperTube.size(); ++i)
  {
    IndexType ind = upperTube[i];
    mass[ind] = massInitial;
    momentum[ind] = momentumInitial;
    pressure[ind] = pressureInitial;
    energy[ind] = energyInitial;
  }

  // Create needed time info
  realsRegistry.addScalar("time", 0.0);
  intsRegistry.addScalar("cycle", 0);

  double dx = 1.0 / static_cast<double>(endTube);
  realsRegistry.addScalar("dx", dx);
  realsRegistry.addScalar("dt", 0.4 * dx);
}

/**
 * \brief Compute F quantities at faces.
 *
 * \details
 *
 * \verbatim
 *  @F   @F0   @F1   @F2
 *  -- = --- + --- + ---
 *  @x   @x    @x    @x
 *  \endverbatim
 *
 *  Calculate F0, F1 and F2 at the face centers.
 *
 */
void ComputeFaceInfo(ShockTubeMesh const& mesh)
{
  using IndexType = ShockTubeMesh::IndexType;

  // Face fields
  RealField& F0 = realsRegistry.getField("F0");
  RealField& F1 = realsRegistry.getField("F1");
  RealField& F2 = realsRegistry.getField("F2");

  // Element fields
  const RealField& mass = realsRegistry.getField("mass");
  const RealField& momentum = realsRegistry.getField("momentum");
  const RealField& energy = realsRegistry.getField("energy");

  // Update face data using element data using the face->elem relation
  ShockTubeMesh::PositionType numFaceElems =
    static_cast<ShockTubeMesh::PositionType>(mesh.faces.size());
  for(ShockTubeMesh::PositionType fIdx = 0; fIdx < numFaceElems; ++fIdx)
  {
    // each face has an upwind and downwind element.
    IndexType upWind = mesh.relationFaceElem[fIdx][UPWIND];      // upwind
                                                                 // element
    IndexType downWind = mesh.relationFaceElem[fIdx][DOWNWIND];  // downwind
                                                                 // element

    // calculate face centered quantities as avg of element centered ones
    double massf = 0.5 * (mass[upWind] + mass[downWind]);
    double momentumf = 0.5 * (momentum[upWind] + momentum[downWind]);
    double energyf = 0.5 * (energy[upWind] + energy[downWind]);
    double pressuref = (gammaa - 1.0) * (energyf - 0.5 * momentumf * momentumf / massf);
    double c = sqrt(gammaa * pressuref / massf);
    double v = momentumf / massf;

    double ev;
    double cLocal;

    // Now that we have the wave speeds, we might want to
    // look for the max wave speed here, and update dt
    // appropriately right before leaving this function.

    // OK, calculate face quantities
    F0[fIdx] = F1[fIdx] = F2[fIdx] = 0.0;

    IndexType contributor = ((v >= 0.0) ? upWind : downWind);
    massf = mass[contributor];
    momentumf = momentum[contributor];
    energyf = energy[contributor];
    pressuref = energyf - 0.5 * momentumf * momentumf / massf;
    ev = v * (gammaa - 1.0);

    F0[fIdx] += ev * massf;
    F1[fIdx] += ev * momentumf;
    F2[fIdx] += ev * (energyf - pressuref);

    contributor = ((v + c >= 0.0) ? upWind : downWind);
    massf = mass[contributor];
    momentumf = momentum[contributor];
    energyf = energy[contributor];
    pressuref = (gammaa - 1.0) * (energyf - 0.5 * momentumf * momentumf / massf);
    ev = 0.5 * (v + c);
    cLocal = sqrt(gammaa * pressuref / massf);

    F0[fIdx] += ev * massf;
    F1[fIdx] += ev * (momentumf + massf * cLocal);
    F2[fIdx] += ev * (energyf + pressuref + momentumf * cLocal);

    contributor = ((v - c >= 0.0) ? upWind : downWind);
    massf = mass[contributor];
    momentumf = momentum[contributor];
    energyf = energy[contributor];
    pressuref = (gammaa - 1.0) * (energyf - 0.5 * momentumf * momentumf / massf);
    ev = 0.5 * (v - c);
    cLocal = sqrt(gammaa * pressuref / massf);

    F0[fIdx] += ev * massf;
    F1[fIdx] += ev * (momentumf - massf * cLocal);
    F2[fIdx] += ev * (energyf + pressuref - momentumf * cLocal);
  }
}

/**************************************************************************
 * Subroutine:  UpdateElemInfo
 * Purpose   :  Q(elem) = Q(elem) + deltaQ(elem)
 *
 *  deltaQ(elem) = - (F(downWindFace) - F(upWindFace)) * dt / dx ;
 *
 *************************************************************************/

void UpdateElemInfo(ShockTubeMesh const& mesh)
{
  // get the element quantities we want to update
  RealField& mass = realsRegistry.getField("mass");
  RealField& momentum = realsRegistry.getField("momentum");
  RealField& energy = realsRegistry.getField("energy");
  RealField& pressure = realsRegistry.getField("pressure");

  // The element update is calculated as the flux between faces
  const RealField& F0 = realsRegistry.getField("F0");
  const RealField& F1 = realsRegistry.getField("F1");
  const RealField& F2 = realsRegistry.getField("F2");

  double dx = realsRegistry.getScalar("dx");
  double dt = realsRegistry.getScalar("dt");
  double& time = realsRegistry.getScalar("time");

  /// Update the element fields based on the face data using the elem->face
  // relation
  ShockTubeMesh::PositionType numTubeElems =
    static_cast<ShockTubeMesh::PositionType>(mesh.tubeElems.size());

  for(ShockTubeMesh::PositionType tPos = 0; tPos < numTubeElems; ++tPos)
  {
    // Relation is over tube elements.
    ShockTubeMesh::IndexType elemIdx = mesh.tubeElems[tPos];

    // Each element inside the tube has an upwind and downwind face
    ShockTubeMesh::IndexType upWind = mesh.relationTubeFace[tPos][UPWIND];      // upwind
                                                                                // face
    ShockTubeMesh::IndexType downWind = mesh.relationTubeFace[tPos][DOWNWIND];  // downwind
                                                                                // face

    mass[elemIdx] -= gammaaInverse * (F0[downWind] - F0[upWind]) * dt / dx;
    momentum[elemIdx] -= gammaaInverse * (F1[downWind] - F1[upWind]) * dt / dx;
    energy[elemIdx] -= gammaaInverse * (F2[downWind] - F2[upWind]) * dt / dx;
    pressure[elemIdx] = (gammaa - 1.0) *
      (energy[elemIdx] - 0.5 * momentum[elemIdx] * momentum[elemIdx] / mass[elemIdx]);
  }

  // update the time
  time += dt;
}

void dumpData(ShockTubeMesh const& mesh)
{
  RealField const& mass = realsRegistry.getField("mass");
  RealField const& momentum = realsRegistry.getField("momentum");
  RealField const& energy = realsRegistry.getField("energy");
  RealField const& pressure = realsRegistry.getField("pressure");

  static const ShockTubeMesh::PositionType MAX_ELEM_DUMP = 20;
  const int maxDumpPerSide = std::min(mesh.elems.size(), MAX_ELEM_DUMP) / 2;

  using ElemSubsetBuilder = ShockTubeMesh::ElemSubset::SetBuilder;
  ShockTubeMesh::ElemSubset begSet, endSet;
  begSet = ElemSubsetBuilder()     //
             .parent(&mesh.elems)  //
             .size(maxDumpPerSide);
  endSet =
    ElemSubsetBuilder().parent(&mesh.elems).size(maxDumpPerSide).offset(mesh.elems.size() - maxDumpPerSide);

  SLIC_ASSERT(begSet.isValid(verboseOutput));
  SLIC_ASSERT(endSet.isValid(verboseOutput));

  // Use subsets to output a few samples from the beginning and end of the tube
  std::stringstream elemStream, mStream, pStream, eStream, rStream;

  elemStream << std::setw(5) << std::setfill(' ');
  mStream << std::setw(5) << std::setfill(' ') << std::setprecision(3);
  pStream << std::setw(5) << std::setfill(' ') << std::setprecision(3);
  eStream << std::setw(5) << std::setfill(' ') << std::setprecision(3);
  rStream << std::setw(5) << std::setfill(' ') << std::setprecision(3);

  for(int i = 0; i < begSet.size(); ++i)
  {
    ShockTubeMesh::IndexType ind = begSet[i];
    if(i == 0)
    {
      elemStream << "IN"
                 << "\t";
    }
    else
    {
      elemStream << ind << "\t";
    }

    mStream << mass[ind] << "\t";
    pStream << momentum[ind] << "\t";
    eStream << energy[ind] << "\t";
    rStream << pressure[ind] << "\t";
  }

  elemStream << "...\t";
  mStream << "...\t";
  pStream << "...\t";
  eStream << "...\t";
  rStream << "...\t";

  for(int i = 0; i < endSet.size(); ++i)
  {
    ShockTubeMesh::IndexType ind = endSet[i];
    if(ind == endSet.parentSet()->size() - 1)
    {
      elemStream << "OUT"
                 << "\t";
    }
    else
    {
      elemStream << ind << "\t";
    }

    mStream << mass[ind] << "\t";
    pStream << momentum[ind] << "\t";
    eStream << energy[ind] << "\t";
    rStream << pressure[ind] << "\t";
  }

  SLIC_INFO("Data dump: \n"
            << "Elem idx: " << elemStream.str() << "\n"
            << "mass:     " << mStream.str() << "\n"
            << "momemtum: " << pStream.str() << "\n"
            << "energy:   " << eStream.str() << "\n"
            << "pressure: " << rStream.str() << "\n");
}

}  // end namespace slamShocktube

/**************************************************************************
 * Subroutine:  main
 * Purpose   :  Simulate a 1D Shock Tube using split flux Euler formulation
 *************************************************************************/

int main(void)
{
  using namespace slamShocktube;
  axom::slic::SimpleLogger logger;

  // We should be able to parallelize pretty easily by
  // adding an MPI_Init() here, modifying the setup slightly,
  // adding a communication subroutine in the main loop,
  // and calling MPI_Finalize() at the end of main()

  ShockTubeMesh mesh;

  GetUserInput();

  CreateShockTubeMesh(&mesh);  // setup sets and relations
  InitializeShockTube(mesh);   // setup fields

  int numTotalCycles = intsRegistry.getScalar("numTotalCycles");
  int dumpInterval = intsRegistry.getScalar("numCyclesPerDump");

  // use the & operation when you want to update the param directly
  auto& currCycle = intsRegistry.getScalar("cycle");
  for(currCycle = 0; currCycle < numTotalCycles; ++currCycle)
  {
    if(currCycle % dumpInterval == 0)
    {
      SLIC_INFO("\tStarting cycle " << currCycle << " at time " << realsRegistry.getScalar("time"));
      dumpData(mesh);
    }

    ComputeFaceInfo(mesh);
    UpdateElemInfo(mesh);
  }

  SLIC_INFO("\tFinished cycle " << currCycle << " at time " << realsRegistry.getScalar("time"));

  dumpData(mesh);

  SLIC_INFO("done.");

  return 0;
}
