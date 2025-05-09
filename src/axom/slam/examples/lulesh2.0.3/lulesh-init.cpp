// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "lulesh.hpp"

#include <math.h>

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif
#ifdef AXOM_USE_OPENMP
  #include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <cstdlib>
#include <algorithm>


namespace slamLulesh {

/////////////////////////////////////////////////////////////////////
  Domain::Domain(Int_t numRanks, Index_t colLoc,
      Index_t rowLoc, Index_t planeLoc,
      Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost)
      : m_e_cut(Real_t(1.0e-7)),
        m_p_cut(Real_t(1.0e-7)),
        m_q_cut(Real_t(1.0e-7)),
        m_v_cut(Real_t(1.0e-10)),
        m_u_cut(Real_t(1.0e-7)),
        m_hgcoef(Real_t(3.0)),
        m_ss4o3(Real_t(4.0) / Real_t(3.0)),
        m_qstop(Real_t(1.0e+12)),
        m_monoq_max_slope(Real_t(1.0)),
        m_monoq_limiter_mult(Real_t(2.0)),
        m_qlc_monoq(Real_t(0.5)),
        m_qqc_monoq(Real_t(2.0) / Real_t(3.0)),
        m_qqc(Real_t(2.0)),
        m_eosvmax(Real_t(1.0e+9)),
        m_eosvmin(Real_t(1.0e-9)),
        m_pmin(Real_t(0.)),
        m_emin(Real_t(-1.0e+15)),
        m_dvovmax(Real_t(0.1)),
        m_refdens(Real_t(1.0))
  {

    Index_t edgeElems = nx;
    Index_t edgeNodes = edgeElems + 1;

    this->cost() = cost;

    m_tp       = tp;
    m_numRanks = numRanks;

    ///////////////////////////////
    //   Initialize Sedov Mesh
    ///////////////////////////////

    // construct a uniform box for this processor

    m_colLoc   =   colLoc;
    m_rowLoc   =   rowLoc;
    m_planeLoc = planeLoc;

    m_sizeX = edgeElems;
    m_sizeY = edgeElems;
    m_sizeZ = edgeElems;
    m_elemSet = ElemSet(edgeElems * edgeElems * edgeElems);
    m_nodeSet = NodeSet(edgeNodes * edgeNodes * edgeNodes);
    m_cornerSet = CornerSet( m_elemSet.size() * 8);

    Int_t facesPerPlane = edgeElems * edgeElems;
    Int_t numExtendedElems = m_elemSet.size()   // local elem
        + 2 * facesPerPlane                     // plane ghosts
        + 2 * facesPerPlane                     // row ghosts
        + 2 * facesPerPlane                     // col ghosts
    ;

    m_extendedElemSet = ExtendedElemSet( numExtendedElems );


    // Elem-centered
    AllocateElemPersistent(numElem());

    // Node-centered
    AllocateNodePersistent(numNode());

    SetupCommBuffers();

    // Basic Field Initialization
    for (Index_t i = 0; i<numElem(); ++i)
    {
      e(i) =  Real_t(0.0);
      p(i) =  Real_t(0.0);
      q(i) =  Real_t(0.0);
      ss(i) = Real_t(0.0);
    }

    // Note - v initializes to 1.0, not 0.0!
    for (Index_t i = 0; i<numElem(); ++i)
    {
      v(i) = Real_t(1.0);
    }

    for (Index_t i = 0; i<numNode(); ++i)
    {
      xd(i) = Real_t(0.0);
      yd(i) = Real_t(0.0);
      zd(i) = Real_t(0.0);
    }

    for (Index_t i = 0; i<numNode(); ++i)
    {
      xdd(i) = Real_t(0.0);
      ydd(i) = Real_t(0.0);
      zdd(i) = Real_t(0.0);
    }

    for (Index_t i = 0; i<numNode(); ++i)
    {
      nodalMass(i) = Real_t(0.0);
    }

    BuildMesh(nx, edgeNodes, edgeElems);

#ifdef AXOM_USE_OPENMP
    SetupThreadSupportStructures();
#endif

    // Setup region index sets. For now, these are constant sized
    // throughout the run, but could be changed every cycle to
    // simulate effects of ALE on the lagrange solver
    CreateRegionIndexSets(nr, balance);

    // Setup symmetry nodesets
    SetupSymmetryPlanes(edgeNodes);

    // Setup element connectivities
    SetupElementConnectivities(edgeElems);

    // Setup symmetry planes and free surface boundary arrays
    SetupBoundaryConditions(edgeElems);


    // Setup defaults

    // These can be changed (requires recompile) if you want to run
    // with a fixed timestep, or to a different end time, but it's
    // probably easier/better to just run a fixed number of timesteps
    // using the -i flag in 2.x

    dtfixed() = Real_t(-1.0e-6); // Negative means use courant condition
    stoptime()  = Real_t(1.0e-2);// *Real_t(edgeElems*tp/45.0) ;

    // Initial conditions
    deltatimemultlb() = Real_t(1.1);
    deltatimemultub() = Real_t(1.2);
    dtcourant()       = Real_t(1.0e+20);
    dthydro()         = Real_t(1.0e+20);
    dtmax()           = Real_t(1.0e-2);
    time()            = Real_t(0.);
    cycle()           = Int_t(0);

    // initialize field data
    for (Index_t i = 0; i<numElem(); ++i)
    {
      Real_t x_local[8], y_local[8], z_local[8];
      ElemNodeSet elemToNode = nodelist(i);
      for( Index_t lnode = 0; lnode< elemToNode.size(); ++lnode )
      {
        Index_t gnode = elemToNode[lnode];
        x_local[lnode] = x(gnode);
        y_local[lnode] = y(gnode);
        z_local[lnode] = z(gnode);
      }

      // volume calculations
      Real_t volume = CalcElemVolume(x_local, y_local, z_local );
      volo(i) = volume;
      elemMass(i) = volume;
      for (Index_t j = 0; j<8; ++j)
      {
        Index_t idx = elemToNode[j];
        nodalMass(idx) += volume / Real_t(8.0);
      }
    }

    // deposit initial energy
    // An energy of 3.948746e+7 is correct for a problem with
    // 45 zones along a side - we need to scale it
    const Real_t ebase = Real_t(3.948746e+7);
    Real_t scale = (nx * m_tp) / Real_t(45.0);
    Real_t einit = ebase * scale * scale * scale;
    if (m_rowLoc + m_colLoc + m_planeLoc == 0)
    {
      // Dump into the first zone (which we know is in the corner)
      // of the domain that sits at the origin
      e(0) = einit;
    }
    //set initial deltatime base on analytic CFL calculation
    deltatime() = (Real_t(.5) * cbrt(volo(0))) / sqrt(Real_t(2.0) * einit);

  } // End constructor


Domain::~Domain()
{
#ifdef AXOM_USE_MPI
    delete [] commDataSend;
    delete [] commDataRecv;
#endif
}

////////////////////////////////////////////////////////////////////////////////
  void
  Domain::BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems)
  {
    Index_t meshEdgeElems = m_tp * nx;

    // initialize nodal coordinates
    Index_t nidx = 0;
    Real_t tz = Real_t(1.125) * Real_t(m_planeLoc * nx) / Real_t(meshEdgeElems);

    for (Index_t plane = 0; plane<edgeNodes; ++plane)
    {
      Real_t ty = Real_t(1.125) * Real_t(m_rowLoc * nx) / Real_t(meshEdgeElems);
      for (Index_t row = 0; row<edgeNodes; ++row)
      {
        Real_t tx = Real_t(1.125) * Real_t(m_colLoc * nx) / Real_t(meshEdgeElems);
        for (Index_t col = 0; col<edgeNodes; ++col)
        {
          x(nidx) = tx;
          y(nidx) = ty;
          z(nidx) = tz;
          ++nidx;
          // tx += ds ; // may accumulate roundoff...
          tx = Real_t(1.125) * Real_t(m_colLoc * nx + col + 1) / Real_t(meshEdgeElems);
        }
        // ty += ds ;  // may accumulate roundoff...
        ty = Real_t(1.125) * Real_t(m_rowLoc * nx + row + 1) / Real_t(meshEdgeElems);
      }
      // tz += ds ;  // may accumulate roundoff...
      tz = Real_t(1.125) * Real_t(m_planeLoc * nx + plane + 1) / Real_t(meshEdgeElems);
    }


    // embed hexahedral elements in nodal point lattice
    // SLAM NOTE: This could be represented using an implicit relation, when available
    // SLAM NOTE: Alternatively, the underlying connectivity should be derivable from an implicit Cartesian grid.
    IntsRegistry::BufferType& local_nodelist = m_intsRegistry.addBuffer("hex_node_boundary", 8 * numElem() );
    Index_t zidx = 0;
    nidx = 0;
    const int planeNodes = edgeNodes * edgeNodes;
    for (Index_t plane = 0; plane<edgeElems; ++plane)
    {
      for (Index_t row = 0; row<edgeElems; ++row)
      {
        for (Index_t col = 0; col<edgeElems; ++col)
        {
          Index_t *localNode = &local_nodelist[Index_t(8) * zidx];

          localNode[0] = nidx;
          localNode[1] = nidx                          + 1;
          localNode[2] = nidx              + edgeNodes + 1;
          localNode[3] = nidx              + edgeNodes;
          localNode[4] = nidx + planeNodes;
          localNode[5] = nidx + planeNodes             + 1;
          localNode[6] = nidx + planeNodes + edgeNodes + 1;
          localNode[7] = nidx + planeNodes + edgeNodes;
          ++zidx;
          ++nidx;
        }
        ++nidx;
      }
      nidx += edgeNodes;
    }

    m_nodelist.bindIndices(local_nodelist.size(), &local_nodelist);
    SLIC_ASSERT( m_nodelist.isValid());
  }


////////////////////////////////////////////////////////////////////////////////
  void
  Domain::SetupThreadSupportStructures()
  {
#ifdef AXOM_USE_OPENMP
    Index_t numthreads = omp_get_max_threads();
#else
    Index_t numthreads = 1;
#endif

    if (numthreads > 1)
    {
      NodeIndexMap nodeCornerCount(&m_nodeSet);         // Keeps track of the number of corners per node

      for (Index_t i = 0; i<numElem(); ++i)               // foreach elem
      {
        ElemNodeSet nl = nodelist(i);                     //   grab the Elem2Node relation for the element
        for (Index_t j = 0; j < nl.size(); ++j)           //   for each node of the element
        {
          ++(nodeCornerCount[ nl[j]] );                   //     increment the count for that node
        }
      }

      // Invert the zones to corner relation
      IntsRegistry::BufferType& local_begins = m_intsRegistry.addBuffer("node_corner_begins", numNode() + 1 );
      IntsRegistry::BufferType& local_relIndices = m_intsRegistry.addBuffer("node_corner_rel_indices", m_cornerSet.size() );

      local_begins[0] = 0;                                    // use the counts array to set the begins array
      for (Index_t i = 1; i <= numNode(); ++i)
      {
        local_begins[i] = local_begins[i - 1] + nodeCornerCount[i - 1];
        nodeCornerCount[i - 1] = 0;
      }

      for (Index_t i = 0; i < numElem(); ++i)                    // foreach elem
      {
        ElemNodeSet nl = nodelist(i);                            //   grab the Elem2Node relation for the elem
        for (Index_t j = 0; j < nl.size(); ++j)                  //   foreach node of the elem
        {
          Index_t m = nl[j];                                     //     m is the node index pointed to by this corner
          Index_t k = i * 8 + j;                                 //     k is the corner index (elem*8+offset)
          Index_t offset = local_begins[m] + nodeCornerCount[m]; //     offset is where this element belongs in the offsets array
          local_relIndices[offset] = k;                             //     this is the offsets array of the node to corner relation
          ++(nodeCornerCount[m]);                                //     increment the count for this node
        }
      }

      // Finally create the relation over these arrays and check validity
      m_nodeCornerRelation = NodeToCornerRelation(&m_nodeSet, &m_cornerSet);
      m_nodeCornerRelation.bindBeginOffsets(numNode(), &local_begins);
      m_nodeCornerRelation.bindIndices(local_relIndices.size(), &local_relIndices);

      SLIC_ASSERT_MSG(m_nodeCornerRelation.isValid(), "Generating Node to Corner relation: Corner index out of range." );
    }
  }


////////////////////////////////////////////////////////////////////////////////
  void
  Domain::SetupCommBuffers()
  {
    // allocate a buffer large enough for nodal ghost data
    Index_t maxEdgeSize = MAX(this->sizeX(), MAX(this->sizeY(), this->sizeZ())) + 1;

    m_maxPlaneSize = CACHE_ALIGN_REAL(maxEdgeSize * maxEdgeSize);
    m_maxEdgeSize = CACHE_ALIGN_REAL(maxEdgeSize);

    // assume communication to 6 neighbors by default
    m_rowMin   = (m_rowLoc == 0)          ? 0 : 1;
    m_rowMax   = (m_rowLoc == m_tp - 1)   ? 0 : 1;
    m_colMin   = (m_colLoc == 0)          ? 0 : 1;
    m_colMax   = (m_colLoc == m_tp - 1)   ? 0 : 1;
    m_planeMin = (m_planeLoc == 0)        ? 0 : 1;
    m_planeMax = (m_planeLoc == m_tp - 1) ? 0 : 1;

#ifdef AXOM_USE_MPI
    // account for face communication
    Index_t comBufSize =
        (m_rowMin + m_rowMax + m_colMin + m_colMax + m_planeMin + m_planeMax) *
        m_maxPlaneSize * MAX_FIELDS_PER_MPI_COMM;

    // account for edge communication
    comBufSize +=
        ((m_rowMin & m_colMin) + (m_rowMin & m_planeMin) + (m_colMin & m_planeMin) +
        (m_rowMax & m_colMax) + (m_rowMax & m_planeMax) + (m_colMax & m_planeMax) +
        (m_rowMax & m_colMin) + (m_rowMin & m_planeMax) + (m_colMin & m_planeMax) +
        (m_rowMin & m_colMax) + (m_rowMax & m_planeMin) + (m_colMax & m_planeMin)) *
        m_maxEdgeSize * MAX_FIELDS_PER_MPI_COMM;

    // account for corner communication
    // factor of 16 is so each buffer has its own cache line
    comBufSize += ((m_rowMin & m_colMin & m_planeMin) +
        (m_rowMin & m_colMin & m_planeMax) +
        (m_rowMin & m_colMax & m_planeMin) +
        (m_rowMin & m_colMax & m_planeMax) +
        (m_rowMax & m_colMin & m_planeMin) +
        (m_rowMax & m_colMin & m_planeMax) +
        (m_rowMax & m_colMax & m_planeMin) +
        (m_rowMax & m_colMax & m_planeMax)) * CACHE_COHERENCE_PAD_REAL;

    this->commDataSend = new Real_t[comBufSize];
    this->commDataRecv = new Real_t[comBufSize];
    // prevent floating point exceptions
    memset( this->commDataSend, 0,  comBufSize * sizeof(Real_t));
    memset( this->commDataRecv, 0,  comBufSize * sizeof(Real_t));
#endif

  }


////////////////////////////////////////////////////////////////////////////////
  void
  Domain::CreateRegionIndexSets(Int_t nr, Int_t balance)
  {
    using RegionToElemDynamicRelation = axom::slam::DynamicVariableRelation<RegionSet, ElemSet>;

#ifdef AXOM_USE_MPI
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    srand(myRank);
#else
    srand(0);
    Index_t myRank = 0;
#endif

    // Create a region set over the number of regions in the mesh -- regions are offset by one
    m_regionSet = RegionSet::SetBuilder().size(nr).offset(1);

    // Generate the region to element map as a dynamic relation.
    // We will linearize this into a static relation below...
    RegionToElemDynamicRelation reg2Elems(&m_regionSet, &m_elemSet);

    //if we only have one region just fill it
    if(numReg() == 1)
    {
      // Create the region material number map over the elements
      const Index_t regMatID = m_regionSet[0];
      m_elemRegNum = ElemIntMap(&m_elemSet, regMatID);
    }
    //If we have more than one region distribute the elements.
    else {
      Int_t regionNum, regionIdx;
      Int_t lastReg = -1;
      Index_t elements;
      Index_t runto = 0;

      // Create the region material map over the elements
      m_elemRegNum = ElemIntMap(&m_elemSet);

      // Setup an unnormalized CDF to randomly choose the region
      RegionIntMap regChooserProbs(&m_regionSet);

      // Determine the relative weights of all the regions as a CDF.
      // This is based off the -b flag.  Balance is the value passed into b.
      Int_t costDenominator = 0;
      for (Index_t i = 0; i< numReg(); ++i)
      {
        costDenominator += pow((i + 1), balance); //Total sum of all regions weights
        regChooserProbs[i] = costDenominator;   //Chance of hitting a given region is (regChooserProbs[i] - regChooserProbs[i-1])/costDenominator
      }

      // Setup an unnormalized CDF to randomly choose the number of elements in each run
      const int sizeChooserProbs[7] = {773, 937, 970, 974, 978, 981, 1000};
      using Builder = axom::slam::RangeSet<>::SetBuilder;
      axom::slam::RangeSet<> sizeChooser[7];
      sizeChooser[0] = Builder().range(1,16);
      sizeChooser[1] = Builder().range(16,32);
      sizeChooser[2] = Builder().range(32,64);
      sizeChooser[3] = Builder().range(64,128);
      sizeChooser[4] = Builder().range(128,256);
      sizeChooser[5] = Builder().range(256,512);
      sizeChooser[6] = Builder().range(512,2048);


      // Assign an element to each region
      Index_t nextIndex = 0;
      while (nextIndex < numElem())
      {
        //pick the region; make sure we don't pick the same region twice in a row
        do {
          Int_t regionVar = rand() % costDenominator;
          Index_t i = 0;
          while(regionVar >= regChooserProbs[i])
            i++;

          // Rotate regions based on MPI rank to give each domain a different region with the highest representation
          regionIdx = (i + myRank) % numReg();
          regionNum = m_regionSet[ regionIdx ];
        } while(regionNum == lastReg);

        // Determine the number of elements in this run using the sizeChooserProbs CDF and sizeChooser sets.
        Int_t binSize = rand() % 1000;
        Index_t binIndex = 0;
        while(binSize >= sizeChooserProbs[binIndex])
          ++binIndex;

        const axom::slam::RangeSet<>& rset = sizeChooser[ binIndex ];
        elements = rset[ rand() % rset.size() ]; // gets a random number from this bin's sizeChooser set

        runto = elements + nextIndex;

        //Store the elements.  If we hit the end before we run out of elements then just stop.
        for(; nextIndex < runto && nextIndex < numElem(); nextIndex++)
        {
          reg2Elems.insert(regionIdx, nextIndex);
          m_elemRegNum[nextIndex] = regionNum;
        }
        lastReg = regionNum;
      }
    }

    SLIC_ASSERT(reg2Elems.isValid());    // Ensure that the dynamic relation is valid

    // Convert from a Dynamic to a Static relation
    IntsRegistry::BufferType& local_begins = m_intsRegistry.addBuffer("reg_elem_begins", numReg() + 1 );
    IntsRegistry::BufferType& local_relIndices = m_intsRegistry.addBuffer("reg_elem_rel_indices", numElem() );
    Index_t curOffIdx = 0;
    for(Index_t regionPos = 0; regionPos < numReg(); ++regionPos)
    {
      local_begins[ regionPos] = curOffIdx;
      for(Index_t elemRelPos = 0; elemRelPos < reg2Elems.size( regionPos); ++elemRelPos)
      {
        local_relIndices[curOffIdx++] = reg2Elems[ regionPos][elemRelPos];
      }
    }
    local_begins[ numReg()] = local_relIndices.size();

    /// We can finally create the region to elem relation
    m_regionElementsRel = RegionToElemRelation(&m_regionSet, &m_elemSet);
    m_regionElementsRel.bindBeginOffsets(numReg(), &local_begins);
    m_regionElementsRel.bindIndices(local_relIndices.size(), &local_relIndices);

    SLIC_ASSERT(m_regionElementsRel.isValid(true));    // Ensure that the relation is valid
  }

/////////////////////////////////////////////////////////////
  void
  Domain::SetupSymmetryPlanes(Int_t edgeNodes)
  {
    Int_t numSymmNodesX = m_colLoc == 0   ? edgeNodes * edgeNodes : 0;
    Int_t numSymmNodesY = m_rowLoc == 0   ? edgeNodes * edgeNodes : 0;
    Int_t numSymmNodesZ = m_planeLoc == 0 ? edgeNodes * edgeNodes : 0;

    // Setup symmetry NodeSets            // TODO: Need to add the mesh's nodes as parent
    m_symmX = SymmNodeSet(numSymmNodesX); // eparentSet = &m_nodeSet
    m_symmY = SymmNodeSet(numSymmNodesY);
    m_symmZ = SymmNodeSet(numSymmNodesZ);

    auto& loc_symmX = m_intsRegistry.addField("symmX", &m_symmX).data();
    auto& loc_symmY = m_intsRegistry.addField("symmY", &m_symmY).data();
    auto& loc_symmZ = m_intsRegistry.addField("symmZ", &m_symmZ).data();

    // SLAM Note: We should be able to compute these directory from a Cartesian product set defining a regular grid.
    Index_t nidx = 0;
    for (Index_t i = 0; i<edgeNodes; ++i)
    {
      Index_t planeInc = i * edgeNodes * edgeNodes;
      Index_t rowInc   = i * edgeNodes;
      for (Index_t j = 0; j<edgeNodes; ++j)
      {
        if (m_planeLoc == 0)
        {
          loc_symmZ[nidx] = rowInc   + j;
        }
        if (m_rowLoc == 0)
        {
          loc_symmY[nidx] = planeInc + j;
        }
        if (m_colLoc == 0)
        {
          loc_symmX[nidx] = planeInc + j * edgeNodes;
        }
        ++nidx;
      }
    }

    m_symmX.ptr() = &loc_symmX;
    m_symmY.ptr() = &loc_symmY;
    m_symmZ.ptr() = &loc_symmZ;

    // Verify validity of the sets.
    SLIC_ASSERT(  m_symmX.isValid() && m_symmX.size() == numSymmNodesX && m_symmX.ptr() == &loc_symmX);
    SLIC_ASSERT(  m_symmY.isValid() && m_symmY.size() == numSymmNodesY && m_symmY.ptr() == &loc_symmY);
    SLIC_ASSERT(  m_symmZ.isValid() && m_symmZ.size() == numSymmNodesZ && m_symmZ.ptr() == &loc_symmZ);
  }



/////////////////////////////////////////////////////////////
  void
  Domain::SetupElementConnectivities(Int_t edgeElems)
  {
    // Create temporary arrays to hold the data
    IntsRegistry::BufferType& indices_xim = m_intsRegistry.addBuffer("zone_face_xi_m", numElem());
    IntsRegistry::BufferType& indices_xip = m_intsRegistry.addBuffer("zone_face_xi_p", numElem());

    // Setup xi face adjacencies
    indices_xim[0] = 0;
    for (Index_t i = 1; i<numElem(); ++i)
    {
      indices_xim[i]   = i - 1;
      indices_xip[i - 1] = i;
    }
    indices_xip[numElem() - 1] = numElem() - 1;
    m_lxim.bindIndices(numElem(), &indices_xim);
    m_lxip.bindIndices(numElem(), &indices_xip);

    // Setup eta face adjacencies
    IntsRegistry::BufferType& indices_etam = m_intsRegistry.addBuffer("zone_face_eta_m", numElem());
    IntsRegistry::BufferType& indices_etap = m_intsRegistry.addBuffer("zone_face_eta_p", numElem());
    for (Index_t i = 0; i<edgeElems; ++i)
    {
      indices_etam[i] = i;
      indices_etap[numElem() - edgeElems + i] = numElem() - edgeElems + i;
    }
    for (Index_t i = edgeElems; i<numElem(); ++i)
    {
      indices_etam[i] = i - edgeElems;
      indices_etap[i - edgeElems] = i;
    }
    m_letam.bindIndices( numElem(), &indices_etam);
    m_letap.bindIndices( numElem(), &indices_etap);

    // Setup zeta face adjacencies
    IntsRegistry::BufferType& indices_zetam = m_intsRegistry.addBuffer("zone_face_zeta_m", numElem());
    IntsRegistry::BufferType& indices_zetap = m_intsRegistry.addBuffer("zone_face_zeta_p", numElem());
    for (Index_t i = 0; i<edgeElems * edgeElems; ++i)
    {
      indices_zetam[i] = i;
      indices_zetap[numElem() - edgeElems * edgeElems + i] = numElem() - edgeElems * edgeElems + i;
    }
    for (Index_t i = edgeElems * edgeElems; i<numElem(); ++i)
    {
      indices_zetam[i] = i - edgeElems * edgeElems;
      indices_zetap[i - edgeElems * edgeElems] = i;
    }
    m_lzetam.bindIndices(numElem(), &indices_zetam);
    m_lzetap.bindIndices(numElem(), &indices_zetap);

    // Ensure that all the indices in the relations are valid
    SLIC_ASSERT(  m_lxim.isValid() );
    SLIC_ASSERT(  m_lxip.isValid() );
    SLIC_ASSERT(  m_letam.isValid() );
    SLIC_ASSERT(  m_letap.isValid() );
    SLIC_ASSERT(  m_lzetam.isValid() );
    SLIC_ASSERT(  m_lzetap.isValid() );
  }

/////////////////////////////////////////////////////////////
  void
  Domain::SetupBoundaryConditions(Int_t edgeElems)
  {
    Index_t ghostIdx[6]; // offsets to ghost locations

    // set up boundary condition information
    for (Index_t i = 0; i<numElem(); ++i)
    {
      elemBC(i) = Int_t(0);
    }

    for (Index_t i = 0; i<6; ++i)
    {
      ghostIdx[i] = INT_MIN;
    }

    Int_t pidx = numElem();
    if (m_planeMin != 0)
    {
      ghostIdx[0] = pidx;
      pidx += sizeX() * sizeY();
    }

    if (m_planeMax != 0)
    {
      ghostIdx[1] = pidx;
      pidx += sizeX() * sizeY();
    }

    if (m_rowMin != 0)
    {
      ghostIdx[2] = pidx;
      pidx += sizeX() * sizeZ();
    }

    if (m_rowMax != 0)
    {
      ghostIdx[3] = pidx;
      pidx += sizeX() * sizeZ();
    }

    if (m_colMin != 0)
    {
      ghostIdx[4] = pidx;
      pidx += sizeY() * sizeZ();
    }

    if (m_colMax != 0)
    {
      ghostIdx[5] = pidx;
    }


    // SLAM HACK: We want to modify the data of a StaticRelation (ElemFaceAdjacencyRelation).
    //            This will be addressed in ATK-927.
    //            For now, grab the underlying array from the registry

    IntsRegistry::BufferType& local_xi_m = m_intsRegistry.getBuffer("zone_face_xi_m");
    IntsRegistry::BufferType& local_xi_p = m_intsRegistry.getBuffer("zone_face_xi_p");
    IntsRegistry::BufferType& local_eta_m = m_intsRegistry.getBuffer("zone_face_eta_m");
    IntsRegistry::BufferType& local_eta_p = m_intsRegistry.getBuffer("zone_face_eta_p");
    IntsRegistry::BufferType& local_zeta_m = m_intsRegistry.getBuffer("zone_face_zeta_m");
    IntsRegistry::BufferType& local_zeta_p = m_intsRegistry.getBuffer("zone_face_zeta_p");

    // symmetry plane or free surface BCs
    for (Index_t i = 0; i<edgeElems; ++i)
    {
      Index_t planeInc = i * edgeElems * edgeElems;
      Index_t rowInc   = i * edgeElems;
      for (Index_t j = 0; j<edgeElems; ++j)
      {
        if (m_planeLoc == 0)
        {
          elemBC(rowInc + j) |= ZETA_M_SYMM;
        }
        else {
          elemBC(rowInc + j) |= ZETA_M_COMM;
          local_zeta_m[rowInc + j] = ghostIdx[0] + rowInc + j;
        }

        if (m_planeLoc == m_tp - 1)
        {
          elemBC(rowInc + j + numElem() - edgeElems * edgeElems) |= ZETA_P_FREE;
        }
        else {
          elemBC(rowInc + j + numElem() - edgeElems * edgeElems) |= ZETA_P_COMM;
          local_zeta_p[rowInc + j + numElem() - edgeElems * edgeElems] = ghostIdx[1] + rowInc + j;
        }

        if (m_rowLoc == 0)
        {
          elemBC(planeInc + j) |= ETA_M_SYMM;
        }
        else {
          elemBC(planeInc + j) |= ETA_M_COMM;
          local_eta_m[planeInc + j] = ghostIdx[2] + rowInc + j;
        }

        if (m_rowLoc == m_tp - 1)
        {
          elemBC(planeInc + j + edgeElems * edgeElems - edgeElems) |= ETA_P_FREE;
        }
        else {
          elemBC(planeInc + j + edgeElems * edgeElems - edgeElems) |= ETA_P_COMM;
          local_eta_p[planeInc + j + edgeElems * edgeElems - edgeElems] = ghostIdx[3] +  rowInc + j;
        }

        if (m_colLoc == 0)
        {
          elemBC(planeInc + j * edgeElems) |= XI_M_SYMM;
        }
        else {
          elemBC(planeInc + j * edgeElems) |= XI_M_COMM;
          local_xi_m[planeInc + j * edgeElems] = ghostIdx[4] + rowInc + j;
        }

        if (m_colLoc == m_tp - 1)
        {
          elemBC(planeInc + j * edgeElems + edgeElems - 1) |= XI_P_FREE;
        }
        else {
          elemBC(planeInc + j * edgeElems + edgeElems - 1) |= XI_P_COMM;
          local_xi_p[planeInc + j * edgeElems + edgeElems - 1] = ghostIdx[5] + rowInc + j;
        }
      }
    }

    // Ensure that all the indices in the element adjacency relations are still valid
    SLIC_ASSERT(  m_lxim.isValid() );
    SLIC_ASSERT(  m_lxip.isValid() );
    SLIC_ASSERT(  m_letam.isValid() );
    SLIC_ASSERT(  m_letap.isValid() );
    SLIC_ASSERT(  m_lzetam.isValid() );
    SLIC_ASSERT(  m_lzetap.isValid() );
  }

///////////////////////////////////////////////////////////////////////////
  void InitMeshDecomp(Int_t numRanks, Int_t myRank,
      Int_t *col, Int_t *row, Int_t *plane, Int_t *side)
  {
    Int_t testProcs;
    Int_t dx, dy, dz;
    Int_t myDom;

    // Assume cube processor layout for now
    testProcs = Int_t(cbrt(Real_t(numRanks)) + 0.5);
    if (testProcs * testProcs * testProcs != numRanks)
    {
      SLIC_WARNING("Num processors must be a cube of an integer (1, 8, 27, ...)");
#ifdef AXOM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
    }
    if (sizeof(Real_t) != 4 && sizeof(Real_t) != 8)
    {
      SLIC_WARNING("MPI operations only support float and double right now...");
#ifdef AXOM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
    }
    if (MAX_FIELDS_PER_MPI_COMM > CACHE_COHERENCE_PAD_REAL)
    {
      SLIC_WARNING("Corner element comm buffers too small.  Fix code.");
#ifdef AXOM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
    }

    dx = testProcs;
    dy = testProcs;
    dz = testProcs;

    // temporary test
    if (dx * dy * dz != numRanks)
    {
      SLIC_WARNING("Error -- must have as many domains as procs.");
#ifdef AXOM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, -1);
#else
      exit(-1);
#endif
    }
    Int_t remainder = dx * dy * dz % numRanks;
    if (myRank < remainder)
    {
      myDom = myRank * ( 1 + (dx * dy * dz / numRanks));
    }
    else {
      myDom = remainder * ( 1 + (dx * dy * dz / numRanks)) +
          (myRank - remainder) * (dx * dy * dz / numRanks);
    }

    *col = myDom % dx;
    *row = (myDom / dx) % dy;
    *plane = myDom / (dx * dy);
    *side = testProcs;

    return;
  }

} // end namespace slamLulesh
