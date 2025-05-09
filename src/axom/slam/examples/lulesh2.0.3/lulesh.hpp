// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic.hpp"
#include "axom/slam.hpp"

#ifdef AXOM_USE_MPI
#include <mpi.h>

/*
   define one of these three symbols:

   SEDOV_SYNC_POS_VEL_NONE
   SEDOV_SYNC_POS_VEL_EARLY
   SEDOV_SYNC_POS_VEL_LATE
 */

#define SEDOV_SYNC_POS_VEL_EARLY 1
#endif // AXOM_USE_MPI

#include <math.h>
#include <vector>

namespace slamLulesh {

//**************************************************
// Allow flexibility for arithmetic representations
//**************************************************

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))


// Precision specification
  typedef float       real4;
  typedef double      real8;
  typedef long double real10;  // 10 bytes on x86

  typedef axom::IndexType  Index_t; // array subscript and loop index
  typedef real8            Real_t; // floating point representation
  typedef axom::IndexType  Int_t; // integer representation

  enum { VolumeError = -1, QStopError = -2 };

  inline real4  SQRT(real4 arg) { return sqrtf(arg); }
  inline real8  SQRT(real8 arg) { return sqrt(arg); }
  inline real10 SQRT(real10 arg) { return sqrtl(arg); }

  inline real4  CBRT(real4 arg) { return cbrtf(arg); }
  inline real8  CBRT(real8 arg) { return cbrt(arg); }
  inline real10 CBRT(real10 arg) { return cbrtl(arg); }

  inline real4  FABS(real4 arg) { return fabsf(arg); }
  inline real8  FABS(real8 arg) { return fabs(arg); }
  inline real10 FABS(real10 arg) { return fabsl(arg); }


// Stuff needed for boundary conditions
// 3 BCs on each of 6 hexahedral faces (18 bits)
  enum { XI_M_SYMM   = 1 << 0,
         XI_M_FREE   = 1 << 1,
         XI_M_COMM   = 1 << 2,
         XI_M        = XI_M_SYMM | XI_M_FREE | XI_M_COMM,

         XI_P_SYMM   = 1 << 3,
         XI_P_FREE   = 1 << 4,
         XI_P_COMM   = 1 << 5,
         XI_P        = XI_P_SYMM | XI_P_FREE | XI_P_COMM,

         ETA_M_SYMM  = 1 << 6,
         ETA_M_FREE  = 1 << 7,
         ETA_M_COMM  = 1 << 8,
         ETA_M       = ETA_M_SYMM | ETA_M_FREE | ETA_M_COMM,

         ETA_P_SYMM  = 1 << 9,
         ETA_P_FREE  = 1 << 10,
         ETA_P_COMM  = 1 << 11,
         ETA_P       = ETA_P_SYMM | ETA_P_FREE | ETA_P_COMM,

         ZETA_M_SYMM = 1 << 12,
         ZETA_M_FREE = 1 << 13,
         ZETA_M_COMM = 1 << 14,
         ZETA_M      = ZETA_M_SYMM | ZETA_M_FREE | ZETA_M_COMM,

         ZETA_P_SYMM = 1 << 15,
         ZETA_P_FREE = 1 << 16,
         ZETA_P_COMM = 1 << 17,
         ZETA_P      = ZETA_P_SYMM | ZETA_P_FREE | ZETA_P_COMM};

// MPI Message Tags
#define MSG_COMM_SBN      1024
#define MSG_SYNC_POS_VEL  2048
#define MSG_MONOQ         3072

#define MAX_FIELDS_PER_MPI_COMM 6

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n) \
  (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL - 1))

#ifdef AXOM_USE_MPI
  struct CommMessageMetadata;
#endif

//////////////////////////////////////////////////////
// Primary data structure
//////////////////////////////////////////////////////

/*
 * The implementation of the data abstraction used for lulesh
 * resides entirely in the Domain class below.  You can change
 * grouping and interleaving of fields here to maximize data layout
 * efficiency for your underlying architecture or compiler.
 *
 * For example, fields can be implemented as STL objects or
 * raw array pointers.  As another example, individual fields
 * m_x, m_y, m_z could be bundled into
 *
 *    struct { Real_t x, y, z ; } *m_coord ;
 *
 * allowing accessor functions such as
 *
 *  "Real_t &x(Index_t idx) { return m_coord[idx].x ; }"
 *  "Real_t &y(Index_t idx) { return m_coord[idx].y ; }"
 *  "Real_t &z(Index_t idx) { return m_coord[idx].z ; }"
 */

  class Domain {

  public:
    using SetBase = axom::slam::Set<>;
    using NullSet = axom::slam::NullSet<>;
    using PositionType = SetBase::PositionType;
    using ElementType = SetBase::ElementType;

    using ElemSet = axom::slam::RangeSet<>;
    using NodeSet = axom::slam::RangeSet<>;
    using CornerSet = axom::slam::RangeSet<>;
    using ExtendedElemSet = axom::slam::RangeSet<>;

    using SymmNodeSet = axom::slam::VectorIndirectionSet<>;

    using STLIndirection = axom::slam::policies::STLVectorIndirection<PositionType, ElementType>;

    enum { NODES_PER_ZONE = 8, FACES_PER_ZONE = 1};
    using ZNStride = axom::slam::policies::CompileTimeStride<PositionType, NODES_PER_ZONE>;
    using ZFStride = axom::slam::policies::CompileTimeStride<PositionType, FACES_PER_ZONE>;

    using ZNCard = axom::slam::policies::ConstantCardinality<PositionType, ZNStride>;
    using ElemToNodeRelation =
        axom::slam::StaticRelation<PositionType, ElementType,ZNCard, STLIndirection, ElemSet, NodeSet>;
    using ElemNodeSet = const ElemToNodeRelation::RelationSubset;

    using ZFCard = axom::slam::policies::ConstantCardinality<PositionType, ZFStride>;
    using ElemFaceAdjacencyRelation =
        axom::slam::StaticRelation<PositionType, ElementType,ZFCard, STLIndirection, ElemSet, ExtendedElemSet>;


    using RegionSet = axom::slam::RangeSet<>;
    using VariableCardinality = axom::slam::policies::VariableCardinality<PositionType, STLIndirection>;
    using RegionToElemRelation = axom::slam::StaticRelation<
                                    PositionType, ElementType,
                                    VariableCardinality,
                                    STLIndirection,
                                    RegionSet,
                                    ElemSet>;
    using RegionElemSet = const RegionToElemRelation::RelationSubset;


    using NodeToCornerRelation = axom::slam::StaticRelation<
                                    PositionType, ElementType,
                                    VariableCardinality,
                                    STLIndirection,
                                    NodeSet,
                                    CornerSet>;
    using NodeCornerSet = const NodeToCornerRelation::RelationSubset ;

    using ElemIndexMap = axom::slam::Map<Index_t>;
    using ElemIntMap = axom::slam::Map<Int_t>;
    //using ElemRealMap = axom::slam::Map<Real_t>;

    using NodeIndexMap = axom::slam::Map<Index_t>;
    //using NodeIntMap = axom::slam::Map<Int_t>;
    //using NodeRealMap = axom::slam::Map<Real_t>;

    //using RegionIndexMap = axom::slam::Map<Index_t>;
    using RegionIntMap = axom::slam::Map<Int_t>;
    //using RegionRealMap = axom::slam::Map<Real_t>;

    //using CornerIndexMap = axom::slam::Map<Index_t>;
    //using CornerIntMap = axom::slam::Map<Int_t>;
    using CornerRealMap = axom::slam::Map<Real_t>;

    using RealsRegistry = axom::slam::FieldRegistry<SetBase, Real_t>;
    using IntsRegistry = axom::slam::FieldRegistry<SetBase, Index_t>;

  public:

    // Constructor
    Domain(Int_t numRanks, Index_t colLoc,
        Index_t rowLoc, Index_t planeLoc,
        Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost);


    // Destructor
    ~Domain();
    
    //
    // ALLOCATION
    //

    void AllocateNodePersistent(Int_t numNode) // Node-centered
    {
      m_x.resize(numNode);  // coordinates
      m_y.resize(numNode);
      m_z.resize(numNode);

      m_xd.resize(numNode); // velocities
      m_yd.resize(numNode);
      m_zd.resize(numNode);

      m_xdd.resize(numNode); // accelerations
      m_ydd.resize(numNode);
      m_zdd.resize(numNode);

      m_fx.resize(numNode);  // forces
      m_fy.resize(numNode);
      m_fz.resize(numNode);

      m_nodalMass.resize(numNode);  // mass
    }

    void AllocateElemPersistent(Int_t numElem) // Elem-centered
    {
      // elem to node incidence relation
      m_nodelist = ElemToNodeRelation(&m_elemSet, &m_nodeSet);

      // elem adjacencies through face
      m_lxim   = ElemFaceAdjacencyRelation(&m_elemSet,&m_extendedElemSet);
      m_lxip   = ElemFaceAdjacencyRelation(&m_elemSet,&m_extendedElemSet);
      m_letam  = ElemFaceAdjacencyRelation(&m_elemSet,&m_extendedElemSet);
      m_letap  = ElemFaceAdjacencyRelation(&m_elemSet,&m_extendedElemSet);
      m_lzetam = ElemFaceAdjacencyRelation(&m_elemSet,&m_extendedElemSet);
      m_lzetap = ElemFaceAdjacencyRelation(&m_elemSet,&m_extendedElemSet);

      m_elemBC = ElemIntMap(&m_elemSet);

      m_e.resize(numElem);
      m_p.resize(numElem);

      m_q.resize(numElem);
      m_ql.resize(numElem);
      m_qq.resize(numElem);

      m_v.resize(numElem);

      m_volo.resize(numElem);
      m_delv.resize(numElem);
      m_vdov.resize(numElem);

      m_arealg.resize(numElem);

      m_ss.resize(numElem);

      m_elemMass.resize(numElem);
    }

    void AllocateGradients(Int_t numElem, Int_t allElem)
    {
      // Position gradients
      m_delx_xi.resize(numElem);
      m_delx_eta.resize(numElem);
      m_delx_zeta.resize(numElem);

      // Velocity gradients
      m_delv_xi.resize(allElem);
      m_delv_eta.resize(allElem);
      m_delv_zeta.resize(allElem);
    }

    void DeallocateGradients()
    {
      m_delx_zeta.clear();
      m_delx_eta.clear();
      m_delx_xi.clear();

      m_delv_zeta.clear();
      m_delv_eta.clear();
      m_delv_xi.clear();
    }

    void AllocateStrains(Int_t numElem)
    {
      m_dxx.resize(numElem);
      m_dyy.resize(numElem);
      m_dzz.resize(numElem);
    }

    void DeallocateStrains()
    {
      m_dzz.clear();
      m_dyy.clear();
      m_dxx.clear();
    }

    //
    // ACCESSORS
    //

    // Node-centered

    // Nodal coordinates
    Real_t& x(Index_t idx)    { return m_x[idx]; }
    Real_t& y(Index_t idx)    { return m_y[idx]; }
    Real_t& z(Index_t idx)    { return m_z[idx]; }

    // Nodal velocities
    Real_t& xd(Index_t idx)   { return m_xd[idx]; }
    Real_t& yd(Index_t idx)   { return m_yd[idx]; }
    Real_t& zd(Index_t idx)   { return m_zd[idx]; }

    // Nodal accelerations
    Real_t& xdd(Index_t idx)  { return m_xdd[idx]; }
    Real_t& ydd(Index_t idx)  { return m_ydd[idx]; }
    Real_t& zdd(Index_t idx)  { return m_zdd[idx]; }

    // Nodal forces
    Real_t& fx(Index_t idx)   { return m_fx[idx]; }
    Real_t& fy(Index_t idx)   { return m_fy[idx]; }
    Real_t& fz(Index_t idx)   { return m_fz[idx]; }

    // Nodal mass
    Real_t& nodalMass(Index_t idx) { return m_nodalMass[idx]; }

    // Nodes on symmetry planes
    Index_t symmX(Index_t idx) { return m_symmX[idx]; }
    Index_t symmY(Index_t idx) { return m_symmY[idx]; }
    Index_t symmZ(Index_t idx) { return m_symmZ[idx]; }
    bool    symmXempty()       { return m_symmX.empty(); }
    bool    symmYempty()       { return m_symmY.empty(); }
    bool    symmZempty()       { return m_symmZ.empty(); }

    //
    // Element-centered
    //
    Index_t         regElemSize(Index_t idx) { return m_regionElementsRel.size(idx); }
    Index_t&        regNumList(Index_t idx) { return m_elemRegNum[idx]; }
    Index_t*        regNumList()            { return &m_elemRegNum[0]; }
    RegionElemSet   regElemlist(Int_t r)    { return m_regionElementsRel[r]; }
    Index_t         regElemlist(Int_t r, Index_t idx) { return m_regionElementsRel[r][idx]; }

    ElemNodeSet     nodelist(Index_t idx)    { return m_nodelist[idx]; }

    // elem adjacencies through face
    const Index_t&  lxim(Index_t idx) { return m_lxim[idx][0]; }
    const Index_t&  lxip(Index_t idx) { return m_lxip[idx][0]; }
    const Index_t&  letam(Index_t idx) { return m_letam[idx][0]; }
    const Index_t&  letap(Index_t idx) { return m_letap[idx][0]; }
    const Index_t&  lzetam(Index_t idx) { return m_lzetam[idx][0]; }
    const Index_t&  lzetap(Index_t idx) { return m_lzetap[idx][0]; }

    // elem face symm/free-surface flag
    Int_t&          elemBC(Index_t idx) { return m_elemBC[idx]; }

    // Principal strains - temporary
    Real_t&         dxx(Index_t idx)  { return m_dxx[idx]; }
    Real_t&         dyy(Index_t idx)  { return m_dyy[idx]; }
    Real_t&         dzz(Index_t idx)  { return m_dzz[idx]; }

    // Velocity gradient - temporary
    Real_t&         delv_xi(Index_t idx)    { return m_delv_xi[idx]; }
    Real_t&         delv_eta(Index_t idx)   { return m_delv_eta[idx]; }
    Real_t&         delv_zeta(Index_t idx)  { return m_delv_zeta[idx]; }

    // Position gradient - temporary
    Real_t&         delx_xi(Index_t idx)    { return m_delx_xi[idx]; }
    Real_t&         delx_eta(Index_t idx)   { return m_delx_eta[idx]; }
    Real_t&         delx_zeta(Index_t idx)  { return m_delx_zeta[idx]; }

    // Energy
    Real_t&         e(Index_t idx)          { return m_e[idx]; }

    // Pressure
    Real_t&         p(Index_t idx)          { return m_p[idx]; }

    // Artificial viscosity
    Real_t&         q(Index_t idx)          { return m_q[idx]; }

    // Linear term for q
    Real_t&         ql(Index_t idx)         { return m_ql[idx]; }
    // Quadratic term for q
    Real_t&         qq(Index_t idx)         { return m_qq[idx]; }

    // Relative volume
    Real_t&         v(Index_t idx)          { return m_v[idx]; }
    Real_t&         delv(Index_t idx)       { return m_delv[idx]; }

    // Reference volume
    Real_t&         volo(Index_t idx)       { return m_volo[idx]; }

    // volume derivative over volume
    Real_t&         vdov(Index_t idx)       { return m_vdov[idx]; }

    // Element characteristic length
    Real_t&         arealg(Index_t idx)     { return m_arealg[idx]; }

    // Sound speed
    Real_t&         ss(Index_t idx)         { return m_ss[idx]; }

    // Element mass
    Real_t&         elemMass(Index_t idx)  { return m_elemMass[idx]; }


    Index_t         nodeElemCount(Index_t idx) { return m_nodeCornerRelation.size(idx); }
    NodeCornerSet   nodeElemCornerList(Index_t idx) { return m_nodeCornerRelation[idx]; }


    /**
     * \brief Returns a const reference to the corner set when threading (omp) is enabled, otherwise, a ref to a NullSet (with no elements).
     */
    const SetBase&      threadingCornerSet()  const {
#ifdef AXOM_USE_OPENMP
      return m_cornerSet;
#else
      static const NullSet s_nullSet;
      return s_nullSet;
#endif
    }

    /**
     * Returns a const reference to the corner set.
     */
    const CornerSet& cornerSet()  const {return m_cornerSet;   }

    // Parameters

    // Cutoffs
    Real_t    u_cut() const { return m_u_cut; }
    Real_t    e_cut() const { return m_e_cut; }
    Real_t    p_cut() const { return m_p_cut; }
    Real_t    q_cut() const { return m_q_cut; }
    Real_t    v_cut() const { return m_v_cut; }

    // Other constants (usually are settable via input file in real codes)
    Real_t    hgcoef() const { return m_hgcoef; }
    Real_t    qstop() const { return m_qstop; }
    Real_t    monoq_max_slope() const { return m_monoq_max_slope; }
    Real_t    monoq_limiter_mult() const { return m_monoq_limiter_mult; }
    Real_t    ss4o3() const { return m_ss4o3; }
    Real_t    qlc_monoq() const { return m_qlc_monoq; }
    Real_t    qqc_monoq() const { return m_qqc_monoq; }
    Real_t    qqc() const { return m_qqc; }

    Real_t    eosvmax() const { return m_eosvmax; }
    Real_t    eosvmin() const { return m_eosvmin; }
    Real_t    pmin() const { return m_pmin; }
    Real_t    emin() const { return m_emin; }
    Real_t    dvovmax() const { return m_dvovmax; }
    Real_t    refdens() const { return m_refdens; }

    // Timestep controls, etc...
    Real_t&   time()                 { return m_time; }
    Real_t&   deltatime()            { return m_deltatime; }
    Real_t&   deltatimemultlb()      { return m_deltatimemultlb; }
    Real_t&   deltatimemultub()      { return m_deltatimemultub; }
    Real_t&   stoptime()             { return m_stoptime; }
    Real_t&   dtcourant()            { return m_dtcourant; }
    Real_t&   dthydro()              { return m_dthydro; }
    Real_t&   dtmax()                { return m_dtmax; }
    Real_t&   dtfixed()              { return m_dtfixed; }

    Int_t&    cycle()                { return m_cycle; }
    Index_t&  numRanks()             { return m_numRanks; }

    Index_t&  colLoc()               { return m_colLoc; }
    Index_t&  rowLoc()               { return m_rowLoc; }
    Index_t&  planeLoc()             { return m_planeLoc; }
    Index_t&  tp()                   { return m_tp; }

    Index_t&  sizeX()                { return m_sizeX; }
    Index_t&  sizeY()                { return m_sizeY; }
    Index_t&  sizeZ()                { return m_sizeZ; }

    Index_t   numReg()               { return m_regionSet.size(); }
    Int_t&    cost()                 { return m_cost; }

    Index_t   numElem()           const { return m_elemSet.size(); }
    Index_t   numNode()           const { return m_nodeSet.size(); }
    Index_t   numElemWithGhosts() const { return m_extendedElemSet.size(); }

    Index_t&  maxPlaneSize()         { return m_maxPlaneSize; }
    Index_t&  maxEdgeSize()          { return m_maxEdgeSize; }

    //
    // MPI-Related additional data
    //

#ifdef AXOM_USE_MPI
    // Communication Work space
    Real_t *commDataSend {nullptr};
    Real_t *commDataRecv {nullptr};

    // Maximum number of block neighbors
    MPI_Request recvRequest[26]; // 6 faces + 12 edges + 8 corners
    MPI_Request sendRequest[26]; // 6 faces + 12 edges + 8 corners
#endif

  private:

    void  BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems);
    void  SetupThreadSupportStructures();
    void  CreateRegionIndexSets(Int_t nreg, Int_t balance);
    void  SetupCommBuffers();
    void  SetupSymmetryPlanes(Int_t edgeNodes);
    void  SetupElementConnectivities(Int_t edgeElems);
    void  SetupBoundaryConditions(Int_t edgeElems);

    //
    // IMPLEMENTATION
    //

    /// Node-centered
    std::vector<Real_t> m_x;  /* coordinates */
    std::vector<Real_t> m_y;
    std::vector<Real_t> m_z;

    std::vector<Real_t> m_xd; /* velocities */
    std::vector<Real_t> m_yd;
    std::vector<Real_t> m_zd;

    std::vector<Real_t> m_xdd; /* accelerations */
    std::vector<Real_t> m_ydd;
    std::vector<Real_t> m_zdd;

    std::vector<Real_t> m_fx;  /* forces */
    std::vector<Real_t> m_fy;
    std::vector<Real_t> m_fz;

    std::vector<Real_t> m_nodalMass;  /* mass */

    SymmNodeSet m_symmX;  /* symmetry plane nodesets */
    SymmNodeSet m_symmY;
    SymmNodeSet m_symmZ;

    /// Element-centered

    // Region information
    //Int_t    m_numReg ;
    RegionSet m_regionSet;

    Int_t m_cost;   //imbalance cost

    ElemIntMap m_elemRegNum;                    // previously m_regNumList
    RegionToElemRelation m_regionElementsRel;   // previously m_regElemSize and m_regElemlist

    ElemToNodeRelation m_nodelist;              // elemToNode connectivity


    ElemFaceAdjacencyRelation m_lxim;  /* element connectivity across each face */
    ElemFaceAdjacencyRelation m_lxip;
    ElemFaceAdjacencyRelation m_letam;
    ElemFaceAdjacencyRelation m_letap;
    ElemFaceAdjacencyRelation m_lzetam;
    ElemFaceAdjacencyRelation m_lzetap;

    ElemIntMap m_elemBC;     /* symmetry/free-surface flags for each elem face */

    std::vector<Real_t> m_dxx;  /* principal strains -- temporary */
    std::vector<Real_t> m_dyy;
    std::vector<Real_t> m_dzz;

    std::vector<Real_t> m_delv_xi;    /* velocity gradient -- temporary */
    std::vector<Real_t> m_delv_eta;
    std::vector<Real_t> m_delv_zeta;

    std::vector<Real_t> m_delx_xi;    /* coordinate gradient -- temporary */
    std::vector<Real_t> m_delx_eta;
    std::vector<Real_t> m_delx_zeta;

    std::vector<Real_t> m_e;   /* energy */

    std::vector<Real_t> m_p;   /* pressure */
    std::vector<Real_t> m_q;   /* q */
    std::vector<Real_t> m_ql;  /* linear term for q */
    std::vector<Real_t> m_qq;  /* quadratic term for q */

    std::vector<Real_t> m_v;     /* relative volume */
    std::vector<Real_t> m_volo;  /* reference volume */
    std::vector<Real_t> m_vnew;  /* new relative volume -- temporary */
    std::vector<Real_t> m_delv;  /* m_vnew - m_v */
    std::vector<Real_t> m_vdov;  /* volume derivative over volume */

    std::vector<Real_t> m_arealg;  /* characteristic length of an element */

    std::vector<Real_t> m_ss;      /* "sound speed" */

    std::vector<Real_t> m_elemMass;  /* mass */

    // Cutoffs (treat as constants)
    const Real_t m_e_cut;              // energy tolerance
    const Real_t m_p_cut;              // pressure tolerance
    const Real_t m_q_cut;              // q tolerance
    const Real_t m_v_cut;              // relative volume tolerance
    const Real_t m_u_cut;              // velocity tolerance

    // Other constants (usually setable, but hardcoded in this proxy app)

    const Real_t m_hgcoef;             // hourglass control
    const Real_t m_ss4o3;
    const Real_t m_qstop;              // excessive q indicator
    const Real_t m_monoq_max_slope;
    const Real_t m_monoq_limiter_mult;
    const Real_t m_qlc_monoq;          // linear term coef for q
    const Real_t m_qqc_monoq;          // quadratic term coef for q
    const Real_t m_qqc;
    const Real_t m_eosvmax;
    const Real_t m_eosvmin;
    const Real_t m_pmin;               // pressure floor
    const Real_t m_emin;               // energy floor
    const Real_t m_dvovmax;            // maximum allowable volume change
    const Real_t m_refdens;            // reference density

    // Variables to keep track of timestep, simulation time, and cycle
    Real_t m_dtcourant;          // courant constraint
    Real_t m_dthydro;            // volume change constraint
    Int_t m_cycle;               // iteration count for simulation
    Real_t m_dtfixed;            // fixed time increment
    Real_t m_time;               // current time
    Real_t m_deltatime;          // variable time increment
    Real_t m_deltatimemultlb;
    Real_t m_deltatimemultub;
    Real_t m_dtmax;              // maximum allowable time increment
    Real_t m_stoptime;           // end time for simulation


    Int_t m_numRanks;

    Index_t m_colLoc;
    Index_t m_rowLoc;
    Index_t m_planeLoc;
    Index_t m_tp;

    Index_t m_sizeX;
    Index_t m_sizeY;
    Index_t m_sizeZ;

    NodeSet m_nodeSet;
    ElemSet m_elemSet;
    CornerSet m_cornerSet;
    ExtendedElemSet m_extendedElemSet;      // Has space for a elements as well as each face on the boundary

    Index_t m_maxPlaneSize;
    Index_t m_maxEdgeSize;

    // OMP hack
    NodeToCornerRelation m_nodeCornerRelation;

    // Used in setup
    Index_t m_rowMin, m_rowMax;
    Index_t m_colMin, m_colMax;
    Index_t m_planeMin, m_planeMax;

    RealsRegistry m_realsRegistry;
    IntsRegistry m_intsRegistry;
  };

  typedef Real_t & (Domain::* Domain_member )(Index_t);

  struct cmdLineOpts
  {
    Int_t its; // -i
    Int_t nx; // -s
    Int_t numReg; // -r
    Int_t numFiles; // -f
    Int_t showProg; // -p
    Int_t quiet; // -q
    Int_t viz; // -v
    Int_t cost; // -c
    Int_t balance; // -b
  };



// Function Prototypes

// lulesh-par
  Real_t CalcElemVolume( const Real_t x[8],
      const Real_t y[8],
      const Real_t z[8]);

// lulesh-util
  void  ParseCommandLineOptions(int argc, char *argv[],
      int myRank, struct cmdLineOpts *opts);
  void  VerifyAndWriteFinalOutput(Real_t elapsed_time,
      Domain& locDom,
      Int_t nx,
      Int_t numRanks);

// lulesh-viz
  void  DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);

// lulesh-comm
  void  CommRecv(Domain& domain, int msgType, Index_t xferFields, Index_t dx, Index_t dy, Index_t dz, bool doRecv, bool planeOnly);
  void  CommSend(Domain& domain, int msgType, Index_t xferFields, Domain_member *fieldData, Index_t dx, Index_t dy, Index_t dz, bool doSend, bool planeOnly);
  void  CommSBN(Domain& domain, int xferFields, Domain_member *fieldData);
  void  CommSyncPosVel(Domain& domain);
  void  CommMonoQ(Domain& domain);

// lulesh-init
  void  InitMeshDecomp(Int_t numRanks, Int_t myRank, Int_t *col, Int_t *row, Int_t *plane, Int_t *side);

} // end namespace slamLulesh
