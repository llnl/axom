// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


// OpenMP will be compiled in if this flag is set to 1 AND the compiler beging
// used supports it (i.e. the _OPENMP symbol is defined)
#define USE_OMP 1

#include "axom/sidre.hpp"

#ifdef AXOM_USE_MPI
#include <mpi.h>

/*
   define one of these three symbols:

   SEDOV_SYNC_POS_VEL_NONE
   SEDOV_SYNC_POS_VEL_EARLY
   SEDOV_SYNC_POS_VEL_LATE
 */

#define SEDOV_SYNC_POS_VEL_EARLY 1
#endif

#include <math.h>
#include <vector>

//**************************************************
// Allow flexibility for arithmetic representations
//**************************************************

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))


// Precision specification
typedef float real4;
typedef double real8;
typedef long double real10;    // 10 bytes on x86

typedef int Index_t;     // array subscript and loop index
typedef real8 Real_t;    // floating point representation
typedef int Int_t;       // integer representation

enum { VolumeError = -1, QStopError = -2 };

inline real4  SQRT(real4 arg) {
  return sqrtf(arg);
}
inline real8  SQRT(real8 arg) {
  return sqrt(arg);
}
inline real10 SQRT(real10 arg) {
  return sqrtl(arg);
}

inline real4  CBRT(real4 arg) {
  return cbrtf(arg);
}
inline real8  CBRT(real8 arg) {
  return cbrt(arg);
}
inline real10 CBRT(real10 arg) {
  return cbrtl(arg);
}

inline real4  FABS(real4 arg) {
  return fabsf(arg);
}
inline real8  FABS(real8 arg) {
  return fabs(arg);
}
inline real10 FABS(real10 arg) {
  return fabsl(arg);
}


// Stuff needed for boundary conditions
// 2 BCs on each of 6 hexahedral faces (12 bits)
#define XI_M        0x00007
#define XI_M_SYMM   0x00001
#define XI_M_FREE   0x00002
#define XI_M_COMM   0x00004

#define XI_P        0x00038
#define XI_P_SYMM   0x00008
#define XI_P_FREE   0x00010
#define XI_P_COMM   0x00020

#define ETA_M       0x001c0
#define ETA_M_SYMM  0x00040
#define ETA_M_FREE  0x00080
#define ETA_M_COMM  0x00100

#define ETA_P       0x00e00
#define ETA_P_SYMM  0x00200
#define ETA_P_FREE  0x00400
#define ETA_P_COMM  0x00800

#define ZETA_M      0x07000
#define ZETA_M_SYMM 0x01000
#define ZETA_M_FREE 0x02000
#define ZETA_M_COMM 0x04000

#define ZETA_P      0x38000
#define ZETA_P_SYMM 0x08000
#define ZETA_P_FREE 0x10000
#define ZETA_P_COMM 0x20000

// MPI Message Tags
#define MSG_COMM_SBN      1024
#define MSG_SYNC_POS_VEL  2048
#define MSG_MONOQ         3072

#define MAX_FIELDS_PER_MPI_COMM 6

// Assume 128 byte coherence
// Assume Real_t is an "integral power of 2" bytes wide
#define CACHE_COHERENCE_PAD_REAL (128 / sizeof(Real_t))

#define CACHE_ALIGN_REAL(n) \
  (((n) + (CACHE_COHERENCE_PAD_REAL - 1)) & ~(CACHE_COHERENCE_PAD_REAL-1))

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
 * m_x, m_y, m_z could be budled into
 *
 *    struct { Real_t x, y, z ; } *m_coord ;
 *
 * allowing accessor functions such as
 *
 *  "Real_t &x(Index_t idx) { return m_coord[idx].x ; }"
 *  "Real_t &y(Index_t idx) { return m_coord[idx].y ; }"
 *  "Real_t &z(Index_t idx) { return m_coord[idx].z ; }"
 */
#define USE_SIDRE 1

class Domain
{

public:

#if USE_SIDRE==1
  typedef Real_t * luleshRealData;
  typedef Int_t * luleshIntData;
  typedef Index_t * luleshIndexData;
#else
  typedef std::vector<Real_t> luleshRealData;
  typedef std::vector<Int_t> luleshIntData;
  typedef std::vector<Index_t> luleshIndexData;
#endif

  // Constructor
  Domain(Int_t numRanks, Index_t colLoc,
         Index_t rowLoc, Index_t planeLoc,
         Index_t nx, Int_t tp, Int_t nr, Int_t balance, Int_t cost);

  // Destructor
  ~Domain();

  //
  // ALLOCATION
  //

  void AllocateNodePersistent(Int_t numNode)  // Node-centered
  {
#if USE_SIDRE==1
    m_x = m_Group->createViewAndAllocate("m_x",axom::sidre::DataType::float64(
                                           numNode))->getData();
    m_y = m_Group->createViewAndAllocate("m_y",axom::sidre::DataType::float64(
                                           numNode))->getData();
    m_z = m_Group->createViewAndAllocate("m_z",axom::sidre::DataType::float64(
                                           numNode))->getData();

    m_xd = m_Group->createViewAndAllocate("m_xd",axom::sidre::DataType::float64(
                                            numNode))->getData();
    m_yd = m_Group->createViewAndAllocate("m_yd",axom::sidre::DataType::float64(
                                            numNode))->getData();
    m_zd = m_Group->createViewAndAllocate("m_zd",axom::sidre::DataType::float64(
                                            numNode))->getData();

    m_xdd = m_Group->createViewAndAllocate("m_xdd",axom::sidre::DataType::float64(
                                             numNode))->getData();
    m_ydd = m_Group->createViewAndAllocate("m_ydd",axom::sidre::DataType::float64(
                                             numNode))->getData();
    m_zdd = m_Group->createViewAndAllocate("m_zdd",axom::sidre::DataType::float64(
                                             numNode))->getData();

    m_fx = m_Group->createViewAndAllocate("m_fx",axom::sidre::DataType::float64(
                                            numNode))->getData();
    m_fy = m_Group->createViewAndAllocate("m_fy",axom::sidre::DataType::float64(
                                            numNode))->getData();
    m_fz = m_Group->createViewAndAllocate("m_fz",axom::sidre::DataType::float64(
                                            numNode))->getData();

    m_nodalMass = m_Group->createViewAndAllocate("m_nodalMass",axom::sidre::DataType::float64(
                                                   numNode))->getData();

#else
    m_x.resize(numNode);    // coordinates
    m_y.resize(numNode);
    m_z.resize(numNode);

    m_xd.resize(numNode);   // velocities
    m_yd.resize(numNode);
    m_zd.resize(numNode);

    m_xdd.resize(numNode);   // accelerations
    m_ydd.resize(numNode);
    m_zdd.resize(numNode);

    m_fx.resize(numNode);    // forces
    m_fy.resize(numNode);
    m_fz.resize(numNode);

    m_nodalMass.resize(numNode);    // mass
#endif
  }

  void AllocateElemPersistent(Int_t numElem)  // Elem-centered
  {
#if USE_SIDRE==1
    m_nodelist = m_Group->createViewAndAllocate("m_nodelist",axom::sidre::DataType::int32(
                                                  8*numElem))->getData();

    // elem connectivities through face
    m_lxim = m_Group->createViewAndAllocate("m_lxim",axom::sidre::DataType::int32(
                                              numElem))->getData();
    m_lxip = m_Group->createViewAndAllocate("m_lxip",axom::sidre::DataType::int32(
                                              numElem))->getData();
    m_letam = m_Group->createViewAndAllocate("m_letam",axom::sidre::DataType::int32(
                                               numElem))->getData();
    m_letap = m_Group->createViewAndAllocate("m_letap",axom::sidre::DataType::int32(
                                               numElem))->getData();
    m_lzetam = m_Group->createViewAndAllocate("m_lzetam",axom::sidre::DataType::int32(
                                                numElem))->getData();
    m_lzetap = m_Group->createViewAndAllocate("m_lzetap",axom::sidre::DataType::int32(
                                                numElem))->getData();



    m_elemBC = m_Group->createViewAndAllocate("m_elemBC",axom::sidre::DataType::int32(
                                                numElem))->getData();

    m_e = m_Group->createViewAndAllocate("m_e",axom::sidre::DataType::float64(
                                           numElem))->getData();
    m_p = m_Group->createViewAndAllocate("m_p",axom::sidre::DataType::float64(
                                           numElem))->getData();

    m_q = m_Group->createViewAndAllocate("m_q",axom::sidre::DataType::float64(
                                           numElem))->getData();
    m_ql = m_Group->createViewAndAllocate("m_ql",axom::sidre::DataType::float64(
                                            numElem))->getData();
    m_qq = m_Group->createViewAndAllocate("m_qq",axom::sidre::DataType::float64(
                                            numElem))->getData();

    m_v = m_Group->createViewAndAllocate("m_v",axom::sidre::DataType::float64(
                                           numElem))->getData();

    m_volo = m_Group->createViewAndAllocate("m_volo",axom::sidre::DataType::float64(
                                              numElem))->getData();
    m_delv = m_Group->createViewAndAllocate("m_delv",axom::sidre::DataType::float64(
                                              numElem))->getData();
    m_vdov = m_Group->createViewAndAllocate("m_vdov",axom::sidre::DataType::float64(
                                              numElem))->getData();





    m_arealg = m_Group->createViewAndAllocate("m_arealg",axom::sidre::DataType::float64(
                                                numElem))->getData();

    m_ss = m_Group->createViewAndAllocate("m_ss",axom::sidre::DataType::float64(
                                            numElem))->getData();

    m_elemMass = m_Group->createViewAndAllocate("m_elemMass",axom::sidre::DataType::float64(
                                                  numElem))->getData();

#else
    m_nodelist.resize(8*numElem);

    // elem connectivities through face
    m_lxim.resize(numElem);
    m_lxip.resize(numElem);
    m_letam.resize(numElem);
    m_letap.resize(numElem);
    m_lzetam.resize(numElem);
    m_lzetap.resize(numElem);

    m_elemBC.resize(numElem);

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
#endif
  }

  void AllocateGradients(Int_t numElem, Int_t allElem)
  {
#if USE_SIDRE==1
    // Position gradients
    m_delx_xi = m_Group->createViewAndAllocate("m_delx_xi",axom::sidre::DataType::float64(
                                                 numElem))->getData();
    m_delx_eta = m_Group->createViewAndAllocate("m_delx_eta",axom::sidre::DataType::float64(
                                                  numElem))->getData();
    m_delx_zeta = m_Group->createViewAndAllocate("m_delx_zeta",axom::sidre::DataType::float64(
                                                   numElem))->getData();

    // Velocity gradients
    m_delv_xi = m_Group->createViewAndAllocate("m_delv_xi",axom::sidre::DataType::float64(
                                                 allElem))->getData();
    m_delv_eta = m_Group->createViewAndAllocate("m_delv_eta",axom::sidre::DataType::float64(
                                                  allElem))->getData();
    m_delv_zeta = m_Group->createViewAndAllocate("m_delv_zeta",axom::sidre::DataType::float64(
                                                   allElem))->getData();
#else
    // Position gradients
    m_delx_xi.resize(numElem);
    m_delx_eta.resize(numElem);
    m_delx_zeta.resize(numElem);

    // Velocity gradients
    m_delv_xi.resize(allElem);
    m_delv_eta.resize(allElem);
    m_delv_zeta.resize(allElem);
#endif
  }

  void DeallocateGradients()
  {
#if USE_SIDRE==1
    m_Group->destroyView("m_delx_zeta");
    m_Group->destroyView("m_delx_eta");
    m_Group->destroyView("m_delx_xi");

    m_Group->destroyView("m_delv_zeta");
    m_Group->destroyView("m_delv_eta");
    m_Group->destroyView("m_delv_xi");
#else
    m_delx_zeta.clear();
    m_delx_eta.clear();
    m_delx_xi.clear();

    m_delv_zeta.clear();
    m_delv_eta.clear();
    m_delv_xi.clear();
#endif
  }

  void AllocateStrains(Int_t numElem)
  {
#if USE_SIDRE==1
    m_dxx = m_Group->createViewAndAllocate("m_dxx",axom::sidre::DataType::float64(
                                             numElem))->getData();
    m_dyy = m_Group->createViewAndAllocate("m_dyy",axom::sidre::DataType::float64(
                                             numElem))->getData();
    m_dzz = m_Group->createViewAndAllocate("m_dzz",axom::sidre::DataType::float64(
                                             numElem))->getData();
#else
    m_dxx.resize(numElem);
    m_dyy.resize(numElem);
    m_dzz.resize(numElem);
#endif
  }

  void DeallocateStrains()
  {
#if USE_SIDRE==1
    m_Group->destroyView("m_dxx");
    m_Group->destroyView("m_dyy");
    m_Group->destroyView("m_dzz");
#else
    m_dzz.clear();
    m_dyy.clear();
    m_dxx.clear();
#endif
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

  // Nodes on symmertry planes
  Index_t symmX(Index_t idx) { return m_symmX[idx]; }
  Index_t symmY(Index_t idx) { return m_symmY[idx]; }
  Index_t symmZ(Index_t idx) { return m_symmZ[idx]; }
  bool symmXempty()          {
    return !(m_Group->hasView("m_symmX")
             && (m_Group->getView("m_symmX")->getNumElements() > 0) );
  }
  bool symmYempty()          {
    return !(m_Group->hasView("m_symmY")
             && (m_Group->getView("m_symmY")->getNumElements() > 0) );
  }
  bool symmZempty()          {
    return !(m_Group->hasView("m_symmZ")
             && (m_Group->getView("m_symmZ")->getNumElements() > 0) );
  }

  //
  // Element-centered
  //
  Index_t&  regElemSize(Index_t idx) { return m_regElemSize[idx]; }
  Index_t&  regNumList(Index_t idx) { return m_regNumList[idx]; }
  Index_t * regNumList()            { return &m_regNumList[0]; }
  Index_t * regElemlist(Int_t r)    { return m_regElemlist[r]; }
  Index_t&  regElemlist(Int_t r, Index_t idx) { return m_regElemlist[r][idx]; }

  Index_t * nodelist(Index_t idx)    { return &m_nodelist[Index_t(8)*idx]; }

  // elem connectivities through face
  Index_t&  lxim(Index_t idx) { return m_lxim[idx]; }
  Index_t&  lxip(Index_t idx) { return m_lxip[idx]; }
  Index_t&  letam(Index_t idx) { return m_letam[idx]; }
  Index_t&  letap(Index_t idx) { return m_letap[idx]; }
  Index_t&  lzetam(Index_t idx) { return m_lzetam[idx]; }
  Index_t&  lzetap(Index_t idx) { return m_lzetap[idx]; }

  // elem face symm/free-surface flag
  Int_t&  elemBC(Index_t idx) { return m_elemBC[idx]; }

  // Principal strains - temporary
  Real_t& dxx(Index_t idx)  { return m_dxx[idx]; }
  Real_t& dyy(Index_t idx)  { return m_dyy[idx]; }
  Real_t& dzz(Index_t idx)  { return m_dzz[idx]; }

  // Velocity gradient - temporary
  Real_t& delv_xi(Index_t idx)    { return m_delv_xi[idx]; }
  Real_t& delv_eta(Index_t idx)   { return m_delv_eta[idx]; }
  Real_t& delv_zeta(Index_t idx)  { return m_delv_zeta[idx]; }

  // Position gradient - temporary
  Real_t& delx_xi(Index_t idx)    { return m_delx_xi[idx]; }
  Real_t& delx_eta(Index_t idx)   { return m_delx_eta[idx]; }
  Real_t& delx_zeta(Index_t idx)  { return m_delx_zeta[idx]; }

  // Energy
  Real_t& e(Index_t idx)          { return m_e[idx]; }

  // Pressure
  Real_t& p(Index_t idx)          { return m_p[idx]; }

  // Artificial viscosity
  Real_t& q(Index_t idx)          { return m_q[idx]; }

  // Linear term for q
  Real_t& ql(Index_t idx)         { return m_ql[idx]; }
  // Quadratic term for q
  Real_t& qq(Index_t idx)         { return m_qq[idx]; }

  // Relative volume
  Real_t& v(Index_t idx)          { return m_v[idx]; }
  Real_t& delv(Index_t idx)       { return m_delv[idx]; }

  // Reference volume
  Real_t& volo(Index_t idx)       { return m_volo[idx]; }

  // volume derivative over volume
  Real_t& vdov(Index_t idx)       { return m_vdov[idx]; }

  // Element characteristic length
  Real_t& arealg(Index_t idx)     { return m_arealg[idx]; }

  // Sound speed
  Real_t& ss(Index_t idx)         { return m_ss[idx]; }

  // Element mass
  Real_t& elemMass(Index_t idx)  { return m_elemMass[idx]; }

  Index_t nodeElemCount(Index_t idx)
  { return m_nodeElemStart[idx+1] - m_nodeElemStart[idx]; }

  Index_t * nodeElemCornerList(Index_t idx)
  { return &m_nodeElemCornerList[m_nodeElemStart[idx]]; }

  // Parameters

  // Cutoffs
  Real_t u_cut() const { return m_u_cut; }
  Real_t e_cut() const { return m_e_cut; }
  Real_t p_cut() const { return m_p_cut; }
  Real_t q_cut() const { return m_q_cut; }
  Real_t v_cut() const { return m_v_cut; }

  // Other constants (usually are settable via input file in real codes)
  Real_t hgcoef() const { return m_hgcoef; }
  Real_t qstop() const { return m_qstop; }
  Real_t monoq_max_slope() const { return m_monoq_max_slope; }
  Real_t monoq_limiter_mult() const { return m_monoq_limiter_mult; }
  Real_t ss4o3() const { return m_ss4o3; }
  Real_t qlc_monoq() const { return m_qlc_monoq; }
  Real_t qqc_monoq() const { return m_qqc_monoq; }
  Real_t qqc() const { return m_qqc; }

  Real_t eosvmax() const { return m_eosvmax; }
  Real_t eosvmin() const { return m_eosvmin; }
  Real_t pmin() const { return m_pmin; }
  Real_t emin() const { return m_emin; }
  Real_t dvovmax() const { return m_dvovmax; }
  Real_t refdens() const { return m_refdens; }

  // Timestep controls, etc...
  Real_t& time()                 { return m_time; }
  Real_t& deltatime()            { return m_deltatime; }
  Real_t& deltatimemultlb()      { return m_deltatimemultlb; }
  Real_t& deltatimemultub()      { return m_deltatimemultub; }
  Real_t& stoptime()             { return m_stoptime; }
  Real_t& dtcourant()            { return m_dtcourant; }
  Real_t& dthydro()              { return m_dthydro; }
  Real_t& dtmax()                { return m_dtmax; }
  Real_t& dtfixed()              { return m_dtfixed; }

  Int_t&  cycle()                { return m_cycle; }
  Index_t&  numRanks()           { return m_numRanks; }

  Index_t&  colLoc()             { return m_colLoc; }
  Index_t&  rowLoc()             { return m_rowLoc; }
  Index_t&  planeLoc()           { return m_planeLoc; }
  Index_t&  tp()                 { return m_tp; }

  Index_t&  sizeX()              { return m_sizeX; }
  Index_t&  sizeY()              { return m_sizeY; }
  Index_t&  sizeZ()              { return m_sizeZ; }
  Index_t&  numReg()             { return m_numReg; }
  Int_t&  cost()                 { return m_cost; }
  Index_t&  numElem()            { return m_numElem; }
  Index_t&  numNode()            { return m_numNode; }

  Index_t&  maxPlaneSize()       { return m_maxPlaneSize; }
  Index_t&  maxEdgeSize()        { return m_maxEdgeSize; }

  //
  // MPI-Related additional data
  //

#ifdef AXOM_USE_MPI
  // Communication Work space
  Real_t * commDataSend {nullptr};
  Real_t * commDataRecv {nullptr};

  // Maximum number of block neighbors
  MPI_Request recvRequest[26];   // 6 faces + 12 edges + 8 corners
  MPI_Request sendRequest[26];   // 6 faces + 12 edges + 8 corners
#endif

private:

  void BuildMesh(Int_t nx, Int_t edgeNodes, Int_t edgeElems);
  void SetupThreadSupportStructures();
  void CreateRegionIndexSets(Int_t nreg, Int_t balance);
  void SetupCommBuffers(Int_t edgeNodes);
  void SetupSymmetryPlanes(Int_t edgeNodes);
  void SetupElementConnectivities(Int_t edgeElems);
  void SetupBoundaryConditions(Int_t edgeElems);

  //
  // IMPLEMENTATION
  //
#if USE_SIDRE==1
  axom::sidre::DataStore m_DataStore;
  axom::sidre::Group * m_Group;

#endif

  /* Node-centered */
  luleshRealData m_x;    /* coordinates */
  luleshRealData m_y;
  luleshRealData m_z;

  luleshRealData m_xd;   /* velocities */
  luleshRealData m_yd;
  luleshRealData m_zd;

  luleshRealData m_xdd;   /* accelerations */
  luleshRealData m_ydd;
  luleshRealData m_zdd;

  luleshRealData m_fx;    /* forces */
  luleshRealData m_fy;
  luleshRealData m_fz;

  luleshRealData m_nodalMass;    /* mass */

  luleshIndexData m_symmX;    /* symmetry plane nodesets */
  luleshIndexData m_symmY;
  luleshIndexData m_symmZ;

  // Element-centered

  // Region information
  Int_t m_numReg;
  Int_t m_cost;     //imbalance cost
  Index_t * m_regElemSize {nullptr};    // Size of region sets
  Index_t * m_regNumList {nullptr};     // Region number per domain element
  Index_t * * m_regElemlist {nullptr};  // region indexset

  luleshIndexData m_nodelist;        /* elemToNode connectivity */

  luleshIndexData m_lxim;     /* element connectivity across each face */
  luleshIndexData m_lxip;
  luleshIndexData m_letam;
  luleshIndexData m_letap;
  luleshIndexData m_lzetam;
  luleshIndexData m_lzetap;

  luleshIntData m_elemBC;       /* symmetry/free-surface flags for each elem face */

  luleshRealData m_dxx;    /* principal strains -- temporary */
  luleshRealData m_dyy;
  luleshRealData m_dzz;

  luleshRealData m_delv_xi;      /* velocity gradient -- temporary */
  luleshRealData m_delv_eta;
  luleshRealData m_delv_zeta;

  luleshRealData m_delx_xi;      /* coordinate gradient -- temporary */
  luleshRealData m_delx_eta;
  luleshRealData m_delx_zeta;

  luleshRealData m_e;     /* energy */

  luleshRealData m_p;     /* pressure */
  luleshRealData m_q;     /* q */
  luleshRealData m_ql;    /* linear term for q */
  luleshRealData m_qq;    /* quadratic term for q */

  luleshRealData m_v;       /* relative volume */
  luleshRealData m_volo;    /* reference volume */
  //luleshRealData m_vnew ;  /* new relative volume -- temporary */
  luleshRealData m_delv;    /* m_vnew - m_v */
  luleshRealData m_vdov;    /* volume derivative over volume */

  luleshRealData m_arealg;    /* characteristic length of an element */

  luleshRealData m_ss;        /* "sound speed" */

  luleshRealData m_elemMass;    /* mass */

  // Cutoffs (treat as constants)
  const Real_t m_e_cut;                // energy tolerance
  const Real_t m_p_cut;                // pressure tolerance
  const Real_t m_q_cut;                // q tolerance
  const Real_t m_v_cut;                // relative volume tolerance
  const Real_t m_u_cut;                // velocity tolerance

  // Other constants (usually setable, but hardcoded in this proxy app)

  const Real_t m_hgcoef;               // hourglass control
  const Real_t m_ss4o3;
  const Real_t m_qstop;                // excessive q indicator
  const Real_t m_monoq_max_slope;
  const Real_t m_monoq_limiter_mult;
  const Real_t m_qlc_monoq;            // linear term coef for q
  const Real_t m_qqc_monoq;            // quadratic term coef for q
  const Real_t m_qqc;
  const Real_t m_eosvmax;
  const Real_t m_eosvmin;
  const Real_t m_pmin;                 // pressure floor
  const Real_t m_emin;                 // energy floor
  const Real_t m_dvovmax;              // maximum allowable volume change
  const Real_t m_refdens;              // reference density

  // Variables to keep track of timestep, simulation time, and cycle
  Real_t m_dtcourant;            // courant constraint
  Real_t m_dthydro;              // volume change constraint
  Int_t m_cycle;                 // iteration count for simulation
  Real_t m_dtfixed;              // fixed time increment
  Real_t m_time;                 // current time
  Real_t m_deltatime;            // variable time increment
  Real_t m_deltatimemultlb;
  Real_t m_deltatimemultub;
  Real_t m_dtmax;                // maximum allowable time increment
  Real_t m_stoptime;             // end time for simulation


  Int_t m_numRanks;

  Index_t m_colLoc;
  Index_t m_rowLoc;
  Index_t m_planeLoc;
  Index_t m_tp;

  Index_t m_sizeX;
  Index_t m_sizeY;
  Index_t m_sizeZ;
  Index_t m_numElem;
  Index_t m_numNode;

  Index_t m_maxPlaneSize;
  Index_t m_maxEdgeSize;

  // OMP hack
  Index_t * m_nodeElemStart;
  Index_t * m_nodeElemCornerList;

  // Used in setup
  Index_t m_rowMin, m_rowMax;
  Index_t m_colMin, m_colMax;
  Index_t m_planeMin, m_planeMax;

};

typedef Real_t & (Domain::* Domain_member )(Index_t);

struct cmdLineOpts
{
  Int_t its;  // -i
  Int_t nx;   // -s
  Int_t numReg;  // -r
  Int_t numFiles;  // -f
  Int_t showProg;  // -p
  Int_t quiet;  // -q
  Int_t viz;  // -v
  Int_t cost;  // -c
  Int_t balance;  // -b
};



// Function Prototypes

// lulesh-par
Real_t CalcElemVolume( const Real_t x[8],
                       const Real_t y[8],
                       const Real_t z[8]);

// lulesh-util
void ParseCommandLineOptions(int argc, char * argv[],
                             Int_t myRank, struct cmdLineOpts * opts);
void VerifyAndWriteFinalOutput(Real_t elapsed_time,
                               Domain& locDom,
                               Int_t nx,
                               Int_t numRanks);

// lulesh-viz
void DumpToVisit(Domain& domain, int numFiles, int myRank, int numRanks);

// lulesh-comm
void CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
              Index_t dx, Index_t dy, Index_t dz,
              bool doRecv, bool planeOnly);
void CommSend(Domain& domain, Int_t msgType,
              Index_t xferFields, Domain_member * fieldData,
              Index_t dx, Index_t dy, Index_t dz,
              bool doSend, bool planeOnly);
void CommSBN(Domain& domain, Int_t xferFields, Domain_member * fieldData);
void CommSyncPosVel(Domain& domain);
void CommMonoQ(Domain& domain);

// lulesh-init
void InitMeshDecomp(Int_t numRanks, Int_t myRank,
                    Int_t * col, Int_t * row, Int_t * plane, Int_t * side);
