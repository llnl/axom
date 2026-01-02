// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"

#ifdef AXOM_USE_MFEM

  #include <string>
  #include <iomanip>      // for setw, setfill
  #include <cstdio>       // for snprintf()
  #include <type_traits>  // for checking layout

  #include "conduit_blueprint.hpp"
  #include "axom/fmt.hpp"

  #include "MFEMSidreDataCollection.hpp"
  #ifdef AXOM_USE_MPI
    #include "axom/sidre/spio/IOManager.hpp"
  #endif

  #include "axom/core/utilities/StringUtilities.hpp"
  #include "axom/core/utilities/Utilities.hpp"

using mfem::Array;
using mfem::GridFunction;
using mfem::Mesh;
using mfem::Ordering;

namespace axom
{
namespace sidre
{
// Initializations for the coordset/topology/field names
// used in the blueprint group
// FIXME: These are still hardcoded in some places
const std::string MFEMSidreDataCollection::s_mesh_topology_name = "mesh";
const std::string MFEMSidreDataCollection::s_boundary_topology_name = "boundary";
const std::string MFEMSidreDataCollection::s_attribute_suffix = "_material_attribute";
const std::string MFEMSidreDataCollection::s_coordset_name = "coords";

namespace detail
{
/**
 * @brief Implements an analogue to mfem::FiniteElementCollection::New
 * for mfem::QuadratureSpaces - basis currently must be of form QF_Default_[ORDER]_[VDIM]
 * @param[in] name The name that encodes the QuadratureSpace data
 * @param[in] mesh The mesh to construct the QuadratureSpace with
 * @param[out] vdim The vector dimension parsed from the basis name - unlike FEColl/FESpace,
 * QSpaces do not contain vdim information, but it may be useful for the caller when constructing
 * the QFunc
 */
mfem::QuadratureSpace* NewQuadratureSpace(const std::string& name, Mesh* mesh, int& vdim)
{
  const auto tokens = utilities::string::rsplitN(name, 4, '_');
  // Uses raw pointers for consistency with MFEM
  mfem::QuadratureSpace* qspace = nullptr;
  if((tokens.size() == 4) && (tokens[0] == "QF"))
  {
    // Right now only the Default type is supported
    if(tokens[1] == "Default")
    {
      if(conduit::utils::string_is_integer(tokens[2]) && conduit::utils::string_is_integer(tokens[3]))
      {
        int order = conduit::utils::string_to_value<int>(tokens[2]);
        vdim = conduit::utils::string_to_value<int>(tokens[3]);
        qspace = new mfem::QuadratureSpace(mesh, order);
      }
    }
  }
  SLIC_ERROR_IF(!qspace, "Unrecognized QuadratureSpace name: " << name);
  return qspace;
}

}  // namespace detail

// Constructor that will automatically create the sidre data store and necessary
// data groups for domain and global data.
MFEMSidreDataCollection::MFEMSidreDataCollection(const std::string& collection_name,
                                                 Mesh* the_mesh,
                                                 bool own_mesh_data)
  : mfem::DataCollection(collection_name, the_mesh)
  , m_owns_datastore(true)
  , m_owns_mesh_data(own_mesh_data)
{
  m_datastore_ptr = new sidre::DataStore();

  sidre::Group* global_grp = m_datastore_ptr->getRoot()->createGroup(collection_name + "_global");
  sidre::Group* domain_grp = m_datastore_ptr->getRoot()->createGroup(collection_name);

  m_bp_grp = domain_grp->createGroup("blueprint");
  // Currently only rank 0 adds anything to bp_index.
  m_bp_index_grp = global_grp->createGroup("blueprint_index/" + name);

  m_named_bufs_grp = domain_grp->createGroup("named_buffers");

  this->own_data = own_mesh_data;

  if(the_mesh)
  {
    MFEMSidreDataCollection::SetMesh(the_mesh);
  }
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  else
  {
    m_comm = MPI_COMM_NULL;
  }
  #endif
}

// Second constructor that allows external code to specify data groups to place
// domain and global data in.

// TODO - Conduit will have the capability to generate a blueprint index group
// in the future.  When this is available, all the blueprint index code can be
// removed from the data collection class.
MFEMSidreDataCollection::MFEMSidreDataCollection(const std::string& collection_name,
                                                 Group* bp_index_grp,
                                                 Group* domain_grp,
                                                 bool own_mesh_data)
  : mfem::DataCollection(collection_name)
  , m_owns_datastore(false)
  , m_owns_mesh_data(own_mesh_data)
  , m_bp_index_grp(bp_index_grp)
{
  m_bp_grp = domain_grp->createGroup("blueprint");

  m_named_bufs_grp = domain_grp->createGroup("named_buffers");

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  m_comm = MPI_COMM_NULL;
  #endif

  this->own_data = false;
}

MFEMSidreDataCollection::~MFEMSidreDataCollection()
{
  attr_map.DeleteData(true);
  if(m_owns_datastore)
  {
    delete m_datastore_ptr;
  }
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
void MFEMSidreDataCollection::SetComm(MPI_Comm comm)
{
  m_comm = comm;
  serial = false;
  appendRankToFileName = true;
  MPI_Comm_rank(m_comm, &myid);
  MPI_Comm_size(m_comm, &num_procs);
}

void MFEMSidreDataCollection::SetNumFiles(int num_files)
{
  SLIC_ASSERT_MSG(num_files > 0, "Output must be to at least 1 file");
  m_num_files = num_files;
}
  #endif

// protected method
sidre::Group* MFEMSidreDataCollection::named_buffers_grp() const
{
  SLIC_ASSERT_MSG(m_named_bufs_grp != nullptr,
                  "No group 'named_buffers' in data collection.  Verify that"
                  " SetMesh was called to set the mesh in the data collection.");
  return m_named_bufs_grp;
}

// protected method
View* MFEMSidreDataCollection::alloc_view(Group* grp, const std::string& view_name)
{
  SLIC_ASSERT_MSG(grp, "Group pointer is nullptr");
  sidre::View* v = nullptr;

  if(!grp->hasView(view_name))
  {
    v = grp->createView(view_name);
    SLIC_ASSERT_MSG(v, "error allocating View " << view_name << " in group " << grp->getPathName());
  }
  else
  {
    v = grp->getView(view_name);
  }

  return v;
}

// protected method
View* MFEMSidreDataCollection::alloc_view(Group* grp,
                                          const std::string& view_name,
                                          const DataType& dtype)
{
  SLIC_ASSERT_MSG(grp, "Group pointer is nullptr");
  sidre::View* v = nullptr;

  if(!grp->hasView(view_name))
  {
    v = grp->createView(view_name, dtype);
    SLIC_ASSERT_MSG(v, "error allocating View " << view_name << " in group " << grp->getPathName());
  }
  else
  {
    v = grp->getView(view_name);
    SLIC_ASSERT_MSG(v->getSchema().dtype().equals(dtype), "");
  }
  return v;
}

// protected method
Group* MFEMSidreDataCollection::alloc_group(Group* grp, const std::string& group_name)
{
  SLIC_ASSERT_MSG(grp, "Group pointer is nullptr");
  sidre::Group* g = nullptr;

  if(!grp->hasGroup(group_name))
  {
    g = grp->createGroup(group_name);
    SLIC_ASSERT_MSG(g, "error allocating Group " << group_name << " in group " << grp->getPathName());
  }
  else
  {
    g = grp->getGroup(group_name);
  }
  return g;
}

// protected method
std::string MFEMSidreDataCollection::get_file_path(const std::string& filename) const
{
  std::stringstream fNameSstr;

  namespace fs = axom::utilities::filesystem;
  fNameSstr << fs::joinPath(prefix_path, filename);

  if(GetCycle() >= 0)
  {
    fNameSstr << "_" << std::setfill('0') << std::setw(pad_digits_cycle) << GetCycle();
  }

  return fNameSstr.str();
}

View* MFEMSidreDataCollection::AllocNamedBuffer(const std::string& buffer_name,
                                                IndexType sz,
                                                TypeID type)
{
  sz = std::max(sz, sidre::IndexType(0));
  sidre::Group* f = named_buffers_grp();
  sidre::View* v = nullptr;

  if(!f->hasView(buffer_name))
  {
    // create a buffer view
    v = f->createViewAndAllocate(buffer_name, type, sz);
  }
  else
  {
    v = f->getView(buffer_name);
    SLIC_ASSERT_MSG(v->getTypeID() == type, "type does not match existing type");

    // Here v is the view holding the buffer in the named_buffers group, so
    // its size is the full size of the buffer.

    // check if we need to resize.
    if(!v->isApplied() || v->getNumElements() < sz)
    {
      // resize, even if the buffer has more than 1 View.
      // v->reallocate(sz); // this will not work for more than 1 view.
      sidre::DataType dtype(v->getSchema().dtype());
      dtype.set_number_of_elements(sz);
      f->destroyViewAndData(buffer_name);
      v = f->createViewAndAllocate(buffer_name, dtype);
    }
  }
  SLIC_ASSERT_MSG(v && v->isApplied(), "allocation failed");
  return v;
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintStubs(bool hasBP)
{
  if(!hasBP)
  {
    m_bp_grp->createGroup("state");
    m_bp_grp->createGroup("coordsets");
    m_bp_grp->createGroup("topologies");
    m_bp_grp->createGroup("fields");
  }

  // If rank is 0, set up blueprint index state group.
  if(myid == 0)
  {
    m_bp_index_grp->createGroup("state");
    m_bp_index_grp->createGroup("coordsets");
    m_bp_index_grp->createGroup("topologies");
    m_bp_index_grp->createGroup("fields");
  }
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintState(bool hasBP)
{
  if(!hasBP)
  {
    // Set up blueprint state group.
    m_bp_grp->createViewScalar("state/cycle", 0);
    m_bp_grp->createViewScalar("state/time", 0.);
    m_bp_grp->createViewScalar("state/domain", myid);
    m_bp_grp->createViewScalar("state/time_step", 0.);
  }

  // If rank is 0, set up blueprint index state group.
  if(myid == 0)
  {
    m_bp_index_grp->createViewScalar("state/cycle", 0);
    m_bp_index_grp->createViewScalar("state/time", 0.);
    m_bp_index_grp->createViewScalar("state/number_of_domains", num_procs);
  }
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintCoordset(bool hasBP)
{
  int dim = mesh->SpaceDimension();
  SLIC_ASSERT_MSG(dim >= 1 && dim <= 3, "invalid mesh dimension");

  // Assuming mfem::Vertex has the layout of a double array.
  const int NUM_COORDS = sizeof(mfem::Vertex) / sizeof(double);

  const int num_vertices = mesh->GetNV();
  const int coordset_len = NUM_COORDS * num_vertices;

  // Add blueprint if not present
  if(!hasBP)
  {
    m_bp_grp->createViewString("coordsets/coords/type", "explicit");

    sidre::DataType dtype = sidre::DataType::c_double(num_vertices);
    const size_t stride = dtype.stride();
    dtype.set_stride(stride * NUM_COORDS);

    // Set up views for x, y, z values
    sidre::View *vx, *vy = nullptr, *vz = nullptr;
    vx = m_bp_grp->createView("coordsets/coords/values/x", dtype);

    if(dim >= 2)
    {
      dtype.set_offset(dtype.offset() + stride);
      vy = m_bp_grp->createView("coordsets/coords/values/y", dtype);
    }
    if(dim >= 3)
    {
      dtype.set_offset(dtype.offset() + stride);
      vz = m_bp_grp->createView("coordsets/coords/values/z", dtype);
    }

    if(m_owns_mesh_data)
    {
      // Allocate buffer for coord values.
      sidre::Buffer* coordbuf = AllocNamedBuffer("vertex_coords", coordset_len)->getBuffer();

      vx->attachBuffer(coordbuf);
      if(dim >= 2)
      {
        vy->attachBuffer(coordbuf);
      }
      if(dim >= 3)
      {
        vz->attachBuffer(coordbuf);
      }
    }
    else
    {
      double* coordbuf = mesh->GetVertex(0);

      vx->setExternalDataPtr(coordbuf);
      if(dim >= 2)
      {
        vy->setExternalDataPtr(coordbuf);
      }
      if(dim >= 3)
      {
        vz->setExternalDataPtr(coordbuf);
      }
    }
  }

  // If rank 0, set up blueprint index for coordinate set.
  if(myid == 0)
  {
    m_bp_index_grp->createViewString("coordsets/coords/path",
                                     m_bp_grp->getPathName() + "/coordsets/coords");

    m_bp_index_grp->getGroup("coordsets/coords")
      ->copyView(m_bp_grp->getView("coordsets/coords/type"));

    m_bp_index_grp->createViewString("coordsets/coords/coord_system/type", "cartesian");

    // These are empty views, their existence in the group tree is used to
    // define the number of dims
    m_bp_index_grp->createView("coordsets/coords/coord_system/axes/x");

    if(dim >= 2)
    {
      m_bp_index_grp->createView("coordsets/coords/coord_system/axes/y");
    }

    if(dim == 3)
    {
      m_bp_index_grp->createView("coordsets/coords/coord_system/axes/z");
    }
  }

  if(m_owns_mesh_data)
  {
    double* coord_values = GetNamedBuffer("vertex_coords")->getData();
    // Change ownership of the mesh vertex data to sidre
    mesh->ChangeVertexDataOwnership(coord_values, coordset_len, hasBP);
  }
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintTopologies(bool hasBP, const std::string& mesh_name)
{
  const bool isBdry = (mesh_name == s_boundary_topology_name);

  const int num_elements = !isBdry ? mesh->GetNE() : mesh->GetNBE();

  const std::string mesh_topo_str = "topologies/" + mesh_name;
  const std::string mesh_attr_str = mesh_name + s_attribute_suffix;

  int num_indices = 0;
  int geom = 0;
  std::string eltTypeStr = "point";

  if(num_elements > 0)
  {
    mfem::Element* elem = !isBdry ? mesh->GetElement(0) : mesh->GetBdrElement(0);
    int element_size = elem->GetNVertices();

    num_indices = num_elements * element_size;

    // Find the element shape
    // Note: Assumes homogeneous elements, so only check the first element
    geom = elem->GetGeometryType();
    eltTypeStr = getElementName(elem->GetType());
  }

  // Create the blueprint "topology" group, if not present
  if(!hasBP)
  {
    sidre::Group* topology_grp = m_bp_grp->createGroup(mesh_topo_str);

    topology_grp->createViewString("type", "unstructured");
    topology_grp->createViewString("elements/shape", eltTypeStr);
    topology_grp->createViewAndAllocate("elements/connectivity", sidre::INT_ID, num_indices);
    topology_grp->createViewString("coordset", "coords");

    // If the mesh has nodes, set the name of the GridFunction holding the
    // mesh nodes in the blueprint group.
    if(!isBdry && mesh->GetNodes() != nullptr)
    {
      topology_grp->createViewString("grid_function", m_meshNodesGFName);
    }
  }

  // Add the mesh's attributes as an attribute field
  RegisterAttributeField(mesh_attr_str, isBdry);

  // Change ownership or copy the element arrays into Sidre
  if(num_elements > 0)
  {
    sidre::View* conn_view = m_bp_grp->getGroup(mesh_topo_str)->getView("elements/connectivity");

    // The MFEMSidreDataCollection always owns these arrays:
    mfem::Array<int> conn_array(conn_view->getData<int*>(), num_indices);
    mfem::Array<int>* attr_array = attr_map.Get(mesh_attr_str);
    if(!isBdry)
    {
      mesh->GetElementData(geom, conn_array, *attr_array);
    }
    else
    {
      mesh->GetBdrElementData(geom, conn_array, *attr_array);
    }
    SLIC_ASSERT_MSG(!conn_array.OwnsData(), "");
    SLIC_ASSERT_MSG(!attr_array->OwnsData(), "");
  }

  // If rank 0, set up blueprint index for topologies group
  if(myid == 0)
  {
    const std::string m_bp_grp_path = m_bp_grp->getPathName();

    if(isBdry)
    {
      // "Shallow" copy the m_bp_grp view into the m_bp_index_grp sub-group.
      // Note that the "topologies/mesh" sub-group has to exist, i.e. this
      // method should be called first with mesh_name = "mesh".
      m_bp_index_grp->getGroup("topologies/mesh")
        ->copyView(m_bp_grp->getView("topologies/mesh/boundary_topology"));
    }

    sidre::Group* bp_index_topo_grp = m_bp_index_grp->createGroup(mesh_topo_str);
    sidre::Group* topology_grp = m_bp_grp->getGroup(mesh_topo_str);

    bp_index_topo_grp->createViewString("path", m_bp_grp_path + "/" + mesh_topo_str);
    bp_index_topo_grp->copyView(topology_grp->getView("type"));
    bp_index_topo_grp->copyView(topology_grp->getView("coordset"));

    // If the mesh has nodes, set the name of the GridFunction holding the
    // mesh nodes in the blueprint_index group.
    if(!isBdry && mesh->GetNodes() != nullptr)
    {
      bp_index_topo_grp->copyView(topology_grp->getView("grid_function"));
    }
  }
}

  // private method
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
void MFEMSidreDataCollection::createMeshBlueprintAdjacencies(bool hasBP)
{
  mfem::ParMesh* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);

  const int GRP_SZ = 25;
  char group_str[GRP_SZ];

  // TODO(JRC): Separate this out into group hierarchy setup and data allocation
  // stages like all of the other "createMeshBlueprint*" functions.
  SLIC_WARNING_IF(hasBP, "The case hasBP == true is not supported yet!");

  sidre::Group* adjset_grp = nullptr;
  if(pmesh->GetNGroups() > 1)
  {
    adjset_grp = m_bp_grp->createGroup("adjsets/mesh");
    adjset_grp->createViewString("association", "vertex");
    adjset_grp->createViewString("topology", "mesh");

    if(myid == 0)
    {
      m_bp_index_grp->createGroup("adjsets");
    }
  }

  int dim = pmesh->SpaceDimension();

  mfem::Array<int> tmp_verts;
  int tmp_edge;
  int tmp_face;
  int tmp_orientation;

  for(int gi = 1; gi < pmesh->GetNGroups(); ++gi)
  {
    // MFEM will build groups with zero shared vertices but still some shared
    // edges or faces, so we need to check all possible elements

    // Shared element views are not created if they would be empty though,
    // as this causes issues during a reload
    int num_gneighbors = pmesh->gtopo.GetGroupSize(gi);
    int num_gvertices = pmesh->GroupNVertices(gi);
    int num_gedges = pmesh->GroupNEdges(gi);
    int num_gtris = pmesh->GroupNTriangles(gi);
    int num_gquads = pmesh->GroupNQuadrilaterals(gi);

    int num_shared_elements = num_gvertices;
    if(dim >= 2)
    {
      num_shared_elements += num_gedges;
      if(dim >= 3)
      {
        num_shared_elements += num_gtris + num_gquads;
      }
    }
    const bool has_shared_elements = num_shared_elements > 0;

    if(has_shared_elements && (num_gneighbors > 1))
    {
      std::snprintf(group_str,
                    GRP_SZ,
                    "groups/g%d_%d",
                    pmesh->gtopo.GetGroupMasterRank(gi),
                    pmesh->gtopo.GetGroupMasterGroup(gi));
      sidre::Group* group_grp = adjset_grp->createGroup(group_str);

      sidre::View* gneighbors_view =
        group_grp->createViewAndAllocate("neighbors", sidre::INT_ID, num_gneighbors - 1);
      int* gneighbors_data = gneighbors_view->getData<int*>();

      // skip local domain when adding Blueprint neighbors
      const int* gneighbors = pmesh->gtopo.GetGroup(gi);
      for(int ni = 0, noff = 0; ni < num_gneighbors; ++ni)
      {
        if(gneighbors[ni] == 0)
        {
          noff++;
        }
        else
        {
          gneighbors_data[ni - noff] = pmesh->gtopo.GetNeighborRank(gneighbors[ni]);
        }
      }

      if(num_gvertices > 0)
      {
        sidre::View* gvertices_view =
          group_grp->createViewAndAllocate("values", sidre::INT_ID, num_gvertices);
        int* gvertices_data = gvertices_view->getData<int*>();

        for(int vi = 0; vi < num_gvertices; ++vi)
        {
          gvertices_data[vi] = pmesh->GroupVertex(gi, vi);
        }
      }

      // For all these higher-order elements, we store a flat list of tuples of
      // vertex indices instead of MFEM's internal edge/face indices for generality
      // This uses more space but technically the Blueprint does not standardize how
      // edge/face indices correspond to vertex indices
      if(dim >= 2)
      {
        // Don't create the group if there are no elements
        if(num_gedges > 0)
        {
          sidre::View* gedges_view =
            group_grp->createViewAndAllocate("edges", sidre::INT_ID, num_gedges * 2);
          int* gedges_data = gedges_view->getData<int*>();

          for(int ei = 0; ei < num_gedges; ++ei)
          {
            pmesh->GroupEdge(gi, ei, tmp_edge, tmp_orientation);
            pmesh->GetEdgeVertices(tmp_edge, tmp_verts);
            // Copy into array such that vertices for a given edge are contiguous
            // that is, v0 of e0, v1 of e0, v0 of e1, v1 of e1, v0 of e2, etc
            std::copy(tmp_verts.begin(), tmp_verts.end(), &gedges_data[2 * ei]);
          }
        }

        if(dim >= 3)
        {
          if(num_gtris > 0)
          {
            sidre::View* gtris_view =
              group_grp->createViewAndAllocate("triangles", sidre::INT_ID, num_gtris * 3);
            int* gtris_data = gtris_view->getData<int*>();
            for(int ti = 0; ti < num_gtris; ++ti)
            {
              pmesh->GroupTriangle(gi, ti, tmp_face, tmp_orientation);
              pmesh->GetFaceVertices(tmp_face, tmp_verts);
              std::copy(tmp_verts.begin(), tmp_verts.end(), &gtris_data[3 * ti]);
            }
          }

          if(num_gquads > 0)
          {
            sidre::View* gquads_view =
              group_grp->createViewAndAllocate("quadrilaterals", sidre::INT_ID, num_gquads * 4);
            int* gquads_data = gquads_view->getData<int*>();
            for(int qi = 0; qi < num_gquads; ++qi)
            {
              pmesh->GroupQuadrilateral(gi, qi, tmp_face, tmp_orientation);
              pmesh->GetFaceVertices(tmp_face, tmp_verts);
              std::copy(tmp_verts.begin(), tmp_verts.end(), &gquads_data[4 * qi]);
            }
          }
        }
      }
    }
  }
}
  #endif

// private method
bool MFEMSidreDataCollection::verifyMeshBlueprint()
{
  conduit::Node mesh_node;
  m_bp_grp->createNativeLayout(mesh_node);

  conduit::Node verify_info;
  bool result = conduit::blueprint::mesh::verify(mesh_node, verify_info);
  SLIC_WARNING_IF(!result,
                  "MFEMSidreDataCollection blueprint verification failed:\n"
                    << verify_info.to_yaml());
  return result;
}

bool MFEMSidreDataCollection::HasBoundaryMesh() const
{
  // check if this rank has any boundary elements
  int hasBndElts = mesh->GetNBE() > 0 ? 1 : 0;

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  // check if any rank has boundary elements
  mfem::ParMesh* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);
  if(pmesh)
  {
    int hasBndElts_g;
    MPI_Allreduce(&hasBndElts, &hasBndElts_g, 1, MPI_INT, MPI_MAX, pmesh->GetComm());

    hasBndElts = hasBndElts_g;
  }
  #endif

  return hasBndElts > 0 ? true : false;
}

void MFEMSidreDataCollection::SetMesh(Mesh* new_mesh)
{
  DataCollection::SetMesh(new_mesh);

  // hasBP is used to indicate if the data currently in the blueprint should be
  // used to replace the data in the mesh.
  bool hasBP = m_bp_grp->getNumViews() > 0 || m_bp_grp->getNumGroups() > 0;
  bool has_bnd_elts = HasBoundaryMesh();

  createMeshBlueprintStubs(hasBP);
  createMeshBlueprintState(hasBP);
  createMeshBlueprintCoordset(hasBP);

  GridFunction* nodes = new_mesh->GetNodes();

  // register the "mesh" topology in the blueprint.
  createMeshBlueprintTopologies(hasBP, s_mesh_topology_name);

  if(has_bnd_elts)
  {
    // Set the "boundary_topology" of "mesh" to "boundary".
    m_bp_grp->createViewString("topologies/" + s_mesh_topology_name + "/boundary_topology",
                               s_boundary_topology_name);

    // register the "boundary" topology in the blueprint.
    createMeshBlueprintTopologies(hasBP, s_boundary_topology_name);
  }

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  mfem::ParMesh* new_pmesh = dynamic_cast<mfem::ParMesh*>(new_mesh);
  m_comm = new_pmesh ? new_pmesh->GetComm() : MPI_COMM_NULL;
  if(new_pmesh)
  {
    createMeshBlueprintAdjacencies(hasBP);
  }
  #endif

  if(nodes)
  {
    // See the comment at the definition of 'hasBP' above.
    if(hasBP)
    {
      // Get the bp mesh nodes name.
      sidre::View* v_bp_nodes_name =
        m_bp_grp->getView("topologies/" + s_mesh_topology_name + "/grid_function");
      std::string bp_nodes_name(v_bp_nodes_name->getString());

      // Check that the names match, e.g. when loading the collection.
      SLIC_WARNING_IF(m_meshNodesGFName == bp_nodes_name,
                      "mismatch of requested and blueprint mesh nodes names");
      // Support renaming bp_nodes_name --> m_meshNodesGFName ?
    }

    if(m_owns_mesh_data)
    {
      // Make sure Sidre owns the data of the new_mesh's Nodes.
      if(!GetNamedBuffer(m_meshNodesGFName))
      {
        int sz = new_mesh->GetNodalFESpace()->GetVSize();
        double* gfData = AllocNamedBuffer(m_meshNodesGFName, sz)->getData();

        // See the comment at the definition of 'hasBP' above.
        if(!hasBP)
        {
          SLIC_ASSERT_MSG(nodes->Size() == sz, "");
          std::memcpy(gfData, nodes->GetData(), sizeof(double) * sz);
        }
      }
      // Since the named buffer for m_meshNodesGFName exists, the call to
      // RegisterField() below will replace the data of the nodes with the
      // data from the named buffer.
    }
    else
    {
      // Make sure Sidre does not have a named buffer for m_meshNodesGFName.
      SLIC_WARNING_IF(GetNamedBuffer(m_meshNodesGFName) != nullptr, "");
    }

    RegisterField(m_meshNodesGFName, nodes);

    if(own_data)
    {
      // Avoid double delete calls (for the nodes gf) when (own_data == true)
      // and the new_mesh owns its Nodes --> take ownership from new_mesh.
      // When new_mesh does not own its Nodes and (own_data == true), we can
      // not take ownership --> verify that does not happen.
      SLIC_WARNING_IF(!new_mesh->OwnsNodes(),
                      "mesh does not own its nodes, "
                      "can not take ownership");
      new_mesh->SetNodesOwner(false);
    }
  }
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
void MFEMSidreDataCollection::SetMesh(MPI_Comm comm, Mesh* new_mesh)
{
  // use MFEMSidreDataCollection's custom SetMesh, then set MPI info
  this->SetMesh(new_mesh);
  this->SetComm(comm);
}
  #endif

void MFEMSidreDataCollection::SetGroupPointers(Group* bp_index_grp, Group* domain_grp)
{
  SLIC_WARNING_IF(!domain_grp->hasGroup("blueprint"),
                  "Domain group does not contain a blueprint group.");

  m_bp_grp = domain_grp->getGroup("blueprint");
  m_bp_index_grp = bp_index_grp;
  m_named_bufs_grp = domain_grp->getGroup("named_buffers");
}

void MFEMSidreDataCollection::Load(const std::string& path, const std::string& protocol)
{
  DataCollection::DeleteAll();
  // Reset DataStore?

  int numInputRanks = 1;
  int numCommRanks = 1;

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(m_comm != MPI_COMM_NULL)
  {
    // The conduit abstraction appears to automatically handle the ".root"
    // suffix, but the IOManager does not, so it gets added here
    using axom::utilities::string::endsWith;
    std::string suffixedPath = endsWith(path, ".root") ? path : path + ".root";

    IOManager reader(m_comm);

    // Get the paths in which the global and domain groups are stored
    // This allows for a datastore structure other than the one used when this
    // class owns the DataStore
    SLIC_ERROR_IF(!m_bp_index_grp || !m_bp_grp, "Blueprint pointers must not be null");
    Group* global_grp = m_bp_index_grp->getParent()->getParent();
    Group* domain_grp = m_bp_grp->getParent();
    const std::string global_path = global_grp->getPath();
    const std::string domain_path = domain_grp->getPath();

    // Create a temporary group to read the data into
    // This is done instead of creating a temp DataStore so the buffers are intact
    Group* temp_root = m_bp_grp->getDataStore()->getRoot()->createGroup("_sidre_tmp_load");
    reader.read(temp_root, suffixedPath);

    // First transfer the global group from the temp group to its correct location
    Group* datastore_root = m_bp_grp->getDataStore()->getRoot();
    Group* new_global_group = temp_root->getGroup(global_grp->getPathName());
    Group* global_group_parent = datastore_root;
    if(!global_path.empty())
    {
      global_group_parent = datastore_root->getGroup(global_path);
    }

    // The group may already exist, but we're going to overwrite it
    global_group_parent->destroyGroup(new_global_group->getName());
    global_group_parent->moveGroup(new_global_group);

    // Then transfer the domain group in the same way
    Group* new_domain_group = temp_root->getGroup(domain_grp->getPathName());
    Group* domain_group_parent = datastore_root;
    if(!domain_path.empty())
    {
      domain_group_parent = datastore_root->getGroup(domain_path);
    }
    domain_group_parent->destroyGroup(new_domain_group->getName());
    domain_group_parent->moveGroup(new_domain_group);

    // Now that the data has been transferred, we can delete the temporary group
    temp_root->getParent()->destroyGroup("_sidre_tmp_load");

    // Get some data in support of error checks below
    numInputRanks = reader.getNumGroupsFromRoot(suffixedPath);
    MPI_Comm_size(m_comm, &numCommRanks);
  }
  else
  #endif
  {
    m_bp_grp->load(path, protocol);
  }

  // If the data collection created the datastore, it knows the layout of where
  // the domain and global groups are, and can restore them after the Load().
  //
  // If the data collection did not create the datastore, the host code must
  // reset these pointers after the load operation and also reset the state variables.
  if(m_owns_datastore)
  {
    // Use the same path format as was used to create the datastore
    SetGroupPointers(m_datastore_ptr->getRoot()->getGroup(name + "_global/blueprint_index/" + name),
                     m_datastore_ptr->getRoot()->getGroup(name));
    SLIC_ERROR_IF(m_bp_grp->getNumGroups() == 0,
                  fmt::format("Loaded datastore is empty, was the datastore created on a "
                              "different number of ranks? Current run has {} ranks, but "
                              "dataset was created using {} ranks",
                              numCommRanks,
                              numInputRanks));

    UpdateStateFromDS();
    UpdateMeshAndFieldsFromDS();
  }
}

void MFEMSidreDataCollection::LoadExternalData(const std::string& filename,
                                               const std::string& group_name)
{
  // Use the user-provided group name or the DataCollection's base group
  Group* grp = m_bp_grp->getDataStore()->getRoot();
  if(!group_name.empty())
  {
    SLIC_ERROR_IF(
      !m_bp_grp->hasGroup(group_name),
      axom::fmt::format("MFEMSidreDataCollection does not have a Sidre Group '{}'", group_name));
    grp = m_bp_grp->getGroup(group_name);
  }

  // Use the user-provided file name or the DataCollection's file name
  std::string path = name;
  if(!filename.empty())
  {
    path = filename;
  }
  path = get_file_path(path);

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(m_comm != MPI_COMM_NULL)
  {
    // The conduit abstraction appears to automatically handle the ".root"
    // suffix, but the IOManager does not, so it gets added here
    using axom::utilities::string::endsWith;
    std::string suffixedPath = endsWith(path, ".root") ? path : path + ".root";
    IOManager reader(m_comm);
    reader.loadExternalData(grp->getDataStore()->getRoot(), grp, suffixedPath);
  }
  else
  #endif
  {
    if(!group_name.empty())
    {
      grp->loadExternalData(path);
    }
    else
    {
      m_bp_grp->loadExternalData(path);
    }
  }
}

void MFEMSidreDataCollection::UpdateStateFromDS()
{
  SetTime(m_bp_grp->getView("state/time")->getData<double>());
  SetCycle(m_bp_grp->getView("state/cycle")->getData<int>());
  SetTimeStep(m_bp_grp->getView("state/time_step")->getData<double>());
}

void MFEMSidreDataCollection::UpdateStateToDS()
{
  SLIC_ASSERT_MSG(mesh != nullptr,
                  "Need to set mesh before updating state in MFEMSidreDataCollection.");

  m_bp_grp->getView("state/cycle")->setScalar(GetCycle());
  m_bp_grp->getView("state/time")->setScalar(GetTime());
  m_bp_grp->getView("state/time_step")->setScalar(GetTimeStep());

  if(myid == 0)
  {
    m_bp_index_grp->getView("state/cycle")->setScalar(GetCycle());
    m_bp_index_grp->getView("state/time")->setScalar(time);
  }
}

void MFEMSidreDataCollection::addMaterialSetToIndex()
{
  if(myid == 0)
  {
    for(auto it = m_matset_associations.begin(); it != m_matset_associations.end(); it++)
    {
      const std::string& mat_prefix = it->first;
      const std::string& matset_name = it->second;

      conduit::Node matnames;
      int nmats = 0;
      // Iterate over the fields in their group order.
      sidre::Group* fields_grp = m_bp_grp->getGroup("fields");
      for(auto& group : fields_grp->groups())
      {
        const std::string& field_name = group.getName();
        const auto tokens = utilities::string::rsplitN(field_name, 2, '_');
        // Expecting [base_field_name, material_id]
        if(tokens.size() != 2)
        {
          // The field_name did not split into 2 tokens matching
          // [base_field_name, material_id]. It's likely not a material
          // so we can skip the field.
          continue;
        }
        if(mat_prefix == tokens[0])
        {
          matnames[tokens[1]] = nmats++;
        }
      }

      // Now we know the material names for the current matset.
      if(nmats > 0)
      {
        Group* bp_matset_index_grp = m_bp_index_grp->createGroup("matsets/" + matset_name);
        bp_matset_index_grp->createViewString("topology", s_mesh_topology_name);
        Group* bp_materials_index_grp = bp_matset_index_grp->createGroup("materials");
        bp_materials_index_grp->importConduitTree(matnames);

        Group* bp_matset_group = m_bp_grp->getGroup("matsets/" + matset_name);
        bp_matset_index_grp->createViewString("path", bp_matset_group->getPathName());
      }
    }
  }
}

void MFEMSidreDataCollection::UpdateMeshAndFieldsFromDS()
{
  // 1. Start by constructing the mesh
  reconstructMesh();

  // 2. Set the mesh nodal grid function when present
  // QUESTION: Does this logic only apply when we own the datastore (i.e., when we know where things are)
  // or should we unconditionally run this logic whenever we Load()?
  const std::string gfPath = fmt::format("topologies/{}/grid_function", s_mesh_topology_name);
  if(GetBPGroup()->hasView(gfPath))
  {
    const std::string nodalGFName = GetBPGroup()->getView(gfPath)->getString();
    Group* fields_grp = m_bp_grp->getGroup("fields");
    if(fields_grp->hasGroup(nodalGFName))
    {
      reconstructField(fields_grp->getGroup(nodalGFName));
      if(this->HasField(nodalGFName))
      {
        this->SetMeshNodesName(nodalGFName);
        mesh->NewNodes(*this->GetField(nodalGFName), false);
      }
    }
  }

  // 3. Now we can finalize the mesh
  const bool fix_orientation = false;
  const bool refine = false;
  mesh->Finalize(refine, fix_orientation);

  // 4. And then finally reconstruct the remaining fields
  reconstructFields();
}

void MFEMSidreDataCollection::PrepareToSave()
{
  verifyMeshBlueprint();
  addMaterialSetToIndex();
  UpdateStateToDS();
}

void MFEMSidreDataCollection::Save()
{
  std::string filename = name;
  std::string protocol = Group::getDefaultIOProtocol();

  Save(filename, protocol);
}

void MFEMSidreDataCollection::Save(const std::string& filename, const std::string& protocol)
{
  PrepareToSave();

  std::string file_path = get_file_path(filename);
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  create_directory(file_path, mesh, myid);
  #else
  create_directory(prefix_path, mesh, myid);
  #endif

  sidre::Group* blueprint_indicies_grp = m_bp_index_grp->getParent();
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(m_comm != MPI_COMM_NULL)
  {
    IOManager writer(m_comm);

    // Create a shallow mirror of the DataStore structure that only includes
    // the data relevant to the calling DataCollection
    SLIC_ERROR_IF(!m_bp_index_grp || !m_bp_grp, "Blueprint pointers must not be null");
    Group* global_grp = m_bp_index_grp->getParent()->getParent();
    Group* domain_grp = m_bp_grp->getParent();
    const std::string global_path = global_grp->getPath();
    const std::string domain_path = domain_grp->getPath();
    Group* temp_root = m_bp_grp->getDataStore()->getRoot()->createGroup("_sidre_tmp_save");

    // First shallow copy the global group into the temporary group
    Group* temp_global_grp = temp_root;
    // We wouldn't need an extra conditional if hasGroup("") returned true
    if(!global_path.empty())
    {
      temp_global_grp = alloc_group(temp_root, global_path);
    }
    temp_global_grp->copyGroup(global_grp);

    // Then shallow copy the domain group in the same way
    Group* temp_domain_grp = temp_root;
    if(!domain_path.empty())
    {
      temp_domain_grp = alloc_group(temp_root, domain_path);
    }
    temp_domain_grp->copyGroup(domain_grp);

    int num_files = num_procs;
    if(m_num_files > 0)
    {
      SLIC_ASSERT_MSG(num_files <= num_procs,
                      "Save output must have num_files less than or equal to "
                      "number of mpi ranks");
      num_files = m_num_files;
      num_files = axom::utilities::clampUpper(num_files, num_procs);
    }

    // Specifying the tree pattern helps visit display our data.
    std::string tree_pattern = "datagroup";
    if(num_files != num_procs)
    {
      tree_pattern = "datagroup_%07d";
    }
    writer.write(temp_root, num_files, file_path, protocol, tree_pattern);

    // Now that we've written the data, we can delete the temporary group
    temp_root->getParent()->destroyGroup("_sidre_tmp_save");

    if(myid == 0)
    {
      // 'state/number_of_domains` may have been previously set to an incorrect
      // value before this object had access to the MPI communicator. Here
      // we reset to the correct value, or create a View to hold the value if
      // one does not already exist.
      if(m_bp_index_grp->hasView("state/number_of_domains"))
      {
        View* num_domains = m_bp_index_grp->getView("state/number_of_domains")->setScalar(num_procs);

        SLIC_ASSERT_MSG(num_domains,
                        "Failed to reset View 'state/number_of_domains' "
                        "in blueprint index to correct number of domains.");
        AXOM_UNUSED_VAR(num_domains);
      }
      else
      {
        m_bp_index_grp->createViewScalar("state/number_of_domains", num_procs);
      }

      if(protocol == "sidre_hdf5")
      {
        writer.writeGroupToRootFile(blueprint_indicies_grp, file_path + ".root");
      }
      // Root file support only available in hdf5.
      else
      {
        blueprint_indicies_grp->save(file_path + ".root", protocol);
      }
    }
  }
  else
  #endif
  {
    // If serial, use sidre group writer.
    m_bp_grp->save(file_path, protocol);

    blueprint_indicies_grp->save(file_path + ".root", protocol);
  }
}

// private method
void MFEMSidreDataCollection::addScalarBasedField(const std::string& field_name,
                                                  mfem::Vector* field,
                                                  const std::string& buffer_name,
                                                  IndexType offset,
                                                  const int num_dofs)
{
  sidre::Group* grp = m_bp_grp->getGroup("fields/" + field_name);
  SLIC_ASSERT_MSG(grp != nullptr, "field " << field_name << " does not exist");

  if(field->GetData() == nullptr)
  {
    AllocNamedBuffer(buffer_name, offset + num_dofs);
    // field->data is set below.
  }

  /*
   *  Mesh blueprint for a scalar-based grid function is of the form
   *    /fields/field_name/basis
   *              -- string value is GridFunction's FEC::Name
   *    /fields/field_name/values
   *              -- array of size num_dofs
   */

  // Make sure we have the View "values".
  sidre::View* vv = alloc_view(grp, "values");

  // Describe and apply the "values" View.
  // If the data store has buffer for field_name (e.g. AllocNamedBuffer was
  // called, or it was loaded from file), use that buffer.
  if(named_buffers_grp()->hasView(buffer_name))
  {
    sidre::View* bv = named_buffers_grp()->getView(buffer_name);
    SLIC_ASSERT_MSG(bv->hasBuffer() && bv->isDescribed(), "");

    // named buffers always have offset 0
    SLIC_ASSERT_MSG(bv->getSchema().dtype().offset() == 0, "");
    SLIC_ASSERT_MSG(bv->getNumElements() >= offset + num_dofs, "");

    if(vv->isEmpty())
    {
      vv->attachBuffer(bv->getBuffer())->apply(sidre::DOUBLE_ID, num_dofs, offset);
    }

    field->NewDataAndSize(vv->getData(), num_dofs);
  }
  else
  {
    // If we are not managing the grid function's data,
    // create a view with the external data
    vv->setExternalDataPtr(sidre::DOUBLE_ID, num_dofs, field->GetData());
  }
  SLIC_ASSERT_MSG(
    (num_dofs > 0 && vv->isApplied()) || (num_dofs == 0 && vv->isEmpty() && vv->isDescribed()),
    "invalid View state");
  SLIC_ASSERT_MSG(num_dofs == 0 || vv->getData() == field->GetData(),
                  "View data is different from GridFunction data");
  SLIC_ASSERT_MSG(vv->getNumElements() == num_dofs, "View size is different from GridFunction size");
}

// private method
void MFEMSidreDataCollection::addVectorBasedGridFunction(const std::string& field_name,
                                                         GridFunction* gf,
                                                         const std::string& buffer_name,
                                                         IndexType offset)
{
  sidre::Group* grp = m_bp_grp->getGroup("fields/" + field_name);
  SLIC_ASSERT_MSG(grp != nullptr, "field " << field_name << " does not exist");

  const int FLD_SZ = 20;
  char fidxName[FLD_SZ];

  int vdim = gf->FESpace()->GetVDim();
  int ndof = gf->FESpace()->GetNDofs();
  Ordering::Type ordering = gf->FESpace()->GetOrdering();

  if(gf->GetData() == nullptr)
  {
    AllocNamedBuffer(buffer_name, offset + vdim * ndof);
    // gf->data is set below.
  }

  /*
   *  Mesh blueprint for a vector-based grid function is of the form
   *    /fields/field_name/basis
   *              -- string value is GridFunction's FEC::Name
   *    /fields/field_name/values/x0
   *    /fields/field_name/values/x1
   *    ...
   *    /fields/field_name/values/xn
   *              -- each coordinate is an array of size ndof
   */

  // Get/create the Group "values".
  sidre::Group* vg = alloc_group(grp, "values");

  // Create the Views "x0", "x1", etc inside the "values" Group, vg.
  // If we have a named buffer for field_name, attach it to the Views;
  // otherwise set the Views to use gf->GetData() as external data.
  sidre::DataType dtype = sidre::DataType::c_double(ndof);
  const int entry_stride = (ordering == Ordering::byNODES ? 1 : vdim);
  const int vdim_stride = (ordering == Ordering::byNODES ? ndof : 1);
  dtype.set_stride(dtype.stride() * entry_stride);

  if(named_buffers_grp()->hasView(buffer_name))
  {
    sidre::View* bv = named_buffers_grp()->getView(buffer_name);
    SLIC_ASSERT_MSG(bv->hasBuffer() && bv->isDescribed(), "");

    // named buffers always have offset 0
    SLIC_ASSERT_MSG(bv->getSchema().dtype().offset() == 0, "");
    dtype.set_offset(dtype.element_bytes() * offset);

    for(int d = 0; d < vdim; d++)
    {
      std::snprintf(fidxName, FLD_SZ, "x%d", d);
      sidre::View* xv = alloc_view(vg, fidxName, dtype);
      xv->attachBuffer(bv->getBuffer());
      dtype.set_offset(dtype.offset() + dtype.element_bytes() * vdim_stride);
    }

    gf->NewDataAndSize(bv->getData<double*>() + offset, vdim * ndof);
  }
  else
  {
    for(int d = 0; d < vdim; d++)
    {
      std::snprintf(fidxName, FLD_SZ, "x%d", d);
      sidre::View* xv = alloc_view(vg, fidxName, dtype);
      xv->setExternalDataPtr(gf->GetData());
      dtype.set_offset(dtype.offset() + dtype.element_bytes() * vdim_stride);
    }
  }

  #ifdef AXOM_DEBUG
  for(int d = 0; d < vdim; d++)
  {
    std::snprintf(fidxName, FLD_SZ, "x%d", d);
    sidre::View* xv = vg->getView(fidxName);
    SLIC_ASSERT_MSG(
      (ndof > 0 && xv->isApplied()) || (ndof == 0 && xv->isEmpty() && xv->isDescribed()),
      "invalid View state");
    SLIC_ASSERT_MSG(ndof == 0 || xv->getData() == gf->GetData() + d * vdim_stride,
                    "View data is different from GridFunction data");
    SLIC_ASSERT_MSG(xv->getNumElements() == ndof, "View size is different from GridFunction size");
  }
  #endif
}

// private method
// Should only be called on mpi rank 0 ( or if serial problem ).
void MFEMSidreDataCollection::RegisterFieldInBPIndex(const std::string& field_name,
                                                     const int number_of_components)
{
  sidre::Group* bp_field_grp = m_bp_grp->getGroup("fields/" + field_name);
  sidre::Group* bp_index_field_grp = m_bp_index_grp->createGroup("fields/" + field_name);

  bp_index_field_grp->createViewString("path", bp_field_grp->getPathName());
  bp_index_field_grp->copyView(bp_field_grp->getView("topology"));
  bp_index_field_grp->copyView(bp_field_grp->getView("basis"));

  bp_index_field_grp->createViewScalar("number_of_components", number_of_components);
}

// private method
// Should only be called on mpi rank 0 ( or if serial problem ).
void MFEMSidreDataCollection::DeregisterFieldInBPIndex(const std::string& field_name)
{
  sidre::Group* fields_grp = m_bp_index_grp->getGroup("fields");
  SLIC_WARNING_IF(!fields_grp->hasGroup(field_name),
                  "No field exists in blueprint index with name " << field_name);

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  fields_grp->destroyGroup(field_name);
}

void MFEMSidreDataCollection::RegisterField(const std::string& field_name,
                                            GridFunction* gf,
                                            const std::string& buffer_name,
                                            IndexType offset)
{
  #ifdef AXOM_DEBUG
  SLIC_WARNING_IF(field_name.empty(), "Name for GridFunction was empty");
  SLIC_WARNING_IF(buffer_name.empty(), "Name for GridFunction destination buffer was empty");

  SLIC_WARNING_IF(gf == nullptr || gf->FESpace() == nullptr,
                  "Field with the name '"
                    << field_name << "' was provided a null GridFunction, so nothing was done.");
  #endif
  if(field_name.empty() || buffer_name.empty() || gf == nullptr || gf->FESpace() == nullptr)
  {
    return;
  }

  SLIC_ASSERT_MSG(mesh != nullptr,
                  "Need to set mesh before registering attributes in "
                  "MFEMSidreDataCollection.");

  // Register field_name in the blueprint group.
  sidre::Group* f = m_bp_grp->getGroup("fields");

  if(f->hasGroup(field_name))
  {
    // There are two possibilities:
    // 1. If HasField(field_name) is true - we are overwriting a field that
    //    was previously registered.
    // 2. Otherwise, the field was loaded from a file, or defined outside of
    //    the data collection.
    if(HasField(field_name))
    {
  #ifdef AXOM_DEBUG
      // Warn about overwriting field.
      // Skip warning when re-registering the nodal grid function
      if(field_name != m_meshNodesGFName)
      {
        SLIC_WARNING("field with the name '" << field_name
                                             << "' is already "
                                                "registered, overwriting the old field");
      }
  #endif
      DeregisterField(field_name);
    }
  }

  // This will return the existing field (if external), otherwise, a new group
  sidre::Group* grp = alloc_group(f, field_name);

  // Set the "basis" string using the gf's finite element space, overwrite if
  // necessary.
  sidre::View* v = alloc_view(grp, "basis");
  v->setString(gf->FESpace()->FEColl()->Name());

  // Set the topology of the GridFunction.
  // This is always 'mesh' except for a special case with the boundary material
  // attributes field.
  v = alloc_view(grp, "topology")->setString("mesh");

  SLIC_ASSERT_MSG(gf->Size() == gf->FESpace()->GetVSize(),
                  "GridFunction size does not match FiniteElementSpace size");

  // Set the data views of the grid function
  // e.g. the number of coefficients per DoF -- either scalar-valued or
  // vector-valued
  bool const isScalarValued = (gf->FESpace()->GetVDim() == 1);
  if(isScalarValued)
  {
    // Set the View "<m_bp_grp>/fields/<field_name>/values"
    addScalarBasedField(field_name, gf, buffer_name, offset, gf->FESpace()->GetVSize());
  }
  else  // vector valued
  {
    // Set the Group "<m_bp_grp>/fields/<field_name>/values"
    addVectorBasedGridFunction(field_name, gf, buffer_name, offset);
  }

  // Register field_name in the blueprint_index group.
  if(myid == 0)
  {
    // Note: The bp index requires GridFunction::VectorDim()
    //       since the GF might be scalar valued and have a vector basis
    //       (e.g. hdiv and hcurl spaces)
    RegisterFieldInBPIndex(field_name, gf->VectorDim());
  }

  // Check if the field contains material metadata
  checkForMaterialSet(field_name);
  checkForSpeciesSet(field_name);
  checkForMaterialDependentField(field_name);

  // Register field_name + gf in field_map.
  DataCollection::RegisterField(field_name, gf);
}

void MFEMSidreDataCollection::RegisterQField(const std::string& field_name,
                                             mfem::QuadratureFunction* qf,
                                             const std::string& buffer_name,
                                             axom::sidre::IndexType offset)
{
  #ifdef AXOM_DEBUG
  SLIC_WARNING_IF(field_name.empty(), "Name for QuadratureFunction was empty");
  SLIC_WARNING_IF(buffer_name.empty(), "Name for QuadratureFunction destination buffer was empty");

  SLIC_WARNING_IF(qf == nullptr || qf->GetSpace() == nullptr,
                  "QField with the name '"
                    << field_name
                    << "' was provided a null QuadratureFunction, so nothing was done.");
  #endif

  if(field_name.empty() || buffer_name.empty() || qf == nullptr || qf->GetSpace() == nullptr)
  {
    return;
  }

  // Register field_name in the blueprint group.
  sidre::Group* f = m_bp_grp->getGroup("fields");

  if(f->hasGroup(field_name))
  {
    // There are two possibilities:
    // 1. If HasQField(field_name) is true - we are overwriting a field that
    //    was previously registered.
    // 2. Otherwise, the field was loaded from a file, or defined outside of
    //    the data collection.
    if(HasQField(field_name))
    {
  #ifdef AXOM_DEBUG
      // Warn about overwriting field.
      // Skip warning when re-registering the nodal grid function
      if(field_name != m_meshNodesGFName)
      {
        SLIC_WARNING("QField with the name '" << field_name
                                              << "' is already "
                                                 "registered, overwriting the old qfield");
      }
  #endif
      DeregisterQField(field_name);
    }
  }

  // This will return the existing field (if external), otherwise, a new group
  sidre::Group* grp = alloc_group(f, field_name);

  // A QField has the following schema:
  // <m_bp_grp>/fields/<field_name>/
  // <m_bp_grp>/fields/<field_name>/topology     = "mesh"
  // <m_bp_grp>/fields/<field_name>/values       = array of size QuadratureFunction::Size()
  // <m_bp_grp>/fields/<field_name>/basis        = string that encodes the QF
  //                                               integration order and vector dimension
  //                                               "QF_Default_[ORDER]_[VDIM]"

  // Set the "basis" string using the qf's order and vdim, overwrite if necessary
  sidre::View* v = alloc_view(grp, "basis");
  // In the future, we should be able to have mfem::QuadratureFunction provide us this name.
  const std::string basis_name =
    fmt::format("QF_Default_{0}_{1}", qf->GetSpace()->GetOrder(), qf->GetVDim());
  v->setString(basis_name);

  // Set the topology of the QuadratureFunction.
  // This is always 'mesh'
  v = alloc_view(grp, "topology")->setString("mesh");

  // Set the View "<m_bp_grp>/fields/<field_name>/values"
  addScalarBasedField(field_name, qf, buffer_name, offset, qf->Size());

  // Register field_name in the blueprint_index group.
  if(myid == 0)
  {
    RegisterFieldInBPIndex(field_name, qf->GetVDim());
  }

  // Register field_name + qf in field_map.
  mfem::DataCollection::RegisterQField(field_name, qf);
}

void MFEMSidreDataCollection::DeregisterQField(const std::string& field_name)
{
  // Deregister field_name from field_map.
  mfem::DataCollection::DeregisterQField(field_name);
  removeField(field_name);
}

void MFEMSidreDataCollection::RegisterAttributeField(const std::string& attr_name, bool is_bdry)
{
  SLIC_ASSERT_MSG(mesh != nullptr,
                  "Need to set mesh before registering attributes in "
                  "MFEMSidreDataCollection.");

  // Register attr_name in the blueprint group.
  sidre::Group* f = m_bp_grp->getGroup("fields");
  if(f->hasGroup(attr_name))
  {
    bool isAttr = attr_map.Has(attr_name);
    bool isFld = field_map.Has(attr_name);

    if(isAttr)
    {
      SLIC_WARNING("field with the name '"
                   << attr_name
                   << "' is already "
                      " registered as an attribute, overwriting old values.");
      DeregisterAttributeField(attr_name);
    }
    else if(isFld)
    {
      SLIC_WARNING("field with the name '"
                   << attr_name
                   << "' is already "
                      " registered as a field, skipping register attribute.");
      return;
    }
  }

  // Generate sidre views and groups for this mesh attribute and allocate space
  addIntegerAttributeField(attr_name, is_bdry);

  if(myid == 0)
  {
    RegisterAttributeFieldInBPIndex(attr_name);
  }

  // Register new attribute array with attr_map
  sidre::View* a = m_bp_grp->getGroup("fields")->getGroup(attr_name)->getView("values");
  mfem::Array<int>* attr =
    new mfem::Array<int>(a->getData<int*>(), static_cast<int>(a->getNumElements()));

  attr_map.Register(attr_name, attr, true);
}

void MFEMSidreDataCollection::RegisterAttributeFieldInBPIndex(const std::string& attr_name)
{
  SLIC_ASSERT_MSG(m_bp_grp->getGroup("fields") != nullptr,
                  "Mesh blueprint does not have 'fields' group");
  SLIC_ASSERT_MSG(m_bp_index_grp->getGroup("fields") != nullptr,
                  "Mesh blueprint index does not have 'fields' group");

  // get the BP attr group
  sidre::Group* attr_grp = m_bp_grp->getGroup("fields")->getGroup(attr_name);

  // create blueprint index for this attribute
  sidre::Group* bp_index_attr_grp = m_bp_index_grp->getGroup("fields")->createGroup(attr_name);

  bp_index_attr_grp->createViewString("path", attr_grp->getPathName());
  bp_index_attr_grp->copyView(attr_grp->getView("association"));
  bp_index_attr_grp->copyView(attr_grp->getView("topology"));
  bp_index_attr_grp->createViewScalar("number_of_components", 1);
}

void MFEMSidreDataCollection::DeregisterAttributeField(const std::string& attr_name)
{
  attr_map.Deregister(name, true);

  sidre::Group* attr_grp = m_bp_grp->getGroup("fields");
  SLIC_WARNING_IF(!attr_grp->hasGroup(attr_name),
                  "No field exists in blueprint with name " << attr_name);

  // Delete attr_name from the blueprint group.

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  attr_grp->destroyGroup(attr_name);

  // Delete field_name from the blueprint_index group.
  if(myid == 0)
  {
    DeregisterAttributeFieldInBPIndex(attr_name);
  }

  // Delete field_name from the named_buffers group, if allocated.
  FreeNamedBuffer(attr_name);
}

void MFEMSidreDataCollection::DeregisterAttributeFieldInBPIndex(const std::string& attr_name)
{
  sidre::Group* fields_grp = m_bp_index_grp->getGroup("fields");
  SLIC_WARNING_IF(!fields_grp->hasGroup(attr_name),
                  "No attribute exists in blueprint index with name " << attr_name);

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  fields_grp->destroyGroup(attr_name);
}

void MFEMSidreDataCollection::addIntegerAttributeField(const std::string& attr_name, bool is_bdry)
{
  sidre::Group* fld_grp = m_bp_grp->getGroup("fields");
  SLIC_ASSERT_MSG(fld_grp != nullptr, "'fields' group does not exist");

  const int num_elem = is_bdry ? mesh->GetNBE() : mesh->GetNE();
  std::string topo_name = is_bdry ? "boundary" : "mesh";

  sidre::Group* attr_grp = fld_grp->createGroup(attr_name);
  attr_grp->createViewString("association", "element");
  attr_grp->createViewAndAllocate("values", sidre::INT_ID, num_elem);
  attr_grp->createViewString("topology", topo_name);
}

void MFEMSidreDataCollection::removeField(const std::string& field_name)
{
  sidre::Group* fields_grp = m_bp_grp->getGroup("fields");
  SLIC_WARNING_IF(!fields_grp->hasGroup(field_name),
                  "No field exists in blueprint with name " << field_name);

  // Delete field_name from the blueprint group.

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  fields_grp->destroyGroup(field_name);

  // Delete field_name from the blueprint_index group.
  if(myid == 0)
  {
    DeregisterFieldInBPIndex(field_name);
  }

  // Delete field_name from the named_buffers group, if allocated.
  FreeNamedBuffer(field_name);
}

void MFEMSidreDataCollection::DeregisterField(const std::string& field_name)
{
  // Deregister field_name from field_map.
  DataCollection::DeregisterField(field_name);
  removeField(field_name);
}

void MFEMSidreDataCollection::AssociateMaterialSet(const std::string& volume_fraction_field_name,
                                                   const std::string& matset_name)
{
  auto iter = m_matset_associations.find(volume_fraction_field_name);
  if(iter != m_matset_associations.end())
  {
    SLIC_WARNING("Volume fraction field "
                 << volume_fraction_field_name
                 << " has already been associated with a material set: " << iter->second);
    return;
  }
  m_matset_associations[volume_fraction_field_name] = matset_name;
  Group* matset_group = m_bp_grp->createGroup("matsets/" + matset_name);
  // Since we're creating the matset, associate it with a topology
  // These will always be associated with the mesh
  matset_group->createViewString("topology", s_mesh_topology_name);
}

void MFEMSidreDataCollection::AssociateSpeciesSet(const std::string& species_field_name,
                                                  const std::string& specset_name,
                                                  const std::string& matset_name,
                                                  const bool volume_dependent)
{
  SLIC_WARNING_IF(!m_bp_grp->hasGroup("matsets/" + matset_name),
                  "The material set '" << matset_name << "' has not been associated with a field");
  auto iter = m_specset_associations.find(species_field_name);
  if(iter != m_specset_associations.end())
  {
    SLIC_WARNING("Species field " << species_field_name
                                  << " has already been associated with a species set: "
                                  << iter->second);
    return;
  }
  m_specset_associations[species_field_name] = specset_name;
  Group* specset_grp = m_bp_grp->createGroup("specsets/" + specset_name);
  specset_grp->createViewScalar("volume_dependent", static_cast<std::int8_t>(volume_dependent));
  // Since we're creating the species set, associate it with a material set
  specset_grp->createViewString("matset", matset_name);
}

void MFEMSidreDataCollection::AssociateMaterialDependentField(
  const std::string& material_dependent_field_name,
  const std::string& matset_name)
{
  SLIC_WARNING_IF(!m_bp_grp->hasGroup("matsets/" + matset_name),
                  "The material set '" << matset_name << "' has not been associated with a field");
  auto iter = m_material_dependent_fields.find(material_dependent_field_name);
  if(iter != m_material_dependent_fields.end())
  {
    SLIC_WARNING("Field " << material_dependent_field_name
                          << " has already been labeled as material-dependent"
                             " and associated with a material set: "
                          << iter->second);
    return;
  }
  m_material_dependent_fields[material_dependent_field_name] = matset_name;
}

View* MFEMSidreDataCollection::getFieldValuesView(const std::string& field_name)
{
  const std::string field_values_name = "fields/" + field_name + "/values";
  View* values_view = nullptr;
  if(m_bp_grp->hasView(field_values_name))
  {
    // Scalar-valued field
    values_view = m_bp_grp->getView(field_values_name);
  }
  else if(m_bp_grp->hasGroup(field_values_name))
  {
    // Vector-valued field
    values_view = m_bp_grp->getGroup(field_values_name)->getView("x0");
  }
  SLIC_WARNING_IF(values_view == nullptr, "Field " << field_name << " was not registered");
  return values_view;
}

void MFEMSidreDataCollection::checkForMaterialSet(const std::string& field_name)
{
  const auto tokens = utilities::string::rsplitN(field_name, 2, '_');
  // Expecting [base_field_name, material_id]
  if(tokens.size() != 2)
  {
    return;
  }

  auto iter = m_matset_associations.find(tokens[0]);
  // If it hasn't been registered as a matset, it's not a volume fraction field
  if(iter == m_matset_associations.end())
  {
    return;
  }
  const std::string matset_name = iter->second;

  View* vol_fractions_view = getFieldValuesView(field_name);

  Group* fractions_group = alloc_group(m_bp_grp, "matsets/" + matset_name + "/volume_fractions");

  View* matset_frac_view = fractions_group->copyView(vol_fractions_view);
  matset_frac_view->rename(tokens[1]);

  // NOTE: The index group is populated with matsets prior to saving once we
  //       know all of the volume fraction fields.
}

void MFEMSidreDataCollection::checkForSpeciesSet(const std::string& field_name)
{
  const auto tokens = utilities::string::rsplitN(field_name, 3, '_');
  // Expecting [base_field_name, material_id, component]
  if(tokens.size() != 3)
  {
    return;
  }

  auto iter = m_specset_associations.find(tokens[0]);
  // If it hasn't been registered as a matset, it's not a species field
  if(iter == m_specset_associations.end())
  {
    return;
  }
  const std::string specset_name = iter->second;

  View* species_values_view = getFieldValuesView(field_name);

  Group* specset_material_group =
    alloc_group(m_bp_grp, "specsets/" + specset_name + "/matset_values/" + tokens[1]);
  View* specset_values = specset_material_group->copyView(species_values_view);
  specset_values->rename(tokens[2]);
  // FIXME: Do we need to add anything to the index group?
}

void MFEMSidreDataCollection::checkForMaterialDependentField(const std::string& field_name)
{
  const auto tokens = utilities::string::rsplitN(field_name, 2, '_');
  // Expecting [base_field_name, material_id]
  if(tokens.size() != 2)
  {
    return;
  }

  auto iter = m_material_dependent_fields.find(tokens[0]);
  // If it hasn't been registered, it's not a material-dependent field
  if(iter == m_material_dependent_fields.end())
  {
    return;
  }
  const std::string matset_name = iter->second;

  View* material_values = getFieldValuesView(field_name);

  Group* field_grp = m_bp_grp->getGroup("fields/" + tokens[0]);
  if(!field_grp->hasView("matset"))
  {
    field_grp->createViewString("matset", matset_name);
  }

  Group* values_grp = alloc_group(field_grp, "matset_values");
  View* added_values = values_grp->copyView(material_values);
  added_values->rename(tokens[1]);
  // FIXME: Do we need to add anything to the index group?
}

// private method
std::string MFEMSidreDataCollection::getElementName(mfem::Element::Type elementEnum)
{
  // Note -- the mapping from Element::Type to string is based on
  //   enum Element::Type { POINT, SEGMENT, TRIANGLE, QUADRILATERAL,
  //                        TETRAHEDRON, HEXAHEDRON };
  // Note: -- the string names are from conduit's blueprint

  switch(elementEnum)
  {
  case mfem::Element::POINT:
    return "point";
  case mfem::Element::SEGMENT:
    return "line";
  case mfem::Element::TRIANGLE:
    return "tri";
  case mfem::Element::QUADRILATERAL:
    return "quad";
  case mfem::Element::TETRAHEDRON:
    return "tet";
  case mfem::Element::HEXAHEDRON:
    return "hex";
  case mfem::Element::WEDGE:
  default:;
  }

  return "unknown";
}

mfem::Geometry::Type MFEMSidreDataCollection::getElementTypeFromName(const std::string& name)
{
  if(name == "point")
  {
    return mfem::Geometry::POINT;
  }
  else if(name == "line")
  {
    return mfem::Geometry::SEGMENT;
  }
  else if(name == "tri")
  {
    return mfem::Geometry::TRIANGLE;
  }
  else if(name == "quad")
  {
    return mfem::Geometry::SQUARE;
  }
  else if(name == "tet")
  {
    return mfem::Geometry::TETRAHEDRON;
  }
  else if(name == "hex")
  {
    return mfem::Geometry::CUBE;
  }
  else
  {
    return mfem::Geometry::INVALID;
  }
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
namespace detail
{
/**
 * @brief Wrapper class to enable the reconstruction of a ParMesh from a
 * Conduit blueprint
 *
 * @note This is only needed because mfem::ParMesh doesn't provide a
 * constructor that allows us to set all the required information.
 */
class SidreParMeshWrapper : public mfem::ParMesh
{
public:
  /**
   * @brief Constructs a new parallel mesh, intended to be assigned
   * to an mfem::ParMesh*
   *
   * @param [in] bp_grp The root group of the Conduit Blueprint structure
   * @param [in] comm The MPI communicator to use with the new mesh
   * @param [in] mesh_topo_name The name of the mesh topology in the Blueprint structure
   *
   * The remaining parameters are used to construct the mfem::Mesh base subobject
   */
  SidreParMeshWrapper(Group* bp_grp,
                      MPI_Comm comm,
                      const std::string& mesh_topo_name,
                      double* vertices,
                      int num_vertices,
                      int* element_indices,
                      mfem::Geometry::Type element_type,
                      int* element_attributes,
                      int num_elements,
                      int* boundary_indices,
                      mfem::Geometry::Type boundary_type,
                      int* boundary_attributes,
                      int num_boundary_elements,
                      int dimension,
                      int space_dimension = -1)
  {
    MyComm = comm;
    gtopo.SetComm(comm);
    ConstructMeshSubObject(vertices,
                           num_vertices,
                           element_indices,
                           element_type,
                           element_attributes,
                           num_elements,
                           boundary_indices,
                           boundary_type,
                           boundary_attributes,
                           num_boundary_elements,
                           dimension,
                           space_dimension);

    SLIC_ERROR_IF(!bp_grp->hasGroup("adjsets"),
                  "Cannot reconstruct a ParMesh without adjacency sets");

    MPI_Comm_size(MyComm, &NRanks);
    MPI_Comm_rank(MyComm, &MyRank);

    // NOTE: The remaining logic is borrowed heavily from the mfem::ParMesh
    // constructors, since it's setting up the internal data structures used
    // to keep track of communication groups and shared elements

    ReduceMeshGen();  // determine the global 'meshgen'

    // FIXME: The word "group" is overloaded in this context
    // Can we do any better with variable naming?
    // A group can refer to a communication group or a Sidre group

    auto mesh_adjset_groups = bp_grp->getGroup("adjsets/" + mesh_topo_name + "/groups");
    auto num_groups = mesh_adjset_groups->getNumGroups();

    auto shared_geoms = GetSharedGeometries(mesh_adjset_groups);

    // Set up the gtopo object
    InitializeGroupTopology(shared_geoms);

    std::size_t total_shared_vertices = 0;
    std::size_t total_shared_edges = 0;
    std::size_t total_shared_faces = 0;

    for(const auto& shared_geom : shared_geoms)
    {
      total_shared_vertices += shared_geom.shared_verts.size();
      total_shared_edges += shared_geom.shared_edges.size();
      total_shared_faces += shared_geom.shared_triangles.size();
      total_shared_faces += shared_geom.shared_quadrilaterals.size();
    }

    // Resizing for shared vertex data
    svert_lvert.SetSize(total_shared_vertices);
    group_svert.SetDims(num_groups, total_shared_vertices);

    if(dimension >= 2)
    {
      // Resizing for shared edge data
      sedge_ledge.SetSize(total_shared_edges);
      shared_edges.SetSize(total_shared_edges);
      group_sedge.SetDims(num_groups, total_shared_edges);
    }
    if(dimension >= 3)
    {
      // Resizing for shared face data
      sface_lface.SetSize(total_shared_faces);
      group_stria.MakeI(num_groups);
      group_squad.MakeI(num_groups);
    }
    else
    {
      group_stria.SetSize(num_groups, 0);  // create empty group_stria
      group_squad.SetSize(num_groups, 0);  // create empty group_squad
    }

    std::size_t comm_group_idx = 1;  // The ID of the current group
    std::size_t total_verts_added = 0;
    std::size_t total_edges_added = 0;

    for(const auto& shared_geom : shared_geoms)
    {
      // Add shared vertices
      auto n_shared_verts = shared_geom.shared_verts.size();
      n_shared_verts += total_verts_added;  // Get the index into the global table
      SLIC_ERROR_IF(n_shared_verts > group_svert.Size_of_connections(),
                    "incorrect number of total_shared_vertices");
      group_svert.GetI()[comm_group_idx] = n_shared_verts;

      for(const auto vert : shared_geom.shared_verts)
      {
        group_svert.GetJ()[total_verts_added] = total_verts_added;
        svert_lvert[total_verts_added] = *vert;
        total_verts_added++;
      }

      if(dimension >= 2)
      {
        auto n_shared_edges = shared_geom.shared_edges.size();
        n_shared_edges += total_edges_added;
        SLIC_ERROR_IF(n_shared_edges > group_sedge.Size_of_connections(),
                      "incorrect number of total_shared_edges");
        group_sedge.GetI()[comm_group_idx] = n_shared_edges;

        for(const auto edge : shared_geom.shared_edges)
        {
          group_sedge.GetJ()[total_edges_added] = total_edges_added;
          shared_edges[total_edges_added] = new mfem::Segment(edge, 1);
          total_edges_added++;
        }

        if(dimension >= 3)
        {
          // Add all the shared faces
          // VectorSpanIterator allows for readable iteration and
          // "registering" of the shared faces
          for(const auto triangle : shared_geom.shared_triangles)
          {
            shared_trias.Append({triangle[0], triangle[1], triangle[2]});
          }

          for(const auto quad : shared_geom.shared_quadrilaterals)
          {
            shared_quads.Append({quad[0], quad[1], quad[2], quad[3]});
          }

          group_stria.AddColumnsInRow(comm_group_idx - 1, shared_geom.shared_triangles.size());
          group_squad.AddColumnsInRow(comm_group_idx - 1, shared_geom.shared_quadrilaterals.size());
        }
      }

      comm_group_idx++;
    }

    if(dimension >= 3)
    {
      SLIC_ERROR_IF(shared_trias.Size() + shared_quads.Size() != sface_lface.Size(),
                    "incorrect number of total_shared_faces");
      // Define the J arrays of group_stria and group_squad -- they just contain
      // consecutive numbers starting from 0 up to shared_trias.Size()-1 and
      // shared_quads.Size()-1, respectively.
      group_stria.MakeJ();
      for(int i = 0; i < shared_trias.Size(); i++)
      {
        group_stria.GetJ()[i] = i;
      }
      group_squad.MakeJ();
      for(int i = 0; i < shared_quads.Size(); i++)
      {
        group_squad.GetJ()[i] = i;
      }
    }

    // Delay the finalization of the mesh until we can set the nodal GridFunction
  }

private:
  /**
   * @brief Clone of mfem::Mesh constructor that allows data fields
   * to be set directly
   *
   * Needed because the mfem::Mesh base subobject cannot be constructed
   * directly, so we need to set the data fields "manually"
   *
   * @see mfem::Mesh::Mesh(double, int, int*, ...)
   */
  void ConstructMeshSubObject(double* _vertices,
                              int num_vertices,
                              int* element_indices,
                              mfem::Geometry::Type element_type,
                              int* element_attributes,
                              int num_elements,
                              int* boundary_indices,
                              mfem::Geometry::Type boundary_type,
                              int* boundary_attributes,
                              int num_boundary_elements,
                              int dimension,
                              int space_dimension)
  {
    // NOTE: This is all copied directly from the mfem::Mesh constructor
    // It needs to be called directly though because we have to "manually"
    // initialize the indirect base subobject - fortunately all relevant
    // variables are "protected"
    if(space_dimension == -1)
    {
      space_dimension = dimension;
    }

    InitMesh(dimension,
             space_dimension,
             /*num_vertices*/ 0,
             num_elements,
             num_boundary_elements);

    int element_index_stride = mfem::Geometry::NumVerts[element_type];
    int boundary_index_stride =
      num_boundary_elements > 0 ? mfem::Geometry::NumVerts[boundary_type] : 0;

    static_assert(std::is_standard_layout<mfem::Vertex>::value,
                  "mfem::Vertex must have standard layout for reinterpret_casting");
    vertices.MakeRef(reinterpret_cast<mfem::Vertex*>(_vertices), num_vertices);
    NumOfVertices = num_vertices;

    for(int i = 0; i < num_elements; i++)
    {
      elements[i] = NewElement(element_type);
      elements[i]->SetVertices(element_indices + i * element_index_stride);
      elements[i]->SetAttribute(element_attributes[i]);
    }
    NumOfElements = num_elements;

    for(int i = 0; i < num_boundary_elements; i++)
    {
      boundary[i] = NewElement(boundary_type);
      boundary[i]->SetVertices(boundary_indices + i * boundary_index_stride);
      boundary[i]->SetAttribute(boundary_attributes[i]);
    }
    NumOfBdrElements = num_boundary_elements;

    FinalizeTopology();
  }

  /**
   * @brief A span over a list of vectors arranged contiguously (not interleaved)
   * Allows for convenient iteration over the list without having to manage
   * sizes and offsets in multiple places
   *
   * \note This is a convenience wrapper over Sidre's "stride" data attribute,
   * perhaps this would be useful as a user-facing utility??
   */
  // TODO: Replace with std::span when C++20 is available
  class VectorSpan
  {
  public:
    /**
       * @brief Constructs a new VectorSpan
       * @param [in] data The data to build a view over
       * @param [in] num_scalars The length of the data array
       * @param [in] num_components The number of components in each vector
       */
    VectorSpan(const int* data, const IndexType num_scalars, const IndexType num_components)
      : m_data(data)
      , m_num_vectors(num_scalars / num_components)
      , m_num_components(num_components)
    {
      SLIC_ERROR_IF(num_scalars % num_components != 0,
                    "VectorSpan number of components does not evenly divide "
                    "length of array");
    }

    VectorSpan() = default;

    /**
     * @brief Returns the number of vectors referenced in the span
     */
    int size() const { return m_num_vectors; }

    /**
     * @brief Helper class to iterate over the vector span
     */
    class VectorSpanIterator
    {
    public:
      VectorSpanIterator(const int* start, const IndexType stride) : m_ptr(start), m_stride(stride)
      { }

      /**
       * @brief Compares two iterators
       */
      bool operator!=(const VectorSpanIterator& other)
      {
        return (m_ptr != other.m_ptr) || (m_stride != other.m_stride);
      }

      /**
       * @brief Advances the iterator to the next vector in the list
       */
      void operator++()
      {
        // Moves the pointer forward by the stride
        m_ptr += m_stride;
      }

      /**
       * @brief Returns the pointer to the start of the current
       * vector - accesses on this pointer are not bounds-checked
       */
      const int* operator*() const { return m_ptr; }

    private:
      const int* m_ptr;
      const IndexType m_stride;
    };

    VectorSpanIterator begin() const { return {m_data, m_num_components}; }
    VectorSpanIterator end() const
    {
      return {m_data + (m_num_vectors * m_num_components), m_num_components};
    }

  private:
    // Use reasonable defaults to avoid extra logic for empty spans
    const int* m_data = nullptr;
    IndexType m_num_vectors = 0;
    IndexType m_num_components = 0;
  };

  /**
   * @brief A data class containing all relevant information needed on geometric
   * communication groups within an mfem::ParMesh - the shared geometric
   * entities and neighboring groups
   */
  struct SharedGeometries
  {
    VectorSpan shared_verts;      // Indices of shared vertices
    VectorSpan shared_edges;      // Pairs of vertex indices corresponding to shared edges
    VectorSpan shared_triangles;  // 3-tuples of vertex indices corresponding to shared triangular faces
    VectorSpan shared_quadrilaterals;  // 4-tuples of vertex indices corresponding to shared quadrilateral faces
    mfem::Array<int> neighbors;  // Rank IDs of the neighboring nodes/processes
  };

  /**
   * @brief Initializes the mfem::GroupTopology object
   * @param [in] shared_geoms The shared geometry information
   * @pre The shared vertex and neighbor fields of each SharedGeometries must be populated
   * @post The gtopo member of the base subobject is initialized
   */
  void InitializeGroupTopology(const std::vector<SharedGeometries>& shared_geoms)
  {
    mfem::ListOfIntegerSets comm_group_integer_sets;

    // The first group always contains only the current rank
    mfem::IntegerSet first_set;
    mfem::Array<int>& first_array = first_set;  // Bind a ref so we can modify it
    first_array.Append(MyRank);
    comm_group_integer_sets.Insert(first_set);

    for(const auto& shared_geom : shared_geoms)
    {
      mfem::IntegerSet integer_set;
      mfem::Array<int>& array = integer_set;  // Bind a ref so we can modify it
      // On a given rank, the group ID corresponding to the current rank
      // is a member of each group
      array.Append(MyRank);
      array.Append(shared_geom.neighbors);
      array.Sort();  // MFEM requires that the sets be sorted
      comm_group_integer_sets.Insert(integer_set);
    }
    // FIXME: 822 or 823?
    // This appears to just be an MPI tag used by MFEM when it creates the
    // group topology object in ParMesh and elsewhere - see
    // https://github.com/mfem/mfem/blob/b0770915511fc63fda3120b13c88db48946e4302/mesh/pmesh.cpp#L206-L207
    // https://github.com/mfem/mfem/blob/b0770915511fc63fda3120b13c88db48946e4302/general/communication.cpp#L324
    gtopo.Create(comm_group_integer_sets, 823);
  }

  /**
   * @brief Retrieves the shared geometry information from the mesh topology
   * adjacency set group
   *
   * @param [in] mesh_adjset_groups The blueprint group corresponding to the
   * adjacency sets, i.e. <blueprint_root>/adjsets/<mesh_topology>/groups
   *
   * @return The set of SharedGeometries, one for each communication group
   */
  std::vector<SharedGeometries> GetSharedGeometries(const Group* mesh_adjset_groups)
  {
    std::vector<SharedGeometries> shared_geoms(mesh_adjset_groups->getNumGroups());

    // Iterate over both the *sidre* groups group
    // to fill in the shared geometry data
    int group_idx = 0;
    for(const Group& grp : mesh_adjset_groups->groups())
    {
      const Group* adjset_grp = &grp;
      auto& shared_geom = shared_geoms[group_idx++];

      // Copy the neighbors array
      const View* group_neighbors = adjset_grp->getView("neighbors");
      shared_geom.neighbors.Append(group_neighbors->getData(), group_neighbors->getNumElements());

      // This group's shared vertices
      if(adjset_grp->hasView("values"))
      {
        auto verts = adjset_grp->getView("values");
        shared_geom.shared_verts = {verts->getData(), verts->getNumElements(), 1};
      }

      // This group's shared edges
      if(adjset_grp->hasView("edges"))
      {
        auto edges = adjset_grp->getView("edges");
        shared_geom.shared_edges = {edges->getData(), edges->getNumElements(), 2};
      }

      // This group's shared triangular faces
      if(adjset_grp->hasView("triangles"))
      {
        auto tris = adjset_grp->getView("triangles");
        shared_geom.shared_triangles = {tris->getData(), tris->getNumElements(), 3};
      }

      // This group's shared quadrilateral faces
      if(adjset_grp->hasView("quadrilaterals"))
      {
        auto quads = adjset_grp->getView("quadrilaterals");
        shared_geom.shared_quadrilaterals = {quads->getData(), quads->getNumElements(), 4};
      }
    }

    return shared_geoms;
  }
};

} /* namespace detail */
  #endif  // defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)

// private method
void MFEMSidreDataCollection::reconstructMesh()
{
  // mfem::ConduitDataCollection::BlueprintMeshToMesh would be useful here, but
  // we need all the parameters to construct a ParMesh (the Mesh base subobject
  // is initialized manually)

  SLIC_ERROR_IF(!verifyMeshBlueprint(),
                "Cannot reconstruct mesh, data does not satisfy Conduit Blueprint");

  SLIC_ERROR_IF(!m_bp_grp->hasView("coordsets/coords/values/x"),
                "Cannot reconstruct a mesh without a Cartesian coordinate set");

  View* vertex_view = m_bp_grp->getView("coordsets/coords/values/x");
  SLIC_ERROR_IF(vertex_view->isExternal(),
                "Cannot reconstruct a mesh that was built using external data");

  // Assumes that Vertex is a standard layout type with only an array member
  static_assert(std::is_standard_layout<mfem::Vertex>::value,
                "mfem::Vertex must have standard layout to determine size of "
                "its array member");
  const std::size_t EXPECTED_STRIDE = sizeof(mfem::Vertex) / sizeof(double);
  if((vertex_view->getTypeID() != DataTypeId::DOUBLE_ID) ||
     (vertex_view->getStride() != EXPECTED_STRIDE))
  {
    SLIC_ERROR("Vertex array must consist of interleaved doubles");
  }
  // Sufficient to always grab the underlying data to the x-coords
  // as the layout of the internal buffer is identical to what MFEM
  // expects, i.e., x0 y0 z0 x1 y1 z1 regardless of dimension
  double* vertices = vertex_view->getData();

  // Use the x to get the number of vertices
  int num_vertices = vertex_view->getNumElements();

  const std::string mesh_elements_path = "topologies/" + s_mesh_topology_name + "/elements";

  SLIC_ERROR_IF(!m_bp_grp->hasGroup(mesh_elements_path),
                "Cannot reconstruct mesh without mesh topology");

  int* element_indices = m_bp_grp->getView(mesh_elements_path + "/connectivity")->getData<int*>();

  // Name of the element type - convert to mfem::Geometry::Type later
  std::string element_name = m_bp_grp->getView(mesh_elements_path + "/shape")->getString();

  const std::string element_attribute_path = "fields/" + s_mesh_topology_name + s_attribute_suffix;
  View* element_attribute_view = m_bp_grp->getView(element_attribute_path + "/values");

  int* element_attributes = element_attribute_view->getData<int*>();
  int num_elements = static_cast<int>(element_attribute_view->getNumElements());

  const std::string boundary_elements_path = "topologies/" + s_boundary_topology_name + "/elements";

  SLIC_ERROR_IF(!m_bp_grp->hasGroup(boundary_elements_path),
                "Cannot reconstruct mesh without boundary topology");

  int* boundary_indices =
    m_bp_grp->getView(boundary_elements_path + "/connectivity")->getData<int*>();

  // Name of the element type - convert to mfem::Geometry::Type later
  std::string bdr_element_name = m_bp_grp->getView(boundary_elements_path + "/shape")->getString();

  const std::string boundary_attribute_path =
    "fields/" + s_boundary_topology_name + s_attribute_suffix;
  View* bdr_attribute_view = m_bp_grp->getView(boundary_attribute_path + "/values");

  int* boundary_attributes = bdr_attribute_view->getData<int*>();
  int num_boundary_elements = bdr_attribute_view->getNumElements();

  const std::string coordset_path = "coordsets/" + s_coordset_name;

  int dimension = 1;

  if(m_bp_grp->hasView(coordset_path + "/values/z"))
  {
    dimension = 3;
  }

  else if(m_bp_grp->hasView(coordset_path + "/values/y"))
  {
    dimension = 2;
  }

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  SLIC_ERROR_IF(m_comm == MPI_COMM_NULL,
                "Must set the communicator with SetComm before a ParMesh can "
                "be reconstructed");
  // If it has an adjacencies group, the reloaded state was a ParMesh
  if(m_bp_grp->hasGroup("adjsets/" + s_mesh_topology_name))
  {
    m_owned_mesh = std::unique_ptr<detail::SidreParMeshWrapper>(
      new detail::SidreParMeshWrapper(m_bp_grp,
                                      m_comm,
                                      s_mesh_topology_name,
                                      vertices,
                                      num_vertices,
                                      element_indices,
                                      getElementTypeFromName(element_name),
                                      element_attributes,
                                      num_elements,
                                      boundary_indices,
                                      getElementTypeFromName(bdr_element_name),
                                      boundary_attributes,
                                      num_boundary_elements,
                                      dimension));
  }
  else
  #endif
  {
    m_owned_mesh =
      std::unique_ptr<mfem::Mesh>(new mfem::Mesh(vertices,
                                                 num_vertices,
                                                 element_indices,
                                                 getElementTypeFromName(element_name),
                                                 element_attributes,
                                                 num_elements,
                                                 boundary_indices,
                                                 getElementTypeFromName(bdr_element_name),
                                                 boundary_attributes,
                                                 num_boundary_elements,
                                                 dimension));
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    // If this is a vacuously parallel run (MPI enabled, but only one rank),
    // we don't know whether the original mesh was a ParMesh.  In case the user
    // is expecting a ParMesh, we create a trivial ParMesh.  Since an
    // mfem::ParMesh is-an mfem::Mesh, this won't affect users who don't need it
    if(num_procs == 1)
    {
      m_owned_mesh = std::unique_ptr<mfem::ParMesh>(new mfem::ParMesh(m_comm, *m_owned_mesh));
    }
  #endif
  }
  // Now that we've initialized an owning pointer, set the base subobject's
  // mesh pointer as a non-owning pointer
  mesh = m_owned_mesh.get();
}

void MFEMSidreDataCollection::reconstructField(Group* field_grp)
{
  // Filter out the non-user-registered attribute fields
  if(!field_grp->hasView("association"))
  {
    int vdim = 1;  // The default
    // GridFunction vs QuadratureFunction
    bool is_gridfunc = true;
    // FiniteElementCollection/QuadratureSpace construction
    const std::string basis_name = field_grp->getView("basis")->getString();
    // Check if it's a QuadratureFunction
    if(basis_name.find("QF") == 0 && (m_quadspaces.count(basis_name) == 0))
    {
      // The vdim is being overwritten here with the value in the basis string so we
      // can correctly construct the QuadratureFunction
      m_quadspaces[basis_name] =
        std::unique_ptr<mfem::QuadratureSpace>(detail::NewQuadratureSpace(basis_name, mesh, vdim));
      is_gridfunc = false;
    }
    // Only need to create a new FEColl if one doesn't already exist
    else if(m_fecolls.count(basis_name) == 0)
    {
      m_fecolls[basis_name] = std::unique_ptr<mfem::FiniteElementCollection>(
        mfem::FiniteElementCollection::New(basis_name.c_str()));
    }

    View* value_view = nullptr;
    auto ordering = mfem::Ordering::byNODES;  // The default
    // Scalar grid function
    if(field_grp->hasView("values"))
    {
      value_view = field_grp->getView("values");
    }

    // Vector grid function
    else if(field_grp->hasGroup("values"))
    {
      // Sufficient to use address of first component as data is interleaved
      value_view = field_grp->getGroup("values")->getView("x0");
      vdim = field_grp->getGroup("values")->getNumViews();
      if(value_view->getStride() == static_cast<axom::sidre::IndexType>(vdim))
      {
        ordering = mfem::Ordering::byVDIM;
      }
    }
    else
    {
      SLIC_ERROR("Cannot reconstruct grid function - field values not found");
    }

    // We cache the FESpaces to avoid reconstructing them when not needed
    // An FESpace is uniquely identified by the basis of its FEColl, ordering, and vdim number
    const std::string fespace_id =
      axom::fmt::format("{0}_{1}_{2}",
                        basis_name,
                        (ordering == mfem::Ordering::byVDIM) ? "vdim" : "nodes",
                        vdim);

    // Only need to create a new FESpace if one doesn't already exist
    if(is_gridfunc && (m_fespaces.count(fespace_id) == 0))
    {
      // FiniteElementSpace - mesh ptr and FEColl ptr
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      auto* parmesh = dynamic_cast<mfem::ParMesh*>(mesh);
      if(parmesh)
      {
        m_fespaces[fespace_id] = std::unique_ptr<mfem::ParFiniteElementSpace>(
          new mfem::ParFiniteElementSpace(parmesh, m_fecolls.at(basis_name).get(), vdim, ordering));
      }
      else
  #endif
      {
        m_fespaces[fespace_id] = std::unique_ptr<mfem::FiniteElementSpace>(
          new mfem::FiniteElementSpace(mesh, m_fecolls.at(basis_name).get(), vdim, ordering));
      }
    }

    double* values = value_view->getData();

    if(is_gridfunc)
    {
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
      auto parfes = dynamic_cast<mfem::ParFiniteElementSpace*>(m_fespaces.at(fespace_id).get());
      if(parfes)
      {
        m_owned_gridfuncs.emplace_back(new mfem::ParGridFunction(parfes, values));
      }
      else
  #endif
      {
        m_owned_gridfuncs.emplace_back(
          new mfem::GridFunction(m_fespaces.at(fespace_id).get(), values));
      }

      // Register a non-owning pointer with the base subobject
      DataCollection::RegisterField(field_grp->getName(), m_owned_gridfuncs.back().get());
    }
    else
    {
      m_owned_quadfuncs.emplace_back(
        new mfem::QuadratureFunction(m_quadspaces.at(basis_name).get(), vdim));

      // FIXME: A bug was fixed in QuadratureFunction::SetSpace after mfem@4.5
      // See: https://github.com/mfem/mfem/pull/3281
      // Update using the following when possible to avoid an unnecessary allocation in the constructor
      // (1) Use default Quadrature function constructor
      // (2) Use qf->SetSpace(qs, values, vdim);
      auto* qf = m_owned_quadfuncs.back().get();
      qf->NewDataAndSize(values, vdim * qf->GetSpace()->GetSize());

      // Register a non-owning pointer with the base subobject
      DataCollection::RegisterQField(field_grp->getName(), qf);
    }
  }
}

void MFEMSidreDataCollection::reconstructFields()
{
  sidre::Group* f = m_bp_grp->getGroup("fields");
  if(f != nullptr)
  {
    for(auto& field_grp : f->groups())
    {
      // The nodal grid function will already have been reconstructed
      if(field_grp.getName() != m_meshNodesGFName)
      {
        reconstructField(&field_grp);
      }
    }
  }
}

}  // namespace sidre
}  // namespace axom

#endif  // AXOM_USE_MFEM
