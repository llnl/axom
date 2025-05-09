// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "axom/mint/mesh/blueprint.hpp"

// Axom includes
#include "axom/core/Types.hpp"  // for nullptr

// Mint includes
#include "axom/mint/config.hpp"          // for AXOM_MINT_USE_SIDRE
#include "axom/mint/mesh/MeshTypes.hpp"  // for mesh types

// Slic includes
#include "axom/slic/interface/slic.hpp"  // for SLIC macros

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre/core/Group.hpp"  // for sidre::Group
  #include "axom/sidre/core/View.hpp"   // for sidre::View
#endif

namespace axom
{
namespace mint
{
namespace blueprint
{
#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
bool isValidRootGroup(const sidre::Group* group)
{
  if(group == nullptr)
  {
    SLIC_WARNING("supplied group is NULL!");
    return false;
  }

  const bool hasCoordsets = group->hasChildGroup("coordsets");
  const bool hasTopologies = group->hasChildGroup("topologies");
  const bool hasFields = group->hasChildGroup("fields");

  SLIC_WARNING_IF(!hasCoordsets,
                  "sidre::Group " << group->getPathName() << " is missing coordsets group!");
  SLIC_WARNING_IF(!hasTopologies,
                  "sidre::Group " << group->getPathName() << " is missing topologies group!");
  SLIC_WARNING_IF(!hasFields, "sidre::Group " << group->getPathName() << " is missing fields group!");

  return ((hasCoordsets && hasTopologies && hasFields));
}

//------------------------------------------------------------------------------
bool isValidTopologyGroup(const sidre::Group* topo)
{
  if(topo == nullptr)
  {
    SLIC_WARNING("supplied topology group is NULL!");
    return false;
  }

  const std::string path = topo->getPathName();

  const bool hasTypeView = topo->hasChildView("type");
  SLIC_WARNING_IF(!hasTypeView, "[" << path << "] is missing 'type' view!");

  const bool isTypeAString = (hasTypeView) ? topo->getView("type")->isString() : false;
  SLIC_WARNING_IF(!isTypeAString, "'type' view in [" << path << "] is not a string");

  const bool hasCoordset = topo->hasChildView("coordset");
  SLIC_WARNING_IF(!hasCoordset, "[" << path << "] is missing 'coordset' view!");

  const bool isCoordsetAString = (hasCoordset) ? topo->getView("coordset")->isString() : false;
  SLIC_WARNING_IF(!isCoordsetAString, "'coordset' view in [" << path << "] is not a string");

  const bool status = hasTypeView && hasCoordset && isTypeAString && isCoordsetAString;

  return (status);
}

//------------------------------------------------------------------------------
bool isValidCoordsetGroup(const sidre::Group* coordset)
{
  if(coordset == nullptr)
  {
    SLIC_WARNING("supplied coordset group is NULL!");
    return false;
  }

  const std::string path = coordset->getPathName();

  const bool hasTypeView = coordset->hasChildView("type");
  SLIC_WARNING_IF(!hasTypeView, "[" << path << "] is missing 'type' view!");

  const bool isTypeAString = (hasTypeView) ? coordset->getView("type")->isString() : false;
  SLIC_WARNING_IF(!isTypeAString, "'type' view in [" << path << "] is not a string");

  return (hasTypeView && isTypeAString);
}

//------------------------------------------------------------------------------
const sidre::Group* getCoordsetGroup(const sidre::Group* group, const sidre::Group* topology)
{
  SLIC_ERROR_IF(!blueprint::isValidRootGroup(group),
                "supplied group does not conform to the blueprint!");
  SLIC_ERROR_IF(topology == nullptr, "supplied topology group is null!");
  SLIC_ERROR_IF(!blueprint::isValidTopologyGroup(topology),
                "supplied topology group does not conform to the blueprint!");

  const sidre::Group* coordsets = group->getGroup("coordsets");

  const char* coordset_name = topology->getView("coordset")->getString();
  SLIC_WARNING_IF(!coordsets->hasChildGroup(coordset_name),
                  "cannot find coordset [" << coordset_name << "] in " << coordsets->getPathName());

  const sidre::Group* coords = coordsets->getGroup(coordset_name);
  SLIC_WARNING_IF(coords == nullptr,
                  "null coordset [" << coordset_name << "] in " << coordsets->getPathName());

  return coords;
}

//------------------------------------------------------------------------------
const sidre::Group* getCoordsetGroup(const sidre::Group* group, const std::string& coords)
{
  SLIC_ERROR_IF(!blueprint::isValidRootGroup(group),
                "supplied group does not conform to the blueprint!");

  const sidre::Group* coordsets = group->getGroup("coordsets");
  const std::string path = coordsets->getPathName();

  // get coordset group
  const sidre::Group* coordset = nullptr;
  if(coords.empty())
  {
    SLIC_ERROR_IF(coordsets->getNumGroups() == 0, "[" << coordsets->getPathName() << "] is empty!");
    SLIC_WARNING_IF(coordsets->getNumGroups() > 1, "multiple coordsets found!  ");
    coordset = coordsets->getGroup(0);
  }
  else
  {
    SLIC_ERROR_IF(!coordsets->hasChildGroup(coords),
                  "[" << path << "] is missing requested coordset group [" << coords << "]");

    coordset = coordsets->getGroup(coords);
  }

  return (coordset);
}

//------------------------------------------------------------------------------
const sidre::Group* getTopologyGroup(const sidre::Group* group, const std::string& topo)
{
  SLIC_ERROR_IF(!blueprint::isValidRootGroup(group),
                "supplied group does not conform to the blueprint!");

  const sidre::Group* topologies = group->getGroup("topologies");
  const std::string path = topologies->getPathName();

  // get topology group
  const sidre::Group* topology = nullptr;
  if(topo.empty())
  {
    SLIC_ERROR_IF(topologies->getNumGroups() == 0, "[" << topologies->getPathName() << "] is empty!");
    SLIC_WARNING_IF(topologies->getNumGroups() > 1, "multiple topologies found!  ");
    topology = topologies->getGroup(0);
  }
  else
  {
    SLIC_ERROR_IF(!topologies->hasChildGroup(topo),
                  "[" << path << "] is missing requested topology group [" << topo << "]");

    topology = topologies->getGroup(topo);
  }

  return (topology);
}

//------------------------------------------------------------------------------
void initializeTopologyGroup(sidre::Group* group,
                             const std::string& topo,
                             const std::string& coordset,
                             const std::string& type)
{
  SLIC_ASSERT(group != nullptr);
  sidre::Group* topo_group = group->getGroup("topologies");
  SLIC_ASSERT(topo_group != nullptr);
  sidre::Group* cur_topo = topo_group->getGroup(topo);
  SLIC_ASSERT(cur_topo != nullptr);

  cur_topo->createView("type")->setString(type);
  cur_topo->createView("coordset")->setString(coordset);
}

//------------------------------------------------------------------------------
void getMeshTypeAndDimension(int& mesh_type,
                             int& dimension,
                             const sidre::Group* group,
                             const std::string& topo)
{
  SLIC_ERROR_IF(!blueprint::isValidRootGroup(group),
                "supplied group does not conform to the blueprint!");

  const sidre::Group* topology = blueprint::getTopologyGroup(group, topo);
  SLIC_ERROR_IF(!blueprint::isValidTopologyGroup(topology),
                "mesh topology does not conform to the blueprint!");

  const sidre::Group* coords = blueprint::getCoordsetGroup(group, topology);
  SLIC_ERROR_IF(!blueprint::isValidCoordsetGroup(coords),
                "mesh coordset does not conform to the blueprint!");

  // get topology type
  const char* topo_type = topology->getView("type")->getString();
  SLIC_ASSERT(topo_type != nullptr);

  // detect mesh type based on the topology type
  if(strcmp(topo_type, "uniform") == 0)
  {
    SLIC_ERROR_IF(
      !coords->hasChildGroup("origin"),
      "missing [origin] group from [" << coords->getPathName() << "], required for a uniform mesh");

    mesh_type = STRUCTURED_UNIFORM_MESH;
    dimension = coords->getGroup("origin")->getNumViews();

  }  // END if UNIFORM MESH
  else if(strcmp(topo_type, "rectilinear") == 0)
  {
    SLIC_ERROR_IF(!coords->hasChildGroup("values"),
                  "missing [values] group from [" << coords->getPathName()
                                                  << "], required for a rectilinear mesh");

    mesh_type = STRUCTURED_RECTILINEAR_MESH;
    dimension = coords->getGroup("values")->getNumViews();

  }  // END if RECTILINEAR_MESH
  else if(strcmp(topo_type, "structured") == 0)
  {
    SLIC_ERROR_IF(!coords->hasChildGroup("values"),
                  "missing [values] group from [" << coords->getPathName()
                                                  << "], required for a structured mesh");

    mesh_type = STRUCTURED_CURVILINEAR_MESH;
    dimension = coords->getGroup("values")->getNumViews();

  }  // END if STRUCTURED_MESH
  else if(strcmp(topo_type, "points") == 0)
  {
    SLIC_ERROR_IF(!coords->hasChildGroup("values"),
                  "missing [values] group from [" << coords->getPathName()
                                                  << "], required for a particle mesh");

    mesh_type = PARTICLE_MESH;
    dimension = coords->getGroup("values")->getNumViews();

  }  // END if PARTICLE_MESH
  else if(strcmp(topo_type, "unstructured") == 0)
  {
    SLIC_ERROR_IF(!coords->hasChildGroup("values"),
                  "missing [values] group from [" << coords->getPathName()
                                                  << "], required for a unstructured mesh");

    // check if this is a particle mesh stored as an unstructured mesh
    const char* shape = topology->getView("elements/shape")->getString();
    mesh_type = (strcmp(shape, "point") == 0) ? PARTICLE_MESH : UNSTRUCTURED_MESH;
    dimension = coords->getGroup("values")->getNumViews();

  }  // END if UNSTRUCTURED_MESH
  else
  {
    mesh_type = UNDEFINED_MESH;
    dimension = -1;
    SLIC_ERROR("invalid mesh topology_type=[" << topo_type << "] ");
  }
}

//------------------------------------------------------------------------------
bool hasMixedCellTypes(const sidre::Group* group, const std::string& topo)
{
  SLIC_ERROR_IF(!blueprint::isValidRootGroup(group),
                "supplied group does not conform to the blueprint!");

  const sidre::Group* topology = blueprint::getTopologyGroup(group, topo);
  SLIC_ERROR_IF(!blueprint::isValidTopologyGroup(topology),
                "mesh topology does not conform to the blueprint!");

  if(topology->getView("type")->getString() != std::string("unstructured"))
  {
    return false;
  }

  SLIC_ERROR_IF(!topology->hasChildGroup("elements"),
                "Unstructured topology has no 'elements' group.");

  const sidre::Group* elems_group = topology->getGroup("elements");

  SLIC_ERROR_IF(!elems_group->hasChildView("shape"), "elements group has no 'shape' view.");

  const sidre::View* shape_view = elems_group->getView("shape");
  SLIC_ERROR_IF(!shape_view->isString(), "'shape' view must hold a string.");

  if(shape_view->getString() == std::string("mixed"))
  {
    return true;
  }
  else
  {
    return false;
  }
}

//------------------------------------------------------------------------------
void getStructuredMeshProperties(int dimension,
                                 IndexType node_dims[3],
                                 int64 node_ext[6],
                                 const sidre::Group* coordset)
{
  SLIC_ERROR_IF(dimension < 1 || dimension > 3, "invalid dimension!");
  SLIC_ERROR_IF(node_dims == nullptr, "supplied extent is null!");
  SLIC_ERROR_IF(node_ext == nullptr, "supplied global extent is null!");
  SLIC_ERROR_IF(!blueprint::isValidCoordsetGroup(coordset), "invalid coordset group!");

  sidre::Group* c = const_cast<sidre::Group*>(coordset);

  const char* dim_names[] = {"dims/i", "dims/j", "dims/k"};
  const char* global_names[] = {"global_ext/i_min",
                                "global_ext/i_max",
                                "global_ext/j_min",
                                "global_ext/j_max",
                                "global_ext/k_min",
                                "global_ext/k_max"};

  for(int dim = 0; dim < dimension; ++dim)
  {
    node_dims[dim] = c->getView(dim_names[dim])->getScalar();
  }  // END for

  for(int i = 0; i < 6; ++i)
  {
    node_ext[i] = c->getView(global_names[i])->getScalar();
  }
}

//------------------------------------------------------------------------------
void setStructuredMeshProperties(int dimension,
                                 const IndexType node_dims[3],
                                 const int64 node_ext[6],
                                 sidre::Group* coordset)
{
  SLIC_ERROR_IF(dimension < 1 || dimension > 3, "invalid dimension!");
  SLIC_ERROR_IF(node_dims == nullptr, "supplied extent is null!");
  SLIC_ERROR_IF(node_ext == nullptr, "supplied global extent is null!");
  SLIC_ERROR_IF(coordset == nullptr, "invalid coordset group!");

  const char* dim_names[] = {"dims/i", "dims/j", "dims/k"};
  const char* global_names[] = {"global_ext/i_min",
                                "global_ext/i_max",
                                "global_ext/j_min",
                                "global_ext/j_max",
                                "global_ext/k_min",
                                "global_ext/k_max"};

  sidre::Group* c = const_cast<sidre::Group*>(coordset);

  for(int dim = 0; dim < dimension; ++dim)
  {
    coordset->createView(dim_names[dim])->setScalar(node_dims[dim]);
  }  // END for

  for(int i = 0; i < 6; ++i)
  {
    c->createView(global_names[i])->setScalar(node_ext[i]);
  }
}

void setExtent(sidre::Group* coordset, const int64 node_ext[6])
{
  SLIC_ERROR_IF(node_ext == nullptr, "supplied global extent is null!");
  SLIC_ERROR_IF(coordset == nullptr, "invalid coordset group!");

  const char* global_names[] = {"global_ext/i_min",
                                "global_ext/i_max",
                                "global_ext/j_min",
                                "global_ext/j_max",
                                "global_ext/k_min",
                                "global_ext/k_max"};

  for(int i = 0; i < 6; ++i)
  {
    coordset->getView(global_names[i])->setScalar(node_ext[i]);
  }
}

//------------------------------------------------------------------------------
void getUniformMeshProperties(int dimension, double* origin, double* spacing, const sidre::Group* coordset)
{
  SLIC_ERROR_IF(dimension < 1 || dimension > 3, "invalid dimension!");
  SLIC_ERROR_IF(origin == nullptr, "supplied null pointer for origin!");
  SLIC_ERROR_IF(spacing == nullptr, "supplied null pointer for spacing!");
  SLIC_ERROR_IF(!blueprint::isValidCoordsetGroup(coordset), "invalid coordset group!");

  sidre::Group* c = const_cast<sidre::Group*>(coordset);

  const char* origin_names[] = {"origin/x", "origin/y", "origin/z"};
  const char* spacing_names[] = {"spacing/dx", "spacing/dy", "spacing/dz"};

  SLIC_ERROR_IF(c->getView("type")->getString() != std::string("uniform"),
                "Mesh is not a UniformMesh.");
  for(int dim = 0; dim < dimension; ++dim)
  {
    origin[dim] = c->getView(origin_names[dim])->getScalar();
    spacing[dim] = c->getView(spacing_names[dim])->getScalar();
  }  // END for
}

//------------------------------------------------------------------------------
void setUniformMeshProperties(int dimension,
                              const double* origin,
                              const double* spacing,
                              sidre::Group* coordset)
{
  SLIC_ERROR_IF(dimension < 1 || dimension > 3, "invalid dimension!");
  SLIC_ERROR_IF(origin == nullptr, "supplied null pointer for origin!");
  SLIC_ERROR_IF(spacing == nullptr, "supplied null pointer for spacing!");

  const char* origin_names[] = {"origin/x", "origin/y", "origin/z"};
  const char* spacing_names[] = {"spacing/dx", "spacing/dy", "spacing/dz"};

  coordset->createView("type")->setString("uniform");
  for(int dim = 0; dim < dimension; ++dim)
  {
    coordset->createView(origin_names[dim])->setScalar(origin[dim]);
    coordset->createView(spacing_names[dim])->setScalar(spacing[dim]);
  }  // END for
}

#endif

} /* namespace blueprint */
} /* namespace mint */
} /* namespace axom */
