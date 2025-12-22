// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Associated header file
#include "Group.hpp"

#include "conduit_relay.hpp"

#ifdef AXOM_USE_HDF5
  #include "conduit_relay_io_hdf5.hpp"
#endif

#include "axom/core/ListCollection.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/MapCollection.hpp"
#include "axom/core/Path.hpp"

// Sidre headers
#include "Buffer.hpp"
#include "DataStore.hpp"

#include "axom/fmt.hpp"

namespace axom
{
namespace sidre
{
// Helper macro for defining a prepend string for sidre::Group log messages
// We are using it to add the pathName() of the group
#ifndef SIDRE_GROUP_LOG_PREPEND
  #define SIDRE_GROUP_LOG_PREPEND \
    "[Group: '" << this->getPathName() << "'" << (this->isRoot() ? " (root)" : "") << "] "
#endif

// Initialization of static path delimiter character for methods that
// support path syntax.
const char Group::s_path_delimiter = '/';

///////////////////////////////////////////////////////////////////////////////
//
// Private utility functions to cast ItemCollections to (named) MapCollections.
//
///////////////////////////////////////////////////////////////////////////////

MapCollection<View>* Group::getNamedViews()
{
  SLIC_ASSERT_MSG(this->isUsingMap(), "Invalid cast: The views in this group do not have names");

  return static_cast<MapCollection<View>*>(m_view_coll);
}

const MapCollection<View>* Group::getNamedViews() const
{
  SLIC_ASSERT_MSG(this->isUsingMap(), "Invalid cast: The views in this group do not have names");

  return static_cast<const MapCollection<View>*>(m_view_coll);
}

MapCollection<Group>* Group::getNamedGroups()
{
  SLIC_ASSERT_MSG(this->isUsingMap(), "Invalid cast: The groups in this group do not have names");

  return static_cast<MapCollection<Group>*>(m_group_coll);
}

const MapCollection<Group>* Group::getNamedGroups() const
{
  SLIC_ASSERT_MSG(this->isUsingMap(), "Invalid cast: The groups in this group do not have names");

  return static_cast<const MapCollection<Group>*>(m_group_coll);
}

////////////////////////////////////////////////////////////////////////
//
// Basic query and accessor methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return path of Group object, not including its name.
 *
 *************************************************************************
 */
std::string Group::getPath() const
{
  const Group* root = getDataStore()->getRoot();
  const Group* curr = getParent();
  std::string thePath = curr->getName();
  curr = curr->getParent();

  while(curr != root)
  {
    thePath = curr->getName() + s_path_delimiter + thePath;
    curr = curr->getParent();
  }

  return thePath;
}

/*
 *************************************************************************
 *
 * Insert information about data associated with Group subtree with this 
 * Group at root of tree (default 'recursive' is true), or for this Group 
 * only ('recursive' is false) in fields of given Conduit Node.
 *
 *************************************************************************
 */
void Group::getDataInfo(Node& n, bool recursive) const
{
  //
  // Initialize Node fields
  //
  IndexType num_groups = 0;
  IndexType num_views = 0;
  IndexType num_views_empty = 0;
  IndexType num_views_buffer = 0;
  IndexType num_views_external = 0;
  IndexType num_views_scalar = 0;
  IndexType num_views_string = 0;
  IndexType num_bytes_assoc_with_views = 0;
  IndexType num_bytes_external = 0;

  n["num_groups"] = num_groups;
  n["num_views"] = num_views;
  n["num_views_empty"] = num_views_empty;
  n["num_views_buffer"] = num_views_buffer;
  n["num_views_external"] = num_views_external;
  n["num_views_scalar"] = num_views_scalar;
  n["num_views_string"] = num_views_string;
  n["num_bytes_assoc_with_views"] = num_bytes_assoc_with_views;
  n["num_bytes_external"] = num_bytes_external;

  std::set<IndexType> buffer_ids;

  getDataInfoHelper(n, buffer_ids, recursive);

  const DataStore* ds = getDataStore();
  IndexType num_bytes_in_buffers = 0;
  for(auto it = buffer_ids.begin(); it != buffer_ids.end(); ++it)
  {
    num_bytes_in_buffers += ds->getBuffer(*it)->getTotalBytes();
  }
  n["num_bytes_in_buffers"] = num_bytes_in_buffers;
}

/*
 *************************************************************************
 *
 * Private helper method to support getDataInfo() method.
 *
 *************************************************************************
 */
void Group::getDataInfoHelper(Node& n, std::set<IndexType>& buffer_ids, bool recursive) const
{
  //
  // Grab Node entries for updating data info for this Group
  //
  IndexType num_groups = n["num_groups"].value();
  IndexType num_views = n["num_views"].value();
  IndexType num_views_empty = n["num_views_empty"].value();
  IndexType num_views_buffer = n["num_views_buffer"].value();
  IndexType num_views_external = n["num_views_external"].value();
  IndexType num_views_scalar = n["num_views_scalar"].value();
  IndexType num_views_string = n["num_views_string"].value();
  IndexType num_bytes_assoc_with_views = n["num_bytes_assoc_with_views"].value();
  IndexType num_bytes_external = n["num_bytes_external"].value();

  num_groups += 1;  // count this group
  num_views += getNumViews();

  //
  // Gather info from Views owned by this Group
  //
  for(auto& view : views())
  {
    if(view.isExternal())
    {
      num_views_external += 1;
      num_bytes_external += view.getTotalBytes();
    }
    else if(view.isScalar())
    {
      num_views_scalar += 1;
      num_bytes_assoc_with_views += view.getTotalBytes();
    }
    else if(view.isString())
    {
      num_views_string += 1;
      num_bytes_assoc_with_views += view.getTotalBytes();
    }
    else if(view.hasBuffer())
    {
      num_views_buffer += 1;
      const Buffer* buf = view.getBuffer();
      if(buf->isAllocated())
      {
        buffer_ids.insert(buf->getIndex());
        num_bytes_assoc_with_views += view.getTotalBytes();
      }
    }
    else
    {
      num_views_empty += 1;
    }
  }

  //
  // Update Node entries with data info for this Group
  //
  n["num_groups"] = num_groups;
  n["num_views"] = num_views;
  n["num_views_empty"] = num_views_empty;
  n["num_views_buffer"] = num_views_buffer;
  n["num_views_external"] = num_views_external;
  n["num_views_scalar"] = num_views_scalar;
  n["num_views_string"] = num_views_string;
  n["num_bytes_assoc_with_views"] = num_bytes_assoc_with_views;
  n["num_bytes_external"] = num_bytes_external;

  //
  // Recursively gather info for Group subtree, if requested
  //
  if(recursive)
  {
    for(auto& group : groups())
    {
      group.getDataInfoHelper(n, buffer_ids, recursive);
    }
  }
}

////////////////////////////////////////////////////////////////////////
//
// View query methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return true if Group owns a View with given name or path; else false.
 *
 *************************************************************************
 */
bool Group::hasView(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
  {
    return false;
  }

  return group->hasChildView(intpath);
}

////////////////////////////////////////////////////////////////////////
//
// View access methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return pointer to non-const View with given name or path if it exists.
 *
 *************************************************************************
 */
View* Group::getView(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path " << path);
    return nullptr;
  }

  SLIC_CHECK_MSG(!intpath.empty() && group->hasChildView(intpath),
                 SIDRE_GROUP_LOG_PREPEND << "No View with name '" << intpath << "'");

  return group->getNamedViews()->getItem(intpath);
}

/*
 *************************************************************************
 *
 * Return pointer to const View with given name or path if it exists.
 *
 *************************************************************************
 */
const View* Group::getView(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path " << path);
    return nullptr;
  }

  SLIC_CHECK_MSG(!intpath.empty() && group->hasChildView(intpath),
                 SIDRE_GROUP_LOG_PREPEND << "No View with name '" << intpath << "'");

  return group->getNamedViews()->getItem(intpath);
}

////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View that has no associated data.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create empty View (i.e., no data description) with given name or path
 * in this Group.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path)
{
  std::string intpath(path);

  Group* group;
  if(intpath.empty())
  {
    SLIC_CHECK_MSG(m_is_list,
                   SIDRE_GROUP_LOG_PREPEND << "Could not create View with empty string "
                                           << "for the path.");
    if(m_is_list)
    {
      group = this;
    }
    else
    {
      return nullptr;
    }
  }
  else
  {
    bool create_groups_in_path = true;
    group = walkPath(intpath, create_groups_in_path);

    if(group == nullptr)
    {
      if(m_is_list)
      {
        SLIC_CHECK_MSG(group != nullptr,
                       SIDRE_GROUP_LOG_PREPEND << "Could not find or create View '" << path << "'."
                                               << "for a Group using the list format.");
      }
      else
      {
        SLIC_CHECK_MSG(group != nullptr,
                       SIDRE_GROUP_LOG_PREPEND << "Could not find or create path '" << path << "'."
                                               << "There is already a View with that name.");
      }
      return nullptr;
    }
    else if(intpath.empty() || group->hasChildView(intpath) || group->hasChildGroup(intpath))
    {
      SLIC_CHECK_MSG(!intpath.empty(),
                     SIDRE_GROUP_LOG_PREPEND << "Cannot create a View with an empty path.");
      SLIC_CHECK_MSG(!group->hasChildView(intpath),
                     SIDRE_GROUP_LOG_PREPEND << "Cannot create View with name '" << intpath << "'. "
                                             << "There is already a View with that name.");
      SLIC_CHECK_MSG(!group->hasChildGroup(intpath),
                     SIDRE_GROUP_LOG_PREPEND << "Cannot create View with name '" << intpath << "'. "
                                             << "There is already has a Group with that name.");
      return nullptr;
    }
  }

  View* view = new(std::nothrow) View(intpath);
  if(view != nullptr)
  {
    group->attachView(view);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, TypeID type, IndexType num_elems)
{
  if(type == NO_TYPE_ID || num_elems < 0)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   SIDRE_GROUP_LOG_PREPEND << "Cannot create View with name '" << path << "'. "
                                           << " Invalid type << " << type);
    SLIC_CHECK_MSG(num_elems >= 0,
                   SIDRE_GROUP_LOG_PREPEND << "Cannot create View with name '" << path << "'. "
                                           << "Number of elements cannot be less than zero.");
    return nullptr;
  }

  View* view = createView(path);
  if(view != nullptr)
  {
    view->describe(type, num_elems);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in this
 * Group.
 *
 *************************************************************************
 */
View* Group::createViewWithShape(const std::string& path, TypeID type, int ndims, const IndexType* shape)
{
  if(type == NO_TYPE_ID || ndims < 0 || shape == nullptr)
  {
    SLIC_CHECK_MSG(type != NO_TYPE_ID,
                   SIDRE_GROUP_LOG_PREPEND << "Problem creating View with name '" << path << "'. "
                                           << " Invalid type: " << type);
    SLIC_CHECK_MSG(ndims >= 0,
                   SIDRE_GROUP_LOG_PREPEND << "Problem creating View with name '" << path << "'. "
                                           << "ndims must be greater than 0.");
    SLIC_CHECK_MSG(shape != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Problem creating View with name '" << path << "'. "
                                           << "shape pointer was null.");
    return nullptr;
  }

  View* view = createView(path);
  if(view != nullptr)
  {
    view->describe(type, ndims, shape);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create View with given name and data described by a conduit Datatype object.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, const DataType& dtype)
{
  View* view = createView(path);
  if(view != nullptr)
  {
    view->describe(dtype);
  }

  return view;
}

////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View and attach a Buffer to it.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create an undescribed View with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, Buffer* buff)
{
  View* view = createView(path);
  if(view != nullptr)
  {
    view->attachBuffer(buff);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, TypeID type, IndexType num_elems, Buffer* buff)
{
  View* view = createView(path, type, num_elems);
  if(view != nullptr)
  {
    view->attachBuffer(buff);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
View* Group::createViewWithShape(const std::string& path,
                                 TypeID type,
                                 int ndims,
                                 const IndexType* shape,
                                 Buffer* buff)
{
  View* view = createViewWithShape(path, type, ndims, shape);
  if(view != nullptr)
  {
    view->attachBuffer(buff);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (Conduit DataType) with given name or path in
 * this Group and attach Buffer to it.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, const DataType& dtype, Buffer* buff)
{
  View* view = createView(path, dtype);
  if(view != nullptr)
  {
    view->attachBuffer(buff);
  }

  return view;
}

////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View and attach external data to it.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create an undescribed View with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, void* external_ptr)
{
  View* view = createView(path);
  if(view != nullptr)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, TypeID type, IndexType num_elems, void* external_ptr)
{
  View* view = createView(path, type, num_elems);
  if(view != nullptr)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
View* Group::createViewWithShape(const std::string& path,
                                 TypeID type,
                                 int ndims,
                                 const IndexType* shape,
                                 void* external_ptr)
{
  View* view = createViewWithShape(path, type, ndims, shape);
  if(view != nullptr)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (Conduit DataType) with given name or path in
 * this Group and attach external data ptr to it.
 *
 *************************************************************************
 */
View* Group::createView(const std::string& path, const DataType& dtype, void* external_ptr)
{
  View* view = createView(path, dtype);
  if(view != nullptr)
  {
    view->setExternalDataPtr(external_ptr);
  }
  return view;
}

////////////////////////////////////////////////////////////////////////
//
//  Methods to create a View and allocate its data.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create described View (type and # elems) with given name or path in
 * this Group and allocate its data.
 *
 *************************************************************************
 */
View* Group::createViewAndAllocate(const std::string& path, TypeID type, IndexType num_elems, int allocID)
{
  View* view = createView(path, type, num_elems);
  if(view != nullptr)
  {
    allocID = getValidArrayAllocatorId(allocID);
    view->allocate(allocID);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (type and shape) with given name or path in
 * this Group and allocate its data.
 *
 *************************************************************************
 */
View* Group::createViewWithShapeAndAllocate(const std::string& path,
                                            TypeID type,
                                            int ndims,
                                            const IndexType* shape,
                                            int allocID)
{
  View* view = createViewWithShape(path, type, ndims, shape);
  if(view != nullptr)
  {
    allocID = getValidArrayAllocatorId(allocID);
    view->allocate(allocID);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create described View (Conduit DataType) with given name or path in
 * this Group and allocate its data.
 *
 *************************************************************************
 */
View* Group::createViewAndAllocate(const std::string& path, const DataType& dtype, int allocID)
{
  View* view = createView(path, dtype);
  if(view != nullptr)
  {
    allocID = getValidArrayAllocatorId(allocID);
    view->allocate(allocID);
  }
  return view;
}

/*
 *************************************************************************
 *
 * Create View with given name or path and set its data to given string.
 *
 *************************************************************************
 */
View* Group::createViewString(const std::string& path, const std::string& value, int allocID)
{
  View* view = createView(path);
  if(view != nullptr)
  {
    view->setString(value, allocID);
  }

  return view;
}

////////////////////////////////////////////////////////////////////////
//
//  Methods for destroying Views and their data.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Destroy View with given name or path and leave its data intact.
 *
 *************************************************************************
 */
void Group::destroyView(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group != nullptr)
  {
    View* view = group->detachView(intpath);
    if(view != nullptr)
    {
      delete view;
    }
  }
}

/*
 *************************************************************************
 *
 * Destroy View with given index and leave its data intact.
 *
 *************************************************************************
 */
void Group::destroyView(IndexType idx)
{
  View* view = detachView(idx);
  if(view != nullptr)
  {
    delete view;
  }
}

/*
 *************************************************************************
 *
 * Destroy all Views in Group and leave their data intact.
 *
 *************************************************************************
 */
void Group::destroyViews()
{
  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    View* view = detachView(vidx);
    if(view != nullptr)
    {
      delete view;
    }

    vidx = getNextValidViewIndex(vidx);
  }

  m_view_coll->removeAllItems();
}

/*
 *************************************************************************
 *
 * Destroy View with given name or path and its data if it's the only View
 * associated with that data.
 *
 *************************************************************************
 */
void Group::destroyViewAndData(const std::string& path) { destroyViewAndData(getView(path)); }

/*
 *************************************************************************
 *
 * Destroy View with given index and its data if it's the only View
 * associated with that data.
 *
 *************************************************************************
 */
void Group::destroyViewAndData(IndexType idx) { destroyViewAndData(getView(idx)); }

/*
 *************************************************************************
 *
 * Destroy all Views in Group as well as the data for each View when it's
 * the only View associated with that data.
 *
 *************************************************************************
 */
void Group::destroyViewsAndData()
{
  IndexType vidx = getFirstValidViewIndex();
  while(indexIsValid(vidx))
  {
    destroyViewAndData(vidx);
    vidx = getNextValidViewIndex(vidx);
  }

  m_view_coll->removeAllItems();
}

////////////////////////////////////////////////////////////////////////
//
//  Methods for moving and copying View objects from one Group to another.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Remove given View from its owning Group and attach to this Group.
 *
 *************************************************************************
 */
View* Group::moveView(View* view)
{
  if(view == nullptr)
  {
    SLIC_CHECK_MSG(view != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no View to move.");
    return nullptr;
  }

  Group* curr_group = view->getOwningGroup();
  if(curr_group == this)
  {
    // this Group already owns the View
    return view;
  }
  else if(hasChildView(view->getName()))
  {
    SLIC_CHECK_MSG(!hasChildView(view->getName()),
                   SIDRE_GROUP_LOG_PREPEND << "Group already has a View named '" << view->getName()
                                           << "' so View move operation cannot happen");
    return nullptr;
  }

  curr_group->detachView(view);
  attachView(view);

  return view;
}

/*
 *************************************************************************
 *
 * Create a copy of given View and attach to this Group.
 *
 * Copying a View does not perform a deep copy of its data Buffer.
 *
 *************************************************************************
 */
View* Group::copyView(View* view)
{
  if(view == nullptr || hasChildView(view->getName()))
  {
    SLIC_CHECK_MSG(view != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no View to copy.");

    if(view != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildView(view->getName()),
                     SIDRE_GROUP_LOG_PREPEND << "Group already has a View named '" << view->getName()
                                             << "' so View copy operation cannot happen");
    }

    return nullptr;
  }

  View* copy = createView(view->getName());
  view->copyView(copy);
  return copy;
}

/*
 *************************************************************************
 *
 * Create a deep copy of given View and attach to this Group.
 *
 * The deep copy performs a copy of all described data.
 *
 *************************************************************************
 */
View* Group::deepCopyView(const View* view, int arrayAllocId, int tupleAllocId)
{
  if(view == nullptr || hasChildView(view->getName()))
  {
    SLIC_CHECK_MSG(view != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no View to copy.");

    if(view != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildView(view->getName()),
                     SIDRE_GROUP_LOG_PREPEND << "Group already has a View named '" << view->getName()
                                             << "' so View copy operation cannot happen");
    }

    return nullptr;
  }

  View* copy = createView(view->getName());
  arrayAllocId = getValidArrayAllocatorId(arrayAllocId);
  tupleAllocId = getValidTupleAllocatorId(tupleAllocId);
  view->deepCopyView(copy, arrayAllocId, tupleAllocId);
  return copy;
}

////////////////////////////////////////////////////////////////////////
//
// Child Group query methods.
//
////////////////////////////////////////////////////////////////////////

/*
 ***********************************************************************
 *
 * Return true if this Group has a descendant Group with given name or path;
 * else false.
 *
 ***********************************************************************
 */
bool Group::hasGroup(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
  {
    return false;
  }
  else
  {
    return group->hasChildGroup(intpath);
  }
}

////////////////////////////////////////////////////////////////////////
//
// Child Group access methods.
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Return pointer to non-const child Group with given name or path
 * if it exists.
 *
 *************************************************************************
 */
Group* Group::getGroup(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path '" << path << "'.");
    return nullptr;
  }

  SLIC_CHECK_MSG(!intpath.empty() && group->hasChildGroup(intpath),
                 SIDRE_GROUP_LOG_PREPEND << "Group has no descendant Group named '" << path << "'.");

  return group->getNamedGroups()->getItem(intpath);
}

/*
 *************************************************************************
 *
 * Return pointer to const child Group with given name or path if it exists.
 *
 *************************************************************************
 */
const Group* Group::getGroup(const std::string& path) const
{
  std::string intpath(path);
  const Group* group = walkPath(intpath);

  if(group == nullptr)
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Non-existent group in path " << path);
    return nullptr;
  }

  SLIC_CHECK_MSG(
    !intpath.empty() && group->hasChildGroup(intpath),
    SIDRE_GROUP_LOG_PREPEND << "Group has no descendant Group with name '" << path << "'.");

  return group->getNamedGroups()->getItem(intpath);
}

////////////////////////////////////////////////////////////////////////
//
//  Methods for managing child Group objects in Group
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * Create Group with given name or path and make it a child of this Group.
 *
 *************************************************************************
 */
Group* Group::createGroup(const std::string& path, bool is_list, bool accept_existing)
{
  std::string intpath(path);
  bool create_groups_in_path = true;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group == nullptr)
  {
    if(m_is_list)
    {
      SLIC_CHECK_MSG(group != nullptr,
                     SIDRE_GROUP_LOG_PREPEND << "Could not find or create Group '" << path << "'."
                                             << "for a Group using the list format.");
    }
    else
    {
      SLIC_CHECK_MSG(group != nullptr,
                     SIDRE_GROUP_LOG_PREPEND << "Could not find or create path '" << path
                                             << "'. There is already a Group with that name");
    }
    return nullptr;
  }
  else if(intpath.empty() || (group->hasChildGroup(intpath) && !accept_existing) ||
          group->hasChildView(intpath))
  {
    SLIC_CHECK_MSG(!intpath.empty(),
                   SIDRE_GROUP_LOG_PREPEND << "Cannot create a Group with an empty path.");
    SLIC_CHECK_MSG(!group->hasChildGroup(intpath),
                   SIDRE_GROUP_LOG_PREPEND << "Cannot create Group with name '" << path
                                           << "'. There is already has a Group with that name.");
    SLIC_CHECK_MSG(!group->hasChildView(intpath),
                   SIDRE_GROUP_LOG_PREPEND << "Cannot create Group with name '" << path
                                           << "'. There is already has a View with that name.");

    return nullptr;
  }

  if(!group->hasGroup(intpath))
  {
    Group* new_group = new(std::nothrow) Group(intpath, group->getDataStore(), is_list);
    if(new_group == nullptr)
    {
      return nullptr;
    }

    SLIC_ASSERT(group->getDefaultArrayAllocatorID() == m_default_array_alloc_id);
    SLIC_ASSERT(group->getDefaultTupleAllocatorID() == m_default_tuple_alloc_id);
    new_group->setDefaultArrayAllocator(m_default_array_alloc_id);
    new_group->setDefaultTupleAllocator(m_default_tuple_alloc_id);
    return group->attachGroup(new_group);
  }

  return group->getGroup(intpath);
}

Group* Group::createUnnamedGroup(bool is_list)
{
  SLIC_CHECK_MSG(m_is_list,
                 SIDRE_GROUP_LOG_PREPEND << "Cannot create an unnamed Group when not using "
                                         << "list format.");

  Group* new_group;
  if(m_is_list)
  {
    new_group = new(std::nothrow) Group("", getDataStore(), is_list);
  }
  else
  {
    new_group = nullptr;
  }

  if(new_group == nullptr)
  {
    return nullptr;
  }

#ifdef AXOM_USE_UMPIRE
  // Why only do this with Umpire?  BTNG.
  new_group->setDefaultArrayAllocator(getDefaultArrayAllocator());
  new_group->setDefaultTupleAllocator(getDefaultTupleAllocator());
#endif
  return attachGroup(new_group);
}

/*
 *************************************************************************
 *
 * Detach child Group with given name or path and destroy it.
 *
 *************************************************************************
 */
void Group::destroyGroup(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group != nullptr)
  {
    Group* targetgroup = group->detachGroup(intpath);
    if(targetgroup != nullptr)
    {
      delete targetgroup;
    }
  }
}

/*
 *************************************************************************
 *
 * Detach child Group with given name or path and destroy it, and destroy
 * any data that becomes disassociated from all Views
 *
 *************************************************************************
 */
void Group::destroyGroupAndData(const std::string& path)
{
  std::string intpath(path);
  bool create_groups_in_path = false;
  Group* group = walkPath(intpath, create_groups_in_path);

  if(group != nullptr)
  {
    Group* targetgroup = group->detachGroup(intpath);
    if(targetgroup != nullptr)
    {
      targetgroup->destroyGroupSubtreeAndData();
      delete targetgroup;
    }
  }
}

/*
 *************************************************************************
 *
 * Detach child Group with given index and destroy it.
 *
 *************************************************************************
 */
void Group::destroyGroup(IndexType idx)
{
  Group* group = detachGroup(idx);
  if(group != nullptr)
  {
    delete group;
  }
}

/*
 *************************************************************************
 *
 * Detach child Group with given index and destroy it, and destroy any data
 * that becomes disassociated from all Views
 *
 *************************************************************************
 */
void Group::destroyGroupAndData(IndexType idx)
{
  Group* group = detachGroup(idx);
  if(group != nullptr)
  {
    group->destroyGroupSubtreeAndData();
    delete group;
  }
}

/*
 *************************************************************************
 *
 * Detach all child Groups and destroy them.
 *
 *************************************************************************
 */
void Group::destroyGroups()
{
  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    Group* group = this->getGroup(gidx);
    delete group;

    gidx = getNextValidGroupIndex(gidx);
  }

  m_group_coll->removeAllItems();
}

/*
 *************************************************************************
 *
 * Detach all child Groups and destroy them, and destroy any data that
 * becomes disassociated with all Views.
 *
 *************************************************************************
 */
void Group::destroyGroupsAndData()
{
  IndexType gidx = getFirstValidGroupIndex();
  while(indexIsValid(gidx))
  {
    destroyGroupAndData(gidx);

    gidx = getNextValidGroupIndex(gidx);
  }

  m_group_coll->removeAllItems();
}

/*
 *************************************************************************
 *
 * Detach all child Groups and Views and destroy them, and destroy any data
 * that becomes disassociated with all Views.
 *
 *************************************************************************
 */
void Group::destroyGroupSubtreeAndData()
{
  destroyViewsAndData();
  destroyGroupsAndData();
}

/*
 *************************************************************************
 *
 * Remove given Group from its owning Group and make it a child of this Group.
 *
 *************************************************************************
 */
Group* Group::moveGroup(Group* group)
{
  if(group == nullptr || hasChildGroup(group->getName()))
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no Group to move.");

    if(group != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildGroup(group->getName()),
                     SIDRE_GROUP_LOG_PREPEND << "Invalid move operation. Group already has "
                                             << "a child named '" << group->getName() << "'.");
    }

    return nullptr;
  }

  Group* curr_group = group->getParent();
  curr_group->detachGroup(group->getName());
  attachGroup(group);
  return group;
}

/*
 *************************************************************************
 *
 * Create a copy of given Group and make it a child of this Group.
 *
 * Copying a Group does not perform a deep copy of any of its Buffers.
 *
 *************************************************************************
 */
Group* Group::copyGroup(Group* group)
{
  if(group == nullptr || hasChildGroup(group->getName()))
  {
    SLIC_CHECK_MSG(group != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no Group to copy.");

    if(group != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildGroup(group->getName()),
                     SIDRE_GROUP_LOG_PREPEND << "Invalid copy operation. Group already has "
                                             << "a child named '" << group->getName() << "'.");
    }

    return nullptr;
  }

  Group* res = createGroup(group->getName());

  // copy child Groups to new Group
  for(auto& grp : group->groups())
  {
    res->copyGroup(&grp);
  }

  // copy Views to new Group
  for(auto& view : group->views())
  {
    res->copyView(&view);
  }

  return res;
}

/*
 *************************************************************************
 *
 * Create a deep copy of given Group and make the copy a child of this Group.
 *
 * The deep copy of a Group will copy the group hierarchy and deep copy
 * all Views within the hierarchy.
 *
 *************************************************************************
 */
Group* Group::deepCopyGroup(const Group* srcGroup, int arrayAllocId, int tupleAllocId)
{
  SLIC_ASSERT(srcGroup != this);
  if(srcGroup == nullptr || hasChildGroup(srcGroup->getName()))
  {
    SLIC_CHECK_MSG(srcGroup != nullptr,
                   SIDRE_GROUP_LOG_PREPEND << "Null pointer provided, no Group to copy.");

    if(srcGroup != nullptr)
    {
      SLIC_CHECK_MSG(!hasChildGroup(srcGroup->getName()),
                     SIDRE_GROUP_LOG_PREPEND << "Invalid copy operation. Group already has "
                                             << "a child named '" << srcGroup->getName() << "'.");
    }

    return nullptr;
  }

  Group* dstGroup = createGroup(srcGroup->getName());

  // copy child Groups to new Group
  for(auto& grp : srcGroup->groups())
  {
    dstGroup->deepCopyGroup(&grp, arrayAllocId, tupleAllocId);
  }

  // copy Views to new Group
  for(auto& view : srcGroup->views())
  {
    dstGroup->deepCopyView(&view, arrayAllocId, tupleAllocId);
  }

  return dstGroup;
}

/*
 *************************************************************************
 *
 * Deep copy another Group into this Group.
 *
 * All contents in this group and their data are destroyed before the copy.
 * The deep copy of a Group will copy the group hierarchy and deep copy
 * all Views within the hierarchy.
 *
 *************************************************************************
 */
Group* Group::deepCopyGroupToSelf(const Group* srcGroup)
{
  SLIC_ASSERT(srcGroup != this);
  SLIC_ERROR_IF(m_is_list && !srcGroup->m_is_list,
                "Group::deepCopyToSelf cannot copy from a list Group '" + srcGroup->getPath() +
                  "' to a non-list Group '" + getPath() + "'");

  destroyGroupsAndData();
  destroyViewsAndData();

  // copy child Groups to new Group
  for(auto& grp : srcGroup->groups())
  {
    deepCopyGroup(&grp);
  }

  // copy Views to new Group
  for(auto& view : srcGroup->views())
  {
    deepCopyView(&view);
  }

  return this;
}

/*
 *************************************************************************
 *
 * Copy Group native layout to given Conduit node.
 *
 *************************************************************************
 */
bool Group::createNativeLayout(Node& n, const Attribute* attr) const
{
  n.set(DataType::object());
  bool hasSavedViews = false;

  // Dump the group's views
  for(const auto& view : this->views())
  {
    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG(
      m_is_list || !hasChildGroup(view.getName()),
      SIDRE_GROUP_LOG_PREPEND << axom::fmt::format("'{}' is the name of both a group and a view",
                                                   view.getName()));

    if(attr == nullptr || view.hasAttributeValue(attr))
    {
      conduit::Node& child_node = m_is_list ? n.append() : n[view.getName()];
      view.createNativeLayout(child_node);
      hasSavedViews = true;
    }
  }

  // Recursively dump the child groups
  for(const auto& group : this->groups())
  {
    conduit::Node& child_node = m_is_list ? n.append() : n[group.getName()];
    if(group.createNativeLayout(child_node, attr))
    {
      hasSavedViews = true;
    }
    else
    {
      if(m_is_list)
      {
        n.remove(group.getName());
      }
      else
      {
        n.remove(n.number_of_children() - 1);
      }
    }
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * Deep-copy Group's native layout to given Conduit node.
 *
 *************************************************************************
 */
bool Group::deepCopyToConduit(Node& dst, int tupleAllocId, int arrayAllocId, const Attribute* attr) const
{
  if(tupleAllocId == INVALID_ALLOCATOR_ID)
  {
    tupleAllocId = ConduitMemory::conduitAllocIdToAxom(dst.allocator());
  }
  if(arrayAllocId == INVALID_ALLOCATOR_ID)
  {
    arrayAllocId = ConduitMemory::conduitAllocIdToAxom(dst.allocator());
  }

  dst.set(DataType::object());
  bool hasSavedViews = false;

  // Dump the group's views
  for(auto& view : views())
  {
    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG(m_is_list || !hasChildGroup(view.getName()),
                   SIDRE_GROUP_LOG_PREPEND << "'" << view.getName()
                                           << "' is the name of both a group and a view.");

    if(attr == nullptr || view.hasAttributeValue(attr))
    {
      conduit::Node& child_node = m_is_list ? dst.append() : dst[view.getName()];
      const int allocId = view.isScalar() || view.isString() ? tupleAllocId : arrayAllocId;
      if(allocId != axom::INVALID_ALLOCATOR_ID)
      {
        const conduit::index_t conduitAllocId = ConduitMemory::axomAllocIdToConduit(allocId);
        child_node.set_allocator(conduitAllocId);
      }
      view.deepCopyToConduit(child_node);
      hasSavedViews = true;
    }
  }

  // Recursively dump the child groups
  for(auto& group : groups())
  {
    conduit::Node& child_node = m_is_list ? dst.append() : dst[group.getName()];
    if(group.deepCopyToConduit(child_node, tupleAllocId, arrayAllocId, attr))
    {
      hasSavedViews = true;
    }
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * Copy Group layout to the given Conduit node, including only the metadata
 * for Views and not the actual data the Views hold.
 *
 *************************************************************************
 */

void Group::createNoDataLayout(Node& n, const Attribute* attr) const
{
  n.set(DataType::object());

  // Dump the group's views
  for(auto& view : views())
  {
    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG(
      m_is_list || !hasChildGroup(view.getName()),
      SIDRE_GROUP_LOG_PREPEND << axom::fmt::format("'{}' is the name of both a group and a view",
                                                   view.getName()));

    if(attr == nullptr || view.hasAttributeValue(attr))
    {
      conduit::Node& child_node = m_is_list ? n.append() : n[view.getName()];
      view.copyMetadataToNode(child_node);
    }
  }

  // Recursively dump the child groups
  for(auto& group : groups())
  {
    conduit::Node& child_node = m_is_list ? n.append() : n[group.getName()];
    group.createNoDataLayout(child_node, attr);
  }
}

/*
 *************************************************************************
 *
 * Adds a conduit node for this group if it has external views,
 * or if any of its children groups has an external view
 *
 *************************************************************************
 * see ATK-736 - Improvements to createNativeLayout and createExternalLayout
 */
bool Group::createExternalLayout(Node& n, const Attribute* attr) const
{
  n.set(DataType::object());

  bool hasExternalViews = false;

  // Dump the group's views
  for(const auto& view : this->views())
  {
    // Check that the view's name is not also a child group name
    SLIC_CHECK_MSG(
      m_is_list || !hasChildGroup(view.getName()),
      SIDRE_GROUP_LOG_PREPEND << axom::fmt::format("'{}' is the name of both a group and a view.",
                                                   view.getName()));

    if(attr == nullptr || view.hasAttributeValue(attr))
    {
      if(view.isExternal())
      {
        if(view.isDescribed())
        {
          conduit::Node& vnode = m_is_list ? n.append() : n[view.getName()];
          view.createNativeLayout(vnode);
        }
        hasExternalViews = true;
      }
    }
  }

  for(const auto& group : this->groups())
  {
    conduit::Node& gnode = m_is_list ? n.append() : n[group.getName()];
    if(group.createExternalLayout(gnode, attr))
    {
      hasExternalViews = true;
    }
    else
    {
      // Remove nodes that do not have any external views
      if(m_is_list)
      {
        n.remove(n.number_of_children() - 1);
      }
      else
      {
        n.remove(group.getName());
      }
    }
  }

  return hasExternalViews;
}

/*
 *************************************************************************
 *
 * Print JSON description of data Group to stdout.
 *
 *************************************************************************
 */
void Group::print() const { print(std::cout); }

/*
 *************************************************************************
 *
 * Print JSON description of data Group to an ostream.
 *
 *************************************************************************
 */
void Group::print(std::ostream& os) const
{
  Node n;
  copyToConduitNode(n);
  n.to_json_stream(os);
}

/*
 *************************************************************************
 *
 * Print sub-tree rooted at this Group to stdout, with level-based padding.
 *
 *************************************************************************
 */
void Group::printTree(const int nlevels, std::ostream& os) const
{
  // lambda to generate the proper number of padding spaces
  auto pad = [](int n) { return axom::fmt::format("{:>{}}", "", 4 * n); };

  os << axom::fmt::format("{}Group {}\n", pad(nlevels), this->getName());

  for(const auto& view : this->views())
  {
    os << axom::fmt::format("{}View {}\n", pad(nlevels + 1), view.getName());
  }

  for(const auto& group : this->groups())
  {
    group.printTree(nlevels + 1, os);
  }
}

void Group::hostPrint(const std::string& indent, std::ostream& os) const
{
  for(const auto& view : this->views())
  {
    os << indent << view.getName() << ':';
    view.hostPrint(os);
    os << std::endl;
  }

  const std::string nextIndent = indent + std::string("  ");
  for(const auto& group : this->groups())
  {
    os << indent << group.getName() << ':' << std::endl;
    group.hostPrint(nextIndent, os);
  }
}

/*
 *************************************************************************
 *
 * Copy data Group description to given Conduit node.
 *
 *************************************************************************
 */
void Group::copyToConduitNode(Node& n) const
{
  n["name"] = m_name;

  for(const auto& view : this->views())
  {
    Node& v = isUsingMap() ? n["views"].fetch(view.getName()) : n["views"].append();
    view.copyToConduitNode(v);
  }

  for(const auto& group : this->groups())
  {
    Node& g = isUsingMap() ? n["groups"].fetch(group.getName()) : n["groups"].append();
    group.copyToConduitNode(g);
  }
}

/*
 *************************************************************************
 *
 * Test this Group for equivalence to another Group.
 *
 *************************************************************************
 */
bool Group::isEquivalentTo(const Group* other, bool checkName) const
{
  // Equality of names
  bool is_equiv = checkName ? m_name == other->m_name : true;

  // Sizes of collections of child items must be equal
  is_equiv = is_equiv && (m_view_coll->getNumItems() == other->m_view_coll->getNumItems()) &&
    (m_group_coll->getNumItems() == other->m_group_coll->getNumItems());

  // Test equivalence of Views
  if(is_equiv)
  {
    for(const auto& view : this->views())
    {
      const std::string& name = view.getName();

      is_equiv = is_equiv && other->hasChildView(name) && view.isEquivalentTo(other->getView(name));
    }
  }

  // Recursively call this method to test equivalence of child Groups
  if(is_equiv)
  {
    for(const auto& group : this->groups())
    {
      const std::string& name = group.getName();

      is_equiv =
        is_equiv && other->hasChildGroup(name) && group.isEquivalentTo(other->getGroup(name));
    }
  }

  return is_equiv;
}

/*
 *************************************************************************
 *
 * Private class used to set Conduit error handler to default (throw exception)
 *
 * Creating an instance of this class will disable the current Conduit error
 * callbacks.  The default Conduit callbacks throw exceptions when they
 * encounter I/O errors.  When the instance is destroyed, the previous
 * callbacks are restored.
 *
 * This class is also a functor.  The operator() method takes a closure,
 * catches any Conduit error, and saves that in the DataStore.
 *
 *************************************************************************
 */

class ConduitErrorSuppressor
{
public:
  using conduit_error_handler = void (*)(const std::string&, const std::string&, int);

  ConduitErrorSuppressor(const DataStore* ds, bool suppress_in_call = true);
  ~ConduitErrorSuppressor();

  void operator()(const std::function<void(void)>& conduitOp)
  {
    if(m_suppress_in_call)
    {
      disable_conduit_error_handlers();
    }

    try
    {
      conduitOp();
    }
    catch(conduit::Error& e)
    {
      m_ds->appendToConduitErrors(e.what());
    }

    if(m_suppress_in_call)
    {
      restore_conduit_error_handlers();
    }
  }

private:
  // saves current error func; set to default
  void disable_conduit_error_handlers()
  {
    m_info_handler = conduit::utils::info_handler();
    m_warning_handler = conduit::utils::warning_handler();
    m_error_handler = conduit::utils::error_handler();

    DataStore::setConduitDefaultMessageHandlers();
  }

  // restores saved error func
  void restore_conduit_error_handlers()
  {
    conduit::utils::set_error_handler(m_error_handler);
    conduit::utils::set_warning_handler(m_warning_handler);
    conduit::utils::set_info_handler(m_info_handler);
  }

  /// The DataStore we report to
  const DataStore* m_ds;

  /// callbacks used for Conduit error interface
  conduit_error_handler m_error_handler, m_warning_handler, m_info_handler;

  /// Suppress error handling in the call operator?
  bool m_suppress_in_call;
};

ConduitErrorSuppressor::ConduitErrorSuppressor(const DataStore* ds, bool suppress_in_call)
  : m_ds(ds)
  , m_error_handler(nullptr)
  , m_warning_handler(nullptr)
  , m_info_handler(nullptr)
  , m_suppress_in_call(suppress_in_call)
{
  if(!m_suppress_in_call)
  {
    disable_conduit_error_handlers();
  }
}

ConduitErrorSuppressor::~ConduitErrorSuppressor()
{
  if(!m_suppress_in_call)
  {
    restore_conduit_error_handlers();
  }
}

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a file
 *
 *************************************************************************
 */

bool Group::save(const std::string& path, const std::string& protocol, const Attribute* attr) const
{
  const DataStore* ds = getDataStore();
  ConduitErrorSuppressor checkConduitCall(ds);
  bool retval = false;

  if(protocol == "sidre_hdf5")
  {
    Node n;
    exportTo(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, "hdf5"); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "sidre_conduit_json")
  {
    Node n;
    exportTo(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, "conduit_json"); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "sidre_json")
  {
    Node n;
    exportTo(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, "json"); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "sidre_layout_json")
  {
    Node n;
    exportWithoutBufferData(n["sidre"], attr);
    ds->saveAttributeLayout(n["sidre/attribute"]);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, "conduit_json"); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "conduit_hdf5")
  {
    Node n;
    createNativeLayout(n, attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, "hdf5"); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "conduit_bin" || protocol == "conduit_json" || protocol == "json")
  {
    Node n;
    createNativeLayout(n, attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, protocol); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "conduit_layout_json")
  {
    Node n;
    createNoDataLayout(n, attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::save(n, path, "conduit_json"); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol << "' for file save.");
    retval = false;
  }

  return retval;
}

/*************************************************************************/

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from a file
 *
 *************************************************************************
 */
bool Group::load(const std::string& path, const std::string& protocol, bool preserve_contents)
{
  std::string new_name;
  return load(path, protocol, preserve_contents, new_name);
}

bool Group::load(const std::string& path,
                 const std::string& protocol,
                 bool preserve_contents,
                 std::string& name_from_file)
{
  ConduitErrorSuppressor checkConduitCall(getDataStore());
  bool retval = false;
  if(protocol == "sidre_hdf5")
  {
    Node n;
    checkConduitCall([&] { conduit::relay::io::load(path, "hdf5", n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      SLIC_ASSERT_MSG(n.has_path("sidre"),
                      SIDRE_GROUP_LOG_PREPEND << "Conduit Node " << n.path() << " does not have sidre "
                                              << "data for this Group " << getPathName() << ".");
      importFrom(n["sidre"], preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else if(protocol == "sidre_conduit_json")
  {
    Node n;
    checkConduitCall([&] { conduit::relay::io::load(path, "conduit_json", n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      SLIC_ASSERT_MSG(n.has_path("sidre"),
                      SIDRE_GROUP_LOG_PREPEND << "Conduit Node " << n.path() << " does not have sidre "
                                              << "data for Group " << getPathName() << ".");
      importFrom(n["sidre"], preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else if(protocol == "sidre_json")
  {
    Node n;
    checkConduitCall([&] { conduit::relay::io::load(path, "json", n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      SLIC_ASSERT_MSG(n.has_path("sidre"),
                      SIDRE_GROUP_LOG_PREPEND << "Conduit Node " << n.path() << " does not have sidre "
                                              << "data for Group " << getPathName() << ".");
      importFrom(n["sidre"], preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else if(protocol == "conduit_hdf5")
  {
    Node n;
    checkConduitCall([&] { conduit::relay::io::load(path, "hdf5", n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      importConduitTree(n, preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else if(protocol == "conduit_bin" || protocol == "conduit_json" || protocol == "json")
  {
    Node n;
    checkConduitCall([&] { conduit::relay::io::load(path, protocol, n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      importConduitTree(n, preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol << "' for file load.");
    retval = false;
  }

  return retval;
}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from a file
 *
 *************************************************************************
 */
Group* Group::createGroupAndLoad(std::string& group_name,
                                 const std::string& path,
                                 const std::string& protocol,
                                 bool& load_success)
{
  load_success = false;
  Group* child = createGroup(group_name);
  if(child != nullptr)
  {
    load_success = child->load(path, protocol, false, group_name);
    if(!load_success)
    {
      destroyGroupAndData(group_name);
      child = nullptr;
    }
  }

  return child;
}

/*
 *************************************************************************
 *
 * Load External Data from a file
 *
 *************************************************************************
 */
bool Group::loadExternalData(const std::string& path)
{
  Node n;
  createExternalLayout(n);
  ConduitErrorSuppressor checkConduitCall(getDataStore());
  bool success;

#ifdef AXOM_USE_HDF5
  // CYRUS'-NOTE, not sure ":" will work with multiple trees per
  // output file
  checkConduitCall([&] { conduit::relay::io::hdf5_read(path + ":sidre/external", n); });

  success = !(getDataStore()->getConduitErrorOccurred());
#else
  AXOM_UNUSED_VAR(path);
  SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND << "External data not loaded. "
                                       << "This function requires hdf5 support. "
                                       << " Please reconfigure with hdf5.");

  success = false;
#endif

  return success;
}

// Functions that directly use the hdf5 API in their signature
#ifdef AXOM_USE_HDF5

/*
 *************************************************************************
 *
 * Save Group (including Views and child Groups) to a hdf5 handle
 *
 *************************************************************************
 */
bool Group::save(const hid_t& h5_id, const std::string& protocol, const Attribute* attr) const
{
  ConduitErrorSuppressor checkConduitCall(getDataStore());
  bool retval = false;

  // supported here:
  // "sidre_hdf5"
  // "conduit_hdf5"
  if(protocol == "sidre_hdf5")
  {
    Node n;
    exportTo(n["sidre"], attr);
    createExternalLayout(n["sidre/external"], attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::hdf5_write(n, h5_id); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else if(protocol == "conduit_hdf5")
  {
    Node n;
    createNativeLayout(n, attr);
    n["sidre_group_name"] = m_name;
    checkConduitCall([&] { conduit::relay::io::hdf5_write(n, h5_id); });
    retval = !(getDataStore()->getConduitErrorOccurred());
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol
                                       << "' for save with hdf5 handle.");
    retval = false;
  }

  return retval;
}

/*
 *************************************************************************
 *
 * Load Group (including Views and child Groups) from an hdf5 handle
 *
 *************************************************************************
 */
bool Group::load(const hid_t& h5_id, const std::string& protocol, bool preserve_contents)
{
  std::string name_from_file;
  return load(h5_id, protocol, preserve_contents, name_from_file);
}

bool Group::load(const hid_t& h5_id,
                 const std::string& protocol,
                 bool preserve_contents,
                 std::string& name_from_file)
{
  ConduitErrorSuppressor checkConduitCall(getDataStore());
  bool retval = false;

  // supported here:
  // "sidre_hdf5"
  // "conduit_hdf5"
  if(protocol == "sidre_hdf5")
  {
    Node n;
    checkConduitCall([&] { conduit::relay::io::hdf5_read(h5_id, n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      SLIC_ASSERT_MSG(n.has_path("sidre"),
                      SIDRE_GROUP_LOG_PREPEND << "Conduit Node " << n.path() << " does not have sidre "
                                              << "data for Group " << getPathName() << ".");
      importFrom(n["sidre"], preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else if(protocol == "conduit_hdf5")
  {
    SLIC_ERROR("Protocol " << protocol << " not yet supported for file load.");
    Node n;
    checkConduitCall([&] { conduit::relay::io::hdf5_read(h5_id, n); });
    if(!getDataStore()->getConduitErrorOccurred())
    {
      importConduitTree(n, preserve_contents);
      if(n.has_path("sidre_group_name"))
      {
        name_from_file = n["sidre_group_name"].as_string();
      }
      retval = true;
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Invalid protocol '" << protocol << "' for file load.");
    retval = false;
  }

  return retval;
}

/*
 *************************************************************************
 *
 * Load External Data from an hdf5 file
 *
 * Note: this ASSUMES usage of the "sidre_hdf5" protocol
 *************************************************************************
 */
bool Group::loadExternalData(const hid_t& h5_id)
{
  Node n;
  createExternalLayout(n);
  ConduitErrorSuppressor checkConduitCall(getDataStore());

  checkConduitCall([&] { conduit::relay::io::hdf5_read(h5_id, "sidre/external", n); });

  return !(getDataStore()->getConduitErrorOccurred());
}

/*
 *************************************************************************
 *
 * Load External Data from a path within an hdf5 file
 *
 * Note: this ASSUMES usage of the "sidre_hdf5" protocol
 *************************************************************************
 */
bool Group::loadExternalData(const hid_t& h5_id, const std::string& group_path)
{
  bool success;
  std::string delim(1, getPathDelimiter());
  if(group_path.empty() || group_path == delim)
  {
    // An empty or trivial path means the load should read from the start of
    // the file.
    success = loadExternalData(h5_id);
  }
  else
  {
    Node n;
    createExternalLayout(n);
    ConduitErrorSuppressor checkConduitCall(getDataStore());

    std::string read_path("sidre/external/" + group_path);

    checkConduitCall([&] { conduit::relay::io::hdf5_read(h5_id, read_path, n); });
    success = !(getDataStore()->getConduitErrorOccurred());
  }
  return success;
}

#endif /* AXOM_USE_HDF5 */

////////////////////////////////////////////////////////////////////////
//
// Private methods below
//
////////////////////////////////////////////////////////////////////////

/*
 *************************************************************************
 *
 * PRIVATE ctor makes Group with given name and make it a child of
 * root Group in datastore.
 *
 *************************************************************************
 */
Group::Group(const std::string& name, DataStore* datastore, bool is_list)
  : m_name(name)
  , m_index(InvalidIndex)
  , m_parent(nullptr)
  , m_datastore(datastore)
  , m_is_list(is_list)
  , m_view_coll(nullptr)
  , m_group_coll(nullptr)
#ifdef AXOM_USE_UMPIRE
  , m_default_array_alloc_id(axom::getDefaultAllocatorID())
  , m_default_tuple_alloc_id(axom::getDefaultAllocatorID())
#else
  , m_default_array_alloc_id(axom::MALLOC_ALLOCATOR_ID)
  , m_default_tuple_alloc_id(axom::MALLOC_ALLOCATOR_ID)
#endif
{
  if(is_list)
  {
    m_view_coll = new ListCollection<View>();
    m_group_coll = new ListCollection<Group>();
  }
  else
  {
    m_view_coll = new MapCollection<View>();
    m_group_coll = new MapCollection<Group>();
  }
}

/*
 *************************************************************************
 *
 * PRIVATE dtor destroys Group and all its contents.
 *
 *************************************************************************
 */
Group::~Group()
{
  destroyViews();
  destroyGroups();
  delete m_view_coll;
  delete m_group_coll;
}

/*
 *************************************************************************
 *
 * PRIVATE method to attach given View to Group.
 *
 *************************************************************************
 */
View* Group::attachView(View* view)
{
  if(view == nullptr || (!view->getName().empty() && hasChildView(view->getName())))
  {
    return nullptr;
  }
  else
  {
    SLIC_ASSERT_MSG(view->m_owning_group == nullptr,
                    SIDRE_GROUP_LOG_PREPEND
                      << "Provided View " << view->getPathName() << " is already "
                      << "attatched to Group " << view->m_owning_group->getPathName()
                      << " and can't be "
                      << "attatched to Group " << getPathName() << ".");
    view->m_owning_group = this;
    view->m_index = m_view_coll->insertItem(view, view->getName());
    return view;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach View with given name from Group.
 *
 *************************************************************************
 */
View* Group::detachView(const std::string& name)
{
  View* view = getNamedViews()->removeItem(name);
  if(view != nullptr)
  {
    view->m_owning_group = nullptr;
    view->m_index = InvalidIndex;
  }

  return view;
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach View with given index from Group.
 *
 *************************************************************************
 */
View* Group::detachView(IndexType idx)
{
  View* view = m_view_coll->removeItem(idx);
  if(view != nullptr)
  {
    view->m_owning_group = nullptr;
    view->m_index = InvalidIndex;
  }

  return view;
}

/*
 *************************************************************************
 *
 * PRIVATE method to destroy View in this Group and its data.
 *
 *************************************************************************
 */
void Group::destroyViewAndData(View* view)
{
  if(view != nullptr)
  {
    Group* group = view->getOwningGroup();
    group->detachView(view->getName());
    Buffer* const buffer = view->detachBuffer();
    if(buffer != nullptr && buffer->getNumViews() == 0)
    {
      getDataStore()->destroyBuffer(buffer);
    }
    delete view;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to make given Group a child of this Group.
 *
 *************************************************************************
 */
Group* Group::attachGroup(Group* group)
{
  if(group == nullptr || (!group->getName().empty() && hasChildGroup(group->getName())))
  {
    return nullptr;
  }
  else
  {
    group->m_parent = this;
    group->m_index = m_group_coll->insertItem(group, group->getName());
    return group;
  }
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach child Group with given name from Group.
 *
 *************************************************************************
 */
Group* Group::detachGroup(const std::string& name)
{
  Group* group = getNamedGroups()->removeItem(name);
  if(group != nullptr)
  {
    group->m_parent = nullptr;
    group->m_index = InvalidIndex;
  }

  return group;
}

/*
 *************************************************************************
 *
 * PRIVATE method to detach child Group with given index from Group.
 *
 *************************************************************************
 */
Group* Group::detachGroup(IndexType idx)
{
  Group* group = m_group_coll->removeItem(idx);
  if(group != nullptr)
  {
    group->m_parent = nullptr;
    group->m_index = InvalidIndex;
  }

  return group;
}

/*
 *************************************************************************
 *
 * Serialize tree identified by a Group into a conduit node.  Include
 * any Buffers attached to Views in that tree.
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */
bool Group::exportTo(conduit::Node& result, const Attribute* attr, bool export_buffer) const
{
  result.set(DataType::object());
  // TODO - This implementation will change in the future.  We want to write
  // out some separate set of conduit nodes:
  // #1 A set of nodes representing the Group and Views (hierarchy), with
  // the data descriptions ( schemas ).
  // #2 A set of nodes for our data ( Buffers, external data, etc ).
  // On a load, we want to be able to create our DataStore tree first,
  // then call allocate ourself, then have conduit load the data directly
  // into our allocated memory areas.  Conduit can do this, as long as the
  // conduit node set is compatible with what's in the file.
  std::set<IndexType> buffer_indices;

  // Tell Group to add itself and all sub-Groups and Views to node.
  // Any Buffers referenced by those Views will be tracked in the
  // buffer_indices
  bool hasSavedViews = exportTo(result, attr, buffer_indices);

  if(!buffer_indices.empty())
  {
    // Now, add all the referenced buffers to the node.
    Node& bnode = result["buffers"];
    for(std::set<IndexType>::iterator s_it = buffer_indices.begin(); s_it != buffer_indices.end();
        ++s_it)
    {
      // Use a dictionary layout here instead of conduit list.
      // Conduit IO HDF5 doesn't support conduit list objects.
      std::ostringstream oss;
      oss << "buffer_id_" << *s_it;
      Node& n_buffer = bnode.fetch(oss.str());
      if(export_buffer)
      {
        getDataStore()->getBuffer(*s_it)->exportTo(n_buffer);
      }
      else
      {
        getDataStore()->getBuffer(*s_it)->exportMetadata(n_buffer);
      }
    }
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * PRIVATE method to copy from Group to given Conduit node using
 * given set of ids to maintain correct association of data Buffers
 * to data Views.
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */
bool Group::exportTo(conduit::Node& result,
                     const Attribute* attr,
                     std::set<IndexType>& buffer_indices) const
{
  result.set(DataType::object());

  bool hasSavedViews = false;
  if(getNumViews() > 0)
  {
    Node& vnode = result["views"];
    for(const auto& view : this->views())
    {
      if(attr == nullptr || view.hasAttributeValue(attr))
      {
        Node& n_view = m_is_list ? vnode.append() : vnode.fetch(view.getName());
        view.exportTo(n_view, buffer_indices);
        hasSavedViews = true;
      }
    }
    if(!hasSavedViews)
    {
      result.remove("views");
    }
  }

  bool hasSavedGroups = false;
  if(getNumGroups() > 0)
  {
    Node& gnode = result["groups"];
    for(const auto& group : this->groups())
    {
      Node& n_group = m_is_list ? gnode.append() : gnode.fetch(group.getName());
      bool hsv = group.exportTo(n_group, attr, buffer_indices);
      hasSavedViews = hasSavedViews || hsv;
      hasSavedGroups = true;
    }
    if(!hasSavedGroups)
    {
      result.remove("groups");
    }
  }

  return hasSavedViews;
}

/*
 *************************************************************************
 *
 * Exports the contents of the Group, excluding the data arrays held
 * by any Buffers attached to the Views.
 *
 *************************************************************************
 */

bool Group::exportWithoutBufferData(conduit::Node& result, const Attribute* attr) const
{
  return exportTo(result, attr, false);
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into this Group.  Includes
 * any Buffers attached to Views in that tree.
 *
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */

void Group::importFrom(conduit::Node& node, bool preserve_contents)
{
  // TODO - May want to put in a little meta-data into these files like a
  // 'version'
  // or tag identifying the data.  We don't want someone giving us a file that
  // doesn't have our full multiView->buffer connectivity in there.

  if(!preserve_contents)
  {
    destroyGroups();
    destroyViews();
  }

  getDataStore()->loadAttributeLayout(node);

  // First - Import Buffers into the DataStore.
  std::map<IndexType, IndexType> buffer_indices_map;

  if(node.has_path("buffers"))
  {
    conduit::NodeIterator buffs_itr = node["buffers"].children();
    while(buffs_itr.has_next())
    {
      Node& n_buffer = buffs_itr.next();
      IndexType old_buffer_id = n_buffer["id"].to_int64();

      Buffer* buffer = getDataStore()->createBuffer();

      // track change of old Buffer id to new Buffer id
      buffer_indices_map[old_buffer_id] = buffer->getIndex();

      // populate the new Buffer's state
      buffer->importFrom(n_buffer);
    }
  }

  // Next - import tree of Groups, sub-Groups, Views into the DataStore.
  // Use the mapping of old to new Buffer ids to connect the Views to the
  // right Buffers.
  importFrom(node, buffer_indices_map);
}

/*
 *************************************************************************
 *
 * PRIVATE method to copy from given Conduit node to this Group using
 * given map of ids to indicate association of Buffer ids in node to
 * those in datastore.
 *
 * Note: This is for the "sidre_{zzz}" protocols.
 *
 *************************************************************************
 */
void Group::importFrom(conduit::Node& node, const std::map<IndexType, IndexType>& buffer_id_map)
{
  if(node.has_path("views"))
  {
    // create the Views
    conduit::NodeIterator views_itr = node["views"].children();
    while(views_itr.has_next())
    {
      Node& n_view = views_itr.next();
      std::string view_name = m_is_list ? "" : views_itr.name();

      View* view = createView(view_name);
      view->importFrom(n_view, buffer_id_map);
    }
  }
  if(node.has_path("groups"))
  {
    // create the child Groups
    conduit::NodeIterator groups_itr = node["groups"].children();
    while(groups_itr.has_next())
    {
      Group* group;
      Node& n_group = groups_itr.next();
      bool create_list = false;
      if(n_group.has_child("views"))
      {
        if(n_group["views"].dtype().is_list())
        {
          create_list = true;
        }
      }
      if(!create_list && n_group.has_child("groups"))
      {
        if(n_group["groups"].dtype().is_list())
        {
          create_list = true;
        }
      }
      std::string group_name = groups_itr.name();
      if(m_is_list)
      {
        group = createUnnamedGroup(create_list);
      }
      else
      {
        std::string group_name = groups_itr.name();
        group = createGroup(group_name, create_list);
      }
      group->importFrom(n_group, buffer_id_map);
    }
  }
}

/*
 *************************************************************************
 *
 * Imports tree from a conduit node into this Group.
 * This takes a generic conduit tree, not one with sidre conventions.
 *
 *************************************************************************
 */

bool Group::importConduitTree(const conduit::Node& node, bool preserve_contents)
{
  bool success = true;
  if(!preserve_contents)
  {
    destroyGroups();
    destroyViews();
  }

  //
  DataType node_dtype = node.dtype();
  if(node_dtype.is_object() || node_dtype.is_list())
  {
    conduit::NodeConstIterator itr = node.children();
    while(itr.has_next())
    {
      const Node& cld_node = itr.next();
      std::string cld_name = m_is_list ? "" : itr.name();
      DataType cld_dtype = cld_node.dtype();

      if(cld_dtype.is_object() || cld_dtype.is_list())
      {
        // create group
        Group* grp = m_is_list ? createUnnamedGroup(cld_dtype.is_list())
                               : createGroup(cld_name, cld_dtype.is_list());
        success = grp->importConduitTree(cld_node, preserve_contents);
      }
      else if(cld_dtype.is_empty())
      {
        //create empty view
        createView(cld_name);
      }
      else if(cld_dtype.is_string())
      {
        if(cld_name != "sidre_group_name")
        {
          //create string view
          createViewString(cld_name, cld_node.as_string());
        }
      }
      else if(cld_dtype.is_number())
      {
        if(cld_dtype.number_of_elements() == 1)
        {
          // create scalar view
          View* view = createView(cld_name);
          view->setScalar(cld_node);
        }
        else
        {
          View* view = createView(cld_name);
          view->importArrayNode(cld_node);
        }
      }
      else
      {
        // All Nodes should have one of the above datatypes, so if
        // we get here something is wrong.
        SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Conduit child Node " << cld_name
                                           << " does not have a recognized datatype."
                                           << " Cannot import into Group " << getPathName());
      }
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Group cannot import non-object Conduit Node");
  }

  return success;
}

bool Group::importConduitTreeExternal(conduit::Node& node, bool preserve_contents)
{
  bool success = true;
  if(!preserve_contents)
  {
    destroyGroups();
    destroyViews();
  }

  DataType node_dtype = node.dtype();
  if(node_dtype.is_object() || node_dtype.is_list())
  {
    conduit::NodeIterator itr = node.children();
    while(itr.has_next())
    {
      Node& cld_node = itr.next();
      std::string cld_name = m_is_list ? "" : itr.name();
      DataType cld_dtype = cld_node.dtype();

      if(cld_dtype.is_object() || cld_dtype.is_list())
      {
        // create group
        Group* grp = m_is_list ? createUnnamedGroup(cld_dtype.is_list())
                               : createGroup(cld_name, cld_dtype.is_list());
        success = grp->importConduitTreeExternal(cld_node, preserve_contents);
      }
      else if(cld_dtype.is_empty())
      {
        //create empty view
        createView(cld_name);
      }
      else if(cld_dtype.is_string())
      {
        if(cld_name != "sidre_group_name")
        {
          //create string view
          // createViewString(cld_name, cld_node.data_ptr());
          int allocId = ConduitMemory::conduitAllocIdToAxom(cld_node.allocator());
          if(allocId == axom::INVALID_ALLOCATOR_ID)
          {
            allocId = axom::MALLOC_ALLOCATOR_ID;
          }
          createViewString(cld_name, cld_node.as_string(), allocId);
        }
      }
      else if(cld_dtype.is_number())
      {
        if(cld_dtype.number_of_elements() == 1)
        {
          // create scalar view
          int allocId = ConduitMemory::conduitAllocIdToAxom(cld_node.allocator());
          if(allocId == axom::INVALID_ALLOCATOR_ID)
          {
            allocId = axom::MALLOC_ALLOCATOR_ID;
          }
          View* view = createView(cld_name);
          view->setScalar(cld_node, allocId);
        }
        else
        {
          void* conduit_ptr = cld_node.data_ptr();
          View* view = createView(cld_name);
          view->setExternalDataPtr(conduit_ptr);
          view->apply(cld_dtype);
        }
      }
      else
      {
        // All Nodes should have one of the above datatypes, so if
        // we get here something is wrong.
        SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Conduit child Node " << cld_name
                                           << " does not have a recognized datatype."
                                           << " Cannot import into Group " << getPathName());
      }
    }
  }
  else
  {
    SLIC_ERROR(SIDRE_GROUP_LOG_PREPEND << "Group cannot import non-object Conduit Node");
  }

  return success;
}

/*
 *************************************************************************
 *
 * PRIVATE method to walk down a path to the next-to-last entry.
 *
 * If an error is encountered, this private function will return nullptr
 *
 *************************************************************************
 */
Group* Group::walkPath(std::string& path, bool create_groups_in_path)
{
  Group* group_ptr = this;

  // Split path into parts
  std::vector<std::string> path_parts = axom::Path(path, s_path_delimiter).parts();

  if(path_parts.size() > 0)
  {
    if(m_is_list && path_parts.size() > 1)
    {
      // A size > 1 indicates that a delimited path string was provided.
      // A path string is invalid when this Group uses the list format.
      SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND << "A delimited path string '" << path
                                           << "' cannot be used as the name of an object "
                                           << "to be created by a Group that uses the list format. "
                                           << "A null pointer will be returned.");
      group_ptr = nullptr;
    }
    else
    {
      // Find stopping point (right before last part of path)
      std::vector<std::string>::const_iterator stop = path_parts.end() - 1;

      // Navigate path down to desired Group
      for(std::vector<std::string>::const_iterator iter = path_parts.begin(); iter < stop; ++iter)
      {
        if(group_ptr->hasChildGroup(*iter))
        {
          group_ptr = group_ptr->getGroup(*iter);
        }
        else if(create_groups_in_path)
        {
          group_ptr = group_ptr->createGroup(*iter);

          if(group_ptr == nullptr)
          {
            iter = stop;
          }
        }
        else
        {
          iter = stop;
          group_ptr = nullptr;
        }
      }
      path = path_parts.back();
    }
  }

  return group_ptr;
}

/*
 *************************************************************************
 *
 * PRIVATE const method to walk down a path to the next-to-last entry.
 *
 * If an error is encountered, this private function will return nullptr
 *
 *************************************************************************
 */
const Group* Group::walkPath(std::string& path) const
{
  const Group* group_ptr = this;

  // Split path into parts
  std::vector<std::string> path_parts = axom::Path(path, s_path_delimiter).parts();

  if(path_parts.size() > 0)
  {
    // Find stopping point (right before last part of path)
    std::vector<std::string>::const_iterator stop = path_parts.end() - 1;

    // Navigate path down to desired Group
    for(std::vector<std::string>::const_iterator iter = path_parts.begin(); iter < stop; ++iter)
    {
      if(group_ptr->hasChildGroup(*iter))
      {
        group_ptr = group_ptr->getGroup(*iter);
      }
      else
      {
        group_ptr = nullptr;
        iter = stop;
      }
    }
    path = path_parts.back();
  }

  return group_ptr;
}

/*
 *************************************************************************
 *
 * Return number of child Groups in a Group object.
 *
 *************************************************************************
 */
IndexType Group::getNumGroups() const { return m_group_coll->getNumItems(); }

/*
 *************************************************************************
 *
 * Return number of Views owned by a Group object.
 *
 *************************************************************************
 */
IndexType Group::getNumViews() const { return m_view_coll->getNumItems(); }

/*
 *************************************************************************
 *
 * Return true if this Group owns a View with given name (not path);
 * else false.
 *
 *************************************************************************
 */
bool Group::hasChildView(const std::string& name) const
{
  return isUsingMap() ? getNamedViews()->hasItem(name) : false;
}

/*
 *************************************************************************
 *
 * Return true if this Group owns a View with given index; else false.
 *
 *************************************************************************
 */
bool Group::hasView(IndexType idx) const { return m_view_coll->hasItem(idx); }

/*
 *************************************************************************
 *
 * Return index of View with given name owned by this Group object.
 *
 * If no such View exists, return sidre::InvalidIndex;
 *
 *************************************************************************
 */
IndexType Group::getViewIndex(const std::string& name) const
{
  SLIC_CHECK_MSG(hasChildView(name),
                 SIDRE_GROUP_LOG_PREPEND << "Group has no View with name '" << name << "'");

  return getNamedViews()->getItemIndex(name);
}

/*
 *************************************************************************
 *
 * Return name of View with given index owned by Group object.
 *
 * If no such View exists, return sidre::InvalidName.
 *
 *************************************************************************
 */
const std::string& Group::getViewName(IndexType idx) const
{
  SLIC_CHECK_MSG(hasView(idx), SIDRE_GROUP_LOG_PREPEND << "Group has no View with index " << idx);

  return getNamedViews()->getItemName(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to non-const View with given index.
 *
 * If no such View exists, nullptr is returned.
 *
 *************************************************************************
 */
View* Group::getView(IndexType idx)
{
  SLIC_CHECK_MSG(hasView(idx), SIDRE_GROUP_LOG_PREPEND << "Group has no View with index " << idx);

  return m_view_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to const View with given index.
 *
 * If no such View exists, nullptr is returned.
 *
 *************************************************************************
 */
const View* Group::getView(IndexType idx) const
{
  SLIC_CHECK_MSG(hasView(idx), SIDRE_GROUP_LOG_PREPEND << "Group has no View with index " << idx);

  return m_view_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return first valid View index in Group object
 *        (i.e., smallest index over all Views).
 *
 * sidre::InvalidIndex is returned if Group has no Views.
 *
 *************************************************************************
 */
IndexType Group::getFirstValidViewIndex() const { return m_view_coll->getFirstValidIndex(); }

/*
 *************************************************************************
 *
 * Return next valid View index in Group object after given index
 *        (i.e., smallest index over all View indices larger than given one).
 *
 * sidre::InvalidIndex is returned if there is no valid index greater
 * than given one.
 *
 *************************************************************************
 */
IndexType Group::getNextValidViewIndex(IndexType idx) const
{
  return m_view_coll->getNextValidIndex(idx);
}

/*!
 * \brief Returns an adaptor to support iterating the collection of views
 */
typename Group::ViewCollection::iterator_adaptor Group::views()
{
  return m_view_coll->getIteratorAdaptor();
}

/*!
 * \brief Returns a const adaptor to support iterating the collection of views
 */
typename Group::ViewCollection::const_iterator_adaptor Group::views() const
{
  return m_view_coll->getIteratorAdaptor();
}

/*!
 * \brief Returns an adaptor to support iterating the collection of groups
 */
typename Group::GroupCollection::iterator_adaptor Group::groups()
{
  return m_group_coll->getIteratorAdaptor();
}

/*!
 * \brief Returns a const adaptor to support iterating the collection of groups
 */
typename Group::GroupCollection::const_iterator_adaptor Group::groups() const
{
  return m_group_coll->getIteratorAdaptor();
}

/*
 *************************************************************************
 *
 * Return true if this Group has a child Group with given
 * name; else false.
 *
 *************************************************************************
 */
bool Group::hasChildGroup(const std::string& name) const
{
  return isUsingMap() ? getNamedGroups()->hasItem(name) : false;
}

/*
 *************************************************************************
 *
 * Return true if Group has an immediate child Group
 * with given index; else false.
 *
 *************************************************************************
 */
bool Group::hasGroup(IndexType idx) const { return m_group_coll->hasItem(idx); }

/*
 *************************************************************************
 * Return the index of immediate child Group with given name.
 *
 * If no such child Group exists, return sidre::InvalidIndex;
 *
 *************************************************************************
 */
IndexType Group::getGroupIndex(const std::string& name) const
{
  SLIC_CHECK_MSG(hasChildGroup(name),
                 SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with name '" << name << "'");

  return getNamedGroups()->getItemIndex(name);
}

/*
 *************************************************************************
 *
 * Return the name of immediate child Group with given index.
 *
 * If no such child Group exists, return sidre::InvalidName.
 *
 *************************************************************************
 */
const std::string& Group::getGroupName(IndexType idx) const
{
  SLIC_CHECK_MSG(hasGroup(idx),
                 SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with index " << idx);

  return getNamedGroups()->getItemName(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to non-const immediate child Group with given index.
 *
 * If no such Group exists, nullptr is returned.
 *
 *************************************************************************
 */
Group* Group::getGroup(IndexType idx)
{
  SLIC_CHECK_MSG(hasGroup(idx),
                 SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with index " << idx);

  return m_group_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return pointer to const immediate child Group with given index.
 *
 * If no such Group exists, nullptr is returned.
 *
 *************************************************************************
 */
const Group* Group::getGroup(IndexType idx) const
{
  SLIC_CHECK_MSG(hasGroup(idx),
                 SIDRE_GROUP_LOG_PREPEND << "Group has no child Group with index " << idx);

  return m_group_coll->getItem(idx);
}

/*
 *************************************************************************
 *
 * Return first valid child Group index (i.e., smallest
 *        index over all child Groups).
 *
 * sidre::InvalidIndex is returned if Group has no child Groups.
 *
 *************************************************************************
 */
IndexType Group::getFirstValidGroupIndex() const { return m_group_coll->getFirstValidIndex(); }

/*
 *************************************************************************
 *
 * Return next valid child Group index after given index
 *        (i.e., smallest index over all child Group indices larger
 *        than given one).
 *
 * sidre::InvalidIndex is returned if there is no valid index greater
 * than given one.
 *
 *************************************************************************
 */
IndexType Group::getNextValidGroupIndex(IndexType idx) const
{
  return m_group_coll->getNextValidIndex(idx);
}

/*
 *************************************************************************
 *
 * Rename this Group with a new string name.
 *
 *************************************************************************
 */
bool Group::rename(const std::string& new_name)
{
  bool do_rename = true;
  if(new_name != m_name)
  {
    if(new_name.empty())
    {
      SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND << "Cannot rename Group to an empty string.");
      do_rename = false;
    }
    else if(new_name.find(s_path_delimiter) != std::string::npos)
    {
      SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND << "Cannot rename Group to path name '" << new_name << "'. "
                                           << "Only strings without path delimiters can "
                                           << "be passed into the rename method.");
      do_rename = false;
    }

    if(do_rename)
    {
      const Group* root = getDataStore()->getRoot();
      Group* parent = getParent();

      //If this is the root group, we don't need
      //to do anything to change the parent's handle to this group.
      if(this != root && parent != nullptr)
      {
        if(parent->hasGroup(new_name) || parent->hasView(new_name))
        {
          SLIC_WARNING(SIDRE_GROUP_LOG_PREPEND << "Parent group " << parent->getPathName()
                                               << " already has a child group named '" << new_name
                                               << "'. Group " << getPathName()
                                               << " will not be renamed.");
          do_rename = false;
        }
        else
        {
          Group* detached_group = parent->detachGroup(m_name);
          SLIC_CHECK_MSG(detached_group == this,
                         SIDRE_GROUP_LOG_PREPEND << "Group detatched from parent '"
                                                 << detached_group->getPathName() << "' is not "
                                                 << "this Group '" << getPathName() << "'.");

          m_name = new_name;

          Group* attached_group = parent->attachGroup(detached_group);
          AXOM_UNUSED_VAR(attached_group);
          SLIC_CHECK_MSG(attached_group == this,
                         SIDRE_GROUP_LOG_PREPEND << "Group attached to parent '"
                                                 << attached_group->getPathName() << "' is not "
                                                 << "this Group '" << getPathName() << "'.");
        }
      }
      else
      {
        m_name = new_name;
      }
    }
  }

  return do_rename;
}

} /* end namespace sidre */
} /* end namespace axom */
