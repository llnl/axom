// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_ConnectivityArrayHelpers_HPP_
#define MINT_ConnectivityArrayHelpers_HPP_

#include "axom/config.hpp"
#include "axom/core/Array.hpp"

#include "axom/mint/config.hpp"
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/utils/ArrayWrapper.hpp"

#include "axom/slic/interface/slic.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

#include <cmath>   /* for std::ceil */
#include <string>  /* for std:string */
#include <cstring> /* for std::srcmp */

namespace axom
{
namespace mint
{
namespace internal
{
#ifdef AXOM_MINT_USE_SIDRE

/*!
 * \brief Initializes the members of a ConnectivityArray instance from a
 *  sidre::Group which already has data.
 *
 * \param [in] group the sidre::Group to create the ConnectivityArray from.
 * \param [out] m_values a pointer to a pointer to the values array.
 * \param [out] m_offsets a pointer to a pointer to the offsets array.
 * \param [out] m_types a pointer to a pointer to the types array.
 *
 * \note if the offset/types pointers aren't given then this method does not
 *  check that the appropriate views exist.
 * \note the given Group must conform to a single Blueprint topology.
 *
 * \pre group != nullptr
 * \pre m_values != nullptr
 */
inline CellType initializeFromGroup(sidre::Group* group,
                                    detail::ArrayWrapper<IndexType>& m_values,
                                    detail::ArrayWrapper<IndexType>& m_offsets,
                                    detail::ArrayWrapper<CellType>& m_types)
{
  SLIC_ERROR_IF(group == nullptr, "sidre::Group pointer must not be null.");

  SLIC_ERROR_IF(
    !group->hasView("coordset"),
    "sidre::Group " << group->getPathName()
                    << " does not conform to mesh blueprint. No child view 'coordset'.");
  SLIC_ERROR_IF(!group->getView("coordset")->isString(),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. Child view 'coordset' "
                                << "does not hold a string.");

  SLIC_ERROR_IF(!group->hasView("type"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. No child view 'type'.");
  sidre::View* type_view = group->getView("type");
  SLIC_ERROR_IF(!type_view->isString(),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. Child view 'type' "
                                << "does not hold a string.");
  SLIC_ERROR_IF(
    std::strcmp(type_view->getString(), "unstructured") != 0,
    "Incorrect type found. Expected 'unstructured' but got '" << type_view->getString() << "'.");

  SLIC_ERROR_IF(!group->hasGroup("elements"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. No 'elements' group "
                                << "found.");
  sidre::Group* elems_group = group->getGroup("elements");

  SLIC_ERROR_IF(!elems_group->hasView("shape"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. The elements group "
                                << "does not have a child view 'shape'.");

  sidre::View* shape_view = elems_group->getView("shape");
  std::string bp_name = shape_view->getString();
  CellType cell_type = UNDEFINED_CELL;

  for(IndexType i = 0; i < NUM_CELL_TYPES; ++i)
  {
    if(cell_info[i].blueprint_name == bp_name)
    {
      cell_type = cell_info[i].cell_type;
      break;
    }
  }

  SLIC_ERROR_IF(!elems_group->hasView("connectivity"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. The elements group "
                                << "does not have a child view 'connectivity'.");

  sidre::View* connec_view = elems_group->getView("connectivity");
  m_values = sidre::Array<IndexType>(connec_view);

  {
    SLIC_ERROR_IF(!elems_group->hasView("offsets"),
                  "sidre::Group " << group->getPathName() << " does not conform to mesh blueprint.");

    sidre::View* offsets_view = elems_group->getView("offsets");
    m_offsets = sidre::Array<IndexType>(offsets_view);

    if(m_offsets.size() == 0)
    {
      m_offsets.push_back(0);
    }
    SLIC_ERROR_IF((*m_offsets)[0] != 0,
                  "The offset of ID 0 must be 0 not " << (*m_offsets)[0] << ".");
  }

  {
    SLIC_ERROR_IF(!elems_group->hasView("types"),
                  "sidre::Group " << group->getPathName() << " does not conform to mesh blueprint.");

    sidre::View* elem_type_view = elems_group->getView("types");
    m_types = sidre::Array<CellType>(elem_type_view);
  }

  return cell_type;
}
/*!
 * \brief Initializes the members of a ConnectivityArray instance from a
 *  sidre::Group which already has data.
 *
 * \param [in] group the sidre::Group to create the ConnectivityArray from.
 * \param [out] m_values a pointer to a pointer to the values array.
 *
 * \note if the offset/types pointers aren't given then this method does not
 *  check that the appropriate views exist.
 * \note the given Group must conform to a single Blueprint topology.
 *
 * \pre group != nullptr
 * \pre m_values != nullptr
 */
inline CellType initializeFromGroup(sidre::Group* group, detail::ArrayWrapper<IndexType, 2>& m_values)
{
  SLIC_ERROR_IF(group == nullptr, "sidre::Group pointer must not be null.");

  SLIC_ERROR_IF(
    !group->hasView("coordset"),
    "sidre::Group " << group->getPathName()
                    << " does not conform to mesh blueprint. No child view 'coordset'.");
  SLIC_ERROR_IF(!group->getView("coordset")->isString(),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. Child view 'coordset' "
                                << "does not hold a string.");

  SLIC_ERROR_IF(!group->hasView("type"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. No child view 'type'.");
  sidre::View* type_view = group->getView("type");
  SLIC_ERROR_IF(!type_view->isString(),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. Child view 'type' "
                                << "does not hold a string.");
  SLIC_ERROR_IF(
    std::strcmp(type_view->getString(), "unstructured") != 0,
    "Incorrect type found. Expected 'unstructured' but got '" << type_view->getString() << "'.");

  SLIC_ERROR_IF(!group->hasGroup("elements"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. No 'elements' group "
                                << "found.");
  sidre::Group* elems_group = group->getGroup("elements");

  SLIC_ERROR_IF(!elems_group->hasView("shape"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. The elements group "
                                << "does not have a child view 'shape'.");

  sidre::View* shape_view = elems_group->getView("shape");
  std::string bp_name = shape_view->getString();
  CellType cell_type = UNDEFINED_CELL;

  for(IndexType i = 0; i < NUM_CELL_TYPES; ++i)
  {
    if(cell_info[i].blueprint_name == bp_name)
    {
      cell_type = cell_info[i].cell_type;
      break;
    }
  }

  SLIC_ERROR_IF(!elems_group->hasView("connectivity"),
                "sidre::Group " << group->getPathName()
                                << " does not conform to mesh blueprint. The elements group "
                                << "does not have a child view 'connectivity'.");

  sidre::View* connec_view = elems_group->getView("connectivity");
  m_values = sidre::MCArray<IndexType>(connec_view);
  return cell_type;
}

/*!
 * \brief Initializes an empty sidre::Group to hold a ConnectivityArray.
 *
 * \param [out] group the empty sidre::Group.
 * \param [in] coordset the name of the Blueprint coordinate set to associate
 *  this ConnectivityArray with.
 * \param [in] create_offsets iff true will create a view for the offsets array.
 * \param [in] create_types iff true will create a view for the types array.
 *
 * \pee group != nullptr
 * \pre group->getNumGroups() == group->getNumViews() == 0
 */
inline void initializeGroup(sidre::Group* group,
                            const std::string& coordset,
                            CellType cell_type,
                            bool create_offsets = false,
                            bool create_types = false)
{
  SLIC_ERROR_IF(group == nullptr, "sidre::Group pointer must not be null.");
  SLIC_ERROR_IF(group->getNumGroups() != 0, "sidre::Group is not empty.");
  SLIC_ERROR_IF(group->getNumViews() != 0, "sidre::Group is not empty.");

  group->createView("coordset")->setString(coordset);
  group->createView("type")->setString("unstructured");

  const std::string bp_name =
    (cell_type == UNDEFINED_CELL) ? "mixed" : getCellInfo(cell_type).blueprint_name;

  sidre::Group* elems_group = group->createGroup("elements");
  elems_group->createView("shape")->setString(bp_name);

  elems_group->createView("connectivity");

  if(create_offsets)
  {
    elems_group->createView("offsets");
  }

  if(create_types)
  {
    elems_group->createView("types");
  }
}

/*!
 * \brief Sets the stride associated with a connectivity array.
 *
 * \param [out] group the group holding the connectivity array.
 * \param [in] stride the stride to set.
 *
 * \pre group != nullptr
 */
inline void setStride(sidre::Group* group, IndexType stride)
{
  SLIC_ERROR_IF(group == nullptr, "sidre::Group pointer must not be null.");

  sidre::Group* elems_group = group->getGroup("elements");
  SLIC_ERROR_IF(elems_group == nullptr, "No group found");
  elems_group->createView("stride")->setScalar(stride);
}

/*!
 * \brief Returns the stride associated with a connectivity array.
 *
 * \param [in] group the group holding the connectivity array.
 *
 * \pee group != nullptr
 */
inline IndexType getStride(const sidre::Group* group)
{
  SLIC_ERROR_IF(group == nullptr, "sidre::Group pointer must not be null.");

  const sidre::Group* elems_group = group->getGroup("elements");
  SLIC_ERROR_IF(elems_group == nullptr, "No group found");
  const sidre::View* stride_view = elems_group->getView("stride");
  SLIC_ERROR_IF(stride_view == nullptr, "No view found");
  return stride_view->getScalar();
}

#endif

/*!
 * \brief Append multiple IDs to the members of a ConnectivityArray.
 *
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] values pointer to the values to append.
 * \param [in] offsets the offsets array of length at least n_IDs + 1.
 * \param [in/out] m_values a pointer to the values array.
 * \param [in/out] m_offsets a pointer to the offsets array.
 *
 * \note The number of values to append is given by
 *  offsets[n_IDs + 1] - offsets[0] and the values array must be at least
 *  this long.
 *
 * \pre n_IDs >= 0
 * \pre values != nullptr
 * \pre offsets != nullptr
 * \pre m_values != nullptr
 * \pre m_offsets != nullptr
 */
inline void append(IndexType n_IDs,
                   const IndexType* values,
                   const IndexType* offsets,
                   detail::ArrayWrapper<IndexType>& m_values,
                   detail::ArrayWrapper<IndexType>& m_offsets)
{
  SLIC_ASSERT(values != nullptr);
  SLIC_ASSERT(offsets != nullptr);

  IndexType n_values_to_add = offsets[n_IDs] - offsets[0];
  IndexType old_n_values = m_values.size();
  IndexType old_n_offsets = m_offsets.size();

  m_offsets.insert(m_offsets.size(), n_IDs, offsets + 1);
  m_values.insert(m_values.size(), n_values_to_add, values);

  /* Correct the appended offsets. */
  IndexType* m_offsets_ptr = m_offsets.data();
  const IndexType correction = old_n_values - offsets[0];
  for(IndexType i = 0; i < n_IDs; ++i)
  {
    m_offsets_ptr[old_n_offsets + i] += correction;
  }
}

/*!
 * \brief Sets the values of multiple IDs starting with the given ID.
 *
 * \param [in] start_ID the ID to start at.
 * \param [in] values pointer to the values to set, of length at least
 *  the sum of the number of values of each ID to set.
 * \param [in] n_IDs the number of IDs to set.
 * \param [in/out] m_values a pointer to the values array.
 * \param [in] m_offsets a pointer to the offsets axom::deprecated::MCArray.
 *
 * \pre start_ID >= 0 && start_ID + n_IDs < getNumberOfIDs()
 * \pre values != nullptr
 * \pre m_values != nullptr
 */
inline void set(IndexType start_ID,
                const IndexType* values,
                IndexType n_IDs,
                detail::ArrayWrapper<IndexType>& m_values,
                detail::ArrayWrapper<IndexType>& m_offsets)
{
  SLIC_ASSERT(start_ID >= 0);
  SLIC_ASSERT(start_ID + n_IDs <= m_offsets.size() - 1);
  SLIC_ASSERT(values != nullptr);

  IndexType offset = (*m_offsets)[start_ID];
  IndexType n_values = (*m_offsets)[start_ID + n_IDs] - offset;
  m_values.set(values, n_values, offset);
}

/*!
 * \brief Insert the values of new IDs before the given ID.
 *
 * \param [in] start_ID the ID to insert at.
 * \param [in] n_IDs the number of IDs to insert.
 * \param [in] values pointer to the values to insert.
 * \param [in] offsets the offsets array of length at least n_IDs + 1.
 * \param [in/out] m_values a pointer to the values axom::deprecated::MCArray.
 * \param [in/out] m_offsets a pointer to the offsets axom::deprecated::MCArray.
 *
 * \note The number of values to insert is given by
 *  offsets[n_IDs + 1] - offsets[0] and the values array must be at least
 *  this long.
 *
 * \pre start_ID >= 0 && start_ID <= getNumberOfIDs()
 * \pre n_IDs >= 0
 * \pre values != nullptr
 * \pre offsets != nullptr
 * \pre m_values != nullptr
 * \pre m_offsets != nullptr
 */
inline void insert(IndexType start_ID,
                   IndexType n_IDs,
                   const IndexType* values,
                   const IndexType* offsets,
                   detail::ArrayWrapper<IndexType>& m_values,
                   detail::ArrayWrapper<IndexType>& m_offsets)
{
  SLIC_ASSERT(start_ID >= 0);
  SLIC_ASSERT(start_ID <= m_offsets.size() - 1);
  SLIC_ASSERT(n_IDs >= 0);
  SLIC_ASSERT(values != nullptr);
  SLIC_ASSERT(offsets != nullptr);

  IndexType n_values = offsets[n_IDs] - offsets[0];
  IndexType* m_offsets_ptr = m_offsets.data();
  IndexType insert_pos = m_offsets_ptr[start_ID];

  /* Increment the offsets after the insertion position. */
  IndexType n_offsets = m_offsets.size();
  for(IndexType i = start_ID + 1; i < n_offsets; ++i)
  {
    m_offsets_ptr[i] += n_values;
  }

  m_offsets.insert(start_ID + 1, n_IDs, offsets + 1);
  m_values.insert(insert_pos, n_values, values);

  /* Correct the inserted offsets. */
  m_offsets_ptr = m_offsets.data();
  const IndexType correction = insert_pos - offsets[0];
  for(IndexType i = 0; i < n_IDs; ++i)
  {
    m_offsets_ptr[start_ID + 1 + i] += correction;
  }
}

/*!
 * \brief Return the value capacity given the number of IDs, the ID_capacity,
 *  and the number of values.
 *
 * \param [in] n_IDs the number of IDs in the ConnectivityArray.
 * \param [in] ID_capacity the ID capacity of the ConnectivityArray.
 * \param [in] n_values the number of values in the ConnectivityArray.
 * \param [in] value_capacity the value capacity of the ConnectivityArray.
 *  If this is not USE_DEFAULT it simply returns the given capacity, otherwise
 *  it calculates the capacity given the other three parameters.
 */
inline IndexType calcValueCapacity(IndexType n_IDs,
                                   IndexType ID_capacity,
                                   IndexType n_values,
                                   IndexType value_capacity)
{
  if(value_capacity == USE_DEFAULT)
  {
    if(n_IDs == 0)
    {
      value_capacity = ID_capacity * MAX_CELL_NODES;
    }
    else
    {
      const double avg_n_vals = double(n_values) / n_IDs;
      value_capacity = static_cast<IndexType>(std::ceil(avg_n_vals * n_IDs));
    }
  }

  return value_capacity;
}

} /* namespace internal */
} /* namespace mint */
} /* namespace axom */

#endif /* MINT_ConnectivityArrayHelpers_HPP_ */
