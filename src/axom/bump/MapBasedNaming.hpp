// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_MAP_BASED_NAMING_HPP_
#define AXOM_BUMP_MAP_BASED_NAMING_HPP_

#include "axom/bump/utilities/utilities.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <map>
#include <set>

namespace axom
{
namespace bump
{

//------------------------------------------------------------------------------
/*!
 * \brief This class implements a naming policy that uses a map to keep track of
 *        naming groups of ids. It is NOT parallel or device-enabled. It can
 *        be used for serial cases or serial debugging.
 *
 * \tparam IndexT The index type that gets hashed.
 * \tparam MAXIDS The max number of ids that get hashed.
 */
template <typename IndexT, int MAXIDS = 14>
class MapBasedNaming
{
public:
  using KeyType = std::uint64_t;
  using IndexType = IndexT;

private:
  using MapType = std::map<std::set<IndexT>, KeyType>;

public:
  /*!
   * \brief A view for making names.
   */
  class View
  {
  public:
    using KeyType = MapBasedNaming::KeyType;

    /*!
     * \brief Make a name from an array of ids.
     *
     * \param p The array of ids.
     * \param n The number of ids in the array.
     *
     * \return The name that describes the array of ids.
     */
    AXOM_HOST_DEVICE
    KeyType makeName(const IndexType *p, int n) const
    {
      std::set<IndexType> ids;
      for(int i = 0; i < n; i++)
      {
        ids.insert(p[i]);
      }
      KeyType name {};
      const auto it = m_map_ptr->find(ids);
      if(it == m_map_ptr->end())
      {
        name = static_cast<KeyType>(m_map_ptr->size());
        m_map_ptr->operator[](ids) = name;
      }
      else
      {
        name = it->second;
      }
      return name;
    }

    /// Set the max number of nodes.
    AXOM_HOST_DEVICE
    void setMaxId(IndexType) { }

    MapType *m_map_ptr {nullptr};
  };

  // Host-callable methods

  /// Make a name from the array of ids.
  KeyType makeName(const IndexType *p, int n) const { return m_view.makeName(p, n); }

  /*!
   * \brief Set the max number of nodes, which can help with id packing/narrowing.
   * \param n The number of nodes.
   */
  void setMaxId(IndexType n) { m_view.setMaxId(n); }

  /// Return a view that can be used on device.
  View view()
  {
    m_view.m_map_ptr = &m_map;
    return m_view;
  }

  /*!
   * \brief Turn name into a string.
   * \param key The name.
   * \return A string that represents the name.
   */
  static std::string toString(KeyType key)
  {
    std::stringstream ss;
    ss << key;
    return ss.str();
  }

  View m_view {};
  MapType m_map;
};

}  // end namespace bump
}  // end namespace axom

#endif
