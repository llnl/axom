// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_VIEWS_MATERIAL_VIEW_HPP_
#define AXOM_BUMP_VIEWS_MATERIAL_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

#include <algorithm>
#include <string>
#include <map>
#include <vector>

namespace axom
{
namespace bump
{
namespace views
{
/*!
 * \brief This object contains information about the materials as provided by a Conduit node.
 *
 * \note This would only be used on the host.
 */
struct Material
{
  int number {};
  std::string name {};
};

using MaterialInformation = std::vector<Material>;

/*!
 * \brief Return a vector of Material from a matset (this is the material_map)
 *
 * \param matset The Conduit node that contains the matset.
 *
 * \return A vector of Material that contains the materials in the material_map.
 */
MaterialInformation materials(const conduit::Node &matset);

//---------------------------------------------------------------------------
// Material views - These objects are meant to wrap Blueprint Matsets behind
//                  an interface that lets us query materials for a single
//                  zone. It is intended that these views will be used in
//                  device kernels.
//---------------------------------------------------------------------------

/*!
 \brief Material view for unibuffer matsets.

 \tparam IndexT The integer type used for material data.
 \tparam FloatT The floating point type used for material data (volume fractions).
 \tparam MAXMATERIALS The maximum number of materials to support.

 \verbatim

matsets:
  matset:
    topology: topology
    material_map:
      a: 1
      b: 2
      c: 0
    material_ids: [0, 1, 2, 2, 2, 0, 1, 0]
    volume_fractions: [0, a0, b2, b1, b0, 0, a1, 0]
    sizes: [2, 2, 1]
    offsets: [0, 2, 4]
    indices: [1, 4, 6, 3, 2]

 \endverbatim
 */
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS = 20>
class UnibufferMaterialView
{
public:
  using MaterialID = IndexT;
  using ZoneIndex = IndexT;
  using IndexType = IndexT;
  using FloatType = FloatT;
  using IDList = StaticArray<MaterialID, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static axom::IndexType MaxMaterials = MAXMATERIALS;

  void set(const axom::ArrayView<IndexType> &material_ids,
           const axom::ArrayView<FloatType> &volume_fractions,
           const axom::ArrayView<IndexType> &sizes,
           const axom::ArrayView<IndexType> &offsets,
           const axom::ArrayView<IndexType> &indices)
  {
#if !defined(AXOM_DEVICE_CODE)
    SLIC_ERROR_IF(material_ids.size() != volume_fractions.size(),
                  "Array views for material_ids, volume_fractions are different sizes.");
    SLIC_ERROR_IF(sizes.size() != offsets.size(),
                  "Array views for sizes, offsets are different sizes.");
#endif
    m_material_ids = material_ids;
    m_volume_fractions = volume_fractions;
    m_sizes = sizes;
    m_offsets = offsets;
    m_indices = indices;
  }

  AXOM_HOST_DEVICE
  inline axom::IndexType numberOfZones() const { return m_sizes.size(); }

  AXOM_HOST_DEVICE
  inline axom::IndexType numberOfMaterials(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));
    return m_sizes[zi];
  }

  AXOM_HOST_DEVICE
  void zoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    ids.clear();
    vfs.clear();

    const auto sz = numberOfMaterials(zi);
    const auto offset = m_offsets[zi];
    for(axom::IndexType i = 0; i < sz; i++)
    {
      const auto idx = m_indices[offset + i];

      ids.push_back(m_material_ids[idx]);
      vfs.push_back(m_volume_fractions[idx]);
    }
  }

  AXOM_HOST_DEVICE
  axom::IndexType zoneMaterials(ZoneIndex zi,
                                axom::ArrayView<IndexType> &ids,
                                axom::ArrayView<FloatType> &vfs) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    const auto sz = numberOfMaterials(zi);
    SLIC_ASSERT(sz <= ids.size());
    const auto offset = m_offsets[zi];
    for(axom::IndexType i = 0; i < sz; i++)
    {
      const auto idx = m_indices[offset + i];

      ids[i] = m_material_ids[idx];
      vfs[i] = m_volume_fractions[idx];
    }
    return sz;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat) const
  {
    FloatType tmp {};
    return zoneContainsMaterial(zi, mat, tmp);
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat, FloatType &vf) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));
    const auto sz = numberOfMaterials(zi);
    const auto offset = m_offsets[zi];
    for(axom::IndexType i = 0; i < sz; i++)
    {
      const auto idx = m_indices[offset + i];

      if(m_material_ids[idx] == mat)
      {
        vf = m_volume_fractions[idx];
        return true;
      }
    }
    vf = 0;
    return false;
  }

  /*!
   * \brief An iterator class for iterating over read-only data in a zone.
   *        The iterator can access material ids and volume fractions for one
   *        material at a time in the associated zone.
   */
  class const_iterator
  {
    // Let the material view call the const_iterator constructor.
    friend class UnibufferMaterialView<IndexT, FloatT, MAXMATERIALS>;

  public:
    /// Get the current material id for the iterator.
    MaterialID AXOM_HOST_DEVICE material_id() const
    {
      SLIC_ASSERT(m_currentIndex < size());
      return m_view->m_material_ids[m_index];
    }
    /// Get the current volume fraction for the iterator.
    FloatType AXOM_HOST_DEVICE volume_fraction() const
    {
      SLIC_ASSERT(m_currentIndex < size());
      return m_view->m_volume_fractions[m_index];
    }
    axom::IndexType AXOM_HOST_DEVICE size() const { return m_view->m_sizes[m_zoneIndex]; }
    ZoneIndex AXOM_HOST_DEVICE zoneIndex() const { return m_zoneIndex; }

    void AXOM_HOST_DEVICE operator++() { advance(true); }
    void AXOM_HOST_DEVICE operator++(int) { advance(true); }
    bool AXOM_HOST_DEVICE operator==(const const_iterator &rhs) const
    {
      return m_currentIndex == rhs.m_currentIndex && m_zoneIndex == rhs.m_zoneIndex &&
        m_view == rhs.m_view;
    }
    bool AXOM_HOST_DEVICE operator!=(const const_iterator &rhs) const
    {
      return m_currentIndex != rhs.m_currentIndex || m_zoneIndex != rhs.m_zoneIndex ||
        m_view != rhs.m_view;
    }

  private:
    DISABLE_DEFAULT_CTOR(const_iterator);

    /// Constructor
    AXOM_HOST_DEVICE const_iterator(const UnibufferMaterialView<IndexT, FloatT, MAXMATERIALS> *view,
                                    ZoneIndex zoneIndex,
                                    axom::IndexType currentIndex = 0)
      : m_view(view)
      , m_zoneIndex(zoneIndex)
      , m_currentIndex(currentIndex)
      , m_index(0)
    { }

    void AXOM_HOST_DEVICE advance(bool doIncrement)
    {
      m_currentIndex += (doIncrement && m_currentIndex < size()) ? 1 : 0;
      const auto idx = m_view->m_offsets[m_zoneIndex] + m_currentIndex;
      if(idx < m_view->m_indices.size())
      {
        m_index = m_view->m_indices[idx];
      }
    }

    const UnibufferMaterialView<IndexT, FloatT, MAXMATERIALS> *m_view;
    ZoneIndex m_zoneIndex;
    axom::IndexType m_currentIndex;
    axom::IndexType m_index;  // not considered in ==, !=
  };
  // Let the const_iterator access members.
  friend class const_iterator;

  /*!
   * \brief Return the iterator for the beginning of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the beginning of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE beginZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    auto it = const_iterator(this, zi, 0);
    it.advance(false);
    return it;
  }

  /*!
   * \brief Return the iterator for the end of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the end of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE endZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    return const_iterator(this, zi, m_sizes[zi]);
  }

private:
  axom::ArrayView<MaterialID> m_material_ids;
  axom::ArrayView<FloatType> m_volume_fractions;
  axom::ArrayView<IndexType> m_sizes;
  axom::ArrayView<IndexType> m_offsets;
  axom::ArrayView<IndexType> m_indices;
};

/*!
 \brief View for multi-buffer matsets.

 \tparam IndexT The integer type used for material data.
 \tparam FloatT The floating point type used for material data (volume fractions).
 \tparam MAXMATERIALS The maximum number of materials to support.

 \verbatim

matsets:
  matset:
    topology: topology
    volume_fractions:
      a:
        values: [0, 0, 0, a1, 0, a0]
        indices: [5, 3]
      b:
        values: [0, b0, b2, b1, 0]
        indices: [1, 3, 2]
    material_map: # (optional)
      a: 0
      b: 1

 \endverbatim
 */
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS = 20>
class MultiBufferMaterialView
{
public:
  using MaterialID = IndexT;
  using ZoneIndex = IndexT;
  using IndexType = IndexT;
  using FloatType = FloatT;
  using IDList = StaticArray<IndexType, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static axom::IndexType MaxMaterials = MAXMATERIALS;
  constexpr static axom::IndexType InvalidIndex = -1;

  void add(MaterialID matno,
           const axom::ArrayView<ZoneIndex> &indices,
           const axom::ArrayView<FloatType> &vfs)
  {
    SLIC_ASSERT(m_size + 1 < MaxMaterials);
#if !defined(AXOM_DEVICE_CODE)
    const auto begin = m_matnos.data();
    const auto end = begin + m_size;
    SLIC_ERROR_IF(std::find(begin, end, matno) != end,
                  "Adding a duplicate material number is not allowed.");
#endif
    m_indices[m_size] = indices;
    m_values[m_size] = vfs;
    m_matnos[m_size] = matno;
    m_size++;
  }

  AXOM_HOST_DEVICE
  axom::IndexType numberOfZones() const
  {
    axom::IndexType nzones = 0;
    for(int i = 0; i < m_size; i++) nzones = axom::utilities::max(nzones, m_indices[i].size());
    return nzones;
  }

  AXOM_HOST_DEVICE
  axom::IndexType numberOfMaterials(ZoneIndex zi) const
  {
    axom::IndexType nmats = 0;
    for(axom::IndexType i = 0; i < m_size; i++)
    {
      const auto &curIndices = m_indices[i];
      const auto &curValues = m_values[i];

      if(zi < static_cast<ZoneIndex>(curIndices.size()))
      {
        const auto idx = curIndices[zi];
        nmats += (curValues[idx] > 0) ? 1 : 0;
      }
    }

    return nmats;
  }

  AXOM_HOST_DEVICE
  void zoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    ids.clear();
    vfs.clear();

    for(axom::IndexType i = 0; i < m_size; i++)
    {
      const auto &curIndices = m_indices[i];
      const auto &curValues = m_values[i];

      if(zi < static_cast<ZoneIndex>(curIndices.size()))
      {
        const auto idx = curIndices[zi];
        if(curValues[idx] > 0)
        {
          ids.push_back(m_matnos[i]);
          vfs.push_back(curValues[idx]);
        }
      }
    }
  }

  AXOM_HOST_DEVICE
  axom::IndexType zoneMaterials(ZoneIndex zi,
                                axom::ArrayView<IndexType> &ids,
                                axom::ArrayView<FloatType> &vfs) const
  {
    axom::IndexType n = 0;
    for(axom::IndexType i = 0; i < m_size; i++)
    {
      const auto &curIndices = m_indices[i];
      const auto &curValues = m_values[i];

      if(zi < static_cast<ZoneIndex>(curIndices.size()))
      {
        const auto idx = curIndices[zi];
        if(curValues[idx] > 0)
        {
          ids[n] = m_matnos[i];
          vfs[n] = curValues[idx];
          n++;
        }
      }
    }
    return n;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat) const
  {
    FloatType tmp {};
    return zoneContainsMaterial(zi, mat, tmp);
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat, FloatType &vf) const
  {
    bool found = false;
    vf = FloatType {};
    axom::IndexType mi = indexOfMaterialID(mat);
    if(mi != InvalidIndex)
    {
      const auto &curIndices = m_indices[mi];
      const auto &curValues = m_values[mi];
      if(zi < static_cast<ZoneIndex>(curIndices.size()))
      {
        const auto idx = curIndices[zi];
        vf = curValues[idx];
        found = curValues[idx] > 0;
      }
    }
    return found;
  }

  /*!
   * \brief An iterator class for iterating over read-only data in a zone.
   *        The iterator can access material ids and volume fractions for one
   *        material at a time in the associated zone.
   */
  class const_iterator
  {
    // Let the material view call the const_iterator constructor.
    friend class MultiBufferMaterialView<IndexT, FloatT, MAXMATERIALS>;

  public:
    /// Get the current material id for the iterator.
    MaterialID AXOM_HOST_DEVICE material_id() const
    {
      SLIC_ASSERT(m_currentIndex < m_view->m_size);
      return m_view->m_matnos[m_currentIndex];
    }
    /// Get the current volume fraction for the iterator.
    FloatType AXOM_HOST_DEVICE volume_fraction() const
    {
      SLIC_ASSERT(m_currentIndex < m_view->m_size);
      const auto &curIndices = m_view->m_indices[m_currentIndex];
      const auto &curValues = m_view->m_values[m_currentIndex];
      const auto idx = curIndices[m_zoneIndex];
      return curValues[idx];
    }
    axom::IndexType AXOM_HOST_DEVICE size() const { return m_view->numberOfMaterials(m_zoneIndex); }
    ZoneIndex AXOM_HOST_DEVICE zoneIndex() const { return m_zoneIndex; }
    void AXOM_HOST_DEVICE operator++()
    {
      m_currentIndex += (m_currentIndex < m_view->m_size) ? 1 : 0;
      advance();
    }
    void AXOM_HOST_DEVICE operator++(int)
    {
      m_currentIndex += (m_currentIndex < m_view->m_size) ? 1 : 0;
      advance();
    }
    bool AXOM_HOST_DEVICE operator==(const const_iterator &rhs) const
    {
      return m_currentIndex == rhs.m_currentIndex && m_zoneIndex == rhs.m_zoneIndex &&
        m_view == rhs.m_view;
    }
    bool AXOM_HOST_DEVICE operator!=(const const_iterator &rhs) const
    {
      return m_currentIndex != rhs.m_currentIndex || m_zoneIndex != rhs.m_zoneIndex ||
        m_view != rhs.m_view;
    }

  private:
    DISABLE_DEFAULT_CTOR(const_iterator);

    /// Constructor
    AXOM_HOST_DEVICE const_iterator(const MultiBufferMaterialView<IndexT, FloatT, MAXMATERIALS> *view,
                                    ZoneIndex zoneIndex,
                                    axom::IndexType currentIndex = 0)
      : m_view(view)
      , m_zoneIndex(zoneIndex)
      , m_currentIndex(currentIndex)
    { }

    /// Advance to the next valid material slot for the zone.
    void AXOM_HOST_DEVICE advance()
    {
      while(m_currentIndex < m_view->m_size)
      {
        const auto &curIndices = m_view->m_indices[m_currentIndex];
        const auto &curValues = m_view->m_values[m_currentIndex];

        if(m_zoneIndex < static_cast<ZoneIndex>(curIndices.size()))
        {
          const auto idx = curIndices[m_zoneIndex];
          if(curValues[idx] > 0)
          {
            break;
          }
        }
        m_currentIndex++;
      }
    }

    const MultiBufferMaterialView<IndexT, FloatT, MAXMATERIALS> *m_view;
    ZoneIndex m_zoneIndex;
    axom::IndexType m_currentIndex;
  };
  // Let the const_iterator access members.
  friend class const_iterator;

  /*!
   * \brief Return the iterator for the beginning of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the beginning of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE beginZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    auto it = const_iterator(this, zi, 0);
    it.advance();
    return it;
  }

  /*!
   * \brief Return the iterator for the end of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the end of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE endZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));
    return const_iterator(this, zi, m_size);
  }

private:
  AXOM_HOST_DEVICE
  axom::IndexType indexOfMaterialID(MaterialID mat) const
  {
    axom::IndexType index = InvalidIndex;
    for(axom::IndexType mi = 0; mi < m_size; mi++)
    {
      if(mat == m_matnos[mi])
      {
        index = mi;
        break;
      }
    }
    return index;
  }

  axom::StackArray<axom::ArrayView<FloatType>, MAXMATERIALS> m_values {};
  axom::StackArray<axom::ArrayView<ZoneIndex>, MAXMATERIALS> m_indices {};
  axom::StackArray<MaterialID, MAXMATERIALS> m_matnos {};
  axom::IndexType m_size {0};
};

/*!
 \brief View for element-dominant matsets.

 \tparam IndexT The integer type used for material data.
 \tparam FloatT The floating point type used for material data (volume fractions).
 \tparam MAXMATERIALS The maximum number of materials to support.

 \verbatim

matsets:
  matset:
    topology: topology
    volume_fractions:
      a: [a0, a1, 0]
      b: [b0, b1, b2]
      c: [0, 0, c2]
    material_map: # (optional)
      a: 0
      b: 1
      c: 2

 \endverbatim
 */
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS = 20>
class ElementDominantMaterialView
{
public:
  using MaterialID = IndexT;
  using ZoneIndex = IndexT;
  using IndexType = IndexT;
  using FloatType = FloatT;
  using IDList = StaticArray<IndexType, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static axom::IndexType MaxMaterials = MAXMATERIALS;
  constexpr static axom::IndexType InvalidIndex = -1;

  void add(MaterialID matno, const axom::ArrayView<FloatType> &vfs)
  {
#if !defined(AXOM_DEVICE_CODE)
    const auto begin = m_matnos.data();
    const auto end = begin + m_volume_fractions.size();
    SLIC_ERROR_IF(std::find(begin, end, matno) != end,
                  "Adding a duplicate material number is not allowed.");
#endif
    if((m_volume_fractions.size() + 1) < m_volume_fractions.capacity())
    {
      m_matnos[m_volume_fractions.size()] = matno;
      m_volume_fractions.push_back(vfs);
    }
#if !defined(AXOM_DEVICE_CODE)
    else
    {
      SLIC_ERROR("Attempted to add more than the maximum number of materials.");
    }
    checkSizes();
#endif
  }

  AXOM_HOST_DEVICE
  axom::IndexType numberOfZones() const
  {
    return (m_volume_fractions.size() > 0) ? m_volume_fractions[0].size() : 0;
  }

  AXOM_HOST_DEVICE
  axom::IndexType numberOfMaterials(ZoneIndex zi) const
  {
    axom::IndexType nmats = 0;
    for(axom::IndexType i = 0; i < m_volume_fractions.size(); i++)
    {
      const auto &currentVF = m_volume_fractions[i];
      SLIC_ASSERT(zi < currentVF.size());
      nmats += currentVF[zi] > 0 ? 1 : 0;
    }
    return nmats;
  }

  AXOM_HOST_DEVICE
  void zoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    ids.clear();
    vfs.clear();

    for(axom::IndexType i = 0; i < m_volume_fractions.size(); i++)
    {
      const auto &currentVF = m_volume_fractions[i];
      SLIC_ASSERT(zi < currentVF.size());
      if(currentVF[zi] > 0)
      {
        ids.push_back(m_matnos[i]);
        vfs.push_back(currentVF[zi]);
      }
    }
  }

  AXOM_HOST_DEVICE
  axom::IndexType zoneMaterials(ZoneIndex zi,
                                axom::ArrayView<IndexType> &ids,
                                axom::ArrayView<FloatType> &vfs) const
  {
    axom::IndexType n = 0;
    for(axom::IndexType i = 0; i < m_volume_fractions.size(); i++)
    {
      const auto &currentVF = m_volume_fractions[i];
      SLIC_ASSERT(zi < currentVF.size());
      if(currentVF[zi] > 0)
      {
        ids[n] = m_matnos[i];
        vfs[n] = currentVF[zi];
        n++;
      }
    }
    return n;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat) const
  {
    FloatType tmp {};
    return zoneContainsMaterial(zi, mat, tmp);
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat, FloatType &vf) const
  {
    bool found = false;
    vf = FloatType {};
    int mi = indexOfMaterialID(mat);
    if(mi != InvalidIndex)
    {
      const auto &currentVF = m_volume_fractions[mi];
      SLIC_ASSERT(zi < currentVF.size());
      vf = currentVF[zi];
      found = vf > 0;
    }
    return found;
  }

  /*!
   * \brief An iterator class for iterating over read-only data in a zone.
   *        The iterator can access material ids and volume fractions for one
   *        material at a time in the associated zone.
   */
  class const_iterator
  {
    // Let the material view call the const_iterator constructor.
    friend class ElementDominantMaterialView<IndexT, FloatT, MAXMATERIALS>;

  public:
    /// Get the current material id for the iterator.
    MaterialID AXOM_HOST_DEVICE material_id() const
    {
      SLIC_ASSERT(m_currentIndex < m_view->m_volume_fractions.size());
      return m_view->m_matnos[m_currentIndex];
    }
    /// Get the current volume fraction for the iterator.
    FloatType AXOM_HOST_DEVICE volume_fraction() const
    {
      SLIC_ASSERT(m_currentIndex < m_view->m_volume_fractions.size());
      return m_view->m_volume_fractions[m_currentIndex][m_zoneIndex];
    }
    axom::IndexType AXOM_HOST_DEVICE size() const { return m_view->numberOfMaterials(m_zoneIndex); }
    ZoneIndex AXOM_HOST_DEVICE zoneIndex() const { return m_zoneIndex; }
    void AXOM_HOST_DEVICE operator++()
    {
      m_currentIndex += (m_currentIndex < m_view->m_volume_fractions.size()) ? 1 : 0;
      advance();
    }
    void AXOM_HOST_DEVICE operator++(int)
    {
      m_currentIndex += (m_currentIndex < m_view->m_volume_fractions.size()) ? 1 : 0;
      advance();
    }
    bool AXOM_HOST_DEVICE operator==(const const_iterator &rhs) const
    {
      return m_currentIndex == rhs.m_currentIndex && m_zoneIndex == rhs.m_zoneIndex &&
        m_view == rhs.m_view;
    }
    bool AXOM_HOST_DEVICE operator!=(const const_iterator &rhs) const
    {
      return m_currentIndex != rhs.m_currentIndex || m_zoneIndex != rhs.m_zoneIndex ||
        m_view != rhs.m_view;
    }

  private:
    DISABLE_DEFAULT_CTOR(const_iterator);

    /// Constructor
    AXOM_HOST_DEVICE const_iterator(const ElementDominantMaterialView<IndexT, FloatT, MAXMATERIALS> *view,
                                    ZoneIndex zoneIndex,
                                    axom::IndexType currentIndex = 0)
      : m_view(view)
      , m_zoneIndex(zoneIndex)
      , m_currentIndex(currentIndex)
    { }

    /// Advance to the next valid material slot for the zone.
    void AXOM_HOST_DEVICE advance()
    {
      while(m_currentIndex < m_view->m_volume_fractions.size())
      {
        if(m_view->m_volume_fractions[m_currentIndex][m_zoneIndex] > 0)
        {
          break;
        }
        m_currentIndex++;
      }
    }

    const ElementDominantMaterialView<IndexT, FloatT, MAXMATERIALS> *m_view;
    ZoneIndex m_zoneIndex;
    axom::IndexType m_currentIndex;
  };
  // Let the const_iterator access members.
  friend class const_iterator;

  /*!
   * \brief Return the iterator for the beginning of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the beginning of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE beginZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    auto it = const_iterator(this, zi, 0);
    it.advance();
    return it;
  }

  /*!
   * \brief Return the iterator for the end of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the end of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE endZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));
    return const_iterator(this, zi, m_volume_fractions.size());
  }

private:
#if !defined(AXOM_DEVICE_CODE)
  void checkSizes() const
  {
    axom::IndexType size = m_volume_fractions.empty() ? 0 : m_volume_fractions[0].size();
    for(axom::IndexType i = 1; i < m_volume_fractions.size(); i++)
    {
      SLIC_ERROR_IF(m_volume_fractions[i].size() != size,
                    "Volume fraction arrays have different sizes.");
    }
  }
#endif

  AXOM_HOST_DEVICE
  axom::IndexType indexOfMaterialID(MaterialID mat) const
  {
    axom::IndexType index = InvalidIndex;
    const auto size = m_volume_fractions.size();
    for(axom::IndexType mi = 0; mi < size; mi++)
    {
      if(mat == m_matnos[mi])
      {
        index = mi;
        break;
      }
    }
    return index;
  }

  axom::StaticArray<axom::ArrayView<FloatType>, MAXMATERIALS> m_volume_fractions {};
  axom::StackArray<MaterialID, MAXMATERIALS> m_matnos {};
};

/*!
 \brief View for material-dominant matsets.

 \tparam IndexT The integer type used for material data.
 \tparam FloatT The floating point type used for material data (volume fractions).
 \tparam MAXMATERIALS The maximum number of materials to support.

 \verbatim

matsets:
  matset:
    topology: topology
    volume_fractions:
      a: [a0, a1]
      b: [b0, b1, b2]
      c: [c2]
    element_ids:
      a: [0, 1]
      b: [0, 1, 2]
      c: [2]
    material_map: # (optional)
      a: 0
      b: 1
      c: 2

 \endverbatim

 \note This matset type does not seem so GPU friendly since there is some work to do for some of the queries.

 */
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS = 20>
class MaterialDominantMaterialView
{
public:
  using MaterialID = IndexT;
  using ZoneIndex = IndexT;
  using IndexType = IndexT;
  using FloatType = FloatT;
  using IDList = StaticArray<IndexType, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static axom::IndexType MaxMaterials = MAXMATERIALS;
  constexpr static axom::IndexType InvalidIndex = -1;

  void add(MaterialID matno,
           const axom::ArrayView<ZoneIndex> &ids,
           const axom::ArrayView<FloatType> &vfs)
  {
#if !defined(AXOM_DEVICE_CODE)
    SLIC_ERROR_IF(ids.size() != vfs.size(), "Array views for ids, vfs have different sizes.");
    SLIC_ERROR_IF(m_size + 1 >= MaxMaterials,
                  "Attempted to add more than the maximum number of materials.");
    const auto begin = m_matnos.data();
    const auto end = begin + m_size;
    SLIC_ERROR_IF(std::find(begin, end, matno) != end,
                  "Adding a duplicate material number is not allowed.");
#endif
    m_element_ids[m_size] = ids;
    m_volume_fractions[m_size] = vfs;
    m_matnos[m_size] = matno;
    m_size++;
  }

  AXOM_HOST_DEVICE
  axom::IndexType numberOfZones() const
  {
    axom::IndexType nzones = -1;
    for(axom::IndexType mi = 0; mi < m_size; mi++)
    {
      const auto &element_ids = m_element_ids[mi];
      const auto sz = element_ids.size();
      for(axom::IndexType i = 0; i < sz; i++)
      {
        const auto ei = static_cast<axom::IndexType>(element_ids[i]);
        nzones = axom::utilities::max(nzones, ei);
      }
    }
    nzones++;
    return nzones;
  }

  AXOM_HOST_DEVICE
  axom::IndexType numberOfMaterials(ZoneIndex zi) const
  {
    axom::IndexType nmats = 0;
    for(axom::IndexType mi = 0; mi < m_size; mi++)
    {
      const auto &element_ids = m_element_ids[mi];
      const auto sz = element_ids.size();
      for(axom::IndexType i = 0; i < sz; i++)
      {
        if(element_ids[i] == zi)
        {
          nmats++;
          break;
        }
      }
    }
    return nmats;
  }

  AXOM_HOST_DEVICE
  void zoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    ids.clear();
    vfs.clear();

    for(axom::IndexType mi = 0; mi < m_size; mi++)
    {
      const auto &element_ids = m_element_ids[mi];
      const auto &volume_fractions = m_volume_fractions[mi];
      const auto sz = element_ids.size();
      for(axom::IndexType i = 0; i < sz; i++)
      {
        if(element_ids[i] == zi)
        {
          ids.push_back(m_matnos[mi]);
          vfs.push_back(volume_fractions[i]);
          break;
        }
      }
    }
  }

  AXOM_HOST_DEVICE
  axom::IndexType zoneMaterials(ZoneIndex zi,
                                axom::ArrayView<IndexType> &ids,
                                axom::ArrayView<FloatType> &vfs) const
  {
    axom::IndexType n = 0;
    for(axom::IndexType mi = 0; mi < m_size; mi++)
    {
      const auto &element_ids = m_element_ids[mi];
      const auto &volume_fractions = m_volume_fractions[mi];
      const auto sz = element_ids.size();
      for(axom::IndexType i = 0; i < sz; i++)
      {
        if(element_ids[i] == zi)
        {
          ids[n] = m_matnos[mi];
          vfs[n] = volume_fractions[i];
          n++;
          break;
        }
      }
    }
    return n;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat) const
  {
    FloatType tmp {};
    return zoneContainsMaterial(zi, mat, tmp);
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialID mat, FloatType &vf) const
  {
    bool found = false;
    vf = FloatType {};
    axom::IndexType mi = indexOfMaterialID(mat);
    if(mi != InvalidIndex)
    {
      const auto &element_ids = m_element_ids[mi];
      const auto &volume_fractions = m_volume_fractions[mi];
      const auto n = element_ids.size();
      for(axom::IndexType i = 0; i < n; i++)
      {
        if(element_ids[i] == zi)
        {
          found = true;
          vf = volume_fractions[i];
          break;
        }
      }
    }
    return found;
  }

  /*!
   * \brief An iterator class for iterating over read-only data in a zone.
   *        The iterator can access material ids and volume fractions for one
   *        material at a time in the associated zone.
   */
  class const_iterator
  {
    // Let the material view call the const_iterator constructor.
    friend class MaterialDominantMaterialView<IndexT, FloatT, MAXMATERIALS>;

  public:
    /// Get the current material id for the iterator.
    MaterialID AXOM_HOST_DEVICE material_id() const
    {
      SLIC_ASSERT(m_miIndex < m_view->m_size);
      return m_view->m_matnos[m_miIndex];
    }
    /// Get the current volume fraction for the iterator.
    FloatType AXOM_HOST_DEVICE volume_fraction() const
    {
      SLIC_ASSERT(m_miIndex < m_view->m_size && m_index < m_view->m_element_ids[m_miIndex].size());
      return m_view->m_volume_fractions[m_miIndex][m_index];
    }
    ZoneIndex AXOM_HOST_DEVICE zoneIndex() const { return m_zoneIndex; }
    axom::IndexType AXOM_HOST_DEVICE size() const { return m_view->numberOfMaterials(m_zoneIndex); }
    void AXOM_HOST_DEVICE operator++() { advance(true); }
    void AXOM_HOST_DEVICE operator++(int) { advance(true); }
    bool AXOM_HOST_DEVICE operator==(const const_iterator &rhs) const
    {
      return m_miIndex == rhs.m_miIndex && m_index == rhs.m_index &&
        m_zoneIndex == rhs.m_zoneIndex && m_view == rhs.m_view;
    }
    bool AXOM_HOST_DEVICE operator!=(const const_iterator &rhs) const
    {
      return m_miIndex != rhs.m_miIndex || m_index != rhs.m_index ||
        m_zoneIndex != rhs.m_zoneIndex || m_view != rhs.m_view;
    }

  private:
    DISABLE_DEFAULT_CTOR(const_iterator);

    /// Constructor
    AXOM_HOST_DEVICE const_iterator(const MaterialDominantMaterialView<IndexT, FloatT, MAXMATERIALS> *view,
                                    ZoneIndex zoneIndex,
                                    axom::IndexType miIndex,
                                    axom::IndexType index)
      : m_view(view)
      , m_zoneIndex(zoneIndex)
      , m_miIndex(miIndex)
      , m_index(index)
    { }

    /// Advance to the next valid material slot for the zone.
    void AXOM_HOST_DEVICE advance(bool doIncrement)
    {
      if(doIncrement)
      {
        if(m_miIndex < m_view->m_size)
        {
          m_index = 0;
          m_miIndex++;
        }
        if(m_miIndex == m_view->m_size)
        {
          m_index = m_view->m_element_ids[m_view->m_size - 1].size();
        }
      }

      // Look for the next m_miIndex,m_index pair that contains material for the selected zone index.
      for(; m_miIndex < m_view->m_size; m_miIndex++)
      {
        const auto &element_ids = m_view->m_element_ids[m_miIndex];
        const auto sz = element_ids.size();
        for(; m_index < sz; m_index++)
        {
          if(element_ids[m_index] == m_zoneIndex)
          {
            return;
          }
        }
        m_index = (m_miIndex + 1 == m_view->m_size) ? m_index : 0;
      }
    }

    const MaterialDominantMaterialView<IndexT, FloatT, MAXMATERIALS> *m_view;
    ZoneIndex m_zoneIndex;
    axom::IndexType m_miIndex;
    axom::IndexType m_index;
  };
  // Let the const_iterator access members.
  friend class const_iterator;

  /*!
   * \brief Return the iterator for the beginning of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the beginning of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE beginZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));

    auto it = const_iterator(this, zi, 0, 0);
    it.advance(false);
    return it;
  }

  /*!
   * \brief Return the iterator for the end of a zone's material data.
   *
   * \param zi The zone index being queried.
   *
   * \return The iterator for the end of a zone's material data.
   */
  const_iterator AXOM_HOST_DEVICE endZone(ZoneIndex zi) const
  {
    SLIC_ASSERT(zi < static_cast<ZoneIndex>(numberOfZones()));
    const axom::IndexType miIndex = m_size;
    const axom::IndexType index = (m_size > 0) ? (m_volume_fractions[m_size - 1].size()) : 0;
    return const_iterator(this, zi, miIndex, index);
  }

private:
  AXOM_HOST_DEVICE
  axom::IndexType indexOfMaterialID(MaterialID mat) const
  {
    axom::IndexType index = InvalidIndex;
    for(axom::IndexType mi = 0; mi < m_size; mi++)
    {
      if(mat == m_matnos[mi])
      {
        index = mi;
        break;
      }
    }
    return index;
  }

  axom::StackArray<axom::ArrayView<IndexType>, MAXMATERIALS> m_element_ids {};
  axom::StackArray<axom::ArrayView<FloatType>, MAXMATERIALS> m_volume_fractions {};
  axom::StackArray<MaterialID, MAXMATERIALS> m_matnos {};
  axom::IndexType m_size {0};
};

}  // end namespace views
}  // end namespace bump
}  // end namespace axom

#endif
