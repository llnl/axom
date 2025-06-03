// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_MATERIAL_VIEW_HPP_
#define AXOM_MIR_VIEWS_MATERIAL_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>

#include <string>
#include <map>
#include <vector>

namespace axom
{
namespace mir
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
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS>
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
    SLIC_ASSERT(material_ids.size() == volume_fractions.size());
    SLIC_ASSERT(sizes.size() == offsets.size());

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
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS>
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
           const axom::ArrayView<ZoneIndex> &ids,
           const axom::ArrayView<FloatType> &vfs)
  {
    SLIC_ASSERT(m_size + 1 < MaxMaterials);

    m_indices[m_size] = ids;
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
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS>
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
    if((m_volume_fractions.size() + 1) < m_volume_fractions.capacity())
    {
      m_matnos[m_volume_fractions.size()] = matno;
      m_volume_fractions.push_back(vfs);
    }
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

private:
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
template <typename IndexT, typename FloatT, axom::IndexType MAXMATERIALS>
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
    SLIC_ASSERT(m_size + 1 < MaxMaterials);

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
}  // end namespace mir
}  // end namespace axom

#endif
