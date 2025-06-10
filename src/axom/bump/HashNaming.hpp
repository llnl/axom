// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_BUMP_HASH_NAMING_HPP_
#define AXOM_BUMP_HASH_NAMING_HPP_

#include "axom/bump/utilities/utilities.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <conduit/conduit.hpp>
#include <conduit/conduit_blueprint.hpp>

#include <cstdint>

namespace axom
{
namespace bump
{

//------------------------------------------------------------------------------
/*!
 * \brief This class implements a naming policy that uses some hashing functions
 *        to produce a "name" for an array of ids.
 *
 * \tparam IndexT The index type that gets hashed.
 * \tparam MAXIDS The max number of ids that get hashed.
 */
template <typename IndexT, int MAXIDS = 14>
class HashNaming
{
public:
  using KeyType = std::uint64_t;
  using IndexType = IndexT;

  // The top 2 bits are reserved for the key type.
  constexpr static KeyType KeyIDSingle = 0;
  constexpr static KeyType KeyIDPair = KeyType(1) << 62;
  constexpr static KeyType KeyIDPack = KeyType(2) << 62;
  constexpr static KeyType KeyIDHash = KeyType(3) << 62;

  // The rest of the bits can be divided in various ways.
  constexpr static KeyType KeyMask = KeyType(3) << 62;
  constexpr static KeyType PayloadMask = ~KeyMask;
  constexpr static KeyType Max15Bit = (KeyType(1) << 15) - 1;
  constexpr static KeyType Max16Bit = (KeyType(1) << 16) - 1;
  constexpr static KeyType Max20Bit = (KeyType(1) << 20) - 1;
  constexpr static KeyType Max31Bit = (KeyType(1) << 31) - 1;
  constexpr static KeyType Max32Bit = (KeyType(1) << 32) - 1;

  /*!
   * \brief A view for making names, suitable for use in device code.
   */
  class View
  {
  public:
    using KeyType = HashNaming::KeyType;

    /*!
     * \brief Make a name from an array of ids.
     *
     * \param p The array of ids.
     * \param n The number of ids in the array.
     *
     * \return The name that describes the array of ids.
     *
     * \note Different make_name_* functions are used because we can skip most
     *       sorting for 1,2 element arrays. Also, there is a small benefit
     *       to some of the other shortcuts for smaller arrays.
     */
    AXOM_HOST_DEVICE
    KeyType makeName(const IndexType *p, int n) const
    {
      KeyType name {};
      if(n == 1)
        name = make_name_1(p[0]);
      else if(n == 2)
        name = make_name_2(p[0], p[1]);
      else
        name = make_name_n(p, n);
      return name;
    }

    /// Set the max number of nodes, which can be useful for packing/narrowing.
    AXOM_HOST_DEVICE
    void setMaxId(IndexType m) { m_maxId = static_cast<KeyType>(m); }

  private:
    /*!
     * \brief Encode a single id as a name.
     * \param p0 The id to encode.
     * \return A name that encodes the id.
     */
    AXOM_HOST_DEVICE
    inline KeyType make_name_1(IndexType p0) const
    {
      SLIC_ASSERT(static_cast<KeyType>(p0) < PayloadMask);
      // Store p0 in the key as a 62-bit integer
      KeyType k0 = (static_cast<KeyType>(p0) & PayloadMask);
      return KeyIDSingle | k0;
    }

    /*!
     * \brief Encode 2 ids as a name.
     * \param p0 The first id to encode.
     * \param p1 The second id to encode.
     * \return A name that encodes the ids.
     */
    AXOM_HOST_DEVICE
    inline KeyType make_name_2(IndexType p0, IndexType p1) const
    {
      SLIC_ASSERT(static_cast<KeyType>(p0) <= Max31Bit && static_cast<KeyType>(p1) <= Max31Bit);
      // Store p0 and p1 both in the 64-bit key as 31-bit integers
      KeyType k0 = (static_cast<KeyType>(axom::utilities::min(p0, p1)) & Max31Bit);
      KeyType k1 = (static_cast<KeyType>(axom::utilities::max(p0, p1)) & Max31Bit);
      return KeyIDPair | (k0 << 31) | k1;
    }

    /*!
     * \brief Encode multiple ids as a name.
     * \param p The ids to encode.
     * \param n The number of ids.
     * \return A name that encodes the ids.
     */
    AXOM_HOST_DEVICE
    KeyType make_name_n(const IndexType *p, int n) const
    {
      KeyType retval {};
      if(n == 3 && m_maxId <= Max20Bit)
      {
        // We can pack 3 values into the id lossless
        IndexType sorted[3];
        sorted[0] = p[0];
        sorted[1] = p[1];
        sorted[2] = p[2];
        axom::utilities::Sorting<IndexType, 3>::sort(sorted, n);

        KeyType k0 = static_cast<KeyType>(sorted[0]) & Max20Bit;
        KeyType k1 = static_cast<KeyType>(sorted[1]) & Max20Bit;
        KeyType k2 = static_cast<KeyType>(sorted[2]) & Max20Bit;
        constexpr KeyType len = KeyType(3 - 1) << 60;
        retval = KeyIDPack | len | (k0 << 40) | (k1 << 20) | k2;
      }
      else if(n == 4 && m_maxId <= Max15Bit)
      {
        // We can pack 4 values into the id lossless
        IndexType sorted[4];
        sorted[0] = p[0];
        sorted[1] = p[1];
        sorted[2] = p[2];
        sorted[3] = p[3];
        axom::utilities::Sorting<IndexType, 4>::sort(sorted, n);

        KeyType k0 = static_cast<KeyType>(sorted[0]) & Max15Bit;
        KeyType k1 = static_cast<KeyType>(sorted[1]) & Max15Bit;
        KeyType k2 = static_cast<KeyType>(sorted[2]) & Max15Bit;
        KeyType k3 = static_cast<KeyType>(sorted[3]) & Max15Bit;
        constexpr KeyType len = KeyType(4 - 1) << 60;
        retval = KeyIDPack | len | (k0 << 45) | (k1 << 30) | (k2 << 15) | k3;
      }
      else if(m_maxId < Max16Bit)
      {
        // Narrow to 16-bit, sort
        std::uint16_t sorted[MAXIDS];
        for(int i = 0; i < n; i++)
        {
          sorted[i] = static_cast<std::uint16_t>(p[i]);
        }
        axom::utilities::Sorting<std::uint16_t, MAXIDS>::sort(sorted, n);

        // Make a hash from the narrowed ids
        void *ptr = static_cast<void *>(sorted);
        KeyType k0 =
          axom::utilities::hash_bytes(static_cast<std::uint8_t *>(ptr), n * sizeof(std::uint16_t));
        retval = KeyIDHash | (k0 & PayloadMask);
      }
      else if(m_maxId < Max32Bit)
      {
        // Narrow to 32-bit, sort
        std::uint32_t sorted[MAXIDS];
        for(int i = 0; i < n; i++)
        {
          sorted[i] = static_cast<std::uint32_t>(p[i]);
        }
        axom::utilities::Sorting<std::uint32_t, MAXIDS>::sort(sorted, n);

        // Make a hash from the narrowed ids
        void *ptr = static_cast<void *>(sorted);
        KeyType k0 =
          axom::utilities::hash_bytes(static_cast<std::uint8_t *>(ptr), n * sizeof(std::uint32_t));
        retval = KeyIDHash | (k0 & PayloadMask);
      }
      else if(n > 0)
      {
        IndexType sorted[MAXIDS];
        for(int i = 0; i < n; i++)
        {
          sorted[i] = p[i];
        }
        axom::utilities::Sorting<IndexType, MAXIDS>::sort(sorted, n);

        // Make a hash from the ids
        void *ptr = static_cast<void *>(sorted);
        KeyType k0 =
          axom::utilities::hash_bytes(static_cast<std::uint8_t *>(ptr), n * sizeof(IndexType));
        retval = KeyIDHash | (k0 & PayloadMask);
      }
      return retval;
    }

    KeyType m_maxId {axom::numeric_limits<KeyType>::max()};
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
  View view() { return m_view; }

  /*!
   * \brief Turn name into a string.
   * \param key The name.
   * \return A string that represents the name.
   */
  static std::string toString(KeyType key)
  {
    std::stringstream ss;
    auto kt = key & KeyMask;
    if(kt == KeyIDSingle)
    {
      auto id = key & PayloadMask;
      ss << "single(" << std::hex << id << ")";
    }
    else if(kt == KeyIDPair)
    {
      auto payload = key & PayloadMask;
      auto p0 = (payload >> 31) & Max31Bit;
      auto p1 = payload & Max31Bit;
      ss << "pair(" << std::hex << p0 << ", " << p1 << ")";
    }
    else if(kt == KeyIDHash)
    {
      ss << "hash(" << std::hex << key << ")";
    }
    else if(kt == KeyIDPack)
    {
      auto npts = ((key >> 60) & 3) + 1;
      if(npts == 3)
      {
        auto p0 = (key >> 40) & Max20Bit;
        auto p1 = (key >> 20) & Max20Bit;
        auto p2 = (key)&Max20Bit;
        ss << "pack(" << std::hex << p0 << ", " << p1 << ", " << p2 << ")";
      }
      else if(npts == 4)
      {
        auto p0 = (key >> 45) & Max15Bit;
        auto p1 = (key >> 30) & Max15Bit;
        auto p2 = (key >> 15) & Max15Bit;
        auto p3 = (key)&Max15Bit;
        ss << "pack(" << std::hex << p0 << ", " << p1 << ", " << p2 << ", " << p3 << ")";
      }
    }
    return ss.str();
  }

  View m_view {};
};

}  // end namespace bump
}  // end namespace axom

#endif
