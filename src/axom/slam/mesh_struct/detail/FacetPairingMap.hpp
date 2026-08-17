// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/**
 * \file FacetPairingMap.hpp
 *
 * \brief Pairs the interior faces of a vertex star during IA connectivity repair.
 *
 * This is used by IAMesh::fixVertexNeighborhood() after a vertex \a apex is inserted
 * and its star -- the new elements incident to \a apex -- is retriangulated. 
 * The job is to recover the element-element adjacencies that are \a internal 
 * to that star. It is a transient pairing table, not a persistent lookup container.
 *
 * ## Why the key is a face of the link
 *
 * Recall that the \a star of \a apex is the collection of elements (triangles in 2D, tetrahedra in 3D)
 * incident to \a apex, and the \a link of \a apex is the boundary of the star:
 * a closed polyline in 2D (or open if \a apex is on the domain boundary)
 * and a triangulated sphere in 3D (or a disk if \a apex is on the domain boundary).
 *
 * Each star element splits its faces into two kinds. 
 * In the IA convention, local face \a i is opposite local vertex \a i-1, 
 * and \a apex is stored at the last local vertex position, 
 * so the face opposite \a apex is local face 0:
 *   - Local face 0 does not contain \a apex. This is the element's contribution 
 *     to the LINK of \a apex (an edge in 2D, a triangle in 3D). Its neighbor lies 
 *     outside the star, so it is repaired directly by fixVertexNeighborhood() and never enters this pairer.
 *   - The other \a TDIM faces are all incident in \a apex and are paired here.
 *
 * Because every incident face contains \a apex, the apex carries no
 * discriminating information. Removing it leaves a face of the link (our hash Key):
 *   - 2D: an incident face is the edge {apex, a}; the key is the single link vertex.
 *   - 3D: an incident face is the triangle {apex, a, b}; the key is the link
 *         edge {a, b}, stored sorted so that both incident tetrahedra produce the same key.
 *
 * Two star elements are adjacent across an incident face if and only if they
 * share that link face, so matching keys is equivalent to reconstructing the
 * star's internal adjacency. The apex is the implicit extra vertex of every key.
 *
 * ## Algorithm: linear-probing open addressing with generational marking
 *
 * **Generational Marking:** Each entry stores a generation counter. Entries with
 * generation != current_generation are treated as empty without actually clearing
 * them. This makes prepareForInsertions() an O(1) "clear" (bump the generation)
 * rather than an O(table) reset between uses.
 *
 * **Generation Wraparound:** On 32-bit wraparound (every ~4 billion
 * prepareForInsertions() calls), we reset all entry generations to 0 and restart
 * at 1. This is the only O(table) operation in the table's lifetime.
 *
 * **Tombstone Markers:** Deleted entries are marked with a special TOMBSTONE value
 * to preserve probe chains. Tombstones can be reused for new insertions.
 *
 * **Thread-Local Storage:** Uses static thread_local storage so the backing array
 * is reused across calls. The table only ever grows: prepareForInsertions() grows
 * it to fit the largest star seen, so after a large pairing pass the table stays
 * oversized for the lifetime of the thread.
 *
 * ## Performance characteristics
 *
 * - **Lookup / insert:** O(1) expected.
 * - **Space:** cold, the table is sized to ~(2*TDIM)x the facet count, i.e. about
 *   12-37% load for typical star sizes. After a large pairing pass it remains at
 *   that high-water mark, so the effective load factor for subsequent small stars is much lower.
 * - **Allocations:** zero after the high-water mark is reached (thread-local reuse).
 *
 * \note This is NOT a general-purpose facet map. The reduced key (a face of the link)
 *  is valid because every key within a single pairing pass shares the same apex.
 *  The global manifold check, IAMesh::isConforming(), has no fixed apex to factor out
 *  and instead keys on the full sorted vertex tuple.
 *
 * \note Host-only by construction (thread_local + std::vector + std::optional).
 *  This is a serial, construction-time helper, not a portable (e.g. device-capable) container.
 *
 * ## Usage pattern
 *
 * ```cpp
 * FacetPairingMap<TDIM> map;
 * map.prepareForInsertions(expected_facet_count);
 *
 * for (each incident face of each star element):
 *   if (auto match = map.findAndExtract(link_face_key))
 *     // Second sighting: wire the two star elements adjacent across this link face
 *     updateAdjacencies(current_slot, match.value());
 *   else
 *     // First sighting: park which star element/face owns this link face
 *     map.insert(link_face_key, current_slot);
 *
 * assert(map.allFacetsPaired());  // every link face matched exactly twice
 * ```
 */

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace axom
{
namespace slam
{
namespace detail
{

/**
 * \brief A face of the link of the inserted vertex (the apex)
 *
 *   - In 2D (triangle star): \a v0 is the single link vertex
 *   - In 3D (tetrahedron star): (v0, v1) are the two endpoints of a link edge, stored in sorted order.
 *
 * The sorted ordering in 3D ensures that the same link edge produces an identical key
 * for both incident tetrahedra sharing in the star of \a apex.
 */
template <int TDIM, typename IndexType>
struct LinkFace
{
  static constexpr IndexType INVALID = static_cast<IndexType>(-1);

  IndexType v0 {INVALID};  ///< Link vertex (2D) or smaller link-edge endpoint (3D)
  IndexType v1 {INVALID};  ///< Unused in 2D; larger link-edge endpoint (3D)

  /// Default constructor - creates invalid key
  LinkFace() = default;

  /// Construct from key components (automatically sorted in 3D)
  /// \param a First link vertex / link-edge endpoint
  /// \param b Second link-edge endpoint (unused in 2D, sorted in 3D)
  LinkFace(IndexType a, IndexType b = INVALID) : v0(a), v1(b)
  {
    // In 3D, ensure consistent ordering: v0 <= v1
    // This allows the same link edge to match from either incident tetrahedron.
    if constexpr(TDIM == 3)
    {
      if(b < a)
      {
        axom::utilities::swap(v0, v1);
      }
    }
  }

  /// Equality comparison
  bool operator==(const LinkFace& other) const { return v0 == other.v0 && v1 == other.v1; }

  /// Inequality comparison
  bool operator!=(const LinkFace& other) const { return !(*this == other); }
};

/**
 * \brief Identifies which star element generated a given link face, and where.
 *
 * On the first sighting of a link face, fixVertexNeighborhood() stores the star
 * element that owns it together with the local face slot it occupies in that element.
 * On the second sighting this is returned so the caller can wire the two
 * star elements that are adjacent across the link face.
 */
template <typename IndexType>
struct StarFacetSlot
{
  static constexpr IndexType INVALID = static_cast<IndexType>(-1);

  IndexType element_idx {INVALID};  ///< Star element owning this link face
  IndexType local_face {INVALID};   ///< Local face slot in that element (0..VERTS_PER_ELEM-1)

  /// Default constructor
  StarFacetSlot() = default;

  /// Construct from element and local face slot
  StarFacetSlot(IndexType elem, IndexType face) : element_idx(elem), local_face(face) { }
};

/**
 * \brief Thread-local hash table for matching facets during mesh connectivity repair.
 *
 * This class implements a specialized hash table optimized for the facet-pairing
 * problem in mesh connectivity repair. It provides:
 * - O(1) expected lookup and insert
 * - Zero allocations after first use (thread-local storage)
 * - Generational marking to avoid clearing
 * - Tombstone deletion to maintain probe chains
 *
 * \tparam TDIM Topological dimension (2 for triangles, 3 for tetrahedra)
 * \tparam IndexType Index type for element and face indices (typically int)
 */
template <int TDIM, typename IndexType = int>
class FacetPairingMap
{
public:
  using KeyType = LinkFace<TDIM, IndexType>;
  using DataType = StarFacetSlot<IndexType>;

  static constexpr IndexType INVALID = static_cast<IndexType>(-1);
  // Emptiness is encoded by the per-entry generation counter (see Entry below),
  // not by a sentinel key, so there is no EMPTY_SLOT. TOMBSTONE_SLOT marks a
  // deleted entry within an otherwise live probe chain.
  static constexpr IndexType TOMBSTONE_SLOT = INVALID - 1;

  /// Golden ratio constant for Knuth's multiplicative hashing
  /// This constant (phi * 2^64) provides excellent bit distribution
  static constexpr std::uint64_t HASH_MULTIPLIER = 0x9e3779b97f4a7c15ULL;

private:
  /// Hash table entry combining key, data, and generation tracking
  struct Entry
  {
    KeyType key;                  ///< Facet key for matching
    DataType data;                ///< Associated element/face data
    unsigned int generation {0};  ///< Generation marker for fast clearing

    Entry() = default;
  };

  // Thread-local storage for zero-allocation operation across calls
  static thread_local std::vector<Entry> s_table;
  static thread_local unsigned int s_generation;

  std::size_t m_mask {0};    ///< Table size - 1 (for fast modulo via bitwise AND)
  int m_unpaired_count {0};  ///< Number of unpaired facets currently in table

public:
  /// Default constructor
  FacetPairingMap() = default;

  /**
   * \brief Prepare the hash table for an expected number of facet insertions.
   *
   * This method sizes the table appropriately and advances the generation counter.
   * It should be called before each batch of insertions.
   *
   * In 2D, each element contributes 2 incident faces (face 0 is on boundary).
   * In 3D, each element contributes 3 incident faces (face 0 is on boundary).
   *
   * \param expected_facet_count Estimated number of facets to insert
   *
   * \note Table size is rounded up to the next power of 2 for fast modulo via bitwise AND.
   */
  void prepareForInsertions(std::size_t expected_facet_count)
  {
    // Cold sizing targets a low load factor: ~(2*TDIM)x the facet count
    // (2D: ~2 incident faces/elem; 3D: ~3 incident faces/elem). The next
    // power-of-two round-up below lowers the realized load further.
    // The table never shrinks, so after a large pairing pass,
    // subsequent small stars run at a much lower load factor.
    const std::size_t target_slots =
      axom::utilities::max<std::size_t>(8, (TDIM == 2 ? 4 : 8) * expected_facet_count);

    // Round up to power of 2 for fast modulo via bitwise AND
    std::size_t table_size = 8;
    while(table_size < target_slots)
    {
      table_size <<= 1;
    }

    // Resize if needed (reuses existing allocation if possible)
    if(s_table.size() < table_size)
    {
      s_table.resize(table_size);
    }

    m_mask = s_table.size() - 1;
    m_unpaired_count = 0;

    // Advance generation counter (handles wraparound)
    advanceGeneration();
  }

  /**
   * \brief Try to find and remove a matching facet from the table.
   *
   * If a match is found, it's removed (marked as tombstone) and its data returned.
   * If no match exists, returns empty optional.
   *
   * \param key The facet key to search for
   * \return Optional containing matching facet data, or empty if not found
   *
   * \note This uses linear probing with generational marking. Empty slots
   *       (including stale generations) terminate the search.
   */
  std::optional<DataType> findAndExtract(const KeyType& key)
  {
    const std::size_t hash = computeHash(key);
    std::size_t slot = hash & m_mask;

    // Linear probe until we find the key or an empty slot
    while(true)
    {
      Entry& entry = s_table[slot];

      // Empty slot (including stale generation) - key not present
      if(entry.generation != s_generation)
      {
        return std::nullopt;
      }

      // Tombstone - continue searching (must maintain probe chain)
      if(entry.data.element_idx == TOMBSTONE_SLOT)
      {
        slot = (slot + 1) & m_mask;
        continue;
      }

      // Found matching key - remove it and return data
      if(entry.key == key)
      {
        DataType result = entry.data;
        entry.data.element_idx = TOMBSTONE_SLOT;  // Mark as deleted
        --m_unpaired_count;
        return result;
      }

      // Key mismatch - continue probing
      slot = (slot + 1) & m_mask;
    }
  }

  /**
   * \brief Insert a new facet into the table.
   *
   * Should only be called after findAndExtract() returns empty for this key.
   * Calling insert() for an already-present key is an error (assertion).
   *
   * \param key The facet key
   * \param data The facet data (element index and face index)
   *
   * \note This uses linear probing with tombstone reuse. If a tombstone is
   *       encountered during probing, it will be reused for insertion.
   */
  void insert(const KeyType& key, const DataType& data)
  {
    const std::size_t hash = computeHash(key);
    std::size_t slot = hash & m_mask;
    std::size_t first_tombstone = s_table.size();  // Track first tombstone for reuse

    while(true)
    {
      Entry& entry = s_table[slot];

      // Empty slot (including stale generation) - insert here (or at tombstone if found)
      if(entry.generation != s_generation)
      {
        Entry& target = (first_tombstone < s_table.size()) ? s_table[first_tombstone] : entry;
        target.key = key;
        target.data = data;
        target.generation = s_generation;
        ++m_unpaired_count;
        return;
      }

      // Tombstone - remember it for potential reuse
      if(entry.data.element_idx == TOMBSTONE_SLOT)
      {
        if(first_tombstone == s_table.size())
        {
          first_tombstone = slot;
        }
        slot = (slot + 1) & m_mask;
        continue;
      }

      // Collision with existing key - this is an error
      SLIC_ASSERT_MSG(entry.key != key,
                      "FacetPairingMap::insert() called with duplicate key. "
                      "Each facet should appear exactly once per element.");

      // Different key - continue probing
      slot = (slot + 1) & m_mask;
    }
  }

  /// \brief Returns the number of link faces seen an odd number of times so far
  /// \return Number of link faces still waiting for their second sighting
  int pendingCount() const { return m_unpaired_count; }

  /**
   * \brief Checks that every link face was paired exactly twice.
   *
   * Should be called after all insertions complete. For an interior apex this is
   * the manifold condition: the link is closed (a loop in 2D, a sphere in 3D),
   * so every link face is interior and is shared by exactly two star elements.
   * A non-zero count indicates a link face seen only once, i.e. a link with boundary
   * or a non-manifold/broken star.
   *
   * \return true if all link faces paired (unpaired count is zero)
   */
  bool allFacetsPaired() const { return m_unpaired_count == 0; }

  /// \brief Returns the current generation counter value (for testing)
  unsigned int currentGeneration() const { return s_generation; }

  /// \brief Returns the table size (for testing/diagnostics)
  std::size_t tableSize() const { return s_table.size(); }

  /// \brief Returns the table mask (for testing/diagnostics)
  std::size_t tableMask() const { return m_mask; }

private:
  /**
   * \brief Compute hash for a link-face key using Knuth's multiplicative method.
   *
   * Uses the golden ratio constant (phi * 2^64) which provides excellent
   * bit distribution. In 3D, the two link-edge endpoints are mixed together.
   *
   * \param key The link-face key to hash
   * \return Hash value
   *
   * \note In 2D the single link vertex is returned unmixed; the table relies on
   *  linear probing and power-of-two masking to spread sequential vertex ids.
   *  With spatially sorted (e.g. BRIO) insertion these ids are not random and can clump,
   *  so applying the same multiplicative mix in 2D may be worth measuring.
   */
  std::size_t computeHash(const KeyType& key) const
  {
    std::size_t seed = static_cast<std::size_t>(key.v0);

    if constexpr(TDIM == 3)
    {
      // Mix in the second link-edge endpoint using bit shifting and XOR.
      // This formula provides good avalanche properties.
      seed ^= static_cast<std::size_t>(key.v1) + HASH_MULTIPLIER + (seed << 6) + (seed >> 2);
    }

    return seed;
  }

  /**
   * \brief Advance the generation counter, handling wraparound.
   *
   * On wraparound (every ~4 billion calls), resets all entries to generation 0
   * and restarts at 1. This is the only O(n) operation in the table's lifetime.
   *
   * Wraparound is extremely rare in practice (would require billions of
   * prepareForInsertions() calls), but this ensures correctness.
   */
  void advanceGeneration()
  {
    if(++s_generation == 0u)
    {
      // Generation counter wrapped around - reset all entries
      // This is rare but ensures correctness
      for(auto& entry : s_table)
      {
        entry.generation = 0u;
      }
      s_generation = 1u;
    }
  }
};

// Static member initialization
template <int TDIM, typename IndexType>
thread_local std::vector<typename FacetPairingMap<TDIM, IndexType>::Entry>
  FacetPairingMap<TDIM, IndexType>::s_table;

template <int TDIM, typename IndexType>
thread_local unsigned int FacetPairingMap<TDIM, IndexType>::s_generation = 0;

}  // namespace detail
}  // namespace slam
}  // namespace axom
