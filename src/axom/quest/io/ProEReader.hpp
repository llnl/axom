// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PROEREADER_HPP_
#define QUEST_PROEREADER_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"

// C/C++ includes
#include <string>  // for std::string
#include <vector>  // for std::vector

namespace axom
{
namespace quest
{
/*!
 * \class ProEReader
 *
 * \brief A simple reader for a Pro/E file encoded in ascii.
 *        Current support for only tetrahedra.
 *
 * Pro/Engineer (also known as Creo) is a modeling application.
 *
 * To read an ASCII Pro/E tet file,
 *   - instantiate a ProEReader,
 *   - set the file name,
 *   - optionally set a tetrahedron predicate to specify a subset of the mesh,
 *   - call \a read(),
 *   - optionally, query the number of nodes or tets,
 *   - retrieve the mesh with \a getMesh().
 *
 * The tetrahedron predicate is optional.  If none is specified using
 * either \a setTetPredFromBoundingBox or \a setTetPred, the reader
 * retains all tets.  In any case, the reader retains all nodes in the
 * Pro/E file, even those not referenced by any tets.
 *
 * The tetrahedron predicate is intended to save memory by retaining tets
 * for which the predicate returns true and discarding all others.  A code
 * can use the convenience function \a setTetPredFromBoundingBox() to
 * discard all tets outside a bounding box, or specify an arbitrary
 * function (\a TetPred) for more complicated decisions.
 *
 * \note Pro/E node IDs start at 1. ProEReader adjusts and stores
 *  node IDs to start at 0 for indexing.
 */
class ProEReader
{
public:
  constexpr static int NUM_NODES_PER_TET = 4;
  constexpr static int NUM_COMPS_PER_NODE = 3;
  using Point3D = primal::Point<double, NUM_COMPS_PER_NODE>;
  using BBox3D = primal::BoundingBox<double, 3>;
  /*! \brief Specify tets to keep.
    *
    * - First argument: Pro/E node IDs for the current tet (1-based)
    * - Second argument: Pro/E tet ID for the current tet (1-based)
    * - Third argument: The node locations, stored interleaved
    *   (x1, y1, z1, x2, y2, z2, ... )
    */
  using TetPred = std::function<bool(int[4], int, std::vector<double>&)>;

public:
  /*!
   * \brief Constructor.
   */
  ProEReader();

  /*!
   * \brief Destructor.
   */
  virtual ~ProEReader();

  /*!
   * \brief Sets the name of the file to read.
   * \param [in] fileName the name of the file to read.
   */
  void setFileName(const std::string& fileName) { m_fileName = fileName; };

  /*!
   * \brief Returns the number of nodes of the mesh.
   * \return numNodes the number of nodes.
   */
  int getNumNodes() const { return static_cast<int>(m_num_nodes); };

  /*!
   * \brief Returns the number of tetrahedra of the mesh
   * \return numTets the number of tetrahedra.
   */
  int getNumTets() const { return static_cast<int>(m_num_tets); };

  /*!
   * \brief Clears all internal data-structures
   */
  void clear();

  /*!
   * \brief Reads in the mesh from a Pro/E file
   * \pre path to input file has been set by calling `setFileName()`
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  virtual int read();

  /// \name Set tetrahedron predicate to read mesh subset
  /// @{
  /*!
   * The reader calls the \a TetPred \a p after reading each tet,
   * retaining tets where p returns true and discarding the rest.
   * If a TetPred is not set, the reader retains all tets.
   */
  virtual void setTetPred(const TetPred& p) { m_tetPredicate = p; }

  /*!
   * Convenience function to set a \a TetPred from a bounding box.
   * 
   * \param box The TetPred will be defined in terms of the provided
   *        bounding box.  If the box is invalid, no TetPred is constructed
   *        and all tets are retained.
   * \param inclusive If true (the default), the TetPred will keep all tets
   *        with at least one node in the box.  If false, the TetPred will
   *        discard any tet with at least one node outside the box.
   */
  virtual void setTetPredFromBoundingBox(BBox3D& box, bool inclusive = true);
  /// @}

  /*!
   * \brief Stores the Pro/E data in the supplied unstructured mesh object.
   * \param [in,out] mesh pointer to the unstructured mesh.
   * \pre mesh != nullptr.
   */
  void getMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh);

protected:
  std::string m_fileName;

  axom::IndexType m_num_nodes;
  axom::IndexType m_num_tets;

  std::vector<double> m_nodes;
  std::vector<int> m_tets;

  TetPred m_tetPredicate;

  /*!
   * \brief Compact internal mesh storage arrays
   * \param [in] elt_count Number of elements stored in array
   * 
   * This method compacts the vertex array and resizes the vertex and the
   * element arrays to their minimum sizes.
   */
  void compact_arrays(int elt_count);

private:
  DISABLE_COPY_AND_ASSIGNMENT(ProEReader);
  DISABLE_MOVE_AND_ASSIGNMENT(ProEReader);
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_PROEREADER_HPP_
