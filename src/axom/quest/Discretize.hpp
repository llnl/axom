// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISCRETIZE_HPP_
#define QUEST_DISCRETIZE_HPP_

// Axom includes
#include "axom/core/Macros.hpp"

// Geometry
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Octahedron.hpp"

namespace axom
{
// forward declare Mesh
namespace mint
{
class Mesh;
}  // end namespace mint

namespace quest
{
/// \name Discretize primitive shapes to linear shapes
/// @{

using SphereType = primal::Sphere<double, 3>;
using OctType = primal::Octahedron<double, 3>;
using Point2D = primal::Point<double, 2>;

/*!
 * \brief Given a primitive sphere and a refinement level, allocate and return
 *   a list of Octahedra approximating the shape.
 * \param [in] s The sphere to approximate
 * \param [in] levels The number of refinements to perform, in addition to
 *   a central level-zero octahedron
 * \param [out] out The newly-initialized Array of octahedra representing \a s
 * \param [out] octcount The number of elements in \a out
 * \return false for invalid input or error in computation; true otherwise
 *
 * This routine generates O(4^level) octahedra.  That's exponential growth.
 * Use appropriate caution.
 *
 * This routine initializes an Array, \a out.
 */
bool discretize(const SphereType& s, int levels, axom::Array<OctType>& out, int& octcount);

/*!
 * \brief Given a 2D polyline revolved around the positive X-axis, allocate
 *   and return a list of Octahedra approximating the shape.
 * \param [in] polyline The polyline to revolve around the X-axis
 * \param [in] len The number of points in \a polyline
 * \param [in] levels The number of refinements to perform, in addition to
 *   a central level-zero octahedron in each segment
 * \param [out] out The newly-initialized Array of octahedra representing the
 *   revolved polyline
 * \param [out] octcount The number of elements in \a out
 * \return false for invalid input or error in computation; true otherwise
 *
 * This routine generates n*O(2^level) octahedra, where n is the number of
 * segments in \a polyline (one less than the length).
 * That's exponential growth.  Use appropriate caution.
 *
 * This routine initializes an Array, \a out.
 */
template <typename ExecSpace>
bool discretize(const axom::ArrayView<Point2D>& polyline,
                int len,
                int levels,
                axom::Array<OctType>& out,
                int& octcount);

/// @}

/// \name Visualize octahedra as tet mesh
/// @{

/*!
 * \brief Produces a mesh of tets from an array of Octahedra
 *
 * \param [in] octs Array of Octahedron objects
 * \param [in] octcount Length of \a octs
 * \param [in] segcount Number of segments in the polyline originating the octs
 * \param [out] mesh Pointer to a mesh object where the mesh will be generated.
 *
 * \return status error code, zero on success
 *
 * This function creates a new Mint mesh out of an array of Octahedron objects.
 * The mesh will be valid for any oct array, but the fields it includes make sense
 * for the output of discretizing a revolved polyline:
 *  - octahedron_volume: octahedron volume calculated by summing tetrahedron volumes.
 *  - oct_as_polyhedron_volume: octahedron volume calculated by changing the oct
 *    into a Polyhedron and reporting the Polyhedron's volume.  Any discrepancy
 *    between this field and octahedron_volume is a bug.
 *  - segment_index: the zero-based index of the segment to which the parent
 *    octahedron belongs.
 *  - octahedron_index: the zero-based index of the parent octahedron within its
 *    segment.
 *  - level_of_refinement: the level of refinement to which parent octahedron belongs.
 *
 * \note Ownership of the mesh object is passed to the caller. Consequently,
 *  the caller is responsible for properly deallocating the mesh object that
 *  the return mesh pointer points to.
 */
int mesh_from_discretized_polyline(const axom::ArrayView<OctType>& octs,
                                   int octcount,
                                   int segcount,
                                   mint::Mesh*& mesh);

/// @}

}  // end namespace quest
}  // end namespace axom

#include "detail/Discretize_detail.hpp"

#endif  // QUEST_DISCRETIZE_HPP_
