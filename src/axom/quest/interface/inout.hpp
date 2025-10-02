// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_INOUT_INTERFACE_HPP_
#define QUEST_INOUT_INTERFACE_HPP_

// Axom includes
#include "axom/config.hpp"

// Quest includes
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"

// C/C++ includes
#include <string>
#include <memory>

/*!
 * \file inout.hpp
 *
 * \brief Defines a C-style interface for evaluating the \a inout query.
 *
 * Given a watertight surface mesh and an arbitrary point in space,
 * the \a inout query determines if the point is contained within the volume
 * enclosed by the surface mesh.
 *
 * The mesh can either be provided via a path to a mesh file, or as a pointer
 * to a \a mint::Mesh object. The interface currently supports reading
 * triangle meshes in the
 * <a href="https://en.wikipedia.org/wiki/STL_(file_format)">STL format</a>.
 *
 * The interface has several parameters that may be set before initializing the
 * query (via \a inout_init() ), and several functions that are available once
 * the query has been initialized.
 */

namespace axom
{
// Forward Mint declarations
namespace mint
{
class Mesh;
}

namespace quest
{
constexpr int QUEST_INOUT_SUCCESS = 0;
constexpr int QUEST_INOUT_FAILED = -1;

/// \name InOut query -- initialization and finalization functions
/// @{

/*!
 * \brief Initializes the quest inout query from a mesh file
 *
 * \param [in] file Path to an STL file containing the surface mesh
 * \param [in] comm The MPI communicator (when running in parallel)
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 *
 * \pre inout_initialized() == false
 * \post inout_initialized() == true, when rc is QUEST_INOUT_SUCCESS
 */
int inout_init(const std::string& file, MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Initialize the inout query using a pre-loaded mesh
 *
 * \param [inout] mesh Pointer to the input mesh. This pointer will
 * be updated during this invocation
 * \param [in] comm The MPI communicator (when running in parallel)
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 *
 * \pre inout_initialized() == false
 * \post inout_initialized() == true, when rc is QUEST_INOUT_SUCCESS
 * \warning The underlying data structure modifies the input mesh (e.g.
 * by welding vertices) and updates the \a mesh pointer. It is the user's
 * responsibility to update any other pointers to this same mesh.
 */
int inout_init(std::shared_ptr<mint::Mesh>& mesh, MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Finalizes the inout query
 *
 * \post inout_initialized() == false
 */
int inout_finalize();

/*!
 * \brief Predicate to test whether the inout query has been initialized
 *
 * \return True if the inout query has been initialized, false otherwise.
 */
bool inout_initialized();

/// @}

/// \name InOut query -- properties and querying functions
/// \note These must be called after \a inout_init()
/// @{

/*!
 * \brief Tests if the point (\a x, \a y, \a z) is inside the contained volume
 *
 * \param [in] x The x-coordinate of the query point
 * \param [in] y The y-coordinate of the query point
 * \param [in] z The z-coordinate of the query point
 * \return True if the point is within the contained volume or on the
 *    surface mesh, false otherwise.
 * \pre inout_initialized() == true
 */
bool inout_evaluate(double x, double y, double z = 0.);

/*!
 * \brief Tests an array of points for containment
 *
 * Upon successful completion, entries in array \a res
 * will have the value 1 for points that are inside
 * and value 0 otherwise.
 *
 * \param [in] x Array of x-coordinates for the query points
 * \param [in] y Array of y-coordinates for the query points
 * \param [in] z Array of z-coordinates for the query points
 * \param [in] npoints The number of points to test
 * \param [out] res An array of results. Each entry has value \a 1
 * if the corresponding point is inside or on the mesh and \a 0 otherwise.
 *
 * \return Return code is QUEST_INOUT_SUCCESS when all preconditions
 * are satisfied and QUEST_INOUT_FAILED otherwise.
 *
 * \pre inout_initialized() == true
 * \pre When \a npoints is greater than zero, arrays \a x, \a y, \a z
 *  and \a res are not \a nullptr and contain sufficient data/space for
 *  \a npoints points.
 */
int inout_evaluate(const double* x, const double* y, const double* z, int npoints, int* res);

/*!
 * \brief Returns the lower coordinates of the mesh's bounding box
 *
 * \param [in] coords A buffer for the returned coordinates
 * \pre \a coords != nullptr and has sufficient storage for the coordinates
 * \pre inout_initialized() == true
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_mesh_min_bounds(double* coords);

/*!
 * \brief Returns the upper coordinates of the mesh's bounding box
 *
 * \param [in] coords A buffer for the returned coordinates
 * \pre \a coords != nullptr and has sufficient storage for the coordinates
 * \pre inout_initialized() == true
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_mesh_max_bounds(double* coords);

/*!
 * \brief Returns the center of mass of the mesh
 *
 * The function computes a discrete center of mass defined by the average of the
 * mesh coordinates rather than a continuous center of mass defined by the mesh
 * faces.
 *
 * \param [in] coords A buffer for the returned coordinates
 * \pre \a coords != nullptr and has sufficient storage for the coordinates
 * \pre inout_initialized() == true
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 */
int inout_mesh_center_of_mass(double* coords);

/*!
 * \brief Gets the spatial dimension of the query
 *
 * \return Returns the spatial dimension for the inout query
 * \note This determines the number of coordinates for the input mesh
 * \note The default dimension is 3
 * \sa inout_set_dimension
 */
int inout_get_dimension();

/// @}

/// \name InOut query -- setup options and parameters
/// \note These must be called before \a inout_init()
/// @{

/*!
 * \brief Sets the spatial dimension of the query
 *
 * \param dim The dimension for the inout query
 * \return Return code is QUEST_INOUT_SUCCESS if successful, QUEST_INOUT_FAILED otherwise
 * \pre dim == 2 || dim == 3
 * \pre inout_initialized() == false
 */
int inout_set_dimension(int dim);

/*!
 * \brief Enables/disables verbose logging output
 *
 * By default, the logging verbosity is set to \a false.
 *
 * \param verbosity True for more verbose, false for less verbose
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 * \pre inout_initialized() == false
 */
int inout_set_verbose(bool verbosity);

/*!
 * \brief Sets the cutoff distance for welding vertices during initialization
 *
 * By default, the welding threshold is 1E-9.
 *
 * The inout query requires the input surface to be watertight so this
 * parameter should be set with care. A welding threshold that is too
 * high could unnecessarily merge vertices and create topological defects,
 * while a value that is too low risks leaving gaps in meshes with tolerances
 * between vertices. The default value tends to work well in practice.
 *
 * \param thresh Cutoff distance for welding vertices
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 * \pre inout_initialized() == false
 * \pre thresh > 0
 */
int inout_set_vertex_weld_threshold(double thresh);

/*!
 * \brief Sets the number of samples for each knot span (2D only)
 *
 * By default, the welding threshold is 25
 *
 * Span intervals are the segments of each curve. This parameter controls
 * the number of samples to use when linearizing each knot span
 * of the contour splines 
 *
 * \param segmentsPerKnotSpan The number of segments for each knot span
 * \return Return code is QUEST_INOUT_SUCCESS if successful
 *  and QUEST_INOUT_FAILED otherwise.
 * \pre inout_initialized() == false
 * \pre segmentsPerKnotSpan >= 1
 */
int inout_set_segments_per_knot_span(int segmentsPerKnotSpan);

/// @}

}  // end namespace quest
}  // end namespace axom

#endif  // QUEST_INOUT_INTERFACE_HPP_
