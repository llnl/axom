// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_SIGNED_DISTANCE_INTERFACE_HPP_
#define QUEST_SIGNED_DISTANCE_INTERFACE_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile-time definitions

// Quest includes
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"  // MPI_COMM_SELF

// C/C++ includes
#include <string>  // for std::string

/*!
 * \file
 *
 * \brief The Signed Distance Query evaluates the signed distance function,
 *  \f$ \phi \f$, at an arbitrary point, \f$ \vec{p} \f$, with respect to an
 *  oriented surface, \f$ \mathcal S\f$
 *
 *  Given a discrete representation of the surface, i.e., a surface mesh, that
 *  is oriented according to an outward-facing unit normal, the signed distance
 *  function, \f$ \phi(\vec{p}) \f$, evaluated at a point, \f$ \vec{p} \f$, is
 *  defined such that:
 *  \f{equation}{
 *    \phi( \vec{p} ) =
 *      \begin{cases}
 *        +d,   & \vec{p} \mbox{ is in front of the surface, } \mathcal{S}  \\
 *        \pm0, & \vec{p} \mbox{ is on the boundary, } \partial\mathcal{S}  \\
 *        -d,   & \vec{p} \mbox{ is behind the surface, } \mathcal{S}
 *       \end{cases}
 *  \f}
 *
 *  \warning Currently, the Signed Distance Query supports:
 *    * reading in surface meshes given in the
 *      <a href="https://en.wikipedia.org/wiki/STL_(file_format)">STL format</a>
 *    * oriented 3D triangular surface meshes
 *
 * <b> Usage Example: </b>
 * \code
 *
 *  int main( int argc, char** argv )
 *  {
 *    MPI_Init( &argc, &argv );
 *    MPI_Comm mpi_comm = MPI_COMM_WORLD;
 *    ...
 *    // STEP 0: set parameters
 *    quest::signed_distance_set_dimension( 3 );
 *    quest::signed_distance_set_max_levels( 25 );
 *    quest::signed_distance_set_max_occupancy( 5 );
 *
 *    // STEP 1: initialize
 *    quest::signed_distance_init( "path/to/stlmesh.stl", mpi_comm );
 *
 *    // STEP 2: evaluate the signed distance at N query points
 *    for ( int i=0; i < N; ++i )
 *    {
 *       ...
 *       phi[ i ] = quest::signed_distance_evaluate( x[i], y[i], z[i] );
 *
 *    } // END for all N query points
 *
 *    // STEP 3: Finalize
 *    quest::signed_distance_finalize( );
 *
 *    MPI_Finalize( );
 *    return 0;
 *  }
 *
 * \endcode
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
// Signed Distance execution policies
enum class SignedDistExec
{
  CPU = 0,
  OpenMP = 1,
  GPU = 2
};

/// \name Signed Distance Query Initialization Methods
/// @{

/*!
 * \brief Initializes the Signed Distance Query with a surface given in an
 *  <a href="https://en.wikipedia.org/wiki/STL_(file_format)">STL</a> formatted
 *  file.
 *
 * \param [in] file the name of the file consisting of the surface mesh.
 * \param [in] comm the MPI communicator (applicable when MPI is available)
 *
 * \return status zero on success, or a non-zero value if an error occurs
 *
 * \note The Signed Distance Query currently only supports reading in meshes
 *  defined in <a href="https://en.wikipedia.org/wiki/STL_(file_format)">STL</a>
 *  file format.
 *
 * \pre file.empty() == false
 * \pre comm != MPI_COMM_NULL (when MPI is available)
 * \pre signed_distance_initialized() == false
 * \post signed_distance_initialized() == true
 */
int signed_distance_init(const std::string& file, MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Initializes the Signed Distance Query with the given surface mesh
 *
 * \param [in] m pointer to the surface mesh object
 * \param [in] comm the MPI communicator (applicable whem MPI is available)
 *
 * \return status zero on success, or a non-zero value if an error occurs.
 *
 * \note The Signed Distance Query currently only supports 3-D triangular
 *  surface meshes.
 *
 * \pre m != nullptr
 * \pre comm != MPI_COMM_NULL (when MPI is available)
 * \pre signed_distance_initialized() == false
 * \post signed_distance_initialized() == true
 *
 * \see mint::Mesh
 */
int signed_distance_init(const mint::Mesh* m, MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Checks if the Signed Distance Query has been initialized
 * \return status true if initialized, else, false.
 */
bool signed_distance_initialized();

/// @}

/// \name Signed Distance Query Options
/// @{

/*!
 * \brief Sets the dimension for the Signed Distance Query.
 * \param [in] dim the dimension, e.g., 2 or 3
 *
 * \warning The Signed Distance function is currently supported in 3D
 * \note Options must be set before initializing the Signed Distance Query.
 */
void signed_distance_set_dimension(int dim);

/*!
 * \brief Indicates whether the input to the signed distance consists of a
 *  water-tight surface mesh, or not.
 *
 * \param [in] status flag indicating whether the input is water-tight
 *
 * \note By default the input type is assumed to be a water-tight surface mesh.
 * \note Options must be set before initializing the Signed Distance Query.
 *
 * \note When the input is not a closed surface mesh, the assumption is that
 *  the surface mesh divides the computational mesh domain into two regions.
 *  Hence, the surface mesh has to span the entire domain of interest, e.g.,
 *  the computational mesh at which the signed distance field is evaluated,
 *  along some plane.
 *
 * \warning The sign of the distance from a given query point is determined by
 *  a pseudo-normal which is computed at the closest point on the surface mesh.
 *  For a non-watertight mesh, the sign of the distance is not defined
 *  everywhere. Specifically, the sign is ambiguous for all points for which
 *  a normal projection onto the surface does not exist.
 *
 */
void signed_distance_set_closed_surface(bool status);

/*!
 * \brief Sets whether the distance query should compute or ignore the sign
 * \param [in] computeSign predicate indicating if sign should be computed
 *
 * \note Options must be set before initializing the Signed Distance Query.
 */
void signed_distance_set_compute_signs(bool computeSign);

/*!
 * \brief Sets the allocator to use for creating internal signed distance query
 *  data structures.
 * \param [in] allocatorID the allocator ID to use
 *
 * \note Allocator should be compatible with the set execution space.
 */
void signed_distance_set_allocator(int allocatorID);

/*!
 * \brief Enables/Disables verbose output for the Signed Distance Query.
 * \param [in] status flag indicating whether to enable/disable verbose output
 *
 * \note Options must be set before initializing the Signed Distance Query.
 *
 * \note Currently, this is only applicable when the Signed Distance Query
 *  initializes the SLIC logging environment.
 */
void signed_distance_set_verbose(bool status);

/*!
 * \brief Enable/Disable the use of MPI-3 on-node shared memory for storing
 *  the surface mesh. By default this option is disabled.
 *
 * \param [in] status flag indicating whether to enable/disable shared memory.
 *
 * \note This option utilities MPI-3 features
 */
void signed_distance_use_shared_memory(bool status);

/*!
 * \brief Set the execution space in which to run signed distance queries. By
 *  default this option is set to SIGNED_DIST_EVAL_CPU.
 *
 * \param [in] exec_space the execution space setting, one of the enum values
 *  in SignedDistExec.
 *
 * \note This function resets the default allocator based on the requested
 *  execution space. If a custom allocator ID is desired, set it in a call to
 *  axom::setDefaultAllocator after this call.
 */
void signed_distance_set_execution_space(SignedDistExec execSpace);

/// @}

/// \name Signed Distance Query Evaluation Methods
/// @{

/*!
 * \brief Evaluates the signed distance function at the given point.
 *
 * \param [in] x the x-coordinate of the point in query
 * \param [in] y the y-coordinate of the point in query
 * \param [in] z the z-coordinate of the point in query
 *
 * \return d the signed distance evaluated at the specified point.
 */
double signed_distance_evaluate(double x, double y, double z = 0.0);

/*!
 * \brief Evaluates the signed distance function at the given 3D point.
 *
 * \param [in] x the x-coordinate of the point in query
 * \param [in] y the y-coordinate of the point in query
 * \param [in] z the z-coordinate of the point in query
 * \param [out] cp_x the x-coordinate of the computed closest point on surface to query point
 * \param [out] cp_y the y-coordinate of the computed closest point on surface to query point
 * \param [out] cp_z the z-coordinate of the computed closest point on surface to query point
 * \param [out] n_x the x-coordinate of the surface normal at the computed closest point
 * \param [out] n_y the y-coordinate of the surface normal at the computed closest point
 * \param [out] n_z the z-coordinate of the surface normal at the computed closest point
 *
 * \return d the signed distance evaluated at the specified point. The closest surface point
 * to the query point and the normal at that point are returned as OUT parameters 
 */
double signed_distance_evaluate(double x,
                                double y,
                                double z,
                                double& cp_x,
                                double& cp_y,
                                double& cp_z,
                                double& n_x,
                                double& n_y,
                                double& n_z);

/*!
 * \brief Evaluates the signed distance function at the given set of points.
 *
 * \param [in] x array consisting of the x-coordinates for each query point
 * \param [in] y array consisting of the y-coordinates for each query point
 * \param [in] z array consisting of the z-coordinates for each query point
 * \param [in] npoints the number of query point
 * \param [out] phi output array storing the signed distance of each point
 *
 * \pre x != nullptr
 * \pre y != nullptr
 * \pre z != nullptr
 * \pre phi != nullptr
 */
void signed_distance_evaluate(const double* x,
                              const double* y,
                              const double* z,
                              int npoints,
                              double* phi);

/*!
 * \brief Computes the bounds of the specified input mesh supplied to the
 *  Signed Distance Query.
 *
 * \param [out] lo buffer to store the lower bound mesh coordinates.
 * \param [out] hi buffer to store the upper bound mesh coordinates.
 *
 * \pre lo != nullptr
 * \pre hi != nullptr
 * \pre hi & lo must point to a buffer that is at least ndims long.
 * \pre signed_distance_initialized() == true
 */
void signed_distance_get_mesh_bounds(double* lo, double* hi);

/// @}

/// \name Signed Distance Query Finalization Methods
/// @{

/*!
 * \brief Finalizes the SignedDistance query
 * \post signed_distance_initialized() == false.
 */
void signed_distance_finalize();

/// @}

}  // end namespace quest
}  // end namespace axom

#endif /* QUEST_SIGNED_DISTANCE_INTERFACE_HPP_ */
