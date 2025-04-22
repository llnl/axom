// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_HELPERS_HPP_
#define QUEST_HELPERS_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/slic.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/quest/interface/internal/mpicomm_wrapper.hpp"
#include "axom/quest/readers/STLReader.hpp"
#include "axom/quest/readers/ProEReader.hpp"

// C/C++ includes
#include <string>

/*!
 * \file QuestHelpers.hpp
 *
 * \brief Helper methods that can be used across the different Quest queries.
 */
namespace axom
{
namespace quest
{
namespace internal
{
constexpr int READ_FAILED = -1;
constexpr int READ_SUCCESS = 0;

/*!
 * \brief Simple RAII-based utility class to update the slic logging level within a scope
 *
 * The original logging level will be restored when this instance goes out of scope
 */
class ScopedLogLevelChanger
{
public:
  ScopedLogLevelChanger(slic::message::Level newLevel)
  {
    if(slic::isInitialized())
    {
      m_previousLevel = slic::getLoggingMsgLevel();
      slic::setLoggingMsgLevel(newLevel);
    }
  }

  ~ScopedLogLevelChanger()
  {
    if(slic::isInitialized())
    {
      slic::setLoggingMsgLevel(m_previousLevel);
    }
  }

private:
  slic::message::Level m_previousLevel {slic::message::Level::Debug};
};

/// \name MPI Helper/Wrapper Methods
/// @{
#ifdef AXOM_USE_MPI

/*!
 * \brief Deallocates the specified MPI window object.
 * \param [in] window handle to the MPI window.
 * \note All buffers attached to the window are also deallocated.
 */
void mpi_win_free(MPI_Win* window);

/*!
 * \brief Deallocates the specified MPI communicator object.
 * \param [in] comm handle to the MPI communicator object.
 */
void mpi_comm_free(MPI_Comm* comm);

/*!
 * \brief Reads the mesh on rank 0 and exchanges the mesh metadata, i.e., the
 *  number of nodes and faces with all other ranks.
 *
 * \param [in] global_rank_id MPI rank w.r.t. the global communicator
 * \param [in] global_comm handle to the global communicator
 * \param [in,out] reader the corresponding STL reader
 * \param [out] mesh_metadata an array consisting of the mesh metadata.
 *
 * \note This method calls read() on the reader on rank 0.
 *
 * \pre global_comm != MPI_COMM_NULL
 * \pre mesh_metadata != nullptr
 */
int read_and_exchange_mesh_metadata(int global_rank_id,
                                    MPI_Comm global_comm,
                                    quest::STLReader& reader,
                                    axom::IndexType mesh_metadata[2]);

#endif /* AXOM_USE_MPI */

#ifdef AXOM_USE_MPI3
/*!
 * \brief Creates inter-node and intra-node communicators from the given global
 *  MPI communicator handle.
 *
 *  The intra-node communicator groups the ranks within the same compute node.
 *  Consequently, all ranks have a global rank ID, w.r.t. the global
 *  communicator, and a corresponding local rank ID, w.r.t. the intra-node
 *  communicator. The global rank ID can span and is unique across multiple
 *  compute nodes, while the local rank ID, is only uniquely defined within the
 *  same compute node and is generally different from the global rank ID.
 *
 *  In contrast, the inter-node communicator groups only a subset of the ranks
 *  defined by the global communicator. Specifically, ranks that have a
 *  local rank ID of zero are included in the inter-node commuinicator and
 *  have a corresponding inter-comm rank ID, w.r.t., the inter-node
 *  communicator. For all other ranks, that are not included in the inter-node
 *  communicator, the inter-comm rank ID is set to "-1".
 *
 * \param [in]  global_comm handle to the global MPI communicator.
 * \param [out] intra_node_comm handle to the intra-node communicator object.
 * \param [out] inter_node_comm handle to the inter-node communicator object.
 * \param [out] global_rank_id rank ID w.r.t. the global communicator
 * \param [out] local_rank_id rank ID within the a compute node
 * \param [out] intercom_rank_id rank ID w.r.t. the inter-node communicator

 * \note The caller must call `MPI_Comm_free` on the corresponding communicator
 *  handles, namely, `intra_node_comm` and `inter_node_comm` which are created
 *  by this routine.
 *
 * \pre global_comm != MPI_COMM_NULL
 * \pre intra_node_comm == MPI_COMM_NULL
 * \pre inter_node_comm == MPI_COMM_NULL
 * \post intra_node_comm != MPI_COMM_NULL
 * \post inter_node_comm != MPI_COMM_NULL
 */
void create_communicators(MPI_Comm global_comm,
                          MPI_Comm& intra_node_comm,
                          MPI_Comm& inter_node_comm,
                          int& global_rank_id,
                          int& local_rank_id,
                          int& intercom_rank_id);
#endif

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)
/*!
 * \brief Allocates a shared memory buffer for the mesh that is shared among
 *  all the ranks within the same compute node.
 *
 * \param [in] intra_node_comm intra-node communicator within a node.
 * \param [in] mesh_metada tuple with the number of nodes/faces on the mesh
 * \param [out] x pointer into the buffer where the x--coordinates are stored.
 * \param [out] y pointer into the buffer where the y--coordinates are stored.
 * \param [out] z pointer into the buffer where the z--coordinates are stored.
 * \param [out] conn pointer into the buffer consisting the cell-connectivity.
 * \param [out] mesh_buffer raw buffer consisting of all the mesh data.
 * \param [out] shared_window MPI window to which the shared buffer is attached.
 *
 * \return bytesize the number of bytes in the raw buffer.
 *
 * \pre intra_node_comm != MPI_COMM_NULL
 * \pre mesh_metadata != nullptr
 * \pre x == nullptr
 * \pre y == nullptr
 * \pre z == nullptr
 * \pre conn == nullptr
 * \pre mesh_buffer == nullptr
 * \pre shared_window == MPI_WIN_NULL
 *
 * \post x != nullptr
 * \post y != nullptr
 * \post z != nullptr
 * \post coon != nullptr
 * \post mesh_buffer != nullptr
 * \post shared_window != MPI_WIN_NULL
 */
MPI_Aint allocate_shared_buffer(int local_rank_id,
                                MPI_Comm intra_node_comm,
                                const axom::IndexType mesh_metadata[2],
                                double*& x,
                                double*& y,
                                double*& z,
                                axom::IndexType*& conn,
                                unsigned char*& mesh_buffer,
                                MPI_Win& shared_window);
#endif

/// @}

/// \name Mesh I/O methods
/// @{

#if defined(AXOM_USE_MPI) && defined(AXOM_USE_MPI3)

/*!
 * \brief Reads in the surface mesh from the specified file into a shared
 *  memory buffer that is attached to the given MPI shared window.
 *
 * \param [in] file the file consisting of the surface mesh
 * \param [in] global_comm handle to the global MPI communicator
 * \param [out] mesh_buffer pointer to the raw mesh buffer
 * \param [out] m pointer to the mesh object
 * \param [out] intra_node_comm handle to the shared MPI communicator.
 * \param [out] shared_window handle to the MPI shared window.
 *
 * \return status set to READ_SUCCESS, or READ_FAILED on error.
 *
 * \note Each rank has a unique mint::Mesh object instance, however, the
 *  mint::Mesh object is constructed using external pointers that point into
 *  the supplied mesh_buffer, an on-node data-structure shared across all
 *  MPI ranks within the same compute node.
 *
 * \pre global_comm != MPI_COMM_NULL
 * \pre mesh_buffer == nullptr
 * \pre m == nullptr
 * \pre intra_node_comm == MPI_COMM_NULL
 * \pre shared_window == MPI_WIN_NULL
 *
 * \post m != nullptr
 * \post m->isExternal() == true
 * \post mesh_buffer != nullptr
 * \post intra_node_comm != MPI_COMM_NULL
 * \post shared_window != MPI_WIN_NULL
 */
int read_stl_mesh_shared(const std::string& file,
                         MPI_Comm global_comm,
                         unsigned char*& mesh_buffer,
                         mint::Mesh*& m,
                         MPI_Comm& intra_node_comm,
                         MPI_Win& shared_window);
#endif

/*!
 * \brief Reads in the surface mesh from the specified file.
 *
 * \param [in] file the file consisting of the surface
 * \param [out] m user-supplied pointer to point to the mesh object.
 * \param [in] comm the MPI communicator, only applicable when MPI is available.
 *
 * \note This method currently expects the surface mesh to be given in STL format.
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function.
 *
 * \return status set to zero on success, or to a non-zero value otherwise.
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::TRIANGLE
 *
 * \see STLReader
 * \see PSTLReader
 */
int read_stl_mesh(const std::string& file, mint::Mesh*& m, MPI_Comm comm = MPI_COMM_SELF);

#ifdef AXOM_USE_C2C
/*!
 * \brief Reads in the contour mesh from the specified file
 *
 * \param [in] file the file consisting of a C2C contour defined by one or more c2c::Piece
 * \param [in] transform A 4x4 matrix that contains a transform to be applied to points.
 * \param [in] segmentsPerPiece number of segments to sample per contour Piece
 * \param [in] vertexWeldThreshold threshold for welding vertices of adjacent curves
 * \param [out] m user-supplied pointer to point to the mesh object
 * \param [out] revolvedVolume An approximation of the revolved volume of the contour
 *                             or 0 if it could not be computed.
 * \param [in] comm the MPI communicator, only applicable when MPI is available
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function
 *
 * \return status set to zero on success, or to a non-zero value otherwise
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::SEGMENT
 * \post revolvedVolume > 0 if it could be computed.
 *
 * \see C2CReader
 * \see PC2CReader
 */
int read_c2c_mesh_uniform(const std::string& file,
                          const numerics::Matrix<double>& transform,
                          int segmentsPerPiece,
                          double vertexWeldThreshold,
                          mint::Mesh*& m,
                          double& revolvedVolume,
                          MPI_Comm comm = MPI_COMM_SELF);

/*!
 * \brief Reads in the contour mesh from the specified file and refines it
 *        according to a percent error.
 *
 * \param [in] file the file consisting of a C2C contour defined by one or more c2c::Piece
 * \param [in] transform A 4x4 matrix that contains a transform to be applied to points.
 * \param [in] percentError An acceptable percent error in the arc length of the curve.
 * \param [in] vertexWeldThreshold threshold for welding vertices of adjacent curves
 * \param [out] m user-supplied pointer to the mesh object
 * \param [out] revolvedVolume An approximation of the revolved volume of the contour
 *                             or 0 if it could not be computed.
 * \param [in] comm the MPI communicator, only applicable when MPI is available
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function
 *
 * \return status set to zero on success, or to a non-zero value otherwise
 *
 * \pre m == nullptr
 * \pre !file.empty()
 * \pre percentError should be in the range (0,1) non-inclusive.
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::SEGMENT
 * \post revolvedVolume > 0 if it could be computed.
 *
 * \see C2CReader
 * \see PC2CReader
 */
int read_c2c_mesh_non_uniform(const std::string& file,
                              const numerics::Matrix<double>& transform,
                              double percentError,
                              double vertexWeldThreshold,
                              mint::Mesh*& m,
                              double& revolvedVolume,
                              MPI_Comm comm = MPI_COMM_SELF);
#endif  // AXOM_USE_C2C

/*!
 * \brief Reads in the Pro/E tetrahedral mesh from the specified file.
 *
 * \param [in] file the file consisting of the Pro/E mesh
 * \param [out] m user-supplied pointer to point to the mesh object.
 * \param [in] comm the MPI communicator, only applicable when MPI is available.
 *
 * \note The caller is responsible for properly de-allocating the mesh object
 *  that is returned by this function.
 *
 * \return zero on success, or a non-zero value otherwise.
 *
 * \pre m == nullptr
 * \pre !file.empty()
 *
 * \post m != nullptr
 * \post m->getMeshType() == mint::UNSTRUCTURED_MESH
 * \post m->hasMixedCellTypes() == false
 * \post m->getCellType() == mint::TET
 *
 * \see ProEReader
 * \see PProEReader
 */
int read_pro_e_mesh(const std::string& file, mint::Mesh*& m, MPI_Comm comm = MPI_COMM_SELF);

/// @}

/// \name Mesh Helper Methods
/// @{

/*!
 * \brief Computes the bounds of the given mesh.
 *
 * \param [in] mesh pointer to the mesh whose bounds will be computed.
 * \param [out] lo buffer to store the lower bound mesh coordinates
 * \param [out] hi buffer to store the upper bound mesh coordinates
 *
 * \pre mesh != nullptr
 * \pre lo != nullptr
 * \pre hi != nullptr
 * \pre hi & lo must point to buffers that are at least N long, where N
 *  corresponds to the mesh dimension.
 */
void compute_mesh_bounds(const mint::Mesh* mesh, double* lo, double* hi);
/// @}

/// \name Logger Initialize/Finalize Methods
/// @{

/*!
 * \brief Helper method to initialize the Slic logger if needed.
 *
 * \param [in,out] isInitialized indicates if Slic is already initialized.
 * \param [out] mustFinalize inidicates if the caller would be responsible
 *  for finalizing the Slic logger.
 * \param [in] verbose flag to control the verbosity
 * \param [in] comm the MPI communicator (applicable when compiled with MPI)
 *
 * \note If Slic is not already initialized, this method will initialize the
 *  Slic Logging environment and set the `isInitialized` flag to true.
 *
 * \note The 'verbose' flag is only applicable when the Slic logging environment
 *  is not already initialized by the calling application. In that case, when
 *  'verbose' is true, all messages will get logged to the console, including,
 *  Info and debug messages. Otherwise, if 'false', only errors will be printed
 *  out.
 *
 *  \see logger_finalize
 */
void logger_init(bool& isInitialized, bool& mustFinalize, bool verbose, MPI_Comm comm);

/*!
 * \brief Finalizes the Slic logger (if needed)
 *
 * \param [in] mustFinalize flag that indicates whether the query is responsible
 *  for finalizing the Slic logger.
 *
 * \see logger_init
 */
void logger_finalize(bool mustFinalize);
/// @}

} /* end namespace internal */
} /* end namespace quest    */
} /* end namespace axom     */

#endif /* QUEST_HELPERS_HPP_ */
