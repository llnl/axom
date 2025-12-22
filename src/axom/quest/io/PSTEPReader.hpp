// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PSTEPREADER_HPP_
#define QUEST_PSTEPREADER_HPP_

#include "axom/config.hpp"

#ifndef AXOM_USE_OPENCASCADE
  #error PSTEPReader should only be included when Axom is configured with opencascade
#endif

#include "axom/core/Macros.hpp"
#include "axom/quest/io/STEPReader.hpp"

#include "mpi.h"

namespace axom
{
namespace quest
{

/*
 * \class PSTEPReader
 *
 * \brief Parallel version of STEPReader.
 *
 * Rank 0 reads and processes the STEP file using STEPReader and then
 * broadcasts mesh-level data (triangle mesh or NURBS patches) to the
 * other ranks in the communicator.
 */
class PSTEPReader : public STEPReader
{
public:
  PSTEPReader() = delete;
  explicit PSTEPReader(MPI_Comm comm);
  ~PSTEPReader() override = default;

  /*!
   * \brief Reads the STEP file on rank 0, and broadcasts success or failure.
   *
   * Rank 0:
   *  - constructs internal::StepFileProcessor,
   *  - validates the model (optionally),
   *  - extracts NURBS patches and trimming curves into m_patches.
   *
   * All ranks:
   *  - receive integer status flag from rank 0.
   *  - receives all patches and their trimming curves from rank 0
   *
   * \param validate_model When true, runs OpenCascade BRep validation on rank 0
   * \return 0 on success on all ranks, non-zero otherwise
   */
  int read(bool validate_model) override;

  /*!
  * \brief Generates a triangulated representation of the STEP file as a Mint unstructured triangle mesh.
  *
  * \param[inout] mesh Pointer to a Mint unstructured mesh that will be populated
  *            with triangular elements approximating the STEP geometry.
  * \param[in] linear_deflection Maximum allowed deviation between the
  *            original geometry and the triangulated approximation.
  * \param[in] angular_deflection Maximum allowed angular deviation (in radians)
  *            between normals of adjacent triangles.
  * \param[in] is_relative When false (default), linear deflection is in mesh units. When true,
              linear deflection is relative to the local edge length of the triangles.
  * \param[in] trimmed If true (default), the triangulation respects trimming curves.
  *            otherwise, we triangulate the untrimmed patches. The latter is mostly to aid 
  *            in understanding the model's patches and is not generally useful.
  *
  * The mesh is constructed on rank 0 and is then broadcast to the other ranks.
  */
  int getTriangleMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* mesh,
                      double linear_deflection = 0.1,
                      double angular_deflection = 0.5,
                      bool is_relative = false,
                      bool trimmed = true) override;

private:
  /// MPI broadcasts an integer from rank 0 and returns the value to all ranks
  /// (other ranks do not need to supply a value)
  int bcast_int(int value = 0);

  /// MPI broadcasts an integer from rank 0 and returns the value to all ranks
  /// (other ranks do not need to supply a value)
  axom::IndexType bcast_index(axom::IndexType value = axom::IndexType {});

  /// MPI broadcasts a bool from rank 0 and returns the value to all ranks
  /// (other ranks do not need to supply a value)
  bool bcast_bool(bool value = false);

  /// MPI broadcasts the values of an ArrayView<ValueType, ARR_DIM>
  /// Assumes all ranks already have the correct size in \a arr
  template <typename ValueType, int DIM>
  void bcast_data(axom::ArrayView<ValueType, DIM> arr)
  {
    static_assert(
      std::is_same_v<ValueType, double> || std::is_same_v<ValueType, int> ||
        std::is_same_v<ValueType, axom::IndexType>,
      "PSTEPReader::bcast_data only supports ValueType of double, int, or axom::IndexType");

    MPI_Bcast(arr.data(), arr.size(), axom::mpi_traits<ValueType>::type, 0, m_comm);
  }

  /// MPI broadcasts an Array<double, ARR_DIM> or Array<PointType, ARR_DIM> for ARR_DIM==1 or ARR_DIM==2
  template <typename ArrayType>
  void bcast_array(ArrayType& arr)
  {
    // Check that we have a 1D or 2D axom::Array
    constexpr int ARR_DIM = axom::detail::ArrayTraits<ArrayType>::dimension;
    static_assert(ARR_DIM == 1 || ARR_DIM == 2);

    // Check that the value_type of the array is either double or Point<double, DIM>
    using value_type = typename ArrayType::value_type;
    static_assert(std::is_same_v<value_type, double> || primal::detail::is_point_v<value_type>);

    const bool is_root = (m_my_rank == 0);

    // first, send/receive the size(s)
    if(is_root)
    {
      if constexpr(ARR_DIM == 1)
      {
        bcast_int(arr.size());
      }
      else if constexpr(ARR_DIM == 2)
      {
        bcast_int(arr.shape()[0]);
        bcast_int(arr.shape()[1]);
      }
    }
    else
    {
      arr.clear();
      if constexpr(ARR_DIM == 1)
      {
        const int size = bcast_int();
        arr.resize(ArrayOptions::Uninitialized {}, size);
      }
      else if constexpr(ARR_DIM == 2)
      {
        const int size0 = bcast_int();
        const int size1 = bcast_int();
        arr.resize(ArrayOptions::Uninitialized {}, size0, size1);
      }
    }

    // then, send/receive the data
    if constexpr(std::is_same_v<value_type, double>)
    {
      // handles Array<double,1> and Array<double,2>
      bcast_data(arr.view());
    }
    else if constexpr(primal::detail::is_point_v<value_type>)
    {
      // handles 1D or 2D array of primal::Point
      // since the data is contiguous, we cast the values to a 1D ArrayView<double>
      constexpr static int POINT_DIM = value_type::DIMENSION;
      const int numVals = arr.size() * POINT_DIM;
      double* dataPtr = nullptr;
      if constexpr(ARR_DIM == 1)  // 1D array of points
      {
        dataPtr = (arr.size() > 0) ? arr[0].data() : nullptr;
      }
      else if constexpr(ARR_DIM == 2)  // 2D array of points
      {
        dataPtr = (arr.size() > 0) ? &arr(0, 0)[0] : nullptr;
      }
      axom::ArrayView<double> flatView(dataPtr, numVals);
      bcast_data(flatView);
    }
    else
    {
      SLIC_ERROR("Unsupported ArrayType dimensionality in PSTEPReader::bcast_array");
    }
  }

private:
  MPI_Comm m_comm {MPI_COMM_NULL};
  int m_my_rank {0};
  int m_num_ranks {-1};

  DISABLE_COPY_AND_ASSIGNMENT(PSTEPReader);
  DISABLE_MOVE_AND_ASSIGNMENT(PSTEPReader);
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_PSTEPREADER_HPP_