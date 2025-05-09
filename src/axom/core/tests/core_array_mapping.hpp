// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/Array.hpp"
#include "axom/core/MDMapping.hpp"
#include "axom/fmt.hpp"

// Google test include
#include "gtest/gtest.h"

// Test strides and permutations.
TEST(core_array_mapping, core_strides_and_permutations)
{
  {
    // 1D
    constexpr int DIM = 1;
    axom::StackArray<axom::IndexType, DIM> lengths {2};

    axom::MDMapping<DIM> colIndexer(lengths, axom::ArrayStrideOrder::COLUMN);
    EXPECT_EQ(colIndexer.getStrideOrder(), axom::ArrayStrideOrder::BOTH);
    axom::StackArray<std::uint16_t, DIM> colSlowestDirs {0};
    axom::StackArray<axom::IndexType, DIM> colStrides {1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(colIndexer.slowestDirs()[d], colSlowestDirs[d]);
      EXPECT_EQ(colIndexer.strides()[d], colStrides[d]);
    }
    EXPECT_TRUE(colIndexer == (axom::MDMapping<DIM>(lengths, colSlowestDirs)));

    axom::MDMapping<DIM> rowIndexer(lengths, axom::ArrayStrideOrder::ROW);
    EXPECT_EQ(rowIndexer.getStrideOrder(), axom::ArrayStrideOrder::BOTH);
    axom::StackArray<std::uint16_t, DIM> rowSlowestDirs {0};
    axom::StackArray<axom::IndexType, DIM> rowStrides {1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(rowIndexer.slowestDirs()[d], rowSlowestDirs[d]);
      EXPECT_EQ(rowIndexer.strides()[d], rowStrides[d]);
    }
    EXPECT_TRUE(rowIndexer == (axom::MDMapping<DIM>(lengths, rowSlowestDirs)));
  }
  {
    // 2D
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};

    axom::MDMapping<DIM> colIndexer(lengths, axom::ArrayStrideOrder::COLUMN);
    EXPECT_EQ(colIndexer.getStrideOrder(), axom::ArrayStrideOrder::COLUMN);
    axom::StackArray<std::uint16_t, DIM> colSlowestDirs {1, 0};
    axom::StackArray<axom::IndexType, DIM> colStrides {1, 3};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(colIndexer.slowestDirs()[d], colSlowestDirs[d]);
      EXPECT_EQ(colIndexer.strides()[d], colStrides[d]);
    }

    axom::MDMapping<DIM> rowIndexer(lengths, axom::ArrayStrideOrder::ROW);
    EXPECT_EQ(rowIndexer.getStrideOrder(), axom::ArrayStrideOrder::ROW);
    axom::StackArray<std::uint16_t, DIM> rowSlowestDirs {0, 1};
    axom::StackArray<axom::IndexType, DIM> rowStrides {2, 1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(rowIndexer.slowestDirs()[d], rowSlowestDirs[d]);
      EXPECT_EQ(rowIndexer.strides()[d], rowStrides[d]);
    }
  }
  {
    // 3D
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};

    axom::MDMapping<DIM> colIndexer(lengths, axom::ArrayStrideOrder::COLUMN);
    EXPECT_EQ(colIndexer.getStrideOrder(), axom::ArrayStrideOrder::COLUMN);
    axom::StackArray<std::uint16_t, DIM> colSlowestDirs {2, 1, 0};
    axom::StackArray<axom::IndexType, DIM> colStrides {1, 5, 15};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(colIndexer.slowestDirs()[d], colSlowestDirs[d]);
      EXPECT_EQ(colIndexer.strides()[d], colStrides[d]);
    }

    axom::MDMapping<DIM> rowIndexer(lengths, axom::ArrayStrideOrder::ROW);
    EXPECT_EQ(rowIndexer.getStrideOrder(), axom::ArrayStrideOrder::ROW);
    axom::StackArray<std::uint16_t, DIM> rowSlowestDirs {0, 1, 2};
    axom::StackArray<axom::IndexType, DIM> rowStrides {6, 2, 1};
    for(int d = 0; d < DIM; ++d)
    {
      EXPECT_EQ(rowIndexer.slowestDirs()[d], rowSlowestDirs[d]);
      EXPECT_EQ(rowIndexer.strides()[d], rowStrides[d]);
    }
  }
}

// Test row-major offsets.
TEST(quest_array_mapping, quest_row_major_offset)
{
  {
    // 1D
    constexpr int DIM = 1;
    axom::StackArray<axom::IndexType, DIM> lengths {3};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::ROW);
    EXPECT_EQ(mdmap.getStrideOrder(), axom::ArrayStrideOrder::BOTH);
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      axom::StackArray<axom::IndexType, DIM> indices {i};
      EXPECT_EQ(mdmap.toFlatIndex(indices), offset);
      EXPECT_EQ(mdmap.toMultiIndex(offset), indices);
      ++offset;
    }
  }
  {
    // 2D
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::ROW);
    EXPECT_EQ(mdmap.getStrideOrder(), axom::ArrayStrideOrder::ROW);
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      for(int j = 0; j < lengths[1]; ++j)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(mdmap.toFlatIndex(indices), offset);
        EXPECT_EQ(mdmap.toMultiIndex(offset), indices);
        ++offset;
      }
    }
  }
  {
    // 3D
    axom::StackArray<axom::IndexType, 3> lengths {5, 3, 2};
    axom::MDMapping<3> mdmap(lengths, axom::ArrayStrideOrder::ROW);
    EXPECT_EQ(mdmap.getStrideOrder(), axom::ArrayStrideOrder::ROW);
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      for(int j = 0; j < lengths[1]; ++j)
      {
        for(int k = 0; k < lengths[2]; ++k)
        {
          axom::StackArray<axom::IndexType, 3> indices {i, j, k};
          EXPECT_EQ(mdmap.toFlatIndex(indices), offset);
          EXPECT_EQ(mdmap.toMultiIndex(offset), indices);
          ++offset;
        }
      }
    }
  }
}

// Test column-major offsets.
TEST(quest_array_mapping, quest_col_major_offset)
{
  {
    // 1D
    constexpr int DIM = 1;
    axom::StackArray<axom::IndexType, DIM> lengths {3};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::COLUMN);
    EXPECT_EQ(mdmap.getStrideOrder(), axom::ArrayStrideOrder::BOTH);
    axom::IndexType offset = 0;
    for(int i = 0; i < lengths[0]; ++i)
    {
      axom::StackArray<axom::IndexType, DIM> indices {i};
      EXPECT_EQ(mdmap.toFlatIndex(indices), offset);
      EXPECT_EQ(mdmap.toMultiIndex(offset), indices);
      ++offset;
    }
  }
  {
    // 2D
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::COLUMN);
    EXPECT_EQ(mdmap.getStrideOrder(), axom::ArrayStrideOrder::COLUMN);
    axom::IndexType offset = 0;
    for(int j = 0; j < lengths[1]; ++j)
    {
      for(int i = 0; i < lengths[0]; ++i)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(mdmap.toFlatIndex(indices), offset);
        EXPECT_EQ(mdmap.toMultiIndex(offset), indices);
        ++offset;
      }
    }
  }
  {
    // 3D
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::COLUMN);
    EXPECT_EQ(mdmap.getStrideOrder(), axom::ArrayStrideOrder::COLUMN);
    axom::IndexType offset = 0;
    for(int k = 0; k < lengths[2]; ++k)
    {
      for(int j = 0; j < lengths[1]; ++j)
      {
        for(int i = 0; i < lengths[0]; ++i)
        {
          axom::StackArray<axom::IndexType, DIM> indices {i, j, k};
          EXPECT_EQ(mdmap.toFlatIndex(indices), offset);
          EXPECT_EQ(mdmap.toMultiIndex(offset), indices);
          ++offset;
        }
      }
    }
  }
}

template <int DIM>
typename std::enable_if<DIM == 1>::type check_arbitrary_strides_nested_loops(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs,
  const axom::MDMapping<DIM>& mdmap)
{
  /*
    We address the array with natural indices, but arange the
    nested loops to increase the offset by 1 each time around.
  */
  axom::StackArray<axom::IndexType, DIM> i;
  auto& n = i[fastestDirs[0]];
  axom::IndexType offset = 0;
  for(n = 0; n < lengths[fastestDirs[0]]; ++n)
  {
    // SLIC_INFO(axom::fmt::format("offset {}   i {}", offset, i));
    EXPECT_EQ(mdmap.toMultiIndex(offset), i);
    EXPECT_EQ(mdmap.toFlatIndex(i), offset);
    ++offset;
  }
}

template <int DIM>
typename std::enable_if<DIM == 2>::type check_arbitrary_strides_nested_loops(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs,
  const axom::MDMapping<DIM>& mdmap)
{
  /*
    We address the array with natural indices, but arange the
    nested loops to increase the offset by 1 each time around.
  */
  axom::StackArray<axom::IndexType, DIM> ij;
  auto& m = ij[fastestDirs[1]];
  auto& n = ij[fastestDirs[0]];
  axom::IndexType offset = 0;
  for(m = 0; m < lengths[fastestDirs[1]]; ++m)
  {
    for(n = 0; n < lengths[fastestDirs[0]]; ++n)
    {
      // SLIC_INFO(axom::fmt::format("offset {}   ij {}", offset, ij));
      EXPECT_EQ(mdmap.toMultiIndex(offset), ij);
      EXPECT_EQ(mdmap.toFlatIndex(ij), offset);
      ++offset;
    }
  }
}

template <int DIM>
typename std::enable_if<DIM == 3>::type check_arbitrary_strides_nested_loops(
  const axom::StackArray<axom::IndexType, DIM>& lengths,
  const axom::StackArray<axom::IndexType, DIM>& fastestDirs,
  const axom::MDMapping<DIM>& mdmap)
{
  /*
    We address the array with natural indices, but arange the
    nested loops to increase the offset by 1 each time around.
  */
  axom::StackArray<axom::IndexType, DIM> ijk;
  auto& l = ijk[fastestDirs[2]];  // Slowest dir
  auto& m = ijk[fastestDirs[1]];
  auto& n = ijk[fastestDirs[0]];  // Fastest dir
  axom::IndexType offset = 0;
  for(l = 0; l < lengths[fastestDirs[2]]; ++l)
  {
    for(m = 0; m < lengths[fastestDirs[1]]; ++m)
    {
      for(n = 0; n < lengths[fastestDirs[0]]; ++n)
      {
        // SLIC_INFO(axom::fmt::format("offset {}   ijk {}", offset, ijk));
        EXPECT_EQ(mdmap.toMultiIndex(offset), ijk);
        EXPECT_EQ(mdmap.toFlatIndex(ijk), offset);
        ++offset;
      }
    }
  }
}

template <int DIM>
void check_arbitrary_strides(const axom::StackArray<axom::IndexType, DIM>& lengths,
                             const axom::StackArray<axom::IndexType, DIM>& fastestDirs)
{
  // fastestDirs should be a permutation.
  // SLIC_INFO(axom::fmt::format("Testing lengths {} with fastestDirs {}", lengths, fastestDirs));

  axom::StackArray<axom::IndexType, DIM> strides;
  axom::IndexType currentStride = 1;
  for(int d = 0; d < DIM; ++d)
  {
    axom::IndexType currentDir = fastestDirs[d];
    strides[currentDir] = currentStride;
    currentStride *= lengths[currentDir];
  }

  axom::MDMapping<DIM> mdmap(strides);

  const auto slowestDirs = mdmap.slowestDirs();
  for(int d = 0; d < DIM; ++d)
  {
    EXPECT_EQ(slowestDirs[d], fastestDirs[DIM - d - 1]);
  }

  check_arbitrary_strides_nested_loops(lengths, fastestDirs, mdmap);
}

// Test arbitrary strides.
TEST(core_array_mapping, core_arbitrary_strides)
{
  {
    constexpr int DIM = 1;
    // MDMapping is overkill for 1D, but think of it as a smoke test.
    axom::StackArray<axom::IndexType, DIM> lengths {3};

    axom::StackArray<axom::IndexType, DIM> fastestDirs {0};
    check_arbitrary_strides(lengths, fastestDirs);
  }
  {
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};

    // Row-major ordering.
    axom::StackArray<axom::IndexType, DIM> rowFastestDirs {0, 1};
    check_arbitrary_strides(lengths, rowFastestDirs);

    // Col-major ordering.
    axom::StackArray<axom::IndexType, DIM> colFastestDirs {1, 0};
    check_arbitrary_strides(lengths, colFastestDirs);

    // A few random orderings.
    axom::StackArray<axom::IndexType, DIM> rand1FastestDirs {0, 1};
    check_arbitrary_strides(lengths, rand1FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand2FastestDirs {1, 0};
    check_arbitrary_strides(lengths, rand2FastestDirs);
  }
  {
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};

    // Row-major ordering.
    axom::StackArray<axom::IndexType, DIM> rowFastestDirs {0, 1, 2};
    check_arbitrary_strides(lengths, rowFastestDirs);

    // Col-major ordering.
    axom::StackArray<axom::IndexType, DIM> colFastestDirs {2, 1, 0};
    check_arbitrary_strides(lengths, colFastestDirs);

    // A few random orderings.
    axom::StackArray<axom::IndexType, DIM> rand1FastestDirs {0, 2, 1};
    check_arbitrary_strides(lengths, rand1FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand2FastestDirs {2, 0, 1};
    check_arbitrary_strides(lengths, rand2FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand3FastestDirs {1, 0, 2};
    check_arbitrary_strides(lengths, rand3FastestDirs);
    axom::StackArray<axom::IndexType, DIM> rand4FastestDirs {1, 2, 0};
    check_arbitrary_strides(lengths, rand4FastestDirs);
  }
}

// Test row- and column-major element offsets with Array's.
TEST(core_array_mapping, core_array_match)
{
  // No test for 1D.  Array provides no non-trivial interface to test.

  {
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::ROW);
    axom::Array<int, DIM> a(lengths, mdmap.getStrideOrder());
    axom::IndexType sequenceNum = 0;
    for(axom::IndexType i = 0; i < lengths[0]; ++i)
    {
      for(axom::IndexType j = 0; j < lengths[1]; ++j)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(mdmap.toFlatIndex(indices), sequenceNum);
        ++sequenceNum;
        EXPECT_EQ(&a(i, j), &a.flatIndex(mdmap.toFlatIndex(indices)));
      }
    }
  }
  {
    constexpr int DIM = 2;
    axom::StackArray<axom::IndexType, DIM> lengths {3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::COLUMN);
    axom::Array<int, DIM> a(lengths, mdmap.getStrideOrder());
    axom::IndexType sequenceNum = 0;
    for(axom::IndexType j = 0; j < lengths[1]; ++j)
    {
      for(axom::IndexType i = 0; i < lengths[0]; ++i)
      {
        axom::StackArray<axom::IndexType, DIM> indices {i, j};
        EXPECT_EQ(mdmap.toFlatIndex(indices), sequenceNum);
        ++sequenceNum;
        EXPECT_EQ(&a(i, j), &a.flatIndex(mdmap.toFlatIndex(indices)));
      }
    }
  }
  {
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::ROW);
    axom::Array<int, DIM> a(lengths, mdmap.getStrideOrder());
    axom::IndexType sequenceNum = 0;
    for(axom::IndexType i = 0; i < lengths[0]; ++i)
    {
      for(axom::IndexType j = 0; j < lengths[1]; ++j)
      {
        for(axom::IndexType k = 0; k < lengths[2]; ++k)
        {
          axom::StackArray<axom::IndexType, DIM> indices {i, j, k};
          EXPECT_EQ(mdmap.toFlatIndex(indices), sequenceNum);
          ++sequenceNum;
          EXPECT_EQ(&a(i, j, k), &a.flatIndex(mdmap.toFlatIndex(indices)));
        }
      }
    }
  }
  {
    constexpr int DIM = 3;
    axom::StackArray<axom::IndexType, DIM> lengths {5, 3, 2};
    axom::MDMapping<DIM> mdmap(lengths, axom::ArrayStrideOrder::COLUMN);
    axom::Array<int, DIM> a(lengths, mdmap.getStrideOrder());
    axom::IndexType sequenceNum = 0;
    for(axom::IndexType k = 0; k < lengths[2]; ++k)
    {
      for(axom::IndexType j = 0; j < lengths[1]; ++j)
      {
        for(axom::IndexType i = 0; i < lengths[0]; ++i)
        {
          axom::StackArray<axom::IndexType, DIM> indices {i, j, k};
          EXPECT_EQ(mdmap.toFlatIndex(indices), sequenceNum);
          ++sequenceNum;
          EXPECT_EQ(&a(i, j, k), &a.flatIndex(mdmap.toFlatIndex(indices)));
        }
      }
    }
  }
}
