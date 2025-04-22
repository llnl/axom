// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/mint/mesh/ConnectivityArray.hpp"
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/config.hpp"
#include "axom/slic.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
#endif

#include <algorithm>

namespace axom
{
namespace mint
{
const char IGNORE_OUTPUT[] = ".*";

constexpr IndexType vertex_stride = 1;
constexpr IndexType hex_stride = 8;

namespace internal
{
/*!
 * \brief Return the new ID capacity for a ConnectivityArray given an
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the amount the number of IDs will increase by
 */
IndexType calc_ID_capacity(const ConnectivityArray<NO_INDIRECTION>& connec, IndexType increase)
{
  IndexType new_n_IDs = connec.getNumberOfIDs() + increase;
  if(new_n_IDs > connec.getIDCapacity())
  {
    IndexType stride = connec.getNumberOfValuesForID();
    IndexType newCapacity =
      axom::utilities::max<IndexType>(connec.getValueCapacity() * connec.getResizeRatio() + 0.5,
                                      new_n_IDs * stride);
    IndexType remainder = newCapacity % stride;
    if(remainder != 0)
    {
      newCapacity += stride - remainder;
    }
    return newCapacity / stride;
  }

  return connec.getIDCapacity();
}

/*!
 * \brief Return the new ID capacity for a ConnectivityArray given an
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the amount the number of IDs will increase by
 */
IndexType calc_ID_capacity(const ConnectivityArray<TYPED_INDIRECTION>& connec, IndexType increase)
{
  IndexType new_n_IDs = connec.getNumberOfIDs() + increase;
  if(new_n_IDs > connec.getIDCapacity())
  {
    return axom::utilities::max<IndexType>(connec.getIDCapacity() * connec.getResizeRatio() + 0.5,
                                           new_n_IDs);
  }

  return connec.getIDCapacity();
}

/*!
 * \brief Return the new value capacity for a ConnectivityArray given an
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the amount the number of values will increase by
 */
IndexType calc_value_capacity(const ConnectivityArray<NO_INDIRECTION>& connec, IndexType increase)
{
  IndexType stride = getCellInfo(connec.getIDType()).num_nodes;
  return calc_ID_capacity(connec, increase / stride) * stride;
}

/*!
 * \brief Return the new value capacity for a ConnectivityArray given an
 *  increase.
 *
 * \param [in] connec the ConnectivityArray in question.
 * \param [in] increase the amount the number of values will increase by
 */
template <ConnectivityType TYPE>
IndexType calc_value_capacity(const ConnectivityArray<TYPE>& connec, IndexType increase)
{
  IndexType new_n_values = connec.getNumberOfValues() + increase;
  if(new_n_values > connec.getValueCapacity())
  {
    return axom::utilities::max<IndexType>(connec.getValueCapacity() * connec.getResizeRatio() + 0.5,
                                           new_n_values);
  }

  return connec.getValueCapacity();
}

/*!
 * \brief Get the type associated with the given ID to be inserted into a
 *  ConnectivityArray.
 *
 * \param [in] ID the ID in question.
 * \param [in] cell_type the CellType of the ConnectivityArray.
 * \param [in] types the array of types.
 */
CellType get_type(IndexType ID, CellType cell_type, const CellType* types)
{
  if(types == nullptr)
  {
    return cell_type;
  }
  else
  {
    return types[ID];
  }
}

/*!
 * \brief Get a pointer to the types starting at the given ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] types the array of types.
 */
const CellType* get_type_ptr(IndexType ID, const CellType* types)
{
  if(types == nullptr)
  {
    return nullptr;
  }
  else
  {
    return types + ID;
  }
}

/*!
 * \brief Return the number of values associated with ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] stride the stride of the ConnectivityArray.
 * \param [in] offsets the offsets array.
 */
IndexType get_n_values(IndexType ID, IndexType stride, const IndexType* offsets)
{
  if(offsets == nullptr)
  {
    return stride;
  }
  else
  {
    return offsets[ID + 1] - offsets[ID];
  }
}

/*!
 * \brief Return the number of values between the given IDs.
 *
 * \param [in] startID the start of the range.
 * \param [in] endID the end of the range.
 * \param [in] stride the stride of the ConnectivityArray.
 * \param [in] offsets pointer to the offsets array for the values in question.
 */
IndexType get_n_values(IndexType startID, IndexType endID, IndexType stride, const IndexType* offsets)
{
  if(offsets == nullptr)
  {
    return stride * (endID - startID);
  }
  else
  {
    return offsets[endID] - offsets[startID];
  }
}

/*!
 * \brief Return a pointer to the values array starting with the given ID.
 *
 * \param [in] ID the ID used to index into the values array.
 * \param [in] stride the stride of the ConnectivityArray.
 * \param [in] values the array of values.
 * \param [in] offsets pointer to the offsets array for the values in question.
 */
const IndexType* get_value_ptr(IndexType ID,
                               IndexType stride,
                               const IndexType* values,
                               const IndexType* offsets)
{
  if(offsets == nullptr)
  {
    return values + ID * stride;
  }
  else
  {
    return values + offsets[ID];
  }
}

/*!
 * \brief Return a pointer to the offsets array starting with the given ID.
 *
 * \param [in] ID the ID in question.
 * \param [in] offsets the offsets array.
 */
const IndexType* get_offset_ptr(IndexType ID, const IndexType* offsets)
{
  if(offsets == nullptr)
  {
    return nullptr;
  }
  else
  {
    return offsets + ID;
  }
}

/*!
 * \brief Check that the given pointers point to the same values or that they
 *  are both null.
 *
 * \param [in] p1 the first pointer to check.
 * \param [in] p2 the second pointer to check.
 * \param [in] n the number of values to check.
 */
template <typename T>
void check_pointers(const T* p1, const T* p2, IndexType n)
{
  if(p1 == nullptr && p2 == nullptr)
  {
    return;
  }

  for(IndexType i = 0; i < n; ++i)
  {
    EXPECT_EQ(p1[i], p2[i]) << "i = " << i;
  }
}

/*!
 * \brief Check that the given ConnectivityArrays are equal.
 *
 * \param [in] lhs the first ConnectivityArray to check.
 * \param [in] rhs the second ConnectivityArray to check.
 */
template <ConnectivityType TYPE>
void check_equality(const ConnectivityArray<TYPE>& lhs, const ConnectivityArray<TYPE>& rhs)
{
  ASSERT_EQ(lhs.getNumberOfIDs(), rhs.getNumberOfIDs());
  ASSERT_EQ(lhs.getNumberOfValues(), rhs.getNumberOfValues());
  EXPECT_EQ(lhs.getIDCapacity(), rhs.getIDCapacity());
  EXPECT_EQ(lhs.isExternal(), rhs.isExternal());
  EXPECT_EQ(lhs.isInSidre(), rhs.isInSidre());

  const IndexType n_IDs = lhs.getNumberOfIDs();
  for(IndexType ID = 0; ID < n_IDs; ++ID)
  {
    const IndexType cur_n_values = lhs.getNumberOfValuesForID(ID);
    ASSERT_EQ(cur_n_values, rhs.getNumberOfValuesForID(ID));
    EXPECT_EQ(lhs.getIDType(ID), rhs.getIDType(ID));
    check_pointers(lhs[ID], rhs[ID], cur_n_values);
  }
}

/*!
 * \brief Checks that the given values were properly appended to the
 *  ConnectivityArray.
 *
 * \param [in/out] connec the ConnectivityArray to check.
 * \param [in] n_IDs the number of IDs appended.
 * \param [in] initial_n_IDs the number of IDs in the ConnectivityArray before
 *  appending.
 * \param [in] values the array of values to check.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID that was appended.
 */
template <ConnectivityType TYPE>
void checkAppend(const ConnectivityArray<TYPE>& connec,
                 IndexType n_IDs,
                 IndexType initial_n_IDs,
                 const IndexType* values,
                 const IndexType* offsets = nullptr,
                 const CellType* types = nullptr)
{
  const CellType cell_type = connec.getIDType();
  const IndexType stride = (cell_type == UNDEFINED_CELL) ? -1 : getCellInfo(cell_type).num_nodes;

  for(IndexType ID = 0; ID < n_IDs; ++ID)
  {
    const IndexType cur_n_values = get_n_values(ID, stride, offsets);
    const IndexType* cur_values = get_value_ptr(ID, stride, values, offsets);
    const CellType cur_type = get_type(ID, cell_type, types);
    EXPECT_EQ(connec.getIDType(initial_n_IDs + ID), cur_type);
    EXPECT_EQ(connec.getNumberOfValuesForID(initial_n_IDs + ID), cur_n_values);
    check_pointers(cur_values, connec[initial_n_IDs + ID], cur_n_values);
  }
}

/*!
 * \brief Check that appending to the given ConnectivityArray functions as
 *  expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] values the array of values to append.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 */
template <ConnectivityType TYPE>
void append(ConnectivityArray<TYPE>& connec,
            IndexType n_IDs,
            const IndexType* values,
            const IndexType* offsets = nullptr,
            const CellType* types = nullptr)
{
  const CellType cell_type = connec.getIDType();
  const IndexType stride = (cell_type == UNDEFINED_CELL) ? -1 : getCellInfo(cell_type).num_nodes;
  const IndexType initial_n_IDs = connec.getNumberOfIDs();

  IndexType half_n_IDs = n_IDs / 2;
  for(IndexType ID = 0; ID < half_n_IDs; ++ID)
  {
    const IndexType cur_n_values = get_n_values(ID, stride, offsets);
    const IndexType* cur_values = get_value_ptr(ID, stride, values, offsets);
    const CellType type = get_type(ID, cell_type, types);

    connec.append(cur_values, cur_n_values, type);
    EXPECT_EQ(connec.getNumberOfIDs(), initial_n_IDs + ID + 1);
    EXPECT_EQ(connec.getIDType(initial_n_IDs + ID), type);
    EXPECT_EQ(connec.getNumberOfValuesForID(initial_n_IDs + ID), cur_n_values);
    check_pointers(cur_values, connec[initial_n_IDs + ID], cur_n_values);
  }

  const IndexType cur_ID = half_n_IDs;
  const IndexType* cur_values = get_value_ptr(cur_ID, stride, values, offsets);
  const IndexType* cur_offsets = get_offset_ptr(cur_ID, offsets);
  const CellType* cur_types = get_type_ptr(cur_ID, types);
  connec.appendM(cur_values, n_IDs - cur_ID, cur_offsets, cur_types);

  EXPECT_EQ(connec.getNumberOfIDs(), initial_n_IDs + n_IDs);
  checkAppend(connec, n_IDs, initial_n_IDs, values, offsets, types);
}

/*!
 * \brief Calls internal::append multiple times and checks that the output is as
 *  expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] max_IDs the max number of IDs to append in a given iteration.
 * \param [in] values the array of values to append.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 */
template <ConnectivityType TYPE>
void testAppend(ConnectivityArray<TYPE>& connec,
                IndexType max_IDs,
                const IndexType* values,
                const IndexType* offsets = nullptr,
                const CellType* types = nullptr)
{
  /* Append the values */
  for(IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs)
  {
    internal::append(connec, n_IDs, values, offsets, types);
  }

  /* Check that the values were appended properly */
  IndexType initial_n_IDs = 0;
  for(IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs)
  {
    checkAppend(connec, n_IDs, initial_n_IDs, values, offsets, types);
    initial_n_IDs += n_IDs;
  }
}

/*!
 * \brief Check that setting the values of the given ConnectivityArray
 *  functions as expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] initial_values the array of values to append.
 * \param [in] values the array of values to set.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 */
template <ConnectivityType TYPE>
void set(ConnectivityArray<TYPE>& connec,
         IndexType n_IDs,
         const IndexType* initial_values,
         const IndexType* values,
         const IndexType* offsets = nullptr,
         const CellType* types = nullptr)
{
  const CellType cell_type = connec.getIDType();
  const IndexType stride = (cell_type == UNDEFINED_CELL) ? -1 : getCellInfo(cell_type).num_nodes;
  const IndexType initial_n_IDs = connec.getNumberOfIDs();

  /* Append the initial values. */
  append(connec, n_IDs, initial_values, offsets, types);

  IndexType half_n_IDs = n_IDs / 2;
  for(IndexType ID = 0; ID < half_n_IDs; ++ID)
  {
    const IndexType cur_n_values = get_n_values(ID, stride, offsets);
    const IndexType* cur_values = get_value_ptr(ID, stride, values, offsets);
    const CellType type = get_type(ID, cell_type, types);

    connec.set(cur_values, initial_n_IDs + ID);
    EXPECT_EQ(connec.getIDType(initial_n_IDs + ID), type);
    EXPECT_EQ(connec.getNumberOfValuesForID(initial_n_IDs + ID), cur_n_values);
    check_pointers(cur_values, connec[initial_n_IDs + ID], cur_n_values);
  }

  const IndexType cur_ID = half_n_IDs;
  const IndexType* cur_values = get_value_ptr(cur_ID, stride, values, offsets);
  connec.setM(cur_values, initial_n_IDs + cur_ID, n_IDs - cur_ID);

  EXPECT_EQ(connec.getNumberOfIDs(), initial_n_IDs + n_IDs);
  checkAppend(connec, n_IDs, initial_n_IDs, values, offsets, types);
}

/*!
 * \brief Calls internal::set multiple times and checks that the output is as
 *  expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to and then set.
 * \param [in] max_IDs the max number of IDs to append in a given iteration.
 * \param [in] initial_values the array of values to append.
 * \param [in] values the array of values to set.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 */
template <ConnectivityType TYPE>
void testSet(ConnectivityArray<TYPE>& connec,
             IndexType max_IDs,
             const IndexType* initial_values,
             const IndexType* values,
             const IndexType* offsets = nullptr,
             const CellType* types = nullptr)
{
  for(IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs)
  {
    internal::set(connec, n_IDs, initial_values, values, offsets, types);
  }

  /* Check that the values were appended and set properly */
  IndexType initial_n_IDs = 0;
  for(IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs)
  {
    checkAppend(connec, n_IDs, initial_n_IDs, values, offsets, types);
    initial_n_IDs += n_IDs;
  }
}

/*!
 * \brief Check that inserting the values into the given ConnectivityArray
 *  functions as expected.
 *
 * \param [in/out] connec the ConnectivityArray to insert into.
 * \param [in] n_IDs the number of IDs to insert.
 * \param [in] values the array of values to set.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 *
 * \pre connec.empty()
 * \pre n_IDs % 6 == 0
 *
 * \note For n_IDs / 6 repetitions inserts an ID to the front, interior, and
 *  rear of the ConnectivityArray. Finally inserts n_IDs / 6 IDs each to the
 *  front, interior, and rear of the ConnectivityArray in one go.
 */
template <ConnectivityType TYPE>
void testInsert(ConnectivityArray<TYPE>& connec,
                IndexType n_IDs,
                const IndexType* values,
                const IndexType* offsets = nullptr,
                const CellType* types = nullptr)
{
  SLIC_ERROR_IF(!connec.empty(), "Insertion test requires an empty ConnectivityArray.");
  SLIC_ERROR_IF(n_IDs % 6 != 0, "Insertion test requires n_IDs to be a multiple of 6.");

  const CellType cell_type = connec.getIDType();
  const IndexType stride = (cell_type == UNDEFINED_CELL) ? -1 : getCellInfo(cell_type).num_nodes;

  const IndexType half_n_IDs = n_IDs / 2;
  const IndexType third_n_IDs = n_IDs / 3;
  const IndexType sixth_n_IDs = n_IDs / 6;

  /* Array of IDs the next ID to insert at the front, middle, and back. */
  IndexType IDs[3] = {third_n_IDs - 1, third_n_IDs, 2 * third_n_IDs};

  /* Array of positions at which to insert for the front, middle, and back. */
  IndexType insert_positions[3] = {0, 1, 2};

  for(IndexType round = 0; round < sixth_n_IDs; ++round)
  {
    for(IndexType i = 0; i < 3; ++i)
    {
      const IndexType cur_ID = IDs[i];
      const IndexType insert_pos = insert_positions[i];
      const IndexType cur_n_values = get_n_values(cur_ID, stride, offsets);
      const IndexType* cur_values = get_value_ptr(cur_ID, stride, values, offsets);
      const CellType type = get_type(cur_ID, cell_type, types);

      connec.insert(cur_values, insert_pos, cur_n_values, type);
      EXPECT_EQ(connec.getNumberOfIDs(), 3 * round + i + 1);
      EXPECT_EQ(connec.getIDType(insert_pos), type);
      EXPECT_EQ(connec.getNumberOfValuesForID(insert_pos), cur_n_values);
      check_pointers(cur_values, connec[insert_pos], cur_n_values);
    }

    /* The middle insertion increases by two, the back by three. */
    insert_positions[1] += 2;
    insert_positions[2] += 3;

    /* The front ID decreases by one while the middle and rear increase by one.
     */
    IDs[0]--;
    IDs[1]++;
    IDs[2]++;
  }

  /* Insert the rest of the front values. */
  const IndexType front_ID = 0;
  const IndexType* front_offsets = get_offset_ptr(front_ID, offsets);
  const IndexType* front_values = get_value_ptr(front_ID, stride, values, offsets);
  const CellType* front_types = get_type_ptr(front_ID, types);
  connec.insertM(front_values, front_ID, sixth_n_IDs, front_offsets, front_types);
  EXPECT_EQ(connec.getNumberOfIDs(), 2 * third_n_IDs);

  /* Insert the rest of the middle values. */
  const IndexType middle_ID = half_n_IDs;
  const IndexType* middle_offsets = get_offset_ptr(middle_ID, offsets);
  const IndexType* middle_values = get_value_ptr(middle_ID, stride, values, offsets);
  const CellType* middle_types = get_type_ptr(middle_ID, types);
  connec.insertM(middle_values, middle_ID, sixth_n_IDs, middle_offsets, middle_types);
  EXPECT_EQ(connec.getNumberOfIDs(), 5 * sixth_n_IDs);

  /* Insert the rest of the back values. */
  const IndexType back_ID = 5 * sixth_n_IDs;
  const IndexType* back_offsets = get_offset_ptr(back_ID, offsets);
  const IndexType* back_values = get_value_ptr(back_ID, stride, values, offsets);
  const CellType* back_types = get_type_ptr(back_ID, types);
  connec.insertM(back_values, back_ID, sixth_n_IDs, back_offsets, back_types);
  EXPECT_EQ(connec.getNumberOfIDs(), n_IDs);

  /* Check that the values were correctly inserted. */
  checkAppend(connec, n_IDs, 0, values, offsets, types);
}

/*!
 * \brief Check that capacity of the given ConnectivityArray functions as
 *  expected.
 *
 * \param [in/out] connec the ConnectivityArray to append to.
 * \param [in] n_IDs the number of IDs to append.
 * \param [in] values the array of values to append.
 * \param [in] offsets the offsets array.
 * \param [in] types the type of each ID to append.
 *
 * \pre connec.empty()
 *
 * \note offsets only need to be specified for ConnectivityArrays that need
 *  indirection and types only needs to be specified for ConnectivityArrays
 *  that have variable ID types. First reserves space for half of n_IDs and half
 *  the values, then appends that many IDs one at a time. Then appends one extra
 *  ID and checks that the resizing occured as expected. Then shrinks the
 *  ConnectivityArray and appends the remaining IDs all at once.
 */
template <ConnectivityType TYPE>
void testCapacity(ConnectivityArray<TYPE>& connec,
                  IndexType n_IDs,
                  const IndexType* values,
                  const IndexType* offsets = nullptr,
                  const CellType* types = nullptr)
{
  SLIC_ERROR_IF(!connec.empty(), "Insertion test requires an empty ConnectivityArray.");
  const CellType cell_type = connec.getIDType();
  const IndexType stride = (cell_type == UNDEFINED_CELL) ? -1 : getCellInfo(cell_type).num_nodes;
  const IndexType half_n_IDs = n_IDs / 2;
  const IndexType n_values = get_n_values(0, n_IDs, stride, offsets);
  const IndexType first_half_n_values = get_n_values(0, half_n_IDs, stride, offsets);
  IndexType cur_ID_size = connec.getNumberOfIDs();
  IndexType cur_value_size = connec.getNumberOfValues();

  EXPECT_EQ(cur_ID_size, 0);
  EXPECT_EQ(cur_value_size, 0);

  connec.reserve(half_n_IDs, first_half_n_values);
  const IndexType* cur_values_ptr = connec.getValuePtr();
  const IndexType* cur_offsets_ptr = connec.getOffsetPtr();
  const CellType* cur_types_ptr = connec.getTypePtr();
  IndexType cur_ID_capacity = connec.getIDCapacity();
  IndexType cur_value_capacity = connec.getValueCapacity();

  EXPECT_NE(cur_values_ptr, nullptr);
  EXPECT_GE(cur_ID_capacity, half_n_IDs);
  EXPECT_GE(cur_value_capacity, first_half_n_values);

  /* Append the first half of the IDs, no resize should occur. */
  for(IndexType ID = 0; ID < half_n_IDs; ++ID)
  {
    const IndexType cur_n_values = get_n_values(ID, stride, offsets);
    const IndexType* cur_values = get_value_ptr(ID, stride, values, offsets);
    const CellType type = get_type(ID, cell_type, types);

    connec.append(cur_values, cur_n_values, type);
    cur_ID_size += 1;
    cur_value_size += cur_n_values;
    EXPECT_EQ(connec.getNumberOfIDs(), cur_ID_size);
    EXPECT_EQ(connec.getIDCapacity(), cur_ID_capacity);
    EXPECT_EQ(connec.getNumberOfValues(), cur_value_size);
    EXPECT_EQ(connec.getValueCapacity(), cur_value_capacity);
    EXPECT_EQ(connec.getValuePtr(), cur_values_ptr);
    EXPECT_EQ(connec.getOffsetPtr(), cur_offsets_ptr);
    EXPECT_EQ(connec.getTypePtr(), cur_types_ptr);
    EXPECT_EQ(connec.getIDType(ID), type);
    EXPECT_EQ(connec.getNumberOfValuesForID(ID), cur_n_values);
    check_pointers(cur_values, connec[ID], cur_n_values);
  }

  /* Append one more value, should trigger a resize. */
  IndexType cur_ID = half_n_IDs;
  IndexType cur_n_values = get_n_values(cur_ID, stride, offsets);
  const IndexType* cur_values = get_value_ptr(cur_ID, stride, values, offsets);
  CellType type = get_type(cur_ID, cell_type, types);

  cur_ID_capacity = calc_ID_capacity(connec, 1);
  cur_value_capacity = calc_value_capacity(connec, cur_n_values);
  connec.append(cur_values, cur_n_values, type);
  cur_ID_size += 1;
  cur_value_size += cur_n_values;

  EXPECT_EQ(connec.getNumberOfIDs(), cur_ID_size);
  EXPECT_EQ(connec.getIDCapacity(), cur_ID_capacity);
  EXPECT_EQ(connec.getNumberOfValues(), cur_value_size);
  EXPECT_EQ(connec.getValueCapacity(), cur_value_capacity);
  EXPECT_EQ(connec.getIDType(cur_ID), type);
  EXPECT_EQ(connec.getNumberOfValuesForID(cur_ID), cur_n_values);
  check_pointers(cur_values, connec[cur_ID], cur_n_values);

  /* Shrink. */
  connec.shrink();
  cur_ID_capacity = cur_ID_size;
  cur_value_capacity = cur_value_size;
  EXPECT_EQ(connec.getIDCapacity(), cur_ID_capacity);
  EXPECT_EQ(connec.getValueCapacity(), cur_value_capacity);

  /* Append the rest of the values all at once, should trigger a resize. */
  cur_ID = half_n_IDs + 1;
  cur_values = get_value_ptr(cur_ID, stride, values, offsets);
  cur_n_values = get_n_values(cur_ID, n_IDs, stride, offsets);
  const IndexType* cur_offsets = get_offset_ptr(cur_ID, offsets);
  const CellType* cur_types = get_type_ptr(cur_ID, types);

  cur_ID_capacity = calc_ID_capacity(connec, n_IDs - cur_ID);
  cur_value_capacity = calc_value_capacity(connec, cur_n_values);
  connec.appendM(cur_values, n_IDs - cur_ID, cur_offsets, cur_types);
  EXPECT_EQ(connec.getNumberOfIDs(), n_IDs);
  EXPECT_EQ(connec.getIDCapacity(), cur_ID_capacity);
  EXPECT_EQ(connec.getNumberOfValues(), n_values);
  EXPECT_EQ(connec.getValueCapacity(), cur_value_capacity);

  for(IndexType ID = 0; ID < n_IDs; ++ID)
  {
    cur_n_values = get_n_values(ID, stride, offsets);
    cur_values = get_value_ptr(ID, stride, values, offsets);
    type = get_type(ID, cell_type, types);
    EXPECT_EQ(connec.getIDType(ID), type);
    check_pointers(cur_values, connec[ID], cur_n_values);
  }

  check_pointers(connec.getValuePtr(), values, n_values);
  check_pointers(connec.getOffsetPtr(), offsets, n_IDs + 1);
  check_pointers(connec.getTypePtr(), types, n_IDs);
}

/*!
 * \brief Create values in the range [0, n_vals) scaled by a factor.
 *
 * \param [in] n_vals the number of values to create.
 * \param [in] factor the scaling factor.
 */
const IndexType* createValues(IndexType n_vals, const IndexType factor = 1)
{
  IndexType* values = new IndexType[n_vals];
  for(IndexType i = 0; i < n_vals; ++i)
  {
    values[i] = factor * i;
  }

  return values;
}

/*!
 * \brief Create offsets and types arrays for the given number of IDs.
 *
 * \param [in] n_IDs the number of values to create.
 * \param [out] offsets where the offsets array is stored.
 * \param [out] types where the types array is stored.
 *
 * \return The number of values accounted for in the offsets array.
 */
IndexType createOffsetsAndTypes(IndexType n_IDs, IndexType*& offsets, CellType*& types)
{
  offsets = new IndexType[n_IDs + 1];
  types = new CellType[n_IDs];
  offsets[0] = 0;
  for(IndexType i = 0; i < n_IDs; ++i)
  {
    if(i % 2 == 0)
    {
      types[i] = VERTEX;
      offsets[i + 1] = offsets[i] + vertex_stride;
    }
    else
    {
      types[i] = HEX;
      offsets[i + 1] = offsets[i] + hex_stride;
    }
  }

  return offsets[n_IDs];
}

/*!
 * \brief Calaculate the total number of IDs and Values to insert given
 *  the maximum number inserted in a given iteration.
 *
 * \param [in] max_IDs the max number of IDs inserted in a single iteration.
 * \param [out] total_IDs the total number of IDs calculated.
 * \param [out] total_values the total number of values calculated.
 */
void calcTotalIDsAndValues(IndexType max_IDs, IndexType& total_IDs, IndexType& total_values)
{
  total_IDs = (max_IDs * (max_IDs + 1)) / 2;
  total_values = 0;
  for(IndexType n_IDs = 1; n_IDs <= max_IDs; ++n_IDs)
  {
    total_values += (n_IDs / 2) * (vertex_stride + hex_stride);
    if(n_IDs % 2 == 1)
    {
      total_values += vertex_stride;
    }
  }
}

/*!
 * \brief Check that the provided ConnectivityArray meets certain criteria.
 *
 * \param [in] connec the ConnectivityArray to check.
 * \param [in] external if the ConnectivityArray is external.
 * \param [in] n_IDs the number of IDs.
 * \param [in] n_values the number of values.
 * \param [in] pointer to the values array.
 * \param [in] pointer to the offsets array.
 * \param [in] pointer to the types array.
 */
template <ConnectivityType TYPE>
void checkConnectivity(const ConnectivityArray<TYPE>& connec,
                       bool external,
                       IndexType n_IDs,
                       IndexType n_values,
                       const IndexType* values = nullptr,
                       const IndexType* offsets = nullptr,
                       const CellType* types = nullptr)
{
  EXPECT_EQ(connec.isExternal(), external);
  EXPECT_EQ(connec.isInSidre(), false);
  EXPECT_EQ(connec.getNumberOfIDs(), n_IDs);
  EXPECT_EQ(connec.getNumberOfValues(), n_values);
  bool empty = n_IDs == 0 && n_values == 0;
  EXPECT_EQ(connec.empty(), empty);

  if(values != nullptr)
  {
    EXPECT_EQ(connec.getValuePtr(), values);
  }

  if(offsets != nullptr)
  {
    EXPECT_EQ(connec.getOffsetPtr(), offsets);
  }

  if(types != nullptr)
  {
    EXPECT_EQ(connec.getTypePtr(), types);
  }

#ifdef AXOM_MINT_USE_SIDRE
  EXPECT_EQ(connec.getGroup(), nullptr);
#endif
}

#ifdef AXOM_MINT_USE_SIDRE
/*!
 * \brief Check that the provided ConnectivityArray meets certain criteria.
 *
 * \param [in] connec the ConnectivityArray to check.
 * \param [in] n_IDs the number of IDs.
 * \param [in] n_values the number of values.
 * \param [in] group the sidr::Group of the ConnectivityArray.
 */
template <ConnectivityType TYPE>
void checkConnectivity(const ConnectivityArray<TYPE>& connec,
                       IndexType n_IDs,
                       IndexType n_values,
                       const sidre::Group* group = nullptr)
{
  EXPECT_EQ(connec.isExternal(), false);
  EXPECT_EQ(connec.isInSidre(), true);
  EXPECT_EQ(connec.getNumberOfIDs(), n_IDs);
  EXPECT_EQ(connec.getNumberOfValues(), n_values);
  bool empty = n_IDs == 0 && n_values == 0;
  EXPECT_EQ(connec.empty(), empty);
  EXPECT_EQ(connec.getGroup(), group);
}

/*!
 * \brief Check that the ConnectivityArray can be reconstructed from
 *  a sidre::Group.
 *
 * \param [in] connec the ConnectivityArray to check.
 */
template <ConnectivityType TYPE>
void checkSidreConstructor(const ConnectivityArray<TYPE>& connec)
{
  ASSERT_TRUE(connec.isInSidre());

  sidre::Group* group = const_cast<sidre::Group*>(connec.getGroup());
  ConnectivityArray<TYPE> cpy(group);
  internal::check_equality(connec, cpy);
}
#endif /* AXOM_MINT_USE_SIDRE */

/*!
 * \brief Check that the ConnectivityArray can be reconstructed from
 *  a external buffers.
 *
 * \param [in] connec the ConnectivityArray to check.
 */
void checkExternalConstructor(const ConnectivityArray<NO_INDIRECTION>& connec)
{
  ASSERT_TRUE(connec.isExternal());

  CellType cell_type = connec.getIDType();
  const IndexType n_IDs = connec.getNumberOfIDs();
  IndexType* values = const_cast<IndexType*>(connec.getValuePtr());
  ConnectivityArray<NO_INDIRECTION> connec_cpy(cell_type, n_IDs, values);
  internal::check_equality(connec, connec_cpy);
}

/*!
 * \brief Check that the ConnectivityArray can be reconstructed from
 *  a external buffers.
 *
 * \param [in] connec the ConnectivityArray to check.
 */
void checkExternalConstructor(const ConnectivityArray<TYPED_INDIRECTION>& connec)
{
  ASSERT_TRUE(connec.isExternal());

  const IndexType n_IDs = connec.getNumberOfIDs();
  IndexType* values = const_cast<IndexType*>(connec.getValuePtr());
  IndexType* offsets = const_cast<IndexType*>(connec.getOffsetPtr());
  CellType* types = const_cast<CellType*>(connec.getTypePtr());
  ConnectivityArray<TYPED_INDIRECTION> connec_cpy(n_IDs, values, offsets, types);
  internal::check_equality(connec, connec_cpy);
}

} /* end namespace internal */

/*******************************************************************************
 *                          Append tests                                       *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionNativeAppend)
{
  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = (max_IDs * (max_IDs + 1)) / 2;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray<NO_INDIRECTION> native_vertex(VERTEX);
  internal::checkConnectivity(native_vertex, false, 0, 0);

  ConnectivityArray<NO_INDIRECTION> native_hex(HEX);
  internal::checkConnectivity(native_vertex, false, 0, 0);

  internal::testAppend(native_vertex, max_IDs, values);
  internal::testAppend(native_hex, max_IDs, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(native_vertex, false, total_IDs, total_IDs * vertex_stride);
  internal::checkConnectivity(native_hex, false, total_IDs, total_IDs * hex_stride);

  delete[] values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, NoIndirectionExternalAppend)
{
  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = (max_IDs * (max_IDs + 1)) / 2;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[total_IDs * vertex_stride];
  ConnectivityArray<NO_INDIRECTION> ext_vertex(VERTEX, 0, ext_vertex_values, total_IDs);
  internal::checkConnectivity(ext_vertex, true, 0, 0, ext_vertex_values);

  IndexType* ext_hex_values = new IndexType[total_IDs * hex_stride];
  ConnectivityArray<NO_INDIRECTION> ext_hex(HEX, 0, ext_hex_values, total_IDs);
  internal::checkConnectivity(ext_hex, true, 0, 0, ext_hex_values);

  /* Check that appending functions properly */
  internal::testAppend(ext_vertex, max_IDs, values);
  internal::testAppend(ext_hex, max_IDs, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(ext_vertex, true, total_IDs, total_IDs * vertex_stride, ext_vertex_values);
  internal::checkConnectivity(ext_hex, true, total_IDs, total_IDs * hex_stride, ext_hex_values);

  /* Check that the external ConnectivityArrays cannot append any more. */
  EXPECT_DEATH_IF_SUPPORTED(ext_vertex.append(values), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(ext_hex.append(values), IGNORE_OUTPUT);

  /* Check that the external constructor functions properly. */
  internal::checkExternalConstructor(ext_vertex);
  internal::checkExternalConstructor(ext_hex);

  delete[] values;
  delete[] ext_vertex_values;
  delete[] ext_hex_values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionSidreAppend)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = (max_IDs * (max_IDs + 1)) / 2;

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Create the sidre storage ConnectivityArrays to be tested. */
  sidre::Group* vertex_group = root->createGroup("vertex");
  ConnectivityArray<NO_INDIRECTION> sidre_vertex(VERTEX, vertex_group, "test");
  internal::checkConnectivity(sidre_vertex, 0, 0, vertex_group);

  sidre::Group* hex_group = root->createGroup("hex");
  ConnectivityArray<NO_INDIRECTION> sidre_hex(HEX, hex_group, "test");
  internal::checkConnectivity(sidre_hex, 0, 0, hex_group);

  internal::testAppend(sidre_vertex, max_IDs, values);
  internal::testAppend(sidre_hex, max_IDs, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(sidre_vertex, total_IDs, total_IDs * vertex_stride, vertex_group);
  internal::checkConnectivity(sidre_hex, total_IDs, total_IDs * hex_stride, hex_group);

  /* Check that restoring from sidre functions properly */
  internal::checkSidreConstructor(sidre_vertex);
  internal::checkSidreConstructor(sidre_hex);

  delete[] values;
#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionNativeAppend)
{
  constexpr IndexType max_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(max_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Calculate the total number of IDs and values. */
  IndexType total_IDs;
  IndexType total_values;
  internal::calcTotalIDsAndValues(max_IDs, total_IDs, total_values);

  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray<TYPED_INDIRECTION> native_mixed;
  internal::checkConnectivity(native_mixed, false, 0, 0);

  /* Append the values */
  internal::testAppend(native_mixed, max_IDs, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(native_mixed, false, total_IDs, total_values);

  delete[] values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, IndirectionExternalAppend)
{
  constexpr IndexType max_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(max_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Calculate the total number of IDs and values. */
  IndexType total_IDs;
  IndexType total_values;
  internal::calcTotalIDsAndValues(max_IDs, total_IDs, total_values);

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_mixed_values = new IndexType[total_IDs * hex_stride];
  IndexType* ext_mixed_offsets = new IndexType[total_IDs + 1];
  CellType* ext_mixed_types = new CellType[total_IDs];
  ConnectivityArray<TYPED_INDIRECTION> ext_mixed(0,
                                                 ext_mixed_values,
                                                 ext_mixed_offsets,
                                                 ext_mixed_types,
                                                 total_IDs,
                                                 total_values);
  internal::checkConnectivity(ext_mixed, true, 0, 0, ext_mixed_values, ext_mixed_offsets, ext_mixed_types);

  /* Check that appending functions properly */
  internal::testAppend(ext_mixed, max_IDs, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(ext_mixed,
                              true,
                              total_IDs,
                              total_values,
                              ext_mixed_values,
                              ext_mixed_offsets,
                              ext_mixed_types);

  /* Check that the external ConnectivityArrays cannot append any more. */
  EXPECT_DEATH_IF_SUPPORTED(ext_mixed.append(values, 1, VERTEX), IGNORE_OUTPUT);

  /* Check that the external constructor functions properly. */
  internal::checkExternalConstructor(ext_mixed);

  delete[] values;
  delete[] offsets;
  delete[] types;
  delete[] ext_mixed_values;
  delete[] ext_mixed_offsets;
  delete[] ext_mixed_types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionSidreAppend)
{
#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType max_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(max_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Calculate the total number of IDs and values. */
  IndexType total_IDs;
  IndexType total_values;
  internal::calcTotalIDsAndValues(max_IDs, total_IDs, total_values);

  /* Create the sidre storage ConnectivityArrays to be tested. */
  sidre::Group* mixed_group = root->createGroup("mixed");
  ConnectivityArray<TYPED_INDIRECTION> sidre_mixed(mixed_group, "test");
  internal::checkConnectivity(sidre_mixed, 0, 0, mixed_group);

  /* Append the values */
  internal::testAppend(sidre_mixed, max_IDs, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(sidre_mixed, total_IDs, total_values, mixed_group);

  /* Check that restoring from sidre functions properly */
  internal::checkSidreConstructor(sidre_mixed);

  delete[] values;
  delete[] offsets;
  delete[] types;

#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

/*******************************************************************************
 *                          Set tests                                          *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionNativeSet)
{
  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = (max_IDs * (max_IDs + 1)) / 2;

  /* Allocate and populate the values buffer. */
  const IndexType* initial_values = internal::createValues(max_values, -1);
  const IndexType* values = internal::createValues(max_values);

  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray<NO_INDIRECTION> native_vertex(VERTEX);
  internal::checkConnectivity(native_vertex, false, 0, 0);

  ConnectivityArray<NO_INDIRECTION> native_hex(HEX);
  internal::checkConnectivity(native_vertex, false, 0, 0);

  /* Append and set the values */
  internal::testSet(native_vertex, max_IDs, initial_values, values);
  internal::testSet(native_hex, max_IDs, initial_values, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(native_vertex, false, total_IDs, total_IDs * vertex_stride);
  internal::checkConnectivity(native_hex, false, total_IDs, total_IDs * hex_stride);

  delete[] values;
  delete[] initial_values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionExternalSet)
{
  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = (max_IDs * (max_IDs + 1)) / 2;

  /* Allocate and populate the values buffer. */
  const IndexType* initial_values = internal::createValues(max_values, -1);
  const IndexType* values = internal::createValues(max_values);

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[total_IDs * vertex_stride];
  ConnectivityArray<NO_INDIRECTION> ext_vertex(VERTEX, 0, ext_vertex_values, total_IDs);
  internal::checkConnectivity(ext_vertex, true, 0, 0, ext_vertex_values);

  IndexType* ext_hex_values = new IndexType[total_IDs * hex_stride];
  ConnectivityArray<NO_INDIRECTION> ext_hex(HEX, 0, ext_hex_values, total_IDs);
  internal::checkConnectivity(ext_hex, true, 0, 0, ext_hex_values);

  /* Append and set the values */
  internal::testSet(ext_vertex, max_IDs, initial_values, values);
  internal::testSet(ext_hex, max_IDs, initial_values, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(ext_vertex, true, total_IDs, total_IDs * vertex_stride, ext_vertex_values);
  internal::checkConnectivity(ext_hex, true, total_IDs, total_IDs * hex_stride, ext_hex_values);

  /* Check that the external constructor functions properly. */
  internal::checkExternalConstructor(ext_vertex);
  internal::checkExternalConstructor(ext_hex);

  delete[] values;
  delete[] initial_values;
  delete[] ext_vertex_values;
  delete[] ext_hex_values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionSidreSet)
{
#ifdef AXOM_MINT_USE_SIDRE
  constexpr IndexType max_IDs = 100;
  constexpr IndexType max_values = hex_stride * max_IDs;
  constexpr IndexType total_IDs = (max_IDs * (max_IDs + 1)) / 2;

  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  /* Allocate and populate the values buffer. */
  const IndexType* initial_values = internal::createValues(max_values, -1);
  const IndexType* values = internal::createValues(max_values);

  /* Create the sidre storage ConnectivityArrays to be tested. */
  sidre::Group* vertex_group = root->createGroup("vertex");
  ConnectivityArray<NO_INDIRECTION> sidre_vertex(VERTEX, vertex_group, "test");
  internal::checkConnectivity(sidre_vertex, 0, 0, vertex_group);

  sidre::Group* hex_group = root->createGroup("hex");
  ConnectivityArray<NO_INDIRECTION> sidre_hex(HEX, hex_group, "test");
  internal::checkConnectivity(sidre_hex, 0, 0, hex_group);

  /* Append and set the values */
  internal::testSet(sidre_vertex, max_IDs, initial_values, values);
  internal::testSet(sidre_hex, max_IDs, initial_values, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(sidre_vertex, total_IDs, total_IDs * vertex_stride, vertex_group);
  internal::checkConnectivity(sidre_hex, total_IDs, total_IDs * hex_stride, hex_group);

  /* Check that restoring from sidre functions properly */
  internal::checkSidreConstructor(sidre_vertex);
  internal::checkSidreConstructor(sidre_hex);

  delete[] values;
  delete[] initial_values;
#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionNativeSet)
{
  constexpr IndexType max_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(max_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* initial_values = internal::createValues(max_values, -1);
  const IndexType* values = internal::createValues(max_values);

  /* Calculate the total number of IDs and values. */
  IndexType total_IDs;
  IndexType total_values;
  internal::calcTotalIDsAndValues(max_IDs, total_IDs, total_values);

  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray<TYPED_INDIRECTION> native_mixed;
  internal::checkConnectivity(native_mixed, false, 0, 0);

  /* Append and set the values */
  internal::testSet(native_mixed, max_IDs, initial_values, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(native_mixed, false, total_IDs, total_values);

  delete[] values;
  delete[] initial_values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionExternalSet)
{
  constexpr IndexType max_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(max_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* initial_values = internal::createValues(max_values, -1);
  const IndexType* values = internal::createValues(max_values);

  /* Calculate the total number of IDs and values. */
  IndexType total_IDs;
  IndexType total_values;
  internal::calcTotalIDsAndValues(max_IDs, total_IDs, total_values);

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_mixed_values = new IndexType[total_IDs * hex_stride];
  IndexType* ext_mixed_offsets = new IndexType[total_IDs + 1];
  CellType* ext_mixed_types = new CellType[total_IDs];
  ConnectivityArray<TYPED_INDIRECTION> ext_mixed(0,
                                                 ext_mixed_values,
                                                 ext_mixed_offsets,
                                                 ext_mixed_types,
                                                 total_IDs,
                                                 total_values);
  internal::checkConnectivity(ext_mixed, true, 0, 0, ext_mixed_values, ext_mixed_offsets, ext_mixed_types);

  /* Append and set the values */
  internal::testSet(ext_mixed, max_IDs, initial_values, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(ext_mixed,
                              true,
                              total_IDs,
                              total_values,
                              ext_mixed_values,
                              ext_mixed_offsets,
                              ext_mixed_types);

  /* Check that the external constructor functions properly. */
  internal::checkExternalConstructor(ext_mixed);

  delete[] initial_values;
  delete[] values;
  delete[] offsets;
  delete[] types;
  delete[] ext_mixed_values;
  delete[] ext_mixed_offsets;
  delete[] ext_mixed_types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionSidreSet)
{
#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType max_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(max_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* initial_values = internal::createValues(max_values, -1);
  const IndexType* values = internal::createValues(max_values);

  /* Calculate the total number of IDs and values. */
  IndexType total_IDs;
  IndexType total_values;
  internal::calcTotalIDsAndValues(max_IDs, total_IDs, total_values);

  /* Create the sidre storage ConnectivityArrays to be tested. */
  sidre::Group* mixed_group = root->createGroup("mixed");
  ConnectivityArray<TYPED_INDIRECTION> sidre_mixed(mixed_group, "test");
  internal::checkConnectivity(sidre_mixed, 0, 0, mixed_group);

  /* Append and set the values */
  internal::testSet(sidre_mixed, max_IDs, initial_values, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(sidre_mixed, total_IDs, total_values, mixed_group);

  /* Check that restoring from sidre functions properly */
  internal::checkSidreConstructor(sidre_mixed);

  delete[] values;
  delete[] initial_values;
  delete[] offsets;
  delete[] types;

#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

/*******************************************************************************
 *                          Insert tests                                       *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionNativeInsert)
{
  constexpr IndexType total_IDs = 500 * 6;
  constexpr IndexType max_values = hex_stride * total_IDs;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray<NO_INDIRECTION> native_vertex(VERTEX);
  internal::checkConnectivity(native_vertex, false, 0, 0);

  ConnectivityArray<NO_INDIRECTION> native_hex(HEX);
  internal::checkConnectivity(native_vertex, false, 0, 0);

  /* Insert the values */
  internal::testInsert(native_vertex, total_IDs, values);
  internal::testInsert(native_hex, total_IDs, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(native_vertex, false, total_IDs, total_IDs * vertex_stride);
  internal::checkConnectivity(native_hex, false, total_IDs, total_IDs * hex_stride);

  delete[] values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, NoIndirectionExternalInsert)
{
  constexpr IndexType total_IDs = 500 * 6;
  constexpr IndexType max_values = hex_stride * total_IDs;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_vertex_values = new IndexType[total_IDs * vertex_stride];
  ConnectivityArray<NO_INDIRECTION> ext_vertex(VERTEX, 0, ext_vertex_values, total_IDs);
  internal::checkConnectivity(ext_vertex, true, 0, 0, ext_vertex_values);

  IndexType* ext_hex_values = new IndexType[total_IDs * hex_stride];
  ConnectivityArray<NO_INDIRECTION> ext_hex(HEX, 0, ext_hex_values, total_IDs);
  internal::checkConnectivity(ext_hex, true, 0, 0, ext_hex_values);

  /* Insert the values */
  internal::testInsert(ext_vertex, total_IDs, values);
  internal::testInsert(ext_hex, total_IDs, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(ext_vertex, true, total_IDs, total_IDs * vertex_stride, ext_vertex_values);
  internal::checkConnectivity(ext_hex, true, total_IDs, total_IDs * hex_stride, ext_hex_values);

  /* Check that the external ConnectivityArrays cannot insert any more. */
  EXPECT_DEATH_IF_SUPPORTED(ext_vertex.insert(values, 0), IGNORE_OUTPUT);
  EXPECT_DEATH_IF_SUPPORTED(ext_hex.insert(values, 0), IGNORE_OUTPUT);

  /* Check that the external constructor functions properly. */
  internal::checkExternalConstructor(ext_vertex);
  internal::checkExternalConstructor(ext_hex);

  delete[] values;
  delete[] ext_vertex_values;
  delete[] ext_hex_values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, NoIndirectionSidreInsert)
{
#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType total_IDs = 500 * 6;
  constexpr IndexType max_values = hex_stride * total_IDs;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  /* Create the sidre storage ConnectivityArrays to be tested. */
  sidre::Group* vertex_group = root->createGroup("vertex");
  ConnectivityArray<NO_INDIRECTION> sidre_vertex(VERTEX, vertex_group, "test");
  internal::checkConnectivity(sidre_vertex, 0, 0, vertex_group);

  sidre::Group* hex_group = root->createGroup("hex");
  ConnectivityArray<NO_INDIRECTION> sidre_hex(HEX, hex_group, "test");
  internal::checkConnectivity(sidre_hex, 0, 0, hex_group);

  /* Insert the values */
  internal::testInsert(sidre_vertex, total_IDs, values);
  internal::testInsert(sidre_hex, total_IDs, values);

  /* Check that the number of IDs is correct. */
  internal::checkConnectivity(sidre_vertex, total_IDs, total_IDs * vertex_stride, vertex_group);
  internal::checkConnectivity(sidre_hex, total_IDs, total_IDs * hex_stride, hex_group);

  /* Check that restoring from sidre functions properly */
  internal::checkSidreConstructor(sidre_vertex);
  internal::checkSidreConstructor(sidre_hex);

  delete[] values;
#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionNativeInsert)
{
  constexpr IndexType total_IDs = 500 * 6;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType total_values = internal::createOffsetsAndTypes(total_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(total_values);

  /* Create the native storage ConnectivityArrays to be tested. */
  ConnectivityArray<TYPED_INDIRECTION> native_mixed;
  internal::checkConnectivity(native_mixed, false, 0, 0);

  /* Insert the values. */
  internal::testInsert(native_mixed, total_IDs, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(native_mixed, false, total_IDs, total_values);

  delete[] values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, IndirectionExternalInsert)
{
  constexpr IndexType total_IDs = 500 * 6;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType total_values = internal::createOffsetsAndTypes(total_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(total_values);

  /* Create the external storage ConnectivityArrays to be tested. */
  IndexType* ext_mixed_values = new IndexType[total_IDs * hex_stride];
  IndexType* ext_mixed_offsets = new IndexType[total_IDs + 1];
  CellType* ext_mixed_types = new CellType[total_IDs];
  ConnectivityArray<TYPED_INDIRECTION> ext_mixed(0,
                                                 ext_mixed_values,
                                                 ext_mixed_offsets,
                                                 ext_mixed_types,
                                                 total_IDs,
                                                 total_values);
  internal::checkConnectivity(ext_mixed, true, 0, 0, ext_mixed_values, ext_mixed_offsets, ext_mixed_types);

  /* Insert the values. */
  internal::testInsert(ext_mixed, total_IDs, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(ext_mixed,
                              true,
                              total_IDs,
                              total_values,
                              ext_mixed_values,
                              ext_mixed_offsets,
                              ext_mixed_types);

  /* Check that the external ConnectivityArrays cannot append any more. */
  EXPECT_DEATH_IF_SUPPORTED(ext_mixed.insert(values, 0, 1, VERTEX), IGNORE_OUTPUT);

  /* Check that the external constructor functions properly. */
  internal::checkExternalConstructor(ext_mixed);

  delete[] values;
  delete[] offsets;
  delete[] types;
  delete[] ext_mixed_values;
  delete[] ext_mixed_offsets;
  delete[] ext_mixed_types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array, IndirectionSidreInsert)
{
#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType total_IDs = 500 * 6;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType total_values = internal::createOffsetsAndTypes(total_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(total_values);

  /* Create the sidre storage ConnectivityArrays to be tested. */
  sidre::Group* mixed_group = root->createGroup("mixed");
  ConnectivityArray<TYPED_INDIRECTION> sidre_mixed(mixed_group, "test");
  internal::checkConnectivity(sidre_mixed, 0, 0, mixed_group);

  /* Insert the values */
  internal::testInsert(sidre_mixed, total_IDs, values, offsets, types);

  /* Check that the number of IDs and values is correct. */
  internal::checkConnectivity(sidre_mixed, total_IDs, total_values, mixed_group);

  /* Check that restoring from sidre functions properly */
  internal::checkSidreConstructor(sidre_mixed);

  delete[] values;
  delete[] offsets;
  delete[] types;
#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

/*******************************************************************************
 *                          Capacity tests                                     *
 ******************************************************************************/

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, NoIndirectionNativeCapacity)
{
  constexpr IndexType n_IDs = 100;
  constexpr IndexType max_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  for(double ratio = 0.75; ratio <= 3.0; ratio += 0.25)
  {
    /* Create the native storage ConnectivityArrays to be tested. */
    ConnectivityArray<NO_INDIRECTION> native_vertex(VERTEX);
    native_vertex.setResizeRatio(ratio);
    internal::checkConnectivity(native_vertex, false, 0, 0);

    ConnectivityArray<NO_INDIRECTION> native_hex(HEX);
    native_hex.setResizeRatio(ratio);
    internal::checkConnectivity(native_vertex, false, 0, 0);

    /* Append the values */
    if(ratio < 1.0)
    {
      EXPECT_DEATH_IF_SUPPORTED(internal::testCapacity(native_vertex, n_IDs, values), IGNORE_OUTPUT);
      EXPECT_DEATH_IF_SUPPORTED(internal::testCapacity(native_hex, n_IDs, values), IGNORE_OUTPUT);
    }
    else
    {
      internal::testCapacity(native_vertex, n_IDs, values);
      internal::testCapacity(native_hex, n_IDs, values);
    }
  }

  delete[] values;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, NoIndirectionSidreCapacity)
{
#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType n_IDs = 100;
  constexpr IndexType max_values = hex_stride * n_IDs;

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  for(double ratio = 0.75; ratio <= 3.0; ratio += 0.25)
  {
    /* Create the sidre storage ConnectivityArrays to be tested. */
    sidre::Group* vertex_group = root->createGroup("vertex");
    ConnectivityArray<NO_INDIRECTION> sidre_vertex(VERTEX, vertex_group, "test");
    sidre_vertex.setResizeRatio(ratio);
    internal::checkConnectivity(sidre_vertex, 0, 0, vertex_group);

    sidre::Group* hex_group = root->createGroup("hex");
    ConnectivityArray<NO_INDIRECTION> sidre_hex(HEX, hex_group, "test");
    sidre_hex.setResizeRatio(ratio);
    internal::checkConnectivity(sidre_hex, 0, 0, hex_group);

    /* Append the values */
    if(ratio < 1.0)
    {
      EXPECT_DEATH_IF_SUPPORTED(internal::testCapacity(sidre_vertex, n_IDs, values), IGNORE_OUTPUT);
      EXPECT_DEATH_IF_SUPPORTED(internal::testCapacity(sidre_hex, n_IDs, values), IGNORE_OUTPUT);
    }
    else
    {
      internal::testCapacity(sidre_vertex, n_IDs, values);
      internal::testCapacity(sidre_hex, n_IDs, values);
    }

    root->destroyGroups();
  }

  delete[] values;
#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, IndirectionNativeCapacity)
{
  constexpr IndexType n_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(n_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  for(double ratio = 0.75; ratio <= 3.0; ratio += 0.25)
  {
    /* Create the native storage ConnectivityArrays to be tested. */
    ConnectivityArray<TYPED_INDIRECTION> native_mixed(0, 0);
    native_mixed.setResizeRatio(ratio);
    internal::checkConnectivity(native_mixed, false, 0, 0);

    /* Append the values */
    if(ratio < 1.0)
    {
      EXPECT_DEATH_IF_SUPPORTED(internal::testCapacity(native_mixed, n_IDs, values, offsets, types),
                                IGNORE_OUTPUT);
    }
    else
    {
      internal::testCapacity(native_mixed, n_IDs, values, offsets, types);
    }
  }

  delete[] values;
  delete[] offsets;
  delete[] types;
}

//------------------------------------------------------------------------------
TEST(mint_connectivity_array_DeathTest, IndirectionSidreCapacity)
{
#ifdef AXOM_MINT_USE_SIDRE
  sidre::DataStore ds;
  sidre::Group* root = ds.getRoot();

  constexpr IndexType n_IDs = 100;

  /* Allocate and populate the types and offsets buffer. */
  IndexType* offsets;
  CellType* types;
  const IndexType max_values = internal::createOffsetsAndTypes(n_IDs, offsets, types);

  /* Allocate and populate the values buffer. */
  const IndexType* values = internal::createValues(max_values);

  for(double ratio = 0.75; ratio <= 3.0; ratio += 0.25)
  {
    /* Create the sidre storage ConnectivityArrays to be tested. */
    sidre::Group* mixed_group = root->createGroup("mixed");
    ConnectivityArray<TYPED_INDIRECTION> sidre_mixed(mixed_group, "test", 0, 0);
    sidre_mixed.setResizeRatio(ratio);
    internal::checkConnectivity(sidre_mixed, 0, 0, mixed_group);

    /* Append the values */
    if(ratio < 1.0)
    {
      EXPECT_DEATH_IF_SUPPORTED(internal::testCapacity(sidre_mixed, n_IDs, values, offsets, types),
                                IGNORE_OUTPUT);
    }
    else
    {
      internal::testCapacity(sidre_mixed, n_IDs, values, offsets, types);
    }

    root->destroyGroups();
  }

  delete[] values;
  delete[] offsets;
  delete[] types;
#else
  EXPECT_TRUE(true);
#endif /* AXOM_MINT_USE_SIDRE */
}

} /* end namespace mint */
} /* end namespace axom */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();
  return result;
}
