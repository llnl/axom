// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SIDRE_ARRAY_HPP_
#define SIDRE_ARRAY_HPP_

#include "axom/core/utilities/Utilities.hpp"  // for memory allocation functions
#include "axom/core/Array.hpp"                // to inherit
#include "axom/core/Types.hpp"

#include "axom/slic/interface/slic.hpp"  // for slic logging macros

#include "View.hpp"    // for View definition
#include "Buffer.hpp"  // for Buffer definition

// C/C++ includes
#include <cstring>  // for std::memcpy

namespace axom
{
namespace sidre
{
/* Provided so that 0 doesn't convert to nullptr and lead to ambiguous
 * constructor calls. */
namespace internal
{
constexpr axom::IndexType ZERO = 0;
}

namespace detail
{
inline void describeViewImpl(TypeID T_type, const StackArray<axom::IndexType, 1>& shape, View* view)
{
  SLIC_ASSERT(view != nullptr);
  IndexType dims[1];
  dims[0] = shape[0];

  view->apply(T_type, 1, dims);
}

inline void describeViewImpl(TypeID T_type, const StackArray<axom::IndexType, 2>& shape, View* view)
{
  SLIC_ASSERT(view != nullptr);

  IndexType dims[2];
  dims[0] = shape[0];
  dims[1] = shape[1];

  view->apply(T_type, 2, dims);
}

template <int DIM>
inline IndexType getViewShapeImpl(int dim, const View* view);

template <>
inline IndexType getViewShapeImpl<1>(int dim, const View* view)
{
  SLIC_ERROR_IF(dim > 0, "Only one dimensional views supported.");
  SLIC_ERROR_IF(view->isEmpty(), "view cannot be empty.");
  SLIC_ERROR_IF(view->getNumDimensions() != 1, "view must have dimension 1.");

  sidre::IndexType dims[2];
  view->getShape(1, dims);
  return static_cast<axom::IndexType>(dims[dim]);
}

template <>
inline IndexType getViewShapeImpl<2>(int dim, const View* view)
{
  SLIC_ERROR_IF(dim > 1, "Only two dimensional views supported.");
  SLIC_ERROR_IF(view->isEmpty(), "view cannot be empty.");
  SLIC_ERROR_IF(view->getNumDimensions() != 2, "view must have dimension 2.");

  sidre::IndexType dims[2];
  view->getShape(2, dims);
  return static_cast<axom::IndexType>(dims[dim]);
}

}  // namespace detail

/*!
 * \class Array
 *
 * \brief Provides a generic multi-component array, contained in Sidre.
 *
 *  This sidre::Array class extends axom::Array by storing
 *  data in a Sidre `DataStore`.  This class provides a generic
 *  multi-component array container with dynamic re-allocation and insertion.
 *  Each element in the array is a tuple consisting of 1 or more components,
 *  which are stored contiguously.
 *
 *  Objects of the sidre::Array class may be constructed from a View.
 *  All array operations can be performed as with the base
 *  axom::Array class.  The size of the Array can grow as needed,
 *  and all memory management is delegated to Sidre.
 *
 *  \note When the Array object is deleted, it does not delete the associated
 *   data in Sidre, since, Sidre owns the data.
 *
 * \warning Reallocations tend to be costly operations in terms of performance.
 *  Use `reserve()` when the number of nodes is known a priori, or opt to
 *  use a constructor that takes an actual size and capacity when possible.
 *
 * \tparam T the type of the values to hold.
 *
 * \see Group
 * \see View
 */
template <typename T, int DIM = 1>
class Array : public axom::Array<T, DIM>
{
public:
  static_assert(DIM <= 2, "Only 1- and 2-dimensional Sidre arrays are permitted");
  /*!
   * \brief Default constructor. Disabled.
   */
  Array() = delete;

  /*!
   * \brief Copy constructor.

   * Deleted because copies would have to reference the same underlying
   * Sidre buffer, which does not match the underlying axom::Array ownership model
   */
  Array(const Array&) = delete;

  /*!
   * \brief Move constructor.
   * \param [in] other The array to move from
   */
  Array(Array&& other);

  /// \name Sidre Array constructors
  /// @{

  /*!
   * \brief Creates an Array instance from a View that already has data.
   *
   * \param [in] view the View that holds this Array's data.
   *
   * \note The Sidre view shape has two dimensions. The first dimension
   *  corresponds to the max capacity of the array and the second corresponds to
   *  the number of components per tuple.
   *
   * \pre view != nullptr
   * \pre view->isEmpty() == false.
   * \pre view->getNumDimensions() == 2
   *
   * \post capacity() == view->getDimension(0)
   * \post numComponents() == view->getDimension(1)
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <int SFINAE = DIM, typename std::enable_if<SFINAE == 1>::type* = nullptr>
  Array(View* view);

  /// \overload
  template <int SFINAE = DIM, typename std::enable_if<SFINAE == 2>::type* = nullptr>
  Array(View* view);

  /*!
   * \brief Creates an Array instance of `num_elements` size
   *  and populates the associated View.
   *
   * \param [in] view the View that will hold this Array's data.
   * \param [in] num_elements the number of values.
   * \param [in] capacity the number of values to allocate space for.
   *
   * \note The capacity argument is optional. If not specified or if less than
   *  num_elements, the capacity of the array will be initialized to
   *  num_elements * DEFAULT_RESIZE_RATIO.
   *
   * \note The non-null view is expected to be empty and will be populated to hold this
   *  Array's data.
   *
   * \note The Sidre view shape has one dimension.
   *
   * \pre view != nullptr
   * \pre view->isEmpty() == true
   * \pre num_elements >= 1
   *
   * \post view->getNumDimensions() == 1
   * \post view->isEmpty() == false
   * \post size() == num_elements.
   * \post numComponents() == num_components
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <int SFINAE = DIM, typename std::enable_if<SFINAE == 1>::type* = nullptr>
  Array(View* view, axom::IndexType num_elements, axom::IndexType capacity = 0);

  /*!
   * \brief Creates an Array instance of `num_tuples` size, where each
   *  tuple consists of `num_components` values and populates the associated
   *  View.
   *
   * \param [in] view the View that will hold this Array's data.
   * \param [in] num_tuples the number of tuples accounted for in the Array.
   * \param [in] num_components the number of values per tuple. If not
   *  specified defaults to 1.
   * \param [in] capacity the number of tuples to allocate space for.
   *
   * \note The capacity argument is optional. If not specified or if less than
   *  num_tuples, the capacity of the Array will be initialized to
   *  num_tuples * DEFAULT_RESIZE_RATIO.
   *
   * \note The non-null view is expected to be empty and will be populated to hold this
   *  Array's data.
   *
   * \note The Sidre view shape has two dimensions. The first dimension
   *  corresponds to the number of tuples and the second corresponds to
   *  the number of components per tuple.
   *
   * \pre view != nullptr
   * \pre view->isEmpty() == true
   * \pre num_tuples >= 1
   * \pre num_components >= 1
   *
   * \post view->getNumDimensions() == 2
   * \post view->isEmpty() == false
   * \post size() == num_tuples.
   * \post numComponents() == num_components
   * \post getResizeRatio() == DEFAULT_RESIZE_RATIO
   */
  template <int SFINAE = DIM, typename std::enable_if<SFINAE == 2>::type* = nullptr>
  Array(View* view,
        axom::IndexType num_tuples,
        axom::IndexType num_components = 1,
        axom::IndexType capacity = 0);

  /// @}

  /*!
   * Destructor.  Frees the associated buffer unless owned by Sidre.
   */
  virtual ~Array();

  /*!
   * \brief Copy assignment.
   * 
   * Deleted because copies would have to reference the same underlying
   * Sidre buffer, which does not match the underlying axom::Array ownership model
   */
  Array& operator=(const Array&) = delete;

  /*!
   * \brief Move assignment.
   * \param [in] other The Array to move from
   */
  Array& operator=(Array&& other);

  /// \name Array methods to query and set attributes
  /// @{

  /*!
   * \brief Return true iff the external buffer constructor was called.
   */
  virtual bool isExternal() const { return false; }

  /*!
   * \brief Return a pointer to the View that this Array wraps.
   */
  const View* getView() const { return m_view; }

  /// @}

protected:
  /*!
   * \brief Update the number of elements.
   *
   * \param [in] new_num_elements the new number of elements.
   */
  virtual void updateNumElements(axom::IndexType new_num_elements);

  /*!
   * \brief Set the number of elements allocated for the data array.
   *
   * \param [in] capacity the new number of elements to allocate.
   */
  virtual void setCapacity(axom::IndexType new_capacity);

  /*!
   * \brief Reallocates the data array when the size exceeds the capacity.
   *
   * \param [in] new_num_elements the number of elements which exceeds the current
   *  capacity.
   */
  virtual void dynamicRealloc(axom::IndexType new_num_elements);

  /*!
   * \brief Return the TypeID corresponding to T. This function
   *  handles when T is not an enum.
   */
  template <typename U = T>
  static constexpr typename std::enable_if<!std::is_enum<U>::value, TypeID>::type sidreTypeId()
  {
    return detail::SidreTT<U>::id;
  }

  /*!
   * \brief Return the TypeID corresponding to T. This function handles
   *  when T is an enum.
   */
  template <typename U = T>
  static constexpr typename std::enable_if<std::is_enum<U>::value, TypeID>::type sidreTypeId()
  {
    return detail::SidreTT<typename std::underlying_type<U>::type>::id;
  }

  /*!
   * \brief Applies this Array's type and dimensions to the sidre View.
   *
   * Calls m_view->apply(sidreTypeId(), 2, {m_num_tuples, m_num_components}).
   */
  void describeView();

  /*!
   * \brief Given a non-empty View of dimension 2, returns the length
   *  of the given dimension.
   *
   * \param [in] view the View to examine.
   * \param [in] dim the dimension (0 or 1) to return the length of.
   *
   * \pre 0 <= dim <= 1
   */
  axom::IndexType getViewShape(int dim) const;

  /*!
   * \brief Allocates space within the Array's View.
   *
   * \param [in] new_capacity the number of elements to allocate.
   */
  void reallocViewData(IndexType new_capacity);

  View* m_view;
};

/// \brief Helper alias for multi-component arrays
template <typename T>
using MCArray = Array<T, 2>;

//------------------------------------------------------------------------------
//                            Array IMPLEMENTATION
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename T, int DIM>
Array<T, DIM>::Array(Array<T, DIM>&& other)
  : axom::Array<T, DIM>(static_cast<axom::Array<T, DIM>&&>(std::move(other)))
  , m_view(other.m_view)
{
  other.m_view = nullptr;
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
Array<T, DIM>& Array<T, DIM>::operator=(Array<T, DIM>&& other)
{
  axom::Array<T, DIM>::operator=(std::move(other));
  m_view = other.m_view;
  other.m_view = nullptr;
  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
template <int SFINAE, typename std::enable_if<SFINAE == 1>::type*>
Array<T, DIM>::Array(View* view) : axom::Array<T, 1>()
                                 , m_view(view)
{
  SLIC_ERROR_IF(m_view == nullptr, "Provided View cannot be null.");
  SLIC_ERROR_IF(m_view->isEmpty(), "Provided View cannot be empty.");

  axom::StackArray<axom::IndexType, 1> newShape {getViewShape(0)};
  this->setShape(newShape);
  this->m_num_elements = newShape[0];

  axom::IndexType buffer_size = m_view->getBuffer()->getNumElements();
  this->m_capacity = buffer_size;

  SLIC_ERROR_IF(this->m_num_elements < 0,
                "Number of tuples (" << this->m_num_elements << ") "
                                     << "cannot be negative.");

  SLIC_ERROR_IF(this->m_num_elements > this->m_capacity,
                "Number of tuples (" << this->m_num_elements << ") "
                                     << "cannot be greater than the tuple capacity "
                                     << "(" << this->m_capacity << ").");

  TypeID view_type = m_view->getTypeID();
  TypeID T_type = sidreTypeId();
  SLIC_ERROR_IF(view_type != T_type,
                "View data type (" << view_type << ")"
                                   << "differs from this Array type (" << T_type << ").");

  this->m_data = static_cast<T*>(m_view->getVoidPtr());
  SLIC_ERROR_IF(this->m_data == nullptr && this->m_capacity > 0,
                "View returned a null pointer when the capacity "
                  << "is greater than zero.");
}

template <typename T, int DIM>
template <int SFINAE, typename std::enable_if<SFINAE == 2>::type*>
Array<T, DIM>::Array(View* view) : axom::Array<T, 2>()
                                 , m_view(view)
{
  SLIC_ERROR_IF(m_view == nullptr, "Provided View cannot be null.");
  SLIC_ERROR_IF(m_view->isEmpty(), "Provided View cannot be empty.");

  axom::StackArray<axom::IndexType, 2> newShape {getViewShape(0), getViewShape(1)};
  this->setShape(newShape);

  axom::IndexType buffer_size = m_view->getBuffer()->getNumElements();
  SLIC_ERROR_IF(buffer_size % this->m_shape[1] != 0,
                "The buffer size (" << buffer_size << ") "
                                    << "is not a multiple of the number of components "
                                    << "(" << this->m_shape[1] << ").");
  this->m_capacity = buffer_size;
  this->m_num_elements = this->m_shape[0] * this->m_shape[1];

  SLIC_ERROR_IF(this->m_shape[0] < 0,
                "Number of tuples (" << this->m_shape[0] << ") "
                                     << "cannot be negative.");

  SLIC_ERROR_IF(this->m_shape[1] <= 0,
                "Number of components (" << this->m_shape[1] << ") "
                                         << "must be greater than 0.");

  SLIC_ERROR_IF((this->m_shape[0] * this->m_shape[1]) > this->m_capacity,
                "Number of elements (" << this->m_shape[0] * this->m_shape[1] << ") "
                                       << "cannot be greater than the element capacity "
                                       << "(" << this->m_capacity << ").");

  TypeID view_type = m_view->getTypeID();
  TypeID T_type = sidreTypeId();
  SLIC_ERROR_IF(view_type != T_type,
                "View data type (" << view_type << ")"
                                   << "differs from this MCArray type (" << T_type << ").");

  this->m_data = static_cast<T*>(m_view->getVoidPtr());
  SLIC_ERROR_IF(this->m_data == nullptr && this->m_capacity > 0,
                "View returned a null pointer when the capacity "
                  << "is greater than zero.");
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
template <int SFINAE, typename std::enable_if<SFINAE == 1>::type*>
Array<T, DIM>::Array(View* view, axom::IndexType num_elements, axom::IndexType capacity)
  : axom::Array<T, 1>()
  , m_view(view)
{
  SLIC_ERROR_IF(m_view == nullptr, "Provided View cannot be null.");
  SLIC_ERROR_IF(!m_view->isEmpty(), "View must be empty.");
  SLIC_ERROR_IF(num_elements < 0,
                "Number of elements (" << num_elements << ") "
                                       << "cannot be negative.");
  this->m_num_elements = num_elements;
  IndexType real_capacity = capacity;
  if(real_capacity < this->m_num_elements)
  {
    real_capacity = this->m_num_elements;
  }
  reallocViewData(real_capacity);

  SLIC_ERROR_IF(this->m_num_elements > this->m_capacity,
                "Number of elements (" << this->m_num_elements << ") "
                                       << "cannot be greater than the capacity "
                                       << "(" << this->m_capacity << ").");
  // sanity checks
  SLIC_ASSERT(this->m_data != nullptr);
  SLIC_ASSERT(this->m_num_elements >= 0);
}

template <typename T, int DIM>
template <int SFINAE, typename std::enable_if<SFINAE == 2>::type*>
Array<T, DIM>::Array(View* view,
                     axom::IndexType num_tuples,
                     axom::IndexType num_components,
                     axom::IndexType capacity)
  : axom::Array<T, DIM>()
  , m_view(view)
{
  SLIC_ERROR_IF(m_view == nullptr, "Provided View cannot be null.");
  SLIC_ERROR_IF(!m_view->isEmpty(), "View must be empty.");
  SLIC_ERROR_IF(num_tuples < 0,
                "Number of tuples (" << num_tuples << ") "
                                     << "cannot be negative.");
  SLIC_ERROR_IF(num_components <= 0,
                "Components per tuple (" << num_components << ") "
                                         << "must be greater than 0.");
  // FIXME: What we probably want is a Sidre allocator (as opposed to the regular host allocator for example)
  // Would something like that even be possible?
  // FIXME: This isn't quite ideal because we do a "regular" host allocation, delete it, then
  // do an allocation within Sidre
  this->m_shape[0] = num_tuples;
  this->m_shape[1] = num_components;
  this->updateStrides();
  this->m_num_elements = num_tuples * num_components;
  IndexType real_capacity = capacity;
  if(real_capacity < this->m_num_elements)
  {
    real_capacity = this->m_num_elements;
  }
  reallocViewData(real_capacity);

  SLIC_ERROR_IF(this->m_num_elements > this->m_capacity,
                "Number of elements (" << this->m_num_elements << ") "
                                       << "cannot be greater than the element capacity "
                                       << "(" << this->m_capacity << ").");

  // sanity checks
  SLIC_ASSERT(this->m_data != nullptr);
  SLIC_ASSERT(this->m_shape[0] >= 0);
  SLIC_ASSERT(this->m_shape[1] >= 1);
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
Array<T, DIM>::~Array()
{
  m_view = nullptr;
  this->m_data = nullptr;
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline void Array<T, DIM>::updateNumElements(axom::IndexType new_num_elements)
{
  SLIC_ASSERT(new_num_elements >= 0);
  SLIC_ASSERT(new_num_elements <= this->m_capacity);
  this->m_num_elements = new_num_elements;
  describeView();
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline void Array<T, DIM>::setCapacity(axom::IndexType new_capacity)
{
  SLIC_ASSERT(new_capacity >= 0);

  if(new_capacity < this->m_num_elements)
  {
    updateNumElements(new_capacity);
  }

  return reallocViewData(new_capacity);
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline void Array<T, DIM>::dynamicRealloc(axom::IndexType new_num_elements)
{
  SLIC_ERROR_IF(this->m_resize_ratio < 1.0,
                "Resize ratio of " << this->m_resize_ratio << " doesn't support dynamic resizing");

  IndexType new_capacity =
    axom::utilities::max<IndexType>(this->capacity() * this->getResizeRatio() + 0.5,
                                    new_num_elements);
  const IndexType block_size = this->blockSize();
  const IndexType remainder = new_capacity % block_size;
  if(remainder != 0)
  {
    new_capacity += block_size - remainder;
  }
  return reallocViewData(new_capacity);
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline void Array<T, DIM>::describeView()
{
  detail::describeViewImpl(sidreTypeId(), this->shape(), m_view);
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline axom::IndexType Array<T, DIM>::getViewShape(int dim) const
{
  return detail::getViewShapeImpl<DIM>(dim, m_view);
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline void Array<T, DIM>::reallocViewData(IndexType new_capacity)
{
  if(m_view->isEmpty())
  {
    constexpr sidre::TypeID T_type = sidreTypeId();
    m_view->allocate(T_type, new_capacity);
  }
  else
  {
    m_view->reallocate(new_capacity);
  }

  this->m_capacity = new_capacity;
  describeView();
  this->m_data = static_cast<T*>(m_view->getVoidPtr());

  SLIC_ERROR_IF(this->m_data == nullptr && this->m_capacity > 0, "Array reallocation failed.");
}

} /* namespace sidre */
} /* namespace axom */

#endif /* SIDRE_ARRAY_HPP_ */
