// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FIELDVARIABLE_HPP_
#define MINT_FIELDVARIABLE_HPP_

#include "axom/mint/mesh/Field.hpp"

// axom includes
#include "axom/core/Macros.hpp"  // for axom Macros
#include "axom/core/Types.hpp"   // for axom types
#include "axom/mint/deprecated/MCArray.hpp"

#include "axom/mint/config.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
  #include "axom/mint/deprecated/SidreMCArray.hpp"
#endif

#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <string>  // for C++ string

namespace axom
{
namespace mint
{
/*!
 * \class FieldVariable
 *
 * \brief Provides the ability to store, access and modify a mesh field.
 *
 *  FieldVariable is a concrete implementation of mint::Field. The FieldVariable
 *  class provides the ability to store, access and modify a field. Most often
 *  a field is associated with entities on a mesh. For example, temperature
 *  or velocity at the nodes of a mesh, or, pressure and mass evaluated at cell
 *  centers.
 *
 *  A FieldVariable may be used to represent a scalar quantity, a vector
 *  quantity or a tensor field. The number of tuples of the field corresponds to
 *  the number of corresponding mesh entities at which the field is evaluated,
 *  e.g., the number of nodes in the mesh for a node-centered field. The number
 *  of components of the FieldVariable may be used to indicate the number of
 *  components for a vector or tensor field. For example, a 3D velocity field
 *  may have \f$ 3 \f$ components for the \f$ x \f$, \f$ y \f$, and \f$ z \f$
 *  velocity components. Similarly, a \f$ 3 \times 3 \f$ tensor field, may be
 *  stored using \f$ 9 \f$ components that correspond to the tensor
 *  components.
 *
 *  The FieldVariable class provides support for a variety of field types
 *  through the use of  templates. A list of supported field types is given in
 *  FieldTypes.hpp
 *
 *  A FieldVariable object may be constructed using (a) native storage, (b)
 *  external storage, or, (c) from Sidre:
 *
 *  * <b> Native Storage </b> <br />
 *
 *    When using native storage, the FieldVariable object owns all associated
 *    memory. The storage can dynamically grow as needed, e.g., when doing
 *    refinement. When the object is deleted, all memory associated with
 *    the given instance is returned to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    A FieldVariable may also be constructed by pointing it to an external,
 *    user-supplied buffer. In this case, the FieldVariable does not own the
 *    memory. Consequently, the FieldVariable cannot be resized, i.e., the
 *    number of tuples and number of components stays fixed for the life-time
 *    of the object.
 *
 *    \warning All calls to shrink(), resize() and reserve() will fail if a
 *     FieldVariable object is constructed from an external buffer.
 *
 *  * <b> Sidre </b> <br />
 *
 *    When Sidre is enabled, a FieldVariable may be constructed from a
 *    corresponding sidre::View. In this case, Sidre has ownership of
 *    the data, but, storage can dynamically grow as needed. All memory
 *    management operations are, therefore, delegated to Sidre.
 *
 *    When the FieldVariable object is deleted, the, associated data in the
 *    Sidre::View remains persistent in Sidre.
 *
 * \tparam T the type of field, e.g., double, int, etc.
 *
 * \see mint::Field
 * \see mint::FieldData
 * \see sidre::View
 */
template <typename T>
class FieldVariable : public Field
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  FieldVariable() = delete;

  /// \name Native Storage Field Variable Constructor
  /// @{

  /*!
   * \brief Creates a FieldVariable instance with the given name,
   *  size and number of components per tuple.
   *
   * \param [in] name the name of this FieldVariable.
   * \param [in] num_tuples the number of tuples.
   * \param [in] num_components num components per tuple (optional).
   * \param [in] capacity number of tuples to allocate space for (optional).
   *
   * \note num_components is set to 1 if not specified.
   * \note If capacity is not explicitly specified, an internal default is
   *  used instead.
   *
   * \pre name.empty() == false
   * \pre num_tuples >= 0
   * \pre num_components >= 1
   */
  FieldVariable(const std::string& name,
                IndexType num_tuples,
                IndexType num_components = 1,
                IndexType capacity = USE_DEFAULT);

  /// @}

  /// \name External Storage FieldVariable Constructors
  /// @{

  /*!
   * \brief Creates a FieldVariable that points to a supplied external buffer.
   *
   * \param [in] name the name of this FieldVariable.
   * \param [in] data pointer to the external buffer.
   * \param [in] num_tuples the number of tuples to allocate.
   * \param [in] num_components the number of components per tuple (optional).
   *
   * \note If num_components is set to 1 if not specified.
   *
   * \note The supplied pointer must point to a buffer that is sufficiently
   *  allocated to hold \f$ num\_tuples \times num\_components \f$ items.
   *
   * \warning When constructing a FieldVariable instance from an external
   *  buffer, all calls to shrink(), resize() and reserve() will fail.
   *
   * \pre name.empty() == false
   * \pre num_tuples >= 0
   * \pre num_components >= 1
   * \pre data != nullptr
   */
  FieldVariable(const std::string& name,
                T* data,
                IndexType num_tuples,
                IndexType num_components = 1,
                IndexType capacity = USE_DEFAULT);

  /// @}

  /// \name Sidre FieldVariable Constructors
  /// @{

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Creates a FieldVariable instance from a sidre::View that has data.
   *
   * \param [in] name the name of this FieldVariable.
   * \param [in] field_view the sidre::View that holds the field data.
   *
   * \pre name.empty() == false
   * \pre field_view != nullptr
   * \pre field_view->isEmpty() == false
   *
   * \post this->getNumTuples() >= 0
   * \post this->getNumComponents() >= 1
   *
   * \see sidre::View
   */
  FieldVariable(const std::string& name, sidre::View* field_view);

  /*!
   * \brief Creates a FieldVariable with the given name, number of tuples and
   *  number of components per tuple, on the supplied, empty, sidre::View.
   *
   * \param [in] field_view pointer to the sidre::View that will hold the field
   * \param [in] name the name associated with this field instance.
   * \param [in] num_tuples the number of tuples of the field.
   * \param [in] num_components the number of components per tuple (optional).
   * \param [in] capacity number of tuples to allocate space for (optional).
   *
   * \note num_components is set to 1 if not specified.
   * \note If capacity is not explicitly specified, an internal default is
   *  used instead.
   *
   * \note The supplied view is expected be be empty and will be populated to
   *  hold the data associated with this FieldVairable.
   *
   * \pre field_view != nullptr
   * \pre field_view->isEmpty() == true
   *
   * \pre num_tuples >= 0
   * \pre num_components >= 1
   *
   * \see sidre::View
   */
  FieldVariable(const std::string& name,
                sidre::View* field_view,
                IndexType num_tuples,
                IndexType num_components = 1,
                IndexType capacity = USE_DEFAULT);
#endif

  /// @}

  /// \name Virtual Methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~FieldVariable() { delete m_field; }

  /*!
   * \brief Returns the number of tuples of this FieldVariable instance.
   * \return N the number of tuples of this FieldVariable.
   * \post N >= 0
   * \see Field::getNumTuples()
   */
  virtual IndexType getNumTuples() const final override { return m_field->size(); }

  /*!
   * \brief Return the number of components per tuple.
   * \return N the number of components per tuple
   * \post N >= 1
   * \see Field::getNumComponents()
   */
  virtual IndexType getNumComponents() const final override { return m_field->numComponents(); };

  /*!
   * \brief Returns the total number of tuples this instance can hold.
   * \return N the capacity of this FieldVariable instance.
   * \post N >= this->getNumTuples()
   * \see Field::getCapacity()
   */
  virtual IndexType getCapacity() const final override { return m_field->capacity(); };

  /*!
   * \brief Resizes the Field such that it can store the given number of tuples.
   * \param [in] newNumTuples the number of tuples of this Field instance.
   * \note Reallocation is done only if the new size exceeds the capacity.
   * \see Field::resize()
   */
  virtual void resize(IndexType newNumTuples) final override { m_field->resize(newNumTuples); }

  /*!
   * \brief Inserts n_tuples with the default value at the given position.
   *
   * \param [in] pos the position of the insert.
   * \param [in] n_tuples the number of tuples to insert.
   *
   * \note The values at pos and above are shifted up and the new tuples
   *  have the default values.
   */
  virtual void emplace(IndexType pos, IndexType num_tuples) final override
  {
    m_field->emplace(num_tuples, pos);
  }

  /*!
   * \brief Increase the Field capacity to hold the given number of tuples.
   * \param [in] newCapacity number of tuples to reserve memory for.
   * \note if newCapacity < getCapacity() this method returns immediately.
   * \see Field::reserve()
   */
  virtual void reserve(IndexType newCapacity) final override { m_field->reserve(newCapacity); }

  /*!
   * \brief Shrinks the field capacity to be equal to the number of tuples.
   * \post getCapacity()==getNumTuple()
   * \see Field::shrink()
   */
  virtual void shrink() final override { m_field->shrink(); }

  /*!
   * \brief Return the resize ratio of this field.
   */
  virtual double getResizeRatio() const final override { return m_field->getResizeRatio(); }

  /*!
   * \brief Set the resize ratio of this field.
   * \param [in] ratio the new resize ratio.
   * \post getResizeRatio() == ratio
   */
  virtual void setResizeRatio(double ratio) final override { m_field->setResizeRatio(ratio); }

  /*!
   * \brief Return true iff the field is stored in an external buffer.
   */
  virtual bool isExternal() const final override { return m_field->isExternal(); }

  /*!
   * \brief Return true iff the field is stored in sidre.
   */
  virtual bool isInSidre() const final override { return m_field->isInSidre(); }

  /// @}

  /// \name Data Access Methods
  /// @{

  /*!
   * \brief Returns pointer to the FieldVariable data.
   * \return ptr pointer to the data associated with this field variable.
   * \post ptr != nullptr
   */
  /// @{

  inline T* getFieldVariablePtr() { return m_field->getData(); }

  inline const T* getFieldVariablePtr() const { return m_field->getData(); }

  /// @}

  /// @}

private:
  axom::deprecated::MCArray<T>* m_field;

  DISABLE_COPY_AND_ASSIGNMENT(FieldVariable);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldVariable);
};

//------------------------------------------------------------------------------
//                FieldVariable IMPLEMENTATION
//------------------------------------------------------------------------------

template <typename T>
FieldVariable<T>::FieldVariable(const std::string& name,
                                IndexType num_tuples,
                                IndexType num_components,
                                IndexType capacity)
  : Field(name, field_traits<T>::type())
{
  m_field = new axom::deprecated::MCArray<T>(num_tuples, num_components, capacity);
  SLIC_ASSERT(m_field != nullptr);
  SLIC_ERROR_IF(m_type == UNDEFINED_FIELD_TYPE, "Undefined field type!");
}

//------------------------------------------------------------------------------
template <typename T>
FieldVariable<T>::FieldVariable(const std::string& name,
                                T* data,
                                IndexType num_tuples,
                                IndexType num_components,
                                IndexType capacity)
  : Field(name, field_traits<T>::type())
{
  m_field = new axom::deprecated::MCArray<T>(data, num_tuples, num_components, capacity);
  SLIC_ASSERT(m_field != nullptr);
  SLIC_ASSERT(m_field->isExternal() == true);
  SLIC_ERROR_IF(m_type == UNDEFINED_FIELD_TYPE, "Undefined field type!");
}

#ifdef AXOM_MINT_USE_SIDRE

//------------------------------------------------------------------------------
template <typename T>
FieldVariable<T>::FieldVariable(const std::string& name, sidre::View* field_view)
  : Field(name, field_traits<T>::type())
{
  m_field = new sidre::deprecated::MCArray<T>(field_view);
  SLIC_ASSERT(m_field != nullptr);
  SLIC_ERROR_IF(m_type == UNDEFINED_FIELD_TYPE, "Undefined field type!");
}

//------------------------------------------------------------------------------
template <typename T>
FieldVariable<T>::FieldVariable(const std::string& name,
                                sidre::View* field_view,
                                IndexType num_tuples,
                                IndexType num_components,
                                IndexType capacity)
  : Field(name, field_traits<T>::type())
{
  m_field = new sidre::deprecated::MCArray<T>(field_view, num_tuples, num_components, capacity);
  SLIC_ASSERT(m_field != nullptr);
  SLIC_ERROR_IF(m_type == UNDEFINED_FIELD_TYPE, "Undefined field type!");
}

#endif

} /* namespace mint */
} /* namespace axom */

#endif /* FIELDVARIABLE_HPP_ */
