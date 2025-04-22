// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FIELD_HPP_
#define MINT_FIELD_HPP_

// axom includes
#include "axom/core/Macros.hpp"  // for axom Macros
#include "axom/core/Types.hpp"   // for nullptr

// mint includes
#include "axom/mint/config.hpp"            // for axom::IndexType
#include "axom/mint/mesh/FieldTypes.hpp"   // for the mint::FieldType enum
#include "axom/mint/fem/FEBasisTypes.hpp"  // for the mint::FEBasisType enum

// slic includes
#include "axom/slic/interface/slic.hpp"  // for slic macros

// C/C++ includes
#include <string>  // for std::string

namespace axom
{
namespace mint
{
// Forward declaration
template <typename T>
class FieldVariable;

/*!
 * \class Field
 *
 * \brief Field is an abstract base class which provides the following:
 *
 *  * A type-agnostic wrapper for FieldVariable< T > objects, which are
 *    templated on the underlying field type, e.g., double, int, etc.
 *
 *  * A unified interface to access Fields regardless of type
 *
 *  * The ability to store fields of different types in a FieldData container.
 *
 * \see FieldVariable
 * \see FieldData
 */
class Field
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  Field() = delete;

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor.
   */
  virtual ~Field() { }

  /*!
   * \brief Returns the number of tuples associated with this field.
   * \return N the number of tuples of this field.
   * \post N >= 0
   */
  virtual IndexType getNumTuples() const = 0;

  /*!
   * \brief Returns the number of components associated with this field.
   * \return N the number of components of this field.
   * \post N >= 1
   */
  virtual IndexType getNumComponents() const = 0;

  /*!
   * \brief Returns the total number of tuples this Field instance can hold.
   * \return N the number of tuples this Field instance can hold.
   * \post getCapacity() >= getNumTuples()
   */
  virtual IndexType getCapacity() const = 0;

  /*!
   * \brief Resizes the Field such that it can store the given number of tuples.
   * \param [in] newNumTuples the number of tuples of this Field instance.
   * \note Reallocation is done only if the new size exceeds the capacity.
   */
  virtual void resize(IndexType newNumTuples) = 0;

  /*!
   * \brief Inserts n_tuples with the default value at the given position.
   *
   * \param [in] pos the position of the insert.
   * \param [in] n_tuples the number of tuples to insert.
   *
   * \note The values at pos and above are shifted up and the new tuples
   *  have the default values.
   *
   * \see FieldVariable
   */
  virtual void emplace(IndexType pos, IndexType n_tuples) = 0;

  /*!
   * \brief Increase the Field capacity to hold the given number of tuples.
   * \param [in] newCapacity number of tuples to reserve memory for.
   * \note if newCapacity < getCapacity() this method returns immediately.
   */
  virtual void reserve(IndexType newCapacity) = 0;

  /*!
   * \brief Shrinks the field capacity to be equal to the number of tuples.
   * \post getCapacity()==getNumTuple()
   */
  virtual void shrink() = 0;

  /*!
   * \brief Return the resize ratio of this field.
   */
  virtual double getResizeRatio() const = 0;

  /*!
   * \brief Set the resize ratio of this field.
   * \param [in] ratio the new resize ratio.
   * \post getResizeRatio() == ratio
   */
  virtual void setResizeRatio(double ratio) = 0;

  /*!
   * \brief Return true iff the field is stored in an external buffer.
   */
  virtual bool isExternal() const = 0;

  /*!
   * \brief Return true iff the field is stored in sidre.
   */
  virtual bool isInSidre() const = 0;

  /// @}

  /// \name Attribute get/set Methods
  /// @{

  /*!
   * \brief Returns the type of the field.
   * \return t the type of the field.
   * \post t < NUMBER_OF_FIELD_TYPES
   * \see FieldTypes
   */
  int getType() const { return m_type; }

  /*!
   * \brief Returns the name of the field.
   * \return name the name of the field.
   */
  const std::string& getName() const { return m_name; }

  /*!
   * \brief Returns the basis associated with this field.
   * \return basis the basis associated with this field.
   * \pre fe_basis >= 0 && fe_basis < MINT_NUM_BASIS_TYPES
   */
  int getBasis() const { return m_basis; };

  /*!
   * \brief Sets the finite element basis associated with this field.
   * \param [in] fe_basis the finite element basis
   *
   * \note A field is not associated with basis type by default. The caller
   *  must explicitly associate a field with a basis by calling `setBasis()`
   *
   * \pre fe_basis >= 0 && fe_basis < MINT_NUM_BASIS_TYPES
   * \see FEBasisTypes for a list of supported basis.
   */
  void setBasis(int fe_basis)
  {
    SLIC_ASSERT(fe_basis >= 0 && fe_basis < MINT_NUM_BASIS_TYPES);
    m_basis = fe_basis;
  }

  /// @}

  /// \name Static Methods
  /// @{

  /*!
   * \brief Returns raw pointer to the field data.
   * \return dataPtr pointer to the field data.
   * \tparam T the field data type, e.g., double, int, etc.
   * \pre field_traits< T >::type() == field->getType()
   * \post dataPtr != nullptr
   */
  /// @{
  template <typename T>
  static inline T* getDataPtr(Field* field);

  template <typename T>
  static inline const T* getDataPtr(const Field* field);
  /// @}

  /// @}

protected:
  /*!
   * \brief Custom constructor to call from a derived class.
   *
   * \param [in] name the name associated with this field instance.
   * \param [in] type the field type.
   *
   * \pre name.empty() == false
   * \pre type >= 0 && type < NUMBER_OF_FIELD_TYPES
   *
   * \see FieldType for a list of supported field types.
   */
  Field(const std::string& name, int type)
    : m_name(name)
    , m_type(type)
    , m_basis(MINT_UNDEFINED_BASIS)
  {
    SLIC_ERROR_IF(m_name.empty(), "Supplied Field name is empty!");
    SLIC_ERROR_IF(m_type == UNDEFINED_FIELD_TYPE,
                  "Supplied field type doesn't map to a supported type!");
  }

  std::string m_name; /*!< the name of the field  */
  int m_type;         /*!< the field type */
  int m_basis;        /*!< associated finite element basis with the field. */

private:
  DISABLE_COPY_AND_ASSIGNMENT(Field);
  DISABLE_MOVE_AND_ASSIGNMENT(Field);
};

//------------------------------------------------------------------------------
//             Field IMPLEMENTATION
//------------------------------------------------------------------------------

template <typename T>
inline T* Field::getDataPtr(Field* field)
{
  SLIC_ASSERT(field != nullptr);

  // check type
  constexpr int type = field_traits<T>::type();
  int ftype = field->getType();
  SLIC_ERROR_IF((type == UNDEFINED_FIELD_TYPE),
                "Template argument to Field::getDataPtr() doesn't map to a supported type");
  SLIC_ERROR_IF((type != ftype),
                "Template argument to Field::getDataPtr() doesn't match the field type");

  FieldVariable<T>* f = static_cast<FieldVariable<T>*>(field);
  SLIC_ASSERT(f != nullptr);
  return (f->getFieldVariablePtr());
}

//------------------------------------------------------------------------------
template <typename T>
inline const T* Field::getDataPtr(const Field* field)
{
  SLIC_ASSERT(field != nullptr);

  // check type
  int type = field_traits<T>::type();
  int ftype = field->getType();
  SLIC_ERROR_IF((type == UNDEFINED_FIELD_TYPE),
                "Template argument to Field::getDataPtr() doesn't map to a supported type");
  SLIC_ERROR_IF((type != ftype),
                "Template argument to Field::getDataPtr() doesn't match the field type");

  const FieldVariable<T>* f = static_cast<const FieldVariable<T>*>(field);
  SLIC_ASSERT(f != nullptr);
  return (f->getFieldVariablePtr());
}

} /* namespace mint */
} /* namespace axom */

#endif /* FIELD_HPP_ */
