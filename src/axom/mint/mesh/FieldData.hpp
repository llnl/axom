// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_FIELDDATA_HPP_
#define MINT_FIELDDATA_HPP_

// Axom includes
#include "axom/core/Macros.hpp"  // for Axom macros
#include "axom/mint/deprecated/MCArray.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre.hpp"
  #include "axom/mint/deprecated/SidreMCArray.hpp"
#endif

// Mint includes
#include "axom/mint/config.hpp"                 // for mint compile time definitions
#include "axom/mint/mesh/Field.hpp"             // for mint::Field definition
#include "axom/mint/mesh/FieldVariable.hpp"     // for mint::FieldVariable
#include "axom/mint/mesh/FieldAssociation.hpp"  // for mint::FieldAssociation

// C/C++ includes
#include <iterator>  // for std::advance
#include <string>    // for std::string
#include <map>       // for std::map
#include <vector>    // for std::vector

namespace axom
{
namespace mint
{
/*!
 * \class FieldData
 *
 * \brief Provides a container for storing fields associated with a specified
 *  mesh topology and methods to create, access and remove fields from the
 *  container.
 *
 *  A FieldData object may store fields associated with a single mesh topology,
 *  e.g., node-, cell- , face- or edge-centered. The FieldData object provides
 *  the means to create, access, modify and remove fields from the container.
 *
 *  Each field in the container is a FieldVariable instance which may
 *  store scalar, vector, tensor fields etc. as explained in more detail in
 *  the FieldVariable documentation.
 *
 *  When Mint is compiled with Sidre support, the FieldData object can be
 *  bound to a particular sidre::Group to store fields in a Sidre hierarchy
 *  according to the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
 *  mesh blueprint </a>.
 *
 *  \note When using Sidre as the back-end storage, the action of creating or
 *   removing fields from the FieldData object is reflected in the associated
 *   Sidre hierarchy by default.
 *
 * \see Field
 * \see FieldVariable
 */
class FieldData
{
public:
  /// \name Constructors
  /// @{

  /*!
   * \brief Default constructor. Disabled.
   */
  FieldData() = delete;

  /*!
   * \brief Creates an empty FieldData instance with the specified mesh
   *  topology association.
   *
   * \param [in] association the mesh topology association, e.g., NODE_CENTERED
   *
   * \note The FieldAssociation enum provides a convenient way to specify
   *  values for the association parameter used in this constructor.
   *
   * \warning ANY_CENTERING, albeit defined in the FieldAssociation enum, it is
   *  used as wild-card and is not a valid value for this constructor.
   *
   * \pre association >= NODE_CENTERED && association < NUM_FIELD_ASSOCIATIONS
   * \post this->getAssociation() == association
   * \post this->empty() == true
   * \post this->hasSidreGroup() == false
   * \post this->getNumFields() == 0
   */
  explicit FieldData(int association);

#ifdef AXOM_MINT_USE_SIDRE
  /*!
   * \brief Constructs a FieldData instance which uses Sidre as the back-end
   *  data-store for storing fields and is bound to the specified group in the
   *  Sidre hierarchy.
   *
   * \param [in] association the mesh topology association, e.g., NODE_CENTERED
   * \param [in] field_group pointer to the Group containing fields.
   *
   * \note The user-supplied Sidre group should conform to the specifications of
   *  outlined in the <a href="http://llnl-conduit.readthedocs.io/en/latest/">
   *  mesh blueprint </a> for storing fields. The key requirements are:
   *  <ul>
   *   <li> Each field is a sub-group under the given group. </li>
   *   <li> Each sub-group consists of the following views:
   *        <ul>
   *          <li> <b> association </b>: encodes the mesh topology </li>
   *          <li> <b> volume_dependent </b>: volume scaling type </li>
   *          <li> <b> values </b>: holds the raw buffer of the field </li>
   *        </ul>
   *   </li>
   *  </ul>
   *
   * \note If the supplied Sidre group is not empty, the resulting FieldData
   *  instance will be populated with the fields that match the specified
   *  association.
   *
   * \pre association >= NODE_CENTERED && association < NUM_FIELD_ASSOCIATIONS
   * \pre field_group != nullptr
   * \post this->getAssociation() == association
   * \post this->hasSidreGroup() == true
   */
  FieldData(int association, sidre::Group* field_group, const std::string& topo);
#endif

  /// @}

  /*!
   * \brief Destructor.
   */
  ~FieldData() { clear(); }

  /// \name Attribute get/set Methods
  /// @{

  /*!
   * \brief Returns the mesh topology association of fields in this container.
   * \return association the field association.
   * \see FieldAssociation.hpp
   */
  inline int getAssociation() const { return m_association; }

  /*!
   * \brief Checks if this FieldData instance is empty.
   * \return status true if empty, else, false.
   */
  inline bool empty() const { return (m_fields.empty()); }

  /*!
   * \brief Returns the number of fields of this FieldData instance.
   * \return N the number of fields in this instance.
   * \post N == 0 \iff this->empty() == true.
   */
  inline int getNumFields() const { return static_cast<int>(m_fields.size()); }

  /*!
   * \brief Checks if the field with the given name exists.
   * \param [in] name the name of the field to check.
   * \return status true if the field exists, else, false.
   */
  inline bool hasField(const std::string& name) const
  {
    return (m_fields.find(name) != m_fields.end());
  }

  /*!
   * \brief Checks if a Sidre Group is associated with this FieldData instance.
   * \return status true if the FieldData is associated with Sidre, else, false.
   */
  inline bool hasSidreGroup() const;

  /*!
   *  \brief Return the resize ratio of this FieldData.
   */
  inline double getResizeRatio() const { return m_resize_ratio; }

  /// @}

  /// \name Methods acting on a single field
  /// @{

  /*!
   * \brief Creates a new field with the given name consisting of the specified
   *  number of tuples and number of components.
   *
   * \param [in] name the user-supplied name to identify this field.
   * \param [in] num_tuples the number of tuples in the field.
   * \param [in] num_components number of components per tuple (optional)
   * \param [in] capacity initial max capacity for the field (optional).
   * \param [in] storeInSidre indicates whether to store the field in the
   *  associated Sidre group of this FieldData instance (optional).
   *
   * \tparam T the underlying data type of the field, e.g., double, int, etc.
   *
   * \note num_components is an optional argument. It defaults to '1' if a value
   *  is not explicitly specified by the caller.
   *
   * \note The 'storeInSidre' boolean argument that indicates whether the field
   *  data should be stored in Sidre. This is an optional argument which
   *  defaults to true and is only relevant when the code is compiled with Sidre
   *  support and the FieldData instance is associated with a Sidre Group. For
   *  persistent data, e.g., state variables etc., Sidre should own the data.
   *  However, for temporary variables it may be desirable to not modify the
   *  Sidre hierarchy.
   *
   * \return ptr pointer to the buffer associated with this field.
   *
   * \pre hasField( name ) == false
   * \pre num_components >= 1
   * \post ptr != nullptr
   *
   */
  template <typename T>
  inline T* createField(const std::string& name,
                        IndexType num_tuples,
                        IndexType num_components = 1,
                        IndexType capacity = USE_DEFAULT,
                        bool storeInSidre = true);

  /*!
   * \brief Creates a new field which will use the supplied external buffer as
   *  storage, consisting of the specified number of tuples and components.
   *
   * \param [in] name the user-supplied name of this field
   * \param [in] data supplied external buffer
   * \param [in] num_tuples the number of tuples in the field.
   * \param [in] num_components the numbere of components per tuple (optional).
   * \param [in] capacity max capacity for the field (optional). If not
   *  specified the capacity defaults to num_tuples.
   *
   * \tparam T the underlying data type of the field, e.g., double, int, etc.
   *
   * \note num_components is an optional argument. It defaults to '1' if a value
   *  is not explicitly specified by the caller.
   *
   * \note The supplied pointer must point to a buffer that is able to hold at
   *  least \f$ num\_tuples \times num\_components \f$
   *
   * \note An external field is not inserted in the Sidre tree hierarchy.
   *
   * \note Once an external field is created, subsequent calls to the FieldData
   *  methods to "resize()" and "reserve()" will fail.
   *
   * \return ptr pointer to the buffer associated with this field.
   *
   * \pre hasField( name ) == false
   * \pre data != nullptr
   * \post ptr != nullptr
   * \post ptr == data
   */
  template <typename T>
  inline T* createField(const std::string& name,
                        T* data,
                        IndexType num_tuples,
                        IndexType num_components = 1,
                        IndexType capacity = USE_DEFAULT);

  /*!
   * \brief Removes the field with the given name.
   *
   * \param [in] name the name of the field to remove.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == true
   *
   */
  void removeField(const std::string& name);

  /*!
   * \brief Removes the ith field from the container.
   *
   * \param [in] i the index of the field to remove from the container
   *
   * \pre i >= 0 && i < getNumFields()
   */
  void removeField(int i);

  /*!
   * \brief Returns the ith field of this FieldData instance.
   *
   * \param [in] i the index of the field in query.
   * \return f pointer to the field in query.
   *
   * \pre i >= 0 && i < this->getNumFields()
   * \post f == nullptr \iff i < 0 || i >= this->getNumberOfFieds()
   */
  /// @{

  inline Field* getField(int i)
  {
    const FieldData* const_this = this;
    return const_cast<Field*>(const_this->getField(i));
  }

  inline const Field* getField(int i) const
  {
    SLIC_ASSERT(i < static_cast<int>(m_fields.size()));
    std::map<std::string, Field*>::const_iterator it = m_fields.begin();
    std::advance(it, i);
    SLIC_ASSERT(it != m_fields.end());
    SLIC_ASSERT(it->second != nullptr);
    return it->second;
  }

  /// @}

  /*!
   * \brief Returns the field with the given name.
   *
   * \param [in] name the name of the field in query.
   * \return f pointer to the field in query.
   *
   * \post f == nullptr \iff this->hasField( name )==false.
   */
  /// @{

  inline Field* getField(const std::string& name)
  {
    const FieldData* const_this = this;
    return const_cast<Field*>(const_this->getField(name));
  }

  inline const Field* getField(const std::string& name) const
  {
    std::map<std::string, Field*>::const_iterator it = m_fields.find(name);
    if(it == m_fields.end())
    {
      return nullptr;
    }
    else
    {
      SLIC_ASSERT(it->second != nullptr);
      return it->second;
    }
  }

  /// @}

  /*!
   * \brief Returns pointer to the buffer of the field with the given name.
   *
   * \param [in] name the name of the field in query.
   * \param [out] num_tuples the number of tuples in the field (optional)
   * \param [out] num_components the number of components per tuple (optional).
   *
   * \return ptr pointer to the buffer of the specified field.
   *
   * \pre hasField( name ) == true
   * \post ptr != nullptr
   */
  /// @{
  template <typename T>
  inline T* getFieldPtr(const std::string& name);

  template <typename T>
  inline T* getFieldPtr(const std::string& name, IndexType& num_tuples);

  template <typename T>
  inline T* getFieldPtr(const std::string& name, IndexType& num_tuples, IndexType& num_components);
  /// @}

  /*!
   * \brief Returns const pointer to the buffer of the field with the given
   * name.
   *
   * \param [in] name the name of the field in query.
   * \param [out] num_tuples the number of tuples in the field (optional)
   * \param [out] num_components the number of components per tuple (optional)
   *
   * \return ptr pointer to the buffer of the specified field.
   *
   * \pre hasField( name ) == true.
   * \post ptr != nullptr
   */
  /// @{
  template <typename T>
  inline const T* getFieldPtr(const std::string& name) const;

  template <typename T>
  inline const T* getFieldPtr(const std::string& name, IndexType& num_tuples) const;

  template <typename T>
  inline const T* getFieldPtr(const std::string& name,
                              IndexType& num_tuples,
                              IndexType& num_components) const;
  /// @}

  /// @}

  /// \name Methods acting on all fields
  /// @{

  /*!
   * \brief Changes the num-tuples of all fields in this FieldData instance.
   *
   * \param [in] newNumTuples the new number of tuples.
   *
   * \warning If the FieldData instance contains a FieldVariable that cannot be
   *  re-sized, e.g., it points to an external buffer, this method will abort
   *  with an error.
   *
   * \see FieldVariable
   */
  void resize(IndexType newNumTuples);

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
  void emplace(IndexType pos, IndexType n_tuples);

  /*!
   * \brief Changes the tuple capacity of all fields in this FieldData instance.
   *
   * \param [in] newCapacity the new max tuple capacity.
   *
   * \warning If the FieldData instance contains a FieldVariable that cannot be
   *  re-sized, e.g., it points to an external buffer, this method will abort
   *  with an error.
   *
   * \see FieldVariable
   */
  void reserve(IndexType newCapacity);

  /*!
   * \brief Shrinks the tuple capacity of all fields in this FieldData instance
   *  to be equal to the actual number of tuples in the field.
   *
   * \see FieldVariable
   */
  void shrink();

  /*!
   *  \brief Set the resize ratio of this FieldData and all associated Fields.
   */
  void setResizeRatio(double ratio);

  /*
   * \brief Return true iff all fields have the given number of tuples, all
   *  non-external fields have the given capacity and the same resize ratio as
   *  this FieldData instance.
   *
   * \param [in] num_tuples the expected number of tuples for each field.
   * \param [in] capacity the expected capacity of the non-external fields.
   */
  bool checkConsistency(IndexType num_tuples, IndexType capacity) const;

  /// @}

private:
  static constexpr int INVALID_FIELD_INDEX = -1;

  /// \name Private helper methods
  /// @{

  /*!
   * \brief Deletes all fields associated with this FieldData instance.
   *
   * \note If the FieldData instance is using Sidre as the back-end data-store,
   *  this method does not remove the associated objects from the Sidre
   *  hierarchy.
   *
   * \post this->empty() == true.
   */
  void clear();

  /*!
   * \brief Returns a string representation of the association name.
   *
   * \return name a string corresponding to the association name.
   *
   * \note By construction, the returned name will never be an empty string.
   *
   * \post name.empty() == false
   */
  inline std::string getAssociationName();

  /*!
   * \brief Removes the field at the given field index.
   *
   * \param [in] i the index of the field to remove.
   *
   * \pre i >= 0 && i < m_fields.size()
   * \pre m_fields[ i ]  != nullptr
   */
  void removeFieldAt(int i);

  /// @}

  int m_association;
  double m_resize_ratio;
  std::map<std::string, Field*> m_fields;

#ifdef AXOM_MINT_USE_SIDRE
  sidre::Group* m_fields_group;
  const std::string m_topology;
#endif

  DISABLE_COPY_AND_ASSIGNMENT(FieldData);
  DISABLE_MOVE_AND_ASSIGNMENT(FieldData);
};

//------------------------------------------------------------------------------
//  IMPLEMENTATION OF TEMPLATE & IN-LINE METHODS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline std::string FieldData::getAssociationName()
{
  // TODO: note, currently edge/face data are not supported by the blueprint

  std::string name = "";
  switch(m_association)
  {
  case NODE_CENTERED:
    name = "vertex";
    break;
  case CELL_CENTERED:
    name = "element";
    break;
  case EDGE_CENTERED:
    name = "edge";
    break;
  case FACE_CENTERED:
    name = "face";
    break;
  default:
    SLIC_ERROR("undefined field association [" << m_association << "]");
  }  // END switch

  return (name);
}

//------------------------------------------------------------------------------
inline bool FieldData::hasSidreGroup() const
{
#ifdef AXOM_MINT_USE_SIDRE
  return (m_fields_group != nullptr);
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
template <typename T>
inline T* FieldData::getFieldPtr(const std::string& name)
{
  IndexType num_tuples = 0;
  IndexType num_components = 0;
  return (getFieldPtr<T>(name, num_tuples, num_components));
}

//------------------------------------------------------------------------------
template <typename T>
inline T* FieldData::getFieldPtr(const std::string& name, IndexType& num_tuples)
{
  IndexType num_components = 0;
  return (getFieldPtr<T>(name, num_tuples, num_components));
}

//------------------------------------------------------------------------------
template <typename T>
inline T* FieldData::getFieldPtr(const std::string& name,
                                 IndexType& num_tuples,
                                 IndexType& num_components)
{
  const FieldData* const_this = this;
  return const_cast<T*>(const_this->getFieldPtr<T>(name, num_tuples, num_components));
}

//------------------------------------------------------------------------------
template <typename T>
inline const T* FieldData::getFieldPtr(const std::string& name) const
{
  IndexType num_tuples = 0;
  IndexType num_components = 0;
  return (getFieldPtr<T>(name, num_tuples, num_components));
}

//------------------------------------------------------------------------------
template <typename T>
inline const T* FieldData::getFieldPtr(const std::string& name, IndexType& num_tuples) const
{
  IndexType num_components = 0;
  return (getFieldPtr<T>(name, num_tuples, num_components));
}

//------------------------------------------------------------------------------
template <typename T>
inline const T* FieldData::getFieldPtr(const std::string& name,
                                       IndexType& num_tuples,
                                       IndexType& num_components) const
{
  SLIC_ERROR_IF(!hasField(name), "field [" << name << "] does not exist!");

  const mint::Field* f = getField(name);
  SLIC_ASSERT(f != nullptr);

  num_tuples = f->getNumTuples();
  num_components = f->getNumComponents();
  SLIC_ASSERT(num_components >= 1);

  return mint::Field::getDataPtr<T>(f);
}

//------------------------------------------------------------------------------
template <typename T>
inline T* FieldData::createField(const std::string& name,
                                 IndexType num_tuples,
                                 IndexType num_components,
                                 IndexType capacity,
                                 bool storeInSidre)
{
  SLIC_ERROR_IF(hasField(name), "Field [" << name << "] already exists!");

  if(capacity == USE_DEFAULT)
  {
    capacity = num_tuples;
  }

  mint::Field* newField = nullptr;
  // create the field on sidre
  if(hasSidreGroup() && storeInSidre)
  {
#ifdef AXOM_MINT_USE_SIDRE

    SLIC_ERROR_IF(m_fields_group->hasGroup(name),
                  "Field [" << name << "] already exists in the Sidre tree!");

    sidre::Group* field = m_fields_group->createGroup(name);
    field->createView("association")->setString(getAssociationName());
    field->createView("volume_dependent")->setString("true");

    // TODO: how should we bind this to the topology?
    field->createView("topology")->setString(m_topology);

    sidre::View* values = field->createView("values");
    newField = new mint::FieldVariable<T>(name, values, num_tuples, num_components, capacity);
#endif
  }  // END if
  else
  {
    newField = new mint::FieldVariable<T>(name, num_tuples, num_components, capacity);
  }  // END else

  SLIC_ASSERT(newField != nullptr);
  newField->setResizeRatio(m_resize_ratio);
  m_fields[name] = newField;

  return (mint::Field::getDataPtr<T>(newField));
}

//------------------------------------------------------------------------------
template <typename T>
inline T* FieldData::createField(const std::string& name,
                                 T* data,
                                 IndexType num_tuples,
                                 IndexType num_components,
                                 IndexType capacity)
{
  SLIC_ERROR_IF(data == nullptr, "supplied buffer is NULL");
  SLIC_ERROR_IF(hasField(name), "Field [" << name << "] already exists!");

  if(capacity == USE_DEFAULT)
  {
    capacity = num_tuples;
  }

  Field* field = new mint::FieldVariable<T>(name, data, num_tuples, num_components, capacity);
  m_fields[name] = field;

  return (mint::Field::getDataPtr<T>(field));
}

} /* namespace mint */
} /* namespace axom */

#endif /* FIELDDATA_HPP_ */
