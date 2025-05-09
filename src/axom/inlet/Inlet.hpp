// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.hpp
 *
 * \brief This file contains the class definition of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#ifndef INLET_INLET_HPP
#define INLET_INLET_HPP

#include <memory>
#include <string>
#include <vector>
#include <functional>

#include "axom/inlet/Container.hpp"
#include "axom/inlet/Field.hpp"
#include "axom/inlet/Proxy.hpp"
#include "axom/inlet/Reader.hpp"

#include "axom/sidre.hpp"

#include "axom/inlet/Writer.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class Inlet
 *
 * \brief This class is the main access point for all Inlet operations from
 *        from defining the schema of the users input file to getting the values
 *        out of the Sidre DataStore.
 *
 * \see Container Field
 *******************************************************************************
 */
class Inlet
{
public:
  /*!
   *****************************************************************************
   * \brief Constructor for the Inlet class.
   *
   * Creates an Inlet class that can then be used with the given Reader and will
   * store data under the given Sidre Group.
   *
   * \param [in] reader Unique (owning) pointer to the input file Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] docEnabled Boolean indicating whether documentation generation
   * is enabled. This also toggles the storing of documentation-specific information.
   * \param [in] reconstruct Whether or not to attempt to reconstruct child Containers
   * and Fields from the data in the sidre Group
   *****************************************************************************
   */
  Inlet(std::unique_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup,
        bool docEnabled = true,
        bool reconstruct = false)
    : m_reader(std::move(reader))
    , m_sidreRootGroup(sidreRootGroup)
    , m_globalContainer("", "", *m_reader, m_sidreRootGroup, m_unexpectedNames, docEnabled, reconstruct)
    , m_docEnabled(docEnabled)
  {
    m_unexpectedNames = m_reader->getAllNames();
  }

  /// \overload
  Inlet(std::unique_ptr<Reader> reader, bool docEnabled = true, bool reconstruct = false)
    : m_reader(std::move(reader))
    , m_datastore(new sidre::DataStore)
    , m_sidreRootGroup(m_datastore->getRoot())
    , m_globalContainer("", "", *m_reader, m_sidreRootGroup, m_unexpectedNames, docEnabled, reconstruct)
    , m_docEnabled(docEnabled)
  {
    m_unexpectedNames = m_reader->getAllNames();
  }

  // Inlet objects must be move only - delete the implicit shallow copy constructor
  Inlet(const Inlet&) = delete;
  Inlet(Inlet&&) = default;

  virtual ~Inlet() = default;

  /*!
   *****************************************************************************
   * \brief Returns the reference to the Reader class.
   *
   * Provides access to the Reader class that is used to access the input file.
   *
   * \return Reference to this instances' Reader class
   *****************************************************************************
   */
  Reader& reader() { return *m_reader; };

  /*!
   *****************************************************************************
   * \brief Returns pointer to the root Sidre Group class for all of Inlet.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for all of Inlet.
   *
   * \return Pointer to the root Sidre Group for Inlet
   *****************************************************************************
   */
  axom::sidre::Group* sidreGroup() { return m_sidreRootGroup; };

  //
  // Functions that define the input file schema
  //

  /*!
   *****************************************************************************
   * \brief Add a structure to the input file schema.
   *
   * Adds a structure/record to the input file schema. Structures can contain
   * fields and/or substructures.  By default, it is not required unless marked with
   * Container::isRequired(). This creates the Sidre Group class with the given name and
   * stores the given description.
   *
   * \param [in] name Name of the struct expected in the input file
   * \param [in] description Description of the struct
   *
   * \return Reference to the created struct, as a Container
   *****************************************************************************
   */
  Container& addStruct(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a Boolean Field to the input file schema.
   *
   * Adds a Boolean Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input file
   * \param [in] description Description of the Field
   *
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addBool(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a Double Field to the input file schema.
   *
   * Adds a Double Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input file
   * \param [in] description Description of the Field
   *
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addDouble(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a Integer Field to the input file schema.
   *
   * Adds a Integer Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input file
   * \param [in] description Description of the Field
   *
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addInt(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a String Field to the input file schema.
   *
   * Adds a String Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Container expected in the input file
   * \param [in] description Description of the Container
   *
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addString(const std::string& name, const std::string& description = "");

  //
  // Functions that get the values out of the datastore
  //

  /*!
   *******************************************************************************
   * \brief Gets a value of arbitrary type out of the datastore
   * 
   * Retrieves a value of primitive or user-defined type.
   * 
   * \param [in] name The name of the subcontainer representing the root of the object
   * \return The retrieved value
   * \tparam The type to retrieve
   * \pre Requires a specialization of FromInlet<T> for user-defined types
   * \note This function does not indicate failure in a way that can be handled
   * by a program - if an object of requested type does not exist at the specified
   * location, the program will terminate
   *******************************************************************************
   */
  template <typename T>
  T get(const std::string& name) const
  {
    return m_globalContainer.get<T>(name);
  }

  /*!
   *****************************************************************************
   * \brief Return whether a subobject with the given name is present in 
   * the datastore.
   *
   * \see Container::contains
   *****************************************************************************
   */
  bool contains(const std::string& name) const { return m_globalContainer.contains(name); }

  /*!
   *****************************************************************************
   * \brief Return whether a subobject with the given name was provided by the user
   *
   * \see Container::isUserProvided
   *****************************************************************************
   */
  bool isUserProvided(const std::string& name) const
  {
    return m_globalContainer.isUserProvided(name);
  }

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the datastore.
   * 
   * \see Container::operator[]
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) const { return m_globalContainer[name]; }

  /*!
   *****************************************************************************
   * \brief Writes input file documentation.
   *
   * This runs the calling Inlet object through the \a writer.
   * 
   * \param [in] writer The writer object to use
   *
   *****************************************************************************
   */
  void write(Writer&& writer);

  /*!
   *****************************************************************************
   * \brief Verifies the contents of the sidreGroup according to Inlet 
   * requirements.
   * \param [in] errors An optional vector of errors to append to in the case
   * of verification failure
   * 
   * Ownership is not taken of @a errors, the raw pointer is only used for its
   * optional reference semantics, as opposed to something like
   * std::optional<std::reference_wrapper<T>>
   *
   * This recursively checks the correctness of each Field and Container in the Sidre
   * Group: ensuring that required Fields are specified, each Field's value 
   * and default value are within the specified range or are equal to a valid 
   * value, and types are consistent. Also ensures that the registered verification
   * functions hold true.
   * 
   * \return true if contents are correct and false if not.
   *
   *****************************************************************************
   */
  bool verify(std::vector<VerificationError>* errors = nullptr) const;

  /*!
   *****************************************************************************
   * \return The global Container.
   *****************************************************************************
   */
  Container& getGlobalContainer() { return m_globalContainer; }

  /*!
   *****************************************************************************
   * \brief Add an array of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addBoolArray(const std::string& name, const std::string& description = "")
  {
    return m_globalContainer.addBoolArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of Integer Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addIntArray(const std::string& name, const std::string& description = "")
  {
    return m_globalContainer.addIntArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of Double Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addDoubleArray(const std::string& name, const std::string& description = "")
  {
    return m_globalContainer.addDoubleArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of String Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addStringArray(const std::string& name, const std::string& description = "")
  {
    return m_globalContainer.addStringArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of user-defined type to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Container& addStructArray(const std::string& name, const std::string& description = "")
  {
    return m_globalContainer.addStructArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Get a function from the input deck
   *
   * \param [in] name        Name of the function
   * \param [in] ret_type    The return type of the function
   * \param [in] arg_types   The argument types of the function
   * \param [in] description Description of the Field
   *
   * \return Reference to the created Function
   *****************************************************************************
   */
  Verifiable<Function>& addFunction(const std::string& name,
                                    const FunctionTag ret_type,
                                    const std::vector<FunctionTag>& arg_types,
                                    const std::string& description = "")
  {
    return m_globalContainer.addFunction(name, ret_type, arg_types, description);
  }
  /*!
   *****************************************************************************
   * \brief Add a dictionary of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Container>& addBoolDictionary(const std::string& name,
                                           const std::string& description = "")
  {
    return m_globalContainer.addBoolDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add a dictionary of Integer Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Container>& addIntDictionary(const std::string& name,
                                          const std::string& description = "")
  {
    return m_globalContainer.addIntDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add a dictionary of Double Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Container>& addDoubleDictionary(const std::string& name,
                                             const std::string& description = "")
  {
    return m_globalContainer.addDoubleDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add a dictionary of String Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Container>& addStringDictionary(const std::string& name,
                                             const std::string& description = "")
  {
    return m_globalContainer.addStringDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an dictionary of user-defined type to the input file schema.
   *
   * \param [in] name Name of the dictionary
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Container& addStructDictionary(const std::string& name, const std::string& description = "")
  {
    return m_globalContainer.addStructDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Returns the global list of unexpected names, i.e., entries
   * in the input file that were not added via an add* call
   *****************************************************************************
   */
  const std::vector<std::string>& unexpectedNames() const { return m_unexpectedNames; }

  // TODO add update value functions
private:
  std::unique_ptr<Reader> m_reader;
  // Used only in the case where the user does not provide an initial root group
  std::unique_ptr<sidre::DataStore> m_datastore;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  Container m_globalContainer;
  bool m_docEnabled;
  std::vector<std::string> m_unexpectedNames;
};

}  // end namespace inlet
}  // end namespace axom

#endif
