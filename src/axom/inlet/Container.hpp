// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Container.hpp
 *
 * \brief This file contains the class definition of Inlet's Container class.
 *******************************************************************************
 */

#ifndef INLET_CONTAINER_HPP
#define INLET_CONTAINER_HPP

#include <memory>
#include <string>
#include <functional>
#include <unordered_map>
#include <tuple>
#include <type_traits>

#include "axom/fmt.hpp"

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Function.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/VariantKey.hpp"
#include "axom/inlet/Verifiable.hpp"

#include "axom/sidre.hpp"

/*!
 *******************************************************************************
 * \brief Prototype for user-defined types wishing to define a by-value read
 * from inlet with axom::inlet::Container::get(const std::string&)
 * 
 * This is the only way of reading in a non-default-constructible type from
 * Inlet.
 * 
 * \see axom::inlet::Container::get(const std::string&)
 *******************************************************************************
 */
template <typename T>
struct FromInlet
{
  // "Poison" the base implementation so it cannot be used
  // and so it can be checked in a SFINAE context
  FromInlet() = delete;
};

namespace axom
{
namespace inlet
{
// Forward declaration for the traits

class Container;

namespace detail
{
/*!
 *******************************************************************************
 * \class is_inlet_primitive
 *
 * \brief A type trait for checking if a type is isomorphic to an Inlet primitive
 * \tparam T The type to check
 * \note An Inlet primitive is any of the following C++ types: int, double, bool,
 * std::string.
 * \see InletType
 *******************************************************************************
 */
template <typename T>
struct is_inlet_primitive
{
  using BaseType = typename std::decay<T>::type;
  static constexpr bool value = std::is_same<BaseType, bool>::value ||
    std::is_same<BaseType, int>::value || std::is_same<BaseType, double>::value ||
    std::is_same<BaseType, std::string>::value;
};

/*!
 *******************************************************************************
 * \class is_inlet_primitive
 *
 * \brief A type trait for checking if a type is isomorphic to an array of Inlet
 * primitives
 * \tparam T The type to check
 *******************************************************************************
 */
template <typename T>
struct is_inlet_primitive_array : std::false_type
{ };

// If it's an unordered map whose value type is a Inlet primitive,
// assume that it's an array
template <typename T>
struct is_inlet_primitive_array<std::unordered_map<int, T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

template <typename T>
struct is_inlet_array : std::false_type
{ };

template <typename T>
struct is_inlet_array<std::unordered_map<int, T>> : std::true_type
{ };

// Determines whether a type is a std::function
template <typename T>
struct is_std_function : std::false_type
{ };

template <typename T>
struct is_std_function<std::function<T>> : std::true_type
{ };

// Extracts the signature of a std::function
template <typename FuncType>
struct std_function_signature;

template <typename FuncType>
struct std_function_signature<std::function<FuncType>>
{
  using type = FuncType;
};

template <typename T>
struct is_inlet_dict : std::false_type
{ };

template <typename T>
struct is_inlet_dict<std::unordered_map<std::string, T>> : std::true_type
{ };

template <typename T>
struct is_inlet_dict<std::unordered_map<VariantKey, T>> : std::true_type
{ };

template <typename T>
struct is_inlet_primitive_dict : std::false_type
{ };

template <typename T>
struct is_inlet_primitive_dict<std::unordered_map<std::string, T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

template <typename T>
struct is_inlet_primitive_dict<std::unordered_map<VariantKey, T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

template <typename T>
struct is_std_vector : std::false_type
{ };

template <typename T>
struct is_std_vector<std::vector<T>> : std::true_type
{ };

template <typename T>
struct is_primitive_std_vector : std::false_type
{ };

template <typename T>
struct is_primitive_std_vector<std::vector<T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

/*!
 *******************************************************************************
 * \class has_FromInlet_specialization
 *
 * \brief A type trait for checking if a type has specialized FromInlet
 * with the required T operator()(axom::inlet::Container&)
 * \tparam T The type to check
 *******************************************************************************
 */
template <typename T, typename SFINAE = void>
struct has_FromInlet_specialization : std::false_type
{ };

template <typename T>
struct has_FromInlet_specialization<
  T,
  typename std::enable_if<
    std::is_same<T, decltype(std::declval<FromInlet<T>&>()(std::declval<const Container&>()))>::value>::type>
  : std::true_type
{ };

/*!
 *******************************************************************************
 * \brief An overloaded utility function for converting a type to a string
 * 
 * \note Needed as std::to_string doesn't implement an identity overload
 * 
 * \param [in] idx The index to convert to string
 *******************************************************************************
 */
inline std::string indexToString(const std::string& idx) { return idx; }
/// \overload
inline std::string indexToString(const int idx) { return std::to_string(idx); }
/// \overload
inline std::string indexToString(const VariantKey& idx)
{
  return idx.type() == InletType::String ? static_cast<std::string>(idx)
                                         : indexToString(static_cast<int>(idx));
}

/*!
 *******************************************************************************
 * \brief An templated utility function for converting an index to the desired type
 * 
 * \param [in] idx The index to convert
 * \tparam From The type of the idx parameter (converting from)
 * \tparam Result The type to convert to
 *******************************************************************************
 */
template <typename Result, typename From>
inline Result toIndex(const From& idx)
{
  // By default, try a static cast
  return static_cast<Result>(idx);
}
/// \overload
template <>
inline int toIndex(const std::string& idx)
{
  SLIC_ERROR_IF(!conduit::utils::string_is_integer(idx),
                fmt::format("[Inlet] Expected an integer, got: {0}", idx));
  return conduit::utils::string_to_value<int>(idx);
}

/*!
 *******************************************************************************
 * \brief Determines whether a variant key is convertible to another type
 * 
 * \tparam Key The type to check the validity of the conversion to
 * 
 * \note That is, returns true if the key holds an integer and Key is int, etc
 *******************************************************************************
 */
template <typename Key>
bool matchesKeyType(const VariantKey& key)
{
  // A VariantKey matches everything
  if(std::is_same<Key, VariantKey>::value)
  {
    return true;
  }
  else if(std::is_same<Key, int>::value && key.type() == InletType::Integer)
  {
    return true;
  }
  else if(std::is_same<Key, std::string>::value && key.type() == InletType::String)
  {
    return true;
  }
  return false;
}

/*!
 *****************************************************************************
 * \brief This is an internal utility intended to be used with arrays/dicts of 
 * user-defined types that returns the indices as strings - integer indices
 * will be converted to strings
 * 
 * \param [in] container The container to retrieve indices from
 * \param [in] trimAbsolute Whether to only return the "basename" if the path
 * is absolute, e.g., an absolute path foo/0/bar will be trimmed to "bar"
 *****************************************************************************
 */
std::vector<VariantKey> collectionIndices(const Container& container, bool trimAbsolute = true);

/*!
 *****************************************************************************
 * \brief This is an internal utility intended to be used with arrays of 
 * user-defined types that returns the a list of pairs, each of which contain
 * an index (a number) and a fully qualified path within the input file to
 * the array element at the corresponding index.
 * 
 * \param [in] container The container to retrieve indices from
 * \param [in] name The name of the array object in the input file
 *****************************************************************************
 */
std::vector<std::pair<std::string, std::string>> collectionIndicesWithPaths(const Container& container,
                                                                            const std::string& name);

/*!
 *******************************************************************************
 * \brief Updates the set of unexpected names to reflect an user-requested access
 * 
 * \param [in] accessedName The path within the input file that will be accessed
 * \param [inout] unexpectedNames The set of input file paths that have not yet
 * been requested by the user
 * 
 * \note To maintain consistency, this function should always be followed by an
 * access to a Reader
 *******************************************************************************
 */
void updateUnexpectedNames(const std::string& accessedName,
                           std::vector<std::string>& unexpectedNames);

}  // namespace detail

class Proxy;
/*!
 *******************************************************************************
 * \class Container
 *
 * \brief Provides functions to help define how individual Container and Field
 *        variables in an input file are expected to behave.  It also holds the
 *        Sidre Group to the individual Container.
 *
 * \see Inlet Field
 *******************************************************************************
 */
class Container : public Verifiable<Container>
{
public:
  /*!
   *****************************************************************************
   * \brief Constructor for the Container class.
   *
   * This class provides functions to define the behavior of the Container
   * data already read and stored in the given Sidre Group. This creates
   * the necessary Sidre Group's and views to store the given name and description.
   *
   * \param [in] name Name of the Container expected in the input file
   * \param [in] description Description of the Container
   * \param [in] reader Reference to the input file Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] unexpectedNames Reference to the global (relative to the Inlet
   * hierarchy) list of unexpected names
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for input feck this Container instance belongs to.
   * \param [in] reconstruct Whether or not to attempt to reconstruct child Containers
   * and Fields from the data in the sidre Group
   *****************************************************************************
   */
  Container(const std::string& name,
            const std::string& description,
            Reader& reader,
            axom::sidre::Group* sidreRootGroup,
            std::vector<std::string>& unexpectedNames,
            bool docEnabled = true,
            bool reconstruct = false);

  // Containers must be move-only - delete the implicit shallow copy constructor
  Container(const Container&) = delete;
  Container(Container&&) = default;

  virtual ~Container() = default;

  /*!
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group class for this Container.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for this Container class.
   *
   * \return Pointer to the Sidre Group for this Container
   *****************************************************************************
   */
  axom::sidre::Group* sidreGroup() { return m_sidreGroup; };
  /// \overload
  const axom::sidre::Group* sidreGroup() const { return m_sidreGroup; };

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
   * \brief Add an array of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addBoolArray(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of Integer Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addIntArray(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of Double Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addDoubleArray(const std::string& name,
                                        const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of String Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Container>& addStringArray(const std::string& name,
                                        const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Container& addStructArray(const std::string& name, const std::string& description = "");
  /*!
   *****************************************************************************
   * \brief Add a dictionary of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the Field
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Container>& addBoolDictionary(const std::string& name,
                                           const std::string& description = "");

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
                                          const std::string& description = "");
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
                                             const std::string& description = "");

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
                                             const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a dictionary of user-defined types to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Container& addStructDictionary(const std::string& name, const std::string& description = "");
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

  /*!
   *****************************************************************************
   * \brief Get a function from the input deck
   *
   * \param [in] name         Name of the function
   * \param [in] ret_type     The return type of the function
   * \param [in] arg_types    The argument types of the function
   * \param [in] description  Description of the function
   * \param [in] pathOverride The path within the input file to read from, if
   * different than the structure of the Sidre datastore
   *
   * \return Reference to the created Function
   *****************************************************************************
   */
  Verifiable<Function>& addFunction(const std::string& name,
                                    const FunctionTag ret_type,
                                    const std::vector<FunctionTag>& arg_types,
                                    const std::string& description = "",
                                    const std::string& pathOverride = "");

  /*!
   *******************************************************************************
   * \brief Returns a stored value of primitive type.
   * 
   * Retrieves a value of primitive type.
   * 
   * \param [in] name Name of the Field value to be gotten
   * \return The retrieved value
   * 
   * \tparam T The primitive type
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_inlet_primitive<T>::value, T>::type get(const std::string& name) const
  {
    if(!hasField(name))
    {
      const std::string msg = fmt::format(
        "[Inlet] Field with specified path "
        "does not exist: {0}",
        name);
      SLIC_ERROR(msg);
    }
    return getField(name).get<T>();
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored value of user-defined type.
   * 
   * Retrieves a value of user-defined type.
   * 
   * \param [in] name The name of the subcontainer representing the root of the object
   * If nothing is passed, the calling container is interpreted as the roof of the object
   * \return The retrieved value
   * \pre Requires a specialization of FromInlet for T
   * 
   * \tparam T The user-defined type
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<!detail::is_inlet_primitive<T>::value && !detail::is_inlet_array<T>::value &&
                            !detail::is_inlet_dict<T>::value && !detail::is_std_vector<T>::value,
                          T>::type
  get(const std::string& name = "") const
  {
    static_assert(detail::has_FromInlet_specialization<T>::value,
                  "To read a user-defined type, specialize FromInlet<T>");
    FromInlet<T> from_inlet;
    if(name.empty())
    {
      return from_inlet(*this);
    }
    else
    {
      if(!hasContainer(name))
      {
        std::string msg = fmt::format("[Inlet] Container with name '{0}' does not exist", name);
        SLIC_ERROR(msg);
      }
      return from_inlet(getContainer(name));
    }
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored collection.
   * 
   * Retrieves a collection of user-defined type.
   * 
   * \return The retrieved collection
   * 
   * \tparam T The collection type, i.e., T = std::unordered_map<K, V>
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_inlet_array<T>::value || detail::is_inlet_dict<T>::value, T>::type
  get() const
  {
    using Key = typename T::key_type;
    using Val = typename T::mapped_type;
    // This needs to work transparently for both references to the underlying
    // internal container and references using the same path as the data file
    if(isCollectionGroup(m_name))
    {
      return getCollection<Key, Val>();
    }
    // If the collection group is a child container, retrieve it and copy its contents
    // into the result
    else
    {
      return getContainer(detail::COLLECTION_GROUP_NAME).getCollection<Key, Val>();
    }
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored collection as a contiguous array.
   * 
   * Retrieves a collection of user-defined type.
   * 
   * \return The values in the retrieved collection
   * 
   * \tparam T The collection type, i.e., T = std::vector<V>
   * 
   * \note Elements in the returned array will be in ascending order by index,
   * regardless of index contiguity or base index
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_std_vector<T>::value, T>::type get() const
  {
    // Only allow retrieval of std::vectors from integer-keyed collections
    using Key = int;
    using Val = typename T::value_type;
    auto map = get<std::unordered_map<Key, Val>>();

    // Retrieve and sort the indices to provide consistent behavior regardless
    // of index contiguity or base index
    std::vector<Key> indices;
    indices.reserve(map.size());

    for(const auto& entry : map)
    {
      indices.push_back(entry.first);
    }
    std::sort(indices.begin(), indices.end());

    std::vector<Val> result;
    result.reserve(map.size());

    for(const Key index : indices)
    {
      // Safe to move from the map as it won't get used afterwards
      result.push_back(std::move(map[index]));
    }

    return result;
  }

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the container for either a Field/Container subobject
   * 
   * Returns a reference via a lightweight proxy object to the element in the 
   * datastore at the index specified by the name.  This can be a field 
   * or a container.
   * 
   * \param [in] name The name of the subobject
   * \return The retrieved array
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) const;

  Container& required(bool isRequired = true) override;

  bool isRequired() const override;

  using Verifiable<Container>::registerVerifier;

  Container& registerVerifier(Verifier lambda) override;

  bool verify(std::vector<VerificationError>* errors = nullptr) const override;

  /*!
   *****************************************************************************
   * \brief Set the strictness of this Container.
   *
   * Set whether this Container is strict, or not - i.e., whether entries other
   * than those added to the schema should be allowed.
   * The default behavior is to not be strict.
   *
   * \param [in] isStrict Boolean value of whether Container is strict
   *
   * \return Reference to this instance of Container
   *****************************************************************************
   */
  Container& strict(bool isStrict = true);

  /*!
   *****************************************************************************
   * \brief Return whether a Container or Field with the given name is present in 
   * this Container's subtree.
   *
   * \return Boolean value indicating whether this Container's subtree contains a
   * Field or Container with the given name.
   *****************************************************************************
   */
  bool contains(const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Returns whether this container or any of its subcontainers exist, 
   * i.e., if they contain a Field or Function that exists
   *****************************************************************************
   */
  bool exists() const;

  /*!
   *****************************************************************************
   * \brief Returns whether this container or any of its subcontainers were
   * provided in the input file, i.e., if they contain a Field or Function that
   * was provided in the input file
   *****************************************************************************
   */
  bool isUserProvided() const;

  /*!
   *****************************************************************************
   * \brief  Return whether a Container or Field with the given name were
   * provided in the input file
   *
   * \param [in] name The path relative to the calling Container to search
   *****************************************************************************
   */
  bool isUserProvided(const std::string& name) const;

  /*!
   *****************************************************************************
   * \return An unordered map from Field names to the child Field pointers for 
   * this Container.
   *****************************************************************************
   */
  const std::unordered_map<std::string, std::unique_ptr<Field>>& getChildFields() const;

  /*!
   *****************************************************************************
   * \return An unordered map from Container names to the child Container pointers for 
   * this Container.
   *****************************************************************************
   */
  const std::unordered_map<std::string, std::unique_ptr<Container>>& getChildContainers() const;

  /*!
   *****************************************************************************
   * \return An unordered map from Function names to the child Function pointers for 
   * this Container.
   *****************************************************************************
   */
  const std::unordered_map<std::string, std::unique_ptr<Function>>& getChildFunctions() const;

  /*!
   *****************************************************************************
   * \return The full name of this Container.
   *****************************************************************************
   */
  std::string name() const;

  /*!
   *****************************************************************************
   * \brief Returns the list of unexpected names "below" the calling container,
   * i.e., entries in the input file structure (e.g., a Lua table or
   * YAML dictionary) corresponding to the calling container that were not
   * requested/retrieved via an add* call
   *****************************************************************************
   */
  std::vector<std::string> unexpectedNames() const;

  /*!
   *****************************************************************************
   * \brief Add a Field to the input file schema.
   *
   * Adds a primitive Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Container expected in the input file
   * \param [in] description Description of the Container
   * \param [in] forArray Whether the primitive is in an array, in which
   * case the provided value should be inserted instead of the one read from
   * the input file
   * \param [in] val A provided value, will be overwritten if found at specified
   * path in input file
   * \param [in] pathOverride The path within the input file to read from, if
   * different than the structure of the Sidre datastore
   *
   * \return Reference to the created Field
   *****************************************************************************
   */

  template <typename T,
            typename SFINAE = typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  VerifiableScalar& addPrimitive(const std::string& name,
                                 const std::string& description = "",
                                 bool forArray = false,
                                 T val = T {},
                                 const std::string& pathOverride = "");

private:
  /*!
   *****************************************************************************
   * \brief Add a Container to the input file schema.
   *
   * Adds a Container to the input file schema. Containers hold a varying amount Fields
   * defined by the user.  By default, it is not required unless marked with
   * Container::isRequired(). This creates the Sidre Group class with the given name and
   * stores the given description.
   *
   * \param [in] name Name of the Container expected in the input file
   * \param [in] description Description of the Container
   *
   * \return Reference to the created Container
   *****************************************************************************
   */
  Container& addContainer(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of primitive Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   * \param [in] isDict Whether to use string-valued keys for the collection
   * \param [in] pathOverride The path within the input file to read from, if
   * different than the structure of the Sidre datastore
   *
   * \return Reference to the created Field
   *****************************************************************************
   */
  template <typename T,
            typename SFINAE = typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  Verifiable<Container>& addPrimitiveArray(const std::string& name,
                                           const std::string& description = "",
                                           const bool isDict = false,
                                           const std::string& pathOverride = "");

  /*!
   *****************************************************************************
   * \brief Return whether a Container with the given name is present in this Container's subtree.
   *
   * \return Boolean value indicating whether the calling Container's subtree
   * contains a Container with the given name.
   *****************************************************************************
   */
  bool hasContainer(const std::string& containerName) const;

  /*!
   *****************************************************************************
   * \brief Return whether a Field with the given name is present in this Container's
   *  subtree.
   *
   * \return Boolean value indicating whether the calling Container's subtree
   * contains a Field with the given name.
   *****************************************************************************
   */
  bool hasField(const std::string& fieldName) const;

  /*!
   *****************************************************************************
   * \brief Return whether a Function with the given name is present in this Container's
   *  subtree.
   *
   * \return Boolean value indicating whether the calling Container's subtree
   * contains a Function with the given name.
   *****************************************************************************
   */
  bool hasFunction(const std::string& fieldName) const;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Container.
   * 
   * \param [in] The string indicating the target name of the Container to be searched for.
   * 
   * \return The Container matching the target name.
   *****************************************************************************
   */
  Container& getContainer(const std::string& containerName) const;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Field.
   * 
   * \param [in] The string indicating the target name of the Field to be searched for.
   * 
   * \return The Field matching the target name.
   *****************************************************************************
   */
  Field& getField(const std::string& fieldName) const;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Function.
   * 
   * \param [in] The string indicating the target name of the Function to be searched for.
   * 
   * \return The Function matching the target name.
   *****************************************************************************
   */
  Function& getFunction(const std::string& funcName) const;

  /*!
   *****************************************************************************
   * \brief Helper method template for adding primitives
   * 
   * Adds the value at the templated type to the sidre group
   * 
   * \param [inout] sidreGroup The group to add the primitive view to
   * \param [in] lookupPath The path within the input file to read from
   * \param [in] forArray Whether the primitive is in an array, in which
   * case the provided value should be inserted instead of the one read from
   * the input file
   * \param [in] val A provided value, will be overwritten if found at specified
   * path in input file
   *
   * \return Type ID for the inserted view
   *****************************************************************************
   */
  template <typename T,
            typename SFINAE = typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  axom::sidre::DataTypeId addPrimitiveHelper(axom::sidre::Group* sidreGroup,
                                             const std::string& lookupPath,
                                             bool forArray,
                                             T val);

  /*!
   *****************************************************************************
   * \brief Helper method template for adding primitives
   * 
   * Reads an array at the provided path into the provided container
   * 
   * \param [inout] container The inlet::Container to add the array to 
   * \param [in] lookupPath The path within the input file to read from
   * \param [in] isDict Whether to use string keys
   * 
   *****************************************************************************
   */
  template <typename T,
            typename SFINAE = typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  void addPrimitiveArrayHelper(Container& container,
                               const std::string& lookupPath,
                               bool isDict = false);

  /*!
   *****************************************************************************
   * \brief Creates the basic Sidre Group for this Container and stores the given
   *        information
   *
   * \return Pointer to the created Sidre Group for this Container
   *****************************************************************************
   */
  axom::sidre::Group* createSidreGroup(const std::string& name, const std::string& description);

  /*!
   *****************************************************************************
   * \brief Adds the Field.
   * 
   * \param [in] The Sidre Group corresponding to the Field that will be added.
   * \param [in] The type ID
   * \param [in] The complete Container sequence for the Container this Field will be added to.
   * \param [in] The Container sequence for the Container this Field will be added to, 
   * relative to this Container.
   * 
   * \return The child Field matching the target name.
   *****************************************************************************
   */
  Field& addField(axom::sidre::Group* sidreGroup,
                  axom::sidre::DataTypeId type,
                  const std::string& fullName,
                  const std::string& name);

  /*!
   *****************************************************************************
   * \brief Adds the Function.
   * 
   * \param [in] The Sidre Group corresponding to the Function that will be added.
   * \param [in] func The actual callable to store
   * \param [in] The complete Container sequence for the Container this Function will be added to.
   * \param [in] The Container sequence for the Container this Function will be added to, 
   * relative to this Container.
   * 
   * \return The child Function matching the target name.
   *****************************************************************************
   */
  Function& addFunctionInternal(axom::sidre::Group* sidreGroup,
                                FunctionVariant&& func,
                                const std::string& fullName,
                                const std::string& name);

  axom::sidre::View* baseGet(const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief This is an internal helper that returns the pointer-to-member for
   * the unordered_map of children of requested type.
   * 
   * \tparam T The type of the child to search for (Field/Container/Function)
   *****************************************************************************
   */
  template <typename T>
  static std::unordered_map<std::string, std::unique_ptr<T>> Container::* getChildren();

  /*!
   *****************************************************************************
   * \brief This is an internal helper. It return whether this Container has a child 
   * with the given name and type.
   * 
   * \tparam T The type of the child to search for (Field/Container/Function)
   * \return Boolean value of whether this Container has the child.
   *****************************************************************************
   */
  template <typename T>
  bool hasChild(const std::string& childName) const;

  /*!
   *****************************************************************************
   * \brief This retrieves the child of requested name and type.
   * 
   * \param [in] The string indicating the target name of the child to be searched for.
   * 
   * \tparam T The type of the child to search for (Field/Container/Function)
   * \return The child matching the target name. If no such child is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  template <typename T>
  T* getChildInternal(const std::string& childName) const;

  /*!
   *****************************************************************************
   * \brief Get a collection represented as an unordered map from the input file
   *****************************************************************************
   */
  template <typename Key, typename Val>
  std::unordered_map<Key, Val> getCollection() const
  {
    std::unordered_map<Key, Val> map;
    for(const auto& indexLabel : detail::collectionIndices(*this))
    {
      if(detail::matchesKeyType<Key>(indexLabel))
      {
        map[detail::toIndex<Key>(indexLabel)] = get<Val>(detail::indexToString(indexLabel));
      }
    }
    return map;
  }

  /*!
   *******************************************************************************
   * \brief Adds a group containing the indices of a collection to the calling 
   * container and optionally a subcontainer for each index
   * 
   * \param [in] indices The indices to add
   * \param [in] description The optional description of the subcontainers
   * \param [in] add_containers Whether to add a subcontainer for each index
   * \tparam Key The type of the indices to add
   *******************************************************************************
   */
  template <typename Key>
  void addIndicesGroup(const std::vector<Key>& indices,
                       const std::string& description = "",
                       const bool add_containers = false);

  /*!
   *****************************************************************************
   * \brief Add a collection of user-defined type to the input file schema.
   *
   * \param [in] name Name of the collection
   * \param [in] description Description of the collection
   *
   * \return Reference to the created collection
   *****************************************************************************
   */
  template <typename Key>
  Container& addStructCollection(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Returns true if the calling object is part of a struct collection,
   * i.e., an array or dictionary of user-defined type
   *****************************************************************************
   */
  bool isStructCollection() const { return m_sidreGroup->hasView(detail::STRUCT_COLLECTION_FLAG); }

  /*!
   *****************************************************************************
   * \brief Calls a function on the subcontainers corresponding to the elements
   * of the collection held by this container
   * 
   * \param [in] func The function to apply to individual collection elements
   * 
   * \pre The calling container must be a struct collection, i.e., isStructCollection()
   * returns true
   * 
   * \pre The function must accept a single argument of type Container&
   * 
   *****************************************************************************
   */
  template <typename Func>
  void forEachCollectionElement(Func&& func) const;

  /*!
   *****************************************************************************
   * \brief Applies a provided function to nested elements of the calling table
   * and stores the result in a range pointed to by an output iterator \a output
   * 
   * \pre The function \a func must accept two arguments of type Table& and
   * const std::string&, respectively.  
   * 
   * This function will pass to the provided function the nested table
   * as the first argument and the path of the nested element in the input file
   * as the second argument, when applicable.
   * 
   * \param [out] output An iterator to the beginning of the output range
   * \param [in] name The name to append to the path described above
   * \param [in] func The function to apply to individual nested elements
   * 
   * \return Whether the calling container had any nested elements (or
   * was a struct container)
   * 
   * \note This function can be thought of as a variant of std::transform that
   * operates on nested elements instead of a provided input range
   *****************************************************************************
   */
  template <typename OutputIt, typename Func>
  bool transformFromNestedElements(OutputIt output, const std::string& name, Func&& func) const;

  std::string m_name;
  Reader& m_reader;
  // Inlet's Root Sidre Group
  axom::sidre::Group* m_sidreRootGroup {nullptr};
  // This Container's Sidre Group
  axom::sidre::Group* m_sidreGroup {nullptr};
  // Hold a reference to the global set of unexpected names so it can be updated when
  // things are added to this Container
  std::vector<std::string>& m_unexpectedNames;
  bool m_docEnabled;
  std::unordered_map<std::string, std::unique_ptr<Container>> m_containerChildren;
  std::unordered_map<std::string, std::unique_ptr<Field>> m_fieldChildren;
  std::unordered_map<std::string, std::unique_ptr<Function>> m_functionChildren;
  Verifier m_verifier;

  // Used for ownership only - need to take ownership of these so children
  // and their aggregates have identical lifetime
  std::vector<AggregateVerifiable<Container>> m_aggregate_containers;
  std::vector<AggregateField> m_aggregate_fields;
  std::vector<AggregateVerifiable<Function>> m_aggregate_funcs;

  // Used when the calling Container is a struct collection within a struct collection
  // Need to delegate schema-defining calls (add*) to the elements of the nested
  // collection
  std::vector<std::reference_wrapper<Container>> m_nested_aggregates;
};

}  // namespace inlet
}  // namespace axom

#endif
