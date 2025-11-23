// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/ndarray.h>
#include <nanobind/make_iterator.h>

#include "axom/core/Types.hpp"
#include "core/SidreTypes.hpp"
#include "core/Buffer.hpp"
#include "core/View.hpp"
#include "core/DataStore.hpp"
#include "core/Group.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace axom
{
namespace sidre
{

// Helper to map TypeID to nanobind dtype
nb::dlpack::dtype typeIDToDtype(DataTypeId id)
{
  switch(id)
  {
  case INT32_ID:
    return nb::dtype<int>();
  case INT64_ID:
    return nb::dtype<int64_t>();

  // DOUBLE_ID also has same value
  case FLOAT64_ID:
    return nb::dtype<double>();
  default:
    SLIC_ERROR("DataTypeId unsupported for numpy");
    return nb::dtype<double>();
  }
}

/*!
 * \brief Returns a View as a numpy array.
 *
 * \note Max dimensions (DMAX) is currently set to 10.
 * \pre data description must have been applied.
 */
nb::ndarray<nb::numpy> viewToNumpyArray(View& self)
{
  // Manually applying offset
  void* data = self.getVoidPtr();
  char* data_with_offset = static_cast<char*>(data) + (self.getOffset() * self.getBytesPerElement());
  data = static_cast<void*>(data_with_offset);

  constexpr int DMAX = 10;

  IndexType shapeOutput[DMAX];
  size_t ndims = self.getShape(DMAX, shapeOutput);
  size_t shape[DMAX];
  for(size_t i = 0; i < ndims; i++)
  {
    shape[i] = static_cast<size_t>(shapeOutput[i]);
  }

  // TODO This is tricky and difficult to understand
  // Delete 'data' when the 'owner' capsule expires
  // nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<char*>(p); });

  // For external memory (numpy owns it), no deletion takes place
  nb::capsule owner(data, [](void*) noexcept { });

  // When stride is not default of 1, guaranteed that shape is 1D.
  int64_t* strides = nullptr;
  int64_t stride_array[1];
  if(self.getStride() != 1)
  {
    stride_array[0] = static_cast<int64_t>(self.getStride());
    strides = stride_array;
  }

  DataTypeId id = self.getTypeID();

  return nb::ndarray<nb::numpy>(
    /* data = */ data,
    /* ndim = */ ndims,
    /* shape = */ shape,
    /* owner = */ owner,
    /* strides = */ strides,
    /* dtype = */ typeIDToDtype(id));
}

/*!
 * \brief Returns a Buffer as a numpy array.
 *
 * \pre data description must have been applied.
 */
nb::ndarray<nb::numpy> bufferToNumpyArray(Buffer& self)
{
  void* data = self.getVoidPtr();

  size_t shape[1] = {static_cast<size_t>(self.getNumElements())};

  // TODO This is tricky and difficult to understand
  // Delete 'data' when the 'owner' capsule expires
  // nb::capsule owner(data, [](void* p) noexcept { delete[] static_cast<char*>(p); });

  // For external memory (numpy owns it), no deletion takes place
  nb::capsule owner(data, [](void*) noexcept { });

  DataTypeId id = self.getTypeID();

  return nb::ndarray<nb::numpy>(
    /* data = */ data,
    /* ndim = */ 1,
    /* shape = */ shape,
    /* owner = */ owner,
    /* strides = */ nullptr,
    /* dtype = */ typeIDToDtype(id));
}

/*!
 * \brief Helper function to bind iterator types
 *
 * \tparam IteratorType The type of the iterator
 */
template <typename IteratorType>
void bindIterator(nb::module_& m, const char* iterator_name)
{
  // Create non-const and const python iterators using nanobind's make_iterator helper
  // We need to specify the ValueType (4th template parameter) as a reference since IndexedCollections
  // have pointer/reference semantics and the underlying type (e.g. sidre::View) are not copyable
  nb::class_<IteratorType>(m, iterator_name)
    .def("__len__", &IteratorType::size)
    .def(
      "__iter__",
      [iterator_name](IteratorType& self) {
        return nb::make_iterator<nb::rv_policy::reference,
                                 decltype(self.begin()),
                                 decltype(self.end()),
                                 typename IteratorType::CollectionType::value_type&>(
          nb::type<IteratorType>(),
          iterator_name,
          self.begin(),
          self.end());
      },
      nb::keep_alive<0, 1>());
}

NB_MODULE(pysidre, m_sidre)
{
  m_sidre.doc() = "A python extension for Axom's Sidre component";

  m_sidre.attr("InvalidIndex") = axom::InvalidIndex;
  m_sidre.attr("InvalidName") = axom::utilities::string::InvalidName;

  m_sidre.def("indexIsValid", &indexIsValid, "Returns true if idx is valid, else false.");
  m_sidre.def("nameIsValid", &nameIsValid, "Returns true if name is valid, else false.");

#if defined(AXOM_USE_HDF5)
  m_sidre.attr("AXOM_USE_HDF5") = true;
#else
  m_sidre.attr("AXOM_USE_HDF5") = false;
#endif

  // Bind the DataTypeId enum (TypeID alias)
  nb::enum_<DataTypeId>(m_sidre, "TypeID")
    .value("NO_TYPE_ID", NO_TYPE_ID)
    .value("INT8_ID", INT8_ID)
    .value("INT16_ID", INT16_ID)
    .value("INT32_ID", INT32_ID)
    .value("INT64_ID", INT64_ID)
    .value("UINT8_ID", UINT8_ID)
    .value("UINT16_ID", UINT16_ID)
    .value("UINT32_ID", UINT32_ID)
    .value("UINT64_ID", UINT64_ID)
    .value("FLOAT32_ID", FLOAT32_ID)
    .value("FLOAT64_ID", FLOAT64_ID)
    .value("CHAR8_STR_ID", CHAR8_STR_ID)
    .value("INT_ID", INT_ID)
    .value("UINT_ID", UINT_ID)
    .value("LONG_ID", LONG_ID)
    .value("ULONG_ID", ULONG_ID)
    .value("FLOAT_ID", FLOAT_ID)
    .value("DOUBLE_ID", DOUBLE_ID)
    .export_values();

  // Bindings to support iterating collections
  using AttributeIterator = axom::IndexedCollection<Attribute>::iterator_adaptor;
  using BufferIterator = axom::IndexedCollection<Buffer>::iterator_adaptor;
  using GroupIterator = axom::IndexedCollection<Group>::iterator_adaptor;
  using ViewIterator = axom::IndexedCollection<View>::iterator_adaptor;

  bindIterator<AttributeIterator>(m_sidre, "AttributeIterator");
  bindIterator<BufferIterator>(m_sidre, "BufferIterator");
  bindIterator<GroupIterator>(m_sidre, "GroupIterator");
  bindIterator<ViewIterator>(m_sidre, "ViewIterator");

  // Bindings for the DataStore class
  nb::class_<DataStore>(m_sidre, "DataStore")
    .def(nb::init<>())
    .def("getRoot",
         nb::overload_cast<>(&DataStore::getRoot),
         nb::rv_policy::reference,
         "Return pointer to the root Group")
    .def("getNumBuffers", &DataStore::getNumBuffers, "Return number of Buffers in the DataStore")
    .def("hasBuffer",
         &DataStore::hasBuffer,
         "Return true if DataStore owns a Buffer with given index; else false")
    .def("getBuffer",
         &DataStore::getBuffer,
         nb::rv_policy::reference,
         "Return pointer to Buffer object with the given index")

    .def("createBuffer",
         nb::overload_cast<>(&DataStore::createBuffer),
         nb::rv_policy::reference,
         "Create an undescribed Buffer object")
    .def("createBuffer",
         nb::overload_cast<TypeID, IndexType>(&DataStore::createBuffer),
         nb::rv_policy::reference,
         "Create a Buffer object with specified type and number of elements")
    .def("destroyBuffer",
         nb::overload_cast<Buffer*>(&DataStore::destroyBuffer),
         "Remove Buffer from the DataStore and destroy it and its data")
    .def("destroyBuffer",
         nb::overload_cast<IndexType>(&DataStore::destroyBuffer),
         "Remove Buffer with given index from the DataStore and destroy it and its data.")
    .def("destroyAllBuffers",
         &DataStore::destroyAllBuffers,
         "Remove all Buffers from the DataStore and destroy them and their data")
    .def("getFirstValidBufferIndex",
         &DataStore::getFirstValidBufferIndex,
         "Return first valid Buffer index")
    .def("getNextValidBufferIndex",
         &DataStore::getNextValidBufferIndex,
         "Return next valid Buffer index after given index")

    .def("generateBlueprintIndex",
         nb::overload_cast<const std::string&, const std::string&, const std::string&, int>(
           &DataStore::generateBlueprintIndex),
         "Generate a Conduit Blueprint index based on a mesh in stored in this DataStore.")
    .def("buffers",
         nb::overload_cast<>(&DataStore::buffers),
         nb::rv_policy::reference,
         "Return an iterator over Buffers")

    .def("getNumAttributes",
         &DataStore::getNumAttributes,
         "Return number of Attributes in the DataStore")
    .def("createAttributeScalar",
         &DataStore::createAttributeScalar<int>,
         nb::rv_policy::reference,
         "Create an Attribute object with a default int scalar value",
         nb::arg("name"),
         nb::arg("default_value").noconvert())
    .def("createAttributeScalar",
         &DataStore::createAttributeScalar<double>,
         nb::rv_policy::reference,
         "Create an Attribute object with a default float (C++ double) scalar value",
         nb::arg("name"),
         nb::arg("default_value").noconvert())
    .def("createAttributeString",
         &DataStore::createAttributeString,
         nb::rv_policy::reference,
         "Create an Attribute object with a default string value")
    .def("hasAttribute",
         nb::overload_cast<const std::string&>(&DataStore::hasAttribute, nb::const_),
         "Return true if DataStore has created attribute name, else false")
    .def("hasAttribute",
         nb::overload_cast<IndexType>(&DataStore::hasAttribute, nb::const_),
         "Return true if DataStore has created attribute with index, else false")
    .def("destroyAttribute",
         nb::overload_cast<const std::string&>(&DataStore::destroyAttribute),
         "Remove Attribute from the DataStore and destroy it and its data")
    .def("destroyAttribute",
         nb::overload_cast<IndexType>(&DataStore::destroyAttribute),
         "Remove Attribute with given index from the DataStore and destroy it and its data")
    .def("destroyAttribute",
         nb::overload_cast<Attribute*>(&DataStore::destroyAttribute),
         "Remove Attribute from the DataStore and destroy it and its data")
    .def("destroyAllAttributes",
         &DataStore::destroyAllAttributes,
         "Remove all Attributes from the DataStore and destroy them and their data")
    .def("getAttribute",
         nb::overload_cast<IndexType>(&DataStore::getAttribute),
         nb::rv_policy::reference,
         "Return pointer to non-const Attribute with given index")
    .def("getAttribute",
         nb::overload_cast<const std::string&>(&DataStore::getAttribute),
         nb::rv_policy::reference,
         "Return pointer to non-const Attribute with given name")

    // Requires conduit::Node information
    // .def("saveAttributeLayout",
    //      &DataStore::saveAttributeLayout,
    //      "Copy Attribute and default value to Conduit node. Return true if attributes were copied.")
    // .def("loadAttributeLayout",
    //      &DataStore::loadAttributeLayout,
    //      "Create attributes from name/value pairs in node['attribute'].")

    .def("getFirstValidAttributeIndex",
         &DataStore::getFirstValidAttributeIndex,
         "Return first valid Attribute index in DataStore object "
         "(i.e., smallest index over all Attributes)")
    .def("getNextValidAttributeIndex",
         &DataStore::getNextValidAttributeIndex,
         "Return next valid Attribute index in DataStore object after given index"
         "(i.e., smallest index over all Attribute indices larger than given one)")
    .def("attributes",
         nb::overload_cast<>(&DataStore::attributes),
         nb::rv_policy::reference,
         "Return an iterator over Attributes")

    // Nanobind fails compilation on blueos
    // #ifdef AXOM_USE_MPI
    //     .def("generateBlueprintIndex",
    //          nb::overload_cast<MPI_Comm, const std::string&, const std::string&, const std::string&>(
    //            &DataStore::generateBlueprintIndex),
    //          "Generate a Conduit Blueprint index from a distributed mesh stored in this Datastore")
    // #endif
    .def("print",
         nb::overload_cast<>(&DataStore::print, nb::const_),
         "Print JSON description of the DataStore");

  // Bindings for the Buffer class
  nb::class_<Buffer>(m_sidre, "Buffer")
    .def("getIndex", &Buffer::getIndex, "Return the unique index of this Buffer object.")
    .def("getNumViews", &Buffer::getNumViews, "Return number of Views this Buffer is attached to.")
    // .def("getVoidPtr", &Buffer::getVoidPtr, "Return void-pointer to data held by Buffer.")
    .def("getDataArray", &bufferToNumpyArray, "Return the data held by the Buffer as a numpy array.")
    .def("getTypeID", &Buffer::getTypeID, "Return type of data owned by this Buffer object.")
    .def("getNumElements",
         &Buffer::getNumElements,
         "Return total number of data elements owned by this Buffer object.")
    .def("getTotalBytes",
         &Buffer::getTotalBytes,
         "Return total number of bytes of data owned by this Buffer object.")
    .def("getBytesPerElement",
         &Buffer::getBytesPerElement,
         "Return the number of bytes per element owned by this Buffer object.")
    .def("isAllocated",
         &Buffer::isAllocated,
         "Return true if Buffer has been (re)allocated with length >= 0, else false")
    .def("isDescribed", &Buffer::isDescribed, "Return true if data description exists")
    .def("describe",
         &Buffer::describe,
         "Describe a Buffer with data type and number of elements.",
         nb::arg("type"),
         nb::arg("num_elems"))
    .def("allocate",
         nb::overload_cast<int>(&Buffer::allocate),
         "Allocate data for a Buffer.",
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("allocate",
         nb::overload_cast<TypeID, IndexType, int>(&Buffer::allocate),
         "Allocate Buffer with data type and number of elements.",
         nb::arg("type"),
         nb::arg("num_elems"),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("reallocate",
         &Buffer::reallocate,
         "Reallocate data to given number of elements.",
         nb::arg("num_elems"))
    .def("print",
         nb::overload_cast<>(&Buffer::print, nb::const_),
         "Print JSON description of Buffer to std::cout.");

  // Bindings for the View class
  nb::class_<View>(m_sidre, "View")
    .def("getIndex", &View::getIndex, "Return the index of the View within its owning Group.")
    .def("getName", &View::getName, "Return the name of the View.")
    .def("getPath", &View::getPath, "Return the path of the View's owning Group object.")
    .def("getPathName",
         &View::getPathName,
         "Return the full path of the View object, including its name.")
    .def("getOwningGroup",
         nb::overload_cast<>(&View::getOwningGroup),
         nb::rv_policy::reference,
         "Return the owning Group of the View.")
    .def("hasBuffer", &View::hasBuffer, "Check if the View has an associated Buffer object.")
    .def("getBuffer",
         nb::overload_cast<>(&View::getBuffer),
         nb::rv_policy::reference,
         "Return the associated Buffer object (non-const).")
    .def("isExternal", &View::isExternal, "Check if the View holds external data.")
    .def("isAllocated", &View::isAllocated, "Check if the View's data is allocated.")
    .def("isApplied", &View::isApplied, "Check if the data description has been applied.")
    .def("isDescribed", &View::isDescribed, "Check if the View has a data description.")
    .def("isEmpty", &View::isEmpty, "Check if the View is empty.")
    .def("isOpaque", &View::isOpaque, "Check if the View is opaque.")
    .def("isScalar", &View::isScalar, "Check if the View contains a scalar value.")
    .def("isString", &View::isString, "Check if the View contains a string value.")
    .def("getTypeID", &View::getTypeID, "Return the type ID of the View's data.")
    .def("getTotalBytes",
         &View::getTotalBytes,
         "Return the total number of bytes described by the View.")
    .def("getNumElements",
         &View::getNumElements,
         "Return the total number of elements described by the View.")
    .def("getBytesPerElement",
         &View::getBytesPerElement,
         "Return the number of bytes per element described by the View.")
    .def("getOffset",
         &View::getOffset,
         "Return the offset in number of elements for the data described by the View.")
    .def("getStride",
         &View::getStride,
         "Return the stride in number of elements for the data described by the View.")
    .def("getNumDimensions",
         &View::getNumDimensions,
         "Return the dimensionality of the View's data.")
    .def(
      "getShape",
      [](View& self, int ndims, nb::ndarray<IndexType>& shape) {
        SLIC_ERROR_IF(static_cast<size_t>(ndims) > shape.size(),
                      "getShape() - shape array size (" << shape.size()
                                                        << ") must be greater or equal to ndims ("
                                                        << static_cast<size_t>(ndims) << ")");
        int ret = self.getShape(ndims, shape.data());
        return nb::make_tuple(ret, shape);
      },
      "Return number of dimensions in data view and shape information"
      " of this data view object."
      " ndims - maximum number of dimensions to return."
      " shape - user supplied numpy 1D array assumed to be ndims long.")

    .def("allocate",
         nb::overload_cast<int>(&View::allocate),
         nb::rv_policy::reference,
         "Allocate data for the View.",
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("allocate",
         nb::overload_cast<TypeID, IndexType, int>(&View::allocate),
         "Allocate data for the View with type and number of elements.",
         nb::rv_policy::reference,
         nb::arg("type"),
         nb::arg("num_elems"),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("reallocate",
         nb::overload_cast<IndexType>(&View::reallocate),
         nb::rv_policy::reference,
         "Reallocate data for the View.")
    .def("attachBuffer",
         nb::overload_cast<Buffer*>(&View::attachBuffer),
         nb::rv_policy::reference,
         "Attach a Buffer object to the View.")
    .def("attachBuffer",
         nb::overload_cast<TypeID, IndexType, Buffer*>(&View::attachBuffer),
         nb::rv_policy::reference,
         "Describe the data view and attach Buffer object.")
    .def("attachBuffer",
         nb::overload_cast<TypeID, int, const IndexType*, Buffer*>(&View::attachBuffer),
         nb::rv_policy::reference,
         "Describe the data view and attach Buffer object")

    .def("clear", &View::clear, "Clear data and metadata from the View.")
    .def("apply", nb::overload_cast<>(&View::apply), "Apply the View's description to its data.")
    .def("apply",
         nb::overload_cast<IndexType, IndexType, IndexType>(&View::apply),
         nb::rv_policy::reference,
         "Apply data description with number of elements, offset, and stride.",
         nb::arg("num_elems"),
         nb::arg("offset") = 0,
         nb::arg("stride") = 1)
    .def("apply",
         nb::overload_cast<TypeID, IndexType, IndexType, IndexType>(&View::apply),
         nb::rv_policy::reference,
         "Apply data description with type, number of elements, offset, and stride.",
         nb::arg("type"),
         nb::arg("num_elems"),
         nb::arg("offset").noconvert() = 0,
         nb::arg("stride").noconvert() = 1)
    .def(
      "apply",
      [](View& self, TypeID type, int ndims, nb::ndarray<IndexType>& shape) {
        return self.apply(type, ndims, shape.data());
      },
      nb::rv_policy::reference,
      "Apply data description with type and numpy shape.")
    .def("setScalar",
         &View::setScalar<int>,
         nb::rv_policy::reference,
         "Set the View to hold a scalar value (int).",
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("setScalar",
         &View::setScalar<double>,
         nb::rv_policy::reference,
         "Set the View to hold a scalar value (python float, C++ double).",
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("setString",
         &View::setString,
         "Set the View to hold a string value.",
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def(
      "setExternalData",
      [](View& self, const nb::ndarray<>& external_ptr) {
        return self.setExternalDataPtr(external_ptr.data());
      },
      nb::rv_policy::reference,
      "Set the View to hold undescribed external data (numpy array).")
    .def(
      "setExternalData",
      [](View& self, TypeID type, IndexType num_elems, const nb::ndarray<>& external_ptr) {
        return self.setExternalDataPtr(type, num_elems, external_ptr.data());
      },
      nb::rv_policy::reference,
      "Set the View to hold described external data  (numpy array).")
    .def(
      "setExternalData",
      [](View& self,
         TypeID type,
         int ndims,
         const nb::ndarray<IndexType>& shape,
         const nb::ndarray<>& external_ptr) {
        return self.setExternalDataPtr(type, ndims, shape.data(), external_ptr.data());
      },
      nb::rv_policy::reference,
      "Set the View to hold described external data (numpy array).")

    .def("getString",
         &View::getString,
         nb::rv_policy::reference,
         "Return the string contained in the View.")
    .def("getDataArray", &viewToNumpyArray, "Return the data held by the View as a numpy array.")

    .def("getDataInt",
         &View::getData<int>,
         nb::rv_policy::reference,
         "Return the scalar data held by the View as an python int type.")
    .def("getDataFloat",
         &View::getData<double>,
         nb::rv_policy::reference,
         "Return the data held by the View as a python float type (C++ double).")
    .def("print",
         nb::overload_cast<>(&View::print, nb::const_),
         "Print JSON description of the View.")
    .def("rename", &View::rename, "Change the name of the View.")

    // Attribute accessors
    .def("getAttribute",
         nb::overload_cast<IndexType>(&View::getAttribute),
         nb::rv_policy::reference,
         "Get Attribute by index")
    .def("getAttribute",
         nb::overload_cast<const std::string&>(&View::getAttribute),
         nb::rv_policy::reference,
         "Get Attribute by name")

    .def("hasAttributeValue",
         nb::overload_cast<IndexType>(&View::hasAttributeValue, nb::const_),
         "Return true if the attribute (by index) has been explicitly set; else false.")
    .def("hasAttributeValue",
         nb::overload_cast<const std::string&>(&View::hasAttributeValue, nb::const_),
         "Return true if the attribute (by name) has been explicitly set; else false.")
    .def("hasAttributeValue",
         nb::overload_cast<const Attribute*>(&View::hasAttributeValue, nb::const_),
         nb::arg("attr").none(),
         "Return true if the attribute (by pointer) has been explicitly set; else false.")

    .def("setAttributeToDefault",
         nb::overload_cast<IndexType>(&View::setAttributeToDefault),
         "Set Attribute (by index) to its default value")
    .def("setAttributeToDefault",
         nb::overload_cast<const std::string&>(&View::setAttributeToDefault),
         "Set Attribute (by name) to its default value")
    .def("setAttributeToDefault",
         nb::overload_cast<const Attribute*>(&View::setAttributeToDefault),
         nb::arg("attr").none(),
         "Set Attribute (by pointer) to its default value")

    // Scalar setters for int and python float (C++ double)
    .def(
      "setAttributeScalar",
      [](View& self, IndexType idx, int value) { return self.setAttributeScalar(idx, value); },
      "Set Attribute (by index) to int value")
    .def(
      "setAttributeScalar",
      [](View& self, IndexType idx, double value) { return self.setAttributeScalar(idx, value); },
      "Set Attribute (by index) to float (C++ double) value")
    .def(
      "setAttributeScalar",
      [](View& self, const std::string& name, int value) {
        return self.setAttributeScalar(name, value);
      },
      "Set Attribute (by name) to int value")
    .def(
      "setAttributeScalar",
      [](View& self, const std::string& name, double value) {
        return self.setAttributeScalar(name, value);
      },
      "Set Attribute (by name) to float (C++ double) value")
    .def(
      "setAttributeScalar",
      [](View& self, const Attribute* attr, int value) {
        return self.setAttributeScalar(attr, value);
      },
      "Set Attribute (by pointer) to int value")
    .def(
      "setAttributeScalar",
      [](View& self, const Attribute* attr, double value) {
        return self.setAttributeScalar(attr, value);
      },
      "Set Attribute (by pointer) to float (C++ double) value")

    // String setters
    .def("setAttributeString",
         nb::overload_cast<IndexType, const std::string&>(&View::setAttributeString),
         "Set Attribute (by index) to string value")
    .def("setAttributeString",
         nb::overload_cast<const std::string&, const std::string&>(&View::setAttributeString),
         "Set Attribute (by name) to string value")
    .def("setAttributeString",
         nb::overload_cast<const Attribute*, const std::string&>(&View::setAttributeString),
         "Set Attribute (by pointer) to string value")

    // Requires conduit::Node information
    // Scalar getters (Node::ConstValue version)
    // .def("getAttributeScalar",
    //      nb::overload_cast<IndexType>(&View::getAttributeScalar, nb::const_),
    //      "Return scalar Attribute value (by index) as Node::ConstValue")
    // .def("getAttributeScalar",
    //      nb::overload_cast<const std::string&>(&View::getAttributeScalar, nb::const_),
    //      "Return scalar Attribute value (by name) as Node::ConstValue")
    // .def("getAttributeScalar",
    //      nb::overload_cast<const Attribute*>(&View::getAttributeScalar, nb::const_),
    //      "Return scalar Attribute value (by pointer) as Node::ConstValue")

    // Scalar getters (templated, explicit for int and float)
    .def(
      "getAttributeScalarInt",
      [](View& self, IndexType idx) { return self.getAttributeScalar<int>(idx); },
      "Return scalar Attribute value (by index) as int")
    .def(
      "getAttributeScalarFloat",
      [](View& self, IndexType idx) { return self.getAttributeScalar<double>(idx); },
      "Return scalar Attribute value (by index) as float (C++ double)")
    .def(
      "getAttributeScalarInt",
      [](View& self, const std::string& name) { return self.getAttributeScalar<int>(name); },
      "Return scalar Attribute value (by name) as int")
    .def(
      "getAttributeScalarFloat",
      [](View& self, const std::string& name) { return self.getAttributeScalar<double>(name); },
      "Return scalar Attribute value (by name) as float (C++ double)")
    .def(
      "getAttributeScalarInt",
      [](View& self, const Attribute* attr) { return self.getAttributeScalar<int>(attr); },
      nb::arg("attr").none(),
      "Return scalar Attribute value (by pointer) as int")
    .def(
      "getAttributeScalarFloat",
      [](View& self, const Attribute* attr) { return self.getAttributeScalar<double>(attr); },
      nb::arg("attr").none(),
      "Return scalar Attribute value (by pointer) as float (C++ double)")

    // String getters
    .def("getAttributeString",
         nb::overload_cast<IndexType>(&View::getAttributeString, nb::const_),
         "Return string Attribute value (by index)")
    .def("getAttributeString",
         nb::overload_cast<const std::string&>(&View::getAttributeString, nb::const_),
         "Return string Attribute value (by name)")
    .def("getAttributeString",
         nb::overload_cast<const Attribute*>(&View::getAttributeString, nb::const_),
         "Return string Attribute value (by pointer)")

    // Requires conduit::Node information
    // Node reference getters
    // .def("getAttributeNodeRef",
    //      nb::overload_cast<IndexType>(&View::getAttributeNodeRef),
    //      nb::rv_policy::reference,
    //      "Return reference to Attribute Node (by index)")
    // .def("getAttributeNodeRef",
    //      nb::overload_cast<const std::string&>(&View::getAttributeNodeRef),
    //      nb::rv_policy::reference,
    //      "Return reference to Attribute Node (by name)")
    // .def("getAttributeNodeRef",
    //      nb::overload_cast<const Attribute*>(&View::getAttributeNodeRef),
    //      nb::rv_policy::reference,
    //      "Return reference to Attribute Node (by pointer)")

    // Attribute index iteration
    .def("getFirstValidAttrValueIndex",
         &View::getFirstValidAttrValueIndex,
         "Return first valid Attribute index for a set Attribute in View object"
         "(i.e., smallest index over all Attributes)")
    .def("getNextValidAttrValueIndex",
         &View::getNextValidAttrValueIndex,
         "Return next valid Attribute index for a set Attribute in View object after given index"
         "(i.e., smallest index over all Attribute indices larger than given one)");

  // Bindings for the Group class
  nb::class_<Group>(m_sidre, "Group")
    .def("getIndex", &Group::getIndex, "Return index of Group object within parent Group.")
    .def("getName", &Group::getName, "Return const reference to name of Group object.")
    .def("getPath", &Group::getPath, "Return path of Group object, not including its name.")
    .def("getPathName", &Group::getPathName, "Return full path of Group object, including its name.")
    .def("getParent",
         nb::overload_cast<>(&Group::getParent, nb::const_),
         nb::rv_policy::reference,
         "Return pointer to non-const parent Group of a Group.")
    .def("getNumGroups", &Group::getNumGroups, "Return number of child Groups in a Group object.")
    .def("getNumViews", &Group::getNumViews, "Return number of Views owned by a Group object.")
    .def("getDataStore",
         nb::overload_cast<>(&Group::getDataStore, nb::const_),
         nb::rv_policy::reference,
         "Return pointer to non-const DataStore object that owns this object.")

    .def("hasView",
         nb::overload_cast<const std::string&>(&Group::hasView, nb::const_),
         "Return true if Group includes a descendant View with given name or path; else false.")
    .def("hasView",
         nb::overload_cast<IndexType>(&Group::hasView, nb::const_),
         "Return true if this Group owns a View with given index; else false")
    .def("hasChildView",
         &Group::hasChildView,
         "Return true if this Group owns a View with given name (not path); else false.")
    .def("getViewIndex",
         &Group::getViewIndex,
         "Return index of View with given name owned by this Group object.")
    .def("getViewName",
         &Group::getViewName,
         "Return name of View with given index owned by Group object.")

    .def("getView",
         nb::overload_cast<const std::string&>(&Group::getView, nb::const_),
         nb::rv_policy::reference,
         "Return pointer to const View with given name or path.")
    .def("getView",
         nb::overload_cast<IndexType>(&Group::getView, nb::const_),
         nb::rv_policy::reference,
         "Return pointer to non-const View with given index.")
    .def("getFirstValidViewIndex",
         &Group::getFirstValidViewIndex,
         "Return first valid View index in Group object.")
    .def("getNextValidViewIndex",
         &Group::getNextValidViewIndex,
         "Return next valid View index in Group object after given index.")

    .def("createView",
         nb::overload_cast<const std::string&>(&Group::createView),
         nb::rv_policy::reference,
         "Create an undescribed (i.e., empty) View object with given name or path in this Group.")
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType>(&Group::createView),
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group that has a data description "
         "with data type and number of elements.")
    .def(
      "createViewWithShape",
      [](Group& self, const std::string& path, TypeID type, int ndims, const nb::ndarray<IndexType>& shape) {
        return self.createViewWithShape(path, type, ndims, shape.data());
      },
      nb::rv_policy::reference,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape.")
    .def("createView",
         nb::overload_cast<const std::string&, Buffer*>(&Group::createView),
         nb::rv_policy::reference,
         "Create an undescribed View object with given name or path in this Group and attach given "
         "Buffer to it.")
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType, Buffer*>(&Group::createView),
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group that has a data description "
         "with data type and number of elements and attach given Buffer to it.")
    .def(
      "createViewWithShape",
      [](Group& self,
         const std::string& path,
         TypeID type,
         int ndims,
         const nb::ndarray<IndexType>& shape,
         Buffer* buffer) {
        return self.createViewWithShape(path, type, ndims, shape.data(), buffer);
      },
      nb::rv_policy::reference,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and attach given Buffer to it.")

    .def(
      "createView",
      [](Group& self, const std::string& path, const nb::ndarray<>& a) {
        return self.createView(path, a.data());
      },
      nb::rv_policy::reference)

    .def(
      "createView",
      [](Group& self, const std::string& path, TypeID id, IndexType num_elems, const nb::ndarray<>& a) {
        return self.createView(path, id, num_elems, a.data());
      },
      nb::rv_policy::reference,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and number of elements and attach externally-owned data to it.")

    .def(
      "createViewWithShape",
      [](Group& self,
         const std::string& path,
         TypeID type,
         int ndims,
         const nb::ndarray<IndexType>& shape,
         const nb::ndarray<>& external_ptr) {
        return self.createViewWithShape(path, type, ndims, shape.data(), external_ptr.data());
      },
      nb::rv_policy::reference,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and attach externally-owned data (numpy array) to it.")
    .def("createViewAndAllocate",
         nb::overload_cast<const std::string&, TypeID, IndexType, int>(&Group::createViewAndAllocate),
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group that has a data description "
         "with data type and number of elements and allocate data for it.",
         nb::arg("path"),
         nb::arg("type"),
         nb::arg("num_elems"),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def(
      "createViewWithShapeAndAllocate",
      [](Group& self, const std::string& path, TypeID type, int ndims, const std::vector<IndexType>& shape) {
        return self.createViewWithShapeAndAllocate(path, type, ndims, shape.data());
      },
      nb::rv_policy::reference,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and allocate data for it.")

    .def("createViewScalar",
         &Group::createViewScalar<int>,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (int).",
         nb::arg("path"),
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("createViewScalar",
         &Group::createViewScalar<double>,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (C++ double, python float).",
         nb::arg("path"),
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("createViewString",
         &Group::createViewString,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given string.",
         nb::arg("path"),
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)

    .def("destroyView",
         nb::overload_cast<const std::string&>(&Group::destroyView),
         "Destroy View with given name or path owned by this Group, but leave its data intact.")
    .def("destroyViewAndData",
         nb::overload_cast<const std::string&>(&Group::destroyViewAndData),
         "Destroy View with given name or path owned by this Group and deallocate")
    .def("destroyViewAndData",
         nb::overload_cast<IndexType>(&Group::destroyViewAndData),
         "Destroy View with given index owned by this Group and deallocate its data if it's the "
         "only View associated with that data.")
    .def("destroyViewsAndData",
         &Group::destroyViewsAndData,
         "Destroy all Views owned by this Group and deallocate "
         "data for each View when it's the only View associated with that data.")

    .def("moveView",
         &Group::moveView,
         nb::rv_policy::reference,
         "Remove given View object from its owning Group and move it to this Group.")
    .def("copyView",
         &Group::copyView,
         nb::rv_policy::reference,
         "Create a (shallow) copy of given View object and add it to this Group.")

    .def("hasGroup",
         nb::overload_cast<const std::string&>(&Group::hasGroup, nb::const_),
         "Return true if this Group has a descendant Group with given name or path; else false.")
    .def("hasGroup",
         nb::overload_cast<IndexType>(&Group::hasGroup, nb::const_),
         "Return true if Group has an immediate child Group with given index; else false.")
    .def("hasChildGroup",
         &Group::hasChildGroup,
         "Return true if this Group has a child Group with given name; else false.")
    .def("getGroupIndex",
         &Group::getGroupIndex,
         "Return the index of immediate child Group with given name.")
    .def("getGroupName",
         &Group::getGroupName,
         "Return the name of immediate child Group with given index.")
    .def("getGroup",
         nb::overload_cast<const std::string&>(&Group::getGroup),
         nb::rv_policy::reference,
         "Return pointer to non-const child Group with given name or path.")
    .def("getGroup",
         nb::overload_cast<IndexType>(&Group::getGroup),
         nb::rv_policy::reference,
         "Return pointer to non-const immediate child Group with given index.")
    .def("views",
         nb::overload_cast<>(&Group::views),
         nb::rv_policy::reference,
         "Return an iterator over Views")
    .def("groups",
         nb::overload_cast<>(&Group::groups),
         nb::rv_policy::reference,
         "Return an iterator over Groups")
    .def("getFirstValidGroupIndex",
         &Group::getFirstValidGroupIndex,
         "Return first valid child Group index (i.e., smallest index over all child Groups).")
    .def("getNextValidGroupIndex",
         &Group::getNextValidGroupIndex,
         "Return next valid child Group index after given index.")
    .def("createGroup",
         &Group::createGroup,
         nb::rv_policy::reference,
         "Create a child Group within this Group with given name or path.",
         nb::arg("path"),
         nb::arg("is_list") = false,
         nb::arg("accept_existing") = false)
    .def("createUnnamedGroup",
         &Group::createUnnamedGroup,
         nb::rv_policy::reference,
         "Create a child Group within this Group with no name.",
         nb::arg("is_list") = false)
    .def("destroyGroup",
         nb::overload_cast<const std::string&>(&Group::destroyGroup),
         "Destroy child Group in this Group with given name or path.")
    .def("destroyGroup",
         nb::overload_cast<IndexType>(&Group::destroyGroup),
         "Destroy child Group within this Group with given index.")
    .def("destroyGroupAndData",
         nb::overload_cast<const std::string&>(&Group::destroyGroupAndData),
         "Destroy child Group at the given path, and destroy data that is "
         "not shared elsewhere.")
    .def("destroyGroupAndData",
         nb::overload_cast<IndexType>(&Group::destroyGroupAndData),
         "Destroy child Group with the given index, and destroy data that "
         "is not shared elsewhere.")
    .def("destroyGroupsAndData",
         &Group::destroyGroupsAndData,
         "Destroy all child Groups held by this Group, and destroy data that "
         "is not shared elsewhere.")
    .def("destroyGroupSubtreeAndData",
         &Group::destroyGroupSubtreeAndData,
         "Destroy the entire subtree of Groups and Views held by this Group, "
         "and destroy data that is not shared elsewhere.")
    .def("destroyGroups", &Group::destroyGroups, "Destroy all child Groups in this Group.")
    .def("moveGroup",
         &Group::moveGroup,
         "Remove given Group object from its parent Group and make it a child of this Group.")
    .def("copyGroup",
         &Group::copyGroup,
         nb::rv_policy::reference,
         "Create a (shallow) copy of Group hierarchy rooted at given "
         "Group and make the copy a child of this Group.")
    .def("deepCopyGroup",
         &Group::deepCopyGroup,
         nb::rv_policy::reference,
         "Create a deep copy of Group hierarchy rooted at given Group and "
         "make the copy a child of this Group.",
         nb::arg("srcGroup"),
         nb::arg("arrayAllocID") = INVALID_ALLOCATOR_ID,
         nb::arg("tupleAllocID") = INVALID_ALLOCATOR_ID)

    .def("print",
         nb::overload_cast<>(&Group::print, nb::const_),
         "Print JSON description of data Group to stdout.")
    .def("isEquivalentTo",
         &Group::isEquivalentTo,
         "Return true if this Group is equivalent to given Group; else false.",
         nb::arg("other"),
         nb::arg("checkName") = true)
    .def("isUsingMap", &Group::isUsingMap, "Return true if this Group holds items in map format.")
    .def("isUsingList", &Group::isUsingList, "Return true if this Group holds items in list format.")
    .def("save",
         nb::overload_cast<const std::string&, const std::string&, const Attribute*>(&Group::save,
                                                                                     nb::const_),
         "Save the Group to a file.",
         nb::arg("path"),
         nb::arg("protocol") = Group::getDefaultIOProtocol(),
         nb::arg("attr") = nullptr)
    .def("load",
         nb::overload_cast<const std::string&, const std::string&, bool>(&Group::load),
         "Load a Group hierarchy from a file into this Group",
         nb::arg("path"),
         nb::arg("protocol"),
         nb::arg("preserve_contents") = false)

    .def("loadExternalData",
         nb::overload_cast<const std::string&>(&Group::loadExternalData),
         "Load data into the Group's external views from a file.")
    .def("rename", &Group::rename, "Change the name of this Group.");

  // Bindings for the Attribute class
  nb::class_<Attribute>(m_sidre, "Attribute")
    .def("getName", &Attribute::getName, "Return the name of the Attribute object.")
    .def("getIndex", &Attribute::getIndex, "Return the unique index of this Attribute object.")
    .def("setDefaultScalar",
         &Attribute::setDefaultScalar<int>,
         "Set default value of Attribute as int. Return true if successfully changed.",
         nb::arg("value").noconvert())
    .def(
      "setDefaultScalar",
      &Attribute::setDefaultScalar<double>,
      "Set default value of Attribute as float (C++ double). Return true if successfully changed.",
      nb::arg("value").noconvert())
    .def("setDefaultString",
         &Attribute::setDefaultString,
         "Set default value of Attribute as string. Return true if successfully changed.")

    // Requires conduit::Node information
    // .def("setDefaultNodeRef",
    //      &Attribute::setDefaultNodeRef,
    //      "Set default value of Attribute as a Node reference.")
    .def("getDefaultNodeRef",
         &Attribute::getDefaultNodeRef,
         nb::rv_policy::reference,
         "Return default value of Attribute as Node reference.")
    .def("getTypeID", &Attribute::getTypeID, "Return type of Attribute.");
}

} /* end namespace sidre */
} /* end namespace axom */
