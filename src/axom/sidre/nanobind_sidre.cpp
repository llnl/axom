// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <nanobind/ndarray.h>

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

// Helper to map TypeID to dtype
nb::dlpack::dtype typeIDToDtype(DataTypeId id)
{
  switch(id)
  {
  case INT32_ID:
    return nb::dtype<int>();
  case INT64_ID:
    return nb::dtype<int64_t>();
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

  int shapeOutput[DMAX];
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
  nb::capsule owner(data, [](void*) noexcept {});

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
  nb::capsule owner(data, [](void*) noexcept {});

  DataTypeId id = self.getTypeID();

  return nb::ndarray<nb::numpy>(
    /* data = */ data,
    /* ndim = */ 1,
    /* shape = */ shape,
    /* owner = */ owner,
    /* strides = */ nullptr,
    /* dtype = */ typeIDToDtype(id));
}

// TODO determine return value policy (nb::rv_policy::reference or not) for ownership between C++ and Python:
// https://nanobind.readthedocs.io/en/latest/ownership.html#return-value-policies
NB_MODULE(pysidre, m_sidre)
{
  m_sidre.doc() = "A python extension for Axom's Sidre component";

  m_sidre.attr("InvalidIndex") = axom::InvalidIndex;
  m_sidre.attr("InvalidName") = axom::utilities::string::InvalidName;

  m_sidre.def("nameIsValid", &nameIsValid, "Returns true if name is valid, else false.");

#if defined(AXOM_USE_HDF5)
  m_sidre.attr("AXOM_USE_HDF5") = true;
#else
  m_sidre.attr("AXOM_USE_HDF5") = false;
#endif

  // Expose IndexType as an alias (bad cast error...)
  // #if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
  // m_sidre.attr("IndexType") = nb::type<std::int64_t>();
  // #else
  // m_sidre.attr("IndexType") = nb::type<int>();
  // #endif

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

  // Bindings for the DataStore class
  nb::class_<DataStore>(m_sidre, "DataStore")
    .def(nb::init<>())
    .def("getRoot",
         nb::overload_cast<>(&DataStore::getRoot),
         nb::rv_policy::reference,
         "Return pointer to the root Group")
    .def("getNumBuffers", &DataStore::getNumBuffers, "Return number of Buffers in the DataStore")
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
         nb::overload_cast<IndexType>(&DataStore::destroyBuffer),
         "Destroy a Buffer object by index")

    .def("generateBlueprintIndex",
         nb::overload_cast<const std::string&, const std::string&, const std::string&, int>(
           &DataStore::generateBlueprintIndex),
         "Generate a Conduit Blueprint index based on a mesh in stored in this DataStore.")

#ifdef AXOM_USE_MPI
    .def("generateBlueprintIndex",
         nb::overload_cast<MPI_Comm, const std::string&, const std::string&, const std::string&>(
           &DataStore::generateBlueprintIndex),
         "Generate a Conduit Blueprint index from a distributed mesh stored in this Datastore")
#endif
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
    // .def("getShape", &View::getShape, "Return the shape of the View's data.")
    .def(
      "getShape",
      [](View& self, int ndims, nb::ndarray<int>& shape) {
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
         nb::arg("offset") = 0,
         nb::arg("stride") = 1)
    .def(
      "apply",
      // nb::overload_cast<TypeID, int, const IndexType*>(&View::apply),
      [](View& self, TypeID type, int ndims, nb::ndarray<int64_t>& shape) {
        std::vector<int> temp(ndims);
        for(int i = 0; i < ndims; i++)
        {
          temp[i] = static_cast<int>(shape.data()[i]);
        }
        return self.apply(type, ndims, temp.data());
      },
      nb::rv_policy::reference,
      "Apply data description with type and numpy shape.")
    .def("setScalar",
         &View::setScalar<int>,
         nb::rv_policy::reference,
         "Set the View to hold a scalar value (int).")
    .def("setScalar",
         &View::setScalar<long>,
         nb::rv_policy::reference,
         "Set the View to hold a scalar value (long).")
    .def("setScalar",
         &View::setScalar<float>,
         nb::rv_policy::reference,
         "Set the View to hold a scalar value (float).")
    .def("setScalar",
         &View::setScalar<double>,
         nb::rv_policy::reference,
         "Set the View to hold a scalar value (double).")

    .def("setString", &View::setString, "Set the View to hold a string value.")
    .def("setExternalDataPtr",
         nb::overload_cast<void*>(&View::setExternalDataPtr),
         "Set the View to hold external data.")
    .def("setExternalDataPtr",
         nb::overload_cast<TypeID, IndexType, void*>(&View::setExternalDataPtr),
         "Set the View to hold described external data.")
    .def("setExternalDataPtr",
         nb::overload_cast<TypeID, int, const IndexType*, void*>(&View::setExternalDataPtr),
         "Set the View to hold described external data.")

    .def("getString",
         &View::getString,
         nb::rv_policy::reference,
         "Return the string contained in the View.")
    .def("getDataArray", &viewToNumpyArray, "Return the data held by the View as a numpy array.")

    .def("getData",
         &View::getData<int>,
         nb::rv_policy::reference,
         "Return the data held by the View (int).")
    .def("getData",
         &View::getData<long>,
         nb::rv_policy::reference,
         "Return the data held by the View (long).")
    .def("getData",
         &View::getData<float>,
         nb::rv_policy::reference,
         "Return the data held by the View (float).")
    .def("getData",
         &View::getData<double>,
         nb::rv_policy::reference,
         "Return the data held by the View (double).")

    // .def("getVoidPtr",
    //      &View::getVoidPtr,
    //      nb::rv_policy::reference,
    //      "Return a void pointer to the View's data.")
    .def("print",
         nb::overload_cast<>(&View::print, nb::const_),
         "Print JSON description of the View.")
    .def("rename", &View::rename, "Change the name of the View.");

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
    .def("createViewWithShape",
         nb::overload_cast<const std::string&, TypeID, int, const IndexType*>(
           &Group::createViewWithShape),
         "Create View object with given name or path in this Group that has a data description "
         "with data type and shape.")
    .def("createView",
         nb::overload_cast<const std::string&, Buffer*>(&Group::createView),
         nb::rv_policy::reference,
         "Create an undescribed View object with given name or path in this Group and attach given "
         "Buffer to it.")
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType, Buffer*>(&Group::createView),
         "Create View object with given name or path in this Group that has a data description "
         "with data type and number of elements and attach given Buffer to it.")
    .def("createViewWithShape",
         nb::overload_cast<const std::string&, TypeID, int, const IndexType*, Buffer*>(
           &Group::createViewWithShape),
         "Create View object with given name or path in this Group that has a data description "
         "with data type and shape and attach given Buffer to it.")

    // TODO
    // .def("createView",
    //      nb::overload_cast<const std::string&, void*>(&Group::createView),
    //      "Create View object with given name with given name or path in this Group and attach "
    //      "external data ptr to it.")

    // Debugging, seeing if can access data
    .def(
      "createView",
      // [](Group& self, const std::string& path, nb::ndarray<nb::c_contig, nb::device::cpu> a) {
      //   if(a.dtype() == nb::dtype<int64_t>() && a.ndim() == 1)
      //   {
      //     auto v = a.view<int64_t, nb::ndim<1>>();  // <-- new!

      //     for(size_t i = 0; i < v.shape(0); ++i)
      //     {
      //       printf("v(%ld) is %ld\n", i, v(i));
      //     }
      //   }
      //   else if(a.dtype() == nb::dtype<double>() && a.ndim() == 1)
      //   {
      //     auto v = a.view<double, nb::ndim<1>>();  // <-- new!

      //     for(size_t i = 0; i < v.shape(0); ++i)
      //     {
      //       printf("v(%ld) is %f\n", i, v(i));
      //     }
      //   }
      //   printf("Array dtype: int16=%i, uint32=%i, float32=%i, double=%i, int64=%i\n",
      //          a.dtype() == nb::dtype<int16_t>(),
      //          a.dtype() == nb::dtype<uint32_t>(),
      //          a.dtype() == nb::dtype<float>(),
      //          a.dtype() == nb::dtype<double>(),
      //          a.dtype() == nb::dtype<int64_t>());
      //   printf("Dtype components are: code=%d, bits=%d, lanes=%d\n",
      //          a.dtype().code,
      //          a.dtype().bits,
      //          a.dtype().lanes);
      //   printf("DIM IS %ld\n", a.ndim());
      //   return self.createView(path, a.data());
      // },

      [](Group& self, const std::string& path, const nb::ndarray<>& a) {
        return self.createView(path, a.data());
      },
      nb::rv_policy::reference)

    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType, void*>(&Group::createView),
         "Create View object with given name or path in this Group that has a data description "
         "with data type and number of elements and attach externally-owned data to it.")
    .def("createViewWithShape",
         nb::overload_cast<const std::string&, TypeID, int, const IndexType*, void*>(
           &Group::createViewWithShape),
         "Create View object with given name or path in this Group that has a data description "
         "with data type and shape and attach externally-owned data to it.")
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
      [](Group& self, const std::string& path, TypeID type, int ndims, const std::vector<int>& shape) {
        return self.createViewWithShapeAndAllocate(path, type, ndims, shape.data());
      },
      nb::rv_policy::reference,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and allocate data for it.")

    .def("createViewScalar",
         &Group::createViewScalar<int>,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (int).")
    .def("createViewScalar",
         &Group::createViewScalar<long>,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (long).")
    .def("createViewScalar",
         &Group::createViewScalar<float>,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (float).")
    .def("createViewScalar",
         &Group::createViewScalar<double>,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (double).")
    .def("createViewString",
         &Group::createViewString,
         nb::rv_policy::reference,
         "Create View object with given name or path in this Group set its data to given string.")

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

    .def("moveView",
         &Group::moveView,
         "Remove given View object from its owning Group and move it to this Group.")
    .def("copyView",
         &Group::copyView,
         nb::rv_policy::reference,
         "Create a (shallow) copy of given View object and add it to this Group.")

    .def("hasGroup",
         nb::overload_cast<const std::string&>(&Group::hasGroup, nb::const_),
         "Return true if this Group has a descendant Group with given name or path; else false.")
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
         nb::arg("is_list") = false)
    .def("destroyGroup",
         nb::overload_cast<const std::string&>(&Group::destroyGroup),
         "Destroy child Group in this Group with given name or path.")
    .def("destroyGroup",
         nb::overload_cast<IndexType>(&Group::destroyGroup),
         "Destroy child Group within this Group with given index.")
    .def("moveGroup",
         &Group::moveGroup,
         "Remove given Group object from its parent Group and make it a child of this Group.")

    .def("print",
         nb::overload_cast<>(&Group::print, nb::const_),
         "Print JSON description of data Group to stdout.")
    .def("isEquivalentTo",
         &Group::isEquivalentTo,
         "Return true if this Group is equivalent to given Group; else false.",
         nb::arg("other"),
         nb::arg("checkName") = true)
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
}

} /* end namespace sidre */
} /* end namespace axom */
