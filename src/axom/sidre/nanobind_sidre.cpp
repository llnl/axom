// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

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

// TODO determine return value policy (nb::rv_policy::reference or not) for ownership between C++ and Python:
// https://nanobind.readthedocs.io/en/latest/ownership.html#return-value-policies
NB_MODULE(pysidre, m_sidre)
{
  m_sidre.doc() = "A python extension for Axom's Sidre component";

  // Bindings for the DataStore class
  nb::class_<DataStore>(m_sidre, "DataStore")
    .def(nb::init<>())  // Bind the default constructor
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
    .def("getVoidPtr", &Buffer::getVoidPtr, "Return void-pointer to data held by Buffer.")
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
    .def("getShape", &View::getShape, "Return the shape of the View's data.")
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
         "Apply data description with number of elements, offset, and stride.")
    .def("apply",
         nb::overload_cast<TypeID, IndexType, IndexType, IndexType>(&View::apply),
         nb::rv_policy::reference,
         "Apply data description with type, number of elements, offset, and stride.")
    .def("apply",
         nb::overload_cast<TypeID, int, const IndexType*>(&View::apply),
         nb::rv_policy::reference,
         "Apply data description with type and shape.")
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
    .def("getData", &View::getData<int>, "Return the data held by the View (int).")
    .def("getData", &View::getData<long>, "Return the data held by the View (long).")
    .def("getData", &View::getData<float>, "Return the data held by the View (float).")
    .def("getData", &View::getData<double>, "Return the data held by the View (double).")

    .def("getVoidPtr",
         &View::getVoidPtr,
         nb::rv_policy::reference,
         "Return a void pointer to the View's data.")
    .def("print",
         nb::overload_cast<>(&View::print, nb::const_),
         "Print JSON description of the View.")
    .def("rename", &View::rename, "Change the name of the View.");

  // Bindings for the Group class
  nb::class_<Group>(m_sidre, "Group")
    .def("getIndex", &Group::getIndex)
    .def("getName", &Group::getName)
    .def("getPath", &Group::getPath)
    .def("getPathName", &Group::getPathName)
    .def("getParent", nb::overload_cast<>(&Group::getParent, nb::const_), nb::rv_policy::reference)
    .def("getNumGroups", &Group::getNumGroups)
    .def("getNumViews", &Group::getNumViews)
    .def("getDataStore", nb::overload_cast<>(&Group::getDataStore, nb::const_), nb::rv_policy::reference)

    .def("hasView", nb::overload_cast<const std::string&>(&Group::hasView, nb::const_))
    .def("hasChildView", &Group::hasChildView)
    .def("getViewIndex", &Group::getViewIndex)
    .def("getViewName", &Group::getViewName)

    .def("getView",
         nb::overload_cast<const std::string&>(&Group::getView, nb::const_),
         nb::rv_policy::reference)
    .def("getView", nb::overload_cast<IndexType>(&Group::getView, nb::const_), nb::rv_policy::reference)
    .def("getFirstValidViewIndex", &Group::getFirstValidViewIndex)
    .def("getNextValidViewIndex", &Group::getNextValidViewIndex)

    .def("createView", nb::overload_cast<const std::string&>(&Group::createView))
    .def("createView", nb::overload_cast<const std::string&, TypeID, IndexType>(&Group::createView))
    .def("createViewWithShape",
         nb::overload_cast<const std::string&, TypeID, int, const IndexType*>(
           &Group::createViewWithShape))
    .def("createView", nb::overload_cast<const std::string&, Buffer*>(&Group::createView))
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType, Buffer*>(&Group::createView))
    .def("createViewWithShape",
         nb::overload_cast<const std::string&, TypeID, int, const IndexType*, Buffer*>(
           &Group::createViewWithShape))
    .def("createView", nb::overload_cast<const std::string&, void*>(&Group::createView))
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType, void*>(&Group::createView))
    .def("createViewWithShape",
         nb::overload_cast<const std::string&, TypeID, int, const IndexType*, void*>(
           &Group::createViewWithShape))
    .def("createViewAndAllocate",
         nb::overload_cast<const std::string&, TypeID, IndexType, int>(&Group::createViewAndAllocate),
         "todo description",
         nb::arg("path"),
         nb::arg("type"),
         nb::arg("num_elems"),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("createViewWithShapeAndAllocate", &Group::createViewWithShapeAndAllocate)

    .def("createViewScalar",
         &Group::createViewScalar<int>,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (int).")
    .def("createViewScalar",
         &Group::createViewScalar<long>,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (long).")
    .def("createViewScalar",
         &Group::createViewScalar<float>,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (float).")
    .def("createViewScalar",
         &Group::createViewScalar<double>,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (double).")
    .def("createViewString", &Group::createViewString)

    .def("destroyView", nb::overload_cast<const std::string&>(&Group::destroyView))
    .def("destroyViewAndData", nb::overload_cast<const std::string&>(&Group::destroyViewAndData))
    .def("destroyViewAndData", nb::overload_cast<IndexType>(&Group::destroyViewAndData))

    .def("moveView", &Group::moveView)
    .def("copyView", &Group::copyView)

    .def("hasGroup", nb::overload_cast<const std::string&>(&Group::hasGroup, nb::const_))
    .def("hasChildGroup", &Group::hasChildGroup)
    .def("getGroupIndex", &Group::getGroupIndex)
    .def("getGroupName", &Group::getGroupName)
    .def("getGroup", nb::overload_cast<const std::string&>(&Group::getGroup), nb::rv_policy::reference)
    .def("getGroup", nb::overload_cast<IndexType>(&Group::getGroup), nb::rv_policy::reference)
    .def("getFirstValidGroupIndex", &Group::getFirstValidGroupIndex)
    .def("getNextValidGroupIndex", &Group::getNextValidGroupIndex)
    .def("createGroup", &Group::createGroup, nb::arg("path"), nb::arg("is_list") = false)
    .def("destroyGroup", nb::overload_cast<const std::string&>(&Group::destroyGroup))
    .def("destroyGroup", nb::overload_cast<IndexType>(&Group::destroyGroup))
    .def("moveGroup", &Group::moveGroup)

    .def("print", nb::overload_cast<>(&Group::print, nb::const_))
    .def("isEquivalentTo",
         &Group::isEquivalentTo,
         "todo description",
         nb::arg("other"),
         nb::arg("checkName") = true)
    .def("save",
         nb::overload_cast<const std::string&, const std::string&, const Attribute*>(&Group::save,
                                                                                     nb::const_))
    .def("load",
         nb::overload_cast<const std::string&, const std::string&, bool>(&Group::load),
         "TODO description",
         nb::arg("path"),
         nb::arg("protocol"),
         nb::arg("preserve_contents") = false)

    .def("loadExternalData", nb::overload_cast<const std::string&>(&Group::loadExternalData))
    .def("rename", &Group::rename);
}

} /* end namespace sidre */
} /* end namespace axom */
