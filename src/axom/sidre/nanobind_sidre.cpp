// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/optional.h>
#include <nanobind/ndarray.h>
#include <nanobind/make_iterator.h>

#include <memory>
#include <optional>
#include <unordered_map>

#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/slic/interface/slic.hpp"

#include "core/SidreTypes.hpp"
#include "core/Buffer.hpp"
#include "core/View.hpp"
#include "core/DataStore.hpp"
#include "core/Group.hpp"
#if defined(AXOM_USE_MPI)
  #include "spio/IOManager.hpp"
#endif

// Separate Conduit header for python functionality
#include "conduit_python.hpp"

namespace nb = nanobind;
using namespace nb::literals;

namespace axom
{
namespace sidre
{

// Helper to map TypeID to nanobind dtype
nb::dlpack::dtype typeIDToDtype(DataTypeId id)
{
  if(id == INT8_ID)
  {
    return nb::dtype<std::int8_t>();
  }
  if(id == INT16_ID)
  {
    return nb::dtype<std::int16_t>();
  }
  if(id == INT32_ID || id == INT_ID)
  {
    return nb::dtype<std::int32_t>();
  }
  if(id == INT64_ID)
  {
    return nb::dtype<std::int64_t>();
  }
  if(id == UINT8_ID)
  {
    return nb::dtype<std::uint8_t>();
  }
  if(id == UINT16_ID)
  {
    return nb::dtype<std::uint16_t>();
  }
  if(id == UINT32_ID || id == UINT_ID)
  {
    return nb::dtype<std::uint32_t>();
  }
  if(id == UINT64_ID)
  {
    return nb::dtype<std::uint64_t>();
  }
  if(id == FLOAT32_ID || id == FLOAT_ID)
  {
    return nb::dtype<float>();
  }
  if(id == FLOAT64_ID || id == DOUBLE_ID)
  {
    return nb::dtype<double>();
  }

  SLIC_ERROR("DataTypeId unsupported for numpy");
  return nb::dtype<double>();
}

/*!
 * \brief Returns a View as a numpy array.
 *
 * \note Max dimensions (DMAX) is currently set to 10.
 * \pre data description must have been applied.
 *
 * \warning The returned array is a zero-copy view into the View's storage; it
 *  does not own the data. The array keeps the Python \a View object (and, via
 *  the View's owner-chain, its Group and DataStore) alive for as long as the
 *  array is referenced. One sharp edge remains: \c View::reallocate (and
 *  \c Buffer::reallocate) can move the underlying storage while a live array
 *  still points at the old allocation, the same hazard as resizing a container
 *  under a numpy view. Re-acquire the array after any reallocation.
 */
nb::ndarray<nb::numpy> viewToNumpyArray(nb::handle_t<View> h)
{
  View& self = nb::cast<View&>(h);

  // Manually applying offset
  void* data = self.getVoidPtr();
  char* data_with_offset = static_cast<char*>(data) + (self.getOffset() * self.getBytesPerElement());
  data = static_cast<void*>(data_with_offset);

  constexpr int DMAX = 10;

  IndexType shapeOutput[DMAX];
  size_t ndims = self.getShape(DMAX, shapeOutput);

  // Guard against buffer overflow if View has more dimensions than our buffer can hold
  SLIC_ERROR_IF(ndims > DMAX,
                "View has " << ndims << " dimensions, exceeds maximum of " << DMAX
                            << ". Cannot convert to numpy array.");

  size_t shape[DMAX];
  for(size_t i = 0; i < ndims; i++)
  {
    shape[i] = static_cast<size_t>(shapeOutput[i]);
  }

  // When stride is not default of 1, guaranteed that shape is 1D.
  int64_t* strides = nullptr;
  int64_t stride_array[1];
  if(self.getStride() != 1)
  {
    stride_array[0] = static_cast<int64_t>(self.getStride());
    strides = stride_array;
  }

  DataTypeId id = self.getTypeID();

  // Pass the Python 'self' object as the array owner: nanobind ties the
  // array's lifetime to it, so the View (and its DataStore) cannot be freed
  // while the array is alive. This resolves the prior no-op owner capsule that
  // left the array dangling once the View was collected.
  return nb::ndarray<nb::numpy>(
    /* data = */ data,
    /* ndim = */ ndims,
    /* shape = */ shape,
    /* owner = */ h,
    /* strides = */ strides,
    /* dtype = */ typeIDToDtype(id));
}

/*!
 * \brief Returns a Buffer as a numpy array.
 *
 * \pre data description must have been applied.
 *
 * \warning The returned array is a zero-copy view into the Buffer's storage; it
 *  does not own the data. The array keeps the Python \a Buffer object (and its
 *  owning DataStore) alive for as long as the array is referenced. As with
 *  \c viewToNumpyArray, \c Buffer::reallocate can move storage under a live
 *  array; re-acquire the array after any reallocation.
 */
nb::ndarray<nb::numpy> bufferToNumpyArray(nb::handle_t<Buffer> h)
{
  Buffer& self = nb::cast<Buffer&>(h);

  void* data = self.getVoidPtr();

  size_t shape[1] = {static_cast<size_t>(self.getNumElements())};

  DataTypeId id = self.getTypeID();

  // Pass the Python 'self' object as the array owner (see viewToNumpyArray).
  return nb::ndarray<nb::numpy>(
    /* data = */ data,
    /* ndim = */ 1,
    /* shape = */ shape,
    /* owner = */ h,
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
          self.end(),
          // Pin each yielded element to the iterator (applies to __next__).
          // The iterator is in turn pinned to the collection by the
          // keep_alive on __iter__ below, and the collection accessor
          // (e.g. Group::views()) is reference_internal, so a harvested
          // element transitively keeps its DataStore alive.
          nb::keep_alive<0, 1>());
      },
      nb::keep_alive<0, 1>());
}

/*!
 * \brief Helper function that converts a C++ Conduit Node into a python nanobind object
 *        that is of type Node in python.
 *
 * \param [in] node C++ Conduit Node
 *
 * \return A python nanobind object
 */
nb::object nodeToNbObject(conduit::Node& node)
{
  // Setup conduit python c api
  if(import_conduit() < 0)
  {
    SLIC_ERROR("Failed to import Conduit Python C-API");
  }

  // 0 - python owns => false
  PyObject* wrapped = PyConduit_Node_Python_Wrap(&node, 0);

  SLIC_ERROR_IF(!wrapped, "PyObject is null");
  SLIC_ERROR_IF(!PyConduit_Node_Check(wrapped), "PyObject is not a Conduit Node");

  // Return nb::object
  return nb::steal<nb::object>(wrapped);
}

/*!
 * \brief Helper function that converts a python nanobind object that is of type Node
 *        in python into a C++ Conduit Node.
 *
 * \param [in] node C++ Conduit Node
 * \pre The python nanobind object must be of type Node in python
 *
 * \return A python nanobind object
 */
conduit::Node& nbObjectToNode(nb::object& o)
{
  // Setup conduit python C API
  if(import_conduit() < 0)
  {
    SLIC_ERROR("Failed to import Conduit Python C-API");
  }

  PyObject* obj = o.ptr();
  SLIC_ERROR_IF(!obj, "PyObject is null");
  SLIC_ERROR_IF(!PyConduit_Node_Check(obj), "PyObject is not a Conduit Node");

  // Turn python PyObject into C++ conduit::Node
  conduit::Node* cpp_node = PyConduit_Node_Get_Node_Ptr(obj);
  SLIC_ERROR_IF(!cpp_node, "Failed to get underlying conduit::Node pointer");

  return *cpp_node;
}

/*!
 * \brief Binding-side owner pinning for external NumPy-backed Sidre views.
 *
 * Sidre external views (View::setExternalDataPtr / Group::createView(..., void*))
 * borrow a raw pointer and do not own the storage. In Python, it is common to
 * pass a temporary NumPy array (e.g. createView("x", np.arange(...))) and then
 * drop it immediately, which can leave Sidre holding a dangling pointer once
 * the ndarray is garbage collected.
 *
 * To keep Sidre's C++ semantics unchanged while making the Python API safe,
 * we maintain a binding-only registry of "pins": copies of the nanobind ndarray
 * wrapper. Copying nb::ndarray increments the underlying ndarray owner's
 * refcount via nanobind's internal handle, so the NumPy storage stays alive for
 * as long as the pin exists. A pin therefore ties the external array's lifetime
 * to the *C++ View's* lifetime, not to any transient Python proxy: the array
 * survives even if the Python View object that set it is discarded, as long as
 * the View still lives in its DataStore (see the lifetime tests).
 *
 * **Per-DataStore scoping.** The registry is keyed by owning DataStore* and,
 * within each DataStore, by View*. This is what makes raw-pointer keys safe:
 *
 *  - When a DataStore's Python object is collected, nanobind destroys the C++
 *    DataStore (the user always holds a DataStore through Python, so the two
 *    lifetimes coincide). We install a weak reference on the DataStore at the
 *    first pin whose callback erases that DataStore's entire sub-map. Pins are
 *    thus released no later than DataStore destruction -- the registry never
 *    grows without bound, even for Views torn down by the C++ DataStore
 *    destructor rather than an explicit destroyView()/destroyGroup().
 *  - A View* is only meaningful within its owning DataStore, and that
 *    DataStore's sub-map is wiped when the DataStore dies. A View* address
 *    reused by a *different* DataStore therefore cannot collide with a stale
 *    pin, and a reused DataStore* address starts from a fresh (empty) sub-map.
 *  - We never look up a pin by a View* that might be stale: copyView/copyGroup
 *    only search for the source pin while the source View is live, then re-pin
 *    the destination from that live ndarray value.
 *
 * Pins are also released eagerly when the external pointer is cleared
 * (View.clear(), setExternalData(None)) and when views/groups are destroyed via
 * the bound destroy* APIs, so memory is reclaimed promptly in the common case
 * rather than waiting for DataStore destruction.
 *
 * \note Thread safety: all access goes through the GIL (see the module
 *  docstring). If the bindings ever release the GIL, this registry needs a mutex.
 *
 * \note Pin scoping assumes a pinned View stays within the DataStore it belonged
 *  to when it was pinned. Sidre reparenting (moveView/moveGroup) stays within a
 *  single DataStore, so a View's owning DataStore is stable for its lifetime and
 *  the DataStore* key never goes stale under a supported operation.
 *  If Sidre ever gained cross-DataStore migration of a live View, that View's pin
 *  would remain under its original DataStore (and be released when that DataStore is collected),
 *  so this invariant would need revisiting.
 */
struct DataStoreExternalPins
{
  std::unordered_map<View*, nb::ndarray<>> pins;
  // Holds the weak reference whose callback clears this sub-map; keeping it here
  // keeps the callback armed for the lifetime of the entries.
  nb::object datastore_weakref;
};

std::unordered_map<DataStore*, DataStoreExternalPins>& externalDataOwnerRegistry()
{
  // Intentionally heap-allocated so any Python-owned references it still holds
  // are not destroyed after interpreter finalization during static shutdown.
  static auto* registry = new std::unordered_map<DataStore*, DataStoreExternalPins>();
  return *registry;
}

//! Return the owning DataStore of a View, or nullptr if it has none yet.
DataStore* owningDataStore(View* view)
{
  if(view == nullptr)
  {
    return nullptr;
  }
  Group* group = view->getOwningGroup();
  return (group != nullptr) ? group->getDataStore() : nullptr;
}

//! Erase all pins recorded for \a ds (called when the DataStore is collected).
void releaseDataStoreExternalPins(DataStore* ds) { externalDataOwnerRegistry().erase(ds); }

/*!
 * \brief Record \a owner as the pin for \a view, scoped to its DataStore.
 *
 * On the first pin into a given DataStore, installs a weak reference on the
 * DataStore's Python object so the sub-map is cleared when the DataStore is
 * destroyed. Re-assigning a View*'s pin releases the previous ndarray wrapper.
 */
void pinExternalDataOwner(View* view, const nb::ndarray<>& owner)
{
  DataStore* ds = owningDataStore(view);
  // Enforce the precondition in debug builds;
  // release builds fall through to the null-safe early return below.
  SLIC_ASSERT_MSG(view == nullptr || ds != nullptr,
                  "pinExternalDataOwner: a non-null View is expected to have an owning DataStore");
  if(view == nullptr || ds == nullptr)
  {
    return;
  }

  DataStoreExternalPins& entry = externalDataOwnerRegistry()[ds];
  if(!entry.datastore_weakref.is_valid())
  {
    // Retrieve the DataStore's existing Python wrapper and attach a weakref
    // whose callback clears this DataStore's pins. nb::find returns a null
    // object if no wrapper exists; in that (unexpected) case we simply skip the
    // weakref -- the eager release paths still apply, and the worst case is the
    // pre-existing process-lifetime retention.
    nb::object ds_obj = nb::find(*ds);
    if(ds_obj.is_valid() && !ds_obj.is_none())
    {
      entry.datastore_weakref =
        nb::weakref(ds_obj, nb::cpp_function([ds](nb::handle) { releaseDataStoreExternalPins(ds); }));
    }
  }

  // Map assignment releases the previous ndarray wrapper if one was present.
  entry.pins[view] = nb::ndarray<>(owner);
}

void releaseExternalDataOwner(View* view)
{
  DataStore* ds = owningDataStore(view);
  if(view == nullptr || ds == nullptr)
  {
    return;
  }
  auto it = externalDataOwnerRegistry().find(ds);
  if(it != externalDataOwnerRegistry().end())
  {
    it->second.pins.erase(view);
  }
}

void releaseExternalDataOwners(Group* group)
{
  if(group == nullptr)
  {
    return;
  }

  for(auto& v : group->views())
  {
    releaseExternalDataOwner(&v);
  }

  for(auto& g : group->groups())
  {
    releaseExternalDataOwners(&g);
  }
}

void releaseExternalDataOwnersOfViews(Group& group)
{
  for(auto& v : group.views())
  {
    releaseExternalDataOwner(&v);
  }
}

/*!
 * \brief Copy the external-data pin from a source View to a destination View.
 *
 * copyView() makes a shallow copy that shares the external pointer, so the
 * destination needs its own pin to keep the NumPy array alive independently of
 * the source. The source View is live for the duration of the copy (the caller
 * holds it), so looking up its pin within its own DataStore's sub-map is safe;
 * we then pin the destination from that live ndarray value. We never search the
 * registry by a View* that could be stale.
 *
 * \param src_view Source View (live; expected to hold an external data pin)
 * \param dst_view Destination View (receives a copy of the pin)
 */
void copyExternalDataOwner(const View* src_view, View* dst_view)
{
  if(src_view == nullptr || dst_view == nullptr)
  {
    return;
  }

  DataStore* src_ds = owningDataStore(const_cast<View*>(src_view));
  if(src_ds == nullptr)
  {
    return;
  }

  auto& registry = externalDataOwnerRegistry();
  auto ds_it = registry.find(src_ds);
  if(ds_it == registry.end())
  {
    return;
  }

  auto pin_it = ds_it->second.pins.find(const_cast<View*>(src_view));
  if(pin_it != ds_it->second.pins.end())
  {
    // Re-pin the destination from the source's live ndarray value (scoped to
    // the destination's own DataStore by pinExternalDataOwner).
    pinExternalDataOwner(dst_view, pin_it->second);
  }
}

/*!
 * \brief Recursively copy external data pins from source Group hierarchy to destination.
 *
 * When copyGroup() creates a shallow copy of a Group hierarchy, all Views with
 * external data in the source hierarchy need their pins copied to the destination.
 *
 * \param src_group Source Group (root of hierarchy to copy from)
 * \param dst_group Destination Group (root of copied hierarchy)
 */
void copyExternalDataOwners(const Group* src_group, Group* dst_group)
{
  if(src_group == nullptr || dst_group == nullptr)
  {
    return;
  }

  // Copy pins for all Views in this Group
  for(auto& src_view : src_group->views())
  {
    if(src_view.isExternal())
    {
      View* dst_view = dst_group->getView(src_view.getName());
      if(dst_view != nullptr)
      {
        copyExternalDataOwner(&src_view, dst_view);
      }
    }
  }

  // Recursively copy pins for all child Groups
  for(auto& src_child : src_group->groups())
  {
    Group* dst_child = dst_group->getGroup(src_child.getName());
    if(dst_child != nullptr)
    {
      copyExternalDataOwners(&src_child, dst_child);
    }
  }
}

View* setExternalDataAndPinOwner(View& view, const nb::ndarray<>& owner)
{
  View* result = view.setExternalDataPtr(owner.data());
  if(result != nullptr && result->isExternal())
  {
    pinExternalDataOwner(result, owner);
  }
  return result;
}

View* setExternalDataAndPinOwner(View& view, TypeID type, IndexType num_elems, const nb::ndarray<>& owner)
{
  View* result = view.setExternalDataPtr(type, num_elems, owner.data());
  if(result != nullptr && result->isExternal())
  {
    pinExternalDataOwner(result, owner);
  }
  return result;
}

View* setExternalDataAndPinOwner(View& view,
                                 TypeID type,
                                 int ndims,
                                 const nb::ndarray<IndexType>& shape,
                                 const nb::ndarray<>& owner)
{
  View* result = view.setExternalDataPtr(type, ndims, shape.data(), owner.data());
  if(result != nullptr && result->isExternal())
  {
    pinExternalDataOwner(result, owner);
  }
  return result;
}

#if defined(AXOM_USE_MPI)
MPI_Comm mpiCommFromObject(nb::object comm)
{
  MPI_Comm c = MPI_COMM_WORLD;
  if(!comm.is_none())
  {
    MPI_Fint handle = nb::cast<MPI_Fint>(comm.attr("py2f")());
    c = MPI_Comm_f2c(handle);
  }
  return c;
}

/*!
 * \brief Python-facing IOManager holder that owns the communicator lifetime.
 *
 * sidre::IOManager stores the MPI_Comm passed to its constructor but does not
 * duplicate or free it. That borrowed-communicator contract works for C++
 * callers, but it is unsafe for mpi4py objects: py2f() exposes the object's
 * current communicator handle, and Python code may later drop or explicitly
 * Free() that object while axom.sidre.IOManager is still alive.
 *
 * PyIOManager keeps the public Python class name as axom.sidre.IOManager while
 * giving the binding its own lifetime boundary. It duplicates the input
 * communicator, constructs sidre::IOManager with that duplicate, destroys the
 * IOManager first, and then frees the duplicate when MPI is still active.
 */
class PyIOManager
{
public:
  PyIOManager(MPI_Comm comm, bool use_scr)
  {
    int err = MPI_Comm_dup(comm, &m_comm);
    SLIC_ERROR_IF(err != MPI_SUCCESS,
                  "Failed to duplicate MPI communicator for axom.sidre.IOManager");
    m_manager = std::make_unique<IOManager>(m_comm, use_scr);
  }

  ~PyIOManager()
  {
    // IOManager may still use m_comm during destruction, so release it before
    // freeing the duplicated communicator.
    m_manager.reset();

    int initialized = 0;
    MPI_Initialized(&initialized);
    if(initialized != 0 && m_comm != MPI_COMM_NULL)
    {
      int finalized = 0;
      MPI_Finalized(&finalized);
      if(finalized == 0)
      {
        MPI_Comm_free(&m_comm);
      }
    }
  }

  PyIOManager(const PyIOManager&) = delete;
  PyIOManager& operator=(const PyIOManager&) = delete;

  IOManager& manager() { return *m_manager; }
  const IOManager& manager() const { return *m_manager; }

private:
  MPI_Comm m_comm {MPI_COMM_NULL};
  std::unique_ptr<IOManager> m_manager;
};
#endif

// The extension installs as ``axom/sidre/_sidre.<tag>.so`` and is re-exported
// by the ``axom.sidre`` package (see src/python/src/axom/sidre/__init__.py).
NB_MODULE(_sidre, m_sidre)
{
  m_sidre.doc() = R"pbdoc(
    A python extension for Axom's Sidre component.

    **Thread Safety:**
    This module relies on Python's Global Interpreter Lock (GIL) for thread safety.

    - **Python threads:** Safe. The GIL serializes all operations.
    - **C++ threads:** Unsafe. Calling Sidre methods from C++ threads without acquiring
      the GIL will cause data races and undefined behavior.
    - **GIL release:** The bindings do not release the GIL during operations, so Python
      threads are always serialized. If future changes add `nb::call_guard<nb::gil_scoped_release>`,
      explicit synchronization (e.g., mutexes) must be added to the external data registry.

    **External Data Lifetime:**
    Views can reference external numpy arrays via setExternalData() or createView().
    The binding automatically pins these arrays to prevent garbage collection while
    the View exists, so an array stays valid even if the Python View object that set
    it is discarded (as long as the View still lives in its DataStore). Pins are
    scoped per-DataStore and are released when:
    - View.clear() or setExternalData(None) is called
    - The View is destroyed via destroyView() / destroyViewAndData()
    - The owning Group hierarchy is destroyed via destroyGroup*() methods
    - The owning DataStore is destroyed (a weak reference on the DataStore clears
      its remaining pins, so the registry never grows without bound -- even for
      Views torn down by the C++ DataStore destructor rather than an explicit
      destroy* call)

    **Reallocation Hazards:**
    Arrays obtained via getDataArray() are zero-copy views into Sidre storage.
    View::reallocate() or Buffer::reallocate() can move the underlying storage,
    leaving existing numpy arrays pointing to freed memory. Always re-acquire
    arrays after any reallocation. The binding cannot prevent this hazard without
    breaking zero-copy semantics.
  )pbdoc";

  // Module version mirrors the Axom release
  m_sidre.attr("__version__") = AXOM_VERSION_FULL;

  m_sidre.attr("InvalidIndex") = axom::InvalidIndex;
  m_sidre.attr("InvalidName") = axom::utilities::string::InvalidName;

  m_sidre.def("indexIsValid", &indexIsValid, "Returns true if idx is valid, else false.");
  m_sidre.def("nameIsValid", &nameIsValid, "Returns true if name is valid, else false.");

#if defined(AXOM_USE_HDF5)
  m_sidre.attr("AXOM_USE_HDF5") = true;
#else
  m_sidre.attr("AXOM_USE_HDF5") = false;
#endif

#if defined(AXOM_USE_MPI)
  m_sidre.attr("AXOM_ENABLE_MPI") = true;
#else
  m_sidre.attr("AXOM_ENABLE_MPI") = false;
#endif

#if defined(AXOM_USE_MPI) && \
  ((NB_VERSION_MAJOR > 2) || (NB_VERSION_MAJOR == 2 && NB_VERSION_MINOR >= 10))
  m_sidre.attr("AXOM_HAS_DISTRIBUTED_BLUEPRINT_INDEX_BINDING") = true;
#else
  m_sidre.attr("AXOM_HAS_DISTRIBUTED_BLUEPRINT_INDEX_BINDING") = false;
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
  // DataStore is weak-referenceable so the external-data registry can attach a
  // weakref callback that releases that DataStore's pins when it is destroyed
  // (see externalDataOwnerRegistry).
  nb::class_<DataStore>(m_sidre, "DataStore", nb::is_weak_referenceable())
    .def(nb::init<>())
    .def("getRoot",
         nb::overload_cast<>(&DataStore::getRoot),
         nb::rv_policy::reference_internal,
         "Return pointer to the root Group")
    .def("getNumBuffers", &DataStore::getNumBuffers, "Return number of Buffers in the DataStore")
    .def("hasBuffer",
         &DataStore::hasBuffer,
         "Return true if DataStore owns a Buffer with given index; else false")
    .def("getBuffer",
         &DataStore::getBuffer,
         nb::rv_policy::reference_internal,
         "Return pointer to Buffer object with the given index")

    .def("createBuffer",
         nb::overload_cast<>(&DataStore::createBuffer),
         nb::rv_policy::reference_internal,
         "Create an undescribed Buffer object")
    .def("createBuffer",
         nb::overload_cast<TypeID, IndexType>(&DataStore::createBuffer),
         nb::rv_policy::reference_internal,
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
         nb::keep_alive<0, 1>(),
         "Return an iterator over Buffers")

    .def("getNumAttributes",
         &DataStore::getNumAttributes,
         "Return number of Attributes in the DataStore")
    .def("createAttributeScalar",
         &DataStore::createAttributeScalar<int>,
         nb::rv_policy::reference_internal,
         "Create an Attribute object with a default int scalar value",
         nb::arg("name"),
         nb::arg("default_value").noconvert())
    .def("createAttributeScalar",
         &DataStore::createAttributeScalar<double>,
         nb::rv_policy::reference_internal,
         "Create an Attribute object with a default float (C++ double) scalar value",
         nb::arg("name"),
         nb::arg("default_value").noconvert())
    .def("createAttributeString",
         &DataStore::createAttributeString,
         nb::rv_policy::reference_internal,
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
         nb::rv_policy::reference_internal,
         "Return pointer to non-const Attribute with given index")
    .def("getAttribute",
         nb::overload_cast<const std::string&>(&DataStore::getAttribute),
         nb::rv_policy::reference_internal,
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
         nb::keep_alive<0, 1>(),
         "Return an iterator over Attributes")

#if defined(AXOM_USE_MPI) && \
  ((NB_VERSION_MAJOR > 2) || (NB_VERSION_MAJOR == 2 && NB_VERSION_MINOR >= 10))
    // Distributed Blueprint index generation builds cleanly under nanobind >= 2.10
    .def(
      "generateBlueprintIndex",
      [](DataStore& self,
         nb::object comm,
         const std::string& domain_path,
         const std::string& mesh_name,
         const std::string& index_path) {
        MPI_Comm c = MPI_COMM_WORLD;
        c = mpiCommFromObject(comm);
        return self.generateBlueprintIndex(c, domain_path, mesh_name, index_path);
      },
      "Generate a Conduit Blueprint index from a distributed mesh stored in this "
      "DataStore. Pass None to use MPI_COMM_WORLD.",
      nb::arg("comm").none(),
      nb::arg("domain_path"),
      nb::arg("mesh_name"),
      nb::arg("index_path"))
#endif

    .def("print",
         nb::overload_cast<>(&DataStore::print, nb::const_),
         "Print JSON description of the DataStore");

  // Bindings for the Buffer class
  nb::class_<Buffer>(m_sidre, "Buffer")
    .def("getIndex", &Buffer::getIndex, "Return the unique index of this Buffer object.")
    .def("getNumViews", &Buffer::getNumViews, "Return number of Views this Buffer is attached to.")
    // .def("getVoidPtr", &Buffer::getVoidPtr, "Return void-pointer to data held by Buffer.")
    .def("getDataArray",
         &bufferToNumpyArray,
         "Return the data held by the Buffer as a numpy array.\n\n"
         "The array is a zero-copy view into the Buffer's storage "
         "and keeps the Buffer (and its DataStore) alive while referenced. "
         "Buffer.reallocate() can move the storage, leaving a previously returned array "
         "pointing at freed memory; re-acquire the array after any reallocation.")
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
    .def("checksum",
         &View::checksum,
         "Return a checksum for the View's name, metadata, and data.",
         nb::arg("includeAttributes") = true)
    .def("getOwningGroup",
         nb::overload_cast<>(&View::getOwningGroup),
         nb::rv_policy::reference_internal,
         "Return the owning Group of the View.")
    .def("hasBuffer", &View::hasBuffer, "Check if the View has an associated Buffer object.")
    .def("getBuffer",
         nb::overload_cast<>(&View::getBuffer),
         nb::rv_policy::reference_internal,
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
         "Attach a Buffer object to the View.",
         nb::arg("buffer").none())
    .def("attachBuffer",
         nb::overload_cast<TypeID, IndexType, Buffer*>(&View::attachBuffer),
         nb::rv_policy::reference,
         "Describe the data view and attach Buffer object.",
         nb::arg("type"),
         nb::arg("num_elems"),
         nb::arg("buffer").none())
    .def("attachBuffer",
         nb::overload_cast<TypeID, int, const IndexType*, Buffer*>(&View::attachBuffer),
         nb::rv_policy::reference,
         "Describe the data view and attach Buffer object",
         nb::arg("type"),
         nb::arg("ndims"),
         nb::arg("shape"),
         nb::arg("buffer").none())

    .def(
      "clear",
      [](View& self) {
        self.clear();
        releaseExternalDataOwner(&self);
      },
      "Clear data and metadata from the View.")
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
      [](View& self, std::optional<nb::ndarray<>> external_ptr) {
        // A single undescribed-data overload covering both None
        // (clear the external pointer and release any pin) and a numpy array (set + pin).
        // Using std::optional<ndarray> lets nanobind reject a non-array argument
        // with a clean "incompatible function arguments" error
        // rather than throwing mid-body from an explicit cast.
        if(!external_ptr.has_value())
        {
          View* result = self.setExternalDataPtr(nullptr);
          releaseExternalDataOwner(&self);
          return result;
        }
        return setExternalDataAndPinOwner(self, *external_ptr);
      },
      nb::rv_policy::reference,
      "Set the View to hold undescribed external data, or clear it when passed None.",
      nb::arg("external_ptr").none())
    .def(
      "setExternalData",
      [](View& self, TypeID type, IndexType num_elems, const nb::ndarray<>& external_ptr) {
        return setExternalDataAndPinOwner(self, type, num_elems, external_ptr);
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
        return setExternalDataAndPinOwner(self, type, ndims, shape, external_ptr);
      },
      nb::rv_policy::reference,
      "Set the View to hold described external data (numpy array).")

    .def("getString",
         &View::getString,
         nb::rv_policy::reference,
         "Return the string contained in the View.")
    .def("getDataArray",
         &viewToNumpyArray,
         "Return the data held by the View as a numpy array.\n\n"
         "The array is a zero-copy view into the View's storage "
         "and keeps the View (and its DataStore) alive while referenced. "
         "View.reallocate() (or reallocating the underlying Buffer) can move the storage, "
         "leaving a previously returned array pointing at freed memory; "
         "re-acquire the array after any reallocation.")

    .def(
      "getDataInt",
      [](View& self) { return self.getData<int>(); },
      nb::rv_policy::reference,
      "Return the scalar data held by the View as an python int type.")
    .def(
      "getDataFloat",
      [](View& self) { return self.getData<double>(); },
      nb::rv_policy::reference,
      "Return the data held by the View as a python float type (C++ double).")
    .def("print",
         nb::overload_cast<>(&View::print, nb::const_),
         "Print JSON description of the View.")
    .def("rename", &View::rename, "Change the name of the View.")

    // Attribute accessors
    .def("getAttribute",
         nb::overload_cast<IndexType>(&View::getAttribute),
         nb::rv_policy::reference_internal,
         "Get Attribute by index")
    .def("getAttribute",
         nb::overload_cast<const std::string&>(&View::getAttribute),
         nb::rv_policy::reference_internal,
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
    .def(
      "getAttributeNodeRef",
      [](View& self, IndexType idx) {
        conduit::Node& node = const_cast<conduit::Node&>(self.getAttributeNodeRef(idx));
        return nodeToNbObject(node);
      },
      nb::rv_policy::reference,
      "Return reference to Attribute Node (by index)")
    .def(
      "getAttributeNodeRef",
      [](View& self, const std::string& name) {
        conduit::Node& node = const_cast<conduit::Node&>(self.getAttributeNodeRef(name));
        return nodeToNbObject(node);
      },
      nb::rv_policy::reference,
      "Return reference to Attribute Node (by name)")
    .def(
      "getAttributeNodeRef",
      [](View& self, const Attribute* attr) {
        conduit::Node& node = const_cast<conduit::Node&>(self.getAttributeNodeRef(attr));
        return nodeToNbObject(node);
      },
      nb::rv_policy::reference,
      "Return reference to Attribute Node (by pointer)")

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
    .def("checksum",
         nb::overload_cast<bool>(&Group::checksum, nb::const_),
         "Return a checksum for the Group's name, child structure, and descendant view data.",
         nb::arg("includeAttributes") = true)
    .def(
      "checksum",
      [](const Group& self, nb::object& o, bool includeAttributes) {
        conduit::Node& cpp_node = nbObjectToNode(o);
        self.checksum(cpp_node, includeAttributes);
      },
      "Populate a Conduit node with checksum metadata for this Group hierarchy.",
      nb::arg("n_checksum"),
      nb::arg("includeAttributes") = true)
    .def("getParent",
         nb::overload_cast<>(&Group::getParent, nb::const_),
         nb::rv_policy::reference_internal,
         "Return pointer to non-const parent Group of a Group.")
    .def("getNumGroups", &Group::getNumGroups, "Return number of child Groups in a Group object.")
    .def("getNumViews", &Group::getNumViews, "Return number of Views owned by a Group object.")
    .def("getDataStore",
         nb::overload_cast<>(&Group::getDataStore, nb::const_),
         nb::rv_policy::reference_internal,
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
         nb::rv_policy::reference_internal,
         "Return pointer to const View with given name or path.")
    .def("getView",
         nb::overload_cast<IndexType>(&Group::getView, nb::const_),
         nb::rv_policy::reference_internal,
         "Return pointer to non-const View with given index.")
    .def("getFirstValidViewIndex",
         &Group::getFirstValidViewIndex,
         "Return first valid View index in Group object.")
    .def("getNextValidViewIndex",
         &Group::getNextValidViewIndex,
         "Return next valid View index in Group object after given index.")

    .def("createView",
         nb::overload_cast<const std::string&>(&Group::createView),
         nb::rv_policy::reference_internal,
         "Create an undescribed (i.e., empty) View object with given name or path in this Group.")
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType>(&Group::createView),
         nb::rv_policy::reference_internal,
         "Create View object with given name or path in this Group that has a data description "
         "with data type and number of elements.")
    .def(
      "createViewWithShape",
      [](Group& self, const std::string& path, TypeID type, int ndims, const nb::ndarray<IndexType>& shape) {
        return self.createViewWithShape(path, type, ndims, shape.data());
      },
      nb::rv_policy::reference_internal,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape.")
    .def("createView",
         nb::overload_cast<const std::string&, Buffer*>(&Group::createView),
         nb::rv_policy::reference_internal,
         "Create an undescribed View object with given name or path in this Group and attach given "
         "Buffer to it.")
    .def("createView",
         nb::overload_cast<const std::string&, TypeID, IndexType, Buffer*>(&Group::createView),
         nb::rv_policy::reference_internal,
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
      nb::rv_policy::reference_internal,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and attach given Buffer to it.")

    .def(
      "createView",
      [](Group& self, const std::string& path, const nb::ndarray<>& a) {
        View* view = self.createView(path, a.data());
        pinExternalDataOwner(view, a);
        return view;
      },
      nb::rv_policy::reference_internal)

    .def(
      "createView",
      [](Group& self, const std::string& path, TypeID id, IndexType num_elems, const nb::ndarray<>& a) {
        View* view = self.createView(path, id, num_elems, a.data());
        pinExternalDataOwner(view, a);
        return view;
      },
      nb::rv_policy::reference_internal,
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
        View* view = self.createViewWithShape(path, type, ndims, shape.data(), external_ptr.data());
        pinExternalDataOwner(view, external_ptr);
        return view;
      },
      nb::rv_policy::reference_internal,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and attach externally-owned data (numpy array) to it.")
    .def("createViewAndAllocate",
         nb::overload_cast<const std::string&, TypeID, IndexType, int>(&Group::createViewAndAllocate),
         nb::rv_policy::reference_internal,
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
      nb::rv_policy::reference_internal,
      "Create View object with given name or path in this Group that has a data description "
      "with data type and shape and allocate data for it.")

    .def("createViewScalar",
         &Group::createViewScalar<int>,
         nb::rv_policy::reference_internal,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (int).",
         nb::arg("path"),
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("createViewScalar",
         &Group::createViewScalar<double>,
         nb::rv_policy::reference_internal,
         "Create View object with given name or path in this Group set its data to given scalar "
         "value (C++ double, python float).",
         nb::arg("path"),
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)
    .def("createViewString",
         &Group::createViewString,
         nb::rv_policy::reference_internal,
         "Create View object with given name or path in this Group set its data to given string.",
         nb::arg("path"),
         nb::arg("value").noconvert(),
         nb::arg("allocID") = INVALID_ALLOCATOR_ID)

    .def(
      "destroyView",
      [](Group& self, const std::string& path) {
        // releaseExternalDataOwner is null-safe; destroyView handles invalid paths gracefully
        releaseExternalDataOwner(self.getView(path));
        self.destroyView(path);
      },
      "Destroy View with given name or path owned by this Group, but leave its data intact.")
    .def(
      "destroyView",
      [](Group& self, IndexType idx) {
        // releaseExternalDataOwner is null-safe; destroyView handles invalid indices gracefully
        releaseExternalDataOwner(self.getView(idx));
        self.destroyView(idx);
      },
      "Destroy View with given index owned by this Group, but leave its data intact.")
    .def(
      "destroyViewAndData",
      [](Group& self, const std::string& path) {
        // releaseExternalDataOwner is null-safe; destroyViewAndData handles invalid paths gracefully
        releaseExternalDataOwner(self.getView(path));
        self.destroyViewAndData(path);
      },
      "Destroy View with given name or path owned by this Group and deallocate")
    .def(
      "destroyViewAndData",
      [](Group& self, IndexType idx) {
        releaseExternalDataOwner(self.getView(idx));
        self.destroyViewAndData(idx);
      },
      "Destroy View with given index owned by this Group and deallocate its data if it's the "
      "only View associated with that data.")
    .def(
      "destroyViewsAndData",
      [](Group& self) {
        releaseExternalDataOwnersOfViews(self);
        self.destroyViewsAndData();
      },
      "Destroy all Views owned by this Group and deallocate "
      "data for each View when it's the only View associated with that data.")

    .def("moveView",
         &Group::moveView,
         nb::rv_policy::reference_internal,
         "Remove given View object from its owning Group and move it to this Group.")
    .def(
      "copyView",
      [](Group& self, View* view) {
        View* copy = self.copyView(view);
        if(copy != nullptr && view != nullptr && view->isExternal())
        {
          // Shallow copy shares external data pointer - copy the pin to prevent
          // the numpy array from being garbage collected
          copyExternalDataOwner(view, copy);
        }
        return copy;
      },
      nb::rv_policy::reference_internal,
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
         nb::rv_policy::reference_internal,
         "Return pointer to non-const child Group with given name or path.")
    .def("getGroup",
         nb::overload_cast<IndexType>(&Group::getGroup),
         nb::rv_policy::reference_internal,
         "Return pointer to non-const immediate child Group with given index.")
    .def("views",
         nb::overload_cast<>(&Group::views),
         nb::keep_alive<0, 1>(),
         "Return an iterator over Views")
    .def("groups",
         nb::overload_cast<>(&Group::groups),
         nb::keep_alive<0, 1>(),
         "Return an iterator over Groups")
    .def("getFirstValidGroupIndex",
         &Group::getFirstValidGroupIndex,
         "Return first valid child Group index (i.e., smallest index over all child Groups).")
    .def("getNextValidGroupIndex",
         &Group::getNextValidGroupIndex,
         "Return next valid child Group index after given index.")
    .def("createGroup",
         &Group::createGroup,
         nb::rv_policy::reference_internal,
         "Create a child Group within this Group with given name or path.",
         nb::arg("path"),
         nb::arg("is_list") = false,
         nb::arg("accept_existing") = false)
    .def("createUnnamedGroup",
         &Group::createUnnamedGroup,
         nb::rv_policy::reference_internal,
         "Create a child Group within this Group with no name.",
         nb::arg("is_list") = false)
    .def(
      "destroyGroup",
      [](Group& self, const std::string& path) {
        // releaseExternalDataOwners is null-safe; destroyGroup handles invalid paths gracefully
        releaseExternalDataOwners(self.getGroup(path));
        self.destroyGroup(path);
      },
      "Destroy child Group in this Group with given name or path.")
    .def(
      "destroyGroup",
      [](Group& self, IndexType idx) {
        // releaseExternalDataOwners is null-safe; destroyGroup handles invalid indices gracefully
        releaseExternalDataOwners(self.getGroup(idx));
        self.destroyGroup(idx);
      },
      "Destroy child Group within this Group with given index.")
    .def(
      "destroyGroupAndData",
      [](Group& self, const std::string& path) {
        // releaseExternalDataOwners is null-safe; destroyGroupAndData handles invalid paths gracefully
        releaseExternalDataOwners(self.getGroup(path));
        self.destroyGroupAndData(path);
      },
      "Destroy child Group at the given path, and destroy data that is "
      "not shared elsewhere.")
    .def(
      "destroyGroupAndData",
      [](Group& self, IndexType idx) {
        releaseExternalDataOwners(self.getGroup(idx));
        self.destroyGroupAndData(idx);
      },
      "Destroy child Group with the given index, and destroy data that "
      "is not shared elsewhere.")
    .def(
      "destroyGroupsAndData",
      [](Group& self) {
        for(auto& g : self.groups())
        {
          releaseExternalDataOwners(&g);
        }
        self.destroyGroupsAndData();
      },
      "Destroy all child Groups held by this Group, and destroy data that "
      "is not shared elsewhere.")
    .def(
      "destroyGroupSubtreeAndData",
      [](Group& self) {
        releaseExternalDataOwners(&self);
        self.destroyGroupSubtreeAndData();
      },
      "Destroy the entire subtree of Groups and Views held by this Group, "
      "and destroy data that is not shared elsewhere.")
    .def(
      "destroyGroups",
      [](Group& self) {
        for(auto& g : self.groups())
        {
          releaseExternalDataOwners(&g);
        }
        self.destroyGroups();
      },
      "Destroy all child Groups in this Group.")
    .def("moveGroup",
         &Group::moveGroup,
         nb::rv_policy::reference_internal,
         "Remove given Group object from its parent Group and make it a child of this Group.")
    .def(
      "copyGroup",
      [](Group& self, Group* group) {
        Group* copy = self.copyGroup(group);
        if(copy != nullptr && group != nullptr)
        {
          // Shallow copy shares external data pointers - recursively copy all pins
          copyExternalDataOwners(group, copy);
        }
        return copy;
      },
      nb::rv_policy::reference_internal,
      "Create a (shallow) copy of Group hierarchy rooted at given "
      "Group and make the copy a child of this Group.")
    .def("deepCopyGroup",
         &Group::deepCopyGroup,
         nb::rv_policy::reference_internal,
         "Create a deep copy of Group hierarchy rooted at given Group and "
         "make the copy a child of this Group.",
         nb::arg("srcGroup"),
         nb::arg("arrayAllocID") = INVALID_ALLOCATOR_ID,
         nb::arg("tupleAllocID") = INVALID_ALLOCATOR_ID)

    .def("print",
         nb::overload_cast<>(&Group::print, nb::const_),
         "Print JSON description of data Group to stdout.")

    .def(
      "createExternalLayout",
      [](Group& self, nb::object& o, const Attribute* attr) {
        conduit::Node& cpp_node = nbObjectToNode(o);
        return self.createExternalLayout(cpp_node, attr);
      },
      "Copy data Group external layout to given Conduit node.",
      nb::arg("n"),
      nb::arg("attr") = nb::none())

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
         nb::arg("protocol") = Group::getDefaultIOProtocol(),
         nb::arg("preserve_contents") = false)

    .def("loadExternalData",
         nb::overload_cast<const std::string&>(&Group::loadExternalData),
         "Load data into the Group's external views from a file.")
    .def_static("getDefaultIOProtocol",
                &Group::getDefaultIOProtocol,
                "Return the default I/O protocol for this Axom build.")
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

    .def(
      "getDefaultNodeRef",
      [](Attribute& self) {
        conduit::Node& node = const_cast<conduit::Node&>(self.getDefaultNodeRef());
        return nodeToNbObject(node);
      },
      nb::rv_policy::reference,
      "Return default value of Attribute as Node reference.")
    .def("getTypeID", &Attribute::getTypeID, "Return type of Attribute.");

#if defined(AXOM_USE_MPI)
  nb::class_<PyIOManager>(m_sidre, "IOManager")
    .def(
      "__init__",
      [](PyIOManager* self, bool use_scr) { new(self) PyIOManager(MPI_COMM_WORLD, use_scr); },
      "Create an IOManager on MPI_COMM_WORLD.",
      nb::arg("use_scr") = false)
    .def(
      "__init__",
      [](PyIOManager* self, nb::object comm, bool use_scr) {
        // Accept an mpi4py communicator without depending on mpi4py's C headers.
        // mpi4py Comm objects expose py2f(), which returns the Fortran integer handle.
        // MPI_Comm_f2c converts that handle back to an MPI_Comm.
        // PyIOManager duplicates the communicator, so callers may drop or free
        // the mpi4py communicator after construction.
        // Passing None uses MPI_COMM_WORLD, while no-argument calls are
        // handled by the bool/use_scr overload above.
        new(self) PyIOManager(mpiCommFromObject(comm), use_scr);
      },
      "Create an IOManager on an mpi4py communicator, or MPI_COMM_WORLD when comm is None.",
      nb::arg("comm"),
      nb::arg("use_scr") = false)
    .def(
      "write",
      [](PyIOManager& self,
         Group* group,
         int num_files,
         const std::string& file_base,
         const std::string& protocol,
         const std::string& tree_pattern) {
        self.manager().write(group, num_files, file_base, protocol, tree_pattern);
      },
      "Write a Group to output files.",
      nb::arg("group"),
      nb::arg("num_files"),
      nb::arg("file_base"),
      nb::arg("protocol"),
      nb::arg("tree_pattern") = "datagroup")
    .def(
      "read",
      [](PyIOManager& self,
         Group* group,
         const std::string& root_file,
         const std::string& protocol,
         bool preserve_contents) {
        self.manager().read(group, root_file, protocol, preserve_contents);
      },
      "Read from input file.",
      nb::arg("group"),
      nb::arg("root_file"),
      nb::arg("protocol"),
      nb::arg("preserve_contents") = false)
    .def(
      "read",
      [](PyIOManager& self, Group* group, const std::string& root_file, bool preserve_contents) {
        self.manager().read(group, root_file, preserve_contents);
      },
      "Read from a root file.",
      nb::arg("group"),
      nb::arg("root_file"),
      nb::arg("preserve_contents") = false)
    .def(
      "loadExternalData",
      [](PyIOManager& self, Group* group, const std::string& root_file) {
        self.manager().loadExternalData(group, root_file);
      },
      "Load external data into a group.",
      nb::arg("group"),
      nb::arg("root_file"))
    .def(
      "loadExternalData",
      [](PyIOManager& self, Group* parent_group, Group* load_group, const std::string& root_file) {
        self.manager().loadExternalData(parent_group, load_group, root_file);
      },
      "Piecewise load of external data into a group.",
      nb::arg("parent_group"),
      nb::arg("load_group"),
      nb::arg("root_file"))
    .def(
      "getNumFilesFromRoot",
      [](PyIOManager& self, const std::string& root_file) {
        return self.manager().getNumFilesFromRoot(root_file);
      },
      "Gets the number of files in the dataset from the specified root file.",
      nb::arg("root_file"))
    .def(
      "getNumGroupsFromRoot",
      [](PyIOManager& self, const std::string& root_file) {
        return self.manager().getNumGroupsFromRoot(root_file);
      },
      "Gets the number of groups in the dataset from the specified root file.",
      nb::arg("root_file"))
    .def_static("correspondingRelayProtocol",
                &IOManager::correspondingRelayProtocol,
                "Finds conduit relay protocol corresponding to a sidre protocol.");
#endif
}

} /* end namespace sidre */
} /* end namespace axom */
