.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _group-label:

==========
Group
==========

Sidre ``Group`` objects are used to define a tree-like hierarchical organization
for application data, such as meshes and fields used in a simulation. Each
group has a name and one parent group (except for the root group, which has no
parent) and contains zero or more child groups and zero or more data views.
A Sidre datastore has exactly one root group that is created when the
datastore is created. The root group's name is initially the empty string.
See :ref:`datastore-label` for more information.

A group hierarchy is constructed by creating child groups of the root group,
children of those groups, and so on. All groups in a subtree rooted at a
particular group are considered descendants of that group. ``View`` objects
can be created in groups to hold or provide access to data of various types.

.. note:: ``Group`` and ``View`` objects can only be created and destroyed using
          ``Group`` class methods. The ``Group`` and ``View`` class
          constructors and destructors are private.

A group or view is owned by the group that created it; i.e., its parent group
or *owning* group, respectively. Groups and views maintain pointers to
their parent/owning group. Thus, one may *walk* up or down a group hierarchy
to access groups and views in it, as needed.

.. note:: * The name (string) of a group or view **must be unique** within its
            parent/owning Group.
          * A group or view has a unique integer identifier within its
            parent/owning group, which is generated when it is created.
          * Views and child groups in a group can be accessed by name or
            integer id.

A group can be moved or copied to another group. When a group is moved
to another group, it is removed from its original parent and the group to
which it is moved becomes its parent. This implies that the entire subtree
of groups and views within the moved group is moved as well and can no longer
be accessed via the original parent group. When a group is copied to another
group, a copy of the entire group subtree rooted at the copied group is added
to the group to which it is copied. A **shallow** copy is performed for the
data in each view. Specifically, a new view objects is created in the
destination, but the data is shared by the original and new view.

.. note:: ``View`` object copy operations perform **shallow** copies of the
          data in a view.

Some methods for creating, destroying, querying, and retrieving groups and
views take a string with *path syntax*, where parent and child group names
are joined with the path separator character '/'. Other methods take the name
of an immediate child of a group. Methods that require the name of an immediate
child are marked with 'Child' in their name, such as ``hasChildView()`` and
``hasChildGroup()``. When a path string is passed to a method that accepts
path syntax, the last item in the path indicates the item to be created,
destroyed, accessed, etc.  For example,::

   View* view = group->createView("foo/bar/baz");

is equivalent to::

   View* view = group->createGroup("foo")->createGroup("bar")->createView("baz");

In particular, intermediate groups "foo" and "bar" will be created in this
case if they don't already exist. The path syntax is similar to a Unix
filesystem, but the path string **may not** contain the parent entry,
such as "../foo", or current group, such as "./bar".

----------------------------
Methods to Operate on Groups
----------------------------

The following lists summarize ``Group`` class methods that support operations
related to groups.

.. note:: * Methods that access groups by index only work with the immediate
            children of the current group because an id has no meaning
            outside of the indexing of the current group. None of these methods
            is marked with 'Child' in its name.
          * When a group is created, destroyed, copied, or moved,
            ids of other views and groups in its parent group may
            become invalid. This is analogous to iterator invalidation for
            containers when the container contents change.

Create, Modify, and Destroy Groups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 * Create a child group given a name (child) or path (other descendant).
   If a path is given, intermediate groups in path are created, if needed.
 * Rename a group.  A group cannot be renamed to the empty string, to
   a string containing the path separator character, or to the name of
   another group or view owned by the same parent.
 * Destroy a descendant group with given id (child), or name/path (child or
   other descendant).
 * Destroy all child groups in a group.

.. note:: When a ``Group`` object is destroyed, all groups and views in the
          subtree rooted at the destroyed group are also destroyed. However,
          the data associated with the views will remain intact.

Group Properties
^^^^^^^^^^^^^^^^^^^^^^^

 * Retrieve the name or id of a group
 * Retrieve the full path name from the root of the tree to a group
 * Get a pointer to the parent group of a group
 * Query the number of child groups of a group
 * Query whether a group has a descendant group with a given name or path
 * Query whether a group has an immediate child group with a given integer id
 * Query the name of a child group with a given id, or the id of a child group
   with a given name
 * Get a pointer to the datastore that owns the hierarchy in which a group
   resides

Group Access
^^^^^^^^^^^^

 * Retrieve an immediate child group with a given name or id, or a descendant
   group with a given path
 * Iterate over the set of child groups of a group.
   One can use the "range-for" syntax or the iterator syntax

  .. code-block:: C++

     // 'range-for' syntax:

     for(auto& grp: someGroup->groups()) { /* ... */ }

     // 'iterator' syntax:
     for(auto it = someGroup->groups().begin(),
           itEnd = someGroup->groups().end(); it != itEnd; ++it)
      {
        auto& grp = *it;
        /* ... */
      }

Move and Copy Groups
^^^^^^^^^^^^^^^^^^^^^^

 * Move a group, and its associated subtree, from its parent group and make it
   a child of another group
 * Create a copy of group subtree rooted at some group and make it a child of
   another group
 * Query whether a group subtree is equivalent to another; i.e., identical
   subtree structures with same names for all groups and views, and views are
   also equivalent (see :ref:`view-interface-label`).

----------------------------
Methods to Operate on Views
----------------------------

``Group`` class methods that support operations related to ``View`` objects are
summarized below. For more details on View concepts and operations, please
see :ref:`view-label`.

.. note:: Methods that access views by index work only with the
          views owned by the current group because an id has no meaning
          outside of the indexing of the current group. None of these methods
          is marked with 'Child' in its name.

Create Views
^^^^^^^^^^^^^

 * Create a view in the group with a name only.
 * Create a view in the group with a name and data description.
 * Create a view in the group with a name and with a Buffer attached. The
   View may or may not have a data description.
 * Create a view in the group with a name and an external data pointer. The
   data may or may not be described.
 * Create a view in the group with a name and data description, and allocate
   the data. Implicitly the data is held in a buffer that is attached to the
   view.
 * Create a view in the group with a name holding a given scalar or string.

Destroy Views
^^^^^^^^^^^^^^

 * Destroy view with given id (child), or name/path (view in the group or some
   descendant group), and leave view data intact.
 * Destroy all views in the group, and leave their data intact.
 * Destroy view with given id, or name/path, and destroy their data.
 * Destroy all views in the group and destroy their data.

View Queries
^^^^^^^^^^^^^^^^

 * Query the number of views in a group.
 * Query whether a group subtree has a view with a given name or path.
 * Query whether a group has a view with a given integer id.
 * Query the name of a view with a given id, or the id of a view with a given
   name.

View Access
^^^^^^^^^^^^^

 * Retrieve a view in the group with a given name or id, or a descendant view
   (somewhere in the subtree) with a given path.
 * Iterate over the set of views owned by the group.
   One can use the "range-for" syntax or the iterator syntax

  .. code-block:: C++

     // 'range-for' syntax:
     for(auto& view: someGroup->views()) { /* ... */ }

     // 'iterator' syntax:
     for(auto it = someGroup->views().begin(),
           itEnd = someGroup->views().end(); it != itEnd; ++it)
      {
        auto& view = *it;
        /* ... */
      }


Move and Copy Views
^^^^^^^^^^^^^^^^^^^^

 * Move a view from its owning group to another group (removed from original
   owning group).
 * Copy a view to another group. Note that this is a **shallow** copy of the
   view data; i.e., it is shared by the original and the new view in the
   destination group.

----------------------------
List Format
----------------------------

The list format is an alternate way for a group to hold its child groups and
views. In list format, any number of child group or view items can be created.
Each can be accessed, in order of creation, using an iterator over groups or
over views. Child groups and views held in list format cannot be accessed by
name.

To create a group that uses the list format, the optional argument ``is_list``
must be set to ``true`` in the call to ``createGroup``.

  .. code-block:: C++

     // list_group will hold its child items in the list format.
     Group* list_group = group->createGroup("my_list", true);

It is recommended but not required that the items held in the list format
be created without names. String names may be assigned to these items,
but the names will not be useful for accessing them from their parent
group, and none of the methods that access child items by name or path will
return a valid pointer. The method ``createUnnamedGroup`` is available to
create an unnamed child group, while unnammed views can be created by passing
an empty string to any of the several ``createView`` methods in the ``Group``
class.

  .. code-block:: C++
          
     Group* g0 = list_group->createUnnamedGroup();
     Group* g1 = list_group->createUnnamedGroup();
     Group* g2 = list_group->createUnnamedGroup();
     View* v0 = list_group->createView("");
     View* v1 = list_group->createViewScalar("", 1.0);
     View* v2 = list_group->createViewString("", "foo");
     View* v3 = list_group->createView("", type, num_elems, buffer);

While it is allowed to pass a non-empty string to be the name of a child
item held in the list format, a string with path syntax, like
``"foo/bar/baz"``, will be considered invalid, and the object creation methods
will return a nullptr if such a string is provided.

  .. code-block:: C++

     // This is valid, but the string name will not be useful for future access.
     View* foo = list_group->createView("foo");
     // This is invalid due to the path syntax, a nullptr will be returned.
     View* baz = list_group->createView("bar/baz");

----------------------------
Group I/O Operations
----------------------------

The group interface provides methods to perform data I/O operations on views
in the group subtree rooted at any group.

 * Copy a description of a group subtree to a ``conduit::Node``.
 * Create native and external data layouts in ``conduit::Node`` hierarchies
   (used mainly for I/O operations)
 * Save and load group subtrees, including data in associated views, to and
   from files. A variety of methods are provided to support different I/O
   operations, different I/O protocols, etc.

I/O methods on the group class use `Conduit <https://github.com/LLNL/conduit>`_
to :ref:`write the data (sub)tree <sidre-serial-io>` rooted in a group to a
file, `HDF5 <https://www.hdfgroup.org/HDF5/>`_ handle, or other
Conduit protocol, or to an in-memory Conduit data structure. Please see
:ref:`sidre-conduit` for more information. An application may
provide an attribute to the method call, so only views with that attribute
explicitly set will be written. See :ref:`spio-core-concepts` for more
information.

