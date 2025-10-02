.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******
Views
******

The BUMP component provides lightweight, device-compatible, view classes that add a C++ interface
for Blueprint data. Blueprint data defines several object protocols represented with arrays of
various data types. Views can simplify the process of writing algorithms to support Blueprint data.
Views do not own their data so they can be easily copied, making them suitable for use in device
kernels.

----------
ArrayView
----------

Axom provides ``axom::ArrayView`` to wrap data in a non-owning data structure that can be passed to
kernels. The BUMP component provides the ``axom::bump::utilities::make_array_view()``
function to help wrap arrays stored in ``conduit::Node`` to ``axom::ArrayView``. To use the
``make_array_view`` function, one must know the type held within the Conduit node. If that is
not the case, then consider using one of the dispatch ''Node_to_ArrayView'' functions.

.. code-block:: cpp

    // Make an axom::ArrayView<float> for X coordinate components.
    auto x = axom::bump::utilities::make_array_view<float>(n_mesh["coordsets/coords/values/x"]);


----------
Coordsets
----------

Blueprint supports multiple coordset *(coordinate set)* types: uniform, rectilinear, explicit. Axom provides functions
to explicitly create coordset views for each of these types.

.. code-block:: cpp

    namespace views = axom::bump::views;
    // Make a 2D uniform coordset view
    auto view1 = views::make_uniform_coordset<2>::view(n_mesh["coordsets/coords"]);
    // Make a 3D uniform coordset view
    auto view2 = views::make_uniform_coordset<3>::view(n_mesh["coordsets/coords"]);
    // Make a 2D rectilinear coordset view with float coordinates
    auto view3 = views::make_rectilinear_coordset<float, 2>::view(n_mesh["coordsets/coords"]);
    // Make a 3D rectilinear coordset view with double coordinates
    auto view4 = views::make_rectilinear_coordset<double, 3>::view(n_mesh["coordsets/coords"]);
    // Make a 2D explicit coordset view with float coordinates
    auto view5 = views::make_explicit_coordset<float, 2>::view(n_mesh["coordsets/coords"]);
    // Make a 3D explicit coordset view with double coordinates
    auto view6 = views::make_explicit_coordset<double, 3>::view(n_mesh["coordsets/coords"]);


----------------
Topology Views
----------------

Topology views provide a layer on top of the Blueprint mesh topology that enables Axom algorithms
to be written while not needing to care specifically about topology types and data types.
Axom provides topology views for structured meshes and unstructured meshes. 

^^^^^^^^^^^^^^^^^^^^^^
Structured Mesh Views
^^^^^^^^^^^^^^^^^^^^^^

The structured mesh topology view, ``axom::bump::views::StructuredTopologyView``, pertains to any of the Blueprint
structured topology types. The ``StructuredTopologyView`` class is a template that takes an indexing policy
as a template argument. The indexing policy computes zone indices and converts to/from
logical/global indices. The ``StridedStructuredIndexingPolicy`` class supports indexing for
strided-structured Blueprint meshes, which are structured meshes that exist over a sub-window
of the overall mesh. There are helper functions for creating structured topology views from
a Conduit node.

.. code-block:: cpp

    namespace views = axom::bump::views;
    conduit::Node &n_topo1 = n_mesh["topologies/mesh2d"];
    conduit::Node &n_topo2 = n_mesh["topologies/mesh3d"];
    conduit::Node &n_topo3 = n_mesh["topologies/mesh2dss"];
    // Make a 2D structured mesh view from the topology.
    auto topologyView1 = views::make_structured_topology<2>::view(n_topo1);
    // Make a 3D structured mesh view from the topology.
    auto topologyView2 = views::make_structured_topology<2>::view(n_topo2);
    // Make a 2D strided-structured mesh view from the topology.
    auto topologyView3 = views::make_strided_structured_topology<2>::view(n_topo3);

^^^^^^^^^^^^^^^^^^^^^^^^
Unstructured Mesh Views
^^^^^^^^^^^^^^^^^^^^^^^^

There are 3 unstructured mesh views, covering single shape meshes, mixed shape meshes, and polyhedral meshes.
The ``axom::bump::views::UnstructuredTopologySingleShapeView`` class wraps a
Blueprint topology that contains a single zone/shape type. The zone type is a template argument
that determines the type of zone that is held within the topology.

.. code-block:: cpp

    // Make a topology view for a tetrahedral mesh with int connectivity.
    namespace views = axom::bump::views;
    namespace utils = axom::bump::utilities;
    const conduit::Node &n_topo = n_mesh["topologies/mesh"];
    const auto connView = utils::make_array_view<int>(n_topo["elements/connectivity"]);
    views::UnstructuredTopologySingleShapeView<views::TetShape<int>> view(connView);

There are multiple shape types defined in ``axom/bump/views/Shapes.hpp`` that can be used with
the ``UnstructuredTopologySingleShapeView`` class: *TriShape*, *QuadShape*, *TetShape*, *PyramidShape*,
*WedgeShape*, and *HexShape*.

Blueprint supports *mixed* topologies that contain multiple shape types. These topologies are
handled using the ``axom::bump::views::UnstructuredTopologyMixedShapeView``. Additional array
views are needed to supply the sizes, offsets, and shapes arrays.

.. code-block:: cpp

    // A shape map helps map values from the values used in the Blueprint topology to
    // the shape ids used in Axom.
    namespace views = axom::bump::views;
    namespace utils = axom::bump::utilities;

    const conduit::Node &n_topo = n_mesh["topologies/mesh"];
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<int> ids, values;
    auto shapeMap = views::buildShapeMap(n_topo, ids, values, allocatorID);

    views::UnstructuredTopologyMixedShapeView<int> view(
      utils::make_array_view(n_topo["elements/connectivity"),
      utils::make_array_view(n_topo["elements/sizes"),
      utils::make_array_view(n_topo["elements/offsets"),
      utils::make_array_view(n_topo["elements/shapes"),
      shapeMap);

The final unstructured topology view is ``axom::bump::views::UnstructuredTopologyPolyhedralView``
and it provides a view interface to polyhedral meshes.

.. literalinclude:: ../../views/dispatch_unstructured_topology.hpp
   :start-after: _bump_views_ph_topoview_begin
   :end-before: _bump_views_ph_topoview_end
   :language: C++


Once a suitable topology view type has wrapped a Blueprint topology, it can be used in
device kernels to obtain zone information.

.. code-block:: cpp

    auto topologyView = ...
    axom::for_all<ExecSpace>(topologyView.numberOfZones(), AXOM_LAMBDA(axom::IndexType zoneIndex)
    {
      // Get the current zone.
      const auto zone = topologyView.zone(zoneIndex);

      // Iterate over this zone's nodes.
      for(const auto &nodeId : zone.getIds())
      {
        // Do something.
      }
    });

----------
Matsets
----------

The BUMP component provides material views to wrap Blueprint matsets behind an interface that
supports queries of the matset data without having to care much about its internal representation.
Blueprint provides 4 flavors of matset, each with a different representation. The 
``axom::bump::views::UnibufferMaterialView`` class wraps unibuffer matsets, which consist of
several arrays that define materials for each zone in the associated topology. The view's
methods allow algorithms to query the list of materials for each zone.

.. literalinclude:: ../../tests/bump_views.cpp
   :start-after: _bump_views_matsetview_begin
   :end-before: _bump_views_matsetview_end
   :language: C++

----------
Dispatch
----------

There are several helper functions that wrap a Conduit node in a specific view type
and **dispatch** the view to a user-supplied lambda for further processing. The lambda
function will typically be instantiated multiple times to handle various data types and
object types (e.g. coordsets). Generic lambdas can be used to process multiple views
and it is possible to nest dispatch functions to handle added complexity.

^^^^^^^^^^^
Array Data
^^^^^^^^^^^

Blueprint data can readily be wrapped in ``axom::ArrayView`` using the ``axom::bump::utilities::make_array_view()``
function. There are dispatch functions for ``conduit::Node`` data arrays that automate the
wrapping to ``axom::ArrayView`` and passing the views to a user-supplied lambda.

To generically wrap any type of datatype supported by Conduit, the ``axom::bump::views::Node_to_ArrayView()``
function can be used. This template function takes a variable number of ``conduit::Node``
arguments and a generic lambda function that accepts the view arguments. The lambda gets
instantiated for every supported Conduit data type.

.. code-block:: cpp

    conduit::Node n; // Assume it contains data values
    axom::bump::views::Node_to_ArrayView(n["foo"], n["bar"], [&](auto fooView, auto barView)
    {
      // Use fooView and barView axom::ArrayView objects to access data.
      // They can have different types.
    });

Using ``axom::bump::views::Node_to_ArrayView`` with multiple data values can instantiate
the supplied lambda many times so be careful. It is more common when wrapping multiple
nodes that they are the same type. The ``axom::bump::views::Node_to_ArrayView_same`` function
ensures that the lambdas get instantiated with views that wrap the Conduit  nodes in
array views that of the same type.

.. code-block:: cpp

    conduit::Node n; // Assume it contains data values
    axom::bump::views::Node_to_ArrayView_same(n["foo"], n["bar"], [&](auto fooView, auto barView)
    {
      // Use fooView and barView axom::ArrayView objects to access data.
      // They have the same types.
    });

When dealing with mesh data structures, it is common to have data that are using only integer
types or only floating-point types. Axom provides functions that limit the lambda instantiation
to only those selected types using the following functions:

 * ``axom::bump::views::IndexNode_to_ArrayView()``
 * ``axom::bump::views::IndexNode_to_ArrayView_same()``
 * ``axom::bump::views::FloatNode_to_ArrayView()``
 * ``axom::bump::views::FloatNode_to_ArrayView_same()``

The "Index" functions limit lambda instantiation to common index types signed/unsigned 32/64-bit
integers. The "Float" functions instantiate lambdas with float32 and float64 types.


^^^^^^^^^^^
Coordsets
^^^^^^^^^^^

The ``axom::bump::views::dispatch_coordset()`` function can wrap Blueprint coordsets in an
appropriate view and pass it to a lambda function.

.. code-block:: cpp

   const conduit::Node &n_coordset = n_mesh["coordsets/coords"];
   axom::bump::views::dispatch_coordset(n_coordset, [&](auto coordsetView) {
     // Get the C++ type of the coordset.
     using CoordsetView = decltype(coordsetView);
     // Implement algorithm using coordsetView.
   });

^^^^^^^^^^^
Topologies
^^^^^^^^^^^

Dispatch functions for topologies enable creation of algorithms that can operate on multiple
topology types through a topology view. These dispatch functions can be called for specific
topology types such as unstructured topologies or they can be called to implement algorithms
that can operate on any topology.

.. code-block:: cpp

    namespace views = axom::bump::views;
    const conduit::Node &n_topo = n_mesh["topologies/mesh"];
    // Handle rectilinear topology type.
    views::dispatch_rectilinear_topology(n_topo, [&](auto topologyView) {
    });
    // Handle structured topology types
    views::dispatch_structured_topology(n_topo, [&](auto topologyView) {
    });
    // Handle unstructured topology types
    views::dispatch_unstructured_topology(n_topo, [&](auto topologyView) {
    });
    // Handle any topology type.
    views::dispatch_topology(n_topo, [&](auto topologyView) {
    });

Nesting dispatch functions permits the calling code to handle both coordset views and
topology views using a compact amount of code. For portability, the actual
algorithm should be placed in a function or class member method when instantiated from
the anonymous lambda function from the dispatch functions.

.. code-block:: cpp

    struct Algorithm
    {
      void execute(const conduit::Node &n_mesh)
      {
        namespace views = axom::bump::views;
        // Handle product of coordset types and topology types.
        views::dispatch_coordset(n_mesh["coordsets/coords"], [&](auto coordsetView)
        {
          views::dispatch_topologies(n_mesh["topologies/mesh"], [&](auto topologyView)
          {
            implementation(coordsetView, topologyView);
          });
        });
      }

      template <typename CoordsetView, typename TopologyView>
      void implementation(CoordsetView coordsetView, TopologyView topologyView) const
      {
        // Algorithm that involves coordsetView and topologyView.
      }
    };

