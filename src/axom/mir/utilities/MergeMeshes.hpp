// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_MERGE_MESHES_HPP_
#define AXOM_MIR_MERGE_MESHES_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/mir/views/Shapes.hpp"
#include "axom/mir/views/dispatch_material.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <string>

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/**
 * \brief A mesh input containing a Blueprint mesh and some mapping array views.
 */
struct MeshInput
{
  conduit::Node *m_input {nullptr};                   //!< Pointer to Blueprint mesh.
  axom::ArrayView<axom::IndexType> m_nodeMapView {};  //!< Map for mesh nodeIds to nodeIds in final mesh.
  axom::ArrayView<axom::IndexType> m_nodeSliceView {};  //!< Node ids to be extracted and added to final mesh.
};

/**
 * \brief Merge multiple unstructured Blueprint meshes through MeshInput.
 *
 * \note The input meshes must currently contain a single coordset/topology/matset.
 */
template <typename ExecSpace>
class MergeMeshes
{
public:
  /**
   * \brief Merge the input Blueprint meshes into a single Blueprint mesh.
   *
   * \param inputs A vector of inputs to be merged.
   * \param options A Node containing algorithm options.
   * \param[out] output The node that will contain the merged mesh.
   *
   * The options node may contain a "topology" string that designates the name of
   * the topology to be merged.
   */
  void execute(const std::vector<MeshInput> &inputs,
               const conduit::Node &options,
               conduit::Node &output) const
  {
    AXOM_ANNOTATE_SCOPE("MergeMeshes");
    bool ok = validInputs(inputs, options);
    if(!ok)
    {
      SLIC_ASSERT_MSG(ok, "Unsupported inputs were provided.");
      return;
    }

    if(inputs.size() == 1)
    {
      singleInput(inputs, output);
    }
    else if(inputs.size() > 1)
    {
      mergeInputs(inputs, output);
    }
  }

// The following members are protected (unless using CUDA)
#if !defined(__CUDACC__)
protected:
#endif

  /**
   * \brief This struct contains information used when merging fields.
   */
  struct FieldInformation
  {
    std::string topology;
    std::string association;
    int dtype;
    std::vector<std::string> components;
  };

  /**
   * \brief Check that the mesh inputs are valid and meet constraints. There must
   *        be 1 coordset/topology/matset. The coordset must be explicit and the
   *        topology must be unstructured and for now, non-polyhedral.
   *
   * \param inputs The mesh inputs.
   *
   * \return True if the inputs appear to be valid; False otherwise.
   */
  bool validInputs(const std::vector<MeshInput> &inputs, const conduit::Node &options) const
  {
    std::string topoName;
    if(options.has_child("topology")) topoName = options.fetch_existing("topology").as_string();

    for(size_t i = 0; i < inputs.size(); i++)
    {
      if(inputs[i].m_input == nullptr) return false;

      // If we did not specify which topology, make sure that there is only 1.
      const char *keys[] = {"coordsets", "topologies", "matsets"};
      if(topoName.empty())
      {
        for(int k = 0; k < 3; k++)
        {
          if(inputs[i].m_input->has_path(keys[k]))
          {
            const conduit::Node &n = inputs[i].m_input->fetch_existing(keys[k]);
            if(n.number_of_children() > 1) return false;
          }
        }
      }

      const conduit::Node &n_coordsets = inputs[i].m_input->fetch_existing("coordsets");
      const conduit::Node &n_coordset = n_coordsets[0];
      if(n_coordset["type"].as_string() != "explicit")
      {
        return false;
      }

      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node *n_topo = nullptr;
      if(topoName.empty())
        n_topo = n_topologies.child_ptr(0);
      else
        n_topo = n_topologies.fetch_ptr(topoName);
      if(n_topo->operator[]("type").as_string() != "unstructured")
      {
        return false;
      }
      // For now
      if(n_topo->operator[]("elements/shape").as_string() == "polyhedral")
      {
        return false;
      }

      // Require no nodeMap/nodeSlice or that they both be present.
      if(!((inputs[i].m_nodeMapView.size() == 0 && inputs[i].m_nodeSliceView.size() == 0) ||
           (inputs[i].m_nodeMapView.size() > 0 && inputs[i].m_nodeSliceView.size() > 0)))
      {
        return false;
      }
    }
    return true;
  }

  /**
   * \brief Merge a single input (copy it to the output).
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the merged mesh. 
   */
  void singleInput(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::copy<ExecSpace>(output, *(inputs[0].m_input));
  }

  /**
   * \brief Merge a multiple inputs.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the merged mesh. 
   */
  void mergeInputs(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    mergeCoordset(inputs, output);
    mergeTopology(inputs, output);
    mergeFields(inputs, output);
    mergeMatset(inputs, output);
  }

  /**
   * \brief Merge multiple coordsets into a single coordset. No node merging takes place.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeCoordset(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("mergeCoordset");
    const axom::IndexType totalNodes = countNodes(inputs);
    conduit::Node &n_newCoordsets = output["coordsets"];
    conduit::Node *n_newValuesPtr = nullptr;
    int nComps = 2;

    axom::IndexType offsets[3] = {0, 0, 0};
    const axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &coordsets = inputs[i].m_input->fetch_existing("coordsets");
      const conduit::Node &n_srcCoordset = coordsets[0];
      const conduit::Node &n_srcValues = n_srcCoordset.fetch_existing("values");

      const auto type = n_srcCoordset.fetch_existing("type").as_string();
      SLIC_ASSERT(type == "explicit");

      // Make all of the components the first time.
      if(i == 0)
      {
        bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

        conduit::Node &n_newCoordset = n_newCoordsets[n_srcCoordset.name()];
        n_newCoordset["type"] = "explicit";
        conduit::Node &n_newValues = n_newCoordset["values"];
        n_newValuesPtr = n_newCoordset.fetch_ptr("values");

        nComps = n_srcValues.number_of_children();
        for(int c = 0; c < nComps; c++)
        {
          const conduit::Node &n_srcComp = n_srcValues[c];
          conduit::Node &n_comp = n_newValues[n_srcComp.name()];
          n_comp.set_allocator(c2a.getConduitAllocatorID());
          n_comp.set(conduit::DataType(n_srcComp.dtype().id(), totalNodes));
        }
      }

      // Copy this input's coordinates into the new coordset.
      axom::mir::views::FloatNode_to_ArrayView(n_srcValues[0], [&](auto comp0) {
        using FloatType = typename decltype(comp0)::value_type;
        for(int c = 0; c < nComps; c++)
        {
          const conduit::Node &n_srcComp = n_srcValues[c];
          conduit::Node &n_comp = n_newValuesPtr->child(c);
          auto srcCompView = bputils::make_array_view<FloatType>(n_srcComp);
          auto compView = bputils::make_array_view<FloatType>(n_comp);

          axom::IndexType size =
            mergeCoordset_copy(inputs[i].m_nodeSliceView, offsets[c], compView, srcCompView);

          offsets[c] += size;
        }
      });
    }
  }

  /*!
   * \brief Assist setting merging coordset data.
   *
   * \param nodeSliceView The view that contains a node slice for the current input mesh.
   * \param offset The current write offset in the new coordset.
   * \param compView The view that exposes the current output coordinate component.
   * \param srcCompView The view that exposes the current output source coordinate component.
   *
   * \return The size of the data copied.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename DataArrayView>
  axom::IndexType mergeCoordset_copy(const axom::ArrayView<axom::IndexType> nodeSliceView,
                                     axom::IndexType offset,
                                     DataArrayView compView,
                                     DataArrayView srcCompView) const
  {
    axom::IndexType size = 0;
    if(nodeSliceView.size() > 0)
    {
      // Pull out specific nodes from the input.
      axom::for_all<ExecSpace>(
        nodeSliceView.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto sliceIndex = nodeSliceView[index];
          compView[offset + index] = srcCompView[sliceIndex];
        });
      size = nodeSliceView.size();
    }
    else
    {
      // Pull out all nodes from the input.
      axom::for_all<ExecSpace>(
        srcCompView.size(),
        AXOM_LAMBDA(axom::IndexType index) { compView[offset + index] = srcCompView[index]; });
      size = srcCompView.size();
    }
    return size;
  }

  /*!
   * \brief Count the number of nodes in the \a index'th input mesh.
   *
   * \param inputs The vector of input meshes.
   * \param index The index of the mesh to count.
   *
   * \return The number of nodes in the \a index input mesh.
   */
  axom::IndexType countNodes(const std::vector<MeshInput> &inputs, size_t index) const
  {
    SLIC_ASSERT(index < inputs.size());

    const conduit::Node &coordsets = inputs[index].m_input->fetch_existing("coordsets");
    const conduit::Node &coordset = coordsets[0];
    const auto type = coordset.fetch_existing("type").as_string();

    axom::IndexType nnodes = 0;
    if(inputs[index].m_nodeSliceView.size() > 0)
    {
      nnodes = inputs[index].m_nodeSliceView.size();
    }
    else
    {
      nnodes = conduit::blueprint::mesh::utils::coordset::length(coordset);
    }
    return nnodes;
  }

  /*!
   * \brief Count the number of nodes in all input meshes.
   *
   * \param inputs The vector of input meshes.
   *
   * \return The total number of nodes in the input meshes.
   */
  axom::IndexType countNodes(const std::vector<MeshInput> &inputs) const
  {
    axom::IndexType nodeTotal = 0;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      const auto nnodes = countNodes(inputs, i);
      nodeTotal += nnodes;
    }
    return nodeTotal;
  }

  /*!
   * \brief Count the number of zones in the \a index'th input mesh.
   *
   * \param inputs The vector of input meshes.
   * \param index The index of the mesh to count.
   *
   * \return The number of zones in the \a index input mesh.
   */
  axom::IndexType countZones(const std::vector<MeshInput> &inputs, size_t index) const
  {
    const conduit::Node &n_topologies = inputs[index].m_input->fetch_existing("topologies");
    const conduit::Node &n_topo = n_topologies[0];

    const conduit::Node &n_size = n_topo.fetch_existing("elements/sizes");
    axom::IndexType nzones = n_size.dtype().number_of_elements();
    return nzones;
  }

  /*!
   * \brief Count the number of nodes in all input meshes.
   *
   * \param inputs The vector of input meshes.
   * \param[out] totalConnLength The total connectivity length for all meshes.
   * \param[out] totalZones The total zones for all meshes.
   */
  void countZones(const std::vector<MeshInput> &inputs,
                  axom::IndexType &totalConnLength,
                  axom::IndexType &totalZones) const
  {
    totalConnLength = 0;
    totalZones = 0;
    axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_topo = n_topologies[0];
      const auto type = n_topo.fetch_existing("type").as_string();
      SLIC_ASSERT(type == "unstructured");

      const conduit::Node &n_size = n_topo.fetch_existing("elements/sizes");
      const auto nzones = n_size.dtype().number_of_elements();
      totalZones += nzones;

      // Instead of counting the total number of elements in the connectivity
      // array, we sum the sizes. This is needed because connectivity might
      // have unused elements in it.
      axom::IndexType connLength = 0;
      axom::mir::views::IndexNode_to_ArrayView(n_size, [&](auto sizesView) {
        connLength = sumArrayView(sizesView);
      });
      totalConnLength += connLength;
    }
  }

  /*!
   * \brief Sum the sizes array view and return the value.
   *
   * \param sizesView The view that contains the sizes.
   *
   * \note This is implemented as a template method because we can't use
   *       axom::for_all in a lambda.
   */
  template <typename ViewType>
  axom::IndexType sumArrayView(ViewType view) const
  {
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    using value_type = typename ViewType::value_type;
    RAJA::ReduceSum<reduce_policy, value_type> sum(0);
    axom::for_all<ExecSpace>(
      view.size(),
      AXOM_LAMBDA(axom::IndexType index) { sum += view[index]; });
    return static_cast<axom::IndexType>(sum.get());
  }

  /**
   * \brief Merge multiple topologies into a single topology.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeTopology(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("mergeTopology");
    axom::IndexType totalConnLen = 0, totalZones = 0;
    countZones(inputs, totalConnLen, totalZones);
    conduit::Node &n_newTopologies = output["topologies"];

    // Check whether there are mixed shapes.
    std::map<std::string, int> shape_map;
    const axom::IndexType n = static_cast<axom::IndexType>(inputs.size());
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_srcTopo = n_topologies[0];
      const auto type = n_srcTopo.fetch_existing("type").as_string();
      const auto shape = n_srcTopo.fetch_existing("elements/shape").as_string();
      SLIC_ASSERT(type == "unstructured");
      if(shape == "mixed")
      {
        const conduit::Node &n_shape_map = n_srcTopo.fetch_existing("elements/shape_map");
        for(int s = 0; s < n_shape_map.number_of_children(); s++)
        {
          const std::string sname = n_shape_map[s].name();
          const auto id = axom::mir::views::shapeNameToID(sname);
          SLIC_ASSERT(id == n_shape_map[s].to_int());
          shape_map[sname] = id;
        }
      }
      else
      {
        shape_map[shape] = axom::mir::views::shapeNameToID(shape);
      }
    }

    // If there are polygon shapes then assume that the rest of the shapes
    // are 2D and should be promoted to polygons.
    if(shape_map.find("polygonal") != shape_map.end() && shape_map.size() > 1)
    {
      shape_map.clear();
      shape_map["polygonal"] = axom::mir::views::shapeNameToID("polygonal");
    }

    conduit::Node *n_newTopoPtr = nullptr;
    axom::IndexType connOffset = 0, sizesOffset = 0, shapesOffset = 0, coordOffset = 0;
    for(axom::IndexType i = 0; i < n; i++)
    {
      const conduit::Node &n_topologies = inputs[i].m_input->fetch_existing("topologies");
      const conduit::Node &n_srcTopo = n_topologies[0];

      const std::string srcShape = n_srcTopo.fetch_existing("elements/shape").as_string();
      const conduit::Node &n_srcConn = n_srcTopo.fetch_existing("elements/connectivity");
      const conduit::Node &n_srcSizes = n_srcTopo.fetch_existing("elements/sizes");
      const conduit::Node &n_srcOffsets = n_srcTopo.fetch_existing("elements/offsets");

      // Make all of the elements the first time.
      if(i == 0)
      {
        bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

        conduit::Node &n_newTopo = n_newTopologies[n_srcTopo.name()];
        n_newTopoPtr = n_newTopologies.fetch_ptr(n_srcTopo.name());
        n_newTopo["type"] = "unstructured";
        n_newTopo["coordset"] = n_srcTopo["coordset"].as_string();

        conduit::Node &n_newConn = n_newTopo["elements/connectivity"];
        n_newConn.set_allocator(c2a.getConduitAllocatorID());
        n_newConn.set(conduit::DataType(n_srcConn.dtype().id(), totalConnLen));

        conduit::Node &n_newSizes = n_newTopo["elements/sizes"];
        n_newSizes.set_allocator(c2a.getConduitAllocatorID());
        n_newSizes.set(conduit::DataType(n_srcSizes.dtype().id(), totalZones));

        conduit::Node &n_newOffsets = n_newTopo["elements/offsets"];
        n_newOffsets.set_allocator(c2a.getConduitAllocatorID());
        n_newOffsets.set(conduit::DataType(n_srcConn.dtype().id(), totalZones));

        if(shape_map.size() > 1)
        {
          n_newTopo["elements/shape"] = "mixed";

          // Build a new shape map in the new topology.
          conduit::Node &n_shape_map = n_newTopo["elements/shape_map"];
          for(auto it = shape_map.begin(); it != shape_map.end(); it++)
            n_shape_map[it->first] = it->second;

          conduit::Node &n_newShapes = n_newTopo["elements/shapes"];
          n_newShapes.set_allocator(c2a.getConduitAllocatorID());
          n_newShapes.set(conduit::DataType(n_srcConn.dtype().id(), totalZones));
        }
        else
        {
          n_newTopo["elements/shape"] = shape_map.begin()->first;
        }
      }

      // Copy this input's connectivity into the new topology.
      axom::mir::views::IndexNode_to_ArrayView_same(
        n_srcConn,
        n_srcSizes,
        n_srcOffsets,
        [&](auto srcConnView, auto srcSizesView, auto srcOffsetsView) {
          using ConnType = typename decltype(srcConnView)::value_type;
          conduit::Node &n_newConn = n_newTopoPtr->fetch_existing("elements/connectivity");
          auto connView = bputils::make_array_view<ConnType>(n_newConn);

          // Copy the relevant connectivity from srcConnView. Also compute how
          // many elements were used.
          mergeTopology_copy(inputs[i].m_nodeMapView,
                             connOffset,
                             coordOffset,
                             connView,
                             srcConnView,
                             srcSizesView,
                             srcOffsetsView);

          connOffset += srcConnView.size();
          coordOffset += countNodes(inputs, static_cast<size_t>(i));
        });

      // Copy this input's sizes into the new topology.
      axom::mir::views::IndexNode_to_ArrayView(n_srcSizes, [&](auto srcSizesView) {
        using ConnType = typename decltype(srcSizesView)::value_type;
        conduit::Node &n_newSizes = n_newTopoPtr->fetch_existing("elements/sizes");
        auto sizesView = bputils::make_array_view<ConnType>(n_newSizes);

        mergeTopology_copy_sizes(sizesOffset, sizesView, srcSizesView);

        sizesOffset += srcSizesView.size();
      });

      // If the shape map contains multiple shapes then we're making an
      // elements/shapes array.
      if(shape_map.size() > 1)
      {
        // Copy shape information if it exists.
        if(n_srcTopo.has_path("elements/shapes"))
        {
          const conduit::Node &n_srcShapes = n_srcTopo.fetch_existing("elements/shapes");

          axom::mir::views::IndexNode_to_ArrayView(n_srcShapes, [&](auto srcShapesView) {
            using ConnType = typename decltype(srcShapesView)::value_type;
            conduit::Node &n_newShapes = n_newTopoPtr->fetch_existing("elements/shapes");
            auto shapesView = bputils::make_array_view<ConnType>(n_newShapes);
            // Copy all sizes from the input.
            mergeTopology_copy_shapes(shapesOffset, shapesView, srcShapesView);
            shapesOffset += srcShapesView.size();
          });
        }
        else
        {
          // Fill in shape information. There is no source shape information. Use
          // sizes to get the number of zones.
          const conduit::Node &n_srcSizes = n_srcTopo.fetch_existing("elements/sizes");
          axom::IndexType nz = n_srcSizes.dtype().number_of_elements();
          conduit::Node &n_newShapes = n_newTopoPtr->fetch_existing("elements/shapes");
          axom::mir::views::IndexNode_to_ArrayView(n_newShapes, [&](auto shapesView) {
            const int shapeId = axom::mir::views::shapeNameToID(srcShape);
            mergeTopology_default_shapes(shapesOffset, shapesView, nz, shapeId);
            shapesOffset += nz;
          });
        }
      }
    }

    // Make new offsets from the sizes.
    conduit::Node &n_newSizes = n_newTopoPtr->fetch_existing("elements/sizes");
    axom::mir::views::IndexNode_to_ArrayView(n_newSizes, [&](auto sizesView) {
      using ConnType = typename decltype(sizesView)::value_type;
      conduit::Node &n_newOffsets = n_newTopoPtr->fetch_existing("elements/offsets");
      auto offsetsView = bputils::make_array_view<ConnType>(n_newOffsets);
      axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);
    });
  }

  /*!
   * \brief Assist copying topology connectivity to the merged topology.
   *
   * \param nodeMapView The node map.
   * \param connOffset The write offset in the new connectivity.
   * \param coordOffset The current mesh's coordinate offset in the new coordinates.
   * \param connView The view that contains the new merged connectivity.
   * \param srcConnView The view that contains the source connectivity.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename ConnectivityView>
  void mergeTopology_copy(axom::ArrayView<axom::IndexType> nodeMapView,
                          axom::IndexType connOffset,
                          axom::IndexType coordOffset,
                          ConnectivityView connView,
                          ConnectivityView srcConnView,
                          ConnectivityView srcSizesView,
                          ConnectivityView srcOffsetsView) const
  {
    using value_type = typename ConnectivityView::value_type;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Compress out any gaps in the offsets by making new offsets from the sizes.
    // This makes the code able to handle meshes that have gaps in its connectivity array.
    axom::Array<value_type> actualOffsets(srcSizesView.size(), srcSizesView.size(), allocatorID);
    auto actualOffsetsView = actualOffsets.view();
    axom::exclusive_scan<ExecSpace>(srcSizesView, actualOffsetsView);

    if(nodeMapView.size() > 0)
    {
      // Copy all zones from the input but map the nodes to new values.
      // The supplied nodeMap is assumed to be a mapping from the current
      // node connectivity to the merged node connectivity.
      axom::for_all<ExecSpace>(
        srcSizesView.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto destOffset = connOffset + actualOffsetsView[index];
          const auto srcOffset = srcOffsetsView[index];
          for(value_type j = 0; j < srcSizesView[index]; j++)
          {
            const auto nodeId = srcConnView[srcOffset + j];
            const auto newNodeId = nodeMapView[nodeId];
            connView[destOffset + j] = newNodeId;
          }
        });
    }
    else
    {
      axom::for_all<ExecSpace>(
        srcSizesView.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto destOffset = connOffset + actualOffsetsView[index];
          const auto srcOffset = srcOffsetsView[index];
          for(value_type j = 0; j < srcSizesView[index]; j++)
          {
            connView[destOffset + j] = coordOffset + srcConnView[srcOffset + j];
          }
        });
    }
  }

  /*!
   * \brief Assist copying topology sizes to the merged topology.
   *
   * \param sizesOffset The write offset for sizes in the new connectivity.
   * \param sizesView The view that contains sizes for the new connectivity.
   * \param srcSizesView The view that contains sizes for the input mesh.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerView>
  void mergeTopology_copy_sizes(axom::IndexType sizesOffset,
                                IntegerView sizesView,
                                IntegerView srcSizesView) const
  {
    axom::for_all<ExecSpace>(
      srcSizesView.size(),
      AXOM_LAMBDA(axom::IndexType index) { sizesView[sizesOffset + index] = srcSizesView[index]; });
  }

  /*!
   * \brief Copy shapes from the source mesh to the merged mesh.
   *
   * \param shapesOffset The write offset for the shapes.
   * \param shapesView The view that exposes shapes for the merged mesh.
   * \param srcShapesView The view that exposes shapes for the source mesh.
   */
  template <typename IntegerView>
  void mergeTopology_copy_shapes(axom::IndexType shapesOffset,
                                 IntegerView shapesView,
                                 IntegerView srcShapesView) const
  {
    axom::for_all<ExecSpace>(
      srcShapesView.size(),
      AXOM_LAMBDA(axom::IndexType index) { shapesView[shapesOffset + index] = srcShapesView[index]; });
  }

  /*!
   * \brief Set shapes in the merged mesh to a specific shape.
   *
   * \param shapesOffset The write offset for the shapes.
   * \param shapesView The view that exposes shapes for the merged mesh.
   * \param srcShapesView The view that exposes shapes for the source mesh.
   */
  template <typename IntegerView>
  void mergeTopology_default_shapes(axom::IndexType shapesOffset,
                                    IntegerView shapesView,
                                    axom::IndexType nzones,
                                    int shapeId) const
  {
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType index) { shapesView[shapesOffset + index] = shapeId; });
  }

  /*!
   * \brief Merge fields that exist on the various mesh inputs. Zero-fill values
   *        where a field does not exist in an input.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  void mergeFields(const std::vector<MeshInput> &inputs, conduit::Node &output) const
  {
    namespace bputils = axom::mir::utilities::blueprint;
    AXOM_ANNOTATE_SCOPE("mergeFields");
    axom::IndexType totalNodes = countNodes(inputs);
    axom::IndexType totalConnLen = 0, totalZones = 0;
    countZones(inputs, totalConnLen, totalZones);
    const axom::IndexType n = static_cast<axom::IndexType>(inputs.size());

    // Determine whether any inputs have fields.
    bool hasFields = false;
    for(axom::IndexType i = 0; i < n; i++)
    {
      hasFields |= inputs[i].m_input->has_child("fields");
    }

    if(hasFields)
    {
      // Make field information in case some inputs do not have the field.
      std::map<std::string, FieldInformation> fieldInfo;
      for(axom::IndexType i = 0; i < n; i++)
      {
        if(inputs[i].m_input->has_child("fields"))
        {
          const conduit::Node &n_fields = inputs[i].m_input->fetch_existing("fields");
          for(conduit::index_t c = 0; c < n_fields.number_of_children(); c++)
          {
            const conduit::Node &n_field = n_fields[c];
            const conduit::Node &n_values = n_field.fetch_existing("values");
            FieldInformation fi;
            fi.topology = n_field.fetch_existing("topology").as_string();
            fi.association = n_field.fetch_existing("association").as_string();
            if(n_values.number_of_children() > 0)
            {
              for(conduit::index_t comp = 0; comp < n_values.number_of_children(); comp++)
              {
                fi.components.push_back(n_values[comp].name());
                fi.dtype = n_values[comp].dtype().id();
              }
            }
            else
            {
              fi.dtype = n_values.dtype().id();
            }
            fieldInfo[n_field.name()] = fi;
          }
        }
      }

      // Make new fields
      bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;
      conduit::Node &n_newFields = output["fields"];
      for(auto it = fieldInfo.begin(); it != fieldInfo.end(); it++)
      {
        conduit::Node &n_newField = n_newFields[it->first];
        n_newField["association"] = it->second.association;
        n_newField["topology"] = it->second.topology;
        conduit::Node &n_values = n_newField["values"];
        if(it->second.components.empty())
        {
          // Scalar
          conduit::Node &n_values = n_newField["values"];
          n_values.set_allocator(c2a.getConduitAllocatorID());
          const std::string srcPath("fields/" + it->first + "/values");
          if(it->second.association == "element")
          {
            n_values.set(conduit::DataType(it->second.dtype, totalZones));
            copyZonal(inputs, n_values, srcPath);
          }
          else if(it->second.association == "vertex")
          {
            n_values.set(conduit::DataType(it->second.dtype, totalNodes));
            copyNodal(inputs, n_values, srcPath);
          }
        }
        else
        {
          // Vector
          for(size_t ci = 0; ci < it->second.components.size(); ci++)
          {
            conduit::Node &n_comp = n_values[it->second.components[ci]];
            n_comp.set_allocator(c2a.getConduitAllocatorID());
            const std::string srcPath("fields/" + it->first + "/values/" + it->second.components[ci]);
            if(it->second.association == "element")
            {
              n_comp.set(conduit::DataType(it->second.dtype, totalZones));
              copyZonal(inputs, n_comp, srcPath);
            }
            else if(it->second.association == "vertex")
            {
              n_comp.set(conduit::DataType(it->second.dtype, totalNodes));
              copyNodal(inputs, n_comp, srcPath);
            }
          }
        }
      }
    }
  }

  /**
   * \brief Copy zonal field data into a Conduit node.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] n_values The node will be populated with data values from the field inputs.
   * \param srcPath The path to the source data in each input node.
   */
  void copyZonal(const std::vector<MeshInput> &inputs,
                 conduit::Node &n_values,
                 const std::string &srcPath) const
  {
    axom::IndexType offset = 0;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      const auto nzones = countZones(inputs, i);

      if(inputs[i].m_input->has_path(srcPath))
      {
        const conduit::Node &n_src_values = inputs[i].m_input->fetch_existing(srcPath);
        axom::mir::views::Node_to_ArrayView(n_values, n_src_values, [&](auto destView, auto srcView) {
          copyZonal_copy(nzones, offset, destView, srcView);
        });
      }
      else
      {
        axom::mir::views::Node_to_ArrayView(n_values, [&](auto destView) {
          fillValues(nzones, offset, destView);
        });
      }
      offset += nzones;
    }
  }

  /*!
   * \brief Copy zonal data from src to dest.
   *
   * \param nzones The number of zones.
   * \param offset The current write offset.
   * \param destView The view that exposes the new merged field.
   * \param srcView The view that exposes the source field.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename DestView, typename SrcView>
  void copyZonal_copy(axom::IndexType nzones,
                      axom::IndexType offset,
                      DestView destView,
                      SrcView srcView) const
  {
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType index) { destView[offset + index] = srcView[index]; });
  }

  /*!
   * \brief Fill data in dest.
   *
   * \param nvalues The number of values.
   * \param offset The current write offset.
   * \param destView The view that exposes the new merged field.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename DestView>
  void fillValues(axom::IndexType nvalues, axom::IndexType offset, DestView destView) const
  {
    axom::for_all<ExecSpace>(
      nvalues,
      AXOM_LAMBDA(axom::IndexType index) { destView[offset + index] = 0; });
  }

  /**
   * \brief Copy nodal field data into a Conduit node.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] n_values The node will be populated with data values from the field inputs.
   * \param srcPath The path to the source data in each input node.
   */
  void copyNodal(const std::vector<MeshInput> &inputs,
                 conduit::Node &n_values,
                 const std::string &srcPath) const
  {
    axom::IndexType offset = 0;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      const auto nnodes = countNodes(inputs, i);

      if(inputs[i].m_input->has_path(srcPath))
      {
        const conduit::Node &n_src_values = inputs[i].m_input->fetch_existing(srcPath);

        axom::mir::views::Node_to_ArrayView(n_src_values, n_values, [&](auto srcView, auto destView) {
          copyNodal_copy(inputs[i].m_nodeSliceView, nnodes, offset, destView, srcView);
        });
      }
      else
      {
        axom::mir::views::Node_to_ArrayView(n_values, [&](auto destView) {
          fillValues(nnodes, offset, destView);
        });
      }
      offset += nnodes;
    }
  }

  /*!
   * \brief Copy nodal data from src to dest.
   *
   * \param nodeSliceView The nodes we're pulling out (if populated).
   * \param nnodes The number of nodes.
   * \param offset The current write offset.
   * \param destView The view that exposes the new merged field.
   * \param srcView The view that exposes the source field.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename SrcViewType, typename DestViewType>
  void copyNodal_copy(axom::ArrayView<axom::IndexType> nodeSliceView,
                      axom::IndexType nnodes,
                      axom::IndexType offset,
                      SrcViewType destView,
                      DestViewType srcView) const
  {
    if(nodeSliceView.empty())
    {
      axom::for_all<ExecSpace>(
        nnodes,
        AXOM_LAMBDA(axom::IndexType index) { destView[offset + index] = srcView[index]; });
    }
    else
    {
      axom::for_all<ExecSpace>(
        nodeSliceView.size(),
        AXOM_LAMBDA(axom::IndexType index) {
          const auto nodeId = nodeSliceView[index];
          destView[offset + index] = srcView[nodeId];
        });
    }
  }

  /**
   * \brief Merge matsets that exist on the various mesh inputs.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  virtual void mergeMatset(const std::vector<MeshInput> &AXOM_UNUSED_PARAM(inputs),
                           conduit::Node &AXOM_UNUSED_PARAM(output)) const
  {
    // Do nothing.
  }
};

/*!
 * \brief Dispatches any type of data as any type of matset.
 */
class DispatchAnyMatset
{
public:
  /*!
   * \brief Takes input matset data as Conduit nodes whose data can be various
   *        types and dispatches the data to array views that are passed to a
   *        lambda function that processes them. These nodes contain data for
   *        the merged material being created.
   *
   * \tparam FuncType A callable object/function/lambda.
   *
   * \param n_material_ids A node containing material ids.
   * \param n_sizes A node containing sizes.
   * \param n_offsets A node containing offsets.
   * \param n_indices A node containing indices.
   * \param n_volume_fractions A node containing volume fractions.
   * \param func The function to invoke on the array views.
   */
  template <typename FuncType>
  void execute(conduit::Node &n_material_ids,
               conduit::Node &n_sizes,
               conduit::Node &n_offsets,
               conduit::Node &n_indices,
               conduit::Node &n_volume_fractions,
               FuncType &&func)
  {
    // Support various types of material data.
    axom::mir::views::IndexNode_to_ArrayView_same(
      n_material_ids,
      n_sizes,
      n_offsets,
      n_indices,
      [&](auto materialIdsView, auto sizesView, auto offsetsView, auto indicesView) {
        axom::mir::views::FloatNode_to_ArrayView(n_volume_fractions, [&](auto volumeFractionsView) {
          // Invoke a function that can use these views.
          func(materialIdsView, sizesView, offsetsView, indicesView, volumeFractionsView);
        });
      });
  }

  /*!
   * \brief Takes a node and turns it into a matset view and sends that view
   *        to the input function.
   *
   * \tparam FuncType A callable object/function/lambda.
   *
   * \param n_matset A node containing any kind of matset.
   * \param func The function to invoke on the matset view.
   */
  template <typename FuncType>
  void dispatchMatset(conduit::Node &n_matset, FuncType &&func)
  {
    axom::mir::views::dispatch_material(n_matset, [&](auto matsetView) { func(matsetView); });
  }
};

/*!
 * \brief Dispatches data of known type to a unibuffer matset.
 */
template <typename IntElement, typename FloatElement, size_t MAXMATERIALS = 20>
class DispatchTypedUnibufferMatset
{
public:
  /*!
   * \brief Takes input matset data as Conduit nodes and turn them into
   *        specific view types since we know the types of data being used.
   *
   * \tparam FuncType A callable object/function/lambda.
   *
   * \param n_material_ids A node containing material ids.
   * \param n_sizes A node containing sizes.
   * \param n_offsets A node containing offsets.
   * \param n_indices A node containing indices.
   * \param n_volume_fractions A node containing volume fractions.
   * \param func The function to invoke on the array views.
   */
  template <typename FuncType>
  void execute(conduit::Node &n_material_ids,
               conduit::Node &n_sizes,
               conduit::Node &n_offsets,
               conduit::Node &n_indices,
               conduit::Node &n_volume_fractions,
               FuncType &&func)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    auto materialIdsView = bputils::make_array_view<IntElement>(n_material_ids);
    auto sizesView = bputils::make_array_view<IntElement>(n_sizes);
    auto offsetsView = bputils::make_array_view<IntElement>(n_offsets);
    auto indicesView = bputils::make_array_view<IntElement>(n_indices);
    auto volumeFractionsView = bputils::make_array_view<FloatElement>(n_volume_fractions);

    func(materialIdsView, sizesView, offsetsView, indicesView, volumeFractionsView);
  }

  /*!
   * \brief Takes a node and turns it into a matset view and sends that view
   *        to the input function.
   *
   * \tparam FuncType A callable object/function/lambda.
   *
   * \param n_matset A node containing any kind of matset.
   * \param func The function to invoke on the matset view.
   */
  template <typename FuncType>
  void dispatchMatset(conduit::Node &n_matset, FuncType &&func)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    // We know the types. Make views explicitly.
    auto material_ids = bputils::make_array_view<IntElement>(n_matset["material_ids"]);
    auto sizes = bputils::make_array_view<IntElement>(n_matset["sizes"]);
    auto offsets = bputils::make_array_view<IntElement>(n_matset["offsets"]);
    auto indices = bputils::make_array_view<IntElement>(n_matset["indices"]);
    auto volume_fractions = bputils::make_array_view<FloatElement>(n_matset["volume_fractions"]);
    // We know we're making a unibuffer matset.
    axom::mir::views::UnibufferMaterialView<IntElement, FloatElement, MAXMATERIALS> matsetView;
    matsetView.set(material_ids, volume_fractions, sizes, offsets, indices);
    // Use it.
    func(matsetView);
  }
};

/**
 * \brief Merge multiple unstructured Blueprint meshes (with matsets) through MeshInput.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 * \tparam MaterialDispatch A policy that helps determine how materials are dispatched.
 *                          The default is to handle any type of matset but that can
 *                          generate a lot of code.
 *
 * \note The input meshes must currently contain a single coordset/topology/matset.
 */
template <typename ExecSpace, typename MaterialDispatch = DispatchAnyMatset>
class MergeMeshesAndMatsets : public MergeMeshes<ExecSpace>
{
public:
// The following members are private (unless using CUDA)
#if !defined(__CUDACC__)
private:
#endif
  /**
   * \brief Merge matsets that exist on the various mesh inputs.
   *
   * \param inputs A vector of inputs to be merged.
   * \param[out] output The node that will contain the output mesh.
   */
  virtual void mergeMatset(const std::vector<MeshInput> &inputs, conduit::Node &output) const override
  {
    AXOM_ANNOTATE_SCOPE("mergeMatset");
    namespace bputils = axom::mir::utilities::blueprint;
    bputils::ConduitAllocateThroughAxom<ExecSpace> c2a;

    // Make a pass through the inputs and make a list of the material names.
    bool hasMatsets = false, defaultMaterial = false;
    int nmats = 0;
    std::map<std::string, int> allMats;
    std::string matsetName, topoName;
    for(size_t i = 0; i < inputs.size(); i++)
    {
      if(inputs[i].m_input->has_path("matsets"))
      {
        conduit::Node &n_matsets = inputs[i].m_input->fetch_existing("matsets");
        conduit::Node &n_matset = n_matsets[0];
        matsetName = n_matset.name();
        topoName = n_matset.fetch_existing("topology").as_string();
        auto matInfo = axom::mir::views::materials(n_matset);
        for(const auto &info : matInfo)
        {
          if(allMats.find(info.name) == allMats.end())
          {
            allMats[info.name] = nmats++;
          }
        }
        hasMatsets = true;
      }
      else
      {
        defaultMaterial = true;
      }
    }

    if(hasMatsets)
    {
      MaterialDispatch disp;
      auto *This = this;

      // One or more inputs did not have a matset.
      if(defaultMaterial) allMats["default"] = nmats++;

      // Make a pass through the matsets to determine the overall storage.
      axom::IndexType totalZones = 0, totalMatCount = 0;
      int itype, ftype;
      {
        AXOM_ANNOTATE_SCOPE("sizes");
        for(size_t i = 0; i < inputs.size(); i++)
        {
          const auto nzones = this->countZones(inputs, i);
          totalZones += nzones;

          if(inputs[i].m_input->has_path("matsets"))
          {
            conduit::Node &n_matsets = inputs[i].m_input->fetch_existing("matsets");
            conduit::Node &n_matset = n_matsets[0];
            axom::IndexType matCount = 0;
            disp.dispatchMatset(n_matset, [&](auto matsetView) {
              // Figure out the types to use for storing the data.
              using IType = typename decltype(matsetView)::IndexType;
              using FType = typename decltype(matsetView)::FloatType;
              itype = bputils::cpp2conduit<IType>::id;
              ftype = bputils::cpp2conduit<FType>::id;

              matCount = This->mergeMatset_count(matsetView, nzones);
            });
            totalMatCount += matCount;
          }
          else
          {
            totalMatCount += nzones;
          }
        }
      }

      // Allocate
      AXOM_ANNOTATE_BEGIN("allocate");
      conduit::Node &n_newMatset = output["matsets/" + matsetName];
      n_newMatset["topology"] = topoName;
      conduit::Node &n_volume_fractions = n_newMatset["volume_fractions"];
      n_volume_fractions.set_allocator(c2a.getConduitAllocatorID());
      n_volume_fractions.set(conduit::DataType(ftype, totalMatCount));

      conduit::Node &n_material_ids = n_newMatset["material_ids"];
      n_material_ids.set_allocator(c2a.getConduitAllocatorID());
      n_material_ids.set(conduit::DataType(itype, totalMatCount));

      conduit::Node &n_sizes = n_newMatset["sizes"];
      n_sizes.set_allocator(c2a.getConduitAllocatorID());
      n_sizes.set(conduit::DataType(itype, totalZones));

      conduit::Node &n_offsets = n_newMatset["offsets"];
      n_offsets.set_allocator(c2a.getConduitAllocatorID());
      n_offsets.set(conduit::DataType(itype, totalZones));

      conduit::Node &n_indices = n_newMatset["indices"];
      n_indices.set_allocator(c2a.getConduitAllocatorID());
      n_indices.set(conduit::DataType(itype, totalMatCount));
      AXOM_ANNOTATE_END("allocate");

      {
        AXOM_ANNOTATE_SCOPE("populate");

        // Make material_map.
        conduit::Node &n_material_map = n_newMatset["material_map"];
        for(auto it = allMats.begin(); it != allMats.end(); it++)
          n_material_map[it->first] = it->second;

        // Populate
        disp.execute(n_material_ids,
                     n_sizes,
                     n_offsets,
                     n_indices,
                     n_volume_fractions,
                     [&](auto materialIdsView,
                         auto sizesView,
                         auto offsetsView,
                         auto indicesView,
                         auto volumeFractionsView) {
                       // Fill in sizes array.
                       axom::IndexType zOffset = 0;
                       for(size_t i = 0; i < inputs.size(); i++)
                       {
                         const auto nzones = This->countZones(inputs, i);

                         if(inputs[i].m_input->has_child("matsets"))
                         {
                           conduit::Node &n_matsets = inputs[i].m_input->fetch_existing("matsets");
                           conduit::Node &n_matset = n_matsets[0];

                           disp.dispatchMatset(n_matset, [&](auto matsetView) {
                             This->mergeMatset_sizes(matsetView, sizesView, nzones, zOffset);
                           });
                         }
                         else
                         {
                           This->mergeMatset_sizes1(sizesView, nzones, zOffset);
                         }
                         zOffset += nzones;
                       }
                       // Make offsets.
                       axom::exclusive_scan<ExecSpace>(sizesView, offsetsView);

                       // Make indices.
                       This->mergeMatset_indices(indicesView, totalMatCount);

                       // Fill in material info.
                       zOffset = 0;
                       for(size_t i = 0; i < inputs.size(); i++)
                       {
                         const auto nzones = This->countZones(inputs, i);

                         if(inputs[i].m_input->has_child("matsets"))
                         {
                           conduit::Node &n_matsets = inputs[i].m_input->fetch_existing("matsets");
                           conduit::Node &n_matset = n_matsets[0];

                           disp.dispatchMatset(n_matset, [&](auto matsetView) {
                             This->mergeMatset_copy(n_matset,
                                                    allMats,
                                                    materialIdsView,
                                                    offsetsView,
                                                    volumeFractionsView,
                                                    matsetView,
                                                    nzones,
                                                    zOffset);
                           });
                         }
                         else
                         {
                           const int dmat = allMats["default"];
                           This->mergeMatset_default(materialIdsView,
                                                     offsetsView,
                                                     volumeFractionsView,
                                                     dmat,
                                                     nzones,
                                                     zOffset);
                         }
                         zOffset += nzones;
                       }
                     });
      }
    }  // if hasMatsets
  }

  /*!
   * \brief Assist in counting the total material elements needed for the input matset.
   *
   * \param matsetView The view that wraps the material data.
   * \param nzones The number of zones.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename MatsetView>
  axom::IndexType mergeMatset_count(MatsetView matsetView, axom::IndexType nzones) const
  {
    using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;
    RAJA::ReduceSum<reduce_policy, axom::IndexType> matCount_reduce(0);
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto nmats = matsetView.numberOfMaterials(zoneIndex);
        matCount_reduce += nmats;
      });
    return matCount_reduce.get();
  }

  /*!
   * \brief Assist in setting sizes for the new matset.
   *
   * \param matsetView The view that wraps the material data.
   * \param sizesView The view that exposes the new matset sizes.
   * \param nzones The number of zones in the current input.
   * \param zOffset The current offset in the merged output.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename MatsetView, typename IntegerArrayView>
  void mergeMatset_sizes(const MatsetView matsetView,
                         IntegerArrayView sizesView,
                         axom::IndexType nzones,
                         axom::IndexType zOffset) const
  {
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        sizesView[zOffset + zoneIndex] = matsetView.numberOfMaterials(zoneIndex);
      });
  }

  /*!
   * \brief Assist in setting sizes for the new matset.
   *
   * \param sizesView The view that exposes the new matset sizes.
   * \param nzones The number of zones in the current input.
   * \param zOffset The current offset in the merged output.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerArrayView>
  void mergeMatset_sizes1(IntegerArrayView sizesView,
                          axom::IndexType nzones,
                          axom::IndexType zOffset) const
  {
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) { sizesView[zOffset + zoneIndex] = 1; });
  }

  /*!
   * \brief Assist in setting indices for the new matset.
   *
   * \param indicesView The view that exposes the new matset indices.
   * \param totalMatCount The total number of elements in the material_ids or volume_fractions array.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerArrayView>
  void mergeMatset_indices(IntegerArrayView indicesView, axom::IndexType totalMatCount) const
  {
    axom::for_all<ExecSpace>(
      totalMatCount,
      AXOM_LAMBDA(axom::IndexType index) { indicesView[index] = index; });
  }

  /*!
   * \brief Assist in copying matset data into the new merged matset.
   *
   * \param n_matset A Node that contains the matset.
   * \param allMats A map of material numbers to merged material numbers.
   * \param materialIdsView The view that exposes the merged material ids.
   * \param offsetsView The view that exposes the merged material offsets.
   * \param volumeFractionsView The view that exposes the merged volume fractions.
   * \param matsetView The matset view that contains the source matset data.
   * \param nzones The number of zones in the current input.
   * \param zOffset The current offset in the merged output.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerView, typename FloatView, typename MatsetView>
  void mergeMatset_copy(const conduit::Node &n_matset,
                        const std::map<std::string, int> &allMats,
                        IntegerView materialIdsView,
                        IntegerView offsetsView,
                        FloatView volumeFractionsView,
                        MatsetView matsetView,
                        axom::IndexType nzones,
                        axom::IndexType zOffset) const
  {
    using IDList = typename decltype(matsetView)::IDList;
    using VFList = typename decltype(matsetView)::VFList;
    using MatID = typename decltype(matsetView)::IndexType;

    // Make some maps for renumbering material numbers.
    const auto localMaterialMap = axom::mir::views::materials(n_matset);
    std::map<MatID, MatID> localToAll;
    for(const auto &info : localMaterialMap)
    {
      const auto it = allMats.find(info.name);
      SLIC_ASSERT(it != allMats.end());
      MatID matno = it->second;
      localToAll[info.number] = matno;
    }
    std::vector<MatID> localVec, allVec;
    for(auto it = localToAll.begin(); it != localToAll.end(); it++)
    {
      localVec.push_back(it->first);
      allVec.push_back(it->second);
    }
    // Put maps on device.
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<MatID> local(localVec.size(), localVec.size(), allocatorID);
    axom::Array<MatID> all(allVec.size(), allVec.size(), allocatorID);
    axom::copy(local.data(), localVec.data(), sizeof(MatID) * local.size());
    axom::copy(all.data(), allVec.data(), sizeof(MatID) * all.size());
    const auto localView = local.view();
    const auto allView = all.view();

    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        // Get this zone's materials.
        IDList ids;
        VFList vfs;
        matsetView.zoneMaterials(zoneIndex, ids, vfs);

        // Store the materials in the new material.
        const auto zoneStart = offsetsView[zOffset + zoneIndex];
        for(axom::IndexType mi = 0; mi < ids.size(); mi++)
        {
          const auto destIndex = zoneStart + mi;
          volumeFractionsView[destIndex] = vfs[mi];

          // Get the index of the material number in the local map.
          const auto mapIndex = axom::mir::utilities::bsearch(ids[mi], localView);
          SLIC_ASSERT(mapIndex != -1);
          // We'll store the all materials number.
          const auto allMatno = allView[mapIndex];
          materialIdsView[destIndex] = allMatno;
        }
      });
  }

  /*!
   * \brief Assist setting default material data when the input mesh has no material.
   *
   * \param materialIdsView The view that exposes the merged material ids.
   * \param offsetsView The view that exposes the merged material offsets.
   * \param volumeFractionsView The view that exposes the merged volume fractions.
   * \param matno The material number to use for these zones.
   * \param nzones The number of zones in the current input.
   * \param zOffset The current offset in the merged output.
   *
   * \note This method was broken out into a template member method since nvcc
   *       would not instantiate the lambda for axom::for_all() from an anonymous
   *       lambda.
   */
  template <typename IntegerView, typename FloatView>
  void mergeMatset_default(IntegerView materialIdsView,
                           IntegerView offsetsView,
                           FloatView volumeFractionsView,
                           int matno,
                           axom::IndexType nzones,
                           axom::IndexType zOffset) const
  {
    axom::for_all<ExecSpace>(
      nzones,
      AXOM_LAMBDA(axom::IndexType zoneIndex) {
        const auto zoneStart = offsetsView[zOffset + zoneIndex];
        volumeFractionsView[zoneStart] = 1;
        materialIdsView[zoneStart] = matno;
      });
  }
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
