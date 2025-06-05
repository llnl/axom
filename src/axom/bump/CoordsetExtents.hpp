// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BUMP_COORDSET_EXTENTS_HPP_
#define AXOM_BUMP_COORDSET_EXTENTS_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/bump/views/UniformCoordsetView.hpp"
#include "axom/bump/views/RectilinearCoordsetView.hpp"

#include <conduit/conduit.hpp>

namespace axom
{
namespace bump
{
namespace detail
{

/*!
 * \brief Base template for computing coordset extents. This algorithm handles
 *        any coordset view.
 */
template <typename ExecSpace, typename CoordsetView, typename DataType, int NDIMS>
struct ComputeCoordsetExtents
{
  /*!
   * \brief Compute the spatial extents of the coordset into a device array.
   *        This implementation iterates over the coordset view and finds
   *        the min/max extents from the points.
   *
   * \param coordsetView The coordset view to use.
   * \param[out] extentsView A view on the device that contains extents,
   *                         stored: xmin, xmax, ymin, ymax, zmin, zmax.
   */
  static void computeExtents(CoordsetView coordsetView, axom::ArrayView<double> extentsView)
  {
    AXOM_ANNOTATE_SCOPE("computeExtents");

    axom::for_all<ExecSpace>(
      CoordsetView::dimension(),
      AXOM_LAMBDA(axom::IndexType dim) {
        double &minValue = extentsView[2 * dim];
        double &maxValue = extentsView[2 * dim + 1];
        minValue = axom::numeric_limits<double>::max();
        maxValue = -axom::numeric_limits<double>::max();
      });

    axom::for_all<ExecSpace>(
      coordsetView.numberOfNodes(),
      AXOM_LAMBDA(axom::IndexType index) {
        const auto pt = coordsetView[index];
        for(int d = 0; d < CoordsetView::dimension(); d++)
        {
          double *minValue = extentsView.data() + 2 * d;
          double *maxValue = minValue + 1;
          const auto value = static_cast<double>(pt[d]);
          axom::atomicMin<ExecSpace>(minValue, value);
          axom::atomicMax<ExecSpace>(maxValue, value);
        }
      });
  }
};

/*!
 * Specialization for UniformCoordsetView that does less work.
 */
template <typename ExecSpace, typename DataType, int NDIMS>
struct ComputeCoordsetExtents<ExecSpace, axom::bump::views::UniformCoordsetView<DataType, NDIMS>, DataType, NDIMS>
{
  using CoordsetView = axom::bump::views::UniformCoordsetView<DataType, NDIMS>;

  /*!
   * \brief Compute the spatial extents of the coordset into a device array.
   *
   * \param coordsetView The coordset view to use.
   * \param[out] extentsView A view on the device that contains extents,
   *                         stored: xmin, xmax, ymin, ymax, zmin, zmax.
   */
  static void computeExtents(CoordsetView coordsetView, axom::ArrayView<double> extentsView)
  {
    AXOM_ANNOTATE_SCOPE("computeExtentsUniform");
    axom::for_all<ExecSpace>(
      NDIMS,
      AXOM_LAMBDA(axom::IndexType d) {
        const auto n = coordsetView.indexing().logicalDimensions()[d] - 1;
        extentsView[2 * d] = coordsetView.origin()[d];
        extentsView[2 * d + 1] = coordsetView.origin()[d] + coordsetView.spacing()[d] * n;
      });
  }
};

/*!
 * Specialization for RectilinearCoordsetView that does less work.
 */
template <typename ExecSpace, typename DataType, int NDIMS>
class ComputeCoordsetExtents<
  ExecSpace,
  typename std::conditional<NDIMS == 3,
                            axom::bump::views::RectilinearCoordsetView3<DataType>,
                            axom::bump::views::RectilinearCoordsetView2<DataType>>::type,
  DataType,
  NDIMS>
{
public:
  using CoordsetView =
    typename std::conditional<NDIMS == 3,
                              axom::bump::views::RectilinearCoordsetView3<DataType>,
                              axom::bump::views::RectilinearCoordsetView2<DataType>>::type;

  /*!
   * \brief Compute the spatial extents of the coordset into a device array.
   *
   * \param coordsetView The coordset view to use.
   * \param[out] extentsView A view on the device that contains extents,
   *                         stored: xmin, xmax, ymin, ymax, zmin, zmax.
   */
  static void computeExtents(CoordsetView coordsetView, axom::ArrayView<double> extentsView)
  {
    AXOM_ANNOTATE_SCOPE("computeExtentsRectilinear");
    axom::for_all<ExecSpace>(
      NDIMS,
      AXOM_LAMBDA(axom::IndexType d) {
        const auto coordsView = coordsetView.getCoordinates(d);
        extentsView[2 * d] = coordsView[0];
        extentsView[2 * d + 1] = coordsView[coordsView.size() - 1];
      });
  }
};

}  // end namespace detail

/**
 * \brief Compute coordset extents.
 *
 * \tparam ExecSpace The execution space where the algorithm runs.
 * \tparam CoordsetView The coordset view type.
 */
template <typename ExecSpace, typename CoordsetView>
class CoordsetExtents
{
public:
  static constexpr int NVALUES = CoordsetView::dimension() * 2;

  /*!
   * \brief Constructor
   *
   * \param coordsetView The coordset view that wraps the coordset to be examined.
   */
  CoordsetExtents(const CoordsetView &coordsetView) : m_coordsetView(coordsetView) { }

  /*!
   * \brief Compute the spatial extents of the coordset and bring results to the host.
   *
   * \param[out] extents The extents on the host.
   */
  void execute(double extents[NVALUES]) const
  {
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    axom::Array<double> deviceExtents(NVALUES, NVALUES, allocatorID);
    auto deviceExtentsView = deviceExtents.view();
    computeExtents(deviceExtentsView);
    axom::copy(extents, deviceExtentsView.data(), sizeof(double) * NVALUES);
  }

  /*!
   * \brief Compute the spatial extents of the coordset into a device array.
   *
   * \param[out] extentsView A view on the device that contains extents,
   *                         stored: xmin, xmax, ymin, ymax, zmin, zmax.
   */
  void computeExtents(axom::ArrayView<double> extentsView) const
  {
    AXOM_ANNOTATE_SCOPE("CoordsetExtents");
    axom::for_all<ExecSpace>(
      extentsView.size(),
      AXOM_LAMBDA(axom::IndexType index) { extentsView[index] = 0.; });
    // Use the appropriate specialization to compute the extents.
    using Implementation = detail::ComputeCoordsetExtents<ExecSpace,
                                                          CoordsetView,
                                                          typename CoordsetView::value_type,
                                                          CoordsetView::dimension()>;
    Implementation::computeExtents(m_coordsetView, extentsView);
  }

  CoordsetView m_coordsetView;
};

}  // end namespace bump
}  // end namespace axom

#endif
