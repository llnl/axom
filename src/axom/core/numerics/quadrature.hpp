// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

#include "axom/core/Array.hpp"
#include "axom/core/memory_management.hpp"

/*!
 * \file quadrature.hpp
 * The functions declared in this header file find the nodes and weights of 
 *   arbitrary order quadrature rules
 */

namespace axom
{
namespace numerics
{

/*!
 * \brief Enumerates the 1D quadrature families implemented in Axom.
 */
enum class QuadratureType : int
{
  Invalid = -1,
  GaussLegendre = 0,
  GaussLobatto = 1,
  OpenUniform = 2,
  ClosedUniform = 3,
  OpenHalfUniform = 4,
  ClosedGL = 5
};

/*!
 * \brief Returns true when the supplied integer corresponds to a valid
 *        `QuadratureType` enumerator.
 */
bool is_valid_quadrature_type(int quadratureType);

/*!
 * \class QuadratureRule
 *
 * \brief Stores fixed views to arrays of 1D quadrature points and weights
 */
class QuadratureRule
{
  // Define friend functions so rules can only be created via get_rule() methods
  friend QuadratureRule get_gauss_legendre(int, int);
  friend QuadratureRule get_gauss_lobatto(int, int);
  friend QuadratureRule get_open_uniform(int, int);
  friend QuadratureRule get_closed_uniform(int, int);
  friend QuadratureRule get_open_half_uniform(int, int);
  friend QuadratureRule get_closed_gl(int, int);

public:
  //! \brief Accessor for the full array of quadrature nodes
  AXOM_HOST_DEVICE
  axom::ArrayView<const double> nodes() const { return axom::ArrayView<const double>(m_nodes); }

  //! \brief Accessor for the full array of quadrature weights
  AXOM_HOST_DEVICE
  axom::ArrayView<const double> weights() const { return axom::ArrayView<const double>(m_weights); }

  //! \brief Accessor for quadrature nodes
  AXOM_HOST_DEVICE
  double node(size_t idx) const { return m_nodes[static_cast<axom::IndexType>(idx)]; };

  //! \brief Accessor for quadrature weights
  AXOM_HOST_DEVICE
  double weight(size_t idx) const { return m_weights[static_cast<axom::IndexType>(idx)]; };

  //! \brief Accessor for the size of the quadrature rule
  AXOM_HOST_DEVICE
  int getNumPoints() const { return static_cast<int>(m_nodes.size()); }

private:
  //! \brief Use a private constructor to avoid creation of an invalid rule
  QuadratureRule(axom::ArrayView<double> nodes, axom::ArrayView<double> weights)
    : m_nodes(nodes)
    , m_weights(weights) { };

  axom::ArrayView<double> m_nodes;
  axom::ArrayView<double> m_weights;
};

/*!
 * \brief Computes a 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * Algorithm adapted from the MFEM implementation in `mfem/fem/intrules.cpp`
 * 
 * \note This method constructs the points by scratch each time, without caching
 * \sa get_gauss_legendre(int)
 */
void compute_gauss_legendre_data(int npts,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of Gauss-Legendre points 
 *
 * \param [in] npts The number of points in the rule
 * 
 * A Gauss-Legendre rule with \a npts points can exactly integrate
 *  polynomials of order 2 * npts - 1
 *
 * \note If this method has already been called for a given order, it will reuse the same quadrature points
 *  without needing to recompute them
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s of stored nodes and weights
 */
QuadratureRule get_gauss_legendre(int npts, int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes a 1D quadrature rule of Gauss-Lobatto points.
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 *
 * A Gauss-Lobatto rule with \a npts points can exactly integrate
 * polynomials of order `2 * npts - 3` for `npts > 1`.
 */
void compute_gauss_lobatto_data(int npts,
                                axom::Array<double>& nodes,
                                axom::Array<double>& weights,
                                int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of
 * Gauss-Lobatto points.
 *
 * \param [in] npts The number of points in the rule
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s
 *         of stored nodes and weights
 */
QuadratureRule get_gauss_lobatto(int npts, int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Returns an Axom quadrature rule by family.
 *
 * \param [in] quadratureType The quadrature family to construct.
 * \param [in] npts The number of quadrature points in the rule.
 *
 * \note `QuadratureType::Invalid` selects Axom's default rule, which is
 *       currently Gauss-Legendre.
 */
QuadratureRule get_quadrature_rule(QuadratureType quadratureType,
                                   int npts,
                                   int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Returns the highest polynomial degree integrated exactly by a 1D
 *        quadrature family with `npts` points.
 *
 * \param [in] quadratureType The quadrature family to query.
 * \param [in] npts The number of quadrature points in the rule.
 *
 * \note `QuadratureType::Invalid` follows Axom's default rule, which is
 *       currently Gauss-Legendre.
 */
int get_exact_degree(QuadratureType quadratureType, int npts);

/*!
 * \brief Computes a 1D quadrature rule of open uniform Newton-Cotes points.
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 *
 * The points are placed at `x_i = (i + 1) / (npts + 1)` for `i = 0, ..., npts - 1`.
 * This matches MFEM's `QuadratureFunctions1D::OpenUniform`.
 *
 * The rule order matches MFEM's convention: `npts - 1 + npts % 2`.
 */
void compute_open_uniform_data(int npts,
                               axom::Array<double>& nodes,
                               axom::Array<double>& weights,
                               int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of open uniform
 * Newton-Cotes points.
 *
 * \param [in] npts The number of points in the rule
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s
 *         of stored nodes and weights
 */
QuadratureRule get_open_uniform(int npts, int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes a 1D quadrature rule of closed uniform Newton-Cotes points.
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 *
 * For `npts > 1`, the points are placed at `x_i = i / (npts - 1)` for
 * `i = 0, ..., npts - 1`. For `npts == 1`, the rule is the midpoint rule at
 * `x = 0.5` with weight `1.0`, matching MFEM's `ClosedUniform`.
 *
 * The rule order matches MFEM's convention: `npts - 1 + npts % 2`.
 */
void compute_closed_uniform_data(int npts,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of closed
 * uniform Newton-Cotes points.
 *
 * \param [in] npts The number of points in the rule
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s
 *         of stored nodes and weights
 */
QuadratureRule get_closed_uniform(int npts, int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes a 1D quadrature rule of open-half uniform Newton-Cotes
 * points.
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 *
 * The points are placed at `x_i = (2 * i + 1) / (2 * npts)` for
 * `i = 0, ..., npts - 1`, matching MFEM's `OpenHalfUniform`.
 *
 * The rule order matches MFEM's convention: `npts - 1 + npts % 2`.
 */
void compute_open_half_uniform_data(int npts,
                                    axom::Array<double>& nodes,
                                    axom::Array<double>& weights,
                                    int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of open-half
 * uniform Newton-Cotes points.
 *
 * \param [in] npts The number of points in the rule
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s
 *         of stored nodes and weights
 */
QuadratureRule get_open_half_uniform(int npts, int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes a 1D quadrature rule of closed Gauss-Legendre points.
 *
 * \param [in] npts The number of points in the rule
 * \param [out] nodes The array of 1D nodes
 * \param [out] weights The array of weights
 *
 * For `npts > 2`, the rule uses the interval endpoints together with the
 * midpoints between adjacent `(npts - 1)`-point Gauss-Legendre nodes,
 * matching MFEM's `ClosedGL`.
 *
 * The rule order matches MFEM's convention: `npts - 1 + npts % 2`.
 */
void compute_closed_gl_data(int npts,
                            axom::Array<double>& nodes,
                            axom::Array<double>& weights,
                            int allocatorID = axom::getDefaultAllocatorID());

/*!
 * \brief Computes or accesses a precomputed 1D quadrature rule of closed
 * Gauss-Legendre points.
 *
 * \param [in] npts The number of points in the rule
 *
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s
 *         of stored nodes and weights
 */
QuadratureRule get_closed_gl(int npts, int allocatorID = axom::getDefaultAllocatorID());

} /* end namespace numerics */
} /* end namespace axom */
