// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/NumericLimits.hpp"
#include "axom/core/numerics/quadrature.hpp"

// For math constants and includes
#include "axom/config.hpp"

#include <cmath>
#include <cassert>
#include <cstdint>
#include <map>
#include <mutex>
#include <string>

namespace axom
{
namespace numerics
{

bool is_valid_quadrature_type(int quadratureType)
{
  switch(static_cast<QuadratureType>(quadratureType))
  {
  case QuadratureType::Invalid:
  case QuadratureType::GaussLegendre:
  case QuadratureType::GaussLobatto:
  case QuadratureType::OpenUniform:
  case QuadratureType::ClosedUniform:
  case QuadratureType::OpenHalfUniform:
  case QuadratureType::ClosedGL:
    return true;
  default:
    return false;
  }
}

const std::map<std::string, QuadratureType>& stringToQuadratureType()
{
  static const std::map<std::string, QuadratureType> quadrature_types {
    {"default", QuadratureType::Invalid},
    {"invalid", QuadratureType::Invalid},
    {"gausslegendre", QuadratureType::GaussLegendre},
    {"gausslobatto", QuadratureType::GaussLobatto},
    {"openuniform", QuadratureType::OpenUniform},
    {"closeduniform", QuadratureType::ClosedUniform},
    {"openhalfuniform", QuadratureType::OpenHalfUniform},
    {"closedgl", QuadratureType::ClosedGL}};

  return quadrature_types;
}

void compute_gauss_legendre_data(int npts,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID);
void compute_gauss_lobatto_data(int npts,
                                axom::Array<double>& nodes,
                                axom::Array<double>& weights,
                                int allocatorID);
void compute_open_half_uniform_data(int npts,
                                    axom::Array<double>& nodes,
                                    axom::Array<double>& weights,
                                    int allocatorID);
void compute_closed_gl_data(int npts,
                            axom::Array<double>& nodes,
                            axom::Array<double>& weights,
                            int allocatorID);

namespace
{
struct RuleStorage
{
  axom::Array<double> nodes;
  axom::Array<double> weights;
};

std::uint64_t make_rule_key(int npts, int allocatorID)
{
  const auto n = static_cast<std::uint64_t>(static_cast<std::uint32_t>(npts));
  const auto a = static_cast<std::uint64_t>(static_cast<std::uint32_t>(allocatorID));
  return (a << 32) | n;
}

template <typename ComputeRuleData>
RuleStorage& get_cached_rule_storage(int npts,
                                     int allocatorID,
                                     axom::FlatMap<std::uint64_t, RuleStorage>& rule_library,
                                     std::mutex& rule_library_mutex,
                                     ComputeRuleData&& computeRuleData)
{
  const std::lock_guard<std::mutex> lock(rule_library_mutex);
  const std::uint64_t key = make_rule_key(npts, allocatorID);

  auto [it, inserted] = rule_library.try_emplace(key);
  if(inserted)
  {
    computeRuleData(npts, it->second.nodes, it->second.weights, allocatorID);
  }

  return it->second;
}

/*!
 * \brief Computes quadrature weights for an interpolatory rule on `[0, 1]`
 *        from its nodes.
 *
 * \param [in] nodes The interpolation nodes that define the Lagrange basis
 * \param [out] weights The quadrature weights corresponding to `nodes`
 * \param [in] allocatorID The allocator used for temporary storage and the
 *             output `weights` array
 *
 * The returned weights are `w_j = integral_0^1 L_j(x) dx`, where `L_j` is the
 * `j`th Lagrange basis polynomial. This makes the rule exact for polynomials
 * up to degree `npts - 1`, and the weights sum to one because the Lagrange
 * basis forms a partition of unity. Weight positivity is not guaranteed for
 * arbitrary nodes.
 *
 * Several quadrature families define only their nodes directly; once the nodes
 * are known, their interpolatory weights can all be computed this way. The
 * temporary Gauss-Legendre rule is used only as an exact integration rule for
 * the degree `npts - 1` basis polynomials.
 */
void compute_interpolatory_weights(const axom::Array<double>& nodes,
                                   axom::Array<double>& weights,
                                   int allocatorID)
{
  const int npts = nodes.size();
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  weights = axom::Array<double>(npts, npts, allocatorID);

  if(npts == 1)
  {
    weights[0] = 1.0;
    return;
  }

  if(npts == 2)
  {
    weights[0] = 0.5;
    weights[1] = 0.5;
    return;
  }

  axom::Array<double> glNodes;
  axom::Array<double> glWeights;
  compute_gauss_legendre_data((npts + 1) / 2, glNodes, glWeights, allocatorID);

  for(int j = 0; j < npts; ++j)
  {
    double wj = 0.0;
    for(int q = 0; q < glNodes.size(); ++q)
    {
      const double x = glNodes[q];
      double basis = 1.0;
      for(int k = 0; k < npts; ++k)
      {
        if(k == j)
        {
          continue;
        }

        basis *= (x - nodes[k]) / (nodes[j] - nodes[k]);
      }
      wj += glWeights[q] * basis;
    }
    weights[j] = wj;
  }
}
}  // namespace

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
                                 int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  weights = axom::Array<double>(npts, npts, allocatorID);

  if(npts == 1)
  {
    nodes[0] = 0.5;
    weights[0] = 1.0;
    return;
  }
  if(npts == 2)
  {
    nodes[0] = 0.21132486540518711775;
    nodes[1] = 0.78867513459481288225;

    weights[0] = weights[1] = 0.5;
    return;
  }
  if(npts == 3)
  {
    nodes[0] = 0.11270166537925831148207345;
    nodes[1] = 0.5;
    nodes[2] = 0.88729833462074168851792655;

    weights[0] = 0.2777777777777777777777778;
    weights[1] = 0.4444444444444444444444444;
    weights[2] = 0.2777777777777777777777778;
    return;
  }

  const int n = npts;
  const int m = (npts + 1) / 2;

  // Nodes are mirrored across x = 0.5, so only need to evaluate half
  for(int i = 1; i <= m; ++i)
  {
    // Each node is the root of a Legendre polynomial,
    //  which are approximately uniformly distributed in arccos(xi).
    // This makes cos a good initial guess for subsequent Newton iterations
    double z = std::cos(M_PI * (i - 0.25) / (n + 0.5));
    double Pp_n, P_n, dz, xi = 0.0;

    bool done = false;
    while(true)
    {
      // Recursively evaluate P_n(z) through the recurrence relation
      //  n * P_n(z) = (2n - 1) * P_{n-1}(z) - (n - 1) * P_{n - 2}(z)
      double P_nm1 = 1.0;  // P_0(z) = 1
      P_n = z;             // P_1(z) = z
      for(int j = 2; j <= n; ++j)
      {
        double P_nm2 = P_nm1;
        P_nm1 = P_n;
        P_n = ((2 * j - 1) * z * P_nm1 - (j - 1) * P_nm2) / j;
      }

      // Evaluate P'_n(z) using another recurrence relation
      //  (z^2 - 1) * P'_n(z) = n * z * P_n(z) - n * P_{n-1}(z)
      Pp_n = n * (z * P_n - P_nm1) / (z * z - 1);

      if(done)
      {
        break;
      }

      // Compute the Newton method step size
      dz = P_n / Pp_n;

      if(std::fabs(dz) < axom::numeric_limits<double>::epsilon())
      {
        done = true;

        // Scale back to [0, 1]
        xi = ((1 - z) + dz) / 2;
      }

      // Take the Newton step, repeat the process
      z -= dz;
    }

    nodes[i - 1] = xi;
    nodes[n - i] = 1.0 - xi;

    // For z \in [-1, 1], w_i = 2 / (1 - z^2) / P'_n(z)^2.
    // For nodes[i] = xi = (1 - z)/2 \in [0, 1], weights[i] = 0.5 * w_i
    weights[i - 1] = weights[n - i] = 1.0 / (4.0 * xi * (1.0 - xi) * Pp_n * Pp_n);
  }
}

void compute_open_uniform_data(int npts,
                               axom::Array<double>& nodes,
                               axom::Array<double>& weights,
                               int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  for(int i = 0; i < npts; ++i)
  {
    nodes[i] = static_cast<double>(i + 1) / static_cast<double>(npts + 1);
  }

  compute_interpolatory_weights(nodes, weights, allocatorID);
}

void compute_gauss_lobatto_data(int npts,
                                axom::Array<double>& nodes,
                                axom::Array<double>& weights,
                                int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  weights = axom::Array<double>(npts, npts, allocatorID);

  if(npts == 1)
  {
    nodes[0] = 0.5;
    weights[0] = 1.0;
    return;
  }

  nodes[0] = 0.0;
  nodes[npts - 1] = 1.0;
  weights[0] = weights[npts - 1] = 1.0 / (npts * (npts - 1.0));

  constexpr int MaxIterations = 16;
  const double tol = axom::numeric_limits<double>::epsilon();

  for(int i = 1; i <= (npts - 1) / 2; ++i)
  {
    double x = std::sin(M_PI * (static_cast<double>(i) / (npts - 1) - 0.5));
    double pNm2 = 1.0;
    double pNm1 = x;

    for(int iter = 0; iter < MaxIterations; ++iter)
    {
      pNm2 = 1.0;
      pNm1 = x;
      for(int l = 1; l < (npts - 1); ++l)
      {
        const double p = ((2 * l + 1) * x * pNm1 - l * pNm2) / (l + 1);
        pNm2 = pNm1;
        pNm1 = p;
      }

      const double dx = (x * pNm1 - pNm2) / (npts * pNm1);
      x -= dx;

      if(std::fabs(dx) <= tol * (1.0 + std::fabs(x)))
      {
        break;
      }

      assert("Gauss-Lobatto Newton iteration did not converge." && iter + 1 < MaxIterations);
    }

    pNm2 = 1.0;
    pNm1 = x;
    for(int l = 1; l < (npts - 1); ++l)
    {
      const double p = ((2 * l + 1) * x * pNm1 - l * pNm2) / (l + 1);
      pNm2 = pNm1;
      pNm1 = p;
    }

    const double node = 0.5 * (1.0 + x);
    const double weight = 1.0 / (npts * (npts - 1.0) * pNm1 * pNm1);

    nodes[i] = node;
    nodes[npts - 1 - i] = 1.0 - node;
    weights[i] = weight;
    weights[npts - 1 - i] = weight;
  }
}

void compute_closed_uniform_data(int npts,
                                 axom::Array<double>& nodes,
                                 axom::Array<double>& weights,
                                 int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  if(npts == 1)
  {
    nodes[0] = 0.5;
    weights = axom::Array<double>(1, 1, allocatorID);
    weights[0] = 1.0;
    return;
  }

  for(int i = 0; i < npts; ++i)
  {
    nodes[i] = static_cast<double>(i) / static_cast<double>(npts - 1);
  }

  compute_interpolatory_weights(nodes, weights, allocatorID);
}

void compute_open_half_uniform_data(int npts,
                                    axom::Array<double>& nodes,
                                    axom::Array<double>& weights,
                                    int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  for(int i = 0; i < npts; ++i)
  {
    nodes[i] = static_cast<double>(2 * i + 1) / static_cast<double>(2 * npts);
  }

  compute_interpolatory_weights(nodes, weights, allocatorID);
}

void compute_closed_gl_data(int npts,
                            axom::Array<double>& nodes,
                            axom::Array<double>& weights,
                            int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  nodes = axom::Array<double>(npts, npts, allocatorID);
  if(npts == 1)
  {
    nodes[0] = 0.5;
    weights = axom::Array<double>(1, 1, allocatorID);
    weights[0] = 1.0;
    return;
  }

  nodes[0] = 0.0;
  nodes[npts - 1] = 1.0;

  if(npts > 2)
  {
    axom::Array<double> glNodes;
    axom::Array<double> glWeights;
    compute_gauss_legendre_data(npts - 1, glNodes, glWeights, allocatorID);

    for(int i = 1; i < npts - 1; ++i)
    {
      nodes[i] = 0.5 * (glNodes[i - 1] + glNodes[i]);
    }
  }

  compute_interpolatory_weights(nodes, weights, allocatorID);
}

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
 * \note This implementation uses a process-wide cache protected by a mutex for thread safety.
 * 
 * \return The `QuadratureRule` object which contains axom::ArrayView<double>'s of stored nodes and weights
 */
QuadratureRule get_gauss_legendre(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  // Store cached rules keyed by (npts, allocatorID).
  static axom::FlatMap<std::uint64_t, RuleStorage> rule_library(64);
  static std::mutex rule_library_mutex;
  auto& storage = get_cached_rule_storage(npts,
                                          allocatorID,
                                          rule_library,
                                          rule_library_mutex,
                                          compute_gauss_legendre_data);
  return QuadratureRule {storage.nodes.view(), storage.weights.view()};
}

QuadratureRule get_gauss_lobatto(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  static axom::FlatMap<std::uint64_t, RuleStorage> rule_library(64);
  static std::mutex rule_library_mutex;
  auto& storage = get_cached_rule_storage(npts,
                                          allocatorID,
                                          rule_library,
                                          rule_library_mutex,
                                          compute_gauss_lobatto_data);
  return QuadratureRule {storage.nodes.view(), storage.weights.view()};
}

QuadratureRule get_quadrature_rule(QuadratureType quadratureType, int npts, int allocatorID)
{
  assert("Invalid Axom quadrature type." &&
         is_valid_quadrature_type(static_cast<int>(quadratureType)));

  switch(quadratureType)
  {
  case QuadratureType::Invalid:
  case QuadratureType::GaussLegendre:
    return get_gauss_legendre(npts, allocatorID);
  case QuadratureType::GaussLobatto:
    return get_gauss_lobatto(npts, allocatorID);
  case QuadratureType::OpenUniform:
    return get_open_uniform(npts, allocatorID);
  case QuadratureType::ClosedUniform:
    return get_closed_uniform(npts, allocatorID);
  case QuadratureType::OpenHalfUniform:
    return get_open_half_uniform(npts, allocatorID);
  case QuadratureType::ClosedGL:
    return get_closed_gl(npts, allocatorID);
  }

  assert("Unhandled Axom quadrature type." && false);
  return get_gauss_legendre(npts, allocatorID);
}

int get_exact_degree(QuadratureType quadratureType, int npts)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));
  assert("Invalid Axom quadrature type." &&
         is_valid_quadrature_type(static_cast<int>(quadratureType)));

  switch(quadratureType)
  {
  case QuadratureType::Invalid:
  case QuadratureType::GaussLegendre:
    return 2 * npts - 1;
  case QuadratureType::GaussLobatto:
    return npts == 1 ? 1 : 2 * npts - 3;
  case QuadratureType::OpenUniform:
  case QuadratureType::ClosedUniform:
  case QuadratureType::OpenHalfUniform:
  case QuadratureType::ClosedGL:
    return npts - 1 + npts % 2;
  }

  assert("Unhandled Axom quadrature type." && false);
  return 2 * npts - 1;
}

QuadratureRule get_open_uniform(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  static axom::FlatMap<std::uint64_t, RuleStorage> rule_library(64);
  static std::mutex rule_library_mutex;
  auto& storage = get_cached_rule_storage(npts,
                                          allocatorID,
                                          rule_library,
                                          rule_library_mutex,
                                          compute_open_uniform_data);
  return QuadratureRule {storage.nodes.view(), storage.weights.view()};
}

QuadratureRule get_closed_uniform(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  static axom::FlatMap<std::uint64_t, RuleStorage> rule_library(64);
  static std::mutex rule_library_mutex;
  auto& storage = get_cached_rule_storage(npts,
                                          allocatorID,
                                          rule_library,
                                          rule_library_mutex,
                                          compute_closed_uniform_data);
  return QuadratureRule {storage.nodes.view(), storage.weights.view()};
}

QuadratureRule get_open_half_uniform(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  static axom::FlatMap<std::uint64_t, RuleStorage> rule_library(64);
  static std::mutex rule_library_mutex;
  auto& storage = get_cached_rule_storage(npts,
                                          allocatorID,
                                          rule_library,
                                          rule_library_mutex,
                                          compute_open_half_uniform_data);
  return QuadratureRule {storage.nodes.view(), storage.weights.view()};
}

QuadratureRule get_closed_gl(int npts, int allocatorID)
{
  assert("Quadrature rules must have >= 1 point" && (npts >= 1));

  static axom::FlatMap<std::uint64_t, RuleStorage> rule_library(64);
  static std::mutex rule_library_mutex;
  auto& storage =
    get_cached_rule_storage(npts, allocatorID, rule_library, rule_library_mutex, compute_closed_gl_data);
  return QuadratureRule {storage.nodes.view(), storage.weights.view()};
}

} /* end namespace numerics */
} /* end namespace axom */
