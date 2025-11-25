// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file KnotVector.hpp
 *
 * \brief A class to represent knot vectors for NURBS
 */

#ifndef AXOM_PRIMAL_KNOTVECTOR_HPP
#define AXOM_PRIMAL_KNOTVECTOR_HPP

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <ostream>
#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T>
class KnotVector;

/*! \brief Overloaded output operator for knot vectors */
template <typename T>
std::ostream& operator<<(std::ostream& os, const KnotVector<T>& kVec);

/*!
 * \class KnotVector
 *
 * \brief Represents the knot vector for a B-Spline/NURBS Curve/Surface
 * \tparam T the coordinate type, e.g., double, float, etc.
 * 
 * Contains methods for knot vector manipulation shared between NURBS curves and surfaces.
 * This class maintains the degree and requirements on a knot vector, namely
 *  > Monotonicity (u_i <= u_i+1)
 *  > Open/Clamped-ness (degree + 1 knots at beginning and end)
 *  > Continuity (maximum internal knot multiplicity = degree)
 */
template <typename T>
class KnotVector
{
public:
  AXOM_STATIC_ASSERT_MSG(std::is_arithmetic<T>::value,
                         "A knot vector must be defined using an arithmetic type");

public:
  ///@{
  /**
   * \name Constructors for KnotVector
   *  
   * The KnotVector class provides constructors that allow initialization from a user-supplied array of knots
   * and a specified degree and enforces several conditions to ensure validity: 
   * - the degree must be at least -1
   * - the knot array must contain at least \a (degree + 1) elements, 
   * - if the knot array is not empty, its data pointer must not be nullptr.
   * 
   * These checks guarantee that the constructed KnotVector adheres to the requirements for B-Spline/NURBS curves,
   * such as monotonicity, clamped ends, and appropriate internal knot multiplicity. The isValid() method is used
   * to verify that the resulting instance meets all necessary criteria for a valid knot span.
   */

  /*!
   * \brief Constructor from a user-supplied knot vector (axom::ArrayView<const T>)
   *
   * \param [in] knots the knot vector
   * \param [in] degree the degree of the curve
   * \pre \a degree >= 1. When degree is less than 0, the KnotVector is invalid
   * \pre The \a knots can be empty when the degree is -1, otherwise knots.data() is not \a nullptr
   * \sa isValid() tests conditions for a valid knot span instance
   */
  KnotVector(axom::ArrayView<const T> knots, int degree) : m_deg(degree)
  {
    SLIC_ASSERT(degree >= -1);
    SLIC_ASSERT(knots.size() >= (degree + 1));
    SLIC_ASSERT(knots.empty() || knots.data() != nullptr);

    if(!knots.empty())
    {
      m_knots = knots;
    }

    SLIC_ASSERT(isValid());
  }

  /*!
   * \brief Constructor from a user-supplied knot vector (axom::ArrayView<T>)
   *
   * \param [in] knots the knot vector
   * \param degree The degree of the knot vector.
   * \overload for ArrayView of non-const T
   */
  KnotVector(axom::ArrayView<T> knots, int degree)
    : KnotVector(axom::ArrayView<const T>(knots.data(), knots.size()), degree)
  { }

  /// \brief Default constructor for an empty (invalid) knot vector
  KnotVector() : m_deg(-1) { }

  /*!
   * \brief Constructor for a normalized, uniform knot vector
   *
   * \param [in] npts the number of control points
   * \param [in] degree the degree
   * 
   * \pre degree is greater than or equal to 0, and deg < npts.
   */
  KnotVector(axom::IndexType npts, int degree) { makeUniform(npts, degree); }

  /*!
   * \brief Constructor from a user-supplied knot vector (C-style array)
   * 
   * \param [in] knots the knot vector
   * \param [in] nkts the length of the knot vector
   * \param [in] degree the degree of the curve
   * 
   * \see KnotVector(axom::ArrayView<const T> knots, int degree)
   */
  KnotVector(const T* knots, axom::IndexType nkts, int degree)
    : KnotVector(axom::ArrayView<const T>(knots, nkts), degree)
  { }

  /*!
   * \brief Constructor from a user-supplied knot vector (axom::Array)
   * 
   * \param [in] knots the knot vector
   * \param [in] degree the degree of the curve
   * 
   * \see KnotVector(axom::ArrayView<const T> knots, int degree)
   */
  KnotVector(const axom::Array<T>& knots, int degree) : KnotVector(knots.view(), degree) { }

  ///@}

  ///@{
  /// \name Query/modify knot vector properties and sizes

  /*!
   * \brief Give the knot vector uniformly spaced internal knots
   *
   * \param [in] npts The number of (implied) control points
   * \param [in] deg The degree of the curve
   * 
   * \pre Requires that npts + deg + 1 > 0 and deg < npts
   */
  void makeUniform(axom::IndexType npts, int deg)
  {
    SLIC_ASSERT(npts + deg + 1 >= 0);
    SLIC_ASSERT(deg < npts);

    // Set the degree
    m_deg = deg;

    // Initialize the knot vector
    m_knots.resize(npts + deg + 1);

    // Knots for clamped-ness
    for(int i = 0; i < deg + 1; ++i)
    {
      m_knots[i] = 0.0;
      m_knots[npts + deg - i] = 1.0;
    }

    // Interior knots (if any)
    for(int i = 0; i < npts - deg - 1; ++i)
    {
      m_knots[deg + 1 + i] = (i + 1.0) / static_cast<T>(npts - deg);
    }
  }

  /// \brief Return the degree of the knot vector
  int getDegree() const { return m_deg; }

  /*!
   * \brief Reset the knot vector for the given degree
   *
   * \param [in] degree The target degree
   * 
   * If the target degree is greater than the current degree, the knot vector
   *  must be expanded to remain valid (clamped)
   * 
   * \warning This method will always replace existing knot values
   */
  void setDegree(int degree)
  {
    SLIC_ASSERT(degree >= 0);

    if(degree > m_deg)
    {
      makeUniform(degree + 1, degree);
    }
    else
    {
      makeUniform(getNumControlPoints(), degree);
    }
  }

  /// \brief Return the number of knots in the knot vector
  axom::IndexType getNumKnots() const { return m_knots.size(); }

  /// \brief Return the number of control points implied by the knot vector
  axom::IndexType getNumControlPoints() const { return m_knots.size() - m_deg - 1; }

  /// \brief Return the value of the smallest knot, or 0 if there are no knots
  T getMinKnot() const { return (m_knots.size() > 0) ? m_knots[0] : T {}; }
  /// \brief Return the value of the largest knot, or 0 if there are no knots
  T getMaxKnot() const { return (m_knots.size() > 0) ? m_knots[m_knots.size() - 1] : T {}; }

  /// \brief Clear the list of knots
  void clear()
  {
    m_deg = -1;
    m_knots.clear();
  }

  /// \brief Return if the knot vector is valid
  bool isValid() const
  {
    // Check degree
    if(m_deg < 0)
    {
      return false;
    }

    // Check for monotonicity
    for(int i = 0; i < m_knots.size() - 1; ++i)
    {
      if(m_knots[i] > m_knots[i + 1])
      {
        return false;
      }
    }

    // Check for clamped-ness
    const auto nkts = m_knots.size();
    const auto minKnot = getMinKnot();
    const auto maxKnot = getMaxKnot();
    for(int i = 0; i < m_deg + 1; ++i)
    {
      if(m_knots[i] != minKnot)
      {
        return false;
      }
      if(m_knots[nkts - 1 - i] != maxKnot)
      {
        return false;
      }
    }

    // Check for continuity
    T this_knot = m_knots[m_deg];
    int this_multiplicity = 1;
    for(int i = m_deg + 1; i < nkts - m_deg - 1; ++i)
    {
      if(m_knots[i] != this_knot)
      {
        this_knot = m_knots[i];
        this_multiplicity = 1;
      }
      else
      {
        this_multiplicity++;
        if(this_multiplicity > m_deg)
        {
          return false;
        }
      }
    }

    return true;
  }

  /*!
   * \brief Check if a knot span is valid 
   *
   * \param [in] span The span to check
   * 
   * A valid span is one that is within the knot vector and has a non-zero length
   * 
   * \return True if the span is valid, false otherwise
   */
  bool isValidSpan(axom::IndexType span) const
  {
    return span >= m_deg && span < m_knots.size() - m_deg - 1 && m_knots[span] != m_knots[span + 1];
  }

  /*!
   * \brief Check if a knot span is valid and contains the parameter value t
   *
   * \param [in] span The span to check
   * \param [in] t The parameter value
   * 
   * A valid span is one that is within the knot vector, has a non-zero length,
   * and contains the parameter value t
   * 
   * \return True if the span is valid, false otherwise
   */
  bool isValidSpan(axom::IndexType span, T t) const
  {
    return isValidSpan(span) && t >= m_knots[span] && t <= m_knots[span + 1];
  }

  ///@}

  ///@{
  /// \name Query/modify/access knots

  /// \brief Accessor for the array of knots
  axom::Array<T>& getArray() { return m_knots; }

  /// \brief Const accessor for the array of knots
  const axom::Array<T>& getArray() const { return m_knots; }

  /*!
   * \brief Getter for the knot vector
   * 
   * \param [in] i The vector index
   * 
   * \return A const reference to the knot value at index i
   */
  const T& operator[](axom::IndexType i) const { return m_knots[i]; }

  /*!
   * \brief Setter for the knot vector
   *
   * \param [in] i The vector index
   * 
   * \pre Assumes that the knot vector will remain valid after the change
   * 
   * \return A reference to the knot value at index i
   */
  T& operator[](axom::IndexType i)
  {
    SLIC_ASSERT(isValid());
    return m_knots[i];
  }

  /// \brief Return an array of the unique knot values
  axom::Array<T> getUniqueKnots() const
  {
    axom::Array<T> unique_knots;
    for(int i = 0; i < m_knots.size(); ++i)
    {
      if(i == 0 || m_knots[i] != m_knots[i - 1])
      {
        unique_knots.push_back(m_knots[i]);
      }
    }

    return unique_knots;
  }

  /// \brief Return the number of valid knot spans
  axom::IndexType getNumKnotSpans() const
  {
    axom::IndexType num_spans = 0;
    for(int i = 0; i < m_knots.size() - 1; ++i)
    {
      if(m_knots[i] != m_knots[i + 1])
      {
        num_spans++;
      }
    }

    return num_spans;
  }

  ///@}

  ///@{
  /// \name Query/modify KnotVector parametrization

  /// \brief Normalize the knot vector to the span of [0, 1]
  void normalize()
  {
    const T min_knot = getMinKnot();
    if(const T span = getMaxKnot() - min_knot; span > 0)
    {
      for(T& knot : m_knots)
      {
        knot = (knot - min_knot) / span;
      }
    }
    else
    {
      for(T& knot : m_knots)
      {
        knot = T {};
      }
    }
  }

  /*!
   * \brief Rescale the knot vector to the span of [a, b]
   * 
   * \param [in] a The lower bound of the new knot vector
   * \param [in] b The upper bound of the new knot vector
   * 
   * \pre Requires a < b
   */
  void rescale(T a, T b)
  {
    SLIC_ASSERT(a < b);

    const T min_knot = getMinKnot();
    if(const T span = getMaxKnot() - min_knot; span > 0)
    {
      for(T& knot : m_knots)
      {
        knot = axom::utilities::lerp(a, b, (knot - min_knot) / span);
      }
    }
    else
    {
      for(T& knot : m_knots)
      {
        knot = a;
      }
    }
  }

  /// \brief Reverse the knot vector
  void reverse()
  {
    const axom::IndexType nkts = m_knots.size();
    const axom::IndexType knot_mid = (nkts + 1) / 2;

    // Flip the vector
    for(int i = 0; i < knot_mid; ++i)
    {
      axom::utilities::swap(m_knots[i], m_knots[nkts - 1 - i]);
    }

    // Replace each knot with sum - knot_value
    const T the_sum = getMinKnot() + getMaxKnot();
    for(int i = 0; i < nkts; ++i)
    {
      m_knots[i] = the_sum - m_knots[i];
    }
  }

  /*!
   * \brief Insert a knot into the vector r times
   * 
   * \param [in] span The span into which the knot will be inserted
   * \param [in] t The value of the knot to insert
   * \param [in] r The number of times to insert the knot
   *  
   * \pre Assumes that the input span is correct for input t, and that
   *  the knot is not present enough times to exceed the degree 
   */
  void insertKnotBySpan(axom::IndexType span, T t, int r)
  {
    SLIC_ASSERT(isValidSpan(span, t));

    for(int i = 0; i < r; ++i)
    {
      m_knots.insert(m_knots.begin() + span + 1, t);
    }

    SLIC_ASSERT(isValid());
  }

  /*!
   * \brief Insert a knot into the vector to have the given multiplicity
   * 
   * \param [in] t The value of the knot to insert
   * \param [in] target_mutliplicity The number of times the knot will be present
   *
   * \pre Assumes that the input t is within the knot vector (up to some tolerance)
   * 
   * \note If the knot is already present, it will be inserted
   *  up to the given multiplicity, or the maximum permitted by the degree
   */
  void insertKnot(T t, int target_multiplicity)
  {
    SLIC_ASSERT(isValidParameter(t));

    int multiplicity;
    const auto span = findSpan(t, multiplicity);

    // Compute how many knots should be inserted
    const int r =
      axom::utilities::clampVal(target_multiplicity - multiplicity, 0, m_deg - multiplicity);

    insertKnotBySpan(span, t, r);
  }

  /// \brief Checks if given parameter is in knot span (to a tolerance)
  bool isValidParameter(T t, T EPS = 1e-5) const
  {
    return (t >= getMinKnot() - EPS) && (t <= getMaxKnot() + EPS);
  }

  /// \brief Checks if given parameter is *interior* to knot span (to a tolerance)
  bool isValidInteriorParameter(T t) const { return (t > getMinKnot()) && (t < getMaxKnot()); }

  ///@}

  ///@{
  /// \name Functions for evaluating knot vector

  /*!
   * \brief Returns the index of the knot span containing parameter t
   * 
   * \param [in] t The parameter value
   * 
   * For knot vector {u_0, ..., u_n}, returns i such that u_i <= t < u_i+1
   *  if t == u_n, returns i such that u_i < t <= u_i+1 (i.e. i = n - degree - 1)
   * 
   * Implementation adapted from Algorithm A2.1 on page 68 of "The NURBS Book"
   * 
   * \pre Assumes that the input t is in the span [u_0, u_n] (up to some tolerance)
   * 
   * \note If t is outside the knot span up to this tolerance, it is clamped to the span
   * 
   * \return The index of the knot span containing t
   */
  axom::IndexType findSpan(T t) const
  {
    SLIC_ASSERT(isValidParameter(t));

    const axom::IndexType nkts = m_knots.size();

    // Handle cases where t is outside the knot span within a tolerance
    //  by implicitly clamping it to the nearest span
    if(t <= getMinKnot())
    {
      return m_deg;
    }

    if(t >= getMaxKnot())
    {
      return nkts - m_deg - 2;
    }

    // perform binary search on the knots,
    axom::IndexType low = m_deg;
    axom::IndexType high = nkts - m_deg - 1;
    axom::IndexType mid = (low + high) / 2;
    while(t < m_knots[mid] || t >= m_knots[mid + 1])
    {
      (t < m_knots[mid]) ? high = mid : low = mid;
      mid = (low + high) / 2;
    }

    return mid;
  }

  /*!
   * \brief Returns the index of the knot span containing parameter t and its multiplicity
   * 
   * \param [in] t The parameter value
   * \param [out] multiplicity The multiplicity of the knot at t
   *
   * For knot vector {u_0, ..., u_n}, returns i such that u_i <= t < u_i+1
   *  if t == u_n, returns i such that u_i < t <= u_i+1 (i.e. i = n - degree - 1)
   * 
   * \pre Assumes that the input t is within the knot vector (up to some tolerance)
   * 
   * \note If t is outside the knot span up to this tolerance, the returned multiplicity
   *  will be equal to the degree + 1 (required for clamped curves)
   * 
   * \return The index of the knot span containing t
   */
  axom::IndexType findSpan(T t, int& multiplicity) const
  {
    SLIC_ASSERT(isValidParameter(t));

    const auto span = findSpan(t);

    // Early exit for known multiplicities
    if(t <= getMinKnot() || t >= getMaxKnot())
    {
      multiplicity = m_deg + 1;
      return span;
    }

    multiplicity = 0;
    for(auto i = span; i >= 0; --i)
    {
      if(m_knots[i] == t)
      {
        multiplicity++;
      }
      else
      {
        break;
      }
    }

    return span;
  }

  /*!
   * \brief Evaluates the NURBS basis functions for span at parameter value t
   * 
   * \param [in] span The span in which to evaluate the basis functions
   * \param [in] t The parameter value
   * 
   * \pre Assumes that the input t is within the correct span
   * Implementation adapted from Algorithm A2.2 on page 70 of "The NURBS Book".
   * 
   * \return An array of the `m_deg + 1` non-zero basis functions evaluated at t
   */
  axom::Array<T> calculateBasisFunctionsBySpan(axom::IndexType span, T t) const
  {
    SLIC_ASSERT(isValidSpan(span, t));

    axom::Array<T> N(m_deg + 1);
    axom::Array<T> left(m_deg + 1);
    axom::Array<T> right(m_deg + 1);

    // Note: This implementation avoids division by zero and redundant computation
    // that might arise from a direct implementation of the recurrence relation
    // for basis functions. See "The NURBS Book" for details.
    N[0] = 1.0;
    for(int j = 1; j <= m_deg; ++j)
    {
      left[j] = t - m_knots[span + 1 - j];
      right[j] = m_knots[span + j] - t;
      T saved = 0.0;
      for(int r = 0; r < j; ++r)
      {
        const T tmp = N[r] / (right[r + 1] + left[j - r]);
        N[r] = saved + right[r + 1] * tmp;
        saved = left[j - r] * tmp;
      }
      N[j] = saved;
    }

    return N;
  }

  /*!
   * \brief Evaluates the NURBS basis functions for parameter value t
   * 
   * \param [in] t The parameter value
   * 
   * \pre Assumes that the input t is within the knot vector (up to a tolerance)
   * 
   * \return An array of the `m_deg + 1` non-zero basis functions evaluated at t
   */
  axom::Array<T> calculateBasisFunctions(T t) const
  {
    SLIC_ASSERT(isValidParameter(t));
    return calculateBasisFunctionsBySpan(findSpan(t), t);
  }

  /*!
   * \brief Evaluates the NURBS basis functions and derivatives for span at parameter value t
   * 
   * \param [in] span The span in which to evaluate the basis functions
   * \param [in] t The parameter value
   * \param [in] n The number of derivatives to compute
   * 
   * Implementation adapted from Algorithm A2.2 on page 70 of "The NURBS Book".
   *
   * \pre Assumes that the input t is within the provided knot span
   * 
   * \return An array of the `n + 1` derivatives evaluated at t
   */
  axom::Array<axom::Array<T>> derivativeBasisFunctionsBySpan(axom::IndexType span, T t, int n) const
  {
    SLIC_ASSERT(isValidSpan(span, t));

    axom::Array<axom::Array<T>> ders(n + 1);

    axom::Array<axom::Array<T>> ndu(m_deg + 1), a(2);
    axom::Array<T> left(m_deg + 1), right(m_deg + 1);
    for(int j = 0; j <= m_deg; j++)
    {
      ndu[j].resize(m_deg + 1);
    }
    for(int j = 0; j <= n; j++)
    {
      ders[j].resize(m_deg + 1);
    }
    a[0].resize(m_deg + 1);
    a[1].resize(m_deg + 1);

    ndu[0][0] = 1.;
    for(int j = 1; j <= m_deg; j++)
    {
      left[j] = t - m_knots[span + 1 - j];
      right[j] = m_knots[span + j] - t;
      T saved = 0.0;
      for(int r = 0; r < j; r++)
      {
        // lower triangle
        ndu[j][r] = right[r + 1] + left[j - r];
        T temp = ndu[r][j - 1] / ndu[j][r];
        // upper triangle
        ndu[r][j] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    // Load basis functions
    for(int j = 0; j <= m_deg; j++)
    {
      ders[0][j] = ndu[j][m_deg];
    }

    // This section computes the derivatives (Eq. [2.9])

    // Loop over function index.
    for(int r = 0; r <= m_deg; r++)
    {
      int s1 = 0, s2 = 1;  // Alternate rows in array a
      a[0][0] = 1.;
      // Loop to compute kth derivative
      for(int k = 1; k <= n; k++)
      {
        T d = 0.;
        int rk = r - k;
        int pk = m_deg - k;
        if(r >= k)
        {
          a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
          d = a[s2][0] * ndu[rk][pk];
        }
        int j1 = (rk >= -1) ? 1 : -rk;
        int j2 = (r - 1 <= pk) ? (k - 1) : (m_deg - r);
        for(int j = j1; j <= j2; j++)
        {
          a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
          d += a[s2][j] * ndu[rk + j][pk];
        }
        if(r <= pk)
        {
          a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
          d += a[s2][k] * ndu[r][pk];
        }
        ders[k][r] = d;
        // Switch rows
        std::swap(s1, s2);
      }
    }

    // Multiply through by the correct factors (Eq. [2.9])
    T r = static_cast<T>(m_deg);
    for(int k = 1; k <= n; k++)
    {
      for(int j = 0; j <= m_deg; j++)
      {
        ders[k][j] *= r;
      }
      r *= static_cast<T>(m_deg - k);
    }

    return ders;
  }

  /*!
   * \brief Evaluates the NURBS basis functions and derivatives for parameter value t
   * 
   * \param [in] t The parameter value
   * \param [in] n The number of derivatives to compute
   * 
   * \pre Assumes that the input t is within the knot vector (up to a tolerance)
   * 
   * \return An array of the `n + 1` derivatives evaluated at t
   */
  axom::Array<axom::Array<T>> derivativeBasisFunctions(T t, int n) const
  {
    SLIC_ASSERT(isValidParameter(t));
    return derivativeBasisFunctionsBySpan(findSpan(t), t, n);
  }

  ///@}

  ///@{
  /// \name Functions dealing with knot vector subdivision

  /*!
   * \brief Split a knot vector at the value at index `span`
   * 
   * \param [in] span The span at which the knot will be split
   * \param [out] k1 The first knot vector
   * \param [out] k2 The second knot vector
   * \param [in] normalize Whether to normalize the output knot vectors
   * 
   * \warning Assumes that the multiplicity of the knot at which 
   *  the vector will be split is equal to the degree, or the
   *  returned knot vectors will be invalid
   */
  void splitBySpan(axom::IndexType span, KnotVector& k1, KnotVector& k2, bool normalize = false) const
  {
    SLIC_ASSERT(isValidSpan(span));

    const auto nkts = getNumKnots();

    // Create a copy of the vector in case k1 or k2 is the same as *this
    KnotVector<T> data(*this);

    // Copy the degree
    k1.m_deg = m_deg;
    k2.m_deg = m_deg;

    // Copy the knots
    k1.m_knots.resize(span + 2);
    k2.m_knots.resize(nkts - span + m_deg);

    k1.m_knots[span + 1] = data[span];
    for(int i = 0; i < span + 1; ++i)
    {
      k1.m_knots[i] = data[i];
    }

    k2.m_knots[0] = data[span];
    for(int i = 0; i < nkts - span + m_deg - 1; ++i)
    {
      k2.m_knots[i + 1] = data[span + i - m_deg + 1];
    }

    if(normalize)
    {
      k1.normalize();
      k2.normalize();
    }

    SLIC_ASSERT(k1.isValid());
    SLIC_ASSERT(k2.isValid());
  }

  /*!
   * \brief Split a knot vector at the value t
   * 
   * \param [in] t The value at which the knot will be split
   * \param [out] k1 The first knot vector
   * \param [out] k2 The second knot vector
   * \param [in] normalize Whether to normalize the output knot vectors
   * 
   * \pre Assumes that the input t is *interior* to the knot vector
   */
  void split(T t, KnotVector& k1, KnotVector& k2, bool normalize = false) const
  {
    SLIC_ASSERT(isValidInteriorParameter(t));

    int multiplicity;
    axom::IndexType span = findSpan(t, multiplicity);

    k1 = *this;
    k1.insertKnotBySpan(span, t, m_deg - multiplicity);

    k1.splitBySpan(span + m_deg - multiplicity, k1, k2, normalize);
  }

  ///@}

  /*!
   * \brief Check equality of two knot vectors
   * 
   * \param [in] lhs The first knot vector
   * \param [in] rhs The second knot vector
   * 
   * \return True if the knot vectors are equal, false otherwise 
   */
  friend inline bool operator==(const KnotVector<T>& lhs, const KnotVector<T>& rhs)
  {
    return lhs.m_deg == rhs.m_deg && lhs.m_knots == rhs.m_knots;
  }

  /*!
   * \brief Check inequality of two knot vectors
   * 
   * \param [in] lhs The first knot vector
   * \param [in] rhs The second knot vector
   * 
   * \return True if the knot vectors are not equal, false otherwise 
   */
  friend inline bool operator!=(const KnotVector<T>& lhs, const KnotVector<T>& rhs)
  {
    return !(lhs == rhs);
  }

  /*!
   * \brief Simple formatted print of a knot vector instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    int nkts = m_knots.size();

    os << "{ degree " << m_deg << " knot vector ";
    for(int i = 0; i < nkts; ++i)
    {
      os << m_knots[i] << (i < nkts - 1 ? ", " : "");
    }
    os << "}";

    return os;
  }

private:
  int m_deg;
  axom::Array<T> m_knots;
};

//------------------------------------------------------------------------------
/// Free functions related to KnotVector
//------------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const KnotVector<T>& kvector)
{
  kvector.print(os);
  return os;
}
}  // namespace primal
}  // namespace axom

/// Overload to format a primal::KnotVector using fmt
template <typename T>
struct axom::fmt::formatter<axom::primal::KnotVector<T>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_KNOTVECTOR_HPP
