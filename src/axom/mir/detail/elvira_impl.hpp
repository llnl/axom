// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// NOTE: This file is meant to be included by ElviraAlgorithm.hpp after its
//       other includes so we do not include much here.

namespace axom
{
namespace mir
{
namespace elvira
{

// NOTE: Many of the functions below were adapted from code originally
//       written by Jeff Grandy. The functions have been templated on
//       "FloatType" instead of using double so we can generate host
//       and device code for HIP without resulting in  multiple function
//       definitions at link time.

enum class Direction
{
  VERTICAL = 0,
  HORIZONTAL = 1
};

enum class Difference
{
  BACKWARD = 0,
  CENTRAL = 1,
  FORWARD = 2
};

/// This struct stores the normal computed for a plane of stencil data.
template <typename FloatType>
struct Result2D
{
  int plane {2};
  Difference difference_used {Difference::BACKWARD};
  FloatType columns[3] {0., 0., 0.};
  FloatType normal[3][2] {{0., 0.}, {0., 0.}, {0., 0.}};
};

/*!
 * \brief Convert int value to Difference.
 * \param value The input value to convert to Difference.
 * \return The converted Difference value.
 */
inline AXOM_HOST_DEVICE Difference IntToDifference(int value)
{
  return (value == 0) ? Difference::BACKWARD
                      : ((value == 2) ? Difference::FORWARD : Difference::CENTRAL);
}

/*!
 * \brief Convert Difference value to int.
 * \param value The input value to convert to int.
 * \return The converted int value.
 */
inline AXOM_HOST_DEVICE int DifferenceToInt(Difference value)
{
  return (value == Difference::BACKWARD) ? 0 : ((value == Difference::FORWARD) ? 2 : 1);
}
/*!
 * \brief Reverse the difference direction.
 * \param value The input value to reverse.
 * \return The reversed direction.
 */
inline AXOM_HOST_DEVICE Difference Reverse(Difference value)
{
  return (value == Difference::BACKWARD)
    ? Difference::FORWARD
    : ((value == Difference::FORWARD) ? Difference::BACKWARD : Difference::CENTRAL);
}

/*!
 * \brief 2nd order 2-d mir function: returns the normal to the slope in 
 *        the direction away from the material.
 *
 * \param[out] result The result object.
 * \param vf   An array containing volume fractions for the current material
 *             arranged logically like this, according to stencil indices:
 *                       0 - lower left neighbor
 *                       1 - lower neighbor
 *                       2 - lower right neighbor
 *                       3 - left neighbor
 *                       4 - zone in question
 *                       5 - right neighbor
 *                       6 - upper left neighbor
 *                       7 - upper neighbor
 *                       8 - upper right neighbor
 *             Note that the actual vf data can be 9 or 27 elements depending
 *             on the dimensionality.
 *
 * \param ivf The indices to use to pull data out of vf.
 * \param direction hv direction to use
 *
 * \return The result object will contain these values:
 *         normal  : normal (x,y) pointing away from the material
 *         err     : the L2 error of the answer returned
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void elvira2d(Result2D<FloatType> &result,
                               const FloatType *vf,
                               const int *ivf,
                               Direction direction)
{
  const FloatType jb = vf[ivf[0]] + vf[ivf[1]] + vf[ivf[2]];  // bottom row
  const FloatType jm = vf[ivf[3]] + vf[ivf[4]] + vf[ivf[5]];  // middle row
  const FloatType jt = vf[ivf[6]] + vf[ivf[7]] + vf[ivf[8]];  // top row

  const FloatType il = vf[ivf[0]] + vf[ivf[3]] + vf[ivf[6]];  // left column
  const FloatType im = vf[ivf[1]] + vf[ivf[4]] + vf[ivf[7]];  // middle column
  const FloatType ir = vf[ivf[2]] + vf[ivf[5]] + vf[ivf[8]];  // right column

  FloatType slope[3] {0., 0., 0.};  // backward, cent, and forward differencing
  FloatType slopev = 0.;            // Sets overall sign of normal.
  int ihv;
  if(direction == Direction::VERTICAL)
  {
    slope[0] = (im - il);        // backward differencing
    slope[1] = (ir - il) / 2.0;  // central differencing
    slope[2] = (ir - im);        // forward differencing
    slopev = (jb - jt);

    result.columns[0] = il;
    result.columns[1] = im;
    result.columns[2] = ir;

    ihv = 0;
  }
  else
  {
    slope[0] = (jm - jb);        // backward differencing
    slope[1] = (jt - jb) / 2.0;  // central differencing
    slope[2] = (jt - jm);        // forward differencing
    slopev = (il - ir);

    result.columns[0] = jb;
    result.columns[1] = jm;
    result.columns[2] = jt;

    ihv = 1;
  }

  const FloatType sign = (slopev > 0.0) ? (1.) : (-1.);

  // return the normal in the direction away from the material
  // i = loop over three difference schemes.
  for(int i = 0; i < 3; i++)
  {
    result.normal[i][1 - ihv] = sign;   // steep direction
    result.normal[i][ihv] = -slope[i];  // shallow
  }
}

/*!
 * \brief Get chi squared for actual and fitted vol fracs.
 *
 * \param vf actual volume fractions.
 * \param vfs trial volume fractions.
 * \param ivf index list of volume fractions to be compared.
 * \param k The number of volume fractions to compare.
 *
 * \return chi squared (measure of error)
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE inline FloatType elvira_chisq(const FloatType *vf,
                                               const FloatType *vfs,
                                               const int *ivf,
                                               int k)
{
  FloatType chisq = 0.0;
  for(; k--;)
  {
    const FloatType c = vf[ivf[k]] - vfs[ivf[k]];
    chisq += c * c;
  }
  return chisq;
}

/*!
 * \brief Normalize 2d vector, append zero Z component.
 *
 * Notes: Assumes n2 is not (0,0) .
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE inline void norm2d(FloatType n2[2], FloatType n3[3])
{
  constexpr FloatType PTINY = 1.e-15;
  // Compute 2D vector magnitude.
  FloatType magnitude = sqrt(n2[0] * n2[0] + n2[1] * n2[1]);
  FloatType magnitude2 = magnitude * magnitude;
  FloatType magnitudeReciprocal = magnitude / (magnitude2 + PTINY);
  // normalize vector, make it 3D
  n3[0] = n2[0] * magnitudeReciprocal;
  n3[1] = n2[1] * magnitudeReciprocal;
  n3[2] = (magnitude2 <= 0.) ? 1.0 : 0.0;
}

/*!
 * \brief Compute 3d normal from two 2d normals.  
 *
 * \param A The first 2d normal.
 * \param B The second 2d normal.
 * \param[out] n3 The resulting 3D normal.
 *
 * \note Off diagonal n2 terms are +-1.  Physical direction not
 *       assigned in this function.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void norm3d(const FloatType A[2], const FloatType B[2], FloatType n3[3])
{
  constexpr FloatType PTINY = 1.e-15;

  n3[0] = A[1] * B[1];
  n3[1] = A[0] * B[0];
  n3[2] = A[1] * B[0];  // "largest" component.

  FloatType mag = sqrt(n3[0] * n3[0] + n3[1] * n3[1] + n3[2] * n3[2]);
  FloatType mag2 = mag * mag;
  // rmag = reciprocal of magnitude
  FloatType rmag = mag / (mag2 + PTINY);

  // check sign and make unit normal
  if((A[1] + B[0]) < 0)
  {
    rmag = -rmag;
  }

  n3[0] *= rmag;
  n3[1] *= rmag;
  n3[2] *= rmag;

  if(mag2 <= 0)
  {
    n3[2] = 1.0;  // test for sum==0.
  }
}

/*!
 * \brief Return sorted, positive components of vector length 3.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE inline void n_sort(const FloatType n[3],  // input vector
                                    FloatType na[3])       // output vector
{
  na[0] = axom::utilities::abs(n[0]);
  na[1] = axom::utilities::abs(n[1]);
  na[2] = axom::utilities::abs(n[2]);

  /* Sort the components. */
  if(na[0] > na[2])
  {
    axom::utilities::swap(na[0], na[2]);
  }
  if(na[1] > na[2])
  {
    axom::utilities::swap(na[1], na[2]);
  }
  if(na[0] > na[1])
  {
    axom::utilities::swap(na[0], na[1]);
  }
}

/*!
 *
 * \brief Compute volume fraction of unit cube at origin below plane.
 *
 * \param d  Plane displacement
 * \param n1 Sorted positive normal components of plane 1
 * \param n2 Sorted positive normal components of plane 2
 * \param n3 Sorted positive normal components of plane 3 (highest)
 *
 * \return Volume fraction
 *
 * \note Plane equation is (n,x) + d = 0.
 *       Cross section formulas from Youngs (1987).
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE FloatType vf_1cube(FloatType d, FloatType n1, FloatType n2, FloatType n3)
{
  FloatType vf = 0.0, dd1, dd2, dd3;
  enum
  {
    cut_triangle,
    cut_quada,
    cut_pent,
    cut_hex,
    cut_quadb
  } xsec;
  constexpr FloatType one6 = 1. / 6.;
  FloatType nbd = n1 + n2 + n3;

  int ifdhi = 0;
  if(d > 0.5 * nbd)
  {
    d = nbd - d;
    ifdhi = 1;  // symmetric at vf = 0.5
  }
  if(d <= 0.0)
  {
    vf = 0.0;
  }  // For numerical safety.
  else
  {
    // Do the vf calculation per Youngs.
    // We have already guaranteed d>0.

    // Select cross section type.
    if(d <= n1)
    {
      xsec = cut_triangle;
    }
    else if(d <= n2)
    {
      xsec = cut_quada;
    }
    else if(d <= n3)
    {
      if(d <= (n1 + n2))
      {
        xsec = cut_pent;
      }
      else
      {
        xsec = cut_quadb;
      }
    }
    else
    {
      xsec = cut_hex;
    }

    // Compute volume fraction based on type.
    switch(xsec)
    {
    case cut_triangle:

      vf = one6 * d * d * d / (n1 * n2 * n3);
      break;

    case cut_quada:

      vf = one6 * (3.0 * d * (d - n1) + n1 * n1) / (n2 * n3);
      break;

    case cut_pent:

      dd1 = d - n1;
      dd2 = d - n2;
      vf = one6 * (d * d * d - dd1 * dd1 * dd1 - dd2 * dd2 * dd2) / (n1 * n2 * n3);
      break;

    case cut_hex:

      dd1 = d - n1;
      dd2 = d - n2;
      dd3 = d - n3;
      vf = one6 * (d * d * d - dd1 * dd1 * dd1 - dd2 * dd2 * dd2 - dd3 * dd3 * dd3) / (n1 * n2 * n3);
      break;

    case cut_quadb:

      vf = (d - 0.5 * (n1 + n2)) / n3;
      break;
    default:
      break;
    }
  }
  if(ifdhi)
  {
    vf = 1.0 - vf;
  }  // vf=0.5 symmetry.
  constexpr FloatType eps = 1.0e-15;
  constexpr FloatType lower = eps;
  constexpr FloatType upper = 1. - eps;
  if(vf < lower) vf = 0.0;  // Truncate.
  if(vf > upper) vf = 1.0;  // Truncate.
  return (vf);
}

/*!
 * \brief Solve a CUBic from 4 equally spaced Points.
 *
 * \param x An array containing 4 x values.[0] and [3] are interval endpoints, [1],[2] unused.
 * \param y An array containing 4 y values.Dependent variable at 4 equally spaced points
 *
 * \return A root of the cubic. 
 *
 * \note  
 *   x[0], x[3] are interval endpoints, with function values y[0], y[3].
 *   y[1], y[2] are function at 2/3x[0] + 1/3x[3] and 1/3x[0] + 2/3x[3].
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE FloatType cub4p(const FloatType *x, const FloatType *y)
{
  FloatType dstar;
  FloatType e0 = 0.0, e1 = 0.0, e2 = 0.0, ep, em, ea, eb;
  FloatType de0, dde0, de0sq, dde0sq;
  FloatType de1, dde1, de1sq, dde1sq, demin, ddemax;
  FloatType y0, y1, y2, y3;
  FloatType w0 = 0.0, w1 = 0.0, w2;
  FloatType b, y30, y12, a;
  FloatType tolsec = 1.0e-28;
  constexpr FloatType one6 = 1. / 6.;
  int ifsec, ifbis;

  y0 = y[0];
  y1 = y[1];
  y2 = y[2];
  y3 = y[3];

  const FloatType tol2 = (y0 * y0 + y1 * y1 + y2 * y2 + y3 * y3) * tolsec;

  /* Initialize secant search.   */
  /* The default covers end cases where vks never 
         changes sign, as can occur for vol fraction 
         very near zero or one.   */
  ifbis = 1;
  if(y2 * y3 < 0.0)
  {
    e0 = one6;
    w0 = y2;
    e1 = 0.5;
    w1 = y3;
  }
  else if(y1 * y2 < 0.0)
  {
    e0 = -one6;
    w0 = y1;
    e1 = one6;
    w1 = y2;
  }
  else if(y0 * y1 < 0.0)
  {
    e0 = -0.5;
    w0 = y0;
    e1 = -one6;
    w1 = y1;
  }

  // The quantities y1, y2, y3 have already been squared.  These
  // test for not equal to 0 (clears compiler warning).
  else if(y3 * y3 <= 0.0)
  {
    e2 = 0.5;
    ifbis = 0;
  }
  else if(y2 * y2 <= 0.0)
  {
    e2 = one6;
    ifbis = 0;
  }
  else if(y1 * y1 <= 0.0)
  {
    e2 = -one6;
    ifbis = 0;
  }
  else
  {
    e2 = -0.5;
    ifbis = 0;
  }

  ifsec = 0;

  /* Initialize with one Newton iteration.   */
  if(ifbis)
  {
    /* Precompute these. */
    y30 = y3 - y0;
    y12 = y1 - y2;
    b = 4.5 * y30 + 13.5 * y12;

    ep = e0 + 0.5;
    em = e0 - 0.5;
    ea = e0 + one6;
    eb = e0 - one6;
    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    de0 = y30 + (ep + em) * a + ep * em * b;
    dde0 = 2.0 * (a + (ep + em) * b);
    de0sq = de0 * de0;
    dde0sq = dde0 * dde0;

    ep = e1 + 0.5;
    em = e1 - 0.5;
    ea = e1 + one6;
    eb = e1 - one6;
    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    de1 = y30 + (ep + em) * a + ep * em * b;
    dde1 = 2.0 * (a + (ep + em) * b);
    de1sq = de1 * de1;
    dde1sq = dde1 * dde1;

    /* This division is safe. */
    e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
  }
  else
  {
    // This clears up compiler warnings for uninitialized vars.
    e0 = e1 = w0 = w1 = 0.0;
    de0 = de0sq = dde0sq = de1 = de1sq = dde1sq = b = y30 = 0.0;
  }

  /* Do the bisection until ok for secant.   */
  while(ifbis)
  {
    ep = e2 + 0.5;
    em = e2 - 0.5;
    ea = e2 + one6;
    eb = e2 - one6;

    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);
    w2 = ep * y3 - em * y0 + ep * em * a;

    // If we have converged, quit.
    if(w2 * w2 < tol2)
    {
      ifbis = ifsec = 0;
    }

    else
    {
      // Else, do the bisection, rotate.
      // Set correct interval to (e0,e1).
      if(w1 * w2 < 0.0)
      {
        e0 = e1;
        w0 = w1;
        de0 = de1;
        dde0 = dde1;
        de0sq = de1sq;
        dde0sq = dde1sq;
      }

      e1 = e2;
      w1 = w2;

      de1 = y30 + (ep + em) * a + ep * em * b;
      dde1 = 2.0 * (a + (ep + em) * b);
      de1sq = de1 * de1;
      dde1sq = dde1 * dde1;

      // Test for secant applicability.
      if(de0 * de1 > 0.0)
      {  // dy/de same sign.

        demin = (de0sq < de1sq) ? de0sq : de1sq;
        ddemax = (dde0sq > dde1sq) ? dde0sq : dde1sq;

        // Test on second derivative.
        // This works since equation is cubic.
        if(demin > (e1 - e0) * (e1 - e0) * ddemax)
        {
          ifbis = 0;
          ifsec = 1;
        }
      }

      // Next estimate.
      if(ifsec)
        e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
      else
        e2 = 0.5 * (e0 + e1);
    }
  }

  // Refine with secant search.
  while(ifsec)
  {
    ep = e2 + 0.5;
    em = e2 - 0.5;
    ea = e2 + one6;
    eb = e2 - one6;

    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    // Evaluate the volume at e2 as w2.
    w2 = ep * y3 - em * y0 + ep * em * a;

    if(w2 * w2 < tol2)
    {
      ifsec = 0;
    }
    else
    {
      e0 = e1;
      e1 = e2;
      w0 = w1;
      w1 = w2;
      e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
    }
  }

  dstar = x[0] + (e2 + 0.5) * (x[3] - x[0]);

  return (dstar);
}

/*!
 *
 * \brief Calculate displacement of plane through unit cube located 
 *        between (1,1,1) and (2,2,2) given volume fraction.  
 *
 * \param n    Plane normal
 * \param vf13 Volume fraction
 *
 * \return Plane displacement, d, such that (n,x)+d=0 defines plane.
 *
 * \note "Lowest d" node depends on signs of n components, hence additional
 *         offsets from the lowest d node to origin.  
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE FloatType d_3cube(const FloatType n[3], FloatType vf13)
{
  int i, ik, ifdhi;
  FloatType d[5], vc[5], v[4], x[4];
  FloatType na[3], n1, n2, n3, nbd, dstar;
  FloatType one3 = 0.3333333333333333;
  FloatType two3 = 2.0 * one3;

  // Get sorted positive normal components.
  n_sort(n, na);
  n1 = na[0];
  n2 = na[1];
  n3 = na[2];
  nbd = n1 + n2 + n3;

  // Remove cases of bad volume fractions.
  if(vf13 <= 0.0)
  {
    dstar = 0.0;
  }
  else if(vf13 >= 1.0)
  {
    dstar = nbd;
  }
  else
  {
    // Set up ascending sequence of d's, where plane passes through
    // cube nodes.
    d[0] = 0.0;
    d[1] = n1;
    d[2] = n2;
    if((n1 + n2) < n3)
    {
      d[3] = n1 + n2;
    }
    else
    {
      d[3] = n3;
    }
    d[4] = 0.5 * nbd;

    vc[0] = 0.0;
    vc[4] = 0.5;

    // Symmetrize: For volume fractions > 0.5, use complements.
    if(vf13 > 0.5)
    {
      ifdhi = 1;
      vc[0] = 1.0;
      for(i = 0; i < 4; i++)
      {
        d[i] = nbd - d[i];
      }
    }
    else
    {
      ifdhi = 0;
    }

    // Get vol fraction for each d, where plane passes through node.
    for(i = 1; i <= 3; i++)
    {
      vc[i] = vf_1cube(d[i], n1, n2, n3);
    }

    // Find volume fraction's interval.
    for(ik = 0; ik < 4; ik++)
    {
      if((vc[ik] - vf13) * (vc[ik + 1] - vf13) <= 0)
      {
        break;
      }
    }

    if(ik == 4) return (d[4]);  // Exact boundary.

    // Cross section is a simple tetrahedron.
    if(ik == 0)
    {
      if(ifdhi)
      {
        dstar = d[0] - pow(6.0 * (1.0 - vf13) * n1 * n2 * n3, one3);
      }
      else
      {
        dstar = pow(6.0 * vf13 * n1 * n2 * n3, one3);
      }
    }
    else
    {
      // More complicated: set up, solve the cubic.
      x[0] = d[ik];
      x[3] = d[ik + 1];
      v[0] = vc[ik] - vf13;
      v[3] = vc[ik + 1] - vf13;
      x[1] = two3 * x[0] + one3 * x[3];
      v[1] = vf_1cube(x[1], n1, n2, n3) - vf13;
      x[2] = one3 * x[0] + two3 * x[3];
      v[2] = vf_1cube(x[2], n1, n2, n3) - vf13;

      dstar = cub4p(x, v);
    }

    if(dstar < 0.0) dstar = 0.0;
    if(dstar > nbd) dstar = nbd;
  }

  dstar += (n[0] + n[1] + n[2]);  // (1,1,1) offset: center zone.

  // Additional offset for lowest d node of cube.
  for(i = 0; i < 3; i++)
  {
    if(n[i] < 0.0)
    {
      dstar += n[i];
    }
  }
  dstar = -dstar;  // conventional (n,x)+d=0.

  return dstar;
}

/*!
 * \brief Given a plane through 3 cubed stencil, find volume fraction 
 *        inside half space away from which normal points ("below" normal)
 *        for the zones that are center, face, edges. 
 *
 * \param n The plane normal (pointing outward).
 * \param pd The -displacement from the lower left corner of the zone.
 * \param[out] vfs An array of volume fractions for the ivf zone ids.
 * \param ivf A list of stencil indices.
 * \param k The number of stencil indices in \a ivf.
 *
 * \note  Plane equation is (n,x) + d = 0.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void vf_3cube(FloatType n[3], FloatType pd, FloatType *vfs, const int *ivf, int k)
{
  // Find node of lowest d
  // as 0 or 1 offset for each coord.
  const FloatType delx0 = (n[0] < 0.0) ? 1.0 : 0.0;
  const FloatType dely0 = (n[1] < 0.0) ? 1.0 : 0.0;
  const FloatType delz0 = (n[2] < 0.0) ? 1.0 : 0.0;

  FloatType na[3];
  n_sort(n, na);

  for(int i = 0; i < k; i++)
  {
    // Compute lower left coordinate of stencil box.
    // The values are 0,1,2
    const FloatType sx = static_cast<FloatType>(ivf[i] % 3);
    const FloatType sy = static_cast<FloatType>((ivf[i] % 9) / 3);
    const FloatType sz = static_cast<FloatType>(ivf[i] / 9);

    const FloatType x0 = sx + delx0;
    const FloatType y0 = sy + dely0;
    const FloatType z0 = sz + delz0;
    const FloatType d = -(pd + (n[0] * x0 + n[1] * y0 + n[2] * z0));

    const FloatType vf = vf_1cube(d, na[0], na[1], na[2]);

    vfs[ivf[i]] = vf;
  }
}

/*!
 * \brief 2nd order 2-d mir function.
 *
 * \param vf The input volume fraction stencil. This is 9 elements.
 * \param[out] n The output normal (x,y,0) pointing away from the material.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void elvira2xy(const FloatType *vf, FloatType n[3])
{
  // These are indices into the volume fractions that pull out values in the XY plane.
  const int ivf[9] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
  // Pick some positions out of ivf.
  constexpr int center = 4;
  constexpr int bottom = 1;
  constexpr int top = 7;
  constexpr int left = 3;
  constexpr int right = 5;

  // Pick a direction.
  FloatType variance[2];
  variance[0] = ((vf[center] - vf[bottom]) * (vf[center] - vf[bottom]) +
                 (vf[center] - vf[top]) * (vf[center] - vf[top]));
  variance[1] = ((vf[center] - vf[left]) * (vf[center] - vf[left]) +
                 (vf[center] - vf[right]) * (vf[center] - vf[right]));
  const auto direction = (variance[0] < variance[1]) ? Direction::HORIZONTAL : Direction::VERTICAL;

  // Compute normals
  Result2D<FloatType> result;
  elvira2d(result, vf, ivf, direction);
  result.plane = 2;

  FloatType chmin = 0.0;
  FloatType vfs[9];
  const FloatType chfac[3] = {1.0, 0.25, 1.0};

  // Choose difference schemes.
  for(int diff = 0; diff < 3; diff++)
  {
    // Turn the normal into normalized 3D vector.
    FloatType n3[3];
    norm2d(result.normal[diff], n3);

    // Compute displacement for center stencil for this normal.
    FloatType d = d_3cube(n3, vf[center]);

    // Compute the vfs for the various boxes in the stencil using current n3.
    vf_3cube(n3, d, vfs, ivf, 9);

    // Compute diffs between computed, actual vf's.
    FloatType chisq = elvira_chisq(vf, vfs, ivf, 9);
    // Prefer central difference (since a low value will lower the error)
    chisq *= chfac[diff];

    // Keep the normal with the least error
    if((diff == 0) || (chisq < chmin))
    {
      chmin = chisq;
      n[0] = n3[0];
      n[1] = n3[1];
      n[2] = n3[2];
    }
  }
}

/*!
 * \brief Compute variance among central zone, 4 face neighbors
 *        in two directions.  Results are used to select planes,
 *        directions for 2d normals.
 *
 * \param vf The input volume fraction stencil. This is 27 elements.
 * \param ivf The volume fraction indices being used to slice the vfs.
 *
 * \return Variance among central zone.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE FloatType det_variance(const FloatType *vf, const int *ivf)
{
  /* "andfn" selects these zones from the ivf stencil.
   *  *---*---*---*
   *  |   | 7 |   |
   *  *---*---*---*
   *  | 3 | 4 | 5 |
   *  *---*---*---*
   *  |   | 1 |   |
   *  *---*---*---*
   */
  const int andfn[5] = {1, 3, 4, 5, 7};
  const int len = 5;

  FloatType avg = 0.0;
  FloatType avg2 = 0.0;
  for(int i = 0; i < len; i++)
  {
    // Pull out certain values from the plane selected by ivf.
    const FloatType x = vf[ivf[andfn[i]]];
    avg += x;
    avg2 += x * x;
  }
  avg /= len;
  avg2 /= len;

  const FloatType var = ((double)len / (len - 1.0)) * (avg2 - avg * avg);
  return var;
}

/*!
 * \brief Select plane normal that best fits stencil zones from trials.
 *
 * \param elv2d Column sums, trial 2d normals in two planes.
 * \param vf    Volume fractions in stencil.
 *
 * \note We require that both planes have central diff, or else both have 
 *       side diff, since the asymmetric case (one central, one side) gives 
 *       poor results for spheres by breaking symmetry of two planes.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void pick_elv(Result2D<FloatType> elv2d[2], const FloatType *vf)
{
  const int ivf[19] = {1, 3, 4, 5, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 21, 22, 23, 25};
  const int idiff[5] = {0, 0, 1, 2, 2};
  const int jdiff[5] = {0, 2, 1, 0, 2};
  const FloatType chfac[5] = {1.0, 1.0, 0.25, 1.0, 1.0};

  const int p0 = elv2d[0].plane;
  const int p1 = elv2d[1].plane;
  const int p2 = 3 - p0 - p1;

  int imin = 0;
  int jmin = 0;
  FloatType chmin = 0.0;

  // Only use differences symmetric over two planes.
  for(int diff = 0; diff < 5; diff++)
  {
    const int i = idiff[diff];
    const int j = jdiff[diff];

    // Get 3d trial slope.
    FloatType n30[3], n3[3];
    norm3d(elv2d[0].normal[i], elv2d[1].normal[j], n30);
    n3[p0] = n30[0];
    n3[p1] = n30[1];
    n3[p2] = n30[2];

    // Calculate plane displacement
    const auto d = d_3cube(n3, vf[13]);

    // Find volume fractions vfs
    FloatType vfs[27];
    vf_3cube(n3, d, vfs, ivf, 19);

    // Compute diffs between computed, actual vf's.
    FloatType chisq = elvira_chisq(vf, vfs, ivf, 19);

    // Preferentially use central difference.
    chisq *= chfac[diff];

    // Select trial normal with lowest error sum.
    if((diff == 0) || (chisq < chmin))
    {
      imin = i;
      jmin = j;
      chmin = chisq;
    }
  }

  // Record which differences were used.
  elv2d[0].difference_used = IntToDifference(imin);
  elv2d[1].difference_used = IntToDifference(jmin);
}

/*!
 * \brief Find estimate of missing volume in side columns.
 *
 * Corners : near : n_x, n_y same sign, far : opposite sign.
 * vsign due to volume vm = 6 h / (n_x n_y) for near corners, and
 *                       vm =-6 h / (n_x n_y) for far corners.
 * h1 = h (height of missing tet).
 *
 * This routine may be called for any two adjacent side columns.
 * Example: if is=1 in crc_sides, c10 is +y and c01 is -x (and vice versa).
 * 
 * \param c00 Center column sum.
 * \param c10 Side column sum
 * \param v10 Old missing volume for c10 
 * \param c01 Side column sum
 * \param v01 Old missing volume for c01
 * \param vma New estimate for missing volume v10.
 * \param useFar Whether to use far corner formulas.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void missvol1(FloatType c00,  // Center column sum.
                               FloatType c10,
                               FloatType v10,
                               FloatType c01,
                               FloatType v01,
                               FloatType &vma,
                               int is_far)
{
  constexpr FloatType one6 = 1. / 6.;

  const FloatType k10 = c10 + v10;
  const FloatType k01 = c01 + v01;
  const FloatType d10 = c00 - k10;
  const FloatType d01 = c00 - k01;
  const FloatType dd1 = 6.0 * d10 * d01;

  FloatType h1, vsign;
  if(is_far)
  {
    h1 = 1.5 * k10 - 0.5 * k01;
    vsign = -1.0;
  }
  else
  {
    h1 = 1.5 * k10 + 0.5 * k01 - c00;
    vsign = 1.0;
  }

  // Max tetrahedron height.
  FloatType hmax;
  {
    const FloatType hmax2 = axom::utilities::abs(d01);
    hmax = axom::utilities::abs(d10);
    hmax = axom::utilities::min(hmax, hmax2);
    hmax = axom::utilities::min(hmax, 1.);
  }

  FloatType lb = 0.0, ub = 0.0;
  bool correction = true;
  if(c10 > c00 + 0.5)
  {
    // Up correction.
    lb = 0.0;
    ub = one6 * hmax;
    h1 -= 3.0;
  }
  else if(c10 < c00 - 0.5)
  {
    // Down correction.
    lb = -one6 * hmax;
    ub = 0.0;
  }
  else
  {
    // No correction.
    correction = false;
  }

  if(correction)
  {
    const FloatType h2 = h1 * h1;
    const FloatType d = 4.5 * h2 - vsign * (dd1 - 6.0 * v10 * d01);

    // d should be between 0 and -6.
    const FloatType va = (d < 0.0) ? ((h2 * h1 - vsign * v10 * dd1) / d) : 0;

    vma -= va;
    vma = axom::utilities::max(vma, lb);
    vma = axom::utilities::min(vma, ub);
  }
}

/*!
 * \brief Solve for side missing volumes, iteratively.
 *
 * \param upcol Column sums, 3 columns each of 2 planes.
 * \param diff Difference formula ID for each plane.
 * \param[out] n2a 2d normal in each plane (output).
 * \param[out] n3a Composite 3d normal (output).
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void crc_sides(const FloatType upcol[2][3],
                                const Difference diff[2],
                                FloatType n2a[2][2],
                                FloatType n3a[3])
{
  constexpr int MAX_ITERATIONS = 100;
  FloatType c00, c10, c01;
  FloatType side[4];

  // central column.
  c00 = upcol[0][1];
  // Cyclical side column sums.
  side[0] = upcol[0][2];
  side[1] = upcol[1][2];
  side[2] = upcol[0][0];
  side[3] = upcol[1][0];

  int ncrc = 1;
  int is[2] = {0, 0};
  if(diff[0] == Difference::FORWARD)
  {
    is[0] = (diff[1] == Difference::FORWARD) ? 0 : 3;
  }
  else if(diff[0] == Difference::BACKWARD)
  {
    is[0] = (diff[1] == Difference::FORWARD) ? 1 : 2;
  }
  else
  {             // Central-Central
    is[0] = 0;  // Central: Apply correction to +x+y and -x-y pairs.
    is[1] = 2;
    ncrc = 2;
  }

  // Loop over one or two corrections.
  for(int i = 0; i < ncrc; i++)
  {
    constexpr FloatType tol = 1.0e-24;
    FloatType del = 1.0;

    // Missing volumes, b ccw from a.
    FloatType vma, vmb, vmaold, vmbold;
    vma = vmb = vmaold = vmbold = 0.0;

    int isp = is[i] + 1;
    if(isp == 4)
    {
      isp = 0;
    }

    // Rotational Orientation.
    c10 = side[is[i]];
    c01 = side[isp];

    // Select near or far corrections.
    const int is_far = ((c10 - c00) * (c01 - c00) > 0) ? 0 : 1;

    int iters = 0;
    while(del > tol)
    {
      missvol1(c00, c10, vmaold, c01, vmbold, vma, is_far);
      missvol1(c00, c01, vmbold, c10, vmaold, vmb, is_far);
      del = (vma - vmaold) * (vma - vmaold) + (vmb - vmbold) * (vmb - vmbold);
      vmaold = vma;
      vmbold = vmb;
      iters++;
      if(iters >= MAX_ITERATIONS)
      {
        del = 0.0;
      }
    }

    side[is[i]] += vma;
    side[isp] += vmb;
  }

  n2a[0][1] = n2a[1][0] = 1.0;

  for(int i = 0; i < 2; i++)
  {
    switch(diff[i])
    {
    case Difference::FORWARD:
      n2a[i][i] = axom::utilities::abs(c00 - side[i]);
      break;
    case Difference::BACKWARD:
      n2a[i][i] = axom::utilities::abs(c00 - side[2 + i]);
      break;
    case Difference::CENTRAL:
      n2a[i][i] = 0.5 * axom::utilities::abs(side[i] - side[2 + i]);
      break;
    default:
      break;
    }
  }

  norm3d(n2a[0], n2a[1], n3a);
}

/*!
 * \brief Do center column corrections by solving cubic.  
 *
 * \param vf Volume fraction of reference (center) zone.
 * \param c0 Center column sum
 * \param cx Side column sum
 * \param cy Side column sum
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE FloatType missvol2(FloatType vf, FloatType c0, FloatType cx, FloatType cy)
{
  FloatType v[4], f[4], vm;
  FloatType cs, sign;
  constexpr FloatType one3 = 1. / 3.;
  constexpr FloatType two3 = 2. / 3.;

  if(vf > 0.5)
  {
    // upper correction - convert to lower.
    sign = -1.0;
    c0 = 3.0 - c0;
    cx = 3.0 - cx;
    cy = 3.0 - cy;
  }
  else
  {
    sign = 1.0;
  }

  cs = 2.0 * c0 - 0.5 * (cx + cy);
  f[3] = cs * cs * cs;

  cs = -0.5 * (cx + cy);
  f[0] = cs * cs * cs + 6.0 * c0 * cx * cy;

  // Satisfy constraints to initiate correction.
  if((c0 < cx) && (c0 < cy) && (f[3] < 0.0) && (f[0] > 0.0))
  {
    cs += two3 * c0;
    f[1] = cs * cs * cs + 4.0 * c0 * (one3 * c0 - cx) * (one3 * c0 - cy);

    cs += two3 * c0;
    f[2] = cs * cs * cs + 2.0 * c0 * (two3 * c0 - cx) * (two3 * c0 - cy);

    v[0] = -c0;
    v[1] = -two3 * c0;
    v[2] = -one3 * c0;
    v[3] = 0.0;

    // Returns a negative value for vm (lower correction).
    vm = cub4p(v, f);

    // Convert to upper correction if vol fraction > 0.5.
    vm = vm * sign;
  }
  else
  {
    // No correction.
    vm = 0.0;
  }

  return vm;
}

/*!
 * \brief Set up center column corrections.
 *
 * \param upcen Volume fractions, zones in central column
 * \param upcol Column sums, 3 columns each of 2 planes.
 * \param diff Difference formula ID for each plane.
 * \param[out] n2a 2d normal in each plane (output).
 * \param[out] n3a Composite 3d normal (output).
 *
 * \note Adapted from code by Jeff Grandy
 *
 */
template <typename FloatType>
AXOM_HOST_DEVICE void crc_cen(const FloatType upcen[3],
                              const FloatType upcol[2][3],
                              const Difference diff[2],
                              FloatType n2a[2][2],
                              FloatType n3a[3])
{
  // No correction for central difference.
  if((diff[0] != Difference::CENTRAL) && (diff[1] != Difference::CENTRAL))
  {
    FloatType cside[2];
    for(int i = 0; i < 2; i++)
    {
      cside[i] = (diff[i] == Difference::FORWARD) ? upcol[i][2] : upcol[i][0];
    }

    FloatType c0 = upcol[0][1];

    // Get the missing volume.
    const FloatType vm = missvol2(upcen[1], c0, cside[0], cside[1]);
    c0 += vm;

    // Find the gradient.
    for(int i = 0; i < 2; i++)
    {
      if(diff[i] == Difference::FORWARD)
      {
        n2a[i][i] = axom::utilities::abs(upcol[i][2] - c0);
      }
      else
      {
        n2a[i][i] = axom::utilities::abs(c0 - upcol[i][0]);
      }
    }
    n2a[0][1] = n2a[1][0] = 1.0;
  }

  norm3d(n2a[0], n2a[1], n3a);
}

/*!
 * \brief Find volume of missing tets (extrapolating planar interface)
 *
 * \param elv2d Results from 2d elvira on two planes.
 * \param vf_cen Volume fractions, three zones in central column
 * \param[out] n3 The 3D normal.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void correct1(Result2D<FloatType> elv2d[2],
                               const FloatType vf_cen[3],
                               FloatType n3[3])
{
  constexpr FloatType one6 = 1. / 6.;
  constexpr FloatType five6 = 5. / 6.;

  const int p0 = elv2d[0].plane;
  const int p1 = elv2d[1].plane;
  const int p2 = 3 - p0 - p1;

  FloatType n30[3];
  norm3d(elv2d[0].normal[DifferenceToInt(elv2d[0].difference_used)],
         elv2d[1].normal[DifferenceToInt(elv2d[1].difference_used)],
         n30);

  // Ensure correct order.
  if((n30[2] * n30[2] >= n30[1] * n30[1]) && (n30[2] * n30[2] >= n30[0] * n30[0]) &&

     (n30[0] * n30[0] > 0.0) &&  // Do not correct if 0 component.
     (n30[1] * n30[1] > 0.0))
  {
    // Load columns, central column vol fractions.
    FloatType upcol[2][3], upcen[3];
    for(int j = 0; j < 3; j++)
    {
      upcen[j] = vf_cen[j];
      for(int i = 0; i < 2; i++)
      {
        upcol[i][j] = elv2d[i].columns[j];
      }
    }

    Difference diff[2];
    for(int i = 0; i < 2; i++)
    {
      diff[i] = elv2d[i].difference_used;
    }

    // Turn over (rotating in K-I plane) if n[2] negative.
    // Corrections apply only to magnitude of normal; signs are
    // correct from mat corner algorithm.
    if(n30[2] < 0.0)
    {
      axom::utilities::swap(upcen[0], upcen[2]);

      diff[0] = Reverse(diff[0]);  // Forward <--> Backward

      axom::utilities::swap(upcol[0][0], upcol[0][2]);
    }

    // Magnitudes of normals.
    FloatType n2a[2][2];
    for(int i = 0; i < 2; i++)
    {
      const FloatType *n2 = elv2d[i].normal[DifferenceToInt(elv2d[i].difference_used)];
      for(int k = 0; k < 2; k++)
      {
        n2a[i][k] = axom::utilities::abs(n2[k]);
      }
    }

    // Signs of components.
    FloatType nsgn[3];
    for(int k = 0; k < 3; k++)
    {
      nsgn[k] = ((n30[k] >= 0) ? 1.0 : -1.0);
    }

    FloatType n30a[3];
    if((vf_cen[1] >= one6) && (vf_cen[1] <= five6))
    {
      crc_sides(upcol, diff, n2a, n30a);
    }
    else
    {
      crc_cen(upcen, upcol, diff, n2a, n30a);
    }

    // Permute and replace signs.
    n3[p0] = n30a[0] * nsgn[0];
    n3[p1] = n30a[1] * nsgn[1];
    n3[p2] = n30a[2] * nsgn[2];
  }
  else
  {
    // No correction done. Only permute.
    n3[p0] = n30[0];
    n3[p1] = n30[1];
    n3[p2] = n30[2];
  }
}

/*!
 * \brief 2nd order 3-d mir function.
 *
 * \param vf The input volume fraction stencil. This is 27 elements.
 * \param[out] n The output normal (x,y,z) pointing away from the material.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void elvira3d(const FloatType *vf, FloatType n[3])
{
  const int stencil2d[3][9] = {{1, 4, 7, 10, 13, 16, 19, 22, 25},     // yz
                               {3, 12, 21, 4, 13, 22, 5, 14, 23},     // zx
                               {9, 10, 11, 12, 13, 14, 15, 16, 17}};  // xy

  // central columns.
  const int ccol[3][3] = {{12, 13, 14}, {10, 13, 16}, {4, 13, 22}};

  FloatType variance[3];
  for(int i = 0; i < 3; i++)
  {
    variance[i] = det_variance(vf, stencil2d[i]);
  }

  // Find dimension with smallest variance.
  int small = 0;
  small = (variance[1] < variance[small]) ? 1 : small;
  small = (variance[2] < variance[small]) ? 2 : small;

  // Get 2d normals to slopes in direction away from material
  // plane 0=yz, 1=zx, 2=xy.
  Result2D<FloatType> elv2d[2];
  int plane = small;
  for(int dir = 0; dir < 2; dir++)
  {
    plane = (plane + 1) % 3;
    elvira2d(elv2d[dir], vf, stencil2d[plane], (dir == 0) ? Direction::VERTICAL : Direction::HORIZONTAL);
    elv2d[dir].plane = plane;
  }

  // Pull out data of the selected central column
  FloatType vf_col[3];
  for(int i = 0; i < 3; i++)
  {
    vf_col[i] = vf[ccol[small][i]];
  }

  // Choose difference schemes.
  pick_elv(elv2d, vf);

  // make corrections, find 3d normal.
  correct1(elv2d, vf_col, n);
}

/*!
 * \brief Multiply the normal times the Jacobian and normalize.
 *
 * \param normal The normal to transform, nx,ny,nz
 * \param jac    The jacobian to multiply by
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename FloatType>
AXOM_HOST_DEVICE void transform(FloatType *normal, const FloatType jac[3][3])
{
  SLIC_ASSERT(normal != nullptr);

  FloatType norm = 0.0;
  FloatType dfvsum = 0.0;  // stores nx^2 + ny^2 + nz^2
  FloatType dfv[3], delfv[3];

  dfv[0] = normal[0];
  dfv[1] = normal[1];
  dfv[2] = normal[2];

  for(int f = 0; f < 3; f++)
  {
    delfv[f] = dfv[0] * jac[0][f] + dfv[1] * jac[1][f] + dfv[2] * jac[2][f];

    norm += delfv[f] * delfv[f];
    dfvsum += dfv[f] * dfv[f];
  }

  // dfvsum tests small variations in vol frac ( see mira1_youngs )
  if((dfvsum > 1.0e-24) && (norm > 0.0))
  {
    norm = 1.0 / sqrt(norm);
    for(int f = 0; f < 3; f++)
    {
      delfv[f] = delfv[f] * norm;
    }
  }
  else
  {
    delfv[0] = 1.0; /* An arbitrary slope. */
    delfv[1] = 0.0;
    delfv[2] = 0.0;
  }

  /* Store the slope.  */
  normal[0] = delfv[0];
  normal[1] = delfv[1];
  normal[2] = delfv[2];
}

/*!
 * \brief Compute the jacobian from an input stencil of x,y,z coordinates derived
 *        from zone centroids.
 *
 * \param xcst A stencil of X coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param ycst A stencil of Y coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param zcst A stencil of Z coordinate values. 9 values if ndims == 2, 27 values if ndims == 3.
 * \param ndims The number of dimensions (2 or 3).
 * \param[out] jac A 3x3 Jacobian matrix.
 *
 * \note  Adapted from code by Jeff Grandy.
 */
template <typename FloatType>
AXOM_HOST_DEVICE void computeJacobian(const FloatType *xcst,
                                      const FloatType *ycst,
                                      const FloatType *zcst,
                                      int ndims,
                                      FloatType jac[3][3])
{
  FloatType del[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}, det;
  int f, f0, f1, f2, g0, g1, g2;
  const int perm1[3] = {1, 2, 0};
  const int perm2[3] = {2, 0, 1};

  // Stencil index pairs in X, X, Y, Y, Z, Z order.
  const int idx_2D[6] = {3, 5, 1, 7, 4, 4};
  const int idx_3D[6] = {12, 14, 10, 16, 4, 22};

  const int *idx = (ndims == 3) ? idx_3D : idx_2D;

  /*
  *  Note face convention : opposite face pairs (01), (23), (45) form
  *  a left handed orientation so that det is negative for uniform
  *  positively oriented meshes.
  */
  for(f = 0; f < ndims; f++)
  {
    f1 = idx[f * 2];
    f2 = idx[f * 2 + 1];

    del[0][f] = 0.5 * (xcst[f2] - xcst[f1]);
    del[1][f] = 0.5 * (ycst[f2] - ycst[f1]);
    if(ndims == 3)
    {
      del[2][f] = 0.5 * (zcst[f2] - zcst[f1]);
    }
  }  // END for all face pairs

  /* Construct the Left Jacobian matrix. */
  /* Neglect the determinant since we renormalize. */
  /* This is jac = inverse(del) = 1/det(del) Adj(del).  */
  /* (Note g0 <--> f0 in jac).    */
  for(f0 = 0; f0 < 3; f0++)
  {
    f1 = perm1[f0];
    f2 = perm2[f0];

    for(g0 = 0; g0 < 3; g0++)
    {
      g1 = perm1[g0];
      g2 = perm2[g0];
      jac[g0][f0] = del[f1][g1] * del[f2][g2] - del[f1][g2] * del[f2][g1];
    }
  }

  /* Check the sign of the determinant. */
  det = del[0][0] * jac[0][0] + del[0][1] * jac[1][0] + del[0][2] * jac[2][0];

  /* We still need the sign of 1/det in Cramer's rule so we
    * change the sign here.  We don't need a normalized jacobian since
    * the mix slope normal is normalized later.    */
  if(det < 0.0)
  {
    for(f0 = 0; f0 < 3; f0++)
    {
      for(int f3 = 0; f3 < 3; f3++)
      {
        jac[f0][f3] *= -1.0;
      }
    }

  }  // END if negative determinant

  // Fix jacobian matrix elements that should be identity for 2D.
  if(ndims == 2)
  {
    jac[2][0] = jac[2][1] = jac[0][2] = jac[1][2] = 0.0;
    jac[2][2] = 1.0;
  }
}

//------------------------------------------------------------------------------
/*!
 * \brief Get the size of the stencil based on dimension.
 *
 * \param NDIMS The stencil size.
 *
 * \return The stencil size for dimension \a NDIMS.
 */
AXOM_HOST_DEVICE
constexpr int getStencilSize(int NDIMS) { return (NDIMS == 3) ? 27 : ((NDIMS == 2) ? 9 : 3); }

/*! 
 * \brief Base template for a class that invokes ELVIRA on various dimension data.
 */
template <int NDIMS>
struct elvira
{ };

/*!
 * \brief 2D specialization that calls elvira2xy to make normals.
 */
template <>
struct elvira<2>
{
  static constexpr int NDIMS = 2;

  /*!
   * \brief Create normals for the material interface fragments in the selected zone.
   *
   * \param matCount The number of materials in the current zone.
   * \param fragmentVFStencilStart The start of the material volume fraction stencil data for all fragments in a zone.
   * \param fragmentVectorsStart The start of the normals for all fragments in a zone.
   * \param iskip The material index to skip.
   *
   * \note Calling this function will update some vectors in the \a fragmentVectorsStart.
   */
  AXOM_HOST_DEVICE
  static void execute(int matCount,
                      const double *fragmentVFStencilStart,
                      double *fragmentVectorsStart,
                      int iskip)
  {
    constexpr int StencilSize = getStencilSize(NDIMS);
    constexpr int numVectorComponents = 3;

    const double *vol_fracs = fragmentVFStencilStart;
    double *normal = fragmentVectorsStart;

    for(int m = 0; m < matCount; m++)
    {
      if(m != iskip)
      {
        elvira2xy(vol_fracs, normal);
      }
      else
      {
        // The last fragment does not use its normal for slicing so fill
        // in a value.
        normal[0] = 1.;
        normal[1] = 0.;
        normal[2] = 0.;
      }

      // Advance to next fragment.
      vol_fracs += StencilSize;
      normal += numVectorComponents;
    }
  }
};

/*!
 * \brief 3D specialization that calls elvira3d to make normals.
 */
template <>
struct elvira<3>
{
  static constexpr int NDIMS = 3;

  /*!
   * \brief Create normals for the material interface fragments in the selected zone.
   *
   * \param matCount The number of materials in the current zone.
   * \param fragmentVFStencilStart The start of the material volume fraction stencil data for all fragments in a zone.
   * \param fragmentVectorsStart The start of the normals for all fragments in a zone.
   * \param iskip The material index to skip.
   *
   * \note Calling this function will update some vectors in the \a fragmentVectorsStart.
   */
  AXOM_HOST_DEVICE
  static void execute(int matCount,
                      const double *fragmentVFStencilStart,
                      double *fragmentVectorsStart,
                      int iskip)
  {
    constexpr int StencilSize = getStencilSize(NDIMS);
    constexpr int numVectorComponents = 3;

    const double *vol_fracs = fragmentVFStencilStart;
    double *normal = fragmentVectorsStart;

    for(int m = 0; m < matCount; m++)
    {
      if(m != iskip)
      {
        elvira3d(vol_fracs, normal);
      }
      else
      {
        normal[0] = 1.;
        normal[1] = 0.;
        normal[2] = 0.;
      }

      // Advance to next fragment.
      vol_fracs += StencilSize;
      normal += numVectorComponents;
    }
  }
};

}  // namespace elvira
}  // namespace mir
}  // namespace axom
