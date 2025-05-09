// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// NOTE: This file is meant to be included by ElviraAlgorithm.hpp after its
//       other includes so we do not include much here.

#include <iostream>

namespace axom
{
namespace mir
{
namespace elvira
{

// NOTE: Many of the functions below were adapted from code originally
//       written by Jeff Grandy. The functions have been templated on
//       "value_type" instead of using double so we can generate host
//       and device code for HIP without resulting in  multiple function
//       definitions at link time.

enum class Direction
{
  VERTICAL = 0,
  HORIZONTAL = 1
};

/// This struct stores the normal computed for a plane of stencil data.
template <typename value_type>
struct Result2D
{
  int plane {2};
  int difference_used {0};
  value_type columns[3] {0., 0., 0.};
  value_type normal[3][2] {{0., 0.}, {0., 0.}, {0., 0.}};
};

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
template <typename value_type>
AXOM_HOST_DEVICE void elvira2d(Result2D<value_type> &result,
                               const value_type *vf,
                               const int *ivf,
                               Direction direction)
{
  const value_type jb = vf[ivf[0]] + vf[ivf[1]] + vf[ivf[2]];  // bottom row
  const value_type jm = vf[ivf[3]] + vf[ivf[4]] + vf[ivf[5]];  // middle row
  const value_type jt = vf[ivf[6]] + vf[ivf[7]] + vf[ivf[8]];  // top row

  const value_type il = vf[ivf[0]] + vf[ivf[3]] + vf[ivf[6]];  // left column
  const value_type im = vf[ivf[1]] + vf[ivf[4]] + vf[ivf[7]];  // middle column
  const value_type ir = vf[ivf[2]] + vf[ivf[5]] + vf[ivf[8]];  // right column

  value_type slope[3] {0., 0., 0.};  // backward, cent, and forward differencing
  value_type slopev = 0.;            // Sets overall sign of normal.
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

  const value_type sign = (slopev > 0.0) ? (1.) : (-1.);

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
template <typename value_type>
AXOM_HOST_DEVICE inline value_type elvira_chisq(const value_type *vf,
                                                const value_type *vfs,
                                                const int *ivf,
                                                int k)
{
  value_type chisq = 0.0;
  for(; k--;)
  {
    const value_type c = vf[ivf[k]] - vfs[ivf[k]];
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
template <typename value_type>
AXOM_HOST_DEVICE inline void norm2d(value_type n2[2], value_type n3[3])
{
  constexpr value_type PTINY = 1.e-15;
  // Compute 2D vector magnitude.
  value_type magnitude = sqrt(n2[0] * n2[0] + n2[1] * n2[1]);
  value_type magnitude2 = magnitude * magnitude;
  value_type magnitudeReciprocal = magnitude / (magnitude2 + PTINY);
  // normalize vector, make it 3D
  n3[0] = n2[0] * magnitudeReciprocal;
  n3[1] = n2[1] * magnitudeReciprocal;
  n3[2] = (magnitude2 <= 0.) ? 1.0 : 0.0;
}

/*!
 * \brief Return sorted, positive components of vector length 3.
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename value_type>
AXOM_HOST_DEVICE inline void n_sort(const value_type n[3],  // input vector
                                    value_type na[3])       // output vector
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
template <typename value_type>
AXOM_HOST_DEVICE value_type vf_1cube(value_type d, value_type n1, value_type n2, value_type n3)
{
  value_type vf = 0.0, dd1, dd2, dd3;
  enum
  {
    cut_triangle,
    cut_quada,
    cut_pent,
    cut_hex,
    cut_quadb
  } xsec;
  constexpr value_type one6 = 1. / 6.;
  value_type nbd = n1 + n2 + n3;

  int ifdhi = 0;
  if(d > 0.5 * nbd)
  {
    d = nbd - d;
    ifdhi = 1; /* symmetric at vf = 0.5 */
  }
  if(d <= 0.0)
  {
    vf = 0.0;
  } /* For numerical safety.   */
  else
  {
    /* Do the vf calculation per Youngs.   */
    /* We have already guaranteed d>0. */

    /* Select cross section type.   */
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

    /* Compute volume fraction based on type.   */
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
  } /* vf=0.5 symmetry.   */
  constexpr value_type eps = 1.0e-15;
  constexpr value_type lower = eps;
  constexpr value_type upper = 1. - eps;
  if(vf < lower) vf = 0.0; /* Truncate. */
  if(vf > upper) vf = 1.0; /* Truncate. */
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
template <typename value_type>
AXOM_HOST_DEVICE value_type cub4p(const value_type *x, const value_type *y)
{
  value_type dstar;
  value_type e0 = 0.0, e1 = 0.0, e2 = 0.0, ep, em, ea, eb;
  value_type de0, dde0, de0sq, dde0sq;
  value_type de1, dde1, de1sq, dde1sq, demin, ddemax;
  value_type y0, y1, y2, y3;
  value_type w0 = 0.0, w1 = 0.0, w2;
  value_type b, y30, y12, a, tol2;
  value_type tolsec = 1.0e-28;
  constexpr value_type one6 = 1. / 6.;
  int ifsec, ifbis;

  y0 = y[0];
  y1 = y[1];
  y2 = y[2];
  y3 = y[3];

  tol2 = (y0 * y0 + y1 * y1 + y2 * y2 + y3 * y3) * tolsec;

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

    /* If we have converged, quit. */
    if(w2 * w2 < tol2)
    {
      ifbis = ifsec = 0;
    }

    else
    {
      /* Else, do the bisection, rotate. */
      /* Set correct interval to (e0,e1). */
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

      /* Test for secant applicability.   */
      if(de0 * de1 > 0.0)
      { /* dy/de same sign. */

        demin = (de0sq < de1sq) ? de0sq : de1sq;
        ddemax = (dde0sq > dde1sq) ? dde0sq : dde1sq;

        /* Test on second derivative.   */
        /* This works since equation is cubic.  */
        if(demin > (e1 - e0) * (e1 - e0) * ddemax)
        {
          ifbis = 0;
          ifsec = 1;
        }
      }

      /* Next estimate.    */
      if(ifsec)
        e2 = e0 + (e1 - e0) * w0 / (w0 - w1);
      else
        e2 = 0.5 * (e0 + e1);
    }
  }

  /* Refine with secant search.   */
  while(ifsec)
  {
    ep = e2 + 0.5;
    em = e2 - 0.5;
    ea = e2 + one6;
    eb = e2 - one6;

    a = 4.5 * (ep * y3 - em * y0) + 13.5 * (eb * y1 - ea * y2);

    /* Evaluate the volume at e2 as w2. */
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
template <typename value_type>
AXOM_HOST_DEVICE value_type d_3cube(const value_type n[3], value_type vf13)
{
  int i, ik, ifdhi;
  value_type d[5], vc[5], v[4], x[4];
  value_type na[3], n1, n2, n3, nbd, dstar;
  value_type one3 = 0.3333333333333333;
  value_type two3 = 2.0 * one3;

  /* Get sorted positive normal components.   */
  n_sort(n, na);
  n1 = na[0];
  n2 = na[1];
  n3 = na[2];
  nbd = n1 + n2 + n3;

  /* Remove cases of bad volume fractions. */
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
    /* Set up ascending sequence of d's, where plane passes through 
      cube nodes.    */
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

    /* Symmetrize: For volume fractions > 0.5, use complements.   */
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

    /* Get vol fraction for each d, where plane passes through node. */
    for(i = 1; i <= 3; i++)
    {
      vc[i] = vf_1cube(d[i], n1, n2, n3);
    }

    /* Find volume fraction's interval.   */
    for(ik = 0; ik < 4; ik++)
    {
      if((vc[ik] - vf13) * (vc[ik + 1] - vf13) <= 0)
      {
        break;
      }
    }

    if(ik == 4) return (d[4]); /* Exact boundary. */

    /* Cross section is a simple tetrahedron.   */
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

    /* More complicated: set up, solve the cubic.   */
    else
    {
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

  dstar += (n[0] + n[1] + n[2]); /* (1,1,1) offset: center zone. */

  /* Additional offset for lowest d node of cube. */
  for(i = 0; i < 3; i++)
    if(n[i] < 0.0)
    {
      dstar += n[i];
    }

  dstar = -dstar; /* conventional (n,x)+d=0.   */
  return (dstar);
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
template <typename value_type>
AXOM_HOST_DEVICE void vf_3cube(value_type n[3], value_type pd, value_type *vfs, const int *ivf, int k)
{
  /* find node of lowest d */
  /* as 0 or 1 offset for each coord.  */
  const value_type delx0 = (n[0] < 0.0) ? 1.0 : 0.0;
  const value_type dely0 = (n[1] < 0.0) ? 1.0 : 0.0;
  const value_type delz0 = (n[2] < 0.0) ? 1.0 : 0.0;

  value_type na[3];
  n_sort(n, na);

  for(int i = 0; i < k; i++)
  {
    // Compute lower left coordinate of stencil box.
    // The values are 0,1,2
    const value_type sx = static_cast<value_type>(ivf[i] % 3);
    const value_type sy = static_cast<value_type>((ivf[i] % 9) / 3);
    const value_type sz = static_cast<value_type>(ivf[i] / 9);

    const value_type x0 = sx + delx0;
    const value_type y0 = sy + dely0;
    const value_type z0 = sz + delz0;
    const value_type d = -(pd + (n[0] * x0 + n[1] * y0 + n[2] * z0));

    const value_type vf = vf_1cube(d, na[0], na[1], na[2]);

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
template <typename value_type>
AXOM_HOST_DEVICE void elvira2xy(const value_type *vf, value_type n[3])
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
  value_type variance[2];
  variance[0] = ((vf[center] - vf[bottom]) * (vf[center] - vf[bottom]) +
                 (vf[center] - vf[top]) * (vf[center] - vf[top]));
  variance[1] = ((vf[center] - vf[left]) * (vf[center] - vf[left]) +
                 (vf[center] - vf[right]) * (vf[center] - vf[right]));
  const auto direction = (variance[0] < variance[1]) ? Direction::HORIZONTAL : Direction::VERTICAL;

  // Compute normals
  Result2D<value_type> result;
  elvira2d(result, vf, ivf, direction);
  result.plane = 2;

  value_type chmin = 0.0;
  value_type vfs[9];
  const value_type chfac[3] = {1.0, 0.25, 1.0};

  // Choose difference schemes.
  for(int diff = 0; diff < 3; diff++)
  {
    // Turn the normal into normalized 3D vector.
    value_type n3[3];
    norm2d(result.normal[diff], n3);

    // Compute displacement for center stencil for this normal.
    value_type d = d_3cube(n3, vf[center]);

    // Compute the vfs for the various boxes in the stencil using current n3.
    vf_3cube(n3, d, vfs, ivf, 9);

    // Compute diffs between computed, actual vf's.
    value_type chisq = elvira_chisq(vf, vfs, ivf, 9);
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
 * \brief Multiply the normal times the Jacobian and normalize.
 *
 * \param normal The normal to transform, nx,ny,nz
 * \param jac    The jacobian to multiply by
 *
 * \note Adapted from code by Jeff Grandy
 */
template <typename value_type>
AXOM_HOST_DEVICE void transform(value_type *normal, const value_type jac[3][3])
{
  SLIC_ASSERT(normal != nullptr);

  value_type norm = 0.0;
  value_type dfvsum = 0.0;  // stores nx^2 + ny^2 + nz^2
  value_type dfv[3], delfv[3];

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
template <typename value_type>
AXOM_HOST_DEVICE void computeJacobian(const value_type *xcst,
                                      const value_type *ycst,
                                      const value_type *zcst,
                                      int ndims,
                                      value_type jac[3][3])
{
  value_type del[3][3] = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 1.}}, det;
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

#if defined(AXOM_ELVIRA_DEBUG) && !defined(AXOM_DEVICE_CODE)
  std::cout << "xcst={";
  for(int i = 0; i < 9; i++)
  {
    std::cout << xcst[i] << ", ";
  }
  std::cout << "}, ycst={";
  for(int i = 0; i < 9; i++)
  {
    std::cout << ycst[i] << ", ";
  }
  std::cout << "}, zcst={";
  for(int i = 0; i < 9; i++)
  {
    std::cout << zcst[i] << ", ";
  }

  std::cout << "}, del={";
  for(int row = 0; row < 3; row++)
  {
    for(int col = 0; col < 3; col++)
    {
      std::cout << del[row][col] << ", ";
    }
  }
  std::cout << "}, jac={";
  for(int row = 0; row < 3; row++)
  {
    for(int col = 0; col < 3; col++)
    {
      std::cout << jac[row][col] << ", ";
    }
  }
  std::cout << "}, ndims=" << ndims << ", det=" << det << std::endl;
#endif

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
  static void execute(int AXOM_UNUSED_PARAM(matCount),
                      const double *AXOM_UNUSED_PARAM(fragmentVFStencilStart),
                      double *AXOM_UNUSED_PARAM(fragmentVectorsStart),
                      int AXOM_UNUSED_PARAM(iskip))
  {
#if 0
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
#endif
  }
};

}  // namespace elvira
}  // namespace mir
}  // namespace axom
