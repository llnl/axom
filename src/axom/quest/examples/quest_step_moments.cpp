// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_step_moments.cpp
 * \brief Example that reads a STEP file and computes geometric moments of the
 *        corresponding trimmed NURBS surface and/or enclosed volume.
 *
 * The example loads a BRep with Quest's STEP reader, then evaluates
 * monomial moments \f$\int x^i y^j z^k \, dS\f$ and/or
 * \f$\int x^i y^j z^k \, dV\f$ for all nonnegative exponent triples with
 * total degree \f$i + j + k \le n\f$, where \a n is user supplied.
 * It also derives centroids, inertia tensors, and volume-based principal-axis
 * fit proxies from the first- and second-order moments.
 *
 * \note Volume moments are only geometrically meaningful when the STEP model
 * represents a closed, consistently oriented boundary.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/quest.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace mint = axom::mint;
namespace primal = axom::primal;
namespace slic = axom::slic;

namespace
{

using Point3D = primal::Point<double, 3>;
using Vector3D = primal::Vector<double, 3>;
using OBB3D = primal::OrientedBoundingBox<double, 3>;
using TriMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
using PatchArray = axom::quest::STEPReader::PatchArray;
using MomentKey = std::tuple<int, int, int>;

enum class IntegralMode
{
  SURFACE,
  VOLUME,
  BOTH
};

const std::map<std::string, IntegralMode> s_validIntegralModes {{"surface", IntegralMode::SURFACE},
                                                                {"volume", IntegralMode::VOLUME},
                                                                {"both", IntegralMode::BOTH}};

constexpr double eps = 1e-14;

const char* integral_mode_name(IntegralMode mode)
{
  switch(mode)
  {
  case IntegralMode::SURFACE:
    return "surface";
  case IntegralMode::VOLUME:
    return "volume";
  case IntegralMode::BOTH:
    return "both";
  }

  return "unknown";
}

bool should_compute_surface(IntegralMode mode)
{
  return mode == IntegralMode::SURFACE || mode == IntegralMode::BOTH;
}

bool should_compute_volume(IntegralMode mode)
{
  return mode == IntegralMode::VOLUME || mode == IntegralMode::BOTH;
}

struct MomentIndex
{
  int degree {};
  int px {};
  int py {};
  int pz {};
};

struct MomentEntry
{
  MomentIndex index {};
  double value {};
};

struct MomentSet
{
  std::vector<MomentEntry> requested_entries;
  std::map<MomentKey, double> values;
};

struct InertiaTensor
{
  double xx {};
  double yy {};
  double zz {};
  double xy {};
  double xz {};
  double yz {};
};

struct MassProperties
{
  bool valid {false};
  double measure {};
  Point3D centroid {};
  InertiaTensor inertia_origin {};
  InertiaTensor inertia_centroid {};
};

struct PrincipalFrame
{
  bool valid {false};
  double measure {};
  Point3D centroid {};
  double principal_inertia[3] {};
  double principal_second_moments[3] {};
  Vector3D axes[3] {};
};

struct EllipsoidFit
{
  bool valid {false};
  Point3D centroid {};
  Vector3D axes[3] {};
  double radii[3] {};
  double volume {};
};

struct ObbFit
{
  bool valid {false};
  OBB3D box {};
  double volume {};
};

std::vector<MomentIndex> enumerate_moment_indices(int max_degree)
{
  std::vector<MomentIndex> moments;
  for(int degree = 0; degree <= max_degree; ++degree)
  {
    for(int px = 0; px <= degree; ++px)
    {
      for(int py = 0; py <= degree - px; ++py)
      {
        const int pz = degree - px - py;
        moments.push_back(MomentIndex {degree, px, py, pz});
      }
    }
  }

  return moments;
}

double ipow(double base, int exponent)
{
  double result = 1.0;
  for(int i = 0; i < exponent; ++i)
  {
    result *= base;
  }
  return result;
}

double evaluate_monomial(const Point3D& x, const MomentIndex& idx)
{
  return ipow(x[0], idx.px) * ipow(x[1], idx.py) * ipow(x[2], idx.pz);
}

template <typename Integrator>
MomentSet compute_moments(int requested_max_degree, int computed_max_degree, Integrator&& integrator)
{
  MomentSet moments;
  const auto indices = enumerate_moment_indices(computed_max_degree);
  moments.requested_entries.reserve(indices.size());

  for(const auto& idx : indices)
  {
    const double value = integrator(idx);
    moments.values[MomentKey {idx.px, idx.py, idx.pz}] = value;

    if(idx.degree <= requested_max_degree)
    {
      moments.requested_entries.push_back(MomentEntry {idx, value});
    }
  }

  return moments;
}

double get_moment(const MomentSet& moments, int px, int py, int pz)
{
  const auto it = moments.values.find(MomentKey {px, py, pz});
  SLIC_ASSERT(it != moments.values.end());
  return it->second;
}

InertiaTensor shift_to_centroid(const InertiaTensor& inertia_origin,
                                double measure,
                                const Point3D& centroid)
{
  InertiaTensor shifted = inertia_origin;

  shifted.xx -= measure * (centroid[1] * centroid[1] + centroid[2] * centroid[2]);
  shifted.yy -= measure * (centroid[0] * centroid[0] + centroid[2] * centroid[2]);
  shifted.zz -= measure * (centroid[0] * centroid[0] + centroid[1] * centroid[1]);
  shifted.xy += measure * centroid[0] * centroid[1];
  shifted.xz += measure * centroid[0] * centroid[2];
  shifted.yz += measure * centroid[1] * centroid[2];

  return shifted;
}

MassProperties compute_mass_properties(const MomentSet& moments)
{
  MassProperties props;
  props.measure = get_moment(moments, 0, 0, 0);
  if(std::abs(props.measure) <= eps)
  {
    return props;
  }

  props.valid = true;
  props.centroid[0] = get_moment(moments, 1, 0, 0) / props.measure;
  props.centroid[1] = get_moment(moments, 0, 1, 0) / props.measure;
  props.centroid[2] = get_moment(moments, 0, 0, 1) / props.measure;

  const double xx = get_moment(moments, 2, 0, 0);
  const double yy = get_moment(moments, 0, 2, 0);
  const double zz = get_moment(moments, 0, 0, 2);
  const double xy = get_moment(moments, 1, 1, 0);
  const double xz = get_moment(moments, 1, 0, 1);
  const double yz = get_moment(moments, 0, 1, 1);

  props.inertia_origin.xx = yy + zz;
  props.inertia_origin.yy = xx + zz;
  props.inertia_origin.zz = xx + yy;
  props.inertia_origin.xy = -xy;
  props.inertia_origin.xz = -xz;
  props.inertia_origin.yz = -yz;

  props.inertia_centroid = shift_to_centroid(props.inertia_origin, props.measure, props.centroid);

  return props;
}

PrincipalFrame compute_principal_frame(const MassProperties& props)
{
  constexpr int ndims = 3;

  PrincipalFrame frame;
  if(!props.valid || std::abs(props.measure) <= eps)
  {
    return frame;
  }

  axom::numerics::Matrix<double> inertia(ndims, ndims);
  inertia(0, 0) = props.inertia_centroid.xx;
  inertia(0, 1) = props.inertia_centroid.xy;
  inertia(0, 2) = props.inertia_centroid.xz;
  inertia(1, 0) = props.inertia_centroid.xy;
  inertia(1, 1) = props.inertia_centroid.yy;
  inertia(1, 2) = props.inertia_centroid.yz;
  inertia(2, 0) = props.inertia_centroid.xz;
  inertia(2, 1) = props.inertia_centroid.yz;
  inertia(2, 2) = props.inertia_centroid.zz;

  axom::numerics::Matrix<double> eigenvectors(ndims, ndims);
  double eigenvalues[ndims] {};
  const int rc = axom::numerics::jacobi_eigensolve(inertia, eigenvectors, eigenvalues);
  if(rc != axom::numerics::JACOBI_EIGENSOLVE_SUCCESS)
  {
    return frame;
  }

  frame.valid = true;
  frame.measure = props.measure;
  frame.centroid = props.centroid;

  double inertia_trace = 0.0;
  for(int i = 0; i < ndims; ++i)
  {
    frame.principal_inertia[i] = eigenvalues[i];
    inertia_trace += eigenvalues[i];

    frame.axes[i][0] = eigenvectors(0, i);
    frame.axes[i][1] = eigenvectors(1, i);
    frame.axes[i][2] = eigenvectors(2, i);
  }

  for(int i = 0; i < ndims; ++i)
  {
    frame.principal_second_moments[i] = 0.5 * (inertia_trace - 2.0 * frame.principal_inertia[i]);
  }

  return frame;
}

bool compute_shape_dimensions(const PrincipalFrame& frame, double factor, double (&dims)[3])
{
  constexpr int ndims = 3;
  constexpr double min_scale = 1.0;

  if(!frame.valid || std::abs(frame.measure) <= eps)
  {
    return false;
  }

  const double scale = std::max({min_scale,
                                 std::abs(frame.measure),
                                 std::abs(frame.principal_second_moments[0]),
                                 std::abs(frame.principal_second_moments[1]),
                                 std::abs(frame.principal_second_moments[2])});
  const double tol = 1e-12 * scale;

  for(int i = 0; i < ndims; ++i)
  {
    const double dim_sq = factor * frame.principal_second_moments[i] / frame.measure;
    if(dim_sq < -tol)
    {
      return false;
    }

    dims[i] = std::sqrt(std::max(dim_sq, 0.0));
  }

  return true;
}

double compute_ellipsoid_volume(const double (&radii)[3])
{
  return (4.0 / 3.0) * M_PI * radii[0] * radii[1] * radii[2];
}

double compute_obb_volume(const Vector3D& extents)
{
  return 8.0 * extents[0] * extents[1] * extents[2];
}

double compute_volume_ratio(double fit_volume, double reference_volume)
{
  if(std::abs(reference_volume) <= eps)
  {
    return 0.0;
  }

  return fit_volume / reference_volume;
}

double compute_effective_density(double mass, double geometric_volume)
{
  if(std::abs(geometric_volume) <= eps)
  {
    return 0.0;
  }

  return mass / geometric_volume;
}

EllipsoidFit compute_inertia_matched_ellipsoid(const PrincipalFrame& frame)
{
  EllipsoidFit fit;
  double radii[3] {};
  if(!compute_shape_dimensions(frame, 5.0, radii))
  {
    return fit;
  }

  fit.valid = true;
  fit.centroid = frame.centroid;
  for(int i = 0; i < 3; ++i)
  {
    fit.axes[i] = frame.axes[i];
    fit.radii[i] = radii[i];
  }
  fit.volume = compute_ellipsoid_volume(fit.radii);

  return fit;
}

EllipsoidFit scale_ellipsoid_to_volume(const EllipsoidFit& fit, double target_volume)
{
  EllipsoidFit scaled;
  if(!fit.valid || std::abs(fit.volume) <= eps || target_volume <= 0.0)
  {
    return scaled;
  }

  const double scale = std::cbrt(target_volume / fit.volume);
  scaled.valid = true;
  scaled.centroid = fit.centroid;
  for(int i = 0; i < 3; ++i)
  {
    scaled.axes[i] = fit.axes[i];
    scaled.radii[i] = scale * fit.radii[i];
  }
  scaled.volume = compute_ellipsoid_volume(scaled.radii);

  return scaled;
}

ObbFit compute_inertia_matched_obb(const PrincipalFrame& frame)
{
  ObbFit fit;
  double extents_data[3] {};
  if(!compute_shape_dimensions(frame, 3.0, extents_data))
  {
    return fit;
  }

  Vector3D extents;
  for(int i = 0; i < 3; ++i)
  {
    extents[i] = extents_data[i];
  }

  fit.box = OBB3D(frame.centroid, frame.axes, extents);
  fit.valid = fit.box.isValid();
  fit.volume = compute_obb_volume(extents);

  return fit;
}

ObbFit scale_obb_to_volume(const ObbFit& fit, double target_volume)
{
  ObbFit scaled;
  if(!fit.valid || std::abs(fit.volume) <= eps || target_volume <= 0.0)
  {
    return scaled;
  }

  const double scale = std::cbrt(target_volume / fit.volume);
  const auto& centroid = fit.box.getCentroid();
  const auto* axes = fit.box.getAxes();
  const auto& extents = fit.box.getExtents();

  Vector3D scaled_axes[3];
  Vector3D scaled_extents;
  for(int i = 0; i < 3; ++i)
  {
    scaled_axes[i] = axes[i];
    scaled_extents[i] = scale * extents[i];
  }

  scaled.box = OBB3D(centroid, scaled_axes, scaled_extents);
  scaled.valid = scaled.box.isValid();
  scaled.volume = compute_obb_volume(scaled_extents);

  return scaled;
}

std::string make_ellipsoid_vtk_filename(const std::string& prefix, const std::string& variant)
{
  return axom::fmt::format("{}_{}.vtk", axom::utilities::string::removeSuffix(prefix, ".vtk"), variant);
}

Point3D transform_unit_sphere_point(const EllipsoidFit& fit, double x, double y, double z)
{
  Point3D point = fit.centroid;
  const double local_coords[3] {x, y, z};

  for(int axis = 0; axis < 3; ++axis)
  {
    for(int dim = 0; dim < 3; ++dim)
    {
      point[dim] += fit.radii[axis] * local_coords[axis] * fit.axes[axis][dim];
    }
  }

  return point;
}

bool write_ellipsoid_vtk(const EllipsoidFit& fit,
                         const std::string& file_path,
                         int theta_resolution = 48,
                         int phi_resolution = 24)
{
  if(!fit.valid || file_path.empty())
  {
    return false;
  }

  SLIC_ASSERT(theta_resolution >= 3);
  SLIC_ASSERT(phi_resolution >= 4);

  const int num_rings = phi_resolution - 2;
  const int total_nodes = 2 + theta_resolution * num_rings;
  const int total_cells = 2 * theta_resolution * (phi_resolution - 2);

  TriMesh mesh(3, mint::TRIANGLE);
  mesh.reserve(total_nodes, total_cells);

  auto append_node = [&mesh, &fit](double x, double y, double z) {
    const Point3D point = transform_unit_sphere_point(fit, x, y, z);
    mesh.appendNode(point[0], point[1], point[2]);
  };

  append_node(0.0, 0.0, 1.0);
  append_node(0.0, 0.0, -1.0);

  for(int i = 0; i < theta_resolution; ++i)
  {
    const double theta = 2.0 * M_PI * static_cast<double>(i) / theta_resolution;
    for(int j = 1; j <= num_rings; ++j)
    {
      const double phi = M_PI * static_cast<double>(j) / (phi_resolution - 1);
      const double sin_phi = std::sin(phi);
      append_node(std::cos(theta) * sin_phi, std::sin(theta) * sin_phi, std::cos(phi));
    }
  }

  auto ring_node = [num_rings](int theta_idx, int ring_idx) {
    return 2 + theta_idx * num_rings + ring_idx;
  };

  axom::IndexType cell[3];

  for(int i = 0; i < theta_resolution; ++i)
  {
    const int next_i = (i + 1) % theta_resolution;
    cell[0] = 0;
    cell[1] = ring_node(next_i, 0);
    cell[2] = ring_node(i, 0);
    mesh.appendCell(cell);
  }

  for(int ring = 0; ring < num_rings - 1; ++ring)
  {
    for(int i = 0; i < theta_resolution; ++i)
    {
      const int next_i = (i + 1) % theta_resolution;
      cell[0] = ring_node(i, ring);
      cell[1] = ring_node(i, ring + 1);
      cell[2] = ring_node(next_i, ring);
      mesh.appendCell(cell);

      cell[0] = ring_node(next_i, ring);
      cell[1] = ring_node(i, ring + 1);
      cell[2] = ring_node(next_i, ring + 1);
      mesh.appendCell(cell);
    }
  }

  const int last_ring = num_rings - 1;
  for(int i = 0; i < theta_resolution; ++i)
  {
    const int next_i = (i + 1) % theta_resolution;
    cell[0] = 1;
    cell[1] = ring_node(i, last_ring);
    cell[2] = ring_node(next_i, last_ring);
    mesh.appendCell(cell);
  }

  return mint::write_vtk(&mesh, file_path) == 0;
}

std::string format_real(double value) { return axom::fmt::format("{:.16e}", value); }

std::string format_point_inline(const Point3D& point)
{
  return axom::fmt::format("[{}, {}, {}]",
                           format_real(point[0]),
                           format_real(point[1]),
                           format_real(point[2]));
}

std::string format_vector_inline(const Vector3D& vector)
{
  return axom::fmt::format("[{}, {}, {}]",
                           format_real(vector[0]),
                           format_real(vector[1]),
                           format_real(vector[2]));
}

std::string format_triplet_inline(double a, double b, double c)
{
  return axom::fmt::format("[{}, {}, {}]", format_real(a), format_real(b), format_real(c));
}

void append_line(std::string& output, int indent, const std::string& line)
{
  output += std::string(indent * 2, ' ');
  output += line;
  output += '\n';
}

void append_inertia_tensor_yaml(std::string& output,
                                int indent,
                                const char* key,
                                const InertiaTensor& tensor)
{
  append_line(output, indent, axom::fmt::format("{}:", key));
  append_line(output, indent + 1, axom::fmt::format("xx: {}", format_real(tensor.xx)));
  append_line(output, indent + 1, axom::fmt::format("yy: {}", format_real(tensor.yy)));
  append_line(output, indent + 1, axom::fmt::format("zz: {}", format_real(tensor.zz)));
  append_line(output, indent + 1, axom::fmt::format("xy: {}", format_real(tensor.xy)));
  append_line(output, indent + 1, axom::fmt::format("xz: {}", format_real(tensor.xz)));
  append_line(output, indent + 1, axom::fmt::format("yz: {}", format_real(tensor.yz)));
}

void append_axes_yaml(std::string& output, int indent, const Vector3D (&axes)[3])
{
  append_line(output, indent, "axes:");
  for(int i = 0; i < 3; ++i)
  {
    append_line(output, indent + 1, axom::fmt::format("- {}", format_vector_inline(axes[i])));
  }
}

void append_axes_yaml(std::string& output, int indent, const Vector3D* axes)
{
  append_line(output, indent, "axes:");
  for(int i = 0; i < 3; ++i)
  {
    append_line(output, indent + 1, axom::fmt::format("- {}", format_vector_inline(axes[i])));
  }
}

void append_moment_entries_yaml(std::string& output, int indent, const MomentSet& moments)
{
  append_line(output, indent, "raw_moments:");
  append_line(output,
              indent + 1,
              axom::fmt::format("measure: {}", format_real(get_moment(moments, 0, 0, 0))));

  if(moments.requested_entries.empty())
  {
    append_line(output, indent + 1, "entries: []");
    return;
  }

  append_line(output, indent + 1, "entries:");
  for(const auto& entry : moments.requested_entries)
  {
    append_line(output, indent + 2, axom::fmt::format("- degree: {}", entry.index.degree));
    append_line(
      output,
      indent + 3,
      axom::fmt::format("exponents: [{}, {}, {}]", entry.index.px, entry.index.py, entry.index.pz));
    append_line(output, indent + 3, axom::fmt::format("value: {}", format_real(entry.value)));
  }
}

void append_mass_properties_yaml(std::string& output, int indent, const MassProperties& props)
{
  append_line(output, indent, "derived:");
  if(!props.valid)
  {
    append_line(output, indent + 1, "available: false");
    append_line(output, indent + 1, "reason: measure_is_zero");
    return;
  }

  append_line(output, indent + 1, "available: true");
  append_line(output,
              indent + 1,
              axom::fmt::format("centroid: {}", format_point_inline(props.centroid)));
  append_inertia_tensor_yaml(output, indent + 1, "inertia_origin", props.inertia_origin);
  append_inertia_tensor_yaml(output, indent + 1, "inertia_centroid", props.inertia_centroid);
}

void append_principal_frame_yaml(std::string& output, int indent, const PrincipalFrame& frame)
{
  append_line(output, indent, "principal_frame:");
  if(!frame.valid)
  {
    append_line(output, indent + 1, "available: false");
    append_line(output, indent + 1, "reason: eigensolve_failed_or_measure_is_zero");
    return;
  }

  append_line(output, indent + 1, "available: true");
  append_line(output,
              indent + 1,
              axom::fmt::format("principal_inertia: {}",
                                format_triplet_inline(frame.principal_inertia[0],
                                                      frame.principal_inertia[1],
                                                      frame.principal_inertia[2])));
  append_line(output,
              indent + 1,
              axom::fmt::format("principal_second_moments: {}",
                                format_triplet_inline(frame.principal_second_moments[0],
                                                      frame.principal_second_moments[1],
                                                      frame.principal_second_moments[2])));
  append_axes_yaml(output, indent + 1, &frame.axes[0]);
}

void append_ellipsoid_fit_yaml(std::string& output,
                               int indent,
                               const char* key,
                               const EllipsoidFit& fit,
                               double reference_volume,
                               const char* assumption,
                               const std::string& vtk_file)
{
  append_line(output, indent, axom::fmt::format("{}:", key));
  if(!fit.valid)
  {
    append_line(output, indent + 1, "available: false");
    append_line(output, indent + 1, "reason: invalid_second_moments");
    return;
  }

  append_line(output, indent + 1, "available: true");
  if(assumption != nullptr)
  {
    append_line(output, indent + 1, axom::fmt::format("assumption: '{}'", assumption));
  }
  if(!vtk_file.empty())
  {
    append_line(output, indent + 1, axom::fmt::format("vtk_file: '{}'", vtk_file));
  }
  append_line(output, indent + 1, axom::fmt::format("center: {}", format_point_inline(fit.centroid)));
  append_line(output,
              indent + 1,
              axom::fmt::format("semiaxes: {}",
                                format_triplet_inline(fit.radii[0], fit.radii[1], fit.radii[2])));
  append_axes_yaml(output, indent + 1, &fit.axes[0]);
  append_line(output, indent + 1, axom::fmt::format("geometric_volume: {}", format_real(fit.volume)));
  append_line(output,
              indent + 1,
              axom::fmt::format("volume_ratio: {}",
                                format_real(compute_volume_ratio(fit.volume, reference_volume))));
  append_line(
    output,
    indent + 1,
    axom::fmt::format("effective_density: {}",
                      format_real(compute_effective_density(reference_volume, fit.volume))));
}

void append_obb_fit_yaml(std::string& output,
                         int indent,
                         const char* key,
                         const ObbFit& fit,
                         double reference_volume,
                         const char* assumption)
{
  append_line(output, indent, axom::fmt::format("{}:", key));
  if(!fit.valid)
  {
    append_line(output, indent + 1, "available: false");
    append_line(output, indent + 1, "reason: invalid_second_moments");
    return;
  }

  const auto& centroid = fit.box.getCentroid();
  const auto& extents = fit.box.getExtents();
  const auto* axes = fit.box.getAxes();

  append_line(output, indent + 1, "available: true");
  if(assumption != nullptr)
  {
    append_line(output, indent + 1, axom::fmt::format("assumption: '{}'", assumption));
  }
  append_line(output, indent + 1, axom::fmt::format("center: {}", format_point_inline(centroid)));
  append_line(
    output,
    indent + 1,
    axom::fmt::format("extents: {}", format_triplet_inline(extents[0], extents[1], extents[2])));
  append_axes_yaml(output, indent + 1, axes);
  append_line(output, indent + 1, axom::fmt::format("geometric_volume: {}", format_real(fit.volume)));
  append_line(output,
              indent + 1,
              axom::fmt::format("volume_ratio: {}",
                                format_real(compute_volume_ratio(fit.volume, reference_volume))));
  append_line(
    output,
    indent + 1,
    axom::fmt::format("effective_density: {}",
                      format_real(compute_effective_density(reference_volume, fit.volume))));
}

void append_surface_yaml(std::string& output,
                         int indent,
                         const MomentSet& moments,
                         const MassProperties& props)
{
  append_line(output, indent, "surface:");
  append_moment_entries_yaml(output, indent + 1, moments);
  append_mass_properties_yaml(output, indent + 1, props);
}

void append_volume_yaml(std::string& output,
                        int indent,
                        const MomentSet& moments,
                        const MassProperties& props,
                        const PrincipalFrame& frame,
                        const EllipsoidFit& inertia_matched_ellipsoid,
                        const EllipsoidFit& same_volume_ellipsoid,
                        const ObbFit& inertia_matched_obb,
                        const ObbFit& same_volume_obb,
                        const std::string& inertia_matched_ellipsoid_vtk_file,
                        const std::string& same_volume_ellipsoid_vtk_file)
{
  append_line(output, indent, "volume:");
  append_line(output, indent + 1, "assumptions:");
  append_line(
    output,
    indent + 2,
    "- 'volume properties assume the STEP model is a closed, consistently oriented boundary'");
  append_line(output,
              indent + 2,
              "- 'inertia_matched fits reproduce the 0th, 1st, and 2nd moments but may not match "
              "the geometric volume at unit density'");
  append_line(output,
              indent + 2,
              "- 'same_volume_scaled fits uniformly scale the inertia_matched shape to match the "
              "geometric volume while preserving principal directions and aspect ratios'");
  append_line(output,
              indent + 2,
              "- 'oriented_boxes reported here are moment-based proxies and are not guaranteed to "
              "bound the model'");

  append_moment_entries_yaml(output, indent + 1, moments);
  append_mass_properties_yaml(output, indent + 1, props);
  append_principal_frame_yaml(output, indent + 1, frame);

  append_line(output, indent + 1, "fit_proxies:");
  append_line(output, indent + 2, "ellipsoids:");
  append_ellipsoid_fit_yaml(output,
                            indent + 3,
                            "inertia_matched",
                            inertia_matched_ellipsoid,
                            props.measure,
                            nullptr,
                            inertia_matched_ellipsoid_vtk_file);
  append_ellipsoid_fit_yaml(output,
                            indent + 3,
                            "same_volume_scaled",
                            same_volume_ellipsoid,
                            props.measure,
                            "uniformly scaled from the inertia_matched ellipsoid",
                            same_volume_ellipsoid_vtk_file);

  append_line(output, indent + 2, "oriented_boxes:");
  append_obb_fit_yaml(output, indent + 3, "inertia_matched", inertia_matched_obb, props.measure, nullptr);
  append_obb_fit_yaml(output,
                      indent + 3,
                      "same_volume_scaled",
                      same_volume_obb,
                      props.measure,
                      "uniformly scaled from the inertia_matched oriented box");
}

std::string build_results_yaml(const std::string& input_file,
                               int patch_count,
                               const std::string& units,
                               int requested_order,
                               int computed_order,
                               int quadrature_order,
                               IntegralMode integral_mode,
                               bool has_surface,
                               const MomentSet& surface_moments,
                               const MassProperties& surface_props,
                               bool has_volume,
                               const MomentSet& volume_moments,
                               const MassProperties& volume_props,
                               const PrincipalFrame& volume_frame,
                               const EllipsoidFit& inertia_matched_ellipsoid,
                               const EllipsoidFit& same_volume_ellipsoid,
                               const ObbFit& inertia_matched_obb,
                               const ObbFit& same_volume_obb,
                               const std::string& inertia_matched_ellipsoid_vtk_file,
                               const std::string& same_volume_ellipsoid_vtk_file)
{
  std::string output;
  append_line(output, 0, "results:");

  append_line(output, 1, "model:");
  append_line(output, 2, axom::fmt::format("file: '{}'", input_file));
  append_line(output, 2, axom::fmt::format("patches: {}", patch_count));
  append_line(output, 2, axom::fmt::format("units: '{}'", units));

  append_line(output, 1, "config:");
  append_line(output, 2, axom::fmt::format("requested_order: {}", requested_order));
  append_line(output, 2, axom::fmt::format("computed_order: {}", computed_order));
  append_line(output, 2, axom::fmt::format("quadrature_order: {}", quadrature_order));
  append_line(output, 2, axom::fmt::format("integral: '{}'", integral_mode_name(integral_mode)));

  append_line(output, 1, "summary:");
  if(has_surface)
  {
    append_line(output,
                2,
                axom::fmt::format("surface_measure: {}", format_real(surface_props.measure)));
  }
  if(has_volume)
  {
    append_line(output, 2, axom::fmt::format("volume_measure: {}", format_real(volume_props.measure)));
  }

  if(has_surface)
  {
    append_surface_yaml(output, 1, surface_moments, surface_props);
  }

  if(has_volume)
  {
    append_volume_yaml(output,
                       1,
                       volume_moments,
                       volume_props,
                       volume_frame,
                       inertia_matched_ellipsoid,
                       same_volume_ellipsoid,
                       inertia_matched_obb,
                       same_volume_obb,
                       inertia_matched_ellipsoid_vtk_file,
                       same_volume_ellipsoid_vtk_file);
  }

  return output;
}

}  // namespace

int main(int argc, char** argv)
{
  std::string input_file;
  int max_degree {1};
  int quadrature_order {0};
  bool verbose {false};
  bool validate_model {false};
  std::string ellipsoid_vtk_prefix;
  std::string annotationMode {"none"};
  IntegralMode integral_mode {IntegralMode::BOTH};

  axom::CLI::App app {
    "Load a STEP model and compute geometric moments, centroids, inertia tensors, and volume-based "
    "fit proxies."};
  app.add_option("-f,--file", input_file)
    ->description("Input STEP file")
    ->required()
    ->check(axom::CLI::ExistingFile);
  app.add_option("-n,--order", max_degree)
    ->description("Maximum total degree n for explicit monomial output x^i y^j z^k with i+j+k<=n")
    ->capture_default_str()
    ->check(axom::CLI::NonNegativeNumber);
  app.add_option("--npts", quadrature_order)
    ->description("Quadrature order for each numerical integration stage; defaults to max(8, n+4)")
    ->capture_default_str()
    ->check(axom::CLI::PositiveNumber);
  app.add_option("--integral", integral_mode)
    ->description("Which measures to compute: 'surface', 'volume', or 'both'")
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(s_validIntegralModes));
  app.add_flag("-v,--verbose", verbose, "Enable verbose output")->capture_default_str();
  app.add_flag("--validate", validate_model, "Run STEP model validation checks")->capture_default_str();
  app.add_option("--ellipsoid-vtk-prefix", ellipsoid_vtk_prefix)
    ->description(
      "Write the two volume ellipsoid fits to '<prefix>_inertia_matched_ellipsoid.vtk' and "
      "'<prefix>_same_volume_scaled_ellipsoid.vtk'");
#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "See Axom's Caliper support for additional modes.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  app.get_formatter()->column_width(44);

  try
  {
    app.parse(argc, argv);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  axom::slic::SimpleLogger logger(axom::slic::message::Info);
#ifdef AXOM_USE_CALIPER
  axom::utilities::raii::AnnotationsWrapper annotation_raii_wrapper(annotationMode);
#endif
  AXOM_ANNOTATE_SCOPE("quest step moments example");

  if(quadrature_order <= 0)
  {
    quadrature_order = axom::utilities::max(8, max_degree + 4);
  }

  const int computed_max_degree = axom::utilities::max(max_degree, 2);

  SLIC_INFO(axom::fmt::format("Reading STEP file '{}'", input_file));

  axom::quest::STEPReader reader;
  reader.setFileName(input_file);
  reader.setVerbosity(verbose);

  {
    AXOM_ANNOTATE_SCOPE("read step");
    const int read_status = reader.read(validate_model);
    if(read_status != 0)
    {
      SLIC_ERROR("Failed to read STEP file.");
      return 1;
    }
  }

  const PatchArray& patches = reader.getPatchArray();
  if(patches.empty())
  {
    SLIC_ERROR("STEP file did not contain any patches.");
    return 1;
  }

  if(verbose)
  {
    SLIC_INFO(axom::fmt::format("STEP file units: '{}'", reader.getFileUnits()));
    SLIC_INFO(reader.getBRepStats());
  }

  MomentSet surface_moments;
  MomentSet volume_moments;
  MassProperties surface_props;
  MassProperties volume_props;
  PrincipalFrame volume_frame;
  EllipsoidFit inertia_matched_ellipsoid;
  EllipsoidFit same_volume_ellipsoid;
  ObbFit inertia_matched_obb;
  ObbFit same_volume_obb;
  std::string inertia_matched_ellipsoid_vtk_file;
  std::string same_volume_ellipsoid_vtk_file;

  if(should_compute_surface(integral_mode))
  {
    AXOM_ANNOTATE_SCOPE("compute surface properties");
    surface_moments = compute_moments(max_degree, computed_max_degree, [&](const MomentIndex& idx) {
      auto integrand = [idx](const Point3D& x) -> double { return evaluate_monomial(x, idx); };
      return primal::evaluate_surface_integral(patches, integrand, quadrature_order);
    });
    surface_props = compute_mass_properties(surface_moments);
  }

  if(should_compute_volume(integral_mode))
  {
    AXOM_ANNOTATE_SCOPE("compute volume properties");
    volume_moments = compute_moments(max_degree, computed_max_degree, [&](const MomentIndex& idx) {
      auto integrand = [idx](const Point3D& x) -> double { return evaluate_monomial(x, idx); };
      return primal::evaluate_volume_integral(patches, integrand, quadrature_order);
    });
    volume_props = compute_mass_properties(volume_moments);

    if(volume_props.valid)
    {
      AXOM_ANNOTATE_SCOPE("compute volume fit proxies");
      volume_frame = compute_principal_frame(volume_props);
      inertia_matched_ellipsoid = compute_inertia_matched_ellipsoid(volume_frame);
      same_volume_ellipsoid =
        scale_ellipsoid_to_volume(inertia_matched_ellipsoid, volume_props.measure);
      inertia_matched_obb = compute_inertia_matched_obb(volume_frame);
      same_volume_obb = scale_obb_to_volume(inertia_matched_obb, volume_props.measure);
    }
  }

  if(!ellipsoid_vtk_prefix.empty())
  {
    AXOM_ANNOTATE_SCOPE("write ellipsoid vtk");
    if(!should_compute_volume(integral_mode))
    {
      SLIC_WARNING("Ellipsoid VTK export requested, but volume integrals were not computed.");
    }
    else if(!inertia_matched_ellipsoid.valid || !same_volume_ellipsoid.valid)
    {
      SLIC_WARNING(
        "Ellipsoid VTK export requested, but the volume ellipsoid fits are unavailable.");
    }
    else
    {
      inertia_matched_ellipsoid_vtk_file =
        make_ellipsoid_vtk_filename(ellipsoid_vtk_prefix, "inertia_matched_ellipsoid");
      same_volume_ellipsoid_vtk_file =
        make_ellipsoid_vtk_filename(ellipsoid_vtk_prefix, "same_volume_scaled_ellipsoid");

      const bool wrote_inertia =
        write_ellipsoid_vtk(inertia_matched_ellipsoid, inertia_matched_ellipsoid_vtk_file);
      const bool wrote_same_volume =
        write_ellipsoid_vtk(same_volume_ellipsoid, same_volume_ellipsoid_vtk_file);

      if(!wrote_inertia || !wrote_same_volume)
      {
        SLIC_WARNING("Failed to write one or more ellipsoid VTK files.");
      }
    }
  }

  {
    AXOM_ANNOTATE_SCOPE("log results");
    std::string summary = "SUMMARY";
    if(should_compute_surface(integral_mode))
    {
      summary += axom::fmt::format(" surface_m0={:.16e}", surface_props.measure);
    }
    if(should_compute_volume(integral_mode))
    {
      summary += axom::fmt::format(" volume_m0={:.16e}", volume_props.measure);
    }
    SLIC_INFO(summary);

    SLIC_INFO(build_results_yaml(input_file,
                                 static_cast<int>(patches.size()),
                                 reader.getFileUnits(),
                                 max_degree,
                                 computed_max_degree,
                                 quadrature_order,
                                 integral_mode,
                                 should_compute_surface(integral_mode),
                                 surface_moments,
                                 surface_props,
                                 should_compute_volume(integral_mode),
                                 volume_moments,
                                 volume_props,
                                 volume_frame,
                                 inertia_matched_ellipsoid,
                                 same_volume_ellipsoid,
                                 inertia_matched_obb,
                                 same_volume_obb,
                                 inertia_matched_ellipsoid_vtk_file,
                                 same_volume_ellipsoid_vtk_file));
  }

  return 0;
}
