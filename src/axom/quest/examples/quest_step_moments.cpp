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

#ifdef AXOM_USE_CONDUIT
  #include "conduit.hpp"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <variant>
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

double clamp_near_zero(double value, double tolerance)
{
  return std::abs(value) <= tolerance ? 0.0 : value;
}

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

  axom::IndexType next_node_id = 0;

  auto append_node = [&mesh, &fit, &next_node_id](double x, double y, double z) {
    // Map a point on the unit sphere into the ellipsoid's principal-axis frame.
    Point3D point = fit.centroid;
    const double local_coords[3] {x, y, z};
    for(int axis = 0; axis < 3; ++axis)
    {
      for(int dim = 0; dim < 3; ++dim)
      {
        point[dim] += fit.radii[axis] * local_coords[axis] * fit.axes[axis][dim];
      }
    }

    mesh.appendNode(point[0], point[1], point[2]);
    return next_node_id++;
  };

  auto append_ring = [&append_node, theta_resolution](std::vector<axom::IndexType>& ring, double phi) {
    ring.clear();

    const double sin_phi = std::sin(phi);
    const double cos_phi = std::cos(phi);
    if(std::abs(sin_phi) <= eps)
    {
      ring.push_back(append_node(0.0, 0.0, cos_phi));
      return;
    }

    ring.reserve(theta_resolution);
    for(int i = 0; i < theta_resolution; ++i)
    {
      const double theta = 2.0 * M_PI * static_cast<double>(i) / theta_resolution;
      ring.push_back(append_node(std::cos(theta) * sin_phi, std::sin(theta) * sin_phi, cos_phi));
    }
  };

  auto append_ring_triangles = [&mesh](const std::vector<axom::IndexType>& lower_ring,
                                       const std::vector<axom::IndexType>& upper_ring) {
    SLIC_ASSERT(!lower_ring.empty());
    SLIC_ASSERT(!upper_ring.empty());

    axom::IndexType cell[3];
    if(lower_ring.size() == 1)
    {
      const axom::IndexType pole = lower_ring.front();
      const int ring_size = static_cast<int>(upper_ring.size());
      for(int i = 0; i < ring_size; ++i)
      {
        const int next_i = (i + 1) % ring_size;
        cell[0] = pole;
        cell[1] = upper_ring[next_i];
        cell[2] = upper_ring[i];
        mesh.appendCell(cell);
      }
      return;
    }

    if(upper_ring.size() == 1)
    {
      const axom::IndexType pole = upper_ring.front();
      const int ring_size = static_cast<int>(lower_ring.size());
      for(int i = 0; i < ring_size; ++i)
      {
        const int next_i = (i + 1) % ring_size;
        cell[0] = pole;
        cell[1] = lower_ring[i];
        cell[2] = lower_ring[next_i];
        mesh.appendCell(cell);
      }
      return;
    }

    SLIC_ASSERT(lower_ring.size() == upper_ring.size());
    const int ring_size = static_cast<int>(lower_ring.size());
    for(int i = 0; i < ring_size; ++i)
    {
      const int next_i = (i + 1) % ring_size;
      cell[0] = lower_ring[i];
      cell[1] = upper_ring[i];
      cell[2] = lower_ring[next_i];
      mesh.appendCell(cell);

      cell[0] = lower_ring[next_i];
      cell[1] = upper_ring[i];
      cell[2] = upper_ring[next_i];
      mesh.appendCell(cell);
    }
  };

  std::vector<axom::IndexType> lower_ring;
  std::vector<axom::IndexType> upper_ring;
  lower_ring.reserve(theta_resolution);
  upper_ring.reserve(theta_resolution);

  append_ring(lower_ring, 0.0);
  for(int j = 1; j <= num_rings; ++j)
  {
    const double phi = M_PI * static_cast<double>(j) / (phi_resolution - 1);
    append_ring(upper_ring, phi);
    append_ring_triangles(lower_ring, upper_ring);
    lower_ring.swap(upper_ring);
  }

  append_ring(upper_ring, M_PI);
  append_ring_triangles(lower_ring, upper_ring);

  return mint::write_vtk(&mesh, file_path) == 0;
}

std::string join_path(const std::string& prefix, const std::string& key)
{
  return prefix.empty() ? key : axom::fmt::format("{}/{}", prefix, key);
}

class ResultsStore
{
public:
  explicit ResultsStore(double zero_tolerance = 0.0) : m_zero_tolerance(zero_tolerance) { }

  void set(const std::string& path, const std::array<double, 3>& value) { set_vec3(path, value); }

  void add(const std::string& prefix, const std::string& key, const std::string& value)
  {
    set_string(join_path(prefix, key), value);
  }

  void add(const std::string& prefix, const std::string& key, const char* value)
  {
    set_string(join_path(prefix, key), value);
  }

  void add(const std::string& prefix, const std::string& key, double value)
  {
    set_real(join_path(prefix, key), value);
  }

  void add(const std::string& prefix, const std::string& key, const std::array<double, 3>& value)
  {
    set_vec3(join_path(prefix, key), value);
  }

  void add_unclamped_real(const std::string& prefix, const std::string& key, double value)
  {
    set_real_unclamped(join_path(prefix, key), value);
  }

  void add(const std::string& prefix, const std::string& key, int value)
  {
    set_integer(join_path(prefix, key), value);
  }

  void add(const std::string& prefix, const std::string& key, bool value)
  {
    set_boolean(join_path(prefix, key), value);
  }

  std::string to_yaml() const
  {
#ifdef AXOM_USE_CONDUIT
    return m_root.to_yaml();
#else
    std::string output;
    append_yaml(output, m_root, 0);
    return output;
#endif
  }

private:
  double sanitize_real(double value) const { return clamp_near_zero(value, m_zero_tolerance); }

  const double m_zero_tolerance;

#ifdef AXOM_USE_CONDUIT
  conduit::Node m_root;

  void set_string(const std::string& path, const std::string& value) { m_root[path] = value; }

  void set_string(const std::string& path, const char* value) { m_root[path].set_string(value); }

  void set_real(const std::string& path, double value) { m_root[path] = sanitize_real(value); }

  void set_vec3(const std::string& path, const std::array<double, 3>& value)
  {
    const std::array<double, 3> sanitized {sanitize_real(value[0]),
                                           sanitize_real(value[1]),
                                           sanitize_real(value[2])};
    m_root[path].set(sanitized.data(), 3);
  }

  void set_real_unclamped(const std::string& path, double value) { m_root[path] = value; }

  void set_integer(const std::string& path, int value) { m_root[path] = value; }

  void set_boolean(const std::string& path, bool value)
  {
    m_root[path].set_string(value ? "true" : "false");
  }
#else
  using ScalarValue = std::variant<std::string, int, double, bool, std::array<double, 3>>;

  struct TreeNode
  {
    bool has_scalar {false};
    ScalarValue scalar_value {std::string {}};
    std::unordered_map<std::string, TreeNode> children;
    std::vector<std::string> child_order;

    TreeNode& fetch_child(const std::string& key)
    {
      auto it = children.find(key);
      if(it == children.end())
      {
        child_order.push_back(key);
        it = children.emplace(key, TreeNode {}).first;
      }

      return it->second;
    }
  };

  TreeNode m_root;

  TreeNode& fetch_path(const std::string& path)
  {
    TreeNode* node = &m_root;

    for(const auto& key : axom::utilities::string::split(path, '/'))
    {
      if(!key.empty())
      {
        node = &node->fetch_child(key);
      }
    }

    return *node;
  }

  void set_string(const std::string& path, const std::string& value)
  {
    TreeNode& node = fetch_path(path);
    node.has_scalar = true;
    node.scalar_value = value;
  }

  void set_string(const std::string& path, const char* value)
  {
    set_string(path, std::string(value));
  }

  void set_real(const std::string& path, double value)
  {
    TreeNode& node = fetch_path(path);
    node.has_scalar = true;
    node.scalar_value = sanitize_real(value);
  }

  void set_vec3(const std::string& path, const std::array<double, 3>& value)
  {
    TreeNode& node = fetch_path(path);
    node.has_scalar = true;
    node.scalar_value = std::array<double, 3> {sanitize_real(value[0]),
                                               sanitize_real(value[1]),
                                               sanitize_real(value[2])};
  }

  void set_real_unclamped(const std::string& path, double value)
  {
    TreeNode& node = fetch_path(path);
    node.has_scalar = true;
    node.scalar_value = value;
  }

  void set_integer(const std::string& path, int value)
  {
    TreeNode& node = fetch_path(path);
    node.has_scalar = true;
    node.scalar_value = value;
  }

  void set_boolean(const std::string& path, bool value)
  {
    TreeNode& node = fetch_path(path);
    node.has_scalar = true;
    node.scalar_value = value;
  }

  static std::string format_scalar(const ScalarValue& value)
  {
    if(const auto* string_value = std::get_if<std::string>(&value))
    {
      std::string escaped;
      escaped.reserve(string_value->size());
      for(char ch : *string_value)
      {
        escaped += ch;
        if(ch == '\'')
        {
          escaped += '\'';
        }
      }

      return axom::fmt::format("'{}'", escaped);
    }
    if(const auto* integer_value = std::get_if<int>(&value))
    {
      return axom::fmt::format("{}", *integer_value);
    }
    if(const auto* real_value = std::get_if<double>(&value))
    {
      return axom::fmt::format("{:.16e}", *real_value);
    }
    if(const auto* vec3_value = std::get_if<std::array<double, 3>>(&value))
    {
      return axom::fmt::format("[{:.16e}, {:.16e}, {:.16e}]",
                               (*vec3_value)[0],
                               (*vec3_value)[1],
                               (*vec3_value)[2]);
    }

    return std::get<bool>(value) ? "true" : "false";
  }

  static void append_yaml(std::string& output, const TreeNode& node, int indent)
  {
    const auto append_line = [&output](int line_indent,
                                       const std::string& key,
                                       const std::string* value = nullptr) {
      axom::fmt::format_to(std::back_inserter(output), "{}{}:", std::string(line_indent * 2, ' '), key);
      if(value != nullptr)
      {
        axom::fmt::format_to(std::back_inserter(output), " {}", *value);
      }
      output += '\n';
    };

    for(const auto& key : node.child_order)
    {
      const auto it = node.children.find(key);
      SLIC_ASSERT(it != node.children.end());
      const TreeNode& child = it->second;

      if(child.children.empty())
      {
        if(child.has_scalar)
        {
          const std::string value = format_scalar(child.scalar_value);
          append_line(indent, key, &value);
        }
        else
        {
          append_line(indent, key);
        }
        continue;
      }

      append_line(indent, key);
      if(child.has_scalar)
      {
        const std::string value = format_scalar(child.scalar_value);
        append_line(indent + 1, "value", &value);
      }
      append_yaml(output, child, indent + 1);
    }
  }
#endif
};

class ResultsNode
{
public:
  ResultsNode(ResultsStore& results, std::string prefix)
    : m_results(results)
    , m_prefix(std::move(prefix))
  { }

  ResultsNode child(const std::string& key) const
  {
    return ResultsNode {m_results, join_path(m_prefix, key)};
  }

  void add(const std::string& key, const std::string& value) const
  {
    m_results.add(m_prefix, key, value);
  }

  void add(const std::string& key, const char* value) const { m_results.add(m_prefix, key, value); }

  void add(const std::string& key, double value) const { m_results.add(m_prefix, key, value); }

  void add(const std::string& key, const std::array<double, 3>& value) const
  {
    m_results.add(m_prefix, key, value);
  }

  void set(const std::array<double, 3>& value) const { m_results.set(m_prefix, value); }

  void add(const std::string& key, int value) const { m_results.add(m_prefix, key, value); }

  void add(const std::string& key, bool value) const { m_results.add(m_prefix, key, value); }

private:
  ResultsStore& m_results;
  std::string m_prefix;
};

void write_xyz_triplet(const ResultsNode& node, double x, double y, double z)
{
  node.set(std::array<double, 3> {x, y, z});
}

void write(const ResultsNode& node, const Point3D& point)
{
  write_xyz_triplet(node, point[0], point[1], point[2]);
}

void write(const ResultsNode& node, const Vector3D& vector)
{
  write_xyz_triplet(node, vector[0], vector[1], vector[2]);
}

void write_axis_scalars(const ResultsNode& node, double a, double b, double c)
{
  node.add("axis_0", a);
  node.add("axis_1", b);
  node.add("axis_2", c);
}

void write_axes(const ResultsNode& node, const Vector3D* axes)
{
  for(int i = 0; i < 3; ++i)
  {
    write(node.child(axom::fmt::format("axis_{}", i)), axes[i]);
  }
}

void write(const ResultsNode& node, const InertiaTensor& tensor)
{
  node.add("xx", tensor.xx);
  node.add("yy", tensor.yy);
  node.add("zz", tensor.zz);
  node.add("xy", tensor.xy);
  node.add("xz", tensor.xz);
  node.add("yz", tensor.yz);
}

void write(const ResultsNode& node, const MomentSet& moments)
{
  node.add("measure", get_moment(moments, 0, 0, 0));

  const ResultsNode entries = node.child("entries");
  entries.add("count", static_cast<int>(moments.requested_entries.size()));
  for(const auto& entry : moments.requested_entries)
  {
    entries.add(axom::fmt::format("m_{}_{}_{}", entry.index.px, entry.index.py, entry.index.pz),
                entry.value);
  }
}

void write(const ResultsNode& node, const MassProperties& props)
{
  if(!props.valid)
  {
    node.add("available", props.valid);
    node.add("reason", "measure_is_zero");
    return;
  }

  write(node.child("centroid"), props.centroid);
  write(node.child("inertia_origin"), props.inertia_origin);
  write(node.child("inertia_centroid"), props.inertia_centroid);
}

void write(const ResultsNode& node, const PrincipalFrame& frame)
{
  if(!frame.valid)
  {
    node.add("available", frame.valid);
    node.add("reason", "eigensolve_failed_or_measure_is_zero");
    return;
  }

  write_axis_scalars(node.child("principal_inertia"),
                     frame.principal_inertia[0],
                     frame.principal_inertia[1],
                     frame.principal_inertia[2]);
  write_axis_scalars(node.child("principal_second_moments"),
                     frame.principal_second_moments[0],
                     frame.principal_second_moments[1],
                     frame.principal_second_moments[2]);
  write_axes(node.child("axes"), &frame.axes[0]);
}

void write_assumptions(const ResultsNode& node, std::initializer_list<const char*> assumptions)
{
  int idx = 0;
  for(const char* assumption : assumptions)
  {
    node.add(axom::fmt::format("assumption_{}", idx++), assumption);
  }
}

void write_ellipsoid_fit(const ResultsNode& node,
                         const EllipsoidFit& fit,
                         double reference_volume,
                         const char* assumption,
                         const std::string& vtk_file)
{
  if(!fit.valid)
  {
    node.add("available", fit.valid);
    node.add("reason", "invalid_second_moments");
    return;
  }

  if(assumption != nullptr)
  {
    node.add("assumption", assumption);
  }
  if(!vtk_file.empty())
  {
    node.add("vtk_file", vtk_file);
  }

  write(node.child("center"), fit.centroid);
  write_axis_scalars(node.child("semiaxes"), fit.radii[0], fit.radii[1], fit.radii[2]);
  write_axes(node.child("axes"), &fit.axes[0]);
  node.add("geometric_volume", fit.volume);
  node.add("volume_ratio", compute_volume_ratio(fit.volume, reference_volume));
  node.add("effective_density", compute_effective_density(reference_volume, fit.volume));
}

void write_obb_fit(const ResultsNode& node,
                   const ObbFit& fit,
                   double reference_volume,
                   const char* assumption)
{
  if(!fit.valid)
  {
    node.add("available", fit.valid);
    node.add("reason", "invalid_second_moments");
    return;
  }

  if(assumption != nullptr)
  {
    node.add("assumption", assumption);
  }

  const auto& centroid = fit.box.getCentroid();
  const auto& extents = fit.box.getExtents();
  const auto* axes = fit.box.getAxes();
  write(node.child("center"), centroid);
  write_xyz_triplet(node.child("extents"), extents[0], extents[1], extents[2]);
  write_axes(node.child("axes"), axes);
  node.add("geometric_volume", fit.volume);
  node.add("volume_ratio", compute_volume_ratio(fit.volume, reference_volume));
  node.add("effective_density", compute_effective_density(reference_volume, fit.volume));
}

void write_results(ResultsStore& results,
                   const std::string& input_file,
                   int patch_count,
                   const std::string& units,
                   int requested_order,
                   int computed_order,
                   int quadrature_order,
                   double zero_tolerance,
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
  const ResultsNode root {results, "results"};

  const ResultsNode model = root.child("model");
  model.add("file", input_file);
  model.add("patches", patch_count);
  model.add("units", units);

  const ResultsNode config = root.child("config");
  config.add("requested_order", requested_order);
  config.add("computed_order", computed_order);
  config.add("quadrature_order", quadrature_order);
  results.add_unclamped_real("results/config", "zero_tolerance", zero_tolerance);
  config.add("integral", integral_mode_name(integral_mode));

  const ResultsNode summary = root.child("summary");
  if(has_surface)
  {
    summary.add("surface_measure", surface_props.measure);
  }
  if(has_volume)
  {
    summary.add("volume_measure", volume_props.measure);
  }

  if(has_surface)
  {
    const ResultsNode surface = root.child("surface");
    write(surface.child("raw_moments"), surface_moments);
    write(surface.child("derived"), surface_props);
  }

  if(has_volume)
  {
    const ResultsNode volume = root.child("volume");
    write_assumptions(
      volume.child("assumptions"),
      {"volume properties assume the STEP model is a closed, consistently oriented boundary",
       "inertia_matched fits reproduce the 0th, 1st, and 2nd moments but may not match the "
       "geometric volume at unit density",
       "same_volume_scaled fits uniformly scale the inertia_matched shape to match the geometric "
       "volume while preserving principal directions and aspect ratios",
       "oriented_boxes reported here are moment-based proxies and are not guaranteed to bound the "
       "model"});
    write(volume.child("raw_moments"), volume_moments);
    write(volume.child("derived"), volume_props);
    write(volume.child("principal_frame"), volume_frame);

    const ResultsNode fit_proxies = volume.child("fit_proxies");
    const ResultsNode ellipsoids = fit_proxies.child("ellipsoids");
    write_ellipsoid_fit(ellipsoids.child("inertia_matched"),
                        inertia_matched_ellipsoid,
                        volume_props.measure,
                        nullptr,
                        inertia_matched_ellipsoid_vtk_file);
    write_ellipsoid_fit(ellipsoids.child("same_volume_scaled"),
                        same_volume_ellipsoid,
                        volume_props.measure,
                        "uniformly scaled from the inertia_matched ellipsoid",
                        same_volume_ellipsoid_vtk_file);

    const ResultsNode oriented_boxes = fit_proxies.child("oriented_boxes");
    write_obb_fit(oriented_boxes.child("inertia_matched"),
                  inertia_matched_obb,
                  volume_props.measure,
                  nullptr);
    write_obb_fit(oriented_boxes.child("same_volume_scaled"),
                  same_volume_obb,
                  volume_props.measure,
                  "uniformly scaled from the inertia_matched oriented box");
  }
}
}  // namespace

int main(int argc, char** argv)
{
  std::string input_file;
  int max_degree {1};
  int quadrature_order {0};
  bool verbose {false};
  bool validate_model {false};
  double zero_tolerance {0.0};
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
  app.add_option("--zero-tolerance", zero_tolerance)
    ->description("Clamp reported floating-point values with abs(value) <= tol to exactly zero")
    ->capture_default_str()
    ->check(axom::CLI::NonNegativeNumber);
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
      summary += axom::fmt::format(" surface_m0={:.16e}",
                                   clamp_near_zero(surface_props.measure, zero_tolerance));
    }
    if(should_compute_volume(integral_mode))
    {
      summary += axom::fmt::format(" volume_m0={:.16e}",
                                   clamp_near_zero(volume_props.measure, zero_tolerance));
    }
    SLIC_INFO(summary);

    ResultsStore results(zero_tolerance);
    write_results(results,
                  input_file,
                  static_cast<int>(patches.size()),
                  reader.getFileUnits(),
                  max_degree,
                  computed_max_degree,
                  quadrature_order,
                  zero_tolerance,
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
                  same_volume_ellipsoid_vtk_file);
    SLIC_INFO(results.to_yaml());
  }

  return 0;
}
