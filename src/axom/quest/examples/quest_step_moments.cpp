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
 * The example loads a STEP BRep with Quest's STEP reader, then evaluates
 * monomial moments \f$\int x^i y^j z^k \, dS\f$ and/or
 * \f$\int x^i y^j z^k \, dV\f$ for all nonnegative exponent triples with
 * total degree \f$i + j + k \le n\f$, where \a n is user supplied.
 * It also derives centroids and inertia tensors from the first- and second-
 * order moments.
 *
 * \note Volume moments are only geometrically meaningful when the STEP model
 * represents a closed, consistently oriented boundary.
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/core/utilities/CommandLineUtilities.hpp"
#include "axom/primal.hpp"
#include "axom/primal/operators/evaluate_integral.hpp"
#include "axom/quest.hpp"
#include "axom/slic.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace primal = axom::primal;
namespace slic = axom::slic;

namespace
{

using Point3D = primal::Point<double, 3>;
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
  constexpr double eps = 1e-14;

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

void log_moment_entries(const char* label, const MomentSet& moments)
{
  for(const auto& entry : moments.requested_entries)
  {
    SLIC_INFO(axom::fmt::format("{}_MOMENT degree={} exponents=({},{},{}) value={:.16e}",
                                label,
                                entry.index.degree,
                                entry.index.px,
                                entry.index.py,
                                entry.index.pz,
                                entry.value));
  }
}

void log_mass_properties(const char* label, const MassProperties& props)
{
  if(!props.valid)
  {
    SLIC_INFO(axom::fmt::format("{}_PROPERTIES unavailable=measure_is_zero", label));
    return;
  }

  SLIC_INFO(axom::fmt::format("{}_CENTROID x={:.16e} y={:.16e} z={:.16e}",
                              label,
                              props.centroid[0],
                              props.centroid[1],
                              props.centroid[2]));
  SLIC_INFO(axom::fmt::format(
    "{}_INERTIA_ORIGIN xx={:.16e} yy={:.16e} zz={:.16e} xy={:.16e} xz={:.16e} yz={:.16e}",
    label,
    props.inertia_origin.xx,
    props.inertia_origin.yy,
    props.inertia_origin.zz,
    props.inertia_origin.xy,
    props.inertia_origin.xz,
    props.inertia_origin.yz));
  SLIC_INFO(axom::fmt::format(
    "{}_INERTIA_CENTROID xx={:.16e} yy={:.16e} zz={:.16e} xy={:.16e} xz={:.16e} yz={:.16e}",
    label,
    props.inertia_centroid.xx,
    props.inertia_centroid.yy,
    props.inertia_centroid.zz,
    props.inertia_centroid.xy,
    props.inertia_centroid.xz,
    props.inertia_centroid.yz));
}

}  // namespace

int main(int argc, char** argv)
{
  std::string input_file;
  int max_degree {1};
  int quadrature_order {0};
  bool verbose {false};
  bool validate_model {false};
  std::string annotationMode {"none"};
  IntegralMode integral_mode {IntegralMode::BOTH};

  axom::CLI::App app {
    "Load a STEP model and compute geometric moments, centroids, and inertia tensors."};
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

  if(should_compute_volume(integral_mode))
  {
    SLIC_WARNING(
      "Volume properties assume the STEP model is a closed, consistently oriented boundary.");
  }

  SLIC_INFO(axom::fmt::format("MODEL file='{}' patches={} units='{}'",
                              input_file,
                              patches.size(),
                              reader.getFileUnits()));
  SLIC_INFO(axom::fmt::format("CONFIG requested_order={} computed_order={} npts={} integral={}",
                              max_degree,
                              computed_max_degree,
                              quadrature_order,
                              integral_mode_name(integral_mode)));

  MomentSet surface_moments;
  MomentSet volume_moments;
  MassProperties surface_props;
  MassProperties volume_props;

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

    if(should_compute_surface(integral_mode))
    {
      log_mass_properties("SURFACE", surface_props);
      log_moment_entries("SURFACE", surface_moments);
    }

    if(should_compute_volume(integral_mode))
    {
      log_mass_properties("VOLUME", volume_props);
      log_moment_entries("VOLUME", volume_moments);
    }
  }

  return 0;
}
