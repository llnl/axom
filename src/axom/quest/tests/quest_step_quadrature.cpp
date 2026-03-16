// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/io/STEPReader.hpp"
#include "axom/primal/operators/evaluate_integral.hpp"
#include "axom/core/numerics/transforms.hpp"

#include "gtest/gtest.h"

#include <cmath>
#include <string>

namespace primal = axom::primal;
namespace quest = axom::quest;

namespace
{
using Point3D = primal::Point<double, 3>;
using Matrix = axom::numerics::Matrix<double>;
using PatchArray = quest::STEPReader::PatchArray;

PatchArray read_step_patches(const std::string& file)
{
  quest::STEPReader reader;
  reader.setFileName(file);
  EXPECT_EQ(reader.read(false), 0);
  return reader.getPatchArray();
}

Matrix embed_linear_transform(const Matrix& linear)
{
  EXPECT_EQ(linear.getNumRows(), 3);
  EXPECT_EQ(linear.getNumColumns(), 3);

  Matrix affine = Matrix::identity(4);
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      affine(i, j) = linear(i, j);
    }
  }

  return affine;
}

void transform_patches(PatchArray& patches, const Matrix& transform)
{
  for(int p = 0; p < patches.size(); ++p)
  {
    auto& control_points = patches[p].getControlPoints();
    for(axom::IndexType ui = 0; ui < control_points.shape()[0]; ++ui)
    {
      for(axom::IndexType vi = 0; vi < control_points.shape()[1]; ++vi)
      {
        control_points(ui, vi) = primal::transform_point(control_points(ui, vi), transform);
      }
    }
  }
}

PatchArray replicate_biquartic_sphere_faces(const PatchArray& base_patches)
{
  EXPECT_EQ(base_patches.size(), 1);

  PatchArray all_patches;
  all_patches.push_back(base_patches[0]);

  const Matrix rot_pos_y = embed_linear_transform(axom::numerics::transforms::xRotation(M_PI / 2.0));
  const Matrix rot_neg_y = embed_linear_transform(axom::numerics::transforms::xRotation(-M_PI / 2.0));
  const Matrix rot_pos_z = embed_linear_transform(axom::numerics::transforms::xRotation(M_PI));
  const Matrix rot_neg_x = embed_linear_transform(axom::numerics::transforms::yRotation(M_PI / 2.0));
  const Matrix rot_pos_x = embed_linear_transform(axom::numerics::transforms::yRotation(-M_PI / 2.0));

  const Matrix rotations[] = {rot_pos_y, rot_neg_y, rot_pos_z, rot_neg_x, rot_pos_x};
  for(const Matrix& rotation : rotations)
  {
    PatchArray rotated(1);
    rotated[0] = base_patches[0];
    transform_patches(rotated, rotation);
    all_patches.push_back(rotated[0]);
  }

  return all_patches;
}

template <typename Lambda>
void expect_matching_surface_values(const PatchArray& patches_a,
                                    const PatchArray& patches_b,
                                    const Lambda& integrand,
                                    int npts,
                                    double expected,
                                    double abs_tol)
{
  const double observed_a = primal::evaluate_surface_integral(patches_a, integrand, npts);
  const double observed_b = primal::evaluate_surface_integral(patches_b, integrand, npts);

  EXPECT_NEAR(observed_a, expected, abs_tol);
  EXPECT_NEAR(observed_b, expected, abs_tol);
  EXPECT_NEAR(observed_a, observed_b, abs_tol);
}

template <typename Lambda>
void expect_matching_volume_values(const PatchArray& patches_a,
                                   const PatchArray& patches_b,
                                   const Lambda& integrand,
                                   int npts,
                                   double expected,
                                   double abs_tol)
{
  const double observed_a = primal::evaluate_volume_integral(patches_a, integrand, npts);
  const double observed_b = primal::evaluate_volume_integral(patches_b, integrand, npts);

  EXPECT_NEAR(observed_a, expected, abs_tol);
  EXPECT_NEAR(observed_b, expected, abs_tol);
  EXPECT_NEAR(observed_a, observed_b, abs_tol);
}

}  // namespace

TEST(quest_step_quadrature, tet_surface_and_volume)
{
  auto patches = read_step_patches(std::string(AXOM_DATA_DIR) + "/quest/step/tet.step");
  EXPECT_EQ(patches.size(), 4);

  auto const_integrand = [](Point3D /*x*/) -> double { return 1.0; };

  constexpr int npts = 6;
  constexpr double abs_tol = 1e-10;

  EXPECT_NEAR(primal::evaluate_surface_integral(patches, const_integrand, npts),
              8.0 * std::sqrt(3.0),
              abs_tol);
  EXPECT_NEAR(primal::evaluate_volume_integral(patches, const_integrand, npts), 8.0 / 3.0, abs_tol);
}

TEST(quest_step_quadrature, sphere_models_match_unit_sphere_moments)
{
  const std::string step_dir = std::string(AXOM_DATA_DIR) + "/quest/step/";

  auto revolved = read_step_patches(step_dir + "revolved_sphere.step");
  auto biquartic = read_step_patches(step_dir + "biquartic_sphere_surface.step");
  biquartic = replicate_biquartic_sphere_faces(biquartic);

  EXPECT_EQ(revolved.size(), 1);
  EXPECT_EQ(biquartic.size(), 6);

  // revolved_sphere.step is a radius-5 sphere, while the biquartic sphere patch is unit-radius.
  transform_patches(revolved, axom::numerics::transforms::scale(1.0 / 5.0, 4));

  constexpr int npts = 16;
  constexpr double abs_tol = 5e-5;

  auto const_integrand = [](Point3D /*x*/) -> double { return 1.0; };
  auto x_integrand = [](Point3D x) -> double { return x[0]; };
  auto y_integrand = [](Point3D x) -> double { return x[1]; };
  auto z_integrand = [](Point3D x) -> double { return x[2]; };
  auto xx_integrand = [](Point3D x) -> double { return x[0] * x[0]; };
  auto yy_integrand = [](Point3D x) -> double { return x[1] * x[1]; };
  auto zz_integrand = [](Point3D x) -> double { return x[2] * x[2]; };

  const double sphere_area = 4.0 * M_PI;
  const double sphere_volume = 4.0 * M_PI / 3.0;
  const double surface_second_moment = 4.0 * M_PI / 3.0;
  const double volume_second_moment = 4.0 * M_PI / 15.0;

  expect_matching_surface_values(revolved, biquartic, const_integrand, npts, sphere_area, abs_tol);
  expect_matching_surface_values(revolved, biquartic, x_integrand, npts, 0.0, abs_tol);
  expect_matching_surface_values(revolved, biquartic, y_integrand, npts, 0.0, abs_tol);
  expect_matching_surface_values(revolved, biquartic, z_integrand, npts, 0.0, abs_tol);
  expect_matching_surface_values(revolved, biquartic, xx_integrand, npts, surface_second_moment, abs_tol);
  expect_matching_surface_values(revolved, biquartic, yy_integrand, npts, surface_second_moment, abs_tol);
  expect_matching_surface_values(revolved, biquartic, zz_integrand, npts, surface_second_moment, abs_tol);

  expect_matching_volume_values(revolved, biquartic, const_integrand, npts, sphere_volume, abs_tol);
  expect_matching_volume_values(revolved, biquartic, x_integrand, npts, 0.0, abs_tol);
  expect_matching_volume_values(revolved, biquartic, y_integrand, npts, 0.0, abs_tol);
  expect_matching_volume_values(revolved, biquartic, z_integrand, npts, 0.0, abs_tol);
  expect_matching_volume_values(revolved, biquartic, xx_integrand, npts, volume_second_moment, abs_tol);
  expect_matching_volume_values(revolved, biquartic, yy_integrand, npts, volume_second_moment, abs_tol);
  expect_matching_volume_values(revolved, biquartic, zz_integrand, npts, volume_second_moment, abs_tol);
}

TEST(quest_step_quadrature, transformed_sphere_models_match_expected_moments)
{
  const std::string step_dir = std::string(AXOM_DATA_DIR) + "/quest/step/";

  auto revolved = read_step_patches(step_dir + "revolved_sphere.step");
  auto biquartic = read_step_patches(step_dir + "biquartic_sphere_surface.step");
  biquartic = replicate_biquartic_sphere_faces(biquartic);

  transform_patches(revolved, axom::numerics::transforms::scale(1.0 / 5.0, 4));

  constexpr double sphere_scale = 1.7;
  const Point3D center {1.25, -0.75, 2.0};

  const Matrix scale = axom::numerics::transforms::scale(sphere_scale, 4);
  const Matrix rotation =
    embed_linear_transform(axom::numerics::transforms::axisRotation(M_PI / 4.0, 1.0, 2.0, 3.0));
  const Matrix translation = axom::numerics::transforms::translate(center[0], center[1], center[2]);

  transform_patches(revolved, scale);
  transform_patches(revolved, rotation);
  transform_patches(revolved, translation);

  transform_patches(biquartic, scale);
  transform_patches(biquartic, rotation);
  transform_patches(biquartic, translation);

  constexpr int npts = 16;
  constexpr double abs_tol = 1e-4;

  auto const_integrand = [](Point3D /*x*/) -> double { return 1.0; };
  auto x_integrand = [](Point3D x) -> double { return x[0]; };
  auto y_integrand = [](Point3D x) -> double { return x[1]; };
  auto z_integrand = [](Point3D x) -> double { return x[2]; };
  auto centered_xx_integrand = [center](Point3D x) -> double {
    const double dx = x[0] - center[0];
    return dx * dx;
  };
  auto centered_yy_integrand = [center](Point3D x) -> double {
    const double dy = x[1] - center[1];
    return dy * dy;
  };
  auto centered_zz_integrand = [center](Point3D x) -> double {
    const double dz = x[2] - center[2];
    return dz * dz;
  };

  const double sphere_area = 4.0 * M_PI * sphere_scale * sphere_scale;
  const double sphere_volume = (4.0 * M_PI / 3.0) * sphere_scale * sphere_scale * sphere_scale;
  const double surface_second_moment = (4.0 * M_PI / 3.0) * std::pow(sphere_scale, 4);
  const double volume_second_moment = (4.0 * M_PI / 15.0) * std::pow(sphere_scale, 5);

  expect_matching_surface_values(revolved, biquartic, const_integrand, npts, sphere_area, abs_tol);
  expect_matching_surface_values(revolved, biquartic, x_integrand, npts, sphere_area * center[0], abs_tol);
  expect_matching_surface_values(revolved, biquartic, y_integrand, npts, sphere_area * center[1], abs_tol);
  expect_matching_surface_values(revolved, biquartic, z_integrand, npts, sphere_area * center[2], abs_tol);
  expect_matching_surface_values(revolved,
                                 biquartic,
                                 centered_xx_integrand,
                                 npts,
                                 surface_second_moment,
                                 abs_tol);
  expect_matching_surface_values(revolved,
                                 biquartic,
                                 centered_yy_integrand,
                                 npts,
                                 surface_second_moment,
                                 abs_tol);
  expect_matching_surface_values(revolved,
                                 biquartic,
                                 centered_zz_integrand,
                                 npts,
                                 surface_second_moment,
                                 abs_tol);

  expect_matching_volume_values(revolved, biquartic, const_integrand, npts, sphere_volume, abs_tol);
  expect_matching_volume_values(revolved, biquartic, x_integrand, npts, sphere_volume * center[0], abs_tol);
  expect_matching_volume_values(revolved, biquartic, y_integrand, npts, sphere_volume * center[1], abs_tol);
  expect_matching_volume_values(revolved, biquartic, z_integrand, npts, sphere_volume * center[2], abs_tol);
  expect_matching_volume_values(revolved,
                                biquartic,
                                centered_xx_integrand,
                                npts,
                                volume_second_moment,
                                abs_tol);
  expect_matching_volume_values(revolved,
                                biquartic,
                                centered_yy_integrand,
                                npts,
                                volume_second_moment,
                                abs_tol);
  expect_matching_volume_values(revolved,
                                biquartic,
                                centered_zz_integrand,
                                npts,
                                volume_second_moment,
                                abs_tol);
}
