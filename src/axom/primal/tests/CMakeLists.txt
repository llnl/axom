# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal unit tests
#------------------------------------------------------------------------------

set( primal_tests
    primal_bezier_curve.cpp
    primal_bezier_intersect.cpp
    primal_bezier_patch.cpp
    primal_boundingbox.cpp
    primal_bounding_box_intersect.cpp
    primal_clip.cpp
    primal_closest_point.cpp
    primal_compute_bounding_box.cpp
    primal_compute_moments.cpp
    primal_cone.cpp
    primal_coord_transform.cpp
    primal_curved_polygon.cpp
    primal_hexahedron.cpp
    primal_in_sphere.cpp
    primal_integral.cpp
    primal_intersect.cpp
    primal_intersect_impl.cpp
    primal_knot_vector.cpp
    primal_nurbs_curve.cpp
    primal_nurbs_patch.cpp
    primal_orientation.cpp
    primal_orientedboundingbox.cpp
    primal_plane.cpp
    primal_point.cpp
    primal_polygon.cpp
    primal_polyhedron.cpp
    primal_rational_bezier.cpp
    primal_ray_intersect.cpp
    primal_segment.cpp
    primal_sphere.cpp
    primal_solid_angle.cpp
    primal_split.cpp
    primal_squared_distance.cpp
    primal_surface_intersect.cpp
    primal_tetrahedron.cpp
    primal_trimming_curves.cpp
    primal_octahedron.cpp
    primal_triangle.cpp
    primal_quadrilateral.cpp
    primal_vector.cpp
    primal_winding_number.cpp
    primal_zip.cpp

    primal_clip_perf.cpp
   )

set(primal_test_depends primal fmt gtest)

foreach ( test ${primal_tests} )

   get_filename_component( test_name ${test} NAME_WE )

   axom_add_executable(
        NAME       ${test_name}_test
        SOURCES    ${test}
        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
        DEPENDS_ON ${primal_test_depends}
        FOLDER     axom/primal/tests
       )

   axom_add_test(
        NAME    ${test_name}
        COMMAND ${test_name}_test
       )

endforeach()
