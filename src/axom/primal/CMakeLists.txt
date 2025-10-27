# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal -- Computational geometry primitives
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       PRIMAL
                        COMPONENTS SLIC)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set( primal_headers

    constants.hpp

    ## geometry
    geometry/BezierCurve.hpp
    geometry/BezierPatch.hpp
    geometry/BoundingBox.hpp
    geometry/CoordinateTransformer.hpp
    geometry/Cone.hpp
    geometry/CurvedPolygon.hpp
    geometry/Hexahedron.hpp
    geometry/KnotVector.hpp
    geometry/Line.hpp
    geometry/OrientedBoundingBox.hpp
    geometry/OrientationResult.hpp
    geometry/NURBSCurve.hpp
    geometry/NURBSPatch.hpp
    geometry/Plane.hpp
    geometry/Point.hpp
    geometry/Polygon.hpp
    geometry/Ray.hpp
    geometry/Segment.hpp
    geometry/Sphere.hpp
    geometry/Polyhedron.hpp
    geometry/Tetrahedron.hpp
    geometry/Octahedron.hpp
    geometry/Triangle.hpp
    geometry/Quadrilateral.hpp
    geometry/Vector.hpp
    geometry/construct.hpp

    ## operators
    operators/clip.hpp
    operators/closest_point.hpp
    operators/evaluate_integral.hpp
    operators/intersect.hpp
    operators/intersection_volume.hpp
    operators/orientation.hpp
    operators/squared_distance.hpp
    operators/compute_bounding_box.hpp
    operators/compute_moments.hpp
    operators/in_curved_polygon.hpp
    operators/in_polyhedron.hpp
    operators/in_polygon.hpp
    operators/in_sphere.hpp
    operators/is_convex.hpp
    operators/split.hpp
    operators/winding_number.hpp

    operators/detail/clip_impl.hpp
    operators/detail/compute_moments_impl.hpp
    operators/detail/evaluate_integral_impl.hpp
    operators/detail/fuzzy_comparators.hpp
    operators/detail/intersect_bezier_impl.hpp
    operators/detail/intersect_bounding_box_impl.hpp
    operators/detail/intersect_patch_impl.hpp
    operators/detail/intersect_impl.hpp
    operators/detail/intersect_ray_impl.hpp
    operators/detail/winding_number_2d_impl.hpp
    operators/detail/winding_number_2d_memoization.hpp
    operators/detail/winding_number_3d_impl.hpp
    operators/detail/winding_number_3d_memoization.hpp

    ## utils
    utils/ZipIndexable.hpp
    utils/ZipBoundingBox.hpp
    utils/ZipPoint.hpp
    utils/ZipRay.hpp
    utils/ZipVector.hpp
   )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set( primal_dependencies
     slic
     )

blt_list_append( TO primal_dependencies ELEMENTS mfem IF MFEM_FOUND )

axom_add_library(
    NAME                  primal
    HEADERS               ${primal_headers}
    FOLDER                axom/primal
    DEPENDS_ON            ${primal_dependencies}
    )

axom_write_unified_header(NAME    primal
                          HEADERS ${primal_headers})

axom_install_component(NAME    primal
                       HEADERS ${primal_headers}
                       )

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()
