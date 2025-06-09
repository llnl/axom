// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MeshTester.hpp"
#include "axom/bump.hpp"

namespace numerics = axom::numerics;
namespace utils = axom::bump::utilities;
using namespace axom::bump::views;

namespace axom
{
namespace bump
{
namespace data
{
//--------------------------------------------------------------------------------
/*!
 * \brief Calculates the percent overlap between the given circle and quad.
 * 
 * \param gridSize  The size of the uniform grid which will be sampled over to check for overlap.
 * \param circleCenter  The center point of the circle.
 * \param circleRadius  The radius of the circle.
 * \param quadP0  The upper left vertex of the quad.
 * \param quadP1  The lower left vertex of the quad.
 * \param quadP2  The lower right vertex of the quad.
 * \param quadP3  The upper right vertex of the quad.
 * 
 * /return The percent value overlap of the circle and the quad between [0, 1].
 */

template <typename PointType>
static axom::float64 calculatePercentOverlapMonteCarlo(int gridSize,
                                                       const PointType& circleCenter,
                                                       axom::float64 circleRadius,
                                                       const PointType& quadP0,
                                                       const PointType& quadP1,
                                                       const PointType& quadP2,
                                                       const PointType& quadP3)
{
  // Check if any of the quad's corners are within the circle
  auto d0Sq = primal::squared_distance(quadP0, circleCenter);
  auto d1Sq = primal::squared_distance(quadP1, circleCenter);
  auto d2Sq = primal::squared_distance(quadP2, circleCenter);
  auto d3Sq = primal::squared_distance(quadP3, circleCenter);
  auto dRSq = circleRadius * circleRadius;

  int inFlags = ((d0Sq < dRSq) ? 1 << 0 : 0) + ((d1Sq < dRSq) ? 1 << 1 : 0) +
    ((d2Sq < dRSq) ? 1 << 2 : 0) + ((d3Sq < dRSq) ? 1 << 3 : 0);
  const int allFlags = 15;
  const int noFlags = 0;

  if(inFlags == allFlags)
  {
    // The entire quad overlaps the circle
    return 1.;
  }
  else if(inFlags == noFlags)
  {
    return 0.;
  }
  else
  {
    // Some of the quad overlaps the circle, so run the Monte Carlo sampling to determine how much
    const int numSamples = std::min(gridSize, 20);
    float delta_x = axom::utilities::abs(quadP2[0] - quadP1[0]) / static_cast<float>(numSamples - 1);
    float delta_y = axom::utilities::abs(quadP0[1] - quadP1[1]) / static_cast<float>(numSamples - 1);
    int countOverlap = 0;
    for(int y = 0; y < numSamples; ++y)
    {
      for(int x = 0; x < numSamples; ++x)
      {
        PointType samplePoint({delta_x * x + quadP1[0], delta_y * y + quadP1[1]});
        if(primal::squared_distance(samplePoint, circleCenter) < dRSq) ++countOverlap;
      }
    }
    return static_cast<axom::float64>(countOverlap) /
      static_cast<axom::float64>(numSamples * numSamples);
  }
}

//--------------------------------------------------------------------------------
template <typename PointType>
int materialAtPoint(const PointType& samplePoint,
                    const PointType& sphereCenter,
                    const std::vector<axom::float64>& sphereRadii2)
{
  const int numSpheres = static_cast<int>(sphereRadii2.size());
  // Default material is always the last index
  const int defaultMaterialID = numSpheres;
  for(int cID = 0; cID < numSpheres; ++cID)
  {
    if(primal::squared_distance(samplePoint, sphereCenter) < sphereRadii2[cID])
    {
      return cID;
    }
  }
  return defaultMaterialID;
}

//--------------------------------------------------------------------------------
template <typename TopologyView, typename CoordsetView>
void generateSphericalVolumeFractions(TopologyView topologyView,
                                      CoordsetView coordsetView,
                                      int gridSize,
                                      int numSpheres,
                                      std::vector<std::vector<axom::float64>>& materialVolumeFractionsData)
{
  using PointType = typename CoordsetView::PointType;
  AXOM_ANNOTATE_SCOPE("generateSphericalVolumeFractions");

  // Generate the element volume fractions with concentric spheres
  int numMaterials = numSpheres + 1;

  // Initialize the radii^2 of the circles
  std::vector<axom::float64> sphereRadii2;
  // Note: The choice of divisors is arbitrary
  axom::float64 maxRadius = gridSize / 2.1;
  axom::float64 minRadius = gridSize / 4.0;

  axom::float64 radiusDelta;
  if(numSpheres <= 1)
    radiusDelta = (maxRadius - minRadius);
  else
    radiusDelta = (maxRadius - minRadius) / static_cast<axom::float64>(numSpheres - 1);

  for(int i = 0; i < numSpheres; ++i)
  {
    auto rad = minRadius + (i * radiusDelta);
    sphereRadii2.push_back(rad * rad);
  }

  // Initialize all material volume fractions to 0
  materialVolumeFractionsData.resize(numMaterials);
  for(int i = 0; i < numMaterials; ++i)
  {
    materialVolumeFractionsData[i].resize(topologyView.numberOfZones(), 0.);
  }

  // all spheres are centered around the same point
  const float c = static_cast<float>(gridSize / 2.0);
  const auto sphereCenter = PointType({c, c, c});
  const int numSamples = 10;
  const axom::float64 numSamples3 = numSamples * numSamples * numSamples;

  // Use the uniform sampling method to generate volume fractions for each material
  axom::Array<int> materialCount(numMaterials, 0);
  for(int eID = 0; eID < topologyView.numberOfZones(); ++eID)
  {
    const auto zone = topologyView.zone(eID);

    // Get the materials at the zone corners.
    int cornerMats[8];
    for(int i = 0; i < 8; i++)
    {
      cornerMats[i] = materialAtPoint(coordsetView[zone.getId(i)], sphereCenter, sphereRadii2);
    }
    // See whether the materials are all the same at the corners.
    bool allSame = true;
    for(int i = 1; i < 8; i++)
    {
      allSame &= cornerMats[i - 1] == cornerMats[i];
    }

    // Run the uniform sampling to determine how much of the current cell is composed of each material
    for(int i = 0; i < numMaterials; ++i)
    {
      materialCount[i] = 0;
    }

    if(allSame)
    {
      // All the samples are the same.
      materialCount[cornerMats[0]] = numSamples3;
    }
    else
    {
      // The zone looks mixed. We have to check various samples.
      const auto v0 = coordsetView[zone.getId(0)];
      const auto v1 = coordsetView[zone.getId(1)];
      const auto v3 = coordsetView[zone.getId(3)];
      const auto v4 = coordsetView[zone.getId(4)];

      float delta_x = axom::utilities::abs(v1[0] - v0[0]) / static_cast<float>(numSamples - 1);
      float delta_y = axom::utilities::abs(v3[1] - v0[1]) / static_cast<float>(numSamples - 1);
      float delta_z = axom::utilities::abs(v4[2] - v0[2]) / static_cast<float>(numSamples - 1);

      for(int z = 0; z < numSamples; ++z)
      {
        for(int y = 0; y < numSamples; ++y)
        {
          for(int x = 0; x < numSamples; ++x)
          {
            const PointType samplePoint({static_cast<float>(delta_x * x + v0[0]),
                                         static_cast<float>(delta_y * y + v0[1]),
                                         static_cast<float>(delta_z * z + v0[2])});

            const int mat = materialAtPoint(samplePoint, sphereCenter, sphereRadii2);
            materialCount[mat]++;
          }
        }
      }
    }

    // Assign the element volume fractions based on the count of the samples in each circle
    const axom::float64 scale = 1. / static_cast<axom::float64>(numSamples * numSamples * numSamples);
    for(int matID = 0; matID < numMaterials; ++matID)
    {
      materialVolumeFractionsData[matID][eID] = materialCount[matID] * scale;
    }
  }
}

//--------------------------------------------------------------------------------
void MeshTester::setStructured(bool structured) { m_structured = structured; }

//--------------------------------------------------------------------------------
void MeshTester::mesh3x3(conduit::Node& mesh)
{
  // clang-format off
  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(std::vector<float>{{
    0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.
  }});
  mesh["coordsets/coords/values/y"].set(std::vector<float>{{
    0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2., 3., 3., 3., 3.
  }});

  if(m_structured)
  {
    mesh["topologies/mesh/type"] = "structured";
    mesh["topologies/mesh/coordset"] = "coords";
    mesh["topologies/mesh/elements/dims/i"] = 3;
    mesh["topologies/mesh/elements/dims/j"] = 3;
  }
  else
  {
    mesh["topologies/mesh/type"] = "unstructured";
    mesh["topologies/mesh/coordset"] = "coords";
    mesh["topologies/mesh/elements/shape"] = "quad";
    mesh["topologies/mesh/elements/connectivity"].set(std::vector<int>{{
      0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10, 8,9,13,12, 9,10,14,13, 10,11,15,14
    }});
    mesh["topologies/mesh/elements/sizes"].set(std::vector<int>{{
      4,4,4, 4,4,4, 4,4,4
    }});
    mesh["topologies/mesh/elements/offsets"].set(std::vector<int>{{
      0, 4, 8, 12, 16, 20, 24, 28, 32
    }});
  }
  // clang-format on
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseOne(conduit::Node& mesh)
{
  mesh3x3(mesh);

  // clang-format off
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/green"] = 0;
  mesh["matsets/mat/material_map/blue"] = 1;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   0,
   0,
   0,
   0,
   0, 1,
   0, 1,
   0, 1,
   1,
   1
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   1.,
   1.,
   1.,
   0.5, 0.5,
   0.2, 0.8,
   0.2, 0.8,
   1.,
   1.
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 1, 1, 1, 2, 2, 2, 1, 1
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    //0, 1, 2, 3, 5, 7, 9, 10, 11
    0, 1, 2, 3, 4, 6, 8, 10, 11
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseTwo(conduit::Node& mesh)
{
  mesh3x3(mesh);

  // clang-format off
  constexpr int BLUE = 0;
  constexpr int RED = 1;
  constexpr int ORANGE = 2;
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/blue"] = BLUE;
  mesh["matsets/mat/material_map/red"] = RED;
  mesh["matsets/mat/material_map/orange"] = ORANGE;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   BLUE,
   BLUE,
   BLUE,
   BLUE,
   BLUE, RED, ORANGE,
   BLUE, RED,
   BLUE, ORANGE,
   RED, ORANGE,
   RED
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   1.,
   1.,
   1.,
   0.5, 0.3, 0.2,
   0.2, 0.8,
   0.2, 0.8,
   0.3, 0.7,
   1.
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 1, 1, 1, 3, 2, 2, 2, 1
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 7, 9, 11, 13
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseThree(conduit::Node& mesh)
{
  // clang-format off
  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(std::vector<float>{{
    1., 0.5, 1.5, 0., 1., 2.,
  }});
  mesh["coordsets/coords/values/y"].set(std::vector<float>{{
    2., 1., 1., 0., 0., 0.
  }});

  // NOTE: structured is not an option here since it makes a triangle mesh.
  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "tri";
  mesh["topologies/mesh/elements/connectivity"].set(std::vector<int>{{
    0,1,2, 1,3,4, 1,4,2, 2,4,5
  }});
  mesh["topologies/mesh/elements/sizes"].set(std::vector<int>{{3,3,3,3}});
  mesh["topologies/mesh/elements/offsets"].set(std::vector<int>{{0,3,6,9}});

  constexpr int BLUE = 0;
  constexpr int RED = 1;
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/blue"] = BLUE;
  mesh["matsets/mat/material_map/red"] = RED;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   RED,
   BLUE, RED,
   BLUE, RED,
   BLUE, RED
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   0.5, 0.5,
   0.8, 0.2,
   0.5, 0.5
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 2, 2, 2
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    0, 1, 3, 5
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------
template <typename TopoView, typename CoordsetView>
static void addCircleMaterial(const TopoView& topoView,
                              const CoordsetView& coordsetView,
                              conduit::Node& mesh,
                              const MeshTester::Point2& circleCenter,
                              axom::float64 circleRadius,
                              int numSamples)
{
  constexpr int GREEN = 0;
  constexpr int BLUE = 1;

  int numElements = topoView.numberOfZones();
  std::vector<float> volFracs[2];
  volFracs[GREEN].resize(numElements);
  volFracs[BLUE].resize(numElements);

  // Generate the element volume fractions for the circle
  axom::ArrayView<float> greenView(volFracs[GREEN].data(), numElements);
  axom::ArrayView<float> blueView(volFracs[BLUE].data(), numElements);
  typename CoordsetView::PointType center;
  center[0] = circleCenter[0];
  center[1] = circleCenter[1];

  const TopoView deviceTopologyView(topoView);
  axom::for_all<axom::SEQ_EXEC>(
    topoView.numberOfZones(),
    AXOM_LAMBDA(axom::IndexType zoneIndex) {
      const auto zone = deviceTopologyView.zone(zoneIndex);
      // NOTE: node ordering shuffled because function takes different node ordering.
      auto vf = calculatePercentOverlapMonteCarlo(numSamples,
                                                  center,
                                                  circleRadius,
                                                  coordsetView[zone.getId(3)],
                                                  coordsetView[zone.getId(0)],
                                                  coordsetView[zone.getId(1)],
                                                  coordsetView[zone.getId(2)]);
      greenView[zoneIndex] = vf;
      blueView[zoneIndex] = 1.0 - vf;
    });

  // Figure out the material buffers from the volume fractions.
  std::vector<int> material_ids, sizes, offsets, indices;
  std::vector<float> volume_fractions;
  for(int i = 0; i < numElements; ++i)
  {
    int nmats = 0;
    offsets.push_back(indices.size());
    if(volFracs[GREEN][i] > 0.)
    {
      material_ids.push_back(GREEN);
      volume_fractions.push_back(volFracs[GREEN][i]);
      indices.push_back(indices.size());
      nmats++;
    }
    if(volFracs[BLUE][i] > 0.)
    {
      material_ids.push_back(BLUE);
      volume_fractions.push_back(volFracs[BLUE][i]);
      indices.push_back(indices.size());
      nmats++;
    }
    sizes.push_back(nmats);
  }

  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/green"] = GREEN;
  mesh["matsets/mat/material_map/blue"] = BLUE;
  mesh["matsets/mat/material_ids"].set(material_ids);
  mesh["matsets/mat/volume_fractions"].set(volume_fractions);
  mesh["matsets/mat/sizes"].set(sizes);
  mesh["matsets/mat/offsets"].set(offsets);
  mesh["matsets/mat/indices"].set(indices);
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseFour(conduit::Node& mesh)
{
  mesh3x3(mesh);

  // Make views
  using CoordsetView = axom::bump::views::ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(utils::make_array_view<float>(mesh["coordsets/coords/values/x"]),
                            utils::make_array_view<float>(mesh["coordsets/coords/values/y"]));
  using TopoView =
    axom::bump::views::UnstructuredTopologySingleShapeView<axom::bump::views::QuadShape<int>>;
  TopoView topoView(utils::make_array_view<int>(mesh["topologies/mesh/elements/connectivity"]));

  // Add material
  const Point2 circleCenter({1.5, 1.5});
  const axom::float64 circleRadius = 1.25;
  const int numSamples = 100;
  addCircleMaterial<TopoView, CoordsetView>(topoView,
                                            coordsetView,
                                            mesh,
                                            circleCenter,
                                            circleRadius,
                                            numSamples);
}

//--------------------------------------------------------------------------------
void MeshTester::createUniformGridTestCaseMesh(int gridSize,
                                               const MeshTester::Point2& circleCenter,
                                               axom::float64 circleRadius,
                                               conduit::Node& mesh)
{
  // Generate the mesh
  generateGrid(gridSize, mesh);

  // Make views
  using CoordsetView = axom::bump::views::ExplicitCoordsetView<float, 2>;
  CoordsetView coordsetView(utils::make_array_view<float>(mesh["coordsets/coords/values/x"]),
                            utils::make_array_view<float>(mesh["coordsets/coords/values/y"]));
  using TopoView =
    axom::bump::views::UnstructuredTopologySingleShapeView<axom::bump::views::QuadShape<int>>;
  TopoView topoView(utils::make_array_view<int>(mesh["topologies/mesh/elements/connectivity"]));

  // Add material
  int numSamples = 100;
  addCircleMaterial<TopoView, CoordsetView>(topoView,
                                            coordsetView,
                                            mesh,
                                            circleCenter,
                                            circleRadius,
                                            numSamples);
}

//--------------------------------------------------------------------------------
void MeshTester::generateGrid(int gridSize, conduit::Node& mesh)
{
  AXOM_ANNOTATE_SCOPE("generateGrid");
  int nx = gridSize + 1;
  int ny = gridSize + 1;
  int nzones = gridSize * gridSize;
  int nnodes = nx * ny;

  std::vector<float> xc, yc;
  xc.reserve(nnodes);
  yc.reserve(nnodes);
  for(int j = 0; j < ny; j++)
  {
    for(int i = 0; i < nx; i++)
    {
      xc.push_back(i);
      yc.push_back(j);
    }
  }

  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(xc);
  mesh["coordsets/coords/values/y"].set(yc);

  if(m_structured)
  {
    mesh["topologies/mesh/type"] = "structured";
    mesh["topologies/mesh/coordset"] = "coords";
    mesh["topologies/mesh/elements/dims/i"] = gridSize;
    mesh["topologies/mesh/elements/dims/j"] = gridSize;
  }
  else
  {
    std::vector<int> conn, sizes, offsets;
    conn.reserve(nzones * 4);
    sizes.reserve(nzones);
    offsets.reserve(nzones);
    for(int j = 0; j < gridSize; j++)
    {
      for(int i = 0; i < gridSize; i++)
      {
        offsets.push_back(offsets.size() * 4);
        sizes.push_back(4);
        conn.push_back(j * nx + i);
        conn.push_back(j * nx + i + 1);
        conn.push_back((j + 1) * nx + i + 1);
        conn.push_back((j + 1) * nx + i);
      }
    }

    mesh["topologies/mesh/type"] = "unstructured";
    mesh["topologies/mesh/coordset"] = "coords";
    mesh["topologies/mesh/elements/shape"] = "quad";
    mesh["topologies/mesh/elements/connectivity"].set(conn);
    mesh["topologies/mesh/elements/sizes"].set(sizes);
    mesh["topologies/mesh/elements/offsets"].set(offsets);
  }
}

//--------------------------------------------------------------------------------
static void addMaterial(axom::IndexType numElements,
                        int numMaterials,
                        const std::vector<std::vector<axom::float64>>& materialVolumeFractionsData,
                        conduit::Node& mesh)
{
  AXOM_ANNOTATE_SCOPE("addMaterial");
  // Figure out the material buffers from the volume fractions.
  std::vector<int> material_ids, sizes, offsets, indices;
  std::vector<float> volume_fractions;
  std::vector<float> sums(numMaterials, 0.);
  for(axom::IndexType i = 0; i < numElements; ++i)
  {
    int nmats = 0;
    offsets.push_back(indices.size());
    for(int mat = 0; mat < numMaterials; mat++)
    {
      if(materialVolumeFractionsData[mat][i] > 0.)
      {
        material_ids.push_back(mat);
        volume_fractions.push_back(materialVolumeFractionsData[mat][i]);
        indices.push_back(indices.size());
        nmats++;

        // Keep a total of the VFs for each material.
        sums[mat] += materialVolumeFractionsData[mat][i];
      }
    }
    sizes.push_back(nmats);
  }

  // Add the material
  mesh["matsets/mat/topology"] = "mesh";
  for(int mat = 0; mat < numMaterials; mat++)
  {
    if(sums[mat] > 0.)
    {
      std::stringstream ss;
      ss << "matsets/mat/material_map/mat" << mat;
      mesh[ss.str()] = mat;
    }
  }
  mesh["matsets/mat/material_ids"].set(material_ids);
  mesh["matsets/mat/volume_fractions"].set(volume_fractions);
  mesh["matsets/mat/sizes"].set(sizes);
  mesh["matsets/mat/offsets"].set(offsets);
  mesh["matsets/mat/indices"].set(indices);
}

//--------------------------------------------------------------------------------
template <typename TopoView, typename CoordsetView>
void addConcentricCircleMaterial(const TopoView& topoView,
                                 const CoordsetView& coordsetView,
                                 const MeshTester::Point2& circleCenter,
                                 std::vector<axom::float64>& circleRadii,
                                 int numSamples,
                                 conduit::Node& mesh)
{
  AXOM_ANNOTATE_SCOPE("addConcentricCircleMaterial");
  // Generate the element volume fractions with concentric circles
  int numMaterials = circleRadii.size() + 1;
  int defaultMaterialID = numMaterials - 1;  // default material is always the last index

  // Initialize all material volume fractions to 0
  std::vector<std::vector<axom::float64>> materialVolumeFractionsData(numMaterials);
  constexpr int MAXMATERIALS = 100;
  axom::StackArray<axom::ArrayView<axom::float64>, MAXMATERIALS> matvfViews;
  for(int i = 0; i < numMaterials; ++i)
  {
    const auto len = topoView.numberOfZones();
    materialVolumeFractionsData[i].resize(len, 0.);
    matvfViews[i] = axom::ArrayView<axom::float64>(materialVolumeFractionsData[i].data(), len);
  }
  // Make a vector of radius^2.
  std::vector<axom::float64> circleRadii2;
  for(const auto& r : circleRadii) circleRadii2.push_back(r * r);

  auto circleRadii2View = axom::ArrayView<axom::float64>(circleRadii2.data(), circleRadii2.size());
  const int numCircles = circleRadii2.size();

  // Use the uniform sampling method to generate volume fractions for each material
  // Note: Assumes that the cell is a parallelogram. This could be modified via biliear interpolation
  const TopoView deviceTopologyView(topoView);
  axom::for_all<axom::SEQ_EXEC>(
    topoView.numberOfZones(),
    AXOM_LAMBDA(axom::IndexType eID) {
      const auto zone = deviceTopologyView.zone(eID);

      const auto v0 = coordsetView[zone.getId(0)];
      const auto v1 = coordsetView[zone.getId(1)];
      const auto v2 = coordsetView[zone.getId(2)];

      // Run the uniform sampling to determine how much of the current cell is composed of each material
      float delta_x = axom::utilities::abs(v1[0] - v0[0]) / (float)(numSamples - 1);
      float delta_y = axom::utilities::abs(v2[1] - v1[1]) / (float)(numSamples - 1);

      // If the corners are all in the same circle then we can skip checking for mix.
      int circle = defaultMaterialID;
      bool sameCircles = true;
      for(int c = 0; c < 4; c++)
      {
        const auto corner = coordsetView[zone.getId(c)];
        const auto dist2 = primal::squared_distance(corner, circleCenter);
        // Check which circle the point is in.
        int currentCircle = defaultMaterialID;
        for(int cID = 0; cID < numCircles; ++cID)
        {
          if(dist2 < circleRadii2View[cID])
          {
            currentCircle = cID;
            break;
          }
        }
        if(c > 0)
        {
          sameCircles &= circle == currentCircle;
        }
        circle = currentCircle;
      }

      if(sameCircles)
      {
        // All of the points were found to be in circle.
        matvfViews[circle][eID] = 1.;
      }
      else
      {
        // There was variation along the edge path so the element is mixed.
        for(int y = 0; y < numSamples; ++y)
        {
          const float yc = static_cast<float>(delta_y * y + v0[1]);
          for(int x = 0; x < numSamples; ++x)
          {
            const float xc = static_cast<float>(delta_x * x + v0[0]);
            bool isPointSampled = false;
            const auto dist2 = primal::squared_distance(MeshTester::Point2({xc, yc}), circleCenter);
            for(int cID = 0; cID < numCircles && !isPointSampled; ++cID)
            {
              if(dist2 < circleRadii2View[cID])
              {
                matvfViews[cID][eID] += 1.;
                isPointSampled = true;
              }
            }
            if(!isPointSampled)
            {
              // The point was not within any of the circles, so increment the count for the default material
              matvfViews[defaultMaterialID][eID] += 1.;
            }
          }
        }

        // Assign the element volume fractions based on the count of the samples in each circle
        const axom::float64 ns2 = static_cast<axom::float64>(numSamples * numSamples);
        for(int matID = 0; matID < numMaterials; ++matID)
        {
          matvfViews[matID][eID] /= ns2;
        }
      }
    });

  addMaterial(topoView.numberOfZones(), numMaterials, materialVolumeFractionsData, mesh);
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseFive(int gridSize, int numCircles, conduit::Node& mesh)
{
  // Generate the mesh topology
  generateGrid(gridSize, mesh);

  Point2 circleCenter({gridSize / 2.f, gridSize / 2.f});  // all circles are centered around the same point

  // Initialize the radii of the circles
  std::vector<axom::float64> circleRadii;
  axom::float64 maxRadius = gridSize / 2.4;  // Note: The choice of divisor is arbitrary
  axom::float64 minRadius = gridSize / 8;    // Note: The choice of divisor is arbitrary

  axom::float64 radiusDelta;
  if(numCircles <= 1)
    radiusDelta = (maxRadius - minRadius);
  else
    radiusDelta = (maxRadius - minRadius) / (double)(numCircles - 1);

  for(int i = 0; i < numCircles; ++i)
  {
    circleRadii.push_back(minRadius + (i * radiusDelta));
  }

  const conduit::Node& n_coordset = mesh["coordsets/coords"];
  const conduit::Node& n_topology = mesh["topologies/mesh"];

  // Make views and add materials.
  auto coordsetView = make_explicit_coordset<float, 2>::view(n_coordset);
  using CoordsetView = decltype(coordsetView);

  if(m_structured)
  {
    auto topologyView = axom::bump::views::make_structured_topology<2>::view(n_topology);
    using TopologyView = decltype(topologyView);

    addConcentricCircleMaterial<TopologyView, CoordsetView>(topologyView,
                                                            coordsetView,
                                                            circleCenter,
                                                            circleRadii,
                                                            100,
                                                            mesh);
  }
  else
  {
    using ShapeType = QuadShape<int>;
    auto topologyView = make_unstructured_single_shape_topology<ShapeType>::view(n_topology);
    using TopologyView = decltype(topologyView);

    // Add the material
    addConcentricCircleMaterial<TopologyView, CoordsetView>(topologyView,
                                                            coordsetView,
                                                            circleCenter,
                                                            circleRadii,
                                                            100,
                                                            mesh);
  }
}

//--------------------------------------------------------------------------------

int MeshTester::circleQuadCornersOverlaps(const MeshTester::Point2& circleCenter,
                                          axom::float64 circleRadius,
                                          const MeshTester::Point2& quadP0,
                                          const MeshTester::Point2& quadP1,
                                          const MeshTester::Point2& quadP2,
                                          const MeshTester::Point2& quadP3)
{
  // Check if any of the quad's corners are within the circle
  auto d0Sq = primal::squared_distance(quadP0, circleCenter);
  auto d1Sq = primal::squared_distance(quadP1, circleCenter);
  auto d2Sq = primal::squared_distance(quadP2, circleCenter);
  auto d3Sq = primal::squared_distance(quadP3, circleCenter);
  auto dRSq = circleRadius * circleRadius;

  int numCorners = 0;

  if(d0Sq < dRSq) numCorners++;
  if(d1Sq < dRSq) numCorners++;
  if(d2Sq < dRSq) numCorners++;
  if(d3Sq < dRSq) numCorners++;

  return numCorners;
}

//--------------------------------------------------------------------------------
void MeshTester::initQuadClippingTestMesh(conduit::Node& mesh)
{
  // Generate the mesh topology
  constexpr int gridSize = 3;
  generateGrid(gridSize, mesh);

  // clang-format off
  mesh["matsets/mat/topology"] = "mesh";
  mesh["matsets/mat/material_map/green"] = 0;
  mesh["matsets/mat/material_map/blue"] = 1;
  mesh["matsets/mat/material_ids"].set(std::vector<int>{{
   0,
   0,
   0,
   0, 1,
   0, 1,
   0, 1,
   1,
   1,
   1
  }});
  mesh["matsets/mat/volume_fractions"].set(std::vector<float>{{
   1.,
   1.,
   1.,
   0.5, 0.5,
   0.5, 0.5,
   0.5, 0.5,
   1.,
   1.,
   1.
  }});
  mesh["matsets/mat/sizes"].set(std::vector<int>{{
    1, 1, 1, 2, 2, 2, 1, 1, 1
  }});
  mesh["matsets/mat/offsets"].set(std::vector<int>{{
    0, 1, 2, 3, 5, 7, 9, 10, 11
  }});
  mesh["matsets/mat/indices"].set(std::vector<int>{{
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
  }});
  // clang-format on
}

//--------------------------------------------------------------------------------
void MeshTester::initTestCaseSix(int gridSize, int numSpheres, conduit::Node& mesh)
{
  // Generate the mesh topology
  generateGrid3D(gridSize, mesh);

  const conduit::Node& n_coordset = mesh["coordsets/coords"];
  const conduit::Node& n_topology = mesh["topologies/mesh"];

  auto coordsetView = make_explicit_coordset<float, 3>::view(n_coordset);

  const int numMaterials = numSpheres + 1;
  std::vector<std::vector<axom::float64>> materialVolumeFractionsData;
  if(m_structured)
  {
    auto topologyView = axom::bump::views::make_structured_topology<3>::view(n_topology);

    generateSphericalVolumeFractions(topologyView,
                                     coordsetView,
                                     gridSize,
                                     numSpheres,
                                     materialVolumeFractionsData);
    addMaterial(topologyView.numberOfZones(), numMaterials, materialVolumeFractionsData, mesh);
  }
  else
  {
    using ShapeType = HexShape<int>;
    auto topologyView = make_unstructured_single_shape_topology<ShapeType>::view(n_topology);

    generateSphericalVolumeFractions(topologyView,
                                     coordsetView,
                                     gridSize,
                                     numSpheres,
                                     materialVolumeFractionsData);
    addMaterial(topologyView.numberOfZones(), numMaterials, materialVolumeFractionsData, mesh);
  }
}

//--------------------------------------------------------------------------------
void MeshTester::generateGrid3D(int gridSize, conduit::Node& mesh)
{
  AXOM_ANNOTATE_SCOPE("generateGrid3D");
  const int nzones = gridSize * gridSize * gridSize;
  const int dims[3] = {gridSize + 1, gridSize + 1, gridSize + 1};
  const int nnodes = dims[0] * dims[1] * dims[2];

  conduit::Node& coordset = mesh["coordsets/coords"];
  conduit::Node& topo = mesh["topologies/mesh"];

  // Make coordset
  coordset["type"] = "explicit";
  conduit::Node& n_x = coordset["values/x"];
  conduit::Node& n_y = coordset["values/y"];
  conduit::Node& n_z = coordset["values/z"];
  n_x.set(conduit::DataType::float32(nnodes));
  n_y.set(conduit::DataType::float32(nnodes));
  n_z.set(conduit::DataType::float32(nnodes));
  auto* x = static_cast<float*>(n_x.data_ptr());
  auto* y = static_cast<float*>(n_y.data_ptr());
  auto* z = static_cast<float*>(n_z.data_ptr());

  for(int k = 0; k < dims[2]; k++)
  {
    for(int j = 0; j < dims[1]; j++)
    {
      for(int i = 0; i < dims[0]; i++)
      {
        *x++ = i;
        *y++ = j;
        *z++ = k;
      }
    }
  }

  // Make topology
  if(m_structured)
  {
    topo["type"] = "structured";
    topo["coordset"] = "coords";
    topo["elements/dims/i"] = gridSize;
    topo["elements/dims/j"] = gridSize;
    topo["elements/dims/k"] = gridSize;
  }
  else
  {
    topo["type"] = "unstructured";
    topo["coordset"] = "coords";
    topo["elements/shape"] = "hex";
    conduit::Node& conn = topo["elements/connectivity"];
    conduit::Node& sizes = topo["elements/sizes"];
    conduit::Node& offsets = topo["elements/offsets"];
    conn.set(conduit::DataType::int32(8 * nzones));
    sizes.set(conduit::DataType::int32(nzones));
    offsets.set(conduit::DataType::int32(nzones));
    auto* conn_ptr = static_cast<int*>(conn.data_ptr());
    auto* sizes_ptr = static_cast<int*>(sizes.data_ptr());
    auto* offsets_ptr = static_cast<int*>(offsets.data_ptr());

    for(int k = 0; k < gridSize; k++)
    {
      const int knxny = k * dims[0] * dims[1];
      const int k1nxny = (k + 1) * dims[0] * dims[1];
      for(int j = 0; j < gridSize; j++)
      {
        const int jnx = j * dims[0];
        const int j1nx = (j + 1) * dims[0];
        for(int i = 0; i < gridSize; i++)
        {
          conn_ptr[0] = knxny + jnx + i;
          conn_ptr[1] = knxny + jnx + i + 1;
          conn_ptr[2] = knxny + j1nx + i + 1;
          conn_ptr[3] = knxny + j1nx + i;
          conn_ptr[4] = k1nxny + jnx + i;
          conn_ptr[5] = k1nxny + jnx + i + 1;
          conn_ptr[6] = k1nxny + j1nx + i + 1;
          conn_ptr[7] = k1nxny + j1nx + i;

          conn_ptr += 8;
        }
      }
    }
    for(int i = 0; i < nzones; i++)
    {
      sizes_ptr[i] = 8;
      offsets_ptr[i] = 8 * i;
    }
  }
}

//--------------------------------------------------------------------------------

}  // namespace data
}  // namespace bump
}  // namespace axom
