// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/mir/reference/MIRMesh.hpp"

#include "axom/core.hpp"
#include "axom/primal.hpp"
#include "axom/fmt.hpp"

#include <sstream>

namespace axom
{
namespace mir
{
MIRMesh::MIRMesh() : m_verts(0), m_elems(0), m_numMaterials(0) { }

//--------------------------------------------------------------------------------

MIRMesh::MIRMesh(const MIRMesh& other)
  : m_verts(other.m_verts)
  , m_elems(other.m_elems)
  , m_materialVolumeFractionsElement(other.m_materialVolumeFractionsElement)
  , m_materialVolumeFractionsVertex(other.m_materialVolumeFractionsVertex)
  , m_numMaterials(other.m_numMaterials)
  , m_meshTopology(other.m_meshTopology)
{
  constructMeshRelations();
  constructVertexPositionMap(other.m_vertexPositions.data());
  constructElementParentMap(other.m_elementParentIDs.data());
  constructElementDominantMaterialMap(other.m_elementDominantMaterials.data());

  m_shapeTypes = IntMap(&m_elems);
  for(auto el : m_elems)
  {
    m_shapeTypes[el] = other.m_shapeTypes[el];
  }
}

MIRMesh& MIRMesh::operator=(const MIRMesh& other)
{
  if(this != &other)
  {
    m_verts = other.m_verts;
    m_elems = other.m_elems;

    m_meshTopology = other.m_meshTopology;

    constructMeshRelations();
    constructVertexPositionMap(other.m_vertexPositions.data());
    constructElementParentMap(other.m_elementParentIDs.data());
    constructElementDominantMaterialMap(other.m_elementDominantMaterials.data());

    m_shapeTypes = IntMap(&m_elems);
    for(auto el : m_elems)
    {
      m_shapeTypes[el] = other.m_shapeTypes[el];
    }

    m_materialVolumeFractionsElement = other.m_materialVolumeFractionsElement;
    m_materialVolumeFractionsVertex = other.m_materialVolumeFractionsVertex;

    m_numMaterials = other.m_numMaterials;
  }
  return *this;
}

//--------------------------------------------------------------------------------

void MIRMesh::initializeMesh(const VertSet _verts,
                             const ElemSet _elems,
                             const int _numMaterials,
                             const CellTopologyData& _topology,
                             const CellMapData& _mapData,
                             const std::vector<std::vector<axom::float64>>& _elementVF)
{
  // Initialize the vertex and element sets
  m_verts = _verts;
  m_elems = _elems;

  m_numMaterials = _numMaterials;

  // Intialize the mesh topology
  m_meshTopology.m_evInds = _topology.m_evInds;
  m_meshTopology.m_evBegins = _topology.m_evBegins;
  m_meshTopology.m_veInds = _topology.m_veInds;
  m_meshTopology.m_veBegins = _topology.m_veBegins;

  // Initialize the mesh relations
  constructMeshRelations();

  // Initialize the mesh's data maps
  constructVertexPositionMap(_mapData.m_vertexPositions);
  constructElementParentMap(_mapData.m_elementParents);
  constructElementDominantMaterialMap(_mapData.m_elementDominantMaterials);
  constructElementShapeTypesMap(_mapData.m_shapeTypes);

  // Initialize the element and vertex volume fraction maps
  if(_elementVF.size() > 0) constructMeshVolumeFractionsMaps(_elementVF);

  // Check validity of the sets
  SLIC_ASSERT_MSG(m_verts.isValid(), "Vertex set is not valid.");
  SLIC_ASSERT_MSG(m_elems.isValid(), "Element set is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::constructMeshRelations()
{
  // construct boundary relation from elements to vertices
  {
    using RelationBuilder = ElemToVertRelation::RelationBuilder;
    m_bdry = RelationBuilder()
               .fromSet(&m_elems)
               .toSet(&m_verts)
               .begins(RelationBuilder::BeginsSetBuilder()
                         .size(m_elems.size())
                         .data(m_meshTopology.m_evBegins.data()))
               .indices(RelationBuilder::IndicesSetBuilder()
                          .size(m_meshTopology.m_evInds.size())
                          .data(m_meshTopology.m_evInds.data()));
  }

  {
    // construct coboundary relation from vertices to elements
    using RelationBuilder = VertToElemRelation::RelationBuilder;
    m_cobdry = RelationBuilder()
                 .fromSet(&m_verts)
                 .toSet(&m_elems)
                 .begins(RelationBuilder::BeginsSetBuilder()
                           .size(m_verts.size())
                           .data(m_meshTopology.m_veBegins.data()))
                 .indices(RelationBuilder::IndicesSetBuilder()
                            .size(m_meshTopology.m_veInds.size())
                            .data(m_meshTopology.m_veInds.data()));
  }

  SLIC_ASSERT_MSG(m_bdry.isValid(), "Boundary relation is not valid.");
  SLIC_ASSERT_MSG(m_cobdry.isValid(), "Coboundary relation is not valid.");

  // SLIC_DEBUG("Elem-Vert relation has size " << m_bdry.totalSize());
  // SLIC_DEBUG("Vert-Elem relation has size " << m_cobdry.totalSize());
}

//--------------------------------------------------------------------------------

void MIRMesh::constructMeshVolumeFractionsMaps(const std::vector<std::vector<axom::float64>>& elementVF)
{
  // Clear the old maps
  m_materialVolumeFractionsElement.clear();
  m_materialVolumeFractionsVertex.clear();

  // Initialize the maps for all of the materials with the input volume fraction data for each material
  for(int matID = 0; matID < m_numMaterials; ++matID)
  {
    // Initialize the map for the current material
    m_materialVolumeFractionsElement.push_back(ScalarMap(&m_elems));

    // Copy the data for the current material
    for(int eID = 0; eID < m_elems.size(); ++eID)
    {
      m_materialVolumeFractionsElement[matID][eID] = elementVF[matID][eID];
    }

    SLIC_ASSERT_MSG(m_materialVolumeFractionsElement[matID].isValid(),
                    "Element volume fraction map is not valid.");
  }

  // Initialize the maps for all of the vertex volume fractions
  for(int matID = 0; matID < m_numMaterials; ++matID)
  {
    // Initialize the new map for the volume fractions
    m_materialVolumeFractionsVertex.push_back(ScalarMap(&m_verts));

    // Calculate the average volume fraction value for the current vertex for the current material
    for(int vID = 0; vID < m_verts.size(); ++vID)
    {
      // Compute the per vertex volume fractions for the green material
      axom::float64 sum = 0;
      auto vertexElements = m_cobdry[vID];

      for(int i = 0; i < vertexElements.size(); ++i)
      {
        auto eID = vertexElements[i];
        sum += m_materialVolumeFractionsElement[matID][eID];
      }

      m_materialVolumeFractionsVertex[matID][vID] = sum / vertexElements.size();
    }

    SLIC_ASSERT_MSG(m_materialVolumeFractionsVertex[matID].isValid(),
                    "Vertex volume fraction map is not valid.");
  }
}

//--------------------------------------------------------------------------------

void MIRMesh::constructMeshVolumeFractionsVertex(const std::vector<std::vector<axom::float64>>& vertexVF)
{
  m_materialVolumeFractionsVertex.clear();

  // Initialize the maps for all of the materials with the input volume fraction data for each vertex
  for(int matID = 0; matID < m_numMaterials; ++matID)
  {
    // Initialize the map for the current material
    m_materialVolumeFractionsVertex.push_back(ScalarMap(&m_verts));

    // Copy the data for the current material
    for(int vID = 0; vID < m_verts.size(); ++vID)
    {
      m_materialVolumeFractionsVertex[matID][vID] = vertexVF[matID][vID];
    }

    SLIC_ASSERT_MSG(m_materialVolumeFractionsVertex[matID].isValid(),
                    "Vertex volume fraction map is not valid.");
  }
}

//--------------------------------------------------------------------------------

void MIRMesh::constructVertexPositionMap(const std::vector<Point2>& data)
{
  // construct the position map on the vertices
  m_vertexPositions = PointMap(&m_verts);

  for(int vID = 0; vID < m_verts.size(); ++vID) m_vertexPositions[vID] = data[vID];

  SLIC_ASSERT_MSG(m_vertexPositions.isValid(), "Position map is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::constructElementParentMap(const std::vector<int>& elementParents)
{
  // Initialize the map for the elements' parent IDs
  m_elementParentIDs = IntMap(&m_elems);

  // Copy the data for the elements
  for(int eID = 0; eID < m_elems.size(); ++eID) m_elementParentIDs[eID] = elementParents[eID];

  SLIC_ASSERT_MSG(m_elementParentIDs.isValid(), "Element parent map is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::constructElementDominantMaterialMap(const std::vector<int>& dominantMaterials)
{
  // Initialize the map for the elements' dominant colors
  m_elementDominantMaterials = IntMap(&m_elems);

  // Copy the dat for the elements
  for(int eID = 0; eID < m_elems.size(); ++eID)
    m_elementDominantMaterials[eID] = dominantMaterials[eID];

  SLIC_ASSERT_MSG(m_elementDominantMaterials.isValid(),
                  "Element dominant materials map is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::constructElementShapeTypesMap(const std::vector<mir::Shape>& shapeTypes)
{
  // Initialize the map for the elements' dominant colors
  m_shapeTypes = IntMap(&m_elems);

  // Copy the data for the elements
  for(int eID = 0; eID < m_elems.size(); ++eID) m_shapeTypes[eID] = shapeTypes[eID];

  SLIC_ASSERT_MSG(m_shapeTypes.isValid(), "Element dominant materials map is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::print()
{
  SLIC_INFO("\n------------------------Printing Mesh Information:------------------------");
  SLIC_INFO(axom::fmt::format("number of vertices: {}", m_verts.size()));
  SLIC_INFO(axom::fmt::format("number of elements: {}", m_elems.size()));
  SLIC_INFO(axom::fmt::format("number of materials: {}", m_numMaterials));

  SLIC_INFO(axom::fmt::format("evInds: {{ {} }}", axom::fmt::join(m_meshTopology.m_evInds, " ")););
  SLIC_INFO(
    axom::fmt::format("evBegins: {{ {} }}", axom::fmt::join(m_meshTopology.m_evBegins, " ")));
  SLIC_INFO(axom::fmt::format("veInds: {{ {} }}", axom::fmt::join(m_meshTopology.m_veInds, " ")));
  SLIC_INFO(
    axom::fmt::format("veBegins: {{ {} }}", axom::fmt::join(m_meshTopology.m_veBegins, " ")));

  SLIC_INFO(axom::fmt::format("vertexPositions: {{ {} }}", axom::fmt::join(m_vertexPositions, " ")));
  SLIC_INFO(
    axom::fmt::format("elementParentIDs: {{ {} }}", axom::fmt::join(m_elementParentIDs, " ")));

  SLIC_INFO(axom::fmt::format("elementDominantMaterials: {{ {} }}",
                              axom::fmt::join(m_elementDominantMaterials, " ")));

  SLIC_INFO(axom::fmt::format("shapeTypes: {{ {} }}", axom::fmt::join(m_shapeTypes, " ")));

  // lambda to convert an array of doubles to an array of formatted strings
  auto formatted_string_array = [](const auto& arr) {
    axom::Array<std::string> result(0, arr.size());
    for(auto v : arr)
    {
      result.emplace_back(axom::fmt::format("{:.3f}", v));
    }
    return result;
  };

  {
    axom::Array<std::string> entries(0, m_materialVolumeFractionsElement.size());
    for(unsigned long i = 0; i < m_materialVolumeFractionsElement.size(); ++i)
    {
      entries.emplace_back(axom::fmt::format(
        "  {{ {} }}",
        axom::fmt::join(formatted_string_array(m_materialVolumeFractionsElement[i]), ", ")));
    }
    SLIC_INFO(axom::fmt::format("elementVolumeFractions: {{ {} }}", axom::fmt::join(entries, " ")));
  }

  {
    axom::Array<std::string> entries(0, m_materialVolumeFractionsVertex.size());
    for(unsigned long i = 0; i < m_materialVolumeFractionsVertex.size(); ++i)
    {
      entries.emplace_back(axom::fmt::format(
        "  {{ {} }}",
        axom::fmt::join(formatted_string_array(m_materialVolumeFractionsVertex[i]), ", ")));
    }
    SLIC_INFO(axom::fmt::format("vertexVolumeFractions: {{ {} }}", axom::fmt::join(entries, " ")));
  }

  SLIC_INFO("--------------------------------------------------------------------------");
}

//--------------------------------------------------------------------------------

void MIRMesh::readMeshFromFile(std::string filename)
{
  printf("Mesh reading functionality not implemented yet. Can't read file: %s", filename.c_str());

  // Read in header

  // Read in POINTS

  // Read in CELLS

  // Read in CELL_TYPES

  // Read in CELL_DATA (element volume fractions)
}

//--------------------------------------------------------------------------------

void MIRMesh::writeMeshToFile(const std::string& dirName,
                              const std::string& fileName,
                              const std::string& separator)
{
  // Ensure that the directory to write to exists
  if(!axom::utilities::filesystem::pathExists(dirName))
  {
    axom::utilities::filesystem::makeDirsForPath(dirName);
  }

  std::string outputLocation = axom::utilities::filesystem::joinPath(dirName, fileName, separator);

  std::ofstream meshfile;
  meshfile.open(outputLocation);
  std::ostream_iterator<PosType> out_it(meshfile, " ");

  // write header
  meshfile << "# vtk DataFile Version 3.0\n"
           << "vtk output\n"
           << "ASCII\n"
           << "DATASET UNSTRUCTURED_GRID\n"
           << "POINTS " << m_verts.size() << " double\n";

  // write positions
  for(int vID = 0; vID < m_verts.size(); ++vID)
  {
    auto& pt = m_vertexPositions[vID];
    for(int i = 0; i < pt.dimension(); ++i)
    {
      meshfile << pt[i] << " ";
    }
    meshfile << "\n";
  }

  meshfile << "\nCELLS " << m_elems.size() << " " << m_meshTopology.m_evInds.size() + m_elems.size();
  for(int i = 0; i < m_elems.size(); ++i)
  {
    int nVerts = m_meshTopology.m_evBegins[i + 1] - m_meshTopology.m_evBegins[i];
    meshfile << "\n" << nVerts;
    for(int j = 0; j < nVerts; ++j)
    {
      int startIndex = m_meshTopology.m_evBegins[i];
      meshfile << " " << m_meshTopology.m_evInds[startIndex + j];
    }
  }

  meshfile << "\n\nCELL_TYPES " << m_elems.size() << "\n";
  for(int i = 0; i < m_elems.size(); ++i)
  {
    if(m_shapeTypes[i] == mir::Shape::Triangle)
      meshfile << "5\n";
    else if(m_shapeTypes[i] == mir::Shape::Quad)
      meshfile << "9\n";
    else if(m_shapeTypes[i] == mir::Shape::Tetrahedron)
      meshfile << "10\n";
    else if(m_shapeTypes[i] == mir::Shape::Triangular_Prism)
      meshfile << "13\n";
    else if(m_shapeTypes[i] == mir::Shape::Pyramid)
      meshfile << "14\n";
    else if(m_shapeTypes[i] == mir::Shape::Hexahedron)
      meshfile << "12\n";
  }

  // write element materials
  meshfile << "\n\nCELL_DATA " << m_elems.size() << "\nSCALARS cellIds int 1"
           << "\nLOOKUP_TABLE default \n";
  for(int i = 0; i < m_elems.size(); ++i)
  {
    meshfile << m_elementDominantMaterials[i] << " ";
  }

  meshfile << "\n";
}

//--------------------------------------------------------------------------------

std::vector<std::vector<axom::float64>> MIRMesh::computeOriginalElementVolumeFractions()
{
  std::map<int, axom::float64> totalAreaOriginalElements;  // the total area of the original elements
  std::map<int, axom::float64> newElementAreas;  // the area of each of the generated child elements
  std::map<int, int> numChildren;  // the number of child elements generated from one of the original elements

  // Compute the total area of each element of the original mesh and also the area of each new element
  for(int eID = 0; eID < m_elems.size(); ++eID)
  {
    int numVertices = m_meshTopology.m_evBegins[eID + 1] - m_meshTopology.m_evBegins[eID];
    if(numVertices == 3)
    {
      Point2 trianglePoints[3];
      for(int i = 0; i < 3; ++i)
        trianglePoints[i] =
          m_vertexPositions[m_meshTopology.m_evInds[m_meshTopology.m_evBegins[eID] + i]];
      newElementAreas[eID] =
        computeTriangleArea(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
    }
    if(numVertices == 4)
    {
      Point2 trianglePoints[4];
      for(int i = 0; i < 4; ++i)
        trianglePoints[i] =
          m_vertexPositions[m_meshTopology.m_evInds[m_meshTopology.m_evBegins[eID] + i]];
      newElementAreas[eID] =
        computeQuadArea(trianglePoints[0], trianglePoints[1], trianglePoints[2], trianglePoints[3]);
    }

    totalAreaOriginalElements[m_elementParentIDs[eID]] += newElementAreas[eID];
    numChildren[m_elementParentIDs[eID]] += .4;
  }

  // Intialize the element volume fraction vectors
  std::vector<std::vector<axom::float64>>
    elementVolumeFractions;  // indexed as: elementVolumeFractions[material][originalElementID] = volumeFraction
  elementVolumeFractions.resize(m_numMaterials);
  for(unsigned long i = 0; i < elementVolumeFractions.size(); ++i)
    elementVolumeFractions[i].resize(totalAreaOriginalElements.size(), 0.0);

  // Compute the volume fractions for each of the original mesh elements
  for(auto itr = newElementAreas.begin(); itr != newElementAreas.end(); itr++)
  {
    int parentElementID = m_elementParentIDs[itr->first];
    int materialID = m_elementDominantMaterials[itr->first];
    elementVolumeFractions[materialID][parentElementID] +=
      (itr->second / totalAreaOriginalElements[parentElementID]);
  }

  return elementVolumeFractions;
}

//--------------------------------------------------------------------------------

axom::float64 MIRMesh::computeTriangleArea(Point2 p0, Point2 p1, Point2 p2)
{
  return primal::Triangle<double, 3>(p0, p1, p2).area();
}

//--------------------------------------------------------------------------------

axom::float64 MIRMesh::computeQuadArea(Point2 p0, Point2 p1, Point2 p2, Point2 p3)
{
  return computeTriangleArea(p0, p1, p2) + computeTriangleArea(p2, p3, p0);
}

//--------------------------------------------------------------------------------

bool MIRMesh::isValid(bool verbose) const
{
  bool bValid = true;

  // Check the topology data
  bValid = bValid && isTopologyValid(verbose);

  // Check the volume fraction data
  bValid = bValid && areVolumeFractionsValid(verbose);

  // Check the map data
  bValid = bValid && areMapsValid(verbose);

  return bValid;
}

bool MIRMesh::isTopologyValid(bool verbose) const
{
  bool bValid = true;

  // check the vertex and element sets
  bValid = bValid && m_verts.isValid(verbose);
  bValid = bValid && m_elems.isValid(verbose);

  // check the relations
  if(!m_verts.empty() && !m_elems.empty())
  {
    bValid = bValid && m_bdry.isValid(verbose);
    bValid = bValid && m_cobdry.isValid(verbose);
  }

  if(!bValid && verbose)
  {
    SLIC_WARNING("MIRMesh invalid: Invalid topology");
  }

  return bValid;
}

bool MIRMesh::areVolumeFractionsValid(bool verbose) const
{
  bool bValid = true;

  int vfElemSize = m_materialVolumeFractionsElement.size();
  if(vfElemSize != m_numMaterials)
  {
    bValid = false;
    if(verbose)
    {
      SLIC_WARNING(
        "MIRMesh invalid: Invalid number of materials in volume fraction "
        "array.");
    }
    return bValid;
  }

  // Check the sizes of the volume fraction arrays
  for(int i = 0; i < m_numMaterials; ++i)
  {
    const auto& mats = m_materialVolumeFractionsElement[i];
    if(mats.size() != m_elems.size())
    {
      bValid = false;
      if(verbose)
      {
        SLIC_WARNING("MIRMesh invalid: Material "
                     << i << " has wrong "
                     << " number of materials in volume fraction array."
                     << " Expected " << m_elems.size() << ", got " << mats.size() << ".");
      }
    }
  }

  if(!bValid)
  {
    return false;
  }

  // Check that the volume fractions sum to 1.0
  for(auto el : m_elems.positions())
  {
    double sum = 0;
    for(int m = 0; m < m_numMaterials; ++m)
    {
      sum += m_materialVolumeFractionsElement[m][el];
    }

    if(!axom::utilities::isNearlyEqual(sum, 1.))
    {
      bValid = false;
      if(verbose)
      {
        std::stringstream sstr;

        for(int m = 0; m < m_numMaterials; ++m)
        {
          sstr << m_materialVolumeFractionsElement[m][el] << "  ";
        }

        SLIC_WARNING("MIRMesh invalid: Sum of materials in element"
                     << el << " was " << sum << ". Expected 1."
                     << " Values were: " << sstr.str());
      }
    }
  }

  return bValid;
}

bool MIRMesh::areMapsValid(bool verbose) const
{
  bool bValid = true;

  // Check the positions map
  if(m_vertexPositions.size() != m_verts.size())
  {
    bValid = false;
    if(verbose)
    {
      SLIC_WARNING("MIRMesh invalid: Incorrect number of vertex positions."
                   << " Expected " << m_verts.size() << ". Got " << m_vertexPositions.size());
    }
  }

  bValid = bValid && m_vertexPositions.isValid(verbose);

  // Check the element maps
  bValid = bValid && m_elementParentIDs.isValid(verbose);
  bValid = bValid && m_elementDominantMaterials.isValid(verbose);
  bValid = bValid && m_shapeTypes.isValid(verbose);

  if(m_elementParentIDs.size() != m_elems.size())
  {
    bValid = false;
    if(verbose)
    {
      SLIC_WARNING("MIRMesh invalid: Incorrect number of elem parent IDs."
                   << " Expected " << m_elems.size() << ". Got " << m_elementParentIDs.size());
    }
  }

  if(m_elementDominantMaterials.size() != m_elems.size())
  {
    bValid = false;
    if(verbose)
    {
      SLIC_WARNING("MIRMesh invalid: Incorrect number of elem dominant mats."
                   << " Expected " << m_elems.size() << ". Got "
                   << m_elementDominantMaterials.size());
    }
  }

  if(m_shapeTypes.size() != m_elems.size())
  {
    bValid = false;
    if(verbose)
    {
      SLIC_WARNING("MIRMesh invalid: Incorrect number of elem shape types."
                   << " Expected " << m_elems.size() << ". Got " << m_shapeTypes.size());
    }
  }

  return bValid;
}

}  // namespace mir
}  // namespace axom
