// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file quest_high_order_remap.cpp
 * \brief Demonstrates conservative field remap on 2D high order meshes
 */

// Axom includes
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"

#ifdef AXOM_USE_MFEM
  #include "mfem.hpp"
#else
  #error "This example requires mfem"
#endif

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <fstream>

namespace primal = axom::primal;
namespace spin = axom::spin;

/**
 * \brief Wrapper for a 2D mfem mesh
 *
 * Helps with conversion to BezierCurves and CurvedPolygons
 *
 * The meshwrapper assumes ownership of the wrapped mesh
 */
class MeshWrapper
{
public:
  constexpr static int DIM = 2;
  using BBox = primal::BoundingBox<double, DIM>;
  using SpacePoint = primal::Point<double, DIM>;

  using CurvedPolygonType = primal::CurvedPolygon<double, DIM>;
  using BCurve = CurvedPolygonType::BezierCurveType;

private:
  /// \brief Checks if the mesh's nodes are in the Bernstein basis
  bool isBernsteinBasis() const
  {
    auto* fec = m_mesh->GetNodalFESpace()->FEColl();

    if(fec == nullptr)
    {
      return false;
    }

    if(const mfem::H1_FECollection* h1Fec =
         dynamic_cast<const mfem::H1_FECollection*>(fec))
    {
      return h1Fec->GetBasisType() == mfem::BasisType::Positive;
    }

    if(const mfem::L2_FECollection* l2Fec =
         dynamic_cast<const mfem::L2_FECollection*>(fec))
    {
      return l2Fec->GetBasisType() == mfem::BasisType::Positive;
    }

    if(dynamic_cast<const mfem::NURBSFECollection*>(fec) ||
       dynamic_cast<const mfem::LinearFECollection*>(fec) ||
       dynamic_cast<const mfem::QuadraticPosFECollection*>(fec))
    {
      return true;
    }

    return false;
  }

public:
  MeshWrapper() : m_mesh(nullptr) { }

  MeshWrapper(mfem::Mesh* mesh) { setMesh(mesh); }

  ~MeshWrapper()
  {
    if(m_mesh != nullptr)
    {
      delete m_mesh;
      m_mesh = nullptr;
    }
  }

  /// Sets the mfem mesh pointer for this MeshWrapper instance
  void setMesh(mfem::Mesh* mesh)
  {
    SLIC_ERROR_IF(mesh == nullptr, "Mesh was null");
    m_mesh = mesh;

    // Check that mesh is high order
    SLIC_ERROR_IF(m_mesh->GetNodes() == nullptr, "The mesh must be high order.");

    // Check that mesh nodes are using Bernstein basis
    SLIC_ERROR_IF(!isBernsteinBasis(),
                  "The mesh must be in the Bernstein basis");

    const double tol = 1E-8;
    computeBoundingBoxes(1 + tol);
  }

  int numVertices() const { return m_mesh->GetNV(); }
  int numElements() const { return m_mesh->GetNE(); }

  bool hasMesh() const { return m_mesh != nullptr; }
  mfem::Mesh* getMesh() { return m_mesh; }
  const mfem::Mesh* getMesh() const { return m_mesh; }

  const BBox& elementBoundingBox(int i) const { return m_boundingBoxes[i]; }
  const BBox& meshBoundingBox() const { return m_meshBBox; }

  /*!
   * \brief Transform the mfem element into a primal::CurvedPolygon
   *
   * \param elemId The index of the element
   * \return The element as a CurvedPolygon composed of BezierCurves
   */
  CurvedPolygonType elemAsCurvedPolygon(int elemId)
  {
    SLIC_ASSERT(elemId >= 0 && elemId < numElements());

    auto* fes = m_mesh->GetNodalFESpace();
    auto* nodes = m_mesh->GetNodes();

    // Get the edge Ids for this element
    mfem::Array<int> edgeIds, edgeOrients;
    m_mesh->GetElementEdges(elemId, edgeIds, edgeOrients);
    const int nEdges = edgeIds.Size();

    CurvedPolygonType poly(nEdges);
    const int order = fes->GetOrder(0);

    // Get the grid function data associated with this edge
    mfem::Array<int> dofIndices;
    BCurve curve(order);
    for(int e = 0; e < nEdges; ++e)
    {
      // get the dof (degree of freedom) indices for this edge
      fes->GetEdgeDofs(edgeIds[e], dofIndices);

      // possibly reverse the dofs, based on the orientation
      // Note: The dofs are ordered by vertices, then by edge
      const bool bReverse = (edgeOrients[e] > 0);
      if(bReverse)
      {
        curve[0] = spacePointFromDof(dofIndices[1], fes, nodes);
        for(int p = 1; p < order; ++p)
        {
          curve[p] = spacePointFromDof(dofIndices[order - (p - 1)], fes, nodes);
        }
        curve[order] = spacePointFromDof(dofIndices[0], fes, nodes);
      }
      else
      {
        curve[0] = spacePointFromDof(dofIndices[0], fes, nodes);
        for(int p = 1; p < order; ++p)
        {
          curve[p] = spacePointFromDof(dofIndices[p + 1], fes, nodes);
        }
        curve[order] = spacePointFromDof(dofIndices[1], fes, nodes);
      }

      // Note: mfem's orientation is reversed w.r.t. primal's
      //SLIC_INFO("Elem " << elemId << " edge " << e << " -- curve: " << curve);
      //curve.reverseOrientation();
      poly[nEdges - e - 1] = curve;
    }
    //  std::cout << poly << std::endl;
    return poly;
  }

private:
  /// Get the coordinates of the point from the dof index
  SpacePoint spacePointFromDof(int idx,
                               const mfem::FiniteElementSpace* fes,
                               const mfem::GridFunction* nodes)
  {
    return SpacePoint {(*nodes)(fes->DofToVDof(idx, 0)),
                       (*nodes)(fes->DofToVDof(idx, 1))};
  }

  /*!
   * \brief Compute the element and mesh bounding boxes
   *
   * \param[in] bboxScaleFac Scale factor to increase the bounding boxes
   * \pre bboxScaleFac >= 1.
   */
  void computeBoundingBoxes(double bboxScaleFac)
  {
    SLIC_ASSERT(bboxScaleFac >= 1.);

    m_boundingBoxes.resize(numElements());
    m_meshBBox.clear();

    auto* nodes = m_mesh->GetNodes();
    auto* fes = m_mesh->GetNodalFESpace();
    mfem::Array<int> dofIndices;

    const int NE = numElements();
    for(int elem = 0; elem < NE; ++elem)
    {
      auto& bbox = m_boundingBoxes[elem];
      bbox.clear();

      // Add each dof of the element to the bbox
      // Note: positivity of Bernstein bases ensures that convex
      //       hull of element nodes contain entire element
      fes->GetElementDofs(elem, dofIndices);
      for(int i = 0; i < dofIndices.Size(); ++i)
      {
        int nIdx = dofIndices[i];
        bbox.addPoint(spacePointFromDof(nIdx, fes, nodes));
      }

      // Slightly scale the bbox to account for numerical noise
      bbox.scale(bboxScaleFac);

      m_meshBBox.addBox(bbox);
    }
  }

private:
  mfem::Mesh* m_mesh;

  std::vector<BBox> m_boundingBoxes;
  BBox m_meshBBox;
};

struct Remapper
{
public:
  constexpr static int DIM = 2;
  using CandidateList = std::vector<int>;

private:
  using GridType = spin::ImplicitGrid<DIM, axom::SEQ_EXEC, int>;

public:
  Remapper() = default;

  ~Remapper()
  {
    fecMap.DeleteData(true);
    fesMap.DeleteData(true);
  }

  mfem::GridFunction* project_to_pos_basis(const mfem::GridFunction* gf,
                                           bool& is_new)
  {
    mfem::GridFunction* out_pos_gf = nullptr;
    is_new = false;

    SLIC_ASSERT(gf != nullptr);

    const mfem::FiniteElementSpace* nodal_fe_space = gf->FESpace();
    SLIC_ERROR_IF(nodal_fe_space == nullptr,
                  "project_to_pos_basis(): nodal_fe_space is NULL!");

    const mfem::FiniteElementCollection* nodal_fe_coll = nodal_fe_space->FEColl();
    SLIC_ERROR_IF(nodal_fe_coll == nullptr,
                  "project_to_pos_basis(): nodal_fe_coll is NULL!");

    mfem::Mesh* gf_mesh = nodal_fe_space->GetMesh();
    SLIC_ERROR_IF(gf_mesh == nullptr,
                  "project_to_pos_basis(): gf_mesh is NULL!");

    int order = nodal_fe_space->GetOrder(0);
    int dim = gf_mesh->Dimension();

    auto* pos_fe_coll =
      new mfem::H1_FECollection(order, dim, mfem::BasisType::Positive);

    // {
    //   mfem::Geometry::Type geom_type = gf_mesh->GetElementBaseGeometry(0);
    //   int map_type = (nodal_fe_coll != nullptr)
    //     ? nodal_fe_coll->FiniteElementForGeometry(geom_type)->GetMapType()
    //     : static_cast<int>(mfem::FiniteElement::VALUE);

    //   fecMap.Register("src_fec", fec, true);
    //   mfem::FiniteElementCollection pos_fe_coll = *nodal_fe_coll;
    //   detail::get_pos_fec(nodal_fe_coll, order, dim, map_type);
    // }

    SLIC_ERROR_IF(pos_fe_coll == nullptr,
                  "Problem generating a positive finite element collection "
                    << "corresponding to the mesh's '" << nodal_fe_coll->Name()
                    << "' finite element collection.");

    //SLIC_INFO("Good so far... pos_fe_coll is not null. Making FESpace and GridFunction.");
    const int dims = nodal_fe_space->GetVDim();

    // Create a positive (Bernstein) grid function for the nodes
    mfem::FiniteElementSpace* pos_fe_space =
      new mfem::FiniteElementSpace(gf_mesh, pos_fe_coll, dims);
    mfem::GridFunction* pos_nodes = new mfem::GridFunction(pos_fe_space);

    // m_pos_nodes takes ownership of pos_fe_coll's memory (and pos_fe_space's memory)
    pos_nodes->MakeOwner(pos_fe_coll);

    // Project the nodal grid function onto this
    pos_nodes->ProjectGridFunction(*gf);

    out_pos_gf = pos_nodes;
    is_new = true;

    SLIC_WARNING_IF(
      out_pos_gf == nullptr,
      "project_to_pos_basis(): Construction failed;  out_pos_gf is NULL!");

    return out_pos_gf;
  }

  /*!
   * \brief Loads a mesh (source or target) and applies some simple transformations
   *
   * \param isSource Determines is we're loading the source mesh (true) or target mesh (false)
   * \param fname File pointing to an mfem mesh. If empty, we'll generate a Cartesian mesh over the unit square
   * \param offset_x Offset for translating the mesh in the x direction
   * \param offset_y Offset for translating the mesh in the y direction
   * \param scale  Factor for uniformly scaling the mesh
   * \param mref Number of uniform refinements to apply to the mesh
   * \param order Polynomial order for the mesh's nodal grid function
   */
  void loadMesh(bool isSource,
                const std::string& fname,
                double offset_x,
                double offset_y,
                double scale,
                int mref,
                int order)
  {
    mfem::Mesh* mesh = nullptr;

    // Load mesh from file or as Cartesian mesh
    if(!fname.empty())
    {
      mesh = new mfem::Mesh(fname.c_str(), 1, 1);
    }
    else
    {
      const auto quadType = mfem::Element::QUADRILATERAL;
      const bool generateEdges = true;
      const int res = 1;
      mesh = new mfem::Mesh(res, res, quadType, generateEdges);
    }

    // Apply uniform refinement
    for(int i = 0; i < mref; ++i)
    {
      mesh->UniformRefinement();
    }

    // Ensure that mesh has high order nodes
    mesh->SetCurvature(order);

    // Scale and offset mesh
    xformMesh(mesh, scale, offset_x, offset_y);

    // dump original mesh
    {
      std::ofstream file;
      file.open(axom::fmt::format("{}_mesh_orig.mfem", isSource ? "src" : "tgt"));
      mesh->Print(file);
    }

    // project nodes to Berstein basis, if necessary
    {
      bool is_mesh_gf_new {false};
      mfem::GridFunction* mesh_nodes = mesh->GetNodes();
      mfem::GridFunction* pos_mesh_nodes_ptr =
        project_to_pos_basis(mesh_nodes, is_mesh_gf_new);

      mfem::GridFunction& pos_mesh_nodes =
        (is_mesh_gf_new ? *pos_mesh_nodes_ptr : *mesh_nodes);
      mesh->NewNodes(pos_mesh_nodes, true);
    }

    // dump modified mesh
    {
      std::ofstream file;
      file.open(axom::fmt::format("{}_mesh_set.mfem", isSource ? "src" : "tgt"));
      mesh->Print(file);
    }

    // set as active source/target mesh
    if(isSource)
    {
      srcMesh.setMesh(mesh);
    }
    else
    {
      tgtMesh.setMesh(mesh);
    }

    SLIC_INFO(axom::fmt::format(
      "Loaded {} mesh from {} w/ {} elements."
      "\n\t(Slightly inflated) mesh bounding box: {}",
      isSource ? "source" : "target",
      fname.empty() ? "Cartesian grid" : fname,
      mesh->GetNE(),
      isSource ? srcMesh.meshBoundingBox() : tgtMesh.meshBoundingBox()));
  }

  /// Setup the implicit grid spatial index over the source mesh
  void initSpatialIndex()
  {
    const int NE = srcMesh.numElements();
    grid.initialize(srcMesh.meshBoundingBox(), nullptr, NE);

    for(int i = 0; i < NE; ++i)
    {
      grid.insert(srcMesh.elementBoundingBox(i), i);
    }
  }

  /*!
   * Computes the overlap areas between the elements of the target mesh
   * to the elements of the source mesh
   */
  double computeOverlapAreas()
  {
    double totalArea = 0.0;
    double correctArea = 0.0;
    const int nTargetElems = tgtMesh.numElements();

    double calcE = 0.0;
    for(int i = 0; i < nTargetElems; ++i)
    {
      // Finds the candidates from the source mesh that can intersect this target element
      auto candidates = getSourceCandidates(i);
      if(candidates.empty())
      {
        continue;
      }

      auto tgtPoly = tgtMesh.elemAsCurvedPolygon(i);
      //SLIC_INFO("Target Element: " << tgtPoly);
      correctArea += primal::area(tgtPoly);
      //SLIC_INFO("Target elem " << i
      //                         << " -- area " << primal::area(tgtPoly)
      //          //<< " -- bbox " << tgtMesh.elementBoundingBox(i)
      //          );

      double A = 0.0;
      for(int srcElem : candidates)
      {
        auto srcPoly = srcMesh.elemAsCurvedPolygon(srcElem);
        //SLIC_INFO("\tSource Element: " << srcPoly);
        //SLIC_INFO(
        //  "\tSource elem \n\t\t" << srcElem << "\n\t\t -- area " << primal::area(srcPoly)
        //              << " -- bbox " <<  srcMesh.elementBoundingBox(srcElem)
        //);

        std::vector<primal::CurvedPolygon<double, 2>> pnew;
        tgtPoly.reverseOrientation();
        srcPoly.reverseOrientation();
        if(primal::intersect(tgtPoly, srcPoly, pnew, 1e-8))
        {
          for(int i = 0; i < static_cast<int>(pnew.size()); ++i)
          {
            A -= primal::area(pnew[i]);
            //   SLIC_INFO("** Intersection area :" << -primal::area(pnew[i])
            //            );
          }
        }
        srcPoly.reverseOrientation();
        tgtPoly.reverseOrientation();
        //  SLIC_INFO("* Calculated area:  " << srcElem
        //                             << " -- area " << A
        //            //<< " -- bbox " << srcMesh.elementBoundingBox(srcElem)
        //            );
      }
      calcE += abs(primal::area(tgtPoly) - A);
      totalArea += A;
      //  SLIC_INFO("Calculated Area :" << A);
    }
    //  std::cout << inclusion << ", ";
    //  std::cout << calcE << ", ";
    //  const double tgt_scale = .999712378102150;
    //  const double tgt_trans = .0001345747181586;
    //  const double tgt_scale = .9999712378102150;
    //  const double tgt_trans = .00001345747181586;
    //  const double tgt_scale = .99999999999712378102150;
    //  const double tgt_trans = .000000000001345747181586;
    //  double trueError = (tgt_scale*tgt_scale-totalArea);
    //  std::cout << trueError << ", ";
    std::cout << "Calculated area (supermesh): " << std::fixed << totalArea
              << std::endl;
    std::cout << "Calculated area (target mesh): " << correctArea << std::endl;
    return (totalArea - correctArea);
  }

  /*!
   * Gets the IDs of the candidate elements from the source mesh
   * that might intersect with \a targetElemID from the target mesh
   */
  CandidateList getSourceCandidates(int targetElemID) const
  {
    SLIC_ASSERT(targetElemID < tgtMesh.numElements());
    using BitsetType = GridType::BitsetType;

    CandidateList filteredCandidates;

    auto& targetBBox = tgtMesh.elementBoundingBox(targetElemID);
    auto candidateBits = grid.getCandidates(targetBBox);

    // Filter the candidates; return as std vec
    for(auto idx = candidateBits.find_first(); idx != BitsetType::npos;
        idx = candidateBits.find_next(idx))
    {
      if(primal::intersect(targetBBox, srcMesh.elementBoundingBox(idx)))
        filteredCandidates.push_back(idx);
    }

    return filteredCandidates;
  }

  void outputAsSVG()
  {
    std::string header = R"html(<svg viewBox='-5 -5 10 10'
           xmlns='http://www.w3.org/2000/svg'>)html";
    std::string footer = "</svg>";

    // lambda to convert a CurvedPolygon to an SVG path string
    auto cpToSVG = [](const MeshWrapper::CurvedPolygonType& cp) {
      axom::fmt::memory_buffer out;
      bool is_first = true;

      for(auto& curve : cp.getEdges())
      {
        // Only write out first point for first edge
        if(is_first)
        {
          axom::fmt::format_to(out, "M {} {} ", curve[0][0], curve[0][1]);
          is_first = false;
        }

        switch(curve.getOrder())
        {
        case 1:
          axom::fmt::format_to(out, "L {} {} ", curve[1][0], curve[1][1]);
          break;
        case 2:
          axom::fmt::format_to(out,
                               "Q {} {}, {} {} ",
                               curve[1][0],
                               curve[1][1],
                               curve[2][0],
                               curve[2][1]);
          break;
        case 3:
          axom::fmt::format_to(out,
                               "C {} {}, {} {}, {} {} ",
                               curve[1][0],
                               curve[1][1],
                               curve[2][0],
                               curve[2][1],
                               curve[3][0],
                               curve[3][1]);
          break;
        default:
          SLIC_WARNING(
            "Unsupported case: can only output up to cubic curves as SVG.");
        }
      }
      return axom::fmt::format("    <path d='{} Z' />\n",
                               axom::fmt::to_string(out));
    };

    std::string srcGroup;
    std::string tgtGroup;
    std::string intersectionGroup;

    // output src mesh
    {
      axom::fmt::memory_buffer out;
      axom::fmt::format_to(
        out,
        "  <g id='source_mesh' stroke='black' stroke-width='.01' "
        "fill='red' fill-opacity='.7'>\n");
      auto& meshWrapper = srcMesh;
      for(int i = 0; i < meshWrapper.numElements(); ++i)
      {
        auto cp = meshWrapper.elemAsCurvedPolygon(i);
        axom::fmt::format_to(out, cpToSVG(cp));
      }
      axom::fmt::format_to(out, "  </g>\n");
      srcGroup = axom::fmt::to_string(out);
    }

    //output tgt mesh
    {
      axom::fmt::memory_buffer out;
      axom::fmt::format_to(
        out,
        "  <g id='target_mesh' stroke='black' stroke-width='.01' "
        "fill='blue' fill-opacity='.7'>\n");
      auto& meshWrapper = tgtMesh;
      for(int i = 0; i < meshWrapper.numElements(); ++i)
      {
        auto cp = meshWrapper.elemAsCurvedPolygon(i);
        axom::fmt::format_to(out, cpToSVG(cp));
      }
      axom::fmt::format_to(out, "  </g>\n");
      tgtGroup = axom::fmt::to_string(out);
    }

    //output intersection elements
    {
      double EPS = 1e-8;
      axom::fmt::memory_buffer out;
      axom::fmt::format_to(
        out,
        "  <g id='intersection_mesh' stroke='black' stroke-width='.01' "
        "fill='green' fill-opacity='.7'>\n");

      // foreach target element, find src intersections and add to group
      const int nTargetElems = tgtMesh.numElements();
      for(int i = 0; i < nTargetElems; ++i)
      {
        auto candidates = getSourceCandidates(i);
        if(candidates.empty())
        {
          continue;
        }

        auto tgtPoly = tgtMesh.elemAsCurvedPolygon(i);
        for(int srcElem : candidates)
        {
          auto srcPoly = srcMesh.elemAsCurvedPolygon(srcElem);

          srcPoly.reverseOrientation();
          tgtPoly.reverseOrientation();

          std::vector<MeshWrapper::CurvedPolygonType> pnew;
          if(primal::intersect(tgtPoly, srcPoly, pnew, EPS))
          {
            for(auto& cp : pnew)
            {
              axom::fmt::format_to(out, cpToSVG(cp));
            }
          }

          srcPoly.reverseOrientation();
          tgtPoly.reverseOrientation();
        }
      }

      axom::fmt::format_to(out, "  </g>\n");
      intersectionGroup = axom::fmt::to_string(out);
    }

    // Write the file
    {
      std::string fname = "high_order_intersections.svg";
      std::ofstream fs(fname);
      fs << header << std::endl;
      fs << srcGroup << std::endl;
      fs << tgtGroup << std::endl;
      fs << intersectionGroup << std::endl;
      fs << footer << std::endl;
    }
  }

public:
  MeshWrapper srcMesh;
  MeshWrapper tgtMesh;
  GridType grid;

private:
  // Named fields map manage memory associated with the
  // finite element collection and spaces
  mfem::NamedFieldsMap<mfem::FiniteElementCollection> fecMap;
  mfem::NamedFieldsMap<mfem::FiniteElementSpace> fesMap;

private:
  // scale and translate the vertices of the given mesh
  void xformMesh(mfem::Mesh* mesh, double sc, double off1, double off2)
  {
    // std::cout << "Transforming Mesh now" << std::endl;
    //for (int v = 0 ; v < NumVerts ; ++v)
    //  {
    //  double pt* = mesh->GetVertex(v);
    //  pt[0] = sc*pt[0] + off1;
    //  pt[1] = sc*pt[1] + off2;
    //  }
    mfem::GridFunction* mesh_nodes = mesh->GetNodes();
    //std::cout << *mesh_nodes[0] << std::endl;
    int NumDofs = mesh_nodes->Size();
    for(int e = 0; e < (NumDofs / 2); ++e)
    {
      (*mesh_nodes)[2 * e] = sc * (*mesh_nodes)[2 * e] + off1;
      (*mesh_nodes)[2 * e + 1] = sc * (*mesh_nodes)[2 * e + 1] + off2;
    }
  }
};

/// Simple struct to hold properties for initializing a mesh instance
struct MeshProps
{
  std::string file;
  std::vector<double> offset;
  double scale {1.};
  int order {2};
  int refinement {0};

  bool isDefault() const { return file.empty() && offset.empty(); }
  bool hasOffset() const { return offset.size() == 2; }

  friend std::ostream& operator<<(std::ostream& os, const MeshProps& props)
  {
    os << axom::fmt::format(
      "{{\n"
      "   file: {} \n"
      "   offset({}): {} \n"
      "   order: {} \n"
      "   scale: {} \n"
      "   refinement: {} \n"
      "}}",
      props.file.empty() ? "<>" : axom::fmt::format("'{}'", props.file),
      props.offset.size(),
      axom::fmt::join(props.offset, " "),
      props.order,
      props.scale,
      props.refinement);
    return os;
  }
};

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;

  SLIC_INFO("The application conservatively maps fields from a source \n"
            << "high order mesh to a target high order mesh!");

  MeshProps srcMesh;
  MeshProps tgtMesh;

  // Setup default meshes if user doesn't pass in data
#ifdef AXOM_DATA_DIR
  // Use a mesh from a file if we have the data directory
  namespace fs = axom::utilities::filesystem;
  MeshProps defaultSrcMesh;
  defaultSrcMesh.file = fs::joinPath(AXOM_DATA_DIR, "mfem/disc-nurbs.mesh");
  defaultSrcMesh.scale = 1.5;
  defaultSrcMesh.offset = {.01517288412347, .02571238506182};
  defaultSrcMesh.order = 3;

  MeshProps defaultTgtMesh;
  defaultTgtMesh.file = fs::joinPath(AXOM_DATA_DIR, "mfem/disc-nurbs.mesh");
  defaultTgtMesh.offset = {.001237586, -.06172376};
  defaultTgtMesh.order = 3;
#else
  // parameters for a Cartesian mesh
  MeshProps defaultSrcMesh;
  defaultSrcMesh.refinement = 1;
  defaultSrcMesh.order = 5;

  // quad mesh partially covering unit square
  MeshProps defaultTgtMesh;
  defaultSrcMesh.refinement = 1;
  defaultTgtMesh.scale = 712378102150;
  defaultTgtMesh.offset = {.1345747181586, .1345747181586};
  defaultTgtMesh.order = 5;
#endif

  // Set up and parse command line args
  axom::CLI::App app {"High order mesh intersection application"};
  {
    app.add_option("--srcFile", srcMesh.file)
      ->description("mfem mesh file for source mesh")
      ->check(axom::CLI::ExistingFile);
    app.add_option("--srcOffset", srcMesh.offset)
      ->description("offset for source mesh")
      ->expected(2);
    app.add_option("--srcScale", srcMesh.scale)
      ->description("scale for source mesh")
      ->capture_default_str();
    app.add_option("--srcOrder", srcMesh.order)
      ->description("polynomial order for source mesh")
      ->capture_default_str();
    app.add_option("--srcRef", srcMesh.refinement)
      ->description("refinement levels for source mesh")
      ->capture_default_str();

    app.add_option("--tgtFile", tgtMesh.file)
      ->description("mfem mesh file for source mesh")
      ->check(axom::CLI::ExistingFile);
    app.add_option("--tgtOffset", tgtMesh.offset)
      ->description("offset for target mesh")
      ->expected(2);
    app.add_option("--tgtScale", tgtMesh.scale)
      ->description("scale for target mesh")
      ->capture_default_str();
    app.add_option("--tgtOrder", tgtMesh.order)
      ->description("polynomial order for target mesh")
      ->capture_default_str();
    app.add_option("--tgtRef", tgtMesh.refinement)
      ->description("refinement levels for target mesh")
      ->capture_default_str();
  }
  CLI11_PARSE(app, argc, argv);

  SLIC_INFO("Running intersection tests...");
  axom::utilities::Timer timer(false);

  Remapper remap;

  // Load the source mesh
  {
    const bool isSource = true;
    if(srcMesh.isDefault())
    {
      srcMesh = defaultSrcMesh;
    }
    if(!srcMesh.hasOffset())
    {
      srcMesh.offset = {0, 0};
    }
    SLIC_INFO("Loading source mesh: " << srcMesh);
    remap.loadMesh(isSource,
                   srcMesh.file,
                   srcMesh.offset[0],
                   srcMesh.offset[1],
                   srcMesh.scale,
                   srcMesh.refinement,
                   srcMesh.order);
  }
  // Load the target mesh
  {
    const bool isSource = false;
    if(tgtMesh.isDefault())
    {
      tgtMesh = defaultTgtMesh;
    }
    if(!tgtMesh.hasOffset())
    {
      tgtMesh.offset = {0, 0};
    }
    SLIC_INFO("Loading target mesh: " << tgtMesh);
    remap.loadMesh(isSource,
                   tgtMesh.file,
                   tgtMesh.offset[0],
                   tgtMesh.offset[1],
                   tgtMesh.scale,
                   tgtMesh.refinement,
                   tgtMesh.order);
  }

  // Set up the spatial index
  remap.initSpatialIndex();

  // Computes the overlaps between elements of the target and source meshes
  timer.start();
  double area = remap.computeOverlapAreas();
  timer.stop();

  remap.outputAsSVG();

  std::cout.precision(16);
  SLIC_INFO(axom::fmt::format("Intersecting meshes: area: {}, time: {}",
                              area,
                              timer.elapsed()));

  return 0;
}
