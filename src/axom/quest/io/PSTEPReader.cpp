#include "axom/quest/io/PSTEPReader.hpp"

namespace axom
{
namespace quest
{
namespace
{
constexpr int READER_SUCCESS = 0;
constexpr int READER_FAILED = -1;

}  // end anonymous namespace

//------------------------------------------------------------------------------
PSTEPReader::PSTEPReader(MPI_Comm comm) : m_comm(comm)
{
  MPI_Comm_rank(m_comm, &m_my_rank);
  MPI_Comm_size(m_comm, &m_num_ranks);
}

//------------------------------------------------------------------------------
int PSTEPReader::read(bool validate_model)
{
  SLIC_ASSERT(m_comm != MPI_COMM_NULL);

  int rc = READER_FAILED;  // return code

  switch(m_my_rank)
  {
  // handle rank 0
  case 0:
    rc = STEPReader::read(validate_model);
    if(m_num_ranks <= 1)
    {
      return rc;
    }

    bcast_int(rc);
    if(rc == READER_SUCCESS)
    {
      // broadcast number of patches, followed by patch data
      bcast_int(m_patches.size());
      for(auto& patch : m_patches)
      {
        // broadcast u- and v- knot vector
        bcast_array(patch.getKnots_u().getArray());

        // broadcast v- knot vector
        bcast_array(patch.getKnots_v().getArray());

        // broadcast control points
        bcast_array(patch.getControlPoints());

        // broadcast rational flag and weights
        const bool isRational = bcast_bool(patch.isRational());
        if(isRational)
        {
          bcast_array(patch.getWeights());
        }

        // broadcast the trimming curves for this patch
        const int numTrimmingCurves = patch.getNumTrimmingCurves();
        bcast_int(numTrimmingCurves);
        for(auto& cur : patch.getTrimmingCurves())
        {
          // broadcast the knot vector
          bcast_array(cur.getKnots().getArray());

          // broadcast the control points
          bcast_array(cur.getControlPoints());

          // broadcast the weights (if rational)
          const bool curIsRational = cur.isRational();
          bcast_bool(curIsRational);
          if(curIsRational)
          {
            bcast_array(cur.getWeights());
          }
        }
      }
    }
    break;
  //handle other ranks
  default:
    m_stepProcessor = nullptr;
    rc = bcast_int();
    if(rc == READER_SUCCESS)
    {
      // receive and reconstruct each NURBSPatch
      const int numPatches = bcast_int();
      m_patches.clear();
      m_patches.reserve(numPatches);
      for(int i = 0; i < numPatches; ++i)
      {
        {
          // receive the u-knotvector
          axom::Array<double> uKnotsArr;
          bcast_array(uKnotsArr);

          // receive the v-knotvector
          axom::Array<double> vKnotsArr;
          bcast_array(vKnotsArr);

          // receive control points data
          NURBSPatch::CoordsMat ctrlPts;
          bcast_array(ctrlPts);

          // receive weights if rational
          const bool isRational = bcast_bool();
          if(isRational)
          {
            NURBSPatch::WeightsMat weights;
            bcast_array(weights);

            m_patches.emplace_back(NURBSPatch {ctrlPts, weights, uKnotsArr, vKnotsArr});
          }
          else
          {
            m_patches.emplace_back(NURBSPatch {ctrlPts, uKnotsArr, vKnotsArr});
          }
        }

        const int numTrimmingCurves = bcast_int();
        for(int j = 0; j < numTrimmingCurves; ++j)
        {
          // receive the knot vector
          axom::Array<double> curKnotsArr;
          bcast_array(curKnotsArr);

          // receive the control points
          NURBSCurve::CoordsVec curControlPoints;
          bcast_array(curControlPoints);

          const bool curIsRational = bcast_bool();
          if(curIsRational)
          {
            NURBSCurve::WeightsVec curWeights;
            bcast_array(curWeights);

            m_patches[i].addTrimmingCurve(NURBSCurve {curControlPoints, curWeights, curKnotsArr});
          }
          else
          {
            m_patches[i].addTrimmingCurve(NURBSCurve {curControlPoints, curKnotsArr});
          }
        }
      }
    }
    break;
  }

  return rc;
}

//------------------------------------------------------------------------------
int PSTEPReader::getTriangleMesh(axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>* mesh,
                                 double linear_deflection,
                                 double angular_deflection,
                                 bool is_relative,
                                 bool trimmed)
{
  SLIC_ASSERT(mesh != nullptr);

  int rc = READER_FAILED;

  // Rank 0 reads the STEP file and generates the triangle mesh
  if(m_my_rank == 0)
  {
    rc =
      STEPReader::getTriangleMesh(mesh, linear_deflection, angular_deflection, is_relative, trimmed);

    if(m_num_ranks <= 1)
    {
      return rc;
    }

    bcast_int(rc);
    if(rc == READER_SUCCESS)
    {
      // Broadcast vertex data
      const axom::IndexType numNodes = mesh->getNumberOfNodes();
      axom::ArrayView<double> coords_x(mesh->getCoordinateArray(mint::X_COORDINATE), numNodes);
      axom::ArrayView<double> coords_y(mesh->getCoordinateArray(mint::Y_COORDINATE), numNodes);
      axom::ArrayView<double> coords_z(mesh->getCoordinateArray(mint::Z_COORDINATE), numNodes);
      bcast_index(numNodes);
      bcast_data(coords_x);
      bcast_data(coords_y);
      bcast_data(coords_z);

      // Broadcast connectivity data
      const axom::IndexType numCells = mesh->getNumberOfCells();
      const axom::IndexType connSize = mesh->getCellNodesSize();
      axom::ArrayView<axom::IndexType> conn(mesh->getCellNodesArray(), connSize);

      bcast_index(numCells);
      bcast_index(connSize);
      bcast_data(conn);

      // broadcast patch_index field
      const bool has_patch_index_field = mesh->hasField("patch_index", mint::CELL_CENTERED);
      bcast_bool(has_patch_index_field);
      if(has_patch_index_field)
      {
        auto* field = mesh->getFieldPtr<int>("patch_index", mint::CELL_CENTERED);
        SLIC_ASSERT(field != nullptr);
        axom::ArrayView<int> patchIdxView(field, numCells);
        bcast_data(patchIdxView);
      }
    }
  }
  else
  {
    rc = bcast_int();
    if(rc == READER_SUCCESS)
    {
      // receive vertex data
      axom::IndexType numCoords = bcast_index();
      axom::Array<double> coords_x(ArrayOptions::Uninitialized {}, numCoords, numCoords);
      axom::Array<double> coords_y(ArrayOptions::Uninitialized {}, numCoords, numCoords);
      axom::Array<double> coords_z(ArrayOptions::Uninitialized {}, numCoords, numCoords);
      bcast_data(coords_x.view());
      bcast_data(coords_y.view());
      bcast_data(coords_z.view());

      // receive connectivity data
      axom::IndexType numCells = bcast_index();
      axom::IndexType connSize = bcast_index();
      axom::Array<axom::IndexType> conn(ArrayOptions::Uninitialized {}, connSize, connSize);
      bcast_data(conn.view());

      // loop through vertices and add coordinates to the input mesh
      mesh->reserveNodes(numCoords);
      for(axom::IndexType i = 0; i < numCoords; ++i)
      {
        mesh->appendNode(coords_x[i], coords_y[i], coords_z[i]);
      }

      // loop through connectivity and add each triple of indices as a cell
      mesh->reserveCells(numCells, connSize);
      for(axom::IndexType c = axom::IndexType {0}; c < numCells; ++c)
      {
        axom::IndexType verts[3] = {conn[3 * c + 0], conn[3 * c + 1], conn[3 * c + 2]};
        mesh->appendCell(verts);
      }

      // receive patch_index field
      const bool has_patch_index_field = bcast_bool();
      if(has_patch_index_field)
      {
        auto* field = mesh->createField<int>("patch_index", mint::CELL_CENTERED, numCells);

        axom::ArrayView<int> patchIdxArr(field, numCells);
        bcast_data(patchIdxArr);
      }
    }
  }

  return rc;
}

/// MPI broadcasts an integer from rank 0
int PSTEPReader::bcast_int(int value)
{
  MPI_Bcast(&value, 1, axom::mpi_traits<int>::type, 0, m_comm);
  return value;
}

/// MPI broadcasts an index from rank 0
axom::IndexType PSTEPReader::bcast_index(axom::IndexType value)
{
  MPI_Bcast(&value, 1, axom::mpi_traits<axom::IndexType>::type, 0, m_comm);
  return value;
}

/// MPI broadcasts a bool from rank 0
bool PSTEPReader::bcast_bool(bool value)
{
  int intValue = value ? 1 : 0;
  MPI_Bcast(&intValue, 1, axom::mpi_traits<int>::type, 0, m_comm);
  return static_cast<bool>(intValue);
}

}  // end namespace quest
}  // end namespace axom
