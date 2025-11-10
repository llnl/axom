// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_BLUEPRINT_TESTING_DATA_HELPERS_HPP_
#define AXOM_BLUEPRINT_TESTING_DATA_HELPERS_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"
#include <conduit/conduit.hpp>
#include <numeric>
#include <string>
#include <vector>

namespace axom
{
namespace blueprint
{
namespace testing
{
namespace data
{
//------------------------------------------------------------------------------

/*!
 * \brief Make a new radial distance field at each of the mesh coordinates and
 *        subtract a \a dist value from it. This makes the zero value of the
 *        field the radius in which we're interested for isosurfacing.
 *
 * \param mesh The node that contains the blueprint mesh and fields.
 * \param dist The radial distance of interest.
 */
void add_distance(conduit::Node &mesh, float dist = 6.5f)
{
  // Make a new distance field.
  const conduit::Node &n_coordset = mesh["coordsets"][0];
  axom::bump::views::dispatch_coordset(n_coordset, [&](auto coordsetView) {
    using PointType = typename decltype(coordsetView)::PointType;
    using SphereType = axom::primal::Sphere<typename PointType::CoordType, PointType::DIMENSION>;
    mesh["fields/distance/topology"] = "mesh";
    mesh["fields/distance/association"] = "vertex";
    conduit::Node &n_values = mesh["fields/distance/values"];
    const auto nnodes = coordsetView.size();
    n_values.set(conduit::DataType::float32(nnodes));
    float *valuesPtr = static_cast<float *>(n_values.data_ptr());
    SphereType s(dist);
    for(int index = 0; index < nnodes; index++)
    {
      valuesPtr[index] = s.computeSignedDistance(coordsetView[index]);
    }
  });
}

/*!
 * \brief Creates a mesh of the desired \a type and dimensions \a dims and puts
 *        the values into the \a mesh node. A distance field is also added.
 *
 * \param type The type of mesh being added. See Conduit's braid function docs.
 * \param dims An array containing the dimensions.
 * \param[out] mesh The node that will contain the new mesh and fields.
 */
template <typename Dimensions>
void braid(const std::string &type, const Dimensions &dims, conduit::Node &mesh)
{
  int d[3] = {0, 0, 0};
  auto n = dims.size();
  for(decltype(n) i = 0; i < n; i++)
  {
    d[i] = dims[i];
  }
  conduit::blueprint::mesh::examples::braid(type, d[0], d[1], d[2], mesh);

  add_distance(mesh);
}

/*!
 * \brief Make a new "unibuffer" matset from the input vectors. The unibuffer
 *        matset is a style of matset in Blueprint that combines material ids
 *        and volume fractions from multiple materials shared arrays.
 *
 * \param vfA The volume fractions for material A over all zones in the mesh.
 * \param vfB The volume fractions for material B over all zones in the mesh.
 * \param vfC The volume fractions for material C over all zones in the mesh.
 * \param matnos The material numbers to use for materials A, B, C.
 * \param[out] matset The node that will contain the matset.
 */
void make_unibuffer(const std::vector<float> &vfA,
                    const std::vector<float> &vfB,
                    const std::vector<float> &vfC,
                    const std::vector<int> &matnos,
                    conduit::Node &matset)
{
  std::vector<int> material_ids;
  std::vector<float> volume_fractions;
  std::vector<int> sizes(vfA.size(), 0);
  std::vector<int> offsets(vfA.size(), 0);
  std::vector<int> indices;

  int offset = 0;
  for(size_t i = 0; i < vfA.size(); i++)
  {
    if(vfA[i] > 0.)
    {
      indices.push_back(static_cast<int>(indices.size()));
      material_ids.push_back(matnos[0]);
      volume_fractions.push_back(vfA[i]);
      sizes[i]++;
    }
    if(vfB[i] > 0.)
    {
      indices.push_back(static_cast<int>(indices.size()));
      material_ids.push_back(matnos[1]);
      volume_fractions.push_back(vfB[i]);
      sizes[i]++;
    }
    if(vfC[i] > 0.)
    {
      indices.push_back(static_cast<int>(indices.size()));
      material_ids.push_back(matnos[2]);
      volume_fractions.push_back(vfC[i]);
      sizes[i]++;
    }

    offsets[i] = offset;
    offset += sizes[i];
  }

  matset["material_ids"].set(material_ids);
  matset["volume_fractions"].set(volume_fractions);
  matset["indices"].set(indices);
  matset["sizes"].set(sizes);
  matset["offsets"].set(offsets);
}

/*!
 * \brief Make a new "multibuffer" matset from the input vectors.
 *
 * \param vfA The volume fractions for material A over all zones in the mesh.
 * \param vfB The volume fractions for material B over all zones in the mesh.
 * \param vfC The volume fractions for material C over all zones in the mesh.
 * \param matnos The material numbers to use for materials A, B, C.
 * \param[out] matset The node that will contain the matset.
 */
void make_multibuffer(const std::vector<float> &vfA,
                      const std::vector<float> &vfB,
                      const std::vector<float> &vfC,
                      const std::vector<int> &AXOM_UNUSED_PARAM(matnos),
                      conduit::Node &matset)
{
  std::vector<int> indices(vfA.size());
  std::iota(indices.begin(), indices.end(), 0);
  matset["volume_fractions/A/values"].set(vfA);
  matset["volume_fractions/A/indices"].set(indices);
  matset["volume_fractions/B/values"].set(vfB);
  matset["volume_fractions/B/indices"].set(indices);
  matset["volume_fractions/C/values"].set(vfC);
  matset["volume_fractions/C/indices"].set(indices);
}

/*!
 * \brief Make a new "element_dominant" matset from the input vectors.
 *
 * \param vfA The volume fractions for material A over all zones in the mesh.
 * \param vfB The volume fractions for material B over all zones in the mesh.
 * \param vfC The volume fractions for material C over all zones in the mesh.
 * \param matnos The material numbers to use for materials A, B, C.
 * \param[out] matset The node that will contain the matset.
 */
void make_element_dominant(const std::vector<float> &vfA,
                           const std::vector<float> &vfB,
                           const std::vector<float> &vfC,
                           const std::vector<int> &AXOM_UNUSED_PARAM(matnos),
                           conduit::Node &matset)
{
  matset["volume_fractions/A"].set(vfA);
  matset["volume_fractions/B"].set(vfB);
  matset["volume_fractions/C"].set(vfC);
}

/*!
 * \brief Make a new "material_dominant" matset from the input vectors.
 *
 * \param vfA The volume fractions for material A over all zones in the mesh.
 * \param vfB The volume fractions for material B over all zones in the mesh.
 * \param vfC The volume fractions for material C over all zones in the mesh.
 * \param matnos The material numbers to use for materials A, B, C.
 * \param[out] matset The node that will contain the matset.
 */
void make_material_dominant(const std::vector<float> &vfA,
                            const std::vector<float> &vfB,
                            const std::vector<float> &vfC,
                            const std::vector<int> &AXOM_UNUSED_PARAM(matnos),
                            conduit::Node &matset)
{
  std::vector<float> svfA, svfB, svfC;
  std::vector<int> ziA, ziB, ziC;
  const size_t n = vfA.size();
  for(size_t zi = 0; zi < n; zi++)
  {
    if(vfA[zi] > 0.f)
    {
      svfA.push_back(vfA[zi]);
      ziA.push_back(zi);
    }
    if(vfB[zi] > 0.f)
    {
      svfB.push_back(vfB[zi]);
      ziB.push_back(zi);
    }
    if(vfC[zi] > 0.f)
    {
      svfC.push_back(vfC[zi]);
      ziC.push_back(zi);
    }
  }
  matset["volume_fractions/A"].set(svfA);
  matset["volume_fractions/B"].set(svfB);
  matset["volume_fractions/C"].set(svfC);
  matset["element_ids/A"].set(ziA);
  matset["element_ids/B"].set(ziB);
  matset["element_ids/C"].set(ziC);
}

/*!
 * \brief Make a new mesh with a matset that has 3 materials.
 *
 * \param type The type of matset to create.
 * \param topoName The name of mesh topology.
 * \param dims The dimensions of the mesh.
 * \param cleanMats Whether to make a matset of only clean zones (1 mat/zone).
 * \param[out] mesh The mesh node to which a matset will be added.
 *
 *   *--------------*
 *   |            ==|
 *   |    A     ==  |
 *   |       ==     |
 *   |=======    C  |
 *   |       ==     |
 *   |    B    ==   |
 *   |           == |
 *   *--------------*
 */
template <typename Dimensions>
void make_matset(const std::string &type,
                 const std::string &topoName,
                 const Dimensions &dims,
                 bool cleanMats,
                 conduit::Node &mesh)
{
  constexpr int sampling = 10;
  int midx = sampling * dims[0] / 2;
  int midy = sampling * dims[1] / 2;

  int ksize = 0;
  int nelements = 1;
  for(int i = 0; i < dims.size(); i++)
  {
    nelements *= dims[i];
    ksize = (i == 0) ? 1 : (ksize * dims[i]);
  }
  std::vector<int> matA(nelements);
  std::vector<int> matB(nelements);
  std::vector<int> matC(nelements);

  int nk = (dims.size() == 3) ? dims[2] : 1;
  for(int k = 0; k < nk; k++)
  {
    for(int j = 0; j < sampling * dims[1]; j++)
    {
      const int jele = j / sampling;
      for(int i = 0; i < sampling * dims[0]; i++)
      {
        const int iele = i / sampling;

        bool gt1 = j >= midy;
        bool gt2 = j >= ((3. / 2.) * (i - midx) + midy);
        bool gt3 = j >= ((-2. / 5.) * (i - midx) + midy);

        int index = k * ksize + jele * dims[0] + iele;

        if(gt1 && gt2)
          matA[index]++;
        else if(!gt1 && !gt3)
          matB[index]++;
        else
          matC[index]++;
      }
    }
  }

  std::vector<float> vfA(nelements);
  std::vector<float> vfB(nelements);
  std::vector<float> vfC(nelements);
  constexpr float s2 = static_cast<float>(sampling * sampling);
  for(int k = 0; k < nk; k++)
  {
    for(int j = 0; j < sampling * dims[1]; j++)
    {
      const int jele = j / sampling;
      for(int i = 0; i < sampling * dims[0]; i++)
      {
        const int iele = i / sampling;
        int index = k * ksize + jele * dims[0] + iele;

        float vfs[3];
        vfs[0] = static_cast<float>(matA[index]) / s2;
        vfs[1] = static_cast<float>(matB[index]) / s2;
        vfs[2] = static_cast<float>(matC[index]) / s2;

        if(cleanMats)
        {
          int maxIdx = (vfs[1] > vfs[0]) ? 1 : 0;
          maxIdx = (vfs[maxIdx] > vfs[2]) ? maxIdx : 2;
          for(int ii = 0; ii < 3; ii++)
          {
            vfs[ii] = (ii == maxIdx) ? 1.f : 0.f;
          }
        }

        vfA[index] = vfs[0];
        vfB[index] = vfs[1];
        vfC[index] = vfs[2];
      }
    }
  }
  matA.clear();
  matB.clear();
  matC.clear();

  // Debugging. Add the vfs as fields.
  mesh["fields/vfA/topology"] = topoName;
  mesh["fields/vfA/association"] = "element";
  mesh["fields/vfA/values"].set(vfA);

  mesh["fields/vfB/topology"] = topoName;
  mesh["fields/vfB/association"] = "element";
  mesh["fields/vfB/values"].set(vfB);

  mesh["fields/vfC/topology"] = topoName;
  mesh["fields/vfC/association"] = "element";
  mesh["fields/vfC/values"].set(vfC);

  const std::vector<int> matnos {{22, 66, 33}};
  conduit::Node &matset = mesh["matsets/mat"];
  matset["topology"] = topoName;
  matset["material_map/A"] = matnos[0];
  matset["material_map/B"] = matnos[1];
  matset["material_map/C"] = matnos[2];

  // produce different material types.
  if(type == "unibuffer")
  {
    make_unibuffer(vfA, vfB, vfC, matnos, matset);
  }
  else if(type == "multibuffer")
  {
    make_multibuffer(vfA, vfB, vfC, matnos, matset);
  }
  else if(type == "element_dominant")
  {
    make_element_dominant(vfA, vfB, vfC, matnos, matset);
  }
  else if(type == "material_dominant")
  {
    make_material_dominant(vfA, vfB, vfC, matnos, matset);
  }
}

/*!
 * \brief Makes a new Blueprint 3D mesh made of mixed cell types.
 *
 * \param[out] mesh The node that will contain the new mesh.
 */
void mixed3d(conduit::Node &mesh)
{
  // clang-format off
  const std::vector<int> conn{{
    // tets
    0,6,1,3,
    3,6,1,9,
    6,7,1,9,
    3,9,1,4,
    9,7,1,4,
    9,7,4,10,
    // pyramids
    1,7,8,2,4,
    11,8,7,10,4,
    2,8,11,5,4,
    // wedges
    6,7,9,12,13,15,
    9,7,10,15,13,16,
    // hex
    7,13,14,8,10,16,17,11
  }};
  const std::vector<int> shapes{{
    0,0,0,0,0,0,
    1,1,1,
    2,2,
    3
  }};
  const std::vector<int> sizes{{
    4,4,4,4,4,4,
    5,5,5,
    6,6,
    8
  }};
  const std::vector<int> offsets{{
    0,4,8,12,16,20,
    24,29,34,
    39,45,
    51
  }};
  constexpr float LOW = -10.f;
  constexpr float MID = 0.f;
  constexpr float HIGH = 10.f;
  const std::vector<float> x{{
    LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH
  }};
  const std::vector<float> y{{
    LOW, LOW, LOW, MID, MID, MID, LOW, LOW, LOW, MID, MID, MID, LOW, LOW, LOW, MID, MID, MID
  }};
  const std::vector<float> z{{
    LOW, LOW, LOW, LOW, LOW, LOW, MID, MID, MID, MID, MID, MID, HIGH, HIGH, HIGH, HIGH, HIGH, HIGH
  }};
  // clang-format off

  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(x);
  mesh["coordsets/coords/values/y"].set(y);
  mesh["coordsets/coords/values/z"].set(z);
  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "mixed";
  mesh["topologies/mesh/elements/connectivity"].set(conn);
  mesh["topologies/mesh/elements/shapes"].set(shapes);
  mesh["topologies/mesh/elements/sizes"].set(sizes);
  mesh["topologies/mesh/elements/offsets"].set(offsets);
  mesh["topologies/mesh/elements/shape_map/tet"] = 0;
  mesh["topologies/mesh/elements/shape_map/pyramid"] = 1;
  mesh["topologies/mesh/elements/shape_map/wedge"] = 2;
  mesh["topologies/mesh/elements/shape_map/hex"] = 3;

  add_distance(mesh, 0.f);
}

/*!
 * \brief Make a Blueprint mesh that contains one hex.
 *
 * \param[out] hostMesh A node that contains the mesh on the host.
 */
void make_one_hex(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(
    std::vector<float> {{0., 1., 1., 0., 0., 1., 1., 0.}});
  hostMesh["coordsets/coords/values/y"].set(
    std::vector<float> {{0., 0., 1., 1., 0., 0., 1., 1.}});
  hostMesh["coordsets/coords/values/z"].set(
    std::vector<float> {{0., 0., 0., 0., 1., 1., 1., 1.}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "hex";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3, 4, 5, 6, 7}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {8});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(
    std::vector<float> {{1., -1., -1., -1., -1., -1., -1., -1.}});
}

/*!
 * \brief Make a Blueprint mesh that contains one tet.
 *
 * \param[out] hostMesh A node that contains the mesh on the host.
 */
void make_one_tet(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(std::vector<float> {{0., 0., 1., 0.}});
  hostMesh["coordsets/coords/values/y"].set(std::vector<float> {{0., 0., 0., 1.}});
  hostMesh["coordsets/coords/values/z"].set(std::vector<float> {{0., 1., 0., 0.}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "tet";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {4});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(std::vector<float> {{-1., -1., -1., 1.}});
}

/*!
 * \brief Make a Blueprint mesh that contains one pyramid.
 *
 * \param[out] hostMesh A node that contains the mesh on the host.
 */

void make_one_pyr(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(
    std::vector<float> {{0., 0., 1., 1., 0.5}});
  hostMesh["coordsets/coords/values/y"].set(
    std::vector<float> {{0., 0., 0., 0., 1.}});
  hostMesh["coordsets/coords/values/z"].set(
    std::vector<float> {{0., 1., 1., 0., 0.5}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "pyramid";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3, 4}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {5});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(
    std::vector<float> {{1., 1., -1., -1., -1.}});
}

/*!
 * \brief Make a Blueprint mesh that contains one wedge.
 *
 * \param[out] hostMesh A node that contains the mesh on the host.
 */
void make_one_wdg(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(
    std::vector<float> {{0., 0., 1., 0., 0., 1}});
  hostMesh["coordsets/coords/values/y"].set(
    std::vector<float> {{0., 0., 0., 1., 1., 1.}});
  hostMesh["coordsets/coords/values/z"].set(
    std::vector<float> {{0., 1., 0., 0., 1., 0.}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "wedge";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3, 4, 5}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {6});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(
    std::vector<float> {{1., 1., -1., -1., -1., -1.}});
}

/*!
 * \brief Make a Blueprint with a strided structured topology.
 *
 * \tparam NDIMS The number of topological mesh dimensions.
 *
 * \param[out] hostMesh A node that contains the mesh on the host.
 */
template <int NDIMS = 2>
void strided_structured(conduit::Node &hostMesh)
{
  // Total padding in each dimension 
  const conduit::index_t total_elt_pad = 4; // two on each end
  const conduit::index_t total_pt_pad = 3;  // two on the low end, one on the high end

  // Size of the window in the larger mesh
  const conduit::index_t npts_x = 4;
  const conduit::index_t npts_y = 3;
  const conduit::index_t npts_z = (NDIMS == 3) ? 3 : 0;
  const conduit::index_t spz = (NDIMS == 3) ? (npts_z + total_pt_pad) : 0;

  const conduit::index_t nelts_x = npts_x - 1;
  const conduit::index_t nelts_y = npts_y - 1;
  const conduit::index_t nelts_z = (NDIMS == 3) ? (npts_z - 1) : 0;
  const conduit::index_t sez = (NDIMS == 3) ? (nelts_z + total_elt_pad) : 0;

  // Origin: where the data starts in the arrays 
  const conduit::index_t origin_x = 2;
  const conduit::index_t origin_y = 2;
  const conduit::index_t origin_z = (NDIMS == 3) ? 2 : 0;

  conduit::Node desc;
  desc["vertex_data/shape"].set(std::vector<conduit::index_t>{{npts_x + total_pt_pad, npts_y + total_pt_pad, spz}});
  desc["vertex_data/origin"].set(std::vector<conduit::index_t>{{origin_x, origin_y, origin_z}});
  desc["element_data/shape"].set(std::vector<conduit::index_t>{{nelts_x + total_elt_pad, nelts_y + total_elt_pad, sez}});
  desc["element_data/origin"].set(std::vector<conduit::index_t>{{origin_x, origin_y, origin_z}});
  conduit::blueprint::mesh::examples::strided_structured(desc, npts_x, npts_y, npts_z, hostMesh);
}

} // end namespace data
} // end namespace testing
} // end namespace blueprint
} // end namespace axom

#endif
