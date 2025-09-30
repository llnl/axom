// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"
#include "axom/fmt.hpp"

#ifndef AXOM_USE_MFEM
  #error This file requires MFEM
#endif

#include "mfem.hpp"
#include "gtest/gtest.h"

#ifdef AXOM_USE_MPI
  #include "mpi.h"
#endif

using axom::sidre::Group;
using axom::sidre::MFEMSidreDataCollection;

constexpr double EPSILON = 1.0e-6;

std::string testName() { return ::testing::UnitTest::GetInstance()->current_test_info()->name(); }

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
  MFEMSidreDataCollection sdc(testName());
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_owning_mesh)
{
  // 1D mesh divided into 10 segments
  auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian1D(10));
  const bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc(testName(), mesh, owns_mesh_data);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_mesh)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  const bool owns_mesh_data = false;
  MFEMSidreDataCollection sdc(testName(), &mesh, owns_mesh_data);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_empty_field)
{
  MFEMSidreDataCollection sdc(testName());
  mfem::GridFunction gf;
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_partial_field)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  mfem::H1_FECollection fec(2);
  mfem::FiniteElementSpace fes(&mesh, &fec);
  mfem::GridFunction gf(&fes);

  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_update_state)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);

  // Arbitrary values for the "state" part of blueprint
  sdc.SetCycle(3);
  sdc.SetTime(1.1);
  sdc.SetTimeStep(0.4);
  // Force refresh from the superclass
  sdc.UpdateStateToDS();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_save)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc.SetComm(MPI_COMM_WORLD);
#endif

  sdc.Save();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_save_single_file)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc.SetComm(MPI_COMM_WORLD);
  sdc.SetNumFiles(1);
#endif

  sdc.Save();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_save_two_files)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc.SetMesh(MPI_COMM_WORLD, &mesh);

  // Try to go from N ranks to 2 files if we are running with more than 2 ranks
  int num_procs;
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  if(num_procs >= 2)
  {
    sdc.SetNumFiles(2);
  }
  else
  {
    sdc.SetNumFiles(1);
  }
#else
  // Just save the single file when not using MPI
  sdc.SetMesh(&mesh);
#endif

  sdc.Save();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_gf)
{
  const std::string field_name = "test_field";
  // 2D mesh divided into triangles
  auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::TRIANGLE));
  auto* fec = new mfem::H1_FECollection(1, mesh->Dimension());
  auto* fes = new mfem::FiniteElementSpace(mesh, fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), mesh, owns_mesh_data);
  auto* gf_write = new mfem::GridFunction(fes, nullptr);
  gf_write->MakeOwner(fec);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, gf_write);
  EXPECT_TRUE(sdc_writer.HasField(field_name));

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write->ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();
  EXPECT_TRUE(sdc_reader.HasField(field_name));

  // No need to reregister, it already exists
  auto* gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_gf_vdim)
{
  const std::string field_name = "test_field";
  const int vdim = 2;
  // 2D mesh divided into triangles
  auto mesh = mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec, vdim, mfem::Ordering::byVDIM);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh_data);
  // the following prevents mfem::DataCollection from deleting the mesh and gfs
  // see 'dc_reload_gf' test for a case where the DataCollection also owns the mesh and gf data
  sdc_writer.SetOwnData(false);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // No need to reregister, it already exists
  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);
  EXPECT_EQ(gf_read->FESpace()->GetOrdering(), mfem::Ordering::byVDIM);
  EXPECT_EQ(gf_read->FESpace()->GetVDim(), vdim);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

// Note: This test use Group::loadExternalData(), which is only valid (implemented)
// for sidre_hdf5 protocol. Let's skip the test in configs w/o hdf5
TEST(sidre_datacollection, dc_reload_externaldata)
{
#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::Group::loadExternalData() is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  const std::string view_name = "external_data";
  const bool owns_mesh_data = true;

  // Create external arrays for writer and reader
  axom::Array<int64_t> writer_data {1, 2, 3, 4};
  axom::Array<int64_t> reader_data {5, 6, 7, 8};

  // Create DC
  {
    auto mesh = mfem::Mesh::MakeCartesian1D(10);
    MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh_data);
    // After creation set owning to false so data doesn't get double free'd by reader and writer
    sdc_writer.SetOwnData(false);
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    sdc_writer.SetComm(MPI_COMM_WORLD);
#endif
    sdc_writer.SetCycle(0);

    axom::sidre::Group* writer_bp_group = sdc_writer.GetBPGroup();
    axom::sidre::View* writer_external_view = writer_bp_group->createView(view_name);
    writer_external_view->setExternalDataPtr(axom::sidre::INT64_ID,
                                             writer_data.size(),
                                             writer_data.data());
    EXPECT_TRUE(writer_bp_group->hasView(view_name));

    sdc_writer.Save();
  }

  // Load DC from file
  {
    MFEMSidreDataCollection sdc_reader(testName());
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
    sdc_reader.SetComm(MPI_COMM_WORLD);
#endif
    // Note: this will recreate the external view but not load the external data yet
    sdc_reader.Load();
    axom::sidre::Group* reader_bp_group = sdc_reader.GetBPGroup();
    EXPECT_TRUE(reader_bp_group->hasView(view_name));
    axom::sidre::View* reader_external_view = reader_bp_group->getView(view_name);

    // Create external buffer with wrong data and load previously saved data into it
    reader_external_view->setExternalDataPtr(reader_data.data());

    sdc_reader.LoadExternalData();
  }

  EXPECT_TRUE(writer_data.size() == reader_data.size());
  SLIC_INFO(axom::fmt::format("~~~~ {}", writer_data.size()));
  for(int i = 0; i < reader_data.size(); ++i)
  {
    SLIC_INFO(axom::fmt::format("~~~~ {} == {}", reader_data[i], writer_data[i]));
    EXPECT_TRUE(reader_data[i] == writer_data[i]);
  }
  SLIC_INFO("~~~ END");
}

TEST(sidre_datacollection, dc_reload_mesh)
{
  const std::string field_name = "test_field";
  // 2D mesh divided into triangles
  auto mesh = mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh_data);
  sdc_writer.SetOwnData(false);
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  // Save some basic info about the mesh
  const int n_verts = sdc_writer.GetMesh()->GetNV();
  const int n_ele = sdc_writer.GetMesh()->GetNE();
  const int n_bdr_ele = sdc_writer.GetMesh()->GetNBE();

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  // No mesh is used here to construct as it will be read in
  MFEMSidreDataCollection sdc_reader(testName());
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_qf)
{
  //Set up a small mesh and a couple of grid function on that mesh
  auto mesh = mfem::Mesh::MakeCartesian2D(2, 3, mfem::Element::QUADRILATERAL, 0, 2., 3.);
  mfem::LinearFECollection fec;
  mfem::FiniteElementSpace fes(&mesh, &fec);

  const int intOrder = 3;
  const int qs_vdim = 1;
  const int qv_vdim = 2;

  mfem::QuadratureSpace qspace(&mesh, intOrder);
  // We want Sidre to allocate the data for us and to own the internal data.
  // If we don't provide a nullptr, the data collection won't save the data structure
  mfem::QuadratureFunction qs(&qspace, qs_vdim);
  qs.NewDataAndSize(nullptr, qs_vdim * qspace.GetSize());

  mfem::QuadratureFunction qv(&qspace, qv_vdim);
  qv.NewDataAndSize(nullptr, qv_vdim * qspace.GetSize());

  int Nq = qs.Size();

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh_data);
  sdc_writer.SetOwnData(false);
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  // sdc_writer owns the quadrature function fields and the underlying data
  sdc_writer.RegisterQField("qs", &qs);
  sdc_writer.RegisterQField("qv", &qv);

  // The data needs to be instantiated before we save it off
  for(int i = 0; i < Nq; ++i)
  {
    qs(i) = double(i);
    qv(2 * i + 0) = double(i);
    qv(2 * i + 1) = double(Nq - i - 1);
  }

  sdc_writer.SetCycle(5);
  sdc_writer.SetTime(8.0);
  sdc_writer.Save();

#ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
#endif

  MFEMSidreDataCollection sdc_reader(testName());
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif
  sdc_reader.Load(sdc_writer.GetCycle());

  ASSERT_TRUE(sdc_reader.HasQField("qs"));
  ASSERT_TRUE(sdc_reader.HasQField("qv"));
  mfem::QuadratureFunction* reader_qs = sdc_reader.GetQField("qs");
  mfem::QuadratureFunction* reader_qv = sdc_reader.GetQField("qv");

  // order_qs should also equal order_qv in this trivial case
  EXPECT_EQ(reader_qs->GetSpace()->GetOrder(), intOrder);
  EXPECT_EQ(reader_qv->GetSpace()->GetOrder(), intOrder);

  EXPECT_EQ(reader_qs->GetVDim(), qs_vdim);
  EXPECT_EQ(reader_qv->GetVDim(), qv_vdim);

  *(reader_qs) -= qs;
  *(reader_qv) -= qv;

  EXPECT_LT(reader_qs->Norml2(), 1e-15);
  EXPECT_LT(reader_qv->Norml2(), 1e-15);
}

// Helper function to check for the existence of two sidre::Views and to check
// that they each refer to the same block of data
void checkReferentialEquality(axom::sidre::Group* grp,
                              const std::string& first,
                              const std::string& second)
{
  const bool has_first = grp->hasView(first);
  const bool has_second = grp->hasView(second);
  EXPECT_TRUE(has_first);
  EXPECT_TRUE(has_second);
  if(has_first && has_second)
  {
    // Conduct a pointer comparison to make sure that the two views actually
    // refer to the same block of memory
    const double* first_data = grp->getView(first)->getData();
    const double* second_data = grp->getView(second)->getData();
    EXPECT_EQ(first_data, second_data);
  }
}

TEST(sidre_datacollection, create_matset)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.AssociateMaterialSet("volume_fraction", "matset");

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 1.0;

  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("matsets"));
  EXPECT_TRUE(bp_grp->hasGroup("matsets/matset"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
}

TEST(sidre_datacollection, create_matset_multi_fraction)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.AssociateMaterialSet("volume_fraction", "matset");

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.7;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.3;
  sdc.RegisterField("volume_fraction_002", &vol_frac_2);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("matsets"));
  EXPECT_TRUE(bp_grp->hasGroup("matsets/matset"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/002",
                           "fields/volume_fraction_002/values");
}

TEST(sidre_datacollection, create_specset)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the specset with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 1.0;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  // Then add the species set
  const bool volume_dependent = false;
  sdc.AssociateSpeciesSet("partial_density", "specset", "matset", volume_dependent);

  mfem::GridFunction partial_density_1_1(&fes);
  partial_density_1_1 = 1.0;

  sdc.RegisterField("partial_density_001_001", &partial_density_1_1);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("specsets"));
  EXPECT_TRUE(bp_grp->hasGroup("specsets/specset"));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/volume_dependent"));
  EXPECT_FALSE(
    static_cast<std::int8_t>(bp_grp->getView("specsets/specset/volume_dependent")->getScalar()));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/matset"));
  EXPECT_EQ(std::string(bp_grp->getView("specsets/specset/matset")->getString()), "matset");
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/matset_values/001/001"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/001/001",
                           "fields/partial_density_001_001/values");
}

TEST(sidre_datacollection, create_specset_multi_fraction)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the specset with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.7;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);
  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.3;
  sdc.RegisterField("volume_fraction_002", &vol_frac_2);

  // Then add the species set
  const bool volume_dependent = false;
  sdc.AssociateSpeciesSet("partial_density", "specset", "matset", volume_dependent);

  mfem::GridFunction partial_density_1_1(&fes);
  partial_density_1_1 = 0.4;
  sdc.RegisterField("partial_density_001_001", &partial_density_1_1);
  mfem::GridFunction partial_density_1_2(&fes);
  partial_density_1_2 = 0.3;
  sdc.RegisterField("partial_density_001_002", &partial_density_1_2);
  mfem::GridFunction partial_density_2_1(&fes);
  partial_density_2_1 = 0.2;
  sdc.RegisterField("partial_density_002_001", &partial_density_2_1);
  mfem::GridFunction partial_density_2_2(&fes);
  partial_density_2_2 = 0.1;
  sdc.RegisterField("partial_density_002_002", &partial_density_2_2);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("specsets"));
  EXPECT_TRUE(bp_grp->hasGroup("specsets/specset"));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/volume_dependent"));
  EXPECT_FALSE(
    static_cast<std::int8_t>(bp_grp->getView("specsets/specset/volume_dependent")->getScalar()));
  EXPECT_TRUE(bp_grp->hasView("specsets/specset/matset"));
  EXPECT_EQ(std::string(bp_grp->getView("specsets/specset/matset")->getString()), "matset");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/002",
                           "fields/volume_fraction_002/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/001/001",
                           "fields/partial_density_001_001/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/001/002",
                           "fields/partial_density_001_002/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/002/001",
                           "fields/partial_density_002_001/values");
  checkReferentialEquality(bp_grp,
                           "specsets/specset/matset_values/002/002",
                           "fields/partial_density_002_002/values");
}

TEST(sidre_datacollection, create_material_dependent_field)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the field with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 1.0;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);

  // Then mark the field as material-dependent
  sdc.AssociateMaterialDependentField("density", "matset");

  mfem::GridFunction density_independent(&fes);
  density_independent = 1.0;
  sdc.RegisterField("density", &density_independent);

  mfem::GridFunction density_dependent(&fes);
  density_dependent = 0.2;
  sdc.RegisterField("density_001", &density_dependent);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("fields/density"));
  EXPECT_TRUE(bp_grp->hasGroup("fields/density/matset_values"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp, "fields/density/matset_values/001", "fields/density_001/values");
}

TEST(sidre_datacollection, create_material_dependent_field_multi_fraction)
{
  // 1D mesh divided into 10 segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // Add a material set so we can associate the field with it
  sdc.AssociateMaterialSet("volume_fraction", "matset");
  mfem::GridFunction vol_frac_1(&fes);
  vol_frac_1 = 0.65;
  sdc.RegisterField("volume_fraction_001", &vol_frac_1);
  mfem::GridFunction vol_frac_2(&fes);
  vol_frac_2 = 0.35;
  sdc.RegisterField("volume_fraction_002", &vol_frac_2);

  // Then mark the field as material-dependent
  sdc.AssociateMaterialDependentField("density", "matset");

  mfem::GridFunction density_indenpendent(&fes);
  density_indenpendent = 1.0;
  sdc.RegisterField("density", &density_indenpendent);

  mfem::GridFunction density_dependent_1(&fes);
  density_dependent_1 = 0.2;
  sdc.RegisterField("density_001", &density_dependent_1);

  mfem::GridFunction density_dependent_2(&fes);
  density_dependent_2 = 0.3;
  sdc.RegisterField("density_002", &density_dependent_2);

  EXPECT_TRUE(sdc.verifyMeshBlueprint());

  const auto bp_grp = sdc.GetBPGroup();
  EXPECT_TRUE(bp_grp->hasGroup("fields/density"));
  EXPECT_TRUE(bp_grp->hasGroup("fields/density/matset_values"));
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/001",
                           "fields/volume_fraction_001/values");
  checkReferentialEquality(bp_grp,
                           "matsets/matset/volume_fractions/002",
                           "fields/volume_fraction_002/values");
  checkReferentialEquality(bp_grp, "fields/density/matset_values/001", "fields/density_001/values");
  checkReferentialEquality(bp_grp, "fields/density/matset_values/002", "fields/density_002/values");
}

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)

TEST(sidre_datacollection, dc_alloc_parmesh)
{
  // case 1: serial and parallel mesh are not pointers
  {
    constexpr bool owns_mesh_data = true;

    // 1D mesh divided into 10 segments
    auto mesh = mfem::Mesh::MakeCartesian1D(10);
    mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
    MFEMSidreDataCollection sdc(testName(), &parmesh, owns_mesh_data);

    // Must set this to avoid mfem's DataCollection from deleting the parmesh
    sdc.SetOwnData(false);

    EXPECT_TRUE(sdc.verifyMeshBlueprint());
  }

  // case 2: serial and parallel meshes are pointers
  // data collection owns and deletes parmesh
  {
    constexpr bool owns_mesh_data = true;

    // 1D mesh divided into 10 segments
    auto* mesh = new mfem::Mesh(mfem::Mesh::MakeCartesian1D(10));
    auto* parmesh = new mfem::ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    MFEMSidreDataCollection sdc(testName(), parmesh, owns_mesh_data);

    EXPECT_TRUE(sdc.verifyMeshBlueprint());
  }

  // case 3: data collection doesn't own mesh/data
  {
    constexpr bool owns_mesh_data = false;

    // 1D mesh divided into 10 segments
    auto mesh = mfem::Mesh::MakeCartesian1D(10);
    mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);

    MFEMSidreDataCollection sdc(testName(), &parmesh, owns_mesh_data);
    EXPECT_TRUE(sdc.verifyMeshBlueprint());
  }
}

struct ParMeshGroupData
{
  int n_verts = 0;
  int n_edges = 0;
  int n_triangles = 0;
  int n_quads = 0;
  bool operator==(const ParMeshGroupData& other) const
  {
    return (n_verts == other.n_verts) && (n_edges == other.n_edges) &&
      (n_triangles == other.n_triangles) && (n_quads == other.n_quads);
  }
};

static std::vector<ParMeshGroupData> getGroupData(const mfem::ParMesh& parmesh)
{
  const int dim = parmesh.Dimension();
  // Zeroth group doesn't matter
  std::vector<ParMeshGroupData> result(parmesh.GetNGroups() - 1);
  // remove when marked const in MFEM
  auto& non_const_parmesh = const_cast<mfem::ParMesh&>(parmesh);

  for(std::size_t i = 1; i <= result.size(); i++)
  {
    result[i - 1].n_verts = non_const_parmesh.GroupNVertices(i);
    if(dim >= 2)
    {
      result[i - 1].n_edges = non_const_parmesh.GroupNEdges(i);
      if(dim >= 3)
      {
        result[i - 1].n_triangles = non_const_parmesh.GroupNTriangles(i);
        result[i - 1].n_quads = non_const_parmesh.GroupNQuadrilaterals(i);
      }
    }
  }

  return result;
}

/**
 * @brief Helper method for testing that a parallel mesh is reconstructed correctly
 * @param [in] base_mesh The serial mesh object to distribute, save, and then reload
 * @param [in] part_method The partitioning method to use - in [0, 5]
 */
static void testParallelMeshReload(mfem::Mesh& base_mesh, const int part_method = 1)
{
  auto* parmesh = new mfem::ParMesh(MPI_COMM_WORLD, base_mesh, nullptr, part_method);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), parmesh, owns_mesh_data);

  // Save some basic info about the mesh
  const int n_verts = sdc_writer.GetMesh()->GetNV();
  const int n_ele = sdc_writer.GetMesh()->GetNE();
  const int n_bdr_ele = sdc_writer.GetMesh()->GetNBE();

  // ParMesh-specific info
  auto writer_pmesh = dynamic_cast<mfem::ParMesh*>(sdc_writer.GetMesh());
  ASSERT_NE(writer_pmesh, nullptr);
  const int n_groups = writer_pmesh->GetNGroups();
  const int n_face_neighbors = writer_pmesh->GetNFaceNeighbors();
  const int n_shared_faces = writer_pmesh->GetNSharedFaces();

  // Group-specific info
  auto writer_group_data = getGroupData(*writer_pmesh);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(testName());

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  // Make sure the ParMesh was reconstructed - even on one rank, this
  // will still be a ParMesh
  auto reader_pmesh = dynamic_cast<mfem::ParMesh*>(sdc_reader.GetMesh());
  ASSERT_NE(reader_pmesh, nullptr);
  EXPECT_EQ(reader_pmesh->GetNGroups(), n_groups);
  EXPECT_EQ(reader_pmesh->GetNFaceNeighbors(), n_face_neighbors);
  EXPECT_EQ(reader_pmesh->GetNSharedFaces(), n_shared_faces);

  auto reader_group_data = getGroupData(*reader_pmesh);
  EXPECT_EQ(writer_group_data, reader_group_data);
  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());

  // Make sure that orientation checks pass
  // These checks assume that the test mesh is orientable and are intended to make sure
  // that Mesh::Nodes was set up correctly - otherwise these can fail for periodic (and HO?) meshes
  // QUESTION: What does this method do for a non-orientable mesh?
  EXPECT_EQ(sdc_reader.GetMesh()->CheckElementOrientation(), 0);
  EXPECT_EQ(sdc_reader.GetMesh()->CheckBdrElementOrientation(), 0);
}

/**
 * @brief Wrapper for testing mesh reconstruction across all partitioning methods
 * @param [in] base_mesh The serial mesh object to distribute, save, and then reload
 */
static void testParallelMeshReloadAllPartitionings(mfem::Mesh& base_mesh)
{
  // MFEM supports partition methods [0, 5]
  constexpr int MAX_PART_METHOD = 5;
  for(int part_method = 0; part_method <= MAX_PART_METHOD; part_method++)
  {
    testParallelMeshReload(base_mesh, part_method);
  }
}

TEST(sidre_datacollection, dc_par_reload_gf)
{
  const std::string field_name = "test_field";
  // 3D tet mesh
  auto mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), &parmesh, owns_mesh_data);
  sdc_writer.SetOwnData(false);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  mfem::ParGridFunction gf_write(&parfes, static_cast<double*>(nullptr));

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  MFEMSidreDataCollection sdc_reader(testName());

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.Load();

  auto gf_read = sdc_reader.GetField(field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(gf_read));

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_par_reload_gf_ordering)
{
  const std::string first_field_name = "test_field_1";
  const std::string second_field_name = "test_field_2";
  const std::string third_field_name = "test_field_3";

  // 3D tet mesh
  auto mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace first_parfes(&parmesh, &fec, 1, mfem::Ordering::byNODES);
  mfem::ParFiniteElementSpace second_parfes(&parmesh, &fec, 3, mfem::Ordering::byVDIM);
  mfem::ParFiniteElementSpace third_parfes(&parmesh, &fec, 3, mfem::Ordering::byNODES);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh_data = true;
  MFEMSidreDataCollection sdc_writer(testName(), &parmesh, owns_mesh_data);
  sdc_writer.SetOwnData(false);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  mfem::ParGridFunction first_gf_write(&first_parfes, static_cast<double*>(nullptr));
  mfem::ParGridFunction second_gf_write(&second_parfes, static_cast<double*>(nullptr));
  mfem::ParGridFunction third_gf_write(&third_parfes, static_cast<double*>(nullptr));

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(first_field_name, &first_gf_write);
  sdc_writer.RegisterField(second_field_name, &second_gf_write);
  sdc_writer.RegisterField(third_field_name, &third_gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  first_gf_write.ProjectCoefficient(three_and_a_half);

  mfem::ConstantCoefficient five_and_a_half(5.5);
  second_gf_write.ProjectCoefficient(five_and_a_half);

  mfem::ConstantCoefficient seven_and_a_half(5.5);
  third_gf_write.ProjectCoefficient(seven_and_a_half);

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  MFEMSidreDataCollection sdc_reader(testName());

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.Load();

  auto first_gf_read = sdc_reader.GetField(first_field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(first_gf_read));
  auto second_gf_read = sdc_reader.GetField(second_field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(second_gf_read));
  auto third_gf_read = sdc_reader.GetField(third_field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(third_gf_read));

  EXPECT_EQ(first_gf_read->FESpace()->GetOrdering(), mfem::Ordering::byNODES);
  EXPECT_EQ(second_gf_read->FESpace()->GetOrdering(), mfem::Ordering::byVDIM);
  EXPECT_EQ(third_gf_read->FESpace()->GetOrdering(), mfem::Ordering::byNODES);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(first_gf_read->ComputeL2Error(three_and_a_half), EPSILON);
  EXPECT_LT(second_gf_read->ComputeL2Error(five_and_a_half), EPSILON);
  EXPECT_LT(third_gf_read->ComputeL2Error(seven_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_par_reload_multi_datastore)
{
  const std::string first_coll_name = testName() + "first";
  const std::string second_coll_name = testName() + "second";
  const std::string field_name = "test_field";
  const std::string useless_view_name = "useless_view";

  // 3D tet mesh
  auto mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::ParMesh first_parmesh(MPI_COMM_WORLD, mesh);
  mfem::ParMesh second_parmesh(MPI_COMM_WORLD, mesh);

  mfem::H1_FECollection fec(1, first_parmesh.Dimension());

  mfem::ParFiniteElementSpace first_parfes(&first_parmesh, &fec);
  mfem::ParFiniteElementSpace second_parfes(&second_parmesh, &fec);

  axom::sidre::DataStore ds_write;

  // We want to make sure this isn't restored
  ds_write.getRoot()->createViewString(useless_view_name, "useless_data");

  auto first_global_grp = ds_write.getRoot()->createGroup(first_coll_name + "_global");
  auto first_bp_index_grp = first_global_grp->createGroup("blueprint_index/" + first_coll_name);
  auto first_domain_grp = ds_write.getRoot()->createGroup(first_coll_name);

  auto second_global_grp = ds_write.getRoot()->createGroup(second_coll_name + "_global");
  auto second_bp_index_grp = second_global_grp->createGroup("blueprint_index/" + second_coll_name);
  auto second_domain_grp = ds_write.getRoot()->createGroup(second_coll_name);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh_data = true;
  MFEMSidreDataCollection first_sdc_writer(first_coll_name,
                                           first_bp_index_grp,
                                           first_domain_grp,
                                           owns_mesh_data);
  first_sdc_writer.SetComm(MPI_COMM_WORLD);
  first_sdc_writer.SetMesh(&first_parmesh);

  MFEMSidreDataCollection second_sdc_writer(second_coll_name,
                                            second_bp_index_grp,
                                            second_domain_grp,
                                            owns_mesh_data);
  second_sdc_writer.SetComm(MPI_COMM_WORLD);
  second_sdc_writer.SetMesh(&second_parmesh);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  mfem::ParGridFunction first_gf_write(&first_parfes, static_cast<double*>(nullptr));
  mfem::ParGridFunction second_gf_write(&second_parfes, static_cast<double*>(nullptr));

  // Register to allocate storage internally, then write to it
  first_sdc_writer.RegisterField(field_name, &first_gf_write);
  second_sdc_writer.RegisterField(field_name, &second_gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  first_gf_write.ProjectCoefficient(three_and_a_half);

  mfem::ConstantCoefficient five_and_a_half(5.5);
  second_gf_write.ProjectCoefficient(five_and_a_half);

  first_sdc_writer.SetCycle(0);
  first_sdc_writer.Save();

  second_sdc_writer.SetCycle(0);
  second_sdc_writer.Save();

  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  axom::sidre::DataStore ds_read;

  first_global_grp = ds_read.getRoot()->createGroup(first_coll_name + "_global");
  first_bp_index_grp = first_global_grp->createGroup("blueprint_index/" + first_coll_name);
  first_domain_grp = ds_read.getRoot()->createGroup(first_coll_name);

  second_global_grp = ds_read.getRoot()->createGroup(second_coll_name + "_global");
  second_bp_index_grp = second_global_grp->createGroup("blueprint_index/" + second_coll_name);
  second_domain_grp = ds_read.getRoot()->createGroup(second_coll_name);

  MFEMSidreDataCollection first_sdc_reader(first_coll_name,
                                           first_bp_index_grp,
                                           first_domain_grp,
                                           owns_mesh_data);
  MFEMSidreDataCollection second_sdc_reader(second_coll_name,
                                            second_bp_index_grp,
                                            second_domain_grp,
                                            owns_mesh_data);

  // Needs to be set "manually" in order for everything to be loaded in properly
  first_sdc_reader.SetComm(MPI_COMM_WORLD);
  first_sdc_reader.Load();

  // Make sure that the useless view wasn't read back in
  EXPECT_FALSE(ds_read.getRoot()->hasView(useless_view_name));

  first_sdc_reader.SetGroupPointers(
    ds_read.getRoot()->getGroup(first_coll_name + "_global/blueprint_index/" + first_coll_name),
    ds_read.getRoot()->getGroup(first_coll_name));

  first_sdc_reader.UpdateStateFromDS();
  first_sdc_reader.UpdateMeshAndFieldsFromDS();

  second_sdc_reader.SetComm(MPI_COMM_WORLD);
  second_sdc_reader.Load();

  second_sdc_reader.SetGroupPointers(
    ds_read.getRoot()->getGroup(second_coll_name + "_global/blueprint_index/" + second_coll_name),
    ds_read.getRoot()->getGroup(second_coll_name));

  second_sdc_reader.UpdateStateFromDS();
  second_sdc_reader.UpdateMeshAndFieldsFromDS();

  auto first_gf_read = first_sdc_reader.GetField(field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(first_gf_read));

  auto second_gf_read = second_sdc_reader.GetField(field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(second_gf_read));

  // Make sure the gridfunction was actually read in
  EXPECT_LT(first_gf_read->ComputeL2Error(three_and_a_half), EPSILON);
  EXPECT_LT(second_gf_read->ComputeL2Error(five_and_a_half), EPSILON);

  EXPECT_TRUE(first_sdc_reader.verifyMeshBlueprint());
  EXPECT_TRUE(second_sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_par_reload_mesh_1D_small)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 1D mesh divided into segments
  auto mesh = mfem::Mesh::MakeCartesian1D(10);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_2D_small)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 2D mesh divided into triangles
  auto mesh = mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::TRIANGLE);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_2D_large)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 2D mesh divided into triangles
  auto mesh = mfem::Mesh::MakeCartesian2D(100, 100, mfem::Element::TRIANGLE);
  testParallelMeshReloadAllPartitionings(mesh);
}

  // The following test requires a function from mfem@4.3
  #if (MFEM_VERSION >= 40300)
TEST(sidre_datacollection, dc_par_reload_mesh_2D_periodic)
{
    #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
    #endif

  // periodic 2D mesh divided into triangles
  auto base_mesh =
    mfem::Mesh::MakeCartesian2D(10, 10, mfem::Element::Type::QUADRILATERAL, false, 1.0, 1.0);
  std::vector<mfem::Vector> translations = {mfem::Vector({1.0, 0.0}), mfem::Vector({0.0, 1.0})};
  auto vertex_map = base_mesh.CreatePeriodicVertexMapping(translations);
  auto mesh = mfem::Mesh::MakePeriodic(base_mesh, vertex_map);
  testParallelMeshReloadAllPartitionings(mesh);
}
  #endif

TEST(sidre_datacollection, dc_par_reload_mesh_3D_small_tet)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 3D mesh divided into tetrahedra
  auto mesh = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_medium_tet)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 3D mesh divided into tetrahedra
  auto mesh = mfem::Mesh::MakeCartesian3D(10, 10, 10, mfem::Element::TETRAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_small_hex)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 3D mesh divided into hexahedra
  auto mesh = mfem::Mesh::MakeCartesian3D(3, 3, 3, mfem::Element::HEXAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_medium_hex)
{
  #ifndef AXOM_USE_HDF5
  SUCCEED() << "sidre::MFEMSidreDataCollection::load(<cycle>) is only implemented "
               "for the 'sidre_hdf5' protocol";
  return;
  #endif

  // 3D mesh divided into hexahedra
  auto mesh = mfem::Mesh::MakeCartesian3D(10, 10, 10, mfem::Element::HEXAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}
#endif  // defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
#endif

  result = RUN_ALL_TESTS();

#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return result;
}
