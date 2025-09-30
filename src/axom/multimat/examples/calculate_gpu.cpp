// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
//
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/multimat/multimat.hpp"

#define ITERMAX 1  // we don't use this, but helper.hpp does
#include "helper.hpp"

#include "axom/CLI11.hpp"

#include <map>

using namespace std::literals;

namespace slic = axom::slic;
namespace mmat = axom::multimat;

enum class RuntimePolicy
{
  seq = 0,
  raja_seq = 1,
  raja_omp = 2,
  raja_cuda = 3,
  raja_hip = 4
};

const std::vector<std::string> policy_strs = {"seq", "raja_seq", "raja_omp", "raja_cuda"};

struct Input
{
  RuntimePolicy policy {RuntimePolicy::seq};
  int num_cells {10000};
  int num_mats {50};

  const static std::map<std::string, RuntimePolicy> s_validPolicies;

  void parse(int argc, char** argv, axom::CLI::App& app);
};

const std::map<std::string, RuntimePolicy> Input::s_validPolicies(
  {{"seq", RuntimePolicy::seq}
#ifdef AXOM_USE_RAJA
   ,
   {"raja_seq", RuntimePolicy::raja_seq}
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
   ,
   {"raja_omp", RuntimePolicy::raja_omp}
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
   ,
   {"raja_cuda", RuntimePolicy::raja_cuda}
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
   ,
   {"raja_hip", RuntimePolicy::raja_hip}
  #endif
#endif
  });

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  std::stringstream pol_sstr;
  pol_sstr << "With \'-m bvh\' or \'-m naive\', set runtime policy. \n"
           << "Set to \'seq\' or 0 to use the sequential algorithm "
           << "(w/o RAJA).";
#ifdef AXOM_USE_RAJA
  pol_sstr << "\nSet to \'raja_seq\' or 1 to use the RAJA sequential policy.";
  #if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  pol_sstr << "\nSet to \'raja_omp\' or 2 to use the RAJA OpenMP policy.";
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  pol_sstr << "\nSet to \'raja_cuda\' or 3 to use the RAJA CUDA policy.";
  #endif
  #if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  pol_sstr << "\nSet to \'raja_hip\' or 4 to use the RAJA HIP policy.";
  #endif
#endif

  app.add_option("-p, --policy", policy, pol_sstr.str())
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(Input::s_validPolicies));

  app.add_option("-c,--num-cells", num_cells, "Sets the number of cells to generate data for.")
    ->capture_default_str();

  app.add_option("-m,--num-mats", num_mats, "Sets the number of materials to generate data for.")
    ->capture_default_str();

  app.get_formatter()->column_width(76);

  // Could throw an exception
  app.parse(argc, argv);

  // Output parsed information
  SLIC_INFO("Using parameter values: " << "\n  execution policy = "s << policy_strs[(int)policy]);
}

int allocator_id;

template <typename T>
using SparseField2D = mmat::MultiMat::SparseField2D<T>;

template <typename ExecSpace>
void avgDensityCompactFlat(mmat::MultiMat& mm)
{
  axom::fmt::print("Running average density compact - flat\n");
  int ncells = mm.getNumberOfCells();

  auto* relationSet = mm.getSparse2dFieldSet(mmat::DataLayout::CELL_DOM);
  int nthreads = relationSet->totalSize();

  const auto density = mm.getSparse2dField<double>("Densityfrac");
  const auto vf = mm.getSparse2dField<double>("Volfrac");

  const auto vol = mm.get1dField<double>("Vol");

  axom::Array<double> densityAvg(ncells, ncells, allocator_id);
  const auto densityAvg_view = densityAvg.view();

  axom::for_all<ExecSpace>(
    nthreads,
    AXOM_LAMBDA(int flatid) {
      int cell_id = relationSet->flatToFirstIndex(flatid);

      double density_avg_slot = density[flatid] * vf[flatid];

      axom::atomicAdd<ExecSpace>(&densityAvg_view[cell_id], density_avg_slot / vol[cell_id]);
    });
}

template <typename ExecSpace>
void avgDensityCompactSubmap(mmat::MultiMat& mm)
{
  axom::fmt::print("Running average density compact - submap\n");
  int ncells = mm.getNumberOfCells();

  const auto density = mm.getSparse2dField<double>("Densityfrac");
  const auto vf = mm.getSparse2dField<double>("Volfrac");

  const auto vol = mm.get1dField<double>("Vol");

  axom::Array<double> densityAvg(ncells, ncells, allocator_id);
  const auto densityAvg_view = densityAvg.view();

  axom::for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(int cellid) {
      double density_avg = 0.0;
      auto density_row = density(cellid);
      auto volfrac_row = vf(cellid);

      for(int slotid = 0; slotid < volfrac_row.size(); slotid++)
      {
        density_avg += density_row(slotid) * volfrac_row(slotid);
      }
      densityAvg_view[cellid] = density_avg / vol[cellid];
    });
}

template <typename ExecSpace>
void avgDensityDirect(mmat::MultiMat& mm)
{
  axom::fmt::print("Running average density - direct\n");
  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();

  const auto density = mm.getDense2dField<double>("Densityfrac");
  const auto vf = mm.getDense2dField<double>("Volfrac");

  const auto vol = mm.get1dField<double>("Vol");

  axom::Array<double> densityAvg(ncells, ncells, allocator_id);
  const auto densityAvg_view = densityAvg.view();

  axom::for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(int cellid) {
      double density_avg = 0.0;

      for(int matid = 0; matid < nmats; matid++)
      {
        density_avg += density(cellid, matid) * vf(cellid, matid);
      }
      densityAvg_view[cellid] = density_avg / vol[cellid];
    });
}

template <typename ExecSpace>
void avgDensitySubmap(mmat::MultiMat& mm)
{
  axom::fmt::print("Running average density - submap\n");
  int ncells = mm.getNumberOfCells();
  int nmats = mm.getNumberOfMaterials();

  const auto density = mm.getDense2dField<double>("Densityfrac");
  const auto vf = mm.getDense2dField<double>("Volfrac");

  const auto vol = mm.get1dField<double>("Vol");

  axom::Array<double> densityAvg(ncells, ncells, allocator_id);
  const auto densityAvg_view = densityAvg.view();

  axom::for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(int cellid) {
      double density_avg = 0.0;
      auto density_row = density(cellid);
      auto volfrac_row = vf(cellid);

      for(int matid = 0; matid < nmats; matid++)
      {
        density_avg += density_row(matid) * volfrac_row(matid);
      }
      densityAvg_view[cellid] = density_avg / vol[cellid];
    });
}

template <typename ExecSpace>
void avgDensityIter(mmat::MultiMat& mm)
{
  axom::fmt::print("Running average density - iterators\n");
  int ncells = mm.getNumberOfCells();

  const auto density = mm.getDense2dField<double>("Densityfrac");
  const auto vf = mm.getDense2dField<double>("Volfrac");

  const auto vol = mm.get1dField<double>("Vol");

  axom::Array<double> densityAvg(ncells, ncells, allocator_id);
  const auto densityAvg_view = densityAvg.view();

  axom::for_all<ExecSpace>(
    ncells,
    AXOM_LAMBDA(int cellid) {
      double density_avg = 0.0;

      auto density_it = density.begin(cellid);
      auto volfrac_it = vf.begin(cellid);

      auto density_end = density.end(cellid);

      while(density_it != density_end)
      {
        density_avg += (*density_it) * (*volfrac_it);
        ++density_it;
        ++volfrac_it;
      }

      densityAvg_view[cellid] = density_avg / vol[cellid];
    });
}

template <typename ExecSpace>
void traverseCells(mmat::MultiMat& mm)
{
  allocator_id = axom::execution_space<ExecSpace>::allocatorID();
  mm.setAllocatorID(allocator_id);

  avgDensityDirect<ExecSpace>(mm);
  avgDensitySubmap<ExecSpace>(mm);
  avgDensityIter<ExecSpace>(mm);

  mm.convertLayoutToSparse();

  avgDensityCompactFlat<ExecSpace>(mm);
  avgDensityCompactSubmap<ExecSpace>(mm);
}

int main(int argc, char** argv)
{
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);
  slic::addStreamToAllMsgLevels(new slic::GenericOutputStream(&std::cout));

  SLIC_INFO("Axom Version:" << " [" << axom::getVersion() << "]");

  // Parse the command line arguments
  Input params;
  axom::CLI::App app {"Multimat GPU calculate example"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  Robey_data data("", params.num_cells, params.num_mats);
  data.set_up_cell_dom_data();

  mmat::MultiMat mm {mmat::DataLayout::CELL_DOM, mmat::SparsityLayout::DENSE};
  mm.setNumberOfMaterials(data.nmats);
  mm.setNumberOfCells(data.ncells);
  mm.setCellMatRel(data.Volfrac_bool, mmat::DataLayout::CELL_DOM);

  //Setting field data in terms of slam
  mm.addField<>("Densityfrac",
                mmat::FieldMapping::PER_CELL_MAT,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.Densityfrac.view());
  mm.addField<>("Vol",
                mmat::FieldMapping::PER_CELL,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.Vol.view());
  mm.addField<>("Volfrac",
                mmat::FieldMapping::PER_CELL_MAT,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.Volfrac.view());
  mm.addField<>("Tempfrac",
                mmat::FieldMapping::PER_CELL_MAT,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.Temperaturefrac.view());
  mm.addField<>("Pressurefrac",
                mmat::FieldMapping::PER_CELL_MAT,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.Pressurefrac.view());
  mm.addField<>("nmatconsts",
                mmat::FieldMapping::PER_MAT,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.nmatconsts.view());
  mm.addField<>("MatDensityAverage",
                mmat::FieldMapping::PER_CELL_MAT,
                mmat::DataLayout::CELL_DOM,
                mmat::SparsityLayout::DENSE,
                data.Pressurefrac.view());

  //printself and check
  mm.isValid(true);

  switch(params.policy)
  {
  case RuntimePolicy::seq:
  case RuntimePolicy::raja_seq:
    traverseCells<axom::SEQ_EXEC>(mm);
    break;
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  case RuntimePolicy::raja_omp:
    traverseCells<axom::OMP_EXEC>(mm);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  case RuntimePolicy::raja_cuda:
    traverseCells<axom::CUDA_EXEC<256>>(mm);
    break;
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  case RuntimePolicy::raja_hip:
    traverseCells<axom::HIP_EXEC<256>>(mm);
    break;
#endif
  default:
    SLIC_ERROR("Unhandled runtime policy case");
    break;
  }

  slic::finalize();

  return 0;
}
