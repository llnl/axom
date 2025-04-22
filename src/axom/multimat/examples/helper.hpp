// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * Set-up multi-material data for examples
 *
 * Set-up code from LANL's MultiMatTest repo:
 *     https://github.com/lanl/MultiMatTest
 * Also defines some helper struct-classes.
 */

#include "axom/core.hpp"
#include "axom/slam.hpp"
#include "axom/fmt.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

namespace slam = axom::slam;

struct Value_Checker
{
  std::vector<double> values;
  void reset() { values.resize(0); }
  void check(std::vector<double>& vec, double EPS = 1e-8)
  {
    if(values.empty())
    {
      values = vec;
    }
    else
    {
      if(values.size() != vec.size())
      {
        SLIC_ERROR(
          axom::fmt::format("Sizes of arrays are different. 'values' has {} "
                            "elements; 'vec' has {} elements",
                            values.size(),
                            vec.size()));
      }
      else
      {
        using axom::utilities::isNearlyEqual;
        bool equivalent = true;
        const int SZ = values.size();
        for(int i = 0; i < SZ; ++i)
        {
          if(!isNearlyEqual(values[i], vec[i], EPS))
          {
            equivalent = false;
          }
        }

        SLIC_ERROR_IF(!equivalent, "Calculated values are not the same!");
      }
    }
  }
  void check(axom::ArrayView<double> vec, double EPS = 1e-8)
  {
    if(values.empty())
    {
      values = std::vector<double>(vec.begin(), vec.end());
    }
    else
    {
      if(static_cast<axom::IndexType>(values.size()) != vec.size())
      {
        SLIC_ERROR(
          axom::fmt::format("Sizes of arrays are different. 'values' has {} "
                            "elements; 'vec' has {} elements",
                            values.size(),
                            vec.size()));
      }
      else
      {
        using axom::utilities::isNearlyEqual;
        bool equivalent = true;
        const int SZ = values.size();
        for(int i = 0; i < SZ; ++i)
        {
          if(!isNearlyEqual(values[i], vec[i], EPS))
          {
            equivalent = false;
          }
        }

        SLIC_ERROR_IF(!equivalent, "Calculated values are not the same!");
      }
    }
  }
};

struct multirun_timer
{
  std::vector<double> time_record;
  axom::utilities::Timer timer;
  double time_sum;

  void reset()
  {
    time_sum = 0;
    time_record.resize(0);
    timer.reset();
  }
  void start() { timer.start(); }

  void record()
  {
    timer.stop();
    double elapsed = timer.elapsed();
    time_sum += elapsed;
    time_record.push_back(elapsed);
  }

  /// Returns the median time
  double get_median()
  {
    auto sz = time_record.size();
    if(sz == 0)
    {
      return 0.;
    }

    std::vector<double> sorted(time_record);
    std::sort(sorted.begin(), sorted.end());

    return (sz % 2) == 0 ? 0.5 * (sorted[(sz >> 1) - 1] + sorted[sz >> 1]) : sorted[sz >> 1];
  }

  double get_average() { return time_sum / (double)time_record.size(); }
};

void make_other_field_data_celldom(int ncells,
                                   int nmats,
                                   std::vector<double>& i_Volfrac_CD,
                                   axom::Array<double>& o_Volfrac,
                                   axom::Array<double>& o_Vol,
                                   axom::Array<double>& o_Densityfrac,
                                   axom::Array<double>& o_Temperaturefrac,
                                   axom::Array<double>& o_Pressurefrac,
                                   std::vector<double>& o_Volfrac_sparse,
                                   std::vector<double>& o_Densityfrac_sparse,
                                   std::vector<double>& o_Temperaturefrac_sparse,
                                   std::vector<double>& o_Pressurefrac_sparse,
                                   std::vector<int>& o_begin_idx,
                                   std::vector<int>& o_col_idx)
{
  //Set up dense data
  o_Vol.resize(ncells, 0);
  o_Volfrac.resize(ncells * nmats, 0);
  o_Densityfrac.resize(ncells * nmats, 0);
  o_Temperaturefrac.resize(ncells * nmats, 0);
  o_Pressurefrac.resize(ncells * nmats, 0);
  int nnz = 0;
  for(int ic = 0; ic < ncells; ic++)
  {
    o_Vol[ic] = 0.0;
    for(int m = 0; m < nmats; m++)
    {
      o_Volfrac[ic * nmats + m] = i_Volfrac_CD[ic * nmats + m];
      if(i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        nnz += 1;
        o_Densityfrac[ic * nmats + m] = 2.0;
        o_Temperaturefrac[ic * nmats + m] = 0.5;
      }
      else
      {
        o_Densityfrac[ic * nmats + m] = 0.0;
        o_Temperaturefrac[ic * nmats + m] = 0.0;
      }
      o_Pressurefrac[ic * nmats + m] = 0.0;

      o_Vol[ic] += i_Volfrac_CD[ic * nmats + m];
    }
  }

  //fill in sparse data
  o_begin_idx.resize(ncells + 1);
  o_col_idx.resize(nnz);
  int ii = 0;
  for(int ic = 0; ic < ncells; ic++)
  {
    o_begin_idx[ic] = ii;
    for(int m = 0; m < nmats; m++)
    {
      if(i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        o_col_idx[ii] = m;
        ii++;
      }
    }
  }
  o_begin_idx[ncells] = ii;
  assert(ii == nnz);

  o_Densityfrac_sparse.resize(nnz);
  o_Volfrac_sparse.resize(nnz);
  o_Temperaturefrac_sparse.resize(nnz);
  o_Pressurefrac_sparse.resize(nnz);

  ii = 0;
  for(int ci = 0; ci < ncells; ci++)
  {
    for(int mi = 0; mi < nmats; mi++)
    {
      if(o_Volfrac[ci * nmats + mi] != 0)
      {
        o_Densityfrac_sparse[ii] = o_Densityfrac[ci * nmats + mi];
        o_Volfrac_sparse[ii] = o_Volfrac[ci * nmats + mi];
        o_Temperaturefrac_sparse[ii] = o_Temperaturefrac[ci * nmats + mi];
        o_Pressurefrac_sparse[ii] = o_Pressurefrac[ci * nmats + mi];
        ii++;
      }
    }
  }
  assert(ii == nnz);
}

void make_other_field_data_matdom(int ncells,
                                  int nmats,
                                  std::vector<double>& i_Volfrac_CD,
                                  axom::Array<double>& o_Volfrac,
                                  axom::Array<double>& o_Vol,
                                  axom::Array<double>& o_Densityfrac,
                                  axom::Array<double>& o_Temperaturefrac,
                                  axom::Array<double>& o_Pressurefrac,
                                  std::vector<double>& o_Volfrac_sparse,
                                  std::vector<double>& o_Densityfrac_sparse,
                                  std::vector<double>& o_Temperaturefrac_sparse,
                                  std::vector<double>& o_Pressurefrac_sparse,
                                  std::vector<int>& o_begin_idx,
                                  std::vector<int>& o_col_idx,
                                  std::vector<int>& o_dense2sparse_idx

)
{
  //Set up dense data
  o_Vol.resize(ncells, 0);
  o_Volfrac.resize(ncells * nmats, 0);
  o_Densityfrac.resize(ncells * nmats, 0);
  o_Temperaturefrac.resize(ncells * nmats, 0);
  o_Pressurefrac.resize(ncells * nmats, 0);
  int nnz = 0;
  for(int ic = 0; ic < ncells; ic++)
  {
    o_Vol[ic] = 0.0;
  }

  for(int m = 0; m < nmats; m++)
  {
    for(int ic = 0; ic < ncells; ic++)
    {
      o_Volfrac[m * ncells + ic] = i_Volfrac_CD[ic * nmats + m];
      if(i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        nnz += 1;
        o_Densityfrac[m * ncells + ic] = 2.0;
        o_Temperaturefrac[m * ncells + ic] = 0.5;
      }
      else
      {
        o_Densityfrac[m * ncells + ic] = 0;
        o_Temperaturefrac[m * ncells + ic] = 0;
      }
      o_Pressurefrac[m * ncells + ic] = 0.0;
      o_Vol[ic] += i_Volfrac_CD[ic * nmats + m];
    }
  }

  //fill in sparse data, CSR layout
  o_begin_idx.resize(nmats + 1);
  o_col_idx.resize(nnz);
  o_dense2sparse_idx.resize(nmats * ncells, -1);
  int ii = 0;
  for(int m = 0; m < nmats; m++)
  {
    o_begin_idx[m] = ii;
    int spidx = 0;
    for(int ic = 0; ic < ncells; ic++)
    {
      if(i_Volfrac_CD[ic * nmats + m] > 0.0)
      {
        o_col_idx[ii] = ic;
        o_dense2sparse_idx[m * ncells + ic] = spidx;
        ii++;
        spidx++;
      }
    }
  }
  o_begin_idx[nmats] = ii;
  assert(ii == nnz);

  o_Densityfrac_sparse.resize(nnz);
  o_Volfrac_sparse.resize(nnz);
  o_Temperaturefrac_sparse.resize(nnz);
  o_Pressurefrac_sparse.resize(nnz);

  ii = 0;
  for(int mi = 0; mi < nmats; mi++)
  {
    for(int ci = 0; ci < ncells; ci++)
    {
      if(o_Volfrac[mi * ncells + ci] != 0)
      {
        o_Densityfrac_sparse[ii] = o_Densityfrac[mi * ncells + ci];
        o_Volfrac_sparse[ii] = o_Volfrac[mi * ncells + ci];
        o_Temperaturefrac_sparse[ii] = o_Temperaturefrac[mi * ncells + ci];
        o_Pressurefrac_sparse[ii] = o_Pressurefrac[mi * ncells + ci];
        ii++;
      }
    }
  }
  assert(ii == nnz);
}

template <typename T1, typename T2>
float ratio_percent(T1 numer, T2 denom)
{
  return static_cast<float>(numer * 100.0) / static_cast<float>(denom);
}

void read_vol_frac_matrix_file(std::string filename,
                               int& ncells,
                               int& nmats,
                               std::vector<double>& Volfrac,
                               float& filled_percentage)
{
  using namespace std;

  ncells = 1000000;

  int status;
  FILE* fp;
  fp = fopen(filename.c_str(), "r");
  if(!fp)
  {
    fprintf(stderr, "unable to read volume fractions from file \"%s\"\n", filename.c_str());
    exit(-1);
  }

  status = fscanf(fp, "%d", &nmats);
  if(status < 0)
  {
    printf("error in read at line %d\n", __LINE__);
    exit(1);
  }

  Volfrac.resize(ncells * nmats, 0.0);

  char matname[256];
  for(int m = 0; m < nmats; m++)
  {
    status = fscanf(fp, "%s", matname);  // read and discard
    if(status < 0)
    {
      printf("error in read at line %d\n", __LINE__);
      exit(1);
    }
  }

  double VolTotal = 0.0;
  int filled_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  int pure_frac_count = 0;
  int pure_cell_count = 0;
  int onematcell = 0;
  int twomatcell = 0;
  int threematcell = 0;
  int fourmatcell = 0;
  int fiveplusmatcell = 0;
  for(int ic = 0; ic < ncells; ic++)
  {
    int mat_count = 0;
    for(int m = 0; m < nmats; m++)
    {
      status = fscanf(fp, "%lf", &(Volfrac[ic * nmats + m]));
      if(status < 0)
      {
        printf("error in read at line %d\n", __LINE__);
        exit(1);
      }
      if(Volfrac[ic * nmats + m] > 0.0)
      {
        filled_count++;
        mat_count++;
      }
      VolTotal += Volfrac[ic * nmats + m];
    }
    if(mat_count >= 2)
    {
      mixed_cell_count++;
      mixed_frac_count += mat_count;
    }
    else
    {
      pure_frac_count++;
    }
    if(mat_count == 1)
    {
      pure_cell_count++;
    }
    if(mat_count == 1)
    {
      onematcell++;
    }
    if(mat_count == 2)
    {
      twomatcell++;
    }
    if(mat_count == 3)
    {
      threematcell++;
    }
    if(mat_count == 4)
    {
      fourmatcell++;
    }
    if(mat_count >= 5)
    {
      fiveplusmatcell++;
    }
  }
  fclose(fp);

  printf("Ratios to Full Data Structure\n");
  filled_percentage = ratio_percent(filled_count, ncells * nmats);
  float sparsity_percentage = ratio_percent(ncells * nmats - filled_count, ncells * nmats);
  printf("Sparsity %lf percent/Filled %lf percent\n\n", sparsity_percentage, filled_percentage);

  printf("Ratios to Number of Cells\n");
  float pure_cell_percentage = ratio_percent(pure_cell_count, ncells);
  float mixed_cell_percentage = ratio_percent(mixed_cell_count, ncells);
  printf("Pure cell %lf percentage/Mixed material %lf percentage\n\n",
         pure_cell_percentage,
         mixed_cell_percentage);

  printf("Ratios to Mixed Material Data Structure\n");
  float mixed_material_sparsity_percentage =
    ratio_percent(mixed_frac_count, mixed_cell_count * nmats);
  float mixed_material_filled_percentage =
    ratio_percent(mixed_cell_count * nmats - mixed_frac_count, mixed_cell_count * nmats);
  printf("Mixed material Sparsity %lf percent/Mixed material Filled %lf percent\n\n",
         mixed_material_sparsity_percentage,
         mixed_material_filled_percentage);

  printf("Vol Total %lf\n", VolTotal);
  printf("%f percent of the cells are filled\n", ratio_percent(filled_count, ncells * nmats));
  printf("%f percent of the cells are mixed\n", ratio_percent(mixed_cell_count, ncells));
  printf("%f percent of the total are mixed\n", ratio_percent(mixed_frac_count, ncells * nmats));
  printf("%f percent of the frac are mixed\n",
         ratio_percent(mixed_frac_count, mixed_cell_count * nmats));
  printf("%f percent sparsity\n", ratio_percent(ncells * nmats - mixed_frac_count, ncells * nmats));
  printf("%f percent of the frac are pure\n", ratio_percent(pure_frac_count, ncells));
  printf("1 matcell %d 2 matcell %d 3 matcell %d 4 matcell %d 5 matcell %d\n\n",
         onematcell,
         twomatcell,
         threematcell,
         fourmatcell,
         fiveplusmatcell);
  printf("Total cells %d\n\n",
         onematcell + 2 * twomatcell + 3 * threematcell + 4 * fourmatcell + 5 * fiveplusmatcell);
}

void get_vol_frac_matrix_rand(int& ncells,
                              int& nmats,
                              std::vector<double>& Volfrac,
                              float& filled_percentage,
                              int default_ncell = 1000000,
                              int default_nmats = 50)
{
  //void get_vol_frac_matrix_rand(double **&Volfrac, float& filled_percentage) {
  using namespace std;

  ncells = default_ncell;
  nmats = default_nmats;

  Volfrac.resize(ncells * nmats, 0.0);

  vector<int> mf_rand(ncells, 0);
  srand(0);
  for(int ic = 0; ic < ncells; ic++)
  {
    mf_rand[ic] = (int)((float)rand() * 1000.0 / (float)((long long)RAND_MAX + 1));
  }

  double VolTotal = 0.0;
  int filled_count = 0;
  int mixed_cell_count = 0;
  int mixed_frac_count = 0;
  int pure_frac_count = 0;
  int pure_cell_count = 0;
  int onematcell = 0;
  int twomatcell = 0;
  int threematcell = 0;
  int fourmatcell = 0;
  int fiveplusmatcell = 0;
  for(int ic = 0; ic < ncells; ic++)
  {
    int m1 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
    if(m1 > nmats - 1)
    {
      m1 = nmats - 1;
    }
    Volfrac[ic * nmats + m1] = 1.0;
    int mf = mf_rand[ic];
    if(mf < 25)
    {
      int m2 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
      while(m2 == m1)
      {
        m2 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
      }
      int m3 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
      while(m3 == m2 || m3 == m1)
      {
        m3 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
      }
      int m4 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
      while(m4 == m3 || m4 == m2 || m4 == m1)
      {
        m4 = (int)((float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1));
      }
      if(m2 > nmats - 1)
      {
        m2 = nmats - 1;
      }
      if(m3 > nmats - 1)
      {
        m3 = nmats - 1;
      }
      if(m4 > nmats - 1)
      {
        m4 = nmats - 1;
      }
      Volfrac[ic * nmats + m1] = 0.4;
      Volfrac[ic * nmats + m2] = 0.3;
      Volfrac[ic * nmats + m3] = 0.2;
      Volfrac[ic * nmats + m4] = 0.1;
    }
    else if(mf < 75)
    {
      int m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while(m2 == m1)
      {
        m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      int m3 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while(m3 == m2 || m3 == m1)
      {
        m3 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      if(m2 > nmats - 1)
      {
        m2 = nmats - 1;
      }
      if(m3 > nmats - 1)
      {
        m3 = nmats - 1;
      }
      Volfrac[ic * nmats + m1] = 0.5;
      Volfrac[ic * nmats + m2] = 0.3;
      Volfrac[ic * nmats + m3] = 0.2;
    }
    else if(mf < 200)
    {
      int m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      while(m2 == m1)
      {
        m2 = (float)rand() * (float)nmats / (float)((long long)RAND_MAX + 1);
      }
      if(m2 > nmats - 1)
      {
        m2 = nmats - 1;
      }
      Volfrac[ic * nmats + m1] = 0.5;
      Volfrac[ic * nmats + m2] = 0.5;
    }
    int mat_count = 0;
    for(int m = 0; m < nmats; m++)
    {
      if(Volfrac[ic * nmats + m] > 0.0)
      {
        filled_count++;
        mat_count++;
      }
      VolTotal += Volfrac[ic * nmats + m];
    }
    if(mat_count >= 2)
    {
      mixed_cell_count++;
      mixed_frac_count += mat_count;
    }
    else
    {
      pure_frac_count++;
    }
    if(mat_count == 1)
    {
      pure_cell_count++;
    }
    if(mat_count == 1)
    {
      onematcell++;
    }
    if(mat_count == 2)
    {
      twomatcell++;
    }
    if(mat_count == 3)
    {
      threematcell++;
    }
    if(mat_count == 4)
    {
      fourmatcell++;
    }
    if(mat_count >= 5)
    {
      fiveplusmatcell++;
    }
  }

  printf("Ratios to Full Data Structure\n");
  filled_percentage = (float)filled_count * 100.0 / (float)(ncells * nmats);
  float sparsity_percentage =
    (float)(ncells * nmats - filled_count) * 100.0 / (float)(ncells * nmats);
  printf("Sparsity %lf percent/Filled %lf percent\n\n", sparsity_percentage, filled_percentage);

  printf("Ratios to Number of Cells\n");
  float pure_cell_percentage = (float)pure_cell_count * 100.0 / (float)ncells;
  float mixed_cell_percentage = (float)mixed_cell_count * 100.0 / (float)ncells;
  printf("Pure cell %lf percentage/Mixed material %lf percentage\n\n",
         pure_cell_percentage,
         mixed_cell_percentage);

  printf("Ratios to Mixed Material Data Structure\n");
  float mixed_material_sparsity_percentage =
    (float)mixed_frac_count * 100.0 / (float)(mixed_cell_count * nmats);
  float mixed_material_filled_percentage =
    (float)(mixed_cell_count * nmats - mixed_frac_count) * 100.0 / (float)(mixed_cell_count * nmats);
  printf("Mixed material Sparsity %lf percent/Mixed material Filled %lf percent\n\n",
         mixed_material_sparsity_percentage,
         mixed_material_filled_percentage);

  printf("Vol Total %lf\n", VolTotal);
  printf("%f percent of the cells are filled\n",
         (float)filled_count * 100.0 / (float)(ncells * nmats));
  printf("%f percent of the cells are mixed\n", (float)mixed_cell_count * 100.0 / (float)ncells);
  printf("%f percent of the total are mixed\n",
         (float)mixed_frac_count * 100.0 / (float)(ncells * nmats));
  printf("%f percent of the frac are mixed\n",
         (float)mixed_frac_count * 100.0 / (float)(mixed_cell_count * nmats));
  printf("%f percent sparsity\n",
         (float)(ncells * nmats - mixed_frac_count) * 100.0 / (float)(ncells * nmats));
  printf("%f percent of the frac are pure\n", (float)pure_frac_count * 100.0 / (float)ncells);
  printf("1 matcell %d 2 matcell %d 3 matcell %d 4 matcell %d 5 matcell %d\n\n",
         onematcell,
         twomatcell,
         threematcell,
         fourmatcell,
         fiveplusmatcell);
  printf("Total cells %d\n\n",
         onematcell + 2 * twomatcell + 3 * threematcell + 4 * fourmatcell + 5 * fiveplusmatcell);
}

void get_neighbors(int ncells, int nnbrs_max, std::vector<int>& num_nbrs, std::vector<int>& nbrs)
{
  int ncells1 = (int)sqrt(ncells);  // assumes ncells is a perfect square
  if(ncells1 * ncells1 != ncells)
  {
    fprintf(stderr, "Number of cells in mesh is not a perfect square");
    exit(-1);
  }
  if(nnbrs_max != 8)
  {
    std::cout << "This neighbor generating code assumes 2D 8-neighbors. "
              << "You need to update this." << std::endl;
    exit(-1);
  }

  for(int i = 0; i < ncells1; i++)
  {
    for(int j = 0; j < ncells1; j++)
    {
      int c = i * ncells1 + j;
      int ilo = i == 0 ? i : i - 1;
      int jlo = j == 0 ? j : j - 1;
      int ihi = i == ncells1 - 1 ? i : i + 1;
      int jhi = j == ncells1 - 1 ? j : j + 1;
      int n = 0;
      for(int i1 = ilo; i1 <= ihi; i1++)
      {
        for(int j1 = jlo; j1 <= jhi; j1++)
        {
          int c2 = i1 * ncells1 + j1;
          if(c2 != c)
          {
            nbrs[c * nnbrs_max + n] = i1 * ncells1 + j1;
            n++;
          }
        }
      }
      num_nbrs[c] = n;
    }
  }
}

void get_centroids(int ncells, std::vector<double>& cen)
{
  int ncells1 = (int)sqrt(ncells);  // assumes ncells is a perfect square
  if(ncells1 * ncells1 != ncells)
  {
    fprintf(stderr, "Number of cells in mesh is not a perfect square");
    exit(-1);
  }

  // Assume domain is a unit square

  double XLO = 0.0, YLO = 0.0, XHI = 1.0, YHI = 1.0;
  double dx = (XHI - XLO) / ncells1, dy = (YHI - YLO) / ncells1;

  for(int i = 0; i < ncells1; i++)
  {
    for(int j = 0; j < ncells1; j++)
    {
      int c = i * ncells1 + j;
      cen[c * 2] = XLO + i * dx;
      cen[c * 2 + 1] = YLO + j * dy;
    }
  }
}

float filled_fraction = -1;

struct Robey_data
{
  /*
  * Data structure to store Robey's data.
  * When constructed, the volume fraction data is created.
  * Then a layout is selected for the other array to be filled in.
  */

  int ncells;
  int nmats;
  float filled_percentage;

  //For creating MultiMat object. Always in dense cell-dominant layout.
  std::vector<bool> Volfrac_bool;
  std::vector<double> Volfrac_CD;  //cell-dominant full volfrac array

  axom::Array<double> Vol;  //per cell, for all layouts

  //Per cellmat
  axom::Array<double> Volfrac;
  axom::Array<double> Densityfrac;
  axom::Array<double> Temperaturefrac;
  axom::Array<double> Pressurefrac;

  int cellmatcount;

  axom::Array<double> nmatconsts;

  int nnbrs_max;            //max number of neighbor = 8 for a 2d structured mesh
  std::vector<int> nnbrs;   //number of neighbors
  std::vector<int> nbrs;    //neighbor element id
  std::vector<double> cen;  //centroids of cells

  std::vector<int> subset2mesh;
  std::vector<int> mesh2subset;
  std::vector<int> nmatscell;  //number of materials in a cell
  std::vector<int> matids;     //material id in a cell, at most 4 materials per cell
  std::vector<int> ncellsmat;  //number of cells a material is in
  std::vector<int> dense2sparse_idx;

  //For CSR layout
  std::vector<int> begin_idx;
  std::vector<int> col_idx;
  std::vector<double> Volfrac_sparse;
  std::vector<double> Densityfrac_sparse;
  std::vector<double> Temperaturefrac_sparse;
  std::vector<double> Pressurefrac_sparse;

  // For slam's sets, relations and maps
  using P = int;
  using ElemSet = slam::PositionSet<P, P>;
  using Ind = slam::policies::STLVectorIndirection<P, P>;
  using Card = slam::policies::VariableCardinality<P, Ind>;
  using NbrRel = slam::StaticRelation<P, P, Card, Ind, ElemSet, ElemSet>;
  using DataInd = slam::policies::STLVectorIndirection<P, double>;
  using CentroidMap = slam::Map<double, ElemSet, DataInd, slam::policies::CompileTimeStride<P, 2>>;

  std::vector<int> slam_nbr_beg;
  std::vector<int> slam_nbr_ind;

  ElemSet slam_elems;
  NbrRel slam_neighbors;
  CentroidMap slam_centroids;

  Robey_data(std::string filename = "", int ncells_in = 100, int nmats_in = 50)
  {
    if(filename != "")
    {
      //read from file... large and takes a long time.
      read_vol_frac_matrix_file(filename, ncells, nmats, Volfrac_CD, filled_percentage);
    }
    else
    {
      //create random data
      //get_vol_frac_matrix_rand(ncells, nmats, Volfrac, filled_percentage);
      get_vol_frac_matrix_rand(ncells,
                               nmats,
                               Volfrac_CD,
                               filled_percentage,
                               ncells_in,
                               nmats_in);  //small version
    }

    filled_fraction = filled_percentage / 100.0f;

    // Some variables on neighbors
    //float L_f = read_from_file_bool ? 0.5 : 1.0;
    //// ave frac of nbrs containing material
    //int nnbrs_ave = 8;  // nearly so; 4000 boundary cells in 1 million cells
    //                    // in 3D, nnbrs_ave would be 26
    nnbrs_max = 8;
    // Build up list of neighbors for each cell
    // Assuming a 2D structured mesh, each cell will have a maximum of 8 nbrs
    nnbrs.resize(ncells);
    nbrs.resize(ncells * nnbrs_max);
    cen.resize(ncells * 2);

    nmatconsts.resize(nmats, 5.0);

    get_neighbors(ncells, nnbrs_max, nnbrs, nbrs);

    // Compute centroids of cells
    get_centroids(ncells, cen);

    //Making data for SLAM cell to mat relation
    Volfrac_bool.resize(ncells * nmats, false);
    cellmatcount = 0;
    for(unsigned int i = 0; i < Volfrac_CD.size(); i++)
    {
      if(Volfrac_CD[i] > 0)
      {
        Volfrac_bool[i] = true;
        cellmatcount++;
      }
    }

    /// Set up slam neighbor relation
    {
      // Initialize the elems set
      slam_elems = ElemSet(ncells);

      // Initialize the begin indices for the nbrs relation
      //via prefix sum on nnbrs array
      slam_nbr_beg.resize(ncells + 1);
      slam_nbr_beg[0] = 0;
      for(auto i : slam_elems.positions())
      {
        slam_nbr_beg[i + 1] = slam_nbr_beg[i] + nnbrs[i];
      }

      // Initialize the indices for the nbrs relation
      // Note -- we're compacting an array of 8 neighbors per element
      //         to one that has the correct size
      const int MAX_NBRS = nnbrs_max;
      slam_nbr_ind.resize(slam_nbr_beg[ncells]);
      int cur = 0;
      for(auto i : slam_elems.positions())
      {
        for(int j = 0; j < nnbrs[i]; ++j)
        {
          slam_nbr_ind[cur++] = nbrs[i * MAX_NBRS + j];
        }
      }

      using NBuilder = typename NbrRel::RelationBuilder;
      using BuildBeg = typename NBuilder::BeginsSetBuilder;
      using BuildInd = typename NBuilder::IndicesSetBuilder;

      // Initialize the nbrs relation
      slam_neighbors = NBuilder()
                         .fromSet(&slam_elems)
                         .toSet(&slam_elems)
                         .begins(BuildBeg().size(slam_nbr_beg.size()).data(&slam_nbr_beg))
                         .indices(BuildInd().size(slam_nbr_ind.size()).data(&slam_nbr_ind));

      // Initialize map over centroids
      slam_centroids = CentroidMap(&slam_elems);
      std::copy(cen.begin(), cen.end(), slam_centroids.data().begin());
    }

  }  //end constructor

  void set_up_cell_dom_data()
  {
    make_other_field_data_celldom(ncells,
                                  nmats,
                                  Volfrac_CD,
                                  Volfrac,
                                  Vol,
                                  Densityfrac,
                                  Temperaturefrac,
                                  Pressurefrac,
                                  Volfrac_sparse,
                                  Densityfrac_sparse,
                                  Temperaturefrac_sparse,
                                  Pressurefrac_sparse,
                                  begin_idx,
                                  col_idx);
  }
  void set_up_mat_dom_data()
  {
    make_other_field_data_matdom(ncells,
                                 nmats,
                                 Volfrac_CD,
                                 Volfrac,
                                 Vol,
                                 Densityfrac,
                                 Temperaturefrac,
                                 Pressurefrac,
                                 Volfrac_sparse,
                                 Densityfrac_sparse,
                                 Temperaturefrac_sparse,
                                 Pressurefrac_sparse,
                                 begin_idx,
                                 col_idx,
                                 dense2sparse_idx);
  }
};

struct Result_Store
{
  using DataLayout = axom::multimat::DataLayout;
  using SparsityLayout = axom::multimat::SparsityLayout;

  const int nAlgo = 3;
  enum Algo
  {
    avg_density,
    neighbor_density,
    pressure_calc
  };
  const char* algo_names[3] = {"Avg density",
                               "Neighbor material density",
                               "Pressure from ideal gas law"};

  static constexpr int nMethod = 15;
  enum Method
  {
    method_csr,
    mm_direct,
    mm_direct_templated_bset,
    mm_direct_templated_full,
    mm_direct_slam,
    mm_direct_slam_tmpl,
    mm_direct_slam_tmpl_stride,
    mm_idxarray,
    mm_submap,
    mm_submap_slam,
    mm_submap_slam_tmpl,
    mm_submap_templated_bset,
    mm_submap_templated_full,
    mm_iter,
    mm_flatiter
  };
  const char* method_names[nMethod] = {"CSR",
                                       "MM-Direct",
                                       "MM-Direct-BSet-Templated",
                                       "MM-Direct-Fully-Templated",
                                       "MM-Direct-Slam-BMap",
                                       "MM-Direct-Slam-BMap-Tmpl",
                                       "MM-Direct-Slam-BMap-Tmpl-StrideOne",
                                       "MM-Index Array",
                                       "MM-Submap",
                                       "MM-Submap-Slam",
                                       "MM-Submap-Slam-Tmpl",
                                       "MM-Submap-BSet-Templated",
                                       "MM-Submap-Fully-Templated",
                                       "MM-Iterator",
                                       "MM-Flat Iterator"};

  const int nLayout = 4;
  const char* const data_layout_str[2] = {"Cell Dominant", "Material Dominant"};
  const char* const sparsity_str[2] = {"Full Matrix", "Compact"};

  Robey_data* robey_data_ptr;
  std::vector<double> result_vec;

  Result_Store()
  {
    robey_data_ptr = nullptr;
    result_vec.resize(nAlgo * nLayout * nMethod, 0.0);
  }

  void init(Robey_data* robey_data_ptr_in) { robey_data_ptr = robey_data_ptr_in; }

  void add_result(Algo algo,
                  DataLayout data_layout,
                  axom::multimat::SparsityLayout sparsity_layout,
                  Method method,
                  double time)
  {
    int data_layout_i = data_layout == DataLayout::CELL_DOM ? 0 : 1;
    int sparsity_layout_i = sparsity_layout == SparsityLayout::DENSE ? 0 : 1;

    int idx = algo * (nLayout * nMethod) + data_layout_i * (nLayout / 2 * nMethod) +
      sparsity_layout_i * (nMethod) + method;

    std::cout << idx << ": " << get_algo_name(idx) << " - " << get_method_name(idx) << std::endl;
    result_vec[idx] = time;
  }

  const char* get_method_name(int index) { return method_names[index % nMethod]; }
  std::string get_algo_name(int index)
  {
    int algo_i = index / (nLayout * nMethod);
    int data_layout_i = (index / (nLayout / 2 * nMethod)) % 2;
    int sparsity_layout_i = (index / (nMethod)) % 2;
    return std::string(algo_names[algo_i]) + "|" + std::string(sparsity_str[sparsity_layout_i]) +
      "|" + std::string(data_layout_str[data_layout_i]);
  }

  void save_to_csv_file(const std::string& filename)
  {
    std::ofstream outputFile;
    outputFile.open(filename.c_str());
    if(!outputFile)
    {
      std::cout << "Writing result to file failed." << std::endl;
      return;
    }

    outputFile << "NCells: " << robey_data_ptr->ncells << " NMats: " << robey_data_ptr->nmats
               << " Sparsity: " << robey_data_ptr->filled_percentage << " NRuns: " << ITERMAX
               << "\n\n";

    outputFile << "Methods";
    for(int i = 0; i < nMethod; i++)
    {
      outputFile << "," << method_names[i];
    }
    outputFile << "\n";

    for(int i = 0; i < (int)result_vec.size() / nMethod; i++)
    {
      int idx = i * nMethod;

      //skip this row if no data is taken
      bool all_zero = true;
      for(int j = 0; j < nMethod; j++)
      {
        if(result_vec[idx + j] != 0.0)
        {
          all_zero = false;
        }
      }
      if(all_zero)
      {
        continue;
      }

      outputFile << get_algo_name(idx) << ",";
      for(int j = 0; j < nMethod; j++)
      {
        outputFile << result_vec[idx + j] << ",";
      }
      outputFile << "\n";
    }

    outputFile.close();

    std::cout << "Result written to csv file." << std::endl;
  }
};
