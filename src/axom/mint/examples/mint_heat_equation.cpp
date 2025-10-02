// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom Includes
#include "axom/core.hpp"
#include "axom/mint.hpp"
#include "axom/slic.hpp"

// C/C++ includes
#include <cmath>   /* for std::exp, std::ciel */
#include <sstream> /* for std::stringstream */
#include <string>  /* for std::string */

namespace axom
{
namespace mint
{
class Gaussian2D
{
public:
  /*!
   * \brief Creates a 2D Guassian with the given amplitude, mean, and
   * covariance.
   * \param [in] amplitude the maximum amplitude of the gaussian \f$ a \f$.
   * \param [in] mean the mean of the gaussian. The format is
   *  \f$ \mu_x, \mu_y \f$.
   * \param [in] covar the covariance of the gaussian. The format is
   *  \f$ \text{var}(x), \text{var}(y), \text{var}(xy) \f$.
   */
  Gaussian2D(double amplitude, const double mean[2], const double covar[3])
  {
    setMean(mean);
    setCovar(covar);
    setAmplitude(amplitude);
  }

  /*!
   * \brief Sets the mean.
   * \param [in] mean the mean of the gaussian. The format is
   *  \f$ \mu_x, \mu_y \f$.
   */
  void setMean(const double mean[2])
  {
    m_mean[0] = mean[0];
    m_mean[1] = mean[1];
  }

  /*!
   * \brief Sets the covariance.
   * \param [in] covar the covariance of the gaussian. The format is
   *  \f$ \text{var}(x), \text{var}(y), \text{var}(xy) \f$.
   * \note The covariance matrix needs to be positive definite.
   * \pre \f$ \text{var}(x) > 0 \f$
   * \pre \f$ \text{var}(y) > 0 \f$
   * \pre \f$ \text{var}(xy)^2 < \text{var}(x) \text{var}(y) \f$
   */
  void setCovar(const double covar[3])
  {
    double det = covar[0] * covar[1] - covar[2] * covar[2];
    if(covar[0] <= 0 || covar[1] <= 0 || det == 0)
    {
      SLIC_ERROR("Invalid covariance: (" << covar[0] << ", " << covar[1] << ", " << covar[2]
                                         << ") must be positive definite.");
    }

    m_half_covar_inv[0] = covar[1] / det / 2.0;
    m_half_covar_inv[1] = covar[0] / det / 2.0;
    m_half_covar_inv[2] = -covar[2] / det / 2.0;
  }

  /*!
   * \brief Sets the amplitude.
   * \param [in] amplitude the maximum amplitude of the gaussian \f$ a \f$.
   */
  void setAmplitude(double amplitude) { m_amplitude = amplitude; }

  /*!
   * \brief Evaluate the gaussian at a point.
   * \param [in] r the position vector at which the gaussian is evaluated.
   * \return \f$ f\left(\vec{r}\right) = ae^{-\frac{1}{2} \left(\vec{r} -
      \vec{\mu}\right)^\top \boldsymbol{\Sigma}^{-1}\left(\vec{r} -
      \vec{\mu}\right)} \f$.
   */
  double evaluate(const double r[2]) const
  {
    double dx = r[0] - m_mean[0];
    double dy = r[1] - m_mean[1];
    double temp = m_half_covar_inv[0] * dx * dx;
    temp += m_half_covar_inv[1] * dy * dy;
    temp += 2 * m_half_covar_inv[2] * dx * dy;
    return m_amplitude * std::exp(-temp);
  }

private:
  double m_mean[2];
  double m_half_covar_inv[3];
  double m_amplitude;
};

class HeatEquationSolver
{
public:
  /*!
   * \brief Sets up a heat equation solver on a uniform mesh.
   * \param [in] h the spacing of the uniform mesh.
   * \param [in] lower_bound the bottom left corner of the bounding box.
   * \param [in] upper_bound the upper right corner of the bounding box.
   * \note the equation we are solving is \f$ \frac{\partial U}{\partial t}
                                          = \alpha \nabla^2 U \f$.
   */
  HeatEquationSolver(double h, const double lower_bound[2], const double upper_bound[2])
    : m_mesh(create_mesh(h, lower_bound, upper_bound))
    , m_h(h)
  { }

  /*!
   * \brief Destroys the heat equation solver by deleting the mesh.
   */
  ~HeatEquationSolver() { delete m_mesh; }

  /*!
   * \brief Set the uniform mesh upon which to solve.
   * \param [in] h the spacing of the uniform mesh.
   * \param [in] lower_bound the bottom left corner of the bounding box.
   * \param [in] upper_bound the upper right corner of the bounding box.
   */
  void setMesh(const double h, const double lower_bound[2], const double upper_bound[2])
  {
    delete m_mesh;
    m_mesh = create_mesh(h, lower_bound, upper_bound);
    m_h = h;
  }

  /*!
   * \brief Apply a two dimensional gaussian to the mesh as an initial
   * condition.
   * \param [in] pulse the pulse to apply.
   */
  void initialize(const Gaussian2D& pulse)
  {
    const double* origin = m_mesh->getOrigin();
    IndexType Ni = m_mesh->getNodeResolution(mint::I_DIRECTION);
    IndexType Nj = m_mesh->getNodeResolution(mint::J_DIRECTION);
    double* t = m_mesh->getFieldPtr<double>("temperature", mint::NODE_CENTERED);

    IndexType idx = 0;
    double node_pos[2] = {origin[0], origin[1]};
    for(IndexType j = 0; j < Nj; ++j)
    {
      for(IndexType i = 0; i < Ni; ++i)
      {
        t[idx++] = pulse.evaluate(node_pos);
        node_pos[0] += m_h;
      }
      node_pos[0] = origin[0];
      node_pos[1] += m_h;
    }
  }

  /*!
   * \brief Solve the heat equation on the uniform mesh.
   * \param [in] alpha the conductivity.
   * \param [in] dt the time step.
   * \param [in] dt the simulation time.
   * \param [in] period the number of cycles between dumps.
   * \param [in] path the base path of the dump files.
   */
  void solve(double alpha, double dt, double t_max, int period, const std::string& path)
  {
    const IndexType num_nodes = m_mesh->getNumberOfNodes();
    double* new_temp = new double[num_nodes];
    double* prev_temp = m_mesh->getFieldPtr<double>("temperature", mint::NODE_CENTERED);

    /* Copy the boundary conditions into new_temp since they won't be copied
       during the time step. */
    copy_boundary(prev_temp, new_temp);

    int cur_dump = 0;
    double cur_time = 0.0;
    const int num_cycles = static_cast<int>(std::ceil(t_max / dt));
    for(int cycle = 0; cycle < num_cycles; ++cycle)
    {
      if(cycle == num_cycles - 1)
      {
        SLIC_INFO("Cycle: " << cycle << " Time: " << cur_time << "\n");
      }
      else
      {
        SLIC_INFO("Cycle: " << cycle << " Time: " << cur_time << "\r");
      }
      slic::flushStreams();

      if(cycle % period == 0)
      {
        write_dump(path, cur_dump++);
      }

      step(alpha, dt, prev_temp, new_temp);
      std::memcpy(prev_temp, new_temp, num_nodes * sizeof(double));

      if(cur_time + dt > t_max)
      {
        dt = t_max - cur_time;
      }
      cur_time += dt;
    }

    delete[] new_temp;
    write_dump(path, cur_dump);
    SLIC_INFO("Finished\n");
  }

private:
  /*!
   * \brief Copy the values on the boundary of the mesh.
   * \param [in] prev_temp the data to copy.
   * \param [in] new_temp the buffer to copy into.
   */
  void copy_boundary(const double* prev_temp, double* new_temp)
  {
    IndexType Ni = m_mesh->getNodeResolution(mint::I_DIRECTION);
    IndexType Nj = m_mesh->getNodeResolution(mint::J_DIRECTION);

    /* Copy the -y side, which is contiguous. */
    const IndexType memcpy_size = Ni * sizeof(double);
    std::memcpy(new_temp, prev_temp, memcpy_size);

    /* Copy the +y side, which is contiguous. */
    const IndexType offset = (Nj - 1) * Ni;
    std::memcpy(new_temp + offset, prev_temp + offset, memcpy_size);

    /* Copy the -x and +x sides which aren't contiguous. */
    for(IndexType idx = Ni; idx < offset; idx += Ni)
    {
      new_temp[idx] = prev_temp[idx];
      new_temp[idx + Ni - 1] = prev_temp[idx + Ni - 1];
    }
  }

  /*!
   * \brief Time step the system.
   * \param [in] alpha the conductivity.
   * \param [in] dt the time step to take.
   * \param [in] prev_temp the temperature at the beginning of the time step.
   * \param [out] new_temp the temperature at the end of the time step.
   * \note The timestepping algorithm is an explicit euler scheme with second
   *  order central difference in space. The update equation is \f$
      T_{i, j}^{n + 1} = \left( 1 - 4 \frac{\alpha \Delta t}{h^2} \right)
      T_{i, j}^n + \frac{\alpha \Delta t}{h^2}  \left( T_{i - 1, j}^n +
      T_{i + 1, j}^n + T_{i, j - 1}^n + T_{i, j + 1}^n \right) \f$.
   */
  void step(double alpha, double dt, const double* prev_temp, double* new_temp)
  {
    const double neighbors_scale = dt * alpha / (m_h * m_h);
    const double self_scale = 1.0 - (4.0 * neighbors_scale);

    /* Since the boundary conditions are fixed we only need to iterate over
       the interior nodes. */
    const IndexType jp = m_mesh->nodeJp();
    const IndexType Ni = m_mesh->getNodeResolution(I_DIRECTION);
    const IndexType Nj = m_mesh->getNodeResolution(J_DIRECTION);
    for(IndexType j = 1; j < Nj; ++j)
    {
      for(IndexType i = 1; i < Ni; ++i)
      {
        const IndexType j_offset = j * jp;
        const IndexType idx = i + j_offset;
        const IndexType north = idx + jp;
        const IndexType south = idx - jp;
        const IndexType east = idx + 1;
        const IndexType west = idx - 1;
        const double neighbors_contrib =
          prev_temp[north] + prev_temp[east] + prev_temp[south] + prev_temp[west];

        new_temp[idx] = neighbors_scale * neighbors_contrib;
        new_temp[idx] += self_scale * prev_temp[idx];
      } /* loop over i */
    } /* loop over j */
  }

  /*!
   * \brief Write out the mesh to a vtk file.
   * \param [in] path the base path of the file to write.
   * \param [in] cur_dump the number of previous dumps.
   */
  void write_dump(const std::string& path, int cur_dump)
  {
    std::stringstream cur_path;
    cur_path << path << "_" << cur_dump << ".vtk";
    if(write_vtk(m_mesh, cur_path.str()) != 0)
    {
      SLIC_WARNING("Unable to write to file: " << cur_path.str() << std::endl);
    }
  }

  /*!
   * \brief Set the uniform mesh upon which to solve.
   * \param [in] h the spacing of the uniform mesh.
   * \param [in] lower_bound the bottom left corner of the bounding box.
   * \param [in] upper_bound the upper right corner of the bounding box.
   * \return a pointer to the new uniform mesh.
   */
  static UniformMesh* create_mesh(const double h,
                                  const double lower_bound[2],
                                  const double upper_bound[2])
  {
    IndexType Ni, Nj;

    const double lx = axom::utilities::abs(upper_bound[0] - lower_bound[0]);
    Ni = static_cast<IndexType>((lx / h) + 1);

    const double ly = axom::utilities::abs(upper_bound[1] - lower_bound[1]);
    Nj = static_cast<IndexType>((ly / h) + 1);

    UniformMesh* mesh = new UniformMesh(lower_bound, upper_bound, Ni, Nj);
    mesh->createField<double>("temperature", mint::NODE_CENTERED);
    return mesh;
  }

  UniformMesh* m_mesh;
  double m_h;
};

} /* end namespace mint */
} /* end namespace axom */

const std::string help_string =
  "\nUsage: ./mint_heat_equation_ex [options]\n"
  "-h, -help\n"
  "\tPrint out this message then exit.\n"
  "-p -path PATH\n"
  "\tThe base bath of the dump files. The files will be written to\n"
  "\tpath_#.vtk where # is the dump number.\n"
  "-s, -spacing FLOAT\n"
  "\tSet the mesh spacing, must be greater than 0.\n."
  "-b FLOAT FLOAT FLOAT FLOAT, -bounds FLOAT FLOAT FLOAT FLOAT\n"
  "\tSet the bounding box for the mesh. Format is lower_x lower_y upper_x\n"
  "\tupper_y.\n"
  "-a FLOAT, -amplitude FLOAT\n"
  "\tSet the amplitude of the gaussian pulse.\n"
  "-m FLOAT FLOAT, -mean FLOAT FLOAT\n"
  "\tSet the mean of the gaussian pulse. Format is x y.\n"
  "-c FLOAT FLOAT FLOAT, -covariance FLOAT FLOAT FLOAT\n"
  "\tSet the covariance of the gaussian. Format is var_x var_y var_xy.\n"
  "-alpha FLOAT\n"
  "\tThe conductivity.\n"
  "-dt FLOAT\n"
  "\tThe time step, must be greater than zero.\n"
  "-t FLOAT\n"
  "\tThe time to run the simulation, must be greater than zero.\n"
  "-d INT, -dumpPeriod INT\n"
  "\tThe number of cycles to wait between writing a dump file.\n";

/// A structure that holds the command line arguments.
struct Arguments
{
  double h;
  double lower_bound[2];
  double upper_bound[2];
  double amplitude;
  double mean[2];
  double covar[3];
  double alpha;
  double dt;
  double t_max;
  int period;
  std::string path;
};

/*!
 * \brief Parses and validates the command line arguments.
 * \param [out] args where the settings are stored.
 * \param [in] argc the number of arguments.
 * \param [in] argv the array of arguments.
 */
void parse_arguments(Arguments& args, int argc, const char** argv)
{
  for(int i = 1; i < argc; ++i)
  {
    if(std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "-help") == 0)
    {
      SLIC_INFO(help_string);
      axom::utilities::processAbort();
    }
    else if(std::strcmp(argv[i], "-p") == 0 || std::strcmp(argv[i], "-path") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.path = argv[i + 1];
      i++;
    }
    else if(std::strcmp(argv[i], "-s") == 0 || std::strcmp(argv[i], "-spacing") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.h = atof(argv[i + 1]);
      i++;
    }
    else if(std::strcmp(argv[i], "-b") == 0 || std::strcmp(argv[i], "-bounds") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 4, "Not enough arguments.");
      args.lower_bound[0] = atof(argv[i + 1]);
      args.lower_bound[1] = atof(argv[i + 2]);
      args.upper_bound[0] = atof(argv[i + 3]);
      args.upper_bound[1] = atof(argv[i + 4]);
      i += 4;
    }
    else if(std::strcmp(argv[i], "-a") == 0 || std::strcmp(argv[i], "-amplitude") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.amplitude = atof(argv[i + 1]);
      i++;
    }
    else if(std::strcmp(argv[i], "-m") == 0 || std::strcmp(argv[i], "-mean") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 2, "Not enough arguments.");
      args.mean[0] = atof(argv[i + 1]);
      args.mean[1] = atof(argv[i + 2]);
      i += 2;
    }
    else if(std::strcmp(argv[i], "-c") == 0 || std::strcmp(argv[i], "-covariance") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 3, "Not enough arguments.");
      args.covar[0] = atof(argv[i + 1]);
      args.covar[1] = atof(argv[i + 2]);
      args.covar[2] = atof(argv[i + 3]);
      i += 3;
    }
    else if(std::strcmp(argv[i], "-alpha") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.alpha = atof(argv[i + 1]);
      i++;
    }
    else if(std::strcmp(argv[i], "-dt") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.dt = atof(argv[i + 1]);
      i++;
    }
    else if(std::strcmp(argv[i], "-t") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.t_max = atof(argv[i + 1]);
      i++;
    }
    else if(std::strcmp(argv[i], "-d") == 0 || std::strcmp(argv[i], "-dumpPeriod") == 0)
    {
      SLIC_ERROR_IF(i >= argc - 1, "Not enough arguments.");
      args.period = atoi(argv[i + 1]);
      i++;
    }
    else
    {
      SLIC_ERROR("Unrecognized argument.");
    }
  }

  /* Validate arguments. */
  if(args.lower_bound[0] >= args.upper_bound[0] || args.lower_bound[1] >= args.upper_bound[1])
  {
    SLIC_ERROR("Invalid bounding box: (" << args.lower_bound[0] << ", " << args.lower_bound[1]
                                         << ") x (" << args.upper_bound[0] << ", "
                                         << args.upper_bound[1] << ")");
  }
  if(args.covar[0] <= 0 || args.covar[1] <= 0 ||
     args.covar[2] * args.covar[2] >= args.covar[0] * args.covar[1])
  {
    SLIC_ERROR("Invalid covariance: (" << args.covar[0] << ", " << args.covar[1] << ", "
                                       << args.covar[2] << ") must be positive definite.");
  }

  SLIC_ERROR_IF(args.h <= 0, "Invalid spacing: " << args.h);
  SLIC_ERROR_IF(args.dt <= 0, "Invalid time step: " << args.dt);
  SLIC_ERROR_IF(args.t_max <= 0, "Invalid run time: " << args.t_max);
  SLIC_ERROR_IF(args.period < 0, "Invalid dump period: " << args.period);

  std::string dir;
  axom::utilities::filesystem::getDirName(dir, args.path);
  if(axom::utilities::filesystem::makeDirsForPath(dir) != 0)
  {
    SLIC_ERROR("Could not make directories for Dump path: " << args.path);
  }

  double dt_max = args.h * args.h / (4 * args.alpha);
  if(args.dt >= dt_max)
  {
    SLIC_WARNING("The chosen time step " << args.dt << " is larger than " << dt_max
                                         << " this will lead to numerical instability.\n");
  }

  SLIC_INFO("Mesh bounds: (" << args.lower_bound[0] << ", " << args.lower_bound[1] << ") x ("
                             << args.upper_bound[0] << ", " << args.upper_bound[1] << ")\n");
  SLIC_INFO("Mesh spacing: " << args.h << std::endl);
  SLIC_INFO("Gaussian amplitude: " << args.amplitude << std::endl);
  SLIC_INFO("Gaussian mean: (" << args.mean[0] << ", " << args.mean[1] << ")\n");
  SLIC_INFO("Gaussian covariance: (" << args.covar[0] << ", " << args.covar[1] << ", "
                                     << args.covar[2] << ")\n");
  SLIC_INFO("Conductivity: " << args.alpha << std::endl);
  SLIC_INFO("Time step: " << args.dt << std::endl);
  SLIC_INFO("Max stable time step: " << dt_max << std::endl);
  SLIC_INFO("CFL Number: " << args.dt / dt_max << std::endl);
  SLIC_INFO("Run time: " << args.t_max << std::endl);
  SLIC_INFO("Number of cycles: " << std::ceil(args.t_max / args.dt) << std::endl);
  SLIC_INFO("Dump period: " << args.period << std::endl);
  SLIC_INFO("Dump path: " + args.path << std::endl);
}

/*!
 * \brief Initializes SLIC.
 */
void init()
{
  axom::slic::initialize();
  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);

  std::string slicFormatStr = "[<LEVEL>] <MESSAGE>";
  axom::slic::GenericOutputStream* defaultStream = new axom::slic::GenericOutputStream(&std::cout);
  axom::slic::GenericOutputStream* compactStream =
    new axom::slic::GenericOutputStream(&std::cout, slicFormatStr);
  axom::slic::addStreamToMsgLevel(defaultStream, axom::slic::message::Error);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Warning);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Info);
  axom::slic::addStreamToMsgLevel(compactStream, axom::slic::message::Debug);
}

/*!
 * \brief Finalizes SLIC.
 */
void finalize() { axom::slic::finalize(); }

int main(int argc, const char** argv)
{
  init();

  /* Default settings */
  Arguments args;
  args.h = 0.25;
  args.lower_bound[0] = -10.0;
  args.lower_bound[1] = -10.0;
  args.upper_bound[0] = 10.0;
  args.upper_bound[1] = 10.0;
  args.amplitude = 1.0;
  args.mean[0] = 0.0;
  args.mean[1] = 0.0;
  args.covar[0] = 1.0;
  args.covar[1] = 1.0;
  args.covar[2] = 0.7;
  args.alpha = 1.0;
  args.dt = 0.01;
  args.t_max = 60.0;
  args.period = 60;
  args.path = "./results/dump";

  /* Read in settings from the command line. */
  parse_arguments(args, argc, argv);

  /* Solve. */
  axom::mint::HeatEquationSolver solver(args.h, args.lower_bound, args.upper_bound);
  axom::mint::Gaussian2D pulse(args.amplitude, args.mean, args.covar);
  SLIC_INFO("Initializing\n");
  solver.initialize(pulse);
  SLIC_INFO("Solving\n");
  solver.solve(args.alpha, args.dt, args.t_max, args.period, args.path);

  finalize();
  return 0;
}
