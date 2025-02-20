#include "libcmaes/cmaes.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <libcmaes/cmastrategy.h>
#include <libcmaes/esoptimizer.h>
#include <libcmaes/llogging.h>

#include "dp_su2f2to2/avgsvcalculator.hpp"
#include "dp_su2f2to2/boltzmann.hpp"
#include "dp_su2f2to2/config.hpp"
#include "dp_su2f2to2/correspondance.hpp"
#include "dp_su2f2to2/leshouchesfrommarty.hpp"
#include "dp_su2f2to2/process.hpp"

using namespace su2f2to2;

using namespace advmath;
using namespace libcmaes;

const std::string PATHPLOTS = std::string(OUTPATH);
constexpr const double omega_h2_target = 0.12;
constexpr const double omega_h2_err = 0.01;

static inline double sq(const double x) { return x * x; }

int main(int argc, char* argv[])
{
  int read_args = 0;
  std::cout << "Running " << argv[read_args++] << '\n';

  constexpr const double sigma = omega_h2_err;
  constexpr const int npar = 1;
  constexpr const int n_required_args = 9;
  constexpr const char* par_names[npar] = {"g_chi"}; // Name of the parameters to be optimised
  constexpr const double m_phi_over_mv = 2.;

  [[maybe_unused]] constexpr const double init_seed_g_chi = 0.3;

  double lbounds[npar] = {1.e-6}, ubounds[npar] = {1.0}; // arrays for lower and upper parameter bounds, respectively

  std::cout << "This function requires " << n_required_args - 1 << " arguments:\n"
            << " - The name of the input file\n"
            << " - The name of the output file\n"
            << " - The seed for the parameter g_f\n"
            << " - The value for the parameter m_chi/m_V\n"
            << " - The initial value for m_V\n"
            << " - The final value for m_V\n"
            << " - The total number of points\n"
            << " - The current point";
  if (argc < n_required_args)
  {
    return 1;
  }

  Param_t input(argv[read_args]);
  std::cout << "The input file is " << argv[read_args++] << '\n';

  std::ofstream outfile{argv[read_args]};
  if (!outfile)
  {
    std::cout << "Impossible to open " << argv[read_args] << '\n';
  }
  std::cout << "Output file: " << argv[read_args++] << '\n';

  char* end;

  // Initialising x0
  std::vector<double> x0;
  x0.reserve(npar);
  x0.emplace_back(std::strtod(argv[read_args++], &end));

  std::cout << "SEED: ";
  for (size_t i = 0; i < x0.size(); i++)
    std::cout << par_names[i] << '=' << x0[i] << ' ';
  std::cout << '\n';

  // Reading parameters from CLI arguments
  const double mchi_over_mv = std::strtod(argv[read_args++], &end);

  const real_t m_V_init = std::strtod(argv[read_args++], &end);
  const real_t m_V_final = std::strtod(argv[read_args++], &end);
  const real_t npoints = std::strtod(argv[read_args++], &end);
  const real_t exponent = std::strtod(argv[read_args++], &end);

  const real_t m_V = m_V_init * std::pow(m_V_final / m_V_init, exponent / npoints);

  std::cout << "m_V=" << m_V << '\n';

  // Function and paramters' first definition
  input.g_f = x0[0];
  input.m_V = m_V;
  input.m_phi = 2. * m_V;
  input.m_chi_dm_1 = mchi_over_mv * m_V;
  input.m_chi_dm_2 = mchi_over_mv * m_V;
  input.refresh();
  BoltzmannSolver boltz(input);

  auto function_to_minimize = [&](const double* parameters, [[maybe_unused]] const int size_array = 2)
  {
    // parameters in input will be
    input.g_f = parameters[0];
    input.refresh();
    boltz.changeInput(input);
    const double omega = boltz.relic_density();
    return std::abs(omega - omega_h2_target);
  };

  // Minimization algorithm

  GenoPheno<pwqBoundStrategy> gp(lbounds, ubounds, npar); // genotype / phenotype transform associated to bounds.

  // -1 for automatically decided lambda, 0 is for random	seeding	of the internal generator.
  CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(x0, sigma, -1, 0, gp);
  cmaparams.set_algo(aCMAES);
  FitFunc f = function_to_minimize;
  CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(f, cmaparams);

  std::cout << "seed is\n";
  for (int i = 0; i < npar; i++)
    std::cout << par_names[i] << " = " << x0[i] << '\n';
  std::cout << "\nbest solution: ";
  cmasols.print(std::cout, 0, gp);
  std::cout << '\n';
  std::cout << "Expected Distance from Minimum: " << cmasols.edm() << '\n';
  std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0 << " seconds\n ";

  Candidate bcand = cmasols.get_best_seen_candidate();

  double fmin = bcand.get_fvalue(); // min objective function value the optimizer converged to
  Eigen::VectorXd bestparameters_dvec =
      gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec()); // Applying pheno to the best candidate output
  double edm = cmasols.edm();                                   // expected distance to the minimum.
  int return_val = cmasols.run_status();

  std::cout << "# m_phi/m_V  m_chi/m_V m_V   g_f   pull-Oh2   EDM  exit_code  \nRESULT=";
  std::cout << m_phi_over_mv << '\t' << mchi_over_mv << '\t' << input.m_V.get() << '\t' << bestparameters_dvec[0]
            << '\t' << fmin << '\t' << edm << '\t' << return_val << '\n';

  if(outfile)
  {
    outfile << "# m_phi/m_V  m_chi/m_V m_V  g_f   pull-Oh2   EDM  exit_code \n";
    outfile << m_phi_over_mv << '\t' << mchi_over_mv << '\t' << input.m_V.get() << '\t' << bestparameters_dvec[0] << '\t'
            << fmin << '\t' << edm << '\t' << return_val << '\n';
  }


  return return_val;
}