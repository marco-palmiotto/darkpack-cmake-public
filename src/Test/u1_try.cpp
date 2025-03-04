#include <filesystem> // Include the filesystem library
#include <fstream>
#include <iomanip>
#include <iostream>

#include "dp_u1f2to2/avgsvcalculator.hpp"
#include "dp_u1f2to2/boltzmann.hpp"
#include "dp_u1f2to2/config.hpp"
#include "dp_u1f2to2/correspondance.hpp"
#include "dp_u1f2to2/leshouchesfrommarty.hpp"
#include "dp_u1f2to2/process.hpp"

using namespace u1f2to2;
using namespace advmath;
namespace fs = std::filesystem; // Create an alias for the filesystem namespace

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
  constexpr const int n_required_args = 2;
  // constexpr const char* par_names[npar] = {"g_chi"}; // Name of the parameters to be optimised
  constexpr const double m_phi_over_mv = 2.;

  [[maybe_unused]] constexpr const double init_seed_g_chi = 0.3;

  double lbounds[npar] = {1.e-6}, ubounds[npar] = {1.0}; // arrays for lower and upper parameter bounds, respectively

  std::cout << "This function requires " << n_required_args - 1 << " arguments:\n"
            << " - The name of the input file\n";
  if (argc < n_required_args)
  {
    return 1;
  }

  std::string input_file = argv[read_args++];
  if (!fs::exists(input_file)) // Check if the file exists
  {
    std::cerr << "Error: The input file " << input_file << " does not exist.\n";
    return 1;
  }

  Param_t input(input_file);
  std::cout << "The input file is " << input_file << '\n';

  input.m_phi = input.m_V_3 * m_phi_over_mv;
  input.refresh();

  input.Print();

  BoltzmannSolver boltzmann(input);

  std::cout << "Oh2=" << boltzmann.relic_density() << '\n';
  return 0;
}