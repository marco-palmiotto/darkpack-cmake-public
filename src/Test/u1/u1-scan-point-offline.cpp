#include <fstream>
#include <iomanip>
#include <iostream>

#include <ranges>

#include "dp_u1f2to2/avgsvcalculator.hpp"
#include "dp_u1f2to2/boltzmann.hpp"
#include "dp_u1f2to2/config.hpp"
#include "dp_u1f2to2/correspondance.hpp"
#include "dp_u1f2to2/leshouchesfrommarty.hpp"
#include "dp_u1f2to2/process.hpp"

using namespace u1f2to2;

using namespace advmath;

static inline double sq(const double x) { return x * x; }

int main(int argc, char* argv[])
{
  int read_args = 0;
  std::cout << "Running " << argv[read_args++] << '\n';

  constexpr const int n_required_args = 2;
  constexpr const double m_phi_over_mv = 2.;

  std::cout << "This function requires " << n_required_args - 1 << " arguments:\n"
            << " - The name of the input file\n";
  if (argc < n_required_args)
  {
    return 1;
  }

  Param_t input(argv[read_args]);
  std::cout << "The input file is " << argv[read_args++] << '\n';

  constexpr const double mchi_over_mv = 0.33;

  // File path
  std::string file_path = "/workspaces/darkpack-cmake/out/scans/results/u1-theta23=0/growing_intervals.txt";

  std::vector<real_t> m_V_jumps;
  std::vector<real_t> g_f_jumps;

  // Open the file
  std::ifstream file(file_path);
  if (!file.is_open())
  {
    std::cerr << "Error: Could not open file " << file_path << std::endl;
    return 1;
  }

  // Read the file line by line
  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream iss(line);
    real_t x, y;

    // Parse x and y values from the line
    if (iss >> x >> y)
    {
      m_V_jumps.push_back(x);
      g_f_jumps.push_back(y);
    }
  }

  // Close the file
  file.close();

  // auto [a, b] = std::views::zip(m_V_jumps, g_f_jumps);

  // Function and paramters' first definition
  {
    const real_t m_V = 10.7978;
    input.g_f = 0.1;
    input.m_V_3 = m_V;
    input.m_phi = 2. * m_V;
    input.m_chi_dm_1 = mchi_over_mv * m_V;
    input.m_chi_dm_2 = mchi_over_mv * m_V;
  }

  input.refresh();
  AvgSvCalculator asv(input);
  size_t n = 0;

  for (auto m_V = m_V_jumps.begin(), g_f = g_f_jumps.begin(); m_V != m_V_jumps.end() || g_f != g_f_jumps.end();
       ++m_V, ++g_f)
  {

    input.m_V_3 = *m_V;
    input.g_f = *g_f;
    input.m_phi = m_phi_over_mv * input.m_V_3;
    input.m_chi_dm_1 = mchi_over_mv * input.m_V_3;
    input.m_chi_dm_2 = mchi_over_mv * input.m_V_3;
    input.refresh();
    asv.changeInput(input);
    std::cout << "m_V = " << *m_V << " g_f = " << *g_f << " 2*m_chi = " << 2. * input.m_chi_dm_1.get() << '\n';
    std::vector<const Process2to2*> list_allowed(asv.compute_contributing_processes());
    if (n < list_allowed.size())
    {
      n = list_allowed.size();
      size_t i = 1;
      for (const auto* proc : list_allowed)
      {
        std::cout << i++ << ": " << proc->getMname() << '\n';
      }
    }

    std::cout << '\n';
  }


  return 0;
}