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

  std::vector<real_t> m_V_jumps = {10.7978, 11.007};
  std::vector<real_t> g_f_jumps = {0.055293, 0.0521988};

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
    BoltzmannSolver boltz(input);

    std::cout << "m_V = " << *m_V << " g_f = " << *g_f << " 2*m_chi = " << 2. * input.m_chi_dm_1.get() << '\n';

    // Create a file named Weff_<value_of_mv>.dat
    std::ostringstream filename;
    filename << "out/scans/results/u1-theta23=0/Weff_" << *m_V << ".dat";
    std::ofstream outfile(filename.str());
    if (!outfile.is_open())
    {
      std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";
      continue;
    }

    // Write the header for the file
    outfile << "# peff sqrts weff\n";

    // Retrieve the contributing processes
    std::vector<const Process2to2*> list_allowed(boltz.get_contributing_processes());
    if (n < list_allowed.size())
    {
      n = list_allowed.size();
      size_t i = 1;
      for (const auto* proc : list_allowed)
      {
        std::cout << i++ << ": " << proc->getMname() << '\n';
      }
    }

    // Generate and write the peff, sqrts, and weff tables
    for (size_t i = 0; i < boltz.getWefftabsize(); ++i)
    {
      real_t peff, sqrts, weff;
      boltz.get_g2_WeffTabElement(i, peff, sqrts, weff);
      outfile << peff << " " << sqrts << " " << weff << "\n";
    }

    outfile.close();
    std::cout << "Weff table written to " << filename.str() << "\n";
    std::cout << '\n';

    // Create a file named sigma_v_<value_of_mv>.dat
    std::ostringstream sigma_v_filename;
    sigma_v_filename << "out/scans/results/u1-theta23=0/sigma_v_" << *m_V << ".dat";
    std::ofstream sigma_v_file(sigma_v_filename.str());
    if (!sigma_v_file.is_open())
    {
      std::cerr << "Error: Could not open file " << sigma_v_filename.str() << " for writing.\n";
      continue;
    }

    boltz.setTfo();
    const real_t Tfo = boltz.getTfo();
    const real_t T_in = 10. * Tfo;
    const real_t T_fin = 2.7;
    const real_t x_fin = 1. / T_fin;

    const unsigned int n_points = 1000;
    const real_t x_factor = std::pow(x_fin * T_in, 1. / n_points);
    real_t x = 1. / T_in;

    for (unsigned int i(0); i < n_points; i++)
    {
      x = std::pow(x_factor, i) / T_in;
      const real_t T = 1. / x;
      const real_t sigma_v = boltz.getAverageSigmav(T);
      const real_t Y_eq = boltz.Yeq(T);
      const real_t sigmavnum = boltz.getAverageSigmav_coan_hightemp_num(T);
      const real_t sigmavden = boltz.getAverageSigmav_coan_hightemp_den(T);
      sigma_v_file << x * boltz.getMassLBSM() << " " << sigma_v << " " << Y_eq << " " << sigmavnum << " " << sigmavden << "\n";
      if( std::abs( (sigmavnum/sigmavden - sigma_v)/(sigmavnum/sigmavden + sigma_v) ) > 1.0e-3)
      {
        std::cerr << "Warning: sigmavnum/sigmavden = " << sigmavnum / sigmavden
                  << " and sigma_v = " << sigma_v << " are not consistent\n";
      }
    }
    sigma_v_file.close();
    std::cout << "Sigma_v table written to " << sigma_v_filename.str() << "\n";
    const real_t relic_density = boltz.relic_density();
    std::cout << "Oh2=" << relic_density << " xfo=" << boltz.getMassLBSM() / Tfo << "\n";
  }

  return 0;
}