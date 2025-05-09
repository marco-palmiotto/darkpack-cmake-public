#include <fstream>
#include <iomanip>
#include <iostream>

// #include <algorithm>
// #include <execution>
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

  std::vector<real_t> m_V_jumps = {10.5926, 10.7978 };
  std::vector<real_t> g_f_jumps = {0.0548977, 0.055293};

  const real_t slope_before_jump = (m_V_jumps[1] - m_V_jumps[0]) / (g_f_jumps[1] - g_f_jumps[0]);
  const std::vector<real_t> m_V_trial_vec({10.9, 10.91});

  auto g_f_trial = [=](const real_t m_V_trial) { return g_f_jumps[1] + 1. / slope_before_jump * (m_V_trial - m_V_jumps[1]); };

  for (auto& elem : m_V_trial_vec)
  {
    m_V_jumps.push_back(elem);
    g_f_jumps.push_back(g_f_trial(elem));
  }

  m_V_jumps.push_back( 11.007);
  g_f_jumps.push_back(0.0521988);   // Function and paramters' first definition
  {
    const real_t m_V = 10.9;
    input.g_f = 0.1;
    input.m_V_3 = m_V;
    input.m_phi = 2. * m_V;
    input.m_chi_dm_1 = mchi_over_mv * m_V;
    input.m_chi_dm_2 = mchi_over_mv * m_V;
  }

  input.refresh();
  size_t n = 0;

  for (auto [m_V, g_f] : std::views::zip(m_V_jumps, g_f_jumps))
  {
    input.m_V_3 = m_V;
    input.g_f = g_f;
    input.m_phi = m_phi_over_mv * input.m_V_3;
    input.m_chi_dm_1 = mchi_over_mv * input.m_V_3;
    input.m_chi_dm_2 = mchi_over_mv * input.m_V_3;
    input.refresh();
    BoltzmannSolver boltz(input);
    // boltz.full_comput = true;

    std::cout << "m_V = " << m_V << " g_f = " << g_f << " 2*m_chi = " << 2. * input.m_chi_dm_1.get() << '\n';
    const real_t relic_density = boltz.relic_density();
    // Create a file named Weff_<value_of_mv>.dat
    std::ostringstream filename;
    filename << "out/scans/results/u1-theta23=0/Weff_" << m_V << ".dat";
    std::ofstream outfile(filename.str());
    if (!outfile.is_open())
    {
      std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";
      continue;
    }

    // Write the header for the file
    outfile << "# peff sqrts weff\n";

    // Retrieve the contributing processes
    std::vector<const Process2to2*> list_allowed(boltz.compute_contributing_processes());
    if (n < list_allowed.size())
    {
      n = list_allowed.size();
      size_t i = 1;
      for (const Process2to2* proc : list_allowed)
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
    sigma_v_filename << "out/scans/results/u1-theta23=0/sigma_v_" << m_V << ".dat";
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

    constexpr const unsigned int n_points = 1000;
    const real_t x_factor = std::pow(x_fin * T_in, 1. / n_points);
    std::array<real_t, n_points + 1> abscissae_sigmav;

    for (auto i : std::views::iota(0u, n_points + 1))
    {
      abscissae_sigmav[i] = std::pow(x_factor, i) / T_in;
    }

    // const auto abscissae_view = std::ranges::iota_view(0u, n_points);

    // std::transform(std::execution::par, abscissae_view.begin(), abscissae_view.end(), abscissae_sigmav.begin(),
    //                [&](const auto i) { return std::pow(x_factor, i) / T_in; });

    for (const real_t& x : abscissae_sigmav)
    {
      const real_t T = 1. / x;
      const real_t sigma_v = boltz.getAverageSigmav(T);
      const real_t Y_eq = boltz.Yeq(T);
      const real_t sigmavnum = boltz.getAverageSigmav_coan_hightemp_num(T);
      const real_t sigmavden = boltz.getAverageSigmav_coan_hightemp_den(T);
      sigma_v_file << x * boltz.getMassLBSM() << " " << sigma_v << " " << Y_eq << " " << sigmavnum << " " << sigmavden
                   << "\n";
      if (std::abs((sigmavnum / sigmavden - sigma_v) / (sigmavnum / sigmavden + sigma_v)) > 1.0e-3)
      {
        std::cerr << "Warning: sigmavnum/sigmavden = " << sigmavnum / sigmavden << " and sigma_v = " << sigma_v
                  << " are not consistent\n";
      }
    }
    sigma_v_file.close();
    std::cout << "Sigma_v table written to " << sigma_v_filename.str() << "\n";

    std::cout << "Oh2=" << relic_density << " xfo=" << boltz.getMassLBSM() / Tfo << "\n";

    std::vector<std::vector<const Process2to2*>> grouped_processes;
    if (!list_allowed.empty())
    {
      for (const Process2to2* proc : list_allowed)
      {
        bool added_to_group = false;
        const real_t threshold_proc = proc->compute_final_state_treshold(input);

        // Check if the process can be added to an existing group
        for (auto& group : grouped_processes)
        {
          const real_t threshold_group = group.front()->compute_final_state_treshold(input);

          if ((threshold_proc == 0. && threshold_group == 0.) ||
              std::abs(threshold_proc - threshold_group) / std::abs(threshold_proc + threshold_group) < 1.0e-3)
          {
            group.push_back(proc);
            added_to_group = true;
            break;
          }
        }

        // If not added to any group, create a new group
        if (!added_to_group)
        {
          grouped_processes.emplace_back(std::vector<const Process2to2*>{proc});
        }
      }
    }

    // Output the grouped processes
    unsigned short int group_index = 1;
    for (const auto& group : grouped_processes)
    {
      std::cout << "Group " << group_index++ << ":\n";
      for (const Process2to2* proc : group)
      {
        std::cout << "  - " << proc->getMname() << " (Threshold: " << proc->compute_final_state_treshold(input)
                  << ")\n";
      }
      
      std::shared_ptr< std::vector< Process2to2 > > group_ptr; 
      
      for( const Process2to2* elem : group )
      {
        if( group_ptr == nullptr )
          group_ptr = std::make_shared< std::vector< Process2to2 > >();
        group_ptr->push_back( *elem );
      }

      AvgSvCalculator avgsv(input, group_ptr);

      // Compute the label for the group based on the threshold of the first process
      const real_t threshold_group = group.front()->compute_final_state_treshold(input);
      std::ostringstream group_label_stream;
      group_label_stream << std::fixed << std::setprecision(3) << threshold_group;
      const std::string group_label = group_label_stream.str();

      // Create a file for the current group
      std::ostringstream group_filename;
      group_filename << "out/scans/results/u1-theta23=0/sigma_vg_" << m_V << "_group_" << group_label << ".dat";
      std::ofstream group_file(group_filename.str());
      if (!group_file.is_open())
      {
          std::cerr << "Error: Could not open file " << group_filename.str() << " for writing.\n";
          continue;
      }

      // Write the header for the file
      group_file << "# x sigmav\n";

      // Compute sigmav for each value in abscissae_sigmav
      real_t sigmav_total = 0.0;
      for (const real_t& x : abscissae_sigmav)
      {
          const real_t T = 1. / x;
        
          sigmav_total += avgsv.getAverageSigmav(T); 

          // Write x and sigmav_total to the file
          group_file << x << " " << sigmav_total << "\n";
      }

      group_file.close();
      std::cout << "Sigma_v table for group " << group_label << " written to " << group_filename.str() << "\n";

      }
  }

  return 0;
}