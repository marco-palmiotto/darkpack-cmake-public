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

  std::vector<real_t> m_V_jumps = {// 0.000510392, 0.000771106, 0.00106936, 0.00118252, 0.00126454,
                                   //  0.00135406,  0.00143791,
                                   0.00152342, 0.00161401, 77.9232};
  std::vector<real_t> g_f_jumps = {// 0.00588459, 0.00591574, 0.00639891, 0.00648346, 0.00650477,
                                   //  0.00649191, 0.00637008,
                                   0.00292935, 0.00594678, 0.12692};

  // size_t n = 0;

  for (auto [m_V, g_f] : std::views::zip(m_V_jumps, g_f_jumps))
  {
    input.m_V_3 = m_V;
    input.g_f = g_f;
    input.m_phi = m_phi_over_mv * input.m_V_3;
    input.m_chi_dm_1 = mchi_over_mv * input.m_V_3;
    input.refresh();
    BoltzmannSolver boltz(input);

    std::cout << "\nThe parameters are:\n"
              << "m_V = " << m_V << " g_f = " << g_f << " 2*m_chi = " << 2. * input.m_chi_dm_1.get() << '\n';

    // Create a file named Weff_<value_of_mv>.dat
    std::ostringstream filename;
    std::ofstream outfile;

    // Create a file named sigma_v_<value_of_mv>.dat
    filename.str("");
    filename.clear();
    filename << "out/scans/results/u1-theta23=0/sigma_v_" << m_V << ".dat";
    std::ofstream sigma_v_file(filename.str());
    if (!sigma_v_file.is_open())
    {
      std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";
      continue;
    }

    boltz.setTfo();
    const real_t Tfo = boltz.getTfo();

    // Retrieve the contributing processes
    std::vector<const Process2to2*> list_allowed(boltz.compute_contributing_processes());
    // if (n < list_allowed.size())
    // {
    // n = list_allowed.size();
    // size_t i = 1;
    // for (const Process2to2* proc : list_allowed)
    // {
    //   std::cout << i++ << ": " << proc->getMname() << '\n';
    // }
    // }

    const real_t relic_density = boltz.relic_density();

    std::cout << "The relic density computed with the full algorithm is " << relic_density << '\n';

    boltz.full_comput = false;
    const real_t relic_density_no_full = boltz.relic_density();

    std::cout << "The relic density computed with the non-full algorithm is " << relic_density_no_full << '\n';


    std::cout << "Oh2 (full)=" << relic_density << " Oh2 (no full)=" << relic_density_no_full
              << " xfo=" << boltz.getMassLBSM() / Tfo << "\n";


    std::cout << "Working on the processes grouped by final state threshold...\n";
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
    unsigned short int group_index(1u);
    for (const auto& group_of_processes : grouped_processes)
    {
      std::shared_ptr<std::vector<Process2to2>> group_ptr(std::make_shared<std::vector<Process2to2>>());
      std::cout << "Group " << group_index++ << ":\n";
      for (const Process2to2* proc : group_of_processes)
      {
        std::cout << "  - " << proc->getMname() << " (Threshold: " << proc->compute_final_state_treshold(input)
                  << ", m_V_min = " << proc->compute_final_state_treshold(input) * 2. / 3. << ")\n";
      }
    }
  }

  // Creating the process list with thresholds also for m_V



  return 0;
}