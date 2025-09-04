#include <fstream>
#include <iomanip>
#include <iostream>
#include <ranges>

#include "dp_u1f2to2/boltzmann.hpp"

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

  // The following 2 vectors are initialised with data from summary.dat
  // and before th jump 
  std::vector<real_t> m_V_jumps = {10.3912, 10.5926};    
  std::vector<real_t> g_f_jumps = {0.0544449, 0.0548977};

  // Using the data to compute the slope
  const real_t delta_g_over_delta_m = (g_f_jumps[1] - g_f_jumps[0])/(m_V_jumps[1] - m_V_jumps[0]);

  // This vector contains the trial m_V values, for which we will lineary interpolate g_f
  const std::vector<real_t> m_V_trial_vec(
      {10.7, 10.8, 10.9, 10.91, 10.92, 10.93, 10.94, 10.95, 10.96, 10.97, 10.98, 10.99, 11.0, 11.01});

  // Lambda function to linearly estrapolate g_f
  auto g_f_trial = [=](const real_t m_V_trial)
  { return g_f_jumps[1] + delta_g_over_delta_m * (m_V_trial - m_V_jumps[1]); };

  // Extending the first two vectors with the trial values
  for (auto& elem : m_V_trial_vec)
  {
    m_V_jumps.push_back(elem);
    g_f_jumps.push_back(g_f_trial(elem));
  }

  const constexpr bool add_points_after_jump = true;
  if (add_points_after_jump)
  {
    // m_V_jumps.push_back(11.007);
    // g_f_jumps.push_back(0.0517103);
    m_V_jumps.push_back(11.2202);
    g_f_jumps.push_back(0.0516542);
  }

  // Initialising the input Param_t structure
  {
    const real_t m_V = m_V_trial_vec [0];
    input.g_f = 0.1;
    input.m_V_3 = m_V;
    input.m_phi = 2. * m_V;
    input.m_chi_dm_1 = mchi_over_mv * m_V;
  }
  input.refresh();

  /* This variable is used to ensure that the list of contributing processes is printed only when 
     a new maximum number of processes is found for a given parameter set */
  size_t last_number_of_contributing_processes = 0;


  /*
    This loop iterates over pairs of (m_V, g_f) values, where:
      - m_V: vector of mediator masses
      - g_f: corresponding vector of coupling values (interpolated where needed)

    For each (m_V, g_f) pair:
      1. Updates the input parameters for the Boltzmann solver.
      2. Computes and writes the sigma_v table to a file "sigma_v_<m_V>.dat"
      3. Computes and writes the Weff table to a file "Weff_<m_V>.dat"
      4. Computes the relic density using both the full and non-full algorithms, and writes results to a file "Oh2_xfo_<m_V>.dat"
         Such a file contains 3 columns: 
          Oh2 (full), Oh2 (non-full), xfo
      5. Groups contributing processes by their final state threshold, and for each group:
         - Computes and writes group-specific sigma_v and Weff tables.
         - Prints process grouping information to stdout.
  */
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

    std::cout << "Generating and writing the sigma_v table for the totality of the processes\n";
    boltz.setTfo();
    const real_t Tfo = boltz.getTfo();

    /* Initial and final values of x = m_LBSM / T, for which the sigma_v table will be computed */
    constexpr const real_t x_in = 50., x_fin = 0.1; // x = m_LBSM / T
    
    constexpr const unsigned int n_points = 1000; // Number of points (minus 1) in the sigma_v table
    
    // Variable to create a logarithmic grid in x
    const real_t x_mult_factor = std::pow(x_fin / x_in, 1. / n_points); 
    
    std::array<real_t, n_points + 1> x_values; // Array to store the x values
    for (auto i : std::views::iota(0u, n_points))
    {
      x_values[i] = std::pow(x_mult_factor, i) * x_in;
    }

    sigma_v_file << "# x sigmav Y_eq sigmavnum sigmavden\n" << std::scientific << std::setprecision(4);
    for (const real_t& x : x_values)
    {
      const real_t T = 1. / x * input.getLightestBSMmass();
      const real_t sigma_v = boltz.getAverageSigmav(T);
      const real_t Y_eq = boltz.Yeq(T);
      const real_t sigmavnum = boltz.getAverageSigmav_coan_hightemp_num(T);
      const real_t sigmavden = boltz.getAverageSigmav_coan_hightemp_den(T);
      sigma_v_file << x << " " << sigma_v << " " << Y_eq << " " << sigmavnum << " " << sigmavden << "\n";
      if (std::abs((sigmavnum / sigmavden - sigma_v) / (sigmavnum / sigmavden + sigma_v)) > 1.0e-3)
      {
        std::cerr << "Warning: sigmavnum/sigmavden = " << sigmavnum / sigmavden << " and sigma_v = " << sigma_v
                  << " are not consistent\n";
      }
    }
    sigma_v_file.close();
    std::cout << "Sigma_v table written to " << filename.str() << "\n";

    filename.str("");
    filename.clear();
    filename << "out/scans/results/u1-theta23=0/Weff_" << m_V << ".dat";
    outfile.open(filename.str());
    if (!outfile.is_open())
    {
      std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";
      continue;
    }

    // Write the header for the file
    outfile << "# peff sqrts weff\n";

    // Retrieve the contributing processes
    std::vector<const Process2to2*> list_allowed(boltz.compute_contributing_processes());
    
    // Print the contributing processes if there are new ones
    if (last_number_of_contributing_processes < list_allowed.size())
    {
      last_number_of_contributing_processes = list_allowed.size();
      size_t i = 1;
      for (const Process2to2* proc_ptr : list_allowed)
      {
        std::cout << i++ << ": " << proc_ptr->getMname() << '\n';
      }
    }

    std::cout << "Generating and writing the peff, sqrts, and weff tables for the totality of the processes\n";
    for (size_t i = 0; i < boltz.getWefftabsize(); ++i)
    {
      real_t peff, sqrts, weff;
      boltz.get_g2_WeffTabElement(i, peff, sqrts, weff);
      outfile << peff << " " << sqrts << " " << weff << "\n";
    }

    outfile.close();
    std::cout << "Weff table written to " << filename.str() << "\n";
    std::cout << '\n';

    boltz.full_comput = true;
    const real_t relic_density = boltz.relic_density();

    std::cout << "The relic density computed with the full algorithm is " << relic_density << '\n';

    boltz.full_comput = false;
    const real_t relic_density_no_full = boltz.relic_density();

    std::cout << "The relic density computed with the non-full algorithm is " << relic_density_no_full << '\n';


    std::cout << "Oh2 (full)=" << relic_density << " Oh2 (no full)=" << relic_density_no_full
              << " xfo=" << boltz.getMassLBSM() / Tfo << "\n";

    // Create a file named Oh2_xfo_<value_of_mv>.dat
    filename.str("");
    filename.clear();
    filename << "out/scans/results/u1-theta23=0/Oh2_xfo_" << m_V << ".dat";
    outfile.open(filename.str());
    if (outfile.is_open())
    {
      outfile << std::scientific << std::setprecision(4) << relic_density << " " << relic_density_no_full << " "
              << boltz.getMassLBSM() / Tfo << "\n";
      outfile.close();
      std::cout << "Oh2_xfo written to " << filename.str() << "\n";
    }
    else
      std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";

    std::cout << "Working on the processes grouped by final state threshold...\n";
    std::vector<std::vector<const Process2to2*>> grouped_processes;
    if (!list_allowed.empty())
    {
      for (const Process2to2* proc_ptr : list_allowed)
      {
        bool added_to_group = false;
        const real_t threshold_proc = proc_ptr->compute_final_state_treshold(input);

        // Check if the process can be added to an existing group
        for (auto& group : grouped_processes)
        {
          const real_t threshold_group = group.front()->compute_final_state_treshold(input);

          if ((threshold_proc == 0. && threshold_group == 0.) or
              std::abs(threshold_proc - threshold_group) / std::abs(threshold_proc + threshold_group) < 1.0e-3)
          {
            group.push_back(proc_ptr);
            added_to_group = true;
            break;
          }
        }

        // If not added to any group, create a new group
        if (!added_to_group)
        {
          grouped_processes.emplace_back(std::vector<const Process2to2*>{proc_ptr});
        }
      }
    }

    // Output the grouped processes
    unsigned short int group_index(1u);
    for (const auto& group_of_processes : grouped_processes)
    {
      std::shared_ptr<std::vector<Process2to2>> group_ptr(std::make_shared<std::vector<Process2to2>>());
      std::cout << "Group " << group_index++ << ":\n";
      for (const Process2to2* proc_ptr : group_of_processes)
      {
        std::cout << "  - " << proc_ptr->getMname() << " (Threshold: " << proc_ptr->compute_final_state_treshold(input)
                  << ")\n";
        group_ptr->push_back(*proc_ptr);
      }

      BoltzmannSolver boltzgroup(input, group_ptr);
      boltzgroup.relic_density();

      // Compute the label for the group based on the threshold of the first process
      const real_t threshold_group = group_of_processes.front()->compute_final_state_treshold(input);
      std::ostringstream group_label_stream;
      group_label_stream << std::fixed << std::setprecision(3) << threshold_group;
      const std::string group_label = group_label_stream.str();

      // Create a file for the current group
      filename.str("");
      filename.clear();
      filename << "out/scans/results/u1-theta23=0/sigma_vg_" << m_V << "_group_" << group_label << ".dat";
      std::ofstream group_file(filename.str());
      if (!group_file.is_open())
      {
        std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";
        continue;
      }

      // Write the header for the file
      group_file << "# x sigmav\n" << std::scientific << std::setprecision(4);

      // Compute sigmav for each value in inverse_T_sigmav
      for (const real_t& x : x_values)
      {
        const real_t T = 1. / x * input.getLightestBSMmass();
        const real_t sigma_v = boltzgroup.getAverageSigmav(T);
        const real_t Y_eq = boltzgroup.Yeq(T);
        const real_t sigmavnum = boltzgroup.getAverageSigmav_coan_hightemp_num(T);
        const real_t sigmavden = boltzgroup.getAverageSigmav_coan_hightemp_den(T);
        group_file << x << " " << sigma_v << " " << Y_eq << " " << sigmavnum << " " << sigmavden << "\n";
        if (std::abs((sigmavnum / sigmavden - sigma_v) / (sigmavnum / sigmavden + sigma_v)) > 1.0e-3)
        {
          std::cerr << "Warning: sigmavnum/sigmavden = " << sigmavnum / sigmavden << " and sigma_v = " << sigma_v
                    << " are not consistent\n";
        }
      }

      group_file.close();
      std::cout << "Sigma_v table for group " << group_label << " written to " << filename.str() << "\n";

      // Create a file named Weff_<value_of_mv>.dat for the current group
      filename.str("");
      filename.clear();
      filename << "out/scans/results/u1-theta23=0/Weff_g_" << m_V << "_group_" << group_label << ".dat";
      outfile.open(filename.str());
      if (!outfile.is_open())
      {
        std::cerr << "Error: Could not open file " << filename.str() << " for writing.\n";
        continue;
      }

      // Write the header for the file
      outfile << "# peff sqrts weff\n";

      // Generate and write the peff, sqrts, and weff tables
      for (size_t i = 0; i < boltzgroup.getWefftabsize(); ++i)
      {
        real_t peff, sqrts, weff;
        boltzgroup.get_g2_WeffTabElement(i, peff, sqrts, weff);
        outfile << peff << " " << sqrts << " " << weff << "\n";
      }

      outfile.close();
      std::cout << "Weff table for group " << group_label << " written to " << filename.str() << "\n";
      
      std::cout << '\n';
    }
  }

  return 0;
}