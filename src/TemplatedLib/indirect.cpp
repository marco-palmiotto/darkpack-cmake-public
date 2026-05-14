#include "indirect.hpp"
#include "copy_to_c_struct.hpp"
// #define DEBUG
//  #define AVG_SV_CALCULATOR
template <typename T> static inline T SQUARE(const T x) { return x * x; }

template <typename T> static inline T SQRT(const T x) { return pow(x, 0.5); }

namespace __SPEC_LIB_NAME__
{
  Indirect_param_t::Indirect_param_t(const Param_t& input_in) : input(input_in), run(input_in)
  {
    // This constructor assumes that the Param_t object input has been refreshed
    // for a proper assignment of the masses vector

    run.RunCharmMass();
    p_ptr = std::make_shared<std::vector<Process2to2>>();
    total_sigma_v = 0;
    find_annihilation_processes();
    compute_sigma_v_at_threshold();
  };

  Indirect_param_t::~Indirect_param_t() = default;

  void Indirect_param_t::find_annihilation_processes()
  {
    p_ptr->clear();
    // Loop over all processes in corr::sumSquaredAmpl and find the ones with
    // two lightest BSM particles in the initial state

    for (auto& it : corr::squaredampl)
    {
      const auto current_process = it.first;
      std::array<Insertion, 4> particle_list;
      bool ap_flag = false;  // To be set true if the loop encounters ANTICHAR
      size_t curr_index = 0; // Index to track the filling of ParticleList

      for (size_t j = 0; j < current_process.size(); j++)
      {
        if (curr_index == 2 && (particle_list[0].field != input.getLightestBSMpart() ||
                                particle_list[1].field != input.getLightestBSMpart()))
        {
          continue; // This process does not have two lightest BSM particles in the initial state
                    // (checked after the first two particles have been filled)
        }

        if (!ap_flag && current_process[j] == ANTICHAR)
        {
          ap_flag = true;
          continue; // Move to the next particle
        }

        if (ap_flag && current_process[j] != ANTICHAR)
        {
          const Insertion Particle(current_process[j] - EMPTYCHAR,
                                   false); // This is the antiparticle of the previous one
          particle_list[curr_index] = Particle;

          ap_flag = false; // Reset the flag
          curr_index++;
          if (curr_index > 3)
          {
            p_ptr->emplace_back(particle_list);
            break; // We have filled all four particles, exit the loop
          }
          continue; // Move to the next particle
        }

        if (!ap_flag && current_process[j] != ANTICHAR)
        {
          const Insertion Particle(current_process[j] - EMPTYCHAR, true);
          particle_list[curr_index] = Particle;
          curr_index++;
          if (curr_index > 3)
          {
            p_ptr->emplace_back(particle_list);
            break; // We have filled all four particles, exit the loop
          }
          continue;
        }
      }
    }
    for (auto& process : *p_ptr)
    {
      process.setRunningData(&run);
      process.setRunningExternal();
    }
  };


  void Indirect_param_t::compute_sigma_v_at_threshold()
  {
    // This function computes the thermally averaged cross section for a 2 to 2 process
    // using a low velocity approximation (assuming processes with two identical particles in the initial state)
    real_t sqrt_s = 2 * input.masses_vector.at(input.getLightestBSMpart()); // s at threshold
    Param_t input_copy = input; // Create a copy of the input to avoid modifying the original one
    run.HandleParamRunning(input_copy, sqrt_s);
#ifdef DEBUG
    std::cerr << "After first running: " << input_copy.m_b << "\n";
#endif
    for (Process2to2& p : *p_ptr)
    {
#ifdef DEBUG
      std::cerr << "Process found: " << p.getName() << "\n";
#endif
      const real_t m1 = p.getMass(1, input_copy);
      const real_t threshold = 2. * m1; // Threshold energy for the process

      if (std::abs((sqrt_s - threshold) / (sqrt_s + threshold) > 1e-3))
      {
        sqrt_s = threshold;
        input_copy = input; // Reset the input copy to the original input

        run.HandleParamRunning(input_copy, sqrt_s); // Run the parameters to the threshold energy
#ifdef DEBUG
        std::cerr << "After next running: " << input_copy.m_b << "\n";
#endif
      }
      const real_t m3 = p.getMass(3, input_copy);
      const real_t m4 = p.getMass(4, input_copy);

      real_t p1, p3, sij[5][5];

      // Returning 0 if the process is kinematically forbidden
      if (!p.setKinematics(input_copy, sqrt_s, 0., p1, p3, sij))
      {
#ifdef DEBUG
        std::cerr << "Process " << p.getName() << " is kinematically forbidden at energy " << sqrt_s
                  << ". Skipping...\n";
#endif
        sigma_v_process.push_back(0.);
        continue;
      }

      /*
      update_kinematics(input_copy, sij);

      cparam_s c_input = copy_to_c_struct(input_copy);

      // Defining a lambda function for the amplitude
      // Everything has to be passed by value because we use multithreading
      auto lambda_opti = [=, this](real_t cosine) mutable
      {
        real_t p1_local, p3_local, sij_local[5][5];
        if (!p.setKinematics(input_copy, sqrt_s, cosine, p1_local, p3_local, sij_local))
          return 0.;
        update_kinematics(input_copy, sij_local);
        real_t m1 = p.getMass(1, input_copy);
        real_t m2 = p.getMass(2, input_copy);
        real_t m3 = p.getMass(3, input_copy);
        real_t m4 = p.getMass(4, input_copy);

        real_t s = sqrt_s * sqrt_s;
        real_t t = m1 * m1 - 2.0 * sij_local[1][3] + m3 * m3;
        real_t u = m1 * m1 - 2.0 * sij_local[1][4] + m4 * m4;

        c_input = update_c_kinematics(input_copy, c_input);
        const real_t SumSquaredAmpl = p.getSumSquaredAmpl_ptr()(&c_input).real;
        return (std::isnormal(SumSquaredAmpl) && SumSquaredAmpl > 0.) ? SumSquaredAmpl : 0.0;
      };

      real_t integral, discrepancy;
      const real_t maxdiscrepancy = 5.0e-3;

      integral = advmath::integrate_gauss_comparative(-1., 1., lambda_opti, maxdiscrepancy, &discrepancy);

      if (discrepancy > maxdiscrepancy)
        integral = advmath::integrate_trap(-1., 1., lambda_opti, maxdiscrepancy, &discrepancy);

      real_t temp_squared_amp = integral / 2.;
      */
      real_t temp_squared_amp = p.getSumSquaredAmpl(input_copy, sqrt_s, 0.);
#ifdef DEBUG
      std::cerr << "Process: " << p.getName() << " at energy " << sqrt_s << " with squared amplitude "
                << temp_squared_amp << "\n";
#endif
      real_t s_v = (std::isnormal(temp_squared_amp)) ? temp_squared_amp / (128. * M_PI * SQUARE(m1)) *
                                                           SQRT(1. - (m3 * m3 + m4 * m4) / (2. * m1 * m1) +
                                                                (pow(m3 * m3 - m4 * m4, 2.)) / (16. * pow(m1, 4.))) *
                                                           3.8937966e+8 * 2.997900e-26 / (p.getSf34())
                                                     : 0.;
      printf("Process %s has a sigma v of %e\n", p.getName().c_str(), s_v);
      sigma_v_process.push_back(s_v);
      total_sigma_v += s_v;
    }
    printf("Total sigma_v = %e\n", total_sigma_v);
    for (int i = 0; i < sigma_v_process.size(); i++)
    {
      if (sigma_v_process[i] > 0.01 * total_sigma_v)
        dof++;
    };
  }

  std::vector<Process_1to2> Indirect_param_t::find_decay_processes(int particle)
  { // Finds all the decay processes for a given particle in the 1 to 2 hash table

    std::vector<Process_1to2> proc_list;

    for (auto& it : corr::squaredampl_1to2)
    {
      const auto current_process = it.first;
      std::array<Insertion, 3> particle_list;
      bool ap_flag = false;  // To be set true if the loop encounters ANTICHAR
      size_t curr_index = 0; // Index to track the filling of ParticleList

      for (size_t j = 0; j < current_process.size(); j++)
      {
        if (curr_index == 1 && particle_list[0].field != particle)
        {
          continue; // This process does not have the desired initial state
        }

        if (!ap_flag && current_process[j] == ANTICHAR)
        {
          ap_flag = true;
          continue; // Move to the next particle
        }

        if (ap_flag && current_process[j] != ANTICHAR)
        {
          const Insertion Particle(current_process[j] - EMPTYCHAR,
                                   false); // This is the antiparticle of the previous one
          particle_list[curr_index] = Particle;

          ap_flag = false; // Reset the flag
          curr_index++;
          if (curr_index > 2)
          {
            proc_list.emplace_back(particle_list);
            break; // We have filled all three particles, exit the loop
          }
          continue; // Move to the next particle
        }

        if (!ap_flag && current_process[j] != ANTICHAR)
        {
          const Insertion Particle(current_process[j] - EMPTYCHAR, true);
          particle_list[curr_index] = Particle;
          curr_index++;
          if (curr_index > 2)
          {
            proc_list.emplace_back(particle_list);
            break; // We have filled all three particles, exit the loop
          }
          continue;
        }
      }
    }
    return proc_list;
  }

  void Indirect_param_t::handle_1to2_decays(real_t sigma_v, int particle, real_t sqrt_s)
  {
    std::vector<real_t> part_widths;
    real_t total_width = 0;
    printf("sqrtS = %e\n", sqrt_s);
    auto proc_list = find_decay_processes(particle);

    Param_t input_state = input;
    real_t initial_mass = input_state.masses_vector.at(particle);
    run.HandleParamRunning(input_state, initial_mass); // Run the parameters to the threshold energy

    for (size_t i = 0; i < proc_list.size(); i++)
    {

      Process_1to2 proc = proc_list[i];
      proc.set_running_external();
      proc.set_running_data(&run);

      real_t branching_ratio = proc.compute_branching_ratio(input_state);

      std::cout << proc.get_name() << " has a BR of " << branching_ratio << std::endl;
      if (branching_ratio == 0.)
      {
        std::cerr << "Process " << proc.get_name() << " has a branching ratio of 0. Skipping...\n";
        continue;
      }

      auto sigma_v_new = sigma_v * branching_ratio;

      // Arrays to store the fields, masses and energies of the final state particles
      std::array<int, 2> final_fields = {proc.get_field(2), proc.get_field(3)};
      std::array<real_t, 2> final_masses = {proc.get_mass(2, input_state), proc.get_mass(3, input_state)};
      printf("Masses of the particles are %f and %f\n", final_masses[0], final_masses[1]);
      std::array<real_t, 2> final_energies = {
          0.5 * (sqrt_s + (SQUARE(final_masses[0]) - SQUARE(final_masses[1])) / sqrt_s),
          sqrt_s - 0.5 * (sqrt_s + (SQUARE(final_masses[0]) - SQUARE(final_masses[1])) / sqrt_s)};

      // Loop over the final state particles, check to see if their energies are already in the table
      for (size_t j = 0; j < 2; j++)
      {
        auto it = std::find(energy_table.begin(), energy_table.end(), final_energies[j]);
        if (it == energy_table.end())
        {
          // Energy not found, add it to the table and update the index
          energy_table.push_back(final_energies[j]);
          it = std::find(energy_table.begin(), energy_table.end(), final_energies[j]);
        }
        size_t energy_index = distance(energy_table.begin(), it);

        if (energy_index == energy_table.size() - 1)
        {
          for (auto& channel : sigma_v_table)
          {
            channel.resize(energy_index + 1);
          }
        }

        switch (final_fields[j])
        {
        case (corr::e):
          sigma_v_table[0][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::mu):
          sigma_v_table[1][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::tau):
          sigma_v_table[2][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::u):
          sigma_v_table[3][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::d):
          sigma_v_table[3][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::s):
          sigma_v_table[3][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::c):
          sigma_v_table[4][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::b):
          sigma_v_table[5][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::t):
          sigma_v_table[6][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::W):
          sigma_v_table[7][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::Z):
          sigma_v_table[8][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::G):
          sigma_v_table[9][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::A):
          sigma_v_table[10][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::nu_e):
          sigma_v_table[11][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::nu_mu):
          sigma_v_table[12][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        case (corr::nu_tau):
          sigma_v_table[13][energy_index] += 0.5 * sigma_v * branching_ratio;
          break;
        default:
        {
          // Calling to handle decays
          handle_1to2_decays(0.5 * sigma_v, final_fields[j], final_energies[j]);
        }
        }
      }
    }
  };

  void Indirect_param_t::fill_spectrum()
  {
    const size_t n_proc = sigma_v_process.size();

    int DM_candidate = input.getLightestBSMpart();
    real_t sqrt_s = 2 * input.masses_vector.at(DM_candidate);
    for (size_t i = 0; i < n_proc; i++)
    {
      if (sigma_v_process[i] < 0.01 * total_sigma_v)
        continue; // Skip processes that do not contribute significantly to the total cross section
      Process2to2 proc = (*p_ptr)[i];
      real_t sigma_v = sigma_v_process[i];

      // Arrays to store the fields, masses and energies of the final state particles
      std::array<int, 2> final_fields = {proc.getField(2), proc.getField(3)};
      std::array<real_t, 2> final_masses = {proc.getMass(3, input), proc.getMass(4, input)};
      std::array<real_t, 2> final_energies = {
          0.5 * (sqrt_s + (SQUARE(final_masses[0]) - SQUARE(final_masses[1])) / sqrt_s),
          sqrt_s - 0.5 * (sqrt_s + (SQUARE(final_masses[0]) - SQUARE(final_masses[1])) / sqrt_s)};

      // Loop over the final state particles, check to see if their energies are already in the table
      for (size_t j = 0; j < 2; j++)
      {
        auto it = std::find(energy_table.begin(), energy_table.end(), final_energies[j]);
        if (it == energy_table.end())
        {
          // Energy not found, add it to the table and update the index
          energy_table.push_back(final_energies[j]);
          it = std::find(energy_table.begin(), energy_table.end(), final_energies[j]);
        }
        size_t energy_index = std::distance(energy_table.begin(), it);
        if (energy_index == energy_table.size() - 1)
        {
          for (auto& channel : sigma_v_table)
          {
            channel.resize(energy_index + 1);
          }
        }

        switch (final_fields[j])
        {
        case (corr::e):
          sigma_v_table[0][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::mu):
          sigma_v_table[1][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::tau):
          sigma_v_table[2][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::u):
          sigma_v_table[3][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::d):
          sigma_v_table[3][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::s):
          sigma_v_table[3][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::c):
          sigma_v_table[4][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::b):
          sigma_v_table[5][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::t):
          sigma_v_table[6][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::W):
          sigma_v_table[7][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::Z):
          sigma_v_table[8][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::G):
          sigma_v_table[9][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::A):
          sigma_v_table[10][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::nu_e):
          sigma_v_table[11][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::nu_mu):
          sigma_v_table[12][energy_index] += 0.5 * sigma_v;
          break;
        case (corr::nu_tau):
          sigma_v_table[13][energy_index] += 0.5 * sigma_v;
          break;
        default:
        {
          // Calling to handle decays
          handle_1to2_decays(0.5 * sigma_v, final_fields[j], final_energies[j]);
        }
        }
      }
    };
  }

  /*
  real_t Indirect_param_t::virtual_particle_integration(int particle, )
  {
    auto proc_list = find_decay_processes(particle);

    for (auto& proc : proc_list)
    {
      real_t sum_squared_amp =
    }

    sVarW sVarStr;
    sVarW* sVar = &sVarStr;

    n_Gauss = N;
    m0_ = m0;
    if (w1 == 0 && w2 == 0)
      return decayPcm(m0, m1, m2);
    else if (w1 == 0)
    {
      if (m1 > m0)
        return 0;
      m1_ = m1;
      m2_p = m2;
      w2_ = w2;
      y_21 = atan(-m2 / w2);
      y_22 = atan(((m0 - m1) * (m0 - m1) - m2 * m2) / (m2 * w2));
      if (n_Gauss)
        return gauss_arg(intDecay2_, sVar, 0, 1, n_Gauss) / M_PI / ME(m0, m1, m2);
      else
        return simpson_arg(intDecay2_, sVar, 0, 1, 1.E-3, NULL) / M_PI / ME(m0, m1, m2);
    }
    else if (w2 == 0)
    {
      if (m2 > m0)
        return 0;
      m1_ = m2;
      m2_p = m1;
      w2_ = w1;
      y_21 = atan(-m1 / w1);
      y_22 = atan(((m0 - m2) * (m0 - m2) - m1 * m1) / (m1 * w1));
      if (n_Gauss)
        return gauss_arg(intDecay2_, sVar, 0, 1, n_Gauss) / M_PI / ME(m0, m1, m2);
      else
        return simpson_arg(intDecay2_, sVar, 0, 1, 1.E-3, NULL) / M_PI / ME(m0, m1, m2);
    }
    else
    {
      if (w1 <= w2)
      {
        w1_ = w1;
        w2_ = w2;
        m1_p = m1;
        m2_p = m2;
      }
      else
      {
        w1_ = w2;
        w2_ = w1;
        m1_p = m2;
        m2_p = m1;
      }
      y_11 = atan(-m1 / w1);
      y_12 = atan((m0 * m0 - m1 * m1) / (m1 * w1));

      //((m0-m2)^2-m1_p*m1_p)/(w1_*m1_p) =tan(y1);

      // m1_=m1_p*m1_p + w1_*m1_p*tan(y1);


      // double y_1m=atan( ((m0-m2+3*w2 )*(m0-m2+3*w2)-m1*m1)/(m1*w1));

      // printf("y_11 y_1m y_12 = %E %E %E\n", y_11, y_1m, y_12);
      // displayPlot("intDecay1_1","y",y_11,y_1m,0,1,"",0,intDecay1,sVar);
      // displayPlot("intDecay1_2","y",y_1m,y_12,0,1,"",0,intDecay1,sVar);

      // y_11=atan( ((m0-m2+3*w2 )*(m0-m2)-m1*m1)/(m1*w1));

      // y_11+pow(x,nn)*(y_12-y_11)=y_1m

      // double x0=pow( (y_1m-y_11)/(y_12-y_11), 1./nn);

      // printf("x0=%e\n",x0);
      // displayPlot("intDecay1_","x",0,1,0,1,"",0,intDecay1_,sVar);
      // double sum=0;
      // for(int i=0;i<10;i++) sum+=intDecay1_((i+0.5)/10.,sVar)/10;
      // printf("=====================\n");
      double res = simpson_arg(intDecay1_, sVar, 0, 1, 1E-5, NULL);
      // printf("sum=%E res=%E\n",sum,res);
      return res / (M_PI * M_PI) / ME(m0, m1, m2);
    }
  }
  */

  void Indirect_param_t::read_file(std::string filename, int ncol, std::vector<std::vector<real_t>>& data)
  {
    data.clear();

    std::fstream file;
    file.open(filename);

    if (!file.is_open())
    {
      std::cerr << "Could not open read_file " << filename << "\n";
      exit(1);
    }

    std::string line;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#')
        continue; // Skip comment lines
      real_t val;
      int count = 0;
      std::vector<real_t> values;
      values.clear();
      std::stringstream ss(line);
      if (!std::isdigit(line[0]) && line[0] != '-' && line[0] != '+')
      {
        std::string dummy;
        ss >> dummy; // Skip non-numeric starting values
      }
      while (ss >> val)
      {
        values.push_back(val);
        count++;
      }

      if (count == ncol)
      {
        data.emplace_back(values);
        continue;
      }
      else
      {
        int ind_start = count - ncol;
        data.emplace_back(values.begin() + ind_start, values.end());
      }
    }
  }


  void Indirect_param_t::fetch_pppc4dmid_data()
  {
    int ind_matrix = 0;
    for (const auto& file : std::filesystem::directory_iterator(FILEPATH_PPPC4DMID))
    {
      DM_masses.emplace_back(std::stod(file.path().stem().string())); // Extract mass from filename
    }
    std::sort(DM_masses.begin(), DM_masses.end());
    for (auto& mass : DM_masses)
    {
      std::string filename = std::format("{}/{:.0f}GeV.txt", FILEPATH_PPPC4DMID, mass);
      read_file(filename, 15, PPPC4DMID.at(ind_matrix));

      ind_matrix++;
    }
  }


  void Indirect_param_t::fermi_energy_bins(const real_t& mass, const std::vector<std::vector<real_t>>& tab,
                                           std::vector<std::vector<real_t>>& final_spec)
  {
    size_t n_cols = tab.at(0).size();

    final_spec.resize(25);


    for (int i = 0; i < 25; i++)
    {
      final_spec.at(i).resize(n_cols);
      final_spec.at(i).at(0) = pow(10., log10(500) + (log10(500000) - log10(500)) / 24 * i) *
                               0.001; // 24 energy bins logarithmically spaced between 500 and 500000 MeV

      if (final_spec.at(i).at(0) > mass)
        for (size_t j = 1; j < n_cols; j++)
          final_spec.at(i).at(j) = 0;

      else
      {
        size_t index = 0;

        while (pow(10., tab[index][0]) * mass <= final_spec[i][0] && index < tab.size())
          index++;

        index--;
        real_t K0 = pow(10., tab[index][0]) * mass;
        real_t K1 = pow(10., tab[index + 1][0]) * mass;
        real_t interpolation_factor = (final_spec[i][0] - K0) / (K1 - K0);

        for (size_t j = 1; j < n_cols; j++)
        {
          final_spec.at(i).at(j) =
              tab.at(index).at(j) + (tab.at(index + 1).at(j) - tab.at(index).at(j)) * interpolation_factor;
        }
      }
    }
  }

  void Indirect_param_t::interpolate_spectrum_gamma(const real_t& mass,
                                                    std::vector<std::vector<real_t>>& interpolated_spectrum)
  {
    size_t n_masses = DM_masses.size();
    interpolated_spectrum.clear();

    interpolated_spectrum.resize(25); // First column is log10(K/m)
    size_t index = 0;

    if (mass < DM_masses.at(0) || mass > DM_masses.back())
      std::cerr << "DM mass out of bounds for PPPC4DMID interpolation" << std::endl;
    while (mass > DM_masses.at(index) && index < n_masses)
      index++;

    index--;

    // prevents the singularities due to the opening of W, Z, or t channels
    if (mass < 90 && mass > 80)
      index++;
    if (mass < 100 && mass > 90)
      index++;
    if (mass < 200 && mass > 173)
      index++;


    if (DM_masses.at(index) == mass)
      fermi_energy_bins(mass, PPPC4DMID[index], interpolated_spectrum);
    else
    {
      // Linear interpolation
      std::vector<std::vector<real_t>> spec_low, spec_high;
      fermi_energy_bins(DM_masses.at(index), PPPC4DMID.at(index), spec_low);
      fermi_energy_bins(DM_masses.at(index + 1), PPPC4DMID.at(index + 1), spec_high);
      real_t interpolation_factor = (mass - DM_masses.at(index)) / (DM_masses.at(index + 1) - DM_masses.at(index));

      for (size_t i = 0; i < spec_low.size(); i++)
      {
        for (size_t j = 0; j < spec_low.at(0).size(); j++)
        {
          if (j == 0)
          {
            interpolated_spectrum.at(i).push_back(spec_low.at(i).at(j)); // Energy bin is the same for both spectra
            continue;
          }
          real_t y = spec_low.at(i).at(j) + (spec_high.at(i).at(j) - spec_low.at(i).at(j)) * interpolation_factor;
          interpolated_spectrum.at(i).push_back(y);
        }
      }
    }
  }

  void Indirect_param_t::spectrum_at_production(std::vector<std::vector<real_t>>& production_spectrum)
  {
    production_spectrum.resize(25); // 24 energy bins

    size_t N = energy_table.size();
    real_t mass_DM = input.getLightestBSMmass();
    if (N == 0)
    {
      fill_spectrum();
      N = energy_table.size();
    }
    std::vector<std::vector<real_t>> tab; // To store the intermediate spectrum
    for (size_t iene = 0; iene < N; iene++)
    {
      real_t energy = energy_table.at(iene);
      interpolate_spectrum_gamma(energy, tab);

      size_t row_index = 0;
      for (auto& row : production_spectrum)
      {
        row.resize(2);                       // First column is energy, second is dN/dE
        row.at(0) = tab.at(row_index).at(0); // Energy of the bin
        for (size_t col_index = 1; col_index < tab.at(0).size(); col_index++)
        {
          row.at(1) +=
              tab.at(row_index).at(col_index) / row.at(0) * sigma_v_table.at(col_index - 1).at(iene); // dN/dE value
        }
        row_index++;
      }
    }

    /*
    for (int i = 0; i < production_spectrum.size(); i++)
    {
      printf("bin %d, prod_spec %e\n", i, production_spectrum.at(i).at(1));
    }
    getchar();
    */
  }

  real_t Indirect_param_t::integrate_spectrum(const std::vector<std::vector<real_t>>& spectrum, const real_t& e_min,
                                              const real_t& e_max)
  {
    real_t integral = 0.0;
    size_t n_bins = spectrum.size();

    for (size_t i = 0; i < n_bins - 1; i++)
    {
      real_t e_low = spectrum[i][0];
      real_t e_high = spectrum[i + 1][0];

      if (e_high < e_min || e_low > e_max)
        continue; // Skip bins outside the integration range

      // Determine the effective bin edges within the integration range
      real_t bin_start = std::max(e_low, e_min);
      real_t bin_end = std::min(e_high, e_max);

      // Linear interpolation of dN/dE at the bin edges
      real_t dN_dE_start;
      real_t dN_dE_end;
      if (bin_start == e_low)
        dN_dE_start = spectrum[i][1] * spectrum[i][0];
      else
        dN_dE_start = spectrum[i][1] * spectrum[i][0] +
                      bin_start * (spectrum[i + 1][1] - spectrum[i][1]) * (bin_start - e_low) / (e_high - e_low);
      if (bin_end == e_high)
        dN_dE_end = spectrum[i + 1][1] * spectrum[i + 1][0];
      else
        dN_dE_end = spectrum[i][1] * spectrum[i][0] +
                    bin_end * (spectrum[i + 1][1] - spectrum[i][1]) * (bin_end - e_low) / (e_high - e_low);

      // Trapezoidal rule for integration over the bin
      integral += 0.5 * (dN_dE_start + dN_dE_end) * (bin_end - bin_start);
    }

    return integral * 1000;
  }

  void Indirect_param_t::read_fermi_data()
  {
    int i = 0;

    // Organizing files in alphabetical order
    std::vector<std::string> file_list;
    for (const auto& file : std::filesystem::directory_iterator(FILEPATH_FERMI))
    {
      std::string filepath = file.path().string();
      file_list.push_back(filepath);
    }
    std::sort(file_list.begin(), file_list.end());

    for (const auto& filepath : file_list)
    {
      std::filesystem::path p(filepath);

      if (p.filename() == "dSphs_list.dat") // This is the file containing the logJ data for the dSph's
      {
        logJ_factors.clear();
        read_file(filepath, 5, logJ_factors);
      }
      else if (p.filename() == "limits_bb.txt" || p.filename() == "limits_tautau.txt")
      {
        continue; // Skip limits files
      }
      else
      {
#ifdef DEBUG
        std::cerr << "Reading Fermi data file: " << std::filesystem::absolute(filepath) << "\n";
#endif
        fermi_data.at(i).clear();
        read_file(filepath, 4, fermi_data.at(i));
        i++;
      }
    }
  }

  real_t Indirect_param_t::interpolate_likelihood(const int& dsph, const int& bin, const real_t& flux)
  {
    const auto data = fermi_data.at(dsph);
    auto total_size = static_cast<int>(data.size());
    int ind_flux = bin * 25; // Each energy bin has 25 flux values
    while (ind_flux < (bin + 1) * 25 - 1 && data.at(ind_flux).at(2) <= flux && ind_flux < total_size)
      ind_flux++;
#ifdef DEBUG
    std::cerr << ind_flux << " " << flux << "\n";
    std::cerr << data.at(ind_flux).at(0) << " " << data.at(ind_flux).at(1) << " " << data.at(ind_flux).at(2) << " "
              << data.at(ind_flux).at(3) << "\n";
#endif
    real_t res = data.at(ind_flux - 1).at(3) + (data.at(ind_flux).at(3) - data.at(ind_flux - 1).at(3)) /
                                                   (data.at(ind_flux).at(2) - data.at(ind_flux - 1).at(2)) *
                                                   (flux - data.at(ind_flux - 1).at(2));

    return res;
  }

  void Indirect_param_t::compute_fluxes(std::vector<std::vector<real_t>>& production_spectrum)
  {
    spectrum_at_production(production_spectrum);


    int DM_candidate = input.getLightestBSMpart();
    real_t DM_mass = input.masses_vector.at(DM_candidate);

    for (int bin = 0; bin < 24; bin++)
    {
      real_t e_min = fermi_data.at(0).at(bin * 25).at(0) * 0.001;
      real_t e_max = fermi_data.at(0).at(bin * 25).at(1) * 0.001; // Converting from MeV to GeV

      fermi_eflux.at(bin) =
          integrate_spectrum(production_spectrum, e_min, e_max) / (4 * M_PI) / (2. * pow(DM_mass, 2.0));
    }
  }

  real_t Indirect_param_t::likelihood_one_dsph(const int& dsph, const real_t& logJ)
  {
    real_t logJ_obs = logJ_factors.at(dsph).at(0);
    real_t sigma_j = logJ_factors.at(dsph).at(1);
    real_t J_obs = pow(10., logJ_obs);
    real_t J = (logJ != 0.0) ? pow(10., logJ) : 0.0;


    real_t log_like = 0.0;
    for (size_t bin = 0; bin < 24; bin++)
    {
      real_t flux_J = J * fermi_eflux.at(bin);
      log_like += interpolate_likelihood(dsph, bin, flux_J);
    }

    if (J != 0.0)
    {
      real_t logJ_term = -0.5 * pow((logJ - logJ_obs) / sigma_j, 2.0);
      log_like += logJ_term;
    }
    return log_like;
  }

  real_t Indirect_param_t::logx_interpol(const std::vector<std::vector<real_t>>& spect, real_t logx)
  /* Interpolates the value of a spectrum at an energy log10(E)=logx using the spectrum tabulated in energy spect */
  {
    int n1 = spect.size();
    if (logx < spect[0][0] || logx > spect[n1 - 1][0])
    {
      printf("logx out of range, %e %e %e\n", logx, spect[0][0], spect[n1 - 1][0]);
      exit(EXIT_FAILURE);
    }
    int index = 0;
    while (spect.at(index).at(0) <= logx)
    {
      index++;
      if (index == n1)
        break;
    }
    index--;

    if (index == n1 - 1)
    {
      return spect.at(index).at(1);
    }

    real_t fact = (logx - spect.at(index).at(0)) / (spect.at(index + 1).at(0) - spect.at(index).at(0));

    return (1. - fact) * spect.at(index).at(1) + fact * spect.at(index + 1).at(1);
  }


}; // namespace __SPEC_LIB_NAME__