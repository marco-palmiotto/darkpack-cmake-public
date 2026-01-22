#include "indirect.hpp"
// #define DEBUG

template <typename T> static inline T SQUARE(const T x) { return x * x; }

template <typename T> static inline T SQRT(const T x) { return pow(x, 0.5); }

namespace __SPEC_LIB_NAME__
{
  Indirectparam_t::Indirectparam_t(const Param_t& input_in) : input(input_in), run(input_in)
  {
    // This constructor assumes that the Param_t object input has been refreshed
    // for a proper assignment of the masses vector
    p_ptr = std::make_shared<std::vector<Process2to2>>();
    total_sigma_v = 0;
    find_annihilation_processes();
    sigma_v(*p_ptr);
  };

  Indirectparam_t::~Indirectparam_t() = default;

  void Indirectparam_t::find_annihilation_processes()
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
  };


  void Indirectparam_t::sigma_v(std::vector<Process2to2> Processes)
  {
    // This function computes the thermally averaged cross section for a 2 to 2 process
    // using a low velocity approximation (assuming processes with two identical particles in the initial state)
    for (size_t i = 0; i < Processes.size(); i++)
    {

      Process2to2 p = Processes[i];
      const real_t m1 = p.getMass(1, input);
      const real_t m3 = p.getMass(3, input);
      const real_t m4 = p.getMass(4, input);

      const real_t sqrt_s = 2 * m1; // s at threshold


      real_t s_v =
          (p.getSumSquaredAmpl(input, sqrt_s, 0.)) / (128. * M_PI * SQUARE(m1)) *
          SQRT(1. - (m3 * m3 + m4 * m4) / (2. * m1 * m1) + (pow(m3 * m3 - m4 * m4, 2.)) / (16. * pow(m1, 4.))) /
          (p.getSf34());
      sigma_v_process.push_back(s_v);
      total_sigma_v += s_v;
    }
  };

  std::vector<Process_1to2> Indirectparam_t::find_decay_processes(int particle)
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

  void Indirectparam_t::handle_1to2_decays(real_t sigma_v, int particle, real_t sqrt_s)
  {

    auto proc_list = find_decay_processes(particle);

    for (size_t i = 0; i < proc_list.size(); i++)
    {
      Process_1to2 proc = proc_list[i];
      auto branching_ratio = proc.get_branching_ratio(input);
      auto sigma_v_new = sigma_v * branching_ratio;

      // Arrays to store the fields, masses and energies of the final state particles
      std::array<int, 2> final_fields = {proc.get_field(3), proc.get_field(4)};
      std::array<real_t, 2> final_masses = {proc.get_mass(3, input), proc.get_mass(4, input)};
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
        }
        size_t energy_index = distance(energy_table.begin(), it);

        switch (final_fields[j])
        {
        case (corr::e):
          sigma_v_table[energy_index][0] += 0.5 * sigma_v;
          break;
        case (corr::mu):
          sigma_v_table[energy_index][1] += 0.5 * sigma_v;
          break;
        case (corr::tau):
          sigma_v_table[energy_index][2] += 0.5 * sigma_v;
          break;
        case (corr::u):
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
          break;
        case (corr::d):
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
          break;
        case (corr::s):
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
          break;
        case (corr::c):
          sigma_v_table[energy_index][4] += 0.5 * sigma_v;
          break;
        case (corr::b):
          sigma_v_table[energy_index][5] += 0.5 * sigma_v;
          break;
        case (corr::t):
          sigma_v_table[energy_index][6] += 0.5 * sigma_v;
          break;
        case (corr::W):
          sigma_v_table[energy_index][7] += 0.5 * sigma_v;
          break;
        case (corr::Z):
          sigma_v_table[energy_index][8] += 0.5 * sigma_v;
          break;
        case (corr::G):
          sigma_v_table[energy_index][9] += 0.5 * sigma_v;
          break;
        case (corr::A):
          sigma_v_table[energy_index][10] += 0.5 * sigma_v;
          break;
        case (corr::nu_e):
          sigma_v_table[energy_index][11] += 0.5 * sigma_v;
          break;
        case (corr::nu_mu):
          sigma_v_table[energy_index][12] += 0.5 * sigma_v;
          break;
        case (corr::nu_tau):
          sigma_v_table[energy_index][13] += 0.5 * sigma_v;
          break;
        default:
        {
          // Calling to handle decays
          handle_1to2_decays(sigma_v, final_fields[j], final_energies[j]);
        }
        }
      }
    }
  };

  void Indirectparam_t::fill_spectrum(std::vector<Process2to2> processes, std::vector<real_t> sigma_v_vector,
                                      real_t sqrt_s)
  {
    const size_t n_proc = sigma_v_vector.size();

    for (size_t i = 0; i < n_proc; i++)
    {
      Process2to2 proc = processes[i];
      real_t sigma_v = sigma_v_vector[i];
      if (sigma_v > 0.01 * total_sigma_v)
        dof++;

      // Arrays to store the fields, masses and energies of the final state particles
      std::array<int, 2> final_fields = {proc.getField(3), proc.getField(4)};
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
        }
        size_t energy_index = distance(energy_table.begin(), it);

        switch (final_fields[j])
        {
        case (corr::e):
          sigma_v_table[energy_index][0] += 0.5 * sigma_v;
          break;
        case (corr::mu):
          sigma_v_table[energy_index][1] += 0.5 * sigma_v;
          break;
        case (corr::tau):
          sigma_v_table[energy_index][2] += 0.5 * sigma_v;
          break;
        case (corr::u):
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
          break;
        case (corr::d):
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
          break;
        case (corr::s):
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
          break;
        case (corr::c):
          sigma_v_table[energy_index][4] += 0.5 * sigma_v;
          break;
        case (corr::b):
          sigma_v_table[energy_index][5] += 0.5 * sigma_v;
          break;
        case (corr::t):
          sigma_v_table[energy_index][6] += 0.5 * sigma_v;
          break;
        case (corr::W):
          sigma_v_table[energy_index][7] += 0.5 * sigma_v;
          break;
        case (corr::Z):
          sigma_v_table[energy_index][8] += 0.5 * sigma_v;
          break;
        case (corr::G):
          sigma_v_table[energy_index][9] += 0.5 * sigma_v;
          break;
        case (corr::A):
          sigma_v_table[energy_index][10] += 0.5 * sigma_v;
          break;
        case (corr::nu_e):
          sigma_v_table[energy_index][11] += 0.5 * sigma_v;
          break;
        case (corr::nu_mu):
          sigma_v_table[energy_index][12] += 0.5 * sigma_v;
          break;
        case (corr::nu_tau):
          sigma_v_table[energy_index][13] += 0.5 * sigma_v;
          break;
        default:
        {
          // Calling to handle decays
          handle_1to2_decays(sigma_v, final_fields[j], final_energies[j]);
        }
        }
      }
    };
  }

  void Indirectparam_t::read_file(std::string filename, int ncol, std::vector<std::vector<real_t>>& data)
  {
    data.clear();

    std::fstream file;
    file.open(filename);

    if (!file.is_open())
    {
      std::cout << "Could not open read_file " << filename << "\n";
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


  void Indirectparam_t::fetch_pppc4dmid_data()
  {
    int ind_matrix = 0;
    for (const auto& file : std::filesystem::directory_iterator(FILEPATH_PPPC4DMID))
    {
      std::string filename = file.path().string();

      read_file(filename, 15, PPPC4DMID.at(ind_matrix));
      DM_masses.emplace_back(std::stod(file.path().stem().string())); // Extract mass from filename
      ind_matrix++;
    }
  }


  void Indirectparam_t::fermi_energy_bins(const real_t& mass, const std::vector<std::vector<real_t>>& tab,
                                          std::vector<std::vector<real_t>>& final_spec)
  {
    size_t n_cols = PPPC4DMID.at(0).at(0).size();


    final_spec.resize(n_cols);


    for (int i = 0; i < 25; i++)
    {
      final_spec.at(i).at(0) = pow(10., log10(500) + (log10(500000) - log10(500)) / 24 * i) *
                               0.001; // 24 energy bins logarithmically spaced between 500 and 500000 MeV

      if (final_spec.at(i).at(0) > mass)
        for (size_t j = 1; j < n_cols; j++)
          final_spec.at(i).at(j) = 0;

      else
      {
        size_t index = 0;

        while (pow(10., tab[index][0] * mass) <= final_spec[i][0] && index < tab.size())
          index++;

        real_t K0 = pow(10., tab[index - 1][0] * mass);
        real_t K1 = pow(10., tab[index][0] * mass);
        real_t interpolation_factor = (final_spec[i][0] - K0) / (K1 - K0);

        for (size_t j = 1; j < n_cols; j++)
        {
          final_spec.at(i).at(j) =
              tab.at(index - 1).at(j) + (tab.at(index).at(j) - tab.at(index - 1).at(j)) * interpolation_factor;
        }
      }
    }
  }

  void Indirectparam_t::interpolate_spectrum_gamma(const real_t& mass,
                                                   std::vector<std::vector<real_t>>& interpolated_spectrum)
  {
    size_t n_cols = PPPC4DMID.at(0).at(0).size();
    size_t n_masses = DM_masses.size();

    interpolated_spectrum.resize(n_cols); // First column is log10(K/m)
    size_t index = 0;

    while (mass < DM_masses.at(index) && index < n_masses)
      index++;

    if (DM_masses.at(index) == mass)
      fermi_energy_bins(mass, PPPC4DMID[index], interpolated_spectrum);

    else if (DM_masses.at(index) > mass)
    {
      // Linear interpolation
      std::vector<std::vector<real_t>> spec_low, spec_high;
      fermi_energy_bins(DM_masses.at(index - 1), PPPC4DMID.at(index - 1), spec_low);
      fermi_energy_bins(DM_masses.at(index), PPPC4DMID.at(index), spec_high);
      real_t interpolation_factor = (mass - DM_masses.at(index - 1)) / (DM_masses.at(index) - DM_masses.at(index - 1));

      for (size_t i = 0; i < spec_low.size(); i++)
      {
        for (size_t j = 0; j < spec_low.at(0).size(); j++)
        {
          real_t y = spec_low.at(i).at(j) + (spec_high.at(i).at(j) - spec_low.at(i).at(j)) * interpolation_factor;
          interpolated_spectrum.at(i).push_back(y);
        }
      }
    }
  }

  void Indirectparam_t::spectrum_at_production(std::vector<std::vector<real_t>>& production_spectrum)
  {
    production_spectrum.resize(25); // 24 energy bins

    std::vector<std::vector<real_t>> tab; // To store the intermediate spectrum
    for (size_t iene = 0; iene < energy_table.size(); iene++)
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
          row.at(1) += tab.at(row_index).at(col_index) / row.at(0) * sigma_v_table.at(iene).at(col_index - 1) /
                       2.; // dN/dE value
        }
        row_index++;
      }
    }
  }

  real_t Indirectparam_t::integrate_spectrum(const std::vector<std::vector<real_t>>& spectrum, const real_t& e_min,
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

    return integral;
  }

  void Indirectparam_t::read_fermi_data()
  {
    std::string path = "/workspaces/darkpack-cmake/src/fermi_data/";
    int i = 0;

    // Organizing files in alphabetical order
    std::vector<std::string> file_list;
    for (const auto& file : std::filesystem::directory_iterator(path))
    {
      std::string filepath = file.path().string();
      file_list.push_back(filepath);
    }
    std::sort(file_list.begin(), file_list.end());

    for (const auto& filepath : file_list)
    {
      if (filepath == "/workspaces/darkpack-cmake/src/fermi_data/dSphs_list.dat")
      {
        logJ_factors.clear();
        read_file(filepath, 5, logJ_factors);
      }
      else if (filepath == "/workspaces/darkpack-cmake/src/fermi_data/limits_bb.txt" ||
               filepath == "/workspaces/darkpack-cmake/src/fermi_data/limits_tautau.txt")
      {
        continue; // Skip limits files
      }
      else
      {
#ifdef DEBUG
        // std::cout << "Reading Fermi data file: " << std::filesystem::absolute(filepath) << "\n";
#endif
        fermi_data.at(i).clear();
        read_file(filepath, 4, fermi_data.at(i));
        i++;
      }
    }
  }

  real_t Indirectparam_t::interpolate_likelihood(const int& dsph, const int& bin, const real_t& flux)
  {
    const auto data = fermi_data.at(dsph);
    auto total_size = static_cast<int>(data.size());
    int ind_flux = bin * 25; // Each energy bin has 25 flux values
    while (ind_flux < (bin + 1) * 25 - 1 && data.at(ind_flux).at(2) < flux && ind_flux < total_size)
      ind_flux++;
#ifdef DEBUG
    std::cout << ind_flux << " " << flux << "\n";
    std::cout << data.at(ind_flux).at(0) << " " << data.at(ind_flux).at(1) << " " << data.at(ind_flux).at(2) << " "
              << data.at(ind_flux).at(3) << "\n";
#endif
    real_t res = data.at(ind_flux - 1).at(3) + (data.at(ind_flux).at(3) - data.at(ind_flux - 1).at(3)) /
                                                   (data.at(ind_flux).at(2) - data.at(ind_flux - 1).at(2)) *
                                                   (flux - data.at(ind_flux - 1).at(2));

    return res;
  }

  void Indirectparam_t::compute_fluxes(std::vector<std::vector<real_t>>& production_spectrum)
  {
    spectrum_at_production(production_spectrum);

    int DM_candidate = input.getLightestBSMpart();
    real_t DM_mass = input.masses_vector.at(DM_candidate);

    for (int bin = 0; bin < 24; bin++)
    {
      real_t logE_min = fermi_data.at(0).at(bin * 25).at(0);
      real_t logE_max = fermi_data.at(0).at(bin * 25).at(1);

      real_t e_min = pow(10., logE_min);
      real_t e_max = pow(10., logE_max);

      fermi_eflux.at(bin) = integrate_spectrum(production_spectrum, e_min, e_max) / (4 * M_PI) / pow(DM_mass, 2.0);
    }
  }

  real_t Indirectparam_t::likelihood_one_dsph(const int& dsph, const real_t& logJ)
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

      if (J != 0.0)
      {
        real_t logJ_term = -0.5 * pow((logJ - logJ_obs) / sigma_j, 2.0);
        log_like += logJ_term;
      }
    }
    return log_like;
  }

  real_t Indirectparam_t::logx_interpol(const std::vector<std::vector<real_t>>& spect, real_t logx)
  /* Interpolates the value of a spectrum at an energy log10(E)=logx using the spectrum tabulated in energy spect */
  {
    int n1 = spect.size();
    if (logx < spect[0][0] || logx > spect[n1 - 1][0])
    {
      printf("logx out of range, %e %e %e\n", logx, spect[0][0], spect[n1 - 1][0]);
      exit(EXIT_FAILURE);
    }

    int index = 0;
    while (spect[index][0] <= logx && index < n1)
      index++;
    index--;

    real_t fact = (logx - spect[index][0]) / (spect[index + 1][0] - spect[index][0]);

    return (1. - fact) * spect[index][1] + fact * spect[index + 1][1];
  }


}; // namespace __SPEC_LIB_NAME__