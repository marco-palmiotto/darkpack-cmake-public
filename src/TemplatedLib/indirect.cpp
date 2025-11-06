#include "indirect.hpp"


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

  // Indirectparam_t::~Indirectparam_t() = default;

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

        if (final_fields[j] == corr::e)
          sigma_v_table[energy_index][0] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::mu)
          sigma_v_table[energy_index][1] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::tau)
          sigma_v_table[energy_index][2] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::u)
          sigma_v_table[energy_index][3] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::d)
          sigma_v_table[energy_index][3] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::s)
          sigma_v_table[energy_index][3] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::c)
          sigma_v_table[energy_index][4] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::b)
          sigma_v_table[energy_index][5] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::t)
          sigma_v_table[energy_index][6] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::W)
          sigma_v_table[energy_index][7] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::Z)
          sigma_v_table[energy_index][8] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::G)
          sigma_v_table[energy_index][9] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::A)
          sigma_v_table[energy_index][10] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::nu_e)
          sigma_v_table[energy_index][11] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::nu_mu)
          sigma_v_table[energy_index][12] += 0.5 * sigma_v_new;
        else if (final_fields[j] == corr::nu_tau)
          sigma_v_table[energy_index][13] += 0.5 * sigma_v_new;
        else
        {
          // Called recursively to handle subsequent decays
          handle_1to2_decays(sigma_v_new, final_fields[j], final_energies[j]);
        }
      }
    }
  }

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

        if (final_fields[j] == corr::e)
          sigma_v_table[energy_index][0] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::mu)
          sigma_v_table[energy_index][1] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::tau)
          sigma_v_table[energy_index][2] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::u)
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::d)
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::s)
          sigma_v_table[energy_index][3] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::c)
          sigma_v_table[energy_index][4] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::b)
          sigma_v_table[energy_index][5] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::t)
          sigma_v_table[energy_index][6] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::W)
          sigma_v_table[energy_index][7] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::Z)
          sigma_v_table[energy_index][8] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::G)
          sigma_v_table[energy_index][9] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::A)
          sigma_v_table[energy_index][10] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::nu_e)
          sigma_v_table[energy_index][11] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::nu_mu)
          sigma_v_table[energy_index][12] += 0.5 * sigma_v;
        else if (final_fields[j] == corr::nu_tau)
          sigma_v_table[energy_index][13] += 0.5 * sigma_v;
        else
        {
          // Calling to handle decays
          handle_1to2_decays(sigma_v, final_fields[j], final_energies[j]);
        }
      }
    }
  };

  std::vector<std::vector<double>> Indirectparam_t::read_file(std::string filename, int ncol)
  {
    std::vector<std::vector<double>> data;
    std::fstream file;
    file.open(filename);

    if(!file.is_open())
    {
      std::cout << "Could not open read_file " << filename << "\n";
      exit(1);
    }

    std::vector<std::vector<double>> data(ncol);
    while(std::getline(file, line))
    {
      double val;
      int count=0;
      std::vector<double> values(ncol);
      std::stringstream ss(line);
      while(ss >> val)
      {
        values.push_back(val);
        count++;
      }

      if(count == ncol){
        data.push_back(values);
      }
      else continue;

    }

  }

  void Indirectparam_t::fetch_pppc4dmid_data(std::string filename, int ncol)
  {
    auto data = read_file(filename, ncol);

    size_t n_rows = data.size();
    size_t n_cols = data[0].size();

    for (size_t i=0; i<PPPC4DMID.size(); i++) 
    {
      PPPC4DMID[i].resize(n_rows);
      for(size_t j = 0; j < n_rows; j++)
      {
        for(size_t k = 1; k < n_cols; k++) // Skip first column (DM mass)
        {
          PPPC4DMID[i][j].push_back(data[j][k]);
        }
      }
    }

    std::vector<real_t> dm_masses; // Temporary vector to store DM masses
    for(size_t i = 0; i < n_rows; i++)
    {
      dm_masses.push_back(data[i][0]);
    }
    int index = 0;

    while(index<n_rows)
    {
      DM_masses.push_back(dm_masses[index]);
      // Counts the number of rows with the same value of the DM mass
      index += std::count(dm_masses.begin()+index,dm_masses.end(), dm_masses[index])
    }
  }


  std::array<std::vector<real_t>> Indirectparam_t::fermi_energy_bins(real_t mass, std::array<std::vector<real_t>> tab)
  {
    size_t n_rows = PPPC4DMID[0].size();
    size_t n_cols = PPPC4DMID[0][0].size();
    size_t n_masses = DM_masses.size();

    std::array<std::vector<real_t>> final_spec(n_cols);
    
    for (int i=0; i<25; i++) 
    {
      final_spec[i][0]=pow(10.,log10(500)+(log10(500000)-log10(500))/24*i)*0.001; // 24 energy bins logarithmically spaced between 500 and 500000 MeV

      if (final_spec[i][0]>mass) for(j=0; j<n_cols-1; j++) energy_bins[i][j]=0;

      else
      {
        size_t index = 0;

        while(pow(10., tab[index][0]*mass)<=final_spec[i][0] && index < tab.size()) index++;

        real_t K0 = pow(10., tab[index-1][0]*mass);
        real_t K1 = pow(10., tab[index][0]*mass);
        real_t interpolation_factor = (final_spec[i][0]-K0)/(K1-K0);

        for (j=1; j<n_cols; j++) 
        {
          final_spec[i][j]=tab[index-1][j]+(tab[index][j] - tab[index-1][j])*interpolation_factor;
        }
      }
    return final_spec;

    }
  }

  std::array<std::vector<real_t>> Indirectparam_t::interpolate_spectrum_gamma(real_t mass)
  {
    size_t n_rows = PPPC4DMID[0].size();
    size_t n_cols = PPPC4DMID[0][0].size();
    size_t n_masses = DM_masses.size();

    std::array<std::vector<real_t>> interpolated_spectrum(n_cols);// First column is log10(K/m)
    size_t index=0;

    while(mass<DM_masses[index] && index<n_masses) index++;

    if(DM_masses[index] == mass) interpolated_spectrum = fermi_energy_bins(mass, PPPC4DMID[index]);

    else if(DM_masses[index] > mass)
    {
        // Linear interpolation
        std::array<std::vector<real_t>> spec_low = fermi_energy_bins(DM_masses[index-1], PPPC4DMID[index-1]);
        std::array<std::vector<real_t>> spec_high = fermi_energy_bins(DM_masses[index], PPPC4DMID[index]);
        real_t interpolation_factor = (mass - DM_masses[index-1])/(DM_masses[index]-DM_masses[index-1]);

        for(size_t i=0; i<spec_low.size(); i++)
        {
          for(size_t j=0; j<spec_low[0].size(); j++)
          {
            real_t y = spec_low[i][j] + (spec_high[i][j]-spec_low[i][j])*interpolation_factor;
            interpolated_spectrum[i].push_back(y);
          }
        }
    }
      
    
    return interpolated_spectrum;
  }

}; // namespace __SPEC_LIB_NAME__