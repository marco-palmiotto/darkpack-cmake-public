#pragma once


#ifndef INDIRECT__SPEC_LIB_NAME__
  #define INDIRECT__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif


  #include "RunningSM.hpp"
  #include "advmath.hpp"
  #include "avgsvcalculator.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "process.hpp"
  #include "process_1to2.hpp"
  #include <array>
  #include <cctype>
  #include <filesystem>
  #include <fstream>
  #include <iostream>
  #include <sstream>
  #include <vector>

namespace __SPEC_LIB_NAME__
{
  class Process_1to2;        //< Forward declaration of Process_1to2 class
  class Propagation_param_t; //< Forward declaration of Propagation_param_t class


  class Indirect_param_t //< Class containing the parameters and methods for indirect detection
  {
private:
    Param_t input; //!< Input parameters.
    RunningSM run; //!< Handles the running of SM parameters.
    std::shared_ptr<std::vector<Process2to2>>
        p_ptr;                        //!< Pointer to the list of processes (running is external to the processes)
    std::vector<real_t> energy_table; //!< Table containing the different energies for outgoing particles
    std::array<std::vector<real_t>, 14>
        sigma_v_table; //!< Table containing the values of &sigma;*v for each channel at different energies
    std::vector<real_t> sigma_v_process; //!< Vector for the thermally averaged cross section of each process
    real_t total_sigma_v; //!< Value of the total thermally averaged cross section to evaluate the no. of d.o.f.
    real_t dof; //!< Number of degrees of freedom (processes for which sigma_v is at least 1% of the total annihilation
                //!< cross section)
    std::array<std::vector<std::vector<real_t>>, 62> PPPC4DMID;  //!< Data from the PPPC4DMID
    std::vector<real_t> DM_masses;                               //!< List of DM masses in PPPC4DMID
    std::array<std::vector<std::vector<real_t>>, 45> fermi_data; //!< Fermi-LAT binned data
    std::array<real_t, 24> fermi_eflux;            //!< Integrated energy fluxes of gamma rays for each energy bin
    std::vector<std::vector<real_t>> logJ_factors; //!< Log J-factors, errors in logJ and sample value for the 45 dSphs

public:
    /**
     * @brief Construct a new Indirect_param_t object
     *
     * @param input Param_t object containing numerical inputs.
     */
    Indirect_param_t(const Param_t& input_in);

    /**
     * @brief Destructor for the Indirect_param_t class.
     *
     */
    ~Indirect_param_t();


    /**
     * @brief Finds the \f$\chi \overline{\chi}\f$ annihilation processes and appends them to the list.
     *
     */
    void find_annihilation_processes();

    /**
     * @brief Finds the thermally averaged cross section for a 2 to 2 process using a low velocity approximation.
     *
     */
    void compute_sigma_v_at_threshold();

    /**
     * @brief Finds the decay processes present in the 1 to 2 hash table for a given particle
     *
     * @param particle The field number for the initial state particle
     *
     */
    std::vector<Process_1to2> find_decay_processes(int particle);

    /**
     * @brief Handles the decay processes for particles whose hadronization spectrum changes in the model.
     *
     * @param sigma_v The thermally averaged cross section of the process that generated the particle
     *
     * @param particle The particle that is decaying
     *
     * @param sqrt_s The energy of the decaying particle
     */
    void handle_1to2_decays(real_t sigma_v, int particle, real_t sqrt_s);


    /**
     * @brief Fills the spectra for the different annihilation channels.
     *
     *
     */
    void fill_spectrum();

    /**
     * @brief Reads an n-column data file. Returns a vector of vectors, one for each column.
     *
     * @param filename The name of the file to be read.
     *
     * @param ncol The number of columns in the file
     *
     * @param data The vector of vectors to store the data.
     *
     */
    void read_file(std::string filename, int ncol, std::vector<std::vector<real_t>>& data);

    /**
     * @brief Fetches the data from the PPPC4DMID data file.
     *
     */
    void fetch_pppc4dmid_data();

    /**
     * @brief Adapts the spectrum of gamma rays to the energy bins used by Fermi-LAT.
     *
     * @param mass The mass of the dark matter particle.
     *
     * @param tab The tabulated spectrum to be adapted.
     *
     * @param final_spec Dummy vector, stores the adapted spectrum.
     */
    void fermi_energy_bins(const real_t& mass, const std::vector<std::vector<real_t>>& tab,
                           std::vector<std::vector<real_t>>& final_spec);

    /**
     * @brief Interpolates the spectrum of gamma rays for any value of the dark matter mass.
     *
     * @param mass The mass of the dark matter particle.
     *
     * @param interpolated_spectrum Dummy vector, stores the interpolated spectrum.
     */
    void interpolate_spectrum_gamma(const real_t& mass, std::vector<std::vector<real_t>>& interpolated_spectrum);

    /**
     * @brief Obtains the spectrum of \f$\gamma\f$-rays at production from the PPPC4DMID data.
     *
     * @param production_spectrum Dummy vector, stores the spectrum at production.
     *
     */
    void spectrum_at_production(std::vector<std::vector<real_t>>& production_spectrum);

    /**
     * @brief Integrates the spectrum of \f$\gamma\f$-rays between two energy values, based
     *
     * @param spectrum The spectrum of \f$\gamma\f$-rays.
     *
     * @param e_min The lower integration bound.
     *
     * @param e_max The upper integration bound.
     *
     */
    real_t integrate_spectrum(const std::vector<std::vector<real_t>>& spectrum, const real_t& e_min,
                              const real_t& e_max);

    /**
     * @brief Reads dSph data from the Fermi-LAT experiment.
     *
     */
    void read_fermi_data();

    /**
     * @brief Interpolates the log likelihood for a given dSph, energy bin and flux value.
     *
     * @param dsph The index of the dSph.
     *
     * @param bin The index of the energy bin.
     *
     * @param flux The value of the energy flux.
     *
     */
    real_t interpolate_likelihood(const int& dsph, const int& bin, const real_t& flux);

    /**
     * @brief Computes the elements of fermi_eflux, which contain the integrated \f$\gamma\f$-ray spectrum for each of
     * the Fermi-LAT energy bins.
     *
     * @param production_spectrum Dummy vector for which spectrum_at_production is called
     *
     */
    void compute_fluxes(std::vector<std::vector<real_t>>& production_spectrum);

    /**
     * @brief Obtains the log-likelihood for one dwarf spheroidal galaxy.
     *
     * @param dsph The index of the dSph in the Fermi-LAT database.
     *
     * @param logJ Logarithm of the true value of the J-factor.
     *
     */
    real_t likelihood_one_dsph(const int& dsph, const real_t& logJ);
    /**
     * @brief Interpolates a given \f$\bar{p}\f$ spectrum to any value of \f$\log_10(E)=\log(x)\f$.
     *
     * @param spect The \f$\bar{p}\f$ spectrum.
     *
     * @param logx The value of \f$\log_10(E)\f$ to interpolate to.
     *
     */
    real_t logx_interpol(const std::vector<std::vector<real_t>>& spect, real_t logx);

    inline real_t get_total_sigma_v() { return total_sigma_v; };

    inline real_t get_dof() { return dof; };

    inline std::vector<Process2to2> get_processes() { return *p_ptr; };

    inline std::vector<real_t> get_energy_table() { return energy_table; };

    inline std::array<std::vector<real_t>, 14> get_sigma_v_table() { return sigma_v_table; };

    inline std::vector<real_t> get_sigma_v_process() { return sigma_v_process; };

    inline std::vector<real_t> get_DM_masses() { return DM_masses; };

    inline std::array<std::vector<std::vector<real_t>>, 45> get_fermi_data() { return fermi_data; };

    inline std::vector<std::vector<real_t>> get_logJ_factors() { return logJ_factors; };

    inline std::array<real_t, 24> get_fermi_eflux() { return fermi_eflux; };


    /**
     * @brief Changes the value of total_sigma_v. Computes the factor by which total_sigma_v changed, and multiplies all
     * the elements in sigma_v_process and sigma_v_table by that same factor.
     *
     * @param new_total_sigma_v New value for the total \f$\langle \sigma v \rangle\f$.
     */
    inline void change_sigma_v(real_t new_total_sigma_v)
    {
      real_t mass_chi = input.masses_vector.at(input.getLightestBSMpart());
      if (energy_table.size() == 0)
        fill_spectrum();

      real_t factor = new_total_sigma_v / total_sigma_v;
      total_sigma_v = new_total_sigma_v;
      for (size_t i = 0; i < sigma_v_process.size(); i++)
        sigma_v_process.at(i) *= factor;

      for (size_t i = 0; i < sigma_v_table.size(); i++)
        for (size_t j = 0; j < sigma_v_table.at(i).size(); j++)
        {
          sigma_v_table.at(i).at(j) *= factor;
        };
    }

    /**
     * @brief Changes the value of the DM mass. Recomputes total_sigma_v and the elements of sigma_v_process,
     * sigma_v_table and energy_table
     *
     * @param new_mass_DM New value for the DM mass.
     */
    inline void change_mass_DM(real_t new_mass_DM)
    {
      sigma_v_process.clear();
      dof = 0;
      total_sigma_v = 0.;
      input.masses_vector.at(input.getLightestBSMpart()) = new_mass_DM;
      input.m_N_1 = new_mass_DM;
      compute_sigma_v_at_threshold();
      energy_table.clear();
      for (auto& row : sigma_v_table)
        row.clear();
      fill_spectrum();
    };

    /**
     * @brief Removes all elements from sigma_v_table, for all but one channel.
     *
     * @param particle The particle whose channel is to be singled out.
     */
    inline void single_out_annihilation_channel(int particle)
    {
      dof = 1;

      real_t mass_DM = input.masses_vector.at(input.getLightestBSMpart());
      auto it = std::find(energy_table.begin(), energy_table.end(), mass_DM);
      if (it == energy_table.end())
      {
        // Energy not found, add it to the table and update the index
        energy_table.push_back(mass_DM);
        it = std::find(energy_table.begin(), energy_table.end(), mass_DM);
      }
      size_t energy_index = std::distance(energy_table.begin(), it);
      if (energy_index == energy_table.size() - 1)
      {
        for (auto& channel : sigma_v_table)
        {
          channel.resize(energy_index + 1);
        }
      }

      for (auto& row : sigma_v_table)
      {
        for (size_t i = 0; i < row.size(); i++)
          row.at(i) = 0;
      }

      switch (particle)
      {
      case (corr::e):
        sigma_v_table[0][energy_index] = total_sigma_v;
        break;
      case (corr::mu):
        sigma_v_table[1][energy_index] = total_sigma_v;
        break;
      case (corr::tau):
        sigma_v_table[2][energy_index] = total_sigma_v;
        break;
      case (corr::u):
        sigma_v_table[3][energy_index] = total_sigma_v;
        break;
      case (corr::d):
        sigma_v_table[3][energy_index] = total_sigma_v;
        break;
      case (corr::s):
        sigma_v_table[3][energy_index] = total_sigma_v;
        break;
      case (corr::c):
        sigma_v_table[4][energy_index] = total_sigma_v;
        break;
      case (corr::b):
        sigma_v_table[5][energy_index] = total_sigma_v;
        break;
      case (corr::t):
        sigma_v_table[6][energy_index] = total_sigma_v;
        break;
      case (corr::W):
        sigma_v_table[7][energy_index] = total_sigma_v;
        break;
      case (corr::Z):
        sigma_v_table[8][energy_index] = total_sigma_v;
        break;
      case (corr::G):
        sigma_v_table[9][energy_index] = total_sigma_v;
        break;
      case (corr::A):
        sigma_v_table[10][energy_index] = total_sigma_v;
        break;
      case (corr::nu_e):
        sigma_v_table[11][energy_index] = total_sigma_v;
        break;
      case (corr::nu_mu):
        sigma_v_table[12][energy_index] = total_sigma_v;
        break;
      case (corr::nu_tau):
        sigma_v_table[13][energy_index] = total_sigma_v;
        break;
      default:
        break;
      }
    };

  }; // class Indirect_param_t

}; // namespace __SPEC_LIB_NAME__

#endif // INDIRECT__SPEC_LIB_NAME__