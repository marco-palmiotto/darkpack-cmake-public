#pragma once


#ifndef INDIRECT__SPEC_LIB_NAME__
  #define INDIRECT__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif


  #include "RunningSM.hpp"
  #include "advmath.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "process.hpp"
  #include "process_1to2.hpp"
  #include <array>
  #include <filesystem>
  #include <fstream>
  #include <iostream>
  #include <sstream>
  #include <vector>

namespace __SPEC_LIB_NAME__
{
  class Process_1to2;        //< Forward declaration of Process_1to2 class
  class Propagation_param_t; //< Forward declaration of Propagation_param_t class


  class Indirectparam_t //< Class containing the parameters and methods for indirect detection
  {
private:
    Param_t input;                                   //!< Input parameters.
    RunningSM run;                                   //!< Handles the running of SM parameters.
    std::shared_ptr<std::vector<Process2to2>> p_ptr; //!< Pointer to the list of processes
    std::vector<real_t> energy_table;                //!< Table containing the different energies for outgoing particles
    std::array<std::vector<real_t>, 14>
        sigma_v_table; //!< Table containing the values of &sigma;*v for each channel at different energies
    std::vector<real_t> sigma_v_process; //!< Vector for the thermally averaged cross section of each process
    real_t total_sigma_v; //!< Value of the total thermally averaged cross section to evaluate the no. of d.o.f.
    real_t dof;           //!< Number of degrees of freedom (processes for which sigma_v is at least 1% of the total)
    std::array<std::vector<std::vector<real_t>>, 62> PPPC4DMID;  //!< Data from the PPPC4DMID
    std::vector<real_t> DM_masses;                               //!< List of DM masses in PPPC4DMID
    std::array<std::vector<std::vector<real_t>>, 47> fermi_data; //!< Fermi-LAT binned data
    std::vector<std::vector<real_t>> logJ_factors;               //!< Log J-factors for the 47 dSphs

public:
    /**
     * @brief Construct a new Indirectparam_t object
     *
     * @param input Param_t object containing numerical inputs.
     */
    Indirectparam_t(const Param_t& input_in);

    /**
     * @brief Destructor for the Indirectparam_t class.
     *
     */
    ~Indirectparam_t();


    /**
     * @brief Finds the \f$\chi \overline{\chi}\f$ annihilation processes and appends them to the list.
     *
     */
    void find_annihilation_processes();

    /**
     * @brief Finds the thermally averaged cross section for a 2 to 2 process using a low velocity approximation.
     *
     */
    void sigma_v(std::vector<Process2to2> Processes);

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
     */
    void handle_1to2_decays(real_t sigma_v, int particle, real_t sqrt_s);


    /**
     * @brief Fills the spectra for the different annihilation channels.
     *
     * @param processes The list of 2 to 2 processes.
     *
     * @param sigma_v_vector The vector containing the thermally averaged cross sections for each process.
     *
     * @param sqrtS The center-of-mass energy.
     *
     */
    void fill_spectrum(std::vector<Process2to2> processes, std::vector<real_t> sigma_v_vector, double sqrt_s);

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
     * @param final_spec The adapted spectrum.
     */
    void fermi_energy_bins(const real_t& mass, const std::vector<std::vector<real_t>>& tab,
                           std::vector<std::vector<real_t>>& final_spec);

    /**
     * @brief Interpolates the spectrum of gamma rays for any value of the dark matter mass.
     *
     * @param mass The mass of the dark matter particle.
     *
     * @param interpolated_spectrum The interpolated spectrum.
     */
    void interpolate_spectrum_gamma(const real_t& mass, std::vector<std::vector<real_t>>& interpolated_spectrum);

    /**
     * @brief Obtains the spectrum of gamma rays at production from the PPPC4DMID data.
     *
     * @param production_spectrum The spectrum at production.
     *
     */
    void spectrum_at_production(std::vector<std::vector<real_t>>& production_spectrum);

    /**
     * @brief Obtains the spectrum of gamma rays at production from the PPPC4DMID data.
     *
     * @param production_spectrum The spectrum at production.
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
     * @brief Obtains the spectrum of gamma rays at production from the PPPC4DMID data.
     *
     * @param production_spectrum The spectrum at production.
     *
     */
    real_t likelihood_one_dsph(const int& dsph, const real_t& logJ);

    inline std::vector<std::vector<real_t>> get_logJ() { return logJ_factors; };



  }; // class Indirectparam_t

}; // namespace __SPEC_LIB_NAME__

#endif // INDIRECT__SPEC_LIB_NAME__