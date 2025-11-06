#pragma once
#ifndef __cplusplus
  #define __cplusplus
#endif

#ifndef INDIRECT__SPEC_LIB_NAME__
  #define INDIRECT__SPEC_LIB_NAME__

  #include "RunningSM.hpp"
  #include "advmath.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "process.hpp"
  #include "process_1to2.hpp"
  #include <array>
  #include <vector>
  #include <fstream>
  #include <sstream>
  #include <iostream>

namespace __SPEC_LIB_NAME__
{
  class Process_1to2;   //< Forward declaration of Process_1to2 class
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
    std::array<std::array<std::vector<real_t>>>(62) PPPC4DMID; //!< Data from the PPPC4DMID
    std::vector<real_t> DM_masses;                               //!< List of DM masses in PPPC4DMID

public:
    /**
     * @brief Construct a new Indirectparam_t object
     *
     * @param input Param_t object containing numerical inputs.
     */
    Indirectparam_t(const Param_t& input);

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
     * @param sqrtS The center-of-mass energy.
     *
     */
    void fill_spectrum(std::vector<Process2to2> processes, std::vector<real_t> sigma_v_vector, double sqrt_s);

    /**
     * @brief Fills the spectra for the different annihilation channels.
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
     */
    std::vector<std::vector<double>> read_file(std::string filename, int ncol);

    /**
     * @brief Creates a list of masses from a data file.
     *
     * @param filename The name of the file to be read.
     * 
     * @param ncol The number of columns in the file
     *
     */
    void fetch_pppc4dmid_data(std::string filename, int ncol);

    /**
     * @brief Adapts the spectrum of gamma rays to the energy bins used by Fermi-LAT.
     * 
     * @param mass The mass of the dark matter particle.
     * 
     * @param tab The tabulated spectrum to be adapted.
     */
    std::array<std::vector<real_t>> fermi_energy_bins(real_t mass, std::array<std::vector<real_t>> tab);

     /**
     * @brief Interpolates the spectrum of gamma rays for any value of the dark matter mass.
     *
     * @param mass The mass of the dark matter particle.
     * 
     */
    std::array<std::vector<real_t>> interpolate_spectrum_gamma(real_t mass);

    

  }; // class Indirectparam_t




}; // namespace __SPEC_LIB_NAME__
#endif // RELICPARAM__SPEC_LIB_NAME__