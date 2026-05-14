#pragma once

#ifndef PROPAGATION__SPEC_LIB_NAME__
  #define PROPAGATION__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif

  #include "RunningSM.hpp"
  #include "advmath.hpp"
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "indirect.hpp"
  #include "process.hpp"
  #include "process_1to2.hpp"
  #include <array>
  #include <cmath>
  #include <fstream>
  #include <iostream>
  #include <sstream>
  #include <vector>

namespace __SPEC_LIB_NAME__
{
  constexpr const int DIM_TAB_PBAR = 300;
  constexpr const int DIM_TAB_PROTON = 500;
  constexpr const int N_BESSEL = 100;
  constexpr const int N_AMS02 = 57;
  constexpr const int DIM_TAB_PBAR_MFGS = 197;
  constexpr const int N_x_pbar_scan = 300;
  constexpr const double x_pbar_scan_min = 1.0e-9;
  constexpr const double x_pbar_scan_max = 1.0;
  constexpr const int NINT_PRODUCTION = 50000;
  constexpr const double R_GAL = 20.0;
  constexpr const double DENSITY_H_DISC = 0.9;  // in cm^{-3}
  constexpr const double DENSITY_HE_DISC = 0.1; // in cm^{-3}
  constexpr const double E_DISC = 0.1;          // in kpc
  constexpr const double RHO_CHI_0 = 1.0;       // in GeV cm^{-3}
  constexpr const double sigma_v_annihilation = 3.e-26;
  constexpr const double RC_SMBH = 0.1;                  // in kpc
  constexpr const double CELERITY_LIGHT = 2.99792458e10; // in cm s^{-1}
  constexpr const double CM_TO_KPC = 3.0856775807e21;
  constexpr const double SEC_TO_MGYR = 3.15581498e13;
  constexpr const double mb_cm2 = 1.e-27;
  constexpr const double PROTON_MASS = 0.938;              // in GeV
  constexpr const double K_BOLTZMANN = 8.617385e-14;       // in GeV K^{-1}
  constexpr const double H_BAR = 6.58211915e-25;           // in GeV s
  constexpr const double MASSE_ELECTRON = 0.510998928e-03; // in GeV
  constexpr const double RADIUS_ELECTRON = 2.818e-13;      // in cm
  constexpr const double V_ION_H = 19.e-9;                 // in GeV
  constexpr const double V_ION_HE = 44.e-9;                // in GeV
  constexpr const double DENSITY_FREE_ELECTRON = 0.033;    // in cm^{-3}
  constexpr const double T_ELECTRONIC = 3.e5;              // in K
  constexpr const double T_PBAR_MIN = 0.1;                 // in GeV
  constexpr const double T_PBAR_MAX = 1.e5;                // in GeV
  constexpr const double E_PROTON_MIN = 7. * PROTON_MASS;
  constexpr const double E_PROTON_MAX = 1.e6;
  /* Valid for the interval: 3.0-2100 GeV/c */
  constexpr const double p_seuil_pp_tot = 2.83;
  constexpr const double A_pp_tot = 48.0; /* error=+/- 0.1   */
  constexpr const double B_pp_tot = 0.0;
  constexpr const double n_pp_tot = 0.0;
  constexpr const double C_pp_tot = 0.522; /* error=+/- 0.005 */
  constexpr const double D_pp_tot = -4.51; /* error=+/- 0.05  */
  /* Valid for the interval: 0.14-3.0 GeV/c */
  constexpr const double P1_pp_tot = 29.038;  /* error=+/- 0.62113E-01 */
  constexpr const double P2_pp_tot = 37.126;  /* error=+/- 0.16886	*/
  constexpr const double P3_pp_tot = 28.260;  /* error=+/- 0.26031	*/
  constexpr const double P4_pp_tot = -52.708; /* error=+/- 0.30967	*/
  constexpr const double P5_pp_tot = -28.658; /* error=+/- 0.30778	*/
  constexpr const double P6_pp_tot = 16.918;  /* error=+/- 0.16909	*/
  constexpr const double P7_pp_tot = 14.320;  /* error=+/- 0.14092	*/

  constexpr const std::array<std::array<double, 3>, 57> AMS02 = {
      {{4.897275e-01, 7.872882e-03, 1.898842e-03}, {6.179301e-01, 6.977685e-03, 1.064817e-03},
       {7.608659e-01, 1.169011e-02, 1.153023e-03}, {9.219596e-01, 1.227228e-02, 9.959445e-04},
       {1.101684e+00, 1.407454e-02, 1.062232e-03}, {1.299103e+00, 1.541921e-02, 1.039032e-03},
       {1.518860e+00, 1.622801e-02, 1.113851e-03}, {1.760836e+00, 1.748933e-02, 1.050304e-03},
       {2.029303e+00, 1.728564e-02, 1.038072e-03}, {2.324807e+00, 1.764431e-02, 1.028262e-03},
       {2.646625e+00, 1.730250e-02, 9.829095e-04}, {2.999779e+00, 1.637098e-02, 8.797087e-04},
       {3.384188e+00, 1.598068e-02, 9.162533e-04}, {3.804385e+00, 1.458807e-02, 7.769190e-04},
       {4.265516e+00, 1.250476e-02, 6.429843e-04}, {4.767542e+00, 1.135444e-02, 5.459422e-04},
       {5.310399e+00, 9.912291e-03, 4.737686e-04}, {5.898723e+00, 8.772939e-03, 4.090835e-04},
       {6.537690e+00, 7.650446e-03, 3.545111e-04}, {7.227260e+00, 6.583575e-03, 3.028366e-04},
       {7.972125e+00, 5.490503e-03, 2.487684e-04}, {8.777476e+00, 4.691915e-03, 2.031865e-04},
       {9.643280e+00, 3.975649e-03, 1.733090e-04}, {1.058854e+01, 3.240747e-03, 1.418919e-04},
       {1.158648e+01, 2.657462e-03, 1.121182e-04}, {1.263257e+01, 2.235346e-03, 9.241645e-05},
       {1.377894e+01, 1.853769e-03, 8.263010e-05}, {1.502560e+01, 1.492579e-03, 6.093289e-05},
       {1.637254e+01, 1.191751e-03, 5.106521e-05}, {1.781975e+01, 9.541937e-04, 3.790240e-05},
       {1.936721e+01, 7.728250e-04, 2.986475e-05}, {2.101491e+01, 6.335786e-04, 2.379145e-05},
       {2.281080e+01, 5.023920e-04, 1.869613e-05}, {2.475898e+01, 4.112741e-04, 1.456993e-05},
       {2.680734e+01, 3.321899e-04, 1.171139e-05}, {2.900394e+01, 2.681316e-04, 8.548199e-06},
       {3.135270e+01, 2.070874e-04, 6.711036e-06}, {3.384972e+01, 1.750636e-04, 5.387123e-06},
       {3.654696e+01, 1.400439e-04, 4.473537e-06}, {3.944441e+01, 1.100297e-04, 3.606524e-06},
       {4.254206e+01, 8.782044e-05, 3.281007e-06}, {4.583991e+01, 7.051418e-05, 2.766420e-06},
       {4.938612e+01, 5.891023e-05, 2.334929e-06}, {5.318441e+01, 4.820724e-05, 1.980196e-06},
       {5.723106e+01, 3.990519e-05, 1.769411e-06}, {6.157796e+01, 3.090348e-05, 1.486774e-06},
       {6.627329e+01, 2.430237e-05, 1.204277e-06}, {7.132072e+01, 2.060174e-05, 1.063104e-06},
       {7.671660e+01, 1.780130e-05, 1.063092e-06}, {8.558923e+01, 1.220072e-05, 6.403500e-07},
       {9.928355e+01, 7.860344e-06, 4.879763e-07}, {1.152528e+02, 5.090166e-06, 3.764429e-07},
       {1.346166e+02, 3.390081e-06, 2.941159e-07}, {1.594513e+02, 2.080036e-06, 2.332421e-07},
       {1.912190e+02, 1.330016e-06, 2.195476e-07}, {2.328302e+02, 9.110073e-07, 1.697453e-07},
       {3.404216e+02, 2.130008e-07, 5.866878e-08}}};
  constexpr const std::array<std::array<double, 11>, 10> chi2_table = {
      {{0.317, 1.0, 2.30, 3.53, 4.72, 5.89, 7.04, 8.18, 9.30, 10.42, 11.54},
       {0.20, 1.64, 3.22, 4.64, 5.99, 7.29, 8.56, 9.80, 11.03, 12.24, 13.44},
       {0.10, 2.71, 4.61, 6.25, 7.78, 9.24, 10.64, 12.02, 13.36, 14.68, 15.99},
       {0.05, 3.84, 5.99, 7.81, 9.49, 11.07, 12.59, 14.07, 15.51, 16.92, 18.31},
       {0.0455, 4.00, 6.18, 8.02, 9.72, 11.31, 12.85, 14.34, 15.79, 17.21, 18.61},
       {0.01, 6.63, 9.21, 11.34, 13.28, 15.09, 16.81, 18.48, 20.09, 21.67, 23.21},
       {0.0027, 9.0, 11.83, 14.16, 16.25, 18.21, 20.06, 21.85, 23.57, 25.26, 26.90},
       {0.001, 10.83, 13.82, 16.27, 18.47, 20.52, 22.46, 24.32, 26.12, 27.88, 29.59},
       {6.e-5, 16.00, 19.33, 22.06, 24.50, 26.77, 28.91, 30.96, 32.93, 34.95, 36.72},
       {3e-7, 25.00, 28.74, 31.81, 34.56, 37.09, 39.49, 41.78, 43.98, 46.12, 48.19}}};

  class Propagation_param_t //< Class containing the parameters and methods for indirect detection
  {
    using input_function_t = real_t (Propagation_param_t::*)(const std::vector<real_t>&, const std::vector<real_t>&,
                                                             const std::vector<std::vector<real_t>>&);
    using density_profile_t = real_t (*)(real_t, real_t);

public:
    /**
     * @brief Enumeration for defining the propagation model to be used, affecting the 5 free parameters of the
     * diffusion equation.
     */
    enum propagation_models : short int
    {
      MAX = 1,
      MED = 2,
      MIN = 3
    };

    /**
     * @brief Enumeration for defining the dark matter density profile to be used.
     */
    enum halo_profiles : short int
    {
      NFW = 1,
      BURKERT = 2,
      EINASTO = 3
    };

    /**
     * @brief Enumeration for defining which sample of dSph's to use.
     */
    enum sample_options : short int
    {
      CONSERVATIVE = -1,
      NOMINAL = 0,
      INCLUSIVE = 1
    };

private:
    enum sample_options sample_option; //!< Specifies whether the conservative, standard or stringent samples are in use
                                       //!< when analyzing the data.
    std::array<real_t, 5> equation_parameters; //!< Parameters for the propagation equation
    real_t r_earth;                            //!< Distance between the Earth and the galactic center
    real_t rho_chi_local;                      //!< Local density of dark matter
    density_profile_t halo_profile;            //!< Pointer to the dark matter halo profile function
    Param_t input;                             //!< Input parameters.
    Indirect_param_t ind_param;                //!< Indirect detection parameters, for storing the dark matter masses
    std::vector<real_t> alpha_i;               //!< Vector containing the zeros of the Bessel function J0
    std::vector<real_t> q_i;                   //!< Bessel coefficients of the primary cosmic ray flux
    std::vector<std::vector<real_t>> bessel_coef_proton;
    std::vector<std::vector<real_t>>
        bessel_coef_helium; //!< Bessel coefficients of the source term for protons and helium
    std::vector<std::vector<real_t>>
        h_on_h_xsection; //!< Table of differential cross sections for p + p -> pbar + X processes
    std::vector<std::vector<real_t>>
        he_on_h_xsection; //!< Table of differential cross sections for He + p -> pbar + X processes
    std::vector<std::vector<real_t>>
        he_on_he_xsection; //!< Table of differential cross sections for He + He -> pbar + X processes
    std::vector<std::vector<real_t>>
        h_on_he_xsection; //!< Table of differential cross sections for p + He -> pbar + X processes
    std::vector<std::vector<real_t>>
        bessel_pri_pbar_spec; //!< Table containing the bessel coefficients for the spectrum of primary antiprotons
    std::vector<std::vector<real_t>>
        bessel_sec_pbar_spec; //!< Table containing the bessel coefficients for the spectrum of secondary antiprotons
    std::vector<std::vector<real_t>>
        bessel_ter_pbar_spec; //!< Table containing the bessel coefficients for the spectrum of tertiary antiprotons
    std::vector<std::vector<real_t>>
        bessel_tot_pbar_spec; //!< Table containing the bessel coefficients for the total spectrum of antiprotons
    std::vector<std::vector<real_t>> table_abar; //!< Table containing the antiproton A_i coefficients, as they are used
                                                 //!< for computing and evolving all three spectra.
    std::vector<std::vector<real_t>>
        xs_uncertainties; //!< Table containing the antiproton energies, along with their respective lower and upper
                          //!< bounds for the secondary spectra (accounting for cross section uncertainties).
    real_t chi2_noDM;     //!< Chi squared value for the background-only hypothesis (AMS-02)
    std::array<real_t, 45>
        logL_noDM; //!< Log-likelihood values for all dSph for the background-only hypothesis (Fermi-LAT)
    std::vector<real_t> optimal_J_factors; //!< Stores the optimal J-factors after maximizing each dSph's likelihood
public:
    /**
     * @brief Construct a new Propagation_param_t object specifying the propagation and halo models.
     *
     * @param input_in Param_t object containing numerical inputs.
     *
     * @param propagation_model Integer that defines the set of values for the equation parameters. It can be set to
     * 1 (MIN), 2 (MED) or 3 (MAX).
     *
     * @param halo_model Integer that defines the dark matter halo profile. It can be set to 1 (NFW), 2 (Burkert) or 3
     * (Einasto).
     *
     */
    Propagation_param_t(const Param_t& input_in, enum propagation_models p_model, enum halo_profiles halo_model);

    /**
     * @brief Destructor for the Propagation_param_t class.
     *
     */
    ~Propagation_param_t();

    /**
     * @brief NFW DM density profile
     *
     * @param r Radial distance from the center of the galaxy.
     *
     * @param z Vertical coordinate relative to the thin galactic disk.
     *
     */
    static real_t nfw(real_t rr, real_t z);

    /**
     * @brief Burkert DM density profile
     *
     * @param r Radial distance from the center of the galaxy.
     *
     * @param z Vertical coordinate relative to the thin galactic disk.
     */
    static real_t burkert(real_t rr, real_t z);

    /**
     * @brief Einasto DM density profile
     *
     * @param r Radial distance from the center of the galaxy.
     *
     * @param z Vertical coordinate relative to the thin galactic disk.
     */
    static real_t einasto(real_t rr, real_t z);

    /**
     * @brief Auxiliary function to make the likelihood maximization compatible with the function optimization methods.
     *
     * @param logJ Vector containing the logarithm of the true value of the J-factor.
     *
     * @param xtra Vector whose first element is the J-factor uncertainty and the second is the dSph index.
     *
     * @param spect Dummy variable, only present for compatibility with optimization methods.
     *
     */
    real_t max_likelihood_aux_func(const std::vector<real_t>& logJ, const std::vector<real_t>& xtra,
                                   const std::vector<std::vector<real_t>>& spect);

    /**
     * @brief Obtains the maximum likelihood for a given dSph by varying the J-factor.
     *
     * @param dsph The index of the dSph in the Fermi-LAT database.
     *
     * @param logJ_factors Vector of vectors containing the logarithm of the J-factor, its error and the sample
     * (conservative, standard or stringent) for all dSphs.
     *
     * @param IDpowell_result Integer that assesses whether IDpowell ran without issues.
     */
    real_t max_likelihood_one_dsph(const int& dsph, std::vector<std::vector<real_t>>& logJ_factors,
                                   int& IDpowell_result);

    /**
     * @brief Sums the maximum likelihood for all dSphs. Can provide a different answer depending on whether the
     * sample_option is set to -1, 0 or 1 (conservative, standard or stringent).
     *
     */
    real_t likelihood_all_dwarfs(int& IDpowell_result);

    /**
     * @brief Searches for the zeroes of the J0 Bessel function.
     *
     * @param tol The tolerance of the zero search.
     *
     */
    void search_zeroes_J0(const real_t tol);

    /**
     * @brief Returns the distribution of primary sources of cosmic rays in the galactic plane.
     *
     * @param u Radial distance from the galactic center, in terms of the galactic radius.
     *
     */
    real_t f_pulsar_distribution(real_t u);

    /**
     * @brief Obtains the bessel coefficients for the flux of primary cosmic rays, q_i.
     *
     */
    void production();

    /**
     * @brief Obtains the differential flux of interstellar protons, with respect to the proton energy.
     *
     * @param E_proton The proton energy.
     *
     */
    real_t flux_proton_EXP(real_t E_proton);

    /**
     * @brief Obtains the differential flux of interstellar helium, with respect to the energy per nucleon.
     *
     * @param E_nucleon The energy per nucleon of the helium nucleus.
     */
    real_t flux_helium_EXP(real_t E_nucleon);

    /**
     * @brief Obtains the bessel coefficients for the fluxes of proton and helium cosmic rays, for a given energy per
     * nucleon.
     *
     * @param E_nucleon The proton energy/energy per nucleon of the helium nucleus.
     *
     * @param bessel_coef_proton_i Dummy array to be filled with the proton bessel coefficients.
     *
     * @param bessel_coef_helium_i Dummy array to be filled with the helium bessel coefficients.
     */
    void calcul_method_BESSEL_i(real_t E_nucleon, std::vector<real_t>& bessel_coef_proton_i,
                                std::vector<real_t>& bessel_coef_helium_i);

    /**
     * @brief Applies calcul_method_BESSEL_i to all values of the nucleon energy.
     */
    void calculation_BESSEL_Ep_i();

    /**
     * @brief Calculates the Bessel coefficients for the primary antiproton flux.
     *
     * @param n_vert Number of vertical grid points for performing the integration.
     *
     * @param n_rad Number of radial grid points for performing the integration.
     *
     * @param primary_source_term The values of the primary cosmic ray source term.
     */
    void calculation_BESSEL_PBAR_PRIMARY(long int n_vert, long int n_rad, std::vector<real_t>& primary_source_term);

    /**
     * @brief Reads the tabulated primary spectra at production for masses m_inf and m_sup such that
     * m_inf<=mass_chi<=m_sup. Returns a vector containing m_inf, m_sup.
     *
     * @param mass_chi Dark matter mass.
     *
     * @param primary_source_terms_inf 2D vector to store the primary source terms for mass m_inf.
     *
     * @param primary_source_terms_sup 2D vector to store the primary source terms for mass m_sup.
     */
    std::vector<real_t> DNPBAR_ON_DTPBAR_read_file(real_t mass_chi,
                                                   std::vector<std::vector<real_t>>& primary_source_terms_inf,
                                                   std::vector<std::vector<real_t>>& primary_source_terms_sup);

    /**
     * @brief Parses the values from the tabulated primary spectra at production with the cross sections per channel
     * from an Indirect_param_t object.
     *
     * @param dNpbar_on_dEpbar_vec Dummy vector to store the value of dNpbar/dEpbar for all values of the antiproton
     * energy.
     */
    int dNpbar_on_dEpbar_calculation(std::vector<real_t>& dNpbar_on_dEpbar_vec);

    /**
     * @brief Finalizes the primary source term calculation.
     *
     * @param mass_chi Dark matter mass.
     *
     * @param dNpbar_on_dEpbar_vec Dummy vector to store the value of dNpbar/dEpbar for all values of the antiproton
     * energy.
     *
     * @param primary_source_term Dummy vector to store the primary source term.
     */
    void primary_source_calculation(real_t mass_chi, std::vector<real_t>& dNpbar_on_dEpbar_vec,
                                    std::vector<real_t>& primary_source_term);

    /**
     * @brief Computes the primary spectra of antiprotons.
     *
     * @param dNpbar_on_dEpbar_vec Dummy vector to store the value of dNpbar/dEpbar for all values of the antiproton
     * energy.
     *
     * @param primary_source_term Dummy vector to store the primary source term.
     *
     * @param PBAR_SPECTRUM 2D vector to store the primary antiproton spectrum.
     */
    void primary_spectra_BCGS_2014(std::vector<real_t>& dNpbar_on_dEpbar_vec, std::vector<real_t>& primary_source_term,
                                   std::vector<std::vector<real_t>>& PBAR_SPECTRUM);


    /**
     * @brief Calculates the Bessel coefficients for the secondary antiproton flux.
     */
    void calculation_BESSEL_PBAR_SECONDARY();

    /**
     * @brief Calculates the lower and upper bounds of the secondary antiproton flux.
     */
    void preliminary_secondary_spectrum_IS_calculation();

    /**
     * @brief Routine used to account for uncertainties in the secondary spectrum from the antineutron production cross
     * section.
     *
     * @param option Can be set to -1 or 1, to compute the lower and upper bounds of the secondary spectrum,
     * respectively
     */
    void calculation_XS_UNCERTAINTIES_ON_BESSEL_PBAR_SECONDARY(int option);

    /**
     * @brief Calculates the background spectrum of antiprotons for the energy values in logE.
     *
     * @param logE The antiproton energies in log space (the tabulated spectra are interpolated to these values).
     *
     * @param A Real number between 0 and 1 such that the total flux is given by (1-A) * lower + A * upper.
     *
     * @param result 2D vector to store the background antiproton spectrum.
     */
    int background_spectAMS(const std::vector<real_t>& logE, real_t A, std::vector<std::vector<real_t>>& result);

    /**
     * @brief Calculates the Bessel coefficients for the tertiary antiproton flux.
     */
    void calculation_BESSEL_PBAR_TERTIARY();

    /**
     * @brief Calculates the sum of Bessel coefficients for the primary, secondary, and tertiary antiproton flux.
     */
    void calculation_BESSEL_PBAR_SUM();

    /**
     * @brief Computes the total spectrum after solar modulation using a force-field approximation.
     *
     * @param spectrum 2D vector containing the antiproton energy values and the spectrum at Earth as given by the
     * equation of transport.
     *
     * @param phi_f The fisk potential, which parametrizes the force-field approximation.
     *
     * @param result 2D vector containing the antiproton energy values and the spectrum at Earth after solar modulation.
     *
     */
    void solar_mod(std::vector<std::vector<real_t>>& spectrum, real_t phi_f, std::vector<std::vector<real_t>>& result);

    /**
     * @brief Returns the flux of cosmic rays at a given position and energy, computed from a set of bessel
     * coefficients.
     *
     * @param r The radial coordinate, relative to the center of the Milky Way.
     *
     * @param z The vertical coordinate, relative to the galactic disk
     *
     * @param energy The energy of the cosmic ray.
     *
     * @param mass The mass of the cosmic ray particle.
     *
     * @param Z_em The electric charge of the cosmic ray particle (in multiples of e).
     *
     * @param bes_coef Vector containing the bessel coefficients.
     *
     */
    real_t GENERIC_FLUX(real_t r, real_t z, real_t energy, real_t mass, real_t Z_em, std::vector<real_t> bes_coef);

    /**
     * @brief Returns the flux of cosmic rays at a given position and energy, computed from a set of bessel
     * coefficients. Modulates the result with Torsten Bringman's coefficients to speed up convergence.
     *
     * @param r The radial coordinate, relative to the center of the Milky Way.
     *
     * @param z The vertical coordinate, relative to the galactic disk
     *
     * @param energy The energy of the cosmic ray.
     *
     * @param mass The mass of the cosmic ray particle.
     *
     * @param Z_em The electric charge of the cosmic ray particle (in multiples of e).
     *
     * @param bes_coef Vector containing the bessel coefficients.
     *
     */
    real_t GENERIC_FLUX_04(real_t r, real_t z, real_t enerçgy, real_t mass, real_t Z_em, std::vector<real_t> bes_coef);

    /**
     * @brief Solves the diffusion equation describing the energy behaviour of the Bessel transforms for all values of
     * the antiproton energy.
     */
    void calculation_BESSEL_PBAR_direct_inversion();

    /**
     * @brief This routine solves the matrix equation [A] * [u] = [r] in the case of a tridiagonal matrix [A]. Pivoting
     * is not used.
     *
     * @param a The upper diagonal of the A matrix.
     *
     * @param b The main diagonal of the A matrix.
     *
     * @param c The lower diagonal of the A matrix.
     *
     * @param u The vector in the left-hand side of the matrix equation.
     *
     * @param r The vector in the right-hand side of the matrix equation.
     */
    void inversion_tridiagonal(std::array<real_t, DIM_TAB_PBAR + 1>& a, std::array<real_t, DIM_TAB_PBAR + 1>& b,
                               std::array<real_t, DIM_TAB_PBAR + 1>& c, std::array<real_t, DIM_TAB_PBAR + 1>& r,
                               std::array<real_t, DIM_TAB_PBAR + 1>& u);

    /**
     * @brief Computes the spatial diffusion term in the diffusion equation.
     *
     * @param energy The energy of the cosmic ray.
     *
     * @param mass The mass of the cosmic ray particle.
     *
     * @param Z_em The electric charge of the cosmic ray particle (in multiples of e).
     */
    real_t K_space_diffusion(real_t energy, real_t mass, real_t Z_em);

    /**
     * @brief Computes the energy diffusion term in the diffusion equation.
     *
     * @param energy The energy of the cosmic ray.
     *
     * @param mass The mass of the cosmic ray particle.
     *
     * @param Z_em The electric charge of the cosmic ray particle (in multiples of e).
     */
    real_t D_energy_diffusion(real_t energy, real_t mass, real_t Z_em);

    /**
     * @brief Computes the term in the diffusion equation which accounts for energy losses from ionization, as well as
     * coulomb friction on the gas of the galactic ridge.
     *
     * @param energy The energy of the cosmic ray.
     *
     * @param mass The mass of the cosmic ray particle.
     *
     * @param Z_em The electric charge of the cosmic ray particle (in multiples of e).
     */
    real_t b_energy_losses(real_t energy, real_t mass, real_t Z_em);

    /**
     * @brief Iteratively applies calculation_BESSEL_PBAR_TOT_direct_inversion_A to compute the effects of energy losses
     * and diffusive reacceleration.
     */
    void ELDR_effect_calculation();

    /**
     * @brief Computes the total cross section for inelastic collisions between proton cosmic rays and interstellar
     * hydrogen.
     *
     * @param E_proton The proton energy.
     */
    real_t sigma_inelastic_pH_TAN_and_NG(real_t E_proton);

    /**
     * @brief Computes the total cross section for collisions between proton cosmic rays and interstellar hydrogen.
     *
     * @param E_proton The proton energy.
     */
    real_t sigma_total_pH(real_t E_proton);

    /**
     * @brief Computes the total cross section for inelastic collisions between antiproton cosmic rays and interstellar
     * hydrogen.
     *
     * @param E_pbar The antiproton energy.
     */
    real_t sigma_inelastic_pbarH_TAN_and_NG(real_t E_pbar);

    /**
     * @brief Computes the total cross section for non-annihilating inelastic collisions between proton cosmic rays and
     * interstellar hydrogen.
     *
     * @param E_pbar The antiproton energy.
     */
    real_t sigma_inelastic_NOANN_pbarH_TAN_and_NG(real_t E_pbar);


    //--------------Optimization Methods------------------//

    /**
     * @brief Optimizes a function using Powell's method.
     *
     * @param n The number of parameters to optimize.
     *
     * @param func Pointer to the member function to be optimized. Must take as input parameters two 1D real_t vectors
     * and one 2D real_t vector, and return a real_t.
     *
     * @param xtra Vector containing extra (constant) parameters to be passed to func.
     *
     * @param spect 2D vector containing antiproton spectral information.
     *
     * @param xlim_min Vector containing the lower limits for each of the parameters.
     *
     * @param xlim_max Vector containing the upper limits for each of the parameters.
     *
     * @param fmin Pointer to store the extremum value of the function.
     *
     * @param xmin Vector to store the values of the parameters that extremize the function.
     *
     * @param ftol Tolerance for the optimization.
     *
     * @param option String to specify whether to minimize or maximize the function.
     */
    int IDpowell(int n, input_function_t func, const std::vector<real_t>& xtra,
                 const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& xlim_min,
                 const std::vector<real_t>& xlim_max, real_t* fmin, std::vector<real_t>& xmin, real_t ftol,
                 const std::string& option);


    /**
     * @brief Called by IDpowell to apply Powell's method.
     *
     * @param n The number of parameters to optimize.
     *
     * @param func Pointer to the member function to be optimized. Must take as input parameters two 1D real_t vectors
     * and one 2D real_t vector, and return a real_t.
     *
     * @param xtra Vector containing extra (constant) parameters to be passed to func.
     *
     * @param spect 2D vector containing antiproton spectral information.
     *
     * @param xlim_min Vector containing the lower limits for each of the parameters.
     *
     * @param xlim_max Vector containing the upper limits for each of the parameters.
     *
     * @param fmin Pointer to store the extremum value of the function.
     *
     * @param xmin Vector to store the values of the parameters that extremize the function.
     *
     * @param ftol Tolerance for the optimization.
     *
     * @param factor Set to 1 for minimization and -1 for maximization.
     */
    int IDpowellaux(int n, input_function_t func, const std::vector<real_t>& xtra,
                    const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& x0,
                    const std::vector<real_t>& xlim_min, const std::vector<real_t>& xlim_max, real_t* fmin,
                    std::vector<real_t>& xmin, real_t ftol, real_t factor);

    /**
     * @brief Applies Brent's method at each step of IDpowellaux to optimize the function along a given direction.
     *
     * @param n The number of parameters to optimize.
     *
     * @param func Pointer to the member function to be optimized. Must take as input parameters two 1D real_t vectors
     * and one 2D real_t vector, and return a real_t.
     *
     * @param xtra Vector containing extra (constant) parameters to be passed to func.
     *
     * @param spect 2D vector containing antiproton spectral information.
     *
     * @param xmin Vector to store the values of the parameters that extremize the function.
     *
     * @param xi Direction along which to optimize the function.
     *
     * @param xlim_min Vector containing the lower limits for each of the parameters.
     *
     * @param xlim_max Vector containing the upper limits for each of the parameters.
     *
     * @param fmin Pointer to store the extremum value of the function.
     *
     * @param factor Set to 1 for minimization and -1 for maximization.
     */
    int IDbrentmethod(int n, input_function_t func, const std::vector<real_t>& xtra,
                      const std::vector<std::vector<real_t>>& spect, std::vector<real_t>& xmin,
                      const std::vector<real_t>& xi, const std::vector<real_t>& xlim_min,
                      const std::vector<real_t>& xlim_max, real_t* fmin, real_t factor);

    /**
     * @brief Auxiliary method called by IDbrentmethod1D. Returns the value of the function along a given direction xi
     * at a given displacement x.
     *
     * @param func Pointer to the member function to be optimized. Must take as input parameters two 1D real_t vectors
     * and one 2D real_t vector, and return a real_t.
     *
     * @param xtra Vector containing extra (constant) parameters to be passed to func.
     *
     * @param spect 2D vector containing antiproton spectral information.
     *
     * @param x Displacement along the direction xi.
     *
     * @param xmini Vector containing the best found values of the parameters so far.
     *
     * @param xi Direction along which to optimize the function.
     *
     * @param n The number of parameters to optimize.
     */
    real_t IDbrentmethod1Dfunc(input_function_t func, const std::vector<real_t>& xtra,
                               const std::vector<std::vector<real_t>>& spect, real_t x,
                               const std::vector<real_t>& xmini, const std::vector<real_t>& xi, int n);


    /**
     * @brief Auxiliary method called by IDbrentmethod to optimize the function along a given direction using Brent's
     * method.
     *
     * @param func Pointer to the member function to be optimized. Must take as input parameters two 1D real_t vectors
     * and one 2D real_t vector, and return a real_t.
     *
     * @param xtra Vector containing extra (constant) parameters to be passed to func.
     *
     * @param spect 2D vector containing antiproton spectral information.
     *
     * @param a0 Lower limit of the interval to search for the extremum.
     *
     * @param x0 Initial guess for the optimal displacement along the direction xi.
     *
     * @param b0 Upper limit of the interval to search for the extremum.
     *
     * @param xmini Vector containing the best found values of the parameters so far.
     *
     * @param xi Direction along which to optimize the function.
     *
     * @param n The number of parameters to optimize.
     */
    int IDbrentmethod1D(input_function_t func, const std::vector<real_t>& xtra,
                        const std::vector<std::vector<real_t>>& spect, real_t a0, real_t x0, real_t b0, real_t* fmin,
                        real_t* xmin, real_t factor, const std::vector<real_t>& xmini, std::vector<real_t> xi, int n);


    /**
     * @brief Finds the optimal values of a0, b0 and x0 to bracket the extremum of the function along a given direction.
     *
     * @param n The number of parameters to optimize.
     *
     * @param func Pointer to the member function to be optimized. Must take as input parameters two 1D real_t vectors
     * and one 2D real_t vector, and return a real_t.
     *
     * @param xtra Vector containing extra (constant) parameters to be passed to func.
     *
     * @param spect 2D vector containing antiproton spectral information.
     *
     * @param a0 Pointer to the value of a0.
     *
     * @param x0 Pointer to the value of x0.
     *
     * @param b0 Pointer to the value of b0.
     *
     * @param xinit Vector containing the initial guess for the parameter values.
     *
     * @param xi Direction along which to optimize the function.
     *
     * @param xlim_min Vector containing the lower limits for each of the parameters.
     *
     * @param xlim_max Vector containing the upper limits for each of the parameters.
     *
     * @param factor Set to 1 for minimization and -1 for maximization.
     */
    void IDbracket(int n, input_function_t func, const std::vector<real_t>& xtra,
                   const std::vector<std::vector<real_t>>& spect, real_t* a0, real_t* x0, real_t* b0,
                   std::vector<real_t> xinit, std::vector<real_t> xi, const std::vector<real_t>& xlim_min,
                   const std::vector<real_t>& xlim_max, real_t factor);
    /**
     * @brief Adds two spectra together. Returns 0 if the spectra are not of the same size.
     *
     * @param a1 The first spectrum to be added.
     *
     * @param a2 The second spectrum to be added. The values in a2 are interpolated to be at the same energies as a1.
     *
     * @param add Dummy vector to store the resulting summed spectrum.
     */
    void add_spectra(const std::vector<std::vector<real_t>>& a1, const std::vector<std::vector<real_t>>& a2,
                     std::vector<std::vector<real_t>>& add);

    /**
     * @brief Returns the chi squared from the AMS-02 experiment in the background-only hypothesis.
     *
     * @param param Vector of nuisance parameters, A and \f$\phi_f\f$ (uncertainties from cross sections and solar modulation,
     * respectively).
     *
     * @param logE The antiproton energies in log space (the tabulated spectra are interpolated to these values).
     *
     * @param spect Dummy variable, only present for compatibility with optimization methods.
     */
    real_t dchi_bkg(const std::vector<real_t>& param, const std::vector<real_t>& logE,
                    const std::vector<std::vector<real_t>>& spect);

    /**
     * @brief Finds the combination of A and \f$\phi_f\f$ that minimizes the chi squared in the background only hypothesis.
     * Stores the result in chi2_noDM.
     *
     * @param IDpowell_result Integer that assesses whether IDpowell ran without issues.
     */
    std::vector<real_t> chi2_bkg(int& IDpowell_result);

    /**
     * @brief Compares a given spectrum with the values tabulated by AMS-02 to return its respective chi squared value.
     *
     * @param spectrum 2D vector containing the antiproton energy values and the spectrum at Earth after solar
     * modulation.
     */
    real_t chi2_AMS(const std::vector<std::vector<real_t>>& spectrum);

    /**
     * @brief Returns the chi squared for the sum of the primary and secondary spectrum.
     *
     * @param param Vector containing the values of A and the Fisk potential.
     *
     * @param logE The antiproton energies of the background spectrum in log space (the primary spectrum is interpolated
     * to these values).
     *
     * @param spect 2D vector containing the antiproton energy values and the primary spectrum at Earth before solar
     * modulation.
     */
    real_t dchi_tot(const std::vector<real_t>& param, const std::vector<real_t>& logE,
                    const std::vector<std::vector<real_t>>& spect);

    /**
     * @brief Finds the values of A and \f$\phi_f\f$ that minimize the \f$\chi^2\f$ for the total spectrum. Subtracts the \f$\chi^2\f$
     * of the background-only hypothesis to return deltachi2_AMS.
     * @param IDpowell_result Integer that assesses whether IDpowell ran without issues.
     */
    real_t deltachi2_AMS(int& IDpowell_result, std::vector<real_t>& xmin);

    /**
     * @brief Computes the delta log-likelihood for the Fermi-LAT data.
     *
     * @param IDpowell_result Integer that assesses whether IDpowell ran without issues.
     */
    real_t deltalikelihood_fermi(int& IDpowell_result);

    inline std::vector<real_t> get_alpha_i() const { return alpha_i; };

    inline std::vector<real_t> get_q_i() const { return q_i; };

    inline std::vector<std::vector<real_t>> get_bessel_coef_proton() const { return bessel_coef_proton; };

    inline std::vector<std::vector<real_t>> get_bessel_coef_helium() const { return bessel_coef_helium; };

    inline std::vector<std::vector<real_t>> get_xs_uncertainties() const { return xs_uncertainties; };

    inline real_t get_chi2_noDM() const { return chi2_noDM; };

    inline std::vector<real_t> get_optimal_J_factors() const { return optimal_J_factors; };

    inline std::array<real_t, 24> get_fermi_eflux() { return ind_param.get_fermi_eflux(); };

    inline std::vector<real_t> get_DM_masses() { return ind_param.get_DM_masses(); };

    inline std::vector<std::vector<real_t>> get_logJ_factors() { return ind_param.get_logJ_factors(); };

    inline real_t get_total_sigma_v() { return ind_param.get_total_sigma_v(); };

    inline int get_dof() { return ind_param.get_dof(); };

    inline std::vector<Process2to2> get_processes() { return ind_param.get_processes(); };

    inline std::vector<real_t> get_sigma_v_process() { return ind_param.get_sigma_v_process(); };

    inline std::array<std::vector<real_t>, 14> get_sigma_v_table() { return ind_param.get_sigma_v_table(); };

    inline std::vector<real_t> get_energy_table() { return ind_param.get_energy_table(); };


    /**
     * @brief Changes the private Indirect_param_t member.
     *
     * @param new_sigma_v New value for the total \f$\langle \sigma v \rangle\f$. If set to zero, this parameter is not changed.
     * 
     * @param new_mass_DM New value for the DM mass. If set to zero, this parameter is not changed.
     */
    inline void change_ind_param(real_t new_sigma_v, real_t new_mass_DM)
    {
      if (new_mass_DM != 0.)
      {
        ind_param.change_mass_DM(new_mass_DM);
        input.masses_vector.at(input.getLightestBSMpart()) = new_mass_DM;
        input.m_N_1 = new_mass_DM;
      }

      auto energy_table = ind_param.get_energy_table();
      auto sigma_v_process = ind_param.get_sigma_v_process();
      auto processes = ind_param.get_processes();
      int DM_candidate = input.getLightestBSMpart();
      real_t mass_chi = input.masses_vector.at(DM_candidate);

      if (new_sigma_v != 0.)
      {
        ind_param.change_sigma_v(new_sigma_v);
      }
    }

    /**
     * @brief Removes all elements from sigma_v_table in the Indirect_param_t member, for all but one channel.
     *
     * @param particle The particle whose channel is to be singled out.
     */
    inline void single_out_annihilation_channel(int particle) { ind_param.single_out_annihilation_channel(particle); };


    real_t test_func(const std::vector<real_t>& x, const std::vector<real_t>& xtra,
                     const std::vector<std::vector<real_t>>& spect);
    real_t test_func_spec(const std::vector<real_t>& x, const std::vector<real_t>& xtra,
                          const std::vector<std::vector<real_t>>& spect);

  }; // class Propagation_param_t
} // namespace __SPEC_LIB_NAME__

#endif // PROPAGATION__SPEC_LIB_NAME__