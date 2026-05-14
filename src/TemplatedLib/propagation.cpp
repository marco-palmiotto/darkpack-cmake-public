#include "propagation.hpp"
// #define DEBUG

using namespace advmath;


template <typename T> static inline T SQUARE(const T x) { return x * x; }

template <typename T> static inline T SQRT(const T x) { return pow(x, 0.5); }

namespace __SPEC_LIB_NAME__
{

  Propagation_param_t::Propagation_param_t(const Param_t& input_in, enum propagation_models p_model,
                                           enum halo_profiles halo_model)
          : input(input_in), ind_param(input_in)
  {
    // Fetch the dark matter masses from the PPPC4DMID data
    ind_param.fetch_pppc4dmid_data();
    ind_param.read_fermi_data();

    // Fetch the cross section data for antiproton production
    ind_param.read_file(FILENAME_H_ON_H, DIM_TAB_PBAR + 1, h_on_h_xsection);
    ind_param.read_file(FILENAME_HE_ON_H, DIM_TAB_PBAR + 1, he_on_h_xsection);
    ind_param.read_file(FILENAME_HE_ON_HE, DIM_TAB_PBAR + 1, he_on_he_xsection);
    ind_param.read_file(FILENAME_H_ON_HE, DIM_TAB_PBAR + 1, h_on_he_xsection);

    switch (p_model)
    {
    case (propagation_models::MAX):
      equation_parameters = {0.0765 * pow(CM_TO_KPC, 2.) / SEC_TO_MGYR, 0.46, 15.0, 5.0 * 1.0e5, 117.6 * 1.0e5};
      break;
    case (propagation_models::MED):
      equation_parameters = {0.0112 * pow(CM_TO_KPC, 2.) / SEC_TO_MGYR, 0.7, 4.0, 12.0 * 1.0e5, 52.9 * 1.0e5};
      break;
    case (propagation_models::MIN):
      equation_parameters = {0.0016 * pow(CM_TO_KPC, 2.) / SEC_TO_MGYR, 0.85, 1.0, 13.5 * 1.0e5, 22.4 * 1.0e5};
      break;
    default:
      std::cerr << "Invalid propagation model specified. Using 'med' parameters by default." << std::endl;
      equation_parameters = {0.0112 * pow(CM_TO_KPC, 2.) / SEC_TO_MGYR, 0.7, 4.0, 12.0 * 1.0e5, 52.9 * 1.0e5};
    }

    switch (halo_model)
    {
    case (halo_profiles::NFW):
      halo_profile = &nfw;
      r_earth = 8.21;
      rho_chi_local = 0.383;
      break;
    case (halo_profiles::BURKERT):
      halo_profile = &burkert;
      r_earth = 7.94;
      rho_chi_local = 0.487;
      break;
    case (halo_profiles::EINASTO):
      halo_profile = &einasto;
      r_earth = 8.25;
      rho_chi_local = 0.386;
      break;
    default:
      std::cerr << "Invalid halo model specified. Using NFW profile by default." << std::endl;
      halo_profile = &nfw;
    }

    sample_option = sample_options::NOMINAL; // default to the nominal dSph sample unless a specific propagation model
                                             // and halo profile is used

    if (p_model == propagation_models::MED && halo_model == halo_profiles::BURKERT)
      sample_option = sample_options::CONSERVATIVE; // Conservative sample is in use
    if (p_model == propagation_models::MED && halo_model == halo_profiles::EINASTO)
      sample_option = sample_options::NOMINAL; // nominal sample is in use
    if (p_model == propagation_models::MAX && halo_model == halo_profiles::EINASTO)
      sample_option = sample_options::INCLUSIVE; // inclusive sample is in use


    alpha_i.resize(N_BESSEL);
    q_i.resize(N_BESSEL);

    search_zeroes_J0(1e-8);
    production();

    bessel_coef_proton.resize(DIM_TAB_PROTON + 1);
    bessel_coef_helium.resize(DIM_TAB_PROTON + 1);
    /* pt_Proton->BESSEL_COEF_Enuc_i is set to zero. */
    for (int i_nuc = 0; i_nuc <= DIM_TAB_PROTON; i_nuc++)
    {
      bessel_coef_proton.at(i_nuc).resize(N_BESSEL);
      bessel_coef_helium.at(i_nuc).resize(N_BESSEL);
    }

    calculation_BESSEL_Ep_i();

    table_abar.resize(DIM_TAB_PBAR + 1);

    bessel_pri_pbar_spec.resize(DIM_TAB_PBAR + 1);
    bessel_sec_pbar_spec.resize(DIM_TAB_PBAR + 1);
    bessel_ter_pbar_spec.resize(DIM_TAB_PBAR + 1);
    bessel_tot_pbar_spec.resize(DIM_TAB_PBAR + 1);

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* bessel_sec_pbar_spec and TABLE_Abar_i are set to zero. */
    {
      table_abar.at(i_pbar).resize(N_BESSEL);

      bessel_pri_pbar_spec.at(i_pbar).resize(N_BESSEL);
      bessel_sec_pbar_spec.at(i_pbar).resize(N_BESSEL);
      bessel_ter_pbar_spec.at(i_pbar).resize(N_BESSEL);
      bessel_tot_pbar_spec.at(i_pbar).resize(N_BESSEL);
    }

    chi2_noDM = 0;
  };


  //-------------------------------DARK MATTER DENSITY PROFILES----------------------------//
  real_t Propagation_param_t::nfw(real_t rr, real_t z)
  {
    real_t r = sqrt(rr * rr + z * z);
    real_t alpha = 1.0; // [NO UNIT]
    real_t beta = 3.0;  // [NO UNIT]
    real_t gamma = 1.0; // [NO UNIT]
    real_t core = 19.6; // [kpc]

    real_t r_earth = 8.21; // [kpc]
    real_t rho_chi_solar = 0.383;

    real_t enhancement = 1.0;

    real_t a2 = 8.0 * gamma * (pow(M_PI, (2.0)) - 9.0 + 6.0 * gamma) / (9.0 * (3.0 - 2.0 * gamma));
    real_t a1 = a2 + (2.0 * gamma);
    real_t x;

    if (r <= RC_SMBH)
    {
      x = r / RC_SMBH;

      if (x <= 1.e-3)
        enhancement = a1 + a2;
      else
        enhancement = (a1 * sin(M_PI * x) / (M_PI * x)) + (a2 * sin(2.0 * M_PI * x) / (2.0 * M_PI * x));

      enhancement += 1.;
      enhancement = sqrt(enhancement);

      r = RC_SMBH;
    }

    // Generic Profile
    real_t rho_chi = (1. + pow((r_earth / core), alpha)) / (1. + pow((r / core), alpha));
    rho_chi = rho_chi_solar * pow((r_earth / r), gamma) * pow(rho_chi, ((beta - gamma) / alpha)); // [GeV cm^{-3}]
    rho_chi *= enhancement;
    rho_chi /= RHO_CHI_0;

    return rho_chi;
  }

  real_t Propagation_param_t::burkert(real_t rr, real_t z)
  {
    real_t core = 9.26;
    real_t r = sqrt(rr * rr + z * z);
    real_t r_earth = 7.94; // [kpc]
    real_t rho_chi_solar = 0.487;
    real_t rho_chi = (1.0 + r / core) * (1.0 + pow(r / core, 2.0));
    rho_chi = rho_chi_solar * (1.0 + r_earth / core) * (1.0 + pow(r_earth / core, 2.0)) / rho_chi;
    rho_chi /= RHO_CHI_0;

    return rho_chi;
  }

  real_t Propagation_param_t::einasto(real_t rr, real_t z)
  {
    real_t alpha = 0.22; // [NO UNIT] 			// Catena & Ullio (2010)
    real_t core = 16.07; // [kpc]				// Catena & Ullio (2010)
    real_t rho_chi;
    real_t r = sqrt(rr * rr + z * z);
    real_t r_earth = 8.25;
    real_t rho_chi_solar = 0.386;
    rho_chi = pow(r_earth / core, alpha) - pow(r / core, alpha);
    rho_chi = rho_chi_solar * exp((2.0 / alpha) * rho_chi);
    rho_chi /= RHO_CHI_0;

    return rho_chi;
  }

  //-------------------------------LIKELIHOOD METHODS FOR FERMI-LAT DSphs----------------------------//

  real_t Propagation_param_t::max_likelihood_aux_func(const std::vector<real_t>& logJ, const std::vector<real_t>& xtra,
                                                      const std::vector<std::vector<real_t>>& spect)
  {
    // Auxiliary function to be used in the maximization of the likelihood (needed redefinition to make the optimization
    // compatible with both fermi-LAT and AMS-02)
    real_t sigma_j = xtra.at(0);
    auto dsph = static_cast<int>(xtra.at(1));

    real_t logJ0 = logJ.at(0);

    return ind_param.likelihood_one_dsph(dsph, logJ0);
  }

  real_t Propagation_param_t::max_likelihood_one_dsph(const int& dsph, std::vector<std::vector<real_t>>& logJ_factors,
                                                      int& IDpowell_result)
  {
    // Calls max_likelihood_aux_func and returns the maximum log likelihood for one dSph by optimizing over logJ
    real_t logJ_obs = logJ_factors.at(dsph).at(0);
    real_t sigma_j = logJ_factors.at(dsph).at(1);

    std::vector<real_t> logJ_opt;
    logJ_opt.resize(1);
    real_t lhmax;

    std::vector<real_t> xtra;
    xtra.resize(2);
    xtra.at(0) = sigma_j;
    xtra.at(1) = dsph;

    std::vector<real_t> xlim_min = {1.0e-3};
    std::vector<real_t> xlim_max = {logJ_obs + 10. * sigma_j};

    std::vector<std::vector<real_t>> spect;

    IDpowell_result = IDpowell(1, &Propagation_param_t::max_likelihood_aux_func, xtra, spect, xlim_min, xlim_max,
                               &lhmax, logJ_opt, 1.0e-3, "max");

    optimal_J_factors.push_back(pow(10., logJ_opt.at(0)));
    return lhmax;
  }

  real_t Propagation_param_t::likelihood_all_dwarfs(int& IDpowell_result)
  /* Sum the maximum log-likelihood of every dSph to this sum is substracted the log-likelihood in the case of no dark
  matter This delta log-likelihood is calculated for a "conservative", "nominal" or "inclusive" case, depending if
  sample_option=-1, 0 or 1 (default is 0)*/
  {
    IDpowell_result = 1;
    int IDpowell_fail = 0; // store the number failures of IDpowell
    std::vector<real_t> l = {0, 0, 0};
    auto logJ_factors = get_logJ_factors();

    optimal_J_factors.clear();
    real_t tmp;
    for (int i = 0; i < 45; i++)
    {
      logL_noDM.at(i) = ind_param.likelihood_one_dsph(i, 0.);
      real_t sigma_j = logJ_factors.at(i).at(1);
      if (sigma_j == -1)
      {
        switch (sample_option)
        {
        case (sample_options::INCLUSIVE):
          logJ_factors.at(i).at(1) = 0.4;
          break;
        case (sample_options::CONSERVATIVE):
          logJ_factors.at(i).at(1) = 0.8;
          break;
        default:
          logJ_factors.at(i).at(1) = 0.6;
        }
      }

      tmp = max_likelihood_one_dsph(i, logJ_factors, IDpowell_result) - logL_noDM.at(i);
      if (IDpowell_result == 0)
        IDpowell_fail++;

      auto sample = (int)logJ_factors.at(i).at(2);
      switch (sample)
      {
      case (sample_options::INCLUSIVE):
        l.at(2) += tmp; // inclusive
        break;
      case (sample_options::CONSERVATIVE):
        l.at(0) += tmp;
        l.at(1) += tmp;
        l.at(2) += tmp; // conservative
        break;
      default:
        l.at(1) += tmp;
        l.at(2) += tmp; // nominal
      }
    }


    for (int i = 1; i <= 2; i++)
    {
      if (l.at(i) < l.at(0))
      {
        tmp = l.at(0);
        l.at(0) = l.at(i);
        l.at(i) = tmp;
      }
    }
    if (l.at(2) < l.at(1))
    {
      tmp = l.at(1);
      l.at(1) = l.at(2);
      l.at(2) = tmp;
    }

    real_t res;

    std::cout << "Sample option is " << sample_option << std::endl;
    switch (sample_option)
    {
    case (sample_options::INCLUSIVE):
      res = l.at(0);
      break;
    case (sample_options::CONSERVATIVE):
      res = l.at(2);
      break;
    default:
      res = l.at(1);
    }

    if (IDpowell_fail != 0)
    {
      IDpowell_result = 0;
      std::cerr << "WARNING: IDpowell failed for " << IDpowell_fail << " dSph's!!!" << std::endl;
    }
    return res;
  }



  void Propagation_param_t::search_zeroes_J0(const real_t tol)
  {
    int n_zero;

    real_t z, dz, zero;



    n_zero = 1;
    z = 2.0;
    dz = 0.01;
    for (long int iz = 1; iz <= 200000 && n_zero <= N_BESSEL; iz++)
    {
      z += dz;
      if ((besselj0(z) * besselj0(z + dz)) <= 0)
      {
        /* The zero of J0 Bessel function is between z and z+dz */
        zero = find_zero_secant(z, z + dz, besselj0, tol);
        alpha_i.at(n_zero - 1) = zero;
        n_zero += 1;
      }
    }
  }

  real_t Propagation_param_t::f_pulsar_distribution(real_t u)
  /*
    This function returns the distribution of primary cosmic-ray sources in the galactic plane.
     This radial profile follows the pulsar distribution. */
  {
    real_t r, r0, rs, fr1, fr2, z0, a, b, resultat;

    r = u * R_GAL;


    r0 = 0.55; //[kpc]
    rs = 8.50; //[kpc]
    a = 1.64;
    b = 4.01;
    z0 = E_DISC; //[kpc]
    fr1 = (r + r0) / (rs + r0);
    fr2 = (r - rs) / (rs + r0);
    resultat = pow(fr1, a) * exp((-1.0) * b * fr2);
    return resultat;
  }


  void Propagation_param_t::production()
  /* computes the flux of primary cosmic rays from the distribution of primary sources f_pulsar_distribution */
  {

    real_t u, du, sum, coefficient, weight_SIMPSON;



    sum = 0.0;
    du = 1. / (real_t)NINT_PRODUCTION;
    u = 0.0;
    for (int i_int = 0; i_int <= NINT_PRODUCTION; i_int++)
    {
      if (i_int == 0 || i_int == NINT_PRODUCTION)
      {
        weight_SIMPSON = 1. / 3.;
      }
      else
      {
        weight_SIMPSON = (1. + (real_t)(i_int % 2)) * 2. / 3.;
      }

      sum += weight_SIMPSON * f_pulsar_distribution(u) * u * du;
      u += du;
    }
    coefficient = (1.0) / (M_PI * pow(CM_TO_KPC * R_GAL, 2)) / sum; // [cm^{-2}]

    /* Calculation of integrals q_i of function f_pulsar_distribution Bessel's transforms. */
    for (int i = 0; i < N_BESSEL; i++)
    {
      sum = 0.0;
      du = 1. / (real_t)NINT_PRODUCTION;
      u = 0.0;
      for (int i_int = 0; i_int <= NINT_PRODUCTION; i_int++)
      {
        if (i_int == 0 || i_int == NINT_PRODUCTION)
        {
          weight_SIMPSON = 1. / 3.;
        }
        else
        {
          weight_SIMPSON = (1. + (real_t)(i_int % 2)) * 2. / 3.;
        }

        sum += weight_SIMPSON * f_pulsar_distribution(u) * besselj0(alpha_i.at(i) * u) * u * du;
        u += du;
      }
      q_i.at(i) = coefficient * sum / pow(besselj1(alpha_i.at(i)), 2); // [cm^{-2}]
    }
  }

  //----------------------------------------PROTON AND HELIUM FLUXES--------------------------------------------//

  real_t Propagation_param_t::flux_proton_EXP(real_t E_proton)
  /* Interstellar proton flux. It is a differential flux with respect to the
     energy, expressed in [cm^{-2} s^{-1} sr^{-1} GeV^{-1}].

     Several options are possible.
     The GAISSER and SCHAEFER proton spectrum is
     resultat = NORMALISATION_FLUX_PROTONS * pow(impulsion,-DELTA) /
     sqrt(1.+pow(P0/impulsion,2));

     The BEREZINSKY proton spectrum is
     resultat = 1.8 * (impulsion/E) * pow(E,-2.75);
     with minimal spectral index of
     resultat = 1.8 * (impulsion/E) * pow(E,-2.70);
     and maximal spectral index of
     resultat = 1.8 * (impulsion/E) * pow(E,-2.80);

     A possible value for the fit is
     resultat = N * beta_lorentz^{-1} * impulsion^{-gamma}
     resultat = N * (E/impulsion) * pow(impulsion,-gamma)
     with
     I_A)  gamma = 2.78        N = 1.8250
     I_B)  gamma = 2.88        N = 1.3700
     I_C)  gamma = 2.69        N = 2.2800

     Another parametrization is
     resultat = N * beta_lorentz * E^{-gamma}
     resultat = N * (impulsion/E) * pow(E,-gamma)
     with
     II_A)  gamma = 2.76        N = 1.5950     MEDIAN FLUX OF THE PAPER
     II_B)  gamma = 2.85        N = 1.2300
     II_C)  gamma = 2.67        N = 1.9600

     Recently, FIORENZA, NICOLAO and SANDRO found less conservative
     fits to the proton spectra.
     resultat = A * beta_lorentz * E^{-gamma}
     resultat = A * (impulsion/E) * pow(E,-gamma)
     with
     III_A)     A = 1.5300     gamma = 2.67
     III_B)     A = 1.6600     gamma = 2.85
     III_C)     A = 1.2300     gamma = 2.61     MINIMAL FLUX OF THE PAPER
     III_D)     A = 1.9600     gamma = 2.89     MAXIMAL FLUX OF THE PAPER

  **********************************************************************************************

    FEBRUARY 2000
    ANALYSIS BY FIORENZA DONATO
    Fitting both helium and proton data -- assuming the same spectral index --
    for the experiments Imax, Caprice, Mass, Leap and BESS, and AMS proton:

    N = 1.3215   gamma = 2.74   median  flux
    N = 1.2750   gamma = 2.75   minimal flux
    N = 1.3680   gamma = 2.73   maximal flux

    where we have chosen a parameterization N * rigidity^{- gamma}
    Please NOTE that the correct variable is now the RIGIDITY.

    MARCH 2005
    FIORENZA HAS SLIGHTLY MODIFIED THE PROTON FLUX WITH

    N = 1.3249   gamma = 2.72   modified median  flux

    where we have chosen the NEW parameterization N * T^{- gamma}
    Please NOTE that the correct variable is now the KINETIC ENERGY. */
  {
    real_t T, impulsion, resultat;

    T = E_proton - PROTON_MASS;
    if (T <= 0.0)
    {
      return (0.0);
    }
    else
    {
      impulsion = sqrt(pow(E_proton, 2) - pow(PROTON_MASS, 2));

      /* AMS and CREAM parametrization with phi_fisk=724MV (2015) */

      real_t C = 2.716e+04, alpha = -0.5115, gamma = -2.885, InvRb = 0.002357, DeltaGamma = 0.242, s = 0.1556;
      real_t R = sqrt(T * T + 2. * PROTON_MASS * T);
      real_t Beta = R / (T + PROTON_MASS);
      real_t dR_On_dT = (T + PROTON_MASS) / R;
      resultat = C * 1e-4 * dR_On_dT * Beta * (1 - exp(R * alpha)) * pow(R, gamma) *
                 pow(1 + pow((R * InvRb), DeltaGamma / s), s);

      return resultat;
    }
  }

  real_t Propagation_param_t::flux_helium_EXP(real_t E_nucleon)
  /* Intersellar flux of helium nuclei. This is a differential flux with respect to
     the energy per nucleon, expressed in [cm^{-2} s^{-1} sr^{-1} (GeV/nucleon)^{-1}].

     FEBRUARY 2000
     ANALYSIS BY FIORENZA DONATO
     Fitting both helium and proton data -- assuming the same spectral index --
     for the experiments Imax, Caprice, Mass, Leap and BESS, and AMS proton:

     N = 0.0645   gamma = 2.70   median  flux

     where we have chosen a parameterization N * T_nucleon^{- gamma}
     Please NOTE that the correct variable is now the KINETIC ENERGY PER NUCLEON.

     MARCH 2005
     FIORENZA HAS SLIGHTLY MODIFIED THE HELIUM FLUX WITH

     N = 0.0721   gamma = 2.74   modified median  flux */
  {
    real_t T_nucleon, P_nucleon, resultat;

    T_nucleon = E_nucleon - PROTON_MASS;
    if (T_nucleon <= 0.0)
    {
      return 0.;
    }
    else
    {
      P_nucleon = sqrt(pow(E_nucleon, 2) - pow(PROTON_MASS, 2));

      /* AMS and CREAM parametrization with phi_fisk=724MV (2015) */
      real_t C = 3564, gamma = -2.765, InvRb = 0.001842, DeltaGamma = 0.2129, s = 0.04654;
      real_t Ttot = 4. * T_nucleon;
      real_t Mtot = 4. * PROTON_MASS;
      real_t R = sqrt(Ttot * (Ttot + 2. * Mtot)) / 2.;
      real_t Beta = sqrt(Ttot * (Ttot + 2. * Mtot)) / (Ttot + Mtot);
      real_t dR_On_dT = (Ttot + Mtot) / (R);
      resultat = C * dR_On_dT * 1e-4 * Beta * pow(R, gamma) * pow(1 + pow(R * InvRb, DeltaGamma / s), s);
      return resultat;
    }
  }

  //----------------------------------------PROTON AND HELIUM BESSEL
  // COEFFICIENTS--------------------------------------------//
  void Propagation_param_t::calcul_method_BESSEL_i(real_t E_nucleon, std::vector<real_t>& bessel_coef_proton_i,
                                                   std::vector<real_t>& bessel_coef_helium_i)
  /* Calculates the flux of nuclei from its Bessel coefficients for a given energy per nucleon E_nucleon */
  {
    real_t momentum_nucleon, v_nucleon, K_proton, K_helium;
    real_t production_E_proton, Si_P, Ai_P;
    real_t production_E_helium, Si_He, Ai_He;
    real_t e_diffus = equation_parameters.at(2);
    real_t v_conv = equation_parameters.at(3);

    bessel_coef_proton_i.resize(N_BESSEL);
    bessel_coef_helium_i.resize(N_BESSEL);


    momentum_nucleon = sqrt(pow(E_nucleon, 2) - pow(PROTON_MASS, 2));
    v_nucleon = CELERITY_LIGHT * momentum_nucleon / E_nucleon;
    K_proton = K_space_diffusion(E_nucleon, PROTON_MASS, 1.0);
    K_helium = K_space_diffusion((4.0 * E_nucleon), (4.0 * PROTON_MASS), 2.0);


    for (int i = 0; i < N_BESSEL; i++)
    {
      Si_P = sqrt(pow(2.0 * alpha_i.at(i) / R_GAL, 2) + pow(v_conv * CM_TO_KPC / K_helium, 2));

      /* Ai is in [cm s^{-1}]. */
      Ai_P = v_conv;
      Ai_P += 2.0 * E_DISC * CM_TO_KPC * (pow(1.0, 2.2 / 3.) * sigma_total_pH(E_nucleon) * v_nucleon * DENSITY_H_DISC);
      Ai_P += K_proton * Si_P / CM_TO_KPC / tanh(Si_P * e_diffus / 2.);


      bessel_coef_proton_i.at(i) = q_i.at(i) / Ai_P; /* expressed in [s^{+1} cm^{-3}]. */

      /* Si is in [kpc^{-1}]. */
      Si_He = sqrt(pow(2.0 * alpha_i.at(i) / R_GAL, 2) + pow(v_conv * CM_TO_KPC / K_helium, 2));

      /* Ai is in [cm s^{-1}]. */
      Ai_He = v_conv;
      Ai_He +=
          2.0 * E_DISC * CM_TO_KPC * (pow(4.0, 2.2 / 3.) * sigma_total_pH(E_nucleon) * v_nucleon * DENSITY_HE_DISC);
      Ai_He += K_helium * Si_He / CM_TO_KPC / tanh(Si_He * e_diffus / 2.);

      bessel_coef_helium_i.at(i) = q_i.at(i) / Ai_He; /* expressed in [s^{+1} cm^{-3}]. */
    }
    production_E_proton =
        flux_proton_EXP(E_nucleon) /
        GENERIC_FLUX(r_earth, 0., E_nucleon, PROTON_MASS, 1., bessel_coef_proton_i); /* [protons s^{-1} GeV^{-1}] */
    production_E_helium = flux_helium_EXP(E_nucleon) /
                          GENERIC_FLUX(r_earth, 0., (4.0 * E_nucleon), (4.0 * PROTON_MASS), 2.0,
                                       bessel_coef_helium_i); /* expressed in  [helions s^{-1} (GeV/nucleon)^{-1}] */
    for (int i = 0; i < N_BESSEL; i++)
    {
      bessel_coef_proton_i.at(i) *= production_E_proton; /* expressed now in [protons cm^{-3} GeV^{-1}]. */
      bessel_coef_helium_i.at(i) *= production_E_helium; /* Expressed now in [helions cm^{-3} (GeV/nucleon)^{-1}]. */
    }
  }

  void Propagation_param_t::calculation_BESSEL_Ep_i()
  /* This routine calculates the Bessel coefficients for proton and helium for all values of E_nucleon, filling the
   * respective vectors*/
  {
    real_t E_nuc;


    for (int i_nuc = 0; i_nuc <= DIM_TAB_PROTON; i_nuc++)
    {
      std::vector<real_t> bessel_coef_proton_i;
      std::vector<real_t> bessel_coef_helium_i;
      E_nuc = E_PROTON_MIN * pow((E_PROTON_MAX / E_PROTON_MIN), ((real_t)i_nuc / (real_t)DIM_TAB_PROTON));
      calcul_method_BESSEL_i(E_nuc, bessel_coef_proton_i, bessel_coef_helium_i);
      for (int i = 0; i < N_BESSEL; i++)
      {
        bessel_coef_proton.at(i_nuc).at(i) = bessel_coef_proton_i.at(i);
        bessel_coef_helium.at(i_nuc).at(i) = bessel_coef_helium_i.at(i);
      }
    }
  }

  //----------------------------------------PRIMARY CONTRIBUTION TO THE ANTIPROTON
  // FLUX--------------------------------------------//
  void Propagation_param_t::calculation_BESSEL_PBAR_PRIMARY(long int n_vert, long int n_rad,
                                                            std::vector<real_t>& primary_source_term)
  /* This routine calculates Bessel coefficients bessel_pri_pbar_spec
     for antiproton kinetic energy T_pbar and Bessel coefficient i. */
  {
    real_t T_pbar, E_pbar, momentum_pbar, v_pbar, K_pbar;
    real_t Si, Abar_i;
    real_t e_diffus = equation_parameters.at(2);
    real_t v_conv = equation_parameters.at(3);

    real_t x_vert, z_vert, dx_vert, weight_SIMPSON_vert;
    real_t x_rad, r_rad, dx_rad, weight_SIMPSON_rad;
    std::array<real_t, 2001> q_pbar_primary_i_z;

    /* n_vert <= 1000 but never greater! */
    if (n_vert > 1000)
    {
      printf(" WARNING! n_vert > 1000 !");
      exit(0);
    }

    /* bessel_pri_pbar_spec is set to zero. */
    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      for (int i = 0; i < N_BESSEL; i++)
      {
        bessel_pri_pbar_spec.at(i_pbar).at(i) = 0.0;
      }
    }

    /* bessel_pri_pbar_spec is filled. */
    for (int i = 0; i < N_BESSEL; i++)
    {
      /* The radial integral q_pbar_primary_i_z is first calculated.
         \beq
         q_pbar_primary_i_z = \frac{2}{J_{1}^{2}(\alpha_{i})} \times
         \int_{0}^{1} \, x \, dx \, J_{0} \left( \alpha_{i} x \right) \,
         \left\{ \frac{\rho_{\chi}}{\rho_0} \right\}^{2} \;\; .
         \eeq

         This integral depends on i and z. Results corresponding to the ith Bessel order are saved
         in array q_pbar_primary_i_z which depends only on z.

         Loop over z followed by a loop over r */
      dx_vert = 1. / (real_t)(2 * n_vert);
      x_vert = 0.;
      for (int i_vert = 0; i_vert <= (2 * n_vert); i_vert++)
      {
        z_vert = e_diffus * x_vert; /* z_vert is the vertical coordinate expressed in [kpc]. */
        q_pbar_primary_i_z.at(i_vert) = 0.0;
        dx_rad = 1. / (real_t)(2 * n_rad);
        x_rad = 0.0;
        for (int i_rad = 0; i_rad <= (2 * n_rad); i_rad++)
        {
          r_rad = R_GAL * x_rad; /* r_rad is the radial coordinate expressed in [kpc]. */
          if (i_rad == 0 || i_rad == (2 * n_rad))
          {
            weight_SIMPSON_rad = 1. / 3.;
          }
          else
          {
            weight_SIMPSON_rad = (1. + (real_t)(i_rad % 2)) * 2. / 3.;
          }

          q_pbar_primary_i_z.at(i_vert) += (x_rad * dx_rad * weight_SIMPSON_rad) * besselj0(alpha_i.at(i) * x_rad) *
                                           pow(halo_profile(r_rad, z_vert), 2.0); //[NO UNIT].


          x_rad += dx_rad;
        }
        q_pbar_primary_i_z.at(i_vert) *= 2. / pow(besselj1(alpha_i.at(i)), 2.0); /* [NO UNIT]. */

        x_vert += dx_vert;
      }

      /* Loop over antiproton energy containing a loop over the vertical coordinate,
         allowing us to calculate calculer bessel_pri_pbar_spec[i_pbar].at(i). */
      for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
      {
        T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
        E_pbar = T_pbar + PROTON_MASS;
        momentum_pbar = sqrt(pow(E_pbar, 2) - pow(PROTON_MASS, 2));
        v_pbar = CELERITY_LIGHT * momentum_pbar / E_pbar;
        K_pbar = K_space_diffusion(E_pbar, PROTON_MASS, 1.0);

        Si = sqrt(pow(2.0 * alpha_i.at(i) / R_GAL, 2) +
                  pow(v_conv * CM_TO_KPC / K_pbar, 2)); /* Si is expressed in [kpc^{-1}]. */

        Abar_i = v_conv; /* Abar_i is expressed in [cm s^{-1}]. */
        Abar_i += 2.0 * E_DISC * CM_TO_KPC *
                  ((sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)) *
                   v_pbar * (DENSITY_H_DISC + pow(4., (2. / 3.)) * 1.0 * DENSITY_HE_DISC));
        Abar_i += K_pbar * Si / CM_TO_KPC / tanh(Si * e_diffus / 2.);

        table_abar.at(i_pbar).at(i) = Abar_i;


        /* Integration over the vertical variable x_vert and over the radial variable x_rad. */
        dx_vert = 1. / (real_t)(2 * n_vert);
        x_vert = 0.;
        bessel_pri_pbar_spec.at(i_pbar).at(i) = 0.0;
        for (int i_vert = 0; i_vert <= (2 * n_vert); i_vert++)
        {
          z_vert = e_diffus * x_vert;

          if (i_vert == 0 || i_vert == (2 * n_vert))
          {
            weight_SIMPSON_vert = 1. / 3.;
          }
          else
          {
            weight_SIMPSON_vert = (1. + (real_t)(i_vert % 2)) * 2. / 3.;
          }

          bessel_pri_pbar_spec.at(i_pbar).at(i) +=
              (dx_vert * weight_SIMPSON_vert) * q_pbar_primary_i_z.at(i_vert) * primary_source_term.at(i_pbar) *
              exp(-v_conv * z_vert * CM_TO_KPC / (2. * K_pbar)) * sinh((Si / 2.) * (e_diffus - z_vert)) /
              sinh((Si / 2.) * e_diffus); //[antiprotons cm^{-3} s^{-1} GeV^{-1}].

          x_vert += dx_vert;
        }
        bessel_pri_pbar_spec.at(i_pbar).at(i) *=
            2. * e_diffus * CM_TO_KPC / Abar_i; /* expressed in [antiprotons cm^{-3} GeV^{-1}]. */
      }
    }
    return;
  }


  //----------------------------------------CALCULATING PRIMARY SOURCE
  // TERM--------------------------------------------//
  std::vector<real_t>
  Propagation_param_t::DNPBAR_ON_DTPBAR_read_file(real_t mass_chi,
                                                  std::vector<std::vector<real_t>>& primary_source_terms_inf,
                                                  std::vector<std::vector<real_t>>& primary_source_terms_sup)
  /* read tabulated primary par spectra at production for DM masses m_inf and m_sup such that m_inf<=mass_chi<=m_sup
   * aand save the spectra in pt_Primary_Source_Term. */
  {
    std::string filename_inf;
    std::string filename_sup;
    real_t mass_min, mass_max, mass_inf, mass_sup;

    std::vector<real_t> mass_list = get_DM_masses();

    const int N_MASS = mass_list.size();



    mass_min = mass_list.at(0);
    mass_max = mass_list.at(N_MASS - 1);

    auto read_file_new = [&](const std::string& filename, std::vector<std::vector<real_t>>& primary_source_terms)
    {
      // Replicate ind_param.read_file, but adapted to sum contributions from different processes that belong to the
      // same channel
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
        if (line.empty() || line.at(0) == '#')
          continue; // Skip comment lines
        real_t val;
        int count = 0;
        std::vector<real_t> values;
        values.clear();
        std::stringstream ss(line);
        if (!std::isdigit(line.at(0)) && line.at(0) != '-' && line.at(0) != '+')
        {
          std::string dummy;
          ss >> dummy; // Skip non-numeric starting values
        }
        while (ss >> val)
        {
          values.push_back(val);
        }

        // Assuming 14 channels
        std::vector<real_t> channels;
        channels.resize(14);
        channels.at(0) = (values.at(1) + values.at(2)) / 2.;
        channels.at(1) = (values.at(3) + values.at(4)) / 2.;
        channels.at(2) = (values.at(5) + values.at(6)) / 2.;
        channels.at(3) = values.at(7);
        channels.at(4) = values.at(8);
        channels.at(5) = values.at(9);
        channels.at(6) = values.at(10);
        channels.at(7) = (2. * values.at(12) + values.at(11)) / 3.;
        channels.at(8) = (2. * values.at(14) + values.at(13)) / 3.;
        channels.at(9) = values.at(15);
        channels.at(10) = values.at(16);
        channels.at(11) = values.at(17);
        channels.at(12) = values.at(18);
        channels.at(13) = values.at(19);
        primary_source_terms.emplace_back(channels);
      }
      std::vector<real_t> dummy_zeros(14, 0);
      primary_source_terms.emplace_back(dummy_zeros);
    };

    if ((mass_chi < mass_min) || (mass_chi > mass_max))
    {
      printf("\n ERROR ! \n Function : 'DNPBAR_ON_DTPBAR_read_file'  \n mass_chi is out of range ! \n");
      exit(0);
    }
    if ((mass_chi >= mass_min) && (mass_chi <= mass_max))
    {
      int ind_mass = 0;
      while (mass_chi < mass_list.at(ind_mass) || mass_chi >= mass_list.at(ind_mass + 1))
      {
        ind_mass++;
      }
      // prevents the singularities due to the opening of W, Z, or t channels
      if (mass_chi < 90 && mass_chi > 80)
        ind_mass++;
      if (mass_chi < 100 && mass_chi > 90)
        ind_mass++;
      if (mass_chi < 200 && mass_chi > 173)
        ind_mass++;
      for (const auto& file : std::filesystem::directory_iterator(PRIMARY_PBAR_SOURCES))
      {
        std::string filename = file.path().string();
        if (std::stod(file.path().stem().string()) == mass_list.at(ind_mass))
        {
          filename_inf = file.path().string();
          read_file_new(filename_inf, primary_source_terms_inf);
          mass_inf = mass_list.at(ind_mass);
        }
        else if (std::stod(file.path().stem().string()) == mass_list.at(ind_mass + 1))
        {
          filename_sup = file.path().string();
          read_file_new(filename_sup, primary_source_terms_sup);
          mass_sup = mass_list.at(ind_mass + 1);
        }
      }
    }

    return {mass_inf, mass_sup};
  }



  int Propagation_param_t::dNpbar_on_dEpbar_calculation(std::vector<real_t>& dNpbar_on_dEpbar_vec)
  /* Calculates the primary pbar spectrum for a DM particle whose mass and annihilation channels ared defined in spec
   * and save it in pt_Primary_Source_Term */
  {
    int i_scan_pbar;
    real_t xi_x_pbar, T_pbar, x_pbar;
    real_t dNpbar_on_dEpbar_xi_m, dNpbar_on_dEpbar_xi_plus_un_m;
    int DM_candidate = input.getLightestBSMpart();
    real_t mass_chi = input.masses_vector.at(DM_candidate);
    auto energy_table = get_energy_table();

    if (energy_table.size() == 0)
    {
      std::cerr << "Energy table is empty, filling the spectrum with the processes and sigma_v_process\n";
      ind_param.fill_spectrum();
      energy_table = get_energy_table();
    }

    auto sigma_v_table = get_sigma_v_table();
    int n1 = sigma_v_table.at(0).size();
    int number_channels = sigma_v_table.size();
    dNpbar_on_dEpbar_vec.resize(DIM_TAB_PBAR + 1);

    std::vector<std::vector<real_t>> temp_source_terms_inf, temp_source_terms_sup;



    for (int i_mass = 0; i_mass < n1; i_mass++)
    {
      real_t mass_chi_new = energy_table.at(i_mass);
      auto inf_sup_masses = DNPBAR_ON_DTPBAR_read_file(mass_chi_new, temp_source_terms_inf, temp_source_terms_sup);
      real_t mass_inf = inf_sup_masses.at(0);
      real_t mass_sup = inf_sup_masses.at(1);
      for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
      {
        T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
        x_pbar = T_pbar / mass_chi_new;

        if (x_pbar < x_pbar_scan_min)
        {
          x_pbar = x_pbar_scan_min;
        }
        if ((x_pbar_scan_min <= x_pbar) && (x_pbar < x_pbar_scan_max))
        {
          xi_x_pbar = ((real_t)N_x_pbar_scan) * log(x_pbar / ((real_t)x_pbar_scan_min)) /
                      log(((real_t)x_pbar_scan_max) / ((real_t)x_pbar_scan_min));
          if (xi_x_pbar > ((real_t)N_x_pbar_scan))
          {
            printf(
                "\n ERROR ! \n Function : 'dNpbar_on_dEpbar_primary_calculation'  \n xi_x_pbar >= N_x_pbar_scan ! \n");
            return 0;
          }
          else if (xi_x_pbar == ((real_t)N_x_pbar_scan))
          {
            for (int channel = 0; channel < number_channels; channel++)
            {
              dNpbar_on_dEpbar_vec.at(i_pbar) +=
                  sigma_v_table.at(channel).at(i_mass) * temp_source_terms_inf.at(N_x_pbar_scan).at(channel);
              dNpbar_on_dEpbar_vec.at(i_pbar) += sigma_v_table.at(channel).at(i_mass) *
                                                 (temp_source_terms_sup.at(N_x_pbar_scan).at(channel) -
                                                  temp_source_terms_inf.at(N_x_pbar_scan).at(channel)) *
                                                 (mass_chi_new - mass_inf) / (mass_sup - mass_inf);
            }
          }
          else
          {
            i_scan_pbar = xi_x_pbar;
            for (int channel = 0; channel < number_channels; channel++)
            {
              dNpbar_on_dEpbar_xi_m = temp_source_terms_inf.at(i_scan_pbar).at(channel);
              dNpbar_on_dEpbar_xi_m += (temp_source_terms_sup.at(i_scan_pbar).at(channel) -
                                        temp_source_terms_inf.at(i_scan_pbar).at(channel)) *
                                       (mass_chi_new - mass_inf) / (mass_sup - mass_inf);


              dNpbar_on_dEpbar_xi_plus_un_m = temp_source_terms_inf.at(i_scan_pbar + 1).at(channel);
              dNpbar_on_dEpbar_xi_plus_un_m += (temp_source_terms_sup.at(i_scan_pbar + 1).at(channel) -
                                                temp_source_terms_inf.at(i_scan_pbar + 1).at(channel)) *
                                               (mass_chi_new - mass_inf) / (mass_sup - mass_inf);


              dNpbar_on_dEpbar_vec.at(i_pbar) += sigma_v_table.at(channel).at(i_mass) * dNpbar_on_dEpbar_xi_m;
              dNpbar_on_dEpbar_vec.at(i_pbar) += sigma_v_table.at(channel).at(i_mass) *
                                                 (dNpbar_on_dEpbar_xi_plus_un_m - dNpbar_on_dEpbar_xi_m) *
                                                 (xi_x_pbar - (real_t)i_scan_pbar);
            }
          }
        }
        if (x_pbar_scan_max <= x_pbar)
        {
          dNpbar_on_dEpbar_vec.at(i_pbar) += 0.0;
        }
      }
    }
    return (1.0);
  }


  void Propagation_param_t::primary_source_calculation(real_t mass_chi, std::vector<real_t>& dNpbar_on_dEpbar_vec,
                                                       std::vector<real_t>& primary_source_term)
  /* modulates the primary pbar spectrum at production with the total annihilation cross section and the DM number
   * density */
  {
    std::cerr << "Calculating primary source term...\n";

    primary_source_term.resize(DIM_TAB_PBAR + 1);

    dNpbar_on_dEpbar_calculation(dNpbar_on_dEpbar_vec);

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      if (dNpbar_on_dEpbar_vec.at(i_pbar) > 1.0e-100)
        primary_source_term.at(i_pbar) = 0.5 * sigma_v_annihilation * pow((RHO_CHI_0 / mass_chi), 2.0) *
                                         dNpbar_on_dEpbar_vec.at(i_pbar); // [#pbar cm{-3} s{-1} GeV{-1}]
    }
  }

  void Propagation_param_t::primary_spectra_BCGS_2014(std::vector<real_t>& dNpbar_on_dEpbar_vec,
                                                      std::vector<real_t>& primary_source_term,
                                                      std::vector<std::vector<real_t>>& PBAR_SPECTRUM)
  /* Calculates the primary pbar spectrum at Earth position for a DM particle defined in structure spec and save the
   * result in array PBAR_SPECTRUM */
  {
    /*	char file_name[1024];*/

    int DM_candidate = input.getLightestBSMpart();
    real_t mass_chi = input.masses_vector.at(DM_candidate);

    real_t T_pbar_IS, E_pbar_IS, flux_antiproton_IS;

    PBAR_SPECTRUM.resize(DIM_TAB_PBAR + 1);

    for (int i = 0; i <= DIM_TAB_PBAR; i++)
    {
      PBAR_SPECTRUM.at(i).resize(2);
    }


    primary_source_calculation(mass_chi, dNpbar_on_dEpbar_vec, primary_source_term);


    /* We define the parameters determined by FIORENZA, DAVID et RICHARD -- hereafter called FDR. */

    /* pt_Pbar->BESSEL_PBAR_SEC_Epbar_i and pt_Pbar->BESSEL_PBAR_TER_Epbar_i are set to zero. */
    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      for (int i = 0; i < N_BESSEL; i++)
      {
        bessel_pri_pbar_spec.at(i_pbar).at(i) = 0.0;
        bessel_sec_pbar_spec.at(i_pbar).at(i) = 0.0;
        bessel_ter_pbar_spec.at(i_pbar).at(i) = 0.0;
      }
    }

    /* CALCULATION OF THE PRIMARY CONTRIBUTION FROM NEUTRALINO ANNIHILATION */
    calculation_BESSEL_PBAR_PRIMARY(100, 500, primary_source_term);

    /* CALCULATION OF SECONDARY CONTRIBUTION FROM THE INTERSTELLAR GAS SPALLATION BY COSMIC-RAY PROTONS AND HELIONS */
    calculation_BESSEL_PBAR_SUM();

    /* CALCULATION OF FINAL ANTIPROTON SPECTRUM */
    for (int i_iteration = 1; i_iteration <= 5; i_iteration++)
    {
      calculation_BESSEL_PBAR_TERTIARY();
      calculation_BESSEL_PBAR_SUM();
    }

    for (int i_iteration = 1; i_iteration <= 5; i_iteration++)
    {
      calculation_BESSEL_PBAR_TERTIARY();
      calculation_BESSEL_PBAR_direct_inversion();
    }

    /* We compute now the antiproton spectrum */
    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      std::vector<real_t> BESSEL_PBARi;
      T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
      PBAR_SPECTRUM.at(i_pbar).at(0) = log10(T_pbar_IS);
      E_pbar_IS = T_pbar_IS + PROTON_MASS;
      for (int i = 0; i < N_BESSEL; i++)
      {
        BESSEL_PBARi.push_back(bessel_tot_pbar_spec.at(i_pbar).at(i));
      }
      flux_antiproton_IS = GENERIC_FLUX_04(r_earth, 0., E_pbar_IS, PROTON_MASS, 1., BESSEL_PBARi);
      PBAR_SPECTRUM.at(i_pbar).at(1) = flux_antiproton_IS * 1.0e+4 / 3.e-26;
    }
  }

  //----------------------------------------SECONDARY CONTRIBUTION--------------------------------------------//
  void Propagation_param_t::calculation_BESSEL_PBAR_SECONDARY()
  /* This routine computes the elements of the array  bessel_sec_pbar_spec
     as a function of the antiproton kinetic energy T_pbar and BESSEL coefficients i. */
  {

    real_t T_pbar, E_pbar, momentum_pbar, v_pbar, K_pbar;
    real_t Si, Abar_i;
    real_t dlog_E_proton, E_proton;
    std::vector<real_t> impulsion_proton, weight_SIMSPON;
    impulsion_proton.resize(DIM_TAB_PROTON + 1);
    weight_SIMSPON.resize(DIM_TAB_PROTON + 1);


    real_t e_diffus = equation_parameters.at(2);
    real_t v_conv = equation_parameters.at(3);

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* bessel_sec_pbar_spec and table_abar are set to zero. */
    {
      for (int i = 0; i < N_BESSEL; i++)
      {
        bessel_sec_pbar_spec.at(i_pbar).at(i) = 0.0;
      }
    }

    dlog_E_proton = pow((E_PROTON_MAX / E_PROTON_MIN), (1. / (real_t)DIM_TAB_PROTON)) - 1.;
    for (int i_proton = 0; i_proton <= DIM_TAB_PROTON; i_proton++)
    /* Loop on the protons in order to calulate E_proton, v_proton and weight_SIMSPON. */
    {
      E_proton = E_PROTON_MIN * pow((E_PROTON_MAX / E_PROTON_MIN), ((real_t)i_proton / (real_t)DIM_TAB_PROTON));
      impulsion_proton.at(i_proton) = sqrt(pow(E_proton, 2) - pow(PROTON_MASS, 2));

      if (i_proton == 0 || i_proton == DIM_TAB_PROTON)
      {
        weight_SIMSPON.at(i_proton) = 1. / 3.;
      }
      else
      {
        weight_SIMSPON.at(i_proton) = (1. + (real_t)(i_proton % 2)) * 2. / 3.;
      }
    }


    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* Array BESSEL_PBAR_SEC_Epbar_i is filled. */
    {
      table_abar.at(i_pbar).resize(N_BESSEL);
      T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
      E_pbar = T_pbar + PROTON_MASS;
      momentum_pbar = sqrt(pow(E_pbar, 2) - pow(PROTON_MASS, 2));
      v_pbar = CELERITY_LIGHT * momentum_pbar / E_pbar;
      K_pbar = K_space_diffusion(E_pbar, PROTON_MASS, 1.0);

      for (int i = 0; i < N_BESSEL; i++)
      {
        /* Si is in [kpc^{-1}]. */
        Si = sqrt(pow(2.0 * alpha_i.at(i) / R_GAL, 2) + pow(v_conv * CM_TO_KPC / K_pbar, 2));
        /* Abar_i in [cm s^{-1}]. */
        Abar_i = v_conv;
        Abar_i += 2.0 * E_DISC * CM_TO_KPC *
                  ((sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar)) *
                   v_pbar * (DENSITY_H_DISC + pow(4., (2. / 3.)) * 1.0 * DENSITY_HE_DISC));
        Abar_i += K_pbar * Si / CM_TO_KPC / tanh(Si * e_diffus / 2.);
        table_abar.at(i_pbar).at(i) = Abar_i;


        for (int i_proton = 0; i_proton <= DIM_TAB_PROTON; i_proton++)
        {
          /* CONTRIBUTION H ON H */
          bessel_sec_pbar_spec.at(i_pbar).at(i) +=
              h_on_h_xsection.at(i_proton).at(i_pbar) * bessel_coef_proton.at(i_proton).at(i) * DENSITY_H_DISC *
              impulsion_proton.at(i_proton) * CELERITY_LIGHT * weight_SIMSPON.at(i_proton) *
              dlog_E_proton; // [antiprotons GeV^{-1} s^{-1} cm^{-3}]

          /* CONTRIBUTION H ON HE */
          bessel_sec_pbar_spec.at(i_pbar).at(i) +=
              h_on_he_xsection.at(i_proton).at(i_pbar) * bessel_coef_proton.at(i_proton).at(i) * DENSITY_HE_DISC *
              impulsion_proton.at(i_proton) * CELERITY_LIGHT * weight_SIMSPON.at(i_proton) *
              dlog_E_proton; // [antiprotons GeV^{-1} s^{-1} cm^{-3}]

          /* CONTRIBUTION HE ON H */
          bessel_sec_pbar_spec.at(i_pbar).at(i) +=
              he_on_h_xsection.at(i_proton).at(i_pbar) * bessel_coef_helium.at(i_proton).at(i) * DENSITY_H_DISC *
              impulsion_proton.at(i_proton) * CELERITY_LIGHT * weight_SIMSPON.at(i_proton) *
              dlog_E_proton; // [antiprotons GeV^{-1} s^{-1} cm^{-3}]

          /* CONTRIBUTION HE ON HE */
          bessel_sec_pbar_spec.at(i_pbar).at(i) +=
              he_on_he_xsection.at(i_proton).at(i_pbar) * bessel_coef_helium.at(i_proton).at(i) * DENSITY_HE_DISC *
              impulsion_proton.at(i_proton) * CELERITY_LIGHT * weight_SIMSPON.at(i_proton) *
              dlog_E_proton; //[antiprotons GeV^{-1} s^{-1} cm^{-3}]
        }
        bessel_sec_pbar_spec.at(i_pbar).at(i) *= 2. * E_DISC * CM_TO_KPC / Abar_i;
        /* expressed in [antiprotons cm^{-3} GeV^{-1}]. */
      }
    }
  }

  void Propagation_param_t::preliminary_secondary_spectrum_IS_calculation()
  /* Routine which computes the secondary antiproton spectrum */
  {
    xs_uncertainties.resize(DIM_TAB_PBAR + 1);

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      xs_uncertainties.at(i_pbar).resize(3); /* Element 0 is energy,
      1 and 2 are the lower and upper bounds of the secondary spectra, respectively*/

      for (int i = 0; i < N_BESSEL; i++)
      {
        bessel_pri_pbar_spec.at(i_pbar).at(i) = 0.0;
        bessel_sec_pbar_spec.at(i_pbar).at(i) = 0.0;
        bessel_ter_pbar_spec.at(i_pbar).at(i) = 0.0;
        bessel_tot_pbar_spec.at(i_pbar).at(i) = 0.0;
      }
    }


    /* calculation of the proton and helium P_i. */
    calculation_BESSEL_Ep_i();

    /* PBAR_IS_spectrum_II_MIN is filled */
    calculation_BESSEL_PBAR_SECONDARY();
    calculation_XS_UNCERTAINTIES_ON_BESSEL_PBAR_SECONDARY(-1);
    calculation_BESSEL_PBAR_SUM();
    ELDR_effect_calculation();

    real_t T_pbar_IS, E_pbar_IS, flux_antiproton_IS;


    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      std::vector<real_t> BESSEL_PBARi_minus;
      T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
      E_pbar_IS = T_pbar_IS + PROTON_MASS;
      for (int i = 0; i < N_BESSEL; i++)
      {
        BESSEL_PBARi_minus.push_back(bessel_tot_pbar_spec.at(i_pbar).at(i));
      }
      flux_antiproton_IS = GENERIC_FLUX_04(r_earth, 0., E_pbar_IS, PROTON_MASS, 1., BESSEL_PBARi_minus);

      xs_uncertainties.at(i_pbar).at(0) = (T_pbar_IS);
      xs_uncertainties.at(i_pbar).at(1) = flux_antiproton_IS * 1.0e4; /* [#pbar cm^{-2} sr^{-1} s^{-1} GeV^{-1}] */
    }

    /* PBAR_IS_spectrum_II_MAX is filled */
    calculation_BESSEL_PBAR_SECONDARY();
    calculation_XS_UNCERTAINTIES_ON_BESSEL_PBAR_SECONDARY(1);
    calculation_BESSEL_PBAR_SUM();
    ELDR_effect_calculation();

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      std::vector<real_t> BESSEL_PBARi_plus;
      T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
      E_pbar_IS = T_pbar_IS + PROTON_MASS;
      for (int i = 0; i < N_BESSEL; i++)
      {
        BESSEL_PBARi_plus.push_back(bessel_tot_pbar_spec[i_pbar].at(i));
      }
      flux_antiproton_IS = GENERIC_FLUX_04(r_earth, 0., E_pbar_IS, PROTON_MASS, 1., BESSEL_PBARi_plus);

      xs_uncertainties.at(i_pbar).at(2) = flux_antiproton_IS * 1.0e+4; /* [#pbar cm^{-2} sr^{-1} s^{-1} GeV^{-1}] */
    }

    return;
  }

  void Propagation_param_t::calculation_XS_UNCERTAINTIES_ON_BESSEL_PBAR_SECONDARY(int option)
  /* Routine which computes the lower (option=-1) or upper (option=1) bound of the secondary pbar spectrum */
  {

    if (option == 0)
    {

      return;
    }
    else
    {


      std::vector<real_t> ENERGY, F_Q_PBAR_SEC_PLUS, F_Q_PBAR_SEC_MINUS;
      ENERGY.resize(DIM_TAB_PBAR_MFGS);
      F_Q_PBAR_SEC_PLUS.resize(DIM_TAB_PBAR_MFGS);
      F_Q_PBAR_SEC_MINUS.resize(DIM_TAB_PBAR_MFGS);
      real_t T_pbar, q_pbar, q_pbar_sec_minus, q_pbar_sec_plus;

      std::vector<std::vector<real_t>> temp_errors;
      ind_param.read_file(XSECTION_ERRORS, 4, temp_errors);

      for (int i_pbar = 0; i_pbar < DIM_TAB_PBAR_MFGS; i_pbar++)
      {

        T_pbar = temp_errors.at(i_pbar).at(0);
        q_pbar = temp_errors.at(i_pbar).at(1);
        q_pbar_sec_minus = temp_errors.at(i_pbar).at(2);
        q_pbar_sec_plus = temp_errors.at(i_pbar).at(3);


        ENERGY.at(i_pbar) = T_pbar;
        F_Q_PBAR_SEC_PLUS.at(i_pbar) = q_pbar_sec_plus / q_pbar;
        F_Q_PBAR_SEC_MINUS.at(i_pbar) = q_pbar_sec_minus / q_pbar;
      }

      for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
      {
        real_t Energy_pbar = pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR)) * T_PBAR_MIN;

        int j_mem = -1;

        if (Energy_pbar < ENERGY.at(0))
        {
          j_mem = 0;
        }
        else if (Energy_pbar > ENERGY.at(DIM_TAB_PBAR_MFGS - 1))
        {
          j_mem = DIM_TAB_PBAR_MFGS - 2;
        }
        else
        {
          int j;
          for (j = 0; j < DIM_TAB_PBAR_MFGS - 1; j++)
          {

            if (Energy_pbar >= ENERGY.at(j) && Energy_pbar <= ENERGY.at(j + 1))
            {
              j_mem = j;
            }
          }
        }
        if (j_mem == -1)
        {
          printf("PB WITH FUNCTION XS_UNCERTAINTIES \t j_mem = %d", j_mem);
        }
        for (int i = 0; i < N_BESSEL; i++)
        {
          if (option == -1)
          {
            bessel_sec_pbar_spec.at(i_pbar).at(i) *=
                ((F_Q_PBAR_SEC_MINUS.at(j_mem + 1) - F_Q_PBAR_SEC_MINUS.at(j_mem)) /
                     (ENERGY.at(j_mem + 1) - ENERGY.at(j_mem)) * (Energy_pbar - ENERGY.at(j_mem)) +
                 F_Q_PBAR_SEC_MINUS.at(j_mem));
          }
          if (option == 1)
          {
            bessel_sec_pbar_spec.at(i_pbar).at(i) *=
                ((F_Q_PBAR_SEC_PLUS.at(j_mem + 1) - F_Q_PBAR_SEC_PLUS.at(j_mem)) /
                     (ENERGY.at(j_mem + 1) - ENERGY.at(j_mem)) * (Energy_pbar - ENERGY.at(j_mem)) +
                 F_Q_PBAR_SEC_PLUS.at(j_mem));
          }
        }
      }
    }

    return;
  }

  int Propagation_param_t::background_spectAMS(const std::vector<real_t>& logE, real_t A,
                                               std::vector<std::vector<real_t>>& result)
  /* Computes secondary pbar spectrum for parameter A  at energies logE */
  {
    int N = xs_uncertainties.size();

    if (N == 0) // Calls the function to compute the cross section uncertainties if the vector has not been filled
    {
      preliminary_secondary_spectrum_IS_calculation();
      N = xs_uncertainties.size();
    }
    int j = 0;

    auto n = logE.size();

    result.resize(n);

    real_t Et;
    real_t result_down, result_up;
    for (size_t i = 0; i < n; i++)
    {
      result.at(i).resize(2);
      Et = pow(10., logE.at(i));
      while (!(xs_uncertainties.at(j).at(0) <= Et && xs_uncertainties.at(j + 1).at(0) >= Et) && j < (N - 1))
        j++;

      if (j == N - 1)
      {
        printf("Error in function background_spectAMS: energy out of range");
        printf("%e %e %e\n", xs_uncertainties.at(j).at(0), Et, xs_uncertainties.at(j + 1).at(0));
        return 0;
      }
      else
      {
        result_down = (xs_uncertainties.at(j).at(1) * xs_uncertainties.at(j + 1).at(0) -
                       xs_uncertainties.at(j).at(0) * xs_uncertainties.at(j + 1).at(1)) /
                          (xs_uncertainties.at(j + 1).at(0) - xs_uncertainties.at(j).at(0)) +
                      (xs_uncertainties.at(j + 1).at(1) - xs_uncertainties.at(j).at(1)) /
                          (xs_uncertainties.at(j + 1).at(0) - xs_uncertainties.at(j).at(0)) * Et;
        result_up = (xs_uncertainties.at(j).at(2) * xs_uncertainties.at(j + 1).at(0) -
                     xs_uncertainties.at(j).at(0) * xs_uncertainties.at(j + 1).at(2)) /
                        (xs_uncertainties.at(j + 1).at(0) - xs_uncertainties.at(j).at(0)) +
                    (xs_uncertainties.at(j + 1).at(2) - xs_uncertainties.at(j).at(2)) /
                        (xs_uncertainties.at(j + 1).at(0) - xs_uncertainties.at(j).at(0)) * Et;
        result.at(i).at(0) = logE.at(i);
        result.at(i).at(1) = result_down + A * (result_up - result_down);
      }
    }
    return 1;
  }

  //----------------------------------------TERTIARY CONTRIBUTION-----------------------------------------//

  void Propagation_param_t::calculation_BESSEL_PBAR_TERTIARY()
  /* This routine computes the elements of the array  bessel_ter_pbar_spec
     as a function of the antiproton kinetic energy T_pbar and BESSEL coefficients i. */
  {

    real_t dlog_T_pbar, T_pbar, E_pbar, momentum_pbar, v_pbar;
    static std::vector<real_t> S_inel_NOANN_times_v_pbar;
    S_inel_NOANN_times_v_pbar.resize(DIM_TAB_PBAR + 1);
    real_t Abar_i, SUM;

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* Loop on the antiprotons in order to calculate the quantity sigma_inelastic_NOANN_pbarH_TAN_and_NG * v_pbar which
     * will be used in the following loops. */
    {
      T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
      E_pbar = T_pbar + PROTON_MASS;
      momentum_pbar = sqrt(pow(E_pbar, 2) - pow(PROTON_MASS, 2));
      v_pbar = CELERITY_LIGHT * momentum_pbar / E_pbar;

      S_inel_NOANN_times_v_pbar.at(i_pbar) = sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar) * v_pbar;
    }

    for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* bessel_ter_pbar_spec is set to zero. */
    {
      for (int i = 0; i < N_BESSEL; i++)
      {
        bessel_ter_pbar_spec.at(i_pbar).at(i) = 0.0;
      }
    }

    dlog_T_pbar = pow((T_PBAR_MAX / T_PBAR_MIN), (1. / (real_t)DIM_TAB_PBAR)) - 1.;
    for (int i = 0; i < N_BESSEL; i++)
    /* bessel_ter_pbar_spec is now filled. */
    {
      SUM = 0.0;
      for (int i_pbar = DIM_TAB_PBAR; i_pbar >= 0; i_pbar--)
      {
        Abar_i = table_abar.at(i_pbar).at(i); /* Abar_i is in [cm s^{-1}]. */
        SUM += dlog_T_pbar * S_inel_NOANN_times_v_pbar.at(i_pbar) * bessel_tot_pbar_spec.at(i_pbar).at(i);
        bessel_ter_pbar_spec.at(i_pbar).at(i) =
            SUM - S_inel_NOANN_times_v_pbar.at(i_pbar) * bessel_tot_pbar_spec.at(i_pbar).at(i);
        /* expressed in [antiprotons GeV^{-1} s^{-1}]. */

        bessel_ter_pbar_spec.at(i_pbar).at(i) *=
            2. * E_DISC * CM_TO_KPC * (DENSITY_H_DISC + pow(4., (2. / 3.)) * 1.0 * DENSITY_HE_DISC) / Abar_i;
        /* expressed now in [antiprotons cm^{-3} GeV^{-1}]. */
      }
    }
  }

  //----------------------------------------SUM OF PRIMARY, SECONDARY AND TERTIARY
  // CONTRIBUTIONS-----------------------------------------//

  void Propagation_param_t::calculation_BESSEL_PBAR_SUM()
  /* This routine calculates the elements ofpt_Pbar->BESSEL_PBAR_TOT_Epbar_i
     and the sum of the primary, secondary and tertiary antiproton contributions
    as a function of the antiproton kinetic energy T_pbar and BESSEL coefficients i. */
  {

    for (long int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* pt_Pbar->BESSEL_PBAR_TOT_Epbar_i is filled. */
    {
      for (long int i = 0; i < N_BESSEL; i++)
      {
        bessel_tot_pbar_spec.at(i_pbar).at(i) = bessel_pri_pbar_spec.at(i_pbar).at(i) +
                                                bessel_sec_pbar_spec.at(i_pbar).at(i) +
                                                bessel_ter_pbar_spec.at(i_pbar).at(i);
      }
    }
  }

  //----------------------------------------SOLAR MODULATION-----------------------------------------//
  void Propagation_param_t::solar_mod(std::vector<std::vector<real_t>>& spectrum, real_t phi_f,
                                      std::vector<std::vector<real_t>>& result)
  /* computes the total spectrum after solar modulation with Fisk potential phi_f */
  {
    real_t Z = 1, A = 1;
    real_t m = 938.2720813e-03;
    real_t e = 1.;

    real_t K, K_prime, logK_prime;
    size_t size = spectrum.size();
    result.resize(size);

    for (size_t i = 0; i < size; i++)
    {
      result.at(i).resize(2);
      K = pow(10., spectrum.at(i).at(0));
      K_prime = K + e * phi_f * Z / A;
      logK_prime = log10(K_prime);

      if (fabs(logK_prime - spectrum.at(size - 1).at(0)) < 1.e-4)
        logK_prime = spectrum.at(size - 1).at(0);

      if (logK_prime >= spectrum.at(0).at(0) && logK_prime <= spectrum.at(size - 1).at(0))
      {
        result.at(i).at(0) = spectrum.at(i).at(0);
        result.at(i).at(1) = ind_param.logx_interpol(spectrum, logK_prime) * (2. * m * K + K * K) /
                             (2. * m * K_prime + K_prime * K_prime);
      }
      else
      {
        result.at(i).at(0) = spectrum.at(i).at(0);
        result.at(i).at(1) = 0;
      }
    }
  }


  //----------------------------------------COSMIC RAY FLUXES FROM BESSEL
  // COEFFICIENTS-----------------------------------------//

  real_t Propagation_param_t::GENERIC_FLUX(real_t r, real_t z, real_t energy, real_t mass, real_t Z_em,
                                           std::vector<real_t> bes_coef)
  /* GENERIC_FLUX is a generic routine which returns a cosmic-ray flux from Bessel coefficients calculated beforehand */
  {
    real_t x, az;
    real_t momentum, velocity, K;
    real_t Si, resultat;

    x = r / R_GAL;
    az = std::abs(z);
    resultat = 0.0;

    real_t e_diffus = equation_parameters.at(2);
    real_t v_conv = equation_parameters.at(3);

    if (az >= e_diffus || x >= 1.0)
    {
      return resultat;
    }
    else
    {
      momentum = sqrt(pow(energy, 2) - pow(mass, 2));
      velocity = CELERITY_LIGHT * momentum / energy;
      K = K_space_diffusion(energy, mass, Z_em);

      for (long int i = 0; i < N_BESSEL; i++)
      {
        /* Si is expressed in [kpc^{-1}]. */
        Si = sqrt(pow(2.0 * alpha_i.at(i) / R_GAL, 2) + pow(v_conv * CM_TO_KPC / K, 2));
        resultat += bes_coef.at(i) * besselj0(alpha_i.at(i) * x) * exp(v_conv * az * CM_TO_KPC / (2. * K)) *
                    sinh((Si / 2.) * (e_diffus - az)) / sinh((Si / 2.) * e_diffus);
      }
      resultat *= (1. / 4. / M_PI) * velocity;
      return resultat; /* [cm^{-2} s^{-1} sr^{-1} GeV^{-1}] */
    }
  }

  real_t Propagation_param_t::GENERIC_FLUX_04(real_t r, real_t z, real_t energy, real_t mass, real_t Z_em,
                                              std::vector<real_t> bes_coef)
  /* GENERIC_FLUX_04 is a generic routine which returns a cosmic-ray flux from Bessel coefficients calculated
  beforehand. The result is modulated by Torsten Bringman's coefficients to speed up convergence. */
  {
    real_t x, az;
    real_t momentum, velocity, K;
    real_t Si, resultat;
    real_t coefficient_torsten = 1.0;

    x = r / R_GAL;
    az = std::abs(z);
    resultat = 0.0;

    real_t e_diffus = equation_parameters.at(2);
    real_t v_conv = equation_parameters.at(3);

    if (az >= e_diffus || x >= 1.0)
    {
      return resultat;
    }
    else
    {
      momentum = sqrt(pow(energy, 2) - pow(mass, 2));
      velocity = CELERITY_LIGHT * momentum / energy;
      K = K_space_diffusion(energy, mass, Z_em);


      for (long int i = 0; i < N_BESSEL; i++)
      {
        /* Si est exprime en [kpc^{-1}]. */
        Si = sqrt(pow(2.0 * alpha_i.at(i) / R_GAL, 2) + pow(v_conv * CM_TO_KPC / K, 2));
        switch (i)
        {
        case (N_BESSEL - 15):
          coefficient_torsten = 624. / 625.;
          break;
        case (N_BESSEL - 14):
          coefficient_torsten = 620. / 625.;
          break;
        case (N_BESSEL - 13):
          coefficient_torsten = 610. / 625.;
          break;
        case (N_BESSEL - 12):
          coefficient_torsten = 590. / 625.;
          break;
        case (N_BESSEL - 11):
          coefficient_torsten = 555. / 625.;
          break;
        case (N_BESSEL - 10):
          coefficient_torsten = 503. / 625.;
          break;
        case (N_BESSEL - 9):
          coefficient_torsten = 435. / 625.;
          break;
        case (N_BESSEL - 8):
          coefficient_torsten = 355. / 625.;
          break;
        case (N_BESSEL - 7):
          coefficient_torsten = 270. / 625.;
          break;
        case (N_BESSEL - 6):
          coefficient_torsten = 190. / 625.;
          break;
        case (N_BESSEL - 5):
          coefficient_torsten = 122. / 625.;
          break;
        case (N_BESSEL - 4):
          coefficient_torsten = 70. / 625.;
          break;
        case (N_BESSEL - 3):
          coefficient_torsten = 35. / 625.;
          break;
        case (N_BESSEL - 2):
          coefficient_torsten = 15. / 625.;
          break;
        case (N_BESSEL - 1):
          coefficient_torsten = 5. / 625.;
          break;
        case (N_BESSEL):
          coefficient_torsten = 1. / 625.;
          break;
        }

        resultat += coefficient_torsten * bes_coef.at(i) * besselj0(alpha_i.at(i) * x) *
                    exp(v_conv * az * CM_TO_KPC / (2. * K)) * sinh((Si / 2.) * (e_diffus - az)) /
                    sinh((Si / 2.) * e_diffus);
      }
      resultat *= (1. / 4. / M_PI) * velocity;

      return resultat; /* [cm^{-2} s^{-1} sr^{-1} GeV^{-1}] */
    }
  }

  /********************************************************************************************/

  void Propagation_param_t::calculation_BESSEL_PBAR_direct_inversion()
  /* We need to solve the differential equation that describes the energy behaviour
  of the BESSEL transforms \bar{\cal{P}}_{i} :

  \bar{A}_{i} \bar{\cal{P}}_{i} \; + \; 2 h \, \frac{\partial}{\partial E} \left\{
  b^{ion} \bar{\cal{P}}_{i} \, - \, D_{E} \, \frac{\partial}{\partial E}
  \bar{\cal{P}}_{i} \right\} \; = \;
  \bar{A}_{i} \, \left\{
  \bar{\cal{P}}_{i}^{pri} + \bar{\cal{P}}_{i}^{sec} + \bar{\cal{P}}_{i}^{ter}
  \right\} \;\; .

  This is a mere diffusion equation. This routine tries to solve that equation
  on the entire kinetic energy range from T_PBAR_MIN up to T_PBAR_MAX. We define
  x[i_pbar] = log(T_pbar/T_PBAR_MIN) so that x[i_pbar] = DELTA_x * (real_t)i_pbar
  with DELTA_x = log(T_PBAR_MAX/T_PBAR_MIN) / (real_t)DIM_TAB_PBAR.

  Then, we define u[i_pbar] = BESSEL_PBAR_TOT_Epbar_i[i_pbar].at(i) for a given BESSEL
  order i. By making it discontinuous, we modify the differential equation to be
  solved into an algebraic relation [A] * [u] = [r]. The matrix [A] is tridiagonal
  so that inversion is a straightforward BUT SOMETIMES HAZARDOUS process! */
  {
    real_t DELTA_x, T_pbar, E_pbar, Abar_i, grand_C_cal;
    std::array<real_t, DIM_TAB_PBAR + 2> a_coeff;
    std::array<real_t, DIM_TAB_PBAR + 1> b_coeff, vec_a, vec_b, vec_c, vec_r, vec_u;

    /* We compute the coefficients a_{j-1/2} = a_coeff.at(j) with j = i_pbar. */
    for (long int i_pbar = 0; i_pbar <= DIM_TAB_PBAR + 1; i_pbar++)
    {
      T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((-0.5 + (real_t)i_pbar) / (real_t)DIM_TAB_PBAR));
      E_pbar = T_pbar + PROTON_MASS;
      a_coeff.at(i_pbar) = D_energy_diffusion(E_pbar, PROTON_MASS, 1.0) / T_pbar; /* in [GeV s^{-1}]. */
    }

    /* We compute the coefficients b_{j} = b_coeff.at(j) with j = i_pbar. */
    for (long int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
      E_pbar = T_pbar + PROTON_MASS;
      b_coeff.at(i_pbar) = b_energy_losses(E_pbar, PROTON_MASS, 1.0); /* in [GeV s^{-1}]. */
    }

    for (long int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    /* BESSEL_PBAR_TOT_Epbar_i is set to zero. */
    {
      for (long int i = 0; i < N_BESSEL; i++)
      {
        bessel_tot_pbar_spec.at(i_pbar).at(i) = 0.0;
      }
    }

    /* BESSEL_PBAR_TOT_Epbar_i is calculated from
    bessel_pri_pbar_spec,  bessel_sec_pbar_spec and
    bessel_ter_pbar_spec by inversing the tridiagonal matrix. */

    DELTA_x = log(T_PBAR_MAX / T_PBAR_MIN) / (real_t)DIM_TAB_PBAR; /* [no unit] */
    for (long int i = 0; i < N_BESSEL; i++)
    {
      for (long int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
      {
        vec_r.at(i_pbar) = bessel_pri_pbar_spec.at(i_pbar).at(i) + bessel_sec_pbar_spec.at(i_pbar).at(i) +
                           bessel_ter_pbar_spec.at(i_pbar).at(i);

        Abar_i = table_abar.at(i_pbar).at(i); /* in [cm s^{-1}]. */
        T_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((real_t)i_pbar / (real_t)DIM_TAB_PBAR));
        grand_C_cal = 2. * E_DISC * CM_TO_KPC / T_pbar / Abar_i;
        /* Fill the tridiagonal matrix before computing the inverse. */
        if (i_pbar == 0)
        {
          vec_a.at(0) = 0.0;
          vec_b.at(0) = 1.0;
          vec_b.at(0) -= grand_C_cal * b_coeff.at(0) / DELTA_x;
          vec_b.at(0) += grand_C_cal * (a_coeff.at(1) - a_coeff.at(0)) / pow(DELTA_x, 2);
          vec_c.at(0) = grand_C_cal * b_coeff.at(1) / DELTA_x;
          vec_c.at(0) -= grand_C_cal * (a_coeff.at(1) - a_coeff.at(0)) / pow(DELTA_x, 2);
        }
        else if (i_pbar == DIM_TAB_PBAR)
        {
          vec_a.at(DIM_TAB_PBAR) = 0.0;
          vec_b.at(DIM_TAB_PBAR) = 1.0;
          vec_c.at(DIM_TAB_PBAR) = 0.0;
        }
        else
        {
          vec_a.at(i_pbar) = -grand_C_cal * b_coeff.at(i_pbar - 1) / 2. / DELTA_x;
          vec_a.at(i_pbar) -= grand_C_cal * a_coeff.at(i_pbar) / pow(DELTA_x, 2);
          vec_b.at(i_pbar) = 1.0;
          vec_b.at(i_pbar) += grand_C_cal * (a_coeff.at(i_pbar) + a_coeff.at(i_pbar + 1)) / pow(DELTA_x, 2);
          vec_c.at(i_pbar) = grand_C_cal * b_coeff.at(i_pbar + 1) / 2. / DELTA_x;
          vec_c.at(i_pbar) -= grand_C_cal * a_coeff.at(i_pbar + 1) / pow(DELTA_x, 2);
        }
      }
      inversion_tridiagonal(vec_a, vec_b, vec_c, vec_r, vec_u);
      for (int i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
      {
        bessel_tot_pbar_spec.at(i_pbar).at(i) = vec_u.at(i_pbar);
      }
    }
  }

  /********************************************************************************************/

  void Propagation_param_t::inversion_tridiagonal(std::array<real_t, DIM_TAB_PBAR + 1>& a,
                                                  std::array<real_t, DIM_TAB_PBAR + 1>& b,
                                                  std::array<real_t, DIM_TAB_PBAR + 1>& c,
                                                  std::array<real_t, DIM_TAB_PBAR + 1>& r,
                                                  std::array<real_t, DIM_TAB_PBAR + 1>& u)
  /* This routine solves the matrix equation [A] * [u] = [r] in the case of a tridiagonal matrix [A]. Pivoting is not
   * used. */
  {
    real_t bet;
    std::array<real_t, DIM_TAB_PBAR + 1> gam;

    if (b.at(0) == 0.0)
    {
      printf(" PROBLEM IN THE TRIDIAGONAL MATRIX INVERSION\n"
             " b.at(0) = %.5e \n",
             b.at(0));
      return;
    }

    bet = b.at(0);
    u.at(0) = r.at(0) / bet;
    for (long int j = 1; j <= DIM_TAB_PBAR; j++)
    {
      gam.at(j) = c.at(j - 1) / bet;
      bet = b.at(j) - a.at(j) * gam.at(j);
      if (bet == 0.0)
      {
        printf(" PROBLEM IN THE TRIDIAGONAL MATRIX INVERSION \n"
               " j = %.ld , beta_j_j = %.5e \n",
               j, bet);
        return;
      }
      u.at(j) = (r.at(j) - a.at(j) * u.at(j - 1)) / bet;
    }
    for (int j = DIM_TAB_PBAR - 1; j >= 0; j--)
    {
      u.at(j) = u.at(j) - gam.at(j + 1) * u.at(j + 1);
    }
  }



  //-----------------------------Terms of the diffusion equation--------------------------//
  real_t Propagation_param_t::K_space_diffusion(real_t energy, real_t mass, real_t Z_em)
  /* Spatial iffusion term in the propagation equation
     K = DIFFUSION_0_GV * beta * pow(rigidity,PUISSANCE_COEFF_DIFF); */
  {
    real_t momentum, rigidity, beta, K;

    momentum = sqrt(pow(energy, 2.) - pow(mass, 2.)); /* [GeV]        */
    rigidity = momentum / Z_em;                       /* [GV]         */
    beta = momentum / energy;                         /* NO DIMENSION */


    /* Spatial diffusion coefficient Strong & al, ApJ 509, 212, 1998. Expressed in [cm^{2} s^{-1}]. */
#ifdef break_in_K
    real_t inv_Rb, Delta_delta, s;

    inv_Rb = 3.205e-3;
    Delta_delta = 0.14;
    s = 0.04;

    K = equation_parameters.at(0) * beta * pow(rigidity, equation_parameters.at(1)) /
        pow(1 + pow((rigidity * inv_Rb), (Delta_delta / (real_t)s)), s);

#else
    K = equation_parameters.at(0) * beta * pow(rigidity, equation_parameters.at(1)); /* [cm^{2} s^{-1}] */
#endif

    return K;
  }

  real_t Propagation_param_t::D_energy_diffusion(real_t energy, real_t mass, real_t Z_em)
  /* This function computes the diffusion coefficient D_EE in energy space that
     is responsible for the diffusive reacceleration. The cosmic ray spectrum gets
     smeared by this process and the average particle energy tends to increase.

     Note that both the TOTAL energy and the TOTAL mass of the cosmic ray species
     are considered here. Whenever one deals with energies per nucleon in the rest
     of the code, conversion into TOTAL energy must be properly performed as this
     function is called.

     Energies and masses are expressed in [GeV]. */
  {
    real_t momentum, rigidity, beta, D_EE;

    momentum = sqrt(pow(energy, 2) - pow(mass, 2)); /* [GeV]        */
    beta = momentum / energy;                       /* NO DIMENSION */

    real_t delta = equation_parameters.at(1);   /* PUISSANCE_COEFF_DIFF */
    real_t v_alpha = equation_parameters.at(4); /* [cm s^{-1}] */
#ifdef DR_energy_drift
    D_EE = 4.0 / 3.0 * pow(v_alpha, 2.0) / K_space_diffusion(energy, mass, Z_em) * pow(beta, 2.0) *
           pow(momentum, 2.0);                                 // [GeV^{2} s^{-1}]
    D_EE /= (delta * (4.0 - pow(delta, 2.0)) * (4.0 - delta)); // [GeV^{2} s^{-1}]
#else
    D_EE = (2. / 9.) * pow(v_alpha, 2) * pow(beta, 4) * pow(energy, 2) /
           K_space_diffusion(energy, mass, Z_em); // [GeV^{2} s^{-1}]
#endif

    return D_EE;
  }

  /********************************************************************************************/

  real_t Propagation_param_t::b_energy_losses(real_t energy, real_t mass, real_t Z_em)
  /* This module computes the energy LOSSES of a given cosmic ray species
     that result from IONIZATION as well as COULOMB FRICTION on the gas of
     the galactic ridge. As the function describes the variation in time of
     the TOTAL energy, note the presence of a MINUS sign.

     Note that both the TOTAL energy and the TOTAL mass of the cosmic ray
     species are considered here. Whenever one deals with energies per nucleon
     in the rest of the code, conversion into TOTAL energy must be properly
     performed whenever this function is called.

     Energies and masses are expressed in [GeV]. */
  {
    real_t momentum, beta, gamma_lorentz;
    real_t Q_max, B1, B2, ionisation;
    real_t x_m, ln_lambda, coulomb, adiabatic;
    real_t dE_dt;
    real_t drift;

    real_t v_conv = equation_parameters.at(3); /* [cm s^{-1}] */
    /* The energy loss dE_dt is expressed in [GeV sec^{-1}]. */

    momentum = sqrt(pow(energy, 2) - pow(mass, 2)); /* [GeV]        */
    beta = momentum / energy;                       /* NO DIMENSION */
    gamma_lorentz = energy / mass;                  /* NO DIMENSION */

    /* IONISATION CONTRIBUTION
       ***********************
       The maximal transfered energy is noted Q_max and expressed in [GeV]. */
    Q_max =
        2. * MASSE_ELECTRON * pow((beta * gamma_lorentz), 2.) / (1.0 + (2. * gamma_lorentz * MASSE_ELECTRON / mass));


    /* B1 and B2 coefficients are dimensionless. */
    B1 = log(2. * MASSE_ELECTRON * pow((beta * gamma_lorentz), 2.) * Q_max / pow(V_ION_H, 2.));
    B1 -= 2. * pow(beta, 2.);
    B2 = log(2. * MASSE_ELECTRON * pow((beta * gamma_lorentz), 2.) * Q_max / pow(V_ION_HE, 2.));
    B2 -= 2. * pow(beta, 2.);

    /* The energy loss from ionisation is expressed in [GeV s^{-1}]. */
    ionisation = (-1.) * 2. * M_PI * pow(RADIUS_ELECTRON, 2) * MASSE_ELECTRON * CELERITY_LIGHT * pow(Z_em, 2) *
                 ((DENSITY_H_DISC * B1) + (DENSITY_HE_DISC * B2)) / beta; /* [GeV s^{-1}] */


    /* COULOMB CONTRIBUTION
     ******************** */
    x_m = pow((3. * sqrt(M_PI) / 4.), (1. / 3.)) * sqrt(2. * K_BOLTZMANN * T_ELECTRONIC / MASSE_ELECTRON);

    ln_lambda = mass / (mass + 2. * gamma_lorentz * MASSE_ELECTRON);
    ln_lambda *= pow((MASSE_ELECTRON / H_BAR / CELERITY_LIGHT), 2.) / RADIUS_ELECTRON; /* [cm^{-3}] */
    ln_lambda *= pow(gamma_lorentz, 2.) * pow(beta, 4.) / M_PI / DENSITY_FREE_ELECTRON;
    ln_lambda = 0.5 * log(ln_lambda);

    coulomb = (-1.) * 4. * M_PI * pow(RADIUS_ELECTRON, 2) * MASSE_ELECTRON * CELERITY_LIGHT * pow(Z_em, 2) *
              DENSITY_FREE_ELECTRON * ln_lambda * pow(beta, 2.) / (pow(x_m, 3.) + pow(beta, 3.)); /* [GeV s^{-1}] */

    /* ADIABATIC CONTRIBUTION
       **********************

       WARNING ! In order to write this term in a form similar to the previous terms, one needs to factorize the term
       2h. hence a term
       \beq
       b_adiab \; = \; - \, E \, \frac{V_{c}}{3h} \;\; .
       \eeq */

    adiabatic = (-1.) * (momentum * momentum / energy) * (v_conv / (3. * E_DISC * CM_TO_KPC)); /* [GeV s^{-1}] */

#ifdef DR_energy_drift
    /* DRIFT TERM
       **********

       This term comes from the formalism based on the phase space distribution function.
       It may be expressed as
       \beq
       b_drift \; = \; \frac{1 \, + \, \beta^{2}}{\beta^{2}} \times \frac{K_EE}{E} \;\; .
       \eeq */
    drift = energy_gain_diff_reac(energy, mass, Z_em);

#else
    drift = 0.0;
#endif

    /* WE SUM UP BOTH CONTRIBUTIONS
     **************************** */
    dE_dt = ionisation + coulomb + adiabatic + drift;
    return dE_dt;
  }

  void Propagation_param_t::ELDR_effect_calculation()
  /* Routines which computes the energy loss and diffusive reacceleration effects */
  {

    for (int i_iteration = 1; i_iteration <= 15; i_iteration++)
    {
      calculation_BESSEL_PBAR_TERTIARY();
      calculation_BESSEL_PBAR_direct_inversion();
    }
    return;
  }

  //-----------------------------------Cosmic ray collision cross sections--------------------------//

  real_t Propagation_param_t::sigma_inelastic_pH_TAN_and_NG(real_t E_proton)
  /* sigma_inelastic_pH_TAN_and_NG is the inelastic scattering cross section in [cm^{2}] of a cosmic-ray proton on a
     hydrogen atom at rest. TAN and NG note this quantity \sigma^{i}_{p}. */
  {
    real_t EK_proton, U, Cp;
    real_t resultat;

    EK_proton = E_proton - PROTON_MASS;
    resultat = 0.0;
    if (EK_proton <= 0.0)
    {
      return resultat;
    }
    else if (EK_proton < 0.3)
    {
      return resultat;
    }
    else if (EK_proton < 3.0)
    {
      U = log(E_proton / 200.0);
      resultat = 1. + 0.0273 * U;
      resultat *= 32.2e-27;
      Cp = 17.9 + 13.8 * log(EK_proton) + 4.41 * log(EK_proton) * log(EK_proton);
      resultat /= 1. + 2.62e-3 * pow(EK_proton, -Cp);
      return resultat;
    }
    else if (E_proton < 200.0)
    {
      U = log(E_proton / 200.0);
      resultat = 1. + 0.0273 * U;
      resultat *= 32.2e-27;
      return resultat;
    }
    else
    {
      U = log(E_proton / 200.0);
      resultat = 1. + 0.0273 * U + 0.01 * U * U;
      resultat *= 32.2e-27;
      return resultat;
    }
  }

  real_t Propagation_param_t::sigma_inelastic_pbarH_TAN_and_NG(real_t E_pbar)
  /* sigma_inelastic_pbarH_TAN_and_NG is the inelastic scattering cross section in [cm^{2}] of a cosmic-ray antiproton
  on a hydrogen atom at rest. sigma_inelastic_pbarH_TAN_and_NG = sigma_ANN_pbarH_TAN_and_NG +
  sigma_inelastic_NOANN_pbarH_TAN_and_NG. TAN et NG write this relation as
  \sigma^{inel}_{pbar} \; = \; \sigma^{an}_{pbar} + \sigma^{i}_{pbar} \;\; .

  sigma_inelastic_pbarH_TAN_and_NG       is noted as \sigma^{inel}_{pbar}
  sigma_ANN_pbarH_TAN_and_NG             is noted as \sigma^{an}_{pbar}
  sigma_inelastic_NOANN_pbarH_TAN_and_NG is noted as \sigma^{i}_{pbar}

  Annihilation dominates at low energy. However, at high energy,  the inelastic interaction becomes non-annihilating.

  Morevover, in the high energy limit, one has the asymptotic behavior
  \sigma^{i}_{pbar} \simeq \sigma^{i}_{p} \;\; . */
  {
    real_t EK_pbar;
    real_t resultat;

    EK_pbar = E_pbar - PROTON_MASS;
    resultat = 0.0;
    if (EK_pbar <= 0.0)
    {
      return resultat;
    }
    else if (EK_pbar < 0.05)
    {
      printf(" WARNING : sigma_inelastic_pbarH_TAN_and_NG NOT DEFINED ! \n");
      return resultat;
    }
    else
    {
      resultat = 1. + 0.584 * pow(EK_pbar, -0.115) + 0.856 * pow(EK_pbar, -0.566);
      resultat *= 24.7e-27;
      return resultat;
    }
  }

  /********************************************************************************************/

  real_t Propagation_param_t::sigma_total_pH(real_t E_proton)
  /*sigma_total_pH is the total cross section of the collision of a cosmic-ray proton on a hydrogen atom at rest. The
   * cross section is in [cm^{2}].*/
  {
    real_t p;
    real_t resultat;

    p = sqrt(pow(E_proton, 2.) - pow(PROTON_MASS, 2.));

    if (p <= p_seuil_pp_tot)
    {
      resultat = P1_pp_tot + P2_pp_tot * log(p) + P3_pp_tot * pow(log(p), 2.) + P4_pp_tot * pow(log(p), 3.) +
                 P5_pp_tot * pow(log(p), 4.) + P6_pp_tot * pow(log(p), 5.) + P7_pp_tot * pow(log(p), 6.);
    }
    else
    {
      resultat = A_pp_tot + B_pp_tot * pow(p, n_pp_tot) + C_pp_tot * pow(log(p), 2.) + D_pp_tot * log(p);
    }
    return (resultat * mb_cm2);
  }

  real_t Propagation_param_t::sigma_inelastic_NOANN_pbarH_TAN_and_NG(real_t E_pbar)
  /* sigma_inelastic_NOANN_pbarH_TAN_and_NG is the non-annihilating inelastic cross section in [cm^{2}]
  of the collision of a cosmic-ray antiproton with an hydrogen atom at rest. */
  {
    real_t EK_pbar, sigma_ANN_pbarH_TAN_and_NG;
    real_t resultat;

    EK_pbar = E_pbar - PROTON_MASS;
    if (EK_pbar <= 0.0)
    {
      sigma_ANN_pbarH_TAN_and_NG = 0.0;
    }
    else if (EK_pbar < 0.05)
    {
      printf(" WARNING : sigma_ANN_pbarH_TAN_and_NG NOT DEFINED ! \n");
      sigma_ANN_pbarH_TAN_and_NG = 0.0;
    }
    else
    {
      sigma_ANN_pbarH_TAN_and_NG = 1. + 0.0115 * pow(EK_pbar, -0.774) - 0.948 * pow(EK_pbar, 0.0151);
      sigma_ANN_pbarH_TAN_and_NG *= 661.0e-27;
    }

    resultat = 0.0;
    if (EK_pbar <= 0.0)
    {
      return resultat;
    }
    else if (EK_pbar < 0.05)
    {
      printf(" WARNING : sigma_inelastic_NOANN_pbarH_TAN_and_NG NOT DEFINED ! \n");
      return resultat;
    }
    else if (EK_pbar <= 13.3)
    {
      resultat = sigma_inelastic_pbarH_TAN_and_NG(E_pbar) - sigma_ANN_pbarH_TAN_and_NG;
      return resultat;
    }
    else
    {
      resultat = sigma_inelastic_pH_TAN_and_NG(E_pbar);
      return resultat;
    }
  }

  //-----------------------------Optimization Methods--------------------------//

  int Propagation_param_t::IDpowell(int n, input_function_t func, const std::vector<real_t>& xtra,
                                    const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& xlim_min,
                                    const std::vector<real_t>& xlim_max, real_t* fmin, std::vector<real_t>& xmin,
                                    real_t ftol, const std::string& option)
  /* finds the minimum (option="min") or maximum (option="max") with a precions ftol of the function func(real_t x[],
   * const std::array<real_t>& xtra) having n free parameters x[] with a lower bound xlim_min and an upper bound
   * xlim_max. The minimum (maximum) value of the function is saved in fmin and the corresponding values of the free
   * parameters in xmin. xtra[] is an array containing the fixed parameters of function func. */
  {
    int nmax = 5;

    srand((unsigned int)time(NULL));
    std::vector<real_t> x0;
    x0.resize(nmax);
    int i;

    real_t factor = 0.;
    if (option == "min")
      factor = 1.;
    else if (option == "max")
      factor = -1.;

    for (i = 0; i < n; i++)
      x0.at(i) = (xlim_min.at(i) + xlim_max.at(i)) / 2.;
    int test = IDpowellaux(n, func, xtra, spect, x0, xlim_min, xlim_max, fmin, xmin, ftol, factor);
    if (option == "max")
      *fmin = -*fmin;

    return test;
  }

  /*----------------------------------------------------------------*/

  int Propagation_param_t::IDpowellaux(int n, input_function_t func, const std::vector<real_t>& xtra,
                                       const std::vector<std::vector<real_t>>& spect, const std::vector<real_t>& x0,
                                       const std::vector<real_t>& xlim_min, const std::vector<real_t>& xlim_max,
                                       real_t* fmin, std::vector<real_t>& xmin, real_t ftol, real_t factor)
  /* finds the minimum (factor=1) or maximum (factor=-1) using a Powell method with a precions ftol of the function
   * func(real_t x[], const std::array<real_t>& xtra) having n free parameters x[] with a lower bound xlim_min and an
   * upper bound xlim_max. The minimum value of the function is saved in fmin and the corresponding values of the free
   * parameters in xmin. xtra[] is an array containing the fixed parameters of function func. The computation starts
   * with the values x0[] of the free parameters x[]. */
  {
    real_t TINY = 1.0e-70;
    int ITMAX = 500;
    int nmax = 5;

    std::vector<std::vector<real_t>> vect; // set of vectors
    vect.resize(nmax);

    /* initializes set of vectors */
    for (int i = 0; i < nmax; i++)
    {
      vect.at(i).resize(nmax);
      for (int j = 0; j < nmax; j++)
      {
        if (i == j)
          vect.at(i).at(j) = 1;
        else
          vect.at(i).at(j) = 0;
      }
    }



    for (int j = 0; j < n; j++)
      xmin.at(j) = x0.at(j); // initializes minimum position

    std::vector<real_t> vecti; // current vector in the iteration
    vecti.resize(nmax);
    real_t fmintemp;
    int imax, j;
    real_t fdeltamax;

    *fmin = factor * (this->*func)(x0, xtra, spect);
    real_t fiter;
    std::vector<real_t> xiter;
    xiter.resize(nmax);
    int testbrent;

    std::vector<real_t> meanvect, xextrapol;
    meanvect.resize(nmax);
    xextrapol.resize(nmax);

    real_t fextrapol;
    real_t t;

    for (int iter = 0; iter <= ITMAX; iter++)
    {
      fiter = *fmin;
      for (j = 0; j < n; j++)
        xiter.at(j) = xmin.at(j);

      fdeltamax = 0;
      for (int i = 0; i < n; i++)
      {
        for (j = 0; j < n; j++)
          vecti.at(j) = vect.at(j).at(i);

        fmintemp = *fmin;
        testbrent = IDbrentmethod(n, func, xtra, spect, xmin, vecti, xlim_min, xlim_max, fmin, factor);
        if (testbrent == 0)
        {
#ifdef DEBUG
          std::cerr << "Brent method unsuccessful\n";
#endif
          return 0;
        }

        if ((fmintemp - *fmin) > fdeltamax)
        {
          imax = i;
          fdeltamax = fmintemp - *fmin;
        }
      }


      if (2. * (fiter - (*fmin)) <= ftol * (fabs(fiter) + fabs(*fmin)) + TINY)
        return 1;

      if (iter == ITMAX)
      {
#ifdef DEBUG
        std::cerr << "Powell exceeding maximum iterations.\n";
#endif
        return 0;
      }

      if (n == 1)
        continue; // no need to check the mean displacement vector in a 1D problem
      for (j = 0; j < n; j++)
      {
        // check if it is worth to keep the mean displacement vector
        meanvect.at(j) = xmin.at(j) - xiter.at(j);
        xextrapol.at(j) = xmin.at(j) + meanvect.at(j);
      }



      fextrapol = factor * (this->*func)(xextrapol, xtra, spect);
      if (fextrapol < fiter)
      {
        t = 2. * (fiter - 2. * (*fmin) + fextrapol) * pow(fiter - (*fmin) - fdeltamax, 2.) -
            fdeltamax * pow(fiter - fextrapol, 2.);


        if (t < 0.)
        {
          IDbrentmethod(n, func, xtra, spect, xmin, meanvect, xlim_min, xlim_max, fmin, factor);

          // vect.at(j).at(n) = meanvect.at(j);
          // Delayed directional update: only update the direction after the next iteration, if the extrapolation is
          // successful. Avoid updating if vect[j][n] has not been used yet, i.e. if iter == 0.
          for (int j = 0; j < n; ++j)
          {
            if (iter > 0)
              vect.at(j).at(imax) = vect.at(j).at(n);
            vect.at(j).at(n) = meanvect.at(j);
          }
        }
      }
    }
    return 1;
  }

  /*--------------------------------------------------------------------*/


  int Propagation_param_t::IDbrentmethod(int n, input_function_t func, const std::vector<real_t>& xtra,
                                         const std::vector<std::vector<real_t>>& spect, std::vector<real_t>& xmin,
                                         const std::vector<real_t>& xi, const std::vector<real_t>& xlim_min,
                                         const std::vector<real_t>& xlim_max, real_t* fmin, real_t factor)
  /* minimizes function func along vector xi */
  {
    real_t a0, b0, x0, xmin1d;
    std::vector<real_t> xitemp;
    xitemp.resize(n);
    int i;

    IDbracket(n, func, xtra, spect, &a0, &x0, &b0, xmin, xi, xlim_min, xlim_max, factor);

    int test = IDbrentmethod1D(func, xtra, spect, a0, x0, b0, fmin, &xmin1d, factor, xmin, xi, n);
    for (i = 0; i < n; i++)
    {
      xmin.at(i) = xmin.at(i) + xmin1d * xi.at(i);
    }
    return test;
  }

  /*--------------------------------------------------------------------*/

  real_t Propagation_param_t::IDbrentmethod1Dfunc(input_function_t func, const std::vector<real_t>& xtra,
                                                  const std::vector<std::vector<real_t>>& spect, real_t x,
                                                  const std::vector<real_t>& xmini, const std::vector<real_t>& xi,
                                                  int n)
  /* function f(x)=func(xmini+x*xi, xtra) */
  {
    int i;
    std::vector<real_t> xitemp;
    xitemp.resize(n);
    for (i = 0; i < n; i++)
      xitemp.at(i) = xmini.at(i) + x * xi.at(i);

    return (this->*func)(xitemp, xtra, spect);
  }

  int Propagation_param_t::IDbrentmethod1D(input_function_t func, const std::vector<real_t>& xtra,
                                           const std::vector<std::vector<real_t>>& spect, real_t a0, real_t x0,
                                           real_t b0, real_t* fmin, real_t* xmin, real_t factor,
                                           const std::vector<real_t>& xmini, std::vector<real_t> xi, int n)
  /* finds the minimum of IDbrentmethod1Dfunc for x between a0 and b0, starting with an initial value x=x0 */
  {
    int brent_iter_max = 1000;
    real_t ZEPS = 1.0e-70;
    real_t CGOLD = 0.1;
    real_t tol = 1.0e-3;

    real_t a, b;

    if (a0 < b0)
    {
      a = a0;
      b = b0;
    }
    else
    {
      a = b0;
      b = a0;
    }

    real_t x, w, v; // points used in the parabolic fit
    real_t u;
    real_t pa1, pa2; // parabola parameters pa2*X^2+pa1*X+pa0
    real_t pamin;    // parabola minimum
    real_t xm;
    real_t fx, fw, fv, fu;
    real_t e = 0, etmp, d;
    int iter = 0;

    x = w = v = x0;
    fx = fw = fv = factor * IDbrentmethod1Dfunc(func, xtra, spect, x0, xmini, xi, n);
    xm = 0.5 * (a + b);

    real_t tol1, tol2;

    do
    {
      // std::cout << "Brent iteration " << iter << ": a= " << a << " b= " << b << " x= " << x << " w= " << w
      //<< " v= " << v << "\n";
      // printf("Brent iteration %d: fx = %f, fw = %f, fv = %f\n", iter, fx, fw, fv);
      tol1 = tol * fabs(x) + ZEPS;
      tol2 = e * tol1;

      if (fabs(e) > tol1)
      {
        // printf("Brent iteration %d: parabolic fit\n", iter);
        pa2 = ((fw - fx) / (w - x) - (fv - fx) / (v - x)) / (w - v); // parabolic fit
        // printf("%f     %f\n", ((fw - fx) / (w - x) - (fv - fx) / (v - x)), (w - v));
        pa1 = (fv - fx) / (v - x) - pa2 * (v + x);
        // printf("%f     %f\n", (fv - fx) / (v - x), (v + x));
        pamin = -0.5 * pa1 / pa2;
        // printf("Brent iteration %d: parabolic fit minimum = %f\n", iter, pamin);
        etmp = e;
        e = d;
        if (fabs(x - pamin) > 0.5 * etmp || pamin < a || pamin > b || !std::isnormal(pamin))
        {
          // printf("Brent iteration %d: parabolic fit not acceptable\n", iter);
          //  cannot use the parabolic fit
          if (x >= xm)
            e = a - x;
          else
            e = b - x;
          d = CGOLD * e;
        }
        else
        {
          // printf("Brent iteration %d: using parabolic fit\n", iter);
          d = pamin - x;
          u = pamin;
          if ((u - a) < tol2 || (b - u) < tol2)
          {
            // pamin too close to borders
            d = fabs(tol1);
            if (xm < x)
              d = -d;
          }
        }
      }
      else
      {
        // printf("Brent iteration %d: parabolic fit not attempted\n", iter);
        if (x >= xm)
          e = a - x;
        else
          e = b - x;
        d = CGOLD * e;
      }

      if (fabs(d) >= tol1)
        u = x + d;
      else
      {
        if (d < 0)
          u = x - fabs(tol1);
        else
          u = x + fabs(tol1);
      }
      fu = factor * IDbrentmethod1Dfunc(func, xtra, spect, u, xmini, xi, n);
      // printf("Brent iteration %d: u = %f, f(u) = %f\n", iter, u, fu);
      if (fu <= fx)
      {
        // printf("Brent iteration %d: new minimum found\n", iter);
        if (u >= x)
          a = x;
        else
          b = x;
        v = w;
        w = x;
        x = u;
        fv = fw;
        fw = fx;
        fx = fu;
      }
      else
      {
        // printf("Brent iteration %d: no improvement\n", iter);
        if (u < x)
          a = u;
        else
          b = u;
        if (fu <= fw || w == x)
        {
          v = w;
          w = u;
          fv = fw;
          fw = fu;
        }
        else if (fu <= fv || v == x || v == w)
        {
          v = u;
          fv = fu;
        }
      }

      xm = 0.5 * (a + b);
      // std::cout << "Brent iteration " << iter << ": a= " << a << " b= " << b << " x= " << x << "\n";

#ifdef DEBUG
      std::vector<real_t> xtemp, utemp;
      xtemp.resize(n);
      utemp.resize(n);
      for (int i = 0; i < n; i++)
      {
        utemp.at(i) = xmini.at(i) + u * xi.at(i);
        xtemp.at(i) = xmini.at(i) + x * xi.at(i);
      }

      std::cerr << "Brent iteration " << iter << ": x = " << x << " f(x) = " << fx << " u = " << u << " f(u) = " << fu
                << "\n";

#endif
      iter++;
      // std::cout << fabs(x - xm) + 0.5 * (b - a) << " " << 2 * tol1 << "\n";
    } while (iter < brent_iter_max && fabs(x - xm) + 0.5 * (b - a) > 2 * tol1);

    *xmin = x;
    *fmin = fx;

    if (iter == brent_iter_max)
    {
      std::cerr << "WARNING: too many iterations in Brent's method\n";
      return 0;
    }
    else
      return 1;
  }

  /*--------------------------------------------------------------------*/

  void Propagation_param_t::IDbracket(int n, input_function_t func, const std::vector<real_t>& xtra,
                                      const std::vector<std::vector<real_t>>& spect, real_t* a0, real_t* x0, real_t* b0,
                                      std::vector<real_t> xinit, std::vector<real_t> xi,
                                      const std::vector<real_t>& xlim_min, const std::vector<real_t>& xlim_max,
                                      real_t factor)
  /* find values  of a0, b0 and x0 required in function IDbrentmethod1D */
  {
    int BRAKETMAX = 50;
    int nmax = 5;

    int i;
    real_t a0tmp, b0tmp;

    real_t xr, frand;
    std::vector<real_t> xrand, xa, xb;
    xrand.resize(nmax);
    xa.resize(nmax);
    xb.resize(nmax);
    real_t fa, fb, fmin;

    int afirst = 1;
    int bfirst = 1;
    bool first = true;

    for (int i = 0; i < n; i++)
    {
      if (xi.at(i) != 0.0)
      {
        real_t l1 = (xlim_min.at(i) - xinit.at(i)) / xi.at(i);
        real_t l2 = (xlim_max.at(i) - xinit.at(i)) / xi.at(i);

        real_t lo = std::min(l1, l2);
        real_t hi = std::max(l1, l2);

        if (first)
        {
          *a0 = lo;
          *b0 = hi;
          first = false;
        }
        else
        {
          *a0 = std::max(*a0, lo);
          *b0 = std::min(*b0, hi);
        }
      }
    }

    if (first || *a0 > *b0)
    {
      *a0 = 0;
      *b0 = 0;
      *x0 = 0;
      return;
    }
    for (i = 0; i < n; i++)
    {
      xa.at(i) = xinit.at(i) + *a0 * xi.at(i);
      xb.at(i) = xinit.at(i) + *b0 * xi.at(i);
    }
    fa = factor * (this->*func)(xa, xtra, spect);
    fb = factor * (this->*func)(xb, xtra, spect);
    // printf("fa = %f, fb = %f\n", fa, fb);
    if (fa < fb)
    {
      *x0 = *a0;
      fmin = fa;
    }
    else
    {
      *x0 = *b0;
      fmin = fb;
    }
    int irand;

    int icount = 0;
    for (irand = 0; irand < BRAKETMAX; irand++)
    {
      xr = *a0 + irand * (*b0 - *a0) / BRAKETMAX;
      for (i = 0; i < n; i++)
        xrand.at(i) = xinit.at(i) + xr * xi.at(i);
      frand = factor * (this->*func)(xrand, xtra, spect);
      if (frand < fmin)
      {
        *x0 = xr;
        fmin = frand;
        icount++;
      }
      if (icount == 5)
        break;
    }
  }


  //-------------------------------CHI SQUARED FUNCTIONS (AMS-02)----------------------------//

  void Propagation_param_t::add_spectra(const std::vector<std::vector<real_t>>& a1,
                                        const std::vector<std::vector<real_t>>& a2,
                                        std::vector<std::vector<real_t>>& add)
  /* sums two spectra a1 and a2 */
  {
    auto a1_size = a1.size();
    auto a2_size = a2.size();
    add.resize(a1_size);

    for (size_t i = 0; i < a1_size; i++)
    {
      add.at(i).resize(2);
      add.at(i).at(0) = a1.at(i).at(0);
      add.at(i).at(1) = a1.at(i).at(1);
      if (a1.at(i).at(0) >= a2.at(0).at(0) && a1.at(i).at(0) <= a2.at(a2_size - 1).at(0))
      {
        add.at(i).at(1) += ind_param.logx_interpol(a2, a1.at(i).at(0));
      }
    }
    return;
  }


  real_t Propagation_param_t::dchi_bkg(const std::vector<real_t>& param, const std::vector<real_t>& logE,
                                       const std::vector<std::vector<real_t>>& spect)
  /* AMS-02 chi^2 with background only */
  {
    std::vector<std::vector<real_t>> back, backS;
    real_t A = param.at(0);
    real_t phi_f = param.at(1);

    if (!background_spectAMS(logE, A, back))
      return -1.e30;
    solar_mod(back, phi_f, backS);
    real_t chi = chi2_AMS(backS);
    return chi;
  }

  std::vector<real_t> Propagation_param_t::chi2_bkg(int& IDpowell_result)
  /* returns AMS-02 delta-chi2 in the hypothesis of no DM */
  {
    std::vector<real_t> logE;
    logE.resize(N_AMS02 + 1);
    for (int i = 0; i < N_AMS02; i++)
      logE.at(i) = log10(AMS02.at(i).at(0));
    logE.at(N_AMS02) = 2.6;

    std::vector<real_t> xlimmin = {0, 0.1};
    std::vector<real_t> xlimmax = {1, 1.1};
    std::vector<real_t> xmin;
    std::vector<std::vector<real_t>> spect;
    xmin.resize(2);

    real_t chi_min;


    IDpowell_result =
        IDpowell(2, &Propagation_param_t::dchi_bkg, logE, spect, xlimmin, xlimmax, &chi_min, xmin, 1.0e-8, "min");

    chi2_noDM = chi_min;

    return xmin;
  }

  real_t Propagation_param_t::chi2_AMS(const std::vector<std::vector<real_t>>& spectrum)
  /* computes the chi^2 between expected spectrum and experimental data */
  {
    real_t chi = 0;

    real_t e, o, sigma;
    for (int i = 0; i < N_AMS02; i++)
    {
      e = spectrum.at(i).at(1);
      o = AMS02.at(i).at(1);
      sigma = AMS02.at(i).at(2);
      chi += (e - o) * (e - o) / (sigma * sigma);
    }

    return chi;
  }

  real_t Propagation_param_t::dchi_tot(const std::vector<real_t>& param, const std::vector<real_t>& logE,
                                       const std::vector<std::vector<real_t>>& spect)
  /* AMS-02 chi^2 for the total (primary+secondary) pbar spectrum */
  {
    std::vector<std::vector<real_t>> tot, back, totS;
    real_t A = param.at(0);
    real_t phi_f = param.at(1);
    if (!background_spectAMS(logE, A, back))
      return -1.e30;
    add_spectra(back, spect, tot);
    solar_mod(tot, phi_f, totS);
    real_t chi2 = chi2_AMS(totS);
    return chi2;
  }

  real_t Propagation_param_t::deltachi2_AMS(int& IDpowell_result, std::vector<real_t>& xmin)
  /* computes AMS-02 delta-chi2 for a DM particle with annihilation spectrum spec and halo and propagation parameters
   * pparam */
  {
    int IDpowell_bkg = 1;
    if (chi2_noDM == 0.)
      auto optimal_params_bkg = chi2_bkg(IDpowell_bkg);
    if (chi2_noDM < -1.e29)
      return -1.e30;

    std::vector<real_t> dEpbar_on_dNpbar_temp, primary_source_temp;
    std::vector<std::vector<real_t>> spect;
    primary_spectra_BCGS_2014(dEpbar_on_dNpbar_temp, primary_source_temp, spect);
    std::cout << "Primary spectrum computed\n";
    real_t chiDM_min = 0;
    std::vector<real_t> logE;
    logE.resize(N_AMS02 + 1);
    for (int i = 0; i < N_AMS02; i++)
      logE.at(i) = log10(AMS02.at(i).at(0));
    logE.at(N_AMS02) = 2.6;

    std::vector<real_t> xlimmin = {0, 0.1};
    std::vector<real_t> xlimmax = {1, 1.1};
    xmin.clear();
    xmin.resize(2);

    IDpowell_result =
        IDpowell(2, &Propagation_param_t::dchi_tot, logE, spect, xlimmin, xlimmax, &chiDM_min, xmin, 1.0e-3, "min");

    if (IDpowell_bkg == 0)
      IDpowell_result = 0;
    return chiDM_min - chi2_noDM;
  }


  //-------------------------------CHI SQUARED FUNCTIONS (FERMI-LAT)----------------------------//

  real_t Propagation_param_t::deltalikelihood_fermi(int& IDpowell_result)
  /* calculates Fermi-LAT dSphs delta log-likelihood for options "conservative", "nominal" or "inclusive" */
  {
    std::vector<real_t> mass_list = get_DM_masses();
    size_t N_MASS = mass_list.size();

    int DM_candidate = input.getLightestBSMpart();
    real_t mass_chi = input.masses_vector.at(DM_candidate);
    if (mass_chi < mass_list.at(0) || mass_chi > mass_list.at(N_MASS - 1))
    {
      printf("Mass out of range, mDM must be between 5 and 1.0e+5 GeV\n");
      return -1;
    }

    std::vector<std::vector<real_t>> production_spectrum;
    ind_param.compute_fluxes(production_spectrum);

    return 2 * likelihood_all_dwarfs(IDpowell_result);
  }

  Propagation_param_t::~Propagation_param_t() = default;



} // namespace __SPEC_LIB_NAME__