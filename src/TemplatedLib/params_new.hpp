/**
 * @brief This file contains the Params_new structure definition
 *
 */

#pragma once

#ifndef PARAMS_NEW__SPEC_LIB_NAME__
  #define PARAMS_NEW__SPEC_LIB_NAME__

  #ifndef __cplusplus
    #define __cplusplus
  #endif
  #include "config.hpp"
  #include "correspondance.hpp"
  #include "cparams.h"
  #include "marty/lha/lha.h"
  #include <chrono>
  #include <cstring>

namespace __SPEC_LIB_NAME__
{

  /**
   * @brief Retrieves a value from an LHA file.
   *
   * @param data The LHA file data.
   * @param block The name of the block in the LHA file.
   * @param N The index of the value within the block.
   * @return real_t The value retrieved from the LHA file.
   * @throws May throw an exception if the value is not found.
   */
  inline real_t get_value(mty::lha::LHAFileData data, const std::string block, const int N)
  {
    return data.getValue(block, N).value();
  }

  /**
   * @brief Retrieves a value from an LHA file with a default fallback.
   *
   * @param data The LHA file data.
   * @param block The name of the block in the LHA file.
   * @param N The index of the value within the block.
   * @param pdgvalue The default value to return if the value is not found in the LHA file.
   * @return real_t The value retrieved from the LHA file or the default pdgvalue.
   */
  inline real_t get_value(mty::lha::LHAFileData data, const std::string block, const int N, const real_t pdgvalue)
  {
    return data.getValue(block, N).value_or(pdgvalue);
  }

  /**
   * @brief Retrieves a matrix diagonal element from an LHA file.
   *
   * @param data The LHA file data.
   * @param block The name of the block in the LHA file.
   * @param N The index of the diagonal element.
   * @return real_t The diagonal element value retrieved from the LHA file.
   * @throws May throw an exception if the value is not found.
   */
  inline real_t get_value_mat(mty::lha::LHAFileData data, const std::string block, const int N)
  {
    return data.getValue(block, N, N).value();
  }

  /**
   * @brief Retrieves a matrix diagonal element from an LHA file with a default fallback.
   *
   * @param data The LHA file data.
   * @param block The name of the block in the LHA file.
   * @param N The index of the diagonal element.
   * @param pdgvalue The default value to return if the value is not found in the LHA file.
   * @return real_t The diagonal element value retrieved from the LHA file or the default pdgvalue.
   */
  inline real_t get_value_mat(mty::lha::LHAFileData data, const std::string block, const int N, const real_t pdgvalue)
  {
    return data.getValue(block, N, N).value_or(pdgvalue);
  }

  namespace mathematicalconstants /// Mathematical constants
  {
    constexpr const real_t pi = M_PI;                     //!< Value of \f$\pi\f$
    constexpr const real_t zeta3 = 3.1415926535897932385; //!< Value of \f$\zeta_3\f$
  } // namespace mathematicalconstants

  namespace units_conversion /// From the exact values of the PDG 2022
  {
    constexpr const real_t c_in_m_s = 299792458.;
    /**< Speed of light in vacuo in m/s.
     */

    constexpr const real_t h_in_J_s = 6.62607015e-34;
    /**< Planck's constant in J*s.
     */

    constexpr const real_t hbar_in_J_s = h_in_J_s / 2. / mathematicalconstants::pi;
    /**< Reduced Planck's constant \f$\hbar\f$ in J*s.
     */

    constexpr const real_t hbar_in_MeV_s = 6.582119569e-22;
    /**< Reduced Planck's constant \f$\hbar\f$ in MeV*s.
     */

    constexpr const real_t hbarc_MeV_fm = 197.3269804;
    /**< \f$\hbar c\f$ in MeV * fm.
     */

    constexpr const real_t hbarc2_Gev2_mbarn = 0.389379372;
    /**< \f$(\hbar c)^2\f$ in \f$\mathrm{GeV}^2 \cdot \mathrm{mbarn}\f$.
     */

    constexpr const real_t e_C = 1.602176634e-19;
    /**< Elementary charge in Coulomb.
     */

    // Energy -> length^-1
    constexpr const real_t MeV_to_fmm1 = 1. / hbarc_MeV_fm;
    /**< Conversion factor from \f$\mathrm{MeV}\f$ to \f$\mathrm{fm}^{-1}\f$.
     */

    constexpr const real_t GeV_to_pmm1 = MeV_to_fmm1;
    /**< Conversion factor from \f$\mathrm{GeV}\f$ to \f$\mathrm{pm}^{-1}\f$.
     */

    constexpr const real_t GeV_to_cmm1 = GeV_to_pmm1 * 1.0e+10;
    /**< Conversion factor from \f$\mathrm{GeV}\f$ to \f$\mathrm{cm}^{-1}\f$.
     */

    constexpr const real_t GeV_to_mm1 = GeV_to_pmm1 * 1.0e+12;
    /**< Conversion factor from \f$\mathrm{GeV}\f$ to \f$\mathrm{m}^{-1}\f$.
     */

    // Energy^-2 -> length^2
    constexpr const real_t GeVm2_to_mbarn = hbarc2_Gev2_mbarn;
    /**< Conversion factor from \f$\mathrm{GeV}^{-2}\f$ to \f$\mathrm{mbarn}\f$.
     */

    constexpr const real_t GeVm2_to_pbarn = GeVm2_to_mbarn * 1.0e+9;
    /**< Conversion factor from \f$\mathrm{GeV}^{-2}\f$ to \f$\mathrm{pbarn}\f$.
     */

    // Energy -> time^-1
    constexpr const real_t GeV_to_sm1 = GeV_to_mm1 / c_in_m_s;
    /**< Conversion factor from \f$\mathrm{GeV}\f$ to \f$\mathrm{s}^{-1}\f$.
     */

    // Need to convert the GeV^-2 into cm^2 (surface) and multiply by c in cm/s to have sigmav in cm^3/s
    constexpr const real_t GeVm2_to_cm3_over_s = c_in_m_s * 1.0e+2 / (GeV_to_cmm1 * GeV_to_cmm1);
    /**< Conversion factor from \f$\mathrm{GeV}^{-2}\f$ to \f$\mathrm{cm}^3/\mathrm{s}\f$.
     */
  } // namespace units_conversion

  namespace pdg2016Value /// Default values from PDG 2016 and their lha block code convention
  {
    constexpr const real_t alpha_em_inv = 137.035999710;
    /**< \f$1/\alpha_{em}\f$ at the scale \f$E = 0\f$
     */

    constexpr const real_t alpha_str_Mz = 0.1181;
    /**< \f$\alpha_s(M_Z)\f$
     *   Corresponsing LHA block: SMINPUTS N=3
     */

    constexpr const real_t mb_mb = 4.18;
    /**< \f$m_b(m_b)\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=5
     */

    constexpr const real_t m_top_pole = 173.34;
    /**< Top quark pole mass.
     *   Corresponding LHA block: SMINPUTS N=6
     *   From the kinematics of ttbar events
     */

    constexpr const real_t m_tau_pole = 1.77686;
    /**< Tau lepton pole mass.
     *   Corresponding LHA block: SMINPUTS N=7
     */

    constexpr const real_t m_el_pole = 510.9989461e-6;
    /**< Electron pole mass in \f$\mathrm{GeV}\f$.
     *   Corresponding LHA block: SMINPUTS N=11
     */

    constexpr const real_t m_mu_pole = 0.105658;
    /**< Muon pole mass.
     *   Corresponding LHA block: SMINPUTS N=13
     */

    constexpr const real_t m_d2GeV = 4.7e-3;
    /**< \f$m_d(2 \mathrm{GeV})\f$: down mass in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=21
     */

    constexpr const real_t m_u2GeV = 2.2e-3;
    /**< \f$m_u(2 \mathrm{GeV})\f$: up mass in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=22
     */

    constexpr const real_t m_s2GeV = 9.6e-2;
    /**< \f$m_s(2\mathrm{GeV})\f$: strange in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=23
     */

    constexpr const real_t m_c_m_c = 1.27;
    /**< \f$m_c(m_c)\f$: charm mass in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=24
     */

    constexpr const real_t m_W = 80.385;
    /**< W boson mass.
     *   Corresponding LHA block: MASS N=24
     */

    constexpr const real_t m_Z = 91.1876;
    /**< Z boson mass.
     *   Corresponding LHA block: SMINPUTS N=4
     */

    constexpr const real_t theta_W = 0.490449;
    /**< Weinberg angle.
     */

    constexpr const real_t delta_wolf = 1.20146;
    /**< CKM CP-violating phase.
     */

    constexpr const real_t V_ud = 0.974339;
    /**< CKM matrix element \f$V_{ud}\f$.
     */

    constexpr const real_t V_us = 0.225057;
    /**< CKM matrix element \f$V_{us}\f$.
     */

    constexpr const real_t V_ub_mod = 0.00366661;
    /**< Magnitude of CKM matrix element \f$|V_{ub}|\f$.
     */

    constexpr const complex_t V_cd = {-0.224922, -0.000136891};
    /**< CKM matrix element \f$V_{cd}\f$.
     */

    constexpr const complex_t V_cs = {0.97351, -3.16198e-05};
    /**< CKM matrix element \f$V_{cs}\f$.
     */

    constexpr const real_t V_cb = 0.0410782;
    /**< CKM matrix element \f$V_{cb}\f$.
     */

    constexpr const complex_t V_td = {0.00795853, -0.00332961};
    /**< CKM matrix element \f$V_{td}\f$.
     */

    constexpr const complex_t V_ts = {-0.0403218, -0.00076909};
    /**< CKM matrix element \f$V_{ts}\f$.
     */

    constexpr const real_t V_tb = 0.999149;
    /**< CKM matrix element \f$V_{tb}\f$.
     */

    constexpr const real_t Gamma_Z = 2.4952;
    /**< Z boson decay width.
     */

    constexpr const real_t Gamma_W = 2.085;
    /**< W boson decay width.
     */

    constexpr const real_t m_e = 510.9989461e-6;
    /**< Electron mass in \f$\mathrm{GeV}\f$.
     */

    constexpr const real_t k_B = 8.617330e-5;
    /**< Boltzmann's constant in \f$\mathrm{GeV}/\mathrm{GK}\f$.
     */

    constexpr const real_t alphaem = 1. / alpha_em_inv;
    /**< Fine-structure constant
     *   (approximately 0.007297353).
     */

    constexpr const real_t sigma_SB = 0.16449340668482282;
    /**< Stefan-Boltzmann constant \f$(\pi^2/60)\f$.
     */

    constexpr const real_t Gn = 6.67428e-8;
    /**< Newton's constant in \f$\mathrm{cm}^3\cdot\mathrm{g}^{-1}\cdot\mathrm{s}^{-2}\f$.
     */

    constexpr const real_t G = 6.708609142443796e-39;
    /**< Newton constant in \f$\mathrm{GeV}^{-2}\f$.
     */

    constexpr const real_t Mplanck = 1.2209102930946623e+19;
    /**< Planck mass in \f$\mathrm{GeV}\f$.
     */

    constexpr const real_t DMpn = 0.0012934;
    /**< Mass difference between neutron and proton in \f$\mathrm{GeV}\f$.
     */

    constexpr const real_t M_u = 1.660538921e-24;
    /**< Unit atomic mass \f$M_u\f$ in \f$\mathrm{g}\f$.
     */

    constexpr const real_t zeta = 1.6103162253325862;
    /**< \f$3k_B/(2c^2M_u)\f$ in \f$\mathrm{GeV}^{-1}\f$.
     */

    constexpr const real_t K_to_eV = k_B;
    /**< Conversion factor: \f$T(10^9 \mathrm{K}) \cdot K\_to\_GeV = T(\mathrm{GeV})\f$ or \f$T(\mathrm{K}) \cdot
     * K\_to\_eV = T(\mathrm{eV})\f$.
     */

    constexpr const real_t g_to_GeV = 5.60958884538932e+23;
    /**< Conversion factor from grams to GeVs:
     *   \f$M(\mathrm{g}) \cdot g\_to\_GeV = M(\mathrm{GeV})\f$.
     */

    constexpr const real_t kg_to_GeV = 5.60958884538932e+26;
    /**< Conversion factor: from kg to GeV
     *   \f$M(\mathrm{kg}) \cdot kg\_to\_GeV = M(\mathrm{GeV})\f$.
     */

    constexpr const real_t m_to_GeV = 5.067730582705779e+15;
    /**< Conversion factor: from meters to \f$\mathrm{GeV}^{-1}\f$
     *   \f$L(\mathrm{m}) \cdot m\_to\_GeV = L(\mathrm{GeV}^{-1})\f$.
     */

    constexpr const real_t cm_to_GeV = 5.067730582705779e+13;
    /**< Conversion factor from cm to \f$\mathrm{GeV}^{-1}\f$:
     *   \f$L(\mathrm{cm}) \cdot cm\_to\_GeV = L(\mathrm{GeV}^{-1})\f$.
     */

    constexpr const real_t s_to_GeV = 1.5192674078711377e+24;
    /**< Conversion factor from s to \f$\mathrm{GeV}^{-1}\f$:
     *   \f$t(\mathrm{s}) \cdot s\_to\_GeV = t(\mathrm{GeV}^{-1})\f$.
     */
  } // namespace pdg2016Value

  namespace pdgValue /// Default values from PDG 2022 and their lha block code convention
  {
    /* In the block SM_INPUT
     * a set of SM parameter is provided in order to compute
     * the spectrum and other quantities
     *
     * In the block MASS
     * the spectrum (pole masses) has to be provided
     *
     * We list here the reference values taken from the PDG 2022
     * They are chosen if no alternative value is provided in the lha file
     */

    constexpr const real_t alpha_em_inv = 137.035999180;
    /**< \f$1/\alpha_{em}(E=0)\f$ at the scale \f$E = 0\f$
     */

    constexpr const real_t alpha_em_Mtau_inv = 133.471;
    /**< \f$1/\alpha_{em}(m_{\tau})\f$ at the scale \f$E = m_{\tau}\f$ in the MS bar scheme.
     */

    constexpr const real_t alpha_em_MZ_inv = 127.951;
    /**< \f$1/\alpha_{em}(M_Z)\f$ at the scale \f$E = M_Z\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=1
     */

    constexpr const real_t Gfermi = 1.1663788e-5;
    /**< Fermi constant in \f$\mathrm{GeV}^{-2}\f$.
     *   Corresponding LHA block: SMINPUTS N=2.
     */

    constexpr const real_t alpha_str_Mz = 0.1179;
    /**< \f$\alpha_s(M_Z)\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=3
     */

    constexpr const real_t m_Z_pole = 91.1876 - 0.034;
    /**< Z pole mass
     *   Corresponding LHA block: SMINPUTS N=4
     */

    constexpr const real_t m_Z_m_Z = 91.1876;
    /**< \f$M_Z(M_Z)\f$ in the MS bar scheme.
     *   Corresponding LHA block: MASS N=23.
     *   Experimental value from Breit-Wigner distribution \f$(M_Z(M_Z))\f$
     */

    constexpr const real_t mb_mb = 4.18;
    /**< \f$m_b(m_b)\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=5
     */

    constexpr const real_t mb_pole = 4.78;
    /**< \f$m_b\f$ pole computed from \f$m_b(m_b)\f$.
     *   Corresponding LHA block: MASS N=5
     */

    constexpr const real_t m_top_pole = 172.69;
    /**< Top quark pole mass.
     *   Corresponding LHA block: SMINPUTS N=6
     *   From the kinematics of ttbar events
     */

    constexpr const real_t m_top_m_top_computed = 163.058;
    /**< \f$m_t(m_t)\f$ computed at 2 loops from the previously defined m_top_pole.
     */

    constexpr const real_t m_tau_pole = 1.77686;
    /**< Tau lepton pole mass.
     *   Corresponding LHA block: SMINPUTS N=7
     */

    constexpr const real_t m_el_pole = 5.1099895000e-4;
    /**< Electron pole mass.
     *   Corresponding LHA block: SMINPUTS N=11
     */

    constexpr const real_t m_mu_pole = 0.1056583755;
    /**< Muon pole mass.
     *   Corresponding LHA block: SMINPUTS N=13
     */

    constexpr const real_t m_d2GeV = 4.67e-3;
    /**< \f$m_d(2 \mathrm{GeV})\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=21
     */

    constexpr const real_t m_u2GeV = 2.16e-3;
    /**< \f$m_u(2 \mathrm{GeV})\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=22
     */

    constexpr const real_t m_s2GeV = 9.34e-2;
    /**< \f$m_s(2 \mathrm{GeV})\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=23
     */

    constexpr const real_t m_c_m_c = 1.27;
    /**< \f$m_c(m_c)\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=24
     */

    constexpr const real_t m_c_pole_computed = 1.66949;
    /**< Pole mass of the c-quarks computed at 2 loops from the previously defined m_c_m_c.
     */

    constexpr const real_t ckm_lambda = 0.22500;
    /**< CKM Wolfenstein parameter \f$\lambda\f$.
     *   Corresponding LHA block: VCKMIN N=1
     */

    constexpr const real_t ckm_A = 0.826;
    /**< CKM Wolfenstein parameter A.
     *   Corresponding LHA block: VCKMIN N=2
     */

    constexpr const real_t ckm_rhobar = 0.159;
    /**< CKM Wolfenstein parameter \f$\bar{\rho}\f$.
     *   Corresponding LHA block: VCKMIN N=3
     */

    constexpr const real_t ckm_etabar = 0.348;
    /**< CKM Wolfenstein parameter \f$\bar{\eta}\f$.
     */

    // CKM matrix assignment from Wolfenstein values
    /***************************************************************/

    const inline static real_t s12 = ckm_lambda;
    /**< CKM standard parameter \f$sin(\theta_{12})\f$.
     */

    const inline static real_t s23 = ckm_A * ckm_lambda * ckm_lambda;
    /**< CKM standard parameter \f$sin(\theta_{23})\f$.
     */

    const inline static complex_t rhopIeta = {ckm_rhobar, ckm_etabar};
    /**< CKM parameter \f$\rho + \marhrm{i} \eta\f$.
     */

    const inline static real_t s13 = std::abs(
        ckm_A * std::pow(ckm_lambda, 3.) * (rhopIeta)*std::sqrt(1. - std::pow(ckm_A * ckm_lambda * ckm_lambda, 2.)) /
        std::sqrt(1. - ckm_lambda * ckm_lambda) / (1. - std::pow(ckm_A * ckm_lambda * ckm_lambda, 2.) * (rhopIeta)));
    /**< CKM standard parameter \f$sin(\theta_{13})\f$.
     */

    const inline static complex_t expid =
        ((complex_t)(ckm_A)*std::pow(ckm_lambda, 3.) *
         (rhopIeta)*std::sqrt(1. - std::pow(ckm_A * ckm_lambda * ckm_lambda, 2.)) /
         std::sqrt(1. - ckm_lambda * ckm_lambda) / (1. - std::pow(ckm_A * ckm_lambda * ckm_lambda, 2.) * (rhopIeta))) /
        s13;
    /**< CKM standard parameter \f$\mathrm{e}^{\mathrm{i}\delta}\f$
     */

    const inline static real_t c12 = std::sqrt(1. - s12 * s12);
    /**< CKM standard parameter \f$cos(\theta_{12})\f$
     */

    const inline static real_t c13 = std::sqrt(1. - s13 * s13);
    /**< CKM standard parameter \f$cos(\theta_{13})\f$
     */

    const inline static real_t c23 = std::sqrt(1. - s23 * s23);
    /**< CKM standard parameter \f$cos(\theta_{23})\f$
     */

    const inline static complex_t CKM[4][4] = {
        {0, 0, 0, 0},
        {0, c12* c13, c12* c13, s13 / expid},
        {0, -s12* c23 - c12* s23* s13* expid, c12* c23 - s12* s23* s13* expid, s23* c13},
        {0, s12* s23 - c12* c23* s13* expid, -c12* s23 - s12* c23* s13* expid, c23* c13}};
    /**< C-style matrix with CKM matrix values.
     *   Indexing starts from 1.
     */


    const real_t V_ud = CKM[1][1].real();
    /**< CKM matrix element \f$V_{ud}\f$.
     */

    const real_t V_us = CKM[1][2].real();
    /**< CKM matrix element \f$V_{us}\f$.
     */

    const real_t V_ub_mod = s13;
    /**< Module of CKM matrix element \f$|V_{ub}|\f$.
     */

    const real_t delta_wolf = std::log(expid).imag();
    /**< CKM phase \f$\delta\f$.
     */

    const complex_t V_cd = CKM[2][1];
    /**< CKM matrix element \f$V_{cd}\f$.
     */

    const complex_t V_cs = CKM[2][2];
    /**< CKM matrix element \f$V_{cs}\f$.
     */
    const real_t V_cb = CKM[2][3].real();
    /**< CKM matrix element \f$V_{cb}\f$.
     */

    const complex_t V_td = CKM[3][1];
    /**< CKM matrix element \f$V_{td}\f$.
     */

    const complex_t V_ts = CKM[3][2];
    /**< CKM matrix element \f$V_{ts}\f$.
     */

    const real_t V_tb = CKM[3][3].real();
    /**< CKM matrix element \f$V_{tb}\f$.
     */
    /**************************************************************/

    constexpr const real_t m_h = 125.25;
    /**< SM Higgs boson mass
     */

    constexpr const real_t Gamma_Z = 2.4952;
    /**< Z boson decay width
     */

    constexpr const real_t Gamma_W = 2.085;
    /**< W boson decay width
     */

    constexpr const real_t m_e = m_el_pole;
    /**< Electron mass (same as electron pole mass)
     */

    constexpr const real_t k_B = 8.617333262e-5;
    /**< Boltzmann constant in \f$\mathrm{eV}/\mathrm{K}\f$
     */

    const real_t alphaem = 1. / alpha_em_inv;
    /**< Fine-structure constant at E=0
     */

    constexpr const real_t Gn = 6.67430e-8;
    /**< Newton's constant in \f$\mathrm{cm}^3\cdot\mathrm{g}^{-1}\cdot\mathrm{s}^{-2}\f$
     */

    constexpr const real_t G = 6.70883e-39;
    /**< Newton constant in \f$\mathrm{GeV}^{-2}\f$
     */

    constexpr const real_t Mplanck = 1.220890e+19;
    /**< Planck mass in \f$\mathrm{GeV}\f$
     */

    constexpr const real_t M_u = 1.6605306660e-24;
    /**< Unit atomic mass \f$M_u\f$ in \f$\mathrm{g}\f$
     */

    constexpr const real_t K_to_eV = k_B;
    /**< Conversion factor:
     * \f$T(10^9 \mathrm{K}) \cdot K\_to\_GeV = T(\mathrm{GeV})\f$ or \f$T(\mathrm{K}) \cdot K\_to\_eV =
     * T(\mathrm{eV})\f$
     */

    constexpr const real_t kg_to_eV = 5.609588603e+35;
    /**< Conversion factor from \f$\mathrm{kg}\f$ to \f$\mathrm{eV}\f$
     */

    constexpr const real_t kg_to_GeV = kg_to_eV * 1.0e-9;
    /**< Conversion factor from \f$\mathrm{kg}\f$ to \f$\mathrm{GeV}\f$
     */

    constexpr const real_t g_to_GeV = 1.0e-3 * kg_to_eV * 1.0e-9;
    /**< Conversion factor from \f$\mathrm{g}\f$ to \f$\mathrm{GeV}\f$
     */

    constexpr const real_t hbarc = 197.3269804;
    /**< \f$\hbar c\f$ in \f$\mathrm{MeV}\cdot\mathrm{fm}\f$
     */

    constexpr const real_t fm_to_eVm1 = 1. / hbarc * 1.0e-6;
    /**< Conversion factor from \f$\mathrm{fm}\f$ to \f$\mathrm{eV}^{-1}\f$
     */

    constexpr const real_t m_to_GeV = 1.0e+15 * fm_to_eVm1 * 1.0e+9;
    /**< Conversion factor from \f$\mathrm{m}\f$ to \f$\mathrm{GeV}^{-1}\f$
     */

    constexpr const real_t cm_to_GeV = 1.0e-2 * m_to_GeV;
    /**< Conversion factor from \f$\mathrm{cm}\f$ to \f$\mathrm{GeV}^{-1}\f$
     */

    constexpr const real_t c_em = 299792458;
    /**< Speed of light in \f$\mathrm{m}/\mathrm{s}\f$
     */

    const real_t s_to_GeV = c_em * m_to_GeV;
    /**< Conversion factor from \f$\mathrm{s}\f$ to \f$\mathrm{GeV}^{-1}\f$
     */

    constexpr const real_t hbar = 6.582119569e-22;
    /**< \f$\hbar\f$ in \f$\mathrm{MeV}\cdot\mathrm{s}\f$
     */

    const real_t sigma_SB = pdg2016Value::sigma_SB;
    /**< Stefan-Boltzmann constant
     */

    const real_t zeta = pdg2016Value::zeta;
    /**< \f$3k_B/(2c^2M_u)\f$ in \f$\mathrm{GeV}^{-1}\f$
     */

    const real_t DMpn = pdg2016Value::DMpn;
    /**< Mass difference between neutron and proton in \f$\mathrm{GeV}\f$
     */
  } // namespace pdgValue

  struct darkpackparam_t
  {
public:
    /* Parameters necessary for the running */
    csl::InitSanitizer<real_t> Running_scale{"Running_scale"};
    /**< Running scale parameter.
     */

    csl::InitSanitizer<real_t> alpha_str_Mz{"alpha_strong(M_Z)"};
    /**< \f$\alpha_s(M_Z)\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=3
     */

    csl::InitSanitizer<real_t> m_Z_pole{"m_Z,pole"};
    /**< Z boson pole mass.
     *   Corresponding LHA block: SMINPUTS N=4
     */

    csl::InitSanitizer<real_t> mbmb{"m_b(m_b)"};
    /**< \f$m_b(m_b)\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=5
     */

    csl::InitSanitizer<real_t> mtop_pole{"m_top,pole"};
    /**< Top quark pole mass.
     *   Corresponding LHA block: SMINPUTS N=6
     */

    /* Parameters as read from lha file */

    csl::InitSanitizer<real_t> alpha_em_MZ_inv{"alpha_em_MZ_inv"};
    /**< \f$1/\alpha_{em}(M_Z)\f$ in the MS bar scheme.
     *   Corresponding LHA block: SMINPUTS N=1
     */

    csl::InitSanitizer<real_t> Gfermi{"G_F"};
    /**< Fermi constant in \f$\mathrm{GeV}^{-2}\f$.
     *   Corresponding LHA block: SMINPUTS N=2
     */

    csl::InitSanitizer<real_t> m_tau_pole{"m_tau_pole"};
    /**< Tau lepton pole mass.
     *   Corresponding LHA block: SMINPUTS N=7
     */

    csl::InitSanitizer<real_t> m_el_pole{"m_el_pole"};
    /**< Electron pole mass.
     *   Corresponding LHA block: SMINPUTS N=11
     */

    csl::InitSanitizer<real_t> m_mu_pole{"m_mu_pole"};
    /**< Muon pole mass.
     *   Corresponding LHA block: SMINPUTS N=13
     */

    csl::InitSanitizer<real_t> m_d2GeV{"m_d2GeV"};
    /**< \f$m_d(2 \mathrm{GeV})\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=21
     */

    csl::InitSanitizer<real_t> m_u2GeV{"m_u2GeV"};
    /**< \f$m_u(2 \mathrm{GeV})\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=22
     */

    csl::InitSanitizer<real_t> m_s2GeV{"m_s2GeV"};
    /**< \f$m_s(2 \mathrm{GeV})\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=23
     */

    csl::InitSanitizer<real_t> m_c_m_c{"m_c_m_c"};
    /**< \f$m_c(m_c)\f$ in the MSbar scheme.
     *   Corresponding LHA block: SMINPUTS N=24
     */

    csl::InitSanitizer<real_t> ckm_lambda{"ckm_lambda"};
    /**< CKM Wolfenstein parameter \f$\lambda\f$.
     *   Corresponding LHA block: VCKMIN N=1
     */

    csl::InitSanitizer<real_t> ckm_A{"ckm_A"};
    /**< CKM Wolfenstein parameter A.
     *   Corresponding LHA block: VCKMIN N=2
     */

    csl::InitSanitizer<real_t> ckm_rhobar{"ckm_rhobar"};
    /**< CKM Wolfenstein parameter \f$\bar{\rho}\f$.
     *   Corresponding LHA block: VCKMIN N=3
     */

    csl::InitSanitizer<real_t> ckm_etabar{"ckm_etabar"};
    /**< CKM Wolfenstein parameter \f$\bar{\eta}\f$.
     *   Corresponding LHA block: VCKMIN N=4
     */

    bool aremassesassigned = false;
    /**< Boolean that determines if the masses_vector array has been initialized.
     */

    bool arewidthsassigned = false;
    /**< Boolean that determines if the width_vector array has been initialized.
     */

    std::array<csl::InitSanitizer<real_t>, corr::TOTAL_PARTICLES + 1> masses_vector;
    /**< Array containing the masses of all particles.
     *   Ordered according to the corr::Part_t enumeration in correspondance.h.
     */

    std::array<csl::InitSanitizer<real_t>, corr::TOTAL_PARTICLES + 1> widths_vector;
    /**< Array containing the widths of all particles.
     *   Ordered according to the corr::Part_t enumeration in correspondance.h.
     */

protected:
    csl::InitSanitizer<int> part_lbsm{"part_lbsm"};
    /**< Integer value labeling the lightest BSM particle.
     *   Uses the corr::Part_t enumeration in correspondance.h.
     */

    csl::InitSanitizer<real_t> m_lbsm{"m_lbsm"};
    /**< Mass of the lightest BSM particle.
     */

    csl::InitSanitizer<real_t> em_charge_lbsm{"q_lbsm"};
    /**< Electric charge of the lightest BSM particle.
     */

    csl::InitSanitizer<int> g_lbsm{"g_lbsm"};
    /**< Spin-helicity degrees of freedom of the lightest BSM particle.
     */

public:
    /** @brief Initializes the object with empty values.
     *  Sets mass and width assignment flags to false.
     *  Assigns names to mass and width variables.
     */
    inline void empty_init()
    {
      aremassesassigned = false;
      arewidthsassigned = false;
      for (size_t i = 0; i < corr::TOTAL_PARTICLES + 1; i++)
      {
        // Here we set the names of the csl::InitSanitizer variables
        // which are contained in arrays
        masses_vector[i].setName("m_" + corr::part_names[i]);
        widths_vector[i].setName("Gamma_" + corr::part_names[i]);
      }
    }

    /** @brief Default constructor.
     *  Calls empty_init() to initialize the object.
     */
    darkpackparam_t() { empty_init(); }

    /** @brief Returns the integer value in the corr::Part_t enumeration
     *         corresponding to the lightest BSM particle.
     *  @return Integer value of the lightest BSM particle.
     */
    inline int getLightestBSMpart() const { return part_lbsm; }

    /** @brief Returns the mass of the lightest BSM particle.
     *  @return Mass of the lightest BSM particle as a real_t.
     */
    inline real_t getLightestBSMmass() const { return m_lbsm; }

    /** @brief Returns the electric charge of the lightest BSM particle.
     *  @return Electric charge of the lightest BSM particle as a real_t.
     */
    inline real_t getLightestBSMem_charge() const { return em_charge_lbsm; }

    /** @brief Returns the spin-helicity degrees of freedom of the lightest BSM particle.
     *  @return Degrees of freedom of the lightest BSM particle as an integer.
     */
    inline int getLightestBSMdof() const { return g_lbsm; }

    /** @brief Compares this object with another darkpackparam_t object.
     *  @param other The darkpackparam_t object to compare with.
     *  @return True if objects are equal, false otherwise.
     */
    bool compare_with_darkpackparam(const darkpackparam_t&) const;

    /** @brief Equality comparison operator.
     *  @param other The darkpackparam_t object to compare with.
     *  @return True if objects are equal, false otherwise.
     */
    virtual bool operator==(const darkpackparam_t&) const;

    /** @brief Inequality comparison operator.
     *  @param other The darkpackparam_t object to compare with.
     *  @return True if objects are not equal, false otherwise.
     */
    virtual bool operator!=(const darkpackparam_t& other) const { return !(*this == other); }

    /** @brief Prints the object's data to the specified output stream.
     *  @param out The output stream to print to (default is std::cout).
     */
    void print(std::ostream& out = std::cout) const;
  };


  struct Param_t : public darkpackparam_t, param_t
  {
public:
    /**  @brief Reads Standard Model parameters from a file.
     *   @param namefile The name of the file containing SM parameters.
     */
    void ReadSMparams(const std::string namefile);

    /** @brief Default constructor.
     *  Initializes the object using the darkpackparam_t constructor.
     */
    Param_t() : darkpackparam_t() {};

    /**  @brief Constructor that reads parameters from an LHA file.
     *   @param filename The name of the LHA file to read from.
     *   @note If leshouchesfrommarty.h/.cpp is not in model-dependent libraries,
     *         the installation script creates an empty library.
     */
    Param_t(const std::string filename);

    /**  @brief Assigns values to the masses_vector array.
     *   Uses corr::getMassFirst to populate the array.
     */
    void AssignMassesVector();

    /**  @brief Determines the lightest BSM particle and assigns relevant data.
     *   Updates the corresponding data member in this instance.
     */
    void setLightestBSMparticle();

    /** @brief Updates the object state after an input change.
     */
    void refresh();

    /**  @brief Modifies the object's content based on another instance.
     *   Checks that SM parameters and DM candidate remain the same.
     *   @param other The Param_t object to copy values from.
     *   @return 0 for success, non-zero for failure.
     */
    unsigned int change_values(const Param_t& other);

    using darkpackparam_t::operator==; // Expose base class method
    using darkpackparam_t::operator!=; // Expose base class method

    /**  @brief Equality comparison operator.
     *   @param other The Param_t object to compare with.
     *   @return True if objects are equal, false otherwise.
     */
    bool operator==(const Param_t&) const;

    /**  @brief Inequality comparison operator.
     *   @param other The Param_t object to compare with.
     *   @return True if objects are not equal, false otherwise.
     */
    bool operator!=(const Param_t& other) const { return !(*this == other); }

    /**  @brief Prints the object's data to the specified output stream.
     *   @param out The output stream to print to (default is std::cout).
     */
    void Print(std::ostream& out = std::cout) const;
  };

  class Timer /// The following class is defined to have an easy way of measuring time
  {
    std::chrono::_V2::system_clock::time_point start, end;

    std::chrono::duration<float> duration;

public:
    Timer() { start = std::chrono::high_resolution_clock::now(); }

    real_t getElapsed()
    {
      end = std::chrono::high_resolution_clock::now();
      duration = end - start;
      return static_cast<real_t>(duration.count());
    }

    ~Timer()
    {
      end = std::chrono::high_resolution_clock::now();
      duration = end - start;
      std::cout << "Timer took " << duration.count() * 1000 << " ms" << std::endl;
    }
  };

} // end of namespace __SPEC_LIB_NAME__

/**  @brief Stream insertion operator for param_t objects.
 *   Outputs the param_t object to the given output stream.
 *   @param out The output stream to write to.
 *   @param par The param_t object to output.
 *   @return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::param_t& par)
{
  par.print(out);
  return out;
}

/**  @brief Stream insertion operator for Param_t objects.
 *   Outputs the Param_t object to the given output stream.
 *   @param out The output stream to write to.
 *   @param par The Param_t object to output.
 *   @return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& out, const __SPEC_LIB_NAME__::Param_t& par)
{
  par.Print(out);
  return out;
}
#endif // __SPEC_LIB_NAME__PARAMS_H