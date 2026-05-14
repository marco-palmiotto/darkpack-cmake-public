#pragma once
#ifndef __cplusplus
  #define __cplusplus
#endif

#include <gtest/gtest.h>

#include "dp_mssm2to2/RunningSM.hpp"
#include "dp_mssm2to2/indirect.hpp"
#include "dp_mssm2to2/process.hpp"
#include "dp_mssm2to2/process_1to2.hpp"
#include "dp_mssm2to2/propagation.hpp"
#include "dp_mssm2to2/relicparam.hpp"

namespace mssm2to2
{

  real_t Propagation_param_t::test_func(const std::vector<real_t>& x, const std::vector<real_t>& xtra,
                                        const std::vector<std::vector<real_t>>& spect)
  {
    real_t a = xtra[0];
    real_t b = xtra[1];
    real_t c = xtra[2];
    real_t d = xtra[3];
    real_t e = xtra[4];

    return a * pow(x[0], 2.) + c * pow(x[1], 2.) + b * x[0] * x[1] + d * x[0] + e * x[1];
  };

  real_t Propagation_param_t::test_func_spec(const std::vector<real_t>& x, const std::vector<real_t>& xtra,
                                             const std::vector<std::vector<real_t>>& spect)
  {
    real_t result = 0;
    for (int i = 0; i < spect.size(); i++)
    {
      for (int j = 0; j < spect[i].size(); j++)
      {
        result += spect[i][j] * (pow(x[0], i + 2));
      }
    }

    result /= xtra[0] * (x[0] - 1);

    return result;
  };
} // namespace mssm2to2

using namespace mssm2to2;
using namespace advmath;

TEST(PropagationParamTest, Changing_mass)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::BURKERT);

  real_t mass_chi_min = 100.;
  real_t mass_chi_max = 10000.;
  real_t dm_mass = mass_chi_min * pow((mass_chi_max / mass_chi_min), ((real_t)(27 - 1) / 100.));
  prop_param.change_ind_param(0., dm_mass);
  int powell;
  real_t delta_likelihood1 = prop_param.deltalikelihood_fermi(powell);
  std::cout << "Delta likelihood 1= " << delta_likelihood1 << std::endl;
}

TEST(IndirectParamTest, Initialization)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  EXPECT_NO_THROW({ Indirect_param_t ind_param(input); });
}

/*
TEST(IndirectParamTest, ReadingFile)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);

  std::vector<std::vector<real_t>> dummy_data;
  ind_param.read_file("/workspaces/darkpack-cmake/src/antiprotons/dummy_file.txt", 4, dummy_data);

  // for (int i = 0; i < dummy_data[0].size(); i++)
  //{
  // std::cout << dummy_data[0][1] << " ";
  //}
  EXPECT_EQ(dummy_data[0].size(), 4);

  EXPECT_EQ(dummy_data[0][0], 0.000000000);
  EXPECT_EQ(dummy_data[19][0], 0.000000019);
  EXPECT_EQ(dummy_data[0][3], 0.000000000);
  EXPECT_EQ(dummy_data[19][3], 0.000019000);
}
*/
TEST(IndirectParamTest, SpectrumIntegration)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);

  std::vector<std::vector<real_t>> spectrum = {{0.1, 10.0}, {1.0, 5.0}, {10.0, 1.0}};

  real_t integral = ind_param.integrate_spectrum(spectrum, 0.5, 5.0);
  // Manually computed integral over the specified bin
  real_t expected_integral = 3.44444444444;

  EXPECT_NEAR(integral, expected_integral, 1e-6);
}

TEST(IndirectParamTest, LikelihoodInterpolation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);

  ind_param.read_fermi_data();

  real_t flux = 2.5e-10;
  int dsph = 0;
  int bin = 5;

  real_t likelihood = ind_param.interpolate_likelihood(dsph, bin, flux);

  // Manually computed expected likelihood value
  real_t expected_likelihood = -0.0017305;

  EXPECT_NEAR(likelihood, expected_likelihood, 1e-6);
}

TEST(IndirectParamTest, FermiReading)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);

  EXPECT_NO_THROW({ ind_param.read_fermi_data(); });

  auto logJ_factors = ind_param.get_logJ_factors();

  EXPECT_EQ(logJ_factors.at(0).at(0), 18.2);
  EXPECT_EQ(logJ_factors.at(18).at(0), 18.3);
  EXPECT_EQ(logJ_factors.at(44).at(0), 18.9);
}

TEST(IndirectParamTest, LogxInterpolation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);

  std::vector<std::vector<real_t>> x_data = {{1.0, 100.0}, {10.0, 10.0}, {100.0, 1.0}};

  real_t x = 5.0;
  real_t y_interp = ind_param.logx_interpol(x_data, x);

  // Manually computed expected interpolated value
  real_t expected_y = 60.;

  EXPECT_NEAR(y_interp, expected_y, 1e-6);
}

TEST(IndirectParamTest, total_sigma_v)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();
  auto total_sigma_v = ind_param.get_total_sigma_v();

  EXPECT_EQ(sigma_v_process.size(), 25);
  EXPECT_NEAR(total_sigma_v, 4.082059e-30,
              4.082059e-30 * 0.1); // Check that total sigma_v is close to the expected value within 10% tolerance
}

TEST(IndirectParamTest, dof_calculation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto dof = ind_param.get_dof();

  EXPECT_EQ(dof, 6); // Check that the number of degrees of freedom is equal to the expected value
}

//------------------------------ TESTING INDIVIDUAL CROSS SECTIONS ------------------------------------------//

TEST(IndirectParamTest, h_z_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(3), 2.05907e-31,
              2.05907e-31 * 0.1); // Check that the thermally averaged cross section for the process at index 3 is close
                                  // to the expected value within 10% tolerance
}

TEST(IndirectParamTest, ubar_u_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(17), 8.04793e-39,
              8.04793e-39 * 0.1); // Check that the thermally averaged cross section for the process at index 17 is
                                  // close to the expected value within 10% tolerance
}

TEST(IndirectParamTest, dbar_d_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(19), 7.71775e-36,
              7.71775e-36 * 0.1); // Check that the thermally averaged cross section for the process at index 19 is
                                  // close to the expected value within 10% tolerance
}

TEST(IndirectParamTest, cbar_c_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(21), 5.86478e-33,
              5.86478e-33 * 0.1); // Check that the thermally averaged cross section for the process at index 21 is
                                  // close to the expected value within 10% tolerance
}

TEST(IndirectParamTest, sbar_s_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(22), 3.02177e-33,
              3.02177e-33 * 0.1); // Check that the thermally averaged cross section for the process at index 22 is
                                  // close to the expected value within 10% tolerance
}

TEST(IndirectParamTest, bbar_b_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(20), 2.87619e-30,
              2.87619e-30 * 0.1); // Check that the thermally averaged cross section for the process at index 20 is
                                  // close to the expected value within 10% tolerance
}


TEST(IndirectParamTest, tbar_t_final_state)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();

  EXPECT_NEAR(sigma_v_process.at(23), 1.72962e-31,
              1.72962e-31 * 0.1); // Check that the thermally averaged cross section for the process at index 23 is
                                  // close to the expected value within 10% tolerance
}


TEST(IndirectParamTest, FillSpectrum)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  int DM_candidate = input.getLightestBSMpart();
  real_t mass_chi = input.masses_vector.at(DM_candidate);
  std::cout << mass_chi << "\n";
  auto sigma_v_process = ind_param.get_sigma_v_process();
  auto processes = ind_param.get_processes();

  ind_param.fill_spectrum();

  auto sigma_v_table = ind_param.get_sigma_v_table();
  auto energy_table = ind_param.get_energy_table();
  std::ofstream output_file("sigma_v_table.txt");

  if (output_file.is_open())
  {
    output_file << "Energy, ee, mumu, tautau, qq, cc, bb, tt, WW, ZZ, gg, aa, nue, numu, nutau \n";
    for (size_t i = 0; i < energy_table.size(); i++)
    {

      // char* fpp;
      // sprintf(fpp, "%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", energy_table[i],
      // sigma_v_table[0][i],
      // sigma_v_table[1][i], sigma_v_table[2][i], sigma_v_table[3][i], sigma_v_table[4][i], sigma_v_table[5][i],
      // sigma_v_table[6][i], sigma_v_table[7][i], sigma_v_table[8][i], sigma_v_table[9][i], sigma_v_table[10][i],
      // sigma_v_table[11][i], sigma_v_table[12][i], sigma_v_table[13][i]);
      // output_file << fpp << std::endl;
      output_file << energy_table[i] << ", " << sigma_v_table[0][i] << ", " << sigma_v_table[1][i] << ", "
                  << sigma_v_table[2][i] << ", " << sigma_v_table[3][i] << ", " << sigma_v_table[4][i] << ", "
                  << sigma_v_table[5][i] << ", " << sigma_v_table[6][i] << ", " << sigma_v_table[7][i] << ", "
                  << sigma_v_table[8][i] << ", " << sigma_v_table[9][i] << ", " << sigma_v_table[10][i] << ", "
                  << sigma_v_table[11][i] << ", " << sigma_v_table[12][i] << ", " << sigma_v_table[13][i] << "\n";
    }
    output_file.close();
  }
  else
  {
    std::cerr << "Unable to open file for writing.\n";
  }

  real_t dof = ind_param.get_dof();

  EXPECT_EQ(dof, 6);
}

TEST(IndirectParamTest, Fluxes)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirect_param_t ind_param(input);
  std::vector<std::vector<real_t>> prod_spec;

  ind_param.compute_fluxes(prod_spec);

  auto fermi_eflux = ind_param.get_fermi_eflux();

  std::ofstream output_file("fluxes_gamma.txt");
  for (int i = 0; i < fermi_eflux.size(); i++)
  {
    real_t bin_energy = log10(500) + (log10(500000) - log10(500)) / 24. * i;
    output_file << "Energy bin " << i << ": " << fermi_eflux[i] << "\n";
  }
  output_file.close();
}

TEST(PropagationParamTest, Initialization)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  EXPECT_NO_THROW({
    Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                   Propagation_param_t::halo_profiles::NFW);
  });
}


TEST(PropagationParamTest, FunctionOptimizationPolin)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);
  std::vector<std::vector<real_t>> dummy_spectrum;
  std::vector<real_t> xlim_min = {-1.0, -1.0};
  std::vector<real_t> xlim_max = {1.0, 1.0};
  real_t ftol = 1.e-8;
  real_t fmin;
  std::vector<real_t> x_min;
  x_min.resize(2);
  std::vector<real_t> xtra = {5.0, 4.0, 3.0, 2.0, 0.5};

  prop_param.IDpowell(2, &Propagation_param_t::test_func, xtra, dummy_spectrum, xlim_min, xlim_max, &fmin, x_min, ftol,
                      "min");

  EXPECT_NEAR(x_min[0], -5.0 / 22.0, 1e-6);
  EXPECT_NEAR(x_min[1], 3.0 / 44.0, 1e-6);
}

TEST(PropagationParamTest, FunctionOptimizationNonLin)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);
  std::vector<std::vector<real_t>> spectrum = {{1.0}};
  std::vector<real_t> xtra = {1.0};
  std::vector<real_t> xlim_min = {1.5};
  std::vector<real_t> xlim_max = {3.};
  real_t ftol = 1.e-8;
  real_t fmin;
  std::vector<real_t> x_min;
  x_min.resize(1);

  prop_param.IDpowell(1, &Propagation_param_t::test_func_spec, xtra, spectrum, xlim_min, xlim_max, &fmin, x_min, ftol,
                      "min");

  EXPECT_NEAR(x_min[0], 2., 1e-6);
}

TEST(PropagationParamTest, BesselZeros)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);


  auto alpha_i = prop_param.get_alpha_i();

  EXPECT_NEAR(alpha_i[0], 2.404826, 1e-6);
  EXPECT_NEAR(alpha_i[9], 30.634606, 1e-6);
  EXPECT_NEAR(alpha_i[99], 313.374266, 1e-6);
}

TEST(PropagationParamTest, f_PSRD)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  real_t result = prop_param.f_pulsar_distribution(0.5);

  EXPECT_NEAR(result, 6.615771e-01, 1e-6);
}

TEST(PropagationParamTest, cosmic_ray_cross_sections)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  std::vector<real_t> E_pbar_vec, E_proton_vec, sigma_in_pbarH_vec, sigma_total_pH_vec, sigma_in_pH_vec,
      sigma_in_noann_pbarH_vec;



  for (size_t i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
  {
    real_t E_pbar = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((double)i_pbar / (double)DIM_TAB_PBAR)) + PROTON_MASS;
    real_t sigma_in_noann_pbarH = prop_param.sigma_inelastic_NOANN_pbarH_TAN_and_NG(E_pbar);
    real_t sigma_in_pbarH = prop_param.sigma_inelastic_pbarH_TAN_and_NG(E_pbar);
    E_pbar_vec.push_back(E_pbar);
    sigma_in_noann_pbarH_vec.push_back(sigma_in_noann_pbarH);
    sigma_in_pbarH_vec.push_back(sigma_in_pbarH);
  }
  for (size_t i_proton = 0; i_proton <= DIM_TAB_PROTON; i_proton++)
  {
    real_t E_proton = E_PROTON_MIN * pow((E_PROTON_MAX / E_PROTON_MIN), ((real_t)i_proton / (real_t)DIM_TAB_PROTON));
    real_t sigma_in_pH = prop_param.sigma_inelastic_pH_TAN_and_NG(E_proton);
    real_t sigma_total_pH = prop_param.sigma_total_pH(E_proton);
    E_proton_vec.push_back(E_proton);
    sigma_in_pH_vec.push_back(sigma_in_pH);
    sigma_total_pH_vec.push_back(sigma_total_pH);
  }
  std::ofstream output_file_pbar("sigma_in_pbarH.txt");
  std::ofstream output_file_proton("sigma_pH.txt");
  if (output_file_pbar.is_open() && output_file_proton.is_open())
  {
    std::cout << "Files opened\n";
    char fpp[500];
    output_file_pbar << "# ENERGY  INELASTIC PBAR-H XSECTION   INELASTIC NO-ANN PBAR-H XSECTION" << std::endl;
    for (size_t i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      sprintf(fpp, "%e    %e    %e", E_pbar_vec[i_pbar], sigma_in_pbarH_vec[i_pbar], sigma_in_noann_pbarH_vec[i_pbar]);
      output_file_pbar << fpp << std::endl;
    }
    output_file_pbar.close();

    output_file_proton << "# ENERGY  TOTAL PROTON-H XSECTION   INELASTIC PROTON-H XSECTION" << std::endl;
    for (size_t i_proton = 0; i_proton <= DIM_TAB_PROTON; i_proton++)
    {
      sprintf(fpp, "%e    %e    %e", E_proton_vec[i_proton], sigma_total_pH_vec[i_proton], sigma_in_pH_vec[i_proton]);
      output_file_proton << fpp << std::endl;
    }
    output_file_proton.close();
  }
}

TEST(PropagationParamTest, primary_CR_flux)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);
  std::vector<real_t> q_i = prop_param.get_q_i();
  std::ofstream output_file("primary_cosmic_flux.txt");


  if (output_file.is_open())
  {
    char fpp[500];
    output_file << "# BESSEL ORDER     BESSEL COEFFICIENT\n";

    for (size_t j = 0; j < N_BESSEL; j++)
    {
      sprintf(fpp, "%d    %e", j + 1, q_i[j]);
      output_file << fpp << std::endl;
    }

    output_file.close();
  }
}


TEST(PropagationParamTest, proton_helium_bessel_coeffs)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  std::vector<real_t> E_nuc_vec;
  for (int i_nuc = 0; i_nuc <= DIM_TAB_PROTON; i_nuc++)
  {
    real_t E_nuc = E_PROTON_MIN * pow((E_PROTON_MAX / E_PROTON_MIN), ((real_t)i_nuc / (real_t)DIM_TAB_PROTON));
    E_nuc_vec.push_back(E_nuc);
  }
  auto bessel_coef_proton = prop_param.get_bessel_coef_proton();
  auto bessel_coef_helium = prop_param.get_bessel_coef_helium();

  std::ofstream output_file_proton("bessel_coef_p.txt");
  std::ofstream output_file_helium("bessel_coef_he.txt");
  if (output_file_proton.is_open() && output_file_helium.is_open())
  {
    char fpp[500];
    output_file_proton << "# ENERGY     BESSEL ORDER     BESSEL COEFFICIENT\n";
    for (size_t i = 0; i < bessel_coef_proton.size(); i++)
    {
      for (size_t j = 0; j < N_BESSEL; j++)
      {
        sprintf(fpp, "%e    %d    %e", E_nuc_vec[i], j + 1, bessel_coef_proton[i][j]);
        output_file_proton << fpp << std::endl;
      }
    }
    output_file_proton.close();

    output_file_helium << "# ENERGY     BESSEL ORDER     BESSEL COEFFICIENT\n";
    for (size_t i = 0; i < bessel_coef_proton.size(); i++)
    {
      for (size_t j = 0; j <= N_BESSEL; j++)
      {
        sprintf(fpp, "%e    %d    %e", E_nuc_vec[i], j, bessel_coef_helium[i][j]);
        output_file_helium << fpp << std::endl;
      }
    }
    output_file_helium.close();
  }
}

TEST(PropagationParamTest, secondary_spectra_uncertainties)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::BURKERT);

  prop_param.preliminary_secondary_spectrum_IS_calculation();

  auto xs_uncertainties = prop_param.get_xs_uncertainties();
  std::ofstream output_file("secondary_spec_xs_uncertainties_low.txt");
  std::ofstream output_file_high("secondary_spec_xs_uncertainties_high.txt");
  if (output_file.is_open())
  {
    char fpp[500];
    char fpp_high[500];
    for (size_t i = 0; i < xs_uncertainties.size(); i++)
    {
      sprintf(fpp, "%e %e", xs_uncertainties[i][0], xs_uncertainties[i][1]);
      sprintf(fpp_high, "%e %e", xs_uncertainties[i][0], xs_uncertainties[i][2]);
      output_file << fpp << std::endl;
      output_file_high << fpp_high << std::endl;
    }
    output_file.close();
    output_file_high.close();
  }
}

TEST(PropagationParamTest, primary_source_calculation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  real_t DM_mass = input.getLightestBSMmass();
  std::vector<real_t> dNpbar_on_dEpbar_vec;
  std::vector<real_t> primary_source_term;
  prop_param.primary_source_calculation(DM_mass, dNpbar_on_dEpbar_vec, primary_source_term);


  std::ofstream output_file("primary_source_output.txt");
  if (output_file.is_open())
  {
    output_file << "# ENERGY   DNPBAR_ON_DEPBAR  PRIMARY_SOURCE_TERM" << std::endl;
    for (size_t i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      real_t T_pbar_IS = T_PBAR_MIN * pow((T_PBAR_MAX / T_PBAR_MIN), ((double)i_pbar / (double)DIM_TAB_PBAR));
      char fpp[500];
      sprintf(fpp, "%e  %e  %e", T_pbar_IS, dNpbar_on_dEpbar_vec[i_pbar], primary_source_term[i_pbar]);
      output_file << fpp << std::endl;
    }
    output_file.close();
  }
}

TEST(PropagationParamTest, primary_flux_calculation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  std::vector<real_t> dNpbar_on_dEpbar_vec;
  std::vector<real_t> primary_source_term;
  std::vector<std::vector<real_t>> pbar_spectrum;
  prop_param.primary_spectra_BCGS_2014(dNpbar_on_dEpbar_vec, primary_source_term, pbar_spectrum);

  std::ofstream output_file("primary_flux_output.txt");
  if (output_file.is_open())
  {
    output_file << "# Log10(E)   PRIMARY CR FLUX" << std::endl;
    for (size_t i_pbar = 0; i_pbar <= DIM_TAB_PBAR; i_pbar++)
    {
      char fpp[500];
      sprintf(fpp, "%e  %e", pbar_spectrum[i_pbar][0], pbar_spectrum[i_pbar][1]);
      output_file << fpp << std::endl;
    }
    output_file.close();
  }
}

TEST(PropagationParamTest, solar_modulation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  std::vector<real_t> logE;
  logE.resize(N_AMS02 + 1);
  for (int i = 0; i < N_AMS02; i++)
    logE.at(i) = log10(AMS02.at(i).at(0));
  logE[N_AMS02] = 2.6;

  std::vector<real_t> param = {0.5, 0.5};
  std::vector<std::vector<real_t>> back, backS, tot, totS, spect;
  prop_param.background_spectAMS(logE, param[0], back);
  prop_param.solar_mod(back, param[1], backS);
  std::ofstream output_file("solar_mod_test.txt");
  if (output_file.is_open())
  {
    output_file << "# Log10(E)   SOLAR_MOD_BKG" << std::endl;
    for (size_t i = 0; i < backS.size() - 1; i++)
    {
      char fpp[500];
      sprintf(fpp, "%e  %e", backS[i][0], backS[i][1]);
      output_file << fpp << std::endl;
    }
    output_file.close();
  }

  std::vector<real_t> dEpbar_on_dNpbar_temp, primary_source_temp;
  prop_param.primary_spectra_BCGS_2014(dEpbar_on_dNpbar_temp, primary_source_temp, spect);
  prop_param.add_spectra(back, spect, tot);
  prop_param.solar_mod(tot, param[1], totS);
  std::ofstream output_file_tot("solar_mod_total.txt");
  if (output_file_tot.is_open())
  {
    output_file_tot << "# Log10(E)   SOLAR_MOD_TOT" << std::endl;
    for (size_t i = 0; i < backS.size() - 1; i++)
    {
      char fpp[500];
      sprintf(fpp, "%e  %e", backS[i][0], backS[i][1]);
      output_file_tot << fpp << std::endl;
    }
    output_file_tot.close();
  }
}

TEST(PropagationParamTest, chi_squared_bg)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW);

  int IDpowell;
  std::vector<real_t> final_params = prop_param.chi2_bkg(IDpowell);
  real_t chi2_noDM = prop_param.get_chi2_noDM();
  std::vector<real_t> param = {1., 0.948228};
  std::cout << "chi2_noDM = " << chi2_noDM << std::endl;
  std::cout << "Final phi_f = " << final_params[1] << " Final A = " << final_params[0] << std::endl;


  std::vector<real_t> logE;
  logE.resize(N_AMS02 + 1);
  for (int i = 0; i < N_AMS02; i++)
    logE.at(i) = log10(AMS02.at(i).at(0));
  logE[N_AMS02] = 2.6;
  std::vector<std::vector<real_t>> spec;
}


TEST(PropagationParamTest, Fermi_Likelihood)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::BURKERT);

  int powell;
  real_t delta_likelihood = prop_param.deltalikelihood_fermi(powell);
  std::cout << "Delta likelihood = " << delta_likelihood << std::endl;
}

/*
class OstreamListener : public testing::EmptyTestEventListener
{
  public:
  explicit OstreamListener(std::ostream& os) : os_(os) {}

  void OnTestEnd(const testing::TestInfo& test_info) override
  {
    os_ << "[TEST] " << test_info.test_suite_name() << "." << test_info.name() << std::endl;
  }

  private:
  std::ostream& os_;
};
*/

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

  /*
  std::ofstream file("test_output.txt");
  auto& listeners = testing::UnitTest::GetInstance()->listeners();
  listeners.Append(new OstreamListener(file));
  */

  return RUN_ALL_TESTS();
}