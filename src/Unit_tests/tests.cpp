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

using namespace mssm2to2;
using namespace advmath;

TEST(IndirectParamTest, Initialization)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  EXPECT_NO_THROW({ Indirectparam_t ind_param(input); });
}

TEST(IndirectParamTest, ReadingFile)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirectparam_t ind_param(input);

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

TEST(IndirectParamTest, SpectrumIntegration)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirectparam_t ind_param(input);

  std::vector<std::vector<real_t>> spectrum = {{0.1, 10.0}, {1.0, 5.0}, {10.0, 1.0}};

  real_t integral = ind_param.integrate_spectrum(spectrum, 0.5, 5.0);
  // Manually computed integral over the specified bin
  real_t expected_integral = 3.44444444444;

  EXPECT_NEAR(integral, expected_integral, 1e-6);
}

TEST(IndirectParamTest, LikelihoodInterpolation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirectparam_t ind_param(input);

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
  Indirectparam_t ind_param(input);

  EXPECT_NO_THROW({ ind_param.read_fermi_data(); });

  auto logJ_factors = ind_param.get_logJ();

  EXPECT_EQ(logJ_factors.at(0).at(0), 18.2);
  EXPECT_EQ(logJ_factors.at(18).at(0), 18.3);
  EXPECT_EQ(logJ_factors.at(44).at(0), 18.9);
}

TEST(IndirectParamTest, LogxInterpolation)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirectparam_t ind_param(input);

  std::vector<std::vector<real_t>> x_data = {{1.0, 100.0}, {10.0, 10.0}, {100.0, 1.0}};

  real_t x = 5.0;
  real_t y_interp = ind_param.logx_interpol(x_data, x);

  // Manually computed expected interpolated value
  real_t expected_y = 60.;

  EXPECT_NEAR(y_interp, expected_y, 1e-6);
}

TEST(IndirectParamTest, sigma_v)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  std::cout << input.getLightestBSMpart() << "\n";
  Indirectparam_t ind_param(input);
  auto sigma_v_process = ind_param.get_sigma_v_process();
  auto total_sigma_v = ind_param.get_total_sigma_v();

  EXPECT_EQ(sigma_v_process.size(), 25);
  EXPECT_NEAR(total_sigma_v, 4.082059e-30, 1e-32);
}

TEST(IndirectParamTest, FillSpectrum)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Indirectparam_t ind_param(input);
  int DM_candidate = input.getLightestBSMpart();
  real_t mass_chi = input.masses_vector.at(DM_candidate);
  std::cout << mass_chi << "\n";
  auto sigma_v_process = ind_param.get_sigma_v_process();
  auto processes = ind_param.get_processes();

  ind_param.fill_spectrum(processes, sigma_v_process, 2. * mass_chi);

  real_t dof = ind_param.get_dof();

  EXPECT_EQ(dof, 6);
}

TEST(PropagationParamTest, Initialization)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  EXPECT_NO_THROW({
    Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                   Propagation_param_t::halo_profiles::NFW, 8.5, 0.3);
  });
}


TEST(PropagationParamTest, FunctionOptimizationPolin)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/mssm.lha");
  Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                 Propagation_param_t::halo_profiles::NFW, 8.5, 0.3);
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
                                 Propagation_param_t::halo_profiles::NFW, 8.5, 0.3);
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
                                 Propagation_param_t::halo_profiles::NFW, 8.5, 0.3);

  prop_param.search_zeros_J0(1e-8);

  auto alpha_i = prop_param.get_alpha_i();

  EXPECT_NEAR(alpha_i[0], 2.404826, 1e-6);
  EXPECT_NEAR(alpha_i[9], 30.634606, 1e-6);
  EXPECT_NEAR(alpha_i[99], 313.374266, 1e-6);
}



int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}