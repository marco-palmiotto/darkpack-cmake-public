#pragma once
#ifndef __cplusplus
  #define __cplusplus
#endif

#include <gtest/gtest.h>

#include "dp_scalar2to2/RunningSM.hpp"
#include "dp_scalar2to2/indirect.hpp"
#include "dp_scalar2to2/process.hpp"
#include "dp_scalar2to2/process_1to2.hpp"
#include "dp_scalar2to2/propagation.hpp"
#include "dp_scalar2to2/relicparam.hpp"

using namespace scalar2to2;
using namespace advmath;

TEST(IndirectParamTest, Initialization)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/scalar.lha");
  EXPECT_NO_THROW({ Indirectparam_t ind_param(input); });
}

TEST(IndirectParamTest, ReadingFile)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/scalar.lha");
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
  struct Param_t input("/workspaces/darkpack-cmake/src/models/scalar.lha");
  Indirectparam_t ind_param(input);

  std::vector<std::vector<real_t>> spectrum = {{0.1, 10.0}, {1.0, 5.0}, {10.0, 1.0}};

  real_t integral = ind_param.integrate_spectrum(spectrum, 0.5, 5.0);
  // Manually computed integral over the specified bin
  real_t expected_integral = 3.44444444444;

  EXPECT_NEAR(integral, expected_integral, 1e-6);
}

TEST(PropagationParamTest, Initialization)
{
  struct Param_t input("/workspaces/darkpack-cmake/src/models/scalar.lha");
  EXPECT_NO_THROW({
    Propagation_param_t prop_param(input, Propagation_param_t::propagation_models::MED,
                                   Propagation_param_t::halo_profiles::NFW, 8.5, 0.3);
  });
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}