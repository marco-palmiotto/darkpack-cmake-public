#pragma once
#ifndef __cplusplus
  #define __cplusplus
#endif

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <optional>
#include <vector>

#include "dp_mssm2to2/RunningSM.hpp"
#include "dp_mssm2to2/indirect.hpp"
#include "dp_mssm2to2/process.hpp"
#include "dp_mssm2to2/process_1to2.hpp"
#include "dp_mssm2to2/propagation.hpp"
#include "dp_mssm2to2/relicparam.hpp"

using namespace mssm2to2;
using namespace advmath;

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    std::cerr << "The program requires three arguments:\n"
              << "1. The path to the model SLHA file\n"
              << "2. the desired propagation model: 0=nominal, 1=conservative, 2=inclusive\n"
              << "3. the desired particle channel: 0=bbar, 1=ttbar, 2=WW, 3=ZZ\n";
    return 1;
  }

  Param_t input(argv[1]);

  int prop_choice = std::stoi(argv[2]);
  int channel_choice = std::stoi(argv[3]);

  std::string filename = "chi2_AMS_";
  std::optional<Propagation_param_t> prop_param;

  switch (prop_choice)
  {
  case 0:
    prop_param.emplace(input, Propagation_param_t::propagation_models::MED,
                       Propagation_param_t::halo_profiles::EINASTO);
    filename += "nominal_";
    break;

  case 1:
    prop_param.emplace(input, Propagation_param_t::propagation_models::MED,
                       Propagation_param_t::halo_profiles::BURKERT);
    filename += "conservative_";
    break;

  default:
    prop_param.emplace(input, Propagation_param_t::propagation_models::MAX,
                       Propagation_param_t::halo_profiles::EINASTO);
    filename += "inclusive_";
    break;
  }

  real_t mass_chi_min;
  real_t mass_chi_max = 10000;


  switch (channel_choice)
  {
  case 0:
    mass_chi_min = 10.;
    filename += "bbar.txt";
    break;

  case 1:
    mass_chi_min = input.masses_vector.at(corr::t);
    filename += "ttbar.txt";
    break;

  case 2:
    mass_chi_min = input.masses_vector.at(corr::W);
    filename += "WW.txt";
    break;

  default:
    mass_chi_min = input.masses_vector.at(corr::Z);
    filename += "ZZ.txt";
    break;
  }

  std::ofstream chi2_file(filename);
  if (!chi2_file)
  {
    std::cerr << "Unable to open output file: " << filename << std::endl;
    return 1;
  }

  std::vector<real_t> optimal_params;
  int powell_min, powell_max;
  real_t chi2_min, chi2_max;

  chi2_file << "# Mass   Total sigma_v (cm^3/s) (95 C.L.)     Delta chi2 AMS02\n";
  chi2_file.flush();

  for (int j = 0; j < 151; j++)
  {
    real_t new_DM_mass = mass_chi_min * std::pow((mass_chi_max / mass_chi_min), ((real_t)(j - 1) / 150.));

    std::cout << "Testing mass: " << new_DM_mass << " GeV" << std::endl;
    std::cout << j << std::endl;

    real_t xnew = 0.;
    real_t chi2_new = 0.;

    if (j != 0)
    {
      prop_param->change_ind_param(0., new_DM_mass);
      switch (channel_choice)
      {
      case 0:
        prop_param->single_out_annihilation_channel(corr::b);
        break;

      case 1:
        prop_param->single_out_annihilation_channel(corr::t);
        break;
      case 2:
        prop_param->single_out_annihilation_channel(corr::W);
        break;
      case 3:
      default:
        prop_param->single_out_annihilation_channel(corr::Z);
        break;
      }
    }
    else
    {
      continue;
    }

    real_t total_sigma_v = prop_param->get_total_sigma_v();

    std::array<std::vector<real_t>, 14> sigma_v_table = prop_param->get_sigma_v_table();
    std::vector<real_t> energy_table = prop_param->get_energy_table();

    std::cout << "Total sigma_v = " << total_sigma_v << std::endl;

    real_t sigma_v_min = total_sigma_v;

    int dof = prop_param->get_dof();

    std::cout << "Degrees of freedom: " << dof << std::endl;
    std::cout << chi2_table.at(4).at(dof) << std::endl;

    chi2_min = -prop_param->deltachi2_AMS(powell_min, optimal_params) + chi2_table.at(4).at(dof);

    int test2 = 0;

    while (powell_min == 0)
    {
      std::cout << "Recalculating chi2_min" << std::endl;

      sigma_v_min *= .9;
      prop_param->change_ind_param(sigma_v_min, 0.);

      chi2_min = -prop_param->deltachi2_AMS(powell_min, optimal_params) + chi2_table.at(4).at(dof);

      std::cout << "New chi_2 min: " << chi2_min << std::endl;

      test2++;
      if (test2 == 10)
        break;
    }

    if (test2 == 10)
    {
      std::cerr << "Unable to find a converging chi_2 after 10 attempts. Skipping this mass point." << std::endl;
      continue;
    }

    real_t sigma_v_max = sigma_v_min * 100.;

    prop_param->change_ind_param(sigma_v_max, 0.);

    std::cout << "New sigma_v = " << prop_param->get_total_sigma_v() << std::endl;

    chi2_max = -prop_param->deltachi2_AMS(powell_max, optimal_params) + chi2_table.at(4).at(dof);

    std::cout << "Chi2 min: " << chi2_min << ", Chi2 max: " << chi2_max << std::endl;

    int test = 0;

    while (powell_max == 0)
    {
      std::cout << "Recalculating chi2_max" << std::endl;

      sigma_v_max *= .9;
      prop_param->change_ind_param(sigma_v_max, 0.);

      chi2_max = -prop_param->deltachi2_AMS(powell_max, optimal_params) + chi2_table.at(4).at(dof);

      std::cout << "New chi_2 max: " << chi2_max << std::endl;

      test2++;
      if (test2 == 10)
        break;
    }

    if (test2 == 10)
    {
      std::cerr << "Unable to find a converging chi_2 after 10 attempts. Skipping this mass point." << std::endl;
      continue;
    }

    while (chi2_min * chi2_max > 0. && test < 20)
    {
      test2 = 0;

      std::cout << "Expanding the range of sigma_v. Attempt: " << test << std::endl;

      chi2_min = chi2_max;
      sigma_v_min = sigma_v_max;

      sigma_v_max *= 2.;
      prop_param->change_ind_param(sigma_v_max, 0.);

      chi2_max = -prop_param->deltachi2_AMS(powell_max, optimal_params) + chi2_table.at(4).at(dof);

      while (powell_max == 0)
      {
        sigma_v_max *= .9;
        prop_param->change_ind_param(sigma_v_max, 0.);

        chi2_max = -prop_param->deltachi2_AMS(powell_max, optimal_params) + chi2_table.at(4).at(dof);

        test2++;
        if (test2 == 7)
          break;
      }

      test++;

      if (test2 == 7)
        break;
    }

    if (test2 == 7)
    {
      std::cerr << "Unable to find a converging chi_2 after 7 attempts. Skipping this mass point." << std::endl;
      continue;
    }

    if (test == 20)
    {
      std::cerr << "Unable to find a valid range for sigma_v after 20 attempts. Skipping this mass point." << std::endl;
      continue;
    }

    const int MAXITER = 50;
    int iter = 0;
    int count = 0;
    real_t xold, chi2_old;
    while (iter < MAXITER)
    {
      std::cout << iter << std::endl;

      xnew = sigma_v_max - chi2_max * (sigma_v_max - sigma_v_min) / (chi2_max - chi2_min);

      if (xnew <= 5.e-28 && channel_choice != 0)
        xnew *= 100.;

      real_t x_init_iter = xnew;

      test2 = 0;

      prop_param->change_ind_param(xnew, 0.);

      chi2_new = -prop_param->deltachi2_AMS(powell_min, optimal_params) + chi2_table.at(4).at(dof);

      if (iter != 0 && chi2_new * chi2_old < 0 && (std::abs(xnew - xold)) / std::abs(xold) < 0.00000001 &&
          std::abs(chi2_new) > 1.)
      {
        std::cout << "Interval collapsed: must be a root\n";
        break;
      }

      real_t chi2_init_iter = chi2_new;

      while (powell_min == 0)
      {
        xnew *= .9;

        if (xnew <= sigma_v_min)
          break;

        prop_param->change_ind_param(xnew, 0.);

        chi2_new = -prop_param->deltachi2_AMS(powell_min, optimal_params) + chi2_table.at(4).at(dof);

        printf("Sigma_v: %e, Chi2: %e\n", xnew, chi2_new);

        test2++;
      }

      if (std::abs(chi2_new) > 100.)
        count++;
      else
        count = 0;

      std::cout << "count = " << count << std::endl;

      if (count == 3)
      {
        std::cout << "Updating xnew to try and solve the problem" << std::endl;

        count = 0;

        if (chi2_new < 0.)
          xnew *= 0.5;

        if (chi2_new > 0.)
          xnew *= 2.;

        prop_param->change_ind_param(xnew, 0.);

        chi2_new = -prop_param->deltachi2_AMS(powell_min, optimal_params) + chi2_table.at(4).at(dof);
      }

      if (xnew <= sigma_v_min)
      {
        xnew = x_init_iter;
        chi2_new = chi2_init_iter;
      }

      printf("Sigma_v: %e, Chi2: %e\n", xnew, chi2_new);

      if (std::abs(chi2_new) < 5.e-5)
        break;

      if (chi2_min * chi2_new < 0.)
      {
        sigma_v_max = xnew;
        chi2_max = chi2_new;
      }
      else
      {
        sigma_v_min = xnew;
        chi2_min = chi2_new;
      }
      chi2_old = chi2_new;
      xold = xnew;
      iter++;
    }

    if (test2 == 10)
    {
      std::cerr << "Unable to find a converging chi_2 after 10 attempts. Skipping this mass point." << std::endl;
      continue;
    }

    if (iter == MAXITER)
      continue;

    chi2_file << new_DM_mass << " " << xnew << " " << chi2_new << "\n";
    chi2_file.flush();
  }

  chi2_file.close();

  return 0;
}