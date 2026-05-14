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

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "The program requires one argument:" << std::endl << "1. The path to the model SLHA file" << std::endl;
    return 1;
  }

  struct Param_t input(argv[1]);
  Propagation_param_t prop_param_std(input, Propagation_param_t::propagation_models::MED,
                                     Propagation_param_t::halo_profiles::EINASTO);
  Propagation_param_t prop_param_str(input, Propagation_param_t::propagation_models::MAX,
                                     Propagation_param_t::halo_profiles::EINASTO);
  Propagation_param_t prop_param_con(input, Propagation_param_t::propagation_models::MED,
                                     Propagation_param_t::halo_profiles::BURKERT);
  // Indirect_param_t ind_param(input);
  // ind_param.read_fermi_data();
  // ind_param.fetch_pppc4dmid_data();
  std::cout << std::scientific << std::setprecision(6);
  /*
  for (int i = 0; i < 45; i++)
  {
    real_t interpolated = ind_param.interpolate_likelihood(i, 0, 0);
    std::cout << "dSph number " << i << " interpolated no DM flux " << interpolated << std::endl;
  }
  */
  // prop_param_con.change_ind_param(0., 10.);

  std::ofstream con_file("spect_conservative_bb_m_10GeV.txt");
  con_file << std::scientific << std::setprecision(6);

  std::ofstream std_file("spect_standard_bb_m_10GeV.txt");
  std_file << std::scientific << std::setprecision(6);

  std::ofstream str_file("spect_stringent_bb_m_10GeV.txt");
  str_file << std::scientific << std::setprecision(6);


  std::vector<real_t> optimal_params_std, optimal_params_str, optimal_params_con;
  int DM_candidate = input.getLightestBSMpart();
  real_t DM_mass = input.masses_vector.at(DM_candidate);

  int powell;
  prop_param_con.change_ind_param(0., 10);
  prop_param_std.change_ind_param(0., 10);
  prop_param_str.change_ind_param(0., 10);
  prop_param_con.single_out_annihilation_channel(corr::b);
  prop_param_std.single_out_annihilation_channel(corr::b);
  prop_param_str.single_out_annihilation_channel(corr::b);

  prop_param_con.change_ind_param(1.39099e-28, 0.);
  real_t chi2_con = -prop_param_con.deltachi2_AMS(powell, optimal_params_con) + chi2_table.at(4).at(6);
  real_t chi2_std = prop_param_std.deltachi2_AMS(powell, optimal_params_std) + chi2_table.at(4).at(6);
  real_t chi2_str = prop_param_str.deltachi2_AMS(powell, optimal_params_str) + chi2_table.at(4).at(6);
  real_t chi2_noDM = prop_param_con.get_chi2_noDM();
  std::cout << "chi2_noDM for conservative: " << chi2_noDM << std::endl;
  std::cout << chi2_con << "\t" << chi2_std << "\t" << chi2_str << std::endl;
  for (int k = 0; k <= 80; ++k)
  {
    real_t sv = 1.24e-29;  // 1e-30 to 1e-24
    real_t sv2 = 1.25e-29; // 1e-30 to 1e-24

    prop_param_con.change_ind_param(sv, 0.);

    int status;
    std::vector<real_t> opt, opt2;

    real_t f = -prop_param_con.deltachi2_AMS(status, opt) + chi2_table.at(4).at(1);
    prop_param_con.change_ind_param(sv2, 0.);
    f = -prop_param_con.deltachi2_AMS(status, opt2) + chi2_table.at(4).at(1);
    std::cout << "A = " << opt[0] << " phi_f = " << opt[1] << "\n";
    std::cout << "A = " << opt2[0] << " phi_f = " << opt2[1] << "\n";
  }
  for (int i = 0; i < 2; i++)
  {
    std::cout << "conservative: " << optimal_params_con[i] << std::endl;
    std::cout << "standard: " << optimal_params_std[i] << std::endl;
    std::cout << "stringent: " << optimal_params_str[i] << std::endl;
  }
  con_file << "# Pbar spectrum with optimal parameters\n";
  con_file.flush();
  if (con_file)
  {
    std::vector<real_t> dEpbar_on_dNpbar_temp, primary_source_temp;
    std::vector<std::vector<real_t>> spec;
    prop_param_con.primary_spectra_BCGS_2014(dEpbar_on_dNpbar_temp, primary_source_temp, spec);

    for (int i = 0; i < spec.size(); i++)
      std::cout << spec.at(i).at(0) << "\t" << spec.at(i).at(1) << std::endl;

    std::vector<real_t> logE;
    logE.resize(N_AMS02 + 1);
    for (int i = 0; i < N_AMS02; i++)
      logE.at(i) = log10(AMS02.at(i).at(0));
    logE[N_AMS02] = 2.6;

    std::vector<std::vector<real_t>> tot, back, totS, specS, backS;
    real_t A = optimal_params_con.at(0);
    real_t phi_f = optimal_params_con.at(1);
    prop_param_con.background_spectAMS(logE, A, back);
    prop_param_con.add_spectra(back, spec, tot);
    prop_param_con.solar_mod(tot, phi_f, totS);

    for (int i = 0; i < totS.size(); i++)
    {
      con_file << totS.at(i).at(0) << "\t" << totS.at(i).at(1) << std::endl;
    }
    prop_param_con.solar_mod(spec, phi_f, specS);
    prop_param_con.solar_mod(back, phi_f, backS);
    con_file << "# Primary pbar spectrum\n";
    for (int i = 0; i < specS.size(); i++)
      con_file << specS.at(i).at(0) << "\t" << specS.at(i).at(1) << std::endl;

    con_file << "# Background pbar spectrum\n";
    for (int i = 0; i < backS.size(); i++)
      con_file << backS.at(i).at(0) << "\t" << backS.at(i).at(1) << std::endl;
  }
  con_file.close();

  std_file << "# Pbar spectrum with optimal parameters\n";
  std_file.flush();
  if (std_file)
  {
    std::vector<real_t> dEpbar_on_dNpbar_temp, primary_source_temp;
    std::vector<std::vector<real_t>> spec;
    prop_param_std.primary_spectra_BCGS_2014(dEpbar_on_dNpbar_temp, primary_source_temp, spec);

    std::vector<real_t> logE;
    logE.resize(N_AMS02 + 1);
    for (int i = 0; i < N_AMS02; i++)
      logE.at(i) = log10(AMS02.at(i).at(0));
    logE[N_AMS02] = 2.6;

    std::vector<std::vector<real_t>> tot, back, totS, backS, specS;
    real_t A = optimal_params_std.at(0);
    real_t phi_f = optimal_params_std.at(1);
    prop_param_std.background_spectAMS(logE, A, back);
    prop_param_std.add_spectra(back, spec, tot);
    prop_param_std.solar_mod(tot, phi_f, totS);

    for (int i = 0; i < totS.size(); i++)
    {
      std_file << totS.at(i).at(0) << "\t" << totS.at(i).at(1) << std::endl;
    }
    prop_param_std.solar_mod(spec, phi_f, specS);
    prop_param_std.solar_mod(back, phi_f, backS);
    std_file << "# Primary pbar spectrum\n";
    for (int i = 0; i < specS.size(); i++)
      std_file << specS.at(i).at(0) << "\t" << specS.at(i).at(1) << std::endl;

    std_file << "# Background pbar spectrum\n";
    for (int i = 0; i < backS.size(); i++)
      std_file << backS.at(i).at(0) << "\t" << backS.at(i).at(1) << std::endl;
  }
  std_file.close();


  str_file << "# Pbar spectrum with optimal parameters\n";
  str_file.flush();
  if (str_file)
  {
    std::vector<real_t> dEpbar_on_dNpbar_temp, primary_source_temp;
    std::vector<std::vector<real_t>> spec;
    prop_param_str.primary_spectra_BCGS_2014(dEpbar_on_dNpbar_temp, primary_source_temp, spec);

    std::vector<real_t> logE;
    logE.resize(N_AMS02 + 1);
    for (int i = 0; i < N_AMS02; i++)
      logE.at(i) = log10(AMS02.at(i).at(0));
    logE[N_AMS02] = 2.6;

    std::vector<std::vector<real_t>> tot, back, totS, specS, backS;
    real_t A = optimal_params_str.at(0);
    real_t phi_f = optimal_params_str.at(1);
    prop_param_str.background_spectAMS(logE, A, back);
    prop_param_str.add_spectra(back, spec, tot);
    prop_param_str.solar_mod(tot, phi_f, totS);

    for (int i = 0; i < totS.size(); i++)
    {
      str_file << totS.at(i).at(0) << "\t" << totS.at(i).at(1) << std::endl;
    }
    prop_param_str.solar_mod(spec, phi_f, specS);
    prop_param_str.solar_mod(back, phi_f, backS);
    str_file << "# Primary pbar spectrum\n";
    for (int i = 0; i < specS.size(); i++)
      str_file << specS.at(i).at(0) << "\t" << specS.at(i).at(1) << std::endl;

    str_file << "# Background pbar spectrum\n";
    for (int i = 0; i < backS.size(); i++)
      str_file << backS.at(i).at(0) << "\t" << backS.at(i).at(1) << std::endl;
  }
  str_file.close();

  return 0;
}
