#include "func_mssm2to2.h"
#include "clooptools.h"
#include "common.h"
#include "marty/core/looptools_init.h"
#include "stdbool.h"
#include <cmath>

#include "clib_mssm2to2.h"
#include "cparams.h"
#include "libcomplexop.h"
#include "params.h"
#include <complex.h>

namespace mssm2to2
{

  complex_t m_Hp(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_A0 = param.m_A0;
    auto res = c_m_Hp(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_h(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.m_A0 = param.m_A0;
    auto res = c_m_h(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_H0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.m_A0 = param.m_A0;
    auto res = c_m_H0(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_sc_L(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_c = param.m_c;
    cparam.beta = param.beta;
    cparam.M_q1L = param.M_q1L;
    cparam.theta_W = param.theta_W;
    auto res = c_m_sc_L(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_su_L(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.M_q1L = param.M_q1L;
    cparam.theta_W = param.theta_W;
    auto res = c_m_su_L(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_sc_R(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_c = param.m_c;
    cparam.beta = param.beta;
    cparam.M_quR = param.M_quR;
    cparam.theta_W = param.theta_W;
    auto res = c_m_sc_R(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_su_R(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.M_quR = param.M_quR;
    cparam.theta_W = param.theta_W;
    auto res = c_m_su_R(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_sd_R(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_d = param.m_d;
    cparam.V_ud = param.V_ud;
    cparam.beta = param.beta;
    cparam.M_qbR = param.M_qbR;
    cparam.M_qdR = param.M_qdR;
    cparam.theta_W = param.theta_W;
    cparam.V_cd = param.V_cd.get().real() + _mty_I * param.V_cd.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    auto res = c_m_sd_R(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_ss_R(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_s = param.m_s;
    cparam.V_us = param.V_us;
    cparam.beta = param.beta;
    cparam.M_qbR = param.M_qbR;
    cparam.M_qdR = param.M_qdR;
    cparam.theta_W = param.theta_W;
    cparam.V_cs = param.V_cs.get().real() + _mty_I * param.V_cs.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    auto res = c_m_ss_R(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_sd_L(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_d = param.m_d;
    cparam.V_ud = param.V_ud;
    cparam.beta = param.beta;
    cparam.M_q1L = param.M_q1L;
    cparam.M_q3L = param.M_q3L;
    cparam.theta_W = param.theta_W;
    cparam.V_cd = param.V_cd.get().real() + _mty_I * param.V_cd.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    auto res = c_m_sd_L(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_ss_L(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_s = param.m_s;
    cparam.V_us = param.V_us;
    cparam.beta = param.beta;
    cparam.M_q1L = param.M_q1L;
    cparam.M_q3L = param.M_q3L;
    cparam.theta_W = param.theta_W;
    cparam.V_cs = param.V_cs.get().real() + _mty_I * param.V_cs.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    auto res = c_m_ss_L(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_se_R(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_e = param.m_e;
    cparam.M_eR = param.M_eR;
    cparam.beta = param.beta;
    cparam.theta_W = param.theta_W;
    auto res = c_m_se_R(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_smu_R(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.M_eR = param.M_eR;
    cparam.beta = param.beta;
    cparam.m_mu = param.m_mu;
    cparam.theta_W = param.theta_W;
    auto res = c_m_smu_R(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_se_L(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_e = param.m_e;
    cparam.M_eL = param.M_eL;
    cparam.beta = param.beta;
    cparam.theta_W = param.theta_W;
    auto res = c_m_se_L(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_smu_L(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.M_eL = param.M_eL;
    cparam.beta = param.beta;
    cparam.m_mu = param.m_mu;
    cparam.theta_W = param.theta_W;
    auto res = c_m_smu_L(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_snu_e(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.M_eL = param.M_eL;
    cparam.beta = param.beta;
    cparam.theta_W = param.theta_W;
    auto res = c_m_snu_e(&cparam);
    return {res.real, res.imag};
  }

  complex_t m_snu_tau(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.M_tauL = param.M_tauL;
    cparam.theta_W = param.theta_W;
    auto res = c_m_snu_tau(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_e_e(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_se_L = param.m_se_L;
    cparam.m_se_R = param.m_se_R;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_el = param.Gamma_el;
    cparam.Gamma_er = param.Gamma_er;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_e_e(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_nu_e_nu_e(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_nuel = param.Gamma_nuel;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_nu_e_nu_e(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_mu_mu(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.m_smu_L = param.m_smu_L;
    cparam.m_smu_R = param.m_smu_R;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_ml = param.Gamma_ml;
    cparam.Gamma_mr = param.Gamma_mr;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_nu_mu_nu_mu(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_numl = param.Gamma_numl;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_nu_mu_nu_mu(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_tau_tau(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_tau = param.m_tau;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_tau1 = param.Gamma_tau1;
    cparam.Gamma_tau2 = param.Gamma_tau2;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    cparam.U_stau_10 = param.U_stau_10.get().real() + _mty_I * param.U_stau_10.get().imag();
    cparam.U_stau_11 = param.U_stau_11.get().real() + _mty_I * param.U_stau_11.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_nu_tau_nu_tau(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.m_snu_tau = param.m_snu_tau;
    cparam.Gamma_nutaul = param.Gamma_nutaul;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_nu_tau_nu_tau(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_u_u(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_su_L = param.m_su_L;
    cparam.m_su_R = param.m_su_R;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_ul = param.Gamma_ul;
    cparam.Gamma_ur = param.Gamma_ur;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_u_u(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_d_d(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_d = param.m_d;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_dl = param.Gamma_dl;
    cparam.Gamma_dr = param.Gamma_dr;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_d_d(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_s_s(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_sl = param.Gamma_sl;
    cparam.Gamma_sr = param.Gamma_sr;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_s_s(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_c_c(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_c = param.m_c;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sc_R = param.m_sc_R;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_cl = param.Gamma_cl;
    cparam.Gamma_cr = param.Gamma_cr;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_c_c(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_b_b(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_b1 = param.Gamma_b1;
    cparam.Gamma_b2 = param.Gamma_b2;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_b_b(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_t_t(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.m_t = param.m_t;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_t1 = param.Gamma_t1;
    cparam.Gamma_t2 = param.Gamma_t2;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_t_t(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_W_W(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_c1 = param.Gamma_c1;
    cparam.Gamma_c2 = param.Gamma_c2;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_W_W(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_Z_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_Z_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_h_h(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_h_h(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_h_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_h_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_Hp_Hp(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_c1 = param.Gamma_c1;
    cparam.Gamma_c2 = param.Gamma_c2;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_Hp_Hp(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_anti_Hp_W(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_c1 = param.Gamma_c1;
    cparam.Gamma_c2 = param.Gamma_c2;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_anti_Hp_W(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_H0_H0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_H0_H0(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_h_H0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_h_H0(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_H0_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_H0_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_A0_A0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_A0_A0(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_h_A0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_h_A0(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_H0_A0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_Z = param.Gamma_Z;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_A0 = param.Gamma_A0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_H0_A0(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_N_1_N_1_to_A0_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.alpha = param.alpha;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.Gamma_h = param.Gamma_h;
    cparam.theta_W = param.theta_W;
    cparam.Gamma_H0 = param.Gamma_H0;
    cparam.Gamma_o1 = param.Gamma_o1;
    cparam.Gamma_o2 = param.Gamma_o2;
    cparam.Gamma_o3 = param.Gamma_o3;
    cparam.Gamma_o4 = param.Gamma_o4;
    cparam.reg_prop = param.reg_prop;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    auto res = c_sumSqAmpl_N_1_N_1_to_A0_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_mu_mu(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_tau_tau(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_s_s(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_s = param.m_s;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_s_s(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_c_c(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_c = param.m_c;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_c_c(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_b_b(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_b = param.m_b;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_b_b(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_t_t(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_t = param.m_t;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_t_t(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_anti_W_W(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_anti_W_W(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_h_to_Z_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_h_to_Z_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_mu_mu(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_tau_tau(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_s_s(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_s = param.m_s;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_s_s(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_c_c(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_c = param.m_c;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_c_c(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_b_b(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_b = param.m_b;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_b_b(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_t_t(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_t = param.m_t;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_t_t(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_anti_W_W(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_anti_W_W(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_Z_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_Z_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_H0_to_h_h(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_H0_to_h_h(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_anti_mu_mu(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.s_23 = param.s_23;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_anti_tau_tau(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_anti_s_s(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_s = param.m_s;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_anti_s_s(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_anti_c_c(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_c = param.m_c;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_anti_c_c(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_anti_b_b(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_b = param.m_b;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_anti_b_b(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_anti_t_t(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.m_t = param.m_t;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_anti_t_t(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_A0_to_h_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_Z = param.m_Z;
    cparam.m_h = param.m_h;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.s_23 = param.s_23;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_A0_to_h_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_Hp_to_W_h(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_Hp = param.m_Hp;
    cparam.s_12 = param.s_12;
    cparam.alpha = param.alpha;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_Hp_to_W_h(&cparam);
    return {res.real, res.imag};
  }

  complex_t sumSqAmpl_Hp_to_W_A0(param_t const& param)
  {
    cparam_t cparam;
    cparam.m_W = param.m_W;
    cparam.e_em = param.e_em;
    cparam.m_Hp = param.m_Hp;
    cparam.s_12 = param.s_12;
    cparam.theta_W = param.theta_W;
    auto res = c_sumSqAmpl_Hp_to_W_A0(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_h(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.A_b = param.A_b;
    cparam.A_t = param.A_t;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_mu = param.m_mu;
    cparam.mu_h = param.mu_h;
    cparam.A_tau = param.A_tau;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_tau = param.m_tau;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sc_R = param.m_sc_R;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.m_se_L = param.m_se_L;
    cparam.m_se_R = param.m_se_R;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.m_su_L = param.m_su_L;
    cparam.m_su_R = param.m_su_R;
    cparam.m_smu_L = param.m_smu_L;
    cparam.m_smu_R = param.m_smu_R;
    cparam.theta_W = param.theta_W;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.m_snu_tau = param.m_snu_tau;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_tb = param.V_tb.get().real() + _mty_I * param.V_tb.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    cparam.U_stau_10 = param.U_stau_10.get().real() + _mty_I * param.U_stau_10.get().imag();
    cparam.U_stau_11 = param.U_stau_11.get().real() + _mty_I * param.U_stau_11.get().imag();
    auto res = c_width_h(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_W(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.V_cb = param.V_cb;
    cparam.V_ud = param.V_ud;
    cparam.V_us = param.V_us;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_mu = param.m_mu;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_tau = param.m_tau;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_se_L = param.m_se_L;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.m_su_L = param.m_su_L;
    cparam.m_smu_L = param.m_smu_L;
    cparam.theta_W = param.theta_W;
    cparam.V_ub_mod = param.V_ub_mod;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.m_snu_tau = param.m_snu_tau;
    cparam.delta_wolf = param.delta_wolf;
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_cd = param.V_cd.get().real() + _mty_I * param.V_cd.get().imag();
    cparam.V_cs = param.V_cs.get().real() + _mty_I * param.V_cs.get().imag();
    cparam.V_tb = param.V_tb.get().real() + _mty_I * param.V_tb.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    auto res = c_width_W(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_Z(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_mu = param.m_mu;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_tau = param.m_tau;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sc_R = param.m_sc_R;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.m_se_L = param.m_se_L;
    cparam.m_se_R = param.m_se_R;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.m_su_L = param.m_su_L;
    cparam.m_su_R = param.m_su_R;
    cparam.m_smu_L = param.m_smu_L;
    cparam.m_smu_R = param.m_smu_R;
    cparam.theta_W = param.theta_W;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.m_snu_tau = param.m_snu_tau;
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    cparam.U_stau_10 = param.U_stau_10.get().real() + _mty_I * param.U_stau_10.get().imag();
    cparam.U_stau_11 = param.U_stau_11.get().real() + _mty_I * param.U_stau_11.get().imag();
    auto res = c_width_Z(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_t(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.g_s = param.g_s;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_d = param.m_d;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_sG = param.m_sG;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.theta_W = param.theta_W;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_tb = param.V_tb.get().real() + _mty_I * param.V_tb.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    auto res = c_width_t(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_H0(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.A_b = param.A_b;
    cparam.A_t = param.A_t;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_mu = param.m_mu;
    cparam.mu_h = param.mu_h;
    cparam.A_tau = param.A_tau;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_tau = param.m_tau;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sc_R = param.m_sc_R;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.m_se_L = param.m_se_L;
    cparam.m_se_R = param.m_se_R;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.m_su_L = param.m_su_L;
    cparam.m_su_R = param.m_su_R;
    cparam.m_smu_L = param.m_smu_L;
    cparam.m_smu_R = param.m_smu_R;
    cparam.theta_W = param.theta_W;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.m_snu_tau = param.m_snu_tau;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_tb = param.V_tb.get().real() + _mty_I * param.V_tb.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    cparam.U_stau_10 = param.U_stau_10.get().real() + _mty_I * param.U_stau_10.get().imag();
    cparam.U_stau_11 = param.U_stau_11.get().real() + _mty_I * param.U_stau_11.get().imag();
    auto res = c_width_H0(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_Hp(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.A_b = param.A_b;
    cparam.A_t = param.A_t;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.V_cb = param.V_cb;
    cparam.V_ud = param.V_ud;
    cparam.V_us = param.V_us;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_mu = param.m_mu;
    cparam.mu_h = param.mu_h;
    cparam.A_tau = param.A_tau;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_tau = param.m_tau;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sc_R = param.m_sc_R;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.m_se_L = param.m_se_L;
    cparam.m_se_R = param.m_se_R;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.m_su_L = param.m_su_L;
    cparam.m_su_R = param.m_su_R;
    cparam.m_smu_L = param.m_smu_L;
    cparam.m_smu_R = param.m_smu_R;
    cparam.theta_W = param.theta_W;
    cparam.V_ub_mod = param.V_ub_mod;
    cparam.m_snu_mu = param.m_snu_mu;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.m_snu_tau = param.m_snu_tau;
    cparam.delta_wolf = param.delta_wolf;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_cd = param.V_cd.get().real() + _mty_I * param.V_cd.get().imag();
    cparam.V_cs = param.V_cs.get().real() + _mty_I * param.V_cs.get().imag();
    cparam.V_tb = param.V_tb.get().real() + _mty_I * param.V_tb.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    cparam.U_stau_10 = param.U_stau_10.get().real() + _mty_I * param.U_stau_10.get().imag();
    cparam.U_stau_11 = param.U_stau_11.get().real() + _mty_I * param.U_stau_11.get().imag();
    auto res = c_width_Hp(&cparam);
    return {res.real, res.imag};
  }

  complex_t width_A0(param_t const& param)
  {
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.A_b = param.A_b;
    cparam.A_t = param.A_t;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.beta = param.beta;
    cparam.e_em = param.e_em;
    cparam.m_A0 = param.m_A0;
    cparam.m_H0 = param.m_H0;
    cparam.m_Hp = param.m_Hp;
    cparam.m_mu = param.m_mu;
    cparam.mu_h = param.mu_h;
    cparam.A_tau = param.A_tau;
    cparam.alpha = param.alpha;
    cparam.m_C_1 = param.m_C_1;
    cparam.m_C_2 = param.m_C_2;
    cparam.m_N_1 = param.m_N_1;
    cparam.m_N_2 = param.m_N_2;
    cparam.m_N_3 = param.m_N_3;
    cparam.m_N_4 = param.m_N_4;
    cparam.m_tau = param.m_tau;
    cparam.m_sb_1 = param.m_sb_1;
    cparam.m_sb_2 = param.m_sb_2;
    cparam.m_sc_L = param.m_sc_L;
    cparam.m_sc_R = param.m_sc_R;
    cparam.m_sd_L = param.m_sd_L;
    cparam.m_sd_R = param.m_sd_R;
    cparam.m_se_L = param.m_se_L;
    cparam.m_se_R = param.m_se_R;
    cparam.m_ss_L = param.m_ss_L;
    cparam.m_ss_R = param.m_ss_R;
    cparam.m_st_1 = param.m_st_1;
    cparam.m_st_2 = param.m_st_2;
    cparam.m_su_L = param.m_su_L;
    cparam.m_su_R = param.m_su_R;
    cparam.m_smu_L = param.m_smu_L;
    cparam.m_smu_R = param.m_smu_R;
    cparam.theta_W = param.theta_W;
    cparam.m_stau_1 = param.m_stau_1;
    cparam.m_stau_2 = param.m_stau_2;
    cparam.N_B1 = param.N_B1.get().real() + _mty_I * param.N_B1.get().imag();
    cparam.N_B2 = param.N_B2.get().real() + _mty_I * param.N_B2.get().imag();
    cparam.N_B3 = param.N_B3.get().real() + _mty_I * param.N_B3.get().imag();
    cparam.N_B4 = param.N_B4.get().real() + _mty_I * param.N_B4.get().imag();
    cparam.N_W1 = param.N_W1.get().real() + _mty_I * param.N_W1.get().imag();
    cparam.N_W2 = param.N_W2.get().real() + _mty_I * param.N_W2.get().imag();
    cparam.N_W3 = param.N_W3.get().real() + _mty_I * param.N_W3.get().imag();
    cparam.N_W4 = param.N_W4.get().real() + _mty_I * param.N_W4.get().imag();
    cparam.N_d1 = param.N_d1.get().real() + _mty_I * param.N_d1.get().imag();
    cparam.N_d2 = param.N_d2.get().real() + _mty_I * param.N_d2.get().imag();
    cparam.N_d3 = param.N_d3.get().real() + _mty_I * param.N_d3.get().imag();
    cparam.N_d4 = param.N_d4.get().real() + _mty_I * param.N_d4.get().imag();
    cparam.N_u1 = param.N_u1.get().real() + _mty_I * param.N_u1.get().imag();
    cparam.N_u2 = param.N_u2.get().real() + _mty_I * param.N_u2.get().imag();
    cparam.N_u3 = param.N_u3.get().real() + _mty_I * param.N_u3.get().imag();
    cparam.N_u4 = param.N_u4.get().real() + _mty_I * param.N_u4.get().imag();
    cparam.U_d1 = param.U_d1.get().real() + _mty_I * param.U_d1.get().imag();
    cparam.U_d2 = param.U_d2.get().real() + _mty_I * param.U_d2.get().imag();
    cparam.V_tb = param.V_tb.get().real() + _mty_I * param.V_tb.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I * param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I * param.V_ts.get().imag();
    cparam.V_u1 = param.V_u1.get().real() + _mty_I * param.V_u1.get().imag();
    cparam.V_u2 = param.V_u2.get().real() + _mty_I * param.V_u2.get().imag();
    cparam.U_Wm1 = param.U_Wm1.get().real() + _mty_I * param.U_Wm1.get().imag();
    cparam.U_Wm2 = param.U_Wm2.get().real() + _mty_I * param.U_Wm2.get().imag();
    cparam.V_Wp1 = param.V_Wp1.get().real() + _mty_I * param.V_Wp1.get().imag();
    cparam.V_Wp2 = param.V_Wp2.get().real() + _mty_I * param.V_Wp2.get().imag();
    cparam.U_sb_00 = param.U_sb_00.get().real() + _mty_I * param.U_sb_00.get().imag();
    cparam.U_sb_01 = param.U_sb_01.get().real() + _mty_I * param.U_sb_01.get().imag();
    cparam.U_sb_10 = param.U_sb_10.get().real() + _mty_I * param.U_sb_10.get().imag();
    cparam.U_sb_11 = param.U_sb_11.get().real() + _mty_I * param.U_sb_11.get().imag();
    cparam.U_st_00 = param.U_st_00.get().real() + _mty_I * param.U_st_00.get().imag();
    cparam.U_st_01 = param.U_st_01.get().real() + _mty_I * param.U_st_01.get().imag();
    cparam.U_st_10 = param.U_st_10.get().real() + _mty_I * param.U_st_10.get().imag();
    cparam.U_st_11 = param.U_st_11.get().real() + _mty_I * param.U_st_11.get().imag();
    cparam.U_stau_00 = param.U_stau_00.get().real() + _mty_I * param.U_stau_00.get().imag();
    cparam.U_stau_01 = param.U_stau_01.get().real() + _mty_I * param.U_stau_01.get().imag();
    cparam.U_stau_10 = param.U_stau_10.get().real() + _mty_I * param.U_stau_10.get().imag();
    cparam.U_stau_11 = param.U_stau_11.get().real() + _mty_I * param.U_stau_11.get().imag();
    auto res = c_width_A0(&cparam);
    return {res.real, res.imag};
  }

} // End of namespace mssm2to2
