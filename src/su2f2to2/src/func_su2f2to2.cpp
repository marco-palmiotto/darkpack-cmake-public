#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "stdbool.h"
#include "func_su2f2to2.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"
#include "cparams.h"
#include "clib_su2f2to2.h"
#include <complex.h>

namespace su2f2to2 {

complex_t m_chi_dm_1(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.m_cdm = param.m_cdm;
    auto res = c_m_chi_dm_1(&cparam);
    return {res.real, res.imag};
}

complex_t m_chi_dm_2(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.m_cdm = param.m_cdm;
    auto res = c_m_chi_dm_2(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_phi(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_phi(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_V(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_V(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_tau_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_tau_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_phi(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_phi(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_tau_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_tau_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_tau_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.Gamma_V = param.Gamma_V;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_phi(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_23 = param.s_23;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_phi(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_V(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_V(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_tau_nu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_23_L = param.theta_23_L;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_tau_nu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_tau(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_tau_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.s_12 = param.s_12;
    cparam.s_13 = param.s_13;
    cparam.s_14 = param.s_14;
    cparam.s_23 = param.s_23;
    cparam.s_24 = param.s_24;
    cparam.s_34 = param.s_34;
    cparam.m_tau = param.m_tau;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_V_3 = param.Gamma_V_3;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_tau_tau(&cparam);
    return {res.real, res.imag};
}

complex_t width_h(
        param_t const &param
        )
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
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    auto res = c_width_h(&cparam);
    return {res.real, res.imag};
}

complex_t width_W(
        param_t const &param
        )
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
    cparam.V_tb = param.V_tb;
    cparam.V_ud = param.V_ud;
    cparam.V_us = param.V_us;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    cparam.V_cd = param.V_cd.get().real() + _mty_I*param.V_cd.get().imag();
    cparam.V_cs = param.V_cs.get().real() + _mty_I*param.V_cs.get().imag();
    cparam.V_td = param.V_td.get().real() + _mty_I*param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I*param.V_ts.get().imag();
    cparam.V_ub = param.V_ub.get().real() + _mty_I*param.V_ub.get().imag();
    auto res = c_width_W(&cparam);
    return {res.real, res.imag};
}

complex_t width_Z(
        param_t const &param
        )
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
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    auto res = c_width_Z(&cparam);
    return {res.real, res.imag};
}

complex_t width_t(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_d = param.m_d;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.V_tb = param.V_tb;
    cparam.e_em = param.e_em;
    cparam.theta_W = param.theta_W;
    cparam.V_td = param.V_td.get().real() + _mty_I*param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I*param.V_ts.get().imag();
    auto res = c_width_t(&cparam);
    return {res.real, res.imag};
}

complex_t width_V(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.m_mu = param.m_mu;
    cparam.m_phi = param.m_phi;
    cparam.m_tau = param.m_tau;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_width_V(&cparam);
    return {res.real, res.imag};
}

complex_t width_V_3(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.g_f = param.g_f;
    cparam.m_V = param.m_V;
    cparam.m_e = param.m_e;
    cparam.m_mu = param.m_mu;
    cparam.m_phi = param.m_phi;
    cparam.m_tau = param.m_tau;
    cparam.m_chi_dm_1 = param.m_chi_dm_1;
    cparam.m_chi_dm_2 = param.m_chi_dm_2;
    cparam.theta_12_L = param.theta_12_L;
    cparam.theta_12_R = param.theta_12_R;
    cparam.theta_13_L = param.theta_13_L;
    cparam.theta_13_R = param.theta_13_R;
    cparam.theta_23_L = param.theta_23_L;
    cparam.theta_23_R = param.theta_23_R;
    auto res = c_width_V_3(&cparam);
    return {res.real, res.imag};
}

} // End of namespace su2f2to2
