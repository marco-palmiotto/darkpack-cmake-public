#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <cmath>
#include "stdbool.h"
#include "func_scalar2to2.h"
#include "common.h"

#include "params.h"
#include "libcomplexop.h"
#include "cparams.h"
#include "clib_scalar2to2.h"
#include <complex.h>

namespace scalar2to2 {

complex_t sumSqAmpl_anti_chi_chi_to_anti_u_u(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_u = param.g_u;
    cparam.m_W = param.m_W;
    cparam.m_u = param.m_u;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_u_u(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_c_c(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_u = param.g_u;
    cparam.m_W = param.m_W;
    cparam.m_c = param.m_c;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_c_c(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_t_t(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_u = param.g_u;
    cparam.m_W = param.m_W;
    cparam.m_t = param.m_t;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_t_t(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_d_d(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_d = param.g_d;
    cparam.m_W = param.m_W;
    cparam.m_d = param.m_d;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_d_d(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_s_s(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_d = param.g_d;
    cparam.m_W = param.m_W;
    cparam.m_s = param.m_s;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_s_s(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_b_b(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_d = param.g_d;
    cparam.m_W = param.m_W;
    cparam.m_b = param.m_b;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_b_b(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_e_e(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_l = param.g_l;
    cparam.m_W = param.m_W;
    cparam.m_e = param.m_e;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_e_e(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_mu_mu(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_l = param.g_l;
    cparam.m_W = param.m_W;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_mu_mu(&cparam);
    return {res.real, res.imag};
}

complex_t sumSqAmpl_anti_chi_chi_to_anti_tau_tau(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.g_l = param.g_l;
    cparam.m_W = param.m_W;
    cparam.e_em = param.e_em;
    cparam.s_12 = param.s_12;
    cparam.s_34 = param.s_34;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    cparam.reg_prop = param.reg_prop;
    cparam.Gamma_phi = param.Gamma_phi;
    auto res = c_sumSqAmpl_anti_chi_chi_to_anti_tau_tau(&cparam);
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
    cparam.g_u = param.g_u;
    cparam.m_W = param.m_W;
    cparam.m_Z = param.m_Z;
    cparam.m_b = param.m_b;
    cparam.m_d = param.m_d;
    cparam.m_h = param.m_h;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.V_tb = param.V_tb;
    cparam.e_em = param.e_em;
    cparam.m_phi = param.m_phi;
    cparam.theta_W = param.theta_W;
    cparam.V_td = param.V_td.get().real() + _mty_I*param.V_td.get().imag();
    cparam.V_ts = param.V_ts.get().real() + _mty_I*param.V_ts.get().imag();
    auto res = c_width_t(&cparam);
    return {res.real, res.imag};
}

complex_t width_phi(
        param_t const &param
        )
{
    cparam_t cparam;
    cparam.pi = param.pi;
    cparam.g_d = param.g_d;
    cparam.g_l = param.g_l;
    cparam.g_u = param.g_u;
    cparam.m_W = param.m_W;
    cparam.m_b = param.m_b;
    cparam.m_c = param.m_c;
    cparam.m_d = param.m_d;
    cparam.m_e = param.m_e;
    cparam.m_s = param.m_s;
    cparam.m_t = param.m_t;
    cparam.m_u = param.m_u;
    cparam.e_em = param.e_em;
    cparam.m_mu = param.m_mu;
    cparam.g_chi = param.g_chi;
    cparam.m_chi = param.m_chi;
    cparam.m_phi = param.m_phi;
    cparam.m_tau = param.m_tau;
    cparam.theta_W = param.theta_W;
    auto res = c_width_phi(&cparam);
    return {res.real, res.imag};
}

} // End of namespace scalar2to2
