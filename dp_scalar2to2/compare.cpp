#include "compare.hpp"
#include "config.hpp"
#include "params_new.hpp"

namespace __SPEC_LIB_NAME__
{

  bool compare_param_t(const param_t& p1, const param_t& p2)
  {
    if (!(p1.pi == p2.pi))
      return false;
    if (!(p1.g_d == p2.g_d))
      return false;
    if (!(p1.g_l == p2.g_l))
      return false;
    if (!(p1.g_u == p2.g_u))
      return false;
    if (!(p1.V_cb == p2.V_cb))
      return false;
    if (!(p1.V_tb == p2.V_tb))
      return false;
    if (!(p1.V_ud == p2.V_ud))
      return false;
    if (!(p1.V_us == p2.V_us))
      return false;
    if (!(p1.e_em == p2.e_em))
      return false;
    if (!(p1.s_12 == p2.s_12))
      return false;
    if (!(p1.s_34 == p2.s_34))
      return false;
    if (!(p1.g_chi == p2.g_chi))
      return false;
    if (!(p1.theta_W == p2.theta_W))
      return false;
    if (!(p1.reg_prop == p2.reg_prop))
      return false;
    if (!(p1.Gamma_phi == p2.Gamma_phi))
      return false;
    if (!(p1.V_cd == p2.V_cd))
      return false;
    if (!(p1.V_cs == p2.V_cs))
      return false;
    if (!(p1.V_td == p2.V_td))
      return false;
    if (!(p1.V_ts == p2.V_ts))
      return false;
    if (!(p1.V_ub == p2.V_ub))
      return false;
    if (!(p1.g_s == p2.g_s))
      return false;
    if (!(p1.m_W == p2.m_W))
      return false;
    if (!(p1.m_Z == p2.m_Z))
      return false;
    if (!(p1.m_b == p2.m_b))
      return false;
    if (!(p1.m_c == p2.m_c))
      return false;
    if (!(p1.m_d == p2.m_d))
      return false;
    if (!(p1.m_e == p2.m_e))
      return false;
    if (!(p1.m_h == p2.m_h))
      return false;
    if (!(p1.m_s == p2.m_s))
      return false;
    if (!(p1.m_t == p2.m_t))
      return false;
    if (!(p1.m_u == p2.m_u))
      return false;
    if (!(p1.m_mu == p2.m_mu))
      return false;
    if (!(p1.m_chi == p2.m_chi))
      return false;
    if (!(p1.m_phi == p2.m_phi))
      return false;
    if (!(p1.m_tau == p2.m_tau))
      return false;
    return true;
  }
} // namespace __SPEC_LIB_NAME__
