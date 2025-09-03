#include "compare.hpp"
#include "config.hpp"
#include "params_new.hpp"

namespace __SPEC_LIB_NAME__
{

  bool compare_param_t(const param_t& p1, const param_t& p2)
  {
    if (!(p1.pi == p2.pi))
      return false;
    if (!(p1.g_f == p2.g_f))
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
    if (!(p1.s_13 == p2.s_13))
      return false;
    if (!(p1.s_14 == p2.s_14))
      return false;
    if (!(p1.s_23 == p2.s_23))
      return false;
    if (!(p1.s_24 == p2.s_24))
      return false;
    if (!(p1.s_34 == p2.s_34))
      return false;
    if (!(p1.m_cdm == p2.m_cdm))
      return false;
    if (!(p1.theta_W == p2.theta_W))
      return false;
    if (!(p1.reg_prop == p2.reg_prop))
      return false;
    if (!(p1.Gamma_V_3 == p2.Gamma_V_3))
      return false;
    if (!(p1.theta_12_L == p2.theta_12_L))
      return false;
    if (!(p1.theta_12_R == p2.theta_12_R))
      return false;
    if (!(p1.theta_13_L == p2.theta_13_L))
      return false;
    if (!(p1.theta_13_R == p2.theta_13_R))
      return false;
    if (!(p1.theta_23_L == p2.theta_23_L))
      return false;
    if (!(p1.theta_23_R == p2.theta_23_R))
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
    if (!(p1.m_V_3 == p2.m_V_3))
      return false;
    if (!(p1.m_phi == p2.m_phi))
      return false;
    if (!(p1.m_tau == p2.m_tau))
      return false;
    if (!(p1.m_chi_dm_1 == p2.m_chi_dm_1))
      return false;
    return true;
  }
} // namespace __SPEC_LIB_NAME__
