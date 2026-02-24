#pragma once

#ifndef __cplusplus
  #define __cplusplus
#endif

#include "clib_mssm2to2.h"
#include "mssm2to2.h"
#include <unordered_map>
namespace mssm2to2
{
// We are using char 32 (space) as possible separator, char 126 (tilde) as antiParticle tag
#define EMPTYCHAR 32
#define ANTICHAR  126
  using Cfptr_t = ccomplex_return_t (*)(cparam_t const*);
  using CXXfptr_t = complex_t (*)(const param_t&);
  namespace corr
  {
    using Entry_t = std::tuple<Cfptr_t, short int, short int>;
    using Entry_t_1to2 = std::tuple<CXXfptr_t, short int, short int>;
    enum Part_t
    {
      G = 1,
      W,
      A,
      Z,
      A0,
      Hp,
      h,
      H0,
      sc_L,
      su_L,
      sc_R,
      su_R,
      sd_R,
      ss_R,
      sd_L,
      ss_L,
      se_R,
      smu_R,
      se_L,
      smu_L,
      snu_e,
      snu_mu,
      snu_tau,
      st_1,
      st_2,
      sb_1,
      sb_2,
      stau_1,
      stau_2,
      nu_e,
      nu_mu,
      nu_tau,
      sG,
      N_1,
      N_2,
      N_3,
      N_4,
      e,
      mu,
      tau,
      u,
      d,
      s,
      c,
      t,
      b,
      C_1,
      C_2
    };

#include "smBsm.hpp"
    inline const std::array<std::string, TOTAL_PARTICLES + 1> part_names = {
        " ",     "G",     "W",      "A",       "Z",    "A0",   "Hp",   "h",    "H0",     "sc_L",
        "su_L",  "sc_R",  "su_R",   "sd_R",    "ss_R", "sd_L", "ss_L", "se_R", "smu_R",  "se_L",
        "smu_L", "snu_e", "snu_mu", "snu_tau", "st_1", "st_2", "sb_1", "sb_2", "stau_1", "stau_2",
        "nu_e",  "nu_mu", "nu_tau", "sG",      "N_1",  "N_2",  "N_3",  "N_4",  "e",      "mu",
        "tau",   "u",     "d",      "s",       "c",    "t",    "b",    "C_1",  "C_2"};

    inline const std::array<bool, TOTAL_PARTICLES + 1> isboson = {
        false, true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
        true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,  true,
        true,  true,  true,  true,  false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, false, false, false};

    inline const std::string mass_codes = {
        EMPTYCHAR,         EMPTYCHAR + G,     EMPTYCHAR + W,      EMPTYCHAR + A,       EMPTYCHAR + Z,
        EMPTYCHAR + A0,    EMPTYCHAR + Hp,    EMPTYCHAR + h,      EMPTYCHAR + H0,      EMPTYCHAR + sc_L,
        EMPTYCHAR + su_L,  EMPTYCHAR + sc_R,  EMPTYCHAR + su_R,   EMPTYCHAR + sd_R,    EMPTYCHAR + ss_R,
        EMPTYCHAR + sd_L,  EMPTYCHAR + ss_L,  EMPTYCHAR + se_R,   EMPTYCHAR + smu_R,   EMPTYCHAR + se_L,
        EMPTYCHAR + smu_L, EMPTYCHAR + snu_e, EMPTYCHAR + snu_mu, EMPTYCHAR + snu_tau, EMPTYCHAR + st_1,
        EMPTYCHAR + st_2,  EMPTYCHAR + sb_1,  EMPTYCHAR + sb_2,   EMPTYCHAR + stau_1,  EMPTYCHAR + stau_2,
        EMPTYCHAR + nu_e,  EMPTYCHAR + nu_mu, EMPTYCHAR + nu_tau, EMPTYCHAR + sG,      EMPTYCHAR + N_1,
        EMPTYCHAR + N_2,   EMPTYCHAR + N_3,   EMPTYCHAR + N_4,    EMPTYCHAR + e,       EMPTYCHAR + mu,
        EMPTYCHAR + tau,   EMPTYCHAR + u,     EMPTYCHAR + d,      EMPTYCHAR + s,       EMPTYCHAR + c,
        EMPTYCHAR + t,     EMPTYCHAR + b,     EMPTYCHAR + C_1,    EMPTYCHAR + C_2};

    inline const std::array<double, TOTAL_PARTICLES + 1> part_charge = {
        0,         0,         1,  0,  0,  0,  1, 0, 0,  0.666667, 0.666667, 0.666667,  0.666667,  -0.333333, -0.333333,
        -0.333333, -0.333333, -1, -1, -1, -1, 0, 0, 0,  0.666667, 0.666667, -0.333333, -0.333333, -1,        -1,
        0,         0,         0,  0,  0,  0,  0, 0, -1, -1,       -1,       0.666667,  -0.333333, -0.333333, 0.666667,
        0.666667,  -0.333333, 1,  1};

    inline const std::array<int, TOTAL_PARTICLES + 1> part_tot_dof = {
        0, 16, 6, 2, 3, 1, 2, 1, 1,  6, 6, 6, 6, 6, 6, 6, 6,  2,  2,  2,  2,  2,  2, 2, 6,
        6, 6,  6, 2, 2, 4, 4, 4, 16, 2, 2, 2, 2, 4, 4, 4, 12, 12, 12, 12, 12, 12, 4, 4};

    inline const std::array<int, TOTAL_PARTICLES + 1> part_hel_dof = {
        0, 16, 3, 2, 3, 1, 1, 1, 1,  3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 3,
        3, 3,  3, 1, 1, 2, 2, 2, 16, 2, 2, 2, 2, 2, 2, 2, 6, 6, 6, 6, 6, 6, 2, 2};

    inline const std::array<bool, TOTAL_PARTICLES + 1> part_isSelfConj = {
        false, true,  false, true,  true,  true,  false, true,  true,  false, false, false, false,
        false, false, false, false, false, false, false, false, false, false, false, false, false,
        false, false, false, false, false, false, false, true,  true,  true,  true,  true,  false,
        false, false, false, false, false, false, false, false, false, false};


    extern const std::unordered_map<std::string, Entry_t> squaredampl;

    extern const std::unordered_map<std::string, Entry_t_1to2> squaredampl_1to2;

    inline double getMassFirst(const int p, const param_t& input)
    {
      switch (p)
      {
      case (G):
        return 0.0;
      case (W):
        return static_cast<double>(input.m_W);
      case (A):
        return 0.0;
      case (Z):
        return static_cast<double>(input.m_Z);
      case (A0):
        return static_cast<double>(input.m_A0);
      case (Hp):
        return static_cast<double>(input.m_Hp);
      case (h):
        return static_cast<double>(input.m_h);
      case (H0):
        return static_cast<double>(input.m_H0);
      case (sc_L):
        return static_cast<double>(input.m_sc_L);
      case (su_L):
        return static_cast<double>(input.m_su_L);
      case (sc_R):
        return static_cast<double>(input.m_sc_R);
      case (su_R):
        return static_cast<double>(input.m_su_R);
      case (sd_R):
        return static_cast<double>(input.m_sd_R);
      case (ss_R):
        return static_cast<double>(input.m_ss_R);
      case (sd_L):
        return static_cast<double>(input.m_sd_L);
      case (ss_L):
        return static_cast<double>(input.m_ss_L);
      case (se_R):
        return static_cast<double>(input.m_se_R);
      case (smu_R):
        return static_cast<double>(input.m_smu_R);
      case (se_L):
        return static_cast<double>(input.m_se_L);
      case (smu_L):
        return static_cast<double>(input.m_smu_L);
      case (snu_e):
        return static_cast<double>(input.m_snu_mu);
      case (snu_mu):
        return static_cast<double>(input.m_snu_mu);
      case (snu_tau):
        return static_cast<double>(input.m_snu_tau);
      case (st_1):
        return static_cast<double>(input.m_st_1);
      case (st_2):
        return static_cast<double>(input.m_st_2);
      case (sb_1):
        return static_cast<double>(input.m_sb_1);
      case (sb_2):
        return static_cast<double>(input.m_sb_2);
      case (stau_1):
        return static_cast<double>(input.m_stau_1);
      case (stau_2):
        return static_cast<double>(input.m_stau_2);
      case (nu_e):
        return 0.0;
      case (nu_mu):
        return 0.0;
      case (nu_tau):
        return 0.0;
      case (sG):
        return static_cast<double>(input.m_sG);
      case (N_1):
        return static_cast<double>(input.m_N_1);
      case (N_2):
        return static_cast<double>(input.m_N_2);
      case (N_3):
        return static_cast<double>(input.m_N_3);
      case (N_4):
        return static_cast<double>(input.m_N_4);
      case (e):
        return static_cast<double>(input.m_e);
      case (mu):
        return static_cast<double>(input.m_mu);
      case (tau):
        return static_cast<double>(input.m_tau);
      case (u):
        return static_cast<double>(input.m_u);
      case (d):
        return static_cast<double>(input.m_d);
      case (s):
        return static_cast<double>(input.m_s);
      case (c):
        return static_cast<double>(input.m_c);
      case (t):
        return static_cast<double>(input.m_t);
      case (b):
        return static_cast<double>(input.m_b);
      case (C_1):
        return static_cast<double>(input.m_C_1);
      case (C_2):
        return static_cast<double>(input.m_C_2);
      default:
        return -1.;
      }
    } // End of getMassFirst
  } // End of namespace corr
} // End of namespace mssm2to2
