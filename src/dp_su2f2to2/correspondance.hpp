#pragma once

#ifndef __cplusplus
  #define __cplusplus
#endif

#include "clib_su2f2to2.h"
#include "su2f2to2.h"
#include <unordered_map>
namespace su2f2to2
{
// We are using char 32 (space) as possible separator, char 126 (tilde) as antiParticle tag
#define EMPTYCHAR 32
#define ANTICHAR  126
  using Cfptr_t = ccomplex_return_t (*)(cparam_t const*);
  namespace corr
  {
    using Entry_t = std::tuple<Cfptr_t, short int, short int>;
    enum Part_t
    {
      G = 1,
      V_3,
      V,
      phi,
      W,
      h,
      A,
      Z,
      nu_e,
      nu_mu,
      nu_tau,
      chi_dm_1,
      chi_dm_2,
      e,
      mu,
      tau,
      u,
      c,
      t,
      d,
      s,
      b
    };

#include "smBsm.hpp"
    inline const std::array<std::string, TOTAL_PARTICLES + 1> part_names = {
        " ",        "G",        "V_3", "V",  "phi", "W", "h", "A", "Z", "nu_e", "nu_mu", "nu_tau",
        "chi_dm_1", "chi_dm_2", "e",   "mu", "tau", "u", "c", "t", "d", "s",    "b"};

    inline const std::array<bool, TOTAL_PARTICLES + 1> isboson = {
        false, true,  true,  true,  true,  true,  true,  true,  true,  false, false, false,
        false, false, false, false, false, false, false, false, false, false, false};

    inline const std::string mass_codes = {
        EMPTYCHAR,         EMPTYCHAR + G,      EMPTYCHAR + V_3,      EMPTYCHAR + V,        EMPTYCHAR + phi,
        EMPTYCHAR + W,     EMPTYCHAR + h,      EMPTYCHAR + A,        EMPTYCHAR + Z,        EMPTYCHAR + nu_e,
        EMPTYCHAR + nu_mu, EMPTYCHAR + nu_tau, EMPTYCHAR + chi_dm_1, EMPTYCHAR + chi_dm_2, EMPTYCHAR + e,
        EMPTYCHAR + mu,    EMPTYCHAR + tau,    EMPTYCHAR + u,        EMPTYCHAR + c,        EMPTYCHAR + t,
        EMPTYCHAR + d,     EMPTYCHAR + s,      EMPTYCHAR + b};

    inline const std::array<double, TOTAL_PARTICLES + 1> part_charge = {
        0, 0, 0,  0,  0,  1,        0,        0,        0,         0,         0,        0,
        0, 0, -1, -1, -1, 0.666667, 0.666667, 0.666667, -0.333333, -0.333333, -0.333333};

    inline const std::array<int, TOTAL_PARTICLES + 1> part_tot_dof = {0, 16, 3, 6, 1, 6,  1,  2,  3,  4,  4, 4,
                                                                      4, 4,  4, 4, 4, 12, 12, 12, 12, 12, 12};

    inline const std::array<int, TOTAL_PARTICLES + 1> part_hel_dof = {0, 16, 3, 3, 1, 3, 1, 2, 3, 2, 2, 2,
                                                                      2, 2,  2, 2, 2, 6, 6, 6, 6, 6, 6};

    inline const std::array<bool, TOTAL_PARTICLES + 1> part_isSelfConj = {
        false, true,  true,  false, true,  false, true,  true,  true,  false, false, false,
        false, false, false, false, false, false, false, false, false, false, false};


    extern const std::unordered_map<std::string, Entry_t> squaredampl;

    inline double getMassFirst(const int p, const param_t& input)
    {
      switch (p)
      {
      case (G):
        return 0.0;
      case (V_3):
        return static_cast<double>(input.m_V);
      case (V):
        return static_cast<double>(input.m_V);
      case (phi):
        return static_cast<double>(input.m_phi);
      case (W):
        return static_cast<double>(input.m_W);
      case (h):
        return static_cast<double>(input.m_h);
      case (A):
        return 0.0;
      case (Z):
        return static_cast<double>(input.m_Z);
      case (nu_e):
        return 0.0;
      case (nu_mu):
        return 0.0;
      case (nu_tau):
        return 0.0;
      case (chi_dm_1):
        return static_cast<double>(input.m_chi_dm_1);
      case (chi_dm_2):
        return static_cast<double>(input.m_chi_dm_2);
      case (e):
        return static_cast<double>(input.m_e);
      case (mu):
        return static_cast<double>(input.m_mu);
      case (tau):
        return static_cast<double>(input.m_tau);
      case (u):
        return static_cast<double>(input.m_u);
      case (c):
        return static_cast<double>(input.m_c);
      case (t):
        return static_cast<double>(input.m_t);
      case (d):
        return static_cast<double>(input.m_d);
      case (s):
        return static_cast<double>(input.m_s);
      case (b):
        return static_cast<double>(input.m_b);
      default:
        return -1.;
      }
    } // End of getMassFirst
  }   // End of namespace corr
} // End of namespace su2f2to2
