#include "correspondance.hpp"
namespace mssm2to2::corr
{
  const std::unordered_map<std::string, Entry_t_1to2> aux_processes{
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::b, EMPTYCHAR + corr::c},
       {
           &partWidth_W_to_b_c,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::c, EMPTYCHAR + corr::d},
       {
           &partWidth_W_to_c_d,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::c, EMPTYCHAR + corr::s},
       {
           &partWidth_W_to_c_s,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::b, EMPTYCHAR + corr::t},
       {
           &partWidth_W_to_b_t,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::d, EMPTYCHAR + corr::t},
       {
           &partWidth_W_to_d_t,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::s, EMPTYCHAR + corr::t},
       {
           &partWidth_W_to_s_t,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::b, EMPTYCHAR + corr::u},
       {
           &partWidth_W_to_b_u,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::d, EMPTYCHAR + corr::u},
       {
           &partWidth_W_to_d_u,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::s, EMPTYCHAR + corr::u},
       {
           &partWidth_W_to_s_u,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::h, EMPTYCHAR + corr::Hp},
       {
           &partWidth_W_to_h_Hp,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::H0, EMPTYCHAR + corr::Hp},
       {
           &partWidth_W_to_H0_Hp,
       }},
      {{EMPTYCHAR + corr::W, EMPTYCHAR + corr::A0, EMPTYCHAR + corr::Hp},
       {
           &partWidth_W_to_A0_Hp,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::W, EMPTYCHAR + corr::W},
       {
           &partWidth_Z_to_anti_W_W,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::b, EMPTYCHAR + corr::b},
       {
           &partWidth_Z_to_anti_b_b,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::c, EMPTYCHAR + corr::c},
       {
           &partWidth_Z_to_anti_c_c,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::d, EMPTYCHAR + corr::d},
       {
           &partWidth_Z_to_anti_d_d,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::e, EMPTYCHAR + corr::e},
       {
           &partWidth_Z_to_anti_e_e,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::s, EMPTYCHAR + corr::s},
       {
           &partWidth_Z_to_anti_s_s,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::t, EMPTYCHAR + corr::t},
       {
           &partWidth_Z_to_anti_t_t,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::u, EMPTYCHAR + corr::u},
       {
           &partWidth_Z_to_anti_u_u,
       }},
      {{EMPTYCHAR + corr::Z, EMPTYCHAR + corr::h, EMPTYCHAR + corr::A0},
       {
           &partWidth_Z_to_h_A0,
       }},
      {{EMPTYCHAR + corr::Z, EMPTYCHAR + corr::A0, EMPTYCHAR + corr::H0},
       {
           &partWidth_Z_to_A0_H0,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::Hp, EMPTYCHAR + corr::Hp},
       {
           &partWidth_Z_to_anti_Hp_Hp,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::mu, EMPTYCHAR + corr::mu},
       {
           &partWidth_Z_to_anti_mu_mu,
       }},
      {{EMPTYCHAR + corr::Z, ANTICHAR, EMPTYCHAR + corr::tau, EMPTYCHAR + corr::tau},
       {
           &partWidth_Z_to_anti_tau_tau,
       }},
      {{EMPTYCHAR + corr::t, EMPTYCHAR + corr::W, EMPTYCHAR + corr::b},
       {
           &partWidth_t_to_W_b,
       }},
      {{EMPTYCHAR + corr::t, EMPTYCHAR + corr::W, EMPTYCHAR + corr::d},
       {
           &partWidth_t_to_W_d,
       }},
      {{EMPTYCHAR + corr::t, EMPTYCHAR + corr::W, EMPTYCHAR + corr::s},
       {
           &partWidth_t_to_W_s,
       }},
      {{EMPTYCHAR + corr::t, EMPTYCHAR + corr::b, EMPTYCHAR + corr::Hp},
       {
           &partWidth_t_to_b_Hp,
       }},
      {{EMPTYCHAR + corr::t, EMPTYCHAR + corr::d, EMPTYCHAR + corr::Hp},
       {
           &partWidth_t_to_d_Hp,
       }},
      {{EMPTYCHAR + corr::t, EMPTYCHAR + corr::s, EMPTYCHAR + corr::Hp},
       {
           &partWidth_t_to_s_Hp,
       }}};
}
