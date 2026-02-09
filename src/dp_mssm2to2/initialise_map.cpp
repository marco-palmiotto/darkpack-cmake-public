#include "correspondance.hpp"
namespace mssm2to2::corr{
const std::unordered_map<std::string, Entry_t> squaredampl{ 
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::e,EMPTYCHAR+corr::e},{  &c_sumSqAmpl_N_1_N_1_to_anti_e_e, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::nu_e,EMPTYCHAR+corr::nu_e},{  &c_sumSqAmpl_N_1_N_1_to_anti_nu_e_nu_e, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::mu,EMPTYCHAR+corr::mu},{  &c_sumSqAmpl_N_1_N_1_to_anti_mu_mu, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::nu_mu,EMPTYCHAR+corr::nu_mu},{  &c_sumSqAmpl_N_1_N_1_to_anti_nu_mu_nu_mu, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::tau,EMPTYCHAR+corr::tau},{  &c_sumSqAmpl_N_1_N_1_to_anti_tau_tau, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::nu_tau,EMPTYCHAR+corr::nu_tau},{  &c_sumSqAmpl_N_1_N_1_to_anti_nu_tau_nu_tau, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::u,EMPTYCHAR+corr::u},{  &c_sumSqAmpl_N_1_N_1_to_anti_u_u, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::d,EMPTYCHAR+corr::d},{  &c_sumSqAmpl_N_1_N_1_to_anti_d_d, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::s,EMPTYCHAR+corr::s},{  &c_sumSqAmpl_N_1_N_1_to_anti_s_s, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::c,EMPTYCHAR+corr::c},{  &c_sumSqAmpl_N_1_N_1_to_anti_c_c, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::b,EMPTYCHAR+corr::b},{  &c_sumSqAmpl_N_1_N_1_to_anti_b_b, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::t,EMPTYCHAR+corr::t},{  &c_sumSqAmpl_N_1_N_1_to_anti_t_t, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::W,EMPTYCHAR+corr::W},{  &c_sumSqAmpl_N_1_N_1_to_anti_W_W, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::Z,EMPTYCHAR+corr::Z},{  &c_sumSqAmpl_N_1_N_1_to_Z_Z, 8,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::h,EMPTYCHAR+corr::h},{  &c_sumSqAmpl_N_1_N_1_to_h_h, 8,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::h,EMPTYCHAR+corr::Z},{  &c_sumSqAmpl_N_1_N_1_to_h_Z, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::Hp,EMPTYCHAR+corr::Hp},{  &c_sumSqAmpl_N_1_N_1_to_anti_Hp_Hp, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,ANTICHAR,EMPTYCHAR+corr::Hp,EMPTYCHAR+corr::W},{  &c_sumSqAmpl_N_1_N_1_to_anti_Hp_W, 4,2} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::H0,EMPTYCHAR+corr::H0},{  &c_sumSqAmpl_N_1_N_1_to_H0_H0, 8,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::h,EMPTYCHAR+corr::H0},{  &c_sumSqAmpl_N_1_N_1_to_h_H0, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::H0,EMPTYCHAR+corr::Z},{  &c_sumSqAmpl_N_1_N_1_to_H0_Z, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::A0,EMPTYCHAR+corr::A0},{  &c_sumSqAmpl_N_1_N_1_to_A0_A0, 8,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::h,EMPTYCHAR+corr::A0},{  &c_sumSqAmpl_N_1_N_1_to_h_A0, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::H0,EMPTYCHAR+corr::A0},{  &c_sumSqAmpl_N_1_N_1_to_H0_A0, 4,1} },
{ {EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::N_1,EMPTYCHAR+corr::A0,EMPTYCHAR+corr::Z},{  &c_sumSqAmpl_N_1_N_1_to_A0_Z, 4,1} }};
}
