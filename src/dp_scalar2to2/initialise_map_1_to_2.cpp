#include "correspondance.hpp"
namespace scalar2to2::corr{
const std::unordered_map<std::string, Entry_t_1to2> squaredampl_1_to_2{ 
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::W,EMPTYCHAR+corr::W},{  &sumSqAmpl_h_to_anti_W_W, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::h,EMPTYCHAR+corr::h},{  &sumSqAmpl_h_to_anti_h_h, 2,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::Z,EMPTYCHAR+corr::Z},{  &sumSqAmpl_h_to_anti_Z_Z, 2,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::u,EMPTYCHAR+corr::u},{  &sumSqAmpl_h_to_anti_u_u, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::c,EMPTYCHAR+corr::c},{  &sumSqAmpl_h_to_anti_c_c, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::t,EMPTYCHAR+corr::t},{  &sumSqAmpl_h_to_anti_t_t, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::d,EMPTYCHAR+corr::d},{  &sumSqAmpl_h_to_anti_d_d, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::s,EMPTYCHAR+corr::s},{  &sumSqAmpl_h_to_anti_s_s, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::b,EMPTYCHAR+corr::b},{  &sumSqAmpl_h_to_anti_b_b, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::e,EMPTYCHAR+corr::e},{  &sumSqAmpl_h_to_anti_e_e, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::mu,EMPTYCHAR+corr::mu},{  &sumSqAmpl_h_to_anti_mu_mu, 1,1} },
{ {EMPTYCHAR+corr::h,ANTICHAR,EMPTYCHAR+corr::tau,EMPTYCHAR+corr::tau},{  &sumSqAmpl_h_to_anti_tau_tau, 1,1} }};
}
