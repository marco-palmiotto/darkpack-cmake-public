#include "correspondance.hpp"
namespace scalar2to2::corr
{
    const std::unordered_map<std::string, Entry_t> squaredampl{
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::u, EMPTYCHAR + corr::u},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_u_u, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::c, EMPTYCHAR + corr::c},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_c_c, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::t, EMPTYCHAR + corr::t},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_t_t, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::d, EMPTYCHAR + corr::d},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_d_d, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::s, EMPTYCHAR + corr::s},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_s_s, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::b, EMPTYCHAR + corr::b},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_b_b, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::e, EMPTYCHAR + corr::e},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_e_e, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::mu, EMPTYCHAR + corr::mu},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_mu_mu, 4, 1}},
        {{ANTICHAR, EMPTYCHAR + corr::chi, EMPTYCHAR + corr::chi, ANTICHAR, EMPTYCHAR + corr::tau,
          EMPTYCHAR + corr::tau},
         {&c_sumSqAmpl_anti_chi_chi_to_anti_tau_tau, 4, 1}}};
}
