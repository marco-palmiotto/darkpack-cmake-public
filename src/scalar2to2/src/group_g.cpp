#include "group_g.h"
namespace scalar2to2 {

const std::array<Callable<complex_t, param_t>, 14> f_G = {
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_u_u", sumSqAmpl_anti_chi_chi_to_anti_u_u},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_c_c", sumSqAmpl_anti_chi_chi_to_anti_c_c},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_t_t", sumSqAmpl_anti_chi_chi_to_anti_t_t},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_d_d", sumSqAmpl_anti_chi_chi_to_anti_d_d},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_s_s", sumSqAmpl_anti_chi_chi_to_anti_s_s},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_b_b", sumSqAmpl_anti_chi_chi_to_anti_b_b},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_e_e", sumSqAmpl_anti_chi_chi_to_anti_e_e},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_mu_mu", sumSqAmpl_anti_chi_chi_to_anti_mu_mu},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_tau_tau", sumSqAmpl_anti_chi_chi_to_anti_tau_tau},
    Callable{"width_h", width_h},
    Callable{"width_W", width_W},
    Callable{"width_Z", width_Z},
    Callable{"width_t", width_t},
    Callable{"width_phi", width_phi},
};

const std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"sumSqAmpl_anti_chi_chi_to_anti_u_u", f_G[0]},
    {"sumSqAmpl_anti_chi_chi_to_anti_c_c", f_G[1]},
    {"sumSqAmpl_anti_chi_chi_to_anti_t_t", f_G[2]},
    {"sumSqAmpl_anti_chi_chi_to_anti_d_d", f_G[3]},
    {"sumSqAmpl_anti_chi_chi_to_anti_s_s", f_G[4]},
    {"sumSqAmpl_anti_chi_chi_to_anti_b_b", f_G[5]},
    {"sumSqAmpl_anti_chi_chi_to_anti_e_e", f_G[6]},
    {"sumSqAmpl_anti_chi_chi_to_anti_mu_mu", f_G[7]},
    {"sumSqAmpl_anti_chi_chi_to_anti_tau_tau", f_G[8]},
    {"width_h", f_G[9]},
    {"width_W", f_G[10]},
    {"width_Z", f_G[11]},
    {"width_t", f_G[12]},
};


} // End of namespace scalar2to2
