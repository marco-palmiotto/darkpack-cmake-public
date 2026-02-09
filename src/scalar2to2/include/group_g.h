#ifndef CSL_LIB_scalar2to2_G_H_INCLUDED
#define CSL_LIB_scalar2to2_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_scalar2to2.h"

namespace scalar2to2 {


inline std::array<Callable<complex_t, param_t>, 26> f_G = {
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_u_u", sumSqAmpl_anti_chi_chi_to_anti_u_u},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_c_c", sumSqAmpl_anti_chi_chi_to_anti_c_c},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_t_t", sumSqAmpl_anti_chi_chi_to_anti_t_t},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_d_d", sumSqAmpl_anti_chi_chi_to_anti_d_d},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_s_s", sumSqAmpl_anti_chi_chi_to_anti_s_s},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_b_b", sumSqAmpl_anti_chi_chi_to_anti_b_b},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_e_e", sumSqAmpl_anti_chi_chi_to_anti_e_e},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_mu_mu", sumSqAmpl_anti_chi_chi_to_anti_mu_mu},
    Callable{"sumSqAmpl_anti_chi_chi_to_anti_tau_tau", sumSqAmpl_anti_chi_chi_to_anti_tau_tau},
    Callable{"sumSqAmpl_h_to_anti_W_W", sumSqAmpl_h_to_anti_W_W},
    Callable{"sumSqAmpl_h_to_anti_h_h", sumSqAmpl_h_to_anti_h_h},
    Callable{"sumSqAmpl_h_to_anti_Z_Z", sumSqAmpl_h_to_anti_Z_Z},
    Callable{"sumSqAmpl_h_to_anti_u_u", sumSqAmpl_h_to_anti_u_u},
    Callable{"sumSqAmpl_h_to_anti_c_c", sumSqAmpl_h_to_anti_c_c},
    Callable{"sumSqAmpl_h_to_anti_t_t", sumSqAmpl_h_to_anti_t_t},
    Callable{"sumSqAmpl_h_to_anti_d_d", sumSqAmpl_h_to_anti_d_d},
    Callable{"sumSqAmpl_h_to_anti_s_s", sumSqAmpl_h_to_anti_s_s},
    Callable{"sumSqAmpl_h_to_anti_b_b", sumSqAmpl_h_to_anti_b_b},
    Callable{"sumSqAmpl_h_to_anti_e_e", sumSqAmpl_h_to_anti_e_e},
    Callable{"sumSqAmpl_h_to_anti_mu_mu", sumSqAmpl_h_to_anti_mu_mu},
    Callable{"sumSqAmpl_h_to_anti_tau_tau", sumSqAmpl_h_to_anti_tau_tau},
    Callable{"width_h", width_h},
    Callable{"width_W", width_W},
    Callable{"width_Z", width_Z},
    Callable{"width_t", width_t},
    Callable{"width_phi", width_phi},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"sumSqAmpl_anti_chi_chi_to_anti_u_u", f_G[0]},
    {"sumSqAmpl_anti_chi_chi_to_anti_c_c", f_G[1]},
    {"sumSqAmpl_anti_chi_chi_to_anti_t_t", f_G[2]},
    {"sumSqAmpl_anti_chi_chi_to_anti_d_d", f_G[3]},
    {"sumSqAmpl_anti_chi_chi_to_anti_s_s", f_G[4]},
    {"sumSqAmpl_anti_chi_chi_to_anti_b_b", f_G[5]},
    {"sumSqAmpl_anti_chi_chi_to_anti_e_e", f_G[6]},
    {"sumSqAmpl_anti_chi_chi_to_anti_mu_mu", f_G[7]},
    {"sumSqAmpl_anti_chi_chi_to_anti_tau_tau", f_G[8]},
    {"sumSqAmpl_h_to_anti_W_W", f_G[9]},
    {"sumSqAmpl_h_to_anti_h_h", f_G[10]},
    {"sumSqAmpl_h_to_anti_Z_Z", f_G[11]},
    {"sumSqAmpl_h_to_anti_u_u", f_G[12]},
    {"sumSqAmpl_h_to_anti_c_c", f_G[13]},
    {"sumSqAmpl_h_to_anti_t_t", f_G[14]},
    {"sumSqAmpl_h_to_anti_d_d", f_G[15]},
    {"sumSqAmpl_h_to_anti_s_s", f_G[16]},
    {"sumSqAmpl_h_to_anti_b_b", f_G[17]},
    {"sumSqAmpl_h_to_anti_e_e", f_G[18]},
    {"sumSqAmpl_h_to_anti_mu_mu", f_G[19]},
    {"sumSqAmpl_h_to_anti_tau_tau", f_G[20]},
    {"width_h", f_G[21]},
    {"width_W", f_G[22]},
    {"width_Z", f_G[23]},
    {"width_t", f_G[24]},
    {"width_phi", f_G[25]},
};


}
 // End of namespace scalar2to2

#endif
