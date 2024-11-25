#ifndef CSL_LIB_su2f2to2_G_H_INCLUDED
#define CSL_LIB_su2f2to2_G_H_INCLUDED

#include <array>
#include "common.h"
#include "librarytensor.h"
#include "callable.h"
#include "csl/initSanitizer.h"
#include "params.h"
#include "func_su2f2to2.h"

namespace su2f2to2 {


inline std::array<Callable<complex_t, param_t>, 51> f_G = {
    Callable{"m_chi_dm_1", m_chi_dm_1},
    Callable{"m_chi_dm_2", m_chi_dm_2},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_phi", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_phi},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_V", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_V},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_e", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_e},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_tau_nu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_tau_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_e", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_e},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_phi", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_phi},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_e", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_e},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_tau_nu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_tau_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_e", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_e},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_mu", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_tau_tau", sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_tau_tau},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_phi", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_phi},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_V", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_V},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_e", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_e},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_mu", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_tau", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_mu", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_tau", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_tau_nu_tau", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_tau_nu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_e", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_e},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_mu", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_mu},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_mu", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_mu},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_tau", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_tau},
    Callable{"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_tau_tau", sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_tau_tau},
    Callable{"width_h", width_h},
    Callable{"width_W", width_W},
    Callable{"width_Z", width_Z},
    Callable{"width_t", width_t},
    Callable{"width_V", width_V},
    Callable{"width_V_3", width_V_3},
};

inline std::map<std::string, Callable<complex_t, param_t>> fmap_G {
    {"m_chi_dm_1", f_G[0]},
    {"m_chi_dm_2", f_G[1]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3", f_G[2]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_phi", f_G[3]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_V", f_G[4]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_e", f_G[5]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_mu", f_G[6]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_e_nu_tau", f_G[7]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_mu", f_G[8]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_mu_nu_tau", f_G[9]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_nu_tau_nu_tau", f_G[10]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_e", f_G[11]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_mu", f_G[12]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_e_tau", f_G[13]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_mu", f_G[14]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_mu_tau", f_G[15]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau", f_G[16]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_phi", f_G[17]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_e", f_G[18]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu", f_G[19]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_tau", f_G[20]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_mu", f_G[21]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_mu_nu_tau", f_G[22]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_tau_nu_tau", f_G[23]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_e", f_G[24]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_mu", f_G[25]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_e_tau", f_G[26]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_mu", f_G[27]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_mu_tau", f_G[28]},
    {"sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_tau_tau", f_G[29]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3", f_G[30]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_phi", f_G[31]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_V", f_G[32]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_e", f_G[33]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_mu", f_G[34]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_e_nu_tau", f_G[35]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_mu", f_G[36]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_mu_nu_tau", f_G[37]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_nu_tau_nu_tau", f_G[38]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_e", f_G[39]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_mu", f_G[40]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau", f_G[41]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_mu", f_G[42]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_mu_tau", f_G[43]},
    {"sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_tau_tau", f_G[44]},
    {"width_h", f_G[45]},
    {"width_W", f_G[46]},
    {"width_Z", f_G[47]},
    {"width_t", f_G[48]},
    {"width_V", f_G[49]},
    {"width_V_3", f_G[50]},
};


}
 // End of namespace su2f2to2

#endif
