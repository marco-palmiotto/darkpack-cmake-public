#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_e_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
    const creal_t m_e = param->m_e;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_tau = param->m_tau;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_V_3 = param->Gamma_V_3;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const creal_t theta_12_L = param->theta_12_L;
    const creal_t theta_12_R = param->theta_12_R;
    const creal_t theta_13_L = param->theta_13_L;
    const creal_t theta_13_R = param->theta_13_R;
    const creal_t theta_23_L = param->theta_23_L;
    const creal_t theta_23_R = param->theta_23_R;
    const ccomplex_t IT_0000 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + (-2)*IT_0000 + -reg_prop + m_V
      *(m_V + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0003 = (-0.5)*IT_0002;
    const ccomplex_t IT_0004 = cos(theta_13_R);
    const ccomplex_t IT_0005 = cos(theta_23_R);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = sin(theta_13_R);
    const ccomplex_t IT_0008 = sin(theta_12_R);
    const ccomplex_t IT_0009 = cos(theta_12_R);
    const ccomplex_t IT_0010 = sin(theta_23_R);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*g_f*(IT_0004*IT_0006*IT_0007
       + (IT_0007*IT_0008 + IT_0004*IT_0009*IT_0010)*(IT_0004*IT_0008 + -IT_0007
      *IT_0009*IT_0010));
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = IT_0001*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = s_14*s_23;
    const ccomplex_t IT_0018 = s_13*s_24;
    const ccomplex_t IT_0019 = s_34*IT_0000;
    const ccomplex_t IT_0020 = 16*IT_0017 + 16*IT_0018 + 16*IT_0019;
    const ccomplex_t IT_0021 = cos(theta_13_L);
    const ccomplex_t IT_0022 = cos(theta_23_L);
    const ccomplex_t IT_0023 = cpow(IT_0022, 2);
    const ccomplex_t IT_0024 = sin(theta_13_L);
    const ccomplex_t IT_0025 = sin(theta_12_L);
    const ccomplex_t IT_0026 = cos(theta_12_L);
    const ccomplex_t IT_0027 = sin(theta_23_L);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*g_f*(IT_0021*IT_0023*IT_0024
       + (IT_0024*IT_0025 + IT_0021*IT_0026*IT_0027)*(IT_0021*IT_0025 + -IT_0024
      *IT_0026*IT_0027));
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0003*IT_0029;
    const ccomplex_t IT_0031 = IT_0001*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (-16)*s_12 + (-32)*IT_0000;
    const ccomplex_t IT_0035 = m_e*m_tau;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = -IT_0036;
    return create_ccomplex_return(IT_0033*(IT_0020*conj(IT_0033) + conj
      (IT_0016)*IT_0037) + IT_0016*(conj(IT_0016)*IT_0020 + conj(IT_0033)
      *IT_0037));
}

