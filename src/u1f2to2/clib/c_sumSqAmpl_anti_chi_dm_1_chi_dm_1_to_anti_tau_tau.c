#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_V_3 = param->m_V_3;
    const creal_t m_tau = param->m_tau;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_V_3 = param->Gamma_V_3;
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const creal_t theta_12_L = param->theta_12_L;
    const creal_t theta_12_R = param->theta_12_R;
    const creal_t theta_13_L = param->theta_13_L;
    const creal_t theta_13_R = param->theta_13_R;
    const creal_t theta_23_L = param->theta_23_L;
    const creal_t theta_23_R = param->theta_23_R;
    const ccomplex_t IT_0000 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0001 = 16*s_12 + 32*IT_0000;
    const ccomplex_t IT_0002 = pow(m_tau, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = cos(theta_13_L);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = cos(theta_23_L);
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = sin(theta_12_L);
    const ccomplex_t IT_0011 = sin(theta_13_L);
    const ccomplex_t IT_0012 = cos(theta_12_L);
    const ccomplex_t IT_0013 = sin(theta_23_L);
    const ccomplex_t IT_0014 = cpow(IT_0010*IT_0011 + IT_0006*IT_0012*IT_0013,
       2);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*g_f*(IT_0007*IT_0009 + 
      -IT_0014);
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = cpow(s_12 + IT_0000 + 0.5*reg_prop + (-0.5)
      *m_V_3*(m_V_3 + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = cos(theta_13_R);
    const ccomplex_t IT_0023 = cpow(IT_0022, 2);
    const ccomplex_t IT_0024 = cos(theta_23_R);
    const ccomplex_t IT_0025 = cpow(IT_0024, 2);
    const ccomplex_t IT_0026 = sin(theta_12_R);
    const ccomplex_t IT_0027 = sin(theta_13_R);
    const ccomplex_t IT_0028 = cos(theta_12_R);
    const ccomplex_t IT_0029 = sin(theta_23_R);
    const ccomplex_t IT_0030 = cpow(IT_0026*IT_0027 + IT_0022*IT_0028*IT_0029,
       2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*g_f*(IT_0023*IT_0025 + 
      -IT_0030);
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0005*IT_0032;
    const ccomplex_t IT_0034 = IT_0018*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = s_14*s_23;
    const ccomplex_t IT_0038 = s_13*s_24;
    const ccomplex_t IT_0039 = s_34*IT_0000;
    const ccomplex_t IT_0040 = 16*IT_0037 + 16*IT_0038 + 16*IT_0039;
    return create_ccomplex_return(IT_0003*(conj(IT_0021)*IT_0036 + IT_0021
      *conj(IT_0036)) + (IT_0021*conj(IT_0021) + IT_0036*conj(IT_0036))*IT_0040);
}

