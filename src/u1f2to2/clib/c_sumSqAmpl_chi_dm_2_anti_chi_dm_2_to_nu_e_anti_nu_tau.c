#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_2_anti_chi_dm_2_to_nu_e_anti_nu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_2_anti_chi_dm_2_to_nu_e_anti_nu_tau(
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
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_V_3 = param->Gamma_V_3;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const creal_t theta_12_L = param->theta_12_L;
    const creal_t theta_13_L = param->theta_13_L;
    const creal_t theta_23_L = param->theta_23_L;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = s_13*s_24;
    const ccomplex_t IT_0002 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0003 = s_34*IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0005 = (-0.5)*IT_0004;
    const ccomplex_t IT_0006 = cos(theta_13_L);
    const ccomplex_t IT_0007 = cos(theta_23_L);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = sin(theta_13_L);
    const ccomplex_t IT_0010 = sin(theta_12_L);
    const ccomplex_t IT_0011 = cos(theta_12_L);
    const ccomplex_t IT_0012 = sin(theta_23_L);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*g_f*(IT_0006*IT_0008*IT_0009
       + (IT_0009*IT_0010 + IT_0006*IT_0011*IT_0012)*(IT_0006*IT_0010 + -IT_0009
      *IT_0011*IT_0012));
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0005*IT_0014;
    const ccomplex_t IT_0016 = cpow(s_12 + IT_0002 + 0.5*reg_prop + (-0.5)
      *m_V_3*(m_V_3 + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    return create_ccomplex_return((16*IT_0000 + 16*IT_0001 + 16*IT_0003)
      *IT_0019*conj(IT_0019));
}

