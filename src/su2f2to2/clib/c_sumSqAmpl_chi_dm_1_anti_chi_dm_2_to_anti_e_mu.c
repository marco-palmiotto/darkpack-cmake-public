#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_1_anti_chi_dm_2_to_anti_e_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_1_anti_chi_dm_2_to_anti_e_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
    const creal_t m_e = param->m_e;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t Gamma_V = param->Gamma_V;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const creal_t theta_12_L = param->theta_12_L;
    const creal_t theta_12_R = param->theta_12_R;
    const creal_t theta_13_L = param->theta_13_L;
    const creal_t theta_13_R = param->theta_13_R;
    const creal_t theta_23_L = param->theta_23_L;
    const creal_t theta_23_R = param->theta_23_R;
    const ccomplex_t IT_0000 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0001 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_12 + m_V*(m_V + (0 + _Complex_I*-1)
      *Gamma_V) + -IT_0000 + -IT_0001 + -reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_12_R);
    const ccomplex_t IT_0006 = cos(theta_23_R);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = sin(theta_13_R);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0005
      *IT_0007*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = s_13*s_24;
    const ccomplex_t IT_0016 = s_34*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0017 = s_14*s_23;
    const ccomplex_t IT_0018 = cos(theta_12_L);
    const ccomplex_t IT_0019 = cos(theta_23_L);
    const ccomplex_t IT_0020 = cpow(IT_0019, 2);
    const ccomplex_t IT_0021 = sin(theta_13_L);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0018
      *IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = IT_0002*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = m_e*m_mu*s_12;
    const ccomplex_t IT_0029 = m_e*m_mu*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0030 = 32*IT_0029;
    const ccomplex_t IT_0031 = 16*IT_0028 + IT_0030;
    return create_ccomplex_return(IT_0027*((16*IT_0015 + 16*IT_0016 + 16
      *IT_0017)*conj(IT_0027) + conj(IT_0014)*IT_0031) + IT_0014*(conj(IT_0014)*
      (16*IT_0015 + 16*IT_0016 + 16*IT_0017) + conj(IT_0027)*IT_0031));
}

