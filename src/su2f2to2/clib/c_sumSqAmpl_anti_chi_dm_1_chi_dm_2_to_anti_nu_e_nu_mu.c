#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_nu_e_nu_mu(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
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
    const creal_t theta_13_L = param->theta_13_L;
    const creal_t theta_23_L = param->theta_23_L;
    const ccomplex_t IT_0000 = s_14*s_23;
    const ccomplex_t IT_0001 = s_34*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0002 = s_13*s_24;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_23_L);
    const ccomplex_t IT_0006 = cos(theta_13_L);
    const ccomplex_t IT_0007 = sin(theta_12_L);
    const ccomplex_t IT_0008 = cos(theta_12_L);
    const ccomplex_t IT_0009 = sin(theta_13_L);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_f*IT_0005*
      (IT_0006*IT_0007 + -IT_0005*IT_0008*IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0014 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + m_V*(m_V + (0 + _Complex_I*-1)
      *Gamma_V) + -IT_0013 + -IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    return create_ccomplex_return(16*(IT_0000 + IT_0001 + IT_0002)*IT_0018
      *conj(IT_0018));
}

