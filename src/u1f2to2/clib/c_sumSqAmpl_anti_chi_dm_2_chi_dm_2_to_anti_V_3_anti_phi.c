#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_anti_phi.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_anti_phi(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_V_3 = param->m_V_3;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_V_3 = param->Gamma_V_3;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const ccomplex_t IT_0000 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0001 = pow(m_V_3, -2);
    const ccomplex_t IT_0002 = s_13*s_23;
    const ccomplex_t IT_0003 = pow(m_V_3, 2);
    const ccomplex_t IT_0004 = s_12*IT_0003;
    const ccomplex_t IT_0005 = (-0.5)*IT_0004;
    const ccomplex_t IT_0006 = IT_0002 + IT_0005;
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*g_f*m_V_3;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = cpow(s_12 + IT_0000 + 0.5*reg_prop + (-0.5)
      *m_V_3*(m_V_3 + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    return create_ccomplex_return(8*(s_12 + 1.5*IT_0000 + IT_0007)*IT_0015
      *conj(IT_0015));
}

