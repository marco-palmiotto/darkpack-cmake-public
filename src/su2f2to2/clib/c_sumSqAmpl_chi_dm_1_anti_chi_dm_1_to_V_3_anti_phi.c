#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_1_anti_chi_dm_1_to_V_3_anti_phi.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_1_anti_chi_dm_1_to_V_3_anti_phi(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_V_3 = param->Gamma_V_3;
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_f*m_V;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0002 = 0.5*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_12 + (-2)*IT_0004 + -reg_prop + m_V
      *(m_V + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_V, -2);
    const ccomplex_t IT_0009 = pow(m_V, 2);
    const ccomplex_t IT_0010 = s_12*IT_0009;
    const ccomplex_t IT_0011 = s_13*s_23;
    const ccomplex_t IT_0012 = (-2)*IT_0011;
    const ccomplex_t IT_0013 = IT_0010 + IT_0012;
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    return create_ccomplex_return(8*IT_0007*conj(IT_0007)*(s_12 + 1.5*IT_0004 
      + (-0.5)*IT_0014));
}

