#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_ss_R_to_s_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_ss_R_to_s_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_s = param->m_s;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0002 = (-0.333333333333333)*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = pow(m_s, 2);
    const ccomplex_t IT_0005 = pow(m_sG, 2);
    const ccomplex_t IT_0006 = pow(m_ss_R, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_12 + IT_0004 + -IT_0005 + -IT_0006 
      + -reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = 0.666666666666667*IT_0001;
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0004 + IT_0005 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (-2)*IT_0009 + -IT_0015;
    const ccomplex_t IT_0017 = s_14*IT_0004;
    const ccomplex_t IT_0018 = s_14*s_34;
    const ccomplex_t IT_0019 = s_34*IT_0005;
    const ccomplex_t IT_0020 = s_13*s_14;
    const ccomplex_t IT_0021 = (-16)*IT_0020;
    const ccomplex_t IT_0022 = 8*IT_0019 + IT_0021;
    const ccomplex_t IT_0023 = 0.333333333333333*IT_0001;
    const ccomplex_t IT_0024 = IT_0014*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (-8)*s_13;
    const ccomplex_t IT_0027 = s_34*IT_0026;
    const ccomplex_t IT_0028 = pow(s_13, 2);
    const ccomplex_t IT_0029 = (-8)*IT_0028;
    const ccomplex_t IT_0030 = 8*IT_0010*((IT_0016 + -conj(IT_0016))*IT_0017 +
       (-4)*IT_0010*IT_0018 + (-0.125)*(IT_0015 + -conj(IT_0015))*IT_0022) +
       IT_0025*(s_14*conj(IT_0015)*IT_0026 + conj(IT_0016)*IT_0027) + (-8)
      *IT_0016*(s_13*IT_0004*conj(IT_0016) + (-0.125)*conj(IT_0025)*IT_0027 + (
      -0.125)*conj(IT_0015)*IT_0029) + (-8)*IT_0015*(s_13*IT_0005*conj(IT_0015) 
      + (-0.125)*s_14*conj(IT_0025)*IT_0026 + (-0.125)*conj(IT_0016)*IT_0029);
    return create_ccomplex_return(IT_0030);
}

