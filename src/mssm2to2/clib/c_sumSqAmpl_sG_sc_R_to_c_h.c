#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sc_R_to_c_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sc_R_to_c_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_c, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_sc_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = cos(alpha);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*m_c*e_em*IT_0005*IT_0006
      *IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = m_c*IT_0015;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -2);
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*((-1.5)*IT_0000*IT_0005
      *IT_0006*IT_0008*IT_0010 + m_W*IT_0009*IT_0019*(IT_0006*IT_0007 + IT_0020
      *IT_0021));
    const ccomplex_t IT_0023 = 0.666666666666667*IT_0022;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (-2)*IT_0017 + IT_0027;
    const ccomplex_t IT_0029 = m_c*s_14;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0005
      *IT_0008*IT_0010*IT_0021;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = m_sG*s_34;
    const ccomplex_t IT_0038 = pow(m_h, 2);
    const ccomplex_t IT_0039 = s_13*IT_0038;
    const ccomplex_t IT_0040 = s_14*s_34;
    const ccomplex_t IT_0041 = m_c*m_sG;
    const ccomplex_t IT_0042 = 8*s_13;
    const ccomplex_t IT_0043 = (-8)*IT_0016*((IT_0028 + -conj(IT_0028))
      *IT_0029 + (IT_0036 + -conj(IT_0036))*IT_0037 + -IT_0016*(IT_0039 + (-2)
      *IT_0040)) + 8*conj(IT_0028)*(IT_0036*IT_0041 + 0.125*IT_0028*IT_0042) + 8
      *conj(IT_0036)*(IT_0028*IT_0041 + 0.125*IT_0036*IT_0042);
    return create_ccomplex_return(IT_0043);
}

