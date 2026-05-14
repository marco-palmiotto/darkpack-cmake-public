#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_ss_R_to_c_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_ss_R_to_c_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_ss_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *IT_0006*IT_0007*IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0005*IT_0013;
    const ccomplex_t IT_0015 = IT_0003*IT_0004*IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0007, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_c*m_s*V_cs
      *e_em*IT_0006*IT_0009*IT_0011*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = IT_0005*IT_0018;
    const ccomplex_t IT_0020 = pow(m_c, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_13 + IT_0001 + IT_0020 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = -IT_0015 + IT_0023;
    const ccomplex_t IT_0025 = 8*s_13;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *mu_h*IT_0006*IT_0011;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = IT_0005*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0001 + IT_0020 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*m_s*V_cs*e_em
      *IT_0006*IT_0008*IT_0011*IT_0016;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0005*IT_0033;
    const ccomplex_t IT_0035 = IT_0003*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = m_sG*IT_0036;
    const ccomplex_t IT_0038 = -IT_0031 + -IT_0037;
    const ccomplex_t IT_0039 = -IT_0036;
    const ccomplex_t IT_0040 = m_sG*s_23;
    const ccomplex_t IT_0041 = m_c*m_sG;
    const ccomplex_t IT_0042 = 8*IT_0024;
    const ccomplex_t IT_0043 = s_12*s_23;
    const ccomplex_t IT_0044 = s_13*IT_0002;
    const ccomplex_t IT_0045 = m_c*s_12;
    const ccomplex_t IT_0046 = 8*conj(IT_0024);
    const ccomplex_t IT_0047 = IT_0024*conj(IT_0024)*IT_0025 + 8*conj(IT_0038)
      *(0.125*IT_0025*IT_0038 + IT_0039*IT_0040 + 0.125*IT_0041*IT_0042) + 8
      *conj(IT_0039)*(IT_0038*IT_0040 + 0.125*IT_0039*(16*IT_0043 + (-8)*IT_0044
      ) + 0.125*IT_0042*IT_0045) + (IT_0038*IT_0041 + IT_0039*IT_0045)*IT_0046;
    return create_ccomplex_return(IT_0047);
}

