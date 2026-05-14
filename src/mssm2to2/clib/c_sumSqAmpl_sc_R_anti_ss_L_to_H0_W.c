#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_ss_L_to_H0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_ss_L_to_H0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = pow(m_sc_R, 2);
    const ccomplex_t IT_0011 = pow(m_ss_L, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0001, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*mu_h*IT_0002*IT_0005*IT_0007*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + IT_0021 
      + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0017*IT_0024;
    const ccomplex_t IT_0026 = IT_0015 + IT_0025;
    const ccomplex_t IT_0027 = pow(m_H0, 2);
    const ccomplex_t IT_0028 = pow(m_W, -2);
    const ccomplex_t IT_0029 = pow(s_34, 2);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0005;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0000
      *IT_0005*IT_0007*IT_0018;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_13 + IT_0010 + IT_0027 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = s_24*s_34*IT_0028;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = s_23 + IT_0040;
    const ccomplex_t IT_0042 = pow(s_24, 2);
    const ccomplex_t IT_0043 = IT_0028*IT_0042;
    return create_ccomplex_return((-3)*IT_0026*(conj(IT_0026)*(IT_0027 + 
      -IT_0030) + conj(IT_0038)*IT_0041) + (-3)*IT_0038*(conj(IT_0026)*IT_0041 +
       conj(IT_0038)*(IT_0011 + -IT_0043)));
}

