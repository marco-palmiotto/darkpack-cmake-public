#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sc_R_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sc_R_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + IT_0004*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = cpow(IT_0004, -1);
    const ccomplex_t IT_0012 = m_c*e_em*mu_h*IT_0002*IT_0007*IT_0010*IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_Z, 2);
    const ccomplex_t IT_0015 = pow(m_sc_L, 2);
    const ccomplex_t IT_0016 = pow(m_sc_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + IT_0014 + -IT_0015 + -IT_0016 
      + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0009*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0006 + (-3)
      *IT_0000*IT_0007);
    const ccomplex_t IT_0022 = 0.333333333333333*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0003
      *IT_0007*IT_0010*IT_0011;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = pow(m_H0, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_23 + IT_0016 + IT_0025 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    const ccomplex_t IT_0030 = e_em*IT_0001*(IT_0003*IT_0004 + -IT_0002
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0031 = m_c*e_em*mu_h*IT_0007*IT_0010;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = IT_0020 + -IT_0029 + IT_0036;
    const ccomplex_t IT_0038 = pow(m_Z, -2);
    const ccomplex_t IT_0039 = pow(s_34, 2);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0006;
    const ccomplex_t IT_0042 = 1.33333333333333*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0015 + IT_0025 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0024*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0042*IT_0045;
    const ccomplex_t IT_0047 = IT_0029 + IT_0046;
    const ccomplex_t IT_0048 = s_24*s_34*IT_0038;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = s_23 + IT_0049;
    const ccomplex_t IT_0051 = pow(s_24, 2);
    const ccomplex_t IT_0052 = IT_0038*IT_0051;
    return create_ccomplex_return((-3)*IT_0037*(conj(IT_0037)*(IT_0025 + 
      -IT_0040) + conj(IT_0047)*IT_0050) + (-3)*IT_0047*(conj(IT_0037)*IT_0050 +
       conj(IT_0047)*(IT_0016 + -IT_0052)));
}

