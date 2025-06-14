#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_R_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_R_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_ss_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm1)*IT_0004;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*V_us
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm2)*IT_0004;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*V_us
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = -IT_0013 + -IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_us
      *e_em*IT_0004*IT_0006*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = IT_0002*IT_0022*IT_0027;
    const ccomplex_t IT_0029 = m_u*IT_0013;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0004*IT_0006*IT_0024;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = IT_0014*IT_0030*IT_0033;
    const ccomplex_t IT_0035 = m_u*IT_0020;
    const ccomplex_t IT_0036 = -IT_0028 + IT_0029 + -IT_0034 + IT_0035;
    const ccomplex_t IT_0037 = m_u*s_13;
    const ccomplex_t IT_0038 = m_u*e_em*IT_0004*IT_0006;
    const ccomplex_t IT_0039 = (-0.5)*IT_0038;
    const ccomplex_t IT_0040 = m_s*e_em*mu_h*IT_0006;
    const ccomplex_t IT_0041 = IT_0004*IT_0008*IT_0023;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0044 = IT_0007*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = 0.5*IT_0042 + (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0039*IT_0046;
    const ccomplex_t IT_0048 = pow(m_Z, 2);
    const ccomplex_t IT_0049 = pow(m_ss_L, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + -IT_0001 + IT_0048 + -IT_0049 
      + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0004*IT_0040;
    const ccomplex_t IT_0054 = IT_0023*IT_0043;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0053 + 0.5*IT_0055;
    const ccomplex_t IT_0057 = m_u*e_em*IT_0004*IT_0006*IT_0007*IT_0024;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + -IT_0001 + -IT_0049 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = -IT_0052 + -IT_0062;
    const ccomplex_t IT_0064 = s_13 + s_14;
    const ccomplex_t IT_0065 = m_u*IT_0064;
    const ccomplex_t IT_0066 = s_13*s_14;
    const ccomplex_t IT_0067 = s_34*IT_0049;
    const ccomplex_t IT_0068 = sin(alpha);
    const ccomplex_t IT_0069 = cos(alpha);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0068 + m_s*e_em*mu_h*IT_0004*IT_0006*IT_0008*IT_0069);
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_u*e_em*IT_0004*IT_0006
      *IT_0024*IT_0069;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = IT_0071*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0001 + -IT_0049 + -reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*-1)*(m_s*e_em*mu_h*IT_0004
      *IT_0006*IT_0008*IT_0068 + (-1.4142135623731)*A_b*V_ts*conj(V_ts)*IT_0069);
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*m_u*e_em*IT_0004*IT_0006
      *IT_0024*IT_0068;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + -IT_0001 + -IT_0049 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0077 + -IT_0085;
    const ccomplex_t IT_0087 = -s_14;
    const ccomplex_t IT_0088 = s_13 + IT_0087;
    const ccomplex_t IT_0089 = m_u*IT_0088;
    const ccomplex_t IT_0090 = s_34 + IT_0000;
    const ccomplex_t IT_0091 = -IT_0000;
    const ccomplex_t IT_0092 = IT_0021*(18*conj(IT_0036)*IT_0037 + (-6)*conj
      (IT_0063)*IT_0065 + conj(IT_0021)*(36*IT_0066 + (-18)*IT_0067) + 6*conj
      (IT_0086)*IT_0089) + IT_0063*((36*s_34 + 36*IT_0000)*conj(IT_0063) + (-6)
      *conj(IT_0021)*IT_0065 + (-6)*conj(IT_0036)*IT_0090) + IT_0086*((36*s_34 +
       (-36)*IT_0000)*conj(IT_0086) + 6*conj(IT_0021)*IT_0089 + 6*conj(IT_0036)*
      (s_34 + IT_0091)) + IT_0036*(18*s_34*conj(IT_0036) + 18*conj(IT_0021)
      *IT_0037 + (-6)*conj(IT_0063)*IT_0090 + 6*conj(IT_0086)*(s_34 + IT_0091));
    return create_ccomplex_return(IT_0092);
}

