#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_R_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_R_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = m_t*e_em*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = m_s*e_em*mu_h*IT_0000;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = IT_0002*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0010 + (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = pow(m_Z, 2);
    const ccomplex_t IT_0017 = pow(m_ss_L, 2);
    const ccomplex_t IT_0018 = pow(m_ss_R, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + IT_0016 + -IT_0017 + -IT_0018 
      + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0002*IT_0005;
    const ccomplex_t IT_0023 = IT_0008*IT_0011;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0022 + 0.5*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0008, -1);
    const ccomplex_t IT_0027 = m_t*e_em*IT_0000*IT_0002*IT_0006*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0017 + -IT_0018 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0021 + -IT_0032;
    const ccomplex_t IT_0034 = pow(m_t, 2);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm1)*IT_0002;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0007;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0018 + IT_0034 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0002;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0007;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0018 + IT_0034 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = -IT_0041 + -IT_0048;
    const ccomplex_t IT_0050 = s_13 + s_14;
    const ccomplex_t IT_0051 = m_t*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u1
      *e_em*IT_0000*IT_0002*IT_0026;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0037*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0056 = IT_0039*IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0000*IT_0002*IT_0026;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0044*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0061 = IT_0046*IT_0059*IT_0060;
    const ccomplex_t IT_0062 = m_t*IT_0041;
    const ccomplex_t IT_0063 = m_t*IT_0048;
    const ccomplex_t IT_0064 = -IT_0056 + -IT_0061 + IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = s_34 + IT_0034;
    const ccomplex_t IT_0066 = sin(alpha);
    const ccomplex_t IT_0067 = cos(alpha);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0066 + m_s*e_em*mu_h*IT_0000*IT_0002*IT_0007*IT_0067);
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_t*e_em*IT_0000*IT_0002
      *IT_0026*IT_0067;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0017 + -IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*-1)*(m_s*e_em*mu_h*IT_0000
      *IT_0002*IT_0007*IT_0066 + (-1.4142135623731)*A_b*V_ts*conj(V_ts)*IT_0067);
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_t*e_em*IT_0000*IT_0002
      *IT_0026*IT_0066;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_12 + -IT_0017 + -IT_0018 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0075 + -IT_0083;
    const ccomplex_t IT_0085 = -s_14;
    const ccomplex_t IT_0086 = s_13 + IT_0085;
    const ccomplex_t IT_0087 = m_t*IT_0086;
    const ccomplex_t IT_0088 = -IT_0034;
    const ccomplex_t IT_0089 = m_t*s_13;
    const ccomplex_t IT_0090 = s_13*s_14;
    const ccomplex_t IT_0091 = s_34*IT_0017;
    const ccomplex_t IT_0092 = (-6)*conj(IT_0033);
    const ccomplex_t IT_0093 = IT_0033*(conj(IT_0033)*(36*s_34 + 36*IT_0034) +
       (-6)*conj(IT_0049)*IT_0051 + (-6)*conj(IT_0064)*IT_0065) + IT_0084*((36
      *s_34 + (-36)*IT_0034)*conj(IT_0084) + 6*conj(IT_0049)*IT_0087 + 6*conj
      (IT_0064)*(s_34 + IT_0088)) + 6*IT_0049*(conj(IT_0084)*IT_0087 + 3*conj
      (IT_0064)*IT_0089 + 1./6*conj(IT_0049)*(36*IT_0090 + (-18)*IT_0091) + 1./6
      *IT_0051*IT_0092) + IT_0064*(18*s_34*conj(IT_0064) + 6*conj(IT_0084)*(s_34
       + IT_0088) + 18*conj(IT_0049)*IT_0089 + IT_0065*IT_0092);
    return create_ccomplex_return(IT_0093);
}

