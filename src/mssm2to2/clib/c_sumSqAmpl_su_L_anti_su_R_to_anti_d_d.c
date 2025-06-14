#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_R_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_R_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
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
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = m_d*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = m_u*e_em*mu_h*IT_0000*IT_0005;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = pow(m_su_L, 2);
    const ccomplex_t IT_0012 = pow(m_su_R, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = m_d*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0003, -1);
    const ccomplex_t IT_0019 = m_u*e_em*mu_h*IT_0000*IT_0001*IT_0005*IT_0018;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = pow(m_Z, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0022 
      + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = -IT_0015 + -IT_0025;
    const ccomplex_t IT_0027 = pow(m_d, 2);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0000*IT_0005*IT_0018;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp2)*IT_0005;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0012 + IT_0027 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0000*IT_0005*IT_0018;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*V_ud*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0012 + IT_0027 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0034 + IT_0041;
    const ccomplex_t IT_0043 = s_13 + s_14;
    const ccomplex_t IT_0044 = m_d*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*V_ud
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = IT_0036*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0049 = IT_0039*IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0029*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0054 = IT_0032*IT_0052*IT_0053;
    const ccomplex_t IT_0055 = m_d*IT_0041;
    const ccomplex_t IT_0056 = m_d*IT_0034;
    const ccomplex_t IT_0057 = IT_0049 + IT_0054 + -IT_0055 + -IT_0056;
    const ccomplex_t IT_0058 = s_34 + IT_0027;
    const ccomplex_t IT_0059 = sin(alpha);
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0000
      *IT_0005*IT_0018*IT_0059;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0061*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = cos(alpha);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0000
      *IT_0005*IT_0018*IT_0068;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = -IT_0067 + -IT_0076;
    const ccomplex_t IT_0078 = -s_14;
    const ccomplex_t IT_0079 = s_13 + IT_0078;
    const ccomplex_t IT_0080 = m_d*IT_0079;
    const ccomplex_t IT_0081 = -IT_0027;
    const ccomplex_t IT_0082 = m_d*s_13;
    const ccomplex_t IT_0083 = s_13*s_14;
    const ccomplex_t IT_0084 = s_34*IT_0011;
    const ccomplex_t IT_0085 = (-6)*conj(IT_0026);
    const ccomplex_t IT_0086 = IT_0026*(conj(IT_0026)*(36*s_34 + 36*IT_0027) +
       (-6)*conj(IT_0042)*IT_0044 + (-6)*conj(IT_0057)*IT_0058) + IT_0077*((36
      *s_34 + (-36)*IT_0027)*conj(IT_0077) + 6*conj(IT_0042)*IT_0080 + 6*conj
      (IT_0057)*(s_34 + IT_0081)) + 6*IT_0042*(conj(IT_0077)*IT_0080 + 3*conj
      (IT_0057)*IT_0082 + 1./6*conj(IT_0042)*(36*IT_0083 + (-18)*IT_0084) + 1./6
      *IT_0044*IT_0085) + IT_0057*(18*s_34*conj(IT_0057) + 6*conj(IT_0077)*(s_34
       + IT_0081) + 18*conj(IT_0042)*IT_0082 + IT_0058*IT_0085);
    return create_ccomplex_return(IT_0086);
}

