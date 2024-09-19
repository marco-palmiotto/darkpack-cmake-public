#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -2);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = pow(m_s, 2);
    const ccomplex_t IT_0009 = cpow(IT_0004, -1);
    const ccomplex_t IT_0010 = cpow(IT_0006, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0002
      *IT_0003 + IT_0004*IT_0005)*IT_0006 + (-3)*IT_0005*IT_0007*IT_0008*IT_0009
      *IT_0010);
    const ccomplex_t IT_0012 = (-0.333333333333333)*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0003, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*m_u*e_em*IT_0002*IT_0007
      *IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = pow(m_ss_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*IT_0001*(IT_0002
      *IT_0004 + -IT_0003*IT_0005)*IT_0006 + (-3)*IT_0002*IT_0007*IT_0008
      *IT_0009*IT_0010);
    const ccomplex_t IT_0022 = 0.333333333333333*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*m_u*e_em*IT_0005*IT_0007
      *IT_0010*IT_0013;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_u, 2);
    const ccomplex_t IT_0031 = cpow(IT_0000, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0031;
    const ccomplex_t IT_0033 = (-0.666666666666667)*IT_0032;
    const ccomplex_t IT_0034 = IT_0006*IT_0031;
    const ccomplex_t IT_0035 = e_em*IT_0034;
    const ccomplex_t IT_0036 = IT_0000*IT_0010;
    const ccomplex_t IT_0037 = e_em*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0035 + (-3)*IT_0037);
    const ccomplex_t IT_0039 = (-0.166666666666667)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0033*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*V_us
      *e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*V_us
      *e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0007*IT_0009*IT_0010;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0052 + -IT_0060;
    const ccomplex_t IT_0062 = IT_0033*IT_0040;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0033*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0017*IT_0030;
    const ccomplex_t IT_0067 = m_u*IT_0033;
    const ccomplex_t IT_0068 = (-0.333333333333333)*IT_0032;
    const ccomplex_t IT_0069 = m_u*IT_0068;
    const ccomplex_t IT_0070 = -(IT_0042 + -IT_0063)*(IT_0067 + -IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = (IT_0042 + -IT_0063)*(IT_0067 + -IT_0069);
    const ccomplex_t IT_0074 = (-18)*IT_0030;
    const ccomplex_t IT_0075 = 18*s_34;
    const ccomplex_t IT_0076 = IT_0074 + IT_0075;
    const ccomplex_t IT_0077 = s_34*IT_0017;
    const ccomplex_t IT_0078 = (-6)*IT_0077;
    const ccomplex_t IT_0079 = s_13*s_14;
    const ccomplex_t IT_0080 = 12*IT_0079;
    const ccomplex_t IT_0081 = IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = (-18)*IT_0077;
    const ccomplex_t IT_0083 = 36*IT_0079;
    const ccomplex_t IT_0084 = IT_0082 + IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0086 = (-2)*IT_0085;
    const ccomplex_t IT_0087 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0086*IT_0089;
    const ccomplex_t IT_0091 = 0.25*IT_0090;
    const ccomplex_t IT_0092 = cpow(IT_0091, 2);
    const ccomplex_t IT_0093 = (-36)*IT_0077;
    const ccomplex_t IT_0094 = (-36)*IT_0066;
    const ccomplex_t IT_0095 = 72*IT_0079;
    const ccomplex_t IT_0096 = IT_0093 + IT_0094 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0098 = 0.666666666666667*IT_0097;
    const ccomplex_t IT_0099 = cpow(IT_0098, 2);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0087*IT_0100;
    const ccomplex_t IT_0102 = (-0.0833333333333333)*IT_0090 + 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0052 + IT_0060;
    const ccomplex_t IT_0104 = m_u*IT_0103;
    const ccomplex_t IT_0105 = (-18)*conj(IT_0072);
    const ccomplex_t IT_0106 = (-18)*IT_0072;
    const ccomplex_t IT_0107 = (-6)*IT_0072;
    const ccomplex_t IT_0108 = (-18)*IT_0066;
    const ccomplex_t IT_0109 = IT_0082 + IT_0083 + IT_0108;
    const ccomplex_t IT_0110 = conj(IT_0044) + conj(IT_0065);
    const ccomplex_t IT_0111 = IT_0044 + IT_0065;
    const ccomplex_t IT_0112 = -IT_0102;
    const ccomplex_t IT_0113 = -IT_0091;
    const ccomplex_t IT_0114 = (-6)*IT_0066;
    const ccomplex_t IT_0115 = IT_0078 + IT_0080 + IT_0114;
    const ccomplex_t IT_0116 = 6*IT_0029;
    const ccomplex_t IT_0117 = -IT_0030;
    const ccomplex_t IT_0118 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0119 = 6*conj(IT_0029);
    const ccomplex_t IT_0120 = 6*IT_0104;
    const ccomplex_t IT_0121 = m_u*s_13;
    const ccomplex_t IT_0122 = 18*conj(IT_0072);
    const ccomplex_t IT_0123 = 3*conj(IT_0104);
    const ccomplex_t IT_0124 = 6*conj(IT_0104);
    const ccomplex_t IT_0125 = 18*conj(IT_0073);
    const ccomplex_t IT_0126 = 18*conj(IT_0065)*IT_0072 + 6*conj(IT_0061)*
      (IT_0073 + 3*IT_0104) + conj(IT_0044)*(18*IT_0073 + IT_0120) + IT_0065
      *IT_0122 + 6*IT_0061*(conj(IT_0073) + IT_0123) + IT_0044*(IT_0124 +
       IT_0125);
    const ccomplex_t IT_0127 = m_u*s_14;
    const ccomplex_t IT_0128 = (-18)*IT_0073;
    const ccomplex_t IT_0129 = (-6)*IT_0104;
    const ccomplex_t IT_0130 = -s_14;
    const ccomplex_t IT_0131 = s_13 + IT_0130;
    const ccomplex_t IT_0132 = m_u*IT_0131;
    const ccomplex_t IT_0133 = 18*IT_0029;
    const ccomplex_t IT_0134 = 18*conj(IT_0029);
    const ccomplex_t IT_0135 = conj(IT_0061)*IT_0116 + IT_0061*IT_0119 + (-12)
      *IT_0091*(IT_0029 + -conj(IT_0029) + -1./2*conj(IT_0072) + 1./2*IT_0073 + 
      -1./2*conj(IT_0073) + 3./2*IT_0104 + (-0.0833333333333333)*IT_0107 + -1./2
      *IT_0123) + IT_0102*((-36)*IT_0029 + 36*conj(IT_0029) + IT_0106 + IT_0122 
      + IT_0124 + IT_0125 + IT_0128 + IT_0129) + IT_0110*IT_0133 + IT_0111
      *IT_0134;
    const ccomplex_t IT_0136 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + (((-18)*conj(IT_0044) + (-6)*conj(IT_0061))*IT_0065 + ((-18)*IT_0044 + 
      (-6)*IT_0061)*conj(IT_0065))*IT_0066 + IT_0072*conj(IT_0072)*IT_0075 + 
      (conj(IT_0029)*(IT_0072 + IT_0073) + IT_0029*(conj(IT_0072) + conj(IT_0073
      )))*IT_0076 + (conj(IT_0044)*IT_0061 + IT_0044*conj(IT_0061))*IT_0081 + 
      (IT_0044*conj(IT_0044) + IT_0061*conj(IT_0061) + IT_0065*conj(IT_0065))
      *IT_0084 + -IT_0092*IT_0096 + (-2)*IT_0102*(((-12)*IT_0066 + (-12)*IT_0077
       + 24*IT_0079)*IT_0091 + 0.5*IT_0096*IT_0102) + 6*IT_0073*
      (0.166666666666667*conj(IT_0073)*IT_0075 + s_34*conj(IT_0104) +
       0.166666666666667*IT_0030*IT_0105) + IT_0030*(conj(IT_0073)*IT_0106 +
       conj(IT_0104)*IT_0107) + IT_0109*(conj(IT_0061)*IT_0091 + IT_0102*IT_0110
       + IT_0111*IT_0112 + IT_0061*IT_0113) + (conj(IT_0061)*IT_0102 + IT_0091
      *IT_0110 + IT_0061*IT_0112 + IT_0111*IT_0113)*IT_0115 + conj(IT_0104)
      *IT_0116*(s_34 + IT_0117) + (s_34*conj(IT_0073) + 0.166666666666667
      *IT_0075*conj(IT_0104) + 0.166666666666667*IT_0030*IT_0118 +
       0.166666666666667*(s_34 + IT_0117)*IT_0119)*IT_0120 + IT_0121*IT_0126 +
       IT_0127*(IT_0065*((-18)*conj(IT_0073) + (-6)*conj(IT_0104)) + IT_0044
      *IT_0105 + conj(IT_0044)*IT_0106 + conj(IT_0061)*IT_0107 + IT_0061*IT_0118
       + conj(IT_0065)*(IT_0128 + IT_0129)) + IT_0132*IT_0135;
    return create_ccomplex_return(IT_0136);
}

