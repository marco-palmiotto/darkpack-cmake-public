#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_R_anti_sc_R_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_R_anti_sc_R_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
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
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = pow(m_c, 2);
    const ccomplex_t IT_0013 = cpow(IT_0010, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0001*IT_0003 + -IT_0010*IT_0011) + 1.5*IT_0000*IT_0005*IT_0011*IT_0012
      *IT_0013);
    const ccomplex_t IT_0015 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = pow(m_sc_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0011;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*(m_W*IT_0004*IT_0009*
      (IT_0003*IT_0010 + IT_0001*IT_0011) + (-1.5)*IT_0000*IT_0003*IT_0005
      *IT_0012*IT_0013);
    const ccomplex_t IT_0024 = 0.666666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_d, 2);
    const ccomplex_t IT_0031 = cpow(IT_0008, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0031;
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = m_d*IT_0033;
    const ccomplex_t IT_0035 = 1.33333333333333*IT_0032;
    const ccomplex_t IT_0036 = m_d*IT_0035;
    const ccomplex_t IT_0037 = IT_0004*IT_0031;
    const ccomplex_t IT_0038 = e_em*IT_0037;
    const ccomplex_t IT_0039 = IT_0005*IT_0008;
    const ccomplex_t IT_0040 = e_em*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0038 + 3*IT_0040);
    const ccomplex_t IT_0042 = (-0.166666666666667)*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = 0.333333333333333*IT_0032;
    const ccomplex_t IT_0047 = IT_0043*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (IT_0034 + -IT_0036)*(IT_0045 + -IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = (-18)*IT_0030;
    const ccomplex_t IT_0052 = 18*s_34;
    const ccomplex_t IT_0053 = IT_0051 + IT_0052;
    const ccomplex_t IT_0054 = IT_0035*IT_0045;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u1
      *e_em*IT_0000*IT_0005*IT_0013;
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u1)*e_em*IT_0000*IT_0005*IT_0013;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0017 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cd)
      *conj(V_u2)*e_em*IT_0000*IT_0005*IT_0013;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*V_u2
      *e_em*IT_0000*IT_0005*IT_0013;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = IT_0064*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0063 + IT_0071;
    const ccomplex_t IT_0073 = s_34*IT_0017;
    const ccomplex_t IT_0074 = (-6)*IT_0073;
    const ccomplex_t IT_0075 = s_13*s_14;
    const ccomplex_t IT_0076 = 12*IT_0075;
    const ccomplex_t IT_0077 = IT_0074 + IT_0076;
    const ccomplex_t IT_0078 = IT_0035*IT_0048;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = (-18)*IT_0073;
    const ccomplex_t IT_0081 = 36*IT_0075;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0084 = (-2)*IT_0083;
    const ccomplex_t IT_0085 = cpow(s_12 + IT_0017 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0083*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0084*IT_0087;
    const ccomplex_t IT_0089 = 0.25*IT_0088;
    const ccomplex_t IT_0090 = cpow(IT_0089, 2);
    const ccomplex_t IT_0091 = (-36)*IT_0073;
    const ccomplex_t IT_0092 = IT_0017*IT_0030;
    const ccomplex_t IT_0093 = (-36)*IT_0092;
    const ccomplex_t IT_0094 = 72*IT_0075;
    const ccomplex_t IT_0095 = IT_0091 + IT_0093 + IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0097 = (-0.333333333333333)*IT_0096;
    const ccomplex_t IT_0098 = (-1.33333333333333)*IT_0096;
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = IT_0085*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = (-0.0833333333333333)*IT_0088 + 0.5*IT_0101;
    const ccomplex_t IT_0103 = (-18)*IT_0050;
    const ccomplex_t IT_0104 = m_d*IT_0072;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = (-6)*IT_0105;
    const ccomplex_t IT_0107 = IT_0103 + IT_0106;
    const ccomplex_t IT_0108 = (-18)*conj(IT_0050);
    const ccomplex_t IT_0109 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = (-18)*IT_0055;
    const ccomplex_t IT_0112 = (-6)*IT_0072;
    const ccomplex_t IT_0113 = IT_0111 + IT_0112;
    const ccomplex_t IT_0114 = (-18)*conj(IT_0055);
    const ccomplex_t IT_0115 = (-6)*conj(IT_0072);
    const ccomplex_t IT_0116 = IT_0114 + IT_0115;
    const ccomplex_t IT_0117 = m_d*s_14;
    const ccomplex_t IT_0118 = -s_14;
    const ccomplex_t IT_0119 = s_13 + IT_0118;
    const ccomplex_t IT_0120 = m_d*IT_0119;
    const ccomplex_t IT_0121 = 6*IT_0029;
    const ccomplex_t IT_0122 = 6*conj(IT_0029);
    const ccomplex_t IT_0123 = 6*conj(IT_0050);
    const ccomplex_t IT_0124 = 18*conj(IT_0105);
    const ccomplex_t IT_0125 = conj(IT_0055) + conj(IT_0079);
    const ccomplex_t IT_0126 = 18*IT_0029;
    const ccomplex_t IT_0127 = IT_0055 + IT_0079;
    const ccomplex_t IT_0128 = 18*conj(IT_0029);
    const ccomplex_t IT_0129 = 6*conj(IT_0105);
    const ccomplex_t IT_0130 = 18*conj(IT_0050);
    const ccomplex_t IT_0131 = conj(IT_0072)*IT_0121 + IT_0072*IT_0122 + (-12)
      *IT_0089*(IT_0029 + -conj(IT_0029) + 1./2*IT_0049 + -1./2*conj(IT_0049) +
       1./2*IT_0050 + 3./2*IT_0105 + (-0.0833333333333333)*IT_0123 + (
      -0.0833333333333333)*IT_0124) + IT_0125*IT_0126 + IT_0127*IT_0128 +
       IT_0102*((-36)*IT_0029 + 36*conj(IT_0029) + (-18)*IT_0049 + 18*conj
      (IT_0049) + IT_0103 + IT_0106 + IT_0129 + IT_0130);
    const ccomplex_t IT_0132 = (-18)*IT_0092;
    const ccomplex_t IT_0133 = IT_0080 + IT_0081 + IT_0132;
    const ccomplex_t IT_0134 = -IT_0102;
    const ccomplex_t IT_0135 = -IT_0089;
    const ccomplex_t IT_0136 = (-6)*IT_0092;
    const ccomplex_t IT_0137 = IT_0074 + IT_0076 + IT_0136;
    const ccomplex_t IT_0138 = -IT_0030;
    const ccomplex_t IT_0139 = 6*IT_0105;
    const ccomplex_t IT_0140 = m_d*s_13;
    const ccomplex_t IT_0141 = IT_0050*(18*conj(IT_0055) + 6*conj(IT_0072)) +
       18*conj(IT_0049)*IT_0079 + 18*IT_0049*conj(IT_0079) + IT_0072*(IT_0123 +
       IT_0124) + IT_0055*(IT_0129 + IT_0130) + (conj(IT_0055) + 3*conj(IT_0072)
      )*IT_0139;
    const ccomplex_t IT_0142 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + IT_0049*conj(IT_0049)*IT_0052 + (conj(IT_0029)*(IT_0049 + IT_0050) +
       IT_0029*(conj(IT_0049) + conj(IT_0050)))*IT_0053 + (conj(IT_0055)*IT_0072
       + IT_0055*conj(IT_0072))*IT_0077 + (IT_0055*conj(IT_0055) + IT_0072*conj
      (IT_0072) + IT_0079*conj(IT_0079))*IT_0082 + -IT_0090*IT_0095 + (-2)
      *IT_0102*(IT_0089*((-12)*IT_0073 + 24*IT_0075 + (-12)*IT_0092) + 0.5
      *IT_0095*IT_0102) + 6*IT_0050*(0.166666666666667*conj(IT_0050)*IT_0052 +
       s_34*conj(IT_0105)) + IT_0030*(conj(IT_0049)*IT_0107 + IT_0049*IT_0110) +
       IT_0092*(conj(IT_0079)*IT_0113 + IT_0079*IT_0116) + (conj(IT_0079)
      *IT_0107 + IT_0079*IT_0110 + conj(IT_0049)*IT_0113 + IT_0049*IT_0116)
      *IT_0117 + IT_0120*IT_0131 + IT_0133*(conj(IT_0072)*IT_0089 + IT_0102
      *IT_0125 + IT_0127*IT_0134 + IT_0072*IT_0135) + (conj(IT_0072)*IT_0102 +
       IT_0089*IT_0125 + IT_0072*IT_0134 + IT_0127*IT_0135)*IT_0137 + conj
      (IT_0105)*IT_0121*(s_34 + IT_0138) + (s_34*conj(IT_0050) +
       0.166666666666667*IT_0052*conj(IT_0105) + 0.166666666666667*IT_0122*(s_34
       + IT_0138))*IT_0139 + IT_0140*IT_0141;
    return create_ccomplex_return(IT_0142);
}

