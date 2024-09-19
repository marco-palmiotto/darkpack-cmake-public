#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_su_R_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_su_R_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_h = param->m_h;
    const creal_t m_u = param->m_u;
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
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t IT_0000 = pow(m_su_R, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = pow(m_u, 2);
    const ccomplex_t IT_0011 = cpow(IT_0006, -1);
    const ccomplex_t IT_0012 = cpow(IT_0008, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0004
      *IT_0005 + -IT_0006*IT_0007)*IT_0008 + 1.5*IT_0007*IT_0009*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0014 = (-0.666666666666667)*IT_0013;
    const ccomplex_t IT_0015 = cpow(IT_0004, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*m_b*e_em*IT_0005*IT_0009
      *IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_b*e_em*IT_0007*IT_0009
      *IT_0012*IT_0015;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0005
      *IT_0006 + IT_0004*IT_0007)*IT_0008 + (-1.5)*IT_0005*IT_0009*IT_0010
      *IT_0011*IT_0012);
    const ccomplex_t IT_0024 = 0.666666666666667*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0020 + -IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = cpow(IT_0002, -1);
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0031;
    const ccomplex_t IT_0033 = 1.33333333333333*IT_0032;
    const ccomplex_t IT_0034 = m_b*IT_0033;
    const ccomplex_t IT_0035 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0036 = m_b*IT_0035;
    const ccomplex_t IT_0037 = 0.333333333333333*IT_0032;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0008*IT_0031;
    const ccomplex_t IT_0042 = e_em*IT_0041;
    const ccomplex_t IT_0043 = IT_0002*IT_0012;
    const ccomplex_t IT_0044 = e_em*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0042 + 3*IT_0044);
    const ccomplex_t IT_0046 = (-0.166666666666667)*IT_0045;
    const ccomplex_t IT_0047 = IT_0038*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (IT_0034 + -IT_0036)*(IT_0040 + -IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = (-18)*IT_0030;
    const ccomplex_t IT_0054 = 18*s_34;
    const ccomplex_t IT_0055 = IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = IT_0033*IT_0048;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0009*IT_0011*IT_0012*IT_0058*V_ub_mod;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0009*IT_0011*IT_0012*IT_0061*V_ub_mod;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0060*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0009*IT_0011*IT_0012*IT_0058*V_ub_mod;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *e_em*IT_0009*IT_0011*IT_0012*IT_0061*V_ub_mod;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0000 + IT_0030 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = IT_0067 + IT_0075;
    const ccomplex_t IT_0077 = s_34*IT_0000;
    const ccomplex_t IT_0078 = (-6)*IT_0077;
    const ccomplex_t IT_0079 = s_13*s_14;
    const ccomplex_t IT_0080 = 12*IT_0079;
    const ccomplex_t IT_0081 = IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = IT_0033*IT_0040;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = (-18)*IT_0077;
    const ccomplex_t IT_0085 = 36*IT_0079;
    const ccomplex_t IT_0086 = IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0088 = (-2)*IT_0087;
    const ccomplex_t IT_0089 = cpow(s_12 + IT_0000 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0088*IT_0091;
    const ccomplex_t IT_0093 = 0.25*IT_0092;
    const ccomplex_t IT_0094 = cpow(IT_0093, 2);
    const ccomplex_t IT_0095 = (-36)*IT_0077;
    const ccomplex_t IT_0096 = IT_0000*IT_0030;
    const ccomplex_t IT_0097 = (-36)*IT_0096;
    const ccomplex_t IT_0098 = 72*IT_0079;
    const ccomplex_t IT_0099 = IT_0095 + IT_0097 + IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0101 = (-0.333333333333333)*IT_0100;
    const ccomplex_t IT_0102 = (-1.33333333333333)*IT_0100;
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = IT_0089*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (-0.0833333333333333)*IT_0092 + 0.5*IT_0105;
    const ccomplex_t IT_0107 = (-6)*IT_0076;
    const ccomplex_t IT_0108 = (-18)*IT_0057 + IT_0107;
    const ccomplex_t IT_0109 = (-18)*conj(IT_0057);
    const ccomplex_t IT_0110 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0111 = IT_0109 + IT_0110;
    const ccomplex_t IT_0112 = m_b*s_14;
    const ccomplex_t IT_0113 = (-18)*IT_0052;
    const ccomplex_t IT_0114 = m_b*IT_0076;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = (-6)*IT_0115;
    const ccomplex_t IT_0117 = IT_0113 + IT_0116;
    const ccomplex_t IT_0118 = (-18)*conj(IT_0052);
    const ccomplex_t IT_0119 = (-6)*conj(IT_0115);
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = m_b*s_13;
    const ccomplex_t IT_0122 = 6*conj(IT_0115);
    const ccomplex_t IT_0123 = 18*conj(IT_0052);
    const ccomplex_t IT_0124 = 6*conj(IT_0052);
    const ccomplex_t IT_0125 = 18*conj(IT_0115);
    const ccomplex_t IT_0126 = 18*conj(IT_0051)*IT_0083 + 18*IT_0051*conj
      (IT_0083) + 6*conj(IT_0076)*(IT_0052 + 3*IT_0115) + conj(IT_0057)*(18
      *IT_0052 + 6*IT_0115) + IT_0057*(IT_0122 + IT_0123) + IT_0076*(IT_0124 +
       IT_0125);
    const ccomplex_t IT_0127 = -s_14;
    const ccomplex_t IT_0128 = s_13 + IT_0127;
    const ccomplex_t IT_0129 = m_b*IT_0128;
    const ccomplex_t IT_0130 = 6*IT_0029;
    const ccomplex_t IT_0131 = 6*conj(IT_0029);
    const ccomplex_t IT_0132 = conj(IT_0057) + conj(IT_0083);
    const ccomplex_t IT_0133 = 18*IT_0029;
    const ccomplex_t IT_0134 = IT_0057 + IT_0083;
    const ccomplex_t IT_0135 = 18*conj(IT_0029);
    const ccomplex_t IT_0136 = IT_0106*((-36)*IT_0029 + 36*conj(IT_0029) + (
      -18)*IT_0051 + 18*conj(IT_0051) + IT_0113 + IT_0116 + IT_0122 + IT_0123) +
       (-12)*IT_0093*(IT_0029 + -conj(IT_0029) + 1./2*IT_0051 + -1./2*conj
      (IT_0051) + 1./2*IT_0052 + 3./2*IT_0115 + (-0.0833333333333333)*IT_0124 + 
      (-0.0833333333333333)*IT_0125) + conj(IT_0076)*IT_0130 + IT_0076*IT_0131 +
       IT_0132*IT_0133 + IT_0134*IT_0135;
    const ccomplex_t IT_0137 = (-18)*IT_0096;
    const ccomplex_t IT_0138 = IT_0084 + IT_0085 + IT_0137;
    const ccomplex_t IT_0139 = -IT_0106;
    const ccomplex_t IT_0140 = -IT_0093;
    const ccomplex_t IT_0141 = (-6)*IT_0096;
    const ccomplex_t IT_0142 = IT_0078 + IT_0080 + IT_0141;
    const ccomplex_t IT_0143 = -IT_0030;
    const ccomplex_t IT_0144 = IT_0029*conj(IT_0029)*(36*s_34 + (-36)*IT_0030)
       + IT_0051*conj(IT_0051)*IT_0054 + (conj(IT_0029)*(IT_0051 + IT_0052) +
       IT_0029*(conj(IT_0051) + conj(IT_0052)))*IT_0055 + (conj(IT_0057)*IT_0076
       + IT_0057*conj(IT_0076))*IT_0081 + (IT_0057*conj(IT_0057) + IT_0076*conj
      (IT_0076) + IT_0083*conj(IT_0083))*IT_0086 + -IT_0094*IT_0099 + (-2)
      *IT_0106*(IT_0093*((-12)*IT_0077 + 24*IT_0079 + (-12)*IT_0096) + 0.5
      *IT_0099*IT_0106) + IT_0096*(conj(IT_0083)*IT_0108 + IT_0083*IT_0111) + 6
      *IT_0052*(0.166666666666667*conj(IT_0052)*IT_0054 + s_34*conj(IT_0115)) +
       IT_0030*(conj(IT_0051)*IT_0117 + IT_0051*IT_0120) + IT_0112*(conj(IT_0051
      )*IT_0108 + IT_0051*IT_0111 + conj(IT_0083)*IT_0117 + IT_0083*IT_0120) +
       IT_0121*IT_0126 + IT_0129*IT_0136 + IT_0138*(conj(IT_0076)*IT_0093 +
       IT_0106*IT_0132 + IT_0134*IT_0139 + IT_0076*IT_0140) + (conj(IT_0076)
      *IT_0106 + IT_0093*IT_0132 + IT_0076*IT_0139 + IT_0134*IT_0140)*IT_0142 +
       conj(IT_0115)*IT_0130*(s_34 + IT_0143) + 6*IT_0115*(s_34*conj(IT_0052) +
       0.166666666666667*IT_0054*conj(IT_0115) + 0.166666666666667*IT_0131*(s_34
       + IT_0143));
    return create_ccomplex_return(IT_0144);
}

