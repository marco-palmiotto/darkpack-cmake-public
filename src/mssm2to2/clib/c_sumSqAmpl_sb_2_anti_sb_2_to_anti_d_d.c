#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = m_d*e_em*IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = m_b*e_em*mu_h*IT_0000*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = m_b*e_em*mu_h*IT_0000*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = cpow(V_tb, 2);
    const ccomplex_t IT_0013 = A_b*IT_0012*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = A_b*IT_0012*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0009 + 0.5*IT_0011 + (-0.5)*IT_0015 +
       0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = pow(m_sb_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + (-2)*IT_0021 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = pow(m_Z, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + (-2)*IT_0021 + IT_0025 + 
      -reg_prop, -1);
    const ccomplex_t IT_0027 = m_d*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0030 = IT_0008*IT_0029;
    const ccomplex_t IT_0031 = IT_0010*IT_0029;
    const ccomplex_t IT_0032 = IT_0001*IT_0013;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = IT_0001*IT_0016;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0030 + 0.5*IT_0031 + 0.5*IT_0033 + (
      -0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = IT_0026*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = -IT_0024 + -IT_0039;
    const ccomplex_t IT_0041 = pow(m_d, 2);
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*((-1.5)*(IT_0005
      *IT_0042 + 1./3*IT_0004*IT_0043)*U_sb_01*conj(U_sb_01) + IT_0004*IT_0043
      *U_sb_11*conj(U_sb_11));
    const ccomplex_t IT_0045 = (-0.666666666666667)*IT_0044;
    const ccomplex_t IT_0046 = IT_0004*IT_0043;
    const ccomplex_t IT_0047 = e_em*IT_0046;
    const ccomplex_t IT_0048 = IT_0005*IT_0042;
    const ccomplex_t IT_0049 = e_em*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0047 + 3*IT_0049);
    const ccomplex_t IT_0051 = (-0.166666666666667)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0045*IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = s_13 + s_14;
    const ccomplex_t IT_0058 = m_d*IT_0057;
    const ccomplex_t IT_0059 = (-0.333333333333333)*IT_0044;
    const ccomplex_t IT_0060 = m_d*IT_0059;
    const ccomplex_t IT_0061 = m_d*IT_0045;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0043;
    const ccomplex_t IT_0063 = 0.333333333333333*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = (IT_0060 + -IT_0061)*(IT_0054 + -IT_0065);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = 18*s_34;
    const ccomplex_t IT_0070 = 18*IT_0041 + IT_0069;
    const ccomplex_t IT_0071 = -IT_0068;
    const ccomplex_t IT_0072 = (-18)*IT_0041;
    const ccomplex_t IT_0073 = (-18)*s_34 + IT_0072;
    const ccomplex_t IT_0074 = cpow(s_12 + IT_0021 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = (-0.333333333333333)*IT_0075;
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = IT_0074*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0082 = (-2)*IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0081;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0080 + (-0.0833333333333333)*IT_0085;
    const ccomplex_t IT_0087 = 0.25*IT_0085;
    const ccomplex_t IT_0088 = s_13*s_14;
    const ccomplex_t IT_0089 = IT_0021*IT_0041;
    const ccomplex_t IT_0090 = s_34*IT_0021;
    const ccomplex_t IT_0091 = (-18)*IT_0090;
    const ccomplex_t IT_0092 = (-18)*IT_0089;
    const ccomplex_t IT_0093 = 36*IT_0088;
    const ccomplex_t IT_0094 = IT_0091 + IT_0092 + IT_0093;
    const ccomplex_t IT_0095 = (-36)*IT_0090;
    const ccomplex_t IT_0096 = (-36)*IT_0089;
    const ccomplex_t IT_0097 = 72*IT_0088;
    const ccomplex_t IT_0098 = IT_0095 + IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = (-6)*IT_0090;
    const ccomplex_t IT_0100 = (-6)*IT_0089;
    const ccomplex_t IT_0101 = 12*IT_0088;
    const ccomplex_t IT_0102 = IT_0099 + IT_0100 + IT_0101;
    const ccomplex_t IT_0103 = cos(alpha);
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = m_b*mu_h*IT_0000;
    const ccomplex_t IT_0107 = sin(alpha);
    const ccomplex_t IT_0108 = -IT_0003*IT_0107;
    const ccomplex_t IT_0109 = IT_0001*IT_0103;
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0042, -2);
    const ccomplex_t IT_0112 = pow(m_b, 2);
    const ccomplex_t IT_0113 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0110*(IT_0005 + 1.
      /3*IT_0004*IT_0111) + IT_0000*IT_0002*IT_0005*IT_0103*IT_0112)*conj
      (U_sb_01) + -1./2*IT_0002*IT_0005*IT_0106*IT_0107*conj(U_sb_11)) + U_sb_11
      *(IT_0002*IT_0005*IT_0106*IT_0107*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0004*IT_0110*IT_0111 + (-3)*IT_0000*IT_0002*IT_0005*IT_0103*IT_0112)
      *conj(U_sb_11));
    const ccomplex_t IT_0114 = (0 + _Complex_I*-1)*(e_em*IT_0113 + (
      -1.4142135623731)*A_b*IT_0012*IT_0103*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = IT_0105*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_12 + (-2)*IT_0021 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0107;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0001*IT_0107;
    const ccomplex_t IT_0124 = IT_0003*IT_0103;
    const ccomplex_t IT_0125 = IT_0123 + IT_0124;
    const ccomplex_t IT_0126 = 2*U_sb_01*((IT_0000*IT_0002*IT_0005*IT_0107
      *IT_0112 + (-0.5)*m_W*(IT_0005 + 1./3*IT_0004*IT_0111)*IT_0125)*conj
      (U_sb_01) + 1./2*IT_0002*IT_0005*IT_0103*IT_0106*conj(U_sb_11)) + U_sb_11*
      (IT_0002*IT_0005*IT_0103*IT_0106*conj(U_sb_01) + 2*(IT_0000*IT_0002
      *IT_0005*IT_0107*IT_0112 + -1./3*m_W*IT_0004*IT_0111*IT_0125)*conj(U_sb_11
      ));
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(e_em*IT_0126 +
       1.4142135623731*A_b*IT_0012*IT_0107*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0122*IT_0128;
    const ccomplex_t IT_0130 = IT_0120*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = -IT_0119 + -IT_0131;
    const ccomplex_t IT_0133 = IT_0069 + IT_0072;
    const ccomplex_t IT_0134 = -s_14;
    const ccomplex_t IT_0135 = s_13 + IT_0134;
    const ccomplex_t IT_0136 = m_d*IT_0135;
    const ccomplex_t IT_0137 = IT_0045*IT_0065;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = (-18)*IT_0071;
    const ccomplex_t IT_0140 = (-6)*IT_0087*(IT_0068 + -conj(IT_0068) +
       IT_0071 + -conj(IT_0071) + 2*IT_0132 + (-2)*conj(IT_0132)) + conj(IT_0132
      )*(18*IT_0056 + 18*IT_0138) + IT_0132*(18*conj(IT_0056) + 18*conj(IT_0138)
      ) + IT_0086*((-18)*IT_0068 + 18*conj(IT_0068) + 18*conj(IT_0071) + (-36)
      *IT_0132 + 36*conj(IT_0132) + IT_0139);
    const ccomplex_t IT_0141 = m_d*s_14;
    const ccomplex_t IT_0142 = (-18)*IT_0056;
    const ccomplex_t IT_0143 = m_d*s_13;
    const ccomplex_t IT_0144 = conj(IT_0132)*IT_0133;
    const ccomplex_t IT_0145 = (-18)*conj(IT_0071);
    const ccomplex_t IT_0146 = (-18)*conj(IT_0056);
    const ccomplex_t IT_0147 = IT_0091 + IT_0093;
    const ccomplex_t IT_0148 = -IT_0087;
    const ccomplex_t IT_0149 = IT_0102*IT_0148;
    const ccomplex_t IT_0150 = IT_0040*(conj(IT_0040)*(36*s_34 + 36*IT_0041) +
       (-18)*conj(IT_0056)*IT_0058 + conj(IT_0068)*IT_0070 + conj(IT_0071)
      *IT_0073) + (-2)*IT_0086*(IT_0087*(24*IT_0088 + (-12)*IT_0089 + (-12)
      *IT_0090) + (-0.5)*conj(IT_0056)*IT_0094 + 0.5*IT_0086*IT_0098) + -IT_0087
      *(IT_0087*IT_0098 + -conj(IT_0056)*IT_0102) + IT_0132*((36*s_34 + (-36)
      *IT_0041)*conj(IT_0132) + (conj(IT_0068) + conj(IT_0071))*IT_0133) +
       IT_0136*IT_0140 + conj(IT_0068)*(IT_0041*IT_0139 + IT_0141*IT_0142) +
       conj(IT_0138)*(18*IT_0040*IT_0058 + IT_0086*IT_0094 + IT_0087*IT_0102 +
       IT_0139*IT_0141 + IT_0089*IT_0142 + 18*IT_0068*IT_0143) + IT_0071*
      (IT_0069*conj(IT_0071) + conj(IT_0040)*IT_0073 + 18*conj(IT_0056)*IT_0143 
      + IT_0144) + IT_0068*(conj(IT_0068)*IT_0069 + conj(IT_0040)*IT_0070 +
       IT_0144 + IT_0041*IT_0145 + IT_0141*IT_0146) + IT_0056*((-18)*conj
      (IT_0040)*IT_0058 + -IT_0086*IT_0094 + 18*conj(IT_0071)*IT_0143 + conj
      (IT_0056)*IT_0147 + IT_0149) + IT_0138*(18*conj(IT_0040)*IT_0058 + 
      -IT_0086*IT_0094 + 18*conj(IT_0068)*IT_0143 + IT_0141*IT_0145 + IT_0089
      *IT_0146 + conj(IT_0138)*IT_0147 + IT_0149);
    return create_ccomplex_return(IT_0150);
}

