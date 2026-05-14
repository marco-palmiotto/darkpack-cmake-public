#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_anti_d_d.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_anti_d_d(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = m_t*mu_h*IT_0000;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0001*IT_0003;
    const ccomplex_t IT_0013 = -IT_0008*IT_0010 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_t, 2);
    const ccomplex_t IT_0017 = (-2)*U_st_00*((0.5*m_W*IT_0013*(IT_0005 + -1./3
      *IT_0004*IT_0015) + IT_0000*IT_0005*IT_0008*IT_0011*IT_0016)*conj(U_st_00)
       + -1./2*IT_0003*IT_0005*IT_0009*IT_0011*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0005*IT_0009*IT_0011*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0004
      *IT_0013*IT_0015 + 1.5*IT_0000*IT_0005*IT_0008*IT_0011*IT_0016)*conj
      (U_st_10));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_t*IT_0008*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = pow(m_st_1, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + (-2)*IT_0021 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_d*e_em*IT_0000*IT_0002
      *IT_0005*IT_0008;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0001*IT_0008;
    const ccomplex_t IT_0028 = IT_0003*IT_0010;
    const ccomplex_t IT_0029 = IT_0027 + IT_0028;
    const ccomplex_t IT_0030 = 2*U_st_00*((IT_0000*IT_0003*IT_0005*IT_0011
      *IT_0016 + (-0.5)*m_W*(IT_0005 + -1./3*IT_0004*IT_0015)*IT_0029)*conj
      (U_st_00) + 1./2*IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_10)) + U_st_10*
      (IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_00) + (-1.33333333333333)*((
      -1.5)*IT_0000*IT_0003*IT_0005*IT_0011*IT_0016 + m_W*IT_0004*IT_0015
      *IT_0029)*conj(U_st_10));
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_t*IT_0003*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0024 + -IT_0036;
    const ccomplex_t IT_0038 = cpow(IT_0014, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0005
      *IT_0014 + -1./3*IT_0004*IT_0038)*U_st_00*conj(U_st_00) + IT_0004*IT_0038
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0040 = 0.666666666666667*IT_0039;
    const ccomplex_t IT_0041 = m_d*IT_0040;
    const ccomplex_t IT_0042 = 1.33333333333333*IT_0039;
    const ccomplex_t IT_0043 = m_d*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0038;
    const ccomplex_t IT_0045 = 0.333333333333333*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0004*IT_0038;
    const ccomplex_t IT_0050 = e_em*IT_0049;
    const ccomplex_t IT_0051 = IT_0005*IT_0014;
    const ccomplex_t IT_0052 = e_em*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0050 + 3*IT_0052);
    const ccomplex_t IT_0054 = (-0.166666666666667)*IT_0053;
    const ccomplex_t IT_0055 = IT_0046*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = (IT_0041 + -IT_0043)*(IT_0048 + -IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -(IT_0041 + -IT_0043)*(IT_0048 + -IT_0056);
    const ccomplex_t IT_0061 = 18*s_34;
    const ccomplex_t IT_0062 = pow(m_d, 2);
    const ccomplex_t IT_0063 = (-18)*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = IT_0042*IT_0056;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0005*IT_0011;
    const ccomplex_t IT_0068 = m_t*V_td*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = V_td*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0074 = IT_0005*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0072 + IT_0077;
    const ccomplex_t IT_0079 = conj(V_td)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0080 = IT_0005*IT_0079;
    const ccomplex_t IT_0081 = m_t*conj(V_td)*conj(V_u1)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0082 = IT_0067*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + (-0.5)*IT_0083);
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0021 + IT_0062 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0078*IT_0088;
    const ccomplex_t IT_0090 = m_t*V_td*V_u2*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0091 = IT_0067*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = V_td*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0096 = IT_0005*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0094 + IT_0099;
    const ccomplex_t IT_0101 = conj(V_td)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0102 = IT_0005*IT_0101;
    const ccomplex_t IT_0103 = m_t*conj(V_td)*conj(V_u2)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0104 = IT_0067*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0102 + (-0.5)*IT_0105);
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0021 + IT_0062 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0100*IT_0110;
    const ccomplex_t IT_0112 = IT_0089 + IT_0111;
    const ccomplex_t IT_0113 = IT_0042*IT_0048;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*conj
      (V_td)*e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = IT_0108*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_td*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1.4142135623731)*m_d*U_d1*conj
      (V_td)*e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0123*IT_0125;
    const ccomplex_t IT_0127 = IT_0086*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0121 + IT_0128;
    const ccomplex_t IT_0130 = s_34*IT_0021;
    const ccomplex_t IT_0131 = (-6)*IT_0130;
    const ccomplex_t IT_0132 = s_13*s_14;
    const ccomplex_t IT_0133 = 12*IT_0132;
    const ccomplex_t IT_0134 = IT_0131 + IT_0133;
    const ccomplex_t IT_0135 = (-18)*IT_0130;
    const ccomplex_t IT_0136 = 36*IT_0132;
    const ccomplex_t IT_0137 = IT_0135 + IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0139 = (-0.333333333333333)*IT_0138;
    const ccomplex_t IT_0140 = (-1.33333333333333)*IT_0138;
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = cpow(s_12 + IT_0021 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0143 = IT_0141*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0146 = (-2)*IT_0145;
    const ccomplex_t IT_0147 = IT_0142*IT_0145;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = IT_0146*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0144 + (-0.0833333333333333)*IT_0149;
    const ccomplex_t IT_0151 = cpow(IT_0150, 2);
    const ccomplex_t IT_0152 = (-36)*IT_0130;
    const ccomplex_t IT_0153 = IT_0021*IT_0062;
    const ccomplex_t IT_0154 = (-36)*IT_0153;
    const ccomplex_t IT_0155 = 72*IT_0132;
    const ccomplex_t IT_0156 = IT_0152 + IT_0154 + IT_0155;
    const ccomplex_t IT_0157 = 0.25*IT_0149;
    const ccomplex_t IT_0158 = IT_0085*IT_0123;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0160 = IT_0086*IT_0158*IT_0159;
    const ccomplex_t IT_0161 = IT_0107*IT_0116;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0163 = IT_0108*IT_0161*IT_0162;
    const ccomplex_t IT_0164 = m_d*IT_0128;
    const ccomplex_t IT_0165 = m_d*IT_0121;
    const ccomplex_t IT_0166 = -IT_0160 + -IT_0163 + -IT_0164 + -IT_0165;
    const ccomplex_t IT_0167 = IT_0078*IT_0125;
    const ccomplex_t IT_0168 = IT_0086*IT_0159*IT_0167;
    const ccomplex_t IT_0169 = IT_0100*IT_0118;
    const ccomplex_t IT_0170 = IT_0108*IT_0162*IT_0169;
    const ccomplex_t IT_0171 = m_d*IT_0072;
    const ccomplex_t IT_0172 = m_d*IT_0076;
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = IT_0171 + IT_0173;
    const ccomplex_t IT_0175 = IT_0088*IT_0174;
    const ccomplex_t IT_0176 = m_d*IT_0094;
    const ccomplex_t IT_0177 = m_d*IT_0098;
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = IT_0110*IT_0179;
    const ccomplex_t IT_0181 = -IT_0168 + -IT_0170 + -IT_0175 + -IT_0180;
    const ccomplex_t IT_0182 = (-18)*IT_0181;
    const ccomplex_t IT_0183 = (-6)*IT_0059;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = (-18)*conj(IT_0181);
    const ccomplex_t IT_0186 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = (-18)*IT_0059;
    const ccomplex_t IT_0189 = (-6)*IT_0181;
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = (-18)*conj(IT_0059);
    const ccomplex_t IT_0192 = (-6)*conj(IT_0181);
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = 6*s_34;
    const ccomplex_t IT_0195 = -IT_0062;
    const ccomplex_t IT_0196 = 6*IT_0037;
    const ccomplex_t IT_0197 = 6*conj(IT_0037);
    const ccomplex_t IT_0198 = m_d*s_13;
    const ccomplex_t IT_0199 = 6*IT_0166;
    const ccomplex_t IT_0200 = 6*conj(IT_0166);
    const ccomplex_t IT_0201 = 6*conj(IT_0059);
    const ccomplex_t IT_0202 = 18*conj(IT_0181);
    const ccomplex_t IT_0203 = 6*conj(IT_0181);
    const ccomplex_t IT_0204 = 18*conj(IT_0059);
    const ccomplex_t IT_0205 = 18*conj(IT_0166);
    const ccomplex_t IT_0206 = 3*conj(IT_0129)*(0.333333333333333*IT_0060 +
       IT_0166) + 3*conj(IT_0112)*(0.333333333333333*IT_0059 + IT_0181) + 3*conj
      (IT_0066)*(IT_0059 + 1./3*IT_0181) + 3*conj(IT_0114)*(IT_0060 +
       0.0555555555555556*IT_0199) + 3*IT_0114*(conj(IT_0060) +
       0.0555555555555556*IT_0200) + 0.166666666666667*IT_0112*(IT_0201 +
       IT_0202) + 0.166666666666667*IT_0066*(IT_0203 + IT_0204) + IT_0129*(conj
      (IT_0060) + 0.166666666666667*IT_0205);
    const ccomplex_t IT_0207 = m_d*e_em*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = m_t*e_em*mu_h*IT_0000*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0210 = IT_0005*IT_0209;
    const ccomplex_t IT_0211 = m_t*e_em*mu_h*IT_0000*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0212 = IT_0005*IT_0211;
    const ccomplex_t IT_0213 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0214 = IT_0001*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0217 = IT_0001*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (-0.5)*IT_0210 + 0.5*IT_0212 + (-0.5)*IT_0215 +
       0.5*IT_0218;
    const ccomplex_t IT_0220 = IT_0208*IT_0219;
    const ccomplex_t IT_0221 = cpow((-2)*s_12 + (-2)*IT_0021 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0222 = IT_0220*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = m_d*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0225 = 0.5*IT_0224;
    const ccomplex_t IT_0226 = IT_0001*IT_0005*IT_0011;
    const ccomplex_t IT_0227 = IT_0209*IT_0226;
    const ccomplex_t IT_0228 = IT_0211*IT_0226;
    const ccomplex_t IT_0229 = IT_0010*IT_0213;
    const ccomplex_t IT_0230 = 1.4142135623731*IT_0229;
    const ccomplex_t IT_0231 = IT_0010*IT_0216;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = 0.5*IT_0227 + (-0.5)*IT_0228 + (-0.5)*IT_0230 +
       0.5*IT_0232;
    const ccomplex_t IT_0234 = IT_0225*IT_0233;
    const ccomplex_t IT_0235 = pow(m_Z, 2);
    const ccomplex_t IT_0236 = cpow((-2)*s_12 + (-2)*IT_0021 + IT_0235 + 
      -reg_prop, -1);
    const ccomplex_t IT_0237 = IT_0234*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = -IT_0223 + -IT_0238;
    const ccomplex_t IT_0240 = s_13 + s_14;
    const ccomplex_t IT_0241 = m_d*IT_0240;
    const ccomplex_t IT_0242 = 18*IT_0062;
    const ccomplex_t IT_0243 = IT_0061 + IT_0242;
    const ccomplex_t IT_0244 = (-18)*s_34 + IT_0063;
    const ccomplex_t IT_0245 = s_34 + IT_0062;
    const ccomplex_t IT_0246 = -s_14;
    const ccomplex_t IT_0247 = s_13 + IT_0246;
    const ccomplex_t IT_0248 = m_d*IT_0247;
    const ccomplex_t IT_0249 = conj(IT_0112) + conj(IT_0129);
    const ccomplex_t IT_0250 = IT_0112 + IT_0129;
    const ccomplex_t IT_0251 = conj(IT_0066) + conj(IT_0114);
    const ccomplex_t IT_0252 = 18*IT_0037;
    const ccomplex_t IT_0253 = IT_0066 + IT_0114;
    const ccomplex_t IT_0254 = 18*conj(IT_0037);
    const ccomplex_t IT_0255 = (-18)*IT_0166;
    const ccomplex_t IT_0256 = (-6)*IT_0166;
    const ccomplex_t IT_0257 = IT_0196*IT_0249 + IT_0197*IT_0250 + IT_0251
      *IT_0252 + IT_0253*IT_0254 + IT_0157*((-12)*IT_0037 + 12*conj(IT_0037) + (
      -6)*IT_0060 + 6*conj(IT_0060) + IT_0182 + IT_0183 + IT_0201 + IT_0202 +
       IT_0205 + IT_0255) + IT_0150*((-36)*IT_0037 + 36*conj(IT_0037) + (-18)
      *IT_0060 + 18*conj(IT_0060) + IT_0188 + IT_0189 + IT_0200 + IT_0203 +
       IT_0204 + IT_0256);
    const ccomplex_t IT_0258 = (-18)*IT_0153;
    const ccomplex_t IT_0259 = IT_0135 + IT_0136 + IT_0258;
    const ccomplex_t IT_0260 = -IT_0157;
    const ccomplex_t IT_0261 = -IT_0150;
    const ccomplex_t IT_0262 = (-6)*IT_0153;
    const ccomplex_t IT_0263 = IT_0131 + IT_0133 + IT_0262;
    const ccomplex_t IT_0264 = 36*s_34;
    const ccomplex_t IT_0265 = IT_0239*(((-18)*conj(IT_0066) + (-6)*conj
      (IT_0112) + 18*conj(IT_0114) + 6*conj(IT_0129))*IT_0241 + conj(IT_0060)
      *IT_0243 + conj(IT_0059)*IT_0244 + (IT_0192 + IT_0200)*IT_0245 + conj
      (IT_0239)*(36*IT_0062 + IT_0264));
    const ccomplex_t IT_0266 = 3*IT_0112;
    const ccomplex_t IT_0267 = IT_0066 + IT_0266;
    const ccomplex_t IT_0268 = (-18)*conj(IT_0066);
    const ccomplex_t IT_0269 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = (-18)*IT_0066;
    const ccomplex_t IT_0272 = (-6)*IT_0112;
    const ccomplex_t IT_0273 = IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = (-0.166666666666667)*conj(IT_0114);
    const ccomplex_t IT_0275 = m_d*s_14;
    const ccomplex_t IT_0276 = (-0.166666666666667)*conj(IT_0129)*IT_0184 + (
      -0.166666666666667)*IT_0129*IT_0187 + (-0.166666666666667)*IT_0114*IT_0193
       + (-0.166666666666667)*conj(IT_0112)*IT_0255 + (-0.166666666666667)*conj
      (IT_0066)*IT_0256 + conj(IT_0166)*IT_0267 + (-0.166666666666667)*IT_0060
      *IT_0270 + (-0.166666666666667)*conj(IT_0060)*IT_0273 + IT_0190*IT_0274;
    const ccomplex_t IT_0277 = (conj(IT_0037)*(IT_0059 + IT_0060) + IT_0037*
      (conj(IT_0059) + conj(IT_0060)))*IT_0064 + (conj(IT_0066)*IT_0112 +
       IT_0066*conj(IT_0112) + conj(IT_0114)*IT_0129 + IT_0114*conj(IT_0129))
      *IT_0134 + (IT_0066*conj(IT_0066) + IT_0112*conj(IT_0112) + IT_0114*conj
      (IT_0114) + IT_0129*conj(IT_0129))*IT_0137 + -IT_0151*IT_0156 + (-2)
      *IT_0157*(IT_0150*((-12)*IT_0130 + 24*IT_0132 + (-12)*IT_0153) + 0.5
      *IT_0156*IT_0157) + IT_0061*(IT_0059*conj(IT_0059) + IT_0060*conj(IT_0060)
       + IT_0166*conj(IT_0166) + IT_0181*conj(IT_0181)) + IT_0062*(conj(IT_0166)
      *IT_0184 + IT_0166*IT_0187 + conj(IT_0060)*IT_0190 + IT_0060*IT_0193) + 
      (conj(IT_0060)*IT_0166 + IT_0060*conj(IT_0166) + conj(IT_0059)*IT_0181 +
       IT_0059*conj(IT_0181))*IT_0194 + (s_34 + IT_0195)*((conj(IT_0166) + conj
      (IT_0181))*IT_0196 + (IT_0166 + IT_0181)*IT_0197) + 6*IT_0198*IT_0206 +
       conj(IT_0239)*(((-18)*IT_0066 + (-6)*IT_0112 + 18*IT_0114 + 6*IT_0129)
      *IT_0241 + IT_0060*IT_0243 + IT_0059*IT_0244 + (IT_0189 + IT_0199)*IT_0245
      ) + IT_0248*IT_0257 + IT_0259*(IT_0157*IT_0249 + IT_0150*IT_0251 + IT_0250
      *IT_0260 + IT_0253*IT_0261) + (IT_0150*IT_0249 + IT_0157*IT_0251 + IT_0253
      *IT_0260 + IT_0250*IT_0261)*IT_0263 + IT_0037*conj(IT_0037)*((-36)*IT_0062
       + IT_0264) + IT_0265 + (-6)*IT_0153*((conj(IT_0066) + 3*conj(IT_0112))
      *IT_0129 + conj(IT_0129)*IT_0267 + (-0.166666666666667)*IT_0114*IT_0270 +
       IT_0273*IT_0274) + (-6)*IT_0275*IT_0276;
    return create_ccomplex_return(IT_0277);
}

