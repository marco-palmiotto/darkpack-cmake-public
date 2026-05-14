#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
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
    const creal_t m_st_2 = param->m_st_2;
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
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
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
    const ccomplex_t IT_0017 = (-2)*U_st_01*((0.5*m_W*IT_0013*(IT_0005 + -1./3
      *IT_0004*IT_0015) + IT_0000*IT_0005*IT_0008*IT_0011*IT_0016)*conj(U_st_01)
       + -1./2*IT_0003*IT_0005*IT_0009*IT_0011*conj(U_st_11)) + U_st_11*(IT_0003
      *IT_0005*IT_0009*IT_0011*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0004
      *IT_0013*IT_0015 + 1.5*IT_0000*IT_0005*IT_0008*IT_0011*IT_0016)*conj
      (U_st_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*-1)*(e_em*IT_0017 + (
      -1.4142135623731)*A_t*IT_0008*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0007*IT_0019;
    const ccomplex_t IT_0021 = pow(m_st_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + (-2)*IT_0021 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_b*e_em*IT_0000*IT_0002
      *IT_0005*IT_0008;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0001*IT_0008;
    const ccomplex_t IT_0028 = IT_0003*IT_0010;
    const ccomplex_t IT_0029 = IT_0027 + IT_0028;
    const ccomplex_t IT_0030 = 2*U_st_01*((IT_0000*IT_0003*IT_0005*IT_0011
      *IT_0016 + (-0.5)*m_W*(IT_0005 + -1./3*IT_0004*IT_0015)*IT_0029)*conj
      (U_st_01) + 1./2*IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_11)) + U_st_11*
      (IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0000*IT_0003*IT_0005*IT_0011*IT_0016 + m_W*IT_0004*IT_0015
      *IT_0029)*conj(U_st_11));
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_t*IT_0003*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
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
      *IT_0014 + -1./3*IT_0004*IT_0038)*U_st_01*conj(U_st_01) + IT_0004*IT_0038
      *U_st_11*conj(U_st_11));
    const ccomplex_t IT_0040 = 0.666666666666667*IT_0039;
    const ccomplex_t IT_0041 = m_b*IT_0040;
    const ccomplex_t IT_0042 = 1.33333333333333*IT_0039;
    const ccomplex_t IT_0043 = m_b*IT_0042;
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
    const ccomplex_t IT_0060 = 18*s_34;
    const ccomplex_t IT_0061 = pow(m_b, 2);
    const ccomplex_t IT_0062 = (-18)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060 + IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_01);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0000*IT_0002*IT_0005*U_st_01;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0021 + IT_0061 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_01);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0000*IT_0002*IT_0005*U_st_01;
    const ccomplex_t IT_0075 = 0.5*IT_0074;
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = cpow((-2)*s_23 + IT_0021 + IT_0061 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0071 + IT_0079;
    const ccomplex_t IT_0081 = IT_0042*IT_0048;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0005*IT_0011;
    const ccomplex_t IT_0084 = m_t*V_tb*V_u1*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = 1.4142135623731*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0090 = IT_0005*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = IT_0088 + IT_0093;
    const ccomplex_t IT_0095 = V_tb*e_em*conj(V_Wp1)*U_st_01;
    const ccomplex_t IT_0096 = IT_0005*IT_0095;
    const ccomplex_t IT_0097 = m_t*V_tb*conj(V_u1)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0098 = IT_0083*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + (-0.5)*IT_0099);
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0077*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0094*IT_0103;
    const ccomplex_t IT_0105 = m_t*V_tb*V_u2*e_em*IT_0000*conj(U_st_11);
    const ccomplex_t IT_0106 = IT_0083*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0111 = IT_0005*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0109 + IT_0114;
    const ccomplex_t IT_0116 = V_tb*e_em*conj(V_Wp2)*U_st_01;
    const ccomplex_t IT_0117 = IT_0005*IT_0116;
    const ccomplex_t IT_0118 = m_t*V_tb*conj(V_u2)*e_em*IT_0000*U_st_11;
    const ccomplex_t IT_0119 = IT_0083*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0117 + (-0.5)*IT_0120);
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = IT_0069*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0115*IT_0124;
    const ccomplex_t IT_0126 = IT_0104 + IT_0125;
    const ccomplex_t IT_0127 = IT_0042*IT_0056;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_34*IT_0021;
    const ccomplex_t IT_0130 = (-6)*IT_0129;
    const ccomplex_t IT_0131 = s_13*s_14;
    const ccomplex_t IT_0132 = 12*IT_0131;
    const ccomplex_t IT_0133 = IT_0130 + IT_0132;
    const ccomplex_t IT_0134 = (-18)*IT_0129;
    const ccomplex_t IT_0135 = 36*IT_0131;
    const ccomplex_t IT_0136 = IT_0134 + IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0138 = (-0.333333333333333)*IT_0137;
    const ccomplex_t IT_0139 = (-1.33333333333333)*IT_0137;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = cpow(s_12 + IT_0021 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0145 = (-2)*IT_0144;
    const ccomplex_t IT_0146 = IT_0141*IT_0144;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0145*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0143 + (-0.0833333333333333)*IT_0148;
    const ccomplex_t IT_0150 = cpow(IT_0149, 2);
    const ccomplex_t IT_0151 = (-36)*IT_0129;
    const ccomplex_t IT_0152 = IT_0021*IT_0061;
    const ccomplex_t IT_0153 = (-36)*IT_0152;
    const ccomplex_t IT_0154 = 72*IT_0131;
    const ccomplex_t IT_0155 = IT_0151 + IT_0153 + IT_0154;
    const ccomplex_t IT_0156 = 0.25*IT_0148;
    const ccomplex_t IT_0157 = IT_0073*IT_0101;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0159 = IT_0077*IT_0157*IT_0158;
    const ccomplex_t IT_0160 = IT_0065*IT_0122;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0162 = IT_0069*IT_0160*IT_0161;
    const ccomplex_t IT_0163 = m_b*IT_0071;
    const ccomplex_t IT_0164 = m_b*IT_0079;
    const ccomplex_t IT_0165 = -IT_0159 + -IT_0162 + -IT_0163 + -IT_0164;
    const ccomplex_t IT_0166 = IT_0075*IT_0094;
    const ccomplex_t IT_0167 = IT_0077*IT_0158*IT_0166;
    const ccomplex_t IT_0168 = IT_0067*IT_0115;
    const ccomplex_t IT_0169 = IT_0069*IT_0161*IT_0168;
    const ccomplex_t IT_0170 = m_b*IT_0088;
    const ccomplex_t IT_0171 = m_b*IT_0092;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = IT_0103*IT_0173;
    const ccomplex_t IT_0175 = m_b*IT_0109;
    const ccomplex_t IT_0176 = m_b*IT_0113;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0175 + IT_0177;
    const ccomplex_t IT_0179 = IT_0124*IT_0178;
    const ccomplex_t IT_0180 = -IT_0167 + -IT_0169 + -IT_0174 + -IT_0179;
    const ccomplex_t IT_0181 = (-18)*IT_0180;
    const ccomplex_t IT_0182 = (-6)*IT_0059;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = (-18)*conj(IT_0180);
    const ccomplex_t IT_0185 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0186 = IT_0184 + IT_0185;
    const ccomplex_t IT_0187 = (-18)*IT_0059;
    const ccomplex_t IT_0188 = (-6)*IT_0180;
    const ccomplex_t IT_0189 = IT_0187 + IT_0188;
    const ccomplex_t IT_0190 = (-18)*conj(IT_0059);
    const ccomplex_t IT_0191 = (-6)*conj(IT_0180);
    const ccomplex_t IT_0192 = IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = (-18)*IT_0126 + (-6)*IT_0128;
    const ccomplex_t IT_0194 = (-18)*IT_0128;
    const ccomplex_t IT_0195 = (-6)*IT_0126;
    const ccomplex_t IT_0196 = IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = (-18)*conj(IT_0128);
    const ccomplex_t IT_0198 = (-6)*conj(IT_0126);
    const ccomplex_t IT_0199 = IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = 6*s_34;
    const ccomplex_t IT_0201 = -IT_0061;
    const ccomplex_t IT_0202 = 6*IT_0037;
    const ccomplex_t IT_0203 = 6*conj(IT_0037);
    const ccomplex_t IT_0204 = m_b*e_em*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = m_t*e_em*mu_h*IT_0000*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0207 = IT_0005*IT_0206;
    const ccomplex_t IT_0208 = m_t*e_em*mu_h*IT_0000*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0209 = IT_0005*IT_0208;
    const ccomplex_t IT_0210 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0211 = IT_0001*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0214 = IT_0001*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0207 + 0.5*IT_0209 + (-0.5)*IT_0212 +
       0.5*IT_0215;
    const ccomplex_t IT_0217 = IT_0205*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_12 + (-2)*IT_0021 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0219 = IT_0217*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = m_b*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = IT_0001*IT_0005*IT_0011;
    const ccomplex_t IT_0224 = IT_0206*IT_0223;
    const ccomplex_t IT_0225 = IT_0208*IT_0223;
    const ccomplex_t IT_0226 = IT_0010*IT_0210;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = IT_0010*IT_0213;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = 0.5*IT_0224 + (-0.5)*IT_0225 + (-0.5)*IT_0227 +
       0.5*IT_0229;
    const ccomplex_t IT_0231 = IT_0222*IT_0230;
    const ccomplex_t IT_0232 = pow(m_Z, 2);
    const ccomplex_t IT_0233 = cpow((-2)*s_12 + (-2)*IT_0021 + IT_0232 + 
      -reg_prop, -1);
    const ccomplex_t IT_0234 = IT_0231*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = -IT_0220 + -IT_0235;
    const ccomplex_t IT_0237 = s_13 + s_14;
    const ccomplex_t IT_0238 = m_b*IT_0237;
    const ccomplex_t IT_0239 = 18*IT_0061;
    const ccomplex_t IT_0240 = IT_0060 + IT_0239;
    const ccomplex_t IT_0241 = (-18)*s_34 + IT_0062;
    const ccomplex_t IT_0242 = s_34 + IT_0061;
    const ccomplex_t IT_0243 = 6*IT_0165;
    const ccomplex_t IT_0244 = m_b*s_14;
    const ccomplex_t IT_0245 = (-18)*IT_0165;
    const ccomplex_t IT_0246 = (-6)*IT_0165;
    const ccomplex_t IT_0247 = (-6)*IT_0152;
    const ccomplex_t IT_0248 = IT_0130 + IT_0132 + IT_0247;
    const ccomplex_t IT_0249 = conj(IT_0080) + conj(IT_0126);
    const ccomplex_t IT_0250 = conj(IT_0082) + conj(IT_0128);
    const ccomplex_t IT_0251 = IT_0082 + IT_0128;
    const ccomplex_t IT_0252 = -IT_0156;
    const ccomplex_t IT_0253 = IT_0080 + IT_0126;
    const ccomplex_t IT_0254 = -IT_0149;
    const ccomplex_t IT_0255 = (-18)*IT_0152;
    const ccomplex_t IT_0256 = IT_0134 + IT_0135 + IT_0255;
    const ccomplex_t IT_0257 = 36*s_34;
    const ccomplex_t IT_0258 = 6*conj(IT_0165);
    const ccomplex_t IT_0259 = (conj(IT_0080) + 3*conj(IT_0082) + -conj
      (IT_0126) + (-3)*conj(IT_0128))*IT_0238 + 0.166666666666667*conj(IT_0058)
      *IT_0240 + 0.166666666666667*conj(IT_0059)*IT_0241 + 0.166666666666667
      *conj(IT_0236)*(36*IT_0061 + IT_0257) + 0.166666666666667*IT_0242*(IT_0191
       + IT_0258);
    const ccomplex_t IT_0260 = m_b*s_13;
    const ccomplex_t IT_0261 = 6*conj(IT_0059);
    const ccomplex_t IT_0262 = 18*conj(IT_0180);
    const ccomplex_t IT_0263 = 6*conj(IT_0180);
    const ccomplex_t IT_0264 = 18*conj(IT_0059);
    const ccomplex_t IT_0265 = 18*conj(IT_0165);
    const ccomplex_t IT_0266 = 18*conj(IT_0058);
    const ccomplex_t IT_0267 = IT_0059*(conj(IT_0126) + 3*conj(IT_0128)) + 1
      *conj(IT_0080)*(IT_0058 + 3*IT_0165) + 0.166666666666667*(18*conj(IT_0126)
       + 6*conj(IT_0128))*IT_0180 + 0.166666666666667*conj(IT_0082)*(18*IT_0058 
      + IT_0243) + 0.166666666666667*IT_0126*(IT_0261 + IT_0262) +
       0.166666666666667*IT_0128*(IT_0263 + IT_0264) + IT_0080*(conj(IT_0058) +
       0.166666666666667*IT_0265) + 0.166666666666667*IT_0082*(IT_0258 + IT_0266);
    const ccomplex_t IT_0268 = -s_14;
    const ccomplex_t IT_0269 = s_13 + IT_0268;
    const ccomplex_t IT_0270 = m_b*IT_0269;
    const ccomplex_t IT_0271 = 18*IT_0037;
    const ccomplex_t IT_0272 = 18*conj(IT_0037);
    const ccomplex_t IT_0273 = IT_0202*IT_0249 + IT_0203*IT_0253 + IT_0156*((
      -12)*IT_0037 + 12*conj(IT_0037) + (-6)*IT_0058 + 6*conj(IT_0058) + IT_0181
       + IT_0182 + IT_0245 + IT_0261 + IT_0262 + IT_0265) + IT_0149*((-36)
      *IT_0037 + 36*conj(IT_0037) + (-18)*IT_0058 + IT_0187 + IT_0188 + IT_0246 
      + IT_0258 + IT_0263 + IT_0264 + IT_0266) + IT_0250*IT_0271 + IT_0251
      *IT_0272;
    const ccomplex_t IT_0274 = (conj(IT_0037)*(IT_0058 + IT_0059) + IT_0037*
      (conj(IT_0058) + conj(IT_0059)))*IT_0063 + (conj(IT_0080)*IT_0082 +
       IT_0080*conj(IT_0082) + conj(IT_0126)*IT_0128 + IT_0126*conj(IT_0128))
      *IT_0133 + (IT_0080*conj(IT_0080) + IT_0082*conj(IT_0082) + IT_0126*conj
      (IT_0126) + IT_0128*conj(IT_0128))*IT_0136 + -IT_0150*IT_0155 + (-2)
      *IT_0156*(IT_0149*((-12)*IT_0129 + 24*IT_0131 + (-12)*IT_0152) + 0.5
      *IT_0155*IT_0156) + IT_0060*(IT_0058*conj(IT_0058) + IT_0059*conj(IT_0059)
       + IT_0165*conj(IT_0165) + IT_0180*conj(IT_0180)) + IT_0061*(conj(IT_0165)
      *IT_0183 + IT_0165*IT_0186 + conj(IT_0058)*IT_0189 + IT_0058*IT_0192) +
       IT_0152*(IT_0080*((-18)*conj(IT_0126) + (-6)*conj(IT_0128)) + conj
      (IT_0080)*IT_0193 + conj(IT_0082)*IT_0196 + IT_0082*IT_0199) + (conj
      (IT_0058)*IT_0165 + IT_0058*conj(IT_0165) + conj(IT_0059)*IT_0180 +
       IT_0059*conj(IT_0180))*IT_0200 + (s_34 + IT_0201)*((conj(IT_0165) + conj
      (IT_0180))*IT_0202 + (IT_0165 + IT_0180)*IT_0203) + 6*conj(IT_0236)*(
      (IT_0080 + 3*IT_0082 + -IT_0126 + (-3)*IT_0128)*IT_0238 +
       0.166666666666667*IT_0058*IT_0240 + 0.166666666666667*IT_0059*IT_0241 +
       0.166666666666667*IT_0242*(IT_0188 + IT_0243)) + IT_0244*(conj(IT_0080)
      *IT_0183 + IT_0080*IT_0186 + conj(IT_0082)*IT_0189 + IT_0082*IT_0192 +
       conj(IT_0165)*IT_0193 + conj(IT_0058)*IT_0196 + IT_0058*IT_0199 + conj
      (IT_0126)*IT_0245 + conj(IT_0128)*IT_0246) + IT_0248*(IT_0149*IT_0249 +
       IT_0156*IT_0250 + IT_0251*IT_0252 + IT_0253*IT_0254) + (IT_0156*IT_0249 +
       IT_0149*IT_0250 + IT_0252*IT_0253 + IT_0251*IT_0254)*IT_0256 + IT_0037
      *conj(IT_0037)*((-36)*IT_0061 + IT_0257) + 6*IT_0236*IT_0259 + 6*IT_0260
      *IT_0267 + IT_0270*IT_0273;
    return create_ccomplex_return(IT_0274);
}

