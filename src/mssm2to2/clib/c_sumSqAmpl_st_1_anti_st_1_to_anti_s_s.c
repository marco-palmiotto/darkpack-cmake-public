#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_anti_s_s.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_anti_s_s(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
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
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = m_t*mu_h*IT_0000;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0001*IT_0003;
    const ccomplex_t IT_0013 = IT_0008*IT_0010;
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_t, 2);
    const ccomplex_t IT_0018 = 2*U_st_00*(((-0.5)*m_W*IT_0014*(IT_0005 + -1./3
      *IT_0004*IT_0016) + IT_0000*IT_0005*IT_0008*IT_0011*IT_0017)*conj(U_st_00)
       + 1./2*IT_0003*IT_0005*IT_0009*IT_0011*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0005*IT_0009*IT_0011*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0004
      *IT_0014*IT_0016 + (-1.5)*IT_0000*IT_0005*IT_0008*IT_0011*IT_0017)*conj
      (U_st_10));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(e_em*IT_0018 +
       1.4142135623731*A_t*IT_0008*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = pow(m_st_1, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + (-2)*IT_0022 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0005*IT_0008;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0001*IT_0008;
    const ccomplex_t IT_0030 = -IT_0003*IT_0010 + IT_0029;
    const ccomplex_t IT_0031 = (-2)*U_st_00*((IT_0000*IT_0003*IT_0005*IT_0011
      *IT_0017 + 0.5*m_W*(IT_0005 + -1./3*IT_0004*IT_0016)*IT_0030)*conj(U_st_00
      ) + -1./2*IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_10)) + U_st_10*
      (IT_0005*IT_0008*IT_0009*IT_0011*conj(U_st_00) + (-1.33333333333333)*(1.5
      *IT_0000*IT_0003*IT_0005*IT_0011*IT_0017 + m_W*IT_0004*IT_0016*IT_0030)
      *conj(U_st_10));
    const ccomplex_t IT_0032 = (0 + _Complex_I*-1)*(e_em*IT_0031 + (
      -1.4142135623731)*A_t*IT_0003*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = IT_0026*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0025 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_s, 2);
    const ccomplex_t IT_0039 = cpow(IT_0015, -1);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0005
      *IT_0015 + -1./3*IT_0004*IT_0039)*U_st_00*conj(U_st_00) + IT_0004*IT_0039
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0041 = 0.666666666666667*IT_0040;
    const ccomplex_t IT_0042 = m_s*IT_0041;
    const ccomplex_t IT_0043 = 1.33333333333333*IT_0040;
    const ccomplex_t IT_0044 = m_s*IT_0043;
    const ccomplex_t IT_0045 = IT_0004*IT_0039;
    const ccomplex_t IT_0046 = e_em*IT_0045;
    const ccomplex_t IT_0047 = IT_0005*IT_0015;
    const ccomplex_t IT_0048 = e_em*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0046 + 3*IT_0048);
    const ccomplex_t IT_0050 = (-0.166666666666667)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0039;
    const ccomplex_t IT_0055 = 0.333333333333333*IT_0054;
    const ccomplex_t IT_0056 = IT_0051*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (IT_0042 + -IT_0044)*(IT_0053 + -IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = 18*s_34;
    const ccomplex_t IT_0061 = (-18)*IT_0038;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = IT_0043*IT_0053;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = conj(V_ts)*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0066 = IT_0005*IT_0065;
    const ccomplex_t IT_0067 = IT_0005*IT_0011;
    const ccomplex_t IT_0068 = m_t*conj(V_ts)*conj(V_u2)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + (-0.5)*IT_0070);
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0022 + IT_0038 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = m_t*V_ts*V_u2*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0077 = IT_0067*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = V_ts*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0082 = IT_0005*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0080 + IT_0085;
    const ccomplex_t IT_0087 = IT_0075*IT_0086;
    const ccomplex_t IT_0088 = m_t*V_ts*V_u1*e_em*IT_0000*conj(U_st_10);
    const ccomplex_t IT_0089 = IT_0067*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = V_ts*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0094 = IT_0005*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = IT_0092 + IT_0097;
    const ccomplex_t IT_0099 = conj(V_ts)*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0100 = IT_0005*IT_0099;
    const ccomplex_t IT_0101 = m_t*conj(V_ts)*conj(V_u1)*e_em*IT_0000*U_st_10;
    const ccomplex_t IT_0102 = IT_0067*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + (-0.5)*IT_0103);
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0022 + IT_0038 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0098*IT_0108;
    const ccomplex_t IT_0110 = IT_0087 + IT_0109;
    const ccomplex_t IT_0111 = IT_0043*IT_0057;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = IT_0114*IT_0116;
    const ccomplex_t IT_0118 = IT_0073*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0000*IT_0002*IT_0005*conj(U_st_00);
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_ts)*e_em*IT_0000*IT_0002*IT_0005*U_st_00;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = IT_0106*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = IT_0119 + IT_0126;
    const ccomplex_t IT_0128 = s_34*IT_0022;
    const ccomplex_t IT_0129 = (-6)*IT_0128;
    const ccomplex_t IT_0130 = s_13*s_14;
    const ccomplex_t IT_0131 = 12*IT_0130;
    const ccomplex_t IT_0132 = IT_0129 + IT_0131;
    const ccomplex_t IT_0133 = (-18)*IT_0128;
    const ccomplex_t IT_0134 = 36*IT_0130;
    const ccomplex_t IT_0135 = IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0137 = (-2)*IT_0136;
    const ccomplex_t IT_0138 = cpow(s_12 + IT_0022 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0139 = IT_0136*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = IT_0137*IT_0140;
    const ccomplex_t IT_0142 = 0.25*IT_0141;
    const ccomplex_t IT_0143 = cpow(IT_0142, 2);
    const ccomplex_t IT_0144 = (-36)*IT_0128;
    const ccomplex_t IT_0145 = IT_0022*IT_0038;
    const ccomplex_t IT_0146 = (-36)*IT_0145;
    const ccomplex_t IT_0147 = 72*IT_0130;
    const ccomplex_t IT_0148 = IT_0144 + IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0150 = (-0.333333333333333)*IT_0149;
    const ccomplex_t IT_0151 = (-1.33333333333333)*IT_0149;
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = IT_0138*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = (-0.0833333333333333)*IT_0141 + 0.5*IT_0154;
    const ccomplex_t IT_0156 = (-18)*IT_0059;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0158 = IT_0086*IT_0116;
    const ccomplex_t IT_0159 = IT_0073*IT_0157*IT_0158;
    const ccomplex_t IT_0160 = m_s*IT_0080;
    const ccomplex_t IT_0161 = m_s*IT_0084;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0160 + IT_0162;
    const ccomplex_t IT_0164 = IT_0075*IT_0163;
    const ccomplex_t IT_0165 = m_s*IT_0092;
    const ccomplex_t IT_0166 = m_s*IT_0096;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = IT_0108*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0171 = IT_0098*IT_0123;
    const ccomplex_t IT_0172 = IT_0106*IT_0170*IT_0171;
    const ccomplex_t IT_0173 = -IT_0159 + -IT_0164 + -IT_0169 + -IT_0172;
    const ccomplex_t IT_0174 = (-6)*IT_0173;
    const ccomplex_t IT_0175 = IT_0156 + IT_0174;
    const ccomplex_t IT_0176 = (-18)*conj(IT_0059);
    const ccomplex_t IT_0177 = (-6)*conj(IT_0173);
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = IT_0105*IT_0121;
    const ccomplex_t IT_0180 = IT_0106*IT_0170*IT_0179;
    const ccomplex_t IT_0181 = IT_0072*IT_0114;
    const ccomplex_t IT_0182 = IT_0073*IT_0157*IT_0181;
    const ccomplex_t IT_0183 = m_s*IT_0126;
    const ccomplex_t IT_0184 = m_s*IT_0119;
    const ccomplex_t IT_0185 = -IT_0180 + -IT_0182 + -IT_0183 + -IT_0184;
    const ccomplex_t IT_0186 = (-18)*IT_0173;
    const ccomplex_t IT_0187 = (-6)*IT_0059;
    const ccomplex_t IT_0188 = IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = (-18)*conj(IT_0173);
    const ccomplex_t IT_0190 = (-6)*conj(IT_0059);
    const ccomplex_t IT_0191 = IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = 6*s_34;
    const ccomplex_t IT_0193 = -IT_0038;
    const ccomplex_t IT_0194 = 6*IT_0037;
    const ccomplex_t IT_0195 = 6*conj(IT_0037);
    const ccomplex_t IT_0196 = m_s*e_em*IT_0000*IT_0002*IT_0005*IT_0010;
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = m_t*e_em*mu_h*IT_0000*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0199 = IT_0005*IT_0198;
    const ccomplex_t IT_0200 = m_t*e_em*mu_h*IT_0000*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0201 = IT_0005*IT_0200;
    const ccomplex_t IT_0202 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0203 = IT_0001*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0206 = IT_0001*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (-0.5)*IT_0199 + 0.5*IT_0201 + (-0.5)*IT_0204 +
       0.5*IT_0207;
    const ccomplex_t IT_0209 = IT_0197*IT_0208;
    const ccomplex_t IT_0210 = cpow((-2)*s_12 + (-2)*IT_0022 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0211 = IT_0209*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = m_s*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0001*IT_0005*IT_0011;
    const ccomplex_t IT_0216 = IT_0198*IT_0215;
    const ccomplex_t IT_0217 = IT_0200*IT_0215;
    const ccomplex_t IT_0218 = IT_0010*IT_0202;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = IT_0010*IT_0205;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = 0.5*IT_0216 + (-0.5)*IT_0217 + (-0.5)*IT_0219 +
       0.5*IT_0221;
    const ccomplex_t IT_0223 = IT_0214*IT_0222;
    const ccomplex_t IT_0224 = pow(m_Z, 2);
    const ccomplex_t IT_0225 = cpow((-2)*s_12 + (-2)*IT_0022 + IT_0224 + 
      -reg_prop, -1);
    const ccomplex_t IT_0226 = IT_0223*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = -IT_0212 + -IT_0227;
    const ccomplex_t IT_0229 = s_13 + s_14;
    const ccomplex_t IT_0230 = m_s*IT_0229;
    const ccomplex_t IT_0231 = 18*IT_0038;
    const ccomplex_t IT_0232 = IT_0060 + IT_0231;
    const ccomplex_t IT_0233 = (-18)*s_34 + IT_0061;
    const ccomplex_t IT_0234 = s_34 + IT_0038;
    const ccomplex_t IT_0235 = 6*conj(IT_0185);
    const ccomplex_t IT_0236 = (36*s_34 + 36*IT_0038)*conj(IT_0228) + ((-18)
      *conj(IT_0064) + (-6)*conj(IT_0110) + 18*conj(IT_0112) + 6*conj(IT_0127))
      *IT_0230 + conj(IT_0058)*IT_0232 + conj(IT_0059)*IT_0233 + IT_0234*
      (IT_0177 + IT_0235);
    const ccomplex_t IT_0237 = (-6)*IT_0145;
    const ccomplex_t IT_0238 = IT_0129 + IT_0131 + IT_0237;
    const ccomplex_t IT_0239 = conj(IT_0064) + conj(IT_0112);
    const ccomplex_t IT_0240 = conj(IT_0110) + conj(IT_0127);
    const ccomplex_t IT_0241 = IT_0110 + IT_0127;
    const ccomplex_t IT_0242 = -IT_0155;
    const ccomplex_t IT_0243 = IT_0064 + IT_0112;
    const ccomplex_t IT_0244 = -IT_0142;
    const ccomplex_t IT_0245 = (-18)*IT_0145;
    const ccomplex_t IT_0246 = IT_0133 + IT_0134 + IT_0245;
    const ccomplex_t IT_0247 = 6*IT_0185;
    const ccomplex_t IT_0248 = 6*IT_0127;
    const ccomplex_t IT_0249 = m_s*s_13;
    const ccomplex_t IT_0250 = conj(IT_0064) + 3*conj(IT_0110);
    const ccomplex_t IT_0251 = 6*conj(IT_0173);
    const ccomplex_t IT_0252 = 18*conj(IT_0059);
    const ccomplex_t IT_0253 = 6*conj(IT_0059);
    const ccomplex_t IT_0254 = 18*conj(IT_0173);
    const ccomplex_t IT_0255 = 18*conj(IT_0185);
    const ccomplex_t IT_0256 = 18*conj(IT_0058);
    const ccomplex_t IT_0257 = IT_0059*(18*conj(IT_0064) + 6*conj(IT_0110)) +
       6*conj(IT_0127)*(IT_0058 + 3*IT_0185) + conj(IT_0112)*(18*IT_0058 +
       IT_0247) + 6*IT_0173*IT_0250 + IT_0064*(IT_0251 + IT_0252) + IT_0110*
      (IT_0253 + IT_0254) + IT_0248*(conj(IT_0058) + 0.166666666666667*IT_0255) 
      + IT_0112*(IT_0235 + IT_0256);
    const ccomplex_t IT_0258 = -s_14;
    const ccomplex_t IT_0259 = s_13 + IT_0258;
    const ccomplex_t IT_0260 = m_s*IT_0259;
    const ccomplex_t IT_0261 = 18*IT_0037;
    const ccomplex_t IT_0262 = 18*conj(IT_0037);
    const ccomplex_t IT_0263 = (-18)*IT_0185;
    const ccomplex_t IT_0264 = (-6)*IT_0185;
    const ccomplex_t IT_0265 = IT_0194*IT_0240 + IT_0195*IT_0241 + IT_0239
      *IT_0261 + IT_0243*IT_0262 + IT_0142*((-12)*IT_0037 + 12*conj(IT_0037) + (
      -6)*IT_0058 + 6*conj(IT_0058) + IT_0186 + IT_0187 + IT_0253 + IT_0254 +
       IT_0255 + IT_0263) + IT_0155*((-36)*IT_0037 + 36*conj(IT_0037) + (-18)
      *IT_0058 + IT_0156 + IT_0174 + IT_0235 + IT_0251 + IT_0252 + IT_0256 +
       IT_0264);
    const ccomplex_t IT_0266 = m_s*s_14;
    const ccomplex_t IT_0267 = IT_0064 + 3*IT_0110;
    const ccomplex_t IT_0268 = (-18)*IT_0064;
    const ccomplex_t IT_0269 = (-6)*IT_0110;
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = (-18)*conj(IT_0064);
    const ccomplex_t IT_0272 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0273 = IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = (-0.166666666666667)*IT_0112;
    const ccomplex_t IT_0275 = (-0.166666666666667)*conj(IT_0112);
    const ccomplex_t IT_0276 = (-0.166666666666667)*conj(IT_0127)*IT_0188 + (
      -0.166666666666667)*IT_0127*IT_0191 + (-0.166666666666667)*conj(IT_0110)
      *IT_0263 + (-0.166666666666667)*conj(IT_0064)*IT_0264 + conj(IT_0185)
      *IT_0267 + (-0.166666666666667)*conj(IT_0058)*IT_0270 + (
      -0.166666666666667)*IT_0058*IT_0273 + IT_0178*IT_0274 + IT_0175*IT_0275;
    const ccomplex_t IT_0277 = IT_0037*conj(IT_0037)*(36*s_34 + (-36)*IT_0038)
       + (conj(IT_0037)*(IT_0058 + IT_0059) + IT_0037*(conj(IT_0058) + conj
      (IT_0059)))*IT_0062 + (conj(IT_0064)*IT_0110 + IT_0064*conj(IT_0110) +
       conj(IT_0112)*IT_0127 + IT_0112*conj(IT_0127))*IT_0132 + (IT_0064*conj
      (IT_0064) + IT_0110*conj(IT_0110) + IT_0112*conj(IT_0112) + IT_0127*conj
      (IT_0127))*IT_0135 + -IT_0143*IT_0148 + (-2)*IT_0155*(IT_0142*((-12)
      *IT_0128 + 24*IT_0130 + (-12)*IT_0145) + 0.5*IT_0148*IT_0155) + IT_0060*
      (IT_0058*conj(IT_0058) + IT_0059*conj(IT_0059) + IT_0173*conj(IT_0173) +
       IT_0185*conj(IT_0185)) + IT_0038*(conj(IT_0058)*IT_0175 + IT_0058*IT_0178
       + conj(IT_0185)*IT_0188 + IT_0185*IT_0191) + (conj(IT_0059)*IT_0173 +
       IT_0059*conj(IT_0173) + conj(IT_0058)*IT_0185 + IT_0058*conj(IT_0185))
      *IT_0192 + (s_34 + IT_0193)*((conj(IT_0173) + conj(IT_0185))*IT_0194 + 
      (IT_0173 + IT_0185)*IT_0195) + IT_0228*IT_0236 + IT_0238*(IT_0142*IT_0239 
      + IT_0155*IT_0240 + IT_0241*IT_0242 + IT_0243*IT_0244) + (IT_0155*IT_0239 
      + IT_0142*IT_0240 + IT_0242*IT_0243 + IT_0241*IT_0244)*IT_0246 + conj
      (IT_0228)*(IT_0058*IT_0232 + IT_0059*IT_0233 + IT_0234*(IT_0174 + IT_0247)
       + IT_0230*((-18)*IT_0064 + (-6)*IT_0110 + 18*IT_0112 + IT_0248)) +
       IT_0249*IT_0257 + IT_0260*IT_0265 + (-6)*IT_0145*(IT_0127*IT_0250 + conj
      (IT_0127)*IT_0267 + IT_0273*IT_0274 + IT_0270*IT_0275) + (-6)*IT_0266
      *IT_0276;
    return create_ccomplex_return(IT_0277);
}

