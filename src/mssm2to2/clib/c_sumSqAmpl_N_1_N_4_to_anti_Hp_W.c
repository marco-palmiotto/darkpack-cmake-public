#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_4_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_4_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = m_N_1*m_N_4;
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_1, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = IT_0009*IT_0017;
    const ccomplex_t IT_0019 = conj(N_u4)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 2*IT_0020;
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = conj(N_B4)*conj(V_u1)*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = conj(N_W4)*conj(V_u1)*e_em;
    const ccomplex_t IT_0029 = IT_0018*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = IT_0027 + IT_0030;
    const ccomplex_t IT_0032 = IT_0021 + IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0016*IT_0034;
    const ccomplex_t IT_0036 = pow(m_Hp, 2);
    const ccomplex_t IT_0037 = pow(m_N_4, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0036 + IT_0037 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0035*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + -IT_0004 + -IT_0037 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0044 = IT_0024*IT_0043;
    const ccomplex_t IT_0045 = sin(beta);
    const ccomplex_t IT_0046 = IT_0023*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0050 = IT_0018*IT_0049;
    const ccomplex_t IT_0051 = IT_0009*IT_0045;
    const ccomplex_t IT_0052 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0055 = IT_0024*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0057 = IT_0046*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0059 = IT_0018*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0061 = IT_0051*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0044 + (-0.5)*IT_0048 + (-0.5)*IT_0050 +
       0.5*IT_0053 + 0.5*IT_0055 + (-0.5)*IT_0057 + (-0.5)*IT_0059 + 0.5*IT_0061;
    const ccomplex_t IT_0063 = e_em*IT_0009;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = IT_0042*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0036 + IT_0037 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0068 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0069 = IT_0009*IT_0068;
    const ccomplex_t IT_0070 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0071 = IT_0009*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-0.5)*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = conj(N_B4)*conj(V_u2)*e_em;
    const ccomplex_t IT_0076 = IT_0024*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W4)*conj(V_u2)*e_em;
    const ccomplex_t IT_0079 = IT_0018*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = IT_0077 + IT_0080;
    const ccomplex_t IT_0082 = conj(N_u4)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0083 = IT_0018*IT_0082;
    const ccomplex_t IT_0084 = 2*IT_0083;
    const ccomplex_t IT_0085 = IT_0081 + IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0074*IT_0087;
    const ccomplex_t IT_0089 = IT_0067*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0004 + -IT_0037 + -reg_prop, -1);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = cos(alpha);
    const ccomplex_t IT_0094 = sin(alpha);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*e_em*IT_0009*(IT_0017
      *IT_0093 + IT_0045*IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0023*IT_0093;
    const ccomplex_t IT_0098 = IT_0043*IT_0097;
    const ccomplex_t IT_0099 = IT_0023*IT_0094;
    const ccomplex_t IT_0100 = IT_0047*IT_0099;
    const ccomplex_t IT_0101 = IT_0054*IT_0097;
    const ccomplex_t IT_0102 = IT_0056*IT_0099;
    const ccomplex_t IT_0103 = IT_0098 + IT_0100 + IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = IT_0009*IT_0093;
    const ccomplex_t IT_0105 = IT_0049*IT_0104;
    const ccomplex_t IT_0106 = IT_0009*IT_0094;
    const ccomplex_t IT_0107 = IT_0052*IT_0106;
    const ccomplex_t IT_0108 = IT_0058*IT_0104;
    const ccomplex_t IT_0109 = IT_0060*IT_0106;
    const ccomplex_t IT_0110 = -IT_0105 + -IT_0107 + -IT_0108 + -IT_0109;
    const ccomplex_t IT_0111 = IT_0103 + IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0096*IT_0113;
    const ccomplex_t IT_0115 = IT_0092*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*e_em*IT_0009*(IT_0045
      *IT_0093 + -IT_0017*IT_0094);
    const ccomplex_t IT_0117 = IT_0047*IT_0097;
    const ccomplex_t IT_0118 = IT_0049*IT_0106;
    const ccomplex_t IT_0119 = IT_0056*IT_0097;
    const ccomplex_t IT_0120 = IT_0058*IT_0106;
    const ccomplex_t IT_0121 = IT_0117 + IT_0118 + IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = IT_0043*IT_0099;
    const ccomplex_t IT_0123 = IT_0052*IT_0104;
    const ccomplex_t IT_0124 = IT_0054*IT_0099;
    const ccomplex_t IT_0125 = IT_0060*IT_0104;
    const ccomplex_t IT_0126 = -IT_0122 + -IT_0123 + -IT_0124 + -IT_0125;
    const ccomplex_t IT_0127 = IT_0121 + IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = IT_0116*IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_12 + -IT_0004 + -IT_0037 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = 2*IT_0040 + -IT_0066 + 2*IT_0090 + -IT_0115 + 
      -IT_0133;
    const ccomplex_t IT_0135 = N_d4*e_em*U_Wm1;
    const ccomplex_t IT_0136 = IT_0051*IT_0135;
    const ccomplex_t IT_0137 = (-2)*IT_0136;
    const ccomplex_t IT_0138 = N_B4*U_d1*e_em;
    const ccomplex_t IT_0139 = IT_0046*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = N_W4*U_d1*e_em;
    const ccomplex_t IT_0142 = IT_0051*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = IT_0140 + IT_0143;
    const ccomplex_t IT_0145 = IT_0137 + IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0149 = IT_0009*IT_0148;
    const ccomplex_t IT_0150 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0151 = IT_0009*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + 0.5*IT_0152);
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = IT_0147*IT_0154;
    const ccomplex_t IT_0156 = IT_0038*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = N_d4*e_em*U_Wm2;
    const ccomplex_t IT_0159 = IT_0051*IT_0158;
    const ccomplex_t IT_0160 = (-2)*IT_0159;
    const ccomplex_t IT_0161 = N_B4*U_d2*e_em;
    const ccomplex_t IT_0162 = IT_0046*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = N_W4*U_d2*e_em;
    const ccomplex_t IT_0165 = IT_0051*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = IT_0163 + IT_0166;
    const ccomplex_t IT_0168 = IT_0160 + IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = 0.5*IT_0169;
    const ccomplex_t IT_0171 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0172 = IT_0009*IT_0171;
    const ccomplex_t IT_0173 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0174 = IT_0009*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0172 + 0.5*IT_0175);
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0170*IT_0177;
    const ccomplex_t IT_0179 = IT_0067*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = 2*IT_0157 + 2*IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0181;
    const ccomplex_t IT_0183 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0184 = IT_0024*IT_0183;
    const ccomplex_t IT_0185 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0186 = IT_0046*IT_0185;
    const ccomplex_t IT_0187 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0188 = IT_0018*IT_0187;
    const ccomplex_t IT_0189 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0190 = IT_0051*IT_0189;
    const ccomplex_t IT_0191 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0192 = IT_0024*IT_0191;
    const ccomplex_t IT_0193 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0194 = IT_0046*IT_0193;
    const ccomplex_t IT_0195 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0196 = IT_0018*IT_0195;
    const ccomplex_t IT_0197 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0198 = IT_0051*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0184 + 0.5*IT_0186 + 0.5*IT_0188 + (
      -0.5)*IT_0190 + (-0.5)*IT_0192 + 0.5*IT_0194 + 0.5*IT_0196 + (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0064*IT_0199;
    const ccomplex_t IT_0201 = IT_0042*IT_0200;
    const ccomplex_t IT_0202 = IT_0097*IT_0185;
    const ccomplex_t IT_0203 = IT_0106*IT_0187;
    const ccomplex_t IT_0204 = IT_0097*IT_0193;
    const ccomplex_t IT_0205 = IT_0106*IT_0195;
    const ccomplex_t IT_0206 = IT_0202 + IT_0203 + IT_0204 + IT_0205;
    const ccomplex_t IT_0207 = IT_0099*IT_0183;
    const ccomplex_t IT_0208 = IT_0104*IT_0189;
    const ccomplex_t IT_0209 = IT_0099*IT_0191;
    const ccomplex_t IT_0210 = IT_0104*IT_0197;
    const ccomplex_t IT_0211 = -IT_0207 + -IT_0208 + -IT_0209 + -IT_0210;
    const ccomplex_t IT_0212 = IT_0206 + IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0116*IT_0214;
    const ccomplex_t IT_0216 = IT_0132*IT_0215;
    const ccomplex_t IT_0217 = IT_0097*IT_0183;
    const ccomplex_t IT_0218 = IT_0099*IT_0185;
    const ccomplex_t IT_0219 = IT_0097*IT_0191;
    const ccomplex_t IT_0220 = IT_0099*IT_0193;
    const ccomplex_t IT_0221 = IT_0217 + IT_0218 + IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = IT_0104*IT_0187;
    const ccomplex_t IT_0223 = IT_0106*IT_0189;
    const ccomplex_t IT_0224 = IT_0104*IT_0195;
    const ccomplex_t IT_0225 = IT_0106*IT_0197;
    const ccomplex_t IT_0226 = -IT_0222 + -IT_0223 + -IT_0224 + -IT_0225;
    const ccomplex_t IT_0227 = IT_0221 + IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = IT_0096*IT_0229;
    const ccomplex_t IT_0231 = IT_0092*IT_0230;
    const ccomplex_t IT_0232 = -IT_0201 + -IT_0216 + -IT_0231;
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = IT_0182 + IT_0233;
    const ccomplex_t IT_0235 = 2*conj(IT_0181);
    const ccomplex_t IT_0236 = (-2)*conj(IT_0232);
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0239 = IT_0051*IT_0238;
    const ccomplex_t IT_0240 = (-2)*IT_0239;
    const ccomplex_t IT_0241 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0242 = IT_0046*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0245 = IT_0051*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = IT_0243 + IT_0246;
    const ccomplex_t IT_0248 = IT_0240 + IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = 0.5*IT_0249;
    const ccomplex_t IT_0251 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0252 = IT_0009*IT_0251;
    const ccomplex_t IT_0253 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0254 = IT_0009*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*(IT_0252 + 0.5*IT_0255);
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0250*IT_0257;
    const ccomplex_t IT_0259 = cpow((-2)*s_13 + IT_0004 + IT_0036 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0260 = IT_0258*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = cpow((-2)*s_13 + IT_0004 + IT_0036 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0263 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0264 = IT_0009*IT_0263;
    const ccomplex_t IT_0265 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0266 = IT_0009*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0264 + 0.5*IT_0267);
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0271 = IT_0051*IT_0270;
    const ccomplex_t IT_0272 = (-2)*IT_0271;
    const ccomplex_t IT_0273 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0274 = IT_0046*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0277 = IT_0051*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = IT_0275 + IT_0278;
    const ccomplex_t IT_0280 = IT_0272 + IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = 0.5*IT_0281;
    const ccomplex_t IT_0283 = IT_0269*IT_0282;
    const ccomplex_t IT_0284 = IT_0262*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = IT_0261 + IT_0285;
    const ccomplex_t IT_0287 = pow(m_W, 2);
    const ccomplex_t IT_0288 = s_14*IT_0287;
    const ccomplex_t IT_0289 = IT_0001*IT_0288;
    const ccomplex_t IT_0290 = (-0.25)*IT_0289;
    const ccomplex_t IT_0291 = s_14 + IT_0290;
    const ccomplex_t IT_0292 = m_N_4*IT_0291;
    const ccomplex_t IT_0293 = (-4)*IT_0292;
    const ccomplex_t IT_0294 = m_N_1*IT_0090;
    const ccomplex_t IT_0295 = m_N_1*IT_0040;
    const ccomplex_t IT_0296 = m_N_4*IT_0285;
    const ccomplex_t IT_0297 = m_N_4*IT_0261;
    const ccomplex_t IT_0298 = IT_0294 + IT_0295 + IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = 2*conj(IT_0298);
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0301 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0302 = IT_0009*IT_0301;
    const ccomplex_t IT_0303 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0304 = IT_0009*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*(IT_0302 + (-0.5)*IT_0305);
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = IT_0250*IT_0307;
    const ccomplex_t IT_0309 = IT_0259*IT_0300*IT_0308;
    const ccomplex_t IT_0310 = IT_0087*IT_0177;
    const ccomplex_t IT_0311 = IT_0067*IT_0300*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0313 = IT_0034*IT_0154;
    const ccomplex_t IT_0314 = IT_0038*IT_0312*IT_0313;
    const ccomplex_t IT_0315 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0316 = IT_0009*IT_0315;
    const ccomplex_t IT_0317 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0318 = IT_0009*IT_0317;
    const ccomplex_t IT_0319 = 1.4142135623731*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*(IT_0316 + (-0.5)*IT_0319);
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0282*IT_0321;
    const ccomplex_t IT_0323 = IT_0262*IT_0312*IT_0322;
    const ccomplex_t IT_0324 = -IT_0309 + -IT_0311 + -IT_0314 + -IT_0323;
    const ccomplex_t IT_0325 = 2*conj(IT_0324);
    const ccomplex_t IT_0326 = IT_0299 + IT_0325;
    const ccomplex_t IT_0327 = s_24*IT_0287;
    const ccomplex_t IT_0328 = IT_0001*IT_0327;
    const ccomplex_t IT_0329 = 0.5*IT_0328;
    const ccomplex_t IT_0330 = s_24 + IT_0329;
    const ccomplex_t IT_0331 = m_N_1*IT_0330;
    const ccomplex_t IT_0332 = 2*IT_0331;
    const ccomplex_t IT_0333 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0334 = IT_0018*IT_0333;
    const ccomplex_t IT_0335 = 2*IT_0334;
    const ccomplex_t IT_0336 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0337 = IT_0024*IT_0336;
    const ccomplex_t IT_0338 = 1.4142135623731*IT_0337;
    const ccomplex_t IT_0339 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0340 = IT_0018*IT_0339;
    const ccomplex_t IT_0341 = 1.4142135623731*IT_0340;
    const ccomplex_t IT_0342 = IT_0338 + IT_0341;
    const ccomplex_t IT_0343 = IT_0335 + IT_0342;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*IT_0343;
    const ccomplex_t IT_0345 = 0.5*IT_0344;
    const ccomplex_t IT_0346 = IT_0269*IT_0345;
    const ccomplex_t IT_0347 = IT_0262*IT_0312*IT_0346;
    const ccomplex_t IT_0348 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0349 = IT_0018*IT_0348;
    const ccomplex_t IT_0350 = 2*IT_0349;
    const ccomplex_t IT_0351 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0352 = IT_0024*IT_0351;
    const ccomplex_t IT_0353 = 1.4142135623731*IT_0352;
    const ccomplex_t IT_0354 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0355 = IT_0018*IT_0354;
    const ccomplex_t IT_0356 = 1.4142135623731*IT_0355;
    const ccomplex_t IT_0357 = IT_0353 + IT_0356;
    const ccomplex_t IT_0358 = IT_0350 + IT_0357;
    const ccomplex_t IT_0359 = (0 + _Complex_I*1)*IT_0358;
    const ccomplex_t IT_0360 = 0.5*IT_0359;
    const ccomplex_t IT_0361 = IT_0257*IT_0360;
    const ccomplex_t IT_0362 = IT_0259*IT_0300*IT_0361;
    const ccomplex_t IT_0363 = IT_0016*IT_0147;
    const ccomplex_t IT_0364 = IT_0038*IT_0312*IT_0363;
    const ccomplex_t IT_0365 = IT_0074*IT_0170;
    const ccomplex_t IT_0366 = IT_0067*IT_0300*IT_0365;
    const ccomplex_t IT_0367 = IT_0347 + IT_0362 + IT_0364 + IT_0366;
    const ccomplex_t IT_0368 = (-2)*conj(IT_0367);
    const ccomplex_t IT_0369 = m_N_1*IT_0157;
    const ccomplex_t IT_0370 = m_N_1*IT_0180;
    const ccomplex_t IT_0371 = IT_0321*IT_0345;
    const ccomplex_t IT_0372 = IT_0262*IT_0371;
    const ccomplex_t IT_0373 = (0 + _Complex_I*1)*IT_0372;
    const ccomplex_t IT_0374 = m_N_4*IT_0373;
    const ccomplex_t IT_0375 = IT_0307*IT_0360;
    const ccomplex_t IT_0376 = IT_0259*IT_0375;
    const ccomplex_t IT_0377 = (0 + _Complex_I*1)*IT_0376;
    const ccomplex_t IT_0378 = m_N_4*IT_0377;
    const ccomplex_t IT_0379 = -IT_0369 + -IT_0370 + -IT_0374 + -IT_0378;
    const ccomplex_t IT_0380 = (-2)*conj(IT_0379);
    const ccomplex_t IT_0381 = IT_0368 + IT_0380;
    const ccomplex_t IT_0382 = IT_0373 + IT_0377;
    const ccomplex_t IT_0383 = -IT_0157 + -IT_0180;
    const ccomplex_t IT_0384 = (-0.25)*IT_0328;
    const ccomplex_t IT_0385 = s_24 + IT_0384;
    const ccomplex_t IT_0386 = m_N_1*IT_0385;
    const ccomplex_t IT_0387 = (-4)*IT_0386;
    const ccomplex_t IT_0388 = 2*conj(IT_0367);
    const ccomplex_t IT_0389 = 2*conj(IT_0379);
    const ccomplex_t IT_0390 = IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = 0.5*IT_0289;
    const ccomplex_t IT_0392 = s_14 + IT_0391;
    const ccomplex_t IT_0393 = m_N_4*IT_0392;
    const ccomplex_t IT_0394 = 2*IT_0393;
    const ccomplex_t IT_0395 = (-2)*conj(IT_0298);
    const ccomplex_t IT_0396 = (-2)*conj(IT_0324);
    const ccomplex_t IT_0397 = IT_0395 + IT_0396;
    const ccomplex_t IT_0398 = -IT_0040 + -IT_0090;
    const ccomplex_t IT_0399 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0400 = -IT_0399;
    const ccomplex_t IT_0401 = s_12 + IT_0400;
    const ccomplex_t IT_0402 = s_12*IT_0401;
    const ccomplex_t IT_0403 = -IT_0066 + -IT_0115 + -IT_0133;
    const ccomplex_t IT_0404 = -IT_0201 + -IT_0216 + -IT_0231 + (-2)*IT_0261 +
       (-2)*IT_0285;
    const ccomplex_t IT_0405 = (-2)*IT_0373 + (-2)*IT_0377;
    const ccomplex_t IT_0406 = 2*IT_0134;
    const ccomplex_t IT_0407 = 2*conj(IT_0134);
    const ccomplex_t IT_0408 = m_N_1*IT_0401;
    const ccomplex_t IT_0409 = conj(IT_0367) + conj(IT_0379);
    const ccomplex_t IT_0410 = 2*IT_0232;
    const ccomplex_t IT_0411 = (-2)*IT_0181;
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = IT_0367 + IT_0379;
    const ccomplex_t IT_0414 = 2*conj(IT_0232);
    const ccomplex_t IT_0415 = (-2)*conj(IT_0181);
    const ccomplex_t IT_0416 = IT_0414 + IT_0415;
    const ccomplex_t IT_0417 = conj(IT_0298) + conj(IT_0324);
    const ccomplex_t IT_0418 = IT_0298 + IT_0324;
    const ccomplex_t IT_0419 = s_12*s_14;
    const ccomplex_t IT_0420 = s_24*IT_0004;
    const ccomplex_t IT_0421 = -IT_0420;
    const ccomplex_t IT_0422 = IT_0419 + IT_0421;
    const ccomplex_t IT_0423 = conj(IT_0286) + conj(IT_0383);
    const ccomplex_t IT_0424 = IT_0286 + IT_0383;
    const ccomplex_t IT_0425 = conj(IT_0382) + conj(IT_0398);
    const ccomplex_t IT_0426 = IT_0382 + IT_0398;
    const ccomplex_t IT_0427 = IT_0000*IT_0401;
    const ccomplex_t IT_0428 = 2*IT_0404;
    const ccomplex_t IT_0429 = 2*conj(IT_0404);
    const ccomplex_t IT_0430 = s_12*s_24;
    const ccomplex_t IT_0431 = s_14*IT_0037;
    const ccomplex_t IT_0432 = -IT_0431;
    const ccomplex_t IT_0433 = IT_0430 + IT_0432;
    const ccomplex_t IT_0434 = 2*IT_0405;
    const ccomplex_t IT_0435 = (-2)*IT_0403;
    const ccomplex_t IT_0436 = IT_0434 + IT_0435;
    const ccomplex_t IT_0437 = 2*conj(IT_0405);
    const ccomplex_t IT_0438 = (-2)*conj(IT_0403);
    const ccomplex_t IT_0439 = IT_0437 + IT_0438;
    const ccomplex_t IT_0440 = m_N_4*IT_0006;
    const ccomplex_t IT_0441 = 2*IT_0367;
    const ccomplex_t IT_0442 = s_12*IT_0287;
    const ccomplex_t IT_0443 = s_14*s_24;
    const ccomplex_t IT_0444 = (-2)*IT_0443;
    const ccomplex_t IT_0445 = IT_0442 + IT_0444;
    const ccomplex_t IT_0446 = IT_0001*IT_0445;
    const ccomplex_t IT_0447 = (-0.5)*IT_0446;
    const ccomplex_t IT_0448 = s_12 + IT_0447;
    const ccomplex_t IT_0449 = 4*IT_0298;
    const ccomplex_t IT_0450 = 6*conj(IT_0367);
    const ccomplex_t IT_0451 = s_12*IT_0006;
    const ccomplex_t IT_0452 = pow(s_24, 2);
    const ccomplex_t IT_0453 = IT_0001*IT_0452;
    const ccomplex_t IT_0454 = -IT_0037;
    const ccomplex_t IT_0455 = IT_0453 + IT_0454;
    const ccomplex_t IT_0456 = IT_0000*IT_0455;
    const ccomplex_t IT_0457 = (-2)*IT_0404;
    const ccomplex_t IT_0458 = (-2)*conj(IT_0404);
    const ccomplex_t IT_0459 = m_N_1*IT_0455;
    const ccomplex_t IT_0460 = s_12*IT_0455;
    const ccomplex_t IT_0461 = m_N_4*IT_0401;
    const ccomplex_t IT_0462 = (-2)*IT_0367;
    const ccomplex_t IT_0463 = IT_0000*IT_0287;
    const ccomplex_t IT_0464 = 3*IT_0463;
    const ccomplex_t IT_0465 = (-3)*IT_0463;
    const ccomplex_t IT_0466 = pow(m_W, 4);
    const ccomplex_t IT_0467 = s_12*IT_0466;
    const ccomplex_t IT_0468 = IT_0001*IT_0467;
    const ccomplex_t IT_0469 = (-0.5)*IT_0468;
    const ccomplex_t IT_0470 = IT_0442 + IT_0444 + IT_0469;
    const ccomplex_t IT_0471 = (-0.25)*IT_0468;
    const ccomplex_t IT_0472 = IT_0443 + IT_0471;
    const ccomplex_t IT_0473 = (-2)*IT_0298;
    const ccomplex_t IT_0474 = (-2)*IT_0324;
    const ccomplex_t IT_0475 = IT_0473 + IT_0474;
    const ccomplex_t IT_0476 = (-0.5)*IT_0394;
    const ccomplex_t IT_0477 = 2*IT_0379;
    const ccomplex_t IT_0478 = IT_0441 + IT_0477;
    const ccomplex_t IT_0479 = (-0.5)*IT_0387;
    const ccomplex_t IT_0480 = (-2)*IT_0379;
    const ccomplex_t IT_0481 = IT_0462 + IT_0480;
    const ccomplex_t IT_0482 = (-0.5)*IT_0332;
    const ccomplex_t IT_0483 = 2*IT_0298;
    const ccomplex_t IT_0484 = 2*IT_0324;
    const ccomplex_t IT_0485 = IT_0483 + IT_0484;
    const ccomplex_t IT_0486 = (-0.5)*IT_0293;
    const ccomplex_t IT_0487 = IT_0007*(conj(IT_0134)*IT_0234 + IT_0134
      *IT_0237) + IT_0286*(IT_0293*IT_0326 + IT_0332*IT_0381) + (IT_0326*IT_0332
       + IT_0293*IT_0381)*IT_0382 + IT_0383*(IT_0387*IT_0390 + IT_0394*IT_0397) 
      + (IT_0390*IT_0394 + IT_0387*IT_0397)*IT_0398 + (-2)*IT_0402*(conj(IT_0134
      )*IT_0403 + IT_0134*conj(IT_0403) + (-0.5)*IT_0237*IT_0404 + (-0.5)
      *IT_0234*conj(IT_0404) + (-0.5)*conj(IT_0405)*IT_0406 + (-0.5)*IT_0405
      *IT_0407) + IT_0408*(IT_0409*IT_0412 + IT_0413*IT_0416 + IT_0406*IT_0417 +
       IT_0407*IT_0418) + IT_0422*(IT_0234*IT_0423 + IT_0237*IT_0424 + IT_0406
      *IT_0425 + IT_0407*IT_0426) + IT_0427*(IT_0237*IT_0405 + IT_0234*conj
      (IT_0405) + conj(IT_0403)*IT_0412 + IT_0403*IT_0416 + conj(IT_0134)
      *IT_0428 + IT_0134*IT_0429) + IT_0433*(IT_0423*IT_0428 + IT_0424*IT_0429 +
       IT_0425*IT_0436 + IT_0426*IT_0439) + IT_0440*(IT_0134*IT_0388 + conj
      (IT_0379)*IT_0406 + IT_0379*IT_0407 + IT_0412*IT_0417 + IT_0416*IT_0418 +
       conj(IT_0134)*IT_0441) + 4*IT_0448*(IT_0409*IT_0413 + IT_0417*(IT_0324 +
       0.25*IT_0449)) + 6*IT_0000*(IT_0413*IT_0417 + IT_0418*(conj(IT_0379) +
       0.166666666666667*IT_0450)) + (conj(IT_0181)*IT_0234 + IT_0181*IT_0236 +
       conj(IT_0134)*IT_0406 + conj(IT_0232)*IT_0410)*IT_0451 + IT_0456*(conj
      (IT_0405)*IT_0428 + IT_0405*IT_0429 + conj(IT_0403)*IT_0457 + IT_0403
      *IT_0458) + (IT_0417*IT_0436 + IT_0418*IT_0439 + IT_0409*IT_0457 + IT_0413
      *IT_0458)*IT_0459 + 2*(IT_0403*conj(IT_0403) + 0.5*conj(IT_0404)*IT_0428 +
       0.5*conj(IT_0405)*IT_0436 + 0.5*IT_0405*IT_0438)*IT_0460 + IT_0461*
      (IT_0368*IT_0403 + IT_0388*IT_0405 + conj(IT_0379)*IT_0436 + IT_0379
      *IT_0439 + conj(IT_0405)*IT_0441 + IT_0417*IT_0457 + IT_0418*IT_0458 +
       conj(IT_0403)*IT_0462) + (-2)*conj(IT_0398)*(IT_0286*IT_0464 + -IT_0383
      *IT_0465 + 2*IT_0398*IT_0470 + (-4)*IT_0382*IT_0472 + IT_0476*IT_0478 +
       IT_0475*IT_0479) + (-2)*conj(IT_0383)*(IT_0382*IT_0464 + -IT_0398*IT_0465
       + 2*IT_0383*IT_0470 + (-4)*IT_0286*IT_0472 + IT_0475*IT_0476 + IT_0478
      *IT_0479) + (-2)*conj(IT_0382)*(IT_0383*IT_0464 + -IT_0286*IT_0465 + 2
      *IT_0382*IT_0470 + (-4)*IT_0398*IT_0472 + IT_0482*IT_0485 + IT_0481
      *IT_0486) + (-2)*conj(IT_0286)*(IT_0398*IT_0464 + -IT_0382*IT_0465 + 2
      *IT_0286*IT_0470 + (-4)*IT_0383*IT_0472 + IT_0481*IT_0482 + IT_0485
      *IT_0486);
    return create_ccomplex_return(IT_0487);
}

