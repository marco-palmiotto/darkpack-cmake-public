#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_H0, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = U_d1*V_Wp1;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = V_u1*U_Wm1;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (IT_0004 + IT_0007)*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0009*IT_0014;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0018 + -IT_0020 + (-2)
      *IT_0022);
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = U_d1*V_Wp2;
    const ccomplex_t IT_0029 = IT_0002*IT_0028;
    const ccomplex_t IT_0030 = V_u2*U_Wm1;
    const ccomplex_t IT_0031 = IT_0005*IT_0030;
    const ccomplex_t IT_0032 = IT_0009*(IT_0029 + IT_0031);
    const ccomplex_t IT_0033 = 1.4142135623731*e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0037 = IT_0016*IT_0036;
    const ccomplex_t IT_0038 = IT_0019*IT_0036;
    const ccomplex_t IT_0039 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0040 = IT_0019*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0038 + (-2)
      *IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = pow(m_C_1, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0000 + IT_0044 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = sin(beta);
    const ccomplex_t IT_0049 = cos(beta);
    const ccomplex_t IT_0050 = e_em*IT_0009*IT_0015*(IT_0002*IT_0048 + 
      -IT_0005*IT_0049);
    const ccomplex_t IT_0051 = IT_0006*IT_0049;
    const ccomplex_t IT_0052 = IT_0003*IT_0048;
    const ccomplex_t IT_0053 = IT_0009*(IT_0051 + IT_0052);
    const ccomplex_t IT_0054 = 1.4142135623731*e_em*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + (-2)*IT_0044 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = e_em*IT_0009*IT_0015*(IT_0005*IT_0048 + IT_0002
      *IT_0049);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0003*IT_0049;
    const ccomplex_t IT_0063 = IT_0006*IT_0048;
    const ccomplex_t IT_0064 = IT_0009*(IT_0062 + -IT_0063);
    const ccomplex_t IT_0065 = 1.4142135623731*e_em*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = pow(m_Z, 2);
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + (-2)*IT_0044 + IT_0067 + 
      -reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0061*IT_0070;
    const ccomplex_t IT_0072 = (-2)*IT_0027 + (-2)*IT_0047 + -IT_0059 + 
      -IT_0071;
    const ccomplex_t IT_0073 = pow(m_Z, -2);
    const ccomplex_t IT_0074 = pow(s_14, 2);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = -IT_0044;
    const ccomplex_t IT_0077 = IT_0075 + IT_0076;
    const ccomplex_t IT_0078 = s_12*IT_0077;
    const ccomplex_t IT_0079 = IT_0027 + IT_0047;
    const ccomplex_t IT_0080 = s_24*IT_0067;
    const ccomplex_t IT_0081 = IT_0073*IT_0080;
    const ccomplex_t IT_0082 = (-0.25)*IT_0081;
    const ccomplex_t IT_0083 = s_24 + IT_0082;
    const ccomplex_t IT_0084 = m_C_1*IT_0083;
    const ccomplex_t IT_0085 = (-4)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0087 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0088 = IT_0016*IT_0087;
    const ccomplex_t IT_0089 = IT_0019*IT_0087;
    const ccomplex_t IT_0090 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0091 = IT_0019*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + -IT_0089 + (-2)
      *IT_0091);
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = IT_0035*IT_0093;
    const ccomplex_t IT_0095 = IT_0045*IT_0086*IT_0094;
    const ccomplex_t IT_0096 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0097 = IT_0002*IT_0096;
    const ccomplex_t IT_0098 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0099 = IT_0005*IT_0098;
    const ccomplex_t IT_0100 = IT_0009*(IT_0097 + IT_0099);
    const ccomplex_t IT_0101 = 1.4142135623731*e_em*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0105 = IT_0016*IT_0104;
    const ccomplex_t IT_0106 = IT_0019*IT_0104;
    const ccomplex_t IT_0107 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0108 = IT_0019*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + -IT_0106 + (-2)
      *IT_0108);
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0103*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0000 + IT_0044 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0086*IT_0111*IT_0112;
    const ccomplex_t IT_0114 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0115 = IT_0002*IT_0114;
    const ccomplex_t IT_0116 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0117 = IT_0005*IT_0116;
    const ccomplex_t IT_0118 = IT_0009*(IT_0115 + IT_0117);
    const ccomplex_t IT_0119 = 1.4142135623731*e_em*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0123 = IT_0016*IT_0122;
    const ccomplex_t IT_0124 = IT_0019*IT_0122;
    const ccomplex_t IT_0125 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0126 = IT_0019*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + -IT_0124 + (-2)
      *IT_0126);
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = IT_0121*IT_0128;
    const ccomplex_t IT_0130 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0132 = IT_0129*IT_0130*IT_0131;
    const ccomplex_t IT_0133 = IT_0013*IT_0128;
    const ccomplex_t IT_0134 = IT_0001*IT_0131*IT_0133;
    const ccomplex_t IT_0135 = m_C_1*IT_0027;
    const ccomplex_t IT_0136 = m_C_1*IT_0047;
    const ccomplex_t IT_0137 = cpow(IT_0014, -4);
    const ccomplex_t IT_0138 = IT_0002*IT_0008*IT_0049*IT_0137;
    const ccomplex_t IT_0139 = IT_0005*IT_0008*IT_0048*IT_0137;
    const ccomplex_t IT_0140 = cpow(IT_0008, -3);
    const ccomplex_t IT_0141 = IT_0002*IT_0049*IT_0140;
    const ccomplex_t IT_0142 = IT_0005*IT_0048*IT_0140;
    const ccomplex_t IT_0143 = cpow(IT_0014, -2);
    const ccomplex_t IT_0144 = IT_0002*IT_0009*IT_0049*IT_0143;
    const ccomplex_t IT_0145 = IT_0005*IT_0009*IT_0048*IT_0143;
    const ccomplex_t IT_0146 = cpow(IT_0008, -2);
    const ccomplex_t IT_0147 = cpow(IT_0143 + IT_0146, -1);
    const ccomplex_t IT_0148 = (IT_0138 + IT_0139 + IT_0141 + IT_0142 + 2
      *IT_0144 + 2*IT_0145)*IT_0147;
    const ccomplex_t IT_0149 = m_W*e_em;
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0128*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0044 + -reg_prop, -1);
    const ccomplex_t IT_0154 = IT_0152*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = IT_0024*IT_0121;
    const ccomplex_t IT_0157 = IT_0130*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = m_C_1*IT_0158;
    const ccomplex_t IT_0160 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0161 = IT_0016*IT_0160;
    const ccomplex_t IT_0162 = IT_0019*IT_0160;
    const ccomplex_t IT_0163 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0164 = IT_0019*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0161 + -IT_0162 + (-2)
      *IT_0164);
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = IT_0103*IT_0166;
    const ccomplex_t IT_0168 = IT_0112*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = m_C_1*IT_0169;
    const ccomplex_t IT_0171 = IT_0095 + IT_0113 + IT_0132 + IT_0134 + IT_0135
       + IT_0136 + -IT_0155 + IT_0159 + IT_0170;
    const ccomplex_t IT_0172 = 2*IT_0171;
    const ccomplex_t IT_0173 = s_14*IT_0067;
    const ccomplex_t IT_0174 = IT_0073*IT_0173;
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = s_14 + IT_0175;
    const ccomplex_t IT_0177 = m_C_1*IT_0176;
    const ccomplex_t IT_0178 = 2*IT_0177;
    const ccomplex_t IT_0179 = IT_0001*IT_0131*IT_0156;
    const ccomplex_t IT_0180 = IT_0025*IT_0130*IT_0131;
    const ccomplex_t IT_0181 = IT_0024*IT_0151;
    const ccomplex_t IT_0182 = IT_0153*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0185 = IT_0002*IT_0184;
    const ccomplex_t IT_0186 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0187 = IT_0005*IT_0186;
    const ccomplex_t IT_0188 = IT_0009*(IT_0185 + IT_0187);
    const ccomplex_t IT_0189 = 1.4142135623731*e_em*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = IT_0042*IT_0191;
    const ccomplex_t IT_0193 = IT_0045*IT_0086*IT_0192;
    const ccomplex_t IT_0194 = U_d2*V_Wp1;
    const ccomplex_t IT_0195 = IT_0002*IT_0194;
    const ccomplex_t IT_0196 = V_u1*U_Wm2;
    const ccomplex_t IT_0197 = IT_0005*IT_0196;
    const ccomplex_t IT_0198 = IT_0009*(IT_0195 + IT_0197);
    const ccomplex_t IT_0199 = 1.4142135623731*e_em*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = IT_0166*IT_0201;
    const ccomplex_t IT_0203 = IT_0086*IT_0112*IT_0202;
    const ccomplex_t IT_0204 = IT_0130*IT_0133;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = m_C_1*IT_0205;
    const ccomplex_t IT_0207 = IT_0110*IT_0201;
    const ccomplex_t IT_0208 = IT_0112*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = m_C_1*IT_0209;
    const ccomplex_t IT_0211 = IT_0001*IT_0129;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = m_C_1*IT_0212;
    const ccomplex_t IT_0214 = IT_0093*IT_0191;
    const ccomplex_t IT_0215 = IT_0045*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = m_C_1*IT_0216;
    const ccomplex_t IT_0218 = -IT_0179 + -IT_0180 + IT_0183 + -IT_0193 + 
      -IT_0203 + -IT_0206 + -IT_0210 + -IT_0213 + -IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0218;
    const ccomplex_t IT_0220 = IT_0205 + IT_0209;
    const ccomplex_t IT_0221 = 0.5*IT_0081;
    const ccomplex_t IT_0222 = s_24 + IT_0221;
    const ccomplex_t IT_0223 = m_C_1*IT_0222;
    const ccomplex_t IT_0224 = 2*IT_0223;
    const ccomplex_t IT_0225 = (-0.25)*IT_0174;
    const ccomplex_t IT_0226 = s_14 + IT_0225;
    const ccomplex_t IT_0227 = m_C_1*IT_0226;
    const ccomplex_t IT_0228 = (-4)*IT_0227;
    const ccomplex_t IT_0229 = m_C_1*IT_0077;
    const ccomplex_t IT_0230 = (-2)*IT_0072;
    const ccomplex_t IT_0231 = -IT_0059 + -IT_0071 + (-2)*IT_0205 + (-2)
      *IT_0209;
    const ccomplex_t IT_0232 = pow(s_24, 2);
    const ccomplex_t IT_0233 = IT_0073*IT_0232;
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0044 + IT_0234;
    const ccomplex_t IT_0236 = s_12*IT_0235;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = 2*IT_0231;
    const ccomplex_t IT_0239 = s_14*s_24*IT_0073;
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = s_12 + IT_0240;
    const ccomplex_t IT_0242 = s_12*IT_0241;
    const ccomplex_t IT_0243 = IT_0049*IT_0114;
    const ccomplex_t IT_0244 = IT_0048*IT_0116;
    const ccomplex_t IT_0245 = IT_0009*(IT_0243 + -IT_0244);
    const ccomplex_t IT_0246 = 1.4142135623731*e_em*IT_0245;
    const ccomplex_t IT_0247 = 0.5*IT_0246;
    const ccomplex_t IT_0248 = IT_0068*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = IT_0061*IT_0249;
    const ccomplex_t IT_0251 = IT_0049*IT_0116;
    const ccomplex_t IT_0252 = IT_0048*IT_0114;
    const ccomplex_t IT_0253 = IT_0009*(IT_0251 + IT_0252);
    const ccomplex_t IT_0254 = 1.4142135623731*e_em*IT_0253;
    const ccomplex_t IT_0255 = (-0.5)*IT_0254;
    const ccomplex_t IT_0256 = IT_0056*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = IT_0050*IT_0257;
    const ccomplex_t IT_0259 = 2*IT_0212 + 2*IT_0216 + -IT_0250 + -IT_0258;
    const ccomplex_t IT_0260 = 2*IT_0259;
    const ccomplex_t IT_0261 = IT_0044*IT_0241;
    const ccomplex_t IT_0262 = IT_0044*IT_0077;
    const ccomplex_t IT_0263 = (-2)*conj(IT_0072);
    const ccomplex_t IT_0264 = 2*conj(IT_0231);
    const ccomplex_t IT_0265 = m_C_1*IT_0241;
    const ccomplex_t IT_0266 = (-2)*conj(IT_0218);
    const ccomplex_t IT_0267 = 2*IT_0158 + 2*IT_0169 + -IT_0250 + -IT_0258;
    const ccomplex_t IT_0268 = 2*IT_0267;
    const ccomplex_t IT_0269 = 2*conj(IT_0267);
    const ccomplex_t IT_0270 = m_C_1*IT_0235;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = (-2)*IT_0259;
    const ccomplex_t IT_0273 = s_12*s_14;
    const ccomplex_t IT_0274 = s_24*IT_0044;
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = IT_0273 + IT_0275;
    const ccomplex_t IT_0277 = -IT_0212 + -IT_0216;
    const ccomplex_t IT_0278 = 2*conj(IT_0072);
    const ccomplex_t IT_0279 = -IT_0158 + -IT_0169;
    const ccomplex_t IT_0280 = (-2)*conj(IT_0259);
    const ccomplex_t IT_0281 = IT_0072*conj(IT_0079) + 0.5*conj(IT_0220)
      *IT_0230 + 0.5*IT_0220*IT_0263 + conj(IT_0259)*IT_0277 + 0.5*IT_0260*conj
      (IT_0277) + 0.5*IT_0079*IT_0278 + 0.5*IT_0272*conj(IT_0279) + 0.5*IT_0279
      *IT_0280;
    const ccomplex_t IT_0282 = 2*IT_0218;
    const ccomplex_t IT_0283 = IT_0044*IT_0235;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = (-2)*IT_0231;
    const ccomplex_t IT_0286 = s_12*s_24;
    const ccomplex_t IT_0287 = s_14*IT_0044;
    const ccomplex_t IT_0288 = -IT_0287;
    const ccomplex_t IT_0289 = IT_0286 + IT_0288;
    const ccomplex_t IT_0290 = (-2)*conj(IT_0231);
    const ccomplex_t IT_0291 = (-0.5)*conj(IT_0220)*IT_0238 + (-0.5)*IT_0220
      *IT_0264 + conj(IT_0267)*IT_0277 + IT_0267*conj(IT_0277) + (-0.5)*IT_0269
      *IT_0279 + (-0.5)*IT_0268*conj(IT_0279) + (-0.5)*conj(IT_0079)*IT_0285 + (
      -0.5)*IT_0079*IT_0290;
    const ccomplex_t IT_0292 = 2*conj(IT_0218);
    const ccomplex_t IT_0293 = (-2)*IT_0171;
    const ccomplex_t IT_0294 = IT_0044*IT_0067;
    const ccomplex_t IT_0295 = (-3)*IT_0294;
    const ccomplex_t IT_0296 = (-2)*IT_0079;
    const ccomplex_t IT_0297 = 3*IT_0294;
    const ccomplex_t IT_0298 = (-2)*IT_0220;
    const ccomplex_t IT_0299 = s_12*IT_0067;
    const ccomplex_t IT_0300 = s_14*s_24;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = IT_0299 + IT_0301;
    const ccomplex_t IT_0303 = IT_0073*IT_0302;
    const ccomplex_t IT_0304 = (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = pow(m_Z, 4);
    const ccomplex_t IT_0306 = s_12*IT_0305;
    const ccomplex_t IT_0307 = IT_0073*IT_0306;
    const ccomplex_t IT_0308 = (-0.5)*IT_0307;
    const ccomplex_t IT_0309 = IT_0299 + IT_0301 + IT_0308;
    const ccomplex_t IT_0310 = (-0.25)*IT_0307;
    const ccomplex_t IT_0311 = IT_0300 + IT_0310;
    const ccomplex_t IT_0312 = (-2)*conj(IT_0079);
    const ccomplex_t IT_0313 = (-2)*conj(IT_0220);
    const ccomplex_t IT_0314 = (-0.25)*IT_0295;
    const ccomplex_t IT_0315 = (-2)*conj(IT_0171);
    const ccomplex_t IT_0316 = (-0.25)*IT_0228;
    const ccomplex_t IT_0317 = (-0.25)*IT_0224;
    const ccomplex_t IT_0318 = 2*conj(IT_0171);
    const ccomplex_t IT_0319 = (-0.25)*IT_0085;
    const ccomplex_t IT_0320 = 2*IT_0072*conj(IT_0072)*IT_0078 + conj(IT_0079)
      *(IT_0085*IT_0172 + IT_0178*IT_0219) + conj(IT_0220)*(IT_0172*IT_0224 +
       IT_0219*IT_0228) + conj(IT_0218)*IT_0229*IT_0230 + conj(IT_0231)*(IT_0237
      *IT_0238 + IT_0230*IT_0242) + conj(IT_0259)*(IT_0078*IT_0260 + IT_0238
      *IT_0261 + IT_0230*IT_0262) + IT_0231*IT_0242*IT_0263 + IT_0259*(IT_0262
      *IT_0263 + IT_0261*IT_0264) + IT_0265*(conj(IT_0218)*IT_0238 + IT_0219
      *conj(IT_0259) + IT_0218*IT_0264 + IT_0259*IT_0266 + conj(IT_0171)*
      (IT_0230 + IT_0268) + IT_0171*(IT_0263 + IT_0269)) + conj(IT_0171)*
      (IT_0238*IT_0271 + IT_0229*IT_0272) + 2*IT_0276*IT_0281 + 2*conj(IT_0267)*
      (IT_0072*IT_0261 + 0.5*IT_0237*IT_0268 + 0.5*IT_0242*IT_0272 + 0.5*IT_0271
      *IT_0282 + 0.5*IT_0284*IT_0285) + (-2)*IT_0289*IT_0291 + IT_0267*(IT_0261
      *IT_0278 + IT_0242*IT_0280 + IT_0284*IT_0290 + IT_0271*IT_0292) + conj
      (IT_0279)*(IT_0224*IT_0282 + IT_0228*IT_0293 + IT_0296*IT_0297 + IT_0295
      *IT_0298) + conj(IT_0277)*(IT_0085*IT_0282 + IT_0178*IT_0293 + IT_0295
      *IT_0296 + IT_0297*IT_0298) + 6*IT_0171*(IT_0044*conj(IT_0218) +
       0.166666666666667*IT_0264*IT_0271 + 0.166666666666667*IT_0229*IT_0280 +
       0.666666666666667*conj(IT_0171)*(s_12 + IT_0304)) + 6*IT_0218*(IT_0044
      *conj(IT_0171) + 0.166666666666667*IT_0229*IT_0263 + 0.666666666666667
      *conj(IT_0218)*(s_12 + IT_0304)) + (-4)*IT_0279*(conj(IT_0279)*IT_0309 + 2
      *conj(IT_0277)*IT_0311 + (-0.25)*IT_0297*IT_0312 + IT_0313*IT_0314 +
       IT_0315*IT_0316 + IT_0292*IT_0317) + (-4)*IT_0220*(conj(IT_0220)*IT_0309 
      + 2*conj(IT_0079)*IT_0311 + IT_0266*IT_0316 + IT_0317*IT_0318) + (-4)
      *IT_0277*(conj(IT_0277)*IT_0309 + 2*conj(IT_0279)*IT_0311 + (-0.25)
      *IT_0297*IT_0313 + IT_0312*IT_0314 + (-0.25)*IT_0178*IT_0315 + IT_0292
      *IT_0319) + (-4)*IT_0079*((-0.25)*IT_0178*IT_0266 + conj(IT_0079)*IT_0309 
      + 2*conj(IT_0220)*IT_0311 + IT_0318*IT_0319);
    return create_ccomplex_return(IT_0320);
}

