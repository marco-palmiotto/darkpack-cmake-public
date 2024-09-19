#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_4_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_4_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t IT_0001 = pow(m_W, -4);
    const ccomplex_t IT_0002 = pow(s_34, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = 1 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*e_em*IT_0007*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0007, -1);
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0009;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0019 = IT_0013*IT_0018;
    const ccomplex_t IT_0020 = IT_0016*IT_0018;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0015 + IT_0017 + 
      -IT_0019 + -IT_0020);
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = m_N_1*IT_0023;
    const ccomplex_t IT_0025 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = IT_0016*IT_0025;
    const ccomplex_t IT_0028 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0029 = IT_0013*IT_0028;
    const ccomplex_t IT_0030 = IT_0016*IT_0028;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0026 + IT_0027 + 
      -IT_0029 + -IT_0030);
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0011*IT_0032;
    const ccomplex_t IT_0034 = m_N_4*IT_0033;
    const ccomplex_t IT_0035 = pow(m_N_1, 2);
    const ccomplex_t IT_0036 = pow(m_N_4, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0035 + -IT_0036 + -reg_prop, -1);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (IT_0024 + -IT_0034)*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = cos(alpha);
    const ccomplex_t IT_0042 = IT_0012*IT_0041;
    const ccomplex_t IT_0043 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = sin(alpha);
    const ccomplex_t IT_0046 = IT_0012*IT_0045;
    const ccomplex_t IT_0047 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0050 = IT_0042*IT_0049;
    const ccomplex_t IT_0051 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = IT_0044 + IT_0048 + IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = IT_0009*IT_0041;
    const ccomplex_t IT_0055 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = IT_0009*IT_0045;
    const ccomplex_t IT_0058 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0061 = IT_0054*IT_0060;
    const ccomplex_t IT_0062 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0063 = IT_0057*IT_0062;
    const ccomplex_t IT_0064 = -IT_0056 + -IT_0059 + -IT_0061 + -IT_0063;
    const ccomplex_t IT_0065 = IT_0053 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = sin(beta);
    const ccomplex_t IT_0069 = cos(beta);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_W*e_em*IT_0009*(IT_0041
      *IT_0068 + -IT_0045*IT_0069);
    const ccomplex_t IT_0071 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0035 + -IT_0036 + -reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0067*IT_0073;
    const ccomplex_t IT_0075 = IT_0042*IT_0047;
    const ccomplex_t IT_0076 = IT_0055*IT_0057;
    const ccomplex_t IT_0077 = IT_0042*IT_0051;
    const ccomplex_t IT_0078 = IT_0057*IT_0060;
    const ccomplex_t IT_0079 = IT_0075 + IT_0076 + IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = IT_0043*IT_0046;
    const ccomplex_t IT_0081 = IT_0054*IT_0058;
    const ccomplex_t IT_0082 = IT_0046*IT_0049;
    const ccomplex_t IT_0083 = IT_0054*IT_0062;
    const ccomplex_t IT_0084 = -IT_0080 + -IT_0081 + -IT_0082 + -IT_0083;
    const ccomplex_t IT_0085 = IT_0079 + IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_W*e_em*IT_0009*(IT_0045
      *IT_0068 + IT_0041*IT_0069);
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + -IT_0035 + -IT_0036 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0087*IT_0091;
    const ccomplex_t IT_0093 = -IT_0074 + -IT_0092;
    const ccomplex_t IT_0094 = m_N_4*IT_0023;
    const ccomplex_t IT_0095 = m_N_1*IT_0033;
    const ccomplex_t IT_0096 = IT_0038*(IT_0094 + -IT_0095);
    const ccomplex_t IT_0097 = (-4)*IT_0096;
    const ccomplex_t IT_0098 = 0.25*IT_0097;
    const ccomplex_t IT_0099 = IT_0093 + IT_0098;
    const ccomplex_t IT_0100 = (-4)*conj(IT_0096);
    const ccomplex_t IT_0101 = 0.25*IT_0100;
    const ccomplex_t IT_0102 = conj(IT_0093) + IT_0101;
    const ccomplex_t IT_0103 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0104 = IT_0042*IT_0103;
    const ccomplex_t IT_0105 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0106 = IT_0046*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0108 = IT_0042*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0110 = IT_0046*IT_0109;
    const ccomplex_t IT_0111 = IT_0104 + IT_0106 + IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0113 = IT_0054*IT_0112;
    const ccomplex_t IT_0114 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0115 = IT_0057*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0117 = IT_0054*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0119 = IT_0057*IT_0118;
    const ccomplex_t IT_0120 = -IT_0113 + -IT_0115 + -IT_0117 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0111 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0073*IT_0123;
    const ccomplex_t IT_0125 = IT_0042*IT_0105;
    const ccomplex_t IT_0126 = IT_0057*IT_0112;
    const ccomplex_t IT_0127 = IT_0042*IT_0109;
    const ccomplex_t IT_0128 = IT_0057*IT_0116;
    const ccomplex_t IT_0129 = IT_0125 + IT_0126 + IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = IT_0046*IT_0103;
    const ccomplex_t IT_0131 = IT_0054*IT_0114;
    const ccomplex_t IT_0132 = IT_0046*IT_0107;
    const ccomplex_t IT_0133 = IT_0054*IT_0118;
    const ccomplex_t IT_0134 = -IT_0130 + -IT_0131 + -IT_0132 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0129 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0091*IT_0137;
    const ccomplex_t IT_0139 = -IT_0124 + -IT_0138;
    const ccomplex_t IT_0140 = (-4)*IT_0093;
    const ccomplex_t IT_0141 = 0.25*IT_0140;
    const ccomplex_t IT_0142 = IT_0096 + IT_0141;
    const ccomplex_t IT_0143 = (-4)*conj(IT_0093);
    const ccomplex_t IT_0144 = 0.25*IT_0143;
    const ccomplex_t IT_0145 = conj(IT_0096) + IT_0144;
    const ccomplex_t IT_0146 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0147 = IT_0009*IT_0146;
    const ccomplex_t IT_0148 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0149 = IT_0009*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 0.5*IT_0150);
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0154 = IT_0009*IT_0153;
    const ccomplex_t IT_0155 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0156 = IT_0009*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + 0.5*IT_0157);
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0152*IT_0159;
    const ccomplex_t IT_0161 = pow(m_W, 2);
    const ccomplex_t IT_0162 = cpow((-2)*s_13 + IT_0035 + IT_0161 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = (-2)*IT_0164;
    const ccomplex_t IT_0166 = pow(m_W, -2);
    const ccomplex_t IT_0167 = pow(s_13, 2);
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = pow(s_14, 2);
    const ccomplex_t IT_0170 = IT_0166*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0173 = -IT_0035 + -IT_0172;
    const ccomplex_t IT_0174 = IT_0171 + IT_0173;
    const ccomplex_t IT_0175 = m_N_4*IT_0174;
    const ccomplex_t IT_0176 = 2*conj(IT_0040);
    const ccomplex_t IT_0177 = (-2)*conj(IT_0139);
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0180 = s_13*s_23*IT_0166;
    const ccomplex_t IT_0181 = s_14*s_24*IT_0166;
    const ccomplex_t IT_0182 = -IT_0180 + -IT_0181;
    const ccomplex_t IT_0183 = s_12 + IT_0179 + IT_0182;
    const ccomplex_t IT_0184 = m_N_1*IT_0183;
    const ccomplex_t IT_0185 = 2*conj(IT_0096);
    const ccomplex_t IT_0186 = (-2)*conj(IT_0093);
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = s_34*IT_0001;
    const ccomplex_t IT_0189 = s_14*s_23;
    const ccomplex_t IT_0190 = s_12*s_34;
    const ccomplex_t IT_0191 = IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = s_13*s_24;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = IT_0188*IT_0194;
    const ccomplex_t IT_0196 = 0.25*IT_0195;
    const ccomplex_t IT_0197 = s_12*IT_0161;
    const ccomplex_t IT_0198 = IT_0166*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = s_12 + IT_0196 + IT_0199;
    const ccomplex_t IT_0201 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0202 = IT_0009*IT_0201;
    const ccomplex_t IT_0203 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0204 = IT_0009*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*(IT_0202 + (-0.5)*IT_0205);
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0209 = IT_0009*IT_0208;
    const ccomplex_t IT_0210 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0211 = IT_0009*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*(IT_0209 + (-0.5)*IT_0212);
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0207*IT_0214;
    const ccomplex_t IT_0216 = cpow((-2)*s_23 + IT_0036 + IT_0161 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0217 = IT_0215*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = cpow((-2)*s_23 + IT_0036 + IT_0161 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0220 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0221 = IT_0009*IT_0220;
    const ccomplex_t IT_0222 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0223 = IT_0009*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + (-0.5)*IT_0224);
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0228 = IT_0009*IT_0227;
    const ccomplex_t IT_0229 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0230 = IT_0009*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + (-0.5)*IT_0231);
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0226*IT_0233;
    const ccomplex_t IT_0235 = IT_0219*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0218 + IT_0236;
    const ccomplex_t IT_0238 = m_N_4*IT_0237;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0241 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0242 = IT_0009*IT_0241;
    const ccomplex_t IT_0243 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0244 = IT_0009*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + 0.5*IT_0245);
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0214*IT_0247;
    const ccomplex_t IT_0249 = IT_0216*IT_0240*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0251 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0252 = IT_0009*IT_0251;
    const ccomplex_t IT_0253 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0254 = IT_0009*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*(IT_0252 + 0.5*IT_0255);
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0233*IT_0257;
    const ccomplex_t IT_0259 = IT_0219*IT_0250*IT_0258;
    const ccomplex_t IT_0260 = IT_0249 + IT_0259;
    const ccomplex_t IT_0261 = 8*IT_0093;
    const ccomplex_t IT_0262 = (-8)*IT_0096;
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 8*conj(IT_0093);
    const ccomplex_t IT_0265 = (-8)*conj(IT_0096);
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = s_23*IT_0035;
    const ccomplex_t IT_0268 = s_12*s_13;
    const ccomplex_t IT_0269 = s_23*IT_0166*IT_0167;
    const ccomplex_t IT_0270 = (-2)*IT_0269;
    const ccomplex_t IT_0271 = IT_0268 + IT_0270;
    const ccomplex_t IT_0272 = s_13*s_14*s_24*IT_0166;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = s_23*IT_0166*IT_0169;
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = s_13*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0277 = 2*IT_0276;
    const ccomplex_t IT_0278 = IT_0267 + IT_0271 + IT_0273 + IT_0275 + IT_0277;
    const ccomplex_t IT_0279 = s_13*IT_0036;
    const ccomplex_t IT_0280 = s_12*s_23;
    const ccomplex_t IT_0281 = pow(s_23, 2);
    const ccomplex_t IT_0282 = s_13*IT_0166*IT_0281;
    const ccomplex_t IT_0283 = (-2)*IT_0282;
    const ccomplex_t IT_0284 = IT_0280 + IT_0283;
    const ccomplex_t IT_0285 = pow(s_24, 2);
    const ccomplex_t IT_0286 = s_13*IT_0166*IT_0285;
    const ccomplex_t IT_0287 = -IT_0286;
    const ccomplex_t IT_0288 = s_14*s_23*s_24*IT_0166;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = s_13*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0291 = 2*IT_0290;
    const ccomplex_t IT_0292 = IT_0279 + IT_0284 + IT_0287 + IT_0289 + IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0010;
    const ccomplex_t IT_0294 = IT_0022*IT_0293;
    const ccomplex_t IT_0295 = IT_0038*IT_0294;
    const ccomplex_t IT_0296 = IT_0023*IT_0038;
    const ccomplex_t IT_0297 = -IT_0296;
    const ccomplex_t IT_0298 = 2*IT_0297;
    const ccomplex_t IT_0299 = conj(IT_0295)*IT_0298;
    const ccomplex_t IT_0300 = 2*conj(IT_0297);
    const ccomplex_t IT_0301 = IT_0295*IT_0300;
    const ccomplex_t IT_0302 = IT_0032*IT_0293;
    const ccomplex_t IT_0303 = IT_0038*IT_0302;
    const ccomplex_t IT_0304 = IT_0033*IT_0038;
    const ccomplex_t IT_0305 = -IT_0304;
    const ccomplex_t IT_0306 = 2*IT_0305;
    const ccomplex_t IT_0307 = conj(IT_0303)*IT_0306;
    const ccomplex_t IT_0308 = 2*conj(IT_0305);
    const ccomplex_t IT_0309 = IT_0303*IT_0308;
    const ccomplex_t IT_0310 = s_14*s_34*IT_0166;
    const ccomplex_t IT_0311 = s_13*s_23;
    const ccomplex_t IT_0312 = (-2)*IT_0311;
    const ccomplex_t IT_0313 = IT_0197 + IT_0312;
    const ccomplex_t IT_0314 = IT_0166*IT_0313;
    const ccomplex_t IT_0315 = (s_13 + -IT_0310)*(s_12 + (-0.5)*IT_0314);
    const ccomplex_t IT_0316 = (-2)*IT_0315;
    const ccomplex_t IT_0317 = s_24*s_34*IT_0166;
    const ccomplex_t IT_0318 = (s_12 + (-0.5)*IT_0314)*(s_23 + -IT_0317);
    const ccomplex_t IT_0319 = (-2)*IT_0318;
    const ccomplex_t IT_0320 = IT_0010*IT_0022;
    const ccomplex_t IT_0321 = IT_0038*IT_0320;
    const ccomplex_t IT_0322 = IT_0298*conj(IT_0321);
    const ccomplex_t IT_0323 = IT_0300*IT_0321;
    const ccomplex_t IT_0324 = IT_0010*IT_0032;
    const ccomplex_t IT_0325 = IT_0038*IT_0324;
    const ccomplex_t IT_0326 = IT_0306*conj(IT_0325);
    const ccomplex_t IT_0327 = IT_0308*IT_0325;
    const ccomplex_t IT_0328 = s_13*IT_0161;
    const ccomplex_t IT_0329 = s_14*s_34;
    const ccomplex_t IT_0330 = -IT_0329;
    const ccomplex_t IT_0331 = IT_0328 + IT_0330;
    const ccomplex_t IT_0332 = IT_0166*IT_0331;
    const ccomplex_t IT_0333 = -IT_0310 + (-0.25)*IT_0332;
    const ccomplex_t IT_0334 = s_13 + IT_0333;
    const ccomplex_t IT_0335 = IT_0000*IT_0334;
    const ccomplex_t IT_0336 = (-4)*IT_0335;
    const ccomplex_t IT_0337 = s_23*IT_0161;
    const ccomplex_t IT_0338 = s_24*s_34;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = IT_0337 + IT_0339;
    const ccomplex_t IT_0341 = IT_0166*IT_0340;
    const ccomplex_t IT_0342 = -IT_0317 + (-0.25)*IT_0341;
    const ccomplex_t IT_0343 = s_23 + IT_0342;
    const ccomplex_t IT_0344 = IT_0000*IT_0343;
    const ccomplex_t IT_0345 = (-4)*IT_0344;
    const ccomplex_t IT_0346 = -IT_0218 + -IT_0236;
    const ccomplex_t IT_0347 = IT_0298*conj(IT_0346);
    const ccomplex_t IT_0348 = conj(IT_0164)*IT_0298;
    const ccomplex_t IT_0349 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0350 = IT_0009*IT_0349;
    const ccomplex_t IT_0351 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0352 = IT_0009*IT_0351;
    const ccomplex_t IT_0353 = 1.4142135623731*IT_0352;
    const ccomplex_t IT_0354 = (0 + _Complex_I*1)*(IT_0350 + 0.5*IT_0353);
    const ccomplex_t IT_0355 = -IT_0354;
    const ccomplex_t IT_0356 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0357 = IT_0009*IT_0356;
    const ccomplex_t IT_0358 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0359 = IT_0009*IT_0358;
    const ccomplex_t IT_0360 = 1.4142135623731*IT_0359;
    const ccomplex_t IT_0361 = (0 + _Complex_I*1)*(IT_0357 + 0.5*IT_0360);
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = IT_0355*IT_0362;
    const ccomplex_t IT_0364 = cpow((-2)*s_13 + IT_0035 + IT_0161 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0365 = IT_0363*IT_0364;
    const ccomplex_t IT_0366 = (0 + _Complex_I*1)*IT_0365;
    const ccomplex_t IT_0367 = IT_0298*conj(IT_0366);
    const ccomplex_t IT_0368 = IT_0300*IT_0346;
    const ccomplex_t IT_0369 = IT_0164*IT_0300;
    const ccomplex_t IT_0370 = IT_0300*IT_0366;
    const ccomplex_t IT_0371 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0372 = IT_0009*IT_0371;
    const ccomplex_t IT_0373 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0374 = IT_0009*IT_0373;
    const ccomplex_t IT_0375 = 1.4142135623731*IT_0374;
    const ccomplex_t IT_0376 = (0 + _Complex_I*1)*(IT_0372 + (-0.5)*IT_0375);
    const ccomplex_t IT_0377 = -IT_0376;
    const ccomplex_t IT_0378 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0379 = IT_0009*IT_0378;
    const ccomplex_t IT_0380 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0381 = IT_0009*IT_0380;
    const ccomplex_t IT_0382 = 1.4142135623731*IT_0381;
    const ccomplex_t IT_0383 = (0 + _Complex_I*1)*(IT_0379 + (-0.5)*IT_0382);
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = IT_0377*IT_0384;
    const ccomplex_t IT_0386 = IT_0162*IT_0385;
    const ccomplex_t IT_0387 = (0 + _Complex_I*1)*IT_0386;
    const ccomplex_t IT_0388 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0389 = IT_0009*IT_0388;
    const ccomplex_t IT_0390 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0391 = IT_0009*IT_0390;
    const ccomplex_t IT_0392 = 1.4142135623731*IT_0391;
    const ccomplex_t IT_0393 = (0 + _Complex_I*1)*(IT_0389 + (-0.5)*IT_0392);
    const ccomplex_t IT_0394 = -IT_0393;
    const ccomplex_t IT_0395 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0396 = IT_0009*IT_0395;
    const ccomplex_t IT_0397 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0398 = IT_0009*IT_0397;
    const ccomplex_t IT_0399 = 1.4142135623731*IT_0398;
    const ccomplex_t IT_0400 = (0 + _Complex_I*1)*(IT_0396 + (-0.5)*IT_0399);
    const ccomplex_t IT_0401 = -IT_0400;
    const ccomplex_t IT_0402 = IT_0394*IT_0401;
    const ccomplex_t IT_0403 = IT_0364*IT_0402;
    const ccomplex_t IT_0404 = (0 + _Complex_I*1)*IT_0403;
    const ccomplex_t IT_0405 = IT_0387 + IT_0404;
    const ccomplex_t IT_0406 = IT_0306*conj(IT_0405);
    const ccomplex_t IT_0407 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0408 = IT_0009*IT_0407;
    const ccomplex_t IT_0409 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0410 = IT_0009*IT_0409;
    const ccomplex_t IT_0411 = 1.4142135623731*IT_0410;
    const ccomplex_t IT_0412 = (0 + _Complex_I*1)*(IT_0408 + 0.5*IT_0411);
    const ccomplex_t IT_0413 = -IT_0412;
    const ccomplex_t IT_0414 = IT_0257*IT_0413;
    const ccomplex_t IT_0415 = IT_0219*IT_0414;
    const ccomplex_t IT_0416 = (0 + _Complex_I*1)*IT_0415;
    const ccomplex_t IT_0417 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0418 = IT_0009*IT_0417;
    const ccomplex_t IT_0419 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0420 = IT_0009*IT_0419;
    const ccomplex_t IT_0421 = 1.4142135623731*IT_0420;
    const ccomplex_t IT_0422 = (0 + _Complex_I*1)*(IT_0418 + 0.5*IT_0421);
    const ccomplex_t IT_0423 = -IT_0422;
    const ccomplex_t IT_0424 = IT_0247*IT_0423;
    const ccomplex_t IT_0425 = IT_0216*IT_0424;
    const ccomplex_t IT_0426 = (0 + _Complex_I*1)*IT_0425;
    const ccomplex_t IT_0427 = -IT_0416 + -IT_0426;
    const ccomplex_t IT_0428 = IT_0306*conj(IT_0427);
    const ccomplex_t IT_0429 = IT_0308*IT_0405;
    const ccomplex_t IT_0430 = IT_0308*IT_0427;
    const ccomplex_t IT_0431 = 2*IT_0218 + 2*IT_0236;
    const ccomplex_t IT_0432 = IT_0166*IT_0281;
    const ccomplex_t IT_0433 = IT_0166*IT_0285;
    const ccomplex_t IT_0434 = IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0436 = -IT_0036 + -IT_0435;
    const ccomplex_t IT_0437 = IT_0434 + IT_0436;
    const ccomplex_t IT_0438 = m_N_1*IT_0437;
    const ccomplex_t IT_0439 = 2*conj(IT_0093);
    const ccomplex_t IT_0440 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0441 = IT_0439 + IT_0440;
    const ccomplex_t IT_0442 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0443 = s_12 + IT_0182 + IT_0442;
    const ccomplex_t IT_0444 = m_N_4*IT_0443;
    const ccomplex_t IT_0445 = 2*conj(IT_0139);
    const ccomplex_t IT_0446 = (-2)*conj(IT_0040);
    const ccomplex_t IT_0447 = IT_0445 + IT_0446;
    const ccomplex_t IT_0448 = s_14*s_24;
    const ccomplex_t IT_0449 = (-2)*IT_0448;
    const ccomplex_t IT_0450 = IT_0197 + IT_0449;
    const ccomplex_t IT_0451 = IT_0166*IT_0450;
    const ccomplex_t IT_0452 = (s_12 + -IT_0180)*(s_12 + (-0.5)*IT_0451);
    const ccomplex_t IT_0453 = (-2)*IT_0452;
    const ccomplex_t IT_0454 = 2*conj(IT_0165);
    const ccomplex_t IT_0455 = s_13*s_23*IT_0001*IT_0169;
    const ccomplex_t IT_0456 = s_14*s_24*IT_0001*IT_0167;
    const ccomplex_t IT_0457 = s_12*IT_0035;
    const ccomplex_t IT_0458 = s_14*s_24*IT_0035*IT_0166;
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = s_13*s_23*IT_0035*IT_0166;
    const ccomplex_t IT_0461 = -IT_0460;
    const ccomplex_t IT_0462 = s_12*s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0463 = -IT_0462;
    const ccomplex_t IT_0464 = IT_0455 + IT_0456 + IT_0457 + IT_0459 + IT_0461
       + IT_0463;
    const ccomplex_t IT_0465 = s_13*s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0466 = IT_0001*IT_0167*IT_0285;
    const ccomplex_t IT_0467 = IT_0035*IT_0036;
    const ccomplex_t IT_0468 = IT_0036*IT_0166*IT_0167;
    const ccomplex_t IT_0469 = -IT_0468;
    const ccomplex_t IT_0470 = IT_0035*IT_0166*IT_0285;
    const ccomplex_t IT_0471 = -IT_0470;
    const ccomplex_t IT_0472 = s_12*s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0473 = -IT_0472;
    const ccomplex_t IT_0474 = IT_0465 + IT_0466 + IT_0467 + IT_0469 + IT_0471
       + IT_0473;
    const ccomplex_t IT_0475 = IT_0464 + IT_0474;
    const ccomplex_t IT_0476 = 2*IT_0010;
    const ccomplex_t IT_0477 = IT_0032*IT_0476;
    const ccomplex_t IT_0478 = IT_0038*IT_0477;
    const ccomplex_t IT_0479 = (-2)*IT_0366 + -IT_0478;
    const ccomplex_t IT_0480 = IT_0306*conj(IT_0479);
    const ccomplex_t IT_0481 = IT_0308*IT_0479;
    const ccomplex_t IT_0482 = 2*IT_0165;
    const ccomplex_t IT_0483 = conj(IT_0305)*IT_0482;
    const ccomplex_t IT_0484 = IT_0305*IT_0454;
    const ccomplex_t IT_0485 = IT_0001*IT_0169*IT_0281;
    const ccomplex_t IT_0486 = IT_0036*IT_0166*IT_0169;
    const ccomplex_t IT_0487 = -IT_0486;
    const ccomplex_t IT_0488 = IT_0035*IT_0166*IT_0281;
    const ccomplex_t IT_0489 = -IT_0488;
    const ccomplex_t IT_0490 = s_12*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0491 = -IT_0490;
    const ccomplex_t IT_0492 = IT_0465 + IT_0467 + IT_0485 + IT_0487 + IT_0489
       + IT_0491;
    const ccomplex_t IT_0493 = s_14*s_24*IT_0001*IT_0281;
    const ccomplex_t IT_0494 = s_13*s_23*IT_0001*IT_0285;
    const ccomplex_t IT_0495 = s_12*IT_0036;
    const ccomplex_t IT_0496 = s_14*s_24*IT_0036*IT_0166;
    const ccomplex_t IT_0497 = -IT_0496;
    const ccomplex_t IT_0498 = s_13*s_23*IT_0036*IT_0166;
    const ccomplex_t IT_0499 = -IT_0498;
    const ccomplex_t IT_0500 = s_12*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0501 = -IT_0500;
    const ccomplex_t IT_0502 = IT_0493 + IT_0494 + IT_0495 + IT_0497 + IT_0499
       + IT_0501;
    const ccomplex_t IT_0503 = IT_0492 + IT_0502;
    const ccomplex_t IT_0504 = IT_0022*IT_0476;
    const ccomplex_t IT_0505 = IT_0038*IT_0504;
    const ccomplex_t IT_0506 = 2*IT_0416 + 2*IT_0426 + -IT_0505;
    const ccomplex_t IT_0507 = IT_0298*conj(IT_0506);
    const ccomplex_t IT_0508 = IT_0300*IT_0506;
    const ccomplex_t IT_0509 = 2*IT_0431;
    const ccomplex_t IT_0510 = conj(IT_0305)*IT_0509;
    const ccomplex_t IT_0511 = 2*conj(IT_0431);
    const ccomplex_t IT_0512 = IT_0305*IT_0511;
    const ccomplex_t IT_0513 = -IT_0035;
    const ccomplex_t IT_0514 = IT_0170 + IT_0513;
    const ccomplex_t IT_0515 = IT_0000*IT_0514;
    const ccomplex_t IT_0516 = -IT_0036;
    const ccomplex_t IT_0517 = IT_0433 + IT_0516;
    const ccomplex_t IT_0518 = IT_0000*IT_0517;
    const ccomplex_t IT_0519 = conj(IT_0297)*IT_0305;
    const ccomplex_t IT_0520 = 6*IT_0519;
    const ccomplex_t IT_0521 = -IT_0310;
    const ccomplex_t IT_0522 = s_13 + IT_0521;
    const ccomplex_t IT_0523 = IT_0000*IT_0522;
    const ccomplex_t IT_0524 = -IT_0317;
    const ccomplex_t IT_0525 = s_23 + IT_0524;
    const ccomplex_t IT_0526 = IT_0000*IT_0525;
    const ccomplex_t IT_0527 = IT_0523 + IT_0526;
    const ccomplex_t IT_0528 = conj(IT_0297)*IT_0325;
    const ccomplex_t IT_0529 = (-6)*IT_0528;
    const ccomplex_t IT_0530 = (-6)*IT_0321;
    const ccomplex_t IT_0531 = conj(IT_0305)*IT_0530;
    const ccomplex_t IT_0532 = (-6)*conj(IT_0321);
    const ccomplex_t IT_0533 = IT_0305*IT_0532;
    const ccomplex_t IT_0534 = s_12*IT_0005;
    const ccomplex_t IT_0535 = (-4)*conj(IT_0139);
    const ccomplex_t IT_0536 = 0.25*IT_0535;
    const ccomplex_t IT_0537 = conj(IT_0040) + IT_0536;
    const ccomplex_t IT_0538 = (-4)*IT_0139;
    const ccomplex_t IT_0539 = 2*IT_0040;
    const ccomplex_t IT_0540 = (-2)*IT_0139;
    const ccomplex_t IT_0541 = IT_0539 + IT_0540;
    const ccomplex_t IT_0542 = 2*IT_0096;
    const ccomplex_t IT_0543 = (-2)*IT_0093;
    const ccomplex_t IT_0544 = IT_0542 + IT_0543;
    const ccomplex_t IT_0545 = -(IT_0035 + -IT_0168)*(s_12 + (-0.5)*IT_0451);
    const ccomplex_t IT_0546 = (-2)*IT_0545;
    const ccomplex_t IT_0547 = (-2)*IT_0165;
    const ccomplex_t IT_0548 = IT_0159*IT_0377;
    const ccomplex_t IT_0549 = IT_0162*IT_0240*IT_0548;
    const ccomplex_t IT_0550 = IT_0362*IT_0394;
    const ccomplex_t IT_0551 = IT_0250*IT_0364*IT_0550;
    const ccomplex_t IT_0552 = m_N_1*IT_0164;
    const ccomplex_t IT_0553 = IT_0549 + IT_0551 + -IT_0552;
    const ccomplex_t IT_0554 = -IT_0190 + -IT_0192;
    const ccomplex_t IT_0555 = IT_0189 + IT_0554;
    const ccomplex_t IT_0556 = IT_0188*IT_0555;
    const ccomplex_t IT_0557 = (-0.25)*IT_0556;
    const ccomplex_t IT_0558 = s_12 + IT_0199 + IT_0557;
    const ccomplex_t IT_0559 = s_23*s_24*s_34*IT_0166;
    const ccomplex_t IT_0560 = IT_0036*IT_0161;
    const ccomplex_t IT_0561 = 2*IT_0285 + (-2)*IT_0559 + -IT_0560;
    const ccomplex_t IT_0562 = IT_0281 + IT_0561;
    const ccomplex_t IT_0563 = IT_0166*IT_0562;
    const ccomplex_t IT_0564 = -IT_0036 + (-0.5)*IT_0563;
    const ccomplex_t IT_0565 = IT_0432 + IT_0564;
    const ccomplex_t IT_0566 = m_N_1*IT_0565;
    const ccomplex_t IT_0567 = (-2)*IT_0566;
    const ccomplex_t IT_0568 = (-0.5)*IT_0314;
    const ccomplex_t IT_0569 = s_12 + IT_0568;
    const ccomplex_t IT_0570 = pow(m_W, 4);
    const ccomplex_t IT_0571 = s_12*IT_0570;
    const ccomplex_t IT_0572 = s_13*s_23*IT_0161;
    const ccomplex_t IT_0573 = s_14*s_24*IT_0161;
    const ccomplex_t IT_0574 = s_13*s_24*s_34;
    const ccomplex_t IT_0575 = (-2)*IT_0572 + (-2)*IT_0573 + 4*IT_0574;
    const ccomplex_t IT_0576 = IT_0571 + IT_0575;
    const ccomplex_t IT_0577 = IT_0001*IT_0576;
    const ccomplex_t IT_0578 = 0.25*IT_0577;
    const ccomplex_t IT_0579 = (-0.5)*IT_0451;
    const ccomplex_t IT_0580 = IT_0569 + IT_0578 + IT_0579;
    const ccomplex_t IT_0581 = 8*IT_0553;
    const ccomplex_t IT_0582 = s_13*s_24*s_34*IT_0166;
    const ccomplex_t IT_0583 = -IT_0311 + (-2)*IT_0448 + 2*IT_0582;
    const ccomplex_t IT_0584 = IT_0197 + IT_0583;
    const ccomplex_t IT_0585 = IT_0166*IT_0584;
    const ccomplex_t IT_0586 = -IT_0180 + (-0.5)*IT_0585;
    const ccomplex_t IT_0587 = s_12 + IT_0586;
    const ccomplex_t IT_0588 = m_N_1*IT_0587;
    const ccomplex_t IT_0589 = (-2)*IT_0588;
    const ccomplex_t IT_0590 = 1.125*IT_0000;
    const ccomplex_t IT_0591 = s_14*s_23*s_34;
    const ccomplex_t IT_0592 = (-2)*IT_0572 + (-2)*IT_0573 + 4*IT_0591;
    const ccomplex_t IT_0593 = IT_0571 + IT_0592;
    const ccomplex_t IT_0594 = IT_0001*IT_0593;
    const ccomplex_t IT_0595 = 0.25*IT_0594;
    const ccomplex_t IT_0596 = IT_0569 + IT_0579 + IT_0595;
    const ccomplex_t IT_0597 = s_14*s_23*s_34*IT_0166;
    const ccomplex_t IT_0598 = -IT_0311 + (-2)*IT_0448 + 2*IT_0597;
    const ccomplex_t IT_0599 = IT_0197 + IT_0598;
    const ccomplex_t IT_0600 = IT_0166*IT_0599;
    const ccomplex_t IT_0601 = -IT_0180 + (-0.5)*IT_0600;
    const ccomplex_t IT_0602 = s_12 + IT_0601;
    const ccomplex_t IT_0603 = m_N_4*IT_0602;
    const ccomplex_t IT_0604 = (-2)*IT_0603;
    const ccomplex_t IT_0605 = (-2)*IT_0002;
    const ccomplex_t IT_0606 = IT_0570 + IT_0605;
    const ccomplex_t IT_0607 = IT_0001*IT_0606;
    const ccomplex_t IT_0608 = 0.25*IT_0607;
    const ccomplex_t IT_0609 = 1 + IT_0608;
    const ccomplex_t IT_0610 = IT_0000*IT_0609;
    const ccomplex_t IT_0611 = s_13*s_14*s_34*IT_0166;
    const ccomplex_t IT_0612 = IT_0035*IT_0161;
    const ccomplex_t IT_0613 = 2*IT_0169 + (-2)*IT_0611 + -IT_0612;
    const ccomplex_t IT_0614 = IT_0167 + IT_0613;
    const ccomplex_t IT_0615 = IT_0166*IT_0614;
    const ccomplex_t IT_0616 = -IT_0035 + (-0.5)*IT_0615;
    const ccomplex_t IT_0617 = IT_0168 + IT_0616;
    const ccomplex_t IT_0618 = m_N_4*IT_0617;
    const ccomplex_t IT_0619 = (-2)*IT_0618;
    const ccomplex_t IT_0620 = (-2)*conj(IT_0165);
    const ccomplex_t IT_0621 = 8*conj(IT_0553);
    const ccomplex_t IT_0622 = 2*IT_0093;
    const ccomplex_t IT_0623 = (-2)*IT_0096;
    const ccomplex_t IT_0624 = IT_0622 + IT_0623;
    const ccomplex_t IT_0625 = 2*IT_0139;
    const ccomplex_t IT_0626 = (-2)*IT_0040;
    const ccomplex_t IT_0627 = IT_0625 + IT_0626;
    const ccomplex_t IT_0628 = -(IT_0036 + -IT_0432)*(s_12 + (-0.5)*IT_0451);
    const ccomplex_t IT_0629 = (-2)*IT_0628;
    const ccomplex_t IT_0630 = (-2)*IT_0431;
    const ccomplex_t IT_0631 = -IT_0311;
    const ccomplex_t IT_0632 = IT_0197 + IT_0631;
    const ccomplex_t IT_0633 = IT_0166*IT_0632;
    const ccomplex_t IT_0634 = -IT_0180 + (-0.25)*IT_0633;
    const ccomplex_t IT_0635 = s_12 + IT_0634;
    const ccomplex_t IT_0636 = m_N_1*IT_0635;
    const ccomplex_t IT_0637 = (-4)*IT_0636;
    const ccomplex_t IT_0638 = s_12*IT_0002;
    const ccomplex_t IT_0639 = (-2)*IT_0572 + (-2)*IT_0573 + 2*IT_0574 + 2
      *IT_0591 + (-2)*IT_0638;
    const ccomplex_t IT_0640 = IT_0571 + IT_0639;
    const ccomplex_t IT_0641 = IT_0001*IT_0640;
    const ccomplex_t IT_0642 = (-0.25)*IT_0641;
    const ccomplex_t IT_0643 = -IT_0181;
    const ccomplex_t IT_0644 = -IT_0180;
    const ccomplex_t IT_0645 = s_12 + IT_0642 + IT_0643 + IT_0644;
    const ccomplex_t IT_0646 = -IT_0560;
    const ccomplex_t IT_0647 = IT_0281 + IT_0646;
    const ccomplex_t IT_0648 = IT_0166*IT_0647;
    const ccomplex_t IT_0649 = IT_0432 + IT_0516;
    const ccomplex_t IT_0650 = (-4)*IT_0649;
    const ccomplex_t IT_0651 = IT_0648 + IT_0650;
    const ccomplex_t IT_0652 = m_N_1*IT_0651;
    const ccomplex_t IT_0653 = (-2)*conj(IT_0431);
    const ccomplex_t IT_0654 = s_34*m_N_4*IT_0166;
    const ccomplex_t IT_0655 = s_13*s_34*IT_0166;
    const ccomplex_t IT_0656 = -IT_0655;
    const ccomplex_t IT_0657 = s_14 + IT_0656;
    const ccomplex_t IT_0658 = IT_0654*IT_0657;
    const ccomplex_t IT_0659 = s_34*m_N_1*IT_0166;
    const ccomplex_t IT_0660 = s_23*s_34*IT_0166;
    const ccomplex_t IT_0661 = -IT_0660;
    const ccomplex_t IT_0662 = s_24 + IT_0661;
    const ccomplex_t IT_0663 = IT_0659*IT_0662;
    const ccomplex_t IT_0664 = s_34*m_N_1*m_N_4*IT_0166;
    const ccomplex_t IT_0665 = IT_0657*IT_0664;
    const ccomplex_t IT_0666 = IT_0662*IT_0664;
    const ccomplex_t IT_0667 = IT_0337 + IT_0338;
    const ccomplex_t IT_0668 = s_23*IT_0002*IT_0166;
    const ccomplex_t IT_0669 = (-2)*IT_0668;
    const ccomplex_t IT_0670 = IT_0667 + IT_0669;
    const ccomplex_t IT_0671 = IT_0166*IT_0670;
    const ccomplex_t IT_0672 = -IT_0317 + 0.5*IT_0671;
    const ccomplex_t IT_0673 = s_23 + IT_0672;
    const ccomplex_t IT_0674 = m_N_1*IT_0673;
    const ccomplex_t IT_0675 = 2*IT_0674;
    const ccomplex_t IT_0676 = 2*IT_0239;
    const ccomplex_t IT_0677 = IT_0002*IT_0166;
    const ccomplex_t IT_0678 = (s_12 + (-0.5)*IT_0314)*(IT_0161 + -IT_0677);
    const ccomplex_t IT_0679 = (-2)*IT_0678;
    const ccomplex_t IT_0680 = 2*IT_0321;
    const ccomplex_t IT_0681 = IT_0328 + IT_0329;
    const ccomplex_t IT_0682 = s_13*IT_0002*IT_0166;
    const ccomplex_t IT_0683 = (-2)*IT_0682;
    const ccomplex_t IT_0684 = IT_0681 + IT_0683;
    const ccomplex_t IT_0685 = IT_0166*IT_0684;
    const ccomplex_t IT_0686 = -IT_0310 + 0.5*IT_0685;
    const ccomplex_t IT_0687 = s_13 + IT_0686;
    const ccomplex_t IT_0688 = m_N_4*IT_0687;
    const ccomplex_t IT_0689 = 2*IT_0688;
    const ccomplex_t IT_0690 = (-2)*IT_0553;
    const ccomplex_t IT_0691 = m_N_4*IT_0334;
    const ccomplex_t IT_0692 = (-4)*IT_0691;
    const ccomplex_t IT_0693 = (-2)*IT_0260;
    const ccomplex_t IT_0694 = 2*conj(IT_0239);
    const ccomplex_t IT_0695 = (-2)*conj(IT_0553);
    const ccomplex_t IT_0696 = (-2)*conj(IT_0260);
    const ccomplex_t IT_0697 = s_24*s_34*IT_0001*IT_0167;
    const ccomplex_t IT_0698 = s_24*s_34*IT_0035*IT_0166;
    const ccomplex_t IT_0699 = -IT_0698;
    const ccomplex_t IT_0700 = s_12*s_13*IT_0001*IT_0002;
    const ccomplex_t IT_0701 = -IT_0700;
    const ccomplex_t IT_0702 = IT_0267 + IT_0271 + IT_0276 + IT_0697 + IT_0699
       + IT_0701;
    const ccomplex_t IT_0703 = (-2)*IT_0387 + (-2)*IT_0404 + -IT_0505;
    const ccomplex_t IT_0704 = 2*IT_0703;
    const ccomplex_t IT_0705 = s_14*s_34*IT_0001*IT_0281;
    const ccomplex_t IT_0706 = s_14*s_34*IT_0036*IT_0166;
    const ccomplex_t IT_0707 = -IT_0706;
    const ccomplex_t IT_0708 = s_12*s_23*IT_0001*IT_0002;
    const ccomplex_t IT_0709 = -IT_0708;
    const ccomplex_t IT_0710 = IT_0279 + IT_0284 + IT_0290 + IT_0705 + IT_0707
       + IT_0709;
    const ccomplex_t IT_0711 = 2*IT_0506;
    const ccomplex_t IT_0712 = m_N_1*IT_0343;
    const ccomplex_t IT_0713 = (-4)*IT_0712;
    const ccomplex_t IT_0714 = 2*IT_0553;
    const ccomplex_t IT_0715 = 2*IT_0260;
    const ccomplex_t IT_0716 = 2*IT_0325;
    const ccomplex_t IT_0717 = -IT_0677;
    const ccomplex_t IT_0718 = IT_0161 + IT_0717;
    const ccomplex_t IT_0719 = IT_0000*IT_0718;
    const ccomplex_t IT_0720 = (-2)*IT_0239;
    const ccomplex_t IT_0721 = IT_0541*IT_0658 + IT_0624*IT_0663 + IT_0482
      *IT_0702 + IT_0665*IT_0704 + IT_0509*IT_0710 + IT_0666*IT_0711 + IT_0713
      *IT_0714 + IT_0675*IT_0715 + IT_0679*IT_0716 + IT_0530*IT_0719 + IT_0692
      *IT_0720;
    const ccomplex_t IT_0722 = 2*conj(IT_0703);
    const ccomplex_t IT_0723 = 2*conj(IT_0506);
    const ccomplex_t IT_0724 = 2*conj(IT_0553);
    const ccomplex_t IT_0725 = 2*conj(IT_0260);
    const ccomplex_t IT_0726 = (-2)*conj(IT_0239);
    const ccomplex_t IT_0727 = IT_0325*(IT_0178*IT_0658 + IT_0441*IT_0663 +
       IT_0454*IT_0702 + IT_0511*IT_0710 + IT_0532*IT_0719 + IT_0665*IT_0722 +
       IT_0666*IT_0723 + IT_0713*IT_0724 + IT_0675*IT_0725 + IT_0692*IT_0726);
    const ccomplex_t IT_0728 = IT_0355*IT_0401;
    const ccomplex_t IT_0729 = IT_0250*IT_0364*IT_0728;
    const ccomplex_t IT_0730 = m_N_1*IT_0404;
    const ccomplex_t IT_0731 = IT_0152*IT_0384;
    const ccomplex_t IT_0732 = IT_0162*IT_0240*IT_0731;
    const ccomplex_t IT_0733 = m_N_1*IT_0387;
    const ccomplex_t IT_0734 = -IT_0729 + IT_0730 + -IT_0732 + IT_0733;
    const ccomplex_t IT_0735 = m_N_4*IT_0635;
    const ccomplex_t IT_0736 = (-4)*IT_0735;
    const ccomplex_t IT_0737 = -IT_0478;
    const ccomplex_t IT_0738 = 2*conj(IT_0737);
    const ccomplex_t IT_0739 = IT_0511 + IT_0738;
    const ccomplex_t IT_0740 = 8*conj(IT_0040);
    const ccomplex_t IT_0741 = (-8)*conj(IT_0139);
    const ccomplex_t IT_0742 = IT_0740 + IT_0741;
    const ccomplex_t IT_0743 = -IT_0612;
    const ccomplex_t IT_0744 = IT_0167 + IT_0743;
    const ccomplex_t IT_0745 = IT_0166*IT_0744;
    const ccomplex_t IT_0746 = IT_0168 + IT_0513;
    const ccomplex_t IT_0747 = (-4)*IT_0746;
    const ccomplex_t IT_0748 = IT_0745 + IT_0747;
    const ccomplex_t IT_0749 = m_N_4*IT_0748;
    const ccomplex_t IT_0750 = (-2)*conj(IT_0479);
    const ccomplex_t IT_0751 = IT_0620 + IT_0750;
    const ccomplex_t IT_0752 = IT_0285 + IT_0646;
    const ccomplex_t IT_0753 = IT_0166*IT_0752;
    const ccomplex_t IT_0754 = (-4)*IT_0517;
    const ccomplex_t IT_0755 = IT_0753 + IT_0754;
    const ccomplex_t IT_0756 = m_N_1*IT_0755;
    const ccomplex_t IT_0757 = IT_0166*IT_0328;
    const ccomplex_t IT_0758 = s_14*IT_0161;
    const ccomplex_t IT_0759 = s_13*s_34;
    const ccomplex_t IT_0760 = (-2)*IT_0759;
    const ccomplex_t IT_0761 = IT_0758 + IT_0760;
    const ccomplex_t IT_0762 = IT_0188*IT_0761;
    const ccomplex_t IT_0763 = (-0.5)*IT_0757 + (-0.25)*IT_0762;
    const ccomplex_t IT_0764 = s_13 + IT_0763;
    const ccomplex_t IT_0765 = m_N_4*IT_0764;
    const ccomplex_t IT_0766 = (-4)*IT_0765;
    const ccomplex_t IT_0767 = 2*conj(IT_0303);
    const ccomplex_t IT_0768 = 2*conj(IT_0325);
    const ccomplex_t IT_0769 = IT_0166*IT_0337;
    const ccomplex_t IT_0770 = s_24*IT_0161;
    const ccomplex_t IT_0771 = IT_0188*IT_0770;
    const ccomplex_t IT_0772 = (-0.5)*IT_0769 + 0.25*IT_0771;
    const ccomplex_t IT_0773 = s_23 + IT_0772;
    const ccomplex_t IT_0774 = m_N_1*IT_0773;
    const ccomplex_t IT_0775 = 4*IT_0774;
    const ccomplex_t IT_0776 = 2*conj(IT_0295);
    const ccomplex_t IT_0777 = (-8)*conj(IT_0260);
    const ccomplex_t IT_0778 = (-8)*conj(IT_0239);
    const ccomplex_t IT_0779 = -IT_0448;
    const ccomplex_t IT_0780 = IT_0197 + IT_0779;
    const ccomplex_t IT_0781 = IT_0166*IT_0780;
    const ccomplex_t IT_0782 = -IT_0181 + (-0.25)*IT_0781;
    const ccomplex_t IT_0783 = s_12 + IT_0782;
    const ccomplex_t IT_0784 = m_N_1*IT_0783;
    const ccomplex_t IT_0785 = (-4)*IT_0784;
    const ccomplex_t IT_0786 = (-2)*conj(IT_0297);
    const ccomplex_t IT_0787 = (-2)*IT_0311 + -IT_0448 + 2*IT_0582;
    const ccomplex_t IT_0788 = IT_0197 + IT_0787;
    const ccomplex_t IT_0789 = IT_0166*IT_0788;
    const ccomplex_t IT_0790 = -IT_0181 + (-0.5)*IT_0789;
    const ccomplex_t IT_0791 = s_12 + IT_0790;
    const ccomplex_t IT_0792 = m_N_4*IT_0791;
    const ccomplex_t IT_0793 = (-2)*IT_0792;
    const ccomplex_t IT_0794 = (-2)*conj(IT_0305);
    const ccomplex_t IT_0795 = (-2)*conj(IT_0506);
    const ccomplex_t IT_0796 = (-2)*conj(IT_0321);
    const ccomplex_t IT_0797 = 16*conj(IT_0553)*IT_0590 + IT_0589*IT_0722 +
       IT_0736*IT_0739 + IT_0558*IT_0742 + IT_0749*IT_0751 + IT_0300*IT_0756 +
       IT_0766*IT_0767 + IT_0689*IT_0768 + IT_0775*IT_0776 + IT_0610*IT_0777 +
       IT_0645*IT_0778 + IT_0785*IT_0786 + IT_0793*IT_0794 + IT_0567*IT_0795 +
       IT_0713*IT_0796;
    const ccomplex_t IT_0798 = s_23*m_N_1;
    const ccomplex_t IT_0799 = IT_0005*IT_0798;
    const ccomplex_t IT_0800 = s_13*m_N_4;
    const ccomplex_t IT_0801 = IT_0005*IT_0800;
    const ccomplex_t IT_0802 = 2*conj(IT_0479);
    const ccomplex_t IT_0803 = IT_0454 + IT_0802;
    const ccomplex_t IT_0804 = s_12*s_14*s_34*IT_0166;
    const ccomplex_t IT_0805 = IT_0267 + IT_0271 + IT_0273 + IT_0275 + IT_0276
       + IT_0697 + IT_0701 + IT_0804;
    const ccomplex_t IT_0806 = s_12*s_24*s_34*IT_0166;
    const ccomplex_t IT_0807 = IT_0279 + IT_0284 + IT_0287 + IT_0289 + IT_0290
       + IT_0705 + IT_0709 + IT_0806;
    const ccomplex_t IT_0808 = IT_0188*IT_0758;
    const ccomplex_t IT_0809 = (-0.5)*IT_0757 + 0.25*IT_0808;
    const ccomplex_t IT_0810 = s_13 + IT_0809;
    const ccomplex_t IT_0811 = m_N_4*IT_0810;
    const ccomplex_t IT_0812 = 4*IT_0811;
    const ccomplex_t IT_0813 = s_23*s_34;
    const ccomplex_t IT_0814 = (-2)*IT_0813;
    const ccomplex_t IT_0815 = IT_0770 + IT_0814;
    const ccomplex_t IT_0816 = IT_0188*IT_0815;
    const ccomplex_t IT_0817 = (-0.5)*IT_0769 + (-0.25)*IT_0816;
    const ccomplex_t IT_0818 = s_23 + IT_0817;
    const ccomplex_t IT_0819 = m_N_1*IT_0818;
    const ccomplex_t IT_0820 = (-4)*IT_0819;
    const ccomplex_t IT_0821 = s_34*IT_0166;
    const ccomplex_t IT_0822 = IT_0189 + IT_0192;
    const ccomplex_t IT_0823 = s_13*s_23*s_34*IT_0166;
    const ccomplex_t IT_0824 = (-2)*IT_0823;
    const ccomplex_t IT_0825 = IT_0822 + IT_0824;
    const ccomplex_t IT_0826 = IT_0821*IT_0825;
    const ccomplex_t IT_0827 = 0.25*IT_0666*IT_0739 + 0.25*IT_0695*IT_0766 +
       IT_0537*IT_0799 + IT_0102*IT_0801 + 0.25*IT_0665*IT_0803 + 0.25*IT_0722
      *IT_0805 + 0.25*IT_0723*IT_0807 + 0.25*IT_0725*IT_0812 + 0.25*IT_0694
      *IT_0820 + 0.25*IT_0796*IT_0826;
    const ccomplex_t IT_0828 = m_N_1*IT_0366;
    const ccomplex_t IT_0829 = 2*IT_0479;
    const ccomplex_t IT_0830 = IT_0482 + IT_0829;
    const ccomplex_t IT_0831 = 8*IT_0828;
    const ccomplex_t IT_0832 = (-8)*IT_0553;
    const ccomplex_t IT_0833 = IT_0831 + IT_0832;
    const ccomplex_t IT_0834 = 8*IT_0096;
    const ccomplex_t IT_0835 = (-8)*IT_0093;
    const ccomplex_t IT_0836 = IT_0834 + IT_0835;
    const ccomplex_t IT_0837 = (-2)*IT_0737;
    const ccomplex_t IT_0838 = IT_0630 + IT_0837;
    const ccomplex_t IT_0839 = 2*IT_0303;
    const ccomplex_t IT_0840 = 2*IT_0295;
    const ccomplex_t IT_0841 = (-8)*IT_0260;
    const ccomplex_t IT_0842 = (-8)*IT_0239;
    const ccomplex_t IT_0843 = (-2)*IT_0297;
    const ccomplex_t IT_0844 = (-2)*IT_0305;
    const ccomplex_t IT_0845 = (-2)*IT_0703;
    const ccomplex_t IT_0846 = (-2)*IT_0325;
    const ccomplex_t IT_0847 = IT_0680*IT_0689 + (-16)*IT_0590*IT_0734 +
       IT_0711*IT_0736 + IT_0306*IT_0756 + IT_0589*IT_0830 + IT_0580*IT_0833 +
       IT_0558*IT_0836 + IT_0567*IT_0838 + IT_0775*IT_0839 + IT_0766*IT_0840 +
       IT_0645*IT_0841 + IT_0610*IT_0842 + IT_0793*IT_0843 + IT_0785*IT_0844 +
       IT_0749*IT_0845 + IT_0713*IT_0846;
    const ccomplex_t IT_0848 = 8*conj(IT_0096);
    const ccomplex_t IT_0849 = (-8)*conj(IT_0093);
    const ccomplex_t IT_0850 = IT_0848 + IT_0849;
    const ccomplex_t IT_0851 = (-2)*conj(IT_0737);
    const ccomplex_t IT_0852 = IT_0653 + IT_0851;
    const ccomplex_t IT_0853 = 2*conj(IT_0321);
    const ccomplex_t IT_0854 = (-8)*conj(IT_0553);
    const ccomplex_t IT_0855 = (-2)*conj(IT_0703);
    const ccomplex_t IT_0856 = (-2)*conj(IT_0325);
    const ccomplex_t IT_0857 = IT_0723*IT_0736 + IT_0308*IT_0756 + IT_0767
      *IT_0775 + IT_0766*IT_0776 + IT_0645*IT_0777 + IT_0610*IT_0778 + IT_0786
      *IT_0793 + IT_0785*IT_0794 + IT_0589*IT_0803 + IT_0558*IT_0850 + IT_0567
      *IT_0852 + IT_0689*IT_0853 + IT_0580*IT_0854 + IT_0749*IT_0855 + IT_0713
      *IT_0856;
    const ccomplex_t IT_0858 = IT_0416 + IT_0426;
    const ccomplex_t IT_0859 = m_N_4*IT_0858;
    const ccomplex_t IT_0860 = IT_0226*IT_0413;
    const ccomplex_t IT_0861 = IT_0219*IT_0250*IT_0860;
    const ccomplex_t IT_0862 = IT_0207*IT_0423;
    const ccomplex_t IT_0863 = IT_0216*IT_0240*IT_0862;
    const ccomplex_t IT_0864 = -IT_0861 + -IT_0863;
    const ccomplex_t IT_0865 = s_13*IT_0570;
    const ccomplex_t IT_0866 = s_14*s_34*IT_0161;
    const ccomplex_t IT_0867 = (-2)*IT_0866;
    const ccomplex_t IT_0868 = IT_0865 + IT_0867;
    const ccomplex_t IT_0869 = IT_0001*IT_0868;
    const ccomplex_t IT_0870 = (-0.666666666666667)*IT_0310 + (
      -0.166666666666667)*IT_0869;
    const ccomplex_t IT_0871 = s_13 + IT_0870;
    const ccomplex_t IT_0872 = m_N_4*IT_0871;
    const ccomplex_t IT_0873 = (-6)*IT_0872;
    const ccomplex_t IT_0874 = (-2)*IT_0329;
    const ccomplex_t IT_0875 = IT_0328 + IT_0874;
    const ccomplex_t IT_0876 = IT_0166*IT_0875;
    const ccomplex_t IT_0877 = 0.166666666666667*IT_0869 + (-0.666666666666667
      )*IT_0876;
    const ccomplex_t IT_0878 = s_13 + IT_0877;
    const ccomplex_t IT_0879 = m_N_4*IT_0878;
    const ccomplex_t IT_0880 = 6*IT_0879;
    const ccomplex_t IT_0881 = s_23*IT_0570;
    const ccomplex_t IT_0882 = IT_0001*IT_0881;
    const ccomplex_t IT_0883 = (-10)*s_23;
    const ccomplex_t IT_0884 = IT_0882 + IT_0883;
    const ccomplex_t IT_0885 = m_N_1*IT_0884;
    const ccomplex_t IT_0886 = 2*IT_0737;
    const ccomplex_t IT_0887 = IT_0509 + IT_0886;
    const ccomplex_t IT_0888 = (-2)*IT_0479;
    const ccomplex_t IT_0889 = IT_0547 + IT_0888;
    const ccomplex_t IT_0890 = m_N_4*IT_0783;
    const ccomplex_t IT_0891 = (-4)*IT_0890;
    const ccomplex_t IT_0892 = (-2)*IT_0311 + -IT_0448 + 2*IT_0597;
    const ccomplex_t IT_0893 = IT_0197 + IT_0892;
    const ccomplex_t IT_0894 = IT_0166*IT_0893;
    const ccomplex_t IT_0895 = -IT_0181 + (-0.5)*IT_0894;
    const ccomplex_t IT_0896 = s_12 + IT_0895;
    const ccomplex_t IT_0897 = m_N_1*IT_0896;
    const ccomplex_t IT_0898 = (-2)*IT_0897;
    const ccomplex_t IT_0899 = 8*IT_0734;
    const ccomplex_t IT_0900 = 0.5*IT_0285 + -IT_0559 + (-0.5)*IT_0560;
    const ccomplex_t IT_0901 = IT_0281 + IT_0900;
    const ccomplex_t IT_0902 = IT_0166*IT_0901;
    const ccomplex_t IT_0903 = -IT_0433;
    const ccomplex_t IT_0904 = IT_0036 + IT_0902 + IT_0903;
    const ccomplex_t IT_0905 = m_N_1*IT_0904;
    const ccomplex_t IT_0906 = (-4)*IT_0305;
    const ccomplex_t IT_0907 = IT_0169 + IT_0743;
    const ccomplex_t IT_0908 = IT_0166*IT_0907;
    const ccomplex_t IT_0909 = (-4)*IT_0514;
    const ccomplex_t IT_0910 = IT_0908 + IT_0909;
    const ccomplex_t IT_0911 = m_N_4*IT_0910;
    const ccomplex_t IT_0912 = (-2)*IT_0303;
    const ccomplex_t IT_0913 = (-2)*IT_0506;
    const ccomplex_t IT_0914 = (-2)*IT_0295;
    const ccomplex_t IT_0915 = IT_0680*IT_0692 + IT_0619*IT_0704 + IT_0645
      *IT_0833 + IT_0200*IT_0836 + IT_0596*IT_0841 + IT_0675*IT_0846 + IT_0590*
      (16*IT_0239 + (-16)*IT_0864) + 2*IT_0405*IT_0873 + 2*IT_0427*IT_0880 + 2
      *IT_0346*IT_0885 + IT_0652*IT_0887 + IT_0637*IT_0889 + IT_0298*IT_0891 +
       IT_0306*IT_0898 + IT_0610*IT_0899 + IT_0905*IT_0906 + IT_0843*IT_0911 +
       IT_0820*IT_0912 + IT_0604*IT_0913 + IT_0812*IT_0914;
    const ccomplex_t IT_0916 = (-8)*IT_0139;
    const ccomplex_t IT_0917 = 0.125*IT_0916;
    const ccomplex_t IT_0918 = IT_0040 + IT_0917;
    const ccomplex_t IT_0919 = (-4)*IT_0297;
    const ccomplex_t IT_0920 = (-2)*IT_0321;
    const ccomplex_t IT_0921 = 16*IT_0260*IT_0590 + IT_0652*IT_0711 + IT_0692
      *IT_0716 + IT_0619*IT_0830 + IT_0610*IT_0833 + IT_0604*IT_0838 + IT_0596
      *IT_0842 + IT_0637*IT_0845 + IT_0306*IT_0891 + IT_0298*IT_0898 + IT_0645
      *IT_0899 + IT_0844*IT_0911 + IT_0812*IT_0912 + IT_0820*IT_0914 + 8*IT_0200
      *IT_0918 + IT_0905*IT_0919 + IT_0675*IT_0920;
    const ccomplex_t IT_0922 = (-4)*IT_0006;
    const ccomplex_t IT_0923 = -conj(IT_0553) + conj(IT_0828);
    const ccomplex_t IT_0924 = -conj(IT_0040);
    const ccomplex_t IT_0925 = -IT_0093;
    const ccomplex_t IT_0926 = -conj(IT_0093);
    const ccomplex_t IT_0927 = -IT_0139;
    const ccomplex_t IT_0928 = -conj(IT_0139);
    const ccomplex_t IT_0929 = conj(IT_0093)*IT_0239 + IT_0093*conj(IT_0239) +
       IT_0139*(conj(IT_0260) + conj(IT_0553)) + -conj(IT_0096)*(IT_0239 + 
      -IT_0734 + -IT_0864) + -IT_0096*(conj(IT_0239) + -conj(IT_0734) + -conj
      (IT_0864)) + -IT_0040*(conj(IT_0260) + -conj(IT_0859) + -IT_0923) + 
      (IT_0260 + IT_0553)*(conj(IT_0139) + IT_0924) + (conj(IT_0734) + conj
      (IT_0864))*IT_0925 + (IT_0734 + IT_0864)*IT_0926 + (conj(IT_0828) + conj
      (IT_0859))*IT_0927 + (IT_0828 + IT_0859)*(conj(IT_0040) + IT_0928);
    const ccomplex_t IT_0930 = 8*conj(IT_0828);
    const ccomplex_t IT_0931 = IT_0854 + IT_0930;
    const ccomplex_t IT_0932 = 8*conj(IT_0734);
    const ccomplex_t IT_0933 = (-4)*conj(IT_0305);
    const ccomplex_t IT_0934 = (-2)*conj(IT_0303);
    const ccomplex_t IT_0935 = (-2)*conj(IT_0295);
    const ccomplex_t IT_0936 = 2*conj(IT_0405);
    const ccomplex_t IT_0937 = 8*IT_0596;
    const ccomplex_t IT_0938 = IT_0619*IT_0722 + IT_0652*IT_0739 + IT_0637
      *IT_0751 + IT_0604*IT_0795 + IT_0200*IT_0850 + IT_0692*IT_0853 + IT_0675
      *IT_0856 + IT_0590*(16*conj(IT_0239) + (-16)*conj(IT_0864)) + 2*conj
      (IT_0427)*IT_0880 + 2*conj(IT_0346)*IT_0885 + IT_0300*IT_0891 + IT_0308
      *IT_0898 + IT_0786*IT_0911 + IT_0645*IT_0931 + IT_0610*IT_0932 + IT_0905
      *IT_0933 + IT_0820*IT_0934 + IT_0812*IT_0935 + IT_0873*IT_0936 + (0.125
      *IT_0777 + conj(IT_0859))*IT_0937;
    const ccomplex_t IT_0939 = (-4)*conj(IT_0297);
    const ccomplex_t IT_0940 = 16*conj(IT_0260)*IT_0590 + IT_0652*IT_0723 +
       IT_0200*IT_0742 + IT_0692*IT_0768 + IT_0675*IT_0796 + IT_0619*IT_0803 +
       IT_0604*IT_0852 + IT_0637*IT_0855 + IT_0308*IT_0891 + IT_0300*IT_0898 +
       IT_0794*IT_0911 + IT_0610*IT_0931 + IT_0645*IT_0932 + IT_0812*IT_0934 +
       IT_0820*IT_0935 + (0.125*IT_0778 + conj(IT_0864))*IT_0937 + IT_0905
      *IT_0939;
    const ccomplex_t IT_0941 = IT_0000*IT_0746;
    const ccomplex_t IT_0942 = s_12 + IT_0644;
    const ccomplex_t IT_0943 = IT_0000*IT_0942;
    const ccomplex_t IT_0944 = 0.166666666666667*IT_0184;
    const ccomplex_t IT_0945 = 0.166666666666667*IT_0175;
    const ccomplex_t IT_0946 = 0.166666666666667*IT_0546;
    const ccomplex_t IT_0947 = 0.166666666666667*IT_0637;
    const ccomplex_t IT_0948 = 0.166666666666667*IT_0619;
    const ccomplex_t IT_0949 = 0.166666666666667*IT_0749;
    const ccomplex_t IT_0950 = 0.166666666666667*IT_0475;
    const ccomplex_t IT_0951 = 0.166666666666667*IT_0702;
    const ccomplex_t IT_0952 = IT_0165*IT_0941 + -IT_0431*IT_0943 + IT_0541
      *IT_0944 + IT_0544*IT_0945 + IT_0845*IT_0946 + IT_0676*IT_0947 + IT_0693
      *IT_0948 + IT_0714*IT_0949 + IT_0298*IT_0950 + IT_0680*IT_0951;
    const ccomplex_t IT_0953 = 6*IT_0703*(conj(IT_0165)*IT_0941 + -conj
      (IT_0431)*IT_0943 + IT_0178*IT_0944 + IT_0187*IT_0945 + IT_0694*IT_0947 +
       IT_0696*IT_0948 + IT_0724*IT_0949 + IT_0300*IT_0950 + IT_0853*IT_0951);
    const ccomplex_t IT_0954 = 0.166666666666667*IT_0453;
    const ccomplex_t IT_0955 = 0.166666666666667*IT_0589;
    const ccomplex_t IT_0956 = IT_0000*IT_0174;
    const ccomplex_t IT_0957 = 0.166666666666667*IT_0956;
    const ccomplex_t IT_0958 = IT_0000*IT_0183;
    const ccomplex_t IT_0959 = 0.166666666666667*IT_0958;
    const ccomplex_t IT_0960 = 0.166666666666667*IT_0665;
    const ccomplex_t IT_0961 = IT_0703*IT_0941 + -IT_0506*IT_0943 + IT_0544
      *IT_0944 + IT_0541*IT_0945 + IT_0889*IT_0946 + IT_0715*IT_0947 + IT_0720
      *IT_0948 + IT_0716*IT_0951 + IT_0509*IT_0954 + IT_0690*IT_0955 + IT_0298
      *IT_0957 + IT_0843*IT_0959 + IT_0680*IT_0960;
    const ccomplex_t IT_0962 = -conj(IT_0506);
    const ccomplex_t IT_0963 = conj(IT_0703)*IT_0941 + IT_0187*IT_0944 +
       IT_0178*IT_0945 + IT_0620*IT_0946 + IT_0725*IT_0947 + IT_0726*IT_0948 +
       IT_0768*IT_0951 + IT_0511*IT_0954 + IT_0695*IT_0955 + IT_0300*IT_0957 +
       IT_0786*IT_0959 + IT_0853*IT_0960 + IT_0943*IT_0962;
    const ccomplex_t IT_0964 = 0.5*IT_0169 + -IT_0611 + (-0.5)*IT_0612;
    const ccomplex_t IT_0965 = IT_0167 + IT_0964;
    const ccomplex_t IT_0966 = IT_0166*IT_0965;
    const ccomplex_t IT_0967 = -IT_0170;
    const ccomplex_t IT_0968 = IT_0035 + IT_0966 + IT_0967;
    const ccomplex_t IT_0969 = m_N_4*IT_0968;
    const ccomplex_t IT_0970 = -IT_0553;
    const ccomplex_t IT_0971 = IT_0828 + IT_0970;
    const ccomplex_t IT_0972 = 0.25*IT_0793;
    const ccomplex_t IT_0973 = 0.25*IT_0785;
    const ccomplex_t IT_0974 = 0.25*IT_0756;
    const ccomplex_t IT_0975 = 0.25*IT_0766;
    const ccomplex_t IT_0976 = 0.25*IT_0589*IT_0704 + 0.25*IT_0689*IT_0716 +
       0.25*IT_0775*IT_0840 + 0.25*IT_0610*IT_0841 + 0.25*IT_0645*IT_0842 + 0.25
      *IT_0736*IT_0887 + 0.25*IT_0749*IT_0889 + 0.25*IT_0580*IT_0899 + 0.25
      *IT_0567*IT_0913 + 2*IT_0558*IT_0918 + 0.25*IT_0713*IT_0920 + IT_0305
      *IT_0969 + (-4)*IT_0590*IT_0971 + IT_0844*IT_0972 + IT_0843*IT_0973 +
       IT_0298*IT_0974 + IT_0839*IT_0975;
    const ccomplex_t IT_0977 = (1 + 0.5*IT_0003)*(IT_0197 + (-2)*IT_0311);
    const ccomplex_t IT_0978 = (-4)*IT_0295;
    const ccomplex_t IT_0979 = (IT_0040 + 0.25*IT_0538)*IT_0799 + IT_0099
      *IT_0801 + 0.25*IT_0704*IT_0805 + 0.25*IT_0711*IT_0807 + 0.25*IT_0715
      *IT_0812 + 0.25*IT_0676*IT_0820 + 0.25*IT_0665*IT_0830 + 0.25*IT_0666
      *IT_0887 + 0.25*IT_0826*IT_0920 + IT_0690*IT_0975 + 0.25*IT_0977*IT_0978;
    const ccomplex_t IT_0980 = 0.25*IT_0438;
    const ccomplex_t IT_0981 = 0.25*IT_0175;
    const ccomplex_t IT_0982 = m_N_4*IT_0183;
    const ccomplex_t IT_0983 = 0.25*IT_0982;
    const ccomplex_t IT_0984 = m_N_1*IT_0443;
    const ccomplex_t IT_0985 = 0.25*IT_0984;
    const ccomplex_t IT_0986 = 0.25*IT_0956;
    const ccomplex_t IT_0987 = 0.25*IT_0958;
    const ccomplex_t IT_0988 = IT_0000*IT_0443;
    const ccomplex_t IT_0989 = 0.25*IT_0988;
    const ccomplex_t IT_0990 = IT_0000*IT_0437;
    const ccomplex_t IT_0991 = 0.25*IT_0990;
    const ccomplex_t IT_0992 = 0.25*IT_0911;
    const ccomplex_t IT_0993 = 0.25*IT_0891;
    const ccomplex_t IT_0994 = 0.25*IT_0898;
    const ccomplex_t IT_0995 = conj(IT_0260)*IT_0905 + IT_0724*IT_0973 +
       IT_0695*IT_0974 + IT_0441*IT_0980 + IT_0178*IT_0981 + IT_0447*IT_0983 +
       IT_0187*IT_0985 + IT_0722*IT_0986 + IT_0855*IT_0987 + IT_0795*IT_0989 +
       IT_0723*IT_0991 + IT_0694*IT_0992 + IT_0726*IT_0993 + IT_0696*IT_0994;
    const ccomplex_t IT_0996 = (s_12 + -IT_0181)*(s_12 + (-0.5)*IT_0314);
    const ccomplex_t IT_0997 = (-2)*IT_0996;
    const ccomplex_t IT_0998 = s_12 + IT_0643;
    const ccomplex_t IT_0999 = IT_0000*IT_0998;
    const ccomplex_t IT_1000 = (-1.5)*conj(IT_0325)*IT_0527 + conj(IT_0239)
      *IT_0905 + IT_0923*IT_0969 + IT_0724*IT_0972 + IT_0447*IT_0980 + IT_0187
      *IT_0981 + IT_0441*IT_0983 + IT_0178*IT_0985 + IT_0454*IT_0986 + IT_0620
      *IT_0987 + IT_0653*IT_0989 + IT_0511*IT_0991 + IT_0725*IT_0992 + IT_0696
      *IT_0993 + IT_0726*IT_0994 + conj(IT_0297)*IT_0997 + 1.5*conj(IT_0305)*
      (IT_0515 + IT_0518 + (-2)*IT_0999);
    const ccomplex_t IT_1001 = -(IT_0035 + -IT_0170)*(s_12 + (-0.5)*IT_0314);
    const ccomplex_t IT_1002 = (-2)*IT_1001;
    const ccomplex_t IT_1003 = -(s_12 + (-0.5)*IT_0314)*(IT_0036 + -IT_0433);
    const ccomplex_t IT_1004 = (-2)*IT_1003;
    const ccomplex_t IT_1005 = IT_1002 + IT_1004;
    const ccomplex_t IT_1006 = 0.25*IT_1005;
    const ccomplex_t IT_1007 = IT_0239*IT_0905 + IT_0969*IT_0971 + IT_0714
      *IT_0972 + IT_0627*IT_0980 + IT_0544*IT_0981 + IT_0624*IT_0983 + IT_0541
      *IT_0985 + IT_0482*IT_0986 + IT_0547*IT_0987 + IT_0630*IT_0989 + IT_0509
      *IT_0991 + IT_0715*IT_0992 + IT_0693*IT_0993 + IT_0720*IT_0994 + (-3)
      *IT_0305*IT_0999 + IT_0843*IT_1006;
    const ccomplex_t IT_1008 = IT_0260*IT_0905 + IT_0734*IT_0969 + IT_0714
      *IT_0973 + IT_0690*IT_0974 + IT_0624*IT_0980 + IT_0541*IT_0981 + IT_0627
      *IT_0983 + IT_0544*IT_0985 + IT_0704*IT_0986 + IT_0845*IT_0987 + IT_0913
      *IT_0989 + IT_0711*IT_0991 + IT_0676*IT_0992 + IT_0720*IT_0993 + IT_0693
      *IT_0994 + IT_0305*IT_0997 + IT_0844*IT_1006;
    const ccomplex_t IT_1009 = pow(m_W, 6);
    const ccomplex_t IT_1010 = s_12*IT_1009;
    const ccomplex_t IT_1011 = s_14*s_24*IT_0570;
    const ccomplex_t IT_1012 = (-2)*IT_1011;
    const ccomplex_t IT_1013 = IT_1010 + IT_1012;
    const ccomplex_t IT_1014 = IT_0001*IT_1013;
    const ccomplex_t IT_1015 = 0.5*IT_1014;
    const ccomplex_t IT_1016 = (-4)*IT_0311;
    const ccomplex_t IT_1017 = IT_0197 + IT_1015 + IT_1016;
    const ccomplex_t IT_1018 = IT_0166*IT_0593;
    const ccomplex_t IT_1019 = -IT_1018;
    const ccomplex_t IT_1020 = IT_1017 + IT_1019;
    const ccomplex_t IT_1021 = (-2)*IT_0572;
    const ccomplex_t IT_1022 = IT_0571 + IT_1021;
    const ccomplex_t IT_1023 = IT_0166*IT_1022;
    const ccomplex_t IT_1024 = 0.5*IT_1023;
    const ccomplex_t IT_1025 = (-2)*IT_0572 + 2*IT_0574 + 2*IT_0591 + (-2)
      *IT_0638;
    const ccomplex_t IT_1026 = IT_0571 + IT_1025;
    const ccomplex_t IT_1027 = IT_0166*IT_1026;
    const ccomplex_t IT_1028 = (-0.5)*IT_1027;
    const ccomplex_t IT_1029 = s_14*s_23*IT_0161;
    const ccomplex_t IT_1030 = s_13*s_24*IT_0161;
    const ccomplex_t IT_1031 = IT_1029 + IT_1030;
    const ccomplex_t IT_1032 = s_12*s_34*IT_0161;
    const ccomplex_t IT_1033 = -IT_1032;
    const ccomplex_t IT_1034 = IT_1031 + IT_1033;
    const ccomplex_t IT_1035 = IT_0188*IT_1034;
    const ccomplex_t IT_1036 = 0.5*IT_1035;
    const ccomplex_t IT_1037 = IT_0197 + IT_0312 + IT_1024 + IT_1028 + IT_1036;
    const ccomplex_t IT_1038 = -IT_0267;
    const ccomplex_t IT_1039 = IT_0268 + IT_1038;
    const ccomplex_t IT_1040 = s_12*s_13*IT_0161;
    const ccomplex_t IT_1041 = s_23*IT_0035*IT_0161;
    const ccomplex_t IT_1042 = s_23*IT_0169;
    const ccomplex_t IT_1043 = s_12*s_14*s_34;
    const ccomplex_t IT_1044 = -IT_1041 + 2*IT_1042 + (-2)*IT_1043;
    const ccomplex_t IT_1045 = IT_1040 + IT_1044;
    const ccomplex_t IT_1046 = IT_0166*IT_1045;
    const ccomplex_t IT_1047 = (-0.5)*IT_1046;
    const ccomplex_t IT_1048 = IT_1039 + IT_1047;
    const ccomplex_t IT_1049 = 2*conj(IT_0864);
    const ccomplex_t IT_1050 = IT_0726 + IT_1049;
    const ccomplex_t IT_1051 = 0.25*IT_1014;
    const ccomplex_t IT_1052 = IT_0197 + IT_1051;
    const ccomplex_t IT_1053 = -IT_0279;
    const ccomplex_t IT_1054 = IT_0280 + IT_1053;
    const ccomplex_t IT_1055 = s_12*s_23*IT_0161;
    const ccomplex_t IT_1056 = s_13*IT_0036*IT_0161;
    const ccomplex_t IT_1057 = s_14*s_23*s_24;
    const ccomplex_t IT_1058 = s_14*s_34*IT_0036;
    const ccomplex_t IT_1059 = -IT_1056 + (-2)*IT_1057 + 2*IT_1058;
    const ccomplex_t IT_1060 = IT_1055 + IT_1059;
    const ccomplex_t IT_1061 = IT_0166*IT_1060;
    const ccomplex_t IT_1062 = (-0.5)*IT_1061;
    const ccomplex_t IT_1063 = IT_1054 + IT_1062;
    const ccomplex_t IT_1064 = (-2)*IT_0338;
    const ccomplex_t IT_1065 = IT_0337 + IT_1064;
    const ccomplex_t IT_1066 = IT_0166*IT_1065;
    const ccomplex_t IT_1067 = 0.5*IT_0769 + (-0.5)*IT_0771 + (-0.5)*IT_1066;
    const ccomplex_t IT_1068 = s_23 + IT_1067;
    const ccomplex_t IT_1069 = m_N_1*IT_1068;
    const ccomplex_t IT_1070 = 2*IT_1069;
    const ccomplex_t IT_1071 = (-0.25)*IT_1070;
    const ccomplex_t IT_1072 = 0.5*IT_0757 + (-0.5)*IT_0808 + (-0.5)*IT_0876;
    const ccomplex_t IT_1073 = s_13 + IT_1072;
    const ccomplex_t IT_1074 = m_N_4*IT_1073;
    const ccomplex_t IT_1075 = 2*IT_1074;
    const ccomplex_t IT_1076 = (-0.25)*IT_1075;
    const ccomplex_t IT_1077 = (-2)*conj(IT_0734);
    const ccomplex_t IT_1078 = s_24*s_34*IT_0161;
    const ccomplex_t IT_1079 = (-2)*IT_1078;
    const ccomplex_t IT_1080 = IT_0881 + IT_1079;
    const ccomplex_t IT_1081 = IT_0001*IT_1080;
    const ccomplex_t IT_1082 = (-0.666666666666667)*IT_0317 + (
      -0.166666666666667)*IT_1081;
    const ccomplex_t IT_1083 = s_23 + IT_1082;
    const ccomplex_t IT_1084 = m_N_1*IT_1083;
    const ccomplex_t IT_1085 = (-6)*IT_1084;
    const ccomplex_t IT_1086 = (-0.25)*IT_1085;
    const ccomplex_t IT_1087 = (-2)*conj(IT_0828);
    const ccomplex_t IT_1088 = IT_0724 + IT_1087;
    const ccomplex_t IT_1089 = m_N_4*IT_0166;
    const ccomplex_t IT_1090 = IT_0166*IT_0865;
    const ccomplex_t IT_1091 = (-2)*IT_0329 + (-0.5)*IT_1090;
    const ccomplex_t IT_1092 = IT_0328 + IT_1091;
    const ccomplex_t IT_1093 = IT_1089*IT_1092;
    const ccomplex_t IT_1094 = (-2)*IT_1093;
    const ccomplex_t IT_1095 = (-0.25)*IT_1094;
    const ccomplex_t IT_1096 = (-0.25)*IT_0880;
    const ccomplex_t IT_1097 = -IT_0002;
    const ccomplex_t IT_1098 = IT_0570 + IT_1097;
    const ccomplex_t IT_1099 = IT_0166*IT_1098;
    const ccomplex_t IT_1100 = -IT_0161 + 0.25*IT_1099;
    const ccomplex_t IT_1101 = IT_0677 + IT_1100;
    const ccomplex_t IT_1102 = IT_0000*IT_1101;
    const ccomplex_t IT_1103 = 4*IT_1102;
    const ccomplex_t IT_1104 = (-0.25)*IT_1103;
    const ccomplex_t IT_1105 = IT_0768*IT_1104;
    const ccomplex_t IT_1106 = s_34*IT_0161;
    const ccomplex_t IT_1107 = IT_0188*IT_1106;
    const ccomplex_t IT_1108 = IT_0166*IT_0606;
    const ccomplex_t IT_1109 = (-0.333333333333333)*IT_1107 + (
      -0.333333333333333)*IT_1108;
    const ccomplex_t IT_1110 = IT_0161 + IT_1109;
    const ccomplex_t IT_1111 = IT_0000*IT_1110;
    const ccomplex_t IT_1112 = 3*IT_1111;
    const ccomplex_t IT_1113 = (-0.25)*IT_1112;
    const ccomplex_t IT_1114 = (-2)*conj(IT_0346);
    const ccomplex_t IT_1115 = IT_0000*IT_0161;
    const ccomplex_t IT_1116 = 9*IT_1115;
    const ccomplex_t IT_1117 = (-0.25)*IT_1116;
    const ccomplex_t IT_1118 = (-4)*conj(IT_0321);
    const ccomplex_t IT_1119 = s_12*IT_0002*IT_0166;
    const ccomplex_t IT_1120 = -IT_1119;
    const ccomplex_t IT_1121 = IT_0571 + IT_0574 + IT_0591;
    const ccomplex_t IT_1122 = (-2)*IT_0572 + -IT_0638;
    const ccomplex_t IT_1123 = IT_1121 + IT_1122;
    const ccomplex_t IT_1124 = IT_0166*IT_1123;
    const ccomplex_t IT_1125 = (-0.5)*IT_1124;
    const ccomplex_t IT_1126 = -IT_0597;
    const ccomplex_t IT_1127 = IT_0197 + IT_0582 + IT_1120 + IT_1125 + IT_1126;
    const ccomplex_t IT_1128 = (-0.25)*IT_1127;
    const ccomplex_t IT_1129 = s_13*s_14*s_24;
    const ccomplex_t IT_1130 = IT_1042 + IT_1129;
    const ccomplex_t IT_1131 = -IT_1041 + -IT_1043;
    const ccomplex_t IT_1132 = IT_1130 + IT_1131;
    const ccomplex_t IT_1133 = IT_0166*IT_1132;
    const ccomplex_t IT_1134 = (-0.5)*IT_1133;
    const ccomplex_t IT_1135 = -IT_0804;
    const ccomplex_t IT_1136 = IT_0267 + IT_0272 + IT_0275 + IT_1134 + IT_1135;
    const ccomplex_t IT_1137 = (-0.25)*IT_1136;
    const ccomplex_t IT_1138 = (-8)*conj(IT_0297);
    const ccomplex_t IT_1139 = s_13*IT_0285;
    const ccomplex_t IT_1140 = IT_1057 + IT_1139;
    const ccomplex_t IT_1141 = s_12*s_24*s_34;
    const ccomplex_t IT_1142 = -IT_1056 + -IT_1141;
    const ccomplex_t IT_1143 = IT_1140 + IT_1142;
    const ccomplex_t IT_1144 = IT_0166*IT_1143;
    const ccomplex_t IT_1145 = (-0.25)*IT_1144;
    const ccomplex_t IT_1146 = (-0.5)*IT_0806;
    const ccomplex_t IT_1147 = 0.5*IT_0286;
    const ccomplex_t IT_1148 = (-0.5)*IT_0288;
    const ccomplex_t IT_1149 = (-0.5)*IT_0279;
    const ccomplex_t IT_1150 = IT_0280 + IT_1145 + IT_1146 + IT_1147 + IT_1148
       + IT_1149;
    const ccomplex_t IT_1151 = (-0.25)*IT_1150;
    const ccomplex_t IT_1152 = conj(IT_0427)*IT_1020 + -conj(IT_0295)*IT_1037 
      + conj(IT_0703)*IT_1048 + (-0.25)*IT_0885*IT_1050 + IT_0936*IT_1052 +
       IT_0962*IT_1063 + IT_0447*IT_1071 + IT_0187*IT_1076 + IT_1077*IT_1086 +
       IT_1088*IT_1095 + IT_0696*IT_1096 + IT_1105 + IT_0767*IT_1113 + IT_1114
      *IT_1117 + IT_1118*IT_1128 + IT_0939*IT_1137 + IT_1138*IT_1151;
    const ccomplex_t IT_1153 = 0.666666666666667*IT_0769 + (-0.166666666666667
      )*IT_0882;
    const ccomplex_t IT_1154 = s_23 + IT_1153;
    const ccomplex_t IT_1155 = m_N_1*IT_1154;
    const ccomplex_t IT_1156 = (-6)*IT_1155;
    const ccomplex_t IT_1157 = (-4)*conj(IT_0325);
    const ccomplex_t IT_1158 = (-8)*conj(IT_0305);
    const ccomplex_t IT_1159 = conj(IT_0346)*IT_1020 + -conj(IT_0303)*IT_1037 
      + (conj(IT_0165) + conj(IT_0479))*IT_1048 + -(conj(IT_0431) + conj(IT_0737
      ))*IT_1063 + IT_0441*IT_1071 + IT_0178*IT_1076 + IT_1086*IT_1088 + IT_1077
      *IT_1095 + IT_1050*IT_1096 + IT_0853*IT_1104 + IT_0776*IT_1113 + IT_0933
      *IT_1137 + (-0.25)*IT_0696*IT_1156 + IT_1128*IT_1157 + IT_1151*IT_1158;
    const ccomplex_t IT_1160 = IT_0166*IT_0576;
    const ccomplex_t IT_1161 = -IT_1160;
    const ccomplex_t IT_1162 = IT_1017 + IT_1161;
    const ccomplex_t IT_1163 = s_24*s_34*IT_0035;
    const ccomplex_t IT_1164 = -IT_1041 + (-2)*IT_1129 + 2*IT_1163;
    const ccomplex_t IT_1165 = IT_1040 + IT_1164;
    const ccomplex_t IT_1166 = IT_0166*IT_1165;
    const ccomplex_t IT_1167 = (-0.5)*IT_1166;
    const ccomplex_t IT_1168 = IT_1039 + IT_1167;
    const ccomplex_t IT_1169 = -IT_1056 + 2*IT_1139 + (-2)*IT_1141;
    const ccomplex_t IT_1170 = IT_1055 + IT_1169;
    const ccomplex_t IT_1171 = IT_0166*IT_1170;
    const ccomplex_t IT_1172 = (-0.5)*IT_1171;
    const ccomplex_t IT_1173 = IT_1054 + IT_1172;
    const ccomplex_t IT_1174 = IT_0001*IT_0865;
    const ccomplex_t IT_1175 = 0.666666666666667*IT_0757 + (-0.166666666666667
      )*IT_1174;
    const ccomplex_t IT_1176 = s_13 + IT_1175;
    const ccomplex_t IT_1177 = m_N_4*IT_1176;
    const ccomplex_t IT_1178 = (-6)*IT_1177;
    const ccomplex_t IT_1179 = (-0.666666666666667)*IT_1066 +
       0.166666666666667*IT_1081;
    const ccomplex_t IT_1180 = s_23 + IT_1179;
    const ccomplex_t IT_1181 = m_N_1*IT_1180;
    const ccomplex_t IT_1182 = 6*IT_1181;
    const ccomplex_t IT_1183 = (-0.25)*IT_1182;
    const ccomplex_t IT_1184 = (-10)*s_13;
    const ccomplex_t IT_1185 = IT_1174 + IT_1184;
    const ccomplex_t IT_1186 = m_N_4*IT_1185;
    const ccomplex_t IT_1187 = (-0.25)*IT_1186;
    const ccomplex_t IT_1188 = m_N_1*IT_0166;
    const ccomplex_t IT_1189 = IT_0166*IT_0881;
    const ccomplex_t IT_1190 = (-2)*IT_0338 + (-0.5)*IT_1189;
    const ccomplex_t IT_1191 = IT_0337 + IT_1190;
    const ccomplex_t IT_1192 = IT_1188*IT_1191;
    const ccomplex_t IT_1193 = (-2)*IT_1192;
    const ccomplex_t IT_1194 = (-0.25)*IT_1193;
    const ccomplex_t IT_1195 = (-0.25)*IT_0873;
    const ccomplex_t IT_1196 = m_N_1*m_N_4*IT_0166;
    const ccomplex_t IT_1197 = (-0.25)*IT_0570;
    const ccomplex_t IT_1198 = IT_0002 + IT_1197;
    const ccomplex_t IT_1199 = IT_1196*IT_1198;
    const ccomplex_t IT_1200 = (-4)*IT_1199;
    const ccomplex_t IT_1201 = (-0.25)*IT_1200;
    const ccomplex_t IT_1202 = -IT_0582;
    const ccomplex_t IT_1203 = IT_0197 + IT_0597 + IT_1120 + IT_1125 + IT_1202;
    const ccomplex_t IT_1204 = (-0.25)*IT_1203;
    const ccomplex_t IT_1205 = (-0.25)*IT_1133;
    const ccomplex_t IT_1206 = (-0.5)*IT_0804;
    const ccomplex_t IT_1207 = (-0.5)*IT_0272;
    const ccomplex_t IT_1208 = 0.5*IT_0274;
    const ccomplex_t IT_1209 = (-0.5)*IT_0267;
    const ccomplex_t IT_1210 = IT_0268 + IT_1205 + IT_1206 + IT_1207 + IT_1208
       + IT_1209;
    const ccomplex_t IT_1211 = (-0.25)*IT_1210;
    const ccomplex_t IT_1212 = (-0.5)*IT_1144;
    const ccomplex_t IT_1213 = -IT_0806;
    const ccomplex_t IT_1214 = IT_0279 + IT_0287 + IT_0288 + IT_1212 + IT_1213;
    const ccomplex_t IT_1215 = (-0.25)*IT_1214;
    const ccomplex_t IT_1216 = -conj(IT_0295)*IT_1037 + 2*conj(IT_0427)
      *IT_1052 + IT_0447*IT_1071 + IT_0187*IT_1076 + IT_1105 + IT_0767*IT_1113 +
       conj(IT_0405)*IT_1162 + -conj(IT_0703)*IT_1168 + conj(IT_0506)*IT_1173 + 
      (-0.25)*IT_1087*IT_1178 + IT_1077*IT_1183 + IT_0724*IT_1187 + IT_1050
      *IT_1194 + IT_0696*IT_1195 + IT_1114*IT_1201 + IT_1118*IT_1204 + IT_1138
      *IT_1211 + IT_0939*IT_1215;
    const ccomplex_t IT_1217 = (-4)*IT_0164;
    const ccomplex_t IT_1218 = (-0.25)*IT_1217;
    const ccomplex_t IT_1219 = IT_0366 + IT_1218;
    const ccomplex_t IT_1220 = 4*IT_0165;
    const ccomplex_t IT_1221 = 4*IT_0479;
    const ccomplex_t IT_1222 = IT_1220 + IT_1221;
    const ccomplex_t IT_1223 = (-4)*IT_0737;
    const ccomplex_t IT_1224 = (-4)*IT_0431;
    const ccomplex_t IT_1225 = IT_1223 + IT_1224;
    const ccomplex_t IT_1226 = (-8)*IT_0346;
    const ccomplex_t IT_1227 = (-2)*IT_0828;
    const ccomplex_t IT_1228 = IT_0714 + IT_1227;
    const ccomplex_t IT_1229 = (-2)*IT_0734;
    const ccomplex_t IT_1230 = 2*IT_0859;
    const ccomplex_t IT_1231 = IT_0693 + IT_1230;
    const ccomplex_t IT_1232 = 2*IT_0864;
    const ccomplex_t IT_1233 = IT_0720 + IT_1232;
    const ccomplex_t IT_1234 = (-2)*IT_0405;
    const ccomplex_t IT_1235 = (-2)*IT_0427;
    const ccomplex_t IT_1236 = (-4)*IT_0325;
    const ccomplex_t IT_1237 = (-8)*IT_0305;
    const ccomplex_t IT_1238 = IT_0303*IT_1037 + -IT_0624*IT_1071 + -IT_0541
      *IT_1076 + -IT_0680*IT_1104 + -IT_0840*IT_1113 + -IT_0906*IT_1215 + 
      -IT_1162*IT_1219 + 0.25*IT_1168*IT_1222 + 0.25*IT_1173*IT_1225 + 0.25
      *IT_1052*IT_1226 + -IT_1183*IT_1228 + -IT_1187*IT_1229 + -IT_1194*IT_1231 
      + -IT_1195*IT_1233 + -IT_1117*IT_1234 + -IT_1201*IT_1235 + -IT_1204
      *IT_1236 + -IT_1211*IT_1237;
    const ccomplex_t IT_1239 = m_N_1*m_N_4*IT_0161;
    const ccomplex_t IT_1240 = IT_0005*IT_1239;
    const ccomplex_t IT_1241 = (-0.25)*IT_0665;
    const ccomplex_t IT_1242 = (-0.25)*IT_0666;
    const ccomplex_t IT_1243 = (-0.25)*IT_0775;
    const ccomplex_t IT_1244 = (-0.25)*IT_0812;
    const ccomplex_t IT_1245 = (-0.25)*IT_0820;
    const ccomplex_t IT_1246 = (-0.25)*IT_0826;
    const ccomplex_t IT_1247 = (-0.25)*IT_0805;
    const ccomplex_t IT_1248 = (-0.25)*IT_0807;
    const ccomplex_t IT_1249 = -IT_0142*IT_0799 + (IT_0040 + -IT_0139)*IT_0801
       + (-0.25)*IT_0978*IT_1240 + IT_0704*IT_1241 + IT_0711*IT_1242 + IT_0690
      *IT_1243 + IT_0676*IT_1244 + IT_0715*IT_1245 + IT_0846*IT_1246 + IT_0830
      *IT_1247 + IT_0887*IT_1248;
    const ccomplex_t IT_1250 = -IT_0145*IT_0799 + (conj(IT_0040) + -conj
      (IT_0139))*IT_0801 + conj(IT_0303)*IT_0977 + conj(IT_0295)*IT_1240 +
       IT_0722*IT_1241 + IT_0723*IT_1242 + IT_0695*IT_1243 + IT_0694*IT_1244 +
       IT_0725*IT_1245 + IT_0856*IT_1246 + IT_0803*IT_1247 + IT_0739*IT_1248;
    const ccomplex_t IT_1251 = IT_0000*IT_0649;
    const ccomplex_t IT_1252 = (-0.166666666666667)*IT_0438;
    const ccomplex_t IT_1253 = (-0.166666666666667)*IT_0444;
    const ccomplex_t IT_1254 = (-0.166666666666667)*IT_0453;
    const ccomplex_t IT_1255 = (-0.166666666666667)*IT_0629;
    const ccomplex_t IT_1256 = (-0.166666666666667)*IT_0652;
    const ccomplex_t IT_1257 = (-0.166666666666667)*IT_0604;
    const ccomplex_t IT_1258 = (-0.166666666666667)*IT_0736;
    const ccomplex_t IT_1259 = (-0.166666666666667)*IT_0710;
    const ccomplex_t IT_1260 = IT_0165*IT_0943 + -IT_0431*IT_1251 + IT_0627
      *IT_1252 + IT_0624*IT_1253 + IT_0704*IT_1254 + IT_0913*IT_1255 + IT_0720
      *IT_1256 + IT_0715*IT_1257 + IT_0690*IT_1258 + IT_0680*IT_1259;
    const ccomplex_t IT_1261 = (-6)*IT_0506*(conj(IT_0165)*IT_0943 + -conj
      (IT_0431)*IT_1251 + IT_0447*IT_1252 + IT_0441*IT_1253 + IT_0722*IT_1254 +
       IT_0726*IT_1256 + IT_0725*IT_1257 + IT_0695*IT_1258 + IT_0853*IT_1259);
    const ccomplex_t IT_1262 = (-0.166666666666667)*IT_0567;
    const ccomplex_t IT_1263 = (-0.166666666666667)*IT_0988;
    const ccomplex_t IT_1264 = (-0.166666666666667)*IT_0990;
    const ccomplex_t IT_1265 = (-0.166666666666667)*IT_0666;
    const ccomplex_t IT_1266 = (-0.166666666666667)*IT_0503;
    const ccomplex_t IT_1267 = IT_0703*IT_0943 + -IT_0506*IT_1251 + IT_0624
      *IT_1252 + IT_0627*IT_1253 + IT_0830*IT_1254 + IT_0838*IT_1255 + IT_0693
      *IT_1256 + IT_0676*IT_1257 + IT_0716*IT_1259 + IT_0714*IT_1262 + IT_0843
      *IT_1263 + IT_0298*IT_1264 + IT_0680*IT_1265 + IT_0306*IT_1266;
    const ccomplex_t IT_1268 = conj(IT_0703)*IT_0943 + IT_0962*IT_1251 +
       IT_0441*IT_1252 + IT_0447*IT_1253 + IT_0803*IT_1254 + IT_0653*IT_1255 +
       IT_0696*IT_1256 + IT_0694*IT_1257 + IT_0768*IT_1259 + IT_0724*IT_1262 +
       IT_0786*IT_1263 + IT_0300*IT_1264 + IT_0853*IT_1265 + IT_0308*IT_1266;
    const ccomplex_t IT_1269 = 0.25*IT_1070;
    const ccomplex_t IT_1270 = IT_0624*IT_1269;
    const ccomplex_t IT_1271 = 0.25*IT_1075;
    const ccomplex_t IT_1272 = IT_0541*IT_1271;
    const ccomplex_t IT_1273 = 0.25*IT_1085;
    const ccomplex_t IT_1274 = 0.25*IT_1094;
    const ccomplex_t IT_1275 = 0.25*IT_0880;
    const ccomplex_t IT_1276 = 0.25*IT_1103;
    const ccomplex_t IT_1277 = 0.25*IT_1112;
    const ccomplex_t IT_1278 = IT_0840*IT_1277;
    const ccomplex_t IT_1279 = 0.25*IT_1127;
    const ccomplex_t IT_1280 = 0.25*IT_1136;
    const ccomplex_t IT_1281 = 0.25*IT_1150;
    const ccomplex_t IT_1282 = IT_0303*IT_1037;
    const ccomplex_t IT_1283 = (IT_0165 + IT_0479)*IT_1048 + -(IT_0431 +
       IT_0737)*IT_1063 + (-0.25)*IT_0693*IT_1156 + -IT_1270 + -IT_1272 + 
      -IT_1228*IT_1273 + -IT_1229*IT_1274 + -IT_1233*IT_1275 + -IT_0680*IT_1276 
      + -IT_1278 + -IT_1236*IT_1279 + -IT_0906*IT_1280 + -IT_1237*IT_1281 + 
      -IT_1282;
    const ccomplex_t IT_1284 = 0.25*IT_1182;
    const ccomplex_t IT_1285 = 0.25*IT_1193;
    const ccomplex_t IT_1286 = 0.25*IT_0873;
    const ccomplex_t IT_1287 = 0.25*IT_1178;
    const ccomplex_t IT_1288 = 0.25*IT_1116;
    const ccomplex_t IT_1289 = 0.25*IT_1200;
    const ccomplex_t IT_1290 = 0.25*IT_1203;
    const ccomplex_t IT_1291 = 0.25*IT_1052;
    const ccomplex_t IT_1292 = 0.25*IT_1168;
    const ccomplex_t IT_1293 = 0.25*IT_1210;
    const ccomplex_t IT_1294 = 0.25*IT_1173;
    const ccomplex_t IT_1295 = 0.25*IT_1214;
    const ccomplex_t IT_1296 = IT_1162*IT_1219 + -IT_1270 + -IT_1272 + 
      -IT_0680*IT_1276 + -IT_1278 + -IT_1282 + -IT_1228*IT_1284 + -IT_1231
      *IT_1285 + -IT_1233*IT_1286 + -IT_1229*IT_1287 + -IT_1234*IT_1288 + 
      -IT_1235*IT_1289 + -IT_1236*IT_1290 + -IT_1226*IT_1291 + -IT_1222*IT_1292 
      + -IT_1237*IT_1293 + -IT_1225*IT_1294 + -IT_0906*IT_1295;
    const ccomplex_t IT_1297 = IT_0441*IT_1269;
    const ccomplex_t IT_1298 = IT_0178*IT_1271;
    const ccomplex_t IT_1299 = IT_1088*IT_1284;
    const ccomplex_t IT_1300 = 0.25*IT_1186;
    const ccomplex_t IT_1301 = 2*conj(IT_0859);
    const ccomplex_t IT_1302 = IT_0696 + IT_1301;
    const ccomplex_t IT_1303 = IT_1285*IT_1302;
    const ccomplex_t IT_1304 = IT_1050*IT_1286;
    const ccomplex_t IT_1305 = IT_0853*IT_1276;
    const ccomplex_t IT_1306 = IT_0776*IT_1277;
    const ccomplex_t IT_1307 = (-2)*conj(IT_0405);
    const ccomplex_t IT_1308 = IT_1288*IT_1307;
    const ccomplex_t IT_1309 = (-2)*conj(IT_0427);
    const ccomplex_t IT_1310 = IT_1289*IT_1309;
    const ccomplex_t IT_1311 = IT_1157*IT_1290;
    const ccomplex_t IT_1312 = (-8)*conj(IT_0346);
    const ccomplex_t IT_1313 = IT_1291*IT_1312;
    const ccomplex_t IT_1314 = 4*conj(IT_0165);
    const ccomplex_t IT_1315 = 4*conj(IT_0479);
    const ccomplex_t IT_1316 = IT_1314 + IT_1315;
    const ccomplex_t IT_1317 = IT_1292*IT_1316;
    const ccomplex_t IT_1318 = IT_1158*IT_1293;
    const ccomplex_t IT_1319 = (-4)*conj(IT_0737);
    const ccomplex_t IT_1320 = (-4)*conj(IT_0431);
    const ccomplex_t IT_1321 = IT_1319 + IT_1320;
    const ccomplex_t IT_1322 = IT_1294*IT_1321;
    const ccomplex_t IT_1323 = IT_0933*IT_1295;
    const ccomplex_t IT_1324 = conj(IT_0303)*IT_1037;
    const ccomplex_t IT_1325 = IT_0716*IT_1276;
    const ccomplex_t IT_1326 = (-2)*IT_0346;
    const ccomplex_t IT_1327 = (-4)*IT_0321;
    const ccomplex_t IT_1328 = (-8)*IT_0297;
    const ccomplex_t IT_1329 = IT_0295*IT_1037;
    const ccomplex_t IT_1330 = IT_0703*IT_1168 + -IT_0506*IT_1173 + IT_0627
      *IT_1269 + IT_0544*IT_1271 + IT_0839*IT_1277 + IT_1229*IT_1284 + IT_1233
      *IT_1285 + IT_0693*IT_1286 + IT_1227*IT_1287 + IT_0919*IT_1295 + IT_0714
      *IT_1300 + IT_1325 + IT_1289*IT_1326 + IT_1290*IT_1327 + IT_1293*IT_1328 +
       IT_1329;
    const ccomplex_t IT_1331 = IT_0703*IT_1048 + -IT_0506*IT_1063 + (-0.25)
      *IT_0885*IT_1233 + -IT_0627*IT_1269 + -IT_0544*IT_1271 + -IT_1229*IT_1273 
      + -IT_1228*IT_1274 + -IT_0693*IT_1275 + -IT_0839*IT_1277 + -IT_0919
      *IT_1280 + -IT_1325 + -IT_1288*IT_1326 + -IT_1279*IT_1327 + -IT_1281
      *IT_1328 + -IT_1329;
    const ccomplex_t IT_1332 = 4*IT_0006*(conj(IT_0040)*IT_0099 + IT_0040
      *IT_0102 + conj(IT_0139)*IT_0142 + IT_0139*IT_0145) + IT_0165*(IT_0175
      *IT_0178 + IT_0184*IT_0187) + (-8)*IT_0200*((conj(IT_0040) + -conj(IT_0139
      ))*IT_0239 + (IT_0040 + -IT_0139)*conj(IT_0239) + (-0.125)*conj(IT_0260)
      *IT_0263 + (-0.125)*IT_0260*IT_0266) + (IT_0278 + IT_0292)*(IT_0299 +
       IT_0301 + IT_0307 + IT_0309) + (IT_0316 + IT_0319)*(IT_0322 + IT_0323 +
       IT_0326 + IT_0327) + (IT_0336 + IT_0345)*(IT_0347 + IT_0348 + IT_0367 +
       IT_0368 + IT_0369 + IT_0370 + IT_0406 + IT_0428 + IT_0429 + IT_0430) +
       IT_0431*(IT_0438*IT_0441 + IT_0444*IT_0447 + IT_0453*IT_0454) + IT_0475*
      (IT_0480 + IT_0481 + IT_0483 + IT_0484) + IT_0503*(IT_0507 + IT_0508 +
       IT_0510 + IT_0512) + (IT_0515 + IT_0518)*IT_0520 + IT_0527*(IT_0529 +
       IT_0531 + IT_0533) + 4*IT_0534*(IT_0093*conj(IT_0093) + IT_0139*conj
      (IT_0139) + 0.25*conj(IT_0096)*IT_0140 + IT_0096*IT_0145 + IT_0040*IT_0537
       + 0.25*conj(IT_0040)*IT_0538) + conj(IT_0165)*(IT_0175*IT_0541 + IT_0184
      *IT_0544 + IT_0546*IT_0547) + conj(IT_0553)*(IT_0263*IT_0558 + IT_0509
      *IT_0567 + IT_0580*IT_0581 + IT_0547*IT_0589) + conj(IT_0239)*((-16)
      *IT_0260*IT_0590 + 8*IT_0239*IT_0596 + IT_0509*IT_0604 + IT_0581*IT_0610 +
       IT_0547*IT_0619) + IT_0553*(IT_0266*IT_0558 + IT_0511*IT_0567 + IT_0589
      *IT_0620) + IT_0239*(IT_0511*IT_0604 + IT_0619*IT_0620 + IT_0610*IT_0621) 
      + conj(IT_0431)*(IT_0453*IT_0482 + IT_0438*IT_0624 + IT_0444*IT_0627 +
       IT_0629*IT_0630) + conj(IT_0260)*((-16)*IT_0239*IT_0590 + 8*IT_0260
      *IT_0596 + IT_0482*IT_0637 + IT_0581*IT_0645 + IT_0630*IT_0652) + IT_0260*
      (IT_0454*IT_0637 + IT_0621*IT_0645 + IT_0652*IT_0653) + conj(IT_0321)*
      (IT_0544*IT_0658 + IT_0627*IT_0663 + IT_0482*IT_0665 + IT_0509*IT_0666 +
       IT_0675*IT_0676 + IT_0679*IT_0680 + IT_0689*IT_0690 + IT_0692*IT_0693) +
       IT_0321*(IT_0187*IT_0658 + IT_0447*IT_0663 + IT_0454*IT_0665 + IT_0511
      *IT_0666 + IT_0675*IT_0694 + IT_0689*IT_0695 + IT_0692*IT_0696) + conj
      (IT_0325)*IT_0721 + IT_0727 + IT_0734*IT_0797 + 4*IT_0295*IT_0827 + conj
      (IT_0828)*IT_0847 + IT_0828*IT_0857 + conj(IT_0859)*IT_0915 + conj(IT_0864
      )*IT_0921 + IT_0922*IT_0929 + IT_0859*IT_0938 + IT_0864*IT_0940 + 6*conj
      (IT_0703)*IT_0952 + IT_0953 + 6*conj(IT_0479)*IT_0961 + 6*IT_0479*IT_0963 
      + 4*conj(IT_0734)*IT_0976 + 4*conj(IT_0295)*IT_0979 + 4*IT_0305*IT_0995 +
       4*IT_0297*IT_1000 + 4*conj(IT_0297)*IT_1007 + 4*conj(IT_0305)*IT_1008 + (
      -4)*IT_0427*IT_1152 + (-4)*IT_0346*IT_1159 + (-4)*IT_0405*IT_1216 + 4*conj
      (IT_0164)*IT_1238 + (-4)*conj(IT_0303)*IT_1249 + (-4)*IT_0303*IT_1250 + (
      -6)*conj(IT_0506)*IT_1260 + IT_1261 + (-6)*conj(IT_0737)*IT_1267 + (-6)
      *IT_0737*IT_1268 + (-4)*conj(IT_0346)*IT_1283 + (-4)*conj(IT_0366)*IT_1296
       + 4*IT_0366*(IT_1077*IT_1287 + IT_1297 + IT_1298 + IT_1299 + IT_1303 +
       IT_1304 + IT_1305 + IT_1306 + IT_1308 + IT_1310 + IT_1311 + IT_1313 +
       IT_1317 + IT_1318 + IT_1322 + IT_1323 + IT_1324) + 4*IT_0164*(IT_1297 +
       IT_1298 + IT_1299 + IT_1077*IT_1300 + IT_1303 + IT_1304 + IT_1305 +
       IT_1306 + IT_1308 + IT_1310 + IT_1311 + IT_1313 + IT_1317 + IT_1318 +
       IT_1322 + IT_1323 + IT_1324) + 4*conj(IT_0405)*IT_1330 + (-4)*conj
      (IT_0427)*IT_1331;
    return create_ccomplex_return(IT_1332);
}

