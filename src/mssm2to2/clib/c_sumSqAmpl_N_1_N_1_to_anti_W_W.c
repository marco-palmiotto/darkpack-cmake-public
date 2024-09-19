#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_W_W(
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, 2);
    const ccomplex_t IT_0005 = s_12*IT_0004;
    const ccomplex_t IT_0006 = s_13*s_23;
    const ccomplex_t IT_0007 = (-2)*IT_0006;
    const ccomplex_t IT_0008 = IT_0005 + IT_0007;
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = (IT_0000 + -IT_0003)*(s_12 + (-0.5)*IT_0009);
    const ccomplex_t IT_0011 = 2*IT_0010;
    const ccomplex_t IT_0012 = pow(s_24, 2);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = (s_12 + (-0.5)*IT_0009)*(IT_0000 + -IT_0013);
    const ccomplex_t IT_0015 = 2*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*IT_0016*IT_0018;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = cpow(IT_0016, -1);
    const ccomplex_t IT_0022 = IT_0017*IT_0021;
    const ccomplex_t IT_0023 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0018;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = IT_0025*IT_0027;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0024 + IT_0026 + 
      -IT_0028 + -IT_0029);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0020*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = (-2)*IT_0035;
    const ccomplex_t IT_0037 = conj(IT_0035)*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0030;
    const ccomplex_t IT_0039 = IT_0020*IT_0038;
    const ccomplex_t IT_0040 = IT_0034*IT_0039;
    const ccomplex_t IT_0041 = (-2)*IT_0040;
    const ccomplex_t IT_0042 = conj(IT_0040)*IT_0041;
    const ccomplex_t IT_0043 = IT_0037 + IT_0042;
    const ccomplex_t IT_0044 = s_12*s_14*s_34*IT_0001;
    const ccomplex_t IT_0045 = s_12*s_23;
    const ccomplex_t IT_0046 = pow(m_W, -4);
    const ccomplex_t IT_0047 = s_13*s_23*s_24*s_34*IT_0046;
    const ccomplex_t IT_0048 = s_13*IT_0000;
    const ccomplex_t IT_0049 = s_12*s_24*s_34*IT_0001;
    const ccomplex_t IT_0050 = pow(s_13, 2);
    const ccomplex_t IT_0051 = s_24*s_34*IT_0046*IT_0050;
    const ccomplex_t IT_0052 = pow(s_23, 2);
    const ccomplex_t IT_0053 = s_14*s_34*IT_0046*IT_0052;
    const ccomplex_t IT_0054 = s_12*s_13;
    const ccomplex_t IT_0055 = s_13*s_14*s_23*s_34*IT_0046;
    const ccomplex_t IT_0056 = s_23*IT_0000;
    const ccomplex_t IT_0057 = IT_0054 + IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = s_13*s_14*s_24*IT_0001;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = s_23*IT_0001*IT_0002;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = s_23*IT_0001*IT_0050;
    const ccomplex_t IT_0063 = (-2)*IT_0062;
    const ccomplex_t IT_0064 = pow(s_34, 2);
    const ccomplex_t IT_0065 = s_12*s_13*IT_0046*IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = s_13*IT_0001*IT_0012;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0001*IT_0052;
    const ccomplex_t IT_0072 = (-2)*IT_0071;
    const ccomplex_t IT_0073 = s_12*s_23*IT_0046*IT_0064;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = (-2)*IT_0019;
    const ccomplex_t IT_0076 = IT_0038*IT_0075;
    const ccomplex_t IT_0077 = IT_0034*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = 2*IT_0019;
    const ccomplex_t IT_0080 = IT_0038*IT_0079;
    const ccomplex_t IT_0081 = IT_0034*IT_0080;
    const ccomplex_t IT_0082 = 2*conj(IT_0081);
    const ccomplex_t IT_0083 = IT_0078*IT_0082;
    const ccomplex_t IT_0084 = s_12*IT_0000;
    const ccomplex_t IT_0085 = s_14*s_24*IT_0046*IT_0052;
    const ccomplex_t IT_0086 = s_13*s_23*IT_0012*IT_0046;
    const ccomplex_t IT_0087 = s_14*s_24*IT_0000*IT_0001;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = s_13*s_23*IT_0000*IT_0001;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = s_12*s_23*s_24*s_34*IT_0046;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0084 + IT_0085 + IT_0086 + IT_0088 + IT_0090
       + IT_0092;
    const ccomplex_t IT_0094 = IT_0002*IT_0046*IT_0052;
    const ccomplex_t IT_0095 = s_13*s_14*s_23*s_24*IT_0046;
    const ccomplex_t IT_0096 = pow(m_N_1, 4);
    const ccomplex_t IT_0097 = IT_0000*IT_0001*IT_0052;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0000*IT_0001*IT_0002;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = s_12*s_14*s_23*s_34*IT_0046;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = IT_0094 + IT_0095 + IT_0096 + IT_0098 + IT_0100
       + IT_0102;
    const ccomplex_t IT_0104 = IT_0031*IT_0079;
    const ccomplex_t IT_0105 = IT_0034*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0107 = IT_0018*IT_0106;
    const ccomplex_t IT_0108 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0109 = IT_0018*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + (-0.5)*IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0114 = IT_0018*IT_0113;
    const ccomplex_t IT_0115 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0116 = IT_0018*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + (-0.5)*IT_0117);
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = IT_0112*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_23 + IT_0000 + IT_0004 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0125 = IT_0018*IT_0124;
    const ccomplex_t IT_0126 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0127 = IT_0018*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-0.5)*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0132 = IT_0018*IT_0131;
    const ccomplex_t IT_0133 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0134 = IT_0018*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + (-0.5)*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0130*IT_0137;
    const ccomplex_t IT_0139 = cpow((-2)*s_23 + IT_0000 + IT_0004 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0105 + 2*IT_0123 + 2*IT_0141;
    const ccomplex_t IT_0143 = 2*IT_0035;
    const ccomplex_t IT_0144 = conj(IT_0142)*IT_0143;
    const ccomplex_t IT_0145 = 2*conj(IT_0035);
    const ccomplex_t IT_0146 = IT_0142*IT_0145;
    const ccomplex_t IT_0147 = 2*IT_0040;
    const ccomplex_t IT_0148 = conj(IT_0081)*IT_0147;
    const ccomplex_t IT_0149 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0150 = IT_0018*IT_0149;
    const ccomplex_t IT_0151 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0152 = IT_0018*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0150 + 0.5*IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0157 = IT_0018*IT_0156;
    const ccomplex_t IT_0158 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0159 = IT_0018*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + 0.5*IT_0160);
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0155*IT_0162;
    const ccomplex_t IT_0164 = IT_0121*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0167 = IT_0018*IT_0166;
    const ccomplex_t IT_0168 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0169 = IT_0018*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0167 + 0.5*IT_0170);
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0174 = IT_0018*IT_0173;
    const ccomplex_t IT_0175 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0176 = IT_0018*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 0.5*IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0172*IT_0179;
    const ccomplex_t IT_0181 = IT_0139*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = 2*IT_0165 + 2*IT_0182;
    const ccomplex_t IT_0184 = IT_0147*conj(IT_0183);
    const ccomplex_t IT_0185 = 2*conj(IT_0040);
    const ccomplex_t IT_0186 = IT_0081*IT_0185;
    const ccomplex_t IT_0187 = IT_0183*IT_0185;
    const ccomplex_t IT_0188 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0189 = (s_12 + (-0.5)*IT_0009)*(s_13 + -IT_0188);
    const ccomplex_t IT_0190 = (-2)*IT_0189;
    const ccomplex_t IT_0191 = s_24*s_34*IT_0001;
    const ccomplex_t IT_0192 = (s_12 + (-0.5)*IT_0009)*(s_23 + -IT_0191);
    const ccomplex_t IT_0193 = (-2)*IT_0192;
    const ccomplex_t IT_0194 = IT_0019*IT_0031;
    const ccomplex_t IT_0195 = IT_0034*IT_0194;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0143*conj(IT_0196);
    const ccomplex_t IT_0198 = IT_0145*IT_0196;
    const ccomplex_t IT_0199 = IT_0019*IT_0038;
    const ccomplex_t IT_0200 = IT_0034*IT_0199;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0147*conj(IT_0201);
    const ccomplex_t IT_0203 = IT_0185*IT_0201;
    const ccomplex_t IT_0204 = s_13*s_23*IT_0002*IT_0046;
    const ccomplex_t IT_0205 = s_14*s_24*IT_0046*IT_0050;
    const ccomplex_t IT_0206 = s_12*s_13*s_14*s_34*IT_0046;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0084 + IT_0088 + IT_0090 + IT_0204 + IT_0205
       + IT_0207;
    const ccomplex_t IT_0209 = IT_0012*IT_0046*IT_0050;
    const ccomplex_t IT_0210 = IT_0000*IT_0001*IT_0012;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0000*IT_0001*IT_0050;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = s_12*s_13*s_24*s_34*IT_0046;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0095 + IT_0096 + IT_0209 + IT_0211 + IT_0213
       + IT_0215;
    const ccomplex_t IT_0217 = cpow((-2)*s_13 + IT_0000 + IT_0004 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0163*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = cpow((-2)*s_13 + IT_0000 + IT_0004 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0221 = IT_0180*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0105 + (-2)*IT_0219 + (-2)*IT_0222;
    const ccomplex_t IT_0224 = IT_0143*conj(IT_0223);
    const ccomplex_t IT_0225 = IT_0145*IT_0223;
    const ccomplex_t IT_0226 = IT_0138*IT_0220;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = IT_0120*IT_0217;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = (-2)*IT_0227 + (-2)*IT_0229;
    const ccomplex_t IT_0231 = IT_0147*conj(IT_0230);
    const ccomplex_t IT_0232 = IT_0185*IT_0230;
    const ccomplex_t IT_0233 = 2*IT_0055;
    const ccomplex_t IT_0234 = IT_0054 + IT_0056 + IT_0059 + IT_0061 + IT_0063
       + IT_0233;
    const ccomplex_t IT_0235 = 2*IT_0047;
    const ccomplex_t IT_0236 = IT_0045 + IT_0048 + IT_0068 + IT_0070 + IT_0072
       + IT_0235;
    const ccomplex_t IT_0237 = IT_0031*IT_0075;
    const ccomplex_t IT_0238 = IT_0034*IT_0237;
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = IT_0143*conj(IT_0239);
    const ccomplex_t IT_0241 = IT_0145*IT_0239;
    const ccomplex_t IT_0242 = conj(IT_0078)*IT_0147;
    const ccomplex_t IT_0243 = IT_0078*IT_0185;
    const ccomplex_t IT_0244 = s_13*IT_0004;
    const ccomplex_t IT_0245 = s_14*s_34;
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = IT_0244 + IT_0246;
    const ccomplex_t IT_0248 = IT_0001*IT_0247;
    const ccomplex_t IT_0249 = -IT_0188 + (-0.25)*IT_0248;
    const ccomplex_t IT_0250 = s_13 + IT_0249;
    const ccomplex_t IT_0251 = IT_0000*IT_0250;
    const ccomplex_t IT_0252 = (-4)*IT_0251;
    const ccomplex_t IT_0253 = s_23*IT_0004;
    const ccomplex_t IT_0254 = s_24*s_34;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = IT_0253 + IT_0255;
    const ccomplex_t IT_0257 = IT_0001*IT_0256;
    const ccomplex_t IT_0258 = -IT_0191 + (-0.25)*IT_0257;
    const ccomplex_t IT_0259 = s_23 + IT_0258;
    const ccomplex_t IT_0260 = IT_0000*IT_0259;
    const ccomplex_t IT_0261 = (-4)*IT_0260;
    const ccomplex_t IT_0262 = IT_0227 + IT_0229;
    const ccomplex_t IT_0263 = IT_0143*conj(IT_0262);
    const ccomplex_t IT_0264 = -IT_0165 + -IT_0182;
    const ccomplex_t IT_0265 = IT_0143*conj(IT_0264);
    const ccomplex_t IT_0266 = IT_0145*IT_0262;
    const ccomplex_t IT_0267 = IT_0145*IT_0264;
    const ccomplex_t IT_0268 = -IT_0123 + -IT_0141;
    const ccomplex_t IT_0269 = IT_0147*conj(IT_0268);
    const ccomplex_t IT_0270 = IT_0147*conj(IT_0222);
    const ccomplex_t IT_0271 = IT_0147*conj(IT_0219);
    const ccomplex_t IT_0272 = IT_0185*IT_0268;
    const ccomplex_t IT_0273 = IT_0185*IT_0222;
    const ccomplex_t IT_0274 = IT_0185*IT_0219;
    const ccomplex_t IT_0275 = IT_0046*IT_0064;
    const ccomplex_t IT_0276 = 0.5*IT_0275;
    const ccomplex_t IT_0277 = 1 + IT_0276;
    const ccomplex_t IT_0278 = s_12*IT_0277;
    const ccomplex_t IT_0279 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0280 = cos(alpha);
    const ccomplex_t IT_0281 = IT_0021*IT_0280;
    const ccomplex_t IT_0282 = IT_0279*IT_0281;
    const ccomplex_t IT_0283 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0284 = sin(alpha);
    const ccomplex_t IT_0285 = IT_0021*IT_0284;
    const ccomplex_t IT_0286 = IT_0283*IT_0285;
    const ccomplex_t IT_0287 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0288 = IT_0018*IT_0280;
    const ccomplex_t IT_0289 = IT_0287*IT_0288;
    const ccomplex_t IT_0290 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0291 = IT_0018*IT_0284;
    const ccomplex_t IT_0292 = IT_0290*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*(IT_0282 + -IT_0286 + 
      -IT_0289 + IT_0292);
    const ccomplex_t IT_0294 = cos(beta);
    const ccomplex_t IT_0295 = sin(beta);
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*m_W*e_em*IT_0018*(IT_0280
      *IT_0294 + IT_0284*IT_0295);
    const ccomplex_t IT_0297 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0298 = IT_0296*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*IT_0298;
    const ccomplex_t IT_0300 = IT_0293*IT_0299;
    const ccomplex_t IT_0301 = m_N_1*IT_0032;
    const ccomplex_t IT_0302 = m_N_1*IT_0039;
    const ccomplex_t IT_0303 = -IT_0302;
    const ccomplex_t IT_0304 = IT_0301 + IT_0303;
    const ccomplex_t IT_0305 = IT_0034*IT_0304;
    const ccomplex_t IT_0306 = IT_0281*IT_0283;
    const ccomplex_t IT_0307 = IT_0279*IT_0285;
    const ccomplex_t IT_0308 = IT_0288*IT_0290;
    const ccomplex_t IT_0309 = IT_0287*IT_0291;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*(IT_0306 + IT_0307 + 
      -IT_0308 + -IT_0309);
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*-1)*m_W*e_em*IT_0018*(IT_0284
      *IT_0294 + -IT_0280*IT_0295);
    const ccomplex_t IT_0313 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0314 = IT_0312*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0311*IT_0315;
    const ccomplex_t IT_0317 = -IT_0300 + -IT_0305 + -IT_0316;
    const ccomplex_t IT_0318 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0319 = IT_0281*IT_0318;
    const ccomplex_t IT_0320 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0321 = IT_0285*IT_0320;
    const ccomplex_t IT_0322 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0323 = IT_0288*IT_0322;
    const ccomplex_t IT_0324 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0325 = IT_0291*IT_0324;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*(IT_0319 + IT_0321 + 
      -IT_0323 + -IT_0325);
    const ccomplex_t IT_0327 = -IT_0326;
    const ccomplex_t IT_0328 = IT_0315*IT_0327;
    const ccomplex_t IT_0329 = IT_0281*IT_0320;
    const ccomplex_t IT_0330 = IT_0285*IT_0318;
    const ccomplex_t IT_0331 = IT_0288*IT_0324;
    const ccomplex_t IT_0332 = IT_0291*IT_0322;
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*(IT_0329 + -IT_0330 + 
      -IT_0331 + IT_0332);
    const ccomplex_t IT_0334 = IT_0299*IT_0333;
    const ccomplex_t IT_0335 = IT_0305 + -IT_0328 + -IT_0334;
    const ccomplex_t IT_0336 = 4*IT_0335;
    const ccomplex_t IT_0337 = -IT_0003;
    const ccomplex_t IT_0338 = IT_0000 + IT_0337;
    const ccomplex_t IT_0339 = IT_0000*IT_0338;
    const ccomplex_t IT_0340 = -IT_0339;
    const ccomplex_t IT_0341 = -IT_0013;
    const ccomplex_t IT_0342 = IT_0000 + IT_0341;
    const ccomplex_t IT_0343 = IT_0000*IT_0342;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = conj(IT_0035)*IT_0040;
    const ccomplex_t IT_0346 = 6*IT_0345;
    const ccomplex_t IT_0347 = -IT_0188;
    const ccomplex_t IT_0348 = s_13 + IT_0347;
    const ccomplex_t IT_0349 = IT_0000*IT_0348;
    const ccomplex_t IT_0350 = -IT_0191;
    const ccomplex_t IT_0351 = s_23 + IT_0350;
    const ccomplex_t IT_0352 = IT_0000*IT_0351;
    const ccomplex_t IT_0353 = IT_0349 + IT_0352;
    const ccomplex_t IT_0354 = conj(IT_0035)*IT_0201;
    const ccomplex_t IT_0355 = (-6)*IT_0354;
    const ccomplex_t IT_0356 = (-6)*IT_0196;
    const ccomplex_t IT_0357 = conj(IT_0040)*IT_0356;
    const ccomplex_t IT_0358 = (-6)*conj(IT_0196);
    const ccomplex_t IT_0359 = IT_0040*IT_0358;
    const ccomplex_t IT_0360 = IT_0000*IT_0277;
    const ccomplex_t IT_0361 = (-4)*IT_0317;
    const ccomplex_t IT_0362 = (-4)*conj(IT_0317);
    const ccomplex_t IT_0363 = s_34*m_N_1*IT_0001;
    const ccomplex_t IT_0364 = s_23*s_34*IT_0001;
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = s_24 + IT_0365;
    const ccomplex_t IT_0367 = IT_0363*IT_0366;
    const ccomplex_t IT_0368 = 2*IT_0317;
    const ccomplex_t IT_0369 = IT_0001*IT_0064;
    const ccomplex_t IT_0370 = (s_12 + (-0.5)*IT_0009)*(IT_0004 + -IT_0369);
    const ccomplex_t IT_0371 = (-2)*IT_0370;
    const ccomplex_t IT_0372 = 2*IT_0196;
    const ccomplex_t IT_0373 = s_13*s_34*IT_0001;
    const ccomplex_t IT_0374 = -IT_0373;
    const ccomplex_t IT_0375 = s_14 + IT_0374;
    const ccomplex_t IT_0376 = IT_0363*IT_0375;
    const ccomplex_t IT_0377 = (-2)*IT_0335;
    const ccomplex_t IT_0378 = 2*conj(IT_0317);
    const ccomplex_t IT_0379 = (-2)*conj(IT_0335);
    const ccomplex_t IT_0380 = 2*IT_0335;
    const ccomplex_t IT_0381 = s_34*IT_0000*IT_0001;
    const ccomplex_t IT_0382 = IT_0366*IT_0381;
    const ccomplex_t IT_0383 = 2*IT_0142;
    const ccomplex_t IT_0384 = 2*IT_0201;
    const ccomplex_t IT_0385 = IT_0375*IT_0381;
    const ccomplex_t IT_0386 = 2*IT_0223;
    const ccomplex_t IT_0387 = -IT_0369;
    const ccomplex_t IT_0388 = IT_0004 + IT_0387;
    const ccomplex_t IT_0389 = IT_0000*IT_0388;
    const ccomplex_t IT_0390 = (-2)*IT_0317;
    const ccomplex_t IT_0391 = m_N_1*IT_0250;
    const ccomplex_t IT_0392 = (-4)*IT_0391;
    const ccomplex_t IT_0393 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0394 = IT_0119*IT_0162;
    const ccomplex_t IT_0395 = IT_0121*IT_0393*IT_0394;
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0397 = IT_0137*IT_0179;
    const ccomplex_t IT_0398 = IT_0139*IT_0396*IT_0397;
    const ccomplex_t IT_0399 = IT_0395 + IT_0398;
    const ccomplex_t IT_0400 = (-2)*IT_0399;
    const ccomplex_t IT_0401 = IT_0244 + IT_0245;
    const ccomplex_t IT_0402 = s_13*IT_0001*IT_0064;
    const ccomplex_t IT_0403 = (-2)*IT_0402;
    const ccomplex_t IT_0404 = IT_0401 + IT_0403;
    const ccomplex_t IT_0405 = IT_0001*IT_0404;
    const ccomplex_t IT_0406 = -IT_0188 + 0.5*IT_0405;
    const ccomplex_t IT_0407 = s_13 + IT_0406;
    const ccomplex_t IT_0408 = m_N_1*IT_0407;
    const ccomplex_t IT_0409 = 2*IT_0408;
    const ccomplex_t IT_0410 = IT_0219 + IT_0222;
    const ccomplex_t IT_0411 = m_N_1*IT_0410;
    const ccomplex_t IT_0412 = -IT_0411;
    const ccomplex_t IT_0413 = (-2)*IT_0412;
    const ccomplex_t IT_0414 = IT_0035*IT_0353 + (-0.166666666666667)*IT_0367
      *IT_0380 + (-0.166666666666667)*IT_0382*IT_0383 + (-0.166666666666667)
      *IT_0371*IT_0384 + (-0.166666666666667)*IT_0385*IT_0386 + (
      -0.166666666666667)*IT_0356*IT_0389 + (-0.166666666666667)*IT_0376*IT_0390
       + (-0.166666666666667)*IT_0392*IT_0400 + (-0.166666666666667)*IT_0409
      *IT_0413;
    const ccomplex_t IT_0415 = 2*conj(IT_0335);
    const ccomplex_t IT_0416 = 2*conj(IT_0142);
    const ccomplex_t IT_0417 = 2*conj(IT_0223);
    const ccomplex_t IT_0418 = (-2)*conj(IT_0317);
    const ccomplex_t IT_0419 = (-2)*conj(IT_0399);
    const ccomplex_t IT_0420 = (-2)*conj(IT_0412);
    const ccomplex_t IT_0421 = -IT_0056;
    const ccomplex_t IT_0422 = IT_0054 + IT_0421;
    const ccomplex_t IT_0423 = s_23*IT_0002;
    const ccomplex_t IT_0424 = s_12*s_14*s_34;
    const ccomplex_t IT_0425 = s_23*IT_0000*IT_0004;
    const ccomplex_t IT_0426 = s_12*s_13*IT_0004;
    const ccomplex_t IT_0427 = -IT_0424 + (-0.5)*IT_0425 + 0.5*IT_0426;
    const ccomplex_t IT_0428 = IT_0423 + IT_0427;
    const ccomplex_t IT_0429 = IT_0001*IT_0428;
    const ccomplex_t IT_0430 = 2*IT_0429;
    const ccomplex_t IT_0431 = (-0.5)*IT_0430;
    const ccomplex_t IT_0432 = IT_0422 + IT_0431;
    const ccomplex_t IT_0433 = -IT_0048;
    const ccomplex_t IT_0434 = IT_0045 + IT_0433;
    const ccomplex_t IT_0435 = s_13*IT_0012;
    const ccomplex_t IT_0436 = s_12*s_24*s_34;
    const ccomplex_t IT_0437 = s_13*IT_0000*IT_0004;
    const ccomplex_t IT_0438 = s_12*s_23*IT_0004;
    const ccomplex_t IT_0439 = -IT_0436 + (-0.5)*IT_0437 + 0.5*IT_0438;
    const ccomplex_t IT_0440 = IT_0435 + IT_0439;
    const ccomplex_t IT_0441 = IT_0001*IT_0440;
    const ccomplex_t IT_0442 = 2*IT_0441;
    const ccomplex_t IT_0443 = (-0.5)*IT_0442;
    const ccomplex_t IT_0444 = IT_0434 + IT_0443;
    const ccomplex_t IT_0445 = IT_0382 + IT_0385;
    const ccomplex_t IT_0446 = 2*conj(IT_0196);
    const ccomplex_t IT_0447 = 2*conj(IT_0239);
    const ccomplex_t IT_0448 = IT_0446 + IT_0447;
    const ccomplex_t IT_0449 = s_14*s_23*s_34*IT_0046;
    const ccomplex_t IT_0450 = s_13*s_23*IT_0001;
    const ccomplex_t IT_0451 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0452 = -IT_0450 + -IT_0451;
    const ccomplex_t IT_0453 = s_12 + IT_0449 + IT_0452;
    const ccomplex_t IT_0454 = m_N_1*IT_0453;
    const ccomplex_t IT_0455 = s_23*s_24*s_34*IT_0046;
    const ccomplex_t IT_0456 = IT_0000 + IT_0455;
    const ccomplex_t IT_0457 = IT_0001*IT_0052;
    const ccomplex_t IT_0458 = -IT_0013 + -IT_0457;
    const ccomplex_t IT_0459 = IT_0456 + IT_0458;
    const ccomplex_t IT_0460 = m_N_1*IT_0459;
    const ccomplex_t IT_0461 = -IT_0460;
    const ccomplex_t IT_0462 = s_14*s_24;
    const ccomplex_t IT_0463 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0464 = -IT_0006 + (-2)*IT_0462 + 2*IT_0463;
    const ccomplex_t IT_0465 = IT_0005 + IT_0464;
    const ccomplex_t IT_0466 = IT_0001*IT_0465;
    const ccomplex_t IT_0467 = -IT_0450 + (-0.5)*IT_0466;
    const ccomplex_t IT_0468 = s_12 + IT_0467;
    const ccomplex_t IT_0469 = m_N_1*IT_0468;
    const ccomplex_t IT_0470 = (-2)*IT_0469;
    const ccomplex_t IT_0471 = 2*conj(IT_0399);
    const ccomplex_t IT_0472 = IT_0000*IT_0004;
    const ccomplex_t IT_0473 = -IT_0472;
    const ccomplex_t IT_0474 = IT_0050 + IT_0473;
    const ccomplex_t IT_0475 = IT_0001*IT_0474;
    const ccomplex_t IT_0476 = IT_0001*IT_0050;
    const ccomplex_t IT_0477 = -IT_0000;
    const ccomplex_t IT_0478 = IT_0476 + IT_0477;
    const ccomplex_t IT_0479 = (-4)*IT_0478;
    const ccomplex_t IT_0480 = IT_0475 + IT_0479;
    const ccomplex_t IT_0481 = m_N_1*IT_0480;
    const ccomplex_t IT_0482 = 2*conj(IT_0412);
    const ccomplex_t IT_0483 = s_13*s_14*s_34*IT_0046;
    const ccomplex_t IT_0484 = IT_0000 + IT_0483;
    const ccomplex_t IT_0485 = -IT_0003 + -IT_0476;
    const ccomplex_t IT_0486 = IT_0484 + IT_0485;
    const ccomplex_t IT_0487 = IT_0000*IT_0486;
    const ccomplex_t IT_0488 = -IT_0487;
    const ccomplex_t IT_0489 = IT_0000*IT_0459;
    const ccomplex_t IT_0490 = -IT_0489;
    const ccomplex_t IT_0491 = IT_0488 + IT_0490;
    const ccomplex_t IT_0492 = s_24*s_34*IT_0000*IT_0001;
    const ccomplex_t IT_0493 = -IT_0492;
    const ccomplex_t IT_0494 = s_14*s_34*IT_0000*IT_0001;
    const ccomplex_t IT_0495 = -IT_0494;
    const ccomplex_t IT_0496 = 2*conj(IT_0201);
    const ccomplex_t IT_0497 = -IT_0457;
    const ccomplex_t IT_0498 = IT_0000 + IT_0497;
    const ccomplex_t IT_0499 = IT_0000*IT_0498;
    const ccomplex_t IT_0500 = -IT_0499;
    const ccomplex_t IT_0501 = 6*conj(IT_0142);
    const ccomplex_t IT_0502 = -IT_0476;
    const ccomplex_t IT_0503 = IT_0000 + IT_0502;
    const ccomplex_t IT_0504 = IT_0000*IT_0503;
    const ccomplex_t IT_0505 = -IT_0504;
    const ccomplex_t IT_0506 = 6*conj(IT_0223);
    const ccomplex_t IT_0507 = -IT_0450;
    const ccomplex_t IT_0508 = s_12 + IT_0507;
    const ccomplex_t IT_0509 = IT_0000*IT_0508;
    const ccomplex_t IT_0510 = (-6)*conj(IT_0142);
    const ccomplex_t IT_0511 = (-6)*conj(IT_0223);
    const ccomplex_t IT_0512 = m_N_1*IT_0486;
    const ccomplex_t IT_0513 = -IT_0512;
    const ccomplex_t IT_0514 = s_13*s_24*s_34*IT_0046;
    const ccomplex_t IT_0515 = s_12 + IT_0452 + IT_0514;
    const ccomplex_t IT_0516 = m_N_1*IT_0515;
    const ccomplex_t IT_0517 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0518 = 2*IT_0002 + -IT_0472 + (-2)*IT_0517;
    const ccomplex_t IT_0519 = IT_0050 + IT_0518;
    const ccomplex_t IT_0520 = IT_0001*IT_0519;
    const ccomplex_t IT_0521 = -IT_0476 + 0.5*IT_0520;
    const ccomplex_t IT_0522 = IT_0000 + IT_0521;
    const ccomplex_t IT_0523 = m_N_1*IT_0522;
    const ccomplex_t IT_0524 = 2*IT_0523;
    const ccomplex_t IT_0525 = IT_0000*IT_0453;
    const ccomplex_t IT_0526 = IT_0000*IT_0515;
    const ccomplex_t IT_0527 = IT_0525 + IT_0526;
    const ccomplex_t IT_0528 = (-2)*conj(IT_0035);
    const ccomplex_t IT_0529 = conj(IT_0264)*IT_0432 + conj(IT_0262)*IT_0444 +
       (-0.25)*IT_0445*IT_0448 + (-0.25)*IT_0378*IT_0454 + (-0.25)*IT_0415
      *IT_0461 + (-0.25)*IT_0470*IT_0471 + (-0.25)*IT_0481*IT_0482 + (-0.25)
      *IT_0145*IT_0491 + (-0.25)*(IT_0045 + IT_0047 + IT_0048 + IT_0051 +
       IT_0053 + IT_0057 + IT_0063 + IT_0066 + IT_0072 + IT_0074 + IT_0493 +
       IT_0495)*IT_0496 + (-0.25)*IT_0500*IT_0501 + (-0.25)*IT_0505*IT_0506 + (
      -0.25)*IT_0509*(IT_0510 + IT_0511) + (-0.25)*IT_0418*IT_0513 + (-0.25)
      *IT_0379*IT_0516 + (-0.25)*IT_0419*IT_0524 + (-0.25)*IT_0527*IT_0528;
    const ccomplex_t IT_0530 = IT_0045 + IT_0047 + IT_0048;
    const ccomplex_t IT_0531 = IT_0053 + IT_0072 + IT_0074 + IT_0495 + IT_0530;
    const ccomplex_t IT_0532 = (-2)*IT_0462;
    const ccomplex_t IT_0533 = IT_0005 + IT_0532;
    const ccomplex_t IT_0534 = IT_0001*IT_0533;
    const ccomplex_t IT_0535 = (s_12 + -IT_0450)*(s_12 + (-0.5)*IT_0534);
    const ccomplex_t IT_0536 = (-2)*IT_0535;
    const ccomplex_t IT_0537 = (IT_0000 + -IT_0457)*(s_12 + (-0.5)*IT_0534);
    const ccomplex_t IT_0538 = 2*IT_0537;
    const ccomplex_t IT_0539 = (-2)*IT_0142;
    const ccomplex_t IT_0540 = IT_0001*IT_0005;
    const ccomplex_t IT_0541 = (-0.5)*IT_0540;
    const ccomplex_t IT_0542 = s_34*IT_0046;
    const ccomplex_t IT_0543 = s_14*s_23;
    const ccomplex_t IT_0544 = s_12*s_34;
    const ccomplex_t IT_0545 = IT_0543 + IT_0544;
    const ccomplex_t IT_0546 = s_13*s_24;
    const ccomplex_t IT_0547 = -IT_0546;
    const ccomplex_t IT_0548 = IT_0545 + IT_0547;
    const ccomplex_t IT_0549 = IT_0542*IT_0548;
    const ccomplex_t IT_0550 = 0.25*IT_0549;
    const ccomplex_t IT_0551 = s_12 + IT_0541 + IT_0550;
    const ccomplex_t IT_0552 = pow(m_W, 4);
    const ccomplex_t IT_0553 = s_12*IT_0552;
    const ccomplex_t IT_0554 = s_13*s_23*IT_0004;
    const ccomplex_t IT_0555 = s_14*s_23*s_34;
    const ccomplex_t IT_0556 = s_13*s_24*s_34;
    const ccomplex_t IT_0557 = s_12*IT_0064;
    const ccomplex_t IT_0558 = s_14*s_24*IT_0004;
    const ccomplex_t IT_0559 = (-2)*IT_0554 + 2*IT_0555 + 2*IT_0556 + (-2)
      *IT_0557 + (-2)*IT_0558;
    const ccomplex_t IT_0560 = IT_0553 + IT_0559;
    const ccomplex_t IT_0561 = IT_0046*IT_0560;
    const ccomplex_t IT_0562 = (-0.25)*IT_0561;
    const ccomplex_t IT_0563 = -IT_0451;
    const ccomplex_t IT_0564 = s_12 + IT_0507 + IT_0562 + IT_0563;
    const ccomplex_t IT_0565 = IT_0253 + IT_0254;
    const ccomplex_t IT_0566 = s_23*IT_0001*IT_0064;
    const ccomplex_t IT_0567 = (-2)*IT_0566;
    const ccomplex_t IT_0568 = IT_0565 + IT_0567;
    const ccomplex_t IT_0569 = IT_0001*IT_0568;
    const ccomplex_t IT_0570 = -IT_0191 + 0.5*IT_0569;
    const ccomplex_t IT_0571 = s_23 + IT_0570;
    const ccomplex_t IT_0572 = m_N_1*IT_0571;
    const ccomplex_t IT_0573 = 2*IT_0572;
    const ccomplex_t IT_0574 = IT_0001*IT_0253;
    const ccomplex_t IT_0575 = s_24*IT_0004;
    const ccomplex_t IT_0576 = s_23*s_34;
    const ccomplex_t IT_0577 = (-2)*IT_0576;
    const ccomplex_t IT_0578 = IT_0575 + IT_0577;
    const ccomplex_t IT_0579 = IT_0542*IT_0578;
    const ccomplex_t IT_0580 = (-0.5)*IT_0574 + (-0.25)*IT_0579;
    const ccomplex_t IT_0581 = s_23 + IT_0580;
    const ccomplex_t IT_0582 = m_N_1*IT_0581;
    const ccomplex_t IT_0583 = (-4)*IT_0582;
    const ccomplex_t IT_0584 = 2*IT_0239;
    const ccomplex_t IT_0585 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0586 = 0.5*IT_0012 + (-0.5)*IT_0472 + -IT_0585;
    const ccomplex_t IT_0587 = IT_0052 + IT_0586;
    const ccomplex_t IT_0588 = IT_0001*IT_0587;
    const ccomplex_t IT_0589 = IT_0000 + IT_0341 + IT_0588;
    const ccomplex_t IT_0590 = m_N_1*IT_0589;
    const ccomplex_t IT_0591 = 4*IT_0035;
    const ccomplex_t IT_0592 = (-2)*IT_0554 + 4*IT_0555 + (-2)*IT_0558;
    const ccomplex_t IT_0593 = IT_0553 + IT_0592;
    const ccomplex_t IT_0594 = IT_0046*IT_0593;
    const ccomplex_t IT_0595 = 0.25*IT_0594;
    const ccomplex_t IT_0596 = (-0.5)*IT_0534;
    const ccomplex_t IT_0597 = (-0.5)*IT_0009;
    const ccomplex_t IT_0598 = s_12 + IT_0595 + IT_0596 + IT_0597;
    const ccomplex_t IT_0599 = 8*IT_0399;
    const ccomplex_t IT_0600 = (-4)*IT_0335;
    const ccomplex_t IT_0601 = (-2)*IT_0006 + -IT_0462 + 2*IT_0463;
    const ccomplex_t IT_0602 = IT_0005 + IT_0601;
    const ccomplex_t IT_0603 = IT_0001*IT_0602;
    const ccomplex_t IT_0604 = -IT_0451 + (-0.5)*IT_0603;
    const ccomplex_t IT_0605 = s_12 + IT_0604;
    const ccomplex_t IT_0606 = m_N_1*IT_0605;
    const ccomplex_t IT_0607 = (-2)*IT_0606;
    const ccomplex_t IT_0608 = IT_0052 + IT_0473;
    const ccomplex_t IT_0609 = IT_0001*IT_0608;
    const ccomplex_t IT_0610 = IT_0457 + IT_0477;
    const ccomplex_t IT_0611 = (-4)*IT_0610;
    const ccomplex_t IT_0612 = IT_0609 + IT_0611;
    const ccomplex_t IT_0613 = m_N_1*IT_0612;
    const ccomplex_t IT_0614 = IT_0317*IT_0551 + IT_0412*IT_0564 + 0.125
      *IT_0372*IT_0573 + 0.125*IT_0583*IT_0584 + 0.125*IT_0590*IT_0591 + 0.125
      *IT_0598*IT_0599 + 0.125*IT_0360*IT_0600 + 0.125*IT_0036*IT_0607 + 0.125
      *IT_0539*IT_0613;
    const ccomplex_t IT_0615 = IT_0049 + IT_0053 + IT_0068 + IT_0070 + IT_0072
       + IT_0074 + IT_0530;
    const ccomplex_t IT_0616 = IT_0044 + IT_0051 + IT_0057 + IT_0059 + IT_0061
       + IT_0063 + IT_0066;
    const ccomplex_t IT_0617 = s_13*m_N_1;
    const ccomplex_t IT_0618 = IT_0277*IT_0617;
    const ccomplex_t IT_0619 = s_23*m_N_1;
    const ccomplex_t IT_0620 = IT_0277*IT_0619;
    const ccomplex_t IT_0621 = (IT_0005 + (-2)*IT_0006)*(1 + 0.5*IT_0275);
    const ccomplex_t IT_0622 = (-4)*IT_0239;
    const ccomplex_t IT_0623 = s_34*IT_0001;
    const ccomplex_t IT_0624 = IT_0543 + IT_0546;
    const ccomplex_t IT_0625 = s_13*s_23*s_34*IT_0001;
    const ccomplex_t IT_0626 = (-2)*IT_0625;
    const ccomplex_t IT_0627 = IT_0624 + IT_0626;
    const ccomplex_t IT_0628 = IT_0623*IT_0627;
    const ccomplex_t IT_0629 = (-2)*IT_0196;
    const ccomplex_t IT_0630 = 0.5*IT_0002 + (-0.5)*IT_0472 + -IT_0517;
    const ccomplex_t IT_0631 = IT_0050 + IT_0630;
    const ccomplex_t IT_0632 = IT_0001*IT_0631;
    const ccomplex_t IT_0633 = IT_0000 + IT_0337 + IT_0632;
    const ccomplex_t IT_0634 = m_N_1*IT_0633;
    const ccomplex_t IT_0635 = s_12 + IT_0563;
    const ccomplex_t IT_0636 = IT_0000*IT_0635;
    const ccomplex_t IT_0637 = IT_0002 + IT_0473;
    const ccomplex_t IT_0638 = IT_0001*IT_0637;
    const ccomplex_t IT_0639 = IT_0003 + IT_0477;
    const ccomplex_t IT_0640 = (-4)*IT_0639;
    const ccomplex_t IT_0641 = IT_0638 + IT_0640;
    const ccomplex_t IT_0642 = m_N_1*IT_0641;
    const ccomplex_t IT_0643 = 2*IT_0399;
    const ccomplex_t IT_0644 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0645 = (-2)*IT_0006 + -IT_0462 + 2*IT_0644;
    const ccomplex_t IT_0646 = IT_0005 + IT_0645;
    const ccomplex_t IT_0647 = IT_0001*IT_0646;
    const ccomplex_t IT_0648 = -IT_0451 + (-0.5)*IT_0647;
    const ccomplex_t IT_0649 = s_12 + IT_0648;
    const ccomplex_t IT_0650 = m_N_1*IT_0649;
    const ccomplex_t IT_0651 = (-2)*IT_0650;
    const ccomplex_t IT_0652 = 2*IT_0412;
    const ccomplex_t IT_0653 = (s_12 + (-0.5)*IT_0009)*(s_12 + -IT_0451);
    const ccomplex_t IT_0654 = (-2)*IT_0653;
    const ccomplex_t IT_0655 = 4*IT_0040;
    const ccomplex_t IT_0656 = -IT_0462;
    const ccomplex_t IT_0657 = IT_0005 + IT_0656;
    const ccomplex_t IT_0658 = IT_0001*IT_0657;
    const ccomplex_t IT_0659 = -IT_0451 + (-0.25)*IT_0658;
    const ccomplex_t IT_0660 = s_12 + IT_0659;
    const ccomplex_t IT_0661 = m_N_1*IT_0660;
    const ccomplex_t IT_0662 = (-4)*IT_0661;
    const ccomplex_t IT_0663 = (-2)*IT_0223;
    const ccomplex_t IT_0664 = (-0.25)*IT_0377*IT_0454 + (-0.25)*IT_0380
      *IT_0461 + (-0.25)*IT_0386*IT_0488 + (-0.25)*IT_0383*IT_0490 + (-0.25)
      *IT_0390*IT_0513 + (-0.25)*IT_0368*IT_0516 + (-0.25)*IT_0525*IT_0539 +
       IT_0412*IT_0634 + (-1.5)*IT_0035*(IT_0340 + IT_0344 + (-2)*IT_0636) + (
      -0.25)*IT_0642*IT_0643 + (-0.25)*IT_0651*IT_0652 + (-0.25)*IT_0654*IT_0655
       + (-0.25)*IT_0400*IT_0662 + (-0.25)*IT_0526*IT_0663;
    const ccomplex_t IT_0665 = -IT_0006;
    const ccomplex_t IT_0666 = IT_0005 + IT_0665;
    const ccomplex_t IT_0667 = IT_0001*IT_0666;
    const ccomplex_t IT_0668 = -IT_0450 + (-0.25)*IT_0667;
    const ccomplex_t IT_0669 = s_12 + IT_0668;
    const ccomplex_t IT_0670 = m_N_1*IT_0669;
    const ccomplex_t IT_0671 = (-4)*IT_0670;
    const ccomplex_t IT_0672 = IT_0217*IT_0393*IT_0394;
    const ccomplex_t IT_0673 = IT_0220*IT_0396*IT_0397;
    const ccomplex_t IT_0674 = -IT_0672 + -IT_0673;
    const ccomplex_t IT_0675 = 2*IT_0081;
    const ccomplex_t IT_0676 = 2*IT_0183;
    const ccomplex_t IT_0677 = IT_0675 + IT_0676;
    const ccomplex_t IT_0678 = 2*conj(IT_0183);
    const ccomplex_t IT_0679 = IT_0082 + IT_0678;
    const ccomplex_t IT_0680 = IT_0130*IT_0172;
    const ccomplex_t IT_0681 = IT_0139*IT_0396*IT_0680;
    const ccomplex_t IT_0682 = -IT_0681;
    const ccomplex_t IT_0683 = IT_0112*IT_0155;
    const ccomplex_t IT_0684 = IT_0121*IT_0393*IT_0683;
    const ccomplex_t IT_0685 = -IT_0684;
    const ccomplex_t IT_0686 = (-2)*IT_0081;
    const ccomplex_t IT_0687 = (-2)*IT_0230;
    const ccomplex_t IT_0688 = IT_0686 + IT_0687;
    const ccomplex_t IT_0689 = (-2)*conj(IT_0081);
    const ccomplex_t IT_0690 = (-2)*conj(IT_0230);
    const ccomplex_t IT_0691 = IT_0689 + IT_0690;
    const ccomplex_t IT_0692 = IT_0123 + IT_0141;
    const ccomplex_t IT_0693 = m_N_1*IT_0692;
    const ccomplex_t IT_0694 = -IT_0693;
    const ccomplex_t IT_0695 = IT_0217*IT_0393*IT_0683;
    const ccomplex_t IT_0696 = IT_0220*IT_0396*IT_0680;
    const ccomplex_t IT_0697 = m_N_1*IT_0229;
    const ccomplex_t IT_0698 = m_N_1*IT_0227;
    const ccomplex_t IT_0699 = -IT_0695 + -IT_0696 + IT_0697 + IT_0698;
    const ccomplex_t IT_0700 = m_N_1*IT_0182;
    const ccomplex_t IT_0701 = m_N_1*IT_0165;
    const ccomplex_t IT_0702 = (-2)*conj(IT_0223);
    const ccomplex_t IT_0703 = IT_0386*conj(IT_0399) + (conj(IT_0081) + conj
      (IT_0183))*IT_0413 + IT_0399*IT_0417 + (IT_0081 + IT_0183)*IT_0420 + conj
      (IT_0674)*IT_0677 + IT_0674*IT_0679 + (conj(IT_0682) + conj(IT_0685))
      *IT_0688 + (IT_0682 + IT_0685)*IT_0691 + 2*(0.5*IT_0082 + conj(IT_0230))
      *IT_0694 + 2*(IT_0230 + 0.5*IT_0675)*conj(IT_0694) + IT_0416*IT_0699 +
       IT_0383*conj(IT_0699) + IT_0663*(conj(IT_0700) + conj(IT_0701)) + 
      (IT_0700 + IT_0701)*IT_0702;
    const ccomplex_t IT_0704 = pow(m_W, 6);
    const ccomplex_t IT_0705 = s_12*IT_0704;
    const ccomplex_t IT_0706 = s_14*s_24*IT_0552;
    const ccomplex_t IT_0707 = (-2)*IT_0706;
    const ccomplex_t IT_0708 = IT_0705 + IT_0707;
    const ccomplex_t IT_0709 = IT_0046*IT_0708;
    const ccomplex_t IT_0710 = 0.5*IT_0709;
    const ccomplex_t IT_0711 = (-4)*IT_0006;
    const ccomplex_t IT_0712 = IT_0005 + IT_0710 + IT_0711;
    const ccomplex_t IT_0713 = (-2)*IT_0554 + 4*IT_0556 + (-2)*IT_0558;
    const ccomplex_t IT_0714 = IT_0553 + IT_0713;
    const ccomplex_t IT_0715 = IT_0001*IT_0714;
    const ccomplex_t IT_0716 = -IT_0715;
    const ccomplex_t IT_0717 = IT_0712 + IT_0716;
    const ccomplex_t IT_0718 = (-4)*IT_0222;
    const ccomplex_t IT_0719 = (-0.25)*IT_0718;
    const ccomplex_t IT_0720 = IT_0219 + IT_0719;
    const ccomplex_t IT_0721 = (-2)*IT_0254;
    const ccomplex_t IT_0722 = IT_0253 + IT_0721;
    const ccomplex_t IT_0723 = IT_0001*IT_0722;
    const ccomplex_t IT_0724 = s_23*IT_0552;
    const ccomplex_t IT_0725 = s_24*s_34*IT_0004;
    const ccomplex_t IT_0726 = (-2)*IT_0725;
    const ccomplex_t IT_0727 = IT_0724 + IT_0726;
    const ccomplex_t IT_0728 = IT_0046*IT_0727;
    const ccomplex_t IT_0729 = (-0.666666666666667)*IT_0723 +
       0.166666666666667*IT_0728;
    const ccomplex_t IT_0730 = s_23 + IT_0729;
    const ccomplex_t IT_0731 = m_N_1*IT_0730;
    const ccomplex_t IT_0732 = 6*IT_0731;
    const ccomplex_t IT_0733 = (-2)*IT_0674;
    const ccomplex_t IT_0734 = IT_0652 + IT_0733;
    const ccomplex_t IT_0735 = m_N_1*IT_0001;
    const ccomplex_t IT_0736 = IT_0001*IT_0724;
    const ccomplex_t IT_0737 = (-2)*IT_0254 + (-0.5)*IT_0736;
    const ccomplex_t IT_0738 = IT_0253 + IT_0737;
    const ccomplex_t IT_0739 = IT_0735*IT_0738;
    const ccomplex_t IT_0740 = (-2)*IT_0739;
    const ccomplex_t IT_0741 = 2*IT_0700;
    const ccomplex_t IT_0742 = 2*IT_0701;
    const ccomplex_t IT_0743 = IT_0741 + IT_0742;
    const ccomplex_t IT_0744 = IT_0400 + IT_0743;
    const ccomplex_t IT_0745 = s_13*IT_0552;
    const ccomplex_t IT_0746 = s_14*s_34*IT_0004;
    const ccomplex_t IT_0747 = (-2)*IT_0746;
    const ccomplex_t IT_0748 = IT_0745 + IT_0747;
    const ccomplex_t IT_0749 = IT_0046*IT_0748;
    const ccomplex_t IT_0750 = (-0.666666666666667)*IT_0188 + (
      -0.166666666666667)*IT_0749;
    const ccomplex_t IT_0751 = s_13 + IT_0750;
    const ccomplex_t IT_0752 = m_N_1*IT_0751;
    const ccomplex_t IT_0753 = (-6)*IT_0752;
    const ccomplex_t IT_0754 = 2*IT_0682;
    const ccomplex_t IT_0755 = 2*IT_0685;
    const ccomplex_t IT_0756 = (-2)*IT_0694;
    const ccomplex_t IT_0757 = IT_0754 + IT_0755 + IT_0756;
    const ccomplex_t IT_0758 = IT_0542*IT_0575;
    const ccomplex_t IT_0759 = 0.5*IT_0574 + (-0.5)*IT_0723 + (-0.5)*IT_0758;
    const ccomplex_t IT_0760 = s_23 + IT_0759;
    const ccomplex_t IT_0761 = m_N_1*IT_0760;
    const ccomplex_t IT_0762 = 2*IT_0761;
    const ccomplex_t IT_0763 = -IT_0064;
    const ccomplex_t IT_0764 = IT_0552 + IT_0763;
    const ccomplex_t IT_0765 = IT_0001*IT_0764;
    const ccomplex_t IT_0766 = -IT_0004 + 0.25*IT_0765;
    const ccomplex_t IT_0767 = IT_0369 + IT_0766;
    const ccomplex_t IT_0768 = IT_0000*IT_0767;
    const ccomplex_t IT_0769 = 4*IT_0768;
    const ccomplex_t IT_0770 = s_34*IT_0004;
    const ccomplex_t IT_0771 = IT_0542*IT_0770;
    const ccomplex_t IT_0772 = (-2)*IT_0064;
    const ccomplex_t IT_0773 = IT_0552 + IT_0772;
    const ccomplex_t IT_0774 = IT_0001*IT_0773;
    const ccomplex_t IT_0775 = (-0.333333333333333)*IT_0771 + (
      -0.333333333333333)*IT_0774;
    const ccomplex_t IT_0776 = IT_0004 + IT_0775;
    const ccomplex_t IT_0777 = IT_0000*IT_0776;
    const ccomplex_t IT_0778 = 3*IT_0777;
    const ccomplex_t IT_0779 = 2*IT_0078;
    const ccomplex_t IT_0780 = s_13*s_14*s_24;
    const ccomplex_t IT_0781 = s_24*s_34*IT_0000;
    const ccomplex_t IT_0782 = 0.5*IT_0425 + (-0.5)*IT_0426 + -IT_0781;
    const ccomplex_t IT_0783 = IT_0780 + IT_0782;
    const ccomplex_t IT_0784 = IT_0001*IT_0783;
    const ccomplex_t IT_0785 = (-2)*IT_0784;
    const ccomplex_t IT_0786 = (-0.5)*IT_0785;
    const ccomplex_t IT_0787 = IT_0422 + IT_0786;
    const ccomplex_t IT_0788 = 4*IT_0223;
    const ccomplex_t IT_0789 = (-2)*IT_0554;
    const ccomplex_t IT_0790 = IT_0553 + IT_0789;
    const ccomplex_t IT_0791 = IT_0001*IT_0790;
    const ccomplex_t IT_0792 = 0.5*IT_0791;
    const ccomplex_t IT_0793 = (-2)*IT_0554 + 2*IT_0555 + 2*IT_0556 + (-2)
      *IT_0557;
    const ccomplex_t IT_0794 = IT_0553 + IT_0793;
    const ccomplex_t IT_0795 = IT_0001*IT_0794;
    const ccomplex_t IT_0796 = (-0.5)*IT_0795;
    const ccomplex_t IT_0797 = s_14*s_23*IT_0004;
    const ccomplex_t IT_0798 = s_13*s_24*IT_0004;
    const ccomplex_t IT_0799 = IT_0797 + IT_0798;
    const ccomplex_t IT_0800 = s_12*s_34*IT_0004;
    const ccomplex_t IT_0801 = -IT_0800;
    const ccomplex_t IT_0802 = IT_0799 + IT_0801;
    const ccomplex_t IT_0803 = IT_0542*IT_0802;
    const ccomplex_t IT_0804 = 0.5*IT_0803;
    const ccomplex_t IT_0805 = IT_0005 + IT_0007 + IT_0792 + IT_0796 + IT_0804;
    const ccomplex_t IT_0806 = 4*IT_0239;
    const ccomplex_t IT_0807 = (-0.5)*IT_0044;
    const ccomplex_t IT_0808 = IT_0423 + IT_0780;
    const ccomplex_t IT_0809 = -IT_0424 + -IT_0425;
    const ccomplex_t IT_0810 = IT_0808 + IT_0809;
    const ccomplex_t IT_0811 = IT_0001*IT_0810;
    const ccomplex_t IT_0812 = (-0.25)*IT_0811;
    const ccomplex_t IT_0813 = (-0.5)*IT_0056;
    const ccomplex_t IT_0814 = (-0.5)*IT_0058;
    const ccomplex_t IT_0815 = 0.5*IT_0060;
    const ccomplex_t IT_0816 = IT_0054 + IT_0807 + IT_0812 + IT_0813 + IT_0814
       + IT_0815;
    const ccomplex_t IT_0817 = (-8)*IT_0035;
    const ccomplex_t IT_0818 = 0.25*IT_0709;
    const ccomplex_t IT_0819 = IT_0005 + IT_0818;
    const ccomplex_t IT_0820 = (-8)*IT_0268;
    const ccomplex_t IT_0821 = -IT_0049;
    const ccomplex_t IT_0822 = s_14*s_23*s_24;
    const ccomplex_t IT_0823 = IT_0435 + IT_0822;
    const ccomplex_t IT_0824 = -IT_0436 + -IT_0437;
    const ccomplex_t IT_0825 = IT_0823 + IT_0824;
    const ccomplex_t IT_0826 = IT_0001*IT_0825;
    const ccomplex_t IT_0827 = (-0.5)*IT_0826;
    const ccomplex_t IT_0828 = IT_0048 + IT_0068 + IT_0069 + IT_0821 + IT_0827;
    const ccomplex_t IT_0829 = (-4)*IT_0035;
    const ccomplex_t IT_0830 = (-4)*IT_0142;
    const ccomplex_t IT_0831 = s_12*IT_0001*IT_0064;
    const ccomplex_t IT_0832 = -IT_0831;
    const ccomplex_t IT_0833 = IT_0553 + IT_0555 + IT_0556;
    const ccomplex_t IT_0834 = (-2)*IT_0554 + -IT_0557;
    const ccomplex_t IT_0835 = IT_0833 + IT_0834;
    const ccomplex_t IT_0836 = IT_0001*IT_0835;
    const ccomplex_t IT_0837 = (-0.5)*IT_0836;
    const ccomplex_t IT_0838 = -IT_0644;
    const ccomplex_t IT_0839 = IT_0005 + IT_0463 + IT_0832 + IT_0837 + IT_0838;
    const ccomplex_t IT_0840 = (-4)*IT_0196;
    const ccomplex_t IT_0841 = IT_0001*IT_0244;
    const ccomplex_t IT_0842 = s_14*IT_0004;
    const ccomplex_t IT_0843 = IT_0542*IT_0842;
    const ccomplex_t IT_0844 = (-2)*IT_0245;
    const ccomplex_t IT_0845 = IT_0244 + IT_0844;
    const ccomplex_t IT_0846 = IT_0001*IT_0845;
    const ccomplex_t IT_0847 = 0.5*IT_0841 + (-0.5)*IT_0843 + (-0.5)*IT_0846;
    const ccomplex_t IT_0848 = s_13 + IT_0847;
    const ccomplex_t IT_0849 = m_N_1*IT_0848;
    const ccomplex_t IT_0850 = 2*IT_0849;
    const ccomplex_t IT_0851 = IT_0046*IT_0745;
    const ccomplex_t IT_0852 = (-10)*s_13;
    const ccomplex_t IT_0853 = IT_0851 + IT_0852;
    const ccomplex_t IT_0854 = m_N_1*IT_0853;
    const ccomplex_t IT_0855 = (-2)*IT_0699;
    const ccomplex_t IT_0856 = 9*IT_0472;
    const ccomplex_t IT_0857 = (-2)*IT_0262;
    const ccomplex_t IT_0858 = IT_0000*IT_0001;
    const ccomplex_t IT_0859 = (-0.25)*IT_0552;
    const ccomplex_t IT_0860 = IT_0064 + IT_0859;
    const ccomplex_t IT_0861 = IT_0858*IT_0860;
    const ccomplex_t IT_0862 = (-4)*IT_0861;
    const ccomplex_t IT_0863 = (-2)*IT_0264;
    const ccomplex_t IT_0864 = IT_0717*IT_0720 + (-0.25)*IT_0732*IT_0734 + (
      -0.25)*IT_0740*IT_0744 + (-0.25)*IT_0753*IT_0757 + (-0.25)*IT_0368*IT_0762
       + (-0.25)*IT_0384*IT_0769 + (-0.25)*IT_0778*IT_0779 + (-0.25)*IT_0787
      *IT_0788 + (-0.25)*IT_0805*IT_0806 + (-0.25)*IT_0816*IT_0817 + (-0.25)
      *IT_0819*IT_0820 + (-0.25)*IT_0828*IT_0829 + (-0.25)*IT_0444*IT_0830 + (
      -0.25)*IT_0839*IT_0840 + (-0.25)*IT_0377*IT_0850 + (-0.25)*IT_0854*IT_0855
       + (-0.25)*IT_0856*IT_0857 + (-0.25)*IT_0862*IT_0863;
    const ccomplex_t IT_0865 = IT_0378*IT_0461;
    const ccomplex_t IT_0866 = -IT_0544 + -IT_0546;
    const ccomplex_t IT_0867 = IT_0543 + IT_0866;
    const ccomplex_t IT_0868 = IT_0542*IT_0867;
    const ccomplex_t IT_0869 = (-0.25)*IT_0868;
    const ccomplex_t IT_0870 = s_12 + IT_0541 + IT_0869;
    const ccomplex_t IT_0871 = 2*IT_0012 + -IT_0472 + (-2)*IT_0585;
    const ccomplex_t IT_0872 = IT_0052 + IT_0871;
    const ccomplex_t IT_0873 = IT_0001*IT_0872;
    const ccomplex_t IT_0874 = -IT_0457 + 0.5*IT_0873;
    const ccomplex_t IT_0875 = IT_0000 + IT_0874;
    const ccomplex_t IT_0876 = m_N_1*IT_0875;
    const ccomplex_t IT_0877 = 2*IT_0876;
    const ccomplex_t IT_0878 = m_N_1*IT_0259;
    const ccomplex_t IT_0879 = (-4)*IT_0878;
    const ccomplex_t IT_0880 = IT_0046*IT_0714;
    const ccomplex_t IT_0881 = 0.25*IT_0880;
    const ccomplex_t IT_0882 = s_12 + IT_0596 + IT_0597 + IT_0881;
    const ccomplex_t IT_0883 = 8*IT_0412;
    const ccomplex_t IT_0884 = IT_0012 + IT_0473;
    const ccomplex_t IT_0885 = IT_0001*IT_0884;
    const ccomplex_t IT_0886 = IT_0013 + IT_0477;
    const ccomplex_t IT_0887 = (-4)*IT_0886;
    const ccomplex_t IT_0888 = IT_0885 + IT_0887;
    const ccomplex_t IT_0889 = m_N_1*IT_0888;
    const ccomplex_t IT_0890 = -IT_0006 + (-2)*IT_0462 + 2*IT_0644;
    const ccomplex_t IT_0891 = IT_0005 + IT_0890;
    const ccomplex_t IT_0892 = IT_0001*IT_0891;
    const ccomplex_t IT_0893 = -IT_0450 + (-0.5)*IT_0892;
    const ccomplex_t IT_0894 = s_12 + IT_0893;
    const ccomplex_t IT_0895 = m_N_1*IT_0894;
    const ccomplex_t IT_0896 = (-2)*IT_0895;
    const ccomplex_t IT_0897 = (-0.5)*IT_0574 + 0.25*IT_0758;
    const ccomplex_t IT_0898 = s_23 + IT_0897;
    const ccomplex_t IT_0899 = m_N_1*IT_0898;
    const ccomplex_t IT_0900 = 4*IT_0899;
    const ccomplex_t IT_0901 = (-2)*IT_0239;
    const ccomplex_t IT_0902 = 0.125*IT_0564*IT_0599 + 0.125*IT_0360*IT_0600 +
       0.125*IT_0143*IT_0662 + IT_0317*IT_0870 + 0.125*IT_0383*IT_0877 + 0.125
      *IT_0372*IT_0879 + 0.125*IT_0882*IT_0883 + 0.125*IT_0036*IT_0889 + 0.125
      *IT_0663*IT_0896 + 0.125*IT_0900*IT_0901;
    const ccomplex_t IT_0903 = 8*conj(IT_0412);
    const ccomplex_t IT_0904 = 1.125*IT_0000;
    const ccomplex_t IT_0905 = (-0.5)*IT_0841 + 0.25*IT_0843;
    const ccomplex_t IT_0906 = s_13 + IT_0905;
    const ccomplex_t IT_0907 = m_N_1*IT_0906;
    const ccomplex_t IT_0908 = 4*IT_0907;
    const ccomplex_t IT_0909 = 2*conj(IT_0078);
    const ccomplex_t IT_0910 = 4*conj(IT_0040);
    const ccomplex_t IT_0911 = (-2)*conj(IT_0040);
    const ccomplex_t IT_0912 = (-2)*conj(IT_0196);
    const ccomplex_t IT_0913 = IT_0046*IT_0773;
    const ccomplex_t IT_0914 = 0.25*IT_0913;
    const ccomplex_t IT_0915 = 1 + IT_0914;
    const ccomplex_t IT_0916 = IT_0000*IT_0915;
    const ccomplex_t IT_0917 = IT_0360*IT_0362 + IT_0416*IT_0470 + 8*conj
      (IT_0335)*IT_0551 + IT_0496*IT_0573 + IT_0145*IT_0642 + IT_0528*IT_0662 +
       IT_0613*IT_0689 + IT_0524*IT_0702 + (-16)*conj(IT_0399)*IT_0904 + IT_0447
      *IT_0908 + IT_0583*IT_0909 + IT_0590*IT_0910 + IT_0607*IT_0911 + IT_0392
      *IT_0912 + IT_0903*IT_0916;
    const ccomplex_t IT_0918 = 0.25*IT_0620;
    const ccomplex_t IT_0919 = 0.25*IT_0628;
    const ccomplex_t IT_0920 = IT_0277*IT_0472;
    const ccomplex_t IT_0921 = (-2)*IT_0201;
    const ccomplex_t IT_0922 = 0.25*IT_0385;
    const ccomplex_t IT_0923 = 0.25*IT_0382;
    const ccomplex_t IT_0924 = 0.25*IT_0908;
    const ccomplex_t IT_0925 = s_13*s_34;
    const ccomplex_t IT_0926 = (-2)*IT_0925;
    const ccomplex_t IT_0927 = IT_0842 + IT_0926;
    const ccomplex_t IT_0928 = IT_0542*IT_0927;
    const ccomplex_t IT_0929 = (-0.5)*IT_0841 + (-0.25)*IT_0928;
    const ccomplex_t IT_0930 = s_13 + IT_0929;
    const ccomplex_t IT_0931 = m_N_1*IT_0930;
    const ccomplex_t IT_0932 = (-4)*IT_0931;
    const ccomplex_t IT_0933 = 0.25*IT_0932;
    const ccomplex_t IT_0934 = IT_0317*IT_0618 + -IT_0078*IT_0621 + 0.25*
      (IT_0044 + IT_0045 + IT_0047 + IT_0048 + IT_0049 + IT_0051 + IT_0053 +
       IT_0057 + IT_0059 + IT_0061 + IT_0063 + IT_0066 + IT_0068 + IT_0070 +
       IT_0072 + IT_0074)*IT_0675 + IT_0600*IT_0918 + 0.25*IT_0622*IT_0920 +
       IT_0919*IT_0921 + IT_0386*IT_0922 + IT_0383*IT_0923 + IT_0643*IT_0924 +
       IT_0413*IT_0933;
    const ccomplex_t IT_0935 = (-4)*conj(IT_0335);
    const ccomplex_t IT_0936 = (-2)*conj(IT_0201);
    const ccomplex_t IT_0937 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0938 = 0.5*IT_0437 + (-0.5)*IT_0438 + -IT_0937;
    const ccomplex_t IT_0939 = IT_0822 + IT_0938;
    const ccomplex_t IT_0940 = IT_0001*IT_0939;
    const ccomplex_t IT_0941 = (-2)*IT_0940;
    const ccomplex_t IT_0942 = (-0.5)*IT_0941;
    const ccomplex_t IT_0943 = IT_0434 + IT_0942;
    const ccomplex_t IT_0944 = (-0.25)*IT_0762;
    const ccomplex_t IT_0945 = (-0.25)*IT_0850;
    const ccomplex_t IT_0946 = IT_0377*IT_0945;
    const ccomplex_t IT_0947 = (-0.25)*IT_0769;
    const ccomplex_t IT_0948 = IT_0384*IT_0947;
    const ccomplex_t IT_0949 = (-0.25)*IT_0778;
    const ccomplex_t IT_0950 = IT_0779*IT_0949;
    const ccomplex_t IT_0951 = IT_0046*IT_0724;
    const ccomplex_t IT_0952 = 0.666666666666667*IT_0574 + (-0.166666666666667
      )*IT_0951;
    const ccomplex_t IT_0953 = s_23 + IT_0952;
    const ccomplex_t IT_0954 = m_N_1*IT_0953;
    const ccomplex_t IT_0955 = (-6)*IT_0954;
    const ccomplex_t IT_0956 = (-0.25)*IT_0955;
    const ccomplex_t IT_0957 = 0.166666666666667*IT_0749 + (-0.666666666666667
      )*IT_0846;
    const ccomplex_t IT_0958 = s_13 + IT_0957;
    const ccomplex_t IT_0959 = m_N_1*IT_0958;
    const ccomplex_t IT_0960 = 6*IT_0959;
    const ccomplex_t IT_0961 = (-0.25)*IT_0960;
    const ccomplex_t IT_0962 = (-0.666666666666667)*IT_0191 + (
      -0.166666666666667)*IT_0728;
    const ccomplex_t IT_0963 = s_23 + IT_0962;
    const ccomplex_t IT_0964 = m_N_1*IT_0963;
    const ccomplex_t IT_0965 = (-6)*IT_0964;
    const ccomplex_t IT_0966 = (-0.25)*IT_0965;
    const ccomplex_t IT_0967 = IT_0001*IT_0745;
    const ccomplex_t IT_0968 = (-2)*IT_0245 + (-0.5)*IT_0967;
    const ccomplex_t IT_0969 = IT_0244 + IT_0968;
    const ccomplex_t IT_0970 = IT_0735*IT_0969;
    const ccomplex_t IT_0971 = (-2)*IT_0970;
    const ccomplex_t IT_0972 = (-0.25)*IT_0971;
    const ccomplex_t IT_0973 = (-10)*s_23;
    const ccomplex_t IT_0974 = IT_0951 + IT_0973;
    const ccomplex_t IT_0975 = m_N_1*IT_0974;
    const ccomplex_t IT_0976 = (-0.25)*IT_0975;
    const ccomplex_t IT_0977 = (-0.25)*IT_0805;
    const ccomplex_t IT_0978 = -IT_0463;
    const ccomplex_t IT_0979 = IT_0005 + IT_0644 + IT_0832 + IT_0837 + IT_0978;
    const ccomplex_t IT_0980 = (-0.25)*IT_0979;
    const ccomplex_t IT_0981 = 0.5*IT_0811;
    const ccomplex_t IT_0982 = IT_0044 + IT_0059 + IT_0060 + IT_0421 + IT_0981;
    const ccomplex_t IT_0983 = (-0.25)*IT_0982;
    const ccomplex_t IT_0984 = (-0.5)*IT_0049;
    const ccomplex_t IT_0985 = (-0.25)*IT_0826;
    const ccomplex_t IT_0986 = (-0.5)*IT_0048;
    const ccomplex_t IT_0987 = 0.5*IT_0067;
    const ccomplex_t IT_0988 = (-0.5)*IT_0069;
    const ccomplex_t IT_0989 = IT_0045 + IT_0984 + IT_0985 + IT_0986 + IT_0987
       + IT_0988;
    const ccomplex_t IT_0990 = (-0.25)*IT_0989;
    const ccomplex_t IT_0991 = IT_0223*IT_0432 + -IT_0142*IT_0943 + IT_0368
      *IT_0944 + IT_0946 + IT_0948 + IT_0950 + IT_0400*IT_0956 + IT_0757*IT_0961
       + IT_0734*IT_0966 + IT_0855*IT_0972 + IT_0743*IT_0976 + IT_0806*IT_0977 +
       IT_0840*IT_0980 + IT_0591*IT_0983 + IT_0817*IT_0990;
    const ccomplex_t IT_0992 = IT_0001*IT_0593;
    const ccomplex_t IT_0993 = -IT_0992;
    const ccomplex_t IT_0994 = IT_0712 + IT_0993;
    const ccomplex_t IT_0995 = 2*conj(IT_0682);
    const ccomplex_t IT_0996 = 2*conj(IT_0685);
    const ccomplex_t IT_0997 = (-2)*conj(IT_0694);
    const ccomplex_t IT_0998 = IT_0995 + IT_0996 + IT_0997;
    const ccomplex_t IT_0999 = (-2)*conj(IT_0674);
    const ccomplex_t IT_1000 = IT_0482 + IT_0999;
    const ccomplex_t IT_1001 = (-2)*conj(IT_0699);
    const ccomplex_t IT_1002 = 2*conj(IT_0700);
    const ccomplex_t IT_1003 = 2*conj(IT_0701);
    const ccomplex_t IT_1004 = IT_1002 + IT_1003;
    const ccomplex_t IT_1005 = 4*conj(IT_0239);
    const ccomplex_t IT_1006 = (-4)*conj(IT_0196);
    const ccomplex_t IT_1007 = 4*conj(IT_0035);
    const ccomplex_t IT_1008 = (-8)*conj(IT_0035);
    const ccomplex_t IT_1009 = conj(IT_0223)*IT_0432 + -conj(IT_0142)*IT_0943 
      + IT_0378*IT_0944 + IT_0379*IT_0945 + IT_0496*IT_0947 + IT_0909*IT_0949 +
       IT_0419*IT_0956 + conj(IT_0268)*IT_0994 + IT_0961*IT_0998 + IT_0966
      *IT_1000 + IT_0972*IT_1001 + IT_0976*IT_1004 + IT_0977*IT_1005 + IT_0980
      *IT_1006 + IT_0983*IT_1007 + IT_0990*IT_1008;
    const ccomplex_t IT_1010 = -IT_0943;
    const ccomplex_t IT_1011 = 4*IT_0081;
    const ccomplex_t IT_1012 = 0.25*IT_1011;
    const ccomplex_t IT_1013 = (-2)*IT_0268;
    const ccomplex_t IT_1014 = (-0.25)*IT_0856;
    const ccomplex_t IT_1015 = 4*IT_0078;
    const ccomplex_t IT_1016 = IT_0977*IT_1015;
    const ccomplex_t IT_1017 = (-4)*IT_0201;
    const ccomplex_t IT_1018 = (-8)*IT_0040;
    const ccomplex_t IT_1019 = IT_0230*IT_0432 + IT_0380*IT_0944 + IT_0390
      *IT_0945 + IT_0372*IT_0947 + IT_0584*IT_0949 + IT_0744*IT_0961 + IT_0855
      *IT_0966 + IT_0734*IT_0972 + IT_0757*IT_0976 + IT_0655*IT_0983 + IT_1010*
      (IT_0183 + IT_1012) + IT_1013*IT_1014 + IT_1016 + IT_0980*IT_1017 +
       IT_0990*IT_1018;
    const ccomplex_t IT_1020 = (-0.25)*IT_0740;
    const ccomplex_t IT_1021 = (-0.25)*IT_0753;
    const ccomplex_t IT_1022 = (-0.25)*IT_0732;
    const ccomplex_t IT_1023 = (-0.25)*IT_0854;
    const ccomplex_t IT_1024 = (-0.25)*IT_0862;
    const ccomplex_t IT_1025 = (-0.25)*IT_0839;
    const ccomplex_t IT_1026 = (-0.25)*IT_0816;
    const ccomplex_t IT_1027 = (-0.25)*IT_0828;
    const ccomplex_t IT_1028 = IT_0717*IT_0720 + (-0.25)*IT_0787*IT_0788 + (
      -0.25)*IT_0819*IT_0820 + (-0.25)*IT_0444*IT_0830 + IT_0368*IT_0944 +
       IT_0946 + IT_0948 + IT_0950 + IT_0806*IT_0977 + IT_0857*IT_1014 + IT_0744
      *IT_1020 + IT_0757*IT_1021 + IT_0734*IT_1022 + IT_0855*IT_1023 + IT_0863
      *IT_1024 + IT_0840*IT_1025 + IT_0817*IT_1026 + IT_0829*IT_1027;
    const ccomplex_t IT_1029 = -IT_0787;
    const ccomplex_t IT_1030 = 0.666666666666667*IT_0841 + (-0.166666666666667
      )*IT_0851;
    const ccomplex_t IT_1031 = s_13 + IT_1030;
    const ccomplex_t IT_1032 = m_N_1*IT_1031;
    const ccomplex_t IT_1033 = (-6)*IT_1032;
    const ccomplex_t IT_1034 = (-0.25)*IT_1033;
    const ccomplex_t IT_1035 = (-4)*IT_0040;
    const ccomplex_t IT_1036 = IT_0183*IT_0444 + IT_0380*IT_0944 + IT_0390
      *IT_0945 + IT_0372*IT_0947 + IT_0584*IT_0949 + IT_1016 + IT_0757*IT_1020 +
       IT_0744*IT_1021 + IT_0855*IT_1022 + IT_0652*IT_1023 + IT_1013*IT_1024 +
       IT_1017*IT_1025 + IT_1018*IT_1026 + (IT_0230 + IT_1012)*IT_1029 + IT_0733
      *IT_1034 + IT_1027*IT_1035;
    const ccomplex_t IT_1037 = 4*conj(IT_0081);
    const ccomplex_t IT_1038 = 0.25*IT_1037;
    const ccomplex_t IT_1039 = IT_0415*IT_0944;
    const ccomplex_t IT_1040 = IT_0418*IT_0945;
    const ccomplex_t IT_1041 = IT_0446*IT_0947;
    const ccomplex_t IT_1042 = IT_0447*IT_0949;
    const ccomplex_t IT_1043 = IT_0419 + IT_1004;
    const ccomplex_t IT_1044 = (-2)*conj(IT_0268);
    const ccomplex_t IT_1045 = 4*conj(IT_0078);
    const ccomplex_t IT_1046 = IT_0977*IT_1045;
    const ccomplex_t IT_1047 = (-4)*conj(IT_0201);
    const ccomplex_t IT_1048 = (-8)*conj(IT_0040);
    const ccomplex_t IT_1049 = (-4)*conj(IT_0040);
    const ccomplex_t IT_1050 = (-4)*conj(IT_0081);
    const ccomplex_t IT_1051 = (-0.25)*IT_1050;
    const ccomplex_t IT_1052 = conj(IT_0262)*IT_0717 + 2*conj(IT_0264)*IT_0819
       + IT_0998*IT_1020 + IT_1001*IT_1022 + IT_0482*IT_1023 + IT_0999*IT_1034 +
       IT_1029*(conj(IT_0230) + IT_1038) + IT_1039 + IT_1040 + IT_1041 + IT_1042
       + IT_1021*IT_1043 + IT_1024*IT_1044 + IT_1046 + IT_1025*IT_1047 + IT_1026
      *IT_1048 + IT_1027*IT_1049 + IT_0444*(conj(IT_0183) + IT_1051);
    const ccomplex_t IT_1053 = 2*conj(IT_0262)*IT_0819 + IT_0910*IT_0983 +
       conj(IT_0264)*IT_0994 + IT_0976*IT_0998 + IT_0972*IT_1000 + IT_0966
      *IT_1001 + IT_1010*(conj(IT_0183) + IT_1038) + IT_1039 + IT_1040 + IT_1041
       + IT_1042 + IT_0961*IT_1043 + IT_1014*IT_1044 + IT_1046 + IT_0980*IT_1047
       + IT_0990*IT_1048 + IT_0432*(conj(IT_0230) + IT_1051);
    const ccomplex_t IT_1054 = 0.25*IT_0607;
    const ccomplex_t IT_1055 = IT_0147*IT_1054;
    const ccomplex_t IT_1056 = 0.25*IT_0662;
    const ccomplex_t IT_1057 = 0.25*IT_0642;
    const ccomplex_t IT_1058 = IT_0036*IT_1057;
    const ccomplex_t IT_1059 = 0.25*IT_0590;
    const ccomplex_t IT_1060 = IT_1035*IT_1059;
    const ccomplex_t IT_1061 = 0.25*IT_0573;
    const ccomplex_t IT_1062 = IT_0921*IT_1061;
    const ccomplex_t IT_1063 = 0.25*IT_0392;
    const ccomplex_t IT_1064 = IT_0372*IT_1063;
    const ccomplex_t IT_1065 = 0.25*IT_0524;
    const ccomplex_t IT_1066 = IT_0386*IT_1065;
    const ccomplex_t IT_1067 = 0.25*IT_0613;
    const ccomplex_t IT_1068 = IT_0677*IT_1067;
    const ccomplex_t IT_1069 = 0.25*IT_0470;
    const ccomplex_t IT_1070 = IT_0539*IT_1069;
    const ccomplex_t IT_1071 = (-2)*IT_0078;
    const ccomplex_t IT_1072 = 0.25*IT_0583;
    const ccomplex_t IT_1073 = IT_1071*IT_1072;
    const ccomplex_t IT_1074 = IT_0901*IT_0924;
    const ccomplex_t IT_1075 = 16*IT_0399;
    const ccomplex_t IT_1076 = (-16)*IT_0700;
    const ccomplex_t IT_1077 = (-16)*IT_0701;
    const ccomplex_t IT_1078 = IT_1075 + IT_1076 + IT_1077;
    const ccomplex_t IT_1079 = 0.25*IT_0904;
    const ccomplex_t IT_1080 = IT_1078*IT_1079;
    const ccomplex_t IT_1081 = 8*IT_0674;
    const ccomplex_t IT_1082 = (-8)*IT_0412;
    const ccomplex_t IT_1083 = IT_1081 + IT_1082;
    const ccomplex_t IT_1084 = 0.25*IT_0916;
    const ccomplex_t IT_1085 = IT_1083*IT_1084;
    const ccomplex_t IT_1086 = (-8)*IT_0335;
    const ccomplex_t IT_1087 = 0.25*IT_0551;
    const ccomplex_t IT_1088 = IT_1086*IT_1087;
    const ccomplex_t IT_1089 = 8*IT_0682;
    const ccomplex_t IT_1090 = (-8)*IT_0694;
    const ccomplex_t IT_1091 = 0.25*IT_0598;
    const ccomplex_t IT_1092 = IT_0317*IT_0360;
    const ccomplex_t IT_1093 = 4*conj(IT_0682)*(IT_1055 + IT_0143*IT_1056 +
       IT_1058 + IT_1060 + IT_1062 + IT_1064 + IT_1066 + IT_1068 + IT_1070 +
       IT_1073 + IT_1074 + IT_1080 + IT_1085 + IT_1088 + (IT_1089 + IT_1090)
      *IT_1091 + IT_1092);
    const ccomplex_t IT_1094 = 0.125*IT_1090;
    const ccomplex_t IT_1095 = 0.125*IT_1089;
    const ccomplex_t IT_1096 = IT_0685 + IT_1094 + IT_1095;
    const ccomplex_t IT_1097 = 0.25*IT_0634;
    const ccomplex_t IT_1098 = 0.25*IT_0889;
    const ccomplex_t IT_1099 = 0.25*IT_0651;
    const ccomplex_t IT_1100 = 0.25*IT_0879;
    const ccomplex_t IT_1101 = 0.25*IT_0409;
    const ccomplex_t IT_1102 = 2*IT_0230;
    const ccomplex_t IT_1103 = IT_0675 + IT_1102;
    const ccomplex_t IT_1104 = 0.25*IT_0896;
    const ccomplex_t IT_1105 = 0.25*IT_0481;
    const ccomplex_t IT_1106 = (-2)*IT_0183;
    const ccomplex_t IT_1107 = IT_0686 + IT_1106;
    const ccomplex_t IT_1108 = 0.25*IT_0877;
    const ccomplex_t IT_1109 = 0.25*IT_0900;
    const ccomplex_t IT_1110 = 8*IT_0700;
    const ccomplex_t IT_1111 = 8*IT_0701;
    const ccomplex_t IT_1112 = (-8)*IT_0399;
    const ccomplex_t IT_1113 = 0.25*IT_0870;
    const ccomplex_t IT_1114 = 8*IT_0699;
    const ccomplex_t IT_1115 = 0.25*IT_0882;
    const ccomplex_t IT_1116 = (IT_0412 + -IT_0674)*IT_0904 + 0.25*IT_0584
      *IT_0933 + 0.25*IT_0041*IT_1056 + 0.25*IT_1092 + 0.5*IT_0564*IT_1096 +
       0.25*IT_0591*IT_1097 + 0.25*IT_0147*IT_1098 + 0.25*IT_0036*IT_1099 + 0.25
      *IT_0921*IT_1100 + 0.25*IT_0372*IT_1101 + 0.25*IT_1103*IT_1104 + 0.25
      *IT_0663*IT_1105 + 0.25*IT_1107*IT_1108 + 0.25*IT_0779*IT_1109 + 0.25
      *IT_1084*(IT_1110 + IT_1111 + IT_1112) + 0.25*IT_1086*IT_1113 + 0.25
      *IT_1114*IT_1115;
    const ccomplex_t IT_1117 = 0.25*IT_0564;
    const ccomplex_t IT_1118 = 0.5*IT_1055 + 0.5*IT_0143*IT_1056 + 0.5*IT_1058
       + 0.5*IT_1060 + 0.5*IT_1062 + 0.5*IT_1064 + 0.5*IT_1066 + 0.5*IT_1068 +
       0.5*IT_1070 + 0.5*IT_1073 + 0.5*IT_1074 + 0.5*IT_1080 + 0.5*IT_1085 + 0.5
      *IT_1088 + 0.5*IT_1092 + IT_0598*IT_1096 + 0.5*IT_1114*IT_1117;
    const ccomplex_t IT_1119 = IT_0185*IT_1054;
    const ccomplex_t IT_1120 = IT_0145*IT_1056;
    const ccomplex_t IT_1121 = IT_0528*IT_1057;
    const ccomplex_t IT_1122 = IT_1049*IT_1059;
    const ccomplex_t IT_1123 = IT_0936*IT_1061;
    const ccomplex_t IT_1124 = IT_0446*IT_1063;
    const ccomplex_t IT_1125 = IT_0417*IT_1065;
    const ccomplex_t IT_1126 = IT_0679*IT_1067;
    const ccomplex_t IT_1127 = (-2)*conj(IT_0142);
    const ccomplex_t IT_1128 = IT_1069*IT_1127;
    const ccomplex_t IT_1129 = (-2)*conj(IT_0078);
    const ccomplex_t IT_1130 = IT_1072*IT_1129;
    const ccomplex_t IT_1131 = (-2)*conj(IT_0239);
    const ccomplex_t IT_1132 = IT_0924*IT_1131;
    const ccomplex_t IT_1133 = 16*conj(IT_0399);
    const ccomplex_t IT_1134 = (-16)*conj(IT_0700);
    const ccomplex_t IT_1135 = (-16)*conj(IT_0701);
    const ccomplex_t IT_1136 = IT_1133 + IT_1134 + IT_1135;
    const ccomplex_t IT_1137 = IT_1079*IT_1136;
    const ccomplex_t IT_1138 = 8*conj(IT_0674);
    const ccomplex_t IT_1139 = (-8)*conj(IT_0412);
    const ccomplex_t IT_1140 = IT_1138 + IT_1139;
    const ccomplex_t IT_1141 = IT_1084*IT_1140;
    const ccomplex_t IT_1142 = (-8)*conj(IT_0335);
    const ccomplex_t IT_1143 = IT_1087*IT_1142;
    const ccomplex_t IT_1144 = 8*conj(IT_0682);
    const ccomplex_t IT_1145 = (-8)*conj(IT_0694);
    const ccomplex_t IT_1146 = IT_1144 + IT_1145;
    const ccomplex_t IT_1147 = conj(IT_0317)*IT_0360;
    const ccomplex_t IT_1148 = 8*conj(IT_0700);
    const ccomplex_t IT_1149 = (-8)*conj(IT_0399);
    const ccomplex_t IT_1150 = 2*conj(IT_0230);
    const ccomplex_t IT_1151 = IT_0082 + IT_1150;
    const ccomplex_t IT_1152 = (-2)*conj(IT_0183);
    const ccomplex_t IT_1153 = IT_0689 + IT_1152;
    const ccomplex_t IT_1154 = (conj(IT_0412) + -conj(IT_0674))*IT_0904 + 0.25
      *IT_0447*IT_0933 + 0.25*IT_0911*IT_1056 + 0.25*IT_1007*IT_1097 + 0.25
      *IT_0185*IT_1098 + 0.25*IT_0528*IT_1099 + 0.25*IT_0936*IT_1100 + 0.25
      *IT_0446*IT_1101 + 0.25*IT_0702*IT_1105 + 0.25*IT_0909*IT_1109 + 0.25
      *IT_1113*IT_1142 + 0.25*IT_1117*IT_1146 + 0.25*IT_1147 + 0.5*IT_0916*(conj
      (IT_0701) + 0.125*IT_1148 + 0.125*IT_1149) + 0.25*IT_1104*IT_1151 + 0.25
      *IT_1108*IT_1153;
    const ccomplex_t IT_1155 = IT_0185*IT_1056;
    const ccomplex_t IT_1156 = IT_0911*IT_1057;
    const ccomplex_t IT_1157 = (-4)*conj(IT_0035);
    const ccomplex_t IT_1158 = IT_0496*IT_1063;
    const ccomplex_t IT_1159 = IT_1065*IT_1151;
    const ccomplex_t IT_1160 = IT_1069*IT_1153;
    const ccomplex_t IT_1161 = IT_0924*IT_1129;
    const ccomplex_t IT_1162 = 16*conj(IT_0694);
    const ccomplex_t IT_1163 = (-8)*conj(IT_0317);
    const ccomplex_t IT_1164 = IT_1148 + IT_1149;
    const ccomplex_t IT_1165 = conj(IT_0335)*IT_0360;
    const ccomplex_t IT_1166 = IT_0145*IT_1054 + IT_0912*IT_1061 + IT_0416
      *IT_1067 + IT_1072*IT_1131 + IT_1117*IT_1140 + IT_1155 + IT_1156 + IT_1059
      *IT_1157 + IT_1158 + IT_1159 + IT_1160 + IT_1161 + IT_1079*IT_1162 +
       IT_1087*IT_1163 + IT_1091*IT_1164 + IT_1165;
    const ccomplex_t IT_1167 = IT_0145*IT_1054 + IT_0912*IT_1061 + IT_0416
      *IT_1067 + IT_1072*IT_1131 + IT_1117*IT_1139 + IT_1091*IT_1149 + IT_1155 +
       IT_1156 + IT_1059*IT_1157 + IT_1158 + IT_1159 + IT_1160 + IT_1161 +
       IT_1079*IT_1162 + IT_1087*IT_1163 + IT_1165;
    const ccomplex_t IT_1168 = IT_0909*IT_0933 + IT_0528*IT_1056 + IT_0910
      *IT_1097 + IT_0145*IT_1098 + IT_0911*IT_1099 + IT_0912*IT_1100 + IT_0496
      *IT_1101 + IT_0417*IT_1104 + IT_0691*IT_1105 + IT_0447*IT_1109 + IT_1108
      *IT_1127 + IT_1115*IT_1139 + IT_1084*IT_1145 + IT_1113*IT_1163 + IT_1117
      *IT_1164 + IT_1165;
    const ccomplex_t IT_1169 = 8*conj(IT_0317);
    const ccomplex_t IT_1170 = IT_0360*IT_0935;
    const ccomplex_t IT_1171 = IT_0372 + IT_0584;
    const ccomplex_t IT_1172 = 6*IT_0142;
    const ccomplex_t IT_1173 = (-6)*IT_0223;
    const ccomplex_t IT_1174 = IT_0380*IT_0461;
    const ccomplex_t IT_1175 = IT_0368*IT_0454 + IT_0143*IT_0490 + IT_0036
      *IT_0525 + IT_0384*IT_0531 + IT_0470*IT_0643 + IT_0613*IT_0756 + IT_0615
      *IT_0779 + IT_0536*IT_1103 + IT_0538*IT_1107 + IT_0382*IT_1171 + IT_0500
      *IT_1172 + IT_0509*IT_1173 + IT_1174;
    const ccomplex_t IT_1176 = IT_0415*IT_0461;
    const ccomplex_t IT_1177 = IT_0382*IT_0448 + IT_0378*IT_0454 + IT_0470
      *IT_0471 + IT_0145*IT_0490 + IT_0500*IT_0501 + IT_0509*IT_0511 + IT_0525
      *IT_0528 + IT_0496*IT_0531 + IT_0538*IT_0689 + IT_0615*IT_0909 + IT_0613
      *IT_0997 + IT_0536*IT_1151 + IT_1176;
    const ccomplex_t IT_1178 = (IT_0000 + -IT_0476)*(s_12 + (-0.5)*IT_0534);
    const ccomplex_t IT_1179 = 2*IT_1178;
    const ccomplex_t IT_1180 = IT_0051 + IT_0057 + IT_0063 + IT_0066 + IT_0493;
    const ccomplex_t IT_1181 = 6*IT_0223;
    const ccomplex_t IT_1182 = (-6)*IT_0142;
    const ccomplex_t IT_1183 = IT_0377*IT_0516;
    const ccomplex_t IT_1184 = IT_0143*IT_0488 + IT_0390*IT_0513 + IT_0400
      *IT_0524 + IT_0036*IT_0526 + IT_0481*IT_0652 + IT_0536*IT_0675 + IT_0616
      *IT_0779 + IT_0385*IT_1171 + IT_0688*IT_1179 + IT_0384*IT_1180 + IT_0505
      *IT_1181 + IT_0509*IT_1182 + IT_1183;
    const ccomplex_t IT_1185 = IT_0368*IT_0454 + IT_0143*IT_0491 + IT_0384*
      (IT_0045 + IT_0047 + IT_0048 + IT_0051 + IT_0053 + IT_0057 + IT_0063 +
       IT_0066 + IT_0072 + IT_0074 + IT_0493 + IT_0495) + IT_0390*IT_0513 +
       IT_0400*IT_0524 + IT_0036*IT_0527 + IT_0470*IT_0643 + IT_0481*IT_0652 +
       IT_0536*IT_1011 + IT_0445*IT_1171 + IT_0500*IT_1172 + IT_1174 + IT_0686*
      (IT_0538 + IT_1179) + IT_0505*IT_1181 + IT_0509*(IT_1173 + IT_1182) +
       IT_1183;
    const ccomplex_t IT_1186 = IT_0418*IT_0513;
    const ccomplex_t IT_1187 = IT_0385*IT_0448 + IT_0481*IT_0482 + IT_0145
      *IT_0488 + IT_0505*IT_0506 + IT_0509*IT_0510 + IT_0379*IT_0516 + IT_0419
      *IT_0524 + IT_0526*IT_0528 + IT_0082*IT_0536 + IT_0616*IT_0909 + IT_0689
      *IT_1179 + IT_0496*IT_1180 + IT_1186;
    const ccomplex_t IT_1188 = IT_0379*IT_0454 + IT_0417*IT_0488 + IT_0416
      *IT_0490 + IT_0378*IT_0516 + (-4)*conj(IT_0412)*IT_0634 + (-12)*conj
      (IT_0035)*IT_0636 + IT_0471*IT_0642 + IT_0482*IT_0651 + IT_0419*IT_0662 +
       IT_0526*IT_0702 + IT_0525*IT_1127 + IT_1176 + IT_1186;
    const ccomplex_t IT_1189 = IT_0377*IT_0513;
    const ccomplex_t IT_1190 = IT_0379*IT_0513;
    const ccomplex_t IT_1191 = IT_0036*IT_0662;
    const ccomplex_t IT_1192 = IT_0360*IT_0361 + IT_0383*IT_0470 + 8*IT_0335
      *IT_0551 + IT_0384*IT_0573 + IT_0041*IT_0607 + IT_0392*IT_0629 + IT_0143
      *IT_0642 + IT_0590*IT_0655 + IT_0524*IT_0663 + IT_0613*IT_0686 + 8*IT_0598
      *IT_0694 + IT_0583*IT_0779 + (-16)*IT_0399*IT_0904 + IT_0584*IT_0908 +
       IT_0883*IT_0916 + IT_1191;
    const ccomplex_t IT_1193 = IT_1110 + IT_1112;
    const ccomplex_t IT_1194 = (-8)*IT_0317;
    const ccomplex_t IT_1195 = IT_0336*IT_0360;
    const ccomplex_t IT_1196 = IT_0384*IT_0409 + IT_0041*IT_0651 + IT_0634
      *IT_0655 + IT_0481*IT_0688 + IT_0539*IT_0877 + IT_0629*IT_0879 + IT_0143
      *IT_0889 + IT_0386*IT_0896 + IT_0584*IT_0900 + IT_0779*IT_0932 + IT_0882
      *IT_1083 + IT_0916*IT_1090 + IT_1191 + IT_0564*IT_1193 + IT_0870*IT_1194 +
       IT_1195;
    const ccomplex_t IT_1197 = IT_0694*IT_0904;
    const ccomplex_t IT_1198 = 16*IT_1197;
    const ccomplex_t IT_1199 = IT_0147*IT_0662;
    const ccomplex_t IT_1200 = IT_0041*IT_0642;
    const ccomplex_t IT_1201 = IT_0590*IT_0829;
    const ccomplex_t IT_1202 = IT_0573*IT_0629;
    const ccomplex_t IT_1203 = IT_0470*IT_1107;
    const ccomplex_t IT_1204 = IT_0583*IT_0901;
    const ccomplex_t IT_1205 = IT_0908*IT_1071;
    const ccomplex_t IT_1206 = IT_0551*IT_1194;
    const ccomplex_t IT_1207 = IT_0384*IT_0392 + IT_0143*IT_0607 + IT_0383
      *IT_0613 + IT_0564*IT_1083 + IT_0524*IT_1103 + IT_0598*(IT_1110 + IT_1111 
      + IT_1112) + IT_1195 + IT_1198 + IT_1199 + IT_1200 + IT_1201 + IT_1202 +
       IT_1203 + IT_1204 + IT_1205 + IT_1206;
    const ccomplex_t IT_1208 = IT_0384*IT_0392 + IT_0143*IT_0607 + IT_0383
      *IT_0613 + IT_0564*IT_1082 + IT_0524*IT_1103 + IT_0598*IT_1193 + IT_1195 +
       IT_1198 + IT_1199 + IT_1200 + IT_1201 + IT_1202 + IT_1203 + IT_1204 +
       IT_1205 + IT_1206;
    const ccomplex_t IT_1209 = IT_0378*IT_0762;
    const ccomplex_t IT_1210 = IT_0379*IT_0850;
    const ccomplex_t IT_1211 = IT_0496*IT_0769;
    const ccomplex_t IT_1212 = IT_0778*IT_0909;
    const ccomplex_t IT_1213 = IT_0740*IT_1043;
    const ccomplex_t IT_1214 = IT_0753*IT_0998;
    const ccomplex_t IT_1215 = IT_0732*IT_1000;
    const ccomplex_t IT_1216 = IT_0854*IT_1001;
    const ccomplex_t IT_1217 = (-2)*conj(IT_0262);
    const ccomplex_t IT_1218 = IT_0856*IT_1217;
    const ccomplex_t IT_1219 = (-2)*conj(IT_0264);
    const ccomplex_t IT_1220 = IT_0862*IT_1219;
    const ccomplex_t IT_1221 = IT_0805*IT_1005;
    const ccomplex_t IT_1222 = IT_0839*IT_1006;
    const ccomplex_t IT_1223 = (-8)*conj(IT_0268);
    const ccomplex_t IT_1224 = IT_0819*IT_1223;
    const ccomplex_t IT_1225 = IT_0816*IT_1008;
    const ccomplex_t IT_1226 = 4*conj(IT_0223);
    const ccomplex_t IT_1227 = IT_0787*IT_1226;
    const ccomplex_t IT_1228 = (-4)*conj(IT_0142);
    const ccomplex_t IT_1229 = IT_0444*IT_1228;
    const ccomplex_t IT_1230 = IT_0828*IT_1157;
    const ccomplex_t IT_1231 = (IT_0011 + IT_0015)*IT_0043 + (IT_0044 +
       IT_0045 + IT_0047 + IT_0048 + IT_0049 + IT_0051 + IT_0053 + IT_0057 +
       IT_0059 + IT_0061 + IT_0063 + IT_0066 + IT_0068 + IT_0070 + IT_0072 +
       IT_0074)*IT_0083 + (IT_0093 + IT_0103)*(IT_0144 + IT_0146 + IT_0148 +
       IT_0184 + IT_0186 + IT_0187) + (IT_0190 + IT_0193)*(IT_0197 + IT_0198 +
       IT_0202 + IT_0203) + (IT_0208 + IT_0216)*(IT_0148 + IT_0186 + IT_0224 +
       IT_0225 + IT_0231 + IT_0232) + (IT_0234 + IT_0236)*(IT_0240 + IT_0241 +
       IT_0242 + IT_0243) + (IT_0252 + IT_0261)*(IT_0263 + IT_0265 + IT_0266 +
       IT_0267 + IT_0269 + IT_0270 + IT_0271 + IT_0272 + IT_0273 + IT_0274) + 4
      *IT_0278*(IT_0317*conj(IT_0317) + 0.25*conj(IT_0335)*IT_0336) + (IT_0340 +
       IT_0344)*IT_0346 + IT_0353*(IT_0355 + IT_0357 + IT_0359) + IT_0360*(conj
      (IT_0335)*IT_0361 + IT_0335*IT_0362) + conj(IT_0196)*(IT_0367*IT_0368 +
       IT_0371*IT_0372 + IT_0376*IT_0377) + IT_0196*(IT_0367*IT_0378 + IT_0376
      *IT_0379) + (-6)*conj(IT_0201)*IT_0414 + IT_0201*(IT_0358*IT_0389 +
       IT_0367*IT_0415 + IT_0382*IT_0416 + IT_0385*IT_0417 + IT_0376*IT_0418 +
       IT_0392*IT_0419 + IT_0409*IT_0420) + (-4)*IT_0081*IT_0529 + conj(IT_0142)
      *(IT_0380*IT_0454 + IT_0368*IT_0461 + IT_0372*IT_0531 + IT_0386*IT_0536 +
       IT_0538*IT_0539) + 8*conj(IT_0399)*IT_0614 + conj(IT_0239)*(IT_0383
      *IT_0615 + IT_0386*IT_0616 + IT_0336*IT_0618 + IT_0361*IT_0620 + IT_0621
      *IT_0622 + IT_0628*IT_0629) + (-4)*conj(IT_0040)*IT_0664 + IT_0671*IT_0703
       + (-4)*conj(IT_0219)*IT_0864 + IT_0142*(IT_0415*IT_0454 + IT_0446*IT_0531
       + IT_0417*IT_0536 + IT_0865) + IT_0902*IT_0903 + IT_0694*IT_0917 + 4
      *IT_0239*(0.25*IT_0416*IT_0615 + 0.25*IT_0417*IT_0616 + conj(IT_0335)
      *IT_0618 + IT_0362*IT_0918 + IT_0912*IT_0919) + 4*conj(IT_0078)*IT_0934 +
       4*IT_0078*(conj(IT_0317)*IT_0618 + -conj(IT_0239)*IT_0920 + IT_0417
      *IT_0922 + IT_0416*IT_0923 + IT_0471*IT_0924 + IT_0420*IT_0933 + IT_0918
      *IT_0935 + IT_0919*IT_0936) + (-4)*conj(IT_0268)*IT_0991 + (-4)*IT_0268
      *IT_1009 + (-4)*conj(IT_0264)*IT_1019 + (-4)*conj(IT_0222)*IT_1028 + (-4)
      *conj(IT_0262)*IT_1036 + (-4)*IT_0262*IT_1052 + (-4)*IT_0264*IT_1053 +
       IT_1093 + 16*conj(IT_0699)*IT_1116 + 8*conj(IT_0685)*IT_1118 + 4*IT_0682*
      (IT_1119 + IT_1120 + IT_1121 + IT_1122 + IT_1123 + IT_1124 + IT_1125 +
       IT_1126 + IT_1128 + IT_1130 + IT_1132 + IT_1137 + IT_1141 + IT_1143 +
       IT_1091*IT_1145 + IT_1147) + 4*IT_0685*(IT_1119 + IT_1120 + IT_1121 +
       IT_1122 + IT_1123 + IT_1124 + IT_1125 + IT_1126 + IT_1128 + IT_1130 +
       IT_1132 + IT_1137 + IT_1141 + IT_1143 + IT_1091*IT_1146 + IT_1147) + 16
      *IT_0699*IT_1154 + 4*IT_0701*IT_1166 + 4*IT_0700*IT_1167 + 4*IT_0674
      *IT_1168 + IT_0399*(IT_0446*IT_0573 + IT_0447*IT_0583 + IT_0528*IT_0607 +
       IT_0590*IT_1007 + IT_0613*IT_1127 + IT_0551*IT_1169 + IT_1170) + IT_0412*
      (IT_0145*IT_0662 + IT_0416*IT_0877 + IT_0446*IT_0879 + IT_0528*IT_0889 +
       IT_0702*IT_0896 + IT_0900*IT_1131 + IT_0870*IT_1169 + IT_1170) + conj
      (IT_0183)*IT_1175 + IT_0183*IT_1177 + conj(IT_0230)*IT_1184 + conj(IT_0081
      )*IT_1185 + IT_0230*IT_1187 + IT_0040*IT_1188 + conj(IT_0035)*(IT_0390
      *IT_0454 + IT_0368*IT_0461 + IT_0380*IT_0516 + IT_0591*IT_0654 + IT_1189) 
      + conj(IT_0223)*(IT_0390*IT_0516 + IT_0663*IT_1179 + IT_0372*IT_1180 +
       IT_1189) + IT_0035*(IT_0418*IT_0454 + IT_0415*IT_0516 + IT_0865 + IT_1190
      ) + IT_0223*(IT_0418*IT_0516 + IT_0446*IT_1180 + IT_1190) + conj(IT_0694)
      *IT_1192 + conj(IT_0674)*IT_1196 + conj(IT_0701)*IT_1207 + conj(IT_0700)
      *IT_1208 + (IT_0219 + IT_0222)*(IT_1209 + IT_1210 + IT_1211 + IT_1212 +
       IT_1213 + IT_1214 + IT_1215 + IT_1216 + IT_1218 + IT_1220 + IT_1221 +
       IT_1222 + IT_1224 + IT_1225 + IT_1227 + IT_1229 + IT_1230);
    return create_ccomplex_return(IT_1231);
}

