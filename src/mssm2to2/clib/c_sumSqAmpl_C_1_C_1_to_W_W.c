#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_C_1_to_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_C_1_to_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + 0.5*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0010 + (-0.5)*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, 2);
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + (-0.5)*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0028*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0043 = IT_0001*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0041 + 0.5*IT_0044);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0048 = IT_0001*IT_0047;
    const ccomplex_t IT_0049 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0050 = IT_0001*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + (-0.5)*IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0046*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0059 = IT_0001*IT_0058;
    const ccomplex_t IT_0060 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0061 = IT_0001*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + 0.5*IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0066 = IT_0001*IT_0065;
    const ccomplex_t IT_0067 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0068 = IT_0001*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-0.5)*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0064*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = (-2)*IT_0021 + (-2)*IT_0039 + (-2)*IT_0057 + (
      -2)*IT_0075;
    const ccomplex_t IT_0077 = pow(m_W, -2);
    const ccomplex_t IT_0078 = s_14*s_24*IT_0077;
    const ccomplex_t IT_0079 = s_12*IT_0017;
    const ccomplex_t IT_0080 = s_13*s_23;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = IT_0077*IT_0082;
    const ccomplex_t IT_0084 = (s_12 + -IT_0078)*(s_12 + (-0.5)*IT_0083);
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0054*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0072*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0016*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0036*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 2*IT_0088 + 2*IT_0091 + 2*IT_0094 + 2*IT_0097;
    const ccomplex_t IT_0099 = 4*conj(IT_0098);
    const ccomplex_t IT_0100 = 1.125*IT_0018;
    const ccomplex_t IT_0101 = IT_0088 + IT_0091 + IT_0094 + IT_0097;
    const ccomplex_t IT_0102 = m_C_1*IT_0101;
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = cpow(IT_0071, 2);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_N_4*IT_0104;
    const ccomplex_t IT_0106 = IT_0089*IT_0105;
    const ccomplex_t IT_0107 = cpow(IT_0035, 2);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_2*IT_0107;
    const ccomplex_t IT_0109 = IT_0095*IT_0108;
    const ccomplex_t IT_0110 = IT_0106 + IT_0109;
    const ccomplex_t IT_0111 = cpow(IT_0015, 2);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_N_1*IT_0111;
    const ccomplex_t IT_0113 = IT_0092*IT_0112;
    const ccomplex_t IT_0114 = cpow(IT_0053, 2);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*m_N_3*IT_0114;
    const ccomplex_t IT_0116 = IT_0086*IT_0115;
    const ccomplex_t IT_0117 = m_C_1*IT_0094;
    const ccomplex_t IT_0118 = m_C_1*IT_0097;
    const ccomplex_t IT_0119 = m_C_1*IT_0088;
    const ccomplex_t IT_0120 = m_C_1*IT_0091;
    const ccomplex_t IT_0121 = -IT_0113 + -IT_0116 + IT_0117 + IT_0118 +
       IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = cpow(IT_0028, 2);
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*m_N_2*IT_0122;
    const ccomplex_t IT_0124 = IT_0037*IT_0123;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = m_C_1*IT_0021;
    const ccomplex_t IT_0127 = cpow(IT_0008, 2);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_N_1*IT_0127;
    const ccomplex_t IT_0129 = IT_0019*IT_0128;
    const ccomplex_t IT_0130 = IT_0126 + -IT_0129;
    const ccomplex_t IT_0131 = m_C_1*IT_0039;
    const ccomplex_t IT_0132 = m_C_1*IT_0057;
    const ccomplex_t IT_0133 = m_C_1*IT_0075;
    const ccomplex_t IT_0134 = cpow(IT_0046, 2);
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_N_3*IT_0134;
    const ccomplex_t IT_0136 = IT_0055*IT_0135;
    const ccomplex_t IT_0137 = cpow(IT_0064, 2);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*m_N_4*IT_0137;
    const ccomplex_t IT_0139 = IT_0073*IT_0138;
    const ccomplex_t IT_0140 = IT_0131 + IT_0132 + IT_0133 + -IT_0136 + 
      -IT_0139;
    const ccomplex_t IT_0141 = conj(IT_0125) + conj(IT_0130) + conj(IT_0140);
    const ccomplex_t IT_0142 = IT_0055*IT_0115;
    const ccomplex_t IT_0143 = -IT_0126 + -IT_0131 + -IT_0132 + -IT_0133 +
       IT_0142;
    const ccomplex_t IT_0144 = 16*IT_0143;
    const ccomplex_t IT_0145 = IT_0019*IT_0112;
    const ccomplex_t IT_0146 = IT_0037*IT_0108;
    const ccomplex_t IT_0147 = IT_0073*IT_0105;
    const ccomplex_t IT_0148 = -IT_0145 + -IT_0146 + -IT_0147;
    const ccomplex_t IT_0149 = (-16)*IT_0148;
    const ccomplex_t IT_0150 = IT_0125 + IT_0130 + IT_0140;
    const ccomplex_t IT_0151 = 16*conj(IT_0143);
    const ccomplex_t IT_0152 = (-16)*conj(IT_0148);
    const ccomplex_t IT_0153 = IT_0089*IT_0138;
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = IT_0095*IT_0123;
    const ccomplex_t IT_0156 = IT_0086*IT_0135;
    const ccomplex_t IT_0157 = -IT_0155 + -IT_0156;
    const ccomplex_t IT_0158 = IT_0092*IT_0128;
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = conj(IT_0157) + conj(IT_0159);
    const ccomplex_t IT_0161 = conj(IT_0154) + IT_0160;
    const ccomplex_t IT_0162 = 16*IT_0110;
    const ccomplex_t IT_0163 = (-16)*IT_0121;
    const ccomplex_t IT_0164 = IT_0157 + IT_0159;
    const ccomplex_t IT_0165 = 16*conj(IT_0110);
    const ccomplex_t IT_0166 = (-16)*conj(IT_0121);
    const ccomplex_t IT_0167 = conj(IT_0103)*((-16)*IT_0110 + 16*IT_0121) +
       IT_0103*((-16)*conj(IT_0110) + 16*conj(IT_0121)) + IT_0141*(IT_0144 +
       IT_0149) + IT_0150*(IT_0151 + IT_0152) + IT_0161*(IT_0162 + IT_0163) + 
      (IT_0154 + IT_0164)*(IT_0165 + IT_0166);
    const ccomplex_t IT_0168 = (-0.5)*IT_0083;
    const ccomplex_t IT_0169 = s_12 + IT_0168;
    const ccomplex_t IT_0170 = s_14*s_24;
    const ccomplex_t IT_0171 = (-2)*IT_0170;
    const ccomplex_t IT_0172 = IT_0079 + IT_0171;
    const ccomplex_t IT_0173 = IT_0077*IT_0172;
    const ccomplex_t IT_0174 = (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = pow(m_W, -4);
    const ccomplex_t IT_0176 = pow(m_W, 4);
    const ccomplex_t IT_0177 = s_12*IT_0176;
    const ccomplex_t IT_0178 = s_13*s_23*IT_0017;
    const ccomplex_t IT_0179 = s_14*s_24*IT_0017;
    const ccomplex_t IT_0180 = s_13*s_24*s_34;
    const ccomplex_t IT_0181 = (-2)*IT_0178 + (-2)*IT_0179 + 4*IT_0180;
    const ccomplex_t IT_0182 = IT_0177 + IT_0181;
    const ccomplex_t IT_0183 = IT_0175*IT_0182;
    const ccomplex_t IT_0184 = 0.25*IT_0183;
    const ccomplex_t IT_0185 = IT_0169 + IT_0174 + IT_0184;
    const ccomplex_t IT_0186 = 8*IT_0148;
    const ccomplex_t IT_0187 = (-8)*IT_0143;
    const ccomplex_t IT_0188 = IT_0186 + IT_0187;
    const ccomplex_t IT_0189 = 8*IT_0143;
    const ccomplex_t IT_0190 = 8*IT_0125;
    const ccomplex_t IT_0191 = 8*conj(IT_0125);
    const ccomplex_t IT_0192 = 8*IT_0130;
    const ccomplex_t IT_0193 = 8*conj(IT_0130);
    const ccomplex_t IT_0194 = (-8)*conj(IT_0143);
    const ccomplex_t IT_0195 = IT_0185*(conj(IT_0148)*IT_0188 + conj(IT_0143)
      *IT_0189 + IT_0141*IT_0190 + IT_0130*IT_0191 + (conj(IT_0130) + conj
      (IT_0140))*IT_0192 + 8*IT_0140*(conj(IT_0140) + 0.125*IT_0191 + 0.125
      *IT_0193) + IT_0148*IT_0194);
    const ccomplex_t IT_0196 = pow(s_34, 2);
    const ccomplex_t IT_0197 = (-2)*IT_0196;
    const ccomplex_t IT_0198 = IT_0176 + IT_0197;
    const ccomplex_t IT_0199 = IT_0175*IT_0198;
    const ccomplex_t IT_0200 = 0.25*IT_0199;
    const ccomplex_t IT_0201 = 1 + IT_0200;
    const ccomplex_t IT_0202 = IT_0018*IT_0201;
    const ccomplex_t IT_0203 = 8*IT_0121;
    const ccomplex_t IT_0204 = (-8)*IT_0110;
    const ccomplex_t IT_0205 = IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = 8*conj(IT_0121);
    const ccomplex_t IT_0207 = (-8)*conj(IT_0110);
    const ccomplex_t IT_0208 = IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = 8*conj(IT_0143);
    const ccomplex_t IT_0210 = 8*conj(IT_0148);
    const ccomplex_t IT_0211 = (-8)*IT_0103;
    const ccomplex_t IT_0212 = (-8)*conj(IT_0103);
    const ccomplex_t IT_0213 = s_14*s_23*s_34;
    const ccomplex_t IT_0214 = (-2)*IT_0178 + (-2)*IT_0179 + 4*IT_0213;
    const ccomplex_t IT_0215 = IT_0177 + IT_0214;
    const ccomplex_t IT_0216 = IT_0175*IT_0215;
    const ccomplex_t IT_0217 = 0.25*IT_0216;
    const ccomplex_t IT_0218 = IT_0169 + IT_0174 + IT_0217;
    const ccomplex_t IT_0219 = 8*IT_0157;
    const ccomplex_t IT_0220 = 8*IT_0154;
    const ccomplex_t IT_0221 = IT_0211 + IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = IT_0211 + IT_0219;
    const ccomplex_t IT_0223 = 8*conj(IT_0157);
    const ccomplex_t IT_0224 = IT_0212 + IT_0223;
    const ccomplex_t IT_0225 = 8*IT_0159;
    const ccomplex_t IT_0226 = 8*conj(IT_0154);
    const ccomplex_t IT_0227 = IT_0103*conj(IT_0103) + IT_0110*conj(IT_0110) +
       0.125*conj(IT_0121)*IT_0205 + 0.125*IT_0121*IT_0207 + 0.125*IT_0157
      *IT_0212 + 0.125*conj(IT_0154)*IT_0221 + 0.125*conj(IT_0157)*IT_0222 +
       0.125*IT_0154*IT_0224 + 0.125*conj(IT_0159)*(IT_0211 + IT_0219 + IT_0220 
      + IT_0225) + 0.125*IT_0159*(IT_0212 + IT_0223 + IT_0226);
    const ccomplex_t IT_0228 = s_13*s_23*IT_0077;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = s_12*IT_0196;
    const ccomplex_t IT_0231 = (-2)*IT_0178 + (-2)*IT_0179 + 2*IT_0180 + 2
      *IT_0213 + (-2)*IT_0230;
    const ccomplex_t IT_0232 = IT_0177 + IT_0231;
    const ccomplex_t IT_0233 = IT_0175*IT_0232;
    const ccomplex_t IT_0234 = (-0.25)*IT_0233;
    const ccomplex_t IT_0235 = -IT_0078;
    const ccomplex_t IT_0236 = s_12 + IT_0229 + IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = conj(IT_0121)*IT_0187 + conj(IT_0110)*IT_0189 +
       conj(IT_0154)*IT_0190 + IT_0154*IT_0191 + conj(IT_0159)*(IT_0190 +
       IT_0192) + IT_0159*(IT_0191 + IT_0193) + IT_0121*IT_0194 + conj(IT_0148)
      *IT_0205 + IT_0148*IT_0208 + IT_0110*IT_0209 + conj(IT_0130)*IT_0221 +
       conj(IT_0125)*IT_0222 + IT_0125*IT_0224 + conj(IT_0140)*(IT_0211 +
       IT_0219 + IT_0220 + IT_0225) + IT_0130*(IT_0212 + IT_0223 + IT_0226) + 8
      *IT_0140*(conj(IT_0159) + 0.125*IT_0212 + 0.125*IT_0223 + 0.125*IT_0226);
    const ccomplex_t IT_0238 = s_12 + IT_0235;
    const ccomplex_t IT_0239 = IT_0018*IT_0238;
    const ccomplex_t IT_0240 = pow(s_24, 2);
    const ccomplex_t IT_0241 = IT_0077*IT_0240;
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = IT_0018 + IT_0242;
    const ccomplex_t IT_0244 = IT_0018*IT_0243;
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = 4*IT_0098;
    const ccomplex_t IT_0247 = (s_12 + (-0.5)*IT_0083)*(IT_0018 + -IT_0241);
    const ccomplex_t IT_0248 = 2*IT_0247;
    const ccomplex_t IT_0249 = (-4)*IT_0076;
    const ccomplex_t IT_0250 = IT_0017*IT_0018;
    const ccomplex_t IT_0251 = s_23*s_24*s_34*IT_0077;
    const ccomplex_t IT_0252 = 0.5*IT_0240 + (-0.5)*IT_0250 + -IT_0251;
    const ccomplex_t IT_0253 = pow(s_23, 2);
    const ccomplex_t IT_0254 = IT_0252 + IT_0253;
    const ccomplex_t IT_0255 = IT_0077*IT_0254;
    const ccomplex_t IT_0256 = IT_0018 + IT_0242 + IT_0255;
    const ccomplex_t IT_0257 = m_C_1*IT_0256;
    const ccomplex_t IT_0258 = conj(IT_0121) + conj(IT_0154) + IT_0160;
    const ccomplex_t IT_0259 = 4*IT_0076;
    const ccomplex_t IT_0260 = IT_0121 + IT_0154 + IT_0164;
    const ccomplex_t IT_0261 = 4*conj(IT_0076);
    const ccomplex_t IT_0262 = conj(IT_0103) + conj(IT_0110);
    const ccomplex_t IT_0263 = IT_0103 + IT_0110;
    const ccomplex_t IT_0264 = (-4)*conj(IT_0076);
    const ccomplex_t IT_0265 = pow(s_14, 2);
    const ccomplex_t IT_0266 = IT_0077*IT_0265;
    const ccomplex_t IT_0267 = -IT_0018;
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = IT_0018*IT_0268;
    const ccomplex_t IT_0270 = -(s_12 + (-0.5)*IT_0083)*(IT_0018 + -IT_0266);
    const ccomplex_t IT_0271 = (-2)*IT_0270;
    const ccomplex_t IT_0272 = (-4)*IT_0098;
    const ccomplex_t IT_0273 = s_13*s_14*s_34*IT_0077;
    const ccomplex_t IT_0274 = (-0.5)*IT_0250 + 0.5*IT_0265 + -IT_0273;
    const ccomplex_t IT_0275 = pow(s_13, 2);
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = IT_0077*IT_0276;
    const ccomplex_t IT_0278 = -IT_0266;
    const ccomplex_t IT_0279 = IT_0018 + IT_0277 + IT_0278;
    const ccomplex_t IT_0280 = m_C_1*IT_0279;
    const ccomplex_t IT_0281 = conj(IT_0125) + conj(IT_0130) + conj(IT_0140) +
       conj(IT_0148);
    const ccomplex_t IT_0282 = IT_0125 + IT_0130 + IT_0140 + IT_0148;
    const ccomplex_t IT_0283 = (-4)*conj(IT_0098);
    const ccomplex_t IT_0284 = s_13*s_24*s_34*IT_0077;
    const ccomplex_t IT_0285 = (-2)*IT_0080 + -IT_0170 + 2*IT_0284;
    const ccomplex_t IT_0286 = IT_0079 + IT_0285;
    const ccomplex_t IT_0287 = IT_0077*IT_0286;
    const ccomplex_t IT_0288 = -IT_0078 + (-0.5)*IT_0287;
    const ccomplex_t IT_0289 = s_12 + IT_0288;
    const ccomplex_t IT_0290 = m_C_1*IT_0289;
    const ccomplex_t IT_0291 = (-2)*IT_0290;
    const ccomplex_t IT_0292 = 2*IT_0076;
    const ccomplex_t IT_0293 = 2*conj(IT_0076);
    const ccomplex_t IT_0294 = (-2)*IT_0076;
    const ccomplex_t IT_0295 = (-2)*conj(IT_0076);
    const ccomplex_t IT_0296 = -IT_0250;
    const ccomplex_t IT_0297 = IT_0240 + IT_0296;
    const ccomplex_t IT_0298 = IT_0077*IT_0297;
    const ccomplex_t IT_0299 = IT_0241 + IT_0267;
    const ccomplex_t IT_0300 = (-4)*IT_0299;
    const ccomplex_t IT_0301 = IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = m_C_1*IT_0301;
    const ccomplex_t IT_0303 = -IT_0170;
    const ccomplex_t IT_0304 = IT_0079 + IT_0303;
    const ccomplex_t IT_0305 = IT_0077*IT_0304;
    const ccomplex_t IT_0306 = -IT_0078 + (-0.25)*IT_0305;
    const ccomplex_t IT_0307 = s_12 + IT_0306;
    const ccomplex_t IT_0308 = m_C_1*IT_0307;
    const ccomplex_t IT_0309 = (-4)*IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0098;
    const ccomplex_t IT_0311 = 2*conj(IT_0098);
    const ccomplex_t IT_0312 = (-2)*IT_0098;
    const ccomplex_t IT_0313 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0314 = s_14*s_23*s_34*IT_0077;
    const ccomplex_t IT_0315 = (-2)*IT_0080 + -IT_0170 + 2*IT_0314;
    const ccomplex_t IT_0316 = IT_0079 + IT_0315;
    const ccomplex_t IT_0317 = IT_0077*IT_0316;
    const ccomplex_t IT_0318 = -IT_0078 + (-0.5)*IT_0317;
    const ccomplex_t IT_0319 = s_12 + IT_0318;
    const ccomplex_t IT_0320 = m_C_1*IT_0319;
    const ccomplex_t IT_0321 = (-2)*IT_0320;
    const ccomplex_t IT_0322 = IT_0265 + IT_0296;
    const ccomplex_t IT_0323 = IT_0077*IT_0322;
    const ccomplex_t IT_0324 = (-4)*IT_0268;
    const ccomplex_t IT_0325 = IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = m_C_1*IT_0325;
    const ccomplex_t IT_0327 = IT_0021 + IT_0039 + IT_0057 + IT_0075;
    const ccomplex_t IT_0328 = 9*IT_0250;
    const ccomplex_t IT_0329 = s_12*s_14*IT_0017;
    const ccomplex_t IT_0330 = s_24*IT_0275;
    const ccomplex_t IT_0331 = s_12*s_13*s_34;
    const ccomplex_t IT_0332 = s_24*IT_0017*IT_0018;
    const ccomplex_t IT_0333 = 2*IT_0330 + (-2)*IT_0331 + -IT_0332;
    const ccomplex_t IT_0334 = IT_0329 + IT_0333;
    const ccomplex_t IT_0335 = IT_0077*IT_0334;
    const ccomplex_t IT_0336 = s_24*IT_0018;
    const ccomplex_t IT_0337 = 2*IT_0336;
    const ccomplex_t IT_0338 = s_12*s_14;
    const ccomplex_t IT_0339 = (-2)*IT_0338;
    const ccomplex_t IT_0340 = IT_0335 + IT_0337 + IT_0339;
    const ccomplex_t IT_0341 = IT_0077*IT_0182;
    const ccomplex_t IT_0342 = pow(m_W, 6);
    const ccomplex_t IT_0343 = s_12*IT_0342;
    const ccomplex_t IT_0344 = s_13*s_23*IT_0176;
    const ccomplex_t IT_0345 = (-2)*IT_0344;
    const ccomplex_t IT_0346 = IT_0343 + IT_0345;
    const ccomplex_t IT_0347 = IT_0175*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0347;
    const ccomplex_t IT_0349 = (-4)*IT_0170;
    const ccomplex_t IT_0350 = (-8)*IT_0327;
    const ccomplex_t IT_0351 = s_12*s_24*IT_0017;
    const ccomplex_t IT_0352 = s_13*s_23*s_24;
    const ccomplex_t IT_0353 = s_14*IT_0017*IT_0018;
    const ccomplex_t IT_0354 = s_13*s_34*IT_0018;
    const ccomplex_t IT_0355 = (-2)*IT_0352 + -IT_0353 + 2*IT_0354;
    const ccomplex_t IT_0356 = IT_0351 + IT_0355;
    const ccomplex_t IT_0357 = IT_0077*IT_0356;
    const ccomplex_t IT_0358 = s_14*IT_0018;
    const ccomplex_t IT_0359 = 2*IT_0358;
    const ccomplex_t IT_0360 = s_12*s_24;
    const ccomplex_t IT_0361 = (-2)*IT_0360;
    const ccomplex_t IT_0362 = IT_0357 + IT_0359 + IT_0361;
    const ccomplex_t IT_0363 = s_24*IT_0176;
    const ccomplex_t IT_0364 = IT_0175*IT_0363;
    const ccomplex_t IT_0365 = s_24*IT_0017;
    const ccomplex_t IT_0366 = IT_0077*IT_0365;
    const ccomplex_t IT_0367 = (-0.166666666666667)*IT_0364 +
       0.666666666666667*IT_0366;
    const ccomplex_t IT_0368 = s_24 + IT_0367;
    const ccomplex_t IT_0369 = m_C_1*IT_0368;
    const ccomplex_t IT_0370 = (-6)*IT_0369;
    const ccomplex_t IT_0371 = (-2)*IT_0148;
    const ccomplex_t IT_0372 = s_14*IT_0176;
    const ccomplex_t IT_0373 = s_13*s_34*IT_0017;
    const ccomplex_t IT_0374 = (-2)*IT_0373;
    const ccomplex_t IT_0375 = IT_0372 + IT_0374;
    const ccomplex_t IT_0376 = IT_0175*IT_0375;
    const ccomplex_t IT_0377 = s_14*IT_0017;
    const ccomplex_t IT_0378 = s_13*s_34;
    const ccomplex_t IT_0379 = (-2)*IT_0378;
    const ccomplex_t IT_0380 = IT_0377 + IT_0379;
    const ccomplex_t IT_0381 = IT_0077*IT_0380;
    const ccomplex_t IT_0382 = 0.166666666666667*IT_0376 + (-0.666666666666667
      )*IT_0381;
    const ccomplex_t IT_0383 = s_14 + IT_0382;
    const ccomplex_t IT_0384 = m_C_1*IT_0383;
    const ccomplex_t IT_0385 = 6*IT_0384;
    const ccomplex_t IT_0386 = 2*IT_0143;
    const ccomplex_t IT_0387 = IT_0371 + IT_0386;
    const ccomplex_t IT_0388 = (-2)*IT_0125;
    const ccomplex_t IT_0389 = (-2)*IT_0130;
    const ccomplex_t IT_0390 = (-2)*IT_0140;
    const ccomplex_t IT_0391 = (-10)*s_24;
    const ccomplex_t IT_0392 = IT_0364 + IT_0391;
    const ccomplex_t IT_0393 = m_C_1*IT_0392;
    const ccomplex_t IT_0394 = m_C_1*IT_0077;
    const ccomplex_t IT_0395 = IT_0077*IT_0372;
    const ccomplex_t IT_0396 = (-2)*IT_0378 + (-0.5)*IT_0395;
    const ccomplex_t IT_0397 = IT_0377 + IT_0396;
    const ccomplex_t IT_0398 = IT_0394*IT_0397;
    const ccomplex_t IT_0399 = (-2)*IT_0398;
    const ccomplex_t IT_0400 = s_23*s_34*IT_0077;
    const ccomplex_t IT_0401 = s_23*s_34*IT_0017;
    const ccomplex_t IT_0402 = (-2)*IT_0401;
    const ccomplex_t IT_0403 = IT_0363 + IT_0402;
    const ccomplex_t IT_0404 = IT_0175*IT_0403;
    const ccomplex_t IT_0405 = (-0.666666666666667)*IT_0400 + (
      -0.166666666666667)*IT_0404;
    const ccomplex_t IT_0406 = s_24 + IT_0405;
    const ccomplex_t IT_0407 = m_C_1*IT_0406;
    const ccomplex_t IT_0408 = (-6)*IT_0407;
    const ccomplex_t IT_0409 = IT_0399 + IT_0408;
    const ccomplex_t IT_0410 = 2*IT_0157;
    const ccomplex_t IT_0411 = 2*IT_0159;
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = 2*IT_0121;
    const ccomplex_t IT_0414 = 2*IT_0154;
    const ccomplex_t IT_0415 = (-2)*IT_0110;
    const ccomplex_t IT_0416 = (-2)*IT_0103;
    const ccomplex_t IT_0417 = IT_0327*IT_0328 + (-0.25)*IT_0246*IT_0340 + (
      -0.25)*(IT_0079 + -IT_0341 + IT_0348 + IT_0349)*IT_0350 + (-0.25)*IT_0249
      *IT_0362 + (-0.25)*IT_0370*IT_0371 + (-0.25)*IT_0385*(IT_0387 + IT_0388 +
       IT_0389 + IT_0390) + (-0.25)*(IT_0386 + IT_0388 + IT_0389 + IT_0390)
      *IT_0393 + (-0.25)*IT_0409*(IT_0412 + IT_0413 + IT_0414 + IT_0415 +
       IT_0416);
    const ccomplex_t IT_0418 = (-2)*conj(IT_0148);
    const ccomplex_t IT_0419 = 2*conj(IT_0143);
    const ccomplex_t IT_0420 = IT_0418 + IT_0419;
    const ccomplex_t IT_0421 = (-2)*conj(IT_0125);
    const ccomplex_t IT_0422 = (-2)*conj(IT_0130);
    const ccomplex_t IT_0423 = (-2)*conj(IT_0140);
    const ccomplex_t IT_0424 = 2*conj(IT_0157);
    const ccomplex_t IT_0425 = 2*conj(IT_0159);
    const ccomplex_t IT_0426 = IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = 2*conj(IT_0121);
    const ccomplex_t IT_0428 = 2*conj(IT_0154);
    const ccomplex_t IT_0429 = (-2)*conj(IT_0110);
    const ccomplex_t IT_0430 = (-2)*conj(IT_0103);
    const ccomplex_t IT_0431 = IT_0099*IT_0340 + IT_0264*IT_0362 + IT_0370
      *IT_0418 + IT_0393*(IT_0419 + IT_0421 + IT_0422 + IT_0423) + IT_0385*
      (IT_0420 + IT_0421 + IT_0422 + IT_0423) + IT_0409*(IT_0426 + IT_0427 +
       IT_0428 + IT_0429 + IT_0430);
    const ccomplex_t IT_0432 = -IT_0088 + -IT_0091 + -IT_0094 + -IT_0097;
    const ccomplex_t IT_0433 = IT_0077*IT_0215;
    const ccomplex_t IT_0434 = -IT_0433;
    const ccomplex_t IT_0435 = IT_0079 + IT_0348 + IT_0349 + IT_0434;
    const ccomplex_t IT_0436 = s_13*s_34*IT_0077;
    const ccomplex_t IT_0437 = (-0.166666666666667)*IT_0376 + (
      -0.666666666666667)*IT_0436;
    const ccomplex_t IT_0438 = s_14 + IT_0437;
    const ccomplex_t IT_0439 = m_C_1*IT_0438;
    const ccomplex_t IT_0440 = (-6)*IT_0439;
    const ccomplex_t IT_0441 = IT_0388 + IT_0389 + IT_0390;
    const ccomplex_t IT_0442 = IT_0077*IT_0363;
    const ccomplex_t IT_0443 = s_23*s_34;
    const ccomplex_t IT_0444 = (-2)*IT_0365 + 4*IT_0443;
    const ccomplex_t IT_0445 = IT_0442 + IT_0444;
    const ccomplex_t IT_0446 = IT_0394*IT_0445;
    const ccomplex_t IT_0447 = (-0.5)*IT_0446;
    const ccomplex_t IT_0448 = (-2)*IT_0443;
    const ccomplex_t IT_0449 = IT_0365 + IT_0448;
    const ccomplex_t IT_0450 = IT_0077*IT_0449;
    const ccomplex_t IT_0451 = 0.166666666666667*IT_0404 + (-0.666666666666667
      )*IT_0450;
    const ccomplex_t IT_0452 = s_24 + IT_0451;
    const ccomplex_t IT_0453 = m_C_1*IT_0452;
    const ccomplex_t IT_0454 = 6*IT_0453;
    const ccomplex_t IT_0455 = (-0.5)*IT_0454;
    const ccomplex_t IT_0456 = IT_0412 + IT_0414 + IT_0416;
    const ccomplex_t IT_0457 = IT_0175*IT_0372;
    const ccomplex_t IT_0458 = (-10)*s_14;
    const ccomplex_t IT_0459 = IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = m_C_1*IT_0459;
    const ccomplex_t IT_0461 = (-0.5)*IT_0460;
    const ccomplex_t IT_0462 = IT_0018*IT_0077;
    const ccomplex_t IT_0463 = (-0.25)*IT_0176;
    const ccomplex_t IT_0464 = IT_0196 + IT_0463;
    const ccomplex_t IT_0465 = IT_0462*IT_0464;
    const ccomplex_t IT_0466 = (-4)*IT_0465;
    const ccomplex_t IT_0467 = (-2)*IT_0327;
    const ccomplex_t IT_0468 = (-0.5)*IT_0466*IT_0467;
    const ccomplex_t IT_0469 = 0.25*IT_0347;
    const ccomplex_t IT_0470 = IT_0079 + IT_0469;
    const ccomplex_t IT_0471 = (-0.5)*IT_0350*IT_0470;
    const ccomplex_t IT_0472 = s_13*s_14*s_23;
    const ccomplex_t IT_0473 = s_23*s_34*IT_0018;
    const ccomplex_t IT_0474 = -IT_0332 + (-2)*IT_0472 + 2*IT_0473;
    const ccomplex_t IT_0475 = IT_0329 + IT_0474;
    const ccomplex_t IT_0476 = IT_0077*IT_0475;
    const ccomplex_t IT_0477 = (-0.5)*IT_0476;
    const ccomplex_t IT_0478 = -IT_0336;
    const ccomplex_t IT_0479 = IT_0338 + IT_0477 + IT_0478;
    const ccomplex_t IT_0480 = (-0.5)*IT_0246*IT_0479;
    const ccomplex_t IT_0481 = s_14*IT_0253;
    const ccomplex_t IT_0482 = s_12*s_23*s_34;
    const ccomplex_t IT_0483 = -IT_0353 + 2*IT_0481 + (-2)*IT_0482;
    const ccomplex_t IT_0484 = IT_0351 + IT_0483;
    const ccomplex_t IT_0485 = IT_0077*IT_0484;
    const ccomplex_t IT_0486 = 0.5*IT_0485;
    const ccomplex_t IT_0487 = -IT_0486;
    const ccomplex_t IT_0488 = -IT_0358;
    const ccomplex_t IT_0489 = IT_0360 + IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = (-0.5)*IT_0249*IT_0489;
    const ccomplex_t IT_0491 = (-2)*conj(IT_0432)*(IT_0328*IT_0432 + 2*IT_0432
      *IT_0435 + (-0.5)*IT_0440*IT_0441 + IT_0387*IT_0447 + (IT_0413 + IT_0415)
      *IT_0455 + IT_0456*IT_0461 + IT_0468 + IT_0471 + IT_0480 + IT_0490);
    const ccomplex_t IT_0492 = IT_0077*IT_0377;
    const ccomplex_t IT_0493 = (-0.166666666666667)*IT_0457 +
       0.666666666666667*IT_0492;
    const ccomplex_t IT_0494 = s_14 + IT_0493;
    const ccomplex_t IT_0495 = m_C_1*IT_0494;
    const ccomplex_t IT_0496 = (-6)*IT_0495;
    const ccomplex_t IT_0497 = IT_0328*IT_0432 + 2*IT_0432*IT_0435 + (-0.5)
      *IT_0387*IT_0440 + IT_0441*IT_0447 + IT_0455*IT_0456 + IT_0413*IT_0461 +
       IT_0468 + IT_0471 + IT_0480 + IT_0490 + (-0.5)*IT_0415*IT_0496;
    const ccomplex_t IT_0498 = IT_0426 + IT_0428 + IT_0430;
    const ccomplex_t IT_0499 = IT_0421 + IT_0422 + IT_0423;
    const ccomplex_t IT_0500 = (-2)*conj(IT_0327);
    const ccomplex_t IT_0501 = IT_0466*IT_0500;
    const ccomplex_t IT_0502 = (-8)*conj(IT_0327);
    const ccomplex_t IT_0503 = IT_0470*IT_0502;
    const ccomplex_t IT_0504 = IT_0099*IT_0479;
    const ccomplex_t IT_0505 = IT_0264*IT_0489;
    const ccomplex_t IT_0506 = IT_0076*IT_0085*IT_0099 + IT_0100*IT_0167 +
       IT_0195 + IT_0202*(IT_0161*IT_0188 + conj(IT_0103)*IT_0189 + IT_0141
      *IT_0205 + IT_0150*IT_0208 + IT_0103*IT_0209 + (IT_0154 + IT_0164)*
      (IT_0194 + IT_0210) + conj(IT_0148)*IT_0211 + IT_0148*IT_0212) + 8*IT_0218
      *IT_0227 + IT_0236*IT_0237 + conj(IT_0076)*((-12)*IT_0098*IT_0239 + 12
      *IT_0076*IT_0245 + IT_0085*IT_0246 + IT_0248*IT_0249) + IT_0257*(IT_0258
      *IT_0259 + IT_0260*IT_0261 + IT_0249*IT_0262 + IT_0263*IT_0264) + conj
      (IT_0098)*((-12)*IT_0076*IT_0239 + 12*IT_0098*IT_0269 + IT_0271*IT_0272) +
       IT_0280*(IT_0099*IT_0143 + conj(IT_0143)*IT_0246 + IT_0272*IT_0281 +
       IT_0282*IT_0283) + IT_0291*(IT_0281*IT_0292 + IT_0282*IT_0293 + conj
      (IT_0143)*IT_0294 + IT_0143*IT_0295) + (conj(IT_0143)*IT_0292 + IT_0143
      *IT_0293 + IT_0281*IT_0294 + IT_0282*IT_0295)*IT_0302 + IT_0309*(IT_0262
      *IT_0292 + IT_0263*IT_0293 + IT_0258*IT_0294 + IT_0260*IT_0295 + IT_0281
      *IT_0310 + IT_0282*IT_0311 + conj(IT_0143)*IT_0312 + IT_0143*IT_0313) + 
      (IT_0262*IT_0310 + IT_0263*IT_0311 + IT_0258*IT_0312 + IT_0260*IT_0313)
      *IT_0321 + (IT_0258*IT_0310 + IT_0260*IT_0311 + IT_0262*IT_0312 + IT_0263
      *IT_0313)*IT_0326 + (-4)*conj(IT_0327)*IT_0417 + IT_0327*IT_0431 + IT_0491
       + (-2)*conj(IT_0432)*IT_0497 + IT_0432*(IT_0420*IT_0446 + (IT_0427 +
       IT_0429)*IT_0454 + IT_0460*IT_0498 + IT_0440*IT_0499 + IT_0501 + IT_0503 
      + IT_0504 + IT_0505) + IT_0432*(IT_0420*IT_0440 + IT_0427*IT_0460 +
       IT_0429*IT_0496 + IT_0454*IT_0498 + IT_0446*IT_0499 + IT_0501 + IT_0503 +
       IT_0504 + IT_0505);
    return create_ccomplex_return(IT_0506);
}

