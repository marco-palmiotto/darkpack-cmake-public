#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_H0_Z(
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
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
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
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = -IT_0000;
    const ccomplex_t IT_0005 = IT_0003 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0014 = sin(alpha);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0008*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0023 = IT_0014*IT_0019;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0012 + -IT_0016 + 
      -IT_0021 + IT_0024);
    const ccomplex_t IT_0026 = IT_0010*IT_0018;
    const ccomplex_t IT_0027 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = IT_0009*IT_0019;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0029*IT_0031;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0028 + IT_0030 + 
      -IT_0032 + -IT_0033);
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = pow(m_H0, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0037 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = cos(beta);
    const ccomplex_t IT_0042 = IT_0010*IT_0041;
    const ccomplex_t IT_0043 = IT_0013*IT_0042;
    const ccomplex_t IT_0044 = sin(beta);
    const ccomplex_t IT_0045 = IT_0010*IT_0044;
    const ccomplex_t IT_0046 = IT_0007*IT_0045;
    const ccomplex_t IT_0047 = IT_0019*IT_0041;
    const ccomplex_t IT_0048 = IT_0022*IT_0047;
    const ccomplex_t IT_0049 = IT_0019*IT_0044;
    const ccomplex_t IT_0050 = IT_0017*IT_0049;
    const ccomplex_t IT_0051 = IT_0043 + -IT_0046 + -IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = -e_em*IT_0010*IT_0019*(IT_0014*IT_0041 + 
      -IT_0008*IT_0044);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0058 = IT_0011*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0060 = IT_0023*IT_0059;
    const ccomplex_t IT_0061 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0062 = IT_0011*IT_0061;
    const ccomplex_t IT_0063 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0064 = IT_0023*IT_0063;
    const ccomplex_t IT_0065 = IT_0058 + IT_0060 + IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0067 = IT_0015*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0069 = IT_0020*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0071 = IT_0015*IT_0070;
    const ccomplex_t IT_0072 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0073 = IT_0020*IT_0072;
    const ccomplex_t IT_0074 = -IT_0067 + -IT_0069 + -IT_0071 + -IT_0073;
    const ccomplex_t IT_0075 = IT_0065 + IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0079 = IT_0026*IT_0078;
    const ccomplex_t IT_0080 = IT_0029*IT_0078;
    const ccomplex_t IT_0081 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0082 = IT_0026*IT_0081;
    const ccomplex_t IT_0083 = IT_0029*IT_0081;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0079 + IT_0080 + 
      -IT_0082 + -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0077*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0091 = IT_0011*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0093 = IT_0023*IT_0092;
    const ccomplex_t IT_0094 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0095 = IT_0011*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0097 = IT_0023*IT_0096;
    const ccomplex_t IT_0098 = IT_0091 + IT_0093 + IT_0095 + IT_0097;
    const ccomplex_t IT_0099 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0100 = IT_0015*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0102 = IT_0020*IT_0101;
    const ccomplex_t IT_0103 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0104 = IT_0015*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0106 = IT_0020*IT_0105;
    const ccomplex_t IT_0107 = -IT_0100 + -IT_0102 + -IT_0104 + -IT_0106;
    const ccomplex_t IT_0108 = IT_0098 + IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0112 = IT_0026*IT_0111;
    const ccomplex_t IT_0113 = IT_0029*IT_0111;
    const ccomplex_t IT_0114 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0115 = IT_0026*IT_0114;
    const ccomplex_t IT_0116 = IT_0029*IT_0114;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*(IT_0112 + IT_0113 + 
      -IT_0115 + -IT_0116);
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = IT_0110*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = IT_0119*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0124 = IT_0011*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0126 = IT_0023*IT_0125;
    const ccomplex_t IT_0127 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0128 = IT_0011*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0130 = IT_0023*IT_0129;
    const ccomplex_t IT_0131 = IT_0124 + IT_0126 + IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0133 = IT_0015*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0135 = IT_0020*IT_0134;
    const ccomplex_t IT_0136 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0137 = IT_0015*IT_0136;
    const ccomplex_t IT_0138 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0139 = IT_0020*IT_0138;
    const ccomplex_t IT_0140 = -IT_0133 + -IT_0135 + -IT_0137 + -IT_0139;
    const ccomplex_t IT_0141 = IT_0131 + IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0145 = IT_0026*IT_0144;
    const ccomplex_t IT_0146 = IT_0029*IT_0144;
    const ccomplex_t IT_0147 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0148 = IT_0026*IT_0147;
    const ccomplex_t IT_0149 = IT_0029*IT_0147;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0145 + IT_0146 + 
      -IT_0148 + -IT_0149);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = IT_0143*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0152*IT_0154;
    const ccomplex_t IT_0156 = IT_0007*IT_0042;
    const ccomplex_t IT_0157 = IT_0013*IT_0045;
    const ccomplex_t IT_0158 = IT_0017*IT_0047;
    const ccomplex_t IT_0159 = IT_0022*IT_0049;
    const ccomplex_t IT_0160 = IT_0156 + IT_0157 + -IT_0158 + -IT_0159;
    const ccomplex_t IT_0161 = e_em*IT_0010*IT_0019*(IT_0008*IT_0041 + IT_0014
      *IT_0044);
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = pow(m_Z, 2);
    const ccomplex_t IT_0165 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0164 + 
      -reg_prop, -1);
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = IT_0163*IT_0166;
    const ccomplex_t IT_0168 = -IT_0040 + (-0.5)*IT_0056 + (-2)*IT_0089 + (-2)
      *IT_0122 + (-2)*IT_0155 + (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0170 = IT_0042*IT_0169;
    const ccomplex_t IT_0171 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0172 = IT_0045*IT_0171;
    const ccomplex_t IT_0173 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0174 = IT_0047*IT_0173;
    const ccomplex_t IT_0175 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0176 = IT_0049*IT_0175;
    const ccomplex_t IT_0177 = -IT_0170 + IT_0172 + IT_0174 + -IT_0176;
    const ccomplex_t IT_0178 = IT_0052*IT_0177;
    const ccomplex_t IT_0179 = IT_0055*IT_0178;
    const ccomplex_t IT_0180 = IT_0011*IT_0171;
    const ccomplex_t IT_0181 = IT_0015*IT_0169;
    const ccomplex_t IT_0182 = IT_0020*IT_0175;
    const ccomplex_t IT_0183 = IT_0023*IT_0173;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*(IT_0180 + -IT_0181 + 
      -IT_0182 + IT_0183);
    const ccomplex_t IT_0185 = IT_0035*IT_0184;
    const ccomplex_t IT_0186 = IT_0039*IT_0185;
    const ccomplex_t IT_0187 = IT_0042*IT_0171;
    const ccomplex_t IT_0188 = IT_0045*IT_0169;
    const ccomplex_t IT_0189 = IT_0047*IT_0175;
    const ccomplex_t IT_0190 = IT_0049*IT_0173;
    const ccomplex_t IT_0191 = -IT_0187 + -IT_0188 + IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = IT_0162*IT_0191;
    const ccomplex_t IT_0193 = IT_0166*IT_0192;
    const ccomplex_t IT_0194 = 0.5*IT_0179 + -IT_0186 + 0.5*IT_0193;
    const ccomplex_t IT_0195 = 2*IT_0194;
    const ccomplex_t IT_0196 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0197 = IT_0011*IT_0196;
    const ccomplex_t IT_0198 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0199 = IT_0023*IT_0198;
    const ccomplex_t IT_0200 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0201 = IT_0011*IT_0200;
    const ccomplex_t IT_0202 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0203 = IT_0023*IT_0202;
    const ccomplex_t IT_0204 = IT_0197 + IT_0199 + IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0206 = IT_0015*IT_0205;
    const ccomplex_t IT_0207 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0208 = IT_0020*IT_0207;
    const ccomplex_t IT_0209 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0210 = IT_0015*IT_0209;
    const ccomplex_t IT_0211 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0212 = IT_0020*IT_0211;
    const ccomplex_t IT_0213 = -IT_0206 + -IT_0208 + -IT_0210 + -IT_0212;
    const ccomplex_t IT_0214 = IT_0204 + IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0218 = IT_0026*IT_0217;
    const ccomplex_t IT_0219 = IT_0029*IT_0217;
    const ccomplex_t IT_0220 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0221 = IT_0026*IT_0220;
    const ccomplex_t IT_0222 = IT_0029*IT_0220;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*(IT_0218 + IT_0219 + 
      -IT_0221 + -IT_0222);
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = IT_0216*IT_0224;
    const ccomplex_t IT_0226 = IT_0154*IT_0225;
    const ccomplex_t IT_0227 = (-0.5)*IT_0034;
    const ccomplex_t IT_0228 = IT_0184*IT_0227;
    const ccomplex_t IT_0229 = IT_0039*IT_0228;
    const ccomplex_t IT_0230 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0231 = IT_0011*IT_0230;
    const ccomplex_t IT_0232 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0233 = IT_0023*IT_0232;
    const ccomplex_t IT_0234 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0235 = IT_0011*IT_0234;
    const ccomplex_t IT_0236 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0237 = IT_0023*IT_0236;
    const ccomplex_t IT_0238 = IT_0231 + IT_0233 + IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0240 = IT_0015*IT_0239;
    const ccomplex_t IT_0241 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0242 = IT_0020*IT_0241;
    const ccomplex_t IT_0243 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0244 = IT_0015*IT_0243;
    const ccomplex_t IT_0245 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0246 = IT_0020*IT_0245;
    const ccomplex_t IT_0247 = -IT_0240 + -IT_0242 + -IT_0244 + -IT_0246;
    const ccomplex_t IT_0248 = IT_0238 + IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = 0.5*IT_0249;
    const ccomplex_t IT_0251 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0252 = IT_0026*IT_0251;
    const ccomplex_t IT_0253 = IT_0029*IT_0251;
    const ccomplex_t IT_0254 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0255 = IT_0026*IT_0254;
    const ccomplex_t IT_0256 = IT_0029*IT_0254;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*(IT_0252 + IT_0253 + 
      -IT_0255 + -IT_0256);
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0250*IT_0258;
    const ccomplex_t IT_0260 = IT_0088*IT_0259;
    const ccomplex_t IT_0261 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0262 = IT_0011*IT_0261;
    const ccomplex_t IT_0263 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0264 = IT_0023*IT_0263;
    const ccomplex_t IT_0265 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0266 = IT_0011*IT_0265;
    const ccomplex_t IT_0267 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0268 = IT_0023*IT_0267;
    const ccomplex_t IT_0269 = IT_0262 + IT_0264 + IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0271 = IT_0015*IT_0270;
    const ccomplex_t IT_0272 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0273 = IT_0020*IT_0272;
    const ccomplex_t IT_0274 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0275 = IT_0015*IT_0274;
    const ccomplex_t IT_0276 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0277 = IT_0020*IT_0276;
    const ccomplex_t IT_0278 = -IT_0271 + -IT_0273 + -IT_0275 + -IT_0277;
    const ccomplex_t IT_0279 = IT_0269 + IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = 0.5*IT_0280;
    const ccomplex_t IT_0282 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0283 = IT_0026*IT_0282;
    const ccomplex_t IT_0284 = IT_0029*IT_0282;
    const ccomplex_t IT_0285 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0286 = IT_0026*IT_0285;
    const ccomplex_t IT_0287 = IT_0029*IT_0285;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*(IT_0283 + IT_0284 + 
      -IT_0286 + -IT_0287);
    const ccomplex_t IT_0289 = (-0.5)*IT_0288;
    const ccomplex_t IT_0290 = IT_0281*IT_0289;
    const ccomplex_t IT_0291 = IT_0121*IT_0290;
    const ccomplex_t IT_0292 = (-0.5)*IT_0179 + (-0.5)*IT_0193 + (-2)*IT_0226 
      + -IT_0229 + (-2)*IT_0260 + (-2)*IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0292;
    const ccomplex_t IT_0294 = IT_0195 + IT_0293;
    const ccomplex_t IT_0295 = 2*conj(IT_0194);
    const ccomplex_t IT_0296 = (-2)*conj(IT_0292);
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = IT_0025*IT_0227;
    const ccomplex_t IT_0299 = IT_0039*IT_0298;
    const ccomplex_t IT_0300 = 0.5*IT_0056 + 0.5*IT_0167 + -IT_0299;
    const ccomplex_t IT_0301 = 2*IT_0292;
    const ccomplex_t IT_0302 = (-2)*IT_0194;
    const ccomplex_t IT_0303 = IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = 2*conj(IT_0292);
    const ccomplex_t IT_0305 = (-2)*conj(IT_0194);
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = s_12 + IT_0308;
    const ccomplex_t IT_0310 = s_12*IT_0309;
    const ccomplex_t IT_0311 = cpow((-2)*s_13 + IT_0037 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*IT_0311;
    const ccomplex_t IT_0313 = IT_0185*IT_0312;
    const ccomplex_t IT_0314 = (-0.5)*IT_0179 + (-0.5)*IT_0193 + IT_0313;
    const ccomplex_t IT_0315 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0316 = (0 + _Complex_I*1)*IT_0315;
    const ccomplex_t IT_0317 = IT_0259*IT_0316;
    const ccomplex_t IT_0318 = IT_0228*IT_0312;
    const ccomplex_t IT_0319 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = IT_0290*IT_0320;
    const ccomplex_t IT_0322 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*IT_0322;
    const ccomplex_t IT_0324 = IT_0225*IT_0323;
    const ccomplex_t IT_0325 = 0.5*IT_0179 + 0.5*IT_0193 + 2*IT_0317 + IT_0318
       + 2*IT_0321 + 2*IT_0324;
    const ccomplex_t IT_0326 = IT_0298*IT_0312;
    const ccomplex_t IT_0327 = (-0.5)*IT_0056 + (-0.5)*IT_0167 + IT_0326;
    const ccomplex_t IT_0328 = 2*IT_0300;
    const ccomplex_t IT_0329 = (-2)*IT_0168;
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = 2*conj(IT_0300);
    const ccomplex_t IT_0332 = (-2)*conj(IT_0168);
    const ccomplex_t IT_0333 = IT_0331 + IT_0332;
    const ccomplex_t IT_0334 = IT_0086*IT_0316;
    const ccomplex_t IT_0335 = IT_0119*IT_0320;
    const ccomplex_t IT_0336 = IT_0152*IT_0323;
    const ccomplex_t IT_0337 = IT_0036*IT_0312;
    const ccomplex_t IT_0338 = 0.5*IT_0056 + 0.5*IT_0167 + 2*IT_0334 + 2
      *IT_0335 + 2*IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = 2*IT_0168;
    const ccomplex_t IT_0340 = (-2)*IT_0300;
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = 2*conj(IT_0168);
    const ccomplex_t IT_0343 = (-2)*conj(IT_0300);
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = s_12*s_14;
    const ccomplex_t IT_0346 = s_24*IT_0000;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = IT_0345 + IT_0347;
    const ccomplex_t IT_0349 = 0.5*IT_0186;
    const ccomplex_t IT_0350 = (-0.5)*IT_0313;
    const ccomplex_t IT_0351 = -IT_0317 + (-0.5)*IT_0318 + -IT_0321 + -IT_0324;
    const ccomplex_t IT_0352 = IT_0226 + 0.5*IT_0229 + IT_0260 + IT_0291;
    const ccomplex_t IT_0353 = conj(IT_0351) + conj(IT_0352);
    const ccomplex_t IT_0354 = IT_0351 + IT_0352;
    const ccomplex_t IT_0355 = 0.5*IT_0299;
    const ccomplex_t IT_0356 = (-0.5)*IT_0326;
    const ccomplex_t IT_0357 = 0.5*IT_0040;
    const ccomplex_t IT_0358 = -IT_0334 + -IT_0335 + -IT_0336 + (-0.5)*IT_0337;
    const ccomplex_t IT_0359 = IT_0089 + IT_0122 + IT_0155;
    const ccomplex_t IT_0360 = IT_0297*(IT_0349 + IT_0350) + IT_0294*(conj
      (IT_0349) + conj(IT_0350)) + IT_0303*IT_0353 + IT_0306*IT_0354 + IT_0333*
      (IT_0355 + IT_0356) + IT_0330*(conj(IT_0355) + conj(IT_0356)) + IT_0344*
      (IT_0357 + IT_0358 + IT_0359) + IT_0341*(conj(IT_0357) + conj(IT_0358) +
       conj(IT_0359));
    const ccomplex_t IT_0361 = IT_0000*IT_0309;
    const ccomplex_t IT_0362 = 2*IT_0325;
    const ccomplex_t IT_0363 = (-2)*IT_0314;
    const ccomplex_t IT_0364 = IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = 2*conj(IT_0325);
    const ccomplex_t IT_0366 = (-2)*conj(IT_0314);
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = 2*IT_0314;
    const ccomplex_t IT_0369 = (-2)*IT_0325;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 2*conj(IT_0314);
    const ccomplex_t IT_0372 = (-2)*conj(IT_0325);
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = pow(s_24, 2);
    const ccomplex_t IT_0375 = IT_0001*IT_0374;
    const ccomplex_t IT_0376 = -IT_0375;
    const ccomplex_t IT_0377 = IT_0000 + IT_0376;
    const ccomplex_t IT_0378 = IT_0000*IT_0377;
    const ccomplex_t IT_0379 = -IT_0378;
    const ccomplex_t IT_0380 = m_N_1*IT_0377;
    const ccomplex_t IT_0381 = -IT_0380;
    const ccomplex_t IT_0382 = m_N_1*IT_0036;
    const ccomplex_t IT_0383 = m_N_1*IT_0298;
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = IT_0382 + IT_0384;
    const ccomplex_t IT_0386 = IT_0312*IT_0385;
    const ccomplex_t IT_0387 = IT_0143*IT_0224;
    const ccomplex_t IT_0388 = m_N_4*IT_0387;
    const ccomplex_t IT_0389 = m_N_1*IT_0152;
    const ccomplex_t IT_0390 = -IT_0389;
    const ccomplex_t IT_0391 = IT_0388 + IT_0390;
    const ccomplex_t IT_0392 = IT_0323*IT_0391;
    const ccomplex_t IT_0393 = cpow(IT_0009, -4);
    const ccomplex_t IT_0394 = IT_0008*IT_0018*IT_0041*IT_0393;
    const ccomplex_t IT_0395 = IT_0014*IT_0018*IT_0044*IT_0393;
    const ccomplex_t IT_0396 = cpow(IT_0018, -3);
    const ccomplex_t IT_0397 = IT_0008*IT_0041*IT_0396;
    const ccomplex_t IT_0398 = IT_0014*IT_0044*IT_0396;
    const ccomplex_t IT_0399 = cpow(IT_0009, -2);
    const ccomplex_t IT_0400 = IT_0008*IT_0019*IT_0041*IT_0399;
    const ccomplex_t IT_0401 = IT_0014*IT_0019*IT_0044*IT_0399;
    const ccomplex_t IT_0402 = cpow(IT_0018, -2);
    const ccomplex_t IT_0403 = cpow(IT_0399 + IT_0402, -1);
    const ccomplex_t IT_0404 = (IT_0394 + IT_0395 + IT_0397 + IT_0398 + 2
      *IT_0400 + 2*IT_0401)*IT_0403;
    const ccomplex_t IT_0405 = m_W*e_em;
    const ccomplex_t IT_0406 = IT_0404*IT_0405;
    const ccomplex_t IT_0407 = (0 + _Complex_I*1)*IT_0406;
    const ccomplex_t IT_0408 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0409 = IT_0407*IT_0408;
    const ccomplex_t IT_0410 = (0 + _Complex_I*1)*IT_0409;
    const ccomplex_t IT_0411 = IT_0227*IT_0410;
    const ccomplex_t IT_0412 = m_N_1*IT_0228;
    const ccomplex_t IT_0413 = m_N_1*IT_0185;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0412 + IT_0414;
    const ccomplex_t IT_0416 = IT_0039*IT_0415;
    const ccomplex_t IT_0417 = m_N_1*IT_0086;
    const ccomplex_t IT_0418 = -IT_0417;
    const ccomplex_t IT_0419 = IT_0077*IT_0258;
    const ccomplex_t IT_0420 = m_N_2*IT_0419;
    const ccomplex_t IT_0421 = IT_0418 + IT_0420;
    const ccomplex_t IT_0422 = IT_0316*IT_0421;
    const ccomplex_t IT_0423 = IT_0110*IT_0289;
    const ccomplex_t IT_0424 = m_N_3*IT_0423;
    const ccomplex_t IT_0425 = m_N_1*IT_0119;
    const ccomplex_t IT_0426 = -IT_0425;
    const ccomplex_t IT_0427 = IT_0424 + IT_0426;
    const ccomplex_t IT_0428 = IT_0320*IT_0427;
    const ccomplex_t IT_0429 = (-0.5)*IT_0386 + IT_0392 + 0.5*IT_0411 + 0.5
      *IT_0416 + IT_0422 + IT_0428;
    const ccomplex_t IT_0430 = IT_0151*IT_0216;
    const ccomplex_t IT_0431 = m_N_4*IT_0430;
    const ccomplex_t IT_0432 = m_N_1*IT_0225;
    const ccomplex_t IT_0433 = -IT_0432;
    const ccomplex_t IT_0434 = IT_0431 + IT_0433;
    const ccomplex_t IT_0435 = IT_0154*IT_0434;
    const ccomplex_t IT_0436 = IT_0035*IT_0410;
    const ccomplex_t IT_0437 = IT_0085*IT_0250;
    const ccomplex_t IT_0438 = m_N_2*IT_0437;
    const ccomplex_t IT_0439 = m_N_1*IT_0259;
    const ccomplex_t IT_0440 = -IT_0439;
    const ccomplex_t IT_0441 = IT_0438 + IT_0440;
    const ccomplex_t IT_0442 = IT_0088*IT_0441;
    const ccomplex_t IT_0443 = IT_0118*IT_0281;
    const ccomplex_t IT_0444 = m_N_3*IT_0443;
    const ccomplex_t IT_0445 = m_N_1*IT_0290;
    const ccomplex_t IT_0446 = -IT_0445;
    const ccomplex_t IT_0447 = IT_0444 + IT_0446;
    const ccomplex_t IT_0448 = IT_0121*IT_0447;
    const ccomplex_t IT_0449 = (-0.5)*IT_0386 + 0.5*IT_0416 + -IT_0435 + (-0.5
      )*IT_0436 + -IT_0442 + -IT_0448;
    const ccomplex_t IT_0450 = conj(IT_0429) + conj(IT_0449);
    const ccomplex_t IT_0451 = IT_0429 + IT_0449;
    const ccomplex_t IT_0452 = IT_0039*IT_0385;
    const ccomplex_t IT_0453 = IT_0316*IT_0441;
    const ccomplex_t IT_0454 = IT_0312*IT_0415;
    const ccomplex_t IT_0455 = IT_0320*IT_0447;
    const ccomplex_t IT_0456 = IT_0323*IT_0434;
    const ccomplex_t IT_0457 = 0.5*IT_0436 + 0.5*IT_0452 + IT_0453 + (-0.5)
      *IT_0454 + IT_0455 + IT_0456;
    const ccomplex_t IT_0458 = IT_0088*IT_0421;
    const ccomplex_t IT_0459 = IT_0121*IT_0427;
    const ccomplex_t IT_0460 = IT_0154*IT_0391;
    const ccomplex_t IT_0461 = (-0.5)*IT_0411 + 0.5*IT_0452 + (-0.5)*IT_0454 +
       -IT_0458 + -IT_0459 + -IT_0460;
    const ccomplex_t IT_0462 = conj(IT_0457) + conj(IT_0461);
    const ccomplex_t IT_0463 = 2*IT_0327;
    const ccomplex_t IT_0464 = (-2)*IT_0338;
    const ccomplex_t IT_0465 = IT_0463 + IT_0464;
    const ccomplex_t IT_0466 = IT_0457 + IT_0461;
    const ccomplex_t IT_0467 = 2*conj(IT_0327);
    const ccomplex_t IT_0468 = (-2)*conj(IT_0338);
    const ccomplex_t IT_0469 = IT_0467 + IT_0468;
    const ccomplex_t IT_0470 = s_24*IT_0164;
    const ccomplex_t IT_0471 = IT_0001*IT_0470;
    const ccomplex_t IT_0472 = (-0.25)*IT_0471;
    const ccomplex_t IT_0473 = s_24 + IT_0472;
    const ccomplex_t IT_0474 = m_N_1*IT_0473;
    const ccomplex_t IT_0475 = (-4)*IT_0474;
    const ccomplex_t IT_0476 = 2*IT_0429;
    const ccomplex_t IT_0477 = 2*IT_0449;
    const ccomplex_t IT_0478 = IT_0476 + IT_0477;
    const ccomplex_t IT_0479 = 2*conj(IT_0429);
    const ccomplex_t IT_0480 = 2*conj(IT_0449);
    const ccomplex_t IT_0481 = IT_0479 + IT_0480;
    const ccomplex_t IT_0482 = conj(IT_0357) + conj(IT_0359);
    const ccomplex_t IT_0483 = 2*IT_0457;
    const ccomplex_t IT_0484 = 2*IT_0461;
    const ccomplex_t IT_0485 = IT_0483 + IT_0484;
    const ccomplex_t IT_0486 = IT_0357 + IT_0359;
    const ccomplex_t IT_0487 = 2*conj(IT_0457);
    const ccomplex_t IT_0488 = 2*conj(IT_0461);
    const ccomplex_t IT_0489 = IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = (-2)*IT_0429;
    const ccomplex_t IT_0491 = (-2)*IT_0449;
    const ccomplex_t IT_0492 = IT_0490 + IT_0491;
    const ccomplex_t IT_0493 = (-2)*conj(IT_0429);
    const ccomplex_t IT_0494 = (-2)*conj(IT_0449);
    const ccomplex_t IT_0495 = IT_0493 + IT_0494;
    const ccomplex_t IT_0496 = (-2)*IT_0457;
    const ccomplex_t IT_0497 = (-2)*IT_0461;
    const ccomplex_t IT_0498 = IT_0496 + IT_0497;
    const ccomplex_t IT_0499 = (-2)*conj(IT_0457);
    const ccomplex_t IT_0500 = (-2)*conj(IT_0461);
    const ccomplex_t IT_0501 = IT_0499 + IT_0500;
    const ccomplex_t IT_0502 = s_14*IT_0164;
    const ccomplex_t IT_0503 = IT_0001*IT_0502;
    const ccomplex_t IT_0504 = (-0.25)*IT_0503;
    const ccomplex_t IT_0505 = s_14 + IT_0504;
    const ccomplex_t IT_0506 = m_N_1*IT_0505;
    const ccomplex_t IT_0507 = (-4)*IT_0506;
    const ccomplex_t IT_0508 = 0.5*IT_0471;
    const ccomplex_t IT_0509 = s_24 + IT_0508;
    const ccomplex_t IT_0510 = m_N_1*IT_0509;
    const ccomplex_t IT_0511 = 2*IT_0510;
    const ccomplex_t IT_0512 = 0.5*IT_0503;
    const ccomplex_t IT_0513 = s_14 + IT_0512;
    const ccomplex_t IT_0514 = m_N_1*IT_0513;
    const ccomplex_t IT_0515 = 2*IT_0514;
    const ccomplex_t IT_0516 = s_12*IT_0164;
    const ccomplex_t IT_0517 = s_14*s_24;
    const ccomplex_t IT_0518 = (-2)*IT_0517;
    const ccomplex_t IT_0519 = IT_0516 + IT_0518;
    const ccomplex_t IT_0520 = IT_0001*IT_0519;
    const ccomplex_t IT_0521 = (-0.5)*IT_0520;
    const ccomplex_t IT_0522 = s_12 + IT_0521;
    const ccomplex_t IT_0523 = 4*IT_0457;
    const ccomplex_t IT_0524 = 6*conj(IT_0429);
    const ccomplex_t IT_0525 = m_N_1*IT_0005;
    const ccomplex_t IT_0526 = s_12*IT_0005;
    const ccomplex_t IT_0527 = m_N_1*IT_0309;
    const ccomplex_t IT_0528 = IT_0297*IT_0429 + IT_0294*conj(IT_0429) + 
      (IT_0328 + IT_0329 + IT_0368 + IT_0369)*IT_0462 + conj(IT_0449)*(IT_0195 +
       IT_0293 + IT_0463 + IT_0464) + (IT_0331 + IT_0332 + IT_0371 + IT_0372)
      *IT_0466 + IT_0449*(IT_0295 + IT_0296 + IT_0467 + IT_0468) + conj(IT_0327)
      *IT_0476 + IT_0327*IT_0479 + conj(IT_0338)*IT_0490 + IT_0338*IT_0493;
    const ccomplex_t IT_0529 = s_12*IT_0377;
    const ccomplex_t IT_0530 = -IT_0529;
    const ccomplex_t IT_0531 = 2*IT_0338;
    const ccomplex_t IT_0532 = s_12*s_24;
    const ccomplex_t IT_0533 = s_14*IT_0000;
    const ccomplex_t IT_0534 = -IT_0533;
    const ccomplex_t IT_0535 = IT_0532 + IT_0534;
    const ccomplex_t IT_0536 = (-2)*IT_0327;
    const ccomplex_t IT_0537 = 2*conj(IT_0338);
    const ccomplex_t IT_0538 = (-2)*conj(IT_0327);
    const ccomplex_t IT_0539 = IT_0353*IT_0364 + IT_0354*IT_0367 + (conj
      (IT_0349) + conj(IT_0350))*IT_0370 + (IT_0349 + IT_0350)*IT_0373 + (conj
      (IT_0355) + conj(IT_0356))*IT_0465 + (IT_0355 + IT_0356)*IT_0469 + (conj
      (IT_0357) + conj(IT_0358) + conj(IT_0359))*(IT_0531 + IT_0536) + (IT_0357 
      + IT_0358 + IT_0359)*(IT_0537 + IT_0538);
    const ccomplex_t IT_0540 = IT_0000*IT_0164;
    const ccomplex_t IT_0541 = 3*IT_0540;
    const ccomplex_t IT_0542 = 2*IT_0349;
    const ccomplex_t IT_0543 = (-2)*IT_0352;
    const ccomplex_t IT_0544 = IT_0542 + IT_0543;
    const ccomplex_t IT_0545 = 2*conj(IT_0349);
    const ccomplex_t IT_0546 = (-2)*conj(IT_0352);
    const ccomplex_t IT_0547 = IT_0545 + IT_0546;
    const ccomplex_t IT_0548 = 2*IT_0351;
    const ccomplex_t IT_0549 = (-2)*IT_0350;
    const ccomplex_t IT_0550 = IT_0548 + IT_0549;
    const ccomplex_t IT_0551 = 2*conj(IT_0351);
    const ccomplex_t IT_0552 = (-2)*conj(IT_0350);
    const ccomplex_t IT_0553 = IT_0551 + IT_0552;
    const ccomplex_t IT_0554 = 2*IT_0355;
    const ccomplex_t IT_0555 = 2*conj(IT_0355);
    const ccomplex_t IT_0556 = 2*IT_0358;
    const ccomplex_t IT_0557 = 2*conj(IT_0358);
    const ccomplex_t IT_0558 = (-2)*IT_0349;
    const ccomplex_t IT_0559 = (-2)*conj(IT_0349);
    const ccomplex_t IT_0560 = (-2)*IT_0351;
    const ccomplex_t IT_0561 = (-2)*conj(IT_0351);
    const ccomplex_t IT_0562 = conj(IT_0356)*IT_0544 + IT_0356*IT_0547 +
       IT_0482*IT_0550 + IT_0486*IT_0553 + conj(IT_0350)*IT_0554 + IT_0350
      *IT_0555 + conj(IT_0352)*IT_0556 + IT_0352*IT_0557 + conj(IT_0358)*IT_0558
       + IT_0358*IT_0559 + conj(IT_0355)*IT_0560 + IT_0355*IT_0561;
    const ccomplex_t IT_0563 = (-3)*IT_0540;
    const ccomplex_t IT_0564 = IT_0482*IT_0544 + IT_0486*IT_0547 + conj
      (IT_0356)*IT_0550 + IT_0356*IT_0553 + conj(IT_0352)*IT_0554 + IT_0352
      *IT_0555 + conj(IT_0350)*IT_0556 + IT_0350*IT_0557 + conj(IT_0355)*IT_0558
       + IT_0355*IT_0559 + conj(IT_0358)*IT_0560 + IT_0358*IT_0561;
    const ccomplex_t IT_0565 = pow(m_Z, 4);
    const ccomplex_t IT_0566 = s_12*IT_0565;
    const ccomplex_t IT_0567 = IT_0001*IT_0566;
    const ccomplex_t IT_0568 = (-0.5)*IT_0567;
    const ccomplex_t IT_0569 = IT_0516 + IT_0518 + IT_0568;
    const ccomplex_t IT_0570 = -IT_0351;
    const ccomplex_t IT_0571 = IT_0350 + IT_0570;
    const ccomplex_t IT_0572 = -conj(IT_0351);
    const ccomplex_t IT_0573 = conj(IT_0350) + IT_0572;
    const ccomplex_t IT_0574 = 4*conj(IT_0355);
    const ccomplex_t IT_0575 = 4*IT_0355;
    const ccomplex_t IT_0576 = (-4)*IT_0359;
    const ccomplex_t IT_0577 = -conj(IT_0352);
    const ccomplex_t IT_0578 = conj(IT_0349) + IT_0577;
    const ccomplex_t IT_0579 = -IT_0352;
    const ccomplex_t IT_0580 = IT_0355*conj(IT_0355) + (IT_0356 + -IT_0358)*
      (conj(IT_0356) + -conj(IT_0358)) + IT_0571*IT_0573 + (-0.25)*IT_0359
      *IT_0574 + IT_0357*(IT_0482 + (-0.25)*IT_0574) + (-0.25)*IT_0482*(IT_0575 
      + IT_0576) + IT_0578*(IT_0349 + IT_0579);
    const ccomplex_t IT_0581 = (-0.25)*IT_0567;
    const ccomplex_t IT_0582 = IT_0517 + IT_0581;
    const ccomplex_t IT_0583 = 8*IT_0358;
    const ccomplex_t IT_0584 = 8*conj(IT_0358);
    const ccomplex_t IT_0585 = conj(IT_0355)*IT_0358 + IT_0355*conj(IT_0358) +
       -conj(IT_0356)*(IT_0355 + -IT_0357 + -IT_0359) + -IT_0356*(conj(IT_0355) 
      + -conj(IT_0357) + -conj(IT_0359)) + -IT_0571*IT_0578 + -IT_0573*(IT_0349 
      + IT_0579) + (-0.125)*IT_0482*IT_0583 + (-0.125)*IT_0486*IT_0584;
    const ccomplex_t IT_0586 = IT_0006*(conj(IT_0168)*IT_0294 + IT_0168
      *IT_0297 + conj(IT_0300)*IT_0303 + IT_0300*IT_0306) + IT_0310*(IT_0297
      *IT_0314 + IT_0294*conj(IT_0314) + IT_0306*IT_0325 + IT_0303*conj(IT_0325)
       + conj(IT_0327)*IT_0330 + IT_0327*IT_0333 + conj(IT_0338)*IT_0341 +
       IT_0338*IT_0344) + IT_0348*IT_0360 + IT_0361*(IT_0306*IT_0327 + IT_0303
      *conj(IT_0327) + IT_0297*IT_0338 + IT_0294*conj(IT_0338) + conj(IT_0300)
      *IT_0364 + IT_0300*IT_0367 + conj(IT_0168)*IT_0370 + IT_0168*IT_0373) + 
      (conj(IT_0327)*IT_0364 + IT_0327*IT_0367 + conj(IT_0338)*IT_0370 + IT_0338
      *IT_0373)*IT_0379 + IT_0381*(IT_0370*IT_0450 + IT_0373*IT_0451 + IT_0462
      *IT_0465 + IT_0466*IT_0469) + IT_0475*(conj(IT_0352)*IT_0478 + IT_0352
      *IT_0481 + IT_0482*IT_0485 + IT_0486*IT_0489 + conj(IT_0349)*IT_0492 +
       IT_0349*IT_0495 + conj(IT_0355)*IT_0498 + IT_0355*IT_0501) + (conj
      (IT_0358)*IT_0478 + IT_0358*IT_0481 + conj(IT_0351)*IT_0485 + IT_0351
      *IT_0489 + conj(IT_0356)*IT_0492 + IT_0356*IT_0495 + conj(IT_0350)*IT_0498
       + IT_0350*IT_0501)*IT_0507 + (conj(IT_0350)*IT_0478 + IT_0350*IT_0481 +
       conj(IT_0356)*IT_0485 + IT_0356*IT_0489 + conj(IT_0351)*IT_0492 + IT_0351
      *IT_0495 + conj(IT_0358)*IT_0498 + IT_0358*IT_0501)*IT_0511 + (conj
      (IT_0355)*IT_0478 + IT_0355*IT_0481 + conj(IT_0349)*IT_0485 + IT_0349
      *IT_0489 + IT_0482*IT_0492 + IT_0486*IT_0495 + conj(IT_0352)*IT_0498 +
       IT_0352*IT_0501)*IT_0515 + 4*IT_0522*(IT_0450*IT_0451 + IT_0462*(IT_0461 
      + 0.25*IT_0523)) + 6*IT_0000*(IT_0451*IT_0462 + IT_0466*(conj(IT_0449) +
       0.166666666666667*IT_0524)) + (IT_0333*IT_0449 + IT_0330*conj(IT_0449) +
       IT_0294*IT_0462 + IT_0297*IT_0466 + conj(IT_0300)*IT_0476 + IT_0300
      *IT_0479 + conj(IT_0168)*IT_0490 + IT_0168*IT_0493)*IT_0525 + (conj
      (IT_0194)*IT_0195 + conj(IT_0292)*IT_0303 + IT_0292*IT_0305 + conj(IT_0300
      )*IT_0328 + conj(IT_0168)*IT_0341 + IT_0168*IT_0343)*IT_0526 + IT_0527
      *IT_0528 + IT_0530*(conj(IT_0325)*IT_0362 + conj(IT_0314)*IT_0370 +
       IT_0314*IT_0372 + conj(IT_0327)*IT_0465 + IT_0327*IT_0468 + conj(IT_0338)
      *IT_0531) + IT_0535*IT_0539 + IT_0541*IT_0562 + IT_0563*IT_0564 + (-4)
      *IT_0569*IT_0580 + (-8)*IT_0582*IT_0585;
    return create_ccomplex_return(IT_0586);
}

