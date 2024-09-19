#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_A0_A0(
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
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
    const ccomplex_t IT_0000 = pow(m_N_2, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0012;
    const ccomplex_t IT_0017 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0006 + -IT_0010 + 
      -IT_0015 + IT_0018);
    const ccomplex_t IT_0020 = cpow(IT_0002, -2);
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = IT_0012 + IT_0021;
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0023, 3);
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0023, 2);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0007
      *IT_0024 + -IT_0026*((IT_0007*IT_0023 + -IT_0001*IT_0025)*IT_0026 + 
      -IT_0001*IT_0027));
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0019*IT_0032;
    const ccomplex_t IT_0034 = IT_0003*IT_0025;
    const ccomplex_t IT_0035 = IT_0009*IT_0034;
    const ccomplex_t IT_0036 = IT_0003*IT_0023;
    const ccomplex_t IT_0037 = IT_0005*IT_0036;
    const ccomplex_t IT_0038 = IT_0012*IT_0025;
    const ccomplex_t IT_0039 = IT_0017*IT_0038;
    const ccomplex_t IT_0040 = IT_0012*IT_0023;
    const ccomplex_t IT_0041 = IT_0014*IT_0040;
    const ccomplex_t IT_0042 = IT_0035 + -IT_0037 + -IT_0039 + IT_0041;
    const ccomplex_t IT_0043 = cpow(IT_0042, 2);
    const ccomplex_t IT_0044 = m_N_2*IT_0043;
    const ccomplex_t IT_0045 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0046 = IT_0034*IT_0045;
    const ccomplex_t IT_0047 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0048 = IT_0036*IT_0047;
    const ccomplex_t IT_0049 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0050 = IT_0038*IT_0049;
    const ccomplex_t IT_0051 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0052 = IT_0040*IT_0051;
    const ccomplex_t IT_0053 = -IT_0046 + IT_0048 + IT_0050 + -IT_0052;
    const ccomplex_t IT_0054 = IT_0042*IT_0053;
    const ccomplex_t IT_0055 = m_N_2*IT_0054;
    const ccomplex_t IT_0056 = IT_0044 + IT_0055;
    const ccomplex_t IT_0057 = pow(m_A0, 2);
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0057 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0056*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0000 + IT_0057 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0064 = IT_0034*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0066 = IT_0036*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0068 = IT_0038*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0070 = IT_0040*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0072 = IT_0034*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0074 = IT_0036*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0076 = IT_0038*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0078 = IT_0040*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0064 + (-0.5)*IT_0066 + (-0.5)*IT_0068 +
       0.5*IT_0070 + 0.5*IT_0072 + (-0.5)*IT_0074 + (-0.5)*IT_0076 + 0.5*IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0079, 2);
    const ccomplex_t IT_0081 = m_N_3*IT_0080;
    const ccomplex_t IT_0082 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0083 = IT_0034*IT_0082;
    const ccomplex_t IT_0084 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0085 = IT_0036*IT_0084;
    const ccomplex_t IT_0086 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0087 = IT_0038*IT_0086;
    const ccomplex_t IT_0088 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0089 = IT_0040*IT_0088;
    const ccomplex_t IT_0090 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0091 = IT_0034*IT_0090;
    const ccomplex_t IT_0092 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0093 = IT_0036*IT_0092;
    const ccomplex_t IT_0094 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0095 = IT_0038*IT_0094;
    const ccomplex_t IT_0096 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0097 = IT_0040*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0083 + 0.5*IT_0085 + 0.5*IT_0087 + (
      -0.5)*IT_0089 + (-0.5)*IT_0091 + 0.5*IT_0093 + 0.5*IT_0095 + (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0079*IT_0098;
    const ccomplex_t IT_0100 = m_N_2*IT_0099;
    const ccomplex_t IT_0101 = IT_0081 + IT_0100;
    const ccomplex_t IT_0102 = IT_0062*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0000 + IT_0057 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0101*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_13 + IT_0000 + IT_0057 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0109 = IT_0034*IT_0108;
    const ccomplex_t IT_0110 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0111 = IT_0036*IT_0110;
    const ccomplex_t IT_0112 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0113 = IT_0038*IT_0112;
    const ccomplex_t IT_0114 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0115 = IT_0040*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0117 = IT_0034*IT_0116;
    const ccomplex_t IT_0118 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0119 = IT_0036*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0121 = IT_0038*IT_0120;
    const ccomplex_t IT_0122 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0123 = IT_0040*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0109 + (-0.5)*IT_0111 + (-0.5)*IT_0113 +
       0.5*IT_0115 + 0.5*IT_0117 + (-0.5)*IT_0119 + (-0.5)*IT_0121 + 0.5*IT_0123;
    const ccomplex_t IT_0125 = cpow(IT_0124, 2);
    const ccomplex_t IT_0126 = m_N_4*IT_0125;
    const ccomplex_t IT_0127 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0128 = IT_0034*IT_0127;
    const ccomplex_t IT_0129 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0130 = IT_0036*IT_0129;
    const ccomplex_t IT_0131 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0132 = IT_0038*IT_0131;
    const ccomplex_t IT_0133 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0134 = IT_0040*IT_0133;
    const ccomplex_t IT_0135 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0136 = IT_0034*IT_0135;
    const ccomplex_t IT_0137 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0138 = IT_0036*IT_0137;
    const ccomplex_t IT_0139 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0140 = IT_0038*IT_0139;
    const ccomplex_t IT_0141 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0142 = IT_0040*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0128 + 0.5*IT_0130 + 0.5*IT_0132 + (
      -0.5)*IT_0134 + (-0.5)*IT_0136 + 0.5*IT_0138 + 0.5*IT_0140 + (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = IT_0124*IT_0143;
    const ccomplex_t IT_0145 = m_N_2*IT_0144;
    const ccomplex_t IT_0146 = IT_0126 + IT_0145;
    const ccomplex_t IT_0147 = IT_0107*IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_23 + IT_0000 + IT_0057 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0146*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0033 + (-0.5)*IT_0060 + -IT_0102 + 
      -IT_0105 + -IT_0147 + -IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_13 + IT_0057 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = cpow(IT_0053, 2);
    const ccomplex_t IT_0155 = m_N_2*IT_0154;
    const ccomplex_t IT_0156 = IT_0055 + IT_0155;
    const ccomplex_t IT_0157 = IT_0153*IT_0156;
    const ccomplex_t IT_0158 = IT_0059*IT_0156;
    const ccomplex_t IT_0159 = IT_0004*IT_0045;
    const ccomplex_t IT_0160 = IT_0008*IT_0047;
    const ccomplex_t IT_0161 = IT_0013*IT_0049;
    const ccomplex_t IT_0162 = IT_0016*IT_0051;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + IT_0160 + 
      -IT_0161 + -IT_0162);
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0001
      *IT_0024 + -IT_0025*(IT_0025*(IT_0001*IT_0023 + IT_0007*IT_0025) + 
      -IT_0007*IT_0027));
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = IT_0164*IT_0169;
    const ccomplex_t IT_0171 = IT_0004*IT_0047;
    const ccomplex_t IT_0172 = IT_0008*IT_0045;
    const ccomplex_t IT_0173 = IT_0013*IT_0051;
    const ccomplex_t IT_0174 = IT_0016*IT_0049;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + -IT_0172 + 
      -IT_0173 + IT_0174);
    const ccomplex_t IT_0176 = IT_0032*IT_0175;
    const ccomplex_t IT_0177 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0178 = IT_0034*IT_0177;
    const ccomplex_t IT_0179 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0180 = IT_0036*IT_0179;
    const ccomplex_t IT_0181 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0182 = IT_0038*IT_0181;
    const ccomplex_t IT_0183 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0184 = IT_0040*IT_0183;
    const ccomplex_t IT_0185 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0186 = IT_0034*IT_0185;
    const ccomplex_t IT_0187 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0188 = IT_0036*IT_0187;
    const ccomplex_t IT_0189 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0190 = IT_0038*IT_0189;
    const ccomplex_t IT_0191 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0192 = IT_0040*IT_0191;
    const ccomplex_t IT_0193 = (-0.5)*IT_0178 + 0.5*IT_0180 + 0.5*IT_0182 + (
      -0.5)*IT_0184 + (-0.5)*IT_0186 + 0.5*IT_0188 + 0.5*IT_0190 + (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = cpow(IT_0193, 2);
    const ccomplex_t IT_0195 = m_N_1*IT_0194;
    const ccomplex_t IT_0196 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0197 = IT_0034*IT_0196;
    const ccomplex_t IT_0198 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0199 = IT_0036*IT_0198;
    const ccomplex_t IT_0200 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0201 = IT_0038*IT_0200;
    const ccomplex_t IT_0202 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0203 = IT_0040*IT_0202;
    const ccomplex_t IT_0204 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0205 = IT_0034*IT_0204;
    const ccomplex_t IT_0206 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0207 = IT_0036*IT_0206;
    const ccomplex_t IT_0208 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0209 = IT_0038*IT_0208;
    const ccomplex_t IT_0210 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0211 = IT_0040*IT_0210;
    const ccomplex_t IT_0212 = 0.5*IT_0197 + (-0.5)*IT_0199 + (-0.5)*IT_0201 +
       0.5*IT_0203 + 0.5*IT_0205 + (-0.5)*IT_0207 + (-0.5)*IT_0209 + 0.5*IT_0211;
    const ccomplex_t IT_0213 = IT_0193*IT_0212;
    const ccomplex_t IT_0214 = m_N_2*IT_0213;
    const ccomplex_t IT_0215 = IT_0195 + IT_0214;
    const ccomplex_t IT_0216 = cpow((-2)*s_13 + IT_0000 + IT_0057 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0215*IT_0217;
    const ccomplex_t IT_0219 = cpow((-2)*s_23 + IT_0000 + IT_0057 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = IT_0215*IT_0220;
    const ccomplex_t IT_0222 = IT_0157 + 0.5*IT_0158 + -IT_0170 + (-0.5)
      *IT_0176 + IT_0218 + IT_0221;
    const ccomplex_t IT_0223 = 2*IT_0222;
    const ccomplex_t IT_0224 = cpow(IT_0098, 2);
    const ccomplex_t IT_0225 = m_N_3*IT_0224;
    const ccomplex_t IT_0226 = IT_0100 + IT_0225;
    const ccomplex_t IT_0227 = IT_0062*IT_0226;
    const ccomplex_t IT_0228 = IT_0104*IT_0226;
    const ccomplex_t IT_0229 = cpow(IT_0143, 2);
    const ccomplex_t IT_0230 = m_N_4*IT_0229;
    const ccomplex_t IT_0231 = IT_0145 + IT_0230;
    const ccomplex_t IT_0232 = IT_0107*IT_0231;
    const ccomplex_t IT_0233 = IT_0149*IT_0231;
    const ccomplex_t IT_0234 = (-0.5)*IT_0158 + 0.5*IT_0176 + -IT_0227 + 
      -IT_0228 + -IT_0232 + -IT_0233;
    const ccomplex_t IT_0235 = (-2)*IT_0234;
    const ccomplex_t IT_0236 = IT_0223 + IT_0235;
    const ccomplex_t IT_0237 = 2*conj(IT_0222);
    const ccomplex_t IT_0238 = (-2)*conj(IT_0234);
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = IT_0004*IT_0009;
    const ccomplex_t IT_0241 = IT_0005*IT_0008;
    const ccomplex_t IT_0242 = IT_0013*IT_0017;
    const ccomplex_t IT_0243 = IT_0014*IT_0016;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0240 + IT_0241 + 
      -IT_0242 + -IT_0243);
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = IT_0169*IT_0245;
    const ccomplex_t IT_0247 = IT_0056*IT_0153;
    const ccomplex_t IT_0248 = cpow(IT_0212, 2);
    const ccomplex_t IT_0249 = m_N_1*IT_0248;
    const ccomplex_t IT_0250 = IT_0214 + IT_0249;
    const ccomplex_t IT_0251 = IT_0217*IT_0250;
    const ccomplex_t IT_0252 = IT_0220*IT_0250;
    const ccomplex_t IT_0253 = (-0.5)*IT_0033 + 0.5*IT_0060 + -IT_0246 +
       IT_0247 + IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 2*IT_0234;
    const ccomplex_t IT_0255 = (-2)*IT_0222;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 2*conj(IT_0234);
    const ccomplex_t IT_0258 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = 2*s_12;
    const ccomplex_t IT_0261 = (-2)*s_12;
    const ccomplex_t IT_0262 = s_24*m_N_2;
    const ccomplex_t IT_0263 = IT_0054*IT_0153;
    const ccomplex_t IT_0264 = IT_0054*IT_0059;
    const ccomplex_t IT_0265 = IT_0213*IT_0217;
    const ccomplex_t IT_0266 = IT_0213*IT_0220;
    const ccomplex_t IT_0267 = -IT_0263 + 0.5*IT_0264 + -IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = IT_0062*IT_0099;
    const ccomplex_t IT_0269 = IT_0099*IT_0104;
    const ccomplex_t IT_0270 = IT_0107*IT_0144;
    const ccomplex_t IT_0271 = IT_0144*IT_0149;
    const ccomplex_t IT_0272 = 0.5*IT_0264 + -IT_0268 + IT_0269 + -IT_0270 +
       IT_0271;
    const ccomplex_t IT_0273 = conj(IT_0267) + conj(IT_0272);
    const ccomplex_t IT_0274 = IT_0267 + IT_0272;
    const ccomplex_t IT_0275 = 2*IT_0151;
    const ccomplex_t IT_0276 = (-2)*IT_0253;
    const ccomplex_t IT_0277 = 2*conj(IT_0151);
    const ccomplex_t IT_0278 = (-2)*conj(IT_0253);
    const ccomplex_t IT_0279 = s_14*m_N_2;
    const ccomplex_t IT_0280 = 2*IT_0253;
    const ccomplex_t IT_0281 = (-2)*IT_0151;
    const ccomplex_t IT_0282 = 2*conj(IT_0253);
    const ccomplex_t IT_0283 = (-2)*conj(IT_0151);
    const ccomplex_t IT_0284 = s_12*IT_0057;
    const ccomplex_t IT_0285 = s_14*s_24;
    const ccomplex_t IT_0286 = (-2)*IT_0285;
    const ccomplex_t IT_0287 = IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0267;
    const ccomplex_t IT_0289 = (-2)*conj(IT_0267);
    const ccomplex_t IT_0290 = (-2)*IT_0272;
    const ccomplex_t IT_0291 = IT_0000*IT_0057;
    const ccomplex_t IT_0292 = IT_0000*(conj(IT_0151)*IT_0236 + IT_0151
      *IT_0239 + conj(IT_0253)*IT_0256 + IT_0253*IT_0259) + (IT_0151*conj
      (IT_0151) + IT_0222*conj(IT_0222) + IT_0234*conj(IT_0234) + IT_0253*conj
      (IT_0253))*IT_0260 + (conj(IT_0222)*IT_0234 + IT_0222*conj(IT_0234) + conj
      (IT_0151)*IT_0253 + IT_0151*conj(IT_0253))*IT_0261 + IT_0262*(IT_0256
      *IT_0273 + IT_0259*IT_0274 + IT_0273*(IT_0275 + IT_0276) + IT_0274*
      (IT_0277 + IT_0278)) + IT_0279*(IT_0236*IT_0273 + IT_0239*IT_0274 +
       IT_0273*(IT_0280 + IT_0281) + IT_0274*(IT_0282 + IT_0283)) + 2*IT_0287*
      (IT_0273*IT_0288 + 1./2*IT_0272*IT_0289 + -IT_0272*(conj(IT_0272) + (-0.5)
      *IT_0289) + 1./2*conj(IT_0272)*IT_0290) + (conj(IT_0267)*IT_0288 + (-2)
      *conj(IT_0272)*(IT_0272 + (-0.5)*IT_0288) + IT_0272*IT_0289 + IT_0274
      *IT_0289 + conj(IT_0272)*(IT_0288 + IT_0290))*IT_0291;
    return create_ccomplex_return(IT_0292);
}

