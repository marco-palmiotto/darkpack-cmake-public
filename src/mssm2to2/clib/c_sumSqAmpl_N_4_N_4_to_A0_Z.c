#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
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
    const ccomplex_t IT_0000 = pow(m_Z, -2);
    const ccomplex_t IT_0001 = pow(m_Z, 2);
    const ccomplex_t IT_0002 = s_14*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = s_14 + IT_0004;
    const ccomplex_t IT_0006 = m_N_4*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cos(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = sin(beta);
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0008*IT_0019;
    const ccomplex_t IT_0021 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0014*IT_0019;
    const ccomplex_t IT_0024 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0029 = IT_0015*IT_0028;
    const ccomplex_t IT_0030 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0033 = IT_0023*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0013 + 0.5*IT_0017 + 0.5*IT_0022 + (
      -0.5)*IT_0025 + (-0.5)*IT_0027 + 0.5*IT_0029 + 0.5*IT_0031 + (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0010*IT_0018;
    const ccomplex_t IT_0036 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = IT_0009*IT_0019;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0041 = IT_0035*IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0040;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0037 + IT_0039 + 
      -IT_0041 + -IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0034*IT_0045;
    const ccomplex_t IT_0047 = pow(m_A0, 2);
    const ccomplex_t IT_0048 = pow(m_N_4, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_23 + IT_0047 + IT_0048 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0046*IT_0050;
    const ccomplex_t IT_0052 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0053 = IT_0011*IT_0052;
    const ccomplex_t IT_0054 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0055 = IT_0015*IT_0054;
    const ccomplex_t IT_0056 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0057 = IT_0020*IT_0056;
    const ccomplex_t IT_0058 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0059 = IT_0023*IT_0058;
    const ccomplex_t IT_0060 = -IT_0053 + IT_0055 + IT_0057 + -IT_0059;
    const ccomplex_t IT_0061 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0062 = IT_0035*IT_0061;
    const ccomplex_t IT_0063 = IT_0038*IT_0061;
    const ccomplex_t IT_0064 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0065 = IT_0038*IT_0064;
    const ccomplex_t IT_0066 = IT_0035*IT_0064;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0062 + IT_0063 + 
      -IT_0065 + -IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0060*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_23 + IT_0047 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0070;
    const ccomplex_t IT_0077 = IT_0060*IT_0076;
    const ccomplex_t IT_0078 = IT_0074*IT_0077;
    const ccomplex_t IT_0079 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0080 = IT_0011*IT_0079;
    const ccomplex_t IT_0081 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0082 = IT_0015*IT_0081;
    const ccomplex_t IT_0083 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0084 = IT_0020*IT_0083;
    const ccomplex_t IT_0085 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0086 = IT_0023*IT_0085;
    const ccomplex_t IT_0087 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0088 = IT_0011*IT_0087;
    const ccomplex_t IT_0089 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0090 = IT_0015*IT_0089;
    const ccomplex_t IT_0091 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0092 = IT_0020*IT_0091;
    const ccomplex_t IT_0093 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0094 = IT_0023*IT_0093;
    const ccomplex_t IT_0095 = (-0.5)*IT_0080 + 0.5*IT_0082 + 0.5*IT_0084 + (
      -0.5)*IT_0086 + (-0.5)*IT_0088 + 0.5*IT_0090 + 0.5*IT_0092 + (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0097 = IT_0038*IT_0096;
    const ccomplex_t IT_0098 = IT_0035*IT_0096;
    const ccomplex_t IT_0099 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0100 = IT_0038*IT_0099;
    const ccomplex_t IT_0101 = IT_0035*IT_0099;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0097 + IT_0098 + 
      -IT_0100 + -IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0047 + IT_0048 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0112 = IT_0011*IT_0111;
    const ccomplex_t IT_0113 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0114 = IT_0015*IT_0113;
    const ccomplex_t IT_0115 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0116 = IT_0020*IT_0115;
    const ccomplex_t IT_0117 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0118 = IT_0023*IT_0117;
    const ccomplex_t IT_0119 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0120 = IT_0011*IT_0119;
    const ccomplex_t IT_0121 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0122 = IT_0015*IT_0121;
    const ccomplex_t IT_0123 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0124 = IT_0020*IT_0123;
    const ccomplex_t IT_0125 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0126 = IT_0023*IT_0125;
    const ccomplex_t IT_0127 = (-0.5)*IT_0112 + 0.5*IT_0114 + 0.5*IT_0116 + (
      -0.5)*IT_0118 + (-0.5)*IT_0120 + 0.5*IT_0122 + 0.5*IT_0124 + (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0129 = IT_0035*IT_0128;
    const ccomplex_t IT_0130 = IT_0038*IT_0128;
    const ccomplex_t IT_0131 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0132 = IT_0038*IT_0131;
    const ccomplex_t IT_0133 = IT_0035*IT_0131;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0129 + IT_0130 + 
      -IT_0132 + -IT_0133);
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = IT_0127*IT_0138;
    const ccomplex_t IT_0140 = cpow((-2)*s_23 + IT_0047 + IT_0048 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0139*IT_0141;
    const ccomplex_t IT_0143 = (-0.5)*IT_0075 + 0.5*IT_0078 + IT_0110 + IT_0142;
    const ccomplex_t IT_0144 = conj(IT_0051) + conj(IT_0143);
    const ccomplex_t IT_0145 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0146 = IT_0011*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0148 = IT_0015*IT_0147;
    const ccomplex_t IT_0149 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0150 = IT_0020*IT_0149;
    const ccomplex_t IT_0151 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0152 = IT_0023*IT_0151;
    const ccomplex_t IT_0153 = IT_0146 + -IT_0148 + -IT_0150 + IT_0152;
    const ccomplex_t IT_0154 = IT_0071*IT_0153;
    const ccomplex_t IT_0155 = m_N_4*IT_0154;
    const ccomplex_t IT_0156 = IT_0074*IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_13 + IT_0047 + IT_0048 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0160 = IT_0035*IT_0159;
    const ccomplex_t IT_0161 = IT_0038*IT_0159;
    const ccomplex_t IT_0162 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0163 = IT_0035*IT_0162;
    const ccomplex_t IT_0164 = IT_0038*IT_0162;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0160 + IT_0161 + 
      -IT_0163 + -IT_0164);
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0095*IT_0168;
    const ccomplex_t IT_0170 = m_N_1*IT_0169;
    const ccomplex_t IT_0171 = m_N_4*IT_0107;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0170 + IT_0172;
    const ccomplex_t IT_0174 = IT_0158*IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_13 + IT_0047 + IT_0048 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0178 = IT_0035*IT_0177;
    const ccomplex_t IT_0179 = IT_0038*IT_0177;
    const ccomplex_t IT_0180 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0181 = IT_0038*IT_0180;
    const ccomplex_t IT_0182 = IT_0035*IT_0180;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0178 + IT_0179 + 
      -IT_0181 + -IT_0182);
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = (-0.5)*IT_0184;
    const ccomplex_t IT_0186 = IT_0127*IT_0185;
    const ccomplex_t IT_0187 = m_N_2*IT_0186;
    const ccomplex_t IT_0188 = m_N_4*IT_0139;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0187 + IT_0189;
    const ccomplex_t IT_0191 = IT_0176*IT_0190;
    const ccomplex_t IT_0192 = cpow((-2)*s_13 + IT_0047 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_N_4*IT_0072;
    const ccomplex_t IT_0195 = m_N_4*IT_0077;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = IT_0193*IT_0197;
    const ccomplex_t IT_0199 = cpow((-2)*s_13 + IT_0047 + IT_0048 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0202 = IT_0035*IT_0201;
    const ccomplex_t IT_0203 = IT_0038*IT_0201;
    const ccomplex_t IT_0204 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0205 = IT_0035*IT_0204;
    const ccomplex_t IT_0206 = IT_0038*IT_0204;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0202 + IT_0203 + 
      -IT_0205 + -IT_0206);
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = IT_0034*IT_0208;
    const ccomplex_t IT_0210 = m_N_3*IT_0209;
    const ccomplex_t IT_0211 = m_N_4*IT_0046;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0210 + IT_0212;
    const ccomplex_t IT_0214 = IT_0200*IT_0213;
    const ccomplex_t IT_0215 = (-0.5)*IT_0156 + -IT_0174 + -IT_0191 + (-0.5)
      *IT_0198 + -IT_0214;
    const ccomplex_t IT_0216 = 2*IT_0215;
    const ccomplex_t IT_0217 = IT_0076*IT_0153;
    const ccomplex_t IT_0218 = m_N_4*IT_0217;
    const ccomplex_t IT_0219 = IT_0074*IT_0218;
    const ccomplex_t IT_0220 = -IT_0155;
    const ccomplex_t IT_0221 = IT_0218 + IT_0220;
    const ccomplex_t IT_0222 = IT_0074*IT_0221;
    const ccomplex_t IT_0223 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0224 = IT_0011*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0226 = IT_0015*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0228 = IT_0020*IT_0227;
    const ccomplex_t IT_0229 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0230 = IT_0023*IT_0229;
    const ccomplex_t IT_0231 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0232 = IT_0011*IT_0231;
    const ccomplex_t IT_0233 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0234 = IT_0015*IT_0233;
    const ccomplex_t IT_0235 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0236 = IT_0020*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0238 = IT_0023*IT_0237;
    const ccomplex_t IT_0239 = 0.5*IT_0224 + (-0.5)*IT_0226 + (-0.5)*IT_0228 +
       0.5*IT_0230 + 0.5*IT_0232 + (-0.5)*IT_0234 + (-0.5)*IT_0236 + 0.5*IT_0238;
    const ccomplex_t IT_0240 = IT_0106*IT_0239;
    const ccomplex_t IT_0241 = m_N_1*IT_0240;
    const ccomplex_t IT_0242 = IT_0168*IT_0239;
    const ccomplex_t IT_0243 = m_N_4*IT_0242;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = IT_0241 + IT_0244;
    const ccomplex_t IT_0246 = IT_0109*IT_0245;
    const ccomplex_t IT_0247 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0248 = IT_0011*IT_0247;
    const ccomplex_t IT_0249 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0250 = IT_0015*IT_0249;
    const ccomplex_t IT_0251 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0252 = IT_0020*IT_0251;
    const ccomplex_t IT_0253 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0254 = IT_0023*IT_0253;
    const ccomplex_t IT_0255 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0256 = IT_0011*IT_0255;
    const ccomplex_t IT_0257 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0258 = IT_0015*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0260 = IT_0020*IT_0259;
    const ccomplex_t IT_0261 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0262 = IT_0023*IT_0261;
    const ccomplex_t IT_0263 = 0.5*IT_0248 + (-0.5)*IT_0250 + (-0.5)*IT_0252 +
       0.5*IT_0254 + 0.5*IT_0256 + (-0.5)*IT_0258 + (-0.5)*IT_0260 + 0.5*IT_0262;
    const ccomplex_t IT_0264 = IT_0138*IT_0263;
    const ccomplex_t IT_0265 = m_N_2*IT_0264;
    const ccomplex_t IT_0266 = IT_0185*IT_0263;
    const ccomplex_t IT_0267 = m_N_4*IT_0266;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = IT_0265 + IT_0268;
    const ccomplex_t IT_0270 = IT_0141*IT_0269;
    const ccomplex_t IT_0271 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0272 = IT_0011*IT_0271;
    const ccomplex_t IT_0273 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0274 = IT_0015*IT_0273;
    const ccomplex_t IT_0275 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0276 = IT_0020*IT_0275;
    const ccomplex_t IT_0277 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0278 = IT_0023*IT_0277;
    const ccomplex_t IT_0279 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0280 = IT_0011*IT_0279;
    const ccomplex_t IT_0281 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0282 = IT_0015*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0284 = IT_0020*IT_0283;
    const ccomplex_t IT_0285 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0286 = IT_0023*IT_0285;
    const ccomplex_t IT_0287 = 0.5*IT_0272 + (-0.5)*IT_0274 + (-0.5)*IT_0276 +
       0.5*IT_0278 + 0.5*IT_0280 + (-0.5)*IT_0282 + (-0.5)*IT_0284 + 0.5*IT_0286;
    const ccomplex_t IT_0288 = IT_0208*IT_0287;
    const ccomplex_t IT_0289 = m_N_4*IT_0288;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = IT_0045*IT_0287;
    const ccomplex_t IT_0292 = m_N_3*IT_0291;
    const ccomplex_t IT_0293 = IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = IT_0050*IT_0293;
    const ccomplex_t IT_0295 = (-0.5)*IT_0198 + 0.5*IT_0219 + 0.5*IT_0222 +
       IT_0246 + IT_0270 + IT_0294;
    const ccomplex_t IT_0296 = 2*IT_0295;
    const ccomplex_t IT_0297 = IT_0216 + IT_0296;
    const ccomplex_t IT_0298 = IT_0051 + IT_0143;
    const ccomplex_t IT_0299 = 2*conj(IT_0215);
    const ccomplex_t IT_0300 = 2*conj(IT_0295);
    const ccomplex_t IT_0301 = IT_0299 + IT_0300;
    const ccomplex_t IT_0302 = IT_0050*IT_0288;
    const ccomplex_t IT_0303 = IT_0074*IT_0217;
    const ccomplex_t IT_0304 = IT_0074*IT_0154;
    const ccomplex_t IT_0305 = IT_0141*IT_0266;
    const ccomplex_t IT_0306 = IT_0109*IT_0242;
    const ccomplex_t IT_0307 = (-0.5)*IT_0303 + 0.5*IT_0304 + IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = conj(IT_0302) + conj(IT_0307);
    const ccomplex_t IT_0309 = IT_0193*IT_0221;
    const ccomplex_t IT_0310 = IT_0200*IT_0293;
    const ccomplex_t IT_0311 = IT_0074*IT_0195;
    const ccomplex_t IT_0312 = IT_0158*IT_0245;
    const ccomplex_t IT_0313 = IT_0176*IT_0269;
    const ccomplex_t IT_0314 = (-0.5)*IT_0309 + -IT_0310 + (-0.5)*IT_0311 + 
      -IT_0312 + -IT_0313;
    const ccomplex_t IT_0315 = 2*IT_0314;
    const ccomplex_t IT_0316 = IT_0074*IT_0194;
    const ccomplex_t IT_0317 = IT_0074*IT_0197;
    const ccomplex_t IT_0318 = IT_0109*IT_0173;
    const ccomplex_t IT_0319 = IT_0141*IT_0190;
    const ccomplex_t IT_0320 = IT_0050*IT_0213;
    const ccomplex_t IT_0321 = (-0.5)*IT_0309 + 0.5*IT_0316 + 0.5*IT_0317 +
       IT_0318 + IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = 2*IT_0321;
    const ccomplex_t IT_0323 = IT_0315 + IT_0322;
    const ccomplex_t IT_0324 = IT_0302 + IT_0307;
    const ccomplex_t IT_0325 = 2*conj(IT_0314);
    const ccomplex_t IT_0326 = 2*conj(IT_0321);
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = pow(s_14, 2);
    const ccomplex_t IT_0329 = IT_0000*IT_0328;
    const ccomplex_t IT_0330 = -IT_0048;
    const ccomplex_t IT_0331 = IT_0329 + IT_0330;
    const ccomplex_t IT_0332 = IT_0048*IT_0331;
    const ccomplex_t IT_0333 = cos(alpha);
    const ccomplex_t IT_0334 = IT_0010*IT_0333;
    const ccomplex_t IT_0335 = IT_0147*IT_0334;
    const ccomplex_t IT_0336 = sin(alpha);
    const ccomplex_t IT_0337 = IT_0010*IT_0336;
    const ccomplex_t IT_0338 = IT_0145*IT_0337;
    const ccomplex_t IT_0339 = IT_0019*IT_0333;
    const ccomplex_t IT_0340 = IT_0151*IT_0339;
    const ccomplex_t IT_0341 = IT_0019*IT_0336;
    const ccomplex_t IT_0342 = IT_0149*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*(IT_0335 + -IT_0338 + 
      -IT_0340 + IT_0342);
    const ccomplex_t IT_0344 = e_em*IT_0010*IT_0019*(IT_0014*IT_0333 + 
      -IT_0008*IT_0336);
    const ccomplex_t IT_0345 = IT_0343*IT_0344;
    const ccomplex_t IT_0346 = cpow((-2)*s_12 + (-2)*IT_0048 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0347 = (0 + _Complex_I*1)*IT_0346;
    const ccomplex_t IT_0348 = IT_0345*IT_0347;
    const ccomplex_t IT_0349 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0048 + -reg_prop, -1);
    const ccomplex_t IT_0350 = (0 + _Complex_I*1)*IT_0349;
    const ccomplex_t IT_0351 = e_em*IT_0010*IT_0019*(IT_0008*IT_0333 + IT_0014
      *IT_0336);
    const ccomplex_t IT_0352 = -IT_0351;
    const ccomplex_t IT_0353 = IT_0145*IT_0334;
    const ccomplex_t IT_0354 = IT_0147*IT_0337;
    const ccomplex_t IT_0355 = IT_0149*IT_0339;
    const ccomplex_t IT_0356 = IT_0151*IT_0341;
    const ccomplex_t IT_0357 = (0 + _Complex_I*1)*(IT_0353 + IT_0354 + 
      -IT_0355 + -IT_0356);
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = IT_0352*IT_0358;
    const ccomplex_t IT_0360 = IT_0350*IT_0359;
    const ccomplex_t IT_0361 = (-2)*IT_0302 + IT_0303 + -IT_0304 + (-2)
      *IT_0305 + (-2)*IT_0306 + (-0.5)*IT_0348 + (-0.5)*IT_0360;
    const ccomplex_t IT_0362 = IT_0054*IT_0334;
    const ccomplex_t IT_0363 = IT_0052*IT_0337;
    const ccomplex_t IT_0364 = IT_0058*IT_0339;
    const ccomplex_t IT_0365 = IT_0056*IT_0341;
    const ccomplex_t IT_0366 = (0 + _Complex_I*1)*(IT_0362 + -IT_0363 + 
      -IT_0364 + IT_0365);
    const ccomplex_t IT_0367 = IT_0344*IT_0366;
    const ccomplex_t IT_0368 = IT_0347*IT_0367;
    const ccomplex_t IT_0369 = IT_0052*IT_0334;
    const ccomplex_t IT_0370 = IT_0054*IT_0337;
    const ccomplex_t IT_0371 = IT_0056*IT_0339;
    const ccomplex_t IT_0372 = IT_0058*IT_0341;
    const ccomplex_t IT_0373 = (0 + _Complex_I*1)*(IT_0369 + IT_0370 + 
      -IT_0371 + -IT_0372);
    const ccomplex_t IT_0374 = -IT_0373;
    const ccomplex_t IT_0375 = IT_0352*IT_0374;
    const ccomplex_t IT_0376 = IT_0350*IT_0375;
    const ccomplex_t IT_0377 = 0.5*IT_0368 + 0.5*IT_0376;
    const ccomplex_t IT_0378 = 2*IT_0377;
    const ccomplex_t IT_0379 = (-2)*IT_0051 + IT_0075 + -IT_0078 + (-2)
      *IT_0110 + (-2)*IT_0142 + (-0.5)*IT_0368 + (-0.5)*IT_0376;
    const ccomplex_t IT_0380 = (-2)*IT_0379;
    const ccomplex_t IT_0381 = IT_0378 + IT_0380;
    const ccomplex_t IT_0382 = 2*conj(IT_0377);
    const ccomplex_t IT_0383 = (-2)*conj(IT_0379);
    const ccomplex_t IT_0384 = IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = 0.5*IT_0348 + 0.5*IT_0360;
    const ccomplex_t IT_0386 = 2*IT_0379;
    const ccomplex_t IT_0387 = (-2)*IT_0377;
    const ccomplex_t IT_0388 = IT_0386 + IT_0387;
    const ccomplex_t IT_0389 = 2*conj(IT_0379);
    const ccomplex_t IT_0390 = (-2)*conj(IT_0377);
    const ccomplex_t IT_0391 = IT_0389 + IT_0390;
    const ccomplex_t IT_0392 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0393 = -IT_0392;
    const ccomplex_t IT_0394 = s_12 + IT_0393;
    const ccomplex_t IT_0395 = s_12*IT_0394;
    const ccomplex_t IT_0396 = IT_0072*IT_0193;
    const ccomplex_t IT_0397 = (-0.5)*IT_0368 + (-0.5)*IT_0376 + IT_0396;
    const ccomplex_t IT_0398 = IT_0193*IT_0217;
    const ccomplex_t IT_0399 = (-0.5)*IT_0348 + (-0.5)*IT_0360 + IT_0398;
    const ccomplex_t IT_0400 = 2*IT_0385;
    const ccomplex_t IT_0401 = (-2)*IT_0361;
    const ccomplex_t IT_0402 = IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = 2*conj(IT_0385);
    const ccomplex_t IT_0404 = (-2)*conj(IT_0361);
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = IT_0107*IT_0158;
    const ccomplex_t IT_0407 = IT_0139*IT_0176;
    const ccomplex_t IT_0408 = IT_0046*IT_0200;
    const ccomplex_t IT_0409 = IT_0077*IT_0193;
    const ccomplex_t IT_0410 = 0.5*IT_0368 + 0.5*IT_0376 + 2*IT_0406 + 2
      *IT_0407 + 2*IT_0408 + IT_0409;
    const ccomplex_t IT_0411 = IT_0176*IT_0266;
    const ccomplex_t IT_0412 = IT_0200*IT_0288;
    const ccomplex_t IT_0413 = IT_0154*IT_0193;
    const ccomplex_t IT_0414 = IT_0158*IT_0242;
    const ccomplex_t IT_0415 = 0.5*IT_0348 + 0.5*IT_0360 + 2*IT_0411 + 2
      *IT_0412 + IT_0413 + 2*IT_0414;
    const ccomplex_t IT_0416 = 2*IT_0361;
    const ccomplex_t IT_0417 = (-2)*IT_0385;
    const ccomplex_t IT_0418 = IT_0416 + IT_0417;
    const ccomplex_t IT_0419 = 2*conj(IT_0361);
    const ccomplex_t IT_0420 = (-2)*conj(IT_0385);
    const ccomplex_t IT_0421 = IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = m_N_4*IT_0331;
    const ccomplex_t IT_0423 = conj(IT_0215) + conj(IT_0295);
    const ccomplex_t IT_0424 = IT_0215 + IT_0295;
    const ccomplex_t IT_0425 = conj(IT_0314) + conj(IT_0321);
    const ccomplex_t IT_0426 = IT_0314 + IT_0321;
    const ccomplex_t IT_0427 = s_12*s_14;
    const ccomplex_t IT_0428 = s_24*IT_0048;
    const ccomplex_t IT_0429 = -IT_0428;
    const ccomplex_t IT_0430 = IT_0427 + IT_0429;
    const ccomplex_t IT_0431 = (-0.5)*IT_0396;
    const ccomplex_t IT_0432 = (-0.5)*IT_0398;
    const ccomplex_t IT_0433 = -IT_0406 + -IT_0407 + -IT_0408 + (-0.5)*IT_0409;
    const ccomplex_t IT_0434 = conj(IT_0051) + conj(IT_0143) + conj(IT_0433);
    const ccomplex_t IT_0435 = IT_0051 + IT_0143 + IT_0433;
    const ccomplex_t IT_0436 = -IT_0411 + -IT_0412 + (-0.5)*IT_0413 + -IT_0414;
    const ccomplex_t IT_0437 = conj(IT_0302) + conj(IT_0307) + conj(IT_0436);
    const ccomplex_t IT_0438 = IT_0302 + IT_0307 + IT_0436;
    const ccomplex_t IT_0439 = IT_0048*IT_0394;
    const ccomplex_t IT_0440 = 2*IT_0410;
    const ccomplex_t IT_0441 = (-2)*IT_0397;
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = 2*conj(IT_0410);
    const ccomplex_t IT_0444 = (-2)*conj(IT_0397);
    const ccomplex_t IT_0445 = IT_0443 + IT_0444;
    const ccomplex_t IT_0446 = 2*IT_0397;
    const ccomplex_t IT_0447 = (-2)*IT_0410;
    const ccomplex_t IT_0448 = IT_0446 + IT_0447;
    const ccomplex_t IT_0449 = 2*conj(IT_0397);
    const ccomplex_t IT_0450 = (-2)*conj(IT_0410);
    const ccomplex_t IT_0451 = IT_0449 + IT_0450;
    const ccomplex_t IT_0452 = pow(s_24, 2);
    const ccomplex_t IT_0453 = IT_0000*IT_0452;
    const ccomplex_t IT_0454 = -IT_0453;
    const ccomplex_t IT_0455 = IT_0048 + IT_0454;
    const ccomplex_t IT_0456 = IT_0048*IT_0455;
    const ccomplex_t IT_0457 = -IT_0456;
    const ccomplex_t IT_0458 = s_24*IT_0001;
    const ccomplex_t IT_0459 = IT_0000*IT_0458;
    const ccomplex_t IT_0460 = 0.5*IT_0459;
    const ccomplex_t IT_0461 = s_24 + IT_0460;
    const ccomplex_t IT_0462 = m_N_4*IT_0461;
    const ccomplex_t IT_0463 = 2*IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0215;
    const ccomplex_t IT_0465 = (-2)*IT_0295;
    const ccomplex_t IT_0466 = IT_0464 + IT_0465;
    const ccomplex_t IT_0467 = (-2)*conj(IT_0215);
    const ccomplex_t IT_0468 = (-2)*conj(IT_0295);
    const ccomplex_t IT_0469 = IT_0467 + IT_0468;
    const ccomplex_t IT_0470 = (-2)*IT_0314;
    const ccomplex_t IT_0471 = (-2)*IT_0321;
    const ccomplex_t IT_0472 = IT_0470 + IT_0471;
    const ccomplex_t IT_0473 = (-2)*conj(IT_0314);
    const ccomplex_t IT_0474 = (-2)*conj(IT_0321);
    const ccomplex_t IT_0475 = IT_0473 + IT_0474;
    const ccomplex_t IT_0476 = (-0.25)*IT_0003;
    const ccomplex_t IT_0477 = s_14 + IT_0476;
    const ccomplex_t IT_0478 = m_N_4*IT_0477;
    const ccomplex_t IT_0479 = (-4)*IT_0478;
    const ccomplex_t IT_0480 = (-0.25)*IT_0459;
    const ccomplex_t IT_0481 = s_24 + IT_0480;
    const ccomplex_t IT_0482 = m_N_4*IT_0481;
    const ccomplex_t IT_0483 = (-4)*IT_0482;
    const ccomplex_t IT_0484 = IT_0001*IT_0048;
    const ccomplex_t IT_0485 = (-3)*IT_0484;
    const ccomplex_t IT_0486 = 2*IT_0433;
    const ccomplex_t IT_0487 = (-2)*IT_0431;
    const ccomplex_t IT_0488 = IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = 2*conj(IT_0433);
    const ccomplex_t IT_0490 = (-2)*conj(IT_0431);
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = (-2)*IT_0051;
    const ccomplex_t IT_0493 = (-2)*IT_0143;
    const ccomplex_t IT_0494 = IT_0492 + IT_0493;
    const ccomplex_t IT_0495 = (-2)*conj(IT_0051);
    const ccomplex_t IT_0496 = (-2)*conj(IT_0143);
    const ccomplex_t IT_0497 = IT_0495 + IT_0496;
    const ccomplex_t IT_0498 = 2*IT_0436;
    const ccomplex_t IT_0499 = 2*conj(IT_0436);
    const ccomplex_t IT_0500 = conj(IT_0433)*IT_0436 + IT_0433*conj(IT_0436) +
       (-0.5)*conj(IT_0432)*IT_0488 + (-0.5)*IT_0432*IT_0491 + (-0.5)*IT_0308
      *IT_0494 + (-0.5)*IT_0324*IT_0497 + (-0.5)*conj(IT_0431)*IT_0498 + (-0.5)
      *IT_0431*IT_0499;
    const ccomplex_t IT_0501 = 3*IT_0484;
    const ccomplex_t IT_0502 = s_12*IT_0001;
    const ccomplex_t IT_0503 = s_14*s_24;
    const ccomplex_t IT_0504 = (-2)*IT_0503;
    const ccomplex_t IT_0505 = IT_0502 + IT_0504;
    const ccomplex_t IT_0506 = IT_0000*IT_0505;
    const ccomplex_t IT_0507 = (-0.5)*IT_0506;
    const ccomplex_t IT_0508 = s_12 + IT_0507;
    const ccomplex_t IT_0509 = 4*IT_0314;
    const ccomplex_t IT_0510 = 6*conj(IT_0215);
    const ccomplex_t IT_0511 = s_12*IT_0331;
    const ccomplex_t IT_0512 = m_N_4*IT_0455;
    const ccomplex_t IT_0513 = -IT_0512;
    const ccomplex_t IT_0514 = 2*IT_0415;
    const ccomplex_t IT_0515 = (-2)*IT_0399;
    const ccomplex_t IT_0516 = IT_0514 + IT_0515;
    const ccomplex_t IT_0517 = 2*conj(IT_0415);
    const ccomplex_t IT_0518 = (-2)*conj(IT_0399);
    const ccomplex_t IT_0519 = IT_0517 + IT_0518;
    const ccomplex_t IT_0520 = m_N_4*IT_0394;
    const ccomplex_t IT_0521 = IT_0299*IT_0361 + IT_0216*conj(IT_0361) + conj
      (IT_0295)*(IT_0416 + IT_0417 + IT_0440 + IT_0441) + conj(IT_0215)*IT_0442 
      + IT_0295*(IT_0419 + IT_0420 + IT_0443 + IT_0444) + IT_0215*IT_0445 + conj
      (IT_0385)*IT_0464 + IT_0385*IT_0467 + IT_0425*(IT_0386 + IT_0387 + IT_0514
       + IT_0515) + IT_0426*(IT_0389 + IT_0390 + IT_0517 + IT_0518);
    const ccomplex_t IT_0522 = s_12*IT_0455;
    const ccomplex_t IT_0523 = -IT_0522;
    const ccomplex_t IT_0524 = 2*IT_0399;
    const ccomplex_t IT_0525 = (-2)*IT_0415;
    const ccomplex_t IT_0526 = IT_0524 + IT_0525;
    const ccomplex_t IT_0527 = (-2)*conj(IT_0415);
    const ccomplex_t IT_0528 = s_12*s_24;
    const ccomplex_t IT_0529 = s_14*IT_0048;
    const ccomplex_t IT_0530 = -IT_0529;
    const ccomplex_t IT_0531 = IT_0528 + IT_0530;
    const ccomplex_t IT_0532 = pow(m_Z, 4);
    const ccomplex_t IT_0533 = s_12*IT_0532;
    const ccomplex_t IT_0534 = IT_0000*IT_0533;
    const ccomplex_t IT_0535 = (-0.5)*IT_0534;
    const ccomplex_t IT_0536 = IT_0502 + IT_0504 + IT_0535;
    const ccomplex_t IT_0537 = (-4)*IT_0051;
    const ccomplex_t IT_0538 = (-4)*IT_0302;
    const ccomplex_t IT_0539 = -IT_0433;
    const ccomplex_t IT_0540 = (-4)*IT_0536*((IT_0432 + -IT_0436)*(conj
      (IT_0432) + -conj(IT_0436)) + IT_0144*(IT_0143 + (-0.25)*IT_0537) +
       IT_0308*(IT_0307 + (-0.25)*IT_0538) + (conj(IT_0431) + -conj(IT_0433))*
      (IT_0431 + IT_0539));
    const ccomplex_t IT_0541 = (-0.25)*IT_0534;
    const ccomplex_t IT_0542 = IT_0503 + IT_0541;
    const ccomplex_t IT_0543 = 8*conj(IT_0433);
    const ccomplex_t IT_0544 = 8*IT_0436;
    const ccomplex_t IT_0545 = 8*conj(IT_0436);
    const ccomplex_t IT_0546 = IT_0007*(IT_0144*IT_0297 + IT_0298*IT_0301 +
       IT_0308*IT_0323 + IT_0324*IT_0327) + IT_0332*(conj(IT_0361)*IT_0381 +
       IT_0361*IT_0384 + conj(IT_0385)*IT_0388 + IT_0385*IT_0391) + IT_0395*
      (IT_0384*IT_0397 + IT_0381*conj(IT_0397) + conj(IT_0399)*IT_0402 + IT_0399
      *IT_0405 + IT_0391*IT_0410 + IT_0388*conj(IT_0410) + conj(IT_0415)*IT_0418
       + IT_0415*IT_0421) + IT_0422*(IT_0388*IT_0423 + IT_0391*IT_0424 + IT_0418
      *IT_0425 + IT_0421*IT_0426) + IT_0430*(IT_0384*IT_0431 + IT_0381*conj
      (IT_0431) + IT_0405*IT_0432 + IT_0402*conj(IT_0432) + IT_0388*IT_0434 +
       IT_0391*IT_0435 + IT_0418*IT_0437 + IT_0421*IT_0438) + IT_0439*(IT_0391
      *IT_0399 + IT_0388*conj(IT_0399) + IT_0384*IT_0415 + IT_0381*conj(IT_0415)
       + conj(IT_0385)*IT_0442 + IT_0385*IT_0445 + conj(IT_0361)*IT_0448 +
       IT_0361*IT_0451) + (conj(IT_0399)*IT_0442 + IT_0399*IT_0445 + conj
      (IT_0415)*IT_0448 + IT_0415*IT_0451)*IT_0457 + IT_0463*(IT_0327*IT_0433 +
       IT_0323*conj(IT_0433) + IT_0301*IT_0436 + IT_0297*conj(IT_0436) + conj
      (IT_0432)*IT_0466 + IT_0432*IT_0469 + conj(IT_0431)*IT_0472 + IT_0431
      *IT_0475) + (IT_0301*IT_0431 + IT_0297*conj(IT_0431) + IT_0327*IT_0432 +
       IT_0323*conj(IT_0432) + conj(IT_0433)*IT_0466 + IT_0433*IT_0469 + conj
      (IT_0436)*IT_0472 + IT_0436*IT_0475)*IT_0479 + (IT_0308*IT_0466 + IT_0324
      *IT_0469 + IT_0144*IT_0472 + IT_0298*IT_0475)*IT_0483 + (-2)*IT_0485
      *IT_0500 + (IT_0308*IT_0488 + IT_0324*IT_0491 + conj(IT_0432)*IT_0494 +
       IT_0432*IT_0497 + IT_0144*IT_0498 + IT_0298*IT_0499)*IT_0501 + 4*IT_0508*
      (IT_0423*IT_0424 + IT_0425*(IT_0321 + 0.25*IT_0509)) + 6*IT_0048*(IT_0424
      *IT_0425 + IT_0426*(conj(IT_0295) + 0.166666666666667*IT_0510)) + (conj
      (IT_0377)*IT_0378 + conj(IT_0379)*IT_0388 + IT_0379*IT_0390 + conj(IT_0385
      )*IT_0400 + conj(IT_0361)*IT_0418 + IT_0361*IT_0420)*IT_0511 + IT_0513*
      (IT_0299*IT_0415 + IT_0216*conj(IT_0415) + IT_0425*IT_0442 + IT_0426
      *IT_0445 + conj(IT_0399)*IT_0464 + IT_0399*IT_0467 + conj(IT_0295)*IT_0516
       + IT_0295*IT_0519) + IT_0520*IT_0521 + IT_0523*(conj(IT_0410)*IT_0440 +
       conj(IT_0397)*IT_0448 + IT_0397*IT_0450 + conj(IT_0415)*IT_0514 + conj
      (IT_0399)*IT_0526 + IT_0399*IT_0527) + (IT_0434*IT_0442 + IT_0435*IT_0445 
      + conj(IT_0431)*IT_0448 + IT_0431*IT_0451 + IT_0437*IT_0516 + IT_0438
      *IT_0519 + conj(IT_0432)*IT_0526 + 2*IT_0432*(conj(IT_0399) + 0.5*IT_0527)
      )*IT_0531 + IT_0540 + (-8)*IT_0542*(IT_0144*(IT_0431 + IT_0539) + IT_0298*
      (conj(IT_0431) + (-0.125)*IT_0543) + IT_0308*(IT_0432 + (-0.125)*IT_0544) 
      + IT_0324*(conj(IT_0432) + (-0.125)*IT_0545));
    return create_ccomplex_return(IT_0546);
}

