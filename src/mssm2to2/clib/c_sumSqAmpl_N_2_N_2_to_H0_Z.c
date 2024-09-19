#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_H0_Z(
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
    const ccomplex_t IT_0000 = s_14*s_24;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_12*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = pow(m_H0, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_13 + IT_0007 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = cos(alpha);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = sin(alpha);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = sin(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0021;
    const ccomplex_t IT_0026 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0015 + -IT_0019 + 
      -IT_0024 + IT_0027);
    const ccomplex_t IT_0029 = IT_0011*IT_0021;
    const ccomplex_t IT_0030 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0020;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = IT_0029*IT_0034;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0031 + IT_0033 + 
      -IT_0035 + -IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0028*IT_0040;
    const ccomplex_t IT_0042 = IT_0009*IT_0041;
    const ccomplex_t IT_0043 = pow(m_N_2, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0007 + IT_0043 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0047 = IT_0013*IT_0046;
    const ccomplex_t IT_0048 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0049 = IT_0025*IT_0048;
    const ccomplex_t IT_0050 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0051 = IT_0013*IT_0050;
    const ccomplex_t IT_0052 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0053 = IT_0025*IT_0052;
    const ccomplex_t IT_0054 = IT_0047 + IT_0049 + IT_0051 + IT_0053;
    const ccomplex_t IT_0055 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0056 = IT_0017*IT_0055;
    const ccomplex_t IT_0057 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0058 = IT_0022*IT_0057;
    const ccomplex_t IT_0059 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0060 = IT_0017*IT_0059;
    const ccomplex_t IT_0061 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0062 = IT_0022*IT_0061;
    const ccomplex_t IT_0063 = -IT_0056 + -IT_0058 + -IT_0060 + -IT_0062;
    const ccomplex_t IT_0064 = IT_0054 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0068 = IT_0032*IT_0067;
    const ccomplex_t IT_0069 = IT_0029*IT_0067;
    const ccomplex_t IT_0070 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0071 = IT_0032*IT_0070;
    const ccomplex_t IT_0072 = IT_0029*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = IT_0045*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0042 + -IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0039;
    const ccomplex_t IT_0079 = IT_0028*IT_0078;
    const ccomplex_t IT_0080 = IT_0009*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0007 + IT_0043 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0085 = IT_0013*IT_0084;
    const ccomplex_t IT_0086 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0087 = IT_0025*IT_0086;
    const ccomplex_t IT_0088 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0089 = IT_0013*IT_0088;
    const ccomplex_t IT_0090 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0091 = IT_0025*IT_0090;
    const ccomplex_t IT_0092 = IT_0085 + IT_0087 + IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0094 = IT_0017*IT_0093;
    const ccomplex_t IT_0095 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0096 = IT_0022*IT_0095;
    const ccomplex_t IT_0097 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0098 = IT_0017*IT_0097;
    const ccomplex_t IT_0099 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0100 = IT_0022*IT_0099;
    const ccomplex_t IT_0101 = -IT_0094 + -IT_0096 + -IT_0098 + -IT_0100;
    const ccomplex_t IT_0102 = IT_0092 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0106 = IT_0032*IT_0105;
    const ccomplex_t IT_0107 = IT_0029*IT_0105;
    const ccomplex_t IT_0108 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0109 = IT_0032*IT_0108;
    const ccomplex_t IT_0110 = IT_0029*IT_0108;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0106 + IT_0107 + 
      -IT_0109 + -IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0104*IT_0113;
    const ccomplex_t IT_0115 = IT_0083*IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_13 + IT_0007 + IT_0043 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0119 = IT_0013*IT_0118;
    const ccomplex_t IT_0120 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0121 = IT_0025*IT_0120;
    const ccomplex_t IT_0122 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0123 = IT_0013*IT_0122;
    const ccomplex_t IT_0124 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0125 = IT_0025*IT_0124;
    const ccomplex_t IT_0126 = IT_0119 + IT_0121 + IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0128 = IT_0017*IT_0127;
    const ccomplex_t IT_0129 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0130 = IT_0022*IT_0129;
    const ccomplex_t IT_0131 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0132 = IT_0017*IT_0131;
    const ccomplex_t IT_0133 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0134 = IT_0022*IT_0133;
    const ccomplex_t IT_0135 = -IT_0128 + -IT_0130 + -IT_0132 + -IT_0134;
    const ccomplex_t IT_0136 = IT_0126 + IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0140 = IT_0029*IT_0139;
    const ccomplex_t IT_0141 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0142 = IT_0029*IT_0141;
    const ccomplex_t IT_0143 = IT_0032*IT_0139;
    const ccomplex_t IT_0144 = IT_0032*IT_0141;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0140 + -IT_0142 +
       IT_0143 + -IT_0144);
    const ccomplex_t IT_0146 = -IT_0145;
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = IT_0138*IT_0147;
    const ccomplex_t IT_0149 = IT_0117*IT_0148;
    const ccomplex_t IT_0150 = -IT_0115 + -IT_0149;
    const ccomplex_t IT_0151 = conj(IT_0081) + conj(IT_0150);
    const ccomplex_t IT_0152 = -conj(IT_0077) + IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_23 + IT_0007 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0041*IT_0154;
    const ccomplex_t IT_0156 = IT_0079*IT_0154;
    const ccomplex_t IT_0157 = cpow((-2)*s_23 + IT_0007 + IT_0043 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*IT_0157;
    const ccomplex_t IT_0159 = IT_0075*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0155 + (-0.5)*IT_0156 + IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_23 + IT_0007 + IT_0043 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0114*IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_23 + IT_0007 + IT_0043 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = IT_0148*IT_0165;
    const ccomplex_t IT_0167 = IT_0163 + IT_0166;
    const ccomplex_t IT_0168 = IT_0160 + -IT_0167;
    const ccomplex_t IT_0169 = IT_0081 + IT_0150;
    const ccomplex_t IT_0170 = -IT_0077;
    const ccomplex_t IT_0171 = IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = -conj(IT_0167);
    const ccomplex_t IT_0173 = conj(IT_0160) + IT_0172;
    const ccomplex_t IT_0174 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0175 = IT_0013*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0177 = IT_0017*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0179 = IT_0022*IT_0178;
    const ccomplex_t IT_0180 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0181 = IT_0025*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0175 + -IT_0177 + 
      -IT_0179 + IT_0181);
    const ccomplex_t IT_0183 = IT_0078*IT_0182;
    const ccomplex_t IT_0184 = IT_0009*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0186 = IT_0013*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0188 = IT_0025*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0190 = IT_0013*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0192 = IT_0025*IT_0191;
    const ccomplex_t IT_0193 = IT_0186 + IT_0188 + IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0195 = IT_0017*IT_0194;
    const ccomplex_t IT_0196 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0197 = IT_0022*IT_0196;
    const ccomplex_t IT_0198 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0199 = IT_0017*IT_0198;
    const ccomplex_t IT_0200 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0201 = IT_0022*IT_0200;
    const ccomplex_t IT_0202 = -IT_0195 + -IT_0197 + -IT_0199 + -IT_0201;
    const ccomplex_t IT_0203 = IT_0193 + IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0204;
    const ccomplex_t IT_0206 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0207 = IT_0029*IT_0206;
    const ccomplex_t IT_0208 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0209 = IT_0032*IT_0208;
    const ccomplex_t IT_0210 = IT_0029*IT_0208;
    const ccomplex_t IT_0211 = IT_0032*IT_0206;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*(IT_0207 + -IT_0209 + 
      -IT_0210 + IT_0211);
    const ccomplex_t IT_0213 = (-0.5)*IT_0212;
    const ccomplex_t IT_0214 = IT_0205*IT_0213;
    const ccomplex_t IT_0215 = IT_0045*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0184 + -IT_0215;
    const ccomplex_t IT_0217 = 0.5*IT_0037;
    const ccomplex_t IT_0218 = IT_0182*IT_0217;
    const ccomplex_t IT_0219 = IT_0009*IT_0218;
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0222 = IT_0013*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0224 = IT_0025*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0226 = IT_0013*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0228 = IT_0025*IT_0227;
    const ccomplex_t IT_0229 = IT_0222 + IT_0224 + IT_0226 + IT_0228;
    const ccomplex_t IT_0230 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0231 = IT_0017*IT_0230;
    const ccomplex_t IT_0232 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0233 = IT_0022*IT_0232;
    const ccomplex_t IT_0234 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0235 = IT_0017*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0237 = IT_0022*IT_0236;
    const ccomplex_t IT_0238 = -IT_0231 + -IT_0233 + -IT_0235 + -IT_0237;
    const ccomplex_t IT_0239 = IT_0229 + IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = 0.5*IT_0240;
    const ccomplex_t IT_0242 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0243 = IT_0029*IT_0242;
    const ccomplex_t IT_0244 = IT_0032*IT_0242;
    const ccomplex_t IT_0245 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0246 = IT_0032*IT_0245;
    const ccomplex_t IT_0247 = IT_0029*IT_0245;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0243 + IT_0244 + 
      -IT_0246 + -IT_0247);
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = IT_0241*IT_0251;
    const ccomplex_t IT_0253 = IT_0083*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0255 = IT_0013*IT_0254;
    const ccomplex_t IT_0256 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0257 = IT_0025*IT_0256;
    const ccomplex_t IT_0258 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0259 = IT_0013*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0261 = IT_0025*IT_0260;
    const ccomplex_t IT_0262 = IT_0255 + IT_0257 + IT_0259 + IT_0261;
    const ccomplex_t IT_0263 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0264 = IT_0017*IT_0263;
    const ccomplex_t IT_0265 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0266 = IT_0022*IT_0265;
    const ccomplex_t IT_0267 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0268 = IT_0017*IT_0267;
    const ccomplex_t IT_0269 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0270 = IT_0022*IT_0269;
    const ccomplex_t IT_0271 = -IT_0264 + -IT_0266 + -IT_0268 + -IT_0270;
    const ccomplex_t IT_0272 = IT_0262 + IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = 0.5*IT_0273;
    const ccomplex_t IT_0275 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0276 = IT_0029*IT_0275;
    const ccomplex_t IT_0277 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0278 = IT_0032*IT_0277;
    const ccomplex_t IT_0279 = IT_0032*IT_0275;
    const ccomplex_t IT_0280 = IT_0029*IT_0277;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*(IT_0276 + -IT_0278 +
       IT_0279 + -IT_0280);
    const ccomplex_t IT_0282 = -IT_0281;
    const ccomplex_t IT_0283 = -IT_0282;
    const ccomplex_t IT_0284 = (-0.5)*IT_0283;
    const ccomplex_t IT_0285 = IT_0274*IT_0284;
    const ccomplex_t IT_0286 = IT_0117*IT_0285;
    const ccomplex_t IT_0287 = -IT_0253 + -IT_0286;
    const ccomplex_t IT_0288 = conj(IT_0220) + conj(IT_0287);
    const ccomplex_t IT_0289 = IT_0154*IT_0218;
    const ccomplex_t IT_0290 = IT_0154*IT_0183;
    const ccomplex_t IT_0291 = IT_0158*IT_0214;
    const ccomplex_t IT_0292 = (-0.5)*IT_0289 + 0.5*IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = IT_0162*IT_0252;
    const ccomplex_t IT_0294 = IT_0165*IT_0285;
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = IT_0292 + -IT_0295;
    const ccomplex_t IT_0297 = IT_0220 + IT_0287;
    const ccomplex_t IT_0298 = -IT_0216;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = conj(IT_0292) + -conj(IT_0295);
    const ccomplex_t IT_0301 = pow(s_14, 2);
    const ccomplex_t IT_0302 = IT_0001*IT_0301;
    const ccomplex_t IT_0303 = -IT_0043;
    const ccomplex_t IT_0304 = IT_0302 + IT_0303;
    const ccomplex_t IT_0305 = IT_0043*IT_0304;
    const ccomplex_t IT_0306 = cos(beta);
    const ccomplex_t IT_0307 = IT_0012*IT_0306;
    const ccomplex_t IT_0308 = IT_0176*IT_0307;
    const ccomplex_t IT_0309 = sin(beta);
    const ccomplex_t IT_0310 = IT_0012*IT_0309;
    const ccomplex_t IT_0311 = IT_0174*IT_0310;
    const ccomplex_t IT_0312 = IT_0021*IT_0306;
    const ccomplex_t IT_0313 = IT_0180*IT_0312;
    const ccomplex_t IT_0314 = IT_0021*IT_0309;
    const ccomplex_t IT_0315 = IT_0178*IT_0314;
    const ccomplex_t IT_0316 = IT_0308 + -IT_0311 + -IT_0313 + IT_0315;
    const ccomplex_t IT_0317 = -e_em*IT_0012*IT_0021*(IT_0016*IT_0306 + 
      -IT_0010*IT_0309);
    const ccomplex_t IT_0318 = IT_0316*IT_0317;
    const ccomplex_t IT_0319 = cpow((-2)*s_12 + (-2)*IT_0043 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = IT_0318*IT_0320;
    const ccomplex_t IT_0322 = (-0.5)*IT_0248;
    const ccomplex_t IT_0323 = IT_0241*IT_0322;
    const ccomplex_t IT_0324 = IT_0162*IT_0323;
    const ccomplex_t IT_0325 = IT_0174*IT_0307;
    const ccomplex_t IT_0326 = IT_0176*IT_0310;
    const ccomplex_t IT_0327 = IT_0178*IT_0312;
    const ccomplex_t IT_0328 = IT_0180*IT_0314;
    const ccomplex_t IT_0329 = IT_0325 + IT_0326 + -IT_0327 + -IT_0328;
    const ccomplex_t IT_0330 = e_em*IT_0012*IT_0021*(IT_0010*IT_0306 + IT_0016
      *IT_0309);
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = IT_0329*IT_0331;
    const ccomplex_t IT_0333 = pow(m_Z, 2);
    const ccomplex_t IT_0334 = cpow((-2)*s_12 + (-2)*IT_0043 + IT_0333 + 
      -reg_prop, -1);
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0332*IT_0335;
    const ccomplex_t IT_0337 = (-2)*IT_0294 + (-0.5)*IT_0321 + (-2)*IT_0324 + 
      (-0.5)*IT_0336;
    const ccomplex_t IT_0338 = IT_0018*IT_0307;
    const ccomplex_t IT_0339 = IT_0014*IT_0310;
    const ccomplex_t IT_0340 = IT_0026*IT_0312;
    const ccomplex_t IT_0341 = IT_0023*IT_0314;
    const ccomplex_t IT_0342 = -IT_0338 + IT_0339 + IT_0340 + -IT_0341;
    const ccomplex_t IT_0343 = IT_0317*IT_0342;
    const ccomplex_t IT_0344 = IT_0320*IT_0343;
    const ccomplex_t IT_0345 = IT_0028*IT_0217;
    const ccomplex_t IT_0346 = IT_0154*IT_0345;
    const ccomplex_t IT_0347 = IT_0014*IT_0307;
    const ccomplex_t IT_0348 = IT_0018*IT_0310;
    const ccomplex_t IT_0349 = IT_0023*IT_0312;
    const ccomplex_t IT_0350 = IT_0026*IT_0314;
    const ccomplex_t IT_0351 = -IT_0347 + -IT_0348 + IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = IT_0331*IT_0351;
    const ccomplex_t IT_0353 = IT_0335*IT_0352;
    const ccomplex_t IT_0354 = IT_0156 + (-2)*IT_0159 + 0.5*IT_0344 + -IT_0346
       + 0.5*IT_0353;
    const ccomplex_t IT_0355 = 2*IT_0354;
    const ccomplex_t IT_0356 = (-2)*IT_0163 + (-2)*IT_0166 + (-0.5)*IT_0344 + 
      (-0.5)*IT_0353;
    const ccomplex_t IT_0357 = (-2)*IT_0356;
    const ccomplex_t IT_0358 = IT_0355 + IT_0357;
    const ccomplex_t IT_0359 = 2*conj(IT_0354);
    const ccomplex_t IT_0360 = (-2)*conj(IT_0356);
    const ccomplex_t IT_0361 = IT_0359 + IT_0360;
    const ccomplex_t IT_0362 = IT_0289 + -IT_0290 + (-2)*IT_0291 + 0.5*IT_0321
       + 0.5*IT_0336;
    const ccomplex_t IT_0363 = 2*IT_0356;
    const ccomplex_t IT_0364 = (-2)*IT_0354;
    const ccomplex_t IT_0365 = IT_0363 + IT_0364;
    const ccomplex_t IT_0366 = 2*conj(IT_0356);
    const ccomplex_t IT_0367 = (-2)*conj(IT_0354);
    const ccomplex_t IT_0368 = IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0370 = -IT_0369;
    const ccomplex_t IT_0371 = s_12 + IT_0370;
    const ccomplex_t IT_0372 = s_12*IT_0371;
    const ccomplex_t IT_0373 = IT_0009*IT_0345;
    const ccomplex_t IT_0374 = 2*IT_0076 + (-0.5)*IT_0344 + (-0.5)*IT_0353 +
       IT_0373;
    const ccomplex_t IT_0375 = IT_0080 + 2*IT_0115 + 2*IT_0149 + 0.5*IT_0344 +
       0.5*IT_0353;
    const ccomplex_t IT_0376 = IT_0184 + 2*IT_0215 + (-0.5)*IT_0321 + (-0.5)
      *IT_0336;
    const ccomplex_t IT_0377 = 2*IT_0362;
    const ccomplex_t IT_0378 = (-2)*IT_0337;
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = 2*conj(IT_0362);
    const ccomplex_t IT_0381 = (-2)*conj(IT_0337);
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = IT_0083*IT_0323;
    const ccomplex_t IT_0384 = IT_0219 + 2*IT_0286 + 0.5*IT_0321 + 0.5*IT_0336
       + 2*IT_0383;
    const ccomplex_t IT_0385 = 2*IT_0337;
    const ccomplex_t IT_0386 = (-2)*IT_0362;
    const ccomplex_t IT_0387 = IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = 2*conj(IT_0337);
    const ccomplex_t IT_0389 = (-2)*conj(IT_0362);
    const ccomplex_t IT_0390 = IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = s_12*s_14;
    const ccomplex_t IT_0392 = s_24*IT_0043;
    const ccomplex_t IT_0393 = -IT_0392;
    const ccomplex_t IT_0394 = IT_0391 + IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0077) + conj(IT_0160);
    const ccomplex_t IT_0396 = IT_0077 + IT_0160;
    const ccomplex_t IT_0397 = conj(IT_0216) + conj(IT_0292);
    const ccomplex_t IT_0398 = IT_0216 + IT_0292;
    const ccomplex_t IT_0399 = IT_0394*((IT_0151 + conj(IT_0167))*IT_0365 + 
      (IT_0167 + IT_0169)*IT_0368 + (IT_0288 + conj(IT_0295))*IT_0387 + (IT_0295
       + IT_0297)*IT_0390 + IT_0358*IT_0395 + IT_0361*IT_0396 + IT_0379*IT_0397 
      + IT_0382*IT_0398);
    const ccomplex_t IT_0400 = IT_0043*IT_0371;
    const ccomplex_t IT_0401 = 2*IT_0375;
    const ccomplex_t IT_0402 = (-2)*IT_0374;
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = 2*conj(IT_0375);
    const ccomplex_t IT_0405 = (-2)*conj(IT_0374);
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = 2*IT_0374;
    const ccomplex_t IT_0408 = (-2)*IT_0375;
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = 2*conj(IT_0374);
    const ccomplex_t IT_0411 = (-2)*conj(IT_0375);
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = pow(s_24, 2);
    const ccomplex_t IT_0414 = IT_0001*IT_0413;
    const ccomplex_t IT_0415 = -IT_0414;
    const ccomplex_t IT_0416 = IT_0043 + IT_0415;
    const ccomplex_t IT_0417 = IT_0043*IT_0416;
    const ccomplex_t IT_0418 = -IT_0417;
    const ccomplex_t IT_0419 = m_N_2*IT_0416;
    const ccomplex_t IT_0420 = -IT_0419;
    const ccomplex_t IT_0421 = cpow(IT_0011, -4);
    const ccomplex_t IT_0422 = IT_0010*IT_0020*IT_0306*IT_0421;
    const ccomplex_t IT_0423 = IT_0016*IT_0020*IT_0309*IT_0421;
    const ccomplex_t IT_0424 = cpow(IT_0020, -3);
    const ccomplex_t IT_0425 = IT_0010*IT_0306*IT_0424;
    const ccomplex_t IT_0426 = IT_0016*IT_0309*IT_0424;
    const ccomplex_t IT_0427 = cpow(IT_0011, -2);
    const ccomplex_t IT_0428 = IT_0010*IT_0021*IT_0306*IT_0427;
    const ccomplex_t IT_0429 = IT_0016*IT_0021*IT_0309*IT_0427;
    const ccomplex_t IT_0430 = cpow(IT_0020, -2);
    const ccomplex_t IT_0431 = cpow(IT_0427 + IT_0430, -1);
    const ccomplex_t IT_0432 = (IT_0422 + IT_0423 + IT_0425 + IT_0426 + 2
      *IT_0428 + 2*IT_0429)*IT_0431;
    const ccomplex_t IT_0433 = m_W*e_em;
    const ccomplex_t IT_0434 = IT_0432*IT_0433;
    const ccomplex_t IT_0435 = (0 + _Complex_I*1)*IT_0434;
    const ccomplex_t IT_0436 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0043 + -reg_prop, -1);
    const ccomplex_t IT_0437 = IT_0435*IT_0436;
    const ccomplex_t IT_0438 = (0 + _Complex_I*1)*IT_0437;
    const ccomplex_t IT_0439 = IT_0078*IT_0438;
    const ccomplex_t IT_0440 = m_N_2*IT_0079;
    const ccomplex_t IT_0441 = m_N_2*IT_0345;
    const ccomplex_t IT_0442 = -IT_0441;
    const ccomplex_t IT_0443 = IT_0440 + IT_0442;
    const ccomplex_t IT_0444 = IT_0154*IT_0443;
    const ccomplex_t IT_0445 = IT_0154*IT_0440;
    const ccomplex_t IT_0446 = m_N_2*IT_0218;
    const ccomplex_t IT_0447 = m_N_2*IT_0183;
    const ccomplex_t IT_0448 = -IT_0447;
    const ccomplex_t IT_0449 = IT_0446 + IT_0448;
    const ccomplex_t IT_0450 = IT_0009*IT_0449;
    const ccomplex_t IT_0451 = IT_0066*IT_0213;
    const ccomplex_t IT_0452 = m_N_1*IT_0451;
    const ccomplex_t IT_0453 = m_N_2*IT_0075;
    const ccomplex_t IT_0454 = -IT_0453;
    const ccomplex_t IT_0455 = IT_0452 + IT_0454;
    const ccomplex_t IT_0456 = IT_0158*IT_0455;
    const ccomplex_t IT_0457 = IT_0113*IT_0241;
    const ccomplex_t IT_0458 = m_N_3*IT_0457;
    const ccomplex_t IT_0459 = m_N_2*IT_0323;
    const ccomplex_t IT_0460 = -IT_0459;
    const ccomplex_t IT_0461 = IT_0458 + IT_0460;
    const ccomplex_t IT_0462 = IT_0083*IT_0461;
    const ccomplex_t IT_0463 = IT_0147*IT_0274;
    const ccomplex_t IT_0464 = m_N_4*IT_0463;
    const ccomplex_t IT_0465 = m_N_2*IT_0285;
    const ccomplex_t IT_0466 = -IT_0465;
    const ccomplex_t IT_0467 = IT_0464 + IT_0466;
    const ccomplex_t IT_0468 = IT_0117*IT_0467;
    const ccomplex_t IT_0469 = IT_0439 + 0.5*IT_0444 + 0.5*IT_0445 + (-0.5)
      *IT_0450 + IT_0456 + IT_0462 + IT_0468;
    const ccomplex_t IT_0470 = IT_0154*IT_0441;
    const ccomplex_t IT_0471 = IT_0074*IT_0205;
    const ccomplex_t IT_0472 = m_N_1*IT_0471;
    const ccomplex_t IT_0473 = m_N_2*IT_0214;
    const ccomplex_t IT_0474 = -IT_0473;
    const ccomplex_t IT_0475 = IT_0472 + IT_0474;
    const ccomplex_t IT_0476 = IT_0045*IT_0475;
    const ccomplex_t IT_0477 = IT_0104*IT_0322;
    const ccomplex_t IT_0478 = m_N_3*IT_0477;
    const ccomplex_t IT_0479 = m_N_2*IT_0114;
    const ccomplex_t IT_0480 = -IT_0479;
    const ccomplex_t IT_0481 = IT_0478 + IT_0480;
    const ccomplex_t IT_0482 = IT_0162*IT_0481;
    const ccomplex_t IT_0483 = IT_0138*IT_0284;
    const ccomplex_t IT_0484 = m_N_4*IT_0483;
    const ccomplex_t IT_0485 = m_N_2*IT_0148;
    const ccomplex_t IT_0486 = -IT_0485;
    const ccomplex_t IT_0487 = IT_0484 + IT_0486;
    const ccomplex_t IT_0488 = IT_0165*IT_0487;
    const ccomplex_t IT_0489 = (-0.5)*IT_0450 + (-0.5)*IT_0470 + -IT_0476 + 
      -IT_0482 + -IT_0488;
    const ccomplex_t IT_0490 = conj(IT_0469) + conj(IT_0489);
    const ccomplex_t IT_0491 = IT_0469 + IT_0489;
    const ccomplex_t IT_0492 = IT_0217*IT_0438;
    const ccomplex_t IT_0493 = IT_0154*IT_0449;
    const ccomplex_t IT_0494 = IT_0154*IT_0446;
    const ccomplex_t IT_0495 = IT_0009*IT_0443;
    const ccomplex_t IT_0496 = IT_0158*IT_0475;
    const ccomplex_t IT_0497 = IT_0083*IT_0481;
    const ccomplex_t IT_0498 = IT_0117*IT_0487;
    const ccomplex_t IT_0499 = IT_0492 + 0.5*IT_0493 + 0.5*IT_0494 + (-0.5)
      *IT_0495 + IT_0496 + IT_0497 + IT_0498;
    const ccomplex_t IT_0500 = IT_0154*IT_0447;
    const ccomplex_t IT_0501 = IT_0045*IT_0455;
    const ccomplex_t IT_0502 = IT_0162*IT_0461;
    const ccomplex_t IT_0503 = IT_0165*IT_0467;
    const ccomplex_t IT_0504 = (-0.5)*IT_0495 + (-0.5)*IT_0500 + -IT_0501 + 
      -IT_0502 + -IT_0503;
    const ccomplex_t IT_0505 = conj(IT_0499) + conj(IT_0504);
    const ccomplex_t IT_0506 = 2*IT_0376;
    const ccomplex_t IT_0507 = (-2)*IT_0384;
    const ccomplex_t IT_0508 = IT_0506 + IT_0507;
    const ccomplex_t IT_0509 = IT_0499 + IT_0504;
    const ccomplex_t IT_0510 = 2*conj(IT_0376);
    const ccomplex_t IT_0511 = (-2)*conj(IT_0384);
    const ccomplex_t IT_0512 = IT_0510 + IT_0511;
    const ccomplex_t IT_0513 = s_24*IT_0333;
    const ccomplex_t IT_0514 = IT_0001*IT_0513;
    const ccomplex_t IT_0515 = (-0.25)*IT_0514;
    const ccomplex_t IT_0516 = s_24 + IT_0515;
    const ccomplex_t IT_0517 = m_N_2*IT_0516;
    const ccomplex_t IT_0518 = (-4)*IT_0517;
    const ccomplex_t IT_0519 = 2*IT_0469;
    const ccomplex_t IT_0520 = 2*IT_0489;
    const ccomplex_t IT_0521 = IT_0519 + IT_0520;
    const ccomplex_t IT_0522 = 2*conj(IT_0469);
    const ccomplex_t IT_0523 = 2*conj(IT_0489);
    const ccomplex_t IT_0524 = IT_0522 + IT_0523;
    const ccomplex_t IT_0525 = 2*IT_0499;
    const ccomplex_t IT_0526 = 2*IT_0504;
    const ccomplex_t IT_0527 = IT_0525 + IT_0526;
    const ccomplex_t IT_0528 = 2*conj(IT_0499);
    const ccomplex_t IT_0529 = 2*conj(IT_0504);
    const ccomplex_t IT_0530 = IT_0528 + IT_0529;
    const ccomplex_t IT_0531 = (-2)*IT_0469;
    const ccomplex_t IT_0532 = (-2)*IT_0489;
    const ccomplex_t IT_0533 = IT_0531 + IT_0532;
    const ccomplex_t IT_0534 = (-2)*conj(IT_0469);
    const ccomplex_t IT_0535 = (-2)*conj(IT_0489);
    const ccomplex_t IT_0536 = IT_0534 + IT_0535;
    const ccomplex_t IT_0537 = (-2)*IT_0499;
    const ccomplex_t IT_0538 = (-2)*IT_0504;
    const ccomplex_t IT_0539 = IT_0537 + IT_0538;
    const ccomplex_t IT_0540 = (-2)*conj(IT_0499);
    const ccomplex_t IT_0541 = (-2)*conj(IT_0504);
    const ccomplex_t IT_0542 = IT_0540 + IT_0541;
    const ccomplex_t IT_0543 = s_14*IT_0333;
    const ccomplex_t IT_0544 = IT_0001*IT_0543;
    const ccomplex_t IT_0545 = 0.5*IT_0544;
    const ccomplex_t IT_0546 = s_14 + IT_0545;
    const ccomplex_t IT_0547 = m_N_2*IT_0546;
    const ccomplex_t IT_0548 = 2*IT_0547;
    const ccomplex_t IT_0549 = (-0.25)*IT_0544;
    const ccomplex_t IT_0550 = s_14 + IT_0549;
    const ccomplex_t IT_0551 = m_N_2*IT_0550;
    const ccomplex_t IT_0552 = (-4)*IT_0551;
    const ccomplex_t IT_0553 = 0.5*IT_0514;
    const ccomplex_t IT_0554 = s_24 + IT_0553;
    const ccomplex_t IT_0555 = m_N_2*IT_0554;
    const ccomplex_t IT_0556 = 2*IT_0555;
    const ccomplex_t IT_0557 = s_12*IT_0333;
    const ccomplex_t IT_0558 = (-2)*IT_0000;
    const ccomplex_t IT_0559 = IT_0557 + IT_0558;
    const ccomplex_t IT_0560 = IT_0001*IT_0559;
    const ccomplex_t IT_0561 = (-0.5)*IT_0560;
    const ccomplex_t IT_0562 = s_12 + IT_0561;
    const ccomplex_t IT_0563 = 4*IT_0499;
    const ccomplex_t IT_0564 = 6*conj(IT_0469);
    const ccomplex_t IT_0565 = m_N_2*IT_0304;
    const ccomplex_t IT_0566 = s_12*IT_0304;
    const ccomplex_t IT_0567 = m_N_2*IT_0371;
    const ccomplex_t IT_0568 = IT_0361*IT_0469 + IT_0358*conj(IT_0469) + 
      (IT_0377 + IT_0378 + IT_0407 + IT_0408)*IT_0505 + conj(IT_0489)*(IT_0355 +
       IT_0357 + IT_0506 + IT_0507) + (IT_0380 + IT_0381 + IT_0410 + IT_0411)
      *IT_0509 + IT_0489*(IT_0359 + IT_0360 + IT_0510 + IT_0511) + conj(IT_0376)
      *IT_0519 + IT_0376*IT_0522 + conj(IT_0384)*IT_0531 + IT_0384*IT_0534;
    const ccomplex_t IT_0569 = s_12*IT_0416;
    const ccomplex_t IT_0570 = -IT_0569;
    const ccomplex_t IT_0571 = 2*IT_0384;
    const ccomplex_t IT_0572 = s_12*s_24;
    const ccomplex_t IT_0573 = s_14*IT_0043;
    const ccomplex_t IT_0574 = -IT_0573;
    const ccomplex_t IT_0575 = IT_0572 + IT_0574;
    const ccomplex_t IT_0576 = (-2)*IT_0376;
    const ccomplex_t IT_0577 = 2*conj(IT_0384);
    const ccomplex_t IT_0578 = (-2)*conj(IT_0376);
    const ccomplex_t IT_0579 = (IT_0151 + conj(IT_0167))*IT_0403 + (IT_0167 +
       IT_0169)*IT_0406 + IT_0395*IT_0409 + IT_0396*IT_0412 + IT_0397*IT_0508 +
       IT_0398*IT_0512 + (IT_0288 + conj(IT_0295))*(IT_0571 + IT_0576) + 
      (IT_0295 + IT_0297)*(IT_0577 + IT_0578);
    const ccomplex_t IT_0580 = IT_0043*IT_0333;
    const ccomplex_t IT_0581 = (-3)*IT_0580;
    const ccomplex_t IT_0582 = 2*IT_0081;
    const ccomplex_t IT_0583 = 2*IT_0150;
    const ccomplex_t IT_0584 = (-2)*IT_0077;
    const ccomplex_t IT_0585 = IT_0582 + IT_0583 + IT_0584;
    const ccomplex_t IT_0586 = 2*conj(IT_0081);
    const ccomplex_t IT_0587 = 2*conj(IT_0150);
    const ccomplex_t IT_0588 = (-2)*conj(IT_0077);
    const ccomplex_t IT_0589 = IT_0586 + IT_0587 + IT_0588;
    const ccomplex_t IT_0590 = 2*IT_0160;
    const ccomplex_t IT_0591 = (-2)*IT_0167;
    const ccomplex_t IT_0592 = IT_0590 + IT_0591;
    const ccomplex_t IT_0593 = 2*conj(IT_0160);
    const ccomplex_t IT_0594 = (-2)*conj(IT_0167);
    const ccomplex_t IT_0595 = IT_0593 + IT_0594;
    const ccomplex_t IT_0596 = 2*IT_0220;
    const ccomplex_t IT_0597 = 2*IT_0287;
    const ccomplex_t IT_0598 = IT_0596 + IT_0597;
    const ccomplex_t IT_0599 = 2*conj(IT_0220);
    const ccomplex_t IT_0600 = 2*conj(IT_0287);
    const ccomplex_t IT_0601 = IT_0599 + IT_0600;
    const ccomplex_t IT_0602 = (-2)*IT_0081;
    const ccomplex_t IT_0603 = (-2)*IT_0150;
    const ccomplex_t IT_0604 = IT_0602 + IT_0603;
    const ccomplex_t IT_0605 = (-2)*conj(IT_0081);
    const ccomplex_t IT_0606 = (-2)*conj(IT_0150);
    const ccomplex_t IT_0607 = IT_0605 + IT_0606;
    const ccomplex_t IT_0608 = 2*IT_0292;
    const ccomplex_t IT_0609 = 2*conj(IT_0292);
    const ccomplex_t IT_0610 = (-2)*IT_0160;
    const ccomplex_t IT_0611 = (-2)*conj(IT_0160);
    const ccomplex_t IT_0612 = conj(IT_0216)*IT_0585 + IT_0216*IT_0589 + conj
      (IT_0295)*IT_0592 + IT_0295*IT_0595 + conj(IT_0077)*IT_0598 + IT_0077
      *IT_0601 + IT_0288*IT_0604 + IT_0297*IT_0607 + conj(IT_0167)*IT_0608 +
       IT_0167*IT_0609 + conj(IT_0292)*IT_0610 + IT_0292*IT_0611;
    const ccomplex_t IT_0613 = 3*IT_0580;
    const ccomplex_t IT_0614 = conj(IT_0295)*IT_0585 + IT_0295*IT_0589 + conj
      (IT_0216)*IT_0592 + IT_0216*IT_0595 + conj(IT_0167)*IT_0598 + IT_0167
      *IT_0601 + conj(IT_0292)*IT_0604 + IT_0292*IT_0607 + conj(IT_0077)*IT_0608
       + IT_0077*IT_0609 + IT_0288*IT_0610 + IT_0297*IT_0611;
    const ccomplex_t IT_0615 = (-0.5)*IT_0004;
    const ccomplex_t IT_0616 = IT_0557 + IT_0558 + IT_0615;
    const ccomplex_t IT_0617 = (-4)*IT_0220;
    const ccomplex_t IT_0618 = (-8)*IT_0006*(IT_0152*IT_0168 + IT_0171*IT_0173
       + -(conj(IT_0216) + -IT_0288)*IT_0296 + IT_0299*IT_0300) + IT_0305*(conj
      (IT_0337)*IT_0358 + IT_0337*IT_0361 + conj(IT_0362)*IT_0365 + IT_0362
      *IT_0368) + IT_0372*(IT_0361*IT_0374 + IT_0358*conj(IT_0374) + IT_0368
      *IT_0375 + IT_0365*conj(IT_0375) + conj(IT_0376)*IT_0379 + IT_0376*IT_0382
       + conj(IT_0384)*IT_0387 + IT_0384*IT_0390) + IT_0399 + IT_0400*(IT_0368
      *IT_0376 + IT_0365*conj(IT_0376) + IT_0361*IT_0384 + IT_0358*conj(IT_0384)
       + conj(IT_0362)*IT_0403 + IT_0362*IT_0406 + conj(IT_0337)*IT_0409 +
       IT_0337*IT_0412) + (conj(IT_0376)*IT_0403 + IT_0376*IT_0406 + conj
      (IT_0384)*IT_0409 + IT_0384*IT_0412)*IT_0418 + IT_0420*(IT_0409*IT_0490 +
       IT_0412*IT_0491 + IT_0505*IT_0508 + IT_0509*IT_0512) + IT_0518*(conj
      (IT_0167)*IT_0521 + IT_0167*IT_0524 + conj(IT_0295)*IT_0527 + IT_0295
      *IT_0530 + conj(IT_0160)*IT_0533 + IT_0160*IT_0536 + conj(IT_0292)*IT_0539
       + IT_0292*IT_0542) + (conj(IT_0292)*IT_0521 + IT_0292*IT_0524 + conj
      (IT_0160)*IT_0527 + IT_0160*IT_0530 + conj(IT_0295)*IT_0533 + IT_0295
      *IT_0536 + conj(IT_0167)*IT_0539 + IT_0167*IT_0542)*IT_0548 + (IT_0288
      *IT_0521 + IT_0297*IT_0524 + IT_0151*IT_0527 + IT_0169*IT_0530 + conj
      (IT_0216)*IT_0533 + IT_0216*IT_0536 + conj(IT_0077)*IT_0539 + IT_0077
      *IT_0542)*IT_0552 + (conj(IT_0077)*IT_0521 + IT_0077*IT_0524 + conj
      (IT_0216)*IT_0527 + IT_0216*IT_0530 + IT_0151*IT_0533 + IT_0169*IT_0536 +
       IT_0288*IT_0539 + IT_0297*IT_0542)*IT_0556 + 4*IT_0562*(IT_0490*IT_0491 +
       IT_0505*(IT_0504 + 0.25*IT_0563)) + 6*IT_0043*(IT_0491*IT_0505 + IT_0509*
      (conj(IT_0489) + 0.166666666666667*IT_0564)) + (IT_0382*IT_0489 + IT_0379
      *conj(IT_0489) + IT_0358*IT_0505 + IT_0361*IT_0509 + conj(IT_0362)*IT_0519
       + IT_0362*IT_0522 + conj(IT_0337)*IT_0531 + IT_0337*IT_0534)*IT_0565 + 
      (conj(IT_0354)*IT_0355 + conj(IT_0356)*IT_0365 + IT_0356*IT_0367 + conj
      (IT_0362)*IT_0377 + conj(IT_0337)*IT_0387 + IT_0337*IT_0389)*IT_0566 +
       IT_0567*IT_0568 + IT_0570*(conj(IT_0375)*IT_0401 + conj(IT_0374)*IT_0409 
      + IT_0374*IT_0411 + conj(IT_0376)*IT_0508 + IT_0376*IT_0511 + conj(IT_0384
      )*IT_0571) + IT_0575*IT_0579 + IT_0581*IT_0612 + IT_0613*IT_0614 + (-4)
      *IT_0616*(IT_0152*IT_0171 + IT_0168*IT_0173 + -conj(IT_0216)*IT_0299 +
       IT_0296*IT_0300 + IT_0288*(IT_0287 + IT_0298 + (-0.25)*IT_0617));
    return create_ccomplex_return(IT_0618);
}

