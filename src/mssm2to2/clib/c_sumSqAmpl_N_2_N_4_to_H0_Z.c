#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_4_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_4_to_H0_Z(
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
    const ccomplex_t IT_0000 = m_N_2*m_N_4;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = pow(m_H0, 2);
    const ccomplex_t IT_0009 = pow(m_N_4, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = cos(alpha);
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(alpha);
    const ccomplex_t IT_0019 = sin(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0027 = IT_0021*IT_0026;
    const ccomplex_t IT_0028 = IT_0017 + IT_0023 + IT_0025 + IT_0027;
    const ccomplex_t IT_0029 = IT_0014*IT_0018;
    const ccomplex_t IT_0030 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0020;
    const ccomplex_t IT_0033 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0036 = IT_0029*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = -IT_0031 + -IT_0034 + -IT_0036 + -IT_0038;
    const ccomplex_t IT_0040 = IT_0028 + IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0013*IT_0020;
    const ccomplex_t IT_0044 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = IT_0014*IT_0019;
    const ccomplex_t IT_0047 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = IT_0043*IT_0047;
    const ccomplex_t IT_0050 = IT_0044*IT_0046;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0045 + -IT_0048 + 
      -IT_0049 + IT_0050);
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0042*IT_0052;
    const ccomplex_t IT_0054 = IT_0011*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0008 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0058 = IT_0015*IT_0057;
    const ccomplex_t IT_0059 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0060 = IT_0029*IT_0059;
    const ccomplex_t IT_0061 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0062 = IT_0032*IT_0061;
    const ccomplex_t IT_0063 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0064 = IT_0021*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0058 + -IT_0060 + 
      -IT_0062 + IT_0064);
    const ccomplex_t IT_0066 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0067 = IT_0046*IT_0066;
    const ccomplex_t IT_0068 = IT_0043*IT_0066;
    const ccomplex_t IT_0069 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0070 = IT_0043*IT_0069;
    const ccomplex_t IT_0071 = IT_0046*IT_0069;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0067 + IT_0068 + 
      -IT_0070 + -IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0065*IT_0076;
    const ccomplex_t IT_0078 = IT_0056*IT_0077;
    const ccomplex_t IT_0079 = (-2)*IT_0054 + (-2)*IT_0078;
    const ccomplex_t IT_0080 = cos(beta);
    const ccomplex_t IT_0081 = IT_0014*IT_0080;
    const ccomplex_t IT_0082 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = sin(beta);
    const ccomplex_t IT_0085 = IT_0014*IT_0084;
    const ccomplex_t IT_0086 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = IT_0020*IT_0080;
    const ccomplex_t IT_0089 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = IT_0020*IT_0084;
    const ccomplex_t IT_0092 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0095 = IT_0081*IT_0094;
    const ccomplex_t IT_0096 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0097 = IT_0085*IT_0096;
    const ccomplex_t IT_0098 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0099 = IT_0088*IT_0098;
    const ccomplex_t IT_0100 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0101 = IT_0091*IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0083 + 0.5*IT_0087 + 0.5*IT_0090 + (
      -0.5)*IT_0093 + (-0.5)*IT_0095 + 0.5*IT_0097 + 0.5*IT_0099 + (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = -e_em*IT_0014*IT_0020*(IT_0018*IT_0080 + 
      -IT_0012*IT_0084);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + -IT_0004 + -IT_0009 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0015*IT_0086;
    const ccomplex_t IT_0111 = IT_0021*IT_0089;
    const ccomplex_t IT_0112 = IT_0015*IT_0096;
    const ccomplex_t IT_0113 = IT_0021*IT_0098;
    const ccomplex_t IT_0114 = IT_0110 + IT_0111 + IT_0112 + IT_0113;
    const ccomplex_t IT_0115 = IT_0029*IT_0082;
    const ccomplex_t IT_0116 = IT_0032*IT_0092;
    const ccomplex_t IT_0117 = IT_0029*IT_0094;
    const ccomplex_t IT_0118 = IT_0032*IT_0100;
    const ccomplex_t IT_0119 = -IT_0115 + -IT_0116 + -IT_0117 + -IT_0118;
    const ccomplex_t IT_0120 = IT_0114 + IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0124 = IT_0043*IT_0123;
    const ccomplex_t IT_0125 = IT_0046*IT_0123;
    const ccomplex_t IT_0126 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0127 = IT_0046*IT_0126;
    const ccomplex_t IT_0128 = IT_0043*IT_0126;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0124 + IT_0125 + 
      -IT_0127 + -IT_0128);
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = IT_0122*IT_0130;
    const ccomplex_t IT_0132 = IT_0109*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0129;
    const ccomplex_t IT_0134 = IT_0122*IT_0133;
    const ccomplex_t IT_0135 = IT_0109*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0137 = IT_0046*IT_0136;
    const ccomplex_t IT_0138 = IT_0043*IT_0136;
    const ccomplex_t IT_0139 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0140 = IT_0046*IT_0139;
    const ccomplex_t IT_0141 = IT_0043*IT_0139;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0137 + IT_0138 + 
      -IT_0140 + -IT_0141);
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0145 = IT_0015*IT_0144;
    const ccomplex_t IT_0146 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0147 = IT_0021*IT_0146;
    const ccomplex_t IT_0148 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0149 = IT_0015*IT_0148;
    const ccomplex_t IT_0150 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0151 = IT_0021*IT_0150;
    const ccomplex_t IT_0152 = IT_0145 + IT_0147 + IT_0149 + IT_0151;
    const ccomplex_t IT_0153 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0154 = IT_0029*IT_0153;
    const ccomplex_t IT_0155 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0156 = IT_0032*IT_0155;
    const ccomplex_t IT_0157 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0158 = IT_0029*IT_0157;
    const ccomplex_t IT_0159 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0160 = IT_0032*IT_0159;
    const ccomplex_t IT_0161 = -IT_0154 + -IT_0156 + -IT_0158 + -IT_0160;
    const ccomplex_t IT_0162 = IT_0152 + IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = 0.5*IT_0163;
    const ccomplex_t IT_0165 = IT_0143*IT_0164;
    const ccomplex_t IT_0166 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = IT_0165*IT_0167;
    const ccomplex_t IT_0169 = pow(m_Z, 2);
    const ccomplex_t IT_0170 = cpow((-2)*s_12 + -IT_0004 + -IT_0009 + IT_0169 
      + -reg_prop, -1);
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = e_em*IT_0014*IT_0020*(IT_0012*IT_0080 + IT_0018
      *IT_0084);
    const ccomplex_t IT_0173 = -IT_0172;
    const ccomplex_t IT_0174 = IT_0081*IT_0086;
    const ccomplex_t IT_0175 = IT_0082*IT_0085;
    const ccomplex_t IT_0176 = IT_0088*IT_0092;
    const ccomplex_t IT_0177 = IT_0089*IT_0091;
    const ccomplex_t IT_0178 = IT_0081*IT_0096;
    const ccomplex_t IT_0179 = IT_0085*IT_0094;
    const ccomplex_t IT_0180 = IT_0088*IT_0100;
    const ccomplex_t IT_0181 = IT_0091*IT_0098;
    const ccomplex_t IT_0182 = (-0.5)*IT_0174 + (-0.5)*IT_0175 + 0.5*IT_0176 +
       0.5*IT_0177 + (-0.5)*IT_0178 + (-0.5)*IT_0179 + 0.5*IT_0180 + 0.5*IT_0181;
    const ccomplex_t IT_0183 = IT_0173*IT_0182;
    const ccomplex_t IT_0184 = IT_0171*IT_0183;
    const ccomplex_t IT_0185 = IT_0107 + IT_0132 + -IT_0135 + (-2)*IT_0168 +
       IT_0184;
    const ccomplex_t IT_0186 = 2*IT_0185;
    const ccomplex_t IT_0187 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0188 = IT_0015*IT_0187;
    const ccomplex_t IT_0189 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0190 = IT_0021*IT_0189;
    const ccomplex_t IT_0191 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0192 = IT_0015*IT_0191;
    const ccomplex_t IT_0193 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0194 = IT_0021*IT_0193;
    const ccomplex_t IT_0195 = IT_0188 + IT_0190 + IT_0192 + IT_0194;
    const ccomplex_t IT_0196 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0197 = IT_0029*IT_0196;
    const ccomplex_t IT_0198 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0199 = IT_0032*IT_0198;
    const ccomplex_t IT_0200 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0201 = IT_0029*IT_0200;
    const ccomplex_t IT_0202 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0203 = IT_0032*IT_0202;
    const ccomplex_t IT_0204 = -IT_0197 + -IT_0199 + -IT_0201 + -IT_0203;
    const ccomplex_t IT_0205 = IT_0195 + IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = 0.5*IT_0206;
    const ccomplex_t IT_0208 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0209 = IT_0046*IT_0208;
    const ccomplex_t IT_0210 = IT_0043*IT_0208;
    const ccomplex_t IT_0211 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0212 = IT_0046*IT_0211;
    const ccomplex_t IT_0213 = IT_0043*IT_0211;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*(IT_0209 + IT_0210 + 
      -IT_0212 + -IT_0213);
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0207*IT_0216;
    const ccomplex_t IT_0218 = IT_0011*IT_0217;
    const ccomplex_t IT_0219 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0220 = IT_0015*IT_0219;
    const ccomplex_t IT_0221 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0222 = IT_0029*IT_0221;
    const ccomplex_t IT_0223 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0224 = IT_0032*IT_0223;
    const ccomplex_t IT_0225 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0226 = IT_0021*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0220 + -IT_0222 + 
      -IT_0224 + IT_0226);
    const ccomplex_t IT_0228 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0229 = IT_0046*IT_0228;
    const ccomplex_t IT_0230 = IT_0043*IT_0228;
    const ccomplex_t IT_0231 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0232 = IT_0043*IT_0231;
    const ccomplex_t IT_0233 = IT_0046*IT_0231;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*(IT_0229 + IT_0230 + 
      -IT_0232 + -IT_0233);
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = (-0.5)*IT_0235;
    const ccomplex_t IT_0237 = IT_0227*IT_0236;
    const ccomplex_t IT_0238 = IT_0056*IT_0237;
    const ccomplex_t IT_0239 = (-2)*IT_0218 + (-2)*IT_0238;
    const ccomplex_t IT_0240 = (-2)*IT_0239;
    const ccomplex_t IT_0241 = IT_0186 + IT_0240;
    const ccomplex_t IT_0242 = 2*conj(IT_0185);
    const ccomplex_t IT_0243 = (-2)*conj(IT_0239);
    const ccomplex_t IT_0244 = IT_0242 + IT_0243;
    const ccomplex_t IT_0245 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0246 = IT_0081*IT_0245;
    const ccomplex_t IT_0247 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0248 = IT_0085*IT_0247;
    const ccomplex_t IT_0249 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0250 = IT_0088*IT_0249;
    const ccomplex_t IT_0251 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0252 = IT_0091*IT_0251;
    const ccomplex_t IT_0253 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0254 = IT_0081*IT_0253;
    const ccomplex_t IT_0255 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0256 = IT_0085*IT_0255;
    const ccomplex_t IT_0257 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0258 = IT_0088*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0260 = IT_0091*IT_0259;
    const ccomplex_t IT_0261 = 0.5*IT_0246 + (-0.5)*IT_0248 + (-0.5)*IT_0250 +
       0.5*IT_0252 + 0.5*IT_0254 + (-0.5)*IT_0256 + (-0.5)*IT_0258 + 0.5*IT_0260;
    const ccomplex_t IT_0262 = IT_0103*IT_0261;
    const ccomplex_t IT_0263 = IT_0106*IT_0262;
    const ccomplex_t IT_0264 = IT_0015*IT_0247;
    const ccomplex_t IT_0265 = IT_0021*IT_0249;
    const ccomplex_t IT_0266 = IT_0015*IT_0255;
    const ccomplex_t IT_0267 = IT_0021*IT_0257;
    const ccomplex_t IT_0268 = IT_0264 + IT_0265 + IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = IT_0029*IT_0245;
    const ccomplex_t IT_0270 = IT_0032*IT_0251;
    const ccomplex_t IT_0271 = IT_0029*IT_0253;
    const ccomplex_t IT_0272 = IT_0032*IT_0259;
    const ccomplex_t IT_0273 = -IT_0269 + -IT_0270 + -IT_0271 + -IT_0272;
    const ccomplex_t IT_0274 = IT_0268 + IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = 0.5*IT_0275;
    const ccomplex_t IT_0277 = IT_0133*IT_0276;
    const ccomplex_t IT_0278 = IT_0109*IT_0277;
    const ccomplex_t IT_0279 = IT_0130*IT_0276;
    const ccomplex_t IT_0280 = IT_0109*IT_0279;
    const ccomplex_t IT_0281 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0282 = IT_0015*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0284 = IT_0021*IT_0283;
    const ccomplex_t IT_0285 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0286 = IT_0015*IT_0285;
    const ccomplex_t IT_0287 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0288 = IT_0021*IT_0287;
    const ccomplex_t IT_0289 = IT_0282 + IT_0284 + IT_0286 + IT_0288;
    const ccomplex_t IT_0290 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0291 = IT_0029*IT_0290;
    const ccomplex_t IT_0292 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0293 = IT_0032*IT_0292;
    const ccomplex_t IT_0294 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0295 = IT_0029*IT_0294;
    const ccomplex_t IT_0296 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0297 = IT_0032*IT_0296;
    const ccomplex_t IT_0298 = -IT_0291 + -IT_0293 + -IT_0295 + -IT_0297;
    const ccomplex_t IT_0299 = IT_0289 + IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = 0.5*IT_0300;
    const ccomplex_t IT_0302 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0303 = IT_0043*IT_0302;
    const ccomplex_t IT_0304 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0305 = IT_0043*IT_0304;
    const ccomplex_t IT_0306 = IT_0046*IT_0302;
    const ccomplex_t IT_0307 = IT_0046*IT_0304;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*(IT_0303 + -IT_0305 +
       IT_0306 + -IT_0307);
    const ccomplex_t IT_0309 = -IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*IT_0309;
    const ccomplex_t IT_0311 = IT_0301*IT_0310;
    const ccomplex_t IT_0312 = IT_0167*IT_0311;
    const ccomplex_t IT_0313 = IT_0081*IT_0247;
    const ccomplex_t IT_0314 = IT_0085*IT_0245;
    const ccomplex_t IT_0315 = IT_0088*IT_0251;
    const ccomplex_t IT_0316 = IT_0091*IT_0249;
    const ccomplex_t IT_0317 = IT_0081*IT_0255;
    const ccomplex_t IT_0318 = IT_0085*IT_0253;
    const ccomplex_t IT_0319 = IT_0088*IT_0259;
    const ccomplex_t IT_0320 = IT_0091*IT_0257;
    const ccomplex_t IT_0321 = 0.5*IT_0313 + 0.5*IT_0314 + (-0.5)*IT_0315 + (
      -0.5)*IT_0316 + 0.5*IT_0317 + 0.5*IT_0318 + (-0.5)*IT_0319 + (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = IT_0173*IT_0321;
    const ccomplex_t IT_0323 = IT_0171*IT_0322;
    const ccomplex_t IT_0324 = IT_0263 + IT_0278 + -IT_0280 + (-2)*IT_0312 +
       IT_0323;
    const ccomplex_t IT_0325 = 2*IT_0239;
    const ccomplex_t IT_0326 = (-2)*IT_0185;
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = 2*conj(IT_0239);
    const ccomplex_t IT_0329 = (-2)*conj(IT_0185);
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0332 = -IT_0331;
    const ccomplex_t IT_0333 = s_12 + IT_0332;
    const ccomplex_t IT_0334 = m_N_2*IT_0333;
    const ccomplex_t IT_0335 = m_N_2*IT_0131;
    const ccomplex_t IT_0336 = IT_0109*IT_0335;
    const ccomplex_t IT_0337 = m_N_2*IT_0134;
    const ccomplex_t IT_0338 = -IT_0337;
    const ccomplex_t IT_0339 = IT_0335 + IT_0338;
    const ccomplex_t IT_0340 = IT_0109*IT_0339;
    const ccomplex_t IT_0341 = m_N_2*IT_0165;
    const ccomplex_t IT_0342 = -IT_0341;
    const ccomplex_t IT_0343 = IT_0164*IT_0310;
    const ccomplex_t IT_0344 = m_N_1*IT_0343;
    const ccomplex_t IT_0345 = IT_0342 + IT_0344;
    const ccomplex_t IT_0346 = IT_0167*IT_0345;
    const ccomplex_t IT_0347 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*IT_0347;
    const ccomplex_t IT_0349 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0350 = IT_0046*IT_0349;
    const ccomplex_t IT_0351 = IT_0043*IT_0349;
    const ccomplex_t IT_0352 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0353 = IT_0046*IT_0352;
    const ccomplex_t IT_0354 = IT_0043*IT_0352;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*(IT_0350 + IT_0351 + 
      -IT_0353 + -IT_0354);
    const ccomplex_t IT_0356 = (-0.5)*IT_0355;
    const ccomplex_t IT_0357 = IT_0276*IT_0356;
    const ccomplex_t IT_0358 = m_N_4*IT_0357;
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = 0.5*IT_0355;
    const ccomplex_t IT_0361 = IT_0276*IT_0360;
    const ccomplex_t IT_0362 = m_N_4*IT_0361;
    const ccomplex_t IT_0363 = IT_0359 + IT_0362;
    const ccomplex_t IT_0364 = IT_0348*IT_0363;
    const ccomplex_t IT_0365 = 0.5*IT_0336 + 0.5*IT_0340 + IT_0346 + (-0.5)
      *IT_0364;
    const ccomplex_t IT_0366 = cpow(IT_0013, -4);
    const ccomplex_t IT_0367 = IT_0012*IT_0019*IT_0080*IT_0366;
    const ccomplex_t IT_0368 = IT_0018*IT_0019*IT_0084*IT_0366;
    const ccomplex_t IT_0369 = cpow(IT_0019, -3);
    const ccomplex_t IT_0370 = IT_0012*IT_0080*IT_0369;
    const ccomplex_t IT_0371 = IT_0018*IT_0084*IT_0369;
    const ccomplex_t IT_0372 = cpow(IT_0013, -2);
    const ccomplex_t IT_0373 = IT_0012*IT_0020*IT_0080*IT_0372;
    const ccomplex_t IT_0374 = IT_0018*IT_0020*IT_0084*IT_0372;
    const ccomplex_t IT_0375 = cpow(IT_0019, -2);
    const ccomplex_t IT_0376 = cpow(IT_0372 + IT_0375, -1);
    const ccomplex_t IT_0377 = (IT_0367 + IT_0368 + IT_0370 + IT_0371 + 2
      *IT_0373 + 2*IT_0374)*IT_0376;
    const ccomplex_t IT_0378 = m_W*e_em;
    const ccomplex_t IT_0379 = IT_0377*IT_0378;
    const ccomplex_t IT_0380 = (0 + _Complex_I*1)*IT_0379;
    const ccomplex_t IT_0381 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0004 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0382 = IT_0380*IT_0381;
    const ccomplex_t IT_0383 = (0 + _Complex_I*1)*IT_0382;
    const ccomplex_t IT_0384 = IT_0076*IT_0383;
    const ccomplex_t IT_0385 = IT_0109*IT_0337;
    const ccomplex_t IT_0386 = IT_0052*IT_0207;
    const ccomplex_t IT_0387 = m_N_3*IT_0386;
    const ccomplex_t IT_0388 = m_N_2*IT_0217;
    const ccomplex_t IT_0389 = -IT_0388;
    const ccomplex_t IT_0390 = IT_0387 + IT_0389;
    const ccomplex_t IT_0391 = IT_0011*IT_0390;
    const ccomplex_t IT_0392 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0393 = (0 + _Complex_I*1)*IT_0392;
    const ccomplex_t IT_0394 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0395 = IT_0015*IT_0394;
    const ccomplex_t IT_0396 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0397 = IT_0021*IT_0396;
    const ccomplex_t IT_0398 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0399 = IT_0015*IT_0398;
    const ccomplex_t IT_0400 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0401 = IT_0021*IT_0400;
    const ccomplex_t IT_0402 = IT_0395 + IT_0397 + IT_0399 + IT_0401;
    const ccomplex_t IT_0403 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0404 = IT_0029*IT_0403;
    const ccomplex_t IT_0405 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0406 = IT_0032*IT_0405;
    const ccomplex_t IT_0407 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0408 = IT_0029*IT_0407;
    const ccomplex_t IT_0409 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0410 = IT_0032*IT_0409;
    const ccomplex_t IT_0411 = -IT_0404 + -IT_0406 + -IT_0408 + -IT_0410;
    const ccomplex_t IT_0412 = IT_0402 + IT_0411;
    const ccomplex_t IT_0413 = (0 + _Complex_I*1)*IT_0412;
    const ccomplex_t IT_0414 = 0.5*IT_0413;
    const ccomplex_t IT_0415 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0416 = IT_0046*IT_0415;
    const ccomplex_t IT_0417 = IT_0043*IT_0415;
    const ccomplex_t IT_0418 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0419 = IT_0046*IT_0418;
    const ccomplex_t IT_0420 = IT_0043*IT_0418;
    const ccomplex_t IT_0421 = (0 + _Complex_I*1)*(IT_0416 + IT_0417 + 
      -IT_0419 + -IT_0420);
    const ccomplex_t IT_0422 = (-0.5)*IT_0421;
    const ccomplex_t IT_0423 = IT_0414*IT_0422;
    const ccomplex_t IT_0424 = m_N_1*IT_0423;
    const ccomplex_t IT_0425 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0426 = IT_0046*IT_0425;
    const ccomplex_t IT_0427 = IT_0043*IT_0425;
    const ccomplex_t IT_0428 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0429 = IT_0046*IT_0428;
    const ccomplex_t IT_0430 = IT_0043*IT_0428;
    const ccomplex_t IT_0431 = (0 + _Complex_I*1)*(IT_0426 + IT_0427 + 
      -IT_0429 + -IT_0430);
    const ccomplex_t IT_0432 = (-0.5)*IT_0431;
    const ccomplex_t IT_0433 = IT_0414*IT_0432;
    const ccomplex_t IT_0434 = m_N_4*IT_0433;
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = IT_0424 + IT_0435;
    const ccomplex_t IT_0437 = IT_0393*IT_0436;
    const ccomplex_t IT_0438 = cpow((-2)*s_13 + IT_0008 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0439 = (0 + _Complex_I*1)*IT_0438;
    const ccomplex_t IT_0440 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0441 = IT_0015*IT_0440;
    const ccomplex_t IT_0442 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0443 = IT_0029*IT_0442;
    const ccomplex_t IT_0444 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0445 = IT_0032*IT_0444;
    const ccomplex_t IT_0446 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0447 = IT_0021*IT_0446;
    const ccomplex_t IT_0448 = (0 + _Complex_I*1)*(IT_0441 + -IT_0443 + 
      -IT_0445 + IT_0447);
    const ccomplex_t IT_0449 = IT_0076*IT_0448;
    const ccomplex_t IT_0450 = m_N_2*IT_0449;
    const ccomplex_t IT_0451 = IT_0236*IT_0448;
    const ccomplex_t IT_0452 = m_N_4*IT_0451;
    const ccomplex_t IT_0453 = -IT_0452;
    const ccomplex_t IT_0454 = IT_0450 + IT_0453;
    const ccomplex_t IT_0455 = IT_0439*IT_0454;
    const ccomplex_t IT_0456 = m_N_2*IT_0237;
    const ccomplex_t IT_0457 = -IT_0456;
    const ccomplex_t IT_0458 = IT_0076*IT_0227;
    const ccomplex_t IT_0459 = m_N_4*IT_0458;
    const ccomplex_t IT_0460 = IT_0457 + IT_0459;
    const ccomplex_t IT_0461 = IT_0056*IT_0460;
    const ccomplex_t IT_0462 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0463 = (0 + _Complex_I*1)*IT_0462;
    const ccomplex_t IT_0464 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0465 = IT_0015*IT_0464;
    const ccomplex_t IT_0466 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0467 = IT_0021*IT_0466;
    const ccomplex_t IT_0468 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0469 = IT_0015*IT_0468;
    const ccomplex_t IT_0470 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0471 = IT_0021*IT_0470;
    const ccomplex_t IT_0472 = IT_0465 + IT_0467 + IT_0469 + IT_0471;
    const ccomplex_t IT_0473 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0474 = IT_0029*IT_0473;
    const ccomplex_t IT_0475 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0476 = IT_0032*IT_0475;
    const ccomplex_t IT_0477 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0478 = IT_0029*IT_0477;
    const ccomplex_t IT_0479 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0480 = IT_0032*IT_0479;
    const ccomplex_t IT_0481 = -IT_0474 + -IT_0476 + -IT_0478 + -IT_0480;
    const ccomplex_t IT_0482 = IT_0472 + IT_0481;
    const ccomplex_t IT_0483 = (0 + _Complex_I*1)*IT_0482;
    const ccomplex_t IT_0484 = 0.5*IT_0483;
    const ccomplex_t IT_0485 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0486 = IT_0046*IT_0485;
    const ccomplex_t IT_0487 = IT_0043*IT_0485;
    const ccomplex_t IT_0488 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0489 = IT_0046*IT_0488;
    const ccomplex_t IT_0490 = IT_0043*IT_0488;
    const ccomplex_t IT_0491 = (0 + _Complex_I*1)*(IT_0486 + IT_0487 + 
      -IT_0489 + -IT_0490);
    const ccomplex_t IT_0492 = (-0.5)*IT_0491;
    const ccomplex_t IT_0493 = IT_0484*IT_0492;
    const ccomplex_t IT_0494 = m_N_3*IT_0493;
    const ccomplex_t IT_0495 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0496 = IT_0046*IT_0495;
    const ccomplex_t IT_0497 = IT_0043*IT_0495;
    const ccomplex_t IT_0498 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0499 = IT_0046*IT_0498;
    const ccomplex_t IT_0500 = IT_0043*IT_0498;
    const ccomplex_t IT_0501 = (0 + _Complex_I*1)*(IT_0496 + IT_0497 + 
      -IT_0499 + -IT_0500);
    const ccomplex_t IT_0502 = (-0.5)*IT_0501;
    const ccomplex_t IT_0503 = IT_0484*IT_0502;
    const ccomplex_t IT_0504 = m_N_4*IT_0503;
    const ccomplex_t IT_0505 = -IT_0504;
    const ccomplex_t IT_0506 = IT_0494 + IT_0505;
    const ccomplex_t IT_0507 = IT_0463*IT_0506;
    const ccomplex_t IT_0508 = (-0.5)*IT_0364 + -IT_0384 + (-0.5)*IT_0385 + 
      -IT_0391 + -IT_0437 + -IT_0455 + -IT_0461 + -IT_0507;
    const ccomplex_t IT_0509 = conj(IT_0365) + conj(IT_0508);
    const ccomplex_t IT_0510 = IT_0365 + IT_0508;
    const ccomplex_t IT_0511 = m_N_2*IT_0277;
    const ccomplex_t IT_0512 = IT_0109*IT_0511;
    const ccomplex_t IT_0513 = m_N_2*IT_0279;
    const ccomplex_t IT_0514 = -IT_0513;
    const ccomplex_t IT_0515 = IT_0511 + IT_0514;
    const ccomplex_t IT_0516 = IT_0109*IT_0515;
    const ccomplex_t IT_0517 = m_N_2*IT_0311;
    const ccomplex_t IT_0518 = -IT_0517;
    const ccomplex_t IT_0519 = IT_0143*IT_0301;
    const ccomplex_t IT_0520 = m_N_1*IT_0519;
    const ccomplex_t IT_0521 = IT_0518 + IT_0520;
    const ccomplex_t IT_0522 = IT_0167*IT_0521;
    const ccomplex_t IT_0523 = IT_0122*IT_0356;
    const ccomplex_t IT_0524 = m_N_4*IT_0523;
    const ccomplex_t IT_0525 = IT_0122*IT_0360;
    const ccomplex_t IT_0526 = m_N_4*IT_0525;
    const ccomplex_t IT_0527 = -IT_0526;
    const ccomplex_t IT_0528 = IT_0524 + IT_0527;
    const ccomplex_t IT_0529 = IT_0348*IT_0528;
    const ccomplex_t IT_0530 = 0.5*IT_0512 + 0.5*IT_0516 + IT_0522 + (-0.5)
      *IT_0529;
    const ccomplex_t IT_0531 = IT_0236*IT_0383;
    const ccomplex_t IT_0532 = IT_0109*IT_0513;
    const ccomplex_t IT_0533 = IT_0042*IT_0216;
    const ccomplex_t IT_0534 = m_N_3*IT_0533;
    const ccomplex_t IT_0535 = m_N_2*IT_0053;
    const ccomplex_t IT_0536 = -IT_0535;
    const ccomplex_t IT_0537 = IT_0534 + IT_0536;
    const ccomplex_t IT_0538 = IT_0011*IT_0537;
    const ccomplex_t IT_0539 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0540 = IT_0015*IT_0539;
    const ccomplex_t IT_0541 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0542 = IT_0021*IT_0541;
    const ccomplex_t IT_0543 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0544 = IT_0015*IT_0543;
    const ccomplex_t IT_0545 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0546 = IT_0021*IT_0545;
    const ccomplex_t IT_0547 = IT_0540 + IT_0542 + IT_0544 + IT_0546;
    const ccomplex_t IT_0548 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0549 = IT_0029*IT_0548;
    const ccomplex_t IT_0550 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0551 = IT_0032*IT_0550;
    const ccomplex_t IT_0552 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0553 = IT_0029*IT_0552;
    const ccomplex_t IT_0554 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0555 = IT_0032*IT_0554;
    const ccomplex_t IT_0556 = -IT_0549 + -IT_0551 + -IT_0553 + -IT_0555;
    const ccomplex_t IT_0557 = IT_0547 + IT_0556;
    const ccomplex_t IT_0558 = (0 + _Complex_I*1)*IT_0557;
    const ccomplex_t IT_0559 = 0.5*IT_0558;
    const ccomplex_t IT_0560 = IT_0432*IT_0559;
    const ccomplex_t IT_0561 = m_N_1*IT_0560;
    const ccomplex_t IT_0562 = IT_0422*IT_0559;
    const ccomplex_t IT_0563 = m_N_4*IT_0562;
    const ccomplex_t IT_0564 = -IT_0563;
    const ccomplex_t IT_0565 = IT_0561 + IT_0564;
    const ccomplex_t IT_0566 = IT_0393*IT_0565;
    const ccomplex_t IT_0567 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0568 = IT_0015*IT_0567;
    const ccomplex_t IT_0569 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0570 = IT_0029*IT_0569;
    const ccomplex_t IT_0571 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0572 = IT_0032*IT_0571;
    const ccomplex_t IT_0573 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0574 = IT_0021*IT_0573;
    const ccomplex_t IT_0575 = (0 + _Complex_I*1)*(IT_0568 + -IT_0570 + 
      -IT_0572 + IT_0574);
    const ccomplex_t IT_0576 = IT_0236*IT_0575;
    const ccomplex_t IT_0577 = m_N_2*IT_0576;
    const ccomplex_t IT_0578 = IT_0076*IT_0575;
    const ccomplex_t IT_0579 = m_N_4*IT_0578;
    const ccomplex_t IT_0580 = -IT_0579;
    const ccomplex_t IT_0581 = IT_0577 + IT_0580;
    const ccomplex_t IT_0582 = IT_0439*IT_0581;
    const ccomplex_t IT_0583 = IT_0065*IT_0236;
    const ccomplex_t IT_0584 = m_N_4*IT_0583;
    const ccomplex_t IT_0585 = m_N_2*IT_0077;
    const ccomplex_t IT_0586 = -IT_0585;
    const ccomplex_t IT_0587 = IT_0584 + IT_0586;
    const ccomplex_t IT_0588 = IT_0056*IT_0587;
    const ccomplex_t IT_0589 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0590 = IT_0015*IT_0589;
    const ccomplex_t IT_0591 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0592 = IT_0021*IT_0591;
    const ccomplex_t IT_0593 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0594 = IT_0015*IT_0593;
    const ccomplex_t IT_0595 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0596 = IT_0021*IT_0595;
    const ccomplex_t IT_0597 = IT_0590 + IT_0592 + IT_0594 + IT_0596;
    const ccomplex_t IT_0598 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0599 = IT_0029*IT_0598;
    const ccomplex_t IT_0600 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0601 = IT_0032*IT_0600;
    const ccomplex_t IT_0602 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0603 = IT_0029*IT_0602;
    const ccomplex_t IT_0604 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0605 = IT_0032*IT_0604;
    const ccomplex_t IT_0606 = -IT_0599 + -IT_0601 + -IT_0603 + -IT_0605;
    const ccomplex_t IT_0607 = IT_0597 + IT_0606;
    const ccomplex_t IT_0608 = (0 + _Complex_I*1)*IT_0607;
    const ccomplex_t IT_0609 = 0.5*IT_0608;
    const ccomplex_t IT_0610 = IT_0502*IT_0609;
    const ccomplex_t IT_0611 = m_N_3*IT_0610;
    const ccomplex_t IT_0612 = IT_0492*IT_0609;
    const ccomplex_t IT_0613 = m_N_4*IT_0612;
    const ccomplex_t IT_0614 = -IT_0613;
    const ccomplex_t IT_0615 = IT_0611 + IT_0614;
    const ccomplex_t IT_0616 = IT_0463*IT_0615;
    const ccomplex_t IT_0617 = (-0.5)*IT_0529 + -IT_0531 + (-0.5)*IT_0532 + 
      -IT_0538 + -IT_0566 + -IT_0582 + -IT_0588 + -IT_0616;
    const ccomplex_t IT_0618 = conj(IT_0530) + conj(IT_0617);
    const ccomplex_t IT_0619 = 2*IT_0324;
    const ccomplex_t IT_0620 = (-2)*IT_0079;
    const ccomplex_t IT_0621 = IT_0619 + IT_0620;
    const ccomplex_t IT_0622 = IT_0530 + IT_0617;
    const ccomplex_t IT_0623 = 2*conj(IT_0324);
    const ccomplex_t IT_0624 = (-2)*conj(IT_0079);
    const ccomplex_t IT_0625 = IT_0623 + IT_0624;
    const ccomplex_t IT_0626 = s_12*IT_0333;
    const ccomplex_t IT_0627 = IT_0348*IT_0525;
    const ccomplex_t IT_0628 = IT_0393*IT_0562;
    const ccomplex_t IT_0629 = IT_0463*IT_0612;
    const ccomplex_t IT_0630 = IT_0439*IT_0578;
    const ccomplex_t IT_0631 = IT_0627 + 2*IT_0628 + 2*IT_0629 + 2*IT_0630;
    const ccomplex_t IT_0632 = IT_0348*IT_0523;
    const ccomplex_t IT_0633 = IT_0107 + IT_0184 + IT_0632;
    const ccomplex_t IT_0634 = IT_0348*IT_0357;
    const ccomplex_t IT_0635 = IT_0439*IT_0451;
    const ccomplex_t IT_0636 = IT_0393*IT_0433;
    const ccomplex_t IT_0637 = IT_0463*IT_0503;
    const ccomplex_t IT_0638 = IT_0634 + 2*IT_0635 + 2*IT_0636 + 2*IT_0637;
    const ccomplex_t IT_0639 = IT_0348*IT_0361;
    const ccomplex_t IT_0640 = IT_0263 + IT_0323 + IT_0639;
    const ccomplex_t IT_0641 = 2*IT_0079;
    const ccomplex_t IT_0642 = (-2)*IT_0324;
    const ccomplex_t IT_0643 = IT_0641 + IT_0642;
    const ccomplex_t IT_0644 = 2*conj(IT_0079);
    const ccomplex_t IT_0645 = (-2)*conj(IT_0324);
    const ccomplex_t IT_0646 = IT_0644 + IT_0645;
    const ccomplex_t IT_0647 = s_12*s_14;
    const ccomplex_t IT_0648 = s_24*IT_0004;
    const ccomplex_t IT_0649 = -IT_0648;
    const ccomplex_t IT_0650 = IT_0647 + IT_0649;
    const ccomplex_t IT_0651 = (-0.5)*IT_0132 + 0.5*IT_0135 + IT_0168;
    const ccomplex_t IT_0652 = (-0.5)*IT_0627 + -IT_0628 + -IT_0629 + -IT_0630;
    const ccomplex_t IT_0653 = conj(IT_0651) + conj(IT_0652);
    const ccomplex_t IT_0654 = IT_0651 + IT_0652;
    const ccomplex_t IT_0655 = (-0.5)*IT_0632;
    const ccomplex_t IT_0656 = (-0.5)*IT_0278 + 0.5*IT_0280 + IT_0312;
    const ccomplex_t IT_0657 = -IT_0637;
    const ccomplex_t IT_0658 = (-0.5)*IT_0634 + -IT_0635 + -IT_0636;
    const ccomplex_t IT_0659 = conj(IT_0656) + conj(IT_0657) + conj(IT_0658);
    const ccomplex_t IT_0660 = IT_0656 + IT_0657 + IT_0658;
    const ccomplex_t IT_0661 = (-0.5)*IT_0639;
    const ccomplex_t IT_0662 = conj(IT_0054) + conj(IT_0078) + conj(IT_0661);
    const ccomplex_t IT_0663 = IT_0054 + IT_0078 + IT_0661;
    const ccomplex_t IT_0664 = IT_0650*(IT_0241*IT_0653 + IT_0244*IT_0654 +
       IT_0330*(IT_0218 + IT_0238 + IT_0655) + IT_0327*(conj(IT_0218) + conj
      (IT_0238) + conj(IT_0655)) + IT_0621*IT_0659 + IT_0625*IT_0660 + IT_0643
      *IT_0662 + IT_0646*IT_0663);
    const ccomplex_t IT_0665 = IT_0000*IT_0333;
    const ccomplex_t IT_0666 = 2*IT_0633;
    const ccomplex_t IT_0667 = (-2)*IT_0631;
    const ccomplex_t IT_0668 = IT_0666 + IT_0667;
    const ccomplex_t IT_0669 = 2*conj(IT_0633);
    const ccomplex_t IT_0670 = (-2)*conj(IT_0631);
    const ccomplex_t IT_0671 = IT_0669 + IT_0670;
    const ccomplex_t IT_0672 = 2*IT_0631;
    const ccomplex_t IT_0673 = (-2)*IT_0633;
    const ccomplex_t IT_0674 = IT_0672 + IT_0673;
    const ccomplex_t IT_0675 = 2*conj(IT_0631);
    const ccomplex_t IT_0676 = (-2)*conj(IT_0633);
    const ccomplex_t IT_0677 = IT_0675 + IT_0676;
    const ccomplex_t IT_0678 = pow(s_24, 2);
    const ccomplex_t IT_0679 = IT_0001*IT_0678;
    const ccomplex_t IT_0680 = -IT_0009;
    const ccomplex_t IT_0681 = IT_0679 + IT_0680;
    const ccomplex_t IT_0682 = IT_0000*IT_0681;
    const ccomplex_t IT_0683 = m_N_2*IT_0681;
    const ccomplex_t IT_0684 = 2*IT_0638;
    const ccomplex_t IT_0685 = (-2)*IT_0640;
    const ccomplex_t IT_0686 = IT_0684 + IT_0685;
    const ccomplex_t IT_0687 = 2*conj(IT_0638);
    const ccomplex_t IT_0688 = (-2)*conj(IT_0640);
    const ccomplex_t IT_0689 = IT_0687 + IT_0688;
    const ccomplex_t IT_0690 = s_24*IT_0169;
    const ccomplex_t IT_0691 = IT_0001*IT_0690;
    const ccomplex_t IT_0692 = 0.5*IT_0691;
    const ccomplex_t IT_0693 = s_24 + IT_0692;
    const ccomplex_t IT_0694 = m_N_2*IT_0693;
    const ccomplex_t IT_0695 = 2*IT_0694;
    const ccomplex_t IT_0696 = 2*IT_0365;
    const ccomplex_t IT_0697 = 2*IT_0508;
    const ccomplex_t IT_0698 = IT_0696 + IT_0697;
    const ccomplex_t IT_0699 = 2*conj(IT_0365);
    const ccomplex_t IT_0700 = 2*conj(IT_0508);
    const ccomplex_t IT_0701 = IT_0699 + IT_0700;
    const ccomplex_t IT_0702 = conj(IT_0657) + conj(IT_0658);
    const ccomplex_t IT_0703 = 2*IT_0530;
    const ccomplex_t IT_0704 = 2*IT_0617;
    const ccomplex_t IT_0705 = IT_0703 + IT_0704;
    const ccomplex_t IT_0706 = IT_0657 + IT_0658;
    const ccomplex_t IT_0707 = 2*conj(IT_0530);
    const ccomplex_t IT_0708 = 2*conj(IT_0617);
    const ccomplex_t IT_0709 = IT_0707 + IT_0708;
    const ccomplex_t IT_0710 = (-2)*IT_0365;
    const ccomplex_t IT_0711 = (-2)*IT_0508;
    const ccomplex_t IT_0712 = IT_0710 + IT_0711;
    const ccomplex_t IT_0713 = (-2)*conj(IT_0365);
    const ccomplex_t IT_0714 = (-2)*conj(IT_0508);
    const ccomplex_t IT_0715 = IT_0713 + IT_0714;
    const ccomplex_t IT_0716 = (-2)*IT_0530;
    const ccomplex_t IT_0717 = (-2)*IT_0617;
    const ccomplex_t IT_0718 = IT_0716 + IT_0717;
    const ccomplex_t IT_0719 = (-2)*conj(IT_0530);
    const ccomplex_t IT_0720 = (-2)*conj(IT_0617);
    const ccomplex_t IT_0721 = IT_0719 + IT_0720;
    const ccomplex_t IT_0722 = (-0.25)*IT_0691;
    const ccomplex_t IT_0723 = s_24 + IT_0722;
    const ccomplex_t IT_0724 = m_N_2*IT_0723;
    const ccomplex_t IT_0725 = (-4)*IT_0724;
    const ccomplex_t IT_0726 = conj(IT_0218) + conj(IT_0238);
    const ccomplex_t IT_0727 = IT_0218 + IT_0238;
    const ccomplex_t IT_0728 = conj(IT_0054) + conj(IT_0078);
    const ccomplex_t IT_0729 = IT_0054 + IT_0078;
    const ccomplex_t IT_0730 = s_14*IT_0169;
    const ccomplex_t IT_0731 = IT_0001*IT_0730;
    const ccomplex_t IT_0732 = (-0.25)*IT_0731;
    const ccomplex_t IT_0733 = s_14 + IT_0732;
    const ccomplex_t IT_0734 = m_N_4*IT_0733;
    const ccomplex_t IT_0735 = (-4)*IT_0734;
    const ccomplex_t IT_0736 = 0.5*IT_0731;
    const ccomplex_t IT_0737 = s_14 + IT_0736;
    const ccomplex_t IT_0738 = m_N_4*IT_0737;
    const ccomplex_t IT_0739 = 2*IT_0738;
    const ccomplex_t IT_0740 = s_12*IT_0169;
    const ccomplex_t IT_0741 = s_14*s_24;
    const ccomplex_t IT_0742 = (-2)*IT_0741;
    const ccomplex_t IT_0743 = IT_0740 + IT_0742;
    const ccomplex_t IT_0744 = IT_0001*IT_0743;
    const ccomplex_t IT_0745 = (-0.5)*IT_0744;
    const ccomplex_t IT_0746 = s_12 + IT_0745;
    const ccomplex_t IT_0747 = 4*IT_0530;
    const ccomplex_t IT_0748 = 6*conj(IT_0365);
    const ccomplex_t IT_0749 = s_12*IT_0006;
    const ccomplex_t IT_0750 = s_12*IT_0681;
    const ccomplex_t IT_0751 = 2*IT_0640;
    const ccomplex_t IT_0752 = s_12*s_24;
    const ccomplex_t IT_0753 = s_14*IT_0009;
    const ccomplex_t IT_0754 = -IT_0753;
    const ccomplex_t IT_0755 = IT_0752 + IT_0754;
    const ccomplex_t IT_0756 = (-2)*IT_0638;
    const ccomplex_t IT_0757 = 2*conj(IT_0640);
    const ccomplex_t IT_0758 = (-2)*conj(IT_0638);
    const ccomplex_t IT_0759 = (conj(IT_0218) + conj(IT_0238) + conj(IT_0655))
      *IT_0668 + (IT_0218 + IT_0238 + IT_0655)*IT_0671 + IT_0653*IT_0674 +
       IT_0654*IT_0677 + IT_0659*IT_0686 + IT_0660*IT_0689 + IT_0662*(IT_0751 +
       IT_0756) + IT_0663*(IT_0757 + IT_0758);
    const ccomplex_t IT_0760 = m_N_4*IT_0006;
    const ccomplex_t IT_0761 = m_N_4*IT_0333;
    const ccomplex_t IT_0762 = IT_0000*IT_0169;
    const ccomplex_t IT_0763 = 3*IT_0762;
    const ccomplex_t IT_0764 = 2*IT_0655;
    const ccomplex_t IT_0765 = (-2)*IT_0652;
    const ccomplex_t IT_0766 = IT_0764 + IT_0765;
    const ccomplex_t IT_0767 = 2*conj(IT_0655);
    const ccomplex_t IT_0768 = (-2)*conj(IT_0652);
    const ccomplex_t IT_0769 = IT_0767 + IT_0768;
    const ccomplex_t IT_0770 = 2*IT_0651;
    const ccomplex_t IT_0771 = (-2)*IT_0238;
    const ccomplex_t IT_0772 = (-2)*IT_0218;
    const ccomplex_t IT_0773 = IT_0770 + IT_0771 + IT_0772;
    const ccomplex_t IT_0774 = 2*conj(IT_0651);
    const ccomplex_t IT_0775 = (-2)*conj(IT_0238);
    const ccomplex_t IT_0776 = (-2)*conj(IT_0218);
    const ccomplex_t IT_0777 = IT_0774 + IT_0775 + IT_0776;
    const ccomplex_t IT_0778 = 2*IT_0661;
    const ccomplex_t IT_0779 = 2*conj(IT_0661);
    const ccomplex_t IT_0780 = 2*IT_0656;
    const ccomplex_t IT_0781 = 2*conj(IT_0656);
    const ccomplex_t IT_0782 = (-2)*IT_0655;
    const ccomplex_t IT_0783 = (-2)*conj(IT_0655);
    const ccomplex_t IT_0784 = (-2)*IT_0651;
    const ccomplex_t IT_0785 = (-2)*conj(IT_0651);
    const ccomplex_t IT_0786 = IT_0728*IT_0766 + IT_0729*IT_0769 + IT_0702
      *IT_0773 + IT_0706*IT_0777 + IT_0726*IT_0778 + IT_0727*IT_0779 + conj
      (IT_0652)*IT_0780 + IT_0652*IT_0781 + conj(IT_0656)*IT_0782 + IT_0656
      *IT_0783 + conj(IT_0661)*IT_0784 + IT_0661*IT_0785;
    const ccomplex_t IT_0787 = (-3)*IT_0762;
    const ccomplex_t IT_0788 = IT_0702*IT_0766 + IT_0706*IT_0769 + IT_0728
      *IT_0773 + IT_0729*IT_0777 + conj(IT_0652)*IT_0778 + IT_0652*IT_0779 +
       IT_0726*IT_0780 + IT_0727*IT_0781 + conj(IT_0661)*IT_0782 + IT_0661
      *IT_0783 + conj(IT_0656)*IT_0784 + IT_0656*IT_0785;
    const ccomplex_t IT_0789 = pow(m_Z, 4);
    const ccomplex_t IT_0790 = s_12*IT_0789;
    const ccomplex_t IT_0791 = IT_0001*IT_0790;
    const ccomplex_t IT_0792 = (-0.25)*IT_0791;
    const ccomplex_t IT_0793 = IT_0741 + IT_0792;
    const ccomplex_t IT_0794 = -IT_0652;
    const ccomplex_t IT_0795 = IT_0655 + IT_0794;
    const ccomplex_t IT_0796 = 8*conj(IT_0655);
    const ccomplex_t IT_0797 = 8*IT_0656;
    const ccomplex_t IT_0798 = 8*conj(IT_0656);
    const ccomplex_t IT_0799 = 8*IT_0661;
    const ccomplex_t IT_0800 = (-8)*IT_0657;
    const ccomplex_t IT_0801 = 8*conj(IT_0661);
    const ccomplex_t IT_0802 = (-8)*conj(IT_0657);
    const ccomplex_t IT_0803 = -conj(IT_0652);
    const ccomplex_t IT_0804 = IT_0651*conj(IT_0655) + conj(IT_0656)*IT_0661 +
       IT_0656*conj(IT_0661) + -(conj(IT_0218) + conj(IT_0238) + -conj(IT_0651))
      *IT_0795 + (-0.125)*IT_0727*IT_0796 + (-0.125)*IT_0702*IT_0797 + (-0.125)
      *IT_0706*IT_0798 + IT_0728*(IT_0658 + (-0.125)*IT_0799 + (-0.125)*IT_0800)
       + IT_0729*(conj(IT_0658) + (-0.125)*IT_0801 + (-0.125)*IT_0802) + -
      (IT_0218 + IT_0238 + -IT_0651)*IT_0803;
    const ccomplex_t IT_0805 = (-0.5)*IT_0791;
    const ccomplex_t IT_0806 = IT_0740 + IT_0742 + IT_0805;
    const ccomplex_t IT_0807 = 4*conj(IT_0651);
    const ccomplex_t IT_0808 = 4*conj(IT_0661);
    const ccomplex_t IT_0809 = 4*conj(IT_0656);
    const ccomplex_t IT_0810 = 4*IT_0651;
    const ccomplex_t IT_0811 = (-4)*IT_0238;
    const ccomplex_t IT_0812 = 4*IT_0661;
    const ccomplex_t IT_0813 = (-4)*IT_0657;
    const ccomplex_t IT_0814 = 4*IT_0656;
    const ccomplex_t IT_0815 = (-4)*IT_0078;
    const ccomplex_t IT_0816 = IT_0651*conj(IT_0651) + IT_0656*conj(IT_0656) +
       IT_0661*conj(IT_0661) + IT_0795*(conj(IT_0655) + IT_0803) + (-0.25)
      *IT_0727*IT_0807 + (-0.25)*IT_0706*IT_0808 + (-0.25)*IT_0729*IT_0809 +
       IT_0726*(IT_0218 + (-0.25)*IT_0810 + (-0.25)*IT_0811) + IT_0702*(IT_0658 
      + (-0.25)*IT_0812 + (-0.25)*IT_0813) + IT_0728*(IT_0054 + (-0.25)*IT_0814 
      + (-0.25)*IT_0815);
    const ccomplex_t IT_0817 = IT_0007*(conj(IT_0079)*IT_0241 + IT_0079
      *IT_0244 + conj(IT_0324)*IT_0327 + IT_0324*IT_0330) + IT_0334*(IT_0241
      *IT_0509 + IT_0244*IT_0510 + IT_0618*IT_0621 + IT_0622*IT_0625) + IT_0626*
      (IT_0244*IT_0631 + IT_0241*conj(IT_0631) + IT_0330*IT_0633 + IT_0327*conj
      (IT_0633) + IT_0625*IT_0638 + IT_0621*conj(IT_0638) + conj(IT_0640)
      *IT_0643 + IT_0640*IT_0646) + IT_0664 + IT_0665*(IT_0330*IT_0638 + IT_0327
      *conj(IT_0638) + IT_0244*IT_0640 + IT_0241*conj(IT_0640) + conj(IT_0324)
      *IT_0668 + IT_0324*IT_0671 + conj(IT_0079)*IT_0674 + IT_0079*IT_0677) + 
      (conj(IT_0638)*IT_0668 + IT_0638*IT_0671 + conj(IT_0640)*IT_0674 + IT_0640
      *IT_0677)*IT_0682 + IT_0683*(IT_0509*IT_0674 + IT_0510*IT_0677 + IT_0618
      *IT_0686 + IT_0622*IT_0689) + IT_0695*(conj(IT_0652)*IT_0698 + IT_0652
      *IT_0701 + IT_0702*IT_0705 + IT_0706*IT_0709 + conj(IT_0655)*IT_0712 +
       IT_0655*IT_0715 + conj(IT_0661)*IT_0718 + IT_0661*IT_0721) + IT_0725*
      (conj(IT_0651)*IT_0712 + IT_0651*IT_0715 + conj(IT_0656)*IT_0718 + IT_0656
      *IT_0721 + IT_0698*IT_0726 + IT_0701*IT_0727 + IT_0705*IT_0728 + IT_0709
      *IT_0729) + (conj(IT_0661)*IT_0698 + IT_0661*IT_0701 + conj(IT_0655)
      *IT_0705 + IT_0655*IT_0709 + IT_0702*IT_0712 + IT_0706*IT_0715 + conj
      (IT_0652)*IT_0718 + IT_0652*IT_0721)*IT_0735 + (conj(IT_0656)*IT_0698 +
       IT_0656*IT_0701 + conj(IT_0651)*IT_0705 + IT_0651*IT_0709 + IT_0718
      *IT_0726 + IT_0721*IT_0727 + IT_0712*IT_0728 + IT_0715*IT_0729)*IT_0739 +
       4*IT_0746*(IT_0509*IT_0510 + IT_0618*(IT_0617 + 0.25*IT_0747)) + 6
      *IT_0000*(IT_0510*IT_0618 + IT_0622*(conj(IT_0508) + 0.166666666666667
      *IT_0748)) + (conj(IT_0185)*IT_0186 + conj(IT_0239)*IT_0327 + IT_0239
      *IT_0329 + conj(IT_0324)*IT_0619 + conj(IT_0079)*IT_0643 + IT_0079*IT_0645
      )*IT_0749 + IT_0750*(conj(IT_0633)*IT_0666 + conj(IT_0631)*IT_0674 +
       IT_0631*IT_0676 + conj(IT_0638)*IT_0686 + IT_0638*IT_0688 + conj(IT_0640)
      *IT_0751) + IT_0755*IT_0759 + (IT_0241*IT_0618 + conj(IT_0508)*IT_0621 +
       IT_0244*IT_0622 + IT_0508*IT_0625 + conj(IT_0324)*IT_0696 + IT_0324
      *IT_0699 + conj(IT_0079)*IT_0710 + IT_0079*IT_0713)*IT_0760 + (IT_0618
      *IT_0674 + IT_0622*IT_0677 + conj(IT_0508)*IT_0686 + IT_0508*IT_0689 +
       conj(IT_0638)*IT_0696 + IT_0638*IT_0699 + conj(IT_0640)*IT_0710 + IT_0640
      *IT_0713)*IT_0761 + IT_0763*IT_0786 + IT_0787*IT_0788 + (-8)*IT_0793
      *IT_0804 + (-4)*IT_0806*IT_0816;
    return create_ccomplex_return(IT_0817);
}

