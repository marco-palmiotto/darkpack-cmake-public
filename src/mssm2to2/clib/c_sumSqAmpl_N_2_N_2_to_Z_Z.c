#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_Z_Z(
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = pow(m_N_2, 2);
    const ccomplex_t IT_0001 = s_14*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(s_23, 2);
    const ccomplex_t IT_0004 = s_14*IT_0003;
    const ccomplex_t IT_0005 = pow(m_Z, 2);
    const ccomplex_t IT_0006 = s_12*s_24*IT_0005;
    const ccomplex_t IT_0007 = s_12*s_23*s_34;
    const ccomplex_t IT_0008 = s_14*IT_0000*IT_0005;
    const ccomplex_t IT_0009 = 0.5*IT_0006 + -IT_0007 + (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = IT_0004 + IT_0009;
    const ccomplex_t IT_0011 = IT_0002*IT_0010;
    const ccomplex_t IT_0012 = s_12*s_24;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = IT_0001 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0005 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0017, -1);
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = IT_0020*IT_0026;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0022 + IT_0025 + 
      -IT_0027 + -IT_0028);
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0029;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = IT_0016*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0000 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0037 = IT_0024*IT_0036;
    const ccomplex_t IT_0038 = IT_0020*IT_0036;
    const ccomplex_t IT_0039 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0040 = IT_0024*IT_0039;
    const ccomplex_t IT_0041 = IT_0020*IT_0039;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0037 + IT_0038 + 
      -IT_0040 + -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0045 = IT_0020*IT_0044;
    const ccomplex_t IT_0046 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0047 = IT_0020*IT_0046;
    const ccomplex_t IT_0048 = IT_0024*IT_0044;
    const ccomplex_t IT_0049 = IT_0024*IT_0046;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0045 + -IT_0047 +
       IT_0048 + -IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0043*IT_0052;
    const ccomplex_t IT_0054 = IT_0035*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0000 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0058 = IT_0024*IT_0057;
    const ccomplex_t IT_0059 = IT_0020*IT_0057;
    const ccomplex_t IT_0060 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0061 = IT_0020*IT_0060;
    const ccomplex_t IT_0062 = IT_0024*IT_0060;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0058 + IT_0059 + 
      -IT_0061 + -IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0067 = IT_0024*IT_0066;
    const ccomplex_t IT_0068 = IT_0020*IT_0066;
    const ccomplex_t IT_0069 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0070 = IT_0020*IT_0069;
    const ccomplex_t IT_0071 = IT_0024*IT_0069;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0067 + IT_0068 + 
      -IT_0070 + -IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0065*IT_0076;
    const ccomplex_t IT_0078 = IT_0056*IT_0077;
    const ccomplex_t IT_0079 = -IT_0033 + (-2)*IT_0054 + (-2)*IT_0078;
    const ccomplex_t IT_0080 = cpow(IT_0031, 2);
    const ccomplex_t IT_0081 = IT_0016*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_13 + IT_0000 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0085 = IT_0020*IT_0084;
    const ccomplex_t IT_0086 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0087 = IT_0024*IT_0086;
    const ccomplex_t IT_0088 = IT_0020*IT_0086;
    const ccomplex_t IT_0089 = IT_0024*IT_0084;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0085 + -IT_0087 + 
      -IT_0088 + IT_0089);
    const ccomplex_t IT_0091 = (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0093 = IT_0024*IT_0092;
    const ccomplex_t IT_0094 = IT_0020*IT_0092;
    const ccomplex_t IT_0095 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0096 = IT_0024*IT_0095;
    const ccomplex_t IT_0097 = IT_0020*IT_0095;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0093 + IT_0094 + 
      -IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0091*IT_0100;
    const ccomplex_t IT_0102 = IT_0083*IT_0101;
    const ccomplex_t IT_0103 = IT_0081 + (-2)*IT_0102;
    const ccomplex_t IT_0104 = conj(IT_0079) + conj(IT_0103);
    const ccomplex_t IT_0105 = cpow((-2)*s_23 + IT_0000 + IT_0005 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0077*IT_0106;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_23 + IT_0005 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = IT_0032*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_23 + IT_0000 + IT_0005 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0101*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0111 + -IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_23 + IT_0000 + IT_0005 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0053*IT_0117;
    const ccomplex_t IT_0119 = -IT_0118;
    const ccomplex_t IT_0120 = cpow(IT_0030, 2);
    const ccomplex_t IT_0121 = IT_0110*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = IT_0119 + IT_0122;
    const ccomplex_t IT_0124 = IT_0108 + IT_0115 + IT_0123;
    const ccomplex_t IT_0125 = IT_0079 + IT_0103;
    const ccomplex_t IT_0126 = conj(IT_0119) + conj(IT_0122);
    const ccomplex_t IT_0127 = conj(IT_0108) + conj(IT_0115) + IT_0126;
    const ccomplex_t IT_0128 = IT_0016*IT_0120;
    const ccomplex_t IT_0129 = (-2)*IT_0102 + IT_0128;
    const ccomplex_t IT_0130 = conj(IT_0079) + conj(IT_0129);
    const ccomplex_t IT_0131 = -IT_0107 + (-0.5)*IT_0111 + -IT_0114;
    const ccomplex_t IT_0132 = IT_0080*IT_0110;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0119 + IT_0133;
    const ccomplex_t IT_0135 = IT_0131 + IT_0134;
    const ccomplex_t IT_0136 = IT_0079 + IT_0129;
    const ccomplex_t IT_0137 = conj(IT_0119) + conj(IT_0133);
    const ccomplex_t IT_0138 = conj(IT_0131) + IT_0137;
    const ccomplex_t IT_0139 = s_13*s_23*s_24;
    const ccomplex_t IT_0140 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0141 = -IT_0008 + (-2)*IT_0139 + 2*IT_0140;
    const ccomplex_t IT_0142 = IT_0006 + IT_0141;
    const ccomplex_t IT_0143 = IT_0002*IT_0142;
    const ccomplex_t IT_0144 = (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = -IT_0001;
    const ccomplex_t IT_0146 = IT_0012 + IT_0144 + IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0033 + IT_0054 + IT_0078;
    const ccomplex_t IT_0148 = (-0.5)*IT_0081;
    const ccomplex_t IT_0149 = IT_0102 + IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = conj(IT_0102) + conj(IT_0147) + conj(IT_0148);
    const ccomplex_t IT_0151 = IT_0102 + (-0.5)*IT_0128;
    const ccomplex_t IT_0152 = IT_0147 + IT_0151;
    const ccomplex_t IT_0153 = conj(IT_0147) + conj(IT_0151);
    const ccomplex_t IT_0154 = pow(m_Z, -4);
    const ccomplex_t IT_0155 = s_34*IT_0154;
    const ccomplex_t IT_0156 = s_23*IT_0005;
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = s_24*IT_0005;
    const ccomplex_t IT_0159 = IT_0002*IT_0158;
    const ccomplex_t IT_0160 = s_23*s_34;
    const ccomplex_t IT_0161 = (-2)*IT_0160;
    const ccomplex_t IT_0162 = IT_0158 + IT_0161;
    const ccomplex_t IT_0163 = IT_0002*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0157 + 0.5*IT_0159 + (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = s_24 + IT_0164;
    const ccomplex_t IT_0166 = m_N_2*IT_0165;
    const ccomplex_t IT_0167 = (-2)*IT_0166;
    const ccomplex_t IT_0168 = cos(alpha);
    const ccomplex_t IT_0169 = IT_0023*IT_0168;
    const ccomplex_t IT_0170 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0171 = IT_0169*IT_0170;
    const ccomplex_t IT_0172 = sin(alpha);
    const ccomplex_t IT_0173 = IT_0023*IT_0172;
    const ccomplex_t IT_0174 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = IT_0019*IT_0168;
    const ccomplex_t IT_0177 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0178 = IT_0176*IT_0177;
    const ccomplex_t IT_0179 = IT_0019*IT_0172;
    const ccomplex_t IT_0180 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0171 + -IT_0175 + 
      -IT_0178 + IT_0181);
    const ccomplex_t IT_0183 = cos(beta);
    const ccomplex_t IT_0184 = cpow(IT_0017, -4);
    const ccomplex_t IT_0185 = IT_0018*IT_0168*IT_0183*IT_0184;
    const ccomplex_t IT_0186 = sin(beta);
    const ccomplex_t IT_0187 = IT_0018*IT_0172*IT_0184*IT_0186;
    const ccomplex_t IT_0188 = cpow(IT_0018, -3);
    const ccomplex_t IT_0189 = IT_0168*IT_0183*IT_0188;
    const ccomplex_t IT_0190 = IT_0172*IT_0186*IT_0188;
    const ccomplex_t IT_0191 = cpow(IT_0017, -2);
    const ccomplex_t IT_0192 = IT_0019*IT_0168*IT_0183*IT_0191;
    const ccomplex_t IT_0193 = IT_0019*IT_0172*IT_0186*IT_0191;
    const ccomplex_t IT_0194 = cpow(IT_0018, -2);
    const ccomplex_t IT_0195 = cpow(IT_0191 + IT_0194, -1);
    const ccomplex_t IT_0196 = (IT_0185 + IT_0187 + IT_0189 + IT_0190 + 2
      *IT_0192 + 2*IT_0193)*IT_0195;
    const ccomplex_t IT_0197 = m_W*e_em;
    const ccomplex_t IT_0198 = IT_0196*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0201 = IT_0199*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0182*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = 2*IT_0204;
    const ccomplex_t IT_0206 = IT_0169*IT_0174;
    const ccomplex_t IT_0207 = IT_0170*IT_0173;
    const ccomplex_t IT_0208 = IT_0176*IT_0180;
    const ccomplex_t IT_0209 = IT_0177*IT_0179;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(IT_0206 + IT_0207 + 
      -IT_0208 + -IT_0209);
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0018*IT_0168*IT_0184*IT_0186;
    const ccomplex_t IT_0213 = IT_0018*IT_0172*IT_0183*IT_0184;
    const ccomplex_t IT_0214 = IT_0168*IT_0186*IT_0188;
    const ccomplex_t IT_0215 = IT_0172*IT_0183*IT_0188;
    const ccomplex_t IT_0216 = IT_0019*IT_0168*IT_0186*IT_0191;
    const ccomplex_t IT_0217 = IT_0019*IT_0172*IT_0183*IT_0191;
    const ccomplex_t IT_0218 = IT_0195*(IT_0212 + -IT_0213 + IT_0214 + 
      -IT_0215 + 2*IT_0216 + (-2)*IT_0217);
    const ccomplex_t IT_0219 = IT_0197*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0222 = IT_0220*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = IT_0211*IT_0223;
    const ccomplex_t IT_0225 = (-0.5)*IT_0203 + -IT_0224;
    const ccomplex_t IT_0226 = (-2)*IT_0225;
    const ccomplex_t IT_0227 = IT_0205 + IT_0226;
    const ccomplex_t IT_0228 = 2*conj(IT_0204);
    const ccomplex_t IT_0229 = (-2)*conj(IT_0225);
    const ccomplex_t IT_0230 = IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0232 = IT_0169*IT_0231;
    const ccomplex_t IT_0233 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0234 = IT_0173*IT_0233;
    const ccomplex_t IT_0235 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0236 = IT_0176*IT_0235;
    const ccomplex_t IT_0237 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0238 = IT_0179*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0232 + -IT_0234 + 
      -IT_0236 + IT_0238);
    const ccomplex_t IT_0240 = IT_0202*IT_0239;
    const ccomplex_t IT_0241 = 0.5*IT_0240;
    const ccomplex_t IT_0242 = 2*IT_0241;
    const ccomplex_t IT_0243 = IT_0169*IT_0233;
    const ccomplex_t IT_0244 = IT_0173*IT_0231;
    const ccomplex_t IT_0245 = IT_0176*IT_0237;
    const ccomplex_t IT_0246 = IT_0179*IT_0235;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*(IT_0243 + IT_0244 + 
      -IT_0245 + -IT_0246);
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = IT_0223*IT_0248;
    const ccomplex_t IT_0250 = (-0.5)*IT_0240 + -IT_0249;
    const ccomplex_t IT_0251 = (-2)*IT_0250;
    const ccomplex_t IT_0252 = IT_0242 + IT_0251;
    const ccomplex_t IT_0253 = 2*conj(IT_0241);
    const ccomplex_t IT_0254 = (-2)*conj(IT_0250);
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = (conj(IT_0119) + conj(IT_0131) + conj(IT_0133) 
      + conj(IT_0147) + conj(IT_0151))*IT_0227 + (IT_0119 + IT_0131 + IT_0133 +
       IT_0147 + IT_0151)*IT_0230 + (conj(IT_0102) + conj(IT_0108) + conj
      (IT_0115) + conj(IT_0119) + conj(IT_0122) + conj(IT_0147) + conj(IT_0148))
      *IT_0252 + (IT_0102 + IT_0108 + IT_0115 + IT_0119 + IT_0122 + IT_0147 +
       IT_0148)*IT_0255;
    const ccomplex_t IT_0257 = s_14*s_23*s_34*IT_0154;
    const ccomplex_t IT_0258 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0259 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0260 = -IT_0258 + -IT_0259;
    const ccomplex_t IT_0261 = s_12 + IT_0257 + IT_0260;
    const ccomplex_t IT_0262 = m_N_2*IT_0261;
    const ccomplex_t IT_0263 = 2*IT_0118 + -IT_0132;
    const ccomplex_t IT_0264 = 2*IT_0107 + IT_0111 + 2*IT_0114;
    const ccomplex_t IT_0265 = conj(IT_0263) + conj(IT_0264);
    const ccomplex_t IT_0266 = IT_0263 + IT_0264;
    const ccomplex_t IT_0267 = 2*IT_0118 + -IT_0121;
    const ccomplex_t IT_0268 = conj(IT_0264) + conj(IT_0267);
    const ccomplex_t IT_0269 = IT_0264 + IT_0267;
    const ccomplex_t IT_0270 = s_13*s_14*s_34*IT_0154;
    const ccomplex_t IT_0271 = IT_0000 + IT_0270;
    const ccomplex_t IT_0272 = pow(s_13, 2);
    const ccomplex_t IT_0273 = IT_0002*IT_0272;
    const ccomplex_t IT_0274 = pow(s_14, 2);
    const ccomplex_t IT_0275 = IT_0002*IT_0274;
    const ccomplex_t IT_0276 = -IT_0273 + -IT_0275;
    const ccomplex_t IT_0277 = IT_0271 + IT_0276;
    const ccomplex_t IT_0278 = m_N_2*IT_0277;
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = s_13*IT_0005;
    const ccomplex_t IT_0281 = IT_0155*IT_0280;
    const ccomplex_t IT_0282 = s_14*IT_0005;
    const ccomplex_t IT_0283 = IT_0002*IT_0282;
    const ccomplex_t IT_0284 = s_13*s_34;
    const ccomplex_t IT_0285 = (-2)*IT_0284;
    const ccomplex_t IT_0286 = IT_0282 + IT_0285;
    const ccomplex_t IT_0287 = IT_0002*IT_0286;
    const ccomplex_t IT_0288 = (-0.5)*IT_0281 + 0.5*IT_0283 + (-0.5)*IT_0287;
    const ccomplex_t IT_0289 = s_14 + IT_0288;
    const ccomplex_t IT_0290 = m_N_2*IT_0289;
    const ccomplex_t IT_0291 = (-2)*IT_0290;
    const ccomplex_t IT_0292 = 2*IT_0225;
    const ccomplex_t IT_0293 = (-2)*IT_0204;
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = 2*conj(IT_0225);
    const ccomplex_t IT_0296 = (-2)*conj(IT_0204);
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = 2*IT_0250;
    const ccomplex_t IT_0299 = (-2)*IT_0241;
    const ccomplex_t IT_0300 = IT_0298 + IT_0299;
    const ccomplex_t IT_0301 = 2*conj(IT_0250);
    const ccomplex_t IT_0302 = (-2)*conj(IT_0241);
    const ccomplex_t IT_0303 = IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = (conj(IT_0102) + conj(IT_0108) + conj(IT_0115) 
      + conj(IT_0119) + conj(IT_0122) + conj(IT_0147) + conj(IT_0148))*IT_0294 +
       (IT_0102 + IT_0108 + IT_0115 + IT_0119 + IT_0122 + IT_0147 + IT_0148)
      *IT_0297 + (conj(IT_0119) + conj(IT_0131) + conj(IT_0133) + conj(IT_0147) 
      + conj(IT_0151))*IT_0300 + (IT_0119 + IT_0131 + IT_0133 + IT_0147 +
       IT_0151)*IT_0303;
    const ccomplex_t IT_0305 = s_23*s_24*s_34*IT_0154;
    const ccomplex_t IT_0306 = IT_0000 + IT_0305;
    const ccomplex_t IT_0307 = IT_0002*IT_0003;
    const ccomplex_t IT_0308 = pow(s_24, 2);
    const ccomplex_t IT_0309 = IT_0002*IT_0308;
    const ccomplex_t IT_0310 = -IT_0307 + -IT_0309;
    const ccomplex_t IT_0311 = IT_0306 + IT_0310;
    const ccomplex_t IT_0312 = m_N_2*IT_0311;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = s_13*s_24*s_34*IT_0154;
    const ccomplex_t IT_0315 = s_12 + IT_0260 + IT_0314;
    const ccomplex_t IT_0316 = m_N_2*IT_0315;
    const ccomplex_t IT_0317 = m_N_2*IT_0002;
    const ccomplex_t IT_0318 = pow(m_Z, 4);
    const ccomplex_t IT_0319 = s_24*IT_0318;
    const ccomplex_t IT_0320 = IT_0002*IT_0319;
    const ccomplex_t IT_0321 = (-2)*IT_0160 + (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = IT_0158 + IT_0321;
    const ccomplex_t IT_0323 = IT_0317*IT_0322;
    const ccomplex_t IT_0324 = (-2)*IT_0323;
    const ccomplex_t IT_0325 = cpow(IT_0100, 2);
    const ccomplex_t IT_0326 = m_N_3*IT_0325;
    const ccomplex_t IT_0327 = IT_0083*IT_0326;
    const ccomplex_t IT_0328 = cpow(IT_0065, 2);
    const ccomplex_t IT_0329 = m_N_4*IT_0328;
    const ccomplex_t IT_0330 = m_N_2*IT_0077;
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = IT_0329 + IT_0331;
    const ccomplex_t IT_0333 = IT_0056*IT_0332;
    const ccomplex_t IT_0334 = IT_0327 + IT_0333;
    const ccomplex_t IT_0335 = 2*IT_0334;
    const ccomplex_t IT_0336 = m_N_2*IT_0120;
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = m_N_2*IT_0032;
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = IT_0016*IT_0339;
    const ccomplex_t IT_0341 = cpow(IT_0043, 2);
    const ccomplex_t IT_0342 = m_N_1*IT_0341;
    const ccomplex_t IT_0343 = m_N_2*IT_0053;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = IT_0342 + IT_0344;
    const ccomplex_t IT_0346 = IT_0035*IT_0345;
    const ccomplex_t IT_0347 = (-0.5)*IT_0340 + IT_0346;
    const ccomplex_t IT_0348 = 2*IT_0347;
    const ccomplex_t IT_0349 = IT_0335 + IT_0348;
    const ccomplex_t IT_0350 = m_N_2*IT_0080;
    const ccomplex_t IT_0351 = -IT_0350;
    const ccomplex_t IT_0352 = IT_0338 + IT_0351;
    const ccomplex_t IT_0353 = IT_0016*IT_0352;
    const ccomplex_t IT_0354 = m_N_2*IT_0101;
    const ccomplex_t IT_0355 = IT_0083*IT_0354;
    const ccomplex_t IT_0356 = 0.5*IT_0353 + IT_0355;
    const ccomplex_t IT_0357 = (-2)*IT_0356;
    const ccomplex_t IT_0358 = IT_0349 + IT_0357;
    const ccomplex_t IT_0359 = 2*conj(IT_0334);
    const ccomplex_t IT_0360 = 2*conj(IT_0347);
    const ccomplex_t IT_0361 = IT_0359 + IT_0360;
    const ccomplex_t IT_0362 = (-2)*conj(IT_0356);
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = cpow(IT_0052, 2);
    const ccomplex_t IT_0365 = m_N_1*IT_0364;
    const ccomplex_t IT_0366 = IT_0344 + IT_0365;
    const ccomplex_t IT_0367 = IT_0035*IT_0366;
    const ccomplex_t IT_0368 = cpow(IT_0091, 2);
    const ccomplex_t IT_0369 = m_N_3*IT_0368;
    const ccomplex_t IT_0370 = IT_0083*IT_0369;
    const ccomplex_t IT_0371 = cpow(IT_0076, 2);
    const ccomplex_t IT_0372 = m_N_4*IT_0371;
    const ccomplex_t IT_0373 = IT_0331 + IT_0372;
    const ccomplex_t IT_0374 = IT_0056*IT_0373;
    const ccomplex_t IT_0375 = (-0.5)*IT_0353 + IT_0367 + IT_0370 + IT_0374;
    const ccomplex_t IT_0376 = 2*IT_0375;
    const ccomplex_t IT_0377 = (-2)*IT_0355;
    const ccomplex_t IT_0378 = -IT_0016*(IT_0336 + -IT_0338);
    const ccomplex_t IT_0379 = 0.5*IT_0378;
    const ccomplex_t IT_0380 = (-2)*IT_0379;
    const ccomplex_t IT_0381 = IT_0376 + IT_0377 + IT_0380;
    const ccomplex_t IT_0382 = 2*conj(IT_0375);
    const ccomplex_t IT_0383 = (-2)*conj(IT_0355);
    const ccomplex_t IT_0384 = (-2)*conj(IT_0379);
    const ccomplex_t IT_0385 = IT_0382 + IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0387 = s_14*IT_0318;
    const ccomplex_t IT_0388 = s_13*s_34*IT_0005;
    const ccomplex_t IT_0389 = (-2)*IT_0388;
    const ccomplex_t IT_0390 = IT_0387 + IT_0389;
    const ccomplex_t IT_0391 = IT_0154*IT_0390;
    const ccomplex_t IT_0392 = (-0.666666666666667)*IT_0386 + (
      -0.166666666666667)*IT_0391;
    const ccomplex_t IT_0393 = s_14 + IT_0392;
    const ccomplex_t IT_0394 = m_N_2*IT_0393;
    const ccomplex_t IT_0395 = (-6)*IT_0394;
    const ccomplex_t IT_0396 = (-0.666666666666667)*IT_0287 +
       0.166666666666667*IT_0391;
    const ccomplex_t IT_0397 = s_14 + IT_0396;
    const ccomplex_t IT_0398 = m_N_2*IT_0397;
    const ccomplex_t IT_0399 = 6*IT_0398;
    const ccomplex_t IT_0400 = s_23*s_34*IT_0005;
    const ccomplex_t IT_0401 = (-2)*IT_0400;
    const ccomplex_t IT_0402 = IT_0319 + IT_0401;
    const ccomplex_t IT_0403 = IT_0154*IT_0402;
    const ccomplex_t IT_0404 = (-0.666666666666667)*IT_0163 +
       0.166666666666667*IT_0403;
    const ccomplex_t IT_0405 = s_24 + IT_0404;
    const ccomplex_t IT_0406 = m_N_2*IT_0405;
    const ccomplex_t IT_0407 = 6*IT_0406;
    const ccomplex_t IT_0408 = IT_0110*IT_0339;
    const ccomplex_t IT_0409 = -IT_0354;
    const ccomplex_t IT_0410 = IT_0326 + IT_0409;
    const ccomplex_t IT_0411 = IT_0113*IT_0410;
    const ccomplex_t IT_0412 = IT_0106*IT_0332;
    const ccomplex_t IT_0413 = 0.5*IT_0408 + -IT_0411 + -IT_0412;
    const ccomplex_t IT_0414 = 2*IT_0413;
    const ccomplex_t IT_0415 = IT_0117*(IT_0342 + -IT_0343);
    const ccomplex_t IT_0416 = (-2)*IT_0415;
    const ccomplex_t IT_0417 = IT_0110*(IT_0338 + -IT_0350);
    const ccomplex_t IT_0418 = (-0.5)*IT_0417;
    const ccomplex_t IT_0419 = (-2)*IT_0418;
    const ccomplex_t IT_0420 = IT_0414 + IT_0416 + IT_0419;
    const ccomplex_t IT_0421 = 2*conj(IT_0413);
    const ccomplex_t IT_0422 = (-2)*conj(IT_0415);
    const ccomplex_t IT_0423 = (-2)*conj(IT_0418);
    const ccomplex_t IT_0424 = IT_0421 + IT_0422 + IT_0423;
    const ccomplex_t IT_0425 = IT_0110*IT_0352;
    const ccomplex_t IT_0426 = IT_0369 + IT_0409;
    const ccomplex_t IT_0427 = IT_0113*IT_0426;
    const ccomplex_t IT_0428 = IT_0106*IT_0373;
    const ccomplex_t IT_0429 = 0.5*IT_0425 + -IT_0427 + -IT_0428;
    const ccomplex_t IT_0430 = 2*IT_0429;
    const ccomplex_t IT_0431 = IT_0117*(IT_0343 + -IT_0365);
    const ccomplex_t IT_0432 = -IT_0431;
    const ccomplex_t IT_0433 = (-2)*IT_0432;
    const ccomplex_t IT_0434 = IT_0110*(IT_0336 + -IT_0338);
    const ccomplex_t IT_0435 = 0.5*IT_0434;
    const ccomplex_t IT_0436 = (-2)*IT_0435;
    const ccomplex_t IT_0437 = IT_0430 + IT_0433 + IT_0436;
    const ccomplex_t IT_0438 = 2*conj(IT_0429);
    const ccomplex_t IT_0439 = (-2)*conj(IT_0432);
    const ccomplex_t IT_0440 = (-2)*conj(IT_0435);
    const ccomplex_t IT_0441 = IT_0438 + IT_0439 + IT_0440;
    const ccomplex_t IT_0442 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0443 = (-0.166666666666667)*IT_0403 + (
      -0.666666666666667)*IT_0442;
    const ccomplex_t IT_0444 = s_24 + IT_0443;
    const ccomplex_t IT_0445 = m_N_2*IT_0444;
    const ccomplex_t IT_0446 = (-6)*IT_0445;
    const ccomplex_t IT_0447 = IT_0002*IT_0387;
    const ccomplex_t IT_0448 = (-2)*IT_0284 + (-0.5)*IT_0447;
    const ccomplex_t IT_0449 = IT_0282 + IT_0448;
    const ccomplex_t IT_0450 = IT_0317*IT_0449;
    const ccomplex_t IT_0451 = (-2)*IT_0450;
    const ccomplex_t IT_0452 = s_12*s_14;
    const ccomplex_t IT_0453 = s_12*s_14*IT_0005;
    const ccomplex_t IT_0454 = s_24*IT_0000*IT_0005;
    const ccomplex_t IT_0455 = s_24*IT_0272;
    const ccomplex_t IT_0456 = s_12*s_13*s_34;
    const ccomplex_t IT_0457 = -IT_0454 + 2*IT_0455 + (-2)*IT_0456;
    const ccomplex_t IT_0458 = IT_0453 + IT_0457;
    const ccomplex_t IT_0459 = IT_0002*IT_0458;
    const ccomplex_t IT_0460 = (-0.5)*IT_0459;
    const ccomplex_t IT_0461 = s_24*IT_0000;
    const ccomplex_t IT_0462 = -IT_0461;
    const ccomplex_t IT_0463 = IT_0452 + IT_0460 + IT_0462;
    const ccomplex_t IT_0464 = (-4)*IT_0267;
    const ccomplex_t IT_0465 = (-4)*IT_0264;
    const ccomplex_t IT_0466 = IT_0464 + IT_0465;
    const ccomplex_t IT_0467 = (-4)*conj(IT_0267);
    const ccomplex_t IT_0468 = (-4)*conj(IT_0264);
    const ccomplex_t IT_0469 = IT_0467 + IT_0468;
    const ccomplex_t IT_0470 = (-4)*IT_0263;
    const ccomplex_t IT_0471 = IT_0465 + IT_0470;
    const ccomplex_t IT_0472 = (-4)*conj(IT_0263);
    const ccomplex_t IT_0473 = IT_0468 + IT_0472;
    const ccomplex_t IT_0474 = IT_0154*IT_0387;
    const ccomplex_t IT_0475 = 0.666666666666667*IT_0283 + (-0.166666666666667
      )*IT_0474;
    const ccomplex_t IT_0476 = s_14 + IT_0475;
    const ccomplex_t IT_0477 = m_N_2*IT_0476;
    const ccomplex_t IT_0478 = (-6)*IT_0477;
    const ccomplex_t IT_0479 = IT_0416 + IT_0419;
    const ccomplex_t IT_0480 = IT_0422 + IT_0423;
    const ccomplex_t IT_0481 = IT_0433 + IT_0436;
    const ccomplex_t IT_0482 = IT_0439 + IT_0440;
    const ccomplex_t IT_0483 = s_13*s_14*s_23;
    const ccomplex_t IT_0484 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0485 = -IT_0454 + (-2)*IT_0483 + 2*IT_0484;
    const ccomplex_t IT_0486 = IT_0453 + IT_0485;
    const ccomplex_t IT_0487 = IT_0002*IT_0486;
    const ccomplex_t IT_0488 = 2*IT_0461;
    const ccomplex_t IT_0489 = (-2)*IT_0452;
    const ccomplex_t IT_0490 = IT_0487 + IT_0488 + IT_0489;
    const ccomplex_t IT_0491 = (-2)*IT_0267;
    const ccomplex_t IT_0492 = (-2)*IT_0264;
    const ccomplex_t IT_0493 = IT_0491 + IT_0492;
    const ccomplex_t IT_0494 = (-2)*conj(IT_0267);
    const ccomplex_t IT_0495 = (-2)*conj(IT_0264);
    const ccomplex_t IT_0496 = IT_0494 + IT_0495;
    const ccomplex_t IT_0497 = (-2)*IT_0263;
    const ccomplex_t IT_0498 = IT_0492 + IT_0497;
    const ccomplex_t IT_0499 = (-2)*conj(IT_0263);
    const ccomplex_t IT_0500 = IT_0495 + IT_0499;
    const ccomplex_t IT_0501 = IT_0154*IT_0319;
    const ccomplex_t IT_0502 = (-10)*s_24;
    const ccomplex_t IT_0503 = IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = m_N_2*IT_0503;
    const ccomplex_t IT_0505 = conj(IT_0102) + conj(IT_0148);
    const ccomplex_t IT_0506 = IT_0102 + IT_0148;
    const ccomplex_t IT_0507 = (-10)*s_14;
    const ccomplex_t IT_0508 = IT_0474 + IT_0507;
    const ccomplex_t IT_0509 = m_N_2*IT_0508;
    const ccomplex_t IT_0510 = conj(IT_0108) + conj(IT_0115);
    const ccomplex_t IT_0511 = IT_0108 + IT_0115;
    const ccomplex_t IT_0512 = s_12*IT_0005;
    const ccomplex_t IT_0513 = s_13*s_23;
    const ccomplex_t IT_0514 = s_14*s_24;
    const ccomplex_t IT_0515 = s_14*s_23*s_34*IT_0002;
    const ccomplex_t IT_0516 = (-2)*IT_0513 + -IT_0514 + 2*IT_0515;
    const ccomplex_t IT_0517 = IT_0512 + IT_0516;
    const ccomplex_t IT_0518 = IT_0002*IT_0517;
    const ccomplex_t IT_0519 = -IT_0259 + (-0.5)*IT_0518;
    const ccomplex_t IT_0520 = s_12 + IT_0519;
    const ccomplex_t IT_0521 = m_N_2*IT_0520;
    const ccomplex_t IT_0522 = (-2)*IT_0521;
    const ccomplex_t IT_0523 = 2*IT_0415;
    const ccomplex_t IT_0524 = 2*IT_0418;
    const ccomplex_t IT_0525 = IT_0523 + IT_0524;
    const ccomplex_t IT_0526 = 2*conj(IT_0415);
    const ccomplex_t IT_0527 = 2*conj(IT_0418);
    const ccomplex_t IT_0528 = IT_0526 + IT_0527;
    const ccomplex_t IT_0529 = 2*IT_0432;
    const ccomplex_t IT_0530 = 2*IT_0435;
    const ccomplex_t IT_0531 = IT_0529 + IT_0530;
    const ccomplex_t IT_0532 = 2*conj(IT_0432);
    const ccomplex_t IT_0533 = 2*conj(IT_0435);
    const ccomplex_t IT_0534 = IT_0532 + IT_0533;
    const ccomplex_t IT_0535 = conj(IT_0432) + conj(IT_0435);
    const ccomplex_t IT_0536 = 2*IT_0267;
    const ccomplex_t IT_0537 = IT_0432 + IT_0435;
    const ccomplex_t IT_0538 = 2*conj(IT_0267);
    const ccomplex_t IT_0539 = IT_0522*(conj(IT_0429)*IT_0493 + IT_0429
      *IT_0496 + conj(IT_0413)*IT_0498 + IT_0413*IT_0500 + IT_0265*IT_0525 +
       IT_0266*IT_0528 + conj(IT_0264)*IT_0531 + IT_0264*IT_0534 + IT_0535
      *IT_0536 + IT_0537*IT_0538);
    const ccomplex_t IT_0540 = IT_0000*IT_0005;
    const ccomplex_t IT_0541 = -IT_0540;
    const ccomplex_t IT_0542 = IT_0308 + IT_0541;
    const ccomplex_t IT_0543 = IT_0002*IT_0542;
    const ccomplex_t IT_0544 = -IT_0000;
    const ccomplex_t IT_0545 = IT_0309 + IT_0544;
    const ccomplex_t IT_0546 = (-4)*IT_0545;
    const ccomplex_t IT_0547 = IT_0543 + IT_0546;
    const ccomplex_t IT_0548 = m_N_2*IT_0547;
    const ccomplex_t IT_0549 = conj(IT_0355) + conj(IT_0379);
    const ccomplex_t IT_0550 = (-2)*IT_0079;
    const ccomplex_t IT_0551 = (-2)*IT_0103;
    const ccomplex_t IT_0552 = IT_0550 + IT_0551;
    const ccomplex_t IT_0553 = IT_0355 + IT_0379;
    const ccomplex_t IT_0554 = (-2)*conj(IT_0079);
    const ccomplex_t IT_0555 = (-2)*conj(IT_0103);
    const ccomplex_t IT_0556 = IT_0554 + IT_0555;
    const ccomplex_t IT_0557 = (-2)*IT_0129;
    const ccomplex_t IT_0558 = IT_0550 + IT_0557;
    const ccomplex_t IT_0559 = (-2)*conj(IT_0129);
    const ccomplex_t IT_0560 = IT_0554 + IT_0559;
    const ccomplex_t IT_0561 = 2*IT_0079;
    const ccomplex_t IT_0562 = 2*conj(IT_0079);
    const ccomplex_t IT_0563 = IT_0548*(IT_0130*IT_0349 + IT_0136*IT_0361 +
       conj(IT_0103)*IT_0376 + IT_0103*IT_0382 + IT_0549*IT_0552 + IT_0553
      *IT_0556 + conj(IT_0356)*IT_0558 + IT_0356*IT_0560 + conj(IT_0375)*IT_0561
       + IT_0375*IT_0562);
    const ccomplex_t IT_0564 = (-2)*IT_0513;
    const ccomplex_t IT_0565 = IT_0512 + IT_0564;
    const ccomplex_t IT_0566 = IT_0002*IT_0565;
    const ccomplex_t IT_0567 = (s_12 + -IT_0259)*(s_12 + (-0.5)*IT_0566);
    const ccomplex_t IT_0568 = (-2)*IT_0567;
    const ccomplex_t IT_0569 = 2*IT_0263;
    const ccomplex_t IT_0570 = 2*conj(IT_0263);
    const ccomplex_t IT_0571 = 2*IT_0264;
    const ccomplex_t IT_0572 = 2*conj(IT_0264);
    const ccomplex_t IT_0573 = IT_0104*IT_0536 + IT_0125*IT_0538 + 2*conj
      (IT_0264)*IT_0561 + 2*IT_0264*IT_0562 + IT_0130*IT_0569 + IT_0136*IT_0570 
      + conj(IT_0103)*IT_0571 + conj(IT_0129)*IT_0571 + IT_0103*IT_0572 +
       IT_0129*IT_0572;
    const ccomplex_t IT_0574 = s_13*s_14*s_34*IT_0002;
    const ccomplex_t IT_0575 = 0.5*IT_0274 + (-0.5)*IT_0540 + -IT_0574;
    const ccomplex_t IT_0576 = IT_0272 + IT_0575;
    const ccomplex_t IT_0577 = IT_0002*IT_0576;
    const ccomplex_t IT_0578 = -IT_0275;
    const ccomplex_t IT_0579 = IT_0000 + IT_0577 + IT_0578;
    const ccomplex_t IT_0580 = m_N_2*IT_0579;
    const ccomplex_t IT_0581 = conj(IT_0334) + conj(IT_0347);
    const ccomplex_t IT_0582 = 4*IT_0267;
    const ccomplex_t IT_0583 = IT_0334 + IT_0347;
    const ccomplex_t IT_0584 = 4*conj(IT_0267);
    const ccomplex_t IT_0585 = IT_0265*IT_0375 + IT_0266*conj(IT_0375) + 0.25
      *conj(IT_0356)*IT_0466 + 0.25*IT_0356*IT_0469 + 0.25*IT_0471*IT_0549 +
       0.25*IT_0473*IT_0553 + IT_0581*(IT_0264 + 0.25*IT_0582) + IT_0583*(conj
      (IT_0264) + 0.25*IT_0584);
    const ccomplex_t IT_0586 = IT_0000 + IT_0578;
    const ccomplex_t IT_0587 = IT_0000*IT_0586;
    const ccomplex_t IT_0588 = -IT_0587;
    const ccomplex_t IT_0589 = 6*conj(IT_0267);
    const ccomplex_t IT_0590 = -IT_0309;
    const ccomplex_t IT_0591 = IT_0000 + IT_0590;
    const ccomplex_t IT_0592 = IT_0000*IT_0591;
    const ccomplex_t IT_0593 = -IT_0592;
    const ccomplex_t IT_0594 = 6*IT_0079;
    const ccomplex_t IT_0595 = 6*conj(IT_0079);
    const ccomplex_t IT_0596 = 1.125*IT_0000;
    const ccomplex_t IT_0597 = (-16)*IT_0415;
    const ccomplex_t IT_0598 = (-16)*IT_0418;
    const ccomplex_t IT_0599 = (-16)*conj(IT_0415);
    const ccomplex_t IT_0600 = (-16)*conj(IT_0418);
    const ccomplex_t IT_0601 = 16*IT_0347;
    const ccomplex_t IT_0602 = (-16)*IT_0356;
    const ccomplex_t IT_0603 = 16*conj(IT_0334);
    const ccomplex_t IT_0604 = 16*conj(IT_0347);
    const ccomplex_t IT_0605 = (-16)*conj(IT_0356);
    const ccomplex_t IT_0606 = ((-16)*conj(IT_0334) + (-16)*conj(IT_0347) + 16
      *conj(IT_0356))*IT_0375 + ((-16)*IT_0347 + 16*IT_0356)*conj(IT_0375) +
       IT_0334*(16*conj(IT_0355) + (-16)*conj(IT_0375) + 16*conj(IT_0379)) + ((
      -16)*conj(IT_0413) + 16*conj(IT_0415) + 16*conj(IT_0418))*IT_0429 + ((-16)
      *IT_0413 + 16*IT_0415 + 16*IT_0418)*conj(IT_0429) + IT_0535*(16*IT_0413 +
       IT_0597 + IT_0598) + IT_0537*(16*conj(IT_0413) + IT_0599 + IT_0600) +
       IT_0549*(IT_0601 + IT_0602) + IT_0553*(IT_0603 + IT_0604 + IT_0605);
    const ccomplex_t IT_0607 = (-2)*IT_0514;
    const ccomplex_t IT_0608 = IT_0512 + IT_0607;
    const ccomplex_t IT_0609 = IT_0002*IT_0608;
    const ccomplex_t IT_0610 = (-0.5)*IT_0609;
    const ccomplex_t IT_0611 = (-0.5)*IT_0566;
    const ccomplex_t IT_0612 = s_12*IT_0318;
    const ccomplex_t IT_0613 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0614 = s_14*s_24*IT_0005;
    const ccomplex_t IT_0615 = s_14*s_23*s_34;
    const ccomplex_t IT_0616 = (-2)*IT_0613 + (-2)*IT_0614 + 4*IT_0615;
    const ccomplex_t IT_0617 = IT_0612 + IT_0616;
    const ccomplex_t IT_0618 = IT_0154*IT_0617;
    const ccomplex_t IT_0619 = 0.25*IT_0618;
    const ccomplex_t IT_0620 = 8*IT_0413;
    const ccomplex_t IT_0621 = (-8)*IT_0415;
    const ccomplex_t IT_0622 = (-8)*IT_0418;
    const ccomplex_t IT_0623 = IT_0620 + IT_0621 + IT_0622;
    const ccomplex_t IT_0624 = 8*IT_0429;
    const ccomplex_t IT_0625 = (-8)*IT_0432;
    const ccomplex_t IT_0626 = (-8)*IT_0435;
    const ccomplex_t IT_0627 = IT_0624 + IT_0625 + IT_0626;
    const ccomplex_t IT_0628 = (-8)*conj(IT_0432);
    const ccomplex_t IT_0629 = (-8)*conj(IT_0435);
    const ccomplex_t IT_0630 = IT_0628 + IT_0629;
    const ccomplex_t IT_0631 = conj(IT_0415) + conj(IT_0418);
    const ccomplex_t IT_0632 = 8*IT_0415;
    const ccomplex_t IT_0633 = 8*conj(IT_0415);
    const ccomplex_t IT_0634 = 8*IT_0418;
    const ccomplex_t IT_0635 = 8*IT_0432;
    const ccomplex_t IT_0636 = 8*conj(IT_0432);
    const ccomplex_t IT_0637 = (-8)*conj(IT_0415);
    const ccomplex_t IT_0638 = (-8)*conj(IT_0418);
    const ccomplex_t IT_0639 = conj(IT_0413)*IT_0623 + conj(IT_0429)*IT_0627 +
       IT_0429*IT_0630 + IT_0631*IT_0632 + IT_0418*IT_0633 + conj(IT_0418)
      *IT_0634 + IT_0535*IT_0635 + 8*IT_0435*(conj(IT_0435) + 0.125*IT_0636) +
       IT_0413*(IT_0637 + IT_0638);
    const ccomplex_t IT_0640 = s_13*s_24*s_34;
    const ccomplex_t IT_0641 = (-2)*IT_0613 + (-2)*IT_0614 + 4*IT_0640;
    const ccomplex_t IT_0642 = IT_0612 + IT_0641;
    const ccomplex_t IT_0643 = IT_0154*IT_0642;
    const ccomplex_t IT_0644 = 0.25*IT_0643;
    const ccomplex_t IT_0645 = 8*IT_0356;
    const ccomplex_t IT_0646 = (-8)*IT_0334;
    const ccomplex_t IT_0647 = (-8)*IT_0347;
    const ccomplex_t IT_0648 = IT_0645 + IT_0646 + IT_0647;
    const ccomplex_t IT_0649 = (-8)*conj(IT_0334);
    const ccomplex_t IT_0650 = (-8)*conj(IT_0347);
    const ccomplex_t IT_0651 = IT_0649 + IT_0650;
    const ccomplex_t IT_0652 = 8*IT_0334;
    const ccomplex_t IT_0653 = 8*conj(IT_0334);
    const ccomplex_t IT_0654 = 8*IT_0347;
    const ccomplex_t IT_0655 = 8*IT_0375;
    const ccomplex_t IT_0656 = 8*IT_0355;
    const ccomplex_t IT_0657 = (-8)*IT_0375;
    const ccomplex_t IT_0658 = (-8)*conj(IT_0375);
    const ccomplex_t IT_0659 = 8*conj(IT_0355);
    const ccomplex_t IT_0660 = conj(IT_0356)*IT_0648 + IT_0356*IT_0651 +
       IT_0581*IT_0652 + IT_0347*IT_0653 + conj(IT_0347)*IT_0654 + conj(IT_0375)
      *IT_0655 + IT_0549*(IT_0656 + IT_0657) + IT_0355*IT_0658 + 8*IT_0379*(conj
      (IT_0379) + 0.125*IT_0658 + 0.125*IT_0659);
    const ccomplex_t IT_0661 = pow(s_34, 2);
    const ccomplex_t IT_0662 = s_12*IT_0661;
    const ccomplex_t IT_0663 = (-2)*IT_0613 + (-2)*IT_0614 + 2*IT_0615 + 2
      *IT_0640 + (-2)*IT_0662;
    const ccomplex_t IT_0664 = IT_0612 + IT_0663;
    const ccomplex_t IT_0665 = IT_0154*IT_0664;
    const ccomplex_t IT_0666 = (-0.25)*IT_0665;
    const ccomplex_t IT_0667 = -IT_0259;
    const ccomplex_t IT_0668 = -IT_0258;
    const ccomplex_t IT_0669 = 8*conj(IT_0413);
    const ccomplex_t IT_0670 = IT_0637 + IT_0638 + IT_0669;
    const ccomplex_t IT_0671 = IT_0415 + IT_0418;
    const ccomplex_t IT_0672 = 8*conj(IT_0347);
    const ccomplex_t IT_0673 = conj(IT_0356)*IT_0623 + conj(IT_0355)*(IT_0625 
      + IT_0626) + conj(IT_0379)*IT_0627 + IT_0553*IT_0630 + conj(IT_0347)
      *IT_0632 + IT_0347*IT_0633 + 8*conj(IT_0375)*(IT_0435 + 0.125*IT_0635) +
       IT_0375*IT_0636 + conj(IT_0413)*(IT_0646 + IT_0647) + IT_0413*IT_0651 +
       IT_0631*IT_0652 + conj(IT_0418)*IT_0654 + conj(IT_0435)*IT_0655 + 8*conj
      (IT_0429)*(IT_0379 + 0.125*IT_0656 + 0.125*IT_0657) + IT_0429*(IT_0658 +
       IT_0659) + IT_0356*IT_0670 + IT_0653*IT_0671 + IT_0418*IT_0672;
    const ccomplex_t IT_0674 = (-2)*IT_0661;
    const ccomplex_t IT_0675 = IT_0318 + IT_0674;
    const ccomplex_t IT_0676 = IT_0154*IT_0675;
    const ccomplex_t IT_0677 = 0.25*IT_0676;
    const ccomplex_t IT_0678 = 1 + IT_0677;
    const ccomplex_t IT_0679 = IT_0000*IT_0678;
    const ccomplex_t IT_0680 = IT_0549*IT_0623 + 8*conj(IT_0356)*(IT_0429 +
       0.125*IT_0625 + 0.125*IT_0626) + IT_0356*IT_0630 + 8*IT_0375*(conj
      (IT_0418) + 0.125*IT_0633) + conj(IT_0375)*(IT_0632 + IT_0634) + conj
      (IT_0429)*IT_0648 + IT_0429*IT_0651 + IT_0535*(IT_0652 + IT_0654) + conj
      (IT_0413)*IT_0657 + IT_0413*IT_0658 + IT_0553*IT_0670 + IT_0537*(IT_0653 +
       IT_0672);
    const ccomplex_t IT_0681 = pow(m_Z, 6);
    const ccomplex_t IT_0682 = s_12*IT_0681;
    const ccomplex_t IT_0683 = s_13*s_23*IT_0318;
    const ccomplex_t IT_0684 = (-2)*IT_0683;
    const ccomplex_t IT_0685 = IT_0682 + IT_0684;
    const ccomplex_t IT_0686 = IT_0154*IT_0685;
    const ccomplex_t IT_0687 = 0.25*IT_0686;
    const ccomplex_t IT_0688 = IT_0512 + IT_0687;
    const ccomplex_t IT_0689 = (-8)*IT_0119;
    const ccomplex_t IT_0690 = (-8)*conj(IT_0119);
    const ccomplex_t IT_0691 = (-8)*IT_0147;
    const ccomplex_t IT_0692 = IT_0108 + IT_0115 + IT_0122;
    const ccomplex_t IT_0693 = (-8)*conj(IT_0147);
    const ccomplex_t IT_0694 = IT_0131 + IT_0133;
    const ccomplex_t IT_0695 = conj(IT_0131) + conj(IT_0133);
    const ccomplex_t IT_0696 = (-8)*IT_0122;
    const ccomplex_t IT_0697 = (-8)*IT_0108;
    const ccomplex_t IT_0698 = (-8)*conj(IT_0122);
    const ccomplex_t IT_0699 = (-8)*conj(IT_0108);
    const ccomplex_t IT_0700 = 2*conj(IT_0147)*IT_0689 + 2*IT_0147*IT_0690 + 
      (conj(IT_0108) + conj(IT_0115) + conj(IT_0122))*IT_0691 + IT_0692*IT_0693 
      + IT_0693*IT_0694 + (-8)*conj(IT_0151)*((-0.125)*IT_0689 + IT_0694) +
       IT_0691*IT_0695 + (-8)*IT_0151*((-0.125)*IT_0690 + IT_0695) + (-8)
      *IT_0505*(IT_0115 + (-0.125)*IT_0689 + (-0.125)*IT_0696 + (-0.125)*IT_0697
      ) + (-8)*IT_0506*(conj(IT_0115) + (-0.125)*IT_0690 + (-0.125)*IT_0698 + (
      -0.125)*IT_0699);
    const ccomplex_t IT_0701 = s_12 + IT_0667;
    const ccomplex_t IT_0702 = IT_0000*IT_0701;
    const ccomplex_t IT_0703 = (-6)*conj(IT_0267);
    const ccomplex_t IT_0704 = (-6)*IT_0079;
    const ccomplex_t IT_0705 = (-6)*conj(IT_0079);
    const ccomplex_t IT_0706 = IT_0154*IT_0661;
    const ccomplex_t IT_0707 = 0.5*IT_0706;
    const ccomplex_t IT_0708 = 1 + IT_0707;
    const ccomplex_t IT_0709 = s_12*IT_0708;
    const ccomplex_t IT_0710 = (-4)*conj(IT_0204);
    const ccomplex_t IT_0711 = 0.25*IT_0710;
    const ccomplex_t IT_0712 = conj(IT_0225) + IT_0711;
    const ccomplex_t IT_0713 = (-4)*IT_0204;
    const ccomplex_t IT_0714 = (-4)*IT_0241;
    const ccomplex_t IT_0715 = (-4)*conj(IT_0241);
    const ccomplex_t IT_0716 = s_23*s_24*s_34*IT_0002;
    const ccomplex_t IT_0717 = 0.5*IT_0308 + (-0.5)*IT_0540 + -IT_0716;
    const ccomplex_t IT_0718 = IT_0003 + IT_0717;
    const ccomplex_t IT_0719 = IT_0002*IT_0718;
    const ccomplex_t IT_0720 = IT_0000 + IT_0590 + IT_0719;
    const ccomplex_t IT_0721 = m_N_2*IT_0720;
    const ccomplex_t IT_0722 = (-4)*IT_0415;
    const ccomplex_t IT_0723 = (-4)*IT_0418;
    const ccomplex_t IT_0724 = (-4)*conj(IT_0415);
    const ccomplex_t IT_0725 = (-4)*conj(IT_0418);
    const ccomplex_t IT_0726 = (-4)*IT_0079;
    const ccomplex_t IT_0727 = (-4)*conj(IT_0079);
    const ccomplex_t IT_0728 = conj(IT_0079)*IT_0429 + IT_0079*conj(IT_0429) +
       IT_0103*(conj(IT_0429) + -IT_0535) + conj(IT_0103)*(IT_0429 + -IT_0537) +
       IT_0130*(IT_0413 + 0.25*IT_0722 + 0.25*IT_0723) + IT_0136*(conj(IT_0413) 
      + 0.25*IT_0724 + 0.25*IT_0725) + 0.25*IT_0535*IT_0726 + 0.25*IT_0537
      *IT_0727;
    const ccomplex_t IT_0729 = 0.5*IT_0686;
    const ccomplex_t IT_0730 = (-4)*IT_0514;
    const ccomplex_t IT_0731 = IT_0512 + IT_0729 + IT_0730;
    const ccomplex_t IT_0732 = IT_0002*IT_0617;
    const ccomplex_t IT_0733 = -IT_0732;
    const ccomplex_t IT_0734 = IT_0731 + IT_0733;
    const ccomplex_t IT_0735 = (-4)*IT_0119;
    const ccomplex_t IT_0736 = (-4)*conj(IT_0119);
    const ccomplex_t IT_0737 = (-4)*IT_0122;
    const ccomplex_t IT_0738 = (-4)*conj(IT_0122);
    const ccomplex_t IT_0739 = (-4)*IT_0133;
    const ccomplex_t IT_0740 = (-4)*IT_0108;
    const ccomplex_t IT_0741 = conj(IT_0119)*IT_0735 + IT_0126*IT_0735 +
       IT_0692*IT_0736 + IT_0694*IT_0736 + conj(IT_0122)*IT_0737 + IT_0511
      *IT_0738 + (-4)*IT_0695*(IT_0131 + (-0.25)*IT_0735 + (-0.25)*IT_0739) + (
      -4)*IT_0510*(IT_0115 + (-0.25)*IT_0735 + (-0.25)*IT_0737 + (-0.25)*IT_0740);
    const ccomplex_t IT_0742 = IT_0002*IT_0642;
    const ccomplex_t IT_0743 = -IT_0742;
    const ccomplex_t IT_0744 = IT_0731 + IT_0743;
    const ccomplex_t IT_0745 = (-4)*IT_0147;
    const ccomplex_t IT_0746 = (-4)*conj(IT_0147);
    const ccomplex_t IT_0747 = (-4)*IT_0102;
    const ccomplex_t IT_0748 = 0.666666666666667*IT_0159 + (-0.166666666666667
      )*IT_0501;
    const ccomplex_t IT_0749 = s_24 + IT_0748;
    const ccomplex_t IT_0750 = m_N_2*IT_0749;
    const ccomplex_t IT_0751 = (-6)*IT_0750;
    const ccomplex_t IT_0752 = -IT_0514;
    const ccomplex_t IT_0753 = IT_0512 + IT_0752;
    const ccomplex_t IT_0754 = IT_0002*IT_0753;
    const ccomplex_t IT_0755 = -IT_0259 + (-0.25)*IT_0754;
    const ccomplex_t IT_0756 = s_12 + IT_0755;
    const ccomplex_t IT_0757 = m_N_2*IT_0756;
    const ccomplex_t IT_0758 = (-4)*IT_0757;
    const ccomplex_t IT_0759 = IT_0536 + IT_0571;
    const ccomplex_t IT_0760 = IT_0538 + IT_0572;
    const ccomplex_t IT_0761 = IT_0569 + IT_0571;
    const ccomplex_t IT_0762 = IT_0570 + IT_0572;
    const ccomplex_t IT_0763 = (-2)*IT_0334;
    const ccomplex_t IT_0764 = (-2)*IT_0347;
    const ccomplex_t IT_0765 = (-2)*conj(IT_0334);
    const ccomplex_t IT_0766 = (-2)*conj(IT_0347);
    const ccomplex_t IT_0767 = (-2)*IT_0375;
    const ccomplex_t IT_0768 = (-2)*conj(IT_0375);
    const ccomplex_t IT_0769 = conj(IT_0375)*IT_0491 + IT_0375*IT_0494 + conj
      (IT_0103)*IT_0525 + IT_0103*IT_0528 + conj(IT_0129)*IT_0531 + IT_0129
      *IT_0534 + conj(IT_0413)*IT_0552 + IT_0413*IT_0556 + conj(IT_0429)*IT_0558
       + IT_0429*IT_0560 + IT_0535*IT_0561 + IT_0537*IT_0562 + IT_0561*IT_0631 +
       IT_0562*IT_0671 + IT_0549*IT_0759 + IT_0553*IT_0760 + conj(IT_0356)
      *IT_0761 + IT_0356*IT_0762 + IT_0265*(IT_0763 + IT_0764) + IT_0266*
      (IT_0765 + IT_0766) + conj(IT_0264)*IT_0767 + IT_0264*IT_0768;
    const ccomplex_t IT_0770 = IT_0274 + IT_0541;
    const ccomplex_t IT_0771 = IT_0002*IT_0770;
    const ccomplex_t IT_0772 = IT_0275 + IT_0544;
    const ccomplex_t IT_0773 = (-4)*IT_0772;
    const ccomplex_t IT_0774 = IT_0771 + IT_0773;
    const ccomplex_t IT_0775 = m_N_2*IT_0774;
    const ccomplex_t IT_0776 = conj(IT_0264)*IT_0479 + IT_0264*IT_0480 + conj
      (IT_0264)*IT_0481 + IT_0264*IT_0482 + IT_0497*IT_0535 + IT_0499*IT_0537 +
       IT_0491*IT_0631 + IT_0494*IT_0671 + conj(IT_0413)*IT_0759 + IT_0413
      *IT_0760 + conj(IT_0429)*IT_0761 + IT_0429*IT_0762;
    const ccomplex_t IT_0777 = (IT_0000 + -IT_0275)*(s_12 + (-0.5)*IT_0566);
    const ccomplex_t IT_0778 = 2*IT_0777;
    const ccomplex_t IT_0779 = s_13*s_24*s_34*IT_0002;
    const ccomplex_t IT_0780 = (-2)*IT_0513 + -IT_0514 + 2*IT_0779;
    const ccomplex_t IT_0781 = IT_0512 + IT_0780;
    const ccomplex_t IT_0782 = IT_0002*IT_0781;
    const ccomplex_t IT_0783 = -IT_0259 + (-0.5)*IT_0782;
    const ccomplex_t IT_0784 = s_12 + IT_0783;
    const ccomplex_t IT_0785 = m_N_2*IT_0784;
    const ccomplex_t IT_0786 = (-2)*IT_0785;
    const ccomplex_t IT_0787 = 2*IT_0129;
    const ccomplex_t IT_0788 = 2*conj(IT_0129);
    const ccomplex_t IT_0789 = conj(IT_0375)*IT_0550 + IT_0375*IT_0554 + conj
      (IT_0356)*IT_0561 + IT_0356*IT_0562 + IT_0550*IT_0581 + IT_0554*IT_0583 +
       2*conj(IT_0103)*(IT_0356 + 0.5*IT_0763 + 0.5*IT_0764) + 2*IT_0103*(conj
      (IT_0356) + 0.5*IT_0765 + 0.5*IT_0766) + conj(IT_0129)*IT_0767 + IT_0129
      *IT_0768 + IT_0549*(IT_0561 + IT_0787) + IT_0553*(IT_0562 + IT_0788);
    const ccomplex_t IT_0790 = (IT_0000 + -IT_0309)*(s_12 + (-0.5)*IT_0566);
    const ccomplex_t IT_0791 = 2*IT_0790;
    const ccomplex_t IT_0792 = IT_0000*IT_0002;
    const ccomplex_t IT_0793 = (-0.25)*IT_0318;
    const ccomplex_t IT_0794 = IT_0661 + IT_0793;
    const ccomplex_t IT_0795 = IT_0792*IT_0794;
    const ccomplex_t IT_0796 = (-4)*IT_0795;
    const ccomplex_t IT_0797 = (-2)*IT_0119;
    const ccomplex_t IT_0798 = (-2)*IT_0122;
    const ccomplex_t IT_0799 = IT_0797 + IT_0798;
    const ccomplex_t IT_0800 = (-2)*conj(IT_0119);
    const ccomplex_t IT_0801 = (-2)*conj(IT_0122);
    const ccomplex_t IT_0802 = IT_0800 + IT_0801;
    const ccomplex_t IT_0803 = (-2)*IT_0133;
    const ccomplex_t IT_0804 = (-2)*IT_0131;
    const ccomplex_t IT_0805 = IT_0797 + IT_0803 + IT_0804;
    const ccomplex_t IT_0806 = (-2)*conj(IT_0133);
    const ccomplex_t IT_0807 = (-2)*conj(IT_0131);
    const ccomplex_t IT_0808 = IT_0800 + IT_0806 + IT_0807;
    const ccomplex_t IT_0809 = (-2)*IT_0147;
    const ccomplex_t IT_0810 = (-2)*IT_0151;
    const ccomplex_t IT_0811 = IT_0809 + IT_0810;
    const ccomplex_t IT_0812 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0813 = (-2)*conj(IT_0151);
    const ccomplex_t IT_0814 = IT_0812 + IT_0813;
    const ccomplex_t IT_0815 = 9*IT_0540;
    const ccomplex_t IT_0816 = IT_0000*IT_0708;
    const ccomplex_t IT_0817 = (-4)*IT_0250;
    const ccomplex_t IT_0818 = (-4)*conj(IT_0250);
    const ccomplex_t IT_0819 = (-4)*IT_0225;
    const ccomplex_t IT_0820 = 0.25*IT_0819;
    const ccomplex_t IT_0821 = (-4)*conj(IT_0225);
    const ccomplex_t IT_0822 = 0.25*IT_0821;
    const ccomplex_t IT_0823 = (IT_0241 + IT_0375 + IT_0432 + IT_0435)*IT_0712
       + (conj(IT_0241) + conj(IT_0375) + conj(IT_0432) + conj(IT_0435))*
      (IT_0225 + 0.25*IT_0713) + 0.25*(conj(IT_0334) + conj(IT_0347) + conj
      (IT_0415) + conj(IT_0418))*IT_0714 + 0.25*(IT_0334 + IT_0347 + IT_0415 +
       IT_0418)*IT_0715 + (conj(IT_0356) + conj(IT_0413))*(IT_0241 + 0.25
      *IT_0817) + (IT_0356 + IT_0413)*(conj(IT_0241) + 0.25*IT_0818) + (conj
      (IT_0355) + conj(IT_0379) + conj(IT_0429))*(IT_0204 + IT_0820) + conj
      (IT_0250)*(IT_0204 + IT_0334 + IT_0347 + IT_0415 + IT_0418 + IT_0820) + 
      (IT_0355 + IT_0379 + IT_0429)*(conj(IT_0204) + IT_0822) + IT_0250*(conj
      (IT_0204) + conj(IT_0334) + conj(IT_0347) + conj(IT_0415) + conj(IT_0418) 
      + IT_0822);
    const ccomplex_t IT_0824 = s_14*s_23;
    const ccomplex_t IT_0825 = s_12*s_34;
    const ccomplex_t IT_0826 = IT_0824 + IT_0825;
    const ccomplex_t IT_0827 = s_13*s_24;
    const ccomplex_t IT_0828 = -IT_0827;
    const ccomplex_t IT_0829 = IT_0826 + IT_0828;
    const ccomplex_t IT_0830 = IT_0155*IT_0829;
    const ccomplex_t IT_0831 = 0.25*IT_0830;
    const ccomplex_t IT_0832 = IT_0002*IT_0512;
    const ccomplex_t IT_0833 = (-0.5)*IT_0832;
    const ccomplex_t IT_0834 = s_12 + IT_0831 + IT_0833;
    const ccomplex_t IT_0835 = 8*IT_0204;
    const ccomplex_t IT_0836 = (-8)*IT_0225;
    const ccomplex_t IT_0837 = IT_0835 + IT_0836;
    const ccomplex_t IT_0838 = 8*conj(IT_0204);
    const ccomplex_t IT_0839 = (-8)*conj(IT_0225);
    const ccomplex_t IT_0840 = IT_0838 + IT_0839;
    const ccomplex_t IT_0841 = 8*conj(IT_0225);
    const ccomplex_t IT_0842 = (-8)*conj(IT_0204);
    const ccomplex_t IT_0843 = IT_0841 + IT_0842;
    const ccomplex_t IT_0844 = 8*IT_0241;
    const ccomplex_t IT_0845 = (-8)*IT_0250;
    const ccomplex_t IT_0846 = IT_0844 + IT_0845;
    const ccomplex_t IT_0847 = 8*conj(IT_0241);
    const ccomplex_t IT_0848 = (-8)*conj(IT_0250);
    const ccomplex_t IT_0849 = IT_0847 + IT_0848;
    const ccomplex_t IT_0850 = 8*IT_0250;
    const ccomplex_t IT_0851 = (-8)*IT_0241;
    const ccomplex_t IT_0852 = IT_0850 + IT_0851;
    const ccomplex_t IT_0853 = 8*conj(IT_0250);
    const ccomplex_t IT_0854 = (-8)*conj(IT_0241);
    const ccomplex_t IT_0855 = IT_0853 + IT_0854;
    const ccomplex_t IT_0856 = 8*IT_0225;
    const ccomplex_t IT_0857 = (-0.125)*IT_0856;
    const ccomplex_t IT_0858 = -IT_0825 + -IT_0827;
    const ccomplex_t IT_0859 = IT_0824 + IT_0858;
    const ccomplex_t IT_0860 = IT_0155*IT_0859;
    const ccomplex_t IT_0861 = (-0.25)*IT_0860;
    const ccomplex_t IT_0862 = s_12 + IT_0833 + IT_0861;
    const ccomplex_t IT_0863 = 4*IT_0014*(IT_0104*IT_0124 + IT_0125*IT_0127 +
       IT_0130*IT_0135 + IT_0136*IT_0138) + 4*IT_0146*(IT_0104*IT_0149 + IT_0125
      *IT_0150 + IT_0130*IT_0152 + IT_0136*IT_0153) + IT_0167*IT_0256 + IT_0262*
      (IT_0227*IT_0265 + IT_0230*IT_0266 + IT_0252*IT_0268 + IT_0255*IT_0269) + 
      (IT_0252*IT_0265 + IT_0255*IT_0266 + IT_0227*IT_0268 + IT_0230*IT_0269)
      *IT_0279 + IT_0291*IT_0304 + (IT_0130*IT_0294 + IT_0136*IT_0297 + IT_0104
      *IT_0300 + IT_0125*IT_0303)*IT_0313 + (IT_0104*IT_0294 + IT_0125*IT_0297 +
       IT_0130*IT_0300 + IT_0136*IT_0303)*IT_0316 + IT_0324*(IT_0138*IT_0358 +
       IT_0135*IT_0363 + IT_0127*IT_0381 + IT_0124*IT_0385) + (IT_0127*IT_0358 +
       IT_0124*IT_0363 + IT_0138*IT_0381 + IT_0135*IT_0385)*IT_0395 + (IT_0150
      *IT_0358 + IT_0149*IT_0363 + IT_0153*IT_0381 + IT_0152*IT_0385)*IT_0399 +
       IT_0407*(IT_0138*IT_0420 + IT_0135*IT_0424 + IT_0127*IT_0437 + IT_0124
      *IT_0441) + (IT_0153*IT_0420 + IT_0152*IT_0424 + IT_0150*IT_0437 + IT_0149
      *IT_0441)*IT_0446 + (IT_0150*IT_0420 + IT_0149*IT_0424 + IT_0153*IT_0437 +
       IT_0152*IT_0441)*IT_0451 + IT_0463*(IT_0150*IT_0466 + IT_0149*IT_0469 +
       IT_0153*IT_0471 + IT_0152*IT_0473) + IT_0478*(IT_0126*IT_0479 + IT_0123
      *IT_0480 + IT_0137*IT_0481 + IT_0134*IT_0482) + IT_0490*(IT_0127*IT_0493 +
       IT_0124*IT_0496 + IT_0138*IT_0498 + IT_0135*IT_0500) + IT_0504*(conj
      (IT_0151)*IT_0349 + conj(IT_0147)*IT_0358 + IT_0151*IT_0361 + IT_0147
      *IT_0363 + conj(IT_0147)*IT_0381 + IT_0147*IT_0385 + IT_0376*IT_0505 +
       IT_0382*IT_0506) + IT_0509*(IT_0127*IT_0414 + IT_0124*IT_0421 + IT_0138
      *IT_0430 + IT_0135*IT_0438 + conj(IT_0131)*IT_0481 + IT_0131*IT_0482 +
       IT_0479*IT_0510 + IT_0480*IT_0511) + IT_0539 + IT_0563 + IT_0568*IT_0573 
      + 4*IT_0580*IT_0585 + 6*IT_0588*(IT_0265*IT_0269 + IT_0266*(conj(IT_0264) 
      + 0.166666666666667*IT_0589)) + 6*IT_0593*(IT_0130*(IT_0103 +
       0.166666666666667*IT_0594) + IT_0136*(conj(IT_0103) + 0.166666666666667
      *IT_0595)) + IT_0596*IT_0606 + (s_12 + IT_0610 + IT_0611 + IT_0619)
      *IT_0639 + (s_12 + IT_0610 + IT_0611 + IT_0644)*IT_0660 + (s_12 + IT_0666 
      + IT_0667 + IT_0668)*IT_0673 + IT_0679*IT_0680 + IT_0688*IT_0700 + (-6)
      *IT_0702*(IT_0130*IT_0269 + IT_0136*(conj(IT_0264) + (-0.166666666666667)
      *IT_0703) + IT_0265*(IT_0103 + (-0.166666666666667)*IT_0704) + IT_0266*
      (conj(IT_0103) + (-0.166666666666667)*IT_0705)) + 4*IT_0709*(IT_0204*conj
      (IT_0204) + IT_0241*conj(IT_0241) + IT_0225*IT_0712 + 0.25*conj(IT_0225)
      *IT_0713 + 0.25*conj(IT_0250)*IT_0714 + IT_0250*(conj(IT_0250) + 0.25
      *IT_0715)) + 4*IT_0721*IT_0728 + IT_0734*IT_0741 + IT_0744*(conj(IT_0147)
      *IT_0745 + (-4)*IT_0153*(IT_0151 + (-0.25)*IT_0745) + IT_0506*IT_0746 + (
      -4)*IT_0505*(IT_0148 + (-0.25)*IT_0745 + (-0.25)*IT_0747)) + (conj(IT_0151
      )*IT_0357 + IT_0151*IT_0362 + (IT_0377 + IT_0380)*IT_0505 + (IT_0383 +
       IT_0384)*IT_0506)*IT_0751 + IT_0758*IT_0769 + IT_0775*IT_0776 + (IT_0268
      *IT_0491 + 2*conj(IT_0264)*IT_0492 + IT_0264*IT_0494 + IT_0265*IT_0497 +
       IT_0264*IT_0499)*IT_0778 + IT_0786*IT_0789 + (IT_0104*IT_0550 + IT_0130
      *IT_0550 + conj(IT_0103)*IT_0551 + IT_0103*IT_0554 + IT_0129*IT_0554 +
       conj(IT_0129)*IT_0557)*IT_0791 + IT_0796*(IT_0153*IT_0799 + IT_0152
      *IT_0802 + IT_0505*IT_0805 + IT_0506*IT_0808 + IT_0138*IT_0809 + IT_0510
      *IT_0811 + IT_0135*IT_0812 + IT_0511*IT_0814) + (IT_0138*IT_0799 + IT_0135
      *IT_0802 + IT_0510*IT_0805 + IT_0511*IT_0808 + IT_0153*IT_0809 + IT_0505
      *IT_0811 + IT_0152*IT_0812 + IT_0506*IT_0814)*IT_0815 + 4*IT_0816*IT_0823 
      + IT_0834*(IT_0631*IT_0837 + IT_0671*IT_0840 + IT_0413*IT_0843 + IT_0535
      *IT_0846 + IT_0537*IT_0849 + conj(IT_0429)*IT_0852 + IT_0429*IT_0855 + (-8
      )*conj(IT_0413)*(IT_0204 + IT_0857)) + (IT_0581*IT_0837 + IT_0583*IT_0840 
      + IT_0356*IT_0843 + conj(IT_0375)*IT_0846 + IT_0375*IT_0849 + IT_0549
      *IT_0852 + IT_0553*IT_0855 + (-8)*conj(IT_0356)*(IT_0204 + IT_0857))
      *IT_0862;
    return create_ccomplex_return(IT_0863);
}

