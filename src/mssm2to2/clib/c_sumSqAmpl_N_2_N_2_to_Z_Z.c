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
    const ccomplex_t IT_0000 = pow(m_Z, -4);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, 2);
    const ccomplex_t IT_0003 = s_13*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = pow(m_Z, -2);
    const ccomplex_t IT_0006 = s_14*IT_0002;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = s_13*s_34;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = IT_0006 + IT_0009;
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0004 + 0.5*IT_0007 + (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = s_14 + IT_0012;
    const ccomplex_t IT_0014 = m_N_2*IT_0013;
    const ccomplex_t IT_0015 = (-2)*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0018, -1);
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0025;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + IT_0024 + 
      -IT_0026 + -IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0031 = IT_0019*IT_0030;
    const ccomplex_t IT_0032 = IT_0023*IT_0030;
    const ccomplex_t IT_0033 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0034 = IT_0019*IT_0033;
    const ccomplex_t IT_0035 = IT_0023*IT_0033;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0031 + IT_0032 + 
      -IT_0034 + -IT_0035);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0029*IT_0037;
    const ccomplex_t IT_0039 = pow(m_N_2, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0002 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0044 = IT_0019*IT_0043;
    const ccomplex_t IT_0045 = IT_0023*IT_0043;
    const ccomplex_t IT_0046 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0047 = IT_0019*IT_0046;
    const ccomplex_t IT_0048 = IT_0023*IT_0046;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0044 + IT_0045 + 
      -IT_0047 + -IT_0048);
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = cpow(IT_0050, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_13 + IT_0002 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = IT_0042 + -IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0038*IT_0057;
    const ccomplex_t IT_0059 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0060 = IT_0019*IT_0059;
    const ccomplex_t IT_0061 = IT_0023*IT_0059;
    const ccomplex_t IT_0062 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0063 = IT_0019*IT_0062;
    const ccomplex_t IT_0064 = IT_0023*IT_0062;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0060 + IT_0061 + 
      -IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0068 = IT_0019*IT_0067;
    const ccomplex_t IT_0069 = IT_0023*IT_0067;
    const ccomplex_t IT_0070 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0071 = IT_0019*IT_0070;
    const ccomplex_t IT_0072 = IT_0023*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0002 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0051*IT_0080;
    const ccomplex_t IT_0082 = -IT_0058 + -IT_0078 + IT_0081;
    const ccomplex_t IT_0083 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0084 = IT_0019*IT_0083;
    const ccomplex_t IT_0085 = IT_0023*IT_0083;
    const ccomplex_t IT_0086 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0019*IT_0086;
    const ccomplex_t IT_0088 = IT_0023*IT_0086;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0084 + IT_0085 + 
      -IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0092 = IT_0019*IT_0091;
    const ccomplex_t IT_0093 = IT_0023*IT_0091;
    const ccomplex_t IT_0094 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0095 = IT_0019*IT_0094;
    const ccomplex_t IT_0096 = IT_0023*IT_0094;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0092 + IT_0093 + 
      -IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0090*IT_0098;
    const ccomplex_t IT_0100 = cpow((-2)*s_13 + IT_0002 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0099*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_13 + IT_0002 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0075*IT_0104;
    const ccomplex_t IT_0106 = IT_0102 + IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0099*IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = conj(IT_0055) + conj(IT_0082) + conj(IT_0106) +
       conj(IT_0110);
    const ccomplex_t IT_0112 = cos(alpha);
    const ccomplex_t IT_0113 = IT_0017*IT_0112;
    const ccomplex_t IT_0114 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = sin(alpha);
    const ccomplex_t IT_0117 = IT_0017*IT_0116;
    const ccomplex_t IT_0118 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0119 = IT_0117*IT_0118;
    const ccomplex_t IT_0120 = IT_0022*IT_0112;
    const ccomplex_t IT_0121 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = IT_0022*IT_0116;
    const ccomplex_t IT_0124 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0115 + -IT_0119 + 
      -IT_0122 + IT_0125);
    const ccomplex_t IT_0127 = cos(beta);
    const ccomplex_t IT_0128 = cpow(IT_0016, -4);
    const ccomplex_t IT_0129 = IT_0018*IT_0112*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = sin(beta);
    const ccomplex_t IT_0131 = IT_0018*IT_0116*IT_0128*IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0018, -3);
    const ccomplex_t IT_0133 = IT_0112*IT_0127*IT_0132;
    const ccomplex_t IT_0134 = IT_0116*IT_0130*IT_0132;
    const ccomplex_t IT_0135 = cpow(IT_0016, -2);
    const ccomplex_t IT_0136 = IT_0022*IT_0112*IT_0127*IT_0135;
    const ccomplex_t IT_0137 = IT_0022*IT_0116*IT_0130*IT_0135;
    const ccomplex_t IT_0138 = cpow(IT_0018, -2);
    const ccomplex_t IT_0139 = cpow(IT_0135 + IT_0138, -1);
    const ccomplex_t IT_0140 = (IT_0129 + IT_0131 + IT_0133 + IT_0134 + 2
      *IT_0136 + 2*IT_0137)*IT_0139;
    const ccomplex_t IT_0141 = m_W*e_em;
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_12 + (-2)*IT_0039 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0126*IT_0146;
    const ccomplex_t IT_0148 = IT_0113*IT_0118;
    const ccomplex_t IT_0149 = IT_0114*IT_0117;
    const ccomplex_t IT_0150 = IT_0120*IT_0124;
    const ccomplex_t IT_0151 = IT_0121*IT_0123;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0148 + IT_0149 + 
      -IT_0150 + -IT_0151);
    const ccomplex_t IT_0153 = -IT_0152;
    const ccomplex_t IT_0154 = IT_0018*IT_0112*IT_0128*IT_0130;
    const ccomplex_t IT_0155 = IT_0018*IT_0116*IT_0127*IT_0128;
    const ccomplex_t IT_0156 = IT_0112*IT_0130*IT_0132;
    const ccomplex_t IT_0157 = IT_0116*IT_0127*IT_0132;
    const ccomplex_t IT_0158 = IT_0022*IT_0112*IT_0130*IT_0135;
    const ccomplex_t IT_0159 = IT_0022*IT_0116*IT_0127*IT_0135;
    const ccomplex_t IT_0160 = IT_0139*(IT_0154 + -IT_0155 + IT_0156 + 
      -IT_0157 + 2*IT_0158 + (-2)*IT_0159);
    const ccomplex_t IT_0161 = IT_0141*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0039 + -reg_prop, -1);
    const ccomplex_t IT_0164 = IT_0162*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = IT_0153*IT_0165;
    const ccomplex_t IT_0167 = (-0.5)*IT_0147 + (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = 2*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0147 + 0.5*IT_0166;
    const ccomplex_t IT_0170 = (-2)*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = IT_0055 + IT_0082 + IT_0106 + IT_0110;
    const ccomplex_t IT_0173 = 2*conj(IT_0167);
    const ccomplex_t IT_0174 = (-2)*conj(IT_0169);
    const ccomplex_t IT_0175 = IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = -IT_0058;
    const ccomplex_t IT_0177 = (-0.5)*IT_0049;
    const ccomplex_t IT_0178 = cpow(IT_0177, 2);
    const ccomplex_t IT_0179 = IT_0053*IT_0178;
    const ccomplex_t IT_0180 = IT_0042 + -IT_0179;
    const ccomplex_t IT_0181 = IT_0080*IT_0178;
    const ccomplex_t IT_0182 = -IT_0078 + IT_0181;
    const ccomplex_t IT_0183 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0184 = IT_0113*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0186 = IT_0117*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0188 = IT_0120*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0190 = IT_0123*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0184 + -IT_0186 + 
      -IT_0188 + IT_0190);
    const ccomplex_t IT_0192 = IT_0146*IT_0191;
    const ccomplex_t IT_0193 = IT_0113*IT_0185;
    const ccomplex_t IT_0194 = IT_0117*IT_0183;
    const ccomplex_t IT_0195 = IT_0120*IT_0189;
    const ccomplex_t IT_0196 = IT_0123*IT_0187;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + IT_0194 + 
      -IT_0195 + -IT_0196);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0165*IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0192 + (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = 2*IT_0200;
    const ccomplex_t IT_0202 = 0.5*IT_0192 + 0.5*IT_0199;
    const ccomplex_t IT_0203 = (-2)*IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = 2*conj(IT_0200);
    const ccomplex_t IT_0206 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0207 = IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_0209 = IT_0039 + IT_0208;
    const ccomplex_t IT_0210 = pow(s_13, 2);
    const ccomplex_t IT_0211 = IT_0005*IT_0210;
    const ccomplex_t IT_0212 = pow(s_14, 2);
    const ccomplex_t IT_0213 = IT_0005*IT_0212;
    const ccomplex_t IT_0214 = -IT_0211 + -IT_0213;
    const ccomplex_t IT_0215 = IT_0209 + IT_0214;
    const ccomplex_t IT_0216 = m_N_2*IT_0215;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = 2*IT_0058 + 2*IT_0078 + (-2)*IT_0081;
    const ccomplex_t IT_0219 = 2*conj(IT_0169);
    const ccomplex_t IT_0220 = (-2)*conj(IT_0167);
    const ccomplex_t IT_0221 = IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0058 + 2*IT_0078 + (-2)*IT_0181;
    const ccomplex_t IT_0223 = 2*IT_0202;
    const ccomplex_t IT_0224 = (-2)*IT_0200;
    const ccomplex_t IT_0225 = IT_0223 + IT_0224;
    const ccomplex_t IT_0226 = 2*conj(IT_0202);
    const ccomplex_t IT_0227 = (-2)*conj(IT_0200);
    const ccomplex_t IT_0228 = IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0230 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0231 = s_14*s_24*IT_0005;
    const ccomplex_t IT_0232 = -IT_0230 + -IT_0231;
    const ccomplex_t IT_0233 = s_12 + IT_0229 + IT_0232;
    const ccomplex_t IT_0234 = m_N_2*IT_0233;
    const ccomplex_t IT_0235 = 2*IT_0169;
    const ccomplex_t IT_0236 = (-2)*IT_0167;
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = s_23*IT_0002;
    const ccomplex_t IT_0239 = IT_0001*IT_0238;
    const ccomplex_t IT_0240 = s_24*IT_0002;
    const ccomplex_t IT_0241 = IT_0005*IT_0240;
    const ccomplex_t IT_0242 = s_23*s_34;
    const ccomplex_t IT_0243 = (-2)*IT_0242;
    const ccomplex_t IT_0244 = IT_0240 + IT_0243;
    const ccomplex_t IT_0245 = IT_0005*IT_0244;
    const ccomplex_t IT_0246 = (-0.5)*IT_0239 + 0.5*IT_0241 + (-0.5)*IT_0245;
    const ccomplex_t IT_0247 = s_24 + IT_0246;
    const ccomplex_t IT_0248 = m_N_2*IT_0247;
    const ccomplex_t IT_0249 = (-2)*IT_0248;
    const ccomplex_t IT_0250 = s_14*s_23;
    const ccomplex_t IT_0251 = s_13*s_24;
    const ccomplex_t IT_0252 = s_12*s_34;
    const ccomplex_t IT_0253 = -IT_0251 + -IT_0252;
    const ccomplex_t IT_0254 = IT_0250 + IT_0253;
    const ccomplex_t IT_0255 = IT_0001*IT_0254;
    const ccomplex_t IT_0256 = (-0.25)*IT_0255;
    const ccomplex_t IT_0257 = s_12*IT_0002;
    const ccomplex_t IT_0258 = IT_0005*IT_0257;
    const ccomplex_t IT_0259 = (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = s_12 + IT_0256 + IT_0259;
    const ccomplex_t IT_0261 = cpow(IT_0098, 2);
    const ccomplex_t IT_0262 = m_N_1*IT_0261;
    const ccomplex_t IT_0263 = m_N_2*IT_0099;
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = IT_0262 + IT_0264;
    const ccomplex_t IT_0266 = IT_0101*IT_0265;
    const ccomplex_t IT_0267 = cpow(IT_0029, 2);
    const ccomplex_t IT_0268 = m_N_4*IT_0267;
    const ccomplex_t IT_0269 = IT_0041*IT_0268;
    const ccomplex_t IT_0270 = cpow(IT_0066, 2);
    const ccomplex_t IT_0271 = m_N_3*IT_0270;
    const ccomplex_t IT_0272 = m_N_2*IT_0075;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = IT_0104*IT_0274;
    const ccomplex_t IT_0276 = IT_0266 + IT_0269 + IT_0275;
    const ccomplex_t IT_0277 = 8*IT_0169;
    const ccomplex_t IT_0278 = (-0.125)*IT_0277;
    const ccomplex_t IT_0279 = IT_0167 + IT_0278;
    const ccomplex_t IT_0280 = IT_0050*IT_0177;
    const ccomplex_t IT_0281 = m_N_2*IT_0280;
    const ccomplex_t IT_0282 = m_N_2*IT_0051;
    const ccomplex_t IT_0283 = -IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = IT_0053*IT_0284;
    const ccomplex_t IT_0286 = m_N_2*IT_0038;
    const ccomplex_t IT_0287 = IT_0041*IT_0286;
    const ccomplex_t IT_0288 = IT_0285 + IT_0287;
    const ccomplex_t IT_0289 = 8*IT_0167;
    const ccomplex_t IT_0290 = (-8)*IT_0169;
    const ccomplex_t IT_0291 = IT_0289 + IT_0290;
    const ccomplex_t IT_0292 = 8*conj(IT_0167);
    const ccomplex_t IT_0293 = (-8)*conj(IT_0169);
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = 8*conj(IT_0169);
    const ccomplex_t IT_0296 = (-8)*conj(IT_0167);
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = m_N_2*IT_0178;
    const ccomplex_t IT_0299 = -IT_0053*(IT_0281 + -IT_0298);
    const ccomplex_t IT_0300 = -IT_0299;
    const ccomplex_t IT_0301 = conj(IT_0287) + conj(IT_0300);
    const ccomplex_t IT_0302 = 8*IT_0200;
    const ccomplex_t IT_0303 = (-8)*IT_0202;
    const ccomplex_t IT_0304 = IT_0302 + IT_0303;
    const ccomplex_t IT_0305 = IT_0287 + IT_0300;
    const ccomplex_t IT_0306 = 8*conj(IT_0200);
    const ccomplex_t IT_0307 = (-8)*conj(IT_0202);
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = cpow(IT_0090, 2);
    const ccomplex_t IT_0310 = m_N_1*IT_0309;
    const ccomplex_t IT_0311 = IT_0264 + IT_0310;
    const ccomplex_t IT_0312 = IT_0101*IT_0311;
    const ccomplex_t IT_0313 = cpow(IT_0037, 2);
    const ccomplex_t IT_0314 = m_N_4*IT_0313;
    const ccomplex_t IT_0315 = IT_0041*IT_0314;
    const ccomplex_t IT_0316 = cpow(IT_0074, 2);
    const ccomplex_t IT_0317 = m_N_3*IT_0316;
    const ccomplex_t IT_0318 = IT_0273 + IT_0317;
    const ccomplex_t IT_0319 = IT_0104*IT_0318;
    const ccomplex_t IT_0320 = IT_0312 + IT_0315 + IT_0319;
    const ccomplex_t IT_0321 = 8*IT_0202;
    const ccomplex_t IT_0322 = (-8)*IT_0200;
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = 8*conj(IT_0202);
    const ccomplex_t IT_0325 = (-8)*conj(IT_0200);
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = conj(IT_0276)*IT_0279 + (-0.125)*conj(IT_0288)
      *IT_0291 + (-0.125)*IT_0288*IT_0294 + (-0.125)*IT_0276*IT_0297 + (-0.125)
      *IT_0301*IT_0304 + (-0.125)*IT_0305*IT_0308 + (-0.125)*conj(IT_0320)
      *IT_0323 + (-0.125)*IT_0320*IT_0326;
    const ccomplex_t IT_0328 = IT_0250 + IT_0252;
    const ccomplex_t IT_0329 = -IT_0251;
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = IT_0001*IT_0330;
    const ccomplex_t IT_0332 = 0.25*IT_0331;
    const ccomplex_t IT_0333 = s_12 + IT_0259 + IT_0332;
    const ccomplex_t IT_0334 = IT_0108*(IT_0262 + -IT_0263);
    const ccomplex_t IT_0335 = IT_0057*(IT_0268 + -IT_0286);
    const ccomplex_t IT_0336 = -IT_0335;
    const ccomplex_t IT_0337 = -IT_0298;
    const ccomplex_t IT_0338 = IT_0281 + IT_0337;
    const ccomplex_t IT_0339 = IT_0080*IT_0338;
    const ccomplex_t IT_0340 = IT_0077*IT_0274;
    const ccomplex_t IT_0341 = IT_0339 + -IT_0340;
    const ccomplex_t IT_0342 = conj(IT_0336) + conj(IT_0341);
    const ccomplex_t IT_0343 = IT_0336 + IT_0341;
    const ccomplex_t IT_0344 = -IT_0080*(IT_0281 + -IT_0282);
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = -IT_0057*(IT_0286 + -IT_0314);
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = -IT_0077*(IT_0272 + -IT_0317);
    const ccomplex_t IT_0349 = -IT_0348;
    const ccomplex_t IT_0350 = conj(IT_0345) + conj(IT_0347) + conj(IT_0349);
    const ccomplex_t IT_0351 = IT_0345 + IT_0347 + IT_0349;
    const ccomplex_t IT_0352 = -IT_0108*(IT_0263 + -IT_0310);
    const ccomplex_t IT_0353 = (-0.125)*IT_0297*IT_0334 + IT_0279*conj(IT_0334
      ) + (-0.125)*IT_0291*IT_0342 + (-0.125)*IT_0294*IT_0343 + (-0.125)*IT_0304
      *IT_0350 + (-0.125)*IT_0308*IT_0351 + (-0.125)*IT_0326*IT_0352 + (-0.125)
      *IT_0323*conj(IT_0352);
    const ccomplex_t IT_0354 = pow(m_Z, 4);
    const ccomplex_t IT_0355 = s_14*IT_0354;
    const ccomplex_t IT_0356 = IT_0000*IT_0355;
    const ccomplex_t IT_0357 = (-10)*s_14;
    const ccomplex_t IT_0358 = IT_0356 + IT_0357;
    const ccomplex_t IT_0359 = m_N_2*IT_0358;
    const ccomplex_t IT_0360 = 2*IT_0336;
    const ccomplex_t IT_0361 = 2*IT_0341;
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = (-2)*IT_0334;
    const ccomplex_t IT_0364 = IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = 2*conj(IT_0336);
    const ccomplex_t IT_0366 = 2*conj(IT_0341);
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = (-2)*conj(IT_0334);
    const ccomplex_t IT_0369 = IT_0367 + IT_0368;
    const ccomplex_t IT_0370 = conj(IT_0176) + conj(IT_0182);
    const ccomplex_t IT_0371 = 2*IT_0345;
    const ccomplex_t IT_0372 = 2*IT_0347;
    const ccomplex_t IT_0373 = 2*IT_0349;
    const ccomplex_t IT_0374 = (-2)*IT_0352;
    const ccomplex_t IT_0375 = IT_0371 + IT_0372 + IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = IT_0176 + IT_0182;
    const ccomplex_t IT_0377 = 2*conj(IT_0345);
    const ccomplex_t IT_0378 = 2*conj(IT_0347);
    const ccomplex_t IT_0379 = 2*conj(IT_0349);
    const ccomplex_t IT_0380 = (-2)*conj(IT_0352);
    const ccomplex_t IT_0381 = IT_0377 + IT_0378 + IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = conj(IT_0110)*IT_0360;
    const ccomplex_t IT_0383 = IT_0110*IT_0365;
    const ccomplex_t IT_0384 = conj(IT_0110)*IT_0371;
    const ccomplex_t IT_0385 = IT_0110*IT_0377;
    const ccomplex_t IT_0386 = conj(IT_0110)*IT_0361;
    const ccomplex_t IT_0387 = IT_0110*IT_0366;
    const ccomplex_t IT_0388 = conj(IT_0110)*IT_0372;
    const ccomplex_t IT_0389 = IT_0110*IT_0378;
    const ccomplex_t IT_0390 = conj(IT_0110)*IT_0373;
    const ccomplex_t IT_0391 = IT_0110*IT_0379;
    const ccomplex_t IT_0392 = s_13*s_23;
    const ccomplex_t IT_0393 = (-2)*IT_0392;
    const ccomplex_t IT_0394 = IT_0257 + IT_0393;
    const ccomplex_t IT_0395 = IT_0005*IT_0394;
    const ccomplex_t IT_0396 = (s_12 + -IT_0231)*(s_12 + (-0.5)*IT_0395);
    const ccomplex_t IT_0397 = (-2)*IT_0396;
    const ccomplex_t IT_0398 = (-2)*IT_0042 + 2*IT_0179;
    const ccomplex_t IT_0399 = 2*IT_0222;
    const ccomplex_t IT_0400 = 2*conj(IT_0222);
    const ccomplex_t IT_0401 = (-2)*IT_0042 + 2*IT_0054;
    const ccomplex_t IT_0402 = 2*conj(IT_0218);
    const ccomplex_t IT_0403 = s_14*IT_0039;
    const ccomplex_t IT_0404 = pow(s_23, 2);
    const ccomplex_t IT_0405 = s_14*IT_0404;
    const ccomplex_t IT_0406 = s_12*s_24*IT_0002;
    const ccomplex_t IT_0407 = s_12*s_23*s_34;
    const ccomplex_t IT_0408 = s_14*IT_0002*IT_0039;
    const ccomplex_t IT_0409 = 0.5*IT_0406 + -IT_0407 + (-0.5)*IT_0408;
    const ccomplex_t IT_0410 = IT_0405 + IT_0409;
    const ccomplex_t IT_0411 = IT_0005*IT_0410;
    const ccomplex_t IT_0412 = s_12*s_24;
    const ccomplex_t IT_0413 = -IT_0412;
    const ccomplex_t IT_0414 = IT_0403 + IT_0411 + IT_0413;
    const ccomplex_t IT_0415 = 4*IT_0398;
    const ccomplex_t IT_0416 = 4*conj(IT_0398);
    const ccomplex_t IT_0417 = 4*IT_0401;
    const ccomplex_t IT_0418 = 4*conj(IT_0401);
    const ccomplex_t IT_0419 = conj(IT_0082) + conj(IT_0176) + conj(IT_0182);
    const ccomplex_t IT_0420 = (-2)*IT_0102 + (-2)*IT_0105;
    const ccomplex_t IT_0421 = 4*IT_0420;
    const ccomplex_t IT_0422 = 4*conj(IT_0420);
    const ccomplex_t IT_0423 = 1.125*IT_0039;
    const ccomplex_t IT_0424 = (-16)*IT_0288;
    const ccomplex_t IT_0425 = (-16)*conj(IT_0288);
    const ccomplex_t IT_0426 = (-16)*IT_0336;
    const ccomplex_t IT_0427 = (-16)*IT_0341;
    const ccomplex_t IT_0428 = (-16)*conj(IT_0336);
    const ccomplex_t IT_0429 = (-16)*conj(IT_0341);
    const ccomplex_t IT_0430 = 16*conj(IT_0288)*IT_0320 + conj(IT_0276)*(16
      *IT_0287 + 16*IT_0300 + (-16)*IT_0320) + 16*IT_0288*conj(IT_0320) +
       IT_0276*(16*conj(IT_0287) + 16*conj(IT_0300) + (-16)*conj(IT_0320)) + (16
      *conj(IT_0336) + 16*conj(IT_0341))*IT_0352 + conj(IT_0334)*(16*IT_0345 +
       16*IT_0347 + 16*IT_0349 + (-16)*IT_0352) + (16*IT_0336 + 16*IT_0341)*conj
      (IT_0352) + IT_0334*(16*conj(IT_0345) + 16*conj(IT_0347) + 16*conj(IT_0349
      ) + (-16)*conj(IT_0352)) + IT_0301*IT_0424 + IT_0305*IT_0425 + IT_0350*
      (IT_0426 + IT_0427) + IT_0351*(IT_0428 + IT_0429);
    const ccomplex_t IT_0431 = pow(s_34, 2);
    const ccomplex_t IT_0432 = (-2)*IT_0431;
    const ccomplex_t IT_0433 = IT_0354 + IT_0432;
    const ccomplex_t IT_0434 = IT_0000*IT_0433;
    const ccomplex_t IT_0435 = 0.25*IT_0434;
    const ccomplex_t IT_0436 = 1 + IT_0435;
    const ccomplex_t IT_0437 = IT_0039*IT_0436;
    const ccomplex_t IT_0438 = 8*IT_0276;
    const ccomplex_t IT_0439 = (-8)*IT_0288;
    const ccomplex_t IT_0440 = IT_0438 + IT_0439;
    const ccomplex_t IT_0441 = 8*conj(IT_0276);
    const ccomplex_t IT_0442 = (-8)*conj(IT_0288);
    const ccomplex_t IT_0443 = IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = 8*IT_0334;
    const ccomplex_t IT_0445 = (-8)*IT_0336;
    const ccomplex_t IT_0446 = (-8)*IT_0341;
    const ccomplex_t IT_0447 = IT_0444 + IT_0445 + IT_0446;
    const ccomplex_t IT_0448 = (-8)*conj(IT_0336);
    const ccomplex_t IT_0449 = (-8)*conj(IT_0341);
    const ccomplex_t IT_0450 = IT_0448 + IT_0449;
    const ccomplex_t IT_0451 = (-8)*conj(IT_0287);
    const ccomplex_t IT_0452 = (-8)*conj(IT_0300);
    const ccomplex_t IT_0453 = IT_0451 + IT_0452;
    const ccomplex_t IT_0454 = (-8)*IT_0345;
    const ccomplex_t IT_0455 = (-8)*IT_0347;
    const ccomplex_t IT_0456 = (-8)*IT_0349;
    const ccomplex_t IT_0457 = IT_0454 + IT_0455 + IT_0456;
    const ccomplex_t IT_0458 = (-8)*conj(IT_0345);
    const ccomplex_t IT_0459 = (-8)*conj(IT_0347);
    const ccomplex_t IT_0460 = (-8)*conj(IT_0349);
    const ccomplex_t IT_0461 = IT_0458 + IT_0459 + IT_0460;
    const ccomplex_t IT_0462 = 8*IT_0288;
    const ccomplex_t IT_0463 = 8*conj(IT_0288);
    const ccomplex_t IT_0464 = 8*IT_0336;
    const ccomplex_t IT_0465 = 8*IT_0341;
    const ccomplex_t IT_0466 = 8*conj(IT_0336);
    const ccomplex_t IT_0467 = 8*conj(IT_0341);
    const ccomplex_t IT_0468 = (-8)*IT_0287;
    const ccomplex_t IT_0469 = (-8)*IT_0300;
    const ccomplex_t IT_0470 = conj(IT_0352)*IT_0440 + IT_0352*IT_0443 + conj
      (IT_0320)*IT_0447 + IT_0320*IT_0450 + IT_0334*IT_0453 + conj(IT_0276)
      *IT_0457 + IT_0276*IT_0461 + IT_0350*IT_0462 + IT_0351*IT_0463 + IT_0301*
      (IT_0464 + IT_0465) + IT_0305*(IT_0466 + IT_0467) + 8*conj(IT_0334)*
      (IT_0320 + 0.125*IT_0468 + 0.125*IT_0469);
    const ccomplex_t IT_0471 = pow(m_Z, 6);
    const ccomplex_t IT_0472 = s_12*IT_0471;
    const ccomplex_t IT_0473 = s_13*s_23*IT_0354;
    const ccomplex_t IT_0474 = (-2)*IT_0473;
    const ccomplex_t IT_0475 = IT_0472 + IT_0474;
    const ccomplex_t IT_0476 = IT_0000*IT_0475;
    const ccomplex_t IT_0477 = 0.25*IT_0476;
    const ccomplex_t IT_0478 = IT_0257 + IT_0477;
    const ccomplex_t IT_0479 = (-8)*IT_0180;
    const ccomplex_t IT_0480 = (-8)*conj(IT_0180);
    const ccomplex_t IT_0481 = s_14*s_24;
    const ccomplex_t IT_0482 = (-2)*IT_0481;
    const ccomplex_t IT_0483 = IT_0257 + IT_0482;
    const ccomplex_t IT_0484 = IT_0005*IT_0483;
    const ccomplex_t IT_0485 = (-0.5)*IT_0484;
    const ccomplex_t IT_0486 = (-0.5)*IT_0395;
    const ccomplex_t IT_0487 = s_12*IT_0354;
    const ccomplex_t IT_0488 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0489 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0490 = s_13*s_24*s_34;
    const ccomplex_t IT_0491 = (-2)*IT_0488 + (-2)*IT_0489 + 4*IT_0490;
    const ccomplex_t IT_0492 = IT_0487 + IT_0491;
    const ccomplex_t IT_0493 = IT_0000*IT_0492;
    const ccomplex_t IT_0494 = 0.25*IT_0493;
    const ccomplex_t IT_0495 = 8*IT_0320;
    const ccomplex_t IT_0496 = IT_0468 + IT_0469 + IT_0495;
    const ccomplex_t IT_0497 = 8*IT_0287;
    const ccomplex_t IT_0498 = 8*conj(IT_0287);
    const ccomplex_t IT_0499 = 8*IT_0300;
    const ccomplex_t IT_0500 = s_14*s_23*s_34;
    const ccomplex_t IT_0501 = s_12*IT_0431;
    const ccomplex_t IT_0502 = (-2)*IT_0488 + (-2)*IT_0489 + 2*IT_0490 + 2
      *IT_0500 + (-2)*IT_0501;
    const ccomplex_t IT_0503 = IT_0487 + IT_0502;
    const ccomplex_t IT_0504 = IT_0000*IT_0503;
    const ccomplex_t IT_0505 = (-0.25)*IT_0504;
    const ccomplex_t IT_0506 = -IT_0231;
    const ccomplex_t IT_0507 = -IT_0230;
    const ccomplex_t IT_0508 = 8*IT_0347;
    const ccomplex_t IT_0509 = 8*conj(IT_0347);
    const ccomplex_t IT_0510 = 8*IT_0349;
    const ccomplex_t IT_0511 = 8*conj(IT_0349);
    const ccomplex_t IT_0512 = conj(IT_0334)*IT_0440 + IT_0334*IT_0443 + conj
      (IT_0276)*(IT_0445 + IT_0446) + IT_0276*IT_0450 + IT_0352*IT_0453 + 8*conj
      (IT_0320)*(IT_0352 + 0.125*IT_0454 + 0.125*IT_0455 + 0.125*IT_0456) +
       IT_0320*IT_0461 + conj(IT_0341)*IT_0462 + IT_0341*IT_0463 + conj(IT_0288)
      *IT_0464 + IT_0288*IT_0466 + conj(IT_0352)*IT_0496 + conj(IT_0349)*IT_0497
       + (IT_0345 + IT_0349)*IT_0498 + conj(IT_0345)*(IT_0497 + IT_0499) + conj
      (IT_0287)*IT_0508 + IT_0305*IT_0509 + 8*conj(IT_0300)*(IT_0345 + 0.125
      *IT_0508 + 0.125*IT_0510) + IT_0300*IT_0511;
    const ccomplex_t IT_0513 = (-2)*IT_0488 + (-2)*IT_0489 + 4*IT_0500;
    const ccomplex_t IT_0514 = IT_0487 + IT_0513;
    const ccomplex_t IT_0515 = IT_0000*IT_0514;
    const ccomplex_t IT_0516 = 0.25*IT_0515;
    const ccomplex_t IT_0517 = conj(IT_0334)*IT_0447 + IT_0334*IT_0450 + conj
      (IT_0352)*IT_0457 + 8*IT_0352*(conj(IT_0352) + 0.125*IT_0458 + 0.125
      *IT_0459 + 0.125*IT_0460) + IT_0342*IT_0464 + conj(IT_0341)*IT_0465 +
       IT_0341*IT_0466 + IT_0350*IT_0508 + IT_0349*IT_0509 + (conj(IT_0345) +
       conj(IT_0349))*IT_0510 + 8*IT_0345*(conj(IT_0345) + 0.125*IT_0509 + 0.125
      *IT_0511);
    const ccomplex_t IT_0518 = -IT_0213;
    const ccomplex_t IT_0519 = IT_0039 + IT_0518;
    const ccomplex_t IT_0520 = IT_0039*IT_0519;
    const ccomplex_t IT_0521 = -IT_0520;
    const ccomplex_t IT_0522 = 6*conj(IT_0222);
    const ccomplex_t IT_0523 = s_12 + IT_0506;
    const ccomplex_t IT_0524 = IT_0039*IT_0523;
    const ccomplex_t IT_0525 = (-6)*conj(IT_0398);
    const ccomplex_t IT_0526 = pow(s_24, 2);
    const ccomplex_t IT_0527 = IT_0005*IT_0526;
    const ccomplex_t IT_0528 = -IT_0527;
    const ccomplex_t IT_0529 = IT_0039 + IT_0528;
    const ccomplex_t IT_0530 = IT_0039*IT_0529;
    const ccomplex_t IT_0531 = -IT_0530;
    const ccomplex_t IT_0532 = 6*conj(IT_0398);
    const ccomplex_t IT_0533 = (-6)*conj(IT_0222);
    const ccomplex_t IT_0534 = IT_0000*IT_0431;
    const ccomplex_t IT_0535 = 0.5*IT_0534;
    const ccomplex_t IT_0536 = 1 + IT_0535;
    const ccomplex_t IT_0537 = s_12*IT_0536;
    const ccomplex_t IT_0538 = (-4)*conj(IT_0167);
    const ccomplex_t IT_0539 = 0.25*IT_0538;
    const ccomplex_t IT_0540 = conj(IT_0169) + IT_0539;
    const ccomplex_t IT_0541 = (-4)*IT_0167;
    const ccomplex_t IT_0542 = (-4)*IT_0200;
    const ccomplex_t IT_0543 = (-4)*conj(IT_0200);
    const ccomplex_t IT_0544 = IT_0002*IT_0039;
    const ccomplex_t IT_0545 = s_13*s_14*s_34*IT_0005;
    const ccomplex_t IT_0546 = 0.5*IT_0212 + (-0.5)*IT_0544 + -IT_0545;
    const ccomplex_t IT_0547 = IT_0210 + IT_0546;
    const ccomplex_t IT_0548 = IT_0005*IT_0547;
    const ccomplex_t IT_0549 = IT_0039 + IT_0518 + IT_0548;
    const ccomplex_t IT_0550 = m_N_2*IT_0549;
    const ccomplex_t IT_0551 = (-4)*conj(IT_0288);
    const ccomplex_t IT_0552 = (-4)*IT_0222;
    const ccomplex_t IT_0553 = (-4)*conj(IT_0222);
    const ccomplex_t IT_0554 = 0.5*IT_0476;
    const ccomplex_t IT_0555 = (-4)*IT_0481;
    const ccomplex_t IT_0556 = IT_0257 + IT_0554 + IT_0555;
    const ccomplex_t IT_0557 = IT_0005*IT_0514;
    const ccomplex_t IT_0558 = -IT_0557;
    const ccomplex_t IT_0559 = IT_0556 + IT_0558;
    const ccomplex_t IT_0560 = (-4)*IT_0176;
    const ccomplex_t IT_0561 = (-4)*conj(IT_0176);
    const ccomplex_t IT_0562 = (-4)*IT_0182;
    const ccomplex_t IT_0563 = (s_12 + (-0.5)*IT_0395)*(IT_0039 + -IT_0527);
    const ccomplex_t IT_0564 = 2*IT_0563;
    const ccomplex_t IT_0565 = (-2)*IT_0398;
    const ccomplex_t IT_0566 = 2*IT_0218;
    const ccomplex_t IT_0567 = (-6)*IT_0222;
    const ccomplex_t IT_0568 = (-4)*IT_0420;
    const ccomplex_t IT_0569 = (-2)*IT_0401;
    const ccomplex_t IT_0570 = 6*conj(IT_0420)*((IT_0398 + IT_0401 + 2*IT_0420
      )*IT_0531 + 0.166666666666667*IT_0397*(IT_0399 + IT_0566) + -IT_0524*
      (IT_0218 + (-0.166666666666667)*IT_0567) + 0.166666666666667*IT_0564*
      (IT_0565 + IT_0568 + IT_0569));
    const ccomplex_t IT_0571 = s_13*s_24*s_34*IT_0005;
    const ccomplex_t IT_0572 = (-2)*IT_0392 + -IT_0481 + 2*IT_0571;
    const ccomplex_t IT_0573 = IT_0257 + IT_0572;
    const ccomplex_t IT_0574 = IT_0005*IT_0573;
    const ccomplex_t IT_0575 = -IT_0231 + (-0.5)*IT_0574;
    const ccomplex_t IT_0576 = s_12 + IT_0575;
    const ccomplex_t IT_0577 = m_N_2*IT_0576;
    const ccomplex_t IT_0578 = (-2)*IT_0577;
    const ccomplex_t IT_0579 = 2*IT_0287;
    const ccomplex_t IT_0580 = 2*IT_0300;
    const ccomplex_t IT_0581 = 2*IT_0288;
    const ccomplex_t IT_0582 = IT_0579 + IT_0580 + IT_0581;
    const ccomplex_t IT_0583 = 2*conj(IT_0287);
    const ccomplex_t IT_0584 = 2*conj(IT_0300);
    const ccomplex_t IT_0585 = 2*conj(IT_0288);
    const ccomplex_t IT_0586 = IT_0583 + IT_0584 + IT_0585;
    const ccomplex_t IT_0587 = 2*IT_0398;
    const ccomplex_t IT_0588 = 2*conj(IT_0398);
    const ccomplex_t IT_0589 = (-2)*IT_0420;
    const ccomplex_t IT_0590 = (-2)*conj(IT_0398);
    const ccomplex_t IT_0591 = (-2)*conj(IT_0420);
    const ccomplex_t IT_0592 = (-2)*conj(IT_0401);
    const ccomplex_t IT_0593 = conj(IT_0401)*IT_0581 + conj(IT_0420)*IT_0582 +
       IT_0401*IT_0585 + IT_0420*IT_0586 + IT_0301*IT_0587 + IT_0305*IT_0588 +
       conj(IT_0320)*(IT_0565 + IT_0589) + conj(IT_0276)*(IT_0569 + IT_0589) +
       IT_0320*(IT_0590 + IT_0591) + IT_0276*(IT_0591 + IT_0592);
    const ccomplex_t IT_0594 = m_N_2*IT_0005;
    const ccomplex_t IT_0595 = IT_0005*IT_0355;
    const ccomplex_t IT_0596 = (-2)*IT_0008 + (-0.5)*IT_0595;
    const ccomplex_t IT_0597 = IT_0006 + IT_0596;
    const ccomplex_t IT_0598 = IT_0594*IT_0597;
    const ccomplex_t IT_0599 = (-2)*IT_0598;
    const ccomplex_t IT_0600 = conj(IT_0106)*IT_0360;
    const ccomplex_t IT_0601 = IT_0106*IT_0365;
    const ccomplex_t IT_0602 = conj(IT_0106)*IT_0371;
    const ccomplex_t IT_0603 = IT_0106*IT_0377;
    const ccomplex_t IT_0604 = conj(IT_0106)*IT_0361;
    const ccomplex_t IT_0605 = IT_0106*IT_0366;
    const ccomplex_t IT_0606 = conj(IT_0106)*IT_0372;
    const ccomplex_t IT_0607 = IT_0106*IT_0378;
    const ccomplex_t IT_0608 = conj(IT_0106)*IT_0373;
    const ccomplex_t IT_0609 = IT_0106*IT_0379;
    const ccomplex_t IT_0610 = conj(IT_0106)*IT_0363;
    const ccomplex_t IT_0611 = IT_0106*IT_0368;
    const ccomplex_t IT_0612 = conj(IT_0106)*IT_0374;
    const ccomplex_t IT_0613 = IT_0106*IT_0380;
    const ccomplex_t IT_0614 = conj(IT_0055)*IT_0364 + IT_0055*IT_0369 + conj
      (IT_0180)*IT_0375 + IT_0180*IT_0381 + IT_0600 + IT_0601 + IT_0602 +
       IT_0603 + IT_0604 + IT_0605 + IT_0606 + IT_0607 + IT_0608 + IT_0609 +
       IT_0610 + IT_0611 + IT_0612 + IT_0613;
    const ccomplex_t IT_0615 = s_23*s_34*IT_0005;
    const ccomplex_t IT_0616 = s_24*IT_0354;
    const ccomplex_t IT_0617 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0618 = (-2)*IT_0617;
    const ccomplex_t IT_0619 = IT_0616 + IT_0618;
    const ccomplex_t IT_0620 = IT_0000*IT_0619;
    const ccomplex_t IT_0621 = (-0.666666666666667)*IT_0615 + (
      -0.166666666666667)*IT_0620;
    const ccomplex_t IT_0622 = s_24 + IT_0621;
    const ccomplex_t IT_0623 = m_N_2*IT_0622;
    const ccomplex_t IT_0624 = (-6)*IT_0623;
    const ccomplex_t IT_0625 = conj(IT_0180)*IT_0364 + IT_0180*IT_0369 + conj
      (IT_0055)*IT_0375 + IT_0055*IT_0381 + IT_0600 + IT_0601 + IT_0602 +
       IT_0603 + IT_0604 + IT_0605 + IT_0606 + IT_0607 + IT_0608 + IT_0609 +
       IT_0610 + IT_0611 + IT_0612 + IT_0613;
    const ccomplex_t IT_0626 = 0.666666666666667*IT_0007 + (-0.166666666666667
      )*IT_0356;
    const ccomplex_t IT_0627 = s_14 + IT_0626;
    const ccomplex_t IT_0628 = m_N_2*IT_0627;
    const ccomplex_t IT_0629 = (-6)*IT_0628;
    const ccomplex_t IT_0630 = conj(IT_0110)*IT_0363;
    const ccomplex_t IT_0631 = IT_0110*IT_0368;
    const ccomplex_t IT_0632 = conj(IT_0110)*IT_0374;
    const ccomplex_t IT_0633 = IT_0110*IT_0380;
    const ccomplex_t IT_0634 = (-0.666666666666667)*IT_0245 +
       0.166666666666667*IT_0620;
    const ccomplex_t IT_0635 = s_24 + IT_0634;
    const ccomplex_t IT_0636 = m_N_2*IT_0635;
    const ccomplex_t IT_0637 = 6*IT_0636;
    const ccomplex_t IT_0638 = IT_0364*IT_0370 + conj(IT_0082)*IT_0375 +
       IT_0369*IT_0376 + IT_0082*IT_0381 + IT_0382 + IT_0383 + IT_0384 + IT_0385
       + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 + IT_0391 + IT_0630 +
       IT_0631 + IT_0632 + IT_0633;
    const ccomplex_t IT_0639 = 9*IT_0544;
    const ccomplex_t IT_0640 = (-2)*IT_0176;
    const ccomplex_t IT_0641 = (-2)*IT_0182;
    const ccomplex_t IT_0642 = IT_0640 + IT_0641;
    const ccomplex_t IT_0643 = (-8)*IT_0055;
    const ccomplex_t IT_0644 = (-4)*IT_0082;
    const ccomplex_t IT_0645 = IT_0005*IT_0039;
    const ccomplex_t IT_0646 = (-0.25)*IT_0354;
    const ccomplex_t IT_0647 = IT_0431 + IT_0646;
    const ccomplex_t IT_0648 = IT_0645*IT_0647;
    const ccomplex_t IT_0649 = (-4)*IT_0648;
    const ccomplex_t IT_0650 = (-2)*IT_0180;
    const ccomplex_t IT_0651 = IT_0000*IT_0616;
    const ccomplex_t IT_0652 = 0.666666666666667*IT_0241 + (-0.166666666666667
      )*IT_0651;
    const ccomplex_t IT_0653 = s_24 + IT_0652;
    const ccomplex_t IT_0654 = m_N_2*IT_0653;
    const ccomplex_t IT_0655 = (-6)*IT_0654;
    const ccomplex_t IT_0656 = (-2)*IT_0287;
    const ccomplex_t IT_0657 = (-2)*IT_0300;
    const ccomplex_t IT_0658 = IT_0656 + IT_0657;
    const ccomplex_t IT_0659 = s_13*s_23*s_24;
    const ccomplex_t IT_0660 = s_13*s_34*IT_0039;
    const ccomplex_t IT_0661 = -IT_0408 + (-2)*IT_0659 + 2*IT_0660;
    const ccomplex_t IT_0662 = IT_0406 + IT_0661;
    const ccomplex_t IT_0663 = IT_0005*IT_0662;
    const ccomplex_t IT_0664 = (-0.5)*IT_0663;
    const ccomplex_t IT_0665 = -IT_0403;
    const ccomplex_t IT_0666 = IT_0412 + IT_0664 + IT_0665;
    const ccomplex_t IT_0667 = IT_0005*IT_0492;
    const ccomplex_t IT_0668 = -IT_0667;
    const ccomplex_t IT_0669 = IT_0556 + IT_0668;
    const ccomplex_t IT_0670 = (-4)*IT_0055;
    const ccomplex_t IT_0671 = s_12*s_14;
    const ccomplex_t IT_0672 = s_12*s_14*IT_0002;
    const ccomplex_t IT_0673 = s_24*IT_0002*IT_0039;
    const ccomplex_t IT_0674 = s_24*IT_0210;
    const ccomplex_t IT_0675 = s_12*s_13*s_34;
    const ccomplex_t IT_0676 = -IT_0673 + 2*IT_0674 + (-2)*IT_0675;
    const ccomplex_t IT_0677 = IT_0672 + IT_0676;
    const ccomplex_t IT_0678 = IT_0005*IT_0677;
    const ccomplex_t IT_0679 = (-0.5)*IT_0678;
    const ccomplex_t IT_0680 = s_24*IT_0039;
    const ccomplex_t IT_0681 = -IT_0680;
    const ccomplex_t IT_0682 = IT_0671 + IT_0679 + IT_0681;
    const ccomplex_t IT_0683 = 2*IT_0109;
    const ccomplex_t IT_0684 = (-4)*IT_0683;
    const ccomplex_t IT_0685 = (-4)*IT_0218;
    const ccomplex_t IT_0686 = (-2)*conj(IT_0176);
    const ccomplex_t IT_0687 = (-2)*conj(IT_0182);
    const ccomplex_t IT_0688 = IT_0686 + IT_0687;
    const ccomplex_t IT_0689 = (-8)*conj(IT_0055);
    const ccomplex_t IT_0690 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0691 = (-2)*conj(IT_0287);
    const ccomplex_t IT_0692 = (-2)*conj(IT_0300);
    const ccomplex_t IT_0693 = IT_0691 + IT_0692;
    const ccomplex_t IT_0694 = (-4)*conj(IT_0683);
    const ccomplex_t IT_0695 = (-4)*conj(IT_0218);
    const ccomplex_t IT_0696 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0697 = s_12 + IT_0232 + IT_0696;
    const ccomplex_t IT_0698 = m_N_2*IT_0697;
    const ccomplex_t IT_0699 = IT_0168*conj(IT_0420);
    const ccomplex_t IT_0700 = IT_0173*IT_0420;
    const ccomplex_t IT_0701 = IT_0201*conj(IT_0420);
    const ccomplex_t IT_0702 = IT_0205*IT_0420;
    const ccomplex_t IT_0703 = IT_0170*conj(IT_0420);
    const ccomplex_t IT_0704 = IT_0174*IT_0420;
    const ccomplex_t IT_0705 = IT_0203*conj(IT_0420);
    const ccomplex_t IT_0706 = IT_0206*IT_0420;
    const ccomplex_t IT_0707 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_0708 = IT_0039 + IT_0707;
    const ccomplex_t IT_0709 = IT_0005*IT_0404;
    const ccomplex_t IT_0710 = -IT_0527 + -IT_0709;
    const ccomplex_t IT_0711 = IT_0708 + IT_0710;
    const ccomplex_t IT_0712 = m_N_2*IT_0711;
    const ccomplex_t IT_0713 = -IT_0712;
    const ccomplex_t IT_0714 = (-8)*IT_0176;
    const ccomplex_t IT_0715 = (-4)*IT_0180;
    const ccomplex_t IT_0716 = (-2)*IT_0288;
    const ccomplex_t IT_0717 = (-8)*conj(IT_0176);
    const ccomplex_t IT_0718 = (-2)*conj(IT_0288);
    const ccomplex_t IT_0719 = -IT_0544;
    const ccomplex_t IT_0720 = IT_0526 + IT_0719;
    const ccomplex_t IT_0721 = IT_0005*IT_0720;
    const ccomplex_t IT_0722 = -IT_0039;
    const ccomplex_t IT_0723 = IT_0527 + IT_0722;
    const ccomplex_t IT_0724 = (-4)*IT_0723;
    const ccomplex_t IT_0725 = IT_0721 + IT_0724;
    const ccomplex_t IT_0726 = m_N_2*IT_0725;
    const ccomplex_t IT_0727 = 2*IT_0420;
    const ccomplex_t IT_0728 = IT_0587 + IT_0727;
    const ccomplex_t IT_0729 = 2*conj(IT_0420);
    const ccomplex_t IT_0730 = IT_0588 + IT_0729;
    const ccomplex_t IT_0731 = 2*IT_0401;
    const ccomplex_t IT_0732 = IT_0727 + IT_0731;
    const ccomplex_t IT_0733 = 2*conj(IT_0401);
    const ccomplex_t IT_0734 = IT_0729 + IT_0733;
    const ccomplex_t IT_0735 = conj(IT_0288)*IT_0565 + IT_0301*IT_0569 +
       IT_0288*IT_0590 + IT_0305*IT_0592 + conj(IT_0420)*(IT_0658 + IT_0716) +
       IT_0420*(IT_0693 + IT_0718) + conj(IT_0276)*IT_0728 + IT_0276*IT_0730 +
       conj(IT_0320)*IT_0732 + IT_0320*IT_0734;
    const ccomplex_t IT_0736 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0737 = (-2)*IT_0736;
    const ccomplex_t IT_0738 = IT_0355 + IT_0737;
    const ccomplex_t IT_0739 = IT_0000*IT_0738;
    const ccomplex_t IT_0740 = (-0.666666666666667)*IT_0011 +
       0.166666666666667*IT_0739;
    const ccomplex_t IT_0741 = s_14 + IT_0740;
    const ccomplex_t IT_0742 = m_N_2*IT_0741;
    const ccomplex_t IT_0743 = 6*IT_0742;
    const ccomplex_t IT_0744 = 2*IT_0276;
    const ccomplex_t IT_0745 = IT_0716 + IT_0744;
    const ccomplex_t IT_0746 = 2*conj(IT_0276);
    const ccomplex_t IT_0747 = IT_0718 + IT_0746;
    const ccomplex_t IT_0748 = 2*IT_0320;
    const ccomplex_t IT_0749 = IT_0656 + IT_0657 + IT_0748;
    const ccomplex_t IT_0750 = 2*conj(IT_0320);
    const ccomplex_t IT_0751 = IT_0691 + IT_0692 + IT_0750;
    const ccomplex_t IT_0752 = conj(IT_0106)*IT_0744;
    const ccomplex_t IT_0753 = IT_0106*IT_0746;
    const ccomplex_t IT_0754 = conj(IT_0106)*IT_0748;
    const ccomplex_t IT_0755 = IT_0106*IT_0750;
    const ccomplex_t IT_0756 = conj(IT_0106)*IT_0656;
    const ccomplex_t IT_0757 = IT_0106*IT_0691;
    const ccomplex_t IT_0758 = conj(IT_0106)*IT_0657;
    const ccomplex_t IT_0759 = IT_0106*IT_0692;
    const ccomplex_t IT_0760 = conj(IT_0106)*IT_0716;
    const ccomplex_t IT_0761 = IT_0106*IT_0718;
    const ccomplex_t IT_0762 = (-10)*s_24;
    const ccomplex_t IT_0763 = IT_0651 + IT_0762;
    const ccomplex_t IT_0764 = m_N_2*IT_0763;
    const ccomplex_t IT_0765 = IT_0005*IT_0616;
    const ccomplex_t IT_0766 = (-2)*IT_0242 + (-0.5)*IT_0765;
    const ccomplex_t IT_0767 = IT_0240 + IT_0766;
    const ccomplex_t IT_0768 = IT_0594*IT_0767;
    const ccomplex_t IT_0769 = (-2)*IT_0768;
    const ccomplex_t IT_0770 = conj(IT_0110)*IT_0744;
    const ccomplex_t IT_0771 = IT_0110*IT_0746;
    const ccomplex_t IT_0772 = conj(IT_0110)*IT_0748;
    const ccomplex_t IT_0773 = IT_0110*IT_0750;
    const ccomplex_t IT_0774 = conj(IT_0110)*IT_0656;
    const ccomplex_t IT_0775 = IT_0110*IT_0691;
    const ccomplex_t IT_0776 = conj(IT_0110)*IT_0657;
    const ccomplex_t IT_0777 = IT_0110*IT_0692;
    const ccomplex_t IT_0778 = conj(IT_0110)*IT_0716;
    const ccomplex_t IT_0779 = IT_0110*IT_0718;
    const ccomplex_t IT_0780 = s_13*s_34*IT_0005;
    const ccomplex_t IT_0781 = (-0.166666666666667)*IT_0739 + (
      -0.666666666666667)*IT_0780;
    const ccomplex_t IT_0782 = s_14 + IT_0781;
    const ccomplex_t IT_0783 = m_N_2*IT_0782;
    const ccomplex_t IT_0784 = (-6)*IT_0783;
    const ccomplex_t IT_0785 = (IT_0039 + -IT_0213)*(s_12 + (-0.5)*IT_0395);
    const ccomplex_t IT_0786 = 2*IT_0785;
    const ccomplex_t IT_0787 = (-2)*IT_0222;
    const ccomplex_t IT_0788 = -IT_0481;
    const ccomplex_t IT_0789 = IT_0257 + IT_0788;
    const ccomplex_t IT_0790 = IT_0005*IT_0789;
    const ccomplex_t IT_0791 = -IT_0231 + (-0.25)*IT_0790;
    const ccomplex_t IT_0792 = s_12 + IT_0791;
    const ccomplex_t IT_0793 = m_N_2*IT_0792;
    const ccomplex_t IT_0794 = (-4)*IT_0793;
    const ccomplex_t IT_0795 = (-2)*IT_0683;
    const ccomplex_t IT_0796 = (-2)*IT_0218;
    const ccomplex_t IT_0797 = IT_0795 + IT_0796;
    const ccomplex_t IT_0798 = (-2)*conj(IT_0683);
    const ccomplex_t IT_0799 = (-2)*conj(IT_0218);
    const ccomplex_t IT_0800 = IT_0798 + IT_0799;
    const ccomplex_t IT_0801 = (-2)*IT_0336;
    const ccomplex_t IT_0802 = (-2)*IT_0341;
    const ccomplex_t IT_0803 = (-2)*conj(IT_0336);
    const ccomplex_t IT_0804 = (-2)*conj(IT_0341);
    const ccomplex_t IT_0805 = (-2)*IT_0345;
    const ccomplex_t IT_0806 = (-2)*IT_0347;
    const ccomplex_t IT_0807 = (-2)*IT_0349;
    const ccomplex_t IT_0808 = (-2)*conj(IT_0345);
    const ccomplex_t IT_0809 = (-2)*conj(IT_0347);
    const ccomplex_t IT_0810 = (-2)*conj(IT_0349);
    const ccomplex_t IT_0811 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0812 = conj(IT_0288)*IT_0399 + IT_0288*IT_0400 +
       IT_0305*IT_0402 + IT_0350*IT_0565 + IT_0301*IT_0566 + IT_0351*IT_0590 +
       IT_0586*IT_0683 + IT_0582*conj(IT_0683) + conj(IT_0352)*IT_0728 + IT_0352
      *IT_0730 + conj(IT_0334)*IT_0732 + IT_0334*IT_0734 + conj(IT_0276)*
      (IT_0787 + IT_0795) + conj(IT_0320)*IT_0797 + IT_0320*IT_0800 + conj
      (IT_0401)*(IT_0801 + IT_0802) + IT_0401*(IT_0803 + IT_0804) + conj(IT_0420
      )*(IT_0801 + IT_0802 + IT_0805 + IT_0806 + IT_0807) + IT_0420*(IT_0803 +
       IT_0804 + IT_0808 + IT_0809 + IT_0810) + IT_0276*(IT_0798 + IT_0811);
    const ccomplex_t IT_0813 = IT_0212 + IT_0719;
    const ccomplex_t IT_0814 = IT_0005*IT_0813;
    const ccomplex_t IT_0815 = IT_0213 + IT_0722;
    const ccomplex_t IT_0816 = (-4)*IT_0815;
    const ccomplex_t IT_0817 = IT_0814 + IT_0816;
    const ccomplex_t IT_0818 = m_N_2*IT_0817;
    const ccomplex_t IT_0819 = (-4)*IT_0288;
    const ccomplex_t IT_0820 = 6*conj(IT_0218)*(0.166666666666667*IT_0225
      *IT_0234 + 0.166666666666667*IT_0217*IT_0237 + IT_0222*IT_0521 + -(IT_0398
       + IT_0420)*IT_0524 + 0.166666666666667*IT_0786*IT_0796 +
       0.666666666666667*IT_0550*(IT_0276 + 0.25*IT_0819));
    const ccomplex_t IT_0821 = s_13*s_14*s_23;
    const ccomplex_t IT_0822 = s_23*s_34*IT_0039;
    const ccomplex_t IT_0823 = -IT_0673 + (-2)*IT_0821 + 2*IT_0822;
    const ccomplex_t IT_0824 = IT_0672 + IT_0823;
    const ccomplex_t IT_0825 = IT_0005*IT_0824;
    const ccomplex_t IT_0826 = 2*IT_0680;
    const ccomplex_t IT_0827 = (-2)*IT_0671;
    const ccomplex_t IT_0828 = IT_0825 + IT_0826 + IT_0827;
    const ccomplex_t IT_0829 = s_14*s_23*s_34*IT_0005;
    const ccomplex_t IT_0830 = (-2)*IT_0392 + -IT_0481 + 2*IT_0829;
    const ccomplex_t IT_0831 = IT_0257 + IT_0830;
    const ccomplex_t IT_0832 = IT_0005*IT_0831;
    const ccomplex_t IT_0833 = -IT_0231 + (-0.5)*IT_0832;
    const ccomplex_t IT_0834 = s_12 + IT_0833;
    const ccomplex_t IT_0835 = m_N_2*IT_0834;
    const ccomplex_t IT_0836 = (-2)*IT_0835;
    const ccomplex_t IT_0837 = 2*IT_0683;
    const ccomplex_t IT_0838 = 2*conj(IT_0683);
    const ccomplex_t IT_0839 = IT_0836*(IT_0342*IT_0787 + IT_0350*IT_0796 +
       IT_0351*IT_0799 + IT_0343*IT_0811 + conj(IT_0334)*(IT_0399 + IT_0837) +
       conj(IT_0352)*(IT_0566 + IT_0837) + IT_0334*(IT_0400 + IT_0838) + IT_0352
      *(IT_0402 + IT_0838));
    const ccomplex_t IT_0840 = 0.125*IT_0417;
    const ccomplex_t IT_0841 = 0.125*IT_0415 + IT_0420 + IT_0840;
    const ccomplex_t IT_0842 = (-4)*IT_0106;
    const ccomplex_t IT_0843 = (-2)*IT_0055;
    const ccomplex_t IT_0844 = IT_0650 + IT_0842 + IT_0843;
    const ccomplex_t IT_0845 = (-2)*IT_0082;
    const ccomplex_t IT_0846 = 0.125*IT_0682*(IT_0552 + IT_0685) + -IT_0478*
      (IT_0082 + IT_0182 + (-0.125)*IT_0714) + -IT_0669*(IT_0106 + (-0.125)
      *IT_0670 + (-0.125)*IT_0715) + IT_0666*IT_0841 + 0.125*IT_0639*IT_0844 +
       0.125*IT_0649*(IT_0640 + IT_0641 + IT_0845);
    const ccomplex_t IT_0847 = IT_0684 + IT_0787 + IT_0796;
    const ccomplex_t IT_0848 = 0.125*IT_0478*((-16)*IT_0106 + IT_0479 +
       IT_0643) + 0.125*IT_0559*(IT_0560 + IT_0562 + IT_0644) + IT_0414*IT_0841 
      + 0.125*IT_0649*IT_0844 + 0.125*IT_0639*(IT_0640 + IT_0641 + IT_0845) +
       0.125*IT_0828*IT_0847;
    const ccomplex_t IT_0849 = 8*conj(IT_0420);
    const ccomplex_t IT_0850 = (-2)*conj(IT_0055);
    const ccomplex_t IT_0851 = (-2)*conj(IT_0082);
    const ccomplex_t IT_0852 = (conj(IT_0055) + conj(IT_0180))*IT_0669 + (
      -0.25)*IT_0682*(IT_0553 + IT_0695) + 2*IT_0478*(conj(IT_0082) + conj
      (IT_0182) + (-0.125)*IT_0717) + (-0.25)*IT_0666*(IT_0416 + IT_0418 +
       IT_0849) + (-0.25)*IT_0639*(IT_0690 + IT_0850) + (-0.25)*IT_0649*(IT_0686
       + IT_0687 + IT_0851);
    const ccomplex_t IT_0853 = (-4)*IT_0559*(conj(IT_0082) + 2*conj(IT_0110) +
       conj(IT_0182) + (-0.25)*IT_0561) + (-16)*IT_0478*(conj(IT_0106) + (
      -0.0625)*IT_0480 + (-0.0625)*IT_0689) + (IT_0694 + IT_0799 + IT_0811)
      *IT_0828 + IT_0414*(IT_0416 + IT_0418 + IT_0849) + (-4)*IT_0649*(conj
      (IT_0106) + (-0.25)*IT_0690 + (-0.25)*IT_0850) + (-4)*IT_0639*(conj
      (IT_0110) + (-0.25)*IT_0686 + (-0.25)*IT_0687 + (-0.25)*IT_0851);
    const ccomplex_t IT_0854 = 0.166666666666667*IT_0524;
    const ccomplex_t IT_0855 = 0.166666666666667*IT_0564;
    const ccomplex_t IT_0856 = -IT_0550;
    const ccomplex_t IT_0857 = IT_0217 + IT_0234;
    const ccomplex_t IT_0858 = 0.25*IT_0857;
    const ccomplex_t IT_0859 = 0.25*IT_0786;
    const ccomplex_t IT_0860 = 0.25*IT_0836;
    const ccomplex_t IT_0861 = 0.25*IT_0818;
    const ccomplex_t IT_0862 = (-0.5)*(IT_0398 + IT_0401 + 2*IT_0420)*IT_0524 
      + (-0.666666666666667)*IT_0106*IT_0682 + IT_0521*(0.5*IT_0218 + 0.5
      *IT_0222 + IT_0683) + 0.333333333333333*IT_0397*(IT_0420 + 0.25*IT_0587 +
       0.25*IT_0731) + (-0.333333333333333)*(IT_0276 + -IT_0305 + IT_0320 + 0.25
      *IT_0819)*IT_0856 + 0.333333333333333*(IT_0223 + IT_0224 + IT_0235 +
       IT_0236)*IT_0858 + 0.333333333333333*IT_0847*IT_0859 + 0.333333333333333*
      (IT_0801 + IT_0802 + IT_0805 + IT_0806 + IT_0807)*IT_0860 +
       0.333333333333333*(IT_0360 + IT_0361 + IT_0371 + IT_0372 + IT_0373)
      *IT_0861;
    const ccomplex_t IT_0863 = (-0.75)*IT_0521*(conj(IT_0218) +
       0.166666666666667*IT_0522) + 0.75*IT_0524*(conj(IT_0401) + 2*conj(IT_0420
      ) + (-0.166666666666667)*IT_0525) + conj(IT_0106)*IT_0682 + (-0.5)*IT_0397
      *(conj(IT_0420) + 0.25*IT_0588 + 0.25*IT_0733) + 0.5*(conj(IT_0276) + 
      -IT_0301 + conj(IT_0320) + 0.25*IT_0551)*IT_0856 + (-0.5)*(IT_0219 +
       IT_0220 + IT_0226 + IT_0227)*IT_0858 + (-0.5)*(IT_0799 + IT_0811)*IT_0859
       + (-0.5)*(IT_0803 + IT_0804 + IT_0808 + IT_0809 + IT_0810)*IT_0860 + (
      -0.5)*(IT_0365 + IT_0366 + IT_0377 + IT_0378 + IT_0379)*IT_0861;
    const ccomplex_t IT_0864 = IT_0039*IT_0536;
    const ccomplex_t IT_0865 = (-4)*IT_0202;
    const ccomplex_t IT_0866 = (-4)*conj(IT_0202);
    const ccomplex_t IT_0867 = (-4)*IT_0169;
    const ccomplex_t IT_0868 = 0.25*IT_0867;
    const ccomplex_t IT_0869 = (-4)*conj(IT_0169);
    const ccomplex_t IT_0870 = 0.25*IT_0869;
    const ccomplex_t IT_0871 = (IT_0200 + IT_0305 + IT_0345 + IT_0347 +
       IT_0349)*IT_0540 + (conj(IT_0200) + IT_0301 + conj(IT_0345) + conj
      (IT_0347) + conj(IT_0349))*(IT_0169 + 0.25*IT_0541) + 0.25*(conj(IT_0288) 
      + conj(IT_0336) + conj(IT_0341))*IT_0542 + 0.25*(IT_0288 + IT_0336 +
       IT_0341)*IT_0543 + (conj(IT_0276) + conj(IT_0334))*(IT_0200 + 0.25
      *IT_0865) + (IT_0276 + IT_0334)*(conj(IT_0200) + 0.25*IT_0866) + (conj
      (IT_0320) + conj(IT_0352))*(IT_0167 + IT_0868) + conj(IT_0202)*(IT_0167 +
       IT_0288 + IT_0336 + IT_0341 + IT_0868) + (IT_0320 + IT_0352)*(conj
      (IT_0167) + IT_0870) + IT_0202*(conj(IT_0167) + conj(IT_0288) + conj
      (IT_0336) + conj(IT_0341) + IT_0870);
    const ccomplex_t IT_0872 = s_23*s_24*s_34*IT_0005;
    const ccomplex_t IT_0873 = 0.5*IT_0526 + (-0.5)*IT_0544 + -IT_0872;
    const ccomplex_t IT_0874 = IT_0404 + IT_0873;
    const ccomplex_t IT_0875 = IT_0005*IT_0874;
    const ccomplex_t IT_0876 = IT_0039 + IT_0528 + IT_0875;
    const ccomplex_t IT_0877 = m_N_2*IT_0876;
    const ccomplex_t IT_0878 = (-0.25)*IT_0568;
    const ccomplex_t IT_0879 = (-4)*conj(IT_0420);
    const ccomplex_t IT_0880 = (-0.25)*IT_0879;
    const ccomplex_t IT_0881 = 0.25*IT_0342*IT_0415 + 0.25*IT_0343*IT_0416 +
       0.25*IT_0350*IT_0417 + 0.25*IT_0351*IT_0418 + (conj(IT_0336) + conj
      (IT_0341) + conj(IT_0345) + conj(IT_0347) + conj(IT_0349))*IT_0420 + 
      (IT_0336 + IT_0341 + IT_0345 + IT_0347 + IT_0349)*conj(IT_0420) + -conj
      (IT_0334)*(IT_0398 + IT_0878) + -conj(IT_0352)*(IT_0401 + IT_0878) + 
      -IT_0334*(conj(IT_0398) + IT_0880) + -IT_0352*(conj(IT_0401) + IT_0880);
    const ccomplex_t IT_0882 = IT_0015*(IT_0111*IT_0171 + IT_0172*IT_0175 + 
      (conj(IT_0106) + conj(IT_0110) + conj(IT_0176) + conj(IT_0180) + conj
      (IT_0182))*IT_0204 + (IT_0106 + IT_0110 + IT_0176 + IT_0180 + IT_0182)
      *IT_0207) + IT_0217*(IT_0218*IT_0221 + conj(IT_0222)*IT_0225 + IT_0222
      *IT_0228) + IT_0234*(IT_0221*IT_0222 + IT_0218*IT_0228 + conj(IT_0222)
      *IT_0237) + ((IT_0106 + IT_0110 + IT_0176 + IT_0180 + IT_0182)*IT_0221 +
       IT_0111*IT_0225 + IT_0172*IT_0228 + (conj(IT_0106) + conj(IT_0110) + conj
      (IT_0176) + conj(IT_0180) + conj(IT_0182))*IT_0237)*IT_0249 + (-8)*IT_0260
      *IT_0327 + (-8)*IT_0333*IT_0353 + IT_0359*(conj(IT_0082)*IT_0364 + IT_0082
      *IT_0369 + IT_0370*IT_0375 + IT_0376*IT_0381 + IT_0382 + IT_0383 + IT_0384
       + IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 + IT_0391) +
       IT_0397*(conj(IT_0398)*IT_0399 + IT_0398*IT_0400 + IT_0401*IT_0402) +
       IT_0414*(IT_0370*IT_0415 + IT_0376*IT_0416 + conj(IT_0082)*IT_0417 +
       IT_0082*IT_0418 + IT_0419*IT_0421 + (IT_0082 + IT_0176 + IT_0182)*IT_0422
      ) + IT_0423*IT_0430 + IT_0437*IT_0470 + IT_0478*(conj(IT_0182)*IT_0479 +
       IT_0182*IT_0480) + (s_12 + IT_0485 + IT_0486 + IT_0494)*(conj(IT_0276)
      *IT_0440 + IT_0276*IT_0442 + IT_0320*IT_0453 + conj(IT_0288)*IT_0462 +
       conj(IT_0320)*IT_0496 + IT_0301*IT_0497 + IT_0300*IT_0498 + conj(IT_0300)
      *IT_0499) + (s_12 + IT_0505 + IT_0506 + IT_0507)*IT_0512 + (s_12 + IT_0485
       + IT_0486 + IT_0516)*IT_0517 + IT_0218*(IT_0521*IT_0522 + IT_0524*IT_0525
      ) + IT_0401*(IT_0531*IT_0532 + IT_0524*IT_0533) + 4*IT_0537*(IT_0167*conj
      (IT_0167) + IT_0200*conj(IT_0200) + IT_0169*IT_0540 + 0.25*conj(IT_0169)
      *IT_0541 + 0.25*conj(IT_0202)*IT_0542 + IT_0202*(conj(IT_0202) + 0.25
      *IT_0543)) + 4*IT_0550*(conj(IT_0222)*IT_0320 + IT_0222*conj(IT_0320) +
       IT_0218*(conj(IT_0276) + 0.25*IT_0551) + 0.25*IT_0301*IT_0552 + 0.25
      *IT_0305*IT_0553) + IT_0559*(IT_0370*IT_0560 + IT_0182*IT_0561 + conj
      (IT_0182)*IT_0562) + conj(IT_0398)*IT_0564*IT_0565 + IT_0570 + IT_0578
      *IT_0593 + IT_0599*IT_0614 + IT_0624*IT_0625 + IT_0629*(IT_0630 + IT_0631 
      + IT_0632 + IT_0633) + IT_0637*IT_0638 + conj(IT_0082)*(IT_0639*IT_0642 +
       IT_0478*IT_0643 + IT_0559*IT_0644 + IT_0649*IT_0650) + conj(IT_0055)*
      (IT_0642*IT_0649 + IT_0639*IT_0650 + IT_0655*IT_0658 + (IT_0417 + IT_0421)
      *IT_0666 + IT_0669*IT_0670 + IT_0682*(IT_0684 + IT_0685)) + IT_0082*
      (IT_0639*IT_0688 + IT_0478*IT_0689 + IT_0649*IT_0690) + IT_0055*((IT_0418 
      + IT_0422)*IT_0666 + IT_0649*IT_0688 + IT_0639*IT_0690 + IT_0655*IT_0693 +
       IT_0682*(IT_0694 + IT_0695)) + IT_0698*(IT_0207*IT_0398 + IT_0204*conj
      (IT_0398) + IT_0175*IT_0401 + IT_0171*conj(IT_0401) + IT_0699 + IT_0700 +
       IT_0701 + IT_0702 + IT_0703 + IT_0704 + IT_0705 + IT_0706) + (IT_0175
      *IT_0398 + IT_0171*conj(IT_0398) + IT_0207*IT_0401 + IT_0204*conj(IT_0401)
       + IT_0699 + IT_0700 + IT_0701 + IT_0702 + IT_0703 + IT_0704 + IT_0705 +
       IT_0706)*IT_0713 + conj(IT_0180)*((IT_0415 + IT_0421)*IT_0666 + IT_0682*
      (IT_0552 + IT_0684) + IT_0478*IT_0714 + IT_0669*IT_0715 + IT_0655*IT_0716)
       + IT_0180*((IT_0416 + IT_0422)*IT_0666 + IT_0682*(IT_0553 + IT_0694) +
       IT_0478*IT_0717 + IT_0655*IT_0718) + IT_0726*IT_0735 + IT_0743*(conj
      (IT_0055)*IT_0745 + IT_0055*IT_0747 + conj(IT_0180)*IT_0749 + IT_0180
      *IT_0751 + IT_0752 + IT_0753 + IT_0754 + IT_0755 + IT_0756 + IT_0757 +
       IT_0758 + IT_0759 + IT_0760 + IT_0761) + (conj(IT_0180)*IT_0744 + IT_0180
      *IT_0746 + conj(IT_0055)*IT_0748 + IT_0055*IT_0750 + IT_0752 + IT_0753 +
       IT_0754 + IT_0755 + IT_0756 + IT_0757 + IT_0758 + IT_0759 + IT_0760 +
       IT_0761)*IT_0764 + IT_0769*(IT_0370*IT_0745 + IT_0376*IT_0747 + conj
      (IT_0082)*IT_0749 + IT_0082*IT_0751 + IT_0770 + IT_0771 + IT_0772 +
       IT_0773 + IT_0774 + IT_0775 + IT_0776 + IT_0777 + IT_0778 + IT_0779) + 
      (conj(IT_0082)*IT_0745 + IT_0082*IT_0747 + IT_0370*IT_0749 + IT_0376
      *IT_0751 + IT_0770 + IT_0771 + IT_0772 + IT_0773 + IT_0774 + IT_0775 +
       IT_0776 + IT_0777 + IT_0778 + IT_0779)*IT_0784 + conj(IT_0222)*IT_0786
      *IT_0787 + IT_0794*IT_0812 + (conj(IT_0218)*IT_0362 + IT_0218*IT_0367 +
       IT_0350*IT_0399 + IT_0351*IT_0400 + conj(IT_0352)*(IT_0787 + IT_0795) +
       conj(IT_0334)*IT_0797 + IT_0334*IT_0800 + IT_0352*(IT_0798 + IT_0811))
      *IT_0818 + IT_0820 + (IT_0370*IT_0787 + IT_0419*IT_0795 + conj(IT_0082)
      *IT_0796 + (IT_0082 + IT_0176 + IT_0182)*IT_0798 + IT_0082*IT_0799 +
       IT_0376*IT_0811)*IT_0828 + IT_0839 + 8*conj(IT_0106)*IT_0846 + 8*conj
      (IT_0110)*IT_0848 + (-4)*IT_0106*IT_0852 + IT_0110*IT_0853 + 6*conj
      (IT_0401)*(IT_0398*IT_0531 + 0.166666666666667*IT_0397*IT_0566 + IT_0567
      *IT_0854 + IT_0569*IT_0855) + IT_0420*(IT_0397*(IT_0400 + IT_0402) + 6
      *IT_0531*(conj(IT_0401) + 0.166666666666667*IT_0532) + 6*IT_0533*IT_0854 +
       6*(IT_0590 + IT_0592)*IT_0855) + 12*conj(IT_0683)*IT_0862 + (-8)*IT_0683
      *IT_0863 + 4*IT_0864*IT_0871 + 4*IT_0877*IT_0881;
    return create_ccomplex_return(IT_0882);
}

