#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_Z_Z(
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
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
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
    const ccomplex_t IT_0014 = m_N_1*IT_0013;
    const ccomplex_t IT_0015 = (-2)*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0018, -1);
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0026 = IT_0019*IT_0025;
    const ccomplex_t IT_0027 = IT_0023*IT_0025;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + IT_0024 + 
      -IT_0026 + -IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0031 = IT_0019*IT_0030;
    const ccomplex_t IT_0032 = IT_0023*IT_0030;
    const ccomplex_t IT_0033 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0034 = IT_0019*IT_0033;
    const ccomplex_t IT_0035 = IT_0023*IT_0033;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0031 + IT_0032 + 
      -IT_0034 + -IT_0035);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0029*IT_0037;
    const ccomplex_t IT_0039 = pow(m_N_1, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0002 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0044 = IT_0019*IT_0043;
    const ccomplex_t IT_0045 = IT_0023*IT_0043;
    const ccomplex_t IT_0046 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0047 = IT_0019*IT_0046;
    const ccomplex_t IT_0048 = IT_0023*IT_0046;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0044 + IT_0045 + 
      -IT_0047 + -IT_0048);
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0052 = IT_0019*IT_0051;
    const ccomplex_t IT_0053 = IT_0023*IT_0051;
    const ccomplex_t IT_0054 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0055 = IT_0019*IT_0054;
    const ccomplex_t IT_0056 = IT_0023*IT_0054;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0052 + IT_0053 + 
      -IT_0055 + -IT_0056);
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0002 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0064 = IT_0019*IT_0063;
    const ccomplex_t IT_0065 = IT_0023*IT_0063;
    const ccomplex_t IT_0066 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0067 = IT_0019*IT_0066;
    const ccomplex_t IT_0068 = IT_0023*IT_0066;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0064 + IT_0065 + 
      -IT_0067 + -IT_0068);
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = cpow(IT_0070, 2);
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0002 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0071*IT_0073;
    const ccomplex_t IT_0075 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0076 = IT_0019*IT_0075;
    const ccomplex_t IT_0077 = IT_0023*IT_0075;
    const ccomplex_t IT_0078 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0079 = IT_0019*IT_0078;
    const ccomplex_t IT_0080 = IT_0023*IT_0078;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0076 + IT_0077 + 
      -IT_0079 + -IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0084 = IT_0019*IT_0083;
    const ccomplex_t IT_0085 = IT_0023*IT_0083;
    const ccomplex_t IT_0086 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0087 = IT_0019*IT_0086;
    const ccomplex_t IT_0088 = IT_0023*IT_0086;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0084 + IT_0085 + 
      -IT_0087 + -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = IT_0082*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_13 + IT_0002 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0091*IT_0093;
    const ccomplex_t IT_0095 = IT_0042 + IT_0062 + -IT_0074 + IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0059*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0038*IT_0100;
    const ccomplex_t IT_0102 = -IT_0098 + -IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0002 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0091*IT_0104;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0002 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0071*IT_0108;
    const ccomplex_t IT_0110 = conj(IT_0106) + conj(IT_0109);
    const ccomplex_t IT_0111 = conj(IT_0095) + conj(IT_0102) + IT_0110;
    const ccomplex_t IT_0112 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0113 = cos(alpha);
    const ccomplex_t IT_0114 = IT_0017*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0117 = sin(alpha);
    const ccomplex_t IT_0118 = IT_0017*IT_0117;
    const ccomplex_t IT_0119 = IT_0116*IT_0118;
    const ccomplex_t IT_0120 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0121 = IT_0022*IT_0113;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0124 = IT_0022*IT_0117;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0115 + -IT_0119 + 
      -IT_0122 + IT_0125);
    const ccomplex_t IT_0127 = cos(beta);
    const ccomplex_t IT_0128 = cpow(IT_0016, -4);
    const ccomplex_t IT_0129 = IT_0018*IT_0113*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = sin(beta);
    const ccomplex_t IT_0131 = IT_0018*IT_0117*IT_0128*IT_0130;
    const ccomplex_t IT_0132 = cpow(IT_0018, -3);
    const ccomplex_t IT_0133 = IT_0113*IT_0127*IT_0132;
    const ccomplex_t IT_0134 = IT_0117*IT_0130*IT_0132;
    const ccomplex_t IT_0135 = cpow(IT_0016, -2);
    const ccomplex_t IT_0136 = IT_0022*IT_0113*IT_0127*IT_0135;
    const ccomplex_t IT_0137 = IT_0022*IT_0117*IT_0130*IT_0135;
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
    const ccomplex_t IT_0148 = IT_0114*IT_0116;
    const ccomplex_t IT_0149 = IT_0112*IT_0118;
    const ccomplex_t IT_0150 = IT_0121*IT_0123;
    const ccomplex_t IT_0151 = IT_0120*IT_0124;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0148 + IT_0149 + 
      -IT_0150 + -IT_0151);
    const ccomplex_t IT_0153 = -IT_0152;
    const ccomplex_t IT_0154 = IT_0018*IT_0113*IT_0128*IT_0130;
    const ccomplex_t IT_0155 = IT_0018*IT_0117*IT_0127*IT_0128;
    const ccomplex_t IT_0156 = IT_0113*IT_0130*IT_0132;
    const ccomplex_t IT_0157 = IT_0117*IT_0127*IT_0132;
    const ccomplex_t IT_0158 = IT_0022*IT_0113*IT_0130*IT_0135;
    const ccomplex_t IT_0159 = IT_0022*IT_0117*IT_0127*IT_0135;
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
    const ccomplex_t IT_0172 = IT_0106 + IT_0109;
    const ccomplex_t IT_0173 = IT_0095 + IT_0102 + IT_0172;
    const ccomplex_t IT_0174 = 2*conj(IT_0167);
    const ccomplex_t IT_0175 = (-2)*conj(IT_0169);
    const ccomplex_t IT_0176 = IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = (-0.5)*IT_0069;
    const ccomplex_t IT_0178 = cpow(IT_0177, 2);
    const ccomplex_t IT_0179 = IT_0073*IT_0178;
    const ccomplex_t IT_0180 = IT_0094 + -IT_0179;
    const ccomplex_t IT_0181 = IT_0108*IT_0178;
    const ccomplex_t IT_0182 = -IT_0105 + IT_0181;
    const ccomplex_t IT_0183 = IT_0042 + IT_0062;
    const ccomplex_t IT_0184 = conj(IT_0102) + conj(IT_0180) + conj(IT_0182) +
       conj(IT_0183);
    const ccomplex_t IT_0185 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0186 = IT_0114*IT_0185;
    const ccomplex_t IT_0187 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0188 = IT_0118*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0190 = IT_0121*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0192 = IT_0124*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0186 + -IT_0188 + 
      -IT_0190 + IT_0192);
    const ccomplex_t IT_0194 = IT_0146*IT_0193;
    const ccomplex_t IT_0195 = IT_0114*IT_0187;
    const ccomplex_t IT_0196 = IT_0118*IT_0185;
    const ccomplex_t IT_0197 = IT_0121*IT_0191;
    const ccomplex_t IT_0198 = IT_0124*IT_0189;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + IT_0196 + 
      -IT_0197 + -IT_0198);
    const ccomplex_t IT_0200 = -IT_0199;
    const ccomplex_t IT_0201 = IT_0165*IT_0200;
    const ccomplex_t IT_0202 = (-0.5)*IT_0194 + (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = 2*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0194 + 0.5*IT_0201;
    const ccomplex_t IT_0205 = (-2)*IT_0204;
    const ccomplex_t IT_0206 = IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = IT_0102 + IT_0180 + IT_0182 + IT_0183;
    const ccomplex_t IT_0208 = 2*conj(IT_0202);
    const ccomplex_t IT_0209 = (-2)*conj(IT_0204);
    const ccomplex_t IT_0210 = IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = s_13*s_14*s_34*IT_0000;
    const ccomplex_t IT_0212 = IT_0039 + IT_0211;
    const ccomplex_t IT_0213 = pow(s_13, 2);
    const ccomplex_t IT_0214 = IT_0005*IT_0213;
    const ccomplex_t IT_0215 = pow(s_14, 2);
    const ccomplex_t IT_0216 = IT_0005*IT_0215;
    const ccomplex_t IT_0217 = -IT_0214 + -IT_0216;
    const ccomplex_t IT_0218 = IT_0212 + IT_0217;
    const ccomplex_t IT_0219 = m_N_1*IT_0218;
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = 2*IT_0105 + (-2)*IT_0109;
    const ccomplex_t IT_0222 = 2*IT_0098 + 2*IT_0101;
    const ccomplex_t IT_0223 = conj(IT_0221) + conj(IT_0222);
    const ccomplex_t IT_0224 = 2*IT_0169;
    const ccomplex_t IT_0225 = (-2)*IT_0167;
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = IT_0221 + IT_0222;
    const ccomplex_t IT_0228 = 2*conj(IT_0169);
    const ccomplex_t IT_0229 = (-2)*conj(IT_0167);
    const ccomplex_t IT_0230 = IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = 2*IT_0105 + (-2)*IT_0181;
    const ccomplex_t IT_0232 = conj(IT_0222) + conj(IT_0231);
    const ccomplex_t IT_0233 = 2*IT_0204;
    const ccomplex_t IT_0234 = (-2)*IT_0202;
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = IT_0222 + IT_0231;
    const ccomplex_t IT_0237 = 2*conj(IT_0204);
    const ccomplex_t IT_0238 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = s_23*IT_0002;
    const ccomplex_t IT_0241 = IT_0001*IT_0240;
    const ccomplex_t IT_0242 = s_24*IT_0002;
    const ccomplex_t IT_0243 = IT_0005*IT_0242;
    const ccomplex_t IT_0244 = s_23*s_34;
    const ccomplex_t IT_0245 = (-2)*IT_0244;
    const ccomplex_t IT_0246 = IT_0242 + IT_0245;
    const ccomplex_t IT_0247 = IT_0005*IT_0246;
    const ccomplex_t IT_0248 = (-0.5)*IT_0241 + 0.5*IT_0243 + (-0.5)*IT_0247;
    const ccomplex_t IT_0249 = s_24 + IT_0248;
    const ccomplex_t IT_0250 = m_N_1*IT_0249;
    const ccomplex_t IT_0251 = (-2)*IT_0250;
    const ccomplex_t IT_0252 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0253 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0254 = s_14*s_24*IT_0005;
    const ccomplex_t IT_0255 = -IT_0253 + -IT_0254;
    const ccomplex_t IT_0256 = s_12 + IT_0252 + IT_0255;
    const ccomplex_t IT_0257 = m_N_1*IT_0256;
    const ccomplex_t IT_0258 = s_12*s_14;
    const ccomplex_t IT_0259 = s_12*s_14*IT_0002;
    const ccomplex_t IT_0260 = s_24*IT_0002*IT_0039;
    const ccomplex_t IT_0261 = s_13*s_14*s_23;
    const ccomplex_t IT_0262 = s_23*s_34*IT_0039;
    const ccomplex_t IT_0263 = -IT_0260 + (-2)*IT_0261 + 2*IT_0262;
    const ccomplex_t IT_0264 = IT_0259 + IT_0263;
    const ccomplex_t IT_0265 = IT_0005*IT_0264;
    const ccomplex_t IT_0266 = (-0.5)*IT_0265;
    const ccomplex_t IT_0267 = s_24*IT_0039;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = IT_0258 + IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = 4*IT_0231;
    const ccomplex_t IT_0271 = 4*IT_0222;
    const ccomplex_t IT_0272 = IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = 4*conj(IT_0231);
    const ccomplex_t IT_0274 = 4*conj(IT_0222);
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 4*IT_0221;
    const ccomplex_t IT_0277 = IT_0271 + IT_0276;
    const ccomplex_t IT_0278 = 4*conj(IT_0221);
    const ccomplex_t IT_0279 = IT_0274 + IT_0278;
    const ccomplex_t IT_0280 = IT_0005*IT_0039;
    const ccomplex_t IT_0281 = pow(s_34, 2);
    const ccomplex_t IT_0282 = pow(m_Z, 4);
    const ccomplex_t IT_0283 = (-0.25)*IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = IT_0280*IT_0284;
    const ccomplex_t IT_0286 = (-4)*IT_0285;
    const ccomplex_t IT_0287 = (-2)*IT_0180;
    const ccomplex_t IT_0288 = (-2)*IT_0183;
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0291 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0292 = IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0042 + (-2)*IT_0062 + (-2)*IT_0094 + 2
      *IT_0179;
    const ccomplex_t IT_0294 = s_23*s_24*s_34*IT_0000;
    const ccomplex_t IT_0295 = IT_0039 + IT_0294;
    const ccomplex_t IT_0296 = pow(s_23, 2);
    const ccomplex_t IT_0297 = IT_0005*IT_0296;
    const ccomplex_t IT_0298 = pow(s_24, 2);
    const ccomplex_t IT_0299 = IT_0005*IT_0298;
    const ccomplex_t IT_0300 = -IT_0297 + -IT_0299;
    const ccomplex_t IT_0301 = IT_0295 + IT_0300;
    const ccomplex_t IT_0302 = m_N_1*IT_0301;
    const ccomplex_t IT_0303 = -IT_0302;
    const ccomplex_t IT_0304 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0305 = s_12 + IT_0255 + IT_0304;
    const ccomplex_t IT_0306 = m_N_1*IT_0305;
    const ccomplex_t IT_0307 = s_12*IT_0002;
    const ccomplex_t IT_0308 = s_13*s_23;
    const ccomplex_t IT_0309 = (-2)*IT_0308;
    const ccomplex_t IT_0310 = IT_0307 + IT_0309;
    const ccomplex_t IT_0311 = IT_0005*IT_0310;
    const ccomplex_t IT_0312 = (s_12 + -IT_0254)*(s_12 + (-0.5)*IT_0311);
    const ccomplex_t IT_0313 = (-2)*IT_0312;
    const ccomplex_t IT_0314 = 2*conj(IT_0231);
    const ccomplex_t IT_0315 = IT_0002*IT_0039;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = IT_0215 + IT_0316;
    const ccomplex_t IT_0318 = IT_0005*IT_0317;
    const ccomplex_t IT_0319 = -IT_0039;
    const ccomplex_t IT_0320 = IT_0216 + IT_0319;
    const ccomplex_t IT_0321 = (-4)*IT_0320;
    const ccomplex_t IT_0322 = IT_0318 + IT_0321;
    const ccomplex_t IT_0323 = m_N_1*IT_0322;
    const ccomplex_t IT_0324 = m_N_1*IT_0071;
    const ccomplex_t IT_0325 = IT_0070*IT_0177;
    const ccomplex_t IT_0326 = m_N_1*IT_0325;
    const ccomplex_t IT_0327 = IT_0108*(IT_0324 + -IT_0326);
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = 2*IT_0328;
    const ccomplex_t IT_0330 = cpow(IT_0090, 2);
    const ccomplex_t IT_0331 = m_N_2*IT_0330;
    const ccomplex_t IT_0332 = m_N_1*IT_0091;
    const ccomplex_t IT_0333 = IT_0104*(IT_0331 + -IT_0332);
    const ccomplex_t IT_0334 = -IT_0333;
    const ccomplex_t IT_0335 = 2*IT_0334;
    const ccomplex_t IT_0336 = IT_0329 + IT_0335;
    const ccomplex_t IT_0337 = 2*conj(IT_0328);
    const ccomplex_t IT_0338 = 2*conj(IT_0334);
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = cpow(IT_0058, 2);
    const ccomplex_t IT_0341 = m_N_4*IT_0340;
    const ccomplex_t IT_0342 = m_N_1*IT_0059;
    const ccomplex_t IT_0343 = -IT_0342;
    const ccomplex_t IT_0344 = IT_0341 + IT_0343;
    const ccomplex_t IT_0345 = IT_0097*IT_0344;
    const ccomplex_t IT_0346 = cpow(IT_0037, 2);
    const ccomplex_t IT_0347 = m_N_3*IT_0346;
    const ccomplex_t IT_0348 = m_N_1*IT_0038;
    const ccomplex_t IT_0349 = -IT_0348;
    const ccomplex_t IT_0350 = IT_0347 + IT_0349;
    const ccomplex_t IT_0351 = IT_0100*IT_0350;
    const ccomplex_t IT_0352 = IT_0345 + IT_0351;
    const ccomplex_t IT_0353 = (-2)*IT_0231;
    const ccomplex_t IT_0354 = (-2)*IT_0222;
    const ccomplex_t IT_0355 = IT_0353 + IT_0354;
    const ccomplex_t IT_0356 = (-2)*conj(IT_0231);
    const ccomplex_t IT_0357 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0358 = IT_0356 + IT_0357;
    const ccomplex_t IT_0359 = cpow(IT_0029, 2);
    const ccomplex_t IT_0360 = m_N_3*IT_0359;
    const ccomplex_t IT_0361 = IT_0349 + IT_0360;
    const ccomplex_t IT_0362 = IT_0100*IT_0361;
    const ccomplex_t IT_0363 = cpow(IT_0050, 2);
    const ccomplex_t IT_0364 = m_N_4*IT_0363;
    const ccomplex_t IT_0365 = IT_0343 + IT_0364;
    const ccomplex_t IT_0366 = IT_0097*IT_0365;
    const ccomplex_t IT_0367 = IT_0362 + IT_0366;
    const ccomplex_t IT_0368 = (-2)*IT_0221;
    const ccomplex_t IT_0369 = IT_0354 + IT_0368;
    const ccomplex_t IT_0370 = (-2)*conj(IT_0221);
    const ccomplex_t IT_0371 = IT_0357 + IT_0370;
    const ccomplex_t IT_0372 = conj(IT_0328) + conj(IT_0334);
    const ccomplex_t IT_0373 = 2*IT_0231;
    const ccomplex_t IT_0374 = IT_0328 + IT_0334;
    const ccomplex_t IT_0375 = m_N_1*IT_0178;
    const ccomplex_t IT_0376 = -IT_0375;
    const ccomplex_t IT_0377 = IT_0326 + IT_0376;
    const ccomplex_t IT_0378 = IT_0108*IT_0377;
    const ccomplex_t IT_0379 = cpow(IT_0082, 2);
    const ccomplex_t IT_0380 = m_N_2*IT_0379;
    const ccomplex_t IT_0381 = -IT_0332;
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = IT_0104*IT_0382;
    const ccomplex_t IT_0384 = IT_0378 + -IT_0383;
    const ccomplex_t IT_0385 = 2*IT_0384;
    const ccomplex_t IT_0386 = 2*conj(IT_0384);
    const ccomplex_t IT_0387 = IT_0323*(conj(IT_0222)*IT_0336 + IT_0222
      *IT_0339 + conj(IT_0352)*IT_0355 + IT_0352*IT_0358 + conj(IT_0367)*IT_0369
       + IT_0367*IT_0371 + IT_0372*IT_0373 + IT_0314*IT_0374 + IT_0223*IT_0385 +
       IT_0227*IT_0386);
    const ccomplex_t IT_0388 = s_14*IT_0282;
    const ccomplex_t IT_0389 = IT_0000*IT_0388;
    const ccomplex_t IT_0390 = (-10)*s_14;
    const ccomplex_t IT_0391 = IT_0389 + IT_0390;
    const ccomplex_t IT_0392 = m_N_1*IT_0391;
    const ccomplex_t IT_0393 = (-2)*IT_0352;
    const ccomplex_t IT_0394 = IT_0336 + IT_0393;
    const ccomplex_t IT_0395 = (-2)*conj(IT_0352);
    const ccomplex_t IT_0396 = IT_0339 + IT_0395;
    const ccomplex_t IT_0397 = (-2)*IT_0367;
    const ccomplex_t IT_0398 = IT_0385 + IT_0397;
    const ccomplex_t IT_0399 = (-2)*conj(IT_0367);
    const ccomplex_t IT_0400 = IT_0386 + IT_0399;
    const ccomplex_t IT_0401 = conj(IT_0102)*IT_0329;
    const ccomplex_t IT_0402 = IT_0102*IT_0337;
    const ccomplex_t IT_0403 = conj(IT_0102)*IT_0385;
    const ccomplex_t IT_0404 = IT_0102*IT_0386;
    const ccomplex_t IT_0405 = conj(IT_0102)*IT_0335;
    const ccomplex_t IT_0406 = IT_0102*IT_0338;
    const ccomplex_t IT_0407 = 2*IT_0293;
    const ccomplex_t IT_0408 = 2*conj(IT_0293);
    const ccomplex_t IT_0409 = (-2)*IT_0042 + (-2)*IT_0062 + 2*IT_0074 + (-2)
      *IT_0094;
    const ccomplex_t IT_0410 = 2*IT_0409;
    const ccomplex_t IT_0411 = 2*conj(IT_0409);
    const ccomplex_t IT_0412 = -IT_0216;
    const ccomplex_t IT_0413 = IT_0039 + IT_0412;
    const ccomplex_t IT_0414 = IT_0039*IT_0413;
    const ccomplex_t IT_0415 = -IT_0414;
    const ccomplex_t IT_0416 = 6*conj(IT_0231);
    const ccomplex_t IT_0417 = s_14*s_24;
    const ccomplex_t IT_0418 = (-2)*IT_0417;
    const ccomplex_t IT_0419 = IT_0307 + IT_0418;
    const ccomplex_t IT_0420 = IT_0005*IT_0419;
    const ccomplex_t IT_0421 = (-0.5)*IT_0420;
    const ccomplex_t IT_0422 = s_12*IT_0282;
    const ccomplex_t IT_0423 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0424 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0425 = s_13*s_24*s_34;
    const ccomplex_t IT_0426 = (-2)*IT_0423 + (-2)*IT_0424 + 4*IT_0425;
    const ccomplex_t IT_0427 = IT_0422 + IT_0426;
    const ccomplex_t IT_0428 = IT_0000*IT_0427;
    const ccomplex_t IT_0429 = 0.25*IT_0428;
    const ccomplex_t IT_0430 = (-0.5)*IT_0311;
    const ccomplex_t IT_0431 = -IT_0061*(IT_0342 + -IT_0364);
    const ccomplex_t IT_0432 = -IT_0431;
    const ccomplex_t IT_0433 = 8*IT_0432;
    const ccomplex_t IT_0434 = -IT_0324;
    const ccomplex_t IT_0435 = IT_0326 + IT_0434;
    const ccomplex_t IT_0436 = IT_0073*IT_0435;
    const ccomplex_t IT_0437 = IT_0093*IT_0332;
    const ccomplex_t IT_0438 = IT_0041*IT_0361;
    const ccomplex_t IT_0439 = IT_0436 + IT_0437 + -IT_0438;
    const ccomplex_t IT_0440 = 8*IT_0439;
    const ccomplex_t IT_0441 = IT_0433 + IT_0440;
    const ccomplex_t IT_0442 = IT_0093*IT_0380;
    const ccomplex_t IT_0443 = 8*IT_0442;
    const ccomplex_t IT_0444 = (-8)*IT_0432;
    const ccomplex_t IT_0445 = (-8)*IT_0439;
    const ccomplex_t IT_0446 = IT_0443 + IT_0444 + IT_0445;
    const ccomplex_t IT_0447 = IT_0093*IT_0331;
    const ccomplex_t IT_0448 = 8*IT_0447;
    const ccomplex_t IT_0449 = IT_0061*(IT_0341 + -IT_0342);
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = (-8)*IT_0450;
    const ccomplex_t IT_0452 = IT_0041*IT_0350;
    const ccomplex_t IT_0453 = IT_0073*IT_0377;
    const ccomplex_t IT_0454 = IT_0437 + -IT_0452 + IT_0453;
    const ccomplex_t IT_0455 = (-8)*IT_0454;
    const ccomplex_t IT_0456 = IT_0448 + IT_0451 + IT_0455;
    const ccomplex_t IT_0457 = (-8)*conj(IT_0450);
    const ccomplex_t IT_0458 = (-8)*conj(IT_0454);
    const ccomplex_t IT_0459 = IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = 8*conj(IT_0439);
    const ccomplex_t IT_0461 = conj(IT_0450) + conj(IT_0454);
    const ccomplex_t IT_0462 = 8*IT_0454;
    const ccomplex_t IT_0463 = 8*conj(IT_0454);
    const ccomplex_t IT_0464 = (-8)*conj(IT_0432);
    const ccomplex_t IT_0465 = (-8)*conj(IT_0439);
    const ccomplex_t IT_0466 = conj(IT_0439)*IT_0440 + conj(IT_0432)*IT_0441 +
       conj(IT_0442)*IT_0446 + conj(IT_0447)*IT_0456 + IT_0447*IT_0459 + IT_0432
      *IT_0460 + IT_0461*IT_0462 + 8*IT_0450*(conj(IT_0450) + 0.125*IT_0463) +
       IT_0442*(IT_0464 + IT_0465);
    const ccomplex_t IT_0467 = s_14*s_23*s_34;
    const ccomplex_t IT_0468 = s_12*IT_0281;
    const ccomplex_t IT_0469 = (-2)*IT_0423 + (-2)*IT_0424 + 2*IT_0425 + 2
      *IT_0467 + (-2)*IT_0468;
    const ccomplex_t IT_0470 = IT_0422 + IT_0469;
    const ccomplex_t IT_0471 = IT_0000*IT_0470;
    const ccomplex_t IT_0472 = (-0.25)*IT_0471;
    const ccomplex_t IT_0473 = -IT_0254;
    const ccomplex_t IT_0474 = -IT_0253;
    const ccomplex_t IT_0475 = 8*IT_0328;
    const ccomplex_t IT_0476 = 8*IT_0334;
    const ccomplex_t IT_0477 = IT_0475 + IT_0476;
    const ccomplex_t IT_0478 = 8*conj(IT_0442);
    const ccomplex_t IT_0479 = IT_0464 + IT_0465 + IT_0478;
    const ccomplex_t IT_0480 = (-8)*conj(IT_0328);
    const ccomplex_t IT_0481 = (-8)*conj(IT_0334);
    const ccomplex_t IT_0482 = IT_0480 + IT_0481;
    const ccomplex_t IT_0483 = 8*IT_0384;
    const ccomplex_t IT_0484 = 8*conj(IT_0384);
    const ccomplex_t IT_0485 = IT_0450 + IT_0454;
    const ccomplex_t IT_0486 = 8*conj(IT_0334);
    const ccomplex_t IT_0487 = (-8)*IT_0384;
    const ccomplex_t IT_0488 = (-8)*conj(IT_0384);
    const ccomplex_t IT_0489 = (-8)*IT_0328;
    const ccomplex_t IT_0490 = (-8)*IT_0334;
    const ccomplex_t IT_0491 = conj(IT_0384)*IT_0440 + conj(IT_0367)*IT_0446 +
       conj(IT_0352)*IT_0456 + IT_0352*IT_0459 + IT_0384*IT_0460 + 8*conj
      (IT_0328)*(IT_0450 + 0.125*IT_0462) + IT_0328*IT_0463 + conj(IT_0454)
      *IT_0476 + conj(IT_0450)*IT_0477 + IT_0367*IT_0479 + IT_0447*IT_0482 +
       conj(IT_0432)*IT_0483 + IT_0432*IT_0484 + IT_0485*IT_0486 + conj(IT_0442)
      *IT_0487 + IT_0442*IT_0488 + 8*conj(IT_0447)*(IT_0352 + 0.125*IT_0489 +
       0.125*IT_0490);
    const ccomplex_t IT_0492 = (-2)*IT_0423 + (-2)*IT_0424 + 4*IT_0467;
    const ccomplex_t IT_0493 = IT_0422 + IT_0492;
    const ccomplex_t IT_0494 = IT_0000*IT_0493;
    const ccomplex_t IT_0495 = 0.25*IT_0494;
    const ccomplex_t IT_0496 = 8*IT_0367;
    const ccomplex_t IT_0497 = IT_0487 + IT_0496;
    const ccomplex_t IT_0498 = IT_0489 + IT_0490;
    const ccomplex_t IT_0499 = conj(IT_0334)*IT_0476 + conj(IT_0328)*IT_0477 +
       8*IT_0352*(conj(IT_0352) + 0.125*IT_0480 + 0.125*IT_0481) + conj(IT_0384)
      *IT_0483 + IT_0328*IT_0486 + IT_0367*IT_0488 + conj(IT_0367)*IT_0497 +
       conj(IT_0352)*IT_0498;
    const ccomplex_t IT_0500 = (-2)*IT_0281;
    const ccomplex_t IT_0501 = IT_0282 + IT_0500;
    const ccomplex_t IT_0502 = IT_0000*IT_0501;
    const ccomplex_t IT_0503 = 0.25*IT_0502;
    const ccomplex_t IT_0504 = 1 + IT_0503;
    const ccomplex_t IT_0505 = IT_0039*IT_0504;
    const ccomplex_t IT_0506 = 8*conj(IT_0432);
    const ccomplex_t IT_0507 = IT_0372*IT_0441 + conj(IT_0352)*IT_0446 + 8
      *conj(IT_0367)*(IT_0447 + 0.125*IT_0451 + 0.125*IT_0455) + IT_0367*IT_0459
       + IT_0352*IT_0479 + IT_0442*IT_0482 + IT_0461*IT_0483 + IT_0484*IT_0485 +
       IT_0447*IT_0488 + conj(IT_0447)*IT_0497 + conj(IT_0442)*IT_0498 + IT_0374
      *(IT_0460 + IT_0506);
    const ccomplex_t IT_0508 = pow(m_Z, 6);
    const ccomplex_t IT_0509 = s_12*IT_0508;
    const ccomplex_t IT_0510 = s_13*s_23*IT_0282;
    const ccomplex_t IT_0511 = (-2)*IT_0510;
    const ccomplex_t IT_0512 = IT_0509 + IT_0511;
    const ccomplex_t IT_0513 = IT_0000*IT_0512;
    const ccomplex_t IT_0514 = 0.25*IT_0513;
    const ccomplex_t IT_0515 = IT_0307 + IT_0514;
    const ccomplex_t IT_0516 = (-8)*IT_0180;
    const ccomplex_t IT_0517 = (-8)*conj(IT_0180);
    const ccomplex_t IT_0518 = (-8)*IT_0095;
    const ccomplex_t IT_0519 = (-8)*conj(IT_0095);
    const ccomplex_t IT_0520 = s_12 + IT_0473;
    const ccomplex_t IT_0521 = IT_0039*IT_0520;
    const ccomplex_t IT_0522 = (-6)*IT_0293;
    const ccomplex_t IT_0523 = (-6)*conj(IT_0293);
    const ccomplex_t IT_0524 = IT_0000*IT_0281;
    const ccomplex_t IT_0525 = 0.5*IT_0524;
    const ccomplex_t IT_0526 = 1 + IT_0525;
    const ccomplex_t IT_0527 = s_12*IT_0526;
    const ccomplex_t IT_0528 = 4*IT_0169;
    const ccomplex_t IT_0529 = (-4)*IT_0167;
    const ccomplex_t IT_0530 = IT_0528 + IT_0529;
    const ccomplex_t IT_0531 = 4*IT_0204;
    const ccomplex_t IT_0532 = (-4)*IT_0202;
    const ccomplex_t IT_0533 = IT_0531 + IT_0532;
    const ccomplex_t IT_0534 = 4*IT_0167;
    const ccomplex_t IT_0535 = 4*IT_0202;
    const ccomplex_t IT_0536 = (-4)*conj(IT_0167);
    const ccomplex_t IT_0537 = (-4)*conj(IT_0202);
    const ccomplex_t IT_0538 = IT_0039*IT_0526;
    const ccomplex_t IT_0539 = conj(IT_0432) + conj(IT_0439);
    const ccomplex_t IT_0540 = 4*conj(IT_0169);
    const ccomplex_t IT_0541 = (-4)*IT_0169;
    const ccomplex_t IT_0542 = 4*conj(IT_0167);
    const ccomplex_t IT_0543 = (-4)*conj(IT_0169);
    const ccomplex_t IT_0544 = IT_0432 + IT_0439;
    const ccomplex_t IT_0545 = 4*conj(IT_0204);
    const ccomplex_t IT_0546 = (-4)*IT_0204;
    const ccomplex_t IT_0547 = 4*conj(IT_0202);
    const ccomplex_t IT_0548 = (-4)*conj(IT_0204);
    const ccomplex_t IT_0549 = (conj(IT_0202) + conj(IT_0328) + conj(IT_0334) 
      + conj(IT_0450) + conj(IT_0454))*IT_0530 + IT_0533*(conj(IT_0384) +
       IT_0539) + (IT_0202 + IT_0328 + IT_0334 + IT_0450 + IT_0454)*(IT_0536 +
       IT_0540) + (conj(IT_0204) + conj(IT_0352) + conj(IT_0447))*(IT_0534 +
       IT_0541) + (IT_0204 + IT_0352 + IT_0447)*(IT_0542 + IT_0543) + (IT_0384 +
       IT_0544)*(IT_0537 + IT_0545) + (conj(IT_0367) + conj(IT_0442))*(IT_0535 +
       IT_0546) + (IT_0367 + IT_0442)*(IT_0547 + IT_0548);
    const ccomplex_t IT_0550 = s_13*s_14*s_34*IT_0005;
    const ccomplex_t IT_0551 = 0.5*IT_0215 + (-0.5)*IT_0315 + -IT_0550;
    const ccomplex_t IT_0552 = IT_0213 + IT_0551;
    const ccomplex_t IT_0553 = IT_0005*IT_0552;
    const ccomplex_t IT_0554 = IT_0039 + IT_0412 + IT_0553;
    const ccomplex_t IT_0555 = m_N_1*IT_0554;
    const ccomplex_t IT_0556 = (-4)*IT_0231;
    const ccomplex_t IT_0557 = (-4)*conj(IT_0231);
    const ccomplex_t IT_0558 = (-4)*IT_0432;
    const ccomplex_t IT_0559 = (-4)*IT_0439;
    const ccomplex_t IT_0560 = (-4)*conj(IT_0432);
    const ccomplex_t IT_0561 = (-4)*conj(IT_0439);
    const ccomplex_t IT_0562 = IT_0279*IT_0442 + IT_0277*conj(IT_0442) +
       IT_0275*IT_0447 + IT_0272*conj(IT_0447) + (-4)*IT_0461*(IT_0222 + (-0.25)
      *IT_0556) + (-4)*IT_0485*(conj(IT_0222) + (-0.25)*IT_0557) + IT_0223*
      (IT_0558 + IT_0559) + IT_0227*(IT_0560 + IT_0561);
    const ccomplex_t IT_0563 = s_23*s_24*s_34*IT_0005;
    const ccomplex_t IT_0564 = 0.5*IT_0298 + (-0.5)*IT_0315 + -IT_0563;
    const ccomplex_t IT_0565 = IT_0296 + IT_0564;
    const ccomplex_t IT_0566 = IT_0005*IT_0565;
    const ccomplex_t IT_0567 = -IT_0299;
    const ccomplex_t IT_0568 = IT_0039 + IT_0566 + IT_0567;
    const ccomplex_t IT_0569 = m_N_1*IT_0568;
    const ccomplex_t IT_0570 = 4*IT_0409;
    const ccomplex_t IT_0571 = 4*conj(IT_0409);
    const ccomplex_t IT_0572 = (-4)*IT_0293;
    const ccomplex_t IT_0573 = (-4)*conj(IT_0293);
    const ccomplex_t IT_0574 = (-4)*IT_0409;
    const ccomplex_t IT_0575 = (-4)*conj(IT_0409);
    const ccomplex_t IT_0576 = conj(IT_0293)*IT_0384 + IT_0293*conj(IT_0384) +
       0.25*IT_0372*IT_0570 + 0.25*IT_0374*IT_0571 + 0.25*conj(IT_0367)*IT_0572 
      + 0.25*IT_0367*IT_0573 + 0.25*conj(IT_0352)*IT_0574 + 0.25*IT_0352*IT_0575;
    const ccomplex_t IT_0577 = s_12*s_24;
    const ccomplex_t IT_0578 = s_12*s_24*IT_0002;
    const ccomplex_t IT_0579 = s_14*IT_0002*IT_0039;
    const ccomplex_t IT_0580 = s_14*IT_0296;
    const ccomplex_t IT_0581 = s_12*s_23*s_34;
    const ccomplex_t IT_0582 = -IT_0579 + 2*IT_0580 + (-2)*IT_0581;
    const ccomplex_t IT_0583 = IT_0578 + IT_0582;
    const ccomplex_t IT_0584 = IT_0005*IT_0583;
    const ccomplex_t IT_0585 = 0.5*IT_0584;
    const ccomplex_t IT_0586 = -IT_0585;
    const ccomplex_t IT_0587 = s_14*IT_0039;
    const ccomplex_t IT_0588 = -IT_0587;
    const ccomplex_t IT_0589 = IT_0577 + IT_0586 + IT_0588;
    const ccomplex_t IT_0590 = 0.5*IT_0513;
    const ccomplex_t IT_0591 = (-4)*IT_0417;
    const ccomplex_t IT_0592 = IT_0307 + IT_0590 + IT_0591;
    const ccomplex_t IT_0593 = IT_0005*IT_0493;
    const ccomplex_t IT_0594 = -IT_0593;
    const ccomplex_t IT_0595 = IT_0592 + IT_0594;
    const ccomplex_t IT_0596 = (-4)*IT_0106;
    const ccomplex_t IT_0597 = (-4)*conj(IT_0106);
    const ccomplex_t IT_0598 = (-4)*IT_0109;
    const ccomplex_t IT_0599 = (-4)*IT_0182;
    const ccomplex_t IT_0600 = s_13*s_34*IT_0005;
    const ccomplex_t IT_0601 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0602 = (-2)*IT_0601;
    const ccomplex_t IT_0603 = IT_0388 + IT_0602;
    const ccomplex_t IT_0604 = IT_0000*IT_0603;
    const ccomplex_t IT_0605 = (-0.666666666666667)*IT_0600 + (
      -0.166666666666667)*IT_0604;
    const ccomplex_t IT_0606 = s_14 + IT_0605;
    const ccomplex_t IT_0607 = m_N_1*IT_0606;
    const ccomplex_t IT_0608 = (-6)*IT_0607;
    const ccomplex_t IT_0609 = 2*IT_0442;
    const ccomplex_t IT_0610 = (-2)*IT_0432;
    const ccomplex_t IT_0611 = (-2)*IT_0439;
    const ccomplex_t IT_0612 = IT_0609 + IT_0610 + IT_0611;
    const ccomplex_t IT_0613 = 2*conj(IT_0442);
    const ccomplex_t IT_0614 = (-2)*conj(IT_0432);
    const ccomplex_t IT_0615 = (-2)*conj(IT_0439);
    const ccomplex_t IT_0616 = IT_0613 + IT_0614 + IT_0615;
    const ccomplex_t IT_0617 = 2*IT_0447;
    const ccomplex_t IT_0618 = (-2)*IT_0450;
    const ccomplex_t IT_0619 = (-2)*IT_0454;
    const ccomplex_t IT_0620 = IT_0617 + IT_0618 + IT_0619;
    const ccomplex_t IT_0621 = 2*conj(IT_0447);
    const ccomplex_t IT_0622 = (-2)*conj(IT_0450);
    const ccomplex_t IT_0623 = (-2)*conj(IT_0454);
    const ccomplex_t IT_0624 = IT_0621 + IT_0622 + IT_0623;
    const ccomplex_t IT_0625 = conj(IT_0102)*IT_0609;
    const ccomplex_t IT_0626 = IT_0102*IT_0613;
    const ccomplex_t IT_0627 = conj(IT_0102)*IT_0617;
    const ccomplex_t IT_0628 = IT_0102*IT_0621;
    const ccomplex_t IT_0629 = conj(IT_0102)*IT_0610;
    const ccomplex_t IT_0630 = IT_0102*IT_0614;
    const ccomplex_t IT_0631 = conj(IT_0102)*IT_0618;
    const ccomplex_t IT_0632 = IT_0102*IT_0622;
    const ccomplex_t IT_0633 = conj(IT_0102)*IT_0611;
    const ccomplex_t IT_0634 = IT_0102*IT_0615;
    const ccomplex_t IT_0635 = conj(IT_0102)*IT_0619;
    const ccomplex_t IT_0636 = IT_0102*IT_0623;
    const ccomplex_t IT_0637 = IT_0608*(IT_0110*IT_0612 + IT_0172*IT_0616 +
       conj(IT_0182)*IT_0620 + IT_0182*IT_0624 + IT_0625 + IT_0626 + IT_0627 +
       IT_0628 + IT_0629 + IT_0630 + IT_0631 + IT_0632 + IT_0633 + IT_0634 +
       IT_0635 + IT_0636);
    const ccomplex_t IT_0638 = m_N_1*IT_0005;
    const ccomplex_t IT_0639 = s_24*IT_0282;
    const ccomplex_t IT_0640 = IT_0005*IT_0639;
    const ccomplex_t IT_0641 = (-2)*IT_0244 + (-0.5)*IT_0640;
    const ccomplex_t IT_0642 = IT_0242 + IT_0641;
    const ccomplex_t IT_0643 = IT_0638*IT_0642;
    const ccomplex_t IT_0644 = (-2)*IT_0643;
    const ccomplex_t IT_0645 = (conj(IT_0182)*IT_0612 + IT_0182*IT_0616 +
       IT_0110*IT_0620 + IT_0172*IT_0624 + IT_0625 + IT_0626 + IT_0627 + IT_0628
       + IT_0629 + IT_0630 + IT_0631 + IT_0632 + IT_0633 + IT_0634 + IT_0635 +
       IT_0636)*IT_0644;
    const ccomplex_t IT_0646 = 0.666666666666667*IT_0007 + (-0.166666666666667
      )*IT_0389;
    const ccomplex_t IT_0647 = s_14 + IT_0646;
    const ccomplex_t IT_0648 = m_N_1*IT_0647;
    const ccomplex_t IT_0649 = (-6)*IT_0648;
    const ccomplex_t IT_0650 = conj(IT_0102)*IT_0393;
    const ccomplex_t IT_0651 = IT_0102*IT_0395;
    const ccomplex_t IT_0652 = conj(IT_0102)*IT_0397;
    const ccomplex_t IT_0653 = IT_0102*IT_0399;
    const ccomplex_t IT_0654 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0655 = (-2)*IT_0654;
    const ccomplex_t IT_0656 = IT_0639 + IT_0655;
    const ccomplex_t IT_0657 = IT_0000*IT_0656;
    const ccomplex_t IT_0658 = (-0.666666666666667)*IT_0247 +
       0.166666666666667*IT_0657;
    const ccomplex_t IT_0659 = s_24 + IT_0658;
    const ccomplex_t IT_0660 = m_N_1*IT_0659;
    const ccomplex_t IT_0661 = 6*IT_0660;
    const ccomplex_t IT_0662 = s_14*s_23*s_34*IT_0005;
    const ccomplex_t IT_0663 = (-2)*IT_0308 + -IT_0417 + 2*IT_0662;
    const ccomplex_t IT_0664 = IT_0307 + IT_0663;
    const ccomplex_t IT_0665 = IT_0005*IT_0664;
    const ccomplex_t IT_0666 = -IT_0254 + (-0.5)*IT_0665;
    const ccomplex_t IT_0667 = s_12 + IT_0666;
    const ccomplex_t IT_0668 = m_N_1*IT_0667;
    const ccomplex_t IT_0669 = (-2)*IT_0668;
    const ccomplex_t IT_0670 = 2*IT_0222;
    const ccomplex_t IT_0671 = IT_0373 + IT_0670;
    const ccomplex_t IT_0672 = 2*conj(IT_0222);
    const ccomplex_t IT_0673 = IT_0314 + IT_0672;
    const ccomplex_t IT_0674 = 2*IT_0221;
    const ccomplex_t IT_0675 = IT_0670 + IT_0674;
    const ccomplex_t IT_0676 = 2*conj(IT_0221);
    const ccomplex_t IT_0677 = IT_0672 + IT_0676;
    const ccomplex_t IT_0678 = (-2)*IT_0384;
    const ccomplex_t IT_0679 = (-2)*conj(IT_0384);
    const ccomplex_t IT_0680 = (-2)*(IT_0222 + (-0.5)*IT_0368)*IT_0372 + (-2)*
      (conj(IT_0222) + (-0.5)*IT_0370)*IT_0374 + IT_0356*IT_0384 + IT_0353*conj
      (IT_0384) + conj(IT_0367)*IT_0671 + IT_0367*IT_0673 + conj(IT_0352)
      *IT_0675 + IT_0352*IT_0677 + conj(IT_0222)*IT_0678 + IT_0222*IT_0679;
    const ccomplex_t IT_0681 = (IT_0039 + -IT_0216)*(s_12 + (-0.5)*IT_0311);
    const ccomplex_t IT_0682 = 2*IT_0681;
    const ccomplex_t IT_0683 = (IT_0039 + -IT_0299)*(s_12 + (-0.5)*IT_0311);
    const ccomplex_t IT_0684 = 2*IT_0683;
    const ccomplex_t IT_0685 = (-2)*IT_0293;
    const ccomplex_t IT_0686 = s_24*IT_0213;
    const ccomplex_t IT_0687 = s_12*s_13*s_34;
    const ccomplex_t IT_0688 = -IT_0260 + 2*IT_0686 + (-2)*IT_0687;
    const ccomplex_t IT_0689 = IT_0259 + IT_0688;
    const ccomplex_t IT_0690 = IT_0005*IT_0689;
    const ccomplex_t IT_0691 = 2*IT_0267;
    const ccomplex_t IT_0692 = (-2)*IT_0258;
    const ccomplex_t IT_0693 = IT_0690 + IT_0691 + IT_0692;
    const ccomplex_t IT_0694 = IT_0005*IT_0388;
    const ccomplex_t IT_0695 = (-2)*IT_0008 + (-0.5)*IT_0694;
    const ccomplex_t IT_0696 = IT_0006 + IT_0695;
    const ccomplex_t IT_0697 = IT_0638*IT_0696;
    const ccomplex_t IT_0698 = (-2)*IT_0697;
    const ccomplex_t IT_0699 = (-0.666666666666667)*IT_0011 +
       0.166666666666667*IT_0604;
    const ccomplex_t IT_0700 = s_14 + IT_0699;
    const ccomplex_t IT_0701 = m_N_1*IT_0700;
    const ccomplex_t IT_0702 = 6*IT_0701;
    const ccomplex_t IT_0703 = s_23*s_34*IT_0005;
    const ccomplex_t IT_0704 = (-0.166666666666667)*IT_0657 + (
      -0.666666666666667)*IT_0703;
    const ccomplex_t IT_0705 = s_24 + IT_0704;
    const ccomplex_t IT_0706 = m_N_1*IT_0705;
    const ccomplex_t IT_0707 = (-6)*IT_0706;
    const ccomplex_t IT_0708 = IT_0000*IT_0639;
    const ccomplex_t IT_0709 = 0.666666666666667*IT_0243 + (-0.166666666666667
      )*IT_0708;
    const ccomplex_t IT_0710 = s_24 + IT_0709;
    const ccomplex_t IT_0711 = m_N_1*IT_0710;
    const ccomplex_t IT_0712 = (-6)*IT_0711;
    const ccomplex_t IT_0713 = IT_0610 + IT_0611;
    const ccomplex_t IT_0714 = (-10)*s_24;
    const ccomplex_t IT_0715 = IT_0708 + IT_0714;
    const ccomplex_t IT_0716 = m_N_1*IT_0715;
    const ccomplex_t IT_0717 = IT_0005*IT_0427;
    const ccomplex_t IT_0718 = -IT_0717;
    const ccomplex_t IT_0719 = IT_0592 + IT_0718;
    const ccomplex_t IT_0720 = (-4)*IT_0180;
    const ccomplex_t IT_0721 = s_13*s_23*s_24;
    const ccomplex_t IT_0722 = s_13*s_34*IT_0039;
    const ccomplex_t IT_0723 = -IT_0579 + (-2)*IT_0721 + 2*IT_0722;
    const ccomplex_t IT_0724 = IT_0578 + IT_0723;
    const ccomplex_t IT_0725 = IT_0005*IT_0724;
    const ccomplex_t IT_0726 = 2*IT_0587;
    const ccomplex_t IT_0727 = (-2)*IT_0577;
    const ccomplex_t IT_0728 = IT_0725 + IT_0726 + IT_0727;
    const ccomplex_t IT_0729 = IT_0614 + IT_0615;
    const ccomplex_t IT_0730 = (-2)*conj(IT_0293);
    const ccomplex_t IT_0731 = 2*(conj(IT_0102) + conj(IT_0182))*IT_0515 + (
      -0.25)*IT_0673*IT_0693 + (-0.25)*IT_0396*IT_0698 + (-0.25)*IT_0624*IT_0702
       + (-0.25)*IT_0400*IT_0707 + (-0.25)*IT_0613*IT_0716 + (conj(IT_0180) +
       conj(IT_0183))*IT_0719 + (-0.25)*IT_0712*IT_0729 + (-0.25)*IT_0728*IT_0730;
    const ccomplex_t IT_0732 = -IT_0417;
    const ccomplex_t IT_0733 = IT_0307 + IT_0732;
    const ccomplex_t IT_0734 = IT_0005*IT_0733;
    const ccomplex_t IT_0735 = -IT_0254 + (-0.25)*IT_0734;
    const ccomplex_t IT_0736 = s_12 + IT_0735;
    const ccomplex_t IT_0737 = m_N_1*IT_0736;
    const ccomplex_t IT_0738 = (-4)*IT_0737;
    const ccomplex_t IT_0739 = 2*IT_0432;
    const ccomplex_t IT_0740 = 2*IT_0439;
    const ccomplex_t IT_0741 = IT_0739 + IT_0740;
    const ccomplex_t IT_0742 = 2*conj(IT_0432);
    const ccomplex_t IT_0743 = 2*conj(IT_0439);
    const ccomplex_t IT_0744 = IT_0742 + IT_0743;
    const ccomplex_t IT_0745 = conj(IT_0352)*IT_0407 + IT_0352*IT_0408 + conj
      (IT_0367)*IT_0410 + IT_0367*IT_0411 + IT_0358*IT_0442 + IT_0355*conj
      (IT_0442) + IT_0371*IT_0447 + IT_0369*conj(IT_0447) + IT_0373*IT_0539 +
       IT_0314*IT_0544 + 2*IT_0461*(IT_0222 + 0.5*IT_0674) + 2*IT_0485*(conj
      (IT_0222) + 0.5*IT_0676) + conj(IT_0409)*IT_0678 + IT_0409*IT_0679 +
       IT_0372*IT_0685 + IT_0374*IT_0730 + conj(IT_0222)*IT_0741 + IT_0222
      *IT_0744;
    const ccomplex_t IT_0746 = s_13*s_24*s_34*IT_0005;
    const ccomplex_t IT_0747 = (-2)*IT_0308 + -IT_0417 + 2*IT_0746;
    const ccomplex_t IT_0748 = IT_0307 + IT_0747;
    const ccomplex_t IT_0749 = IT_0005*IT_0748;
    const ccomplex_t IT_0750 = -IT_0254 + (-0.5)*IT_0749;
    const ccomplex_t IT_0751 = s_12 + IT_0750;
    const ccomplex_t IT_0752 = m_N_1*IT_0751;
    const ccomplex_t IT_0753 = (-2)*IT_0752;
    const ccomplex_t IT_0754 = (-2)*IT_0409;
    const ccomplex_t IT_0755 = (-2)*conj(IT_0409);
    const ccomplex_t IT_0756 = IT_0298 + IT_0316;
    const ccomplex_t IT_0757 = IT_0005*IT_0756;
    const ccomplex_t IT_0758 = IT_0299 + IT_0319;
    const ccomplex_t IT_0759 = (-4)*IT_0758;
    const ccomplex_t IT_0760 = IT_0757 + IT_0759;
    const ccomplex_t IT_0761 = m_N_1*IT_0760;
    const ccomplex_t IT_0762 = 9*IT_0315;
    const ccomplex_t IT_0763 = (-2)*IT_0182;
    const ccomplex_t IT_0764 = (-2)*conj(IT_0182);
    const ccomplex_t IT_0765 = 1.125*IT_0039;
    const ccomplex_t IT_0766 = (-16)*IT_0439;
    const ccomplex_t IT_0767 = (-16)*conj(IT_0432);
    const ccomplex_t IT_0768 = (-16)*conj(IT_0439);
    const ccomplex_t IT_0769 = (-16)*IT_0384;
    const ccomplex_t IT_0770 = (-16)*conj(IT_0384);
    const ccomplex_t IT_0771 = 16*conj(IT_0447);
    const ccomplex_t IT_0772 = (16*conj(IT_0328) + 16*conj(IT_0334) + (-16)
      *conj(IT_0352))*IT_0367 + (16*IT_0328 + 16*IT_0334 + (-16)*IT_0352)*conj
      (IT_0367) + 16*conj(IT_0352)*IT_0384 + 16*IT_0352*conj(IT_0384) + (16*conj
      (IT_0432) + 16*conj(IT_0439) + (-16)*conj(IT_0442))*IT_0447 + IT_0442*((
      -16)*conj(IT_0447) + 16*conj(IT_0450) + 16*conj(IT_0454)) + IT_0461
      *IT_0766 + IT_0485*(16*conj(IT_0442) + IT_0767 + IT_0768) + IT_0372
      *IT_0769 + IT_0374*IT_0770 + IT_0432*((-16)*conj(IT_0450) + (-16)*conj
      (IT_0454) + IT_0771) + IT_0439*IT_0771;
    const ccomplex_t IT_0773 = (-0.5)*IT_0269;
    const ccomplex_t IT_0774 = (-0.5)*IT_0589;
    const ccomplex_t IT_0775 = IT_0286*IT_0289 + IT_0515*(IT_0516 + IT_0518) +
       IT_0595*(IT_0596 + IT_0598 + IT_0599) + (-2)*IT_0762*(IT_0172 + (-0.5)
      *IT_0763) + (-2)*(IT_0270 + 2*IT_0271 + IT_0276)*IT_0773 + (-2)*(IT_0572 +
       IT_0574)*IT_0774;
    const ccomplex_t IT_0776 = (-0.5)*IT_0764;
    const ccomplex_t IT_0777 = 2*conj(IT_0106)*IT_0515 + (-0.25)*IT_0677
      *IT_0693 + (-0.25)*IT_0400*IT_0698 + (-0.25)*IT_0616*IT_0702 + (-0.25)
      *IT_0396*IT_0707 + (-0.25)*(IT_0622 + IT_0623)*IT_0712 + (-0.25)*IT_0621
      *IT_0716 + conj(IT_0095)*IT_0719 + (-0.25)*IT_0728*IT_0755 + (-0.25)
      *IT_0292*IT_0762 + 0.5*IT_0286*(conj(IT_0102) + IT_0776);
    const ccomplex_t IT_0778 = 2*conj(IT_0102);
    const ccomplex_t IT_0779 = 0.25*IT_0286*(conj(IT_0095) + (-0.5)*IT_0290 + 
      (-0.5)*IT_0291) + IT_0515*(conj(IT_0183) + (-0.125)*IT_0517 + (-0.125)
      *IT_0519) + 0.25*(IT_0273 + 2*IT_0274 + IT_0278)*IT_0773 + 0.25*(IT_0573 +
       IT_0575)*IT_0774 + 0.5*IT_0595*(conj(IT_0109) + conj(IT_0182) + (-0.25)
      *IT_0597 + IT_0778) + 0.25*IT_0762*(IT_0110 + IT_0776 + IT_0778);
    const ccomplex_t IT_0780 = IT_0039 + IT_0567;
    const ccomplex_t IT_0781 = IT_0039*IT_0780;
    const ccomplex_t IT_0782 = -IT_0781;
    const ccomplex_t IT_0783 = (-0.166666666666667)*IT_0306;
    const ccomplex_t IT_0784 = (-0.166666666666667)*IT_0313;
    const ccomplex_t IT_0785 = (-0.125)*IT_0712;
    const ccomplex_t IT_0786 = (-0.125)*IT_0702;
    const ccomplex_t IT_0787 = (-0.125)*IT_0707;
    const ccomplex_t IT_0788 = (-0.125)*IT_0698;
    const ccomplex_t IT_0789 = (-0.125)*IT_0716;
    const ccomplex_t IT_0790 = (-0.125)*IT_0693;
    const ccomplex_t IT_0791 = (-0.125)*IT_0728;
    const ccomplex_t IT_0792 = (-8)*conj(IT_0183)*(IT_0182*IT_0515 + (-0.125)
      *IT_0719*IT_0720 + IT_0713*IT_0785 + IT_0620*IT_0786 + IT_0398*IT_0787 +
       IT_0394*IT_0788 + IT_0609*IT_0789 + IT_0671*IT_0790 + IT_0685*IT_0791);
    const ccomplex_t IT_0793 = IT_0106*IT_0515 + (-0.125)*IT_0289*IT_0762 + (
      -0.125)*IT_0286*IT_0763 + (IT_0618 + IT_0619)*IT_0785 + IT_0612*IT_0786 +
       IT_0394*IT_0787 + IT_0398*IT_0788 + IT_0617*IT_0789 + IT_0675*IT_0790 +
       IT_0754*IT_0791;
    const ccomplex_t IT_0794 = IT_0005*IT_0307;
    const ccomplex_t IT_0795 = (-0.5)*IT_0794;
    const ccomplex_t IT_0796 = s_14*s_23;
    const ccomplex_t IT_0797 = s_13*s_24;
    const ccomplex_t IT_0798 = s_12*s_34;
    const ccomplex_t IT_0799 = -IT_0797 + -IT_0798;
    const ccomplex_t IT_0800 = IT_0796 + IT_0799;
    const ccomplex_t IT_0801 = IT_0001*IT_0800;
    const ccomplex_t IT_0802 = (-0.25)*IT_0801;
    const ccomplex_t IT_0803 = s_12 + IT_0795 + IT_0802;
    const ccomplex_t IT_0804 = 8*IT_0167;
    const ccomplex_t IT_0805 = (-8)*IT_0169;
    const ccomplex_t IT_0806 = IT_0804 + IT_0805;
    const ccomplex_t IT_0807 = 8*conj(IT_0167);
    const ccomplex_t IT_0808 = (-8)*conj(IT_0169);
    const ccomplex_t IT_0809 = IT_0807 + IT_0808;
    const ccomplex_t IT_0810 = 8*conj(IT_0169);
    const ccomplex_t IT_0811 = (-8)*conj(IT_0167);
    const ccomplex_t IT_0812 = IT_0810 + IT_0811;
    const ccomplex_t IT_0813 = 8*IT_0202;
    const ccomplex_t IT_0814 = (-8)*IT_0204;
    const ccomplex_t IT_0815 = IT_0813 + IT_0814;
    const ccomplex_t IT_0816 = 8*conj(IT_0202);
    const ccomplex_t IT_0817 = (-8)*conj(IT_0204);
    const ccomplex_t IT_0818 = IT_0816 + IT_0817;
    const ccomplex_t IT_0819 = 8*IT_0204;
    const ccomplex_t IT_0820 = (-8)*IT_0202;
    const ccomplex_t IT_0821 = IT_0819 + IT_0820;
    const ccomplex_t IT_0822 = 8*conj(IT_0204);
    const ccomplex_t IT_0823 = (-8)*conj(IT_0202);
    const ccomplex_t IT_0824 = IT_0822 + IT_0823;
    const ccomplex_t IT_0825 = 8*IT_0169;
    const ccomplex_t IT_0826 = (-0.125)*IT_0825;
    const ccomplex_t IT_0827 = IT_0796 + IT_0798;
    const ccomplex_t IT_0828 = -IT_0797;
    const ccomplex_t IT_0829 = IT_0827 + IT_0828;
    const ccomplex_t IT_0830 = IT_0001*IT_0829;
    const ccomplex_t IT_0831 = 0.25*IT_0830;
    const ccomplex_t IT_0832 = s_12 + IT_0795 + IT_0831;
    const ccomplex_t IT_0833 = IT_0015*(IT_0111*IT_0171 + IT_0173*IT_0176 +
       IT_0184*IT_0206 + IT_0207*IT_0210) + IT_0220*(IT_0223*IT_0226 + IT_0227
      *IT_0230 + IT_0232*IT_0235 + IT_0236*IT_0239) + (IT_0184*IT_0226 + IT_0207
      *IT_0230 + IT_0111*IT_0235 + IT_0173*IT_0239)*IT_0251 + (IT_0226*IT_0232 +
       IT_0223*IT_0235 + IT_0230*IT_0236 + IT_0227*IT_0239)*IT_0257 + IT_0269*
      (conj(IT_0182)*IT_0272 + IT_0182*IT_0275 + IT_0110*IT_0277 + IT_0172
      *IT_0279) + IT_0286*(IT_0110*IT_0289 + IT_0172*IT_0292) + IT_0293*(IT_0176
      *IT_0303 + IT_0210*IT_0306 + IT_0313*IT_0314) + IT_0387 + IT_0392*(conj
      (IT_0182)*IT_0394 + IT_0182*IT_0396 + IT_0110*IT_0398 + IT_0172*IT_0400 +
       IT_0401 + IT_0402 + IT_0403 + IT_0404 + IT_0405 + IT_0406) + IT_0313*
      (conj(IT_0222)*IT_0407 + IT_0222*IT_0408 + conj(IT_0222)*IT_0410 + IT_0222
      *IT_0411) + 6*IT_0415*(IT_0223*IT_0236 + IT_0227*(conj(IT_0222) +
       0.166666666666667*IT_0416)) + (s_12 + IT_0421 + IT_0429 + IT_0430)
      *IT_0466 + (s_12 + IT_0472 + IT_0473 + IT_0474)*IT_0491 + (s_12 + IT_0421 
      + IT_0430 + IT_0495)*IT_0499 + IT_0505*IT_0507 + IT_0515*(conj(IT_0182)
      *IT_0516 + IT_0182*IT_0517 + conj(IT_0109)*IT_0518 + IT_0109*IT_0519) +
       IT_0521*(IT_0223*IT_0522 + IT_0227*IT_0523) + IT_0527*(conj(IT_0169)
      *IT_0530 + conj(IT_0204)*IT_0533 + conj(IT_0167)*IT_0534 + conj(IT_0202)
      *IT_0535 + IT_0169*IT_0536 + IT_0204*IT_0537) + IT_0538*IT_0549 + IT_0555
      *IT_0562 + 4*IT_0569*IT_0576 + (conj(IT_0182)*IT_0572 + IT_0182*IT_0573 +
       IT_0110*IT_0574 + IT_0172*IT_0575)*IT_0589 + IT_0595*(IT_0110*IT_0596 +
       IT_0109*IT_0597 + conj(IT_0109)*IT_0598 + conj(IT_0182)*IT_0599) +
       IT_0637 + IT_0645 + IT_0649*(IT_0650 + IT_0651 + IT_0652 + IT_0653) + 
      (IT_0110*IT_0394 + IT_0172*IT_0396 + conj(IT_0182)*IT_0398 + IT_0182
      *IT_0400 + IT_0401 + IT_0402 + IT_0403 + IT_0404 + IT_0405 + IT_0406 +
       IT_0650 + IT_0651 + IT_0652 + IT_0653)*IT_0661 + IT_0669*IT_0680 + 
      (IT_0232*IT_0353 + 2*conj(IT_0222)*IT_0354 + IT_0222*IT_0356 + IT_0223
      *IT_0368 + IT_0222*IT_0370)*IT_0682 + conj(IT_0293)*(IT_0171*IT_0303 +
       IT_0206*IT_0306 + IT_0313*IT_0373 + IT_0684*IT_0685) + conj(IT_0180)*
      (IT_0671*IT_0693 + IT_0394*IT_0698 + IT_0620*IT_0702 + IT_0398*IT_0707 +
       IT_0712*IT_0713 + IT_0609*IT_0716 + IT_0719*IT_0720 + IT_0685*IT_0728) +
       IT_0180*(IT_0673*IT_0693 + IT_0396*IT_0698 + IT_0624*IT_0702 + IT_0400
      *IT_0707 + IT_0613*IT_0716 + IT_0712*IT_0729 + IT_0728*IT_0730) + (-4)
      *IT_0183*IT_0731 + IT_0738*IT_0745 + IT_0753*(IT_0407*IT_0461 + IT_0408
      *IT_0485 + conj(IT_0447)*IT_0685 + IT_0447*IT_0730 + conj(IT_0409)*IT_0741
       + IT_0409*IT_0744 + conj(IT_0442)*IT_0754 + IT_0442*IT_0755) + (IT_0408
      *IT_0442 + IT_0407*conj(IT_0442) + IT_0411*IT_0447 + IT_0410*conj(IT_0447)
       + IT_0539*IT_0685 + IT_0544*IT_0730 + IT_0461*IT_0754 + IT_0485*IT_0755)
      *IT_0761 + IT_0762*(IT_0110*IT_0763 + IT_0172*IT_0764) + IT_0765*IT_0772 +
       conj(IT_0102)*IT_0775 + (-4)*IT_0095*IT_0777 + (-8)*IT_0102*IT_0779 + 6
      *conj(IT_0409)*(0.166666666666667*IT_0206*IT_0303 + -IT_0236*IT_0521 +
       0.166666666666667*IT_0684*IT_0754 + IT_0293*IT_0782 + -IT_0171*IT_0783 + 
      -IT_0674*IT_0784) + 6*IT_0409*(0.166666666666667*IT_0210*IT_0303 + 
      -IT_0232*IT_0521 + conj(IT_0293)*IT_0782 + -IT_0176*IT_0783 + -IT_0676
      *IT_0784) + IT_0792 + (-8)*conj(IT_0095)*IT_0793 + IT_0803*(IT_0539
      *IT_0806 + IT_0544*IT_0809 + IT_0442*IT_0812 + IT_0461*IT_0815 + IT_0485
      *IT_0818 + conj(IT_0447)*IT_0821 + IT_0447*IT_0824 + (-8)*conj(IT_0442)*
      (IT_0167 + IT_0826)) + (conj(IT_0384)*IT_0806 + IT_0384*IT_0809 + IT_0367
      *IT_0812 + IT_0372*IT_0815 + IT_0374*IT_0818 + conj(IT_0352)*IT_0821 +
       IT_0352*IT_0824 + (-8)*conj(IT_0367)*(IT_0167 + IT_0826))*IT_0832;
    return create_ccomplex_return(IT_0833);
}

