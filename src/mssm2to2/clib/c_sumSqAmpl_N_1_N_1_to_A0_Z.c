#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_A0_Z(
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
    const ccomplex_t IT_0000 = s_14*s_24;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_12*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = -IT_0012 + IT_0016 + IT_0021 + -IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0017;
    const ccomplex_t IT_0027 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = IT_0008*IT_0018;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0029*IT_0031;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0028 + IT_0030 + 
      -IT_0032 + -IT_0033);
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = pow(m_A0, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0037 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0043 = IT_0010*IT_0042;
    const ccomplex_t IT_0044 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0045 = IT_0014*IT_0044;
    const ccomplex_t IT_0046 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0047 = IT_0019*IT_0046;
    const ccomplex_t IT_0048 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0049 = IT_0022*IT_0048;
    const ccomplex_t IT_0050 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0051 = IT_0010*IT_0050;
    const ccomplex_t IT_0052 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0053 = IT_0014*IT_0052;
    const ccomplex_t IT_0054 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0055 = IT_0019*IT_0054;
    const ccomplex_t IT_0056 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0057 = IT_0022*IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0043 + 0.5*IT_0045 + 0.5*IT_0047 + (
      -0.5)*IT_0049 + (-0.5)*IT_0051 + 0.5*IT_0053 + 0.5*IT_0055 + (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0060 = IT_0026*IT_0059;
    const ccomplex_t IT_0061 = IT_0029*IT_0059;
    const ccomplex_t IT_0062 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0063 = IT_0026*IT_0062;
    const ccomplex_t IT_0064 = IT_0029*IT_0062;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0060 + IT_0061 + 
      -IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = pow(m_N_1, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_13 + IT_0037 + IT_0068 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0067*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0034;
    const ccomplex_t IT_0073 = IT_0025*IT_0072;
    const ccomplex_t IT_0074 = IT_0039*IT_0073;
    const ccomplex_t IT_0075 = -IT_0071 + (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0079 = IT_0014*IT_0078;
    const ccomplex_t IT_0080 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0081 = IT_0019*IT_0080;
    const ccomplex_t IT_0082 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0083 = IT_0022*IT_0082;
    const ccomplex_t IT_0084 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0085 = IT_0010*IT_0084;
    const ccomplex_t IT_0086 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0087 = IT_0014*IT_0086;
    const ccomplex_t IT_0088 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0089 = IT_0019*IT_0088;
    const ccomplex_t IT_0090 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0091 = IT_0022*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0077 + 0.5*IT_0079 + 0.5*IT_0081 + (
      -0.5)*IT_0083 + (-0.5)*IT_0085 + 0.5*IT_0087 + 0.5*IT_0089 + (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0094 = IT_0026*IT_0093;
    const ccomplex_t IT_0095 = IT_0029*IT_0093;
    const ccomplex_t IT_0096 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0097 = IT_0026*IT_0096;
    const ccomplex_t IT_0098 = IT_0029*IT_0096;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0094 + IT_0095 + 
      -IT_0097 + -IT_0098);
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0092*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_13 + IT_0037 + IT_0068 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0106 = IT_0010*IT_0105;
    const ccomplex_t IT_0107 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0108 = IT_0014*IT_0107;
    const ccomplex_t IT_0109 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0110 = IT_0019*IT_0109;
    const ccomplex_t IT_0111 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0112 = IT_0022*IT_0111;
    const ccomplex_t IT_0113 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0114 = IT_0010*IT_0113;
    const ccomplex_t IT_0115 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0116 = IT_0014*IT_0115;
    const ccomplex_t IT_0117 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0118 = IT_0019*IT_0117;
    const ccomplex_t IT_0119 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0120 = IT_0022*IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0106 + 0.5*IT_0108 + 0.5*IT_0110 + (
      -0.5)*IT_0112 + (-0.5)*IT_0114 + 0.5*IT_0116 + 0.5*IT_0118 + (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0123 = IT_0026*IT_0122;
    const ccomplex_t IT_0124 = IT_0029*IT_0122;
    const ccomplex_t IT_0125 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0126 = IT_0026*IT_0125;
    const ccomplex_t IT_0127 = IT_0029*IT_0125;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0123 + IT_0124 + 
      -IT_0126 + -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = IT_0121*IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0037 + IT_0068 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = -IT_0104 + -IT_0133;
    const ccomplex_t IT_0135 = conj(IT_0075) + conj(IT_0134);
    const ccomplex_t IT_0136 = -conj(IT_0041) + IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_23 + IT_0037 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = -IT_0036;
    const ccomplex_t IT_0140 = IT_0073 + IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_23 + IT_0037 + IT_0068 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0101*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_23 + IT_0037 + IT_0068 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0067*IT_0147;
    const ccomplex_t IT_0149 = cpow((-2)*s_23 + IT_0037 + IT_0068 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = IT_0130*IT_0150;
    const ccomplex_t IT_0152 = IT_0145 + IT_0148 + IT_0151;
    const ccomplex_t IT_0153 = IT_0142 + -IT_0152;
    const ccomplex_t IT_0154 = IT_0075 + IT_0134;
    const ccomplex_t IT_0155 = -IT_0041;
    const ccomplex_t IT_0156 = IT_0154 + IT_0155;
    const ccomplex_t IT_0157 = -conj(IT_0152);
    const ccomplex_t IT_0158 = conj(IT_0142) + IT_0157;
    const ccomplex_t IT_0159 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0160 = IT_0010*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0162 = IT_0014*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0164 = IT_0019*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0166 = IT_0022*IT_0165;
    const ccomplex_t IT_0167 = IT_0160 + -IT_0162 + -IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = IT_0072*IT_0167;
    const ccomplex_t IT_0169 = IT_0039*IT_0168;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0035*IT_0167;
    const ccomplex_t IT_0172 = IT_0039*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0174 = IT_0010*IT_0173;
    const ccomplex_t IT_0175 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0176 = IT_0014*IT_0175;
    const ccomplex_t IT_0177 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0178 = IT_0019*IT_0177;
    const ccomplex_t IT_0179 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0180 = IT_0022*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0182 = IT_0010*IT_0181;
    const ccomplex_t IT_0183 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0184 = IT_0014*IT_0183;
    const ccomplex_t IT_0185 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0186 = IT_0019*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0188 = IT_0022*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0174 + (-0.5)*IT_0176 + (-0.5)*IT_0178 +
       0.5*IT_0180 + 0.5*IT_0182 + (-0.5)*IT_0184 + (-0.5)*IT_0186 + 0.5*IT_0188;
    const ccomplex_t IT_0190 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0191 = IT_0026*IT_0190;
    const ccomplex_t IT_0192 = IT_0029*IT_0190;
    const ccomplex_t IT_0193 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0194 = IT_0026*IT_0193;
    const ccomplex_t IT_0195 = IT_0029*IT_0193;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0191 + IT_0192 + 
      -IT_0194 + -IT_0195);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0189*IT_0197;
    const ccomplex_t IT_0199 = IT_0070*IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0172 + -IT_0199;
    const ccomplex_t IT_0201 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0202 = IT_0010*IT_0201;
    const ccomplex_t IT_0203 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0204 = IT_0014*IT_0203;
    const ccomplex_t IT_0205 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0206 = IT_0019*IT_0205;
    const ccomplex_t IT_0207 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0208 = IT_0022*IT_0207;
    const ccomplex_t IT_0209 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0210 = IT_0010*IT_0209;
    const ccomplex_t IT_0211 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0212 = IT_0014*IT_0211;
    const ccomplex_t IT_0213 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0214 = IT_0019*IT_0213;
    const ccomplex_t IT_0215 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0216 = IT_0022*IT_0215;
    const ccomplex_t IT_0217 = 0.5*IT_0202 + (-0.5)*IT_0204 + (-0.5)*IT_0206 +
       0.5*IT_0208 + 0.5*IT_0210 + (-0.5)*IT_0212 + (-0.5)*IT_0214 + 0.5*IT_0216;
    const ccomplex_t IT_0218 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0219 = IT_0026*IT_0218;
    const ccomplex_t IT_0220 = IT_0029*IT_0218;
    const ccomplex_t IT_0221 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0222 = IT_0026*IT_0221;
    const ccomplex_t IT_0223 = IT_0029*IT_0221;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*(IT_0219 + IT_0220 + 
      -IT_0222 + -IT_0223);
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = IT_0217*IT_0225;
    const ccomplex_t IT_0227 = IT_0132*IT_0226;
    const ccomplex_t IT_0228 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0229 = IT_0010*IT_0228;
    const ccomplex_t IT_0230 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0231 = IT_0014*IT_0230;
    const ccomplex_t IT_0232 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0233 = IT_0019*IT_0232;
    const ccomplex_t IT_0234 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0235 = IT_0022*IT_0234;
    const ccomplex_t IT_0236 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0237 = IT_0010*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0239 = IT_0014*IT_0238;
    const ccomplex_t IT_0240 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0241 = IT_0019*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0243 = IT_0022*IT_0242;
    const ccomplex_t IT_0244 = 0.5*IT_0229 + (-0.5)*IT_0231 + (-0.5)*IT_0233 +
       0.5*IT_0235 + 0.5*IT_0237 + (-0.5)*IT_0239 + (-0.5)*IT_0241 + 0.5*IT_0243;
    const ccomplex_t IT_0245 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0246 = IT_0026*IT_0245;
    const ccomplex_t IT_0247 = IT_0029*IT_0245;
    const ccomplex_t IT_0248 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0249 = IT_0026*IT_0248;
    const ccomplex_t IT_0250 = IT_0029*IT_0248;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*(IT_0246 + IT_0247 + 
      -IT_0249 + -IT_0250);
    const ccomplex_t IT_0252 = (-0.5)*IT_0251;
    const ccomplex_t IT_0253 = IT_0244*IT_0252;
    const ccomplex_t IT_0254 = IT_0103*IT_0253;
    const ccomplex_t IT_0255 = -IT_0227 + -IT_0254;
    const ccomplex_t IT_0256 = conj(IT_0200) + conj(IT_0255);
    const ccomplex_t IT_0257 = -IT_0168;
    const ccomplex_t IT_0258 = IT_0171 + IT_0257;
    const ccomplex_t IT_0259 = IT_0138*IT_0258;
    const ccomplex_t IT_0260 = (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = IT_0150*IT_0226;
    const ccomplex_t IT_0262 = IT_0147*IT_0198;
    const ccomplex_t IT_0263 = IT_0144*IT_0253;
    const ccomplex_t IT_0264 = IT_0261 + IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = IT_0260 + -IT_0264;
    const ccomplex_t IT_0266 = IT_0200 + IT_0255;
    const ccomplex_t IT_0267 = -IT_0170;
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = conj(IT_0260) + -conj(IT_0264);
    const ccomplex_t IT_0270 = pow(s_14, 2);
    const ccomplex_t IT_0271 = IT_0001*IT_0270;
    const ccomplex_t IT_0272 = -IT_0068;
    const ccomplex_t IT_0273 = IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = IT_0068*IT_0273;
    const ccomplex_t IT_0275 = cos(alpha);
    const ccomplex_t IT_0276 = IT_0009*IT_0275;
    const ccomplex_t IT_0277 = IT_0161*IT_0276;
    const ccomplex_t IT_0278 = sin(alpha);
    const ccomplex_t IT_0279 = IT_0009*IT_0278;
    const ccomplex_t IT_0280 = IT_0159*IT_0279;
    const ccomplex_t IT_0281 = IT_0018*IT_0275;
    const ccomplex_t IT_0282 = IT_0165*IT_0281;
    const ccomplex_t IT_0283 = IT_0018*IT_0278;
    const ccomplex_t IT_0284 = IT_0163*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*(IT_0277 + -IT_0280 + 
      -IT_0282 + IT_0284);
    const ccomplex_t IT_0286 = e_em*IT_0009*IT_0018*(IT_0013*IT_0275 + 
      -IT_0007*IT_0278);
    const ccomplex_t IT_0287 = IT_0285*IT_0286;
    const ccomplex_t IT_0288 = cpow((-2)*s_12 + (-2)*IT_0068 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0287*IT_0289;
    const ccomplex_t IT_0291 = IT_0159*IT_0276;
    const ccomplex_t IT_0292 = IT_0161*IT_0279;
    const ccomplex_t IT_0293 = IT_0163*IT_0281;
    const ccomplex_t IT_0294 = IT_0165*IT_0283;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*(IT_0291 + IT_0292 + 
      -IT_0293 + -IT_0294);
    const ccomplex_t IT_0296 = -IT_0295;
    const ccomplex_t IT_0297 = e_em*IT_0009*IT_0018*(IT_0007*IT_0275 + IT_0013
      *IT_0278);
    const ccomplex_t IT_0298 = -IT_0297;
    const ccomplex_t IT_0299 = IT_0296*IT_0298;
    const ccomplex_t IT_0300 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0068 + -reg_prop, -1);
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = IT_0299*IT_0301;
    const ccomplex_t IT_0303 = (-2)*IT_0261 + (-2)*IT_0262 + (-2)*IT_0263 +
       0.5*IT_0290 + 0.5*IT_0302;
    const ccomplex_t IT_0304 = IT_0073*IT_0138;
    const ccomplex_t IT_0305 = IT_0015*IT_0276;
    const ccomplex_t IT_0306 = IT_0011*IT_0279;
    const ccomplex_t IT_0307 = IT_0023*IT_0281;
    const ccomplex_t IT_0308 = IT_0020*IT_0283;
    const ccomplex_t IT_0309 = (0 + _Complex_I*1)*(IT_0305 + -IT_0306 + 
      -IT_0307 + IT_0308);
    const ccomplex_t IT_0310 = IT_0286*IT_0309;
    const ccomplex_t IT_0311 = IT_0289*IT_0310;
    const ccomplex_t IT_0312 = IT_0011*IT_0276;
    const ccomplex_t IT_0313 = IT_0015*IT_0279;
    const ccomplex_t IT_0314 = IT_0020*IT_0281;
    const ccomplex_t IT_0315 = IT_0023*IT_0283;
    const ccomplex_t IT_0316 = (0 + _Complex_I*1)*(IT_0312 + IT_0313 + 
      -IT_0314 + -IT_0315);
    const ccomplex_t IT_0317 = -IT_0316;
    const ccomplex_t IT_0318 = IT_0298*IT_0317;
    const ccomplex_t IT_0319 = IT_0301*IT_0318;
    const ccomplex_t IT_0320 = IT_0036*IT_0138;
    const ccomplex_t IT_0321 = IT_0304 + (-0.5)*IT_0311 + (-0.5)*IT_0319 + 
      -IT_0320;
    const ccomplex_t IT_0322 = 2*IT_0321;
    const ccomplex_t IT_0323 = (-2)*IT_0145 + (-2)*IT_0148 + (-2)*IT_0151 +
       0.5*IT_0311 + 0.5*IT_0319;
    const ccomplex_t IT_0324 = (-2)*IT_0323;
    const ccomplex_t IT_0325 = IT_0322 + IT_0324;
    const ccomplex_t IT_0326 = 2*conj(IT_0321);
    const ccomplex_t IT_0327 = (-2)*conj(IT_0323);
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = IT_0138*IT_0168;
    const ccomplex_t IT_0330 = IT_0138*IT_0171;
    const ccomplex_t IT_0331 = (-0.5)*IT_0290 + (-0.5)*IT_0302 + -IT_0329 +
       IT_0330;
    const ccomplex_t IT_0332 = 2*IT_0323;
    const ccomplex_t IT_0333 = (-2)*IT_0321;
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = 2*conj(IT_0323);
    const ccomplex_t IT_0336 = (-2)*conj(IT_0321);
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = s_12 + IT_0339;
    const ccomplex_t IT_0341 = s_12*IT_0340;
    const ccomplex_t IT_0342 = IT_0040 + 0.5*IT_0311 + 0.5*IT_0319;
    const ccomplex_t IT_0343 = 2*IT_0071 + IT_0074 + 2*IT_0104 + 2*IT_0133 + (
      -0.5)*IT_0311 + (-0.5)*IT_0319;
    const ccomplex_t IT_0344 = IT_0169 + 0.5*IT_0290 + 0.5*IT_0302;
    const ccomplex_t IT_0345 = 2*IT_0331;
    const ccomplex_t IT_0346 = (-2)*IT_0303;
    const ccomplex_t IT_0347 = IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = 2*conj(IT_0331);
    const ccomplex_t IT_0349 = (-2)*conj(IT_0303);
    const ccomplex_t IT_0350 = IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = IT_0172 + 2*IT_0199 + 2*IT_0227 + 2*IT_0254 + (
      -0.5)*IT_0290 + (-0.5)*IT_0302;
    const ccomplex_t IT_0352 = 2*IT_0303;
    const ccomplex_t IT_0353 = (-2)*IT_0331;
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = 2*conj(IT_0303);
    const ccomplex_t IT_0356 = (-2)*conj(IT_0331);
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = s_12*s_14;
    const ccomplex_t IT_0359 = s_24*IT_0068;
    const ccomplex_t IT_0360 = -IT_0359;
    const ccomplex_t IT_0361 = IT_0358 + IT_0360;
    const ccomplex_t IT_0362 = conj(IT_0041) + conj(IT_0142);
    const ccomplex_t IT_0363 = IT_0041 + IT_0142;
    const ccomplex_t IT_0364 = conj(IT_0170) + conj(IT_0260);
    const ccomplex_t IT_0365 = IT_0170 + IT_0260;
    const ccomplex_t IT_0366 = IT_0361*((IT_0135 + conj(IT_0152))*IT_0334 + 
      (IT_0152 + IT_0154)*IT_0337 + (IT_0256 + conj(IT_0264))*IT_0354 + (IT_0264
       + IT_0266)*IT_0357 + IT_0325*IT_0362 + IT_0328*IT_0363 + IT_0347*IT_0364 
      + IT_0350*IT_0365);
    const ccomplex_t IT_0367 = IT_0068*IT_0340;
    const ccomplex_t IT_0368 = 2*IT_0343;
    const ccomplex_t IT_0369 = (-2)*IT_0342;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 2*conj(IT_0343);
    const ccomplex_t IT_0372 = (-2)*conj(IT_0342);
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = 2*IT_0342;
    const ccomplex_t IT_0375 = (-2)*IT_0343;
    const ccomplex_t IT_0376 = IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = 2*conj(IT_0342);
    const ccomplex_t IT_0378 = (-2)*conj(IT_0343);
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = pow(s_24, 2);
    const ccomplex_t IT_0381 = IT_0001*IT_0380;
    const ccomplex_t IT_0382 = -IT_0381;
    const ccomplex_t IT_0383 = IT_0068 + IT_0382;
    const ccomplex_t IT_0384 = IT_0068*IT_0383;
    const ccomplex_t IT_0385 = -IT_0384;
    const ccomplex_t IT_0386 = m_N_1*IT_0383;
    const ccomplex_t IT_0387 = -IT_0386;
    const ccomplex_t IT_0388 = IT_0100*IT_0244;
    const ccomplex_t IT_0389 = m_N_2*IT_0388;
    const ccomplex_t IT_0390 = m_N_1*IT_0253;
    const ccomplex_t IT_0391 = -IT_0390;
    const ccomplex_t IT_0392 = IT_0389 + IT_0391;
    const ccomplex_t IT_0393 = IT_0103*IT_0392;
    const ccomplex_t IT_0394 = m_N_1*IT_0073;
    const ccomplex_t IT_0395 = IT_0138*IT_0394;
    const ccomplex_t IT_0396 = m_N_1*IT_0036;
    const ccomplex_t IT_0397 = -IT_0396;
    const ccomplex_t IT_0398 = IT_0394 + IT_0397;
    const ccomplex_t IT_0399 = IT_0138*IT_0398;
    const ccomplex_t IT_0400 = m_N_1*IT_0171;
    const ccomplex_t IT_0401 = m_N_1*IT_0168;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = IT_0039*IT_0403;
    const ccomplex_t IT_0405 = IT_0066*IT_0189;
    const ccomplex_t IT_0406 = m_N_3*IT_0405;
    const ccomplex_t IT_0407 = m_N_1*IT_0198;
    const ccomplex_t IT_0408 = -IT_0407;
    const ccomplex_t IT_0409 = IT_0406 + IT_0408;
    const ccomplex_t IT_0410 = IT_0070*IT_0409;
    const ccomplex_t IT_0411 = IT_0129*IT_0217;
    const ccomplex_t IT_0412 = m_N_4*IT_0411;
    const ccomplex_t IT_0413 = m_N_1*IT_0226;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0412 + IT_0414;
    const ccomplex_t IT_0416 = IT_0132*IT_0415;
    const ccomplex_t IT_0417 = IT_0393 + 0.5*IT_0395 + 0.5*IT_0399 + (-0.5)
      *IT_0404 + IT_0410 + IT_0416;
    const ccomplex_t IT_0418 = IT_0138*IT_0396;
    const ccomplex_t IT_0419 = IT_0092*IT_0252;
    const ccomplex_t IT_0420 = m_N_2*IT_0419;
    const ccomplex_t IT_0421 = m_N_1*IT_0101;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0420 + IT_0422;
    const ccomplex_t IT_0424 = IT_0144*IT_0423;
    const ccomplex_t IT_0425 = IT_0058*IT_0197;
    const ccomplex_t IT_0426 = m_N_3*IT_0425;
    const ccomplex_t IT_0427 = m_N_1*IT_0067;
    const ccomplex_t IT_0428 = -IT_0427;
    const ccomplex_t IT_0429 = IT_0426 + IT_0428;
    const ccomplex_t IT_0430 = IT_0147*IT_0429;
    const ccomplex_t IT_0431 = m_N_1*IT_0130;
    const ccomplex_t IT_0432 = -IT_0431;
    const ccomplex_t IT_0433 = IT_0121*IT_0225;
    const ccomplex_t IT_0434 = m_N_4*IT_0433;
    const ccomplex_t IT_0435 = IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = IT_0150*IT_0435;
    const ccomplex_t IT_0437 = (-0.5)*IT_0404 + (-0.5)*IT_0418 + -IT_0424 + 
      -IT_0430 + -IT_0436;
    const ccomplex_t IT_0438 = conj(IT_0417) + conj(IT_0437);
    const ccomplex_t IT_0439 = IT_0417 + IT_0437;
    const ccomplex_t IT_0440 = IT_0070*IT_0429;
    const ccomplex_t IT_0441 = IT_0138*IT_0400;
    const ccomplex_t IT_0442 = IT_0138*IT_0403;
    const ccomplex_t IT_0443 = IT_0039*IT_0398;
    const ccomplex_t IT_0444 = IT_0103*IT_0423;
    const ccomplex_t IT_0445 = IT_0132*IT_0435;
    const ccomplex_t IT_0446 = IT_0440 + 0.5*IT_0441 + 0.5*IT_0442 + (-0.5)
      *IT_0443 + IT_0444 + IT_0445;
    const ccomplex_t IT_0447 = IT_0138*IT_0401;
    const ccomplex_t IT_0448 = IT_0144*IT_0392;
    const ccomplex_t IT_0449 = IT_0147*IT_0409;
    const ccomplex_t IT_0450 = IT_0150*IT_0415;
    const ccomplex_t IT_0451 = (-0.5)*IT_0443 + (-0.5)*IT_0447 + -IT_0448 + 
      -IT_0449 + -IT_0450;
    const ccomplex_t IT_0452 = conj(IT_0446) + conj(IT_0451);
    const ccomplex_t IT_0453 = 2*IT_0344;
    const ccomplex_t IT_0454 = (-2)*IT_0351;
    const ccomplex_t IT_0455 = IT_0453 + IT_0454;
    const ccomplex_t IT_0456 = IT_0446 + IT_0451;
    const ccomplex_t IT_0457 = 2*conj(IT_0344);
    const ccomplex_t IT_0458 = (-2)*conj(IT_0351);
    const ccomplex_t IT_0459 = IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = pow(m_Z, 2);
    const ccomplex_t IT_0461 = s_24*IT_0460;
    const ccomplex_t IT_0462 = IT_0001*IT_0461;
    const ccomplex_t IT_0463 = (-0.25)*IT_0462;
    const ccomplex_t IT_0464 = s_24 + IT_0463;
    const ccomplex_t IT_0465 = m_N_1*IT_0464;
    const ccomplex_t IT_0466 = (-4)*IT_0465;
    const ccomplex_t IT_0467 = 2*IT_0417;
    const ccomplex_t IT_0468 = 2*IT_0437;
    const ccomplex_t IT_0469 = IT_0467 + IT_0468;
    const ccomplex_t IT_0470 = 2*conj(IT_0417);
    const ccomplex_t IT_0471 = 2*conj(IT_0437);
    const ccomplex_t IT_0472 = IT_0470 + IT_0471;
    const ccomplex_t IT_0473 = 2*IT_0446;
    const ccomplex_t IT_0474 = 2*IT_0451;
    const ccomplex_t IT_0475 = IT_0473 + IT_0474;
    const ccomplex_t IT_0476 = 2*conj(IT_0446);
    const ccomplex_t IT_0477 = 2*conj(IT_0451);
    const ccomplex_t IT_0478 = IT_0476 + IT_0477;
    const ccomplex_t IT_0479 = (-2)*IT_0417;
    const ccomplex_t IT_0480 = (-2)*IT_0437;
    const ccomplex_t IT_0481 = IT_0479 + IT_0480;
    const ccomplex_t IT_0482 = (-2)*conj(IT_0417);
    const ccomplex_t IT_0483 = (-2)*conj(IT_0437);
    const ccomplex_t IT_0484 = IT_0482 + IT_0483;
    const ccomplex_t IT_0485 = (-2)*IT_0446;
    const ccomplex_t IT_0486 = (-2)*IT_0451;
    const ccomplex_t IT_0487 = IT_0485 + IT_0486;
    const ccomplex_t IT_0488 = (-2)*conj(IT_0446);
    const ccomplex_t IT_0489 = (-2)*conj(IT_0451);
    const ccomplex_t IT_0490 = IT_0488 + IT_0489;
    const ccomplex_t IT_0491 = s_14*IT_0460;
    const ccomplex_t IT_0492 = IT_0001*IT_0491;
    const ccomplex_t IT_0493 = 0.5*IT_0492;
    const ccomplex_t IT_0494 = s_14 + IT_0493;
    const ccomplex_t IT_0495 = m_N_1*IT_0494;
    const ccomplex_t IT_0496 = 2*IT_0495;
    const ccomplex_t IT_0497 = (-0.25)*IT_0492;
    const ccomplex_t IT_0498 = s_14 + IT_0497;
    const ccomplex_t IT_0499 = m_N_1*IT_0498;
    const ccomplex_t IT_0500 = (-4)*IT_0499;
    const ccomplex_t IT_0501 = 0.5*IT_0462;
    const ccomplex_t IT_0502 = s_24 + IT_0501;
    const ccomplex_t IT_0503 = m_N_1*IT_0502;
    const ccomplex_t IT_0504 = 2*IT_0503;
    const ccomplex_t IT_0505 = s_12*IT_0460;
    const ccomplex_t IT_0506 = (-2)*IT_0000;
    const ccomplex_t IT_0507 = IT_0505 + IT_0506;
    const ccomplex_t IT_0508 = IT_0001*IT_0507;
    const ccomplex_t IT_0509 = (-0.5)*IT_0508;
    const ccomplex_t IT_0510 = s_12 + IT_0509;
    const ccomplex_t IT_0511 = 4*IT_0446;
    const ccomplex_t IT_0512 = 6*conj(IT_0417);
    const ccomplex_t IT_0513 = m_N_1*IT_0273;
    const ccomplex_t IT_0514 = s_12*IT_0273;
    const ccomplex_t IT_0515 = m_N_1*IT_0340;
    const ccomplex_t IT_0516 = IT_0328*IT_0417 + IT_0325*conj(IT_0417) + 
      (IT_0345 + IT_0346 + IT_0374 + IT_0375)*IT_0452 + conj(IT_0437)*(IT_0322 +
       IT_0324 + IT_0453 + IT_0454) + (IT_0348 + IT_0349 + IT_0377 + IT_0378)
      *IT_0456 + IT_0437*(IT_0326 + IT_0327 + IT_0457 + IT_0458) + conj(IT_0344)
      *IT_0467 + IT_0344*IT_0470 + conj(IT_0351)*IT_0479 + IT_0351*IT_0482;
    const ccomplex_t IT_0517 = s_12*IT_0383;
    const ccomplex_t IT_0518 = -IT_0517;
    const ccomplex_t IT_0519 = 2*IT_0351;
    const ccomplex_t IT_0520 = s_12*s_24;
    const ccomplex_t IT_0521 = s_14*IT_0068;
    const ccomplex_t IT_0522 = -IT_0521;
    const ccomplex_t IT_0523 = IT_0520 + IT_0522;
    const ccomplex_t IT_0524 = (-2)*IT_0344;
    const ccomplex_t IT_0525 = 2*conj(IT_0351);
    const ccomplex_t IT_0526 = (-2)*conj(IT_0344);
    const ccomplex_t IT_0527 = (IT_0135 + conj(IT_0152))*IT_0370 + (IT_0152 +
       IT_0154)*IT_0373 + IT_0362*IT_0376 + IT_0363*IT_0379 + IT_0364*IT_0455 +
       IT_0365*IT_0459 + (IT_0256 + conj(IT_0264))*(IT_0519 + IT_0524) + 
      (IT_0264 + IT_0266)*(IT_0525 + IT_0526);
    const ccomplex_t IT_0528 = IT_0068*IT_0460;
    const ccomplex_t IT_0529 = (-3)*IT_0528;
    const ccomplex_t IT_0530 = 2*IT_0075;
    const ccomplex_t IT_0531 = 2*IT_0134;
    const ccomplex_t IT_0532 = (-2)*IT_0041;
    const ccomplex_t IT_0533 = IT_0530 + IT_0531 + IT_0532;
    const ccomplex_t IT_0534 = 2*conj(IT_0075);
    const ccomplex_t IT_0535 = 2*conj(IT_0134);
    const ccomplex_t IT_0536 = (-2)*conj(IT_0041);
    const ccomplex_t IT_0537 = IT_0534 + IT_0535 + IT_0536;
    const ccomplex_t IT_0538 = 2*IT_0142;
    const ccomplex_t IT_0539 = (-2)*IT_0152;
    const ccomplex_t IT_0540 = IT_0538 + IT_0539;
    const ccomplex_t IT_0541 = 2*conj(IT_0142);
    const ccomplex_t IT_0542 = (-2)*conj(IT_0152);
    const ccomplex_t IT_0543 = IT_0541 + IT_0542;
    const ccomplex_t IT_0544 = 2*IT_0200;
    const ccomplex_t IT_0545 = 2*IT_0255;
    const ccomplex_t IT_0546 = IT_0544 + IT_0545;
    const ccomplex_t IT_0547 = 2*conj(IT_0200);
    const ccomplex_t IT_0548 = 2*conj(IT_0255);
    const ccomplex_t IT_0549 = IT_0547 + IT_0548;
    const ccomplex_t IT_0550 = (-2)*IT_0075;
    const ccomplex_t IT_0551 = (-2)*IT_0134;
    const ccomplex_t IT_0552 = IT_0550 + IT_0551;
    const ccomplex_t IT_0553 = (-2)*conj(IT_0075);
    const ccomplex_t IT_0554 = (-2)*conj(IT_0134);
    const ccomplex_t IT_0555 = IT_0553 + IT_0554;
    const ccomplex_t IT_0556 = 2*IT_0260;
    const ccomplex_t IT_0557 = 2*conj(IT_0260);
    const ccomplex_t IT_0558 = (-2)*IT_0142;
    const ccomplex_t IT_0559 = (-2)*conj(IT_0142);
    const ccomplex_t IT_0560 = conj(IT_0170)*IT_0533 + IT_0170*IT_0537 + conj
      (IT_0264)*IT_0540 + IT_0264*IT_0543 + conj(IT_0041)*IT_0546 + IT_0041
      *IT_0549 + IT_0256*IT_0552 + IT_0266*IT_0555 + conj(IT_0152)*IT_0556 +
       IT_0152*IT_0557 + conj(IT_0260)*IT_0558 + IT_0260*IT_0559;
    const ccomplex_t IT_0561 = 3*IT_0528;
    const ccomplex_t IT_0562 = conj(IT_0264)*IT_0533 + IT_0264*IT_0537 + conj
      (IT_0170)*IT_0540 + IT_0170*IT_0543 + conj(IT_0152)*IT_0546 + IT_0152
      *IT_0549 + conj(IT_0260)*IT_0552 + IT_0260*IT_0555 + conj(IT_0041)*IT_0556
       + IT_0041*IT_0557 + IT_0256*IT_0558 + IT_0266*IT_0559;
    const ccomplex_t IT_0563 = (-0.5)*IT_0004;
    const ccomplex_t IT_0564 = IT_0505 + IT_0506 + IT_0563;
    const ccomplex_t IT_0565 = (-4)*IT_0200;
    const ccomplex_t IT_0566 = (-8)*IT_0006*(IT_0136*IT_0153 + IT_0156*IT_0158
       + -(conj(IT_0170) + -IT_0256)*IT_0265 + IT_0268*IT_0269) + IT_0274*(conj
      (IT_0303)*IT_0325 + IT_0303*IT_0328 + conj(IT_0331)*IT_0334 + IT_0331
      *IT_0337) + IT_0341*(IT_0328*IT_0342 + IT_0325*conj(IT_0342) + IT_0337
      *IT_0343 + IT_0334*conj(IT_0343) + conj(IT_0344)*IT_0347 + IT_0344*IT_0350
       + conj(IT_0351)*IT_0354 + IT_0351*IT_0357) + IT_0366 + IT_0367*(IT_0337
      *IT_0344 + IT_0334*conj(IT_0344) + IT_0328*IT_0351 + IT_0325*conj(IT_0351)
       + conj(IT_0331)*IT_0370 + IT_0331*IT_0373 + conj(IT_0303)*IT_0376 +
       IT_0303*IT_0379) + (conj(IT_0344)*IT_0370 + IT_0344*IT_0373 + conj
      (IT_0351)*IT_0376 + IT_0351*IT_0379)*IT_0385 + IT_0387*(IT_0376*IT_0438 +
       IT_0379*IT_0439 + IT_0452*IT_0455 + IT_0456*IT_0459) + IT_0466*(conj
      (IT_0152)*IT_0469 + IT_0152*IT_0472 + conj(IT_0264)*IT_0475 + IT_0264
      *IT_0478 + conj(IT_0142)*IT_0481 + IT_0142*IT_0484 + conj(IT_0260)*IT_0487
       + IT_0260*IT_0490) + (conj(IT_0260)*IT_0469 + IT_0260*IT_0472 + conj
      (IT_0142)*IT_0475 + IT_0142*IT_0478 + conj(IT_0264)*IT_0481 + IT_0264
      *IT_0484 + conj(IT_0152)*IT_0487 + IT_0152*IT_0490)*IT_0496 + (IT_0256
      *IT_0469 + IT_0266*IT_0472 + IT_0135*IT_0475 + IT_0154*IT_0478 + conj
      (IT_0170)*IT_0481 + IT_0170*IT_0484 + conj(IT_0041)*IT_0487 + IT_0041
      *IT_0490)*IT_0500 + (conj(IT_0041)*IT_0469 + IT_0041*IT_0472 + conj
      (IT_0170)*IT_0475 + IT_0170*IT_0478 + IT_0135*IT_0481 + IT_0154*IT_0484 +
       IT_0256*IT_0487 + IT_0266*IT_0490)*IT_0504 + 4*IT_0510*(IT_0438*IT_0439 +
       IT_0452*(IT_0451 + 0.25*IT_0511)) + 6*IT_0068*(IT_0439*IT_0452 + IT_0456*
      (conj(IT_0437) + 0.166666666666667*IT_0512)) + (IT_0350*IT_0437 + IT_0347
      *conj(IT_0437) + IT_0325*IT_0452 + IT_0328*IT_0456 + conj(IT_0331)*IT_0467
       + IT_0331*IT_0470 + conj(IT_0303)*IT_0479 + IT_0303*IT_0482)*IT_0513 + 
      (conj(IT_0321)*IT_0322 + conj(IT_0323)*IT_0334 + IT_0323*IT_0336 + conj
      (IT_0331)*IT_0345 + conj(IT_0303)*IT_0354 + IT_0303*IT_0356)*IT_0514 +
       IT_0515*IT_0516 + IT_0518*(conj(IT_0343)*IT_0368 + conj(IT_0342)*IT_0376 
      + IT_0342*IT_0378 + conj(IT_0344)*IT_0455 + IT_0344*IT_0458 + conj(IT_0351
      )*IT_0519) + IT_0523*IT_0527 + IT_0529*IT_0560 + IT_0561*IT_0562 + (-4)
      *IT_0564*(IT_0136*IT_0156 + IT_0153*IT_0158 + -conj(IT_0170)*IT_0268 +
       IT_0265*IT_0269 + IT_0256*(IT_0255 + IT_0267 + (-0.25)*IT_0565));
    return create_ccomplex_return(IT_0566);
}

