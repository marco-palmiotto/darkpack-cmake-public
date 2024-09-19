#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_4_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_4_to_H0_Z(
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
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = s_12*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(m_Z, 4);
    const ccomplex_t IT_0004 = s_12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = s_14*s_24;
    const ccomplex_t IT_0008 = (-2)*IT_0007;
    const ccomplex_t IT_0009 = IT_0001 + IT_0006 + IT_0008;
    const ccomplex_t IT_0010 = cos(alpha);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = sin(alpha);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = sin(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = IT_0016*IT_0021;
    const ccomplex_t IT_0026 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0015 + -IT_0019 + 
      -IT_0024 + IT_0027);
    const ccomplex_t IT_0029 = IT_0012*IT_0020;
    const ccomplex_t IT_0030 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = IT_0011*IT_0021;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0035 = IT_0029*IT_0034;
    const ccomplex_t IT_0036 = IT_0032*IT_0034;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0031 + IT_0033 + 
      -IT_0035 + -IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0028*IT_0038;
    const ccomplex_t IT_0040 = pow(m_H0, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0040 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0039*IT_0042;
    const ccomplex_t IT_0044 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0045 = IT_0013*IT_0044;
    const ccomplex_t IT_0046 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0047 = IT_0025*IT_0046;
    const ccomplex_t IT_0048 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0049 = IT_0013*IT_0048;
    const ccomplex_t IT_0050 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0051 = IT_0025*IT_0050;
    const ccomplex_t IT_0052 = IT_0045 + IT_0047 + IT_0049 + IT_0051;
    const ccomplex_t IT_0053 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0054 = IT_0017*IT_0053;
    const ccomplex_t IT_0055 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0056 = IT_0022*IT_0055;
    const ccomplex_t IT_0057 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0058 = IT_0017*IT_0057;
    const ccomplex_t IT_0059 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0060 = IT_0022*IT_0059;
    const ccomplex_t IT_0061 = -IT_0054 + -IT_0056 + -IT_0058 + -IT_0060;
    const ccomplex_t IT_0062 = IT_0052 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0066 = IT_0029*IT_0065;
    const ccomplex_t IT_0067 = IT_0032*IT_0065;
    const ccomplex_t IT_0068 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0069 = IT_0029*IT_0068;
    const ccomplex_t IT_0070 = IT_0032*IT_0068;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + IT_0067 + 
      -IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0064*IT_0072;
    const ccomplex_t IT_0074 = pow(m_N_4, 2);
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0040 + IT_0074 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0073*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0071;
    const ccomplex_t IT_0079 = IT_0064*IT_0078;
    const ccomplex_t IT_0080 = IT_0076*IT_0079;
    const ccomplex_t IT_0081 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0082 = IT_0013*IT_0081;
    const ccomplex_t IT_0083 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0084 = IT_0025*IT_0083;
    const ccomplex_t IT_0085 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0086 = IT_0013*IT_0085;
    const ccomplex_t IT_0087 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0088 = IT_0025*IT_0087;
    const ccomplex_t IT_0089 = IT_0082 + IT_0084 + IT_0086 + IT_0088;
    const ccomplex_t IT_0090 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0091 = IT_0017*IT_0090;
    const ccomplex_t IT_0092 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0093 = IT_0022*IT_0092;
    const ccomplex_t IT_0094 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0095 = IT_0017*IT_0094;
    const ccomplex_t IT_0096 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0097 = IT_0022*IT_0096;
    const ccomplex_t IT_0098 = -IT_0091 + -IT_0093 + -IT_0095 + -IT_0097;
    const ccomplex_t IT_0099 = IT_0089 + IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0103 = IT_0029*IT_0102;
    const ccomplex_t IT_0104 = IT_0032*IT_0102;
    const ccomplex_t IT_0105 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0106 = IT_0029*IT_0105;
    const ccomplex_t IT_0107 = IT_0032*IT_0105;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0103 + IT_0104 + 
      -IT_0106 + -IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = IT_0101*IT_0109;
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0040 + IT_0074 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0115 = IT_0013*IT_0114;
    const ccomplex_t IT_0116 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0117 = IT_0025*IT_0116;
    const ccomplex_t IT_0118 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0119 = IT_0013*IT_0118;
    const ccomplex_t IT_0120 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0121 = IT_0025*IT_0120;
    const ccomplex_t IT_0122 = IT_0115 + IT_0117 + IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0124 = IT_0017*IT_0123;
    const ccomplex_t IT_0125 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0126 = IT_0022*IT_0125;
    const ccomplex_t IT_0127 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0128 = IT_0017*IT_0127;
    const ccomplex_t IT_0129 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0130 = IT_0022*IT_0129;
    const ccomplex_t IT_0131 = -IT_0124 + -IT_0126 + -IT_0128 + -IT_0130;
    const ccomplex_t IT_0132 = IT_0122 + IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0136 = IT_0029*IT_0135;
    const ccomplex_t IT_0137 = IT_0032*IT_0135;
    const ccomplex_t IT_0138 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0139 = IT_0029*IT_0138;
    const ccomplex_t IT_0140 = IT_0032*IT_0138;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0136 + IT_0137 + 
      -IT_0139 + -IT_0140);
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = IT_0134*IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_23 + IT_0040 + IT_0074 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = -IT_0043 + (-0.5)*IT_0077 + 0.5*IT_0080 + 
      -IT_0113 + -IT_0146;
    const ccomplex_t IT_0148 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0149 = IT_0013*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0151 = IT_0025*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0153 = IT_0013*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0155 = IT_0025*IT_0154;
    const ccomplex_t IT_0156 = IT_0149 + IT_0151 + IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0158 = IT_0017*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0160 = IT_0022*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0162 = IT_0017*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0164 = IT_0022*IT_0163;
    const ccomplex_t IT_0165 = -IT_0158 + -IT_0160 + -IT_0162 + -IT_0164;
    const ccomplex_t IT_0166 = IT_0156 + IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0170 = IT_0029*IT_0169;
    const ccomplex_t IT_0171 = IT_0032*IT_0169;
    const ccomplex_t IT_0172 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0173 = IT_0029*IT_0172;
    const ccomplex_t IT_0174 = IT_0032*IT_0172;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0170 + IT_0171 + 
      -IT_0173 + -IT_0174);
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = IT_0168*IT_0176;
    const ccomplex_t IT_0178 = IT_0145*IT_0177;
    const ccomplex_t IT_0179 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0180 = IT_0013*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0182 = IT_0017*IT_0181;
    const ccomplex_t IT_0183 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0184 = IT_0022*IT_0183;
    const ccomplex_t IT_0185 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0186 = IT_0025*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0180 + -IT_0182 + 
      -IT_0184 + IT_0186);
    const ccomplex_t IT_0188 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0189 = IT_0029*IT_0188;
    const ccomplex_t IT_0190 = IT_0032*IT_0188;
    const ccomplex_t IT_0191 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0192 = IT_0029*IT_0191;
    const ccomplex_t IT_0193 = IT_0032*IT_0191;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0189 + IT_0190 + 
      -IT_0192 + -IT_0193);
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = IT_0187*IT_0195;
    const ccomplex_t IT_0197 = IT_0042*IT_0196;
    const ccomplex_t IT_0198 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0199 = IT_0013*IT_0198;
    const ccomplex_t IT_0200 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0201 = IT_0025*IT_0200;
    const ccomplex_t IT_0202 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0203 = IT_0013*IT_0202;
    const ccomplex_t IT_0204 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0205 = IT_0025*IT_0204;
    const ccomplex_t IT_0206 = IT_0199 + IT_0201 + IT_0203 + IT_0205;
    const ccomplex_t IT_0207 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0208 = IT_0017*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0210 = IT_0022*IT_0209;
    const ccomplex_t IT_0211 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0212 = IT_0017*IT_0211;
    const ccomplex_t IT_0213 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0214 = IT_0022*IT_0213;
    const ccomplex_t IT_0215 = -IT_0208 + -IT_0210 + -IT_0212 + -IT_0214;
    const ccomplex_t IT_0216 = IT_0206 + IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = 0.5*IT_0217;
    const ccomplex_t IT_0219 = IT_0078*IT_0218;
    const ccomplex_t IT_0220 = IT_0076*IT_0219;
    const ccomplex_t IT_0221 = IT_0072*IT_0218;
    const ccomplex_t IT_0222 = IT_0076*IT_0221;
    const ccomplex_t IT_0223 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0224 = IT_0013*IT_0223;
    const ccomplex_t IT_0225 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0226 = IT_0025*IT_0225;
    const ccomplex_t IT_0227 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0228 = IT_0013*IT_0227;
    const ccomplex_t IT_0229 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0230 = IT_0025*IT_0229;
    const ccomplex_t IT_0231 = IT_0224 + IT_0226 + IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0233 = IT_0017*IT_0232;
    const ccomplex_t IT_0234 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0235 = IT_0022*IT_0234;
    const ccomplex_t IT_0236 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0237 = IT_0017*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0239 = IT_0022*IT_0238;
    const ccomplex_t IT_0240 = -IT_0233 + -IT_0235 + -IT_0237 + -IT_0239;
    const ccomplex_t IT_0241 = IT_0231 + IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = 0.5*IT_0242;
    const ccomplex_t IT_0244 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0245 = IT_0029*IT_0244;
    const ccomplex_t IT_0246 = IT_0032*IT_0244;
    const ccomplex_t IT_0247 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0248 = IT_0029*IT_0247;
    const ccomplex_t IT_0249 = IT_0032*IT_0247;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0245 + IT_0246 + 
      -IT_0248 + -IT_0249);
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = IT_0243*IT_0251;
    const ccomplex_t IT_0253 = IT_0112*IT_0252;
    const ccomplex_t IT_0254 = -IT_0178 + -IT_0197 + (-0.5)*IT_0220 + 0.5
      *IT_0222 + -IT_0253;
    const ccomplex_t IT_0255 = pow(m_N_1, 2);
    const ccomplex_t IT_0256 = cpow((-2)*s_13 + IT_0040 + IT_0255 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0259 = IT_0029*IT_0258;
    const ccomplex_t IT_0260 = IT_0032*IT_0258;
    const ccomplex_t IT_0261 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0262 = IT_0029*IT_0261;
    const ccomplex_t IT_0263 = IT_0032*IT_0261;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0259 + IT_0260 + 
      -IT_0262 + -IT_0263);
    const ccomplex_t IT_0265 = (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = IT_0064*IT_0265;
    const ccomplex_t IT_0267 = IT_0257*IT_0266;
    const ccomplex_t IT_0268 = (-0.5)*IT_0267;
    const ccomplex_t IT_0269 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0270 = IT_0013*IT_0269;
    const ccomplex_t IT_0271 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0272 = IT_0025*IT_0271;
    const ccomplex_t IT_0273 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0274 = IT_0013*IT_0273;
    const ccomplex_t IT_0275 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0276 = IT_0025*IT_0275;
    const ccomplex_t IT_0277 = IT_0270 + IT_0272 + IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0279 = IT_0017*IT_0278;
    const ccomplex_t IT_0280 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0281 = IT_0022*IT_0280;
    const ccomplex_t IT_0282 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0283 = IT_0017*IT_0282;
    const ccomplex_t IT_0284 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0285 = IT_0022*IT_0284;
    const ccomplex_t IT_0286 = -IT_0279 + -IT_0281 + -IT_0283 + -IT_0285;
    const ccomplex_t IT_0287 = IT_0277 + IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*IT_0287;
    const ccomplex_t IT_0289 = 0.5*IT_0288;
    const ccomplex_t IT_0290 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0291 = IT_0029*IT_0290;
    const ccomplex_t IT_0292 = IT_0032*IT_0290;
    const ccomplex_t IT_0293 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0294 = IT_0029*IT_0293;
    const ccomplex_t IT_0295 = IT_0032*IT_0293;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*(IT_0291 + IT_0292 + 
      -IT_0294 + -IT_0295);
    const ccomplex_t IT_0297 = (-0.5)*IT_0296;
    const ccomplex_t IT_0298 = IT_0289*IT_0297;
    const ccomplex_t IT_0299 = cpow((-2)*s_13 + IT_0040 + IT_0255 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0298*IT_0300;
    const ccomplex_t IT_0302 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0303 = IT_0013*IT_0302;
    const ccomplex_t IT_0304 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0305 = IT_0025*IT_0304;
    const ccomplex_t IT_0306 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0307 = IT_0013*IT_0306;
    const ccomplex_t IT_0308 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0309 = IT_0025*IT_0308;
    const ccomplex_t IT_0310 = IT_0303 + IT_0305 + IT_0307 + IT_0309;
    const ccomplex_t IT_0311 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0312 = IT_0017*IT_0311;
    const ccomplex_t IT_0313 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0314 = IT_0022*IT_0313;
    const ccomplex_t IT_0315 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0316 = IT_0017*IT_0315;
    const ccomplex_t IT_0317 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0318 = IT_0022*IT_0317;
    const ccomplex_t IT_0319 = -IT_0312 + -IT_0314 + -IT_0316 + -IT_0318;
    const ccomplex_t IT_0320 = IT_0310 + IT_0319;
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*IT_0320;
    const ccomplex_t IT_0322 = 0.5*IT_0321;
    const ccomplex_t IT_0323 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0324 = IT_0029*IT_0323;
    const ccomplex_t IT_0325 = IT_0032*IT_0323;
    const ccomplex_t IT_0326 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0327 = IT_0029*IT_0326;
    const ccomplex_t IT_0328 = IT_0032*IT_0326;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*(IT_0324 + IT_0325 + 
      -IT_0327 + -IT_0328);
    const ccomplex_t IT_0330 = (-0.5)*IT_0329;
    const ccomplex_t IT_0331 = IT_0322*IT_0330;
    const ccomplex_t IT_0332 = cpow((-2)*s_13 + IT_0040 + IT_0255 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*IT_0332;
    const ccomplex_t IT_0334 = IT_0331*IT_0333;
    const ccomplex_t IT_0335 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0336 = IT_0013*IT_0335;
    const ccomplex_t IT_0337 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0338 = IT_0017*IT_0337;
    const ccomplex_t IT_0339 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0340 = IT_0022*IT_0339;
    const ccomplex_t IT_0341 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0342 = IT_0025*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*(IT_0336 + -IT_0338 + 
      -IT_0340 + IT_0342);
    const ccomplex_t IT_0344 = IT_0195*IT_0343;
    const ccomplex_t IT_0345 = cpow((-2)*s_13 + IT_0040 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*IT_0345;
    const ccomplex_t IT_0347 = IT_0344*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0264;
    const ccomplex_t IT_0349 = IT_0064*IT_0348;
    const ccomplex_t IT_0350 = IT_0257*IT_0349;
    const ccomplex_t IT_0351 = -IT_0301 + -IT_0334 + -IT_0347 + (-0.5)*IT_0350;
    const ccomplex_t IT_0352 = -IT_0351;
    const ccomplex_t IT_0353 = IT_0268 + IT_0352;
    const ccomplex_t IT_0354 = -conj(IT_0351);
    const ccomplex_t IT_0355 = conj(IT_0268) + IT_0354;
    const ccomplex_t IT_0356 = IT_0218*IT_0348;
    const ccomplex_t IT_0357 = IT_0257*IT_0356;
    const ccomplex_t IT_0358 = (-0.5)*IT_0357;
    const ccomplex_t IT_0359 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0360 = IT_0013*IT_0359;
    const ccomplex_t IT_0361 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0362 = IT_0025*IT_0361;
    const ccomplex_t IT_0363 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0364 = IT_0013*IT_0363;
    const ccomplex_t IT_0365 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0366 = IT_0025*IT_0365;
    const ccomplex_t IT_0367 = IT_0360 + IT_0362 + IT_0364 + IT_0366;
    const ccomplex_t IT_0368 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0369 = IT_0017*IT_0368;
    const ccomplex_t IT_0370 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0371 = IT_0022*IT_0370;
    const ccomplex_t IT_0372 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0373 = IT_0017*IT_0372;
    const ccomplex_t IT_0374 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0375 = IT_0022*IT_0374;
    const ccomplex_t IT_0376 = -IT_0369 + -IT_0371 + -IT_0373 + -IT_0375;
    const ccomplex_t IT_0377 = IT_0367 + IT_0376;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*IT_0377;
    const ccomplex_t IT_0379 = 0.5*IT_0378;
    const ccomplex_t IT_0380 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0381 = IT_0029*IT_0380;
    const ccomplex_t IT_0382 = IT_0032*IT_0380;
    const ccomplex_t IT_0383 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0384 = IT_0029*IT_0383;
    const ccomplex_t IT_0385 = IT_0032*IT_0383;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*(IT_0381 + IT_0382 + 
      -IT_0384 + -IT_0385);
    const ccomplex_t IT_0387 = (-0.5)*IT_0386;
    const ccomplex_t IT_0388 = IT_0379*IT_0387;
    const ccomplex_t IT_0389 = IT_0300*IT_0388;
    const ccomplex_t IT_0390 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0391 = IT_0013*IT_0390;
    const ccomplex_t IT_0392 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0393 = IT_0025*IT_0392;
    const ccomplex_t IT_0394 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0395 = IT_0013*IT_0394;
    const ccomplex_t IT_0396 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0397 = IT_0025*IT_0396;
    const ccomplex_t IT_0398 = IT_0391 + IT_0393 + IT_0395 + IT_0397;
    const ccomplex_t IT_0399 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0400 = IT_0017*IT_0399;
    const ccomplex_t IT_0401 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0402 = IT_0022*IT_0401;
    const ccomplex_t IT_0403 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0404 = IT_0017*IT_0403;
    const ccomplex_t IT_0405 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0406 = IT_0022*IT_0405;
    const ccomplex_t IT_0407 = -IT_0400 + -IT_0402 + -IT_0404 + -IT_0406;
    const ccomplex_t IT_0408 = IT_0398 + IT_0407;
    const ccomplex_t IT_0409 = (0 + _Complex_I*1)*IT_0408;
    const ccomplex_t IT_0410 = 0.5*IT_0409;
    const ccomplex_t IT_0411 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0412 = IT_0029*IT_0411;
    const ccomplex_t IT_0413 = IT_0032*IT_0411;
    const ccomplex_t IT_0414 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0415 = IT_0029*IT_0414;
    const ccomplex_t IT_0416 = IT_0032*IT_0414;
    const ccomplex_t IT_0417 = (0 + _Complex_I*1)*(IT_0412 + IT_0413 + 
      -IT_0415 + -IT_0416);
    const ccomplex_t IT_0418 = (-0.5)*IT_0417;
    const ccomplex_t IT_0419 = IT_0410*IT_0418;
    const ccomplex_t IT_0420 = IT_0333*IT_0419;
    const ccomplex_t IT_0421 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0422 = IT_0013*IT_0421;
    const ccomplex_t IT_0423 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0424 = IT_0017*IT_0423;
    const ccomplex_t IT_0425 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0426 = IT_0022*IT_0425;
    const ccomplex_t IT_0427 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0428 = IT_0025*IT_0427;
    const ccomplex_t IT_0429 = (0 + _Complex_I*1)*(IT_0422 + -IT_0424 + 
      -IT_0426 + IT_0428);
    const ccomplex_t IT_0430 = IT_0038*IT_0429;
    const ccomplex_t IT_0431 = IT_0346*IT_0430;
    const ccomplex_t IT_0432 = IT_0218*IT_0265;
    const ccomplex_t IT_0433 = IT_0257*IT_0432;
    const ccomplex_t IT_0434 = -IT_0389 + -IT_0420 + -IT_0431 + (-0.5)*IT_0433;
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = IT_0358 + IT_0435;
    const ccomplex_t IT_0437 = -conj(IT_0434);
    const ccomplex_t IT_0438 = conj(IT_0358) + IT_0437;
    const ccomplex_t IT_0439 = m_N_1*m_N_4;
    const ccomplex_t IT_0440 = pow(s_24, 2);
    const ccomplex_t IT_0441 = IT_0002*IT_0440;
    const ccomplex_t IT_0442 = -IT_0074;
    const ccomplex_t IT_0443 = IT_0441 + IT_0442;
    const ccomplex_t IT_0444 = IT_0439*IT_0443;
    const ccomplex_t IT_0445 = 2*IT_0389 + 2*IT_0420 + 2*IT_0431 + IT_0433;
    const ccomplex_t IT_0446 = cos(beta);
    const ccomplex_t IT_0447 = IT_0012*IT_0446;
    const ccomplex_t IT_0448 = IT_0044*IT_0447;
    const ccomplex_t IT_0449 = sin(beta);
    const ccomplex_t IT_0450 = IT_0012*IT_0449;
    const ccomplex_t IT_0451 = IT_0053*IT_0450;
    const ccomplex_t IT_0452 = IT_0021*IT_0446;
    const ccomplex_t IT_0453 = IT_0055*IT_0452;
    const ccomplex_t IT_0454 = IT_0021*IT_0449;
    const ccomplex_t IT_0455 = IT_0046*IT_0454;
    const ccomplex_t IT_0456 = IT_0048*IT_0447;
    const ccomplex_t IT_0457 = IT_0057*IT_0450;
    const ccomplex_t IT_0458 = IT_0059*IT_0452;
    const ccomplex_t IT_0459 = IT_0050*IT_0454;
    const ccomplex_t IT_0460 = (-0.5)*IT_0448 + (-0.5)*IT_0451 + 0.5*IT_0453 +
       0.5*IT_0455 + (-0.5)*IT_0456 + (-0.5)*IT_0457 + 0.5*IT_0458 + 0.5*IT_0459;
    const ccomplex_t IT_0461 = e_em*IT_0012*IT_0021*(IT_0010*IT_0446 + IT_0016
      *IT_0449);
    const ccomplex_t IT_0462 = -IT_0461;
    const ccomplex_t IT_0463 = IT_0460*IT_0462;
    const ccomplex_t IT_0464 = cpow((-2)*s_12 + IT_0000 + -IT_0074 + -IT_0255 
      + -reg_prop, -1);
    const ccomplex_t IT_0465 = (0 + _Complex_I*1)*IT_0464;
    const ccomplex_t IT_0466 = IT_0463*IT_0465;
    const ccomplex_t IT_0467 = IT_0053*IT_0447;
    const ccomplex_t IT_0468 = IT_0044*IT_0450;
    const ccomplex_t IT_0469 = IT_0046*IT_0452;
    const ccomplex_t IT_0470 = IT_0055*IT_0454;
    const ccomplex_t IT_0471 = IT_0057*IT_0447;
    const ccomplex_t IT_0472 = IT_0048*IT_0450;
    const ccomplex_t IT_0473 = IT_0050*IT_0452;
    const ccomplex_t IT_0474 = IT_0059*IT_0454;
    const ccomplex_t IT_0475 = (-0.5)*IT_0467 + 0.5*IT_0468 + 0.5*IT_0469 + (
      -0.5)*IT_0470 + (-0.5)*IT_0471 + 0.5*IT_0472 + 0.5*IT_0473 + (-0.5)*IT_0474;
    const ccomplex_t IT_0476 = -e_em*IT_0012*IT_0021*(IT_0016*IT_0446 + 
      -IT_0010*IT_0449);
    const ccomplex_t IT_0477 = IT_0475*IT_0476;
    const ccomplex_t IT_0478 = cpow((-2)*s_12 + -IT_0074 + -IT_0255 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0479 = (0 + _Complex_I*1)*IT_0478;
    const ccomplex_t IT_0480 = IT_0477*IT_0479;
    const ccomplex_t IT_0481 = IT_0267 + IT_0466 + IT_0480;
    const ccomplex_t IT_0482 = 2*IT_0301 + 2*IT_0334 + 2*IT_0347 + IT_0350;
    const ccomplex_t IT_0483 = (-2)*IT_0482;
    const ccomplex_t IT_0484 = 0.5*IT_0483;
    const ccomplex_t IT_0485 = IT_0481 + IT_0484;
    const ccomplex_t IT_0486 = (-2)*conj(IT_0482);
    const ccomplex_t IT_0487 = 0.5*IT_0486;
    const ccomplex_t IT_0488 = conj(IT_0481) + IT_0487;
    const ccomplex_t IT_0489 = IT_0198*IT_0447;
    const ccomplex_t IT_0490 = IT_0207*IT_0450;
    const ccomplex_t IT_0491 = IT_0209*IT_0452;
    const ccomplex_t IT_0492 = IT_0200*IT_0454;
    const ccomplex_t IT_0493 = IT_0202*IT_0447;
    const ccomplex_t IT_0494 = IT_0211*IT_0450;
    const ccomplex_t IT_0495 = IT_0213*IT_0452;
    const ccomplex_t IT_0496 = IT_0204*IT_0454;
    const ccomplex_t IT_0497 = 0.5*IT_0489 + 0.5*IT_0490 + (-0.5)*IT_0491 + (
      -0.5)*IT_0492 + 0.5*IT_0493 + 0.5*IT_0494 + (-0.5)*IT_0495 + (-0.5)*IT_0496;
    const ccomplex_t IT_0498 = IT_0462*IT_0497;
    const ccomplex_t IT_0499 = IT_0465*IT_0498;
    const ccomplex_t IT_0500 = IT_0207*IT_0447;
    const ccomplex_t IT_0501 = IT_0198*IT_0450;
    const ccomplex_t IT_0502 = IT_0200*IT_0452;
    const ccomplex_t IT_0503 = IT_0209*IT_0454;
    const ccomplex_t IT_0504 = IT_0211*IT_0447;
    const ccomplex_t IT_0505 = IT_0202*IT_0450;
    const ccomplex_t IT_0506 = IT_0204*IT_0452;
    const ccomplex_t IT_0507 = IT_0213*IT_0454;
    const ccomplex_t IT_0508 = 0.5*IT_0500 + (-0.5)*IT_0501 + (-0.5)*IT_0502 +
       0.5*IT_0503 + 0.5*IT_0504 + (-0.5)*IT_0505 + (-0.5)*IT_0506 + 0.5*IT_0507;
    const ccomplex_t IT_0509 = IT_0476*IT_0508;
    const ccomplex_t IT_0510 = IT_0479*IT_0509;
    const ccomplex_t IT_0511 = IT_0357 + IT_0499 + IT_0510;
    const ccomplex_t IT_0512 = (-2)*IT_0481;
    const ccomplex_t IT_0513 = 0.5*IT_0512;
    const ccomplex_t IT_0514 = IT_0482 + IT_0513;
    const ccomplex_t IT_0515 = (-2)*conj(IT_0481);
    const ccomplex_t IT_0516 = 0.5*IT_0515;
    const ccomplex_t IT_0517 = conj(IT_0482) + IT_0516;
    const ccomplex_t IT_0518 = 2*IT_0043 + IT_0077 + -IT_0080 + 2*IT_0113 + 2
      *IT_0146 + IT_0466 + IT_0480;
    const ccomplex_t IT_0519 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0520 = -IT_0519;
    const ccomplex_t IT_0521 = s_12 + IT_0520;
    const ccomplex_t IT_0522 = s_12*IT_0521;
    const ccomplex_t IT_0523 = IT_0439*IT_0521;
    const ccomplex_t IT_0524 = m_N_1*IT_0521;
    const ccomplex_t IT_0525 = m_N_1*IT_0039;
    const ccomplex_t IT_0526 = IT_0042*IT_0525;
    const ccomplex_t IT_0527 = m_N_1*IT_0073;
    const ccomplex_t IT_0528 = IT_0076*IT_0527;
    const ccomplex_t IT_0529 = m_N_1*IT_0079;
    const ccomplex_t IT_0530 = -IT_0529;
    const ccomplex_t IT_0531 = IT_0527 + IT_0530;
    const ccomplex_t IT_0532 = IT_0076*IT_0531;
    const ccomplex_t IT_0533 = m_N_1*IT_0110;
    const ccomplex_t IT_0534 = IT_0112*IT_0533;
    const ccomplex_t IT_0535 = m_N_1*IT_0143;
    const ccomplex_t IT_0536 = IT_0145*IT_0535;
    const ccomplex_t IT_0537 = m_N_4*IT_0356;
    const ccomplex_t IT_0538 = m_N_4*IT_0432;
    const ccomplex_t IT_0539 = -IT_0538;
    const ccomplex_t IT_0540 = IT_0537 + IT_0539;
    const ccomplex_t IT_0541 = IT_0257*IT_0540;
    const ccomplex_t IT_0542 = IT_0526 + 0.5*IT_0528 + 0.5*IT_0532 + IT_0534 +
       IT_0536 + (-0.5)*IT_0541;
    const ccomplex_t IT_0543 = cpow(IT_0011, -4);
    const ccomplex_t IT_0544 = IT_0010*IT_0020*IT_0446*IT_0543;
    const ccomplex_t IT_0545 = IT_0016*IT_0020*IT_0449*IT_0543;
    const ccomplex_t IT_0546 = cpow(IT_0020, -3);
    const ccomplex_t IT_0547 = IT_0010*IT_0446*IT_0546;
    const ccomplex_t IT_0548 = IT_0016*IT_0449*IT_0546;
    const ccomplex_t IT_0549 = cpow(IT_0011, -2);
    const ccomplex_t IT_0550 = IT_0010*IT_0021*IT_0446*IT_0549;
    const ccomplex_t IT_0551 = IT_0016*IT_0021*IT_0449*IT_0549;
    const ccomplex_t IT_0552 = cpow(IT_0020, -2);
    const ccomplex_t IT_0553 = cpow(IT_0549 + IT_0552, -1);
    const ccomplex_t IT_0554 = (IT_0544 + IT_0545 + IT_0547 + IT_0548 + 2
      *IT_0550 + 2*IT_0551)*IT_0553;
    const ccomplex_t IT_0555 = m_W*e_em;
    const ccomplex_t IT_0556 = IT_0554*IT_0555;
    const ccomplex_t IT_0557 = (0 + _Complex_I*1)*IT_0556;
    const ccomplex_t IT_0558 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0074 + -IT_0255 + -reg_prop, -1);
    const ccomplex_t IT_0559 = IT_0557*IT_0558;
    const ccomplex_t IT_0560 = (0 + _Complex_I*1)*IT_0559;
    const ccomplex_t IT_0561 = IT_0195*IT_0560;
    const ccomplex_t IT_0562 = IT_0076*IT_0529;
    const ccomplex_t IT_0563 = IT_0101*IT_0251;
    const ccomplex_t IT_0564 = m_N_2*IT_0563;
    const ccomplex_t IT_0565 = IT_0112*IT_0564;
    const ccomplex_t IT_0566 = IT_0134*IT_0176;
    const ccomplex_t IT_0567 = m_N_3*IT_0566;
    const ccomplex_t IT_0568 = IT_0145*IT_0567;
    const ccomplex_t IT_0569 = IT_0195*IT_0429;
    const ccomplex_t IT_0570 = m_N_1*IT_0569;
    const ccomplex_t IT_0571 = m_N_4*IT_0430;
    const ccomplex_t IT_0572 = -IT_0571;
    const ccomplex_t IT_0573 = IT_0570 + IT_0572;
    const ccomplex_t IT_0574 = IT_0346*IT_0573;
    const ccomplex_t IT_0575 = IT_0028*IT_0195;
    const ccomplex_t IT_0576 = m_N_4*IT_0575;
    const ccomplex_t IT_0577 = IT_0042*IT_0576;
    const ccomplex_t IT_0578 = IT_0297*IT_0379;
    const ccomplex_t IT_0579 = m_N_2*IT_0578;
    const ccomplex_t IT_0580 = m_N_4*IT_0388;
    const ccomplex_t IT_0581 = -IT_0580;
    const ccomplex_t IT_0582 = IT_0579 + IT_0581;
    const ccomplex_t IT_0583 = IT_0300*IT_0582;
    const ccomplex_t IT_0584 = IT_0330*IT_0410;
    const ccomplex_t IT_0585 = m_N_3*IT_0584;
    const ccomplex_t IT_0586 = m_N_4*IT_0419;
    const ccomplex_t IT_0587 = -IT_0586;
    const ccomplex_t IT_0588 = IT_0585 + IT_0587;
    const ccomplex_t IT_0589 = IT_0333*IT_0588;
    const ccomplex_t IT_0590 = (-0.5)*IT_0541 + -IT_0561 + (-0.5)*IT_0562 + 
      -IT_0565 + -IT_0568 + -IT_0574 + -IT_0577 + -IT_0583 + -IT_0589;
    const ccomplex_t IT_0591 = IT_0542 + IT_0590;
    const ccomplex_t IT_0592 = pow(s_14, 2);
    const ccomplex_t IT_0593 = IT_0002*IT_0592;
    const ccomplex_t IT_0594 = -IT_0255;
    const ccomplex_t IT_0595 = IT_0593 + IT_0594;
    const ccomplex_t IT_0596 = m_N_4*IT_0595;
    const ccomplex_t IT_0597 = m_N_1*IT_0177;
    const ccomplex_t IT_0598 = IT_0145*IT_0597;
    const ccomplex_t IT_0599 = m_N_1*IT_0219;
    const ccomplex_t IT_0600 = IT_0076*IT_0599;
    const ccomplex_t IT_0601 = m_N_1*IT_0221;
    const ccomplex_t IT_0602 = -IT_0601;
    const ccomplex_t IT_0603 = IT_0599 + IT_0602;
    const ccomplex_t IT_0604 = IT_0076*IT_0603;
    const ccomplex_t IT_0605 = m_N_1*IT_0252;
    const ccomplex_t IT_0606 = IT_0112*IT_0605;
    const ccomplex_t IT_0607 = m_N_4*IT_0266;
    const ccomplex_t IT_0608 = m_N_4*IT_0349;
    const ccomplex_t IT_0609 = -IT_0608;
    const ccomplex_t IT_0610 = IT_0607 + IT_0609;
    const ccomplex_t IT_0611 = IT_0257*IT_0610;
    const ccomplex_t IT_0612 = m_N_1*IT_0196;
    const ccomplex_t IT_0613 = IT_0042*IT_0612;
    const ccomplex_t IT_0614 = IT_0598 + 0.5*IT_0600 + 0.5*IT_0604 + IT_0606 +
       (-0.5)*IT_0611 + IT_0613;
    const ccomplex_t IT_0615 = IT_0038*IT_0560;
    const ccomplex_t IT_0616 = IT_0076*IT_0601;
    const ccomplex_t IT_0617 = IT_0109*IT_0243;
    const ccomplex_t IT_0618 = m_N_2*IT_0617;
    const ccomplex_t IT_0619 = IT_0112*IT_0618;
    const ccomplex_t IT_0620 = IT_0142*IT_0168;
    const ccomplex_t IT_0621 = m_N_3*IT_0620;
    const ccomplex_t IT_0622 = IT_0145*IT_0621;
    const ccomplex_t IT_0623 = IT_0038*IT_0343;
    const ccomplex_t IT_0624 = m_N_1*IT_0623;
    const ccomplex_t IT_0625 = m_N_4*IT_0344;
    const ccomplex_t IT_0626 = -IT_0625;
    const ccomplex_t IT_0627 = IT_0624 + IT_0626;
    const ccomplex_t IT_0628 = IT_0346*IT_0627;
    const ccomplex_t IT_0629 = IT_0038*IT_0187;
    const ccomplex_t IT_0630 = m_N_4*IT_0629;
    const ccomplex_t IT_0631 = IT_0042*IT_0630;
    const ccomplex_t IT_0632 = IT_0289*IT_0387;
    const ccomplex_t IT_0633 = m_N_2*IT_0632;
    const ccomplex_t IT_0634 = m_N_4*IT_0298;
    const ccomplex_t IT_0635 = -IT_0634;
    const ccomplex_t IT_0636 = IT_0633 + IT_0635;
    const ccomplex_t IT_0637 = IT_0300*IT_0636;
    const ccomplex_t IT_0638 = IT_0322*IT_0418;
    const ccomplex_t IT_0639 = m_N_3*IT_0638;
    const ccomplex_t IT_0640 = m_N_4*IT_0331;
    const ccomplex_t IT_0641 = -IT_0640;
    const ccomplex_t IT_0642 = IT_0639 + IT_0641;
    const ccomplex_t IT_0643 = IT_0333*IT_0642;
    const ccomplex_t IT_0644 = (-0.5)*IT_0611 + -IT_0615 + (-0.5)*IT_0616 + 
      -IT_0619 + -IT_0622 + -IT_0628 + -IT_0631 + -IT_0637 + -IT_0643;
    const ccomplex_t IT_0645 = IT_0614 + IT_0644;
    const ccomplex_t IT_0646 = s_12*IT_0595;
    const ccomplex_t IT_0647 = conj(IT_0542) + conj(IT_0590);
    const ccomplex_t IT_0648 = conj(IT_0614) + conj(IT_0644);
    const ccomplex_t IT_0649 = IT_0000*IT_0439;
    const ccomplex_t IT_0650 = 3*IT_0649;
    const ccomplex_t IT_0651 = (-0.25)*IT_0005;
    const ccomplex_t IT_0652 = IT_0007 + IT_0651;
    const ccomplex_t IT_0653 = s_14*IT_0000;
    const ccomplex_t IT_0654 = IT_0002*IT_0653;
    const ccomplex_t IT_0655 = 0.5*IT_0654;
    const ccomplex_t IT_0656 = s_14 + IT_0655;
    const ccomplex_t IT_0657 = m_N_4*IT_0656;
    const ccomplex_t IT_0658 = 2*IT_0657;
    const ccomplex_t IT_0659 = s_24*IT_0000;
    const ccomplex_t IT_0660 = IT_0002*IT_0659;
    const ccomplex_t IT_0661 = (-0.25)*IT_0660;
    const ccomplex_t IT_0662 = s_24 + IT_0661;
    const ccomplex_t IT_0663 = m_N_1*IT_0662;
    const ccomplex_t IT_0664 = (-4)*IT_0663;
    const ccomplex_t IT_0665 = (-2)*IT_0614;
    const ccomplex_t IT_0666 = (-2)*IT_0644;
    const ccomplex_t IT_0667 = IT_0665 + IT_0666;
    const ccomplex_t IT_0668 = (-0.25)*IT_0654;
    const ccomplex_t IT_0669 = s_14 + IT_0668;
    const ccomplex_t IT_0670 = m_N_4*IT_0669;
    const ccomplex_t IT_0671 = (-4)*IT_0670;
    const ccomplex_t IT_0672 = (-2)*IT_0542;
    const ccomplex_t IT_0673 = (-2)*IT_0590;
    const ccomplex_t IT_0674 = IT_0672 + IT_0673;
    const ccomplex_t IT_0675 = (-2)*conj(IT_0542);
    const ccomplex_t IT_0676 = (-2)*conj(IT_0590);
    const ccomplex_t IT_0677 = IT_0675 + IT_0676;
    const ccomplex_t IT_0678 = (-2)*conj(IT_0614);
    const ccomplex_t IT_0679 = (-2)*conj(IT_0644);
    const ccomplex_t IT_0680 = IT_0678 + IT_0679;
    const ccomplex_t IT_0681 = 2*IT_0671*(conj(IT_0358)*IT_0591 + conj(IT_0268
      )*IT_0645 + IT_0358*IT_0647 + IT_0268*IT_0648 + 0.5*conj(IT_0351)*IT_0667 
      + 0.5*conj(IT_0434)*IT_0674 + 0.5*IT_0434*IT_0677 + 0.5*IT_0351*IT_0680);
    const ccomplex_t IT_0682 = IT_0001 + IT_0008;
    const ccomplex_t IT_0683 = IT_0002*IT_0682;
    const ccomplex_t IT_0684 = (-0.5)*IT_0683;
    const ccomplex_t IT_0685 = s_12 + IT_0684;
    const ccomplex_t IT_0686 = 4*IT_0614;
    const ccomplex_t IT_0687 = 6*conj(IT_0542);
    const ccomplex_t IT_0688 = (-2)*IT_0518;
    const ccomplex_t IT_0689 = (-2)*conj(IT_0518);
    const ccomplex_t IT_0690 = 2*IT_0178 + 2*IT_0197 + IT_0220 + -IT_0222 + 2
      *IT_0253 + IT_0499 + IT_0510;
    const ccomplex_t IT_0691 = (-2)*IT_0690;
    const ccomplex_t IT_0692 = (-2)*conj(IT_0690);
    const ccomplex_t IT_0693 = s_12*IT_0443;
    const ccomplex_t IT_0694 = (-2)*conj(IT_0511);
    const ccomplex_t IT_0695 = 0.5*IT_0694;
    const ccomplex_t IT_0696 = conj(IT_0445) + IT_0695;
    const ccomplex_t IT_0697 = (-2)*IT_0511;
    const ccomplex_t IT_0698 = m_N_1*IT_0443;
    const ccomplex_t IT_0699 = m_N_4*IT_0521;
    const ccomplex_t IT_0700 = conj(IT_0445)*IT_0591 + (conj(IT_0482) + 0.5
      *IT_0515)*IT_0645 + IT_0445*IT_0647 + IT_0514*IT_0648 + 0.5*conj(IT_0511)
      *IT_0672 + 0.5*IT_0511*IT_0675 + 0.5*IT_0590*IT_0694 + 0.5*conj(IT_0590)
      *IT_0697;
    const ccomplex_t IT_0701 = (-2)*IT_0147;
    const ccomplex_t IT_0702 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0703 = (-3)*IT_0649;
    const ccomplex_t IT_0704 = (-2)*conj(IT_0268);
    const ccomplex_t IT_0705 = 0.5*IT_0704;
    const ccomplex_t IT_0706 = conj(IT_0351) + IT_0705;
    const ccomplex_t IT_0707 = s_12*s_24;
    const ccomplex_t IT_0708 = s_14*IT_0074;
    const ccomplex_t IT_0709 = -IT_0708;
    const ccomplex_t IT_0710 = IT_0707 + IT_0709;
    const ccomplex_t IT_0711 = IT_0147 + IT_0351;
    const ccomplex_t IT_0712 = conj(IT_0147) + conj(IT_0351);
    const ccomplex_t IT_0713 = conj(IT_0254) + conj(IT_0434);
    const ccomplex_t IT_0714 = IT_0254 + IT_0434;
    const ccomplex_t IT_0715 = (conj(IT_0254) + -conj(IT_0358) + -IT_0437)
      *IT_0445 + (IT_0254 + -IT_0358 + -IT_0435)*conj(IT_0445) + conj(IT_0268)
      *IT_0485 + IT_0268*IT_0488 + conj(IT_0358)*IT_0511 + IT_0358*conj(IT_0511)
       + IT_0517*IT_0711 + IT_0514*IT_0712 + 0.5*IT_0697*IT_0713 + 0.5*IT_0694
      *IT_0714;
    const ccomplex_t IT_0716 = 0.5*IT_0660;
    const ccomplex_t IT_0717 = s_24 + IT_0716;
    const ccomplex_t IT_0718 = m_N_1*IT_0717;
    const ccomplex_t IT_0719 = 2*IT_0718;
    const ccomplex_t IT_0720 = 0.5*conj(IT_0268);
    const ccomplex_t IT_0721 = 0.5*IT_0358;
    const ccomplex_t IT_0722 = 0.5*conj(IT_0358);
    const ccomplex_t IT_0723 = s_12*s_14;
    const ccomplex_t IT_0724 = s_24*IT_0255;
    const ccomplex_t IT_0725 = -IT_0724;
    const ccomplex_t IT_0726 = IT_0723 + IT_0725;
    const ccomplex_t IT_0727 = IT_0439*IT_0595;
    const ccomplex_t IT_0728 = 0.5*IT_0727;
    const ccomplex_t IT_0729 = 0.5*IT_0664;
    const ccomplex_t IT_0730 = (-4)*IT_0009*(IT_0147*conj(IT_0147) + IT_0254
      *conj(IT_0254) + IT_0353*IT_0355 + IT_0436*IT_0438) + 2*IT_0444*(conj
      (IT_0445)*IT_0485 + IT_0445*IT_0488 + conj(IT_0511)*IT_0514 + IT_0511
      *IT_0517) + 2*conj(IT_0518)*(IT_0482*IT_0522 + IT_0511*IT_0523 + IT_0524
      *IT_0591 + IT_0596*IT_0645) + 2*IT_0518*(conj(IT_0482)*IT_0522 + conj
      (IT_0511)*IT_0523 + conj(IT_0518)*IT_0646 + IT_0524*IT_0647 + IT_0596
      *IT_0648) + (-2)*conj(IT_0254)*(IT_0353*IT_0650 + 4*IT_0436*IT_0652 + 
      -IT_0591*IT_0658 + (-0.5)*IT_0664*IT_0667) + IT_0681 + 4*IT_0685*(IT_0591
      *IT_0647 + IT_0648*(IT_0644 + 0.25*IT_0686)) + 6*IT_0439*(IT_0591*IT_0648 
      + IT_0645*(conj(IT_0590) + 0.166666666666667*IT_0687)) + IT_0523*(conj
      (IT_0445)*IT_0688 + IT_0445*IT_0689) + IT_0522*(conj(IT_0481)*IT_0688 +
       IT_0481*IT_0689 + conj(IT_0511)*IT_0691 + IT_0511*IT_0692) + 2*IT_0693*
      (IT_0481*conj(IT_0481) + IT_0511*conj(IT_0511) + 0.5*conj(IT_0482)*IT_0512
       + IT_0482*IT_0517 + IT_0445*IT_0696 + 0.5*conj(IT_0445)*IT_0697) + 2*
      (IT_0517*IT_0591 + IT_0514*IT_0647 + IT_0645*IT_0696 + IT_0648*(IT_0445 +
       0.5*IT_0697))*IT_0698 + 2*IT_0699*IT_0700 + IT_0650*(conj(IT_0358)
      *IT_0701 + IT_0358*IT_0702) + 2*IT_0703*(IT_0355*IT_0434 + -IT_0353
      *IT_0438 + 0.5*conj(IT_0254)*IT_0701 + 0.5*IT_0254*IT_0702 + IT_0358
      *IT_0706) + 2*IT_0710*IT_0715 + 2*IT_0719*(conj(IT_0351)*IT_0591 + conj
      (IT_0434)*IT_0645 + IT_0351*IT_0647 + IT_0434*IT_0648 + 0.5*IT_0268
      *IT_0677 + IT_0674*IT_0720 + IT_0680*IT_0721 + IT_0667*IT_0722) + 2*(0.5
      *IT_0268*IT_0689 + conj(IT_0518)*IT_0711 + IT_0518*IT_0712 + IT_0690
      *IT_0713 + conj(IT_0690)*IT_0714 + IT_0688*IT_0720 + IT_0692*IT_0721 +
       IT_0691*IT_0722)*IT_0726 + 2*conj(IT_0690)*(IT_0445*IT_0522 + IT_0485
      *IT_0523 + IT_0591*IT_0596 + IT_0524*IT_0645 + IT_0688*IT_0728) + 2
      *IT_0690*(conj(IT_0445)*IT_0522 + IT_0488*IT_0523 + IT_0596*IT_0647 +
       IT_0524*IT_0648 + IT_0646*conj(IT_0690) + IT_0689*IT_0728) + 2*conj
      (IT_0147)*(IT_0434*IT_0650 + (-4)*IT_0353*IT_0652 + IT_0645*IT_0658 +
       IT_0674*IT_0729) + 2*IT_0147*(conj(IT_0434)*IT_0650 + (-4)*IT_0355
      *IT_0652 + IT_0648*IT_0658 + IT_0677*IT_0729) + (-8)*IT_0254*(IT_0438
      *IT_0652 + (-0.25)*IT_0647*IT_0658 + (-0.25)*IT_0650*IT_0706 + (-0.25)
      *IT_0680*IT_0729);
    return create_ccomplex_return(IT_0730);
}

