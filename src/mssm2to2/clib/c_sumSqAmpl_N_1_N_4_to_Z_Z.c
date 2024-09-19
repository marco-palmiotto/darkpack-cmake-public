#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_4_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_4_to_Z_Z(
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
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0013 = IT_0003*IT_0012;
    const ccomplex_t IT_0014 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = IT_0005 + IT_0011 + IT_0013 + IT_0015;
    const ccomplex_t IT_0017 = IT_0002*IT_0006;
    const ccomplex_t IT_0018 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0008;
    const ccomplex_t IT_0021 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0016 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = cos(beta);
    const ccomplex_t IT_0032 = cpow(IT_0001, -4);
    const ccomplex_t IT_0033 = IT_0000*IT_0007*IT_0031*IT_0032;
    const ccomplex_t IT_0034 = sin(beta);
    const ccomplex_t IT_0035 = IT_0006*IT_0007*IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0007, -3);
    const ccomplex_t IT_0037 = IT_0000*IT_0031*IT_0036;
    const ccomplex_t IT_0038 = IT_0006*IT_0034*IT_0036;
    const ccomplex_t IT_0039 = cpow(IT_0001, -2);
    const ccomplex_t IT_0040 = IT_0000*IT_0008*IT_0031*IT_0039;
    const ccomplex_t IT_0041 = IT_0006*IT_0008*IT_0034*IT_0039;
    const ccomplex_t IT_0042 = cpow(IT_0007, -2);
    const ccomplex_t IT_0043 = cpow(IT_0039 + IT_0042, -1);
    const ccomplex_t IT_0044 = (IT_0033 + IT_0035 + IT_0037 + IT_0038 + 2
      *IT_0040 + 2*IT_0041)*IT_0043;
    const ccomplex_t IT_0045 = m_W*e_em;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = pow(m_N_1, 2);
    const ccomplex_t IT_0049 = pow(m_N_4, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + -IT_0048 + -IT_0049 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0047*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0030*IT_0052;
    const ccomplex_t IT_0054 = IT_0003*IT_0018;
    const ccomplex_t IT_0055 = IT_0004*IT_0017;
    const ccomplex_t IT_0056 = IT_0003*IT_0023;
    const ccomplex_t IT_0057 = IT_0012*IT_0017;
    const ccomplex_t IT_0058 = IT_0054 + IT_0055 + IT_0056 + IT_0057;
    const ccomplex_t IT_0059 = IT_0010*IT_0020;
    const ccomplex_t IT_0060 = IT_0009*IT_0021;
    const ccomplex_t IT_0061 = IT_0014*IT_0020;
    const ccomplex_t IT_0062 = IT_0009*IT_0025;
    const ccomplex_t IT_0063 = -IT_0059 + -IT_0060 + -IT_0061 + -IT_0062;
    const ccomplex_t IT_0064 = IT_0058 + IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0000*IT_0007*IT_0032*IT_0034;
    const ccomplex_t IT_0068 = IT_0006*IT_0007*IT_0031*IT_0032;
    const ccomplex_t IT_0069 = IT_0000*IT_0034*IT_0036;
    const ccomplex_t IT_0070 = IT_0006*IT_0031*IT_0036;
    const ccomplex_t IT_0071 = IT_0000*IT_0008*IT_0034*IT_0039;
    const ccomplex_t IT_0072 = IT_0006*IT_0008*IT_0031*IT_0039;
    const ccomplex_t IT_0073 = IT_0043*(IT_0067 + -IT_0068 + IT_0069 + 
      -IT_0070 + 2*IT_0071 + (-2)*IT_0072);
    const ccomplex_t IT_0074 = IT_0045*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0048 + -IT_0049 + -reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = IT_0066*IT_0078;
    const ccomplex_t IT_0080 = -IT_0053 + -IT_0079;
    const ccomplex_t IT_0081 = pow(m_Z, -4);
    const ccomplex_t IT_0082 = pow(s_34, 2);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = 1 + IT_0084;
    const ccomplex_t IT_0086 = s_12*IT_0085;
    const ccomplex_t IT_0087 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0088 = IT_0003*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0090 = IT_0009*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0094 = IT_0009*IT_0093;
    const ccomplex_t IT_0095 = IT_0088 + IT_0090 + IT_0092 + IT_0094;
    const ccomplex_t IT_0096 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0097 = IT_0017*IT_0096;
    const ccomplex_t IT_0098 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0099 = IT_0020*IT_0098;
    const ccomplex_t IT_0100 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0101 = IT_0017*IT_0100;
    const ccomplex_t IT_0102 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0103 = IT_0020*IT_0102;
    const ccomplex_t IT_0104 = -IT_0097 + -IT_0099 + -IT_0101 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0095 + IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0106;
    const ccomplex_t IT_0108 = IT_0052*IT_0107;
    const ccomplex_t IT_0109 = IT_0003*IT_0096;
    const ccomplex_t IT_0110 = IT_0017*IT_0087;
    const ccomplex_t IT_0111 = IT_0003*IT_0100;
    const ccomplex_t IT_0112 = IT_0017*IT_0091;
    const ccomplex_t IT_0113 = IT_0109 + IT_0110 + IT_0111 + IT_0112;
    const ccomplex_t IT_0114 = IT_0020*IT_0089;
    const ccomplex_t IT_0115 = IT_0009*IT_0098;
    const ccomplex_t IT_0116 = IT_0020*IT_0093;
    const ccomplex_t IT_0117 = IT_0009*IT_0102;
    const ccomplex_t IT_0118 = -IT_0114 + -IT_0115 + -IT_0116 + -IT_0117;
    const ccomplex_t IT_0119 = IT_0113 + IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0120;
    const ccomplex_t IT_0122 = IT_0078*IT_0121;
    const ccomplex_t IT_0123 = -IT_0108 + -IT_0122;
    const ccomplex_t IT_0124 = m_N_1*m_N_4;
    const ccomplex_t IT_0125 = IT_0085*IT_0124;
    const ccomplex_t IT_0126 = pow(m_Z, -2);
    const ccomplex_t IT_0127 = s_14*s_24*IT_0126;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = s_12 + IT_0128;
    const ccomplex_t IT_0130 = IT_0124*IT_0129;
    const ccomplex_t IT_0131 = IT_0002*IT_0007;
    const ccomplex_t IT_0132 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = IT_0001*IT_0008;
    const ccomplex_t IT_0135 = IT_0132*IT_0134;
    const ccomplex_t IT_0136 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0137 = IT_0131*IT_0136;
    const ccomplex_t IT_0138 = IT_0134*IT_0136;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0133 + IT_0135 + 
      -IT_0137 + -IT_0138);
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0142 = IT_0131*IT_0141;
    const ccomplex_t IT_0143 = IT_0134*IT_0141;
    const ccomplex_t IT_0144 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0145 = IT_0131*IT_0144;
    const ccomplex_t IT_0146 = IT_0134*IT_0144;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0142 + IT_0143 + 
      -IT_0145 + -IT_0146);
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0140*IT_0148;
    const ccomplex_t IT_0150 = pow(m_Z, 2);
    const ccomplex_t IT_0151 = cpow((-2)*s_13 + IT_0048 + IT_0150 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = IT_0149*IT_0152;
    const ccomplex_t IT_0154 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0155 = IT_0131*IT_0154;
    const ccomplex_t IT_0156 = IT_0134*IT_0154;
    const ccomplex_t IT_0157 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0158 = IT_0131*IT_0157;
    const ccomplex_t IT_0159 = IT_0134*IT_0157;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*(IT_0155 + IT_0156 + 
      -IT_0158 + -IT_0159);
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0163 = IT_0131*IT_0162;
    const ccomplex_t IT_0164 = IT_0134*IT_0162;
    const ccomplex_t IT_0165 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0166 = IT_0131*IT_0165;
    const ccomplex_t IT_0167 = IT_0134*IT_0165;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0163 + IT_0164 + 
      -IT_0166 + -IT_0167);
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0161*IT_0169;
    const ccomplex_t IT_0171 = cpow((-2)*s_13 + IT_0048 + IT_0150 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0175 = IT_0131*IT_0174;
    const ccomplex_t IT_0176 = IT_0134*IT_0174;
    const ccomplex_t IT_0177 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0178 = IT_0131*IT_0177;
    const ccomplex_t IT_0179 = IT_0134*IT_0177;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0175 + IT_0176 + 
      -IT_0178 + -IT_0179);
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0183 = IT_0131*IT_0182;
    const ccomplex_t IT_0184 = IT_0134*IT_0182;
    const ccomplex_t IT_0185 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0186 = IT_0131*IT_0185;
    const ccomplex_t IT_0187 = IT_0134*IT_0185;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0183 + IT_0184 + 
      -IT_0186 + -IT_0187);
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = IT_0181*IT_0189;
    const ccomplex_t IT_0191 = cpow((-2)*s_13 + IT_0048 + IT_0150 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0190*IT_0192;
    const ccomplex_t IT_0194 = 0.5*IT_0180;
    const ccomplex_t IT_0195 = IT_0189*IT_0194;
    const ccomplex_t IT_0196 = IT_0192*IT_0195;
    const ccomplex_t IT_0197 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0198 = IT_0131*IT_0197;
    const ccomplex_t IT_0199 = IT_0134*IT_0197;
    const ccomplex_t IT_0200 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0201 = IT_0131*IT_0200;
    const ccomplex_t IT_0202 = IT_0134*IT_0200;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*(IT_0198 + IT_0199 + 
      -IT_0201 + -IT_0202);
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = IT_0189*IT_0204;
    const ccomplex_t IT_0206 = cpow((-2)*s_13 + IT_0150 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = IT_0205*IT_0207;
    const ccomplex_t IT_0209 = 2*IT_0153 + 2*IT_0173 + IT_0193 + -IT_0196 + 2
      *IT_0208;
    const ccomplex_t IT_0210 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0211 = IT_0131*IT_0210;
    const ccomplex_t IT_0212 = IT_0134*IT_0210;
    const ccomplex_t IT_0213 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0214 = IT_0131*IT_0213;
    const ccomplex_t IT_0215 = IT_0134*IT_0213;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0211 + IT_0212 + 
      -IT_0214 + -IT_0215);
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = IT_0181*IT_0217;
    const ccomplex_t IT_0219 = cpow((-2)*s_23 + IT_0150 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = IT_0218*IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0221;
    const ccomplex_t IT_0223 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0224 = IT_0131*IT_0223;
    const ccomplex_t IT_0225 = IT_0134*IT_0223;
    const ccomplex_t IT_0226 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0227 = IT_0131*IT_0226;
    const ccomplex_t IT_0228 = IT_0134*IT_0226;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0224 + IT_0225 + 
      -IT_0227 + -IT_0228);
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0232 = IT_0131*IT_0231;
    const ccomplex_t IT_0233 = IT_0134*IT_0231;
    const ccomplex_t IT_0234 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0235 = IT_0131*IT_0234;
    const ccomplex_t IT_0236 = IT_0134*IT_0234;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0232 + IT_0233 + 
      -IT_0235 + -IT_0236);
    const ccomplex_t IT_0238 = (-0.5)*IT_0237;
    const ccomplex_t IT_0239 = IT_0230*IT_0238;
    const ccomplex_t IT_0240 = cpow((-2)*s_23 + IT_0049 + IT_0150 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = IT_0239*IT_0241;
    const ccomplex_t IT_0243 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0244 = IT_0131*IT_0243;
    const ccomplex_t IT_0245 = IT_0134*IT_0243;
    const ccomplex_t IT_0246 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0247 = IT_0131*IT_0246;
    const ccomplex_t IT_0248 = IT_0134*IT_0246;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*(IT_0244 + IT_0245 + 
      -IT_0247 + -IT_0248);
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0252 = IT_0131*IT_0251;
    const ccomplex_t IT_0253 = IT_0134*IT_0251;
    const ccomplex_t IT_0254 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0255 = IT_0131*IT_0254;
    const ccomplex_t IT_0256 = IT_0134*IT_0254;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*(IT_0252 + IT_0253 + 
      -IT_0255 + -IT_0256);
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0250*IT_0258;
    const ccomplex_t IT_0260 = cpow((-2)*s_23 + IT_0049 + IT_0150 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = IT_0259*IT_0261;
    const ccomplex_t IT_0263 = IT_0204*IT_0217;
    const ccomplex_t IT_0264 = cpow((-2)*s_23 + IT_0049 + IT_0150 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = IT_0263*IT_0265;
    const ccomplex_t IT_0267 = 0.5*IT_0203;
    const ccomplex_t IT_0268 = IT_0217*IT_0267;
    const ccomplex_t IT_0269 = IT_0265*IT_0268;
    const ccomplex_t IT_0270 = (-2)*IT_0242 + (-2)*IT_0262 + IT_0266 + -IT_0269;
    const ccomplex_t IT_0271 = IT_0222 + IT_0270;
    const ccomplex_t IT_0272 = conj(IT_0222) + conj(IT_0270);
    const ccomplex_t IT_0273 = IT_0152*IT_0239;
    const ccomplex_t IT_0274 = IT_0172*IT_0259;
    const ccomplex_t IT_0275 = IT_0194*IT_0217;
    const ccomplex_t IT_0276 = IT_0192*IT_0275;
    const ccomplex_t IT_0277 = IT_0192*IT_0218;
    const ccomplex_t IT_0278 = IT_0207*IT_0268;
    const ccomplex_t IT_0279 = 2*IT_0273 + 2*IT_0274 + IT_0276 + -IT_0277 + 2
      *IT_0278;
    const ccomplex_t IT_0280 = IT_0195*IT_0220;
    const ccomplex_t IT_0281 = 2*IT_0280;
    const ccomplex_t IT_0282 = IT_0149*IT_0241;
    const ccomplex_t IT_0283 = IT_0170*IT_0261;
    const ccomplex_t IT_0284 = IT_0189*IT_0267;
    const ccomplex_t IT_0285 = IT_0265*IT_0284;
    const ccomplex_t IT_0286 = IT_0205*IT_0265;
    const ccomplex_t IT_0287 = (-2)*IT_0282 + (-2)*IT_0283 + IT_0285 + -IT_0286;
    const ccomplex_t IT_0288 = IT_0281 + IT_0287;
    const ccomplex_t IT_0289 = conj(IT_0281) + conj(IT_0287);
    const ccomplex_t IT_0290 = s_12*IT_0150;
    const ccomplex_t IT_0291 = s_13*s_23;
    const ccomplex_t IT_0292 = (-2)*IT_0291;
    const ccomplex_t IT_0293 = IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = IT_0126*IT_0293;
    const ccomplex_t IT_0295 = (s_12 + -IT_0127)*(s_12 + (-0.5)*IT_0294);
    const ccomplex_t IT_0296 = (-2)*IT_0295;
    const ccomplex_t IT_0297 = 2*IT_0222;
    const ccomplex_t IT_0298 = 2*IT_0270;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = 2*conj(IT_0222);
    const ccomplex_t IT_0301 = 2*conj(IT_0270);
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0281;
    const ccomplex_t IT_0304 = 2*IT_0287;
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = 2*conj(IT_0281);
    const ccomplex_t IT_0307 = 2*conj(IT_0287);
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = s_24*IT_0150;
    const ccomplex_t IT_0310 = s_23*s_34;
    const ccomplex_t IT_0311 = (-2)*IT_0310;
    const ccomplex_t IT_0312 = IT_0309 + IT_0311;
    const ccomplex_t IT_0313 = IT_0126*IT_0312;
    const ccomplex_t IT_0314 = pow(m_Z, 4);
    const ccomplex_t IT_0315 = s_24*IT_0314;
    const ccomplex_t IT_0316 = s_23*s_34*IT_0150;
    const ccomplex_t IT_0317 = (-2)*IT_0316;
    const ccomplex_t IT_0318 = IT_0315 + IT_0317;
    const ccomplex_t IT_0319 = IT_0081*IT_0318;
    const ccomplex_t IT_0320 = (-0.666666666666667)*IT_0313 +
       0.166666666666667*IT_0319;
    const ccomplex_t IT_0321 = s_24 + IT_0320;
    const ccomplex_t IT_0322 = m_N_1*IT_0321;
    const ccomplex_t IT_0323 = 6*IT_0322;
    const ccomplex_t IT_0324 = -IT_0221;
    const ccomplex_t IT_0325 = 0.5*IT_0269;
    const ccomplex_t IT_0326 = IT_0242 + (-0.5)*IT_0266;
    const ccomplex_t IT_0327 = conj(IT_0262) + conj(IT_0324) + conj(IT_0325) +
       conj(IT_0326);
    const ccomplex_t IT_0328 = m_N_1*IT_0218;
    const ccomplex_t IT_0329 = m_N_4*IT_0190;
    const ccomplex_t IT_0330 = IT_0220*(IT_0328 + -IT_0329);
    const ccomplex_t IT_0331 = 2*IT_0330;
    const ccomplex_t IT_0332 = m_N_1*IT_0239;
    const ccomplex_t IT_0333 = IT_0140*IT_0238;
    const ccomplex_t IT_0334 = m_N_2*IT_0333;
    const ccomplex_t IT_0335 = IT_0241*(IT_0332 + -IT_0334);
    const ccomplex_t IT_0336 = (-2)*IT_0335;
    const ccomplex_t IT_0337 = m_N_1*IT_0263;
    const ccomplex_t IT_0338 = m_N_1*IT_0268;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = IT_0337 + IT_0339;
    const ccomplex_t IT_0341 = IT_0265*IT_0340;
    const ccomplex_t IT_0342 = IT_0161*IT_0258;
    const ccomplex_t IT_0343 = m_N_3*IT_0342;
    const ccomplex_t IT_0344 = m_N_1*IT_0259;
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = IT_0343 + IT_0345;
    const ccomplex_t IT_0347 = IT_0261*IT_0346;
    const ccomplex_t IT_0348 = -IT_0341 + -IT_0347;
    const ccomplex_t IT_0349 = (-2)*IT_0348;
    const ccomplex_t IT_0350 = IT_0331 + IT_0336 + IT_0349;
    const ccomplex_t IT_0351 = IT_0262 + IT_0324 + IT_0325 + IT_0326;
    const ccomplex_t IT_0352 = 2*conj(IT_0330);
    const ccomplex_t IT_0353 = (-2)*conj(IT_0335);
    const ccomplex_t IT_0354 = (-2)*conj(IT_0348);
    const ccomplex_t IT_0355 = IT_0352 + IT_0353 + IT_0354;
    const ccomplex_t IT_0356 = -IT_0280;
    const ccomplex_t IT_0357 = IT_0282 + (-0.5)*IT_0285 + 0.5*IT_0286;
    const ccomplex_t IT_0358 = conj(IT_0283) + conj(IT_0356) + conj(IT_0357);
    const ccomplex_t IT_0359 = m_N_1*IT_0195;
    const ccomplex_t IT_0360 = m_N_4*IT_0275;
    const ccomplex_t IT_0361 = IT_0220*(IT_0359 + -IT_0360);
    const ccomplex_t IT_0362 = 2*IT_0361;
    const ccomplex_t IT_0363 = m_N_1*IT_0170;
    const ccomplex_t IT_0364 = IT_0169*IT_0250;
    const ccomplex_t IT_0365 = m_N_3*IT_0364;
    const ccomplex_t IT_0366 = IT_0261*(IT_0363 + -IT_0365);
    const ccomplex_t IT_0367 = -IT_0366;
    const ccomplex_t IT_0368 = -IT_0367;
    const ccomplex_t IT_0369 = (-2)*IT_0368;
    const ccomplex_t IT_0370 = m_N_1*IT_0284;
    const ccomplex_t IT_0371 = m_N_1*IT_0205;
    const ccomplex_t IT_0372 = -IT_0371;
    const ccomplex_t IT_0373 = IT_0370 + IT_0372;
    const ccomplex_t IT_0374 = IT_0265*IT_0373;
    const ccomplex_t IT_0375 = IT_0148*IT_0230;
    const ccomplex_t IT_0376 = m_N_2*IT_0375;
    const ccomplex_t IT_0377 = m_N_1*IT_0149;
    const ccomplex_t IT_0378 = -IT_0377;
    const ccomplex_t IT_0379 = IT_0376 + IT_0378;
    const ccomplex_t IT_0380 = IT_0241*IT_0379;
    const ccomplex_t IT_0381 = -IT_0374 + -IT_0380;
    const ccomplex_t IT_0382 = (-2)*IT_0381;
    const ccomplex_t IT_0383 = IT_0362 + IT_0369 + IT_0382;
    const ccomplex_t IT_0384 = IT_0283 + IT_0356 + IT_0357;
    const ccomplex_t IT_0385 = 2*conj(IT_0361);
    const ccomplex_t IT_0386 = (-2)*conj(IT_0368);
    const ccomplex_t IT_0387 = (-2)*conj(IT_0381);
    const ccomplex_t IT_0388 = IT_0385 + IT_0386 + IT_0387;
    const ccomplex_t IT_0389 = m_N_4*IT_0126;
    const ccomplex_t IT_0390 = s_14*IT_0150;
    const ccomplex_t IT_0391 = s_13*s_34;
    const ccomplex_t IT_0392 = s_14*IT_0314;
    const ccomplex_t IT_0393 = IT_0126*IT_0392;
    const ccomplex_t IT_0394 = (-2)*IT_0391 + (-0.5)*IT_0393;
    const ccomplex_t IT_0395 = IT_0390 + IT_0394;
    const ccomplex_t IT_0396 = IT_0389*IT_0395;
    const ccomplex_t IT_0397 = (-2)*IT_0396;
    const ccomplex_t IT_0398 = -IT_0173 + (-0.5)*IT_0193 + 0.5*IT_0196;
    const ccomplex_t IT_0399 = -IT_0153;
    const ccomplex_t IT_0400 = -IT_0208;
    const ccomplex_t IT_0401 = conj(IT_0398) + conj(IT_0399) + conj(IT_0400);
    const ccomplex_t IT_0402 = IT_0398 + IT_0399 + IT_0400;
    const ccomplex_t IT_0403 = -IT_0274;
    const ccomplex_t IT_0404 = -IT_0273 + (-0.5)*IT_0276;
    const ccomplex_t IT_0405 = 0.5*IT_0277;
    const ccomplex_t IT_0406 = -IT_0278;
    const ccomplex_t IT_0407 = conj(IT_0403) + conj(IT_0404) + conj(IT_0405) +
       conj(IT_0406);
    const ccomplex_t IT_0408 = IT_0403 + IT_0404;
    const ccomplex_t IT_0409 = IT_0405 + IT_0406 + IT_0408;
    const ccomplex_t IT_0410 = s_23*s_34*IT_0126;
    const ccomplex_t IT_0411 = (-0.166666666666667)*IT_0319 + (
      -0.666666666666667)*IT_0410;
    const ccomplex_t IT_0412 = s_24 + IT_0411;
    const ccomplex_t IT_0413 = m_N_1*IT_0412;
    const ccomplex_t IT_0414 = (-6)*IT_0413;
    const ccomplex_t IT_0415 = m_N_1*IT_0126;
    const ccomplex_t IT_0416 = IT_0126*IT_0315;
    const ccomplex_t IT_0417 = (-2)*IT_0310 + (-0.5)*IT_0416;
    const ccomplex_t IT_0418 = IT_0309 + IT_0417;
    const ccomplex_t IT_0419 = IT_0415*IT_0418;
    const ccomplex_t IT_0420 = (-2)*IT_0419;
    const ccomplex_t IT_0421 = m_N_4*IT_0149;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0334 + IT_0422;
    const ccomplex_t IT_0424 = IT_0152*IT_0423;
    const ccomplex_t IT_0425 = m_N_4*IT_0170;
    const ccomplex_t IT_0426 = -IT_0425;
    const ccomplex_t IT_0427 = IT_0343 + IT_0426;
    const ccomplex_t IT_0428 = IT_0172*IT_0427;
    const ccomplex_t IT_0429 = m_N_4*IT_0195;
    const ccomplex_t IT_0430 = -IT_0429;
    const ccomplex_t IT_0431 = IT_0329 + IT_0430;
    const ccomplex_t IT_0432 = IT_0192*IT_0431;
    const ccomplex_t IT_0433 = -IT_0424 + -IT_0428 + IT_0432;
    const ccomplex_t IT_0434 = 2*IT_0433;
    const ccomplex_t IT_0435 = m_N_4*IT_0205;
    const ccomplex_t IT_0436 = IT_0207*(IT_0337 + -IT_0435);
    const ccomplex_t IT_0437 = -IT_0436;
    const ccomplex_t IT_0438 = 2*IT_0437;
    const ccomplex_t IT_0439 = IT_0434 + IT_0438;
    const ccomplex_t IT_0440 = 2*conj(IT_0433);
    const ccomplex_t IT_0441 = 2*conj(IT_0437);
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = m_N_4*IT_0259;
    const ccomplex_t IT_0444 = IT_0172*(IT_0365 + -IT_0443);
    const ccomplex_t IT_0445 = -IT_0444;
    const ccomplex_t IT_0446 = 2*IT_0445;
    const ccomplex_t IT_0447 = m_N_4*IT_0218;
    const ccomplex_t IT_0448 = -IT_0447;
    const ccomplex_t IT_0449 = IT_0360 + IT_0448;
    const ccomplex_t IT_0450 = IT_0192*IT_0449;
    const ccomplex_t IT_0451 = m_N_4*IT_0239;
    const ccomplex_t IT_0452 = -IT_0451;
    const ccomplex_t IT_0453 = IT_0376 + IT_0452;
    const ccomplex_t IT_0454 = IT_0152*IT_0453;
    const ccomplex_t IT_0455 = IT_0450 + -IT_0454;
    const ccomplex_t IT_0456 = 2*IT_0455;
    const ccomplex_t IT_0457 = m_N_4*IT_0268;
    const ccomplex_t IT_0458 = IT_0207*(IT_0370 + -IT_0457);
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = 2*IT_0459;
    const ccomplex_t IT_0461 = IT_0446 + IT_0456 + IT_0460;
    const ccomplex_t IT_0462 = 2*conj(IT_0445);
    const ccomplex_t IT_0463 = 2*conj(IT_0455);
    const ccomplex_t IT_0464 = 2*conj(IT_0459);
    const ccomplex_t IT_0465 = IT_0462 + IT_0463 + IT_0464;
    const ccomplex_t IT_0466 = s_13*s_34*IT_0126;
    const ccomplex_t IT_0467 = s_13*s_34*IT_0150;
    const ccomplex_t IT_0468 = (-2)*IT_0467;
    const ccomplex_t IT_0469 = IT_0392 + IT_0468;
    const ccomplex_t IT_0470 = IT_0081*IT_0469;
    const ccomplex_t IT_0471 = (-0.666666666666667)*IT_0466 + (
      -0.166666666666667)*IT_0470;
    const ccomplex_t IT_0472 = s_14 + IT_0471;
    const ccomplex_t IT_0473 = m_N_4*IT_0472;
    const ccomplex_t IT_0474 = (-6)*IT_0473;
    const ccomplex_t IT_0475 = (-2)*IT_0391;
    const ccomplex_t IT_0476 = IT_0390 + IT_0475;
    const ccomplex_t IT_0477 = IT_0126*IT_0476;
    const ccomplex_t IT_0478 = 0.166666666666667*IT_0470 + (-0.666666666666667
      )*IT_0477;
    const ccomplex_t IT_0479 = s_14 + IT_0478;
    const ccomplex_t IT_0480 = m_N_4*IT_0479;
    const ccomplex_t IT_0481 = 6*IT_0480;
    const ccomplex_t IT_0482 = IT_0081*IT_0315;
    const ccomplex_t IT_0483 = (-10)*s_24;
    const ccomplex_t IT_0484 = IT_0482 + IT_0483;
    const ccomplex_t IT_0485 = m_N_1*IT_0484;
    const ccomplex_t IT_0486 = s_12*s_14*IT_0150;
    const ccomplex_t IT_0487 = s_13*s_14*s_23;
    const ccomplex_t IT_0488 = s_24*IT_0048*IT_0150;
    const ccomplex_t IT_0489 = s_23*s_34*IT_0048;
    const ccomplex_t IT_0490 = (-2)*IT_0487 + -IT_0488 + 2*IT_0489;
    const ccomplex_t IT_0491 = IT_0486 + IT_0490;
    const ccomplex_t IT_0492 = IT_0126*IT_0491;
    const ccomplex_t IT_0493 = s_24*IT_0048;
    const ccomplex_t IT_0494 = 2*IT_0493;
    const ccomplex_t IT_0495 = s_12*s_14;
    const ccomplex_t IT_0496 = (-2)*IT_0495;
    const ccomplex_t IT_0497 = IT_0492 + IT_0494 + IT_0496;
    const ccomplex_t IT_0498 = (-2)*IT_0222;
    const ccomplex_t IT_0499 = (-2)*IT_0270;
    const ccomplex_t IT_0500 = IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = (-2)*conj(IT_0222);
    const ccomplex_t IT_0502 = (-2)*conj(IT_0270);
    const ccomplex_t IT_0503 = IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = (-2)*IT_0281;
    const ccomplex_t IT_0505 = (-2)*IT_0287;
    const ccomplex_t IT_0506 = IT_0504 + IT_0505;
    const ccomplex_t IT_0507 = (-2)*conj(IT_0281);
    const ccomplex_t IT_0508 = (-2)*conj(IT_0287);
    const ccomplex_t IT_0509 = IT_0507 + IT_0508;
    const ccomplex_t IT_0510 = s_14*s_24;
    const ccomplex_t IT_0511 = -IT_0510;
    const ccomplex_t IT_0512 = IT_0290 + IT_0511;
    const ccomplex_t IT_0513 = IT_0126*IT_0512;
    const ccomplex_t IT_0514 = -IT_0127 + (-0.25)*IT_0513;
    const ccomplex_t IT_0515 = s_12 + IT_0514;
    const ccomplex_t IT_0516 = m_N_1*IT_0515;
    const ccomplex_t IT_0517 = (-4)*IT_0516;
    const ccomplex_t IT_0518 = conj(IT_0433) + conj(IT_0437);
    const ccomplex_t IT_0519 = IT_0433 + IT_0437;
    const ccomplex_t IT_0520 = conj(IT_0445) + conj(IT_0455) + conj(IT_0459);
    const ccomplex_t IT_0521 = IT_0445 + IT_0455 + IT_0459;
    const ccomplex_t IT_0522 = s_34*IT_0081;
    const ccomplex_t IT_0523 = s_13*IT_0150;
    const ccomplex_t IT_0524 = IT_0522*IT_0523;
    const ccomplex_t IT_0525 = IT_0126*IT_0390;
    const ccomplex_t IT_0526 = (-0.5)*IT_0477 + (-0.5)*IT_0524 + 0.5*IT_0525;
    const ccomplex_t IT_0527 = s_14 + IT_0526;
    const ccomplex_t IT_0528 = m_N_4*IT_0527;
    const ccomplex_t IT_0529 = (-2)*IT_0528;
    const ccomplex_t IT_0530 = 2*IT_0123;
    const ccomplex_t IT_0531 = 2*conj(IT_0123);
    const ccomplex_t IT_0532 = 2*IT_0080;
    const ccomplex_t IT_0533 = 2*conj(IT_0080);
    const ccomplex_t IT_0534 = (conj(IT_0262) + conj(IT_0324) + conj(IT_0325) 
      + conj(IT_0326) + conj(IT_0403) + conj(IT_0404) + conj(IT_0405) + conj
      (IT_0406))*IT_0530 + (IT_0262 + IT_0324 + IT_0325 + IT_0326 + IT_0403 +
       IT_0404 + IT_0405 + IT_0406)*IT_0531 + (conj(IT_0283) + conj(IT_0356) +
       conj(IT_0357) + conj(IT_0398) + conj(IT_0399) + conj(IT_0400))*IT_0532 + 
      (IT_0283 + IT_0356 + IT_0357 + IT_0398 + IT_0399 + IT_0400)*IT_0533;
    const ccomplex_t IT_0535 = pow(s_23, 2);
    const ccomplex_t IT_0536 = IT_0126*IT_0535;
    const ccomplex_t IT_0537 = pow(s_24, 2);
    const ccomplex_t IT_0538 = IT_0126*IT_0537;
    const ccomplex_t IT_0539 = IT_0536 + IT_0538;
    const ccomplex_t IT_0540 = s_23*s_24*s_34*IT_0081;
    const ccomplex_t IT_0541 = -IT_0049 + -IT_0540;
    const ccomplex_t IT_0542 = IT_0539 + IT_0541;
    const ccomplex_t IT_0543 = m_N_1*IT_0542;
    const ccomplex_t IT_0544 = s_13*s_24*s_34*IT_0081;
    const ccomplex_t IT_0545 = s_13*s_23*IT_0126;
    const ccomplex_t IT_0546 = -IT_0127 + -IT_0545;
    const ccomplex_t IT_0547 = s_12 + IT_0544 + IT_0546;
    const ccomplex_t IT_0548 = m_N_4*IT_0547;
    const ccomplex_t IT_0549 = pow(s_14, 2);
    const ccomplex_t IT_0550 = IT_0048*IT_0150;
    const ccomplex_t IT_0551 = -IT_0550;
    const ccomplex_t IT_0552 = IT_0549 + IT_0551;
    const ccomplex_t IT_0553 = IT_0126*IT_0552;
    const ccomplex_t IT_0554 = IT_0126*IT_0549;
    const ccomplex_t IT_0555 = -IT_0048;
    const ccomplex_t IT_0556 = IT_0554 + IT_0555;
    const ccomplex_t IT_0557 = (-4)*IT_0556;
    const ccomplex_t IT_0558 = IT_0553 + IT_0557;
    const ccomplex_t IT_0559 = m_N_4*IT_0558;
    const ccomplex_t IT_0560 = conj(IT_0368) + conj(IT_0381);
    const ccomplex_t IT_0561 = IT_0368 + IT_0381;
    const ccomplex_t IT_0562 = conj(IT_0335) + conj(IT_0348);
    const ccomplex_t IT_0563 = IT_0335 + IT_0348;
    const ccomplex_t IT_0564 = IT_0559*(IT_0289*IT_0331 + IT_0288*IT_0352 +
       IT_0300*IT_0361 + IT_0297*conj(IT_0361) + conj(IT_0270)*IT_0362 + IT_0270
      *IT_0385 + IT_0500*IT_0560 + IT_0503*IT_0561 + IT_0506*IT_0562 + IT_0509
      *IT_0563);
    const ccomplex_t IT_0565 = IT_0081*IT_0392;
    const ccomplex_t IT_0566 = (-10)*s_14;
    const ccomplex_t IT_0567 = IT_0565 + IT_0566;
    const ccomplex_t IT_0568 = m_N_4*IT_0567;
    const ccomplex_t IT_0569 = IT_0369 + IT_0382;
    const ccomplex_t IT_0570 = IT_0386 + IT_0387;
    const ccomplex_t IT_0571 = conj(IT_0283) + conj(IT_0357);
    const ccomplex_t IT_0572 = IT_0283 + IT_0357;
    const ccomplex_t IT_0573 = IT_0049*IT_0150;
    const ccomplex_t IT_0574 = -IT_0573;
    const ccomplex_t IT_0575 = IT_0537 + IT_0574;
    const ccomplex_t IT_0576 = IT_0126*IT_0575;
    const ccomplex_t IT_0577 = -IT_0049;
    const ccomplex_t IT_0578 = IT_0538 + IT_0577;
    const ccomplex_t IT_0579 = (-4)*IT_0578;
    const ccomplex_t IT_0580 = IT_0576 + IT_0579;
    const ccomplex_t IT_0581 = m_N_1*IT_0580;
    const ccomplex_t IT_0582 = 2*IT_0209;
    const ccomplex_t IT_0583 = 2*conj(IT_0209);
    const ccomplex_t IT_0584 = 2*IT_0279;
    const ccomplex_t IT_0585 = 2*conj(IT_0279);
    const ccomplex_t IT_0586 = s_13*s_14*s_34*IT_0126;
    const ccomplex_t IT_0587 = 0.5*IT_0549 + (-0.5)*IT_0550 + -IT_0586;
    const ccomplex_t IT_0588 = pow(s_13, 2);
    const ccomplex_t IT_0589 = IT_0587 + IT_0588;
    const ccomplex_t IT_0590 = IT_0126*IT_0589;
    const ccomplex_t IT_0591 = -IT_0554;
    const ccomplex_t IT_0592 = IT_0048 + IT_0590 + IT_0591;
    const ccomplex_t IT_0593 = m_N_4*IT_0592;
    const ccomplex_t IT_0594 = 4*IT_0281;
    const ccomplex_t IT_0595 = 4*IT_0287;
    const ccomplex_t IT_0596 = 4*conj(IT_0281);
    const ccomplex_t IT_0597 = 4*conj(IT_0287);
    const ccomplex_t IT_0598 = 4*IT_0222;
    const ccomplex_t IT_0599 = 4*IT_0270;
    const ccomplex_t IT_0600 = 4*conj(IT_0222);
    const ccomplex_t IT_0601 = 4*conj(IT_0270);
    const ccomplex_t IT_0602 = s_14*IT_0049;
    const ccomplex_t IT_0603 = s_14*IT_0535;
    const ccomplex_t IT_0604 = s_12*s_24*IT_0150;
    const ccomplex_t IT_0605 = s_12*s_23*s_34;
    const ccomplex_t IT_0606 = s_14*IT_0049*IT_0150;
    const ccomplex_t IT_0607 = 0.5*IT_0604 + -IT_0605 + (-0.5)*IT_0606;
    const ccomplex_t IT_0608 = IT_0603 + IT_0607;
    const ccomplex_t IT_0609 = IT_0126*IT_0608;
    const ccomplex_t IT_0610 = s_12*s_24;
    const ccomplex_t IT_0611 = -IT_0610;
    const ccomplex_t IT_0612 = IT_0602 + IT_0609 + IT_0611;
    const ccomplex_t IT_0613 = 4*IT_0209;
    const ccomplex_t IT_0614 = 4*conj(IT_0209);
    const ccomplex_t IT_0615 = 4*IT_0279;
    const ccomplex_t IT_0616 = 4*conj(IT_0279);
    const ccomplex_t IT_0617 = s_13*s_23*s_24;
    const ccomplex_t IT_0618 = s_13*s_34*IT_0049;
    const ccomplex_t IT_0619 = -IT_0606 + (-2)*IT_0617 + 2*IT_0618;
    const ccomplex_t IT_0620 = IT_0604 + IT_0619;
    const ccomplex_t IT_0621 = IT_0126*IT_0620;
    const ccomplex_t IT_0622 = (-0.5)*IT_0621;
    const ccomplex_t IT_0623 = -IT_0602;
    const ccomplex_t IT_0624 = IT_0610 + IT_0622 + IT_0623;
    const ccomplex_t IT_0625 = IT_0124*IT_0556;
    const ccomplex_t IT_0626 = 6*conj(IT_0222);
    const ccomplex_t IT_0627 = (-2)*IT_0510;
    const ccomplex_t IT_0628 = IT_0290 + IT_0627;
    const ccomplex_t IT_0629 = IT_0126*IT_0628;
    const ccomplex_t IT_0630 = (-0.5)*IT_0629;
    const ccomplex_t IT_0631 = (-0.5)*IT_0294;
    const ccomplex_t IT_0632 = s_12*IT_0314;
    const ccomplex_t IT_0633 = s_13*s_23*IT_0150;
    const ccomplex_t IT_0634 = s_14*s_24*IT_0150;
    const ccomplex_t IT_0635 = s_13*s_24*s_34;
    const ccomplex_t IT_0636 = (-2)*IT_0633 + (-2)*IT_0634 + 4*IT_0635;
    const ccomplex_t IT_0637 = IT_0632 + IT_0636;
    const ccomplex_t IT_0638 = IT_0081*IT_0637;
    const ccomplex_t IT_0639 = 0.25*IT_0638;
    const ccomplex_t IT_0640 = 8*IT_0433;
    const ccomplex_t IT_0641 = 8*conj(IT_0433);
    const ccomplex_t IT_0642 = 8*IT_0437;
    const ccomplex_t IT_0643 = 8*IT_0445;
    const ccomplex_t IT_0644 = 8*conj(IT_0445);
    const ccomplex_t IT_0645 = 8*IT_0455;
    const ccomplex_t IT_0646 = 8*conj(IT_0455);
    const ccomplex_t IT_0647 = (s_12 + IT_0630 + IT_0631 + IT_0639)*(IT_0518
      *IT_0640 + IT_0437*IT_0641 + conj(IT_0437)*IT_0642 + IT_0520*IT_0643 +
       IT_0455*IT_0644 + (conj(IT_0455) + conj(IT_0459))*IT_0645 + 8*IT_0459*
      (conj(IT_0459) + 0.125*IT_0644 + 0.125*IT_0646));
    const ccomplex_t IT_0648 = 1.125*IT_0124;
    const ccomplex_t IT_0649 = 16*conj(IT_0361);
    const ccomplex_t IT_0650 = (-16)*IT_0335;
    const ccomplex_t IT_0651 = (-16)*IT_0348;
    const ccomplex_t IT_0652 = 16*conj(IT_0330);
    const ccomplex_t IT_0653 = (-16)*conj(IT_0335);
    const ccomplex_t IT_0654 = (-16)*conj(IT_0348);
    const ccomplex_t IT_0655 = (-16)*IT_0433;
    const ccomplex_t IT_0656 = (-16)*IT_0437;
    const ccomplex_t IT_0657 = (-16)*conj(IT_0433);
    const ccomplex_t IT_0658 = (-16)*conj(IT_0437);
    const ccomplex_t IT_0659 = ((-16)*conj(IT_0330) + 16*conj(IT_0335) + 16
      *conj(IT_0348))*IT_0361 + IT_0330*((-16)*conj(IT_0361) + 16*conj(IT_0368) 
      + 16*conj(IT_0381)) + IT_0563*IT_0649 + IT_0560*(IT_0650 + IT_0651) +
       IT_0561*(IT_0652 + IT_0653 + IT_0654) + IT_0520*(IT_0655 + IT_0656) +
       IT_0521*(IT_0657 + IT_0658);
    const ccomplex_t IT_0660 = IT_0126*IT_0290;
    const ccomplex_t IT_0661 = (-0.5)*IT_0660;
    const ccomplex_t IT_0662 = s_14*s_23;
    const ccomplex_t IT_0663 = s_12*s_34;
    const ccomplex_t IT_0664 = IT_0662 + IT_0663;
    const ccomplex_t IT_0665 = s_13*s_24;
    const ccomplex_t IT_0666 = -IT_0665;
    const ccomplex_t IT_0667 = IT_0664 + IT_0666;
    const ccomplex_t IT_0668 = IT_0522*IT_0667;
    const ccomplex_t IT_0669 = 0.25*IT_0668;
    const ccomplex_t IT_0670 = s_12 + IT_0661 + IT_0669;
    const ccomplex_t IT_0671 = (-8)*IT_0123;
    const ccomplex_t IT_0672 = (-8)*conj(IT_0123);
    const ccomplex_t IT_0673 = (-8)*IT_0080;
    const ccomplex_t IT_0674 = (-8)*conj(IT_0080);
    const ccomplex_t IT_0675 = 8*IT_0670*(conj(IT_0080)*IT_0330 + IT_0080*conj
      (IT_0330) + conj(IT_0123)*IT_0361 + IT_0123*conj(IT_0361) + 0.125*IT_0560
      *IT_0671 + 0.125*IT_0561*IT_0672 + 0.125*IT_0562*IT_0673 + 0.125*IT_0563
      *IT_0674);
    const ccomplex_t IT_0676 = -IT_0663 + -IT_0665;
    const ccomplex_t IT_0677 = IT_0662 + IT_0676;
    const ccomplex_t IT_0678 = IT_0522*IT_0677;
    const ccomplex_t IT_0679 = (-0.25)*IT_0678;
    const ccomplex_t IT_0680 = s_12 + IT_0661 + IT_0679;
    const ccomplex_t IT_0681 = (-2)*IT_0082;
    const ccomplex_t IT_0682 = IT_0314 + IT_0681;
    const ccomplex_t IT_0683 = IT_0081*IT_0682;
    const ccomplex_t IT_0684 = 0.25*IT_0683;
    const ccomplex_t IT_0685 = 1 + IT_0684;
    const ccomplex_t IT_0686 = IT_0124*IT_0685;
    const ccomplex_t IT_0687 = 8*IT_0335;
    const ccomplex_t IT_0688 = 8*IT_0348;
    const ccomplex_t IT_0689 = (-8)*IT_0330;
    const ccomplex_t IT_0690 = IT_0687 + IT_0688 + IT_0689;
    const ccomplex_t IT_0691 = 8*conj(IT_0437);
    const ccomplex_t IT_0692 = 8*conj(IT_0335);
    const ccomplex_t IT_0693 = 8*conj(IT_0348);
    const ccomplex_t IT_0694 = (-8)*conj(IT_0330);
    const ccomplex_t IT_0695 = (-8)*IT_0361;
    const ccomplex_t IT_0696 = (-8)*conj(IT_0361);
    const ccomplex_t IT_0697 = s_14*s_23*s_34;
    const ccomplex_t IT_0698 = s_12*IT_0082;
    const ccomplex_t IT_0699 = (-2)*IT_0633 + (-2)*IT_0634 + 2*IT_0635 + 2
      *IT_0697 + (-2)*IT_0698;
    const ccomplex_t IT_0700 = IT_0632 + IT_0699;
    const ccomplex_t IT_0701 = IT_0081*IT_0700;
    const ccomplex_t IT_0702 = (-0.25)*IT_0701;
    const ccomplex_t IT_0703 = -IT_0545;
    const ccomplex_t IT_0704 = IT_0687 + IT_0689;
    const ccomplex_t IT_0705 = IT_0692 + IT_0694;
    const ccomplex_t IT_0706 = 8*IT_0368;
    const ccomplex_t IT_0707 = IT_0695 + IT_0706;
    const ccomplex_t IT_0708 = 8*conj(IT_0368);
    const ccomplex_t IT_0709 = 8*IT_0381;
    const ccomplex_t IT_0710 = conj(IT_0348)*IT_0640 + IT_0348*IT_0641 + conj
      (IT_0368)*IT_0643 + IT_0561*IT_0644 + 8*conj(IT_0381)*(IT_0459 + 0.125
      *IT_0643 + 0.125*IT_0645) + IT_0381*IT_0646 + conj(IT_0437)*IT_0690 +
       IT_0437*(IT_0692 + IT_0693 + IT_0694) + conj(IT_0445)*IT_0695 + IT_0521
      *IT_0696 + conj(IT_0433)*IT_0704 + IT_0433*IT_0705 + conj(IT_0455)*IT_0707
       + (IT_0455 + IT_0459)*IT_0708 + conj(IT_0459)*(IT_0695 + IT_0706 +
       IT_0709);
    const ccomplex_t IT_0711 = (-2)*IT_0633 + (-2)*IT_0634 + 4*IT_0697;
    const ccomplex_t IT_0712 = IT_0632 + IT_0711;
    const ccomplex_t IT_0713 = IT_0081*IT_0712;
    const ccomplex_t IT_0714 = 0.25*IT_0713;
    const ccomplex_t IT_0715 = IT_0330*conj(IT_0330) + IT_0361*conj(IT_0361) +
       0.125*conj(IT_0348)*IT_0690 + 0.125*IT_0335*IT_0694 + 0.125*IT_0368
      *IT_0696 + 0.125*conj(IT_0335)*IT_0704 + 0.125*IT_0348*IT_0705 + 0.125
      *conj(IT_0368)*IT_0707 + 0.125*IT_0381*(IT_0696 + IT_0708) + 0.125*conj
      (IT_0381)*(IT_0695 + IT_0706 + IT_0709);
    const ccomplex_t IT_0716 = pow(m_Z, 6);
    const ccomplex_t IT_0717 = s_12*IT_0716;
    const ccomplex_t IT_0718 = s_13*s_23*IT_0314;
    const ccomplex_t IT_0719 = (-2)*IT_0718;
    const ccomplex_t IT_0720 = IT_0717 + IT_0719;
    const ccomplex_t IT_0721 = IT_0081*IT_0720;
    const ccomplex_t IT_0722 = 0.25*IT_0721;
    const ccomplex_t IT_0723 = IT_0290 + IT_0722;
    const ccomplex_t IT_0724 = (-8)*IT_0356;
    const ccomplex_t IT_0725 = (-8)*conj(IT_0356);
    const ccomplex_t IT_0726 = conj(IT_0325) + conj(IT_0326);
    const ccomplex_t IT_0727 = (-8)*IT_0403;
    const ccomplex_t IT_0728 = IT_0325 + IT_0326;
    const ccomplex_t IT_0729 = (-8)*conj(IT_0403);
    const ccomplex_t IT_0730 = (-8)*IT_0324;
    const ccomplex_t IT_0731 = (-8)*IT_0262;
    const ccomplex_t IT_0732 = (-8)*conj(IT_0324);
    const ccomplex_t IT_0733 = (-8)*conj(IT_0262);
    const ccomplex_t IT_0734 = (-8)*IT_0398;
    const ccomplex_t IT_0735 = (-8)*conj(IT_0398);
    const ccomplex_t IT_0736 = conj(IT_0399) + conj(IT_0400);
    const ccomplex_t IT_0737 = (-8)*IT_0283;
    const ccomplex_t IT_0738 = IT_0399 + IT_0400;
    const ccomplex_t IT_0739 = (-8)*conj(IT_0283);
    const ccomplex_t IT_0740 = conj(IT_0262) + conj(IT_0325) + conj(IT_0326);
    const ccomplex_t IT_0741 = (-8)*IT_0405;
    const ccomplex_t IT_0742 = IT_0262 + IT_0325 + IT_0326;
    const ccomplex_t IT_0743 = (-8)*conj(IT_0405);
    const ccomplex_t IT_0744 = conj(IT_0404) + conj(IT_0406);
    const ccomplex_t IT_0745 = (-8)*IT_0325;
    const ccomplex_t IT_0746 = IT_0404 + IT_0406;
    const ccomplex_t IT_0747 = (-8)*conj(IT_0325);
    const ccomplex_t IT_0748 = conj(IT_0398)*IT_0724 + IT_0398*IT_0725 +
       IT_0726*IT_0727 + IT_0728*IT_0729 + conj(IT_0405)*IT_0730 + conj(IT_0403)
      *(IT_0730 + IT_0731) + IT_0405*IT_0732 + IT_0403*(IT_0732 + IT_0733) +
       IT_0571*IT_0734 + IT_0572*IT_0735 + (-8)*IT_0736*(IT_0357 + (-0.125)
      *IT_0724 + (-0.125)*IT_0737) + (-8)*IT_0738*(conj(IT_0357) + (-0.125)
      *IT_0725 + (-0.125)*IT_0739) + IT_0740*IT_0741 + IT_0742*IT_0743 + (-8)
      *IT_0744*(IT_0326 + (-0.125)*IT_0730 + (-0.125)*IT_0731 + (-0.125)*IT_0745
      ) + (-8)*IT_0746*(conj(IT_0326) + (-0.125)*IT_0732 + (-0.125)*IT_0733 + (
      -0.125)*IT_0747);
    const ccomplex_t IT_0749 = s_24*IT_0588;
    const ccomplex_t IT_0750 = s_12*s_13*s_34;
    const ccomplex_t IT_0751 = -IT_0488 + 2*IT_0749 + (-2)*IT_0750;
    const ccomplex_t IT_0752 = IT_0486 + IT_0751;
    const ccomplex_t IT_0753 = IT_0126*IT_0752;
    const ccomplex_t IT_0754 = (-0.5)*IT_0753;
    const ccomplex_t IT_0755 = -IT_0493;
    const ccomplex_t IT_0756 = IT_0495 + IT_0754 + IT_0755;
    const ccomplex_t IT_0757 = (-4)*IT_0281;
    const ccomplex_t IT_0758 = (-4)*IT_0287;
    const ccomplex_t IT_0759 = (-4)*conj(IT_0281);
    const ccomplex_t IT_0760 = (-4)*conj(IT_0287);
    const ccomplex_t IT_0761 = (-4)*IT_0222;
    const ccomplex_t IT_0762 = (-4)*IT_0270;
    const ccomplex_t IT_0763 = (-4)*conj(IT_0222);
    const ccomplex_t IT_0764 = (-4)*conj(IT_0270);
    const ccomplex_t IT_0765 = s_23*s_24*s_34*IT_0126;
    const ccomplex_t IT_0766 = 0.5*IT_0537 + (-0.5)*IT_0573 + -IT_0765;
    const ccomplex_t IT_0767 = IT_0535 + IT_0766;
    const ccomplex_t IT_0768 = IT_0126*IT_0767;
    const ccomplex_t IT_0769 = -IT_0538;
    const ccomplex_t IT_0770 = IT_0049 + IT_0768 + IT_0769;
    const ccomplex_t IT_0771 = m_N_1*IT_0770;
    const ccomplex_t IT_0772 = (-4)*IT_0209;
    const ccomplex_t IT_0773 = (-4)*conj(IT_0209);
    const ccomplex_t IT_0774 = (-4)*IT_0279;
    const ccomplex_t IT_0775 = (-4)*conj(IT_0279);
    const ccomplex_t IT_0776 = 4*IT_0771*(conj(IT_0279)*IT_0330 + IT_0279*conj
      (IT_0330) + conj(IT_0209)*IT_0361 + IT_0209*conj(IT_0361) + 0.25*IT_0560
      *IT_0772 + 0.25*IT_0561*IT_0773 + 0.25*IT_0562*IT_0774 + 0.25*IT_0563
      *IT_0775);
    const ccomplex_t IT_0777 = 0.5*IT_0721;
    const ccomplex_t IT_0778 = (-4)*IT_0510;
    const ccomplex_t IT_0779 = IT_0290 + IT_0777 + IT_0778;
    const ccomplex_t IT_0780 = IT_0126*IT_0637;
    const ccomplex_t IT_0781 = -IT_0780;
    const ccomplex_t IT_0782 = IT_0779 + IT_0781;
    const ccomplex_t IT_0783 = (-4)*conj(IT_0403);
    const ccomplex_t IT_0784 = (-4)*IT_0398;
    const ccomplex_t IT_0785 = (-4)*conj(IT_0398);
    const ccomplex_t IT_0786 = (-4)*IT_0399;
    const ccomplex_t IT_0787 = (-4)*IT_0405;
    const ccomplex_t IT_0788 = (-4)*IT_0403;
    const ccomplex_t IT_0789 = (-4)*IT_0404;
    const ccomplex_t IT_0790 = (-4)*conj(IT_0405);
    const ccomplex_t IT_0791 = IT_0746*IT_0783 + conj(IT_0398)*IT_0784 +
       IT_0738*IT_0785 + (-4)*IT_0736*(IT_0400 + (-0.25)*IT_0784 + (-0.25)
      *IT_0786) + conj(IT_0405)*IT_0787 + conj(IT_0403)*(IT_0787 + IT_0788) + (
      -4)*IT_0744*(IT_0406 + (-0.25)*IT_0787 + (-0.25)*IT_0788 + (-0.25)*IT_0789
      ) + (IT_0406 + IT_0408)*IT_0790;
    const ccomplex_t IT_0792 = IT_0126*IT_0712;
    const ccomplex_t IT_0793 = -IT_0792;
    const ccomplex_t IT_0794 = IT_0779 + IT_0793;
    const ccomplex_t IT_0795 = (-4)*IT_0356;
    const ccomplex_t IT_0796 = (-4)*conj(IT_0356);
    const ccomplex_t IT_0797 = (-4)*conj(IT_0262);
    const ccomplex_t IT_0798 = (-4)*IT_0324;
    const ccomplex_t IT_0799 = (-4)*IT_0262;
    const ccomplex_t IT_0800 = (-4)*conj(IT_0324);
    const ccomplex_t IT_0801 = (-4)*IT_0283;
    const ccomplex_t IT_0802 = (-4)*IT_0325;
    const ccomplex_t IT_0803 = conj(IT_0356)*IT_0795 + IT_0572*IT_0796 +
       IT_0728*IT_0797 + conj(IT_0324)*IT_0798 + conj(IT_0262)*(IT_0798 +
       IT_0799) + IT_0742*IT_0800 + (-4)*IT_0571*(IT_0357 + (-0.25)*IT_0795 + (
      -0.25)*IT_0801) + (-4)*IT_0726*(IT_0326 + (-0.25)*IT_0798 + (-0.25)
      *IT_0799 + (-0.25)*IT_0802);
    const ccomplex_t IT_0804 = s_23*IT_0150;
    const ccomplex_t IT_0805 = IT_0522*IT_0804;
    const ccomplex_t IT_0806 = IT_0126*IT_0309;
    const ccomplex_t IT_0807 = (-0.5)*IT_0313 + (-0.5)*IT_0805 + 0.5*IT_0806;
    const ccomplex_t IT_0808 = s_24 + IT_0807;
    const ccomplex_t IT_0809 = m_N_1*IT_0808;
    const ccomplex_t IT_0810 = (-2)*IT_0809;
    const ccomplex_t IT_0811 = (-2)*IT_0123;
    const ccomplex_t IT_0812 = (-2)*conj(IT_0123);
    const ccomplex_t IT_0813 = (-2)*IT_0080;
    const ccomplex_t IT_0814 = (-2)*conj(IT_0080);
    const ccomplex_t IT_0815 = (conj(IT_0283) + conj(IT_0356) + conj(IT_0357) 
      + conj(IT_0398) + conj(IT_0399) + conj(IT_0400))*IT_0811 + (IT_0283 +
       IT_0356 + IT_0357 + IT_0398 + IT_0399 + IT_0400)*IT_0812 + (conj(IT_0262)
       + conj(IT_0324) + conj(IT_0325) + conj(IT_0326) + conj(IT_0403) + conj
      (IT_0404) + conj(IT_0405) + conj(IT_0406))*IT_0813 + (IT_0262 + IT_0324 +
       IT_0325 + IT_0326 + IT_0403 + IT_0404 + IT_0405 + IT_0406)*IT_0814;
    const ccomplex_t IT_0816 = s_14*s_23*s_34*IT_0081;
    const ccomplex_t IT_0817 = s_12 + IT_0546 + IT_0816;
    const ccomplex_t IT_0818 = m_N_1*IT_0817;
    const ccomplex_t IT_0819 = IT_0126*IT_0588;
    const ccomplex_t IT_0820 = IT_0554 + IT_0819;
    const ccomplex_t IT_0821 = s_13*s_14*s_34*IT_0081;
    const ccomplex_t IT_0822 = -IT_0048 + -IT_0821;
    const ccomplex_t IT_0823 = IT_0820 + IT_0822;
    const ccomplex_t IT_0824 = m_N_4*IT_0823;
    const ccomplex_t IT_0825 = -(s_12 + (-0.5)*IT_0294)*(IT_0048 + -IT_0554);
    const ccomplex_t IT_0826 = (-2)*IT_0825;
    const ccomplex_t IT_0827 = s_14*s_23*s_34*IT_0126;
    const ccomplex_t IT_0828 = (-2)*IT_0291 + -IT_0510 + 2*IT_0827;
    const ccomplex_t IT_0829 = IT_0290 + IT_0828;
    const ccomplex_t IT_0830 = IT_0126*IT_0829;
    const ccomplex_t IT_0831 = -IT_0127 + (-0.5)*IT_0830;
    const ccomplex_t IT_0832 = s_12 + IT_0831;
    const ccomplex_t IT_0833 = m_N_1*IT_0832;
    const ccomplex_t IT_0834 = (-2)*IT_0833;
    const ccomplex_t IT_0835 = (-2)*IT_0330;
    const ccomplex_t IT_0836 = (-2)*conj(IT_0330);
    const ccomplex_t IT_0837 = (-2)*IT_0361;
    const ccomplex_t IT_0838 = (-2)*conj(IT_0361);
    const ccomplex_t IT_0839 = conj(IT_0330)*IT_0498 + IT_0330*IT_0501 + conj
      (IT_0361)*IT_0504 + IT_0361*IT_0507 + IT_0305*IT_0560 + IT_0308*IT_0561 +
       IT_0299*IT_0562 + IT_0302*IT_0563 + conj(IT_0270)*IT_0835 + IT_0270
      *IT_0836 + conj(IT_0287)*IT_0837 + IT_0287*IT_0838;
    const ccomplex_t IT_0840 = m_N_4*IT_0515;
    const ccomplex_t IT_0841 = (-4)*IT_0840;
    const ccomplex_t IT_0842 = 0.666666666666667*IT_0525 + (-0.166666666666667
      )*IT_0565;
    const ccomplex_t IT_0843 = s_14 + IT_0842;
    const ccomplex_t IT_0844 = m_N_4*IT_0843;
    const ccomplex_t IT_0845 = (-6)*IT_0844;
    const ccomplex_t IT_0846 = s_13*s_24*s_34*IT_0126;
    const ccomplex_t IT_0847 = (-2)*IT_0291 + -IT_0510 + 2*IT_0846;
    const ccomplex_t IT_0848 = IT_0290 + IT_0847;
    const ccomplex_t IT_0849 = IT_0126*IT_0848;
    const ccomplex_t IT_0850 = -IT_0127 + (-0.5)*IT_0849;
    const ccomplex_t IT_0851 = s_12 + IT_0850;
    const ccomplex_t IT_0852 = m_N_4*IT_0851;
    const ccomplex_t IT_0853 = (-2)*IT_0852;
    const ccomplex_t IT_0854 = (-2)*IT_0209;
    const ccomplex_t IT_0855 = (-2)*conj(IT_0209);
    const ccomplex_t IT_0856 = (-2)*IT_0279;
    const ccomplex_t IT_0857 = (-2)*conj(IT_0279);
    const ccomplex_t IT_0858 = IT_0124*IT_0150;
    const ccomplex_t IT_0859 = 9*IT_0858;
    const ccomplex_t IT_0860 = (-2)*IT_0356;
    const ccomplex_t IT_0861 = (-2)*IT_0283;
    const ccomplex_t IT_0862 = (-2)*IT_0357;
    const ccomplex_t IT_0863 = IT_0860 + IT_0861 + IT_0862;
    const ccomplex_t IT_0864 = (-2)*conj(IT_0356);
    const ccomplex_t IT_0865 = (-2)*conj(IT_0283);
    const ccomplex_t IT_0866 = (-2)*conj(IT_0357);
    const ccomplex_t IT_0867 = IT_0864 + IT_0865 + IT_0866;
    const ccomplex_t IT_0868 = (-2)*IT_0398;
    const ccomplex_t IT_0869 = (-2)*IT_0399;
    const ccomplex_t IT_0870 = (-2)*IT_0400;
    const ccomplex_t IT_0871 = IT_0868 + IT_0869 + IT_0870;
    const ccomplex_t IT_0872 = (-2)*conj(IT_0398);
    const ccomplex_t IT_0873 = (-2)*conj(IT_0399);
    const ccomplex_t IT_0874 = (-2)*conj(IT_0400);
    const ccomplex_t IT_0875 = IT_0872 + IT_0873 + IT_0874;
    const ccomplex_t IT_0876 = (-2)*IT_0324;
    const ccomplex_t IT_0877 = (-2)*conj(IT_0324);
    const ccomplex_t IT_0878 = m_N_1*m_N_4*IT_0126;
    const ccomplex_t IT_0879 = (-0.25)*IT_0314;
    const ccomplex_t IT_0880 = IT_0082 + IT_0879;
    const ccomplex_t IT_0881 = IT_0878*IT_0880;
    const ccomplex_t IT_0882 = (-4)*IT_0881;
    const ccomplex_t IT_0883 = (-4)*IT_0125;
    const ccomplex_t IT_0884 = -IT_0123;
    const ccomplex_t IT_0885 = -conj(IT_0123);
    const ccomplex_t IT_0886 = -IT_0080;
    const ccomplex_t IT_0887 = -conj(IT_0080);
    const ccomplex_t IT_0888 = conj(IT_0123)*IT_0330 + IT_0123*conj(IT_0330) +
       conj(IT_0080)*IT_0361 + IT_0080*conj(IT_0361) + (conj(IT_0433) + conj
      (IT_0437) + IT_0562)*IT_0884 + (IT_0433 + IT_0437 + IT_0563)*IT_0885 + 
      (conj(IT_0445) + conj(IT_0455) + conj(IT_0459) + IT_0560)*IT_0886 + 
      (IT_0445 + IT_0455 + IT_0459 + IT_0561)*IT_0887;
    const ccomplex_t IT_0889 = IT_0124*IT_0578;
    const ccomplex_t IT_0890 = 0.166666666666667*IT_0543;
    const ccomplex_t IT_0891 = -(s_12 + (-0.5)*IT_0294)*(IT_0049 + -IT_0538);
    const ccomplex_t IT_0892 = (-2)*IT_0891;
    const ccomplex_t IT_0893 = 0.166666666666667*IT_0892;
    const ccomplex_t IT_0894 = 4*IT_0123*(IT_0086*conj(IT_0123) + -conj
      (IT_0080)*IT_0125) + 4*IT_0080*(conj(IT_0080)*IT_0086 + -conj(IT_0123)
      *IT_0125) + (-6)*IT_0130*(conj(IT_0209)*IT_0271 + IT_0209*IT_0272 + conj
      (IT_0279)*IT_0288 + IT_0279*IT_0289) + IT_0296*(conj(IT_0279)*IT_0299 +
       IT_0279*IT_0302 + conj(IT_0209)*IT_0305 + IT_0209*IT_0308) + IT_0323*
      (IT_0327*IT_0350 + IT_0351*IT_0355 + IT_0358*IT_0383 + IT_0384*IT_0388) +
       IT_0397*(IT_0350*IT_0401 + IT_0355*IT_0402 + IT_0383*IT_0407 + IT_0388
      *IT_0409) + (IT_0383*IT_0401 + IT_0388*IT_0402 + IT_0350*IT_0407 + IT_0355
      *IT_0409)*IT_0414 + IT_0420*(IT_0327*IT_0439 + IT_0351*IT_0442 + IT_0358
      *IT_0461 + IT_0384*IT_0465) + (IT_0358*IT_0439 + IT_0384*IT_0442 + IT_0327
      *IT_0461 + IT_0351*IT_0465)*IT_0474 + (IT_0401*IT_0439 + IT_0402*IT_0442 +
       IT_0407*IT_0461 + IT_0409*IT_0465)*IT_0481 + (IT_0407*IT_0439 + IT_0409
      *IT_0442 + IT_0401*IT_0461 + IT_0402*IT_0465)*IT_0485 + IT_0497*(IT_0327
      *IT_0500 + IT_0351*IT_0503 + IT_0358*IT_0506 + IT_0384*IT_0509) + IT_0517*
      (IT_0500*IT_0518 + IT_0503*IT_0519 + IT_0506*IT_0520 + IT_0509*IT_0521) +
       IT_0529*IT_0534 + IT_0279*(IT_0533*IT_0543 + IT_0531*IT_0548) + IT_0209*
      (IT_0531*IT_0543 + IT_0533*IT_0548) + IT_0564 + IT_0568*(IT_0355*IT_0356 +
       IT_0350*conj(IT_0356) + IT_0327*IT_0362 + IT_0351*IT_0385 + conj(IT_0324)
      *IT_0569 + IT_0324*IT_0570 + IT_0331*IT_0571 + IT_0352*IT_0572) + IT_0581*
      (IT_0520*IT_0582 + IT_0521*IT_0583 + IT_0518*IT_0584 + IT_0519*IT_0585) +
       IT_0593*(IT_0518*(IT_0594 + IT_0595) + IT_0519*(IT_0596 + IT_0597) +
       IT_0520*(IT_0598 + IT_0599) + IT_0521*(IT_0600 + IT_0601)) + IT_0612*
      (IT_0358*IT_0613 + IT_0384*IT_0614 + IT_0327*IT_0615 + IT_0351*IT_0616) + 
      (IT_0401*IT_0613 + IT_0402*IT_0614 + IT_0407*IT_0615 + IT_0409*IT_0616)
      *IT_0624 + 6*IT_0625*(IT_0271*IT_0289 + IT_0288*(conj(IT_0270) +
       0.166666666666667*IT_0626)) + IT_0647 + IT_0648*IT_0659 + IT_0675 + 
      (IT_0520*IT_0671 + IT_0521*IT_0672 + IT_0518*IT_0673 + IT_0519*IT_0674)
      *IT_0680 + IT_0686*(IT_0560*(IT_0640 + IT_0642) + IT_0520*IT_0690 +
       IT_0561*(IT_0641 + IT_0691) + IT_0521*(IT_0692 + IT_0693 + IT_0694) +
       IT_0518*IT_0695 + IT_0519*IT_0696) + (s_12 + IT_0128 + IT_0702 + IT_0703)
      *IT_0710 + 8*(s_12 + IT_0630 + IT_0631 + IT_0714)*IT_0715 + IT_0723
      *IT_0748 + IT_0756*(IT_0401*(IT_0757 + IT_0758) + IT_0402*(IT_0759 +
       IT_0760) + IT_0407*(IT_0761 + IT_0762) + IT_0409*(IT_0763 + IT_0764)) +
       IT_0776 + IT_0782*IT_0791 + IT_0794*IT_0803 + IT_0810*IT_0815 + (IT_0289
      *IT_0811 + IT_0288*IT_0812 + IT_0272*IT_0813 + IT_0271*IT_0814)*IT_0818 + 
      (IT_0272*IT_0811 + IT_0271*IT_0812 + IT_0289*IT_0813 + IT_0288*IT_0814)
      *IT_0824 + (IT_0272*IT_0498 + conj(IT_0270)*IT_0499 + IT_0270*IT_0501 +
       IT_0289*IT_0504 + conj(IT_0287)*IT_0505 + IT_0287*IT_0507)*IT_0826 +
       IT_0834*IT_0839 + (IT_0562*IT_0582 + IT_0563*IT_0583 + IT_0560*IT_0584 +
       IT_0561*IT_0585 + conj(IT_0209)*IT_0835 + IT_0209*IT_0836 + conj(IT_0279)
      *IT_0837 + IT_0279*IT_0838)*IT_0841 + ((IT_0336 + IT_0349)*IT_0571 + 
      (IT_0353 + IT_0354)*IT_0572 + IT_0569*IT_0740 + IT_0570*IT_0742)*IT_0845 +
       IT_0853*(IT_0518*IT_0854 + IT_0519*IT_0855 + IT_0520*IT_0856 + IT_0521
      *IT_0857) + IT_0859*(IT_0740*IT_0863 + IT_0742*IT_0867 + IT_0407*IT_0871 +
       IT_0409*IT_0875 + IT_0358*IT_0876 + IT_0384*IT_0877) + (IT_0407*IT_0863 +
       IT_0409*IT_0867 + IT_0740*IT_0871 + IT_0742*IT_0875 + IT_0401*IT_0876 +
       IT_0402*IT_0877)*IT_0882 + IT_0883*IT_0888 + 6*conj(IT_0209)*
      (0.166666666666667*IT_0532*IT_0548 + IT_0279*IT_0889 + IT_0530*IT_0890 +
       IT_0854*IT_0893) + 6*conj(IT_0279)*(0.166666666666667*IT_0530*IT_0548 +
       IT_0209*IT_0889 + IT_0532*IT_0890 + IT_0856*IT_0893);
    return create_ccomplex_return(IT_0894);
}

