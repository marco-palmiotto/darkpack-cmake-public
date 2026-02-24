#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_h_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0000 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0016 = IT_0007*IT_0012;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0005 + IT_0009 + 
      -IT_0014 + -IT_0017);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0021 = IT_0004*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0027 = IT_0016*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + IT_0023 + 
      -IT_0025 + -IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0019*IT_0029;
    const ccomplex_t IT_0031 = pow(m_h, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0031 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0031 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = IT_0030*IT_0037;
    const ccomplex_t IT_0039 = pow(m_N_1, 2);
    const ccomplex_t IT_0040 = IT_0031*IT_0039;
    const ccomplex_t IT_0041 = s_12*IT_0031;
    const ccomplex_t IT_0042 = s_14*s_24;
    const ccomplex_t IT_0043 = (-2)*IT_0042;
    const ccomplex_t IT_0044 = IT_0041 + IT_0043;
    const ccomplex_t IT_0045 = IT_0030*IT_0035;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0031 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0051 = IT_0008*IT_0050;
    const ccomplex_t IT_0052 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0053 = IT_0004*IT_0052;
    const ccomplex_t IT_0054 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0055 = IT_0008*IT_0054;
    const ccomplex_t IT_0056 = IT_0049 + IT_0051 + IT_0053 + IT_0055;
    const ccomplex_t IT_0057 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0058 = IT_0013*IT_0057;
    const ccomplex_t IT_0059 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0060 = IT_0016*IT_0059;
    const ccomplex_t IT_0061 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0062 = IT_0013*IT_0061;
    const ccomplex_t IT_0063 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0064 = IT_0016*IT_0063;
    const ccomplex_t IT_0065 = -IT_0058 + -IT_0060 + -IT_0062 + -IT_0064;
    const ccomplex_t IT_0066 = IT_0056 + IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0070 = IT_0004*IT_0069;
    const ccomplex_t IT_0071 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0072 = IT_0008*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0074 = IT_0004*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0076 = IT_0008*IT_0075;
    const ccomplex_t IT_0077 = IT_0070 + IT_0072 + IT_0074 + IT_0076;
    const ccomplex_t IT_0078 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0079 = IT_0013*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0081 = IT_0016*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0083 = IT_0013*IT_0082;
    const ccomplex_t IT_0084 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0085 = IT_0016*IT_0084;
    const ccomplex_t IT_0086 = -IT_0079 + -IT_0081 + -IT_0083 + -IT_0085;
    const ccomplex_t IT_0087 = IT_0077 + IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0068*IT_0089;
    const ccomplex_t IT_0091 = IT_0047*IT_0090;
    const ccomplex_t IT_0092 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0093 = IT_0004*IT_0092;
    const ccomplex_t IT_0094 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0095 = IT_0008*IT_0094;
    const ccomplex_t IT_0096 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0097 = IT_0004*IT_0096;
    const ccomplex_t IT_0098 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0099 = IT_0008*IT_0098;
    const ccomplex_t IT_0100 = IT_0093 + IT_0095 + IT_0097 + IT_0099;
    const ccomplex_t IT_0101 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0102 = IT_0013*IT_0101;
    const ccomplex_t IT_0103 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0104 = IT_0016*IT_0103;
    const ccomplex_t IT_0105 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0106 = IT_0013*IT_0105;
    const ccomplex_t IT_0107 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0108 = IT_0016*IT_0107;
    const ccomplex_t IT_0109 = -IT_0102 + -IT_0104 + -IT_0106 + -IT_0108;
    const ccomplex_t IT_0110 = IT_0100 + IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0114 = IT_0004*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0116 = IT_0008*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0118 = IT_0004*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0120 = IT_0008*IT_0119;
    const ccomplex_t IT_0121 = IT_0114 + IT_0116 + IT_0118 + IT_0120;
    const ccomplex_t IT_0122 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0123 = IT_0013*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0125 = IT_0016*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0127 = IT_0013*IT_0126;
    const ccomplex_t IT_0128 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0129 = IT_0016*IT_0128;
    const ccomplex_t IT_0130 = -IT_0123 + -IT_0125 + -IT_0127 + -IT_0129;
    const ccomplex_t IT_0131 = IT_0121 + IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = IT_0112*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0031 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = IT_0134*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_13 + IT_0031 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = IT_0134*IT_0139;
    const ccomplex_t IT_0141 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0142 = IT_0004*IT_0141;
    const ccomplex_t IT_0143 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0144 = IT_0008*IT_0143;
    const ccomplex_t IT_0145 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0146 = IT_0004*IT_0145;
    const ccomplex_t IT_0147 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0148 = IT_0008*IT_0147;
    const ccomplex_t IT_0149 = IT_0142 + IT_0144 + IT_0146 + IT_0148;
    const ccomplex_t IT_0150 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0151 = IT_0013*IT_0150;
    const ccomplex_t IT_0152 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0153 = IT_0016*IT_0152;
    const ccomplex_t IT_0154 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0155 = IT_0013*IT_0154;
    const ccomplex_t IT_0156 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0157 = IT_0016*IT_0156;
    const ccomplex_t IT_0158 = -IT_0151 + -IT_0153 + -IT_0155 + -IT_0157;
    const ccomplex_t IT_0159 = IT_0149 + IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0163 = IT_0004*IT_0162;
    const ccomplex_t IT_0164 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0165 = IT_0008*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0167 = IT_0004*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0169 = IT_0008*IT_0168;
    const ccomplex_t IT_0170 = IT_0163 + IT_0165 + IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0172 = IT_0013*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0174 = IT_0016*IT_0173;
    const ccomplex_t IT_0175 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0176 = IT_0013*IT_0175;
    const ccomplex_t IT_0177 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0178 = IT_0016*IT_0177;
    const ccomplex_t IT_0179 = -IT_0172 + -IT_0174 + -IT_0176 + -IT_0178;
    const ccomplex_t IT_0180 = IT_0170 + IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0161*IT_0182;
    const ccomplex_t IT_0184 = cpow((-2)*s_23 + IT_0031 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0183*IT_0185;
    const ccomplex_t IT_0187 = cpow((-2)*s_13 + IT_0031 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = IT_0183*IT_0188;
    const ccomplex_t IT_0190 = cpow((-2)*s_23 + IT_0031 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = IT_0090*IT_0191;
    const ccomplex_t IT_0193 = (-0.5)*IT_0045 + -IT_0091 + IT_0137 + -IT_0140 
      + IT_0186 + -IT_0189 + IT_0192;
    const ccomplex_t IT_0194 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0195 = IT_0044*IT_0194;
    const ccomplex_t IT_0196 = IT_0040*IT_0194;
    const ccomplex_t IT_0197 = IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = s_14*m_N_1;
    const ccomplex_t IT_0199 = cpow((-2)*s_12 + (-2)*IT_0039 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0200 = cpow(IT_0002, -2);
    const ccomplex_t IT_0201 = IT_0011*IT_0200;
    const ccomplex_t IT_0202 = IT_0012 + IT_0201;
    const ccomplex_t IT_0203 = sin(beta);
    const ccomplex_t IT_0204 = cpow(IT_0007, 3);
    const ccomplex_t IT_0205 = cos(beta);
    const ccomplex_t IT_0206 = cpow(IT_0001, 2);
    const ccomplex_t IT_0207 = cpow(IT_0007, 2);
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*m_W*e_em*IT_0202*(IT_0203
      *IT_0204 + (-5)*IT_0001*(IT_0001*IT_0007*IT_0203 + (-0.2)*IT_0205*(IT_0206
       + (-5)*IT_0207)));
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = IT_0199*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = IT_0004*IT_0006;
    const ccomplex_t IT_0213 = IT_0000*IT_0008;
    const ccomplex_t IT_0214 = IT_0013*IT_0015;
    const ccomplex_t IT_0215 = IT_0010*IT_0016;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0212 + -IT_0213 + 
      -IT_0214 + IT_0215);
    const ccomplex_t IT_0217 = IT_0211*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0039 + -reg_prop, -1);
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*m_W*e_em*IT_0202*(IT_0204
      *IT_0205 + -IT_0001*(IT_0001*(IT_0001*IT_0203 + IT_0007*IT_0205) + 
      -IT_0203*IT_0207));
    const ccomplex_t IT_0220 = 1.5*IT_0219;
    const ccomplex_t IT_0221 = IT_0218*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = IT_0019*IT_0222;
    const ccomplex_t IT_0224 = cpow(IT_0019, 2);
    const ccomplex_t IT_0225 = m_N_1*IT_0224;
    const ccomplex_t IT_0226 = m_N_1*IT_0030;
    const ccomplex_t IT_0227 = IT_0225 + IT_0226;
    const ccomplex_t IT_0228 = IT_0033*IT_0227;
    const ccomplex_t IT_0229 = IT_0035*IT_0227;
    const ccomplex_t IT_0230 = (-0.5)*IT_0217 + (-0.5)*IT_0223 + IT_0228 + 0.5
      *IT_0229;
    const ccomplex_t IT_0231 = IT_0004*IT_0022;
    const ccomplex_t IT_0232 = IT_0008*IT_0020;
    const ccomplex_t IT_0233 = IT_0013*IT_0026;
    const ccomplex_t IT_0234 = IT_0016*IT_0024;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*(IT_0231 + -IT_0232 + 
      -IT_0233 + IT_0234);
    const ccomplex_t IT_0236 = IT_0211*IT_0235;
    const ccomplex_t IT_0237 = IT_0029*IT_0222;
    const ccomplex_t IT_0238 = cpow(IT_0029, 2);
    const ccomplex_t IT_0239 = m_N_1*IT_0238;
    const ccomplex_t IT_0240 = IT_0226 + IT_0239;
    const ccomplex_t IT_0241 = IT_0033*IT_0240;
    const ccomplex_t IT_0242 = IT_0035*IT_0240;
    const ccomplex_t IT_0243 = (-0.5)*IT_0236 + (-0.5)*IT_0237 + IT_0241 + 0.5
      *IT_0242;
    const ccomplex_t IT_0244 = IT_0230 + IT_0243;
    const ccomplex_t IT_0245 = 2*conj(IT_0038);
    const ccomplex_t IT_0246 = s_24*m_N_1;
    const ccomplex_t IT_0247 = cpow(IT_0112, 2);
    const ccomplex_t IT_0248 = m_N_3*IT_0247;
    const ccomplex_t IT_0249 = m_N_1*IT_0134;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = IT_0136*IT_0250;
    const ccomplex_t IT_0252 = IT_0139*IT_0250;
    const ccomplex_t IT_0253 = cpow(IT_0161, 2);
    const ccomplex_t IT_0254 = m_N_4*IT_0253;
    const ccomplex_t IT_0255 = m_N_1*IT_0183;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = IT_0185*IT_0256;
    const ccomplex_t IT_0258 = IT_0188*IT_0256;
    const ccomplex_t IT_0259 = cpow(IT_0068, 2);
    const ccomplex_t IT_0260 = m_N_2*IT_0259;
    const ccomplex_t IT_0261 = m_N_1*IT_0090;
    const ccomplex_t IT_0262 = IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = IT_0191*IT_0262;
    const ccomplex_t IT_0264 = IT_0047*IT_0262;
    const ccomplex_t IT_0265 = 0.5*IT_0217 + 0.5*IT_0223 + (-0.5)*IT_0229 + 
      -IT_0251 + -IT_0252 + -IT_0257 + -IT_0258 + -IT_0263 + -IT_0264;
    const ccomplex_t IT_0266 = cpow(IT_0133, 2);
    const ccomplex_t IT_0267 = m_N_3*IT_0266;
    const ccomplex_t IT_0268 = IT_0249 + IT_0267;
    const ccomplex_t IT_0269 = IT_0139*IT_0268;
    const ccomplex_t IT_0270 = cpow(IT_0182, 2);
    const ccomplex_t IT_0271 = m_N_4*IT_0270;
    const ccomplex_t IT_0272 = IT_0255 + IT_0271;
    const ccomplex_t IT_0273 = IT_0188*IT_0272;
    const ccomplex_t IT_0274 = cpow(IT_0089, 2);
    const ccomplex_t IT_0275 = m_N_2*IT_0274;
    const ccomplex_t IT_0276 = IT_0261 + IT_0275;
    const ccomplex_t IT_0277 = IT_0191*IT_0276;
    const ccomplex_t IT_0278 = IT_0047*IT_0276;
    const ccomplex_t IT_0279 = IT_0185*IT_0272;
    const ccomplex_t IT_0280 = IT_0136*IT_0268;
    const ccomplex_t IT_0281 = 0.5*IT_0236 + 0.5*IT_0237 + (-0.5)*IT_0242 + 
      -IT_0269 + -IT_0273 + -IT_0277 + -IT_0278 + -IT_0279 + -IT_0280;
    const ccomplex_t IT_0282 = 2*IT_0038;
    const ccomplex_t IT_0283 = IT_0265 + IT_0281;
    const ccomplex_t IT_0284 = (-2)*IT_0038;
    const ccomplex_t IT_0285 = IT_0194*IT_0198;
    const ccomplex_t IT_0286 = 2*s_12;
    const ccomplex_t IT_0287 = IT_0198*IT_0282;
    const ccomplex_t IT_0288 = 2*IT_0039;
    const ccomplex_t IT_0289 = 2*IT_0198;
    const ccomplex_t IT_0290 = (-2)*s_12;
    const ccomplex_t IT_0291 = (-2)*IT_0039;
    const ccomplex_t IT_0292 = (-2)*IT_0246;
    const ccomplex_t IT_0293 = 2*IT_0246;
    const ccomplex_t IT_0294 = IT_0198*IT_0284;
    const ccomplex_t IT_0295 = (-2)*IT_0198;
    const ccomplex_t IT_0296 = IT_0044*IT_0284;
    const ccomplex_t IT_0297 = IT_0040*IT_0284;
    const ccomplex_t IT_0298 = (-2)*IT_0044;
    const ccomplex_t IT_0299 = (-2)*IT_0040;
    const ccomplex_t IT_0300 = (-4)*IT_0038*conj(IT_0038)*(IT_0040 + IT_0044) 
      + 2*IT_0193*IT_0197 + IT_0198*IT_0244*IT_0245 + IT_0246*(IT_0194*IT_0244 +
       (conj(IT_0265) + conj(IT_0281))*IT_0282 + IT_0245*IT_0283 + (conj(IT_0230
      ) + conj(IT_0243))*IT_0284) + IT_0283*IT_0285 + conj(IT_0243)*(IT_0243
      *IT_0286 + IT_0287 + IT_0265*IT_0288 + IT_0193*IT_0289 + IT_0281*IT_0290 +
       IT_0230*IT_0291 + IT_0193*IT_0292) + conj(IT_0230)*(IT_0230*IT_0286 +
       IT_0287 + IT_0281*IT_0288 + IT_0193*IT_0289 + IT_0265*IT_0290 + IT_0243
      *IT_0291 + IT_0193*IT_0292) + conj(IT_0281)*(IT_0281*IT_0286 + IT_0230
      *IT_0288 + IT_0243*IT_0290 + IT_0265*IT_0291 + IT_0193*IT_0293 + IT_0294 +
       IT_0193*IT_0295) + conj(IT_0265)*(IT_0265*IT_0286 + IT_0243*IT_0288 +
       IT_0230*IT_0290 + IT_0281*IT_0291 + IT_0193*IT_0293 + IT_0294 + IT_0193
      *IT_0295) + conj(IT_0193)*(IT_0230*IT_0289 + IT_0243*IT_0292 + IT_0281
      *IT_0293 + IT_0265*IT_0295 + IT_0296 + IT_0297 + IT_0193*IT_0298 + IT_0193
      *IT_0299) + conj(IT_0193)*(IT_0243*IT_0289 + IT_0230*IT_0292 + IT_0265
      *IT_0293 + IT_0281*IT_0295 + IT_0296 + IT_0297 + IT_0193*IT_0298 + IT_0193
      *IT_0299);
    return create_ccomplex_return(IT_0300);
}

