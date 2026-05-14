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
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = IT_0030*IT_0037;
    const ccomplex_t IT_0039 = pow(m_N_1, 2);
    const ccomplex_t IT_0040 = IT_0031*IT_0039;
    const ccomplex_t IT_0041 = s_12*IT_0031;
    const ccomplex_t IT_0042 = s_14*s_24;
    const ccomplex_t IT_0043 = (-2)*IT_0042;
    const ccomplex_t IT_0044 = IT_0041 + IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + (-2)*IT_0039 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0046 = cpow(IT_0002, -2);
    const ccomplex_t IT_0047 = IT_0011*IT_0046;
    const ccomplex_t IT_0048 = IT_0012 + IT_0047;
    const ccomplex_t IT_0049 = sin(beta);
    const ccomplex_t IT_0050 = cpow(IT_0007, 3);
    const ccomplex_t IT_0051 = cos(beta);
    const ccomplex_t IT_0052 = cpow(IT_0001, 2);
    const ccomplex_t IT_0053 = cpow(IT_0007, 2);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_W*e_em*IT_0048*(IT_0049
      *IT_0050 + (-5)*IT_0001*(IT_0001*IT_0007*IT_0049 + (-0.2)*IT_0051*(IT_0052
       + (-5)*IT_0053)));
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0045*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0004*IT_0006;
    const ccomplex_t IT_0059 = IT_0000*IT_0008;
    const ccomplex_t IT_0060 = IT_0013*IT_0015;
    const ccomplex_t IT_0061 = IT_0010*IT_0016;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + -IT_0059 + 
      -IT_0060 + IT_0061);
    const ccomplex_t IT_0063 = IT_0057*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0039 + -reg_prop, -1);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_W*e_em*IT_0048*(IT_0050
      *IT_0051 + -IT_0001*(IT_0001*(IT_0001*IT_0049 + IT_0007*IT_0051) + 
      -IT_0049*IT_0053));
    const ccomplex_t IT_0066 = 1.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0064*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0019*IT_0068;
    const ccomplex_t IT_0070 = m_N_1*IT_0030;
    const ccomplex_t IT_0071 = cpow(IT_0019, 2);
    const ccomplex_t IT_0072 = m_N_1*IT_0071;
    const ccomplex_t IT_0073 = IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = IT_0033*IT_0073;
    const ccomplex_t IT_0075 = IT_0035*IT_0073;
    const ccomplex_t IT_0076 = (-0.5)*IT_0063 + -IT_0069 + IT_0074 + IT_0075;
    const ccomplex_t IT_0077 = IT_0029*IT_0068;
    const ccomplex_t IT_0078 = cpow(IT_0029, 2);
    const ccomplex_t IT_0079 = m_N_1*IT_0078;
    const ccomplex_t IT_0080 = IT_0070 + IT_0079;
    const ccomplex_t IT_0081 = IT_0035*IT_0080;
    const ccomplex_t IT_0082 = IT_0004*IT_0022;
    const ccomplex_t IT_0083 = IT_0008*IT_0020;
    const ccomplex_t IT_0084 = IT_0013*IT_0026;
    const ccomplex_t IT_0085 = IT_0016*IT_0024;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + -IT_0083 + 
      -IT_0084 + IT_0085);
    const ccomplex_t IT_0087 = IT_0057*IT_0086;
    const ccomplex_t IT_0088 = IT_0033*IT_0080;
    const ccomplex_t IT_0089 = -IT_0077 + IT_0081 + (-0.5)*IT_0087 + IT_0088;
    const ccomplex_t IT_0090 = conj(IT_0076) + conj(IT_0089);
    const ccomplex_t IT_0091 = s_14*m_N_1;
    const ccomplex_t IT_0092 = 2*IT_0038;
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = s_24*m_N_1;
    const ccomplex_t IT_0095 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0096 = IT_0004*IT_0095;
    const ccomplex_t IT_0097 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0098 = IT_0008*IT_0097;
    const ccomplex_t IT_0099 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0100 = IT_0004*IT_0099;
    const ccomplex_t IT_0101 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0102 = IT_0008*IT_0101;
    const ccomplex_t IT_0103 = IT_0096 + IT_0098 + IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0105 = IT_0013*IT_0104;
    const ccomplex_t IT_0106 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0107 = IT_0016*IT_0106;
    const ccomplex_t IT_0108 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0109 = IT_0013*IT_0108;
    const ccomplex_t IT_0110 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0111 = IT_0016*IT_0110;
    const ccomplex_t IT_0112 = -IT_0105 + -IT_0107 + -IT_0109 + -IT_0111;
    const ccomplex_t IT_0113 = IT_0103 + IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = cpow(IT_0115, 2);
    const ccomplex_t IT_0117 = m_N_2*IT_0116;
    const ccomplex_t IT_0118 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0119 = IT_0004*IT_0118;
    const ccomplex_t IT_0120 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0121 = IT_0008*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0123 = IT_0004*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0125 = IT_0008*IT_0124;
    const ccomplex_t IT_0126 = IT_0119 + IT_0121 + IT_0123 + IT_0125;
    const ccomplex_t IT_0127 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0128 = IT_0013*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0130 = IT_0016*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0132 = IT_0013*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0134 = IT_0016*IT_0133;
    const ccomplex_t IT_0135 = -IT_0128 + -IT_0130 + -IT_0132 + -IT_0134;
    const ccomplex_t IT_0136 = IT_0126 + IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = IT_0115*IT_0138;
    const ccomplex_t IT_0140 = m_N_1*IT_0139;
    const ccomplex_t IT_0141 = IT_0117 + IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_23 + IT_0031 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = IT_0141*IT_0143;
    const ccomplex_t IT_0145 = cpow((-2)*s_23 + IT_0031 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0148 = IT_0004*IT_0147;
    const ccomplex_t IT_0149 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0150 = IT_0008*IT_0149;
    const ccomplex_t IT_0151 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0152 = IT_0004*IT_0151;
    const ccomplex_t IT_0153 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0154 = IT_0008*IT_0153;
    const ccomplex_t IT_0155 = IT_0148 + IT_0150 + IT_0152 + IT_0154;
    const ccomplex_t IT_0156 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0157 = IT_0013*IT_0156;
    const ccomplex_t IT_0158 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0159 = IT_0016*IT_0158;
    const ccomplex_t IT_0160 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0161 = IT_0013*IT_0160;
    const ccomplex_t IT_0162 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0163 = IT_0016*IT_0162;
    const ccomplex_t IT_0164 = -IT_0157 + -IT_0159 + -IT_0161 + -IT_0163;
    const ccomplex_t IT_0165 = IT_0155 + IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = cpow(IT_0167, 2);
    const ccomplex_t IT_0169 = m_N_3*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0171 = IT_0004*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0173 = IT_0008*IT_0172;
    const ccomplex_t IT_0174 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0175 = IT_0004*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0177 = IT_0008*IT_0176;
    const ccomplex_t IT_0178 = IT_0171 + IT_0173 + IT_0175 + IT_0177;
    const ccomplex_t IT_0179 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0180 = IT_0013*IT_0179;
    const ccomplex_t IT_0181 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0182 = IT_0016*IT_0181;
    const ccomplex_t IT_0183 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0184 = IT_0013*IT_0183;
    const ccomplex_t IT_0185 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0186 = IT_0016*IT_0185;
    const ccomplex_t IT_0187 = -IT_0180 + -IT_0182 + -IT_0184 + -IT_0186;
    const ccomplex_t IT_0188 = IT_0178 + IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0167*IT_0190;
    const ccomplex_t IT_0192 = m_N_1*IT_0191;
    const ccomplex_t IT_0193 = IT_0169 + IT_0192;
    const ccomplex_t IT_0194 = IT_0146*IT_0193;
    const ccomplex_t IT_0195 = cpow((-2)*s_13 + IT_0031 + IT_0039 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = IT_0193*IT_0196;
    const ccomplex_t IT_0198 = cpow((-2)*s_23 + IT_0031 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0201 = IT_0004*IT_0200;
    const ccomplex_t IT_0202 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0203 = IT_0008*IT_0202;
    const ccomplex_t IT_0204 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0205 = IT_0004*IT_0204;
    const ccomplex_t IT_0206 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0207 = IT_0008*IT_0206;
    const ccomplex_t IT_0208 = IT_0201 + IT_0203 + IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0210 = IT_0013*IT_0209;
    const ccomplex_t IT_0211 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0212 = IT_0016*IT_0211;
    const ccomplex_t IT_0213 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0214 = IT_0013*IT_0213;
    const ccomplex_t IT_0215 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0216 = IT_0016*IT_0215;
    const ccomplex_t IT_0217 = -IT_0210 + -IT_0212 + -IT_0214 + -IT_0216;
    const ccomplex_t IT_0218 = IT_0208 + IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = cpow(IT_0220, 2);
    const ccomplex_t IT_0222 = m_N_4*IT_0221;
    const ccomplex_t IT_0223 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0224 = IT_0004*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0226 = IT_0008*IT_0225;
    const ccomplex_t IT_0227 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0228 = IT_0004*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0230 = IT_0008*IT_0229;
    const ccomplex_t IT_0231 = IT_0224 + IT_0226 + IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0233 = IT_0013*IT_0232;
    const ccomplex_t IT_0234 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0235 = IT_0016*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0237 = IT_0013*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0239 = IT_0016*IT_0238;
    const ccomplex_t IT_0240 = -IT_0233 + -IT_0235 + -IT_0237 + -IT_0239;
    const ccomplex_t IT_0241 = IT_0231 + IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = IT_0220*IT_0243;
    const ccomplex_t IT_0245 = m_N_1*IT_0244;
    const ccomplex_t IT_0246 = IT_0222 + IT_0245;
    const ccomplex_t IT_0247 = IT_0199*IT_0246;
    const ccomplex_t IT_0248 = cpow((-2)*s_13 + IT_0031 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = IT_0246*IT_0249;
    const ccomplex_t IT_0251 = cpow((-2)*s_13 + IT_0031 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = IT_0141*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0063 + -IT_0144 + -IT_0194 + -IT_0197 + 
      -IT_0247 + -IT_0250 + -IT_0253;
    const ccomplex_t IT_0255 = cpow(IT_0138, 2);
    const ccomplex_t IT_0256 = m_N_2*IT_0255;
    const ccomplex_t IT_0257 = IT_0140 + IT_0256;
    const ccomplex_t IT_0258 = IT_0252*IT_0257;
    const ccomplex_t IT_0259 = IT_0143*IT_0257;
    const ccomplex_t IT_0260 = cpow(IT_0190, 2);
    const ccomplex_t IT_0261 = m_N_3*IT_0260;
    const ccomplex_t IT_0262 = IT_0192 + IT_0261;
    const ccomplex_t IT_0263 = IT_0146*IT_0262;
    const ccomplex_t IT_0264 = IT_0196*IT_0262;
    const ccomplex_t IT_0265 = cpow(IT_0243, 2);
    const ccomplex_t IT_0266 = m_N_4*IT_0265;
    const ccomplex_t IT_0267 = IT_0245 + IT_0266;
    const ccomplex_t IT_0268 = IT_0199*IT_0267;
    const ccomplex_t IT_0269 = IT_0249*IT_0267;
    const ccomplex_t IT_0270 = 0.5*IT_0087 + -IT_0258 + -IT_0259 + -IT_0263 + 
      -IT_0264 + -IT_0268 + -IT_0269;
    const ccomplex_t IT_0271 = conj(IT_0254) + conj(IT_0270);
    const ccomplex_t IT_0272 = (-2)*IT_0038;
    const ccomplex_t IT_0273 = IT_0091*IT_0272;
    const ccomplex_t IT_0274 = IT_0146*IT_0191;
    const ccomplex_t IT_0275 = IT_0191*IT_0196;
    const ccomplex_t IT_0276 = IT_0199*IT_0244;
    const ccomplex_t IT_0277 = IT_0244*IT_0249;
    const ccomplex_t IT_0278 = IT_0139*IT_0252;
    const ccomplex_t IT_0279 = IT_0139*IT_0143;
    const ccomplex_t IT_0280 = IT_0274 + -IT_0275 + IT_0276 + -IT_0277 + 
      -IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = IT_0044*IT_0272*conj(IT_0280);
    const ccomplex_t IT_0282 = IT_0040*IT_0272*conj(IT_0280);
    const ccomplex_t IT_0283 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0284 = IT_0044*IT_0280*IT_0283;
    const ccomplex_t IT_0285 = IT_0040*IT_0280*IT_0283;
    const ccomplex_t IT_0286 = (-2)*IT_0044*IT_0280*conj(IT_0280);
    const ccomplex_t IT_0287 = 2*IT_0094;
    const ccomplex_t IT_0288 = (-2)*IT_0091;
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = 2*IT_0091;
    const ccomplex_t IT_0291 = (-2)*IT_0094;
    const ccomplex_t IT_0292 = IT_0290 + IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0040;
    const ccomplex_t IT_0294 = 2*s_12;
    const ccomplex_t IT_0295 = 2*IT_0039;
    const ccomplex_t IT_0296 = (-2)*s_12;
    const ccomplex_t IT_0297 = (-2)*IT_0039;
    const ccomplex_t IT_0298 = IT_0094*IT_0283;
    const ccomplex_t IT_0299 = 2*conj(IT_0038);
    const ccomplex_t IT_0300 = IT_0091*IT_0299;
    const ccomplex_t IT_0301 = conj(IT_0280)*IT_0292;
    const ccomplex_t IT_0302 = IT_0091*IT_0283;
    const ccomplex_t IT_0303 = conj(IT_0280)*IT_0289;
    const ccomplex_t IT_0304 = (-4)*IT_0038*conj(IT_0038)*(IT_0040 + IT_0044) 
      + IT_0090*IT_0093 + IT_0094*(IT_0092*IT_0271 + IT_0090*IT_0272) + IT_0271
      *IT_0273 + 2*IT_0281 + 2*IT_0282 + 2*IT_0284 + 2*IT_0285 + 2*IT_0286 +
       IT_0280*(IT_0271*IT_0289 + IT_0090*IT_0292 + 2*conj(IT_0280)*IT_0293) +
       IT_0089*(conj(IT_0089)*IT_0294 + conj(IT_0254)*IT_0295 + conj(IT_0270)
      *IT_0296 + conj(IT_0076)*IT_0297 + IT_0298 + IT_0300 + IT_0301) + IT_0076*
      (conj(IT_0076)*IT_0294 + conj(IT_0270)*IT_0295 + conj(IT_0254)*IT_0296 +
       conj(IT_0089)*IT_0297 + IT_0298 + IT_0300 + IT_0301) + IT_0270*(conj
      (IT_0270)*IT_0294 + conj(IT_0076)*IT_0295 + conj(IT_0089)*IT_0296 + conj
      (IT_0254)*IT_0297 + IT_0094*IT_0299 + IT_0302 + IT_0303) + IT_0254*(conj
      (IT_0254)*IT_0294 + conj(IT_0089)*IT_0295 + conj(IT_0076)*IT_0296 + conj
      (IT_0270)*IT_0297 + IT_0094*IT_0299 + IT_0302 + IT_0303);
    return create_ccomplex_return(IT_0304);
}

