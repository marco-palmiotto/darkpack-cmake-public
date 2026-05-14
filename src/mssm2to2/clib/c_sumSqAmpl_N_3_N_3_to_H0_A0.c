#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_3_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_3_to_H0_A0(
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
    const creal_t s_23 = param->s_23;
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
    const ccomplex_t IT_0000 = pow(m_N_3, 2);
    const ccomplex_t IT_0001 = pow(m_H0, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0001 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = IT_0006*IT_0010;
    const ccomplex_t IT_0012 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0004*IT_0015;
    const ccomplex_t IT_0017 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0010*IT_0015;
    const ccomplex_t IT_0020 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0009 + -IT_0013 + -IT_0018 + IT_0021;
    const ccomplex_t IT_0023 = cos(alpha);
    const ccomplex_t IT_0024 = IT_0006*IT_0023;
    const ccomplex_t IT_0025 = IT_0012*IT_0024;
    const ccomplex_t IT_0026 = sin(alpha);
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = IT_0008*IT_0027;
    const ccomplex_t IT_0029 = IT_0015*IT_0023;
    const ccomplex_t IT_0030 = IT_0020*IT_0029;
    const ccomplex_t IT_0031 = IT_0015*IT_0026;
    const ccomplex_t IT_0032 = IT_0017*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0025 + -IT_0028 + 
      -IT_0030 + IT_0032);
    const ccomplex_t IT_0034 = IT_0022*IT_0033;
    const ccomplex_t IT_0035 = m_N_3*IT_0034;
    const ccomplex_t IT_0036 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0037 = IT_0024*IT_0036;
    const ccomplex_t IT_0038 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0039 = IT_0027*IT_0038;
    const ccomplex_t IT_0040 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0041 = IT_0029*IT_0040;
    const ccomplex_t IT_0042 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0043 = IT_0031*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0037 + -IT_0039 + 
      -IT_0041 + IT_0043);
    const ccomplex_t IT_0045 = IT_0022*IT_0044;
    const ccomplex_t IT_0046 = m_N_3*IT_0045;
    const ccomplex_t IT_0047 = IT_0035 + IT_0046;
    const ccomplex_t IT_0048 = IT_0003*IT_0047;
    const ccomplex_t IT_0049 = IT_0007*IT_0012;
    const ccomplex_t IT_0050 = IT_0008*IT_0011;
    const ccomplex_t IT_0051 = IT_0016*IT_0020;
    const ccomplex_t IT_0052 = IT_0017*IT_0019;
    const ccomplex_t IT_0053 = IT_0049 + IT_0050 + -IT_0051 + -IT_0052;
    const ccomplex_t IT_0054 = cpow(IT_0005, -2);
    const ccomplex_t IT_0055 = IT_0014*IT_0054;
    const ccomplex_t IT_0056 = IT_0015 + IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*IT_0010*
      (IT_0004*IT_0023 + -IT_0010*IT_0026)*IT_0056;
    const ccomplex_t IT_0058 = pow(m_Z, 2);
    const ccomplex_t IT_0059 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0058 + 
      -reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0066 = IT_0007*IT_0065;
    const ccomplex_t IT_0067 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0068 = IT_0011*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0070 = IT_0016*IT_0069;
    const ccomplex_t IT_0071 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0072 = IT_0019*IT_0071;
    const ccomplex_t IT_0073 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0074 = IT_0007*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0076 = IT_0011*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0078 = IT_0016*IT_0077;
    const ccomplex_t IT_0079 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0080 = IT_0019*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0066 + (-0.5)*IT_0068 + (-0.5)*IT_0070 +
       0.5*IT_0072 + 0.5*IT_0074 + (-0.5)*IT_0076 + (-0.5)*IT_0078 + 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0024*IT_0067;
    const ccomplex_t IT_0083 = IT_0031*IT_0069;
    const ccomplex_t IT_0084 = IT_0024*IT_0075;
    const ccomplex_t IT_0085 = IT_0031*IT_0077;
    const ccomplex_t IT_0086 = IT_0082 + IT_0083 + IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = IT_0027*IT_0065;
    const ccomplex_t IT_0088 = IT_0029*IT_0071;
    const ccomplex_t IT_0089 = IT_0027*IT_0073;
    const ccomplex_t IT_0090 = IT_0029*IT_0079;
    const ccomplex_t IT_0091 = -IT_0087 + -IT_0088 + -IT_0089 + -IT_0090;
    const ccomplex_t IT_0092 = IT_0086 + IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0081*IT_0094;
    const ccomplex_t IT_0096 = m_N_2*IT_0095;
    const ccomplex_t IT_0097 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0098 = IT_0024*IT_0097;
    const ccomplex_t IT_0099 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0100 = IT_0031*IT_0099;
    const ccomplex_t IT_0101 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0102 = IT_0024*IT_0101;
    const ccomplex_t IT_0103 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0104 = IT_0031*IT_0103;
    const ccomplex_t IT_0105 = IT_0098 + IT_0100 + IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0107 = IT_0027*IT_0106;
    const ccomplex_t IT_0108 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0109 = IT_0029*IT_0108;
    const ccomplex_t IT_0110 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0111 = IT_0027*IT_0110;
    const ccomplex_t IT_0112 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0113 = IT_0029*IT_0112;
    const ccomplex_t IT_0114 = -IT_0107 + -IT_0109 + -IT_0111 + -IT_0113;
    const ccomplex_t IT_0115 = IT_0105 + IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0081*IT_0117;
    const ccomplex_t IT_0119 = m_N_3*IT_0118;
    const ccomplex_t IT_0120 = IT_0096 + IT_0119;
    const ccomplex_t IT_0121 = IT_0064*IT_0120;
    const ccomplex_t IT_0122 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = e_em*IT_0006*IT_0015*(IT_0010*IT_0023 + 
      -IT_0004*IT_0026);
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0006*IT_0014;
    const ccomplex_t IT_0127 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = IT_0005*IT_0015;
    const ccomplex_t IT_0130 = IT_0127*IT_0129;
    const ccomplex_t IT_0131 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0132 = IT_0126*IT_0131;
    const ccomplex_t IT_0133 = IT_0129*IT_0131;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0128 + IT_0130 + 
      -IT_0132 + -IT_0133);
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = IT_0125*IT_0135;
    const ccomplex_t IT_0137 = m_N_3*IT_0136;
    const ccomplex_t IT_0138 = 0.5*IT_0134;
    const ccomplex_t IT_0139 = IT_0125*IT_0138;
    const ccomplex_t IT_0140 = m_N_3*IT_0139;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0137 + IT_0141;
    const ccomplex_t IT_0143 = IT_0123*IT_0142;
    const ccomplex_t IT_0144 = cpow(IT_0010, 3);
    const ccomplex_t IT_0145 = -IT_0004;
    const ccomplex_t IT_0146 = cpow(IT_0010, 2);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*m_W*e_em*IT_0056*(IT_0026
      *IT_0144 + IT_0145*((IT_0004*IT_0023 + -IT_0010*IT_0026)*IT_0145 + IT_0023
      *IT_0146));
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = IT_0022*IT_0151;
    const ccomplex_t IT_0153 = cpow((-2)*s_23 + IT_0001 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0047*IT_0154;
    const ccomplex_t IT_0156 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0159 = IT_0007*IT_0158;
    const ccomplex_t IT_0160 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0161 = IT_0011*IT_0160;
    const ccomplex_t IT_0162 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0163 = IT_0016*IT_0162;
    const ccomplex_t IT_0164 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0165 = IT_0019*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0167 = IT_0007*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0169 = IT_0011*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0171 = IT_0016*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0173 = IT_0019*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0159 + (-0.5)*IT_0161 + (-0.5)*IT_0163 +
       0.5*IT_0165 + 0.5*IT_0167 + (-0.5)*IT_0169 + (-0.5)*IT_0171 + 0.5*IT_0173;
    const ccomplex_t IT_0175 = IT_0024*IT_0160;
    const ccomplex_t IT_0176 = IT_0031*IT_0162;
    const ccomplex_t IT_0177 = IT_0024*IT_0168;
    const ccomplex_t IT_0178 = IT_0031*IT_0170;
    const ccomplex_t IT_0179 = IT_0175 + IT_0176 + IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = IT_0027*IT_0158;
    const ccomplex_t IT_0181 = IT_0029*IT_0164;
    const ccomplex_t IT_0182 = IT_0027*IT_0166;
    const ccomplex_t IT_0183 = IT_0029*IT_0172;
    const ccomplex_t IT_0184 = -IT_0180 + -IT_0181 + -IT_0182 + -IT_0183;
    const ccomplex_t IT_0185 = IT_0179 + IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = IT_0174*IT_0187;
    const ccomplex_t IT_0189 = m_N_1*IT_0188;
    const ccomplex_t IT_0190 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0191 = IT_0024*IT_0190;
    const ccomplex_t IT_0192 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0193 = IT_0031*IT_0192;
    const ccomplex_t IT_0194 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0195 = IT_0024*IT_0194;
    const ccomplex_t IT_0196 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0197 = IT_0031*IT_0196;
    const ccomplex_t IT_0198 = IT_0191 + IT_0193 + IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0200 = IT_0027*IT_0199;
    const ccomplex_t IT_0201 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0202 = IT_0029*IT_0201;
    const ccomplex_t IT_0203 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0204 = IT_0027*IT_0203;
    const ccomplex_t IT_0205 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0206 = IT_0029*IT_0205;
    const ccomplex_t IT_0207 = -IT_0200 + -IT_0202 + -IT_0204 + -IT_0206;
    const ccomplex_t IT_0208 = IT_0198 + IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = 0.5*IT_0209;
    const ccomplex_t IT_0211 = IT_0174*IT_0210;
    const ccomplex_t IT_0212 = m_N_3*IT_0211;
    const ccomplex_t IT_0213 = IT_0189 + IT_0212;
    const ccomplex_t IT_0214 = IT_0157*IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = IT_0213*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = IT_0120*IT_0219;
    const ccomplex_t IT_0221 = 0.5*IT_0048 + (-0.5)*IT_0062 + IT_0121 +
       IT_0143 + (-0.5)*IT_0152 + IT_0155 + IT_0214 + IT_0217 + IT_0220;
    const ccomplex_t IT_0222 = IT_0007*IT_0036;
    const ccomplex_t IT_0223 = IT_0011*IT_0038;
    const ccomplex_t IT_0224 = IT_0016*IT_0040;
    const ccomplex_t IT_0225 = IT_0019*IT_0042;
    const ccomplex_t IT_0226 = -IT_0222 + -IT_0223 + IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = IT_0061*IT_0226;
    const ccomplex_t IT_0228 = IT_0007*IT_0038;
    const ccomplex_t IT_0229 = IT_0011*IT_0036;
    const ccomplex_t IT_0230 = IT_0016*IT_0042;
    const ccomplex_t IT_0231 = IT_0019*IT_0040;
    const ccomplex_t IT_0232 = -IT_0228 + IT_0229 + IT_0230 + -IT_0231;
    const ccomplex_t IT_0233 = IT_0151*IT_0232;
    const ccomplex_t IT_0234 = IT_0044*IT_0232;
    const ccomplex_t IT_0235 = m_N_3*IT_0234;
    const ccomplex_t IT_0236 = IT_0033*IT_0232;
    const ccomplex_t IT_0237 = m_N_3*IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = IT_0003*IT_0238;
    const ccomplex_t IT_0240 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0243 = IT_0007*IT_0242;
    const ccomplex_t IT_0244 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0245 = IT_0011*IT_0244;
    const ccomplex_t IT_0246 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0247 = IT_0016*IT_0246;
    const ccomplex_t IT_0248 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0249 = IT_0019*IT_0248;
    const ccomplex_t IT_0250 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0251 = IT_0007*IT_0250;
    const ccomplex_t IT_0252 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0253 = IT_0011*IT_0252;
    const ccomplex_t IT_0254 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0255 = IT_0016*IT_0254;
    const ccomplex_t IT_0256 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0257 = IT_0019*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0243 + 0.5*IT_0245 + 0.5*IT_0247 + (
      -0.5)*IT_0249 + (-0.5)*IT_0251 + 0.5*IT_0253 + 0.5*IT_0255 + (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0024*IT_0244;
    const ccomplex_t IT_0260 = IT_0031*IT_0246;
    const ccomplex_t IT_0261 = IT_0024*IT_0252;
    const ccomplex_t IT_0262 = IT_0031*IT_0254;
    const ccomplex_t IT_0263 = IT_0259 + IT_0260 + IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = IT_0027*IT_0242;
    const ccomplex_t IT_0265 = IT_0029*IT_0248;
    const ccomplex_t IT_0266 = IT_0027*IT_0250;
    const ccomplex_t IT_0267 = IT_0029*IT_0256;
    const ccomplex_t IT_0268 = -IT_0264 + -IT_0265 + -IT_0266 + -IT_0267;
    const ccomplex_t IT_0269 = IT_0263 + IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = 0.5*IT_0270;
    const ccomplex_t IT_0272 = IT_0258*IT_0271;
    const ccomplex_t IT_0273 = m_N_4*IT_0272;
    const ccomplex_t IT_0274 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0275 = IT_0024*IT_0274;
    const ccomplex_t IT_0276 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0277 = IT_0031*IT_0276;
    const ccomplex_t IT_0278 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0279 = IT_0024*IT_0278;
    const ccomplex_t IT_0280 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0281 = IT_0031*IT_0280;
    const ccomplex_t IT_0282 = IT_0275 + IT_0277 + IT_0279 + IT_0281;
    const ccomplex_t IT_0283 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0284 = IT_0027*IT_0283;
    const ccomplex_t IT_0285 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0286 = IT_0029*IT_0285;
    const ccomplex_t IT_0287 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0288 = IT_0027*IT_0287;
    const ccomplex_t IT_0289 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0290 = IT_0029*IT_0289;
    const ccomplex_t IT_0291 = -IT_0284 + -IT_0286 + -IT_0288 + -IT_0290;
    const ccomplex_t IT_0292 = IT_0282 + IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*IT_0292;
    const ccomplex_t IT_0294 = 0.5*IT_0293;
    const ccomplex_t IT_0295 = IT_0258*IT_0294;
    const ccomplex_t IT_0296 = m_N_3*IT_0295;
    const ccomplex_t IT_0297 = IT_0273 + IT_0296;
    const ccomplex_t IT_0298 = IT_0241*IT_0297;
    const ccomplex_t IT_0299 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0297*IT_0300;
    const ccomplex_t IT_0302 = 0.5*IT_0227 + 0.5*IT_0233 + (-0.5)*IT_0239 + 
      -IT_0298 + -IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0302;
    const ccomplex_t IT_0304 = IT_0007*IT_0199;
    const ccomplex_t IT_0305 = IT_0011*IT_0190;
    const ccomplex_t IT_0306 = IT_0016*IT_0192;
    const ccomplex_t IT_0307 = IT_0019*IT_0201;
    const ccomplex_t IT_0308 = IT_0007*IT_0203;
    const ccomplex_t IT_0309 = IT_0011*IT_0194;
    const ccomplex_t IT_0310 = IT_0016*IT_0196;
    const ccomplex_t IT_0311 = IT_0019*IT_0205;
    const ccomplex_t IT_0312 = (-0.5)*IT_0304 + 0.5*IT_0305 + 0.5*IT_0306 + (
      -0.5)*IT_0307 + (-0.5)*IT_0308 + 0.5*IT_0309 + 0.5*IT_0310 + (-0.5)*IT_0311;
    const ccomplex_t IT_0313 = IT_0210*IT_0312;
    const ccomplex_t IT_0314 = m_N_1*IT_0313;
    const ccomplex_t IT_0315 = IT_0187*IT_0312;
    const ccomplex_t IT_0316 = m_N_3*IT_0315;
    const ccomplex_t IT_0317 = IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = IT_0216*IT_0317;
    const ccomplex_t IT_0319 = IT_0157*IT_0317;
    const ccomplex_t IT_0320 = IT_0007*IT_0106;
    const ccomplex_t IT_0321 = IT_0011*IT_0097;
    const ccomplex_t IT_0322 = IT_0016*IT_0099;
    const ccomplex_t IT_0323 = IT_0019*IT_0108;
    const ccomplex_t IT_0324 = IT_0007*IT_0110;
    const ccomplex_t IT_0325 = IT_0011*IT_0101;
    const ccomplex_t IT_0326 = IT_0016*IT_0103;
    const ccomplex_t IT_0327 = IT_0019*IT_0112;
    const ccomplex_t IT_0328 = (-0.5)*IT_0320 + 0.5*IT_0321 + 0.5*IT_0322 + (
      -0.5)*IT_0323 + (-0.5)*IT_0324 + 0.5*IT_0325 + 0.5*IT_0326 + (-0.5)*IT_0327;
    const ccomplex_t IT_0329 = IT_0117*IT_0328;
    const ccomplex_t IT_0330 = m_N_2*IT_0329;
    const ccomplex_t IT_0331 = IT_0094*IT_0328;
    const ccomplex_t IT_0332 = m_N_3*IT_0331;
    const ccomplex_t IT_0333 = IT_0330 + IT_0332;
    const ccomplex_t IT_0334 = IT_0219*IT_0333;
    const ccomplex_t IT_0335 = IT_0064*IT_0333;
    const ccomplex_t IT_0336 = IT_0154*IT_0238;
    const ccomplex_t IT_0337 = -IT_0143 + (-0.5)*IT_0227 + (-0.5)*IT_0233 +
       0.5*IT_0239 + IT_0318 + IT_0319 + IT_0334 + IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = (-2)*IT_0337;
    const ccomplex_t IT_0339 = IT_0303 + IT_0338;
    const ccomplex_t IT_0340 = 2*conj(IT_0302);
    const ccomplex_t IT_0341 = (-2)*conj(IT_0337);
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = IT_0007*IT_0283;
    const ccomplex_t IT_0344 = IT_0011*IT_0274;
    const ccomplex_t IT_0345 = IT_0016*IT_0276;
    const ccomplex_t IT_0346 = IT_0019*IT_0285;
    const ccomplex_t IT_0347 = IT_0007*IT_0287;
    const ccomplex_t IT_0348 = IT_0011*IT_0278;
    const ccomplex_t IT_0349 = IT_0016*IT_0280;
    const ccomplex_t IT_0350 = IT_0019*IT_0289;
    const ccomplex_t IT_0351 = 0.5*IT_0343 + (-0.5)*IT_0344 + (-0.5)*IT_0345 +
       0.5*IT_0346 + 0.5*IT_0347 + (-0.5)*IT_0348 + (-0.5)*IT_0349 + 0.5*IT_0350;
    const ccomplex_t IT_0352 = IT_0271*IT_0351;
    const ccomplex_t IT_0353 = m_N_3*IT_0352;
    const ccomplex_t IT_0354 = IT_0294*IT_0351;
    const ccomplex_t IT_0355 = m_N_4*IT_0354;
    const ccomplex_t IT_0356 = IT_0353 + IT_0355;
    const ccomplex_t IT_0357 = IT_0241*IT_0356;
    const ccomplex_t IT_0358 = IT_0300*IT_0356;
    const ccomplex_t IT_0359 = (-0.5)*IT_0048 + 0.5*IT_0062 + 0.5*IT_0152 + 
      -IT_0357 + -IT_0358;
    const ccomplex_t IT_0360 = 2*IT_0337;
    const ccomplex_t IT_0361 = (-2)*IT_0302;
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = 2*conj(IT_0337);
    const ccomplex_t IT_0364 = (-2)*conj(IT_0302);
    const ccomplex_t IT_0365 = IT_0363 + IT_0364;
    const ccomplex_t IT_0366 = 2*s_12;
    const ccomplex_t IT_0367 = (-2)*s_12;
    const ccomplex_t IT_0368 = s_13*m_N_3;
    const ccomplex_t IT_0369 = IT_0003*IT_0045;
    const ccomplex_t IT_0370 = IT_0241*IT_0295;
    const ccomplex_t IT_0371 = IT_0300*IT_0352;
    const ccomplex_t IT_0372 = 0.5*IT_0369 + -IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = IT_0154*IT_0236;
    const ccomplex_t IT_0374 = IT_0219*IT_0331;
    const ccomplex_t IT_0375 = IT_0064*IT_0118;
    const ccomplex_t IT_0376 = IT_0124*IT_0138;
    const ccomplex_t IT_0377 = IT_0123*IT_0376;
    const ccomplex_t IT_0378 = IT_0157*IT_0315;
    const ccomplex_t IT_0379 = IT_0211*IT_0216;
    const ccomplex_t IT_0380 = 0.5*IT_0369 + -IT_0373 + -IT_0374 + IT_0375 + 
      -IT_0377 + -IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = conj(IT_0372) + conj(IT_0380);
    const ccomplex_t IT_0382 = IT_0372 + IT_0380;
    const ccomplex_t IT_0383 = IT_0003*IT_0236;
    const ccomplex_t IT_0384 = IT_0241*IT_0352;
    const ccomplex_t IT_0385 = IT_0295*IT_0300;
    const ccomplex_t IT_0386 = 0.5*IT_0383 + -IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = IT_0045*IT_0154;
    const ccomplex_t IT_0388 = IT_0118*IT_0219;
    const ccomplex_t IT_0389 = IT_0064*IT_0331;
    const ccomplex_t IT_0390 = IT_0124*IT_0135;
    const ccomplex_t IT_0391 = IT_0123*IT_0390;
    const ccomplex_t IT_0392 = IT_0157*IT_0211;
    const ccomplex_t IT_0393 = IT_0216*IT_0315;
    const ccomplex_t IT_0394 = 0.5*IT_0383 + -IT_0387 + -IT_0388 + IT_0389 + 
      -IT_0391 + -IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0386) + conj(IT_0394);
    const ccomplex_t IT_0396 = 2*IT_0221;
    const ccomplex_t IT_0397 = (-2)*IT_0359;
    const ccomplex_t IT_0398 = IT_0386 + IT_0394;
    const ccomplex_t IT_0399 = 2*conj(IT_0221);
    const ccomplex_t IT_0400 = (-2)*conj(IT_0359);
    const ccomplex_t IT_0401 = s_23*m_N_3;
    const ccomplex_t IT_0402 = 2*IT_0359;
    const ccomplex_t IT_0403 = (-2)*IT_0221;
    const ccomplex_t IT_0404 = 2*conj(IT_0359);
    const ccomplex_t IT_0405 = (-2)*conj(IT_0221);
    const ccomplex_t IT_0406 = s_12*IT_0001;
    const ccomplex_t IT_0407 = s_13*s_23;
    const ccomplex_t IT_0408 = (-2)*IT_0407;
    const ccomplex_t IT_0409 = IT_0406 + IT_0408;
    const ccomplex_t IT_0410 = (-2)*IT_0372;
    const ccomplex_t IT_0411 = (-2)*conj(IT_0372);
    const ccomplex_t IT_0412 = (-2)*IT_0386;
    const ccomplex_t IT_0413 = (-2)*conj(IT_0386);
    const ccomplex_t IT_0414 = (-2)*IT_0394;
    const ccomplex_t IT_0415 = IT_0000*IT_0001;
    const ccomplex_t IT_0416 = IT_0000*(conj(IT_0221)*IT_0339 + IT_0221
      *IT_0342 + conj(IT_0359)*IT_0362 + IT_0359*IT_0365) + (IT_0221*conj
      (IT_0221) + IT_0302*conj(IT_0302) + IT_0337*conj(IT_0337) + IT_0359*conj
      (IT_0359))*IT_0366 + (conj(IT_0302)*IT_0337 + IT_0302*conj(IT_0337) + conj
      (IT_0221)*IT_0359 + IT_0221*conj(IT_0359))*IT_0367 + IT_0368*(IT_0362
      *IT_0381 + IT_0365*IT_0382 + IT_0395*(IT_0396 + IT_0397) + IT_0398*
      (IT_0399 + IT_0400)) + IT_0401*(IT_0339*IT_0395 + IT_0342*IT_0398 +
       IT_0381*(IT_0402 + IT_0403) + IT_0382*(IT_0404 + IT_0405)) + IT_0409*
      (IT_0381*IT_0410 + (-2)*IT_0380*(conj(IT_0380) + (-0.5)*IT_0411) + IT_0395
      *IT_0412 + IT_0394*IT_0413 + conj(IT_0394)*IT_0414) + (conj(IT_0386)
      *IT_0410 + (-2)*conj(IT_0394)*(IT_0380 + (-0.5)*IT_0410) + IT_0398*IT_0411
       + IT_0380*IT_0413 + conj(IT_0380)*(IT_0412 + IT_0414))*IT_0415;
    return create_ccomplex_return(IT_0416);
}

