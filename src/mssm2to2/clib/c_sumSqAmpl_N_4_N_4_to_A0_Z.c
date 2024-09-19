#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_A0_Z(
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
    const ccomplex_t IT_0011 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = -IT_0012 + IT_0016 + IT_0021 + -IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0017;
    const ccomplex_t IT_0027 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = IT_0008*IT_0018;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0029*IT_0031;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0028 + IT_0030 + 
      -IT_0032 + -IT_0033);
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = pow(m_A0, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0037 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0036*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0034;
    const ccomplex_t IT_0043 = IT_0025*IT_0042;
    const ccomplex_t IT_0044 = IT_0039*IT_0043;
    const ccomplex_t IT_0045 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0046 = IT_0010*IT_0045;
    const ccomplex_t IT_0047 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0048 = IT_0014*IT_0047;
    const ccomplex_t IT_0049 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0050 = IT_0019*IT_0049;
    const ccomplex_t IT_0051 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0052 = IT_0022*IT_0051;
    const ccomplex_t IT_0053 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0054 = IT_0010*IT_0053;
    const ccomplex_t IT_0055 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0056 = IT_0014*IT_0055;
    const ccomplex_t IT_0057 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0058 = IT_0019*IT_0057;
    const ccomplex_t IT_0059 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0060 = IT_0022*IT_0059;
    const ccomplex_t IT_0061 = (-0.5)*IT_0046 + 0.5*IT_0048 + 0.5*IT_0050 + (
      -0.5)*IT_0052 + (-0.5)*IT_0054 + 0.5*IT_0056 + 0.5*IT_0058 + (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0063 = IT_0026*IT_0062;
    const ccomplex_t IT_0064 = IT_0029*IT_0062;
    const ccomplex_t IT_0065 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0066 = IT_0026*IT_0065;
    const ccomplex_t IT_0067 = IT_0029*IT_0065;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0063 + IT_0064 + 
      -IT_0066 + -IT_0067);
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0061*IT_0069;
    const ccomplex_t IT_0071 = pow(m_N_4, 2);
    const ccomplex_t IT_0072 = cpow((-2)*s_13 + IT_0037 + IT_0071 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0070*IT_0073;
    const ccomplex_t IT_0075 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0076 = IT_0010*IT_0075;
    const ccomplex_t IT_0077 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0078 = IT_0014*IT_0077;
    const ccomplex_t IT_0079 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0080 = IT_0019*IT_0079;
    const ccomplex_t IT_0081 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0082 = IT_0022*IT_0081;
    const ccomplex_t IT_0083 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0084 = IT_0010*IT_0083;
    const ccomplex_t IT_0085 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0086 = IT_0014*IT_0085;
    const ccomplex_t IT_0087 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0088 = IT_0019*IT_0087;
    const ccomplex_t IT_0089 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0090 = IT_0022*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0076 + 0.5*IT_0078 + 0.5*IT_0080 + (
      -0.5)*IT_0082 + (-0.5)*IT_0084 + 0.5*IT_0086 + 0.5*IT_0088 + (-0.5)*IT_0090;
    const ccomplex_t IT_0092 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0093 = IT_0026*IT_0092;
    const ccomplex_t IT_0094 = IT_0029*IT_0092;
    const ccomplex_t IT_0095 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0096 = IT_0026*IT_0095;
    const ccomplex_t IT_0097 = IT_0029*IT_0095;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0093 + IT_0094 + 
      -IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0091*IT_0099;
    const ccomplex_t IT_0101 = cpow((-2)*s_13 + IT_0037 + IT_0071 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100*IT_0102;
    const ccomplex_t IT_0104 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0105 = IT_0010*IT_0104;
    const ccomplex_t IT_0106 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0107 = IT_0014*IT_0106;
    const ccomplex_t IT_0108 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0109 = IT_0019*IT_0108;
    const ccomplex_t IT_0110 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0111 = IT_0022*IT_0110;
    const ccomplex_t IT_0112 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0113 = IT_0010*IT_0112;
    const ccomplex_t IT_0114 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0115 = IT_0014*IT_0114;
    const ccomplex_t IT_0116 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0117 = IT_0019*IT_0116;
    const ccomplex_t IT_0118 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0119 = IT_0022*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0105 + 0.5*IT_0107 + 0.5*IT_0109 + (
      -0.5)*IT_0111 + (-0.5)*IT_0113 + 0.5*IT_0115 + 0.5*IT_0117 + (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0122 = IT_0026*IT_0121;
    const ccomplex_t IT_0123 = IT_0029*IT_0121;
    const ccomplex_t IT_0124 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0125 = IT_0026*IT_0124;
    const ccomplex_t IT_0126 = IT_0029*IT_0124;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0122 + IT_0123 + 
      -IT_0125 + -IT_0126);
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = IT_0120*IT_0128;
    const ccomplex_t IT_0130 = cpow((-2)*s_13 + IT_0037 + IT_0071 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0129*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0044 + -IT_0074 + -IT_0103 + -IT_0132;
    const ccomplex_t IT_0134 = conj(IT_0041) + -conj(IT_0133);
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0037 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = IT_0036*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_23 + IT_0037 + IT_0071 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = IT_0129*IT_0139;
    const ccomplex_t IT_0141 = cpow((-2)*s_23 + IT_0037 + IT_0071 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = IT_0070*IT_0142;
    const ccomplex_t IT_0144 = (-0.5)*IT_0137 + IT_0140 + IT_0143;
    const ccomplex_t IT_0145 = IT_0043*IT_0136;
    const ccomplex_t IT_0146 = cpow((-2)*s_23 + IT_0037 + IT_0071 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = IT_0100*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0145 + IT_0148;
    const ccomplex_t IT_0150 = IT_0144 + IT_0149;
    const ccomplex_t IT_0151 = IT_0041 + -IT_0133;
    const ccomplex_t IT_0152 = conj(IT_0144) + conj(IT_0149);
    const ccomplex_t IT_0153 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0154 = IT_0010*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0156 = IT_0014*IT_0155;
    const ccomplex_t IT_0157 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0158 = IT_0019*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0160 = IT_0022*IT_0159;
    const ccomplex_t IT_0161 = IT_0154 + -IT_0156 + -IT_0158 + IT_0160;
    const ccomplex_t IT_0162 = IT_0042*IT_0161;
    const ccomplex_t IT_0163 = IT_0039*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0166 = IT_0010*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0168 = IT_0014*IT_0167;
    const ccomplex_t IT_0169 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0170 = IT_0019*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0172 = IT_0022*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0174 = IT_0010*IT_0173;
    const ccomplex_t IT_0175 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0176 = IT_0014*IT_0175;
    const ccomplex_t IT_0177 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0178 = IT_0019*IT_0177;
    const ccomplex_t IT_0179 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0180 = IT_0022*IT_0179;
    const ccomplex_t IT_0181 = 0.5*IT_0166 + (-0.5)*IT_0168 + (-0.5)*IT_0170 +
       0.5*IT_0172 + 0.5*IT_0174 + (-0.5)*IT_0176 + (-0.5)*IT_0178 + 0.5*IT_0180;
    const ccomplex_t IT_0182 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0183 = IT_0026*IT_0182;
    const ccomplex_t IT_0184 = IT_0029*IT_0182;
    const ccomplex_t IT_0185 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0186 = IT_0026*IT_0185;
    const ccomplex_t IT_0187 = IT_0029*IT_0185;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0183 + IT_0184 + 
      -IT_0186 + -IT_0187);
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = IT_0181*IT_0189;
    const ccomplex_t IT_0191 = IT_0073*IT_0190;
    const ccomplex_t IT_0192 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0193 = IT_0010*IT_0192;
    const ccomplex_t IT_0194 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0195 = IT_0014*IT_0194;
    const ccomplex_t IT_0196 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0197 = IT_0019*IT_0196;
    const ccomplex_t IT_0198 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0199 = IT_0022*IT_0198;
    const ccomplex_t IT_0200 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0201 = IT_0010*IT_0200;
    const ccomplex_t IT_0202 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0203 = IT_0014*IT_0202;
    const ccomplex_t IT_0204 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0205 = IT_0019*IT_0204;
    const ccomplex_t IT_0206 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0207 = IT_0022*IT_0206;
    const ccomplex_t IT_0208 = 0.5*IT_0193 + (-0.5)*IT_0195 + (-0.5)*IT_0197 +
       0.5*IT_0199 + 0.5*IT_0201 + (-0.5)*IT_0203 + (-0.5)*IT_0205 + 0.5*IT_0207;
    const ccomplex_t IT_0209 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0210 = IT_0026*IT_0209;
    const ccomplex_t IT_0211 = IT_0029*IT_0209;
    const ccomplex_t IT_0212 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0213 = IT_0026*IT_0212;
    const ccomplex_t IT_0214 = IT_0029*IT_0212;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0210 + IT_0211 + 
      -IT_0213 + -IT_0214);
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0208*IT_0216;
    const ccomplex_t IT_0218 = IT_0131*IT_0217;
    const ccomplex_t IT_0219 = IT_0035*IT_0161;
    const ccomplex_t IT_0220 = IT_0039*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0222 = IT_0010*IT_0221;
    const ccomplex_t IT_0223 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0224 = IT_0014*IT_0223;
    const ccomplex_t IT_0225 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0226 = IT_0019*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0228 = IT_0022*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0230 = IT_0010*IT_0229;
    const ccomplex_t IT_0231 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0232 = IT_0014*IT_0231;
    const ccomplex_t IT_0233 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0234 = IT_0019*IT_0233;
    const ccomplex_t IT_0235 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0236 = IT_0022*IT_0235;
    const ccomplex_t IT_0237 = 0.5*IT_0222 + (-0.5)*IT_0224 + (-0.5)*IT_0226 +
       0.5*IT_0228 + 0.5*IT_0230 + (-0.5)*IT_0232 + (-0.5)*IT_0234 + 0.5*IT_0236;
    const ccomplex_t IT_0238 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0239 = IT_0026*IT_0238;
    const ccomplex_t IT_0240 = IT_0029*IT_0238;
    const ccomplex_t IT_0241 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0242 = IT_0026*IT_0241;
    const ccomplex_t IT_0243 = IT_0029*IT_0241;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0239 + IT_0240 + 
      -IT_0242 + -IT_0243);
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = IT_0237*IT_0245;
    const ccomplex_t IT_0247 = IT_0102*IT_0246;
    const ccomplex_t IT_0248 = -IT_0191 + -IT_0218 + (-0.5)*IT_0220 + -IT_0247;
    const ccomplex_t IT_0249 = IT_0136*IT_0162;
    const ccomplex_t IT_0250 = IT_0139*IT_0217;
    const ccomplex_t IT_0251 = IT_0142*IT_0190;
    const ccomplex_t IT_0252 = (-0.5)*IT_0249 + IT_0250 + IT_0251;
    const ccomplex_t IT_0253 = IT_0136*IT_0219;
    const ccomplex_t IT_0254 = IT_0147*IT_0246;
    const ccomplex_t IT_0255 = 0.5*IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = IT_0252 + IT_0255;
    const ccomplex_t IT_0257 = IT_0164 + -IT_0248;
    const ccomplex_t IT_0258 = conj(IT_0252) + conj(IT_0255);
    const ccomplex_t IT_0259 = pow(m_Z, 2);
    const ccomplex_t IT_0260 = s_14*IT_0259;
    const ccomplex_t IT_0261 = IT_0001*IT_0260;
    const ccomplex_t IT_0262 = 0.5*IT_0261;
    const ccomplex_t IT_0263 = s_14 + IT_0262;
    const ccomplex_t IT_0264 = m_N_4*IT_0263;
    const ccomplex_t IT_0265 = 2*IT_0264;
    const ccomplex_t IT_0266 = m_N_4*IT_0036;
    const ccomplex_t IT_0267 = IT_0136*IT_0266;
    const ccomplex_t IT_0268 = m_N_4*IT_0162;
    const ccomplex_t IT_0269 = m_N_4*IT_0219;
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = IT_0268 + IT_0270;
    const ccomplex_t IT_0272 = IT_0039*IT_0271;
    const ccomplex_t IT_0273 = IT_0061*IT_0189;
    const ccomplex_t IT_0274 = m_N_3*IT_0273;
    const ccomplex_t IT_0275 = m_N_4*IT_0070;
    const ccomplex_t IT_0276 = -IT_0275;
    const ccomplex_t IT_0277 = IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = IT_0142*IT_0277;
    const ccomplex_t IT_0279 = m_N_4*IT_0043;
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = IT_0266 + IT_0280;
    const ccomplex_t IT_0282 = IT_0136*IT_0281;
    const ccomplex_t IT_0283 = IT_0120*IT_0216;
    const ccomplex_t IT_0284 = m_N_1*IT_0283;
    const ccomplex_t IT_0285 = m_N_4*IT_0129;
    const ccomplex_t IT_0286 = -IT_0285;
    const ccomplex_t IT_0287 = IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = IT_0139*IT_0287;
    const ccomplex_t IT_0289 = m_N_4*IT_0100;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = IT_0091*IT_0245;
    const ccomplex_t IT_0292 = m_N_2*IT_0291;
    const ccomplex_t IT_0293 = IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = IT_0147*IT_0293;
    const ccomplex_t IT_0295 = 0.5*IT_0267 + (-0.5)*IT_0272 + IT_0278 + 0.5
      *IT_0282 + IT_0288 + IT_0294;
    const ccomplex_t IT_0296 = 2*IT_0295;
    const ccomplex_t IT_0297 = IT_0099*IT_0237;
    const ccomplex_t IT_0298 = m_N_2*IT_0297;
    const ccomplex_t IT_0299 = m_N_4*IT_0246;
    const ccomplex_t IT_0300 = -IT_0299;
    const ccomplex_t IT_0301 = IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = IT_0102*IT_0301;
    const ccomplex_t IT_0303 = IT_0136*IT_0279;
    const ccomplex_t IT_0304 = IT_0128*IT_0208;
    const ccomplex_t IT_0305 = m_N_1*IT_0304;
    const ccomplex_t IT_0306 = m_N_4*IT_0217;
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = IT_0305 + IT_0307;
    const ccomplex_t IT_0309 = IT_0131*IT_0308;
    const ccomplex_t IT_0310 = IT_0069*IT_0181;
    const ccomplex_t IT_0311 = m_N_3*IT_0310;
    const ccomplex_t IT_0312 = m_N_4*IT_0190;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = IT_0311 + IT_0313;
    const ccomplex_t IT_0315 = IT_0073*IT_0314;
    const ccomplex_t IT_0316 = (-0.5)*IT_0272 + -IT_0302 + (-0.5)*IT_0303 + 
      -IT_0309 + -IT_0315;
    const ccomplex_t IT_0317 = 2*IT_0316;
    const ccomplex_t IT_0318 = IT_0296 + IT_0317;
    const ccomplex_t IT_0319 = 2*conj(IT_0295);
    const ccomplex_t IT_0320 = 2*conj(IT_0316);
    const ccomplex_t IT_0321 = IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = IT_0136*IT_0268;
    const ccomplex_t IT_0323 = IT_0039*IT_0281;
    const ccomplex_t IT_0324 = IT_0142*IT_0314;
    const ccomplex_t IT_0325 = IT_0136*IT_0271;
    const ccomplex_t IT_0326 = IT_0139*IT_0308;
    const ccomplex_t IT_0327 = IT_0147*IT_0301;
    const ccomplex_t IT_0328 = 0.5*IT_0322 + (-0.5)*IT_0323 + IT_0324 + 0.5
      *IT_0325 + IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = 2*IT_0328;
    const ccomplex_t IT_0330 = IT_0136*IT_0269;
    const ccomplex_t IT_0331 = IT_0073*IT_0277;
    const ccomplex_t IT_0332 = IT_0131*IT_0287;
    const ccomplex_t IT_0333 = IT_0102*IT_0293;
    const ccomplex_t IT_0334 = (-0.5)*IT_0323 + (-0.5)*IT_0330 + -IT_0331 + 
      -IT_0332 + -IT_0333;
    const ccomplex_t IT_0335 = 2*IT_0334;
    const ccomplex_t IT_0336 = IT_0329 + IT_0335;
    const ccomplex_t IT_0337 = 2*conj(IT_0328);
    const ccomplex_t IT_0338 = 2*conj(IT_0334);
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = pow(s_14, 2);
    const ccomplex_t IT_0341 = IT_0001*IT_0340;
    const ccomplex_t IT_0342 = -IT_0071;
    const ccomplex_t IT_0343 = IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = IT_0071*IT_0343;
    const ccomplex_t IT_0345 = cos(alpha);
    const ccomplex_t IT_0346 = IT_0009*IT_0345;
    const ccomplex_t IT_0347 = IT_0153*IT_0346;
    const ccomplex_t IT_0348 = sin(alpha);
    const ccomplex_t IT_0349 = IT_0009*IT_0348;
    const ccomplex_t IT_0350 = IT_0155*IT_0349;
    const ccomplex_t IT_0351 = IT_0018*IT_0345;
    const ccomplex_t IT_0352 = IT_0157*IT_0351;
    const ccomplex_t IT_0353 = IT_0018*IT_0348;
    const ccomplex_t IT_0354 = IT_0159*IT_0353;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*(IT_0347 + IT_0350 + 
      -IT_0352 + -IT_0354);
    const ccomplex_t IT_0356 = -IT_0355;
    const ccomplex_t IT_0357 = e_em*IT_0009*IT_0018*(IT_0007*IT_0345 + IT_0013
      *IT_0348);
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = IT_0356*IT_0358;
    const ccomplex_t IT_0360 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0071 + -reg_prop, -1);
    const ccomplex_t IT_0361 = (0 + _Complex_I*1)*IT_0360;
    const ccomplex_t IT_0362 = IT_0359*IT_0361;
    const ccomplex_t IT_0363 = IT_0155*IT_0346;
    const ccomplex_t IT_0364 = IT_0153*IT_0349;
    const ccomplex_t IT_0365 = IT_0159*IT_0351;
    const ccomplex_t IT_0366 = IT_0157*IT_0353;
    const ccomplex_t IT_0367 = (0 + _Complex_I*1)*(IT_0363 + -IT_0364 + 
      -IT_0365 + IT_0366);
    const ccomplex_t IT_0368 = e_em*IT_0009*IT_0018*(IT_0013*IT_0345 + 
      -IT_0007*IT_0348);
    const ccomplex_t IT_0369 = IT_0367*IT_0368;
    const ccomplex_t IT_0370 = cpow((-2)*s_12 + (-2)*IT_0071 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0371 = (0 + _Complex_I*1)*IT_0370;
    const ccomplex_t IT_0372 = IT_0369*IT_0371;
    const ccomplex_t IT_0373 = 0.5*IT_0362 + 0.5*IT_0372;
    const ccomplex_t IT_0374 = IT_0011*IT_0346;
    const ccomplex_t IT_0375 = IT_0015*IT_0349;
    const ccomplex_t IT_0376 = IT_0020*IT_0351;
    const ccomplex_t IT_0377 = IT_0023*IT_0353;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*(IT_0374 + IT_0375 + 
      -IT_0376 + -IT_0377);
    const ccomplex_t IT_0379 = -IT_0378;
    const ccomplex_t IT_0380 = IT_0358*IT_0379;
    const ccomplex_t IT_0381 = IT_0361*IT_0380;
    const ccomplex_t IT_0382 = IT_0015*IT_0346;
    const ccomplex_t IT_0383 = IT_0011*IT_0349;
    const ccomplex_t IT_0384 = IT_0023*IT_0351;
    const ccomplex_t IT_0385 = IT_0020*IT_0353;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*(IT_0382 + -IT_0383 + 
      -IT_0384 + IT_0385);
    const ccomplex_t IT_0387 = IT_0368*IT_0386;
    const ccomplex_t IT_0388 = IT_0371*IT_0387;
    const ccomplex_t IT_0389 = IT_0137 + (-2)*IT_0140 + (-2)*IT_0143 + 
      -IT_0145 + (-2)*IT_0148 + (-0.5)*IT_0381 + (-0.5)*IT_0388;
    const ccomplex_t IT_0390 = 2*IT_0389;
    const ccomplex_t IT_0391 = 0.5*IT_0381 + 0.5*IT_0388;
    const ccomplex_t IT_0392 = (-2)*IT_0391;
    const ccomplex_t IT_0393 = IT_0390 + IT_0392;
    const ccomplex_t IT_0394 = 2*conj(IT_0389);
    const ccomplex_t IT_0395 = (-2)*conj(IT_0391);
    const ccomplex_t IT_0396 = IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = IT_0249 + (-2)*IT_0250 + (-2)*IT_0251 + 
      -IT_0253 + (-2)*IT_0254 + (-0.5)*IT_0362 + (-0.5)*IT_0372;
    const ccomplex_t IT_0398 = 2*IT_0391;
    const ccomplex_t IT_0399 = (-2)*IT_0389;
    const ccomplex_t IT_0400 = IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = 2*conj(IT_0391);
    const ccomplex_t IT_0402 = (-2)*conj(IT_0389);
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = s_12*s_14;
    const ccomplex_t IT_0405 = s_24*IT_0071;
    const ccomplex_t IT_0406 = -IT_0405;
    const ccomplex_t IT_0407 = IT_0404 + IT_0406;
    const ccomplex_t IT_0408 = conj(IT_0133) + IT_0152;
    const ccomplex_t IT_0409 = IT_0133 + IT_0150;
    const ccomplex_t IT_0410 = conj(IT_0248) + IT_0258;
    const ccomplex_t IT_0411 = 2*IT_0397;
    const ccomplex_t IT_0412 = (-2)*IT_0373;
    const ccomplex_t IT_0413 = IT_0411 + IT_0412;
    const ccomplex_t IT_0414 = IT_0248 + IT_0256;
    const ccomplex_t IT_0415 = 2*conj(IT_0397);
    const ccomplex_t IT_0416 = (-2)*conj(IT_0373);
    const ccomplex_t IT_0417 = IT_0415 + IT_0416;
    const ccomplex_t IT_0418 = 2*IT_0373;
    const ccomplex_t IT_0419 = (-2)*IT_0397;
    const ccomplex_t IT_0420 = IT_0418 + IT_0419;
    const ccomplex_t IT_0421 = 2*conj(IT_0373);
    const ccomplex_t IT_0422 = (-2)*conj(IT_0397);
    const ccomplex_t IT_0423 = IT_0421 + IT_0422;
    const ccomplex_t IT_0424 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0425 = -IT_0424;
    const ccomplex_t IT_0426 = s_12 + IT_0425;
    const ccomplex_t IT_0427 = s_12*IT_0426;
    const ccomplex_t IT_0428 = IT_0044 + 2*IT_0074 + 2*IT_0103 + 2*IT_0132 +
       0.5*IT_0381 + 0.5*IT_0388;
    const ccomplex_t IT_0429 = IT_0040 + (-0.5)*IT_0381 + (-0.5)*IT_0388;
    const ccomplex_t IT_0430 = 2*IT_0191 + 2*IT_0218 + IT_0220 + 2*IT_0247 +
       0.5*IT_0362 + 0.5*IT_0372;
    const ccomplex_t IT_0431 = IT_0163 + (-0.5)*IT_0362 + (-0.5)*IT_0372;
    const ccomplex_t IT_0432 = IT_0071*IT_0426;
    const ccomplex_t IT_0433 = 2*IT_0429;
    const ccomplex_t IT_0434 = (-2)*IT_0428;
    const ccomplex_t IT_0435 = IT_0433 + IT_0434;
    const ccomplex_t IT_0436 = 2*conj(IT_0429);
    const ccomplex_t IT_0437 = (-2)*conj(IT_0428);
    const ccomplex_t IT_0438 = IT_0436 + IT_0437;
    const ccomplex_t IT_0439 = 2*IT_0428;
    const ccomplex_t IT_0440 = (-2)*IT_0429;
    const ccomplex_t IT_0441 = IT_0439 + IT_0440;
    const ccomplex_t IT_0442 = 2*conj(IT_0428);
    const ccomplex_t IT_0443 = (-2)*conj(IT_0429);
    const ccomplex_t IT_0444 = IT_0442 + IT_0443;
    const ccomplex_t IT_0445 = pow(s_24, 2);
    const ccomplex_t IT_0446 = IT_0001*IT_0445;
    const ccomplex_t IT_0447 = -IT_0446;
    const ccomplex_t IT_0448 = IT_0071 + IT_0447;
    const ccomplex_t IT_0449 = IT_0071*IT_0448;
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = m_N_4*IT_0448;
    const ccomplex_t IT_0452 = -IT_0451;
    const ccomplex_t IT_0453 = conj(IT_0295) + conj(IT_0316);
    const ccomplex_t IT_0454 = IT_0295 + IT_0316;
    const ccomplex_t IT_0455 = conj(IT_0328) + conj(IT_0334);
    const ccomplex_t IT_0456 = 2*IT_0430;
    const ccomplex_t IT_0457 = (-2)*IT_0431;
    const ccomplex_t IT_0458 = IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = IT_0328 + IT_0334;
    const ccomplex_t IT_0460 = 2*conj(IT_0430);
    const ccomplex_t IT_0461 = (-2)*conj(IT_0431);
    const ccomplex_t IT_0462 = IT_0460 + IT_0461;
    const ccomplex_t IT_0463 = s_24*IT_0259;
    const ccomplex_t IT_0464 = IT_0001*IT_0463;
    const ccomplex_t IT_0465 = 0.5*IT_0464;
    const ccomplex_t IT_0466 = s_24 + IT_0465;
    const ccomplex_t IT_0467 = m_N_4*IT_0466;
    const ccomplex_t IT_0468 = 2*IT_0467;
    const ccomplex_t IT_0469 = (-2)*IT_0295;
    const ccomplex_t IT_0470 = (-2)*IT_0316;
    const ccomplex_t IT_0471 = IT_0469 + IT_0470;
    const ccomplex_t IT_0472 = (-2)*conj(IT_0295);
    const ccomplex_t IT_0473 = (-2)*conj(IT_0316);
    const ccomplex_t IT_0474 = IT_0472 + IT_0473;
    const ccomplex_t IT_0475 = (-2)*IT_0328;
    const ccomplex_t IT_0476 = (-2)*IT_0334;
    const ccomplex_t IT_0477 = IT_0475 + IT_0476;
    const ccomplex_t IT_0478 = (-2)*conj(IT_0328);
    const ccomplex_t IT_0479 = (-2)*conj(IT_0334);
    const ccomplex_t IT_0480 = IT_0478 + IT_0479;
    const ccomplex_t IT_0481 = (-0.25)*IT_0261;
    const ccomplex_t IT_0482 = s_14 + IT_0481;
    const ccomplex_t IT_0483 = m_N_4*IT_0482;
    const ccomplex_t IT_0484 = (-4)*IT_0483;
    const ccomplex_t IT_0485 = (-0.25)*IT_0464;
    const ccomplex_t IT_0486 = s_24 + IT_0485;
    const ccomplex_t IT_0487 = m_N_4*IT_0486;
    const ccomplex_t IT_0488 = (-4)*IT_0487;
    const ccomplex_t IT_0489 = IT_0071*IT_0259;
    const ccomplex_t IT_0490 = 3*IT_0489;
    const ccomplex_t IT_0491 = (-2)*IT_0041;
    const ccomplex_t IT_0492 = 0.5*IT_0491;
    const ccomplex_t IT_0493 = IT_0133 + IT_0492;
    const ccomplex_t IT_0494 = (-2)*conj(IT_0041);
    const ccomplex_t IT_0495 = 0.5*IT_0494;
    const ccomplex_t IT_0496 = conj(IT_0133) + IT_0495;
    const ccomplex_t IT_0497 = (-2)*IT_0144;
    const ccomplex_t IT_0498 = (-2)*IT_0149;
    const ccomplex_t IT_0499 = IT_0497 + IT_0498;
    const ccomplex_t IT_0500 = (-2)*conj(IT_0144);
    const ccomplex_t IT_0501 = (-2)*conj(IT_0149);
    const ccomplex_t IT_0502 = IT_0500 + IT_0501;
    const ccomplex_t IT_0503 = (-3)*IT_0489;
    const ccomplex_t IT_0504 = s_12*s_24;
    const ccomplex_t IT_0505 = s_14*IT_0071;
    const ccomplex_t IT_0506 = -IT_0505;
    const ccomplex_t IT_0507 = IT_0504 + IT_0506;
    const ccomplex_t IT_0508 = 2*IT_0431;
    const ccomplex_t IT_0509 = -IT_0164*(conj(IT_0430) + -conj(IT_0431)) + 0.5
      *conj(IT_0041)*IT_0435 + 0.5*IT_0041*IT_0438 + 0.5*IT_0408*IT_0441 + 0.5
      *IT_0409*IT_0444 + 0.5*IT_0410*IT_0458 + 0.5*IT_0414*IT_0462 + -conj
      (IT_0164)*(IT_0430 + (-0.5)*IT_0508);
    const ccomplex_t IT_0510 = s_12*IT_0259;
    const ccomplex_t IT_0511 = (-2)*IT_0000;
    const ccomplex_t IT_0512 = IT_0510 + IT_0511;
    const ccomplex_t IT_0513 = IT_0001*IT_0512;
    const ccomplex_t IT_0514 = (-0.5)*IT_0513;
    const ccomplex_t IT_0515 = s_12 + IT_0514;
    const ccomplex_t IT_0516 = 4*IT_0328;
    const ccomplex_t IT_0517 = 6*conj(IT_0295);
    const ccomplex_t IT_0518 = (-0.5)*IT_0004;
    const ccomplex_t IT_0519 = IT_0510 + IT_0511 + IT_0518;
    const ccomplex_t IT_0520 = (-4)*IT_0144;
    const ccomplex_t IT_0521 = (-4)*IT_0252;
    const ccomplex_t IT_0522 = m_N_4*IT_0343;
    const ccomplex_t IT_0523 = s_12*IT_0343;
    const ccomplex_t IT_0524 = m_N_4*IT_0426;
    const ccomplex_t IT_0525 = conj(IT_0295)*IT_0393 + IT_0295*IT_0396 +
       IT_0319*IT_0430 + IT_0296*conj(IT_0430) + (IT_0411 + IT_0412 + IT_0439 +
       IT_0440)*IT_0455 + conj(IT_0316)*(IT_0390 + IT_0392 + IT_0456 + IT_0457) 
      + (IT_0415 + IT_0416 + IT_0442 + IT_0443)*IT_0459 + IT_0316*(IT_0394 +
       IT_0395 + IT_0460 + IT_0461) + conj(IT_0431)*IT_0469 + IT_0431*IT_0472;
    const ccomplex_t IT_0526 = s_12*IT_0448;
    const ccomplex_t IT_0527 = -IT_0526;
    const ccomplex_t IT_0528 = (-8)*IT_0006*(IT_0134*IT_0150 + IT_0151*IT_0152
       + (conj(IT_0164) + -conj(IT_0248))*IT_0256 + IT_0257*IT_0258) + IT_0265*
      (IT_0258*IT_0318 + IT_0256*IT_0321 + IT_0152*IT_0336 + IT_0150*IT_0339) +
       IT_0344*(conj(IT_0373)*IT_0393 + IT_0373*IT_0396 + conj(IT_0397)*IT_0400 
      + IT_0397*IT_0403) + IT_0407*(conj(IT_0041)*IT_0400 + IT_0041*IT_0403 +
       IT_0393*IT_0408 + IT_0396*IT_0409 + IT_0410*IT_0413 + IT_0414*IT_0417 +
       conj(IT_0164)*IT_0420 + IT_0164*IT_0423) + IT_0427*(IT_0396*IT_0428 +
       IT_0393*conj(IT_0428) + IT_0403*IT_0429 + IT_0400*conj(IT_0429) + IT_0417
      *IT_0430 + IT_0413*conj(IT_0430) + IT_0423*IT_0431 + IT_0420*conj(IT_0431)
      ) + IT_0432*(IT_0403*IT_0430 + IT_0400*conj(IT_0430) + IT_0396*IT_0431 +
       IT_0393*conj(IT_0431) + conj(IT_0397)*IT_0435 + IT_0397*IT_0438 + conj
      (IT_0373)*IT_0441 + IT_0373*IT_0444) + (conj(IT_0430)*IT_0435 + IT_0430
      *IT_0438 + conj(IT_0431)*IT_0441 + IT_0431*IT_0444)*IT_0450 + IT_0452*
      (IT_0441*IT_0453 + IT_0444*IT_0454 + IT_0455*IT_0458 + IT_0459*IT_0462) +
       IT_0468*(conj(IT_0133)*IT_0318 + IT_0133*IT_0321 + conj(IT_0248)*IT_0336 
      + IT_0248*IT_0339 + conj(IT_0041)*IT_0471 + IT_0041*IT_0474 + conj(IT_0164
      )*IT_0477 + IT_0164*IT_0480) + (conj(IT_0164)*IT_0318 + IT_0164*IT_0321 +
       conj(IT_0041)*IT_0336 + IT_0041*IT_0339 + conj(IT_0248)*IT_0471 + IT_0248
      *IT_0474 + conj(IT_0133)*IT_0477 + IT_0133*IT_0480)*IT_0484 + (IT_0152
      *IT_0471 + IT_0150*IT_0474 + IT_0258*IT_0477 + IT_0256*IT_0480)*IT_0488 +
       2*IT_0490*(IT_0152*IT_0248 + IT_0150*conj(IT_0248) + IT_0258*IT_0493 +
       IT_0256*IT_0496 + 0.5*conj(IT_0164)*IT_0499 + 0.5*IT_0164*IT_0502) + 2*
      (IT_0134*IT_0248 + IT_0151*conj(IT_0248) + conj(IT_0164)*IT_0493 + IT_0164
      *IT_0496 + 0.5*IT_0258*IT_0499 + 0.5*IT_0256*IT_0502)*IT_0503 + 2*IT_0507
      *IT_0509 + 4*IT_0515*(IT_0453*IT_0454 + IT_0455*(IT_0334 + 0.25*IT_0516)) 
      + 6*IT_0071*(IT_0454*IT_0455 + IT_0459*(conj(IT_0316) + 0.166666666666667
      *IT_0517)) + (-4)*IT_0519*(IT_0134*IT_0151 + (conj(IT_0164) + -conj
      (IT_0248))*IT_0257 + IT_0152*(IT_0149 + (-0.25)*IT_0520) + IT_0258*
      (IT_0255 + (-0.25)*IT_0521)) + (IT_0319*IT_0397 + IT_0296*conj(IT_0397) +
       conj(IT_0316)*IT_0413 + IT_0316*IT_0417 + IT_0393*IT_0455 + IT_0396
      *IT_0459 + conj(IT_0373)*IT_0469 + IT_0373*IT_0472)*IT_0522 + (conj
      (IT_0389)*IT_0390 + conj(IT_0391)*IT_0400 + IT_0391*IT_0402 + conj(IT_0397
      )*IT_0411 + conj(IT_0373)*IT_0420 + IT_0373*IT_0422)*IT_0523 + IT_0524
      *IT_0525 + (conj(IT_0429)*IT_0433 + conj(IT_0428)*IT_0441 + IT_0428
      *IT_0443 + conj(IT_0430)*IT_0458 + IT_0430*IT_0461 + conj(IT_0431)*IT_0508
      )*IT_0527;
    return create_ccomplex_return(IT_0528);
}

