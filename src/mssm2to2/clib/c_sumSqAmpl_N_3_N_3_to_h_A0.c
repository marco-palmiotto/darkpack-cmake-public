#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_3_to_h_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_3_to_h_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
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
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0012;
    const ccomplex_t IT_0017 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0006 + IT_0010 + -IT_0015 + -IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0002, -2);
    const ccomplex_t IT_0021 = IT_0011*IT_0020;
    const ccomplex_t IT_0022 = IT_0012 + IT_0021;
    const ccomplex_t IT_0023 = cos(alpha);
    const ccomplex_t IT_0024 = sin(alpha);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*IT_0007
      *IT_0022*(IT_0007*IT_0023 + IT_0001*IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = pow(m_Z, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0027 + 
      -reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0019*IT_0030;
    const ccomplex_t IT_0032 = pow(m_h, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0032 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0004*IT_0009;
    const ccomplex_t IT_0036 = IT_0005*IT_0008;
    const ccomplex_t IT_0037 = IT_0013*IT_0017;
    const ccomplex_t IT_0038 = IT_0014*IT_0016;
    const ccomplex_t IT_0039 = IT_0035 + -IT_0036 + -IT_0037 + IT_0038;
    const ccomplex_t IT_0040 = IT_0003*IT_0023;
    const ccomplex_t IT_0041 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = IT_0003*IT_0024;
    const ccomplex_t IT_0044 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = IT_0012*IT_0023;
    const ccomplex_t IT_0047 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = IT_0012*IT_0024;
    const ccomplex_t IT_0050 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0042 + IT_0045 + 
      -IT_0048 + -IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0039*IT_0053;
    const ccomplex_t IT_0055 = m_N_3*IT_0054;
    const ccomplex_t IT_0056 = IT_0009*IT_0040;
    const ccomplex_t IT_0057 = IT_0005*IT_0043;
    const ccomplex_t IT_0058 = IT_0017*IT_0046;
    const ccomplex_t IT_0059 = IT_0014*IT_0049;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0056 + IT_0057 + 
      -IT_0058 + -IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0039*IT_0061;
    const ccomplex_t IT_0063 = m_N_3*IT_0062;
    const ccomplex_t IT_0064 = IT_0055 + IT_0063;
    const ccomplex_t IT_0065 = IT_0034*IT_0064;
    const ccomplex_t IT_0066 = cpow(IT_0007, 3);
    const ccomplex_t IT_0067 = cpow(IT_0007, 2);
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0023
      *IT_0066 + -IT_0001*(IT_0001*(IT_0007*IT_0023 + IT_0001*IT_0024) + 
      -IT_0024*IT_0067));
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0039*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_13 + IT_0000 + IT_0032 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0077 = IT_0004*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0079 = IT_0008*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0081 = IT_0013*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0083 = IT_0016*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0085 = IT_0004*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0089 = IT_0013*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0091 = IT_0016*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0077 + (-0.5)*IT_0079 + (-0.5)*IT_0081 +
       0.5*IT_0083 + 0.5*IT_0085 + (-0.5)*IT_0087 + (-0.5)*IT_0089 + 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0040*IT_0076;
    const ccomplex_t IT_0094 = IT_0043*IT_0078;
    const ccomplex_t IT_0095 = IT_0040*IT_0084;
    const ccomplex_t IT_0096 = IT_0043*IT_0086;
    const ccomplex_t IT_0097 = IT_0093 + IT_0094 + IT_0095 + IT_0096;
    const ccomplex_t IT_0098 = IT_0046*IT_0080;
    const ccomplex_t IT_0099 = IT_0049*IT_0082;
    const ccomplex_t IT_0100 = IT_0046*IT_0088;
    const ccomplex_t IT_0101 = IT_0049*IT_0090;
    const ccomplex_t IT_0102 = -IT_0098 + -IT_0099 + -IT_0100 + -IT_0101;
    const ccomplex_t IT_0103 = IT_0097 + IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = IT_0092*IT_0105;
    const ccomplex_t IT_0107 = m_N_4*IT_0106;
    const ccomplex_t IT_0108 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0109 = IT_0040*IT_0108;
    const ccomplex_t IT_0110 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0111 = IT_0043*IT_0110;
    const ccomplex_t IT_0112 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0113 = IT_0040*IT_0112;
    const ccomplex_t IT_0114 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0115 = IT_0043*IT_0114;
    const ccomplex_t IT_0116 = IT_0109 + IT_0111 + IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0118 = IT_0046*IT_0117;
    const ccomplex_t IT_0119 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0120 = IT_0049*IT_0119;
    const ccomplex_t IT_0121 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0122 = IT_0046*IT_0121;
    const ccomplex_t IT_0123 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0124 = IT_0049*IT_0123;
    const ccomplex_t IT_0125 = -IT_0118 + -IT_0120 + -IT_0122 + -IT_0124;
    const ccomplex_t IT_0126 = IT_0116 + IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = (-0.5)*IT_0127;
    const ccomplex_t IT_0129 = IT_0092*IT_0128;
    const ccomplex_t IT_0130 = m_N_3*IT_0129;
    const ccomplex_t IT_0131 = IT_0107 + IT_0130;
    const ccomplex_t IT_0132 = IT_0075*IT_0131;
    const ccomplex_t IT_0133 = cpow((-2)*s_23 + IT_0000 + IT_0032 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0131*IT_0134;
    const ccomplex_t IT_0136 = 0.5*IT_0031 + (-0.5)*IT_0065 + 0.5*IT_0073 + 
      -IT_0132 + -IT_0135;
    const ccomplex_t IT_0137 = IT_0004*IT_0041;
    const ccomplex_t IT_0138 = IT_0008*IT_0044;
    const ccomplex_t IT_0139 = IT_0013*IT_0047;
    const ccomplex_t IT_0140 = IT_0016*IT_0050;
    const ccomplex_t IT_0141 = -IT_0137 + IT_0138 + IT_0139 + -IT_0140;
    const ccomplex_t IT_0142 = IT_0072*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_13 + IT_0032 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0053*IT_0141;
    const ccomplex_t IT_0146 = m_N_3*IT_0145;
    const ccomplex_t IT_0147 = IT_0061*IT_0141;
    const ccomplex_t IT_0148 = m_N_3*IT_0147;
    const ccomplex_t IT_0149 = IT_0146 + IT_0148;
    const ccomplex_t IT_0150 = IT_0144*IT_0149;
    const ccomplex_t IT_0151 = IT_0034*IT_0149;
    const ccomplex_t IT_0152 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0153 = IT_0004*IT_0152;
    const ccomplex_t IT_0154 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0155 = IT_0008*IT_0154;
    const ccomplex_t IT_0156 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0157 = IT_0013*IT_0156;
    const ccomplex_t IT_0158 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0159 = IT_0016*IT_0158;
    const ccomplex_t IT_0160 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0161 = IT_0004*IT_0160;
    const ccomplex_t IT_0162 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0163 = IT_0008*IT_0162;
    const ccomplex_t IT_0164 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0165 = IT_0013*IT_0164;
    const ccomplex_t IT_0166 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0167 = IT_0016*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0153 + 0.5*IT_0155 + 0.5*IT_0157 + (
      -0.5)*IT_0159 + (-0.5)*IT_0161 + 0.5*IT_0163 + 0.5*IT_0165 + (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0040*IT_0152;
    const ccomplex_t IT_0170 = IT_0043*IT_0154;
    const ccomplex_t IT_0171 = IT_0040*IT_0160;
    const ccomplex_t IT_0172 = IT_0043*IT_0162;
    const ccomplex_t IT_0173 = IT_0169 + IT_0170 + IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = IT_0046*IT_0156;
    const ccomplex_t IT_0175 = IT_0049*IT_0158;
    const ccomplex_t IT_0176 = IT_0046*IT_0164;
    const ccomplex_t IT_0177 = IT_0049*IT_0166;
    const ccomplex_t IT_0178 = -IT_0174 + -IT_0175 + -IT_0176 + -IT_0177;
    const ccomplex_t IT_0179 = IT_0173 + IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0168*IT_0181;
    const ccomplex_t IT_0183 = m_N_1*IT_0182;
    const ccomplex_t IT_0184 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0185 = IT_0040*IT_0184;
    const ccomplex_t IT_0186 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0187 = IT_0043*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0189 = IT_0040*IT_0188;
    const ccomplex_t IT_0190 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0191 = IT_0043*IT_0190;
    const ccomplex_t IT_0192 = IT_0185 + IT_0187 + IT_0189 + IT_0191;
    const ccomplex_t IT_0193 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0194 = IT_0046*IT_0193;
    const ccomplex_t IT_0195 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0196 = IT_0049*IT_0195;
    const ccomplex_t IT_0197 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0198 = IT_0046*IT_0197;
    const ccomplex_t IT_0199 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0200 = IT_0049*IT_0199;
    const ccomplex_t IT_0201 = -IT_0194 + -IT_0196 + -IT_0198 + -IT_0200;
    const ccomplex_t IT_0202 = IT_0192 + IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = IT_0168*IT_0204;
    const ccomplex_t IT_0206 = m_N_3*IT_0205;
    const ccomplex_t IT_0207 = IT_0183 + IT_0206;
    const ccomplex_t IT_0208 = cpow((-2)*s_23 + IT_0000 + IT_0032 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = IT_0207*IT_0209;
    const ccomplex_t IT_0211 = cpow((-2)*s_13 + IT_0000 + IT_0032 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0214 = IT_0004*IT_0213;
    const ccomplex_t IT_0215 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0216 = IT_0008*IT_0215;
    const ccomplex_t IT_0217 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0218 = IT_0013*IT_0217;
    const ccomplex_t IT_0219 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0220 = IT_0016*IT_0219;
    const ccomplex_t IT_0221 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0222 = IT_0004*IT_0221;
    const ccomplex_t IT_0223 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0224 = IT_0008*IT_0223;
    const ccomplex_t IT_0225 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0226 = IT_0013*IT_0225;
    const ccomplex_t IT_0227 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0228 = IT_0016*IT_0227;
    const ccomplex_t IT_0229 = (-0.5)*IT_0214 + 0.5*IT_0216 + 0.5*IT_0218 + (
      -0.5)*IT_0220 + (-0.5)*IT_0222 + 0.5*IT_0224 + 0.5*IT_0226 + (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = IT_0040*IT_0213;
    const ccomplex_t IT_0231 = IT_0043*IT_0215;
    const ccomplex_t IT_0232 = IT_0040*IT_0221;
    const ccomplex_t IT_0233 = IT_0043*IT_0223;
    const ccomplex_t IT_0234 = IT_0230 + IT_0231 + IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = IT_0046*IT_0217;
    const ccomplex_t IT_0236 = IT_0049*IT_0219;
    const ccomplex_t IT_0237 = IT_0046*IT_0225;
    const ccomplex_t IT_0238 = IT_0049*IT_0227;
    const ccomplex_t IT_0239 = -IT_0235 + -IT_0236 + -IT_0237 + -IT_0238;
    const ccomplex_t IT_0240 = IT_0234 + IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = (-0.5)*IT_0241;
    const ccomplex_t IT_0243 = IT_0229*IT_0242;
    const ccomplex_t IT_0244 = m_N_2*IT_0243;
    const ccomplex_t IT_0245 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0246 = IT_0040*IT_0245;
    const ccomplex_t IT_0247 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0248 = IT_0043*IT_0247;
    const ccomplex_t IT_0249 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0250 = IT_0040*IT_0249;
    const ccomplex_t IT_0251 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0252 = IT_0043*IT_0251;
    const ccomplex_t IT_0253 = IT_0246 + IT_0248 + IT_0250 + IT_0252;
    const ccomplex_t IT_0254 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0255 = IT_0046*IT_0254;
    const ccomplex_t IT_0256 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0257 = IT_0049*IT_0256;
    const ccomplex_t IT_0258 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0259 = IT_0046*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0261 = IT_0049*IT_0260;
    const ccomplex_t IT_0262 = -IT_0255 + -IT_0257 + -IT_0259 + -IT_0261;
    const ccomplex_t IT_0263 = IT_0253 + IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = IT_0229*IT_0265;
    const ccomplex_t IT_0267 = m_N_3*IT_0266;
    const ccomplex_t IT_0268 = IT_0244 + IT_0267;
    const ccomplex_t IT_0269 = IT_0212*IT_0268;
    const ccomplex_t IT_0270 = cpow((-2)*s_23 + IT_0000 + IT_0032 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0268*IT_0271;
    const ccomplex_t IT_0273 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = e_em*IT_0003*IT_0012*(IT_0001*IT_0023 + IT_0007
      *IT_0024);
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0003*IT_0011;
    const ccomplex_t IT_0278 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0279 = IT_0277*IT_0278;
    const ccomplex_t IT_0280 = IT_0002*IT_0012;
    const ccomplex_t IT_0281 = IT_0278*IT_0280;
    const ccomplex_t IT_0282 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0283 = IT_0277*IT_0282;
    const ccomplex_t IT_0284 = IT_0280*IT_0282;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*(IT_0279 + IT_0281 + 
      -IT_0283 + -IT_0284);
    const ccomplex_t IT_0286 = 0.5*IT_0285;
    const ccomplex_t IT_0287 = IT_0276*IT_0286;
    const ccomplex_t IT_0288 = m_N_3*IT_0287;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = (-0.5)*IT_0285;
    const ccomplex_t IT_0291 = IT_0276*IT_0290;
    const ccomplex_t IT_0292 = m_N_3*IT_0291;
    const ccomplex_t IT_0293 = IT_0289 + IT_0292;
    const ccomplex_t IT_0294 = IT_0274*IT_0293;
    const ccomplex_t IT_0295 = cpow((-2)*s_13 + IT_0000 + IT_0032 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = IT_0207*IT_0296;
    const ccomplex_t IT_0298 = IT_0004*IT_0044;
    const ccomplex_t IT_0299 = IT_0008*IT_0041;
    const ccomplex_t IT_0300 = IT_0013*IT_0050;
    const ccomplex_t IT_0301 = IT_0016*IT_0047;
    const ccomplex_t IT_0302 = -IT_0298 + -IT_0299 + IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0030*IT_0302;
    const ccomplex_t IT_0304 = (-0.5)*IT_0142 + IT_0150 + 0.5*IT_0151 +
       IT_0210 + IT_0269 + IT_0272 + -IT_0294 + IT_0297 + (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0304;
    const ccomplex_t IT_0306 = IT_0004*IT_0108;
    const ccomplex_t IT_0307 = IT_0008*IT_0110;
    const ccomplex_t IT_0308 = IT_0013*IT_0117;
    const ccomplex_t IT_0309 = IT_0016*IT_0119;
    const ccomplex_t IT_0310 = IT_0004*IT_0112;
    const ccomplex_t IT_0311 = IT_0008*IT_0114;
    const ccomplex_t IT_0312 = IT_0013*IT_0121;
    const ccomplex_t IT_0313 = IT_0016*IT_0123;
    const ccomplex_t IT_0314 = (-0.5)*IT_0306 + 0.5*IT_0307 + 0.5*IT_0308 + (
      -0.5)*IT_0309 + (-0.5)*IT_0310 + 0.5*IT_0311 + 0.5*IT_0312 + (-0.5)*IT_0313;
    const ccomplex_t IT_0315 = IT_0105*IT_0314;
    const ccomplex_t IT_0316 = m_N_3*IT_0315;
    const ccomplex_t IT_0317 = IT_0128*IT_0314;
    const ccomplex_t IT_0318 = m_N_4*IT_0317;
    const ccomplex_t IT_0319 = IT_0316 + IT_0318;
    const ccomplex_t IT_0320 = IT_0075*IT_0319;
    const ccomplex_t IT_0321 = IT_0134*IT_0319;
    const ccomplex_t IT_0322 = 0.5*IT_0142 + (-0.5)*IT_0151 + 0.5*IT_0303 + 
      -IT_0320 + -IT_0321;
    const ccomplex_t IT_0323 = (-2)*IT_0322;
    const ccomplex_t IT_0324 = IT_0305 + IT_0323;
    const ccomplex_t IT_0325 = 2*conj(IT_0304);
    const ccomplex_t IT_0326 = (-2)*conj(IT_0322);
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = IT_0004*IT_0245;
    const ccomplex_t IT_0329 = IT_0008*IT_0247;
    const ccomplex_t IT_0330 = IT_0013*IT_0254;
    const ccomplex_t IT_0331 = IT_0016*IT_0256;
    const ccomplex_t IT_0332 = IT_0004*IT_0249;
    const ccomplex_t IT_0333 = IT_0008*IT_0251;
    const ccomplex_t IT_0334 = IT_0013*IT_0258;
    const ccomplex_t IT_0335 = IT_0016*IT_0260;
    const ccomplex_t IT_0336 = 0.5*IT_0328 + (-0.5)*IT_0329 + (-0.5)*IT_0330 +
       0.5*IT_0331 + 0.5*IT_0332 + (-0.5)*IT_0333 + (-0.5)*IT_0334 + 0.5*IT_0335;
    const ccomplex_t IT_0337 = IT_0242*IT_0336;
    const ccomplex_t IT_0338 = m_N_3*IT_0337;
    const ccomplex_t IT_0339 = IT_0265*IT_0336;
    const ccomplex_t IT_0340 = m_N_2*IT_0339;
    const ccomplex_t IT_0341 = IT_0338 + IT_0340;
    const ccomplex_t IT_0342 = IT_0212*IT_0341;
    const ccomplex_t IT_0343 = IT_0271*IT_0341;
    const ccomplex_t IT_0344 = IT_0004*IT_0184;
    const ccomplex_t IT_0345 = IT_0008*IT_0186;
    const ccomplex_t IT_0346 = IT_0013*IT_0193;
    const ccomplex_t IT_0347 = IT_0016*IT_0195;
    const ccomplex_t IT_0348 = IT_0004*IT_0188;
    const ccomplex_t IT_0349 = IT_0008*IT_0190;
    const ccomplex_t IT_0350 = IT_0013*IT_0197;
    const ccomplex_t IT_0351 = IT_0016*IT_0199;
    const ccomplex_t IT_0352 = 0.5*IT_0344 + (-0.5)*IT_0345 + (-0.5)*IT_0346 +
       0.5*IT_0347 + 0.5*IT_0348 + (-0.5)*IT_0349 + (-0.5)*IT_0350 + 0.5*IT_0351;
    const ccomplex_t IT_0353 = IT_0181*IT_0352;
    const ccomplex_t IT_0354 = m_N_3*IT_0353;
    const ccomplex_t IT_0355 = IT_0204*IT_0352;
    const ccomplex_t IT_0356 = m_N_1*IT_0355;
    const ccomplex_t IT_0357 = IT_0354 + IT_0356;
    const ccomplex_t IT_0358 = IT_0296*IT_0357;
    const ccomplex_t IT_0359 = IT_0064*IT_0144;
    const ccomplex_t IT_0360 = IT_0209*IT_0357;
    const ccomplex_t IT_0361 = (-0.5)*IT_0031 + 0.5*IT_0065 + (-0.5)*IT_0073 +
       IT_0294 + IT_0342 + IT_0343 + IT_0358 + IT_0359 + IT_0360;
    const ccomplex_t IT_0362 = 2*IT_0322;
    const ccomplex_t IT_0363 = (-2)*IT_0304;
    const ccomplex_t IT_0364 = IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = 2*conj(IT_0322);
    const ccomplex_t IT_0366 = (-2)*conj(IT_0304);
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = 2*s_12;
    const ccomplex_t IT_0369 = (-2)*s_12;
    const ccomplex_t IT_0370 = s_23*m_N_3;
    const ccomplex_t IT_0371 = IT_0209*IT_0353;
    const ccomplex_t IT_0372 = IT_0212*IT_0266;
    const ccomplex_t IT_0373 = IT_0271*IT_0337;
    const ccomplex_t IT_0374 = IT_0144*IT_0147;
    const ccomplex_t IT_0375 = -IT_0275;
    const ccomplex_t IT_0376 = IT_0290*IT_0375;
    const ccomplex_t IT_0377 = IT_0274*IT_0376;
    const ccomplex_t IT_0378 = IT_0034*IT_0054;
    const ccomplex_t IT_0379 = IT_0205*IT_0296;
    const ccomplex_t IT_0380 = -IT_0371 + IT_0372 + -IT_0373 + IT_0374 + 
      -IT_0377 + (-0.5)*IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = IT_0075*IT_0315;
    const ccomplex_t IT_0382 = IT_0129*IT_0134;
    const ccomplex_t IT_0383 = (-0.5)*IT_0378 + IT_0381 + -IT_0382;
    const ccomplex_t IT_0384 = conj(IT_0380) + conj(IT_0383);
    const ccomplex_t IT_0385 = IT_0380 + IT_0383;
    const ccomplex_t IT_0386 = IT_0296*IT_0353;
    const ccomplex_t IT_0387 = IT_0205*IT_0209;
    const ccomplex_t IT_0388 = IT_0212*IT_0337;
    const ccomplex_t IT_0389 = IT_0266*IT_0271;
    const ccomplex_t IT_0390 = IT_0034*IT_0147;
    const ccomplex_t IT_0391 = IT_0054*IT_0144;
    const ccomplex_t IT_0392 = IT_0286*IT_0375;
    const ccomplex_t IT_0393 = IT_0274*IT_0392;
    const ccomplex_t IT_0394 = IT_0386 + -IT_0387 + IT_0388 + -IT_0389 + (-0.5
      )*IT_0390 + IT_0391 + -IT_0393;
    const ccomplex_t IT_0395 = IT_0075*IT_0129;
    const ccomplex_t IT_0396 = IT_0134*IT_0315;
    const ccomplex_t IT_0397 = (-0.5)*IT_0390 + IT_0395 + -IT_0396;
    const ccomplex_t IT_0398 = conj(IT_0394) + conj(IT_0397);
    const ccomplex_t IT_0399 = 2*IT_0136;
    const ccomplex_t IT_0400 = (-2)*IT_0361;
    const ccomplex_t IT_0401 = IT_0394 + IT_0397;
    const ccomplex_t IT_0402 = 2*conj(IT_0136);
    const ccomplex_t IT_0403 = (-2)*conj(IT_0361);
    const ccomplex_t IT_0404 = s_13*m_N_3;
    const ccomplex_t IT_0405 = 2*IT_0361;
    const ccomplex_t IT_0406 = (-2)*IT_0136;
    const ccomplex_t IT_0407 = 2*conj(IT_0361);
    const ccomplex_t IT_0408 = (-2)*conj(IT_0136);
    const ccomplex_t IT_0409 = s_12*IT_0032;
    const ccomplex_t IT_0410 = s_13*s_23;
    const ccomplex_t IT_0411 = (-2)*IT_0410;
    const ccomplex_t IT_0412 = IT_0409 + IT_0411;
    const ccomplex_t IT_0413 = (-2)*IT_0380;
    const ccomplex_t IT_0414 = (-2)*conj(IT_0380);
    const ccomplex_t IT_0415 = (-2)*IT_0394;
    const ccomplex_t IT_0416 = (-2)*conj(IT_0394);
    const ccomplex_t IT_0417 = (-2)*IT_0397;
    const ccomplex_t IT_0418 = IT_0000*IT_0032;
    const ccomplex_t IT_0419 = IT_0000*(conj(IT_0136)*IT_0324 + IT_0136
      *IT_0327 + conj(IT_0361)*IT_0364 + IT_0361*IT_0367) + (IT_0136*conj
      (IT_0136) + IT_0304*conj(IT_0304) + IT_0322*conj(IT_0322) + IT_0361*conj
      (IT_0361))*IT_0368 + (conj(IT_0304)*IT_0322 + IT_0304*conj(IT_0322) + conj
      (IT_0136)*IT_0361 + IT_0136*conj(IT_0361))*IT_0369 + IT_0370*(IT_0364
      *IT_0384 + IT_0367*IT_0385 + IT_0398*(IT_0399 + IT_0400) + IT_0401*
      (IT_0402 + IT_0403)) + IT_0404*(IT_0324*IT_0398 + IT_0327*IT_0401 +
       IT_0384*(IT_0405 + IT_0406) + IT_0385*(IT_0407 + IT_0408)) + IT_0412*
      (IT_0384*IT_0413 + (-2)*IT_0383*(conj(IT_0383) + (-0.5)*IT_0414) + IT_0398
      *IT_0415 + IT_0397*IT_0416 + conj(IT_0397)*IT_0417) + (conj(IT_0394)
      *IT_0413 + (-2)*conj(IT_0397)*(IT_0383 + (-0.5)*IT_0413) + IT_0401*IT_0414
       + IT_0383*IT_0416 + conj(IT_0383)*(IT_0415 + IT_0417))*IT_0418;
    return create_ccomplex_return(IT_0419);
}

