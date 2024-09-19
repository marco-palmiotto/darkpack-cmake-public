#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_3_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_3_to_A0_Z(
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
    const ccomplex_t IT_0000 = m_N_2*m_N_3;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_2, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = pow(m_A0, 2);
    const ccomplex_t IT_0009 = pow(m_N_3, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = sin(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0012*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = IT_0018*IT_0023;
    const ccomplex_t IT_0028 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0031 = IT_0015*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0033 = IT_0019*IT_0032;
    const ccomplex_t IT_0034 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0035 = IT_0024*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0017 + (-0.5)*IT_0021 + (-0.5)*IT_0026 +
       0.5*IT_0029 + 0.5*IT_0031 + (-0.5)*IT_0033 + (-0.5)*IT_0035 + 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0013*IT_0023;
    const ccomplex_t IT_0040 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = IT_0014*IT_0022;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = IT_0039*IT_0044;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0041 + IT_0043 + 
      -IT_0045 + -IT_0046);
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = IT_0038*IT_0048;
    const ccomplex_t IT_0050 = IT_0011*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0054 = IT_0015*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0056 = IT_0019*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0058 = IT_0024*IT_0057;
    const ccomplex_t IT_0059 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0060 = IT_0027*IT_0059;
    const ccomplex_t IT_0061 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0062 = IT_0015*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0064 = IT_0019*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0066 = IT_0024*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0068 = IT_0027*IT_0067;
    const ccomplex_t IT_0069 = 0.5*IT_0054 + (-0.5)*IT_0056 + (-0.5)*IT_0058 +
       0.5*IT_0060 + 0.5*IT_0062 + (-0.5)*IT_0064 + (-0.5)*IT_0066 + 0.5*IT_0068;
    const ccomplex_t IT_0070 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0071 = IT_0042*IT_0070;
    const ccomplex_t IT_0072 = IT_0039*IT_0070;
    const ccomplex_t IT_0073 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0074 = IT_0039*IT_0073;
    const ccomplex_t IT_0075 = IT_0042*IT_0073;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0071 + IT_0072 + 
      -IT_0074 + -IT_0075);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0069*IT_0080;
    const ccomplex_t IT_0082 = IT_0052*IT_0081;
    const ccomplex_t IT_0083 = -IT_0050 + (-2)*IT_0082;
    const ccomplex_t IT_0084 = cos(alpha);
    const ccomplex_t IT_0085 = IT_0014*IT_0084;
    const ccomplex_t IT_0086 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = sin(alpha);
    const ccomplex_t IT_0089 = IT_0014*IT_0088;
    const ccomplex_t IT_0090 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0093 = IT_0085*IT_0092;
    const ccomplex_t IT_0094 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0095 = IT_0089*IT_0094;
    const ccomplex_t IT_0096 = IT_0087 + IT_0091 + IT_0093 + IT_0095;
    const ccomplex_t IT_0097 = IT_0023*IT_0084;
    const ccomplex_t IT_0098 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = IT_0023*IT_0088;
    const ccomplex_t IT_0101 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0104 = IT_0097*IT_0103;
    const ccomplex_t IT_0105 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0106 = IT_0100*IT_0105;
    const ccomplex_t IT_0107 = -IT_0099 + -IT_0102 + -IT_0104 + -IT_0106;
    const ccomplex_t IT_0108 = IT_0096 + IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = e_em*IT_0014*IT_0023*(IT_0012*IT_0084 + IT_0018
      *IT_0088);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0004 + -IT_0009 + -reg_prop, -1);
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_12 + -IT_0004 + -IT_0009 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0085*IT_0090;
    const ccomplex_t IT_0120 = IT_0098*IT_0100;
    const ccomplex_t IT_0121 = IT_0085*IT_0094;
    const ccomplex_t IT_0122 = IT_0100*IT_0103;
    const ccomplex_t IT_0123 = IT_0119 + IT_0120 + IT_0121 + IT_0122;
    const ccomplex_t IT_0124 = IT_0086*IT_0089;
    const ccomplex_t IT_0125 = IT_0097*IT_0101;
    const ccomplex_t IT_0126 = IT_0089*IT_0092;
    const ccomplex_t IT_0127 = IT_0097*IT_0105;
    const ccomplex_t IT_0128 = -IT_0124 + -IT_0125 + -IT_0126 + -IT_0127;
    const ccomplex_t IT_0129 = IT_0123 + IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = e_em*IT_0014*IT_0023*(IT_0018*IT_0084 + 
      -IT_0012*IT_0088);
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = IT_0118*IT_0133;
    const ccomplex_t IT_0135 = IT_0015*IT_0086;
    const ccomplex_t IT_0136 = IT_0019*IT_0090;
    const ccomplex_t IT_0137 = IT_0024*IT_0098;
    const ccomplex_t IT_0138 = IT_0027*IT_0101;
    const ccomplex_t IT_0139 = IT_0015*IT_0092;
    const ccomplex_t IT_0140 = IT_0019*IT_0094;
    const ccomplex_t IT_0141 = IT_0024*IT_0103;
    const ccomplex_t IT_0142 = IT_0027*IT_0105;
    const ccomplex_t IT_0143 = (-0.5)*IT_0135 + 0.5*IT_0136 + 0.5*IT_0137 + (
      -0.5)*IT_0138 + (-0.5)*IT_0139 + 0.5*IT_0140 + 0.5*IT_0141 + (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = IT_0048*IT_0143;
    const ccomplex_t IT_0145 = IT_0011*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_23 + IT_0008 + IT_0009 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0149 = IT_0015*IT_0148;
    const ccomplex_t IT_0150 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0151 = IT_0019*IT_0150;
    const ccomplex_t IT_0152 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0153 = IT_0024*IT_0152;
    const ccomplex_t IT_0154 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0155 = IT_0027*IT_0154;
    const ccomplex_t IT_0156 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0157 = IT_0015*IT_0156;
    const ccomplex_t IT_0158 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0159 = IT_0019*IT_0158;
    const ccomplex_t IT_0160 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0161 = IT_0024*IT_0160;
    const ccomplex_t IT_0162 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0163 = IT_0027*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0149 + 0.5*IT_0151 + 0.5*IT_0153 + (
      -0.5)*IT_0155 + (-0.5)*IT_0157 + 0.5*IT_0159 + 0.5*IT_0161 + (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0166 = IT_0042*IT_0165;
    const ccomplex_t IT_0167 = IT_0039*IT_0165;
    const ccomplex_t IT_0168 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0169 = IT_0042*IT_0168;
    const ccomplex_t IT_0170 = IT_0039*IT_0168;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0166 + IT_0167 + 
      -IT_0169 + -IT_0170);
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = IT_0164*IT_0172;
    const ccomplex_t IT_0174 = IT_0147*IT_0173;
    const ccomplex_t IT_0175 = cpow((-2)*s_23 + IT_0008 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0178 = IT_0015*IT_0177;
    const ccomplex_t IT_0179 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0180 = IT_0019*IT_0179;
    const ccomplex_t IT_0181 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0182 = IT_0024*IT_0181;
    const ccomplex_t IT_0183 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0184 = IT_0027*IT_0183;
    const ccomplex_t IT_0185 = -IT_0178 + IT_0180 + IT_0182 + -IT_0184;
    const ccomplex_t IT_0186 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0187 = IT_0042*IT_0186;
    const ccomplex_t IT_0188 = IT_0039*IT_0186;
    const ccomplex_t IT_0189 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0190 = IT_0042*IT_0189;
    const ccomplex_t IT_0191 = IT_0039*IT_0189;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0187 + IT_0188 + 
      -IT_0190 + -IT_0191);
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = IT_0185*IT_0194;
    const ccomplex_t IT_0196 = IT_0176*IT_0195;
    const ccomplex_t IT_0197 = -IT_0116 + -IT_0134 + -IT_0145 + (-2)*IT_0174 +
       2*IT_0196;
    const ccomplex_t IT_0198 = 2*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0047;
    const ccomplex_t IT_0200 = IT_0143*IT_0199;
    const ccomplex_t IT_0201 = IT_0011*IT_0200;
    const ccomplex_t IT_0202 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0203 = IT_0015*IT_0202;
    const ccomplex_t IT_0204 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0205 = IT_0019*IT_0204;
    const ccomplex_t IT_0206 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0207 = IT_0024*IT_0206;
    const ccomplex_t IT_0208 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0209 = IT_0027*IT_0208;
    const ccomplex_t IT_0210 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0211 = IT_0015*IT_0210;
    const ccomplex_t IT_0212 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0213 = IT_0019*IT_0212;
    const ccomplex_t IT_0214 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0215 = IT_0024*IT_0214;
    const ccomplex_t IT_0216 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0217 = IT_0027*IT_0216;
    const ccomplex_t IT_0218 = (-0.5)*IT_0203 + 0.5*IT_0205 + 0.5*IT_0207 + (
      -0.5)*IT_0209 + (-0.5)*IT_0211 + 0.5*IT_0213 + 0.5*IT_0215 + (-0.5)*IT_0217;
    const ccomplex_t IT_0219 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0220 = IT_0042*IT_0219;
    const ccomplex_t IT_0221 = IT_0039*IT_0219;
    const ccomplex_t IT_0222 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0223 = IT_0039*IT_0222;
    const ccomplex_t IT_0224 = IT_0042*IT_0222;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0220 + IT_0221 + 
      -IT_0223 + -IT_0224);
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = IT_0218*IT_0227;
    const ccomplex_t IT_0229 = IT_0052*IT_0228;
    const ccomplex_t IT_0230 = -IT_0201 + (-2)*IT_0229;
    const ccomplex_t IT_0231 = (-2)*IT_0230;
    const ccomplex_t IT_0232 = IT_0198 + IT_0231;
    const ccomplex_t IT_0233 = 2*conj(IT_0197);
    const ccomplex_t IT_0234 = (-2)*conj(IT_0230);
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = IT_0020*IT_0085;
    const ccomplex_t IT_0237 = IT_0025*IT_0100;
    const ccomplex_t IT_0238 = IT_0032*IT_0085;
    const ccomplex_t IT_0239 = IT_0034*IT_0100;
    const ccomplex_t IT_0240 = IT_0236 + IT_0237 + IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = IT_0016*IT_0089;
    const ccomplex_t IT_0242 = IT_0028*IT_0097;
    const ccomplex_t IT_0243 = IT_0030*IT_0089;
    const ccomplex_t IT_0244 = IT_0036*IT_0097;
    const ccomplex_t IT_0245 = -IT_0241 + -IT_0242 + -IT_0243 + -IT_0244;
    const ccomplex_t IT_0246 = IT_0240 + IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = 0.5*IT_0247;
    const ccomplex_t IT_0249 = IT_0132*IT_0248;
    const ccomplex_t IT_0250 = IT_0118*IT_0249;
    const ccomplex_t IT_0251 = IT_0016*IT_0085;
    const ccomplex_t IT_0252 = IT_0020*IT_0089;
    const ccomplex_t IT_0253 = IT_0030*IT_0085;
    const ccomplex_t IT_0254 = IT_0032*IT_0089;
    const ccomplex_t IT_0255 = IT_0251 + IT_0252 + IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = IT_0025*IT_0097;
    const ccomplex_t IT_0257 = IT_0028*IT_0100;
    const ccomplex_t IT_0258 = IT_0034*IT_0097;
    const ccomplex_t IT_0259 = IT_0036*IT_0100;
    const ccomplex_t IT_0260 = -IT_0256 + -IT_0257 + -IT_0258 + -IT_0259;
    const ccomplex_t IT_0261 = IT_0255 + IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = (-0.5)*IT_0262;
    const ccomplex_t IT_0264 = IT_0112*IT_0263;
    const ccomplex_t IT_0265 = IT_0115*IT_0264;
    const ccomplex_t IT_0266 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0267 = IT_0015*IT_0266;
    const ccomplex_t IT_0268 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0269 = IT_0019*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0271 = IT_0024*IT_0270;
    const ccomplex_t IT_0272 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0273 = IT_0027*IT_0272;
    const ccomplex_t IT_0274 = IT_0267 + -IT_0269 + -IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0276 = IT_0039*IT_0275;
    const ccomplex_t IT_0277 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0278 = IT_0042*IT_0277;
    const ccomplex_t IT_0279 = IT_0039*IT_0277;
    const ccomplex_t IT_0280 = IT_0042*IT_0275;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*(IT_0276 + -IT_0278 + 
      -IT_0279 + IT_0280);
    const ccomplex_t IT_0282 = (-0.5)*IT_0281;
    const ccomplex_t IT_0283 = IT_0274*IT_0282;
    const ccomplex_t IT_0284 = IT_0176*IT_0283;
    const ccomplex_t IT_0285 = IT_0038*IT_0199;
    const ccomplex_t IT_0286 = IT_0011*IT_0285;
    const ccomplex_t IT_0287 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0288 = IT_0015*IT_0287;
    const ccomplex_t IT_0289 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0290 = IT_0019*IT_0289;
    const ccomplex_t IT_0291 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0292 = IT_0024*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0294 = IT_0027*IT_0293;
    const ccomplex_t IT_0295 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0296 = IT_0015*IT_0295;
    const ccomplex_t IT_0297 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0298 = IT_0019*IT_0297;
    const ccomplex_t IT_0299 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0300 = IT_0024*IT_0299;
    const ccomplex_t IT_0301 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0302 = IT_0027*IT_0301;
    const ccomplex_t IT_0303 = 0.5*IT_0288 + (-0.5)*IT_0290 + (-0.5)*IT_0292 +
       0.5*IT_0294 + 0.5*IT_0296 + (-0.5)*IT_0298 + (-0.5)*IT_0300 + 0.5*IT_0302;
    const ccomplex_t IT_0304 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0305 = IT_0039*IT_0304;
    const ccomplex_t IT_0306 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0307 = IT_0039*IT_0306;
    const ccomplex_t IT_0308 = IT_0042*IT_0304;
    const ccomplex_t IT_0309 = IT_0042*IT_0306;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*(IT_0305 + -IT_0307 +
       IT_0308 + -IT_0309);
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = (-0.5)*IT_0311;
    const ccomplex_t IT_0313 = IT_0303*IT_0312;
    const ccomplex_t IT_0314 = IT_0147*IT_0313;
    const ccomplex_t IT_0315 = -IT_0250 + -IT_0265 + 2*IT_0284 + -IT_0286 + (
      -2)*IT_0314;
    const ccomplex_t IT_0316 = 2*IT_0230;
    const ccomplex_t IT_0317 = (-2)*IT_0197;
    const ccomplex_t IT_0318 = IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = 2*conj(IT_0230);
    const ccomplex_t IT_0320 = (-2)*conj(IT_0197);
    const ccomplex_t IT_0321 = IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0323 = -IT_0322;
    const ccomplex_t IT_0324 = s_12 + IT_0323;
    const ccomplex_t IT_0325 = m_N_2*IT_0324;
    const ccomplex_t IT_0326 = m_N_2*IT_0200;
    const ccomplex_t IT_0327 = m_N_2*IT_0144;
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = IT_0326 + IT_0328;
    const ccomplex_t IT_0330 = IT_0011*IT_0329;
    const ccomplex_t IT_0331 = IT_0164*IT_0312;
    const ccomplex_t IT_0332 = m_N_1*IT_0331;
    const ccomplex_t IT_0333 = m_N_2*IT_0173;
    const ccomplex_t IT_0334 = -IT_0333;
    const ccomplex_t IT_0335 = IT_0332 + IT_0334;
    const ccomplex_t IT_0336 = IT_0147*IT_0335;
    const ccomplex_t IT_0337 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0340 = IT_0042*IT_0339;
    const ccomplex_t IT_0341 = IT_0039*IT_0339;
    const ccomplex_t IT_0342 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0343 = IT_0042*IT_0342;
    const ccomplex_t IT_0344 = IT_0039*IT_0342;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*(IT_0340 + IT_0341 + 
      -IT_0343 + -IT_0344);
    const ccomplex_t IT_0346 = (-0.5)*IT_0345;
    const ccomplex_t IT_0347 = IT_0038*IT_0346;
    const ccomplex_t IT_0348 = m_N_3*IT_0347;
    const ccomplex_t IT_0349 = -IT_0348;
    const ccomplex_t IT_0350 = 0.5*IT_0345;
    const ccomplex_t IT_0351 = IT_0038*IT_0350;
    const ccomplex_t IT_0352 = m_N_3*IT_0351;
    const ccomplex_t IT_0353 = IT_0349 + IT_0352;
    const ccomplex_t IT_0354 = IT_0338*IT_0353;
    const ccomplex_t IT_0355 = m_N_2*IT_0195;
    const ccomplex_t IT_0356 = IT_0176*IT_0355;
    const ccomplex_t IT_0357 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*IT_0357;
    const ccomplex_t IT_0359 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0360 = IT_0015*IT_0359;
    const ccomplex_t IT_0361 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0362 = IT_0019*IT_0361;
    const ccomplex_t IT_0363 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0364 = IT_0024*IT_0363;
    const ccomplex_t IT_0365 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0366 = IT_0027*IT_0365;
    const ccomplex_t IT_0367 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0368 = IT_0015*IT_0367;
    const ccomplex_t IT_0369 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0370 = IT_0019*IT_0369;
    const ccomplex_t IT_0371 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0372 = IT_0024*IT_0371;
    const ccomplex_t IT_0373 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0374 = IT_0027*IT_0373;
    const ccomplex_t IT_0375 = 0.5*IT_0360 + (-0.5)*IT_0362 + (-0.5)*IT_0364 +
       0.5*IT_0366 + 0.5*IT_0368 + (-0.5)*IT_0370 + (-0.5)*IT_0372 + 0.5*IT_0374;
    const ccomplex_t IT_0376 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0377 = IT_0042*IT_0376;
    const ccomplex_t IT_0378 = IT_0039*IT_0376;
    const ccomplex_t IT_0379 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0380 = IT_0042*IT_0379;
    const ccomplex_t IT_0381 = IT_0039*IT_0379;
    const ccomplex_t IT_0382 = (0 + _Complex_I*1)*(IT_0377 + IT_0378 + 
      -IT_0380 + -IT_0381);
    const ccomplex_t IT_0383 = (-0.5)*IT_0382;
    const ccomplex_t IT_0384 = IT_0375*IT_0383;
    const ccomplex_t IT_0385 = m_N_4*IT_0384;
    const ccomplex_t IT_0386 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0387 = IT_0042*IT_0386;
    const ccomplex_t IT_0388 = IT_0039*IT_0386;
    const ccomplex_t IT_0389 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0390 = IT_0042*IT_0389;
    const ccomplex_t IT_0391 = IT_0039*IT_0389;
    const ccomplex_t IT_0392 = (0 + _Complex_I*1)*(IT_0387 + IT_0388 + 
      -IT_0390 + -IT_0391);
    const ccomplex_t IT_0393 = (-0.5)*IT_0392;
    const ccomplex_t IT_0394 = IT_0375*IT_0393;
    const ccomplex_t IT_0395 = m_N_3*IT_0394;
    const ccomplex_t IT_0396 = -IT_0395;
    const ccomplex_t IT_0397 = IT_0385 + IT_0396;
    const ccomplex_t IT_0398 = IT_0358*IT_0397;
    const ccomplex_t IT_0399 = 0.5*IT_0330 + IT_0336 + (-0.5)*IT_0354 +
       IT_0356 + IT_0398;
    const ccomplex_t IT_0400 = IT_0185*IT_0282;
    const ccomplex_t IT_0401 = m_N_3*IT_0400;
    const ccomplex_t IT_0402 = IT_0176*IT_0401;
    const ccomplex_t IT_0403 = cpow((-2)*s_13 + IT_0004 + IT_0008 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0404 = (0 + _Complex_I*1)*IT_0403;
    const ccomplex_t IT_0405 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0406 = IT_0015*IT_0405;
    const ccomplex_t IT_0407 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0408 = IT_0019*IT_0407;
    const ccomplex_t IT_0409 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0410 = IT_0024*IT_0409;
    const ccomplex_t IT_0411 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0412 = IT_0027*IT_0411;
    const ccomplex_t IT_0413 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0414 = IT_0015*IT_0413;
    const ccomplex_t IT_0415 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0416 = IT_0019*IT_0415;
    const ccomplex_t IT_0417 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0418 = IT_0024*IT_0417;
    const ccomplex_t IT_0419 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0420 = IT_0027*IT_0419;
    const ccomplex_t IT_0421 = 0.5*IT_0406 + (-0.5)*IT_0408 + (-0.5)*IT_0410 +
       0.5*IT_0412 + 0.5*IT_0414 + (-0.5)*IT_0416 + (-0.5)*IT_0418 + 0.5*IT_0420;
    const ccomplex_t IT_0422 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0423 = IT_0042*IT_0422;
    const ccomplex_t IT_0424 = IT_0039*IT_0422;
    const ccomplex_t IT_0425 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0426 = IT_0042*IT_0425;
    const ccomplex_t IT_0427 = IT_0039*IT_0425;
    const ccomplex_t IT_0428 = (0 + _Complex_I*1)*(IT_0423 + IT_0424 + 
      -IT_0426 + -IT_0427);
    const ccomplex_t IT_0429 = (-0.5)*IT_0428;
    const ccomplex_t IT_0430 = IT_0421*IT_0429;
    const ccomplex_t IT_0431 = m_N_3*IT_0430;
    const ccomplex_t IT_0432 = -IT_0431;
    const ccomplex_t IT_0433 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0434 = IT_0042*IT_0433;
    const ccomplex_t IT_0435 = IT_0039*IT_0433;
    const ccomplex_t IT_0436 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0437 = IT_0042*IT_0436;
    const ccomplex_t IT_0438 = IT_0039*IT_0436;
    const ccomplex_t IT_0439 = (0 + _Complex_I*1)*(IT_0434 + IT_0435 + 
      -IT_0437 + -IT_0438);
    const ccomplex_t IT_0440 = (-0.5)*IT_0439;
    const ccomplex_t IT_0441 = IT_0421*IT_0440;
    const ccomplex_t IT_0442 = m_N_1*IT_0441;
    const ccomplex_t IT_0443 = IT_0432 + IT_0442;
    const ccomplex_t IT_0444 = IT_0404*IT_0443;
    const ccomplex_t IT_0445 = cpow((-2)*s_13 + IT_0008 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0446 = (0 + _Complex_I*1)*IT_0445;
    const ccomplex_t IT_0447 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0448 = IT_0015*IT_0447;
    const ccomplex_t IT_0449 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0450 = IT_0019*IT_0449;
    const ccomplex_t IT_0451 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0452 = IT_0024*IT_0451;
    const ccomplex_t IT_0453 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0454 = IT_0027*IT_0453;
    const ccomplex_t IT_0455 = IT_0448 + -IT_0450 + -IT_0452 + IT_0454;
    const ccomplex_t IT_0456 = IT_0282*IT_0455;
    const ccomplex_t IT_0457 = m_N_2*IT_0456;
    const ccomplex_t IT_0458 = IT_0194*IT_0455;
    const ccomplex_t IT_0459 = m_N_3*IT_0458;
    const ccomplex_t IT_0460 = -IT_0459;
    const ccomplex_t IT_0461 = IT_0457 + IT_0460;
    const ccomplex_t IT_0462 = IT_0446*IT_0461;
    const ccomplex_t IT_0463 = IT_0080*IT_0218;
    const ccomplex_t IT_0464 = m_N_4*IT_0463;
    const ccomplex_t IT_0465 = m_N_2*IT_0228;
    const ccomplex_t IT_0466 = -IT_0465;
    const ccomplex_t IT_0467 = IT_0464 + IT_0466;
    const ccomplex_t IT_0468 = IT_0052*IT_0467;
    const ccomplex_t IT_0469 = 0.5*IT_0330 + (-0.5)*IT_0354 + -IT_0402 + 
      -IT_0444 + -IT_0462 + -IT_0468;
    const ccomplex_t IT_0470 = conj(IT_0399) + conj(IT_0469);
    const ccomplex_t IT_0471 = IT_0399 + IT_0469;
    const ccomplex_t IT_0472 = m_N_2*IT_0049;
    const ccomplex_t IT_0473 = m_N_2*IT_0285;
    const ccomplex_t IT_0474 = -IT_0473;
    const ccomplex_t IT_0475 = IT_0472 + IT_0474;
    const ccomplex_t IT_0476 = IT_0011*IT_0475;
    const ccomplex_t IT_0477 = IT_0172*IT_0303;
    const ccomplex_t IT_0478 = m_N_1*IT_0477;
    const ccomplex_t IT_0479 = m_N_2*IT_0313;
    const ccomplex_t IT_0480 = -IT_0479;
    const ccomplex_t IT_0481 = IT_0478 + IT_0480;
    const ccomplex_t IT_0482 = IT_0147*IT_0481;
    const ccomplex_t IT_0483 = IT_0143*IT_0346;
    const ccomplex_t IT_0484 = m_N_3*IT_0483;
    const ccomplex_t IT_0485 = IT_0143*IT_0350;
    const ccomplex_t IT_0486 = m_N_3*IT_0485;
    const ccomplex_t IT_0487 = -IT_0486;
    const ccomplex_t IT_0488 = IT_0484 + IT_0487;
    const ccomplex_t IT_0489 = IT_0338*IT_0488;
    const ccomplex_t IT_0490 = m_N_2*IT_0283;
    const ccomplex_t IT_0491 = IT_0176*IT_0490;
    const ccomplex_t IT_0492 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0493 = IT_0015*IT_0492;
    const ccomplex_t IT_0494 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0495 = IT_0019*IT_0494;
    const ccomplex_t IT_0496 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0497 = IT_0024*IT_0496;
    const ccomplex_t IT_0498 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0499 = IT_0027*IT_0498;
    const ccomplex_t IT_0500 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0501 = IT_0015*IT_0500;
    const ccomplex_t IT_0502 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0503 = IT_0019*IT_0502;
    const ccomplex_t IT_0504 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0505 = IT_0024*IT_0504;
    const ccomplex_t IT_0506 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0507 = IT_0027*IT_0506;
    const ccomplex_t IT_0508 = (-0.5)*IT_0493 + 0.5*IT_0495 + 0.5*IT_0497 + (
      -0.5)*IT_0499 + (-0.5)*IT_0501 + 0.5*IT_0503 + 0.5*IT_0505 + (-0.5)*IT_0507;
    const ccomplex_t IT_0509 = IT_0393*IT_0508;
    const ccomplex_t IT_0510 = m_N_4*IT_0509;
    const ccomplex_t IT_0511 = IT_0383*IT_0508;
    const ccomplex_t IT_0512 = m_N_3*IT_0511;
    const ccomplex_t IT_0513 = -IT_0512;
    const ccomplex_t IT_0514 = IT_0510 + IT_0513;
    const ccomplex_t IT_0515 = IT_0358*IT_0514;
    const ccomplex_t IT_0516 = 0.5*IT_0476 + IT_0482 + (-0.5)*IT_0489 +
       IT_0491 + IT_0515;
    const ccomplex_t IT_0517 = IT_0194*IT_0274;
    const ccomplex_t IT_0518 = m_N_3*IT_0517;
    const ccomplex_t IT_0519 = IT_0176*IT_0518;
    const ccomplex_t IT_0520 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0521 = IT_0015*IT_0520;
    const ccomplex_t IT_0522 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0523 = IT_0019*IT_0522;
    const ccomplex_t IT_0524 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0525 = IT_0024*IT_0524;
    const ccomplex_t IT_0526 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0527 = IT_0027*IT_0526;
    const ccomplex_t IT_0528 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0529 = IT_0015*IT_0528;
    const ccomplex_t IT_0530 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0531 = IT_0019*IT_0530;
    const ccomplex_t IT_0532 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0533 = IT_0024*IT_0532;
    const ccomplex_t IT_0534 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0535 = IT_0027*IT_0534;
    const ccomplex_t IT_0536 = (-0.5)*IT_0521 + 0.5*IT_0523 + 0.5*IT_0525 + (
      -0.5)*IT_0527 + (-0.5)*IT_0529 + 0.5*IT_0531 + 0.5*IT_0533 + (-0.5)*IT_0535;
    const ccomplex_t IT_0537 = IT_0440*IT_0536;
    const ccomplex_t IT_0538 = m_N_3*IT_0537;
    const ccomplex_t IT_0539 = -IT_0538;
    const ccomplex_t IT_0540 = IT_0429*IT_0536;
    const ccomplex_t IT_0541 = m_N_1*IT_0540;
    const ccomplex_t IT_0542 = IT_0539 + IT_0541;
    const ccomplex_t IT_0543 = IT_0404*IT_0542;
    const ccomplex_t IT_0544 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0545 = IT_0015*IT_0544;
    const ccomplex_t IT_0546 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0547 = IT_0019*IT_0546;
    const ccomplex_t IT_0548 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0549 = IT_0024*IT_0548;
    const ccomplex_t IT_0550 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0551 = IT_0027*IT_0550;
    const ccomplex_t IT_0552 = -IT_0545 + IT_0547 + IT_0549 + -IT_0551;
    const ccomplex_t IT_0553 = IT_0194*IT_0552;
    const ccomplex_t IT_0554 = m_N_2*IT_0553;
    const ccomplex_t IT_0555 = IT_0282*IT_0552;
    const ccomplex_t IT_0556 = m_N_3*IT_0555;
    const ccomplex_t IT_0557 = -IT_0556;
    const ccomplex_t IT_0558 = IT_0554 + IT_0557;
    const ccomplex_t IT_0559 = IT_0446*IT_0558;
    const ccomplex_t IT_0560 = IT_0069*IT_0227;
    const ccomplex_t IT_0561 = m_N_4*IT_0560;
    const ccomplex_t IT_0562 = m_N_2*IT_0081;
    const ccomplex_t IT_0563 = -IT_0562;
    const ccomplex_t IT_0564 = IT_0561 + IT_0563;
    const ccomplex_t IT_0565 = IT_0052*IT_0564;
    const ccomplex_t IT_0566 = 0.5*IT_0476 + (-0.5)*IT_0489 + -IT_0519 + 
      -IT_0543 + -IT_0559 + -IT_0565;
    const ccomplex_t IT_0567 = conj(IT_0516) + conj(IT_0566);
    const ccomplex_t IT_0568 = 2*IT_0315;
    const ccomplex_t IT_0569 = (-2)*IT_0083;
    const ccomplex_t IT_0570 = IT_0568 + IT_0569;
    const ccomplex_t IT_0571 = IT_0516 + IT_0566;
    const ccomplex_t IT_0572 = 2*conj(IT_0315);
    const ccomplex_t IT_0573 = (-2)*conj(IT_0083);
    const ccomplex_t IT_0574 = IT_0572 + IT_0573;
    const ccomplex_t IT_0575 = s_12*IT_0324;
    const ccomplex_t IT_0576 = IT_0446*IT_0555;
    const ccomplex_t IT_0577 = IT_0338*IT_0485;
    const ccomplex_t IT_0578 = IT_0404*IT_0537;
    const ccomplex_t IT_0579 = 2*IT_0576 + IT_0577 + 2*IT_0578;
    const ccomplex_t IT_0580 = IT_0338*IT_0483;
    const ccomplex_t IT_0581 = IT_0358*IT_0511;
    const ccomplex_t IT_0582 = -IT_0116 + -IT_0134 + IT_0580 + 2*IT_0581;
    const ccomplex_t IT_0583 = IT_0338*IT_0347;
    const ccomplex_t IT_0584 = IT_0404*IT_0430;
    const ccomplex_t IT_0585 = IT_0446*IT_0458;
    const ccomplex_t IT_0586 = IT_0583 + 2*IT_0584 + 2*IT_0585;
    const ccomplex_t IT_0587 = IT_0338*IT_0351;
    const ccomplex_t IT_0588 = IT_0358*IT_0394;
    const ccomplex_t IT_0589 = -IT_0250 + -IT_0265 + IT_0587 + 2*IT_0588;
    const ccomplex_t IT_0590 = 2*IT_0083;
    const ccomplex_t IT_0591 = (-2)*IT_0315;
    const ccomplex_t IT_0592 = IT_0590 + IT_0591;
    const ccomplex_t IT_0593 = 2*conj(IT_0083);
    const ccomplex_t IT_0594 = (-2)*conj(IT_0315);
    const ccomplex_t IT_0595 = IT_0593 + IT_0594;
    const ccomplex_t IT_0596 = s_12*s_14;
    const ccomplex_t IT_0597 = s_24*IT_0004;
    const ccomplex_t IT_0598 = -IT_0597;
    const ccomplex_t IT_0599 = IT_0596 + IT_0598;
    const ccomplex_t IT_0600 = 0.5*IT_0145 + IT_0174 + -IT_0196;
    const ccomplex_t IT_0601 = -IT_0576 + (-0.5)*IT_0577 + -IT_0578;
    const ccomplex_t IT_0602 = 0.5*IT_0201 + IT_0229;
    const ccomplex_t IT_0603 = -IT_0581;
    const ccomplex_t IT_0604 = (-0.5)*IT_0580;
    const ccomplex_t IT_0605 = conj(IT_0603) + conj(IT_0604);
    const ccomplex_t IT_0606 = conj(IT_0602) + IT_0605;
    const ccomplex_t IT_0607 = IT_0603 + IT_0604;
    const ccomplex_t IT_0608 = IT_0602 + IT_0607;
    const ccomplex_t IT_0609 = -IT_0284 + 0.5*IT_0286 + IT_0314;
    const ccomplex_t IT_0610 = (-0.5)*IT_0583 + -IT_0584 + -IT_0585;
    const ccomplex_t IT_0611 = conj(IT_0609) + conj(IT_0610);
    const ccomplex_t IT_0612 = IT_0609 + IT_0610;
    const ccomplex_t IT_0613 = 0.5*IT_0050;
    const ccomplex_t IT_0614 = -IT_0588;
    const ccomplex_t IT_0615 = (-0.5)*IT_0587;
    const ccomplex_t IT_0616 = conj(IT_0614) + conj(IT_0615);
    const ccomplex_t IT_0617 = IT_0614 + IT_0615;
    const ccomplex_t IT_0618 = IT_0235*(IT_0600 + IT_0601) + IT_0232*(conj
      (IT_0600) + conj(IT_0601)) + IT_0318*IT_0606 + IT_0321*IT_0608 + IT_0570
      *IT_0611 + IT_0574*IT_0612 + IT_0592*(conj(IT_0082) + conj(IT_0613) +
       IT_0616) + IT_0595*(IT_0082 + IT_0613 + IT_0617);
    const ccomplex_t IT_0619 = IT_0000*IT_0324;
    const ccomplex_t IT_0620 = 2*IT_0582;
    const ccomplex_t IT_0621 = (-2)*IT_0579;
    const ccomplex_t IT_0622 = IT_0620 + IT_0621;
    const ccomplex_t IT_0623 = 2*conj(IT_0582);
    const ccomplex_t IT_0624 = (-2)*conj(IT_0579);
    const ccomplex_t IT_0625 = IT_0623 + IT_0624;
    const ccomplex_t IT_0626 = 2*IT_0579;
    const ccomplex_t IT_0627 = (-2)*IT_0582;
    const ccomplex_t IT_0628 = IT_0626 + IT_0627;
    const ccomplex_t IT_0629 = 2*conj(IT_0579);
    const ccomplex_t IT_0630 = (-2)*conj(IT_0582);
    const ccomplex_t IT_0631 = IT_0629 + IT_0630;
    const ccomplex_t IT_0632 = pow(s_24, 2);
    const ccomplex_t IT_0633 = IT_0001*IT_0632;
    const ccomplex_t IT_0634 = -IT_0009;
    const ccomplex_t IT_0635 = IT_0633 + IT_0634;
    const ccomplex_t IT_0636 = IT_0000*IT_0635;
    const ccomplex_t IT_0637 = m_N_2*IT_0635;
    const ccomplex_t IT_0638 = 2*IT_0586;
    const ccomplex_t IT_0639 = (-2)*IT_0589;
    const ccomplex_t IT_0640 = IT_0638 + IT_0639;
    const ccomplex_t IT_0641 = 2*conj(IT_0586);
    const ccomplex_t IT_0642 = (-2)*conj(IT_0589);
    const ccomplex_t IT_0643 = IT_0641 + IT_0642;
    const ccomplex_t IT_0644 = pow(m_Z, 2);
    const ccomplex_t IT_0645 = s_24*IT_0644;
    const ccomplex_t IT_0646 = IT_0001*IT_0645;
    const ccomplex_t IT_0647 = (-0.25)*IT_0646;
    const ccomplex_t IT_0648 = s_24 + IT_0647;
    const ccomplex_t IT_0649 = m_N_2*IT_0648;
    const ccomplex_t IT_0650 = (-4)*IT_0649;
    const ccomplex_t IT_0651 = 2*IT_0399;
    const ccomplex_t IT_0652 = 2*IT_0469;
    const ccomplex_t IT_0653 = IT_0651 + IT_0652;
    const ccomplex_t IT_0654 = 2*conj(IT_0399);
    const ccomplex_t IT_0655 = 2*conj(IT_0469);
    const ccomplex_t IT_0656 = IT_0654 + IT_0655;
    const ccomplex_t IT_0657 = conj(IT_0082) + conj(IT_0613);
    const ccomplex_t IT_0658 = 2*IT_0516;
    const ccomplex_t IT_0659 = 2*IT_0566;
    const ccomplex_t IT_0660 = IT_0658 + IT_0659;
    const ccomplex_t IT_0661 = IT_0082 + IT_0613;
    const ccomplex_t IT_0662 = 2*conj(IT_0516);
    const ccomplex_t IT_0663 = 2*conj(IT_0566);
    const ccomplex_t IT_0664 = IT_0662 + IT_0663;
    const ccomplex_t IT_0665 = (-2)*IT_0399;
    const ccomplex_t IT_0666 = (-2)*IT_0469;
    const ccomplex_t IT_0667 = IT_0665 + IT_0666;
    const ccomplex_t IT_0668 = (-2)*conj(IT_0399);
    const ccomplex_t IT_0669 = (-2)*conj(IT_0469);
    const ccomplex_t IT_0670 = IT_0668 + IT_0669;
    const ccomplex_t IT_0671 = (-2)*IT_0516;
    const ccomplex_t IT_0672 = (-2)*IT_0566;
    const ccomplex_t IT_0673 = IT_0671 + IT_0672;
    const ccomplex_t IT_0674 = (-2)*conj(IT_0516);
    const ccomplex_t IT_0675 = (-2)*conj(IT_0566);
    const ccomplex_t IT_0676 = IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = s_14*IT_0644;
    const ccomplex_t IT_0678 = IT_0001*IT_0677;
    const ccomplex_t IT_0679 = (-0.25)*IT_0678;
    const ccomplex_t IT_0680 = s_14 + IT_0679;
    const ccomplex_t IT_0681 = m_N_3*IT_0680;
    const ccomplex_t IT_0682 = (-4)*IT_0681;
    const ccomplex_t IT_0683 = 0.5*IT_0678;
    const ccomplex_t IT_0684 = s_14 + IT_0683;
    const ccomplex_t IT_0685 = m_N_3*IT_0684;
    const ccomplex_t IT_0686 = 2*IT_0685;
    const ccomplex_t IT_0687 = 0.5*IT_0646;
    const ccomplex_t IT_0688 = s_24 + IT_0687;
    const ccomplex_t IT_0689 = m_N_2*IT_0688;
    const ccomplex_t IT_0690 = 2*IT_0689;
    const ccomplex_t IT_0691 = s_12*IT_0644;
    const ccomplex_t IT_0692 = s_14*s_24;
    const ccomplex_t IT_0693 = (-2)*IT_0692;
    const ccomplex_t IT_0694 = IT_0691 + IT_0693;
    const ccomplex_t IT_0695 = IT_0001*IT_0694;
    const ccomplex_t IT_0696 = (-0.5)*IT_0695;
    const ccomplex_t IT_0697 = s_12 + IT_0696;
    const ccomplex_t IT_0698 = 4*IT_0516;
    const ccomplex_t IT_0699 = 6*conj(IT_0399);
    const ccomplex_t IT_0700 = pow(m_Z, 4);
    const ccomplex_t IT_0701 = s_12*IT_0700;
    const ccomplex_t IT_0702 = IT_0001*IT_0701;
    const ccomplex_t IT_0703 = (-0.25)*IT_0702;
    const ccomplex_t IT_0704 = IT_0692 + IT_0703;
    const ccomplex_t IT_0705 = -conj(IT_0601);
    const ccomplex_t IT_0706 = IT_0605 + IT_0705;
    const ccomplex_t IT_0707 = -IT_0602;
    const ccomplex_t IT_0708 = IT_0600 + IT_0707;
    const ccomplex_t IT_0709 = -IT_0601 + IT_0607;
    const ccomplex_t IT_0710 = -conj(IT_0602);
    const ccomplex_t IT_0711 = conj(IT_0600) + IT_0710;
    const ccomplex_t IT_0712 = -IT_0610;
    const ccomplex_t IT_0713 = IT_0617 + IT_0712;
    const ccomplex_t IT_0714 = 8*IT_0614;
    const ccomplex_t IT_0715 = 8*IT_0615;
    const ccomplex_t IT_0716 = 8*conj(IT_0614);
    const ccomplex_t IT_0717 = 8*conj(IT_0615);
    const ccomplex_t IT_0718 = -IT_0609*(conj(IT_0610) + -IT_0616) + IT_0706
      *IT_0708 + IT_0709*IT_0711 + conj(IT_0609)*IT_0713 + IT_0657*(IT_0610 + (
      -0.125)*IT_0714 + (-0.125)*IT_0715) + IT_0661*(conj(IT_0610) + (-0.125)
      *IT_0716 + (-0.125)*IT_0717);
    const ccomplex_t IT_0719 = s_12*IT_0006;
    const ccomplex_t IT_0720 = s_12*IT_0635;
    const ccomplex_t IT_0721 = 2*IT_0589;
    const ccomplex_t IT_0722 = s_12*s_24;
    const ccomplex_t IT_0723 = s_14*IT_0009;
    const ccomplex_t IT_0724 = -IT_0723;
    const ccomplex_t IT_0725 = IT_0722 + IT_0724;
    const ccomplex_t IT_0726 = (-2)*IT_0586;
    const ccomplex_t IT_0727 = 2*conj(IT_0589);
    const ccomplex_t IT_0728 = (-2)*conj(IT_0586);
    const ccomplex_t IT_0729 = IT_0606*IT_0622 + IT_0608*IT_0625 + (conj
      (IT_0600) + conj(IT_0601))*IT_0628 + (IT_0600 + IT_0601)*IT_0631 + IT_0611
      *IT_0640 + IT_0612*IT_0643 + (conj(IT_0082) + conj(IT_0613) + IT_0616)*
      (IT_0721 + IT_0726) + (IT_0082 + IT_0613 + IT_0617)*(IT_0727 + IT_0728);
    const ccomplex_t IT_0730 = m_N_3*IT_0006;
    const ccomplex_t IT_0731 = m_N_3*IT_0324;
    const ccomplex_t IT_0732 = IT_0000*IT_0644;
    const ccomplex_t IT_0733 = (-3)*IT_0732;
    const ccomplex_t IT_0734 = 2*IT_0603;
    const ccomplex_t IT_0735 = 2*IT_0604;
    const ccomplex_t IT_0736 = (-2)*IT_0601;
    const ccomplex_t IT_0737 = IT_0734 + IT_0735 + IT_0736;
    const ccomplex_t IT_0738 = 2*conj(IT_0603);
    const ccomplex_t IT_0739 = 2*conj(IT_0604);
    const ccomplex_t IT_0740 = (-2)*conj(IT_0601);
    const ccomplex_t IT_0741 = IT_0738 + IT_0739 + IT_0740;
    const ccomplex_t IT_0742 = 2*IT_0600;
    const ccomplex_t IT_0743 = (-2)*IT_0602;
    const ccomplex_t IT_0744 = IT_0742 + IT_0743;
    const ccomplex_t IT_0745 = 2*conj(IT_0600);
    const ccomplex_t IT_0746 = (-2)*conj(IT_0602);
    const ccomplex_t IT_0747 = IT_0745 + IT_0746;
    const ccomplex_t IT_0748 = 2*IT_0614;
    const ccomplex_t IT_0749 = 2*IT_0615;
    const ccomplex_t IT_0750 = IT_0748 + IT_0749;
    const ccomplex_t IT_0751 = 2*conj(IT_0614);
    const ccomplex_t IT_0752 = 2*conj(IT_0615);
    const ccomplex_t IT_0753 = IT_0751 + IT_0752;
    const ccomplex_t IT_0754 = (-2)*IT_0603;
    const ccomplex_t IT_0755 = (-2)*IT_0604;
    const ccomplex_t IT_0756 = IT_0754 + IT_0755;
    const ccomplex_t IT_0757 = (-2)*conj(IT_0603);
    const ccomplex_t IT_0758 = (-2)*conj(IT_0604);
    const ccomplex_t IT_0759 = IT_0757 + IT_0758;
    const ccomplex_t IT_0760 = 2*IT_0609;
    const ccomplex_t IT_0761 = 2*conj(IT_0609);
    const ccomplex_t IT_0762 = (-2)*IT_0600;
    const ccomplex_t IT_0763 = (-2)*conj(IT_0600);
    const ccomplex_t IT_0764 = conj(IT_0610)*IT_0737 + IT_0610*IT_0741 +
       IT_0657*IT_0744 + IT_0661*IT_0747 + conj(IT_0601)*IT_0750 + IT_0601
      *IT_0753 + IT_0616*IT_0756 + IT_0617*IT_0759 + conj(IT_0602)*IT_0760 +
       IT_0602*IT_0761 + conj(IT_0609)*IT_0762 + IT_0609*IT_0763;
    const ccomplex_t IT_0765 = 3*IT_0732;
    const ccomplex_t IT_0766 = IT_0657*IT_0737 + IT_0661*IT_0741 + conj
      (IT_0610)*IT_0744 + IT_0610*IT_0747 + conj(IT_0602)*IT_0750 + IT_0602
      *IT_0753 + conj(IT_0609)*IT_0756 + IT_0609*IT_0759 + conj(IT_0601)*IT_0760
       + IT_0601*IT_0761 + IT_0616*IT_0762 + IT_0617*IT_0763;
    const ccomplex_t IT_0767 = (-0.5)*IT_0702;
    const ccomplex_t IT_0768 = IT_0691 + IT_0693 + IT_0767;
    const ccomplex_t IT_0769 = 4*conj(IT_0609);
    const ccomplex_t IT_0770 = 4*IT_0609;
    const ccomplex_t IT_0771 = (-4)*IT_0082;
    const ccomplex_t IT_0772 = (-4)*IT_0614;
    const ccomplex_t IT_0773 = IT_0609*conj(IT_0609) + IT_0706*IT_0709 +
       IT_0708*IT_0711 + -conj(IT_0610)*IT_0713 + (-0.25)*IT_0661*IT_0769 +
       IT_0657*(IT_0613 + (-0.25)*IT_0770 + (-0.25)*IT_0771) + IT_0616*(IT_0615 
      + IT_0712 + (-0.25)*IT_0772);
    const ccomplex_t IT_0774 = IT_0007*(conj(IT_0083)*IT_0232 + IT_0083
      *IT_0235 + conj(IT_0315)*IT_0318 + IT_0315*IT_0321) + IT_0325*(IT_0232
      *IT_0470 + IT_0235*IT_0471 + IT_0567*IT_0570 + IT_0571*IT_0574) + IT_0575*
      (IT_0235*IT_0579 + IT_0232*conj(IT_0579) + IT_0321*IT_0582 + IT_0318*conj
      (IT_0582) + IT_0574*IT_0586 + IT_0570*conj(IT_0586) + conj(IT_0589)
      *IT_0592 + IT_0589*IT_0595) + IT_0599*IT_0618 + IT_0619*(IT_0321*IT_0586 +
       IT_0318*conj(IT_0586) + IT_0235*IT_0589 + IT_0232*conj(IT_0589) + conj
      (IT_0315)*IT_0622 + IT_0315*IT_0625 + conj(IT_0083)*IT_0628 + IT_0083
      *IT_0631) + (conj(IT_0586)*IT_0622 + IT_0586*IT_0625 + conj(IT_0589)
      *IT_0628 + IT_0589*IT_0631)*IT_0636 + IT_0637*(IT_0470*IT_0628 + IT_0471
      *IT_0631 + IT_0567*IT_0640 + IT_0571*IT_0643) + IT_0650*(conj(IT_0602)
      *IT_0653 + IT_0602*IT_0656 + IT_0657*IT_0660 + IT_0661*IT_0664 + conj
      (IT_0600)*IT_0667 + IT_0600*IT_0670 + conj(IT_0609)*IT_0673 + IT_0609
      *IT_0676) + (IT_0616*IT_0653 + IT_0617*IT_0656 + IT_0605*IT_0660 + IT_0607
      *IT_0664 + conj(IT_0610)*IT_0667 + IT_0610*IT_0670 + conj(IT_0601)*IT_0673
       + IT_0601*IT_0676)*IT_0682 + (conj(IT_0609)*IT_0653 + IT_0609*IT_0656 +
       conj(IT_0600)*IT_0660 + IT_0600*IT_0664 + IT_0657*IT_0667 + IT_0661
      *IT_0670 + conj(IT_0602)*IT_0673 + IT_0602*IT_0676)*IT_0686 + (conj
      (IT_0601)*IT_0653 + IT_0601*IT_0656 + conj(IT_0610)*IT_0660 + IT_0610
      *IT_0664 + IT_0605*IT_0667 + IT_0607*IT_0670 + IT_0616*IT_0673 + IT_0617
      *IT_0676)*IT_0690 + 4*IT_0697*(IT_0470*IT_0471 + IT_0567*(IT_0566 + 0.25
      *IT_0698)) + 6*IT_0000*(IT_0471*IT_0567 + IT_0571*(conj(IT_0469) +
       0.166666666666667*IT_0699)) + (-8)*IT_0704*IT_0718 + (conj(IT_0197)
      *IT_0198 + conj(IT_0230)*IT_0318 + IT_0230*IT_0320 + conj(IT_0315)*IT_0568
       + conj(IT_0083)*IT_0592 + IT_0083*IT_0594)*IT_0719 + IT_0720*(conj
      (IT_0582)*IT_0620 + conj(IT_0579)*IT_0628 + IT_0579*IT_0630 + conj(IT_0586
      )*IT_0640 + IT_0586*IT_0642 + conj(IT_0589)*IT_0721) + IT_0725*IT_0729 + 
      (IT_0232*IT_0567 + conj(IT_0469)*IT_0570 + IT_0235*IT_0571 + IT_0469
      *IT_0574 + conj(IT_0315)*IT_0651 + IT_0315*IT_0654 + conj(IT_0083)*IT_0665
       + IT_0083*IT_0668)*IT_0730 + (IT_0567*IT_0628 + IT_0571*IT_0631 + conj
      (IT_0469)*IT_0640 + IT_0469*IT_0643 + conj(IT_0586)*IT_0651 + IT_0586
      *IT_0654 + conj(IT_0589)*IT_0665 + IT_0589*IT_0668)*IT_0731 + IT_0733
      *IT_0764 + IT_0765*IT_0766 + (-4)*IT_0768*IT_0773;
    return create_ccomplex_return(IT_0774);
}

