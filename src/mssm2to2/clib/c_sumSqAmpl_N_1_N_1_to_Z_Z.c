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
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = cos(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -4);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = IT_0002*IT_0003*IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = sin(alpha);
    const ccomplex_t IT_0010 = IT_0005*IT_0006*IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0006, -3);
    const ccomplex_t IT_0012 = IT_0002*IT_0003*IT_0011;
    const ccomplex_t IT_0013 = IT_0008*IT_0009*IT_0011;
    const ccomplex_t IT_0014 = cpow(IT_0004, -2);
    const ccomplex_t IT_0015 = cpow(IT_0006, -1);
    const ccomplex_t IT_0016 = IT_0002*IT_0003*IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0008*IT_0009*IT_0014*IT_0015;
    const ccomplex_t IT_0018 = cpow(IT_0006, -2);
    const ccomplex_t IT_0019 = cpow(IT_0014 + IT_0018, -1);
    const ccomplex_t IT_0020 = (IT_0007 + IT_0010 + IT_0012 + IT_0013 + 2
      *IT_0016 + 2*IT_0017)*IT_0019;
    const ccomplex_t IT_0021 = m_W*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0027 = cpow(IT_0004, -1);
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0031 = IT_0009*IT_0027;
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0034 = IT_0003*IT_0015;
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0037 = IT_0009*IT_0015;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0029 + -IT_0032 + 
      -IT_0035 + IT_0038);
    const ccomplex_t IT_0040 = IT_0025*IT_0039;
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = pow(m_Z, 2);
    const ccomplex_t IT_0043 = cpow((-2)*s_13 + IT_0042 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0006*IT_0027;
    const ccomplex_t IT_0046 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = IT_0004*IT_0015;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0051 = IT_0045*IT_0050;
    const ccomplex_t IT_0052 = IT_0048*IT_0050;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0047 + IT_0049 + 
      -IT_0051 + -IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0053;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = m_N_1*IT_0056;
    const ccomplex_t IT_0058 = cpow(IT_0055, 2);
    const ccomplex_t IT_0059 = m_N_1*IT_0058;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = IT_0057 + IT_0060;
    const ccomplex_t IT_0062 = IT_0044*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_13 + IT_0000 + IT_0042 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0066 = IT_0045*IT_0065;
    const ccomplex_t IT_0067 = IT_0048*IT_0065;
    const ccomplex_t IT_0068 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0069 = IT_0045*IT_0068;
    const ccomplex_t IT_0070 = IT_0048*IT_0068;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + IT_0067 + 
      -IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = cpow(IT_0072, 2);
    const ccomplex_t IT_0074 = m_N_2*IT_0073;
    const ccomplex_t IT_0075 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0076 = IT_0045*IT_0075;
    const ccomplex_t IT_0077 = IT_0048*IT_0075;
    const ccomplex_t IT_0078 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0079 = IT_0045*IT_0078;
    const ccomplex_t IT_0080 = IT_0048*IT_0078;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0076 + IT_0077 + 
      -IT_0079 + -IT_0080);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = IT_0072*IT_0082;
    const ccomplex_t IT_0084 = m_N_1*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0074 + IT_0085;
    const ccomplex_t IT_0087 = IT_0064*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0000 + IT_0042 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0091 = IT_0045*IT_0090;
    const ccomplex_t IT_0092 = IT_0048*IT_0090;
    const ccomplex_t IT_0093 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0094 = IT_0045*IT_0093;
    const ccomplex_t IT_0095 = IT_0048*IT_0093;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0091 + IT_0092 + 
      -IT_0094 + -IT_0095);
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = cpow(IT_0097, 2);
    const ccomplex_t IT_0099 = m_N_3*IT_0098;
    const ccomplex_t IT_0100 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0101 = IT_0045*IT_0100;
    const ccomplex_t IT_0102 = IT_0048*IT_0100;
    const ccomplex_t IT_0103 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0104 = IT_0045*IT_0103;
    const ccomplex_t IT_0105 = IT_0048*IT_0103;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(IT_0101 + IT_0102 + 
      -IT_0104 + -IT_0105);
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = IT_0097*IT_0107;
    const ccomplex_t IT_0109 = m_N_1*IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0099 + IT_0110;
    const ccomplex_t IT_0112 = IT_0089*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_13 + IT_0000 + IT_0042 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0116 = IT_0045*IT_0115;
    const ccomplex_t IT_0117 = IT_0048*IT_0115;
    const ccomplex_t IT_0118 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0119 = IT_0045*IT_0118;
    const ccomplex_t IT_0120 = IT_0048*IT_0118;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0116 + IT_0117 + 
      -IT_0119 + -IT_0120);
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = cpow(IT_0122, 2);
    const ccomplex_t IT_0124 = m_N_4*IT_0123;
    const ccomplex_t IT_0125 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0126 = IT_0045*IT_0125;
    const ccomplex_t IT_0127 = IT_0048*IT_0125;
    const ccomplex_t IT_0128 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0129 = IT_0045*IT_0128;
    const ccomplex_t IT_0130 = IT_0048*IT_0128;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0126 + IT_0127 + 
      -IT_0129 + -IT_0130);
    const ccomplex_t IT_0132 = (-0.5)*IT_0131;
    const ccomplex_t IT_0133 = IT_0122*IT_0132;
    const ccomplex_t IT_0134 = m_N_1*IT_0133;
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = IT_0124 + IT_0135;
    const ccomplex_t IT_0137 = IT_0114*IT_0136;
    const ccomplex_t IT_0138 = -IT_0062 + IT_0087 + IT_0112 + IT_0137;
    const ccomplex_t IT_0139 = conj(IT_0041) + conj(IT_0138);
    const ccomplex_t IT_0140 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0141 = IT_0028*IT_0140;
    const ccomplex_t IT_0142 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0143 = IT_0031*IT_0142;
    const ccomplex_t IT_0144 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0145 = IT_0034*IT_0144;
    const ccomplex_t IT_0146 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0147 = IT_0037*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0141 + -IT_0143 + 
      -IT_0145 + IT_0147);
    const ccomplex_t IT_0149 = IT_0025*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = pow(m_Z, -4);
    const ccomplex_t IT_0152 = pow(s_34, 2);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = 1 + IT_0154;
    const ccomplex_t IT_0156 = IT_0000*IT_0155;
    const ccomplex_t IT_0157 = (-4)*IT_0156;
    const ccomplex_t IT_0158 = IT_0150*IT_0157;
    const ccomplex_t IT_0159 = cpow((-2)*s_23 + IT_0042 + (0 + _Complex_I*1)
      *m_N_1*Gamma_o1 + reg_prop, -1);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0056*IT_0160;
    const ccomplex_t IT_0162 = IT_0058*IT_0160;
    const ccomplex_t IT_0163 = cpow((-2)*s_23 + IT_0000 + IT_0042 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0083*IT_0164;
    const ccomplex_t IT_0166 = cpow((-2)*s_23 + IT_0000 + IT_0042 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = IT_0133*IT_0167;
    const ccomplex_t IT_0169 = cpow((-2)*s_23 + IT_0000 + IT_0042 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0108*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0161 + 0.5*IT_0162 + -IT_0165 + 
      -IT_0168 + -IT_0171;
    const ccomplex_t IT_0173 = IT_0044*IT_0056;
    const ccomplex_t IT_0174 = cpow(IT_0054, 2);
    const ccomplex_t IT_0175 = IT_0044*IT_0174;
    const ccomplex_t IT_0176 = IT_0114*IT_0133;
    const ccomplex_t IT_0177 = IT_0089*IT_0108;
    const ccomplex_t IT_0178 = IT_0064*IT_0083;
    const ccomplex_t IT_0179 = 0.5*IT_0173 + (-0.5)*IT_0175 + IT_0176 +
       IT_0177 + IT_0178;
    const ccomplex_t IT_0180 = conj(IT_0172) + conj(IT_0179);
    const ccomplex_t IT_0181 = s_34*IT_0151;
    const ccomplex_t IT_0182 = s_23*IT_0042;
    const ccomplex_t IT_0183 = IT_0181*IT_0182;
    const ccomplex_t IT_0184 = pow(m_Z, -2);
    const ccomplex_t IT_0185 = s_24*IT_0042;
    const ccomplex_t IT_0186 = IT_0184*IT_0185;
    const ccomplex_t IT_0187 = s_23*s_34;
    const ccomplex_t IT_0188 = (-2)*IT_0187;
    const ccomplex_t IT_0189 = IT_0185 + IT_0188;
    const ccomplex_t IT_0190 = IT_0184*IT_0189;
    const ccomplex_t IT_0191 = (-0.5)*IT_0183 + 0.5*IT_0186 + (-0.5)*IT_0190;
    const ccomplex_t IT_0192 = s_24 + IT_0191;
    const ccomplex_t IT_0193 = m_N_1*IT_0192;
    const ccomplex_t IT_0194 = (-2)*IT_0193;
    const ccomplex_t IT_0195 = 2*IT_0194;
    const ccomplex_t IT_0196 = IT_0041*IT_0195;
    const ccomplex_t IT_0197 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0198 = IT_0003*IT_0005*IT_0006*IT_0008;
    const ccomplex_t IT_0199 = IT_0002*IT_0005*IT_0006*IT_0009;
    const ccomplex_t IT_0200 = IT_0003*IT_0008*IT_0011;
    const ccomplex_t IT_0201 = IT_0002*IT_0009*IT_0011;
    const ccomplex_t IT_0202 = IT_0003*IT_0008*IT_0014*IT_0015;
    const ccomplex_t IT_0203 = IT_0002*IT_0009*IT_0014*IT_0015;
    const ccomplex_t IT_0204 = IT_0019*(IT_0198 + -IT_0199 + IT_0200 + 
      -IT_0201 + 2*IT_0202 + (-2)*IT_0203);
    const ccomplex_t IT_0205 = IT_0021*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = IT_0197*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = IT_0028*IT_0030;
    const ccomplex_t IT_0210 = IT_0026*IT_0031;
    const ccomplex_t IT_0211 = IT_0034*IT_0036;
    const ccomplex_t IT_0212 = IT_0033*IT_0037;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*(IT_0209 + IT_0210 + 
      -IT_0211 + -IT_0212);
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0208*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0040 + -IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0194;
    const ccomplex_t IT_0218 = IT_0216*IT_0217;
    const ccomplex_t IT_0219 = IT_0160*IT_0174;
    const ccomplex_t IT_0220 = (-0.5)*IT_0161 + -IT_0165 + -IT_0168 + -IT_0171
       + 0.5*IT_0219;
    const ccomplex_t IT_0221 = IT_0044*IT_0058;
    const ccomplex_t IT_0222 = 0.5*IT_0173 + IT_0176 + IT_0177 + IT_0178 + (
      -0.5)*IT_0221;
    const ccomplex_t IT_0223 = conj(IT_0220) + conj(IT_0222);
    const ccomplex_t IT_0224 = s_13*IT_0042;
    const ccomplex_t IT_0225 = IT_0181*IT_0224;
    const ccomplex_t IT_0226 = s_14*IT_0042;
    const ccomplex_t IT_0227 = IT_0184*IT_0226;
    const ccomplex_t IT_0228 = s_13*s_34;
    const ccomplex_t IT_0229 = (-2)*IT_0228;
    const ccomplex_t IT_0230 = IT_0226 + IT_0229;
    const ccomplex_t IT_0231 = IT_0184*IT_0230;
    const ccomplex_t IT_0232 = (-0.5)*IT_0225 + 0.5*IT_0227 + (-0.5)*IT_0231;
    const ccomplex_t IT_0233 = s_14 + IT_0232;
    const ccomplex_t IT_0234 = m_N_1*IT_0233;
    const ccomplex_t IT_0235 = (-2)*IT_0234;
    const ccomplex_t IT_0236 = 2*IT_0235;
    const ccomplex_t IT_0237 = IT_0216*IT_0236;
    const ccomplex_t IT_0238 = (-2)*IT_0235;
    const ccomplex_t IT_0239 = IT_0041*IT_0238;
    const ccomplex_t IT_0240 = conj(IT_0041)*IT_0195;
    const ccomplex_t IT_0241 = IT_0028*IT_0142;
    const ccomplex_t IT_0242 = IT_0031*IT_0140;
    const ccomplex_t IT_0243 = IT_0034*IT_0146;
    const ccomplex_t IT_0244 = IT_0037*IT_0144;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*(IT_0241 + IT_0242 + 
      -IT_0243 + -IT_0244);
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = IT_0208*IT_0246;
    const ccomplex_t IT_0248 = (-0.5)*IT_0149 + -IT_0247;
    const ccomplex_t IT_0249 = IT_0236*conj(IT_0248);
    const ccomplex_t IT_0250 = cpow(IT_0082, 2);
    const ccomplex_t IT_0251 = m_N_2*IT_0250;
    const ccomplex_t IT_0252 = IT_0085 + IT_0251;
    const ccomplex_t IT_0253 = IT_0164*IT_0252;
    const ccomplex_t IT_0254 = cpow(IT_0107, 2);
    const ccomplex_t IT_0255 = m_N_3*IT_0254;
    const ccomplex_t IT_0256 = IT_0110 + IT_0255;
    const ccomplex_t IT_0257 = IT_0170*IT_0256;
    const ccomplex_t IT_0258 = cpow(IT_0132, 2);
    const ccomplex_t IT_0259 = m_N_4*IT_0258;
    const ccomplex_t IT_0260 = IT_0135 + IT_0259;
    const ccomplex_t IT_0261 = IT_0167*IT_0260;
    const ccomplex_t IT_0262 = m_N_1*IT_0174;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = IT_0057 + IT_0263;
    const ccomplex_t IT_0265 = IT_0160*IT_0264;
    const ccomplex_t IT_0266 = -IT_0253 + -IT_0257 + -IT_0261 + IT_0265;
    const ccomplex_t IT_0267 = m_N_1*IT_0184;
    const ccomplex_t IT_0268 = pow(m_Z, 4);
    const ccomplex_t IT_0269 = s_14*IT_0268;
    const ccomplex_t IT_0270 = IT_0184*IT_0269;
    const ccomplex_t IT_0271 = (-2)*IT_0228 + (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = IT_0226 + IT_0271;
    const ccomplex_t IT_0273 = IT_0267*IT_0272;
    const ccomplex_t IT_0274 = (-2)*IT_0273;
    const ccomplex_t IT_0275 = 2*IT_0274;
    const ccomplex_t IT_0276 = IT_0061*IT_0160;
    const ccomplex_t IT_0277 = IT_0086*IT_0164;
    const ccomplex_t IT_0278 = IT_0111*IT_0170;
    const ccomplex_t IT_0279 = IT_0136*IT_0167;
    const ccomplex_t IT_0280 = IT_0276 + -IT_0277 + -IT_0278 + -IT_0279;
    const ccomplex_t IT_0281 = s_23*s_34*IT_0184;
    const ccomplex_t IT_0282 = s_24*IT_0268;
    const ccomplex_t IT_0283 = s_23*s_34*IT_0042;
    const ccomplex_t IT_0284 = (-2)*IT_0283;
    const ccomplex_t IT_0285 = IT_0282 + IT_0284;
    const ccomplex_t IT_0286 = IT_0151*IT_0285;
    const ccomplex_t IT_0287 = (-0.666666666666667)*IT_0281 + (
      -0.166666666666667)*IT_0286;
    const ccomplex_t IT_0288 = s_24 + IT_0287;
    const ccomplex_t IT_0289 = m_N_1*IT_0288;
    const ccomplex_t IT_0290 = (-6)*IT_0289;
    const ccomplex_t IT_0291 = 2*IT_0290;
    const ccomplex_t IT_0292 = IT_0044*IT_0264;
    const ccomplex_t IT_0293 = IT_0064*IT_0252;
    const ccomplex_t IT_0294 = IT_0089*IT_0256;
    const ccomplex_t IT_0295 = IT_0114*IT_0260;
    const ccomplex_t IT_0296 = -IT_0292 + IT_0293 + IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_13*s_34*IT_0042;
    const ccomplex_t IT_0298 = (-2)*IT_0297;
    const ccomplex_t IT_0299 = IT_0269 + IT_0298;
    const ccomplex_t IT_0300 = IT_0151*IT_0299;
    const ccomplex_t IT_0301 = (-0.666666666666667)*IT_0231 +
       0.166666666666667*IT_0300;
    const ccomplex_t IT_0302 = s_14 + IT_0301;
    const ccomplex_t IT_0303 = m_N_1*IT_0302;
    const ccomplex_t IT_0304 = 6*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0304;
    const ccomplex_t IT_0306 = IT_0151*IT_0282;
    const ccomplex_t IT_0307 = (-10)*s_24;
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = m_N_1*IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0309;
    const ccomplex_t IT_0311 = -IT_0173 + IT_0175 + (-2)*IT_0176 + (-2)
      *IT_0177 + (-2)*IT_0178;
    const ccomplex_t IT_0312 = s_12*s_24;
    const ccomplex_t IT_0313 = s_12*s_24*IT_0042;
    const ccomplex_t IT_0314 = s_14*IT_0000*IT_0042;
    const ccomplex_t IT_0315 = s_13*s_23*s_24;
    const ccomplex_t IT_0316 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0317 = -IT_0314 + (-2)*IT_0315 + 2*IT_0316;
    const ccomplex_t IT_0318 = IT_0313 + IT_0317;
    const ccomplex_t IT_0319 = IT_0184*IT_0318;
    const ccomplex_t IT_0320 = (-0.5)*IT_0319;
    const ccomplex_t IT_0321 = s_14*IT_0000;
    const ccomplex_t IT_0322 = -IT_0321;
    const ccomplex_t IT_0323 = IT_0312 + IT_0320 + IT_0322;
    const ccomplex_t IT_0324 = 4*IT_0323;
    const ccomplex_t IT_0325 = s_12*IT_0042;
    const ccomplex_t IT_0326 = pow(m_Z, 6);
    const ccomplex_t IT_0327 = s_12*IT_0326;
    const ccomplex_t IT_0328 = s_13*s_23*IT_0268;
    const ccomplex_t IT_0329 = (-2)*IT_0328;
    const ccomplex_t IT_0330 = IT_0327 + IT_0329;
    const ccomplex_t IT_0331 = IT_0151*IT_0330;
    const ccomplex_t IT_0332 = 0.25*IT_0331;
    const ccomplex_t IT_0333 = IT_0325 + IT_0332;
    const ccomplex_t IT_0334 = (-8)*IT_0333;
    const ccomplex_t IT_0335 = IT_0161 + -IT_0162 + 2*IT_0165 + 2*IT_0168 + 2
      *IT_0171;
    const ccomplex_t IT_0336 = s_12*s_14;
    const ccomplex_t IT_0337 = s_12*s_14*IT_0042;
    const ccomplex_t IT_0338 = s_24*IT_0000*IT_0042;
    const ccomplex_t IT_0339 = pow(s_13, 2);
    const ccomplex_t IT_0340 = s_24*IT_0339;
    const ccomplex_t IT_0341 = s_12*s_13*s_34;
    const ccomplex_t IT_0342 = -IT_0338 + 2*IT_0340 + (-2)*IT_0341;
    const ccomplex_t IT_0343 = IT_0337 + IT_0342;
    const ccomplex_t IT_0344 = IT_0184*IT_0343;
    const ccomplex_t IT_0345 = (-0.5)*IT_0344;
    const ccomplex_t IT_0346 = s_24*IT_0000;
    const ccomplex_t IT_0347 = -IT_0346;
    const ccomplex_t IT_0348 = IT_0336 + IT_0345 + IT_0347;
    const ccomplex_t IT_0349 = (-4)*IT_0348;
    const ccomplex_t IT_0350 = s_12*IT_0268;
    const ccomplex_t IT_0351 = s_13*s_23*IT_0042;
    const ccomplex_t IT_0352 = s_14*s_24*IT_0042;
    const ccomplex_t IT_0353 = s_13*s_24*s_34;
    const ccomplex_t IT_0354 = (-2)*IT_0351 + (-2)*IT_0352 + 4*IT_0353;
    const ccomplex_t IT_0355 = IT_0350 + IT_0354;
    const ccomplex_t IT_0356 = IT_0184*IT_0355;
    const ccomplex_t IT_0357 = 0.5*IT_0331;
    const ccomplex_t IT_0358 = s_14*s_24;
    const ccomplex_t IT_0359 = (-4)*IT_0358;
    const ccomplex_t IT_0360 = (-4)*IT_0325 + 4*IT_0356 + (-4)*IT_0357 + (-4)
      *IT_0359;
    const ccomplex_t IT_0361 = IT_0000*IT_0042;
    const ccomplex_t IT_0362 = 9*IT_0361;
    const ccomplex_t IT_0363 = (-2)*IT_0362;
    const ccomplex_t IT_0364 = conj(IT_0216)*IT_0217;
    const ccomplex_t IT_0365 = conj(IT_0150)*IT_0238;
    const ccomplex_t IT_0366 = IT_0000*IT_0184;
    const ccomplex_t IT_0367 = (-0.25)*IT_0268;
    const ccomplex_t IT_0368 = IT_0152 + IT_0367;
    const ccomplex_t IT_0369 = IT_0366*IT_0368;
    const ccomplex_t IT_0370 = (-4)*IT_0369;
    const ccomplex_t IT_0371 = (-2)*IT_0370;
    const ccomplex_t IT_0372 = IT_0240 + IT_0249 + conj(IT_0266)*IT_0275 +
       conj(IT_0280)*IT_0291 + conj(IT_0296)*IT_0305 + conj(IT_0138)*IT_0310 +
       conj(IT_0311)*IT_0324 + conj(IT_0172)*IT_0334 + conj(IT_0335)*IT_0349 +
       conj(IT_0179)*IT_0360 + conj(IT_0222)*IT_0363 + IT_0364 + IT_0365 + conj
      (IT_0220)*IT_0371;
    const ccomplex_t IT_0373 = conj(IT_0150)*IT_0195;
    const ccomplex_t IT_0374 = conj(IT_0216)*IT_0236;
    const ccomplex_t IT_0375 = -IT_0173 + (-2)*IT_0176 + (-2)*IT_0177 + (-2)
      *IT_0178 + IT_0221;
    const ccomplex_t IT_0376 = IT_0161 + 2*IT_0165 + 2*IT_0168 + 2*IT_0171 + 
      -IT_0219;
    const ccomplex_t IT_0377 = IT_0217*conj(IT_0248);
    const ccomplex_t IT_0378 = conj(IT_0041)*IT_0238;
    const ccomplex_t IT_0379 = IT_0275*conj(IT_0280) + conj(IT_0266)*IT_0291 +
       conj(IT_0138)*IT_0305 + conj(IT_0296)*IT_0310 + conj(IT_0220)*IT_0334 +
       conj(IT_0222)*IT_0360 + conj(IT_0179)*IT_0363 + conj(IT_0172)*IT_0371 +
       IT_0373 + IT_0374 + IT_0324*conj(IT_0375) + IT_0349*conj(IT_0376) +
       IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = s_12*IT_0155;
    const ccomplex_t IT_0381 = 4*IT_0156;
    const ccomplex_t IT_0382 = IT_0157*conj(IT_0248);
    const ccomplex_t IT_0383 = (-4)*IT_0380;
    const ccomplex_t IT_0384 = 0.25*conj(IT_0041);
    const ccomplex_t IT_0385 = s_14*s_23;
    const ccomplex_t IT_0386 = s_13*s_24;
    const ccomplex_t IT_0387 = s_12*s_34;
    const ccomplex_t IT_0388 = -IT_0386 + -IT_0387;
    const ccomplex_t IT_0389 = IT_0385 + IT_0388;
    const ccomplex_t IT_0390 = IT_0181*IT_0389;
    const ccomplex_t IT_0391 = (-0.25)*IT_0390;
    const ccomplex_t IT_0392 = IT_0184*IT_0325;
    const ccomplex_t IT_0393 = (-0.5)*IT_0392;
    const ccomplex_t IT_0394 = s_12 + IT_0391 + IT_0393;
    const ccomplex_t IT_0395 = (-8)*IT_0394;
    const ccomplex_t IT_0396 = 0.25*conj(IT_0138);
    const ccomplex_t IT_0397 = s_13*s_14*s_34*IT_0151;
    const ccomplex_t IT_0398 = IT_0000 + IT_0397;
    const ccomplex_t IT_0399 = IT_0184*IT_0339;
    const ccomplex_t IT_0400 = pow(s_14, 2);
    const ccomplex_t IT_0401 = IT_0184*IT_0400;
    const ccomplex_t IT_0402 = -IT_0399 + -IT_0401;
    const ccomplex_t IT_0403 = IT_0398 + IT_0402;
    const ccomplex_t IT_0404 = m_N_1*IT_0403;
    const ccomplex_t IT_0405 = -IT_0404;
    const ccomplex_t IT_0406 = (-2)*IT_0405;
    const ccomplex_t IT_0407 = 0.25*conj(IT_0376);
    const ccomplex_t IT_0408 = 0.25*conj(IT_0266);
    const ccomplex_t IT_0409 = s_14*s_23*s_34*IT_0151;
    const ccomplex_t IT_0410 = s_13*s_23*IT_0184;
    const ccomplex_t IT_0411 = s_14*s_24*IT_0184;
    const ccomplex_t IT_0412 = -IT_0410 + -IT_0411;
    const ccomplex_t IT_0413 = s_12 + IT_0409 + IT_0412;
    const ccomplex_t IT_0414 = m_N_1*IT_0413;
    const ccomplex_t IT_0415 = (-2)*IT_0414;
    const ccomplex_t IT_0416 = 0.25*conj(IT_0335);
    const ccomplex_t IT_0417 = IT_0385 + IT_0387;
    const ccomplex_t IT_0418 = -IT_0386;
    const ccomplex_t IT_0419 = IT_0417 + IT_0418;
    const ccomplex_t IT_0420 = IT_0181*IT_0419;
    const ccomplex_t IT_0421 = 0.25*IT_0420;
    const ccomplex_t IT_0422 = s_12 + IT_0393 + IT_0421;
    const ccomplex_t IT_0423 = 8*IT_0422;
    const ccomplex_t IT_0424 = 0.25*conj(IT_0280);
    const ccomplex_t IT_0425 = s_23*s_24*s_34*IT_0151;
    const ccomplex_t IT_0426 = IT_0000 + IT_0425;
    const ccomplex_t IT_0427 = pow(s_23, 2);
    const ccomplex_t IT_0428 = IT_0184*IT_0427;
    const ccomplex_t IT_0429 = pow(s_24, 2);
    const ccomplex_t IT_0430 = IT_0184*IT_0429;
    const ccomplex_t IT_0431 = -IT_0428 + -IT_0430;
    const ccomplex_t IT_0432 = IT_0426 + IT_0431;
    const ccomplex_t IT_0433 = m_N_1*IT_0432;
    const ccomplex_t IT_0434 = -IT_0433;
    const ccomplex_t IT_0435 = 2*IT_0434;
    const ccomplex_t IT_0436 = 0.25*conj(IT_0311);
    const ccomplex_t IT_0437 = s_13*s_24*s_34*IT_0151;
    const ccomplex_t IT_0438 = s_12 + IT_0412 + IT_0437;
    const ccomplex_t IT_0439 = m_N_1*IT_0438;
    const ccomplex_t IT_0440 = 2*IT_0439;
    const ccomplex_t IT_0441 = 0.25*conj(IT_0375);
    const ccomplex_t IT_0442 = conj(IT_0216)*IT_0380 + 0.25*(conj(IT_0150) +
       conj(IT_0296))*IT_0381 + 0.25*IT_0382 + IT_0383*IT_0384 + IT_0395*IT_0396
       + IT_0406*IT_0407 + IT_0157*IT_0408 + IT_0415*IT_0416 + IT_0423*IT_0424 +
       IT_0435*IT_0436 + IT_0440*IT_0441;
    const ccomplex_t IT_0443 = conj(IT_0150)*IT_0157;
    const ccomplex_t IT_0444 = 0.25*conj(IT_0216);
    const ccomplex_t IT_0445 = 8*IT_0394;
    const ccomplex_t IT_0446 = 2*IT_0405;
    const ccomplex_t IT_0447 = 2*IT_0414;
    const ccomplex_t IT_0448 = (-8)*IT_0422;
    const ccomplex_t IT_0449 = (-2)*IT_0434;
    const ccomplex_t IT_0450 = (-2)*IT_0439;
    const ccomplex_t IT_0451 = 0.25*conj(IT_0296);
    const ccomplex_t IT_0452 = conj(IT_0041)*IT_0380 + 0.25*(conj(IT_0248) +
       conj(IT_0266))*IT_0381 + 0.25*IT_0443 + IT_0383*IT_0444 + IT_0396*IT_0445
       + IT_0407*IT_0446 + IT_0416*IT_0447 + IT_0424*IT_0448 + IT_0436*IT_0449 +
       IT_0441*IT_0450 + IT_0157*IT_0451;
    const ccomplex_t IT_0453 = s_23*s_24*s_34*IT_0184;
    const ccomplex_t IT_0454 = (-0.5)*IT_0361 + 0.5*IT_0429 + -IT_0453;
    const ccomplex_t IT_0455 = IT_0427 + IT_0454;
    const ccomplex_t IT_0456 = IT_0184*IT_0455;
    const ccomplex_t IT_0457 = -IT_0430;
    const ccomplex_t IT_0458 = IT_0000 + IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = m_N_1*IT_0458;
    const ccomplex_t IT_0460 = s_14*IT_0427;
    const ccomplex_t IT_0461 = s_12*s_23*s_34;
    const ccomplex_t IT_0462 = 0.5*IT_0313 + (-0.5)*IT_0314 + -IT_0461;
    const ccomplex_t IT_0463 = IT_0460 + IT_0462;
    const ccomplex_t IT_0464 = IT_0184*IT_0463;
    const ccomplex_t IT_0465 = -IT_0312;
    const ccomplex_t IT_0466 = IT_0321 + IT_0464 + IT_0465;
    const ccomplex_t IT_0467 = 0.25*conj(IT_0248);
    const ccomplex_t IT_0468 = 0.25*conj(IT_0150);
    const ccomplex_t IT_0469 = -IT_0361;
    const ccomplex_t IT_0470 = IT_0429 + IT_0469;
    const ccomplex_t IT_0471 = IT_0184*IT_0470;
    const ccomplex_t IT_0472 = -IT_0000;
    const ccomplex_t IT_0473 = IT_0430 + IT_0472;
    const ccomplex_t IT_0474 = (-4)*IT_0473;
    const ccomplex_t IT_0475 = IT_0471 + IT_0474;
    const ccomplex_t IT_0476 = m_N_1*IT_0475;
    const ccomplex_t IT_0477 = 2*IT_0476;
    const ccomplex_t IT_0478 = -IT_0411;
    const ccomplex_t IT_0479 = s_12 + IT_0478;
    const ccomplex_t IT_0480 = IT_0000*IT_0479;
    const ccomplex_t IT_0481 = (-6)*IT_0480;
    const ccomplex_t IT_0482 = -IT_0358;
    const ccomplex_t IT_0483 = IT_0325 + IT_0482;
    const ccomplex_t IT_0484 = IT_0184*IT_0483;
    const ccomplex_t IT_0485 = -IT_0411 + (-0.25)*IT_0484;
    const ccomplex_t IT_0486 = s_12 + IT_0485;
    const ccomplex_t IT_0487 = m_N_1*IT_0486;
    const ccomplex_t IT_0488 = (-4)*IT_0487;
    const ccomplex_t IT_0489 = (-2)*IT_0488;
    const ccomplex_t IT_0490 = s_13*s_23;
    const ccomplex_t IT_0491 = (-2)*IT_0490;
    const ccomplex_t IT_0492 = IT_0325 + IT_0491;
    const ccomplex_t IT_0493 = IT_0184*IT_0492;
    const ccomplex_t IT_0494 = (s_12 + -IT_0411)*(s_12 + (-0.5)*IT_0493);
    const ccomplex_t IT_0495 = (-2)*IT_0494;
    const ccomplex_t IT_0496 = 2*IT_0495;
    const ccomplex_t IT_0497 = (IT_0000 + -IT_0430)*(s_12 + (-0.5)*IT_0493);
    const ccomplex_t IT_0498 = 2*IT_0497;
    const ccomplex_t IT_0499 = (-2)*IT_0498;
    const ccomplex_t IT_0500 = IT_0000 + IT_0457;
    const ccomplex_t IT_0501 = IT_0000*IT_0500;
    const ccomplex_t IT_0502 = -IT_0501;
    const ccomplex_t IT_0503 = 6*IT_0502;
    const ccomplex_t IT_0504 = s_13*s_24*s_34*IT_0184;
    const ccomplex_t IT_0505 = -IT_0358 + (-2)*IT_0490 + 2*IT_0504;
    const ccomplex_t IT_0506 = IT_0325 + IT_0505;
    const ccomplex_t IT_0507 = IT_0184*IT_0506;
    const ccomplex_t IT_0508 = -IT_0411 + (-0.5)*IT_0507;
    const ccomplex_t IT_0509 = s_12 + IT_0508;
    const ccomplex_t IT_0510 = m_N_1*IT_0509;
    const ccomplex_t IT_0511 = (-2)*IT_0510;
    const ccomplex_t IT_0512 = (-2)*IT_0511;
    const ccomplex_t IT_0513 = 0.25*conj(IT_0179);
    const ccomplex_t IT_0514 = IT_0435*IT_0444 + IT_0384*IT_0449 + conj
      (IT_0280)*IT_0459 + conj(IT_0172)*IT_0466 + IT_0440*IT_0467 + IT_0450
      *IT_0468 + IT_0396*IT_0477 + IT_0407*IT_0481 + IT_0408*IT_0489 + IT_0416
      *IT_0496 + IT_0436*IT_0499 + IT_0441*IT_0503 + IT_0451*IT_0512 + IT_0324
      *IT_0513;
    const ccomplex_t IT_0515 = -IT_0401;
    const ccomplex_t IT_0516 = s_13*s_14*s_34*IT_0184;
    const ccomplex_t IT_0517 = (-0.5)*IT_0361 + 0.5*IT_0400 + -IT_0516;
    const ccomplex_t IT_0518 = IT_0339 + IT_0517;
    const ccomplex_t IT_0519 = IT_0184*IT_0518;
    const ccomplex_t IT_0520 = IT_0000 + IT_0515 + IT_0519;
    const ccomplex_t IT_0521 = m_N_1*IT_0520;
    const ccomplex_t IT_0522 = IT_0000 + IT_0515;
    const ccomplex_t IT_0523 = IT_0000*IT_0522;
    const ccomplex_t IT_0524 = -IT_0523;
    const ccomplex_t IT_0525 = 6*IT_0524;
    const ccomplex_t IT_0526 = IT_0400 + IT_0469;
    const ccomplex_t IT_0527 = IT_0184*IT_0526;
    const ccomplex_t IT_0528 = IT_0401 + IT_0472;
    const ccomplex_t IT_0529 = (-4)*IT_0528;
    const ccomplex_t IT_0530 = IT_0527 + IT_0529;
    const ccomplex_t IT_0531 = m_N_1*IT_0530;
    const ccomplex_t IT_0532 = 2*IT_0531;
    const ccomplex_t IT_0533 = (IT_0000 + -IT_0401)*(s_12 + (-0.5)*IT_0493);
    const ccomplex_t IT_0534 = 2*IT_0533;
    const ccomplex_t IT_0535 = (-2)*IT_0534;
    const ccomplex_t IT_0536 = s_14*s_23*s_34*IT_0184;
    const ccomplex_t IT_0537 = -IT_0358 + (-2)*IT_0490 + 2*IT_0536;
    const ccomplex_t IT_0538 = IT_0325 + IT_0537;
    const ccomplex_t IT_0539 = IT_0184*IT_0538;
    const ccomplex_t IT_0540 = -IT_0411 + (-0.5)*IT_0539;
    const ccomplex_t IT_0541 = s_12 + IT_0540;
    const ccomplex_t IT_0542 = m_N_1*IT_0541;
    const ccomplex_t IT_0543 = (-2)*IT_0542;
    const ccomplex_t IT_0544 = (-2)*IT_0543;
    const ccomplex_t IT_0545 = s_13*s_14*s_23;
    const ccomplex_t IT_0546 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0547 = -IT_0338 + (-2)*IT_0545 + 2*IT_0546;
    const ccomplex_t IT_0548 = IT_0337 + IT_0547;
    const ccomplex_t IT_0549 = IT_0184*IT_0548;
    const ccomplex_t IT_0550 = 2*IT_0346;
    const ccomplex_t IT_0551 = (-2)*IT_0336;
    const ccomplex_t IT_0552 = IT_0549 + IT_0550 + IT_0551;
    const ccomplex_t IT_0553 = (-2)*IT_0552;
    const ccomplex_t IT_0554 = 0.25*conj(IT_0172);
    const ccomplex_t IT_0555 = IT_0415*IT_0444 + IT_0384*IT_0447 + IT_0406
      *IT_0467 + IT_0446*IT_0468 + IT_0441*IT_0481 + IT_0396*IT_0489 + IT_0436
      *IT_0496 + IT_0349*IT_0513 + conj(IT_0296)*IT_0521 + IT_0407*IT_0525 +
       IT_0408*IT_0532 + IT_0416*IT_0535 + IT_0424*IT_0544 + IT_0553*IT_0554;
    const ccomplex_t IT_0556 = s_14*s_23*s_34;
    const ccomplex_t IT_0557 = s_12*IT_0152;
    const ccomplex_t IT_0558 = (-2)*IT_0351 + (-2)*IT_0352 + 2*IT_0353 + 2
      *IT_0556 + (-2)*IT_0557;
    const ccomplex_t IT_0559 = IT_0350 + IT_0558;
    const ccomplex_t IT_0560 = IT_0151*IT_0559;
    const ccomplex_t IT_0561 = (-0.25)*IT_0560;
    const ccomplex_t IT_0562 = -IT_0410;
    const ccomplex_t IT_0563 = s_12 + IT_0478 + IT_0561 + IT_0562;
    const ccomplex_t IT_0564 = (-8)*IT_0563;
    const ccomplex_t IT_0565 = 1.125*IT_0000;
    const ccomplex_t IT_0566 = (-16)*IT_0565;
    const ccomplex_t IT_0567 = (-2)*IT_0351 + (-2)*IT_0352 + 4*IT_0556;
    const ccomplex_t IT_0568 = IT_0350 + IT_0567;
    const ccomplex_t IT_0569 = IT_0151*IT_0568;
    const ccomplex_t IT_0570 = (-2)*IT_0358;
    const ccomplex_t IT_0571 = IT_0325 + IT_0570;
    const ccomplex_t IT_0572 = IT_0184*IT_0571;
    const ccomplex_t IT_0573 = (-0.5)*IT_0493;
    const ccomplex_t IT_0574 = s_12 + IT_0573;
    const ccomplex_t IT_0575 = 2*IT_0569 + (-4)*IT_0572 + 8*IT_0574;
    const ccomplex_t IT_0576 = (-2)*IT_0152;
    const ccomplex_t IT_0577 = IT_0268 + IT_0576;
    const ccomplex_t IT_0578 = IT_0151*IT_0577;
    const ccomplex_t IT_0579 = 0.25*IT_0578;
    const ccomplex_t IT_0580 = 1 + IT_0579;
    const ccomplex_t IT_0581 = IT_0000*IT_0580;
    const ccomplex_t IT_0582 = (-8)*IT_0581;
    const ccomplex_t IT_0583 = IT_0151*IT_0269;
    const ccomplex_t IT_0584 = (-10)*s_14;
    const ccomplex_t IT_0585 = IT_0583 + IT_0584;
    const ccomplex_t IT_0586 = m_N_1*IT_0585;
    const ccomplex_t IT_0587 = 2*IT_0586;
    const ccomplex_t IT_0588 = 0.25*conj(IT_0220);
    const ccomplex_t IT_0589 = (-0.666666666666667)*IT_0190 +
       0.166666666666667*IT_0286;
    const ccomplex_t IT_0590 = s_24 + IT_0589;
    const ccomplex_t IT_0591 = m_N_1*IT_0590;
    const ccomplex_t IT_0592 = 6*IT_0591;
    const ccomplex_t IT_0593 = 2*IT_0592;
    const ccomplex_t IT_0594 = 0.25*conj(IT_0222);
    const ccomplex_t IT_0595 = 0.25*IT_0382 + IT_0423*IT_0444 + IT_0384
      *IT_0448 + conj(IT_0311)*IT_0459 + IT_0381*IT_0468 + IT_0441*IT_0489 +
       IT_0291*IT_0513 + IT_0407*IT_0532 + IT_0416*IT_0544 + IT_0396*IT_0564 +
       IT_0408*IT_0566 + IT_0424*IT_0575 + IT_0451*IT_0582 + IT_0587*IT_0588 +
       IT_0554*IT_0593 + IT_0275*IT_0594;
    const ccomplex_t IT_0596 = IT_0381*IT_0384 + IT_0157*IT_0444 + conj
      (IT_0375)*IT_0459 + IT_0423*IT_0467 + IT_0448*IT_0468 + IT_0436*IT_0489 +
       IT_0275*IT_0513 + IT_0416*IT_0532 + IT_0407*IT_0544 + IT_0451*IT_0564 +
       IT_0424*IT_0566 + IT_0408*IT_0575 + IT_0396*IT_0582 + IT_0554*IT_0587 +
       IT_0588*IT_0593 + IT_0291*IT_0594;
    const ccomplex_t IT_0597 = IT_0151*IT_0355;
    const ccomplex_t IT_0598 = 0.25*IT_0597;
    const ccomplex_t IT_0599 = (-0.5)*IT_0572;
    const ccomplex_t IT_0600 = IT_0574 + IT_0598 + IT_0599;
    const ccomplex_t IT_0601 = 8*IT_0600;
    const ccomplex_t IT_0602 = s_13*s_34*IT_0184;
    const ccomplex_t IT_0603 = (-0.166666666666667)*IT_0300 + (
      -0.666666666666667)*IT_0602;
    const ccomplex_t IT_0604 = s_14 + IT_0603;
    const ccomplex_t IT_0605 = m_N_1*IT_0604;
    const ccomplex_t IT_0606 = (-6)*IT_0605;
    const ccomplex_t IT_0607 = 2*IT_0606;
    const ccomplex_t IT_0608 = IT_0184*IT_0282;
    const ccomplex_t IT_0609 = (-2)*IT_0187 + (-0.5)*IT_0608;
    const ccomplex_t IT_0610 = IT_0185 + IT_0609;
    const ccomplex_t IT_0611 = IT_0267*IT_0610;
    const ccomplex_t IT_0612 = (-2)*IT_0611;
    const ccomplex_t IT_0613 = 2*IT_0612;
    const ccomplex_t IT_0614 = 0.25*IT_0443 + IT_0395*IT_0444 + IT_0384
      *IT_0445 + IT_0381*IT_0467 + IT_0436*IT_0477 + IT_0416*IT_0489 + IT_0441
      *IT_0512 + IT_0310*IT_0513 + conj(IT_0376)*IT_0521 + IT_0424*IT_0564 +
       IT_0451*IT_0566 + IT_0408*IT_0582 + IT_0305*IT_0594 + IT_0396*IT_0601 +
       IT_0588*IT_0607 + IT_0554*IT_0613;
    const ccomplex_t IT_0615 = IT_0157*IT_0384 + IT_0381*IT_0444 + IT_0395
      *IT_0467 + IT_0445*IT_0468 + IT_0441*IT_0477 + IT_0407*IT_0489 + IT_0436
      *IT_0512 + IT_0305*IT_0513 + conj(IT_0335)*IT_0521 + IT_0408*IT_0564 +
       IT_0396*IT_0566 + IT_0424*IT_0582 + IT_0310*IT_0594 + IT_0451*IT_0601 +
       IT_0554*IT_0607 + IT_0588*IT_0613;
    const ccomplex_t IT_0616 = IT_0440*IT_0444 + IT_0384*IT_0450 + conj
      (IT_0266)*IT_0459 + conj(IT_0220)*IT_0466 + IT_0435*IT_0467 + IT_0449
      *IT_0468 + IT_0451*IT_0477 + IT_0416*IT_0481 + IT_0424*IT_0489 + IT_0407
      *IT_0496 + IT_0441*IT_0499 + IT_0436*IT_0503 + IT_0396*IT_0512 + IT_0324
      *IT_0594;
    const ccomplex_t IT_0617 = IT_0184*IT_0568;
    const ccomplex_t IT_0618 = -IT_0617;
    const ccomplex_t IT_0619 = IT_0325 + IT_0357 + IT_0359 + IT_0618;
    const ccomplex_t IT_0620 = (-4)*IT_0619;
    const ccomplex_t IT_0621 = 0.25*IT_0373 + 0.25*IT_0374 + 0.25*IT_0377 +
       0.25*IT_0378 + conj(IT_0375)*IT_0466 + IT_0371*IT_0513 + IT_0407*IT_0553 
      + IT_0363*IT_0554 + IT_0424*IT_0587 + IT_0408*IT_0593 + IT_0334*IT_0594 +
       IT_0396*IT_0607 + IT_0451*IT_0613 + IT_0588*IT_0620;
    const ccomplex_t IT_0622 = IT_0406*IT_0444 + IT_0384*IT_0446 + IT_0415
      *IT_0467 + IT_0447*IT_0468 + IT_0436*IT_0481 + IT_0451*IT_0489 + IT_0441
      *IT_0496 + conj(IT_0138)*IT_0521 + IT_0416*IT_0525 + IT_0424*IT_0532 +
       IT_0407*IT_0535 + IT_0408*IT_0544 + IT_0553*IT_0588 + IT_0349*IT_0594;
    const ccomplex_t IT_0623 = 0.25*IT_0240 + 0.25*IT_0249 + 0.25*IT_0364 +
       0.25*IT_0365 + conj(IT_0311)*IT_0466 + IT_0334*IT_0513 + IT_0416*IT_0553 
      + IT_0408*IT_0587 + IT_0363*IT_0588 + IT_0424*IT_0593 + IT_0371*IT_0594 +
       IT_0451*IT_0607 + IT_0396*IT_0613 + IT_0554*IT_0620;
    const ccomplex_t IT_0624 = conj(IT_0216) + conj(IT_0280);
    const ccomplex_t IT_0625 = 0.25*IT_0624;
    const ccomplex_t IT_0626 = 0.25*IT_0180;
    const ccomplex_t IT_0627 = 0.25*IT_0223;
    const ccomplex_t IT_0628 = conj(IT_0150)*IT_0380 + IT_0416*IT_0446 +
       IT_0407*IT_0447 + IT_0408*IT_0448 + IT_0441*IT_0449 + IT_0436*IT_0450 +
       IT_0445*IT_0451 + IT_0383*IT_0467 + IT_0381*IT_0625 + IT_0238*IT_0626 +
       IT_0195*IT_0627;
    const ccomplex_t IT_0629 = conj(IT_0248)*IT_0380 + 0.25*IT_0139*IT_0381 +
       IT_0407*IT_0415 + IT_0406*IT_0416 + IT_0408*IT_0423 + IT_0436*IT_0440 +
       IT_0435*IT_0441 + IT_0395*IT_0451 + IT_0383*IT_0468 + IT_0157*IT_0625 +
       IT_0236*IT_0626 + IT_0217*IT_0627;
    const ccomplex_t IT_0630 = IT_0139*IT_0158 + IT_0180*(IT_0196 + IT_0218) +
       IT_0223*(IT_0237 + IT_0239) + IT_0179*IT_0372 + IT_0222*IT_0379 + 4
      *IT_0216*IT_0442 + 4*IT_0041*IT_0452 + 4*IT_0311*IT_0514 + 4*IT_0335
      *IT_0555 + 4*IT_0280*IT_0595 + 4*IT_0266*IT_0596 + 4*IT_0138*IT_0614 + 4
      *IT_0296*IT_0615 + 4*IT_0375*IT_0616 + 4*IT_0220*IT_0621 + 4*IT_0376
      *IT_0622 + 4*IT_0172*IT_0623 + 4*IT_0150*IT_0628 + 4*IT_0248*IT_0629;
    return create_ccomplex_return(IT_0630);
}

