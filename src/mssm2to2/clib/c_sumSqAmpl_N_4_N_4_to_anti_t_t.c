#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = m_t*m_N_4;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = pow(m_t, 2);
    const ccomplex_t IT_0006 = pow(m_N_4, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = sin(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_t*conj(N_u4)*e_em*IT_0021*U_st_10;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0012 + 3*IT_0017 + 3
      *IT_0024);
    const ccomplex_t IT_0026 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_0027 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = m_t*N_u4*e_em*IT_0021*conj(U_st_10);
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0029 + IT_0032 +
       0.333333333333333*IT_0035);
    const ccomplex_t IT_0037 = 3*IT_0036;
    const ccomplex_t IT_0038 = 0.166666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0026*IT_0038;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0043 = N_B4*e_em*conj(U_st_01);
    const ccomplex_t IT_0044 = IT_0009*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = N_W4*e_em*conj(U_st_01);
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = m_t*N_u4*e_em*IT_0021*conj(U_st_11);
    const ccomplex_t IT_0050 = IT_0020*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0045 + 3*IT_0048 + 3
      *IT_0051);
    const ccomplex_t IT_0053 = 0.166666666666667*IT_0052;
    const ccomplex_t IT_0054 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0055 = IT_0009*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0058 = IT_0014*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = m_t*conj(N_u4)*e_em*IT_0021*U_st_11;
    const ccomplex_t IT_0061 = IT_0020*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0056 + 3*IT_0059 + 3
      *IT_0062);
    const ccomplex_t IT_0064 = 0.166666666666667*IT_0063;
    const ccomplex_t IT_0065 = IT_0053*IT_0064;
    const ccomplex_t IT_0066 = IT_0042*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0041 + -IT_0067;
    const ccomplex_t IT_0069 = cos(alpha);
    const ccomplex_t IT_0070 = IT_0009*IT_0069;
    const ccomplex_t IT_0071 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = sin(alpha);
    const ccomplex_t IT_0074 = IT_0009*IT_0073;
    const ccomplex_t IT_0075 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = IT_0014*IT_0069;
    const ccomplex_t IT_0078 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = IT_0014*IT_0073;
    const ccomplex_t IT_0081 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0072 + IT_0076 + 
      -IT_0079 + -IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_t*e_em*IT_0014*IT_0019
      *IT_0021*IT_0069;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0084*IT_0089;
    const ccomplex_t IT_0091 = IT_0070*IT_0075;
    const ccomplex_t IT_0092 = IT_0071*IT_0074;
    const ccomplex_t IT_0093 = IT_0077*IT_0081;
    const ccomplex_t IT_0094 = IT_0078*IT_0080;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + -IT_0092 + 
      -IT_0093 + IT_0094);
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_t*e_em*IT_0014*IT_0019
      *IT_0021*IT_0073;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = cpow((-2)*s_12 + (-2)*IT_0006 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0095*IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0090 + (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = 3*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0090 + 0.5*IT_0101;
    const ccomplex_t IT_0105 = (-3)*IT_0104;
    const ccomplex_t IT_0106 = IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = 3*conj(IT_0102);
    const ccomplex_t IT_0108 = (-3)*conj(IT_0104);
    const ccomplex_t IT_0109 = IT_0107 + IT_0108;
    const ccomplex_t IT_0110 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0111 = IT_0009*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = m_t*N_u4*e_em*IT_0021*U_st_00;
    const ccomplex_t IT_0114 = IT_0020*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0112 + (-0.75)*IT_0115);
    const ccomplex_t IT_0117 = 0.666666666666667*IT_0116;
    const ccomplex_t IT_0118 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0119 = IT_0009*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = m_t*conj(N_u4)*e_em*IT_0021*conj(U_st_00);
    const ccomplex_t IT_0122 = IT_0020*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + (-0.75)*IT_0123);
    const ccomplex_t IT_0125 = 0.666666666666667*IT_0124;
    const ccomplex_t IT_0126 = IT_0117*IT_0125;
    const ccomplex_t IT_0127 = IT_0007*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B4)*e_em*conj(U_st_11);
    const ccomplex_t IT_0130 = IT_0009*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = m_t*conj(N_u4)*e_em*IT_0021*conj(U_st_01);
    const ccomplex_t IT_0133 = IT_0020*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + (-0.75)*IT_0134);
    const ccomplex_t IT_0136 = 0.666666666666667*IT_0135;
    const ccomplex_t IT_0137 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0138 = IT_0009*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = m_t*N_u4*e_em*IT_0021*U_st_01;
    const ccomplex_t IT_0141 = IT_0020*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + (-0.75)*IT_0142);
    const ccomplex_t IT_0144 = 0.666666666666667*IT_0143;
    const ccomplex_t IT_0145 = IT_0136*IT_0144;
    const ccomplex_t IT_0146 = IT_0042*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = -IT_0128 + -IT_0147;
    const ccomplex_t IT_0149 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0150 = IT_0070*IT_0149;
    const ccomplex_t IT_0151 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0152 = IT_0074*IT_0151;
    const ccomplex_t IT_0153 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0154 = IT_0077*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0156 = IT_0080*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0150 + IT_0152 + 
      -IT_0154 + -IT_0156);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0089*IT_0158;
    const ccomplex_t IT_0160 = IT_0070*IT_0151;
    const ccomplex_t IT_0161 = IT_0074*IT_0149;
    const ccomplex_t IT_0162 = IT_0077*IT_0155;
    const ccomplex_t IT_0163 = IT_0080*IT_0153;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + -IT_0161 + 
      -IT_0162 + IT_0163);
    const ccomplex_t IT_0165 = IT_0100*IT_0164;
    const ccomplex_t IT_0166 = (-0.5)*IT_0159 + (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = 3*IT_0166;
    const ccomplex_t IT_0168 = 0.5*IT_0159 + 0.5*IT_0165;
    const ccomplex_t IT_0169 = (-3)*IT_0168;
    const ccomplex_t IT_0170 = IT_0167 + IT_0169;
    const ccomplex_t IT_0171 = 3*conj(IT_0166);
    const ccomplex_t IT_0172 = (-3)*conj(IT_0168);
    const ccomplex_t IT_0173 = IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = -s_24;
    const ccomplex_t IT_0175 = s_23 + IT_0174;
    const ccomplex_t IT_0176 = IT_0000*IT_0175;
    const ccomplex_t IT_0177 = (-2)*IT_0176;
    const ccomplex_t IT_0178 = 4*IT_0003;
    const ccomplex_t IT_0179 = IT_0009*IT_0013;
    const ccomplex_t IT_0180 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = IT_0008*IT_0014;
    const ccomplex_t IT_0183 = IT_0180*IT_0182;
    const ccomplex_t IT_0184 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0185 = IT_0182*IT_0184;
    const ccomplex_t IT_0186 = IT_0179*IT_0184;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0181 + IT_0183 + 
      -IT_0185 + -IT_0186);
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0013;
    const ccomplex_t IT_0193 = (-0.666666666666667)*IT_0192;
    const ccomplex_t IT_0194 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0195 = IT_0193*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = IT_0191*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0197;
    const ccomplex_t IT_0199 = e_em*IT_0179;
    const ccomplex_t IT_0200 = e_em*IT_0182;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0199 + (-3)*IT_0200);
    const ccomplex_t IT_0202 = (-0.166666666666667)*IT_0201;
    const ccomplex_t IT_0203 = IT_0194*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = IT_0191*IT_0204;
    const ccomplex_t IT_0206 = 0.5*IT_0205;
    const ccomplex_t IT_0207 = (-0.5)*IT_0190;
    const ccomplex_t IT_0208 = IT_0204*IT_0207;
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = IT_0196*IT_0207;
    const ccomplex_t IT_0211 = 0.5*IT_0210;
    const ccomplex_t IT_0212 = (-4)*IT_0176;
    const ccomplex_t IT_0213 = (-2)*IT_0003;
    const ccomplex_t IT_0214 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0215 = IT_0126*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0145*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = IT_0216 + IT_0219;
    const ccomplex_t IT_0221 = 3*IT_0104;
    const ccomplex_t IT_0222 = (-3)*IT_0102;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = 3*conj(IT_0104);
    const ccomplex_t IT_0225 = (-3)*conj(IT_0102);
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = IT_0039*IT_0214;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = IT_0065*IT_0217;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = 3*IT_0168;
    const ccomplex_t IT_0233 = (-3)*IT_0166;
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = 3*conj(IT_0168);
    const ccomplex_t IT_0236 = (-3)*conj(IT_0166);
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = 2*IT_0176;
    const ccomplex_t IT_0239 = (-4)*IT_0003;
    const ccomplex_t IT_0240 = (-0.5)*IT_0208;
    const ccomplex_t IT_0241 = (-0.5)*IT_0210;
    const ccomplex_t IT_0242 = conj(IT_0240) + conj(IT_0241);
    const ccomplex_t IT_0243 = IT_0240 + IT_0241;
    const ccomplex_t IT_0244 = (-0.5)*IT_0197;
    const ccomplex_t IT_0245 = (-0.5)*IT_0205;
    const ccomplex_t IT_0246 = conj(IT_0244) + conj(IT_0245);
    const ccomplex_t IT_0247 = IT_0244 + IT_0245;
    const ccomplex_t IT_0248 = 4*IT_0176;
    const ccomplex_t IT_0249 = s_12*IT_0005;
    const ccomplex_t IT_0250 = s_13*s_24;
    const ccomplex_t IT_0251 = s_12*s_34;
    const ccomplex_t IT_0252 = s_14*s_23;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = IT_0249 + IT_0250 + IT_0251 + IT_0253;
    const ccomplex_t IT_0255 = IT_0038*IT_0117;
    const ccomplex_t IT_0256 = IT_0007*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = IT_0053*IT_0144;
    const ccomplex_t IT_0259 = IT_0042*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = IT_0257 + IT_0260;
    const ccomplex_t IT_0262 = cos(beta);
    const ccomplex_t IT_0263 = IT_0009*IT_0262;
    const ccomplex_t IT_0264 = IT_0075*IT_0263;
    const ccomplex_t IT_0265 = IT_0009*IT_0018;
    const ccomplex_t IT_0266 = IT_0071*IT_0265;
    const ccomplex_t IT_0267 = IT_0014*IT_0262;
    const ccomplex_t IT_0268 = IT_0081*IT_0267;
    const ccomplex_t IT_0269 = IT_0014*IT_0018;
    const ccomplex_t IT_0270 = IT_0078*IT_0269;
    const ccomplex_t IT_0271 = -IT_0264 + -IT_0266 + IT_0268 + IT_0270;
    const ccomplex_t IT_0272 = m_t*e_em*IT_0014*IT_0021;
    const ccomplex_t IT_0273 = (-0.5)*IT_0272;
    const ccomplex_t IT_0274 = pow(m_Z, 2);
    const ccomplex_t IT_0275 = cpow((-2)*s_12 + (-2)*IT_0006 + IT_0274 + 
      -reg_prop, -1);
    const ccomplex_t IT_0276 = IT_0273*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = IT_0271*IT_0277;
    const ccomplex_t IT_0279 = IT_0071*IT_0263;
    const ccomplex_t IT_0280 = IT_0075*IT_0265;
    const ccomplex_t IT_0281 = IT_0078*IT_0267;
    const ccomplex_t IT_0282 = IT_0081*IT_0269;
    const ccomplex_t IT_0283 = -IT_0279 + IT_0280 + IT_0281 + -IT_0282;
    const ccomplex_t IT_0284 = m_t*e_em*IT_0014*IT_0019*IT_0021*IT_0262;
    const ccomplex_t IT_0285 = (-0.5)*IT_0284;
    const ccomplex_t IT_0286 = cpow((-2)*s_12 + (-2)*IT_0006 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0287 = IT_0285*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*IT_0287;
    const ccomplex_t IT_0289 = IT_0283*IT_0288;
    const ccomplex_t IT_0290 = 0.5*IT_0278 + 0.5*IT_0289;
    const ccomplex_t IT_0291 = 6*IT_0290;
    const ccomplex_t IT_0292 = (-0.5)*IT_0278 + (-0.5)*IT_0289;
    const ccomplex_t IT_0293 = (-6)*IT_0292;
    const ccomplex_t IT_0294 = IT_0291 + IT_0293;
    const ccomplex_t IT_0295 = 6*conj(IT_0290);
    const ccomplex_t IT_0296 = (-6)*conj(IT_0292);
    const ccomplex_t IT_0297 = IT_0295 + IT_0296;
    const ccomplex_t IT_0298 = IT_0026*IT_0125;
    const ccomplex_t IT_0299 = IT_0007*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0064*IT_0136;
    const ccomplex_t IT_0302 = IT_0042*IT_0301;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*IT_0302;
    const ccomplex_t IT_0304 = IT_0300 + IT_0303;
    const ccomplex_t IT_0305 = IT_0149*IT_0263;
    const ccomplex_t IT_0306 = IT_0151*IT_0265;
    const ccomplex_t IT_0307 = IT_0153*IT_0267;
    const ccomplex_t IT_0308 = IT_0155*IT_0269;
    const ccomplex_t IT_0309 = IT_0305 + -IT_0306 + -IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = IT_0288*IT_0309;
    const ccomplex_t IT_0311 = IT_0151*IT_0263;
    const ccomplex_t IT_0312 = IT_0149*IT_0265;
    const ccomplex_t IT_0313 = IT_0155*IT_0267;
    const ccomplex_t IT_0314 = IT_0153*IT_0269;
    const ccomplex_t IT_0315 = IT_0311 + IT_0312 + -IT_0313 + -IT_0314;
    const ccomplex_t IT_0316 = IT_0277*IT_0315;
    const ccomplex_t IT_0317 = (-0.5)*IT_0310 + (-0.5)*IT_0316;
    const ccomplex_t IT_0318 = 6*IT_0317;
    const ccomplex_t IT_0319 = 0.5*IT_0310 + 0.5*IT_0316;
    const ccomplex_t IT_0320 = (-6)*IT_0319;
    const ccomplex_t IT_0321 = IT_0318 + IT_0320;
    const ccomplex_t IT_0322 = 6*conj(IT_0317);
    const ccomplex_t IT_0323 = (-6)*conj(IT_0319);
    const ccomplex_t IT_0324 = IT_0322 + IT_0323;
    const ccomplex_t IT_0325 = -IT_0250;
    const ccomplex_t IT_0326 = IT_0249 + IT_0251 + IT_0252 + IT_0325;
    const ccomplex_t IT_0327 = IT_0214*IT_0255;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = IT_0217*IT_0258;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*IT_0329;
    const ccomplex_t IT_0331 = -IT_0328 + -IT_0330;
    const ccomplex_t IT_0332 = 6*IT_0292;
    const ccomplex_t IT_0333 = (-6)*IT_0290;
    const ccomplex_t IT_0334 = IT_0332 + IT_0333;
    const ccomplex_t IT_0335 = 6*conj(IT_0292);
    const ccomplex_t IT_0336 = (-6)*conj(IT_0290);
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = IT_0214*IT_0298;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*IT_0338;
    const ccomplex_t IT_0340 = IT_0217*IT_0301;
    const ccomplex_t IT_0341 = (0 + _Complex_I*1)*IT_0340;
    const ccomplex_t IT_0342 = -IT_0339 + -IT_0341;
    const ccomplex_t IT_0343 = 6*IT_0319;
    const ccomplex_t IT_0344 = (-6)*IT_0317;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = 6*conj(IT_0319);
    const ccomplex_t IT_0347 = (-6)*conj(IT_0317);
    const ccomplex_t IT_0348 = IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = s_23 + s_24;
    const ccomplex_t IT_0350 = IT_0000*IT_0349;
    const ccomplex_t IT_0351 = conj(IT_0209) + conj(IT_0245);
    const ccomplex_t IT_0352 = 12*IT_0292;
    const ccomplex_t IT_0353 = (-12)*IT_0290;
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = IT_0209 + IT_0245;
    const ccomplex_t IT_0356 = 12*conj(IT_0292);
    const ccomplex_t IT_0357 = (-12)*conj(IT_0290);
    const ccomplex_t IT_0358 = IT_0356 + IT_0357;
    const ccomplex_t IT_0359 = conj(IT_0211) + conj(IT_0244);
    const ccomplex_t IT_0360 = 12*IT_0290;
    const ccomplex_t IT_0361 = (-12)*IT_0292;
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = IT_0211 + IT_0244;
    const ccomplex_t IT_0364 = 12*conj(IT_0290);
    const ccomplex_t IT_0365 = (-12)*conj(IT_0292);
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = conj(IT_0206) + conj(IT_0240);
    const ccomplex_t IT_0368 = 12*IT_0317;
    const ccomplex_t IT_0369 = (-12)*IT_0319;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = IT_0206 + IT_0240;
    const ccomplex_t IT_0372 = 12*conj(IT_0317);
    const ccomplex_t IT_0373 = (-12)*conj(IT_0319);
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = conj(IT_0198) + conj(IT_0241);
    const ccomplex_t IT_0376 = 12*IT_0319;
    const ccomplex_t IT_0377 = (-12)*IT_0317;
    const ccomplex_t IT_0378 = IT_0376 + IT_0377;
    const ccomplex_t IT_0379 = IT_0198 + IT_0241;
    const ccomplex_t IT_0380 = 12*conj(IT_0319);
    const ccomplex_t IT_0381 = (-12)*conj(IT_0317);
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = conj(IT_0148)*IT_0294 + IT_0148*IT_0297 + conj
      (IT_0068)*IT_0321 + IT_0068*IT_0324 + conj(IT_0231)*IT_0334 + IT_0231
      *IT_0337 + conj(IT_0220)*IT_0345 + IT_0220*IT_0348 + IT_0351*IT_0354 +
       IT_0355*IT_0358 + IT_0359*IT_0362 + IT_0363*IT_0366 + IT_0367*IT_0370 +
       IT_0371*IT_0374 + IT_0375*IT_0378 + IT_0379*IT_0382;
    const ccomplex_t IT_0384 = -IT_0251;
    const ccomplex_t IT_0385 = IT_0250 + IT_0252 + IT_0384;
    const ccomplex_t IT_0386 = 6*IT_0261;
    const ccomplex_t IT_0387 = 6*conj(IT_0261);
    const ccomplex_t IT_0388 = s_13 + s_14;
    const ccomplex_t IT_0389 = IT_0000*IT_0388;
    const ccomplex_t IT_0390 = 6*IT_0148;
    const ccomplex_t IT_0391 = 6*conj(IT_0148);
    const ccomplex_t IT_0392 = 6*IT_0231;
    const ccomplex_t IT_0393 = 6*conj(IT_0231);
    const ccomplex_t IT_0394 = conj(IT_0068)*IT_0294 + IT_0068*IT_0297 + conj
      (IT_0220)*IT_0334 + IT_0220*IT_0337 + IT_0362*IT_0367 + IT_0359*IT_0370 +
       IT_0366*IT_0371 + IT_0363*IT_0374 + IT_0354*IT_0375 + IT_0351*IT_0378 +
       IT_0358*IT_0379 + IT_0355*IT_0382 + (-6)*conj(IT_0317)*(IT_0231 + (
      -0.166666666666667)*IT_0390) + (-6)*IT_0317*(conj(IT_0231) + (
      -0.166666666666667)*IT_0391) + (-6)*conj(IT_0319)*(IT_0148 + (
      -0.166666666666667)*IT_0392) + (-6)*IT_0319*(conj(IT_0148) + (
      -0.166666666666667)*IT_0393);
    const ccomplex_t IT_0395 = m_t*s_14*m_N_4;
    const ccomplex_t IT_0396 = 24*IT_0209;
    const ccomplex_t IT_0397 = 24*IT_0245;
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = 12*IT_0231;
    const ccomplex_t IT_0400 = IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = 24*conj(IT_0209);
    const ccomplex_t IT_0402 = 24*conj(IT_0245);
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = 12*conj(IT_0231);
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = 12*IT_0331;
    const ccomplex_t IT_0407 = IT_0386 + IT_0406;
    const ccomplex_t IT_0408 = 12*conj(IT_0331);
    const ccomplex_t IT_0409 = IT_0387 + IT_0408;
    const ccomplex_t IT_0410 = 12*IT_0261;
    const ccomplex_t IT_0411 = 12*conj(IT_0261);
    const ccomplex_t IT_0412 = 12*IT_0304;
    const ccomplex_t IT_0413 = 12*conj(IT_0304);
    const ccomplex_t IT_0414 = conj(IT_0304)*IT_0392 + IT_0304*IT_0393 + conj
      (IT_0342)*IT_0400 + IT_0342*IT_0405 + conj(IT_0220)*IT_0407 + IT_0220
      *IT_0409 + IT_0375*(24*IT_0331 + IT_0410) + IT_0379*(24*conj(IT_0331) +
       IT_0411) + IT_0351*IT_0412 + IT_0355*IT_0413;
    const ccomplex_t IT_0415 = m_t*s_13*m_N_4;
    const ccomplex_t IT_0416 = 12*IT_0148;
    const ccomplex_t IT_0417 = 24*IT_0211;
    const ccomplex_t IT_0418 = 24*IT_0244;
    const ccomplex_t IT_0419 = IT_0416 + IT_0417 + IT_0418;
    const ccomplex_t IT_0420 = 12*conj(IT_0148);
    const ccomplex_t IT_0421 = 24*conj(IT_0211);
    const ccomplex_t IT_0422 = 24*conj(IT_0244);
    const ccomplex_t IT_0423 = IT_0420 + IT_0421 + IT_0422;
    const ccomplex_t IT_0424 = 24*IT_0206;
    const ccomplex_t IT_0425 = 24*IT_0240;
    const ccomplex_t IT_0426 = IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = 24*conj(IT_0206);
    const ccomplex_t IT_0428 = 24*conj(IT_0240);
    const ccomplex_t IT_0429 = IT_0427 + IT_0428;
    const ccomplex_t IT_0430 = 6*IT_0331;
    const ccomplex_t IT_0431 = IT_0410 + IT_0430;
    const ccomplex_t IT_0432 = 6*conj(IT_0331);
    const ccomplex_t IT_0433 = IT_0411 + IT_0432;
    const ccomplex_t IT_0434 = 12*IT_0342;
    const ccomplex_t IT_0435 = 12*conj(IT_0342);
    const ccomplex_t IT_0436 = conj(IT_0342)*IT_0390 + IT_0342*IT_0391 +
       IT_0367*IT_0406 + IT_0371*IT_0408 + conj(IT_0304)*IT_0419 + IT_0304
      *IT_0423 + conj(IT_0261)*IT_0426 + IT_0261*IT_0429 + conj(IT_0068)*IT_0431
       + IT_0068*IT_0433 + IT_0359*IT_0434 + IT_0363*IT_0435;
    const ccomplex_t IT_0437 = m_t*s_24*m_N_4;
    const ccomplex_t IT_0438 = 12*IT_0068;
    const ccomplex_t IT_0439 = IT_0426 + IT_0438;
    const ccomplex_t IT_0440 = 6*IT_0068;
    const ccomplex_t IT_0441 = 6*conj(IT_0068);
    const ccomplex_t IT_0442 = 12*conj(IT_0068);
    const ccomplex_t IT_0443 = IT_0359*(24*IT_0261 + IT_0406) + IT_0363*(24
      *conj(IT_0261) + IT_0408) + conj(IT_0148)*IT_0431 + IT_0148*IT_0433 +
       IT_0367*IT_0434 + IT_0371*IT_0435 + conj(IT_0304)*IT_0439 + conj(IT_0342)
      *IT_0440 + IT_0342*IT_0441 + IT_0304*(IT_0429 + IT_0442);
    const ccomplex_t IT_0444 = s_34*IT_0006;
    const ccomplex_t IT_0445 = IT_0396 + IT_0397 + IT_0399;
    const ccomplex_t IT_0446 = IT_0401 + IT_0402 + IT_0404;
    const ccomplex_t IT_0447 = 12*IT_0220;
    const ccomplex_t IT_0448 = 12*conj(IT_0220);
    const ccomplex_t IT_0449 = conj(IT_0220)*IT_0390 + IT_0220*IT_0391 + conj
      (IT_0068)*IT_0392 + IT_0068*IT_0393 + conj(IT_0206)*(IT_0396 + IT_0399) +
       IT_0206*(IT_0401 + IT_0404) + IT_0375*IT_0419 + IT_0379*IT_0423 + conj
      (IT_0245)*IT_0424 + IT_0245*IT_0427 + IT_0351*IT_0438 + IT_0355*IT_0442 +
       conj(IT_0240)*IT_0445 + IT_0240*IT_0446 + IT_0359*IT_0447 + IT_0363
      *IT_0448;
    const ccomplex_t IT_0450 = 48*IT_0206;
    const ccomplex_t IT_0451 = 24*IT_0068;
    const ccomplex_t IT_0452 = IT_0450 + IT_0451;
    const ccomplex_t IT_0453 = 48*IT_0240;
    const ccomplex_t IT_0454 = IT_0452 + IT_0453;
    const ccomplex_t IT_0455 = IT_0416 + IT_0417;
    const ccomplex_t IT_0456 = 48*IT_0211;
    const ccomplex_t IT_0457 = 48*conj(IT_0211);
    const ccomplex_t IT_0458 = 48*conj(IT_0206);
    const ccomplex_t IT_0459 = conj(IT_0261)*IT_0410 + conj(IT_0304)*IT_0412 +
       conj(IT_0068)*IT_0439 + conj(IT_0206)*IT_0452 + conj(IT_0240)*IT_0454 +
       conj(IT_0148)*IT_0455 + IT_0359*(24*IT_0148 + IT_0456) + IT_0244*(24*conj
      (IT_0148) + 48*conj(IT_0244) + IT_0457) + IT_0240*IT_0458;
    const ccomplex_t IT_0460 = -s_34;
    const ccomplex_t IT_0461 = IT_0005 + IT_0460;
    const ccomplex_t IT_0462 = s_12*IT_0461;
    const ccomplex_t IT_0463 = (-24)*conj(IT_0168);
    const ccomplex_t IT_0464 = 24*conj(IT_0166);
    const ccomplex_t IT_0465 = IT_0463 + IT_0464;
    const ccomplex_t IT_0466 = (-24)*IT_0104;
    const ccomplex_t IT_0467 = s_34 + IT_0005;
    const ccomplex_t IT_0468 = s_12*IT_0467;
    const ccomplex_t IT_0469 = (-24)*IT_0292;
    const ccomplex_t IT_0470 = (-24)*conj(IT_0292);
    const ccomplex_t IT_0471 = 6*conj(IT_0102);
    const ccomplex_t IT_0472 = (-6)*conj(IT_0104);
    const ccomplex_t IT_0473 = IT_0471 + IT_0472;
    const ccomplex_t IT_0474 = 6*IT_0166;
    const ccomplex_t IT_0475 = (-6)*IT_0168;
    const ccomplex_t IT_0476 = IT_0474 + IT_0475;
    const ccomplex_t IT_0477 = (-6)*IT_0104;
    const ccomplex_t IT_0478 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0479 = conj(IT_0231)*IT_0390 + IT_0231*IT_0391 +
       IT_0351*IT_0416 + IT_0355*IT_0420 + IT_0375*(IT_0424 + IT_0425 + IT_0438)
       + conj(IT_0220)*IT_0440 + IT_0220*IT_0441 + IT_0379*(IT_0427 + IT_0428 +
       IT_0442) + IT_0359*IT_0445 + IT_0363*IT_0446 + IT_0367*IT_0447 + IT_0371
      *IT_0448;
    const ccomplex_t IT_0480 = m_t*s_23*m_N_4;
    const ccomplex_t IT_0481 = 24*IT_0241;
    const ccomplex_t IT_0482 = 24*IT_0198 + IT_0447 + IT_0481;
    const ccomplex_t IT_0483 = 6*conj(IT_0304);
    const ccomplex_t IT_0484 = conj(IT_0220)*IT_0304 + 0.166666666666667*conj
      (IT_0331)*IT_0398 + 0.166666666666667*IT_0331*IT_0403 + 0.166666666666667
      *conj(IT_0231)*IT_0407 + 0.166666666666667*IT_0231*IT_0409 +
       0.166666666666667*IT_0351*IT_0410 + 0.166666666666667*IT_0355*IT_0411 +
       0.166666666666667*IT_0375*IT_0412 + 0.166666666666667*IT_0379*IT_0413 +
       0.166666666666667*IT_0342*(24*conj(IT_0198) + 24*conj(IT_0241) + IT_0448)
       + 0.166666666666667*conj(IT_0342)*IT_0482 + 0.166666666666667*IT_0220
      *IT_0483;
    const ccomplex_t IT_0485 = IT_0006*IT_0461;
    const ccomplex_t IT_0486 = 6*conj(IT_0104);
    const ccomplex_t IT_0487 = (-6)*conj(IT_0102);
    const ccomplex_t IT_0488 = IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = 6*IT_0168;
    const ccomplex_t IT_0490 = (-6)*IT_0166;
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = 6*conj(IT_0168);
    const ccomplex_t IT_0493 = (-6)*conj(IT_0166);
    const ccomplex_t IT_0494 = IT_0492 + IT_0493;
    const ccomplex_t IT_0495 = 6*IT_0104;
    const ccomplex_t IT_0496 = conj(IT_0104)*((-24)*IT_0166 + 24*IT_0168) +
       conj(IT_0102)*(24*IT_0166 + (-24)*IT_0168) + 24*IT_0104*conj(IT_0168) + 6
      *conj(IT_0166)*(IT_0261 + 0.166666666666667*IT_0466) + IT_0304*IT_0473 +
       conj(IT_0261)*IT_0476 + conj(IT_0304)*IT_0477 + IT_0261*IT_0478 + (-6)
      *IT_0102*(conj(IT_0342) + (-0.166666666666667)*IT_0465 + (
      -0.166666666666667)*IT_0483) + IT_0342*IT_0488 + conj(IT_0331)*IT_0491 +
       IT_0331*IT_0494 + conj(IT_0342)*IT_0495;
    const ccomplex_t IT_0497 = IT_0005*IT_0006;
    const ccomplex_t IT_0498 = 48*IT_0209;
    const ccomplex_t IT_0499 = 24*IT_0231;
    const ccomplex_t IT_0500 = IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = 48*conj(IT_0198);
    const ccomplex_t IT_0502 = 48*conj(IT_0209);
    const ccomplex_t IT_0503 = 48*IT_0198;
    const ccomplex_t IT_0504 = 48*IT_0241;
    const ccomplex_t IT_0505 = 48*IT_0245;
    const ccomplex_t IT_0506 = conj(IT_0220)*IT_0400 + IT_0220*IT_0405 + conj
      (IT_0342)*IT_0407 + IT_0342*IT_0409 + IT_0068*(IT_0420 + IT_0421) + conj
      (IT_0148)*IT_0426 + IT_0148*IT_0429 + conj(IT_0304)*IT_0431 + IT_0304
      *IT_0433 + conj(IT_0244)*IT_0454 + conj(IT_0068)*IT_0455 + IT_0367*IT_0456
       + IT_0371*IT_0457 + IT_0244*(24*conj(IT_0068) + 48*conj(IT_0240) +
       IT_0458) + conj(IT_0198)*IT_0500 + IT_0245*IT_0501 + IT_0379*(24*conj
      (IT_0231) + IT_0502) + conj(IT_0245)*(IT_0503 + IT_0504) + conj(IT_0241)*
      (IT_0498 + IT_0499 + IT_0505);
    const ccomplex_t IT_0507 = 24*IT_0220;
    const ccomplex_t IT_0508 = IT_0503 + IT_0507;
    const ccomplex_t IT_0509 = conj(IT_0231)*IT_0400 + conj(IT_0331)*IT_0406 +
       conj(IT_0342)*IT_0434 + conj(IT_0220)*IT_0482 + conj(IT_0209)*IT_0500 +
       IT_0241*IT_0501 + IT_0245*IT_0502 + conj(IT_0245)*(IT_0498 + IT_0499 +
       IT_0505) + conj(IT_0198)*IT_0508 + conj(IT_0241)*(IT_0504 + IT_0508);
    const ccomplex_t IT_0510 = (-0.166666666666667)*IT_0304;
    const ccomplex_t IT_0511 = (-0.166666666666667)*conj(IT_0304);
    const ccomplex_t IT_0512 = IT_0006*IT_0467;
    const ccomplex_t IT_0513 = IT_0337*IT_0342 + IT_0334*conj(IT_0342) + (-6)
      *conj(IT_0319)*(IT_0261 + 4*IT_0290 + (-4)*IT_0292 + (-0.166666666666667)
      *IT_0430) + (-6)*IT_0319*(conj(IT_0261) + 4*conj(IT_0290) + (-4)*conj
      (IT_0292) + (-0.166666666666667)*IT_0432) + 24*conj(IT_0317)*(IT_0290 + 
      -1./4*IT_0331 + 0.0416666666666667*IT_0386 + 0.0416666666666667*IT_0469) +
       24*IT_0317*(conj(IT_0290) + -1./4*conj(IT_0331) + 0.0416666666666667
      *IT_0387 + 0.0416666666666667*IT_0470) + (-6)*IT_0297*IT_0510 + (-6)
      *IT_0294*IT_0511;
    const ccomplex_t IT_0514 = IT_0004*(conj(IT_0068)*IT_0106 + IT_0068
      *IT_0109 + conj(IT_0148)*IT_0170 + IT_0148*IT_0173) + (IT_0109*IT_0148 +
       IT_0106*conj(IT_0148) + conj(IT_0068)*IT_0170 + IT_0068*IT_0173)*IT_0177 
      + IT_0178*(IT_0109*(IT_0198 + IT_0206) + IT_0106*(conj(IT_0198) + conj
      (IT_0206)) + IT_0173*(IT_0209 + IT_0211) + IT_0170*(conj(IT_0209) + conj
      (IT_0211))) + (IT_0173*(IT_0198 + IT_0206) + IT_0170*(conj(IT_0198) + conj
      (IT_0206)) + IT_0109*(IT_0209 + IT_0211) + IT_0106*(conj(IT_0209) + conj
      (IT_0211)))*IT_0212 + IT_0213*(conj(IT_0220)*IT_0223 + IT_0220*IT_0226 +
       conj(IT_0231)*IT_0234 + IT_0231*IT_0237) + (IT_0226*IT_0231 + IT_0223
      *conj(IT_0231) + conj(IT_0220)*IT_0234 + IT_0220*IT_0237)*IT_0238 +
       IT_0239*(IT_0223*IT_0242 + IT_0226*IT_0243 + IT_0234*IT_0246 + IT_0237
      *IT_0247) + (IT_0234*IT_0242 + IT_0237*IT_0243 + IT_0223*IT_0246 + IT_0226
      *IT_0247)*IT_0248 + IT_0254*(conj(IT_0261)*IT_0294 + IT_0261*IT_0297 +
       conj(IT_0304)*IT_0321 + IT_0304*IT_0324) + IT_0326*(conj(IT_0331)*IT_0334
       + IT_0331*IT_0337 + conj(IT_0342)*IT_0345 + IT_0342*IT_0348) + IT_0350
      *IT_0383 + 6*IT_0385*(conj(IT_0304)*IT_0342 + IT_0304*conj(IT_0342) +
       0.166666666666667*conj(IT_0331)*IT_0386 + 0.166666666666667*IT_0331
      *IT_0387) + IT_0389*IT_0394 + IT_0395*IT_0414 + IT_0415*IT_0436 + IT_0437
      *IT_0443 + IT_0444*IT_0449 + IT_0250*IT_0459 + IT_0462*(24*conj(IT_0102)
      *IT_0104 + IT_0102*((-24)*conj(IT_0102) + 24*conj(IT_0104)) + IT_0166*((
      -24)*conj(IT_0166) + 24*conj(IT_0168)) + IT_0168*IT_0465 + conj(IT_0104)
      *IT_0466) + IT_0468*(24*IT_0292*conj(IT_0292) + IT_0319*((-24)*conj
      (IT_0317) + 24*conj(IT_0319)) + IT_0317*(24*conj(IT_0317) + (-24)*conj
      (IT_0319)) + conj(IT_0290)*IT_0469 + IT_0290*(24*conj(IT_0290) + IT_0470))
       + (IT_0249 + IT_0250 + IT_0253 + IT_0384)*(IT_0331*IT_0473 + conj(IT_0342
      )*IT_0476 + 6*conj(IT_0331)*(IT_0102 + 0.166666666666667*IT_0477) + 6
      *IT_0342*(conj(IT_0166) + 0.166666666666667*IT_0478)) + IT_0249*IT_0479 +
       6*IT_0480*IT_0484 + IT_0485*IT_0496 + IT_0497*IT_0506 + IT_0252*IT_0509 +
       (IT_0249 + IT_0252 + IT_0325 + IT_0384)*(IT_0261*IT_0488 + (-6)*conj
      (IT_0261)*(IT_0102 + (-0.166666666666667)*IT_0495) + (-6)*IT_0494*IT_0510 
      + (-6)*IT_0491*IT_0511) + IT_0512*IT_0513;
    return create_ccomplex_return(IT_0514);
}

