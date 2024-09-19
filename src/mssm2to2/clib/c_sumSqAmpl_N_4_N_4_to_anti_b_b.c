#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_b_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
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
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = m_b*m_N_4;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = pow(m_b, 2);
    const ccomplex_t IT_0006 = pow(m_N_4, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = pow(m_W, -1);
    const ccomplex_t IT_0022 = m_b*N_d4*e_em*IT_0021*conj(U_sb_10);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0012 + (-3)*IT_0017 + 3
      *IT_0024);
    const ccomplex_t IT_0026 = 0.166666666666667*IT_0025;
    const ccomplex_t IT_0027 = m_b*conj(N_d4)*e_em*IT_0021*U_sb_10;
    const ccomplex_t IT_0028 = IT_0020*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0031 = IT_0014*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0034 = IT_0009*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0029 + -IT_0032 +
       0.333333333333333*IT_0035);
    const ccomplex_t IT_0037 = 3*IT_0036;
    const ccomplex_t IT_0038 = 0.166666666666667*IT_0037;
    const ccomplex_t IT_0039 = IT_0026*IT_0038;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0043 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0044 = IT_0009*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0047 = IT_0014*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = m_b*N_d4*e_em*IT_0021*conj(U_sb_11);
    const ccomplex_t IT_0050 = IT_0020*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0045 + (-3)*IT_0048 + 3
      *IT_0051);
    const ccomplex_t IT_0053 = 0.166666666666667*IT_0052;
    const ccomplex_t IT_0054 = m_b*conj(N_d4)*e_em*IT_0021*U_sb_11;
    const ccomplex_t IT_0055 = IT_0020*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0058 = IT_0014*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0061 = IT_0009*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0056 + -IT_0059 +
       0.333333333333333*IT_0062);
    const ccomplex_t IT_0064 = 3*IT_0063;
    const ccomplex_t IT_0065 = 0.166666666666667*IT_0064;
    const ccomplex_t IT_0066 = IT_0053*IT_0065;
    const ccomplex_t IT_0067 = IT_0042*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0041 + -IT_0068;
    const ccomplex_t IT_0070 = cos(alpha);
    const ccomplex_t IT_0071 = IT_0009*IT_0070;
    const ccomplex_t IT_0072 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = sin(alpha);
    const ccomplex_t IT_0075 = IT_0009*IT_0074;
    const ccomplex_t IT_0076 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = IT_0014*IT_0070;
    const ccomplex_t IT_0079 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = IT_0014*IT_0074;
    const ccomplex_t IT_0082 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0073 + -IT_0077 + 
      -IT_0080 + IT_0083);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_b*e_em*IT_0014*IT_0019
      *IT_0021*IT_0070;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_12 + (-2)*IT_0006 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0084*IT_0089;
    const ccomplex_t IT_0091 = IT_0071*IT_0076;
    const ccomplex_t IT_0092 = IT_0072*IT_0075;
    const ccomplex_t IT_0093 = IT_0078*IT_0082;
    const ccomplex_t IT_0094 = IT_0079*IT_0081;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + IT_0092 + 
      -IT_0093 + -IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_b*e_em*IT_0014*IT_0019
      *IT_0021*IT_0074;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0096*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0090 + (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = 3*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0090 + 0.5*IT_0102;
    const ccomplex_t IT_0106 = (-3)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = 3*conj(IT_0103);
    const ccomplex_t IT_0109 = (-3)*conj(IT_0105);
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0112 = IT_0009*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = m_b*conj(N_d4)*e_em*IT_0021*conj(U_sb_00);
    const ccomplex_t IT_0115 = IT_0020*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*(IT_0113 + 1.5*IT_0116);
    const ccomplex_t IT_0118 = (-0.333333333333333)*IT_0117;
    const ccomplex_t IT_0119 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0120 = IT_0009*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = m_b*N_d4*e_em*IT_0021*U_sb_00;
    const ccomplex_t IT_0123 = IT_0020*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 1.5*IT_0124);
    const ccomplex_t IT_0126 = (-0.333333333333333)*IT_0125;
    const ccomplex_t IT_0127 = IT_0118*IT_0126;
    const ccomplex_t IT_0128 = IT_0007*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0131 = IT_0009*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = m_b*N_d4*e_em*IT_0021*U_sb_01;
    const ccomplex_t IT_0134 = IT_0020*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 1.5*IT_0135);
    const ccomplex_t IT_0137 = (-0.333333333333333)*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0139 = IT_0009*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = m_b*conj(N_d4)*e_em*IT_0021*conj(U_sb_01);
    const ccomplex_t IT_0142 = IT_0020*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + 1.5*IT_0143);
    const ccomplex_t IT_0145 = (-0.333333333333333)*IT_0144;
    const ccomplex_t IT_0146 = IT_0137*IT_0145;
    const ccomplex_t IT_0147 = IT_0042*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = -IT_0129 + -IT_0148;
    const ccomplex_t IT_0150 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0151 = IT_0071*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0153 = IT_0075*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0155 = IT_0078*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0157 = IT_0081*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0151 + -IT_0153 + 
      -IT_0155 + IT_0157);
    const ccomplex_t IT_0159 = IT_0089*IT_0158;
    const ccomplex_t IT_0160 = IT_0071*IT_0152;
    const ccomplex_t IT_0161 = IT_0075*IT_0150;
    const ccomplex_t IT_0162 = IT_0078*IT_0156;
    const ccomplex_t IT_0163 = IT_0081*IT_0154;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + IT_0161 + 
      -IT_0162 + -IT_0163);
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = IT_0101*IT_0165;
    const ccomplex_t IT_0167 = (-0.5)*IT_0159 + (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = 3*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0159 + 0.5*IT_0166;
    const ccomplex_t IT_0170 = (-3)*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = 3*conj(IT_0167);
    const ccomplex_t IT_0173 = (-3)*conj(IT_0169);
    const ccomplex_t IT_0174 = IT_0172 + IT_0173;
    const ccomplex_t IT_0175 = -s_24;
    const ccomplex_t IT_0176 = s_23 + IT_0175;
    const ccomplex_t IT_0177 = IT_0000*IT_0176;
    const ccomplex_t IT_0178 = (-2)*IT_0177;
    const ccomplex_t IT_0179 = 4*IT_0003;
    const ccomplex_t IT_0180 = IT_0009*IT_0013;
    const ccomplex_t IT_0181 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0182 = IT_0180*IT_0181;
    const ccomplex_t IT_0183 = IT_0008*IT_0014;
    const ccomplex_t IT_0184 = IT_0181*IT_0183;
    const ccomplex_t IT_0185 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0186 = IT_0180*IT_0185;
    const ccomplex_t IT_0187 = IT_0183*IT_0185;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0182 + IT_0184 + 
      -IT_0186 + -IT_0187);
    const ccomplex_t IT_0189 = 0.5*IT_0188;
    const ccomplex_t IT_0190 = e_em*IT_0180;
    const ccomplex_t IT_0191 = e_em*IT_0183;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0190 + 3*IT_0191);
    const ccomplex_t IT_0193 = (-0.166666666666667)*IT_0192;
    const ccomplex_t IT_0194 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0195 = IT_0193*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = IT_0189*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0013;
    const ccomplex_t IT_0200 = 0.333333333333333*IT_0199;
    const ccomplex_t IT_0201 = IT_0194*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0189*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = (-0.5)*IT_0188;
    const ccomplex_t IT_0206 = IT_0196*IT_0205;
    const ccomplex_t IT_0207 = 0.5*IT_0206;
    const ccomplex_t IT_0208 = IT_0202*IT_0205;
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = (-4)*IT_0177;
    const ccomplex_t IT_0211 = (-2)*IT_0003;
    const ccomplex_t IT_0212 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0213 = IT_0127*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0146*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0214 + IT_0217;
    const ccomplex_t IT_0219 = 3*IT_0105;
    const ccomplex_t IT_0220 = (-3)*IT_0103;
    const ccomplex_t IT_0221 = IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = 3*conj(IT_0105);
    const ccomplex_t IT_0223 = (-3)*conj(IT_0103);
    const ccomplex_t IT_0224 = IT_0222 + IT_0223;
    const ccomplex_t IT_0225 = IT_0066*IT_0215;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0039*IT_0212;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = IT_0226 + IT_0228;
    const ccomplex_t IT_0230 = 3*IT_0169;
    const ccomplex_t IT_0231 = (-3)*IT_0167;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = 3*conj(IT_0169);
    const ccomplex_t IT_0234 = (-3)*conj(IT_0167);
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = 2*IT_0177;
    const ccomplex_t IT_0237 = 4*IT_0177;
    const ccomplex_t IT_0238 = (-0.5)*IT_0197;
    const ccomplex_t IT_0239 = (-0.5)*IT_0203;
    const ccomplex_t IT_0240 = conj(IT_0238) + conj(IT_0239);
    const ccomplex_t IT_0241 = IT_0238 + IT_0239;
    const ccomplex_t IT_0242 = (-0.5)*IT_0206;
    const ccomplex_t IT_0243 = (-0.5)*IT_0208;
    const ccomplex_t IT_0244 = conj(IT_0242) + conj(IT_0243);
    const ccomplex_t IT_0245 = IT_0242 + IT_0243;
    const ccomplex_t IT_0246 = (-4)*IT_0003;
    const ccomplex_t IT_0247 = s_12*IT_0005;
    const ccomplex_t IT_0248 = s_13*s_24;
    const ccomplex_t IT_0249 = s_12*s_34;
    const ccomplex_t IT_0250 = s_14*s_23;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0247 + IT_0248 + IT_0249 + IT_0251;
    const ccomplex_t IT_0253 = IT_0026*IT_0126;
    const ccomplex_t IT_0254 = IT_0007*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = IT_0053*IT_0137;
    const ccomplex_t IT_0257 = IT_0042*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = IT_0255 + IT_0258;
    const ccomplex_t IT_0260 = IT_0009*IT_0018;
    const ccomplex_t IT_0261 = IT_0076*IT_0260;
    const ccomplex_t IT_0262 = sin(beta);
    const ccomplex_t IT_0263 = IT_0009*IT_0262;
    const ccomplex_t IT_0264 = IT_0072*IT_0263;
    const ccomplex_t IT_0265 = IT_0014*IT_0018;
    const ccomplex_t IT_0266 = IT_0082*IT_0265;
    const ccomplex_t IT_0267 = IT_0014*IT_0262;
    const ccomplex_t IT_0268 = IT_0079*IT_0267;
    const ccomplex_t IT_0269 = -IT_0261 + IT_0264 + IT_0266 + -IT_0268;
    const ccomplex_t IT_0270 = m_b*e_em*IT_0014*IT_0019*IT_0021*IT_0262;
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = cpow((-2)*s_12 + (-2)*IT_0006 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0273 = IT_0271*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = IT_0269*IT_0274;
    const ccomplex_t IT_0276 = IT_0072*IT_0260;
    const ccomplex_t IT_0277 = IT_0076*IT_0263;
    const ccomplex_t IT_0278 = IT_0079*IT_0265;
    const ccomplex_t IT_0279 = IT_0082*IT_0267;
    const ccomplex_t IT_0280 = -IT_0276 + -IT_0277 + IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = m_b*e_em*IT_0014*IT_0021;
    const ccomplex_t IT_0282 = 0.5*IT_0281;
    const ccomplex_t IT_0283 = pow(m_Z, 2);
    const ccomplex_t IT_0284 = cpow((-2)*s_12 + (-2)*IT_0006 + IT_0283 + 
      -reg_prop, -1);
    const ccomplex_t IT_0285 = IT_0282*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0280*IT_0286;
    const ccomplex_t IT_0288 = 0.5*IT_0275 + 0.5*IT_0287;
    const ccomplex_t IT_0289 = 6*IT_0288;
    const ccomplex_t IT_0290 = (-0.5)*IT_0275 + (-0.5)*IT_0287;
    const ccomplex_t IT_0291 = (-6)*IT_0290;
    const ccomplex_t IT_0292 = IT_0289 + IT_0291;
    const ccomplex_t IT_0293 = 6*conj(IT_0288);
    const ccomplex_t IT_0294 = (-6)*conj(IT_0290);
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = IT_0065*IT_0145;
    const ccomplex_t IT_0297 = IT_0042*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = IT_0038*IT_0118;
    const ccomplex_t IT_0300 = IT_0007*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = IT_0298 + IT_0301;
    const ccomplex_t IT_0303 = IT_0150*IT_0260;
    const ccomplex_t IT_0304 = IT_0152*IT_0263;
    const ccomplex_t IT_0305 = IT_0154*IT_0265;
    const ccomplex_t IT_0306 = IT_0156*IT_0267;
    const ccomplex_t IT_0307 = IT_0303 + IT_0304 + -IT_0305 + -IT_0306;
    const ccomplex_t IT_0308 = IT_0286*IT_0307;
    const ccomplex_t IT_0309 = IT_0152*IT_0260;
    const ccomplex_t IT_0310 = IT_0150*IT_0263;
    const ccomplex_t IT_0311 = IT_0156*IT_0265;
    const ccomplex_t IT_0312 = IT_0154*IT_0267;
    const ccomplex_t IT_0313 = IT_0309 + -IT_0310 + -IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = IT_0274*IT_0313;
    const ccomplex_t IT_0315 = (-0.5)*IT_0308 + (-0.5)*IT_0314;
    const ccomplex_t IT_0316 = 6*IT_0315;
    const ccomplex_t IT_0317 = 0.5*IT_0308 + 0.5*IT_0314;
    const ccomplex_t IT_0318 = (-6)*IT_0317;
    const ccomplex_t IT_0319 = IT_0316 + IT_0318;
    const ccomplex_t IT_0320 = 6*conj(IT_0315);
    const ccomplex_t IT_0321 = (-6)*conj(IT_0317);
    const ccomplex_t IT_0322 = IT_0320 + IT_0321;
    const ccomplex_t IT_0323 = -IT_0248;
    const ccomplex_t IT_0324 = IT_0247 + IT_0249 + IT_0250 + IT_0323;
    const ccomplex_t IT_0325 = IT_0215*IT_0256;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*IT_0325;
    const ccomplex_t IT_0327 = IT_0212*IT_0253;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = -IT_0326 + -IT_0328;
    const ccomplex_t IT_0330 = 6*IT_0290;
    const ccomplex_t IT_0331 = (-6)*IT_0288;
    const ccomplex_t IT_0332 = IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = 6*conj(IT_0290);
    const ccomplex_t IT_0334 = (-6)*conj(IT_0288);
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = IT_0212*IT_0299;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*IT_0336;
    const ccomplex_t IT_0338 = IT_0215*IT_0296;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*IT_0338;
    const ccomplex_t IT_0340 = -IT_0337 + -IT_0339;
    const ccomplex_t IT_0341 = 6*IT_0317;
    const ccomplex_t IT_0342 = (-6)*IT_0315;
    const ccomplex_t IT_0343 = IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = 6*conj(IT_0317);
    const ccomplex_t IT_0345 = (-6)*conj(IT_0315);
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = s_23 + s_24;
    const ccomplex_t IT_0348 = IT_0000*IT_0347;
    const ccomplex_t IT_0349 = conj(IT_0207) + conj(IT_0238);
    const ccomplex_t IT_0350 = 12*IT_0290;
    const ccomplex_t IT_0351 = (-12)*IT_0288;
    const ccomplex_t IT_0352 = IT_0350 + IT_0351;
    const ccomplex_t IT_0353 = IT_0207 + IT_0238;
    const ccomplex_t IT_0354 = 12*conj(IT_0290);
    const ccomplex_t IT_0355 = (-12)*conj(IT_0288);
    const ccomplex_t IT_0356 = IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = conj(IT_0209) + conj(IT_0239);
    const ccomplex_t IT_0358 = 12*IT_0288;
    const ccomplex_t IT_0359 = (-12)*IT_0290;
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = IT_0209 + IT_0239;
    const ccomplex_t IT_0362 = 12*conj(IT_0288);
    const ccomplex_t IT_0363 = (-12)*conj(IT_0290);
    const ccomplex_t IT_0364 = IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = conj(IT_0198) + conj(IT_0242);
    const ccomplex_t IT_0366 = 12*IT_0315;
    const ccomplex_t IT_0367 = (-12)*IT_0317;
    const ccomplex_t IT_0368 = IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = IT_0198 + IT_0242;
    const ccomplex_t IT_0370 = 12*conj(IT_0315);
    const ccomplex_t IT_0371 = (-12)*conj(IT_0317);
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = conj(IT_0204) + conj(IT_0243);
    const ccomplex_t IT_0374 = 12*IT_0317;
    const ccomplex_t IT_0375 = (-12)*IT_0315;
    const ccomplex_t IT_0376 = IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = IT_0204 + IT_0243;
    const ccomplex_t IT_0378 = 12*conj(IT_0317);
    const ccomplex_t IT_0379 = (-12)*conj(IT_0315);
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = conj(IT_0149)*IT_0292 + IT_0149*IT_0295 + conj
      (IT_0069)*IT_0319 + IT_0069*IT_0322 + conj(IT_0229)*IT_0332 + IT_0229
      *IT_0335 + conj(IT_0218)*IT_0343 + IT_0218*IT_0346 + IT_0349*IT_0352 +
       IT_0353*IT_0356 + IT_0357*IT_0360 + IT_0361*IT_0364 + IT_0365*IT_0368 +
       IT_0369*IT_0372 + IT_0373*IT_0376 + IT_0377*IT_0380;
    const ccomplex_t IT_0382 = -IT_0249;
    const ccomplex_t IT_0383 = IT_0248 + IT_0250 + IT_0382;
    const ccomplex_t IT_0384 = 6*IT_0259;
    const ccomplex_t IT_0385 = 6*conj(IT_0259);
    const ccomplex_t IT_0386 = s_13 + s_14;
    const ccomplex_t IT_0387 = IT_0000*IT_0386;
    const ccomplex_t IT_0388 = 6*IT_0149;
    const ccomplex_t IT_0389 = 6*conj(IT_0149);
    const ccomplex_t IT_0390 = 6*IT_0229;
    const ccomplex_t IT_0391 = 6*conj(IT_0229);
    const ccomplex_t IT_0392 = conj(IT_0069)*IT_0292 + IT_0069*IT_0295 + conj
      (IT_0218)*IT_0332 + IT_0218*IT_0335 + IT_0360*IT_0365 + IT_0357*IT_0368 +
       IT_0364*IT_0369 + IT_0361*IT_0372 + IT_0352*IT_0373 + IT_0349*IT_0376 +
       IT_0356*IT_0377 + IT_0353*IT_0380 + (-6)*conj(IT_0315)*(IT_0229 + (
      -0.166666666666667)*IT_0388) + (-6)*IT_0315*(conj(IT_0229) + (
      -0.166666666666667)*IT_0389) + (-6)*conj(IT_0317)*(IT_0149 + (
      -0.166666666666667)*IT_0390) + (-6)*IT_0317*(conj(IT_0149) + (
      -0.166666666666667)*IT_0391);
    const ccomplex_t IT_0393 = m_b*s_14*m_N_4;
    const ccomplex_t IT_0394 = 12*IT_0229;
    const ccomplex_t IT_0395 = 24*IT_0207;
    const ccomplex_t IT_0396 = 24*IT_0238;
    const ccomplex_t IT_0397 = IT_0394 + IT_0395 + IT_0396;
    const ccomplex_t IT_0398 = 12*conj(IT_0229);
    const ccomplex_t IT_0399 = 24*conj(IT_0207);
    const ccomplex_t IT_0400 = 24*conj(IT_0238);
    const ccomplex_t IT_0401 = IT_0398 + IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = 24*IT_0243;
    const ccomplex_t IT_0403 = 24*IT_0204;
    const ccomplex_t IT_0404 = IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = 24*conj(IT_0243);
    const ccomplex_t IT_0406 = 24*conj(IT_0204);
    const ccomplex_t IT_0407 = IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = 12*IT_0329;
    const ccomplex_t IT_0409 = IT_0384 + IT_0408;
    const ccomplex_t IT_0410 = 12*conj(IT_0329);
    const ccomplex_t IT_0411 = IT_0385 + IT_0410;
    const ccomplex_t IT_0412 = 12*IT_0259;
    const ccomplex_t IT_0413 = 12*conj(IT_0259);
    const ccomplex_t IT_0414 = 12*IT_0302;
    const ccomplex_t IT_0415 = 12*conj(IT_0302);
    const ccomplex_t IT_0416 = conj(IT_0302)*IT_0390 + IT_0302*IT_0391 + conj
      (IT_0340)*IT_0397 + IT_0340*IT_0401 + conj(IT_0329)*IT_0404 + IT_0329
      *IT_0407 + conj(IT_0218)*IT_0409 + IT_0218*IT_0411 + IT_0373*IT_0412 +
       IT_0377*IT_0413 + IT_0349*IT_0414 + IT_0353*IT_0415;
    const ccomplex_t IT_0417 = m_b*s_13*m_N_4;
    const ccomplex_t IT_0418 = 24*IT_0239;
    const ccomplex_t IT_0419 = 24*IT_0209;
    const ccomplex_t IT_0420 = IT_0418 + IT_0419;
    const ccomplex_t IT_0421 = 12*IT_0149;
    const ccomplex_t IT_0422 = IT_0420 + IT_0421;
    const ccomplex_t IT_0423 = 24*conj(IT_0239);
    const ccomplex_t IT_0424 = 24*conj(IT_0209);
    const ccomplex_t IT_0425 = IT_0423 + IT_0424;
    const ccomplex_t IT_0426 = 12*conj(IT_0149);
    const ccomplex_t IT_0427 = IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = 24*IT_0198;
    const ccomplex_t IT_0429 = 24*IT_0242;
    const ccomplex_t IT_0430 = IT_0428 + IT_0429;
    const ccomplex_t IT_0431 = 24*conj(IT_0198);
    const ccomplex_t IT_0432 = 24*conj(IT_0242);
    const ccomplex_t IT_0433 = IT_0431 + IT_0432;
    const ccomplex_t IT_0434 = 6*IT_0329;
    const ccomplex_t IT_0435 = IT_0412 + IT_0434;
    const ccomplex_t IT_0436 = 6*conj(IT_0329);
    const ccomplex_t IT_0437 = IT_0413 + IT_0436;
    const ccomplex_t IT_0438 = 12*IT_0340;
    const ccomplex_t IT_0439 = 12*conj(IT_0340);
    const ccomplex_t IT_0440 = conj(IT_0340)*IT_0388 + IT_0340*IT_0389 +
       IT_0365*IT_0408 + IT_0369*IT_0410 + conj(IT_0302)*IT_0422 + IT_0302
      *IT_0427 + conj(IT_0259)*IT_0430 + IT_0259*IT_0433 + conj(IT_0069)*IT_0435
       + IT_0069*IT_0437 + IT_0357*IT_0438 + IT_0361*IT_0439;
    const ccomplex_t IT_0441 = m_b*s_24*m_N_4;
    const ccomplex_t IT_0442 = 12*IT_0069;
    const ccomplex_t IT_0443 = IT_0430 + IT_0442;
    const ccomplex_t IT_0444 = 6*IT_0069;
    const ccomplex_t IT_0445 = 6*conj(IT_0069);
    const ccomplex_t IT_0446 = 12*conj(IT_0069);
    const ccomplex_t IT_0447 = IT_0357*IT_0408 + IT_0361*IT_0410 + conj
      (IT_0259)*IT_0420 + IT_0259*IT_0425 + conj(IT_0149)*IT_0435 + IT_0149
      *IT_0437 + IT_0365*IT_0438 + IT_0369*IT_0439 + conj(IT_0302)*IT_0443 +
       conj(IT_0340)*IT_0444 + IT_0340*IT_0445 + IT_0302*(IT_0433 + IT_0446);
    const ccomplex_t IT_0448 = m_b*s_23*m_N_4;
    const ccomplex_t IT_0449 = 12*IT_0218;
    const ccomplex_t IT_0450 = 12*conj(IT_0218);
    const ccomplex_t IT_0451 = conj(IT_0218)*IT_0302 + IT_0218*conj(IT_0302) +
       0.166666666666667*conj(IT_0229)*IT_0409 + 0.166666666666667*IT_0229
      *IT_0411 + 4*IT_0349*(IT_0329 + 0.0416666666666667*IT_0412) + 4*IT_0353*
      (conj(IT_0329) + 0.0416666666666667*IT_0413) + 0.166666666666667*IT_0373
      *IT_0414 + 0.166666666666667*IT_0377*IT_0415 + 4*conj(IT_0340)*(IT_0377 +
       0.0416666666666667*IT_0449) + 4*IT_0340*(IT_0373 + 0.0416666666666667
      *IT_0450);
    const ccomplex_t IT_0452 = s_34*IT_0006;
    const ccomplex_t IT_0453 = IT_0404 + IT_0449;
    const ccomplex_t IT_0454 = conj(IT_0218)*IT_0388 + IT_0218*IT_0389 + conj
      (IT_0069)*IT_0390 + IT_0069*IT_0391 + conj(IT_0198)*(IT_0394 + IT_0395) +
       conj(IT_0242)*IT_0397 + IT_0198*(IT_0398 + IT_0399) + IT_0242*IT_0401 +
       IT_0373*IT_0421 + IT_0377*IT_0426 + conj(IT_0238)*IT_0428 + IT_0238
      *IT_0431 + IT_0349*IT_0442 + IT_0353*IT_0446 + IT_0361*(IT_0407 + IT_0450)
       + IT_0357*IT_0453;
    const ccomplex_t IT_0455 = -s_34;
    const ccomplex_t IT_0456 = IT_0005 + IT_0455;
    const ccomplex_t IT_0457 = s_12*IT_0456;
    const ccomplex_t IT_0458 = (-24)*conj(IT_0169);
    const ccomplex_t IT_0459 = 24*conj(IT_0167);
    const ccomplex_t IT_0460 = IT_0458 + IT_0459;
    const ccomplex_t IT_0461 = (-24)*IT_0105;
    const ccomplex_t IT_0462 = s_34 + IT_0005;
    const ccomplex_t IT_0463 = s_12*IT_0462;
    const ccomplex_t IT_0464 = (-24)*IT_0290;
    const ccomplex_t IT_0465 = (-24)*conj(IT_0290);
    const ccomplex_t IT_0466 = IT_0006*IT_0456;
    const ccomplex_t IT_0467 = 6*conj(IT_0103);
    const ccomplex_t IT_0468 = (-6)*conj(IT_0105);
    const ccomplex_t IT_0469 = IT_0467 + IT_0468;
    const ccomplex_t IT_0470 = 6*conj(IT_0105);
    const ccomplex_t IT_0471 = (-6)*conj(IT_0103);
    const ccomplex_t IT_0472 = IT_0470 + IT_0471;
    const ccomplex_t IT_0473 = 6*IT_0167;
    const ccomplex_t IT_0474 = (-6)*IT_0169;
    const ccomplex_t IT_0475 = IT_0473 + IT_0474;
    const ccomplex_t IT_0476 = 6*IT_0169;
    const ccomplex_t IT_0477 = (-6)*IT_0167;
    const ccomplex_t IT_0478 = IT_0476 + IT_0477;
    const ccomplex_t IT_0479 = 6*conj(IT_0169);
    const ccomplex_t IT_0480 = (-6)*conj(IT_0167);
    const ccomplex_t IT_0481 = IT_0479 + IT_0480;
    const ccomplex_t IT_0482 = 6*IT_0105;
    const ccomplex_t IT_0483 = (-6)*IT_0105;
    const ccomplex_t IT_0484 = (-6)*conj(IT_0169);
    const ccomplex_t IT_0485 = conj(IT_0105)*((-24)*IT_0167 + 24*IT_0169) +
       conj(IT_0103)*(24*IT_0167 + (-24)*IT_0169) + 24*IT_0105*conj(IT_0169) + 6
      *IT_0103*(conj(IT_0302) + -conj(IT_0340) + 0.166666666666667*IT_0460) + 6
      *conj(IT_0167)*(IT_0259 + 0.166666666666667*IT_0461) + IT_0302*IT_0469 +
       IT_0340*IT_0472 + conj(IT_0259)*IT_0475 + conj(IT_0329)*IT_0478 + IT_0329
      *IT_0481 + conj(IT_0340)*IT_0482 + conj(IT_0302)*IT_0483 + IT_0259*IT_0484;
    const ccomplex_t IT_0486 = 48*IT_0207;
    const ccomplex_t IT_0487 = 24*IT_0229;
    const ccomplex_t IT_0488 = IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = 48*IT_0204;
    const ccomplex_t IT_0490 = 24*IT_0218;
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = 48*conj(IT_0207);
    const ccomplex_t IT_0493 = 48*conj(IT_0204);
    const ccomplex_t IT_0494 = 48*IT_0243;
    const ccomplex_t IT_0495 = conj(IT_0229)*IT_0397 + conj(IT_0329)*IT_0408 +
       conj(IT_0340)*IT_0438 + conj(IT_0218)*IT_0453 + conj(IT_0238)*(48*IT_0238
       + IT_0486 + IT_0487) + conj(IT_0207)*IT_0488 + conj(IT_0204)*IT_0491 +
       IT_0238*IT_0492 + IT_0243*IT_0493 + conj(IT_0243)*(IT_0491 + IT_0494);
    const ccomplex_t IT_0496 = conj(IT_0229)*IT_0388 + IT_0229*IT_0389 +
       IT_0357*IT_0397 + IT_0361*IT_0401 + IT_0349*IT_0421 + IT_0353*IT_0426 +
       IT_0373*(IT_0428 + IT_0429 + IT_0442) + conj(IT_0218)*IT_0444 + IT_0218
      *IT_0445 + IT_0377*(IT_0431 + IT_0432 + IT_0446) + IT_0365*IT_0449 +
       IT_0369*IT_0450;
    const ccomplex_t IT_0497 = IT_0005*IT_0006;
    const ccomplex_t IT_0498 = 48*IT_0209;
    const ccomplex_t IT_0499 = 24*IT_0149;
    const ccomplex_t IT_0500 = IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = 48*conj(IT_0209);
    const ccomplex_t IT_0502 = 48*conj(IT_0198);
    const ccomplex_t IT_0503 = 48*IT_0198;
    const ccomplex_t IT_0504 = 48*IT_0242;
    const ccomplex_t IT_0505 = 48*IT_0239;
    const ccomplex_t IT_0506 = conj(IT_0218)*IT_0397 + IT_0218*IT_0401 + conj
      (IT_0229)*IT_0404 + conj(IT_0340)*IT_0409 + IT_0340*IT_0411 + conj(IT_0069
      )*IT_0422 + IT_0069*IT_0427 + conj(IT_0149)*IT_0430 + conj(IT_0302)
      *IT_0435 + IT_0302*IT_0437 + conj(IT_0243)*(48*IT_0238 + IT_0486 + IT_0487
      ) + conj(IT_0204)*IT_0488 + IT_0377*IT_0492 + IT_0238*IT_0493 + conj
      (IT_0238)*(IT_0489 + IT_0494) + conj(IT_0198)*IT_0500 + IT_0369*IT_0501 +
       IT_0239*IT_0502 + conj(IT_0239)*(IT_0503 + IT_0504) + conj(IT_0242)*
      (IT_0498 + IT_0499 + IT_0505);
    const ccomplex_t IT_0507 = 24*IT_0069;
    const ccomplex_t IT_0508 = IT_0503 + IT_0507;
    const ccomplex_t IT_0509 = conj(IT_0259)*IT_0412 + conj(IT_0302)*IT_0414 +
       conj(IT_0149)*IT_0422 + conj(IT_0069)*IT_0443 + conj(IT_0209)*IT_0500 +
       IT_0239*IT_0501 + IT_0242*IT_0502 + conj(IT_0239)*(IT_0498 + IT_0499 +
       IT_0505) + conj(IT_0198)*IT_0508 + conj(IT_0242)*(IT_0504 + IT_0508);
    const ccomplex_t IT_0510 = (-0.166666666666667)*IT_0302;
    const ccomplex_t IT_0511 = (-0.166666666666667)*conj(IT_0302);
    const ccomplex_t IT_0512 = IT_0006*IT_0462;
    const ccomplex_t IT_0513 = IT_0335*IT_0340 + IT_0332*conj(IT_0340) + (-6)
      *conj(IT_0317)*(IT_0259 + 4*IT_0288 + (-4)*IT_0290 + (-0.166666666666667)
      *IT_0434) + (-6)*IT_0317*(conj(IT_0259) + 4*conj(IT_0288) + (-4)*conj
      (IT_0290) + (-0.166666666666667)*IT_0436) + 24*conj(IT_0315)*(IT_0288 + 
      -1./4*IT_0329 + 0.0416666666666667*IT_0384 + 0.0416666666666667*IT_0464) +
       24*IT_0315*(conj(IT_0288) + -1./4*conj(IT_0329) + 0.0416666666666667
      *IT_0385 + 0.0416666666666667*IT_0465) + (-6)*IT_0295*IT_0510 + (-6)
      *IT_0292*IT_0511;
    const ccomplex_t IT_0514 = IT_0004*(conj(IT_0069)*IT_0107 + IT_0069
      *IT_0110 + conj(IT_0149)*IT_0171 + IT_0149*IT_0174) + (IT_0110*IT_0149 +
       IT_0107*conj(IT_0149) + conj(IT_0069)*IT_0171 + IT_0069*IT_0174)*IT_0178 
      + IT_0179*(IT_0110*(IT_0198 + IT_0204) + IT_0107*(conj(IT_0198) + conj
      (IT_0204)) + IT_0174*(IT_0207 + IT_0209) + IT_0171*(conj(IT_0207) + conj
      (IT_0209))) + (IT_0174*(IT_0198 + IT_0204) + IT_0171*(conj(IT_0198) + conj
      (IT_0204)) + IT_0110*(IT_0207 + IT_0209) + IT_0107*(conj(IT_0207) + conj
      (IT_0209)))*IT_0210 + IT_0211*(conj(IT_0218)*IT_0221 + IT_0218*IT_0224 +
       conj(IT_0229)*IT_0232 + IT_0229*IT_0235) + (IT_0224*IT_0229 + IT_0221
      *conj(IT_0229) + conj(IT_0218)*IT_0232 + IT_0218*IT_0235)*IT_0236 +
       IT_0237*(IT_0221*IT_0240 + IT_0224*IT_0241 + IT_0232*IT_0244 + IT_0235
      *IT_0245) + (IT_0232*IT_0240 + IT_0235*IT_0241 + IT_0221*IT_0244 + IT_0224
      *IT_0245)*IT_0246 + IT_0252*(conj(IT_0259)*IT_0292 + IT_0259*IT_0295 +
       conj(IT_0302)*IT_0319 + IT_0302*IT_0322) + IT_0324*(conj(IT_0329)*IT_0332
       + IT_0329*IT_0335 + conj(IT_0340)*IT_0343 + IT_0340*IT_0346) + IT_0348
      *IT_0381 + 6*IT_0383*(conj(IT_0302)*IT_0340 + IT_0302*conj(IT_0340) +
       0.166666666666667*conj(IT_0329)*IT_0384 + 0.166666666666667*IT_0329
      *IT_0385) + IT_0387*IT_0392 + IT_0393*IT_0416 + IT_0417*IT_0440 + IT_0441
      *IT_0447 + 6*IT_0448*IT_0451 + IT_0452*IT_0454 + IT_0457*(24*conj(IT_0103)
      *IT_0105 + IT_0103*((-24)*conj(IT_0103) + 24*conj(IT_0105)) + IT_0167*((
      -24)*conj(IT_0167) + 24*conj(IT_0169)) + IT_0169*IT_0460 + conj(IT_0105)
      *IT_0461) + IT_0463*(24*IT_0290*conj(IT_0290) + IT_0317*((-24)*conj
      (IT_0315) + 24*conj(IT_0317)) + IT_0315*(24*conj(IT_0315) + (-24)*conj
      (IT_0317)) + conj(IT_0288)*IT_0464 + IT_0288*(24*conj(IT_0288) + IT_0465))
       + (IT_0247 + IT_0248 + IT_0251 + IT_0382)*(IT_0329*IT_0469 + conj(IT_0340
      )*IT_0475 + 6*conj(IT_0329)*(IT_0103 + 0.166666666666667*IT_0483) + 6
      *IT_0340*(conj(IT_0167) + 0.166666666666667*IT_0484)) + IT_0466*IT_0485 +
       IT_0250*IT_0495 + IT_0247*IT_0496 + IT_0497*IT_0506 + IT_0248*IT_0509 + 
      (IT_0247 + IT_0250 + IT_0323 + IT_0382)*(IT_0259*IT_0472 + (-6)*conj
      (IT_0259)*(IT_0103 + (-0.166666666666667)*IT_0482) + (-6)*IT_0481*IT_0510 
      + (-6)*IT_0478*IT_0511) + IT_0512*IT_0513;
    return create_ccomplex_return(IT_0514);
}

