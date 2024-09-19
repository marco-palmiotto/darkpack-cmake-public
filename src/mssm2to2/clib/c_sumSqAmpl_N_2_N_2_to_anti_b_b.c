#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_anti_b_b(
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
    const creal_t m_N_2 = param->m_N_2;
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = m_b*m_N_2;
    const ccomplex_t IT_0001 = -s_24;
    const ccomplex_t IT_0002 = s_23 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = pow(m_b, 2);
    const ccomplex_t IT_0006 = pow(m_N_2, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = m_b*N_d2*e_em*IT_0013*conj(U_sb_11);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = N_W2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0018 = IT_0011*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = N_B2*e_em*conj(U_sb_01);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0016 + -IT_0019 +
       0.333333333333333*IT_0024);
    const ccomplex_t IT_0026 = 3*IT_0025;
    const ccomplex_t IT_0027 = 0.166666666666667*IT_0026;
    const ccomplex_t IT_0028 = m_b*conj(N_d2)*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0029 = IT_0012*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0032 = IT_0011*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0035 = IT_0021*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0030 + -IT_0033 +
       0.333333333333333*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = (-3)*IT_0038;
    const ccomplex_t IT_0040 = 0.166666666666667*IT_0039;
    const ccomplex_t IT_0041 = IT_0027*IT_0040;
    const ccomplex_t IT_0042 = IT_0007*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0045 = m_b*N_d2*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0046 = IT_0012*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0049 = IT_0011*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0052 = IT_0021*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0047 + -IT_0050 +
       0.333333333333333*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = (-3)*IT_0057;
    const ccomplex_t IT_0059 = 0.166666666666667*IT_0058;
    const ccomplex_t IT_0060 = m_b*conj(N_d2)*e_em*IT_0013*U_sb_10;
    const ccomplex_t IT_0061 = IT_0012*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0064 = IT_0021*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0067 = IT_0011*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0062 + 0.333333333333333
      *IT_0065 + -IT_0068);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = 3*IT_0073;
    const ccomplex_t IT_0075 = 0.333333333333333*IT_0074;
    const ccomplex_t IT_0076 = 3*IT_0075;
    const ccomplex_t IT_0077 = (-0.333333333333333)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = 3*IT_0078;
    const ccomplex_t IT_0080 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0081 = IT_0059*IT_0080;
    const ccomplex_t IT_0082 = IT_0044*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0043 + IT_0083;
    const ccomplex_t IT_0085 = cos(alpha);
    const ccomplex_t IT_0086 = IT_0021*IT_0085;
    const ccomplex_t IT_0087 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = sin(alpha);
    const ccomplex_t IT_0090 = IT_0021*IT_0089;
    const ccomplex_t IT_0091 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = IT_0011*IT_0085;
    const ccomplex_t IT_0094 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = IT_0011*IT_0089;
    const ccomplex_t IT_0097 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0088 + -IT_0092 + 
      -IT_0095 + IT_0098);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_b*e_em*IT_0009*IT_0011
      *IT_0013*IT_0085;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + (-2)*IT_0006 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0099*IT_0104;
    const ccomplex_t IT_0106 = IT_0086*IT_0091;
    const ccomplex_t IT_0107 = IT_0087*IT_0090;
    const ccomplex_t IT_0108 = IT_0093*IT_0097;
    const ccomplex_t IT_0109 = IT_0094*IT_0096;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + IT_0107 + 
      -IT_0108 + -IT_0109);
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*m_b*e_em*IT_0009*IT_0011
      *IT_0013*IT_0089;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0111*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0105 + 0.5*IT_0117;
    const ccomplex_t IT_0119 = 3*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0105 + (-0.5)*IT_0117;
    const ccomplex_t IT_0121 = (-3)*IT_0120;
    const ccomplex_t IT_0122 = IT_0119 + IT_0121;
    const ccomplex_t IT_0123 = 3*conj(IT_0118);
    const ccomplex_t IT_0124 = (-3)*conj(IT_0120);
    const ccomplex_t IT_0125 = IT_0123 + IT_0124;
    const ccomplex_t IT_0126 = m_b*conj(N_d2)*e_em*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0127 = IT_0012*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0130 = IT_0021*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + 0.666666666666667
      *IT_0131);
    const ccomplex_t IT_0133 = 1.5*IT_0132;
    const ccomplex_t IT_0134 = 0.666666666666667*IT_0133;
    const ccomplex_t IT_0135 = 1.5*IT_0134;
    const ccomplex_t IT_0136 = 0.666666666666667*IT_0135;
    const ccomplex_t IT_0137 = 1.5*IT_0136;
    const ccomplex_t IT_0138 = 0.666666666666667*IT_0137;
    const ccomplex_t IT_0139 = 1.5*IT_0138;
    const ccomplex_t IT_0140 = (-0.333333333333333)*IT_0139;
    const ccomplex_t IT_0141 = m_b*N_d2*e_em*IT_0013*U_sb_00;
    const ccomplex_t IT_0142 = IT_0012*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0145 = IT_0021*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + 0.666666666666667
      *IT_0146);
    const ccomplex_t IT_0148 = 1.5*IT_0147;
    const ccomplex_t IT_0149 = 0.666666666666667*IT_0148;
    const ccomplex_t IT_0150 = 1.5*IT_0149;
    const ccomplex_t IT_0151 = 0.666666666666667*IT_0150;
    const ccomplex_t IT_0152 = 1.5*IT_0151;
    const ccomplex_t IT_0153 = 0.666666666666667*IT_0152;
    const ccomplex_t IT_0154 = 1.5*IT_0153;
    const ccomplex_t IT_0155 = 0.666666666666667*IT_0154;
    const ccomplex_t IT_0156 = 1.5*IT_0155;
    const ccomplex_t IT_0157 = 0.666666666666667*IT_0156;
    const ccomplex_t IT_0158 = 1.5*IT_0157;
    const ccomplex_t IT_0159 = (-0.333333333333333)*IT_0158;
    const ccomplex_t IT_0160 = IT_0140*IT_0159;
    const ccomplex_t IT_0161 = IT_0044*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = conj(N_B2)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0164 = IT_0021*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = m_b*conj(N_d2)*e_em*IT_0013*conj(U_sb_01);
    const ccomplex_t IT_0167 = IT_0012*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0165 + 1.5*IT_0168);
    const ccomplex_t IT_0170 = (-0.333333333333333)*IT_0169;
    const ccomplex_t IT_0171 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0172 = IT_0021*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = m_b*N_d2*e_em*IT_0013*U_sb_01;
    const ccomplex_t IT_0175 = IT_0012*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0173 + 1.5*IT_0176);
    const ccomplex_t IT_0178 = (-0.333333333333333)*IT_0177;
    const ccomplex_t IT_0179 = IT_0170*IT_0178;
    const ccomplex_t IT_0180 = IT_0007*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = IT_0162 + IT_0181;
    const ccomplex_t IT_0183 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0184 = IT_0086*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0186 = IT_0090*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0188 = IT_0093*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0190 = IT_0096*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0184 + -IT_0186 + 
      -IT_0188 + IT_0190);
    const ccomplex_t IT_0192 = IT_0104*IT_0191;
    const ccomplex_t IT_0193 = IT_0086*IT_0185;
    const ccomplex_t IT_0194 = IT_0090*IT_0183;
    const ccomplex_t IT_0195 = IT_0093*IT_0189;
    const ccomplex_t IT_0196 = IT_0096*IT_0187;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + IT_0194 + 
      -IT_0195 + -IT_0196);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0116*IT_0198;
    const ccomplex_t IT_0200 = 0.5*IT_0192 + 0.5*IT_0199;
    const ccomplex_t IT_0201 = 3*IT_0200;
    const ccomplex_t IT_0202 = (-0.5)*IT_0192 + (-0.5)*IT_0199;
    const ccomplex_t IT_0203 = (-3)*IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = 3*conj(IT_0200);
    const ccomplex_t IT_0206 = (-3)*conj(IT_0202);
    const ccomplex_t IT_0207 = IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = -s_14;
    const ccomplex_t IT_0209 = s_13 + IT_0208;
    const ccomplex_t IT_0210 = IT_0000*IT_0209;
    const ccomplex_t IT_0211 = (-2)*IT_0210;
    const ccomplex_t IT_0212 = (-4)*IT_0210;
    const ccomplex_t IT_0213 = IT_0010*IT_0021;
    const ccomplex_t IT_0214 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0215 = IT_0213*IT_0214;
    const ccomplex_t IT_0216 = IT_0011*IT_0020;
    const ccomplex_t IT_0217 = IT_0214*IT_0216;
    const ccomplex_t IT_0218 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0219 = IT_0213*IT_0218;
    const ccomplex_t IT_0220 = IT_0216*IT_0218;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0215 + IT_0217 + 
      -IT_0219 + -IT_0220);
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = e_em*IT_0213;
    const ccomplex_t IT_0224 = e_em*IT_0216;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0223 + 3*IT_0224);
    const ccomplex_t IT_0226 = (-0.166666666666667)*IT_0225;
    const ccomplex_t IT_0227 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0228 = IT_0226*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = IT_0222*IT_0229;
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0021;
    const ccomplex_t IT_0233 = 0.333333333333333*IT_0232;
    const ccomplex_t IT_0234 = IT_0227*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = IT_0222*IT_0235;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = conj(IT_0231) + conj(IT_0237);
    const ccomplex_t IT_0239 = IT_0231 + IT_0237;
    const ccomplex_t IT_0240 = 0.5*IT_0221;
    const ccomplex_t IT_0241 = IT_0235*IT_0240;
    const ccomplex_t IT_0242 = (-0.5)*IT_0241;
    const ccomplex_t IT_0243 = IT_0229*IT_0240;
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = 4*IT_0003;
    const ccomplex_t IT_0246 = (-2)*IT_0003;
    const ccomplex_t IT_0247 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0248 = IT_0160*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0251 = IT_0179*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = -IT_0249 + -IT_0252;
    const ccomplex_t IT_0254 = 3*IT_0120;
    const ccomplex_t IT_0255 = (-3)*IT_0118;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 3*conj(IT_0120);
    const ccomplex_t IT_0258 = (-3)*conj(IT_0118);
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = IT_0081*IT_0247;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = IT_0041*IT_0250;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = -IT_0261 + -IT_0263;
    const ccomplex_t IT_0265 = 3*IT_0202;
    const ccomplex_t IT_0266 = (-3)*IT_0200;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = 3*conj(IT_0202);
    const ccomplex_t IT_0269 = (-3)*conj(IT_0200);
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = 2*IT_0210;
    const ccomplex_t IT_0272 = (-4)*IT_0003;
    const ccomplex_t IT_0273 = 0.5*IT_0230;
    const ccomplex_t IT_0274 = 0.5*IT_0236;
    const ccomplex_t IT_0275 = conj(IT_0273) + conj(IT_0274);
    const ccomplex_t IT_0276 = IT_0273 + IT_0274;
    const ccomplex_t IT_0277 = 0.5*IT_0243;
    const ccomplex_t IT_0278 = 0.5*IT_0241;
    const ccomplex_t IT_0279 = conj(IT_0277) + conj(IT_0278);
    const ccomplex_t IT_0280 = IT_0277 + IT_0278;
    const ccomplex_t IT_0281 = 4*IT_0210;
    const ccomplex_t IT_0282 = s_12*IT_0005;
    const ccomplex_t IT_0283 = s_14*s_23;
    const ccomplex_t IT_0284 = s_12*s_34;
    const ccomplex_t IT_0285 = s_13*s_24;
    const ccomplex_t IT_0286 = -IT_0285;
    const ccomplex_t IT_0287 = IT_0282 + IT_0283 + IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = IT_0027*IT_0178;
    const ccomplex_t IT_0289 = IT_0007*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*IT_0289;
    const ccomplex_t IT_0291 = (-3)*IT_0055;
    const ccomplex_t IT_0292 = 0.166666666666667*IT_0291;
    const ccomplex_t IT_0293 = IT_0159*IT_0292;
    const ccomplex_t IT_0294 = IT_0044*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = -IT_0290 + -IT_0295;
    const ccomplex_t IT_0297 = IT_0008*IT_0021;
    const ccomplex_t IT_0298 = IT_0087*IT_0297;
    const ccomplex_t IT_0299 = sin(beta);
    const ccomplex_t IT_0300 = IT_0021*IT_0299;
    const ccomplex_t IT_0301 = IT_0091*IT_0300;
    const ccomplex_t IT_0302 = IT_0008*IT_0011;
    const ccomplex_t IT_0303 = IT_0094*IT_0302;
    const ccomplex_t IT_0304 = IT_0011*IT_0299;
    const ccomplex_t IT_0305 = IT_0097*IT_0304;
    const ccomplex_t IT_0306 = -IT_0298 + -IT_0301 + IT_0303 + IT_0305;
    const ccomplex_t IT_0307 = m_b*e_em*IT_0011*IT_0013;
    const ccomplex_t IT_0308 = 0.5*IT_0307;
    const ccomplex_t IT_0309 = pow(m_Z, 2);
    const ccomplex_t IT_0310 = cpow((-2)*s_12 + (-2)*IT_0006 + IT_0309 + 
      -reg_prop, -1);
    const ccomplex_t IT_0311 = IT_0308*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*IT_0311;
    const ccomplex_t IT_0313 = IT_0306*IT_0312;
    const ccomplex_t IT_0314 = IT_0091*IT_0297;
    const ccomplex_t IT_0315 = IT_0087*IT_0300;
    const ccomplex_t IT_0316 = IT_0097*IT_0302;
    const ccomplex_t IT_0317 = IT_0094*IT_0304;
    const ccomplex_t IT_0318 = -IT_0314 + IT_0315 + IT_0316 + -IT_0317;
    const ccomplex_t IT_0319 = m_b*e_em*IT_0009*IT_0011*IT_0013*IT_0299;
    const ccomplex_t IT_0320 = (-0.5)*IT_0319;
    const ccomplex_t IT_0321 = cpow((-2)*s_12 + (-2)*IT_0006 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0322 = IT_0320*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*IT_0322;
    const ccomplex_t IT_0324 = IT_0318*IT_0323;
    const ccomplex_t IT_0325 = (-0.5)*IT_0313 + (-0.5)*IT_0324;
    const ccomplex_t IT_0326 = 6*IT_0325;
    const ccomplex_t IT_0327 = 0.5*IT_0313 + 0.5*IT_0324;
    const ccomplex_t IT_0328 = (-6)*IT_0327;
    const ccomplex_t IT_0329 = IT_0326 + IT_0328;
    const ccomplex_t IT_0330 = 6*conj(IT_0325);
    const ccomplex_t IT_0331 = (-6)*conj(IT_0327);
    const ccomplex_t IT_0332 = IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = IT_0080*IT_0140;
    const ccomplex_t IT_0334 = IT_0044*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0040*IT_0170;
    const ccomplex_t IT_0337 = IT_0007*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = -IT_0335 + -IT_0338;
    const ccomplex_t IT_0340 = IT_0183*IT_0297;
    const ccomplex_t IT_0341 = IT_0185*IT_0300;
    const ccomplex_t IT_0342 = IT_0187*IT_0302;
    const ccomplex_t IT_0343 = IT_0189*IT_0304;
    const ccomplex_t IT_0344 = IT_0340 + IT_0341 + -IT_0342 + -IT_0343;
    const ccomplex_t IT_0345 = IT_0312*IT_0344;
    const ccomplex_t IT_0346 = IT_0185*IT_0297;
    const ccomplex_t IT_0347 = IT_0183*IT_0300;
    const ccomplex_t IT_0348 = IT_0189*IT_0302;
    const ccomplex_t IT_0349 = IT_0187*IT_0304;
    const ccomplex_t IT_0350 = IT_0346 + -IT_0347 + -IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = IT_0323*IT_0350;
    const ccomplex_t IT_0352 = 0.5*IT_0345 + 0.5*IT_0351;
    const ccomplex_t IT_0353 = 6*IT_0352;
    const ccomplex_t IT_0354 = (-0.5)*IT_0345 + (-0.5)*IT_0351;
    const ccomplex_t IT_0355 = (-6)*IT_0354;
    const ccomplex_t IT_0356 = IT_0353 + IT_0355;
    const ccomplex_t IT_0357 = 6*conj(IT_0352);
    const ccomplex_t IT_0358 = (-6)*conj(IT_0354);
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = -IT_0283;
    const ccomplex_t IT_0361 = IT_0282 + IT_0284 + IT_0285 + IT_0360;
    const ccomplex_t IT_0362 = IT_0059*IT_0159;
    const ccomplex_t IT_0363 = IT_0247*IT_0362;
    const ccomplex_t IT_0364 = (0 + _Complex_I*1)*IT_0363;
    const ccomplex_t IT_0365 = IT_0250*IT_0288;
    const ccomplex_t IT_0366 = (0 + _Complex_I*1)*IT_0365;
    const ccomplex_t IT_0367 = IT_0364 + IT_0366;
    const ccomplex_t IT_0368 = 6*IT_0327;
    const ccomplex_t IT_0369 = (-6)*IT_0325;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 6*conj(IT_0327);
    const ccomplex_t IT_0372 = (-6)*conj(IT_0325);
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = IT_0247*IT_0333;
    const ccomplex_t IT_0375 = (0 + _Complex_I*1)*IT_0374;
    const ccomplex_t IT_0376 = IT_0250*IT_0336;
    const ccomplex_t IT_0377 = (0 + _Complex_I*1)*IT_0376;
    const ccomplex_t IT_0378 = IT_0375 + IT_0377;
    const ccomplex_t IT_0379 = 6*IT_0354;
    const ccomplex_t IT_0380 = (-6)*IT_0352;
    const ccomplex_t IT_0381 = IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = 6*conj(IT_0354);
    const ccomplex_t IT_0383 = (-6)*conj(IT_0352);
    const ccomplex_t IT_0384 = IT_0382 + IT_0383;
    const ccomplex_t IT_0385 = s_13 + s_14;
    const ccomplex_t IT_0386 = IT_0000*IT_0385;
    const ccomplex_t IT_0387 = conj(IT_0231) + conj(IT_0277);
    const ccomplex_t IT_0388 = 12*IT_0327;
    const ccomplex_t IT_0389 = (-12)*IT_0325;
    const ccomplex_t IT_0390 = IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = IT_0231 + IT_0277;
    const ccomplex_t IT_0392 = 12*conj(IT_0327);
    const ccomplex_t IT_0393 = (-12)*conj(IT_0325);
    const ccomplex_t IT_0394 = IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0237) + conj(IT_0278);
    const ccomplex_t IT_0396 = 12*IT_0325;
    const ccomplex_t IT_0397 = (-12)*IT_0327;
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = IT_0237 + IT_0278;
    const ccomplex_t IT_0400 = 12*conj(IT_0325);
    const ccomplex_t IT_0401 = (-12)*conj(IT_0327);
    const ccomplex_t IT_0402 = IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = conj(IT_0244) + conj(IT_0273);
    const ccomplex_t IT_0404 = 12*IT_0352;
    const ccomplex_t IT_0405 = (-12)*IT_0354;
    const ccomplex_t IT_0406 = IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = IT_0244 + IT_0273;
    const ccomplex_t IT_0408 = 12*conj(IT_0352);
    const ccomplex_t IT_0409 = (-12)*conj(IT_0354);
    const ccomplex_t IT_0410 = IT_0408 + IT_0409;
    const ccomplex_t IT_0411 = conj(IT_0242) + conj(IT_0274);
    const ccomplex_t IT_0412 = 12*IT_0354;
    const ccomplex_t IT_0413 = (-12)*IT_0352;
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = IT_0242 + IT_0274;
    const ccomplex_t IT_0416 = 12*conj(IT_0354);
    const ccomplex_t IT_0417 = (-12)*conj(IT_0352);
    const ccomplex_t IT_0418 = IT_0416 + IT_0417;
    const ccomplex_t IT_0419 = conj(IT_0182)*IT_0329 + IT_0182*IT_0332 + conj
      (IT_0084)*IT_0356 + IT_0084*IT_0359 + conj(IT_0264)*IT_0370 + IT_0264
      *IT_0373 + conj(IT_0253)*IT_0381 + IT_0253*IT_0384 + IT_0387*IT_0390 +
       IT_0391*IT_0394 + IT_0395*IT_0398 + IT_0399*IT_0402 + IT_0403*IT_0406 +
       IT_0407*IT_0410 + IT_0411*IT_0414 + IT_0415*IT_0418;
    const ccomplex_t IT_0420 = -IT_0284;
    const ccomplex_t IT_0421 = IT_0283 + IT_0285 + IT_0420;
    const ccomplex_t IT_0422 = 6*IT_0296;
    const ccomplex_t IT_0423 = 6*conj(IT_0296);
    const ccomplex_t IT_0424 = s_23 + s_24;
    const ccomplex_t IT_0425 = IT_0000*IT_0424;
    const ccomplex_t IT_0426 = 6*IT_0182;
    const ccomplex_t IT_0427 = 6*conj(IT_0182);
    const ccomplex_t IT_0428 = 6*IT_0264;
    const ccomplex_t IT_0429 = 6*conj(IT_0264);
    const ccomplex_t IT_0430 = conj(IT_0084)*IT_0329 + IT_0084*IT_0332 + conj
      (IT_0253)*IT_0370 + IT_0253*IT_0373 + IT_0398*IT_0403 + IT_0395*IT_0406 +
       IT_0402*IT_0407 + IT_0399*IT_0410 + IT_0390*IT_0411 + IT_0387*IT_0414 +
       IT_0394*IT_0415 + IT_0391*IT_0418 + (-6)*conj(IT_0352)*(IT_0264 + (
      -0.166666666666667)*IT_0426) + (-6)*IT_0352*(conj(IT_0264) + (
      -0.166666666666667)*IT_0427) + (-6)*conj(IT_0354)*(IT_0182 + (
      -0.166666666666667)*IT_0428) + (-6)*IT_0354*(conj(IT_0182) + (
      -0.166666666666667)*IT_0429);
    const ccomplex_t IT_0431 = m_b*s_24*m_N_2;
    const ccomplex_t IT_0432 = 12*IT_0264;
    const ccomplex_t IT_0433 = 24*IT_0277;
    const ccomplex_t IT_0434 = 24*IT_0231;
    const ccomplex_t IT_0435 = IT_0432 + IT_0433 + IT_0434;
    const ccomplex_t IT_0436 = 12*conj(IT_0264);
    const ccomplex_t IT_0437 = 24*conj(IT_0231);
    const ccomplex_t IT_0438 = 24*conj(IT_0277);
    const ccomplex_t IT_0439 = IT_0436 + IT_0437 + IT_0438;
    const ccomplex_t IT_0440 = 12*IT_0296;
    const ccomplex_t IT_0441 = 24*IT_0367;
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = 12*conj(IT_0296);
    const ccomplex_t IT_0444 = 24*conj(IT_0367);
    const ccomplex_t IT_0445 = IT_0443 + IT_0444;
    const ccomplex_t IT_0446 = 12*IT_0367;
    const ccomplex_t IT_0447 = IT_0422 + IT_0446;
    const ccomplex_t IT_0448 = 12*conj(IT_0367);
    const ccomplex_t IT_0449 = IT_0423 + IT_0448;
    const ccomplex_t IT_0450 = 12*IT_0339;
    const ccomplex_t IT_0451 = 12*conj(IT_0339);
    const ccomplex_t IT_0452 = IT_0431*(conj(IT_0339)*IT_0428 + IT_0339
      *IT_0429 + conj(IT_0378)*IT_0435 + IT_0378*IT_0439 + IT_0411*IT_0442 +
       IT_0415*IT_0445 + conj(IT_0253)*IT_0447 + IT_0253*IT_0449 + IT_0387
      *IT_0450 + IT_0391*IT_0451);
    const ccomplex_t IT_0453 = m_b*s_23*m_N_2;
    const ccomplex_t IT_0454 = 12*IT_0182;
    const ccomplex_t IT_0455 = 24*IT_0237;
    const ccomplex_t IT_0456 = 24*IT_0278;
    const ccomplex_t IT_0457 = IT_0454 + IT_0455 + IT_0456;
    const ccomplex_t IT_0458 = 12*conj(IT_0182);
    const ccomplex_t IT_0459 = 24*conj(IT_0237);
    const ccomplex_t IT_0460 = 24*conj(IT_0278);
    const ccomplex_t IT_0461 = IT_0458 + IT_0459 + IT_0460;
    const ccomplex_t IT_0462 = 24*IT_0244;
    const ccomplex_t IT_0463 = 24*IT_0273;
    const ccomplex_t IT_0464 = IT_0462 + IT_0463;
    const ccomplex_t IT_0465 = 24*conj(IT_0244);
    const ccomplex_t IT_0466 = 24*conj(IT_0273);
    const ccomplex_t IT_0467 = IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = 6*IT_0367;
    const ccomplex_t IT_0469 = IT_0440 + IT_0468;
    const ccomplex_t IT_0470 = 6*conj(IT_0367);
    const ccomplex_t IT_0471 = IT_0443 + IT_0470;
    const ccomplex_t IT_0472 = 12*IT_0378;
    const ccomplex_t IT_0473 = 12*conj(IT_0378);
    const ccomplex_t IT_0474 = conj(IT_0378)*IT_0426 + IT_0378*IT_0427 +
       IT_0403*IT_0446 + IT_0407*IT_0448 + conj(IT_0339)*IT_0457 + IT_0339
      *IT_0461 + conj(IT_0296)*IT_0464 + IT_0296*IT_0467 + conj(IT_0084)*IT_0469
       + IT_0084*IT_0471 + IT_0395*IT_0472 + IT_0399*IT_0473;
    const ccomplex_t IT_0475 = m_b*s_14*m_N_2;
    const ccomplex_t IT_0476 = 12*IT_0084;
    const ccomplex_t IT_0477 = IT_0464 + IT_0476;
    const ccomplex_t IT_0478 = 12*conj(IT_0084);
    const ccomplex_t IT_0479 = IT_0467 + IT_0478;
    const ccomplex_t IT_0480 = 6*IT_0084;
    const ccomplex_t IT_0481 = 6*conj(IT_0084);
    const ccomplex_t IT_0482 = IT_0395*(24*IT_0296 + IT_0446) + IT_0399*(24
      *conj(IT_0296) + IT_0448) + conj(IT_0182)*IT_0469 + IT_0182*IT_0471 +
       IT_0403*IT_0472 + IT_0407*IT_0473 + conj(IT_0339)*IT_0477 + IT_0339
      *IT_0479 + conj(IT_0378)*IT_0480 + IT_0378*IT_0481;
    const ccomplex_t IT_0483 = s_34*IT_0006;
    const ccomplex_t IT_0484 = IT_0462 + IT_0463 + IT_0476;
    const ccomplex_t IT_0485 = IT_0465 + IT_0466 + IT_0478;
    const ccomplex_t IT_0486 = IT_0436 + IT_0437;
    const ccomplex_t IT_0487 = 12*IT_0253;
    const ccomplex_t IT_0488 = 12*conj(IT_0253);
    const ccomplex_t IT_0489 = conj(IT_0253)*IT_0426 + IT_0253*IT_0427 + conj
      (IT_0084)*IT_0428 + IT_0084*IT_0429 + IT_0403*IT_0432 + IT_0411*IT_0457 +
       IT_0415*IT_0461 + conj(IT_0231)*IT_0476 + IT_0231*IT_0479 + conj(IT_0277)
      *IT_0484 + IT_0277*IT_0485 + IT_0407*IT_0486 + IT_0395*IT_0487 + IT_0399
      *IT_0488;
    const ccomplex_t IT_0490 = conj(IT_0264)*IT_0426 + IT_0264*IT_0427 +
       IT_0395*IT_0432 + conj(IT_0231)*IT_0454 + conj(IT_0277)*IT_0457 + IT_0391
      *IT_0461 + conj(IT_0253)*IT_0480 + IT_0253*IT_0481 + IT_0411*IT_0484 +
       IT_0415*IT_0485 + IT_0399*IT_0486 + IT_0403*IT_0487 + IT_0407*IT_0488;
    const ccomplex_t IT_0491 = 48*IT_0244;
    const ccomplex_t IT_0492 = 24*IT_0084;
    const ccomplex_t IT_0493 = IT_0491 + IT_0492;
    const ccomplex_t IT_0494 = 48*IT_0273;
    const ccomplex_t IT_0495 = IT_0493 + IT_0494;
    const ccomplex_t IT_0496 = IT_0454 + IT_0455;
    const ccomplex_t IT_0497 = 48*IT_0237;
    const ccomplex_t IT_0498 = 48*conj(IT_0237);
    const ccomplex_t IT_0499 = 48*conj(IT_0244);
    const ccomplex_t IT_0500 = conj(IT_0296)*IT_0440 + conj(IT_0339)*IT_0450 +
       conj(IT_0084)*IT_0477 + conj(IT_0244)*IT_0493 + conj(IT_0273)*IT_0495 +
       conj(IT_0182)*IT_0496 + IT_0395*(24*IT_0182 + IT_0497) + IT_0278*(24*conj
      (IT_0182) + 48*conj(IT_0278) + IT_0498) + IT_0273*IT_0499;
    const ccomplex_t IT_0501 = -s_34;
    const ccomplex_t IT_0502 = IT_0005 + IT_0501;
    const ccomplex_t IT_0503 = s_12*IT_0502;
    const ccomplex_t IT_0504 = (-24)*conj(IT_0202);
    const ccomplex_t IT_0505 = 24*conj(IT_0200);
    const ccomplex_t IT_0506 = IT_0504 + IT_0505;
    const ccomplex_t IT_0507 = (-24)*IT_0120;
    const ccomplex_t IT_0508 = s_34 + IT_0005;
    const ccomplex_t IT_0509 = s_12*IT_0508;
    const ccomplex_t IT_0510 = (-24)*IT_0327;
    const ccomplex_t IT_0511 = (-24)*conj(IT_0327);
    const ccomplex_t IT_0512 = 6*conj(IT_0118);
    const ccomplex_t IT_0513 = (-6)*conj(IT_0120);
    const ccomplex_t IT_0514 = IT_0512 + IT_0513;
    const ccomplex_t IT_0515 = 6*IT_0200;
    const ccomplex_t IT_0516 = (-6)*IT_0202;
    const ccomplex_t IT_0517 = IT_0515 + IT_0516;
    const ccomplex_t IT_0518 = (-6)*IT_0120;
    const ccomplex_t IT_0519 = (-6)*conj(IT_0202);
    const ccomplex_t IT_0520 = m_b*s_13*m_N_2;
    const ccomplex_t IT_0521 = 24*IT_0274;
    const ccomplex_t IT_0522 = 24*IT_0242;
    const ccomplex_t IT_0523 = IT_0487 + IT_0521 + IT_0522;
    const ccomplex_t IT_0524 = 6*conj(IT_0339);
    const ccomplex_t IT_0525 = conj(IT_0253)*IT_0339 + 0.166666666666667
      *IT_0387*IT_0442 + 0.166666666666667*IT_0391*IT_0445 + 0.166666666666667
      *conj(IT_0264)*IT_0447 + 0.166666666666667*IT_0264*IT_0449 +
       0.166666666666667*IT_0411*IT_0450 + 0.166666666666667*IT_0415*IT_0451 +
       0.166666666666667*IT_0378*(24*conj(IT_0242) + 24*conj(IT_0274) + IT_0488)
       + 0.166666666666667*conj(IT_0378)*IT_0523 + 0.166666666666667*IT_0253
      *IT_0524;
    const ccomplex_t IT_0526 = IT_0006*IT_0502;
    const ccomplex_t IT_0527 = 6*conj(IT_0120);
    const ccomplex_t IT_0528 = (-6)*conj(IT_0118);
    const ccomplex_t IT_0529 = IT_0527 + IT_0528;
    const ccomplex_t IT_0530 = 6*IT_0202;
    const ccomplex_t IT_0531 = (-6)*IT_0200;
    const ccomplex_t IT_0532 = IT_0530 + IT_0531;
    const ccomplex_t IT_0533 = 6*conj(IT_0202);
    const ccomplex_t IT_0534 = (-6)*conj(IT_0200);
    const ccomplex_t IT_0535 = IT_0533 + IT_0534;
    const ccomplex_t IT_0536 = 6*IT_0120;
    const ccomplex_t IT_0537 = conj(IT_0120)*((-24)*IT_0200 + 24*IT_0202) +
       conj(IT_0118)*(24*IT_0200 + (-24)*IT_0202) + 24*IT_0120*conj(IT_0202) + 6
      *conj(IT_0200)*(IT_0296 + 0.166666666666667*IT_0507) + IT_0339*IT_0514 +
       conj(IT_0296)*IT_0517 + conj(IT_0339)*IT_0518 + IT_0296*IT_0519 + (-6)
      *IT_0118*(conj(IT_0378) + (-0.166666666666667)*IT_0506 + (
      -0.166666666666667)*IT_0524) + IT_0378*IT_0529 + conj(IT_0367)*IT_0532 +
       IT_0367*IT_0535 + conj(IT_0378)*IT_0536;
    const ccomplex_t IT_0538 = IT_0005*IT_0006;
    const ccomplex_t IT_0539 = 48*IT_0231;
    const ccomplex_t IT_0540 = 24*IT_0264;
    const ccomplex_t IT_0541 = IT_0539 + IT_0540;
    const ccomplex_t IT_0542 = 48*conj(IT_0242);
    const ccomplex_t IT_0543 = 48*IT_0242;
    const ccomplex_t IT_0544 = 48*IT_0274;
    const ccomplex_t IT_0545 = 48*conj(IT_0231);
    const ccomplex_t IT_0546 = 24*conj(IT_0264);
    const ccomplex_t IT_0547 = conj(IT_0253)*IT_0435 + IT_0253*IT_0439 + conj
      (IT_0378)*IT_0447 + IT_0378*IT_0449 + IT_0084*(IT_0458 + IT_0459) + conj
      (IT_0182)*IT_0464 + IT_0182*IT_0467 + conj(IT_0339)*IT_0469 + IT_0339
      *IT_0471 + conj(IT_0278)*IT_0495 + conj(IT_0084)*IT_0496 + IT_0403*IT_0497
       + IT_0407*IT_0498 + IT_0278*(24*conj(IT_0084) + 48*conj(IT_0273) +
       IT_0499) + conj(IT_0274)*(48*IT_0277 + IT_0539 + IT_0540) + conj(IT_0242)
      *IT_0541 + IT_0277*IT_0542 + conj(IT_0277)*(IT_0543 + IT_0544) + IT_0415*
      (IT_0545 + IT_0546);
    const ccomplex_t IT_0548 = 24*IT_0253;
    const ccomplex_t IT_0549 = IT_0543 + IT_0548;
    const ccomplex_t IT_0550 = conj(IT_0264)*(IT_0432 + IT_0434) + conj
      (IT_0367)*IT_0446 + conj(IT_0378)*IT_0472 + conj(IT_0253)*IT_0523 +
       IT_0387*IT_0541 + IT_0274*IT_0542 + IT_0277*(48*conj(IT_0277) + IT_0545 +
       IT_0546) + conj(IT_0242)*IT_0549 + conj(IT_0274)*(IT_0544 + IT_0549);
    const ccomplex_t IT_0551 = (-0.166666666666667)*IT_0339;
    const ccomplex_t IT_0552 = (-0.166666666666667)*conj(IT_0339);
    const ccomplex_t IT_0553 = IT_0006*IT_0508;
    const ccomplex_t IT_0554 = IT_0373*IT_0378 + IT_0370*conj(IT_0378) + (-6)
      *conj(IT_0354)*(IT_0296 + 4*IT_0325 + (-4)*IT_0327 + (-0.166666666666667)
      *IT_0468) + (-6)*IT_0354*(conj(IT_0296) + 4*conj(IT_0325) + (-4)*conj
      (IT_0327) + (-0.166666666666667)*IT_0470) + 24*conj(IT_0352)*(IT_0325 + 
      -1./4*IT_0367 + 0.0416666666666667*IT_0422 + 0.0416666666666667*IT_0510) +
       24*IT_0352*(conj(IT_0325) + -1./4*conj(IT_0367) + 0.0416666666666667
      *IT_0423 + 0.0416666666666667*IT_0511) + (-6)*IT_0332*IT_0551 + (-6)
      *IT_0329*IT_0552;
    const ccomplex_t IT_0555 = IT_0004*(conj(IT_0084)*IT_0122 + IT_0084
      *IT_0125 + conj(IT_0182)*IT_0204 + IT_0182*IT_0207) + (IT_0125*IT_0182 +
       IT_0122*conj(IT_0182) + conj(IT_0084)*IT_0204 + IT_0084*IT_0207)*IT_0211 
      + IT_0212*(IT_0122*IT_0238 + IT_0125*IT_0239 + IT_0207*(IT_0242 + IT_0244)
       + IT_0204*(conj(IT_0242) + conj(IT_0244))) + (IT_0204*IT_0238 + IT_0207
      *IT_0239 + IT_0125*(IT_0242 + IT_0244) + IT_0122*(conj(IT_0242) + conj
      (IT_0244)))*IT_0245 + IT_0246*(conj(IT_0253)*IT_0256 + IT_0253*IT_0259 +
       conj(IT_0264)*IT_0267 + IT_0264*IT_0270) + (IT_0259*IT_0264 + IT_0256
      *conj(IT_0264) + conj(IT_0253)*IT_0267 + IT_0253*IT_0270)*IT_0271 +
       IT_0272*(IT_0256*IT_0275 + IT_0259*IT_0276 + IT_0267*IT_0279 + IT_0270
      *IT_0280) + (IT_0267*IT_0275 + IT_0270*IT_0276 + IT_0256*IT_0279 + IT_0259
      *IT_0280)*IT_0281 + IT_0287*(conj(IT_0296)*IT_0329 + IT_0296*IT_0332 +
       conj(IT_0339)*IT_0356 + IT_0339*IT_0359) + IT_0361*(conj(IT_0367)*IT_0370
       + IT_0367*IT_0373 + conj(IT_0378)*IT_0381 + IT_0378*IT_0384) + IT_0386
      *IT_0419 + 6*IT_0421*(conj(IT_0339)*IT_0378 + IT_0339*conj(IT_0378) +
       0.166666666666667*conj(IT_0367)*IT_0422 + 0.166666666666667*IT_0367
      *IT_0423) + IT_0425*IT_0430 + IT_0452 + IT_0453*IT_0474 + IT_0475*IT_0482 
      + IT_0483*IT_0489 + IT_0282*IT_0490 + IT_0283*IT_0500 + IT_0503*(24*conj
      (IT_0118)*IT_0120 + IT_0118*((-24)*conj(IT_0118) + 24*conj(IT_0120)) +
       IT_0200*((-24)*conj(IT_0200) + 24*conj(IT_0202)) + IT_0202*IT_0506 + conj
      (IT_0120)*IT_0507) + IT_0509*(24*IT_0327*conj(IT_0327) + IT_0354*((-24)
      *conj(IT_0352) + 24*conj(IT_0354)) + IT_0352*(24*conj(IT_0352) + (-24)
      *conj(IT_0354)) + conj(IT_0325)*IT_0510 + IT_0325*(24*conj(IT_0325) +
       IT_0511)) + (IT_0282 + IT_0283 + IT_0286 + IT_0420)*(IT_0367*IT_0514 +
       conj(IT_0378)*IT_0517 + 6*conj(IT_0367)*(IT_0118 + 0.166666666666667
      *IT_0518) + 6*IT_0378*(conj(IT_0200) + 0.166666666666667*IT_0519)) + 6
      *IT_0520*IT_0525 + IT_0526*IT_0537 + IT_0538*IT_0547 + IT_0285*IT_0550 + 
      (IT_0282 + IT_0285 + IT_0360 + IT_0420)*(IT_0296*IT_0529 + (-6)*conj
      (IT_0296)*(IT_0118 + (-0.166666666666667)*IT_0536) + (-6)*IT_0535*IT_0551 
      + (-6)*IT_0532*IT_0552) + IT_0553*IT_0554;
    return create_ccomplex_return(IT_0555);
}

