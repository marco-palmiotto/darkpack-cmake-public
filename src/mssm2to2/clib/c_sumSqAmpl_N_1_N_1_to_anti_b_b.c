#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_b_b(
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
    const creal_t m_N_1 = param->m_N_1;
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
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = m_b*m_N_1;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = pow(m_b, 2);
    const ccomplex_t IT_0006 = pow(m_N_1, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = m_b*N_d1*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0018 = IT_0011*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0016 + -IT_0019 +
       0.333333333333333*IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = (-3)*IT_0026;
    const ccomplex_t IT_0028 = 0.166666666666667*IT_0027;
    const ccomplex_t IT_0029 = m_b*conj(N_d1)*e_em*IT_0013*U_sb_10;
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0033 = IT_0011*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0036 = IT_0021*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*(IT_0031 + -IT_0034 +
       0.333333333333333*IT_0037);
    const ccomplex_t IT_0039 = -IT_0038;
    const ccomplex_t IT_0040 = (-3)*IT_0039;
    const ccomplex_t IT_0041 = (-0.333333333333333)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = 3*IT_0042;
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = IT_0028*IT_0044;
    const ccomplex_t IT_0046 = IT_0007*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0049 = m_b*conj(N_d1)*e_em*IT_0013*U_sb_11;
    const ccomplex_t IT_0050 = IT_0012*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0053 = IT_0011*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0056 = IT_0021*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0051 + -IT_0054 +
       0.333333333333333*IT_0057);
    const ccomplex_t IT_0059 = 3*IT_0058;
    const ccomplex_t IT_0060 = 0.166666666666667*IT_0059;
    const ccomplex_t IT_0061 = N_W1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0062 = IT_0011*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = m_b*N_d1*e_em*IT_0013*conj(U_sb_11);
    const ccomplex_t IT_0065 = IT_0012*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_B1*e_em*conj(U_sb_01);
    const ccomplex_t IT_0068 = IT_0021*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0063 + -IT_0066 + (
      -0.333333333333333)*IT_0069);
    const ccomplex_t IT_0071 = (-3)*IT_0070;
    const ccomplex_t IT_0072 = 0.166666666666667*IT_0071;
    const ccomplex_t IT_0073 = IT_0060*IT_0072;
    const ccomplex_t IT_0074 = IT_0048*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = -IT_0047 + -IT_0075;
    const ccomplex_t IT_0077 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0078 = cos(alpha);
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0082 = sin(alpha);
    const ccomplex_t IT_0083 = IT_0021*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0086 = IT_0011*IT_0078;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0089 = IT_0011*IT_0082;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0080 + -IT_0084 + 
      -IT_0087 + IT_0090);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*m_b*e_em*IT_0009*IT_0011
      *IT_0013*IT_0078;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_12 + (-2)*IT_0006 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = IT_0091*IT_0096;
    const ccomplex_t IT_0098 = IT_0079*IT_0081;
    const ccomplex_t IT_0099 = IT_0077*IT_0083;
    const ccomplex_t IT_0100 = IT_0086*IT_0088;
    const ccomplex_t IT_0101 = IT_0085*IT_0089;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + IT_0099 + 
      -IT_0100 + -IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_b*e_em*IT_0009*IT_0011
      *IT_0013*IT_0082;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0105*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0103*IT_0108;
    const ccomplex_t IT_0110 = (-0.5)*IT_0097 + (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = 3*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0097 + 0.5*IT_0109;
    const ccomplex_t IT_0113 = (-3)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = 3*conj(IT_0110);
    const ccomplex_t IT_0116 = (-3)*conj(IT_0112);
    const ccomplex_t IT_0117 = IT_0115 + IT_0116;
    const ccomplex_t IT_0118 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0119 = IT_0021*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = m_b*N_d1*e_em*IT_0013*U_sb_00;
    const ccomplex_t IT_0122 = IT_0012*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + 1.5*IT_0123);
    const ccomplex_t IT_0125 = (-0.333333333333333)*IT_0124;
    const ccomplex_t IT_0126 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0127 = IT_0021*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = m_b*conj(N_d1)*e_em*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0130 = IT_0012*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + 1.5*IT_0131);
    const ccomplex_t IT_0133 = (-0.333333333333333)*IT_0132;
    const ccomplex_t IT_0134 = IT_0125*IT_0133;
    const ccomplex_t IT_0135 = IT_0007*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = conj(N_B1)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0138 = IT_0021*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = m_b*conj(N_d1)*e_em*IT_0013*conj(U_sb_01);
    const ccomplex_t IT_0141 = IT_0012*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 1.5*IT_0142);
    const ccomplex_t IT_0144 = (-0.333333333333333)*IT_0143;
    const ccomplex_t IT_0145 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0146 = IT_0021*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_b*N_d1*e_em*IT_0013*U_sb_01;
    const ccomplex_t IT_0149 = IT_0012*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 1.5*IT_0150);
    const ccomplex_t IT_0152 = (-0.333333333333333)*IT_0151;
    const ccomplex_t IT_0153 = IT_0144*IT_0152;
    const ccomplex_t IT_0154 = IT_0048*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = -IT_0136 + -IT_0155;
    const ccomplex_t IT_0157 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0158 = IT_0079*IT_0157;
    const ccomplex_t IT_0159 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0160 = IT_0083*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0162 = IT_0086*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0164 = IT_0089*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0158 + -IT_0160 + 
      -IT_0162 + IT_0164);
    const ccomplex_t IT_0166 = IT_0096*IT_0165;
    const ccomplex_t IT_0167 = IT_0079*IT_0159;
    const ccomplex_t IT_0168 = IT_0083*IT_0157;
    const ccomplex_t IT_0169 = IT_0086*IT_0163;
    const ccomplex_t IT_0170 = IT_0089*IT_0161;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0167 + IT_0168 + 
      -IT_0169 + -IT_0170);
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0108*IT_0172;
    const ccomplex_t IT_0174 = (-0.5)*IT_0166 + (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = 3*IT_0174;
    const ccomplex_t IT_0176 = 0.5*IT_0166 + 0.5*IT_0173;
    const ccomplex_t IT_0177 = (-3)*IT_0176;
    const ccomplex_t IT_0178 = IT_0175 + IT_0177;
    const ccomplex_t IT_0179 = 3*conj(IT_0174);
    const ccomplex_t IT_0180 = (-3)*conj(IT_0176);
    const ccomplex_t IT_0181 = IT_0179 + IT_0180;
    const ccomplex_t IT_0182 = -s_24;
    const ccomplex_t IT_0183 = s_23 + IT_0182;
    const ccomplex_t IT_0184 = IT_0000*IT_0183;
    const ccomplex_t IT_0185 = (-2)*IT_0184;
    const ccomplex_t IT_0186 = 4*IT_0003;
    const ccomplex_t IT_0187 = IT_0010*IT_0021;
    const ccomplex_t IT_0188 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = IT_0011*IT_0020;
    const ccomplex_t IT_0191 = IT_0188*IT_0190;
    const ccomplex_t IT_0192 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0193 = IT_0187*IT_0192;
    const ccomplex_t IT_0194 = IT_0190*IT_0192;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*(IT_0189 + IT_0191 + 
      -IT_0193 + -IT_0194);
    const ccomplex_t IT_0196 = 0.5*IT_0195;
    const ccomplex_t IT_0197 = e_em*IT_0187;
    const ccomplex_t IT_0198 = e_em*IT_0190;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0197 + 3*IT_0198);
    const ccomplex_t IT_0200 = (-0.166666666666667)*IT_0199;
    const ccomplex_t IT_0201 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0202 = IT_0200*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = IT_0196*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*e_em*IT_0010*IT_0021;
    const ccomplex_t IT_0207 = 0.333333333333333*IT_0206;
    const ccomplex_t IT_0208 = IT_0201*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = IT_0196*IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0210;
    const ccomplex_t IT_0212 = (-0.5)*IT_0195;
    const ccomplex_t IT_0213 = IT_0203*IT_0212;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0209*IT_0212;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = (-4)*IT_0184;
    const ccomplex_t IT_0218 = (-2)*IT_0003;
    const ccomplex_t IT_0219 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0220 = IT_0134*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0223 = IT_0153*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0221 + IT_0224;
    const ccomplex_t IT_0226 = 3*IT_0112;
    const ccomplex_t IT_0227 = (-3)*IT_0110;
    const ccomplex_t IT_0228 = IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = 3*conj(IT_0112);
    const ccomplex_t IT_0230 = (-3)*conj(IT_0110);
    const ccomplex_t IT_0231 = IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = IT_0045*IT_0219;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = IT_0073*IT_0222;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = IT_0233 + IT_0235;
    const ccomplex_t IT_0237 = 3*IT_0176;
    const ccomplex_t IT_0238 = (-3)*IT_0174;
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = 3*conj(IT_0176);
    const ccomplex_t IT_0241 = (-3)*conj(IT_0174);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = 2*IT_0184;
    const ccomplex_t IT_0244 = (-4)*IT_0003;
    const ccomplex_t IT_0245 = (-0.5)*IT_0213;
    const ccomplex_t IT_0246 = (-0.5)*IT_0215;
    const ccomplex_t IT_0247 = conj(IT_0245) + conj(IT_0246);
    const ccomplex_t IT_0248 = IT_0245 + IT_0246;
    const ccomplex_t IT_0249 = (-0.5)*IT_0210;
    const ccomplex_t IT_0250 = (-0.5)*IT_0204;
    const ccomplex_t IT_0251 = conj(IT_0249) + conj(IT_0250);
    const ccomplex_t IT_0252 = IT_0249 + IT_0250;
    const ccomplex_t IT_0253 = 4*IT_0184;
    const ccomplex_t IT_0254 = s_12*IT_0005;
    const ccomplex_t IT_0255 = s_13*s_24;
    const ccomplex_t IT_0256 = s_12*s_34;
    const ccomplex_t IT_0257 = s_14*s_23;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = IT_0254 + IT_0255 + IT_0256 + IT_0258;
    const ccomplex_t IT_0260 = IT_0028*IT_0125;
    const ccomplex_t IT_0261 = IT_0007*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = IT_0072*IT_0152;
    const ccomplex_t IT_0264 = IT_0048*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = IT_0262 + IT_0265;
    const ccomplex_t IT_0267 = IT_0008*IT_0021;
    const ccomplex_t IT_0268 = IT_0081*IT_0267;
    const ccomplex_t IT_0269 = sin(beta);
    const ccomplex_t IT_0270 = IT_0021*IT_0269;
    const ccomplex_t IT_0271 = IT_0077*IT_0270;
    const ccomplex_t IT_0272 = IT_0008*IT_0011;
    const ccomplex_t IT_0273 = IT_0088*IT_0272;
    const ccomplex_t IT_0274 = IT_0011*IT_0269;
    const ccomplex_t IT_0275 = IT_0085*IT_0274;
    const ccomplex_t IT_0276 = -IT_0268 + IT_0271 + IT_0273 + -IT_0275;
    const ccomplex_t IT_0277 = m_b*e_em*IT_0009*IT_0011*IT_0013*IT_0269;
    const ccomplex_t IT_0278 = (-0.5)*IT_0277;
    const ccomplex_t IT_0279 = cpow((-2)*s_12 + (-2)*IT_0006 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0280 = IT_0278*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = IT_0276*IT_0281;
    const ccomplex_t IT_0283 = IT_0077*IT_0267;
    const ccomplex_t IT_0284 = IT_0081*IT_0270;
    const ccomplex_t IT_0285 = IT_0085*IT_0272;
    const ccomplex_t IT_0286 = IT_0088*IT_0274;
    const ccomplex_t IT_0287 = -IT_0283 + -IT_0284 + IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = m_b*e_em*IT_0011*IT_0013;
    const ccomplex_t IT_0289 = 0.5*IT_0288;
    const ccomplex_t IT_0290 = pow(m_Z, 2);
    const ccomplex_t IT_0291 = cpow((-2)*s_12 + (-2)*IT_0006 + IT_0290 + 
      -reg_prop, -1);
    const ccomplex_t IT_0292 = IT_0289*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*IT_0292;
    const ccomplex_t IT_0294 = IT_0287*IT_0293;
    const ccomplex_t IT_0295 = 0.5*IT_0282 + 0.5*IT_0294;
    const ccomplex_t IT_0296 = 6*IT_0295;
    const ccomplex_t IT_0297 = (-0.5)*IT_0282 + (-0.5)*IT_0294;
    const ccomplex_t IT_0298 = (-6)*IT_0297;
    const ccomplex_t IT_0299 = IT_0296 + IT_0298;
    const ccomplex_t IT_0300 = 6*conj(IT_0295);
    const ccomplex_t IT_0301 = (-6)*conj(IT_0297);
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0060*IT_0144;
    const ccomplex_t IT_0304 = IT_0048*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*IT_0304;
    const ccomplex_t IT_0306 = IT_0044*IT_0133;
    const ccomplex_t IT_0307 = IT_0007*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = IT_0305 + IT_0308;
    const ccomplex_t IT_0310 = IT_0159*IT_0267;
    const ccomplex_t IT_0311 = IT_0157*IT_0270;
    const ccomplex_t IT_0312 = IT_0163*IT_0272;
    const ccomplex_t IT_0313 = IT_0161*IT_0274;
    const ccomplex_t IT_0314 = IT_0310 + -IT_0311 + -IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = IT_0281*IT_0314;
    const ccomplex_t IT_0316 = IT_0157*IT_0267;
    const ccomplex_t IT_0317 = IT_0159*IT_0270;
    const ccomplex_t IT_0318 = IT_0161*IT_0272;
    const ccomplex_t IT_0319 = IT_0163*IT_0274;
    const ccomplex_t IT_0320 = IT_0316 + IT_0317 + -IT_0318 + -IT_0319;
    const ccomplex_t IT_0321 = IT_0293*IT_0320;
    const ccomplex_t IT_0322 = (-0.5)*IT_0315 + (-0.5)*IT_0321;
    const ccomplex_t IT_0323 = 6*IT_0322;
    const ccomplex_t IT_0324 = 0.5*IT_0315 + 0.5*IT_0321;
    const ccomplex_t IT_0325 = (-6)*IT_0324;
    const ccomplex_t IT_0326 = IT_0323 + IT_0325;
    const ccomplex_t IT_0327 = 6*conj(IT_0322);
    const ccomplex_t IT_0328 = (-6)*conj(IT_0324);
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = -IT_0255;
    const ccomplex_t IT_0331 = IT_0254 + IT_0256 + IT_0257 + IT_0330;
    const ccomplex_t IT_0332 = IT_0219*IT_0260;
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*IT_0332;
    const ccomplex_t IT_0334 = IT_0222*IT_0263;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = -IT_0333 + -IT_0335;
    const ccomplex_t IT_0337 = 6*IT_0297;
    const ccomplex_t IT_0338 = (-6)*IT_0295;
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = 6*conj(IT_0297);
    const ccomplex_t IT_0341 = (-6)*conj(IT_0295);
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = IT_0222*IT_0303;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*IT_0343;
    const ccomplex_t IT_0345 = IT_0219*IT_0306;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*IT_0345;
    const ccomplex_t IT_0347 = -IT_0344 + -IT_0346;
    const ccomplex_t IT_0348 = 6*IT_0324;
    const ccomplex_t IT_0349 = (-6)*IT_0322;
    const ccomplex_t IT_0350 = IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = 6*conj(IT_0324);
    const ccomplex_t IT_0352 = (-6)*conj(IT_0322);
    const ccomplex_t IT_0353 = IT_0351 + IT_0352;
    const ccomplex_t IT_0354 = s_23 + s_24;
    const ccomplex_t IT_0355 = IT_0000*IT_0354;
    const ccomplex_t IT_0356 = conj(IT_0214) + conj(IT_0250);
    const ccomplex_t IT_0357 = 12*IT_0297;
    const ccomplex_t IT_0358 = (-12)*IT_0295;
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = IT_0214 + IT_0250;
    const ccomplex_t IT_0361 = 12*conj(IT_0297);
    const ccomplex_t IT_0362 = (-12)*conj(IT_0295);
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = conj(IT_0216) + conj(IT_0249);
    const ccomplex_t IT_0365 = 12*IT_0295;
    const ccomplex_t IT_0366 = (-12)*IT_0297;
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = IT_0216 + IT_0249;
    const ccomplex_t IT_0369 = 12*conj(IT_0295);
    const ccomplex_t IT_0370 = (-12)*conj(IT_0297);
    const ccomplex_t IT_0371 = IT_0369 + IT_0370;
    const ccomplex_t IT_0372 = conj(IT_0205) + conj(IT_0245);
    const ccomplex_t IT_0373 = 12*IT_0322;
    const ccomplex_t IT_0374 = (-12)*IT_0324;
    const ccomplex_t IT_0375 = IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = IT_0205 + IT_0245;
    const ccomplex_t IT_0377 = 12*conj(IT_0322);
    const ccomplex_t IT_0378 = (-12)*conj(IT_0324);
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = conj(IT_0211) + conj(IT_0246);
    const ccomplex_t IT_0381 = 12*IT_0324;
    const ccomplex_t IT_0382 = (-12)*IT_0322;
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = IT_0211 + IT_0246;
    const ccomplex_t IT_0385 = 12*conj(IT_0324);
    const ccomplex_t IT_0386 = (-12)*conj(IT_0322);
    const ccomplex_t IT_0387 = IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = conj(IT_0156)*IT_0299 + IT_0156*IT_0302 + conj
      (IT_0076)*IT_0326 + IT_0076*IT_0329 + conj(IT_0236)*IT_0339 + IT_0236
      *IT_0342 + conj(IT_0225)*IT_0350 + IT_0225*IT_0353 + IT_0356*IT_0359 +
       IT_0360*IT_0363 + IT_0364*IT_0367 + IT_0368*IT_0371 + IT_0372*IT_0375 +
       IT_0376*IT_0379 + IT_0380*IT_0383 + IT_0384*IT_0387;
    const ccomplex_t IT_0389 = -IT_0256;
    const ccomplex_t IT_0390 = IT_0255 + IT_0257 + IT_0389;
    const ccomplex_t IT_0391 = 6*IT_0266;
    const ccomplex_t IT_0392 = 6*conj(IT_0266);
    const ccomplex_t IT_0393 = s_13 + s_14;
    const ccomplex_t IT_0394 = IT_0000*IT_0393;
    const ccomplex_t IT_0395 = 6*IT_0156;
    const ccomplex_t IT_0396 = 6*conj(IT_0156);
    const ccomplex_t IT_0397 = 6*IT_0236;
    const ccomplex_t IT_0398 = 6*conj(IT_0236);
    const ccomplex_t IT_0399 = conj(IT_0076)*IT_0299 + IT_0076*IT_0302 + conj
      (IT_0225)*IT_0339 + IT_0225*IT_0342 + IT_0367*IT_0372 + IT_0364*IT_0375 +
       IT_0371*IT_0376 + IT_0368*IT_0379 + IT_0359*IT_0380 + IT_0356*IT_0383 +
       IT_0363*IT_0384 + IT_0360*IT_0387 + (-6)*conj(IT_0322)*(IT_0236 + (
      -0.166666666666667)*IT_0395) + (-6)*IT_0322*(conj(IT_0236) + (
      -0.166666666666667)*IT_0396) + (-6)*conj(IT_0324)*(IT_0156 + (
      -0.166666666666667)*IT_0397) + (-6)*IT_0324*(conj(IT_0156) + (
      -0.166666666666667)*IT_0398);
    const ccomplex_t IT_0400 = m_b*s_14*m_N_1;
    const ccomplex_t IT_0401 = 12*IT_0236;
    const ccomplex_t IT_0402 = 24*IT_0214;
    const ccomplex_t IT_0403 = 24*IT_0250;
    const ccomplex_t IT_0404 = IT_0401 + IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = 12*conj(IT_0236);
    const ccomplex_t IT_0406 = 24*conj(IT_0214);
    const ccomplex_t IT_0407 = 24*conj(IT_0250);
    const ccomplex_t IT_0408 = IT_0405 + IT_0406 + IT_0407;
    const ccomplex_t IT_0409 = 24*IT_0246;
    const ccomplex_t IT_0410 = 24*IT_0211;
    const ccomplex_t IT_0411 = IT_0409 + IT_0410;
    const ccomplex_t IT_0412 = 24*conj(IT_0246);
    const ccomplex_t IT_0413 = 24*conj(IT_0211);
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = 12*IT_0336;
    const ccomplex_t IT_0416 = IT_0391 + IT_0415;
    const ccomplex_t IT_0417 = 12*conj(IT_0336);
    const ccomplex_t IT_0418 = IT_0392 + IT_0417;
    const ccomplex_t IT_0419 = 12*IT_0266;
    const ccomplex_t IT_0420 = 12*conj(IT_0266);
    const ccomplex_t IT_0421 = 12*IT_0309;
    const ccomplex_t IT_0422 = 12*conj(IT_0309);
    const ccomplex_t IT_0423 = conj(IT_0309)*IT_0397 + IT_0309*IT_0398 + conj
      (IT_0347)*IT_0404 + IT_0347*IT_0408 + conj(IT_0336)*IT_0411 + IT_0336
      *IT_0414 + conj(IT_0225)*IT_0416 + IT_0225*IT_0418 + IT_0380*IT_0419 +
       IT_0384*IT_0420 + IT_0356*IT_0421 + IT_0360*IT_0422;
    const ccomplex_t IT_0424 = m_b*s_24*m_N_1;
    const ccomplex_t IT_0425 = 24*IT_0205;
    const ccomplex_t IT_0426 = 24*IT_0245;
    const ccomplex_t IT_0427 = IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = 12*IT_0076;
    const ccomplex_t IT_0429 = IT_0427 + IT_0428;
    const ccomplex_t IT_0430 = 6*IT_0336;
    const ccomplex_t IT_0431 = IT_0419 + IT_0430;
    const ccomplex_t IT_0432 = 6*conj(IT_0336);
    const ccomplex_t IT_0433 = IT_0420 + IT_0432;
    const ccomplex_t IT_0434 = 6*IT_0076;
    const ccomplex_t IT_0435 = 6*conj(IT_0076);
    const ccomplex_t IT_0436 = 12*IT_0347;
    const ccomplex_t IT_0437 = 12*conj(IT_0347);
    const ccomplex_t IT_0438 = 24*conj(IT_0205);
    const ccomplex_t IT_0439 = 24*conj(IT_0245);
    const ccomplex_t IT_0440 = IT_0438 + IT_0439;
    const ccomplex_t IT_0441 = 12*conj(IT_0076);
    const ccomplex_t IT_0442 = IT_0364*(24*IT_0266 + IT_0415) + IT_0368*(24
      *conj(IT_0266) + IT_0417) + conj(IT_0309)*IT_0429 + conj(IT_0156)*IT_0431 
      + IT_0156*IT_0433 + conj(IT_0347)*IT_0434 + IT_0347*IT_0435 + IT_0372
      *IT_0436 + IT_0376*IT_0437 + IT_0309*(IT_0440 + IT_0441);
    const ccomplex_t IT_0443 = m_b*s_23*m_N_1;
    const ccomplex_t IT_0444 = 12*IT_0225;
    const ccomplex_t IT_0445 = 12*conj(IT_0225);
    const ccomplex_t IT_0446 = conj(IT_0225)*IT_0309 + IT_0225*conj(IT_0309) +
       0.166666666666667*conj(IT_0236)*IT_0416 + 0.166666666666667*IT_0236
      *IT_0418 + 4*IT_0356*(IT_0336 + 0.0416666666666667*IT_0419) + 4*IT_0360*
      (conj(IT_0336) + 0.0416666666666667*IT_0420) + 0.166666666666667*IT_0380
      *IT_0421 + 0.166666666666667*IT_0384*IT_0422 + 4*conj(IT_0347)*(IT_0384 +
       0.0416666666666667*IT_0444) + 4*IT_0347*(IT_0380 + 0.0416666666666667
      *IT_0445);
    const ccomplex_t IT_0447 = s_34*IT_0006;
    const ccomplex_t IT_0448 = IT_0411 + IT_0444;
    const ccomplex_t IT_0449 = 12*IT_0156;
    const ccomplex_t IT_0450 = 12*conj(IT_0156);
    const ccomplex_t IT_0451 = conj(IT_0225)*IT_0395 + IT_0225*IT_0396 + conj
      (IT_0076)*IT_0397 + IT_0076*IT_0398 + conj(IT_0205)*(IT_0401 + IT_0402) +
       conj(IT_0245)*IT_0404 + IT_0205*(IT_0405 + IT_0406) + IT_0245*IT_0408 +
       conj(IT_0250)*IT_0425 + IT_0356*IT_0428 + IT_0250*IT_0438 + IT_0360
      *IT_0441 + IT_0368*(IT_0414 + IT_0445) + IT_0364*IT_0448 + IT_0380*IT_0449
       + IT_0384*IT_0450;
    const ccomplex_t IT_0452 = 48*IT_0205;
    const ccomplex_t IT_0453 = 24*IT_0076;
    const ccomplex_t IT_0454 = IT_0452 + IT_0453;
    const ccomplex_t IT_0455 = 48*IT_0245;
    const ccomplex_t IT_0456 = IT_0454 + IT_0455;
    const ccomplex_t IT_0457 = 24*IT_0216;
    const ccomplex_t IT_0458 = IT_0449 + IT_0457;
    const ccomplex_t IT_0459 = 48*IT_0216;
    const ccomplex_t IT_0460 = 48*conj(IT_0216);
    const ccomplex_t IT_0461 = 48*conj(IT_0205);
    const ccomplex_t IT_0462 = conj(IT_0266)*IT_0419 + conj(IT_0309)*IT_0421 +
       conj(IT_0076)*IT_0429 + conj(IT_0205)*IT_0454 + conj(IT_0245)*IT_0456 +
       conj(IT_0156)*IT_0458 + IT_0364*(24*IT_0156 + IT_0459) + IT_0249*(24*conj
      (IT_0156) + 48*conj(IT_0249) + IT_0460) + IT_0245*IT_0461;
    const ccomplex_t IT_0463 = -s_34;
    const ccomplex_t IT_0464 = IT_0005 + IT_0463;
    const ccomplex_t IT_0465 = s_12*IT_0464;
    const ccomplex_t IT_0466 = (-24)*conj(IT_0176);
    const ccomplex_t IT_0467 = 24*conj(IT_0174);
    const ccomplex_t IT_0468 = IT_0466 + IT_0467;
    const ccomplex_t IT_0469 = (-24)*IT_0112;
    const ccomplex_t IT_0470 = s_34 + IT_0005;
    const ccomplex_t IT_0471 = s_12*IT_0470;
    const ccomplex_t IT_0472 = (-24)*IT_0297;
    const ccomplex_t IT_0473 = (-24)*conj(IT_0297);
    const ccomplex_t IT_0474 = IT_0006*IT_0464;
    const ccomplex_t IT_0475 = 6*conj(IT_0110);
    const ccomplex_t IT_0476 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0477 = IT_0475 + IT_0476;
    const ccomplex_t IT_0478 = 6*conj(IT_0112);
    const ccomplex_t IT_0479 = (-6)*conj(IT_0110);
    const ccomplex_t IT_0480 = IT_0478 + IT_0479;
    const ccomplex_t IT_0481 = 6*IT_0174;
    const ccomplex_t IT_0482 = (-6)*IT_0176;
    const ccomplex_t IT_0483 = IT_0481 + IT_0482;
    const ccomplex_t IT_0484 = 6*IT_0176;
    const ccomplex_t IT_0485 = (-6)*IT_0174;
    const ccomplex_t IT_0486 = IT_0484 + IT_0485;
    const ccomplex_t IT_0487 = 6*conj(IT_0176);
    const ccomplex_t IT_0488 = (-6)*conj(IT_0174);
    const ccomplex_t IT_0489 = IT_0487 + IT_0488;
    const ccomplex_t IT_0490 = 6*IT_0112;
    const ccomplex_t IT_0491 = (-6)*IT_0112;
    const ccomplex_t IT_0492 = (-6)*conj(IT_0176);
    const ccomplex_t IT_0493 = conj(IT_0112)*((-24)*IT_0174 + 24*IT_0176) +
       conj(IT_0110)*(24*IT_0174 + (-24)*IT_0176) + 24*IT_0112*conj(IT_0176) + 6
      *IT_0110*(conj(IT_0309) + -conj(IT_0347) + 0.166666666666667*IT_0468) + 6
      *conj(IT_0174)*(IT_0266 + 0.166666666666667*IT_0469) + IT_0309*IT_0477 +
       IT_0347*IT_0480 + conj(IT_0266)*IT_0483 + conj(IT_0336)*IT_0486 + IT_0336
      *IT_0489 + conj(IT_0347)*IT_0490 + conj(IT_0309)*IT_0491 + IT_0266*IT_0492;
    const ccomplex_t IT_0494 = conj(IT_0236)*IT_0395 + IT_0236*IT_0396 +
       IT_0364*IT_0404 + IT_0368*IT_0408 + IT_0380*(IT_0425 + IT_0426 + IT_0428)
       + conj(IT_0225)*IT_0434 + IT_0225*IT_0435 + IT_0384*(IT_0438 + IT_0439 +
       IT_0441) + IT_0372*IT_0444 + IT_0376*IT_0445 + IT_0356*IT_0449 + IT_0360
      *IT_0450;
    const ccomplex_t IT_0495 = m_b*s_13*m_N_1;
    const ccomplex_t IT_0496 = 24*conj(IT_0216);
    const ccomplex_t IT_0497 = conj(IT_0347)*IT_0395 + IT_0347*IT_0396 +
       IT_0372*IT_0415 + IT_0376*IT_0417 + conj(IT_0266)*IT_0427 + conj(IT_0076)
      *IT_0431 + IT_0076*IT_0433 + IT_0364*IT_0436 + IT_0368*IT_0437 + IT_0266
      *IT_0440 + conj(IT_0309)*(24*IT_0249 + IT_0449 + IT_0457) + IT_0309*(24
      *conj(IT_0249) + IT_0450 + IT_0496);
    const ccomplex_t IT_0498 = IT_0005*IT_0006;
    const ccomplex_t IT_0499 = 48*IT_0214;
    const ccomplex_t IT_0500 = 24*IT_0236;
    const ccomplex_t IT_0501 = IT_0499 + IT_0500;
    const ccomplex_t IT_0502 = 48*conj(IT_0214);
    const ccomplex_t IT_0503 = 48*conj(IT_0211);
    const ccomplex_t IT_0504 = 48*IT_0211;
    const ccomplex_t IT_0505 = 48*IT_0246;
    const ccomplex_t IT_0506 = 48*IT_0250;
    const ccomplex_t IT_0507 = conj(IT_0225)*IT_0404 + IT_0225*IT_0408 + conj
      (IT_0236)*IT_0411 + conj(IT_0347)*IT_0416 + IT_0347*IT_0418 + conj(IT_0156
      )*IT_0427 + conj(IT_0309)*IT_0431 + IT_0309*IT_0433 + IT_0156*IT_0440 +
       conj(IT_0249)*IT_0456 + conj(IT_0076)*IT_0458 + IT_0372*IT_0459 + IT_0376
      *IT_0460 + IT_0249*(24*conj(IT_0076) + 48*conj(IT_0245) + IT_0461) +
       IT_0076*(IT_0450 + IT_0496) + conj(IT_0211)*IT_0501 + IT_0384*IT_0502 +
       IT_0250*IT_0503 + conj(IT_0250)*(IT_0504 + IT_0505) + conj(IT_0246)*
      (IT_0499 + IT_0500 + IT_0506);
    const ccomplex_t IT_0508 = 24*IT_0225;
    const ccomplex_t IT_0509 = IT_0504 + IT_0508;
    const ccomplex_t IT_0510 = conj(IT_0236)*IT_0404 + conj(IT_0336)*IT_0415 +
       conj(IT_0347)*IT_0436 + conj(IT_0225)*IT_0448 + conj(IT_0214)*IT_0501 +
       IT_0250*IT_0502 + IT_0246*IT_0503 + conj(IT_0250)*(IT_0499 + IT_0500 +
       IT_0506) + conj(IT_0211)*IT_0509 + conj(IT_0246)*(IT_0505 + IT_0509);
    const ccomplex_t IT_0511 = (-0.166666666666667)*IT_0309;
    const ccomplex_t IT_0512 = (-0.166666666666667)*conj(IT_0309);
    const ccomplex_t IT_0513 = IT_0006*IT_0470;
    const ccomplex_t IT_0514 = IT_0342*IT_0347 + IT_0339*conj(IT_0347) + (-6)
      *conj(IT_0324)*(IT_0266 + 4*IT_0295 + (-4)*IT_0297 + (-0.166666666666667)
      *IT_0430) + (-6)*IT_0324*(conj(IT_0266) + 4*conj(IT_0295) + (-4)*conj
      (IT_0297) + (-0.166666666666667)*IT_0432) + 24*conj(IT_0322)*(IT_0295 + 
      -1./4*IT_0336 + 0.0416666666666667*IT_0391 + 0.0416666666666667*IT_0472) +
       24*IT_0322*(conj(IT_0295) + -1./4*conj(IT_0336) + 0.0416666666666667
      *IT_0392 + 0.0416666666666667*IT_0473) + (-6)*IT_0302*IT_0511 + (-6)
      *IT_0299*IT_0512;
    const ccomplex_t IT_0515 = IT_0004*(conj(IT_0076)*IT_0114 + IT_0076
      *IT_0117 + conj(IT_0156)*IT_0178 + IT_0156*IT_0181) + (IT_0117*IT_0156 +
       IT_0114*conj(IT_0156) + conj(IT_0076)*IT_0178 + IT_0076*IT_0181)*IT_0185 
      + IT_0186*(IT_0117*(IT_0205 + IT_0211) + IT_0114*(conj(IT_0205) + conj
      (IT_0211)) + IT_0181*(IT_0214 + IT_0216) + IT_0178*(conj(IT_0214) + conj
      (IT_0216))) + (IT_0181*(IT_0205 + IT_0211) + IT_0178*(conj(IT_0205) + conj
      (IT_0211)) + IT_0117*(IT_0214 + IT_0216) + IT_0114*(conj(IT_0214) + conj
      (IT_0216)))*IT_0217 + IT_0218*(conj(IT_0225)*IT_0228 + IT_0225*IT_0231 +
       conj(IT_0236)*IT_0239 + IT_0236*IT_0242) + (IT_0231*IT_0236 + IT_0228
      *conj(IT_0236) + conj(IT_0225)*IT_0239 + IT_0225*IT_0242)*IT_0243 +
       IT_0244*(IT_0228*IT_0247 + IT_0231*IT_0248 + IT_0239*IT_0251 + IT_0242
      *IT_0252) + (IT_0239*IT_0247 + IT_0242*IT_0248 + IT_0228*IT_0251 + IT_0231
      *IT_0252)*IT_0253 + IT_0259*(conj(IT_0266)*IT_0299 + IT_0266*IT_0302 +
       conj(IT_0309)*IT_0326 + IT_0309*IT_0329) + IT_0331*(conj(IT_0336)*IT_0339
       + IT_0336*IT_0342 + conj(IT_0347)*IT_0350 + IT_0347*IT_0353) + IT_0355
      *IT_0388 + 6*IT_0390*(conj(IT_0309)*IT_0347 + IT_0309*conj(IT_0347) +
       0.166666666666667*conj(IT_0336)*IT_0391 + 0.166666666666667*IT_0336
      *IT_0392) + IT_0394*IT_0399 + IT_0400*IT_0423 + IT_0424*IT_0442 + 6
      *IT_0443*IT_0446 + IT_0447*IT_0451 + IT_0255*IT_0462 + IT_0465*(24*conj
      (IT_0110)*IT_0112 + IT_0110*((-24)*conj(IT_0110) + 24*conj(IT_0112)) +
       IT_0174*((-24)*conj(IT_0174) + 24*conj(IT_0176)) + IT_0176*IT_0468 + conj
      (IT_0112)*IT_0469) + IT_0471*(24*IT_0297*conj(IT_0297) + IT_0324*((-24)
      *conj(IT_0322) + 24*conj(IT_0324)) + IT_0322*(24*conj(IT_0322) + (-24)
      *conj(IT_0324)) + conj(IT_0295)*IT_0472 + IT_0295*(24*conj(IT_0295) +
       IT_0473)) + (IT_0254 + IT_0255 + IT_0258 + IT_0389)*(IT_0336*IT_0477 +
       conj(IT_0347)*IT_0483 + 6*conj(IT_0336)*(IT_0110 + 0.166666666666667
      *IT_0491) + 6*IT_0347*(conj(IT_0174) + 0.166666666666667*IT_0492)) +
       IT_0474*IT_0493 + IT_0254*IT_0494 + IT_0495*IT_0497 + IT_0498*IT_0507 +
       IT_0257*IT_0510 + (IT_0254 + IT_0257 + IT_0330 + IT_0389)*(IT_0266
      *IT_0480 + (-6)*conj(IT_0266)*(IT_0110 + (-0.166666666666667)*IT_0490) + (
      -6)*IT_0489*IT_0511 + (-6)*IT_0486*IT_0512) + IT_0513*IT_0514;
    return create_ccomplex_return(IT_0515);
}

