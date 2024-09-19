#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_3_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_3_to_anti_b_b(
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
    const creal_t m_N_3 = param->m_N_3;
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
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = m_b*m_N_3;
    const ccomplex_t IT_0001 = -s_24;
    const ccomplex_t IT_0002 = s_23 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 2*IT_0003;
    const ccomplex_t IT_0005 = pow(m_b, 2);
    const ccomplex_t IT_0006 = pow(m_N_3, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_b*conj(N_d3)*e_em*IT_0018*conj(U_sb_00);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0012 + 1.5*IT_0021);
    const ccomplex_t IT_0023 = 0.666666666666667*IT_0022;
    const ccomplex_t IT_0024 = 1.5*IT_0023;
    const ccomplex_t IT_0025 = (-0.333333333333333)*IT_0024;
    const ccomplex_t IT_0026 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = m_b*N_d3*e_em*IT_0018*U_sb_00;
    const ccomplex_t IT_0030 = IT_0017*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + 1.5*IT_0031);
    const ccomplex_t IT_0033 = 0.666666666666667*IT_0032;
    const ccomplex_t IT_0034 = 1.5*IT_0033;
    const ccomplex_t IT_0035 = (-0.333333333333333)*IT_0034;
    const ccomplex_t IT_0036 = IT_0025*IT_0035;
    const ccomplex_t IT_0037 = IT_0007*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_23 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0040 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0041 = IT_0009*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_b*conj(N_d3)*e_em*IT_0018*conj(U_sb_01);
    const ccomplex_t IT_0044 = IT_0017*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 1.5*IT_0045);
    const ccomplex_t IT_0047 = (-0.333333333333333)*IT_0046;
    const ccomplex_t IT_0048 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0049 = IT_0009*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = m_b*N_d3*e_em*IT_0018*U_sb_01;
    const ccomplex_t IT_0052 = IT_0017*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + 1.5*IT_0053);
    const ccomplex_t IT_0055 = 0.666666666666667*IT_0054;
    const ccomplex_t IT_0056 = 1.5*IT_0055;
    const ccomplex_t IT_0057 = 0.666666666666667*IT_0056;
    const ccomplex_t IT_0058 = 1.5*IT_0057;
    const ccomplex_t IT_0059 = 0.666666666666667*IT_0058;
    const ccomplex_t IT_0060 = 1.5*IT_0059;
    const ccomplex_t IT_0061 = (-0.333333333333333)*IT_0060;
    const ccomplex_t IT_0062 = IT_0047*IT_0061;
    const ccomplex_t IT_0063 = IT_0039*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0038 + IT_0064;
    const ccomplex_t IT_0066 = cos(alpha);
    const ccomplex_t IT_0067 = IT_0009*IT_0066;
    const ccomplex_t IT_0068 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = sin(alpha);
    const ccomplex_t IT_0071 = IT_0009*IT_0070;
    const ccomplex_t IT_0072 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = IT_0016*IT_0066;
    const ccomplex_t IT_0075 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = IT_0016*IT_0070;
    const ccomplex_t IT_0078 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0069 + -IT_0073 + 
      -IT_0076 + IT_0079);
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_b*e_em*IT_0014*IT_0016
      *IT_0018*IT_0066;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + (-2)*IT_0006 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0080*IT_0085;
    const ccomplex_t IT_0087 = IT_0067*IT_0072;
    const ccomplex_t IT_0088 = IT_0068*IT_0071;
    const ccomplex_t IT_0089 = IT_0074*IT_0078;
    const ccomplex_t IT_0090 = IT_0075*IT_0077;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + IT_0088 + 
      -IT_0089 + -IT_0090);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_b*e_em*IT_0014*IT_0016
      *IT_0018*IT_0070;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0092*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0086 + 0.5*IT_0098;
    const ccomplex_t IT_0100 = 3*IT_0099;
    const ccomplex_t IT_0101 = (-0.5)*IT_0086 + (-0.5)*IT_0098;
    const ccomplex_t IT_0102 = (-3)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = 3*conj(IT_0099);
    const ccomplex_t IT_0105 = (-3)*conj(IT_0101);
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = m_b*N_d3*e_em*IT_0018*conj(U_sb_10);
    const ccomplex_t IT_0108 = IT_0017*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0111 = IT_0009*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0114 = IT_0016*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0109 + 0.333333333333333
      *IT_0112 + -IT_0115);
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = (-3)*IT_0117;
    const ccomplex_t IT_0119 = 0.166666666666667*IT_0118;
    const ccomplex_t IT_0120 = m_b*conj(N_d3)*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0121 = IT_0017*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0124 = IT_0009*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0127 = IT_0016*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0122 + 0.333333333333333
      *IT_0125 + -IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = 3*IT_0133;
    const ccomplex_t IT_0135 = (-0.333333333333333)*IT_0134;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = 3*IT_0136;
    const ccomplex_t IT_0138 = (-0.333333333333333)*IT_0137;
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = 3*IT_0143;
    const ccomplex_t IT_0145 = (-0.333333333333333)*IT_0144;
    const ccomplex_t IT_0146 = (-3)*IT_0145;
    const ccomplex_t IT_0147 = 0.333333333333333*IT_0146;
    const ccomplex_t IT_0148 = 3*IT_0147;
    const ccomplex_t IT_0149 = 0.166666666666667*IT_0148;
    const ccomplex_t IT_0150 = IT_0119*IT_0149;
    const ccomplex_t IT_0151 = IT_0007*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_b*conj(N_d3)*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0154 = IT_0017*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0157 = IT_0009*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0160 = IT_0016*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0155 + 0.333333333333333
      *IT_0158 + -IT_0161);
    const ccomplex_t IT_0163 = 3*IT_0162;
    const ccomplex_t IT_0164 = 0.166666666666667*IT_0163;
    const ccomplex_t IT_0165 = m_b*N_d3*e_em*IT_0018*conj(U_sb_11);
    const ccomplex_t IT_0166 = IT_0017*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0169 = IT_0016*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0172 = IT_0009*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + -IT_0170 +
       0.333333333333333*IT_0173);
    const ccomplex_t IT_0175 = 3*IT_0174;
    const ccomplex_t IT_0176 = 0.166666666666667*IT_0175;
    const ccomplex_t IT_0177 = IT_0164*IT_0176;
    const ccomplex_t IT_0178 = IT_0039*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = IT_0152 + IT_0179;
    const ccomplex_t IT_0181 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0182 = IT_0067*IT_0181;
    const ccomplex_t IT_0183 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0184 = IT_0071*IT_0183;
    const ccomplex_t IT_0185 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0186 = IT_0074*IT_0185;
    const ccomplex_t IT_0187 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0188 = IT_0077*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0182 + IT_0184 + 
      -IT_0186 + -IT_0188);
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = IT_0097*IT_0190;
    const ccomplex_t IT_0192 = IT_0067*IT_0183;
    const ccomplex_t IT_0193 = IT_0071*IT_0181;
    const ccomplex_t IT_0194 = IT_0074*IT_0187;
    const ccomplex_t IT_0195 = IT_0077*IT_0185;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0192 + -IT_0193 + 
      -IT_0194 + IT_0195);
    const ccomplex_t IT_0197 = IT_0085*IT_0196;
    const ccomplex_t IT_0198 = 0.5*IT_0191 + 0.5*IT_0197;
    const ccomplex_t IT_0199 = 3*IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0191 + (-0.5)*IT_0197;
    const ccomplex_t IT_0201 = (-3)*IT_0200;
    const ccomplex_t IT_0202 = IT_0199 + IT_0201;
    const ccomplex_t IT_0203 = 3*conj(IT_0198);
    const ccomplex_t IT_0204 = (-3)*conj(IT_0200);
    const ccomplex_t IT_0205 = IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = -s_14;
    const ccomplex_t IT_0207 = s_13 + IT_0206;
    const ccomplex_t IT_0208 = IT_0000*IT_0207;
    const ccomplex_t IT_0209 = (-2)*IT_0208;
    const ccomplex_t IT_0210 = (-4)*IT_0208;
    const ccomplex_t IT_0211 = IT_0009*IT_0015;
    const ccomplex_t IT_0212 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0213 = IT_0211*IT_0212;
    const ccomplex_t IT_0214 = IT_0008*IT_0016;
    const ccomplex_t IT_0215 = IT_0212*IT_0214;
    const ccomplex_t IT_0216 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0217 = IT_0211*IT_0216;
    const ccomplex_t IT_0218 = IT_0214*IT_0216;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0213 + IT_0215 + 
      -IT_0217 + -IT_0218);
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = e_em*IT_0211;
    const ccomplex_t IT_0222 = e_em*IT_0214;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*(IT_0221 + 3*IT_0222);
    const ccomplex_t IT_0224 = (-0.166666666666667)*IT_0223;
    const ccomplex_t IT_0225 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0226 = IT_0224*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = IT_0220*IT_0227;
    const ccomplex_t IT_0229 = (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0015;
    const ccomplex_t IT_0231 = 0.333333333333333*IT_0230;
    const ccomplex_t IT_0232 = IT_0225*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = IT_0220*IT_0233;
    const ccomplex_t IT_0235 = (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = (-0.5)*IT_0219;
    const ccomplex_t IT_0237 = IT_0233*IT_0236;
    const ccomplex_t IT_0238 = (-0.5)*IT_0237;
    const ccomplex_t IT_0239 = IT_0227*IT_0236;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = 4*IT_0003;
    const ccomplex_t IT_0242 = (-2)*IT_0003;
    const ccomplex_t IT_0243 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0244 = IT_0036*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = cpow((-2)*s_13 + IT_0005 + IT_0006 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0247 = IT_0062*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = -IT_0245 + -IT_0248;
    const ccomplex_t IT_0250 = 3*IT_0200;
    const ccomplex_t IT_0251 = (-3)*IT_0198;
    const ccomplex_t IT_0252 = IT_0250 + IT_0251;
    const ccomplex_t IT_0253 = 3*conj(IT_0200);
    const ccomplex_t IT_0254 = (-3)*conj(IT_0198);
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = IT_0150*IT_0243;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = IT_0177*IT_0246;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = -IT_0257 + -IT_0259;
    const ccomplex_t IT_0261 = 3*IT_0101;
    const ccomplex_t IT_0262 = (-3)*IT_0099;
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 3*conj(IT_0101);
    const ccomplex_t IT_0265 = (-3)*conj(IT_0099);
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = 2*IT_0208;
    const ccomplex_t IT_0268 = 4*IT_0208;
    const ccomplex_t IT_0269 = 0.5*IT_0228;
    const ccomplex_t IT_0270 = 0.5*IT_0234;
    const ccomplex_t IT_0271 = conj(IT_0269) + conj(IT_0270);
    const ccomplex_t IT_0272 = IT_0269 + IT_0270;
    const ccomplex_t IT_0273 = 0.5*IT_0239;
    const ccomplex_t IT_0274 = 0.5*IT_0237;
    const ccomplex_t IT_0275 = conj(IT_0273) + conj(IT_0274);
    const ccomplex_t IT_0276 = IT_0273 + IT_0274;
    const ccomplex_t IT_0277 = (-4)*IT_0003;
    const ccomplex_t IT_0278 = s_12*IT_0005;
    const ccomplex_t IT_0279 = s_14*s_23;
    const ccomplex_t IT_0280 = s_12*s_34;
    const ccomplex_t IT_0281 = s_13*s_24;
    const ccomplex_t IT_0282 = -IT_0281;
    const ccomplex_t IT_0283 = IT_0278 + IT_0279 + IT_0280 + IT_0282;
    const ccomplex_t IT_0284 = IT_0035*IT_0119;
    const ccomplex_t IT_0285 = IT_0007*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0061*IT_0176;
    const ccomplex_t IT_0288 = IT_0039*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0286 + -IT_0289;
    const ccomplex_t IT_0291 = IT_0009*IT_0013;
    const ccomplex_t IT_0292 = IT_0183*IT_0291;
    const ccomplex_t IT_0293 = sin(beta);
    const ccomplex_t IT_0294 = IT_0009*IT_0293;
    const ccomplex_t IT_0295 = IT_0181*IT_0294;
    const ccomplex_t IT_0296 = IT_0013*IT_0016;
    const ccomplex_t IT_0297 = IT_0187*IT_0296;
    const ccomplex_t IT_0298 = IT_0016*IT_0293;
    const ccomplex_t IT_0299 = IT_0185*IT_0298;
    const ccomplex_t IT_0300 = -IT_0292 + -IT_0295 + IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = m_b*e_em*IT_0016*IT_0018;
    const ccomplex_t IT_0302 = 0.5*IT_0301;
    const ccomplex_t IT_0303 = pow(m_Z, 2);
    const ccomplex_t IT_0304 = cpow((-2)*s_12 + (-2)*IT_0006 + IT_0303 + 
      -reg_prop, -1);
    const ccomplex_t IT_0305 = IT_0302*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = IT_0300*IT_0306;
    const ccomplex_t IT_0308 = IT_0181*IT_0291;
    const ccomplex_t IT_0309 = IT_0183*IT_0294;
    const ccomplex_t IT_0310 = IT_0185*IT_0296;
    const ccomplex_t IT_0311 = IT_0187*IT_0298;
    const ccomplex_t IT_0312 = -IT_0308 + IT_0309 + IT_0310 + -IT_0311;
    const ccomplex_t IT_0313 = m_b*e_em*IT_0014*IT_0016*IT_0018*IT_0293;
    const ccomplex_t IT_0314 = (-0.5)*IT_0313;
    const ccomplex_t IT_0315 = cpow((-2)*s_12 + (-2)*IT_0006 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0316 = IT_0314*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = IT_0312*IT_0317;
    const ccomplex_t IT_0319 = (-0.5)*IT_0307 + (-0.5)*IT_0318;
    const ccomplex_t IT_0320 = 6*IT_0319;
    const ccomplex_t IT_0321 = 0.5*IT_0307 + 0.5*IT_0318;
    const ccomplex_t IT_0322 = (-6)*IT_0321;
    const ccomplex_t IT_0323 = IT_0320 + IT_0322;
    const ccomplex_t IT_0324 = 6*conj(IT_0319);
    const ccomplex_t IT_0325 = (-6)*conj(IT_0321);
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = IT_0025*IT_0149;
    const ccomplex_t IT_0328 = IT_0007*IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*IT_0328;
    const ccomplex_t IT_0330 = IT_0047*IT_0164;
    const ccomplex_t IT_0331 = IT_0039*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = -IT_0329 + -IT_0332;
    const ccomplex_t IT_0334 = IT_0072*IT_0291;
    const ccomplex_t IT_0335 = IT_0068*IT_0294;
    const ccomplex_t IT_0336 = IT_0078*IT_0296;
    const ccomplex_t IT_0337 = IT_0075*IT_0298;
    const ccomplex_t IT_0338 = IT_0334 + -IT_0335 + -IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = IT_0317*IT_0338;
    const ccomplex_t IT_0340 = IT_0068*IT_0291;
    const ccomplex_t IT_0341 = IT_0072*IT_0294;
    const ccomplex_t IT_0342 = IT_0075*IT_0296;
    const ccomplex_t IT_0343 = IT_0078*IT_0298;
    const ccomplex_t IT_0344 = IT_0340 + IT_0341 + -IT_0342 + -IT_0343;
    const ccomplex_t IT_0345 = IT_0306*IT_0344;
    const ccomplex_t IT_0346 = 0.5*IT_0339 + 0.5*IT_0345;
    const ccomplex_t IT_0347 = 6*IT_0346;
    const ccomplex_t IT_0348 = (-0.5)*IT_0339 + (-0.5)*IT_0345;
    const ccomplex_t IT_0349 = (-6)*IT_0348;
    const ccomplex_t IT_0350 = IT_0347 + IT_0349;
    const ccomplex_t IT_0351 = 6*conj(IT_0346);
    const ccomplex_t IT_0352 = (-6)*conj(IT_0348);
    const ccomplex_t IT_0353 = IT_0351 + IT_0352;
    const ccomplex_t IT_0354 = -IT_0279;
    const ccomplex_t IT_0355 = IT_0278 + IT_0280 + IT_0281 + IT_0354;
    const ccomplex_t IT_0356 = IT_0243*IT_0284;
    const ccomplex_t IT_0357 = (0 + _Complex_I*1)*IT_0356;
    const ccomplex_t IT_0358 = IT_0246*IT_0287;
    const ccomplex_t IT_0359 = (0 + _Complex_I*1)*IT_0358;
    const ccomplex_t IT_0360 = IT_0357 + IT_0359;
    const ccomplex_t IT_0361 = 6*IT_0321;
    const ccomplex_t IT_0362 = (-6)*IT_0319;
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = 6*conj(IT_0321);
    const ccomplex_t IT_0365 = (-6)*conj(IT_0319);
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = IT_0243*IT_0327;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*IT_0367;
    const ccomplex_t IT_0369 = IT_0246*IT_0330;
    const ccomplex_t IT_0370 = (0 + _Complex_I*1)*IT_0369;
    const ccomplex_t IT_0371 = IT_0368 + IT_0370;
    const ccomplex_t IT_0372 = 6*IT_0348;
    const ccomplex_t IT_0373 = (-6)*IT_0346;
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = 6*conj(IT_0348);
    const ccomplex_t IT_0376 = (-6)*conj(IT_0346);
    const ccomplex_t IT_0377 = IT_0375 + IT_0376;
    const ccomplex_t IT_0378 = s_13 + s_14;
    const ccomplex_t IT_0379 = IT_0000*IT_0378;
    const ccomplex_t IT_0380 = conj(IT_0240) + conj(IT_0269);
    const ccomplex_t IT_0381 = 12*IT_0321;
    const ccomplex_t IT_0382 = (-12)*IT_0319;
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = IT_0240 + IT_0269;
    const ccomplex_t IT_0385 = 12*conj(IT_0321);
    const ccomplex_t IT_0386 = (-12)*conj(IT_0319);
    const ccomplex_t IT_0387 = IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = conj(IT_0238) + conj(IT_0270);
    const ccomplex_t IT_0389 = 12*IT_0319;
    const ccomplex_t IT_0390 = (-12)*IT_0321;
    const ccomplex_t IT_0391 = IT_0389 + IT_0390;
    const ccomplex_t IT_0392 = IT_0238 + IT_0270;
    const ccomplex_t IT_0393 = 12*conj(IT_0319);
    const ccomplex_t IT_0394 = (-12)*conj(IT_0321);
    const ccomplex_t IT_0395 = IT_0393 + IT_0394;
    const ccomplex_t IT_0396 = conj(IT_0229) + conj(IT_0273);
    const ccomplex_t IT_0397 = 12*IT_0346;
    const ccomplex_t IT_0398 = (-12)*IT_0348;
    const ccomplex_t IT_0399 = IT_0397 + IT_0398;
    const ccomplex_t IT_0400 = IT_0229 + IT_0273;
    const ccomplex_t IT_0401 = 12*conj(IT_0346);
    const ccomplex_t IT_0402 = (-12)*conj(IT_0348);
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = conj(IT_0235) + conj(IT_0274);
    const ccomplex_t IT_0405 = 12*IT_0348;
    const ccomplex_t IT_0406 = (-12)*IT_0346;
    const ccomplex_t IT_0407 = IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = IT_0235 + IT_0274;
    const ccomplex_t IT_0409 = 12*conj(IT_0348);
    const ccomplex_t IT_0410 = (-12)*conj(IT_0346);
    const ccomplex_t IT_0411 = IT_0409 + IT_0410;
    const ccomplex_t IT_0412 = conj(IT_0065)*IT_0323 + IT_0065*IT_0326 + conj
      (IT_0180)*IT_0350 + IT_0180*IT_0353 + conj(IT_0260)*IT_0363 + IT_0260
      *IT_0366 + conj(IT_0249)*IT_0374 + IT_0249*IT_0377 + IT_0380*IT_0383 +
       IT_0384*IT_0387 + IT_0388*IT_0391 + IT_0392*IT_0395 + IT_0396*IT_0399 +
       IT_0400*IT_0403 + IT_0404*IT_0407 + IT_0408*IT_0411;
    const ccomplex_t IT_0413 = s_23 + s_24;
    const ccomplex_t IT_0414 = IT_0000*IT_0413;
    const ccomplex_t IT_0415 = 6*IT_0065;
    const ccomplex_t IT_0416 = 6*conj(IT_0065);
    const ccomplex_t IT_0417 = 6*IT_0260;
    const ccomplex_t IT_0418 = 6*conj(IT_0260);
    const ccomplex_t IT_0419 = conj(IT_0180)*IT_0323 + IT_0180*IT_0326 + conj
      (IT_0249)*IT_0363 + IT_0249*IT_0366 + IT_0391*IT_0396 + IT_0388*IT_0399 +
       IT_0395*IT_0400 + IT_0392*IT_0403 + IT_0383*IT_0404 + IT_0380*IT_0407 +
       IT_0387*IT_0408 + IT_0384*IT_0411 + (-6)*conj(IT_0346)*(IT_0260 + (
      -0.166666666666667)*IT_0415) + (-6)*IT_0346*(conj(IT_0260) + (
      -0.166666666666667)*IT_0416) + (-6)*conj(IT_0348)*(IT_0065 + (
      -0.166666666666667)*IT_0417) + (-6)*IT_0348*(conj(IT_0065) + (
      -0.166666666666667)*IT_0418);
    const ccomplex_t IT_0420 = m_b*s_24*m_N_3;
    const ccomplex_t IT_0421 = 12*IT_0260;
    const ccomplex_t IT_0422 = 24*IT_0240;
    const ccomplex_t IT_0423 = 24*IT_0269;
    const ccomplex_t IT_0424 = IT_0421 + IT_0422 + IT_0423;
    const ccomplex_t IT_0425 = 12*conj(IT_0260);
    const ccomplex_t IT_0426 = 24*conj(IT_0240);
    const ccomplex_t IT_0427 = 24*conj(IT_0269);
    const ccomplex_t IT_0428 = IT_0425 + IT_0426 + IT_0427;
    const ccomplex_t IT_0429 = 24*IT_0274;
    const ccomplex_t IT_0430 = 24*IT_0235;
    const ccomplex_t IT_0431 = IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 24*conj(IT_0274);
    const ccomplex_t IT_0433 = 24*conj(IT_0235);
    const ccomplex_t IT_0434 = IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = 6*IT_0290;
    const ccomplex_t IT_0436 = 12*IT_0360;
    const ccomplex_t IT_0437 = IT_0435 + IT_0436;
    const ccomplex_t IT_0438 = 6*conj(IT_0290);
    const ccomplex_t IT_0439 = 12*conj(IT_0360);
    const ccomplex_t IT_0440 = IT_0438 + IT_0439;
    const ccomplex_t IT_0441 = 12*IT_0290;
    const ccomplex_t IT_0442 = 12*conj(IT_0290);
    const ccomplex_t IT_0443 = 12*IT_0333;
    const ccomplex_t IT_0444 = 12*conj(IT_0333);
    const ccomplex_t IT_0445 = conj(IT_0333)*IT_0417 + IT_0333*IT_0418 + conj
      (IT_0371)*IT_0424 + IT_0371*IT_0428 + conj(IT_0360)*IT_0431 + IT_0360
      *IT_0434 + conj(IT_0249)*IT_0437 + IT_0249*IT_0440 + IT_0404*IT_0441 +
       IT_0408*IT_0442 + IT_0380*IT_0443 + IT_0384*IT_0444;
    const ccomplex_t IT_0446 = m_b*s_23*m_N_3;
    const ccomplex_t IT_0447 = 12*IT_0065;
    const ccomplex_t IT_0448 = 24*IT_0270;
    const ccomplex_t IT_0449 = 24*IT_0238;
    const ccomplex_t IT_0450 = IT_0447 + IT_0448 + IT_0449;
    const ccomplex_t IT_0451 = 12*conj(IT_0065);
    const ccomplex_t IT_0452 = 24*conj(IT_0270);
    const ccomplex_t IT_0453 = 24*conj(IT_0238);
    const ccomplex_t IT_0454 = IT_0451 + IT_0452 + IT_0453;
    const ccomplex_t IT_0455 = 24*IT_0229;
    const ccomplex_t IT_0456 = 24*IT_0273;
    const ccomplex_t IT_0457 = IT_0455 + IT_0456;
    const ccomplex_t IT_0458 = 24*conj(IT_0229);
    const ccomplex_t IT_0459 = 24*conj(IT_0273);
    const ccomplex_t IT_0460 = IT_0458 + IT_0459;
    const ccomplex_t IT_0461 = 6*IT_0360;
    const ccomplex_t IT_0462 = IT_0441 + IT_0461;
    const ccomplex_t IT_0463 = 6*conj(IT_0360);
    const ccomplex_t IT_0464 = IT_0442 + IT_0463;
    const ccomplex_t IT_0465 = 12*IT_0371;
    const ccomplex_t IT_0466 = 12*conj(IT_0371);
    const ccomplex_t IT_0467 = conj(IT_0371)*IT_0415 + IT_0371*IT_0416 +
       IT_0396*IT_0436 + IT_0400*IT_0439 + conj(IT_0333)*IT_0450 + IT_0333
      *IT_0454 + conj(IT_0290)*IT_0457 + IT_0290*IT_0460 + conj(IT_0180)*IT_0462
       + IT_0180*IT_0464 + IT_0388*IT_0465 + IT_0392*IT_0466;
    const ccomplex_t IT_0468 = m_b*s_14*m_N_3;
    const ccomplex_t IT_0469 = 12*IT_0180;
    const ccomplex_t IT_0470 = IT_0457 + IT_0469;
    const ccomplex_t IT_0471 = 6*IT_0180;
    const ccomplex_t IT_0472 = 6*conj(IT_0180);
    const ccomplex_t IT_0473 = 12*conj(IT_0180);
    const ccomplex_t IT_0474 = IT_0388*(24*IT_0290 + IT_0436) + IT_0392*(24
      *conj(IT_0290) + IT_0439) + conj(IT_0065)*IT_0462 + IT_0065*IT_0464 +
       IT_0396*IT_0465 + IT_0400*IT_0466 + conj(IT_0333)*IT_0470 + conj(IT_0371)
      *IT_0471 + IT_0371*IT_0472 + IT_0333*(IT_0460 + IT_0473);
    const ccomplex_t IT_0475 = m_b*s_13*m_N_3;
    const ccomplex_t IT_0476 = 12*IT_0249;
    const ccomplex_t IT_0477 = 12*conj(IT_0249);
    const ccomplex_t IT_0478 = conj(IT_0249)*IT_0333 + IT_0249*conj(IT_0333) +
       0.166666666666667*conj(IT_0260)*IT_0437 + 0.166666666666667*IT_0260
      *IT_0440 + 4*IT_0380*(IT_0360 + 0.0416666666666667*IT_0441) + 4*IT_0384*
      (conj(IT_0360) + 0.0416666666666667*IT_0442) + 0.166666666666667*IT_0404
      *IT_0443 + 0.166666666666667*IT_0408*IT_0444 + 4*conj(IT_0371)*(IT_0408 +
       0.0416666666666667*IT_0476) + 4*IT_0371*(IT_0404 + 0.0416666666666667
      *IT_0477);
    const ccomplex_t IT_0479 = s_34*IT_0006;
    const ccomplex_t IT_0480 = IT_0431 + IT_0476;
    const ccomplex_t IT_0481 = conj(IT_0249)*IT_0415 + IT_0249*IT_0416 + conj
      (IT_0180)*IT_0417 + IT_0180*IT_0418 + conj(IT_0229)*(IT_0421 + IT_0422) +
       conj(IT_0273)*IT_0424 + IT_0229*(IT_0425 + IT_0426) + IT_0273*IT_0428 +
       IT_0404*IT_0447 + IT_0408*IT_0451 + conj(IT_0269)*IT_0455 + IT_0269
      *IT_0458 + IT_0380*IT_0469 + IT_0384*IT_0473 + IT_0392*(IT_0434 + IT_0477)
       + IT_0388*IT_0480;
    const ccomplex_t IT_0482 = s_34 + IT_0005;
    const ccomplex_t IT_0483 = s_12*IT_0482;
    const ccomplex_t IT_0484 = (-24)*IT_0321;
    const ccomplex_t IT_0485 = (-24)*conj(IT_0321);
    const ccomplex_t IT_0486 = -s_34;
    const ccomplex_t IT_0487 = IT_0005 + IT_0486;
    const ccomplex_t IT_0488 = s_12*IT_0487;
    const ccomplex_t IT_0489 = (-24)*IT_0198;
    const ccomplex_t IT_0490 = 24*IT_0200;
    const ccomplex_t IT_0491 = IT_0489 + IT_0490;
    const ccomplex_t IT_0492 = (-24)*IT_0200;
    const ccomplex_t IT_0493 = 24*IT_0198;
    const ccomplex_t IT_0494 = IT_0006*IT_0487;
    const ccomplex_t IT_0495 = 6*IT_0198;
    const ccomplex_t IT_0496 = (-6)*IT_0200;
    const ccomplex_t IT_0497 = IT_0495 + IT_0496;
    const ccomplex_t IT_0498 = 6*conj(IT_0198);
    const ccomplex_t IT_0499 = (-6)*conj(IT_0200);
    const ccomplex_t IT_0500 = IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = 6*IT_0200;
    const ccomplex_t IT_0502 = (-6)*IT_0198;
    const ccomplex_t IT_0503 = IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = 6*conj(IT_0200);
    const ccomplex_t IT_0505 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0506 = IT_0504 + IT_0505;
    const ccomplex_t IT_0507 = 6*conj(IT_0101);
    const ccomplex_t IT_0508 = (-6)*conj(IT_0099);
    const ccomplex_t IT_0509 = IT_0507 + IT_0508;
    const ccomplex_t IT_0510 = 6*IT_0101;
    const ccomplex_t IT_0511 = (-6)*IT_0101;
    const ccomplex_t IT_0512 = (-6)*conj(IT_0101);
    const ccomplex_t IT_0513 = IT_0101*((-24)*conj(IT_0198) + 24*conj(IT_0200)
      ) + IT_0099*(24*conj(IT_0198) + (-24)*conj(IT_0200) + 6*conj(IT_0290) + (
      -6)*conj(IT_0360)) + conj(IT_0101)*IT_0491 + 6*conj(IT_0099)*(IT_0290 +
       0.166666666666667*IT_0492 + 0.166666666666667*IT_0493) + conj(IT_0333)
      *IT_0497 + IT_0333*IT_0500 + conj(IT_0371)*IT_0503 + IT_0371*IT_0506 +
       IT_0360*IT_0509 + conj(IT_0360)*IT_0510 + conj(IT_0290)*IT_0511 + IT_0290
      *IT_0512;
    const ccomplex_t IT_0514 = 48*IT_0240;
    const ccomplex_t IT_0515 = 24*IT_0260;
    const ccomplex_t IT_0516 = IT_0514 + IT_0515;
    const ccomplex_t IT_0517 = 48*IT_0235;
    const ccomplex_t IT_0518 = 24*IT_0249;
    const ccomplex_t IT_0519 = IT_0517 + IT_0518;
    const ccomplex_t IT_0520 = 48*conj(IT_0240);
    const ccomplex_t IT_0521 = 48*conj(IT_0235);
    const ccomplex_t IT_0522 = 48*IT_0274;
    const ccomplex_t IT_0523 = conj(IT_0260)*IT_0424 + conj(IT_0360)*IT_0436 +
       conj(IT_0371)*IT_0465 + conj(IT_0249)*IT_0480 + conj(IT_0269)*(48*IT_0269
       + IT_0514 + IT_0515) + conj(IT_0240)*IT_0516 + conj(IT_0235)*IT_0519 +
       IT_0269*IT_0520 + IT_0274*IT_0521 + conj(IT_0274)*(IT_0519 + IT_0522);
    const ccomplex_t IT_0524 = conj(IT_0260)*IT_0415 + IT_0260*IT_0416 +
       IT_0388*IT_0424 + IT_0392*IT_0428 + IT_0380*IT_0447 + IT_0384*IT_0451 +
       IT_0404*(IT_0455 + IT_0456 + IT_0469) + conj(IT_0249)*IT_0471 + IT_0249
      *IT_0472 + IT_0408*(IT_0458 + IT_0459 + IT_0473) + IT_0396*IT_0476 +
       IT_0400*IT_0477;
    const ccomplex_t IT_0525 = IT_0005*IT_0006;
    const ccomplex_t IT_0526 = 48*IT_0238;
    const ccomplex_t IT_0527 = 24*IT_0065;
    const ccomplex_t IT_0528 = IT_0526 + IT_0527;
    const ccomplex_t IT_0529 = 48*conj(IT_0238);
    const ccomplex_t IT_0530 = 48*conj(IT_0229);
    const ccomplex_t IT_0531 = 48*IT_0229;
    const ccomplex_t IT_0532 = 48*IT_0273;
    const ccomplex_t IT_0533 = 48*IT_0270;
    const ccomplex_t IT_0534 = conj(IT_0249)*IT_0424 + IT_0249*IT_0428 + conj
      (IT_0260)*IT_0431 + conj(IT_0371)*IT_0437 + IT_0371*IT_0440 + conj(IT_0180
      )*IT_0450 + IT_0180*IT_0454 + conj(IT_0065)*IT_0457 + conj(IT_0333)
      *IT_0462 + IT_0333*IT_0464 + conj(IT_0274)*(48*IT_0269 + IT_0514 + IT_0515
      ) + conj(IT_0235)*IT_0516 + IT_0408*IT_0520 + IT_0269*IT_0521 + conj
      (IT_0269)*(IT_0517 + IT_0522) + conj(IT_0229)*IT_0528 + IT_0400*IT_0529 +
       IT_0270*IT_0530 + conj(IT_0270)*(IT_0531 + IT_0532) + conj(IT_0273)*
      (IT_0526 + IT_0527 + IT_0533);
    const ccomplex_t IT_0535 = conj(IT_0290)*IT_0441 + conj(IT_0333)*IT_0443 +
       conj(IT_0065)*IT_0450 + conj(IT_0180)*IT_0470 + conj(IT_0238)*IT_0528 +
       IT_0270*IT_0529 + IT_0273*IT_0530 + conj(IT_0229)*(24*IT_0180 + IT_0531) 
      + conj(IT_0273)*(24*IT_0180 + IT_0531 + IT_0532) + conj(IT_0270)*(IT_0526 
      + IT_0527 + IT_0533);
    const ccomplex_t IT_0536 = -IT_0280;
    const ccomplex_t IT_0537 = 0.166666666666667*conj(IT_0360);
    const ccomplex_t IT_0538 = IT_0279 + IT_0281 + IT_0536;
    const ccomplex_t IT_0539 = (-0.166666666666667)*IT_0333;
    const ccomplex_t IT_0540 = IT_0006*IT_0482;
    const ccomplex_t IT_0541 = IT_0323*conj(IT_0333) + IT_0366*IT_0371 +
       IT_0363*conj(IT_0371) + (-6)*conj(IT_0348)*(IT_0290 + 4*IT_0319 + (-4)
      *IT_0321 + (-0.166666666666667)*IT_0461) + (-6)*IT_0348*(conj(IT_0290) + 4
      *conj(IT_0319) + (-4)*conj(IT_0321) + (-0.166666666666667)*IT_0463) + 24
      *conj(IT_0346)*(IT_0319 + -1./4*IT_0360 + 0.0416666666666667*IT_0435 +
       0.0416666666666667*IT_0484) + 24*IT_0346*(conj(IT_0319) + -1./4*conj
      (IT_0360) + 0.0416666666666667*IT_0438 + 0.0416666666666667*IT_0485) + (-6
      )*IT_0326*IT_0539;
    const ccomplex_t IT_0542 = IT_0004*(conj(IT_0065)*IT_0103 + IT_0065
      *IT_0106 + conj(IT_0180)*IT_0202 + IT_0180*IT_0205) + (IT_0106*IT_0180 +
       IT_0103*conj(IT_0180) + conj(IT_0065)*IT_0202 + IT_0065*IT_0205)*IT_0209 
      + IT_0210*(IT_0106*(IT_0229 + IT_0235) + IT_0103*(conj(IT_0229) + conj
      (IT_0235)) + IT_0205*(IT_0238 + IT_0240) + IT_0202*(conj(IT_0238) + conj
      (IT_0240))) + (IT_0205*(IT_0229 + IT_0235) + IT_0202*(conj(IT_0229) + conj
      (IT_0235)) + IT_0106*(IT_0238 + IT_0240) + IT_0103*(conj(IT_0238) + conj
      (IT_0240)))*IT_0241 + IT_0242*(conj(IT_0249)*IT_0252 + IT_0249*IT_0255 +
       conj(IT_0260)*IT_0263 + IT_0260*IT_0266) + (IT_0255*IT_0260 + IT_0252
      *conj(IT_0260) + conj(IT_0249)*IT_0263 + IT_0249*IT_0266)*IT_0267 +
       IT_0268*(IT_0252*IT_0271 + IT_0255*IT_0272 + IT_0263*IT_0275 + IT_0266
      *IT_0276) + (IT_0263*IT_0271 + IT_0266*IT_0272 + IT_0252*IT_0275 + IT_0255
      *IT_0276)*IT_0277 + IT_0283*(conj(IT_0290)*IT_0323 + IT_0290*IT_0326 +
       conj(IT_0333)*IT_0350 + IT_0333*IT_0353) + IT_0355*(conj(IT_0360)*IT_0363
       + IT_0360*IT_0366 + conj(IT_0371)*IT_0374 + IT_0371*IT_0377) + IT_0379
      *IT_0412 + IT_0414*IT_0419 + IT_0420*IT_0445 + IT_0446*IT_0467 + IT_0468
      *IT_0474 + 6*IT_0475*IT_0478 + IT_0479*IT_0481 + IT_0483*(24*IT_0321*conj
      (IT_0321) + IT_0348*((-24)*conj(IT_0346) + 24*conj(IT_0348)) + IT_0346*(24
      *conj(IT_0346) + (-24)*conj(IT_0348)) + conj(IT_0319)*IT_0484 + IT_0319*
      (24*conj(IT_0319) + IT_0485)) + IT_0488*(IT_0099*((-24)*conj(IT_0099) + 24
      *conj(IT_0101)) + IT_0101*(24*conj(IT_0099) + (-24)*conj(IT_0101)) + conj
      (IT_0198)*IT_0491 + conj(IT_0200)*(IT_0492 + IT_0493)) + IT_0494*IT_0513 +
       IT_0281*IT_0523 + IT_0278*IT_0524 + IT_0525*IT_0534 + IT_0279*IT_0535 + 
      (IT_0278 + IT_0279 + IT_0282 + IT_0536)*(IT_0360*IT_0500 + 6*conj(IT_0371)
      *(IT_0099 + 0.166666666666667*IT_0511) + 6*IT_0371*(conj(IT_0099) +
       0.166666666666667*IT_0512) + 6*IT_0497*IT_0537) + 6*(conj(IT_0333)
      *IT_0371 + IT_0333*conj(IT_0371) + 0.166666666666667*IT_0360*IT_0438 +
       IT_0435*IT_0537)*IT_0538 + (IT_0278 + IT_0281 + IT_0354 + IT_0536)*(conj
      (IT_0290)*IT_0503 + IT_0290*IT_0506 + (-6)*conj(IT_0333)*(IT_0099 + (
      -0.166666666666667)*IT_0510) + (-6)*IT_0509*IT_0539) + IT_0540*IT_0541;
    return create_ccomplex_return(IT_0542);
}

