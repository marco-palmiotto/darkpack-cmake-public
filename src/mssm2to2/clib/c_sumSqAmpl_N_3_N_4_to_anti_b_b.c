#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_N_4_to_anti_b_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_N_4_to_anti_b_b(
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
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = (-0.5)*IT_0005 + (-0.5)*IT_0009 + 0.5*IT_0014 +
       0.5*IT_0017 + (-0.5)*IT_0019 + (-0.5)*IT_0021 + 0.5*IT_0023 + 0.5*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = m_b*e_em*IT_0011*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = pow(m_Z, 2);
    const ccomplex_t IT_0031 = pow(m_N_3, 2);
    const ccomplex_t IT_0032 = pow(m_N_4, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + IT_0030 + -IT_0031 + -IT_0032 
      + -reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = IT_0003*IT_0008;
    const ccomplex_t IT_0038 = IT_0004*IT_0007;
    const ccomplex_t IT_0039 = IT_0012*IT_0016;
    const ccomplex_t IT_0040 = IT_0013*IT_0015;
    const ccomplex_t IT_0041 = IT_0003*IT_0020;
    const ccomplex_t IT_0042 = IT_0007*IT_0018;
    const ccomplex_t IT_0043 = IT_0012*IT_0024;
    const ccomplex_t IT_0044 = IT_0015*IT_0022;
    const ccomplex_t IT_0045 = (-0.5)*IT_0037 + 0.5*IT_0038 + 0.5*IT_0039 + (
      -0.5)*IT_0040 + (-0.5)*IT_0041 + 0.5*IT_0042 + 0.5*IT_0043 + (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = cpow(IT_0000, -1);
    const ccomplex_t IT_0047 = m_b*e_em*IT_0006*IT_0011*IT_0027*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0031 + -IT_0032 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = -IT_0036 + -IT_0052;
    const ccomplex_t IT_0054 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0055 = IT_0003*IT_0054;
    const ccomplex_t IT_0056 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0057 = IT_0007*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0059 = IT_0012*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0061 = IT_0015*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0065 = IT_0007*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0067 = IT_0012*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0069 = IT_0015*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0055 + 0.5*IT_0057 + (-0.5)*IT_0059 + (
      -0.5)*IT_0061 + 0.5*IT_0063 + 0.5*IT_0065 + (-0.5)*IT_0067 + (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0035*IT_0070;
    const ccomplex_t IT_0072 = IT_0003*IT_0056;
    const ccomplex_t IT_0073 = IT_0007*IT_0054;
    const ccomplex_t IT_0074 = IT_0012*IT_0060;
    const ccomplex_t IT_0075 = IT_0015*IT_0058;
    const ccomplex_t IT_0076 = IT_0003*IT_0064;
    const ccomplex_t IT_0077 = IT_0007*IT_0062;
    const ccomplex_t IT_0078 = IT_0012*IT_0068;
    const ccomplex_t IT_0079 = IT_0015*IT_0066;
    const ccomplex_t IT_0080 = 0.5*IT_0072 + (-0.5)*IT_0073 + (-0.5)*IT_0074 +
       0.5*IT_0075 + 0.5*IT_0076 + (-0.5)*IT_0077 + (-0.5)*IT_0078 + 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0051*IT_0080;
    const ccomplex_t IT_0082 = -IT_0071 + -IT_0081;
    const ccomplex_t IT_0083 = pow(m_b, 2);
    const ccomplex_t IT_0084 = s_34 + IT_0083;
    const ccomplex_t IT_0085 = s_12*IT_0084;
    const ccomplex_t IT_0086 = cos(alpha);
    const ccomplex_t IT_0087 = IT_0002*IT_0086;
    const ccomplex_t IT_0088 = IT_0004*IT_0087;
    const ccomplex_t IT_0089 = sin(alpha);
    const ccomplex_t IT_0090 = IT_0011*IT_0089;
    const ccomplex_t IT_0091 = IT_0016*IT_0090;
    const ccomplex_t IT_0092 = IT_0018*IT_0087;
    const ccomplex_t IT_0093 = IT_0024*IT_0090;
    const ccomplex_t IT_0094 = IT_0088 + IT_0091 + IT_0092 + IT_0093;
    const ccomplex_t IT_0095 = IT_0002*IT_0089;
    const ccomplex_t IT_0096 = IT_0008*IT_0095;
    const ccomplex_t IT_0097 = IT_0011*IT_0086;
    const ccomplex_t IT_0098 = IT_0013*IT_0097;
    const ccomplex_t IT_0099 = IT_0020*IT_0095;
    const ccomplex_t IT_0100 = IT_0022*IT_0097;
    const ccomplex_t IT_0101 = -IT_0096 + -IT_0098 + -IT_0099 + -IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_b*e_em*IT_0011*IT_0027
      *IT_0046*IT_0086;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_12 + -IT_0031 + -IT_0032 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0104*IT_0109;
    const ccomplex_t IT_0111 = IT_0008*IT_0087;
    const ccomplex_t IT_0112 = IT_0004*IT_0095;
    const ccomplex_t IT_0113 = IT_0020*IT_0087;
    const ccomplex_t IT_0114 = IT_0018*IT_0095;
    const ccomplex_t IT_0115 = IT_0111 + IT_0112 + IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = IT_0016*IT_0097;
    const ccomplex_t IT_0117 = IT_0013*IT_0090;
    const ccomplex_t IT_0118 = IT_0024*IT_0097;
    const ccomplex_t IT_0119 = IT_0022*IT_0090;
    const ccomplex_t IT_0120 = -IT_0116 + -IT_0117 + -IT_0118 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0115 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*m_b*e_em*IT_0011*IT_0027
      *IT_0046*IT_0089;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0031 + -IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0123*IT_0128;
    const ccomplex_t IT_0130 = -IT_0110 + -IT_0129;
    const ccomplex_t IT_0131 = cpow((-2)*s_13 + IT_0031 + IT_0083 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0132 = N_B4*e_em*U_sb_10;
    const ccomplex_t IT_0133 = IT_0002*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = IT_0011*IT_0046;
    const ccomplex_t IT_0136 = m_b*N_d4*e_em*IT_0027*U_sb_00;
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0134 + 1.5*IT_0138);
    const ccomplex_t IT_0140 = (-0.333333333333333)*IT_0139;
    const ccomplex_t IT_0141 = m_b*N_d3*e_em*IT_0027*conj(U_sb_10);
    const ccomplex_t IT_0142 = IT_0135*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0145 = IT_0011*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0148 = IT_0002*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0143 + -IT_0146 +
       0.333333333333333*IT_0149);
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = (-3)*IT_0151;
    const ccomplex_t IT_0153 = 0.166666666666667*IT_0152;
    const ccomplex_t IT_0154 = IT_0140*IT_0153;
    const ccomplex_t IT_0155 = IT_0131*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = cpow((-2)*s_13 + IT_0031 + IT_0083 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0158 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0159 = IT_0002*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = m_b*N_d4*e_em*IT_0027*U_sb_01;
    const ccomplex_t IT_0162 = IT_0135*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + 1.5*IT_0163);
    const ccomplex_t IT_0165 = 0.666666666666667*IT_0164;
    const ccomplex_t IT_0166 = 1.5*IT_0165;
    const ccomplex_t IT_0167 = (-0.333333333333333)*IT_0166;
    const ccomplex_t IT_0168 = m_b*N_d3*e_em*IT_0027*conj(U_sb_11);
    const ccomplex_t IT_0169 = IT_0135*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = N_W3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0172 = IT_0011*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = N_B3*e_em*conj(U_sb_01);
    const ccomplex_t IT_0175 = IT_0002*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0170 + -IT_0173 +
       0.333333333333333*IT_0176);
    const ccomplex_t IT_0178 = 3*IT_0177;
    const ccomplex_t IT_0179 = 0.166666666666667*IT_0178;
    const ccomplex_t IT_0180 = IT_0167*IT_0179;
    const ccomplex_t IT_0181 = IT_0157*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = IT_0156 + IT_0182;
    const ccomplex_t IT_0184 = s_12*IT_0083;
    const ccomplex_t IT_0185 = s_14*s_23;
    const ccomplex_t IT_0186 = s_13*s_24;
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = s_12*s_34;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0184 + IT_0185 + IT_0187 + IT_0189;
    const ccomplex_t IT_0191 = cpow((-2)*s_23 + IT_0032 + IT_0083 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0192 = m_b*N_d4*e_em*IT_0027*conj(U_sb_10);
    const ccomplex_t IT_0193 = IT_0135*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0196 = IT_0011*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0199 = IT_0002*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*(IT_0194 + -IT_0197 +
       0.333333333333333*IT_0200);
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = (-3)*IT_0202;
    const ccomplex_t IT_0204 = 0.166666666666667*IT_0203;
    const ccomplex_t IT_0205 = N_B3*e_em*U_sb_10;
    const ccomplex_t IT_0206 = IT_0002*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = m_b*N_d3*e_em*IT_0027*U_sb_00;
    const ccomplex_t IT_0209 = IT_0135*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0207 + 1.5*IT_0210);
    const ccomplex_t IT_0212 = 0.666666666666667*IT_0211;
    const ccomplex_t IT_0213 = 1.5*IT_0212;
    const ccomplex_t IT_0214 = (-0.333333333333333)*IT_0213;
    const ccomplex_t IT_0215 = IT_0204*IT_0214;
    const ccomplex_t IT_0216 = IT_0191*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_23 + IT_0032 + IT_0083 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0219 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0220 = IT_0002*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = m_b*N_d3*e_em*IT_0027*U_sb_01;
    const ccomplex_t IT_0223 = IT_0135*IT_0222;
    const ccomplex_t IT_0224 = 1.4142135623731*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*(IT_0221 + 1.5*IT_0224);
    const ccomplex_t IT_0226 = 0.666666666666667*IT_0225;
    const ccomplex_t IT_0227 = 1.5*IT_0226;
    const ccomplex_t IT_0228 = 0.666666666666667*IT_0227;
    const ccomplex_t IT_0229 = 1.5*IT_0228;
    const ccomplex_t IT_0230 = 0.666666666666667*IT_0229;
    const ccomplex_t IT_0231 = 1.5*IT_0230;
    const ccomplex_t IT_0232 = (-0.333333333333333)*IT_0231;
    const ccomplex_t IT_0233 = N_W4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0234 = IT_0011*IT_0233;
    const ccomplex_t IT_0235 = 1.4142135623731*IT_0234;
    const ccomplex_t IT_0236 = m_b*N_d4*e_em*IT_0027*conj(U_sb_11);
    const ccomplex_t IT_0237 = IT_0135*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = N_B4*e_em*conj(U_sb_01);
    const ccomplex_t IT_0240 = IT_0002*IT_0239;
    const ccomplex_t IT_0241 = 1.4142135623731*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*(IT_0235 + -IT_0238 + (
      -0.333333333333333)*IT_0241);
    const ccomplex_t IT_0243 = (-3)*IT_0242;
    const ccomplex_t IT_0244 = 0.166666666666667*IT_0243;
    const ccomplex_t IT_0245 = IT_0232*IT_0244;
    const ccomplex_t IT_0246 = IT_0218*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = -IT_0217 + -IT_0247;
    const ccomplex_t IT_0249 = -IT_0185;
    const ccomplex_t IT_0250 = IT_0184 + IT_0186 + IT_0189 + IT_0249;
    const ccomplex_t IT_0251 = -s_34;
    const ccomplex_t IT_0252 = IT_0083 + IT_0251;
    const ccomplex_t IT_0253 = s_12*IT_0252;
    const ccomplex_t IT_0254 = IT_0054*IT_0087;
    const ccomplex_t IT_0255 = IT_0060*IT_0090;
    const ccomplex_t IT_0256 = IT_0062*IT_0087;
    const ccomplex_t IT_0257 = IT_0068*IT_0090;
    const ccomplex_t IT_0258 = IT_0254 + IT_0255 + IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = IT_0056*IT_0095;
    const ccomplex_t IT_0260 = IT_0058*IT_0097;
    const ccomplex_t IT_0261 = IT_0064*IT_0095;
    const ccomplex_t IT_0262 = IT_0066*IT_0097;
    const ccomplex_t IT_0263 = -IT_0259 + -IT_0260 + -IT_0261 + -IT_0262;
    const ccomplex_t IT_0264 = IT_0258 + IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = 0.5*IT_0265;
    const ccomplex_t IT_0267 = IT_0109*IT_0266;
    const ccomplex_t IT_0268 = IT_0056*IT_0087;
    const ccomplex_t IT_0269 = IT_0054*IT_0095;
    const ccomplex_t IT_0270 = IT_0064*IT_0087;
    const ccomplex_t IT_0271 = IT_0062*IT_0095;
    const ccomplex_t IT_0272 = IT_0268 + IT_0269 + IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = IT_0060*IT_0097;
    const ccomplex_t IT_0274 = IT_0058*IT_0090;
    const ccomplex_t IT_0275 = IT_0068*IT_0097;
    const ccomplex_t IT_0276 = IT_0066*IT_0090;
    const ccomplex_t IT_0277 = -IT_0273 + -IT_0274 + -IT_0275 + -IT_0276;
    const ccomplex_t IT_0278 = IT_0272 + IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = (-0.5)*IT_0279;
    const ccomplex_t IT_0281 = IT_0128*IT_0280;
    const ccomplex_t IT_0282 = -IT_0267 + -IT_0281;
    const ccomplex_t IT_0283 = conj(N_B4)*e_em*U_sb_00;
    const ccomplex_t IT_0284 = IT_0002*IT_0283;
    const ccomplex_t IT_0285 = 1.4142135623731*IT_0284;
    const ccomplex_t IT_0286 = m_b*conj(N_d4)*e_em*IT_0027*U_sb_10;
    const ccomplex_t IT_0287 = IT_0135*IT_0286;
    const ccomplex_t IT_0288 = 1.4142135623731*IT_0287;
    const ccomplex_t IT_0289 = conj(N_W4)*e_em*U_sb_00;
    const ccomplex_t IT_0290 = IT_0011*IT_0289;
    const ccomplex_t IT_0291 = 1.4142135623731*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*(IT_0285 + 3*IT_0288 + (-3)
      *IT_0291);
    const ccomplex_t IT_0293 = 0.333333333333333*IT_0292;
    const ccomplex_t IT_0294 = 3*IT_0293;
    const ccomplex_t IT_0295 = 0.166666666666667*IT_0294;
    const ccomplex_t IT_0296 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0297 = IT_0002*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = m_b*conj(N_d3)*e_em*IT_0027*conj(U_sb_00);
    const ccomplex_t IT_0300 = IT_0135*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*(IT_0298 + 1.5*IT_0301);
    const ccomplex_t IT_0303 = 0.666666666666667*IT_0302;
    const ccomplex_t IT_0304 = 1.5*IT_0303;
    const ccomplex_t IT_0305 = (-0.333333333333333)*IT_0304;
    const ccomplex_t IT_0306 = IT_0295*IT_0305;
    const ccomplex_t IT_0307 = IT_0131*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = conj(N_B3)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0310 = IT_0002*IT_0309;
    const ccomplex_t IT_0311 = 1.4142135623731*IT_0310;
    const ccomplex_t IT_0312 = m_b*conj(N_d3)*e_em*IT_0027*conj(U_sb_01);
    const ccomplex_t IT_0313 = IT_0135*IT_0312;
    const ccomplex_t IT_0314 = 1.4142135623731*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*(IT_0311 + 1.5*IT_0314);
    const ccomplex_t IT_0316 = (-0.333333333333333)*IT_0315;
    const ccomplex_t IT_0317 = m_b*conj(N_d4)*e_em*IT_0027*U_sb_11;
    const ccomplex_t IT_0318 = IT_0135*IT_0317;
    const ccomplex_t IT_0319 = 1.4142135623731*IT_0318;
    const ccomplex_t IT_0320 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0321 = IT_0002*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0324 = IT_0011*IT_0323;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*(IT_0319 + 0.333333333333333
      *IT_0322 + -IT_0325);
    const ccomplex_t IT_0327 = 3*IT_0326;
    const ccomplex_t IT_0328 = 0.166666666666667*IT_0327;
    const ccomplex_t IT_0329 = IT_0316*IT_0328;
    const ccomplex_t IT_0330 = IT_0157*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = IT_0308 + IT_0331;
    const ccomplex_t IT_0333 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0334 = IT_0002*IT_0333;
    const ccomplex_t IT_0335 = 1.4142135623731*IT_0334;
    const ccomplex_t IT_0336 = m_b*conj(N_d4)*e_em*IT_0027*conj(U_sb_00);
    const ccomplex_t IT_0337 = IT_0135*IT_0336;
    const ccomplex_t IT_0338 = 1.4142135623731*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*(IT_0335 + 1.5*IT_0338);
    const ccomplex_t IT_0340 = 0.666666666666667*IT_0339;
    const ccomplex_t IT_0341 = 1.5*IT_0340;
    const ccomplex_t IT_0342 = (-0.333333333333333)*IT_0341;
    const ccomplex_t IT_0343 = m_b*conj(N_d3)*e_em*IT_0027*U_sb_10;
    const ccomplex_t IT_0344 = IT_0135*IT_0343;
    const ccomplex_t IT_0345 = 1.4142135623731*IT_0344;
    const ccomplex_t IT_0346 = conj(N_B3)*e_em*U_sb_00;
    const ccomplex_t IT_0347 = IT_0002*IT_0346;
    const ccomplex_t IT_0348 = 1.4142135623731*IT_0347;
    const ccomplex_t IT_0349 = conj(N_W3)*e_em*U_sb_00;
    const ccomplex_t IT_0350 = IT_0011*IT_0349;
    const ccomplex_t IT_0351 = 1.4142135623731*IT_0350;
    const ccomplex_t IT_0352 = (0 + _Complex_I*1)*(IT_0345 + 0.333333333333333
      *IT_0348 + -IT_0351);
    const ccomplex_t IT_0353 = -IT_0352;
    const ccomplex_t IT_0354 = -IT_0353;
    const ccomplex_t IT_0355 = -IT_0354;
    const ccomplex_t IT_0356 = -IT_0355;
    const ccomplex_t IT_0357 = 3*IT_0356;
    const ccomplex_t IT_0358 = (-0.333333333333333)*IT_0357;
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = 3*IT_0359;
    const ccomplex_t IT_0361 = (-0.333333333333333)*IT_0360;
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = -IT_0363;
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = -IT_0365;
    const ccomplex_t IT_0367 = 3*IT_0366;
    const ccomplex_t IT_0368 = (-0.333333333333333)*IT_0367;
    const ccomplex_t IT_0369 = (-3)*IT_0368;
    const ccomplex_t IT_0370 = 0.333333333333333*IT_0369;
    const ccomplex_t IT_0371 = 3*IT_0370;
    const ccomplex_t IT_0372 = 0.166666666666667*IT_0371;
    const ccomplex_t IT_0373 = IT_0342*IT_0372;
    const ccomplex_t IT_0374 = IT_0191*IT_0373;
    const ccomplex_t IT_0375 = (0 + _Complex_I*1)*IT_0374;
    const ccomplex_t IT_0376 = m_b*conj(N_d3)*e_em*IT_0027*U_sb_11;
    const ccomplex_t IT_0377 = IT_0135*IT_0376;
    const ccomplex_t IT_0378 = 1.4142135623731*IT_0377;
    const ccomplex_t IT_0379 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0380 = IT_0002*IT_0379;
    const ccomplex_t IT_0381 = 1.4142135623731*IT_0380;
    const ccomplex_t IT_0382 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0383 = IT_0011*IT_0382;
    const ccomplex_t IT_0384 = 1.4142135623731*IT_0383;
    const ccomplex_t IT_0385 = (0 + _Complex_I*1)*(IT_0378 + 0.333333333333333
      *IT_0381 + -IT_0384);
    const ccomplex_t IT_0386 = 3*IT_0385;
    const ccomplex_t IT_0387 = 0.166666666666667*IT_0386;
    const ccomplex_t IT_0388 = conj(N_B4)*e_em*conj(U_sb_11);
    const ccomplex_t IT_0389 = IT_0002*IT_0388;
    const ccomplex_t IT_0390 = 1.4142135623731*IT_0389;
    const ccomplex_t IT_0391 = m_b*conj(N_d4)*e_em*IT_0027*conj(U_sb_01);
    const ccomplex_t IT_0392 = IT_0135*IT_0391;
    const ccomplex_t IT_0393 = 1.4142135623731*IT_0392;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*(IT_0390 + 1.5*IT_0393);
    const ccomplex_t IT_0395 = (-0.333333333333333)*IT_0394;
    const ccomplex_t IT_0396 = IT_0387*IT_0395;
    const ccomplex_t IT_0397 = IT_0218*IT_0396;
    const ccomplex_t IT_0398 = (0 + _Complex_I*1)*IT_0397;
    const ccomplex_t IT_0399 = -IT_0375 + -IT_0398;
    const ccomplex_t IT_0400 = m_N_3*m_N_4;
    const ccomplex_t IT_0401 = -IT_0083;
    const ccomplex_t IT_0402 = s_34 + IT_0401;
    const ccomplex_t IT_0403 = IT_0400*IT_0402;
    const ccomplex_t IT_0404 = (-8)*IT_0403;
    const ccomplex_t IT_0405 = 3*IT_0130;
    const ccomplex_t IT_0406 = 3*conj(IT_0130);
    const ccomplex_t IT_0407 = m_b*m_N_4;
    const ccomplex_t IT_0408 = -s_14;
    const ccomplex_t IT_0409 = s_13 + IT_0408;
    const ccomplex_t IT_0410 = IT_0407*IT_0409;
    const ccomplex_t IT_0411 = 2*IT_0410;
    const ccomplex_t IT_0412 = IT_0153*IT_0295;
    const ccomplex_t IT_0413 = IT_0131*IT_0412;
    const ccomplex_t IT_0414 = (0 + _Complex_I*1)*IT_0413;
    const ccomplex_t IT_0415 = IT_0179*IT_0328;
    const ccomplex_t IT_0416 = IT_0157*IT_0415;
    const ccomplex_t IT_0417 = (0 + _Complex_I*1)*IT_0416;
    const ccomplex_t IT_0418 = -IT_0414 + -IT_0417;
    const ccomplex_t IT_0419 = IT_0140*IT_0305;
    const ccomplex_t IT_0420 = IT_0131*IT_0419;
    const ccomplex_t IT_0421 = (0 + _Complex_I*1)*IT_0420;
    const ccomplex_t IT_0422 = IT_0167*IT_0316;
    const ccomplex_t IT_0423 = IT_0157*IT_0422;
    const ccomplex_t IT_0424 = (0 + _Complex_I*1)*IT_0423;
    const ccomplex_t IT_0425 = -IT_0421 + -IT_0424;
    const ccomplex_t IT_0426 = 3*IT_0282;
    const ccomplex_t IT_0427 = 3*conj(IT_0282);
    const ccomplex_t IT_0428 = m_b*m_N_3;
    const ccomplex_t IT_0429 = -s_24;
    const ccomplex_t IT_0430 = s_23 + IT_0429;
    const ccomplex_t IT_0431 = IT_0428*IT_0430;
    const ccomplex_t IT_0432 = (-2)*IT_0431;
    const ccomplex_t IT_0433 = IT_0185 + IT_0186 + IT_0189;
    const ccomplex_t IT_0434 = 6*IT_0183;
    const ccomplex_t IT_0435 = 6*conj(IT_0183);
    const ccomplex_t IT_0436 = m_b*s_13*m_N_4;
    const ccomplex_t IT_0437 = 6*IT_0248;
    const ccomplex_t IT_0438 = 12*IT_0183;
    const ccomplex_t IT_0439 = IT_0437 + IT_0438;
    const ccomplex_t IT_0440 = 6*conj(IT_0248);
    const ccomplex_t IT_0441 = 12*conj(IT_0183);
    const ccomplex_t IT_0442 = IT_0440 + IT_0441;
    const ccomplex_t IT_0443 = IT_0002*IT_0010;
    const ccomplex_t IT_0444 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0445 = IT_0443*IT_0444;
    const ccomplex_t IT_0446 = IT_0001*IT_0011;
    const ccomplex_t IT_0447 = IT_0444*IT_0446;
    const ccomplex_t IT_0448 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0449 = IT_0443*IT_0448;
    const ccomplex_t IT_0450 = IT_0446*IT_0448;
    const ccomplex_t IT_0451 = (0 + _Complex_I*1)*(IT_0445 + IT_0447 + 
      -IT_0449 + -IT_0450);
    const ccomplex_t IT_0452 = (-0.5)*IT_0451;
    const ccomplex_t IT_0453 = e_em*IT_0443;
    const ccomplex_t IT_0454 = e_em*IT_0446;
    const ccomplex_t IT_0455 = (0 + _Complex_I*1)*(IT_0453 + 3*IT_0454);
    const ccomplex_t IT_0456 = (-0.166666666666667)*IT_0455;
    const ccomplex_t IT_0457 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0031 + -IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0458 = IT_0456*IT_0457;
    const ccomplex_t IT_0459 = (0 + _Complex_I*1)*IT_0458;
    const ccomplex_t IT_0460 = IT_0452*IT_0459;
    const ccomplex_t IT_0461 = -IT_0460;
    const ccomplex_t IT_0462 = 12*IT_0248;
    const ccomplex_t IT_0463 = 24*IT_0183;
    const ccomplex_t IT_0464 = IT_0462 + IT_0463;
    const ccomplex_t IT_0465 = 12*conj(IT_0248);
    const ccomplex_t IT_0466 = 24*conj(IT_0183);
    const ccomplex_t IT_0467 = IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0469 = IT_0443*IT_0468;
    const ccomplex_t IT_0470 = IT_0446*IT_0468;
    const ccomplex_t IT_0471 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0472 = IT_0443*IT_0471;
    const ccomplex_t IT_0473 = IT_0446*IT_0471;
    const ccomplex_t IT_0474 = (0 + _Complex_I*1)*(IT_0469 + IT_0470 + 
      -IT_0472 + -IT_0473);
    const ccomplex_t IT_0475 = (-0.5)*IT_0474;
    const ccomplex_t IT_0476 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0010;
    const ccomplex_t IT_0477 = 0.333333333333333*IT_0476;
    const ccomplex_t IT_0478 = IT_0457*IT_0477;
    const ccomplex_t IT_0479 = (0 + _Complex_I*1)*IT_0478;
    const ccomplex_t IT_0480 = IT_0475*IT_0479;
    const ccomplex_t IT_0481 = -IT_0480;
    const ccomplex_t IT_0482 = 12*IT_0399;
    const ccomplex_t IT_0483 = 24*IT_0332;
    const ccomplex_t IT_0484 = IT_0482 + IT_0483;
    const ccomplex_t IT_0485 = 12*conj(IT_0399);
    const ccomplex_t IT_0486 = 24*conj(IT_0332);
    const ccomplex_t IT_0487 = IT_0485 + IT_0486;
    const ccomplex_t IT_0488 = 6*IT_0425;
    const ccomplex_t IT_0489 = 6*conj(IT_0425);
    const ccomplex_t IT_0490 = 12*IT_0425;
    const ccomplex_t IT_0491 = 12*conj(IT_0425);
    const ccomplex_t IT_0492 = IT_0436*(conj(IT_0418)*IT_0439 + IT_0418
      *IT_0442 + conj(IT_0461)*IT_0464 + IT_0461*IT_0467 + conj(IT_0481)*IT_0484
       + IT_0481*IT_0487 + conj(IT_0399)*IT_0488 + IT_0399*IT_0489 + conj
      (IT_0332)*IT_0490 + IT_0332*IT_0491);
    const ccomplex_t IT_0493 = m_b*s_14*m_N_4;
    const ccomplex_t IT_0494 = IT_0214*IT_0342;
    const ccomplex_t IT_0495 = IT_0191*IT_0494;
    const ccomplex_t IT_0496 = (0 + _Complex_I*1)*IT_0495;
    const ccomplex_t IT_0497 = IT_0232*IT_0395;
    const ccomplex_t IT_0498 = IT_0218*IT_0497;
    const ccomplex_t IT_0499 = (0 + _Complex_I*1)*IT_0498;
    const ccomplex_t IT_0500 = IT_0496 + IT_0499;
    const ccomplex_t IT_0501 = IT_0434 + IT_0462;
    const ccomplex_t IT_0502 = IT_0435 + IT_0465;
    const ccomplex_t IT_0503 = IT_0452*IT_0479;
    const ccomplex_t IT_0504 = -IT_0503;
    const ccomplex_t IT_0505 = 24*IT_0248;
    const ccomplex_t IT_0506 = IT_0438 + IT_0505;
    const ccomplex_t IT_0507 = 24*conj(IT_0248);
    const ccomplex_t IT_0508 = IT_0441 + IT_0507;
    const ccomplex_t IT_0509 = IT_0459*IT_0475;
    const ccomplex_t IT_0510 = -IT_0509;
    const ccomplex_t IT_0511 = 12*IT_0332;
    const ccomplex_t IT_0512 = 24*IT_0399;
    const ccomplex_t IT_0513 = IT_0511 + IT_0512;
    const ccomplex_t IT_0514 = 12*conj(IT_0332);
    const ccomplex_t IT_0515 = 24*conj(IT_0399);
    const ccomplex_t IT_0516 = IT_0514 + IT_0515;
    const ccomplex_t IT_0517 = IT_0204*IT_0372;
    const ccomplex_t IT_0518 = IT_0191*IT_0517;
    const ccomplex_t IT_0519 = (0 + _Complex_I*1)*IT_0518;
    const ccomplex_t IT_0520 = IT_0244*IT_0387;
    const ccomplex_t IT_0521 = IT_0218*IT_0520;
    const ccomplex_t IT_0522 = (0 + _Complex_I*1)*IT_0521;
    const ccomplex_t IT_0523 = IT_0519 + IT_0522;
    const ccomplex_t IT_0524 = 6*IT_0523;
    const ccomplex_t IT_0525 = 6*conj(IT_0523);
    const ccomplex_t IT_0526 = 12*IT_0523;
    const ccomplex_t IT_0527 = 12*conj(IT_0523);
    const ccomplex_t IT_0528 = IT_0493*(conj(IT_0500)*IT_0501 + IT_0500
      *IT_0502 + conj(IT_0504)*IT_0506 + IT_0504*IT_0508 + conj(IT_0510)*IT_0513
       + IT_0510*IT_0516 + conj(IT_0332)*IT_0524 + IT_0332*IT_0525 + conj
      (IT_0399)*IT_0526 + IT_0399*IT_0527);
    const ccomplex_t IT_0529 = m_b*s_24*m_N_3;
    const ccomplex_t IT_0530 = 6*IT_0418;
    const ccomplex_t IT_0531 = 6*conj(IT_0418);
    const ccomplex_t IT_0532 = 12*IT_0418;
    const ccomplex_t IT_0533 = 12*conj(IT_0418);
    const ccomplex_t IT_0534 = IT_0529*(conj(IT_0425)*IT_0439 + IT_0425
      *IT_0442 + IT_0467*IT_0481 + IT_0464*conj(IT_0481) + conj(IT_0461)*IT_0484
       + IT_0461*IT_0487 + conj(IT_0399)*IT_0530 + IT_0399*IT_0531 + conj
      (IT_0332)*IT_0532 + IT_0332*IT_0533);
    const ccomplex_t IT_0535 = 24*IT_0481;
    const ccomplex_t IT_0536 = IT_0490 + IT_0535;
    const ccomplex_t IT_0537 = 24*conj(IT_0481);
    const ccomplex_t IT_0538 = IT_0491 + IT_0537;
    const ccomplex_t IT_0539 = 24*IT_0461;
    const ccomplex_t IT_0540 = IT_0532 + IT_0539;
    const ccomplex_t IT_0541 = 24*conj(IT_0461);
    const ccomplex_t IT_0542 = IT_0533 + IT_0541;
    const ccomplex_t IT_0543 = 12*IT_0500;
    const ccomplex_t IT_0544 = 12*conj(IT_0500);
    const ccomplex_t IT_0545 = IT_0489*IT_0523 + IT_0488*conj(IT_0523) + conj
      (IT_0481)*IT_0526 + IT_0481*IT_0527 + conj(IT_0500)*IT_0530 + IT_0500
      *IT_0531 + conj(IT_0510)*IT_0536 + IT_0510*IT_0538 + conj(IT_0504)*IT_0540
       + IT_0504*IT_0542 + conj(IT_0461)*IT_0543 + IT_0461*IT_0544;
    const ccomplex_t IT_0546 = s_34*m_N_3*m_N_4;
    const ccomplex_t IT_0547 = IT_0489*IT_0500 + IT_0488*conj(IT_0500) + conj
      (IT_0418)*IT_0524 + IT_0418*IT_0525 + conj(IT_0461)*IT_0526 + IT_0461
      *IT_0527 + conj(IT_0504)*IT_0536 + IT_0504*IT_0538 + conj(IT_0510)*IT_0540
       + IT_0510*IT_0542 + conj(IT_0481)*IT_0543 + IT_0481*IT_0544;
    const ccomplex_t IT_0548 = m_b*s_23*m_N_3;
    const ccomplex_t IT_0549 = conj(IT_0332)*IT_0500 + IT_0332*conj(IT_0500) +
       0.166666666666667*IT_0508*IT_0510 + 0.166666666666667*IT_0506*conj
      (IT_0510) + 0.166666666666667*conj(IT_0504)*IT_0513 + 0.166666666666667
      *IT_0504*IT_0516 + 0.166666666666667*IT_0502*IT_0523 + 0.166666666666667
      *IT_0501*conj(IT_0523) + 0.166666666666667*conj(IT_0399)*IT_0543 +
       0.166666666666667*IT_0399*IT_0544;
    const ccomplex_t IT_0550 = 24*IT_0500;
    const ccomplex_t IT_0551 = 48*IT_0504;
    const ccomplex_t IT_0552 = IT_0550 + IT_0551;
    const ccomplex_t IT_0553 = 24*IT_0523;
    const ccomplex_t IT_0554 = 24*conj(IT_0523);
    const ccomplex_t IT_0555 = 24*conj(IT_0500);
    const ccomplex_t IT_0556 = s_13 + s_14;
    const ccomplex_t IT_0557 = IT_0407*IT_0556;
    const ccomplex_t IT_0558 = 6*IT_0053;
    const ccomplex_t IT_0559 = 6*conj(IT_0053);
    const ccomplex_t IT_0560 = 12*IT_0053;
    const ccomplex_t IT_0561 = 12*conj(IT_0053);
    const ccomplex_t IT_0562 = 12*IT_0082;
    const ccomplex_t IT_0563 = 12*conj(IT_0082);
    const ccomplex_t IT_0564 = (-12)*conj(IT_0053);
    const ccomplex_t IT_0565 = (-12)*IT_0082;
    const ccomplex_t IT_0566 = (-12)*conj(IT_0082);
    const ccomplex_t IT_0567 = (-6)*IT_0053;
    const ccomplex_t IT_0568 = (-6)*conj(IT_0053);
    const ccomplex_t IT_0569 = (-12)*IT_0053*conj(IT_0461) + (-6)*conj(IT_0082
      )*((-0.166666666666667)*IT_0488 + IT_0523) + (-6)*IT_0082*((
      -0.166666666666667)*IT_0489 + conj(IT_0523)) + conj(IT_0500)*IT_0558 +
       IT_0500*IT_0559 + conj(IT_0504)*IT_0560 + IT_0504*IT_0561 + conj(IT_0481)
      *IT_0562 + IT_0481*IT_0563 + IT_0461*IT_0564 + conj(IT_0510)*IT_0565 +
       IT_0510*IT_0566 + conj(IT_0418)*IT_0567 + IT_0418*IT_0568;
    const ccomplex_t IT_0570 = IT_0184 + IT_0186 + IT_0188 + IT_0249;
    const ccomplex_t IT_0571 = 6*IT_0082;
    const ccomplex_t IT_0572 = 6*conj(IT_0082);
    const ccomplex_t IT_0573 = IT_0084*IT_0400;
    const ccomplex_t IT_0574 = 4*conj(IT_0082)*(IT_0053 + 1./4*IT_0248 + (
      -0.0416666666666667)*IT_0434) + 4*IT_0082*(conj(IT_0053) + 1./4*conj
      (IT_0248) + (-0.0416666666666667)*IT_0435) + (-0.166666666666667)*conj
      (IT_0399)*IT_0558 + (-0.166666666666667)*IT_0399*IT_0559 + (
      -0.166666666666667)*conj(IT_0332)*IT_0567 + (-0.166666666666667)*IT_0332
      *IT_0568;
    const ccomplex_t IT_0575 = s_23 + s_24;
    const ccomplex_t IT_0576 = IT_0428*IT_0575;
    const ccomplex_t IT_0577 = (-6)*IT_0082;
    const ccomplex_t IT_0578 = (-6)*conj(IT_0082);
    const ccomplex_t IT_0579 = (-12)*IT_0053*conj(IT_0481) + conj(IT_0523)
      *IT_0558 + IT_0523*IT_0559 + conj(IT_0510)*IT_0560 + IT_0510*IT_0561 +
       conj(IT_0461)*IT_0562 + IT_0461*IT_0563 + IT_0481*IT_0564 + conj(IT_0504)
      *IT_0565 + IT_0504*IT_0566 + conj(IT_0425)*IT_0567 + IT_0425*IT_0568 +
       conj(IT_0418)*IT_0571 + IT_0418*IT_0572 + conj(IT_0500)*IT_0577 + IT_0500
      *IT_0578;
    const ccomplex_t IT_0580 = IT_0184 + IT_0185 + IT_0187 + IT_0188;
    const ccomplex_t IT_0581 = (-2)*IT_0410;
    const ccomplex_t IT_0582 = (-3)*IT_0130;
    const ccomplex_t IT_0583 = (-3)*conj(IT_0130);
    const ccomplex_t IT_0584 = (-3)*IT_0282;
    const ccomplex_t IT_0585 = (-3)*conj(IT_0282);
    const ccomplex_t IT_0586 = 2*IT_0431;
    const ccomplex_t IT_0587 = (-2)*IT_0403;
    const ccomplex_t IT_0588 = 4*IT_0431;
    const ccomplex_t IT_0589 = conj(IT_0481) + conj(IT_0510);
    const ccomplex_t IT_0590 = IT_0481 + IT_0510;
    const ccomplex_t IT_0591 = conj(IT_0461) + conj(IT_0504);
    const ccomplex_t IT_0592 = IT_0461 + IT_0504;
    const ccomplex_t IT_0593 = (-4)*IT_0410;
    const ccomplex_t IT_0594 = m_N_3*m_N_4*IT_0083;
    const ccomplex_t IT_0595 = 24*conj(IT_0425);
    const ccomplex_t IT_0596 = 24*IT_0418;
    const ccomplex_t IT_0597 = 24*conj(IT_0418);
    const ccomplex_t IT_0598 = 48*IT_0481;
    const ccomplex_t IT_0599 = conj(IT_0332)*IT_0439 + IT_0332*IT_0442 + conj
      (IT_0418)*IT_0490 + IT_0418*IT_0491 + conj(IT_0399)*IT_0501 + IT_0399
      *IT_0502 + conj(IT_0500)*IT_0526 + IT_0500*IT_0527 + conj(IT_0510)*IT_0552
       + conj(IT_0504)*(48*IT_0510 + IT_0553) + IT_0504*IT_0554 + IT_0510
      *IT_0555 + IT_0461*(48*conj(IT_0481) + IT_0595) + conj(IT_0481)*IT_0596 +
       IT_0481*IT_0597 + conj(IT_0461)*(24*IT_0425 + IT_0598);
    const ccomplex_t IT_0600 = IT_0186*(conj(IT_0183)*IT_0438 + conj(IT_0425)
      *IT_0490 + conj(IT_0332)*IT_0511 + conj(IT_0418)*IT_0532 + IT_0481*IT_0595
       + conj(IT_0461)*IT_0596 + IT_0461*(48*conj(IT_0461) + IT_0597) + conj
      (IT_0481)*(24*IT_0425 + IT_0598));
    const ccomplex_t IT_0601 = (24*IT_0053*conj(IT_0053) + 24*IT_0082*conj
      (IT_0082))*IT_0085 + (-6)*conj(IT_0130)*(IT_0183*IT_0190 + -IT_0248
      *IT_0250) + (-6)*IT_0130*(conj(IT_0183)*IT_0190 + -conj(IT_0248)*IT_0250 +
       4*conj(IT_0130)*IT_0253) + (-6)*conj(IT_0282)*(IT_0190*IT_0332 + -IT_0250
      *IT_0399 + (-0.166666666666667)*IT_0404*IT_0405) + IT_0282*((-24)*IT_0253
      *conj(IT_0282) + (-6)*IT_0190*conj(IT_0332) + 6*IT_0250*conj(IT_0399) +
       IT_0404*IT_0406) + IT_0411*(IT_0406*IT_0418 + IT_0405*conj(IT_0418) +
       conj(IT_0425)*IT_0426 + IT_0425*IT_0427) + (IT_0406*IT_0425 + IT_0405
      *conj(IT_0425) + conj(IT_0418)*IT_0426 + IT_0418*IT_0427)*IT_0432 + 6
      *IT_0433*(conj(IT_0332)*IT_0399 + IT_0332*conj(IT_0399) +
       0.166666666666667*conj(IT_0248)*IT_0434 + 0.166666666666667*IT_0248
      *IT_0435) + IT_0492 + IT_0528 + IT_0534 + IT_0184*IT_0545 + IT_0546
      *IT_0547 + 6*IT_0548*IT_0549 + IT_0185*(conj(IT_0248)*IT_0462 + conj
      (IT_0399)*IT_0482 + conj(IT_0523)*IT_0526 + conj(IT_0500)*IT_0543 + conj
      (IT_0504)*IT_0552 + conj(IT_0510)*IT_0553 + IT_0510*(48*conj(IT_0510) +
       IT_0554) + IT_0504*IT_0555) + IT_0557*IT_0569 + IT_0570*(conj(IT_0183)
      *IT_0567 + IT_0183*IT_0568 + conj(IT_0332)*IT_0571 + IT_0332*IT_0572) + (
      -6)*IT_0573*IT_0574 + IT_0576*IT_0579 + (conj(IT_0248)*IT_0558 + IT_0248
      *IT_0559 + conj(IT_0399)*IT_0577 + IT_0399*IT_0578)*IT_0580 + IT_0581*
      (conj(IT_0500)*IT_0582 + IT_0500*IT_0583 + conj(IT_0523)*IT_0584 + IT_0523
      *IT_0585) + (conj(IT_0523)*IT_0582 + IT_0523*IT_0583 + conj(IT_0500)
      *IT_0584 + IT_0500*IT_0585)*IT_0586 + (conj(IT_0332)*IT_0405 + IT_0332
      *IT_0406 + conj(IT_0183)*IT_0426 + IT_0183*IT_0427 + conj(IT_0399)*IT_0582
       + IT_0399*IT_0583 + conj(IT_0248)*IT_0584 + IT_0248*IT_0585)*IT_0587 +
       IT_0588*(IT_0582*IT_0589 + IT_0583*IT_0590 + IT_0584*IT_0591 + IT_0585
      *IT_0592) + (IT_0584*IT_0589 + IT_0585*IT_0590 + IT_0582*IT_0591 + IT_0583
      *IT_0592)*IT_0593 + IT_0594*IT_0599 + IT_0600;
    return create_ccomplex_return(IT_0601);
}

