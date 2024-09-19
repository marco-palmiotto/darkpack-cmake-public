#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_st_1_to_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_st_1_to_b_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_st_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 0.5*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0004*IT_0012*IT_0014*U_st_00;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0021 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + 0.5*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0004*IT_0012*IT_0014*U_st_00;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = cos(theta_W);
    const ccomplex_t IT_0036 = cpow(IT_0035, -1);
    const ccomplex_t IT_0037 = N_B1*e_em*U_sb_10;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = IT_0004*IT_0014;
    const ccomplex_t IT_0041 = m_b*N_d1*e_em*IT_0012*U_sb_00;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0039 + 1.5*IT_0043);
    const ccomplex_t IT_0045 = (-0.333333333333333)*IT_0044;
    const ccomplex_t IT_0046 = pow(m_N_1, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0000 + IT_0046 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0034*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0054 = IT_0036*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = m_b*N_d1*e_em*IT_0012*U_sb_01;
    const ccomplex_t IT_0057 = IT_0040*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + 1.5*IT_0058);
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0000 + IT_0046 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = IT_0052*IT_0063;
    const ccomplex_t IT_0065 = (-2)*IT_0019 + (-2)*IT_0032 + IT_0050 + IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -2);
    const ccomplex_t IT_0067 = pow(s_14, 2);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = -IT_0046;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = s_13*IT_0070;
    const ccomplex_t IT_0072 = m_b*conj(N_d1)*e_em*IT_0012*U_sb_11;
    const ccomplex_t IT_0073 = IT_0040*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0076 = IT_0004*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0079 = IT_0036*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0074 + -IT_0077 +
       0.333333333333333*IT_0080);
    const ccomplex_t IT_0082 = 3*IT_0081;
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0082;
    const ccomplex_t IT_0084 = IT_0061*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0052*IT_0085;
    const ccomplex_t IT_0087 = m_b*conj(N_d1)*e_em*IT_0012*U_sb_10;
    const ccomplex_t IT_0088 = IT_0040*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W1)*e_em*U_sb_00;
    const ccomplex_t IT_0091 = IT_0004*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B1)*e_em*U_sb_00;
    const ccomplex_t IT_0094 = IT_0036*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0089 + -IT_0092 +
       0.333333333333333*IT_0095);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = (-3)*IT_0097;
    const ccomplex_t IT_0099 = (-0.333333333333333)*IT_0098;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = 3*IT_0100;
    const ccomplex_t IT_0102 = 0.166666666666667*IT_0101;
    const ccomplex_t IT_0103 = IT_0047*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0034*IT_0104;
    const ccomplex_t IT_0106 = -IT_0086 + -IT_0105;
    const ccomplex_t IT_0107 = m_b*m_N_1;
    const ccomplex_t IT_0108 = IT_0070*IT_0107;
    const ccomplex_t IT_0109 = IT_0019 + IT_0032;
    const ccomplex_t IT_0110 = s_13*s_14;
    const ccomplex_t IT_0111 = s_34*IT_0046;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0110 + IT_0112;
    const ccomplex_t IT_0114 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0115 = IT_0004*IT_0114;
    const ccomplex_t IT_0116 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0117 = IT_0004*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + (-0.5)*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0122 = IT_0004*IT_0121;
    const ccomplex_t IT_0123 = sin(beta);
    const ccomplex_t IT_0124 = cpow(IT_0123, -1);
    const ccomplex_t IT_0125 = IT_0004*IT_0124;
    const ccomplex_t IT_0126 = m_t*V_tb*conj(V_u2)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0122 + (-0.5)*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0120*IT_0130;
    const ccomplex_t IT_0132 = IT_0002*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0135 = IT_0004*IT_0134;
    const ccomplex_t IT_0136 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0137 = IT_0004*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0135 + (-0.5)*IT_0138);
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0142 = IT_0004*IT_0141;
    const ccomplex_t IT_0143 = m_t*V_tb*conj(V_u1)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0144 = IT_0125*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*(IT_0142 + (-0.5)*IT_0145);
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = IT_0140*IT_0147;
    const ccomplex_t IT_0149 = IT_0020*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = -IT_0133 + -IT_0150;
    const ccomplex_t IT_0152 = pow(m_W, 2);
    const ccomplex_t IT_0153 = IT_0107*IT_0152;
    const ccomplex_t IT_0154 = (-3)*IT_0153;
    const ccomplex_t IT_0155 = 3*IT_0153;
    const ccomplex_t IT_0156 = cpow((-2)*s_12 + -IT_0001 + -IT_0046 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0157 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0158 = IT_0036*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0161 = IT_0004*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = m_t*conj(N_u1)*e_em*IT_0012*U_st_10;
    const ccomplex_t IT_0164 = IT_0125*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0159 + 3*IT_0162 + 3
      *IT_0165);
    const ccomplex_t IT_0167 = 0.166666666666667*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0004;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = IT_0167*IT_0169;
    const ccomplex_t IT_0171 = IT_0156*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = -IT_0050 + -IT_0064;
    const ccomplex_t IT_0174 = s_14*IT_0000;
    const ccomplex_t IT_0175 = s_13*s_34;
    const ccomplex_t IT_0176 = -IT_0175;
    const ccomplex_t IT_0177 = IT_0174 + IT_0176;
    const ccomplex_t IT_0178 = pow(m_W, 4);
    const ccomplex_t IT_0179 = s_13*IT_0178;
    const ccomplex_t IT_0180 = IT_0066*IT_0179;
    const ccomplex_t IT_0181 = s_14*s_34;
    const ccomplex_t IT_0182 = 24*IT_0181;
    const ccomplex_t IT_0183 = s_13*IT_0152;
    const ccomplex_t IT_0184 = (-12)*IT_0183;
    const ccomplex_t IT_0185 = 6*IT_0180 + IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = s_14*IT_0152;
    const ccomplex_t IT_0187 = IT_0066*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = s_14 + IT_0188;
    const ccomplex_t IT_0190 = m_b*IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0193 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0194 = IT_0036*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = m_t*N_u1*e_em*IT_0012*U_st_00;
    const ccomplex_t IT_0197 = IT_0125*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + (-0.75)*IT_0198);
    const ccomplex_t IT_0200 = 0.666666666666667*IT_0199;
    const ccomplex_t IT_0201 = IT_0169*IT_0200;
    const ccomplex_t IT_0202 = IT_0156*IT_0192*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0204 = IT_0011*IT_0130;
    const ccomplex_t IT_0205 = IT_0002*IT_0203*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0207 = IT_0027*IT_0147;
    const ccomplex_t IT_0208 = IT_0020*IT_0206*IT_0207;
    const ccomplex_t IT_0209 = -IT_0202 + -IT_0205 + -IT_0208;
    const ccomplex_t IT_0210 = IT_0133 + IT_0150;
    const ccomplex_t IT_0211 = m_N_1*IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = s_34*IT_0152;
    const ccomplex_t IT_0214 = IT_0066*IT_0213;
    const ccomplex_t IT_0215 = (-0.25)*IT_0214;
    const ccomplex_t IT_0216 = s_34 + IT_0215;
    const ccomplex_t IT_0217 = m_N_1*IT_0216;
    const ccomplex_t IT_0218 = (-4)*IT_0217;
    const ccomplex_t IT_0219 = m_N_1*IT_0032;
    const ccomplex_t IT_0220 = m_N_1*IT_0019;
    const ccomplex_t IT_0221 = m_b*IT_0172;
    const ccomplex_t IT_0222 = IT_0219 + IT_0220 + -IT_0221;
    const ccomplex_t IT_0223 = IT_0029*IT_0140;
    const ccomplex_t IT_0224 = IT_0020*IT_0206*IT_0223;
    const ccomplex_t IT_0225 = IT_0016*IT_0120;
    const ccomplex_t IT_0226 = IT_0002*IT_0203*IT_0225;
    const ccomplex_t IT_0227 = -IT_0224 + -IT_0226;
    const ccomplex_t IT_0228 = IT_0222 + IT_0227;
    const ccomplex_t IT_0229 = 6*conj(IT_0109)*(IT_0065*IT_0113 + -IT_0151
      *IT_0154 + IT_0155*IT_0172 + IT_0173*IT_0177 + 0.166666666666667*IT_0109
      *IT_0185 + IT_0191*IT_0212 + IT_0218*IT_0228);
    const ccomplex_t IT_0230 = (-18)*conj(IT_0209);
    const ccomplex_t IT_0231 = s_14*s_34*IT_0066;
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = s_13 + IT_0232;
    const ccomplex_t IT_0234 = s_13*IT_0233;
    const ccomplex_t IT_0235 = 2*IT_0133 + 2*IT_0150;
    const ccomplex_t IT_0236 = IT_0086 + IT_0105 + 2*IT_0172;
    const ccomplex_t IT_0237 = IT_0107*IT_0233;
    const ccomplex_t IT_0238 = m_N_1*IT_0233;
    const ccomplex_t IT_0239 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0240 = 0.5*IT_0214;
    const ccomplex_t IT_0241 = s_34 + IT_0240;
    const ccomplex_t IT_0242 = m_N_1*IT_0241;
    const ccomplex_t IT_0243 = 2*IT_0242;
    const ccomplex_t IT_0244 = (-2)*IT_0181;
    const ccomplex_t IT_0245 = IT_0183 + IT_0244;
    const ccomplex_t IT_0246 = IT_0066*IT_0245;
    const ccomplex_t IT_0247 = (-6)*IT_0246;
    const ccomplex_t IT_0248 = 12*s_13;
    const ccomplex_t IT_0249 = IT_0247 + IT_0248;
    const ccomplex_t IT_0250 = m_b*IT_0233;
    const ccomplex_t IT_0251 = (-6)*IT_0173;
    const ccomplex_t IT_0252 = pow(s_34, 2);
    const ccomplex_t IT_0253 = IT_0066*IT_0252;
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = IT_0000 + IT_0254;
    const ccomplex_t IT_0256 = IT_0107*IT_0255;
    const ccomplex_t IT_0257 = m_N_1*IT_0255;
    const ccomplex_t IT_0258 = (-6)*conj(IT_0173);
    const ccomplex_t IT_0259 = (-6)*IT_0180;
    const ccomplex_t IT_0260 = IT_0182 + IT_0259;
    const ccomplex_t IT_0261 = (-6)*IT_0106;
    const ccomplex_t IT_0262 = 0.166666666666667*IT_0261;
    const ccomplex_t IT_0263 = IT_0235 + IT_0262;
    const ccomplex_t IT_0264 = (-0.25)*IT_0187;
    const ccomplex_t IT_0265 = s_14 + IT_0264;
    const ccomplex_t IT_0266 = m_b*IT_0265;
    const ccomplex_t IT_0267 = (-4)*IT_0266;
    const ccomplex_t IT_0268 = (-6)*IT_0236;
    const ccomplex_t IT_0269 = IT_0151 + IT_0172;
    const ccomplex_t IT_0270 = (-6)*conj(IT_0236);
    const ccomplex_t IT_0271 = (-6)*IT_0209;
    const ccomplex_t IT_0272 = (-6)*conj(IT_0151)*(IT_0109*IT_0154 + (
      -0.166666666666667)*IT_0151*IT_0185 + IT_0191*IT_0228 + (
      -0.166666666666667)*IT_0172*IT_0260 + -IT_0113*IT_0263 + (
      -0.166666666666667)*IT_0177*IT_0268 + IT_0218*(IT_0211 + (
      -0.166666666666667)*IT_0271));
    const ccomplex_t IT_0273 = (-6)*conj(IT_0209);
    const ccomplex_t IT_0274 = (-6)*conj(IT_0235);
    const ccomplex_t IT_0275 = m_b*IT_0070;
    const ccomplex_t IT_0276 = 6*IT_0275;
    const ccomplex_t IT_0277 = 6*IT_0109;
    const ccomplex_t IT_0278 = 6*IT_0172;
    const ccomplex_t IT_0279 = (-6)*IT_0151*IT_0218 + IT_0107*((-18)*IT_0222 +
       (-18)*IT_0227) + IT_0212*IT_0249 + IT_0250*IT_0251 + 6*IT_0238*IT_0263 +
       IT_0257*IT_0268 + 6*IT_0065*IT_0275 + IT_0191*IT_0277 + IT_0243*IT_0278;
    const ccomplex_t IT_0280 = (-6)*IT_0065;
    const ccomplex_t IT_0281 = (-6)*IT_0235;
    const ccomplex_t IT_0282 = (-6)*IT_0151*IT_0191 + IT_0107*((-18)*IT_0209 +
       (-18)*IT_0211) + IT_0228*IT_0249 + IT_0251*IT_0257 + IT_0250*IT_0268 +
       IT_0218*IT_0277 + IT_0267*IT_0278 + IT_0238*IT_0280 + IT_0276*(IT_0106 +
       0.166666666666667*IT_0281);
    const ccomplex_t IT_0283 = 0.166666666666667*IT_0237;
    const ccomplex_t IT_0284 = s_13*IT_0255;
    const ccomplex_t IT_0285 = 0.166666666666667*IT_0284;
    const ccomplex_t IT_0286 = (-6)*conj(IT_0106);
    const ccomplex_t IT_0287 = IT_0113*IT_0269;
    const ccomplex_t IT_0288 = 6*conj(IT_0065)*(IT_0065*IT_0071 + IT_0106
      *IT_0108 + IT_0109*IT_0113) + IT_0229 + IT_0107*IT_0228*IT_0230 + (IT_0173
      *IT_0234 + IT_0108*IT_0235 + IT_0236*IT_0237 + IT_0228*IT_0238)*IT_0239 +
       6*conj(IT_0209)*(IT_0109*IT_0191 + IT_0235*IT_0238 + IT_0172*IT_0243 +
       0.166666666666667*IT_0212*IT_0249 + 0.166666666666667*IT_0250*IT_0251) + 
      (IT_0106*IT_0237 + IT_0212*IT_0250 + IT_0236*IT_0256 + IT_0228*IT_0257)
      *IT_0258 + 6*conj(IT_0172)*(IT_0109*IT_0155 + 0.166666666666667*IT_0172
      *IT_0185 + IT_0212*IT_0243 + 0.166666666666667*IT_0151*IT_0260 + IT_0113
      *IT_0263 + IT_0228*IT_0267 + 0.166666666666667*IT_0177*IT_0268) + (IT_0228
      *IT_0250 + IT_0211*IT_0257 + IT_0177*IT_0269)*IT_0270 + IT_0272 + (IT_0151
      *IT_0218 + IT_0106*IT_0238 + IT_0236*IT_0257)*IT_0273 + (IT_0065*conj
      (IT_0209) + conj(IT_0065)*IT_0212 + IT_0228*(conj(IT_0106) +
       0.166666666666667*IT_0274))*IT_0276 + conj(IT_0211)*IT_0279 + (conj
      (IT_0222) + conj(IT_0227))*IT_0282 + 6*conj(IT_0106)*(IT_0071*IT_0106 +
       IT_0065*IT_0108 + 0.166666666666667*IT_0238*IT_0271 + IT_0251*IT_0283) +
       6*conj(IT_0173)*(IT_0109*IT_0177 + IT_0235*IT_0237 + 0.166666666666667
      *IT_0234*IT_0280 + IT_0251*IT_0285) + 6*conj(IT_0236)*(0.166666666666667
      *IT_0251*IT_0256 + IT_0234*IT_0263 + 0.166666666666667*IT_0257*IT_0271 +
       IT_0280*IT_0283 + IT_0268*IT_0285) + IT_0286*(IT_0071*IT_0235 + IT_0234
      *IT_0236 + IT_0211*IT_0238 + IT_0287) + 6*conj(IT_0235)*(IT_0234*IT_0236 +
       IT_0173*IT_0237 + IT_0212*IT_0238 + IT_0071*IT_0263 + 0.166666666666667
      *IT_0108*IT_0280 + IT_0287);
    return create_ccomplex_return(IT_0288);
}

