#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_st_1_to_b_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_st_1_to_b_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
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
    const ccomplex_t IT_0000 = m_b*m_N_2;
    const ccomplex_t IT_0001 = pow(m_b, 2);
    const ccomplex_t IT_0002 = pow(m_N_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cpow(IT_0005, -1);
    const ccomplex_t IT_0008 = cpow(IT_0004, -1);
    const ccomplex_t IT_0009 = m_t*conj(U_sb_00)*(m_t*IT_0004*IT_0007*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0005*IT_0008*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0007*IT_0008*U_st_10));
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0011*((-0.5)*IT_0006
      *IT_0009 + m_W*IT_0004*IT_0005*conj(U_sb_00)*U_st_00) + A_t*IT_0004*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0013 = IT_0008*IT_0011;
    const ccomplex_t IT_0014 = m_b*conj(N_d2)*e_em*IT_0006*U_sb_10;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = conj(N_B2)*e_em*U_sb_00;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W2)*e_em*U_sb_00;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0016 + 0.333333333333333
      *IT_0021 + -IT_0024);
    const ccomplex_t IT_0026 = 3*IT_0025;
    const ccomplex_t IT_0027 = 0.166666666666667*IT_0026;
    const ccomplex_t IT_0028 = IT_0012*IT_0027;
    const ccomplex_t IT_0029 = IT_0003*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0001 + IT_0002 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0032 = m_t*conj(U_sb_01)*(m_t*IT_0004*IT_0007*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0005*IT_0008*conj(U_sb_01)*U_st_00 + conj
      (U_sb_11)*(mu_h*U_st_00 + m_t*IT_0007*IT_0008*U_st_10));
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0005*conj
      (U_sb_11)*U_st_00 + (-1.4142135623731)*e_em*IT_0011*((-0.5)*IT_0006
      *IT_0032 + m_W*IT_0004*IT_0005*conj(U_sb_01)*U_st_00) + A_t*IT_0004*conj
      (U_sb_01)*U_st_10);
    const ccomplex_t IT_0034 = m_b*conj(N_d2)*e_em*IT_0006*U_sb_11;
    const ccomplex_t IT_0035 = IT_0013*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0038 = IT_0011*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0041 = IT_0018*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0036 + -IT_0039 +
       0.333333333333333*IT_0042);
    const ccomplex_t IT_0044 = 3*IT_0043;
    const ccomplex_t IT_0045 = 0.166666666666667*IT_0044;
    const ccomplex_t IT_0046 = IT_0033*IT_0045;
    const ccomplex_t IT_0047 = IT_0031*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0005*IT_0011;
    const ccomplex_t IT_0050 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (-2)*IT_0051;
    const ccomplex_t IT_0053 = IT_0005*IT_0018;
    const ccomplex_t IT_0054 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0058 = IT_0049*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = IT_0056 + IT_0059;
    const ccomplex_t IT_0061 = IT_0052 + IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = V_tb*e_em*conj(V_Wp1)*U_st_00;
    const ccomplex_t IT_0065 = IT_0011*IT_0064;
    const ccomplex_t IT_0066 = IT_0007*IT_0011;
    const ccomplex_t IT_0067 = m_t*V_tb*conj(V_u1)*e_em*IT_0006*U_st_10;
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0065 + (-0.5)*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0063*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0074 = pow(m_st_1, 2);
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0001 + IT_0074 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0072*IT_0073*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0004*IT_0006*IT_0007*IT_0011;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B2)*e_em*U_st_00;
    const ccomplex_t IT_0081 = IT_0018*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W2)*e_em*U_st_00;
    const ccomplex_t IT_0084 = IT_0011*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = m_t*conj(N_u2)*e_em*IT_0006*U_st_10;
    const ccomplex_t IT_0087 = IT_0066*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0082 + 3*IT_0085 + 3
      *IT_0088);
    const ccomplex_t IT_0090 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0091 = IT_0079*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_12 + -IT_0002 + -IT_0074 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0093 = IT_0077*IT_0091*IT_0092;
    const ccomplex_t IT_0094 = N_B2*e_em*U_st_10;
    const ccomplex_t IT_0095 = IT_0018*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = m_t*N_u2*e_em*IT_0006*U_st_00;
    const ccomplex_t IT_0098 = IT_0066*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0096 + (-0.75)*IT_0099);
    const ccomplex_t IT_0101 = 0.666666666666667*IT_0100;
    const ccomplex_t IT_0102 = IT_0079*IT_0101;
    const ccomplex_t IT_0103 = IT_0092*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = m_N_2*IT_0104;
    const ccomplex_t IT_0106 = -IT_0030 + -IT_0048 + IT_0076 + IT_0093 + 
      -IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_b*U_d1*V_tb
      *e_em*IT_0006*IT_0008*IT_0011*U_st_00;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0063*IT_0108;
    const ccomplex_t IT_0110 = IT_0075*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = m_b*IT_0111;
    const ccomplex_t IT_0113 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0114 = IT_0049*IT_0113;
    const ccomplex_t IT_0115 = (-2)*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0117 = IT_0053*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0120 = IT_0049*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = IT_0118 + IT_0121;
    const ccomplex_t IT_0123 = IT_0115 + IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = V_tb*e_em*conj(V_Wp2)*U_st_00;
    const ccomplex_t IT_0127 = IT_0011*IT_0126;
    const ccomplex_t IT_0128 = m_t*V_tb*conj(V_u2)*e_em*IT_0006*U_st_10;
    const ccomplex_t IT_0129 = IT_0066*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + (-0.5)*IT_0130);
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0125*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0001 + IT_0074 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0136 = IT_0133*IT_0134*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_b*U_d2*V_tb
      *e_em*IT_0006*IT_0008*IT_0011*U_st_00;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = IT_0125*IT_0138;
    const ccomplex_t IT_0140 = IT_0135*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = m_b*IT_0141;
    const ccomplex_t IT_0143 = -IT_0112 + -IT_0136 + -IT_0142;
    const ccomplex_t IT_0144 = N_B2*e_em*U_sb_10;
    const ccomplex_t IT_0145 = IT_0018*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = m_b*N_d2*e_em*IT_0006*U_sb_00;
    const ccomplex_t IT_0148 = IT_0013*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + 1.5*IT_0149);
    const ccomplex_t IT_0151 = (-0.333333333333333)*IT_0150;
    const ccomplex_t IT_0152 = IT_0012*IT_0151;
    const ccomplex_t IT_0153 = IT_0003*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0004*IT_0011;
    const ccomplex_t IT_0156 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = 2*IT_0157;
    const ccomplex_t IT_0159 = IT_0004*IT_0018;
    const ccomplex_t IT_0160 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0161 = IT_0159*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0164 = IT_0155*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = IT_0162 + IT_0165;
    const ccomplex_t IT_0167 = IT_0158 + IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = IT_0138*IT_0169;
    const ccomplex_t IT_0171 = IT_0134*IT_0135*IT_0170;
    const ccomplex_t IT_0172 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0173 = IT_0018*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = m_b*N_d2*e_em*IT_0006*U_sb_01;
    const ccomplex_t IT_0176 = IT_0013*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 1.5*IT_0177);
    const ccomplex_t IT_0179 = (-0.333333333333333)*IT_0178;
    const ccomplex_t IT_0180 = IT_0033*IT_0179;
    const ccomplex_t IT_0181 = IT_0031*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0005*IT_0006*IT_0008*IT_0011;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0101*IT_0184;
    const ccomplex_t IT_0186 = IT_0077*IT_0092*IT_0185;
    const ccomplex_t IT_0187 = IT_0090*IT_0184;
    const ccomplex_t IT_0188 = IT_0092*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = m_N_2*IT_0189;
    const ccomplex_t IT_0191 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0192 = IT_0159*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0195 = IT_0155*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = IT_0193 + IT_0196;
    const ccomplex_t IT_0198 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0199 = IT_0155*IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = IT_0197 + IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0202;
    const ccomplex_t IT_0204 = IT_0108*IT_0203;
    const ccomplex_t IT_0205 = IT_0073*IT_0075*IT_0204;
    const ccomplex_t IT_0206 = IT_0154 + -IT_0171 + IT_0182 + -IT_0186 +
       IT_0190 + -IT_0205;
    const ccomplex_t IT_0207 = IT_0132*IT_0169;
    const ccomplex_t IT_0208 = IT_0135*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = IT_0071*IT_0203;
    const ccomplex_t IT_0211 = IT_0075*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = IT_0209 + IT_0212;
    const ccomplex_t IT_0214 = m_b*IT_0213;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0206 + IT_0215;
    const ccomplex_t IT_0217 = IT_0106 + -IT_0143;
    const ccomplex_t IT_0218 = m_b*m_N_2*IT_0074;
    const ccomplex_t IT_0219 = -IT_0111 + -IT_0141;
    const ccomplex_t IT_0220 = -IT_0104;
    const ccomplex_t IT_0221 = IT_0213 + IT_0220;
    const ccomplex_t IT_0222 = IT_0189 + IT_0219;
    const ccomplex_t IT_0223 = conj(IT_0213) + conj(IT_0220);
    const ccomplex_t IT_0224 = m_b*s_12;
    const ccomplex_t IT_0225 = conj(IT_0189) + conj(IT_0219);
    const ccomplex_t IT_0226 = (-6)*IT_0143;
    const ccomplex_t IT_0227 = IT_0106 + 0.166666666666667*IT_0226;
    const ccomplex_t IT_0228 = (-6)*conj(IT_0143);
    const ccomplex_t IT_0229 = 0.166666666666667*IT_0228;
    const ccomplex_t IT_0230 = conj(IT_0106) + IT_0229;
    const ccomplex_t IT_0231 = -conj(IT_0214);
    const ccomplex_t IT_0232 = conj(IT_0206) + IT_0231;
    const ccomplex_t IT_0233 = s_23*m_N_2;
    const ccomplex_t IT_0234 = s_13*IT_0074;
    const ccomplex_t IT_0235 = (-6)*IT_0234;
    const ccomplex_t IT_0236 = s_12*s_23;
    const ccomplex_t IT_0237 = 12*IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = 6*s_13;
    const ccomplex_t IT_0240 = (-6)*s_13;
    const ccomplex_t IT_0241 = 6*IT_0000*((conj(IT_0106) + -conj(IT_0143))
      *IT_0216 + (conj(IT_0206) + -conj(IT_0214))*IT_0217) + 6*IT_0218*((conj
      (IT_0189) + conj(IT_0219))*IT_0221 + IT_0222*IT_0223) + 6*IT_0224*(IT_0216
      *IT_0223 + IT_0225*IT_0227 + IT_0222*IT_0230 + IT_0221*IT_0232) + 6*
      (IT_0216*IT_0225 + IT_0223*IT_0227 + IT_0221*IT_0230 + IT_0222*IT_0232)
      *IT_0233 + (IT_0221*IT_0223 + IT_0222*IT_0225)*IT_0238 + (IT_0106*conj
      (IT_0106) + IT_0143*conj(IT_0143) + IT_0206*conj(IT_0206) + IT_0214*conj
      (IT_0214))*IT_0239 + (conj(IT_0106)*IT_0143 + IT_0106*conj(IT_0143) + conj
      (IT_0206)*IT_0214 + IT_0206*conj(IT_0214))*IT_0240;
    return create_ccomplex_return(IT_0241);
}

