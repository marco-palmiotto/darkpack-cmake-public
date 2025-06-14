#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_sb_1_to_anti_b_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_sb_1_to_anti_b_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_sb_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = m_b*conj(N_d1)*e_em*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0007 + 1.5*IT_0016);
    const ccomplex_t IT_0018 = (-0.333333333333333)*IT_0017;
    const ccomplex_t IT_0019 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0022 = IT_0011*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0020 + 0.5*IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0018*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0032 = IT_0011*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + 0.5*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0037 = IT_0004*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = m_b*conj(N_d3)*e_em*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + 1.5*IT_0041);
    const ccomplex_t IT_0043 = (-0.333333333333333)*IT_0042;
    const ccomplex_t IT_0044 = IT_0035*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0049 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0050 = IT_0011*IT_0049;
    const ccomplex_t IT_0051 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0052 = IT_0011*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + 0.5*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = m_b*conj(N_d4)*e_em*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0060 = IT_0012*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 1.5*IT_0061);
    const ccomplex_t IT_0063 = (-0.333333333333333)*IT_0062;
    const ccomplex_t IT_0064 = IT_0055*IT_0063;
    const ccomplex_t IT_0065 = IT_0048*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0068 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0069 = IT_0011*IT_0068;
    const ccomplex_t IT_0070 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0071 = IT_0011*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 0.5*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0076 = IT_0004*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = m_b*conj(N_d2)*e_em*IT_0013*conj(U_sb_00);
    const ccomplex_t IT_0079 = IT_0012*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 1.5*IT_0080);
    const ccomplex_t IT_0082 = (-0.333333333333333)*IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = IT_0067*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0028 + IT_0047 + IT_0066 + IT_0085;
    const ccomplex_t IT_0087 = m_C_2*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0089 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0090 = IT_0011*IT_0089;
    const ccomplex_t IT_0091 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0092 = IT_0011*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + (-0.5)*IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0018*IT_0095;
    const ccomplex_t IT_0097 = IT_0002*IT_0088*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0099 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0100 = IT_0011*IT_0099;
    const ccomplex_t IT_0101 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0102 = IT_0011*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + (-0.5)*IT_0103);
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0082*IT_0105;
    const ccomplex_t IT_0107 = IT_0067*IT_0098*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0109 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0110 = IT_0011*IT_0109;
    const ccomplex_t IT_0111 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0112 = IT_0011*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + (-0.5)*IT_0113);
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0043*IT_0115;
    const ccomplex_t IT_0117 = IT_0045*IT_0108*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0119 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0120 = IT_0011*IT_0119;
    const ccomplex_t IT_0121 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0122 = IT_0011*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + (-0.5)*IT_0123);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = IT_0063*IT_0125;
    const ccomplex_t IT_0127 = IT_0048*IT_0118*IT_0126;
    const ccomplex_t IT_0128 = -IT_0097 + -IT_0107 + -IT_0117 + -IT_0127;
    const ccomplex_t IT_0129 = pow(m_W, -2);
    const ccomplex_t IT_0130 = pow(m_W, 2);
    const ccomplex_t IT_0131 = s_13*IT_0130;
    const ccomplex_t IT_0132 = s_14*s_34;
    const ccomplex_t IT_0133 = (-2)*IT_0132;
    const ccomplex_t IT_0134 = IT_0131 + IT_0133;
    const ccomplex_t IT_0135 = IT_0129*IT_0134;
    const ccomplex_t IT_0136 = (-6)*IT_0135;
    const ccomplex_t IT_0137 = 12*s_13;
    const ccomplex_t IT_0138 = IT_0136 + IT_0137;
    const ccomplex_t IT_0139 = pow(s_14, 2);
    const ccomplex_t IT_0140 = IT_0129*IT_0139;
    const ccomplex_t IT_0141 = pow(m_C_2, 2);
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = m_b*IT_0143;
    const ccomplex_t IT_0145 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0146 = IT_0011*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = m_b*N_d2*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0149 = IT_0012*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0152 = IT_0004*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0147 + -IT_0150 + (
      -0.333333333333333)*IT_0153);
    const ccomplex_t IT_0155 = (-3)*IT_0154;
    const ccomplex_t IT_0156 = 0.166666666666667*IT_0155;
    const ccomplex_t IT_0157 = IT_0105*IT_0156;
    const ccomplex_t IT_0158 = IT_0067*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0161 = IT_0004*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0164 = IT_0011*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = m_b*N_d1*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0167 = IT_0012*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0162 + (-3)*IT_0165 + 3
      *IT_0168);
    const ccomplex_t IT_0170 = 0.166666666666667*IT_0169;
    const ccomplex_t IT_0171 = IT_0095*IT_0170;
    const ccomplex_t IT_0172 = IT_0002*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0175 = IT_0004*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0178 = IT_0011*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = m_b*N_d3*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0181 = IT_0012*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0176 + (-3)*IT_0179 + 3
      *IT_0182);
    const ccomplex_t IT_0184 = 0.166666666666667*IT_0183;
    const ccomplex_t IT_0185 = IT_0115*IT_0184;
    const ccomplex_t IT_0186 = IT_0045*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0189 = IT_0004*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0192 = IT_0011*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_b*N_d4*e_em*IT_0013*conj(U_sb_10);
    const ccomplex_t IT_0195 = IT_0012*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0190 + (-3)*IT_0193 + 3
      *IT_0196);
    const ccomplex_t IT_0198 = 0.166666666666667*IT_0197;
    const ccomplex_t IT_0199 = IT_0125*IT_0198;
    const ccomplex_t IT_0200 = IT_0048*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = 2*IT_0159 + 2*IT_0173 + 2*IT_0187 + 2*IT_0201;
    const ccomplex_t IT_0203 = s_14*s_34*IT_0129;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = s_13 + IT_0204;
    const ccomplex_t IT_0206 = m_b*IT_0205;
    const ccomplex_t IT_0207 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0208 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0011
      *IT_0207*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = cpow((-2)*s_13 + IT_0000 + IT_0141 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0211 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0011*IT_0211
      *V_ub_mod;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = IT_0210*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0209*IT_0215;
    const ccomplex_t IT_0217 = cpow((-2)*s_13 + IT_0000 + IT_0141 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0220 = IT_0011*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0011*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0222*IT_0224;
    const ccomplex_t IT_0226 = sin(beta);
    const ccomplex_t IT_0227 = cpow(IT_0226, -1);
    const ccomplex_t IT_0228 = IT_0011*IT_0227;
    const ccomplex_t IT_0229 = m_t*V_tb*V_u2*e_em*IT_0013*conj(U_st_10);
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = IT_0224*IT_0233;
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = IT_0225 + IT_0235;
    const ccomplex_t IT_0237 = IT_0218*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0011*conj(U_sb_00);
    const ccomplex_t IT_0239 = -IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*V_cb*e_em*V_Wp2*IT_0011;
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = cpow((-2)*s_13 + IT_0000 + IT_0141 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0243 = IT_0241*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = IT_0239*IT_0244;
    const ccomplex_t IT_0246 = cpow((-2)*s_13 + IT_0000 + IT_0141 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0249 = IT_0011*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0011*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = IT_0251*IT_0253;
    const ccomplex_t IT_0255 = m_t*V_tb*V_u2*e_em*IT_0013*conj(U_st_11);
    const ccomplex_t IT_0256 = IT_0228*IT_0255;
    const ccomplex_t IT_0257 = 1.4142135623731*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = IT_0253*IT_0259;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0254 + IT_0261;
    const ccomplex_t IT_0263 = IT_0247*IT_0262;
    const ccomplex_t IT_0264 = IT_0216 + IT_0237 + IT_0245 + IT_0263;
    const ccomplex_t IT_0265 = m_C_2*IT_0205;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0009*IT_0011*IT_0013*conj(U_st_01);
    const ccomplex_t IT_0267 = (-0.5)*IT_0266;
    const ccomplex_t IT_0268 = IT_0246*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = IT_0253*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0009*IT_0011*IT_0013*conj(U_st_00);
    const ccomplex_t IT_0272 = (-0.5)*IT_0271;
    const ccomplex_t IT_0273 = IT_0217*IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = IT_0224*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_cb*e_em*IT_0009*IT_0011*IT_0013;
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = IT_0242*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0239*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *e_em*IT_0009*IT_0011*IT_0013*IT_0211*V_ub_mod;
    const ccomplex_t IT_0282 = (-0.5)*IT_0281;
    const ccomplex_t IT_0283 = IT_0210*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0209*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0028 + 2*IT_0047 + 2*IT_0066 + 2*IT_0085 +
       IT_0270 + IT_0275 + IT_0280 + IT_0285;
    const ccomplex_t IT_0287 = pow(s_34, 2);
    const ccomplex_t IT_0288 = IT_0129*IT_0287;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = IT_0000 + IT_0289;
    const ccomplex_t IT_0291 = m_C_2*IT_0290;
    const ccomplex_t IT_0292 = -IT_0270 + -IT_0275 + -IT_0280 + -IT_0285;
    const ccomplex_t IT_0293 = s_14*IT_0130;
    const ccomplex_t IT_0294 = IT_0129*IT_0293;
    const ccomplex_t IT_0295 = 0.5*IT_0294;
    const ccomplex_t IT_0296 = s_14 + IT_0295;
    const ccomplex_t IT_0297 = m_b*IT_0296;
    const ccomplex_t IT_0298 = 2*IT_0297;
    const ccomplex_t IT_0299 = -IT_0159 + -IT_0173 + -IT_0187 + -IT_0201;
    const ccomplex_t IT_0300 = 6*IT_0087;
    const ccomplex_t IT_0301 = s_34*IT_0130;
    const ccomplex_t IT_0302 = IT_0129*IT_0301;
    const ccomplex_t IT_0303 = (-0.25)*IT_0302;
    const ccomplex_t IT_0304 = s_34 + IT_0303;
    const ccomplex_t IT_0305 = m_C_2*IT_0304;
    const ccomplex_t IT_0306 = (-4)*IT_0305;
    const ccomplex_t IT_0307 = -IT_0028 + -IT_0047 + -IT_0066 + -IT_0085;
    const ccomplex_t IT_0308 = -IT_0216 + -IT_0237 + -IT_0245 + -IT_0263;
    const ccomplex_t IT_0309 = (-6)*IT_0087;
    const ccomplex_t IT_0310 = s_13*s_14;
    const ccomplex_t IT_0311 = s_34*IT_0141;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = IT_0310 + IT_0312;
    const ccomplex_t IT_0314 = 6*IT_0313;
    const ccomplex_t IT_0315 = 6*conj(IT_0087);
    const ccomplex_t IT_0316 = 6*IT_0144;
    const ccomplex_t IT_0317 = s_13*IT_0143;
    const ccomplex_t IT_0318 = 6*IT_0317;
    const ccomplex_t IT_0319 = s_13*IT_0205;
    const ccomplex_t IT_0320 = 6*IT_0319;
    const ccomplex_t IT_0321 = m_b*m_C_2;
    const ccomplex_t IT_0322 = IT_0205*IT_0321;
    const ccomplex_t IT_0323 = 6*IT_0322;
    const ccomplex_t IT_0324 = (-6)*IT_0317;
    const ccomplex_t IT_0325 = IT_0159 + IT_0173 + IT_0187 + IT_0201;
    const ccomplex_t IT_0326 = m_C_2*IT_0325;
    const ccomplex_t IT_0327 = -IT_0326;
    const ccomplex_t IT_0328 = IT_0035*IT_0184;
    const ccomplex_t IT_0329 = IT_0045*IT_0108*IT_0328;
    const ccomplex_t IT_0330 = IT_0074*IT_0156;
    const ccomplex_t IT_0331 = IT_0067*IT_0098*IT_0330;
    const ccomplex_t IT_0332 = IT_0025*IT_0170;
    const ccomplex_t IT_0333 = IT_0002*IT_0088*IT_0332;
    const ccomplex_t IT_0334 = IT_0055*IT_0198;
    const ccomplex_t IT_0335 = IT_0048*IT_0118*IT_0334;
    const ccomplex_t IT_0336 = IT_0329 + IT_0331 + IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = conj(IT_0327) + conj(IT_0336);
    const ccomplex_t IT_0338 = (-6)*IT_0265;
    const ccomplex_t IT_0339 = IT_0143*IT_0321;
    const ccomplex_t IT_0340 = (-6)*IT_0339;
    const ccomplex_t IT_0341 = 6*IT_0265;
    const ccomplex_t IT_0342 = 6*IT_0339;
    const ccomplex_t IT_0343 = (-6)*IT_0319;
    const ccomplex_t IT_0344 = (-6)*IT_0322;
    const ccomplex_t IT_0345 = 6*IT_0206;
    const ccomplex_t IT_0346 = 6*IT_0291;
    const ccomplex_t IT_0347 = 6*IT_0298;
    const ccomplex_t IT_0348 = (-18)*IT_0321;
    const ccomplex_t IT_0349 = (-6)*IT_0144;
    const ccomplex_t IT_0350 = (-6)*IT_0306;
    const ccomplex_t IT_0351 = pow(m_W, 4);
    const ccomplex_t IT_0352 = s_13*IT_0351;
    const ccomplex_t IT_0353 = IT_0129*IT_0352;
    const ccomplex_t IT_0354 = 6*IT_0353;
    const ccomplex_t IT_0355 = 24*IT_0132;
    const ccomplex_t IT_0356 = (-12)*IT_0131;
    const ccomplex_t IT_0357 = IT_0354 + IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = 6*IT_0306;
    const ccomplex_t IT_0359 = s_14*IT_0000;
    const ccomplex_t IT_0360 = s_13*s_34;
    const ccomplex_t IT_0361 = -IT_0360;
    const ccomplex_t IT_0362 = IT_0359 + IT_0361;
    const ccomplex_t IT_0363 = (-6)*IT_0362;
    const ccomplex_t IT_0364 = IT_0130*IT_0321;
    const ccomplex_t IT_0365 = (-3)*IT_0364;
    const ccomplex_t IT_0366 = (-6)*IT_0365;
    const ccomplex_t IT_0367 = conj(IT_0308)*IT_0313 + (-0.166666666666667)
      *conj(IT_0202)*IT_0314 + (-0.166666666666667)*IT_0298*IT_0315 + (
      -0.166666666666667)*conj(IT_0128)*IT_0347 + (-0.166666666666667)*conj
      (IT_0299)*IT_0357 + (-0.166666666666667)*IT_0337*IT_0358 + (
      -0.166666666666667)*conj(IT_0264)*IT_0363 + (-0.166666666666667)*conj
      (IT_0307)*IT_0366;
    const ccomplex_t IT_0368 = 6*IT_0362;
    const ccomplex_t IT_0369 = (-6)*conj(IT_0087);
    const ccomplex_t IT_0370 = (-6)*IT_0298;
    const ccomplex_t IT_0371 = IT_0292*IT_0345;
    const ccomplex_t IT_0372 = (-18)*IT_0087*IT_0321;
    const ccomplex_t IT_0373 = IT_0286*IT_0349;
    const ccomplex_t IT_0374 = IT_0138*IT_0327;
    const ccomplex_t IT_0375 = (-0.166666666666667)*conj(IT_0308);
    const ccomplex_t IT_0376 = conj(IT_0299)*IT_0313 + (-0.166666666666667)
      *conj(IT_0202)*IT_0324 + (-0.166666666666667)*IT_0337*IT_0341 + (
      -0.166666666666667)*conj(IT_0286)*IT_0342 + (-0.166666666666667)*conj
      (IT_0264)*IT_0343 + (-0.166666666666667)*conj(IT_0292)*IT_0344 + (
      -0.166666666666667)*conj(IT_0128)*IT_0349 + (-0.166666666666667)*IT_0144
      *IT_0369 + IT_0318*IT_0375;
    const ccomplex_t IT_0377 = s_13*IT_0290;
    const ccomplex_t IT_0378 = IT_0290*IT_0321;
    const ccomplex_t IT_0379 = (-6)*IT_0378;
    const ccomplex_t IT_0380 = (-0.166666666666667)*IT_0206*IT_0315 + (
      -0.166666666666667)*conj(IT_0202)*IT_0320 + (-0.166666666666667)*conj
      (IT_0286)*IT_0344 + (-0.166666666666667)*conj(IT_0128)*IT_0345 + (
      -0.166666666666667)*IT_0337*IT_0346 + (-0.166666666666667)*conj(IT_0299)
      *IT_0363 + IT_0343*IT_0375 + conj(IT_0264)*IT_0377 + (-0.166666666666667)
      *conj(IT_0292)*IT_0379;
    const ccomplex_t IT_0381 = (-0.166666666666667)*IT_0291*IT_0315 + (
      -0.166666666666667)*conj(IT_0202)*IT_0323 + (-0.166666666666667)*conj
      (IT_0286)*IT_0343 + (-0.166666666666667)*conj(IT_0128)*IT_0346 + (
      -0.166666666666667)*conj(IT_0307)*IT_0368 + IT_0344*IT_0375 + conj(IT_0292
      )*IT_0377 + (-0.166666666666667)*conj(IT_0264)*IT_0379;
    const ccomplex_t IT_0382 = conj(IT_0292)*IT_0345;
    const ccomplex_t IT_0383 = conj(IT_0308)*IT_0341;
    const ccomplex_t IT_0384 = conj(IT_0264)*IT_0346;
    const ccomplex_t IT_0385 = (-18)*conj(IT_0087)*IT_0321;
    const ccomplex_t IT_0386 = conj(IT_0286)*IT_0349;
    const ccomplex_t IT_0387 = conj(IT_0202)*IT_0338;
    const ccomplex_t IT_0388 = conj(IT_0307)*IT_0370;
    const ccomplex_t IT_0389 = conj(IT_0128)*IT_0348;
    const ccomplex_t IT_0390 = conj(IT_0299)*IT_0358;
    const ccomplex_t IT_0391 = IT_0087*(conj(IT_0087) + conj(IT_0128))*IT_0138
       + (IT_0144*conj(IT_0202) + IT_0206*conj(IT_0264) + IT_0265*conj(IT_0286) 
      + IT_0291*conj(IT_0292) + IT_0298*conj(IT_0299))*IT_0300 + (IT_0306*conj
      (IT_0307) + IT_0144*conj(IT_0308))*IT_0309 + IT_0202*(conj(IT_0299)
      *IT_0314 + IT_0144*IT_0315 + conj(IT_0128)*IT_0316 + conj(IT_0202)*IT_0318
       + conj(IT_0264)*IT_0320 + conj(IT_0292)*IT_0323 + conj(IT_0308)*IT_0324 +
       IT_0337*IT_0338 + conj(IT_0286)*IT_0340) + IT_0286*(conj(IT_0307)*IT_0314
       + IT_0265*IT_0315 + conj(IT_0286)*IT_0318 + conj(IT_0202)*IT_0340 + conj
      (IT_0128)*IT_0341 + conj(IT_0308)*IT_0342 + conj(IT_0292)*IT_0343 + conj
      (IT_0264)*IT_0344) + IT_0128*((conj(IT_0087) + conj(IT_0128))*IT_0138 +
       conj(IT_0202)*IT_0316 + conj(IT_0286)*IT_0341 + conj(IT_0264)*IT_0345 +
       conj(IT_0292)*IT_0346 + conj(IT_0299)*IT_0347 + IT_0337*IT_0348 + conj
      (IT_0308)*IT_0349 + conj(IT_0307)*IT_0350) + (-6)*IT_0299*IT_0367 +
       IT_0307*(conj(IT_0286)*IT_0314 + conj(IT_0128)*IT_0350 + conj(IT_0307)
      *IT_0357 + conj(IT_0299)*IT_0366 + conj(IT_0292)*IT_0368 + IT_0306*IT_0369
       + IT_0337*IT_0370) + IT_0337*(IT_0371 + IT_0372 + IT_0373 + IT_0374) + (
      -6)*IT_0308*IT_0376 + (-6)*IT_0264*IT_0380 + (-6)*IT_0292*IT_0381 +
       IT_0327*(IT_0382 + IT_0383 + IT_0384 + IT_0385 + IT_0386 + IT_0387 +
       IT_0388 + IT_0389 + IT_0390) + IT_0336*(IT_0138*IT_0337 + IT_0382 +
       IT_0383 + IT_0384 + IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389 +
       IT_0390);
    return create_ccomplex_return(IT_0391);
}

