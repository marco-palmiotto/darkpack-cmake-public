#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_b_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_b_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_C_2, 2);
    const ccomplex_t IT_0001 = pow(m_sb_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*e_em
      *IT_0003*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0003*IT_0006*IT_0008*U_sb_01;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = pow(m_b, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0004*IT_0008;
    const ccomplex_t IT_0019 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 2*IT_0020;
    const ccomplex_t IT_0022 = cos(theta_W);
    const ccomplex_t IT_0023 = cpow(IT_0022, -1);
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0029 = IT_0018*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = IT_0027 + IT_0030;
    const ccomplex_t IT_0032 = IT_0021 + IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = conj(N_B1)*e_em*U_sb_01;
    const ccomplex_t IT_0036 = IT_0023*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = conj(N_W1)*e_em*U_sb_01;
    const ccomplex_t IT_0039 = IT_0008*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = cpow(IT_0004, -1);
    const ccomplex_t IT_0042 = IT_0008*IT_0041;
    const ccomplex_t IT_0043 = m_b*conj(N_d1)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0037 + (-3)*IT_0040 + 3
      *IT_0045);
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0034*IT_0047;
    const ccomplex_t IT_0049 = IT_0017*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0052 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0053 = IT_0018*IT_0052;
    const ccomplex_t IT_0054 = 2*IT_0053;
    const ccomplex_t IT_0055 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0056 = IT_0024*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0059 = IT_0018*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = IT_0057 + IT_0060;
    const ccomplex_t IT_0062 = IT_0054 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B2)*e_em*U_sb_01;
    const ccomplex_t IT_0066 = IT_0023*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W2)*e_em*U_sb_01;
    const ccomplex_t IT_0069 = IT_0008*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_b*conj(N_d2)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0072 = IT_0042*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + (-3)*IT_0070 + 3
      *IT_0073);
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = IT_0051*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0080 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0081 = IT_0018*IT_0080;
    const ccomplex_t IT_0082 = 2*IT_0081;
    const ccomplex_t IT_0083 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0084 = IT_0024*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0087 = IT_0018*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = IT_0082 + IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = conj(N_B3)*e_em*U_sb_01;
    const ccomplex_t IT_0094 = IT_0023*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W3)*e_em*U_sb_01;
    const ccomplex_t IT_0097 = IT_0008*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = m_b*conj(N_d3)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0100 = IT_0042*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0095 + (-3)*IT_0098 + 3
      *IT_0101);
    const ccomplex_t IT_0103 = 0.166666666666667*IT_0102;
    const ccomplex_t IT_0104 = IT_0092*IT_0103;
    const ccomplex_t IT_0105 = IT_0079*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0001 + IT_0016 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0108 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0109 = IT_0018*IT_0108;
    const ccomplex_t IT_0110 = 2*IT_0109;
    const ccomplex_t IT_0111 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0112 = IT_0024*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0115 = IT_0018*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = IT_0113 + IT_0116;
    const ccomplex_t IT_0118 = IT_0110 + IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B4)*e_em*U_sb_01;
    const ccomplex_t IT_0122 = IT_0023*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W4)*e_em*U_sb_01;
    const ccomplex_t IT_0125 = IT_0008*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = m_b*conj(N_d4)*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0128 = IT_0042*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0123 + (-3)*IT_0126 + 3
      *IT_0129);
    const ccomplex_t IT_0131 = 0.166666666666667*IT_0130;
    const ccomplex_t IT_0132 = IT_0120*IT_0131;
    const ccomplex_t IT_0133 = IT_0107*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = pow(m_c, 2);
    const ccomplex_t IT_0136 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + IT_0135 
      + -reg_prop, -1);
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0003*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*V_u2
      *e_em*IT_0003*IT_0006*IT_0008*U_sb_01;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = IT_0136*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = pow(m_u, 2);
    const ccomplex_t IT_0145 = cpow((-2)*s_12 + -IT_0000 + -IT_0001 + IT_0144 
      + -reg_prop, -1);
    const ccomplex_t IT_0146 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_u*e_em
      *IT_0003*IT_0004*IT_0006*IT_0008*IT_0146*V_ub_mod;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0006*IT_0008*IT_0149*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0148*IT_0151;
    const ccomplex_t IT_0153 = IT_0145*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = -IT_0015 + -IT_0050 + -IT_0078 + -IT_0106 + 
      -IT_0134 + -IT_0143 + -IT_0154;
    const ccomplex_t IT_0156 = m_b*m_C_2*IT_0001;
    const ccomplex_t IT_0157 = N_B1*e_em*U_sb_11;
    const ccomplex_t IT_0158 = IT_0023*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = m_b*N_d1*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0161 = IT_0042*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + 1.5*IT_0162);
    const ccomplex_t IT_0164 = (-0.333333333333333)*IT_0163;
    const ccomplex_t IT_0165 = IT_0005*IT_0008;
    const ccomplex_t IT_0166 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0167 = IT_0165*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = IT_0005*IT_0023;
    const ccomplex_t IT_0170 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0171 = IT_0169*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = IT_0168 + IT_0172;
    const ccomplex_t IT_0174 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0175 = IT_0165*IT_0174;
    const ccomplex_t IT_0176 = (-2)*IT_0175;
    const ccomplex_t IT_0177 = IT_0173 + IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = IT_0164*IT_0179;
    const ccomplex_t IT_0181 = IT_0017*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = V_cb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0184 = IT_0008*IT_0183;
    const ccomplex_t IT_0185 = m_b*conj(U_d2)*V_cb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0186 = IT_0042*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0184 + (-0.5)*IT_0187);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1.4142135623731)*m_b*V_cb*e_em
      *IT_0003*IT_0005*IT_0008*IT_0041;
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = IT_0188*IT_0190;
    const ccomplex_t IT_0192 = IT_0136*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0195 = IT_0008*IT_0194;
    const ccomplex_t IT_0196 = m_b*conj(U_d2)*V_tb*e_em*IT_0003*U_sb_11;
    const ccomplex_t IT_0197 = IT_0042*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + (-0.5)*IT_0198);
    const ccomplex_t IT_0200 = (0 + _Complex_I*1.4142135623731)*m_b*V_tb*e_em
      *IT_0003*IT_0005*IT_0008*IT_0041;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = IT_0199*IT_0201;
    const ccomplex_t IT_0203 = IT_0002*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = IT_0008*IT_0149;
    const ccomplex_t IT_0206 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0207 = IT_0205*IT_0206;
    const ccomplex_t IT_0208 = IT_0008*IT_0041*IT_0149;
    const ccomplex_t IT_0209 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0210 = IT_0208*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*(IT_0207 + (-0.5)*IT_0211);
    const ccomplex_t IT_0213 = (0 + _Complex_I*1.4142135623731)*m_b*e_em
      *IT_0003*IT_0005*IT_0008*IT_0041*IT_0146*V_ub_mod;
    const ccomplex_t IT_0214 = 0.5*IT_0213;
    const ccomplex_t IT_0215 = IT_0212*IT_0214;
    const ccomplex_t IT_0216 = IT_0145*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0219 = IT_0165*IT_0218;
    const ccomplex_t IT_0220 = (-2)*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0222 = IT_0169*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0225 = IT_0165*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = IT_0223 + IT_0226;
    const ccomplex_t IT_0228 = IT_0220 + IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = 0.5*IT_0229;
    const ccomplex_t IT_0231 = N_B4*e_em*U_sb_11;
    const ccomplex_t IT_0232 = IT_0023*IT_0231;
    const ccomplex_t IT_0233 = 1.4142135623731*IT_0232;
    const ccomplex_t IT_0234 = m_b*N_d4*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0235 = IT_0042*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*(IT_0233 + 1.5*IT_0236);
    const ccomplex_t IT_0238 = (-0.333333333333333)*IT_0237;
    const ccomplex_t IT_0239 = IT_0230*IT_0238;
    const ccomplex_t IT_0240 = IT_0107*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0243 = IT_0165*IT_0242;
    const ccomplex_t IT_0244 = (-2)*IT_0243;
    const ccomplex_t IT_0245 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0246 = IT_0169*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0249 = IT_0165*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = IT_0247 + IT_0250;
    const ccomplex_t IT_0252 = IT_0244 + IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = N_B2*e_em*U_sb_11;
    const ccomplex_t IT_0256 = IT_0023*IT_0255;
    const ccomplex_t IT_0257 = 1.4142135623731*IT_0256;
    const ccomplex_t IT_0258 = m_b*N_d2*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0259 = IT_0042*IT_0258;
    const ccomplex_t IT_0260 = 1.4142135623731*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*(IT_0257 + 1.5*IT_0260);
    const ccomplex_t IT_0262 = (-0.333333333333333)*IT_0261;
    const ccomplex_t IT_0263 = IT_0254*IT_0262;
    const ccomplex_t IT_0264 = IT_0051*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0267 = IT_0165*IT_0266;
    const ccomplex_t IT_0268 = (-2)*IT_0267;
    const ccomplex_t IT_0269 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0270 = IT_0169*IT_0269;
    const ccomplex_t IT_0271 = 1.4142135623731*IT_0270;
    const ccomplex_t IT_0272 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0273 = IT_0165*IT_0272;
    const ccomplex_t IT_0274 = 1.4142135623731*IT_0273;
    const ccomplex_t IT_0275 = IT_0271 + IT_0274;
    const ccomplex_t IT_0276 = IT_0268 + IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = 0.5*IT_0277;
    const ccomplex_t IT_0279 = N_B3*e_em*U_sb_11;
    const ccomplex_t IT_0280 = IT_0023*IT_0279;
    const ccomplex_t IT_0281 = 1.4142135623731*IT_0280;
    const ccomplex_t IT_0282 = m_b*N_d3*e_em*IT_0003*U_sb_01;
    const ccomplex_t IT_0283 = IT_0042*IT_0282;
    const ccomplex_t IT_0284 = 1.4142135623731*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*(IT_0281 + 1.5*IT_0284);
    const ccomplex_t IT_0286 = (-0.333333333333333)*IT_0285;
    const ccomplex_t IT_0287 = IT_0278*IT_0286;
    const ccomplex_t IT_0288 = IT_0079*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = -IT_0182 + IT_0193 + IT_0204 + IT_0217 + 
      -IT_0241 + -IT_0265 + -IT_0289;
    const ccomplex_t IT_0291 = m_b*s_12;
    const ccomplex_t IT_0292 = 6*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0294 = IT_0151*IT_0214;
    const ccomplex_t IT_0295 = IT_0145*IT_0293*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0297 = IT_0034*IT_0164;
    const ccomplex_t IT_0298 = IT_0017*IT_0296*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0300 = IT_0092*IT_0286;
    const ccomplex_t IT_0301 = IT_0079*IT_0299*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0303 = IT_0012*IT_0201;
    const ccomplex_t IT_0304 = IT_0002*IT_0302*IT_0303;
    const ccomplex_t IT_0305 = m_b*IT_0050;
    const ccomplex_t IT_0306 = m_b*IT_0078;
    const ccomplex_t IT_0307 = m_b*IT_0106;
    const ccomplex_t IT_0308 = m_b*IT_0134;
    const ccomplex_t IT_0309 = m_C_2*IT_0193;
    const ccomplex_t IT_0310 = m_C_2*IT_0204;
    const ccomplex_t IT_0311 = m_C_2*IT_0217;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0313 = IT_0064*IT_0262;
    const ccomplex_t IT_0314 = IT_0051*IT_0312*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0316 = IT_0120*IT_0238;
    const ccomplex_t IT_0317 = IT_0107*IT_0315*IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0319 = IT_0140*IT_0190;
    const ccomplex_t IT_0320 = IT_0136*IT_0318*IT_0319;
    const ccomplex_t IT_0321 = -IT_0295 + -IT_0298 + -IT_0301 + -IT_0304 +
       IT_0305 + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311 + 
      -IT_0314 + -IT_0317 + -IT_0320;
    const ccomplex_t IT_0322 = 0.166666666666667*conj(IT_0321);
    const ccomplex_t IT_0323 = s_23*m_C_2;
    const ccomplex_t IT_0324 = 6*IT_0323;
    const ccomplex_t IT_0325 = IT_0103*IT_0278;
    const ccomplex_t IT_0326 = IT_0079*IT_0299*IT_0325;
    const ccomplex_t IT_0327 = IT_0047*IT_0179;
    const ccomplex_t IT_0328 = IT_0017*IT_0296*IT_0327;
    const ccomplex_t IT_0329 = IT_0131*IT_0230;
    const ccomplex_t IT_0330 = IT_0107*IT_0315*IT_0329;
    const ccomplex_t IT_0331 = IT_0138*IT_0188;
    const ccomplex_t IT_0332 = IT_0136*IT_0318*IT_0331;
    const ccomplex_t IT_0333 = IT_0010*IT_0199;
    const ccomplex_t IT_0334 = IT_0002*IT_0302*IT_0333;
    const ccomplex_t IT_0335 = m_b*IT_0182;
    const ccomplex_t IT_0336 = IT_0075*IT_0254;
    const ccomplex_t IT_0337 = IT_0051*IT_0312*IT_0336;
    const ccomplex_t IT_0338 = m_b*IT_0289;
    const ccomplex_t IT_0339 = m_b*IT_0241;
    const ccomplex_t IT_0340 = m_C_2*IT_0154;
    const ccomplex_t IT_0341 = IT_0148*IT_0212;
    const ccomplex_t IT_0342 = IT_0145*IT_0293*IT_0341;
    const ccomplex_t IT_0343 = m_b*IT_0265;
    const ccomplex_t IT_0344 = m_C_2*IT_0143;
    const ccomplex_t IT_0345 = m_C_2*IT_0015;
    const ccomplex_t IT_0346 = -IT_0326 + -IT_0328 + -IT_0330 + IT_0332 +
       IT_0334 + IT_0335 + -IT_0337 + IT_0338 + IT_0339 + -IT_0340 + IT_0342 +
       IT_0343 + -IT_0344 + -IT_0345;
    const ccomplex_t IT_0347 = 0.166666666666667*conj(IT_0346);
    const ccomplex_t IT_0348 = s_13*IT_0001;
    const ccomplex_t IT_0349 = (-6)*IT_0348;
    const ccomplex_t IT_0350 = s_12*s_23;
    const ccomplex_t IT_0351 = 12*IT_0350;
    const ccomplex_t IT_0352 = IT_0349 + IT_0351;
    const ccomplex_t IT_0353 = 0.166666666666667*conj(IT_0155);
    const ccomplex_t IT_0354 = 0.166666666666667*conj(IT_0290);
    const ccomplex_t IT_0355 = m_b*m_C_2;
    const ccomplex_t IT_0356 = 6*IT_0355;
    const ccomplex_t IT_0357 = 6*IT_0155*(IT_0156*conj(IT_0290) + IT_0292
      *IT_0322 + IT_0324*IT_0347 + IT_0352*IT_0353) + 6*IT_0290*(conj(IT_0155)
      *IT_0156 + IT_0322*IT_0324 + IT_0292*IT_0347 + IT_0352*IT_0354) + 6
      *IT_0346*(s_13*conj(IT_0346) + IT_0324*IT_0353 + IT_0292*IT_0354 + IT_0322
      *IT_0356) + 6*IT_0321*(s_13*conj(IT_0321) + IT_0292*IT_0353 + IT_0324
      *IT_0354 + IT_0347*IT_0356);
    return create_ccomplex_return(IT_0357);
}

