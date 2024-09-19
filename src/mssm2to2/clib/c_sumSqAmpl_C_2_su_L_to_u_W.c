#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_su_L_to_u_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_su_L_to_u_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_sl = param->Gamma_sl;
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
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = pow(m_u, 2);
    const ccomplex_t IT_0010 = pow(m_C_2, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0001;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0023*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0001*IT_0004*IT_0006*IT_0026*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0025*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0023*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0001*IT_0004*IT_0006*IT_0026*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0034*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0044 = IT_0001*IT_0043;
    const ccomplex_t IT_0045 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + (-0.5)*IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = IT_0042*IT_0049;
    const ccomplex_t IT_0051 = pow(m_su_L, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0009 + IT_0051 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0009 + IT_0051 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0059 = IT_0001*IT_0058;
    const ccomplex_t IT_0060 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0061 = IT_0001*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + (-0.5)*IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0057*IT_0064;
    const ccomplex_t IT_0066 = IT_0055*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0009 + IT_0051 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0069 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0072 = IT_0001*IT_0071;
    const ccomplex_t IT_0073 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0074 = IT_0001*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + (-0.5)*IT_0075);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0070*IT_0077;
    const ccomplex_t IT_0079 = IT_0068*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0001*IT_0004*IT_0006;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0084 = IT_0001*IT_0083;
    const ccomplex_t IT_0085 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0086 = IT_0001*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0084 + (-0.5)*IT_0087);
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = IT_0082*IT_0089;
    const ccomplex_t IT_0091 = cpow((-2)*s_23 + IT_0009 + IT_0051 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = IT_0014 + IT_0022 + IT_0032 + IT_0040 + 2
      *IT_0054 + 2*IT_0067 + 2*IT_0080 + 2*IT_0093;
    const ccomplex_t IT_0095 = pow(m_W, -2);
    const ccomplex_t IT_0096 = pow(s_14, 2);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = -IT_0010;
    const ccomplex_t IT_0099 = IT_0097 + IT_0098;
    const ccomplex_t IT_0100 = s_13*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0102 = cos(theta_W);
    const ccomplex_t IT_0103 = cpow(IT_0102, -1);
    const ccomplex_t IT_0104 = conj(N_B1)*e_em;
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = conj(N_W1)*e_em;
    const ccomplex_t IT_0108 = IT_0001*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + 3*IT_0109);
    const ccomplex_t IT_0111 = 0.166666666666667*IT_0110;
    const ccomplex_t IT_0112 = IT_0049*IT_0111;
    const ccomplex_t IT_0113 = IT_0052*IT_0101*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0115 = conj(N_B2)*e_em;
    const ccomplex_t IT_0116 = IT_0103*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W2)*e_em;
    const ccomplex_t IT_0119 = IT_0001*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0117 + 3*IT_0120);
    const ccomplex_t IT_0122 = 0.166666666666667*IT_0121;
    const ccomplex_t IT_0123 = IT_0064*IT_0122;
    const ccomplex_t IT_0124 = IT_0055*IT_0114*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0126 = conj(N_B3)*e_em;
    const ccomplex_t IT_0127 = IT_0103*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W3)*e_em;
    const ccomplex_t IT_0130 = IT_0001*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*(IT_0128 + 3*IT_0131);
    const ccomplex_t IT_0133 = 0.166666666666667*IT_0132;
    const ccomplex_t IT_0134 = IT_0077*IT_0133;
    const ccomplex_t IT_0135 = IT_0068*IT_0125*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0137 = conj(N_B4)*e_em;
    const ccomplex_t IT_0138 = IT_0103*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = conj(N_W4)*e_em;
    const ccomplex_t IT_0141 = IT_0001*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 3*IT_0142);
    const ccomplex_t IT_0144 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_0145 = IT_0089*IT_0144;
    const ccomplex_t IT_0146 = IT_0091*IT_0136*IT_0145;
    const ccomplex_t IT_0147 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0148 = IT_0001*IT_0147;
    const ccomplex_t IT_0149 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0150 = IT_0001*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0148 + 0.5*IT_0151);
    const ccomplex_t IT_0153 = -IT_0152;
    const ccomplex_t IT_0154 = IT_0111*IT_0153;
    const ccomplex_t IT_0155 = IT_0052*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = m_C_2*IT_0156;
    const ccomplex_t IT_0158 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0159 = IT_0001*IT_0158;
    const ccomplex_t IT_0160 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0161 = IT_0001*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + 0.5*IT_0162);
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0122*IT_0164;
    const ccomplex_t IT_0166 = IT_0055*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = m_C_2*IT_0167;
    const ccomplex_t IT_0169 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0170 = IT_0001*IT_0169;
    const ccomplex_t IT_0171 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0172 = IT_0001*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0170 + 0.5*IT_0173);
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = IT_0133*IT_0175;
    const ccomplex_t IT_0177 = IT_0068*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = m_C_2*IT_0178;
    const ccomplex_t IT_0180 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0183 = IT_0001*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0181 + 0.5*IT_0184);
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = IT_0144*IT_0186;
    const ccomplex_t IT_0188 = IT_0091*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = m_C_2*IT_0189;
    const ccomplex_t IT_0191 = -IT_0113 + -IT_0124 + -IT_0135 + -IT_0146 +
       IT_0157 + IT_0168 + IT_0179 + IT_0190;
    const ccomplex_t IT_0192 = m_u*IT_0099;
    const ccomplex_t IT_0193 = -IT_0054 + -IT_0067 + -IT_0080 + -IT_0093;
    const ccomplex_t IT_0194 = s_13*s_14;
    const ccomplex_t IT_0195 = s_34*IT_0010;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = -IT_0014 + -IT_0022 + -IT_0032 + -IT_0040;
    const ccomplex_t IT_0199 = s_14*IT_0009;
    const ccomplex_t IT_0200 = s_13*s_34;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0199 + IT_0201;
    const ccomplex_t IT_0203 = s_14*s_34*IT_0095;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = s_13 + IT_0204;
    const ccomplex_t IT_0206 = s_13*IT_0205;
    const ccomplex_t IT_0207 = (-6)*conj(IT_0094);
    const ccomplex_t IT_0208 = m_u*IT_0205;
    const ccomplex_t IT_0209 = (-6)*IT_0198;
    const ccomplex_t IT_0210 = (-6)*IT_0094;
    const ccomplex_t IT_0211 = pow(s_34, 2);
    const ccomplex_t IT_0212 = IT_0095*IT_0211;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = IT_0009 + IT_0213;
    const ccomplex_t IT_0215 = s_13*IT_0214;
    const ccomplex_t IT_0216 = pow(m_W, 2);
    const ccomplex_t IT_0217 = s_14*IT_0216;
    const ccomplex_t IT_0218 = IT_0095*IT_0217;
    const ccomplex_t IT_0219 = 0.5*IT_0218;
    const ccomplex_t IT_0220 = s_14 + IT_0219;
    const ccomplex_t IT_0221 = m_u*IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0221;
    const ccomplex_t IT_0223 = IT_0156 + IT_0167 + IT_0178 + IT_0189;
    const ccomplex_t IT_0224 = s_34*IT_0216;
    const ccomplex_t IT_0225 = IT_0095*IT_0224;
    const ccomplex_t IT_0226 = (-0.25)*IT_0225;
    const ccomplex_t IT_0227 = s_34 + IT_0226;
    const ccomplex_t IT_0228 = m_C_2*IT_0227;
    const ccomplex_t IT_0229 = (-4)*IT_0228;
    const ccomplex_t IT_0230 = 12*s_13;
    const ccomplex_t IT_0231 = s_13*IT_0216;
    const ccomplex_t IT_0232 = s_14*s_34;
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = IT_0231 + IT_0233;
    const ccomplex_t IT_0235 = IT_0095*IT_0234;
    const ccomplex_t IT_0236 = (-6)*IT_0235;
    const ccomplex_t IT_0237 = IT_0230 + IT_0236;
    const ccomplex_t IT_0238 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0239 = IT_0057*IT_0164;
    const ccomplex_t IT_0240 = IT_0055*IT_0114*IT_0239;
    const ccomplex_t IT_0241 = IT_0042*IT_0153;
    const ccomplex_t IT_0242 = IT_0052*IT_0101*IT_0241;
    const ccomplex_t IT_0243 = IT_0082*IT_0186;
    const ccomplex_t IT_0244 = IT_0091*IT_0136*IT_0243;
    const ccomplex_t IT_0245 = m_C_2*IT_0054;
    const ccomplex_t IT_0246 = m_C_2*IT_0093;
    const ccomplex_t IT_0247 = IT_0070*IT_0175;
    const ccomplex_t IT_0248 = IT_0068*IT_0125*IT_0247;
    const ccomplex_t IT_0249 = m_C_2*IT_0067;
    const ccomplex_t IT_0250 = m_C_2*IT_0080;
    const ccomplex_t IT_0251 = IT_0240 + IT_0242 + IT_0244 + -IT_0245 + 
      -IT_0246 + IT_0248 + -IT_0249 + -IT_0250;
    const ccomplex_t IT_0252 = m_u*m_C_2;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0254 = IT_0019*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = IT_0016*IT_0255;
    const ccomplex_t IT_0257 = IT_0001*IT_0026;
    const ccomplex_t IT_0258 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0259 = IT_0257*IT_0258;
    const ccomplex_t IT_0260 = cos(beta);
    const ccomplex_t IT_0261 = cpow(IT_0260, -1);
    const ccomplex_t IT_0262 = IT_0001*IT_0026*IT_0261;
    const ccomplex_t IT_0263 = m_b*conj(U_d2)*e_em*IT_0004*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0264 = IT_0262*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*(IT_0259 + (-0.5)*IT_0265);
    const ccomplex_t IT_0267 = IT_0029*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = IT_0025*IT_0268;
    const ccomplex_t IT_0270 = m_b*conj(U_d2)*e_em*IT_0004*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0271 = IT_0262*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0274 = IT_0257*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*(IT_0272 + (-2)*IT_0274);
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0037*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = IT_0034*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0281 = IT_0011*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = IT_0003*IT_0282;
    const ccomplex_t IT_0284 = (-2)*IT_0156 + (-2)*IT_0167 + (-2)*IT_0178 + (
      -2)*IT_0189 + -IT_0256 + -IT_0269 + -IT_0279 + -IT_0283;
    const ccomplex_t IT_0285 = IT_0256 + IT_0269 + IT_0279 + IT_0283;
    const ccomplex_t IT_0286 = m_C_2*IT_0205;
    const ccomplex_t IT_0287 = m_C_2*IT_0214;
    const ccomplex_t IT_0288 = 6*IT_0222*conj(IT_0223) + 6*conj(IT_0193)
      *IT_0229 + IT_0237*conj(IT_0251) + (-18)*conj(IT_0191)*IT_0252 + 6*IT_0192
      *conj(IT_0284) + (-6)*IT_0208*conj(IT_0285) + IT_0207*IT_0286 + IT_0238
      *IT_0287;
    const ccomplex_t IT_0289 = IT_0099*IT_0252;
    const ccomplex_t IT_0290 = IT_0205*IT_0252;
    const ccomplex_t IT_0291 = (-6)*IT_0191;
    const ccomplex_t IT_0292 = IT_0216*IT_0252;
    const ccomplex_t IT_0293 = (-3)*IT_0292;
    const ccomplex_t IT_0294 = pow(m_W, 4);
    const ccomplex_t IT_0295 = s_13*IT_0294;
    const ccomplex_t IT_0296 = IT_0095*IT_0295;
    const ccomplex_t IT_0297 = 24*IT_0232;
    const ccomplex_t IT_0298 = (-12)*IT_0231;
    const ccomplex_t IT_0299 = 6*IT_0296 + IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = 6*IT_0193;
    const ccomplex_t IT_0301 = 6*IT_0223;
    const ccomplex_t IT_0302 = (-6)*IT_0285;
    const ccomplex_t IT_0303 = (-6)*conj(IT_0191);
    const ccomplex_t IT_0304 = 0.166666666666667*IT_0290;
    const ccomplex_t IT_0305 = IT_0214*IT_0252;
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0305;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0287;
    const ccomplex_t IT_0308 = 6*IT_0094*(conj(IT_0094)*IT_0100 + conj(IT_0191
      )*IT_0192 + conj(IT_0193)*IT_0197) + 6*IT_0198*(conj(IT_0193)*IT_0202 +
       0.166666666666667*IT_0206*IT_0207) + conj(IT_0191)*IT_0208*IT_0209 + conj
      (IT_0198)*(IT_0206*IT_0210 + IT_0209*IT_0215) + 6*IT_0191*(conj(IT_0094)
      *IT_0192 + conj(IT_0193)*IT_0222 + conj(IT_0223)*IT_0229 +
       0.166666666666667*conj(IT_0191)*IT_0237 + 0.166666666666667*IT_0208
      *IT_0238) + IT_0251*IT_0288 + 6*conj(IT_0284)*(IT_0197*IT_0223 + -IT_0206
      *IT_0285 + IT_0094*IT_0289 + 0.166666666666667*IT_0209*IT_0290 +
       0.166666666666667*IT_0286*IT_0291) + (conj(IT_0094)*IT_0197 + conj
      (IT_0198)*IT_0202 + conj(IT_0191)*IT_0222 + conj(IT_0223)*IT_0293 +
       0.166666666666667*conj(IT_0193)*IT_0299)*IT_0300 + (conj(IT_0191)*IT_0229
       + conj(IT_0193)*IT_0293 + 0.166666666666667*conj(IT_0223)*IT_0299)
      *IT_0301 + conj(IT_0251)*((-18)*IT_0191*IT_0252 + 6*IT_0192*IT_0284 +
       IT_0210*IT_0286 + IT_0209*IT_0287 + IT_0229*IT_0300 + IT_0222*IT_0301 +
       IT_0208*IT_0302) + 6*IT_0284*(IT_0197*conj(IT_0223) + IT_0100*conj
      (IT_0284) + -IT_0206*conj(IT_0285) + conj(IT_0094)*IT_0289 +
       0.166666666666667*IT_0286*IT_0303 + IT_0238*IT_0304) + 6*conj(IT_0285)*
      (IT_0202*IT_0223 + 0.166666666666667*IT_0215*IT_0302 + IT_0210*IT_0304 +
       IT_0209*IT_0306 + IT_0291*IT_0307) + 6*IT_0285*(IT_0202*conj(IT_0223) +
       IT_0207*IT_0304 + IT_0238*IT_0306 + IT_0303*IT_0307);
    return create_ccomplex_return(IT_0308);
}

