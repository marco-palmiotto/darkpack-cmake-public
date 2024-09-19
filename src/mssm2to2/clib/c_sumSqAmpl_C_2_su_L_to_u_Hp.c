#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_su_L_to_u_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_su_L_to_u_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
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
    const ccomplex_t IT_0000 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = pow(m_u, 2);
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = cpow(IT_0001, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0000
      *IT_0008*(m_W*IT_0001*IT_0002*conj(U_sb_01) + (-0.5)*IT_0003*(IT_0001
      *IT_0004*IT_0005*conj(U_sb_01) + m_b*(m_b*IT_0002*IT_0006*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0005*IT_0008*IT_0011*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0004 + IT_0015 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0008;
    const ccomplex_t IT_0023 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 2*IT_0024;
    const ccomplex_t IT_0026 = cos(theta_W);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = IT_0001*IT_0027;
    const ccomplex_t IT_0029 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0033 = IT_0022*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = IT_0031 + IT_0034;
    const ccomplex_t IT_0036 = IT_0025 + IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0021*IT_0038;
    const ccomplex_t IT_0040 = pow(m_su_L, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0004 + IT_0040 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0019*IT_0039*IT_0041;
    const ccomplex_t IT_0043 = pow(m_s, 2);
    const ccomplex_t IT_0044 = (0 + _Complex_I*(-2.82842712474619))*V_us*e_em
      *IT_0008*(IT_0001*(m_W*IT_0002 + (-0.5)*IT_0003*IT_0004*IT_0005) + -1./2
      *IT_0002*IT_0003*IT_0006*IT_0043);
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0047 = 0.5*IT_0046;
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0004 + IT_0015 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0056 = IT_0022*IT_0055;
    const ccomplex_t IT_0057 = 2*IT_0056;
    const ccomplex_t IT_0058 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0059 = IT_0028*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0062 = IT_0022*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0060 + IT_0063;
    const ccomplex_t IT_0065 = IT_0057 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = IT_0054*IT_0067;
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0004 + IT_0040 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0052*IT_0068*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0004 + IT_0040 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0073 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0076 = IT_0022*IT_0075;
    const ccomplex_t IT_0077 = 2*IT_0076;
    const ccomplex_t IT_0078 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0079 = IT_0028*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0082 = IT_0022*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = IT_0080 + IT_0083;
    const ccomplex_t IT_0085 = IT_0077 + IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0074*IT_0087;
    const ccomplex_t IT_0089 = IT_0071*IT_0072*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_13 + IT_0004 + IT_0015 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0091 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0000
      *IT_0008*(m_W*IT_0001*IT_0002*conj(U_sb_00) + (-0.5)*IT_0003*(IT_0001
      *IT_0004*IT_0005*conj(U_sb_00) + m_b*(m_b*IT_0002*IT_0006*conj(U_sb_00) +
       mu_h*conj(U_sb_10))))*V_ub_mod;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0005*IT_0008*IT_0011*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = IT_0090*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = conj(N_B1)*e_em;
    const ccomplex_t IT_0099 = IT_0027*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = conj(N_W1)*e_em;
    const ccomplex_t IT_0102 = IT_0008*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + 3*IT_0103);
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = IT_0038*IT_0105;
    const ccomplex_t IT_0107 = IT_0041*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = m_u*IT_0108;
    const ccomplex_t IT_0110 = pow(m_d, 2);
    const ccomplex_t IT_0111 = (0 + _Complex_I*(-2.82842712474619))*V_ud*e_em
      *IT_0008*(IT_0001*(m_W*IT_0002 + (-0.5)*IT_0003*IT_0004*IT_0005) + -1./2
      *IT_0002*IT_0003*IT_0006*IT_0110);
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_13 + IT_0004 + IT_0015 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0004 + IT_0040 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0003*IT_0005*IT_0008;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0124 = IT_0022*IT_0123;
    const ccomplex_t IT_0125 = 2*IT_0124;
    const ccomplex_t IT_0126 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0127 = IT_0028*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0130 = IT_0022*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = IT_0128 + IT_0131;
    const ccomplex_t IT_0133 = IT_0125 + IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = IT_0122*IT_0135;
    const ccomplex_t IT_0137 = IT_0119*IT_0120*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B2)*e_em;
    const ccomplex_t IT_0139 = IT_0027*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = conj(N_W2)*e_em;
    const ccomplex_t IT_0142 = IT_0008*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + 3*IT_0143);
    const ccomplex_t IT_0145 = 0.166666666666667*IT_0144;
    const ccomplex_t IT_0146 = IT_0135*IT_0145;
    const ccomplex_t IT_0147 = IT_0120*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = m_u*IT_0148;
    const ccomplex_t IT_0150 = conj(N_B3)*e_em;
    const ccomplex_t IT_0151 = IT_0027*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = conj(N_W3)*e_em;
    const ccomplex_t IT_0154 = IT_0008*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0152 + 3*IT_0155);
    const ccomplex_t IT_0157 = 0.166666666666667*IT_0156;
    const ccomplex_t IT_0158 = IT_0067*IT_0157;
    const ccomplex_t IT_0159 = IT_0069*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = m_u*IT_0160;
    const ccomplex_t IT_0162 = conj(N_B4)*e_em;
    const ccomplex_t IT_0163 = IT_0027*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W4)*e_em;
    const ccomplex_t IT_0166 = IT_0008*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*(IT_0164 + 3*IT_0167);
    const ccomplex_t IT_0169 = 0.166666666666667*IT_0168;
    const ccomplex_t IT_0170 = IT_0087*IT_0169;
    const ccomplex_t IT_0171 = IT_0072*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = m_u*IT_0172;
    const ccomplex_t IT_0174 = IT_0018 + -IT_0042 + IT_0051 + -IT_0070 + 
      -IT_0089 + IT_0097 + IT_0109 + IT_0118 + -IT_0137 + IT_0149 + IT_0161 +
       IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm2)*IT_0008;
    const ccomplex_t IT_0176 = IT_0045*IT_0175;
    const ccomplex_t IT_0177 = IT_0049*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0002*IT_0008;
    const ccomplex_t IT_0180 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (-2)*IT_0181;
    const ccomplex_t IT_0183 = IT_0002*IT_0027;
    const ccomplex_t IT_0184 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0185 = IT_0183*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0188 = IT_0179*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = IT_0186 + IT_0189;
    const ccomplex_t IT_0191 = IT_0182 + IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = 0.5*IT_0192;
    const ccomplex_t IT_0194 = IT_0145*IT_0193;
    const ccomplex_t IT_0195 = IT_0119*IT_0120*IT_0194;
    const ccomplex_t IT_0196 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0197 = IT_0179*IT_0196;
    const ccomplex_t IT_0198 = (-2)*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0200 = IT_0183*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0203 = IT_0179*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = IT_0201 + IT_0204;
    const ccomplex_t IT_0206 = IT_0198 + IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = 0.5*IT_0207;
    const ccomplex_t IT_0209 = IT_0157*IT_0208;
    const ccomplex_t IT_0210 = IT_0052*IT_0069*IT_0209;
    const ccomplex_t IT_0211 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0212 = IT_0179*IT_0211;
    const ccomplex_t IT_0213 = (-2)*IT_0212;
    const ccomplex_t IT_0214 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0215 = IT_0183*IT_0214;
    const ccomplex_t IT_0216 = 1.4142135623731*IT_0215;
    const ccomplex_t IT_0217 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0218 = IT_0179*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = IT_0216 + IT_0219;
    const ccomplex_t IT_0221 = IT_0213 + IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0222;
    const ccomplex_t IT_0224 = IT_0169*IT_0223;
    const ccomplex_t IT_0225 = IT_0071*IT_0072*IT_0224;
    const ccomplex_t IT_0226 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0227 = IT_0179*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0230 = IT_0183*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = IT_0228 + IT_0231;
    const ccomplex_t IT_0233 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0234 = IT_0179*IT_0233;
    const ccomplex_t IT_0235 = (-2)*IT_0234;
    const ccomplex_t IT_0236 = IT_0232 + IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = IT_0021*IT_0238;
    const ccomplex_t IT_0240 = IT_0041*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = m_u*IT_0241;
    const ccomplex_t IT_0243 = IT_0105*IT_0238;
    const ccomplex_t IT_0244 = IT_0019*IT_0041*IT_0243;
    const ccomplex_t IT_0245 = IT_0008*IT_0011;
    const ccomplex_t IT_0246 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0247 = IT_0245*IT_0246;
    const ccomplex_t IT_0248 = IT_0006*IT_0008*IT_0011;
    const ccomplex_t IT_0249 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0250 = IT_0248*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*(IT_0247 + (-0.5)*IT_0251);
    const ccomplex_t IT_0253 = IT_0010*IT_0252;
    const ccomplex_t IT_0254 = IT_0016*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0003*IT_0008;
    const ccomplex_t IT_0257 = 0.5*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0259 = (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = IT_0257*IT_0259;
    const ccomplex_t IT_0261 = cpow((-2)*s_13 + IT_0004 + IT_0015 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0262 = IT_0260*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *mu_h*IT_0003*IT_0008;
    const ccomplex_t IT_0265 = 0.5*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0267 = (-0.5)*IT_0266;
    const ccomplex_t IT_0268 = IT_0265*IT_0267;
    const ccomplex_t IT_0269 = cpow((-2)*s_13 + IT_0004 + IT_0015 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0270 = IT_0268*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0008;
    const ccomplex_t IT_0273 = IT_0112*IT_0272;
    const ccomplex_t IT_0274 = IT_0116*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0277 = IT_0248*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0280 = IT_0245*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*(IT_0278 + (-2)*IT_0280);
    const ccomplex_t IT_0282 = (-0.5)*IT_0281;
    const ccomplex_t IT_0283 = IT_0092*IT_0282;
    const ccomplex_t IT_0284 = IT_0090*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = IT_0122*IT_0193;
    const ccomplex_t IT_0287 = IT_0120*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*IT_0287;
    const ccomplex_t IT_0289 = m_u*IT_0288;
    const ccomplex_t IT_0290 = IT_0054*IT_0208;
    const ccomplex_t IT_0291 = IT_0069*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = m_u*IT_0292;
    const ccomplex_t IT_0294 = IT_0074*IT_0223;
    const ccomplex_t IT_0295 = IT_0072*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = m_u*IT_0296;
    const ccomplex_t IT_0298 = -IT_0178 + -IT_0195 + -IT_0210 + -IT_0225 +
       IT_0242 + -IT_0244 + -IT_0255 + -IT_0263 + -IT_0271 + -IT_0275 + -IT_0285
       + IT_0289 + IT_0293 + IT_0297;
    const ccomplex_t IT_0299 = 6*s_13;
    const ccomplex_t IT_0300 = m_u*m_C_2;
    const ccomplex_t IT_0301 = 6*IT_0174;
    const ccomplex_t IT_0302 = -IT_0108 + -IT_0148 + -IT_0160 + -IT_0172;
    const ccomplex_t IT_0303 = m_u*s_12;
    const ccomplex_t IT_0304 = -IT_0241 + -IT_0288 + -IT_0292 + -IT_0296;
    const ccomplex_t IT_0305 = s_23*m_C_2;
    const ccomplex_t IT_0306 = 6*conj(IT_0174);
    const ccomplex_t IT_0307 = 6*conj(IT_0298);
    const ccomplex_t IT_0308 = m_u*m_C_2*IT_0040;
    const ccomplex_t IT_0309 = s_13*IT_0040;
    const ccomplex_t IT_0310 = (-6)*IT_0309;
    const ccomplex_t IT_0311 = s_12*s_23;
    const ccomplex_t IT_0312 = 12*IT_0311;
    const ccomplex_t IT_0313 = IT_0310 + IT_0312;
    const ccomplex_t IT_0314 = 6*IT_0298;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0303;
    const ccomplex_t IT_0316 = 0.166666666666667*IT_0305;
    const ccomplex_t IT_0317 = (IT_0174*conj(IT_0174) + IT_0298*conj(IT_0298))
      *IT_0299 + conj(IT_0298)*IT_0300*IT_0301 + (IT_0298*IT_0300 + IT_0302
      *IT_0303 + IT_0304*IT_0305)*IT_0306 + (IT_0303*IT_0304 + IT_0302*IT_0305)
      *IT_0307 + 6*conj(IT_0304)*(IT_0302*IT_0308 + 0.166666666666667*IT_0304
      *IT_0313 + IT_0314*IT_0315 + IT_0301*IT_0316) + 6*conj(IT_0302)*(IT_0304
      *IT_0308 + 0.166666666666667*IT_0302*IT_0313 + IT_0301*IT_0315 + IT_0314
      *IT_0316);
    return create_ccomplex_return(IT_0317);
}

