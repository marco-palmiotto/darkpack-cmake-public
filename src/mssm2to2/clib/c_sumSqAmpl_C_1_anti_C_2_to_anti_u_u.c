#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_u_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_u_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = e_em*IT_0003;
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + (-3)*IT_0007);
    const ccomplex_t IT_0009 = (-0.166666666666667)*IT_0008;
    const ccomplex_t IT_0010 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0010;
    const ccomplex_t IT_0013 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + -IT_0012 + (-2)
      *IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = pow(m_C_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0018 + -IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0025 = (-0.666666666666667)*IT_0024;
    const ccomplex_t IT_0026 = IT_0016*IT_0025;
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0031 = cos(alpha);
    const ccomplex_t IT_0032 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = sin(alpha);
    const ccomplex_t IT_0035 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0005*(IT_0033 + IT_0036);
    const ccomplex_t IT_0038 = 1.4142135623731*e_em*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = pow(m_W, -1);
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*m_u*e_em*IT_0005*IT_0034
      *IT_0041*IT_0043;
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0040*IT_0045;
    const ccomplex_t IT_0047 = IT_0030*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0018 + -IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0031*IT_0035;
    const ccomplex_t IT_0051 = IT_0032*IT_0034;
    const ccomplex_t IT_0052 = IT_0005*(IT_0050 + -IT_0051);
    const ccomplex_t IT_0053 = 1.4142135623731*e_em*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*m_u*e_em*IT_0005*IT_0031
      *IT_0041*IT_0043;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = IT_0049*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = -IT_0048 + -IT_0060;
    const ccomplex_t IT_0062 = m_u*m_C_1;
    const ccomplex_t IT_0063 = -s_24;
    const ccomplex_t IT_0064 = s_23 + IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = 4*IT_0065;
    const ccomplex_t IT_0067 = 3*IT_0066;
    const ccomplex_t IT_0068 = conj(IT_0061)*IT_0067;
    const ccomplex_t IT_0069 = U_d2*V_Wp1;
    const ccomplex_t IT_0070 = IT_0031*IT_0069;
    const ccomplex_t IT_0071 = V_u1*U_Wm2;
    const ccomplex_t IT_0072 = IT_0034*IT_0071;
    const ccomplex_t IT_0073 = IT_0005*(IT_0070 + IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0045*IT_0076;
    const ccomplex_t IT_0078 = IT_0030*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = IT_0031*IT_0071;
    const ccomplex_t IT_0081 = IT_0034*IT_0069;
    const ccomplex_t IT_0082 = IT_0005*(IT_0080 + -IT_0081);
    const ccomplex_t IT_0083 = 1.4142135623731*e_em*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0057*IT_0085;
    const ccomplex_t IT_0087 = IT_0049*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = -IT_0079 + -IT_0088;
    const ccomplex_t IT_0090 = m_u*m_C_2;
    const ccomplex_t IT_0091 = -s_14;
    const ccomplex_t IT_0092 = s_13 + IT_0091;
    const ccomplex_t IT_0093 = IT_0090*IT_0092;
    const ccomplex_t IT_0094 = (-4)*IT_0093;
    const ccomplex_t IT_0095 = 3*IT_0094;
    const ccomplex_t IT_0096 = conj(IT_0089)*IT_0095;
    const ccomplex_t IT_0097 = IT_0068 + IT_0096;
    const ccomplex_t IT_0098 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0099 = IT_0003*IT_0098;
    const ccomplex_t IT_0100 = IT_0006*IT_0098;
    const ccomplex_t IT_0101 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0102 = IT_0006*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + -IT_0100 + (-2)
      *IT_0102);
    const ccomplex_t IT_0104 = 0.5*IT_0103;
    const ccomplex_t IT_0105 = IT_0025*IT_0104;
    const ccomplex_t IT_0106 = IT_0020*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 3*conj(IT_0107);
    const ccomplex_t IT_0109 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0110 = pow(m_u, 2);
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_us
      *e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_us*e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = IT_0113*IT_0115;
    const ccomplex_t IT_0117 = IT_0111*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0121 = 0.5*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_ud*e_em*IT_0005*IT_0041*IT_0043;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0121*IT_0123;
    const ccomplex_t IT_0125 = IT_0119*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0128 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0005*IT_0041*IT_0043*IT_0128*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0005*IT_0041*IT_0043*IT_0131*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0130*IT_0133;
    const ccomplex_t IT_0135 = IT_0127*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0138 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *e_em*IT_0005*IT_0041*IT_0043*IT_0128*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0139 = 0.5*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0005*IT_0041*IT_0043*IT_0131*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0139*IT_0141;
    const ccomplex_t IT_0143 = IT_0137*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0118 + IT_0126 + IT_0136 + IT_0144;
    const ccomplex_t IT_0146 = s_23 + s_24;
    const ccomplex_t IT_0147 = IT_0062*IT_0146;
    const ccomplex_t IT_0148 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0149 = cos(beta);
    const ccomplex_t IT_0150 = IT_0071*IT_0149;
    const ccomplex_t IT_0151 = IT_0042*IT_0069;
    const ccomplex_t IT_0152 = IT_0005*(IT_0150 + IT_0151);
    const ccomplex_t IT_0153 = 1.4142135623731*e_em*IT_0152;
    const ccomplex_t IT_0154 = 0.5*IT_0153;
    const ccomplex_t IT_0155 = m_u*e_em*IT_0005*IT_0041*IT_0043*IT_0149;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0154*IT_0156;
    const ccomplex_t IT_0158 = IT_0148*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = pow(m_Z, 2);
    const ccomplex_t IT_0161 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + IT_0160 
      + -reg_prop, -1);
    const ccomplex_t IT_0162 = IT_0069*IT_0149;
    const ccomplex_t IT_0163 = IT_0042*IT_0071;
    const ccomplex_t IT_0164 = IT_0005*(IT_0162 + -IT_0163);
    const ccomplex_t IT_0165 = 1.4142135623731*e_em*IT_0164;
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = m_u*e_em*IT_0005*IT_0041;
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0166*IT_0168;
    const ccomplex_t IT_0170 = IT_0161*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = -IT_0159 + -IT_0171;
    const ccomplex_t IT_0173 = m_u*s_24*m_C_1;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*V_us*e_em*U_Wm2*IT_0005;
    const ccomplex_t IT_0175 = IT_0113*IT_0174;
    const ccomplex_t IT_0176 = IT_0111*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0005*IT_0128;
    const ccomplex_t IT_0179 = e_em*U_Wm2*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0180 = IT_0178*IT_0179;
    const ccomplex_t IT_0181 = cpow(IT_0149, -1);
    const ccomplex_t IT_0182 = IT_0005*IT_0128*IT_0181;
    const ccomplex_t IT_0183 = m_b*U_d2*e_em*IT_0041*conj(U_sb_10)*V_ub_mod;
    const ccomplex_t IT_0184 = IT_0182*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0180 + (-0.5)*IT_0185);
    const ccomplex_t IT_0187 = IT_0133*IT_0186;
    const ccomplex_t IT_0188 = IT_0127*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*V_ud*e_em*U_Wm2*IT_0005;
    const ccomplex_t IT_0191 = IT_0121*IT_0190;
    const ccomplex_t IT_0192 = IT_0119*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = e_em*U_Wm2*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0195 = IT_0178*IT_0194;
    const ccomplex_t IT_0196 = m_b*U_d2*e_em*IT_0041*conj(U_sb_11)*V_ub_mod;
    const ccomplex_t IT_0197 = IT_0182*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + (-0.5)*IT_0198);
    const ccomplex_t IT_0200 = IT_0141*IT_0199;
    const ccomplex_t IT_0201 = IT_0137*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = -IT_0177 + -IT_0189 + -IT_0193 + -IT_0202;
    const ccomplex_t IT_0204 = s_12*IT_0110;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0005*IT_0041*IT_0181;
    const ccomplex_t IT_0206 = (-0.5)*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1.4142135623731)*m_d*U_d2*V_ud
      *e_em*IT_0005*IT_0041*IT_0181;
    const ccomplex_t IT_0208 = (-0.5)*IT_0207;
    const ccomplex_t IT_0209 = IT_0206*IT_0208;
    const ccomplex_t IT_0210 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0211 = IT_0209*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0005*IT_0041*IT_0181;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*V_us
      *e_em*IT_0005*IT_0041*IT_0181;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0214*IT_0216;
    const ccomplex_t IT_0218 = cpow((-2)*s_23 + IT_0019 + IT_0110 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0219 = IT_0217*IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0222 = IT_0174*IT_0221;
    const ccomplex_t IT_0223 = IT_0111*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0005*IT_0131;
    const ccomplex_t IT_0226 = e_em*conj(U_Wm1)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0227 = IT_0225*IT_0226;
    const ccomplex_t IT_0228 = IT_0005*IT_0131*IT_0181;
    const ccomplex_t IT_0229 = m_b*conj(U_d1)*e_em*IT_0041*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0227 + (-0.5)*IT_0231);
    const ccomplex_t IT_0233 = IT_0186*IT_0232;
    const ccomplex_t IT_0234 = IT_0127*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0005;
    const ccomplex_t IT_0237 = IT_0190*IT_0236;
    const ccomplex_t IT_0238 = IT_0119*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = e_em*conj(U_Wm1)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0241 = IT_0225*IT_0240;
    const ccomplex_t IT_0242 = m_b*conj(U_d1)*e_em*IT_0041*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0243 = IT_0228*IT_0242;
    const ccomplex_t IT_0244 = 1.4142135623731*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*(IT_0241 + (-0.5)*IT_0244);
    const ccomplex_t IT_0246 = IT_0199*IT_0245;
    const ccomplex_t IT_0247 = IT_0137*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0212 + IT_0220 + IT_0224 + IT_0235 + IT_0239
       + IT_0248;
    const ccomplex_t IT_0250 = m_u*s_13*m_C_2;
    const ccomplex_t IT_0251 = IT_0123*IT_0236;
    const ccomplex_t IT_0252 = IT_0119*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = IT_0139*IT_0245;
    const ccomplex_t IT_0255 = IT_0137*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = IT_0115*IT_0221;
    const ccomplex_t IT_0258 = IT_0111*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = IT_0130*IT_0232;
    const ccomplex_t IT_0261 = IT_0127*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = -IT_0253 + -IT_0256 + -IT_0259 + -IT_0262;
    const ccomplex_t IT_0264 = 12*conj(IT_0107);
    const ccomplex_t IT_0265 = m_C_1*m_C_2*IT_0110;
    const ccomplex_t IT_0266 = m_u*s_23*m_C_1;
    const ccomplex_t IT_0267 = s_14*s_23;
    const ccomplex_t IT_0268 = m_u*s_14*m_C_2;
    const ccomplex_t IT_0269 = IT_0009*IT_0104;
    const ccomplex_t IT_0270 = IT_0020*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = 24*conj(IT_0271);
    const ccomplex_t IT_0273 = IT_0035*IT_0149;
    const ccomplex_t IT_0274 = IT_0032*IT_0042;
    const ccomplex_t IT_0275 = IT_0005*(IT_0273 + IT_0274);
    const ccomplex_t IT_0276 = 1.4142135623731*e_em*IT_0275;
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = IT_0156*IT_0277;
    const ccomplex_t IT_0279 = IT_0148*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = IT_0032*IT_0149;
    const ccomplex_t IT_0282 = IT_0035*IT_0042;
    const ccomplex_t IT_0283 = IT_0005*(IT_0281 + -IT_0282);
    const ccomplex_t IT_0284 = 1.4142135623731*e_em*IT_0283;
    const ccomplex_t IT_0285 = 0.5*IT_0284;
    const ccomplex_t IT_0286 = IT_0168*IT_0285;
    const ccomplex_t IT_0287 = IT_0161*IT_0286;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*IT_0287;
    const ccomplex_t IT_0289 = -IT_0280 + -IT_0288;
    const ccomplex_t IT_0290 = s_13 + s_14;
    const ccomplex_t IT_0291 = IT_0090*IT_0290;
    const ccomplex_t IT_0292 = s_12*s_34;
    const ccomplex_t IT_0293 = s_13*s_24;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = IT_0204 + IT_0267 + IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = m_C_1*m_C_2;
    const ccomplex_t IT_0297 = s_34 + IT_0110;
    const ccomplex_t IT_0298 = IT_0296*IT_0297;
    const ccomplex_t IT_0299 = 6*IT_0298;
    const ccomplex_t IT_0300 = 12*IT_0147;
    const ccomplex_t IT_0301 = s_12*IT_0297;
    const ccomplex_t IT_0302 = 24*IT_0301;
    const ccomplex_t IT_0303 = (-24)*IT_0298;
    const ccomplex_t IT_0304 = (-12)*IT_0147;
    const ccomplex_t IT_0305 = IT_0145*IT_0147 + 0.166666666666667*((-12)
      *IT_0107 + 6*IT_0249 + 12*IT_0271)*IT_0291 + IT_0263*IT_0295 +
       0.166666666666667*IT_0203*IT_0299 + 0.166666666666667*IT_0029*IT_0300 +
       0.166666666666667*IT_0289*IT_0302 + 0.166666666666667*IT_0172*IT_0303 +
       0.166666666666667*IT_0023*IT_0304;
    const ccomplex_t IT_0306 = -IT_0110;
    const ccomplex_t IT_0307 = s_34 + IT_0306;
    const ccomplex_t IT_0308 = IT_0296*IT_0307;
    const ccomplex_t IT_0309 = (-2)*IT_0308;
    const ccomplex_t IT_0310 = 3*IT_0309;
    const ccomplex_t IT_0311 = 12*IT_0267;
    const ccomplex_t IT_0312 = 12*IT_0266;
    const ccomplex_t IT_0313 = 12*IT_0268;
    const ccomplex_t IT_0314 = 12*IT_0265;
    const ccomplex_t IT_0315 = 12*IT_0250;
    const ccomplex_t IT_0316 = 12*IT_0107;
    const ccomplex_t IT_0317 = 24*IT_0271;
    const ccomplex_t IT_0318 = 24*IT_0268;
    const ccomplex_t IT_0319 = -IT_0292;
    const ccomplex_t IT_0320 = -IT_0267;
    const ccomplex_t IT_0321 = IT_0204 + IT_0293 + IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = (-6)*IT_0321;
    const ccomplex_t IT_0323 = IT_0172*IT_0295 + (-0.166666666666667)*IT_0289
      *IT_0299 + (-0.166666666666667)*IT_0061*IT_0310 + (-0.166666666666667)
      *IT_0203*IT_0311 + (-0.166666666666667)*IT_0249*IT_0312 + (
      -0.166666666666667)*IT_0145*IT_0313 + (-0.166666666666667)*IT_0263*IT_0314
       + (-0.166666666666667)*IT_0023*IT_0315 + (-0.166666666666667)*IT_0173
      *IT_0316 + (-0.166666666666667)*IT_0266*IT_0317 + (-0.166666666666667)
      *IT_0029*IT_0318 + (-0.166666666666667)*IT_0089*IT_0322;
    const ccomplex_t IT_0324 = 12*IT_0291;
    const ccomplex_t IT_0325 = (-12)*IT_0291;
    const ccomplex_t IT_0326 = (-6)*IT_0147;
    const ccomplex_t IT_0327 = (-6)*IT_0291;
    const ccomplex_t IT_0328 = (-6)*IT_0298;
    const ccomplex_t IT_0329 = IT_0203*IT_0295 + (-0.166666666666667)*IT_0172
      *IT_0302 + (-0.166666666666667)*IT_0289*IT_0303 + (-0.166666666666667)
      *IT_0147*((-12)*IT_0271 + IT_0316) + (-0.166666666666667)*IT_0023*IT_0324 
      + (-0.166666666666667)*IT_0029*IT_0325 + (-0.166666666666667)*IT_0249
      *IT_0326 + (-0.166666666666667)*IT_0145*IT_0327 + (-0.166666666666667)
      *IT_0263*IT_0328;
    const ccomplex_t IT_0330 = 48*IT_0107;
    const ccomplex_t IT_0331 = 24*IT_0109;
    const ccomplex_t IT_0332 = 48*IT_0265;
    const ccomplex_t IT_0333 = 0.166666666666667*IT_0145;
    const ccomplex_t IT_0334 = (-2)*IT_0093;
    const ccomplex_t IT_0335 = 3*IT_0334;
    const ccomplex_t IT_0336 = 0.166666666666667*IT_0089;
    const ccomplex_t IT_0337 = 2*IT_0065;
    const ccomplex_t IT_0338 = 3*IT_0337;
    const ccomplex_t IT_0339 = 0.166666666666667*IT_0061;
    const ccomplex_t IT_0340 = 0.166666666666667*IT_0172;
    const ccomplex_t IT_0341 = 0.166666666666667*IT_0203;
    const ccomplex_t IT_0342 = 0.166666666666667*IT_0249;
    const ccomplex_t IT_0343 = 0.166666666666667*IT_0263;
    const ccomplex_t IT_0344 = 12*IT_0204;
    const ccomplex_t IT_0345 = 0.166666666666667*IT_0023;
    const ccomplex_t IT_0346 = 24*IT_0267;
    const ccomplex_t IT_0347 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0348 = IT_0147*IT_0289 + 0.166666666666667*IT_0109
      *IT_0316 + 0.166666666666667*IT_0265*IT_0317 + IT_0311*IT_0333 + IT_0335
      *IT_0336 + IT_0338*IT_0339 + IT_0327*IT_0340 + IT_0313*IT_0341 + IT_0314
      *IT_0342 + IT_0312*IT_0343 + IT_0344*IT_0345 + IT_0346*IT_0347;
    const ccomplex_t IT_0349 = 12*IT_0173;
    const ccomplex_t IT_0350 = 24*IT_0266;
    const ccomplex_t IT_0351 = IT_0289*IT_0295 + 0.166666666666667*IT_0250
      *IT_0316 + 0.166666666666667*IT_0268*IT_0317 + IT_0312*IT_0333 + IT_0310
      *IT_0336 + IT_0322*IT_0339 + IT_0328*IT_0340 + IT_0314*IT_0341 + IT_0313
      *IT_0342 + IT_0311*IT_0343 + IT_0345*IT_0349 + IT_0347*IT_0350;
    const ccomplex_t IT_0352 = 12*IT_0109;
    const ccomplex_t IT_0353 = 24*IT_0265;
    const ccomplex_t IT_0354 = IT_0289*IT_0291 + 0.166666666666667*IT_0204
      *IT_0316 + 0.166666666666667*IT_0267*IT_0317 + IT_0314*IT_0333 + IT_0336
      *IT_0338 + IT_0335*IT_0339 + IT_0326*IT_0340 + IT_0312*IT_0341 + IT_0311
      *IT_0342 + IT_0313*IT_0343 + IT_0345*IT_0352 + IT_0347*IT_0353;
    const ccomplex_t IT_0355 = -s_34;
    const ccomplex_t IT_0356 = IT_0110 + IT_0355;
    const ccomplex_t IT_0357 = s_12*IT_0356;
    const ccomplex_t IT_0358 = (-24)*IT_0357;
    const ccomplex_t IT_0359 = 0.333333333333333*IT_0089;
    const ccomplex_t IT_0360 = (-8)*IT_0308;
    const ccomplex_t IT_0361 = 3*IT_0360;
    const ccomplex_t IT_0362 = 0.333333333333333*IT_0061;
    const ccomplex_t IT_0363 = 0.333333333333333*IT_0203;
    const ccomplex_t IT_0364 = 0.333333333333333*IT_0249;
    const ccomplex_t IT_0365 = 0.333333333333333*IT_0263;
    const ccomplex_t IT_0366 = 3*IT_0107;
    const ccomplex_t IT_0367 = 0.333333333333333*IT_0366;
    const ccomplex_t IT_0368 = IT_0145*IT_0335 + 3*IT_0358*IT_0359 + 3*IT_0361
      *IT_0362 + 3*IT_0322*IT_0363 + 3*IT_0338*IT_0364 + 3*IT_0310*IT_0365 + 3
      *IT_0066*(IT_0271 + IT_0367);
    const ccomplex_t IT_0369 = IT_0145*IT_0338 + 3*IT_0359*IT_0361 + 3*IT_0358
      *IT_0362 + 3*IT_0310*IT_0363 + 3*IT_0335*IT_0364 + 3*IT_0322*IT_0365 + 3
      *IT_0094*(IT_0271 + IT_0367);
    const ccomplex_t IT_0370 = 24*IT_0204;
    const ccomplex_t IT_0371 = IT_0089*IT_0095;
    const ccomplex_t IT_0372 = IT_0061*IT_0067;
    const ccomplex_t IT_0373 = 48*IT_0023*IT_0293 + IT_0289*IT_0304 + IT_0203
      *IT_0315 + IT_0109*IT_0317 + IT_0172*IT_0324 + IT_0265*IT_0330 + IT_0145
      *IT_0344 + IT_0263*IT_0349 + IT_0249*IT_0352 + IT_0029*IT_0370 + IT_0371 +
       IT_0372;
    const ccomplex_t IT_0374 = 48*IT_0029*IT_0267 + 48*IT_0265*IT_0271 +
       IT_0289*IT_0300 + IT_0203*IT_0318 + IT_0172*IT_0325 + IT_0107*IT_0331 +
       IT_0145*IT_0346 + IT_0263*IT_0350 + IT_0249*IT_0353 + IT_0023*IT_0370 +
       IT_0371 + IT_0372;
    const ccomplex_t IT_0375 = (IT_0023 + IT_0029)*IT_0097 + (IT_0066*IT_0089 
      + IT_0061*IT_0094)*IT_0108 + (IT_0109*IT_0145 + IT_0147*IT_0172 + IT_0173
      *IT_0203 + IT_0204*IT_0249 + IT_0250*IT_0263)*IT_0264 + (IT_0023*IT_0109 +
       IT_0107*IT_0204 + IT_0145*IT_0265 + IT_0203*IT_0266 + IT_0249*IT_0267 +
       IT_0263*IT_0268)*IT_0272 + 6*conj(IT_0289)*IT_0305 + (-6)*conj(IT_0203)
      *IT_0323 + (-6)*conj(IT_0172)*IT_0329 + conj(IT_0107)*((-12)*IT_0289
      *IT_0291 + IT_0204*IT_0317 + IT_0293*IT_0330 + IT_0029*IT_0331 + IT_0023
      *IT_0332) + 3*conj(IT_0271)*(IT_0066*IT_0089 + IT_0061*IT_0094 + (-4)
      *IT_0147*IT_0172 + 16*IT_0267*IT_0271 + 4*IT_0289*IT_0291 + 1./3*IT_0029
      *IT_0332) + 6*conj(IT_0145)*IT_0348 + 6*conj(IT_0263)*IT_0351 + 6*conj
      (IT_0249)*IT_0354 + conj(IT_0089)*IT_0368 + conj(IT_0061)*IT_0369 + conj
      (IT_0023)*IT_0373 + conj(IT_0029)*IT_0374;
    return create_ccomplex_return(IT_0375);
}

