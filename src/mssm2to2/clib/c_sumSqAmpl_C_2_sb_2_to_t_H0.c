#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_t_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_t_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
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
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sb_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = U_d1*V_Wp2;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = V_u2*U_Wm1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (IT_0005 + IT_0008)*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*e_em*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = m_b*conj(U_d1)*V_tb*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = V_tb*e_em*conj(U_Wm1)*U_sb_01;
    const ccomplex_t IT_0023 = IT_0010*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + (-2)*IT_0023);
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0030 = U_d2*V_Wp2;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = V_u2*U_Wm2;
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = IT_0010*(IT_0031 + IT_0033);
    const ccomplex_t IT_0035 = 1.4142135623731*e_em*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = m_b*conj(U_d2)*V_tb*e_em*IT_0018*U_sb_11;
    const ccomplex_t IT_0039 = IT_0017*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0042 = IT_0010*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0040 + (-2)*IT_0042);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0037*IT_0044;
    const ccomplex_t IT_0046 = IT_0029*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = pow(m_C_2, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0001 + -IT_0048 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0050 = sin(beta);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*m_t*e_em*IT_0006*IT_0010
      *IT_0018*IT_0051;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0010*IT_0018*IT_0051*U_sb_01;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = IT_0049*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0028 + -IT_0047 + -IT_0058;
    const ccomplex_t IT_0060 = pow(m_H0, 2);
    const ccomplex_t IT_0061 = m_t*m_C_2*IT_0060;
    const ccomplex_t IT_0062 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0063 = IT_0003*IT_0062;
    const ccomplex_t IT_0064 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0065 = IT_0006*IT_0064;
    const ccomplex_t IT_0066 = IT_0010*(IT_0063 + IT_0065);
    const ccomplex_t IT_0067 = 1.4142135623731*e_em*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u1
      *e_em*IT_0010*IT_0018*IT_0051*U_sb_01;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = IT_0002*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0044*IT_0053;
    const ccomplex_t IT_0076 = IT_0049*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0079 = IT_0003*IT_0078;
    const ccomplex_t IT_0080 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0081 = IT_0006*IT_0080;
    const ccomplex_t IT_0082 = IT_0010*(IT_0079 + IT_0081);
    const ccomplex_t IT_0083 = 1.4142135623731*e_em*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = IT_0055*IT_0085;
    const ccomplex_t IT_0087 = IT_0029*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0074 + IT_0077 + IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0091 = IT_0037*IT_0055;
    const ccomplex_t IT_0092 = IT_0029*IT_0090*IT_0091;
    const ccomplex_t IT_0093 = m_t*IT_0058;
    const ccomplex_t IT_0094 = m_C_2*IT_0088;
    const ccomplex_t IT_0095 = cpow((-2)*s_13 + IT_0000 + IT_0048 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0096 = cpow(V_tb, 2);
    const ccomplex_t IT_0097 = m_b*mu_h*IT_0018;
    const ccomplex_t IT_0098 = -IT_0006*IT_0050;
    const ccomplex_t IT_0099 = IT_0003*IT_0015;
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = cos(theta_W);
    const ccomplex_t IT_0102 = cpow(IT_0101, -2);
    const ccomplex_t IT_0103 = pow(m_b, 2);
    const ccomplex_t IT_0104 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0100*(IT_0010 + 1.
      /3*IT_0009*IT_0102) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0103)*conj
      (U_sb_01) + -1./2*IT_0006*IT_0010*IT_0016*IT_0097*conj(U_sb_11)) + U_sb_11
      *(IT_0006*IT_0010*IT_0016*IT_0097*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0009*IT_0100*IT_0102 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0103)
      *conj(U_sb_11));
    const ccomplex_t IT_0105 = (0 + _Complex_I*-1)*(e_em*IT_0104 + (
      -1.4142135623731)*A_b*IT_0003*IT_0096*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0055*IT_0106;
    const ccomplex_t IT_0108 = IT_0095*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_13 + IT_0000 + IT_0048 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*m_t*V_ts*V_u2
      *e_em*IT_0010*IT_0018*IT_0051;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0003*U_sb_11;
    const ccomplex_t IT_0114 = (-0.5)*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = IT_0110*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = cpow((-2)*s_13 + IT_0000 + IT_0048 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_t*V_td*V_u2
      *e_em*IT_0010*IT_0018*IT_0051;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_11;
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0120*IT_0122;
    const ccomplex_t IT_0124 = IT_0118*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = cpow((-2)*s_13 + IT_0000 + IT_0048 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0010*IT_0018*IT_0051*U_sb_00;
    const ccomplex_t IT_0128 = 0.5*IT_0127;
    const ccomplex_t IT_0129 = (-2)*conj(U_sb_00)*(((-0.5)*m_W*IT_0100*
      (IT_0010 + 1./3*IT_0009*IT_0102) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0103
      )*U_sb_01 + -1./2*IT_0006*IT_0010*IT_0016*IT_0097*U_sb_11) + conj(U_sb_10)
      *(IT_0006*IT_0010*IT_0016*IT_0097*U_sb_01 + 0.666666666666667*(m_W*IT_0009
      *IT_0100*IT_0102 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0103)*U_sb_11);
    const ccomplex_t IT_0130 = (0 + _Complex_I*-1)*(e_em*IT_0129 + (
      -1.4142135623731)*A_b*IT_0003*IT_0096*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0128*IT_0131;
    const ccomplex_t IT_0133 = IT_0126*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0136 = IT_0014*IT_0071;
    const ccomplex_t IT_0137 = IT_0002*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = m_C_2*IT_0074;
    const ccomplex_t IT_0139 = -IT_0092 + (-2)*IT_0093 + -IT_0094 + IT_0109 +
       IT_0117 + IT_0125 + IT_0134 + -IT_0137 + -IT_0138;
    const ccomplex_t IT_0140 = m_t*s_14;
    const ccomplex_t IT_0141 = 6*IT_0140;
    const ccomplex_t IT_0142 = s_34*m_C_2;
    const ccomplex_t IT_0143 = 6*IT_0142;
    const ccomplex_t IT_0144 = cpow((-2)*s_13 + IT_0000 + IT_0048 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0145 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_ts*e_em*IT_0010*IT_0016*IT_0018;
    const ccomplex_t IT_0146 = (-0.5)*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_ts)*IT_0003*U_sb_01;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0146*IT_0148;
    const ccomplex_t IT_0150 = IT_0144*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_13 + IT_0000 + IT_0048 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0153 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_td*e_em*IT_0010*IT_0016*IT_0018;
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0003*U_sb_01;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0154*IT_0156;
    const ccomplex_t IT_0158 = IT_0152*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0100*(IT_0010 + 1.
      /3*IT_0009*IT_0102) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0103)*conj
      (U_sb_01) + -1./2*IT_0006*IT_0010*IT_0016*IT_0097*conj(U_sb_11)) + U_sb_11
      *(IT_0006*IT_0010*IT_0016*IT_0097*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0009*IT_0100*IT_0102 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0103)
      *conj(U_sb_11));
    const ccomplex_t IT_0161 = (0 + _Complex_I*-1)*(e_em*IT_0160 + (
      -1.4142135623731)*A_b*IT_0003*IT_0096*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0044*IT_0162;
    const ccomplex_t IT_0164 = IT_0095*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*V_ts*e_em*conj(U_Wm2)*IT_0010;
    const ccomplex_t IT_0167 = IT_0114*IT_0166;
    const ccomplex_t IT_0168 = IT_0110*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*V_td*e_em*conj(U_Wm2)*IT_0010;
    const ccomplex_t IT_0171 = IT_0122*IT_0170;
    const ccomplex_t IT_0172 = IT_0118*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = m_b*conj(U_d2)*V_tb*e_em*IT_0018*U_sb_10;
    const ccomplex_t IT_0175 = IT_0017*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0178 = IT_0010*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0176 + (-2)*IT_0178);
    const ccomplex_t IT_0180 = (-0.5)*IT_0179;
    const ccomplex_t IT_0181 = (-2)*conj(U_sb_00)*(((-0.5)*m_W*IT_0100*
      (IT_0010 + 1./3*IT_0009*IT_0102) + IT_0003*IT_0010*IT_0016*IT_0018*IT_0103
      )*U_sb_01 + -1./2*IT_0006*IT_0010*IT_0016*IT_0097*U_sb_11) + conj(U_sb_10)
      *(IT_0006*IT_0010*IT_0016*IT_0097*U_sb_01 + 0.666666666666667*(m_W*IT_0009
      *IT_0100*IT_0102 + (-3)*IT_0003*IT_0010*IT_0016*IT_0018*IT_0103)*U_sb_11);
    const ccomplex_t IT_0182 = (0 + _Complex_I*-1)*(e_em*IT_0181 + (
      -1.4142135623731)*A_b*IT_0003*IT_0096*(U_sb_01*conj(U_sb_10) + conj
      (U_sb_00)*U_sb_11));
    const ccomplex_t IT_0183 = (-0.5)*IT_0182;
    const ccomplex_t IT_0184 = IT_0180*IT_0183;
    const ccomplex_t IT_0185 = IT_0126*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = IT_0025*IT_0069;
    const ccomplex_t IT_0188 = IT_0002*IT_0135*IT_0187;
    const ccomplex_t IT_0189 = IT_0044*IT_0085;
    const ccomplex_t IT_0190 = IT_0029*IT_0090*IT_0189;
    const ccomplex_t IT_0191 = m_t*IT_0077;
    const ccomplex_t IT_0192 = m_C_2*IT_0028;
    const ccomplex_t IT_0193 = m_C_2*IT_0047;
    const ccomplex_t IT_0194 = -IT_0151 + -IT_0159 + -IT_0165 + -IT_0169 + 
      -IT_0173 + -IT_0186 + IT_0188 + IT_0190 + 2*IT_0191 + IT_0192 + IT_0193;
    const ccomplex_t IT_0195 = 0.166666666666667*conj(IT_0194);
    const ccomplex_t IT_0196 = s_13*IT_0060;
    const ccomplex_t IT_0197 = s_14*s_34;
    const ccomplex_t IT_0198 = (-2)*IT_0197;
    const ccomplex_t IT_0199 = IT_0196 + IT_0198;
    const ccomplex_t IT_0200 = (-6)*IT_0199;
    const ccomplex_t IT_0201 = 0.166666666666667*conj(IT_0059);
    const ccomplex_t IT_0202 = m_t*m_C_2;
    const ccomplex_t IT_0203 = 6*IT_0202;
    const ccomplex_t IT_0204 = 0.166666666666667*conj(IT_0089);
    const ccomplex_t IT_0205 = 6*IT_0059*(IT_0061*conj(IT_0089) +
       0.166666666666667*conj(IT_0139)*IT_0141 + IT_0143*IT_0195 + IT_0200
      *IT_0201) + 6*IT_0194*(s_13*conj(IT_0194) + IT_0143*IT_0201 +
       0.166666666666667*conj(IT_0139)*IT_0203 + IT_0141*IT_0204) + 6*IT_0139*
      (s_13*conj(IT_0139) + IT_0141*IT_0201 + IT_0195*IT_0203 + IT_0143*IT_0204)
       + 6*IT_0089*(conj(IT_0059)*IT_0061 + 0.166666666666667*conj(IT_0139)
      *IT_0143 + IT_0141*IT_0195 + IT_0200*IT_0204);
    return create_ccomplex_return(IT_0205);
}

