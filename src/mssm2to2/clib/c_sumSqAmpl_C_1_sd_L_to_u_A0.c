#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sd_L_to_u_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sd_L_to_u_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = V_u1*U_Wm2;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = U_d2*V_Wp1;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = IT_0001*(IT_0005 + IT_0008);
    const ccomplex_t IT_0010 = 1.4142135623731*e_em*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0011;
    const ccomplex_t IT_0013 = pow(m_u, 2);
    const ccomplex_t IT_0014 = pow(m_sd_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + IT_0013 + -IT_0014 + -IT_0018 
      + -reg_prop, -1);
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = cpow(IT_0006, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0001*IT_0020*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = m_u*e_em*IT_0001*IT_0003*IT_0020*IT_0021;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0030 = V_u1*U_Wm1;
    const ccomplex_t IT_0031 = IT_0003*IT_0030;
    const ccomplex_t IT_0032 = U_d1*V_Wp1;
    const ccomplex_t IT_0033 = IT_0006*IT_0032;
    const ccomplex_t IT_0034 = IT_0001*(IT_0031 + IT_0033);
    const ccomplex_t IT_0035 = 1.4142135623731*e_em*IT_0034;
    const ccomplex_t IT_0036 = 0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0029*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0017 + -IT_0028 + -IT_0040;
    const ccomplex_t IT_0042 = pow(m_A0, 2);
    const ccomplex_t IT_0043 = m_u*m_C_1*IT_0042;
    const ccomplex_t IT_0044 = IT_0025*IT_0029;
    const ccomplex_t IT_0045 = IT_0019*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0001*IT_0020*IT_0021;
    const ccomplex_t IT_0048 = 0.5*IT_0047;
    const ccomplex_t IT_0049 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0050 = IT_0003*IT_0049;
    const ccomplex_t IT_0051 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0052 = IT_0006*IT_0051;
    const ccomplex_t IT_0053 = IT_0001*(IT_0050 + IT_0052);
    const ccomplex_t IT_0054 = 1.4142135623731*e_em*IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0048*IT_0055;
    const ccomplex_t IT_0057 = IT_0015*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0062 = IT_0006*IT_0061;
    const ccomplex_t IT_0063 = IT_0001*(IT_0060 + IT_0062);
    const ccomplex_t IT_0064 = 1.4142135623731*e_em*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0023*IT_0065;
    const ccomplex_t IT_0067 = IT_0038*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0046 + IT_0058 + IT_0068;
    const ccomplex_t IT_0070 = s_13*IT_0042;
    const ccomplex_t IT_0071 = s_14*s_34;
    const ccomplex_t IT_0072 = (-2)*IT_0071;
    const ccomplex_t IT_0073 = IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = (-6)*IT_0073;
    const ccomplex_t IT_0075 = 0.166666666666667*conj(IT_0041);
    const ccomplex_t IT_0076 = m_u*s_14;
    const ccomplex_t IT_0077 = 6*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0079 = IT_0023*IT_0036;
    const ccomplex_t IT_0080 = IT_0038*IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0082 = IT_0011*IT_0048;
    const ccomplex_t IT_0083 = IT_0015*IT_0081*IT_0082;
    const ccomplex_t IT_0084 = cpow((-2)*s_13 + IT_0013 + IT_0018 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0085 = 1.4142135623731*A_b*V_tb*V_td*IT_0006*conj
      (U_sb_10);
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0001*IT_0020*IT_0021*IT_0087*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = IT_0086*IT_0089;
    const ccomplex_t IT_0091 = IT_0084*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = m_C_1*IT_0068;
    const ccomplex_t IT_0094 = 1.4142135623731*A_b*V_tb*V_td*IT_0006*conj
      (U_sb_11);
    const ccomplex_t IT_0095 = 0.5*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0001*IT_0020*IT_0021*IT_0087*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0097 = 0.5*IT_0096;
    const ccomplex_t IT_0098 = IT_0095*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_13 + IT_0013 + IT_0018 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_C_1*IT_0058;
    const ccomplex_t IT_0103 = -IT_0080 + -IT_0083 + IT_0092 + -IT_0093 +
       IT_0101 + -IT_0102;
    const ccomplex_t IT_0104 = 0.166666666666667*conj(IT_0103);
    const ccomplex_t IT_0105 = s_34*m_C_1;
    const ccomplex_t IT_0106 = 6*IT_0105;
    const ccomplex_t IT_0107 = IT_0001*IT_0087;
    const ccomplex_t IT_0108 = e_em*conj(U_Wm1)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = cpow(IT_0003, -1);
    const ccomplex_t IT_0111 = IT_0001*IT_0087*IT_0110;
    const ccomplex_t IT_0112 = m_b*conj(U_d1)*e_em*IT_0020*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0113 = IT_0111*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0109 + (-0.5)*IT_0114);
    const ccomplex_t IT_0116 = IT_0086*IT_0115;
    const ccomplex_t IT_0117 = IT_0084*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0029*IT_0065;
    const ccomplex_t IT_0120 = IT_0038*IT_0078*IT_0119;
    const ccomplex_t IT_0121 = IT_0002*IT_0055;
    const ccomplex_t IT_0122 = IT_0015*IT_0081*IT_0121;
    const ccomplex_t IT_0123 = m_C_1*IT_0040;
    const ccomplex_t IT_0124 = 1.4142135623731*A_b*V_td*conj(V_ts)*IT_0006;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0001*IT_0020*IT_0110;
    const ccomplex_t IT_0127 = (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = IT_0125*IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_13 + IT_0013 + IT_0018 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = e_em*conj(U_Wm1)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0133 = IT_0107*IT_0132;
    const ccomplex_t IT_0134 = m_b*conj(U_d1)*e_em*IT_0020*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0135 = IT_0111*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + (-0.5)*IT_0136);
    const ccomplex_t IT_0138 = IT_0095*IT_0137;
    const ccomplex_t IT_0139 = IT_0099*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_d*e_em*mu_h*IT_0020;
    const ccomplex_t IT_0142 = IT_0001*IT_0141;
    const ccomplex_t IT_0143 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0144 = IT_0006*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = 0.5*IT_0142 + 0.5*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0001*IT_0020*IT_0110;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0146*IT_0148;
    const ccomplex_t IT_0150 = cpow((-2)*s_13 + IT_0013 + IT_0018 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_C_1*IT_0017;
    const ccomplex_t IT_0154 = -IT_0118 + IT_0120 + IT_0122 + IT_0123 + 
      -IT_0131 + -IT_0140 + -IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = 0.166666666666667*conj(IT_0154);
    const ccomplex_t IT_0156 = m_u*m_C_1;
    const ccomplex_t IT_0157 = 6*IT_0156;
    const ccomplex_t IT_0158 = 0.166666666666667*conj(IT_0069);
    const ccomplex_t IT_0159 = 6*IT_0041*(IT_0043*conj(IT_0069) + IT_0074
      *IT_0075 + IT_0077*IT_0104 + IT_0106*IT_0155) + 6*IT_0069*(conj(IT_0041)
      *IT_0043 + IT_0104*IT_0106 + IT_0077*IT_0155 + IT_0074*IT_0158) + 6
      *IT_0154*(IT_0075*IT_0106 + s_13*conj(IT_0154) + IT_0104*IT_0157 + IT_0077
      *IT_0158) + 6*IT_0103*(IT_0075*IT_0077 + s_13*conj(IT_0103) + IT_0155
      *IT_0157 + IT_0106*IT_0158);
    return create_ccomplex_return(IT_0159);
}

