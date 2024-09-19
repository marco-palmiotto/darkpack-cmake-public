#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sd_R_to_u_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sd_R_to_u_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*A_b*V_td*conj
      (V_ts)*IT_0000;
    const ccomplex_t IT_0002 = (-0.5)*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = pow(m_u, 2);
    const ccomplex_t IT_0012 = pow(m_C_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_00);
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0005*IT_0007*IT_0018*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*V_td
      *IT_0000*conj(U_sb_01);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*e_em
      *IT_0003*IT_0005*IT_0007*IT_0018*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cos(beta);
    const ccomplex_t IT_0034 = cpow(IT_0033, -1);
    const ccomplex_t IT_0035 = sin(alpha);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0000 + -m_d*e_em*mu_h*IT_0003*IT_0007*IT_0034*IT_0035);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_ud
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0015 + IT_0024 + IT_0032 + IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0046 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0047 = IT_0000*IT_0046;
    const ccomplex_t IT_0048 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0049 = IT_0035*IT_0048;
    const ccomplex_t IT_0050 = IT_0007*(IT_0047 + IT_0049);
    const ccomplex_t IT_0051 = 1.4142135623731*e_em*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0003*IT_0007*IT_0034;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = pow(m_sd_R, 2);
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0011 + IT_0057 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0045*IT_0056*IT_0058;
    const ccomplex_t IT_0060 = cpow((-2)*s_12 + IT_0011 + -IT_0012 + -IT_0057 
      + -reg_prop, -1);
    const ccomplex_t IT_0061 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0003*IT_0007*IT_0034;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*m_u*e_em*IT_0003*IT_0005
      *IT_0007*IT_0035;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = IT_0060*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = m_u*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*V_ud*e_em*conj(U_Wm2)*IT_0007;
    const ccomplex_t IT_0070 = IT_0037*IT_0069;
    const ccomplex_t IT_0071 = IT_0041*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0011 + IT_0057 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0075 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0076 = IT_0000*IT_0075;
    const ccomplex_t IT_0077 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0078 = IT_0035*IT_0077;
    const ccomplex_t IT_0079 = IT_0007*(IT_0076 + IT_0078);
    const ccomplex_t IT_0080 = 1.4142135623731*e_em*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = IT_0073*IT_0074*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm2)*IT_0007;
    const ccomplex_t IT_0086 = IT_0002*IT_0085;
    const ccomplex_t IT_0087 = IT_0013*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_13 + IT_0011 + IT_0012 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0090 = cos(theta_W);
    const ccomplex_t IT_0091 = cpow(IT_0090, -2);
    const ccomplex_t IT_0092 = pow(m_d, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*e_em*(m_W*IT_0006*(IT_0000
      *IT_0033 + -IT_0004*IT_0035)*IT_0091 + (-3)*IT_0000*IT_0003*IT_0007
      *IT_0034*IT_0092);
    const ccomplex_t IT_0094 = 0.333333333333333*IT_0093;
    const ccomplex_t IT_0095 = IT_0062*IT_0094;
    const ccomplex_t IT_0096 = IT_0089*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0007*IT_0018;
    const ccomplex_t IT_0099 = e_em*conj(U_Wm2)*U_sb_01*V_ub_mod;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = IT_0007*IT_0018*IT_0034;
    const ccomplex_t IT_0102 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_11*V_ub_mod;
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0100 + (-0.5)*IT_0104);
    const ccomplex_t IT_0106 = IT_0026*IT_0105;
    const ccomplex_t IT_0107 = IT_0030*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = e_em*conj(U_Wm2)*U_sb_00*V_ub_mod;
    const ccomplex_t IT_0110 = IT_0098*IT_0109;
    const ccomplex_t IT_0111 = m_b*conj(U_d2)*e_em*IT_0003*U_sb_10*V_ub_mod;
    const ccomplex_t IT_0112 = IT_0101*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + (-0.5)*IT_0113);
    const ccomplex_t IT_0115 = IT_0017*IT_0114;
    const ccomplex_t IT_0116 = IT_0022*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = U_d1*V_Wp2;
    const ccomplex_t IT_0119 = IT_0000*IT_0118;
    const ccomplex_t IT_0120 = V_u2*U_Wm1;
    const ccomplex_t IT_0121 = IT_0035*IT_0120;
    const ccomplex_t IT_0122 = IT_0007*(IT_0119 + IT_0121);
    const ccomplex_t IT_0123 = 1.4142135623731*e_em*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = 0.5*IT_0124;
    const ccomplex_t IT_0126 = IT_0055*IT_0125;
    const ccomplex_t IT_0127 = IT_0058*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = m_C_2*IT_0128;
    const ccomplex_t IT_0130 = U_d2*V_Wp2;
    const ccomplex_t IT_0131 = IT_0000*IT_0130;
    const ccomplex_t IT_0132 = V_u2*U_Wm2;
    const ccomplex_t IT_0133 = IT_0035*IT_0132;
    const ccomplex_t IT_0134 = IT_0007*(IT_0131 + IT_0133);
    const ccomplex_t IT_0135 = 1.4142135623731*e_em*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0062*IT_0137;
    const ccomplex_t IT_0139 = IT_0074*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_C_2*IT_0140;
    const ccomplex_t IT_0142 = IT_0059 + 2*IT_0068 + -IT_0072 + IT_0084 + 
      -IT_0088 + -IT_0097 + -IT_0108 + -IT_0117 + IT_0129 + IT_0141;
    const ccomplex_t IT_0143 = 6*s_13;
    const ccomplex_t IT_0144 = m_u*m_C_2;
    const ccomplex_t IT_0145 = 6*IT_0044;
    const ccomplex_t IT_0146 = -IT_0128 + -IT_0140;
    const ccomplex_t IT_0147 = m_u*s_14;
    const ccomplex_t IT_0148 = s_34*m_C_2;
    const ccomplex_t IT_0149 = 6*conj(IT_0044);
    const ccomplex_t IT_0150 = 6*conj(IT_0142);
    const ccomplex_t IT_0151 = pow(m_H0, 2);
    const ccomplex_t IT_0152 = m_u*m_C_2*IT_0151;
    const ccomplex_t IT_0153 = s_13*IT_0151;
    const ccomplex_t IT_0154 = s_14*s_34;
    const ccomplex_t IT_0155 = (-2)*IT_0154;
    const ccomplex_t IT_0156 = IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = 6*IT_0142;
    const ccomplex_t IT_0158 = 0.166666666666667*IT_0148;
    const ccomplex_t IT_0159 = (IT_0044*conj(IT_0044) + IT_0142*conj(IT_0142))
      *IT_0143 + conj(IT_0142)*IT_0144*IT_0145 + (IT_0142*IT_0144 + IT_0146
      *IT_0147 + IT_0067*IT_0148)*IT_0149 + (IT_0067*IT_0147 + IT_0146*IT_0148)
      *IT_0150 + 6*conj(IT_0067)*(IT_0146*IT_0152 + -IT_0067*IT_0156 +
       0.166666666666667*IT_0147*IT_0157 + IT_0145*IT_0158) + 6*conj(IT_0146)*
      (0.166666666666667*IT_0145*IT_0147 + IT_0067*IT_0152 + -IT_0146*IT_0156 +
       IT_0157*IT_0158);
    return create_ccomplex_return(IT_0159);
}

