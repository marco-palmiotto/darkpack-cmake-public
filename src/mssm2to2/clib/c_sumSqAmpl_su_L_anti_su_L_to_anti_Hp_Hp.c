#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_L_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_L_to_anti_Hp_Hp(
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
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = pow(m_W, -2);
    const ccomplex_t IT_0010 = pow(m_d, 2);
    const ccomplex_t IT_0011 = pow(V_ud, 2);
    const ccomplex_t IT_0012 = pow(m_s, 2);
    const ccomplex_t IT_0013 = pow(V_us, 2);
    const ccomplex_t IT_0014 = pow(m_b, 2);
    const ccomplex_t IT_0015 = pow(V_ub_mod, 2);
    const ccomplex_t IT_0016 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0017 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0018 = cpow(IT_0001, -2);
    const ccomplex_t IT_0019 = IT_0002*(IT_0004 + 3*IT_0006) + -IT_0008*
      (IT_0004 + 3*IT_0006*(1 + (-2)*IT_0009*(IT_0010*IT_0011 + IT_0012*IT_0013 
      + IT_0014*IT_0015*IT_0016*IT_0017)*IT_0018));
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0000*IT_0019;
    const ccomplex_t IT_0021 = (-0.0833333333333333)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0023 = (-0.666666666666667)*IT_0022;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = pow(m_su_L, 2);
    const ccomplex_t IT_0026 = s_12 + IT_0025;
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = 2*IT_0022;
    const ccomplex_t IT_0029 = IT_0023*IT_0028;
    const ccomplex_t IT_0030 = s_13 + s_23;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (-1.33333333333333)*IT_0022;
    const ccomplex_t IT_0033 = IT_0022*IT_0032;
    const ccomplex_t IT_0034 = IT_0026*IT_0033;
    const ccomplex_t IT_0035 = IT_0028*IT_0032;
    const ccomplex_t IT_0036 = s_23*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0031 + (-0.5)*IT_0034 + 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0027 + IT_0037;
    const ccomplex_t IT_0039 = cpow(s_12 + IT_0025 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = cpow(IT_0003, -1);
    const ccomplex_t IT_0043 = cpow(IT_0005, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*(IT_0005*IT_0042 + 
      -IT_0003*IT_0043);
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*e_em*(IT_0005*IT_0042 + (-3)
      *IT_0003*IT_0043);
    const ccomplex_t IT_0046 = 0.166666666666667*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*s_13 + (-0.5)*s_23;
    const ccomplex_t IT_0049 = IT_0026 + IT_0048;
    const ccomplex_t IT_0050 = IT_0047*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0044;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = IT_0026*IT_0052;
    const ccomplex_t IT_0054 = 0.333333333333333*IT_0045;
    const ccomplex_t IT_0055 = IT_0051*IT_0054;
    const ccomplex_t IT_0056 = IT_0026*IT_0055;
    const ccomplex_t IT_0057 = IT_0044*IT_0054;
    const ccomplex_t IT_0058 = -s_23;
    const ccomplex_t IT_0059 = IT_0026 + IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = -IT_0053 + 0.5*IT_0056 + (-0.5)*IT_0060;
    const ccomplex_t IT_0062 = IT_0050 + IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0025 + -reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = pow(m_W, -1);
    const ccomplex_t IT_0067 = pow(m_u, 2);
    const ccomplex_t IT_0068 = cpow(IT_0007, -1);
    const ccomplex_t IT_0069 = sin(alpha);
    const ccomplex_t IT_0070 = cos(alpha);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*e_em*(IT_0043*IT_0066
      *IT_0067*IT_0068*IT_0069 + 0.166666666666667*m_W*(IT_0004*IT_0005 + (-3)
      *IT_0043)*(IT_0007*IT_0069 + -IT_0001*IT_0070));
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = -IT_0043;
    const ccomplex_t IT_0074 = IT_0004*IT_0005;
    const ccomplex_t IT_0075 = IT_0073 + IT_0074;
    const ccomplex_t IT_0076 = cpow(IT_0007, 3);
    const ccomplex_t IT_0077 = -IT_0001;
    const ccomplex_t IT_0078 = 3*IT_0043 + IT_0074;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_W*e_em*(IT_0069*IT_0075
      *IT_0076 + IT_0077*(IT_0008*IT_0070*IT_0078 + IT_0077*(IT_0001*IT_0070
      *IT_0075 + -IT_0007*IT_0069*IT_0078)));
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0072*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_12 + (-2)*IT_0025 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0043*IT_0066;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = cpow(IT_0086, 2);
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = pow(m_Hp, 2);
    const ccomplex_t IT_0090 = cpow((-2)*s_23 + IT_0025 + IT_0089 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *mu_h*IT_0043*IT_0066;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = cpow(IT_0093, 2);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_23 + IT_0025 + IT_0089 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = cpow(IT_0001, -1);
    const ccomplex_t IT_0099 = (0 + _Complex_I*(-2.82842712474619))*V_ud*e_em
      *IT_0043*(IT_0001*(m_W*IT_0007 + (-0.5)*IT_0066*IT_0067*IT_0068) + -1./2
      *IT_0007*IT_0010*IT_0066*IT_0098);
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = cpow(IT_0100, 2);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0025 + IT_0089 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*(-2.82842712474619))*V_us*e_em
      *IT_0043*(IT_0001*(m_W*IT_0007 + (-0.5)*IT_0066*IT_0067*IT_0068) + -1./2
      *IT_0007*IT_0012*IT_0066*IT_0098);
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = cpow(IT_0106, 2);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_23 + IT_0025 + IT_0089 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016
      *IT_0043*(m_W*IT_0001*IT_0007*conj(U_sb_01) + (-0.5)*IT_0066*(IT_0001
      *IT_0067*IT_0068*conj(U_sb_01) + m_b*(m_b*IT_0007*IT_0098*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0017
      *IT_0043*(m_W*IT_0001*IT_0007*U_sb_01 + (-0.5)*IT_0066*(IT_0001*IT_0067
      *IT_0068*U_sb_01 + m_b*(m_b*IT_0007*IT_0098*U_sb_01 + mu_h*U_sb_11)))
      *V_ub_mod;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = cpow((-2)*s_23 + IT_0025 + IT_0089 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0025 + -reg_prop, -1);
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*e_em*(IT_0043*IT_0066
      *IT_0067*IT_0068*IT_0070 + 0.166666666666667*m_W*(IT_0004*IT_0005 + (-3)
      *IT_0043)*(IT_0001*IT_0069 + IT_0007*IT_0070));
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = -IT_0043 + IT_0074;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*m_W*e_em*(IT_0070*IT_0076
      *IT_0122 + IT_0001*(IT_0008*IT_0069*IT_0078 + -IT_0001*(IT_0007*IT_0070
      *IT_0078 + IT_0001*IT_0069*IT_0122)));
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = IT_0121*IT_0124;
    const ccomplex_t IT_0126 = IT_0119*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = cpow((-2)*s_23 + IT_0025 + IT_0089 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0017
      *IT_0043*(m_W*IT_0001*IT_0007*U_sb_00 + (-0.5)*IT_0066*(IT_0001*IT_0067
      *IT_0068*U_sb_00 + m_b*(m_b*IT_0007*IT_0098*U_sb_00 + mu_h*U_sb_10)))
      *V_ub_mod;
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0016
      *IT_0043*(m_W*IT_0001*IT_0007*conj(U_sb_00) + (-0.5)*IT_0066*(IT_0001
      *IT_0067*IT_0068*conj(U_sb_00) + m_b*(m_b*IT_0007*IT_0098*conj(U_sb_00) +
       mu_h*conj(U_sb_10))))*V_ub_mod;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0130*IT_0132;
    const ccomplex_t IT_0134 = IT_0128*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0021 + IT_0041 + 2*IT_0065 + -IT_0084 +
       IT_0091 + IT_0097 + IT_0104 + IT_0110 + IT_0118 + -IT_0127 + IT_0135;
    return create_ccomplex_return(3*IT_0136*conj(IT_0136));
}

