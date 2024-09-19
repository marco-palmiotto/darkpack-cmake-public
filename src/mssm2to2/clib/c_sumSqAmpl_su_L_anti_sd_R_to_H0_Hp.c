#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_R_to_H0_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_R_to_H0_Hp(
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
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = pow(m_d, 2);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = pow(m_u, 2);
    const ccomplex_t IT_0006 = cpow(IT_0004, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em*
      (IT_0000*IT_0001*IT_0003*IT_0004 + (-2)*IT_0002*(m_W*IT_0004 + (-0.5)
      *IT_0000*IT_0005*IT_0006))*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = cos(alpha);
    const ccomplex_t IT_0012 = sin(alpha);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_td
      *conj(V_td)*IT_0011 + -m_d*e_em*mu_h*IT_0000*IT_0003*IT_0008*IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0014;
    const ccomplex_t IT_0016 = pow(m_H0, 2);
    const ccomplex_t IT_0017 = pow(m_sd_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = pow(m_s, 2);
    const ccomplex_t IT_0022 = (0 + _Complex_I*(-2.82842712474619))*V_us*e_em
      *IT_0008*(IT_0002*(m_W*IT_0004 + (-0.5)*IT_0000*IT_0005*IT_0006) + -1./2
      *IT_0000*IT_0003*IT_0004*IT_0021);
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*A_b*conj(V_td)
      *V_ts*IT_0011;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0008
      *IT_0030*(m_W*IT_0002*IT_0004*conj(U_sb_01) + (-0.5)*IT_0000*(IT_0002
      *IT_0005*IT_0006*conj(U_sb_01) + m_b*(m_b*IT_0003*IT_0004*conj(U_sb_01) +
       mu_h*conj(U_sb_11))))*V_ub_mod;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_01;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0000*IT_0008;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -2);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0005
      *IT_0006*IT_0008*IT_0012 + 0.5*m_W*(IT_0002*IT_0011 + -IT_0004*IT_0012)*
      (IT_0008 + -1./3*IT_0007*IT_0042));
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0040*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_13 + IT_0016 + (0 + _Complex_I*1)
      *m_su_L*Gamma_ul + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*-3)*e_em*(IT_0000*IT_0001
      *IT_0003*IT_0008*IT_0011 + -1./3*m_W*IT_0007*(IT_0002*IT_0011 + -IT_0004
      *IT_0012)*IT_0042);
    const ccomplex_t IT_0050 = 0.333333333333333*IT_0049;
    const ccomplex_t IT_0051 = IT_0040*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0016 + (0 + _Complex_I*1)
      *m_sd_R*Gamma_dr + reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0016 + IT_0017 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0008
      *IT_0030*(m_W*IT_0002*IT_0004*conj(U_sb_00) + (-0.5)*IT_0000*(IT_0002
      *IT_0005*IT_0006*conj(U_sb_00) + m_b*(m_b*IT_0003*IT_0004*conj(U_sb_00) +
       mu_h*conj(U_sb_10))))*V_ub_mod;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0011*U_sb_00;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = IT_0055*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = cpow(IT_0004, 3);
    const ccomplex_t IT_0064 = cpow(IT_0002, 2);
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_W*e_em*(IT_0008*IT_0011
      *IT_0063 + (-2)*IT_0002*(IT_0004*(IT_0002*IT_0011 + -IT_0004*IT_0012)*(0.5
      *IT_0008 + IT_0007*IT_0042) + 1./2*IT_0008*IT_0012*IT_0064));
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0000*IT_0003*IT_0004*IT_0008;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = pow(m_W, 2);
    const ccomplex_t IT_0071 = pow(m_su_L, 2);
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + -IT_0017 + IT_0070 + -IT_0071 
      + -reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0069*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_12 + -IT_0017 + -IT_0071 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0076 = -IT_0008;
    const ccomplex_t IT_0077 = IT_0007*IT_0042;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = -IT_0002;
    const ccomplex_t IT_0080 = 3*IT_0008 + IT_0077;
    const ccomplex_t IT_0081 = cpow(IT_0004, 2);
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*m_W*e_em*(IT_0012*IT_0063
      *IT_0078 + IT_0079*(IT_0079*(IT_0002*IT_0011*IT_0078 + -IT_0004*IT_0012
      *IT_0080) + IT_0011*IT_0080*IT_0081));
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0040*IT_0083;
    const ccomplex_t IT_0085 = IT_0075*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0000*IT_0003*IT_0006*IT_0008;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0000
      *IT_0006*IT_0008*IT_0011;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_13 + IT_0016 + IT_0071 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0020 + IT_0029 + IT_0038 + IT_0048 + IT_0054
       + IT_0062 + -IT_0074 + -IT_0086 + IT_0094;
    return create_ccomplex_return(3*IT_0095*conj(IT_0095));
}

