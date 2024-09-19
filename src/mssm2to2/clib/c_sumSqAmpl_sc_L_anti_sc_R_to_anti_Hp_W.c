#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sc_R_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sc_R_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t V_cb = param->V_cb;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = cpow(IT_0001, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0000
      *IT_0005*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = pow(m_sc_L, 2);
    const ccomplex_t IT_0012 = pow(m_sc_R, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_c*e_em*mu_h*IT_0003
      *IT_0005*IT_0007*IT_0008;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = e_em*IT_0005;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = m_c*e_em*mu_h*IT_0005*IT_0007;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = -IT_0016 + -IT_0024 + -IT_0032;
    const ccomplex_t IT_0034 = pow(m_W, -2);
    const ccomplex_t IT_0035 = pow(s_14, 2);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*IT_0005;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1.4142135623731)*m_c*V_cd*e_em
      *mu_h*IT_0005*IT_0007;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = pow(m_Hp, 2);
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0012 + IT_0041 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0038*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0005;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*e_em
      *mu_h*IT_0005*IT_0007;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_23 + IT_0012 + IT_0041 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0047*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0005*conj(U_sb_00);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = cpow(IT_0002, -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0005*IT_0007*(mu_h*U_sb_00 + m_b*IT_0008*IT_0056*U_sb_10);
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0012 + IT_0041 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0055*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0005*conj(U_sb_01);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1.4142135623731)*m_c*V_cb*e_em
      *IT_0005*IT_0007*(mu_h*U_sb_01 + m_b*IT_0008*IT_0056*U_sb_11);
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_23 + IT_0012 + IT_0041 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0064*IT_0069;
    const ccomplex_t IT_0071 = -IT_0045 + -IT_0053 + -IT_0062 + -IT_0070;
    const ccomplex_t IT_0072 = s_14*s_34*IT_0034;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = s_13 + IT_0073;
    const ccomplex_t IT_0075 = -IT_0016 + -IT_0024 + -IT_0032 + IT_0045 +
       IT_0053 + IT_0062 + IT_0070;
    const ccomplex_t IT_0076 = s_14*s_24*IT_0034;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = s_12 + IT_0077;
    const ccomplex_t IT_0079 = (-3)*conj(IT_0033);
    const ccomplex_t IT_0080 = pow(s_34, 2);
    const ccomplex_t IT_0081 = IT_0034*IT_0080;
    const ccomplex_t IT_0082 = s_24*s_34*IT_0034;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = s_23 + IT_0083;
    const ccomplex_t IT_0085 = pow(s_24, 2);
    const ccomplex_t IT_0086 = IT_0034*IT_0085;
    const ccomplex_t IT_0087 = (-3)*IT_0075;
    const ccomplex_t IT_0088 = (-3)*IT_0033*(conj(IT_0033)*(IT_0011 + -IT_0036
      ) + conj(IT_0071)*IT_0074 + conj(IT_0075)*IT_0078) + IT_0075*IT_0078
      *IT_0079 + (-3)*IT_0071*((-0.333333333333333)*IT_0074*IT_0079 + conj
      (IT_0071)*(IT_0041 + -IT_0081) + conj(IT_0075)*IT_0084) + (conj(IT_0071)
      *IT_0084 + conj(IT_0075)*(IT_0012 + -IT_0086))*IT_0087;
    return create_ccomplex_return(IT_0088);
}

