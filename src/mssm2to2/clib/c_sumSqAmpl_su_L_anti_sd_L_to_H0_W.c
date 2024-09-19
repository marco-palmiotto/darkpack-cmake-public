#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_sd_L_to_H0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_sd_L_to_H0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001;
    const ccomplex_t IT_0003 = (-0.5)*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*(IT_0004
      *IT_0005 + IT_0006*IT_0007);
    const ccomplex_t IT_0009 = pow(m_sd_L, 2);
    const ccomplex_t IT_0010 = pow(m_su_L, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0009 + -IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -2);
    const ccomplex_t IT_0017 = pow(s_14, 2);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0019*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0005*U_sb_11;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = pow(m_H0, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0009 + IT_0024 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0019*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*A_b*V_tb*conj
      (V_td)*IT_0005*U_sb_10;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0009 + IT_0024 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0030*IT_0035;
    const ccomplex_t IT_0037 = -IT_0002;
    const ccomplex_t IT_0038 = pow(m_W, -1);
    const ccomplex_t IT_0039 = pow(m_d, 2);
    const ccomplex_t IT_0040 = cpow(IT_0004, -1);
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -2);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0005
      *IT_0038*IT_0039*IT_0040 + (-0.5)*m_W*(IT_0004*IT_0005 + -IT_0006*IT_0007)
      *(IT_0001 + 1./3*IT_0000*IT_0042));
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0024 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0037*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0004
      *IT_0005 + IT_0006*IT_0007);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = pow(m_u, 2);
    const ccomplex_t IT_0052 = cpow(IT_0004, 2);
    const ccomplex_t IT_0053 = cpow(IT_0006, 2);
    const ccomplex_t IT_0054 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001*(IT_0038*(IT_0039 + -IT_0051) + -m_W*(IT_0052 + -IT_0053));
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = pow(m_W, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + IT_0056 
      + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0050*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0005
      *IT_0006 + -IT_0004*IT_0007);
    const ccomplex_t IT_0062 = cpow(IT_0006, -1);
    const ccomplex_t IT_0063 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0001*(IT_0006*IT_0038*IT_0039*IT_0040 + (-2)*IT_0004*(m_W*IT_0006 + (
      -0.5)*IT_0038*IT_0051*IT_0062));
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + -IT_0009 + -IT_0010 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0061*IT_0067;
    const ccomplex_t IT_0069 = -IT_0028 + -IT_0036 + -IT_0048 + IT_0060 +
       IT_0068;
    const ccomplex_t IT_0070 = s_14*s_34*IT_0016;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = s_13 + IT_0071;
    const ccomplex_t IT_0073 = -IT_0037;
    const ccomplex_t IT_0074 = IT_0003 + IT_0073;
    const ccomplex_t IT_0075 = IT_0013*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*e_em*(0.5*m_W*(IT_0004
      *IT_0005 + -IT_0006*IT_0007)*(IT_0001 + -1./3*IT_0000*IT_0042) + IT_0001
      *IT_0007*IT_0038*IT_0051*IT_0062);
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0024 + (0 + _Complex_I*1)
      *m_su_L*Gamma_ul + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0037*IT_0080;
    const ccomplex_t IT_0082 = IT_0028 + IT_0036 + IT_0048 + -IT_0075 + IT_0081;
    const ccomplex_t IT_0083 = s_14*s_24*IT_0016;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = s_12 + IT_0084;
    const ccomplex_t IT_0086 = (-3)*conj(IT_0015);
    const ccomplex_t IT_0087 = pow(s_34, 2);
    const ccomplex_t IT_0088 = IT_0016*IT_0087;
    const ccomplex_t IT_0089 = s_24*s_34*IT_0016;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = s_23 + IT_0090;
    const ccomplex_t IT_0092 = pow(s_24, 2);
    const ccomplex_t IT_0093 = IT_0016*IT_0092;
    const ccomplex_t IT_0094 = (-3)*IT_0082;
    const ccomplex_t IT_0095 = (-3)*IT_0015*(conj(IT_0015)*(IT_0010 + -IT_0018
      ) + conj(IT_0069)*IT_0072 + conj(IT_0082)*IT_0085) + IT_0082*IT_0085
      *IT_0086 + (-3)*IT_0069*((-0.333333333333333)*IT_0072*IT_0086 + conj
      (IT_0069)*(IT_0024 + -IT_0088) + conj(IT_0082)*IT_0091) + (conj(IT_0069)
      *IT_0091 + conj(IT_0082)*(IT_0009 + -IT_0093))*IT_0094;
    return create_ccomplex_return(IT_0095);
}

