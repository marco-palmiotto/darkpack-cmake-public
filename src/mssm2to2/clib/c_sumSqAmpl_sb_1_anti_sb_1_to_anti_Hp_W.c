#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_1_anti_sb_1_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_1_anti_sb_1_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
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
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_sb_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0002 = cpow(V_tb, 2);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_b*mu_h*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0003*IT_0006;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = pow(m_b, 2);
    const ccomplex_t IT_0018 = 2*U_sb_00*(((-0.5)*m_W*IT_0014*(IT_0010 + 1./3
      *IT_0009*IT_0016) + IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj(U_sb_00)
       + 1./2*IT_0005*IT_0007*IT_0008*IT_0010*conj(U_sb_10)) + U_sb_10*(IT_0005
      *IT_0007*IT_0008*IT_0010*conj(U_sb_00) + (-0.666666666666667)*(m_W*IT_0009
      *IT_0014*IT_0016 + (-3)*IT_0003*IT_0004*IT_0007*IT_0010*IT_0017)*conj
      (U_sb_10));
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(e_em*IT_0018 +
       1.4142135623731*A_b*IT_0002*IT_0003*(conj(U_sb_00)*U_sb_10 + U_sb_00*conj
      (U_sb_10)));
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*e_em*IT_0010*(IT_0006
      *IT_0008 + IT_0003*IT_0012);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0027 = -IT_0003*IT_0012;
    const ccomplex_t IT_0028 = IT_0006*IT_0008;
    const ccomplex_t IT_0029 = IT_0027 + IT_0028;
    const ccomplex_t IT_0030 = (-2)*U_sb_00*((IT_0004*IT_0007*IT_0008*IT_0010
      *IT_0017 + (-0.5)*m_W*(IT_0010 + 1./3*IT_0009*IT_0016)*IT_0029)*conj
      (U_sb_00) + -1./2*IT_0003*IT_0005*IT_0007*IT_0010*conj(U_sb_10)) + U_sb_10
      *(IT_0003*IT_0005*IT_0007*IT_0010*conj(U_sb_00) + (-2)*(IT_0004*IT_0007
      *IT_0008*IT_0010*IT_0017 + -1./3*m_W*IT_0009*IT_0016*IT_0029)*conj(U_sb_10
      ));
    const ccomplex_t IT_0031 = (0 + _Complex_I*-1)*(e_em*IT_0030 + (
      -1.4142135623731)*A_b*IT_0002*IT_0008*(conj(U_sb_00)*U_sb_10 + U_sb_00
      *conj(U_sb_10)));
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0026*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*-1)*e_em*IT_0010*(IT_0003
      *IT_0006 + -IT_0008*IT_0012);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = pow(m_Hp, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0039 = pow(m_c, 2);
    const ccomplex_t IT_0040 = cpow(IT_0012, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010*(m_W*IT_0006*IT_0012*U_sb_00 + (-0.5)*IT_0004*(IT_0006*IT_0039
      *IT_0040*U_sb_00 + m_b*(m_b*IT_0007*IT_0012*U_sb_00 + mu_h*U_sb_10)));
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0038*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0010*conj(U_sb_00);
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0049 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0050 = pow(m_u, 2);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0010
      *IT_0049*(m_W*IT_0006*IT_0012*U_sb_00 + (-0.5)*IT_0004*(IT_0006*IT_0040
      *IT_0050*U_sb_00 + m_b*(m_b*IT_0007*IT_0012*U_sb_00 + mu_h*U_sb_10)))
      *V_ub_mod;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = IT_0048*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0010
      *IT_0055*conj(U_sb_00)*V_ub_mod;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0054*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0060 = m_t*U_sb_00*(m_t*IT_0006*IT_0040*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0007*IT_0012*U_sb_00*conj(U_st_00) +
       U_sb_10*(mu_h*conj(U_st_00) + m_t*IT_0007*IT_0040*conj(U_st_10)));
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0012*U_sb_10
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0010*((-0.5)*IT_0004*IT_0060 +
       m_W*IT_0006*IT_0012*U_sb_00*conj(U_st_00)) + A_t*IT_0006*U_sb_00*conj
      (U_st_10));
    const ccomplex_t IT_0062 = IT_0059*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0010*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = IT_0063*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_13 + IT_0000 + IT_0037 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0068 = m_t*U_sb_00*(m_t*IT_0006*IT_0040*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0007*IT_0012*U_sb_00*conj(U_st_01) +
       U_sb_10*(mu_h*conj(U_st_01) + m_t*IT_0007*IT_0040*conj(U_st_11)));
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0012*U_sb_10
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0010*((-0.5)*IT_0004*IT_0068 +
       m_W*IT_0006*IT_0012*U_sb_00*conj(U_st_01)) + A_t*IT_0006*U_sb_00*conj
      (U_st_11));
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0010*conj(U_sb_00)*U_st_01;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = IT_0071*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0076 = m_b*e_em*mu_h*IT_0004*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0077 = IT_0010*IT_0076;
    const ccomplex_t IT_0078 = m_b*e_em*mu_h*IT_0004*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0079 = IT_0010*IT_0078;
    const ccomplex_t IT_0080 = A_b*IT_0002*conj(U_sb_00)*U_sb_10;
    const ccomplex_t IT_0081 = IT_0012*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = A_b*IT_0002*U_sb_00*conj(U_sb_10);
    const ccomplex_t IT_0084 = IT_0012*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0077 + 0.5*IT_0079 + (-0.5)*IT_0082 +
       0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0075*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = e_em*IT_0010;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = IT_0088*IT_0090;
    const ccomplex_t IT_0092 = -IT_0025 + -IT_0036 + IT_0047 + IT_0058 +
       IT_0066 + IT_0074 + -IT_0091;
    const ccomplex_t IT_0093 = pow(m_W, -2);
    const ccomplex_t IT_0094 = pow(s_24, 2);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = -IT_0025 + -IT_0036 + -IT_0091;
    const ccomplex_t IT_0097 = s_14*s_24*IT_0093;
    const ccomplex_t IT_0098 = (-3)*s_12 + 3*IT_0097;
    const ccomplex_t IT_0099 = pow(s_14, 2);
    const ccomplex_t IT_0100 = IT_0093*IT_0099;
    return create_ccomplex_return((-3)*IT_0092*(conj(IT_0092)*(IT_0000 + 
      -IT_0095) + (-0.333333333333333)*conj(IT_0096)*IT_0098) + (-3)*IT_0096*((
      -0.333333333333333)*conj(IT_0092)*IT_0098 + conj(IT_0096)*(IT_0000 + 
      -IT_0100)));
}

