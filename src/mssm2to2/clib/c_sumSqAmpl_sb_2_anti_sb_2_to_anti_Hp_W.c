#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sb_2_anti_sb_2_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sb_2_anti_sb_2_to_anti_Hp_W(
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
    const creal_t m_sb_2 = param->m_sb_2;
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
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_b*e_em*mu_h*IT_0004*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = m_b*e_em*mu_h*IT_0004*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = cpow(V_tb, 2);
    const ccomplex_t IT_0011 = A_b*IT_0010*conj(U_sb_01)*U_sb_11;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = A_b*IT_0010*U_sb_01*conj(U_sb_11);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0006 + 0.5*IT_0008 + (-0.5)*IT_0013 +
       0.5*IT_0016;
    const ccomplex_t IT_0018 = pow(m_sb_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + (-2)*IT_0018 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0003*IT_0021;
    const ccomplex_t IT_0023 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0023*conj(U_sb_01)*V_ub_mod;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0027 = cos(beta);
    const ccomplex_t IT_0028 = pow(m_u, 2);
    const ccomplex_t IT_0029 = cpow(IT_0009, -1);
    const ccomplex_t IT_0030 = cpow(IT_0027, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0026*(m_W*IT_0009*IT_0027*U_sb_01 + (-0.5)*IT_0004*(IT_0027*IT_0028
      *IT_0029*U_sb_01 + m_b*(m_b*IT_0009*IT_0030*U_sb_01 + mu_h*U_sb_11)))
      *V_ub_mod;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = pow(m_Hp, 2);
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0018 + IT_0033 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0025*IT_0036;
    const ccomplex_t IT_0038 = cos(alpha);
    const ccomplex_t IT_0039 = sin(alpha);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0009
      *IT_0038 + -IT_0027*IT_0039);
    const ccomplex_t IT_0041 = m_b*mu_h*IT_0004;
    const ccomplex_t IT_0042 = -IT_0009*IT_0039;
    const ccomplex_t IT_0043 = IT_0027*IT_0038;
    const ccomplex_t IT_0044 = IT_0042 + IT_0043;
    const ccomplex_t IT_0045 = cos(theta_W);
    const ccomplex_t IT_0046 = cpow(IT_0045, -2);
    const ccomplex_t IT_0047 = pow(m_b, 2);
    const ccomplex_t IT_0048 = (-2)*U_sb_01*(((-0.5)*m_W*IT_0044*(IT_0001 + 1.
      /3*IT_0000*IT_0046) + IT_0001*IT_0004*IT_0030*IT_0038*IT_0047)*conj
      (U_sb_01) + -1./2*IT_0001*IT_0030*IT_0039*IT_0041*conj(U_sb_11)) + U_sb_11
      *(IT_0001*IT_0030*IT_0039*IT_0041*conj(U_sb_01) + 0.666666666666667*(m_W
      *IT_0000*IT_0044*IT_0046 + (-3)*IT_0001*IT_0004*IT_0030*IT_0038*IT_0047)
      *conj(U_sb_11));
    const ccomplex_t IT_0049 = (0 + _Complex_I*-1)*(e_em*IT_0048 + (
      -1.4142135623731)*A_b*IT_0010*IT_0038*(conj(U_sb_01)*U_sb_11 + U_sb_01
      *conj(U_sb_11)));
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0040*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0027
      *IT_0038 + IT_0009*IT_0039);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0027*IT_0039;
    const ccomplex_t IT_0058 = IT_0009*IT_0038;
    const ccomplex_t IT_0059 = IT_0057 + IT_0058;
    const ccomplex_t IT_0060 = 2*U_sb_01*((IT_0001*IT_0004*IT_0030*IT_0039
      *IT_0047 + (-0.5)*m_W*(IT_0001 + 1./3*IT_0000*IT_0046)*IT_0059)*conj
      (U_sb_01) + 1./2*IT_0001*IT_0030*IT_0038*IT_0041*conj(U_sb_11)) + U_sb_11*
      (IT_0001*IT_0030*IT_0038*IT_0041*conj(U_sb_01) + 2*(IT_0001*IT_0004
      *IT_0030*IT_0039*IT_0047 + -1./3*m_W*IT_0000*IT_0046*IT_0059)*conj(U_sb_11
      ));
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(e_em*IT_0060 +
       1.4142135623731*A_b*IT_0010*IT_0039*(conj(U_sb_01)*U_sb_11 + U_sb_01*conj
      (U_sb_11)));
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0056*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_00;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = m_t*U_sb_01*(m_t*IT_0027*IT_0029*conj(U_st_00) 
      + mu_h*conj(U_st_10)) + m_b*(m_b*IT_0009*IT_0030*U_sb_01*conj(U_st_00) +
       U_sb_11*(mu_h*conj(U_st_00) + m_t*IT_0029*IT_0030*conj(U_st_10)));
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0009*U_sb_11
      *conj(U_st_00) + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0004*IT_0069 +
       m_W*IT_0009*IT_0027*U_sb_01*conj(U_st_00)) + A_t*IT_0027*U_sb_01*conj
      (U_st_10));
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0018 + IT_0033 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0068*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = m_t*U_sb_01*(m_t*IT_0027*IT_0029*conj(U_st_01) 
      + mu_h*conj(U_st_11)) + m_b*(m_b*IT_0009*IT_0030*U_sb_01*conj(U_st_01) +
       U_sb_11*(mu_h*conj(U_st_01) + m_t*IT_0029*IT_0030*conj(U_st_11)));
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0009*U_sb_11
      *conj(U_st_01) + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0004*IT_0077 +
       m_W*IT_0009*IT_0027*U_sb_01*conj(U_st_01)) + A_t*IT_0027*U_sb_01*conj
      (U_st_11));
    const ccomplex_t IT_0079 = cpow((-2)*s_13 + IT_0018 + IT_0033 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = IT_0076*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0001*conj(U_sb_01);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = pow(m_c, 2);
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0001*(m_W*IT_0009*IT_0027*U_sb_01 + (-0.5)*IT_0004*(IT_0027*IT_0029
      *IT_0085*U_sb_01 + m_b*(m_b*IT_0009*IT_0030*U_sb_01 + mu_h*U_sb_11)));
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0018 + IT_0033 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0084*IT_0090;
    const ccomplex_t IT_0092 = -IT_0022 + IT_0037 + -IT_0054 + -IT_0066 +
       IT_0074 + IT_0082 + IT_0091;
    const ccomplex_t IT_0093 = pow(m_W, -2);
    const ccomplex_t IT_0094 = pow(s_24, 2);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = -IT_0022 + -IT_0054 + -IT_0066;
    const ccomplex_t IT_0097 = s_14*s_24*IT_0093;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = pow(s_14, 2);
    const ccomplex_t IT_0100 = IT_0093*IT_0099;
    return create_ccomplex_return((-3)*IT_0092*(conj(IT_0092)*(IT_0018 + 
      -IT_0095) + conj(IT_0096)*(s_12 + IT_0098)) + (-3)*IT_0096*(conj(IT_0092)*
      (s_12 + IT_0098) + conj(IT_0096)*(IT_0018 + -IT_0100)));
}

