#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_L_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_L_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t m_su_L = param->m_su_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = pow(m_d, 2);
    const ccomplex_t IT_0009 = cpow(IT_0002, -1);
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0005
      *IT_0007*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0000*IT_0002 + -IT_0001*IT_0003)
      *(IT_0005 + 1./3*IT_0004*IT_0011));
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_sd_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + (-2)*IT_0014 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0005
      *IT_0007*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0000*IT_0001 + IT_0002*IT_0003)*
      (IT_0005 + 1./3*IT_0004*IT_0011));
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0020*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0005;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = pow(m_u, 2);
    const ccomplex_t IT_0029 = cpow(IT_0001, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0005*(IT_0001*IT_0007*IT_0008*IT_0009 + (-2)*IT_0002*(m_W*IT_0001 + (
      -0.5)*IT_0007*IT_0028*IT_0029));
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = pow(m_Hp, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0014 + IT_0032 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = IT_0027*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0005*U_st_00;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = A_t*IT_0002*conj(U_st_10) + (-1.4142135623731)
      *e_em*IT_0005*(m_W*IT_0001*IT_0002*conj(U_st_00) + (-0.5)*IT_0007*(IT_0001
      *IT_0008*IT_0009*conj(U_st_00) + m_t*(m_t*IT_0002*IT_0029*conj(U_st_00) +
       mu_h*conj(U_st_10))));
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*V_td*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0014 + IT_0032 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*conj(V_td)
      *e_em*IT_0005*U_st_01;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = A_t*IT_0002*conj(U_st_11) + (-1.4142135623731)
      *e_em*IT_0005*(m_W*IT_0001*IT_0002*conj(U_st_01) + (-0.5)*IT_0007*(IT_0001
      *IT_0008*IT_0009*conj(U_st_01) + m_t*(m_t*IT_0002*IT_0029*conj(U_st_01) +
       mu_h*conj(U_st_11))));
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*V_td*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0014 + IT_0032 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*IT_0005;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = pow(m_c, 2);
    const ccomplex_t IT_0056 = (0 + _Complex_I*1.4142135623731)*V_cd*e_em
      *IT_0005*(IT_0001*IT_0007*IT_0008*IT_0009 + (-2)*IT_0002*(m_W*IT_0001 + (
      -0.5)*IT_0007*IT_0029*IT_0055));
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_13 + IT_0014 + IT_0032 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0054*IT_0060;
    const ccomplex_t IT_0062 = -IT_0018 + -IT_0025 + IT_0036 + IT_0044 +
       IT_0052 + IT_0061;
    const ccomplex_t IT_0063 = pow(m_W, -2);
    const ccomplex_t IT_0064 = pow(s_24, 2);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = -IT_0018 + -IT_0025;
    const ccomplex_t IT_0067 = s_14*s_24*IT_0063;
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = pow(s_14, 2);
    const ccomplex_t IT_0070 = IT_0063*IT_0069;
    return create_ccomplex_return((-3)*IT_0062*(conj(IT_0062)*(IT_0014 + 
      -IT_0065) + conj(IT_0066)*(s_12 + IT_0068)) + (-3)*IT_0066*(conj(IT_0062)*
      (s_12 + IT_0068) + conj(IT_0066)*(IT_0014 + -IT_0070)));
}

