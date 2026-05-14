#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_sb_1_to_Hp_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_sb_1_to_Hp_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_01*conj(U_st_01) + IT_0001*IT_0002*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = pow(m_W, -1);
    const ccomplex_t IT_0009 = cpow(IT_0007, -1);
    const ccomplex_t IT_0010 = cpow(IT_0006, -1);
    const ccomplex_t IT_0011 = m_t*conj(U_sb_00)*(m_t*IT_0006*IT_0009*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0007*IT_0010*conj(U_sb_00)*U_st_01 + conj
      (U_sb_10)*(mu_h*U_st_01 + m_t*IT_0009*IT_0010*U_st_11));
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*conj
      (U_sb_10)*U_st_01 + (-1.4142135623731)*e_em*IT_0003*((-0.5)*IT_0008
      *IT_0011 + m_W*IT_0006*IT_0007*conj(U_sb_00)*U_st_01) + A_t*IT_0006*conj
      (U_sb_00)*U_st_11);
    const ccomplex_t IT_0013 = pow(m_Hp, 2);
    const ccomplex_t IT_0014 = pow(m_sb_1, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0005*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*conj(U_st_00)*U_st_01 + IT_0001*IT_0002*conj
      (U_st_10)*U_st_11);
    const ccomplex_t IT_0020 = 1.33333333333333*IT_0019;
    const ccomplex_t IT_0021 = m_t*conj(U_sb_00)*(m_t*IT_0006*IT_0009*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0007*IT_0010*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0009*IT_0010*U_st_10));
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0003*((-0.5)*IT_0008
      *IT_0021 + m_W*IT_0006*IT_0007*conj(U_sb_00)*U_st_00) + A_t*IT_0006*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*e_em*(IT_0001*IT_0002 + 
      -IT_0000*IT_0003);
    const ccomplex_t IT_0028 = pow(m_st_2, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + -IT_0014 + -IT_0028 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0027*IT_0031;
    const ccomplex_t IT_0033 = -IT_0018 + -IT_0026 + IT_0032;
    const ccomplex_t IT_0034 = pow(m_Z, -2);
    const ccomplex_t IT_0035 = pow(s_34, 2);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*U_sb_00*conj(U_sb_00) + IT_0001*IT_0002
      *U_sb_10*conj(U_sb_10));
    const ccomplex_t IT_0038 = (-0.666666666666667)*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0013 + IT_0028 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0002 + 3*IT_0000*IT_0003)*conj(U_sb_00)*U_sb_01 + IT_0001*IT_0002*conj
      (U_sb_10)*U_sb_11);
    const ccomplex_t IT_0044 = (-0.666666666666667)*IT_0043;
    const ccomplex_t IT_0045 = m_t*conj(U_sb_01)*(m_t*IT_0006*IT_0009*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0007*IT_0010*conj(U_sb_01)*U_st_01 + conj
      (U_sb_11)*(mu_h*U_st_01 + m_t*IT_0009*IT_0010*U_st_11));
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0007*conj
      (U_sb_11)*U_st_01 + (-1.4142135623731)*e_em*IT_0003*((-0.5)*IT_0008
      *IT_0045 + m_W*IT_0006*IT_0007*conj(U_sb_01)*U_st_01) + A_t*IT_0006*conj
      (U_sb_01)*U_st_11);
    const ccomplex_t IT_0047 = cpow((-2)*s_13 + IT_0013 + IT_0028 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0044*IT_0049;
    const ccomplex_t IT_0051 = IT_0018 + IT_0026 + IT_0042 + IT_0050;
    const ccomplex_t IT_0052 = s_24*s_34*IT_0034;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = s_23 + IT_0053;
    const ccomplex_t IT_0055 = pow(s_24, 2);
    const ccomplex_t IT_0056 = IT_0034*IT_0055;
    return create_ccomplex_return((-3)*IT_0033*(conj(IT_0033)*(IT_0013 + 
      -IT_0036) + conj(IT_0051)*IT_0054) + (-3)*IT_0051*(conj(IT_0033)*IT_0054 +
       conj(IT_0051)*(IT_0014 + -IT_0056)));
}

