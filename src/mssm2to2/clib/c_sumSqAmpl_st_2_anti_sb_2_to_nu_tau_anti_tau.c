#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_sb_2_to_nu_tau_anti_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_sb_2_to_nu_tau_anti_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_tau = param->m_tau;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_01)*U_st_01;
    const ccomplex_t IT_0003 = (-0.5)*IT_0002;
    const ccomplex_t IT_0004 = m_tau*IT_0003;
    const ccomplex_t IT_0005 = -IT_0002;
    const ccomplex_t IT_0006 = pow(m_W, -2);
    const ccomplex_t IT_0007 = m_tau*IT_0006;
    const ccomplex_t IT_0008 = pow(m_sb_2, 2);
    const ccomplex_t IT_0009 = s_12 + IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0011 = IT_0004 + IT_0010;
    const ccomplex_t IT_0012 = m_tau*IT_0005;
    const ccomplex_t IT_0013 = pow(m_st_2, 2);
    const ccomplex_t IT_0014 = 0.5*IT_0008 + 0.5*IT_0013;
    const ccomplex_t IT_0015 = s_12 + IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0007*IT_0015;
    const ccomplex_t IT_0017 = -IT_0012 + (-2)*IT_0016;
    const ccomplex_t IT_0018 = IT_0011 + IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0008 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = sin(beta);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0001*IT_0025*IT_0027*IT_0028;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = cpow(IT_0028, -1);
    const ccomplex_t IT_0032 = m_t*conj(U_sb_01)*(m_t*IT_0026*IT_0031*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0027*IT_0028*conj(U_sb_01)*U_st_01 + conj
      (U_sb_11)*(mu_h*U_st_01 + m_t*IT_0027*IT_0031*U_st_11));
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0028*conj
      (U_sb_11)*U_st_01 + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0025
      *IT_0032 + m_W*IT_0026*IT_0028*conj(U_sb_01)*U_st_01) + A_t*IT_0026*conj
      (U_sb_01)*U_st_11);
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + -IT_0008 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0024 + -IT_0037;
    const ccomplex_t IT_0039 = IT_0005*IT_0023;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = s_13*m_tau;
    const ccomplex_t IT_0042 = s_13*s_14;
    const ccomplex_t IT_0043 = s_34*IT_0013;
    return create_ccomplex_return(6*IT_0038*(s_34*conj(IT_0038) + -conj
      (IT_0040)*IT_0041) + (-6)*IT_0040*(conj(IT_0038)*IT_0041 + (
      -0.166666666666667)*conj(IT_0040)*(12*IT_0042 + (-6)*IT_0043)));
}

