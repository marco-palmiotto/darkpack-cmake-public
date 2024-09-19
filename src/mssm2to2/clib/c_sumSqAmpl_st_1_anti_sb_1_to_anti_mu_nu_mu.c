#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_sb_1_to_anti_mu_nu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_sb_1_to_anti_mu_nu_mu(
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
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = pow(m_sb_1, 2);
    const ccomplex_t IT_0005 = pow(m_st_1, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0004 + -IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0001*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = m_mu*IT_0010;
    const ccomplex_t IT_0012 = -IT_0009;
    const ccomplex_t IT_0013 = m_mu*IT_0012;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = IT_0011 + IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + -IT_0004 + -IT_0005 + IT_0017 
      + -reg_prop, -1);
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0001*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = cos(beta);
    const ccomplex_t IT_0023 = sin(beta);
    const ccomplex_t IT_0024 = cpow(IT_0022, -1);
    const ccomplex_t IT_0025 = cpow(IT_0023, -1);
    const ccomplex_t IT_0026 = pow(m_b, 2);
    const ccomplex_t IT_0027 = pow(m_t, 2);
    const ccomplex_t IT_0028 = cpow(IT_0022, 2);
    const ccomplex_t IT_0029 = cpow(IT_0023, 2);
    const ccomplex_t IT_0030 = m_b*mu_h*IT_0019*IT_0023*IT_0024*conj(U_sb_10)
      *U_st_00 + -conj(U_sb_00)*((IT_0019*(IT_0026 + -IT_0027) + -m_W*(IT_0028 +
       -IT_0029))*U_st_00 + m_t*mu_h*IT_0019*IT_0022*IT_0025*U_st_10);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0001*IT_0030 + A_b*IT_0022*conj(U_sb_10)*U_st_00 + -A_t*IT_0023
      *conj(U_sb_00)*U_st_10);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = IT_0021*IT_0032;
    const ccomplex_t IT_0034 = IT_0018*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0001*IT_0019*IT_0023*IT_0024;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = m_t*conj(U_sb_00)*(m_t*IT_0022*IT_0025*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0023*IT_0024*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0024*IT_0025*U_st_10));
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0023*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0001*((-0.5)*IT_0019
      *IT_0038 + m_W*IT_0022*IT_0023*conj(U_sb_00)*U_st_00) + A_t*IT_0022*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + -IT_0004 + -IT_0005 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0016 + -IT_0035 + -IT_0043;
    const ccomplex_t IT_0045 = IT_0008*IT_0012;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = m_mu*s_14;
    const ccomplex_t IT_0048 = s_13*s_14;
    const ccomplex_t IT_0049 = s_34*IT_0005;
    return create_ccomplex_return(6*IT_0044*(s_34*conj(IT_0044) + -conj
      (IT_0046)*IT_0047) + (-6)*IT_0046*(conj(IT_0044)*IT_0047 + (
      -0.166666666666667)*conj(IT_0046)*(12*IT_0048 + (-6)*IT_0049)));
}

