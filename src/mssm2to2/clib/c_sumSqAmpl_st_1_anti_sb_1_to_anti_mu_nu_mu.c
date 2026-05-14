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
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0003, -1);
    const ccomplex_t IT_0009 = m_t*conj(U_sb_00)*(m_t*IT_0001*IT_0008*U_st_00 
      + mu_h*U_st_10) + m_b*(m_b*IT_0002*IT_0003*conj(U_sb_00)*U_st_00 + conj
      (U_sb_10)*(mu_h*U_st_00 + m_t*IT_0002*IT_0008*U_st_10));
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0003*conj
      (U_sb_10)*U_st_00 + (-1.4142135623731)*e_em*IT_0005*((-0.5)*IT_0000
      *IT_0009 + m_W*IT_0001*IT_0003*conj(U_sb_00)*U_st_00) + A_t*IT_0001*conj
      (U_sb_00)*U_st_10);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = pow(m_sb_1, 2);
    const ccomplex_t IT_0013 = pow(m_st_1, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*V_tb*e_em
      *IT_0005*conj(U_sb_00)*U_st_00;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = m_mu*IT_0023;
    const ccomplex_t IT_0025 = -IT_0022;
    const ccomplex_t IT_0026 = m_mu*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu
      *IT_0000*IT_0005;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = pow(m_b, 2);
    const ccomplex_t IT_0033 = pow(m_t, 2);
    const ccomplex_t IT_0034 = cpow(IT_0001, 2);
    const ccomplex_t IT_0035 = cpow(IT_0003, 2);
    const ccomplex_t IT_0036 = m_b*mu_h*IT_0000*IT_0002*IT_0003*conj(U_sb_10)
      *U_st_00 + -conj(U_sb_00)*((IT_0000*(IT_0032 + -IT_0033) + -m_W*(IT_0034 +
       -IT_0035))*U_st_00 + m_t*mu_h*IT_0000*IT_0001*IT_0008*U_st_10);
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*V_tb*((-0.707106781186548)
      *e_em*IT_0005*IT_0036 + A_b*IT_0001*conj(U_sb_10)*U_st_00 + -A_t*IT_0003
      *conj(U_sb_00)*U_st_10);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = pow(m_W, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0040 
      + -reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = -IT_0016 + IT_0029 + -IT_0043;
    const ccomplex_t IT_0045 = IT_0021*IT_0025;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = m_mu*s_14;
    const ccomplex_t IT_0048 = s_13*s_14;
    const ccomplex_t IT_0049 = s_34*IT_0013;
    return create_ccomplex_return(6*IT_0044*(s_34*conj(IT_0044) + -conj
      (IT_0046)*IT_0047) + (-6)*IT_0046*(conj(IT_0044)*IT_0047 + (
      -0.166666666666667)*conj(IT_0046)*(12*IT_0048 + (-6)*IT_0049)));
}

