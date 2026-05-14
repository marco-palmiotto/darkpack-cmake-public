#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_sb_2_to_Hp_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_sb_2_to_Hp_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cpow(IT_0003, -1);
    const ccomplex_t IT_0006 = cpow(IT_0002, -1);
    const ccomplex_t IT_0007 = m_t*conj(U_sb_01)*(m_t*IT_0002*IT_0005*U_st_01 
      + mu_h*U_st_11) + m_b*(m_b*IT_0003*IT_0006*conj(U_sb_01)*U_st_01 + conj
      (U_sb_11)*(mu_h*U_st_01 + m_t*IT_0005*IT_0006*U_st_11));
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*V_tb*(A_b*IT_0003*conj
      (U_sb_11)*U_st_01 + (-1.4142135623731)*e_em*IT_0009*((-0.5)*IT_0004
      *IT_0007 + m_W*IT_0002*IT_0003*conj(U_sb_01)*U_st_01) + A_t*IT_0002*conj
      (U_sb_01)*U_st_11);
    const ccomplex_t IT_0011 = pow(m_Hp, 2);
    const ccomplex_t IT_0012 = pow(m_sb_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = pow(m_st_2, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0011 + IT_0016 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0015 + IT_0019;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0015;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = (-4)*s_23;
    const ccomplex_t IT_0025 = -IT_0000;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (-4)*s_24;
    const ccomplex_t IT_0029 = (-4)*s_34;
    const ccomplex_t IT_0030 = (-4)*IT_0021*(IT_0012*conj(IT_0021) + (-0.25)
      *conj(IT_0023)*IT_0024 + (-0.25)*conj(IT_0027)*IT_0028) + IT_0027*(conj
      (IT_0021)*IT_0028 + conj(IT_0023)*IT_0029) + (-4)*IT_0023*(IT_0011*conj
      (IT_0023) + (-0.25)*conj(IT_0021)*IT_0024 + (-0.25)*conj(IT_0027)*IT_0029);
    return create_ccomplex_return(IT_0030);
}

