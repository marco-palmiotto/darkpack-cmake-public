#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_A0_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_A0_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_t*e_em*mu_h*IT_0004*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = m_t*e_em*mu_h*IT_0004*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0006 + 0.5*IT_0008 + (-0.5)*IT_0012 +
       0.5*IT_0015;
    const ccomplex_t IT_0017 = pow(m_A0, 2);
    const ccomplex_t IT_0018 = pow(m_st_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0017 + IT_0018 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = pow(m_st_1, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0017 + IT_0022 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0016*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0021 + IT_0025;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = IT_0001*IT_0021;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-4)*s_23;
    const ccomplex_t IT_0031 = -IT_0000;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (-4)*s_24;
    const ccomplex_t IT_0035 = (-4)*s_34;
    const ccomplex_t IT_0036 = (-4)*IT_0027*(IT_0018*conj(IT_0027) + (-0.25)
      *conj(IT_0029)*IT_0030 + (-0.25)*conj(IT_0033)*IT_0034) + IT_0033*(conj
      (IT_0027)*IT_0034 + conj(IT_0029)*IT_0035) + (-4)*IT_0029*(IT_0017*conj
      (IT_0029) + (-0.25)*conj(IT_0027)*IT_0030 + (-0.25)*conj(IT_0033)*IT_0035);
    return create_ccomplex_return(IT_0036);
}

