#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_R_to_H0_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_R_to_H0_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0002
      *IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = pow(m_H0, 2);
    const ccomplex_t IT_0011 = pow(m_se_R, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0010 + IT_0011 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = pow(m_se_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0010 + IT_0015 + -m_se_R*
      (m_se_R + (0 + _Complex_I*-1)*Gamma_er) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0014 + IT_0018;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = IT_0001*IT_0014;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = -s_23;
    const ccomplex_t IT_0024 = IT_0000*IT_0019;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = -s_24;
    const ccomplex_t IT_0027 = -s_34;
    const ccomplex_t IT_0028 = -IT_0020*(IT_0011*conj(IT_0020) + -conj(IT_0022
      )*IT_0023 + -conj(IT_0025)*IT_0026) + IT_0025*(conj(IT_0020)*IT_0026 +
       conj(IT_0022)*IT_0027) + -IT_0022*(IT_0010*conj(IT_0022) + -conj(IT_0020)
      *IT_0023 + -conj(IT_0025)*IT_0027);
    return create_ccomplex_return(IT_0028);
}

