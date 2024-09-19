#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_snu_mu_to_A0_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_snu_mu_to_A0_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = pow(m_W, -2);
    const ccomplex_t IT_0005 = pow(m_mu, 2);
    const ccomplex_t IT_0006 = cpow(IT_0001, -2);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0000*(IT_0002 + -IT_0003*
      (IT_0003 + -IT_0003*IT_0004*IT_0005*IT_0006))*IT_0008;
    const ccomplex_t IT_0010 = (-0.25)*IT_0009;
    const ccomplex_t IT_0011 = s_13 + s_23;
    const ccomplex_t IT_0012 = cpow(IT_0007, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = e_em*IT_0012;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = IT_0011*IT_0017;
    const ccomplex_t IT_0019 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0020 = s_12 + IT_0019;
    const ccomplex_t IT_0021 = -IT_0013;
    const ccomplex_t IT_0022 = (-0.5)*IT_0015;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = IT_0018 + IT_0024;
    const ccomplex_t IT_0026 = IT_0016*IT_0021;
    const ccomplex_t IT_0027 = s_23*IT_0026;
    const ccomplex_t IT_0028 = IT_0014*IT_0022;
    const ccomplex_t IT_0029 = pow(m_smu_L, 2);
    const ccomplex_t IT_0030 = 0.5*IT_0019 + 0.5*IT_0029;
    const ccomplex_t IT_0031 = s_12 + IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = -IT_0027 + (-2)*IT_0032;
    const ccomplex_t IT_0034 = IT_0025 + IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0019 + -IT_0029 + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = pow(m_W, -1);
    const ccomplex_t IT_0039 = e_em*m_mu*mu_h*IT_0012*IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0012*IT_0038;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = pow(m_A0, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0029 + IT_0044 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = pow(m_W, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0019 + -IT_0029 + IT_0048 
      + -reg_prop, -1);
    const ccomplex_t IT_0050 = cpow(IT_0003, 2);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0012*
      (IT_0005*IT_0038 + -m_W*(IT_0002 + -IT_0050));
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = m_W*e_em*IT_0012;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = IT_0049*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0010 + IT_0037 + IT_0047 + -IT_0057;
    return create_ccomplex_return(IT_0058*conj(IT_0058));
}

