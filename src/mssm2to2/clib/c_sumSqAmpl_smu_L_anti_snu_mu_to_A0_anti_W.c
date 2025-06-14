#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_snu_mu_to_A0_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_snu_mu_to_A0_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = e_em*IT_0001;
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = pow(m_mu, 2);
    const ccomplex_t IT_0006 = cpow(IT_0003, -1);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001*
      (m_W*IT_0003 + (-0.5)*IT_0004*IT_0005*IT_0006)*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = pow(m_smu_L, 2);
    const ccomplex_t IT_0011 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = IT_0002*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, -2);
    const ccomplex_t IT_0018 = pow(s_14, 2);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = s_14*s_24*IT_0017;
    const ccomplex_t IT_0021 = pow(s_24, 2);
    const ccomplex_t IT_0022 = IT_0017*IT_0021;
    return create_ccomplex_return((-2)*IT_0016*conj(IT_0016)*(s_12 + 0.5
      *IT_0010 + 0.5*IT_0011 + (-0.5)*IT_0019 + -IT_0020 + (-0.5)*IT_0022));
}

