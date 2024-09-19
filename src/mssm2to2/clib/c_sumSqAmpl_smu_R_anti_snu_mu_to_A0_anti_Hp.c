#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_snu_mu_to_A0_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_snu_mu_to_A0_anti_Hp(
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
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = e_em*m_mu*mu_h*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-0.5)*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = pow(m_mu, 2);
    const ccomplex_t IT_0007 = cpow(IT_0005, -1);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0002*
      (m_W*IT_0005 + (-0.5)*IT_0000*IT_0006*IT_0007)*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = pow(m_A0, 2);
    const ccomplex_t IT_0013 = pow(m_smu_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = m_W*e_em*IT_0002;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0000*IT_0002*IT_0007*IT_0008;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, 2);
    const ccomplex_t IT_0023 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + -IT_0013 + IT_0022 + -IT_0023 
      + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0016 + -IT_0026;
    return create_ccomplex_return(IT_0027*conj(IT_0027));
}

