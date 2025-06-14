#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_snu_mu_to_anti_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_snu_mu_to_anti_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t m_mu = param->m_mu;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*e_em*m_mu*mu_h
      *IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = pow(m_Hp, 2);
    const ccomplex_t IT_0008 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_23 + IT_0007 + IT_0008 + -m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = -IT_0011;
    const ccomplex_t IT_0014 = pow(m_smu_R, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + -IT_0008 + -IT_0014 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0013 + IT_0017;
    const ccomplex_t IT_0019 = IT_0001*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0018;
    const ccomplex_t IT_0021 = -s_24;
    const ccomplex_t IT_0022 = -s_34;
    const ccomplex_t IT_0023 = -IT_0012*(IT_0008*conj(IT_0012) + s_23*conj
      (IT_0019) + -conj(IT_0020)*IT_0021) + IT_0020*(conj(IT_0012)*IT_0021 +
       conj(IT_0019)*IT_0022) + -IT_0019*(s_23*conj(IT_0012) + IT_0007*conj
      (IT_0019) + -conj(IT_0020)*IT_0022);
    return create_ccomplex_return(IT_0023);
}

