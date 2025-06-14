#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_snu_e_to_anti_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_snu_e_to_anti_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = pow(m_e, 2);
    const ccomplex_t IT_0005 = cpow(IT_0002, -1);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*e_em*(m_W
      *IT_0002 + (-0.5)*IT_0003*IT_0004*IT_0005)*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = pow(m_Hp, 2);
    const ccomplex_t IT_0012 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = -IT_0015;
    const ccomplex_t IT_0018 = pow(m_se_L, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + -IT_0012 + -IT_0018 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0010*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0017 + IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = IT_0000*IT_0022;
    const ccomplex_t IT_0025 = -s_24;
    const ccomplex_t IT_0026 = -s_34;
    const ccomplex_t IT_0027 = -IT_0016*(IT_0012*conj(IT_0016) + s_23*conj
      (IT_0023) + -conj(IT_0024)*IT_0025) + IT_0024*(conj(IT_0016)*IT_0025 +
       conj(IT_0023)*IT_0026) + -IT_0023*(s_23*conj(IT_0016) + IT_0011*conj
      (IT_0023) + -conj(IT_0024)*IT_0026);
    return create_ccomplex_return(IT_0027);
}

