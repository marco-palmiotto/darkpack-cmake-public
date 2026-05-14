#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_L_anti_snu_mu_to_anti_nu_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_L_anti_snu_mu_to_anti_nu_tau_tau(
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
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = pow(m_mu, 2);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0003
      *IT_0005*(m_W*IT_0001 + (-0.5)*IT_0000*IT_0002*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = pow(m_smu_L, 2);
    const ccomplex_t IT_0013 = pow(m_snu_mu, 2);
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
    const ccomplex_t IT_0022 = (-0.5)*IT_0017;
    const ccomplex_t IT_0023 = m_tau*IT_0022;
    const ccomplex_t IT_0024 = -IT_0017;
    const ccomplex_t IT_0025 = m_tau*IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0023 + IT_0026;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0000*IT_0005;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = cpow(IT_0001, 2);
    const ccomplex_t IT_0032 = cpow(IT_0003, 2);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005*
      (IT_0000*IT_0008 + -m_W*(IT_0031 + -IT_0032));
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0030*IT_0034;
    const ccomplex_t IT_0036 = pow(m_W, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0036 
      + -reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = -IT_0016 + -IT_0028 + -IT_0039;
    const ccomplex_t IT_0041 = IT_0021*IT_0024;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = s_13*m_tau;
    const ccomplex_t IT_0044 = s_13*s_14;
    const ccomplex_t IT_0045 = s_34*IT_0012;
    return create_ccomplex_return(2*IT_0040*(s_34*conj(IT_0040) + conj(IT_0042
      )*IT_0043) + 2*IT_0042*(conj(IT_0040)*IT_0043 + 2*conj(IT_0042)*(IT_0044 +
       (-0.5)*IT_0045)));
}

