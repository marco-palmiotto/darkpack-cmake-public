#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_mu_sd_R_to_mu_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_mu_sd_R_to_mu_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em*m_mu
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_mu, 2);
    const ccomplex_t IT_0011 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em*m_mu
      *IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = IT_0015*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0014 + IT_0022;
    const ccomplex_t IT_0024 = s_13*s_14;
    const ccomplex_t IT_0025 = s_34*IT_0011;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0004;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0029 = IT_0012*IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0004;
    const ccomplex_t IT_0032 = IT_0019*IT_0031;
    const ccomplex_t IT_0033 = IT_0015*IT_0030*IT_0032;
    const ccomplex_t IT_0034 = IT_0029 + IT_0033;
    const ccomplex_t IT_0035 = m_mu*IT_0023;
    const ccomplex_t IT_0036 = m_mu*s_14;
    const ccomplex_t IT_0037 = (-6)*s_34;
    const ccomplex_t IT_0038 = IT_0023*conj(IT_0023)*(12*IT_0024 + (-6)
      *IT_0025) + 6*s_34*(IT_0034*conj(IT_0034) + IT_0035*conj(IT_0035)) + 6*
      (conj(IT_0023)*(IT_0034 + -IT_0035) + IT_0023*(conj(IT_0034) + -conj
      (IT_0035)))*IT_0036 + (conj(IT_0034)*IT_0035 + IT_0034*conj(IT_0035))
      *IT_0037;
    return create_ccomplex_return(IT_0038);
}

