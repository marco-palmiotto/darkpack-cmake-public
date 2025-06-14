#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_snu_e_to_A0_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_snu_e_to_A0_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_e*e_em*mu_h*IT_0001*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_A0, 2);
    const ccomplex_t IT_0008 = pow(m_se_R, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = e_em*IT_0001;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0001*IT_0004;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + -IT_0008 + -IT_0016 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0013*IT_0019;
    const ccomplex_t IT_0021 = IT_0012 + -IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -2);
    const ccomplex_t IT_0023 = pow(s_24, 2);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = -IT_0020;
    const ccomplex_t IT_0026 = s_14*s_24*IT_0022;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = pow(s_14, 2);
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    return create_ccomplex_return(-IT_0021*(conj(IT_0021)*(IT_0016 + -IT_0024)
       + conj(IT_0025)*(s_12 + IT_0027)) + -IT_0025*(conj(IT_0021)*(s_12 +
       IT_0027) + conj(IT_0025)*(IT_0008 + -IT_0029)));
}

