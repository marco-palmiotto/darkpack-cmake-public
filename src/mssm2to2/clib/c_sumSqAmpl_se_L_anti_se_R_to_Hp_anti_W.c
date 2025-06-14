#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_R_to_Hp_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_R_to_Hp_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = e_em*IT_0001;
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = m_e*e_em*mu_h*IT_0001*IT_0003;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = pow(m_se_L, 2);
    const ccomplex_t IT_0007 = pow(m_se_R, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0002*IT_0010;
    const ccomplex_t IT_0012 = cos(alpha);
    const ccomplex_t IT_0013 = sin(beta);
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = sin(alpha);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0012
      *IT_0013 + -IT_0014*IT_0015);
    const ccomplex_t IT_0017 = cpow(IT_0014, -1);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0001
      *IT_0003*IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + -IT_0006 + -IT_0007 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0012
      *IT_0014 + IT_0013*IT_0015);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0001
      *IT_0003*IT_0012*IT_0017;
    const ccomplex_t IT_0027 = 0.5*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0006 + -IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = IT_0025*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0001*IT_0003;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = pow(m_Hp, 2);
    const ccomplex_t IT_0037 = cpow((-2)*s_23 + IT_0007 + IT_0036 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0038 = IT_0035*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = IT_0033*IT_0039;
    const ccomplex_t IT_0041 = IT_0011 + IT_0023 + IT_0031 + -IT_0040;
    const ccomplex_t IT_0042 = pow(m_W, -2);
    const ccomplex_t IT_0043 = pow(s_34, 2);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = s_24*s_34*IT_0042;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = s_23 + IT_0046;
    const ccomplex_t IT_0048 = pow(s_24, 2);
    const ccomplex_t IT_0049 = IT_0042*IT_0048;
    return create_ccomplex_return(-IT_0041*(conj(IT_0041)*(IT_0036 + -IT_0044)
       + conj(IT_0040)*IT_0047) + -IT_0040*(conj(IT_0041)*IT_0047 + conj(IT_0040
      )*(IT_0007 + -IT_0049)));
}

