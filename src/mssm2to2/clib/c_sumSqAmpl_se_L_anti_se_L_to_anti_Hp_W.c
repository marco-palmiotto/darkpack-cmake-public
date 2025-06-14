#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_anti_se_L_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_anti_se_L_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = pow(m_e, 2);
    const ccomplex_t IT_0009 = cpow(IT_0002, -1);
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0005
      *IT_0007*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0000*IT_0002 + -IT_0001*IT_0003)
      *(IT_0005 + -IT_0004*IT_0011));
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_se_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + (-2)*IT_0014 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0003*IT_0005
      *IT_0007*IT_0008*IT_0009 + (-0.5)*m_W*(IT_0000*IT_0001 + IT_0002*IT_0003)*
      (IT_0005 + -IT_0004*IT_0011));
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0020*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0005*(m_W*IT_0002 + (-0.5)*IT_0007*IT_0008*IT_0009);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = pow(m_Hp, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_13 + IT_0014 + IT_0030 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0033;
    const ccomplex_t IT_0035 = -IT_0018 + -IT_0025 + IT_0034;
    const ccomplex_t IT_0036 = pow(m_W, -2);
    const ccomplex_t IT_0037 = pow(s_24, 2);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = -IT_0018 + -IT_0025;
    const ccomplex_t IT_0040 = s_14*s_24*IT_0036;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = pow(s_14, 2);
    const ccomplex_t IT_0043 = IT_0036*IT_0042;
    return create_ccomplex_return(-IT_0035*(conj(IT_0035)*(IT_0014 + -IT_0038)
       + conj(IT_0039)*(s_12 + IT_0041)) + -IT_0039*(conj(IT_0035)*(s_12 +
       IT_0041) + conj(IT_0039)*(IT_0014 + -IT_0043)));
}

