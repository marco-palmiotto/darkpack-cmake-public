#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_anti_snu_e_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_anti_snu_e_to_anti_Hp_W(
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
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0001 + 
      -IT_0002*IT_0003)*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = IT_0005 + IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0002 + 
      -IT_0001*IT_0003)*IT_0010;
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + (-2)*IT_0013 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0006*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*e_em*(IT_0000*IT_0002 +
       IT_0001*IT_0003)*IT_0005;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0001 +
       IT_0002*IT_0003)*IT_0010;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = -IT_0017 + -IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -2);
    const ccomplex_t IT_0028 = pow(s_14, 2);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, -1);
    const ccomplex_t IT_0033 = pow(m_e, 2);
    const ccomplex_t IT_0034 = cpow(IT_0002, -1);
    const ccomplex_t IT_0035 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0005*(m_W*IT_0002 + (-0.5)*IT_0032*IT_0033*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = pow(m_Hp, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0013 + IT_0037 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0031*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = s_14*s_34*IT_0027;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = s_13 + IT_0044;
    const ccomplex_t IT_0046 = -IT_0017 + -IT_0025 + IT_0041;
    const ccomplex_t IT_0047 = s_14*s_24*IT_0027;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = s_12 + IT_0048;
    const ccomplex_t IT_0050 = -IT_0026;
    const ccomplex_t IT_0051 = -conj(IT_0026);
    const ccomplex_t IT_0052 = pow(s_34, 2);
    const ccomplex_t IT_0053 = IT_0027*IT_0052;
    const ccomplex_t IT_0054 = s_24*s_34*IT_0027;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = s_23 + IT_0055;
    const ccomplex_t IT_0057 = pow(s_24, 2);
    const ccomplex_t IT_0058 = IT_0027*IT_0057;
    const ccomplex_t IT_0059 = -IT_0046;
    const ccomplex_t IT_0060 = (conj(IT_0026)*(IT_0013 + -IT_0029) + conj
      (IT_0042)*IT_0045 + conj(IT_0046)*IT_0049)*IT_0050 + IT_0046*IT_0049
      *IT_0051 + IT_0042*(IT_0045*IT_0051 + -conj(IT_0042)*(IT_0037 + -IT_0053) 
      + -conj(IT_0046)*IT_0056) + (conj(IT_0042)*IT_0056 + conj(IT_0046)*
      (IT_0013 + -IT_0058))*IT_0059;
    return create_ccomplex_return(IT_0060);
}

