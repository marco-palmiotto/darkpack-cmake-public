#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_snu_e_to_H0_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_snu_e_to_H0_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_se_L = param->m_se_L;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 3);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -2);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0002
      *IT_0004 + -IT_0005*(IT_0004*IT_0006*IT_0007 + 2*IT_0001*(IT_0000*IT_0005 
      + -IT_0001*IT_0007)*(0.5*IT_0004 + IT_0003*IT_0009)));
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cpow(IT_0005, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0001*IT_0004*IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, 2);
    const ccomplex_t IT_0018 = pow(m_se_R, 2);
    const ccomplex_t IT_0019 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + IT_0017 + -IT_0018 + -IT_0019 
      + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_e*e_em*mu_h
      *IT_0004*IT_0012;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = -IT_0004;
    const ccomplex_t IT_0026 = IT_0003*IT_0009;
    const ccomplex_t IT_0027 = IT_0025 + IT_0026;
    const ccomplex_t IT_0028 = -IT_0005;
    const ccomplex_t IT_0029 = 3*IT_0004 + IT_0026;
    const ccomplex_t IT_0030 = cpow(IT_0001, 2);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_W*e_em*(IT_0002*IT_0007
      *IT_0027 + IT_0028*(IT_0028*(IT_0000*IT_0005*IT_0027 + -IT_0001*IT_0007
      *IT_0029) + IT_0000*IT_0029*IT_0030));
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0024*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_12 + -IT_0018 + -IT_0019 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*m_e*e_em*mu_h*IT_0004
      *IT_0007*IT_0012*IT_0013;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = pow(m_e, 2);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *IT_0004*(m_W*IT_0005 + (-0.5)*IT_0012*IT_0013*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = pow(m_H0, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_13 + IT_0018 + IT_0043 + -m_se_L*
      (m_se_L + (0 + _Complex_I*-1)*Gamma_el) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0000
      *IT_0005 + -IT_0001*IT_0007)*IT_0009 + -IT_0000*IT_0004*IT_0012*IT_0013
      *IT_0039);
    const ccomplex_t IT_0048 = IT_0024*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_13 + IT_0043 + (0 + _Complex_I*1)
      *m_se_R*Gamma_er + reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_23 + IT_0043 + reg_prop + (0 +
       _Complex_I*1)*m_snu_mu*Gamma_nuel, -1);
    const ccomplex_t IT_0053 = IT_0004 + IT_0026;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_W*e_em*(IT_0000*IT_0005 + 
      -IT_0001*IT_0007)*IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = IT_0024*IT_0055;
    const ccomplex_t IT_0057 = IT_0052*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0022 + -IT_0036 + IT_0046 + IT_0051 +
       IT_0058;
    return create_ccomplex_return(IT_0059*conj(IT_0059));
}

