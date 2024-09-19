#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_snu_tau_to_h_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_snu_tau_to_h_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = sin(alpha);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -2);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0000*IT_0001*IT_0003
      *IT_0004*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0005;
    const ccomplex_t IT_0010 = cpow(IT_0005, -1);
    const ccomplex_t IT_0011 = IT_0009 + IT_0010;
    const ccomplex_t IT_0012 = sin(beta);
    const ccomplex_t IT_0013 = cpow(IT_0004, 3);
    const ccomplex_t IT_0014 = cos(beta);
    const ccomplex_t IT_0015 = cpow(IT_0001, 2);
    const ccomplex_t IT_0016 = cpow(IT_0004, 2);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0012
      *IT_0013 + (-5)*IT_0001*(IT_0001*IT_0004*IT_0012 + (-0.2)*IT_0014*(IT_0015
       + (-5)*IT_0016)));
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0001
      *IT_0012 + IT_0004*IT_0014);
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0022 + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = pow(m_h, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0026 + reg_prop + (0 +
       _Complex_I*1)*m_snu_tau*Gamma_nutaul, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*-1)*m_W*e_em*IT_0011*(IT_0004
      *IT_0012 + -IT_0001*IT_0014);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0020*IT_0029;
    const ccomplex_t IT_0031 = IT_0027*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0026 + reg_prop + (0 +
       _Complex_I*1)*m_snu_tau*Gamma_nutaul, -1);
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_W*e_em*IT_0011*(IT_0013
      *IT_0014 + 5*IT_0001*((-0.2)*IT_0001*(IT_0001*IT_0012 + 5*IT_0004*IT_0014)
       + IT_0012*IT_0016));
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0029*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + (-2)*IT_0022 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0008 + -IT_0025 + IT_0032 + IT_0035 + 
      -IT_0041;
    return create_ccomplex_return(IT_0042*conj(IT_0042));
}

