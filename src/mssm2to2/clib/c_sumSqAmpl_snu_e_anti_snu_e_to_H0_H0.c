#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_anti_snu_e_to_H0_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_anti_snu_e_to_H0_H0(
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
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -2);
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = cos(alpha);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0000*IT_0005*(IT_0007 + 
      -IT_0009);
    const ccomplex_t IT_0011 = (-0.25)*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0003;
    const ccomplex_t IT_0013 = cpow(IT_0003, -1);
    const ccomplex_t IT_0014 = IT_0012 + IT_0013;
    const ccomplex_t IT_0015 = cos(beta);
    const ccomplex_t IT_0016 = cpow(IT_0008, 3);
    const ccomplex_t IT_0017 = sin(beta);
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_W*e_em*IT_0014*(IT_0015
      *IT_0016 + 5*IT_0006*(IT_0009*IT_0017 + (-1)*IT_0006*(IT_0008*IT_0015 + 1.
      /5*IT_0006*IT_0017)));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_W*e_em*IT_0014*(IT_0008
      *IT_0015 + IT_0006*IT_0017);
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*m_W*e_em*IT_0014*(IT_0016
      *IT_0017 + IT_0006*((IT_0007 + -IT_0009)*IT_0015 + -IT_0006*IT_0008
      *IT_0017));
    const ccomplex_t IT_0028 = (-1.5)*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_W*e_em*IT_0014*(IT_0006
      *IT_0015 + -IT_0008*IT_0017);
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = pow(m_H0, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0035 + reg_prop + (0 +
       _Complex_I*1)*m_snu_mu*Gamma_nuel, -1);
    const ccomplex_t IT_0037 = cpow(IT_0030, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0035 + reg_prop + (0 +
       _Complex_I*1)*m_snu_mu*Gamma_nuel, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = IT_0011 + -IT_0026 + -IT_0034 + IT_0039 +
       IT_0041;
    return create_ccomplex_return(IT_0042*conj(IT_0042));
}

