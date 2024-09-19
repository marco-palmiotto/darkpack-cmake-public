#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_smu_R_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_smu_R_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0000, -2);
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = pow(m_mu, 2);
    const ccomplex_t IT_0012 = cpow(IT_0004, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*(IT_0002*IT_0004 +
       -IT_0003*IT_0005)*IT_0006*IT_0009 + -IT_0002*IT_0007*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0014 = pow(m_smu_R, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + (-2)*IT_0014 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0008*IT_0017;
    const ccomplex_t IT_0019 = e_em*IT_0001*(IT_0002*IT_0004 + IT_0003*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(m_W*(IT_0002*IT_0003 +
       IT_0004*IT_0005)*IT_0006*IT_0009 + -IT_0005*IT_0007*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0018 + -IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, -2);
    const ccomplex_t IT_0029 = pow(s_14, 2);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = s_14*s_24*IT_0028;
    const ccomplex_t IT_0032 = pow(s_24, 2);
    const ccomplex_t IT_0033 = IT_0028*IT_0032;
    return create_ccomplex_return((-2)*IT_0027*conj(IT_0027)*(s_12 + IT_0014 +
       (-0.5)*IT_0030 + -IT_0031 + (-0.5)*IT_0033));
}

