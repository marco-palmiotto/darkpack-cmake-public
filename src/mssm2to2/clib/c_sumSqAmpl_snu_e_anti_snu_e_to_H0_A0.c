#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_e_anti_snu_e_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_e_anti_snu_e_to_H0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t m_snu_mu = param->m_snu_mu;
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
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0007;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = s_13 + s_23;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = pow(m_snu_mu, 2);
    const ccomplex_t IT_0015 = s_12 + IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0008;
    const ccomplex_t IT_0017 = -IT_0009;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = IT_0013 + IT_0019;
    const ccomplex_t IT_0021 = IT_0008*IT_0017;
    const ccomplex_t IT_0022 = s_23*IT_0021;
    const ccomplex_t IT_0023 = IT_0010*IT_0016;
    const ccomplex_t IT_0024 = IT_0015*IT_0023;
    const ccomplex_t IT_0025 = -IT_0022 + (-2)*IT_0024;
    const ccomplex_t IT_0026 = IT_0020 + IT_0025;
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0014 + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    return create_ccomplex_return(IT_0029*conj(IT_0029));
}

