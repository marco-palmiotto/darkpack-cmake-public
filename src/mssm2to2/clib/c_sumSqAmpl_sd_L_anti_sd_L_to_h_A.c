#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_L_to_h_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_L_to_h_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = pow(m_d, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0003
      *IT_0005*IT_0006*IT_0008 + (-0.5)*m_W*(IT_0004*IT_0006 + IT_0009*IT_0010)*
      (IT_0008 + 1./3*IT_0007*IT_0012));
    const ccomplex_t IT_0014 = pow(m_h, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0014 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0014 + (0 + _Complex_I*1)
      *m_sd_L*Gamma_dl + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0013*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017 + IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = pow(m_sd_L, 2);
    const ccomplex_t IT_0024 = IT_0001*IT_0017;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (-3)*s_23;
    const ccomplex_t IT_0027 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0028 = IT_0021*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-3)*s_24;
    const ccomplex_t IT_0031 = (-3)*s_34;
    const ccomplex_t IT_0032 = (-3)*IT_0022*(conj(IT_0022)*IT_0023 + (
      -0.333333333333333)*conj(IT_0025)*IT_0026 + (-0.333333333333333)*conj
      (IT_0029)*IT_0030) + IT_0029*(conj(IT_0022)*IT_0030 + conj(IT_0025)
      *IT_0031) + (-3)*IT_0025*(IT_0014*conj(IT_0025) + (-0.333333333333333)
      *conj(IT_0022)*IT_0026 + (-0.333333333333333)*conj(IT_0029)*IT_0031);
    return create_ccomplex_return(IT_0032);
}

