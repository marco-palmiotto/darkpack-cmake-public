#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_R_anti_se_R_to_h_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_R_anti_se_R_to_h_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_se_R = param->m_se_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_er = param->Gamma_er;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -2);
    const ccomplex_t IT_0004 = cos(alpha);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = pow(m_e, 2);
    const ccomplex_t IT_0011 = cpow(IT_0006, -1);
    const ccomplex_t IT_0012 = cpow(IT_0008, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(m_W*IT_0003*(IT_0004
      *IT_0005 + IT_0006*IT_0007)*IT_0008 + -IT_0007*IT_0009*IT_0010*IT_0011
      *IT_0012);
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = pow(m_h, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0015 + (0 + _Complex_I*1)
      *m_se_R*Gamma_er + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0015 + (0 + _Complex_I*1)
      *m_se_R*Gamma_er + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0018 + IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = pow(m_se_R, 2);
    const ccomplex_t IT_0025 = IT_0001*IT_0021;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = -s_23;
    const ccomplex_t IT_0028 = IT_0000*IT_0022;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = -s_24;
    const ccomplex_t IT_0031 = -s_34;
    const ccomplex_t IT_0032 = -IT_0023*(conj(IT_0023)*IT_0024 + -conj(IT_0026
      )*IT_0027 + -conj(IT_0029)*IT_0030) + IT_0029*(conj(IT_0023)*IT_0030 +
       conj(IT_0026)*IT_0031) + -IT_0026*(IT_0015*conj(IT_0026) + -conj(IT_0023)
      *IT_0027 + -conj(IT_0029)*IT_0031);
    return create_ccomplex_return(IT_0032);
}

