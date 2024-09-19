#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_R_to_H0_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_R_to_H0_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_su_L = param->m_su_L;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(alpha);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*m_u*e_em*mu_h*IT_0002
      *IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = pow(m_H0, 2);
    const ccomplex_t IT_0011 = pow(m_su_R, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0010 + IT_0011 + -m_su_L*
      (m_su_L + (0 + _Complex_I*-1)*Gamma_ul) + reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = pow(m_su_L, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0010 + IT_0015 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0014 + IT_0018;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = IT_0001*IT_0014;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = (-4)*s_23;
    const ccomplex_t IT_0024 = -IT_0000;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = (-4)*s_24;
    const ccomplex_t IT_0028 = (-4)*s_34;
    const ccomplex_t IT_0029 = (-4)*IT_0020*(IT_0011*conj(IT_0020) + (-0.25)
      *conj(IT_0022)*IT_0023 + (-0.25)*conj(IT_0026)*IT_0027) + IT_0026*(conj
      (IT_0020)*IT_0027 + conj(IT_0022)*IT_0028) + (-4)*IT_0022*(IT_0010*conj
      (IT_0022) + (-0.25)*conj(IT_0020)*IT_0023 + (-0.25)*conj(IT_0026)*IT_0028);
    return create_ccomplex_return(IT_0029);
}

