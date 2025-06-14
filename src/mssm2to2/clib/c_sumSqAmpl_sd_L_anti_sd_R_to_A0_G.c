#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_L_anti_sd_R_to_A0_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_L_anti_sd_R_to_A0_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_td = param->V_td;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_d*e_em*mu_h*IT_0004;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = A_b*V_td*conj(V_td);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0006 + 0.5*IT_0010;
    const ccomplex_t IT_0012 = pow(m_A0, 2);
    const ccomplex_t IT_0013 = pow(m_sd_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_23 + IT_0012 + IT_0013 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_sd_L, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0012 + IT_0017 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0016 + IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0016;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = (-4)*s_23;
    const ccomplex_t IT_0026 = -IT_0000;
    const ccomplex_t IT_0027 = IT_0021*IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = (-4)*s_24;
    const ccomplex_t IT_0030 = (-4)*s_34;
    const ccomplex_t IT_0031 = (-4)*IT_0022*(IT_0013*conj(IT_0022) + (-0.25)
      *conj(IT_0024)*IT_0025 + (-0.25)*conj(IT_0028)*IT_0029) + IT_0028*(conj
      (IT_0022)*IT_0029 + conj(IT_0024)*IT_0030) + (-4)*IT_0024*(IT_0012*conj
      (IT_0024) + (-0.25)*conj(IT_0022)*IT_0025 + (-0.25)*conj(IT_0028)*IT_0030);
    return create_ccomplex_return(IT_0031);
}

