#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sd_R_to_Hp_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sd_R_to_Hp_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = cos(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_d*m_u*V_ud
      *e_em*IT_0002*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = pow(m_Hp, 2);
    const ccomplex_t IT_0012 = pow(m_sd_R, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_su_R*
      (m_su_R + (0 + _Complex_I*-1)*Gamma_ur) + reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = pow(m_su_R, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0011 + IT_0016 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0015 + IT_0019;
    const ccomplex_t IT_0021 = IT_0001*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0015;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = (-4)*s_23;
    const ccomplex_t IT_0025 = -IT_0000;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = (-4)*s_24;
    const ccomplex_t IT_0029 = (-4)*s_34;
    const ccomplex_t IT_0030 = (-4)*IT_0021*(IT_0012*conj(IT_0021) + (-0.25)
      *conj(IT_0023)*IT_0024 + (-0.25)*conj(IT_0027)*IT_0028) + IT_0027*(conj
      (IT_0021)*IT_0028 + conj(IT_0023)*IT_0029) + (-4)*IT_0023*(IT_0011*conj
      (IT_0023) + (-0.25)*conj(IT_0021)*IT_0024 + (-0.25)*conj(IT_0027)*IT_0029);
    return create_ccomplex_return(IT_0030);
}

