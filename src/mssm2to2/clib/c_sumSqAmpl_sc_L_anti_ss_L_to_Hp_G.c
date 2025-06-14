#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_L_to_Hp_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_L_to_Hp_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = pow(m_s, 2);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = pow(m_c, 2);
    const ccomplex_t IT_0008 = cpow(IT_0006, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*(IT_0002*IT_0003*IT_0005*IT_0006 + (-2)*IT_0004*(m_W*IT_0006 + (-0.5
      )*IT_0002*IT_0007*IT_0008))*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = pow(m_Hp, 2);
    const ccomplex_t IT_0014 = pow(m_ss_L, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_sc_L, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0013 + IT_0018 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0012*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0017 + IT_0021;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = IT_0001*IT_0017;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (-4)*s_23;
    const ccomplex_t IT_0027 = -IT_0000;
    const ccomplex_t IT_0028 = IT_0022*IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = (-4)*s_24;
    const ccomplex_t IT_0031 = (-4)*s_34;
    const ccomplex_t IT_0032 = (-4)*IT_0023*(IT_0014*conj(IT_0023) + (-0.25)
      *conj(IT_0025)*IT_0026 + (-0.25)*conj(IT_0029)*IT_0030) + IT_0029*(conj
      (IT_0023)*IT_0030 + conj(IT_0025)*IT_0031) + (-4)*IT_0025*(IT_0013*conj
      (IT_0025) + (-0.25)*conj(IT_0023)*IT_0026 + (-0.25)*conj(IT_0029)*IT_0031);
    return create_ccomplex_return(IT_0032);
}

