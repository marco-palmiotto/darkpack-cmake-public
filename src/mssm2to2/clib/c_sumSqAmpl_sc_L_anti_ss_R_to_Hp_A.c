#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_ss_R_to_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_ss_R_to_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*mu_h*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = pow(m_Hp, 2);
    const ccomplex_t IT_0008 = pow(m_sc_L, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_13 + IT_0007 + IT_0008 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = pow(m_ss_R, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0008 + -IT_0013 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0006*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0007 + IT_0013 + -m_sc_L*
      (m_sc_L + (0 + _Complex_I*-1)*Gamma_cl) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0006*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = -IT_0012 + IT_0017 + -IT_0022;
    const ccomplex_t IT_0024 = 2*IT_0000;
    const ccomplex_t IT_0025 = IT_0016*IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = (-3)*s_14;
    const ccomplex_t IT_0028 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0031 = IT_0011*IT_0030;
    const ccomplex_t IT_0032 = -IT_0025 + IT_0029 + IT_0031;
    const ccomplex_t IT_0033 = (-3)*s_24;
    const ccomplex_t IT_0034 = -IT_0029;
    const ccomplex_t IT_0035 = (-3)*s_34;
    const ccomplex_t IT_0036 = (-3)*s_13;
    const ccomplex_t IT_0037 = (-3)*s_23;
    const ccomplex_t IT_0038 = (-0.333333333333333)*conj(IT_0034);
    const ccomplex_t IT_0039 = IT_0023*(conj(IT_0026)*IT_0027 + conj(IT_0032)
      *IT_0033 + conj(IT_0034)*IT_0035) + (-3)*IT_0034*(IT_0007*conj(IT_0034) + 
      (-0.333333333333333)*conj(IT_0023)*IT_0035 + (-0.333333333333333)*conj
      (IT_0026)*IT_0036 + (-0.333333333333333)*conj(IT_0032)*IT_0037) + (-3)
      *IT_0026*(IT_0008*conj(IT_0026) + (-0.333333333333333)*conj(IT_0023)
      *IT_0027 + s_12*conj(IT_0032) + IT_0036*IT_0038) + (-3)*IT_0032*(s_12*conj
      (IT_0026) + IT_0013*conj(IT_0032) + (-0.333333333333333)*conj(IT_0023)
      *IT_0033 + IT_0037*IT_0038);
    return create_ccomplex_return(IT_0039);
}

