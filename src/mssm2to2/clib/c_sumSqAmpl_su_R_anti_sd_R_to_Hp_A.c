#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_sd_R_to_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_sd_R_to_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = (-0.666666666666667)*IT_0000;
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
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = pow(m_su_R, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + -IT_0012 + -IT_0017 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0010*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0000*IT_0020;
    const ccomplex_t IT_0022 = 0.333333333333333*IT_0000;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0011 + IT_0017 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0010*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = -IT_0016 + IT_0021 + -IT_0026;
    const ccomplex_t IT_0028 = 2*IT_0000;
    const ccomplex_t IT_0029 = IT_0020*IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (-3)*s_14;
    const ccomplex_t IT_0032 = 0.666666666666667*IT_0000;
    const ccomplex_t IT_0033 = IT_0025*IT_0032;
    const ccomplex_t IT_0034 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0035 = IT_0015*IT_0034;
    const ccomplex_t IT_0036 = -IT_0029 + IT_0033 + IT_0035;
    const ccomplex_t IT_0037 = (-3)*s_24;
    const ccomplex_t IT_0038 = -IT_0035;
    const ccomplex_t IT_0039 = (-3)*s_34;
    const ccomplex_t IT_0040 = (-3)*s_13;
    const ccomplex_t IT_0041 = (-3)*s_23;
    const ccomplex_t IT_0042 = (-0.333333333333333)*conj(IT_0038);
    const ccomplex_t IT_0043 = IT_0027*(conj(IT_0030)*IT_0031 + conj(IT_0036)
      *IT_0037 + conj(IT_0038)*IT_0039) + (-3)*IT_0038*(IT_0011*conj(IT_0038) + 
      (-0.333333333333333)*conj(IT_0027)*IT_0039 + (-0.333333333333333)*conj
      (IT_0030)*IT_0040 + (-0.333333333333333)*conj(IT_0036)*IT_0041) + (-3)
      *IT_0030*(IT_0017*conj(IT_0030) + (-0.333333333333333)*conj(IT_0027)
      *IT_0031 + s_12*conj(IT_0036) + IT_0040*IT_0042) + (-3)*IT_0036*(s_12*conj
      (IT_0030) + IT_0012*conj(IT_0036) + (-0.333333333333333)*conj(IT_0027)
      *IT_0037 + IT_0041*IT_0042);
    return create_ccomplex_return(IT_0043);
}

