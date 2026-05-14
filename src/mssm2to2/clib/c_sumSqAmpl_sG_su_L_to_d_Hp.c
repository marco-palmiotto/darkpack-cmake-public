#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_su_L_to_d_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_su_L_to_d_Hp(
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
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_sd_L = param->m_sd_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_dl = param->Gamma_dl;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = pow(m_su_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0005 = pow(m_W, -1);
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *IT_0005*IT_0007*IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_d*V_ud*e_em
      *mu_h*IT_0005*IT_0010;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0004*IT_0017;
    const ccomplex_t IT_0019 = pow(m_d, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0001 + IT_0019 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = m_sG*IT_0015;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = m_sG*s_23;
    const ccomplex_t IT_0026 = cpow(IT_0008, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0010*(IT_0005*IT_0007*IT_0008*IT_0019 + (-2)*IT_0006*(m_W*IT_0008 + (
      -0.5)*IT_0000*IT_0005*IT_0026));
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_13 + IT_0001 + IT_0019 + -m_sd_L*
      (m_sd_L + (0 + _Complex_I*-1)*Gamma_dl) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_u*V_ud*e_em
      *IT_0005*IT_0006*IT_0010*IT_0026;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0037 = IT_0003*IT_0035*IT_0036;
    const ccomplex_t IT_0038 = -IT_0032 + IT_0037;
    const ccomplex_t IT_0039 = m_d*s_12;
    const ccomplex_t IT_0040 = s_12*s_23;
    const ccomplex_t IT_0041 = s_13*IT_0002;
    const ccomplex_t IT_0042 = m_d*m_sG;
    const ccomplex_t IT_0043 = 8*s_13;
    const ccomplex_t IT_0044 = (-8)*IT_0015*((IT_0024 + -conj(IT_0024))
      *IT_0025 + (IT_0038 + -conj(IT_0038))*IT_0039 + 0.125*IT_0015*(16*IT_0040 
      + (-8)*IT_0041)) + 8*conj(IT_0024)*(IT_0038*IT_0042 + 0.125*IT_0024
      *IT_0043) + 8*conj(IT_0038)*(IT_0024*IT_0042 + 0.125*IT_0038*IT_0043);
    return create_ccomplex_return(IT_0044);
}

