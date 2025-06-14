#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_ss_R_to_s_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_ss_R_to_s_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = m_s*e_em*mu_h*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = A_b*V_ts*conj(V_ts);
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0004 + (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = pow(m_s, 2);
    const ccomplex_t IT_0013 = pow(m_sG, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_13 + IT_0012 + IT_0013 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = m_s*e_em*IT_0001*IT_0002*IT_0005*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = pow(m_ss_R, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + IT_0012 + -IT_0013 + -IT_0023 
      + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = m_sG*s_34;
    const ccomplex_t IT_0029 = pow(m_A0, 2);
    const ccomplex_t IT_0030 = s_13*IT_0029;
    const ccomplex_t IT_0031 = s_14*s_34;
    return create_ccomplex_return(8*s_13*IT_0017*conj(IT_0017) + 8*IT_0027*(
      (IT_0017 + conj(IT_0017))*IT_0028 + -IT_0027*(IT_0030 + (-2)*IT_0031)));
}

