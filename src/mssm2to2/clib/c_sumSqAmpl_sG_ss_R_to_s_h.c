#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_ss_R_to_s_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_ss_R_to_s_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_b = param->A_b;
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t m_s = param->m_s;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_s*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = pow(m_s, 2);
    const ccomplex_t IT_0011 = pow(m_sG, 2);
    const ccomplex_t IT_0012 = pow(m_ss_R, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + IT_0010 + -IT_0011 + -IT_0012 
      + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = cos(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -2);
    const ccomplex_t IT_0019 = cos(alpha);
    const ccomplex_t IT_0020 = sin(beta);
    const ccomplex_t IT_0021 = (0 + _Complex_I*-3)*e_em*(IT_0000*IT_0002
      *IT_0003*IT_0005*IT_0010 + -1./3*m_W*IT_0004*IT_0018*(IT_0001*IT_0003 +
       IT_0019*IT_0020));
    const ccomplex_t IT_0022 = (-0.333333333333333)*IT_0021;
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_ss_R*
      (m_ss_R + (0 + _Complex_I*-1)*Gamma_sr) + reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = m_s*IT_0015;
    const ccomplex_t IT_0028 = IT_0026 + (-2)*IT_0027;
    const ccomplex_t IT_0029 = m_s*s_14;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(1.4142135623731*A_b*V_ts
      *conj(V_ts)*IT_0003 + m_s*e_em*mu_h*IT_0000*IT_0002*IT_0005*IT_0019);
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0008*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0010 + IT_0011 + -m_ss_L*
      (m_ss_L + (0 + _Complex_I*-1)*Gamma_sl) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = m_sG*s_34;
    const ccomplex_t IT_0038 = pow(m_h, 2);
    const ccomplex_t IT_0039 = s_13*IT_0038;
    const ccomplex_t IT_0040 = s_14*s_34;
    const ccomplex_t IT_0041 = m_s*m_sG;
    const ccomplex_t IT_0042 = 8*s_13;
    const ccomplex_t IT_0043 = (-8)*IT_0016*((IT_0028 + -conj(IT_0028))
      *IT_0029 + (IT_0036 + -conj(IT_0036))*IT_0037 + -IT_0016*(IT_0039 + (-2)
      *IT_0040)) + 8*conj(IT_0028)*(IT_0036*IT_0041 + 0.125*IT_0028*IT_0042) + 8
      *conj(IT_0036)*(IT_0028*IT_0041 + 0.125*IT_0036*IT_0042);
    return create_ccomplex_return(IT_0043);
}

