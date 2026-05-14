#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_sc_R_to_anti_s_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_sc_R_to_anti_s_G(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_cr = param->Gamma_cr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_c*V_cs*V_u1
      *e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = pow(m_s, 2);
    const ccomplex_t IT_0010 = pow(m_C_1, 2);
    const ccomplex_t IT_0011 = pow(m_sc_R, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + IT_0009 + -IT_0010 + -IT_0011 
      + -reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = s_14*s_34;
    const ccomplex_t IT_0016 = s_14*IT_0009;
    const ccomplex_t IT_0017 = (-2)*IT_0000;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0009 + IT_0010 + -m_sc_R*
      (m_sc_R + (0 + _Complex_I*-1)*Gamma_cr) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = 2*IT_0014 + IT_0021;
    const ccomplex_t IT_0023 = -IT_0021;
    const ccomplex_t IT_0024 = s_34*IT_0010;
    const ccomplex_t IT_0025 = s_13*s_14;
    const ccomplex_t IT_0026 = (-16)*IT_0025;
    const ccomplex_t IT_0027 = 8*IT_0024 + IT_0026;
    const ccomplex_t IT_0028 = -IT_0000;
    const ccomplex_t IT_0029 = -IT_0020*(IT_0017 + -IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (-8)*s_13;
    const ccomplex_t IT_0032 = s_34*IT_0031;
    const ccomplex_t IT_0033 = pow(s_13, 2);
    const ccomplex_t IT_0034 = (-8)*IT_0033;
    const ccomplex_t IT_0035 = IT_0014*(32*conj(IT_0014)*IT_0015 + (-8)
      *IT_0016*conj(IT_0022) + conj(IT_0023)*IT_0027) + IT_0030*(s_14*conj
      (IT_0023)*IT_0031 + conj(IT_0022)*IT_0032) + IT_0023*(conj(IT_0014)
      *IT_0027 + (IT_0010*conj(IT_0023) + s_14*conj(IT_0030))*IT_0031 + conj
      (IT_0022)*IT_0034) + (-8)*IT_0022*(conj(IT_0014)*IT_0016 + s_13*IT_0009
      *conj(IT_0022) + (-0.125)*conj(IT_0030)*IT_0032 + (-0.125)*conj(IT_0023)
      *IT_0034);
    return create_ccomplex_return(IT_0035);
}

