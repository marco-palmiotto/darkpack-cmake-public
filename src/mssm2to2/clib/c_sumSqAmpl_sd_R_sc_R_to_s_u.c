#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_sc_R_to_s_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_sc_R_to_s_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_R = param->m_sc_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = sin(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u1)*e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0001*IT_0005*IT_0009;
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = pow(m_s, 2);
    const ccomplex_t IT_0014 = pow(m_sc_R, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0000*IT_0012*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0001*IT_0005*IT_0009;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *conj(V_u2)*e_em*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0023 = cpow((-2)*s_23 + IT_0013 + IT_0014 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0021*IT_0022*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    return create_ccomplex_return(s_34*(IT_0025*((-18)*conj(IT_0016) + 18*conj
      (IT_0025)) + IT_0016*(18*conj(IT_0016) + (-18)*conj(IT_0025))));
}

