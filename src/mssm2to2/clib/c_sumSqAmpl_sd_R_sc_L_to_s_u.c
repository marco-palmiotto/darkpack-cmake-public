#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_sc_L_to_s_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_sc_L_to_s_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0001 = pow(m_s, 2);
    const ccomplex_t IT_0002 = pow(m_sc_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp1)
      *IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = IT_0000*IT_0003*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*conj(V_cs)*e_em*conj(V_Wp2)
      *IT_0005;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0019 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0018*IT_0019;
    const ccomplex_t IT_0021 = IT_0013 + IT_0020;
    const ccomplex_t IT_0022 = 18*s_34;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*conj
      (V_cs)*e_em*IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0011*IT_0024;
    const ccomplex_t IT_0026 = IT_0003*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*conj
      (V_cs)*e_em*IT_0005*IT_0007*IT_0009;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0015*IT_0029;
    const ccomplex_t IT_0031 = IT_0019*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0027 + IT_0032;
    const ccomplex_t IT_0034 = m_u*IT_0033;
    const ccomplex_t IT_0035 = -IT_0027 + -IT_0032;
    const ccomplex_t IT_0036 = m_u*s_13;
    const ccomplex_t IT_0037 = m_s*m_u;
    const ccomplex_t IT_0038 = 18*IT_0021;
    const ccomplex_t IT_0039 = s_13*s_14;
    const ccomplex_t IT_0040 = pow(m_sd_R, 2);
    const ccomplex_t IT_0041 = s_34*IT_0040;
    const ccomplex_t IT_0042 = m_s*s_14;
    const ccomplex_t IT_0043 = 18*conj(IT_0021);
    const ccomplex_t IT_0044 = IT_0021*conj(IT_0021)*IT_0022 + conj(IT_0034)*
      (IT_0022*IT_0034 + 18*IT_0035*IT_0036 + IT_0037*IT_0038) + conj(IT_0035)*
      (18*IT_0034*IT_0036 + IT_0035*(36*IT_0039 + (-18)*IT_0041) + IT_0038
      *IT_0042) + (IT_0034*IT_0037 + IT_0035*IT_0042)*IT_0043;
    return create_ccomplex_return(IT_0044);
}

