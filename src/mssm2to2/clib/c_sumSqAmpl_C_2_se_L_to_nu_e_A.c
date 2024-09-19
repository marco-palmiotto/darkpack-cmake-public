#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_se_L_to_nu_e_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_se_L_to_nu_e_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_el = param->Gamma_el;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = pow(m_se_L, 2);
    const ccomplex_t IT_0006 = cpow(s_23 + (-0.5)*IT_0005 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = s_14*s_34;
    const ccomplex_t IT_0011 = 2*IT_0003;
    const ccomplex_t IT_0012 = pow(m_C_2, 2);
    const ccomplex_t IT_0013 = cpow(s_13 + (-0.5)*IT_0012 + 0.5*m_se_L*(m_se_L
       + (0 + _Complex_I*-1)*Gamma_el) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0002*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0011*IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = s_13*s_34;
    const ccomplex_t IT_0019 = IT_0003*IT_0015;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = (-2)*s_13;
    const ccomplex_t IT_0022 = s_34*IT_0021;
    const ccomplex_t IT_0023 = IT_0008 + 0.5*IT_0016;
    const ccomplex_t IT_0024 = s_34*IT_0012;
    const ccomplex_t IT_0025 = pow(s_13, 2);
    const ccomplex_t IT_0026 = (-2)*IT_0025;
    const ccomplex_t IT_0027 = 8*IT_0009*conj(IT_0009)*IT_0010 + (-4)*conj
      (IT_0017)*(IT_0009*IT_0018 + (-0.25)*IT_0020*IT_0022) + (-4)*IT_0017*(conj
      (IT_0009)*IT_0018 + (-0.25)*conj(IT_0020)*IT_0022) + (-2)*conj(IT_0023)*((
      -0.5)*IT_0021*(s_14*IT_0020 + IT_0012*IT_0023) + IT_0009*IT_0024 + (-0.5)
      *IT_0017*IT_0026) + (-2)*IT_0023*((-0.5)*s_14*conj(IT_0020)*IT_0021 + conj
      (IT_0009)*IT_0024 + (-0.5)*conj(IT_0017)*IT_0026);
    return create_ccomplex_return(IT_0027);
}

