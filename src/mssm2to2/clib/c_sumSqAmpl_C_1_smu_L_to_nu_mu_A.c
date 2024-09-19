#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_smu_L_to_nu_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_smu_L_to_nu_mu_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_smu_L = param->m_smu_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_ml = param->Gamma_ml;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0002;
    const ccomplex_t IT_0004 = pow(m_C_1, 2);
    const ccomplex_t IT_0005 = cpow(s_13 + (-0.5)*IT_0004 + 0.5*m_smu_L*
      (m_smu_L + (0 + _Complex_I*-1)*Gamma_ml) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = 2*IT_0000;
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = IT_0000*IT_0003;
    const ccomplex_t IT_0013 = pow(m_smu_L, 2);
    const ccomplex_t IT_0014 = cpow(s_23 + (-0.5)*IT_0013 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = 0.5*IT_0011 + IT_0016;
    const ccomplex_t IT_0018 = (-2)*s_13;
    const ccomplex_t IT_0019 = (-0.5)*IT_0011;
    const ccomplex_t IT_0020 = s_34*IT_0018;
    const ccomplex_t IT_0021 = s_13*s_34;
    const ccomplex_t IT_0022 = (-0.5)*IT_0016;
    const ccomplex_t IT_0023 = pow(s_13, 2);
    const ccomplex_t IT_0024 = (-2)*IT_0023;
    const ccomplex_t IT_0025 = s_14*s_34;
    const ccomplex_t IT_0026 = s_34*IT_0004;
    const ccomplex_t IT_0027 = (-2)*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*(s_14*conj(IT_0017)*IT_0018 + conj
      (IT_0019)*IT_0020) + (-4)*IT_0019*((-0.25)*conj(IT_0009)*IT_0020 + IT_0021
      *conj(IT_0022) + (-0.25)*conj(IT_0017)*IT_0024) + (-4)*IT_0022*(conj
      (IT_0019)*IT_0021 + (-2)*conj(IT_0022)*IT_0025 + (-0.25)*conj(IT_0017)
      *IT_0027) + IT_0017*((s_14*conj(IT_0009) + IT_0004*conj(IT_0017))*IT_0018 
      + conj(IT_0019)*IT_0024 + conj(IT_0022)*IT_0027);
    return create_ccomplex_return(IT_0028);
}

