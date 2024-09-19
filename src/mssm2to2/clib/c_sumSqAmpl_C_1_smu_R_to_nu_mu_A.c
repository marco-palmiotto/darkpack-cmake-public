#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_smu_R_to_nu_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_smu_R_to_nu_mu_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_mu = param->m_mu;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_C_1, 2);
    const ccomplex_t IT_0009 = cpow(s_13 + (-0.5)*IT_0008 + 0.5*m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = 2*IT_0000;
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = IT_0000*IT_0007;
    const ccomplex_t IT_0017 = pow(m_smu_R, 2);
    const ccomplex_t IT_0018 = cpow(s_23 + (-0.5)*IT_0017 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0015 + IT_0020;
    const ccomplex_t IT_0022 = (-2)*s_13;
    const ccomplex_t IT_0023 = (-0.5)*IT_0015;
    const ccomplex_t IT_0024 = s_34*IT_0022;
    const ccomplex_t IT_0025 = s_13*s_34;
    const ccomplex_t IT_0026 = (-0.5)*IT_0020;
    const ccomplex_t IT_0027 = pow(s_13, 2);
    const ccomplex_t IT_0028 = (-2)*IT_0027;
    const ccomplex_t IT_0029 = s_14*s_34;
    const ccomplex_t IT_0030 = s_34*IT_0008;
    const ccomplex_t IT_0031 = (-2)*IT_0030;
    const ccomplex_t IT_0032 = IT_0013*(s_14*conj(IT_0021)*IT_0022 + conj
      (IT_0023)*IT_0024) + (-4)*IT_0023*((-0.25)*conj(IT_0013)*IT_0024 + IT_0025
      *conj(IT_0026) + (-0.25)*conj(IT_0021)*IT_0028) + (-4)*IT_0026*(conj
      (IT_0023)*IT_0025 + (-2)*conj(IT_0026)*IT_0029 + (-0.25)*conj(IT_0021)
      *IT_0031) + IT_0021*((s_14*conj(IT_0013) + IT_0008*conj(IT_0021))*IT_0022 
      + conj(IT_0023)*IT_0028 + conj(IT_0026)*IT_0031);
    return create_ccomplex_return(IT_0032);
}

