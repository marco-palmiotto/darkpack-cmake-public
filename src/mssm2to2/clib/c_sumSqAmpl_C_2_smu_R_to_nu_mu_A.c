#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_smu_R_to_nu_mu_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_smu_R_to_nu_mu_A(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_mr = param->Gamma_mr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = pow(m_W, -1);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0001*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-0.5)*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0007;
    const ccomplex_t IT_0009 = pow(m_smu_R, 2);
    const ccomplex_t IT_0010 = cpow(s_23 + (-0.5)*IT_0009 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = s_14*s_34;
    const ccomplex_t IT_0015 = 2*IT_0000;
    const ccomplex_t IT_0016 = pow(m_C_2, 2);
    const ccomplex_t IT_0017 = cpow(s_13 + (-0.5)*IT_0016 + 0.5*m_smu_R*
      (m_smu_R + (0 + _Complex_I*-1)*Gamma_mr) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0007*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = s_13*s_34;
    const ccomplex_t IT_0023 = IT_0000*IT_0019;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = (-2)*s_13;
    const ccomplex_t IT_0026 = s_34*IT_0025;
    const ccomplex_t IT_0027 = IT_0012 + 0.5*IT_0020;
    const ccomplex_t IT_0028 = s_34*IT_0016;
    const ccomplex_t IT_0029 = pow(s_13, 2);
    const ccomplex_t IT_0030 = (-2)*IT_0029;
    const ccomplex_t IT_0031 = 8*IT_0013*conj(IT_0013)*IT_0014 + (-4)*conj
      (IT_0021)*(IT_0013*IT_0022 + (-0.25)*IT_0024*IT_0026) + (-4)*IT_0021*(conj
      (IT_0013)*IT_0022 + (-0.25)*conj(IT_0024)*IT_0026) + (-2)*conj(IT_0027)*((
      -0.5)*IT_0025*(s_14*IT_0024 + IT_0016*IT_0027) + IT_0013*IT_0028 + (-0.5)
      *IT_0021*IT_0030) + (-2)*IT_0027*((-0.5)*s_14*conj(IT_0024)*IT_0025 + conj
      (IT_0013)*IT_0028 + (-0.5)*conj(IT_0021)*IT_0030);
    return create_ccomplex_return(IT_0031);
}

