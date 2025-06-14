#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_smu_R_anti_stau_1_to_nu_mu_anti_nu_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_smu_R_anti_stau_1_to_nu_mu_anti_nu_tau(
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
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t m_smu_R = param->m_smu_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*conj(U_d2)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = e_em*U_Wm2*conj(U_stau_00);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = IT_0002*IT_0004;
    const ccomplex_t IT_0010 = U_d2*e_em*m_tau*IT_0000*conj(U_stau_10);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0012);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = pow(m_smu_R, 2);
    const ccomplex_t IT_0016 = cpow(s_13 + (-0.5)*IT_0015 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = e_em*U_Wm1*conj(U_stau_00);
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = U_d1*e_em*m_tau*IT_0000*conj(U_stau_10);
    const ccomplex_t IT_0024 = IT_0009*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0022 + (-0.5)*IT_0025);
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = cpow(s_13 + (-0.5)*IT_0015 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0018 + 0.5*IT_0030;
    const ccomplex_t IT_0032 = s_13*s_14;
    const ccomplex_t IT_0033 = s_34*IT_0015;
    return create_ccomplex_return(4*IT_0031*conj(IT_0031)*(IT_0032 + (-0.5)
      *IT_0033));
}

