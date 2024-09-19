#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_se_L_snu_tau_to_nu_e_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_se_L_snu_tau_to_nu_e_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0001;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0006 = pow(m_se_L, 2);
    const ccomplex_t IT_0007 = cpow(s_13 + (-0.5)*IT_0006 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0001;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0014 = cpow(s_13 + (-0.5)*IT_0006 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em
      *m_tau*IT_0001*IT_0016*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = IT_0014*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = m_tau*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*U_d1*e_em
      *m_tau*IT_0001*IT_0016*IT_0018;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = IT_0007*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = m_tau*IT_0029;
    const ccomplex_t IT_0031 = (-0.5)*IT_0015 + 0.5*IT_0024 + 0.5*IT_0030;
    const ccomplex_t IT_0032 = 0.5*IT_0023 + 0.5*IT_0029;
    const ccomplex_t IT_0033 = s_23*m_tau;
    const ccomplex_t IT_0034 = s_23*s_24;
    const ccomplex_t IT_0035 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0036 = s_34*IT_0035;
    const ccomplex_t IT_0037 = (-2)*s_34;
    const ccomplex_t IT_0038 = 2*s_34*(IT_0009*conj(IT_0009) + IT_0031*conj
      (IT_0031)) + 2*((conj(IT_0009) + -conj(IT_0031))*IT_0032 + (IT_0009 + 
      -IT_0031)*conj(IT_0032))*IT_0033 + 4*IT_0032*conj(IT_0032)*(IT_0034 + (
      -0.5)*IT_0036) + (conj(IT_0009)*IT_0031 + IT_0009*conj(IT_0031))*IT_0037;
    return create_ccomplex_return(IT_0038);
}

