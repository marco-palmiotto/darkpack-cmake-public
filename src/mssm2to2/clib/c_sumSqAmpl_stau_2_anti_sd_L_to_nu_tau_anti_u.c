#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_sd_L_to_nu_tau_anti_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_sd_L_to_nu_tau_anti_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0004*IT_0010;
    const ccomplex_t IT_0012 = conj(U_d1)*e_em*m_tau*IT_0000*U_stau_11;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0008 + (-0.5)*IT_0014);
    const ccomplex_t IT_0016 = IT_0006*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0018 = pow(m_stau_2, 2);
    const ccomplex_t IT_0019 = cpow(s_13 + (-0.5)*IT_0018 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0022 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = conj(U_d2)*e_em*m_tau*IT_0000*U_stau_11;
    const ccomplex_t IT_0025 = IT_0011*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + (-0.5)*IT_0026);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_u*conj(V_u2)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = cpow(s_13 + (-0.5)*IT_0018 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0021*IT_0030*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0020 + 0.5*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*V_ud*e_em*U_Wm1*IT_0004;
    const ccomplex_t IT_0035 = IT_0015*IT_0034;
    const ccomplex_t IT_0036 = IT_0019*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*V_ud*e_em*U_Wm2*IT_0004;
    const ccomplex_t IT_0039 = IT_0027*IT_0038;
    const ccomplex_t IT_0040 = IT_0031*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0037 + 0.5*IT_0041;
    const ccomplex_t IT_0043 = m_u*s_13;
    const ccomplex_t IT_0044 = s_13*s_14;
    const ccomplex_t IT_0045 = s_34*IT_0018;
    return create_ccomplex_return(6*IT_0033*(s_34*conj(IT_0033) + -conj
      (IT_0042)*IT_0043) + (-6)*IT_0042*(conj(IT_0033)*IT_0043 + (
      -0.166666666666667)*conj(IT_0042)*(12*IT_0044 + (-6)*IT_0045)));
}

