#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_snu_tau_to_A0_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_snu_tau_to_A0_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_24 = param->s_24;
    const creal_t A_tau = param->A_tau;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *U_stau_01;
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = e_em*mu_h*m_tau*IT_0004*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = e_em*mu_h*m_tau*IT_0004*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0006 + 0.5*IT_0008 + (-0.5)*IT_0012 +
       0.5*IT_0015;
    const ccomplex_t IT_0017 = pow(m_A0, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0017 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *U_stau_00;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = e_em*mu_h*m_tau*IT_0004*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = e_em*mu_h*m_tau*IT_0004*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0027 = IT_0001*IT_0026;
    const ccomplex_t IT_0028 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0029 = IT_0009*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0032 = IT_0009*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0025 + (-0.5)*IT_0027 + 0.5*IT_0030 + (
      -0.5)*IT_0033;
    const ccomplex_t IT_0035 = pow(m_stau_2, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0017 + IT_0035 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0023*IT_0038;
    const ccomplex_t IT_0040 = e_em*IT_0001;
    const ccomplex_t IT_0041 = cos(beta);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(A_tau*IT_0009*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0001*(m_W*IT_0009*IT_0041*U_stau_01 + (-0.5)
      *m_tau*IT_0004*(m_tau*IT_0009*IT_0042*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0044 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_12 + -IT_0035 + -IT_0044 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0040*IT_0047;
    const ccomplex_t IT_0049 = IT_0021 + IT_0039 + -IT_0048;
    const ccomplex_t IT_0050 = pow(m_W, -2);
    const ccomplex_t IT_0051 = pow(s_24, 2);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = -IT_0048;
    const ccomplex_t IT_0054 = s_14*s_24*IT_0050;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = pow(s_14, 2);
    const ccomplex_t IT_0057 = IT_0050*IT_0056;
    return create_ccomplex_return(-IT_0049*(conj(IT_0049)*(IT_0044 + -IT_0052)
       + conj(IT_0053)*(s_12 + IT_0055)) + -IT_0053*(conj(IT_0049)*(s_12 +
       IT_0055) + conj(IT_0053)*(IT_0035 + -IT_0057)));
}

