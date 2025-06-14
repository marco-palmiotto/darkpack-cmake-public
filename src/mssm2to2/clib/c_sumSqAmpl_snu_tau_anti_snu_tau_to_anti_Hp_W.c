#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_snu_tau_anti_snu_tau_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_snu_tau_anti_snu_tau_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
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
      *conj(U_stau_00);
    const ccomplex_t IT_0003 = -IT_0002;
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cpow(IT_0005, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(A_tau*IT_0004*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0001*(m_W*IT_0004*IT_0005*U_stau_00 + (-0.5)
      *m_tau*IT_0006*(m_tau*IT_0004*IT_0007*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0009 = pow(m_Hp, 2);
    const ccomplex_t IT_0010 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_23 + IT_0009 + IT_0010 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001
      *conj(U_stau_01);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(A_tau*IT_0004*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0001*(m_W*IT_0004*IT_0005*U_stau_01 + (-0.5)
      *m_tau*IT_0006*(m_tau*IT_0004*IT_0007*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0018 = cpow((-2)*s_23 + IT_0009 + IT_0010 + reg_prop +
       -m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = -IT_0014 + -IT_0021;
    const ccomplex_t IT_0023 = pow(m_W, -2);
    const ccomplex_t IT_0024 = pow(s_34, 2);
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = cos(alpha);
    const ccomplex_t IT_0027 = sin(alpha);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0004
      *IT_0026 + -IT_0005*IT_0027);
    const ccomplex_t IT_0029 = cos(theta_W);
    const ccomplex_t IT_0030 = cpow(IT_0029, -2);
    const ccomplex_t IT_0031 = IT_0000*IT_0030;
    const ccomplex_t IT_0032 = IT_0001 + IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_W*e_em*(IT_0005*IT_0026 + 
      -IT_0004*IT_0027)*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = cpow((-2)*s_12 + (-2)*IT_0010 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0028*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*e_em*IT_0001*(IT_0005
      *IT_0026 + IT_0004*IT_0027);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*m_W*e_em*(IT_0004*IT_0026 +
       IT_0005*IT_0027)*IT_0032;
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0040*IT_0045;
    const ccomplex_t IT_0047 = IT_0014 + IT_0021 + -IT_0038 + -IT_0046;
    const ccomplex_t IT_0048 = s_24*s_34*IT_0023;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = s_23 + IT_0049;
    const ccomplex_t IT_0051 = s_14*s_34*IT_0023;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = s_13 + IT_0052;
    const ccomplex_t IT_0054 = -IT_0038 + -IT_0046;
    const ccomplex_t IT_0055 = -conj(IT_0054);
    const ccomplex_t IT_0056 = pow(s_14, 2);
    const ccomplex_t IT_0057 = IT_0023*IT_0056;
    const ccomplex_t IT_0058 = s_14*s_24*IT_0023;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = s_12 + IT_0059;
    const ccomplex_t IT_0061 = -conj(IT_0022);
    const ccomplex_t IT_0062 = pow(s_24, 2);
    const ccomplex_t IT_0063 = IT_0023*IT_0062;
    const ccomplex_t IT_0064 = -IT_0022*(conj(IT_0022)*(IT_0009 + -IT_0025) +
       conj(IT_0047)*IT_0050 + -IT_0053*IT_0055) + -IT_0054*(conj(IT_0054)*
      (IT_0010 + -IT_0057) + conj(IT_0047)*IT_0060 + -IT_0053*IT_0061) + IT_0047
      *(IT_0055*IT_0060 + IT_0050*IT_0061 + -conj(IT_0047)*(IT_0010 + -IT_0063));
    return create_ccomplex_return(IT_0064);
}

