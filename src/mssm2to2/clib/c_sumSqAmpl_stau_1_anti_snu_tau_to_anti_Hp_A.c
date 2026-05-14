#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_snu_tau_to_anti_Hp_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_snu_tau_to_anti_Hp_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = sin(beta);
    const ccomplex_t IT_0002 = cos(beta);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cpow(IT_0002, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(A_tau*IT_0001*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0006*(m_W*IT_0001*IT_0002*U_stau_00 + (-0.5)
      *m_tau*IT_0003*(m_tau*IT_0001*IT_0004*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0008 = pow(m_Hp, 2);
    const ccomplex_t IT_0009 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0010 = cpow((-2)*s_23 + IT_0008 + IT_0009 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_stau_1, 2);
    const ccomplex_t IT_0015 = cpow((-2)*s_12 + -IT_0009 + -IT_0014 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = IT_0013 + IT_0017;
    const ccomplex_t IT_0019 = IT_0000*IT_0018;
    const ccomplex_t IT_0020 = 2*IT_0000;
    const ccomplex_t IT_0021 = IT_0012*IT_0020;
    const ccomplex_t IT_0022 = -s_24;
    const ccomplex_t IT_0023 = IT_0018*IT_0020;
    const ccomplex_t IT_0024 = -s_34;
    const ccomplex_t IT_0025 = -s_23;
    const ccomplex_t IT_0026 = -conj(IT_0019);
    const ccomplex_t IT_0027 = IT_0019*(conj(IT_0021)*IT_0022 + conj(IT_0023)
      *IT_0024) + -IT_0021*(IT_0009*conj(IT_0021) + -conj(IT_0023)*IT_0025 +
       IT_0022*IT_0026) + -IT_0023*(IT_0008*conj(IT_0023) + -conj(IT_0021)
      *IT_0025 + IT_0024*IT_0026);
    return create_ccomplex_return(IT_0027);
}

