#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_snu_tau_to_e_anti_nu_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_snu_tau_to_e_anti_nu_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0000*IT_0002*IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(A_tau*IT_0003*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0005*(m_W*IT_0001*IT_0003*U_stau_01 + (-0.5)
      *m_tau*IT_0000*(m_tau*IT_0002*IT_0003*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = pow(m_stau_2, 2);
    const ccomplex_t IT_0011 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0010 + -IT_0011 + -reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0005
      *U_stau_01;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = m_e*IT_0021;
    const ccomplex_t IT_0023 = -IT_0020;
    const ccomplex_t IT_0024 = m_e*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0022 + IT_0025;
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_e*e_em
      *IT_0000*IT_0005;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = pow(m_tau, 2);
    const ccomplex_t IT_0031 = cpow(IT_0001, 2);
    const ccomplex_t IT_0032 = cpow(IT_0003, 2);
    const ccomplex_t IT_0033 = A_tau*IT_0001*U_stau_11 + 0.707106781186548
      *e_em*IT_0005*((IT_0000*IT_0030 + -m_W*(IT_0031 + -IT_0032))*U_stau_01 + 
      -mu_h*m_tau*IT_0000*IT_0002*IT_0003*U_stau_11);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0029*IT_0035;
    const ccomplex_t IT_0037 = pow(m_W, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + -IT_0010 + -IT_0011 + IT_0037 
      + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0014 + -IT_0027 + -IT_0040;
    const ccomplex_t IT_0042 = IT_0019*IT_0023;
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = m_e*s_14;
    const ccomplex_t IT_0045 = s_13*s_14;
    const ccomplex_t IT_0046 = s_34*IT_0010;
    return create_ccomplex_return(2*IT_0041*(s_34*conj(IT_0041) + conj(IT_0043
      )*IT_0044) + 2*IT_0043*(conj(IT_0041)*IT_0044 + 2*conj(IT_0043)*(IT_0045 +
       (-0.5)*IT_0046)));
}

