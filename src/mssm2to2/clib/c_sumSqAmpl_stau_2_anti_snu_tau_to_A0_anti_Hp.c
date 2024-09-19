#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_snu_tau_to_A0_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_snu_tau_to_A0_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t A_tau = param->A_tau;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = pow(m_W, -2);
    const ccomplex_t IT_0005 = pow(m_tau, 2);
    const ccomplex_t IT_0006 = cpow(IT_0001, -2);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -2);
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0000*(IT_0002 + -IT_0003*
      (IT_0003 + -IT_0003*IT_0004*IT_0005*IT_0006))*IT_0008*U_stau_01;
    const ccomplex_t IT_0010 = (-0.25)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = cpow(IT_0001, -1);
    const ccomplex_t IT_0013 = cpow(IT_0007, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(A_tau*IT_0003*U_stau_10 + (
      -1.4142135623731)*e_em*IT_0013*(m_W*IT_0001*IT_0003*U_stau_00 + (-0.5)
      *m_tau*IT_0011*(m_tau*IT_0003*IT_0012*U_stau_00 + mu_h*U_stau_10)));
    const ccomplex_t IT_0015 = e_em*mu_h*m_tau*IT_0011*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0016 = IT_0013*IT_0015;
    const ccomplex_t IT_0017 = e_em*mu_h*m_tau*IT_0011*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0018 = IT_0013*IT_0017;
    const ccomplex_t IT_0019 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0020 = IT_0003*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = 0.5*IT_0016 + (-0.5)*IT_0018 + 0.5*IT_0021 + (
      -0.5)*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*IT_0025;
    const ccomplex_t IT_0027 = pow(m_A0, 2);
    const ccomplex_t IT_0028 = pow(m_stau_2, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0027 + IT_0028 + reg_prop +
       -m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0030 = IT_0026*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = e_em*mu_h*m_tau*IT_0011*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0033 = IT_0013*IT_0032;
    const ccomplex_t IT_0034 = e_em*mu_h*m_tau*IT_0011*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0035 = IT_0013*IT_0034;
    const ccomplex_t IT_0036 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0040 = IT_0003*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0033 + 0.5*IT_0035 + (-0.5)*IT_0038 +
       0.5*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(A_tau*IT_0003*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0013*(m_W*IT_0001*IT_0003*U_stau_01 + (-0.5)
      *m_tau*IT_0011*(m_tau*IT_0003*IT_0012*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = cpow((-2)*s_13 + IT_0027 + reg_prop + (0 +
       _Complex_I*1)*m_stau_2*Gamma_tau2, -1);
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = pow(m_snu_tau, 2);
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0028 + -IT_0048 + -reg_prop, -1);
    const ccomplex_t IT_0050 = s_13 + s_23;
    const ccomplex_t IT_0051 = e_em*IT_0013;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0013
      *U_stau_01;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = s_12 + IT_0048;
    const ccomplex_t IT_0058 = (-0.5)*IT_0051;
    const ccomplex_t IT_0059 = -IT_0053;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = IT_0057*IT_0060;
    const ccomplex_t IT_0062 = IT_0056 + IT_0061;
    const ccomplex_t IT_0063 = IT_0052*IT_0059;
    const ccomplex_t IT_0064 = s_23*IT_0063;
    const ccomplex_t IT_0065 = IT_0054*IT_0058;
    const ccomplex_t IT_0066 = 0.5*IT_0028 + 0.5*IT_0048;
    const ccomplex_t IT_0067 = s_12 + IT_0066;
    const ccomplex_t IT_0068 = IT_0065*IT_0067;
    const ccomplex_t IT_0069 = -IT_0064 + (-2)*IT_0068;
    const ccomplex_t IT_0070 = IT_0062 + IT_0069;
    const ccomplex_t IT_0071 = IT_0049*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = pow(m_W, 2);
    const ccomplex_t IT_0074 = cpow((-2)*s_12 + -IT_0028 + -IT_0048 + IT_0073 
      + -reg_prop, -1);
    const ccomplex_t IT_0075 = cpow(IT_0003, 2);
    const ccomplex_t IT_0076 = A_tau*IT_0001*U_stau_11 + 0.707106781186548
      *e_em*IT_0013*((IT_0005*IT_0011 + -m_W*(IT_0002 + -IT_0075))*U_stau_01 + 
      -mu_h*m_tau*IT_0003*IT_0011*IT_0012*U_stau_11);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = m_W*e_em*IT_0013;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = IT_0078*IT_0080;
    const ccomplex_t IT_0082 = IT_0074*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = IT_0010 + IT_0031 + IT_0047 + IT_0072 + -IT_0083;
    return create_ccomplex_return(IT_0084*conj(IT_0084));
}

