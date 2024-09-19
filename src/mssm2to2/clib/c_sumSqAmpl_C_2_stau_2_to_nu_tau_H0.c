#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_stau_2_to_nu_tau_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_stau_2_to_nu_tau_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = conj(U_d2)*e_em*m_tau*IT_0007*U_stau_10;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0010);
    const ccomplex_t IT_0012 = cos(alpha);
    const ccomplex_t IT_0013 = mu_h*m_tau*IT_0007;
    const ccomplex_t IT_0014 = sin(alpha);
    const ccomplex_t IT_0015 = sin(beta);
    const ccomplex_t IT_0016 = -IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0004*IT_0012;
    const ccomplex_t IT_0018 = IT_0016 + IT_0017;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -2);
    const ccomplex_t IT_0021 = pow(m_tau, 2);
    const ccomplex_t IT_0022 = (-2)*conj(U_stau_00)*(((-0.5)*m_W*IT_0018*
      (IT_0001 + -IT_0000*IT_0020) + IT_0001*IT_0005*IT_0007*IT_0012*IT_0021)
      *U_stau_01 + -1./2*IT_0001*IT_0005*IT_0013*IT_0014*U_stau_11) + conj
      (U_stau_10)*(IT_0001*IT_0005*IT_0013*IT_0014*U_stau_01 + 2*(m_W*IT_0000
      *IT_0018*IT_0020 + -IT_0001*IT_0005*IT_0007*IT_0012*IT_0021)*U_stau_11);
    const ccomplex_t IT_0023 = (0 + _Complex_I*-1)*(e_em*IT_0022 + (
      -1.4142135623731)*A_tau*IT_0012*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0011*IT_0024;
    const ccomplex_t IT_0026 = pow(m_C_2, 2);
    const ccomplex_t IT_0027 = cpow(s_13 + (-0.5)*IT_0026 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0031 = IT_0012*IT_0030;
    const ccomplex_t IT_0032 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0033 = IT_0014*IT_0032;
    const ccomplex_t IT_0034 = IT_0001*(IT_0031 + IT_0033);
    const ccomplex_t IT_0035 = 1.4142135623731*e_em*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0039 = IT_0001*IT_0038;
    const ccomplex_t IT_0040 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + (-0.5)*IT_0042);
    const ccomplex_t IT_0044 = IT_0037*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0046 = pow(m_stau_2, 2);
    const ccomplex_t IT_0047 = cpow(s_23 + (-0.5)*IT_0046 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0044*IT_0045*IT_0047;
    const ccomplex_t IT_0049 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0050 = IT_0012*IT_0049;
    const ccomplex_t IT_0051 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0052 = IT_0014*IT_0051;
    const ccomplex_t IT_0053 = IT_0001*(IT_0050 + IT_0052);
    const ccomplex_t IT_0054 = 1.4142135623731*e_em*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0058 = IT_0001*IT_0057;
    const ccomplex_t IT_0059 = conj(U_d2)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0060 = IT_0006*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-0.5)*IT_0061);
    const ccomplex_t IT_0063 = IT_0056*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0065 = cpow(s_23 + (-0.5)*IT_0046 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0063*IT_0064*IT_0065;
    const ccomplex_t IT_0067 = U_d1*V_Wp2;
    const ccomplex_t IT_0068 = IT_0012*IT_0067;
    const ccomplex_t IT_0069 = V_u2*U_Wm1;
    const ccomplex_t IT_0070 = IT_0014*IT_0069;
    const ccomplex_t IT_0071 = IT_0001*(IT_0068 + IT_0070);
    const ccomplex_t IT_0072 = 1.4142135623731*e_em*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0043*IT_0074;
    const ccomplex_t IT_0076 = IT_0047*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = m_C_2*IT_0077;
    const ccomplex_t IT_0079 = U_d2*V_Wp2;
    const ccomplex_t IT_0080 = IT_0012*IT_0079;
    const ccomplex_t IT_0081 = V_u2*U_Wm2;
    const ccomplex_t IT_0082 = IT_0014*IT_0081;
    const ccomplex_t IT_0083 = IT_0001*(IT_0080 + IT_0082);
    const ccomplex_t IT_0084 = 1.4142135623731*e_em*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0062*IT_0086;
    const ccomplex_t IT_0088 = IT_0065*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = m_C_2*IT_0089;
    const ccomplex_t IT_0091 = cpow(s_13 + (-0.5)*IT_0026 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0092 = (-2)*U_stau_01*(((-0.5)*m_W*IT_0018*(IT_0001 + 
      -IT_0000*IT_0020) + IT_0001*IT_0005*IT_0007*IT_0012*IT_0021)*conj
      (U_stau_01) + -1./2*IT_0001*IT_0005*IT_0013*IT_0014*conj(U_stau_11)) +
       U_stau_11*(IT_0001*IT_0005*IT_0013*IT_0014*conj(U_stau_01) + 2*(m_W
      *IT_0000*IT_0018*IT_0020 + -IT_0001*IT_0005*IT_0007*IT_0012*IT_0021)*conj
      (U_stau_11));
    const ccomplex_t IT_0093 = (0 + _Complex_I*-1)*(e_em*IT_0092 + (
      -1.4142135623731)*A_tau*IT_0012*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0062*IT_0094;
    const ccomplex_t IT_0096 = IT_0091*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0029 + (-0.5)*IT_0048 + (-0.5)*IT_0066 +
       (-0.5)*IT_0078 + (-0.5)*IT_0090 + 0.5*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0077 + 0.5*IT_0089;
    const ccomplex_t IT_0100 = s_34*m_C_2;
    const ccomplex_t IT_0101 = pow(m_H0, 2);
    const ccomplex_t IT_0102 = s_13*IT_0101;
    const ccomplex_t IT_0103 = s_14*s_34;
    return create_ccomplex_return(2*IT_0098*(s_13*conj(IT_0098) + conj(IT_0099
      )*IT_0100) + 2*IT_0099*(conj(IT_0098)*IT_0100 + -conj(IT_0099)*(IT_0102 + 
      (-2)*IT_0103)));
}

