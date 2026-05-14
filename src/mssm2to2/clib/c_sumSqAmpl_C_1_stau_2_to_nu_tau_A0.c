#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_stau_2_to_nu_tau_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_stau_2_to_nu_tau_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = V_u1*U_Wm1;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = U_d1*V_Wp1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0000, -1);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = pow(m_W, -1);
    const ccomplex_t IT_0016 = conj(U_d1)*e_em*m_tau*IT_0015*U_stau_11;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0012 + (-0.5)*IT_0018);
    const ccomplex_t IT_0020 = IT_0010*IT_0019;
    const ccomplex_t IT_0021 = pow(m_stau_2, 2);
    const ccomplex_t IT_0022 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = V_u1*U_Wm2;
    const ccomplex_t IT_0026 = IT_0000*IT_0025;
    const ccomplex_t IT_0027 = U_d2*V_Wp1;
    const ccomplex_t IT_0028 = IT_0003*IT_0027;
    const ccomplex_t IT_0029 = IT_0007*(IT_0026 + IT_0028);
    const ccomplex_t IT_0030 = 1.4142135623731*e_em*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = conj(U_d2)*e_em*m_tau*IT_0015*U_stau_11;
    const ccomplex_t IT_0035 = IT_0014*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + (-0.5)*IT_0036);
    const ccomplex_t IT_0038 = IT_0031*IT_0037;
    const ccomplex_t IT_0039 = cpow(s_23 + (-0.5)*IT_0021 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0024 + 0.5*IT_0041;
    const ccomplex_t IT_0043 = pow(m_A0, 2);
    const ccomplex_t IT_0044 = s_13*IT_0043;
    const ccomplex_t IT_0045 = s_14*s_34;
    const ccomplex_t IT_0046 = e_em*mu_h*m_tau*IT_0015*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0047 = IT_0007*IT_0046;
    const ccomplex_t IT_0048 = e_em*mu_h*m_tau*IT_0015*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0049 = IT_0007*IT_0048;
    const ccomplex_t IT_0050 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0051 = IT_0003*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (-0.5)*IT_0047 + 0.5*IT_0049 + (-0.5)*IT_0052 +
       0.5*IT_0055;
    const ccomplex_t IT_0057 = IT_0019*IT_0056;
    const ccomplex_t IT_0058 = pow(m_C_1, 2);
    const ccomplex_t IT_0059 = cpow(s_13 + (-0.5)*IT_0058 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0060 = IT_0057*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0063 = IT_0007*IT_0062;
    const ccomplex_t IT_0064 = conj(U_d1)*e_em*m_tau*IT_0015*U_stau_10;
    const ccomplex_t IT_0065 = IT_0014*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + (-0.5)*IT_0066);
    const ccomplex_t IT_0068 = e_em*mu_h*m_tau*IT_0015*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0069 = IT_0007*IT_0068;
    const ccomplex_t IT_0070 = e_em*mu_h*m_tau*IT_0015*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0071 = IT_0007*IT_0070;
    const ccomplex_t IT_0072 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0073 = IT_0003*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0076 = IT_0003*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0069 + (-0.5)*IT_0071 + 0.5*IT_0074 + (
      -0.5)*IT_0077;
    const ccomplex_t IT_0079 = IT_0067*IT_0078;
    const ccomplex_t IT_0080 = cpow(s_13 + (-0.5)*IT_0058 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0081 = IT_0079*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*IT_0081;
    const ccomplex_t IT_0083 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0084 = IT_0000*IT_0083;
    const ccomplex_t IT_0085 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0086 = IT_0003*IT_0085;
    const ccomplex_t IT_0087 = IT_0007*(IT_0084 + IT_0086);
    const ccomplex_t IT_0088 = 1.4142135623731*e_em*IT_0087;
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0019*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0092 = IT_0022*IT_0090*IT_0091;
    const ccomplex_t IT_0093 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0094 = IT_0000*IT_0093;
    const ccomplex_t IT_0095 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = IT_0007*(IT_0094 + IT_0096);
    const ccomplex_t IT_0098 = 1.4142135623731*e_em*IT_0097;
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0037*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0102 = IT_0039*IT_0100*IT_0101;
    const ccomplex_t IT_0103 = m_C_1*IT_0024;
    const ccomplex_t IT_0104 = m_C_1*IT_0041;
    const ccomplex_t IT_0105 = 0.5*IT_0061 + 0.5*IT_0082 + (-0.5)*IT_0092 + (
      -0.5)*IT_0102 + (-0.5)*IT_0103 + (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = s_34*m_C_1;
    const ccomplex_t IT_0107 = 2*IT_0106;
    return create_ccomplex_return(2*IT_0105*(s_13*conj(IT_0105) + 0.5*conj
      (IT_0042)*IT_0107) + (-2)*IT_0042*(conj(IT_0042)*(IT_0044 + (-2)*IT_0045) 
      + (-0.5)*conj(IT_0105)*IT_0107));
}

