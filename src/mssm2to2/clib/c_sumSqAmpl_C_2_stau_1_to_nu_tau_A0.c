#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_stau_1_to_nu_tau_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_stau_1_to_nu_tau_A0(
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = conj(U_d2)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0010);
    const ccomplex_t IT_0012 = e_em*mu_h*m_tau*IT_0007*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = e_em*mu_h*m_tau*IT_0007*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = sin(beta);
    const ccomplex_t IT_0017 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0013 + 0.5*IT_0015 + (-0.5)*IT_0019 +
       0.5*IT_0022;
    const ccomplex_t IT_0024 = IT_0011*IT_0023;
    const ccomplex_t IT_0025 = pow(m_C_2, 2);
    const ccomplex_t IT_0026 = cpow(s_13 + (-0.5)*IT_0025 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = IT_0001*(IT_0030 + IT_0032);
    const ccomplex_t IT_0034 = 1.4142135623731*e_em*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_10;
    const ccomplex_t IT_0039 = IT_0006*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = IT_0035*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0044 = pow(m_stau_1, 2);
    const ccomplex_t IT_0045 = cpow(s_23 + (-0.5)*IT_0044 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0042*IT_0043*IT_0045;
    const ccomplex_t IT_0047 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0048 = IT_0004*IT_0047;
    const ccomplex_t IT_0049 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0050 = IT_0016*IT_0049;
    const ccomplex_t IT_0051 = IT_0001*(IT_0048 + IT_0050);
    const ccomplex_t IT_0052 = 1.4142135623731*e_em*IT_0051;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = e_em*conj(U_Wm2)*U_stau_00;
    const ccomplex_t IT_0055 = IT_0001*IT_0054;
    const ccomplex_t IT_0056 = conj(U_d2)*e_em*m_tau*IT_0007*U_stau_10;
    const ccomplex_t IT_0057 = IT_0006*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + (-0.5)*IT_0058);
    const ccomplex_t IT_0060 = IT_0053*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0062 = cpow(s_23 + (-0.5)*IT_0044 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0060*IT_0061*IT_0062;
    const ccomplex_t IT_0064 = V_u2*U_Wm1;
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = U_d1*V_Wp2;
    const ccomplex_t IT_0067 = IT_0016*IT_0066;
    const ccomplex_t IT_0068 = IT_0001*(IT_0065 + IT_0067);
    const ccomplex_t IT_0069 = 1.4142135623731*e_em*IT_0068;
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = IT_0041*IT_0070;
    const ccomplex_t IT_0072 = IT_0045*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = m_C_2*IT_0073;
    const ccomplex_t IT_0075 = V_u2*U_Wm2;
    const ccomplex_t IT_0076 = IT_0004*IT_0075;
    const ccomplex_t IT_0077 = U_d2*V_Wp2;
    const ccomplex_t IT_0078 = IT_0016*IT_0077;
    const ccomplex_t IT_0079 = IT_0001*(IT_0076 + IT_0078);
    const ccomplex_t IT_0080 = 1.4142135623731*e_em*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0059*IT_0081;
    const ccomplex_t IT_0083 = IT_0062*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = m_C_2*IT_0084;
    const ccomplex_t IT_0086 = e_em*mu_h*m_tau*IT_0007*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = e_em*mu_h*m_tau*IT_0007*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0089 = IT_0001*IT_0088;
    const ccomplex_t IT_0090 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0091 = IT_0016*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0094 = IT_0016*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0087 + 0.5*IT_0089 + (-0.5)*IT_0092 +
       0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0059*IT_0096;
    const ccomplex_t IT_0098 = cpow(s_13 + (-0.5)*IT_0025 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0028 + (-0.5)*IT_0046 + (-0.5)*IT_0063 +
       (-0.5)*IT_0074 + (-0.5)*IT_0085 + 0.5*IT_0100;
    const ccomplex_t IT_0102 = 0.5*IT_0073 + 0.5*IT_0084;
    const ccomplex_t IT_0103 = s_34*m_C_2;
    const ccomplex_t IT_0104 = pow(m_A0, 2);
    const ccomplex_t IT_0105 = s_13*IT_0104;
    const ccomplex_t IT_0106 = s_14*s_34;
    return create_ccomplex_return(2*IT_0101*(s_13*conj(IT_0101) + conj(IT_0102
      )*IT_0103) + 2*IT_0102*(conj(IT_0101)*IT_0103 + -conj(IT_0102)*(IT_0105 + 
      (-2)*IT_0106)));
}

