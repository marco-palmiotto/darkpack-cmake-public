#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_stau_2_to_nu_tau_anti_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_stau_2_to_nu_tau_anti_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t A_tau = param->A_tau;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = conj(N_u2)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = 2*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0000*IT_0008;
    const ccomplex_t IT_0010 = conj(N_B2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = conj(N_W2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = IT_0012 + IT_0015;
    const ccomplex_t IT_0017 = IT_0006 + IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0021 = IT_0002*IT_0020;
    const ccomplex_t IT_0022 = cpow(IT_0000, -1);
    const ccomplex_t IT_0023 = IT_0002*IT_0022;
    const ccomplex_t IT_0024 = pow(m_W, -1);
    const ccomplex_t IT_0025 = conj(U_d1)*e_em*m_tau*IT_0024*U_stau_11;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0021 + (-0.5)*IT_0027);
    const ccomplex_t IT_0029 = IT_0019*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0031 = pow(m_stau_2, 2);
    const ccomplex_t IT_0032 = cpow(s_23 + (-0.5)*IT_0031 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0029*IT_0030*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B2)*e_em;
    const ccomplex_t IT_0035 = IT_0008*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = conj(N_W2)*e_em;
    const ccomplex_t IT_0038 = IT_0002*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + -IT_0039);
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = sin(beta);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(A_tau*IT_0042*U_stau_11 + (
      -1.4142135623731)*e_em*IT_0002*(m_W*IT_0000*IT_0042*U_stau_01 + (-0.5)
      *m_tau*IT_0024*(m_tau*IT_0022*IT_0042*U_stau_01 + mu_h*U_stau_11)));
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = pow(m_N_2, 2);
    const ccomplex_t IT_0046 = cpow(s_13 + (-0.5)*IT_0045 + (-0.5)*reg_prop +
       0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)*Gamma_nutaul), -1);
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = pow(m_tau, 2);
    const ccomplex_t IT_0050 = cpow(s_12 + 0.5*IT_0031 + 0.5*IT_0045 + (-0.5)
      *IT_0049 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0051 = (0 + _Complex_I*1.4142135623731)*e_em*m_tau
      *IT_0002*IT_0022*IT_0024*IT_0042;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0054 = IT_0008*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_d2*e_em*m_tau*IT_0024*U_stau_01;
    const ccomplex_t IT_0057 = IT_0023*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + 0.5*IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = IT_0052*IT_0060;
    const ccomplex_t IT_0062 = IT_0050*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = m_N_2*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0066 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0067 = IT_0008*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0070 = IT_0002*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_d2)*e_em*m_tau*IT_0024*U_stau_11;
    const ccomplex_t IT_0073 = IT_0023*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + IT_0071 + 
      -IT_0074);
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = IT_0052*IT_0076;
    const ccomplex_t IT_0078 = IT_0050*IT_0065*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0080 = cpow(s_23 + (-0.5)*IT_0031 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0081 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0082 = IT_0002*IT_0081;
    const ccomplex_t IT_0083 = conj(U_d2)*e_em*m_tau*IT_0024*U_stau_11;
    const ccomplex_t IT_0084 = IT_0023*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-0.5)*IT_0085);
    const ccomplex_t IT_0087 = conj(N_u2)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0088 = IT_0003*IT_0087;
    const ccomplex_t IT_0089 = 2*IT_0088;
    const ccomplex_t IT_0090 = conj(N_B2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0091 = IT_0009*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0094 = IT_0003*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = IT_0092 + IT_0095;
    const ccomplex_t IT_0097 = IT_0089 + IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = 0.5*IT_0098;
    const ccomplex_t IT_0100 = IT_0086*IT_0099;
    const ccomplex_t IT_0101 = IT_0079*IT_0080*IT_0100;
    const ccomplex_t IT_0102 = (-0.5)*IT_0033 + 0.5*IT_0048 + 0.5*IT_0064 + (
      -0.5)*IT_0078 + (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0002*IT_0042;
    const ccomplex_t IT_0104 = N_d2*e_em*U_Wm1;
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (-2)*IT_0105;
    const ccomplex_t IT_0107 = IT_0008*IT_0042;
    const ccomplex_t IT_0108 = N_B2*U_d1*e_em;
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = N_W2*U_d1*e_em;
    const ccomplex_t IT_0112 = IT_0103*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = IT_0110 + IT_0113;
    const ccomplex_t IT_0115 = IT_0106 + IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0028*IT_0117;
    const ccomplex_t IT_0119 = IT_0032*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = N_d2*e_em*U_Wm2;
    const ccomplex_t IT_0122 = IT_0103*IT_0121;
    const ccomplex_t IT_0123 = (-2)*IT_0122;
    const ccomplex_t IT_0124 = N_B2*U_d2*e_em;
    const ccomplex_t IT_0125 = IT_0107*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = N_W2*U_d2*e_em;
    const ccomplex_t IT_0128 = IT_0103*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = IT_0126 + IT_0129;
    const ccomplex_t IT_0131 = IT_0123 + IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0086*IT_0133;
    const ccomplex_t IT_0135 = IT_0080*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0063 + (-0.5)*IT_0120 + (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = s_23*m_N_2;
    const ccomplex_t IT_0139 = s_12*s_23;
    const ccomplex_t IT_0140 = s_13*IT_0031;
    return create_ccomplex_return(2*IT_0102*(s_13*conj(IT_0102) + conj(IT_0137
      )*IT_0138) + 2*IT_0137*(conj(IT_0102)*IT_0138 + 2*conj(IT_0137)*(IT_0139 +
       (-0.5)*IT_0140)));
}

