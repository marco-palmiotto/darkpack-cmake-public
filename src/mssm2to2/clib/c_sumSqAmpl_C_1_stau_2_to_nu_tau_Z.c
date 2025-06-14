#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_stau_2_to_nu_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_stau_2_to_nu_tau_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0000*IT_0013;
    const ccomplex_t IT_0015 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0012;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0016 + -IT_0018 + (-2)
      *IT_0020);
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0025 = pow(m_stau_2, 2);
    const ccomplex_t IT_0026 = cpow(s_23 + (-0.5)*IT_0025 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0024*IT_0026;
    const ccomplex_t IT_0028 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0029 = IT_0014*IT_0028;
    const ccomplex_t IT_0030 = IT_0017*IT_0028;
    const ccomplex_t IT_0031 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0032 = IT_0017*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*(IT_0029 + -IT_0030 + (-2)
      *IT_0032);
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0011*IT_0034;
    const ccomplex_t IT_0036 = IT_0026*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = m_C_1*IT_0037;
    const ccomplex_t IT_0039 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0040 = IT_0001*IT_0039;
    const ccomplex_t IT_0041 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0042 = IT_0006*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + (-0.5)*IT_0043);
    const ccomplex_t IT_0045 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0046 = IT_0014*IT_0045;
    const ccomplex_t IT_0047 = IT_0017*IT_0045;
    const ccomplex_t IT_0048 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0049 = IT_0017*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + -IT_0047 + (-2)
      *IT_0049);
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0044*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0054 = cpow(s_23 + (-0.5)*IT_0025 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0057 = IT_0014*IT_0056;
    const ccomplex_t IT_0058 = IT_0017*IT_0056;
    const ccomplex_t IT_0059 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0060 = IT_0017*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + -IT_0058 + (-2)
      *IT_0060);
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0044*IT_0062;
    const ccomplex_t IT_0064 = IT_0054*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = m_C_1*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0027 + 0.5*IT_0038 + 0.5*IT_0055 + 0.5
      *IT_0066;
    const ccomplex_t IT_0068 = pow(m_Z, -2);
    const ccomplex_t IT_0069 = pow(m_Z, 2);
    const ccomplex_t IT_0070 = s_13*IT_0069;
    const ccomplex_t IT_0071 = s_14*s_34;
    const ccomplex_t IT_0072 = (-2)*IT_0071;
    const ccomplex_t IT_0073 = IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = IT_0068*IT_0073;
    const ccomplex_t IT_0075 = e_em*IT_0014;
    const ccomplex_t IT_0076 = e_em*IT_0017;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*(IT_0075 + IT_0076);
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0044*IT_0078;
    const ccomplex_t IT_0080 = pow(m_C_1, 2);
    const ccomplex_t IT_0081 = cpow(s_12 + 0.5*IT_0025 + 0.5*IT_0080 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0012 + -IT_0000*IT_0013)*conj(U_stau_00)*U_stau_01 + IT_0000*IT_0013
      *conj(U_stau_10)*U_stau_11);
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0087 = IT_0001*IT_0086;
    const ccomplex_t IT_0088 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_10;
    const ccomplex_t IT_0089 = IT_0006*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-0.5)*IT_0090);
    const ccomplex_t IT_0092 = cpow(s_13 + (-0.5)*IT_0080 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0085*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0012 + -IT_0000*IT_0013)*U_stau_01*conj(U_stau_01) + IT_0000*IT_0013
      *U_stau_11*conj(U_stau_11));
    const ccomplex_t IT_0097 = (-2)*IT_0096;
    const ccomplex_t IT_0098 = cpow(s_13 + (-0.5)*IT_0080 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0099 = IT_0044*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0097*IT_0100;
    const ccomplex_t IT_0102 = -IT_0083 + (-0.5)*IT_0095 + (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = pow(s_34, 2);
    const ccomplex_t IT_0104 = m_C_1*IT_0068*IT_0103;
    const ccomplex_t IT_0105 = 2*IT_0104;
    const ccomplex_t IT_0106 = (-0.5)*IT_0083;
    const ccomplex_t IT_0107 = s_34*IT_0069;
    const ccomplex_t IT_0108 = IT_0068*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = s_34 + IT_0109;
    const ccomplex_t IT_0111 = m_C_1*IT_0110;
    const ccomplex_t IT_0112 = 2*IT_0111;
    const ccomplex_t IT_0113 = 2*IT_0112;
    const ccomplex_t IT_0114 = 0.5*IT_0037 + 0.5*IT_0065;
    const ccomplex_t IT_0115 = (-0.25)*IT_0108;
    const ccomplex_t IT_0116 = s_34 + IT_0115;
    const ccomplex_t IT_0117 = m_C_1*IT_0116;
    const ccomplex_t IT_0118 = (-4)*IT_0117;
    const ccomplex_t IT_0119 = 2*IT_0118;
    const ccomplex_t IT_0120 = -IT_0037 + -IT_0065 + 0.5*IT_0095 + 0.5*IT_0101;
    const ccomplex_t IT_0121 = s_14*s_34*IT_0068;
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = s_13 + IT_0122;
    const ccomplex_t IT_0124 = m_C_1*IT_0123;
    const ccomplex_t IT_0125 = (-2)*IT_0124;
    const ccomplex_t IT_0126 = s_13*s_34;
    const ccomplex_t IT_0127 = 2*IT_0126;
    const ccomplex_t IT_0128 = pow(m_Z, 4);
    const ccomplex_t IT_0129 = s_13*IT_0128;
    const ccomplex_t IT_0130 = IT_0068*IT_0129;
    const ccomplex_t IT_0131 = (-0.25)*IT_0130;
    const ccomplex_t IT_0132 = IT_0071 + IT_0131;
    const ccomplex_t IT_0133 = (-8)*IT_0132;
    const ccomplex_t IT_0134 = IT_0070 + IT_0072 + (-0.5)*IT_0130;
    const ccomplex_t IT_0135 = (-4)*IT_0134;
    const ccomplex_t IT_0136 = s_13*s_14;
    const ccomplex_t IT_0137 = s_34*IT_0080;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = conj(IT_0120)*IT_0139;
    const ccomplex_t IT_0141 = pow(s_14, 2);
    const ccomplex_t IT_0142 = IT_0068*IT_0141;
    const ccomplex_t IT_0143 = -IT_0080;
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = s_13*IT_0144;
    const ccomplex_t IT_0146 = s_13*IT_0123;
    const ccomplex_t IT_0147 = (-2)*IT_0146;
    const ccomplex_t IT_0148 = 0.5*conj(IT_0102);
    const ccomplex_t IT_0149 = 0.5*conj(IT_0067);
    const ccomplex_t IT_0150 = (-2)*IT_0126;
    const ccomplex_t IT_0151 = 0.5*conj(IT_0114);
    const ccomplex_t IT_0152 = s_13*IT_0068*IT_0103;
    const ccomplex_t IT_0153 = 4*IT_0067*(conj(IT_0067)*(s_13 + (-0.5)*IT_0074
      ) + 0.25*conj(IT_0102)*IT_0105 + 0.25*conj(IT_0106)*IT_0113 + 0.25*conj
      (IT_0114)*IT_0119 + 0.25*conj(IT_0120)*IT_0125) + IT_0106*(conj(IT_0067)
      *IT_0113 + conj(IT_0102)*IT_0127 + conj(IT_0114)*IT_0133 + conj(IT_0106)
      *IT_0135 + (-2)*IT_0140) + (-2)*IT_0120*((conj(IT_0106) + -conj(IT_0114))
      *IT_0139 + -conj(IT_0120)*IT_0145 + -IT_0147*IT_0148 + -IT_0125*IT_0149) +
       IT_0114*(conj(IT_0106)*IT_0133 + 2*IT_0140 + 2*IT_0119*IT_0149 + 2
      *IT_0148*IT_0150 + 2*IT_0135*IT_0151) + 2*IT_0102*(0.5*conj(IT_0106)
      *IT_0127 + 0.5*conj(IT_0120)*IT_0147 + IT_0105*IT_0149 + IT_0150*IT_0151 +
       conj(IT_0102)*IT_0152);
    return create_ccomplex_return(IT_0153);
}

