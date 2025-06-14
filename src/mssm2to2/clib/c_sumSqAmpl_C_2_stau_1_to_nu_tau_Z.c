#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_stau_1_to_nu_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_stau_1_to_nu_tau_Z(
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
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0000*IT_0013;
    const ccomplex_t IT_0015 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0001*IT_0012;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0016 + -IT_0018 + (-2)
      *IT_0020);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = pow(m_stau_1, 2);
    const ccomplex_t IT_0025 = cpow(s_23 + (-0.5)*IT_0024 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0012 + -IT_0000*IT_0013)*U_stau_00*conj(U_stau_00) + IT_0000*IT_0013
      *U_stau_10*conj(U_stau_10));
    const ccomplex_t IT_0029 = (-2)*IT_0028;
    const ccomplex_t IT_0030 = pow(m_C_2, 2);
    const ccomplex_t IT_0031 = cpow(s_13 + (-0.5)*IT_0030 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0032 = IT_0011*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0012 + -IT_0000*IT_0013)*U_stau_00*conj(U_stau_01) + IT_0000*IT_0013
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0036 = (-2)*IT_0035;
    const ccomplex_t IT_0037 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0038 = IT_0001*IT_0037;
    const ccomplex_t IT_0039 = conj(U_d2)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0040 = IT_0006*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + (-0.5)*IT_0041);
    const ccomplex_t IT_0043 = cpow(s_13 + (-0.5)*IT_0030 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0036*IT_0045;
    const ccomplex_t IT_0047 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0048 = IT_0001*IT_0047;
    const ccomplex_t IT_0049 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_10;
    const ccomplex_t IT_0050 = IT_0006*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + (-0.5)*IT_0051);
    const ccomplex_t IT_0053 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0054 = IT_0014*IT_0053;
    const ccomplex_t IT_0055 = IT_0017*IT_0053;
    const ccomplex_t IT_0056 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0057 = IT_0017*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0054 + -IT_0055 + (-2)
      *IT_0057);
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = IT_0052*IT_0059;
    const ccomplex_t IT_0061 = cpow(s_23 + (-0.5)*IT_0024 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = -IT_0027 + 0.5*IT_0034 + 0.5*IT_0046 + -IT_0063;
    const ccomplex_t IT_0065 = pow(m_Z, -2);
    const ccomplex_t IT_0066 = pow(s_14, 2);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = -IT_0030;
    const ccomplex_t IT_0069 = IT_0067 + IT_0068;
    const ccomplex_t IT_0070 = s_13*IT_0069;
    const ccomplex_t IT_0071 = e_em*IT_0014;
    const ccomplex_t IT_0072 = e_em*IT_0017;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0071 + IT_0072);
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = IT_0011*IT_0074;
    const ccomplex_t IT_0076 = cpow(s_12 + 0.5*IT_0024 + 0.5*IT_0030 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0034 + (-0.5)*IT_0046 + -IT_0078;
    const ccomplex_t IT_0080 = pow(s_34, 2);
    const ccomplex_t IT_0081 = s_13*IT_0065*IT_0080;
    const ccomplex_t IT_0082 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0083 = IT_0014*IT_0082;
    const ccomplex_t IT_0084 = IT_0017*IT_0082;
    const ccomplex_t IT_0085 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0086 = IT_0017*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(IT_0083 + -IT_0084 + (-2)
      *IT_0086);
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0011*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0091 = IT_0025*IT_0089*IT_0090;
    const ccomplex_t IT_0092 = m_C_2*IT_0063;
    const ccomplex_t IT_0093 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0094 = IT_0014*IT_0093;
    const ccomplex_t IT_0095 = IT_0017*IT_0093;
    const ccomplex_t IT_0096 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0097 = IT_0017*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + -IT_0095 + (-2)
      *IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0052*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0102 = IT_0061*IT_0100*IT_0101;
    const ccomplex_t IT_0103 = m_C_2*IT_0027;
    const ccomplex_t IT_0104 = 0.5*IT_0091 + 0.5*IT_0092 + 0.5*IT_0102 + 0.5
      *IT_0103;
    const ccomplex_t IT_0105 = m_C_2*IT_0065*IT_0080;
    const ccomplex_t IT_0106 = (-0.5)*IT_0078;
    const ccomplex_t IT_0107 = s_13*s_34;
    const ccomplex_t IT_0108 = 2*IT_0079;
    const ccomplex_t IT_0109 = pow(m_Z, 2);
    const ccomplex_t IT_0110 = s_34*IT_0109;
    const ccomplex_t IT_0111 = IT_0065*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = s_34 + IT_0112;
    const ccomplex_t IT_0114 = m_C_2*IT_0113;
    const ccomplex_t IT_0115 = 2*IT_0114;
    const ccomplex_t IT_0116 = 2*IT_0104;
    const ccomplex_t IT_0117 = s_14*s_34*IT_0065;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = s_13 + IT_0118;
    const ccomplex_t IT_0120 = s_13*IT_0119;
    const ccomplex_t IT_0121 = m_C_2*IT_0119;
    const ccomplex_t IT_0122 = s_13*s_14;
    const ccomplex_t IT_0123 = s_34*IT_0030;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = (-2)*IT_0064;
    const ccomplex_t IT_0127 = (-2)*conj(IT_0064);
    const ccomplex_t IT_0128 = s_13*IT_0109;
    const ccomplex_t IT_0129 = s_14*s_34;
    const ccomplex_t IT_0130 = (-2)*IT_0129;
    const ccomplex_t IT_0131 = IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = IT_0065*IT_0131;
    const ccomplex_t IT_0133 = 2*conj(IT_0079);
    const ccomplex_t IT_0134 = pow(m_Z, 4);
    const ccomplex_t IT_0135 = s_13*IT_0134;
    const ccomplex_t IT_0136 = IT_0065*IT_0135;
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = IT_0128 + IT_0130 + IT_0137;
    const ccomplex_t IT_0139 = 2*conj(IT_0104);
    const ccomplex_t IT_0140 = (-0.25)*IT_0136;
    const ccomplex_t IT_0141 = IT_0129 + IT_0140;
    const ccomplex_t IT_0142 = 0.5*IT_0027 + 0.5*IT_0063;
    const ccomplex_t IT_0143 = 2*conj(IT_0142);
    const ccomplex_t IT_0144 = (-0.25)*IT_0111;
    const ccomplex_t IT_0145 = s_34 + IT_0144;
    const ccomplex_t IT_0146 = m_C_2*IT_0145;
    const ccomplex_t IT_0147 = (-4)*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = 2*IT_0064*conj(IT_0064)*IT_0070 + (conj(IT_0079
      )*IT_0081 + conj(IT_0104)*IT_0105 + conj(IT_0106)*IT_0107)*IT_0108 + conj
      (IT_0106)*IT_0115*IT_0116 + (conj(IT_0079)*IT_0120 + conj(IT_0104)*IT_0121
       + conj(IT_0106)*IT_0125)*IT_0126 + IT_0079*IT_0120*IT_0127 + 4*IT_0104*
      (0.25*IT_0121*IT_0127 + conj(IT_0104)*(s_13 + (-0.5)*IT_0132) + 0.25
      *IT_0105*IT_0133) + (-4)*IT_0106*((-0.25)*IT_0125*IT_0127 + (-0.25)
      *IT_0107*IT_0133 + conj(IT_0106)*IT_0138 + (-0.25)*IT_0115*IT_0139 +
       IT_0141*IT_0143) + (-2)*conj(IT_0142)*(IT_0079*IT_0107 + -IT_0064*IT_0125
       + IT_0116*IT_0148) + (-2)*IT_0142*(conj(IT_0079)*IT_0107 + -conj(IT_0064)
      *IT_0125 + 4*conj(IT_0106)*IT_0141 + IT_0138*IT_0143 + IT_0139*IT_0148);
    return create_ccomplex_return(IT_0149);
}

