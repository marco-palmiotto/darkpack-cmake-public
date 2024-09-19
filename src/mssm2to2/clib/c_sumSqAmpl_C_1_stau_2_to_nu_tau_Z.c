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
    const ccomplex_t IT_0002 = e_em*conj(U_Wm1)*U_stau_01;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0010);
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0000*IT_0013;
    const ccomplex_t IT_0015 = e_em*IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0012;
    const ccomplex_t IT_0017 = e_em*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0015 + IT_0017);
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0019;
    const ccomplex_t IT_0021 = pow(m_C_1, 2);
    const ccomplex_t IT_0022 = pow(m_stau_2, 2);
    const ccomplex_t IT_0023 = cpow(s_12 + 0.5*IT_0021 + 0.5*IT_0022 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0012 + -IT_0000*IT_0013)*conj(U_stau_00)*U_stau_01 + IT_0000*IT_0013
      *conj(U_stau_10)*U_stau_11);
    const ccomplex_t IT_0028 = (-2)*IT_0027;
    const ccomplex_t IT_0029 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = conj(U_d1)*e_em*m_tau*IT_0007*U_stau_10;
    const ccomplex_t IT_0032 = IT_0006*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + (-0.5)*IT_0033);
    const ccomplex_t IT_0035 = cpow(s_13 + (-0.5)*IT_0021 + (-0.5)*reg_prop +
       0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1), -1);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0028*IT_0037;
    const ccomplex_t IT_0039 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0040 = IT_0014*IT_0039;
    const ccomplex_t IT_0041 = IT_0016*IT_0039;
    const ccomplex_t IT_0042 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0043 = IT_0016*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + -IT_0041 + (-2)
      *IT_0043);
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = IT_0011*IT_0045;
    const ccomplex_t IT_0047 = cpow(s_23 + (-0.5)*IT_0022 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = e_em*conj(U_Wm2)*U_stau_01;
    const ccomplex_t IT_0051 = IT_0001*IT_0050;
    const ccomplex_t IT_0052 = conj(U_d2)*e_em*m_tau*IT_0007*U_stau_11;
    const ccomplex_t IT_0053 = IT_0006*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + (-0.5)*IT_0054);
    const ccomplex_t IT_0056 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0057 = IT_0014*IT_0056;
    const ccomplex_t IT_0058 = IT_0016*IT_0056;
    const ccomplex_t IT_0059 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0060 = IT_0016*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + -IT_0058 + (-2)
      *IT_0060);
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = IT_0055*IT_0062;
    const ccomplex_t IT_0064 = cpow(s_23 + (-0.5)*IT_0022 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0001
      *IT_0012 + -IT_0000*IT_0013)*U_stau_01*conj(U_stau_01) + IT_0000*IT_0013
      *U_stau_11*conj(U_stau_11));
    const ccomplex_t IT_0068 = (-2)*IT_0067;
    const ccomplex_t IT_0069 = cpow(s_13 + (-0.5)*IT_0021 + (-0.5)*reg_prop +
       0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2), -1);
    const ccomplex_t IT_0070 = IT_0011*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0068*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0038 + -IT_0049 + -IT_0066 + 0.5*IT_0072;
    const ccomplex_t IT_0074 = s_13*s_14;
    const ccomplex_t IT_0075 = s_34*IT_0021;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0074 + IT_0076;
    const ccomplex_t IT_0078 = -IT_0025 + (-0.5)*IT_0038 + (-0.5)*IT_0072;
    const ccomplex_t IT_0079 = s_13*s_34;
    const ccomplex_t IT_0080 = 2*IT_0079;
    const ccomplex_t IT_0081 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0082 = IT_0014*IT_0081;
    const ccomplex_t IT_0083 = IT_0016*IT_0081;
    const ccomplex_t IT_0084 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0085 = IT_0016*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + -IT_0083 + (-2)
      *IT_0085);
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0055*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0090 = IT_0064*IT_0088*IT_0089;
    const ccomplex_t IT_0091 = m_C_1*IT_0066;
    const ccomplex_t IT_0092 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0093 = IT_0014*IT_0092;
    const ccomplex_t IT_0094 = IT_0016*IT_0092;
    const ccomplex_t IT_0095 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0096 = IT_0016*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + -IT_0094 + (-2)
      *IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = IT_0011*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0101 = IT_0047*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = m_C_1*IT_0049;
    const ccomplex_t IT_0103 = 0.5*IT_0090 + 0.5*IT_0091 + 0.5*IT_0101 + 0.5
      *IT_0102;
    const ccomplex_t IT_0104 = pow(m_Z, -2);
    const ccomplex_t IT_0105 = pow(m_Z, 2);
    const ccomplex_t IT_0106 = s_34*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = s_34 + IT_0108;
    const ccomplex_t IT_0110 = m_C_1*IT_0109;
    const ccomplex_t IT_0111 = 2*IT_0110;
    const ccomplex_t IT_0112 = 2*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0049 + 0.5*IT_0066;
    const ccomplex_t IT_0114 = s_14*s_34;
    const ccomplex_t IT_0115 = pow(m_Z, 4);
    const ccomplex_t IT_0116 = s_13*IT_0115;
    const ccomplex_t IT_0117 = IT_0104*IT_0116;
    const ccomplex_t IT_0118 = (-0.25)*IT_0117;
    const ccomplex_t IT_0119 = IT_0114 + IT_0118;
    const ccomplex_t IT_0120 = (-8)*IT_0119;
    const ccomplex_t IT_0121 = s_13*IT_0105;
    const ccomplex_t IT_0122 = (-2)*IT_0114;
    const ccomplex_t IT_0123 = (-0.5)*IT_0117 + IT_0121 + IT_0122;
    const ccomplex_t IT_0124 = (-4)*IT_0123;
    const ccomplex_t IT_0125 = (-0.25)*IT_0107;
    const ccomplex_t IT_0126 = s_34 + IT_0125;
    const ccomplex_t IT_0127 = m_C_1*IT_0126;
    const ccomplex_t IT_0128 = (-4)*IT_0127;
    const ccomplex_t IT_0129 = 2*IT_0128;
    const ccomplex_t IT_0130 = IT_0121 + IT_0122;
    const ccomplex_t IT_0131 = IT_0104*IT_0130;
    const ccomplex_t IT_0132 = pow(s_34, 2);
    const ccomplex_t IT_0133 = m_C_1*IT_0104*IT_0132;
    const ccomplex_t IT_0134 = 2*IT_0133;
    const ccomplex_t IT_0135 = s_14*s_34*IT_0104;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = s_13 + IT_0136;
    const ccomplex_t IT_0138 = m_C_1*IT_0137;
    const ccomplex_t IT_0139 = (-2)*IT_0138;
    const ccomplex_t IT_0140 = s_13*IT_0104*IT_0132;
    const ccomplex_t IT_0141 = s_13*IT_0137;
    const ccomplex_t IT_0142 = (-2)*IT_0141;
    const ccomplex_t IT_0143 = 0.5*conj(IT_0103);
    const ccomplex_t IT_0144 = pow(s_14, 2);
    const ccomplex_t IT_0145 = IT_0104*IT_0144;
    const ccomplex_t IT_0146 = -IT_0021;
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = s_13*IT_0147;
    const ccomplex_t IT_0149 = (-2)*IT_0026*(conj(IT_0073)*IT_0077 + (-0.5)
      *conj(IT_0078)*IT_0080 + (-0.5)*conj(IT_0103)*IT_0112 + (-0.5)*conj
      (IT_0113)*IT_0120 + (-0.5)*conj(IT_0026)*IT_0124) + 2*IT_0113*(conj
      (IT_0073)*IT_0077 + -conj(IT_0078)*IT_0079 + 0.5*conj(IT_0026)*IT_0120 +
       0.5*conj(IT_0113)*IT_0124 + 0.5*conj(IT_0103)*IT_0129) + 4*IT_0103*(0.25
      *conj(IT_0026)*IT_0112 + 0.25*conj(IT_0113)*IT_0129 + conj(IT_0103)*(s_13 
      + (-0.5)*IT_0131) + 0.25*conj(IT_0078)*IT_0134 + 0.25*conj(IT_0073)
      *IT_0139) + 2*IT_0078*(0.5*conj(IT_0026)*IT_0080 + -IT_0079*conj(IT_0113) 
      + conj(IT_0078)*IT_0140 + 0.5*conj(IT_0073)*IT_0142 + IT_0134*IT_0143) + (
      -2)*IT_0073*(IT_0077*(conj(IT_0026) + -conj(IT_0113)) + (-0.5)*conj
      (IT_0078)*IT_0142 + -IT_0139*IT_0143 + -conj(IT_0073)*IT_0148);
    return create_ccomplex_return(IT_0149);
}

