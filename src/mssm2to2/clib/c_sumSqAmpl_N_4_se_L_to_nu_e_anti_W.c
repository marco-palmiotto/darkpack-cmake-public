#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_se_L_to_nu_e_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_se_L_to_nu_e_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_e = param->m_e;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_se_L = param->m_se_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em*conj(U_Wm2)*IT_0001;
    const ccomplex_t IT_0003 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*(IT_0004 + 0.5*IT_0007);
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = pow(m_se_L, 2);
    const ccomplex_t IT_0012 = cpow(s_23 + (-0.5)*IT_0011 + 0.5*m_C_2*(m_C_2 +
       (0 + _Complex_I*-1)*Gamma_c2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = m_N_4*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0001;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = pow(m_W, -1);
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*m_e*N_d4*e_em
      *IT_0001*IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = pow(m_e, 2);
    const ccomplex_t IT_0026 = pow(m_N_4, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + -IT_0011 + IT_0025 + -IT_0026 
      + -reg_prop, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_e;
    const ccomplex_t IT_0029 = IT_0024*IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*e_em*conj(U_Wm1)*IT_0001;
    const ccomplex_t IT_0031 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + (-0.5)*IT_0035);
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0030*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0040 = cpow(s_23 + (-0.5)*IT_0011 + 0.5*m_C_1*(m_C_1 +
       (0 + _Complex_I*-1)*Gamma_c1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0039*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0043 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0044 = IT_0001*IT_0043;
    const ccomplex_t IT_0045 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0046 = IT_0001*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + (-0.5)*IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = IT_0002*IT_0049;
    const ccomplex_t IT_0051 = IT_0012*IT_0042*IT_0050;
    const ccomplex_t IT_0052 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0053 = IT_0001*IT_0052;
    const ccomplex_t IT_0054 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0055 = IT_0001*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 0.5*IT_0056);
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0030*IT_0058;
    const ccomplex_t IT_0060 = IT_0040*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = m_N_4*IT_0061;
    const ccomplex_t IT_0063 = -IT_0029 + 0.5*IT_0041 + 0.5*IT_0051 + (-0.5)
      *IT_0062;
    const ccomplex_t IT_0064 = pow(m_W, -2);
    const ccomplex_t IT_0065 = pow(m_W, 2);
    const ccomplex_t IT_0066 = s_13*IT_0065;
    const ccomplex_t IT_0067 = s_14*s_34;
    const ccomplex_t IT_0068 = (-2)*IT_0067;
    const ccomplex_t IT_0069 = IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = IT_0064*IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = s_13 + IT_0071;
    const ccomplex_t IT_0073 = conj(IT_0016) + conj(IT_0063);
    const ccomplex_t IT_0074 = 0.5*IT_0014;
    const ccomplex_t IT_0075 = (-0.5)*IT_0061;
    const ccomplex_t IT_0076 = pow(m_W, 4);
    const ccomplex_t IT_0077 = s_13*IT_0076;
    const ccomplex_t IT_0078 = IT_0064*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = IT_0066 + IT_0068 + IT_0079;
    const ccomplex_t IT_0081 = -conj(IT_0075);
    const ccomplex_t IT_0082 = conj(IT_0074) + IT_0081;
    const ccomplex_t IT_0083 = cos(theta_W);
    const ccomplex_t IT_0084 = cpow(IT_0083, -1);
    const ccomplex_t IT_0085 = conj(N_B4)*e_em;
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W4)*e_em;
    const ccomplex_t IT_0089 = IT_0001*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + IT_0090);
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = IT_0018*IT_0092;
    const ccomplex_t IT_0094 = IT_0027*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = (-0.25)*IT_0078;
    const ccomplex_t IT_0097 = IT_0067 + IT_0096;
    const ccomplex_t IT_0098 = s_34*IT_0065;
    const ccomplex_t IT_0099 = IT_0064*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = s_34 + IT_0100;
    const ccomplex_t IT_0102 = m_N_4*IT_0101;
    const ccomplex_t IT_0103 = 2*IT_0102;
    const ccomplex_t IT_0104 = 2*conj(IT_0016);
    const ccomplex_t IT_0105 = 2*conj(IT_0063);
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = s_14*s_34*IT_0064;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = s_13 + IT_0108;
    const ccomplex_t IT_0110 = s_13*IT_0109;
    const ccomplex_t IT_0111 = -IT_0017;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0087 + -IT_0090);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = cpow(s_13 + (-0.5)*IT_0026 + (-0.5)*reg_prop +
       0.5*m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = IT_0111*IT_0116;
    const ccomplex_t IT_0118 = 2*IT_0095 + (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = -IT_0014;
    const ccomplex_t IT_0120 = 2*IT_0119;
    const ccomplex_t IT_0121 = IT_0061 + 0.5*IT_0117;
    const ccomplex_t IT_0122 = (-2)*IT_0121;
    const ccomplex_t IT_0123 = IT_0120 + IT_0122;
    const ccomplex_t IT_0124 = 2*conj(IT_0119);
    const ccomplex_t IT_0125 = (-2)*conj(IT_0121);
    const ccomplex_t IT_0126 = IT_0124 + IT_0125;
    const ccomplex_t IT_0127 = m_N_4*IT_0109;
    const ccomplex_t IT_0128 = pow(s_34, 2);
    const ccomplex_t IT_0129 = s_13*IT_0064*IT_0128;
    const ccomplex_t IT_0130 = 2*IT_0118;
    const ccomplex_t IT_0131 = m_N_4*IT_0064*IT_0128;
    const ccomplex_t IT_0132 = 2*conj(IT_0118);
    const ccomplex_t IT_0133 = s_13*s_34;
    const ccomplex_t IT_0134 = conj(IT_0074) + conj(IT_0095);
    const ccomplex_t IT_0135 = IT_0074 + IT_0095;
    const ccomplex_t IT_0136 = (-0.25)*IT_0099;
    const ccomplex_t IT_0137 = s_34 + IT_0136;
    const ccomplex_t IT_0138 = m_N_4*IT_0137;
    const ccomplex_t IT_0139 = (-4)*IT_0138;
    const ccomplex_t IT_0140 = 2*IT_0063;
    const ccomplex_t IT_0141 = pow(s_14, 2);
    const ccomplex_t IT_0142 = IT_0064*IT_0141;
    const ccomplex_t IT_0143 = -IT_0026;
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = s_13*IT_0144;
    const ccomplex_t IT_0146 = 2*IT_0121;
    const ccomplex_t IT_0147 = (-0.5)*IT_0146;
    const ccomplex_t IT_0148 = IT_0119 + IT_0147;
    const ccomplex_t IT_0149 = (-2)*conj(IT_0119);
    const ccomplex_t IT_0150 = s_13*s_14;
    const ccomplex_t IT_0151 = s_34*IT_0026;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = IT_0150 + IT_0152;
    const ccomplex_t IT_0154 = 4*(IT_0016 + IT_0063)*IT_0072*IT_0073 + (-4)*
      (IT_0074 + -IT_0075)*IT_0080*IT_0082 + 8*IT_0095*(IT_0082*IT_0097 + 0.125
      *IT_0103*IT_0106) + IT_0110*(conj(IT_0118)*IT_0123 + IT_0118*IT_0126) + 
      (IT_0073*IT_0123 + (IT_0016 + IT_0063)*IT_0126)*IT_0127 + conj(IT_0118)
      *IT_0129*IT_0130 + IT_0131*(IT_0073*IT_0130 + (IT_0016 + IT_0063)*IT_0132)
       + (-2)*IT_0133*(conj(IT_0075)*IT_0118 + IT_0075*conj(IT_0118) + (-0.5)
      *IT_0130*IT_0134 + (-0.5)*IT_0132*IT_0135) + (-2)*IT_0139*(IT_0073*IT_0074
       + IT_0063*conj(IT_0074) + IT_0016*IT_0082 + (-0.5)*IT_0075*IT_0106 + (
      -0.5)*conj(IT_0075)*IT_0140) + (-4)*conj(IT_0095)*(IT_0080*IT_0095 + (-2)*
      (IT_0074 + -IT_0075)*IT_0097 + (-0.5)*IT_0103*(IT_0016 + 0.5*IT_0140)) + (
      -2)*IT_0145*((-0.5)*conj(IT_0119)*IT_0120 + conj(IT_0121)*IT_0148 + (-0.5)
      *IT_0121*IT_0149) + (IT_0123*IT_0134 + IT_0126*IT_0135 + 2*IT_0081*IT_0148
       + 2*IT_0075*(conj(IT_0121) + 0.5*IT_0149))*IT_0153;
    return create_ccomplex_return(IT_0154);
}

