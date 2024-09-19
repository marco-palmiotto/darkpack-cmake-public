#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_sd_R_to_u_anti_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_sd_R_to_u_anti_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_d = param->m_d;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t IT_0000 = pow(m_d, 2);
    const ccomplex_t IT_0001 = pow(m_N_2, 2);
    const ccomplex_t IT_0002 = pow(m_sd_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0006;
    const ccomplex_t IT_0008 = (-0.333333333333333)*IT_0007;
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = IT_0008*IT_0012;
    const ccomplex_t IT_0014 = IT_0003*IT_0004*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0016 = pow(m_u, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_23 + IT_0002 + IT_0016 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0018 = N_W2*e_em*conj(V_Wp1);
    const ccomplex_t IT_0019 = IT_0010*IT_0018;
    const ccomplex_t IT_0020 = N_u2*conj(V_u1)*e_em;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0019 + (-0.5)*IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0010*IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = IT_0015*IT_0017*IT_0030;
    const ccomplex_t IT_0032 = conj(N_W2)*e_em*U_Wm1;
    const ccomplex_t IT_0033 = IT_0010*IT_0032;
    const ccomplex_t IT_0034 = conj(N_d2)*U_d1*e_em;
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + 0.5*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = IT_0017*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = m_N_2*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0044 = IT_0010*IT_0043;
    const ccomplex_t IT_0045 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0046 = IT_0010*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + 0.5*IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0010*IT_0025*IT_0027;
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0002 + IT_0016 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = m_N_2*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0058 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0059 = IT_0010*IT_0058;
    const ccomplex_t IT_0060 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0061 = IT_0010*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + (-0.5)*IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = IT_0051*IT_0064;
    const ccomplex_t IT_0066 = IT_0053*IT_0057*IT_0065;
    const ccomplex_t IT_0067 = -IT_0014 + -IT_0031 + IT_0042 + IT_0056 + 
      -IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1.4142135623731)*m_d*conj(N_d2)
      *e_em*IT_0010*IT_0025*IT_0027;
    const ccomplex_t IT_0069 = 0.5*IT_0068;
    const ccomplex_t IT_0070 = IT_0012*IT_0069;
    const ccomplex_t IT_0071 = IT_0003*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = m_u*IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = m_u*m_N_2;
    const ccomplex_t IT_0076 = (-2)*IT_0041 + (-2)*IT_0055;
    const ccomplex_t IT_0077 = pow(m_W, -2);
    const ccomplex_t IT_0078 = s_14*s_34*IT_0077;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = s_13 + IT_0079;
    const ccomplex_t IT_0081 = m_N_2*IT_0080;
    const ccomplex_t IT_0082 = 2*IT_0072;
    const ccomplex_t IT_0083 = pow(s_34, 2);
    const ccomplex_t IT_0084 = IT_0077*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0016 + IT_0085;
    const ccomplex_t IT_0087 = m_N_2*IT_0086;
    const ccomplex_t IT_0088 = pow(m_W, 2);
    const ccomplex_t IT_0089 = s_34*IT_0088;
    const ccomplex_t IT_0090 = IT_0077*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = s_34 + IT_0091;
    const ccomplex_t IT_0093 = m_N_2*IT_0092;
    const ccomplex_t IT_0094 = 2*IT_0093;
    const ccomplex_t IT_0095 = IT_0041 + IT_0055;
    const ccomplex_t IT_0096 = (-0.25)*IT_0090;
    const ccomplex_t IT_0097 = s_34 + IT_0096;
    const ccomplex_t IT_0098 = m_N_2*IT_0097;
    const ccomplex_t IT_0099 = (-4)*IT_0098;
    const ccomplex_t IT_0100 = s_13*IT_0088;
    const ccomplex_t IT_0101 = s_14*s_34;
    const ccomplex_t IT_0102 = (-2)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = IT_0077*IT_0103;
    const ccomplex_t IT_0105 = (-6)*IT_0104;
    const ccomplex_t IT_0106 = 12*s_13 + IT_0105;
    const ccomplex_t IT_0107 = pow(s_14, 2);
    const ccomplex_t IT_0108 = IT_0077*IT_0107;
    const ccomplex_t IT_0109 = -IT_0001;
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = s_13*IT_0110;
    const ccomplex_t IT_0112 = m_u*IT_0110;
    const ccomplex_t IT_0113 = s_13*s_14;
    const ccomplex_t IT_0114 = s_34*IT_0001;
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = 6*IT_0076;
    const ccomplex_t IT_0118 = s_14*IT_0088;
    const ccomplex_t IT_0119 = IT_0077*IT_0118;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = s_14 + IT_0120;
    const ccomplex_t IT_0122 = m_u*IT_0121;
    const ccomplex_t IT_0123 = 2*IT_0122;
    const ccomplex_t IT_0124 = 6*IT_0074;
    const ccomplex_t IT_0125 = (-6)*conj(IT_0067);
    const ccomplex_t IT_0126 = s_14*IT_0016;
    const ccomplex_t IT_0127 = s_13*s_34;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = s_13*IT_0080;
    const ccomplex_t IT_0131 = (-6)*conj(IT_0076);
    const ccomplex_t IT_0132 = m_u*IT_0080;
    const ccomplex_t IT_0133 = (-6)*IT_0082;
    const ccomplex_t IT_0134 = (-6)*IT_0076;
    const ccomplex_t IT_0135 = s_13*IT_0086;
    const ccomplex_t IT_0136 = (-0.25)*IT_0119;
    const ccomplex_t IT_0137 = s_14 + IT_0136;
    const ccomplex_t IT_0138 = m_u*IT_0137;
    const ccomplex_t IT_0139 = (-4)*IT_0138;
    const ccomplex_t IT_0140 = 6*conj(IT_0076);
    const ccomplex_t IT_0141 = (-6)*conj(IT_0082);
    const ccomplex_t IT_0142 = pow(m_W, 4);
    const ccomplex_t IT_0143 = s_13*IT_0142;
    const ccomplex_t IT_0144 = IT_0077*IT_0143;
    const ccomplex_t IT_0145 = 6*IT_0144;
    const ccomplex_t IT_0146 = 24*IT_0101;
    const ccomplex_t IT_0147 = (-12)*IT_0100;
    const ccomplex_t IT_0148 = IT_0145 + IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = 6*conj(IT_0074);
    const ccomplex_t IT_0150 = (-24)*IT_0101;
    const ccomplex_t IT_0151 = IT_0145 + IT_0150;
    const ccomplex_t IT_0152 = 0.166666666666667*conj(IT_0072);
    const ccomplex_t IT_0153 = 0.166666666666667*IT_0116;
    const ccomplex_t IT_0154 = IT_0067*((-18)*conj(IT_0074)*IT_0075 + (-6)
      *conj(IT_0076)*IT_0081 + (-6)*conj(IT_0082)*IT_0087 + 6*conj(IT_0072)
      *IT_0094 + 6*conj(IT_0095)*IT_0099 + conj(IT_0067)*IT_0106) + (conj
      (IT_0076)*IT_0111 + conj(IT_0074)*IT_0112 + conj(IT_0095)*IT_0116)*IT_0117
       + conj(IT_0095)*IT_0123*IT_0124 + IT_0076*IT_0081*IT_0125 + 6*IT_0082*
      (0.166666666666667*IT_0087*IT_0125 + conj(IT_0095)*IT_0129 +
       0.166666666666667*IT_0130*IT_0131) + conj(IT_0074)*IT_0132*IT_0133 + conj
      (IT_0082)*(IT_0130*IT_0134 + IT_0133*IT_0135) + conj(IT_0072)*(IT_0129
      *IT_0133 + IT_0116*IT_0134 + IT_0124*IT_0139) + IT_0074*((-18)*conj
      (IT_0067)*IT_0075 + conj(IT_0074)*IT_0106 + IT_0112*IT_0140 + IT_0132
      *IT_0141) + 6*IT_0072*(conj(IT_0067)*IT_0094 + 0.166666666666667*IT_0129
      *IT_0141 + 0.166666666666667*IT_0139*IT_0149 + 0.166666666666667*conj
      (IT_0095)*IT_0151 + IT_0148*IT_0152 + IT_0131*IT_0153) + 6*IT_0095*(conj
      (IT_0067)*IT_0099 + conj(IT_0082)*IT_0129 + 0.166666666666667*conj(IT_0095
      )*IT_0148 + 0.166666666666667*IT_0123*IT_0149 + IT_0151*IT_0152 + IT_0140
      *IT_0153);
    return create_ccomplex_return(IT_0154);
}

