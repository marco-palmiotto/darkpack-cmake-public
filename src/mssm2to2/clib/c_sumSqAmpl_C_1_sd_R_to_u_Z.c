#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_sd_R_to_u_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_sd_R_to_u_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
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
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0003*IT_0008;
    const ccomplex_t IT_0010 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = IT_0004*IT_0007;
    const ccomplex_t IT_0013 = IT_0010*IT_0012;
    const ccomplex_t IT_0014 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0011 + -IT_0013 + (-2)
      *IT_0015);
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = IT_0006*IT_0017;
    const ccomplex_t IT_0019 = pow(m_u, 2);
    const ccomplex_t IT_0020 = pow(m_sd_R, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_ud*e_em*IT_0000*IT_0002*IT_0004;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = IT_0012*IT_0026;
    const ccomplex_t IT_0029 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + -IT_0028 + (-2)
      *IT_0030);
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0025*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0019 + IT_0020 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = -IT_0023 + -IT_0036;
    const ccomplex_t IT_0038 = pow(m_Z, -2);
    const ccomplex_t IT_0039 = pow(m_Z, 4);
    const ccomplex_t IT_0040 = s_13*IT_0039;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = 6*IT_0041;
    const ccomplex_t IT_0043 = s_14*s_34;
    const ccomplex_t IT_0044 = 24*IT_0043;
    const ccomplex_t IT_0045 = pow(m_Z, 2);
    const ccomplex_t IT_0046 = s_13*IT_0045;
    const ccomplex_t IT_0047 = (-12)*IT_0046;
    const ccomplex_t IT_0048 = IT_0042 + IT_0044 + IT_0047;
    const ccomplex_t IT_0049 = e_em*IT_0009;
    const ccomplex_t IT_0050 = e_em*IT_0012;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0049 + (-3)*IT_0050);
    const ccomplex_t IT_0052 = (-0.166666666666667)*IT_0051;
    const ccomplex_t IT_0053 = IT_0025*IT_0052;
    const ccomplex_t IT_0054 = pow(m_C_1, 2);
    const ccomplex_t IT_0055 = cpow((-2)*s_12 + IT_0019 + -IT_0020 + -IT_0054 
      + -reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (-24)*IT_0043;
    const ccomplex_t IT_0059 = IT_0042 + IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0008;
    const ccomplex_t IT_0061 = (-0.666666666666667)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_13 + IT_0019 + IT_0054 + -m_sd_R*
      (m_sd_R + (0 + _Complex_I*-1)*Gamma_dr) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0025*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = IT_0061*IT_0064;
    const ccomplex_t IT_0066 = 2*IT_0023 + 2*IT_0036 + -IT_0065;
    const ccomplex_t IT_0067 = s_13*s_14;
    const ccomplex_t IT_0068 = s_34*IT_0054;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0067 + IT_0069;
    const ccomplex_t IT_0071 = 6*IT_0070;
    const ccomplex_t IT_0072 = 2*IT_0057 + IT_0065;
    const ccomplex_t IT_0073 = s_14*IT_0019;
    const ccomplex_t IT_0074 = s_13*s_34;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0073 + IT_0075;
    const ccomplex_t IT_0077 = 6*IT_0076;
    const ccomplex_t IT_0078 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0079 = IT_0009*IT_0078;
    const ccomplex_t IT_0080 = IT_0012*IT_0078;
    const ccomplex_t IT_0081 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0082 = IT_0012*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + -IT_0080 + (-2)
      *IT_0082);
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = IT_0006*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0087 = IT_0021*IT_0085*IT_0086;
    const ccomplex_t IT_0088 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0089 = IT_0009*IT_0088;
    const ccomplex_t IT_0090 = IT_0012*IT_0088;
    const ccomplex_t IT_0091 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0092 = IT_0012*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + -IT_0090 + (-2)
      *IT_0092);
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0025*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0097 = IT_0034*IT_0095*IT_0096;
    const ccomplex_t IT_0098 = m_C_1*IT_0023;
    const ccomplex_t IT_0099 = m_C_1*IT_0036;
    const ccomplex_t IT_0100 = -IT_0087 + -IT_0097 + -IT_0098 + -IT_0099;
    const ccomplex_t IT_0101 = s_34*IT_0045;
    const ccomplex_t IT_0102 = IT_0038*IT_0101;
    const ccomplex_t IT_0103 = (-0.25)*IT_0102;
    const ccomplex_t IT_0104 = s_34 + IT_0103;
    const ccomplex_t IT_0105 = m_C_1*IT_0104;
    const ccomplex_t IT_0106 = (-4)*IT_0105;
    const ccomplex_t IT_0107 = 6*IT_0106;
    const ccomplex_t IT_0108 = IT_0025*IT_0061;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0110 = IT_0055*IT_0108*IT_0109;
    const ccomplex_t IT_0111 = m_u*IT_0057;
    const ccomplex_t IT_0112 = IT_0110 + -IT_0111;
    const ccomplex_t IT_0113 = s_14*IT_0045;
    const ccomplex_t IT_0114 = IT_0038*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = s_14 + IT_0115;
    const ccomplex_t IT_0117 = m_u*IT_0116;
    const ccomplex_t IT_0118 = 2*IT_0117;
    const ccomplex_t IT_0119 = 6*IT_0118;
    const ccomplex_t IT_0120 = 12*s_13;
    const ccomplex_t IT_0121 = (-2)*IT_0043;
    const ccomplex_t IT_0122 = IT_0046 + IT_0121;
    const ccomplex_t IT_0123 = IT_0038*IT_0122;
    const ccomplex_t IT_0124 = (-6)*IT_0123;
    const ccomplex_t IT_0125 = IT_0120 + IT_0124;
    const ccomplex_t IT_0126 = 0.5*IT_0102;
    const ccomplex_t IT_0127 = s_34 + IT_0126;
    const ccomplex_t IT_0128 = m_C_1*IT_0127;
    const ccomplex_t IT_0129 = 2*IT_0128;
    const ccomplex_t IT_0130 = 6*conj(IT_0057);
    const ccomplex_t IT_0131 = m_u*m_C_1;
    const ccomplex_t IT_0132 = (-18)*IT_0131;
    const ccomplex_t IT_0133 = s_14*s_34*IT_0038;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = s_13 + IT_0134;
    const ccomplex_t IT_0136 = m_C_1*IT_0135;
    const ccomplex_t IT_0137 = (-6)*IT_0136;
    const ccomplex_t IT_0138 = pow(s_34, 2);
    const ccomplex_t IT_0139 = IT_0038*IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0019 + IT_0140;
    const ccomplex_t IT_0142 = m_C_1*IT_0141;
    const ccomplex_t IT_0143 = (-6)*IT_0142;
    const ccomplex_t IT_0144 = (-0.25)*IT_0114;
    const ccomplex_t IT_0145 = s_14 + IT_0144;
    const ccomplex_t IT_0146 = m_u*IT_0145;
    const ccomplex_t IT_0147 = (-4)*IT_0146;
    const ccomplex_t IT_0148 = pow(s_14, 2);
    const ccomplex_t IT_0149 = IT_0038*IT_0148;
    const ccomplex_t IT_0150 = -IT_0054;
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = m_u*IT_0151;
    const ccomplex_t IT_0153 = 6*IT_0152;
    const ccomplex_t IT_0154 = m_u*IT_0135;
    const ccomplex_t IT_0155 = (-6)*IT_0154;
    const ccomplex_t IT_0156 = s_13*IT_0141;
    const ccomplex_t IT_0157 = (-6)*conj(IT_0057);
    const ccomplex_t IT_0158 = s_13*IT_0135;
    const ccomplex_t IT_0159 = (-6)*IT_0158;
    const ccomplex_t IT_0160 = 0.166666666666667*conj(IT_0037);
    const ccomplex_t IT_0161 = s_13*IT_0151;
    const ccomplex_t IT_0162 = IT_0037*(conj(IT_0037)*IT_0048 + conj(IT_0057)
      *IT_0059 + conj(IT_0066)*IT_0071 + conj(IT_0072)*IT_0077 + conj(IT_0100)
      *IT_0107 + conj(IT_0112)*IT_0119) + IT_0100*(conj(IT_0037)*IT_0107 + conj
      (IT_0100)*IT_0125 + IT_0129*IT_0130 + conj(IT_0112)*IT_0132 + conj(IT_0066
      )*IT_0137 + conj(IT_0072)*IT_0143) + IT_0112*(conj(IT_0037)*IT_0119 + conj
      (IT_0112)*IT_0125 + conj(IT_0100)*IT_0132 + IT_0130*IT_0147 + conj(IT_0066
      )*IT_0153 + conj(IT_0072)*IT_0155) + (-6)*IT_0072*((-0.166666666666667)
      *conj(IT_0037)*IT_0077 + (-0.166666666666667)*conj(IT_0100)*IT_0143 + (
      -0.166666666666667)*conj(IT_0112)*IT_0155 + conj(IT_0072)*IT_0156 + (
      -0.166666666666667)*IT_0076*IT_0157 + (-0.166666666666667)*conj(IT_0066)
      *IT_0159) + 6*IT_0057*(0.166666666666667*IT_0048*conj(IT_0057) + -conj
      (IT_0066)*IT_0070 + -conj(IT_0072)*IT_0076 + conj(IT_0100)*IT_0129 + conj
      (IT_0112)*IT_0147 + IT_0059*IT_0160) + 6*IT_0066*(0.166666666666667*conj
      (IT_0100)*IT_0137 + 0.166666666666667*conj(IT_0112)*IT_0153 +
       0.166666666666667*IT_0070*IT_0157 + 0.166666666666667*conj(IT_0072)
      *IT_0159 + IT_0071*IT_0160 + conj(IT_0066)*IT_0161);
    return create_ccomplex_return(IT_0162);
}

