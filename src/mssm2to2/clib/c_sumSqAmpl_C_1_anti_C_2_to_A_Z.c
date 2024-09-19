#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_A_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0003, -1);
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = IT_0005*IT_0008;
    const ccomplex_t IT_0010 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0006 + -IT_0009 + (-2)
      *IT_0011);
    const ccomplex_t IT_0013 = (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = IT_0000*IT_0013;
    const ccomplex_t IT_0015 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = pow(m_Z, -2);
    const ccomplex_t IT_0020 = pow(s_14, 2);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = pow(m_C_1, 2);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0021 + IT_0023;
    const ccomplex_t IT_0025 = s_12*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0017;
    const ccomplex_t IT_0027 = s_12*s_14;
    const ccomplex_t IT_0028 = s_24*IT_0022;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0027 + IT_0029;
    const ccomplex_t IT_0031 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = IT_0008*IT_0031;
    const ccomplex_t IT_0034 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0035 = IT_0008*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + -IT_0033 + (-2)
      *IT_0035);
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = IT_0000*IT_0037;
    const ccomplex_t IT_0039 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = m_C_1*m_C_2;
    const ccomplex_t IT_0044 = s_14*s_24*IT_0019;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = s_12 + IT_0045;
    const ccomplex_t IT_0047 = IT_0043*IT_0046;
    const ccomplex_t IT_0048 = IT_0015*IT_0038;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = m_C_1*IT_0049;
    const ccomplex_t IT_0051 = m_C_2*IT_0017;
    const ccomplex_t IT_0052 = (-0.5)*IT_0050 + (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = pow(s_24, 2);
    const ccomplex_t IT_0054 = IT_0019*IT_0053;
    const ccomplex_t IT_0055 = pow(m_C_2, 2);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0054 + IT_0056;
    const ccomplex_t IT_0058 = m_C_1*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0041;
    const ccomplex_t IT_0060 = s_12*s_24;
    const ccomplex_t IT_0061 = s_14*IT_0055;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = IT_0060 + IT_0062;
    const ccomplex_t IT_0064 = m_C_2*IT_0046;
    const ccomplex_t IT_0065 = IT_0014*IT_0039;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = m_C_1*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = m_C_2*IT_0041;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0068 + IT_0070;
    const ccomplex_t IT_0072 = s_12*IT_0057;
    const ccomplex_t IT_0073 = conj(IT_0068) + conj(IT_0070);
    const ccomplex_t IT_0074 = s_24*m_C_1;
    const ccomplex_t IT_0075 = m_C_2*IT_0049;
    const ccomplex_t IT_0076 = m_C_1*IT_0017;
    const ccomplex_t IT_0077 = 0.5*IT_0075 + 0.5*IT_0076;
    const ccomplex_t IT_0078 = (-0.5)*IT_0066;
    const ccomplex_t IT_0079 = (-0.5)*IT_0049;
    const ccomplex_t IT_0080 = 12*IT_0052;
    const ccomplex_t IT_0081 = m_C_2*IT_0066;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = (-24)*IT_0082;
    const ccomplex_t IT_0084 = m_C_1*IT_0041;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = (-24)*IT_0085;
    const ccomplex_t IT_0087 = IT_0080 + IT_0083 + IT_0086;
    const ccomplex_t IT_0088 = 12*conj(IT_0052);
    const ccomplex_t IT_0089 = (-24)*conj(IT_0082);
    const ccomplex_t IT_0090 = (-24)*conj(IT_0085);
    const ccomplex_t IT_0091 = IT_0088 + IT_0089 + IT_0090;
    const ccomplex_t IT_0092 = (-12)*conj(IT_0026)*IT_0077 + (-12)*IT_0026
      *conj(IT_0077) + ((-24)*conj(IT_0068) + (-24)*conj(IT_0070) + 12*conj
      (IT_0077))*IT_0078 + ((-24)*IT_0068 + (-24)*IT_0070 + 12*IT_0077)*conj
      (IT_0078) + (-12)*conj(IT_0052)*IT_0079 + (-12)*IT_0052*conj(IT_0079) +
       conj(IT_0059)*IT_0087 + IT_0059*IT_0091;
    const ccomplex_t IT_0093 = pow(m_Z, 2);
    const ccomplex_t IT_0094 = s_14*IT_0093;
    const ccomplex_t IT_0095 = IT_0019*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = s_14 + IT_0096;
    const ccomplex_t IT_0098 = m_C_2*IT_0097;
    const ccomplex_t IT_0099 = 2*IT_0098;
    const ccomplex_t IT_0100 = 8*conj(IT_0052);
    const ccomplex_t IT_0101 = s_12*IT_0093;
    const ccomplex_t IT_0102 = s_14*s_24;
    const ccomplex_t IT_0103 = (-2)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = IT_0019*IT_0104;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = conj(IT_0082) + conj(IT_0085);
    const ccomplex_t IT_0108 = 8*IT_0082;
    const ccomplex_t IT_0109 = 8*conj(IT_0082);
    const ccomplex_t IT_0110 = 8*IT_0052;
    const ccomplex_t IT_0111 = 8*IT_0085;
    const ccomplex_t IT_0112 = 8*IT_0068;
    const ccomplex_t IT_0113 = 8*conj(IT_0068);
    const ccomplex_t IT_0114 = 8*IT_0077;
    const ccomplex_t IT_0115 = 8*conj(IT_0077);
    const ccomplex_t IT_0116 = s_14*m_C_2;
    const ccomplex_t IT_0117 = 12*IT_0085;
    const ccomplex_t IT_0118 = 12*conj(IT_0082);
    const ccomplex_t IT_0119 = 12*conj(IT_0085);
    const ccomplex_t IT_0120 = IT_0073*(12*IT_0059 + (-12)*IT_0079) + IT_0071*
      (12*conj(IT_0059) + (-12)*conj(IT_0079)) + ((-12)*conj(IT_0026) + 12*conj
      (IT_0078))*IT_0082 + (-12)*conj(IT_0026)*IT_0085 + IT_0026*((-12)*conj
      (IT_0082) + (-12)*conj(IT_0085)) + conj(IT_0078)*IT_0117 + IT_0078*
      (IT_0118 + IT_0119);
    const ccomplex_t IT_0121 = (-8)*IT_0042;
    const ccomplex_t IT_0122 = (-8)*conj(IT_0042);
    const ccomplex_t IT_0123 = IT_0082 + IT_0085;
    const ccomplex_t IT_0124 = m_C_1*m_C_2*IT_0093;
    const ccomplex_t IT_0125 = (-4)*IT_0102;
    const ccomplex_t IT_0126 = IT_0101 + IT_0125;
    const ccomplex_t IT_0127 = (-4)*IT_0018;
    const ccomplex_t IT_0128 = m_C_1*IT_0046;
    const ccomplex_t IT_0129 = 8*IT_0018;
    const ccomplex_t IT_0130 = 8*conj(IT_0018);
    const ccomplex_t IT_0131 = (-4)*conj(IT_0018);
    const ccomplex_t IT_0132 = IT_0043*IT_0057;
    const ccomplex_t IT_0133 = s_12*IT_0046;
    const ccomplex_t IT_0134 = IT_0047*conj(IT_0049) + (-0.5)*conj(IT_0066)
      *IT_0072 + IT_0058*(conj(IT_0068) + conj(IT_0070) + (-0.5)*conj(IT_0077)) 
      + 0.5*IT_0063*(conj(IT_0026) + -conj(IT_0078)) + IT_0064*(conj(IT_0052) + 
      -1./2*conj(IT_0082) + -1./2*conj(IT_0085)) + -conj(IT_0042)*IT_0132 + (
      -0.125)*IT_0131*IT_0133;
    const ccomplex_t IT_0135 = m_C_2*IT_0024;
    const ccomplex_t IT_0136 = IT_0043*IT_0093;
    const ccomplex_t IT_0137 = (-3)*IT_0136;
    const ccomplex_t IT_0138 = (-4)*IT_0042;
    const ccomplex_t IT_0139 = IT_0024*IT_0043;
    const ccomplex_t IT_0140 = (-4)*conj(IT_0042);
    const ccomplex_t IT_0141 = IT_0025*conj(IT_0049) + (-2)*IT_0047*conj
      (IT_0066) + -IT_0030*(conj(IT_0059) + -conj(IT_0079)) + -(conj(IT_0052) + 
      (-2)*conj(IT_0085) + (-0.25)*IT_0109)*IT_0128 + -(IT_0073 + (-0.25)
      *IT_0115)*IT_0135 + 0.25*IT_0130*IT_0139 + 0.25*IT_0133*IT_0140;
    const ccomplex_t IT_0142 = (-4)*IT_0026;
    const ccomplex_t IT_0143 = -IT_0058*(IT_0068 + IT_0070 + (-0.5)*IT_0077) +
       -IT_0064*(IT_0052 + -1./2*IT_0082 + -1./2*IT_0085) + IT_0042*IT_0132 +
       0.125*IT_0127*IT_0133 + 0.5*IT_0063*(IT_0078 + 0.25*IT_0142);
    const ccomplex_t IT_0144 = (-4)*IT_0059;
    const ccomplex_t IT_0145 = (-4)*(IT_0052 + (-0.25)*IT_0108 + (-0.25)
      *IT_0111)*IT_0128 + (-4)*(IT_0071 + (-0.25)*IT_0114)*IT_0135 + IT_0133
      *IT_0138 + IT_0129*IT_0139 + 4*IT_0030*(IT_0079 + 0.25*IT_0144);
    const ccomplex_t IT_0146 = 4*conj(IT_0018)*(IT_0018*IT_0025 + IT_0026
      *IT_0030) + (-8)*conj(IT_0042)*(IT_0018*IT_0047 + (-0.5)*IT_0052*IT_0058 +
       (-0.5)*IT_0059*IT_0063 + (-0.5)*IT_0064*IT_0071) + (-8)*IT_0042*(conj
      (IT_0018)*IT_0047 + (-0.5)*conj(IT_0052)*IT_0058 + (-0.5)*conj(IT_0059)
      *IT_0063 + (-0.5)*conj(IT_0042)*IT_0072 + (-0.5)*IT_0064*IT_0073) +
       IT_0074*IT_0092 + IT_0026*IT_0099*IT_0100 + (s_12 + IT_0106)*(IT_0107
      *IT_0108 + IT_0085*IT_0109 + conj(IT_0052)*IT_0110 + conj(IT_0085)*IT_0111
       + IT_0073*IT_0112 + 8*IT_0070*(conj(IT_0070) + 0.125*IT_0113) + conj
      (IT_0077)*IT_0114) + IT_0046*(IT_0085*IT_0100 + conj(IT_0052)*IT_0108 +
       IT_0052*IT_0109 + conj(IT_0085)*IT_0110 + conj(IT_0077)*IT_0112 + IT_0077
      *IT_0113 + conj(IT_0070)*IT_0114 + IT_0070*IT_0115) + IT_0043*(IT_0073
      *IT_0087 + IT_0071*IT_0091 + conj(IT_0077)*((-24)*IT_0052 + 12*IT_0082 +
       IT_0117) + IT_0077*((-24)*conj(IT_0052) + IT_0118 + IT_0119)) + IT_0116
      *IT_0120 + IT_0064*(conj(IT_0077)*IT_0121 + IT_0077*IT_0122) + IT_0058*
      (IT_0107*IT_0121 + IT_0122*IT_0123) + conj(IT_0078)*((-12)*IT_0026*IT_0101
       + (-24)*IT_0059*IT_0124 + (-4)*IT_0078*IT_0126 + IT_0030*IT_0127) +
       IT_0078*((-24)*conj(IT_0059)*IT_0124 + IT_0030*IT_0131) + IT_0128*(conj
      (IT_0077)*IT_0127 + IT_0073*IT_0129 + IT_0071*IT_0130 + IT_0077*IT_0131) +
       (-8)*IT_0066*IT_0134 + (IT_0107*IT_0127 + conj(IT_0052)*IT_0129 + IT_0052
      *IT_0130 + IT_0123*IT_0131)*IT_0135 + conj(IT_0079)*((-12)*IT_0059*IT_0101
       + IT_0099*IT_0114 + (-4)*IT_0079*IT_0126 + 8*IT_0026*IT_0137 + IT_0063
      *IT_0138) + IT_0079*((-12)*conj(IT_0059)*IT_0101 + IT_0099*IT_0115 +
       IT_0063*IT_0140) + 4*IT_0049*IT_0141 + 4*conj(IT_0026)*(IT_0018*IT_0030 +
       (-3)*IT_0078*IT_0101 + 1./4*IT_0099*IT_0110 + 2*IT_0079*IT_0137 + 1./4
      *IT_0126*IT_0142) + 8*conj(IT_0066)*IT_0143 + conj(IT_0059)*IT_0126
      *IT_0144 + conj(IT_0049)*IT_0145;
    return create_ccomplex_return(IT_0146);
}

