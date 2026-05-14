#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_A_Z(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0004, -1);
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + -IT_0010 + (-2)
      *IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = IT_0000*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_2, 2);
    const ccomplex_t IT_0019 = pow(m_Z, -2);
    const ccomplex_t IT_0020 = pow(s_24, 2);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = IT_0018 + IT_0022;
    const ccomplex_t IT_0024 = IT_0018*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0027 = IT_0005*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0026;
    const ccomplex_t IT_0029 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0030 = IT_0009*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0027 + -IT_0028 + (-2)
      *IT_0030);
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = IT_0000*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = s_12*IT_0023;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = 4*IT_0017;
    const ccomplex_t IT_0040 = pow(s_14, 2);
    const ccomplex_t IT_0041 = IT_0019*IT_0040;
    const ccomplex_t IT_0042 = -IT_0018;
    const ccomplex_t IT_0043 = IT_0041 + IT_0042;
    const ccomplex_t IT_0044 = m_C_2*IT_0043;
    const ccomplex_t IT_0045 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0015*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0033*IT_0045;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0047 + IT_0049;
    const ccomplex_t IT_0051 = m_C_2*IT_0050;
    const ccomplex_t IT_0052 = 0.5*IT_0051;
    const ccomplex_t IT_0053 = 8*conj(IT_0049);
    const ccomplex_t IT_0054 = s_12*IT_0043;
    const ccomplex_t IT_0055 = (-0.5)*IT_0049;
    const ccomplex_t IT_0056 = s_12*s_14;
    const ccomplex_t IT_0057 = s_24*IT_0018;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = s_14*s_24*IT_0019;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = s_12 + IT_0061;
    const ccomplex_t IT_0063 = IT_0018*IT_0062;
    const ccomplex_t IT_0064 = (-8)*IT_0017;
    const ccomplex_t IT_0065 = (-8)*conj(IT_0017);
    const ccomplex_t IT_0066 = m_C_2*IT_0023;
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = 4*conj(IT_0017);
    const ccomplex_t IT_0069 = IT_0017 + IT_0035;
    const ccomplex_t IT_0070 = m_C_2*IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = -IT_0047;
    const ccomplex_t IT_0073 = IT_0018*IT_0043;
    const ccomplex_t IT_0074 = s_12*IT_0062;
    const ccomplex_t IT_0075 = (-4)*conj(IT_0017);
    const ccomplex_t IT_0076 = s_12*s_24;
    const ccomplex_t IT_0077 = s_14*IT_0018;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = IT_0076 + IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0017;
    const ccomplex_t IT_0081 = 0.5*IT_0047;
    const ccomplex_t IT_0082 = m_C_2*IT_0062;
    const ccomplex_t IT_0083 = 4*IT_0082;
    const ccomplex_t IT_0084 = 8*IT_0062;
    const ccomplex_t IT_0085 = 8*IT_0082;
    const ccomplex_t IT_0086 = pow(m_Z, 2);
    const ccomplex_t IT_0087 = s_12*IT_0086;
    const ccomplex_t IT_0088 = s_14*s_24;
    const ccomplex_t IT_0089 = (-2)*IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = IT_0019*IT_0090;
    const ccomplex_t IT_0092 = (-0.5)*IT_0091;
    const ccomplex_t IT_0093 = s_12 + IT_0092;
    const ccomplex_t IT_0094 = 8*IT_0093;
    const ccomplex_t IT_0095 = (-0.5)*IT_0051;
    const ccomplex_t IT_0096 = 12*IT_0018;
    const ccomplex_t IT_0097 = 0.5*IT_0035;
    const ccomplex_t IT_0098 = s_14*m_C_2;
    const ccomplex_t IT_0099 = 12*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0070;
    const ccomplex_t IT_0101 = (-24)*IT_0018;
    const ccomplex_t IT_0102 = s_24*m_C_2;
    const ccomplex_t IT_0103 = (-24)*IT_0102;
    const ccomplex_t IT_0104 = (-12)*IT_0098;
    const ccomplex_t IT_0105 = (-4)*IT_0049;
    const ccomplex_t IT_0106 = IT_0064*IT_0067 + IT_0036*IT_0083 + IT_0052
      *IT_0084 + IT_0072*IT_0085 + IT_0071*IT_0094 + IT_0095*IT_0096 + IT_0097
      *IT_0099 + IT_0100*IT_0101 + IT_0080*IT_0103 + IT_0055*IT_0104 + IT_0044
      *IT_0105;
    const ccomplex_t IT_0107 = 4*IT_0067;
    const ccomplex_t IT_0108 = s_14*IT_0086;
    const ccomplex_t IT_0109 = IT_0019*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = s_14 + IT_0110;
    const ccomplex_t IT_0112 = m_C_2*IT_0111;
    const ccomplex_t IT_0113 = 2*IT_0112;
    const ccomplex_t IT_0114 = 8*IT_0113;
    const ccomplex_t IT_0115 = 8*IT_0044;
    const ccomplex_t IT_0116 = 12*IT_0102;
    const ccomplex_t IT_0117 = (-12)*IT_0102;
    const ccomplex_t IT_0118 = IT_0094*IT_0095 + IT_0071*IT_0096 + IT_0084
      *IT_0100 + IT_0052*IT_0101 + IT_0082*(IT_0064 + IT_0105) + IT_0036*IT_0107
       + IT_0081*IT_0114 + IT_0072*IT_0115 + IT_0097*IT_0116 + IT_0055*IT_0117;
    const ccomplex_t IT_0119 = (-4)*conj(IT_0049);
    const ccomplex_t IT_0120 = 4*IT_0079;
    const ccomplex_t IT_0121 = (-8)*IT_0063;
    const ccomplex_t IT_0122 = (-8)*IT_0082;
    const ccomplex_t IT_0123 = (-8)*IT_0067;
    const ccomplex_t IT_0124 = (-4)*IT_0079;
    const ccomplex_t IT_0125 = IT_0017*IT_0025 + 0.5*IT_0036*IT_0038 + 0.125
      *IT_0071*IT_0083 + 0.125*IT_0074*IT_0105 + 0.125*IT_0095*IT_0107 + 0.125
      *IT_0097*IT_0120 + 0.125*IT_0072*IT_0121 + 0.125*IT_0052*IT_0122 + 0.125
      *IT_0100*IT_0123 + 0.125*IT_0055*IT_0124;
    const ccomplex_t IT_0126 = IT_0018*IT_0086;
    const ccomplex_t IT_0127 = (-24)*IT_0126;
    const ccomplex_t IT_0128 = (-12)*IT_0087;
    const ccomplex_t IT_0129 = (-4)*IT_0059;
    const ccomplex_t IT_0130 = (-4)*IT_0088;
    const ccomplex_t IT_0131 = IT_0087 + IT_0130;
    const ccomplex_t IT_0132 = (-4)*IT_0131;
    const ccomplex_t IT_0133 = (-3)*IT_0126;
    const ccomplex_t IT_0134 = 8*IT_0133;
    const ccomplex_t IT_0135 = IT_0049*IT_0059 + 0.25*IT_0071*IT_0104 + 0.25
      *IT_0052*IT_0114 + 0.25*IT_0095*IT_0117 + 0.25*IT_0036*IT_0124 + 0.25
      *IT_0097*IT_0128 + 0.25*IT_0055*IT_0132 + 0.25*IT_0081*IT_0134;
    const ccomplex_t IT_0136 = 8*IT_0049;
    const ccomplex_t IT_0137 = (-4)*IT_0082;
    const ccomplex_t IT_0138 = IT_0039*IT_0067 + IT_0071*IT_0084 + IT_0052
      *IT_0094 + IT_0096*IT_0100 + IT_0095*IT_0101 + IT_0055*IT_0114 + IT_0080
      *IT_0116 + IT_0081*IT_0117 + IT_0036*IT_0122 + IT_0044*IT_0136 + IT_0072
      *IT_0137;
    const ccomplex_t IT_0139 = (-4)*IT_0044;
    const ccomplex_t IT_0140 = IT_0084*IT_0095 + IT_0052*IT_0096 + IT_0080
      *IT_0099 + IT_0094*IT_0100 + IT_0071*IT_0101 + IT_0097*IT_0103 + IT_0081
      *IT_0104 + IT_0036*IT_0123 + IT_0082*(IT_0039 + IT_0136) + IT_0072*IT_0139;
    const ccomplex_t IT_0141 = (-0.25)*IT_0100;
    const ccomplex_t IT_0142 = (-0.25)*IT_0095;
    const ccomplex_t IT_0143 = (-0.25)*IT_0052;
    const ccomplex_t IT_0144 = (-0.25)*IT_0080;
    const ccomplex_t IT_0145 = 4*IT_0059;
    const ccomplex_t IT_0146 = (-0.25)*IT_0081;
    const ccomplex_t IT_0147 = -IT_0054*IT_0072 + IT_0017*IT_0074 + (-0.25)
      *IT_0071*IT_0085 + (-0.25)*IT_0036*IT_0121 + (-0.25)*IT_0073*IT_0136 +
       IT_0139*IT_0141 + IT_0115*IT_0142 + IT_0137*IT_0143 + IT_0129*IT_0144 +
       IT_0145*IT_0146;
    const ccomplex_t IT_0148 = 8*conj(IT_0017)*(IT_0025*IT_0036 + 0.125
      *IT_0038*IT_0039) + IT_0044*IT_0052*IT_0053 + 4*conj(IT_0049)*(IT_0049
      *IT_0054 + IT_0055*IT_0059 + 0.25*IT_0063*IT_0064) + IT_0049*IT_0063
      *IT_0065 + IT_0067*(IT_0052*IT_0068 + IT_0065*IT_0071) + IT_0072*(IT_0053
      *IT_0073 + IT_0074*IT_0075) + IT_0079*(IT_0068*IT_0080 + IT_0075*IT_0081) 
      + conj(IT_0071)*IT_0106 + conj(IT_0095)*IT_0118 + (IT_0044*IT_0071 +
       IT_0036*IT_0074 + IT_0059*IT_0097)*IT_0119 + IT_0082*((IT_0053 + IT_0068)
      *IT_0100 + IT_0095*(IT_0065 + IT_0119)) + 8*conj(IT_0036)*IT_0125 + conj
      (IT_0080)*(IT_0039*IT_0079 + IT_0099*IT_0100 + IT_0071*IT_0103 + IT_0052
      *IT_0116 + IT_0097*IT_0127 + IT_0081*IT_0128 + IT_0072*IT_0129 + IT_0080
      *IT_0132) + conj(IT_0097)*(IT_0071*IT_0099 + IT_0100*IT_0103 + IT_0059
      *IT_0105 + IT_0095*IT_0116 + IT_0036*IT_0120 + IT_0080*IT_0127 + IT_0055
      *IT_0128 + IT_0097*IT_0132) + 4*conj(IT_0055)*IT_0135 + conj(IT_0052)
      *IT_0138 + conj(IT_0100)*IT_0140 + (-4)*conj(IT_0081)*(IT_0017*IT_0079 + (
      -0.25)*IT_0055*IT_0134 + IT_0104*IT_0141 + IT_0114*IT_0142 + IT_0117
      *IT_0143 + IT_0128*IT_0144 + (-0.25)*IT_0072*IT_0145 + IT_0132*IT_0146) + 
      (-4)*conj(IT_0072)*IT_0147;
    return create_ccomplex_return(IT_0148);
}

