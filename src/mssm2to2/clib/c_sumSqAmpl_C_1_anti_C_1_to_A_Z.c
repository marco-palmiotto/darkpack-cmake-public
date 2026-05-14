#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_A_Z(
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
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0004, -1);
    const ccomplex_t IT_0009 = IT_0002*IT_0008;
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*(IT_0007 + -IT_0010 + (-2)
      *IT_0012);
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = IT_0000*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0015*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0024 = IT_0005*IT_0023;
    const ccomplex_t IT_0025 = IT_0009*IT_0023;
    const ccomplex_t IT_0026 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + -IT_0025 + (-2)
      *IT_0027);
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = IT_0000*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0019*IT_0030;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = pow(m_Z, 2);
    const ccomplex_t IT_0038 = s_12*IT_0037;
    const ccomplex_t IT_0039 = pow(m_C_1, 2);
    const ccomplex_t IT_0040 = IT_0037*IT_0039;
    const ccomplex_t IT_0041 = -IT_0021;
    const ccomplex_t IT_0042 = pow(m_Z, -2);
    const ccomplex_t IT_0043 = pow(s_14, 2);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = -IT_0039;
    const ccomplex_t IT_0046 = IT_0044 + IT_0045;
    const ccomplex_t IT_0047 = s_12*IT_0046;
    const ccomplex_t IT_0048 = -IT_0032;
    const ccomplex_t IT_0049 = s_14*s_24*IT_0042;
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = s_12 + IT_0050;
    const ccomplex_t IT_0052 = IT_0039*IT_0051;
    const ccomplex_t IT_0053 = pow(s_24, 2);
    const ccomplex_t IT_0054 = IT_0042*IT_0053;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0039 + IT_0055;
    const ccomplex_t IT_0057 = IT_0039*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = (-3)*IT_0040;
    const ccomplex_t IT_0060 = s_14*m_C_1;
    const ccomplex_t IT_0061 = m_C_1*IT_0017;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = m_C_1*IT_0032;
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0062 + IT_0064;
    const ccomplex_t IT_0066 = conj(IT_0062) + conj(IT_0064);
    const ccomplex_t IT_0067 = (-0.5)*IT_0063;
    const ccomplex_t IT_0068 = (-0.5)*IT_0061;
    const ccomplex_t IT_0069 = IT_0067 + IT_0068;
    const ccomplex_t IT_0070 = conj(IT_0067) + conj(IT_0068);
    const ccomplex_t IT_0071 = IT_0060*((12*conj(IT_0018) + (-12)*conj(IT_0022
      ))*IT_0065 + (12*IT_0018 + (-12)*IT_0022)*IT_0066 + (12*conj(IT_0033) + (
      -12)*conj(IT_0036))*IT_0069 + (12*IT_0033 + (-12)*IT_0036)*IT_0070);
    const ccomplex_t IT_0072 = s_12*IT_0056;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = 4*IT_0048;
    const ccomplex_t IT_0075 = m_C_1*IT_0035;
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = 8*IT_0076;
    const ccomplex_t IT_0078 = 8*conj(IT_0076);
    const ccomplex_t IT_0079 = m_C_1*IT_0021;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = 8*IT_0062;
    const ccomplex_t IT_0082 = 8*conj(IT_0062);
    const ccomplex_t IT_0083 = 8*IT_0080;
    const ccomplex_t IT_0084 = 8*conj(IT_0080);
    const ccomplex_t IT_0085 = 0.5*IT_0079;
    const ccomplex_t IT_0086 = 0.5*IT_0075;
    const ccomplex_t IT_0087 = conj(IT_0085) + conj(IT_0086);
    const ccomplex_t IT_0088 = 8*IT_0067;
    const ccomplex_t IT_0089 = IT_0085 + IT_0086;
    const ccomplex_t IT_0090 = 8*conj(IT_0067);
    const ccomplex_t IT_0091 = 8*IT_0085;
    const ccomplex_t IT_0092 = 8*conj(IT_0085);
    const ccomplex_t IT_0093 = 8*IT_0068;
    const ccomplex_t IT_0094 = IT_0066*IT_0077 + IT_0065*IT_0078 + conj
      (IT_0080)*IT_0081 + IT_0080*IT_0082 + conj(IT_0064)*IT_0083 + IT_0064
      *IT_0084 + IT_0087*IT_0088 + IT_0089*IT_0090 + 8*conj(IT_0068)*(IT_0086 +
       0.125*IT_0091) + IT_0068*IT_0092 + conj(IT_0086)*IT_0093;
    const ccomplex_t IT_0095 = s_14*s_24;
    const ccomplex_t IT_0096 = (-2)*IT_0095;
    const ccomplex_t IT_0097 = IT_0038 + IT_0096;
    const ccomplex_t IT_0098 = IT_0042*IT_0097;
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = conj(IT_0076) + conj(IT_0080);
    const ccomplex_t IT_0101 = 8*IT_0064;
    const ccomplex_t IT_0102 = 8*IT_0086;
    const ccomplex_t IT_0103 = IT_0078*IT_0080 + IT_0066*IT_0081 + IT_0064
      *IT_0082 + conj(IT_0080)*IT_0083 + IT_0070*IT_0088 + IT_0068*IT_0090 +
       IT_0087*IT_0091 + IT_0086*IT_0092 + conj(IT_0068)*IT_0093 + IT_0077
      *IT_0100 + conj(IT_0064)*IT_0101 + conj(IT_0086)*IT_0102;
    const ccomplex_t IT_0104 = s_14*IT_0037;
    const ccomplex_t IT_0105 = IT_0042*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = s_14 + IT_0106;
    const ccomplex_t IT_0108 = m_C_1*IT_0107;
    const ccomplex_t IT_0109 = 2*IT_0108;
    const ccomplex_t IT_0110 = 8*conj(IT_0086);
    const ccomplex_t IT_0111 = s_12*IT_0051;
    const ccomplex_t IT_0112 = (-4)*IT_0041;
    const ccomplex_t IT_0113 = (-4)*conj(IT_0041);
    const ccomplex_t IT_0114 = m_C_1*IT_0046;
    const ccomplex_t IT_0115 = 8*IT_0041;
    const ccomplex_t IT_0116 = IT_0076 + IT_0080;
    const ccomplex_t IT_0117 = 8*conj(IT_0041);
    const ccomplex_t IT_0118 = (-4)*conj(IT_0035)*(IT_0069 + (-0.25)*IT_0091 +
       (-0.25)*IT_0102) + (-4)*IT_0035*(IT_0070 + (-0.25)*IT_0092 + (-0.25)
      *IT_0110) + IT_0066*IT_0112 + IT_0065*IT_0113 + IT_0100*IT_0115 + IT_0116
      *IT_0117;
    const ccomplex_t IT_0119 = IT_0039*IT_0046;
    const ccomplex_t IT_0120 = (-4)*IT_0048;
    const ccomplex_t IT_0121 = (-4)*conj(IT_0048);
    const ccomplex_t IT_0122 = s_12*s_24;
    const ccomplex_t IT_0123 = s_14*IT_0039;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = 4*conj(IT_0048);
    const ccomplex_t IT_0127 = (-4)*IT_0022;
    const ccomplex_t IT_0128 = 4*IT_0125*(IT_0017*(conj(IT_0018) + -conj
      (IT_0022)) + 0.25*conj(IT_0033)*IT_0074 + 0.25*conj(IT_0036)*IT_0120 +
       0.25*IT_0036*IT_0121 + 0.25*IT_0033*IT_0126 + conj(IT_0017)*(IT_0018 +
       0.25*IT_0127));
    const ccomplex_t IT_0129 = (-4)*IT_0095;
    const ccomplex_t IT_0130 = IT_0038 + IT_0129;
    const ccomplex_t IT_0131 = (-4)*IT_0033;
    const ccomplex_t IT_0132 = s_12*s_14;
    const ccomplex_t IT_0133 = s_24*IT_0039;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = IT_0132 + IT_0134;
    const ccomplex_t IT_0136 = (-4)*(IT_0035*(conj(IT_0033) + -conj(IT_0036)) 
      + -conj(IT_0022)*IT_0041 + -IT_0022*conj(IT_0041) + (-0.25)*conj(IT_0018)
      *IT_0112 + (-0.25)*IT_0018*IT_0113 + -conj(IT_0035)*(IT_0036 + 0.25
      *IT_0131))*IT_0135;
    const ccomplex_t IT_0137 = (-24)*conj(IT_0062);
    const ccomplex_t IT_0138 = (-24)*conj(IT_0064);
    const ccomplex_t IT_0139 = 12*conj(IT_0080);
    const ccomplex_t IT_0140 = 12*conj(IT_0076);
    const ccomplex_t IT_0141 = IT_0137 + IT_0138 + IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = (-24)*IT_0080;
    const ccomplex_t IT_0143 = (-24)*conj(IT_0076);
    const ccomplex_t IT_0144 = (-24)*conj(IT_0080);
    const ccomplex_t IT_0145 = (-24)*IT_0062;
    const ccomplex_t IT_0146 = (-24)*IT_0064;
    const ccomplex_t IT_0147 = 12*IT_0080;
    const ccomplex_t IT_0148 = IT_0076*(12*conj(IT_0067) + 12*conj(IT_0068) + 
      (-24)*conj(IT_0085) + (-24)*conj(IT_0086)) + IT_0069*IT_0141 + IT_0087*(12
      *IT_0062 + 12*IT_0064 + IT_0142) + IT_0089*(12*conj(IT_0062) + 12*conj
      (IT_0064) + IT_0143 + IT_0144) + IT_0070*(IT_0145 + IT_0146 + IT_0147);
    const ccomplex_t IT_0149 = s_24*m_C_1;
    const ccomplex_t IT_0150 = conj(IT_0036)*((-12)*IT_0076 + (-12)*IT_0080) +
       IT_0036*((-12)*conj(IT_0076) + (-12)*conj(IT_0080)) + conj(IT_0018)*((-24
      )*IT_0067 + (-24)*IT_0068 + 12*IT_0085 + 12*IT_0086) + conj(IT_0022)*((-12
      )*IT_0085 + (-12)*IT_0086) + IT_0018*((-24)*conj(IT_0067) + (-24)*conj
      (IT_0068) + 12*conj(IT_0085) + 12*conj(IT_0086)) + IT_0022*((-12)*conj
      (IT_0085) + (-12)*conj(IT_0086)) + IT_0033*IT_0141 + conj(IT_0033)*(12
      *IT_0076 + IT_0145 + IT_0146 + IT_0147);
    const ccomplex_t IT_0151 = m_C_1*IT_0056;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = (-2)*conj(IT_0064);
    const ccomplex_t IT_0154 = (-2)*conj(IT_0048)*(IT_0062 + IT_0064 + -1./2
      *IT_0076 + -1./2*IT_0080) + (-2)*conj(IT_0017)*(IT_0067 + IT_0068 + (-0.5)
      *IT_0085 + (-0.5)*IT_0086) + (-2)*IT_0017*(conj(IT_0067) + conj(IT_0068) +
       (-0.5)*conj(IT_0085) + (-0.5)*conj(IT_0086)) + (-2)*IT_0048*(conj(IT_0062
      ) + -1./2*conj(IT_0076) + -1./2*conj(IT_0080) + -1./2*IT_0153);
    const ccomplex_t IT_0155 = m_C_1*IT_0051;
    const ccomplex_t IT_0156 = (-8)*IT_0048;
    const ccomplex_t IT_0157 = (-8)*conj(IT_0048);
    const ccomplex_t IT_0158 = conj(IT_0017)*((-0.5)*IT_0062 + (-0.5)*IT_0064 
      + IT_0076 + IT_0080) + IT_0017*((-0.5)*conj(IT_0062) + (-0.5)*conj(IT_0064
      ) + conj(IT_0076) + conj(IT_0080)) + (-0.125)*IT_0070*(IT_0074 + IT_0115) 
      + 0.5*conj(IT_0035)*((-0.25)*IT_0081 + (-0.25)*IT_0101 + IT_0116) + (
      -0.125)*IT_0069*(IT_0117 + IT_0126) + 0.5*IT_0035*((-0.25)*IT_0082 +
       IT_0100 + IT_0153) + (-0.125)*IT_0087*(IT_0112 + IT_0156) + (-0.125)
      *IT_0089*(IT_0113 + IT_0157);
    const ccomplex_t IT_0159 = ((-12)*conj(IT_0018)*IT_0022 + (-12)*IT_0018
      *conj(IT_0022) + (-12)*conj(IT_0033)*IT_0036 + (-12)*IT_0033*conj(IT_0036)
      )*IT_0038 + ((-24)*conj(IT_0018)*IT_0033 + (-24)*IT_0018*conj(IT_0033))
      *IT_0040 + 4*IT_0041*(conj(IT_0041)*IT_0047 + (-2)*conj(IT_0048)*IT_0052) 
      + (-8)*IT_0048*(conj(IT_0041)*IT_0052 + -conj(IT_0017)*IT_0058) + 8*(conj
      (IT_0022)*IT_0036 + IT_0022*conj(IT_0036))*IT_0059 + IT_0071 + conj
      (IT_0048)*IT_0073*IT_0074 + IT_0051*IT_0094 + (s_12 + IT_0099)*IT_0103 +
       IT_0109*(conj(IT_0022)*(IT_0077 + IT_0083) + IT_0022*(IT_0078 + IT_0084) 
      + conj(IT_0036)*(IT_0091 + IT_0102) + IT_0036*(IT_0092 + IT_0110)) + conj
      (IT_0017)*IT_0111*IT_0112 + (-8)*IT_0017*(conj(IT_0035)*IT_0052 + -conj
      (IT_0048)*IT_0058 + (-0.5)*conj(IT_0017)*IT_0073 + (-0.125)*IT_0111
      *IT_0113) + IT_0114*IT_0118 + conj(IT_0035)*(IT_0115*IT_0119 + IT_0111
      *IT_0120) + 4*IT_0035*(conj(IT_0035)*IT_0047 + (-2)*conj(IT_0017)*IT_0052 
      + 0.25*IT_0117*IT_0119 + 0.25*IT_0111*IT_0121) + IT_0128 + (-4)*IT_0130*
      (IT_0018*conj(IT_0018) + IT_0036*conj(IT_0036) + (-0.25)*conj(IT_0022)
      *IT_0127 + (-0.25)*conj(IT_0033)*IT_0131) + IT_0136 + IT_0039*IT_0148 +
       IT_0149*IT_0150 + 4*IT_0152*IT_0154 + (-8)*IT_0155*IT_0158;
    return create_ccomplex_return(IT_0159);
}

