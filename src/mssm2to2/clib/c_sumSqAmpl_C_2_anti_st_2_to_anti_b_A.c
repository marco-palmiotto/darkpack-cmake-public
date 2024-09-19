#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_st_2_to_anti_b_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_st_2_to_anti_b_A(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_C_2, 2);
    const ccomplex_t IT_0001 = s_34*IT_0000;
    const ccomplex_t IT_0002 = pow(m_b, 2);
    const ccomplex_t IT_0003 = pow(m_st_2, 2);
    const ccomplex_t IT_0004 = cpow((-2)*s_23 + IT_0002 + IT_0003 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0006 = pow(m_W, -1);
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0006*IT_0008*IT_0010*conj(U_st_01);
    const ccomplex_t IT_0012 = (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = IT_0004*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = (-1.33333333333333)*IT_0005;
    const ccomplex_t IT_0018 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = 6*IT_0021;
    const ccomplex_t IT_0023 = 2*IT_0015;
    const ccomplex_t IT_0024 = (-6)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022 + IT_0024;
    const ccomplex_t IT_0026 = 6*conj(IT_0021);
    const ccomplex_t IT_0027 = (-6)*conj(IT_0023);
    const ccomplex_t IT_0028 = IT_0026 + IT_0027;
    const ccomplex_t IT_0029 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0030 = IT_0010*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = sin(beta);
    const ccomplex_t IT_0034 = cpow(IT_0033, -1);
    const ccomplex_t IT_0035 = IT_0010*IT_0034;
    const ccomplex_t IT_0036 = m_t*V_tb*V_u2*e_em*IT_0006*conj(U_st_11);
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (-0.5)*IT_0039;
    const ccomplex_t IT_0041 = IT_0004*IT_0005;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = (IT_0032 + -IT_0040)*IT_0042;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0047 = IT_0017*IT_0040;
    const ccomplex_t IT_0048 = IT_0017*IT_0032;
    const ccomplex_t IT_0049 = IT_0046*(IT_0047 + -IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = 6*IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0043;
    const ccomplex_t IT_0054 = (-6)*IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = 6*conj(IT_0051);
    const ccomplex_t IT_0057 = (-6)*conj(IT_0053);
    const ccomplex_t IT_0058 = IT_0056 + IT_0057;
    const ccomplex_t IT_0059 = m_b*s_14*m_C_2;
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0005;
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + -IT_0000 + IT_0002 + -IT_0003 
      + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = (IT_0032 + -IT_0040)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0064;
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0012*IT_0060;
    const ccomplex_t IT_0068 = IT_0061*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = (-0.166666666666667)*IT_0028*IT_0045 + (
      -0.166666666666667)*IT_0025*conj(IT_0045) + (-0.166666666666667)*conj
      (IT_0016)*IT_0055 + (-0.166666666666667)*IT_0016*IT_0058 + (conj(IT_0021) 
      + -conj(IT_0023))*IT_0066 + (IT_0021 + -IT_0023)*conj(IT_0066) + (conj
      (IT_0051) + -conj(IT_0053))*IT_0070 + (IT_0051 + -IT_0053)*conj(IT_0070);
    const ccomplex_t IT_0072 = s_13*s_14;
    const ccomplex_t IT_0073 = (-2)*IT_0072;
    const ccomplex_t IT_0074 = IT_0001 + IT_0073;
    const ccomplex_t IT_0075 = s_14*s_34;
    const ccomplex_t IT_0076 = 24*IT_0016;
    const ccomplex_t IT_0077 = (-0.666666666666667)*IT_0005;
    const ccomplex_t IT_0078 = IT_0040*IT_0077;
    const ccomplex_t IT_0079 = IT_0032*IT_0077;
    const ccomplex_t IT_0080 = -IT_0046*(IT_0047 + -IT_0048 + -IT_0078 +
       IT_0079);
    const ccomplex_t IT_0081 = 6*s_13;
    const ccomplex_t IT_0082 = s_14*IT_0081;
    const ccomplex_t IT_0083 = (-6)*s_13;
    const ccomplex_t IT_0084 = s_14*IT_0083;
    const ccomplex_t IT_0085 = -IT_0032;
    const ccomplex_t IT_0086 = IT_0040 + IT_0085;
    const ccomplex_t IT_0087 = IT_0063*IT_0086;
    const ccomplex_t IT_0088 = -IT_0047;
    const ccomplex_t IT_0089 = IT_0048 + IT_0088;
    const ccomplex_t IT_0090 = IT_0046*IT_0089;
    const ccomplex_t IT_0091 = (-2)*IT_0087 + IT_0090;
    const ccomplex_t IT_0092 = s_34*IT_0083;
    const ccomplex_t IT_0093 = -IT_0020*(IT_0017 + -IT_0077);
    const ccomplex_t IT_0094 = -IT_0021 + (-2)*IT_0069;
    const ccomplex_t IT_0095 = IT_0000*IT_0081;
    const ccomplex_t IT_0096 = pow(s_13, 2);
    const ccomplex_t IT_0097 = 6*IT_0096;
    const ccomplex_t IT_0098 = IT_0000*IT_0083;
    const ccomplex_t IT_0099 = (-6)*IT_0002;
    const ccomplex_t IT_0100 = (-6)*IT_0096;
    const ccomplex_t IT_0101 = m_b*m_C_2;
    const ccomplex_t IT_0102 = 6*s_14;
    const ccomplex_t IT_0103 = (-6)*s_14;
    const ccomplex_t IT_0104 = (-6)*s_34;
    const ccomplex_t IT_0105 = 6*IT_0000;
    const ccomplex_t IT_0106 = (-6)*IT_0000;
    const ccomplex_t IT_0107 = conj(IT_0091)*(IT_0023*IT_0081 + IT_0021
      *IT_0083 + IT_0094*IT_0099 + IT_0093*IT_0104) + IT_0091*(conj(IT_0023)
      *IT_0081 + conj(IT_0021)*IT_0083 + conj(IT_0094)*IT_0099 + conj(IT_0093)
      *IT_0104) + conj(IT_0080)*(IT_0023*IT_0102 + IT_0021*IT_0103 + IT_0094
      *IT_0104) + IT_0080*(conj(IT_0023)*IT_0102 + conj(IT_0021)*IT_0103 + conj
      (IT_0094)*IT_0104) + conj(IT_0051)*(IT_0083*IT_0094 + IT_0093*IT_0103 +
       IT_0023*IT_0105 + IT_0021*IT_0106) + IT_0051*(IT_0083*conj(IT_0094) +
       conj(IT_0093)*IT_0103 + conj(IT_0023)*IT_0105 + conj(IT_0021)*IT_0106) +
       conj(IT_0053)*(IT_0081*IT_0094 + IT_0093*IT_0102 + IT_0021*IT_0105 +
       IT_0023*IT_0106) + IT_0053*(IT_0081*conj(IT_0094) + conj(IT_0093)*IT_0102
       + conj(IT_0021)*IT_0105 + conj(IT_0023)*IT_0106);
    const ccomplex_t IT_0108 = s_14*IT_0002;
    const ccomplex_t IT_0109 = s_13*s_34;
    const ccomplex_t IT_0110 = (-2)*IT_0109;
    const ccomplex_t IT_0111 = IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = (-6)*IT_0091;
    const ccomplex_t IT_0113 = (-6)*conj(IT_0091);
    const ccomplex_t IT_0114 = (-6)*IT_0094;
    const ccomplex_t IT_0115 = (-6)*conj(IT_0094);
    const ccomplex_t IT_0116 = m_b*s_34*m_C_2;
    const ccomplex_t IT_0117 = 6*(conj(IT_0016)*IT_0091 + IT_0016*conj(IT_0091
      ) + conj(IT_0045)*IT_0094 + IT_0045*conj(IT_0094) + 0.166666666666667*conj
      (IT_0070)*IT_0112 + 0.166666666666667*IT_0070*IT_0113 + 0.166666666666667
      *conj(IT_0066)*IT_0114 + 0.166666666666667*IT_0066*IT_0115)*IT_0116;
    const ccomplex_t IT_0118 = IT_0001*(conj(IT_0016)*IT_0025 + IT_0016
      *IT_0028 + conj(IT_0045)*IT_0055 + IT_0045*IT_0058) + (-6)*IT_0059*IT_0071
       + (IT_0058*IT_0066 + IT_0055*conj(IT_0066) + IT_0028*IT_0070 + IT_0025
      *conj(IT_0070))*IT_0074 + IT_0075*((IT_0045 + IT_0066)*(24*conj(IT_0045) +
       24*conj(IT_0066)) + (conj(IT_0016) + conj(IT_0070))*(24*IT_0070 + IT_0076
      )) + conj(IT_0080)*(IT_0053*IT_0082 + IT_0051*IT_0084 + IT_0091*IT_0092) +
       conj(IT_0093)*(IT_0023*IT_0082 + IT_0021*IT_0084 + IT_0092*IT_0094) +
       conj(IT_0023)*(IT_0082*IT_0093 + IT_0021*IT_0095 + IT_0094*IT_0097 +
       IT_0023*IT_0098) + conj(IT_0053)*(IT_0080*IT_0082 + IT_0051*IT_0095 +
       IT_0091*IT_0097 + IT_0053*IT_0098) + conj(IT_0094)*(IT_0092*IT_0093 +
       IT_0023*IT_0097 + s_13*IT_0094*IT_0099 + IT_0021*IT_0100) + conj(IT_0091)
      *(IT_0080*IT_0092 + IT_0053*IT_0097 + s_13*IT_0091*IT_0099 + IT_0051
      *IT_0100) + conj(IT_0051)*(IT_0080*IT_0084 + IT_0053*IT_0095 + IT_0051
      *IT_0098 + IT_0091*IT_0100) + conj(IT_0021)*(IT_0084*IT_0093 + IT_0023
      *IT_0095 + IT_0021*IT_0098 + IT_0094*IT_0100) + IT_0101*IT_0107 + IT_0111*
      (conj(IT_0045)*IT_0112 + IT_0045*IT_0113 + conj(IT_0016)*IT_0114 + IT_0016
      *IT_0115) + IT_0108*(conj(IT_0066)*IT_0112 + IT_0066*IT_0113 + conj
      (IT_0070)*IT_0114 + IT_0070*IT_0115) + IT_0117;
    return create_ccomplex_return(IT_0118);
}

