#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_A_Z(
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
    const creal_t s_34 = param->s_34;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0001, -2);
    const ccomplex_t IT_0005 = cpow(IT_0003, -2);
    const ccomplex_t IT_0006 = tan(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = cpow(1 + IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0000*IT_0002*IT_0003
      *IT_0008*(0.25*(IT_0004 + (-3)*IT_0005)*U_st_00*conj(U_st_01) + IT_0004
      *U_st_10*conj(U_st_11));
    const ccomplex_t IT_0010 = (-0.888888888888889)*IT_0009;
    const ccomplex_t IT_0011 = pow(m_st_1, 2);
    const ccomplex_t IT_0012 = pow(m_Z, -2);
    const ccomplex_t IT_0013 = pow(s_14, 2);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0011 + IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0018 = (-0.666666666666667)*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0003, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0002*IT_0003 
      + (-3)*IT_0001*IT_0019)*U_st_00*conj(U_st_01) + IT_0002*IT_0003*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0021 = 0.666666666666667*IT_0020;
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_st_2
      *Gamma_t2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_st_1
      *Gamma_t1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (-1.33333333333333)*IT_0017;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = -IT_0022;
    const ccomplex_t IT_0031 = IT_0029 + IT_0030;
    const ccomplex_t IT_0032 = IT_0027*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0025 + 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0027*IT_0029;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = 3*s_13;
    const ccomplex_t IT_0037 = IT_0035*IT_0036;
    const ccomplex_t IT_0038 = IT_0024*IT_0029;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = 3*s_23;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = IT_0037 + IT_0041;
    const ccomplex_t IT_0043 = conj(IT_0035)*IT_0036;
    const ccomplex_t IT_0044 = conj(IT_0039)*IT_0040;
    const ccomplex_t IT_0045 = IT_0043 + IT_0044;
    const ccomplex_t IT_0046 = 3*s_12;
    const ccomplex_t IT_0047 = conj(IT_0035)*IT_0046;
    const ccomplex_t IT_0048 = IT_0035*IT_0046;
    const ccomplex_t IT_0049 = pow(m_st_2, 2);
    const ccomplex_t IT_0050 = 3*IT_0049;
    const ccomplex_t IT_0051 = IT_0039*IT_0050;
    const ccomplex_t IT_0052 = 3*IT_0011;
    const ccomplex_t IT_0053 = IT_0035*IT_0052;
    const ccomplex_t IT_0054 = 3*IT_0010;
    const ccomplex_t IT_0055 = 3*conj(IT_0010);
    const ccomplex_t IT_0056 = s_14*s_24*IT_0012;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = s_12 + IT_0057;
    const ccomplex_t IT_0059 = 1.33333333333333*IT_0020;
    const ccomplex_t IT_0060 = IT_0028*IT_0059;
    const ccomplex_t IT_0061 = IT_0024*(IT_0029 + -IT_0060);
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = conj(IT_0033)*IT_0040;
    const ccomplex_t IT_0064 = conj(IT_0039)*IT_0050;
    const ccomplex_t IT_0065 = IT_0047 + IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = IT_0027*(IT_0029 + -IT_0060);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0039*IT_0046;
    const ccomplex_t IT_0069 = IT_0033*IT_0036;
    const ccomplex_t IT_0070 = IT_0053 + IT_0054 + IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = conj(IT_0039)*IT_0046;
    const ccomplex_t IT_0072 = conj(IT_0033)*IT_0036;
    const ccomplex_t IT_0073 = conj(IT_0035)*IT_0052;
    const ccomplex_t IT_0074 = IT_0055 + IT_0071 + IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = IT_0018*IT_0059;
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0022 + IT_0076;
    const ccomplex_t IT_0078 = IT_0024*IT_0077;
    const ccomplex_t IT_0079 = IT_0029 + IT_0075;
    const ccomplex_t IT_0080 = -IT_0022 + -IT_0060;
    const ccomplex_t IT_0081 = IT_0079 + IT_0080;
    const ccomplex_t IT_0082 = IT_0027*IT_0081;
    const ccomplex_t IT_0083 = (-0.5)*IT_0078 + 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0033*IT_0040;
    const ccomplex_t IT_0085 = s_14*s_34*IT_0012;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = s_13 + IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0061;
    const ccomplex_t IT_0089 = 0.5*IT_0034;
    const ccomplex_t IT_0090 = (-0.5)*IT_0032 + 0.5*IT_0078;
    const ccomplex_t IT_0091 = s_24*s_34*IT_0012;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = s_23 + IT_0092;
    const ccomplex_t IT_0094 = IT_0036*IT_0067;
    const ccomplex_t IT_0095 = IT_0040*IT_0062;
    const ccomplex_t IT_0096 = IT_0094 + IT_0095;
    const ccomplex_t IT_0097 = IT_0036*conj(IT_0067);
    const ccomplex_t IT_0098 = IT_0040*conj(IT_0062);
    const ccomplex_t IT_0099 = IT_0097 + IT_0098;
    const ccomplex_t IT_0100 = IT_0052*IT_0067;
    const ccomplex_t IT_0101 = IT_0046*IT_0067;
    const ccomplex_t IT_0102 = IT_0050*IT_0062;
    const ccomplex_t IT_0103 = IT_0046*conj(IT_0067);
    const ccomplex_t IT_0104 = IT_0055*IT_0083 + IT_0054*conj(IT_0083) + 
      (IT_0046*conj(IT_0062) + IT_0052*conj(IT_0067) + IT_0036*conj(IT_0083))
      *IT_0089 + conj(IT_0090)*IT_0096 + IT_0090*IT_0099 + conj(IT_0089)*
      (IT_0046*IT_0062 + IT_0036*IT_0083 + IT_0100) + conj(IT_0088)*(IT_0054 +
       IT_0040*IT_0083 + IT_0101 + IT_0102) + IT_0088*(IT_0055 + IT_0050*conj
      (IT_0062) + IT_0040*conj(IT_0083) + IT_0103);
    const ccomplex_t IT_0105 = pow(s_24, 2);
    const ccomplex_t IT_0106 = IT_0012*IT_0105;
    const ccomplex_t IT_0107 = -IT_0106;
    const ccomplex_t IT_0108 = IT_0049 + IT_0107;
    const ccomplex_t IT_0109 = pow(s_34, 2);
    const ccomplex_t IT_0110 = IT_0012*IT_0109;
    const ccomplex_t IT_0111 = (-3)*s_13;
    const ccomplex_t IT_0112 = (-3)*s_23;
    const ccomplex_t IT_0113 = (-0.333333333333333)*conj(IT_0090);
    const ccomplex_t IT_0114 = IT_0010*conj(IT_0090) + IT_0090*(conj(IT_0010) 
      + (-0.333333333333333)*conj(IT_0089)*IT_0111 + (-0.333333333333333)*conj
      (IT_0088)*IT_0112) + IT_0089*(s_12*conj(IT_0088) + IT_0011*conj(IT_0089) +
       IT_0111*IT_0113) + IT_0088*(IT_0049*conj(IT_0088) + s_12*conj(IT_0089) +
       IT_0112*IT_0113);
    const ccomplex_t IT_0115 = 9*IT_0010*conj(IT_0010) + IT_0016*(conj(IT_0033
      )*IT_0042 + IT_0033*IT_0045 + IT_0039*IT_0047 + conj(IT_0039)*(IT_0048 +
       IT_0051) + conj(IT_0035)*(IT_0053 + IT_0054) + IT_0035*IT_0055) + IT_0058
      *(conj(IT_0039)*IT_0054 + IT_0039*IT_0055 + IT_0062*IT_0065 + conj(IT_0067
      )*IT_0070 + IT_0067*IT_0074 + IT_0045*IT_0083 + IT_0042*conj(IT_0083) +
       conj(IT_0062)*(IT_0048 + IT_0051 + IT_0084)) + IT_0087*(conj(IT_0033)
      *IT_0054 + IT_0033*IT_0055 + IT_0065*IT_0088 + (IT_0048 + IT_0051 +
       IT_0084)*conj(IT_0088) + IT_0074*IT_0089 + IT_0070*conj(IT_0089) +
       IT_0045*IT_0090 + IT_0042*conj(IT_0090)) + IT_0093*IT_0104 + (conj
      (IT_0083)*IT_0096 + IT_0083*IT_0099 + conj(IT_0067)*IT_0100 + conj(IT_0062
      )*(IT_0054 + IT_0101 + IT_0102) + IT_0062*(IT_0055 + IT_0103))*IT_0108 + (
      -3)*IT_0110*IT_0114;
    return create_ccomplex_return(IT_0115);
}

