#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_2_to_t_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_2_to_t_A(
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = V_tb*e_em*conj(U_Wm2)*U_sb_01;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = pow(m_W, -1);
    const ccomplex_t IT_0008 = m_b*conj(U_d2)*V_tb*e_em*IT_0007*U_sb_11;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0010);
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0013 = 0.666666666666667*IT_0012;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = pow(m_sb_2, 2);
    const ccomplex_t IT_0017 = cpow((-2)*s_12 + -IT_0015 + -IT_0016 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0018 = IT_0014*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0012;
    const ccomplex_t IT_0021 = pow(m_t, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0016 + IT_0021 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = s_14*s_34;
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0015 + IT_0021 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0011*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0013*IT_0028;
    const ccomplex_t IT_0030 = (-2)*IT_0024 + -IT_0029;
    const ccomplex_t IT_0031 = s_34*IT_0015;
    const ccomplex_t IT_0032 = 6*IT_0031;
    const ccomplex_t IT_0033 = s_13*s_14;
    const ccomplex_t IT_0034 = (-12)*IT_0033;
    const ccomplex_t IT_0035 = IT_0032 + IT_0034;
    const ccomplex_t IT_0036 = 0.333333333333333*IT_0012;
    const ccomplex_t IT_0037 = IT_0028*IT_0036;
    const ccomplex_t IT_0038 = (-6)*s_13;
    const ccomplex_t IT_0039 = s_14*IT_0038;
    const ccomplex_t IT_0040 = sin(beta);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0001*IT_0007*IT_0041*U_sb_01;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = IT_0026*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0036*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = m_t*m_C_2;
    const ccomplex_t IT_0049 = (-6)*s_14*IT_0048;
    const ccomplex_t IT_0050 = (-6)*IT_0024;
    const ccomplex_t IT_0051 = IT_0012*IT_0043;
    const ccomplex_t IT_0052 = IT_0022*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0013*IT_0045;
    const ccomplex_t IT_0055 = 2*IT_0053 + IT_0054;
    const ccomplex_t IT_0056 = IT_0013*IT_0043;
    const ccomplex_t IT_0057 = IT_0017*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = m_t*s_14*m_C_2;
    const ccomplex_t IT_0061 = (-6)*IT_0019;
    const ccomplex_t IT_0062 = IT_0050 + IT_0061;
    const ccomplex_t IT_0063 = -IT_0053;
    const ccomplex_t IT_0064 = (-6)*IT_0031;
    const ccomplex_t IT_0065 = (-6)*conj(IT_0024);
    const ccomplex_t IT_0066 = (-6)*conj(IT_0019);
    const ccomplex_t IT_0067 = IT_0065 + IT_0066;
    const ccomplex_t IT_0068 = IT_0015*IT_0038;
    const ccomplex_t IT_0069 = 2*IT_0019 + IT_0029;
    const ccomplex_t IT_0070 = IT_0038*IT_0048;
    const ccomplex_t IT_0071 = -IT_0054 + (-2)*IT_0058;
    const ccomplex_t IT_0072 = pow(s_13, 2);
    const ccomplex_t IT_0073 = (-6)*IT_0072;
    const ccomplex_t IT_0074 = (-6)*IT_0015*IT_0048;
    const ccomplex_t IT_0075 = conj(IT_0059) + conj(IT_0063);
    const ccomplex_t IT_0076 = (-6)*IT_0060;
    const ccomplex_t IT_0077 = m_t*s_34*m_C_2;
    const ccomplex_t IT_0078 = s_14*IT_0021;
    const ccomplex_t IT_0079 = s_13*s_34;
    const ccomplex_t IT_0080 = (-2)*IT_0079;
    const ccomplex_t IT_0081 = IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = 6*IT_0081;
    const ccomplex_t IT_0083 = s_34*IT_0038;
    const ccomplex_t IT_0084 = (-6)*s_34*IT_0048;
    const ccomplex_t IT_0085 = (-6)*IT_0078;
    const ccomplex_t IT_0086 = (-6)*IT_0021;
    const ccomplex_t IT_0087 = IT_0048*IT_0086;
    const ccomplex_t IT_0088 = 24*IT_0025;
    const ccomplex_t IT_0089 = (-24)*IT_0025;
    const ccomplex_t IT_0090 = conj(IT_0030)*IT_0076;
    const ccomplex_t IT_0091 = (-6)*IT_0077;
    const ccomplex_t IT_0092 = conj(IT_0069)*IT_0091;
    const ccomplex_t IT_0093 = 0.166666666666667*IT_0078;
    const ccomplex_t IT_0094 = 0.166666666666667*conj(IT_0055)*IT_0070 +
       0.166666666666667*conj(IT_0030)*IT_0073 + conj(IT_0024)*IT_0081 +
       0.166666666666667*conj(IT_0037)*IT_0083 + 0.166666666666667*conj(IT_0047)
      *IT_0084 + 0.166666666666667*s_13*conj(IT_0069)*IT_0086 +
       0.166666666666667*conj(IT_0071)*IT_0087 + 0.166666666666667*IT_0075
      *IT_0091 + IT_0066*IT_0093;
    const ccomplex_t IT_0095 = (IT_0024*((-24)*conj(IT_0019) + 24*conj(IT_0024
      )) + IT_0019*(24*conj(IT_0019) + (-24)*conj(IT_0024)))*IT_0025 + conj
      (IT_0030)*(IT_0019*IT_0035 + IT_0037*IT_0039 + IT_0047*IT_0049 + IT_0031
      *IT_0050) + conj(IT_0055)*(IT_0039*IT_0047 + IT_0037*IT_0049 + IT_0035
      *IT_0059 + IT_0060*IT_0062 + IT_0063*IT_0064) + IT_0055*(IT_0039*conj
      (IT_0047) + conj(IT_0037)*IT_0049 + IT_0035*conj(IT_0059) + conj(IT_0063)
      *IT_0064 + IT_0060*IT_0067 + conj(IT_0055)*IT_0068 + conj(IT_0069)*IT_0070
       + conj(IT_0071)*IT_0073 + conj(IT_0030)*IT_0074) + IT_0030*(conj(IT_0019)
      *IT_0035 + conj(IT_0037)*IT_0039 + conj(IT_0047)*IT_0049 + IT_0031*IT_0065
       + conj(IT_0030)*IT_0068 + IT_0070*conj(IT_0071) + conj(IT_0069)*IT_0073 +
       conj(IT_0055)*IT_0074 + IT_0075*IT_0076) + conj(IT_0071)*(IT_0062*IT_0077
       + IT_0063*IT_0082 + IT_0047*IT_0083 + IT_0037*IT_0084 + IT_0059*IT_0085) 
      + IT_0071*(conj(IT_0030)*IT_0070 + conj(IT_0055)*IT_0073 + IT_0067*IT_0077
       + conj(IT_0063)*IT_0082 + conj(IT_0047)*IT_0083 + conj(IT_0037)*IT_0084 +
       conj(IT_0059)*IT_0085 + s_13*conj(IT_0071)*IT_0086 + conj(IT_0069)
      *IT_0087) + IT_0063*(conj(IT_0063)*IT_0088 + conj(IT_0059)*IT_0089 +
       IT_0090 + IT_0092) + IT_0059*(conj(IT_0059)*IT_0088 + conj(IT_0063)
      *IT_0089 + IT_0090 + IT_0092) + 6*conj(IT_0069)*(IT_0024*IT_0081 +
       0.166666666666667*IT_0037*IT_0083 + 0.166666666666667*IT_0047*IT_0084 +
       IT_0061*IT_0093) + 6*IT_0069*IT_0094;
    return create_ccomplex_return(IT_0095);
}

