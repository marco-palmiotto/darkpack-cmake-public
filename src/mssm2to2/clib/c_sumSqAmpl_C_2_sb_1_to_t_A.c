#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_sb_1_to_t_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_sb_1_to_t_A(
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
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = pow(m_t, 2);
    const ccomplex_t IT_0001 = pow(m_sb_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = m_b*conj(U_d2)*V_tb*e_em*IT_0009*U_sb_10;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = V_tb*e_em*conj(U_Wm2)*U_sb_00;
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0012 + (-2)*IT_0014);
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_C_2, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + -IT_0001 + -IT_0020 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0022 = 0.666666666666667*IT_0003;
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = s_14*s_34;
    const ccomplex_t IT_0027 = cpow((-2)*s_13 + IT_0000 + IT_0020 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0016*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = (-2)*IT_0019 + -IT_0030;
    const ccomplex_t IT_0032 = s_34*IT_0020;
    const ccomplex_t IT_0033 = 6*IT_0032;
    const ccomplex_t IT_0034 = s_13*s_14;
    const ccomplex_t IT_0035 = (-12)*IT_0034;
    const ccomplex_t IT_0036 = IT_0033 + IT_0035;
    const ccomplex_t IT_0037 = 0.333333333333333*IT_0003;
    const ccomplex_t IT_0038 = IT_0029*IT_0037;
    const ccomplex_t IT_0039 = (-6)*s_13;
    const ccomplex_t IT_0040 = s_14*IT_0039;
    const ccomplex_t IT_0041 = sin(beta);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*m_t*V_tb*V_u2
      *e_em*IT_0007*IT_0009*IT_0042*U_sb_00;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0027*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0037*IT_0046;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = m_t*m_C_2;
    const ccomplex_t IT_0050 = (-6)*s_14*IT_0049;
    const ccomplex_t IT_0051 = (-6)*IT_0019;
    const ccomplex_t IT_0052 = IT_0022*IT_0046;
    const ccomplex_t IT_0053 = IT_0003*IT_0044;
    const ccomplex_t IT_0054 = IT_0002*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0052 + 2*IT_0055;
    const ccomplex_t IT_0057 = IT_0022*IT_0044;
    const ccomplex_t IT_0058 = IT_0021*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = m_t*s_14*m_C_2;
    const ccomplex_t IT_0062 = (-6)*IT_0025;
    const ccomplex_t IT_0063 = IT_0051 + IT_0062;
    const ccomplex_t IT_0064 = -IT_0055;
    const ccomplex_t IT_0065 = (-6)*IT_0032;
    const ccomplex_t IT_0066 = (-6)*conj(IT_0025);
    const ccomplex_t IT_0067 = (-6)*conj(IT_0019);
    const ccomplex_t IT_0068 = IT_0066 + IT_0067;
    const ccomplex_t IT_0069 = IT_0020*IT_0039;
    const ccomplex_t IT_0070 = 2*IT_0025 + IT_0030;
    const ccomplex_t IT_0071 = IT_0039*IT_0049;
    const ccomplex_t IT_0072 = -IT_0052 + (-2)*IT_0059;
    const ccomplex_t IT_0073 = pow(s_13, 2);
    const ccomplex_t IT_0074 = (-6)*IT_0073;
    const ccomplex_t IT_0075 = (-6)*IT_0020*IT_0049;
    const ccomplex_t IT_0076 = conj(IT_0060) + conj(IT_0064);
    const ccomplex_t IT_0077 = (-6)*IT_0061;
    const ccomplex_t IT_0078 = m_t*s_34*m_C_2;
    const ccomplex_t IT_0079 = s_14*IT_0000;
    const ccomplex_t IT_0080 = s_13*s_34;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = 6*IT_0082;
    const ccomplex_t IT_0084 = s_34*IT_0039;
    const ccomplex_t IT_0085 = (-6)*s_34*IT_0049;
    const ccomplex_t IT_0086 = (-6)*IT_0079;
    const ccomplex_t IT_0087 = (-6)*IT_0000;
    const ccomplex_t IT_0088 = IT_0049*IT_0087;
    const ccomplex_t IT_0089 = 24*IT_0026;
    const ccomplex_t IT_0090 = (-24)*IT_0026;
    const ccomplex_t IT_0091 = conj(IT_0031)*IT_0077;
    const ccomplex_t IT_0092 = (-6)*IT_0078;
    const ccomplex_t IT_0093 = conj(IT_0070)*IT_0092;
    const ccomplex_t IT_0094 = 0.166666666666667*IT_0079;
    const ccomplex_t IT_0095 = 0.166666666666667*conj(IT_0056)*IT_0071 +
       0.166666666666667*conj(IT_0031)*IT_0074 + conj(IT_0019)*IT_0082 +
       0.166666666666667*conj(IT_0038)*IT_0084 + 0.166666666666667*conj(IT_0048)
      *IT_0085 + 0.166666666666667*s_13*conj(IT_0070)*IT_0087 +
       0.166666666666667*conj(IT_0072)*IT_0088 + 0.166666666666667*IT_0076
      *IT_0092 + IT_0066*IT_0094;
    const ccomplex_t IT_0096 = (IT_0025*((-24)*conj(IT_0019) + 24*conj(IT_0025
      )) + IT_0019*(24*conj(IT_0019) + (-24)*conj(IT_0025)))*IT_0026 + conj
      (IT_0031)*(IT_0025*IT_0036 + IT_0038*IT_0040 + IT_0048*IT_0050 + IT_0032
      *IT_0051) + conj(IT_0056)*(IT_0040*IT_0048 + IT_0038*IT_0050 + IT_0036
      *IT_0060 + IT_0061*IT_0063 + IT_0064*IT_0065) + IT_0056*(IT_0040*conj
      (IT_0048) + conj(IT_0038)*IT_0050 + IT_0036*conj(IT_0060) + conj(IT_0064)
      *IT_0065 + IT_0061*IT_0068 + conj(IT_0056)*IT_0069 + conj(IT_0070)*IT_0071
       + conj(IT_0072)*IT_0074 + conj(IT_0031)*IT_0075) + IT_0031*(conj(IT_0025)
      *IT_0036 + conj(IT_0038)*IT_0040 + conj(IT_0048)*IT_0050 + IT_0032*IT_0067
       + conj(IT_0031)*IT_0069 + IT_0071*conj(IT_0072) + conj(IT_0070)*IT_0074 +
       conj(IT_0056)*IT_0075 + IT_0076*IT_0077) + conj(IT_0072)*(IT_0063*IT_0078
       + IT_0064*IT_0083 + IT_0048*IT_0084 + IT_0038*IT_0085 + IT_0060*IT_0086) 
      + IT_0072*(conj(IT_0031)*IT_0071 + conj(IT_0056)*IT_0074 + IT_0068*IT_0078
       + conj(IT_0064)*IT_0083 + conj(IT_0048)*IT_0084 + conj(IT_0038)*IT_0085 +
       conj(IT_0060)*IT_0086 + s_13*conj(IT_0072)*IT_0087 + conj(IT_0070)
      *IT_0088) + IT_0064*(conj(IT_0064)*IT_0089 + conj(IT_0060)*IT_0090 +
       IT_0091 + IT_0093) + IT_0060*(conj(IT_0060)*IT_0089 + conj(IT_0064)
      *IT_0090 + IT_0091 + IT_0093) + 6*conj(IT_0070)*(IT_0019*IT_0082 +
       0.166666666666667*IT_0038*IT_0084 + 0.166666666666667*IT_0048*IT_0085 +
       IT_0062*IT_0094) + 6*IT_0070*IT_0095;
    return create_ccomplex_return(IT_0096);
}

