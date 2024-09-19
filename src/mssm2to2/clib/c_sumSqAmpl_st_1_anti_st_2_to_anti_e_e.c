#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = m_e*e_em*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = m_t*e_em*mu_h*IT_0000*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0002*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = IT_0005*IT_0009;
    const ccomplex_t IT_0011 = m_t*e_em*mu_h*IT_0000*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = 1.4142135623731*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0010 + (-0.5)*IT_0012 + (-0.5)*IT_0015 +
       0.5*IT_0018;
    const ccomplex_t IT_0020 = IT_0004*IT_0019;
    const ccomplex_t IT_0021 = pow(m_Z, 2);
    const ccomplex_t IT_0022 = pow(m_st_1, 2);
    const ccomplex_t IT_0023 = pow(m_st_2, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + IT_0021 + -IT_0022 + -IT_0023 
      + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0020*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = cpow(IT_0006, -1);
    const ccomplex_t IT_0028 = m_e*e_em*IT_0000*IT_0002*IT_0007*IT_0027;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0002*IT_0005;
    const ccomplex_t IT_0031 = IT_0002*IT_0011;
    const ccomplex_t IT_0032 = IT_0006*IT_0013;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = IT_0006*IT_0016;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (-0.5)*IT_0030 + 0.5*IT_0031 + (-0.5)*IT_0033 +
       0.5*IT_0035;
    const ccomplex_t IT_0037 = IT_0029*IT_0036;
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = -IT_0026 + -IT_0040;
    const ccomplex_t IT_0042 = cos(theta_W);
    const ccomplex_t IT_0043 = cpow(IT_0042, -1);
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0002
      *IT_0042 + -1./3*IT_0001*IT_0043)*U_st_00*conj(U_st_01) + IT_0001*IT_0043
      *U_st_10*conj(U_st_11));
    const ccomplex_t IT_0045 = 1.33333333333333*IT_0044;
    const ccomplex_t IT_0046 = IT_0001*IT_0043;
    const ccomplex_t IT_0047 = e_em*IT_0046;
    const ccomplex_t IT_0048 = IT_0002*IT_0042;
    const ccomplex_t IT_0049 = e_em*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0047 + -IT_0049);
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0022 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0045*IT_0054;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0043;
    const ccomplex_t IT_0058 = IT_0052*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0045*IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = s_13 + s_14;
    const ccomplex_t IT_0063 = m_e*IT_0062;
    const ccomplex_t IT_0064 = pow(m_e, 2);
    const ccomplex_t IT_0065 = s_34 + IT_0064;
    const ccomplex_t IT_0066 = m_e*IT_0045;
    const ccomplex_t IT_0067 = 0.666666666666667*IT_0044;
    const ccomplex_t IT_0068 = m_e*IT_0067;
    const ccomplex_t IT_0069 = -(IT_0054 + -IT_0059)*(IT_0066 + -IT_0068);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = (IT_0054 + -IT_0059)*(IT_0066 + -IT_0068);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = (-6)*conj(IT_0074);
    const ccomplex_t IT_0076 = 12*s_34;
    const ccomplex_t IT_0077 = cos(alpha);
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0027*IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0078;
    const ccomplex_t IT_0080 = sin(alpha);
    const ccomplex_t IT_0081 = m_t*mu_h*IT_0000;
    const ccomplex_t IT_0082 = IT_0006*IT_0077;
    const ccomplex_t IT_0083 = -IT_0007*IT_0080 + IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0042, -2);
    const ccomplex_t IT_0085 = pow(m_t, 2);
    const ccomplex_t IT_0086 = (-2)*U_st_00*((0.5*m_W*IT_0083*(IT_0002 + -1./3
      *IT_0001*IT_0084) + IT_0000*IT_0002*IT_0008*IT_0080*IT_0085)*conj(U_st_01)
       + -1./2*IT_0002*IT_0008*IT_0077*IT_0081*conj(U_st_11)) + U_st_10*(IT_0002
      *IT_0008*IT_0077*IT_0081*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0001
      *IT_0083*IT_0084 + 1.5*IT_0000*IT_0002*IT_0008*IT_0080*IT_0085)*conj
      (U_st_11));
    const ccomplex_t IT_0087 = (0 + _Complex_I*-1)*(e_em*IT_0086 + (
      -1.4142135623731)*A_t*IT_0080*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0079*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_12 + -IT_0022 + -IT_0023 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0027*IT_0080;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0006*IT_0080;
    const ccomplex_t IT_0096 = IT_0007*IT_0077;
    const ccomplex_t IT_0097 = IT_0095 + IT_0096;
    const ccomplex_t IT_0098 = 2*U_st_00*((IT_0000*IT_0002*IT_0008*IT_0077
      *IT_0085 + (-0.5)*m_W*(IT_0002 + -1./3*IT_0001*IT_0084)*IT_0097)*conj
      (U_st_01) + 1./2*IT_0002*IT_0008*IT_0080*IT_0081*conj(U_st_11)) + U_st_10*
      (IT_0002*IT_0008*IT_0080*IT_0081*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0000*IT_0002*IT_0008*IT_0077*IT_0085 + m_W*IT_0001*IT_0084
      *IT_0097)*conj(U_st_11));
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(e_em*IT_0098 +
       1.4142135623731*A_t*IT_0077*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0094*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0022 + -IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0092 + -IT_0104;
    const ccomplex_t IT_0106 = -s_14;
    const ccomplex_t IT_0107 = s_13 + IT_0106;
    const ccomplex_t IT_0108 = m_e*IT_0107;
    const ccomplex_t IT_0109 = -IT_0064;
    const ccomplex_t IT_0110 = s_34 + IT_0109;
    const ccomplex_t IT_0111 = m_e*s_13;
    const ccomplex_t IT_0112 = IT_0022*IT_0064;
    const ccomplex_t IT_0113 = m_e*s_14;
    const ccomplex_t IT_0114 = s_34*IT_0022;
    const ccomplex_t IT_0115 = (-6)*IT_0114;
    const ccomplex_t IT_0116 = s_13*s_14;
    const ccomplex_t IT_0117 = 12*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = conj(IT_0105)*IT_0108;
    const ccomplex_t IT_0120 = (-6)*conj(IT_0056);
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0113;
    const ccomplex_t IT_0122 = 6*s_34;
    const ccomplex_t IT_0123 = conj(IT_0105)*IT_0110;
    const ccomplex_t IT_0124 = (-6)*IT_0074;
    const ccomplex_t IT_0125 = conj(IT_0041)*IT_0065;
    const ccomplex_t IT_0126 = (-6)*IT_0041*((conj(IT_0056) + -conj(IT_0061))
      *IT_0063 + -IT_0065*(conj(IT_0071) + 0.166666666666667*IT_0075) + (
      -0.166666666666667)*conj(IT_0041)*(12*IT_0064 + IT_0076)) + 6*IT_0105*
      (0.166666666666667*((-12)*IT_0064 + IT_0076)*conj(IT_0105) + (conj(IT_0056
      ) + conj(IT_0061))*IT_0108 + (conj(IT_0071) + conj(IT_0074))*IT_0110) + (
      -6)*IT_0056*(conj(IT_0041)*IT_0063 + -conj(IT_0074)*IT_0111 + conj(IT_0061
      )*IT_0112 + conj(IT_0071)*IT_0113 + (-0.166666666666667)*conj(IT_0056)
      *IT_0118 + -IT_0119) + 6*IT_0061*(conj(IT_0041)*IT_0063 + conj(IT_0071)
      *IT_0111 + 0.166666666666667*conj(IT_0061)*IT_0118 + IT_0119 +
       0.166666666666667*IT_0112*IT_0120 + IT_0075*IT_0121) + 6*IT_0074*(conj
      (IT_0056)*IT_0111 + 0.166666666666667*conj(IT_0074)*IT_0122 + IT_0123) +
       IT_0124*(IT_0064*conj(IT_0071) + conj(IT_0061)*IT_0113 + IT_0125) + 6
      *IT_0071*(0.166666666666667*IT_0064*IT_0075 + conj(IT_0061)*IT_0111 +
       IT_0120*IT_0121 + 0.166666666666667*conj(IT_0071)*IT_0122 + IT_0123 +
       IT_0125);
    return create_ccomplex_return(IT_0126);
}

