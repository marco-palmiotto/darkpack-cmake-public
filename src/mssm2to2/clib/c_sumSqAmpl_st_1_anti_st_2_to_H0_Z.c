#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_00*conj(U_st_01) + IT_0001*IT_0002*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0006 = cos(beta);
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = cpow(IT_0000, -4);
    const ccomplex_t IT_0009 = IT_0002*IT_0006*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = IT_0002*IT_0008*IT_0010*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0002, -3);
    const ccomplex_t IT_0014 = IT_0006*IT_0007*IT_0013;
    const ccomplex_t IT_0015 = IT_0010*IT_0011*IT_0013;
    const ccomplex_t IT_0016 = cpow(IT_0000, -2);
    const ccomplex_t IT_0017 = IT_0003*IT_0006*IT_0007*IT_0016;
    const ccomplex_t IT_0018 = IT_0003*IT_0010*IT_0011*IT_0016;
    const ccomplex_t IT_0019 = cpow(IT_0002, -2);
    const ccomplex_t IT_0020 = cpow(IT_0016 + IT_0019, -1);
    const ccomplex_t IT_0021 = (IT_0009 + IT_0012 + IT_0014 + IT_0015 + 2
      *IT_0017 + 2*IT_0018)*IT_0020;
    const ccomplex_t IT_0022 = m_W*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = pow(m_st_1, 2);
    const ccomplex_t IT_0026 = pow(m_st_2, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0025 + -IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0024*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0029;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = pow(m_Z, -2);
    const ccomplex_t IT_0033 = pow(s_14, 2);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_00*conj(U_st_00) + IT_0001*IT_0002*U_st_10
      *conj(U_st_10));
    const ccomplex_t IT_0036 = 1.33333333333333*IT_0035;
    const ccomplex_t IT_0037 = pow(m_W, -1);
    const ccomplex_t IT_0038 = m_t*mu_h*IT_0037;
    const ccomplex_t IT_0039 = cpow(IT_0010, -1);
    const ccomplex_t IT_0040 = IT_0006*IT_0007;
    const ccomplex_t IT_0041 = -IT_0010*IT_0011 + IT_0040;
    const ccomplex_t IT_0042 = pow(m_t, 2);
    const ccomplex_t IT_0043 = (-2)*U_st_00*((0.5*m_W*(IT_0003 + -1./3*IT_0002
      *IT_0016)*IT_0041 + IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj(U_st_01)
       + -1./2*IT_0003*IT_0007*IT_0038*IT_0039*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0007*IT_0038*IT_0039*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0016*IT_0041 + 1.5*IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj
      (U_st_11));
    const ccomplex_t IT_0044 = (0 + _Complex_I*-1)*(e_em*IT_0043 + (
      -1.4142135623731)*A_t*IT_0011*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = pow(m_H0, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0026 + IT_0046 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0036*IT_0049;
    const ccomplex_t IT_0051 = e_em*IT_0001*IT_0003*(IT_0007*IT_0010 + 
      -IT_0006*IT_0011);
    const ccomplex_t IT_0052 = m_t*e_em*mu_h*IT_0037*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0053 = IT_0003*IT_0052;
    const ccomplex_t IT_0054 = m_t*e_em*mu_h*IT_0037*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0055 = IT_0003*IT_0054;
    const ccomplex_t IT_0056 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0057 = IT_0006*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0060 = IT_0006*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0053 + 0.5*IT_0055 + (-0.5)*IT_0058 +
       0.5*IT_0061;
    const ccomplex_t IT_0063 = cpow((-2)*s_12 + -IT_0025 + -IT_0026 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0051*IT_0065;
    const ccomplex_t IT_0067 = e_em*IT_0001*IT_0003*(IT_0006*IT_0007 + IT_0010
      *IT_0011);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0003*IT_0006*IT_0039;
    const ccomplex_t IT_0070 = IT_0052*IT_0069;
    const ccomplex_t IT_0071 = IT_0054*IT_0069;
    const ccomplex_t IT_0072 = IT_0010*IT_0056;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = IT_0010*IT_0059;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0070 + (-0.5)*IT_0071 + (-0.5)*IT_0073 +
       0.5*IT_0075;
    const ccomplex_t IT_0077 = pow(m_Z, 2);
    const ccomplex_t IT_0078 = cpow((-2)*s_12 + -IT_0025 + -IT_0026 + IT_0077 
      + -reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0068*IT_0080;
    const ccomplex_t IT_0082 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0083 = (-2)*U_st_01*((0.5*m_W*(IT_0003 + -1./3*IT_0002
      *IT_0016)*IT_0041 + IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj(U_st_01)
       + -1./2*IT_0003*IT_0007*IT_0038*IT_0039*conj(U_st_11)) + U_st_11*(IT_0003
      *IT_0007*IT_0038*IT_0039*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0016*IT_0041 + 1.5*IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj
      (U_st_11));
    const ccomplex_t IT_0084 = (0 + _Complex_I*-1)*(e_em*IT_0083 + (
      -1.4142135623731)*A_t*IT_0011*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0046 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0082*IT_0088;
    const ccomplex_t IT_0090 = -IT_0050 + IT_0066 + IT_0081 + -IT_0089;
    const ccomplex_t IT_0091 = s_14*s_34*IT_0032;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = s_13 + IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_01*conj(U_st_01) + IT_0001*IT_0002*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0095 = 1.33333333333333*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0025 + IT_0046 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0045*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0095*IT_0098;
    const ccomplex_t IT_0100 = (-2)*U_st_00*((0.5*m_W*(IT_0003 + -1./3*IT_0002
      *IT_0016)*IT_0041 + IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj(U_st_00)
       + -1./2*IT_0003*IT_0007*IT_0038*IT_0039*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0007*IT_0038*IT_0039*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0016*IT_0041 + 1.5*IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj
      (U_st_10));
    const ccomplex_t IT_0101 = (0 + _Complex_I*-1)*(e_em*IT_0100 + (
      -1.4142135623731)*A_t*IT_0011*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0102 = (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_13 + IT_0046 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0082*IT_0105;
    const ccomplex_t IT_0107 = (-2)*U_st_01*((0.5*m_W*(IT_0003 + -1./3*IT_0002
      *IT_0016)*IT_0041 + IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj(U_st_01)
       + -1./2*IT_0003*IT_0007*IT_0038*IT_0039*conj(U_st_11)) + U_st_11*(IT_0003
      *IT_0007*IT_0038*IT_0039*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0016*IT_0041 + 1.5*IT_0003*IT_0011*IT_0037*IT_0039*IT_0042)*conj
      (U_st_11));
    const ccomplex_t IT_0108 = (0 + _Complex_I*-1)*(e_em*IT_0107 + (
      -1.4142135623731)*A_t*IT_0011*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = IT_0086*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0082*IT_0111;
    const ccomplex_t IT_0113 = -IT_0082;
    const ccomplex_t IT_0114 = IT_0005 + IT_0113;
    const ccomplex_t IT_0115 = IT_0029*IT_0114;
    const ccomplex_t IT_0116 = IT_0050 + IT_0099 + IT_0106 + IT_0112 + -IT_0115;
    const ccomplex_t IT_0117 = s_14*s_24*IT_0032;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = s_12 + IT_0118;
    const ccomplex_t IT_0120 = (-3)*conj(IT_0031);
    const ccomplex_t IT_0121 = pow(s_34, 2);
    const ccomplex_t IT_0122 = IT_0032*IT_0121;
    const ccomplex_t IT_0123 = s_24*s_34*IT_0032;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = s_23 + IT_0124;
    const ccomplex_t IT_0126 = pow(s_24, 2);
    const ccomplex_t IT_0127 = IT_0032*IT_0126;
    const ccomplex_t IT_0128 = (-3)*IT_0116;
    const ccomplex_t IT_0129 = (-3)*IT_0031*(conj(IT_0031)*(IT_0025 + -IT_0034
      ) + conj(IT_0090)*IT_0093 + conj(IT_0116)*IT_0119) + IT_0116*IT_0119
      *IT_0120 + (-3)*IT_0090*((-0.333333333333333)*IT_0093*IT_0120 + conj
      (IT_0090)*(IT_0046 + -IT_0122) + conj(IT_0116)*IT_0125) + (conj(IT_0090)
      *IT_0125 + conj(IT_0116)*(IT_0026 + -IT_0127))*IT_0128;
    return create_ccomplex_return(IT_0129);
}

