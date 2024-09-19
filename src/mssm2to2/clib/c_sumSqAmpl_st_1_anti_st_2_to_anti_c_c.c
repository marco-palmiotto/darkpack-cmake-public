#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_anti_c_c.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_anti_c_c(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_c = param->m_c;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0004 = (-0.666666666666667)*IT_0003;
    const ccomplex_t IT_0005 = pow(m_st_1, 2);
    const ccomplex_t IT_0006 = pow(m_st_2, 2);
    const ccomplex_t IT_0007 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0005 + -IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = IT_0001*IT_0002;
    const ccomplex_t IT_0011 = e_em*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0002, -1);
    const ccomplex_t IT_0013 = IT_0000*IT_0012;
    const ccomplex_t IT_0014 = e_em*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + (-3)*IT_0014);
    const ccomplex_t IT_0016 = (-0.166666666666667)*IT_0015;
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0012)*U_st_00*conj(U_st_01) + IT_0001*IT_0002*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0020 = 1.33333333333333*IT_0019;
    const ccomplex_t IT_0021 = m_c*IT_0020;
    const ccomplex_t IT_0022 = 0.666666666666667*IT_0019;
    const ccomplex_t IT_0023 = m_c*IT_0022;
    const ccomplex_t IT_0024 = (IT_0009 + -IT_0018)*(IT_0021 + -IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = pow(m_c, 2);
    const ccomplex_t IT_0028 = -(IT_0009 + -IT_0018)*(IT_0021 + -IT_0023);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (-18)*IT_0030;
    const ccomplex_t IT_0032 = m_c*s_14;
    const ccomplex_t IT_0033 = IT_0018*IT_0020;
    const ccomplex_t IT_0034 = -IT_0033;
    const ccomplex_t IT_0035 = (-18)*IT_0034;
    const ccomplex_t IT_0036 = IT_0009*IT_0020;
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = IT_0005*IT_0027;
    const ccomplex_t IT_0039 = m_c*s_13;
    const ccomplex_t IT_0040 = pow(m_Z, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_12 + -IT_0005 + -IT_0006 + IT_0040 
      + -reg_prop, -1);
    const ccomplex_t IT_0042 = pow(m_W, -1);
    const ccomplex_t IT_0043 = m_c*e_em*IT_0012*IT_0042;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = m_t*e_em*mu_h*IT_0042*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0046 = cos(beta);
    const ccomplex_t IT_0047 = sin(beta);
    const ccomplex_t IT_0048 = cpow(IT_0047, -1);
    const ccomplex_t IT_0049 = IT_0012*IT_0046*IT_0048;
    const ccomplex_t IT_0050 = IT_0045*IT_0049;
    const ccomplex_t IT_0051 = m_t*e_em*mu_h*IT_0042*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0054 = IT_0047*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0057 = IT_0047*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0050 + (-0.5)*IT_0052 + (-0.5)*IT_0055 +
       0.5*IT_0058;
    const ccomplex_t IT_0060 = IT_0044*IT_0059;
    const ccomplex_t IT_0061 = IT_0041*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0012*IT_0045;
    const ccomplex_t IT_0064 = IT_0012*IT_0051;
    const ccomplex_t IT_0065 = IT_0046*IT_0053;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = IT_0046*IT_0056;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0063 + 0.5*IT_0064 + (-0.5)*IT_0066 +
       0.5*IT_0068;
    const ccomplex_t IT_0070 = m_c*e_em*IT_0012*IT_0042*IT_0046*IT_0048;
    const ccomplex_t IT_0071 = (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_12 + -IT_0005 + -IT_0006 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = -IT_0062 + -IT_0075;
    const ccomplex_t IT_0077 = s_13 + s_14;
    const ccomplex_t IT_0078 = m_c*IT_0077;
    const ccomplex_t IT_0079 = s_34*IT_0005;
    const ccomplex_t IT_0080 = (-18)*IT_0079;
    const ccomplex_t IT_0081 = s_13*s_14;
    const ccomplex_t IT_0082 = 36*IT_0081;
    const ccomplex_t IT_0083 = IT_0080 + IT_0082;
    const ccomplex_t IT_0084 = sin(alpha);
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*m_c*e_em*IT_0012*IT_0042
      *IT_0048*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = m_t*mu_h*IT_0042;
    const ccomplex_t IT_0088 = cos(alpha);
    const ccomplex_t IT_0089 = IT_0046*IT_0088;
    const ccomplex_t IT_0090 = -IT_0047*IT_0084 + IT_0089;
    const ccomplex_t IT_0091 = cpow(IT_0000, -2);
    const ccomplex_t IT_0092 = pow(m_t, 2);
    const ccomplex_t IT_0093 = (-2)*U_st_00*((0.5*m_W*IT_0090*(IT_0012 + -1./3
      *IT_0002*IT_0091) + IT_0012*IT_0042*IT_0048*IT_0084*IT_0092)*conj(U_st_01)
       + -1./2*IT_0012*IT_0048*IT_0087*IT_0088*conj(U_st_11)) + U_st_10*(IT_0012
      *IT_0048*IT_0087*IT_0088*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0090*IT_0091 + 1.5*IT_0012*IT_0042*IT_0048*IT_0084*IT_0092)*conj
      (U_st_11));
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*(e_em*IT_0093 + (
      -1.4142135623731)*A_t*IT_0084*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0086*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_12 + -IT_0005 + -IT_0006 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_c*e_em*IT_0012*IT_0042
      *IT_0048*IT_0088;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0046*IT_0084;
    const ccomplex_t IT_0103 = IT_0047*IT_0088;
    const ccomplex_t IT_0104 = IT_0102 + IT_0103;
    const ccomplex_t IT_0105 = 2*U_st_00*((IT_0012*IT_0042*IT_0048*IT_0088
      *IT_0092 + (-0.5)*m_W*(IT_0012 + -1./3*IT_0002*IT_0091)*IT_0104)*conj
      (U_st_01) + 1./2*IT_0012*IT_0048*IT_0084*IT_0087*conj(U_st_11)) + U_st_10*
      (IT_0012*IT_0048*IT_0084*IT_0087*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0012*IT_0042*IT_0048*IT_0088*IT_0092 + m_W*IT_0002*IT_0091
      *IT_0104)*conj(U_st_11));
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*(e_em*IT_0105 +
       1.4142135623731*A_t*IT_0088*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = IT_0101*IT_0107;
    const ccomplex_t IT_0109 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0005 + -IT_0006 + -reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0108*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0099 + -IT_0111;
    const ccomplex_t IT_0113 = -s_14;
    const ccomplex_t IT_0114 = s_13 + IT_0113;
    const ccomplex_t IT_0115 = m_c*IT_0114;
    const ccomplex_t IT_0116 = 18*conj(IT_0112)*IT_0115;
    const ccomplex_t IT_0117 = (-18)*conj(IT_0030);
    const ccomplex_t IT_0118 = (-18)*conj(IT_0034);
    const ccomplex_t IT_0119 = 18*s_34;
    const ccomplex_t IT_0120 = 18*IT_0027 + IT_0119;
    const ccomplex_t IT_0121 = (-18)*IT_0027;
    const ccomplex_t IT_0122 = (-18)*s_34 + IT_0121;
    const ccomplex_t IT_0123 = 18*conj(IT_0037);
    const ccomplex_t IT_0124 = IT_0119 + IT_0121;
    const ccomplex_t IT_0125 = 18*conj(IT_0034);
    const ccomplex_t IT_0126 = conj(IT_0112)*IT_0124;
    const ccomplex_t IT_0127 = conj(IT_0026)*(IT_0027*IT_0031 + IT_0032
      *IT_0035) + conj(IT_0037)*(IT_0031*IT_0032 + IT_0035*IT_0038) + IT_0034*
      (18*conj(IT_0030)*IT_0039 + (-18)*conj(IT_0076)*IT_0078 + conj(IT_0034)
      *IT_0083 + IT_0116) + IT_0037*(18*conj(IT_0026)*IT_0039 + 18*conj(IT_0076)
      *IT_0078 + conj(IT_0037)*IT_0083 + IT_0116 + IT_0032*IT_0117 + IT_0038
      *IT_0118) + IT_0076*((36*s_34 + 36*IT_0027)*conj(IT_0076) + conj(IT_0026)
      *IT_0120 + conj(IT_0030)*IT_0122 + IT_0078*((-18)*conj(IT_0034) + IT_0123)
      ) + IT_0112*((36*s_34 + (-36)*IT_0027)*conj(IT_0112) + (conj(IT_0026) +
       conj(IT_0030))*IT_0124 + IT_0115*(IT_0123 + IT_0125)) + IT_0026*(IT_0027
      *IT_0117 + IT_0032*IT_0118 + conj(IT_0026)*IT_0119 + conj(IT_0076)*IT_0120
       + IT_0039*IT_0123 + IT_0126) + IT_0030*(conj(IT_0030)*IT_0119 + conj
      (IT_0076)*IT_0122 + IT_0039*IT_0125 + IT_0126);
    return create_ccomplex_return(IT_0127);
}

