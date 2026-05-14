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
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_00*conj(U_st_01) + IT_0001*IT_0002*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0002;
    const ccomplex_t IT_0007 = e_em*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0003;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0007 + (-3)*IT_0009);
    const ccomplex_t IT_0011 = (-0.166666666666667)*IT_0010;
    const ccomplex_t IT_0012 = pow(m_st_1, 2);
    const ccomplex_t IT_0013 = pow(m_st_2, 2);
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0020 = m_c*IT_0019;
    const ccomplex_t IT_0021 = m_c*IT_0005;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0023 = (-0.666666666666667)*IT_0022;
    const ccomplex_t IT_0024 = IT_0014*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = (IT_0020 + -IT_0021)*(IT_0016 + -IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = m_c*s_13;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = m_t*e_em*mu_h*IT_0031*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = m_t*e_em*mu_h*IT_0031*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0035 = IT_0003*IT_0034;
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0041 = IT_0036*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0033 + 0.5*IT_0035 + (-0.5)*IT_0039 +
       0.5*IT_0042;
    const ccomplex_t IT_0044 = sin(beta);
    const ccomplex_t IT_0045 = cpow(IT_0044, -1);
    const ccomplex_t IT_0046 = m_c*e_em*IT_0003*IT_0031*IT_0036*IT_0045;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = IT_0043*IT_0047;
    const ccomplex_t IT_0049 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = m_c*e_em*IT_0003*IT_0031;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0003*IT_0036*IT_0045;
    const ccomplex_t IT_0055 = IT_0032*IT_0054;
    const ccomplex_t IT_0056 = IT_0034*IT_0054;
    const ccomplex_t IT_0057 = IT_0037*IT_0044;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = IT_0040*IT_0044;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0055 + (-0.5)*IT_0056 + (-0.5)*IT_0058 +
       0.5*IT_0060;
    const ccomplex_t IT_0062 = IT_0053*IT_0061;
    const ccomplex_t IT_0063 = pow(m_Z, 2);
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + IT_0063 
      + -reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = -IT_0051 + -IT_0066;
    const ccomplex_t IT_0068 = s_13 + s_14;
    const ccomplex_t IT_0069 = m_c*IT_0068;
    const ccomplex_t IT_0070 = sin(alpha);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_c*e_em*IT_0003*IT_0031
      *IT_0045*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = m_t*mu_h*IT_0031;
    const ccomplex_t IT_0074 = cos(alpha);
    const ccomplex_t IT_0075 = IT_0036*IT_0074;
    const ccomplex_t IT_0076 = -IT_0044*IT_0070 + IT_0075;
    const ccomplex_t IT_0077 = cpow(IT_0000, -2);
    const ccomplex_t IT_0078 = pow(m_t, 2);
    const ccomplex_t IT_0079 = (-2)*U_st_00*((0.5*m_W*IT_0076*(IT_0003 + -1./3
      *IT_0002*IT_0077) + IT_0003*IT_0031*IT_0045*IT_0070*IT_0078)*conj(U_st_01)
       + -1./2*IT_0003*IT_0045*IT_0073*IT_0074*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0045*IT_0073*IT_0074*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0076*IT_0077 + 1.5*IT_0003*IT_0031*IT_0045*IT_0070*IT_0078)*conj
      (U_st_11));
    const ccomplex_t IT_0080 = (0 + _Complex_I*-1)*(e_em*IT_0079 + (
      -1.4142135623731)*A_t*IT_0070*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = IT_0072*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + -IT_0012 + -IT_0013 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*m_c*e_em*IT_0003*IT_0031
      *IT_0045*IT_0074;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0036*IT_0070;
    const ccomplex_t IT_0089 = IT_0044*IT_0074;
    const ccomplex_t IT_0090 = IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = 2*U_st_00*((IT_0003*IT_0031*IT_0045*IT_0074
      *IT_0078 + (-0.5)*m_W*(IT_0003 + -1./3*IT_0002*IT_0077)*IT_0090)*conj
      (U_st_01) + 1./2*IT_0003*IT_0045*IT_0070*IT_0073*conj(U_st_11)) + U_st_10*
      (IT_0003*IT_0045*IT_0070*IT_0073*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0003*IT_0031*IT_0045*IT_0074*IT_0078 + m_W*IT_0002*IT_0077
      *IT_0090)*conj(U_st_11));
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(e_em*IT_0091 +
       1.4142135623731*A_t*IT_0074*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = IT_0087*IT_0093;
    const ccomplex_t IT_0095 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0012 + -IT_0013 + -reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = -IT_0085 + -IT_0097;
    const ccomplex_t IT_0099 = -s_14;
    const ccomplex_t IT_0100 = s_13 + IT_0099;
    const ccomplex_t IT_0101 = m_c*IT_0100;
    const ccomplex_t IT_0102 = s_34*IT_0012;
    const ccomplex_t IT_0103 = (-18)*IT_0102;
    const ccomplex_t IT_0104 = s_13*s_14;
    const ccomplex_t IT_0105 = 36*IT_0104;
    const ccomplex_t IT_0106 = IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = IT_0005*IT_0025;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = pow(m_c, 2);
    const ccomplex_t IT_0110 = IT_0012*IT_0109;
    const ccomplex_t IT_0111 = (-18)*IT_0018;
    const ccomplex_t IT_0112 = m_c*s_14;
    const ccomplex_t IT_0113 = (-18)*IT_0029;
    const ccomplex_t IT_0114 = (-18)*conj(IT_0018);
    const ccomplex_t IT_0115 = (-18)*conj(IT_0029);
    const ccomplex_t IT_0116 = 18*s_34;
    const ccomplex_t IT_0117 = 18*IT_0109 + IT_0116;
    const ccomplex_t IT_0118 = (-18)*IT_0109;
    const ccomplex_t IT_0119 = (-18)*s_34 + IT_0118;
    const ccomplex_t IT_0120 = 18*IT_0108;
    const ccomplex_t IT_0121 = IT_0116 + IT_0118;
    const ccomplex_t IT_0122 = 18*IT_0018;
    const ccomplex_t IT_0123 = IT_0098*IT_0121;
    const ccomplex_t IT_0124 = conj(IT_0018)*(18*IT_0029*IT_0030 + (-18)
      *IT_0067*IT_0069 + 18*IT_0098*IT_0101 + IT_0018*IT_0106) + conj(IT_0108)*
      (18*IT_0026*IT_0030 + 18*IT_0067*IT_0069 + 18*IT_0098*IT_0101 + IT_0106
      *IT_0108 + IT_0110*IT_0111 + IT_0112*IT_0113) + IT_0026*(IT_0112*IT_0114 +
       IT_0109*IT_0115) + IT_0108*(IT_0110*IT_0114 + IT_0112*IT_0115) + conj
      (IT_0067)*(IT_0067*(36*s_34 + 36*IT_0109) + IT_0026*IT_0117 + IT_0029
      *IT_0119 + IT_0069*((-18)*IT_0018 + IT_0120)) + conj(IT_0098)*(IT_0098*(36
      *s_34 + (-36)*IT_0109) + (IT_0026 + IT_0029)*IT_0121 + IT_0101*(IT_0120 +
       IT_0122)) + conj(IT_0026)*(IT_0111*IT_0112 + IT_0109*IT_0113 + IT_0026
      *IT_0116 + IT_0067*IT_0117 + IT_0030*IT_0120 + IT_0123) + conj(IT_0029)*
      (IT_0029*IT_0116 + IT_0067*IT_0119 + IT_0030*IT_0122 + IT_0123);
    return create_ccomplex_return(IT_0124);
}

