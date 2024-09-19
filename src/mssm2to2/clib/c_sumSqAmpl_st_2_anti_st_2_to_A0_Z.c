#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
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
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = pow(m_W, -1);
    const ccomplex_t IT_0010 = m_t*mu_h*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0003, -1);
    const ccomplex_t IT_0012 = IT_0002*IT_0004;
    const ccomplex_t IT_0013 = -IT_0003*IT_0005 + IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0000, -2);
    const ccomplex_t IT_0015 = pow(m_t, 2);
    const ccomplex_t IT_0016 = (-2)*U_st_01*((0.5*m_W*IT_0013*(IT_0007 + -1./3
      *IT_0006*IT_0014) + IT_0005*IT_0007*IT_0009*IT_0011*IT_0015)*conj(U_st_01)
       + -1./2*IT_0002*IT_0007*IT_0010*IT_0011*conj(U_st_11)) + U_st_11*(IT_0002
      *IT_0007*IT_0010*IT_0011*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0006
      *IT_0013*IT_0014 + 1.5*IT_0005*IT_0007*IT_0009*IT_0011*IT_0015)*conj
      (U_st_11));
    const ccomplex_t IT_0017 = (0 + _Complex_I*-1)*(e_em*IT_0016 + (
      -1.4142135623731)*A_t*IT_0005*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = pow(m_st_2, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_12 + (-2)*IT_0019 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0008*IT_0022;
    const ccomplex_t IT_0024 = e_em*IT_0001*(IT_0002*IT_0004 + IT_0003*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0004*IT_0005;
    const ccomplex_t IT_0027 = IT_0002*IT_0003;
    const ccomplex_t IT_0028 = IT_0026 + IT_0027;
    const ccomplex_t IT_0029 = 2*U_st_01*((IT_0002*IT_0007*IT_0009*IT_0011
      *IT_0015 + (-0.5)*m_W*(IT_0007 + -1./3*IT_0006*IT_0014)*IT_0028)*conj
      (U_st_01) + 1./2*IT_0005*IT_0007*IT_0010*IT_0011*conj(U_st_11)) + U_st_11*
      (IT_0005*IT_0007*IT_0010*IT_0011*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0002*IT_0007*IT_0009*IT_0011*IT_0015 + m_W*IT_0006*IT_0014
      *IT_0028)*conj(U_st_11));
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(e_em*IT_0029 +
       1.4142135623731*A_t*IT_0002*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0019 + -reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0031*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0025*IT_0034;
    const ccomplex_t IT_0036 = -IT_0023 + -IT_0035;
    const ccomplex_t IT_0037 = pow(m_Z, -2);
    const ccomplex_t IT_0038 = pow(s_14, 2);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0006 
      + (-3)*IT_0000*IT_0007)*conj(U_st_00)*U_st_01 + IT_0001*IT_0006*conj
      (U_st_10)*U_st_11);
    const ccomplex_t IT_0041 = 1.33333333333333*IT_0040;
    const ccomplex_t IT_0042 = m_t*e_em*mu_h*IT_0009*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = m_t*e_em*mu_h*IT_0009*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0047 = IT_0004*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0050 = IT_0004*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0043 + 0.5*IT_0045 + (-0.5)*IT_0048 +
       0.5*IT_0051;
    const ccomplex_t IT_0053 = pow(m_A0, 2);
    const ccomplex_t IT_0054 = cpow((-2)*s_23 + IT_0019 + IT_0053 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = IT_0041*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0006 
      + (-3)*IT_0000*IT_0007)*U_st_01*conj(U_st_01) + IT_0001*IT_0006*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0059 = 1.33333333333333*IT_0058;
    const ccomplex_t IT_0060 = m_t*e_em*mu_h*IT_0009*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0061 = IT_0007*IT_0060;
    const ccomplex_t IT_0062 = m_t*e_em*mu_h*IT_0009*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0063 = IT_0007*IT_0062;
    const ccomplex_t IT_0064 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0065 = IT_0004*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0068 = IT_0004*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0061 + 0.5*IT_0063 + (-0.5)*IT_0066 +
       0.5*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0053 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0059*IT_0073;
    const ccomplex_t IT_0075 = -IT_0057 + -IT_0074;
    const ccomplex_t IT_0076 = s_14*s_34*IT_0037;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = s_13 + IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0006 
      + (-3)*IT_0000*IT_0007)*U_st_00*conj(U_st_01) + IT_0001*IT_0006*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0080 = 1.33333333333333*IT_0079;
    const ccomplex_t IT_0081 = m_t*e_em*mu_h*IT_0009*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0082 = IT_0007*IT_0081;
    const ccomplex_t IT_0083 = m_t*e_em*mu_h*IT_0009*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0084 = IT_0007*IT_0083;
    const ccomplex_t IT_0085 = A_t*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0086 = IT_0004*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = A_t*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0089 = IT_0004*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0082 + (-0.5)*IT_0084 + 0.5*IT_0087 + (
      -0.5)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_13 + IT_0019 + IT_0053 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0080*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0053 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0070*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0059*IT_0098;
    const ccomplex_t IT_0100 = -IT_0023 + -IT_0035 + IT_0057 + IT_0074 +
       IT_0095 + IT_0099;
    const ccomplex_t IT_0101 = s_14*s_24*IT_0037;
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = s_12 + IT_0102;
    const ccomplex_t IT_0104 = (-3)*conj(IT_0036);
    const ccomplex_t IT_0105 = pow(s_34, 2);
    const ccomplex_t IT_0106 = IT_0037*IT_0105;
    const ccomplex_t IT_0107 = s_24*s_34*IT_0037;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = s_23 + IT_0108;
    const ccomplex_t IT_0110 = pow(s_24, 2);
    const ccomplex_t IT_0111 = IT_0037*IT_0110;
    const ccomplex_t IT_0112 = (-3)*IT_0100;
    const ccomplex_t IT_0113 = (-3)*IT_0036*(conj(IT_0036)*(IT_0019 + -IT_0039
      ) + conj(IT_0075)*IT_0078 + conj(IT_0100)*IT_0103) + IT_0100*IT_0103
      *IT_0104 + (-3)*IT_0075*((-0.333333333333333)*IT_0078*IT_0104 + conj
      (IT_0075)*(IT_0053 + -IT_0106) + conj(IT_0100)*IT_0109) + (conj(IT_0075)
      *IT_0109 + conj(IT_0100)*(IT_0019 + -IT_0111))*IT_0112;
    return create_ccomplex_return(IT_0113);
}

