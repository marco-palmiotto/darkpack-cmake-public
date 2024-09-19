#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_Z_Z(
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
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*IT_0001;
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0004*IT_0005 
      + (-3)*IT_0003*IT_0006)*U_st_01*conj(U_st_01) + IT_0004*IT_0005*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0008 = 1.33333333333333*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = pow(m_st_2, 2);
    const ccomplex_t IT_0012 = cpow((-2)*s_13 + IT_0000 + IT_0011 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0004*IT_0005 
      + (-3)*IT_0003*IT_0006)*conj(U_st_00)*U_st_01 + IT_0004*IT_0005*conj
      (U_st_10)*U_st_11);
    const ccomplex_t IT_0015 = 1.33333333333333*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0004*IT_0005 
      + (-3)*IT_0003*IT_0006)*U_st_00*conj(U_st_01) + IT_0004*IT_0005*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0017 = 1.33333333333333*IT_0016;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = IT_0013*IT_0018;
    const ccomplex_t IT_0020 = -IT_0010 + -IT_0019;
    const ccomplex_t IT_0021 = cpow((-2)*s_23 + IT_0000 + IT_0011 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0018*IT_0022;
    const ccomplex_t IT_0024 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = IT_0023 + IT_0026;
    const ccomplex_t IT_0028 = pow(m_Z, -2);
    const ccomplex_t IT_0029 = s_13*s_23*IT_0028;
    const ccomplex_t IT_0030 = pow(s_24, 2);
    const ccomplex_t IT_0031 = IT_0028*IT_0030;
    const ccomplex_t IT_0032 = -(s_12 + -IT_0029)*(IT_0011 + -IT_0031);
    const ccomplex_t IT_0033 = pow(e_em, 2);
    const ccomplex_t IT_0034 = cpow(IT_0003, -2);
    const ccomplex_t IT_0035 = U_st_01*conj(U_st_01);
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = cpow(IT_0003, -4);
    const ccomplex_t IT_0038 = cpow(IT_0005, 2);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = IT_0035*IT_0039;
    const ccomplex_t IT_0041 = U_st_11*conj(U_st_11);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = cpow(IT_0005, -2);
    const ccomplex_t IT_0044 = IT_0035*IT_0043;
    const ccomplex_t IT_0045 = tan(theta_W);
    const ccomplex_t IT_0046 = cpow(IT_0045, 2);
    const ccomplex_t IT_0047 = cpow(1 + IT_0046, -1);
    const ccomplex_t IT_0048 = (IT_0036 + (-0.166666666666667)*IT_0040 + (
      -2.66666666666667)*IT_0042 + (-1.5)*IT_0044)*IT_0047;
    const ccomplex_t IT_0049 = IT_0033*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = (-0.333333333333333)*IT_0050;
    const ccomplex_t IT_0052 = m_W*e_em;
    const ccomplex_t IT_0053 = cos(alpha);
    const ccomplex_t IT_0054 = sin(beta);
    const ccomplex_t IT_0055 = IT_0005*IT_0037*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = cos(beta);
    const ccomplex_t IT_0057 = sin(alpha);
    const ccomplex_t IT_0058 = IT_0005*IT_0037*IT_0056*IT_0057;
    const ccomplex_t IT_0059 = cpow(IT_0005, -3);
    const ccomplex_t IT_0060 = IT_0053*IT_0054*IT_0059;
    const ccomplex_t IT_0061 = IT_0056*IT_0057*IT_0059;
    const ccomplex_t IT_0062 = IT_0006*IT_0034*IT_0053*IT_0054;
    const ccomplex_t IT_0063 = IT_0006*IT_0034*IT_0056*IT_0057;
    const ccomplex_t IT_0064 = cpow(IT_0034 + IT_0043, -1);
    const ccomplex_t IT_0065 = (IT_0055 + -IT_0058 + IT_0060 + -IT_0061 + 2
      *IT_0062 + (-2)*IT_0063)*IT_0064;
    const ccomplex_t IT_0066 = IT_0052*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = pow(m_W, -1);
    const ccomplex_t IT_0069 = m_t*mu_h*IT_0068;
    const ccomplex_t IT_0070 = cpow(IT_0054, -1);
    const ccomplex_t IT_0071 = IT_0056*IT_0057;
    const ccomplex_t IT_0072 = IT_0053*IT_0054;
    const ccomplex_t IT_0073 = IT_0071 + IT_0072;
    const ccomplex_t IT_0074 = pow(m_t, 2);
    const ccomplex_t IT_0075 = 2*U_st_01*(((-0.5)*m_W*(IT_0006 + -1./3*IT_0005
      *IT_0034)*IT_0073 + IT_0006*IT_0053*IT_0068*IT_0070*IT_0074)*conj(U_st_01)
       + 1./2*IT_0006*IT_0057*IT_0069*IT_0070*conj(U_st_11)) + U_st_11*(IT_0006
      *IT_0057*IT_0069*IT_0070*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0005
      *IT_0034*IT_0073 + (-1.5)*IT_0006*IT_0053*IT_0068*IT_0070*IT_0074)*conj
      (U_st_11));
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(e_em*IT_0075 +
       1.4142135623731*A_t*IT_0053*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0077 = (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = IT_0067*IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0011 + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = IT_0005*IT_0037*IT_0053*IT_0056;
    const ccomplex_t IT_0083 = IT_0005*IT_0037*IT_0054*IT_0057;
    const ccomplex_t IT_0084 = IT_0053*IT_0056*IT_0059;
    const ccomplex_t IT_0085 = IT_0054*IT_0057*IT_0059;
    const ccomplex_t IT_0086 = IT_0006*IT_0034*IT_0053*IT_0056;
    const ccomplex_t IT_0087 = IT_0006*IT_0034*IT_0054*IT_0057;
    const ccomplex_t IT_0088 = IT_0064*(IT_0082 + IT_0083 + IT_0084 + IT_0085 
      + 2*IT_0086 + 2*IT_0087);
    const ccomplex_t IT_0089 = IT_0052*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0053*IT_0056;
    const ccomplex_t IT_0092 = -IT_0054*IT_0057 + IT_0091;
    const ccomplex_t IT_0093 = (-2)*U_st_01*((IT_0006*IT_0057*IT_0068*IT_0070
      *IT_0074 + 0.5*m_W*(IT_0006 + -1./3*IT_0005*IT_0034)*IT_0092)*conj(U_st_01
      ) + -1./2*IT_0006*IT_0053*IT_0069*IT_0070*conj(U_st_11)) + U_st_11*
      (IT_0006*IT_0053*IT_0069*IT_0070*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0006*IT_0057*IT_0068*IT_0070*IT_0074 + m_W*IT_0005*IT_0034*IT_0092)
      *conj(U_st_11));
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*(e_em*IT_0093 + (
      -1.4142135623731)*A_t*IT_0057*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0090*IT_0095;
    const ccomplex_t IT_0097 = cpow((-2)*s_12 + (-2)*IT_0011 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0051 + -IT_0081 + -IT_0099;
    const ccomplex_t IT_0101 = pow(m_Z, -4);
    const ccomplex_t IT_0102 = s_13*s_24*s_34*IT_0101;
    const ccomplex_t IT_0103 = s_14*s_24*IT_0028;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = s_12 + -IT_0029 + IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = pow(s_34, 2);
    const ccomplex_t IT_0107 = IT_0101*IT_0106;
    const ccomplex_t IT_0108 = -IT_0023 + -IT_0026;
    const ccomplex_t IT_0109 = s_34*IT_0028;
    const ccomplex_t IT_0110 = s_23*s_34*IT_0028;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = s_24 + IT_0111;
    const ccomplex_t IT_0113 = IT_0109*IT_0112;
    const ccomplex_t IT_0114 = pow(s_23, 2);
    const ccomplex_t IT_0115 = IT_0028*IT_0114;
    const ccomplex_t IT_0116 = -IT_0011;
    const ccomplex_t IT_0117 = s_23*s_24*s_34*IT_0101;
    const ccomplex_t IT_0118 = -IT_0117;
    const ccomplex_t IT_0119 = IT_0031 + IT_0115 + IT_0116 + IT_0118;
    const ccomplex_t IT_0120 = pow(s_13, 2);
    const ccomplex_t IT_0121 = IT_0028*IT_0120;
    const ccomplex_t IT_0122 = -(IT_0011 + -IT_0031)*(IT_0011 + -IT_0121);
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = s_24*s_34*IT_0028;
    const ccomplex_t IT_0125 = (s_12 + -IT_0029)*(s_23 + -IT_0124);
    const ccomplex_t IT_0126 = 3*IT_0020;
    const ccomplex_t IT_0127 = IT_0028*IT_0106;
    const ccomplex_t IT_0128 = -(IT_0011 + -IT_0115)*(IT_0000 + -IT_0127);
    const ccomplex_t IT_0129 = (-3)*conj(IT_0100);
    const ccomplex_t IT_0130 = -(IT_0011 + -IT_0115)*(s_23 + -IT_0124);
    const ccomplex_t IT_0131 = -conj(IT_0027);
    const ccomplex_t IT_0132 = (IT_0031 + IT_0116)*(IT_0115 + IT_0116);
    const ccomplex_t IT_0133 = (-3)*IT_0020*(conj(IT_0027)*IT_0032 + -conj
      (IT_0100)*IT_0105) + 6*IT_0100*(0.5*conj(IT_0020)*IT_0105 + conj(IT_0100)*
      (1 + 0.5*IT_0107) + (-0.5)*conj(IT_0108)*IT_0113 + (-0.5)*conj(IT_0027)
      *IT_0119) + (conj(IT_0020)*IT_0123 + conj(IT_0108)*IT_0125)*IT_0126 + 3
      *IT_0108*(conj(IT_0020)*IT_0125 + -conj(IT_0108)*IT_0128 +
       0.333333333333333*IT_0113*IT_0129 + IT_0130*IT_0131) + (-3)*IT_0027*(conj
      (IT_0020)*IT_0032 + (-0.333333333333333)*IT_0119*IT_0129 + conj(IT_0108)
      *IT_0130 + IT_0131*IT_0132);
    return create_ccomplex_return(IT_0133);
}

