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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_00*conj(U_st_01) + IT_0001*IT_0002*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = pow(m_st_1, 2);
    const ccomplex_t IT_0008 = pow(m_st_2, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0006*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = -IT_0012;
    const ccomplex_t IT_0014 = pow(m_e, 2);
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = m_e*s_14;
    const ccomplex_t IT_0017 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0018 = m_e*IT_0017;
    const ccomplex_t IT_0019 = m_e*IT_0005;
    const ccomplex_t IT_0020 = IT_0001*IT_0002;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0000*IT_0003;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023);
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = (IT_0018 + -IT_0019)*(IT_0011 + -IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0005*IT_0027;
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = (-6)*conj(IT_0031);
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = m_e*e_em*IT_0003*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = m_t*e_em*mu_h*IT_0033*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0037 = cos(beta);
    const ccomplex_t IT_0038 = sin(beta);
    const ccomplex_t IT_0039 = cpow(IT_0038, -1);
    const ccomplex_t IT_0040 = IT_0003*IT_0037*IT_0039;
    const ccomplex_t IT_0041 = IT_0036*IT_0040;
    const ccomplex_t IT_0042 = m_t*e_em*mu_h*IT_0033*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0045 = IT_0038*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0048 = IT_0038*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = 0.5*IT_0041 + (-0.5)*IT_0043 + (-0.5)*IT_0046 +
       0.5*IT_0049;
    const ccomplex_t IT_0051 = IT_0035*IT_0050;
    const ccomplex_t IT_0052 = pow(m_Z, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + IT_0052 
      + -reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = cpow(IT_0037, -1);
    const ccomplex_t IT_0057 = m_e*e_em*IT_0003*IT_0033*IT_0038*IT_0056;
    const ccomplex_t IT_0058 = (-0.5)*IT_0057;
    const ccomplex_t IT_0059 = IT_0003*IT_0036;
    const ccomplex_t IT_0060 = IT_0003*IT_0042;
    const ccomplex_t IT_0061 = IT_0037*IT_0044;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = IT_0037*IT_0047;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0059 + 0.5*IT_0060 + (-0.5)*IT_0062 +
       0.5*IT_0064;
    const ccomplex_t IT_0066 = IT_0058*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = -IT_0055 + -IT_0069;
    const ccomplex_t IT_0071 = s_13 + s_14;
    const ccomplex_t IT_0072 = m_e*IT_0071;
    const ccomplex_t IT_0073 = s_34 + IT_0014;
    const ccomplex_t IT_0074 = -IT_0029;
    const ccomplex_t IT_0075 = (-6)*IT_0074;
    const ccomplex_t IT_0076 = 12*s_34;
    const ccomplex_t IT_0077 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0078 = cos(alpha);
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0033
      *IT_0056*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = sin(alpha);
    const ccomplex_t IT_0082 = m_t*mu_h*IT_0033;
    const ccomplex_t IT_0083 = IT_0037*IT_0078;
    const ccomplex_t IT_0084 = -IT_0038*IT_0081 + IT_0083;
    const ccomplex_t IT_0085 = cpow(IT_0000, -2);
    const ccomplex_t IT_0086 = pow(m_t, 2);
    const ccomplex_t IT_0087 = (-2)*U_st_00*((0.5*m_W*IT_0084*(IT_0003 + -1./3
      *IT_0002*IT_0085) + IT_0003*IT_0033*IT_0039*IT_0081*IT_0086)*conj(U_st_01)
       + -1./2*IT_0003*IT_0039*IT_0078*IT_0082*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0039*IT_0078*IT_0082*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0084*IT_0085 + 1.5*IT_0003*IT_0033*IT_0039*IT_0081*IT_0086)*conj
      (U_st_11));
    const ccomplex_t IT_0088 = (0 + _Complex_I*-1)*(e_em*IT_0087 + (
      -1.4142135623731)*A_t*IT_0081*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = IT_0080*IT_0089;
    const ccomplex_t IT_0091 = IT_0077*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0033
      *IT_0056*IT_0081;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = IT_0037*IT_0081;
    const ccomplex_t IT_0096 = IT_0038*IT_0078;
    const ccomplex_t IT_0097 = IT_0095 + IT_0096;
    const ccomplex_t IT_0098 = 2*U_st_00*((IT_0003*IT_0033*IT_0039*IT_0078
      *IT_0086 + (-0.5)*m_W*(IT_0003 + -1./3*IT_0002*IT_0085)*IT_0097)*conj
      (U_st_01) + 1./2*IT_0003*IT_0039*IT_0081*IT_0082*conj(U_st_11)) + U_st_10*
      (IT_0003*IT_0039*IT_0081*IT_0082*conj(U_st_01) + (-1.33333333333333)*((
      -1.5)*IT_0003*IT_0033*IT_0039*IT_0078*IT_0086 + m_W*IT_0002*IT_0085
      *IT_0097)*conj(U_st_11));
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(e_em*IT_0098 +
       1.4142135623731*A_t*IT_0078*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0094*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0092 + -IT_0104;
    const ccomplex_t IT_0106 = -s_14;
    const ccomplex_t IT_0107 = s_13 + IT_0106;
    const ccomplex_t IT_0108 = m_e*IT_0107;
    const ccomplex_t IT_0109 = -IT_0014;
    const ccomplex_t IT_0110 = s_34 + IT_0109;
    const ccomplex_t IT_0111 = (-6)*conj(IT_0074);
    const ccomplex_t IT_0112 = IT_0070*IT_0073;
    const ccomplex_t IT_0113 = m_e*s_13;
    const ccomplex_t IT_0114 = s_34*IT_0007;
    const ccomplex_t IT_0115 = (-6)*IT_0114;
    const ccomplex_t IT_0116 = s_13*s_14;
    const ccomplex_t IT_0117 = 12*IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = IT_0105*IT_0108;
    const ccomplex_t IT_0120 = (-6)*IT_0031;
    const ccomplex_t IT_0121 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0122 = 6*s_34;
    const ccomplex_t IT_0123 = IT_0105*IT_0110;
    const ccomplex_t IT_0124 = (IT_0013*IT_0015 + IT_0016*IT_0029)*IT_0032 + 6
      *conj(IT_0070)*((IT_0013 + -IT_0031)*IT_0072 + IT_0073*(IT_0029 +
       0.166666666666667*IT_0075) + 0.166666666666667*IT_0070*(12*IT_0014 +
       IT_0076)) + 6*conj(IT_0105)*(0.166666666666667*((-12)*IT_0014 + IT_0076)
      *IT_0105 + (IT_0013 + IT_0031)*IT_0108 + (IT_0029 + IT_0074)*IT_0110) +
       IT_0111*(IT_0013*IT_0016 + IT_0014*IT_0029 + IT_0112) + (-6)*conj(IT_0031
      )*(IT_0070*IT_0072 + -IT_0074*IT_0113 + (-0.166666666666667)*IT_0031
      *IT_0118 + -IT_0119) + 6*conj(IT_0013)*(IT_0070*IT_0072 + IT_0029*IT_0113 
      + 0.166666666666667*IT_0013*IT_0118 + IT_0119 + 0.166666666666667*IT_0015
      *IT_0120 + IT_0075*IT_0121) + 6*conj(IT_0029)*(0.166666666666667*IT_0014
      *IT_0075 + IT_0112 + IT_0013*IT_0113 + IT_0120*IT_0121 + 0.166666666666667
      *IT_0029*IT_0122 + IT_0123) + 6*conj(IT_0074)*(IT_0031*IT_0113 +
       0.166666666666667*IT_0074*IT_0122 + IT_0123);
    return create_ccomplex_return(IT_0124);
}

