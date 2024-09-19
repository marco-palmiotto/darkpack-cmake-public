#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_anti_e_e(
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_00*conj(U_st_00) + IT_0001*IT_0002*U_st_10
      *conj(U_st_10));
    const ccomplex_t IT_0005 = 1.33333333333333*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0007 = pow(m_st_1, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = IT_0005*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = pow(m_e, 2);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = m_e*s_14;
    const ccomplex_t IT_0016 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0017 = m_e*IT_0016;
    const ccomplex_t IT_0018 = m_e*IT_0005;
    const ccomplex_t IT_0019 = IT_0001*IT_0002;
    const ccomplex_t IT_0020 = e_em*IT_0019;
    const ccomplex_t IT_0021 = IT_0000*IT_0003;
    const ccomplex_t IT_0022 = e_em*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0020 + -IT_0022);
    const ccomplex_t IT_0024 = 0.5*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = (IT_0017 + -IT_0018)*(IT_0010 + -IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0005*IT_0026;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (-6)*IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, -1);
    const ccomplex_t IT_0033 = cos(beta);
    const ccomplex_t IT_0034 = cpow(IT_0033, -1);
    const ccomplex_t IT_0035 = cos(alpha);
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0032
      *IT_0034*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = sin(alpha);
    const ccomplex_t IT_0039 = m_t*mu_h*IT_0032;
    const ccomplex_t IT_0040 = sin(beta);
    const ccomplex_t IT_0041 = cpow(IT_0040, -1);
    const ccomplex_t IT_0042 = IT_0033*IT_0035;
    const ccomplex_t IT_0043 = -IT_0038*IT_0040 + IT_0042;
    const ccomplex_t IT_0044 = cpow(IT_0000, -2);
    const ccomplex_t IT_0045 = pow(m_t, 2);
    const ccomplex_t IT_0046 = (-2)*U_st_00*((0.5*m_W*IT_0043*(IT_0003 + -1./3
      *IT_0002*IT_0044) + IT_0003*IT_0032*IT_0038*IT_0041*IT_0045)*conj(U_st_00)
       + -1./2*IT_0003*IT_0035*IT_0039*IT_0041*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0035*IT_0039*IT_0041*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0043*IT_0044 + 1.5*IT_0003*IT_0032*IT_0038*IT_0041*IT_0045)*conj
      (U_st_10));
    const ccomplex_t IT_0047 = (0 + _Complex_I*-1)*(e_em*IT_0046 + (
      -1.4142135623731)*A_t*IT_0038*(conj(U_st_00)*U_st_10 + U_st_00*conj
      (U_st_10)));
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = IT_0037*IT_0048;
    const ccomplex_t IT_0050 = cpow((-2)*s_12 + (-2)*IT_0007 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_e*e_em*IT_0003*IT_0032
      *IT_0034*IT_0038;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0033*IT_0038;
    const ccomplex_t IT_0056 = IT_0035*IT_0040;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = 2*U_st_00*((IT_0003*IT_0032*IT_0035*IT_0041
      *IT_0045 + (-0.5)*m_W*(IT_0003 + -1./3*IT_0002*IT_0044)*IT_0057)*conj
      (U_st_00) + 1./2*IT_0003*IT_0038*IT_0039*IT_0041*conj(U_st_10)) + U_st_10*
      (IT_0003*IT_0038*IT_0039*IT_0041*conj(U_st_00) + (-1.33333333333333)*((
      -1.5)*IT_0003*IT_0032*IT_0035*IT_0041*IT_0045 + m_W*IT_0002*IT_0044
      *IT_0057)*conj(U_st_10));
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(e_em*IT_0058 +
       1.4142135623731*A_t*IT_0035*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0054*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0052 + -IT_0064;
    const ccomplex_t IT_0066 = -s_14;
    const ccomplex_t IT_0067 = s_13 + IT_0066;
    const ccomplex_t IT_0068 = m_e*IT_0067;
    const ccomplex_t IT_0069 = -IT_0028;
    const ccomplex_t IT_0070 = -IT_0013;
    const ccomplex_t IT_0071 = s_34 + IT_0070;
    const ccomplex_t IT_0072 = 12*s_34;
    const ccomplex_t IT_0073 = m_e*e_em*IT_0003*IT_0032;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = m_t*e_em*mu_h*IT_0032*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0076 = IT_0003*IT_0033*IT_0041;
    const ccomplex_t IT_0077 = IT_0075*IT_0076;
    const ccomplex_t IT_0078 = m_t*e_em*mu_h*IT_0032*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0079 = IT_0076*IT_0078;
    const ccomplex_t IT_0080 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0081 = IT_0040*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0084 = IT_0040*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0077 + (-0.5)*IT_0079 + (-0.5)*IT_0082 +
       0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0074*IT_0086;
    const ccomplex_t IT_0088 = pow(m_Z, 2);
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + (-2)*IT_0007 + IT_0088 + 
      -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = m_e*e_em*IT_0003*IT_0032*IT_0034*IT_0040;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = IT_0003*IT_0075;
    const ccomplex_t IT_0095 = IT_0003*IT_0078;
    const ccomplex_t IT_0096 = IT_0033*IT_0080;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = IT_0033*IT_0083;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0094 + 0.5*IT_0095 + (-0.5)*IT_0097 +
       0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0093*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_12 + (-2)*IT_0007 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = -IT_0091 + -IT_0104;
    const ccomplex_t IT_0106 = s_34 + IT_0013;
    const ccomplex_t IT_0107 = (-6)*conj(IT_0069);
    const ccomplex_t IT_0108 = s_13 + s_14;
    const ccomplex_t IT_0109 = m_e*IT_0108;
    const ccomplex_t IT_0110 = -conj(IT_0012);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = (-1.33333333333333)*IT_0111;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = cpow(s_12 + IT_0007 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0117;
    const ccomplex_t IT_0119 = s_13*s_14;
    const ccomplex_t IT_0120 = s_34*IT_0007;
    const ccomplex_t IT_0121 = (-6)*IT_0069;
    const ccomplex_t IT_0122 = (-6)*IT_0120;
    const ccomplex_t IT_0123 = (-6)*IT_0014;
    const ccomplex_t IT_0124 = 12*IT_0119;
    const ccomplex_t IT_0125 = IT_0122 + IT_0123 + IT_0124;
    const ccomplex_t IT_0126 = IT_0118*((-12)*IT_0014 + 24*IT_0119 + (-12)
      *IT_0120) + 6*IT_0068*(IT_0028 + -conj(IT_0028) + 2*IT_0065 + (-2)*conj
      (IT_0065) + -conj(IT_0069) + (-0.166666666666667)*IT_0121) + (IT_0012 +
       IT_0030 + -conj(IT_0030) + IT_0110)*IT_0125;
    const ccomplex_t IT_0127 = conj(IT_0105)*IT_0106;
    const ccomplex_t IT_0128 = m_e*s_13;
    const ccomplex_t IT_0129 = IT_0122 + IT_0124;
    const ccomplex_t IT_0130 = conj(IT_0065)*IT_0068;
    const ccomplex_t IT_0131 = (-6)*conj(IT_0030);
    const ccomplex_t IT_0132 = 0.166666666666667*IT_0015;
    const ccomplex_t IT_0133 = 6*s_34;
    const ccomplex_t IT_0134 = conj(IT_0065)*IT_0071;
    const ccomplex_t IT_0135 = (conj(IT_0012)*IT_0014 + IT_0015*conj(IT_0028))
      *IT_0031 + 6*IT_0065*((conj(IT_0012) + conj(IT_0030))*IT_0068 + (conj
      (IT_0028) + conj(IT_0069))*IT_0071 + 0.166666666666667*conj(IT_0065)*((-12
      )*IT_0013 + IT_0072)) + 6*IT_0105*(0.166666666666667*(12*IT_0013 + IT_0072
      )*conj(IT_0105) + IT_0106*(conj(IT_0028) + 0.166666666666667*IT_0107) + 
      -IT_0109*(conj(IT_0030) + IT_0110)) + -IT_0118*IT_0126 + IT_0121*(conj
      (IT_0012)*IT_0015 + IT_0013*conj(IT_0028) + IT_0127) + (-6)*IT_0030*(conj
      (IT_0105)*IT_0109 + -conj(IT_0069)*IT_0128 + (-0.166666666666667)*conj
      (IT_0030)*IT_0129 + -IT_0130) + 6*IT_0012*(conj(IT_0105)*IT_0109 + conj
      (IT_0028)*IT_0128 + 0.166666666666667*conj(IT_0012)*IT_0129 + IT_0130 +
       0.166666666666667*IT_0014*IT_0131 + IT_0107*IT_0132) + 6*IT_0028*
      (0.166666666666667*IT_0013*IT_0107 + IT_0127 + conj(IT_0012)*IT_0128 +
       IT_0131*IT_0132 + 0.166666666666667*conj(IT_0028)*IT_0133 + IT_0134) + 6
      *IT_0069*(conj(IT_0030)*IT_0128 + 0.166666666666667*conj(IT_0069)*IT_0133 
      + IT_0134);
    return create_ccomplex_return(IT_0135);
}

