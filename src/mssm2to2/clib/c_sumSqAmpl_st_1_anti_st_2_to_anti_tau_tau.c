#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_anti_tau_tau(
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
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_tau = param->m_tau;
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
    const ccomplex_t IT_0014 = pow(m_tau, 2);
    const ccomplex_t IT_0015 = IT_0007*IT_0014;
    const ccomplex_t IT_0016 = s_14*m_tau;
    const ccomplex_t IT_0017 = 0.666666666666667*IT_0004;
    const ccomplex_t IT_0018 = m_tau*IT_0017;
    const ccomplex_t IT_0019 = m_tau*IT_0005;
    const ccomplex_t IT_0020 = IT_0001*IT_0002;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = IT_0000*IT_0003;
    const ccomplex_t IT_0023 = e_em*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0021 + -IT_0023);
    const ccomplex_t IT_0025 = 0.5*IT_0024;
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = -(IT_0018 + -IT_0019)*(IT_0011 + -IT_0027);
    const ccomplex_t IT_0029 = IT_0005*IT_0027;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = (-6)*IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, -1);
    const ccomplex_t IT_0033 = e_em*m_tau*IT_0003*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = m_t*e_em*mu_h*IT_0032*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = sin(beta);
    const ccomplex_t IT_0038 = cpow(IT_0037, -1);
    const ccomplex_t IT_0039 = IT_0003*IT_0036*IT_0038;
    const ccomplex_t IT_0040 = IT_0035*IT_0039;
    const ccomplex_t IT_0041 = m_t*e_em*mu_h*IT_0032*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0042 = IT_0039*IT_0041;
    const ccomplex_t IT_0043 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0044 = IT_0037*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0047 = IT_0037*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0040 + (-0.5)*IT_0042 + (-0.5)*IT_0045 +
       0.5*IT_0048;
    const ccomplex_t IT_0050 = IT_0034*IT_0049;
    const ccomplex_t IT_0051 = pow(m_Z, 2);
    const ccomplex_t IT_0052 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + IT_0051 
      + -reg_prop, -1);
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow(IT_0036, -1);
    const ccomplex_t IT_0056 = e_em*m_tau*IT_0003*IT_0032*IT_0037*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0056;
    const ccomplex_t IT_0058 = IT_0003*IT_0035;
    const ccomplex_t IT_0059 = IT_0003*IT_0041;
    const ccomplex_t IT_0060 = IT_0036*IT_0043;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = IT_0036*IT_0046;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0058 + 0.5*IT_0059 + (-0.5)*IT_0061 +
       0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0057*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = -IT_0054 + -IT_0068;
    const ccomplex_t IT_0070 = s_13 + s_14;
    const ccomplex_t IT_0071 = m_tau*IT_0070;
    const ccomplex_t IT_0072 = s_34 + IT_0014;
    const ccomplex_t IT_0073 = -IT_0028;
    const ccomplex_t IT_0074 = (-6)*conj(IT_0073);
    const ccomplex_t IT_0075 = 12*s_34;
    const ccomplex_t IT_0076 = sin(alpha);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*e_em*m_tau*IT_0003*IT_0032
      *IT_0055*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = cos(alpha);
    const ccomplex_t IT_0080 = m_t*mu_h*IT_0032;
    const ccomplex_t IT_0081 = IT_0036*IT_0076;
    const ccomplex_t IT_0082 = IT_0037*IT_0079;
    const ccomplex_t IT_0083 = IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = cpow(IT_0000, -2);
    const ccomplex_t IT_0085 = pow(m_t, 2);
    const ccomplex_t IT_0086 = 2*U_st_00*(((-0.5)*m_W*IT_0083*(IT_0003 + -1./3
      *IT_0002*IT_0084) + IT_0003*IT_0032*IT_0038*IT_0079*IT_0085)*conj(U_st_01)
       + 1./2*IT_0003*IT_0038*IT_0076*IT_0080*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0038*IT_0076*IT_0080*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0083*IT_0084 + (-1.5)*IT_0003*IT_0032*IT_0038*IT_0079*IT_0085)*conj
      (U_st_11));
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*(e_em*IT_0086 +
       1.4142135623731*A_t*IT_0079*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0078*IT_0088;
    const ccomplex_t IT_0090 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*e_em*m_tau*IT_0003*IT_0032
      *IT_0055*IT_0079;
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0036*IT_0079;
    const ccomplex_t IT_0097 = -IT_0037*IT_0076 + IT_0096;
    const ccomplex_t IT_0098 = (-2)*U_st_00*((IT_0003*IT_0032*IT_0038*IT_0076
      *IT_0085 + 0.5*m_W*(IT_0003 + -1./3*IT_0002*IT_0084)*IT_0097)*conj(U_st_01
      ) + -1./2*IT_0003*IT_0038*IT_0079*IT_0080*conj(U_st_11)) + U_st_10*
      (IT_0003*IT_0038*IT_0079*IT_0080*conj(U_st_01) + (-1.33333333333333)*(1.5
      *IT_0003*IT_0032*IT_0038*IT_0076*IT_0085 + m_W*IT_0002*IT_0084*IT_0097)
      *conj(U_st_11));
    const ccomplex_t IT_0099 = (0 + _Complex_I*-1)*(e_em*IT_0098 + (
      -1.4142135623731)*A_t*IT_0076*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0095*IT_0100;
    const ccomplex_t IT_0102 = IT_0093*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = -IT_0092 + -IT_0103;
    const ccomplex_t IT_0105 = -s_14;
    const ccomplex_t IT_0106 = s_13 + IT_0105;
    const ccomplex_t IT_0107 = m_tau*IT_0106;
    const ccomplex_t IT_0108 = -IT_0014;
    const ccomplex_t IT_0109 = s_34 + IT_0108;
    const ccomplex_t IT_0110 = s_13*m_tau;
    const ccomplex_t IT_0111 = s_34*IT_0007;
    const ccomplex_t IT_0112 = (-6)*IT_0111;
    const ccomplex_t IT_0113 = s_13*s_14;
    const ccomplex_t IT_0114 = 12*IT_0113;
    const ccomplex_t IT_0115 = IT_0112 + IT_0114;
    const ccomplex_t IT_0116 = conj(IT_0104)*IT_0107;
    const ccomplex_t IT_0117 = (-6)*conj(IT_0030);
    const ccomplex_t IT_0118 = 0.166666666666667*IT_0016;
    const ccomplex_t IT_0119 = 6*s_34;
    const ccomplex_t IT_0120 = conj(IT_0104)*IT_0109;
    const ccomplex_t IT_0121 = (-6)*IT_0073;
    const ccomplex_t IT_0122 = conj(IT_0069)*IT_0072;
    const ccomplex_t IT_0123 = (conj(IT_0013)*IT_0015 + IT_0016*conj(IT_0028))
      *IT_0031 + 6*IT_0069*((conj(IT_0013) + -conj(IT_0030))*IT_0071 + IT_0072*
      (conj(IT_0028) + 0.166666666666667*IT_0074) + 0.166666666666667*conj
      (IT_0069)*(12*IT_0014 + IT_0075)) + 6*IT_0104*(0.166666666666667*((-12)
      *IT_0014 + IT_0075)*conj(IT_0104) + (conj(IT_0013) + conj(IT_0030))
      *IT_0107 + (conj(IT_0028) + conj(IT_0073))*IT_0109) + (-6)*IT_0030*(conj
      (IT_0069)*IT_0071 + -conj(IT_0073)*IT_0110 + (-0.166666666666667)*conj
      (IT_0030)*IT_0115 + -IT_0116) + 6*IT_0013*(conj(IT_0069)*IT_0071 + conj
      (IT_0028)*IT_0110 + 0.166666666666667*conj(IT_0013)*IT_0115 + IT_0116 +
       0.166666666666667*IT_0015*IT_0117 + IT_0074*IT_0118) + 6*IT_0073*(conj
      (IT_0030)*IT_0110 + 0.166666666666667*conj(IT_0073)*IT_0119 + IT_0120) +
       IT_0121*(conj(IT_0013)*IT_0016 + IT_0014*conj(IT_0028) + IT_0122) + 6
      *IT_0028*(0.166666666666667*IT_0014*IT_0074 + conj(IT_0013)*IT_0110 +
       IT_0117*IT_0118 + 0.166666666666667*conj(IT_0028)*IT_0119 + IT_0120 +
       IT_0122);
    return create_ccomplex_return(IT_0123);
}

