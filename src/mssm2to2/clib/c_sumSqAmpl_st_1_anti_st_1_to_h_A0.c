#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_h_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_h_A0(
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
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
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
    const ccomplex_t IT_0000 = pow(m_h, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = m_t*e_em*mu_h*IT_0004*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = m_t*e_em*mu_h*IT_0004*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0006 + 0.5*IT_0008 + (-0.5)*IT_0012 +
       0.5*IT_0015;
    const ccomplex_t IT_0017 = cos(alpha);
    const ccomplex_t IT_0018 = m_t*mu_h*IT_0004;
    const ccomplex_t IT_0019 = sin(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = IT_0009*IT_0021;
    const ccomplex_t IT_0023 = IT_0017*IT_0019;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = cos(theta_W);
    const ccomplex_t IT_0026 = cpow(IT_0025, -2);
    const ccomplex_t IT_0027 = pow(m_t, 2);
    const ccomplex_t IT_0028 = 2*U_st_00*(((-0.5)*m_W*IT_0024*(IT_0003 + -1./3
      *IT_0002*IT_0026) + IT_0003*IT_0004*IT_0017*IT_0020*IT_0027)*conj(U_st_00)
       + 1./2*IT_0003*IT_0018*IT_0020*IT_0021*conj(U_st_10)) + U_st_10*(IT_0003
      *IT_0018*IT_0020*IT_0021*conj(U_st_00) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0024*IT_0026 + (-1.5)*IT_0003*IT_0004*IT_0017*IT_0020*IT_0027)*conj
      (U_st_10));
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(e_em*IT_0028 +
       1.4142135623731*A_t*IT_0017*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = IT_0016*IT_0030;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0002*IT_0026;
    const ccomplex_t IT_0035 = IT_0003 + IT_0034;
    const ccomplex_t IT_0036 = cpow(IT_0019, 3);
    const ccomplex_t IT_0037 = cpow(IT_0019, 2);
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_W*e_em*IT_0035*(IT_0017
      *IT_0036 + -IT_0009*(IT_0009*(IT_0017*IT_0019 + IT_0009*IT_0021) + 
      -IT_0021*IT_0037));
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0016*IT_0039;
    const ccomplex_t IT_0041 = pow(m_st_1, 2);
    const ccomplex_t IT_0042 = cpow((-2)*s_12 + (-2)*IT_0041 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*m_W*e_em*IT_0009*IT_0019*
      (IT_0017*IT_0019 + IT_0009*IT_0021)*IT_0035;
    const ccomplex_t IT_0046 = -IT_0045;
    const ccomplex_t IT_0047 = IT_0003*IT_0009*IT_0020;
    const ccomplex_t IT_0048 = IT_0005*IT_0047;
    const ccomplex_t IT_0049 = IT_0007*IT_0047;
    const ccomplex_t IT_0050 = IT_0010*IT_0019;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = IT_0013*IT_0019;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0048 + (-0.5)*IT_0049 + (-0.5)*IT_0051 +
       0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0046*IT_0054;
    const ccomplex_t IT_0056 = pow(m_Z, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + (-2)*IT_0041 + IT_0056 + 
      -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = cpow(IT_0025, -1);
    const ccomplex_t IT_0061 = e_em*IT_0003*(IT_0009*IT_0017 + IT_0019*IT_0021
      )*IT_0060;
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0003
      *IT_0025 + -1./3*IT_0002*IT_0060)*U_st_00*conj(U_st_00) + IT_0002*IT_0060
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0064 = 0.666666666666667*IT_0063;
    const ccomplex_t IT_0065 = IT_0062*IT_0064;
    const ccomplex_t IT_0066 = s_12 + IT_0041;
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = -IT_0061;
    const ccomplex_t IT_0069 = IT_0064*IT_0068;
    const ccomplex_t IT_0070 = s_13 + s_23;
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = 1.33333333333333*IT_0063;
    const ccomplex_t IT_0073 = IT_0062*IT_0072;
    const ccomplex_t IT_0074 = IT_0066*IT_0073;
    const ccomplex_t IT_0075 = IT_0068*IT_0072;
    const ccomplex_t IT_0076 = s_23*IT_0075;
    const ccomplex_t IT_0077 = (-0.5)*IT_0071 + (-0.5)*IT_0074 + 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0067 + IT_0077;
    const ccomplex_t IT_0079 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0041 + -reg_prop, -1);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0031*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_13 + IT_0000 + IT_0041 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0086 = m_t*e_em*mu_h*IT_0004*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0087 = IT_0003*IT_0086;
    const ccomplex_t IT_0088 = m_t*e_em*mu_h*IT_0004*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = A_t*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0091 = IT_0009*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = A_t*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0094 = IT_0009*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0087 + (-0.5)*IT_0089 + 0.5*IT_0092 + (
      -0.5)*IT_0095;
    const ccomplex_t IT_0097 = 2*U_st_00*(((-0.5)*m_W*IT_0024*(IT_0003 + -1./3
      *IT_0002*IT_0026) + IT_0003*IT_0004*IT_0017*IT_0020*IT_0027)*conj(U_st_01)
       + 1./2*IT_0003*IT_0018*IT_0020*IT_0021*conj(U_st_11)) + U_st_10*(IT_0003
      *IT_0018*IT_0020*IT_0021*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0002
      *IT_0024*IT_0026 + (-1.5)*IT_0003*IT_0004*IT_0017*IT_0020*IT_0027)*conj
      (U_st_11));
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(e_em*IT_0097 +
       1.4142135623731*A_t*IT_0017*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0096*IT_0099;
    const ccomplex_t IT_0101 = IT_0085*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0000 + IT_0041 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0104 = m_t*e_em*mu_h*IT_0004*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0105 = IT_0003*IT_0104;
    const ccomplex_t IT_0106 = m_t*e_em*mu_h*IT_0004*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0107 = IT_0003*IT_0106;
    const ccomplex_t IT_0108 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0109 = IT_0009*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0112 = IT_0009*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0105 + 0.5*IT_0107 + (-0.5)*IT_0110 +
       0.5*IT_0113;
    const ccomplex_t IT_0115 = 2*conj(U_st_00)*(((-0.5)*m_W*IT_0024*(IT_0003 +
       -1./3*IT_0002*IT_0026) + IT_0003*IT_0004*IT_0017*IT_0020*IT_0027)*U_st_01
       + 1./2*IT_0003*IT_0018*IT_0020*IT_0021*U_st_11) + conj(U_st_10)*(IT_0003
      *IT_0018*IT_0020*IT_0021*U_st_01 + (-1.33333333333333)*(m_W*IT_0002
      *IT_0024*IT_0026 + (-1.5)*IT_0003*IT_0004*IT_0017*IT_0020*IT_0027)*U_st_11);
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(e_em*IT_0115 +
       1.4142135623731*A_t*IT_0017*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = IT_0114*IT_0117;
    const ccomplex_t IT_0119 = IT_0103*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = IT_0033 + -IT_0044 + -IT_0059 + (-2)*IT_0081 +
       IT_0084 + IT_0102 + IT_0120;
    return create_ccomplex_return(3*IT_0121*conj(IT_0121));
}

