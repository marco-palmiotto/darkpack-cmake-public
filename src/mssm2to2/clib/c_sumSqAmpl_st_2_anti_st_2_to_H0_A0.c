#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_H0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_H0_A0(
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0004 + 0.5*IT_0006 + (-0.5)*IT_0010 +
       0.5*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -2);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = IT_0001 + IT_0017;
    const ccomplex_t IT_0019 = sin(beta);
    const ccomplex_t IT_0020 = cpow(IT_0019, 3);
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = -IT_0007;
    const ccomplex_t IT_0023 = cos(alpha);
    const ccomplex_t IT_0024 = cpow(IT_0019, 2);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_W*e_em*IT_0018*(IT_0020
      *IT_0021 + -IT_0022*(IT_0022*(IT_0019*IT_0021 + -IT_0007*IT_0023) + 
      -IT_0023*IT_0024));
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0014*IT_0026;
    const ccomplex_t IT_0028 = pow(m_st_2, 2);
    const ccomplex_t IT_0029 = cpow((-2)*s_12 + (-2)*IT_0028 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = m_t*e_em*mu_h*IT_0002*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = m_t*e_em*mu_h*IT_0002*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = A_t*U_st_01*conj(U_st_10);
    const ccomplex_t IT_0037 = IT_0007*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = A_t*conj(U_st_00)*U_st_11;
    const ccomplex_t IT_0040 = IT_0007*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = 0.5*IT_0033 + (-0.5)*IT_0035 + 0.5*IT_0038 + (
      -0.5)*IT_0041;
    const ccomplex_t IT_0043 = m_t*mu_h*IT_0002;
    const ccomplex_t IT_0044 = cpow(IT_0019, -1);
    const ccomplex_t IT_0045 = IT_0007*IT_0023;
    const ccomplex_t IT_0046 = -IT_0019*IT_0021 + IT_0045;
    const ccomplex_t IT_0047 = pow(m_t, 2);
    const ccomplex_t IT_0048 = (-2)*U_st_00*((0.5*m_W*(IT_0001 + -1./3*IT_0000
      *IT_0016)*IT_0046 + IT_0001*IT_0002*IT_0021*IT_0044*IT_0047)*conj(U_st_01)
       + -1./2*IT_0001*IT_0023*IT_0043*IT_0044*conj(U_st_11)) + U_st_10*(IT_0001
      *IT_0023*IT_0043*IT_0044*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0000
      *IT_0016*IT_0046 + 1.5*IT_0001*IT_0002*IT_0021*IT_0044*IT_0047)*conj
      (U_st_11));
    const ccomplex_t IT_0049 = (0 + _Complex_I*-1)*(e_em*IT_0048 + (
      -1.4142135623731)*A_t*IT_0021*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = IT_0042*IT_0050;
    const ccomplex_t IT_0052 = pow(m_H0, 2);
    const ccomplex_t IT_0053 = cpow((-2)*s_23 + IT_0028 + IT_0052 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = pow(m_Z, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + (-2)*IT_0028 + IT_0056 + 
      -reg_prop, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*-1)*m_W*e_em*IT_0007*IT_0018
      *IT_0019*(IT_0019*IT_0021 + -IT_0007*IT_0023);
    const ccomplex_t IT_0059 = IT_0001*IT_0007*IT_0044;
    const ccomplex_t IT_0060 = IT_0003*IT_0059;
    const ccomplex_t IT_0061 = IT_0005*IT_0059;
    const ccomplex_t IT_0062 = IT_0008*IT_0019;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = IT_0011*IT_0019;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0060 + (-0.5)*IT_0061 + (-0.5)*IT_0063 +
       0.5*IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = IT_0057*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0028 + -reg_prop, -1);
    const ccomplex_t IT_0071 = s_12 + IT_0028;
    const ccomplex_t IT_0072 = cpow(IT_0015, -1);
    const ccomplex_t IT_0073 = -e_em*IT_0001*(IT_0007*IT_0021 + -IT_0019
      *IT_0023)*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0001
      *IT_0015 + -1./3*IT_0000*IT_0072)*U_st_01*conj(U_st_01) + IT_0000*IT_0072
      *U_st_11*conj(U_st_11));
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = IT_0071*IT_0077;
    const ccomplex_t IT_0079 = 1.33333333333333*IT_0075;
    const ccomplex_t IT_0080 = IT_0073*IT_0079;
    const ccomplex_t IT_0081 = s_23*IT_0080;
    const ccomplex_t IT_0082 = (-2)*IT_0078 + -IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0079;
    const ccomplex_t IT_0084 = IT_0071*IT_0083;
    const ccomplex_t IT_0085 = s_13 + s_23;
    const ccomplex_t IT_0086 = IT_0073*IT_0076;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0082 + IT_0088;
    const ccomplex_t IT_0090 = IT_0070*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_13 + IT_0052 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0093 = (-2)*U_st_01*((0.5*m_W*(IT_0001 + -1./3*IT_0000
      *IT_0016)*IT_0046 + IT_0001*IT_0002*IT_0021*IT_0044*IT_0047)*conj(U_st_01)
       + -1./2*IT_0001*IT_0023*IT_0043*IT_0044*conj(U_st_11)) + U_st_11*(IT_0001
      *IT_0023*IT_0043*IT_0044*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0000
      *IT_0016*IT_0046 + 1.5*IT_0001*IT_0002*IT_0021*IT_0044*IT_0047)*conj
      (U_st_11));
    const ccomplex_t IT_0094 = (0 + _Complex_I*-1)*(e_em*IT_0093 + (
      -1.4142135623731)*A_t*IT_0021*(conj(U_st_01)*U_st_11 + U_st_01*conj
      (U_st_11)));
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = IT_0014*IT_0095;
    const ccomplex_t IT_0097 = IT_0092*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = cpow((-2)*s_23 + IT_0052 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0100 = IT_0096*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_13 + IT_0028 + IT_0052 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0103 = m_t*e_em*mu_h*IT_0002*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0104 = IT_0001*IT_0103;
    const ccomplex_t IT_0105 = m_t*e_em*mu_h*IT_0002*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0106 = IT_0001*IT_0105;
    const ccomplex_t IT_0107 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0108 = IT_0007*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0104 + 0.5*IT_0106 + (-0.5)*IT_0109 +
       0.5*IT_0112;
    const ccomplex_t IT_0114 = (-2)*conj(U_st_00)*((0.5*m_W*(IT_0001 + -1./3
      *IT_0000*IT_0016)*IT_0046 + IT_0001*IT_0002*IT_0021*IT_0044*IT_0047)
      *U_st_01 + -1./2*IT_0001*IT_0023*IT_0043*IT_0044*U_st_11) + conj(U_st_10)*
      (IT_0001*IT_0023*IT_0043*IT_0044*U_st_01 + (-1.33333333333333)*(m_W
      *IT_0000*IT_0016*IT_0046 + 1.5*IT_0001*IT_0002*IT_0021*IT_0044*IT_0047)
      *U_st_11);
    const ccomplex_t IT_0115 = (0 + _Complex_I*-1)*(e_em*IT_0114 + (
      -1.4142135623731)*A_t*IT_0021*(U_st_01*conj(U_st_10) + conj(U_st_00)
      *U_st_11));
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0113*IT_0116;
    const ccomplex_t IT_0118 = IT_0102*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0031 + IT_0055 + -IT_0069 + IT_0091 +
       IT_0098 + IT_0101 + IT_0119;
    return create_ccomplex_return(3*IT_0120*conj(IT_0120));
}

