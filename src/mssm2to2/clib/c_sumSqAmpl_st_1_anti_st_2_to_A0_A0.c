#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_A0_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_A0_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t A_t = param->A_t;
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
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
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -2);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -2);
    const ccomplex_t IT_0009 = pow(m_t, 2);
    const ccomplex_t IT_0010 = cpow(IT_0003, -2);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -2);
    const ccomplex_t IT_0013 = (-0.25)*(IT_0004*(IT_0002 + (-3)*IT_0012) +
       IT_0007*(IT_0002 + (-3)*(1 + (-2)*IT_0008*IT_0009*IT_0010)*IT_0012))
      *U_st_00*conj(U_st_01) + (IT_0002*IT_0004 + IT_0007*(IT_0002 + (-1.5)
      *IT_0008*IT_0009*IT_0010*IT_0012))*U_st_10*conj(U_st_11);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0000*IT_0013;
    const ccomplex_t IT_0015 = (-0.333333333333333)*IT_0014;
    const ccomplex_t IT_0016 = IT_0002*IT_0011;
    const ccomplex_t IT_0017 = cpow(IT_0011, -1);
    const ccomplex_t IT_0018 = IT_0016 + IT_0017;
    const ccomplex_t IT_0019 = cos(alpha);
    const ccomplex_t IT_0020 = cpow(IT_0003, 3);
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_W*e_em*IT_0018*(IT_0019
      *IT_0020 + IT_0005*(IT_0004*IT_0021 + -IT_0005*(IT_0003*IT_0019 + IT_0005
      *IT_0021)));
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = pow(m_W, -1);
    const ccomplex_t IT_0025 = m_t*mu_h*IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0003, -1);
    const ccomplex_t IT_0027 = IT_0005*IT_0021;
    const ccomplex_t IT_0028 = IT_0003*IT_0019;
    const ccomplex_t IT_0029 = IT_0027 + IT_0028;
    const ccomplex_t IT_0030 = 2*U_st_00*((IT_0009*IT_0017*IT_0019*IT_0024
      *IT_0026 + 0.166666666666667*m_W*(IT_0002*IT_0011 + (-3)*IT_0017)*IT_0029)
      *conj(U_st_01) + 1./2*IT_0017*IT_0021*IT_0025*IT_0026*conj(U_st_11)) +
       U_st_10*(IT_0017*IT_0021*IT_0025*IT_0026*conj(U_st_01) + (
      -1.33333333333333)*((-1.5)*IT_0009*IT_0017*IT_0019*IT_0024*IT_0026 + m_W
      *IT_0002*IT_0011*IT_0029)*conj(U_st_11));
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(e_em*IT_0030 +
       1.4142135623731*A_t*IT_0019*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0023*IT_0032;
    const ccomplex_t IT_0034 = pow(m_st_1, 2);
    const ccomplex_t IT_0035 = pow(m_st_2, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0034 + -IT_0035 + -reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0005;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em*IT_0018*(IT_0020
      *IT_0021 + IT_0039*(IT_0004*IT_0019 + (IT_0005*IT_0019 + -IT_0003*IT_0021)
      *IT_0039));
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = IT_0005*IT_0019;
    const ccomplex_t IT_0043 = -IT_0003*IT_0021 + IT_0042;
    const ccomplex_t IT_0044 = (-2)*U_st_00*((IT_0009*IT_0017*IT_0021*IT_0024
      *IT_0026 + (-0.166666666666667)*m_W*(IT_0002*IT_0011 + (-3)*IT_0017)
      *IT_0043)*conj(U_st_01) + -1./2*IT_0017*IT_0019*IT_0025*IT_0026*conj
      (U_st_11)) + U_st_10*(IT_0017*IT_0019*IT_0025*IT_0026*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0009*IT_0017*IT_0021*IT_0024*IT_0026 + m_W
      *IT_0002*IT_0011*IT_0043)*conj(U_st_11));
    const ccomplex_t IT_0045 = (0 + _Complex_I*-1)*(e_em*IT_0044 + (
      -1.4142135623731)*A_t*IT_0021*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0041*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_12 + -IT_0034 + -IT_0035 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = m_t*e_em*mu_h*IT_0024*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0052 = IT_0017*IT_0051;
    const ccomplex_t IT_0053 = m_t*e_em*mu_h*IT_0024*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0054 = IT_0017*IT_0053;
    const ccomplex_t IT_0055 = A_t*conj(U_st_00)*U_st_10;
    const ccomplex_t IT_0056 = IT_0005*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = A_t*U_st_00*conj(U_st_10);
    const ccomplex_t IT_0059 = IT_0005*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (-0.5)*IT_0052 + 0.5*IT_0054 + (-0.5)*IT_0057 +
       0.5*IT_0060;
    const ccomplex_t IT_0062 = m_t*e_em*mu_h*IT_0024*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = m_t*e_em*mu_h*IT_0024*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0065 = IT_0017*IT_0064;
    const ccomplex_t IT_0066 = A_t*conj(U_st_01)*U_st_10;
    const ccomplex_t IT_0067 = IT_0005*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = A_t*U_st_00*conj(U_st_11);
    const ccomplex_t IT_0070 = IT_0005*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (-0.5)*IT_0063 + 0.5*IT_0065 + (-0.5)*IT_0068 +
       0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0061*IT_0072;
    const ccomplex_t IT_0074 = pow(m_A0, 2);
    const ccomplex_t IT_0075 = cpow((-2)*s_13 + IT_0074 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0035 + IT_0074 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0073*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = m_t*e_em*mu_h*IT_0024*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0082 = IT_0017*IT_0081;
    const ccomplex_t IT_0083 = m_t*e_em*mu_h*IT_0024*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0084 = IT_0017*IT_0083;
    const ccomplex_t IT_0085 = A_t*conj(U_st_01)*U_st_11;
    const ccomplex_t IT_0086 = IT_0005*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = A_t*U_st_01*conj(U_st_11);
    const ccomplex_t IT_0089 = IT_0005*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (-0.5)*IT_0082 + 0.5*IT_0084 + (-0.5)*IT_0087 +
       0.5*IT_0090;
    const ccomplex_t IT_0092 = IT_0072*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_23 + IT_0074 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_13 + IT_0034 + IT_0074 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0092*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0015 + -IT_0038 + -IT_0050 + IT_0077 +
       IT_0080 + IT_0095 + IT_0098;
    return create_ccomplex_return(3*IT_0099*conj(IT_0099));
}

