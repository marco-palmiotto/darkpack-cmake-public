#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_h_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_h_h(
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
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, -2);
    const ccomplex_t IT_0009 = pow(m_t, 2);
    const ccomplex_t IT_0010 = sin(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -2);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -2);
    const ccomplex_t IT_0014 = (-0.25)*(IT_0004*(IT_0002 + (-3)*IT_0013) +
       IT_0007*(IT_0002 + (-3)*(1 + (-2)*IT_0008*IT_0009*IT_0011)*IT_0013))
      *U_st_00*conj(U_st_01) + (IT_0002*IT_0004 + IT_0007*(IT_0002 + (-1.5)
      *IT_0008*IT_0009*IT_0011*IT_0013))*U_st_10*conj(U_st_11);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0000*IT_0014;
    const ccomplex_t IT_0016 = (-0.333333333333333)*IT_0015;
    const ccomplex_t IT_0017 = pow(m_st_1, 2);
    const ccomplex_t IT_0018 = pow(m_st_2, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + -IT_0017 + -IT_0018 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0002*IT_0012;
    const ccomplex_t IT_0021 = cpow(IT_0012, -1);
    const ccomplex_t IT_0022 = IT_0020 + IT_0021;
    const ccomplex_t IT_0023 = cpow(IT_0003, 3);
    const ccomplex_t IT_0024 = cos(beta);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0010
      *IT_0023 + (-5)*IT_0005*(IT_0003*IT_0005*IT_0010 + 1*(IT_0004 + -1./5
      *IT_0006)*IT_0024));
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = m_t*mu_h*IT_0027;
    const ccomplex_t IT_0029 = cpow(IT_0010, -1);
    const ccomplex_t IT_0030 = IT_0005*IT_0024;
    const ccomplex_t IT_0031 = -IT_0003*IT_0010 + IT_0030;
    const ccomplex_t IT_0032 = (-2)*U_st_00*((IT_0003*IT_0009*IT_0021*IT_0027
      *IT_0029 + (-0.166666666666667)*m_W*(IT_0002*IT_0012 + (-3)*IT_0021)
      *IT_0031)*conj(U_st_01) + -1./2*IT_0005*IT_0021*IT_0028*IT_0029*conj
      (U_st_11)) + U_st_10*(IT_0005*IT_0021*IT_0028*IT_0029*conj(U_st_01) + (
      -1.33333333333333)*(1.5*IT_0003*IT_0009*IT_0021*IT_0027*IT_0029 + m_W
      *IT_0002*IT_0012*IT_0031)*conj(U_st_11));
    const ccomplex_t IT_0033 = (0 + _Complex_I*-1)*(e_em*IT_0032 + (
      -1.4142135623731)*A_t*IT_0003*(conj(U_st_01)*U_st_10 + U_st_00*conj
      (U_st_11)));
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = IT_0026*IT_0034;
    const ccomplex_t IT_0036 = IT_0019*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*m_W*e_em*IT_0022*(IT_0023
      *IT_0024 + IT_0005*(IT_0004*IT_0010 + -IT_0005*(IT_0005*IT_0010 + IT_0003
      *IT_0024)));
    const ccomplex_t IT_0039 = 1.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0003*IT_0024;
    const ccomplex_t IT_0041 = IT_0005*IT_0010;
    const ccomplex_t IT_0042 = IT_0040 + IT_0041;
    const ccomplex_t IT_0043 = 2*U_st_00*((IT_0005*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.166666666666667*m_W*(IT_0002*IT_0012 + (-3)*IT_0021)*IT_0042)
      *conj(U_st_01) + 1./2*IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_11)) +
       U_st_10*(IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_01) + (
      -1.33333333333333)*((-1.5)*IT_0005*IT_0009*IT_0021*IT_0027*IT_0029 + m_W
      *IT_0002*IT_0012*IT_0042)*conj(U_st_11));
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(e_em*IT_0043 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = IT_0039*IT_0045;
    const ccomplex_t IT_0047 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0017 + -IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = pow(m_h, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0050 + (0 + _Complex_I*1)
      *m_st_1*Gamma_t1 + reg_prop, -1);
    const ccomplex_t IT_0052 = 2*U_st_00*((IT_0005*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.166666666666667*m_W*(IT_0002*IT_0012 + (-3)*IT_0021)*IT_0042)
      *conj(U_st_00) + 1./2*IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_10)) +
       U_st_10*(IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_00) + (
      -1.33333333333333)*((-1.5)*IT_0005*IT_0009*IT_0021*IT_0027*IT_0029 + m_W
      *IT_0002*IT_0012*IT_0042)*conj(U_st_10));
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(e_em*IT_0052 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_00)*U_st_10 + U_st_00*conj(U_st_10
      )));
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0045*IT_0054;
    const ccomplex_t IT_0056 = IT_0051*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0018 + IT_0050 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0055*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_13 + IT_0017 + IT_0050 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0062 = 2*U_st_01*((IT_0005*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.166666666666667*m_W*(IT_0002*IT_0012 + (-3)*IT_0021)*IT_0042)
      *conj(U_st_01) + 1./2*IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_11)) +
       U_st_11*(IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_01) + (
      -1.33333333333333)*((-1.5)*IT_0005*IT_0009*IT_0021*IT_0027*IT_0029 + m_W
      *IT_0002*IT_0012*IT_0042)*conj(U_st_11));
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(e_em*IT_0062 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_01)*U_st_11 + U_st_01*conj(U_st_11
      )));
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = 2*U_st_00*((IT_0005*IT_0009*IT_0021*IT_0027
      *IT_0029 + 0.166666666666667*m_W*(IT_0002*IT_0012 + (-3)*IT_0021)*IT_0042)
      *conj(U_st_01) + 1./2*IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_11)) +
       U_st_10*(IT_0003*IT_0021*IT_0028*IT_0029*conj(U_st_01) + (
      -1.33333333333333)*((-1.5)*IT_0005*IT_0009*IT_0021*IT_0027*IT_0029 + m_W
      *IT_0002*IT_0012*IT_0042)*conj(U_st_11));
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(e_em*IT_0065 +
       1.4142135623731*A_t*IT_0005*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0064*IT_0067;
    const ccomplex_t IT_0069 = IT_0061*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = cpow((-2)*s_23 + IT_0050 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0068*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0016 + -IT_0037 + -IT_0049 + IT_0057 +
       IT_0060 + IT_0070 + IT_0073;
    return create_ccomplex_return(3*IT_0074*conj(IT_0074));
}

