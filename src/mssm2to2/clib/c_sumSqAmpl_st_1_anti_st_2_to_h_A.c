#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_2_to_h_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_2_to_h_A(
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
    const creal_t mu_h = param->mu_h;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = (-1.33333333333333)*IT_0000;
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = m_t*mu_h*IT_0003;
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = IT_0002*IT_0005;
    const ccomplex_t IT_0013 = IT_0011 + IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -2);
    const ccomplex_t IT_0016 = pow(m_t, 2);
    const ccomplex_t IT_0017 = 2*U_st_00*(((-0.5)*m_W*IT_0013*(IT_0009 + -1./3
      *IT_0008*IT_0015) + IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)*conj(U_st_01)
       + 1./2*IT_0004*IT_0006*IT_0007*IT_0009*conj(U_st_11)) + U_st_10*(IT_0004
      *IT_0006*IT_0007*IT_0009*conj(U_st_01) + (-1.33333333333333)*(m_W*IT_0008
      *IT_0013*IT_0015 + (-1.5)*IT_0002*IT_0003*IT_0006*IT_0009*IT_0016)*conj
      (U_st_11));
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(e_em*IT_0017 +
       1.4142135623731*A_t*IT_0002*(conj(U_st_01)*U_st_10 + U_st_00*conj(U_st_11
      )));
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = pow(m_h, 2);
    const ccomplex_t IT_0021 = pow(m_st_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_23 + IT_0020 + IT_0021 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = pow(m_st_1, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0020 + IT_0025 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0024 + IT_0028;
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = IT_0001*IT_0024;
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (-3)*s_23;
    const ccomplex_t IT_0034 = (-0.666666666666667)*IT_0000;
    const ccomplex_t IT_0035 = IT_0029*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = (-3)*s_24;
    const ccomplex_t IT_0038 = (-3)*s_34;
    const ccomplex_t IT_0039 = (-3)*IT_0030*(IT_0021*conj(IT_0030) + (
      -0.333333333333333)*conj(IT_0032)*IT_0033 + (-0.333333333333333)*conj
      (IT_0036)*IT_0037) + IT_0036*(conj(IT_0030)*IT_0037 + conj(IT_0032)
      *IT_0038) + (-3)*IT_0032*(IT_0020*conj(IT_0032) + (-0.333333333333333)
      *conj(IT_0030)*IT_0033 + (-0.333333333333333)*conj(IT_0036)*IT_0038);
    return create_ccomplex_return(IT_0039);
}

