#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_Z_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_Z_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_s*e_em*(0.25*(IT_0001
      *IT_0002 + (-3)*IT_0000*IT_0003)*U_st_01*conj(U_st_01) + IT_0001*IT_0002
      *U_st_11*conj(U_st_11));
    const ccomplex_t IT_0005 = (-1.33333333333333)*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_st_2, 2);
    const ccomplex_t IT_0008 = pow(m_Z, -2);
    const ccomplex_t IT_0009 = pow(s_13, 2);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = IT_0007 + IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0014 = -IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0001*IT_0002 
      + (-3)*IT_0000*IT_0003)*U_st_01*conj(U_st_01) + IT_0001*IT_0002*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0016 = 1.33333333333333*IT_0015;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = pow(m_Z, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0018 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (-2)*IT_0013;
    const ccomplex_t IT_0023 = IT_0016*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = 4*s_24;
    const ccomplex_t IT_0027 = conj(IT_0025)*IT_0026;
    const ccomplex_t IT_0028 = 4*IT_0007;
    const ccomplex_t IT_0029 = s_13*s_23*IT_0008;
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = s_12 + IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0018 + (0 + _Complex_I*1)
      *m_st_2*Gamma_t2 + reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0023*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = 4*s_23;
    const ccomplex_t IT_0037 = (-4)*IT_0005;
    const ccomplex_t IT_0038 = IT_0028*conj(IT_0034);
    const ccomplex_t IT_0039 = IT_0017*IT_0033;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = IT_0026*conj(IT_0040);
    const ccomplex_t IT_0042 = conj(IT_0035)*IT_0036 + IT_0037 + IT_0038 +
       IT_0041;
    const ccomplex_t IT_0043 = 4*s_34;
    const ccomplex_t IT_0044 = conj(IT_0035)*IT_0043;
    const ccomplex_t IT_0045 = IT_0026*conj(IT_0034) + IT_0044;
    const ccomplex_t IT_0046 = 4*IT_0005;
    const ccomplex_t IT_0047 = pow(s_23, 2);
    const ccomplex_t IT_0048 = IT_0008*IT_0047;
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = IT_0007 + IT_0049;
    const ccomplex_t IT_0051 = s_13*s_34*IT_0008;
    const ccomplex_t IT_0052 = -IT_0051;
    const ccomplex_t IT_0053 = s_14 + IT_0052;
    const ccomplex_t IT_0054 = s_23*s_34*IT_0008;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = s_24 + IT_0055;
    const ccomplex_t IT_0057 = (-12)*IT_0006 + IT_0012*(IT_0021*IT_0027 +
       IT_0025*(conj(IT_0021)*IT_0026 + conj(IT_0025)*IT_0028)) + IT_0031*
      (IT_0027*IT_0040 + IT_0025*IT_0042 + conj(IT_0021)*(IT_0026*IT_0034 +
       IT_0035*IT_0043) + IT_0021*IT_0045 + conj(IT_0025)*(IT_0028*IT_0034 +
       IT_0035*IT_0036 + IT_0046)) + (IT_0034*IT_0042 + 4*IT_0035*(IT_0018*conj
      (IT_0035) + 0.25*conj(IT_0034)*IT_0036 + 0.25*conj(IT_0040)*IT_0043) +
       IT_0040*IT_0045 + conj(IT_0034)*IT_0046)*IT_0050 + (IT_0021*IT_0037 +
       conj(IT_0021)*IT_0046)*IT_0053 + (IT_0037*IT_0040 + conj(IT_0040)*IT_0046
      )*IT_0056;
    return create_ccomplex_return(IT_0057);
}

