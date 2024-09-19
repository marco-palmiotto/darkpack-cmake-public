#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_2_anti_st_2_to_A_A.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_2_anti_st_2_to_A_A(
        cparam_t const *param
        )
{
    clearcache();
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
    const ccomplex_t IT_0000 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_st_2
      *Gamma_t2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*IT_0000;
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0003 = (-0.666666666666667)*IT_0002;
    const ccomplex_t IT_0004 = (-1.33333333333333)*IT_0002;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0003, 2);
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = IT_0005 + IT_0007;
    const ccomplex_t IT_0009 = IT_0001*IT_0008;
    const ccomplex_t IT_0010 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_st_2
      *Gamma_t2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*IT_0010;
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0009 + (-0.5)*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0004, 2);
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = 3*s_13;
    const ccomplex_t IT_0018 = s_24*IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0014;
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = 3*s_23;
    const ccomplex_t IT_0022 = s_24*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0019;
    const ccomplex_t IT_0024 = s_34*IT_0021;
    const ccomplex_t IT_0025 = IT_0001*IT_0005;
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = 3*s_12;
    const ccomplex_t IT_0028 = s_24*IT_0027;
    const ccomplex_t IT_0029 = s_34*IT_0027;
    const ccomplex_t IT_0030 = IT_0005*IT_0011;
    const ccomplex_t IT_0031 = -IT_0005;
    const ccomplex_t IT_0032 = IT_0014 + IT_0031;
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0030 + (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0030;
    const ccomplex_t IT_0036 = s_34*IT_0017;
    const ccomplex_t IT_0037 = pow(m_st_2, 2);
    const ccomplex_t IT_0038 = 3*s_24;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = pow(e_em, 2);
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = tan(theta_W);
    const ccomplex_t IT_0044 = cpow(IT_0043, 2);
    const ccomplex_t IT_0045 = cpow(IT_0041 + IT_0041*IT_0044, -1);
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0040*IT_0042*IT_0045;
    const ccomplex_t IT_0047 = 0.888888888888889*IT_0046;
    const ccomplex_t IT_0048 = 0.5*IT_0030;
    const ccomplex_t IT_0049 = 3*s_34;
    const ccomplex_t IT_0050 = (-3)*s_34*IT_0013 + (-3)*s_12*IT_0016 + (-3)
      *s_14*(IT_0026 + -conj(IT_0026)) + conj(IT_0016)*IT_0027 + (-3)*(IT_0020 +
       -conj(IT_0020))*IT_0037 + (-12)*IT_0047 + (-3)*s_24*IT_0048 + IT_0038
      *conj(IT_0048) + conj(IT_0013)*IT_0049;
    const ccomplex_t IT_0051 = IT_0037*IT_0049;
    const ccomplex_t IT_0052 = IT_0017*IT_0037;
    const ccomplex_t IT_0053 = pow(m_st_2, 4);
    const ccomplex_t IT_0054 = 3*IT_0053;
    const ccomplex_t IT_0055 = IT_0023 + IT_0034;
    const ccomplex_t IT_0056 = IT_0021*IT_0037;
    const ccomplex_t IT_0057 = pow(s_23, 2);
    const ccomplex_t IT_0058 = 3*IT_0057;
    const ccomplex_t IT_0059 = (-3)*s_23*IT_0047;
    const ccomplex_t IT_0060 = s_23*IT_0016;
    const ccomplex_t IT_0061 = IT_0021*IT_0047;
    const ccomplex_t IT_0062 = IT_0020*IT_0056;
    const ccomplex_t IT_0063 = conj(IT_0013)*(IT_0016*IT_0018 + IT_0020
      *IT_0022 + IT_0023*IT_0024) + conj(IT_0026)*(IT_0020*IT_0028 + IT_0023
      *IT_0029 + IT_0018*IT_0034 + IT_0035*IT_0036 + IT_0016*IT_0039) + IT_0047
      *IT_0050 + conj(IT_0048)*(IT_0016*IT_0028 + IT_0022*IT_0034 + IT_0024
      *IT_0035 + IT_0020*IT_0039 + IT_0023*IT_0051) + conj(IT_0016)*(IT_0013
      *IT_0018 + s_23*IT_0017*IT_0035 + IT_0027*(s_23*IT_0023 + IT_0020*IT_0037)
       + IT_0026*IT_0039 + IT_0028*IT_0048 + IT_0034*IT_0052 + IT_0016*IT_0054) 
      + conj(IT_0020)*(IT_0013*IT_0022 + IT_0026*IT_0028 + IT_0016*IT_0027
      *IT_0037 + IT_0039*IT_0048 + IT_0020*IT_0054 + IT_0055*IT_0056 + IT_0035
      *IT_0058) + IT_0055*IT_0059 + conj(IT_0035)*(IT_0026*IT_0036 + IT_0024
      *IT_0048 + IT_0020*IT_0058 + IT_0017*IT_0060) + conj(IT_0034)*(IT_0018
      *IT_0026 + IT_0022*IT_0048 + IT_0016*IT_0052 + IT_0023*IT_0058 + IT_0061 +
       IT_0062) + conj(IT_0023)*(IT_0013*IT_0024 + IT_0026*IT_0029 + IT_0048
      *IT_0051 + IT_0034*IT_0058 + IT_0027*IT_0060 + IT_0061 + IT_0062);
    return create_ccomplex_return(IT_0063);
}

