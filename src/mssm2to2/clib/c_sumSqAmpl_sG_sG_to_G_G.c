#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sG_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sG_to_G_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = -g_s;
    const ccomplex_t IT_0001 = pow(m_sG, 2);
    const ccomplex_t IT_0002 = cpow(s_12 + IT_0001 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*IT_0003;
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = g_s*IT_0004;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = 1152*s_12;
    const ccomplex_t IT_0010 = s_24*IT_0008*IT_0009;
    const ccomplex_t IT_0011 = 2304*s_14*IT_0008;
    const ccomplex_t IT_0012 = 2304*s_24*IT_0008;
    const ccomplex_t IT_0013 = s_14*IT_0001;
    const ccomplex_t IT_0014 = (g_s + -IT_0000)*IT_0004;
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = cpow(IT_0000, 2);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(s_23 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = 1152*IT_0015*IT_0020;
    const ccomplex_t IT_0022 = pow(s_13 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0019;
    const ccomplex_t IT_0026 = 1152*IT_0025;
    const ccomplex_t IT_0027 = (-1152)*IT_0008;
    const ccomplex_t IT_0028 = s_14*s_24;
    const ccomplex_t IT_0029 = cpow(IT_0020, 2);
    const ccomplex_t IT_0030 = cpow(IT_0024, 2);
    const ccomplex_t IT_0031 = (-1152)*IT_0015*IT_0024;
    const ccomplex_t IT_0032 = s_24*IT_0001;
    const ccomplex_t IT_0033 = s_12*s_14;
    const ccomplex_t IT_0034 = (-4)*IT_0015;
    const ccomplex_t IT_0035 = s_12*s_24;
    const ccomplex_t IT_0036 = (-288)*IT_0035;
    const ccomplex_t IT_0037 = 144*IT_0013;
    const ccomplex_t IT_0038 = IT_0036 + IT_0037;
    const ccomplex_t IT_0039 = 576*s_12*IT_0001;
    const ccomplex_t IT_0040 = 1152*IT_0001;
    const ccomplex_t IT_0041 = IT_0039 + IT_0001*IT_0040;
    const ccomplex_t IT_0042 = 288*IT_0015;
    const ccomplex_t IT_0043 = (-288)*IT_0033;
    const ccomplex_t IT_0044 = (-288)*IT_0015;
    const ccomplex_t IT_0045 = (-288)*IT_0006;
    const ccomplex_t IT_0046 = pow(s_12, 2);
    const ccomplex_t IT_0047 = 576*IT_0046;
    const ccomplex_t IT_0048 = IT_0020*IT_0038 + (-0.25)*IT_0023*IT_0041 + (
      -0.25)*IT_0013*(288*IT_0006 + 1152*IT_0024 + IT_0042) + (-0.25)*IT_0032*
      (1152*IT_0020 + 2304*IT_0024 + IT_0044 + IT_0045) + (-0.25)*IT_0025*
      (IT_0001*IT_0009 + IT_0047);
    const ccomplex_t IT_0049 = 1152*IT_0046;
    const ccomplex_t IT_0050 = 2304*s_12;
    const ccomplex_t IT_0051 = (-576)*s_12;
    const ccomplex_t IT_0052 = s_14*(IT_0006*IT_0009 + (576*s_12 + 1152
      *IT_0001)*IT_0025) + IT_0024*(576*IT_0013 + (-1728)*IT_0032 + (-1152)
      *IT_0035) + IT_0020*(1728*IT_0013 + 576*IT_0032 + (-4)*IT_0043) + IT_0008*
      (s_12*IT_0040 + IT_0049 + IT_0001*(2304*IT_0001 + IT_0050)) + s_24*IT_0023
      *((-1152)*IT_0001 + IT_0051);
    const ccomplex_t IT_0053 = s_14*IT_0008*IT_0009 + IT_0010 + (-1728)
      *IT_0013*IT_0024 + ((-1152)*IT_0006 + (-2304)*IT_0015 + IT_0027)*IT_0028 +
       1728*IT_0020*IT_0032 + IT_0025*(IT_0039 + IT_0047) + IT_0015*IT_0049 +
       IT_0001*(IT_0011 + IT_0012 + (-1152)*s_12*IT_0023 + IT_0001*(2304*IT_0015
       + (-1152)*IT_0023) + IT_0025*(IT_0009 + IT_0040) + IT_0015*IT_0050) +
       s_12*(IT_0015*IT_0040 + IT_0023*((-576)*IT_0001 + IT_0051));
    const ccomplex_t IT_0054 = IT_0006*(IT_0010 + IT_0001*(IT_0011 + IT_0012))
       + IT_0013*(IT_0021 + IT_0024*(IT_0026 + IT_0027)) + IT_0028*(IT_0021 +
       IT_0006*IT_0027 + (-2304)*IT_0029 + (-2304)*IT_0030 + IT_0031) + IT_0031
      *IT_0032 + IT_0024*((-144)*IT_0032 + 288*IT_0033)*IT_0034 + IT_0020*
      (IT_0026*IT_0032 + IT_0034*IT_0038) + IT_0025*(IT_0025*IT_0041 + IT_0013*
      (288*IT_0006 + 2304*IT_0020 + IT_0042) + (-4)*IT_0024*(144*IT_0032 +
       IT_0043) + IT_0032*(IT_0044 + IT_0045)) + (-4)*IT_0023*IT_0048 + IT_0008
      *IT_0052 + IT_0015*IT_0053;
    return create_ccomplex_return(IT_0054);
}

