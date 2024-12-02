#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_V_3_V_3.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_V_3_V_3(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0001 = (-0.5)*IT_0000;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0002;
    const ccomplex_t IT_0004 = pow(m_V, 2);
    const ccomplex_t IT_0005 = cpow((-2)*s_23 + IT_0004 + reg_prop, -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = (-2)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_V, -2);
    const ccomplex_t IT_0009 = pow(s_14, 2);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0012 = s_12*IT_0004;
    const ccomplex_t IT_0013 = s_13*s_23;
    const ccomplex_t IT_0014 = (-2)*IT_0013;
    const ccomplex_t IT_0015 = IT_0012 + IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = (IT_0010 + -IT_0011)*(s_12 + 1.5*IT_0011 + (
      -0.5)*IT_0016);
    const ccomplex_t IT_0018 = s_14*s_24*IT_0008;
    const ccomplex_t IT_0019 = (s_12 + 1.5*IT_0011 + (-0.5)*IT_0016)*(s_12 + 
      -IT_0018);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0004 + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0003*IT_0020;
    const ccomplex_t IT_0022 = 2*IT_0021;
    const ccomplex_t IT_0023 = s_12*s_14*IT_0004;
    const ccomplex_t IT_0024 = s_13*s_14*s_23;
    const ccomplex_t IT_0025 = s_24*IT_0004*IT_0011;
    const ccomplex_t IT_0026 = s_23*s_34*IT_0011;
    const ccomplex_t IT_0027 = (-2)*IT_0024 + -IT_0025 + 2*IT_0026;
    const ccomplex_t IT_0028 = IT_0023 + IT_0027;
    const ccomplex_t IT_0029 = IT_0008*IT_0028;
    const ccomplex_t IT_0030 = s_12*s_14;
    const ccomplex_t IT_0031 = s_24*IT_0011;
    const ccomplex_t IT_0032 = IT_0030 + -IT_0031;
    const ccomplex_t IT_0033 = -IT_0021;
    const ccomplex_t IT_0034 = pow(s_13, 2);
    const ccomplex_t IT_0035 = s_24*IT_0034;
    const ccomplex_t IT_0036 = s_12*s_13*s_34;
    const ccomplex_t IT_0037 = -IT_0025 + 2*IT_0035 + (-2)*IT_0036;
    const ccomplex_t IT_0038 = IT_0023 + IT_0037;
    const ccomplex_t IT_0039 = IT_0008*IT_0038;
    const ccomplex_t IT_0040 = pow(s_24, 2);
    const ccomplex_t IT_0041 = IT_0008*IT_0040;
    const ccomplex_t IT_0042 = (s_12 + 1.5*IT_0011 + (-0.5)*IT_0016)*(IT_0011 
      + -IT_0041);
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = s_12*s_24*IT_0004;
    const ccomplex_t IT_0045 = s_14*IT_0004*IT_0011;
    const ccomplex_t IT_0046 = s_13*s_23*s_24;
    const ccomplex_t IT_0047 = s_13*s_34*IT_0011;
    const ccomplex_t IT_0048 = -IT_0045 + (-2)*IT_0046 + 2*IT_0047;
    const ccomplex_t IT_0049 = IT_0044 + IT_0048;
    const ccomplex_t IT_0050 = IT_0008*IT_0049;
    const ccomplex_t IT_0051 = s_12*s_24;
    const ccomplex_t IT_0052 = s_14*IT_0011;
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0051 + IT_0053;
    const ccomplex_t IT_0055 = pow(s_23, 2);
    const ccomplex_t IT_0056 = s_14*IT_0055;
    const ccomplex_t IT_0057 = s_12*s_23*s_34;
    const ccomplex_t IT_0058 = 0.5*IT_0044 + (-0.5)*IT_0045 + -IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = IT_0008*IT_0059;
    const ccomplex_t IT_0061 = IT_0004*IT_0011;
    const ccomplex_t IT_0062 = pow(s_34, 2);
    const ccomplex_t IT_0063 = IT_0008*IT_0011*IT_0062;
    const ccomplex_t IT_0064 = pow(m_V, -4);
    const ccomplex_t IT_0065 = pow(m_V, 6);
    const ccomplex_t IT_0066 = s_12*IT_0065;
    const ccomplex_t IT_0067 = pow(m_V, 4);
    const ccomplex_t IT_0068 = s_13*s_23*IT_0067;
    const ccomplex_t IT_0069 = (-2)*IT_0068;
    const ccomplex_t IT_0070 = IT_0066 + IT_0069;
    const ccomplex_t IT_0071 = IT_0064*IT_0070;
    const ccomplex_t IT_0072 = s_12*IT_0067;
    const ccomplex_t IT_0073 = s_13*s_23*IT_0004;
    const ccomplex_t IT_0074 = s_14*s_24*IT_0004;
    const ccomplex_t IT_0075 = s_14*s_23*s_34;
    const ccomplex_t IT_0076 = (-2)*IT_0073 + (-2)*IT_0074 + 4*IT_0075;
    const ccomplex_t IT_0077 = IT_0072 + IT_0076;
    const ccomplex_t IT_0078 = IT_0008*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0071;
    const ccomplex_t IT_0080 = 4.5*IT_0061;
    const ccomplex_t IT_0081 = s_14*s_24;
    const ccomplex_t IT_0082 = (-4)*IT_0081;
    const ccomplex_t IT_0083 = cpow(IT_0033, 2);
    const ccomplex_t IT_0084 = s_13*s_24*s_34;
    const ccomplex_t IT_0085 = (-2)*IT_0073 + (-2)*IT_0074 + 4*IT_0084;
    const ccomplex_t IT_0086 = IT_0072 + IT_0085;
    const ccomplex_t IT_0087 = IT_0008*IT_0086;
    const ccomplex_t IT_0088 = (-8)*IT_0007*(IT_0007*IT_0017 + (-2)*IT_0019
      *IT_0022 + -IT_0006*(IT_0029 + (-2)*IT_0032) + (-2)*IT_0033*(IT_0032 + -1.
      /2*IT_0039)) + (-8)*IT_0022*(IT_0022*IT_0043 + -IT_0033*(IT_0050 + (-2)
      *IT_0054) + (-2)*IT_0006*(IT_0054 + -IT_0060)) + (-2)*IT_0006*(IT_0033*((
      -16)*IT_0012 + (-4)*IT_0061 + 16*IT_0063 + (-4)*IT_0071) + (-4)*IT_0006*
      (IT_0012 + -IT_0078 + IT_0079 + IT_0080 + IT_0082)) + 8*IT_0083*(IT_0012 +
       IT_0079 + IT_0080 + IT_0082 + -IT_0087);
    return create_ccomplex_return(IT_0088);
}

