#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_1_chi_dm_1_to_anti_V_3_V_3(
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
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const ccomplex_t IT_0000 = pow(m_V, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_23 + IT_0000 + reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = (-2)*IT_0006;
    const ccomplex_t IT_0008 = pow(m_V, -2);
    const ccomplex_t IT_0009 = pow(s_14, 2);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0012 = s_12*IT_0000;
    const ccomplex_t IT_0013 = s_13*s_23;
    const ccomplex_t IT_0014 = (-2)*IT_0013;
    const ccomplex_t IT_0015 = IT_0012 + IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = (IT_0010 + -IT_0011)*(s_12 + 1.5*IT_0011 + (
      -0.5)*IT_0016);
    const ccomplex_t IT_0018 = s_14*s_24*IT_0008;
    const ccomplex_t IT_0019 = (s_12 + 1.5*IT_0011 + (-0.5)*IT_0016)*(s_12 + 
      -IT_0018);
    const ccomplex_t IT_0020 = cpow((-2)*s_13 + IT_0000 + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0005*IT_0020;
    const ccomplex_t IT_0022 = 2*IT_0021;
    const ccomplex_t IT_0023 = s_12*s_14;
    const ccomplex_t IT_0024 = s_24*IT_0011;
    const ccomplex_t IT_0025 = IT_0023 + -IT_0024;
    const ccomplex_t IT_0026 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0027 = pow(s_13, 2);
    const ccomplex_t IT_0028 = s_24*IT_0027;
    const ccomplex_t IT_0029 = s_12*s_13*s_34;
    const ccomplex_t IT_0030 = s_24*IT_0000*IT_0011;
    const ccomplex_t IT_0031 = 2*IT_0028 + (-2)*IT_0029 + -IT_0030;
    const ccomplex_t IT_0032 = IT_0026 + IT_0031;
    const ccomplex_t IT_0033 = IT_0008*IT_0032;
    const ccomplex_t IT_0034 = -IT_0021;
    const ccomplex_t IT_0035 = pow(s_24, 2);
    const ccomplex_t IT_0036 = IT_0008*IT_0035;
    const ccomplex_t IT_0037 = (s_12 + 1.5*IT_0011 + (-0.5)*IT_0016)*(IT_0011 
      + -IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0040 = s_13*s_23*s_24;
    const ccomplex_t IT_0041 = s_14*IT_0000*IT_0011;
    const ccomplex_t IT_0042 = s_13*s_34*IT_0011;
    const ccomplex_t IT_0043 = (-2)*IT_0040 + -IT_0041 + 2*IT_0042;
    const ccomplex_t IT_0044 = IT_0039 + IT_0043;
    const ccomplex_t IT_0045 = IT_0008*IT_0044;
    const ccomplex_t IT_0046 = s_12*s_24;
    const ccomplex_t IT_0047 = s_14*IT_0011;
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0046 + IT_0048;
    const ccomplex_t IT_0050 = cpow(IT_0034, 2);
    const ccomplex_t IT_0051 = pow(m_V, 4);
    const ccomplex_t IT_0052 = s_12*IT_0051;
    const ccomplex_t IT_0053 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0054 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0055 = s_13*s_24*s_34;
    const ccomplex_t IT_0056 = (-2)*IT_0053 + (-2)*IT_0054 + 4*IT_0055;
    const ccomplex_t IT_0057 = IT_0052 + IT_0056;
    const ccomplex_t IT_0058 = IT_0008*IT_0057;
    const ccomplex_t IT_0059 = s_14*s_24;
    const ccomplex_t IT_0060 = IT_0012 + (-4)*IT_0059;
    const ccomplex_t IT_0061 = pow(m_V, -4);
    const ccomplex_t IT_0062 = pow(m_V, 6);
    const ccomplex_t IT_0063 = s_12*IT_0062;
    const ccomplex_t IT_0064 = s_13*s_23*IT_0051;
    const ccomplex_t IT_0065 = (-2)*IT_0064;
    const ccomplex_t IT_0066 = IT_0063 + IT_0065;
    const ccomplex_t IT_0067 = IT_0061*IT_0066;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0000*IT_0011;
    const ccomplex_t IT_0070 = 4.5*IT_0069;
    const ccomplex_t IT_0071 = s_13*s_14*s_23;
    const ccomplex_t IT_0072 = s_23*s_34*IT_0011;
    const ccomplex_t IT_0073 = -IT_0030 + (-2)*IT_0071 + 2*IT_0072;
    const ccomplex_t IT_0074 = IT_0026 + IT_0073;
    const ccomplex_t IT_0075 = IT_0008*IT_0074;
    const ccomplex_t IT_0076 = pow(s_34, 2);
    const ccomplex_t IT_0077 = IT_0008*IT_0011*IT_0076;
    const ccomplex_t IT_0078 = pow(s_23, 2);
    const ccomplex_t IT_0079 = s_14*IT_0078;
    const ccomplex_t IT_0080 = s_12*s_23*s_34;
    const ccomplex_t IT_0081 = -IT_0041 + 2*IT_0079 + (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0039 + IT_0081;
    const ccomplex_t IT_0083 = IT_0008*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = s_14*s_23*s_34;
    const ccomplex_t IT_0086 = (-2)*IT_0053 + (-2)*IT_0054 + 4*IT_0085;
    const ccomplex_t IT_0087 = IT_0052 + IT_0086;
    const ccomplex_t IT_0088 = IT_0008*IT_0087;
    const ccomplex_t IT_0089 = (-16)*IT_0007*(IT_0025 + (-0.5)*IT_0075) + (-2)
      *IT_0034*((-16)*IT_0012 + (-4)*IT_0067 + (-4)*IT_0069 + 16*IT_0077) + 16
      *IT_0022*(IT_0049 + -IT_0084) + 8*IT_0006*(IT_0060 + IT_0068 + IT_0070 + 
      -IT_0088);
    const ccomplex_t IT_0090 = (-8)*IT_0007*(IT_0007*IT_0017 + (-2)*IT_0019
      *IT_0022 + (-2)*(IT_0025 + (-0.5)*IT_0033)*IT_0034) + (-8)*IT_0022*
      (IT_0022*IT_0038 + -IT_0034*(IT_0045 + (-2)*IT_0049)) + (-8)*IT_0050*
      (IT_0058 + -IT_0060 + -IT_0068 + -IT_0070) + IT_0006*IT_0089;
    return create_ccomplex_return(IT_0090);
}

