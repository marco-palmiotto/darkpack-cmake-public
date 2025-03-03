#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_2_chi_dm_2_to_anti_V_3_V_3(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_V_3 = param->m_V_3;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0001 = (-0.5)*IT_0000;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0002;
    const ccomplex_t IT_0004 = pow(m_V_3, 2);
    const ccomplex_t IT_0005 = cpow(s_13 + (-0.5)*IT_0004 + (-0.5)*reg_prop, 
      -1);
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = s_12*IT_0004;
    const ccomplex_t IT_0009 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0010 = IT_0004*IT_0009;
    const ccomplex_t IT_0011 = pow(s_34, 2);
    const ccomplex_t IT_0012 = pow(m_V_3, -2);
    const ccomplex_t IT_0013 = IT_0009*IT_0011*IT_0012;
    const ccomplex_t IT_0014 = pow(m_V_3, -4);
    const ccomplex_t IT_0015 = pow(m_V_3, 4);
    const ccomplex_t IT_0016 = s_13*s_23*IT_0015;
    const ccomplex_t IT_0017 = pow(m_V_3, 6);
    const ccomplex_t IT_0018 = s_12*IT_0017;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = IT_0014*IT_0020;
    const ccomplex_t IT_0022 = cpow(s_23 + (-0.5)*IT_0004 + (-0.5)*reg_prop, 
      -1);
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = s_14*s_24;
    const ccomplex_t IT_0026 = s_13*s_24*s_34;
    const ccomplex_t IT_0027 = s_13*s_23*IT_0004;
    const ccomplex_t IT_0028 = s_14*s_24*IT_0004;
    const ccomplex_t IT_0029 = s_12*IT_0015;
    const ccomplex_t IT_0030 = (-0.5)*IT_0027 + (-0.5)*IT_0028 + 0.25*IT_0029;
    const ccomplex_t IT_0031 = IT_0026 + IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0031;
    const ccomplex_t IT_0033 = -IT_0021;
    const ccomplex_t IT_0034 = 4.5*IT_0010;
    const ccomplex_t IT_0035 = cpow(IT_0024, 2);
    const ccomplex_t IT_0036 = IT_0027 + IT_0028;
    const ccomplex_t IT_0037 = s_14*s_23*s_34;
    const ccomplex_t IT_0038 = (-0.5)*IT_0029 + (-2)*IT_0037;
    const ccomplex_t IT_0039 = IT_0036 + IT_0038;
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = pow(s_14, 2);
    const ccomplex_t IT_0043 = IT_0012*IT_0042;
    const ccomplex_t IT_0044 = s_13*s_23;
    const ccomplex_t IT_0045 = (-0.5)*IT_0008;
    const ccomplex_t IT_0046 = IT_0044 + IT_0045;
    const ccomplex_t IT_0047 = IT_0012*IT_0046;
    const ccomplex_t IT_0048 = -(IT_0009 + -IT_0043)*(s_12 + 1.5*IT_0009 +
       IT_0047);
    const ccomplex_t IT_0049 = s_14*s_24*IT_0012;
    const ccomplex_t IT_0050 = (s_12 + 1.5*IT_0009 + IT_0047)*(s_12 + -IT_0049);
    const ccomplex_t IT_0051 = -IT_0006;
    const ccomplex_t IT_0052 = s_12*s_14*IT_0004;
    const ccomplex_t IT_0053 = s_24*IT_0004*IT_0009;
    const ccomplex_t IT_0054 = s_13*s_14*s_23;
    const ccomplex_t IT_0055 = s_23*s_34*IT_0009;
    const ccomplex_t IT_0056 = -IT_0053 + (-2)*IT_0054 + 2*IT_0055;
    const ccomplex_t IT_0057 = IT_0052 + IT_0056;
    const ccomplex_t IT_0058 = IT_0012*IT_0057;
    const ccomplex_t IT_0059 = (-0.5)*IT_0058;
    const ccomplex_t IT_0060 = s_12*s_14;
    const ccomplex_t IT_0061 = s_24*IT_0009;
    const ccomplex_t IT_0062 = IT_0060 + -IT_0061;
    const ccomplex_t IT_0063 = pow(s_13, 2);
    const ccomplex_t IT_0064 = s_24*IT_0063;
    const ccomplex_t IT_0065 = s_12*s_13*s_34;
    const ccomplex_t IT_0066 = 0.5*IT_0052 + (-0.5)*IT_0053 + -IT_0065;
    const ccomplex_t IT_0067 = IT_0064 + IT_0066;
    const ccomplex_t IT_0068 = IT_0012*IT_0067;
    const ccomplex_t IT_0069 = 2*IT_0007;
    const ccomplex_t IT_0070 = pow(s_24, 2);
    const ccomplex_t IT_0071 = IT_0012*IT_0070;
    const ccomplex_t IT_0072 = (s_12 + 1.5*IT_0009 + IT_0047)*(IT_0009 + 
      -IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = s_12*s_24*IT_0004;
    const ccomplex_t IT_0075 = s_14*IT_0004*IT_0009;
    const ccomplex_t IT_0076 = pow(s_23, 2);
    const ccomplex_t IT_0077 = s_14*IT_0076;
    const ccomplex_t IT_0078 = s_12*s_23*s_34;
    const ccomplex_t IT_0079 = -IT_0075 + 2*IT_0077 + (-2)*IT_0078;
    const ccomplex_t IT_0080 = IT_0074 + IT_0079;
    const ccomplex_t IT_0081 = IT_0012*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = s_12*s_24;
    const ccomplex_t IT_0084 = s_14*IT_0009;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = s_13*s_23*s_24;
    const ccomplex_t IT_0088 = s_13*s_34*IT_0009;
    const ccomplex_t IT_0089 = (-0.5)*IT_0074 + 0.5*IT_0075 + -IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = IT_0012*IT_0090;
    const ccomplex_t IT_0092 = (-2)*IT_0007*(((-16)*IT_0008 + (-4)*IT_0010 +
       16*IT_0013 + 8*IT_0021)*IT_0024 + (-4)*IT_0007*(IT_0008 + (-4)*IT_0025 + 
      (-4)*IT_0032 + IT_0033 + IT_0034)) + 8*IT_0035*(IT_0008 + (-4)*IT_0025 +
       IT_0033 + IT_0034 + (-4)*IT_0041) + (-8)*IT_0023*(IT_0023*IT_0048 + (-2)
      *IT_0050*IT_0051 + 2*IT_0024*(IT_0059 + IT_0062) + -(IT_0062 + -IT_0068)
      *IT_0069) + (-8)*IT_0051*(IT_0051*IT_0073 + 2*IT_0024*(IT_0082 + -IT_0086)
       + IT_0069*(IT_0086 + IT_0091));
    return create_ccomplex_return(IT_0092);
}

