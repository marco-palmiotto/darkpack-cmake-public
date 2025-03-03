#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_1_anti_chi_dm_1_to_V_3_V_3.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_1_anti_chi_dm_1_to_V_3_V_3(
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
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const ccomplex_t IT_0000 = pow(m_V_3, 2);
    const ccomplex_t IT_0001 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*reg_prop, 
      -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0003, 2);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*IT_0004;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = -IT_0006;
    const ccomplex_t IT_0008 = pow(s_14, 2);
    const ccomplex_t IT_0009 = pow(m_V_3, -2);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0012 = s_13*s_23;
    const ccomplex_t IT_0013 = s_12*IT_0000;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = IT_0012 + IT_0014;
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = (IT_0010 + -IT_0011)*(s_12 + 1.5*IT_0011 +
       IT_0016);
    const ccomplex_t IT_0018 = s_12*s_14;
    const ccomplex_t IT_0019 = s_24*IT_0011;
    const ccomplex_t IT_0020 = IT_0018 + -IT_0019;
    const ccomplex_t IT_0021 = pow(s_13, 2);
    const ccomplex_t IT_0022 = s_24*IT_0021;
    const ccomplex_t IT_0023 = s_12*s_13*s_34;
    const ccomplex_t IT_0024 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0025 = s_24*IT_0000*IT_0011;
    const ccomplex_t IT_0026 = -IT_0023 + 0.5*IT_0024 + (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0022 + IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0027;
    const ccomplex_t IT_0029 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*reg_prop, 
      -1);
    const ccomplex_t IT_0030 = IT_0005*IT_0029;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = s_13*s_14*s_23;
    const ccomplex_t IT_0033 = s_23*s_34*IT_0011;
    const ccomplex_t IT_0034 = -IT_0025 + (-2)*IT_0032 + 2*IT_0033;
    const ccomplex_t IT_0035 = IT_0024 + IT_0034;
    const ccomplex_t IT_0036 = IT_0009*IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0006;
    const ccomplex_t IT_0039 = pow(s_24, 2);
    const ccomplex_t IT_0040 = IT_0009*IT_0039;
    const ccomplex_t IT_0041 = (s_12 + 1.5*IT_0011 + IT_0016)*(IT_0011 + 
      -IT_0040);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = s_14*s_24*IT_0009;
    const ccomplex_t IT_0044 = (s_12 + 1.5*IT_0011 + IT_0016)*(s_12 + -IT_0043);
    const ccomplex_t IT_0045 = s_13*s_23*s_24;
    const ccomplex_t IT_0046 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0047 = s_13*s_34*IT_0011;
    const ccomplex_t IT_0048 = s_14*IT_0000*IT_0011;
    const ccomplex_t IT_0049 = (-0.5)*IT_0046 + -IT_0047 + 0.5*IT_0048;
    const ccomplex_t IT_0050 = IT_0045 + IT_0049;
    const ccomplex_t IT_0051 = IT_0009*IT_0050;
    const ccomplex_t IT_0052 = s_12*s_24;
    const ccomplex_t IT_0053 = s_14*IT_0011;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = pow(s_23, 2);
    const ccomplex_t IT_0057 = s_14*IT_0056;
    const ccomplex_t IT_0058 = s_12*s_23*s_34;
    const ccomplex_t IT_0059 = -IT_0048 + 2*IT_0057 + (-2)*IT_0058;
    const ccomplex_t IT_0060 = IT_0046 + IT_0059;
    const ccomplex_t IT_0061 = IT_0009*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = cpow(IT_0038, 2);
    const ccomplex_t IT_0064 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0065 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0066 = IT_0064 + IT_0065;
    const ccomplex_t IT_0067 = pow(m_V_3, 4);
    const ccomplex_t IT_0068 = s_12*IT_0067;
    const ccomplex_t IT_0069 = s_14*s_23*s_34;
    const ccomplex_t IT_0070 = (-0.5)*IT_0068 + (-2)*IT_0069;
    const ccomplex_t IT_0071 = IT_0066 + IT_0070;
    const ccomplex_t IT_0072 = IT_0009*IT_0071;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = s_14*s_24;
    const ccomplex_t IT_0075 = IT_0013 + (-4)*IT_0074;
    const ccomplex_t IT_0076 = pow(m_V_3, -4);
    const ccomplex_t IT_0077 = s_13*s_23*IT_0067;
    const ccomplex_t IT_0078 = pow(m_V_3, 6);
    const ccomplex_t IT_0079 = s_12*IT_0078;
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0077 + IT_0080;
    const ccomplex_t IT_0082 = IT_0076*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = IT_0000*IT_0011;
    const ccomplex_t IT_0085 = 4.5*IT_0084;
    const ccomplex_t IT_0086 = pow(s_34, 2);
    const ccomplex_t IT_0087 = IT_0009*IT_0011*IT_0086;
    const ccomplex_t IT_0088 = s_13*s_24*s_34;
    const ccomplex_t IT_0089 = (-0.5)*IT_0064 + (-0.5)*IT_0065 + 0.25*IT_0068;
    const ccomplex_t IT_0090 = IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = IT_0009*IT_0090;
    const ccomplex_t IT_0092 = (-8)*IT_0007*(IT_0007*IT_0017 + (-2)*(IT_0020 +
       -IT_0028)*IT_0031 + 2*(IT_0020 + IT_0037)*IT_0038) + (-8)*IT_0030*
      (IT_0030*IT_0042 + (-2)*IT_0007*IT_0044 + 2*IT_0031*(IT_0051 + IT_0055) + 
      (-2)*IT_0038*(IT_0055 + -IT_0062)) + (-32)*IT_0063*(IT_0073 + (-0.25)
      *IT_0075 + (-0.25)*IT_0083 + (-0.25)*IT_0085) + (-2)*IT_0031*(IT_0038*((
      -16)*IT_0013 + 8*IT_0082 + (-4)*IT_0084 + 16*IT_0087) + 16*IT_0031*((-0.25
      )*IT_0075 + (-0.25)*IT_0083 + (-0.25)*IT_0085 + IT_0091));
    return create_ccomplex_return(IT_0092);
}

