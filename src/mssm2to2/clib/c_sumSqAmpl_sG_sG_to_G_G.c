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
    const ccomplex_t IT_0000 = pow(s_13 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = -g_s;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*IT_0002;
    const ccomplex_t IT_0004 = IT_0000*IT_0003;
    const ccomplex_t IT_0005 = (-0.5)*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = pow(m_sG, 2);
    const ccomplex_t IT_0008 = 576*s_12*IT_0007;
    const ccomplex_t IT_0009 = pow(s_23 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = pow(s_12, 2);
    const ccomplex_t IT_0013 = 576*IT_0012;
    const ccomplex_t IT_0014 = cpow(s_12 + IT_0007 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0015 = IT_0001*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = (g_s + -IT_0001)*IT_0016;
    const ccomplex_t IT_0018 = 0.25*IT_0017;
    const ccomplex_t IT_0019 = 1152*IT_0012*IT_0018;
    const ccomplex_t IT_0020 = g_s*IT_0016;
    const ccomplex_t IT_0021 = (-0.25)*IT_0020;
    const ccomplex_t IT_0022 = IT_0001*IT_0016;
    const ccomplex_t IT_0023 = 0.25*IT_0022;
    const ccomplex_t IT_0024 = 4608*s_12;
    const ccomplex_t IT_0025 = s_24*IT_0023*IT_0024;
    const ccomplex_t IT_0026 = s_14*s_24;
    const ccomplex_t IT_0027 = (-0.25)*IT_0010;
    const ccomplex_t IT_0028 = cpow(IT_0027, 2);
    const ccomplex_t IT_0029 = 0.25*IT_0004;
    const ccomplex_t IT_0030 = cpow(IT_0029, 2);
    const ccomplex_t IT_0031 = (-0.25)*IT_0004;
    const ccomplex_t IT_0032 = 2304*IT_0018;
    const ccomplex_t IT_0033 = IT_0027*IT_0032;
    const ccomplex_t IT_0034 = (-4608)*IT_0023;
    const ccomplex_t IT_0035 = 0.25*IT_0010;
    const ccomplex_t IT_0036 = (-2304)*IT_0018;
    const ccomplex_t IT_0037 = IT_0029*IT_0036;
    const ccomplex_t IT_0038 = (-2304)*IT_0028 + (-2304)*IT_0030 + IT_0031*
      (4608*IT_0029 + (-2304)*IT_0031 + IT_0032) + IT_0033 + IT_0018*((-9216)
      *IT_0018 + (-4608)*IT_0021 + IT_0034) + IT_0021*IT_0034 + IT_0035*(4608
      *IT_0027 + (-2304)*IT_0035 + IT_0036) + IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0004;
    const ccomplex_t IT_0040 = 1152*IT_0007;
    const ccomplex_t IT_0041 = 4608*IT_0012;
    const ccomplex_t IT_0042 = (-0.5)*IT_0010;
    const ccomplex_t IT_0043 = IT_0005 + IT_0042;
    const ccomplex_t IT_0044 = (-1152)*s_12*IT_0007;
    const ccomplex_t IT_0045 = (-1152)*IT_0012;
    const ccomplex_t IT_0046 = (-576)*IT_0012;
    const ccomplex_t IT_0047 = s_12*s_14;
    const ccomplex_t IT_0048 = (-288)*IT_0047;
    const ccomplex_t IT_0049 = s_24*IT_0007;
    const ccomplex_t IT_0050 = 144*IT_0049;
    const ccomplex_t IT_0051 = IT_0048 + IT_0050;
    const ccomplex_t IT_0052 = (-144)*IT_0049;
    const ccomplex_t IT_0053 = 288*IT_0047 + IT_0052;
    const ccomplex_t IT_0054 = s_14*IT_0007;
    const ccomplex_t IT_0055 = 576*IT_0021;
    const ccomplex_t IT_0056 = 576*IT_0018;
    const ccomplex_t IT_0057 = IT_0055 + IT_0056;
    const ccomplex_t IT_0058 = (-576)*IT_0021;
    const ccomplex_t IT_0059 = (-576)*IT_0018;
    const ccomplex_t IT_0060 = IT_0058 + IT_0059;
    const ccomplex_t IT_0061 = 1152*IT_0011;
    const ccomplex_t IT_0062 = 1152*IT_0005;
    const ccomplex_t IT_0063 = IT_0061 + IT_0062;
    const ccomplex_t IT_0064 = (-3456)*IT_0018;
    const ccomplex_t IT_0065 = (-2304)*IT_0023;
    const ccomplex_t IT_0066 = (-1152)*IT_0042;
    const ccomplex_t IT_0067 = (-1152)*IT_0039;
    const ccomplex_t IT_0068 = 1152*IT_0042;
    const ccomplex_t IT_0069 = 1152*IT_0039;
    const ccomplex_t IT_0070 = IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = 2304*IT_0023;
    const ccomplex_t IT_0072 = 3456*IT_0018;
    const ccomplex_t IT_0073 = (-1152)*IT_0011;
    const ccomplex_t IT_0074 = (-1152)*IT_0005;
    const ccomplex_t IT_0075 = IT_0033 + IT_0023*(3456*IT_0027 + (-3456)
      *IT_0035) + IT_0035*IT_0036 + IT_0011*(2304*IT_0027 + (-2304)*IT_0035 +
       IT_0055 + IT_0056) + IT_0005*IT_0057 + IT_0042*((-576)*IT_0018 + (-2304)
      *IT_0027 + 2304*IT_0035 + IT_0058) + IT_0039*IT_0060 + IT_0029*(IT_0063 +
       IT_0064 + IT_0065 + IT_0066 + IT_0067) + IT_0031*(IT_0070 + IT_0071 +
       IT_0072 + IT_0073 + IT_0074);
    const ccomplex_t IT_0076 = s_12*s_24;
    const ccomplex_t IT_0077 = 288*IT_0076;
    const ccomplex_t IT_0078 = (-144)*IT_0054;
    const ccomplex_t IT_0079 = IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = (-8)*IT_0018;
    const ccomplex_t IT_0081 = (-4)*IT_0005;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = (-8)*IT_0023;
    const ccomplex_t IT_0084 = (-288)*IT_0076;
    const ccomplex_t IT_0085 = 144*IT_0054;
    const ccomplex_t IT_0086 = IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = (-4)*IT_0039;
    const ccomplex_t IT_0088 = (-4)*IT_0011;
    const ccomplex_t IT_0089 = 2304*IT_0039;
    const ccomplex_t IT_0090 = IT_0037 + IT_0029*(2304*IT_0005 + (-3456)
      *IT_0023 + (-2304)*IT_0039) + IT_0042*(576*IT_0018 + IT_0055) + IT_0039
      *IT_0057 + (IT_0005 + IT_0011)*IT_0060 + IT_0027*(IT_0063 + IT_0066 +
       IT_0067 + IT_0071 + IT_0072) + IT_0035*(IT_0064 + IT_0065 + IT_0070 +
       IT_0073 + IT_0074) + IT_0031*((-2304)*IT_0005 + 3456*IT_0023 + IT_0032 +
       IT_0089);
    const ccomplex_t IT_0091 = 9216*s_14*IT_0023;
    const ccomplex_t IT_0092 = 9216*s_24*IT_0023;
    const ccomplex_t IT_0093 = 1152*s_12;
    const ccomplex_t IT_0094 = 2304*s_12*IT_0018;
    const ccomplex_t IT_0095 = (-2304)*IT_0007;
    const ccomplex_t IT_0096 = (-1152)*s_12;
    const ccomplex_t IT_0097 = 9216*s_12;
    const ccomplex_t IT_0098 = (-2304)*s_12;
    const ccomplex_t IT_0099 = IT_0005*IT_0096;
    const ccomplex_t IT_0100 = IT_0006*IT_0040 + IT_0021*(IT_0091 + IT_0092) +
       IT_0011*(IT_0011*IT_0040 + IT_0005*IT_0093 + IT_0094 + IT_0042*IT_0095) +
       IT_0039*(IT_0039*IT_0040 + IT_0042*IT_0093 + IT_0094 + IT_0005*IT_0095 +
       IT_0011*IT_0096) + IT_0023*(s_14*(2304*IT_0011 + (-2304)*IT_0042) + s_24*
      ((-2304)*IT_0005 + IT_0089) + IT_0023*(9216*IT_0007 + IT_0097)) + IT_0018*
      (IT_0007*(2304*IT_0011 + 9216*IT_0018 + IT_0089) + IT_0091 + IT_0092 +
       IT_0018*IT_0097 + IT_0043*(IT_0095 + IT_0098)) + IT_0042*(IT_0040*IT_0042
       + IT_0099);
    const ccomplex_t IT_0101 = IT_0006*IT_0008 + IT_0011*(IT_0008*IT_0011 +
       IT_0005*IT_0013 + IT_0019) + IT_0021*IT_0025 + IT_0026*IT_0038 + IT_0018*
      (s_14*IT_0023*IT_0024 + IT_0025 + s_12*(4608*IT_0007*IT_0018 + (IT_0011 +
       IT_0039)*IT_0040) + IT_0018*IT_0041 + IT_0043*(IT_0044 + IT_0045)) +
       IT_0039*(IT_0019 + IT_0008*IT_0039 + IT_0005*IT_0044 + IT_0011*IT_0046) +
       (-4)*IT_0042*((-0.25)*IT_0013*IT_0039 + (-0.25)*IT_0008*IT_0042 + (-0.25)
      *IT_0011*IT_0044 + (-0.25)*IT_0005*IT_0046 + IT_0031*IT_0051 + IT_0029
      *IT_0053) + IT_0054*IT_0075 + IT_0027*(IT_0051*IT_0083 + IT_0082*IT_0086 +
       IT_0079*IT_0087) + IT_0035*(IT_0079*IT_0082 + IT_0053*IT_0083 + IT_0086
      *IT_0087) + IT_0029*(IT_0053*IT_0080 + IT_0079*IT_0083 + IT_0051*IT_0088) 
      + IT_0031*(IT_0051*IT_0080 + IT_0083*IT_0086 + IT_0053*IT_0088) + IT_0049
      *IT_0090 + IT_0023*(IT_0023*(4608*s_12*IT_0007 + IT_0041) + s_14*(IT_0021
      *IT_0024 + IT_0011*IT_0093 + IT_0042*IT_0096) + s_24*(IT_0039*IT_0093 +
       IT_0099)) + IT_0007*IT_0100;
    return create_ccomplex_return(IT_0101);
}

