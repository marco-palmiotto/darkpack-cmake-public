#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_anti_su_L_to_A_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_anti_su_L_to_A_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_ul = param->Gamma_ul;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = pow(e_em, 2);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0001, -2);
    const ccomplex_t IT_0005 = cpow(IT_0003, -2);
    const ccomplex_t IT_0006 = tan(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, 2);
    const ccomplex_t IT_0008 = cpow(1 + IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0000*IT_0002*IT_0003*
      (IT_0004 + (-3)*IT_0005)*IT_0008;
    const ccomplex_t IT_0010 = (-0.222222222222222)*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0010, 2);
    const ccomplex_t IT_0012 = pow(m_su_L, 2);
    const ccomplex_t IT_0013 = pow(m_Z, -2);
    const ccomplex_t IT_0014 = pow(s_14, 2);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0012 + IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0019 = (-0.666666666666667)*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0003, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*e_em*(IT_0002*IT_0003 + (-3)
      *IT_0001*IT_0020);
    const ccomplex_t IT_0022 = 0.166666666666667*IT_0021;
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_su_L
      *Gamma_ul + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_su_L
      *Gamma_ul + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = (-1.33333333333333)*IT_0018;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = -IT_0023;
    const ccomplex_t IT_0032 = IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = IT_0028*IT_0032;
    const ccomplex_t IT_0034 = (-0.5)*IT_0026 + 0.5*IT_0033;
    const ccomplex_t IT_0035 = IT_0028*IT_0030;
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = 3*s_13;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = IT_0025*IT_0030;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = 3*s_23;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = IT_0038 + IT_0042;
    const ccomplex_t IT_0044 = conj(IT_0036)*IT_0037;
    const ccomplex_t IT_0045 = conj(IT_0040)*IT_0041;
    const ccomplex_t IT_0046 = IT_0044 + IT_0045;
    const ccomplex_t IT_0047 = 3*s_12;
    const ccomplex_t IT_0048 = conj(IT_0036)*IT_0047;
    const ccomplex_t IT_0049 = IT_0036*IT_0047;
    const ccomplex_t IT_0050 = 3*IT_0012;
    const ccomplex_t IT_0051 = IT_0040*IT_0050;
    const ccomplex_t IT_0052 = IT_0036*IT_0050;
    const ccomplex_t IT_0053 = 3*IT_0010;
    const ccomplex_t IT_0054 = (-3)*IT_0010;
    const ccomplex_t IT_0055 = s_14*s_24*IT_0013;
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = s_12 + IT_0056;
    const ccomplex_t IT_0058 = 0.333333333333333*IT_0021;
    const ccomplex_t IT_0059 = IT_0029*IT_0058;
    const ccomplex_t IT_0060 = -IT_0025*(IT_0030 + -IT_0059);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = conj(IT_0034)*IT_0041;
    const ccomplex_t IT_0064 = conj(IT_0040)*IT_0050;
    const ccomplex_t IT_0065 = IT_0048 + IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = -IT_0028*(IT_0030 + -IT_0059);
    const ccomplex_t IT_0067 = -IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0040*IT_0047;
    const ccomplex_t IT_0070 = IT_0034*IT_0037;
    const ccomplex_t IT_0071 = IT_0052 + IT_0053 + IT_0069 + IT_0070;
    const ccomplex_t IT_0072 = conj(IT_0040)*IT_0047;
    const ccomplex_t IT_0073 = conj(IT_0034)*IT_0037;
    const ccomplex_t IT_0074 = conj(IT_0036)*IT_0050;
    const ccomplex_t IT_0075 = IT_0054 + IT_0072 + IT_0073 + IT_0074;
    const ccomplex_t IT_0076 = IT_0019*IT_0058;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0023 + IT_0077;
    const ccomplex_t IT_0079 = IT_0025*IT_0078;
    const ccomplex_t IT_0080 = IT_0030 + IT_0076;
    const ccomplex_t IT_0081 = -IT_0023 + -IT_0059;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = IT_0028*IT_0082;
    const ccomplex_t IT_0084 = (-0.5)*IT_0079 + 0.5*IT_0083;
    const ccomplex_t IT_0085 = IT_0034*IT_0041;
    const ccomplex_t IT_0086 = s_14*s_34*IT_0013;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = s_13 + IT_0087;
    const ccomplex_t IT_0089 = IT_0025*(IT_0030 + -IT_0059);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0035;
    const ccomplex_t IT_0092 = (-0.5)*IT_0033 + 0.5*IT_0079;
    const ccomplex_t IT_0093 = s_24*s_34*IT_0013;
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = s_23 + IT_0094;
    const ccomplex_t IT_0096 = IT_0037*IT_0068;
    const ccomplex_t IT_0097 = IT_0041*IT_0062;
    const ccomplex_t IT_0098 = IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = IT_0037*conj(IT_0068);
    const ccomplex_t IT_0100 = IT_0041*conj(IT_0062);
    const ccomplex_t IT_0101 = IT_0099 + IT_0100;
    const ccomplex_t IT_0102 = IT_0050*IT_0068;
    const ccomplex_t IT_0103 = IT_0047*IT_0068;
    const ccomplex_t IT_0104 = IT_0050*IT_0062;
    const ccomplex_t IT_0105 = IT_0047*conj(IT_0068);
    const ccomplex_t IT_0106 = IT_0054*IT_0084 + IT_0053*conj(IT_0084) + 
      (IT_0047*conj(IT_0062) + IT_0050*conj(IT_0068) + IT_0037*conj(IT_0084))
      *IT_0091 + conj(IT_0092)*IT_0098 + IT_0092*IT_0101 + conj(IT_0091)*
      (IT_0047*IT_0062 + IT_0037*IT_0084 + IT_0102) + conj(IT_0090)*(IT_0053 +
       IT_0041*IT_0084 + IT_0103 + IT_0104) + IT_0090*(IT_0054 + IT_0050*conj
      (IT_0062) + IT_0041*conj(IT_0084) + IT_0105);
    const ccomplex_t IT_0107 = pow(s_24, 2);
    const ccomplex_t IT_0108 = IT_0013*IT_0107;
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = IT_0012 + IT_0109;
    const ccomplex_t IT_0111 = pow(s_34, 2);
    const ccomplex_t IT_0112 = IT_0013*IT_0111;
    const ccomplex_t IT_0113 = (-3)*s_13;
    const ccomplex_t IT_0114 = (-3)*s_23;
    const ccomplex_t IT_0115 = (-3)*IT_0012;
    const ccomplex_t IT_0116 = (-0.333333333333333)*IT_0092;
    const ccomplex_t IT_0117 = IT_0053*IT_0092 + conj(IT_0092)*(IT_0054 +
       IT_0091*IT_0113 + IT_0090*IT_0114) + (-3)*conj(IT_0091)*(s_12*IT_0090 + (
      -0.333333333333333)*IT_0091*IT_0115 + IT_0113*IT_0116) + (-3)*conj(IT_0090
      )*(s_12*IT_0091 + (-0.333333333333333)*IT_0090*IT_0115 + IT_0114*IT_0116);
    const ccomplex_t IT_0118 = (-9)*IT_0011 + IT_0017*(conj(IT_0034)*IT_0043 +
       IT_0034*IT_0046 + IT_0040*IT_0048 + conj(IT_0040)*(IT_0049 + IT_0051) +
       conj(IT_0036)*(IT_0052 + IT_0053) + IT_0036*IT_0054) + IT_0057*(conj
      (IT_0040)*IT_0053 + IT_0040*IT_0054 + IT_0062*IT_0065 + conj(IT_0068)
      *IT_0071 + IT_0068*IT_0075 + IT_0046*IT_0084 + IT_0043*conj(IT_0084) +
       conj(IT_0062)*(IT_0049 + IT_0051 + IT_0085)) + IT_0088*(conj(IT_0034)
      *IT_0053 + IT_0034*IT_0054 + IT_0065*IT_0090 + (IT_0049 + IT_0051 +
       IT_0085)*conj(IT_0090) + IT_0075*IT_0091 + IT_0071*conj(IT_0091) +
       IT_0046*IT_0092 + IT_0043*conj(IT_0092)) + IT_0095*IT_0106 + (conj
      (IT_0084)*IT_0098 + IT_0084*IT_0101 + conj(IT_0068)*IT_0102 + conj(IT_0062
      )*(IT_0053 + IT_0103 + IT_0104) + IT_0062*(IT_0054 + IT_0105))*IT_0110 +
       IT_0112*IT_0117;
    return create_ccomplex_return(IT_0118);
}

