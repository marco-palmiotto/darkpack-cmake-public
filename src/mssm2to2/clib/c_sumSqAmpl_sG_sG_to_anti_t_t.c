#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sG_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sG_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_t = param->m_t;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = -g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_sG, 2);
    const ccomplex_t IT_0004 = cpow(s_12 + IT_0003 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = 0.25*IT_0007;
    const ccomplex_t IT_0009 = (-0.25)*IT_0007;
    const ccomplex_t IT_0010 = s_13*s_24;
    const ccomplex_t IT_0011 = s_14*s_23;
    const ccomplex_t IT_0012 = pow(m_t, 2);
    const ccomplex_t IT_0013 = s_12*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_00);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = cpow((-2)*s_13 + IT_0003 + IT_0012 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_01);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0003 + IT_0012 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0019 + IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_10);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = IT_0017*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1.4142135623731)*g_s*conj
      (U_st_11);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = IT_0023*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = IT_0031 + IT_0036;
    const ccomplex_t IT_0038 = m_t*m_sG;
    const ccomplex_t IT_0039 = 2*s_24;
    const ccomplex_t IT_0040 = s_23 + IT_0039;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = 0.5*s_14;
    const ccomplex_t IT_0043 = s_13 + IT_0042;
    const ccomplex_t IT_0044 = IT_0038*IT_0043;
    const ccomplex_t IT_0045 = IT_0015*IT_0027;
    const ccomplex_t IT_0046 = IT_0017*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = IT_0021*IT_0032;
    const ccomplex_t IT_0049 = IT_0023*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = -IT_0047 + -IT_0050;
    const ccomplex_t IT_0052 = IT_0014*IT_0028;
    const ccomplex_t IT_0053 = IT_0017*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = IT_0020*IT_0033;
    const ccomplex_t IT_0056 = IT_0023*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = -IT_0054 + -IT_0057;
    const ccomplex_t IT_0059 = 0.5*s_34;
    const ccomplex_t IT_0060 = IT_0012 + IT_0059;
    const ccomplex_t IT_0061 = IT_0003*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0003 + IT_0012 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0016*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0003 + IT_0012 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0022*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = -IT_0064 + -IT_0067;
    const ccomplex_t IT_0069 = IT_0029*IT_0062;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0034*IT_0065;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = -IT_0070 + -IT_0072;
    const ccomplex_t IT_0074 = 0.5*s_24;
    const ccomplex_t IT_0075 = s_23 + IT_0074;
    const ccomplex_t IT_0076 = IT_0038*IT_0075;
    const ccomplex_t IT_0077 = 2*s_14;
    const ccomplex_t IT_0078 = s_13 + IT_0077;
    const ccomplex_t IT_0079 = IT_0038*IT_0078;
    const ccomplex_t IT_0080 = 21.3333333333333*IT_0013;
    const ccomplex_t IT_0081 = 42.6666666666667*IT_0010;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = 85.3333333333333*IT_0013;
    const ccomplex_t IT_0084 = 85.3333333333333*IT_0011;
    const ccomplex_t IT_0085 = 85.3333333333333*IT_0010;
    const ccomplex_t IT_0086 = IT_0083 + IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = (-8)*IT_0086;
    const ccomplex_t IT_0088 = (-16)*s_34 + (-32)*IT_0012;
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = 21.3333333333333*IT_0089;
    const ccomplex_t IT_0091 = 16*s_34 + 32*IT_0012;
    const ccomplex_t IT_0092 = IT_0003*IT_0091;
    const ccomplex_t IT_0093 = (-21.3333333333333)*IT_0092;
    const ccomplex_t IT_0094 = IT_0087 + IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = (-16)*IT_0009*((-10.6666666666667)*IT_0010 + (
      -10.6666666666667)*IT_0011 + (-10.6666666666667)*IT_0013) + ((
      -42.6666666666667)*IT_0026 + 42.6666666666667*conj(IT_0026) + (
      -42.6666666666667)*IT_0037 + 42.6666666666667*conj(IT_0037))*IT_0041 + ((
      -85.3333333333333)*IT_0026 + 85.3333333333333*conj(IT_0026) + (
      -85.3333333333333)*IT_0037 + 85.3333333333333*conj(IT_0037))*IT_0044 + 
      (341.333333333333*IT_0009 + (-85.3333333333333)*IT_0051 + 85.3333333333333
      *conj(IT_0051) + (-85.3333333333333)*IT_0058 + 85.3333333333333*conj
      (IT_0058))*IT_0061 + (10.6666666666667*IT_0068 + (-10.6666666666667)*conj
      (IT_0068) + 10.6666666666667*IT_0073 + (-10.6666666666667)*conj(IT_0073))
      *IT_0076 + 5.33333333333333*(IT_0068 + -conj(IT_0068) + IT_0073 + -conj
      (IT_0073))*IT_0079 + (-2)*(IT_0051 + -conj(IT_0051) + IT_0058 + -conj
      (IT_0058))*IT_0082 + IT_0008*IT_0094;
    const ccomplex_t IT_0096 = IT_0052*IT_0062;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0055*IT_0065;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0097 + IT_0099;
    const ccomplex_t IT_0101 = IT_0045*IT_0062;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = IT_0048*IT_0065;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = 42.6666666666667*IT_0011;
    const ccomplex_t IT_0107 = IT_0080 + IT_0106;
    const ccomplex_t IT_0108 = (IT_0026 + -conj(IT_0026) + IT_0037 + -conj
      (IT_0037))*IT_0041 + 0.1875*(10.6666666666667*IT_0026 + (-10.6666666666667
      )*conj(IT_0026) + 10.6666666666667*IT_0037 + (-10.6666666666667)*conj
      (IT_0037))*IT_0044 + 0.1875*((-85.3333333333333)*IT_0068 +
       85.3333333333333*conj(IT_0068) + (-85.3333333333333)*IT_0073 +
       85.3333333333333*conj(IT_0073))*IT_0076 + 0.1875*((-42.6666666666667)
      *IT_0068 + 42.6666666666667*conj(IT_0068) + (-42.6666666666667)*IT_0073 +
       42.6666666666667*conj(IT_0073))*IT_0079 + 0.1875*IT_0009*IT_0094 + 0.1875
      *IT_0061*((-85.3333333333333)*IT_0100 + 85.3333333333333*conj(IT_0100) + (
      -85.3333333333333)*IT_0105 + 85.3333333333333*conj(IT_0105)) + (-0.375)*
      (IT_0100 + -conj(IT_0100) + IT_0105 + -conj(IT_0105))*IT_0107;
    const ccomplex_t IT_0109 = m_t*m_sG*s_23;
    const ccomplex_t IT_0110 = (-1.33333333333333)*IT_0037;
    const ccomplex_t IT_0111 = 21.3333333333333*IT_0073;
    const ccomplex_t IT_0112 = IT_0110 + IT_0111;
    const ccomplex_t IT_0113 = (-1.33333333333333)*conj(IT_0037);
    const ccomplex_t IT_0114 = 21.3333333333333*conj(IT_0073);
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = 21.3333333333333*IT_0068;
    const ccomplex_t IT_0117 = (-1.33333333333333)*IT_0026;
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = 21.3333333333333*conj(IT_0068);
    const ccomplex_t IT_0120 = (-1.33333333333333)*conj(IT_0026);
    const ccomplex_t IT_0121 = IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = m_t*m_sG*s_14;
    const ccomplex_t IT_0123 = 0.5*IT_0013;
    const ccomplex_t IT_0124 = IT_0010 + IT_0123;
    const ccomplex_t IT_0125 = 10.6666666666667*IT_0009*IT_0051;
    const ccomplex_t IT_0126 = 10.6666666666667*IT_0009*IT_0058;
    const ccomplex_t IT_0127 = IT_0125 + IT_0126;
    const ccomplex_t IT_0128 = (-10.6666666666667)*IT_0009*conj(IT_0051);
    const ccomplex_t IT_0129 = (-10.6666666666667)*IT_0009*conj(IT_0058);
    const ccomplex_t IT_0130 = IT_0011 + IT_0123;
    const ccomplex_t IT_0131 = 10.6666666666667*IT_0008*IT_0100;
    const ccomplex_t IT_0132 = 10.6666666666667*IT_0008*IT_0105;
    const ccomplex_t IT_0133 = (-10.6666666666667)*IT_0008*conj(IT_0100);
    const ccomplex_t IT_0134 = (-10.6666666666667)*IT_0008*conj(IT_0105);
    const ccomplex_t IT_0135 = s_34*IT_0003;
    const ccomplex_t IT_0136 = (-1.33333333333333)*IT_0051;
    const ccomplex_t IT_0137 = (-1.33333333333333)*conj(IT_0051);
    const ccomplex_t IT_0138 = (-1.33333333333333)*IT_0058;
    const ccomplex_t IT_0139 = (-1.33333333333333)*conj(IT_0058);
    const ccomplex_t IT_0140 = s_12*s_34;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0010 + IT_0011 + IT_0141;
    const ccomplex_t IT_0143 = m_t*m_sG*s_24;
    const ccomplex_t IT_0144 = 21.3333333333333*IT_0058;
    const ccomplex_t IT_0145 = 21.3333333333333*IT_0026;
    const ccomplex_t IT_0146 = 21.3333333333333*conj(IT_0026);
    const ccomplex_t IT_0147 = 21.3333333333333*IT_0051;
    const ccomplex_t IT_0148 = m_t*m_sG*s_13;
    const ccomplex_t IT_0149 = 21.3333333333333*conj(IT_0051);
    const ccomplex_t IT_0150 = 21.3333333333333*IT_0100;
    const ccomplex_t IT_0151 = IT_0003*IT_0012;
    const ccomplex_t IT_0152 = 21.3333333333333*conj(IT_0100)*IT_0105 + conj
      (IT_0068)*IT_0110 + IT_0068*IT_0113 + conj(IT_0073)*IT_0118 + IT_0073
      *IT_0121 + conj(IT_0037)*IT_0145 + IT_0037*IT_0146 + conj(IT_0058)*IT_0147
       + IT_0058*IT_0149 + conj(IT_0105)*IT_0150;
    const ccomplex_t IT_0153 = IT_0008*IT_0095 + 5.33333333333333*IT_0009
      *IT_0108 + IT_0109*(conj(IT_0105)*IT_0112 + IT_0105*IT_0115 + conj(IT_0100
      )*IT_0118 + IT_0100*IT_0121) + (conj(IT_0100)*IT_0112 + IT_0100*IT_0115 +
       conj(IT_0105)*IT_0118 + IT_0105*IT_0121)*IT_0122 + IT_0124*(IT_0127 +
       IT_0128 + IT_0129) + IT_0130*(IT_0131 + IT_0132 + IT_0133 + IT_0134) +
       IT_0061*(IT_0127 + IT_0128 + IT_0129 + IT_0131 + IT_0132 + IT_0133 +
       IT_0134) + IT_0135*(conj(IT_0105)*IT_0136 + IT_0105*IT_0137 + conj
      (IT_0100)*IT_0138 + IT_0100*IT_0139) + IT_0013*(conj(IT_0100)*IT_0136 +
       IT_0100*IT_0137 + conj(IT_0105)*IT_0138 + IT_0105*IT_0139) + (conj
      (IT_0073)*IT_0110 + IT_0073*IT_0113 + conj(IT_0068)*IT_0117 + IT_0068
      *IT_0120)*IT_0142 + IT_0143*(21.3333333333333*IT_0037*conj(IT_0058) + conj
      (IT_0068)*IT_0136 + IT_0068*IT_0137 + conj(IT_0073)*IT_0138 + IT_0073
      *IT_0139 + conj(IT_0037)*IT_0144 + conj(IT_0051)*IT_0145 + IT_0051*IT_0146
      ) + IT_0010*(21.3333333333333*IT_0037*conj(IT_0037) + conj(IT_0058)
      *IT_0144 + conj(IT_0026)*IT_0145 + conj(IT_0051)*IT_0147) + IT_0148*(conj
      (IT_0073)*IT_0136 + IT_0073*IT_0137 + conj(IT_0068)*IT_0138 + IT_0068
      *IT_0139 + conj(IT_0058)*IT_0145 + IT_0058*IT_0146 + conj(IT_0037)*IT_0147
       + IT_0037*IT_0149) + IT_0011*(21.3333333333333*IT_0105*conj(IT_0105) +
       conj(IT_0073)*IT_0111 + conj(IT_0068)*IT_0116 + conj(IT_0100)*IT_0150) +
       IT_0151*IT_0152;
    return create_ccomplex_return(IT_0153);
}

