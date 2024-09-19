#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_st_2_to_t_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_st_2_to_t_G(
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
    const creal_t s_34 = param->s_34;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = -g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_11;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_t, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_st_2, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0009;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = cpow((-2)*s_12 + -IT_0004 + -IT_0005 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_01;
    const ccomplex_t IT_0019 = IT_0000*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0022;
    const ccomplex_t IT_0025 = IT_0012*IT_0018;
    const ccomplex_t IT_0026 = IT_0014*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0010*(96*conj(IT_0010) + (-96)*conj(IT_0011)
       + 48*conj(IT_0017)) + IT_0011*((-96)*conj(IT_0010) + 96*conj(IT_0011) + (
      -48)*conj(IT_0017)) + IT_0017*(48*conj(IT_0010) + (-48)*conj(IT_0011) +
       42.6666666666667*conj(IT_0017)) + IT_0023*(96*conj(IT_0023) + (-96)*conj
      (IT_0024) + 48*conj(IT_0027)) + IT_0024*((-96)*conj(IT_0023) + 96*conj
      (IT_0024) + (-48)*conj(IT_0027)) + IT_0027*(48*conj(IT_0023) + (-48)*conj
      (IT_0024) + 42.6666666666667*conj(IT_0027));
    const ccomplex_t IT_0029 = s_14*s_34;
    const ccomplex_t IT_0030 = -IT_0009;
    const ccomplex_t IT_0031 = (-2)*IT_0012;
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0031*IT_0034;
    const ccomplex_t IT_0036 = -IT_0022;
    const ccomplex_t IT_0037 = IT_0018*IT_0032;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = conj(IT_0017)*(IT_0009 + -IT_0030 + IT_0035) +
       IT_0017*(conj(IT_0009) + -conj(IT_0030) + conj(IT_0035)) + -conj(IT_0027)
      *(IT_0022 + -IT_0036 + -IT_0040) + -IT_0027*(conj(IT_0022) + -conj(IT_0036
      ) + -conj(IT_0040));
    const ccomplex_t IT_0042 = s_13*s_14;
    const ccomplex_t IT_0043 = s_34*IT_0004;
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = (-10.6666666666667)*IT_0043;
    const ccomplex_t IT_0047 = 21.3333333333333*IT_0042;
    const ccomplex_t IT_0048 = IT_0046 + IT_0047;
    const ccomplex_t IT_0049 = 10.6666666666667*IT_0043;
    const ccomplex_t IT_0050 = (-21.3333333333333)*IT_0042;
    const ccomplex_t IT_0051 = IT_0049 + IT_0050;
    const ccomplex_t IT_0052 = -IT_0035;
    const ccomplex_t IT_0053 = (-2)*IT_0016;
    const ccomplex_t IT_0054 = 2*IT_0027;
    const ccomplex_t IT_0055 = s_13*s_34;
    const ccomplex_t IT_0056 = (-21.3333333333333)*IT_0055;
    const ccomplex_t IT_0057 = s_14*IT_0003;
    const ccomplex_t IT_0058 = 10.6666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = 12*IT_0036;
    const ccomplex_t IT_0061 = 12*IT_0040;
    const ccomplex_t IT_0062 = (-24)*IT_0022;
    const ccomplex_t IT_0063 = 2*IT_0060 + IT_0061 + IT_0062;
    const ccomplex_t IT_0064 = 12*conj(IT_0036);
    const ccomplex_t IT_0065 = 12*conj(IT_0040);
    const ccomplex_t IT_0066 = (-24)*conj(IT_0022);
    const ccomplex_t IT_0067 = 2*IT_0064 + IT_0065 + IT_0066;
    const ccomplex_t IT_0068 = 12*IT_0035;
    const ccomplex_t IT_0069 = 24*IT_0009;
    const ccomplex_t IT_0070 = (-24)*IT_0030;
    const ccomplex_t IT_0071 = IT_0068 + IT_0069 + IT_0070;
    const ccomplex_t IT_0072 = 12*conj(IT_0035);
    const ccomplex_t IT_0073 = 24*conj(IT_0009);
    const ccomplex_t IT_0074 = (-24)*conj(IT_0030);
    const ccomplex_t IT_0075 = IT_0072 + IT_0073 + IT_0074;
    const ccomplex_t IT_0076 = 24*IT_0022;
    const ccomplex_t IT_0077 = (-12)*IT_0036;
    const ccomplex_t IT_0078 = (-12)*IT_0040;
    const ccomplex_t IT_0079 = IT_0076 + 2*IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = 24*conj(IT_0022);
    const ccomplex_t IT_0081 = (-12)*conj(IT_0036);
    const ccomplex_t IT_0082 = (-12)*conj(IT_0040);
    const ccomplex_t IT_0083 = IT_0080 + 2*IT_0081 + IT_0082;
    const ccomplex_t IT_0084 = 24*IT_0030;
    const ccomplex_t IT_0085 = (-24)*IT_0009;
    const ccomplex_t IT_0086 = (-12)*IT_0035;
    const ccomplex_t IT_0087 = IT_0084 + IT_0085 + IT_0086;
    const ccomplex_t IT_0088 = 24*conj(IT_0030);
    const ccomplex_t IT_0089 = (-24)*conj(IT_0009);
    const ccomplex_t IT_0090 = (-12)*conj(IT_0035);
    const ccomplex_t IT_0091 = IT_0088 + IT_0089 + IT_0090;
    const ccomplex_t IT_0092 = m_t*m_sG*s_14;
    const ccomplex_t IT_0093 = conj(IT_0027)*(12*IT_0009 + (-12)*IT_0030 +
       1.33333333333333*IT_0035) + IT_0027*(12*conj(IT_0009) + (-12)*conj
      (IT_0030) + 1.33333333333333*conj(IT_0035)) + conj(IT_0017)*((-12)*IT_0022
       + 12*IT_0036 + 1.33333333333333*IT_0040) + IT_0017*((-12)*conj(IT_0022) +
       12*conj(IT_0036) + 1.33333333333333*conj(IT_0040)) + conj(IT_0011)
      *IT_0063 + IT_0011*IT_0067 + conj(IT_0024)*IT_0071 + IT_0024*IT_0075 +
       conj(IT_0010)*IT_0079 + IT_0010*IT_0083 + conj(IT_0023)*IT_0087 + IT_0023
      *IT_0091;
    const ccomplex_t IT_0094 = m_t*m_sG*s_34;
    const ccomplex_t IT_0095 = 1.33333333333333*IT_0052;
    const ccomplex_t IT_0096 = (-10.6666666666667)*IT_0053;
    const ccomplex_t IT_0097 = IT_0095 + IT_0096;
    const ccomplex_t IT_0098 = 1.33333333333333*conj(IT_0052);
    const ccomplex_t IT_0099 = (-10.6666666666667)*conj(IT_0053);
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = 1.33333333333333*IT_0039;
    const ccomplex_t IT_0102 = (-10.6666666666667)*IT_0054;
    const ccomplex_t IT_0103 = IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = 1.33333333333333*conj(IT_0039);
    const ccomplex_t IT_0105 = (-10.6666666666667)*conj(IT_0054);
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = ((-12)*conj(IT_0010) + 12*conj(IT_0011))
      *IT_0039 + ((-12)*IT_0010 + 12*IT_0011)*conj(IT_0039) + ((-12)*conj
      (IT_0023) + 12*conj(IT_0024))*IT_0052 + ((-12)*IT_0023 + 12*IT_0024)*conj
      (IT_0052) + (12*conj(IT_0023) + (-12)*conj(IT_0024))*IT_0053 + (12*IT_0023
       + (-12)*IT_0024)*conj(IT_0053) + (12*conj(IT_0010) + (-12)*conj(IT_0011))
      *IT_0054 + (12*IT_0010 + (-12)*IT_0011)*conj(IT_0054) + conj(IT_0027)
      *IT_0097 + IT_0027*IT_0100 + conj(IT_0017)*IT_0103 + IT_0017*IT_0106;
    const ccomplex_t IT_0108 = 2*IT_0036 + IT_0040;
    const ccomplex_t IT_0109 = 12*IT_0004;
    const ccomplex_t IT_0110 = s_13*IT_0109;
    const ccomplex_t IT_0111 = -IT_0012;
    const ccomplex_t IT_0112 = IT_0038*IT_0111;
    const ccomplex_t IT_0113 = 12*s_13;
    const ccomplex_t IT_0114 = s_14*IT_0113;
    const ccomplex_t IT_0115 = pow(s_13, 2);
    const ccomplex_t IT_0116 = 12*IT_0115;
    const ccomplex_t IT_0117 = (-24)*IT_0004;
    const ccomplex_t IT_0118 = s_13*IT_0117;
    const ccomplex_t IT_0119 = (-12)*IT_0115;
    const ccomplex_t IT_0120 = IT_0034*IT_0111;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = 24*IT_0004;
    const ccomplex_t IT_0123 = s_13*IT_0122;
    const ccomplex_t IT_0124 = (-12)*s_13*IT_0004;
    const ccomplex_t IT_0125 = (-12)*s_13;
    const ccomplex_t IT_0126 = s_14*IT_0125;
    const ccomplex_t IT_0127 = 21.3333333333333*IT_0055;
    const ccomplex_t IT_0128 = (-10.6666666666667)*IT_0057;
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = 1.33333333333333*IT_0003;
    const ccomplex_t IT_0131 = s_13*IT_0130;
    const ccomplex_t IT_0132 = 1.33333333333333*s_13;
    const ccomplex_t IT_0133 = s_34*IT_0132;
    const ccomplex_t IT_0134 = IT_0036 + IT_0040;
    const ccomplex_t IT_0135 = 1.33333333333333*IT_0115;
    const ccomplex_t IT_0136 = 10.6666666666667*IT_0115;
    const ccomplex_t IT_0137 = (-10.6666666666667)*IT_0003;
    const ccomplex_t IT_0138 = s_13*IT_0137;
    const ccomplex_t IT_0139 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0140 = s_14*IT_0139;
    const ccomplex_t IT_0141 = s_34*IT_0139;
    const ccomplex_t IT_0142 = (-10.6666666666667)*IT_0115;
    const ccomplex_t IT_0143 = (-10.6666666666667)*IT_0004;
    const ccomplex_t IT_0144 = s_13*IT_0143;
    const ccomplex_t IT_0145 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0146 = s_14*IT_0145;
    const ccomplex_t IT_0147 = (-1.33333333333333)*IT_0115;
    const ccomplex_t IT_0148 = m_t*m_sG;
    const ccomplex_t IT_0149 = IT_0004*(conj(IT_0035)*IT_0036 + IT_0035*conj
      (IT_0036) + (-0.75)*conj(IT_0009)*((-24)*IT_0036 + (-12)*IT_0040) + (-0.75
      )*IT_0009*((-24)*conj(IT_0036) + (-12)*conj(IT_0040)));
    const ccomplex_t IT_0150 = 12*s_14;
    const ccomplex_t IT_0151 = (-12)*s_14;
    const ccomplex_t IT_0152 = conj(IT_0036) + conj(IT_0040);
    const ccomplex_t IT_0153 = IT_0053*IT_0132;
    const ccomplex_t IT_0154 = conj(IT_0053)*IT_0132;
    const ccomplex_t IT_0155 = 1.33333333333333*s_34;
    const ccomplex_t IT_0156 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0157 = (-10.6666666666667)*s_34;
    const ccomplex_t IT_0158 = (-1.33333333333333)*s_14;
    const ccomplex_t IT_0159 = (-0.75)*(conj(IT_0022)*IT_0035 + IT_0022*conj
      (IT_0035) + 2*IT_0030*(conj(IT_0036) + 1./2*conj(IT_0040)) + conj(IT_0030)
      *IT_0108)*IT_0109 + (-0.75)*(conj(IT_0030)*IT_0039 + IT_0030*conj(IT_0039)
       + conj(IT_0009)*IT_0054 + IT_0009*conj(IT_0054))*IT_0113 + (-0.75)*(conj
      (IT_0009)*IT_0039 + IT_0009*conj(IT_0039) + conj(IT_0022)*IT_0053 +
       IT_0022*conj(IT_0053) + conj(IT_0030)*IT_0054 + IT_0030*conj(IT_0054))
      *IT_0125 + IT_0149 + (-0.75)*conj(IT_0030)*(IT_0022*IT_0117 + IT_0112
      *IT_0150) + (-0.75)*IT_0030*(conj(IT_0022)*IT_0117 + conj(IT_0112)*IT_0150
      ) + (-0.75)*conj(IT_0009)*(IT_0022*IT_0122 + IT_0112*IT_0151) + (-0.75)
      *IT_0009*(conj(IT_0022)*IT_0122 + conj(IT_0112)*IT_0151) + (-0.75)*IT_0152
      *IT_0153 + (-0.75)*IT_0134*IT_0154 + (-0.75)*conj(IT_0053)*
      (10.6666666666667*s_13*IT_0036 + IT_0039*IT_0130 + IT_0054*IT_0137 +
       IT_0112*IT_0155) + (-0.75)*IT_0053*(10.6666666666667*s_13*conj(IT_0036) +
       conj(IT_0039)*IT_0130 + conj(IT_0054)*IT_0137 + conj(IT_0112)*IT_0155) + 
      (-0.75)*conj(IT_0035)*(IT_0054*IT_0132 + IT_0039*IT_0139 + IT_0134*IT_0143
       + IT_0112*IT_0156) + (-0.75)*IT_0035*(conj(IT_0054)*IT_0132 + conj
      (IT_0039)*IT_0139 + IT_0143*IT_0152 + conj(IT_0112)*IT_0156) + (-0.75)
      *conj(IT_0052)*(IT_0022*IT_0113 + IT_0054*IT_0130 + IT_0039*IT_0137 +
       IT_0134*IT_0139 + IT_0036*IT_0145 + IT_0112*IT_0157) + (-0.75)*IT_0052*
      (conj(IT_0022)*IT_0113 + conj(IT_0054)*IT_0130 + conj(IT_0039)*IT_0137 +
       conj(IT_0036)*IT_0145 + IT_0139*IT_0152 + conj(IT_0112)*IT_0157) + (-0.75
      )*conj(IT_0121)*(IT_0022*IT_0150 + IT_0054*IT_0155 + IT_0134*IT_0156 +
       IT_0039*IT_0157 + IT_0036*IT_0158) + (-0.75)*IT_0121*(conj(IT_0022)
      *IT_0150 + conj(IT_0054)*IT_0155 + IT_0152*IT_0156 + conj(IT_0039)*IT_0157
       + conj(IT_0036)*IT_0158);
    const ccomplex_t IT_0160 = IT_0036*IT_0144;
    const ccomplex_t IT_0161 = (-1.33333333333333)*s_13*IT_0004;
    const ccomplex_t IT_0162 = IT_0022*IT_0110;
    const ccomplex_t IT_0163 = conj(IT_0011)*IT_0052;
    const ccomplex_t IT_0164 = conj(IT_0024)*IT_0039;
    const ccomplex_t IT_0165 = (-2)*IT_0055;
    const ccomplex_t IT_0166 = IT_0057 + IT_0165;
    const ccomplex_t IT_0167 = IT_0024*conj(IT_0039) + IT_0011*conj(IT_0052) +
       (conj(IT_0010) + -conj(IT_0011))*IT_0053 + (IT_0010 + -IT_0011)*conj
      (IT_0053) + (conj(IT_0023) + -conj(IT_0024))*IT_0054 + (IT_0023 + -IT_0024
      )*conj(IT_0054) + (-0.75)*conj(IT_0010)*IT_0095 + (-0.75)*IT_0010*IT_0098 
      + (-0.75)*conj(IT_0023)*IT_0101 + (-0.75)*IT_0023*IT_0104 + IT_0163 +
       IT_0164;
    const ccomplex_t IT_0168 = IT_0028*IT_0029 + 2.66666666666667*IT_0041
      *IT_0045 + conj(IT_0027)*(IT_0036*IT_0048 + IT_0022*IT_0051) + conj
      (IT_0017)*(IT_0009*IT_0048 + IT_0030*IT_0051) + (conj(IT_0023)*IT_0039 +
       conj(IT_0010)*IT_0052 + conj(IT_0011)*IT_0053 + conj(IT_0024)*IT_0054)
      *IT_0059 + IT_0043*(conj(IT_0024)*IT_0063 + IT_0024*IT_0067 + conj(IT_0011
      )*IT_0071 + IT_0011*IT_0075 + conj(IT_0023)*IT_0079 + IT_0023*IT_0083 +
       conj(IT_0010)*IT_0087 + IT_0010*IT_0091) + IT_0092*IT_0093 + IT_0057*
      (conj(IT_0017)*IT_0097 + IT_0017*IT_0100 + conj(IT_0027)*IT_0103 + IT_0027
      *IT_0106) + IT_0094*IT_0107 + conj(IT_0022)*(IT_0027*IT_0051 + IT_0108
      *IT_0110 + IT_0112*IT_0114 + IT_0039*IT_0116 + IT_0022*IT_0118 + IT_0054
      *IT_0119) + conj(IT_0030)*(IT_0017*IT_0051 + IT_0035*IT_0110 + IT_0052
      *IT_0116 + IT_0030*IT_0118 + IT_0053*IT_0119 + IT_0114*IT_0121 + IT_0009
      *IT_0123) + conj(IT_0009)*(IT_0017*IT_0048 + IT_0053*IT_0116 + IT_0009
      *IT_0118 + IT_0052*IT_0119 + IT_0030*IT_0123 + IT_0035*IT_0124 + IT_0121
      *IT_0126) + conj(IT_0053)*(IT_0011*IT_0059 + IT_0009*IT_0116 + IT_0030
      *IT_0119 + IT_0010*IT_0129 + IT_0052*IT_0131 + IT_0121*IT_0133 + IT_0035
      *IT_0135 + IT_0053*IT_0138) + conj(IT_0054)*(IT_0024*IT_0059 + IT_0022
      *IT_0119 + IT_0023*IT_0129 + IT_0039*IT_0131 + IT_0112*IT_0133 + IT_0134
      *IT_0135 + IT_0036*IT_0136 + IT_0054*IT_0138) + conj(IT_0121)*(IT_0030
      *IT_0114 + IT_0009*IT_0126 + IT_0053*IT_0133 + IT_0035*IT_0140 + IT_0052
      *IT_0141) + conj(IT_0052)*(IT_0010*IT_0059 + IT_0030*IT_0116 + IT_0009
      *IT_0119 + IT_0011*IT_0129 + IT_0053*IT_0131 + IT_0052*IT_0138 + IT_0121
      *IT_0141 + IT_0035*IT_0142) + conj(IT_0035)*(IT_0030*IT_0110 + IT_0009
      *IT_0124 + IT_0053*IT_0135 + IT_0121*IT_0140 + IT_0052*IT_0142 + IT_0035
      *IT_0144) + conj(IT_0112)*(IT_0022*IT_0114 + IT_0054*IT_0133 + IT_0134
      *IT_0140 + IT_0039*IT_0141 + IT_0036*IT_0146) + conj(IT_0039)*(IT_0023
      *IT_0059 + IT_0022*IT_0116 + IT_0024*IT_0129 + IT_0054*IT_0131 + IT_0039
      *IT_0138 + IT_0112*IT_0141 + IT_0134*IT_0142 + IT_0036*IT_0147) + (
      -1.33333333333333)*IT_0148*IT_0159 + IT_0152*(IT_0054*IT_0135 + IT_0112
      *IT_0140 + IT_0039*IT_0142 + IT_0040*IT_0144 + IT_0160 + IT_0036*IT_0161 +
       IT_0162) + conj(IT_0036)*(IT_0027*IT_0048 + IT_0054*IT_0136 + IT_0036
      *IT_0144 + IT_0112*IT_0146 + IT_0039*IT_0147 + IT_0134*IT_0161 + IT_0162) 
      + IT_0129*(conj(IT_0010)*IT_0053 + conj(IT_0023)*IT_0054 + IT_0163 +
       IT_0164) + (-1.33333333333333)*IT_0166*IT_0167;
    return create_ccomplex_return(IT_0168);
}

