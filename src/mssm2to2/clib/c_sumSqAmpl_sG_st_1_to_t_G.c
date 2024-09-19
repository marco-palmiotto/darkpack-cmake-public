#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_st_1_to_t_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_st_1_to_t_G(
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
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t IT_0000 = -g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_10;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_t, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_st_1, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + -IT_0004 + -IT_0005 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = (-0.5)*IT_0009;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1.4142135623731)*g_s*U_st_00;
    const ccomplex_t IT_0019 = IT_0000*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = (-0.5)*IT_0022;
    const ccomplex_t IT_0025 = IT_0011*IT_0018;
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0016*((-48)*conj(IT_0010) + 42.6666666666667
      *conj(IT_0016) + 48*conj(IT_0017)) + IT_0017*((-96)*conj(IT_0010) + 48
      *conj(IT_0016) + 96*conj(IT_0017)) + IT_0010*(96*conj(IT_0010) + (-48)
      *conj(IT_0016) + (-96)*conj(IT_0017)) + IT_0023*(96*conj(IT_0023) + (-96)
      *conj(IT_0024) + 48*conj(IT_0027)) + IT_0024*((-96)*conj(IT_0023) + 96
      *conj(IT_0024) + (-48)*conj(IT_0027)) + IT_0027*(48*conj(IT_0023) + (-48)
      *conj(IT_0024) + 42.6666666666667*conj(IT_0027));
    const ccomplex_t IT_0029 = s_14*s_34;
    const ccomplex_t IT_0030 = -IT_0009;
    const ccomplex_t IT_0031 = (-2)*IT_0011;
    const ccomplex_t IT_0032 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0001*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0031*IT_0034;
    const ccomplex_t IT_0036 = -IT_0022;
    const ccomplex_t IT_0037 = IT_0018*IT_0032;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = conj(IT_0016)*(IT_0009 + -IT_0030 + IT_0035) +
       IT_0016*(conj(IT_0009) + -conj(IT_0030) + conj(IT_0035)) + -conj(IT_0027)
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
    const ccomplex_t IT_0052 = 2*IT_0027;
    const ccomplex_t IT_0053 = -IT_0035;
    const ccomplex_t IT_0054 = (-2)*IT_0015;
    const ccomplex_t IT_0055 = s_13*s_34;
    const ccomplex_t IT_0056 = (-21.3333333333333)*IT_0055;
    const ccomplex_t IT_0057 = s_14*IT_0003;
    const ccomplex_t IT_0058 = 10.6666666666667*IT_0057;
    const ccomplex_t IT_0059 = IT_0056 + IT_0058;
    const ccomplex_t IT_0060 = 21.3333333333333*IT_0055;
    const ccomplex_t IT_0061 = (-10.6666666666667)*IT_0057;
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = 12*IT_0036;
    const ccomplex_t IT_0064 = 12*IT_0040;
    const ccomplex_t IT_0065 = (-24)*IT_0022;
    const ccomplex_t IT_0066 = 2*IT_0063 + IT_0064 + IT_0065;
    const ccomplex_t IT_0067 = 12*conj(IT_0036);
    const ccomplex_t IT_0068 = 12*conj(IT_0040);
    const ccomplex_t IT_0069 = (-24)*conj(IT_0022);
    const ccomplex_t IT_0070 = 2*IT_0067 + IT_0068 + IT_0069;
    const ccomplex_t IT_0071 = 12*IT_0035;
    const ccomplex_t IT_0072 = 24*IT_0009;
    const ccomplex_t IT_0073 = (-24)*IT_0030;
    const ccomplex_t IT_0074 = IT_0071 + IT_0072 + IT_0073;
    const ccomplex_t IT_0075 = 12*conj(IT_0035);
    const ccomplex_t IT_0076 = 24*conj(IT_0009);
    const ccomplex_t IT_0077 = (-24)*conj(IT_0030);
    const ccomplex_t IT_0078 = IT_0075 + IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = 24*IT_0022;
    const ccomplex_t IT_0080 = (-12)*IT_0036;
    const ccomplex_t IT_0081 = (-12)*IT_0040;
    const ccomplex_t IT_0082 = IT_0079 + 2*IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = 24*conj(IT_0022);
    const ccomplex_t IT_0084 = (-12)*conj(IT_0036);
    const ccomplex_t IT_0085 = (-12)*conj(IT_0040);
    const ccomplex_t IT_0086 = IT_0083 + 2*IT_0084 + IT_0085;
    const ccomplex_t IT_0087 = 24*IT_0030;
    const ccomplex_t IT_0088 = (-24)*IT_0009;
    const ccomplex_t IT_0089 = (-12)*IT_0035;
    const ccomplex_t IT_0090 = IT_0087 + IT_0088 + IT_0089;
    const ccomplex_t IT_0091 = 24*conj(IT_0030);
    const ccomplex_t IT_0092 = (-24)*conj(IT_0009);
    const ccomplex_t IT_0093 = (-12)*conj(IT_0035);
    const ccomplex_t IT_0094 = IT_0091 + IT_0092 + IT_0093;
    const ccomplex_t IT_0095 = m_t*m_sG*s_14;
    const ccomplex_t IT_0096 = conj(IT_0027)*(12*IT_0009 + (-12)*IT_0030 +
       1.33333333333333*IT_0035) + IT_0027*(12*conj(IT_0009) + (-12)*conj
      (IT_0030) + 1.33333333333333*conj(IT_0035)) + conj(IT_0016)*((-12)*IT_0022
       + 12*IT_0036 + 1.33333333333333*IT_0040) + IT_0016*((-12)*conj(IT_0022) +
       12*conj(IT_0036) + 1.33333333333333*conj(IT_0040)) + conj(IT_0010)
      *IT_0066 + IT_0010*IT_0070 + conj(IT_0024)*IT_0074 + IT_0024*IT_0078 +
       conj(IT_0017)*IT_0082 + IT_0017*IT_0086 + conj(IT_0023)*IT_0090 + IT_0023
      *IT_0094;
    const ccomplex_t IT_0097 = m_t*m_sG*s_34;
    const ccomplex_t IT_0098 = 1.33333333333333*IT_0053;
    const ccomplex_t IT_0099 = (-10.6666666666667)*IT_0054;
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = 1.33333333333333*conj(IT_0053);
    const ccomplex_t IT_0102 = (-10.6666666666667)*conj(IT_0054);
    const ccomplex_t IT_0103 = IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = 1.33333333333333*IT_0039;
    const ccomplex_t IT_0105 = (-10.6666666666667)*IT_0052;
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = 1.33333333333333*conj(IT_0039);
    const ccomplex_t IT_0108 = (-10.6666666666667)*conj(IT_0052);
    const ccomplex_t IT_0109 = IT_0107 + IT_0108;
    const ccomplex_t IT_0110 = (12*conj(IT_0010) + (-12)*conj(IT_0017))
      *IT_0039 + (12*IT_0010 + (-12)*IT_0017)*conj(IT_0039) + ((-12)*conj
      (IT_0010) + 12*conj(IT_0017))*IT_0052 + ((-12)*IT_0010 + 12*IT_0017)*conj
      (IT_0052) + conj(IT_0023)*((-12)*IT_0053 + 12*IT_0054) + conj(IT_0024)*(12
      *IT_0053 + (-12)*IT_0054) + IT_0023*((-12)*conj(IT_0053) + 12*conj(IT_0054
      )) + IT_0024*(12*conj(IT_0053) + (-12)*conj(IT_0054)) + conj(IT_0027)
      *IT_0100 + IT_0027*IT_0103 + conj(IT_0016)*IT_0106 + IT_0016*IT_0109;
    const ccomplex_t IT_0111 = 2*IT_0036 + IT_0040;
    const ccomplex_t IT_0112 = 12*IT_0004;
    const ccomplex_t IT_0113 = s_13*IT_0112;
    const ccomplex_t IT_0114 = -IT_0011;
    const ccomplex_t IT_0115 = IT_0038*IT_0114;
    const ccomplex_t IT_0116 = 12*s_13;
    const ccomplex_t IT_0117 = s_14*IT_0116;
    const ccomplex_t IT_0118 = pow(s_13, 2);
    const ccomplex_t IT_0119 = 12*IT_0118;
    const ccomplex_t IT_0120 = (-24)*IT_0004;
    const ccomplex_t IT_0121 = s_13*IT_0120;
    const ccomplex_t IT_0122 = (-12)*IT_0118;
    const ccomplex_t IT_0123 = IT_0034*IT_0114;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = 24*IT_0004;
    const ccomplex_t IT_0126 = s_13*IT_0125;
    const ccomplex_t IT_0127 = (-12)*s_13*IT_0004;
    const ccomplex_t IT_0128 = (-12)*s_13;
    const ccomplex_t IT_0129 = s_14*IT_0128;
    const ccomplex_t IT_0130 = (-2)*IT_0055;
    const ccomplex_t IT_0131 = IT_0057 + IT_0130;
    const ccomplex_t IT_0132 = conj(IT_0024)*IT_0039 + IT_0024*conj(IT_0039) +
       (conj(IT_0023) + -conj(IT_0024))*IT_0052 + (IT_0023 + -IT_0024)*conj
      (IT_0052) + conj(IT_0010)*IT_0053 + IT_0010*conj(IT_0053) + -(conj(IT_0010
      ) + -conj(IT_0017))*IT_0054 + -(IT_0010 + -IT_0017)*conj(IT_0054) + (-0.75
      )*conj(IT_0017)*IT_0098 + (-0.75)*IT_0017*IT_0101 + (-0.75)*conj(IT_0023)
      *IT_0104 + (-0.75)*IT_0023*IT_0107;
    const ccomplex_t IT_0133 = 1.33333333333333*IT_0003;
    const ccomplex_t IT_0134 = s_13*IT_0133;
    const ccomplex_t IT_0135 = 1.33333333333333*s_13;
    const ccomplex_t IT_0136 = s_34*IT_0135;
    const ccomplex_t IT_0137 = IT_0036 + IT_0040;
    const ccomplex_t IT_0138 = 1.33333333333333*IT_0118;
    const ccomplex_t IT_0139 = 10.6666666666667*IT_0118;
    const ccomplex_t IT_0140 = (-10.6666666666667)*IT_0003;
    const ccomplex_t IT_0141 = s_13*IT_0140;
    const ccomplex_t IT_0142 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0143 = s_14*IT_0142;
    const ccomplex_t IT_0144 = s_34*IT_0142;
    const ccomplex_t IT_0145 = (-10.6666666666667)*IT_0118;
    const ccomplex_t IT_0146 = (-10.6666666666667)*IT_0004;
    const ccomplex_t IT_0147 = s_13*IT_0146;
    const ccomplex_t IT_0148 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0149 = s_14*IT_0148;
    const ccomplex_t IT_0150 = (-1.33333333333333)*IT_0118;
    const ccomplex_t IT_0151 = m_t*m_sG;
    const ccomplex_t IT_0152 = IT_0004*(conj(IT_0035)*IT_0036 + IT_0035*conj
      (IT_0036) + (-0.75)*conj(IT_0009)*((-24)*IT_0036 + (-12)*IT_0040) + (-0.75
      )*IT_0009*((-24)*conj(IT_0036) + (-12)*conj(IT_0040)));
    const ccomplex_t IT_0153 = 12*s_14;
    const ccomplex_t IT_0154 = (-12)*s_14;
    const ccomplex_t IT_0155 = conj(IT_0036) + conj(IT_0040);
    const ccomplex_t IT_0156 = IT_0054*IT_0135;
    const ccomplex_t IT_0157 = conj(IT_0054)*IT_0135;
    const ccomplex_t IT_0158 = 1.33333333333333*s_34;
    const ccomplex_t IT_0159 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0160 = (-10.6666666666667)*s_34;
    const ccomplex_t IT_0161 = (-1.33333333333333)*s_14;
    const ccomplex_t IT_0162 = (-0.75)*(conj(IT_0022)*IT_0035 + IT_0022*conj
      (IT_0035) + 2*IT_0030*(conj(IT_0036) + 1./2*conj(IT_0040)) + conj(IT_0030)
      *IT_0111)*IT_0112 + (-0.75)*(conj(IT_0030)*IT_0039 + IT_0030*conj(IT_0039)
       + conj(IT_0009)*IT_0052 + IT_0009*conj(IT_0052))*IT_0116 + (-0.75)*(conj
      (IT_0009)*IT_0039 + IT_0009*conj(IT_0039) + conj(IT_0030)*IT_0052 +
       IT_0030*conj(IT_0052) + conj(IT_0022)*IT_0054 + IT_0022*conj(IT_0054))
      *IT_0128 + IT_0152 + (-0.75)*conj(IT_0030)*(IT_0022*IT_0120 + IT_0115
      *IT_0153) + (-0.75)*IT_0030*(conj(IT_0022)*IT_0120 + conj(IT_0115)*IT_0153
      ) + (-0.75)*conj(IT_0009)*(IT_0022*IT_0125 + IT_0115*IT_0154) + (-0.75)
      *IT_0009*(conj(IT_0022)*IT_0125 + conj(IT_0115)*IT_0154) + (-0.75)*IT_0155
      *IT_0156 + (-0.75)*IT_0137*IT_0157 + (-0.75)*conj(IT_0054)*
      (10.6666666666667*s_13*IT_0036 + IT_0039*IT_0133 + IT_0052*IT_0140 +
       IT_0115*IT_0158) + (-0.75)*IT_0054*(10.6666666666667*s_13*conj(IT_0036) +
       conj(IT_0039)*IT_0133 + conj(IT_0052)*IT_0140 + conj(IT_0115)*IT_0158) + 
      (-0.75)*conj(IT_0035)*(IT_0052*IT_0135 + IT_0039*IT_0142 + IT_0137*IT_0146
       + IT_0115*IT_0159) + (-0.75)*IT_0035*(conj(IT_0052)*IT_0135 + conj
      (IT_0039)*IT_0142 + IT_0146*IT_0155 + conj(IT_0115)*IT_0159) + (-0.75)
      *conj(IT_0053)*(IT_0022*IT_0116 + IT_0052*IT_0133 + IT_0039*IT_0140 +
       IT_0137*IT_0142 + IT_0036*IT_0148 + IT_0115*IT_0160) + (-0.75)*IT_0053*
      (conj(IT_0022)*IT_0116 + conj(IT_0052)*IT_0133 + conj(IT_0039)*IT_0140 +
       conj(IT_0036)*IT_0148 + IT_0142*IT_0155 + conj(IT_0115)*IT_0160) + (-0.75
      )*conj(IT_0124)*(IT_0022*IT_0153 + IT_0052*IT_0158 + IT_0137*IT_0159 +
       IT_0039*IT_0160 + IT_0036*IT_0161) + (-0.75)*IT_0124*(conj(IT_0022)
      *IT_0153 + conj(IT_0052)*IT_0158 + IT_0155*IT_0159 + conj(IT_0039)*IT_0160
       + conj(IT_0036)*IT_0161);
    const ccomplex_t IT_0163 = IT_0036*IT_0147;
    const ccomplex_t IT_0164 = (-1.33333333333333)*s_13*IT_0004;
    const ccomplex_t IT_0165 = IT_0022*IT_0113;
    const ccomplex_t IT_0166 = IT_0028*IT_0029 + 2.66666666666667*IT_0041
      *IT_0045 + conj(IT_0027)*(IT_0036*IT_0048 + IT_0022*IT_0051) + conj
      (IT_0016)*(IT_0009*IT_0048 + IT_0030*IT_0051) + (conj(IT_0023)*IT_0039 +
       conj(IT_0024)*IT_0052 + conj(IT_0017)*IT_0053 + conj(IT_0010)*IT_0054)
      *IT_0059 + (conj(IT_0024)*IT_0039 + conj(IT_0023)*IT_0052 + conj(IT_0010)
      *IT_0053 + conj(IT_0017)*IT_0054)*IT_0062 + IT_0043*(conj(IT_0024)*IT_0066
       + IT_0024*IT_0070 + conj(IT_0010)*IT_0074 + IT_0010*IT_0078 + conj
      (IT_0023)*IT_0082 + IT_0023*IT_0086 + conj(IT_0017)*IT_0090 + IT_0017
      *IT_0094) + IT_0095*IT_0096 + IT_0057*(conj(IT_0016)*IT_0100 + IT_0016
      *IT_0103 + conj(IT_0027)*IT_0106 + IT_0027*IT_0109) + IT_0097*IT_0110 +
       conj(IT_0022)*(IT_0027*IT_0051 + IT_0111*IT_0113 + IT_0115*IT_0117 +
       IT_0039*IT_0119 + IT_0022*IT_0121 + IT_0052*IT_0122) + conj(IT_0030)*
      (IT_0016*IT_0051 + IT_0035*IT_0113 + IT_0053*IT_0119 + IT_0030*IT_0121 +
       IT_0054*IT_0122 + IT_0117*IT_0124 + IT_0009*IT_0126) + conj(IT_0009)*
      (IT_0016*IT_0048 + IT_0054*IT_0119 + IT_0009*IT_0121 + IT_0053*IT_0122 +
       IT_0030*IT_0126 + IT_0035*IT_0127 + IT_0124*IT_0129) + (-1.33333333333333
      )*IT_0131*IT_0132 + conj(IT_0052)*(IT_0024*IT_0059 + IT_0023*IT_0062 +
       IT_0022*IT_0122 + IT_0039*IT_0134 + IT_0115*IT_0136 + IT_0137*IT_0138 +
       IT_0036*IT_0139 + IT_0052*IT_0141) + conj(IT_0054)*(IT_0010*IT_0059 +
       IT_0017*IT_0062 + IT_0009*IT_0119 + IT_0030*IT_0122 + IT_0053*IT_0134 +
       IT_0124*IT_0136 + IT_0035*IT_0138 + IT_0054*IT_0141) + conj(IT_0124)*
      (IT_0030*IT_0117 + IT_0009*IT_0129 + IT_0054*IT_0136 + IT_0035*IT_0143 +
       IT_0053*IT_0144) + conj(IT_0053)*(IT_0017*IT_0059 + IT_0010*IT_0062 +
       IT_0030*IT_0119 + IT_0009*IT_0122 + IT_0054*IT_0134 + IT_0053*IT_0141 +
       IT_0124*IT_0144 + IT_0035*IT_0145) + conj(IT_0035)*(IT_0030*IT_0113 +
       IT_0009*IT_0127 + IT_0054*IT_0138 + IT_0124*IT_0143 + IT_0053*IT_0145 +
       IT_0035*IT_0147) + conj(IT_0115)*(IT_0022*IT_0117 + IT_0052*IT_0136 +
       IT_0137*IT_0143 + IT_0039*IT_0144 + IT_0036*IT_0149) + conj(IT_0039)*
      (IT_0023*IT_0059 + IT_0024*IT_0062 + IT_0022*IT_0119 + IT_0052*IT_0134 +
       IT_0039*IT_0141 + IT_0115*IT_0144 + IT_0137*IT_0145 + IT_0036*IT_0150) + 
      (-1.33333333333333)*IT_0151*IT_0162 + IT_0155*(IT_0052*IT_0138 + IT_0115
      *IT_0143 + IT_0039*IT_0145 + IT_0040*IT_0147 + IT_0163 + IT_0036*IT_0164 +
       IT_0165) + conj(IT_0036)*(IT_0027*IT_0048 + IT_0052*IT_0139 + IT_0036
      *IT_0147 + IT_0115*IT_0149 + IT_0039*IT_0150 + IT_0137*IT_0164 + IT_0165);
    return create_ccomplex_return(IT_0166);
}

