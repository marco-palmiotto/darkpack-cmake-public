#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_1_to_b_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_1_to_b_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_b = param->m_b;
    const creal_t m_sG = param->m_sG;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_sb_1 = param->m_sb_1;
    const creal_t Gamma_b1 = param->Gamma_b1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_10;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_b, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sb_1, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0005 
      + -reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = -IT_0008;
    const ccomplex_t IT_0010 = -g_s;
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = s_13*s_14;
    const ccomplex_t IT_0018 = (-21.3333333333333)*IT_0017;
    const ccomplex_t IT_0019 = s_34*IT_0004;
    const ccomplex_t IT_0020 = 10.6666666666667*IT_0019;
    const ccomplex_t IT_0021 = IT_0018 + IT_0020;
    const ccomplex_t IT_0022 = 21.3333333333333*IT_0017;
    const ccomplex_t IT_0023 = (-10.6666666666667)*IT_0019;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_00;
    const ccomplex_t IT_0026 = IT_0000*IT_0025;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = IT_0010*IT_0025;
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = (-2)*IT_0008;
    const ccomplex_t IT_0035 = 0.5*IT_0015;
    const ccomplex_t IT_0036 = 2*IT_0028;
    const ccomplex_t IT_0037 = (-0.5)*IT_0032;
    const ccomplex_t IT_0038 = (-2)*IT_0000;
    const ccomplex_t IT_0039 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sb_1*
      (m_sb_1 + (0 + _Complex_I*-1)*Gamma_b1) + reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0025*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0032;
    const ccomplex_t IT_0044 = IT_0001*IT_0039;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0038*IT_0045;
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = (-0.5)*IT_0015;
    const ccomplex_t IT_0049 = s_13*s_34;
    const ccomplex_t IT_0050 = (-21.3333333333333)*IT_0049;
    const ccomplex_t IT_0051 = s_14*IT_0003;
    const ccomplex_t IT_0052 = 10.6666666666667*IT_0051;
    const ccomplex_t IT_0053 = IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = 21.3333333333333*IT_0049;
    const ccomplex_t IT_0055 = (-10.6666666666667)*IT_0051;
    const ccomplex_t IT_0056 = IT_0054 + IT_0055;
    const ccomplex_t IT_0057 = 1.33333333333333*IT_0047;
    const ccomplex_t IT_0058 = (-10.6666666666667)*IT_0034;
    const ccomplex_t IT_0059 = IT_0057 + IT_0058;
    const ccomplex_t IT_0060 = 1.33333333333333*conj(IT_0047);
    const ccomplex_t IT_0061 = (-10.6666666666667)*conj(IT_0034);
    const ccomplex_t IT_0062 = IT_0060 + IT_0061;
    const ccomplex_t IT_0063 = 1.33333333333333*IT_0042;
    const ccomplex_t IT_0064 = (-10.6666666666667)*IT_0036;
    const ccomplex_t IT_0065 = IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = 1.33333333333333*conj(IT_0042);
    const ccomplex_t IT_0067 = (-10.6666666666667)*conj(IT_0036);
    const ccomplex_t IT_0068 = IT_0066 + IT_0067;
    const ccomplex_t IT_0069 = 12*IT_0016*IT_0019*conj(IT_0048);
    const ccomplex_t IT_0070 = 12*conj(IT_0016)*IT_0019*IT_0048;
    const ccomplex_t IT_0071 = m_b*m_sG*s_14;
    const ccomplex_t IT_0072 = 12*IT_0032*conj(IT_0048)*IT_0071;
    const ccomplex_t IT_0073 = 12*conj(IT_0032)*IT_0048*IT_0071;
    const ccomplex_t IT_0074 = s_14*s_34;
    const ccomplex_t IT_0075 = (-48)*IT_0035*conj(IT_0048)*IT_0074;
    const ccomplex_t IT_0076 = (-48)*conj(IT_0035);
    const ccomplex_t IT_0077 = 96*IT_0035*conj(IT_0035) + IT_0028*
      (42.6666666666667*conj(IT_0028) + (-48)*conj(IT_0037) + 48*conj(IT_0043)) 
      + IT_0043*(48*conj(IT_0028) + (-96)*conj(IT_0037) + 96*conj(IT_0043)) +
       IT_0037*((-48)*conj(IT_0028) + 96*conj(IT_0037) + (-96)*conj(IT_0043)) +
       conj(IT_0009)*(42.6666666666667*IT_0009 + (-48)*IT_0035 + 48*IT_0048) +
       10.6666666666667*IT_0048*conj(IT_0048) + IT_0009*(48*conj(IT_0048) +
       IT_0076);
    const ccomplex_t IT_0078 = IT_0048*IT_0074*IT_0076;
    const ccomplex_t IT_0079 = (-2)*IT_0017;
    const ccomplex_t IT_0080 = IT_0019 + IT_0079;
    const ccomplex_t IT_0081 = 0.75*IT_0021;
    const ccomplex_t IT_0082 = IT_0080 + IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0014 + IT_0046;
    const ccomplex_t IT_0084 = 12*s_13*IT_0004;
    const ccomplex_t IT_0085 = -IT_0000;
    const ccomplex_t IT_0086 = IT_0045*IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = 12*s_13;
    const ccomplex_t IT_0089 = s_14*IT_0088;
    const ccomplex_t IT_0090 = pow(s_13, 2);
    const ccomplex_t IT_0091 = 12*IT_0090;
    const ccomplex_t IT_0092 = (-24)*IT_0004;
    const ccomplex_t IT_0093 = s_13*IT_0092;
    const ccomplex_t IT_0094 = (-12)*IT_0090;
    const ccomplex_t IT_0095 = (0 + _Complex_I*-1)*IT_0031 + -IT_0042;
    const ccomplex_t IT_0096 = IT_0041*IT_0085;
    const ccomplex_t IT_0097 = 42.6666666666667*IT_0048*conj(IT_0048)*IT_0074;
    const ccomplex_t IT_0098 = 1.33333333333333*IT_0003;
    const ccomplex_t IT_0099 = s_13*IT_0098;
    const ccomplex_t IT_0100 = 1.33333333333333*s_13;
    const ccomplex_t IT_0101 = s_34*IT_0100;
    const ccomplex_t IT_0102 = 1.33333333333333*IT_0090;
    const ccomplex_t IT_0103 = 10.6666666666667*IT_0090;
    const ccomplex_t IT_0104 = (-10.6666666666667)*IT_0003;
    const ccomplex_t IT_0105 = s_13*IT_0104;
    const ccomplex_t IT_0106 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0107 = s_14*IT_0106;
    const ccomplex_t IT_0108 = s_34*IT_0106;
    const ccomplex_t IT_0109 = (-1.33333333333333)*s_13;
    const ccomplex_t IT_0110 = s_14*IT_0109;
    const ccomplex_t IT_0111 = (-10.6666666666667)*IT_0090;
    const ccomplex_t IT_0112 = (-1.33333333333333)*IT_0090;
    const ccomplex_t IT_0113 = (-1.33333333333333)*IT_0083;
    const ccomplex_t IT_0114 = (-1.33333333333333)*conj(IT_0083);
    const ccomplex_t IT_0115 = (-1.33333333333333)*IT_0095;
    const ccomplex_t IT_0116 = (-1.33333333333333)*conj(IT_0095);
    const ccomplex_t IT_0117 = (-2)*IT_0049;
    const ccomplex_t IT_0118 = IT_0051 + IT_0117;
    const ccomplex_t IT_0119 = 1.33333333333333*IT_0036;
    const ccomplex_t IT_0120 = (-1.33333333333333)*IT_0042;
    const ccomplex_t IT_0121 = 1.33333333333333*conj(IT_0036);
    const ccomplex_t IT_0122 = (-1.33333333333333)*conj(IT_0042);
    const ccomplex_t IT_0123 = conj(IT_0035)*IT_0047 + IT_0035*conj(IT_0047) +
       -conj(IT_0034)*(IT_0035 + -IT_0048) + -IT_0034*(conj(IT_0035) + -conj
      (IT_0048)) + (-0.75)*conj(IT_0048)*IT_0057 + (-0.75)*IT_0048*IT_0060 +
       conj(IT_0043)*(IT_0036 + (-0.75)*IT_0063) + IT_0043*(conj(IT_0036) + (
      -0.75)*IT_0066) + (-0.75)*conj(IT_0037)*(IT_0119 + IT_0120) + (-0.75)
      *IT_0037*(IT_0121 + IT_0122);
    const ccomplex_t IT_0124 = m_b*m_sG*s_34;
    const ccomplex_t IT_0125 = conj(IT_0035)*((-12)*IT_0036 + 12*IT_0042) +
       IT_0035*((-12)*conj(IT_0036) + 12*conj(IT_0042)) + conj(IT_0034)*((-12)
      *IT_0037 + 12*IT_0043) + IT_0034*((-12)*conj(IT_0037) + 12*conj(IT_0043)) 
      + (12*conj(IT_0037) + (-12)*conj(IT_0043))*IT_0047 + (12*IT_0037 + (-12)
      *IT_0043)*conj(IT_0047) + conj(IT_0028)*IT_0059 + IT_0028*IT_0062 + conj
      (IT_0009)*IT_0065 + IT_0009*IT_0068 + conj(IT_0048)*(10.6666666666667
      *IT_0036 + (-10.6666666666667)*IT_0042 + IT_0119 + IT_0120) + IT_0048*
      (10.6666666666667*conj(IT_0036) + (-10.6666666666667)*conj(IT_0042) +
       IT_0121 + IT_0122);
    const ccomplex_t IT_0126 = m_b*m_sG;
    const ccomplex_t IT_0127 = 12*s_14;
    const ccomplex_t IT_0128 = (-12)*s_13;
    const ccomplex_t IT_0129 = 1.33333333333333*s_34;
    const ccomplex_t IT_0130 = 10.6666666666667*s_13;
    const ccomplex_t IT_0131 = (-10.6666666666667)*IT_0004;
    const ccomplex_t IT_0132 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0133 = (-10.6666666666667)*s_34;
    const ccomplex_t IT_0134 = (-1.33333333333333)*s_14;
    const ccomplex_t IT_0135 = 12*conj(IT_0032);
    const ccomplex_t IT_0136 = (12*IT_0032 + (-1.33333333333333)*IT_0033)*conj
      (IT_0083) + conj(IT_0016)*(12*IT_0033 + 12*IT_0095) + conj(IT_0015)*(12
      *IT_0032 + (-1.33333333333333)*IT_0095) + IT_0016*(12*conj(IT_0033) + 12
      *conj(IT_0095)) + (-1.33333333333333)*IT_0083*(conj(IT_0033) + (-0.75)
      *IT_0135) + (-1.33333333333333)*IT_0015*(conj(IT_0095) + (-0.75)*IT_0135);
    const ccomplex_t IT_0137 = conj(IT_0016)*(IT_0042*IT_0088 + IT_0032
      *IT_0092 + IT_0096*IT_0127 + IT_0036*IT_0128) + IT_0016*(conj(IT_0042)
      *IT_0088 + conj(IT_0032)*IT_0092 + conj(IT_0096)*IT_0127 + conj(IT_0036)
      *IT_0128) + conj(IT_0034)*(IT_0042*IT_0098 + IT_0095*IT_0100 + IT_0036
      *IT_0104 + IT_0032*IT_0128 + IT_0096*IT_0129 + IT_0033*IT_0130) + IT_0034*
      (conj(IT_0042)*IT_0098 + conj(IT_0095)*IT_0100 + conj(IT_0036)*IT_0104 +
       conj(IT_0032)*IT_0128 + conj(IT_0096)*IT_0129 + conj(IT_0033)*IT_0130) +
       conj(IT_0083)*(IT_0036*IT_0100 + IT_0042*IT_0106 + IT_0095*IT_0131 +
       IT_0096*IT_0132) + IT_0083*(conj(IT_0036)*IT_0100 + conj(IT_0042)*IT_0106
       + conj(IT_0095)*IT_0131 + conj(IT_0096)*IT_0132) + conj(IT_0047)*(IT_0032
      *IT_0088 + IT_0036*IT_0098 + IT_0042*IT_0104 + IT_0095*IT_0106 + IT_0033
      *IT_0109 + IT_0096*IT_0133) + IT_0047*(conj(IT_0032)*IT_0088 + conj
      (IT_0036)*IT_0098 + conj(IT_0042)*IT_0104 + conj(IT_0095)*IT_0106 + conj
      (IT_0033)*IT_0109 + conj(IT_0096)*IT_0133) + conj(IT_0087)*(IT_0032
      *IT_0127 + IT_0036*IT_0129 + IT_0095*IT_0132 + IT_0042*IT_0133 + IT_0033
      *IT_0134) + IT_0087*(conj(IT_0032)*IT_0127 + conj(IT_0036)*IT_0129 + conj
      (IT_0095)*IT_0132 + conj(IT_0042)*IT_0133 + conj(IT_0033)*IT_0134) + conj
      (IT_0015)*(IT_0042*IT_0109 + IT_0036*IT_0130 + IT_0033*IT_0131 + IT_0096
      *IT_0134) + IT_0015*(conj(IT_0042)*IT_0109 + conj(IT_0036)*IT_0130 + conj
      (IT_0033)*IT_0131 + conj(IT_0096)*IT_0134) + IT_0004*IT_0136;
    const ccomplex_t IT_0138 = 12*IT_0015;
    const ccomplex_t IT_0139 = 12*IT_0083;
    const ccomplex_t IT_0140 = (-24)*IT_0016;
    const ccomplex_t IT_0141 = IT_0138 + IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = 12*conj(IT_0015);
    const ccomplex_t IT_0143 = 12*conj(IT_0083);
    const ccomplex_t IT_0144 = (-24)*conj(IT_0016);
    const ccomplex_t IT_0145 = IT_0142 + IT_0143 + IT_0144;
    const ccomplex_t IT_0146 = 12*IT_0033;
    const ccomplex_t IT_0147 = 12*IT_0095;
    const ccomplex_t IT_0148 = (-24)*IT_0032;
    const ccomplex_t IT_0149 = IT_0146 + IT_0147 + IT_0148;
    const ccomplex_t IT_0150 = 12*conj(IT_0033);
    const ccomplex_t IT_0151 = 12*conj(IT_0095);
    const ccomplex_t IT_0152 = (-24)*conj(IT_0032);
    const ccomplex_t IT_0153 = IT_0150 + IT_0151 + IT_0152;
    const ccomplex_t IT_0154 = (-12)*conj(IT_0015);
    const ccomplex_t IT_0155 = (-12)*IT_0033;
    const ccomplex_t IT_0156 = (-12)*conj(IT_0033);
    const ccomplex_t IT_0157 = conj(IT_0048)*((-12)*IT_0015 + (
      -10.6666666666667)*IT_0083 + IT_0113) + conj(IT_0035)*IT_0141 + IT_0035
      *IT_0145 + conj(IT_0037)*IT_0149 + IT_0037*IT_0153 + IT_0048*((
      -10.6666666666667)*conj(IT_0083) + IT_0114 + IT_0154) + conj(IT_0043)*(24
      *IT_0032 + (-12)*IT_0095 + IT_0155) + IT_0043*(24*conj(IT_0032) + (-12)
      *conj(IT_0095) + IT_0156);
    const ccomplex_t IT_0158 = conj(IT_0043)*((-12)*IT_0015 + 24*IT_0016 + (
      -12)*IT_0083) + conj(IT_0028)*(10.6666666666667*IT_0015 + (-12)*IT_0016 +
       1.33333333333333*IT_0083) + IT_0028*(10.6666666666667*conj(IT_0015) + (
      -12)*conj(IT_0016) + 1.33333333333333*conj(IT_0083)) + conj(IT_0009)*((-12
      )*IT_0032 + 10.6666666666667*IT_0033 + 1.33333333333333*IT_0095) + IT_0009
      *((-12)*conj(IT_0032) + 10.6666666666667*conj(IT_0033) + 1.33333333333333
      *conj(IT_0095)) + conj(IT_0037)*IT_0141 + IT_0037*IT_0145 + conj(IT_0035)
      *IT_0149 + IT_0035*IT_0153 + IT_0043*(24*conj(IT_0016) + (-12)*conj
      (IT_0083) + IT_0154) + conj(IT_0048)*((-10.6666666666667)*IT_0095 +
       IT_0115 + IT_0155) + IT_0048*((-10.6666666666667)*conj(IT_0095) + IT_0116
       + IT_0156);
    const ccomplex_t IT_0159 = s_13*IT_0131;
    const ccomplex_t IT_0160 = (-1.33333333333333)*s_13*IT_0004;
    const ccomplex_t IT_0161 = IT_0016*IT_0084;
    const ccomplex_t IT_0162 = IT_0032*IT_0084;
    const ccomplex_t IT_0163 = conj(IT_0009)*(IT_0016*IT_0021 + IT_0015
      *IT_0024) + conj(IT_0028)*(IT_0021*IT_0032 + IT_0024*IT_0033) + (IT_0034
      *conj(IT_0035) + IT_0036*conj(IT_0037) + IT_0042*conj(IT_0043) + IT_0047
      *conj(IT_0048))*IT_0053 + (conj(IT_0037)*IT_0042 + IT_0036*conj(IT_0043) +
       conj(IT_0035)*IT_0047 + IT_0034*conj(IT_0048))*IT_0056 + IT_0051*(conj
      (IT_0009)*IT_0059 + IT_0009*IT_0062 + conj(IT_0028)*IT_0065 + IT_0028
      *IT_0068) + 2*IT_0069 + 2*IT_0070 + 2*IT_0072 + 2*IT_0073 + 2*IT_0075 +
       IT_0074*IT_0077 + 2*IT_0078 + 1.33333333333333*conj(IT_0016)*(IT_0009
      *IT_0082 + 0.75*(IT_0015 + IT_0083)*IT_0084 + 0.75*IT_0087*IT_0089 + 0.75
      *IT_0047*IT_0091 + 0.75*IT_0016*IT_0093 + 0.75*IT_0034*IT_0094) +
       1.33333333333333*conj(IT_0032)*(IT_0028*IT_0082 + 0.75*IT_0042*IT_0091 +
       0.75*IT_0032*IT_0093 + 0.75*IT_0036*IT_0094 + 0.75*IT_0084*(IT_0033 +
       IT_0095) + 0.75*IT_0089*IT_0096) + 2*IT_0097 + conj(IT_0034)*(IT_0035
      *IT_0053 + IT_0048*IT_0056 + IT_0016*IT_0094 + IT_0047*IT_0099 + IT_0087
      *IT_0101 + IT_0083*IT_0102 + IT_0015*IT_0103 + IT_0034*IT_0105) + conj
      (IT_0036)*(IT_0037*IT_0053 + IT_0043*IT_0056 + IT_0032*IT_0094 + IT_0042
      *IT_0099 + IT_0096*IT_0101 + IT_0095*IT_0102 + IT_0033*IT_0103 + IT_0036
      *IT_0105) + conj(IT_0087)*(IT_0016*IT_0089 + IT_0034*IT_0101 + IT_0083
      *IT_0107 + IT_0047*IT_0108 + IT_0015*IT_0110) + conj(IT_0096)*(IT_0032
      *IT_0089 + IT_0036*IT_0101 + IT_0095*IT_0107 + IT_0042*IT_0108 + IT_0033
      *IT_0110) + conj(IT_0047)*(IT_0048*IT_0053 + IT_0035*IT_0056 + IT_0016
      *IT_0091 + IT_0034*IT_0099 + IT_0047*IT_0105 + IT_0087*IT_0108 + IT_0083
      *IT_0111 + IT_0015*IT_0112) + conj(IT_0042)*(IT_0043*IT_0053 + IT_0037
      *IT_0056 + IT_0032*IT_0091 + IT_0036*IT_0099 + IT_0042*IT_0105 + IT_0096
      *IT_0108 + IT_0095*IT_0111 + IT_0033*IT_0112) + 1.33333333333333*IT_0080*
      (conj(IT_0009)*(IT_0016 + 0.75*IT_0113) + 0.75*IT_0009*IT_0114 + conj
      (IT_0028)*(IT_0032 + 0.75*IT_0115) + 0.75*IT_0028*IT_0116) + (
      -1.33333333333333)*IT_0118*IT_0123 + IT_0124*IT_0125 + IT_0126*IT_0137 +
       IT_0019*IT_0157 + IT_0071*IT_0158 + conj(IT_0083)*(IT_0034*IT_0102 +
       IT_0087*IT_0107 + IT_0047*IT_0111 + IT_0083*IT_0159 + IT_0015*IT_0160 +
       IT_0161) + conj(IT_0015)*(IT_0009*IT_0024 + IT_0034*IT_0103 + IT_0087
      *IT_0110 + IT_0047*IT_0112 + IT_0015*IT_0159 + IT_0083*IT_0160 + IT_0161) 
      + conj(IT_0095)*(IT_0036*IT_0102 + IT_0096*IT_0107 + IT_0042*IT_0111 +
       IT_0095*IT_0159 + IT_0033*IT_0160 + IT_0162) + conj(IT_0033)*(IT_0024
      *IT_0028 + IT_0036*IT_0103 + IT_0096*IT_0110 + IT_0042*IT_0112 + IT_0033
      *IT_0159 + IT_0095*IT_0160 + IT_0162);
    return create_ccomplex_return(IT_0163);
}

