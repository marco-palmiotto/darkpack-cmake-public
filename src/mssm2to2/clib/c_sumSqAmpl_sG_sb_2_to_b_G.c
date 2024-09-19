#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sG_sb_2_to_b_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sG_sb_2_to_b_G(
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
    const creal_t m_sb_2 = param->m_sb_2;
    const creal_t Gamma_b2 = param->Gamma_b2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_sb_01 = param->U_sb_01;
    const ccomplex_t U_sb_11 = param->U_sb_11;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_11;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = pow(m_b, 2);
    const ccomplex_t IT_0004 = pow(m_sG, 2);
    const ccomplex_t IT_0005 = pow(m_sb_2, 2);
    const ccomplex_t IT_0006 = cpow((-2)*s_12 + IT_0003 + -IT_0004 + -IT_0005 
      + -reg_prop, -1);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = -g_s;
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = cpow((-2)*s_23 + IT_0003 + -IT_0004 + IT_0005 +
       reg_prop, -1);
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = s_14*IT_0003;
    const ccomplex_t IT_0018 = s_13*s_34;
    const ccomplex_t IT_0019 = (-2)*IT_0018;
    const ccomplex_t IT_0020 = IT_0017 + IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0015;
    const ccomplex_t IT_0022 = (-21.3333333333333)*IT_0018;
    const ccomplex_t IT_0023 = 10.6666666666667*IT_0017;
    const ccomplex_t IT_0024 = IT_0022 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1.4142135623731)*g_s*U_sb_01;
    const ccomplex_t IT_0026 = IT_0000*IT_0025;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 2*IT_0028;
    const ccomplex_t IT_0030 = IT_0010*IT_0025;
    const ccomplex_t IT_0031 = IT_0012*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0033;
    const ccomplex_t IT_0036 = s_34*IT_0004;
    const ccomplex_t IT_0037 = s_13*s_14;
    const ccomplex_t IT_0038 = IT_0036 + (-2)*IT_0037;
    const ccomplex_t IT_0039 = -IT_0008;
    const ccomplex_t IT_0040 = (-2)*IT_0000;
    const ccomplex_t IT_0041 = cpow((-2)*s_13 + IT_0003 + IT_0004 + -m_sb_2*
      (m_sb_2 + (0 + _Complex_I*-1)*Gamma_b2) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0001*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0040*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0014 + IT_0044;
    const ccomplex_t IT_0046 = IT_0015 + IT_0045;
    const ccomplex_t IT_0047 = conj(IT_0015) + conj(IT_0045);
    const ccomplex_t IT_0048 = -IT_0033;
    const ccomplex_t IT_0049 = IT_0025*IT_0041;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = IT_0040*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*-1)*IT_0032 + -IT_0051;
    const ccomplex_t IT_0053 = IT_0048 + IT_0052;
    const ccomplex_t IT_0054 = conj(IT_0048) + conj(IT_0052);
    const ccomplex_t IT_0055 = -IT_0044;
    const ccomplex_t IT_0056 = -IT_0000;
    const ccomplex_t IT_0057 = IT_0043*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0056;
    const ccomplex_t IT_0060 = s_13*(conj(IT_0009)*IT_0046 + IT_0009*IT_0047 +
       conj(IT_0029)*IT_0053 + IT_0029*IT_0054) + IT_0003*(conj(IT_0029)*IT_0051
       + IT_0029*conj(IT_0051) + conj(IT_0009)*IT_0055 + IT_0009*conj(IT_0055)) 
      + s_34*(conj(IT_0009)*IT_0058 + IT_0009*conj(IT_0058) + conj(IT_0029)
      *IT_0059 + IT_0029*conj(IT_0059));
    const ccomplex_t IT_0061 = -IT_0015;
    const ccomplex_t IT_0062 = 2.66666666666667*s_13*IT_0004;
    const ccomplex_t IT_0063 = (-21.3333333333333)*IT_0037;
    const ccomplex_t IT_0064 = 10.6666666666667*IT_0036;
    const ccomplex_t IT_0065 = IT_0063 + IT_0064;
    const ccomplex_t IT_0066 = 2.66666666666667*s_13*s_14;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*m_sG*IT_0014;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0067;
    const ccomplex_t IT_0070 = 18.6666666666667*s_13;
    const ccomplex_t IT_0071 = 21.3333333333333*s_13*IT_0068*conj(IT_0068);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_sG*IT_0032;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = 0.5*IT_0072;
    const ccomplex_t IT_0075 = IT_0073 + IT_0074;
    const ccomplex_t IT_0076 = conj(IT_0073) + conj(IT_0074);
    const ccomplex_t IT_0077 = 37.3333333333333*s_13;
    const ccomplex_t IT_0078 = m_b*m_sG;
    const ccomplex_t IT_0079 = (-37.3333333333333)*conj(IT_0073);
    const ccomplex_t IT_0080 = IT_0068*IT_0078*IT_0079;
    const ccomplex_t IT_0081 = (-42.6666666666667)*IT_0004;
    const ccomplex_t IT_0082 = s_13*IT_0081;
    const ccomplex_t IT_0083 = pow(s_13, 2);
    const ccomplex_t IT_0084 = (-21.3333333333333)*IT_0083;
    const ccomplex_t IT_0085 = (-10.6666666666667)*IT_0003;
    const ccomplex_t IT_0086 = s_13*IT_0085;
    const ccomplex_t IT_0087 = (-10.6666666666667)*IT_0004;
    const ccomplex_t IT_0088 = s_13*IT_0087;
    const ccomplex_t IT_0089 = IT_0048*IT_0088;
    const ccomplex_t IT_0090 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0091 = s_14*IT_0090;
    const ccomplex_t IT_0092 = 2.66666666666667*IT_0083;
    const ccomplex_t IT_0093 = s_34*IT_0090;
    const ccomplex_t IT_0094 = (-10.6666666666667)*IT_0083;
    const ccomplex_t IT_0095 = (-10.6666666666667)*IT_0009;
    const ccomplex_t IT_0096 = (-10.6666666666667)*conj(IT_0009);
    const ccomplex_t IT_0097 = (-10.6666666666667)*IT_0029;
    const ccomplex_t IT_0098 = (-10.6666666666667)*conj(IT_0029);
    const ccomplex_t IT_0099 = (-2.66666666666667)*s_13;
    const ccomplex_t IT_0100 = s_14*s_34;
    const ccomplex_t IT_0101 = (-21.3333333333333)*conj(IT_0021);
    const ccomplex_t IT_0102 = (-21.3333333333333)*conj(IT_0035);
    const ccomplex_t IT_0103 = IT_0035*((-85.3333333333333)*conj(IT_0028) + (
      -21.3333333333333)*conj(IT_0034) + 170.666666666667*conj(IT_0035)) +
       IT_0028*(42.6666666666667*conj(IT_0028) + 10.6666666666667*conj(IT_0034) 
      + (-85.3333333333333)*conj(IT_0035)) + IT_0039*(10.6666666666667*conj
      (IT_0016) + (-85.3333333333333)*conj(IT_0021) + 42.6666666666667*conj
      (IT_0039)) + IT_0021*((-21.3333333333333)*conj(IT_0016) + 170.666666666667
      *conj(IT_0021) + (-85.3333333333333)*conj(IT_0039)) + IT_0016*
      (170.666666666667*conj(IT_0016) + 10.6666666666667*conj(IT_0039) + IT_0101
      ) + IT_0034*(10.6666666666667*conj(IT_0028) + 170.666666666667*conj
      (IT_0034) + IT_0102);
    const ccomplex_t IT_0104 = m_b*m_sG*s_34;
    const ccomplex_t IT_0105 = 0.375*(2.66666666666667*IT_0016 + (
      -21.3333333333333)*IT_0021)*conj(IT_0029) + 0.375*conj(IT_0009)*
      (2.66666666666667*IT_0034 + (-21.3333333333333)*IT_0035) + (-8)*(conj
      (IT_0016) + -1./8*conj(IT_0021) + (-0.0625)*conj(IT_0039))*IT_0051 + (-8)*
      (IT_0016 + -1./8*IT_0021 + (-0.0625)*IT_0039)*conj(IT_0051) + (0.5*conj
      (IT_0028) + (-8)*conj(IT_0034) + conj(IT_0035))*IT_0055 + (0.5*IT_0028 + (
      -8)*IT_0034 + IT_0035)*conj(IT_0055) + 0.375*conj(IT_0028)*IT_0095 + 0.375
      *IT_0028*IT_0096 + 0.375*conj(IT_0039)*IT_0097 + 0.375*IT_0039*IT_0098 +
       IT_0029*(conj(IT_0016) + 0.375*IT_0101) + IT_0009*(conj(IT_0034) + 0.375
      *IT_0102);
    const ccomplex_t IT_0106 = 2*IT_0033;
    const ccomplex_t IT_0107 = 2*conj(IT_0033);
    const ccomplex_t IT_0108 = m_b*IT_0004;
    const ccomplex_t IT_0109 = (-18.6666666666667)*IT_0061;
    const ccomplex_t IT_0110 = (-9.33333333333333)*IT_0015;
    const ccomplex_t IT_0111 = IT_0109 + IT_0110;
    const ccomplex_t IT_0112 = (-18.6666666666667)*conj(IT_0061);
    const ccomplex_t IT_0113 = (-9.33333333333333)*conj(IT_0015);
    const ccomplex_t IT_0114 = IT_0112 + IT_0113;
    const ccomplex_t IT_0115 = (-9.33333333333333)*IT_0073;
    const ccomplex_t IT_0116 = (-9.33333333333333)*IT_0074;
    const ccomplex_t IT_0117 = IT_0115 + IT_0116;
    const ccomplex_t IT_0118 = (-9.33333333333333)*conj(IT_0073);
    const ccomplex_t IT_0119 = (-9.33333333333333)*conj(IT_0074);
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = (-21.3333333333333)*IT_0033;
    const ccomplex_t IT_0122 = (-21.3333333333333)*conj(IT_0033);
    const ccomplex_t IT_0123 = (-9.33333333333333)*IT_0069;
    const ccomplex_t IT_0124 = (-9.33333333333333)*conj(IT_0069);
    const ccomplex_t IT_0125 = (-18.6666666666667)*conj(IT_0033);
    const ccomplex_t IT_0126 = (-9.33333333333333)*conj(IT_0048);
    const ccomplex_t IT_0127 = IT_0125 + IT_0126;
    const ccomplex_t IT_0128 = 0.375*IT_0069;
    const ccomplex_t IT_0129 = (-18.6666666666667)*IT_0033;
    const ccomplex_t IT_0130 = (-9.33333333333333)*IT_0048;
    const ccomplex_t IT_0131 = IT_0129 + IT_0130;
    const ccomplex_t IT_0132 = 0.375*conj(IT_0069);
    const ccomplex_t IT_0133 = IT_0076*IT_0111 + IT_0075*IT_0114 + conj
      (IT_0045)*IT_0117 + IT_0045*IT_0120 + 2.66666666666667*conj(IT_0068)*
      (IT_0033 + (-3.5)*IT_0048 + (-3.5)*IT_0052 + 0.375*IT_0121) +
       2.66666666666667*IT_0068*(conj(IT_0033) + (-3.5)*conj(IT_0048) + (-3.5)
      *conj(IT_0052) + 0.375*IT_0122) + conj(IT_0052)*IT_0123 + IT_0052*IT_0124 
      + 2.66666666666667*IT_0127*IT_0128 + 2.66666666666667*IT_0131*IT_0132;
    const ccomplex_t IT_0134 = (-74.6666666666667)*IT_0069;
    const ccomplex_t IT_0135 = (-74.6666666666667)*conj(IT_0069);
    const ccomplex_t IT_0136 = (-37.3333333333333)*IT_0073;
    const ccomplex_t IT_0137 = (-37.3333333333333)*IT_0074;
    const ccomplex_t IT_0138 = (-37.3333333333333)*conj(IT_0074);
    const ccomplex_t IT_0139 = (-21.3333333333333)*s_13;
    const ccomplex_t IT_0140 = IT_0015*IT_0087;
    const ccomplex_t IT_0141 = conj(IT_0015)*IT_0087;
    const ccomplex_t IT_0142 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0143 = (-10.6666666666667)*s_34;
    const ccomplex_t IT_0144 = 0.5*conj(IT_0029)*IT_0046 + 0.5*IT_0029*IT_0047
       + 0.5*conj(IT_0009)*IT_0053 + 0.5*IT_0009*IT_0054 + conj(IT_0051)*IT_0061
       + IT_0051*conj(IT_0061) + 0.5*conj(IT_0055)*IT_0106 + 0.5*IT_0055*IT_0107;
    const ccomplex_t IT_0145 = (-74.6666666666667)*conj(IT_0073);
    const ccomplex_t IT_0146 = (-74.6666666666667)*IT_0073;
    const ccomplex_t IT_0147 = IT_0003*(conj(IT_0009)*IT_0051 + IT_0009*conj
      (IT_0051) + conj(IT_0029)*IT_0055 + IT_0029*conj(IT_0055)) + s_34*(conj
      (IT_0029)*IT_0058 + IT_0029*conj(IT_0058) + conj(IT_0009)*IT_0059 +
       IT_0009*conj(IT_0059)) + 2*IT_0004*(conj(IT_0033)*IT_0046 + IT_0033
      *IT_0047 + IT_0054*IT_0061 + IT_0053*conj(IT_0061)) + 2*s_14*(conj(IT_0033
      )*IT_0058 + IT_0033*conj(IT_0058) + conj(IT_0059)*IT_0061 + IT_0059*conj
      (IT_0061)) + 0.75*(conj(IT_0033)*IT_0061 + IT_0033*conj(IT_0061))*IT_0081 
      + 0.75*(conj(IT_0009)*IT_0029 + IT_0009*conj(IT_0029) + conj(IT_0051)
      *IT_0055 + IT_0051*conj(IT_0055))*IT_0085 + 0.75*(conj(IT_0045)*IT_0053 +
       IT_0045*IT_0054)*IT_0087 + 0.75*(IT_0047*IT_0051 + IT_0046*conj(IT_0051) 
      + IT_0054*IT_0055 + IT_0053*conj(IT_0055))*IT_0090 + 0.75*conj(IT_0074)
      *IT_0134 + 0.75*IT_0074*IT_0135 + 0.75*conj(IT_0068)*((-37.3333333333333)
      *IT_0074 + 2*IT_0136 + IT_0137) + 0.75*IT_0068*((-37.3333333333333)*conj
      (IT_0074) + IT_0138) + 0.75*(conj(IT_0009)*IT_0033 + IT_0009*conj(IT_0033)
       + conj(IT_0029)*IT_0061 + IT_0029*conj(IT_0061))*IT_0139 + 0.75*IT_0054
      *IT_0140 + 0.75*IT_0053*IT_0141 + 0.75*(IT_0054*IT_0058 + IT_0053*conj
      (IT_0058) + IT_0047*IT_0059 + IT_0046*conj(IT_0059))*IT_0142 + 0.75*(conj
      (IT_0051)*IT_0058 + IT_0051*conj(IT_0058) + conj(IT_0055)*IT_0059 +
       IT_0055*conj(IT_0059))*IT_0143 + 2*s_13*IT_0144 + 2*IT_0128*IT_0145 + 2
      *IT_0132*IT_0146;
    const ccomplex_t IT_0148 = (-21.3333333333333)*IT_0015;
    const ccomplex_t IT_0149 = (-21.3333333333333)*IT_0045;
    const ccomplex_t IT_0150 = 0.1875*IT_0149;
    const ccomplex_t IT_0151 = IT_0061 + 0.1875*IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = 2.66666666666667*IT_0048;
    const ccomplex_t IT_0153 = 2.66666666666667*IT_0052;
    const ccomplex_t IT_0154 = (-42.6666666666667)*IT_0033;
    const ccomplex_t IT_0155 = IT_0152 + IT_0153 + IT_0154;
    const ccomplex_t IT_0156 = 2.66666666666667*conj(IT_0048);
    const ccomplex_t IT_0157 = 2.66666666666667*conj(IT_0052);
    const ccomplex_t IT_0158 = (-42.6666666666667)*conj(IT_0033);
    const ccomplex_t IT_0159 = IT_0156 + IT_0157 + IT_0158;
    const ccomplex_t IT_0160 = 5.33333333333333*conj(IT_0061);
    const ccomplex_t IT_0161 = (-21.3333333333333)*conj(IT_0015);
    const ccomplex_t IT_0162 = (-21.3333333333333)*conj(IT_0045);
    const ccomplex_t IT_0163 = IT_0160 + IT_0161 + IT_0162;
    const ccomplex_t IT_0164 = 5.33333333333333*IT_0033;
    const ccomplex_t IT_0165 = (-21.3333333333333)*IT_0048;
    const ccomplex_t IT_0166 = (-21.3333333333333)*IT_0052;
    const ccomplex_t IT_0167 = IT_0164 + IT_0165 + IT_0166;
    const ccomplex_t IT_0168 = 5.33333333333333*conj(IT_0033);
    const ccomplex_t IT_0169 = (-21.3333333333333)*conj(IT_0048);
    const ccomplex_t IT_0170 = (-21.3333333333333)*conj(IT_0052);
    const ccomplex_t IT_0171 = IT_0168 + IT_0169 + IT_0170;
    const ccomplex_t IT_0172 = (-42.6666666666667)*IT_0061;
    const ccomplex_t IT_0173 = (-42.6666666666667)*conj(IT_0061);
    const ccomplex_t IT_0174 = 0.375*IT_0173;
    const ccomplex_t IT_0175 = conj(IT_0016)*IT_0151 + 0.1875*conj(IT_0035)
      *IT_0155 + 0.1875*IT_0035*IT_0159 + 0.1875*IT_0016*IT_0163 + 0.1875*conj
      (IT_0034)*IT_0167 + 0.1875*IT_0034*IT_0171 + 0.5*conj(IT_0021)*(IT_0046 +
       0.375*IT_0172) + 0.5*IT_0021*(IT_0047 + IT_0174);
    const ccomplex_t IT_0176 = m_b*m_sG*s_14;
    const ccomplex_t IT_0177 = (-21.3333333333333)*IT_0061;
    const ccomplex_t IT_0178 = (-21.3333333333333)*conj(IT_0061);
    const ccomplex_t IT_0179 = 0.25*conj(IT_0039)*(IT_0053 + 0.75*IT_0121) +
       0.25*IT_0039*(IT_0054 + 0.75*IT_0122) + conj(IT_0034)*IT_0151 + 0.1875
      *conj(IT_0021)*IT_0155 + 0.1875*IT_0021*IT_0159 + 0.1875*IT_0034*IT_0163 +
       0.1875*conj(IT_0016)*IT_0167 + 0.1875*IT_0016*IT_0171 + 0.5*conj(IT_0035)
      *(IT_0046 + 0.375*IT_0172) + 0.5*IT_0035*(IT_0047 + IT_0174) + 0.25*conj
      (IT_0028)*(IT_0046 + 0.75*IT_0177) + 0.25*IT_0028*(IT_0047 + 0.75*IT_0178);
    const ccomplex_t IT_0180 = 0.75*IT_0055;
    const ccomplex_t IT_0181 = (-2.66666666666667)*conj(IT_0055)*(IT_0020
      *IT_0021 + (-0.75)*IT_0016*IT_0024 + -1./2*IT_0017*IT_0039 + (-0.375)
      *IT_0061*IT_0092 + (-0.375)*IT_0058*IT_0093 + (-0.375)*IT_0046*IT_0094 + 
      -1./2*IT_0086*IT_0180);
    const ccomplex_t IT_0182 = 0.75*IT_0051;
    const ccomplex_t IT_0183 = 1.33333333333333*conj(IT_0051)*(IT_0017*IT_0028
       + 1.5*IT_0024*IT_0034 + (-2)*IT_0020*IT_0035 + 0.75*IT_0033*IT_0092 +
       0.75*IT_0059*IT_0093 + 0.75*IT_0053*IT_0094 + IT_0086*IT_0182);
    const ccomplex_t IT_0184 = m_b*s_13;
    const ccomplex_t IT_0185 = (-9.33333333333333)*IT_0009;
    const ccomplex_t IT_0186 = (-9.33333333333333)*conj(IT_0009);
    const ccomplex_t IT_0187 = (-9.33333333333333)*conj(IT_0029);
    const ccomplex_t IT_0188 = 0.75*IT_0069;
    const ccomplex_t IT_0189 = (-9.33333333333333)*IT_0029;
    const ccomplex_t IT_0190 = 0.75*conj(IT_0069);
    const ccomplex_t IT_0191 = 0.75*conj(IT_0055);
    const ccomplex_t IT_0192 = 0.75*conj(IT_0051);
    const ccomplex_t IT_0193 = conj(IT_0068)*(IT_0029 + (-7)*IT_0051 + 0.75
      *IT_0097) + IT_0068*(conj(IT_0029) + (-7)*conj(IT_0051) + 0.75*IT_0098) +
       IT_0120*IT_0180 + IT_0124*IT_0182 + 0.75*IT_0076*IT_0185 + 0.75*IT_0075
      *IT_0186 + IT_0187*IT_0188 + IT_0189*IT_0190 + IT_0117*IT_0191 + IT_0123
      *IT_0192;
    const ccomplex_t IT_0194 = m_sG*IT_0003;
    const ccomplex_t IT_0195 = conj(IT_0068)*(IT_0009 + (-7)*IT_0055 + 0.75
      *IT_0095) + IT_0068*(conj(IT_0009) + (-7)*conj(IT_0055) + 0.75*IT_0096) +
       0.75*(conj(IT_0029) + conj(IT_0051))*IT_0115 + 0.75*(IT_0029 + IT_0051)
      *IT_0118 + IT_0124*IT_0180 + IT_0119*IT_0182 + 0.75*IT_0074*IT_0187 +
       IT_0186*IT_0188 + 0.75*conj(IT_0074)*IT_0189 + IT_0185*IT_0190 + IT_0123
      *IT_0191 + IT_0116*IT_0192;
    const ccomplex_t IT_0196 = m_sG*s_34;
    const ccomplex_t IT_0197 = conj(IT_0016) + conj(IT_0021);
    const ccomplex_t IT_0198 = IT_0016 + IT_0021;
    const ccomplex_t IT_0199 = 10.6666666666667*conj(IT_0068);
    const ccomplex_t IT_0200 = (-85.3333333333333)*conj(IT_0068);
    const ccomplex_t IT_0201 = (-74.6666666666667)*IT_0074;
    const ccomplex_t IT_0202 = (-74.6666666666667)*conj(IT_0074);
    const ccomplex_t IT_0203 = (-9.33333333333333)*IT_0068;
    const ccomplex_t IT_0204 = (-9.33333333333333)*conj(IT_0068);
    const ccomplex_t IT_0205 = 18.6666666666667*conj(IT_0039)*IT_0069 +
       18.6666666666667*IT_0039*conj(IT_0069) + conj(IT_0028)*(18.6666666666667
      *IT_0073 + 18.6666666666667*IT_0074) + IT_0028*(18.6666666666667*conj
      (IT_0073) + 18.6666666666667*conj(IT_0074)) + conj(IT_0059)*IT_0117 +
       IT_0059*IT_0120 + conj(IT_0058)*IT_0123 + IT_0058*IT_0124 + IT_0134
      *IT_0197 + IT_0198*(IT_0135 + IT_0199 + IT_0200) + (conj(IT_0034) + conj
      (IT_0035))*(IT_0146 + IT_0201) + (IT_0034 + IT_0035)*(IT_0145 + IT_0202) +
       8*(conj(IT_0016) + conj(IT_0021) + -1./4*conj(IT_0039) + 1./8*conj
      (IT_0058))*IT_0203 + (-2)*(IT_0039 + -1./2*IT_0058)*IT_0204;
    const ccomplex_t IT_0206 = m_b*s_14;
    const ccomplex_t IT_0207 = 37.3333333333333*IT_0069;
    const ccomplex_t IT_0208 = 37.3333333333333*IT_0073;
    const ccomplex_t IT_0209 = 37.3333333333333*IT_0074;
    const ccomplex_t IT_0210 = 37.3333333333333*conj(IT_0073);
    const ccomplex_t IT_0211 = 37.3333333333333*conj(IT_0074);
    const ccomplex_t IT_0212 = 37.3333333333333*conj(IT_0069);
    const ccomplex_t IT_0213 = 42.6666666666667*conj(IT_0068);
    const ccomplex_t IT_0214 = (-5.33333333333333)*conj(IT_0068);
    const ccomplex_t IT_0215 = (-37.3333333333333)*conj(IT_0028)*IT_0069 + (
      -37.3333333333333)*IT_0028*conj(IT_0069) + conj(IT_0058)*IT_0117 + IT_0058
      *IT_0120 + conj(IT_0059)*IT_0123 + IT_0059*IT_0124 + conj(IT_0039)*
      (IT_0136 + IT_0137) + IT_0039*(IT_0079 + IT_0138) + 4*(conj(IT_0028) + 
      -conj(IT_0034) + -conj(IT_0035) + 1./4*conj(IT_0059))*IT_0203 + 4*(IT_0028
       + 1./4*IT_0059)*IT_0204 + (conj(IT_0034) + conj(IT_0035))*IT_0207 +
       IT_0197*(IT_0208 + IT_0209) + IT_0198*(IT_0210 + IT_0211) + (IT_0034 +
       IT_0035)*(IT_0212 + IT_0213 + IT_0214);
    const ccomplex_t IT_0216 = m_sG*s_13;
    const ccomplex_t IT_0217 = (-18.6666666666667)*conj(IT_0033)*IT_0073 + (
      -18.6666666666667)*IT_0033*conj(IT_0073) + conj(IT_0069)*IT_0111 + IT_0069
      *IT_0114 + IT_0054*IT_0115 + conj(IT_0052)*IT_0116 + IT_0053*IT_0118 +
       IT_0052*IT_0119 + conj(IT_0045)*IT_0123 + IT_0045*IT_0124 + IT_0074
      *IT_0127 + conj(IT_0074)*IT_0131 + (conj(IT_0015) + conj(IT_0045) + (
      -0.285714285714286)*conj(IT_0061) + (-0.107142857142857)*IT_0178)*IT_0203 
      + (IT_0015 + IT_0045 + (-0.285714285714286)*IT_0061 + (-0.107142857142857)
      *IT_0177)*IT_0204;
    const ccomplex_t IT_0218 = (-2.66666666666667)*IT_0009*(conj(IT_0016)
      *IT_0020 + (-0.75)*conj(IT_0021)*IT_0024) + (-2.66666666666667)*IT_0029*
      (IT_0020*conj(IT_0034) + (-0.75)*IT_0024*conj(IT_0035)) + (
      -1.33333333333333)*IT_0038*(conj(IT_0039)*IT_0046 + IT_0039*IT_0047 + conj
      (IT_0028)*IT_0053 + IT_0028*IT_0054) + 1.33333333333333*s_13*IT_0060 + 
      (conj(IT_0033)*IT_0053 + IT_0033*IT_0054 + IT_0047*IT_0061 + IT_0046*conj
      (IT_0061))*IT_0062 + 2*IT_0061*(conj(IT_0039)*IT_0065 + 0.5*conj(IT_0058)
      *IT_0066) + 2*IT_0068*conj(IT_0069)*IT_0070 + 2*IT_0071 + (IT_0069*conj
      (IT_0069) + IT_0075*IT_0076)*IT_0077 + 2*IT_0080 + 2*conj(IT_0061)*
      (IT_0039*IT_0065 + 0.5*IT_0058*IT_0066 + 0.5*IT_0061*IT_0082 + 0.5*IT_0009
      *IT_0084) + (-2.66666666666667)*conj(IT_0009)*(IT_0016*IT_0020 + (-0.75)
      *IT_0021*IT_0024 + (-0.375)*IT_0061*IT_0084 + (-0.375)*IT_0009*IT_0086) + 
      (-2.66666666666667)*conj(IT_0029)*(IT_0020*IT_0034 + (-0.75)*IT_0024
      *IT_0035 + (-0.375)*IT_0033*IT_0084 + (-0.375)*IT_0029*IT_0086) + (IT_0046
      *IT_0047 + IT_0052*IT_0054)*IT_0088 + IT_0054*IT_0089 + (IT_0047*IT_0058 +
       IT_0046*conj(IT_0058) + IT_0054*IT_0059 + IT_0053*conj(IT_0059))*IT_0091 
      + (-2.66666666666667)*IT_0055*(IT_0020*conj(IT_0021) + (-0.75)*conj
      (IT_0016)*IT_0024 + -1./2*IT_0017*conj(IT_0039) + (-0.375)*conj(IT_0061)
      *IT_0092 + (-0.375)*conj(IT_0058)*IT_0093 + (-0.375)*IT_0047*IT_0094) +
       1.33333333333333*IT_0051*(IT_0017*conj(IT_0028) + 1.5*IT_0024*conj
      (IT_0034) + (-2)*IT_0020*conj(IT_0035) + 0.75*conj(IT_0033)*IT_0092 + 0.75
      *conj(IT_0059)*IT_0093 + 0.75*IT_0054*IT_0094) + IT_0017*(conj(IT_0039)
      *IT_0095 + IT_0039*IT_0096 + conj(IT_0028)*IT_0097 + IT_0028*IT_0098) + 2
      *conj(IT_0068)*(IT_0069*IT_0070 + IT_0068*IT_0099) + IT_0100*IT_0103 +
       2.66666666666667*IT_0104*IT_0105 + (conj(IT_0028)*IT_0065 + 0.5*conj
      (IT_0059)*IT_0066)*IT_0106 + (IT_0028*IT_0065 + 0.5*IT_0059*IT_0066 + 0.5
      *IT_0033*IT_0082 + 0.5*IT_0029*IT_0084)*IT_0107 + IT_0108*IT_0133 +
       1.33333333333333*IT_0078*IT_0147 + 5.33333333333333*IT_0036*IT_0175 +
       5.33333333333333*IT_0176*IT_0179 + IT_0181 + IT_0183 + 1.33333333333333
      *IT_0184*IT_0193 + 1.33333333333333*IT_0194*IT_0195 + IT_0196*IT_0205 +
       IT_0206*IT_0215 + IT_0216*IT_0217;
    return create_ccomplex_return(IT_0218);
}

