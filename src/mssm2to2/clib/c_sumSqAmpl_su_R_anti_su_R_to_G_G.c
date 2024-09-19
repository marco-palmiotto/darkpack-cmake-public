#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_R_anti_su_R_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_R_anti_su_R_to_G_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_su_R = param->m_su_R;
    const creal_t Gamma_ur = param->Gamma_ur;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = -IT_0000;
    const ccomplex_t IT_0002 = (-2)*IT_0000;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_su_R
      *Gamma_ur + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*IT_0004;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0001, 2);
    const ccomplex_t IT_0009 = IT_0005*(IT_0003 + -IT_0008);
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = conj(IT_0007) + conj(IT_0010);
    const ccomplex_t IT_0012 = pow(g_s, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_su_R, 2);
    const ccomplex_t IT_0015 = cpow(s_12 + IT_0014 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0016 = s_13 + s_23;
    const ccomplex_t IT_0017 = IT_0001*IT_0016;
    const ccomplex_t IT_0018 = g_s*IT_0017;
    const ccomplex_t IT_0019 = -g_s;
    const ccomplex_t IT_0020 = IT_0002*IT_0019;
    const ccomplex_t IT_0021 = s_12 + IT_0014;
    const ccomplex_t IT_0022 = -s_23;
    const ccomplex_t IT_0023 = IT_0021 + IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0023;
    const ccomplex_t IT_0025 = IT_0001*IT_0019;
    const ccomplex_t IT_0026 = (-0.5)*s_13 + (-0.5)*s_23;
    const ccomplex_t IT_0027 = IT_0021 + IT_0026;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = g_s*s_23;
    const ccomplex_t IT_0030 = IT_0002*IT_0029;
    const ccomplex_t IT_0031 = -IT_0024 + 2*IT_0028 + -IT_0030;
    const ccomplex_t IT_0032 = IT_0018 + IT_0031;
    const ccomplex_t IT_0033 = IT_0015*IT_0032;
    const ccomplex_t IT_0034 = IT_0013 + 0.5*IT_0033;
    const ccomplex_t IT_0035 = 5.33333333333333*s_13;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0013 + (-0.5)*IT_0033;
    const ccomplex_t IT_0038 = (-0.666666666666667)*s_13;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_su_R
      *Gamma_ur + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (IT_0003 + -IT_0008)*IT_0041;
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = IT_0003*IT_0041;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = s_23*IT_0035;
    const ccomplex_t IT_0047 = cpow(IT_0002, 2);
    const ccomplex_t IT_0048 = IT_0041*(IT_0003 + -IT_0047);
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = pow(s_23, 2);
    const ccomplex_t IT_0051 = 5.33333333333333*IT_0050;
    const ccomplex_t IT_0052 = IT_0005*(IT_0003 + -IT_0047);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = s_23*IT_0038;
    const ccomplex_t IT_0055 = (-0.5)*IT_0006;
    const ccomplex_t IT_0056 = pow(s_13, 2);
    const ccomplex_t IT_0057 = (-0.666666666666667)*IT_0056;
    const ccomplex_t IT_0058 = IT_0005*(IT_0003 + (-0.5)*IT_0008 + (-0.5)
      *IT_0047);
    const ccomplex_t IT_0059 = 5.33333333333333*s_12;
    const ccomplex_t IT_0060 = s_13*IT_0059;
    const ccomplex_t IT_0061 = IT_0014*IT_0059;
    const ccomplex_t IT_0062 = IT_0014*IT_0035;
    const ccomplex_t IT_0063 = (-0.5)*IT_0009;
    const ccomplex_t IT_0064 = 5.33333333333333*IT_0056;
    const ccomplex_t IT_0065 = pow(m_su_R, 4);
    const ccomplex_t IT_0066 = 5.33333333333333*IT_0065;
    const ccomplex_t IT_0067 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0068 = IT_0014*IT_0067;
    const ccomplex_t IT_0069 = pow(s_12, 2);
    const ccomplex_t IT_0070 = (-0.666666666666667)*IT_0069;
    const ccomplex_t IT_0071 = IT_0008*IT_0041;
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0014*IT_0038;
    const ccomplex_t IT_0074 = IT_0014*IT_0037 + 1.5*conj(IT_0043)*IT_0057 +
       1.5*conj(IT_0058)*IT_0060 + 1.5*conj(IT_0053)*IT_0061 + 1.5*IT_0011
      *IT_0062 + 1.5*conj(IT_0063)*IT_0064 + 1.5*conj(IT_0055)*IT_0066 + 1.5
      *conj(IT_0045)*IT_0068 + 1.5*conj(IT_0049)*IT_0070 + 1.5*conj(IT_0072)
      *IT_0073;
    const ccomplex_t IT_0075 = (-0.5)*IT_0048;
    const ccomplex_t IT_0076 = 0.5*IT_0042;
    const ccomplex_t IT_0077 = conj(IT_0075) + conj(IT_0076);
    const ccomplex_t IT_0078 = 5.33333333333333*s_23;
    const ccomplex_t IT_0079 = IT_0037*IT_0078;
    const ccomplex_t IT_0080 = s_13*IT_0067;
    const ccomplex_t IT_0081 = IT_0055*IT_0080;
    const ccomplex_t IT_0082 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0083 = IT_0034*IT_0082;
    const ccomplex_t IT_0084 = (-0.666666666666667)*IT_0050;
    const ccomplex_t IT_0085 = 0.666666666666667*s_23;
    const ccomplex_t IT_0086 = g_s*IT_0001;
    const ccomplex_t IT_0087 = IT_0015*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = 1.33333333333333*s_13;
    const ccomplex_t IT_0090 = 5.33333333333333*IT_0014;
    const ccomplex_t IT_0091 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0092 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0093 = IT_0037*((-21.3333333333333)*IT_0037 + conj
      (IT_0045)*IT_0059 + conj(IT_0053)*IT_0067 + IT_0035*conj(IT_0072) + conj
      (IT_0058)*IT_0082 + IT_0058*IT_0085 + IT_0088*IT_0089 + conj(IT_0049)
      *IT_0090 + IT_0072*IT_0091 + conj(IT_0055)*IT_0092);
    const ccomplex_t IT_0094 = g_s*IT_0002;
    const ccomplex_t IT_0095 = IT_0015*(IT_0025 + -IT_0094);
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = 0.666666666666667*s_12;
    const ccomplex_t IT_0098 = s_23*IT_0097;
    const ccomplex_t IT_0099 = IT_0014*IT_0097;
    const ccomplex_t IT_0100 = IT_0049*IT_0099;
    const ccomplex_t IT_0101 = 0.666666666666667*s_13;
    const ccomplex_t IT_0102 = s_23*IT_0101;
    const ccomplex_t IT_0103 = IT_0043*IT_0102;
    const ccomplex_t IT_0104 = IT_0014*IT_0101;
    const ccomplex_t IT_0105 = IT_0014*IT_0085;
    const ccomplex_t IT_0106 = 0.666666666666667*IT_0065;
    const ccomplex_t IT_0107 = -IT_0015*(IT_0086 + -IT_0094);
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = 1.33333333333333*s_12;
    const ccomplex_t IT_0110 = s_23*IT_0109;
    const ccomplex_t IT_0111 = IT_0037*IT_0109;
    const ccomplex_t IT_0112 = -IT_0015*(IT_0020 + -IT_0025 + -IT_0094);
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = 1.33333333333333*IT_0069;
    const ccomplex_t IT_0115 = 0.5*IT_0087;
    const ccomplex_t IT_0116 = 1.33333333333333*s_23;
    const ccomplex_t IT_0117 = IT_0014*IT_0116;
    const ccomplex_t IT_0118 = s_23*IT_0059;
    const ccomplex_t IT_0119 = conj(IT_0055)*IT_0061;
    const ccomplex_t IT_0120 = 5.33333333333333*IT_0069;
    const ccomplex_t IT_0121 = IT_0046*conj(IT_0063);
    const ccomplex_t IT_0122 = IT_0014*IT_0078;
    const ccomplex_t IT_0123 = IT_0015*IT_0025;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0126 = s_13*IT_0125;
    const ccomplex_t IT_0127 = IT_0124*IT_0126;
    const ccomplex_t IT_0128 = 0.5*IT_0107;
    const ccomplex_t IT_0129 = s_23*IT_0125;
    const ccomplex_t IT_0130 = IT_0034*IT_0125;
    const ccomplex_t IT_0131 = (-0.5)*IT_0112;
    const ccomplex_t IT_0132 = (-10.6666666666667)*IT_0069;
    const ccomplex_t IT_0133 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0134 = IT_0014*IT_0133;
    const ccomplex_t IT_0135 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0136 = s_13*IT_0135;
    const ccomplex_t IT_0137 = s_23*IT_0135;
    const ccomplex_t IT_0138 = IT_0014*IT_0135;
    const ccomplex_t IT_0139 = IT_0055*IT_0138;
    const ccomplex_t IT_0140 = (-5.33333333333333)*IT_0069;
    const ccomplex_t IT_0141 = s_23*IT_0091;
    const ccomplex_t IT_0142 = IT_0063*IT_0141;
    const ccomplex_t IT_0143 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0144 = IT_0014*IT_0143;
    const ccomplex_t IT_0145 = (-5.33333333333333)*IT_0065;
    const ccomplex_t IT_0146 = s_23*IT_0067;
    const ccomplex_t IT_0147 = conj(IT_0049)*IT_0068;
    const ccomplex_t IT_0148 = conj(IT_0043)*IT_0054;
    const ccomplex_t IT_0149 = IT_0014*IT_0082;
    const ccomplex_t IT_0150 = (-0.666666666666667)*IT_0065;
    const ccomplex_t IT_0151 = conj(IT_0007)*IT_0060 + IT_0073*conj(IT_0075) +
       IT_0076*IT_0098 + IT_0100 + IT_0103 + IT_0075*IT_0104 + IT_0072*IT_0105 +
       IT_0045*IT_0106 + IT_0108*IT_0110 + IT_0111 + IT_0113*IT_0114 + IT_0115
      *IT_0117 + conj(IT_0058)*IT_0118 + IT_0119 + conj(IT_0053)*IT_0120 +
       IT_0121 + conj(IT_0010)*IT_0122 + IT_0127 + IT_0128*IT_0129 + IT_0130 +
       IT_0131*IT_0132 + IT_0088*IT_0134 + IT_0007*IT_0136 + IT_0058*IT_0137 +
       IT_0139 + IT_0053*IT_0140 + IT_0142 + IT_0010*IT_0144 + IT_0096*IT_0145 +
       conj(IT_0076)*IT_0146 + IT_0147 + IT_0148 + conj(IT_0072)*IT_0149 + conj
      (IT_0045)*IT_0150;
    const ccomplex_t IT_0152 = IT_0075 + IT_0076;
    const ccomplex_t IT_0153 = (-5.33333333333333)*IT_0014;
    const ccomplex_t IT_0154 = IT_0089*IT_0115;
    const ccomplex_t IT_0155 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0156 = IT_0088*IT_0155;
    const ccomplex_t IT_0157 = (-21.3333333333333)*IT_0034 + 5.33333333333333
      *IT_0037 + 0.666666666666667*IT_0014*IT_0049 + conj(IT_0053)*IT_0059 +
       conj(IT_0045)*IT_0067 + IT_0038*conj(IT_0072) + conj(IT_0058)*IT_0078 +
       conj(IT_0055)*IT_0090 + conj(IT_0049)*IT_0092 + IT_0045*IT_0097 + IT_0072
      *IT_0101 + IT_0053*IT_0135 + IT_0058*IT_0143 + IT_0085*IT_0152 + IT_0055
      *IT_0153 + IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0123;
    const ccomplex_t IT_0159 = IT_0014*IT_0109;
    const ccomplex_t IT_0160 = IT_0134 + IT_0159;
    const ccomplex_t IT_0161 = s_13*IT_0097;
    const ccomplex_t IT_0162 = IT_0045*IT_0099;
    const ccomplex_t IT_0163 = 0.666666666666667*IT_0069;
    const ccomplex_t IT_0164 = IT_0055*IT_0104;
    const ccomplex_t IT_0165 = 0.666666666666667*IT_0056;
    const ccomplex_t IT_0166 = 0.5*IT_0095;
    const ccomplex_t IT_0167 = 1.33333333333333*IT_0056;
    const ccomplex_t IT_0168 = 1.33333333333333*IT_0014;
    const ccomplex_t IT_0169 = IT_0037*IT_0168;
    const ccomplex_t IT_0170 = 1.33333333333333*IT_0065;
    const ccomplex_t IT_0171 = conj(IT_0049)*IT_0118;
    const ccomplex_t IT_0172 = s_13*IT_0109;
    const ccomplex_t IT_0173 = IT_0014*IT_0125;
    const ccomplex_t IT_0174 = IT_0037*IT_0155;
    const ccomplex_t IT_0175 = (-10.6666666666667)*IT_0056;
    const ccomplex_t IT_0176 = IT_0049*IT_0137;
    const ccomplex_t IT_0177 = IT_0014*IT_0091;
    const ccomplex_t IT_0178 = (-5.33333333333333)*IT_0056;
    const ccomplex_t IT_0179 = conj(IT_0045)*IT_0068;
    const ccomplex_t IT_0180 = conj(IT_0055)*IT_0073;
    const ccomplex_t IT_0181 = s_23*IT_0128;
    const ccomplex_t IT_0182 = s_23*IT_0155;
    const ccomplex_t IT_0183 = IT_0115 + IT_0158;
    const ccomplex_t IT_0184 = (-10.6666666666667)*IT_0014;
    const ccomplex_t IT_0185 = IT_0034*IT_0184;
    const ccomplex_t IT_0186 = (conj(IT_0010) + conj(IT_0043))*IT_0057 + conj
      (IT_0045)*IT_0060 + conj(IT_0058)*(IT_0054 + IT_0060) + IT_0011*IT_0062 +
       conj(IT_0055)*IT_0066 + conj(IT_0049)*IT_0070 + IT_0064*(conj(IT_0063) +
       conj(IT_0072)) + conj(IT_0072)*IT_0073 + conj(IT_0075)*IT_0080 + conj
      (IT_0076)*(IT_0046 + IT_0080) + IT_0045*IT_0136 + IT_0058*(IT_0102 +
       IT_0136) + IT_0053*(IT_0105 + IT_0138) + IT_0076*IT_0141 + conj(IT_0053)*
      (IT_0061 + IT_0149) + IT_0145*(IT_0055 + IT_0158) + IT_0113*IT_0160 +
       IT_0152*IT_0161 + IT_0162 + IT_0049*IT_0163 + IT_0164 + IT_0043*IT_0165 +
       IT_0159*IT_0166 + IT_0088*IT_0167 + IT_0169 + IT_0124*IT_0170 + IT_0171 +
       IT_0131*(IT_0117 + IT_0173) + IT_0096*(IT_0172 + IT_0173) + IT_0174 +
       IT_0115*IT_0175 + IT_0176 + IT_0007*IT_0177 + IT_0010*(IT_0165 + IT_0177)
       + IT_0063*IT_0178 + IT_0072*(IT_0104 + IT_0178) + IT_0179 + IT_0180 +
       IT_0089*(IT_0034 + IT_0181) + IT_0108*(IT_0172 + IT_0182) + IT_0014*((
      -21.3333333333333)*s_13*IT_0124 + IT_0156 + IT_0089*IT_0183) + IT_0185;
    const ccomplex_t IT_0187 = IT_0055*IT_0099;
    const ccomplex_t IT_0188 = conj(IT_0049)*IT_0061;
    const ccomplex_t IT_0189 = conj(IT_0043)*IT_0046;
    const ccomplex_t IT_0190 = IT_0049*IT_0138;
    const ccomplex_t IT_0191 = IT_0043*IT_0141;
    const ccomplex_t IT_0192 = conj(IT_0055)*IT_0068;
    const ccomplex_t IT_0193 = IT_0054*conj(IT_0063);
    const ccomplex_t IT_0194 = IT_0037*IT_0125;
    const ccomplex_t IT_0195 = IT_0060*conj(IT_0072) + conj(IT_0058)*IT_0073 +
       IT_0062*conj(IT_0076) + conj(IT_0010)*IT_0080 + IT_0063*IT_0102 + IT_0058
      *IT_0104 + IT_0007*IT_0105 + IT_0053*IT_0106 + IT_0034*IT_0109 + conj
      (IT_0075)*IT_0118 + conj(IT_0045)*IT_0120 + IT_0115*IT_0126 + IT_0072
      *IT_0136 + IT_0075*IT_0137 + IT_0045*IT_0140 + IT_0113*IT_0145 + conj
      (IT_0007)*IT_0149 + conj(IT_0053)*IT_0150 + IT_0014*(IT_0089*IT_0128 +
       IT_0108*IT_0155) + IT_0010*IT_0161 + IT_0131*IT_0170 + IT_0088*IT_0172 +
       IT_0076*IT_0177 + IT_0187 + IT_0188 + IT_0189 + IT_0190 + IT_0191 +
       IT_0192 + IT_0193 + IT_0194;
    const ccomplex_t IT_0196 = -IT_0015*(IT_0020 + -IT_0025);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0102 + IT_0137;
    const ccomplex_t IT_0199 = (-10.6666666666667)*IT_0050;
    const ccomplex_t IT_0200 = IT_0117 + IT_0199;
    const ccomplex_t IT_0201 = IT_0117 + IT_0145;
    const ccomplex_t IT_0202 = 1.33333333333333*IT_0050;
    const ccomplex_t IT_0203 = IT_0134 + IT_0202;
    const ccomplex_t IT_0204 = IT_0055*IT_0161;
    const ccomplex_t IT_0205 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0206 = IT_0113 + IT_0166;
    const ccomplex_t IT_0207 = IT_0034*IT_0116;
    const ccomplex_t IT_0208 = IT_0088 + IT_0124;
    const ccomplex_t IT_0209 = IT_0113 + IT_0208;
    const ccomplex_t IT_0210 = IT_0037*IT_0133;
    const ccomplex_t IT_0211 = IT_0049*IT_0144;
    const ccomplex_t IT_0212 = IT_0063 + IT_0076;
    const ccomplex_t IT_0213 = (-5.33333333333333)*IT_0050;
    const ccomplex_t IT_0214 = conj(IT_0055)*IT_0080;
    const ccomplex_t IT_0215 = conj(IT_0053)*IT_0146;
    const ccomplex_t IT_0216 = conj(IT_0010)*IT_0054;
    const ccomplex_t IT_0217 = conj(IT_0058)*IT_0084;
    const ccomplex_t IT_0218 = conj(IT_0045)*IT_0062;
    const ccomplex_t IT_0219 = conj(IT_0053)*IT_0061 + IT_0051*(conj(IT_0063) 
      + conj(IT_0076)) + conj(IT_0043)*IT_0084 + IT_0089*(s_23*IT_0088 + IT_0014
      *IT_0096) + IT_0049*IT_0106 + IT_0011*IT_0118 + conj(IT_0055)*IT_0120 + 
      (conj(IT_0049) + conj(IT_0058))*IT_0122 + IT_0114*IT_0124 + IT_0127 +
       IT_0007*IT_0137 + IT_0053*(IT_0098 + IT_0138) + IT_0055*IT_0140 + IT_0072
      *(IT_0098 + IT_0141) + conj(IT_0072)*(IT_0046 + IT_0146) + IT_0077*IT_0149
       + conj(IT_0049)*IT_0150 + IT_0105*IT_0152 + IT_0162 + IT_0014*IT_0155
      *IT_0166 + IT_0169 + IT_0158*(IT_0110 + IT_0132 + IT_0172) + IT_0096
      *IT_0173 + IT_0131*(IT_0110 + IT_0173) + IT_0045*IT_0177 + IT_0179 +
       IT_0115*(IT_0110 + IT_0182) + IT_0185 + IT_0010*IT_0198 + IT_0108*IT_0200
       + IT_0197*IT_0201 + IT_0128*IT_0203 + IT_0204 + IT_0043*IT_0205 + IT_0058
      *(IT_0144 + IT_0205) + IT_0159*IT_0206 + IT_0207 + IT_0129*IT_0209 +
       IT_0210 + IT_0211 + IT_0212*IT_0213 + IT_0214 + IT_0215 + IT_0216 +
       IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = IT_0046*conj(IT_0075);
    const ccomplex_t IT_0221 = conj(IT_0053)*IT_0080;
    const ccomplex_t IT_0222 = conj(IT_0007)*IT_0057;
    const ccomplex_t IT_0223 = conj(IT_0045)*IT_0122;
    const ccomplex_t IT_0224 = conj(IT_0058)*IT_0146;
    const ccomplex_t IT_0225 = conj(IT_0053)*IT_0070;
    const ccomplex_t IT_0226 = conj(IT_0007)*IT_0080;
    const ccomplex_t IT_0227 = conj(IT_0010)*IT_0149;
    const ccomplex_t IT_0228 = conj(IT_0045)*IT_0066;
    const ccomplex_t IT_0229 = conj(IT_0043)*IT_0046 + conj(IT_0049)*IT_0061 +
       IT_0062*conj(IT_0075) + IT_0088*IT_0105 + conj(IT_0076)*IT_0118 + conj
      (IT_0072)*IT_0122 + IT_0098*IT_0128 + IT_0037*IT_0135 + IT_0108*IT_0137 +
       IT_0115*IT_0144 + IT_0192 + IT_0193 + IT_0224 + IT_0225 + IT_0226 +
       IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = IT_0124*IT_0172;
    const ccomplex_t IT_0231 = IT_0126*IT_0158;
    const ccomplex_t IT_0232 = IT_0062*conj(IT_0075) + IT_0058*IT_0098 +
       IT_0063*IT_0102 + IT_0010*IT_0105 + IT_0034*IT_0109 + IT_0088*IT_0117 +
       conj(IT_0076)*IT_0118 + conj(IT_0072)*IT_0122 + IT_0110*IT_0128 + IT_0108
      *IT_0129 + IT_0114*IT_0131 + IT_0113*IT_0132 + IT_0115*IT_0134 + IT_0076
      *IT_0137 + IT_0072*IT_0144 + IT_0007*IT_0161 + IT_0053*IT_0163 + IT_0145*
      (IT_0045 + IT_0166) + IT_0096*IT_0170 + IT_0075*IT_0177 + IT_0187 +
       IT_0188 + IT_0189 + IT_0190 + IT_0191 + IT_0192 + IT_0193 + IT_0194 +
       IT_0224 + IT_0225 + IT_0226 + IT_0227 + IT_0228 + IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = conj(IT_0045)*IT_0118;
    const ccomplex_t IT_0234 = IT_0051*conj(IT_0075);
    const ccomplex_t IT_0235 = conj(IT_0007)*IT_0054;
    const ccomplex_t IT_0236 = conj(IT_0053)*IT_0073;
    const ccomplex_t IT_0237 = conj(IT_0049)*IT_0122;
    const ccomplex_t IT_0238 = IT_0037*IT_0143;
    const ccomplex_t IT_0239 = conj(IT_0045)*IT_0061;
    const ccomplex_t IT_0240 = conj(IT_0007)*IT_0146;
    const ccomplex_t IT_0241 = conj(IT_0010)*IT_0146;
    const ccomplex_t IT_0242 = conj(IT_0053)*IT_0068;
    const ccomplex_t IT_0243 = conj(IT_0055)*IT_0070;
    const ccomplex_t IT_0244 = conj(IT_0058)*IT_0149;
    const ccomplex_t IT_0245 = conj(IT_0063)*IT_0084;
    const ccomplex_t IT_0246 = conj(IT_0049)*IT_0066;
    const ccomplex_t IT_0247 = 0.5*IT_0196;
    const ccomplex_t IT_0248 = IT_0007 + IT_0010;
    const ccomplex_t IT_0249 = IT_0049*IT_0105;
    const ccomplex_t IT_0250 = IT_0037*IT_0116;
    const ccomplex_t IT_0251 = conj(IT_0055)*IT_0060;
    const ccomplex_t IT_0252 = IT_0034*IT_0133;
    const ccomplex_t IT_0253 = IT_0037*IT_0184;
    const ccomplex_t IT_0254 = IT_0055*IT_0136;
    const ccomplex_t IT_0255 = IT_0045*IT_0138;
    const ccomplex_t IT_0256 = conj(IT_0049)*IT_0149;
    const ccomplex_t IT_0257 = IT_0034*IT_0168;
    const ccomplex_t IT_0258 = conj(IT_0010)*IT_0046 + IT_0051*(conj(IT_0043) 
      + conj(IT_0058)) + IT_0054*conj(IT_0072) + conj(IT_0045)*IT_0073 + conj
      (IT_0076)*IT_0084 + IT_0045*IT_0104 + (conj(IT_0053) + conj(IT_0072))
      *IT_0118 + IT_0077*IT_0122 + IT_0124*IT_0132 + IT_0053*(IT_0099 + IT_0137)
       + IT_0010*IT_0141 + IT_0049*IT_0145 + IT_0144*IT_0152 + IT_0114*IT_0158 +
       (IT_0096 + IT_0131)*IT_0159 + IT_0055*IT_0163 + IT_0014*(IT_0096*IT_0155 
      + IT_0089*IT_0166) + IT_0129*(IT_0131 + IT_0183) + IT_0170*IT_0197 + s_23*
      (IT_0154 + IT_0156 + (-21.3333333333333)*IT_0014*IT_0197) + IT_0072
      *IT_0198 + IT_0128*IT_0200 + IT_0108*IT_0203 + IT_0173*IT_0206 + IT_0110
      *IT_0209 + IT_0205*IT_0212 + IT_0043*IT_0213 + IT_0058*(IT_0105 + IT_0213)
       + IT_0230 + IT_0231 + IT_0239 + IT_0240 + IT_0241 + IT_0242 + IT_0243 +
       IT_0244 + IT_0245 + IT_0246 + IT_0201*IT_0247 + IT_0098*IT_0248 + IT_0249
       + IT_0250 + IT_0251 + IT_0252 + IT_0253 + IT_0254 + IT_0255 + IT_0256 +
       IT_0257;
    const ccomplex_t IT_0259 = IT_0049*IT_0098;
    const ccomplex_t IT_0260 = conj(IT_0007)*IT_0064;
    const ccomplex_t IT_0261 = IT_0055*IT_0177;
    const ccomplex_t IT_0262 = conj(IT_0049)*IT_0146;
    const ccomplex_t IT_0263 = conj(IT_0055)*IT_0062;
    const ccomplex_t IT_0264 = conj(IT_0045)*IT_0149;
    const ccomplex_t IT_0265 = IT_0054*conj(IT_0075);
    const ccomplex_t IT_0266 = IT_0037*IT_0101;
    const ccomplex_t IT_0267 = IT_0034*IT_0091;
    const ccomplex_t IT_0268 = conj(IT_0053)*IT_0062;
    const ccomplex_t IT_0269 = conj(IT_0045)*IT_0146;
    const ccomplex_t IT_0270 = conj(IT_0075)*IT_0084;
    const ccomplex_t IT_0271 = conj(IT_0053)*IT_0122;
    const ccomplex_t IT_0272 = conj(IT_0045)*IT_0080;
    const ccomplex_t IT_0273 = IT_0054*conj(IT_0076);
    const ccomplex_t IT_0274 = IT_0057*conj(IT_0063) + conj(IT_0043)*IT_0064 +
       (IT_0057 + IT_0062)*conj(IT_0072) + conj(IT_0007)*IT_0073 + conj(IT_0010)
      *(IT_0064 + IT_0073) + IT_0060*IT_0077 + conj(IT_0058)*(IT_0046 + IT_0080)
       + IT_0055*IT_0106 + conj(IT_0049)*IT_0120 + IT_0075*IT_0136 + IT_0076*
      (IT_0102 + IT_0136) + IT_0049*IT_0140 + IT_0053*(IT_0099 + IT_0144) +
       IT_0124*IT_0145 + conj(IT_0055)*(IT_0062 + IT_0150) + IT_0096*IT_0159 +
       IT_0131*IT_0160 + IT_0045*IT_0161 + IT_0058*(IT_0141 + IT_0161) + 
      (IT_0063 + IT_0072)*IT_0165 + IT_0115*IT_0167 + IT_0166*IT_0173 + IT_0113*
      (IT_0117 + IT_0173) + IT_0088*IT_0175 + IT_0072*IT_0177 + (IT_0010 +
       IT_0043)*IT_0178 + IT_0155*(IT_0034 + IT_0014*IT_0115 + IT_0181) +
       IT_0089*(IT_0037 + s_23*IT_0108 + IT_0014*IT_0208) + IT_0239 + IT_0242 +
       IT_0104*IT_0248 + IT_0253 + IT_0255 + IT_0257 + IT_0259 + IT_0261 +
       IT_0262 + IT_0271 + IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = IT_0046*conj(IT_0058) + conj(IT_0010)*IT_0064 +
       IT_0057*conj(IT_0072) + IT_0102*IT_0108 + IT_0128*IT_0141 + IT_0115
      *IT_0165 + IT_0088*IT_0178 + IT_0262 + IT_0263 + IT_0266 + IT_0267 +
       IT_0271 + IT_0272 + IT_0273;
    const ccomplex_t IT_0276 = conj(IT_0053)*IT_0066;
    const ccomplex_t IT_0277 = conj(IT_0007)*IT_0122;
    const ccomplex_t IT_0278 = conj(IT_0010)*IT_0060;
    const ccomplex_t IT_0279 = conj(IT_0058)*IT_0062;
    const ccomplex_t IT_0280 = conj(IT_0045)*IT_0070;
    const ccomplex_t IT_0281 = conj(IT_0075)*IT_0146;
    const ccomplex_t IT_0282 = conj(IT_0072)*IT_0080;
    const ccomplex_t IT_0283 = IT_0073*conj(IT_0076);
    const ccomplex_t IT_0284 = conj(IT_0043)*IT_0054 + IT_0046*conj(IT_0063) +
       conj(IT_0049)*IT_0068 + IT_0037*IT_0097 + IT_0104*IT_0108 + IT_0119 +
       IT_0088*IT_0136 + IT_0115*IT_0161 + IT_0128*IT_0177 + IT_0276 + IT_0277 +
       IT_0278 + IT_0279 + IT_0280 + IT_0281 + IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = IT_0075*IT_0098 + IT_0100 + IT_0103 + IT_0076
      *IT_0104 + IT_0111 + IT_0119 + IT_0121 + IT_0088*IT_0126 + IT_0130 +
       IT_0010*IT_0136 + IT_0139 + IT_0142 + IT_0007*IT_0144 + (IT_0053 +
       IT_0131)*IT_0145 + IT_0147 + IT_0148 + IT_0014*(IT_0089*IT_0108 + IT_0128
      *IT_0155) + IT_0072*IT_0161 + IT_0045*IT_0163 + IT_0115*IT_0172 + IT_0058
      *IT_0177 + IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280 + IT_0281 +
       IT_0282 + IT_0283;
    const ccomplex_t IT_0286 = IT_0011*(IT_0036 + IT_0039) + IT_0043*(conj
      (IT_0045)*IT_0046 + conj(IT_0049)*IT_0051 + conj(IT_0053)*IT_0054 + conj
      (IT_0055)*IT_0057) + 0.666666666666667*IT_0055*IT_0074 + IT_0077*(IT_0079 
      + IT_0081 + IT_0083) + IT_0063*(IT_0046*conj(IT_0053) + conj(IT_0045)
      *IT_0054 + conj(IT_0055)*IT_0064 + conj(IT_0049)*IT_0084) + IT_0093 +
       IT_0096*IT_0151 + IT_0034*IT_0157 + IT_0158*IT_0186 + IT_0113*IT_0195 +
       IT_0197*IT_0219 + IT_0072*(IT_0171 + IT_0180 + IT_0220 + IT_0221 +
       IT_0222 + IT_0223) + IT_0115*(IT_0164 + IT_0171 + IT_0174 + IT_0176 +
       IT_0180 + IT_0220 + IT_0221 + IT_0222 + IT_0223) + IT_0045*IT_0229 +
       IT_0166*IT_0232 + IT_0108*(IT_0127 + IT_0204 + IT_0207 + IT_0210 +
       IT_0211 + IT_0214 + IT_0233 + IT_0234 + IT_0235 + IT_0236 + IT_0237) +
       IT_0075*(IT_0046*conj(IT_0072) + IT_0051*conj(IT_0076) + IT_0088*IT_0102 
      + IT_0115*IT_0141 + IT_0128*IT_0205 + IT_0108*IT_0213 + IT_0214 + IT_0215 
      + IT_0216 + IT_0217 + IT_0218 + IT_0237 + IT_0238) + IT_0076*(IT_0214 +
       IT_0233 + IT_0234 + IT_0235 + IT_0236 + IT_0237 + IT_0238) + IT_0049*
      (conj(IT_0043)*IT_0051 + conj(IT_0072)*IT_0118 + IT_0077*IT_0122 + IT_0037
      *IT_0153 + IT_0239 + IT_0240 + IT_0241 + IT_0242 + IT_0243 + IT_0244 +
       IT_0245 + IT_0246) + IT_0247*IT_0258 + IT_0088*(conj(IT_0053)*IT_0060 +
       IT_0259 + IT_0260 + IT_0261 + IT_0262 + IT_0263 + IT_0264 + IT_0265) +
       IT_0010*(conj(IT_0053)*IT_0060 + IT_0260 + IT_0262 + IT_0263 + IT_0264 +
       IT_0265 + IT_0266 + IT_0267) + IT_0058*(conj(IT_0007)*IT_0046 + conj
      (IT_0049)*IT_0149 + IT_0251 + IT_0268 + IT_0269 + IT_0270) + IT_0128*(conj
      (IT_0007)*IT_0046 + IT_0230 + IT_0231 + IT_0249 + IT_0250 + IT_0251 +
       IT_0252 + IT_0254 + IT_0256 + IT_0268 + IT_0269 + IT_0270) + IT_0124
      *IT_0274 + IT_0007*IT_0275 + IT_0053*IT_0284 + IT_0131*IT_0285;
    return create_ccomplex_return(IT_0286);
}

