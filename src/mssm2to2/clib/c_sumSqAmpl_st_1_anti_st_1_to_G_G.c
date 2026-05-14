#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_st_1_anti_st_1_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_st_1_anti_st_1_to_G_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = -IT_0000;
    const ccomplex_t IT_0002 = (-2)*IT_0000;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_st_1
      *Gamma_t1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*IT_0004;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = cpow(IT_0001, 2);
    const ccomplex_t IT_0009 = IT_0005*(IT_0003 + -IT_0008);
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = conj(IT_0007) + conj(IT_0010);
    const ccomplex_t IT_0012 = pow(g_s, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_st_1, 2);
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
    const ccomplex_t IT_0040 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_st_1
      *Gamma_t1 + (-0.5)*reg_prop, -1);
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
    const ccomplex_t IT_0058 = (-0.5)*IT_0048;
    const ccomplex_t IT_0059 = 0.5*IT_0042;
    const ccomplex_t IT_0060 = conj(IT_0058) + conj(IT_0059);
    const ccomplex_t IT_0061 = 5.33333333333333*s_23;
    const ccomplex_t IT_0062 = IT_0037*IT_0061;
    const ccomplex_t IT_0063 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0064 = s_13*IT_0063;
    const ccomplex_t IT_0065 = IT_0055*IT_0064;
    const ccomplex_t IT_0066 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0067 = IT_0034*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0009;
    const ccomplex_t IT_0069 = 5.33333333333333*IT_0056;
    const ccomplex_t IT_0070 = (-0.666666666666667)*IT_0050;
    const ccomplex_t IT_0071 = -IT_0015*(IT_0020 + -IT_0025);
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = 0.666666666666667*s_12;
    const ccomplex_t IT_0074 = s_23*IT_0073;
    const ccomplex_t IT_0075 = IT_0014*IT_0073;
    const ccomplex_t IT_0076 = pow(s_12, 2);
    const ccomplex_t IT_0077 = 0.666666666666667*IT_0076;
    const ccomplex_t IT_0078 = IT_0005*(IT_0003 + (-0.5)*IT_0008 + (-0.5)
      *IT_0047);
    const ccomplex_t IT_0079 = 0.666666666666667*s_23;
    const ccomplex_t IT_0080 = IT_0014*IT_0079;
    const ccomplex_t IT_0081 = 0.666666666666667*IT_0050;
    const ccomplex_t IT_0082 = g_s*IT_0001;
    const ccomplex_t IT_0083 = IT_0015*IT_0082;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = 1.33333333333333*s_12;
    const ccomplex_t IT_0086 = s_23*IT_0085;
    const ccomplex_t IT_0087 = IT_0015*IT_0025;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = IT_0014*IT_0085;
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = 1.33333333333333*IT_0076;
    const ccomplex_t IT_0092 = IT_0088*IT_0091;
    const ccomplex_t IT_0093 = 1.33333333333333*IT_0014;
    const ccomplex_t IT_0094 = IT_0034*IT_0093;
    const ccomplex_t IT_0095 = IT_0008*IT_0041;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = 5.33333333333333*s_12;
    const ccomplex_t IT_0098 = s_23*IT_0097;
    const ccomplex_t IT_0099 = IT_0014*IT_0097;
    const ccomplex_t IT_0100 = conj(IT_0045)*IT_0099;
    const ccomplex_t IT_0101 = IT_0014*IT_0061;
    const ccomplex_t IT_0102 = pow(m_st_1, 4);
    const ccomplex_t IT_0103 = 5.33333333333333*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0083;
    const ccomplex_t IT_0105 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0106 = s_23*IT_0105;
    const ccomplex_t IT_0107 = IT_0088*IT_0106;
    const ccomplex_t IT_0108 = g_s*IT_0002;
    const ccomplex_t IT_0109 = -IT_0015*(IT_0082 + -IT_0108);
    const ccomplex_t IT_0110 = (-0.5)*IT_0109;
    const ccomplex_t IT_0111 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0112 = IT_0014*IT_0111;
    const ccomplex_t IT_0113 = (-0.5)*IT_0071;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (-10.6666666666667)*IT_0014;
    const ccomplex_t IT_0116 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0117 = s_23*IT_0116;
    const ccomplex_t IT_0118 = IT_0014*IT_0116;
    const ccomplex_t IT_0119 = IT_0045*IT_0118;
    const ccomplex_t IT_0120 = IT_0058 + IT_0059;
    const ccomplex_t IT_0121 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0122 = IT_0014*IT_0121;
    const ccomplex_t IT_0123 = (-5.33333333333333)*IT_0050;
    const ccomplex_t IT_0124 = (-5.33333333333333)*IT_0102;
    const ccomplex_t IT_0125 = s_23*IT_0063;
    const ccomplex_t IT_0126 = conj(IT_0007)*IT_0125;
    const ccomplex_t IT_0127 = conj(IT_0010)*IT_0125;
    const ccomplex_t IT_0128 = IT_0014*IT_0063;
    const ccomplex_t IT_0129 = conj(IT_0053)*IT_0128;
    const ccomplex_t IT_0130 = (-0.666666666666667)*IT_0076;
    const ccomplex_t IT_0131 = conj(IT_0055)*IT_0130;
    const ccomplex_t IT_0132 = IT_0014*IT_0066;
    const ccomplex_t IT_0133 = conj(IT_0078)*IT_0132;
    const ccomplex_t IT_0134 = conj(IT_0068)*IT_0070;
    const ccomplex_t IT_0135 = conj(IT_0043)*IT_0051 + (IT_0007 + IT_0010)
      *IT_0074 + IT_0053*IT_0075 + IT_0055*IT_0077 + IT_0078*IT_0080 + IT_0068
      *IT_0081 + IT_0084*IT_0086 + IT_0090 + IT_0092 + IT_0094 + conj(IT_0096)
      *IT_0098 + IT_0100 + IT_0060*IT_0101 + conj(IT_0049)*IT_0103 + IT_0104
      *IT_0106 + IT_0107 + IT_0110*IT_0112 + IT_0114 + IT_0037*IT_0115 + IT_0096
      *IT_0117 + IT_0119 + IT_0120*IT_0122 + IT_0043*IT_0123 + (IT_0049 +
       IT_0072)*IT_0124 + IT_0126 + IT_0127 + IT_0129 + IT_0131 + IT_0133 +
       IT_0134;
    const ccomplex_t IT_0136 = IT_0014*IT_0034 + 1.5*conj(IT_0043)*IT_0051 +
       1.5*conj(IT_0096)*IT_0098 + 1.5*IT_0100 + 1.5*IT_0060*IT_0101 + 1.5*conj
      (IT_0049)*IT_0103 + 1.5*IT_0126 + 1.5*IT_0127 + 1.5*IT_0129 + 1.5*IT_0131 
      + 1.5*IT_0133 + 1.5*IT_0134;
    const ccomplex_t IT_0137 = IT_0015*IT_0020;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = 0.666666666666667*s_13;
    const ccomplex_t IT_0140 = s_23*IT_0139;
    const ccomplex_t IT_0141 = IT_0014*IT_0139;
    const ccomplex_t IT_0142 = IT_0045*IT_0141;
    const ccomplex_t IT_0143 = IT_0049*IT_0080;
    const ccomplex_t IT_0144 = IT_0059*IT_0081;
    const ccomplex_t IT_0145 = 0.5*IT_0087;
    const ccomplex_t IT_0146 = s_13*IT_0085;
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = 1.33333333333333*s_13;
    const ccomplex_t IT_0149 = s_23*IT_0148;
    const ccomplex_t IT_0150 = 1.33333333333333*s_23;
    const ccomplex_t IT_0151 = IT_0014*IT_0150;
    const ccomplex_t IT_0152 = IT_0072*IT_0151;
    const ccomplex_t IT_0153 = IT_0037*IT_0150;
    const ccomplex_t IT_0154 = 1.33333333333333*IT_0050;
    const ccomplex_t IT_0155 = s_13*IT_0097;
    const ccomplex_t IT_0156 = conj(IT_0055)*IT_0155;
    const ccomplex_t IT_0157 = conj(IT_0010)*IT_0046;
    const ccomplex_t IT_0158 = s_13*IT_0105;
    const ccomplex_t IT_0159 = IT_0088*IT_0158;
    const ccomplex_t IT_0160 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0161 = s_23*IT_0160;
    const ccomplex_t IT_0162 = IT_0014*IT_0160;
    const ccomplex_t IT_0163 = IT_0088*IT_0162;
    const ccomplex_t IT_0164 = IT_0034*IT_0111;
    const ccomplex_t IT_0165 = 0.5*IT_0109;
    const ccomplex_t IT_0166 = (-10.6666666666667)*IT_0050;
    const ccomplex_t IT_0167 = s_13*IT_0116;
    const ccomplex_t IT_0168 = IT_0055*IT_0167;
    const ccomplex_t IT_0169 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0170 = s_23*IT_0169;
    const ccomplex_t IT_0171 = IT_0054*conj(IT_0096);
    const ccomplex_t IT_0172 = IT_0014*IT_0038;
    const ccomplex_t IT_0173 = conj(IT_0045)*IT_0172;
    const ccomplex_t IT_0174 = conj(IT_0049)*IT_0132;
    const ccomplex_t IT_0175 = conj(IT_0059)*IT_0070;
    const ccomplex_t IT_0176 = IT_0051*conj(IT_0078) + conj(IT_0053)*IT_0098 +
       IT_0114 + IT_0053*IT_0117 + IT_0078*IT_0123 + IT_0096*IT_0140 + IT_0142 +
       IT_0143 + IT_0144 + IT_0147 + IT_0104*IT_0149 + IT_0152 + IT_0153 +
       IT_0110*IT_0154 + IT_0156 + IT_0157 + IT_0159 + IT_0084*IT_0161 + IT_0163
       + IT_0164 + IT_0165*IT_0166 + IT_0168 + IT_0010*IT_0170 + IT_0171 +
       IT_0173 + IT_0174 + IT_0175;
    const ccomplex_t IT_0177 = 5.33333333333333*IT_0014;
    const ccomplex_t IT_0178 = (-5.33333333333333)*IT_0014;
    const ccomplex_t IT_0179 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0180 = 5.33333333333333*IT_0034 + (-21.3333333333333)
      *IT_0037 + 0.666666666666667*IT_0014*IT_0055 + conj(IT_0053)*IT_0063 +
       IT_0066*conj(IT_0078) + IT_0078*IT_0079 + IT_0035*conj(IT_0096) + conj
      (IT_0045)*IT_0097 + IT_0120*IT_0121 + IT_0096*IT_0169 + conj(IT_0049)
      *IT_0177 + IT_0049*IT_0178 + conj(IT_0055)*IT_0179;
    const ccomplex_t IT_0181 = IT_0015*IT_0108;
    const ccomplex_t IT_0182 = 0.5*IT_0181;
    const ccomplex_t IT_0183 = IT_0059*IT_0074;
    const ccomplex_t IT_0184 = IT_0049*IT_0075;
    const ccomplex_t IT_0185 = IT_0043*IT_0140;
    const ccomplex_t IT_0186 = IT_0058*IT_0141;
    const ccomplex_t IT_0187 = IT_0080*IT_0096;
    const ccomplex_t IT_0188 = 0.666666666666667*IT_0102;
    const ccomplex_t IT_0189 = IT_0045*IT_0188;
    const ccomplex_t IT_0190 = IT_0088*IT_0146;
    const ccomplex_t IT_0191 = IT_0089*IT_0145;
    const ccomplex_t IT_0192 = IT_0037*IT_0085;
    const ccomplex_t IT_0193 = conj(IT_0007)*IT_0155;
    const ccomplex_t IT_0194 = conj(IT_0078)*IT_0098;
    const ccomplex_t IT_0195 = conj(IT_0055)*IT_0099;
    const ccomplex_t IT_0196 = 5.33333333333333*IT_0076;
    const ccomplex_t IT_0197 = IT_0046*conj(IT_0068);
    const ccomplex_t IT_0198 = IT_0014*IT_0105;
    const ccomplex_t IT_0199 = IT_0088*IT_0198;
    const ccomplex_t IT_0200 = IT_0034*IT_0105;
    const ccomplex_t IT_0201 = IT_0007*IT_0167;
    const ccomplex_t IT_0202 = IT_0078*IT_0117;
    const ccomplex_t IT_0203 = IT_0055*IT_0118;
    const ccomplex_t IT_0204 = (-5.33333333333333)*IT_0076;
    const ccomplex_t IT_0205 = IT_0068*IT_0170;
    const ccomplex_t IT_0206 = IT_0010*IT_0122;
    const ccomplex_t IT_0207 = conj(IT_0059)*IT_0125;
    const ccomplex_t IT_0208 = conj(IT_0049)*IT_0128;
    const ccomplex_t IT_0209 = conj(IT_0043)*IT_0054;
    const ccomplex_t IT_0210 = conj(IT_0058)*IT_0172;
    const ccomplex_t IT_0211 = conj(IT_0096)*IT_0132;
    const ccomplex_t IT_0212 = (-0.666666666666667)*IT_0102;
    const ccomplex_t IT_0213 = conj(IT_0045)*IT_0212;
    const ccomplex_t IT_0214 = IT_0072*IT_0089 + conj(IT_0010)*IT_0101 +
       IT_0086*IT_0110 + IT_0084*IT_0112 + IT_0104*IT_0151 + IT_0138*IT_0162 +
       IT_0106*IT_0165 + IT_0124*IT_0182 + IT_0183 + IT_0184 + IT_0185 + IT_0186
       + IT_0187 + IT_0189 + IT_0190 + IT_0191 + IT_0192 + IT_0193 + IT_0194 +
       IT_0195 + conj(IT_0053)*IT_0196 + IT_0197 + IT_0199 + IT_0200 + IT_0201 +
       IT_0202 + IT_0203 + IT_0053*IT_0204 + IT_0205 + IT_0206 + IT_0207 +
       IT_0208 + IT_0209 + IT_0210 + IT_0211 + IT_0213;
    const ccomplex_t IT_0215 = -IT_0015*(IT_0020 + -IT_0025 + -IT_0108);
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = s_13*IT_0073;
    const ccomplex_t IT_0218 = IT_0055*IT_0075;
    const ccomplex_t IT_0219 = IT_0068*IT_0140;
    const ccomplex_t IT_0220 = IT_0034*IT_0085;
    const ccomplex_t IT_0221 = IT_0014*IT_0148;
    const ccomplex_t IT_0222 = conj(IT_0049)*IT_0099;
    const ccomplex_t IT_0223 = conj(IT_0043)*IT_0046;
    const ccomplex_t IT_0224 = IT_0014*IT_0035;
    const ccomplex_t IT_0225 = IT_0037*IT_0105;
    const ccomplex_t IT_0226 = IT_0088*IT_0112;
    const ccomplex_t IT_0227 = IT_0049*IT_0118;
    const ccomplex_t IT_0228 = IT_0043*IT_0170;
    const ccomplex_t IT_0229 = IT_0014*IT_0169;
    const ccomplex_t IT_0230 = conj(IT_0055)*IT_0128;
    const ccomplex_t IT_0231 = IT_0054*conj(IT_0068);
    const ccomplex_t IT_0232 = IT_0072*IT_0198;
    const ccomplex_t IT_0233 = conj(IT_0010)*IT_0064 + IT_0007*IT_0080 +
       IT_0090 + IT_0092 + conj(IT_0058)*IT_0098 + IT_0058*IT_0117 + conj
      (IT_0007)*IT_0132 + IT_0086*IT_0138 + IT_0078*IT_0141 + IT_0084*IT_0146 +
       conj(IT_0096)*IT_0155 + IT_0104*IT_0158 + IT_0110*IT_0162 + IT_0096
      *IT_0167 + conj(IT_0078)*IT_0172 + IT_0091*IT_0182 + IT_0053*IT_0188 +
       conj(IT_0045)*IT_0196 + IT_0045*IT_0204 + conj(IT_0053)*IT_0212 + IT_0124
      *IT_0216 + IT_0010*IT_0217 + IT_0218 + IT_0219 + IT_0220 + IT_0165*IT_0221
       + IT_0222 + IT_0223 + conj(IT_0059)*IT_0224 + IT_0225 + IT_0226 + IT_0227
       + IT_0228 + IT_0059*IT_0229 + IT_0230 + IT_0231 + IT_0232;
    const ccomplex_t IT_0234 = IT_0198 + IT_0221;
    const ccomplex_t IT_0235 = IT_0055*IT_0217;
    const ccomplex_t IT_0236 = IT_0045*IT_0075;
    const ccomplex_t IT_0237 = IT_0086*IT_0088;
    const ccomplex_t IT_0238 = IT_0088*IT_0221;
    const ccomplex_t IT_0239 = IT_0034*IT_0150;
    const ccomplex_t IT_0240 = 1.33333333333333*IT_0102;
    const ccomplex_t IT_0241 = IT_0106*IT_0145;
    const ccomplex_t IT_0242 = (-0.5)*IT_0215;
    const ccomplex_t IT_0243 = IT_0015*(IT_0025 + -IT_0108);
    const ccomplex_t IT_0244 = 0.5*IT_0243;
    const ccomplex_t IT_0245 = IT_0037*IT_0111;
    const ccomplex_t IT_0246 = IT_0034*IT_0115;
    const ccomplex_t IT_0247 = IT_0049*IT_0122;
    const ccomplex_t IT_0248 = conj(IT_0055)*IT_0064;
    const ccomplex_t IT_0249 = conj(IT_0053)*IT_0125;
    const ccomplex_t IT_0250 = conj(IT_0045)*IT_0128;
    const ccomplex_t IT_0251 = conj(IT_0010)*IT_0054;
    const ccomplex_t IT_0252 = IT_0070*conj(IT_0078);
    const ccomplex_t IT_0253 = conj(IT_0053)*IT_0099;
    const ccomplex_t IT_0254 = conj(IT_0045)*IT_0224;
    const ccomplex_t IT_0255 = IT_0051*(conj(IT_0059) + conj(IT_0068)) + conj
      (IT_0043)*IT_0070 + IT_0043*IT_0081 + IT_0037*IT_0093 + IT_0011*IT_0098 + 
      (conj(IT_0049) + conj(IT_0078))*IT_0101 + IT_0007*IT_0117 + IT_0053*
      (IT_0074 + IT_0118) + IT_0080*IT_0120 + IT_0078*(IT_0081 + IT_0122) + 
      (IT_0059 + IT_0068)*IT_0123 + conj(IT_0096)*(IT_0046 + IT_0125) + IT_0060
      *IT_0132 + IT_0010*(IT_0117 + IT_0140) + IT_0084*(IT_0106 + IT_0149) +
       IT_0113*(IT_0124 + IT_0151) + IT_0104*(IT_0086 + IT_0161) + IT_0154
      *IT_0165 + IT_0110*(IT_0151 + IT_0166) + IT_0096*(IT_0074 + IT_0170) +
       IT_0049*IT_0188 + IT_0190 + conj(IT_0055)*IT_0196 + IT_0199 + IT_0055
      *IT_0204 + conj(IT_0049)*IT_0212 + (IT_0089 + IT_0106)*IT_0216 + IT_0045
      *IT_0229 + IT_0182*IT_0234 + IT_0235 + IT_0236 + IT_0237 + IT_0238 +
       IT_0239 + IT_0072*IT_0240 + IT_0241 + (IT_0086 + IT_0198)*IT_0242 + 
      (IT_0089 + IT_0162)*IT_0244 + IT_0245 + IT_0246 + IT_0247 + IT_0248 +
       IT_0249 + IT_0250 + IT_0251 + IT_0252 + IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = conj(IT_0049)*IT_0098;
    const ccomplex_t IT_0257 = IT_0046*conj(IT_0058);
    const ccomplex_t IT_0258 = conj(IT_0053)*IT_0064;
    const ccomplex_t IT_0259 = conj(IT_0055)*IT_0172;
    const ccomplex_t IT_0260 = conj(IT_0007)*IT_0057;
    const ccomplex_t IT_0261 = conj(IT_0045)*IT_0101;
    const ccomplex_t IT_0262 = IT_0055*IT_0141;
    const ccomplex_t IT_0263 = 0.666666666666667*IT_0056;
    const ccomplex_t IT_0264 = IT_0034*IT_0148;
    const ccomplex_t IT_0265 = IT_0037*IT_0160;
    const ccomplex_t IT_0266 = IT_0049*IT_0117;
    const ccomplex_t IT_0267 = conj(IT_0078)*IT_0125;
    const ccomplex_t IT_0268 = conj(IT_0053)*IT_0130;
    const ccomplex_t IT_0269 = conj(IT_0007)*IT_0064;
    const ccomplex_t IT_0270 = conj(IT_0010)*IT_0132;
    const ccomplex_t IT_0271 = conj(IT_0045)*IT_0103;
    const ccomplex_t IT_0272 = conj(IT_0043)*IT_0046 + IT_0034*IT_0073 +
       IT_0080*IT_0084 + conj(IT_0059)*IT_0098 + conj(IT_0049)*IT_0099 + conj
      (IT_0096)*IT_0101 + IT_0037*IT_0116 + IT_0110*IT_0117 + IT_0104*IT_0122 +
       conj(IT_0058)*IT_0224 + IT_0230 + IT_0231 + IT_0267 + IT_0268 + IT_0269 +
       IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = (-10.6666666666667)*IT_0076;
    const ccomplex_t IT_0274 = IT_0053*IT_0077 + IT_0074*IT_0078 + IT_0010
      *IT_0080 + IT_0090 + conj(IT_0059)*IT_0098 + conj(IT_0096)*IT_0101 +
       IT_0106*IT_0110 + IT_0104*IT_0112 + IT_0059*IT_0117 + IT_0096*IT_0122 +
       IT_0147 + IT_0084*IT_0151 + IT_0159 + IT_0086*IT_0165 + IT_0007*IT_0217 +
       IT_0218 + IT_0219 + IT_0220 + IT_0138*IT_0221 + IT_0222 + IT_0223 + conj
      (IT_0058)*IT_0224 + IT_0225 + IT_0227 + IT_0228 + IT_0058*IT_0229 +
       IT_0230 + IT_0231 + IT_0232 + IT_0182*IT_0240 + IT_0091*IT_0242 + IT_0124
      *(IT_0045 + IT_0244) + IT_0267 + IT_0268 + IT_0269 + IT_0270 + IT_0271 +
       IT_0216*IT_0273;
    const ccomplex_t IT_0275 = IT_0034*IT_0079;
    const ccomplex_t IT_0276 = conj(IT_0045)*IT_0098;
    const ccomplex_t IT_0277 = IT_0051*conj(IT_0058);
    const ccomplex_t IT_0278 = conj(IT_0007)*IT_0054;
    const ccomplex_t IT_0279 = conj(IT_0053)*IT_0172;
    const ccomplex_t IT_0280 = conj(IT_0049)*IT_0101;
    const ccomplex_t IT_0281 = IT_0037*IT_0139;
    const ccomplex_t IT_0282 = conj(IT_0007)*IT_0069;
    const ccomplex_t IT_0283 = IT_0034*IT_0169;
    const ccomplex_t IT_0284 = conj(IT_0049)*IT_0125;
    const ccomplex_t IT_0285 = conj(IT_0055)*IT_0224;
    const ccomplex_t IT_0286 = conj(IT_0045)*IT_0132;
    const ccomplex_t IT_0287 = IT_0054*conj(IT_0058);
    const ccomplex_t IT_0288 = IT_0049*IT_0074;
    const ccomplex_t IT_0289 = IT_0037*IT_0148;
    const ccomplex_t IT_0290 = IT_0034*IT_0160;
    const ccomplex_t IT_0291 = IT_0055*IT_0229;
    const ccomplex_t IT_0292 = (-5.33333333333333)*IT_0056;
    const ccomplex_t IT_0293 = conj(IT_0053)*IT_0224;
    const ccomplex_t IT_0294 = conj(IT_0045)*IT_0125;
    const ccomplex_t IT_0295 = conj(IT_0058)*IT_0070;
    const ccomplex_t IT_0296 = conj(IT_0007)*IT_0046 + IT_0045*IT_0074 +
       IT_0058*IT_0081 + IT_0107 + IT_0114 + IT_0143 + IT_0147 + IT_0152 +
       IT_0153 + IT_0156 + IT_0159 + IT_0164 + IT_0168 + IT_0007*IT_0170 +
       IT_0174 + IT_0053*IT_0229 + IT_0293 + IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = IT_0140 + IT_0167;
    const ccomplex_t IT_0298 = 1.33333333333333*IT_0056;
    const ccomplex_t IT_0299 = (-10.6666666666667)*IT_0056;
    const ccomplex_t IT_0300 = IT_0146 + IT_0221;
    const ccomplex_t IT_0301 = conj(IT_0045)*IT_0064;
    const ccomplex_t IT_0302 = IT_0054*conj(IT_0059);
    const ccomplex_t IT_0303 = IT_0057*conj(IT_0068) + conj(IT_0043)*IT_0069 +
       (IT_0046 + IT_0051 + IT_0064)*conj(IT_0078) + IT_0094 + IT_0100 + conj
      (IT_0053)*(IT_0098 + IT_0101) + IT_0037*IT_0115 + IT_0119 + IT_0053*
      (IT_0075 + IT_0117 + IT_0122) + IT_0129 + IT_0007*IT_0141 + IT_0142 +
       IT_0143 + IT_0144 + IT_0153 + IT_0060*IT_0155 + IT_0156 + IT_0157 +
       IT_0113*(IT_0091 + IT_0112 + IT_0158) + IT_0110*(IT_0149 + IT_0154 +
       IT_0158) + IT_0164 + IT_0165*(IT_0161 + IT_0166) + IT_0058*IT_0167 +
       IT_0168 + IT_0171 + conj(IT_0007)*IT_0172 + conj(IT_0010)*(IT_0069 +
       IT_0172) + IT_0173 + IT_0174 + IT_0175 + (IT_0158 + IT_0162)*IT_0182 +
       IT_0055*IT_0188 + conj(IT_0049)*IT_0196 + IT_0049*IT_0204 + (IT_0086 +
       IT_0151 + IT_0198)*IT_0216 + IT_0078*(IT_0123 + IT_0170 + IT_0217) +
       IT_0045*IT_0217 + conj(IT_0096)*(IT_0057 + IT_0224) + conj(IT_0055)*
      (IT_0212 + IT_0224) + IT_0088*(s_13*((-21.3333333333333)*s_12 + (-32)
      *IT_0014) + IT_0240) + IT_0112*IT_0242 + IT_0234*IT_0244 + IT_0096*
      (IT_0140 + IT_0229 + IT_0263) + IT_0068*IT_0263 + IT_0072*(IT_0086 +
       IT_0151 + IT_0273) + IT_0284 + IT_0288 + IT_0289 + IT_0290 + IT_0291 +
       IT_0010*(IT_0141 + IT_0170 + IT_0292) + IT_0043*IT_0292 + IT_0059*IT_0297
       + IT_0104*(IT_0149 + IT_0162 + IT_0298) + IT_0084*(IT_0161 + IT_0221 +
       IT_0299) + IT_0145*(IT_0124 + IT_0300) + IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = conj(IT_0053)*IT_0103;
    const ccomplex_t IT_0305 = conj(IT_0078)*IT_0224;
    const ccomplex_t IT_0306 = conj(IT_0007)*IT_0101;
    const ccomplex_t IT_0307 = conj(IT_0010)*IT_0155;
    const ccomplex_t IT_0308 = conj(IT_0045)*IT_0130;
    const ccomplex_t IT_0309 = conj(IT_0058)*IT_0125;
    const ccomplex_t IT_0310 = IT_0064*conj(IT_0096);
    const ccomplex_t IT_0311 = conj(IT_0059)*IT_0172;
    const ccomplex_t IT_0312 = conj(IT_0043)*IT_0054 + IT_0046*conj(IT_0068) +
       IT_0037*IT_0073 + IT_0034*IT_0116 + conj(IT_0049)*IT_0128 + IT_0110
      *IT_0141 + IT_0084*IT_0167 + IT_0195 + IT_0104*IT_0217 + IT_0304 + IT_0305
       + IT_0306 + IT_0307 + IT_0308 + IT_0309 + IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = IT_0058*IT_0074 + IT_0045*IT_0077 + IT_0072
      *IT_0089 + IT_0007*IT_0122 + IT_0106*IT_0138 + IT_0059*IT_0141 + IT_0104
      *IT_0146 + IT_0084*IT_0158 + IT_0162*IT_0165 + IT_0010*IT_0167 + IT_0184 +
       IT_0185 + IT_0191 + IT_0192 + IT_0195 + IT_0197 + IT_0199 + IT_0200 +
       IT_0203 + IT_0205 + IT_0208 + IT_0209 + IT_0096*IT_0217 + IT_0110*IT_0221
       + IT_0078*IT_0229 + IT_0216*IT_0240 + IT_0241 + IT_0124*(IT_0053 +
       IT_0242) + IT_0182*IT_0273 + IT_0304 + IT_0305 + IT_0306 + IT_0307 +
       IT_0308 + IT_0309 + IT_0310 + IT_0311;
    const ccomplex_t IT_0314 = (-10.6666666666667)*IT_0102;
    const ccomplex_t IT_0315 = conj(IT_0055)*IT_0103;
    const ccomplex_t IT_0316 = conj(IT_0049)*IT_0130;
    const ccomplex_t IT_0317 = conj(IT_0096)*IT_0172;
    const ccomplex_t IT_0318 = conj(IT_0043)*IT_0057;
    const ccomplex_t IT_0319 = conj(IT_0058)*IT_0064 + conj(IT_0059)*(IT_0046 
      + IT_0064) + IT_0049*IT_0077 + IT_0037*IT_0093 + IT_0069*(conj(IT_0068) +
       conj(IT_0096)) + IT_0055*IT_0124 + conj(IT_0045)*IT_0155 + conj(IT_0078)*
      (IT_0054 + IT_0155) + IT_0110*IT_0161 + IT_0149*IT_0165 + IT_0045*IT_0167 
      + IT_0059*IT_0170 + IT_0183 + IT_0184 + IT_0185 + IT_0186 + IT_0187 +
       IT_0189 + IT_0191 + IT_0192 + IT_0193 + IT_0194 + IT_0195 + conj(IT_0053)
      *(IT_0099 + IT_0132 + IT_0196) + IT_0197 + IT_0200 + IT_0201 + IT_0202 +
       IT_0203 + IT_0053*(IT_0080 + IT_0118 + IT_0204) + IT_0205 + IT_0206 +
       IT_0207 + IT_0208 + IT_0209 + IT_0210 + IT_0211 + IT_0213 + IT_0120
      *IT_0217 + conj(IT_0010)*(IT_0057 + IT_0101 + IT_0224) + conj(IT_0007)
      *IT_0224 + IT_0007*IT_0229 + IT_0236 + IT_0240*IT_0244 + IT_0246 + IT_0250
       + IT_0256 + IT_0259 + IT_0262 + IT_0043*IT_0263 + IT_0010*(IT_0229 +
       IT_0263) + IT_0264 + IT_0265 + IT_0266 + IT_0113*IT_0273 + IT_0242*
      (IT_0151 + IT_0273) + IT_0068*IT_0292 + IT_0096*(IT_0141 + IT_0292) +
       IT_0078*IT_0297 + IT_0084*IT_0298 + IT_0104*(IT_0151 + IT_0299) + IT_0088
      *(IT_0198 + IT_0300 + IT_0314) + IT_0182*IT_0314 + IT_0315 + IT_0316 +
       IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = IT_0011*(IT_0036 + IT_0039) + IT_0043*(conj
      (IT_0045)*IT_0046 + conj(IT_0049)*IT_0051 + conj(IT_0053)*IT_0054 + conj
      (IT_0055)*IT_0057) + IT_0060*(IT_0062 + IT_0065 + IT_0067) + IT_0068*
      (IT_0046*conj(IT_0053) + conj(IT_0045)*IT_0054 + conj(IT_0055)*IT_0069 +
       conj(IT_0049)*IT_0070) + IT_0072*IT_0135 + 0.666666666666667*IT_0049
      *IT_0136 + IT_0138*IT_0176 + IT_0034*((-21.3333333333333)*IT_0034 + conj
      (IT_0045)*IT_0063 + IT_0061*conj(IT_0078) + IT_0038*conj(IT_0096) + conj
      (IT_0053)*IT_0097 + IT_0078*IT_0121 + IT_0096*IT_0139 + conj(IT_0055)
      *IT_0177 + conj(IT_0049)*IT_0179) + IT_0037*IT_0180 + IT_0182*IT_0214 +
       IT_0216*IT_0233 + IT_0113*IT_0255 + IT_0096*(IT_0256 + IT_0257 + IT_0258 
      + IT_0259 + IT_0260 + IT_0261) + IT_0104*(IT_0058*IT_0170 + IT_0238 +
       IT_0256 + IT_0257 + IT_0258 + IT_0259 + IT_0260 + IT_0261 + IT_0262 +
       IT_0007*IT_0263 + IT_0264 + IT_0265 + IT_0266) + IT_0045*IT_0272 +
       IT_0244*IT_0274 + IT_0058*(IT_0051*conj(IT_0059) + IT_0046*conj(IT_0096) 
      + IT_0248 + IT_0249 + IT_0251 + IT_0252 + IT_0254 + IT_0275 + IT_0280) +
       IT_0110*(IT_0058*IT_0123 + IT_0007*IT_0140 + IT_0190 + IT_0235 + IT_0237 
      + IT_0239 + IT_0245 + IT_0247 + IT_0248 + IT_0276 + IT_0277 + IT_0278 +
       IT_0279 + IT_0280) + IT_0059*(IT_0248 + IT_0275 + IT_0276 + IT_0277 +
       IT_0278 + IT_0279 + IT_0280) + IT_0010*(conj(IT_0053)*IT_0155 + IT_0281 +
       IT_0282 + IT_0283 + IT_0284 + IT_0285 + IT_0286 + IT_0287) + IT_0084*
      (IT_0058*IT_0140 + conj(IT_0053)*IT_0155 + IT_0163 + IT_0226 + IT_0282 +
       IT_0284 + IT_0285 + IT_0286 + IT_0287 + IT_0288 + IT_0289 + IT_0290 +
       IT_0291 + IT_0007*IT_0292) + IT_0078*(conj(IT_0007)*IT_0046 + conj
      (IT_0049)*IT_0132 + IT_0156 + IT_0293 + IT_0294 + IT_0295) + IT_0165
      *IT_0296 + IT_0007*(conj(IT_0010)*IT_0069 + IT_0046*conj(IT_0078) +
       IT_0057*conj(IT_0096) + conj(IT_0053)*IT_0101 + IT_0281 + IT_0283 +
       IT_0284 + IT_0285 + IT_0301 + IT_0302) + IT_0145*IT_0303 + IT_0053
      *IT_0312 + IT_0242*IT_0313 + IT_0055*(conj(IT_0068)*IT_0069 + conj(IT_0045
      )*IT_0128 + conj(IT_0078)*IT_0155 + IT_0034*IT_0178 + IT_0011*IT_0224 +
       IT_0253 + IT_0315 + IT_0316 + IT_0317 + IT_0318) + IT_0088*IT_0319;
    return create_ccomplex_return(IT_0320);
}

