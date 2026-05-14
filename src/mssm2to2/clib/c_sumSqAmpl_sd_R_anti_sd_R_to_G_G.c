#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sd_R_anti_sd_R_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sd_R_anti_sd_R_to_G_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_sd_R = param->m_sd_R;
    const creal_t Gamma_dr = param->Gamma_dr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_sd_R
      *Gamma_dr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = 0.5*IT_0005;
    const ccomplex_t IT_0007 = -IT_0000;
    const ccomplex_t IT_0008 = IT_0001*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0007, 2);
    const ccomplex_t IT_0010 = IT_0004*(IT_0008 + -IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = 5.33333333333333*s_13;
    const ccomplex_t IT_0013 = s_24*IT_0012;
    const ccomplex_t IT_0014 = IT_0004*(IT_0002 + -IT_0008);
    const ccomplex_t IT_0015 = (-0.5)*IT_0014;
    const ccomplex_t IT_0016 = pow(m_sd_R, 2);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = IT_0004*IT_0008;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = 5.33333333333333*s_24;
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_sd_R
      *Gamma_dr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = (-0.666666666666667)*s_13;
    const ccomplex_t IT_0027 = s_23*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0023;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = s_24*IT_0026;
    const ccomplex_t IT_0031 = s_34*IT_0026;
    const ccomplex_t IT_0032 = g_s*IT_0007;
    const ccomplex_t IT_0033 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = IT_0002*IT_0023;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = 0.666666666666667*s_12;
    const ccomplex_t IT_0039 = s_23*IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0024;
    const ccomplex_t IT_0041 = s_34*IT_0038;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = 0.666666666666667*s_13;
    const ccomplex_t IT_0044 = s_23*IT_0043;
    const ccomplex_t IT_0045 = s_34*IT_0043;
    const ccomplex_t IT_0046 = s_23*IT_0012;
    const ccomplex_t IT_0047 = s_34*IT_0012;
    const ccomplex_t IT_0048 = 5.33333333333333*s_23;
    const ccomplex_t IT_0049 = IT_0016*IT_0048;
    const ccomplex_t IT_0050 = 5.33333333333333*s_34;
    const ccomplex_t IT_0051 = IT_0016*IT_0050;
    const ccomplex_t IT_0052 = -g_s;
    const ccomplex_t IT_0053 = IT_0007*IT_0052;
    const ccomplex_t IT_0054 = IT_0033*(IT_0032 + IT_0053);
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0057 = IT_0016*IT_0056;
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0060 = s_23*IT_0059;
    const ccomplex_t IT_0061 = s_34*IT_0059;
    const ccomplex_t IT_0062 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0063 = IT_0016*IT_0062;
    const ccomplex_t IT_0064 = (-5.33333333333333)*s_34;
    const ccomplex_t IT_0065 = IT_0016*IT_0064;
    const ccomplex_t IT_0066 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0067 = s_23*IT_0066;
    const ccomplex_t IT_0068 = s_34*IT_0066;
    const ccomplex_t IT_0069 = conj(IT_0040)*IT_0068;
    const ccomplex_t IT_0070 = conj(IT_0029)*IT_0031 + IT_0037*IT_0039 +
       IT_0027*conj(IT_0040) + IT_0042 + IT_0040*IT_0044 + IT_0029*IT_0045 +
       conj(IT_0015)*IT_0046 + conj(IT_0011)*IT_0047 + conj(IT_0006)*IT_0049 +
       conj(IT_0019)*IT_0051 + IT_0058 + IT_0015*IT_0060 + IT_0011*IT_0061 +
       IT_0006*IT_0063 + IT_0019*IT_0065 + conj(IT_0037)*IT_0067 + IT_0069;
    const ccomplex_t IT_0071 = IT_0001*IT_0052;
    const ccomplex_t IT_0072 = IT_0033*(IT_0032 + IT_0053 + -IT_0071);
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = s_14*IT_0038;
    const ccomplex_t IT_0075 = pow(s_12, 2);
    const ccomplex_t IT_0076 = 0.666666666666667*IT_0075;
    const ccomplex_t IT_0077 = 0.666666666666667*s_14;
    const ccomplex_t IT_0078 = s_23*IT_0077;
    const ccomplex_t IT_0079 = 1.33333333333333*s_12;
    const ccomplex_t IT_0080 = s_13*IT_0079;
    const ccomplex_t IT_0081 = IT_0033*(IT_0053 + -IT_0071);
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = s_24*IT_0079;
    const ccomplex_t IT_0084 = g_s*IT_0001;
    const ccomplex_t IT_0085 = IT_0033*(IT_0032 + IT_0053 + -IT_0071 + 
      -IT_0084);
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0016*IT_0079;
    const ccomplex_t IT_0088 = pow(g_s, 2);
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = s_13 + s_23;
    const ccomplex_t IT_0091 = IT_0053*IT_0090;
    const ccomplex_t IT_0092 = s_14 + s_24;
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = g_s*IT_0093;
    const ccomplex_t IT_0095 = IT_0091 + IT_0094;
    const ccomplex_t IT_0096 = g_s*s_24;
    const ccomplex_t IT_0097 = IT_0001*IT_0096;
    const ccomplex_t IT_0098 = s_23*IT_0071;
    const ccomplex_t IT_0099 = -IT_0097 + -IT_0098;
    const ccomplex_t IT_0100 = IT_0095 + IT_0099;
    const ccomplex_t IT_0101 = IT_0033*IT_0100;
    const ccomplex_t IT_0102 = IT_0089 + (-0.5)*IT_0101;
    const ccomplex_t IT_0103 = IT_0079*IT_0102;
    const ccomplex_t IT_0104 = IT_0033*(IT_0032 + IT_0053 + -IT_0084);
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = 1.33333333333333*IT_0075;
    const ccomplex_t IT_0107 = IT_0033*(IT_0032 + -IT_0084);
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = 1.33333333333333*s_13;
    const ccomplex_t IT_0110 = IT_0016*IT_0109;
    const ccomplex_t IT_0111 = IT_0033*IT_0053;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = 1.33333333333333*s_24;
    const ccomplex_t IT_0114 = IT_0016*IT_0113;
    const ccomplex_t IT_0115 = IT_0033*IT_0071;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = pow(m_sd_R, 4);
    const ccomplex_t IT_0118 = 1.33333333333333*IT_0117;
    const ccomplex_t IT_0119 = 5.33333333333333*s_12;
    const ccomplex_t IT_0120 = IT_0016*IT_0119;
    const ccomplex_t IT_0121 = conj(IT_0037)*IT_0120;
    const ccomplex_t IT_0122 = conj(IT_0025)*IT_0046;
    const ccomplex_t IT_0123 = 0.5*IT_0036;
    const ccomplex_t IT_0124 = 5.33333333333333*s_14;
    const ccomplex_t IT_0125 = s_23*IT_0124;
    const ccomplex_t IT_0126 = s_23*IT_0119;
    const ccomplex_t IT_0127 = IT_0016*IT_0124;
    const ccomplex_t IT_0128 = 0.5*IT_0034;
    const ccomplex_t IT_0129 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0130 = s_13*IT_0129;
    const ccomplex_t IT_0131 = 0.5*IT_0081;
    const ccomplex_t IT_0132 = s_24*IT_0129;
    const ccomplex_t IT_0133 = IT_0089 + 0.5*IT_0101;
    const ccomplex_t IT_0134 = IT_0129*IT_0133;
    const ccomplex_t IT_0135 = 0.5*IT_0107;
    const ccomplex_t IT_0136 = 0.5*IT_0111;
    const ccomplex_t IT_0137 = (-10.6666666666667)*s_24;
    const ccomplex_t IT_0138 = IT_0016*IT_0137;
    const ccomplex_t IT_0139 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0140 = IT_0016*IT_0139;
    const ccomplex_t IT_0141 = IT_0037*IT_0140;
    const ccomplex_t IT_0142 = IT_0025*IT_0060;
    const ccomplex_t IT_0143 = IT_0016*IT_0059;
    const ccomplex_t IT_0144 = (-5.33333333333333)*s_14;
    const ccomplex_t IT_0145 = s_23*IT_0144;
    const ccomplex_t IT_0146 = s_23*IT_0139;
    const ccomplex_t IT_0147 = IT_0016*IT_0144;
    const ccomplex_t IT_0148 = (-5.33333333333333)*IT_0117;
    const ccomplex_t IT_0149 = s_14*IT_0066;
    const ccomplex_t IT_0150 = (-0.666666666666667)*IT_0075;
    const ccomplex_t IT_0151 = (-0.666666666666667)*s_14;
    const ccomplex_t IT_0152 = s_23*IT_0151;
    const ccomplex_t IT_0153 = IT_0015*IT_0039 + conj(IT_0015)*IT_0067 +
       IT_0019*IT_0074 + IT_0006*IT_0076 + IT_0011*IT_0078 + IT_0035*IT_0080 +
       IT_0082*IT_0083 + IT_0086*IT_0087 + IT_0103 + IT_0105*IT_0106 + IT_0108
      *IT_0110 + IT_0112*IT_0114 + IT_0116*IT_0118 + IT_0121 + IT_0122 + IT_0017
      *conj(IT_0123) + conj(IT_0029)*IT_0125 + conj(IT_0040)*(IT_0126 + IT_0127)
       + IT_0128*IT_0130 + IT_0131*IT_0132 + IT_0134 + IT_0057*IT_0135 + IT_0136
      *IT_0138 + IT_0141 + IT_0142 + IT_0123*IT_0143 + IT_0029*IT_0145 + IT_0040
      *(IT_0146 + IT_0147) + IT_0073*IT_0148 + conj(IT_0019)*IT_0149 + conj
      (IT_0006)*IT_0150 + conj(IT_0011)*IT_0152;
    const ccomplex_t IT_0154 = IT_0039*IT_0040;
    const ccomplex_t IT_0155 = IT_0016*IT_0038;
    const ccomplex_t IT_0156 = IT_0037*IT_0155;
    const ccomplex_t IT_0157 = IT_0025*IT_0044;
    const ccomplex_t IT_0158 = IT_0016*IT_0043;
    const ccomplex_t IT_0159 = IT_0123*IT_0158;
    const ccomplex_t IT_0160 = IT_0016*IT_0077;
    const ccomplex_t IT_0161 = IT_0079*IT_0133;
    const ccomplex_t IT_0162 = s_14*IT_0119;
    const ccomplex_t IT_0163 = 5.33333333333333*IT_0075;
    const ccomplex_t IT_0164 = IT_0102*IT_0129;
    const ccomplex_t IT_0165 = s_14*IT_0139;
    const ccomplex_t IT_0166 = (-5.33333333333333)*IT_0075;
    const ccomplex_t IT_0167 = conj(IT_0040)*IT_0067;
    const ccomplex_t IT_0168 = IT_0016*IT_0066;
    const ccomplex_t IT_0169 = conj(IT_0037)*IT_0168;
    const ccomplex_t IT_0170 = conj(IT_0025)*IT_0027;
    const ccomplex_t IT_0171 = IT_0016*IT_0026;
    const ccomplex_t IT_0172 = conj(IT_0123)*IT_0171;
    const ccomplex_t IT_0173 = IT_0016*IT_0151;
    const ccomplex_t IT_0174 = IT_0029*IT_0078 + IT_0057*IT_0108 + conj
      (IT_0011)*IT_0125 + conj(IT_0015)*IT_0126 + IT_0080*IT_0128 + IT_0035
      *IT_0130 + IT_0083*IT_0131 + IT_0082*IT_0132 + IT_0110*IT_0135 + IT_0114
      *IT_0136 + IT_0112*IT_0138 + IT_0011*IT_0145 + IT_0015*IT_0146 + IT_0116
      *IT_0148 + conj(IT_0029)*IT_0152 + IT_0154 + IT_0156 + IT_0157 + IT_0159 +
       IT_0040*IT_0160 + IT_0161 + conj(IT_0019)*IT_0162 + conj(IT_0006)*IT_0163
       + IT_0164 + IT_0019*IT_0165 + IT_0006*IT_0166 + IT_0167 + IT_0169 +
       IT_0170 + IT_0172 + conj(IT_0040)*IT_0173;
    const ccomplex_t IT_0175 = pow(s_23, 2);
    const ccomplex_t IT_0176 = 5.33333333333333*IT_0175;
    const ccomplex_t IT_0177 = IT_0051 + IT_0176;
    const ccomplex_t IT_0178 = s_24*IT_0066;
    const ccomplex_t IT_0179 = IT_0171 + IT_0178;
    const ccomplex_t IT_0180 = s_34*IT_0048;
    const ccomplex_t IT_0181 = conj(IT_0037)*IT_0049;
    const ccomplex_t IT_0182 = IT_0062*IT_0133;
    const ccomplex_t IT_0183 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0184 = s_24*IT_0183;
    const ccomplex_t IT_0185 = s_24*IT_0048;
    const ccomplex_t IT_0186 = IT_0016*IT_0183;
    const ccomplex_t IT_0187 = IT_0048*IT_0133;
    const ccomplex_t IT_0188 = IT_0102*IT_0183;
    const ccomplex_t IT_0189 = s_34*IT_0183;
    const ccomplex_t IT_0190 = (-0.666666666666667)*IT_0175;
    const ccomplex_t IT_0191 = 0.666666666666667*s_23;
    const ccomplex_t IT_0192 = IT_0016*IT_0191;
    const ccomplex_t IT_0193 = 0.666666666666667*s_34;
    const ccomplex_t IT_0194 = IT_0016*IT_0193;
    const ccomplex_t IT_0195 = s_34*IT_0119;
    const ccomplex_t IT_0196 = conj(IT_0040)*IT_0195;
    const ccomplex_t IT_0197 = s_34*IT_0139;
    const ccomplex_t IT_0198 = IT_0040*IT_0197;
    const ccomplex_t IT_0199 = (-0.666666666666667)*s_34;
    const ccomplex_t IT_0200 = IT_0016*IT_0199;
    const ccomplex_t IT_0201 = conj(IT_0015)*IT_0027 + conj(IT_0011)*IT_0031 +
       IT_0015*IT_0044 + IT_0011*IT_0045 + conj(IT_0040)*IT_0046 + conj(IT_0029)
      *IT_0047 + IT_0040*IT_0060 + IT_0029*IT_0061 + conj(IT_0037)*IT_0126 +
       IT_0037*IT_0146 + conj(IT_0006)*IT_0186 + IT_0006*IT_0192 + IT_0019
      *IT_0194 + IT_0196 + IT_0198 + conj(IT_0019)*IT_0200;
    const ccomplex_t IT_0202 = s_34*IT_0191;
    const ccomplex_t IT_0203 = 0.666666666666667*IT_0175;
    const ccomplex_t IT_0204 = IT_0040*IT_0194;
    const ccomplex_t IT_0205 = s_34*IT_0062;
    const ccomplex_t IT_0206 = (-5.33333333333333)*IT_0175;
    const ccomplex_t IT_0207 = conj(IT_0040)*IT_0200;
    const ccomplex_t IT_0208 = IT_0058 + conj(IT_0006)*IT_0126 + IT_0006
      *IT_0146 + conj(IT_0015)*IT_0176 + conj(IT_0011)*IT_0180 + conj(IT_0037)
      *IT_0186 + conj(IT_0029)*IT_0189 + conj(IT_0040)*IT_0190 + IT_0037*IT_0192
       + conj(IT_0019)*IT_0195 + IT_0019*IT_0197 + IT_0029*IT_0202 + IT_0040
      *IT_0203 + IT_0204 + IT_0011*IT_0205 + IT_0015*IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = s_24*IT_0191;
    const ccomplex_t IT_0210 = 0.666666666666667*s_24;
    const ccomplex_t IT_0211 = IT_0016*IT_0210;
    const ccomplex_t IT_0212 = s_24*IT_0124;
    const ccomplex_t IT_0213 = 1.33333333333333*s_14;
    const ccomplex_t IT_0214 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0215 = -conj(IT_0040);
    const ccomplex_t IT_0216 = IT_0040 + IT_0215;
    const ccomplex_t IT_0217 = 1.33333333333333*s_23;
    const ccomplex_t IT_0218 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0219 = s_24*IT_0144;
    const ccomplex_t IT_0220 = (-5.33333333333333)*s_24;
    const ccomplex_t IT_0221 = (-0.666666666666667)*s_24;
    const ccomplex_t IT_0222 = IT_0016*IT_0221;
    const ccomplex_t IT_0223 = conj(IT_0015)*IT_0051 + IT_0015*IT_0065 + conj
      (IT_0006)*IT_0127 + IT_0113*IT_0133 + IT_0102*IT_0137 + IT_0006*IT_0147 +
       conj(IT_0123)*IT_0184 + conj(IT_0025)*IT_0189 + IT_0025*IT_0202 + IT_0204
       + IT_0207 + IT_0123*IT_0209 + IT_0037*IT_0211 + conj(IT_0019)*IT_0212 +
       s_23*(IT_0128*IT_0213 + IT_0035*IT_0214) + 0.666666666666667*s_24*(s_24
      *IT_0216 + 1.5*IT_0135*IT_0217 + 1.5*IT_0108*IT_0218) + IT_0019*IT_0219 +
       s_34*(conj(IT_0011)*IT_0020 + IT_0029*IT_0210 + IT_0011*IT_0220 + conj
      (IT_0029)*IT_0221) + conj(IT_0037)*IT_0222;
    const ccomplex_t IT_0224 = (-0.5)*IT_0104;
    const ccomplex_t IT_0225 = s_24*IT_0043;
    const ccomplex_t IT_0226 = 0.666666666666667*IT_0117;
    const ccomplex_t IT_0227 = s_14*IT_0079;
    const ccomplex_t IT_0228 = s_23*IT_0079;
    const ccomplex_t IT_0229 = IT_0055*IT_0106;
    const ccomplex_t IT_0230 = IT_0016*IT_0213;
    const ccomplex_t IT_0231 = IT_0016*IT_0217;
    const ccomplex_t IT_0232 = s_24*IT_0119;
    const ccomplex_t IT_0233 = s_14*IT_0129;
    const ccomplex_t IT_0234 = s_23*IT_0129;
    const ccomplex_t IT_0235 = (-0.5)*IT_0054;
    const ccomplex_t IT_0236 = IT_0016*IT_0129;
    const ccomplex_t IT_0237 = (-10.6666666666667)*IT_0075;
    const ccomplex_t IT_0238 = IT_0016*IT_0214;
    const ccomplex_t IT_0239 = IT_0016*IT_0218;
    const ccomplex_t IT_0240 = s_24*IT_0059;
    const ccomplex_t IT_0241 = s_24*IT_0139;
    const ccomplex_t IT_0242 = (-0.666666666666667)*IT_0117;
    const ccomplex_t IT_0243 = IT_0013*conj(IT_0029) + conj(IT_0011)*IT_0030 +
       IT_0086*IT_0087 + IT_0103 + IT_0106*IT_0116 + IT_0105*IT_0118 + IT_0121 +
       IT_0122 + conj(IT_0123)*IT_0126 + IT_0134 + IT_0141 + IT_0142 + IT_0123
      *IT_0146 + IT_0015*IT_0158 + conj(IT_0015)*IT_0171 + IT_0019*IT_0211 +
       conj(IT_0019)*IT_0222 + IT_0148*IT_0224 + IT_0011*IT_0225 + IT_0006
      *IT_0226 + IT_0112*IT_0227 + IT_0108*IT_0228 + IT_0229 + IT_0082*IT_0230 +
       IT_0035*IT_0231 + conj(IT_0040)*(IT_0017 + IT_0232) + IT_0136*IT_0233 +
       IT_0135*IT_0234 + IT_0235*IT_0236 + IT_0073*IT_0237 + IT_0131*IT_0238 +
       IT_0128*IT_0239 + IT_0029*IT_0240 + IT_0040*(IT_0143 + IT_0241) + conj
      (IT_0006)*IT_0242;
    const ccomplex_t IT_0244 = IT_0040 + IT_0123;
    const ccomplex_t IT_0245 = 1.33333333333333*IT_0016;
    const ccomplex_t IT_0246 = IT_0133*IT_0245;
    const ccomplex_t IT_0247 = (-10.6666666666667)*IT_0016;
    const ccomplex_t IT_0248 = s_24*IT_0062;
    const ccomplex_t IT_0249 = conj(IT_0040) + conj(IT_0123);
    const ccomplex_t IT_0250 = conj(IT_0015)*IT_0049 + IT_0015*IT_0063 + conj
      (IT_0006)*IT_0120 + IT_0114*IT_0131 + IT_0112*IT_0132 + IT_0083*IT_0136 +
       IT_0082*IT_0138 + IT_0006*IT_0140 + IT_0086*IT_0148 + conj(IT_0029)
      *IT_0184 + conj(IT_0011)*IT_0185 + conj(IT_0025)*IT_0190 + IT_0025*IT_0203
       + IT_0029*IT_0209 + IT_0040*IT_0211 + conj(IT_0040)*IT_0222 + IT_0037
      *IT_0226 + IT_0128*IT_0228 + IT_0135*IT_0231 + conj(IT_0019)*IT_0232 +
       IT_0035*IT_0234 + IT_0116*IT_0236 + IT_0108*IT_0239 + IT_0019*IT_0241 +
       conj(IT_0037)*IT_0242 + IT_0192*IT_0244 + IT_0246 + IT_0102*IT_0247 +
       IT_0011*IT_0248 + IT_0186*IT_0249;
    const ccomplex_t IT_0251 = s_34*IT_0077;
    const ccomplex_t IT_0252 = s_24*IT_0109;
    const ccomplex_t IT_0253 = s_34*IT_0124;
    const ccomplex_t IT_0254 = IT_0019 + -conj(IT_0019);
    const ccomplex_t IT_0255 = s_24*IT_0056;
    const ccomplex_t IT_0256 = s_34*IT_0144;
    const ccomplex_t IT_0257 = s_34*IT_0151;
    const ccomplex_t IT_0258 = 1.5*IT_0006;
    const ccomplex_t IT_0259 = IT_0015*IT_0041 + conj(IT_0025)*IT_0047 +
       IT_0025*IT_0061 + conj(IT_0015)*IT_0068 + IT_0013*conj(IT_0123) + conj
      (IT_0006)*IT_0149 + IT_0196 + IT_0198 + conj(IT_0040)*IT_0212 + IT_0102
      *IT_0213 + IT_0133*IT_0214 + IT_0040*IT_0219 + conj(IT_0037)*IT_0232 +
       IT_0123*IT_0240 + IT_0037*IT_0241 + IT_0011*IT_0251 + IT_0108*IT_0252 +
       conj(IT_0029)*IT_0253 + 0.666666666666667*s_14*(1.5*IT_0035*IT_0109 + 1.5
      *IT_0056*IT_0128 + s_14*IT_0254) + IT_0135*IT_0255 + IT_0029*IT_0256 +
       conj(IT_0011)*IT_0257 + 0.666666666666667*IT_0074*IT_0258;
    const ccomplex_t IT_0260 = conj(IT_0040)*IT_0051;
    const ccomplex_t IT_0261 = IT_0040*IT_0065;
    const ccomplex_t IT_0262 = (-0.1875)*IT_0102;
    const ccomplex_t IT_0263 = (-0.1875)*IT_0006;
    const ccomplex_t IT_0264 = s_24*IT_0077;
    const ccomplex_t IT_0265 = (-0.1875)*IT_0019;
    const ccomplex_t IT_0266 = s_24*IT_0151;
    const ccomplex_t IT_0267 = (-0.1875)*conj(IT_0019);
    const ccomplex_t IT_0268 = (-0.1875)*IT_0015;
    const ccomplex_t IT_0269 = (-0.1875)*conj(IT_0015);
    const ccomplex_t IT_0270 = IT_0016*IT_0220;
    const ccomplex_t IT_0271 = (-0.1875)*IT_0037;
    const ccomplex_t IT_0272 = (-0.1875)*conj(IT_0037);
    const ccomplex_t IT_0273 = (-0.1875)*IT_0108;
    const ccomplex_t IT_0274 = (-0.1875)*IT_0135;
    const ccomplex_t IT_0275 = IT_0133*IT_0137 + conj(IT_0006)*IT_0173 + conj
      (IT_0025)*IT_0180 + conj(IT_0123)*IT_0185 + IT_0025*IT_0205 + s_23*
      (IT_0035*IT_0213 + IT_0128*IT_0214) + s_34*(IT_0020*conj(IT_0029) +
       IT_0011*IT_0210 + IT_0029*IT_0220 + conj(IT_0011)*IT_0221) + IT_0123
      *IT_0248 + IT_0260 + IT_0261 + (-5.33333333333333)*IT_0113*IT_0262 + (
      -5.33333333333333)*IT_0160*IT_0263 + (-5.33333333333333)*IT_0264*IT_0265 +
       (-5.33333333333333)*IT_0266*IT_0267 + (-5.33333333333333)*IT_0194*IT_0268
       + (-5.33333333333333)*IT_0200*IT_0269 + (-5.33333333333333)*IT_0270
      *IT_0271 + (-5.33333333333333)*IT_0021*IT_0272 + (-5.33333333333333)*s_24*
      (s_24*IT_0216 + IT_0217*IT_0273 + IT_0218*IT_0274);
    const ccomplex_t IT_0276 = s_13*IT_0038;
    const ccomplex_t IT_0277 = (-0.5)*IT_0085;
    const ccomplex_t IT_0278 = s_13*IT_0119;
    const ccomplex_t IT_0279 = (-10.6666666666667)*IT_0117;
    const ccomplex_t IT_0280 = s_13*IT_0139;
    const ccomplex_t IT_0281 = s_14*IT_0059;
    const ccomplex_t IT_0282 = s_13*IT_0066;
    const ccomplex_t IT_0283 = pow(s_13, 2);
    const ccomplex_t IT_0284 = (IT_0025 + -conj(IT_0025))*IT_0283;
    const ccomplex_t IT_0285 = s_14*(conj(IT_0011)*IT_0012 + IT_0026*conj
      (IT_0029) + IT_0029*IT_0043) + IT_0037*IT_0076 + IT_0029*IT_0078 + conj
      (IT_0011)*IT_0125 + (IT_0035 + IT_0108)*IT_0130 + (IT_0080 + IT_0110)*
      (IT_0128 + IT_0135) + conj(IT_0037)*IT_0150 + conj(IT_0029)*IT_0152 +
       IT_0154 + IT_0156 + IT_0157 + IT_0159 + IT_0161 + conj(IT_0019)*(IT_0127 
      + IT_0162) + conj(IT_0006)*(IT_0120 + IT_0163) + IT_0164 + IT_0019*
      (IT_0147 + IT_0165) + IT_0167 + IT_0169 + IT_0170 + IT_0172 + IT_0131*
      (IT_0083 + IT_0227) + IT_0136*(IT_0114 + IT_0230) + IT_0082*(IT_0132 +
       IT_0233) + IT_0118*(IT_0073 + IT_0235) + (IT_0055 + IT_0086 + IT_0105 +
       IT_0116)*IT_0236 + (IT_0086 + IT_0105)*IT_0237 + IT_0112*(IT_0138 +
       IT_0238) + IT_0246 + IT_0102*IT_0247 + 0.666666666666667*(IT_0140 +
       IT_0166)*IT_0258 + IT_0040*(IT_0074 + IT_0160 + IT_0276) + IT_0123
      *IT_0276 + IT_0087*(IT_0073 + IT_0224 + IT_0235 + IT_0277) + conj(IT_0015)
      *(IT_0126 + IT_0278) + (IT_0055 + IT_0116)*IT_0279 + IT_0015*(IT_0146 +
       IT_0280) + IT_0011*(IT_0145 + IT_0281) + conj(IT_0040)*(IT_0149 + IT_0173
       + IT_0282) + conj(IT_0123)*IT_0282 + 0.666666666666667*IT_0284;
    const ccomplex_t IT_0286 = IT_0133*IT_0247;
    const ccomplex_t IT_0287 = conj(IT_0006)*IT_0168;
    const ccomplex_t IT_0288 = (-0.1875)*IT_0128;
    const ccomplex_t IT_0289 = (-0.1875)*IT_0035;
    const ccomplex_t IT_0290 = s_14*(conj(IT_0011)*IT_0026 + IT_0012*conj
      (IT_0029) + IT_0011*IT_0043) + IT_0086*IT_0106 + IT_0087*(IT_0105 +
       IT_0116) + conj(IT_0040)*IT_0162 + IT_0040*IT_0165 + IT_0082*IT_0227 +
       IT_0112*IT_0230 + IT_0131*IT_0233 + IT_0148*IT_0235 + IT_0073*IT_0236 +
       IT_0136*IT_0238 + (-5.33333333333333)*IT_0245*IT_0262 + (
      -5.33333333333333)*IT_0155*IT_0263 + (-5.33333333333333)*IT_0160*IT_0265 +
       (-5.33333333333333)*IT_0173*IT_0267 + (-5.33333333333333)*IT_0166*IT_0271
       + (-5.33333333333333)*IT_0163*IT_0272 + (-5.33333333333333)*IT_0080
      *IT_0273 + (-5.33333333333333)*IT_0130*IT_0274 + (-5.33333333333333)
      *IT_0268*IT_0276 + IT_0249*IT_0278 + IT_0244*IT_0280 + IT_0029*IT_0281 + (
      -5.33333333333333)*IT_0269*IT_0282 + (-5.33333333333333)*IT_0284 + IT_0286
       + IT_0287 + (-5.33333333333333)*IT_0057*IT_0288 + (-5.33333333333333)
      *IT_0110*IT_0289;
    const ccomplex_t IT_0291 = m_sd_R*(IT_0037 + -conj(IT_0037));
    const ccomplex_t IT_0292 = conj(IT_0015)*IT_0048 + conj(IT_0011)*IT_0050 +
       IT_0035*IT_0056 + IT_0015*IT_0062 + IT_0011*IT_0064 + (-21.3333333333333)
      *IT_0102 + conj(IT_0006)*IT_0119 + conj(IT_0019)*IT_0124 + IT_0109*IT_0128
       + 5.33333333333333*IT_0133 + IT_0006*IT_0139 + IT_0019*IT_0144 + IT_0029
      *IT_0193 + conj(IT_0029)*IT_0199 + IT_0040*IT_0210 + IT_0135*IT_0217 +
       IT_0108*IT_0218 + conj(IT_0040)*IT_0221 + IT_0191*IT_0244 +
       0.666666666666667*m_sd_R*IT_0291;
    const ccomplex_t IT_0293 = IT_0006*IT_0038 + IT_0020*conj(IT_0040) + conj
      (IT_0029)*IT_0050 + IT_0029*IT_0064 + conj(IT_0006)*IT_0066 + IT_0019
      *IT_0077 + IT_0035*IT_0109 + IT_0056*IT_0128 + (-21.3333333333333)*IT_0133
       + conj(IT_0019)*IT_0151 + conj(IT_0015)*IT_0183 + IT_0015*IT_0191 +
       IT_0011*IT_0193 + conj(IT_0011)*IT_0199 + IT_0108*IT_0217 + IT_0135
      *IT_0218 + IT_0040*IT_0220 + (-5.33333333333333)*m_sd_R*IT_0291;
    const ccomplex_t IT_0294 = 5.33333333333333*IT_0117;
    const ccomplex_t IT_0295 = conj(IT_0015)*IT_0017;
    const ccomplex_t IT_0296 = IT_0067*conj(IT_0123);
    const ccomplex_t IT_0297 = s_24*IT_0038;
    const ccomplex_t IT_0298 = conj(IT_0011)*IT_0013 + conj(IT_0019)*IT_0021 +
       conj(IT_0029)*IT_0030 + IT_0039*IT_0123 + IT_0015*IT_0143 + (IT_0006 +
       IT_0105)*IT_0148 + IT_0156 + IT_0157 + IT_0161 + IT_0164 + IT_0169 +
       IT_0170 + conj(IT_0040)*IT_0179 + IT_0029*IT_0225 + IT_0136*IT_0227 +
       IT_0135*IT_0228 + IT_0131*IT_0230 + IT_0128*IT_0231 + IT_0112*IT_0233 +
       IT_0108*IT_0234 + IT_0086*IT_0236 + IT_0116*IT_0237 + IT_0082*IT_0238 +
       IT_0035*IT_0239 + IT_0011*IT_0240 + IT_0019*IT_0270 + conj(IT_0006)
      *IT_0294 + IT_0295 + IT_0296 + IT_0040*(IT_0158 + IT_0297);
    const ccomplex_t IT_0299 = conj(IT_0025)*IT_0031;
    const ccomplex_t IT_0300 = IT_0042 + IT_0025*IT_0045 + IT_0069 + IT_0030
      *conj(IT_0123) + conj(IT_0006)*IT_0162 + IT_0006*IT_0165 + conj(IT_0037)
      *IT_0178 + conj(IT_0015)*IT_0195 + IT_0015*IT_0197 + IT_0133*IT_0213 +
       IT_0102*IT_0214 + IT_0123*IT_0225 + IT_0029*IT_0251 + IT_0135*IT_0252 +
       conj(IT_0011)*IT_0253 + IT_0108*IT_0255 + IT_0011*IT_0256 + conj(IT_0029)
      *IT_0257 + IT_0040*IT_0264 + conj(IT_0040)*IT_0266 + (-5.33333333333333)
      *s_14*(s_14*IT_0254 + IT_0109*IT_0288 + IT_0056*IT_0289) + IT_0037*IT_0297
       + IT_0299;
    const ccomplex_t IT_0301 = conj(IT_0015)*IT_0186;
    const ccomplex_t IT_0302 = conj(IT_0019)*IT_0178;
    const ccomplex_t IT_0303 = conj(IT_0011)*IT_0184;
    const ccomplex_t IT_0304 = conj(IT_0037)*IT_0294;
    const ccomplex_t IT_0305 = conj(IT_0025)*IT_0176;
    const ccomplex_t IT_0306 = conj(IT_0029)*IT_0185;
    const ccomplex_t IT_0307 = IT_0021*conj(IT_0040) + IT_0083*IT_0112 +
       IT_0082*IT_0114 + IT_0087*(IT_0105 + IT_0116) + IT_0086*IT_0118 + IT_0132
      *IT_0136 + IT_0131*IT_0138 + IT_0006*IT_0155 + IT_0015*IT_0192 + IT_0025
      *IT_0206 + IT_0011*IT_0209 + IT_0035*IT_0228 + IT_0229 + IT_0108*IT_0231 +
       IT_0128*IT_0234 + (IT_0073 + IT_0224)*IT_0236 + IT_0235*IT_0237 + IT_0135
      *IT_0239 + IT_0063*IT_0244 + IT_0102*IT_0245 + IT_0029*IT_0248 + IT_0049
      *IT_0249 + IT_0040*IT_0270 + IT_0148*(IT_0037 + IT_0277) + IT_0286 +
       IT_0287 + IT_0019*IT_0297 + IT_0301 + IT_0302 + IT_0303 + IT_0304 +
       IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = conj(IT_0019)*IT_0068;
    const ccomplex_t IT_0309 = conj(IT_0011)*IT_0189;
    const ccomplex_t IT_0310 = conj(IT_0015)*IT_0190;
    const ccomplex_t IT_0311 = conj(IT_0006)*IT_0067;
    const ccomplex_t IT_0312 = conj(IT_0029)*IT_0180;
    const ccomplex_t IT_0313 = IT_0006*IT_0039 + IT_0019*IT_0041 + IT_0037
      *IT_0063 + conj(IT_0040)*IT_0176 + IT_0181 + IT_0011*IT_0202 + IT_0015
      *IT_0203 + IT_0029*IT_0205 + IT_0040*IT_0206 + IT_0260 + IT_0261 + IT_0308
       + IT_0309 + IT_0310 + IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = conj(IT_0006)*(IT_0011*IT_0013 + IT_0015
      *IT_0017 + IT_0019*IT_0021 + IT_0025*IT_0027 + IT_0029*IT_0030) + conj
      (IT_0019)*(IT_0013*IT_0015 + IT_0025*IT_0031) + IT_0035*IT_0070 + IT_0073
      *IT_0153 + IT_0116*IT_0174 + IT_0040*(conj(IT_0019)*IT_0030 + IT_0021*conj
      (IT_0037) + conj(IT_0123)*IT_0177 + conj(IT_0006)*IT_0179 + conj(IT_0025)
      *IT_0180 + IT_0181 + IT_0182 + conj(IT_0015)*IT_0184) + conj(IT_0037)*
      (IT_0025*IT_0176 + IT_0019*IT_0178 + IT_0011*IT_0184 + IT_0029*IT_0185 +
       IT_0015*IT_0186) + conj(IT_0040)*(IT_0019*IT_0030 + 10.6666666666667*s_23
      *s_24*IT_0040 + IT_0037*(IT_0021 + IT_0049) + IT_0123*IT_0177 + IT_0006
      *IT_0179 + IT_0025*IT_0180 + IT_0015*IT_0184 + IT_0187 + IT_0188) + conj
      (IT_0123)*(IT_0019*IT_0068 + IT_0029*IT_0180 + IT_0187 + IT_0188 + IT_0011
      *IT_0189 + IT_0015*IT_0190) + IT_0128*IT_0201 + IT_0108*IT_0208 + IT_0082
      *IT_0223 + IT_0224*IT_0243 + IT_0086*IT_0250 + IT_0136*IT_0259 + IT_0131
      *IT_0275 + IT_0055*IT_0285 + IT_0235*IT_0290 + IT_0102*IT_0292 + IT_0133
      *IT_0293 + IT_0006*(conj(IT_0011)*IT_0013 + conj(IT_0019)*IT_0021 + conj
      (IT_0025)*IT_0027 + conj(IT_0029)*IT_0030 + conj(IT_0037)*IT_0168 + conj
      (IT_0006)*IT_0294 + IT_0295 + IT_0296) + IT_0105*IT_0298 + IT_0019*
      (IT_0013*conj(IT_0015) + IT_0299) + IT_0112*IT_0300 + IT_0037*(IT_0049
      *conj(IT_0123) + IT_0287 + IT_0301 + IT_0302 + IT_0303 + IT_0304 + IT_0305
       + IT_0306) + IT_0277*IT_0307 + IT_0123*(IT_0181 + IT_0182 + IT_0308 +
       IT_0309 + IT_0310 + IT_0311 + IT_0312) + IT_0135*IT_0313;
    return create_ccomplex_return(IT_0314);
}

