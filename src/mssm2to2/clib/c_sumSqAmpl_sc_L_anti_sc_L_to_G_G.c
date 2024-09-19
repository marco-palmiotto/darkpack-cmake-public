#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_sc_L_anti_sc_L_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_sc_L_anti_sc_L_to_G_G(
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t Gamma_cl = param->Gamma_cl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_sc_L
      *Gamma_cl + (-0.5)*reg_prop, -1);
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
    const ccomplex_t IT_0016 = pow(m_sc_L, 2);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = IT_0004*IT_0008;
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = 5.33333333333333*s_24;
    const ccomplex_t IT_0021 = IT_0016*IT_0020;
    const ccomplex_t IT_0022 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_sc_L
      *Gamma_cl + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = (-0.5)*IT_0024;
    const ccomplex_t IT_0026 = (-0.666666666666667)*s_13;
    const ccomplex_t IT_0027 = s_23*IT_0026;
    const ccomplex_t IT_0028 = IT_0009*IT_0023;
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = s_24*IT_0026;
    const ccomplex_t IT_0031 = s_34*IT_0026;
    const ccomplex_t IT_0032 = IT_0002*IT_0023;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0035 = s_24*IT_0034;
    const ccomplex_t IT_0036 = conj(IT_0033)*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0032;
    const ccomplex_t IT_0038 = s_34*IT_0034;
    const ccomplex_t IT_0039 = conj(IT_0025)*IT_0031;
    const ccomplex_t IT_0040 = 5.33333333333333*s_23;
    const ccomplex_t IT_0041 = pow(s_23, 2);
    const ccomplex_t IT_0042 = 5.33333333333333*IT_0041;
    const ccomplex_t IT_0043 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0044 = s_24*IT_0043;
    const ccomplex_t IT_0045 = s_34*IT_0043;
    const ccomplex_t IT_0046 = g_s*IT_0007;
    const ccomplex_t IT_0047 = g_s*IT_0001;
    const ccomplex_t IT_0048 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0049 = (IT_0046 + -IT_0047)*IT_0048;
    const ccomplex_t IT_0050 = (-0.5)*IT_0049;
    const ccomplex_t IT_0051 = 0.666666666666667*s_23;
    const ccomplex_t IT_0052 = s_34*IT_0051;
    const ccomplex_t IT_0053 = IT_0016*IT_0051;
    const ccomplex_t IT_0054 = 5.33333333333333*s_12;
    const ccomplex_t IT_0055 = s_23*IT_0054;
    const ccomplex_t IT_0056 = s_34*IT_0054;
    const ccomplex_t IT_0057 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0058 = s_23*IT_0057;
    const ccomplex_t IT_0059 = s_34*IT_0057;
    const ccomplex_t IT_0060 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0061 = s_34*IT_0060;
    const ccomplex_t IT_0062 = (-5.33333333333333)*IT_0041;
    const ccomplex_t IT_0063 = IT_0016*conj(IT_0033)*IT_0043;
    const ccomplex_t IT_0064 = s_34*conj(IT_0011)*IT_0040 + conj(IT_0015)
      *IT_0042 + conj(IT_0029)*IT_0045 + IT_0029*IT_0052 + IT_0033*IT_0053 +
       conj(IT_0006)*IT_0055 + conj(IT_0019)*IT_0056 + IT_0006*IT_0058 + IT_0019
      *IT_0059 + IT_0011*IT_0061 + IT_0015*IT_0062 + IT_0063;
    const ccomplex_t IT_0065 = (-0.666666666666667)*IT_0041;
    const ccomplex_t IT_0066 = -g_s;
    const ccomplex_t IT_0067 = IT_0007*IT_0066;
    const ccomplex_t IT_0068 = IT_0001*IT_0066;
    const ccomplex_t IT_0069 = IT_0048*(IT_0046 + -IT_0047 + IT_0067 + 
      -IT_0068);
    const ccomplex_t IT_0070 = 0.5*IT_0069;
    const ccomplex_t IT_0071 = s_24*IT_0051;
    const ccomplex_t IT_0072 = 0.666666666666667*IT_0041;
    const ccomplex_t IT_0073 = pow(m_sc_L, 4);
    const ccomplex_t IT_0074 = 0.666666666666667*IT_0073;
    const ccomplex_t IT_0075 = IT_0046*IT_0048;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = 1.33333333333333*s_12;
    const ccomplex_t IT_0078 = s_23*IT_0077;
    const ccomplex_t IT_0079 = IT_0048*IT_0067;
    const ccomplex_t IT_0080 = 0.5*IT_0079;
    const ccomplex_t IT_0081 = s_24*IT_0077;
    const ccomplex_t IT_0082 = IT_0048*(IT_0046 + -IT_0047 + IT_0067);
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = IT_0016*IT_0077;
    const ccomplex_t IT_0085 = 0.5*IT_0049;
    const ccomplex_t IT_0086 = 1.33333333333333*s_23;
    const ccomplex_t IT_0087 = IT_0016*IT_0086;
    const ccomplex_t IT_0088 = IT_0048*(IT_0067 + -IT_0068);
    const ccomplex_t IT_0089 = 0.5*IT_0088;
    const ccomplex_t IT_0090 = 1.33333333333333*s_24;
    const ccomplex_t IT_0091 = IT_0016*IT_0090;
    const ccomplex_t IT_0092 = pow(g_s, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = s_13 + s_23;
    const ccomplex_t IT_0095 = IT_0067*IT_0094;
    const ccomplex_t IT_0096 = s_14 + s_24;
    const ccomplex_t IT_0097 = IT_0007*IT_0096;
    const ccomplex_t IT_0098 = g_s*IT_0097;
    const ccomplex_t IT_0099 = IT_0095 + IT_0098;
    const ccomplex_t IT_0100 = g_s*s_24;
    const ccomplex_t IT_0101 = IT_0001*IT_0100;
    const ccomplex_t IT_0102 = s_23*IT_0068;
    const ccomplex_t IT_0103 = -IT_0101 + -IT_0102;
    const ccomplex_t IT_0104 = IT_0099 + IT_0103;
    const ccomplex_t IT_0105 = IT_0048*IT_0104;
    const ccomplex_t IT_0106 = IT_0093 + 0.5*IT_0105;
    const ccomplex_t IT_0107 = 1.33333333333333*IT_0016;
    const ccomplex_t IT_0108 = s_24*IT_0054;
    const ccomplex_t IT_0109 = IT_0016*IT_0054;
    const ccomplex_t IT_0110 = conj(IT_0006)*IT_0109;
    const ccomplex_t IT_0111 = (-0.5)*IT_0075;
    const ccomplex_t IT_0112 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0113 = s_23*IT_0112;
    const ccomplex_t IT_0114 = (-0.5)*IT_0079;
    const ccomplex_t IT_0115 = s_24*IT_0112;
    const ccomplex_t IT_0116 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0117 = IT_0016*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0088;
    const ccomplex_t IT_0119 = (-10.6666666666667)*s_24;
    const ccomplex_t IT_0120 = IT_0016*IT_0119;
    const ccomplex_t IT_0121 = IT_0093 + (-0.5)*IT_0105;
    const ccomplex_t IT_0122 = (-10.6666666666667)*IT_0016;
    const ccomplex_t IT_0123 = IT_0121*IT_0122;
    const ccomplex_t IT_0124 = s_24*IT_0057;
    const ccomplex_t IT_0125 = IT_0016*IT_0057;
    const ccomplex_t IT_0126 = IT_0006*IT_0125;
    const ccomplex_t IT_0127 = s_24*IT_0060;
    const ccomplex_t IT_0128 = IT_0016*IT_0060;
    const ccomplex_t IT_0129 = (-5.33333333333333)*IT_0073;
    const ccomplex_t IT_0130 = (-0.666666666666667)*IT_0073;
    const ccomplex_t IT_0131 = (s_24*conj(IT_0011) + conj(IT_0015)*IT_0016)
      *IT_0040 + IT_0016*conj(IT_0037)*IT_0043 + conj(IT_0029)*IT_0044 + IT_0037
      *IT_0053 + conj(IT_0025)*IT_0065 + IT_0029*IT_0071 + IT_0025*IT_0072 +
       IT_0033*IT_0074 + IT_0076*IT_0078 + IT_0080*IT_0081 + IT_0083*IT_0084 +
       IT_0085*IT_0087 + IT_0089*IT_0091 + IT_0106*IT_0107 + conj(IT_0019)
      *IT_0108 + IT_0110 + IT_0111*IT_0113 + IT_0114*IT_0115 + IT_0050*IT_0117 +
       IT_0118*IT_0120 + IT_0123 + IT_0019*IT_0124 + IT_0126 + IT_0011*IT_0127 +
       IT_0015*IT_0128 + IT_0070*IT_0129 + conj(IT_0033)*IT_0130;
    const ccomplex_t IT_0132 = 0.666666666666667*s_12;
    const ccomplex_t IT_0133 = s_14*IT_0132;
    const ccomplex_t IT_0134 = 0.5*IT_0024;
    const ccomplex_t IT_0135 = IT_0056*conj(IT_0134);
    const ccomplex_t IT_0136 = s_34*IT_0012;
    const ccomplex_t IT_0137 = 5.33333333333333*s_14;
    const ccomplex_t IT_0138 = s_34*IT_0137;
    const ccomplex_t IT_0139 = IT_0059*IT_0134;
    const ccomplex_t IT_0140 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0141 = s_34*IT_0140;
    const ccomplex_t IT_0142 = (-0.666666666666667)*s_14;
    const ccomplex_t IT_0143 = s_34*IT_0142;
    const ccomplex_t IT_0144 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0145 = 1.5*IT_0106;
    const ccomplex_t IT_0146 = 1.33333333333333*s_14;
    const ccomplex_t IT_0147 = 1.5*IT_0121;
    const ccomplex_t IT_0148 = 1.5*IT_0033;
    const ccomplex_t IT_0149 = 1.5*conj(IT_0033);
    const ccomplex_t IT_0150 = s_24*IT_0140;
    const ccomplex_t IT_0151 = 1.5*IT_0037;
    const ccomplex_t IT_0152 = 1.5*conj(IT_0037);
    const ccomplex_t IT_0153 = (-5.33333333333333)*s_14;
    const ccomplex_t IT_0154 = s_34*IT_0153;
    const ccomplex_t IT_0155 = 1.5*IT_0029;
    const ccomplex_t IT_0156 = s_14*IT_0034;
    const ccomplex_t IT_0157 = 1.5*conj(IT_0006);
    const ccomplex_t IT_0158 = 0.666666666666667*s_14;
    const ccomplex_t IT_0159 = s_34*IT_0158;
    const ccomplex_t IT_0160 = 1.5*IT_0011;
    const ccomplex_t IT_0161 = -conj(IT_0019);
    const ccomplex_t IT_0162 = IT_0019 + IT_0161;
    const ccomplex_t IT_0163 = s_14*IT_0162;
    const ccomplex_t IT_0164 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0165 = 1.5*IT_0076;
    const ccomplex_t IT_0166 = 1.33333333333333*s_13;
    const ccomplex_t IT_0167 = 1.5*IT_0111;
    const ccomplex_t IT_0168 = s_24*IT_0164;
    const ccomplex_t IT_0169 = 1.5*IT_0085;
    const ccomplex_t IT_0170 = s_24*IT_0166;
    const ccomplex_t IT_0171 = 1.5*IT_0050;
    const ccomplex_t IT_0172 = s_34*IT_0132;
    const ccomplex_t IT_0173 = 1.5*IT_0015;
    const ccomplex_t IT_0174 = 1.5*conj(IT_0015);
    const ccomplex_t IT_0175 = IT_0006*IT_0133 + IT_0135 + conj(IT_0025)
      *IT_0136 + conj(IT_0029)*IT_0138 + IT_0139 + IT_0025*IT_0141 + conj
      (IT_0011)*IT_0143 + 0.666666666666667*IT_0144*IT_0145 + 0.666666666666667
      *IT_0146*IT_0147 + 0.666666666666667*IT_0124*IT_0148 + 0.666666666666667
      *IT_0108*IT_0149 + 0.666666666666667*IT_0150*IT_0151 + 0.666666666666667
      *IT_0013*IT_0152 + 0.666666666666667*IT_0154*IT_0155 + 0.666666666666667
      *IT_0156*IT_0157 + 0.666666666666667*IT_0159*IT_0160 + 0.666666666666667
      *s_14*(IT_0163 + IT_0164*IT_0165 + IT_0166*IT_0167) + 0.666666666666667
      *IT_0168*IT_0169 + 0.666666666666667*IT_0170*IT_0171 + 0.666666666666667
      *IT_0172*IT_0173 + 0.666666666666667*IT_0038*IT_0174;
    const ccomplex_t IT_0176 = s_24*IT_0132;
    const ccomplex_t IT_0177 = IT_0134*IT_0172;
    const ccomplex_t IT_0178 = 0.666666666666667*s_13;
    const ccomplex_t IT_0179 = s_24*IT_0178;
    const ccomplex_t IT_0180 = s_34*IT_0178;
    const ccomplex_t IT_0181 = s_14*IT_0054;
    const ccomplex_t IT_0182 = s_14*IT_0057;
    const ccomplex_t IT_0183 = IT_0038*conj(IT_0134);
    const ccomplex_t IT_0184 = (-0.1875)*IT_0076;
    const ccomplex_t IT_0185 = (-0.1875)*IT_0111;
    const ccomplex_t IT_0186 = IT_0036 + IT_0030*conj(IT_0037) + IT_0039 +
       conj(IT_0015)*IT_0056 + IT_0015*IT_0059 + conj(IT_0011)*IT_0138 + conj
      (IT_0029)*IT_0143 + IT_0121*IT_0144 + IT_0106*IT_0146 + IT_0011*IT_0154 +
       IT_0029*IT_0159 + IT_0050*IT_0168 + IT_0085*IT_0170 + IT_0033*IT_0176 +
       IT_0177 + IT_0037*IT_0179 + IT_0025*IT_0180 + conj(IT_0006)*IT_0181 +
       IT_0006*IT_0182 + IT_0183 + (-5.33333333333333)*s_14*(IT_0163 + IT_0166
      *IT_0184 + IT_0164*IT_0185);
    const ccomplex_t IT_0187 = IT_0048*(IT_0046 + IT_0067);
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = s_14*IT_0077;
    const ccomplex_t IT_0190 = IT_0016*IT_0153;
    const ccomplex_t IT_0191 = pow(s_12, 2);
    const ccomplex_t IT_0192 = 0.666666666666667*IT_0191;
    const ccomplex_t IT_0193 = (-0.666666666666667)*IT_0191;
    const ccomplex_t IT_0194 = s_13*IT_0132;
    const ccomplex_t IT_0195 = s_13*IT_0034;
    const ccomplex_t IT_0196 = s_14*IT_0178;
    const ccomplex_t IT_0197 = IT_0016*IT_0137;
    const ccomplex_t IT_0198 = 1.5*conj(IT_0019);
    const ccomplex_t IT_0199 = s_14*IT_0140;
    const ccomplex_t IT_0200 = 1.5*IT_0188;
    const ccomplex_t IT_0201 = 1.5*IT_0083;
    const ccomplex_t IT_0202 = IT_0016*IT_0166;
    const ccomplex_t IT_0203 = IT_0016*IT_0164;
    const ccomplex_t IT_0204 = (-10.6666666666667)*IT_0191;
    const ccomplex_t IT_0205 = 1.5*IT_0070;
    const ccomplex_t IT_0206 = 1.33333333333333*IT_0191;
    const ccomplex_t IT_0207 = (-0.5)*IT_0069;
    const ccomplex_t IT_0208 = 1.5*IT_0207;
    const ccomplex_t IT_0209 = s_13*IT_0077;
    const ccomplex_t IT_0210 = s_13*IT_0112;
    const ccomplex_t IT_0211 = s_14*IT_0112;
    const ccomplex_t IT_0212 = 1.5*IT_0118;
    const ccomplex_t IT_0213 = s_13*IT_0057;
    const ccomplex_t IT_0214 = s_13*IT_0054;
    const ccomplex_t IT_0215 = IT_0016*IT_0146;
    const ccomplex_t IT_0216 = 1.5*IT_0080;
    const ccomplex_t IT_0217 = IT_0016*IT_0144;
    const ccomplex_t IT_0218 = 1.5*IT_0114;
    const ccomplex_t IT_0219 = pow(s_13, 2);
    const ccomplex_t IT_0220 = (IT_0025 + -conj(IT_0025))*IT_0219;
    const ccomplex_t IT_0221 = s_14*(conj(IT_0011)*IT_0012 + IT_0026*conj
      (IT_0029)) + IT_0110 + IT_0123 + IT_0126 + 0.666666666666667*IT_0107
      *IT_0145 + IT_0089*IT_0189 + IT_0019*IT_0190 + 0.666666666666667*IT_0148
      *IT_0192 + 0.666666666666667*IT_0149*IT_0193 + 0.666666666666667*IT_0151
      *IT_0194 + 0.666666666666667*IT_0152*IT_0195 + 0.666666666666667*IT_0155
      *IT_0196 + 0.666666666666667*IT_0197*IT_0198 + 0.666666666666667*IT_0160
      *IT_0199 + 0.666666666666667*IT_0129*IT_0200 + 0.666666666666667*IT_0084
      *IT_0201 + 0.666666666666667*IT_0165*IT_0202 + 0.666666666666667*IT_0167
      *IT_0203 + 0.666666666666667*IT_0204*IT_0205 + 0.666666666666667*IT_0206
      *IT_0208 + 0.666666666666667*IT_0169*IT_0209 + 0.666666666666667*IT_0171
      *IT_0210 + 0.666666666666667*IT_0211*IT_0212 + 0.666666666666667*IT_0173
      *IT_0213 + 0.666666666666667*IT_0174*IT_0214 + 0.666666666666667*IT_0215
      *IT_0216 + 0.666666666666667*IT_0217*IT_0218 + 0.666666666666667*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0187;
    const ccomplex_t IT_0223 = IT_0016*IT_0132;
    const ccomplex_t IT_0224 = IT_0006*IT_0223;
    const ccomplex_t IT_0225 = IT_0016*IT_0158;
    const ccomplex_t IT_0226 = IT_0107*IT_0121;
    const ccomplex_t IT_0227 = 1.33333333333333*IT_0073;
    const ccomplex_t IT_0228 = 5.33333333333333*IT_0191;
    const ccomplex_t IT_0229 = IT_0016*IT_0112;
    const ccomplex_t IT_0230 = (-5.33333333333333)*IT_0191;
    const ccomplex_t IT_0231 = IT_0016*IT_0034;
    const ccomplex_t IT_0232 = conj(IT_0006)*IT_0231;
    const ccomplex_t IT_0233 = IT_0016*IT_0142;
    const ccomplex_t IT_0234 = s_14*(conj(IT_0011)*IT_0026 + IT_0012*conj
      (IT_0029)) + IT_0106*IT_0122 + IT_0118*IT_0189 + IT_0015*IT_0194 + conj
      (IT_0015)*IT_0195 + IT_0011*IT_0196 + IT_0029*IT_0199 + (-5.33333333333333
      )*IT_0185*IT_0202 + (-5.33333333333333)*IT_0184*IT_0203 + IT_0070*IT_0206 
      + IT_0204*IT_0207 + IT_0050*IT_0209 + IT_0085*IT_0210 + IT_0089*IT_0211 +
       IT_0037*IT_0213 + conj(IT_0037)*IT_0214 + IT_0114*IT_0215 + IT_0080
      *IT_0217 + (-5.33333333333333)*IT_0220 + IT_0129*IT_0222 + IT_0224 +
       IT_0019*IT_0225 + IT_0226 + IT_0188*IT_0227 + conj(IT_0033)*IT_0228 +
       IT_0083*IT_0229 + IT_0033*IT_0230 + IT_0232 + conj(IT_0019)*IT_0233;
    const ccomplex_t IT_0235 = s_23*IT_0012;
    const ccomplex_t IT_0236 = s_24*IT_0137;
    const ccomplex_t IT_0237 = IT_0048*(IT_0046 + IT_0067 + -IT_0068);
    const ccomplex_t IT_0238 = (-0.5)*IT_0237;
    const ccomplex_t IT_0239 = IT_0033 + IT_0207;
    const ccomplex_t IT_0240 = IT_0037 + IT_0085;
    const ccomplex_t IT_0241 = 5.33333333333333*s_34;
    const ccomplex_t IT_0242 = IT_0016*IT_0241;
    const ccomplex_t IT_0243 = s_24*IT_0142;
    const ccomplex_t IT_0244 = 0.5*IT_0237;
    const ccomplex_t IT_0245 = s_23*IT_0034;
    const ccomplex_t IT_0246 = (-0.666666666666667)*s_24;
    const ccomplex_t IT_0247 = IT_0016*IT_0246;
    const ccomplex_t IT_0248 = (-0.666666666666667)*s_34;
    const ccomplex_t IT_0249 = IT_0016*IT_0248;
    const ccomplex_t IT_0250 = (-8)*IT_0089;
    const ccomplex_t IT_0251 = IT_0118 + IT_0250;
    const ccomplex_t IT_0252 = s_24*IT_0251;
    const ccomplex_t IT_0253 = IT_0019*IT_0030 + IT_0015*IT_0044 + IT_0020
      *IT_0106 + IT_0083*(IT_0017 + IT_0108) + IT_0027*IT_0111 + IT_0188*
      (IT_0156 + IT_0195) + (IT_0181 + IT_0214)*IT_0222 + IT_0076*IT_0235 +
       IT_0080*IT_0236 + (IT_0055 + IT_0197)*IT_0238 + IT_0021*IT_0239 + IT_0040
      *(s_34*IT_0025 + IT_0106 + IT_0016*IT_0239) + IT_0042*IT_0240 + (IT_0089 +
       IT_0240)*IT_0242 + IT_0114*IT_0243 + IT_0244*(IT_0233 + IT_0245) +
       IT_0121*(IT_0043 + IT_0246) + IT_0070*(IT_0016*IT_0043 + IT_0247) +
       IT_0118*IT_0249 + IT_0050*(IT_0065 + IT_0249) + (-0.666666666666667)*s_24
      *IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0082;
    const ccomplex_t IT_0255 = IT_0016*IT_0178;
    const ccomplex_t IT_0256 = s_23*IT_0132;
    const ccomplex_t IT_0257 = IT_0051*IT_0121;
    const ccomplex_t IT_0258 = (-5.33333333333333)*s_24;
    const ccomplex_t IT_0259 = 0.666666666666667*s_24;
    const ccomplex_t IT_0260 = IT_0016*IT_0026;
    const ccomplex_t IT_0261 = IT_0016*IT_0140;
    const ccomplex_t IT_0262 = s_23*IT_0140;
    const ccomplex_t IT_0263 = s_23*IT_0178;
    const ccomplex_t IT_0264 = IT_0016*IT_0259;
    const ccomplex_t IT_0265 = IT_0016*IT_0258;
    const ccomplex_t IT_0266 = 0.666666666666667*s_34;
    const ccomplex_t IT_0267 = IT_0016*IT_0266;
    const ccomplex_t IT_0268 = s_24*IT_0153;
    const ccomplex_t IT_0269 = s_24*IT_0158;
    const ccomplex_t IT_0270 = s_34*conj(IT_0025)*IT_0040 + 0.666666666666667*
      (IT_0021 + IT_0016*IT_0040)*IT_0149 + 0.666666666666667*IT_0062*IT_0169 +
       0.666666666666667*IT_0044*IT_0174 + 0.666666666666667*IT_0030*IT_0198 +
       0.666666666666667*(IT_0133 + IT_0194)*IT_0200 + (IT_0182 + IT_0213)
      *IT_0222 + (IT_0058 + IT_0190)*IT_0238 + 0.666666666666667*IT_0152*
      (IT_0042 + IT_0242) + 0.666666666666667*s_24*(16*s_23*conj(IT_0134) +
       IT_0252) + IT_0254*(IT_0176 + IT_0255) + IT_0244*(IT_0225 + IT_0256) +
       IT_0257 + 0.666666666666667*IT_0145*(IT_0060 + IT_0258) +
       0.666666666666667*IT_0147*IT_0259 + 0.666666666666667*IT_0157*(IT_0035 +
       IT_0260) + 0.666666666666667*IT_0201*(IT_0124 + IT_0261) +
       0.666666666666667*IT_0165*IT_0262 + 0.666666666666667*IT_0167*IT_0263 +
       0.666666666666667*IT_0205*(IT_0053 + IT_0264) + 0.666666666666667*IT_0208
      *(IT_0128 + IT_0265) + 0.666666666666667*IT_0212*IT_0267 +
       0.666666666666667*IT_0171*(IT_0072 + IT_0267) + 0.666666666666667*IT_0216
      *IT_0268 + 0.666666666666667*IT_0218*IT_0269;
    const ccomplex_t IT_0271 = (-5.33333333333333)*s_34;
    const ccomplex_t IT_0272 = IT_0016*IT_0271;
    const ccomplex_t IT_0273 = conj(IT_0006)*IT_0016;
    const ccomplex_t IT_0274 = conj(IT_0029)*IT_0031 + IT_0006*IT_0128 + conj
      (IT_0011)*IT_0136 + IT_0011*IT_0141 + IT_0177 + IT_0029*IT_0180 + IT_0183 
      + conj(IT_0015)*IT_0235 + conj(IT_0019)*IT_0242 + conj(IT_0033)*IT_0245 +
       IT_0033*IT_0256 + IT_0015*IT_0262 + IT_0019*IT_0272 + IT_0040*IT_0273;
    const ccomplex_t IT_0275 = conj(IT_0015)*IT_0027 + conj(IT_0011)*IT_0031 +
       IT_0006*IT_0053 + conj(IT_0033)*IT_0055 + IT_0033*IT_0058 + IT_0135 +
       conj(IT_0029)*IT_0136 + IT_0139 + IT_0029*IT_0141 + IT_0011*IT_0180 +
       conj(IT_0019)*IT_0249 + IT_0015*IT_0263 + IT_0019*IT_0267 + IT_0043
      *IT_0273;
    const ccomplex_t IT_0276 = IT_0070 + IT_0188;
    const ccomplex_t IT_0277 = IT_0077*IT_0121;
    const ccomplex_t IT_0278 = conj(IT_0025)*IT_0235;
    const ccomplex_t IT_0279 = IT_0207 + IT_0222;
    const ccomplex_t IT_0280 = IT_0106*IT_0112;
    const ccomplex_t IT_0281 = IT_0025*IT_0262;
    const ccomplex_t IT_0282 = IT_0033*IT_0125;
    const ccomplex_t IT_0283 = conj(IT_0033)*IT_0109;
    const ccomplex_t IT_0284 = IT_0017*conj(IT_0037) + IT_0091*IT_0114 +
       IT_0089*IT_0115 + IT_0081*IT_0118 + IT_0080*IT_0120 + IT_0019*IT_0133 +
       conj(IT_0019)*IT_0156 + s_23*(conj(IT_0029)*IT_0137 + conj(IT_0011)
      *IT_0142 + IT_0029*IT_0153 + IT_0011*IT_0158) + IT_0006*IT_0192 + conj
      (IT_0006)*IT_0193 + IT_0050*IT_0202 + IT_0085*IT_0203 + IT_0083*IT_0204 +
       IT_0111*IT_0209 + IT_0076*IT_0210 + IT_0129*IT_0238 + IT_0227*IT_0244 +
       conj(IT_0015)*IT_0245 + IT_0015*IT_0256 + IT_0037*IT_0261 + IT_0084
      *IT_0276 + IT_0277 + IT_0278 + IT_0229*IT_0279 + IT_0280 + IT_0281 +
       IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = IT_0013*conj(IT_0029) + conj(IT_0011)*IT_0030 +
       conj(IT_0037)*IT_0055 + IT_0037*IT_0058 + IT_0006*IT_0074 + IT_0050
      *IT_0078 + IT_0087*IT_0111 + IT_0085*IT_0113 + IT_0076*IT_0117 + IT_0083
      *IT_0129 + conj(IT_0006)*IT_0130 + IT_0029*IT_0150 + IT_0011*IT_0179 +
       IT_0114*IT_0189 + IT_0080*IT_0211 + IT_0118*IT_0215 + IT_0089*IT_0217 +
       conj(IT_0019)*IT_0247 + IT_0015*IT_0255 + conj(IT_0015)*IT_0260 + IT_0019
      *IT_0264 + IT_0277 + IT_0278 + IT_0280 + IT_0281 + IT_0282 + IT_0283;
    const ccomplex_t IT_0286 = IT_0025*IT_0263;
    const ccomplex_t IT_0287 = conj(IT_0033)*IT_0231;
    const ccomplex_t IT_0288 = conj(IT_0025)*IT_0027;
    const ccomplex_t IT_0289 = IT_0077*IT_0106;
    const ccomplex_t IT_0290 = IT_0112*IT_0121;
    const ccomplex_t IT_0291 = IT_0033*IT_0223;
    const ccomplex_t IT_0292 = IT_0229*IT_0276;
    const ccomplex_t IT_0293 = IT_0084*IT_0279;
    const ccomplex_t IT_0294 = conj(IT_0015)*IT_0055 + IT_0015*IT_0058 +
       IT_0081*IT_0089 + IT_0080*IT_0091 + IT_0115*IT_0118 + IT_0114*IT_0120 +
       s_23*(conj(IT_0011)*IT_0137 + conj(IT_0029)*IT_0142 + IT_0011*IT_0153 +
       IT_0029*IT_0158) + conj(IT_0019)*IT_0181 + IT_0019*IT_0182 + IT_0085
      *IT_0202 + IT_0050*IT_0203 + IT_0083*IT_0206 + IT_0076*IT_0209 + IT_0111
      *IT_0210 + conj(IT_0006)*IT_0228 + IT_0006*IT_0230 + IT_0129*IT_0244 +
       IT_0037*IT_0255 + conj(IT_0037)*IT_0260 + IT_0286 + IT_0287 + IT_0288 +
       IT_0289 + IT_0290 + IT_0291 + IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = IT_0134*IT_0272;
    const ccomplex_t IT_0296 = conj(IT_0037)*IT_0040;
    const ccomplex_t IT_0297 = IT_0085*IT_0116;
    const ccomplex_t IT_0298 = IT_0050*IT_0086;
    const ccomplex_t IT_0299 = IT_0021*conj(IT_0033) + IT_0025*IT_0061 +
       IT_0106*IT_0119 + IT_0090*IT_0121 + IT_0037*IT_0127 + s_23*(IT_0076
      *IT_0144 + IT_0111*IT_0146) + IT_0006*IT_0225 + conj(IT_0006)*IT_0233 +
       conj(IT_0019)*IT_0243 + conj(IT_0015)*IT_0249 + s_34*(IT_0020*conj
      (IT_0029) + conj(IT_0025)*IT_0040 + conj(IT_0011)*IT_0246 + IT_0029
      *IT_0258 + IT_0011*IT_0259) + IT_0033*IT_0265 + IT_0015*IT_0267 + IT_0019
      *IT_0269 + IT_0295 + s_24*(IT_0296 + IT_0297 + IT_0298);
    const ccomplex_t IT_0300 = IT_0033 + -conj(IT_0033);
    const ccomplex_t IT_0301 = conj(IT_0015)*IT_0043;
    const ccomplex_t IT_0302 = conj(IT_0006)*IT_0034 + IT_0015*IT_0051 +
       IT_0037*IT_0060 + (-21.3333333333333)*IT_0106 + 5.33333333333333*IT_0121 
      + IT_0006*IT_0132 + conj(IT_0019)*IT_0142 + IT_0019*IT_0158 + IT_0076
      *IT_0164 + IT_0111*IT_0166 + conj(IT_0029)*IT_0241 + conj(IT_0011)*IT_0248
       + IT_0011*IT_0266 + IT_0029*IT_0271 + IT_0296 + IT_0297 + IT_0298 + (
      -5.33333333333333)*IT_0016*IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = IT_0085*IT_0086;
    const ccomplex_t IT_0304 = IT_0050*IT_0116;
    const ccomplex_t IT_0305 = conj(IT_0015)*IT_0040 + conj(IT_0037)*IT_0043 +
       conj(IT_0006)*IT_0054 + IT_0006*IT_0057 + IT_0015*IT_0060 + (
      -21.3333333333333)*IT_0121 + conj(IT_0019)*IT_0137 + IT_0019*IT_0153 +
       IT_0111*IT_0164 + IT_0076*IT_0166 + conj(IT_0011)*IT_0241 + conj(IT_0029)
      *IT_0248 + IT_0029*IT_0266 + IT_0011*IT_0271 + 0.666666666666667*IT_0016
      *IT_0300 + IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = conj(IT_0037)*IT_0044 + conj(IT_0025)*IT_0045 +
       IT_0025*IT_0052 + IT_0037*IT_0071 + IT_0090*IT_0106 + IT_0119*IT_0121 +
       s_23*(IT_0111*IT_0144 + IT_0076*IT_0146) + IT_0006*IT_0190 + conj(IT_0006
      )*IT_0197 + conj(IT_0019)*IT_0236 + conj(IT_0015)*IT_0242 + conj(IT_0033)
      *IT_0247 + s_34*(conj(IT_0011)*IT_0020 + conj(IT_0029)*IT_0246 + IT_0011
      *IT_0258 + IT_0029*IT_0259) + IT_0033*IT_0264 + IT_0019*IT_0268 + IT_0015
      *IT_0272 + s_24*(IT_0303 + IT_0304);
    const ccomplex_t IT_0307 = s_24*conj(IT_0029)*IT_0040;
    const ccomplex_t IT_0308 = IT_0016*(conj(IT_0037)*IT_0040 + conj(IT_0015)
      *IT_0043);
    const ccomplex_t IT_0309 = 5.33333333333333*IT_0073;
    const ccomplex_t IT_0310 = conj(IT_0033)*IT_0309;
    const ccomplex_t IT_0311 = conj(IT_0019)*IT_0035;
    const ccomplex_t IT_0312 = conj(IT_0011)*IT_0044;
    const ccomplex_t IT_0313 = conj(IT_0025)*IT_0042;
    const ccomplex_t IT_0314 = IT_0015*IT_0053 + IT_0025*IT_0062 + IT_0011
      *IT_0071 + IT_0050*IT_0087 + IT_0078*IT_0111 + IT_0076*IT_0113 + IT_0081
      *IT_0114 + IT_0080*IT_0115 + IT_0085*IT_0117 + IT_0091*IT_0118 + IT_0089
      *IT_0120 + IT_0106*IT_0122 + IT_0029*IT_0127 + IT_0037*IT_0128 + IT_0019
      *IT_0176 + IT_0224 + IT_0226 + IT_0070*IT_0227 + IT_0083*IT_0229 + IT_0232
       + IT_0129*IT_0239 + IT_0307 + IT_0308 + IT_0310 + IT_0311 + IT_0312 +
       IT_0313;
    const ccomplex_t IT_0315 = IT_0035*conj(IT_0134);
    const ccomplex_t IT_0316 = conj(IT_0029)*IT_0030;
    const ccomplex_t IT_0317 = conj(IT_0134)*IT_0260;
    const ccomplex_t IT_0318 = conj(IT_0019)*IT_0021;
    const ccomplex_t IT_0319 = conj(IT_0011)*IT_0013;
    const ccomplex_t IT_0320 = conj(IT_0015)*IT_0017;
    const ccomplex_t IT_0321 = IT_0078*IT_0085 + IT_0076*IT_0087 + IT_0050
      *IT_0113 + IT_0111*IT_0117 + IT_0011*IT_0150 + IT_0029*IT_0179 + IT_0080
      *IT_0189 + IT_0114*IT_0211 + IT_0089*IT_0215 + IT_0118*IT_0217 + IT_0083
      *IT_0227 + IT_0206*IT_0238 + IT_0204*IT_0244 + conj(IT_0037)*IT_0245 +
       IT_0129*(IT_0006 + IT_0254) + IT_0037*IT_0256 + IT_0015*IT_0261 + IT_0019
      *IT_0265 + IT_0286 + IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291 +
       IT_0292 + IT_0293 + conj(IT_0006)*IT_0309 + IT_0315 + IT_0316 + IT_0317 +
       IT_0318 + IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = s_34*conj(IT_0029);
    const ccomplex_t IT_0323 = IT_0016*conj(IT_0033) + IT_0322;
    const ccomplex_t IT_0324 = IT_0040*IT_0323;
    const ccomplex_t IT_0325 = conj(IT_0006)*IT_0245;
    const ccomplex_t IT_0326 = conj(IT_0019)*IT_0038;
    const ccomplex_t IT_0327 = conj(IT_0011)*IT_0045;
    const ccomplex_t IT_0328 = conj(IT_0015)*IT_0065;
    const ccomplex_t IT_0329 = conj(IT_0006)*(IT_0011*IT_0013 + IT_0015
      *IT_0017 + IT_0019*IT_0021 + IT_0025*IT_0027 + IT_0029*IT_0030) + conj
      (IT_0019)*IT_0025*IT_0031 + IT_0019*(IT_0013*conj(IT_0015) + IT_0036 +
       conj(IT_0037)*IT_0038 + IT_0039) + conj(IT_0033)*(s_24*IT_0029*IT_0040 +
       IT_0025*IT_0042 + IT_0011*IT_0044) + conj(IT_0037)*(s_34*IT_0029*IT_0040 
      + IT_0011*IT_0045) + IT_0050*IT_0064 + IT_0015*(IT_0013*conj(IT_0019) +
       IT_0063 + conj(IT_0037)*IT_0065) + IT_0070*IT_0131 + IT_0080*IT_0175 +
       IT_0114*IT_0186 + IT_0188*IT_0221 + IT_0222*IT_0234 + conj(IT_0134)
      *IT_0253 + IT_0134*IT_0270 + IT_0111*IT_0274 + IT_0076*IT_0275 + IT_0238
      *IT_0284 + IT_0083*IT_0285 + IT_0244*IT_0294 + IT_0089*IT_0299 + IT_0106
      *IT_0302 + IT_0121*IT_0305 + IT_0118*IT_0306 + IT_0033*(conj(IT_0006)
      *IT_0231 + IT_0307 + IT_0308 + IT_0310 + IT_0311 + IT_0312 + IT_0313) +
       IT_0207*IT_0314 + IT_0006*(conj(IT_0037)*IT_0245 + IT_0287 + IT_0288 +
       conj(IT_0006)*IT_0309 + IT_0315 + IT_0316 + IT_0317 + IT_0318 + IT_0319 +
       IT_0320) + IT_0254*IT_0321 + IT_0037*(IT_0257 + IT_0324 + IT_0325 +
       IT_0326 + IT_0327 + IT_0328) + IT_0085*(IT_0011*IT_0052 + IT_0029*IT_0061
       + IT_0015*IT_0072 + IT_0033*IT_0128 + IT_0019*IT_0172 + IT_0006*IT_0256 +
       IT_0295 + IT_0324 + IT_0325 + IT_0326 + IT_0327 + IT_0328);
    return create_ccomplex_return(IT_0329);
}

