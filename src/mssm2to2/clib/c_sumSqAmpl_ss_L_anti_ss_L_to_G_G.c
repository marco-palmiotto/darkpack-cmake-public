#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_L_anti_ss_L_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_L_anti_ss_L_to_G_G(
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
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t Gamma_sl = param->Gamma_sl;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = (-2)*IT_0000;
    const ccomplex_t IT_0002 = cpow(IT_0001, 2);
    const ccomplex_t IT_0003 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_ss_L
      *Gamma_sl + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = -IT_0000;
    const ccomplex_t IT_0008 = cpow(IT_0007, 2);
    const ccomplex_t IT_0009 = IT_0004*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = IT_0001*IT_0007;
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0012;
    const ccomplex_t IT_0015 = 0.5*IT_0005;
    const ccomplex_t IT_0016 = pow(m_ss_L, 2);
    const ccomplex_t IT_0017 = conj(IT_0013) + conj(IT_0015);
    const ccomplex_t IT_0018 = 5.33333333333333*s_23;
    const ccomplex_t IT_0019 = (s_24*(conj(IT_0006)*IT_0010 + IT_0006*conj
      (IT_0010)) + s_34*(conj(IT_0013)*IT_0014 + IT_0013*conj(IT_0014) + conj
      (IT_0010)*IT_0015 + IT_0010*conj(IT_0015)) + IT_0016*(conj(IT_0006)*
      (IT_0013 + IT_0015) + IT_0006*IT_0017))*IT_0018;
    const ccomplex_t IT_0020 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_ss_L
      *Gamma_sl + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = IT_0011*IT_0021;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = (-0.666666666666667)*s_13;
    const ccomplex_t IT_0025 = s_34*IT_0024;
    const ccomplex_t IT_0026 = (IT_0002 + -IT_0011)*IT_0021;
    const ccomplex_t IT_0027 = (-0.5)*IT_0026;
    const ccomplex_t IT_0028 = 5.33333333333333*s_13;
    const ccomplex_t IT_0029 = s_24*IT_0028;
    const ccomplex_t IT_0030 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0031 = s_24*IT_0030;
    const ccomplex_t IT_0032 = conj(IT_0006)*IT_0031;
    const ccomplex_t IT_0033 = s_34*IT_0030;
    const ccomplex_t IT_0034 = s_24*IT_0024;
    const ccomplex_t IT_0035 = conj(IT_0014)*IT_0025;
    const ccomplex_t IT_0036 = IT_0002*IT_0021;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = -(IT_0008 + -IT_0011)*IT_0021;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0016*IT_0028;
    const ccomplex_t IT_0041 = 5.33333333333333*s_24;
    const ccomplex_t IT_0042 = IT_0016*IT_0041;
    const ccomplex_t IT_0043 = pow(m_ss_L, 4);
    const ccomplex_t IT_0044 = 5.33333333333333*IT_0043;
    const ccomplex_t IT_0045 = s_23*IT_0030;
    const ccomplex_t IT_0046 = IT_0016*IT_0030;
    const ccomplex_t IT_0047 = s_23*IT_0024;
    const ccomplex_t IT_0048 = IT_0016*IT_0024;
    const ccomplex_t IT_0049 = -g_s;
    const ccomplex_t IT_0050 = IT_0007*IT_0049;
    const ccomplex_t IT_0051 = cpow(s_12 + IT_0016 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = 0.666666666666667*s_12;
    const ccomplex_t IT_0055 = s_24*IT_0054;
    const ccomplex_t IT_0056 = s_34*IT_0054;
    const ccomplex_t IT_0057 = IT_0013*IT_0056;
    const ccomplex_t IT_0058 = IT_0016*IT_0054;
    const ccomplex_t IT_0059 = IT_0037*IT_0058;
    const ccomplex_t IT_0060 = 0.666666666666667*s_13;
    const ccomplex_t IT_0061 = s_24*IT_0060;
    const ccomplex_t IT_0062 = s_34*IT_0060;
    const ccomplex_t IT_0063 = 0.666666666666667*s_14;
    const ccomplex_t IT_0064 = IT_0016*IT_0063;
    const ccomplex_t IT_0065 = IT_0023*IT_0064;
    const ccomplex_t IT_0066 = g_s*IT_0001;
    const ccomplex_t IT_0067 = IT_0051*(IT_0050 + -IT_0066);
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = g_s*IT_0007;
    const ccomplex_t IT_0070 = IT_0001*IT_0049;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0069 + IT_0071;
    const ccomplex_t IT_0073 = IT_0051*IT_0072;
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = IT_0068 + IT_0074;
    const ccomplex_t IT_0076 = 1.33333333333333*s_12;
    const ccomplex_t IT_0077 = s_24*IT_0076;
    const ccomplex_t IT_0078 = -conj(IT_0023);
    const ccomplex_t IT_0079 = IT_0023 + IT_0078;
    const ccomplex_t IT_0080 = IT_0051*IT_0069;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = 1.33333333333333*s_13;
    const ccomplex_t IT_0083 = pow(g_s, 2);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = s_13 + s_23;
    const ccomplex_t IT_0086 = IT_0050*IT_0085;
    const ccomplex_t IT_0087 = s_14 + s_24;
    const ccomplex_t IT_0088 = IT_0007*IT_0087;
    const ccomplex_t IT_0089 = g_s*IT_0088;
    const ccomplex_t IT_0090 = IT_0086 + IT_0089;
    const ccomplex_t IT_0091 = g_s*s_24;
    const ccomplex_t IT_0092 = IT_0001*IT_0091;
    const ccomplex_t IT_0093 = s_23*IT_0070;
    const ccomplex_t IT_0094 = -IT_0092 + -IT_0093;
    const ccomplex_t IT_0095 = IT_0090 + IT_0094;
    const ccomplex_t IT_0096 = IT_0051*IT_0095;
    const ccomplex_t IT_0097 = IT_0084 + 0.5*IT_0096;
    const ccomplex_t IT_0098 = 1.33333333333333*s_14;
    const ccomplex_t IT_0099 = IT_0051*(IT_0050 + IT_0069 + -IT_0070);
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = 1.33333333333333*s_24;
    const ccomplex_t IT_0102 = IT_0016*IT_0101;
    const ccomplex_t IT_0103 = IT_0084 + (-0.5)*IT_0096;
    const ccomplex_t IT_0104 = 1.33333333333333*IT_0016;
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = 5.33333333333333*s_12;
    const ccomplex_t IT_0107 = s_13*IT_0106;
    const ccomplex_t IT_0108 = conj(IT_0015)*IT_0107;
    const ccomplex_t IT_0109 = s_14*IT_0106;
    const ccomplex_t IT_0110 = conj(IT_0013)*IT_0109;
    const ccomplex_t IT_0111 = pow(s_12, 2);
    const ccomplex_t IT_0112 = 5.33333333333333*IT_0111;
    const ccomplex_t IT_0113 = conj(IT_0006)*IT_0112;
    const ccomplex_t IT_0114 = pow(s_13, 2);
    const ccomplex_t IT_0115 = 5.33333333333333*conj(IT_0014)*IT_0114;
    const ccomplex_t IT_0116 = -IT_0051*(IT_0066 + -IT_0069);
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0119 = s_13*IT_0118;
    const ccomplex_t IT_0120 = -IT_0066 + -IT_0070;
    const ccomplex_t IT_0121 = IT_0069 + IT_0120;
    const ccomplex_t IT_0122 = IT_0051*IT_0121;
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = s_24*IT_0118;
    const ccomplex_t IT_0125 = (-0.5)*IT_0116;
    const ccomplex_t IT_0126 = s_13*IT_0076;
    const ccomplex_t IT_0127 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0128 = (-10.6666666666667)*s_14;
    const ccomplex_t IT_0129 = IT_0051*(IT_0050 + IT_0069);
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = 1.33333333333333*IT_0043;
    const ccomplex_t IT_0132 = IT_0016*IT_0128;
    const ccomplex_t IT_0133 = IT_0016*IT_0127;
    const ccomplex_t IT_0134 = (-10.6666666666667)*s_24;
    const ccomplex_t IT_0135 = IT_0016*IT_0134;
    const ccomplex_t IT_0136 = (-10.6666666666667)*IT_0016;
    const ccomplex_t IT_0137 = IT_0097*IT_0136;
    const ccomplex_t IT_0138 = s_24*IT_0063;
    const ccomplex_t IT_0139 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0140 = s_13*IT_0139;
    const ccomplex_t IT_0141 = IT_0015*IT_0140;
    const ccomplex_t IT_0142 = s_14*IT_0139;
    const ccomplex_t IT_0143 = IT_0013*IT_0142;
    const ccomplex_t IT_0144 = s_13*IT_0054;
    const ccomplex_t IT_0145 = s_34*IT_0139;
    const ccomplex_t IT_0146 = (-5.33333333333333)*IT_0111;
    const ccomplex_t IT_0147 = IT_0006*IT_0146;
    const ccomplex_t IT_0148 = s_34*IT_0063;
    const ccomplex_t IT_0149 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0150 = s_14*IT_0149;
    const ccomplex_t IT_0151 = IT_0014*IT_0114;
    const ccomplex_t IT_0152 = (-5.33333333333333)*IT_0151;
    const ccomplex_t IT_0153 = s_14*IT_0060;
    const ccomplex_t IT_0154 = (-5.33333333333333)*s_14;
    const ccomplex_t IT_0155 = s_34*IT_0154;
    const ccomplex_t IT_0156 = IT_0016*IT_0098;
    const ccomplex_t IT_0157 = (-5.33333333333333)*IT_0043;
    const ccomplex_t IT_0158 = s_34*IT_0106;
    const ccomplex_t IT_0159 = s_13*IT_0030;
    const ccomplex_t IT_0160 = conj(IT_0037)*IT_0046;
    const ccomplex_t IT_0161 = 5.33333333333333*s_14;
    const ccomplex_t IT_0162 = s_34*IT_0161;
    const ccomplex_t IT_0163 = s_14*IT_0024;
    const ccomplex_t IT_0164 = IT_0033 + IT_0107;
    const ccomplex_t IT_0165 = (-0.666666666666667)*s_14;
    const ccomplex_t IT_0166 = s_24*IT_0165;
    const ccomplex_t IT_0167 = s_14*IT_0028;
    const ccomplex_t IT_0168 = s_34*IT_0165;
    const ccomplex_t IT_0169 = IT_0016*IT_0165;
    const ccomplex_t IT_0170 = conj(IT_0023)*IT_0169;
    const ccomplex_t IT_0171 = IT_0032 + conj(IT_0015)*IT_0034 + IT_0035 +
       IT_0006*IT_0055 + IT_0057 + IT_0059 + IT_0015*IT_0061 + IT_0014*IT_0062 +
       IT_0065 + IT_0075*IT_0077 + (-5.33333333333333)*s_14*(s_14*IT_0079 + (
      -0.1875)*IT_0081*IT_0082) + IT_0097*IT_0098 + IT_0100*IT_0102 + IT_0105 +
       IT_0108 + conj(IT_0037)*IT_0109 + IT_0110 + IT_0113 + IT_0115 + IT_0117*
      (s_24*IT_0082 + IT_0119) + IT_0123*IT_0124 + IT_0125*(IT_0126 + s_24
      *IT_0127) + IT_0103*IT_0128 + IT_0130*(IT_0131 + IT_0132) + IT_0081*
      (IT_0133 + IT_0135) + IT_0137 + IT_0013*(IT_0138 + IT_0140) + IT_0141 +
       IT_0037*IT_0142 + IT_0143 + IT_0027*(IT_0144 + IT_0145) + IT_0147 +
       IT_0010*(IT_0148 + IT_0150) + IT_0152 + IT_0039*(IT_0153 + IT_0155) +
       IT_0053*(IT_0156 + IT_0157) + conj(IT_0027)*(IT_0158 + IT_0159) + IT_0160
       + conj(IT_0039)*(IT_0162 + IT_0163) + conj(IT_0013)*(IT_0164 + IT_0166) +
       conj(IT_0010)*(IT_0167 + IT_0168) + IT_0170;
    const ccomplex_t IT_0172 = pow(s_23, 2);
    const ccomplex_t IT_0173 = 5.33333333333333*IT_0172;
    const ccomplex_t IT_0174 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0175 = s_24*IT_0174;
    const ccomplex_t IT_0176 = s_34*IT_0174;
    const ccomplex_t IT_0177 = 0.666666666666667*s_23;
    const ccomplex_t IT_0178 = 0.666666666666667*s_24;
    const ccomplex_t IT_0179 = 0.666666666666667*s_34;
    const ccomplex_t IT_0180 = 1.33333333333333*s_23;
    const ccomplex_t IT_0181 = 5.33333333333333*s_34;
    const ccomplex_t IT_0182 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0183 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0184 = (-5.33333333333333)*s_34;
    const ccomplex_t IT_0185 = (-0.666666666666667)*s_24;
    const ccomplex_t IT_0186 = (-0.666666666666667)*s_34;
    const ccomplex_t IT_0187 = 0.666666666666667*(IT_0006 + -conj(IT_0006))
      *IT_0016 + IT_0018*conj(IT_0027) + (-21.3333333333333)*IT_0103 + conj
      (IT_0037)*IT_0106 + IT_0037*IT_0139 + IT_0023*IT_0154 + conj(IT_0023)
      *IT_0161 + conj(IT_0015)*IT_0174 + IT_0015*IT_0177 + IT_0013*(IT_0177 +
       IT_0178) + IT_0010*IT_0179 + IT_0117*IT_0180 + conj(IT_0039)*IT_0181 +
       IT_0125*IT_0182 + IT_0027*IT_0183 + IT_0039*IT_0184 + conj(IT_0013)*
      (IT_0174 + IT_0185) + conj(IT_0010)*IT_0186;
    const ccomplex_t IT_0188 = (-5.33333333333333)*s_24;
    const ccomplex_t IT_0189 = (-5.33333333333333)*(IT_0006 + -conj(IT_0006))
      *IT_0016 + conj(IT_0015)*IT_0018 + IT_0030*conj(IT_0037) + conj(IT_0013)*
      (IT_0018 + IT_0041) + IT_0037*IT_0054 + IT_0023*IT_0063 + (
      -21.3333333333333)*IT_0097 + 5.33333333333333*IT_0103 + conj(IT_0023)
      *IT_0165 + conj(IT_0027)*IT_0174 + IT_0027*IT_0177 + IT_0039*IT_0179 +
       IT_0125*IT_0180 + conj(IT_0010)*IT_0181 + IT_0117*IT_0182 + IT_0015
      *IT_0183 + IT_0010*IT_0184 + conj(IT_0039)*IT_0186 + IT_0013*(IT_0183 +
       IT_0188);
    const ccomplex_t IT_0190 = s_24*IT_0177;
    const ccomplex_t IT_0191 = IT_0016*IT_0177;
    const ccomplex_t IT_0192 = IT_0013*IT_0191;
    const ccomplex_t IT_0193 = 0.666666666666667*IT_0172;
    const ccomplex_t IT_0194 = IT_0014*IT_0193;
    const ccomplex_t IT_0195 = IT_0016*IT_0178;
    const ccomplex_t IT_0196 = 0.666666666666667*IT_0043;
    const ccomplex_t IT_0197 = s_23*IT_0076;
    const ccomplex_t IT_0198 = IT_0081*IT_0197;
    const ccomplex_t IT_0199 = IT_0051*(IT_0050 + -IT_0066 + IT_0069);
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = IT_0016*IT_0076;
    const ccomplex_t IT_0202 = IT_0100*IT_0201;
    const ccomplex_t IT_0203 = 1.33333333333333*IT_0111;
    const ccomplex_t IT_0204 = IT_0053*IT_0203;
    const ccomplex_t IT_0205 = IT_0016*IT_0180;
    const ccomplex_t IT_0206 = IT_0097*IT_0104;
    const ccomplex_t IT_0207 = s_24*IT_0106;
    const ccomplex_t IT_0208 = conj(IT_0023)*IT_0207;
    const ccomplex_t IT_0209 = IT_0016*IT_0106;
    const ccomplex_t IT_0210 = conj(IT_0037)*IT_0209;
    const ccomplex_t IT_0211 = 0.5*IT_0199;
    const ccomplex_t IT_0212 = IT_0016*IT_0118;
    const ccomplex_t IT_0213 = IT_0081*IT_0212;
    const ccomplex_t IT_0214 = IT_0051*(IT_0050 + -IT_0070);
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = IT_0212*IT_0215;
    const ccomplex_t IT_0217 = (-10.6666666666667)*IT_0111;
    const ccomplex_t IT_0218 = IT_0130*IT_0217;
    const ccomplex_t IT_0219 = IT_0016*IT_0182;
    const ccomplex_t IT_0220 = (-0.5)*IT_0214;
    const ccomplex_t IT_0221 = IT_0103*IT_0136;
    const ccomplex_t IT_0222 = s_24*IT_0139;
    const ccomplex_t IT_0223 = IT_0023*IT_0222;
    const ccomplex_t IT_0224 = IT_0016*IT_0139;
    const ccomplex_t IT_0225 = IT_0037*IT_0224;
    const ccomplex_t IT_0226 = s_24*IT_0183;
    const ccomplex_t IT_0227 = IT_0039*IT_0226;
    const ccomplex_t IT_0228 = IT_0016*IT_0183;
    const ccomplex_t IT_0229 = IT_0027*IT_0228;
    const ccomplex_t IT_0230 = IT_0016*IT_0174;
    const ccomplex_t IT_0231 = (-0.666666666666667)*IT_0172;
    const ccomplex_t IT_0232 = conj(IT_0014)*IT_0231;
    const ccomplex_t IT_0233 = IT_0016*IT_0185;
    const ccomplex_t IT_0234 = (-0.666666666666667)*IT_0043;
    const ccomplex_t IT_0235 = IT_0018*(IT_0016*conj(IT_0027) + s_24*conj
      (IT_0039)) + IT_0123*IT_0157 + conj(IT_0010)*IT_0175 + IT_0010*IT_0190 +
       IT_0015*IT_0191 + IT_0192 + IT_0194 + IT_0013*IT_0195 + IT_0006*IT_0196 +
       IT_0198 + IT_0200*IT_0201 + IT_0202 + IT_0204 + IT_0117*IT_0205 + IT_0206
       + IT_0208 + IT_0210 + IT_0211*IT_0212 + IT_0213 + IT_0216 + IT_0218 +
       IT_0125*IT_0219 + IT_0135*IT_0220 + IT_0221 + IT_0223 + IT_0225 + IT_0227
       + IT_0229 + conj(IT_0015)*IT_0230 + IT_0232 + conj(IT_0013)*(IT_0230 +
       IT_0233) + conj(IT_0006)*IT_0234;
    const ccomplex_t IT_0236 = s_14*IT_0076;
    const ccomplex_t IT_0237 = IT_0053*IT_0236;
    const ccomplex_t IT_0238 = 1.5*IT_0097;
    const ccomplex_t IT_0239 = 1.5*IT_0103;
    const ccomplex_t IT_0240 = 1.5*IT_0006;
    const ccomplex_t IT_0241 = 1.5*conj(IT_0006);
    const ccomplex_t IT_0242 = IT_0016*IT_0179;
    const ccomplex_t IT_0243 = 1.5*IT_0013;
    const ccomplex_t IT_0244 = IT_0016*IT_0186;
    const ccomplex_t IT_0245 = 1.5*conj(IT_0013);
    const ccomplex_t IT_0246 = 1.5*IT_0015;
    const ccomplex_t IT_0247 = 1.5*conj(IT_0015);
    const ccomplex_t IT_0248 = s_34*IT_0177;
    const ccomplex_t IT_0249 = 1.5*IT_0014;
    const ccomplex_t IT_0250 = 1.5*conj(IT_0014);
    const ccomplex_t IT_0251 = IT_0016*IT_0154;
    const ccomplex_t IT_0252 = 1.5*IT_0037;
    const ccomplex_t IT_0253 = IT_0016*IT_0161;
    const ccomplex_t IT_0254 = 1.5*conj(IT_0037);
    const ccomplex_t IT_0255 = IT_0016*IT_0184;
    const ccomplex_t IT_0256 = 1.5*IT_0027;
    const ccomplex_t IT_0257 = IT_0016*IT_0181;
    const ccomplex_t IT_0258 = 1.5*conj(IT_0027);
    const ccomplex_t IT_0259 = s_24*IT_0154;
    const ccomplex_t IT_0260 = 1.5*IT_0023;
    const ccomplex_t IT_0261 = s_24*IT_0161;
    const ccomplex_t IT_0262 = 1.5*conj(IT_0023);
    const ccomplex_t IT_0263 = s_14*IT_0118;
    const ccomplex_t IT_0264 = 1.5*IT_0130;
    const ccomplex_t IT_0265 = -conj(IT_0013);
    const ccomplex_t IT_0266 = 1.5*IT_0117;
    const ccomplex_t IT_0267 = 1.5*IT_0125;
    const ccomplex_t IT_0268 = s_34*(conj(IT_0039)*IT_0041 + IT_0010*IT_0178 +
       conj(IT_0010)*IT_0185 + IT_0039*IT_0188) + IT_0237 + 0.666666666666667
      *IT_0101*IT_0238 + 0.666666666666667*IT_0134*IT_0239 + 0.666666666666667
      *IT_0195*IT_0240 + 0.666666666666667*IT_0233*IT_0241 + 0.666666666666667
      *IT_0242*IT_0243 + 0.666666666666667*IT_0244*IT_0245 + 0.666666666666667
      *IT_0190*IT_0246 + 0.666666666666667*IT_0175*IT_0247 + 0.666666666666667
      *IT_0248*IT_0249 + 0.666666666666667*IT_0176*IT_0250 + 0.666666666666667
      *IT_0251*IT_0252 + 0.666666666666667*IT_0253*IT_0254 + 0.666666666666667
      *IT_0255*IT_0256 + 0.666666666666667*IT_0257*IT_0258 + 0.666666666666667
      *IT_0259*IT_0260 + 0.666666666666667*IT_0261*IT_0262 + 0.666666666666667
      *IT_0263*IT_0264 + 0.666666666666667*s_24*(s_24*(IT_0013 + IT_0265) +
       IT_0180*IT_0266 + IT_0182*IT_0267);
    const ccomplex_t IT_0269 = 0.666666666666667*IT_0111;
    const ccomplex_t IT_0270 = (-0.666666666666667)*IT_0111;
    const ccomplex_t IT_0271 = s_14*IT_0054;
    const ccomplex_t IT_0272 = s_14*IT_0030;
    const ccomplex_t IT_0273 = 1.5*IT_0010;
    const ccomplex_t IT_0274 = 1.5*conj(IT_0010);
    const ccomplex_t IT_0275 = 1.5*IT_0039;
    const ccomplex_t IT_0276 = 1.5*conj(IT_0039);
    const ccomplex_t IT_0277 = 1.5*IT_0206;
    const ccomplex_t IT_0278 = 1.5*IT_0210;
    const ccomplex_t IT_0279 = 1.5*IT_0221;
    const ccomplex_t IT_0280 = 1.5*IT_0225;
    const ccomplex_t IT_0281 = (IT_0014 + -conj(IT_0014))*IT_0114 + IT_0144
      *IT_0246 + IT_0159*IT_0247 + IT_0140*IT_0256 + IT_0107*IT_0258 + IT_0251
      *IT_0260 + IT_0253*IT_0262 + IT_0157*IT_0264 + IT_0126*IT_0266 + IT_0119
      *IT_0267 + IT_0240*IT_0269 + IT_0241*IT_0270 + IT_0243*(IT_0144 + IT_0271)
       + IT_0245*(IT_0159 + IT_0272) + IT_0153*IT_0273 + IT_0163*IT_0274 +
       IT_0150*IT_0275 + IT_0167*IT_0276 + IT_0277 + IT_0278 + IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = 0.5*IT_0052;
    const ccomplex_t IT_0283 = IT_0077 + IT_0157;
    const ccomplex_t IT_0284 = (-10.6666666666667)*IT_0043;
    const ccomplex_t IT_0285 = s_24*IT_0149;
    const ccomplex_t IT_0286 = s_34*IT_0149;
    const ccomplex_t IT_0287 = s_34*IT_0028;
    const ccomplex_t IT_0288 = s_23*IT_0118;
    const ccomplex_t IT_0289 = (-0.5)*IT_0080;
    const ccomplex_t IT_0290 = 1.5*IT_0289;
    const ccomplex_t IT_0291 = 1.5*IT_0081;
    const ccomplex_t IT_0292 = IT_0053*(IT_0124 + IT_0132) + IT_0192 + IT_0194
       + IT_0198 + IT_0202 + IT_0204 + IT_0208 + IT_0213 + 2.66666666666667*s_24
      *(0.375*IT_0018*conj(IT_0039) + 0.375*IT_0082*IT_0125 + 0.375*IT_0117
      *IT_0127 + IT_0016*IT_0215) + IT_0216 + IT_0218 + IT_0135*(IT_0100 +
       IT_0220) + IT_0223 + IT_0227 + IT_0229 + IT_0232 + IT_0211*(IT_0212 +
       IT_0236) + 0.666666666666667*IT_0128*IT_0238 + 0.666666666666667*IT_0098
      *IT_0239 + 0.666666666666667*(IT_0196 + IT_0222)*IT_0240 +
       0.666666666666667*(IT_0207 + IT_0234)*IT_0241 + 0.666666666666667*
      (IT_0029 + IT_0230)*IT_0247 + 0.666666666666667*IT_0056*IT_0256 +
       0.666666666666667*(IT_0016*IT_0018 + IT_0033)*IT_0258 + 0.666666666666667
      *IT_0243*(IT_0145 + IT_0195 + IT_0259) + 0.666666666666667*IT_0245*
      (IT_0158 + IT_0230 + IT_0233 + IT_0261) + IT_0200*(IT_0201 + IT_0263) +
       0.666666666666667*IT_0156*IT_0264 + 0.666666666666667*IT_0205*IT_0266 +
       0.666666666666667*IT_0219*IT_0267 + 0.666666666666667*IT_0252*IT_0271 +
       0.666666666666667*IT_0254*IT_0272 + 0.666666666666667*(IT_0155 + IT_0190)
      *IT_0273 + 0.666666666666667*(IT_0162 + IT_0175)*IT_0274 +
       0.666666666666667*IT_0148*IT_0275 + 0.666666666666667*IT_0168*IT_0276 +
       0.666666666666667*IT_0277 + 0.666666666666667*IT_0278 + 0.666666666666667
      *IT_0279 + 0.666666666666667*IT_0280 + IT_0282*IT_0283 + IT_0123*(IT_0077 
      + IT_0284) + 0.666666666666667*IT_0246*(IT_0191 + IT_0285) +
       0.666666666666667*IT_0249*IT_0286 + 0.666666666666667*IT_0250*IT_0287 +
       0.666666666666667*(IT_0132 + IT_0203 + IT_0288)*IT_0290 +
       0.666666666666667*IT_0102*IT_0291 + 0.666666666666667*s_14*(s_14*IT_0079 
      + IT_0082*IT_0290 + IT_0127*IT_0291);
    const ccomplex_t IT_0293 = s_23*IT_0054;
    const ccomplex_t IT_0294 = IT_0016*IT_0082;
    const ccomplex_t IT_0295 = s_23*IT_0060;
    const ccomplex_t IT_0296 = s_23*IT_0149;
    const ccomplex_t IT_0297 = s_23*IT_0028;
    const ccomplex_t IT_0298 = s_23*IT_0128;
    const ccomplex_t IT_0299 = IT_0016*IT_0018*conj(IT_0037) + conj(IT_0006)
      *IT_0045 + IT_0057 + IT_0059 + IT_0065 + IT_0082*IT_0097 + IT_0105 +
       IT_0108 + IT_0110 + IT_0113 + IT_0115 + IT_0117*IT_0119 + IT_0125*IT_0126
       + IT_0103*IT_0127 + IT_0081*(IT_0119 + IT_0133) + IT_0130*(IT_0131 +
       IT_0133) + IT_0137 + IT_0141 + IT_0143 + IT_0147 + IT_0010*(IT_0062 +
       IT_0150) + IT_0152 + IT_0160 + conj(IT_0013)*(IT_0047 + IT_0164) + conj
      (IT_0010)*(IT_0025 + IT_0167) + IT_0170 + IT_0100*(IT_0126 + IT_0212) +
       IT_0200*(IT_0205 + IT_0212) + IT_0075*(IT_0197 + IT_0217) + IT_0211
      *IT_0219 + IT_0037*IT_0228 + IT_0023*IT_0255 + conj(IT_0023)*IT_0257 +
       IT_0039*(IT_0153 + IT_0286) + conj(IT_0039)*(IT_0163 + IT_0287) + IT_0123
      *(IT_0203 + IT_0288) + IT_0006*IT_0293 + IT_0289*(IT_0157 + IT_0294) +
       IT_0053*(s_14*IT_0127 + IT_0156 + IT_0284 + IT_0294) + IT_0013*(IT_0140 +
       IT_0295) + IT_0027*(IT_0144 + IT_0296) + conj(IT_0027)*(IT_0159 + IT_0297
      ) + IT_0220*(IT_0236 + IT_0298);
    const ccomplex_t IT_0300 = s_23*IT_0063;
    const ccomplex_t IT_0301 = IT_0130*IT_0201;
    const ccomplex_t IT_0302 = IT_0076*IT_0103;
    const ccomplex_t IT_0303 = conj(IT_0014)*IT_0297;
    const ccomplex_t IT_0304 = s_23*IT_0161;
    const ccomplex_t IT_0305 = s_23*IT_0106;
    const ccomplex_t IT_0306 = IT_0053*IT_0212;
    const ccomplex_t IT_0307 = IT_0097*IT_0118;
    const ccomplex_t IT_0308 = IT_0014*IT_0296;
    const ccomplex_t IT_0309 = IT_0016*IT_0149;
    const ccomplex_t IT_0310 = s_23*IT_0154;
    const ccomplex_t IT_0311 = s_23*IT_0139;
    const ccomplex_t IT_0312 = s_23*IT_0165;
    const ccomplex_t IT_0313 = IT_0006*IT_0224;
    const ccomplex_t IT_0314 = conj(IT_0006)*IT_0209;
    const ccomplex_t IT_0315 = conj(IT_0015)*IT_0040 + conj(IT_0027)*IT_0045 +
       IT_0117*IT_0133 + IT_0100*IT_0157 + IT_0077*IT_0220 + IT_0037*IT_0269 +
       conj(IT_0037)*IT_0270 + IT_0023*IT_0271 + conj(IT_0023)*IT_0272 + IT_0027
      *IT_0293 + IT_0125*IT_0294 + IT_0039*IT_0300 + IT_0301 + IT_0302 + IT_0303
       + conj(IT_0010)*IT_0304 + conj(IT_0013)*(IT_0253 + IT_0305) + IT_0306 +
       IT_0307 + IT_0308 + IT_0015*IT_0309 + IT_0010*IT_0310 + IT_0013*(IT_0251 
      + IT_0311) + conj(IT_0039)*IT_0312 + IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0016*IT_0060;
    const ccomplex_t IT_0317 = conj(IT_0010)*IT_0029 + IT_0034*conj(IT_0039) +
       conj(IT_0027)*IT_0048 + IT_0039*IT_0061 + IT_0023*IT_0195 + IT_0037
      *IT_0196 + IT_0125*IT_0197 + IT_0157*IT_0200 + conj(IT_0013)*(IT_0040 +
       IT_0207) + IT_0100*IT_0217 + IT_0156*IT_0220 + conj(IT_0023)*IT_0233 +
       conj(IT_0037)*IT_0234 + IT_0237 + IT_0010*IT_0285 + IT_0117*IT_0288 +
       IT_0301 + IT_0302 + IT_0303 + conj(IT_0015)*IT_0305 + IT_0306 + IT_0307 +
       IT_0308 + IT_0013*(IT_0222 + IT_0309) + IT_0015*IT_0311 + IT_0313 +
       IT_0314 + IT_0027*IT_0316;
    const ccomplex_t IT_0318 = conj(IT_0006)*IT_0042;
    const ccomplex_t IT_0319 = conj(IT_0006)*IT_0230;
    const ccomplex_t IT_0320 = s_34*IT_0183;
    const ccomplex_t IT_0321 = (-5.33333333333333)*IT_0172;
    const ccomplex_t IT_0322 = s_34*IT_0018*conj(IT_0039) + IT_0023*IT_0145 +
       conj(IT_0023)*IT_0158 + conj(IT_0027)*IT_0173 + conj(IT_0010)*IT_0176 +
       IT_0006*IT_0191 + IT_0013*(IT_0193 + IT_0242) + conj(IT_0013)*(IT_0231 +
       IT_0244) + IT_0010*IT_0248 + conj(IT_0037)*IT_0305 + IT_0037*IT_0311 +
       IT_0319 + IT_0039*IT_0320 + IT_0027*IT_0321;
    const ccomplex_t IT_0323 = conj(IT_0027)*IT_0230;
    const ccomplex_t IT_0324 = conj(IT_0013)*IT_0031;
    const ccomplex_t IT_0325 = conj(IT_0006)*IT_0046;
    const ccomplex_t IT_0326 = conj(IT_0014)*IT_0047;
    const ccomplex_t IT_0327 = conj(IT_0010)*IT_0034;
    const ccomplex_t IT_0328 = conj(IT_0013)*IT_0048;
    const ccomplex_t IT_0329 = conj(IT_0015)*IT_0045;
    const ccomplex_t IT_0330 = conj(IT_0027)*IT_0040;
    const ccomplex_t IT_0331 = conj(IT_0023)*IT_0042;
    const ccomplex_t IT_0332 = IT_0029*conj(IT_0039);
    const ccomplex_t IT_0333 = IT_0006*IT_0058;
    const ccomplex_t IT_0334 = IT_0014*IT_0295;
    const ccomplex_t IT_0335 = IT_0201*IT_0289;
    const ccomplex_t IT_0336 = IT_0053*IT_0201;
    const ccomplex_t IT_0337 = IT_0076*IT_0097;
    const ccomplex_t IT_0338 = IT_0053*IT_0263;
    const ccomplex_t IT_0339 = IT_0130*IT_0212;
    const ccomplex_t IT_0340 = IT_0103*IT_0118;
    const ccomplex_t IT_0341 = IT_0016*IT_0188;
    const ccomplex_t IT_0342 = conj(IT_0037)*IT_0044 + IT_0010*IT_0061 +
       IT_0117*IT_0197 + IT_0131*IT_0200 + IT_0100*IT_0203 + IT_0157*(IT_0037 +
       IT_0211) + IT_0132*IT_0220 + IT_0039*IT_0285 + IT_0125*IT_0288 + IT_0015
      *IT_0293 + IT_0027*IT_0309 + IT_0013*(IT_0055 + IT_0316) + IT_0324 +
       IT_0325 + IT_0326 + IT_0327 + IT_0328 + IT_0329 + IT_0330 + IT_0331 +
       IT_0332 + IT_0333 + IT_0334 + IT_0335 + IT_0336 + IT_0337 + IT_0338 +
       IT_0339 + IT_0340 + IT_0023*IT_0341;
    const ccomplex_t IT_0343 = conj(IT_0013)*IT_0173;
    const ccomplex_t IT_0344 = conj(IT_0013)*IT_0257;
    const ccomplex_t IT_0345 = conj(IT_0027)*IT_0231;
    const ccomplex_t IT_0346 = conj(IT_0023)*IT_0033;
    const ccomplex_t IT_0347 = conj(IT_0039)*IT_0176;
    const ccomplex_t IT_0348 = (s_34*conj(IT_0010) + conj(IT_0006)*IT_0016)
      *IT_0018 + conj(IT_0037)*IT_0045 + IT_0023*IT_0056 + IT_0027*IT_0193 +
       IT_0006*IT_0228 + IT_0039*IT_0248 + IT_0037*IT_0293 + IT_0010*IT_0320 +
       IT_0013*(IT_0255 + IT_0321) + IT_0343 + IT_0344 + IT_0345 + IT_0346 +
       IT_0347;
    const ccomplex_t IT_0349 = IT_0023*IT_0055;
    const ccomplex_t IT_0350 = IT_0039*IT_0190;
    const ccomplex_t IT_0351 = IT_0081*IT_0201;
    const ccomplex_t IT_0352 = IT_0201*IT_0215;
    const ccomplex_t IT_0353 = IT_0201*IT_0211;
    const ccomplex_t IT_0354 = IT_0130*IT_0203;
    const ccomplex_t IT_0355 = IT_0131*IT_0282;
    const ccomplex_t IT_0356 = conj(IT_0014)*IT_0173;
    const ccomplex_t IT_0357 = conj(IT_0013)*IT_0042;
    const ccomplex_t IT_0358 = conj(IT_0006)*IT_0044;
    const ccomplex_t IT_0359 = IT_0081*IT_0288;
    const ccomplex_t IT_0360 = IT_0124*IT_0282;
    const ccomplex_t IT_0361 = IT_0100*IT_0212;
    const ccomplex_t IT_0362 = IT_0200*IT_0212;
    const ccomplex_t IT_0363 = IT_0053*IT_0217;
    const ccomplex_t IT_0364 = IT_0135*IT_0215;
    const ccomplex_t IT_0365 = IT_0015*IT_0228;
    const ccomplex_t IT_0366 = IT_0014*IT_0321;
    const ccomplex_t IT_0367 = IT_0013*IT_0341;
    const ccomplex_t IT_0368 = conj(IT_0023)*IT_0031;
    const ccomplex_t IT_0369 = conj(IT_0039)*IT_0175;
    const ccomplex_t IT_0370 = IT_0013*IT_0228;
    const ccomplex_t IT_0371 = IT_0010*IT_0226;
    const ccomplex_t IT_0372 = (conj(IT_0013) + conj(IT_0015))*IT_0016;
    const ccomplex_t IT_0373 = s_24*conj(IT_0010) + IT_0372;
    const ccomplex_t IT_0374 = IT_0018*IT_0373;
    const ccomplex_t IT_0375 = IT_0027*IT_0191;
    const ccomplex_t IT_0376 = IT_0102*IT_0220;
    const ccomplex_t IT_0377 = IT_0117*IT_0219;
    const ccomplex_t IT_0378 = IT_0125*IT_0205;
    const ccomplex_t IT_0379 = IT_0059 + IT_0105 + IT_0123*IT_0131 + IT_0137 +
       (IT_0006 + IT_0068)*IT_0157 + IT_0160 + IT_0074*IT_0284 + IT_0323 +
       IT_0349 + IT_0350 + IT_0351 + IT_0352 + IT_0353 + IT_0354 + IT_0355 +
       IT_0356 + IT_0357 + IT_0358 + IT_0359 + IT_0360 + IT_0361 + IT_0362 +
       IT_0363 + IT_0364 + IT_0365 + IT_0366 + IT_0367 + IT_0368 + IT_0369 +
       IT_0370 + IT_0371 + IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = IT_0059 + IT_0105 + IT_0123*IT_0131 + IT_0137 +
       (IT_0006 + IT_0074)*IT_0157 + IT_0160 + IT_0323 + IT_0349 + IT_0350 +
       IT_0351 + IT_0352 + IT_0353 + IT_0354 + IT_0355 + IT_0356 + IT_0357 +
       IT_0358 + IT_0359 + IT_0360 + IT_0361 + IT_0362 + IT_0363 + IT_0364 +
       IT_0365 + IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370 + IT_0371 +
       IT_0374 + IT_0375 + IT_0376 + IT_0377 + IT_0378;
    const ccomplex_t IT_0381 = IT_0013*IT_0293;
    const ccomplex_t IT_0382 = IT_0015*IT_0316;
    const ccomplex_t IT_0383 = conj(IT_0023)*IT_0109;
    const ccomplex_t IT_0384 = conj(IT_0027)*IT_0305;
    const ccomplex_t IT_0385 = conj(IT_0037)*IT_0112;
    const ccomplex_t IT_0386 = IT_0023*IT_0142;
    const ccomplex_t IT_0387 = IT_0027*IT_0311;
    const ccomplex_t IT_0388 = IT_0037*IT_0146;
    const ccomplex_t IT_0389 = conj(IT_0013)*IT_0045;
    const ccomplex_t IT_0390 = conj(IT_0015)*IT_0048;
    const ccomplex_t IT_0391 = IT_0125*IT_0133;
    const ccomplex_t IT_0392 = conj(IT_0027)*IT_0047 + IT_0082*IT_0103 +
       IT_0097*IT_0127 + IT_0100*(IT_0119 + IT_0131) + IT_0081*(IT_0126 +
       IT_0157) + IT_0037*IT_0191 + IT_0211*(IT_0205 + IT_0217) + IT_0200*
      (IT_0203 + IT_0219) + (s_23*IT_0098 + IT_0124)*IT_0220 + conj(IT_0037)
      *IT_0230 + IT_0023*IT_0242 + conj(IT_0023)*IT_0244 + (IT_0117 + IT_0130)
      *IT_0294 + IT_0027*IT_0295 + IT_0013*(IT_0064 + IT_0145 + IT_0296) + conj
      (IT_0013)*(IT_0158 + IT_0169 + IT_0297) + IT_0010*(IT_0286 + IT_0300) +
       conj(IT_0039)*(IT_0025 + IT_0304) + conj(IT_0006)*IT_0305 + IT_0039*
      (IT_0062 + IT_0310) + IT_0006*IT_0311 + conj(IT_0010)*(IT_0287 + IT_0312) 
      + IT_0325 + IT_0326 + IT_0333 + IT_0334 + IT_0335 + IT_0336 + IT_0337 +
       IT_0339 + IT_0340 + IT_0381 + IT_0382 + IT_0383 + IT_0384 + IT_0385 +
       IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = IT_0037*IT_0064 + IT_0101*IT_0103 + IT_0102
      *IT_0123 + IT_0100*(IT_0124 + IT_0131) + IT_0097*IT_0134 + IT_0053*IT_0135
       + IT_0023*IT_0138 + conj(IT_0023)*IT_0166 + conj(IT_0037)*IT_0169 + s_34*
      (conj(IT_0014)*IT_0018 + conj(IT_0010)*IT_0041 + IT_0039*IT_0178 + conj
      (IT_0039)*IT_0185 + IT_0010*IT_0188) + IT_0200*(IT_0132 + IT_0203) +
       IT_0211*(IT_0156 + IT_0217) + IT_0124*IT_0220 + IT_0015*IT_0226 + IT_0130
      *IT_0236 + IT_0027*IT_0242 + conj(IT_0027)*IT_0244 + IT_0013*(IT_0064 +
       IT_0255) + conj(IT_0013)*(IT_0169 + IT_0257) + s_24*(conj(IT_0015)
      *IT_0018 + IT_0125*IT_0180 + IT_0117*IT_0182 + (-5.33333333333333)*s_24*
      (IT_0013 + IT_0265)) + IT_0215*IT_0283 + (s_23*IT_0098 + IT_0119 + IT_0263
      )*IT_0289 + IT_0117*IT_0294 + IT_0081*(IT_0077 + IT_0126 + IT_0284 +
       IT_0298) + IT_0010*IT_0300 + conj(IT_0039)*IT_0304 + IT_0039*IT_0310 +
       conj(IT_0010)*IT_0312 + IT_0318 + IT_0014*IT_0320 + IT_0325 + IT_0326 +
       IT_0333 + IT_0334 + IT_0335 + IT_0336 + IT_0337 + IT_0338 + IT_0339 +
       IT_0340 + IT_0006*IT_0341 + IT_0381 + IT_0382 + IT_0383 + IT_0384 +
       IT_0385 + IT_0386 + IT_0387 + IT_0388 + IT_0389 + IT_0390 + IT_0391;
    const ccomplex_t IT_0394 = IT_0019 + IT_0014*conj(IT_0023)*IT_0025 +
       IT_0023*(conj(IT_0027)*IT_0029 + IT_0032 + conj(IT_0015)*IT_0033 + conj
      (IT_0013)*IT_0034 + IT_0035) + conj(IT_0037)*(IT_0010*IT_0034 + IT_0029
      *IT_0039 + IT_0027*IT_0040 + IT_0023*IT_0042 + IT_0037*IT_0044 + IT_0015
      *IT_0045 + IT_0006*IT_0046 + IT_0014*IT_0047 + IT_0013*(IT_0031 + IT_0048)
      ) + IT_0053*IT_0171 + conj(IT_0006)*(IT_0014*IT_0173 + IT_0039*IT_0175) +
       conj(IT_0015)*IT_0039*IT_0176 + IT_0103*IT_0187 + IT_0097*IT_0189 +
       IT_0123*IT_0235 + IT_0220*IT_0268 + 0.666666666666667*IT_0130*IT_0281 +
       IT_0282*IT_0292 + IT_0289*IT_0299 + IT_0100*IT_0315 + IT_0200*IT_0317 +
       IT_0013*(10.6666666666667*s_23*s_24*conj(IT_0013) + conj(IT_0023)*IT_0034
       + conj(IT_0027)*IT_0175 + conj(IT_0015)*(IT_0173 + IT_0257) + IT_0318) +
       IT_0027*(conj(IT_0023)*IT_0029 + conj(IT_0013)*IT_0175 + conj(IT_0015)
      *IT_0231 + IT_0319) + IT_0125*IT_0322 + IT_0006*(conj(IT_0023)*IT_0031 +
       conj(IT_0013)*IT_0042 + conj(IT_0006)*IT_0044 + conj(IT_0014)*IT_0173 +
       conj(IT_0039)*IT_0175 + IT_0323) + IT_0037*(IT_0324 + IT_0325 + IT_0326 +
       IT_0327 + IT_0328 + IT_0329 + IT_0330 + IT_0331 + IT_0332) + IT_0211
      *IT_0342 + IT_0015*(IT_0343 + IT_0344 + IT_0345 + IT_0346 + IT_0347) +
       IT_0117*IT_0348 + IT_0068*IT_0379 + IT_0074*IT_0380 + IT_0081*IT_0392 +
       IT_0215*IT_0393;
    return create_ccomplex_return(IT_0394);
}

