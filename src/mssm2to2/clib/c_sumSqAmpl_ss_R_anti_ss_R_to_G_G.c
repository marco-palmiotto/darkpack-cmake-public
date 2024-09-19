#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_ss_R_anti_ss_R_to_G_G.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_ss_R_anti_ss_R_to_G_G(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t Gamma_sr = param->Gamma_sr;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*g_s;
    const ccomplex_t IT_0001 = -IT_0000;
    const ccomplex_t IT_0002 = (-2)*IT_0000;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0001, 2);
    const ccomplex_t IT_0005 = cpow(s_23 + (0 + _Complex_I*(-0.5))*m_ss_R
      *Gamma_sr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*IT_0005;
    const ccomplex_t IT_0007 = (IT_0003 + -IT_0004)*IT_0006;
    const ccomplex_t IT_0008 = 0.5*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0002, 2);
    const ccomplex_t IT_0010 = IT_0006*(IT_0003 + -IT_0009);
    const ccomplex_t IT_0011 = (-0.5)*IT_0010;
    const ccomplex_t IT_0012 = pow(g_s, 2);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = pow(m_ss_R, 2);
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
    const ccomplex_t IT_0034 = IT_0013 + (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = 5.33333333333333*s_23;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0013 + 0.5*IT_0033;
    const ccomplex_t IT_0038 = (-0.666666666666667)*s_23;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = IT_0036 + IT_0039;
    const ccomplex_t IT_0041 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_ss_R
      *Gamma_sr + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = (IT_0003 + -IT_0004)*IT_0042;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0003*IT_0042;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = 5.33333333333333*s_13;
    const ccomplex_t IT_0048 = IT_0014*IT_0047;
    const ccomplex_t IT_0049 = IT_0046*IT_0048;
    const ccomplex_t IT_0050 = IT_0037*IT_0047;
    const ccomplex_t IT_0051 = 0.5*IT_0010;
    const ccomplex_t IT_0052 = (-0.666666666666667)*s_12;
    const ccomplex_t IT_0053 = s_23*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = (-0.666666666666667)*s_13;
    const ccomplex_t IT_0056 = IT_0034*IT_0055;
    const ccomplex_t IT_0057 = (-0.5)*IT_0007;
    const ccomplex_t IT_0058 = IT_0003*IT_0006;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = s_23*IT_0047;
    const ccomplex_t IT_0061 = pow(s_23, 2);
    const ccomplex_t IT_0062 = 5.33333333333333*IT_0061;
    const ccomplex_t IT_0063 = (IT_0003 + -IT_0009)*IT_0042;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = s_23*IT_0055;
    const ccomplex_t IT_0066 = pow(s_13, 2);
    const ccomplex_t IT_0067 = (-0.666666666666667)*IT_0066;
    const ccomplex_t IT_0068 = 0.5*IT_0045;
    const ccomplex_t IT_0069 = (IT_0003 + (-0.5)*IT_0004 + (-0.5)*IT_0009)
      *IT_0042;
    const ccomplex_t IT_0070 = g_s*IT_0002;
    const ccomplex_t IT_0071 = g_s*IT_0001;
    const ccomplex_t IT_0072 = IT_0015*(IT_0070 + -IT_0071);
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0015*IT_0071;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = 5.33333333333333*IT_0066;
    const ccomplex_t IT_0077 = IT_0004*IT_0006;
    const ccomplex_t IT_0078 = (-0.5)*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0074;
    const ccomplex_t IT_0080 = (-0.5)*IT_0043;
    const ccomplex_t IT_0081 = (-0.666666666666667)*IT_0061;
    const ccomplex_t IT_0082 = 5.33333333333333*s_12;
    const ccomplex_t IT_0083 = s_13*IT_0082;
    const ccomplex_t IT_0084 = conj(IT_0046)*IT_0083;
    const ccomplex_t IT_0085 = IT_0048*conj(IT_0064);
    const ccomplex_t IT_0086 = IT_0053*conj(IT_0059);
    const ccomplex_t IT_0087 = IT_0014*IT_0038;
    const ccomplex_t IT_0088 = conj(IT_0051)*IT_0087;
    const ccomplex_t IT_0089 = conj(IT_0011)*IT_0081;
    const ccomplex_t IT_0090 = 0.666666666666667*s_23;
    const ccomplex_t IT_0091 = IT_0014*IT_0090;
    const ccomplex_t IT_0092 = IT_0051*IT_0091;
    const ccomplex_t IT_0093 = 1.33333333333333*s_23;
    const ccomplex_t IT_0094 = IT_0034*IT_0093;
    const ccomplex_t IT_0095 = IT_0015*IT_0025;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = (-10.6666666666667)*s_12;
    const ccomplex_t IT_0098 = s_13*IT_0097;
    const ccomplex_t IT_0099 = IT_0096*IT_0098;
    const ccomplex_t IT_0100 = (-10.6666666666667)*s_23;
    const ccomplex_t IT_0101 = IT_0037*IT_0100;
    const ccomplex_t IT_0102 = (-5.33333333333333)*s_12;
    const ccomplex_t IT_0103 = s_13*IT_0102;
    const ccomplex_t IT_0104 = IT_0046*IT_0103;
    const ccomplex_t IT_0105 = IT_0015*(IT_0025 + -IT_0070);
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = 0.666666666666667*s_12;
    const ccomplex_t IT_0108 = s_23*IT_0107;
    const ccomplex_t IT_0109 = IT_0014*IT_0107;
    const ccomplex_t IT_0110 = IT_0051*IT_0109;
    const ccomplex_t IT_0111 = 0.666666666666667*s_13;
    const ccomplex_t IT_0112 = IT_0014*IT_0111;
    const ccomplex_t IT_0113 = pow(m_ss_R, 4);
    const ccomplex_t IT_0114 = 0.666666666666667*IT_0113;
    const ccomplex_t IT_0115 = 1.33333333333333*s_12;
    const ccomplex_t IT_0116 = s_13*IT_0115;
    const ccomplex_t IT_0117 = IT_0096*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0072;
    const ccomplex_t IT_0119 = s_23*IT_0115;
    const ccomplex_t IT_0120 = IT_0034*IT_0115;
    const ccomplex_t IT_0121 = -IT_0015*(IT_0020 + -IT_0025 + -IT_0070);
    const ccomplex_t IT_0122 = 0.5*IT_0121;
    const ccomplex_t IT_0123 = pow(s_12, 2);
    const ccomplex_t IT_0124 = 1.33333333333333*IT_0123;
    const ccomplex_t IT_0125 = IT_0014*IT_0093;
    const ccomplex_t IT_0126 = 0.5*IT_0105;
    const ccomplex_t IT_0127 = 1.33333333333333*IT_0113;
    const ccomplex_t IT_0128 = s_23*IT_0082;
    const ccomplex_t IT_0129 = IT_0014*IT_0082;
    const ccomplex_t IT_0130 = conj(IT_0046)*IT_0129;
    const ccomplex_t IT_0131 = 5.33333333333333*IT_0123;
    const ccomplex_t IT_0132 = IT_0060*conj(IT_0080);
    const ccomplex_t IT_0133 = IT_0014*IT_0035;
    const ccomplex_t IT_0134 = 0.5*IT_0095;
    const ccomplex_t IT_0135 = IT_0098*IT_0134;
    const ccomplex_t IT_0136 = s_23*IT_0097;
    const ccomplex_t IT_0137 = IT_0014*IT_0097;
    const ccomplex_t IT_0138 = IT_0096*IT_0137;
    const ccomplex_t IT_0139 = IT_0037*IT_0097;
    const ccomplex_t IT_0140 = (-0.5)*IT_0121;
    const ccomplex_t IT_0141 = (-10.6666666666667)*IT_0123;
    const ccomplex_t IT_0142 = IT_0014*IT_0100;
    const ccomplex_t IT_0143 = s_23*IT_0102;
    const ccomplex_t IT_0144 = IT_0014*IT_0102;
    const ccomplex_t IT_0145 = IT_0046*IT_0144;
    const ccomplex_t IT_0146 = (-5.33333333333333)*IT_0123;
    const ccomplex_t IT_0147 = (-5.33333333333333)*s_13;
    const ccomplex_t IT_0148 = s_23*IT_0147;
    const ccomplex_t IT_0149 = IT_0080*IT_0148;
    const ccomplex_t IT_0150 = (-5.33333333333333)*s_23;
    const ccomplex_t IT_0151 = IT_0014*IT_0150;
    const ccomplex_t IT_0152 = (-5.33333333333333)*IT_0113;
    const ccomplex_t IT_0153 = IT_0014*IT_0052;
    const ccomplex_t IT_0154 = conj(IT_0051)*IT_0153;
    const ccomplex_t IT_0155 = conj(IT_0057)*IT_0065;
    const ccomplex_t IT_0156 = IT_0014*IT_0055;
    const ccomplex_t IT_0157 = (-0.666666666666667)*IT_0113;
    const ccomplex_t IT_0158 = conj(IT_0008)*IT_0053 + conj(IT_0068)*IT_0083 +
       conj(IT_0078)*IT_0087 + IT_0078*IT_0091 + IT_0068*IT_0103 + IT_0008
      *IT_0108 + IT_0110 + s_23*IT_0057*IT_0111 + IT_0011*IT_0112 + IT_0059
      *IT_0114 + IT_0117 + IT_0118*IT_0119 + IT_0120 + IT_0122*IT_0124 + IT_0079
      *IT_0125 + IT_0126*IT_0127 + conj(IT_0069)*IT_0128 + IT_0130 + conj
      (IT_0064)*IT_0131 + IT_0132 + conj(IT_0044)*IT_0133 + IT_0135 + IT_0073
      *IT_0136 + IT_0138 + IT_0139 + IT_0140*IT_0141 + IT_0075*IT_0142 + IT_0069
      *IT_0143 + IT_0145 + IT_0064*IT_0146 + IT_0149 + IT_0044*IT_0151 + IT_0106
      *IT_0152 + IT_0154 + IT_0155 + conj(IT_0011)*IT_0156 + conj(IT_0059)
      *IT_0157;
    const ccomplex_t IT_0159 = IT_0051*IT_0108;
    const ccomplex_t IT_0160 = conj(IT_0046)*IT_0048;
    const ccomplex_t IT_0161 = (-10.6666666666667)*s_13*IT_0037;
    const ccomplex_t IT_0162 = IT_0014*IT_0147;
    const ccomplex_t IT_0163 = IT_0046*IT_0162;
    const ccomplex_t IT_0164 = conj(IT_0051)*IT_0053;
    const ccomplex_t IT_0165 = conj(IT_0059)*IT_0087;
    const ccomplex_t IT_0166 = conj(IT_0011)*IT_0065;
    const ccomplex_t IT_0167 = IT_0037*IT_0147;
    const ccomplex_t IT_0168 = 0.666666666666667*IT_0066;
    const ccomplex_t IT_0169 = IT_0060*conj(IT_0069);
    const ccomplex_t IT_0170 = conj(IT_0044)*IT_0076;
    const ccomplex_t IT_0171 = (-5.33333333333333)*IT_0066;
    const ccomplex_t IT_0172 = s_13*IT_0052;
    const ccomplex_t IT_0173 = conj(IT_0059)*IT_0172;
    const ccomplex_t IT_0174 = conj(IT_0008)*IT_0065;
    const ccomplex_t IT_0175 = IT_0067*conj(IT_0078);
    const ccomplex_t IT_0176 = conj(IT_0064)*IT_0133;
    const ccomplex_t IT_0177 = s_13*IT_0107;
    const ccomplex_t IT_0178 = IT_0148 + IT_0177;
    const ccomplex_t IT_0179 = IT_0112 + IT_0171;
    const ccomplex_t IT_0180 = 1.33333333333333*s_13;
    const ccomplex_t IT_0181 = IT_0014*IT_0180;
    const ccomplex_t IT_0182 = (-10.6666666666667)*IT_0066;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = IT_0152 + IT_0181;
    const ccomplex_t IT_0185 = IT_0064*IT_0109;
    const ccomplex_t IT_0186 = IT_0014*IT_0115;
    const ccomplex_t IT_0187 = 1.33333333333333*IT_0066;
    const ccomplex_t IT_0188 = 1.33333333333333*IT_0014;
    const ccomplex_t IT_0189 = IT_0037*IT_0188;
    const ccomplex_t IT_0190 = conj(IT_0008) + conj(IT_0011);
    const ccomplex_t IT_0191 = conj(IT_0059)*IT_0129;
    const ccomplex_t IT_0192 = (-10.6666666666667)*IT_0014;
    const ccomplex_t IT_0193 = conj(IT_0064)*IT_0153;
    const ccomplex_t IT_0194 = conj(IT_0057)*IT_0076 + IT_0048*conj(IT_0078) +
       IT_0067*conj(IT_0080) + s_13*((-10.6666666666667)*s_23*IT_0073 + IT_0014*
      ((-10.6666666666667)*IT_0079 + (-21.3333333333333)*IT_0096)) + (IT_0008 +
       IT_0011)*IT_0103 + IT_0068*IT_0112 + IT_0046*IT_0114 + IT_0073*IT_0116 +
       IT_0096*IT_0127 + conj(IT_0051)*IT_0131 + IT_0126*(IT_0116 + IT_0137) +
       IT_0122*(IT_0125 + IT_0137) + IT_0140*IT_0142 + IT_0051*IT_0146 + IT_0064
      *IT_0151 + (conj(IT_0044) + conj(IT_0068))*IT_0156 + conj(IT_0046)*IT_0157
       + IT_0159 + IT_0160 + IT_0161 + IT_0163 + IT_0164 + IT_0080*IT_0168 +
       IT_0078*(IT_0162 + IT_0168) + IT_0169 + IT_0170 + IT_0057*IT_0171 + conj
      (IT_0069)*IT_0172 + IT_0173 + IT_0174 + IT_0175 + IT_0176 + IT_0059*
      (IT_0144 + IT_0177) + IT_0069*IT_0178 + IT_0044*IT_0179 + s_23*(IT_0008
      *IT_0111 + IT_0118*IT_0180) + IT_0075*IT_0183 + IT_0134*IT_0184 + IT_0185 
      + (IT_0106 + IT_0140)*IT_0186 + IT_0079*IT_0187 + IT_0189 + IT_0083
      *IT_0190 + IT_0191 + IT_0034*(IT_0180 + IT_0192) + IT_0193;
    const ccomplex_t IT_0195 = conj(IT_0059)*IT_0133;
    const ccomplex_t IT_0196 = conj(IT_0064)*IT_0172;
    const ccomplex_t IT_0197 = conj(IT_0051)*IT_0128;
    const ccomplex_t IT_0198 = IT_0051*IT_0143;
    const ccomplex_t IT_0199 = conj(IT_0046)*IT_0156;
    const ccomplex_t IT_0200 = 5.33333333333333*IT_0113;
    const ccomplex_t IT_0201 = conj(IT_0078)*IT_0172;
    const ccomplex_t IT_0202 = conj(IT_0011)*IT_0053;
    const ccomplex_t IT_0203 = (-0.666666666666667)*IT_0123;
    const ccomplex_t IT_0204 = conj(IT_0059)*IT_0203;
    const ccomplex_t IT_0205 = conj(IT_0008)*IT_0156;
    const ccomplex_t IT_0206 = conj(IT_0044)*IT_0083;
    const ccomplex_t IT_0207 = IT_0048*conj(IT_0069);
    const ccomplex_t IT_0208 = IT_0060*conj(IT_0080) + IT_0037*IT_0102 +
       IT_0075*IT_0103 + IT_0112*IT_0118 + IT_0130 + conj(IT_0068)*IT_0133 +
       IT_0154 + IT_0155 + IT_0073*IT_0162 + IT_0079*IT_0177 + conj(IT_0064)
      *IT_0200 + IT_0201 + IT_0202 + IT_0204 + IT_0205 + IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = 0.666666666666667*IT_0123;
    const ccomplex_t IT_0210 = (-10.6666666666667)*s_13;
    const ccomplex_t IT_0211 = IT_0075*IT_0098 + IT_0044*IT_0103 + IT_0011
      *IT_0108 + IT_0110 + s_23*IT_0057*IT_0111 + IT_0008*IT_0112 + IT_0079
      *IT_0116 + IT_0120 + IT_0122*IT_0127 + IT_0130 + IT_0132 + conj(IT_0068)
      *IT_0133 + IT_0138 + IT_0139 + IT_0145 + IT_0149 + IT_0068*IT_0151 + 
      (IT_0064 + IT_0140)*IT_0152 + IT_0154 + IT_0155 + IT_0069*IT_0162 +
       IT_0078*IT_0177 + IT_0118*IT_0181 + conj(IT_0064)*IT_0200 + IT_0201 +
       IT_0202 + IT_0204 + IT_0205 + IT_0206 + IT_0207 + IT_0059*IT_0209 +
       IT_0014*IT_0073*IT_0210;
    const ccomplex_t IT_0212 = conj(IT_0059)*IT_0200;
    const ccomplex_t IT_0213 = conj(IT_0011)*IT_0048;
    const ccomplex_t IT_0214 = conj(IT_0044)*IT_0087;
    const ccomplex_t IT_0215 = conj(IT_0064)*IT_0203;
    const ccomplex_t IT_0216 = IT_0053*conj(IT_0069);
    const ccomplex_t IT_0217 = conj(IT_0068)*IT_0172;
    const ccomplex_t IT_0218 = conj(IT_0057)*IT_0060 + IT_0065*conj(IT_0080) +
       IT_0075*IT_0091 + IT_0037*IT_0107 + IT_0073*IT_0108 + conj(IT_0008)
      *IT_0128 + conj(IT_0051)*IT_0129 + conj(IT_0078)*IT_0133 + IT_0118*IT_0143
       + IT_0079*IT_0151 + conj(IT_0046)*IT_0153 + IT_0212 + IT_0213 + IT_0214 +
       IT_0215 + IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = conj(IT_0011)*IT_0172;
    const ccomplex_t IT_0220 = conj(IT_0059)*IT_0153;
    const ccomplex_t IT_0221 = conj(IT_0051)*IT_0203;
    const ccomplex_t IT_0222 = conj(IT_0057)*IT_0067;
    const ccomplex_t IT_0223 = IT_0076*conj(IT_0080);
    const ccomplex_t IT_0224 = conj(IT_0064)*IT_0129;
    const ccomplex_t IT_0225 = IT_0037*IT_0192;
    const ccomplex_t IT_0226 = IT_0064*IT_0144;
    const ccomplex_t IT_0227 = IT_0044 + IT_0068;
    const ccomplex_t IT_0228 = conj(IT_0044)*(IT_0048 + IT_0067) + IT_0048
      *conj(IT_0068) + IT_0065*conj(IT_0069) + (conj(IT_0059) + conj(IT_0069))
      *IT_0083 + conj(IT_0064)*IT_0087 + IT_0064*IT_0091 + IT_0069*IT_0103 +
       IT_0059*(IT_0103 + IT_0109) + s_13*((-10.6666666666667)*IT_0034 + (
      -10.6666666666667)*IT_0014*IT_0075 + (-10.6666666666667)*s_23*IT_0118) +
       IT_0125*IT_0140 + IT_0122*IT_0142 + IT_0046*(IT_0112 + IT_0152) + conj
      (IT_0078)*(IT_0076 + IT_0156) + (IT_0044 + IT_0057)*IT_0168 + IT_0080
      *IT_0171 + conj(IT_0008)*(IT_0060 + IT_0172) + IT_0011*IT_0177 + IT_0008
      *IT_0178 + IT_0078*IT_0179 + IT_0037*IT_0180 + s_23*(IT_0069*IT_0111 +
       IT_0073*IT_0180) + IT_0079*IT_0183 + IT_0096*IT_0184 + IT_0075*IT_0187 +
       IT_0034*IT_0188 + IT_0197 + IT_0198 + IT_0199 + conj(IT_0046)*IT_0200 +
       IT_0051*IT_0209 + IT_0219 + IT_0220 + IT_0221 + IT_0222 + IT_0223 +
       IT_0224 + IT_0225 + IT_0226 + IT_0162*IT_0227;
    const ccomplex_t IT_0229 = -IT_0015*(IT_0020 + -IT_0025);
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0108 + IT_0148;
    const ccomplex_t IT_0232 = (-5.33333333333333)*IT_0061;
    const ccomplex_t IT_0233 = IT_0091 + IT_0232;
    const ccomplex_t IT_0234 = (-10.6666666666667)*IT_0061;
    const ccomplex_t IT_0235 = IT_0125 + IT_0234;
    const ccomplex_t IT_0236 = IT_0125 + IT_0152;
    const ccomplex_t IT_0237 = 1.33333333333333*IT_0061;
    const ccomplex_t IT_0238 = IT_0142 + IT_0237;
    const ccomplex_t IT_0239 = IT_0065 + IT_0128;
    const ccomplex_t IT_0240 = IT_0053 + IT_0060;
    const ccomplex_t IT_0241 = IT_0081 + IT_0133;
    const ccomplex_t IT_0242 = IT_0046*IT_0177;
    const ccomplex_t IT_0243 = 0.666666666666667*IT_0061;
    const ccomplex_t IT_0244 = IT_0079 + IT_0096 + IT_0140;
    const ccomplex_t IT_0245 = IT_0037*IT_0093;
    const ccomplex_t IT_0246 = IT_0048*conj(IT_0059);
    const ccomplex_t IT_0247 = conj(IT_0008)*IT_0062;
    const ccomplex_t IT_0248 = IT_0034*IT_0100;
    const ccomplex_t IT_0249 = IT_0051*IT_0151;
    const ccomplex_t IT_0250 = conj(IT_0046)*IT_0172;
    const ccomplex_t IT_0251 = IT_0053*conj(IT_0064);
    const ccomplex_t IT_0252 = IT_0075*IT_0180;
    const ccomplex_t IT_0253 = IT_0062*conj(IT_0080) + conj(IT_0057)*IT_0081 +
       IT_0011*IT_0091 + IT_0064*IT_0108 + IT_0051*IT_0114 + IT_0117 + conj
      (IT_0068)*IT_0128 + conj(IT_0046)*IT_0131 + IT_0124*IT_0134 + IT_0135 + 
      (IT_0075 + IT_0122 + IT_0134)*IT_0136 + IT_0137*IT_0140 + IT_0096*IT_0141 
      + IT_0046*IT_0146 + IT_0069*IT_0151 + conj(IT_0051)*(IT_0133 + IT_0157) +
       IT_0059*(IT_0109 + IT_0162) + IT_0106*(IT_0137 + IT_0181) + IT_0122
      *IT_0186 + IT_0034*IT_0188 + IT_0087*IT_0190 + IT_0126*(IT_0186 + IT_0014
      *IT_0210) + IT_0220 + IT_0224 + IT_0225 + IT_0226 + IT_0143*IT_0227 +
       IT_0078*IT_0231 + IT_0080*IT_0232 + IT_0008*IT_0233 + IT_0118*IT_0235 +
       IT_0230*IT_0236 + IT_0073*IT_0238 + conj(IT_0044)*IT_0239 + conj(IT_0078)
      *IT_0240 + conj(IT_0069)*IT_0241 + IT_0242 + (IT_0057 + IT_0069)*IT_0243 +
       IT_0119*IT_0244 + IT_0245 + IT_0246 + IT_0247 + IT_0248 + IT_0249 +
       IT_0250 + IT_0251 + s_23*(IT_0044*IT_0111 + IT_0079*IT_0210 + IT_0252);
    const ccomplex_t IT_0254 = 0.666666666666667*IT_0014;
    const ccomplex_t IT_0255 = 5.33333333333333*IT_0014;
    const ccomplex_t IT_0256 = (-5.33333333333333)*IT_0014;
    const ccomplex_t IT_0257 = (-0.666666666666667)*IT_0014;
    const ccomplex_t IT_0258 = (-21.3333333333333)*IT_0034 + 5.33333333333333
      *IT_0037 + IT_0052*conj(IT_0064) + IT_0038*conj(IT_0069) + IT_0047*conj
      (IT_0078) + conj(IT_0059)*IT_0082 + IT_0069*IT_0090 + IT_0059*IT_0102 +
       IT_0064*IT_0107 + IT_0078*IT_0147 + IT_0079*IT_0210 + IT_0111*IT_0227 +
       IT_0252 + IT_0046*IT_0254 + conj(IT_0051)*IT_0255 + IT_0051*IT_0256 +
       conj(IT_0046)*IT_0257;
    const ccomplex_t IT_0259 = IT_0078*IT_0111;
    const ccomplex_t IT_0260 = IT_0079*IT_0180;
    const ccomplex_t IT_0261 = (-21.3333333333333)*IT_0037 + IT_0052*conj
      (IT_0059) + IT_0035*conj(IT_0069) + IT_0055*conj(IT_0078) + conj(IT_0064)
      *IT_0082 + IT_0069*IT_0150 + IT_0051*IT_0254 + conj(IT_0046)*IT_0255 +
       IT_0046*IT_0256 + conj(IT_0051)*IT_0257 + IT_0259 + IT_0260;
    const ccomplex_t IT_0262 = IT_0046*IT_0109;
    const ccomplex_t IT_0263 = IT_0096*IT_0186;
    const ccomplex_t IT_0264 = conj(IT_0051)*IT_0129;
    const ccomplex_t IT_0265 = conj(IT_0057)*IT_0060;
    const ccomplex_t IT_0266 = IT_0051*IT_0144;
    const ccomplex_t IT_0267 = IT_0057*IT_0148;
    const ccomplex_t IT_0268 = conj(IT_0046)*IT_0153;
    const ccomplex_t IT_0269 = IT_0065*conj(IT_0080);
    const ccomplex_t IT_0270 = s_23*IT_0080*IT_0111;
    const ccomplex_t IT_0271 = IT_0034*IT_0097;
    const ccomplex_t IT_0272 = IT_0037*IT_0115;
    const ccomplex_t IT_0273 = conj(IT_0008)*IT_0048 + conj(IT_0078)*IT_0083 +
       conj(IT_0068)*IT_0087 + IT_0068*IT_0091 + IT_0079*IT_0098 + IT_0078
      *IT_0103 + IT_0069*IT_0112 + IT_0064*IT_0114 + IT_0075*IT_0116 + conj
      (IT_0011)*IT_0128 + conj(IT_0059)*IT_0131 + IT_0011*IT_0143 + IT_0059
      *IT_0146 + IT_0122*IT_0152 + conj(IT_0069)*IT_0156 + conj(IT_0064)*IT_0157
       + IT_0008*IT_0162 + conj(IT_0044)*IT_0172 + IT_0044*IT_0177 + IT_0073
      *IT_0181 + IT_0014*IT_0118*IT_0210 + IT_0262 + IT_0263 + IT_0264 + IT_0265
       + IT_0266 + IT_0267 + IT_0268 + IT_0269 + IT_0270 + IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = IT_0044*IT_0091 + IT_0099 + IT_0069*IT_0108 +
       IT_0073*IT_0119 + IT_0075*IT_0125 + conj(IT_0008)*IT_0128 + conj(IT_0078)
      *IT_0133 + IT_0118*IT_0136 + IT_0124*IT_0140 + IT_0122*IT_0141 + IT_0079
      *IT_0142 + IT_0008*IT_0143 + IT_0078*IT_0151 + (IT_0059 + IT_0126)*IT_0152
       + IT_0011*IT_0162 + IT_0068*IT_0177 + IT_0064*IT_0209 + IT_0212 + IT_0213
       + IT_0214 + IT_0215 + IT_0216 + IT_0217 + IT_0262 + IT_0263 + IT_0264 +
       IT_0265 + IT_0266 + IT_0267 + IT_0268 + IT_0269 + IT_0270 + IT_0271 +
       IT_0272;
    const ccomplex_t IT_0275 = IT_0037*IT_0090;
    const ccomplex_t IT_0276 = IT_0034*IT_0150;
    const ccomplex_t IT_0277 = conj(IT_0051)*IT_0133;
    const ccomplex_t IT_0278 = IT_0011*(conj(IT_0044)*IT_0065 + IT_0060*conj
      (IT_0078) + conj(IT_0069)*IT_0081 + s_23*IT_0075*IT_0111 + IT_0079*IT_0148
       + IT_0073*IT_0243 + IT_0246 + IT_0247 + IT_0250 + IT_0251 + IT_0275 +
       IT_0276 + IT_0277);
    const ccomplex_t IT_0279 = conj(IT_0059)*IT_0128;
    const ccomplex_t IT_0280 = conj(IT_0011)*IT_0062;
    const ccomplex_t IT_0281 = conj(IT_0064)*IT_0156;
    const ccomplex_t IT_0282 = conj(IT_0046)*IT_0203;
    const ccomplex_t IT_0283 = conj(IT_0069)*IT_0087;
    const ccomplex_t IT_0284 = conj(IT_0080)*IT_0081;
    const ccomplex_t IT_0285 = conj(IT_0051)*IT_0200;
    const ccomplex_t IT_0286 = 0.5*IT_0229;
    const ccomplex_t IT_0287 = IT_0053*conj(IT_0068) + IT_0062*(conj(IT_0057) 
      + conj(IT_0069)) + IT_0084 + IT_0088 + IT_0092 + IT_0094 + IT_0099 +
       IT_0101 + IT_0104 + IT_0068*IT_0108 + IT_0119*(IT_0075 + IT_0122) +
       IT_0096*IT_0124 + conj(IT_0064)*IT_0128 + conj(IT_0011)*IT_0133 + IT_0122
      *IT_0137 + IT_0134*(IT_0116 + IT_0119 + IT_0141) + (IT_0064 + IT_0078)
      *IT_0143 + IT_0059*(IT_0112 + IT_0144) + IT_0011*IT_0151 + IT_0051*IT_0152
       + conj(IT_0059)*IT_0156 + IT_0126*(IT_0137 + IT_0181) + IT_0185 + IT_0140
      *IT_0186 + IT_0189 + IT_0191 + IT_0034*IT_0192 + IT_0193 + IT_0046*IT_0209
       + IT_0106*(IT_0186 + IT_0014*IT_0210) + IT_0127*IT_0230 + IT_0044*IT_0231
       + IT_0057*IT_0232 + IT_0069*IT_0233 + IT_0073*IT_0235 + IT_0118*IT_0238 +
       conj(IT_0078)*IT_0239 + conj(IT_0044)*IT_0240 + conj(IT_0008)*IT_0241 +
       IT_0080*IT_0243 + IT_0008*(IT_0151 + IT_0243) + IT_0136*IT_0244 + s_23*
      (IT_0075*IT_0210 + (-21.3333333333333)*IT_0014*IT_0230 + IT_0259 + IT_0260
      ) + IT_0282 + IT_0283 + IT_0284 + IT_0285 + IT_0236*IT_0286;
    const ccomplex_t IT_0288 = (conj(IT_0008) + conj(IT_0011))*IT_0040 + conj
      (IT_0044)*(IT_0049 + IT_0050 + IT_0054 + IT_0056) + IT_0057*(conj(IT_0059)
      *IT_0060 + conj(IT_0051)*IT_0062 + conj(IT_0064)*IT_0065 + conj(IT_0046)
      *IT_0067) + conj(IT_0068)*(IT_0049 + IT_0050 + IT_0054 + IT_0056 + IT_0008
      *IT_0065 + IT_0060*(IT_0069 + IT_0073) + (IT_0044 + IT_0075)*IT_0076 +
       IT_0067*(IT_0078 + IT_0079)) + IT_0080*(IT_0060*conj(IT_0064) + conj
      (IT_0059)*IT_0065 + conj(IT_0046)*IT_0076 + conj(IT_0051)*IT_0081) +
       IT_0069*(IT_0084 + IT_0085 + IT_0086 + IT_0088 + IT_0089) + IT_0073*
      (IT_0084 + IT_0085 + IT_0086 + IT_0088 + IT_0089 + IT_0092 + IT_0094 +
       IT_0099 + IT_0101 + IT_0104) + IT_0106*IT_0158 + IT_0075*(conj(IT_0064)
      *IT_0083 + IT_0159 + IT_0160 + IT_0161 + IT_0163 + IT_0164 + IT_0165 +
       IT_0166) + IT_0044*(conj(IT_0064)*IT_0083 + IT_0160 + IT_0164 + IT_0165 +
       IT_0166 + IT_0167) + IT_0068*(IT_0073*IT_0148 + IT_0160 + IT_0164 +
       IT_0167 + IT_0079*IT_0168 + IT_0169 + IT_0170 + IT_0075*IT_0171 + IT_0173
       + IT_0174 + IT_0175 + IT_0176) + IT_0134*IT_0194 + IT_0078*(conj(IT_0011)
      *IT_0060 + conj(IT_0051)*IT_0128 + conj(IT_0046)*IT_0156 + IT_0195 +
       IT_0196) + IT_0079*(conj(IT_0011)*IT_0060 + IT_0195 + IT_0196 + IT_0197 +
       IT_0198 + IT_0199) + IT_0064*IT_0208 + IT_0140*IT_0211 + IT_0059*IT_0218 
      + IT_0046*(conj(IT_0069)*IT_0083 + IT_0079*IT_0112 + conj(IT_0064)*IT_0129
       + conj(IT_0078)*IT_0156 + conj(IT_0008)*IT_0172 + conj(IT_0046)*IT_0200 +
       IT_0219 + IT_0220 + IT_0221 + IT_0222 + IT_0223) + IT_0096*IT_0228 +
       IT_0230*IT_0253 + IT_0034*IT_0258 + IT_0037*IT_0261 + IT_0122*IT_0273 +
       IT_0126*IT_0274 + IT_0278 + IT_0118*(IT_0065*conj(IT_0068) + s_23*IT_0068
      *IT_0111 + IT_0117 + IT_0135 + IT_0011*IT_0232 + IT_0242 + IT_0245 +
       IT_0248 + IT_0249 + IT_0250 + IT_0277 + IT_0279 + IT_0280 + IT_0281) +
       IT_0008*(IT_0250 + IT_0275 + IT_0276 + IT_0277 + IT_0279 + IT_0280 +
       IT_0281) + IT_0051*(conj(IT_0057)*IT_0062 + conj(IT_0078)*IT_0128 + conj
      (IT_0064)*IT_0153 + IT_0133*IT_0190 + IT_0191 + IT_0282 + IT_0283 +
       IT_0284 + IT_0285) + IT_0286*IT_0287;
    return create_ccomplex_return(IT_0288);
}

