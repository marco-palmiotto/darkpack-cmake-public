#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_1_anti_chi_dm_1_to_V_anti_V.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_1_anti_chi_dm_1_to_V_anti_V(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_f = param->g_f;
    const creal_t m_V = param->m_V;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_V_3 = param->Gamma_V_3;
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const ccomplex_t IT_0000 = pow(m_V, 2);
    const ccomplex_t IT_0001 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0002 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = cpow(IT_0008, 2);
    const ccomplex_t IT_0010 = s_12*IT_0000;
    const ccomplex_t IT_0011 = s_14*s_24;
    const ccomplex_t IT_0012 = pow(m_V, -4);
    const ccomplex_t IT_0013 = pow(m_V, 6);
    const ccomplex_t IT_0014 = s_12*IT_0013;
    const ccomplex_t IT_0015 = pow(m_V, 4);
    const ccomplex_t IT_0016 = s_13*s_23*IT_0015;
    const ccomplex_t IT_0017 = (-2)*IT_0016;
    const ccomplex_t IT_0018 = IT_0014 + IT_0017;
    const ccomplex_t IT_0019 = IT_0012*IT_0018;
    const ccomplex_t IT_0020 = pow(m_V, -2);
    const ccomplex_t IT_0021 = s_12*IT_0015;
    const ccomplex_t IT_0022 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0023 = s_13*s_24*s_34;
    const ccomplex_t IT_0024 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0025 = (-2)*IT_0022 + 4*IT_0023 + (-2)*IT_0024;
    const ccomplex_t IT_0026 = IT_0021 + IT_0025;
    const ccomplex_t IT_0027 = IT_0020*IT_0026;
    const ccomplex_t IT_0028 = IT_0000*IT_0001;
    const ccomplex_t IT_0029 = 9*IT_0028;
    const ccomplex_t IT_0030 = (-2)*IT_0008;
    const ccomplex_t IT_0031 = pow(s_24, 2);
    const ccomplex_t IT_0032 = IT_0020*IT_0031;
    const ccomplex_t IT_0033 = s_13*s_23;
    const ccomplex_t IT_0034 = (-2)*IT_0033;
    const ccomplex_t IT_0035 = IT_0010 + IT_0034;
    const ccomplex_t IT_0036 = IT_0020*IT_0035;
    const ccomplex_t IT_0037 = (IT_0001 + -IT_0032)*(s_12 + (-0.5)*IT_0036);
    const ccomplex_t IT_0038 = 2*IT_0037;
    const ccomplex_t IT_0039 = -IT_0032;
    const ccomplex_t IT_0040 = IT_0001 + IT_0039;
    const ccomplex_t IT_0041 = IT_0001*IT_0040;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = s_12*s_24;
    const ccomplex_t IT_0044 = s_14*IT_0001;
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = s_13*s_23*s_24;
    const ccomplex_t IT_0047 = s_14*IT_0000*IT_0001;
    const ccomplex_t IT_0048 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0049 = s_13*s_34*IT_0001;
    const ccomplex_t IT_0050 = 0.5*IT_0047 + (-0.5)*IT_0048 + -IT_0049;
    const ccomplex_t IT_0051 = IT_0046 + IT_0050;
    const ccomplex_t IT_0052 = IT_0020*IT_0051;
    const ccomplex_t IT_0053 = (-2)*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = IT_0043 + IT_0045 + IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + (-2)*IT_0001 + -reg_prop + m_V
      *(m_V + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0058 = 0.5*IT_0057;
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*IT_0059;
    const ccomplex_t IT_0061 = IT_0057*IT_0060;
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = s_13*s_14*s_23*s_24*IT_0012;
    const ccomplex_t IT_0064 = pow(s_13, 2);
    const ccomplex_t IT_0065 = IT_0012*IT_0031*IT_0064;
    const ccomplex_t IT_0066 = pow(m_chi_dm_1, 4);
    const ccomplex_t IT_0067 = s_13*s_24*s_34*IT_0012;
    const ccomplex_t IT_0068 = s_13*s_23*IT_0020;
    const ccomplex_t IT_0069 = s_14*s_24*IT_0020;
    const ccomplex_t IT_0070 = -IT_0068 + -IT_0069;
    const ccomplex_t IT_0071 = s_12 + IT_0067 + IT_0070;
    const ccomplex_t IT_0072 = IT_0001*IT_0071;
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = IT_0001*IT_0020*IT_0031;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0001*IT_0020*IT_0064;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = s_12*s_13*s_24*s_34*IT_0012;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = IT_0063 + IT_0065 + IT_0066 + IT_0073 + IT_0075
       + IT_0077 + IT_0079;
    const ccomplex_t IT_0081 = (-2)*IT_0057;
    const ccomplex_t IT_0082 = IT_0060*IT_0081;
    const ccomplex_t IT_0083 = -IT_0082;
    const ccomplex_t IT_0084 = (s_12 + 1.5*IT_0001 + (-0.5)*IT_0036)*(s_12 + 
      -IT_0069);
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = (IT_0001 + -IT_0032)*(s_12 + 1.5*IT_0001 + (
      -0.5)*IT_0036);
    const ccomplex_t IT_0087 = 2*IT_0086;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = pow(s_34, 2);
    const ccomplex_t IT_0091 = s_12*s_24*IT_0012*IT_0090;
    const ccomplex_t IT_0092 = pow(s_23, 2);
    const ccomplex_t IT_0093 = s_14*IT_0020*IT_0092;
    const ccomplex_t IT_0094 = s_13*s_23*s_24*IT_0020;
    const ccomplex_t IT_0095 = s_14*IT_0020*IT_0031;
    const ccomplex_t IT_0096 = 2*IT_0095;
    const ccomplex_t IT_0097 = -IT_0043;
    const ccomplex_t IT_0098 = IT_0093 + IT_0094 + IT_0096 + IT_0097;
    const ccomplex_t IT_0099 = s_23*s_34*IT_0020;
    const ccomplex_t IT_0100 = s_24*IT_0012*IT_0090;
    const ccomplex_t IT_0101 = -IT_0100;
    const ccomplex_t IT_0102 = IT_0099 + IT_0101;
    const ccomplex_t IT_0103 = IT_0001*IT_0102;
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = s_13*s_34*IT_0012*IT_0031;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = s_12*s_23*s_34*IT_0020;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = s_14*s_23*s_24*s_34*IT_0012;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0045 + IT_0091 + IT_0098 + IT_0104 + IT_0106
       + IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = -IT_0057;
    const ccomplex_t IT_0113 = IT_0060*(IT_0057 + -IT_0112);
    const ccomplex_t IT_0114 = -conj(IT_0113);
    const ccomplex_t IT_0115 = IT_0113 + IT_0114;
    const ccomplex_t IT_0116 = s_13*s_34*IT_0001*IT_0020;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = -IT_0091;
    const ccomplex_t IT_0119 = (-2)*IT_0095;
    const ccomplex_t IT_0120 = IT_0043 + IT_0044 + IT_0103 + IT_0105 + IT_0109
       + IT_0117 + IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = 2*IT_0057;
    const ccomplex_t IT_0122 = IT_0060*(IT_0057 + -IT_0121);
    const ccomplex_t IT_0123 = -conj(IT_0122);
    const ccomplex_t IT_0124 = IT_0122 + IT_0123;
    const ccomplex_t IT_0125 = s_12*IT_0001;
    const ccomplex_t IT_0126 = s_14*s_24*IT_0012*IT_0092;
    const ccomplex_t IT_0127 = s_13*s_23*IT_0012*IT_0031;
    const ccomplex_t IT_0128 = s_23*s_24*s_34*IT_0012;
    const ccomplex_t IT_0129 = IT_0001 + IT_0128;
    const ccomplex_t IT_0130 = IT_0020*IT_0092;
    const ccomplex_t IT_0131 = -IT_0032 + -IT_0130;
    const ccomplex_t IT_0132 = IT_0129 + IT_0131;
    const ccomplex_t IT_0133 = IT_0001*IT_0132;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = s_14*s_24*IT_0001*IT_0020;
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = s_13*s_23*IT_0001*IT_0020;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = s_12*s_23*s_24*s_34*IT_0012;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = IT_0125 + IT_0126 + IT_0127 + IT_0134 + IT_0136
       + IT_0138 + IT_0140;
    const ccomplex_t IT_0142 = 4*IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0062)*IT_0142;
    const ccomplex_t IT_0144 = s_14*s_23*s_34;
    const ccomplex_t IT_0145 = s_12*IT_0090;
    const ccomplex_t IT_0146 = (-2)*IT_0022 + 2*IT_0023 + 2*IT_0144 + (-2)
      *IT_0145;
    const ccomplex_t IT_0147 = IT_0021 + IT_0146;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = IT_0001*IT_0020;
    const ccomplex_t IT_0150 = (-2)*IT_0090;
    const ccomplex_t IT_0151 = IT_0015 + IT_0150;
    const ccomplex_t IT_0152 = IT_0149*IT_0151;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = s_34*IT_0000;
    const ccomplex_t IT_0155 = s_34*IT_0001*IT_0012;
    const ccomplex_t IT_0156 = IT_0154*IT_0155;
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = (-2)*IT_0022;
    const ccomplex_t IT_0159 = IT_0021 + IT_0158;
    const ccomplex_t IT_0160 = IT_0020*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = s_34*IT_0012;
    const ccomplex_t IT_0163 = s_14*s_23*IT_0000;
    const ccomplex_t IT_0164 = s_13*s_24*IT_0000;
    const ccomplex_t IT_0165 = IT_0163 + IT_0164;
    const ccomplex_t IT_0166 = s_12*s_34*IT_0000;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0165 + IT_0167;
    const ccomplex_t IT_0169 = IT_0162*IT_0168;
    const ccomplex_t IT_0170 = 0.5*IT_0169;
    const ccomplex_t IT_0171 = 1.5*IT_0028;
    const ccomplex_t IT_0172 = (-2)*IT_0011;
    const ccomplex_t IT_0173 = IT_0010 + (-0.5)*IT_0148 + IT_0153 + IT_0157 +
       IT_0161 + IT_0170 + IT_0171 + IT_0172;
    const ccomplex_t IT_0174 = 8*IT_0008;
    const ccomplex_t IT_0175 = IT_0020*IT_0064;
    const ccomplex_t IT_0176 = IT_0010 + IT_0172;
    const ccomplex_t IT_0177 = IT_0020*IT_0176;
    const ccomplex_t IT_0178 = (IT_0001 + -IT_0175)*(s_12 + 1.5*IT_0001 + (
      -0.5)*IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = (s_12 + -IT_0068)*(s_12 + 1.5*IT_0001 + (-0.5)
      *IT_0177);
    const ccomplex_t IT_0181 = (IT_0001 + -IT_0130)*(s_12 + 1.5*IT_0001 + (
      -0.5)*IT_0177);
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = s_12*s_14;
    const ccomplex_t IT_0184 = s_13*s_14*s_24*s_34*IT_0012;
    const ccomplex_t IT_0185 = s_24*IT_0001;
    const ccomplex_t IT_0186 = pow(s_14, 2);
    const ccomplex_t IT_0187 = s_24*IT_0020*IT_0186;
    const ccomplex_t IT_0188 = (-2)*IT_0187;
    const ccomplex_t IT_0189 = s_24*IT_0020*IT_0064;
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = s_13*s_14*s_23*IT_0020;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = IT_0183 + 2*IT_0184 + IT_0185 + IT_0188 +
       IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = 4*IT_0193;
    const ccomplex_t IT_0195 = (-2)*IT_0109;
    const ccomplex_t IT_0196 = IT_0045 + IT_0098 + IT_0195;
    const ccomplex_t IT_0197 = (-4)*IT_0196;
    const ccomplex_t IT_0198 = IT_0194 + IT_0197;
    const ccomplex_t IT_0199 = s_13*s_34*IT_0020;
    const ccomplex_t IT_0200 = (s_12 + 1.5*IT_0001 + (-0.5)*IT_0177)*(s_14 + 
      -IT_0199);
    const ccomplex_t IT_0201 = (-8)*IT_0200;
    const ccomplex_t IT_0202 = (s_24 + -IT_0099)*(s_12 + 1.5*IT_0001 + (-0.5)
      *IT_0177);
    const ccomplex_t IT_0203 = (-8)*IT_0202;
    const ccomplex_t IT_0204 = IT_0201 + IT_0203;
    const ccomplex_t IT_0205 = IT_0012*IT_0092*IT_0186;
    const ccomplex_t IT_0206 = s_14*s_23*s_34*IT_0012;
    const ccomplex_t IT_0207 = s_12 + IT_0070 + IT_0206;
    const ccomplex_t IT_0208 = IT_0001*IT_0207;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0001*IT_0020*IT_0092;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0001*IT_0020*IT_0186;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = s_12*s_14*s_23*s_34*IT_0012;
    const ccomplex_t IT_0215 = -IT_0214;
    const ccomplex_t IT_0216 = IT_0063 + IT_0066 + IT_0205 + IT_0209 + IT_0211
       + IT_0213 + IT_0215;
    const ccomplex_t IT_0217 = s_13*s_23*IT_0012*IT_0186;
    const ccomplex_t IT_0218 = s_14*s_24*IT_0012*IT_0064;
    const ccomplex_t IT_0219 = s_13*s_14*s_34*IT_0012;
    const ccomplex_t IT_0220 = IT_0001 + IT_0219;
    const ccomplex_t IT_0221 = IT_0020*IT_0186;
    const ccomplex_t IT_0222 = -IT_0175 + -IT_0221;
    const ccomplex_t IT_0223 = IT_0220 + IT_0222;
    const ccomplex_t IT_0224 = IT_0001*IT_0223;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = s_12*s_13*s_14*s_34*IT_0012;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0125 + IT_0136 + IT_0138 + IT_0217 + IT_0218
       + IT_0225 + IT_0227;
    const ccomplex_t IT_0229 = 4*IT_0228;
    const ccomplex_t IT_0230 = 0.25*IT_0229;
    const ccomplex_t IT_0231 = IT_0080 + IT_0216 + IT_0230;
    const ccomplex_t IT_0232 = -IT_0199;
    const ccomplex_t IT_0233 = s_14 + IT_0232;
    const ccomplex_t IT_0234 = IT_0001*IT_0233;
    const ccomplex_t IT_0235 = s_14*IT_0000;
    const ccomplex_t IT_0236 = s_13*s_34;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = IT_0149*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = s_13*s_14*s_23;
    const ccomplex_t IT_0242 = s_24*IT_0064;
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = s_12*s_13*s_34;
    const ccomplex_t IT_0245 = s_24*IT_0000*IT_0001;
    const ccomplex_t IT_0246 = -IT_0244 + -IT_0245;
    const ccomplex_t IT_0247 = IT_0243 + IT_0246;
    const ccomplex_t IT_0248 = IT_0020*IT_0247;
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = s_12*s_13*s_34*IT_0020;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0185 + IT_0190 + IT_0191 + 2*IT_0234 +
       IT_0240 + IT_0249 + IT_0251;
    const ccomplex_t IT_0253 = 2*IT_0043;
    const ccomplex_t IT_0254 = -IT_0099;
    const ccomplex_t IT_0255 = s_24 + IT_0254;
    const ccomplex_t IT_0256 = IT_0001*IT_0255;
    const ccomplex_t IT_0257 = 2*IT_0256;
    const ccomplex_t IT_0258 = -IT_0094;
    const ccomplex_t IT_0259 = s_14*IT_0092;
    const ccomplex_t IT_0260 = IT_0046 + IT_0259;
    const ccomplex_t IT_0261 = s_12*s_23*s_34;
    const ccomplex_t IT_0262 = -IT_0047 + -IT_0261;
    const ccomplex_t IT_0263 = IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = IT_0020*IT_0263;
    const ccomplex_t IT_0265 = (-0.5)*IT_0264;
    const ccomplex_t IT_0266 = s_24*IT_0000;
    const ccomplex_t IT_0267 = s_23*s_34;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = IT_0266 + IT_0268;
    const ccomplex_t IT_0270 = IT_0149*IT_0269;
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = IT_0045 + IT_0093 + IT_0108 + IT_0252 + IT_0253
       + IT_0257 + IT_0258 + IT_0265 + IT_0271;
    const ccomplex_t IT_0273 = (-0.5)*IT_0030*(IT_0080 + IT_0141) + conj
      (IT_0062)*(IT_0179 + (-2)*IT_0180 + IT_0182) + 0.125*conj(IT_0113)*IT_0198
       + 0.125*conj(IT_0122)*IT_0204 + 0.5*conj(IT_0083)*(0.25*IT_0142 + IT_0231
      ) + 0.125*IT_0174*IT_0272;
    const ccomplex_t IT_0274 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0275 = -IT_0244 + (-0.5)*IT_0245 + 0.5*IT_0274;
    const ccomplex_t IT_0276 = IT_0242 + IT_0275;
    const ccomplex_t IT_0277 = IT_0020*IT_0276;
    const ccomplex_t IT_0278 = 2*IT_0277;
    const ccomplex_t IT_0279 = (-0.5)*IT_0278;
    const ccomplex_t IT_0280 = -IT_0185;
    const ccomplex_t IT_0281 = IT_0183 + IT_0279 + IT_0280;
    const ccomplex_t IT_0282 = s_14*s_23*s_34*IT_0020;
    const ccomplex_t IT_0283 = IT_0020*IT_0090;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = IT_0000 + IT_0284;
    const ccomplex_t IT_0286 = IT_0001*IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0286;
    const ccomplex_t IT_0288 = s_12*IT_0020*IT_0090;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = -IT_0090;
    const ccomplex_t IT_0291 = IT_0015 + IT_0290;
    const ccomplex_t IT_0292 = IT_0149*IT_0291;
    const ccomplex_t IT_0293 = (-0.5)*IT_0292;
    const ccomplex_t IT_0294 = IT_0021 + IT_0023 + IT_0144;
    const ccomplex_t IT_0295 = (-2)*IT_0022 + -IT_0145;
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = IT_0020*IT_0296;
    const ccomplex_t IT_0298 = (-0.5)*IT_0297;
    const ccomplex_t IT_0299 = s_13*s_24*s_34*IT_0020;
    const ccomplex_t IT_0300 = -IT_0299;
    const ccomplex_t IT_0301 = IT_0010 + IT_0282 + IT_0287 + IT_0289 + IT_0293
       + IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = (-8)*IT_0008;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*IT_0006*m_chi_dm_2;
    const ccomplex_t IT_0304 = IT_0003*IT_0303;
    const ccomplex_t IT_0305 = IT_0008*m_chi_dm_1;
    const ccomplex_t IT_0306 = IT_0304 + -IT_0305;
    const ccomplex_t IT_0307 = -IT_0001;
    const ccomplex_t IT_0308 = IT_0175 + IT_0307;
    const ccomplex_t IT_0309 = (-2)*IT_0011 + -IT_0033 + 2*IT_0299;
    const ccomplex_t IT_0310 = IT_0010 + IT_0309;
    const ccomplex_t IT_0311 = IT_0020*IT_0310;
    const ccomplex_t IT_0312 = -IT_0028;
    const ccomplex_t IT_0313 = IT_0064 + IT_0312;
    const ccomplex_t IT_0314 = IT_0020*IT_0313;
    const ccomplex_t IT_0315 = -IT_0068 + (-2)*IT_0308 + (-0.5)*IT_0311 + 0.5
      *IT_0314;
    const ccomplex_t IT_0316 = s_12 + IT_0315;
    const ccomplex_t IT_0317 = IT_0316*m_chi_dm_1;
    const ccomplex_t IT_0318 = (-2)*IT_0317;
    const ccomplex_t IT_0319 = -IT_0033;
    const ccomplex_t IT_0320 = IT_0010 + IT_0319;
    const ccomplex_t IT_0321 = IT_0020*IT_0320;
    const ccomplex_t IT_0322 = s_23*s_24*s_34*IT_0020;
    const ccomplex_t IT_0323 = -IT_0028 + 2*IT_0031 + (-2)*IT_0322;
    const ccomplex_t IT_0324 = IT_0092 + IT_0323;
    const ccomplex_t IT_0325 = IT_0020*IT_0324;
    const ccomplex_t IT_0326 = -IT_0068;
    const ccomplex_t IT_0327 = s_12 + IT_0326;
    const ccomplex_t IT_0328 = -IT_0130 + (-0.5)*IT_0321 + 0.5*IT_0325 + 2
      *IT_0327;
    const ccomplex_t IT_0329 = IT_0001 + IT_0328;
    const ccomplex_t IT_0330 = IT_0329*m_chi_dm_1;
    const ccomplex_t IT_0331 = (-2)*IT_0330;
    const ccomplex_t IT_0332 = s_13*s_14*s_34*IT_0020;
    const ccomplex_t IT_0333 = (-0.5)*IT_0028 + 0.5*IT_0186 + -IT_0332;
    const ccomplex_t IT_0334 = IT_0064 + IT_0333;
    const ccomplex_t IT_0335 = IT_0020*IT_0334;
    const ccomplex_t IT_0336 = IT_0001 + IT_0335;
    const ccomplex_t IT_0337 = -IT_0011;
    const ccomplex_t IT_0338 = IT_0010 + IT_0337;
    const ccomplex_t IT_0339 = IT_0020*IT_0338;
    const ccomplex_t IT_0340 = -IT_0069;
    const ccomplex_t IT_0341 = s_12 + IT_0340;
    const ccomplex_t IT_0342 = -IT_0221 + (-0.5)*IT_0339 + 2*IT_0341;
    const ccomplex_t IT_0343 = IT_0336 + IT_0342;
    const ccomplex_t IT_0344 = IT_0343*m_chi_dm_1;
    const ccomplex_t IT_0345 = (-2)*IT_0344;
    const ccomplex_t IT_0346 = -IT_0011 + (-2)*IT_0033 + 2*IT_0299;
    const ccomplex_t IT_0347 = IT_0010 + IT_0346;
    const ccomplex_t IT_0348 = IT_0020*IT_0347;
    const ccomplex_t IT_0349 = IT_0031 + IT_0312;
    const ccomplex_t IT_0350 = IT_0020*IT_0349;
    const ccomplex_t IT_0351 = IT_0032 + IT_0307;
    const ccomplex_t IT_0352 = -IT_0069 + (-0.5)*IT_0348 + 0.5*IT_0350 + (-2)
      *IT_0351;
    const ccomplex_t IT_0353 = s_12 + IT_0352;
    const ccomplex_t IT_0354 = IT_0353*m_chi_dm_1;
    const ccomplex_t IT_0355 = (-2)*IT_0354;
    const ccomplex_t IT_0356 = -IT_0069 + (-0.5)*IT_0348;
    const ccomplex_t IT_0357 = s_12 + IT_0356;
    const ccomplex_t IT_0358 = IT_0357*m_chi_dm_1;
    const ccomplex_t IT_0359 = (-2)*IT_0358;
    const ccomplex_t IT_0360 = (-4)*IT_0351;
    const ccomplex_t IT_0361 = IT_0350 + IT_0360;
    const ccomplex_t IT_0362 = IT_0361*m_chi_dm_1;
    const ccomplex_t IT_0363 = 1.125*IT_0001;
    const ccomplex_t IT_0364 = IT_0012*IT_0026;
    const ccomplex_t IT_0365 = s_13*IT_0000;
    const ccomplex_t IT_0366 = IT_0162*IT_0365;
    const ccomplex_t IT_0367 = IT_0020*IT_0235;
    const ccomplex_t IT_0368 = IT_0020*IT_0266;
    const ccomplex_t IT_0369 = s_23*IT_0000;
    const ccomplex_t IT_0370 = s_24*s_34;
    const ccomplex_t IT_0371 = (-2)*IT_0370;
    const ccomplex_t IT_0372 = IT_0369 + IT_0371;
    const ccomplex_t IT_0373 = IT_0162*IT_0372;
    const ccomplex_t IT_0374 = -s_24 + 0.25*IT_0366 + (-0.5)*IT_0367 + 0.5
      *IT_0368 + 0.25*IT_0373;
    const ccomplex_t IT_0375 = s_14 + IT_0374;
    const ccomplex_t IT_0376 = IT_0375*m_chi_dm_1;
    const ccomplex_t IT_0377 = 4*IT_0376;
    const ccomplex_t IT_0378 = IT_0020*IT_0238;
    const ccomplex_t IT_0379 = IT_0266 + IT_0267;
    const ccomplex_t IT_0380 = s_24*IT_0020*IT_0090;
    const ccomplex_t IT_0381 = (-2)*IT_0380;
    const ccomplex_t IT_0382 = IT_0379 + IT_0381;
    const ccomplex_t IT_0383 = IT_0020*IT_0382;
    const ccomplex_t IT_0384 = -IT_0099 + 2*IT_0233 + (-0.5)*IT_0378 + 0.5
      *IT_0383;
    const ccomplex_t IT_0385 = s_24 + IT_0384;
    const ccomplex_t IT_0386 = IT_0385*m_chi_dm_1;
    const ccomplex_t IT_0387 = (-2)*IT_0386;
    const ccomplex_t IT_0388 = s_14*IT_0015;
    const ccomplex_t IT_0389 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0390 = (-2)*IT_0389;
    const ccomplex_t IT_0391 = IT_0388 + IT_0390;
    const ccomplex_t IT_0392 = IT_0012*IT_0391;
    const ccomplex_t IT_0393 = (-2)*IT_0236;
    const ccomplex_t IT_0394 = IT_0235 + IT_0393;
    const ccomplex_t IT_0395 = IT_0020*IT_0394;
    const ccomplex_t IT_0396 = 0.166666666666667*IT_0392 + (-0.666666666666667
      )*IT_0395;
    const ccomplex_t IT_0397 = s_14 + IT_0396;
    const ccomplex_t IT_0398 = IT_0397*m_chi_dm_1;
    const ccomplex_t IT_0399 = 6*IT_0398;
    const ccomplex_t IT_0400 = s_24*IT_0015;
    const ccomplex_t IT_0401 = IT_0012*IT_0400;
    const ccomplex_t IT_0402 = (-10)*s_24;
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = IT_0403*m_chi_dm_1;
    const ccomplex_t IT_0405 = (IT_0062 + -conj(IT_0062))*(IT_0318 + IT_0331) 
      + -(IT_0083 + -conj(IT_0083))*(IT_0345 + IT_0355) + (-2)*IT_0030*(IT_0359 
      + -IT_0362) + (-0.25)*IT_0306*((-16)*s_12 + 8*IT_0036 + 8*IT_0177 + 32
      *IT_0363 + (-4)*IT_0364) + -IT_0115*IT_0377 + IT_0124*IT_0387 + (-0.25)
      *IT_0302*(IT_0399 + IT_0404);
    const ccomplex_t IT_0406 = (s_12 + 1.5*IT_0001 + (-0.5)*IT_0036)*(IT_0001 
      + -IT_0221);
    const ccomplex_t IT_0407 = -IT_0406;
    const ccomplex_t IT_0408 = -IT_0086;
    const ccomplex_t IT_0409 = s_23*s_34*IT_0012*IT_0186;
    const ccomplex_t IT_0410 = s_14*IT_0012*IT_0090;
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = IT_0199 + IT_0411;
    const ccomplex_t IT_0413 = IT_0001*IT_0412;
    const ccomplex_t IT_0414 = s_12*s_14*IT_0012*IT_0090;
    const ccomplex_t IT_0415 = -IT_0414;
    const ccomplex_t IT_0416 = s_23*s_34*IT_0001*IT_0020;
    const ccomplex_t IT_0417 = -IT_0416;
    const ccomplex_t IT_0418 = IT_0183 + IT_0184 + IT_0185 + IT_0188 + IT_0409
       + IT_0413 + IT_0415 + IT_0417;
    const ccomplex_t IT_0419 = IT_0120 + IT_0418;
    const ccomplex_t IT_0420 = -IT_0111;
    const ccomplex_t IT_0421 = (-0.5)*IT_0030*IT_0089 + 0.125*IT_0143 + 0.5
      *conj(IT_0062)*IT_0231 + 0.125*IT_0174*IT_0281 + 0.125*IT_0055*IT_0302 + (
      -2)*conj(IT_0083)*(IT_0084 + -1./2*IT_0407 + -1./2*IT_0408) + 0.5*conj
      (IT_0122)*IT_0419 + 0.5*conj(IT_0113)*(IT_0183 + IT_0184 + IT_0185 +
       IT_0188 + IT_0190 + IT_0192 + IT_0250 + IT_0409 + IT_0413 + IT_0415 +
       IT_0420);
    const ccomplex_t IT_0422 = IT_0012*IT_0090;
    const ccomplex_t IT_0423 = (IT_0010 + (-2)*IT_0011 + IT_0028)*(1 + 0.5
      *IT_0422);
    const ccomplex_t IT_0424 = s_14*s_24*IT_0012*IT_0090;
    const ccomplex_t IT_0425 = (-2)*IT_0424;
    const ccomplex_t IT_0426 = IT_0282 + IT_0299 + IT_0425;
    const ccomplex_t IT_0427 = (-4)*IT_0426;
    const ccomplex_t IT_0428 = (-0.5)*conj(IT_0083);
    const ccomplex_t IT_0429 = (-0.125)*conj(IT_0062);
    const ccomplex_t IT_0430 = (-0.125)*IT_0173*IT_0302 + conj(IT_0113)
      *IT_0423 + (-0.125)*conj(IT_0122)*IT_0427 + (IT_0183 + IT_0184 + IT_0185 +
       IT_0188 + IT_0190 + IT_0192 + IT_0250 + IT_0409 + IT_0413 + IT_0415 +
       IT_0420)*IT_0428 + IT_0198*IT_0429;
    const ccomplex_t IT_0431 = (s_12 + 1.5*IT_0001 + (-0.5)*IT_0177)*(IT_0000 
      + -IT_0283);
    const ccomplex_t IT_0432 = 8*IT_0009*(IT_0010 + (-4)*IT_0011 + 0.5*IT_0019
       + -IT_0027 + 0.5*IT_0029) + 4*IT_0030*(IT_0030*(IT_0038 + (-3)*IT_0042) +
       (-4)*IT_0008*IT_0055 + conj(IT_0062)*IT_0080 + conj(IT_0083)*IT_0089 +
       IT_0111*IT_0115 + -IT_0120*IT_0124 + 0.25*IT_0143) + (IT_0055*conj
      (IT_0083) + conj(IT_0113)*IT_0173)*IT_0174 + 8*IT_0062*IT_0273 + (conj
      (IT_0062)*IT_0272 + conj(IT_0083)*IT_0281 + conj(IT_0122)*IT_0301)*IT_0302
       + (-4)*IT_0306*IT_0405 + 8*IT_0083*IT_0421 + (-8)*IT_0113*IT_0430 + (-8)
      *IT_0122*((-0.125)*IT_0174*IT_0301 + (-0.125)*conj(IT_0113)*IT_0427 +
       IT_0419*IT_0428 + IT_0204*IT_0429 + conj(IT_0122)*IT_0431);
    return create_ccomplex_return(IT_0432);
}

