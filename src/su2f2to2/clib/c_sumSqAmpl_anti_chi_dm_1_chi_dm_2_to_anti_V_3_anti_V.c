#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_3_anti_V.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_anti_chi_dm_1_chi_dm_2_to_anti_V_3_anti_V(
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
    const creal_t Gamma_V = param->Gamma_V;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_chi_dm_1 = param->m_chi_dm_1;
    const creal_t m_chi_dm_2 = param->m_chi_dm_2;
    const ccomplex_t IT_0000 = pow(m_V, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_13 + IT_0000 + reg_prop, -1);
    const ccomplex_t IT_0002 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0003 = 0.5*IT_0002;
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = IT_0003*IT_0005;
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = -IT_0004;
    const ccomplex_t IT_0010 = 2*IT_0009;
    const ccomplex_t IT_0011 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0012 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_V*(m_V + (0 + _Complex_I*-1)
      *Gamma_V) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0003*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0010*IT_0015;
    const ccomplex_t IT_0017 = (-2)*IT_0008 + -IT_0016;
    const ccomplex_t IT_0018 = 2*conj(IT_0017);
    const ccomplex_t IT_0019 = IT_0009*IT_0015;
    const ccomplex_t IT_0020 = pow(m_V, -2);
    const ccomplex_t IT_0021 = s_34*IT_0020*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0022 = s_13*s_34*IT_0020;
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = s_14 + IT_0023;
    const ccomplex_t IT_0025 = IT_0021*IT_0024;
    const ccomplex_t IT_0026 = 2*IT_0025;
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = s_23*IT_0000;
    const ccomplex_t IT_0029 = IT_0020*IT_0028;
    const ccomplex_t IT_0030 = s_24*IT_0000;
    const ccomplex_t IT_0031 = pow(m_V, -4);
    const ccomplex_t IT_0032 = s_34*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = s_24*s_34;
    const ccomplex_t IT_0035 = (-2)*IT_0034;
    const ccomplex_t IT_0036 = IT_0028 + IT_0035;
    const ccomplex_t IT_0037 = IT_0020*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0029 + (-0.5)*IT_0033 + (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = s_23 + IT_0038;
    const ccomplex_t IT_0040 = IT_0039*m_chi_dm_1;
    const ccomplex_t IT_0041 = 2*IT_0040;
    const ccomplex_t IT_0042 = IT_0004*m_chi_dm_1;
    const ccomplex_t IT_0043 = 2*IT_0004;
    const ccomplex_t IT_0044 = IT_0043*m_chi_dm_2;
    const ccomplex_t IT_0045 = IT_0042 + IT_0044;
    const ccomplex_t IT_0046 = IT_0004*m_chi_dm_2;
    const ccomplex_t IT_0047 = IT_0043*m_chi_dm_1;
    const ccomplex_t IT_0048 = -IT_0046 + -IT_0047;
    const ccomplex_t IT_0049 = IT_0045 + IT_0048;
    const ccomplex_t IT_0050 = IT_0015*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0053 = s_34*IT_0000;
    const ccomplex_t IT_0054 = IT_0032*IT_0053;
    const ccomplex_t IT_0055 = pow(m_V, 4);
    const ccomplex_t IT_0056 = pow(s_34, 2);
    const ccomplex_t IT_0057 = (-2)*IT_0056;
    const ccomplex_t IT_0058 = IT_0055 + IT_0057;
    const ccomplex_t IT_0059 = IT_0020*IT_0058;
    const ccomplex_t IT_0060 = (-0.333333333333333)*IT_0054 + (
      -0.333333333333333)*IT_0059;
    const ccomplex_t IT_0061 = IT_0000 + IT_0060;
    const ccomplex_t IT_0062 = IT_0052*IT_0061;
    const ccomplex_t IT_0063 = 3*IT_0062;
    const ccomplex_t IT_0064 = IT_0015*IT_0043;
    const ccomplex_t IT_0065 = IT_0020*IT_0056;
    const ccomplex_t IT_0066 = -IT_0056;
    const ccomplex_t IT_0067 = IT_0055 + IT_0066;
    const ccomplex_t IT_0068 = IT_0020*IT_0067;
    const ccomplex_t IT_0069 = -IT_0065 + (-0.25)*IT_0068;
    const ccomplex_t IT_0070 = IT_0000 + IT_0069;
    const ccomplex_t IT_0071 = IT_0052*IT_0070;
    const ccomplex_t IT_0072 = (-4)*IT_0071;
    const ccomplex_t IT_0073 = (-2)*IT_0004;
    const ccomplex_t IT_0074 = IT_0015*(IT_0009 + -IT_0073);
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = s_14*s_34*IT_0020;
    const ccomplex_t IT_0077 = s_13*IT_0000;
    const ccomplex_t IT_0078 = s_14*s_34;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = IT_0077 + IT_0079;
    const ccomplex_t IT_0081 = IT_0020*IT_0080;
    const ccomplex_t IT_0082 = -IT_0076 + (-0.25)*IT_0081;
    const ccomplex_t IT_0083 = s_13 + IT_0082;
    const ccomplex_t IT_0084 = IT_0052*IT_0083;
    const ccomplex_t IT_0085 = (-4)*IT_0084;
    const ccomplex_t IT_0086 = s_24*s_34*IT_0020;
    const ccomplex_t IT_0087 = -IT_0034;
    const ccomplex_t IT_0088 = IT_0028 + IT_0087;
    const ccomplex_t IT_0089 = IT_0020*IT_0088;
    const ccomplex_t IT_0090 = -IT_0086 + (-0.25)*IT_0089;
    const ccomplex_t IT_0091 = s_23 + IT_0090;
    const ccomplex_t IT_0092 = IT_0052*IT_0091;
    const ccomplex_t IT_0093 = (-4)*IT_0092;
    const ccomplex_t IT_0094 = IT_0085 + IT_0093;
    const ccomplex_t IT_0095 = s_12*s_13;
    const ccomplex_t IT_0096 = s_23*IT_0011;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = IT_0095 + IT_0097;
    const ccomplex_t IT_0099 = s_13*s_14*s_24;
    const ccomplex_t IT_0100 = s_23*IT_0000*IT_0011;
    const ccomplex_t IT_0101 = s_12*s_13*IT_0000;
    const ccomplex_t IT_0102 = s_24*s_34*IT_0011;
    const ccomplex_t IT_0103 = 0.5*IT_0100 + (-0.5)*IT_0101 + -IT_0102;
    const ccomplex_t IT_0104 = IT_0099 + IT_0103;
    const ccomplex_t IT_0105 = IT_0020*IT_0104;
    const ccomplex_t IT_0106 = (-2)*IT_0105;
    const ccomplex_t IT_0107 = (-0.5)*IT_0106;
    const ccomplex_t IT_0108 = IT_0098 + IT_0107;
    const ccomplex_t IT_0109 = 4*IT_0008;
    const ccomplex_t IT_0110 = s_12*IT_0000;
    const ccomplex_t IT_0111 = s_12*IT_0055;
    const ccomplex_t IT_0112 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0113 = (-2)*IT_0112;
    const ccomplex_t IT_0114 = IT_0111 + IT_0113;
    const ccomplex_t IT_0115 = IT_0020*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = s_14*s_23*s_34;
    const ccomplex_t IT_0118 = s_13*s_24*s_34;
    const ccomplex_t IT_0119 = s_12*IT_0056;
    const ccomplex_t IT_0120 = (-2)*IT_0112 + 2*IT_0117 + 2*IT_0118 + (-2)
      *IT_0119;
    const ccomplex_t IT_0121 = IT_0111 + IT_0120;
    const ccomplex_t IT_0122 = IT_0020*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = s_14*s_23*IT_0000;
    const ccomplex_t IT_0125 = s_13*s_24*IT_0000;
    const ccomplex_t IT_0126 = IT_0124 + IT_0125;
    const ccomplex_t IT_0127 = s_12*s_34*IT_0000;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = IT_0032*IT_0129;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = s_13*s_23;
    const ccomplex_t IT_0133 = (-2)*IT_0132;
    const ccomplex_t IT_0134 = IT_0110 + IT_0116 + IT_0123 + IT_0131 + IT_0133;
    const ccomplex_t IT_0135 = 8*IT_0008;
    const ccomplex_t IT_0136 = s_14*s_23*s_24*IT_0020;
    const ccomplex_t IT_0137 = s_13*IT_0012;
    const ccomplex_t IT_0138 = pow(s_24, 2);
    const ccomplex_t IT_0139 = s_13*IT_0020*IT_0138;
    const ccomplex_t IT_0140 = -IT_0139;
    const ccomplex_t IT_0141 = s_14*s_23*s_24;
    const ccomplex_t IT_0142 = s_13*IT_0138;
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = s_12*s_24*s_34;
    const ccomplex_t IT_0145 = s_13*IT_0000*IT_0012;
    const ccomplex_t IT_0146 = -IT_0144 + -IT_0145;
    const ccomplex_t IT_0147 = IT_0143 + IT_0146;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = s_12*s_24*s_34*IT_0020;
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = IT_0136 + IT_0137 + IT_0140 + IT_0149 + IT_0151;
    const ccomplex_t IT_0153 = s_14*s_23*s_34*IT_0020;
    const ccomplex_t IT_0154 = s_12*IT_0020*IT_0056;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = s_13*s_24*s_34*IT_0020;
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = IT_0111 + IT_0117 + IT_0118;
    const ccomplex_t IT_0159 = (-2)*IT_0112 + -IT_0119;
    const ccomplex_t IT_0160 = IT_0158 + IT_0159;
    const ccomplex_t IT_0161 = IT_0020*IT_0160;
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = (-8)*IT_0008;
    const ccomplex_t IT_0164 = IT_0020*IT_0077;
    const ccomplex_t IT_0165 = s_14*IT_0000;
    const ccomplex_t IT_0166 = IT_0032*IT_0165;
    const ccomplex_t IT_0167 = (-2)*IT_0078;
    const ccomplex_t IT_0168 = IT_0077 + IT_0167;
    const ccomplex_t IT_0169 = IT_0020*IT_0168;
    const ccomplex_t IT_0170 = 0.5*IT_0164 + (-0.5)*IT_0166 + (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = s_13 + IT_0170;
    const ccomplex_t IT_0172 = IT_0171*m_chi_dm_2;
    const ccomplex_t IT_0173 = 2*IT_0172;
    const ccomplex_t IT_0174 = (-4)*IT_0008;
    const ccomplex_t IT_0175 = cpow((-2)*s_23 + IT_0000 + reg_prop, -1);
    const ccomplex_t IT_0176 = (-0.5)*IT_0004;
    const ccomplex_t IT_0177 = IT_0003*IT_0176;
    const ccomplex_t IT_0178 = IT_0175*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = -IT_0016 + 2*IT_0179;
    const ccomplex_t IT_0181 = 2*conj(IT_0180);
    const ccomplex_t IT_0182 = s_23*s_34*IT_0020;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = s_24 + IT_0183;
    const ccomplex_t IT_0185 = IT_0021*IT_0184;
    const ccomplex_t IT_0186 = 2*IT_0185;
    const ccomplex_t IT_0187 = IT_0019*IT_0186;
    const ccomplex_t IT_0188 = s_12*s_23;
    const ccomplex_t IT_0189 = -IT_0137;
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = s_12*s_23*IT_0000;
    const ccomplex_t IT_0192 = -IT_0144 + (-0.5)*IT_0145 + 0.5*IT_0191;
    const ccomplex_t IT_0193 = IT_0142 + IT_0192;
    const ccomplex_t IT_0194 = IT_0020*IT_0193;
    const ccomplex_t IT_0195 = 2*IT_0194;
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0190 + IT_0196;
    const ccomplex_t IT_0198 = IT_0174*IT_0197;
    const ccomplex_t IT_0199 = IT_0110 + IT_0133;
    const ccomplex_t IT_0200 = IT_0020*IT_0199;
    const ccomplex_t IT_0201 = (IT_0000 + -IT_0065)*(s_12 + (-0.5)*IT_0200);
    const ccomplex_t IT_0202 = (-2)*IT_0201;
    const ccomplex_t IT_0203 = -IT_0065;
    const ccomplex_t IT_0204 = IT_0000 + IT_0203;
    const ccomplex_t IT_0205 = IT_0052*IT_0204;
    const ccomplex_t IT_0206 = (s_13 + -IT_0076)*(s_12 + (-0.5)*IT_0200);
    const ccomplex_t IT_0207 = (-2)*IT_0206;
    const ccomplex_t IT_0208 = 4*IT_0207;
    const ccomplex_t IT_0209 = (s_23 + -IT_0086)*(s_12 + (-0.5)*IT_0200);
    const ccomplex_t IT_0210 = (-2)*IT_0209;
    const ccomplex_t IT_0211 = 4*IT_0210;
    const ccomplex_t IT_0212 = -IT_0086;
    const ccomplex_t IT_0213 = s_23 + IT_0212;
    const ccomplex_t IT_0214 = IT_0052*IT_0213;
    const ccomplex_t IT_0215 = (-12)*IT_0214;
    const ccomplex_t IT_0216 = -IT_0076;
    const ccomplex_t IT_0217 = s_13 + IT_0216;
    const ccomplex_t IT_0218 = IT_0052*IT_0217;
    const ccomplex_t IT_0219 = (-12)*IT_0218;
    const ccomplex_t IT_0220 = IT_0208 + IT_0211 + IT_0215 + IT_0219;
    const ccomplex_t IT_0221 = s_34*IT_0020*m_chi_dm_1;
    const ccomplex_t IT_0222 = IT_0184*IT_0221;
    const ccomplex_t IT_0223 = 4*IT_0222;
    const ccomplex_t IT_0224 = s_34*IT_0020*m_chi_dm_2;
    const ccomplex_t IT_0225 = IT_0024*IT_0224;
    const ccomplex_t IT_0226 = (-4)*IT_0225;
    const ccomplex_t IT_0227 = IT_0223 + IT_0226;
    const ccomplex_t IT_0228 = s_13*s_14*s_23*s_34*IT_0031;
    const ccomplex_t IT_0229 = pow(s_13, 2);
    const ccomplex_t IT_0230 = s_24*s_34*IT_0031*IT_0229;
    const ccomplex_t IT_0231 = s_23*IT_0020*IT_0229;
    const ccomplex_t IT_0232 = (-2)*IT_0231;
    const ccomplex_t IT_0233 = IT_0095 + IT_0232;
    const ccomplex_t IT_0234 = s_12*s_13*IT_0031*IT_0056;
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = s_24*s_34*IT_0011*IT_0020;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = IT_0096 + IT_0228 + IT_0230 + IT_0233 + IT_0235
       + IT_0237;
    const ccomplex_t IT_0239 = 2*IT_0238;
    const ccomplex_t IT_0240 = s_13*s_23*s_24*s_34*IT_0031;
    const ccomplex_t IT_0241 = pow(s_23, 2);
    const ccomplex_t IT_0242 = s_14*s_34*IT_0031*IT_0241;
    const ccomplex_t IT_0243 = s_12*s_23*IT_0031*IT_0056;
    const ccomplex_t IT_0244 = s_14*s_34*IT_0012*IT_0020;
    const ccomplex_t IT_0245 = s_13*IT_0020*IT_0241;
    const ccomplex_t IT_0246 = (-2)*IT_0245;
    const ccomplex_t IT_0247 = 2*IT_0137 + 2*IT_0188 + 2*IT_0240 + 2*IT_0242 +
       (-2)*IT_0243 + (-2)*IT_0244 + 2*IT_0246;
    const ccomplex_t IT_0248 = s_34*IT_0020;
    const ccomplex_t IT_0249 = s_14*s_23;
    const ccomplex_t IT_0250 = s_13*s_24;
    const ccomplex_t IT_0251 = IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = s_13*s_23*s_34*IT_0020;
    const ccomplex_t IT_0253 = (-2)*IT_0252;
    const ccomplex_t IT_0254 = IT_0251 + IT_0253;
    const ccomplex_t IT_0255 = IT_0248*IT_0254;
    const ccomplex_t IT_0256 = (-4)*IT_0255;
    const ccomplex_t IT_0257 = 0.25*IT_0135*(IT_0110 + IT_0153 + IT_0155 +
       IT_0157 + IT_0162) + 0.25*IT_0072*IT_0174 + conj(IT_0019)*(IT_0202 + (-3)
      *IT_0205) + 0.25*conj(IT_0075)*IT_0220 + 0.25*conj(IT_0051)*IT_0227 + 0.25
      *IT_0018*IT_0239 + 0.25*IT_0181*IT_0247 + 0.25*conj(IT_0064)*IT_0256;
    const ccomplex_t IT_0258 = IT_0031*IT_0056;
    const ccomplex_t IT_0259 = 0.5*IT_0258;
    const ccomplex_t IT_0260 = 1 + IT_0259;
    const ccomplex_t IT_0261 = s_12*IT_0260;
    const ccomplex_t IT_0262 = IT_0052*IT_0260;
    const ccomplex_t IT_0263 = IT_0020*IT_0241;
    const ccomplex_t IT_0264 = IT_0020*IT_0138;
    const ccomplex_t IT_0265 = IT_0263 + IT_0264;
    const ccomplex_t IT_0266 = s_23*s_24*s_34*IT_0031;
    const ccomplex_t IT_0267 = -IT_0012 + -IT_0266;
    const ccomplex_t IT_0268 = IT_0265 + IT_0267;
    const ccomplex_t IT_0269 = IT_0268*m_chi_dm_1;
    const ccomplex_t IT_0270 = s_13*s_24*s_34*IT_0031;
    const ccomplex_t IT_0271 = s_13*s_23*IT_0020;
    const ccomplex_t IT_0272 = s_14*s_24*IT_0020;
    const ccomplex_t IT_0273 = -IT_0271 + -IT_0272;
    const ccomplex_t IT_0274 = s_12 + IT_0270 + IT_0273;
    const ccomplex_t IT_0275 = IT_0274*m_chi_dm_2;
    const ccomplex_t IT_0276 = IT_0020*IT_0229;
    const ccomplex_t IT_0277 = pow(s_14, 2);
    const ccomplex_t IT_0278 = IT_0020*IT_0277;
    const ccomplex_t IT_0279 = IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = s_13*s_14*s_34*IT_0031;
    const ccomplex_t IT_0281 = -IT_0011 + -IT_0280;
    const ccomplex_t IT_0282 = IT_0279 + IT_0281;
    const ccomplex_t IT_0283 = IT_0282*m_chi_dm_2;
    const ccomplex_t IT_0284 = (-4)*IT_0283;
    const ccomplex_t IT_0285 = 0.25*IT_0284;
    const ccomplex_t IT_0286 = s_14*s_23*s_34*IT_0031;
    const ccomplex_t IT_0287 = s_12 + IT_0273 + IT_0286;
    const ccomplex_t IT_0288 = IT_0287*m_chi_dm_1;
    const ccomplex_t IT_0289 = (-4)*IT_0288;
    const ccomplex_t IT_0290 = 0.25*IT_0289;
    const ccomplex_t IT_0291 = IT_0269 + IT_0275 + IT_0285 + IT_0290;
    const ccomplex_t IT_0292 = s_13*m_chi_dm_2;
    const ccomplex_t IT_0293 = IT_0260*IT_0292;
    const ccomplex_t IT_0294 = 8*IT_0293;
    const ccomplex_t IT_0295 = s_23*m_chi_dm_1;
    const ccomplex_t IT_0296 = IT_0260*IT_0295;
    const ccomplex_t IT_0297 = (-8)*IT_0296;
    const ccomplex_t IT_0298 = IT_0294 + IT_0297;
    const ccomplex_t IT_0299 = IT_0287*m_chi_dm_2;
    const ccomplex_t IT_0300 = 2*IT_0299;
    const ccomplex_t IT_0301 = 2*IT_0269;
    const ccomplex_t IT_0302 = IT_0274*m_chi_dm_1;
    const ccomplex_t IT_0303 = (-2)*IT_0302;
    const ccomplex_t IT_0304 = (-2)*IT_0283;
    const ccomplex_t IT_0305 = 0.125*IT_0109*IT_0173 + 0.125*IT_0041*IT_0174 +
       0.125*conj(IT_0019)*IT_0227 + conj(IT_0051)*(IT_0261 + -IT_0262) + 0.5
      *conj(IT_0075)*IT_0291 + 0.125*conj(IT_0064)*IT_0298 + 0.125*IT_0181*
      (IT_0300 + IT_0301) + 0.125*IT_0018*(IT_0303 + IT_0304);
    const ccomplex_t IT_0306 = s_23*IT_0020*IT_0277;
    const ccomplex_t IT_0307 = s_13*s_14*s_24*IT_0020;
    const ccomplex_t IT_0308 = s_23*IT_0277;
    const ccomplex_t IT_0309 = IT_0099 + IT_0308;
    const ccomplex_t IT_0310 = s_12*s_14*s_34;
    const ccomplex_t IT_0311 = -IT_0100 + -IT_0310;
    const ccomplex_t IT_0312 = IT_0309 + IT_0311;
    const ccomplex_t IT_0313 = IT_0020*IT_0312;
    const ccomplex_t IT_0314 = (-0.25)*IT_0313;
    const ccomplex_t IT_0315 = s_12*s_14*s_34*IT_0020;
    const ccomplex_t IT_0316 = (-0.5)*IT_0315;
    const ccomplex_t IT_0317 = (-0.5)*IT_0096;
    const ccomplex_t IT_0318 = IT_0095 + 0.5*IT_0306 + (-0.5)*IT_0307 +
       IT_0314 + IT_0316 + IT_0317;
    const ccomplex_t IT_0319 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0320 = (-2)*IT_0112 + 4*IT_0118 + (-2)*IT_0319;
    const ccomplex_t IT_0321 = IT_0111 + IT_0320;
    const ccomplex_t IT_0322 = IT_0020*IT_0321;
    const ccomplex_t IT_0323 = (-4)*IT_0132;
    const ccomplex_t IT_0324 = pow(m_V, 6);
    const ccomplex_t IT_0325 = s_12*IT_0324;
    const ccomplex_t IT_0326 = s_14*s_24*IT_0055;
    const ccomplex_t IT_0327 = (-2)*IT_0326;
    const ccomplex_t IT_0328 = IT_0325 + IT_0327;
    const ccomplex_t IT_0329 = IT_0031*IT_0328;
    const ccomplex_t IT_0330 = 0.5*IT_0329;
    const ccomplex_t IT_0331 = IT_0000*IT_0052;
    const ccomplex_t IT_0332 = 9*IT_0331;
    const ccomplex_t IT_0333 = 0.5*IT_0332;
    const ccomplex_t IT_0334 = -IT_0179;
    const ccomplex_t IT_0335 = (-0.25)*IT_0148;
    const ccomplex_t IT_0336 = (-0.5)*IT_0150;
    const ccomplex_t IT_0337 = (-0.5)*IT_0137;
    const ccomplex_t IT_0338 = 0.5*IT_0139;
    const ccomplex_t IT_0339 = (-0.5)*IT_0136;
    const ccomplex_t IT_0340 = IT_0188 + IT_0335 + IT_0336 + IT_0337 + IT_0338
       + IT_0339;
    const ccomplex_t IT_0341 = (-4)*IT_0340;
    const ccomplex_t IT_0342 = -IT_0307;
    const ccomplex_t IT_0343 = 0.5*IT_0313;
    const ccomplex_t IT_0344 = IT_0097 + IT_0306 + IT_0315 + IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = 2*IT_0344;
    const ccomplex_t IT_0346 = IT_0094 + IT_0341 + IT_0345;
    const ccomplex_t IT_0347 = -IT_0153;
    const ccomplex_t IT_0348 = (-0.5)*IT_0100 + 0.5*IT_0101 + -IT_0310;
    const ccomplex_t IT_0349 = IT_0308 + IT_0348;
    const ccomplex_t IT_0350 = IT_0020*IT_0349;
    const ccomplex_t IT_0351 = 2*IT_0350;
    const ccomplex_t IT_0352 = (-0.5)*IT_0351;
    const ccomplex_t IT_0353 = IT_0098 + IT_0352;
    const ccomplex_t IT_0354 = 4*IT_0353;
    const ccomplex_t IT_0355 = s_14*s_34*IT_0012;
    const ccomplex_t IT_0356 = 0.5*IT_0145 + (-0.5)*IT_0191 + -IT_0355;
    const ccomplex_t IT_0357 = IT_0141 + IT_0356;
    const ccomplex_t IT_0358 = IT_0020*IT_0357;
    const ccomplex_t IT_0359 = (-2)*IT_0358;
    const ccomplex_t IT_0360 = (-0.5)*IT_0359;
    const ccomplex_t IT_0361 = IT_0190 + IT_0360;
    const ccomplex_t IT_0362 = 4*IT_0361;
    const ccomplex_t IT_0363 = IT_0020*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0364 = (-0.25)*IT_0055;
    const ccomplex_t IT_0365 = IT_0056 + IT_0364;
    const ccomplex_t IT_0366 = IT_0363*IT_0365;
    const ccomplex_t IT_0367 = (-4)*IT_0366;
    const ccomplex_t IT_0368 = (-4)*IT_0353;
    const ccomplex_t IT_0369 = (-4)*IT_0361;
    const ccomplex_t IT_0370 = (-2)*IT_0112 + 4*IT_0117 + (-2)*IT_0319;
    const ccomplex_t IT_0371 = IT_0111 + IT_0370;
    const ccomplex_t IT_0372 = IT_0020*IT_0371;
    const ccomplex_t IT_0373 = -conj(IT_0051);
    const ccomplex_t IT_0374 = (-4)*(IT_0064 + -conj(IT_0064))*(IT_0063 + 2
      *IT_0134) + 32*IT_0008*(IT_0110 + 0.25*IT_0329) + 4*conj(IT_0075)*IT_0346 
      + 8*(IT_0019 + -conj(IT_0019))*((-0.5)*IT_0072 + IT_0110 + IT_0155 +
       IT_0156 + IT_0162 + IT_0347) + 2*IT_0017*IT_0354 + IT_0181*IT_0362 +
       IT_0135*IT_0367 + IT_0018*IT_0368 + 2*IT_0180*IT_0369 + 8*IT_0334*
      (IT_0110 + IT_0323 + IT_0330 + IT_0333 + -IT_0372) + (-4)*(IT_0041 + 
      -IT_0173)*(IT_0051 + IT_0373);
    const ccomplex_t IT_0375 = (IT_0110 + (-2)*IT_0132)*(1 + 0.5*IT_0258);
    const ccomplex_t IT_0376 = IT_0000*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0377 = IT_0260*IT_0376;
    const ccomplex_t IT_0378 = 2*IT_0240;
    const ccomplex_t IT_0379 = -IT_0136;
    const ccomplex_t IT_0380 = IT_0137 + IT_0140 + IT_0188 + IT_0246 + IT_0378
       + IT_0379;
    const ccomplex_t IT_0381 = 0.5*IT_0380;
    const ccomplex_t IT_0382 = -IT_0306;
    const ccomplex_t IT_0383 = IT_0233 + IT_0342 + IT_0382;
    const ccomplex_t IT_0384 = 0.5*IT_0383;
    const ccomplex_t IT_0385 = 0.5*IT_0096 + IT_0228 + IT_0381 + IT_0384;
    const ccomplex_t IT_0386 = IT_0096 + IT_0228 + IT_0230 + IT_0235 + IT_0315
       + IT_0383;
    const ccomplex_t IT_0387 = (-0.25)*IT_0018;
    const ccomplex_t IT_0388 = -IT_0243;
    const ccomplex_t IT_0389 = IT_0137 + IT_0140 + IT_0150 + IT_0188 + IT_0240
       + IT_0242 + IT_0246 + IT_0379 + IT_0388;
    const ccomplex_t IT_0390 = (-0.25)*IT_0181;
    const ccomplex_t IT_0391 = (-0.125)*IT_0134*IT_0163 + (-0.125)*IT_0063
      *IT_0174 + (-0.125)*conj(IT_0019)*IT_0256 + (-0.125)*conj(IT_0051)*IT_0298
       + conj(IT_0064)*(IT_0375 + IT_0377) + -conj(IT_0075)*IT_0385 + (0.5
      *IT_0026 + IT_0386)*IT_0387 + (0.5*IT_0186 + IT_0389)*IT_0390;
    const ccomplex_t IT_0392 = -(s_12 + (-0.5)*IT_0200)*(IT_0011 + -IT_0278);
    const ccomplex_t IT_0393 = (-2)*IT_0392;
    const ccomplex_t IT_0394 = (s_12 + (-0.5)*IT_0200)*(s_12 + -IT_0272);
    const ccomplex_t IT_0395 = (-2)*IT_0394;
    const ccomplex_t IT_0396 = -IT_0011;
    const ccomplex_t IT_0397 = IT_0278 + IT_0396;
    const ccomplex_t IT_0398 = IT_0052*IT_0397;
    const ccomplex_t IT_0399 = -IT_0272;
    const ccomplex_t IT_0400 = s_12 + IT_0399;
    const ccomplex_t IT_0401 = IT_0052*IT_0400;
    const ccomplex_t IT_0402 = -(s_12 + (-0.5)*IT_0200)*(IT_0012 + -IT_0264);
    const ccomplex_t IT_0403 = (-2)*IT_0402;
    const ccomplex_t IT_0404 = -IT_0012;
    const ccomplex_t IT_0405 = IT_0264 + IT_0404;
    const ccomplex_t IT_0406 = IT_0052*IT_0405;
    const ccomplex_t IT_0407 = s_13*s_23*IT_0031*IT_0277;
    const ccomplex_t IT_0408 = s_14*s_24*IT_0031*IT_0229;
    const ccomplex_t IT_0409 = s_12*IT_0011;
    const ccomplex_t IT_0410 = s_14*s_24*IT_0011*IT_0020;
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = s_13*s_23*IT_0011*IT_0020;
    const ccomplex_t IT_0413 = -IT_0412;
    const ccomplex_t IT_0414 = s_12*s_13*s_14*s_34*IT_0031;
    const ccomplex_t IT_0415 = -IT_0414;
    const ccomplex_t IT_0416 = IT_0407 + IT_0408 + IT_0409 + IT_0411 + IT_0413
       + IT_0415;
    const ccomplex_t IT_0417 = 2*IT_0416;
    const ccomplex_t IT_0418 = s_13*s_14*s_23*s_24*IT_0031;
    const ccomplex_t IT_0419 = IT_0011*IT_0012;
    const ccomplex_t IT_0420 = IT_0031*IT_0138*IT_0229;
    const ccomplex_t IT_0421 = IT_0012*IT_0020*IT_0229;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0011*IT_0020*IT_0138;
    const ccomplex_t IT_0424 = -IT_0423;
    const ccomplex_t IT_0425 = s_12*s_13*s_24*s_34*IT_0031;
    const ccomplex_t IT_0426 = -IT_0425;
    const ccomplex_t IT_0427 = IT_0418 + IT_0419 + IT_0420 + IT_0422 + IT_0424
       + IT_0426;
    const ccomplex_t IT_0428 = 2*IT_0427;
    const ccomplex_t IT_0429 = IT_0052*IT_0282;
    const ccomplex_t IT_0430 = 2*IT_0429;
    const ccomplex_t IT_0431 = IT_0052*IT_0274;
    const ccomplex_t IT_0432 = (-2)*IT_0431;
    const ccomplex_t IT_0433 = IT_0031*IT_0241*IT_0277;
    const ccomplex_t IT_0434 = IT_0012*IT_0020*IT_0277;
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = IT_0011*IT_0020*IT_0241;
    const ccomplex_t IT_0437 = -IT_0436;
    const ccomplex_t IT_0438 = s_12*s_14*s_23*s_34*IT_0031;
    const ccomplex_t IT_0439 = -IT_0438;
    const ccomplex_t IT_0440 = IT_0418 + IT_0419 + IT_0433 + IT_0435 + IT_0437
       + IT_0439;
    const ccomplex_t IT_0441 = 2*IT_0440;
    const ccomplex_t IT_0442 = s_14*s_24*IT_0031*IT_0241;
    const ccomplex_t IT_0443 = s_13*s_23*IT_0031*IT_0138;
    const ccomplex_t IT_0444 = s_12*IT_0012;
    const ccomplex_t IT_0445 = s_14*s_24*IT_0012*IT_0020;
    const ccomplex_t IT_0446 = -IT_0445;
    const ccomplex_t IT_0447 = s_13*s_23*IT_0012*IT_0020;
    const ccomplex_t IT_0448 = -IT_0447;
    const ccomplex_t IT_0449 = s_12*s_23*s_24*s_34*IT_0031;
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = IT_0442 + IT_0443 + IT_0444 + IT_0446 + IT_0448
       + IT_0450;
    const ccomplex_t IT_0452 = 2*IT_0451;
    const ccomplex_t IT_0453 = IT_0052*IT_0268;
    const ccomplex_t IT_0454 = 2*IT_0453;
    const ccomplex_t IT_0455 = IT_0052*IT_0287;
    const ccomplex_t IT_0456 = (-2)*IT_0455;
    const ccomplex_t IT_0457 = (-0.25)*IT_0135*IT_0152 + (-0.25)*IT_0094
      *IT_0174 + (-0.25)*conj(IT_0019)*IT_0220 + (-4)*IT_0008*IT_0318 + IT_0334
      *IT_0346 + IT_0291*IT_0373 + (-2)*conj(IT_0064)*IT_0385 + conj(IT_0075)*
      (IT_0393 + (-2)*IT_0395 + (-3)*IT_0398 + 6*IT_0401 + IT_0403 + (-3)
      *IT_0406) + IT_0387*(IT_0417 + IT_0428 + IT_0430 + IT_0432) + IT_0390*
      (IT_0441 + IT_0452 + IT_0454 + IT_0456);
    const ccomplex_t IT_0458 = 0.5*conj(IT_0051);
    const ccomplex_t IT_0459 = IT_0303*IT_0458;
    const ccomplex_t IT_0460 = s_14*s_24;
    const ccomplex_t IT_0461 = (-2)*IT_0460;
    const ccomplex_t IT_0462 = IT_0110 + IT_0461;
    const ccomplex_t IT_0463 = IT_0020*IT_0462;
    const ccomplex_t IT_0464 = -(IT_0011 + -IT_0276)*(s_12 + (-0.5)*IT_0463);
    const ccomplex_t IT_0465 = (-2)*IT_0464;
    const ccomplex_t IT_0466 = (-2)*IT_0465;
    const ccomplex_t IT_0467 = 0.5*conj(IT_0017);
    const ccomplex_t IT_0468 = IT_0276 + IT_0396;
    const ccomplex_t IT_0469 = IT_0052*IT_0468;
    const ccomplex_t IT_0470 = 6*IT_0469;
    const ccomplex_t IT_0471 = -IT_0271;
    const ccomplex_t IT_0472 = s_12 + IT_0471;
    const ccomplex_t IT_0473 = IT_0052*IT_0472;
    const ccomplex_t IT_0474 = (-6)*IT_0473;
    const ccomplex_t IT_0475 = 0.5*conj(IT_0180);
    const ccomplex_t IT_0476 = (s_12 + -IT_0271)*(s_12 + (-0.5)*IT_0463);
    const ccomplex_t IT_0477 = (-2)*IT_0476;
    const ccomplex_t IT_0478 = 2*IT_0477;
    const ccomplex_t IT_0479 = conj(IT_0075)*IT_0417;
    const ccomplex_t IT_0480 = 0.5*IT_0479;
    const ccomplex_t IT_0481 = conj(IT_0075)*IT_0428;
    const ccomplex_t IT_0482 = 0.5*IT_0481;
    const ccomplex_t IT_0483 = conj(IT_0019)*IT_0239;
    const ccomplex_t IT_0484 = 0.5*IT_0483;
    const ccomplex_t IT_0485 = IT_0026*conj(IT_0064);
    const ccomplex_t IT_0486 = 0.5*IT_0485;
    const ccomplex_t IT_0487 = conj(IT_0019)*IT_0026;
    const ccomplex_t IT_0488 = 0.5*IT_0487;
    const ccomplex_t IT_0489 = conj(IT_0075)*IT_0430;
    const ccomplex_t IT_0490 = 0.5*IT_0489;
    const ccomplex_t IT_0491 = IT_0108*IT_0174;
    const ccomplex_t IT_0492 = 0.5*IT_0491;
    const ccomplex_t IT_0493 = conj(IT_0051)*IT_0304;
    const ccomplex_t IT_0494 = 0.5*IT_0493;
    const ccomplex_t IT_0495 = conj(IT_0075)*IT_0432;
    const ccomplex_t IT_0496 = 0.5*IT_0495;
    const ccomplex_t IT_0497 = conj(IT_0064)*IT_0386;
    const ccomplex_t IT_0498 = 2*IT_0017*(IT_0459 + IT_0466*IT_0467 + IT_0467
      *IT_0470 + IT_0474*IT_0475 + IT_0475*IT_0478 + IT_0480 + IT_0482 + IT_0484
       + IT_0486 + IT_0488 + IT_0490 + IT_0492 + IT_0494 + IT_0496 + IT_0497);
    const ccomplex_t IT_0499 = -(IT_0012 + -IT_0263)*(s_12 + (-0.5)*IT_0463);
    const ccomplex_t IT_0500 = (-2)*IT_0499;
    const ccomplex_t IT_0501 = (-2)*IT_0500;
    const ccomplex_t IT_0502 = IT_0263 + IT_0404;
    const ccomplex_t IT_0503 = IT_0052*IT_0502;
    const ccomplex_t IT_0504 = 6*IT_0503;
    const ccomplex_t IT_0505 = conj(IT_0075)*IT_0441;
    const ccomplex_t IT_0506 = 0.5*IT_0505;
    const ccomplex_t IT_0507 = conj(IT_0075)*IT_0452;
    const ccomplex_t IT_0508 = 0.5*IT_0507;
    const ccomplex_t IT_0509 = conj(IT_0019)*IT_0247;
    const ccomplex_t IT_0510 = 0.5*IT_0509;
    const ccomplex_t IT_0511 = conj(IT_0051)*IT_0301;
    const ccomplex_t IT_0512 = 0.5*IT_0511;
    const ccomplex_t IT_0513 = conj(IT_0064)*IT_0186;
    const ccomplex_t IT_0514 = 0.5*IT_0513;
    const ccomplex_t IT_0515 = conj(IT_0019)*IT_0186;
    const ccomplex_t IT_0516 = 0.5*IT_0515;
    const ccomplex_t IT_0517 = conj(IT_0075)*IT_0454;
    const ccomplex_t IT_0518 = 0.5*IT_0517;
    const ccomplex_t IT_0519 = IT_0109*IT_0197;
    const ccomplex_t IT_0520 = 0.5*IT_0519;
    const ccomplex_t IT_0521 = conj(IT_0075)*IT_0456;
    const ccomplex_t IT_0522 = 0.5*IT_0521;
    const ccomplex_t IT_0523 = conj(IT_0064)*IT_0389;
    const ccomplex_t IT_0524 = 2*IT_0180*(IT_0300*IT_0458 + IT_0467*IT_0474 +
       IT_0467*IT_0478 + IT_0475*IT_0501 + IT_0475*IT_0504 + IT_0506 + IT_0508 +
       IT_0510 + IT_0512 + IT_0514 + IT_0516 + IT_0518 + IT_0520 + IT_0522 +
       IT_0523);
    const ccomplex_t IT_0525 = IT_0018*IT_0027 + (IT_0041*conj(IT_0051) +
       IT_0063*conj(IT_0064) + conj(IT_0019)*IT_0072 + conj(IT_0075)*IT_0094 +
       IT_0018*IT_0108)*IT_0109 + conj(IT_0064)*IT_0134*IT_0135 + (conj(IT_0075)
      *IT_0152 + conj(IT_0019)*(IT_0110 + IT_0153 + IT_0155 + IT_0157 + IT_0162)
      )*IT_0163 + conj(IT_0051)*IT_0173*IT_0174 + IT_0181*(IT_0187 + IT_0198) +
       4*IT_0019*IT_0257 + 8*IT_0051*IT_0305 + (-16)*IT_0008*(conj(IT_0075)
      *IT_0318 + (-0.5)*IT_0008*(IT_0110 + -IT_0322 + IT_0323 + IT_0330 +
       IT_0333)) + IT_0334*IT_0374 + (-8)*IT_0064*IT_0391 + (-4)*IT_0075*IT_0457
       + 2*IT_0498 + 2*IT_0524;
    return create_ccomplex_return(IT_0525);
}

