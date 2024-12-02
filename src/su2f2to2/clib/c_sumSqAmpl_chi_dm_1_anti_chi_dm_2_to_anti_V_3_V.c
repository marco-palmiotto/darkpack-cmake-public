#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_1_anti_chi_dm_2_to_anti_V_3_V.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_1_anti_chi_dm_2_to_anti_V_3_V(
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
    const ccomplex_t IT_0000 = pow(m_V, -2);
    const ccomplex_t IT_0001 = pow(m_V, 2);
    const ccomplex_t IT_0002 = s_23*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = pow(m_V, -4);
    const ccomplex_t IT_0005 = s_34*IT_0004;
    const ccomplex_t IT_0006 = s_24*IT_0001;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = s_24*s_34;
    const ccomplex_t IT_0009 = (-2)*IT_0008;
    const ccomplex_t IT_0010 = IT_0002 + IT_0009;
    const ccomplex_t IT_0011 = IT_0000*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0003 + (-0.5)*IT_0007 + (-0.5)*IT_0011;
    const ccomplex_t IT_0013 = s_23 + IT_0012;
    const ccomplex_t IT_0014 = IT_0013*m_chi_dm_1;
    const ccomplex_t IT_0015 = 2*IT_0014;
    const ccomplex_t IT_0016 = cpow((-2)*s_13 + IT_0001 + reg_prop, -1);
    const ccomplex_t IT_0017 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0020 = 0.5*IT_0019;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_chi_dm_1, 2);
    const ccomplex_t IT_0026 = pow(m_chi_dm_2, 2);
    const ccomplex_t IT_0027 = cpow((-2)*s_12 + m_V*(m_V + (0 + _Complex_I*-1)
      *Gamma_V) + -IT_0025 + -IT_0026 + -reg_prop, -1);
    const ccomplex_t IT_0028 = IT_0018*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = -IT_0019;
    const ccomplex_t IT_0031 = IT_0030*m_chi_dm_1;
    const ccomplex_t IT_0032 = (-2)*IT_0019;
    const ccomplex_t IT_0033 = IT_0032*m_chi_dm_2;
    const ccomplex_t IT_0034 = IT_0031 + IT_0033;
    const ccomplex_t IT_0035 = IT_0030*m_chi_dm_2;
    const ccomplex_t IT_0036 = IT_0032*m_chi_dm_1;
    const ccomplex_t IT_0037 = -IT_0035 + -IT_0036;
    const ccomplex_t IT_0038 = IT_0034 + IT_0037;
    const ccomplex_t IT_0039 = IT_0029*IT_0038;
    const ccomplex_t IT_0040 = 2*IT_0039;
    const ccomplex_t IT_0041 = IT_0015*IT_0024*IT_0040;
    const ccomplex_t IT_0042 = pow(s_13, 2);
    const ccomplex_t IT_0043 = IT_0000*IT_0042;
    const ccomplex_t IT_0044 = pow(s_14, 2);
    const ccomplex_t IT_0045 = IT_0000*IT_0044;
    const ccomplex_t IT_0046 = IT_0043 + IT_0045;
    const ccomplex_t IT_0047 = s_13*s_14*s_34*IT_0004;
    const ccomplex_t IT_0048 = -IT_0025 + -IT_0047;
    const ccomplex_t IT_0049 = IT_0046 + IT_0048;
    const ccomplex_t IT_0050 = IT_0049*m_chi_dm_2;
    const ccomplex_t IT_0051 = 2*IT_0019;
    const ccomplex_t IT_0052 = IT_0029*(IT_0019 + -IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0040*IT_0050*conj(IT_0053);
    const ccomplex_t IT_0055 = s_14*s_23*s_34*IT_0004;
    const ccomplex_t IT_0056 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0057 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0058 = -IT_0056 + -IT_0057;
    const ccomplex_t IT_0059 = s_12 + IT_0055 + IT_0058;
    const ccomplex_t IT_0060 = IT_0059*m_chi_dm_1;
    const ccomplex_t IT_0061 = IT_0040*conj(IT_0053)*IT_0060;
    const ccomplex_t IT_0062 = s_34*IT_0000*m_chi_dm_2;
    const ccomplex_t IT_0063 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = s_14 + IT_0064;
    const ccomplex_t IT_0066 = IT_0062*IT_0065;
    const ccomplex_t IT_0067 = IT_0019*IT_0029;
    const ccomplex_t IT_0068 = IT_0040*IT_0066*conj(IT_0067);
    const ccomplex_t IT_0069 = s_13*IT_0001;
    const ccomplex_t IT_0070 = IT_0000*IT_0069;
    const ccomplex_t IT_0071 = s_14*IT_0001;
    const ccomplex_t IT_0072 = IT_0005*IT_0071;
    const ccomplex_t IT_0073 = s_14*s_34;
    const ccomplex_t IT_0074 = (-2)*IT_0073;
    const ccomplex_t IT_0075 = IT_0069 + IT_0074;
    const ccomplex_t IT_0076 = IT_0000*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0070 + (-0.5)*IT_0072 + (-0.5)*IT_0076;
    const ccomplex_t IT_0078 = s_13 + IT_0077;
    const ccomplex_t IT_0079 = IT_0078*m_chi_dm_2;
    const ccomplex_t IT_0080 = 2*IT_0079;
    const ccomplex_t IT_0081 = 2*conj(IT_0039);
    const ccomplex_t IT_0082 = IT_0024*IT_0080*IT_0081;
    const ccomplex_t IT_0083 = IT_0050*IT_0053*IT_0081;
    const ccomplex_t IT_0084 = IT_0053*IT_0060*IT_0081;
    const ccomplex_t IT_0085 = IT_0066*IT_0067*IT_0081;
    const ccomplex_t IT_0086 = pow(s_34, 2);
    const ccomplex_t IT_0087 = IT_0004*IT_0086;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = 1 + IT_0088;
    const ccomplex_t IT_0090 = s_12*IT_0089;
    const ccomplex_t IT_0091 = 4*IT_0039*conj(IT_0039)*IT_0090;
    const ccomplex_t IT_0092 = s_23*m_chi_dm_1;
    const ccomplex_t IT_0093 = IT_0089*IT_0092;
    const ccomplex_t IT_0094 = IT_0029*IT_0032;
    const ccomplex_t IT_0095 = 4*conj(IT_0039)*IT_0093*IT_0094;
    const ccomplex_t IT_0096 = m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0097 = s_34*IT_0001;
    const ccomplex_t IT_0098 = IT_0005*IT_0097;
    const ccomplex_t IT_0099 = pow(m_V, 4);
    const ccomplex_t IT_0100 = (-2)*IT_0086;
    const ccomplex_t IT_0101 = IT_0099 + IT_0100;
    const ccomplex_t IT_0102 = IT_0000*IT_0101;
    const ccomplex_t IT_0103 = (-0.333333333333333)*IT_0098 + (
      -0.333333333333333)*IT_0102;
    const ccomplex_t IT_0104 = IT_0001 + IT_0103;
    const ccomplex_t IT_0105 = IT_0096*IT_0104;
    const ccomplex_t IT_0106 = 3*IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0001 + reg_prop, -1);
    const ccomplex_t IT_0108 = (-0.5)*IT_0019;
    const ccomplex_t IT_0109 = IT_0018*IT_0108;
    const ccomplex_t IT_0110 = IT_0107*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = 4*IT_0111;
    const ccomplex_t IT_0113 = s_12*IT_0001;
    const ccomplex_t IT_0114 = s_12*IT_0099;
    const ccomplex_t IT_0115 = s_13*s_23*IT_0001;
    const ccomplex_t IT_0116 = (-2)*IT_0115;
    const ccomplex_t IT_0117 = IT_0114 + IT_0116;
    const ccomplex_t IT_0118 = IT_0000*IT_0117;
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = s_14*s_23*s_34;
    const ccomplex_t IT_0121 = s_13*s_24*s_34;
    const ccomplex_t IT_0122 = s_12*IT_0086;
    const ccomplex_t IT_0123 = (-2)*IT_0115 + 2*IT_0120 + 2*IT_0121 + (-2)
      *IT_0122;
    const ccomplex_t IT_0124 = IT_0114 + IT_0123;
    const ccomplex_t IT_0125 = IT_0000*IT_0124;
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = s_14*s_23*IT_0001;
    const ccomplex_t IT_0128 = s_13*s_24*IT_0001;
    const ccomplex_t IT_0129 = IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = s_12*s_34*IT_0001;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = IT_0129 + IT_0131;
    const ccomplex_t IT_0133 = IT_0005*IT_0132;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = s_13*s_23;
    const ccomplex_t IT_0136 = (-2)*IT_0135;
    const ccomplex_t IT_0137 = IT_0113 + IT_0119 + IT_0126 + IT_0134 + IT_0136;
    const ccomplex_t IT_0138 = 8*IT_0111;
    const ccomplex_t IT_0139 = s_12*s_23;
    const ccomplex_t IT_0140 = s_14*s_23*s_24*IT_0000;
    const ccomplex_t IT_0141 = pow(s_24, 2);
    const ccomplex_t IT_0142 = s_13*IT_0000*IT_0141;
    const ccomplex_t IT_0143 = s_12*s_24*s_34*IT_0000;
    const ccomplex_t IT_0144 = (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = s_14*s_23*s_24;
    const ccomplex_t IT_0146 = s_13*IT_0141;
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = s_12*s_24*s_34;
    const ccomplex_t IT_0149 = s_13*IT_0001*IT_0026;
    const ccomplex_t IT_0150 = -IT_0148 + -IT_0149;
    const ccomplex_t IT_0151 = IT_0147 + IT_0150;
    const ccomplex_t IT_0152 = IT_0000*IT_0151;
    const ccomplex_t IT_0153 = (-0.25)*IT_0152;
    const ccomplex_t IT_0154 = s_13*IT_0026;
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = IT_0139 + (-0.5)*IT_0140 + 0.5*IT_0142 +
       IT_0144 + IT_0153 + IT_0155;
    const ccomplex_t IT_0157 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0158 = -IT_0073;
    const ccomplex_t IT_0159 = IT_0069 + IT_0158;
    const ccomplex_t IT_0160 = IT_0000*IT_0159;
    const ccomplex_t IT_0161 = -IT_0157 + (-0.25)*IT_0160;
    const ccomplex_t IT_0162 = s_13 + IT_0161;
    const ccomplex_t IT_0163 = IT_0096*IT_0162;
    const ccomplex_t IT_0164 = (-4)*IT_0163;
    const ccomplex_t IT_0165 = s_24*s_34*IT_0000;
    const ccomplex_t IT_0166 = -IT_0008;
    const ccomplex_t IT_0167 = IT_0002 + IT_0166;
    const ccomplex_t IT_0168 = IT_0000*IT_0167;
    const ccomplex_t IT_0169 = -IT_0165 + (-0.25)*IT_0168;
    const ccomplex_t IT_0170 = s_23 + IT_0169;
    const ccomplex_t IT_0171 = IT_0096*IT_0170;
    const ccomplex_t IT_0172 = (-4)*IT_0171;
    const ccomplex_t IT_0173 = IT_0164 + IT_0172;
    const ccomplex_t IT_0174 = -IT_0142;
    const ccomplex_t IT_0175 = -IT_0143;
    const ccomplex_t IT_0176 = (-0.5)*IT_0152;
    const ccomplex_t IT_0177 = IT_0140 + IT_0154 + IT_0174 + IT_0175 + IT_0176;
    const ccomplex_t IT_0178 = (-2)*IT_0177;
    const ccomplex_t IT_0179 = s_12*s_13;
    const ccomplex_t IT_0180 = s_23*IT_0000*IT_0044;
    const ccomplex_t IT_0181 = s_13*s_14*s_24*IT_0000;
    const ccomplex_t IT_0182 = s_23*IT_0044;
    const ccomplex_t IT_0183 = s_13*s_14*s_24;
    const ccomplex_t IT_0184 = IT_0182 + IT_0183;
    const ccomplex_t IT_0185 = s_12*s_14*s_34;
    const ccomplex_t IT_0186 = s_23*IT_0001*IT_0025;
    const ccomplex_t IT_0187 = -IT_0185 + -IT_0186;
    const ccomplex_t IT_0188 = IT_0184 + IT_0187;
    const ccomplex_t IT_0189 = IT_0000*IT_0188;
    const ccomplex_t IT_0190 = s_12*s_14*s_34*IT_0000;
    const ccomplex_t IT_0191 = (-0.5)*IT_0190;
    const ccomplex_t IT_0192 = s_23*IT_0025;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0179 + 0.5*IT_0180 + (-0.5)*IT_0181 + (-0.25
      )*IT_0189 + IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = (-4)*IT_0194;
    const ccomplex_t IT_0196 = IT_0173 + IT_0178 + IT_0195;
    const ccomplex_t IT_0197 = IT_0113 + IT_0136;
    const ccomplex_t IT_0198 = IT_0000*IT_0197;
    const ccomplex_t IT_0199 = (s_13 + -IT_0157)*(s_12 + (-0.5)*IT_0198);
    const ccomplex_t IT_0200 = (-2)*IT_0199;
    const ccomplex_t IT_0201 = 4*IT_0200;
    const ccomplex_t IT_0202 = (s_23 + -IT_0165)*(s_12 + (-0.5)*IT_0198);
    const ccomplex_t IT_0203 = (-2)*IT_0202;
    const ccomplex_t IT_0204 = 4*IT_0203;
    const ccomplex_t IT_0205 = -IT_0165;
    const ccomplex_t IT_0206 = s_23 + IT_0205;
    const ccomplex_t IT_0207 = IT_0096*IT_0206;
    const ccomplex_t IT_0208 = (-12)*IT_0207;
    const ccomplex_t IT_0209 = -IT_0157;
    const ccomplex_t IT_0210 = s_13 + IT_0209;
    const ccomplex_t IT_0211 = IT_0096*IT_0210;
    const ccomplex_t IT_0212 = (-12)*IT_0211;
    const ccomplex_t IT_0213 = IT_0201 + IT_0204 + IT_0208 + IT_0212;
    const ccomplex_t IT_0214 = 2*IT_0030;
    const ccomplex_t IT_0215 = IT_0029*IT_0214;
    const ccomplex_t IT_0216 = 2*IT_0023 + IT_0215;
    const ccomplex_t IT_0217 = s_13*s_23*IT_0004*IT_0044;
    const ccomplex_t IT_0218 = s_14*s_24*IT_0004*IT_0042;
    const ccomplex_t IT_0219 = s_12*s_13*s_14*s_34*IT_0004;
    const ccomplex_t IT_0220 = s_12*IT_0025;
    const ccomplex_t IT_0221 = s_14*s_24*IT_0000*IT_0025;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = s_13*s_23*IT_0000*IT_0025;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0217 + IT_0218 + -IT_0219 + IT_0220 +
       IT_0222 + IT_0224;
    const ccomplex_t IT_0226 = 4*IT_0225;
    const ccomplex_t IT_0227 = s_13*s_14*s_23*s_24*IT_0004;
    const ccomplex_t IT_0228 = IT_0025*IT_0026;
    const ccomplex_t IT_0229 = IT_0004*IT_0042*IT_0141;
    const ccomplex_t IT_0230 = IT_0000*IT_0026*IT_0042;
    const ccomplex_t IT_0231 = -IT_0230;
    const ccomplex_t IT_0232 = IT_0000*IT_0025*IT_0141;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = s_12*s_13*s_24*s_34*IT_0004;
    const ccomplex_t IT_0235 = -IT_0234;
    const ccomplex_t IT_0236 = IT_0227 + IT_0228 + IT_0229 + IT_0231 + IT_0233
       + IT_0235;
    const ccomplex_t IT_0237 = 4*IT_0236;
    const ccomplex_t IT_0238 = IT_0049*IT_0096;
    const ccomplex_t IT_0239 = 4*IT_0238;
    const ccomplex_t IT_0240 = s_13*s_24*s_34*IT_0004;
    const ccomplex_t IT_0241 = s_12 + IT_0058 + IT_0240;
    const ccomplex_t IT_0242 = IT_0096*IT_0241;
    const ccomplex_t IT_0243 = (-4)*IT_0242;
    const ccomplex_t IT_0244 = IT_0226 + IT_0237 + IT_0239 + IT_0243;
    const ccomplex_t IT_0245 = (-2)*IT_0111 + IT_0215;
    const ccomplex_t IT_0246 = pow(s_23, 2);
    const ccomplex_t IT_0247 = IT_0004*IT_0044*IT_0246;
    const ccomplex_t IT_0248 = IT_0000*IT_0026*IT_0044;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0000*IT_0025*IT_0246;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = s_12*s_14*s_23*s_34*IT_0004;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = IT_0227 + IT_0228 + IT_0247 + IT_0249 + IT_0251
       + IT_0253;
    const ccomplex_t IT_0255 = 4*IT_0254;
    const ccomplex_t IT_0256 = s_14*s_24*IT_0004*IT_0246;
    const ccomplex_t IT_0257 = s_13*s_23*IT_0004*IT_0141;
    const ccomplex_t IT_0258 = s_12*IT_0026;
    const ccomplex_t IT_0259 = s_14*s_24*IT_0000*IT_0026;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = s_13*s_23*IT_0000*IT_0026;
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = s_12*s_23*s_24*s_34*IT_0004;
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = IT_0256 + IT_0257 + IT_0258 + IT_0260 + IT_0262
       + IT_0264;
    const ccomplex_t IT_0266 = 4*IT_0265;
    const ccomplex_t IT_0267 = IT_0000*IT_0246;
    const ccomplex_t IT_0268 = IT_0000*IT_0141;
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = s_23*s_24*s_34*IT_0004;
    const ccomplex_t IT_0271 = -IT_0026 + -IT_0270;
    const ccomplex_t IT_0272 = IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = IT_0096*IT_0272;
    const ccomplex_t IT_0274 = 4*IT_0273;
    const ccomplex_t IT_0275 = IT_0059*IT_0096;
    const ccomplex_t IT_0276 = (-4)*IT_0275;
    const ccomplex_t IT_0277 = IT_0255 + IT_0266 + IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = s_13*s_14*s_23*s_34*IT_0004;
    const ccomplex_t IT_0279 = s_23*IT_0000*IT_0042;
    const ccomplex_t IT_0280 = (-2)*IT_0279;
    const ccomplex_t IT_0281 = IT_0179 + IT_0280;
    const ccomplex_t IT_0282 = -IT_0181;
    const ccomplex_t IT_0283 = -IT_0180;
    const ccomplex_t IT_0284 = IT_0281 + IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = 4*IT_0192 + 8*IT_0278 + 4*IT_0284;
    const ccomplex_t IT_0286 = s_13*s_23*s_24*s_34*IT_0004;
    const ccomplex_t IT_0287 = s_13*IT_0000*IT_0246;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = IT_0139 + IT_0288;
    const ccomplex_t IT_0290 = -IT_0140;
    const ccomplex_t IT_0291 = IT_0174 + IT_0289 + IT_0290;
    const ccomplex_t IT_0292 = 4*IT_0154 + 8*IT_0286 + 4*IT_0291;
    const ccomplex_t IT_0293 = IT_0285 + IT_0292;
    const ccomplex_t IT_0294 = -IT_0192;
    const ccomplex_t IT_0295 = 0.5*IT_0189;
    const ccomplex_t IT_0296 = IT_0180 + IT_0190 + IT_0282 + IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = IT_0000*IT_0086;
    const ccomplex_t IT_0298 = -IT_0086;
    const ccomplex_t IT_0299 = IT_0099 + IT_0298;
    const ccomplex_t IT_0300 = IT_0000*IT_0299;
    const ccomplex_t IT_0301 = -IT_0001 + 0.25*IT_0300;
    const ccomplex_t IT_0302 = IT_0297 + IT_0301;
    const ccomplex_t IT_0303 = IT_0096*IT_0302;
    const ccomplex_t IT_0304 = 4*IT_0303;
    const ccomplex_t IT_0305 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0306 = IT_0114 + IT_0120 + IT_0121;
    const ccomplex_t IT_0307 = (-2)*IT_0115 + -IT_0122;
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = IT_0000*IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*IT_0309;
    const ccomplex_t IT_0311 = IT_0113 + IT_0310;
    const ccomplex_t IT_0312 = s_12*IT_0000*IT_0086;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = IT_0305 + IT_0311 + IT_0313 + IT_0315;
    const ccomplex_t IT_0317 = (-8)*IT_0111;
    const ccomplex_t IT_0318 = IT_0241*m_chi_dm_1;
    const ccomplex_t IT_0319 = IT_0050 + IT_0318;
    const ccomplex_t IT_0320 = IT_0179 + IT_0294;
    const ccomplex_t IT_0321 = s_12*s_13*IT_0001;
    const ccomplex_t IT_0322 = -IT_0185 + (-0.5)*IT_0186 + 0.5*IT_0321;
    const ccomplex_t IT_0323 = IT_0182 + IT_0322;
    const ccomplex_t IT_0324 = IT_0000*IT_0323;
    const ccomplex_t IT_0325 = 2*IT_0324;
    const ccomplex_t IT_0326 = (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = IT_0320 + IT_0326;
    const ccomplex_t IT_0328 = s_13*m_chi_dm_2;
    const ccomplex_t IT_0329 = IT_0089*IT_0328;
    const ccomplex_t IT_0330 = (-4)*IT_0039*conj(IT_0094)*IT_0329;
    const ccomplex_t IT_0331 = (-4)*conj(IT_0039);
    const ccomplex_t IT_0332 = IT_0094*IT_0329*IT_0331;
    const ccomplex_t IT_0333 = -IT_0026;
    const ccomplex_t IT_0334 = IT_0268 + IT_0333;
    const ccomplex_t IT_0335 = IT_0096*IT_0334;
    const ccomplex_t IT_0336 = -(IT_0025 + -IT_0045)*(s_12 + (-0.5)*IT_0198);
    const ccomplex_t IT_0337 = (-2)*IT_0336;
    const ccomplex_t IT_0338 = (s_12 + -IT_0057)*(s_12 + (-0.5)*IT_0198);
    const ccomplex_t IT_0339 = (-2)*IT_0338;
    const ccomplex_t IT_0340 = -IT_0025;
    const ccomplex_t IT_0341 = IT_0045 + IT_0340;
    const ccomplex_t IT_0342 = IT_0096*IT_0341;
    const ccomplex_t IT_0343 = -IT_0057;
    const ccomplex_t IT_0344 = s_12 + IT_0343;
    const ccomplex_t IT_0345 = IT_0096*IT_0344;
    const ccomplex_t IT_0346 = -(s_12 + (-0.5)*IT_0198)*(IT_0026 + -IT_0268);
    const ccomplex_t IT_0347 = (-2)*IT_0346;
    const ccomplex_t IT_0348 = (-4)*IT_0111;
    const ccomplex_t IT_0349 = 16*IT_0111*IT_0156 + (-4)*IT_0024*IT_0196 +
       conj(IT_0067)*IT_0213 + conj(IT_0216)*IT_0244 + conj(IT_0245)*IT_0277 +
       conj(IT_0094)*IT_0293 + IT_0296*IT_0317 + conj(IT_0053)*(12*IT_0335 + (-4
      )*IT_0337 + 8*IT_0339 + 12*IT_0342 + (-24)*IT_0345 + (-4)*IT_0347) +
       IT_0173*IT_0348;
    const ccomplex_t IT_0350 = (1 + 0.5*IT_0087)*(IT_0113 + (-2)*IT_0135);
    const ccomplex_t IT_0351 = IT_0001*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0352 = IT_0089*IT_0351;
    const ccomplex_t IT_0353 = s_24*s_34*IT_0004*IT_0042;
    const ccomplex_t IT_0354 = s_12*s_13*IT_0004*IT_0086;
    const ccomplex_t IT_0355 = -IT_0354;
    const ccomplex_t IT_0356 = s_34*IT_0000*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0357 = IT_0065*IT_0356;
    const ccomplex_t IT_0358 = 4*IT_0357;
    const ccomplex_t IT_0359 = 0.25*IT_0358;
    const ccomplex_t IT_0360 = IT_0190 + IT_0192 + IT_0278 + IT_0284 + IT_0353
       + IT_0355 + IT_0359;
    const ccomplex_t IT_0361 = s_14*s_34*IT_0004*IT_0246;
    const ccomplex_t IT_0362 = s_12*s_23*IT_0004*IT_0086;
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = s_24 + IT_0365;
    const ccomplex_t IT_0367 = IT_0356*IT_0366;
    const ccomplex_t IT_0368 = 4*IT_0367;
    const ccomplex_t IT_0369 = 0.25*IT_0368;
    const ccomplex_t IT_0370 = IT_0143 + IT_0154 + IT_0286 + IT_0291 + IT_0361
       + IT_0363 + IT_0369;
    const ccomplex_t IT_0371 = s_34*IT_0000;
    const ccomplex_t IT_0372 = s_14*s_23;
    const ccomplex_t IT_0373 = s_13*s_24;
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = s_13*s_23*s_34*IT_0000;
    const ccomplex_t IT_0376 = (-2)*IT_0375;
    const ccomplex_t IT_0377 = IT_0374 + IT_0376;
    const ccomplex_t IT_0378 = IT_0371*IT_0377;
    const ccomplex_t IT_0379 = (-4)*IT_0378;
    const ccomplex_t IT_0380 = (s_12 + (-0.5)*IT_0198)*(IT_0001 + -IT_0297);
    const ccomplex_t IT_0381 = (-2)*IT_0380;
    const ccomplex_t IT_0382 = -IT_0297;
    const ccomplex_t IT_0383 = IT_0001 + IT_0382;
    const ccomplex_t IT_0384 = IT_0096*IT_0383;
    const ccomplex_t IT_0385 = s_24*s_34*IT_0000*IT_0025;
    const ccomplex_t IT_0386 = 4*IT_0192 + 4*IT_0278 + 4*IT_0281 + 4*IT_0353 +
       (-4)*IT_0354 + (-4)*IT_0385;
    const ccomplex_t IT_0387 = IT_0358 + IT_0386;
    const ccomplex_t IT_0388 = s_14*s_34*IT_0000*IT_0026;
    const ccomplex_t IT_0389 = 4*IT_0154 + 4*IT_0286 + 4*IT_0289 + 4*IT_0361 +
       (-4)*IT_0362 + (-4)*IT_0388;
    const ccomplex_t IT_0390 = IT_0368 + IT_0389;
    const ccomplex_t IT_0391 = 4*IT_0067*(0.25*IT_0138*IT_0316 + 0.25*IT_0304
      *IT_0348 + 0.25*conj(IT_0094)*IT_0379 + conj(IT_0067)*(IT_0381 + (-3)
      *IT_0384) + 0.25*conj(IT_0216)*IT_0387 + 0.25*conj(IT_0245)*IT_0390);
    const ccomplex_t IT_0392 = -IT_0154;
    const ccomplex_t IT_0393 = IT_0139 + IT_0392;
    const ccomplex_t IT_0394 = s_12*s_23*IT_0001;
    const ccomplex_t IT_0395 = s_14*s_34*IT_0026;
    const ccomplex_t IT_0396 = 0.5*IT_0149 + (-0.5)*IT_0394 + -IT_0395;
    const ccomplex_t IT_0397 = IT_0145 + IT_0396;
    const ccomplex_t IT_0398 = IT_0000*IT_0397;
    const ccomplex_t IT_0399 = (-2)*IT_0398;
    const ccomplex_t IT_0400 = (-0.5)*IT_0399;
    const ccomplex_t IT_0401 = IT_0393 + IT_0400;
    const ccomplex_t IT_0402 = IT_0059*m_chi_dm_2;
    const ccomplex_t IT_0403 = (-2)*IT_0039;
    const ccomplex_t IT_0404 = IT_0024*IT_0080*IT_0403;
    const ccomplex_t IT_0405 = IT_0241*m_chi_dm_2;
    const ccomplex_t IT_0406 = conj(IT_0053)*IT_0403*IT_0405;
    const ccomplex_t IT_0407 = IT_0272*m_chi_dm_1;
    const ccomplex_t IT_0408 = IT_0403*IT_0407;
    const ccomplex_t IT_0409 = conj(IT_0053)*IT_0408;
    const ccomplex_t IT_0410 = conj(IT_0245)*IT_0408;
    const ccomplex_t IT_0411 = s_34*IT_0000*m_chi_dm_1;
    const ccomplex_t IT_0412 = IT_0366*IT_0411;
    const ccomplex_t IT_0413 = conj(IT_0067)*IT_0403*IT_0412;
    const ccomplex_t IT_0414 = (-2)*conj(IT_0039);
    const ccomplex_t IT_0415 = IT_0015*IT_0024*IT_0414;
    const ccomplex_t IT_0416 = IT_0053*IT_0405*IT_0414;
    const ccomplex_t IT_0417 = IT_0407*IT_0414;
    const ccomplex_t IT_0418 = IT_0053*IT_0417;
    const ccomplex_t IT_0419 = IT_0245*IT_0417;
    const ccomplex_t IT_0420 = IT_0067*IT_0412*IT_0414;
    const ccomplex_t IT_0421 = 2*IT_0111;
    const ccomplex_t IT_0422 = IT_0089*IT_0096;
    const ccomplex_t IT_0423 = (-2)*IT_0111;
    const ccomplex_t IT_0424 = IT_0001*IT_0096;
    const ccomplex_t IT_0425 = 9*IT_0424;
    const ccomplex_t IT_0426 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0427 = (-2)*IT_0115 + 4*IT_0120 + (-2)*IT_0426;
    const ccomplex_t IT_0428 = IT_0114 + IT_0427;
    const ccomplex_t IT_0429 = IT_0000*IT_0428;
    const ccomplex_t IT_0430 = (-4)*IT_0135;
    const ccomplex_t IT_0431 = pow(m_V, 6);
    const ccomplex_t IT_0432 = s_12*IT_0431;
    const ccomplex_t IT_0433 = s_14*s_24*IT_0099;
    const ccomplex_t IT_0434 = (-2)*IT_0433;
    const ccomplex_t IT_0435 = IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = IT_0004*IT_0435;
    const ccomplex_t IT_0437 = 0.5*IT_0436;
    const ccomplex_t IT_0438 = s_14*s_24;
    const ccomplex_t IT_0439 = (-2)*IT_0438;
    const ccomplex_t IT_0440 = IT_0113 + IT_0439;
    const ccomplex_t IT_0441 = IT_0000*IT_0440;
    const ccomplex_t IT_0442 = -(IT_0025 + -IT_0043)*(s_12 + (-0.5)*IT_0441);
    const ccomplex_t IT_0443 = (-2)*IT_0442;
    const ccomplex_t IT_0444 = IT_0043 + IT_0340;
    const ccomplex_t IT_0445 = IT_0096*IT_0444;
    const ccomplex_t IT_0446 = (s_12 + -IT_0056)*(s_12 + (-0.5)*IT_0441);
    const ccomplex_t IT_0447 = (-2)*IT_0446;
    const ccomplex_t IT_0448 = 4*IT_0447;
    const ccomplex_t IT_0449 = -IT_0056;
    const ccomplex_t IT_0450 = s_12 + IT_0449;
    const ccomplex_t IT_0451 = IT_0096*IT_0450;
    const ccomplex_t IT_0452 = (-12)*IT_0451;
    const ccomplex_t IT_0453 = IT_0448 + IT_0452;
    const ccomplex_t IT_0454 = -conj(IT_0094);
    const ccomplex_t IT_0455 = -(IT_0026 + -IT_0267)*(s_12 + (-0.5)*IT_0441);
    const ccomplex_t IT_0456 = (-2)*IT_0455;
    const ccomplex_t IT_0457 = IT_0267 + IT_0333;
    const ccomplex_t IT_0458 = IT_0096*IT_0457;
    const ccomplex_t IT_0459 = s_24*s_34*IT_0025;
    const ccomplex_t IT_0460 = 0.5*IT_0186 + (-0.5)*IT_0321 + -IT_0459;
    const ccomplex_t IT_0461 = IT_0183 + IT_0460;
    const ccomplex_t IT_0462 = IT_0000*IT_0461;
    const ccomplex_t IT_0463 = (-2)*IT_0462;
    const ccomplex_t IT_0464 = (-0.5)*IT_0463;
    const ccomplex_t IT_0465 = IT_0320 + IT_0464;
    const ccomplex_t IT_0466 = -IT_0148 + (-0.5)*IT_0149 + 0.5*IT_0394;
    const ccomplex_t IT_0467 = IT_0146 + IT_0466;
    const ccomplex_t IT_0468 = IT_0000*IT_0467;
    const ccomplex_t IT_0469 = 2*IT_0468;
    const ccomplex_t IT_0470 = (-0.5)*IT_0469;
    const ccomplex_t IT_0471 = IT_0393 + IT_0470;
    const ccomplex_t IT_0472 = IT_0000*m_chi_dm_1*m_chi_dm_2;
    const ccomplex_t IT_0473 = (-0.25)*IT_0099;
    const ccomplex_t IT_0474 = IT_0086 + IT_0473;
    const ccomplex_t IT_0475 = IT_0472*IT_0474;
    const ccomplex_t IT_0476 = (-4)*IT_0475;
    const ccomplex_t IT_0477 = -IT_0305;
    const ccomplex_t IT_0478 = (-2)*IT_0115 + 4*IT_0121 + (-2)*IT_0426;
    const ccomplex_t IT_0479 = IT_0114 + IT_0478;
    const ccomplex_t IT_0480 = IT_0000*IT_0479;
    const ccomplex_t IT_0481 = 32*IT_0111*(IT_0113 + 0.25*IT_0436) + (-4)*
      (IT_0106 + 2*IT_0137)*(IT_0094 + IT_0454) + (-8)*(IT_0216 + -conj(IT_0216)
      )*IT_0465 + 8*(IT_0245 + -conj(IT_0245))*IT_0471 + IT_0138*IT_0476 + 8*
      (IT_0067 + -conj(IT_0067))*((-0.5)*IT_0304 + IT_0311 + IT_0313 + IT_0314 +
       IT_0477) + 8*IT_0024*(IT_0113 + 0.5*IT_0425 + IT_0430 + IT_0437 + 
      -IT_0480);
    const ccomplex_t IT_0482 = 2*IT_0041 + 2*IT_0054 + 2*IT_0061 + 2*IT_0068 +
       2*IT_0082 + 2*IT_0083 + 2*IT_0084 + 2*IT_0085 + 2*IT_0091 + 2*IT_0095 +
       conj(IT_0094)*(IT_0106*IT_0112 + IT_0137*IT_0138) + conj(IT_0053)*((-16)
      *IT_0111*IT_0156 + IT_0112*IT_0173 + 4*IT_0024*IT_0196 + IT_0067*IT_0213 +
       IT_0216*IT_0244 + IT_0245*IT_0277 + IT_0094*IT_0293 + IT_0138*IT_0296) +
       conj(IT_0067)*(IT_0112*IT_0304 + IT_0316*IT_0317) + 2*conj(IT_0216)*
      (IT_0040*IT_0319 + 0.5*IT_0317*IT_0327) + 2*IT_0330 + 2*IT_0332 + IT_0053
      *IT_0349 + (-8)*IT_0094*((-0.125)*IT_0137*IT_0317 + (-0.125)*IT_0106
      *IT_0348 + conj(IT_0094)*(IT_0350 + IT_0352) + (-0.5)*conj(IT_0216)
      *IT_0360 + (-0.5)*conj(IT_0245)*IT_0370 + (-0.125)*conj(IT_0067)*IT_0379) 
      + IT_0391 + conj(IT_0245)*(IT_0138*IT_0401 + 2*IT_0402*IT_0403) + 2
      *IT_0404 + 2*IT_0406 + 2*IT_0409 + 2*IT_0410 + 2*IT_0413 + 2*IT_0415 + 2
      *IT_0416 + 2*IT_0418 + 2*IT_0419 + 2*IT_0420 + conj(IT_0039)*(IT_0080
      *IT_0421 + IT_0015*IT_0423) + 8*IT_0039*(IT_0093*conj(IT_0094) + 0.125
      *IT_0015*IT_0421 + (-0.5)*(conj(IT_0039) + (-0.25)*IT_0331)*IT_0422 +
       0.125*IT_0080*IT_0423) + IT_0111*(IT_0080*(IT_0081 + IT_0403) + IT_0015*
      (IT_0040 + IT_0414) + 8*IT_0111*(IT_0113 + 0.5*IT_0425 + -IT_0429 +
       IT_0430 + IT_0437)) + (-4)*IT_0216*((-0.5)*IT_0081*IT_0319 + (-0.25)
      *IT_0138*IT_0327 + (-0.25)*conj(IT_0067)*IT_0387 + conj(IT_0216)*(IT_0443 
      + (-3)*IT_0445) + (-0.25)*conj(IT_0245)*IT_0453 + IT_0360*IT_0454) + (-4)
      *IT_0245*((-0.25)*conj(IT_0067)*IT_0390 + (-0.25)*IT_0317*IT_0401 + (-0.5)
      *IT_0402*IT_0414 + (-0.25)*conj(IT_0216)*IT_0453 + IT_0370*IT_0454 + conj
      (IT_0245)*(IT_0456 + (-3)*IT_0458)) + IT_0024*IT_0481;
    return create_ccomplex_return(IT_0482);
}

