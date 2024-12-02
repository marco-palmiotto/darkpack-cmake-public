#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_chi_dm_2_anti_chi_dm_2_to_anti_V_V.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_chi_dm_2_anti_chi_dm_2_to_anti_V_V(
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
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + -IT_0001 + IT_0002 +
       reg_prop, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1.4142135623731)*g_f;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = cpow(IT_0005, 2);
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*IT_0006*m_chi_dm_1;
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*IT_0006;
    const ccomplex_t IT_0010 = IT_0003*IT_0009;
    const ccomplex_t IT_0011 = IT_0010*m_chi_dm_2;
    const ccomplex_t IT_0012 = IT_0008 + -IT_0011;
    const ccomplex_t IT_0013 = pow(m_V, -2);
    const ccomplex_t IT_0014 = IT_0000*IT_0002;
    const ccomplex_t IT_0015 = pow(s_14, 2);
    const ccomplex_t IT_0016 = s_13*s_14*s_34*IT_0013;
    const ccomplex_t IT_0017 = (-0.5)*IT_0014 + 0.5*IT_0015 + -IT_0016;
    const ccomplex_t IT_0018 = pow(s_13, 2);
    const ccomplex_t IT_0019 = IT_0017 + IT_0018;
    const ccomplex_t IT_0020 = IT_0013*IT_0019;
    const ccomplex_t IT_0021 = IT_0002 + IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0015;
    const ccomplex_t IT_0023 = s_12*IT_0000;
    const ccomplex_t IT_0024 = s_14*s_24;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0023 + IT_0025;
    const ccomplex_t IT_0027 = IT_0013*IT_0026;
    const ccomplex_t IT_0028 = s_14*s_24*IT_0013;
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = s_12 + IT_0029;
    const ccomplex_t IT_0031 = -IT_0022 + (-0.5)*IT_0027 + 2*IT_0030;
    const ccomplex_t IT_0032 = IT_0021 + IT_0031;
    const ccomplex_t IT_0033 = IT_0032*m_chi_dm_2;
    const ccomplex_t IT_0034 = (-2)*IT_0033;
    const ccomplex_t IT_0035 = s_13*s_23;
    const ccomplex_t IT_0036 = s_13*s_24*s_34*IT_0013;
    const ccomplex_t IT_0037 = -IT_0024 + (-2)*IT_0035 + 2*IT_0036;
    const ccomplex_t IT_0038 = IT_0023 + IT_0037;
    const ccomplex_t IT_0039 = IT_0013*IT_0038;
    const ccomplex_t IT_0040 = pow(s_24, 2);
    const ccomplex_t IT_0041 = -IT_0014;
    const ccomplex_t IT_0042 = IT_0040 + IT_0041;
    const ccomplex_t IT_0043 = IT_0013*IT_0042;
    const ccomplex_t IT_0044 = IT_0013*IT_0040;
    const ccomplex_t IT_0045 = -IT_0002;
    const ccomplex_t IT_0046 = IT_0044 + IT_0045;
    const ccomplex_t IT_0047 = -IT_0028 + (-0.5)*IT_0039 + 0.5*IT_0043 + (-2)
      *IT_0046;
    const ccomplex_t IT_0048 = s_12 + IT_0047;
    const ccomplex_t IT_0049 = IT_0048*m_chi_dm_2;
    const ccomplex_t IT_0050 = (-2)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + (-2)*IT_0002 + -reg_prop + m_V
      *(m_V + (0 + _Complex_I*-1)*Gamma_V_3), -1);
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*g_f;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = IT_0052*IT_0055;
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = s_13*s_23*IT_0013;
    const ccomplex_t IT_0059 = IT_0013*IT_0018;
    const ccomplex_t IT_0060 = IT_0045 + IT_0059;
    const ccomplex_t IT_0061 = (-2)*IT_0024 + -IT_0035 + 2*IT_0036;
    const ccomplex_t IT_0062 = IT_0023 + IT_0061;
    const ccomplex_t IT_0063 = IT_0013*IT_0062;
    const ccomplex_t IT_0064 = IT_0018 + IT_0041;
    const ccomplex_t IT_0065 = IT_0013*IT_0064;
    const ccomplex_t IT_0066 = -IT_0058 + (-2)*IT_0060 + (-0.5)*IT_0063 + 0.5
      *IT_0065;
    const ccomplex_t IT_0067 = s_12 + IT_0066;
    const ccomplex_t IT_0068 = IT_0067*m_chi_dm_2;
    const ccomplex_t IT_0069 = (-2)*IT_0068;
    const ccomplex_t IT_0070 = pow(s_23, 2);
    const ccomplex_t IT_0071 = IT_0013*IT_0070;
    const ccomplex_t IT_0072 = -IT_0035;
    const ccomplex_t IT_0073 = IT_0023 + IT_0072;
    const ccomplex_t IT_0074 = IT_0013*IT_0073;
    const ccomplex_t IT_0075 = s_23*s_24*s_34*IT_0013;
    const ccomplex_t IT_0076 = -IT_0014 + 2*IT_0040 + (-2)*IT_0075;
    const ccomplex_t IT_0077 = IT_0070 + IT_0076;
    const ccomplex_t IT_0078 = IT_0013*IT_0077;
    const ccomplex_t IT_0079 = -IT_0058;
    const ccomplex_t IT_0080 = s_12 + IT_0079;
    const ccomplex_t IT_0081 = -IT_0071 + (-0.5)*IT_0074 + 0.5*IT_0078 + 2
      *IT_0080;
    const ccomplex_t IT_0082 = IT_0002 + IT_0081;
    const ccomplex_t IT_0083 = IT_0082*m_chi_dm_2;
    const ccomplex_t IT_0084 = (-2)*IT_0083;
    const ccomplex_t IT_0085 = (-2)*IT_0052;
    const ccomplex_t IT_0086 = IT_0055*IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = 1.125*IT_0002;
    const ccomplex_t IT_0089 = (-2)*IT_0035;
    const ccomplex_t IT_0090 = IT_0023 + IT_0089;
    const ccomplex_t IT_0091 = IT_0013*IT_0090;
    const ccomplex_t IT_0092 = (-2)*IT_0024;
    const ccomplex_t IT_0093 = IT_0023 + IT_0092;
    const ccomplex_t IT_0094 = IT_0013*IT_0093;
    const ccomplex_t IT_0095 = pow(m_V, -4);
    const ccomplex_t IT_0096 = pow(m_V, 4);
    const ccomplex_t IT_0097 = s_12*IT_0096;
    const ccomplex_t IT_0098 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0099 = s_13*s_24*s_34;
    const ccomplex_t IT_0100 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0101 = (-2)*IT_0098 + 4*IT_0099 + (-2)*IT_0100;
    const ccomplex_t IT_0102 = IT_0097 + IT_0101;
    const ccomplex_t IT_0103 = IT_0095*IT_0102;
    const ccomplex_t IT_0104 = (-4)*IT_0060;
    const ccomplex_t IT_0105 = IT_0065 + IT_0104;
    const ccomplex_t IT_0106 = IT_0105*m_chi_dm_2;
    const ccomplex_t IT_0107 = -IT_0058 + (-0.5)*IT_0063;
    const ccomplex_t IT_0108 = s_12 + IT_0107;
    const ccomplex_t IT_0109 = IT_0108*m_chi_dm_2;
    const ccomplex_t IT_0110 = (-2)*IT_0109;
    const ccomplex_t IT_0111 = 2*IT_0010;
    const ccomplex_t IT_0112 = s_13*IT_0000;
    const ccomplex_t IT_0113 = IT_0013*IT_0112;
    const ccomplex_t IT_0114 = s_13*IT_0096;
    const ccomplex_t IT_0115 = IT_0095*IT_0114;
    const ccomplex_t IT_0116 = 0.666666666666667*IT_0113 + (-0.166666666666667
      )*IT_0115;
    const ccomplex_t IT_0117 = s_13 + IT_0116;
    const ccomplex_t IT_0118 = IT_0117*m_chi_dm_2;
    const ccomplex_t IT_0119 = (-6)*IT_0118;
    const ccomplex_t IT_0120 = s_23*IT_0096;
    const ccomplex_t IT_0121 = s_24*s_34*IT_0000;
    const ccomplex_t IT_0122 = (-2)*IT_0121;
    const ccomplex_t IT_0123 = IT_0120 + IT_0122;
    const ccomplex_t IT_0124 = IT_0095*IT_0123;
    const ccomplex_t IT_0125 = s_23*IT_0000;
    const ccomplex_t IT_0126 = s_24*s_34;
    const ccomplex_t IT_0127 = (-2)*IT_0126;
    const ccomplex_t IT_0128 = IT_0125 + IT_0127;
    const ccomplex_t IT_0129 = IT_0013*IT_0128;
    const ccomplex_t IT_0130 = 0.166666666666667*IT_0124 + (-0.666666666666667
      )*IT_0129;
    const ccomplex_t IT_0131 = s_23 + IT_0130;
    const ccomplex_t IT_0132 = IT_0131*m_chi_dm_2;
    const ccomplex_t IT_0133 = 6*IT_0132;
    const ccomplex_t IT_0134 = -IT_0010;
    const ccomplex_t IT_0135 = s_24*IT_0000;
    const ccomplex_t IT_0136 = s_34*IT_0095;
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = IT_0013*IT_0125;
    const ccomplex_t IT_0139 = s_14*IT_0000;
    const ccomplex_t IT_0140 = s_13*s_34;
    const ccomplex_t IT_0141 = (-2)*IT_0140;
    const ccomplex_t IT_0142 = IT_0139 + IT_0141;
    const ccomplex_t IT_0143 = IT_0136*IT_0142;
    const ccomplex_t IT_0144 = -s_23 + (-0.5)*IT_0113 + (-0.25)*IT_0137 + 0.5
      *IT_0138 + (-0.25)*IT_0143;
    const ccomplex_t IT_0145 = s_13 + IT_0144;
    const ccomplex_t IT_0146 = IT_0145*m_chi_dm_2;
    const ccomplex_t IT_0147 = (-4)*IT_0146;
    const ccomplex_t IT_0148 = -IT_0052;
    const ccomplex_t IT_0149 = IT_0055*(IT_0052 + -IT_0148);
    const ccomplex_t IT_0150 = -conj(IT_0149);
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = s_14*s_34*IT_0013;
    const ccomplex_t IT_0153 = s_24*s_34*IT_0013;
    const ccomplex_t IT_0154 = -IT_0153;
    const ccomplex_t IT_0155 = s_23 + IT_0154;
    const ccomplex_t IT_0156 = -IT_0126;
    const ccomplex_t IT_0157 = IT_0125 + IT_0156;
    const ccomplex_t IT_0158 = IT_0013*IT_0157;
    const ccomplex_t IT_0159 = s_14*s_34;
    const ccomplex_t IT_0160 = IT_0112 + IT_0159;
    const ccomplex_t IT_0161 = pow(s_34, 2);
    const ccomplex_t IT_0162 = s_13*IT_0013*IT_0161;
    const ccomplex_t IT_0163 = (-2)*IT_0162;
    const ccomplex_t IT_0164 = IT_0160 + IT_0163;
    const ccomplex_t IT_0165 = IT_0013*IT_0164;
    const ccomplex_t IT_0166 = -IT_0152 + 2*IT_0155 + (-0.5)*IT_0158 + 0.5
      *IT_0165;
    const ccomplex_t IT_0167 = s_13 + IT_0166;
    const ccomplex_t IT_0168 = IT_0167*m_chi_dm_2;
    const ccomplex_t IT_0169 = (-2)*IT_0168;
    const ccomplex_t IT_0170 = 2*IT_0052;
    const ccomplex_t IT_0171 = IT_0055*(IT_0052 + -IT_0170);
    const ccomplex_t IT_0172 = IT_0171 + -conj(IT_0171);
    const ccomplex_t IT_0173 = (IT_0034 + IT_0050)*(IT_0057 + -conj(IT_0057)) 
      + -(IT_0069 + IT_0084)*(IT_0087 + -conj(IT_0087)) + 0.25*IT_0012*((-16)
      *s_12 + 32*IT_0088 + 8*IT_0091 + 8*IT_0094 + (-4)*IT_0103) + 2*(IT_0106 + 
      -IT_0110)*IT_0111 + 2*(IT_0119 + IT_0133)*IT_0134 + -IT_0147*IT_0151 +
       IT_0169*IT_0172;
    const ccomplex_t IT_0174 = pow(m_V, 6);
    const ccomplex_t IT_0175 = s_12*IT_0174;
    const ccomplex_t IT_0176 = s_14*s_24*IT_0096;
    const ccomplex_t IT_0177 = (-2)*IT_0176;
    const ccomplex_t IT_0178 = IT_0175 + IT_0177;
    const ccomplex_t IT_0179 = IT_0095*IT_0178;
    const ccomplex_t IT_0180 = IT_0013*IT_0102;
    const ccomplex_t IT_0181 = 9*IT_0014;
    const ccomplex_t IT_0182 = s_12*s_13;
    const ccomplex_t IT_0183 = s_23*IT_0013*IT_0015;
    const ccomplex_t IT_0184 = s_13*s_14*s_24*IT_0013;
    const ccomplex_t IT_0185 = s_23*IT_0015;
    const ccomplex_t IT_0186 = s_13*s_14*s_24;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = s_12*s_14*s_34;
    const ccomplex_t IT_0189 = s_23*IT_0000*IT_0002;
    const ccomplex_t IT_0190 = -IT_0188 + -IT_0189;
    const ccomplex_t IT_0191 = IT_0187 + IT_0190;
    const ccomplex_t IT_0192 = IT_0013*IT_0191;
    const ccomplex_t IT_0193 = -IT_0152;
    const ccomplex_t IT_0194 = s_13 + IT_0193;
    const ccomplex_t IT_0195 = IT_0002*IT_0194;
    const ccomplex_t IT_0196 = IT_0002*IT_0013;
    const ccomplex_t IT_0197 = -IT_0159;
    const ccomplex_t IT_0198 = IT_0112 + IT_0197;
    const ccomplex_t IT_0199 = IT_0196*IT_0198;
    const ccomplex_t IT_0200 = (-0.25)*IT_0199;
    const ccomplex_t IT_0201 = s_12*s_14*s_34*IT_0013;
    const ccomplex_t IT_0202 = (-0.5)*IT_0201;
    const ccomplex_t IT_0203 = s_23*IT_0002;
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = IT_0182 + 0.5*IT_0183 + (-0.5)*IT_0184 + (-0.25
      )*IT_0192 + IT_0195 + IT_0200 + IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = s_14*s_23*s_24*IT_0013;
    const ccomplex_t IT_0207 = s_13*IT_0002;
    const ccomplex_t IT_0208 = s_13*IT_0013*IT_0040;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0002*IT_0155;
    const ccomplex_t IT_0211 = 2*IT_0210;
    const ccomplex_t IT_0212 = IT_0157*IT_0196;
    const ccomplex_t IT_0213 = (-0.5)*IT_0212;
    const ccomplex_t IT_0214 = s_14*s_23*s_24;
    const ccomplex_t IT_0215 = s_13*IT_0040;
    const ccomplex_t IT_0216 = IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = s_12*s_24*s_34;
    const ccomplex_t IT_0218 = s_13*IT_0000*IT_0002;
    const ccomplex_t IT_0219 = -IT_0217 + -IT_0218;
    const ccomplex_t IT_0220 = IT_0216 + IT_0219;
    const ccomplex_t IT_0221 = IT_0013*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = s_12*s_24*s_34*IT_0013;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0206 + IT_0207 + IT_0209 + IT_0211 + IT_0213
       + IT_0222 + IT_0224;
    const ccomplex_t IT_0226 = s_14*s_23*IT_0000;
    const ccomplex_t IT_0227 = s_13*s_24*IT_0000;
    const ccomplex_t IT_0228 = IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = s_12*s_34*IT_0000;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = IT_0136*IT_0231;
    const ccomplex_t IT_0233 = (-2)*IT_0161;
    const ccomplex_t IT_0234 = IT_0096 + IT_0233;
    const ccomplex_t IT_0235 = IT_0196*IT_0234;
    const ccomplex_t IT_0236 = (-0.5)*IT_0235;
    const ccomplex_t IT_0237 = s_34*IT_0000;
    const ccomplex_t IT_0238 = s_34*IT_0002*IT_0095;
    const ccomplex_t IT_0239 = IT_0237*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = (-2)*IT_0098;
    const ccomplex_t IT_0242 = IT_0097 + IT_0241;
    const ccomplex_t IT_0243 = IT_0013*IT_0242;
    const ccomplex_t IT_0244 = 0.5*IT_0243;
    const ccomplex_t IT_0245 = s_14*s_23*s_34;
    const ccomplex_t IT_0246 = s_12*IT_0161;
    const ccomplex_t IT_0247 = (-2)*IT_0098 + 2*IT_0099 + 2*IT_0245 + (-2)
      *IT_0246;
    const ccomplex_t IT_0248 = IT_0097 + IT_0247;
    const ccomplex_t IT_0249 = IT_0013*IT_0248;
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = 1.5*IT_0014;
    const ccomplex_t IT_0252 = IT_0023 + (-2)*IT_0035 + 0.5*IT_0232 + IT_0236 
      + IT_0240 + IT_0244 + IT_0250 + IT_0251;
    const ccomplex_t IT_0253 = s_14*s_23*s_34*IT_0013;
    const ccomplex_t IT_0254 = IT_0013*IT_0161;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = IT_0000 + IT_0255;
    const ccomplex_t IT_0257 = IT_0002*IT_0256;
    const ccomplex_t IT_0258 = 2*IT_0257;
    const ccomplex_t IT_0259 = s_12*IT_0013*IT_0161;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = -IT_0161;
    const ccomplex_t IT_0262 = IT_0096 + IT_0261;
    const ccomplex_t IT_0263 = IT_0196*IT_0262;
    const ccomplex_t IT_0264 = (-0.5)*IT_0263;
    const ccomplex_t IT_0265 = IT_0097 + IT_0099 + IT_0245;
    const ccomplex_t IT_0266 = (-2)*IT_0098 + -IT_0246;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = IT_0013*IT_0267;
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = -IT_0036;
    const ccomplex_t IT_0271 = IT_0023 + IT_0253 + IT_0258 + IT_0260 + IT_0264
       + IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = s_12*s_13*IT_0000;
    const ccomplex_t IT_0273 = s_24*s_34*IT_0002;
    const ccomplex_t IT_0274 = 0.5*IT_0189 + (-0.5)*IT_0272 + -IT_0273;
    const ccomplex_t IT_0275 = IT_0186 + IT_0274;
    const ccomplex_t IT_0276 = IT_0013*IT_0275;
    const ccomplex_t IT_0277 = (-2)*IT_0276;
    const ccomplex_t IT_0278 = (-0.5)*IT_0277;
    const ccomplex_t IT_0279 = -IT_0203;
    const ccomplex_t IT_0280 = IT_0182 + IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = s_12*s_23;
    const ccomplex_t IT_0282 = s_12*s_23*IT_0000;
    const ccomplex_t IT_0283 = -IT_0217 + (-0.5)*IT_0218 + 0.5*IT_0282;
    const ccomplex_t IT_0284 = IT_0215 + IT_0283;
    const ccomplex_t IT_0285 = IT_0013*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0285;
    const ccomplex_t IT_0287 = (-0.5)*IT_0286;
    const ccomplex_t IT_0288 = -IT_0207;
    const ccomplex_t IT_0289 = IT_0281 + IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = IT_0280 + -IT_0289;
    const ccomplex_t IT_0291 = IT_0134*(IT_0023 + (-4)*IT_0035 + 0.5*IT_0179 +
       -IT_0180 + 0.5*IT_0181) + 0.125*conj(IT_0057)*((-16)*IT_0205 + (-8)
      *IT_0225) + -IT_0151*IT_0252 + IT_0172*IT_0271 + conj(IT_0087)*IT_0290;
    const ccomplex_t IT_0292 = (IT_0002 + -IT_0059)*(s_12 + (-0.5)*IT_0094);
    const ccomplex_t IT_0293 = 2*IT_0292;
    const ccomplex_t IT_0294 = -IT_0059;
    const ccomplex_t IT_0295 = IT_0002 + IT_0294;
    const ccomplex_t IT_0296 = IT_0002*IT_0295;
    const ccomplex_t IT_0297 = -IT_0296;
    const ccomplex_t IT_0298 = s_13*s_14*s_23*s_24*IT_0095;
    const ccomplex_t IT_0299 = pow(m_chi_dm_2, 4);
    const ccomplex_t IT_0300 = IT_0018*IT_0040*IT_0095;
    const ccomplex_t IT_0301 = s_13*s_24*s_34*IT_0095;
    const ccomplex_t IT_0302 = -IT_0028 + -IT_0058;
    const ccomplex_t IT_0303 = s_12 + IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = IT_0002*IT_0303;
    const ccomplex_t IT_0305 = -IT_0304;
    const ccomplex_t IT_0306 = IT_0002*IT_0013*IT_0040;
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = IT_0002*IT_0013*IT_0018;
    const ccomplex_t IT_0309 = -IT_0308;
    const ccomplex_t IT_0310 = s_12*s_13*s_24*s_34*IT_0095;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = IT_0298 + IT_0299 + IT_0300 + IT_0305 + IT_0307
       + IT_0309 + IT_0311;
    const ccomplex_t IT_0313 = (IT_0002 + -IT_0059)*(s_12 + 1.5*IT_0002 + (
      -0.5)*IT_0094);
    const ccomplex_t IT_0314 = 2*IT_0313;
    const ccomplex_t IT_0315 = (s_12 + -IT_0058)*(s_12 + 1.5*IT_0002 + (-0.5)
      *IT_0094);
    const ccomplex_t IT_0316 = (-2)*IT_0315;
    const ccomplex_t IT_0317 = IT_0314 + -IT_0316;
    const ccomplex_t IT_0318 = s_24*s_34*IT_0018*IT_0095;
    const ccomplex_t IT_0319 = s_13*IT_0095*IT_0161;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = IT_0152 + IT_0320;
    const ccomplex_t IT_0322 = IT_0002*IT_0321;
    const ccomplex_t IT_0323 = s_13*s_14*s_23*s_34*IT_0095;
    const ccomplex_t IT_0324 = IT_0182 + IT_0203 + IT_0323;
    const ccomplex_t IT_0325 = -IT_0184;
    const ccomplex_t IT_0326 = -IT_0183;
    const ccomplex_t IT_0327 = s_23*IT_0013*IT_0018;
    const ccomplex_t IT_0328 = (-2)*IT_0327;
    const ccomplex_t IT_0329 = s_12*s_13*IT_0095*IT_0161;
    const ccomplex_t IT_0330 = -IT_0329;
    const ccomplex_t IT_0331 = IT_0318 + IT_0322 + IT_0324;
    const ccomplex_t IT_0332 = s_24*s_34*IT_0002*IT_0013;
    const ccomplex_t IT_0333 = -IT_0332;
    const ccomplex_t IT_0334 = s_13*s_23*IT_0015*IT_0095;
    const ccomplex_t IT_0335 = s_14*s_24*IT_0018*IT_0095;
    const ccomplex_t IT_0336 = s_12*IT_0002;
    const ccomplex_t IT_0337 = s_13*s_14*s_34*IT_0095;
    const ccomplex_t IT_0338 = IT_0002 + IT_0337;
    const ccomplex_t IT_0339 = -IT_0022 + -IT_0059;
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = IT_0002*IT_0340;
    const ccomplex_t IT_0342 = -IT_0341;
    const ccomplex_t IT_0343 = s_14*s_24*IT_0002*IT_0013;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = s_13*s_23*IT_0002*IT_0013;
    const ccomplex_t IT_0346 = -IT_0345;
    const ccomplex_t IT_0347 = s_12*s_13*s_14*s_34*IT_0095;
    const ccomplex_t IT_0348 = -IT_0347;
    const ccomplex_t IT_0349 = IT_0334 + IT_0335 + IT_0336 + IT_0342 + IT_0344
       + IT_0346 + IT_0348;
    const ccomplex_t IT_0350 = 4*IT_0349;
    const ccomplex_t IT_0351 = conj(IT_0057)*IT_0350;
    const ccomplex_t IT_0352 = (-4)*IT_0134*IT_0280 + IT_0111*(IT_0293 + (-3)
      *IT_0297) + conj(IT_0057)*IT_0312 + -conj(IT_0087)*IT_0317 + -IT_0151*
      (IT_0201 + IT_0318 + IT_0322 + IT_0324 + IT_0325 + IT_0326 + IT_0328 +
       IT_0330) + -IT_0172*(IT_0328 + IT_0330 + IT_0331 + IT_0333) + 0.25*IT_0351;
    const ccomplex_t IT_0353 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0091)*(IT_0000 
      + -IT_0254);
    const ccomplex_t IT_0354 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0091)*(s_13 + 
      -IT_0152);
    const ccomplex_t IT_0355 = (-8)*IT_0354;
    const ccomplex_t IT_0356 = (s_12 + 1.5*IT_0002 + (-0.5)*IT_0091)*(s_23 + 
      -IT_0153);
    const ccomplex_t IT_0357 = (-8)*IT_0356;
    const ccomplex_t IT_0358 = IT_0355 + IT_0357;
    const ccomplex_t IT_0359 = s_13*s_23*s_24*s_34*IT_0095;
    const ccomplex_t IT_0360 = s_14*s_34*IT_0070*IT_0095;
    const ccomplex_t IT_0361 = s_23*IT_0095*IT_0161;
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = IT_0153 + IT_0362;
    const ccomplex_t IT_0364 = IT_0002*IT_0363;
    const ccomplex_t IT_0365 = s_13*IT_0013*IT_0070;
    const ccomplex_t IT_0366 = (-2)*IT_0365;
    const ccomplex_t IT_0367 = IT_0281 + IT_0366;
    const ccomplex_t IT_0368 = s_14*s_34*IT_0002*IT_0013;
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = s_12*s_23*IT_0095*IT_0161;
    const ccomplex_t IT_0371 = -IT_0370;
    const ccomplex_t IT_0372 = IT_0207 + IT_0359 + IT_0360 + IT_0364 + IT_0367
       + IT_0369 + IT_0371;
    const ccomplex_t IT_0373 = s_13*s_23*IT_0095*IT_0161;
    const ccomplex_t IT_0374 = (-2)*IT_0373;
    const ccomplex_t IT_0375 = IT_0036 + IT_0253 + IT_0374;
    const ccomplex_t IT_0376 = (-4)*IT_0375;
    const ccomplex_t IT_0377 = IT_0095*IT_0161;
    const ccomplex_t IT_0378 = (IT_0014 + IT_0023 + (-2)*IT_0035)*(1 + 0.5
      *IT_0377);
    const ccomplex_t IT_0379 = -IT_0206;
    const ccomplex_t IT_0380 = IT_0209 + IT_0367 + IT_0379;
    const ccomplex_t IT_0381 = 2*IT_0359;
    const ccomplex_t IT_0382 = IT_0207 + IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = 4*IT_0382;
    const ccomplex_t IT_0384 = 2*IT_0323;
    const ccomplex_t IT_0385 = IT_0182 + IT_0203 + IT_0325 + IT_0326 + IT_0328
       + 0.25*IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = (-0.125)*conj(IT_0171)*IT_0376 + conj(IT_0149)
      *IT_0378 + (-0.5)*conj(IT_0087)*(IT_0201 + IT_0207 + IT_0223 + IT_0318 +
       IT_0322 + IT_0324 + IT_0325 + IT_0326 + IT_0328 + IT_0330 + IT_0359 +
       IT_0360 + IT_0364 + IT_0371 + IT_0380) + (-0.5)*conj(IT_0057)*IT_0385;
    const ccomplex_t IT_0387 = (IT_0002 + -IT_0022)*(s_12 + 1.5*IT_0002 + (
      -0.5)*IT_0091);
    const ccomplex_t IT_0388 = -IT_0387;
    const ccomplex_t IT_0389 = (s_12 + -IT_0028)*(s_12 + 1.5*IT_0002 + (-0.5)
      *IT_0091);
    const ccomplex_t IT_0390 = (IT_0002 + -IT_0044)*(s_12 + 1.5*IT_0002 + (
      -0.5)*IT_0091);
    const ccomplex_t IT_0391 = -IT_0390;
    const ccomplex_t IT_0392 = IT_0015*IT_0070*IT_0095;
    const ccomplex_t IT_0393 = s_14*s_23*s_34*IT_0095;
    const ccomplex_t IT_0394 = s_12 + IT_0302 + IT_0393;
    const ccomplex_t IT_0395 = IT_0002*IT_0394;
    const ccomplex_t IT_0396 = -IT_0395;
    const ccomplex_t IT_0397 = IT_0002*IT_0013*IT_0070;
    const ccomplex_t IT_0398 = -IT_0397;
    const ccomplex_t IT_0399 = IT_0002*IT_0013*IT_0015;
    const ccomplex_t IT_0400 = -IT_0399;
    const ccomplex_t IT_0401 = s_12*s_14*s_23*s_34*IT_0095;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = IT_0298 + IT_0299 + IT_0392 + IT_0396 + IT_0398
       + IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = IT_0312 + IT_0403;
    const ccomplex_t IT_0405 = s_14*s_24*IT_0070*IT_0095;
    const ccomplex_t IT_0406 = s_13*s_23*IT_0040*IT_0095;
    const ccomplex_t IT_0407 = s_23*s_24*s_34*IT_0095;
    const ccomplex_t IT_0408 = IT_0002 + IT_0407;
    const ccomplex_t IT_0409 = -IT_0044 + -IT_0071;
    const ccomplex_t IT_0410 = IT_0408 + IT_0409;
    const ccomplex_t IT_0411 = IT_0002*IT_0410;
    const ccomplex_t IT_0412 = -IT_0411;
    const ccomplex_t IT_0413 = s_12*s_23*s_24*s_34*IT_0095;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0336 + IT_0344 + IT_0346 + IT_0405 + IT_0406
       + IT_0412 + IT_0414;
    const ccomplex_t IT_0416 = 4*IT_0415;
    const ccomplex_t IT_0417 = 0.25*IT_0416;
    const ccomplex_t IT_0418 = 0.5*conj(IT_0149);
    const ccomplex_t IT_0419 = 2*IT_0134*(IT_0205 + 0.5*IT_0225) + (-0.5)
      *IT_0111*(IT_0312 + IT_0349) + 0.125*conj(IT_0171)*IT_0358 + conj(IT_0057)
      *(IT_0388 + (-2)*IT_0389 + IT_0391) + 0.5*conj(IT_0087)*(0.25*IT_0350 +
       IT_0404 + IT_0417) + IT_0385*IT_0418;
    const ccomplex_t IT_0420 = -IT_0313;
    const ccomplex_t IT_0421 = (IT_0002 + -IT_0071)*(s_12 + 1.5*IT_0002 + (
      -0.5)*IT_0094);
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = -IT_0134*IT_0290 + 0.5*IT_0111*IT_0317 + 0.125
      *IT_0351 + 0.5*conj(IT_0171)*(IT_0328 + IT_0330 + IT_0331 + IT_0333 +
       IT_0372) + 0.5*conj(IT_0057)*(IT_0404 + IT_0417) + (IT_0201 + IT_0207 +
       IT_0223 + IT_0318 + IT_0322 + IT_0324 + IT_0325 + IT_0326 + IT_0328 +
       IT_0330 + IT_0359 + IT_0360 + IT_0364 + IT_0371 + IT_0380)*IT_0418 + (-2)
      *conj(IT_0087)*(IT_0315 + -1./2*IT_0420 + -1./2*IT_0422);
    const ccomplex_t IT_0424 = 4*IT_0012*IT_0173 + 8*IT_0134*IT_0291 + 4
      *IT_0111*IT_0352 + (-8)*IT_0171*(conj(IT_0171)*IT_0353 + (-0.125)*conj
      (IT_0057)*IT_0358 + (-0.5)*conj(IT_0087)*(IT_0328 + IT_0330 + IT_0331 +
       IT_0333 + IT_0372) + (-0.125)*conj(IT_0149)*IT_0376) + (-8)*IT_0149
      *IT_0386 + 8*IT_0057*IT_0419 + 8*IT_0087*IT_0423;
    return create_ccomplex_return(IT_0424);
}

