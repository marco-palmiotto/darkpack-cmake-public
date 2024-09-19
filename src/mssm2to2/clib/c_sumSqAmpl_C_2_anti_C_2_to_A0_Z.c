#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_A0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = e_em*IT_0001*(IT_0002*IT_0003 + -IT_0004
      *IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = IT_0007*(IT_0010 + IT_0012);
    const ccomplex_t IT_0014 = 1.4142135623731*e_em*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = 0.5*IT_0015;
    const ccomplex_t IT_0017 = pow(m_C_2, 2);
    const ccomplex_t IT_0018 = cpow((-2)*s_12 + (-2)*IT_0017 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = IT_0008*IT_0020;
    const ccomplex_t IT_0022 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = IT_0007*(IT_0023 + IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0001*IT_0006;
    const ccomplex_t IT_0030 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = IT_0000*IT_0007;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0031 + -IT_0033 + (-2)
      *IT_0035);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0028*IT_0037;
    const ccomplex_t IT_0039 = pow(m_A0, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0017 + IT_0039 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = e_em*IT_0001*(IT_0002*IT_0004 + IT_0003*IT_0005
      )*IT_0007;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0002*IT_0011;
    const ccomplex_t IT_0046 = IT_0005*IT_0009;
    const ccomplex_t IT_0047 = IT_0007*(IT_0045 + -IT_0046);
    const ccomplex_t IT_0048 = 1.4142135623731*e_em*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0017 + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0044*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0004*IT_0011;
    const ccomplex_t IT_0057 = IT_0003*IT_0009;
    const ccomplex_t IT_0058 = IT_0007*(IT_0056 + IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0062 = IT_0029*IT_0061;
    const ccomplex_t IT_0063 = IT_0032*IT_0061;
    const ccomplex_t IT_0064 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0065 = IT_0032*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + -IT_0063 + (-2)
      *IT_0065);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0060*IT_0067;
    const ccomplex_t IT_0069 = IT_0055*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0021 + 2*IT_0042 + IT_0054 + 2*IT_0070;
    const ccomplex_t IT_0072 = -IT_0042 + -IT_0070;
    const ccomplex_t IT_0073 = s_12*s_24;
    const ccomplex_t IT_0074 = s_14*IT_0017;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = IT_0073 + IT_0075;
    const ccomplex_t IT_0077 = 2*IT_0076;
    const ccomplex_t IT_0078 = pow(m_Z, -2);
    const ccomplex_t IT_0079 = pow(s_24, 2);
    const ccomplex_t IT_0080 = IT_0078*IT_0079;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = IT_0017 + IT_0081;
    const ccomplex_t IT_0083 = s_12*IT_0082;
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = 2*IT_0084;
    const ccomplex_t IT_0086 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0087 = IT_0029*IT_0086;
    const ccomplex_t IT_0088 = IT_0032*IT_0086;
    const ccomplex_t IT_0089 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0090 = IT_0032*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + -IT_0088 + (-2)
      *IT_0090);
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0060*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0095 = IT_0055*IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0097 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0098 = IT_0029*IT_0097;
    const ccomplex_t IT_0099 = IT_0032*IT_0097;
    const ccomplex_t IT_0100 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0101 = IT_0032*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + -IT_0099 + (-2)
      *IT_0101);
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0028*IT_0103;
    const ccomplex_t IT_0105 = IT_0040*IT_0096*IT_0104;
    const ccomplex_t IT_0106 = V_u1*U_Wm2;
    const ccomplex_t IT_0107 = IT_0004*IT_0106;
    const ccomplex_t IT_0108 = U_d2*V_Wp1;
    const ccomplex_t IT_0109 = IT_0003*IT_0108;
    const ccomplex_t IT_0110 = IT_0007*(IT_0107 + IT_0109);
    const ccomplex_t IT_0111 = 1.4142135623731*e_em*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0114 = IT_0029*IT_0113;
    const ccomplex_t IT_0115 = IT_0032*IT_0113;
    const ccomplex_t IT_0116 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0117 = IT_0032*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + -IT_0115 + (-2)
      *IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0112*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_23 + IT_0017 + IT_0039 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_C_2*IT_0123;
    const ccomplex_t IT_0125 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0126 = IT_0029*IT_0125;
    const ccomplex_t IT_0127 = IT_0032*IT_0125;
    const ccomplex_t IT_0128 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0129 = IT_0032*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + -IT_0127 + (-2)
      *IT_0129);
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0112*IT_0131;
    const ccomplex_t IT_0133 = IT_0096*IT_0121*IT_0132;
    const ccomplex_t IT_0134 = cpow((-2)*s_23 + IT_0039 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0135 = V_u2*U_Wm2;
    const ccomplex_t IT_0136 = IT_0004*IT_0135;
    const ccomplex_t IT_0137 = U_d2*V_Wp2;
    const ccomplex_t IT_0138 = IT_0003*IT_0137;
    const ccomplex_t IT_0139 = IT_0007*(IT_0136 + IT_0138);
    const ccomplex_t IT_0140 = 1.4142135623731*e_em*IT_0139;
    const ccomplex_t IT_0141 = 0.5*IT_0140;
    const ccomplex_t IT_0142 = IT_0092*IT_0141;
    const ccomplex_t IT_0143 = IT_0094*IT_0134*IT_0142;
    const ccomplex_t IT_0144 = m_C_2*IT_0070;
    const ccomplex_t IT_0145 = m_C_2*IT_0042;
    const ccomplex_t IT_0146 = IT_0067*IT_0141;
    const ccomplex_t IT_0147 = IT_0134*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = m_C_2*IT_0148;
    const ccomplex_t IT_0150 = IT_0095 + IT_0105 + IT_0124 + IT_0133 + IT_0143
       + IT_0144 + IT_0145 + IT_0149;
    const ccomplex_t IT_0151 = s_14*s_24*IT_0078;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = s_12 + IT_0152;
    const ccomplex_t IT_0154 = m_C_2*IT_0153;
    const ccomplex_t IT_0155 = 2*IT_0154;
    const ccomplex_t IT_0156 = IT_0002*IT_0137;
    const ccomplex_t IT_0157 = IT_0005*IT_0135;
    const ccomplex_t IT_0158 = IT_0007*(IT_0156 + IT_0157);
    const ccomplex_t IT_0159 = 1.4142135623731*e_em*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = IT_0018*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = IT_0008*IT_0163;
    const ccomplex_t IT_0165 = IT_0002*IT_0135;
    const ccomplex_t IT_0166 = IT_0005*IT_0137;
    const ccomplex_t IT_0167 = IT_0007*(IT_0165 + -IT_0166);
    const ccomplex_t IT_0168 = 1.4142135623731*e_em*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = IT_0051*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = IT_0044*IT_0172;
    const ccomplex_t IT_0174 = (-2)*IT_0123 + (-2)*IT_0148 + IT_0164 + IT_0173;
    const ccomplex_t IT_0175 = IT_0017*IT_0153;
    const ccomplex_t IT_0176 = 2*IT_0175;
    const ccomplex_t IT_0177 = IT_0055*IT_0094*IT_0146;
    const ccomplex_t IT_0178 = V_u2*U_Wm1;
    const ccomplex_t IT_0179 = IT_0004*IT_0178;
    const ccomplex_t IT_0180 = U_d1*V_Wp2;
    const ccomplex_t IT_0181 = IT_0003*IT_0180;
    const ccomplex_t IT_0182 = IT_0007*(IT_0179 + IT_0181);
    const ccomplex_t IT_0183 = 1.4142135623731*e_em*IT_0182;
    const ccomplex_t IT_0184 = 0.5*IT_0183;
    const ccomplex_t IT_0185 = IT_0037*IT_0184;
    const ccomplex_t IT_0186 = IT_0040*IT_0096*IT_0185;
    const ccomplex_t IT_0187 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0188 = IT_0004*IT_0187;
    const ccomplex_t IT_0189 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0190 = IT_0003*IT_0189;
    const ccomplex_t IT_0191 = IT_0007*(IT_0188 + IT_0190);
    const ccomplex_t IT_0192 = 1.4142135623731*e_em*IT_0191;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0119*IT_0193;
    const ccomplex_t IT_0195 = IT_0096*IT_0121*IT_0194;
    const ccomplex_t IT_0196 = IT_0131*IT_0193;
    const ccomplex_t IT_0197 = IT_0121*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = m_C_2*IT_0198;
    const ccomplex_t IT_0200 = IT_0103*IT_0184;
    const ccomplex_t IT_0201 = IT_0040*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = m_C_2*IT_0202;
    const ccomplex_t IT_0204 = IT_0068*IT_0094*IT_0134;
    const ccomplex_t IT_0205 = IT_0093*IT_0134;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = m_C_2*IT_0206;
    const ccomplex_t IT_0208 = IT_0055*IT_0142;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = m_C_2*IT_0209;
    const ccomplex_t IT_0211 = -IT_0177 + -IT_0186 + -IT_0195 + -IT_0199 + 
      -IT_0203 + -IT_0204 + -IT_0207 + -IT_0210;
    const ccomplex_t IT_0212 = m_C_2*IT_0082;
    const ccomplex_t IT_0213 = -IT_0212;
    const ccomplex_t IT_0214 = 2*IT_0213;
    const ccomplex_t IT_0215 = -IT_0198 + -IT_0206;
    const ccomplex_t IT_0216 = (-2)*IT_0076;
    const ccomplex_t IT_0217 = IT_0021 + IT_0054 + 2*IT_0198 + 2*IT_0206;
    const ccomplex_t IT_0218 = s_12*IT_0153;
    const ccomplex_t IT_0219 = (-2)*IT_0218;
    const ccomplex_t IT_0220 = IT_0164 + IT_0173 + (-2)*IT_0202 + (-2)*IT_0209;
    const ccomplex_t IT_0221 = IT_0017*IT_0082;
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = (-2)*IT_0222;
    const ccomplex_t IT_0224 = IT_0202 + IT_0209;
    const ccomplex_t IT_0225 = IT_0123 + IT_0148;
    const ccomplex_t IT_0226 = pow(m_Z, 2);
    const ccomplex_t IT_0227 = s_24*IT_0226;
    const ccomplex_t IT_0228 = IT_0078*IT_0227;
    const ccomplex_t IT_0229 = (-0.25)*IT_0228;
    const ccomplex_t IT_0230 = s_24 + IT_0229;
    const ccomplex_t IT_0231 = m_C_2*IT_0230;
    const ccomplex_t IT_0232 = (-4)*IT_0231;
    const ccomplex_t IT_0233 = 2*IT_0232;
    const ccomplex_t IT_0234 = 0.5*IT_0228;
    const ccomplex_t IT_0235 = s_24 + IT_0234;
    const ccomplex_t IT_0236 = m_C_2*IT_0235;
    const ccomplex_t IT_0237 = 2*IT_0236;
    const ccomplex_t IT_0238 = 2*IT_0237;
    const ccomplex_t IT_0239 = s_12*IT_0226;
    const ccomplex_t IT_0240 = s_14*s_24;
    const ccomplex_t IT_0241 = (-2)*IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = IT_0078*IT_0242;
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = s_12 + IT_0244;
    const ccomplex_t IT_0246 = 4*IT_0245;
    const ccomplex_t IT_0247 = 6*IT_0017;
    const ccomplex_t IT_0248 = pow(s_14, 2);
    const ccomplex_t IT_0249 = IT_0078*IT_0248;
    const ccomplex_t IT_0250 = -IT_0017;
    const ccomplex_t IT_0251 = IT_0249 + IT_0250;
    const ccomplex_t IT_0252 = m_C_2*IT_0251;
    const ccomplex_t IT_0253 = (-2)*IT_0252;
    const ccomplex_t IT_0254 = (-2)*IT_0154;
    const ccomplex_t IT_0255 = s_14*IT_0226;
    const ccomplex_t IT_0256 = IT_0078*IT_0255;
    const ccomplex_t IT_0257 = 0.5*IT_0256;
    const ccomplex_t IT_0258 = s_14 + IT_0257;
    const ccomplex_t IT_0259 = m_C_2*IT_0258;
    const ccomplex_t IT_0260 = 2*IT_0259;
    const ccomplex_t IT_0261 = (-2)*IT_0260;
    const ccomplex_t IT_0262 = (-0.25)*IT_0256;
    const ccomplex_t IT_0263 = s_14 + IT_0262;
    const ccomplex_t IT_0264 = m_C_2*IT_0263;
    const ccomplex_t IT_0265 = (-4)*IT_0264;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = IT_0211*(conj(IT_0071)*IT_0214 + IT_0155*conj
      (IT_0220) + conj(IT_0215)*IT_0233 + conj(IT_0072)*IT_0238 + conj(IT_0211)
      *IT_0246 + conj(IT_0150)*IT_0247 + conj(IT_0174)*IT_0253 + conj(IT_0217)
      *IT_0254 + conj(IT_0225)*IT_0261 + conj(IT_0224)*IT_0266);
    const ccomplex_t IT_0268 = IT_0150*(conj(IT_0071)*IT_0155 + IT_0214*conj
      (IT_0220) + conj(IT_0225)*IT_0233 + conj(IT_0224)*IT_0238 + conj(IT_0150)
      *IT_0246 + conj(IT_0211)*IT_0247 + conj(IT_0217)*IT_0253 + conj(IT_0174)
      *IT_0254 + conj(IT_0215)*IT_0261 + conj(IT_0072)*IT_0266);
    const ccomplex_t IT_0269 = s_12*s_14;
    const ccomplex_t IT_0270 = s_24*IT_0017;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = 2*IT_0272;
    const ccomplex_t IT_0274 = pow(m_Z, 4);
    const ccomplex_t IT_0275 = s_12*IT_0274;
    const ccomplex_t IT_0276 = IT_0078*IT_0275;
    const ccomplex_t IT_0277 = (-0.25)*IT_0276;
    const ccomplex_t IT_0278 = IT_0240 + IT_0277;
    const ccomplex_t IT_0279 = (-8)*IT_0278;
    const ccomplex_t IT_0280 = (-0.5)*IT_0276;
    const ccomplex_t IT_0281 = IT_0239 + IT_0241 + IT_0280;
    const ccomplex_t IT_0282 = (-4)*IT_0281;
    const ccomplex_t IT_0283 = IT_0017*IT_0226;
    const ccomplex_t IT_0284 = 3*IT_0283;
    const ccomplex_t IT_0285 = (-2)*IT_0284;
    const ccomplex_t IT_0286 = (-3)*IT_0283;
    const ccomplex_t IT_0287 = (-2)*IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0272;
    const ccomplex_t IT_0289 = s_12*IT_0251;
    const ccomplex_t IT_0290 = IT_0017*IT_0251;
    const ccomplex_t IT_0291 = (-2)*IT_0290;
    const ccomplex_t IT_0292 = 0.5*conj(IT_0220);
    const ccomplex_t IT_0293 = 0.5*conj(IT_0211);
    const ccomplex_t IT_0294 = 0.5*conj(IT_0071);
    const ccomplex_t IT_0295 = 0.5*conj(IT_0150);
    const ccomplex_t IT_0296 = 2*IT_0217*(0.5*conj(IT_0215)*IT_0273 + 0.5*conj
      (IT_0072)*IT_0288 + conj(IT_0217)*IT_0289 + 0.5*conj(IT_0174)*IT_0291 +
       IT_0176*IT_0292 + IT_0254*IT_0293 + IT_0219*IT_0294 + IT_0253*IT_0295);
    const ccomplex_t IT_0297 = 2*IT_0174*(0.5*conj(IT_0225)*IT_0273 + 0.5*conj
      (IT_0224)*IT_0288 + conj(IT_0174)*IT_0289 + 0.5*conj(IT_0217)*IT_0291 +
       IT_0219*IT_0292 + IT_0253*IT_0293 + IT_0176*IT_0294 + IT_0254*IT_0295);
    const ccomplex_t IT_0298 = IT_0071*(conj(IT_0072)*IT_0077 + conj(IT_0071)
      *IT_0085 + conj(IT_0150)*IT_0155 + conj(IT_0174)*IT_0176 + conj(IT_0211)
      *IT_0214 + conj(IT_0215)*IT_0216 + conj(IT_0217)*IT_0219 + conj(IT_0220)
      *IT_0223) + IT_0220*(IT_0155*conj(IT_0211) + conj(IT_0150)*IT_0214 +
       IT_0176*conj(IT_0217) + conj(IT_0174)*IT_0219 + IT_0085*conj(IT_0220) +
       conj(IT_0071)*IT_0223 + IT_0077*conj(IT_0224) + IT_0216*conj(IT_0225)) +
       IT_0267 + IT_0268 + IT_0225*(IT_0216*conj(IT_0220) + conj(IT_0150)
      *IT_0233 + conj(IT_0211)*IT_0261 + conj(IT_0174)*IT_0273 + conj(IT_0224)
      *IT_0279 + conj(IT_0225)*IT_0282 + conj(IT_0072)*IT_0285 + conj(IT_0215)
      *IT_0287) + IT_0215*(conj(IT_0071)*IT_0216 + conj(IT_0211)*IT_0233 + conj
      (IT_0150)*IT_0261 + conj(IT_0217)*IT_0273 + conj(IT_0072)*IT_0279 + conj
      (IT_0215)*IT_0282 + conj(IT_0224)*IT_0285 + conj(IT_0225)*IT_0287) +
       IT_0224*(IT_0077*conj(IT_0220) + conj(IT_0150)*IT_0238 + conj(IT_0211)
      *IT_0266 + conj(IT_0225)*IT_0279 + conj(IT_0224)*IT_0282 + conj(IT_0215)
      *IT_0285 + conj(IT_0072)*IT_0287 + conj(IT_0174)*IT_0288) + IT_0072*(conj
      (IT_0071)*IT_0077 + conj(IT_0211)*IT_0238 + conj(IT_0150)*IT_0266 + conj
      (IT_0215)*IT_0279 + conj(IT_0072)*IT_0282 + conj(IT_0225)*IT_0285 + conj
      (IT_0224)*IT_0287 + conj(IT_0217)*IT_0288) + IT_0296 + IT_0297;
    return create_ccomplex_return(IT_0298);
}

