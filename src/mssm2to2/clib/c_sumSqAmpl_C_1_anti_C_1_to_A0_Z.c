#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_A0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_A0_Z(
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
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = V_u2*U_Wm1;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(beta);
    const ccomplex_t IT_0004 = U_d1*V_Wp2;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = 0.5*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0006*IT_0012;
    const ccomplex_t IT_0014 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0011;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0015 + -IT_0017 + (-2)
      *IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = pow(m_A0, 2);
    const ccomplex_t IT_0024 = pow(m_C_1, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cos(alpha);
    const ccomplex_t IT_0029 = sin(alpha);
    const ccomplex_t IT_0030 = e_em*IT_0007*IT_0012*(IT_0003*IT_0028 + 
      -IT_0000*IT_0029);
    const ccomplex_t IT_0031 = U_d1*V_Wp1;
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = V_u1*U_Wm1;
    const ccomplex_t IT_0034 = IT_0029*IT_0033;
    const ccomplex_t IT_0035 = IT_0007*(IT_0032 + IT_0034);
    const ccomplex_t IT_0036 = 1.4142135623731*e_em*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = cpow((-2)*s_12 + (-2)*IT_0024 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0040 = IT_0038*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = IT_0030*IT_0041;
    const ccomplex_t IT_0043 = e_em*IT_0007*IT_0012*(IT_0000*IT_0028 + IT_0003
      *IT_0029);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0028*IT_0033;
    const ccomplex_t IT_0046 = IT_0029*IT_0031;
    const ccomplex_t IT_0047 = IT_0007*(IT_0045 + -IT_0046);
    const ccomplex_t IT_0048 = 1.4142135623731*e_em*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0044*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0023 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0000*IT_0033;
    const ccomplex_t IT_0057 = IT_0003*IT_0031;
    const ccomplex_t IT_0058 = IT_0007*(IT_0056 + IT_0057);
    const ccomplex_t IT_0059 = 1.4142135623731*e_em*IT_0058;
    const ccomplex_t IT_0060 = 0.5*IT_0059;
    const ccomplex_t IT_0061 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0062 = IT_0013*IT_0061;
    const ccomplex_t IT_0063 = IT_0016*IT_0061;
    const ccomplex_t IT_0064 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0065 = IT_0016*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*(IT_0062 + -IT_0063 + (-2)
      *IT_0065);
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = IT_0060*IT_0067;
    const ccomplex_t IT_0069 = IT_0055*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (-2)*IT_0027 + IT_0042 + IT_0054 + (-2)*IT_0070;
    const ccomplex_t IT_0072 = pow(m_Z, -2);
    const ccomplex_t IT_0073 = pow(s_14, 2);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = -IT_0024;
    const ccomplex_t IT_0076 = IT_0074 + IT_0075;
    const ccomplex_t IT_0077 = s_12*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0079 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0080 = IT_0000*IT_0079;
    const ccomplex_t IT_0081 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = IT_0007*(IT_0080 + IT_0082);
    const ccomplex_t IT_0084 = 1.4142135623731*e_em*IT_0083;
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0067*IT_0085;
    const ccomplex_t IT_0087 = IT_0055*IT_0078*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0023 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0068*IT_0078*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0091 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0092 = IT_0000*IT_0091;
    const ccomplex_t IT_0093 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0094 = IT_0003*IT_0093;
    const ccomplex_t IT_0095 = IT_0007*(IT_0092 + IT_0094);
    const ccomplex_t IT_0096 = 1.4142135623731*e_em*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0021*IT_0097;
    const ccomplex_t IT_0099 = IT_0025*IT_0090*IT_0098;
    const ccomplex_t IT_0100 = V_u1*U_Wm2;
    const ccomplex_t IT_0101 = IT_0000*IT_0100;
    const ccomplex_t IT_0102 = U_d2*V_Wp1;
    const ccomplex_t IT_0103 = IT_0003*IT_0102;
    const ccomplex_t IT_0104 = IT_0007*(IT_0101 + IT_0103);
    const ccomplex_t IT_0105 = 1.4142135623731*e_em*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0105;
    const ccomplex_t IT_0107 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0108 = IT_0013*IT_0107;
    const ccomplex_t IT_0109 = IT_0016*IT_0107;
    const ccomplex_t IT_0110 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0111 = IT_0016*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + -IT_0109 + (-2)
      *IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0090*IT_0114*IT_0115;
    const ccomplex_t IT_0117 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0118 = IT_0013*IT_0117;
    const ccomplex_t IT_0119 = IT_0016*IT_0117;
    const ccomplex_t IT_0120 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0121 = IT_0016*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + -IT_0119 + (-2)
      *IT_0121);
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0085*IT_0123;
    const ccomplex_t IT_0125 = IT_0055*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = m_C_1*IT_0126;
    const ccomplex_t IT_0128 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0129 = IT_0013*IT_0128;
    const ccomplex_t IT_0130 = IT_0016*IT_0128;
    const ccomplex_t IT_0131 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0132 = IT_0016*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + -IT_0130 + (-2)
      *IT_0132);
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = IT_0097*IT_0134;
    const ccomplex_t IT_0136 = IT_0025*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = m_C_1*IT_0137;
    const ccomplex_t IT_0139 = IT_0060*IT_0123;
    const ccomplex_t IT_0140 = IT_0088*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = m_C_1*IT_0141;
    const ccomplex_t IT_0143 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0144 = IT_0013*IT_0143;
    const ccomplex_t IT_0145 = IT_0016*IT_0143;
    const ccomplex_t IT_0146 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0147 = IT_0016*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0144 + -IT_0145 + (-2)
      *IT_0147);
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = IT_0106*IT_0149;
    const ccomplex_t IT_0151 = IT_0115*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = m_C_1*IT_0152;
    const ccomplex_t IT_0154 = -IT_0087 + -IT_0089 + -IT_0099 + -IT_0116 + 
      -IT_0127 + -IT_0138 + -IT_0142 + -IT_0153;
    const ccomplex_t IT_0155 = m_C_1*IT_0076;
    const ccomplex_t IT_0156 = (-2)*IT_0071;
    const ccomplex_t IT_0157 = IT_0042 + IT_0054 + (-2)*IT_0141 + (-2)*IT_0152;
    const ccomplex_t IT_0158 = pow(s_24, 2);
    const ccomplex_t IT_0159 = IT_0072*IT_0158;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0024 + IT_0160;
    const ccomplex_t IT_0162 = s_12*IT_0161;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = 2*IT_0157;
    const ccomplex_t IT_0165 = s_14*s_24*IT_0072;
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = s_12 + IT_0166;
    const ccomplex_t IT_0168 = s_12*IT_0167;
    const ccomplex_t IT_0169 = IT_0028*IT_0079;
    const ccomplex_t IT_0170 = IT_0029*IT_0081;
    const ccomplex_t IT_0171 = IT_0007*(IT_0169 + -IT_0170);
    const ccomplex_t IT_0172 = 1.4142135623731*e_em*IT_0171;
    const ccomplex_t IT_0173 = 0.5*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = IT_0051*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = IT_0044*IT_0176;
    const ccomplex_t IT_0178 = IT_0028*IT_0081;
    const ccomplex_t IT_0179 = IT_0029*IT_0079;
    const ccomplex_t IT_0180 = IT_0007*(IT_0178 + IT_0179);
    const ccomplex_t IT_0181 = 1.4142135623731*e_em*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = IT_0039*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0030*IT_0185;
    const ccomplex_t IT_0187 = 2*IT_0126 + 2*IT_0137 + IT_0177 + IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0187;
    const ccomplex_t IT_0189 = IT_0024*IT_0167;
    const ccomplex_t IT_0190 = IT_0024*IT_0076;
    const ccomplex_t IT_0191 = IT_0141 + IT_0152;
    const ccomplex_t IT_0192 = pow(m_Z, 2);
    const ccomplex_t IT_0193 = s_24*IT_0192;
    const ccomplex_t IT_0194 = IT_0072*IT_0193;
    const ccomplex_t IT_0195 = 0.5*IT_0194;
    const ccomplex_t IT_0196 = s_24 + IT_0195;
    const ccomplex_t IT_0197 = m_C_1*IT_0196;
    const ccomplex_t IT_0198 = 2*IT_0197;
    const ccomplex_t IT_0199 = IT_0055*IT_0078*IT_0139;
    const ccomplex_t IT_0200 = IT_0078*IT_0088*IT_0124;
    const ccomplex_t IT_0201 = IT_0010*IT_0134;
    const ccomplex_t IT_0202 = IT_0025*IT_0090*IT_0201;
    const ccomplex_t IT_0203 = m_C_1*IT_0027;
    const ccomplex_t IT_0204 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0205 = IT_0000*IT_0204;
    const ccomplex_t IT_0206 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0207 = IT_0003*IT_0206;
    const ccomplex_t IT_0208 = IT_0007*(IT_0205 + IT_0207);
    const ccomplex_t IT_0209 = 1.4142135623731*e_em*IT_0208;
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0149*IT_0210;
    const ccomplex_t IT_0212 = IT_0090*IT_0115*IT_0211;
    const ccomplex_t IT_0213 = IT_0086*IT_0088;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = m_C_1*IT_0214;
    const ccomplex_t IT_0216 = m_C_1*IT_0070;
    const ccomplex_t IT_0217 = IT_0113*IT_0210;
    const ccomplex_t IT_0218 = IT_0115*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = m_C_1*IT_0219;
    const ccomplex_t IT_0221 = IT_0199 + IT_0200 + IT_0202 + IT_0203 + IT_0212
       + IT_0215 + IT_0216 + IT_0220;
    const ccomplex_t IT_0222 = 2*IT_0221;
    const ccomplex_t IT_0223 = s_12*s_24;
    const ccomplex_t IT_0224 = s_14*IT_0024;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = IT_0223 + IT_0225;
    const ccomplex_t IT_0227 = s_14*IT_0192;
    const ccomplex_t IT_0228 = IT_0072*IT_0227;
    const ccomplex_t IT_0229 = (-0.25)*IT_0228;
    const ccomplex_t IT_0230 = s_14 + IT_0229;
    const ccomplex_t IT_0231 = m_C_1*IT_0230;
    const ccomplex_t IT_0232 = (-4)*IT_0231;
    const ccomplex_t IT_0233 = (-2)*IT_0154;
    const ccomplex_t IT_0234 = s_12*s_14;
    const ccomplex_t IT_0235 = s_24*IT_0024;
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = IT_0234 + IT_0236;
    const ccomplex_t IT_0238 = (-2)*conj(IT_0071);
    const ccomplex_t IT_0239 = 2*conj(IT_0157);
    const ccomplex_t IT_0240 = m_C_1*IT_0167;
    const ccomplex_t IT_0241 = (-2)*conj(IT_0154);
    const ccomplex_t IT_0242 = IT_0177 + IT_0186 + 2*IT_0214 + 2*IT_0219;
    const ccomplex_t IT_0243 = 2*IT_0242;
    const ccomplex_t IT_0244 = 2*conj(IT_0242);
    const ccomplex_t IT_0245 = m_C_1*IT_0161;
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = (-2)*IT_0187;
    const ccomplex_t IT_0248 = 2*IT_0154;
    const ccomplex_t IT_0249 = IT_0024*IT_0161;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = (-2)*IT_0157;
    const ccomplex_t IT_0252 = -IT_0126 + -IT_0137;
    const ccomplex_t IT_0253 = 2*conj(IT_0154);
    const ccomplex_t IT_0254 = 2*conj(IT_0071);
    const ccomplex_t IT_0255 = (-2)*conj(IT_0187);
    const ccomplex_t IT_0256 = (-2)*conj(IT_0157);
    const ccomplex_t IT_0257 = (-0.25)*IT_0194;
    const ccomplex_t IT_0258 = s_24 + IT_0257;
    const ccomplex_t IT_0259 = m_C_1*IT_0258;
    const ccomplex_t IT_0260 = (-4)*IT_0259;
    const ccomplex_t IT_0261 = 0.5*IT_0228;
    const ccomplex_t IT_0262 = s_14 + IT_0261;
    const ccomplex_t IT_0263 = m_C_1*IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0263;
    const ccomplex_t IT_0265 = (-2)*IT_0221;
    const ccomplex_t IT_0266 = IT_0024*IT_0192;
    const ccomplex_t IT_0267 = 3*IT_0266;
    const ccomplex_t IT_0268 = (-2)*IT_0191;
    const ccomplex_t IT_0269 = -IT_0214 + -IT_0219;
    const ccomplex_t IT_0270 = (-3)*IT_0266;
    const ccomplex_t IT_0271 = s_12*IT_0192;
    const ccomplex_t IT_0272 = s_14*s_24;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = IT_0072*IT_0274;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0027 + IT_0070;
    const ccomplex_t IT_0278 = pow(m_Z, 4);
    const ccomplex_t IT_0279 = s_12*IT_0278;
    const ccomplex_t IT_0280 = IT_0072*IT_0279;
    const ccomplex_t IT_0281 = (-0.25)*IT_0280;
    const ccomplex_t IT_0282 = IT_0272 + IT_0281;
    const ccomplex_t IT_0283 = 2*conj(IT_0221);
    const ccomplex_t IT_0284 = (-0.5)*IT_0280;
    const ccomplex_t IT_0285 = IT_0271 + IT_0273 + IT_0284;
    const ccomplex_t IT_0286 = 2*conj(IT_0277);
    const ccomplex_t IT_0287 = (-0.5)*IT_0237*IT_0254 + (-0.5)*IT_0226*IT_0256
       + (-0.5)*IT_0241*IT_0264 + IT_0267*conj(IT_0269) + conj(IT_0252)*IT_0270 
      + 4*conj(IT_0191)*IT_0282 + (-0.5)*IT_0260*IT_0283 + IT_0285*IT_0286;
    const ccomplex_t IT_0288 = (-2)*conj(IT_0277);
    const ccomplex_t IT_0289 = (-2)*conj(IT_0191);
    const ccomplex_t IT_0290 = (-2)*conj(IT_0221);
    const ccomplex_t IT_0291 = 0.5*IT_0264;
    const ccomplex_t IT_0292 = 0.5*IT_0260;
    const ccomplex_t IT_0293 = conj(IT_0187)*IT_0237 + -IT_0226*conj(IT_0242) 
      + (-4)*conj(IT_0269)*IT_0282 + (-2)*conj(IT_0252)*IT_0285 + 0.5*IT_0270
      *IT_0288 + 0.5*IT_0267*IT_0289 + IT_0290*IT_0291 + IT_0253*IT_0292;
    const ccomplex_t IT_0294 = (-0.25)*IT_0232;
    const ccomplex_t IT_0295 = (-0.25)*IT_0198;
    const ccomplex_t IT_0296 = (-0.25)*IT_0237;
    const ccomplex_t IT_0297 = (-0.25)*IT_0226;
    const ccomplex_t IT_0298 = (-8)*IT_0269*(conj(IT_0252)*IT_0282 + 1./2*conj
      (IT_0269)*IT_0285 + (-0.125)*IT_0267*IT_0288 + (-0.125)*IT_0270*IT_0289 +
       1./2*IT_0290*IT_0294 + 1./2*IT_0253*IT_0295 + 1./2*IT_0255*IT_0296 + 1./2
      *IT_0244*IT_0297);
    const ccomplex_t IT_0299 = 2*IT_0071*conj(IT_0071)*IT_0077 + conj(IT_0154)
      *IT_0155*IT_0156 + conj(IT_0157)*(IT_0163*IT_0164 + IT_0156*IT_0168) +
       conj(IT_0187)*(IT_0077*IT_0188 + IT_0164*IT_0189 + IT_0156*IT_0190) +
       conj(IT_0191)*(IT_0198*IT_0222 + IT_0164*IT_0226 + IT_0232*IT_0233 +
       IT_0156*IT_0237) + IT_0157*IT_0168*IT_0238 + IT_0187*(IT_0190*IT_0238 +
       IT_0189*IT_0239) + IT_0240*(conj(IT_0154)*IT_0164 + conj(IT_0187)*IT_0233
       + IT_0154*IT_0239 + IT_0187*IT_0241 + conj(IT_0221)*(IT_0156 + IT_0243) +
       IT_0221*(IT_0238 + IT_0244)) + conj(IT_0221)*(IT_0164*IT_0246 + IT_0155
      *IT_0247) + 2*conj(IT_0242)*(IT_0071*IT_0189 + 0.5*IT_0163*IT_0243 + 0.5
      *IT_0168*IT_0247 + 0.5*IT_0246*IT_0248 + 0.5*IT_0250*IT_0251) + (-2)
      *IT_0242*(IT_0226*conj(IT_0252) + (-0.5)*IT_0246*IT_0253 + (-0.5)*IT_0189
      *IT_0254 + (-0.5)*IT_0168*IT_0255 + (-0.5)*IT_0250*IT_0256) + conj(IT_0252
      )*(IT_0188*IT_0237 + IT_0248*IT_0260 + IT_0264*IT_0265 + IT_0267*IT_0268) 
      + conj(IT_0269)*(IT_0226*IT_0243 + IT_0237*IT_0247 + IT_0198*IT_0248 +
       IT_0232*IT_0265 + IT_0268*IT_0270) + 6*IT_0154*(IT_0024*conj(IT_0221) +
       0.166666666666667*IT_0155*IT_0238 + 0.666666666666667*conj(IT_0154)*(s_12
       + IT_0276)) + 6*IT_0221*(IT_0024*conj(IT_0154) + 0.166666666666667
      *IT_0239*IT_0246 + 0.166666666666667*IT_0155*IT_0255 + 0.666666666666667
      *conj(IT_0221)*(s_12 + IT_0276)) + (-2)*IT_0277*IT_0287 + IT_0286*(IT_0071
      *IT_0237 + 0.5*IT_0226*IT_0251 + IT_0233*IT_0291 + IT_0222*IT_0292) + 2
      *IT_0252*IT_0293 + (-4)*IT_0191*(conj(IT_0191)*IT_0285 + IT_0282*IT_0286 +
       IT_0241*IT_0294 + IT_0283*IT_0295 + IT_0238*IT_0296 + IT_0239*IT_0297) +
       IT_0298;
    return create_ccomplex_return(IT_0299);
}

