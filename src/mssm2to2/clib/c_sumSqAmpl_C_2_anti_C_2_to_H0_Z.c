#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_H0_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = U_d1*V_Wp2;
    const ccomplex_t IT_0002 = IT_0000*IT_0001;
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = V_u2*U_Wm1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (IT_0002 + IT_0005)*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = IT_0006*IT_0013;
    const ccomplex_t IT_0015 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0007*IT_0012;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0016 + -IT_0018 + (-2)
      *IT_0020);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = IT_0011*IT_0022;
    const ccomplex_t IT_0024 = pow(m_H0, 2);
    const ccomplex_t IT_0025 = pow(m_C_2, 2);
    const ccomplex_t IT_0026 = cpow((-2)*s_13 + IT_0024 + IT_0025 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0027 = IT_0023*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow((-2)*s_13 + IT_0024 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0030 = U_d2*V_Wp2;
    const ccomplex_t IT_0031 = IT_0000*IT_0030;
    const ccomplex_t IT_0032 = V_u2*U_Wm2;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = IT_0007*(IT_0031 + IT_0033);
    const ccomplex_t IT_0035 = 1.4142135623731*e_em*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = 0.5*IT_0036;
    const ccomplex_t IT_0038 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0039 = IT_0014*IT_0038;
    const ccomplex_t IT_0040 = IT_0017*IT_0038;
    const ccomplex_t IT_0041 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + -IT_0040 + (-2)
      *IT_0042);
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0037*IT_0044;
    const ccomplex_t IT_0046 = IT_0029*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = pow(m_Z, 2);
    const ccomplex_t IT_0049 = s_12*IT_0048;
    const ccomplex_t IT_0050 = pow(m_Z, -2);
    const ccomplex_t IT_0051 = pow(m_Z, 4);
    const ccomplex_t IT_0052 = s_12*IT_0051;
    const ccomplex_t IT_0053 = IT_0050*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = s_14*s_24;
    const ccomplex_t IT_0056 = (-2)*IT_0055;
    const ccomplex_t IT_0057 = IT_0049 + IT_0054 + IT_0056;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0024 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0059 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0060 = IT_0014*IT_0059;
    const ccomplex_t IT_0061 = IT_0017*IT_0059;
    const ccomplex_t IT_0062 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0063 = IT_0017*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*(IT_0060 + -IT_0061 + (-2)
      *IT_0063);
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = IT_0037*IT_0065;
    const ccomplex_t IT_0067 = IT_0058*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = U_d2*V_Wp1;
    const ccomplex_t IT_0070 = IT_0000*IT_0069;
    const ccomplex_t IT_0071 = V_u1*U_Wm2;
    const ccomplex_t IT_0072 = IT_0003*IT_0071;
    const ccomplex_t IT_0073 = IT_0007*(IT_0070 + IT_0072);
    const ccomplex_t IT_0074 = 1.4142135623731*e_em*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0078 = IT_0014*IT_0077;
    const ccomplex_t IT_0079 = IT_0017*IT_0077;
    const ccomplex_t IT_0080 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0081 = IT_0017*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + -IT_0079 + (-2)
      *IT_0081);
    const ccomplex_t IT_0083 = (-0.5)*IT_0082;
    const ccomplex_t IT_0084 = IT_0076*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_23 + IT_0024 + IT_0025 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0068 + IT_0087;
    const ccomplex_t IT_0089 = (-0.25)*IT_0053;
    const ccomplex_t IT_0090 = IT_0055 + IT_0089;
    const ccomplex_t IT_0091 = s_24*IT_0048;
    const ccomplex_t IT_0092 = IT_0050*IT_0091;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = s_24 + IT_0093;
    const ccomplex_t IT_0095 = m_C_2*IT_0094;
    const ccomplex_t IT_0096 = 2*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0098 = IT_0045*IT_0058*IT_0097;
    const ccomplex_t IT_0099 = cos(beta);
    const ccomplex_t IT_0100 = cpow(IT_0012, -4);
    const ccomplex_t IT_0101 = IT_0000*IT_0006*IT_0099*IT_0100;
    const ccomplex_t IT_0102 = sin(beta);
    const ccomplex_t IT_0103 = IT_0003*IT_0006*IT_0100*IT_0102;
    const ccomplex_t IT_0104 = cpow(IT_0006, -3);
    const ccomplex_t IT_0105 = IT_0000*IT_0099*IT_0104;
    const ccomplex_t IT_0106 = IT_0003*IT_0102*IT_0104;
    const ccomplex_t IT_0107 = cpow(IT_0012, -2);
    const ccomplex_t IT_0108 = IT_0000*IT_0007*IT_0099*IT_0107;
    const ccomplex_t IT_0109 = IT_0003*IT_0007*IT_0102*IT_0107;
    const ccomplex_t IT_0110 = cpow(IT_0006, -2);
    const ccomplex_t IT_0111 = cpow(IT_0107 + IT_0110, -1);
    const ccomplex_t IT_0112 = (IT_0101 + IT_0103 + IT_0105 + IT_0106 + 2
      *IT_0108 + 2*IT_0109)*IT_0111;
    const ccomplex_t IT_0113 = m_W*e_em;
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = IT_0044*IT_0115;
    const ccomplex_t IT_0117 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0025 + -reg_prop, -1);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0121 = IT_0014*IT_0120;
    const ccomplex_t IT_0122 = IT_0017*IT_0120;
    const ccomplex_t IT_0123 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0124 = IT_0017*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + -IT_0122 + (-2)
      *IT_0124);
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = IT_0076*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0129 = IT_0085*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0131 = IT_0000*IT_0130;
    const ccomplex_t IT_0132 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0133 = IT_0003*IT_0132;
    const ccomplex_t IT_0134 = IT_0007*(IT_0131 + IT_0133);
    const ccomplex_t IT_0135 = 1.4142135623731*e_em*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0044*IT_0137;
    const ccomplex_t IT_0139 = IT_0029*IT_0097*IT_0138;
    const ccomplex_t IT_0140 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0141 = IT_0000*IT_0140;
    const ccomplex_t IT_0142 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0143 = IT_0003*IT_0142;
    const ccomplex_t IT_0144 = IT_0007*(IT_0141 + IT_0143);
    const ccomplex_t IT_0145 = 1.4142135623731*e_em*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = 0.5*IT_0146;
    const ccomplex_t IT_0148 = IT_0022*IT_0147;
    const ccomplex_t IT_0149 = IT_0026*IT_0128*IT_0148;
    const ccomplex_t IT_0150 = m_C_2*IT_0068;
    const ccomplex_t IT_0151 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0152 = IT_0014*IT_0151;
    const ccomplex_t IT_0153 = IT_0017*IT_0151;
    const ccomplex_t IT_0154 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0155 = IT_0017*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0152 + -IT_0153 + (-2)
      *IT_0155);
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = IT_0147*IT_0157;
    const ccomplex_t IT_0159 = IT_0026*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = m_C_2*IT_0160;
    const ccomplex_t IT_0162 = m_C_2*IT_0087;
    const ccomplex_t IT_0163 = IT_0065*IT_0137;
    const ccomplex_t IT_0164 = IT_0029*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = m_C_2*IT_0165;
    const ccomplex_t IT_0167 = IT_0098 + -IT_0119 + IT_0129 + IT_0139 +
       IT_0149 + IT_0150 + IT_0161 + IT_0162 + IT_0166;
    const ccomplex_t IT_0168 = -e_em*IT_0007*IT_0013*(IT_0003*IT_0099 + 
      -IT_0000*IT_0102);
    const ccomplex_t IT_0169 = IT_0032*IT_0099;
    const ccomplex_t IT_0170 = IT_0030*IT_0102;
    const ccomplex_t IT_0171 = IT_0007*(IT_0169 + IT_0170);
    const ccomplex_t IT_0172 = 1.4142135623731*e_em*IT_0171;
    const ccomplex_t IT_0173 = 0.5*IT_0172;
    const ccomplex_t IT_0174 = cpow((-2)*s_12 + (-2)*IT_0025 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = IT_0168*IT_0176;
    const ccomplex_t IT_0178 = e_em*IT_0007*IT_0013*(IT_0000*IT_0099 + IT_0003
      *IT_0102);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0030*IT_0099;
    const ccomplex_t IT_0181 = IT_0032*IT_0102;
    const ccomplex_t IT_0182 = IT_0007*(IT_0180 + -IT_0181);
    const ccomplex_t IT_0183 = 1.4142135623731*e_em*IT_0182;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = cpow((-2)*s_12 + (-2)*IT_0025 + IT_0048 + 
      -reg_prop, -1);
    const ccomplex_t IT_0186 = IT_0184*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = IT_0179*IT_0187;
    const ccomplex_t IT_0189 = (-2)*IT_0028 + (-2)*IT_0047 + -IT_0177 + 
      -IT_0188;
    const ccomplex_t IT_0190 = s_12*s_24;
    const ccomplex_t IT_0191 = s_14*IT_0025;
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = IT_0190 + IT_0192;
    const ccomplex_t IT_0194 = 2*conj(IT_0028);
    const ccomplex_t IT_0195 = 2*conj(IT_0047);
    const ccomplex_t IT_0196 = IT_0025*IT_0048;
    const ccomplex_t IT_0197 = (-3)*IT_0196;
    const ccomplex_t IT_0198 = -IT_0160 + -IT_0165;
    const ccomplex_t IT_0199 = s_14*IT_0048;
    const ccomplex_t IT_0200 = IT_0050*IT_0199;
    const ccomplex_t IT_0201 = (-0.25)*IT_0200;
    const ccomplex_t IT_0202 = s_14 + IT_0201;
    const ccomplex_t IT_0203 = m_C_2*IT_0202;
    const ccomplex_t IT_0204 = (-4)*IT_0203;
    const ccomplex_t IT_0205 = IT_0058*IT_0097*IT_0163;
    const ccomplex_t IT_0206 = IT_0029*IT_0066*IT_0097;
    const ccomplex_t IT_0207 = IT_0065*IT_0115;
    const ccomplex_t IT_0208 = IT_0117*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = IT_0011*IT_0157;
    const ccomplex_t IT_0211 = IT_0026*IT_0128*IT_0210;
    const ccomplex_t IT_0212 = m_C_2*IT_0028;
    const ccomplex_t IT_0213 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0214 = IT_0000*IT_0213;
    const ccomplex_t IT_0215 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0216 = IT_0003*IT_0215;
    const ccomplex_t IT_0217 = IT_0007*(IT_0214 + IT_0216);
    const ccomplex_t IT_0218 = 1.4142135623731*e_em*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = IT_0083*IT_0220;
    const ccomplex_t IT_0222 = IT_0085*IT_0128*IT_0221;
    const ccomplex_t IT_0223 = IT_0126*IT_0220;
    const ccomplex_t IT_0224 = IT_0085*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = m_C_2*IT_0225;
    const ccomplex_t IT_0227 = IT_0058*IT_0138;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = m_C_2*IT_0228;
    const ccomplex_t IT_0230 = m_C_2*IT_0047;
    const ccomplex_t IT_0231 = -IT_0205 + -IT_0206 + IT_0209 + -IT_0211 + 
      -IT_0212 + -IT_0222 + -IT_0226 + -IT_0229 + -IT_0230;
    const ccomplex_t IT_0232 = 3*IT_0196;
    const ccomplex_t IT_0233 = -IT_0225 + -IT_0228;
    const ccomplex_t IT_0234 = (-2)*IT_0068 + (-2)*IT_0087 + -IT_0177 + 
      -IT_0188;
    const ccomplex_t IT_0235 = s_12*s_14;
    const ccomplex_t IT_0236 = s_24*IT_0025;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = (-2)*conj(IT_0028);
    const ccomplex_t IT_0240 = (-2)*conj(IT_0047);
    const ccomplex_t IT_0241 = IT_0099*IT_0132;
    const ccomplex_t IT_0242 = IT_0102*IT_0130;
    const ccomplex_t IT_0243 = IT_0007*(IT_0241 + IT_0242);
    const ccomplex_t IT_0244 = 1.4142135623731*e_em*IT_0243;
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = IT_0174*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = IT_0168*IT_0247;
    const ccomplex_t IT_0249 = IT_0099*IT_0130;
    const ccomplex_t IT_0250 = IT_0102*IT_0132;
    const ccomplex_t IT_0251 = IT_0007*(IT_0249 + -IT_0250);
    const ccomplex_t IT_0252 = 1.4142135623731*e_em*IT_0251;
    const ccomplex_t IT_0253 = 0.5*IT_0252;
    const ccomplex_t IT_0254 = IT_0185*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = IT_0179*IT_0255;
    const ccomplex_t IT_0257 = 2*IT_0160 + 2*IT_0165 + -IT_0248 + -IT_0256;
    const ccomplex_t IT_0258 = 2*IT_0193;
    const ccomplex_t IT_0259 = pow(s_24, 2);
    const ccomplex_t IT_0260 = IT_0050*IT_0259;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0025 + IT_0261;
    const ccomplex_t IT_0263 = s_12*IT_0262;
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = 2*IT_0264;
    const ccomplex_t IT_0266 = s_14*s_24*IT_0050;
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = s_12 + IT_0267;
    const ccomplex_t IT_0269 = IT_0025*IT_0268;
    const ccomplex_t IT_0270 = 2*IT_0269;
    const ccomplex_t IT_0271 = m_C_2*IT_0268;
    const ccomplex_t IT_0272 = 2*IT_0271;
    const ccomplex_t IT_0273 = m_C_2*IT_0262;
    const ccomplex_t IT_0274 = -IT_0273;
    const ccomplex_t IT_0275 = 2*IT_0274;
    const ccomplex_t IT_0276 = (-2)*IT_0193;
    const ccomplex_t IT_0277 = 2*IT_0225 + 2*IT_0228 + -IT_0248 + -IT_0256;
    const ccomplex_t IT_0278 = s_12*IT_0268;
    const ccomplex_t IT_0279 = (-2)*IT_0278;
    const ccomplex_t IT_0280 = IT_0025*IT_0262;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = (-2)*IT_0281;
    const ccomplex_t IT_0283 = (-2)*IT_0028;
    const ccomplex_t IT_0284 = (-2)*IT_0047;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0096;
    const ccomplex_t IT_0287 = (-0.25)*IT_0092;
    const ccomplex_t IT_0288 = s_24 + IT_0287;
    const ccomplex_t IT_0289 = m_C_2*IT_0288;
    const ccomplex_t IT_0290 = (-4)*IT_0289;
    const ccomplex_t IT_0291 = 2*IT_0290;
    const ccomplex_t IT_0292 = IT_0049 + IT_0056;
    const ccomplex_t IT_0293 = IT_0050*IT_0292;
    const ccomplex_t IT_0294 = 4*s_12 + (-2)*IT_0293;
    const ccomplex_t IT_0295 = 6*IT_0025;
    const ccomplex_t IT_0296 = pow(s_14, 2);
    const ccomplex_t IT_0297 = IT_0050*IT_0296;
    const ccomplex_t IT_0298 = -IT_0025;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = m_C_2*IT_0299;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = (-2)*IT_0271;
    const ccomplex_t IT_0303 = 0.5*IT_0200;
    const ccomplex_t IT_0304 = s_14 + IT_0303;
    const ccomplex_t IT_0305 = m_C_2*IT_0304;
    const ccomplex_t IT_0306 = 2*IT_0305;
    const ccomplex_t IT_0307 = (-2)*IT_0306;
    const ccomplex_t IT_0308 = conj(IT_0231)*(IT_0189*IT_0272 + IT_0257
      *IT_0275 + IT_0204*IT_0285 + IT_0198*IT_0286 + IT_0233*IT_0291 + IT_0231
      *IT_0294 + IT_0167*IT_0295 + IT_0234*IT_0301 + IT_0277*IT_0302 + IT_0088
      *IT_0307);
    const ccomplex_t IT_0309 = 2*IT_0238;
    const ccomplex_t IT_0310 = (-8)*IT_0090;
    const ccomplex_t IT_0311 = (-4)*IT_0057;
    const ccomplex_t IT_0312 = (-2)*IT_0197;
    const ccomplex_t IT_0313 = (-2)*IT_0238;
    const ccomplex_t IT_0314 = (-2)*IT_0204;
    const ccomplex_t IT_0315 = (-2)*IT_0232;
    const ccomplex_t IT_0316 = (IT_0028 + IT_0047)*IT_0090 + (-0.125)*IT_0189
      *IT_0276 + (-0.125)*IT_0167*IT_0291 + (-0.125)*IT_0231*IT_0307 + (-0.125)
      *IT_0234*IT_0309 + (-0.125)*IT_0088*IT_0311 + (-0.125)*IT_0233*IT_0312 + (
      -0.125)*IT_0198*IT_0315;
    const ccomplex_t IT_0317 = 2*IT_0047;
    const ccomplex_t IT_0318 = 0.5*IT_0317;
    const ccomplex_t IT_0319 = IT_0028 + IT_0318;
    const ccomplex_t IT_0320 = 0.5*IT_0234;
    const ccomplex_t IT_0321 = 0.5*IT_0088;
    const ccomplex_t IT_0322 = 0.5*IT_0189;
    const ccomplex_t IT_0323 = 0.5*IT_0277;
    const ccomplex_t IT_0324 = 0.5*IT_0257;
    const ccomplex_t IT_0325 = 0.5*IT_0231;
    const ccomplex_t IT_0326 = 0.5*IT_0167;
    const ccomplex_t IT_0327 = s_12*IT_0299;
    const ccomplex_t IT_0328 = IT_0025*IT_0299;
    const ccomplex_t IT_0329 = (-2)*IT_0328;
    const ccomplex_t IT_0330 = 0.5*IT_0198;
    const ccomplex_t IT_0331 = 0.5*IT_0233;
    const ccomplex_t IT_0332 = 2*conj(IT_0167)*(IT_0096*IT_0319 + IT_0302
      *IT_0320 + IT_0291*IT_0321 + IT_0275*IT_0322 + IT_0301*IT_0323 + IT_0272
      *IT_0324 + IT_0295*IT_0325 + IT_0294*IT_0326 + IT_0314*IT_0330 + IT_0307
      *IT_0331);
    const ccomplex_t IT_0333 = (-4)*(conj(IT_0028) + conj(IT_0047))*((IT_0028 
      + IT_0047)*IT_0057 + 2*IT_0088*IT_0090) + (IT_0096*IT_0167 + IT_0189
      *IT_0193)*(IT_0194 + IT_0195) + (IT_0197*IT_0198 + IT_0204*IT_0231 +
       IT_0232*IT_0233 + IT_0234*IT_0238)*(IT_0239 + IT_0240) + conj(IT_0257)*
      (IT_0198*IT_0258 + IT_0257*IT_0265 + IT_0234*IT_0270 + IT_0167*IT_0272 +
       IT_0231*IT_0275 + IT_0233*IT_0276 + IT_0277*IT_0279 + IT_0189*IT_0282) +
       IT_0308 + conj(IT_0233)*(IT_0257*IT_0276 + IT_0232*IT_0285 + IT_0231
      *IT_0291 + IT_0167*IT_0307 + IT_0277*IT_0309 + IT_0198*IT_0310 + IT_0233
      *IT_0311 + IT_0088*IT_0312) + conj(IT_0198)*(IT_0257*IT_0258 + IT_0197
      *IT_0285 + IT_0231*IT_0286 + IT_0233*IT_0310 + IT_0198*IT_0311 + IT_0277
      *IT_0313 + IT_0167*IT_0314 + IT_0088*IT_0315) + (-8)*conj(IT_0088)*IT_0316
       + 2*conj(IT_0189)*(IT_0193*IT_0319 + IT_0279*IT_0320 + IT_0276*IT_0321 +
       IT_0265*IT_0322 + IT_0270*IT_0323 + IT_0282*IT_0324 + IT_0272*IT_0325 +
       IT_0275*IT_0326) + 2*conj(IT_0234)*(0.5*IT_0238*IT_0285 + IT_0309*IT_0321
       + IT_0279*IT_0322 + IT_0270*IT_0324 + IT_0301*IT_0325 + IT_0302*IT_0326 +
       IT_0234*IT_0327 + IT_0323*IT_0329) + 2*conj(IT_0277)*(IT_0270*IT_0322 +
       IT_0279*IT_0324 + IT_0302*IT_0325 + IT_0301*IT_0326 + IT_0277*IT_0327 +
       IT_0320*IT_0329 + IT_0313*IT_0330 + IT_0309*IT_0331) + IT_0332;
    return create_ccomplex_return(IT_0333);
}

