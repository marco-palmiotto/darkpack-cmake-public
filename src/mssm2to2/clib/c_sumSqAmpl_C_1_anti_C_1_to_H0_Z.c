#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_H0_Z(
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
    const ccomplex_t IT_0000 = pow(m_H0, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = U_d1*V_Wp1;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = V_u1*U_Wm1;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (IT_0004 + IT_0007)*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*IT_0011;
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0009*IT_0014;
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0018 + -IT_0020 + (-2)
      *IT_0022);
    const ccomplex_t IT_0024 = (-0.5)*IT_0023;
    const ccomplex_t IT_0025 = IT_0013*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = U_d1*V_Wp2;
    const ccomplex_t IT_0029 = IT_0002*IT_0028;
    const ccomplex_t IT_0030 = V_u2*U_Wm1;
    const ccomplex_t IT_0031 = IT_0005*IT_0030;
    const ccomplex_t IT_0032 = IT_0009*(IT_0029 + IT_0031);
    const ccomplex_t IT_0033 = 1.4142135623731*e_em*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0037 = IT_0016*IT_0036;
    const ccomplex_t IT_0038 = IT_0019*IT_0036;
    const ccomplex_t IT_0039 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0040 = IT_0019*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0038 + (-2)
      *IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = pow(m_C_1, 2);
    const ccomplex_t IT_0045 = cpow((-2)*s_23 + IT_0000 + IT_0044 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0046 = IT_0043*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = sin(beta);
    const ccomplex_t IT_0049 = cos(beta);
    const ccomplex_t IT_0050 = e_em*IT_0009*IT_0015*(IT_0002*IT_0048 + 
      -IT_0005*IT_0049);
    const ccomplex_t IT_0051 = IT_0006*IT_0049;
    const ccomplex_t IT_0052 = IT_0003*IT_0048;
    const ccomplex_t IT_0053 = IT_0009*(IT_0051 + IT_0052);
    const ccomplex_t IT_0054 = 1.4142135623731*e_em*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_12 + (-2)*IT_0044 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0050*IT_0058;
    const ccomplex_t IT_0060 = e_em*IT_0009*IT_0015*(IT_0005*IT_0048 + IT_0002
      *IT_0049);
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = IT_0003*IT_0049;
    const ccomplex_t IT_0063 = IT_0006*IT_0048;
    const ccomplex_t IT_0064 = IT_0009*(IT_0062 + -IT_0063);
    const ccomplex_t IT_0065 = 1.4142135623731*e_em*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = pow(m_Z, 2);
    const ccomplex_t IT_0068 = cpow((-2)*s_12 + (-2)*IT_0044 + IT_0067 + 
      -reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = IT_0061*IT_0070;
    const ccomplex_t IT_0072 = (-2)*IT_0027 + (-2)*IT_0047 + -IT_0059 + 
      -IT_0071;
    const ccomplex_t IT_0073 = pow(m_Z, -2);
    const ccomplex_t IT_0074 = pow(s_14, 2);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = -IT_0044;
    const ccomplex_t IT_0077 = IT_0075 + IT_0076;
    const ccomplex_t IT_0078 = s_12*IT_0077;
    const ccomplex_t IT_0079 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0080 = IT_0002*IT_0079;
    const ccomplex_t IT_0081 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0082 = IT_0005*IT_0081;
    const ccomplex_t IT_0083 = IT_0009*(IT_0080 + IT_0082);
    const ccomplex_t IT_0084 = 1.4142135623731*e_em*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0088 = IT_0016*IT_0087;
    const ccomplex_t IT_0089 = IT_0019*IT_0087;
    const ccomplex_t IT_0090 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0091 = IT_0019*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + -IT_0089 + (-2)
      *IT_0091);
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = IT_0086*IT_0093;
    const ccomplex_t IT_0095 = IT_0045*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0098 = IT_0049*IT_0097;
    const ccomplex_t IT_0099 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0100 = IT_0048*IT_0099;
    const ccomplex_t IT_0101 = IT_0009*(IT_0098 + IT_0100);
    const ccomplex_t IT_0102 = 1.4142135623731*e_em*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = IT_0056*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0050*IT_0105;
    const ccomplex_t IT_0107 = IT_0002*IT_0099;
    const ccomplex_t IT_0108 = IT_0005*IT_0097;
    const ccomplex_t IT_0109 = IT_0009*(IT_0107 + IT_0108);
    const ccomplex_t IT_0110 = 1.4142135623731*e_em*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = 0.5*IT_0111;
    const ccomplex_t IT_0113 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0114 = IT_0016*IT_0113;
    const ccomplex_t IT_0115 = IT_0019*IT_0113;
    const ccomplex_t IT_0116 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0117 = IT_0019*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0114 + -IT_0115 + (-2)
      *IT_0117);
    const ccomplex_t IT_0119 = 0.5*IT_0118;
    const ccomplex_t IT_0120 = IT_0112*IT_0119;
    const ccomplex_t IT_0121 = IT_0001*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0049*IT_0099;
    const ccomplex_t IT_0124 = IT_0048*IT_0097;
    const ccomplex_t IT_0125 = IT_0009*(IT_0123 + -IT_0124);
    const ccomplex_t IT_0126 = 1.4142135623731*e_em*IT_0125;
    const ccomplex_t IT_0127 = 0.5*IT_0126;
    const ccomplex_t IT_0128 = IT_0068*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0061*IT_0129;
    const ccomplex_t IT_0131 = 2*IT_0096 + -IT_0106 + 2*IT_0122 + -IT_0130;
    const ccomplex_t IT_0132 = s_14*s_24*IT_0073;
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = s_12 + IT_0133;
    const ccomplex_t IT_0135 = IT_0044*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0137 = IT_0013*IT_0119;
    const ccomplex_t IT_0138 = IT_0001*IT_0136*IT_0137;
    const ccomplex_t IT_0139 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0140 = IT_0002*IT_0139;
    const ccomplex_t IT_0141 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0142 = IT_0005*IT_0141;
    const ccomplex_t IT_0143 = IT_0009*(IT_0140 + IT_0142);
    const ccomplex_t IT_0144 = 1.4142135623731*e_em*IT_0143;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*IT_0144;
    const ccomplex_t IT_0146 = 0.5*IT_0145;
    const ccomplex_t IT_0147 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0148 = IT_0016*IT_0147;
    const ccomplex_t IT_0149 = IT_0019*IT_0147;
    const ccomplex_t IT_0150 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0151 = IT_0019*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0148 + -IT_0149 + (-2)
      *IT_0151);
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = IT_0146*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0156 = cpow((-2)*s_13 + IT_0000 + IT_0044 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0157 = IT_0154*IT_0155*IT_0156;
    const ccomplex_t IT_0158 = IT_0035*IT_0093;
    const ccomplex_t IT_0159 = IT_0045*IT_0155*IT_0158;
    const ccomplex_t IT_0160 = m_C_1*IT_0027;
    const ccomplex_t IT_0161 = m_C_1*IT_0047;
    const ccomplex_t IT_0162 = cpow(IT_0014, -4);
    const ccomplex_t IT_0163 = IT_0002*IT_0008*IT_0049*IT_0162;
    const ccomplex_t IT_0164 = IT_0005*IT_0008*IT_0048*IT_0162;
    const ccomplex_t IT_0165 = cpow(IT_0008, -3);
    const ccomplex_t IT_0166 = IT_0002*IT_0049*IT_0165;
    const ccomplex_t IT_0167 = IT_0005*IT_0048*IT_0165;
    const ccomplex_t IT_0168 = cpow(IT_0014, -2);
    const ccomplex_t IT_0169 = IT_0002*IT_0009*IT_0049*IT_0168;
    const ccomplex_t IT_0170 = IT_0005*IT_0009*IT_0048*IT_0168;
    const ccomplex_t IT_0171 = cpow(IT_0008, -2);
    const ccomplex_t IT_0172 = cpow(IT_0168 + IT_0171, -1);
    const ccomplex_t IT_0173 = (IT_0163 + IT_0164 + IT_0166 + IT_0167 + 2
      *IT_0169 + 2*IT_0170)*IT_0172;
    const ccomplex_t IT_0174 = m_W*e_em;
    const ccomplex_t IT_0175 = IT_0173*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = IT_0119*IT_0176;
    const ccomplex_t IT_0178 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0044 + -reg_prop, -1);
    const ccomplex_t IT_0179 = IT_0177*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0182 = IT_0024*IT_0112;
    const ccomplex_t IT_0183 = IT_0181*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = m_C_1*IT_0184;
    const ccomplex_t IT_0186 = IT_0120*IT_0136*IT_0181;
    const ccomplex_t IT_0187 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0188 = IT_0016*IT_0187;
    const ccomplex_t IT_0189 = IT_0019*IT_0187;
    const ccomplex_t IT_0190 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0191 = IT_0019*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0188 + -IT_0189 + (-2)
      *IT_0191);
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0146*IT_0193;
    const ccomplex_t IT_0195 = IT_0156*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = m_C_1*IT_0196;
    const ccomplex_t IT_0198 = IT_0138 + IT_0157 + IT_0159 + IT_0160 + IT_0161
       + -IT_0180 + IT_0185 + IT_0186 + IT_0197;
    const ccomplex_t IT_0199 = pow(s_24, 2);
    const ccomplex_t IT_0200 = IT_0073*IT_0199;
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = IT_0044 + IT_0201;
    const ccomplex_t IT_0203 = m_C_1*IT_0202;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = U_d2*V_Wp1;
    const ccomplex_t IT_0206 = IT_0002*IT_0205;
    const ccomplex_t IT_0207 = V_u1*U_Wm2;
    const ccomplex_t IT_0208 = IT_0005*IT_0207;
    const ccomplex_t IT_0209 = IT_0009*(IT_0206 + IT_0208);
    const ccomplex_t IT_0210 = 1.4142135623731*e_em*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = IT_0153*IT_0212;
    const ccomplex_t IT_0214 = IT_0156*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0137*IT_0181;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = -IT_0059 + -IT_0071 + (-2)*IT_0215 + (-2)
      *IT_0217;
    const ccomplex_t IT_0219 = 2*conj(IT_0218);
    const ccomplex_t IT_0220 = IT_0001*IT_0136*IT_0182;
    const ccomplex_t IT_0221 = IT_0025*IT_0136*IT_0181;
    const ccomplex_t IT_0222 = IT_0024*IT_0176;
    const ccomplex_t IT_0223 = IT_0178*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0042*IT_0086;
    const ccomplex_t IT_0226 = IT_0045*IT_0155*IT_0225;
    const ccomplex_t IT_0227 = IT_0193*IT_0212;
    const ccomplex_t IT_0228 = IT_0155*IT_0156*IT_0227;
    const ccomplex_t IT_0229 = m_C_1*IT_0122;
    const ccomplex_t IT_0230 = m_C_1*IT_0215;
    const ccomplex_t IT_0231 = m_C_1*IT_0096;
    const ccomplex_t IT_0232 = m_C_1*IT_0217;
    const ccomplex_t IT_0233 = -IT_0220 + -IT_0221 + IT_0224 + -IT_0226 + 
      -IT_0228 + -IT_0229 + -IT_0230 + -IT_0231 + -IT_0232;
    const ccomplex_t IT_0234 = -IT_0106 + -IT_0130 + 2*IT_0184 + 2*IT_0196;
    const ccomplex_t IT_0235 = 2*conj(IT_0234);
    const ccomplex_t IT_0236 = s_12*IT_0202;
    const ccomplex_t IT_0237 = -IT_0236;
    const ccomplex_t IT_0238 = 2*IT_0218;
    const ccomplex_t IT_0239 = s_12*IT_0134;
    const ccomplex_t IT_0240 = (-2)*IT_0072;
    const ccomplex_t IT_0241 = 2*IT_0131;
    const ccomplex_t IT_0242 = IT_0044*IT_0077;
    const ccomplex_t IT_0243 = m_C_1*IT_0077;
    const ccomplex_t IT_0244 = (-2)*conj(IT_0072);
    const ccomplex_t IT_0245 = IT_0027 + IT_0047;
    const ccomplex_t IT_0246 = s_24*IT_0067;
    const ccomplex_t IT_0247 = IT_0073*IT_0246;
    const ccomplex_t IT_0248 = (-0.25)*IT_0247;
    const ccomplex_t IT_0249 = s_24 + IT_0248;
    const ccomplex_t IT_0250 = m_C_1*IT_0249;
    const ccomplex_t IT_0251 = (-4)*IT_0250;
    const ccomplex_t IT_0252 = 2*IT_0198;
    const ccomplex_t IT_0253 = s_14*IT_0067;
    const ccomplex_t IT_0254 = IT_0073*IT_0253;
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = s_14 + IT_0255;
    const ccomplex_t IT_0257 = m_C_1*IT_0256;
    const ccomplex_t IT_0258 = 2*IT_0257;
    const ccomplex_t IT_0259 = (-2)*IT_0233;
    const ccomplex_t IT_0260 = IT_0215 + IT_0217;
    const ccomplex_t IT_0261 = 0.5*IT_0247;
    const ccomplex_t IT_0262 = s_24 + IT_0261;
    const ccomplex_t IT_0263 = m_C_1*IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0263;
    const ccomplex_t IT_0265 = (-0.25)*IT_0254;
    const ccomplex_t IT_0266 = s_14 + IT_0265;
    const ccomplex_t IT_0267 = m_C_1*IT_0266;
    const ccomplex_t IT_0268 = (-4)*IT_0267;
    const ccomplex_t IT_0269 = (-2)*conj(IT_0131);
    const ccomplex_t IT_0270 = s_12*s_14;
    const ccomplex_t IT_0271 = s_24*IT_0044;
    const ccomplex_t IT_0272 = -IT_0271;
    const ccomplex_t IT_0273 = IT_0270 + IT_0272;
    const ccomplex_t IT_0274 = -IT_0096 + -IT_0122;
    const ccomplex_t IT_0275 = 2*conj(IT_0072);
    const ccomplex_t IT_0276 = -IT_0184 + -IT_0196;
    const ccomplex_t IT_0277 = (-2)*IT_0131;
    const ccomplex_t IT_0278 = IT_0072*conj(IT_0245) + 0.5*IT_0244*IT_0260 +
       0.5*IT_0240*conj(IT_0260) + conj(IT_0131)*IT_0274 + 0.5*IT_0241*conj
      (IT_0274) + 0.5*IT_0245*IT_0275 + 0.5*IT_0269*IT_0276 + 0.5*conj(IT_0276)
      *IT_0277;
    const ccomplex_t IT_0279 = m_C_1*IT_0134;
    const ccomplex_t IT_0280 = 2*IT_0234;
    const ccomplex_t IT_0281 = IT_0044*IT_0202;
    const ccomplex_t IT_0282 = -IT_0281;
    const ccomplex_t IT_0283 = (-2)*IT_0218;
    const ccomplex_t IT_0284 = s_12*s_24;
    const ccomplex_t IT_0285 = s_14*IT_0044;
    const ccomplex_t IT_0286 = -IT_0285;
    const ccomplex_t IT_0287 = IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = (-2)*conj(IT_0218);
    const ccomplex_t IT_0289 = (-0.5)*IT_0219*IT_0260 + (-0.5)*IT_0238*conj
      (IT_0260) + conj(IT_0234)*IT_0274 + IT_0234*conj(IT_0274) + (-0.5)*IT_0235
      *IT_0276 + (-0.5)*conj(IT_0276)*IT_0280 + (-0.5)*conj(IT_0245)*IT_0283 + (
      -0.5)*IT_0245*IT_0288;
    const ccomplex_t IT_0290 = 2*IT_0233;
    const ccomplex_t IT_0291 = IT_0044*IT_0067;
    const ccomplex_t IT_0292 = (-3)*IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0245;
    const ccomplex_t IT_0294 = (-2)*IT_0198;
    const ccomplex_t IT_0295 = 3*IT_0291;
    const ccomplex_t IT_0296 = (-2)*IT_0260;
    const ccomplex_t IT_0297 = s_12*IT_0067;
    const ccomplex_t IT_0298 = s_14*s_24;
    const ccomplex_t IT_0299 = (-2)*IT_0298;
    const ccomplex_t IT_0300 = IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = IT_0073*IT_0300;
    const ccomplex_t IT_0302 = (-0.5)*IT_0301;
    const ccomplex_t IT_0303 = 0.166666666666667*IT_0204;
    const ccomplex_t IT_0304 = pow(m_Z, 4);
    const ccomplex_t IT_0305 = s_12*IT_0304;
    const ccomplex_t IT_0306 = IT_0073*IT_0305;
    const ccomplex_t IT_0307 = (-0.5)*IT_0306;
    const ccomplex_t IT_0308 = IT_0297 + IT_0299 + IT_0307;
    const ccomplex_t IT_0309 = (-0.25)*IT_0306;
    const ccomplex_t IT_0310 = IT_0298 + IT_0309;
    const ccomplex_t IT_0311 = (-2)*conj(IT_0245);
    const ccomplex_t IT_0312 = (-2)*conj(IT_0260);
    const ccomplex_t IT_0313 = (-0.25)*IT_0292;
    const ccomplex_t IT_0314 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0315 = (-0.25)*IT_0268;
    const ccomplex_t IT_0316 = 2*conj(IT_0233);
    const ccomplex_t IT_0317 = (-0.25)*IT_0264;
    const ccomplex_t IT_0318 = (-2)*conj(IT_0233);
    const ccomplex_t IT_0319 = 2*conj(IT_0198);
    const ccomplex_t IT_0320 = (-0.25)*IT_0258;
    const ccomplex_t IT_0321 = (-0.25)*IT_0251;
    const ccomplex_t IT_0322 = 2*IT_0072*conj(IT_0072)*IT_0078 + (IT_0131
      *IT_0135 + IT_0198*IT_0204)*IT_0219 + IT_0204*IT_0233*IT_0235 + conj
      (IT_0218)*(IT_0237*IT_0238 + IT_0239*IT_0240) + conj(IT_0131)*(IT_0135
      *IT_0238 + IT_0078*IT_0241 + IT_0240*IT_0242) + (IT_0218*IT_0239 + IT_0131
      *IT_0242 + IT_0233*IT_0243)*IT_0244 + conj(IT_0245)*(IT_0251*IT_0252 +
       IT_0258*IT_0259) + conj(IT_0260)*(IT_0252*IT_0264 + IT_0259*IT_0268) +
       IT_0198*IT_0243*IT_0269 + 2*IT_0273*IT_0278 + IT_0279*(IT_0198*(IT_0235 +
       IT_0244) + IT_0233*(IT_0219 + IT_0269) + conj(IT_0233)*(IT_0238 + IT_0277
      ) + conj(IT_0198)*(IT_0240 + IT_0280)) + 2*conj(IT_0234)*(IT_0072*IT_0135 
      + 0.5*IT_0239*IT_0277 + 0.5*IT_0237*IT_0280 + 0.5*IT_0282*IT_0283) +
       IT_0234*(IT_0239*IT_0269 + IT_0135*IT_0275 + IT_0282*IT_0288) + (-2)
      *IT_0287*IT_0289 + conj(IT_0276)*(IT_0264*IT_0290 + IT_0268*IT_0294 +
       IT_0293*IT_0295 + IT_0292*IT_0296) + conj(IT_0274)*(IT_0251*IT_0290 +
       IT_0292*IT_0293 + IT_0258*IT_0294 + IT_0295*IT_0296) + 6*conj(IT_0198)*
      (IT_0044*IT_0233 + 0.166666666666667*IT_0243*IT_0277 + 0.666666666666667
      *IT_0198*(s_12 + IT_0302) + IT_0238*IT_0303) + 6*conj(IT_0233)*(IT_0044
      *IT_0198 + 0.166666666666667*IT_0240*IT_0243 + 0.666666666666667*IT_0233*
      (s_12 + IT_0302) + IT_0280*IT_0303) + (-4)*IT_0276*(conj(IT_0276)*IT_0308 
      + 2*conj(IT_0274)*IT_0310 + (-0.25)*IT_0295*IT_0311 + IT_0312*IT_0313 +
       IT_0314*IT_0315 + IT_0316*IT_0317) + (-4)*IT_0260*(conj(IT_0260)*IT_0308 
      + 2*conj(IT_0245)*IT_0310 + IT_0315*IT_0318 + IT_0317*IT_0319) + (-4)
      *IT_0274*(conj(IT_0274)*IT_0308 + 2*conj(IT_0276)*IT_0310 + (-0.25)
      *IT_0295*IT_0312 + IT_0311*IT_0313 + IT_0314*IT_0320 + IT_0316*IT_0321) + 
      (-4)*IT_0245*(conj(IT_0245)*IT_0308 + 2*conj(IT_0260)*IT_0310 + IT_0318
      *IT_0320 + IT_0319*IT_0321);
    return create_ccomplex_return(IT_0322);
}

