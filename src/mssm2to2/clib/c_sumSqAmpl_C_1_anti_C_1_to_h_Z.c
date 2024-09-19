#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_h_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
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
    const ccomplex_t IT_0000 = pow(m_h, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0002 = cos(alpha);
    const ccomplex_t IT_0003 = V_u1*U_Wm1;
    const ccomplex_t IT_0004 = IT_0002*IT_0003;
    const ccomplex_t IT_0005 = sin(alpha);
    const ccomplex_t IT_0006 = U_d1*V_Wp1;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (IT_0004 + -IT_0007)*IT_0009;
    const ccomplex_t IT_0011 = 1.4142135623731*e_em*IT_0010;
    const ccomplex_t IT_0012 = 0.5*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
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
    const ccomplex_t IT_0028 = V_u2*U_Wm1;
    const ccomplex_t IT_0029 = IT_0002*IT_0028;
    const ccomplex_t IT_0030 = U_d1*V_Wp2;
    const ccomplex_t IT_0031 = IT_0005*IT_0030;
    const ccomplex_t IT_0032 = IT_0009*(IT_0029 + -IT_0031);
    const ccomplex_t IT_0033 = 1.4142135623731*e_em*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
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
    const ccomplex_t IT_0048 = cos(beta);
    const ccomplex_t IT_0049 = sin(beta);
    const ccomplex_t IT_0050 = e_em*IT_0009*IT_0015*(IT_0002*IT_0048 + IT_0005
      *IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0003*IT_0048;
    const ccomplex_t IT_0053 = IT_0006*IT_0049;
    const ccomplex_t IT_0054 = IT_0009*(IT_0052 + IT_0053);
    const ccomplex_t IT_0055 = 1.4142135623731*e_em*IT_0054;
    const ccomplex_t IT_0056 = 0.5*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + (-2)*IT_0044 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0051*IT_0059;
    const ccomplex_t IT_0061 = -e_em*IT_0009*IT_0015*(IT_0005*IT_0048 + 
      -IT_0002*IT_0049);
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = IT_0006*IT_0048;
    const ccomplex_t IT_0064 = IT_0003*IT_0049;
    const ccomplex_t IT_0065 = IT_0009*(IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = 1.4142135623731*e_em*IT_0065;
    const ccomplex_t IT_0067 = (-0.5)*IT_0066;
    const ccomplex_t IT_0068 = pow(m_Z, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_12 + (-2)*IT_0044 + IT_0068 + 
      -reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0062*IT_0071;
    const ccomplex_t IT_0073 = (-2)*IT_0027 + (-2)*IT_0047 + -IT_0060 + 
      -IT_0072;
    const ccomplex_t IT_0074 = pow(m_Z, -2);
    const ccomplex_t IT_0075 = pow(s_14, 2);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = -IT_0044;
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = s_12*IT_0078;
    const ccomplex_t IT_0080 = 2*IT_0073;
    const ccomplex_t IT_0081 = s_24*IT_0068;
    const ccomplex_t IT_0082 = IT_0074*IT_0081;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = s_24 + IT_0083;
    const ccomplex_t IT_0085 = m_C_1*IT_0084;
    const ccomplex_t IT_0086 = 2*IT_0085;
    const ccomplex_t IT_0087 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0088 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0089 = IT_0002*IT_0088;
    const ccomplex_t IT_0090 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0091 = IT_0005*IT_0090;
    const ccomplex_t IT_0092 = IT_0009*(IT_0089 + -IT_0091);
    const ccomplex_t IT_0093 = 1.4142135623731*e_em*IT_0092;
    const ccomplex_t IT_0094 = 0.5*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0024*IT_0095;
    const ccomplex_t IT_0097 = IT_0087*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0100 = IT_0002*IT_0099;
    const ccomplex_t IT_0101 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0102 = IT_0005*IT_0101;
    const ccomplex_t IT_0103 = IT_0009*(IT_0100 + -IT_0102);
    const ccomplex_t IT_0104 = 1.4142135623731*e_em*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0108 = IT_0016*IT_0107;
    const ccomplex_t IT_0109 = IT_0019*IT_0107;
    const ccomplex_t IT_0110 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0111 = IT_0019*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + -IT_0109 + (-2)
      *IT_0111);
    const ccomplex_t IT_0113 = (-0.5)*IT_0112;
    const ccomplex_t IT_0114 = IT_0106*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0000 + IT_0044 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = -IT_0098 + -IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0120 = IT_0001*IT_0096*IT_0119;
    const ccomplex_t IT_0121 = IT_0025*IT_0087*IT_0119;
    const ccomplex_t IT_0122 = m_W*e_em;
    const ccomplex_t IT_0123 = cpow(IT_0014, -4);
    const ccomplex_t IT_0124 = IT_0002*IT_0008*IT_0049*IT_0123;
    const ccomplex_t IT_0125 = IT_0005*IT_0008*IT_0048*IT_0123;
    const ccomplex_t IT_0126 = cpow(IT_0008, -3);
    const ccomplex_t IT_0127 = IT_0002*IT_0049*IT_0126;
    const ccomplex_t IT_0128 = IT_0005*IT_0048*IT_0126;
    const ccomplex_t IT_0129 = cpow(IT_0014, -2);
    const ccomplex_t IT_0130 = IT_0002*IT_0009*IT_0049*IT_0129;
    const ccomplex_t IT_0131 = IT_0005*IT_0009*IT_0048*IT_0129;
    const ccomplex_t IT_0132 = cpow(IT_0008, -2);
    const ccomplex_t IT_0133 = cpow(IT_0129 + IT_0132, -1);
    const ccomplex_t IT_0134 = (IT_0124 + -IT_0125 + IT_0127 + -IT_0128 + 2
      *IT_0130 + (-2)*IT_0131)*IT_0133;
    const ccomplex_t IT_0135 = IT_0122*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = IT_0024*IT_0136;
    const ccomplex_t IT_0138 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0044 + -reg_prop, -1);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0142 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0143 = IT_0002*IT_0142;
    const ccomplex_t IT_0144 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0145 = IT_0005*IT_0144;
    const ccomplex_t IT_0146 = IT_0009*(IT_0143 + -IT_0145);
    const ccomplex_t IT_0147 = 1.4142135623731*e_em*IT_0146;
    const ccomplex_t IT_0148 = 0.5*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = IT_0042*IT_0149;
    const ccomplex_t IT_0151 = IT_0045*IT_0141*IT_0150;
    const ccomplex_t IT_0152 = V_u1*U_Wm2;
    const ccomplex_t IT_0153 = IT_0002*IT_0152;
    const ccomplex_t IT_0154 = U_d2*V_Wp1;
    const ccomplex_t IT_0155 = IT_0005*IT_0154;
    const ccomplex_t IT_0156 = IT_0009*(IT_0153 + -IT_0155);
    const ccomplex_t IT_0157 = 1.4142135623731*e_em*IT_0156;
    const ccomplex_t IT_0158 = 0.5*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0113*IT_0159;
    const ccomplex_t IT_0161 = IT_0115*IT_0141*IT_0160;
    const ccomplex_t IT_0162 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0163 = IT_0016*IT_0162;
    const ccomplex_t IT_0164 = IT_0019*IT_0162;
    const ccomplex_t IT_0165 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0166 = IT_0019*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + -IT_0164 + (-2)
      *IT_0166);
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = IT_0095*IT_0168;
    const ccomplex_t IT_0170 = IT_0001*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_C_1*IT_0171;
    const ccomplex_t IT_0173 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0174 = IT_0016*IT_0173;
    const ccomplex_t IT_0175 = IT_0019*IT_0173;
    const ccomplex_t IT_0176 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0177 = IT_0019*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + -IT_0175 + (-2)
      *IT_0177);
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = IT_0159*IT_0179;
    const ccomplex_t IT_0181 = IT_0115*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = m_C_1*IT_0182;
    const ccomplex_t IT_0184 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0185 = IT_0016*IT_0184;
    const ccomplex_t IT_0186 = IT_0019*IT_0184;
    const ccomplex_t IT_0187 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0188 = IT_0019*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0185 + -IT_0186 + (-2)
      *IT_0188);
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = IT_0149*IT_0190;
    const ccomplex_t IT_0192 = IT_0045*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_C_1*IT_0193;
    const ccomplex_t IT_0195 = IT_0013*IT_0168;
    const ccomplex_t IT_0196 = IT_0087*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = m_C_1*IT_0197;
    const ccomplex_t IT_0199 = -IT_0120 + -IT_0121 + IT_0140 + -IT_0151 + 
      -IT_0161 + -IT_0172 + -IT_0183 + -IT_0194 + -IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = 2*conj(IT_0199);
    const ccomplex_t IT_0202 = IT_0182 + IT_0197;
    const ccomplex_t IT_0203 = IT_0001*IT_0119*IT_0195;
    const ccomplex_t IT_0204 = IT_0106*IT_0179;
    const ccomplex_t IT_0205 = IT_0115*IT_0141*IT_0204;
    const ccomplex_t IT_0206 = IT_0035*IT_0190;
    const ccomplex_t IT_0207 = IT_0045*IT_0141*IT_0206;
    const ccomplex_t IT_0208 = m_C_1*IT_0027;
    const ccomplex_t IT_0209 = m_C_1*IT_0047;
    const ccomplex_t IT_0210 = IT_0136*IT_0168;
    const ccomplex_t IT_0211 = IT_0138*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = m_C_1*IT_0098;
    const ccomplex_t IT_0214 = IT_0087*IT_0119*IT_0169;
    const ccomplex_t IT_0215 = m_C_1*IT_0117;
    const ccomplex_t IT_0216 = IT_0203 + IT_0205 + IT_0207 + IT_0208 + IT_0209
       + -IT_0212 + IT_0213 + IT_0214 + IT_0215;
    const ccomplex_t IT_0217 = 2*IT_0216;
    const ccomplex_t IT_0218 = 2*conj(IT_0216);
    const ccomplex_t IT_0219 = (-0.25)*IT_0082;
    const ccomplex_t IT_0220 = s_24 + IT_0219;
    const ccomplex_t IT_0221 = m_C_1*IT_0220;
    const ccomplex_t IT_0222 = (-4)*IT_0221;
    const ccomplex_t IT_0223 = -IT_0171 + -IT_0193;
    const ccomplex_t IT_0224 = conj(IT_0027) + conj(IT_0047);
    const ccomplex_t IT_0225 = IT_0027 + IT_0047;
    const ccomplex_t IT_0226 = IT_0048*IT_0088;
    const ccomplex_t IT_0227 = IT_0049*IT_0090;
    const ccomplex_t IT_0228 = IT_0009*(IT_0226 + IT_0227);
    const ccomplex_t IT_0229 = 1.4142135623731*e_em*IT_0228;
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0057*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = IT_0051*IT_0232;
    const ccomplex_t IT_0234 = IT_0048*IT_0090;
    const ccomplex_t IT_0235 = IT_0049*IT_0088;
    const ccomplex_t IT_0236 = IT_0009*(IT_0234 + -IT_0235);
    const ccomplex_t IT_0237 = 1.4142135623731*e_em*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = IT_0069*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = IT_0062*IT_0240;
    const ccomplex_t IT_0242 = 2*IT_0171 + 2*IT_0193 + -IT_0233 + -IT_0241;
    const ccomplex_t IT_0243 = s_14*s_24*IT_0074;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = s_12 + IT_0244;
    const ccomplex_t IT_0246 = IT_0044*IT_0245;
    const ccomplex_t IT_0247 = pow(s_24, 2);
    const ccomplex_t IT_0248 = IT_0074*IT_0247;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0044 + IT_0249;
    const ccomplex_t IT_0251 = m_C_1*IT_0250;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = -IT_0060 + -IT_0072 + (-2)*IT_0182 + (-2)
      *IT_0197;
    const ccomplex_t IT_0254 = 2*conj(IT_0253);
    const ccomplex_t IT_0255 = 2*IT_0098 + 2*IT_0117 + -IT_0233 + -IT_0241;
    const ccomplex_t IT_0256 = 2*conj(IT_0255);
    const ccomplex_t IT_0257 = s_14*s_24;
    const ccomplex_t IT_0258 = pow(m_Z, 4);
    const ccomplex_t IT_0259 = s_12*IT_0258;
    const ccomplex_t IT_0260 = IT_0074*IT_0259;
    const ccomplex_t IT_0261 = (-0.25)*IT_0260;
    const ccomplex_t IT_0262 = IT_0257 + IT_0261;
    const ccomplex_t IT_0263 = (-8)*conj(IT_0202);
    const ccomplex_t IT_0264 = s_12*IT_0068;
    const ccomplex_t IT_0265 = (-0.5)*IT_0260;
    const ccomplex_t IT_0266 = (-2)*IT_0257;
    const ccomplex_t IT_0267 = IT_0264 + IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = (-4)*IT_0027;
    const ccomplex_t IT_0269 = s_12*IT_0250;
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = 2*IT_0253;
    const ccomplex_t IT_0272 = s_12*IT_0245;
    const ccomplex_t IT_0273 = (-2)*IT_0073;
    const ccomplex_t IT_0274 = 2*IT_0242;
    const ccomplex_t IT_0275 = IT_0044*IT_0078;
    const ccomplex_t IT_0276 = m_C_1*IT_0078;
    const ccomplex_t IT_0277 = (-2)*conj(IT_0073);
    const ccomplex_t IT_0278 = s_14*IT_0068;
    const ccomplex_t IT_0279 = IT_0074*IT_0278;
    const ccomplex_t IT_0280 = 0.5*IT_0279;
    const ccomplex_t IT_0281 = s_14 + IT_0280;
    const ccomplex_t IT_0282 = m_C_1*IT_0281;
    const ccomplex_t IT_0283 = 2*IT_0282;
    const ccomplex_t IT_0284 = (-2)*IT_0199;
    const ccomplex_t IT_0285 = (-2)*conj(IT_0199);
    const ccomplex_t IT_0286 = (-2)*IT_0216;
    const ccomplex_t IT_0287 = (-2)*conj(IT_0216);
    const ccomplex_t IT_0288 = (-0.25)*IT_0279;
    const ccomplex_t IT_0289 = s_14 + IT_0288;
    const ccomplex_t IT_0290 = m_C_1*IT_0289;
    const ccomplex_t IT_0291 = (-4)*IT_0290;
    const ccomplex_t IT_0292 = (-2)*conj(IT_0242);
    const ccomplex_t IT_0293 = s_12*s_14;
    const ccomplex_t IT_0294 = s_24*IT_0044;
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0293 + IT_0295;
    const ccomplex_t IT_0297 = 2*conj(IT_0073);
    const ccomplex_t IT_0298 = (-2)*IT_0242;
    const ccomplex_t IT_0299 = 0.5*IT_0080*IT_0224 + IT_0223*conj(IT_0242) +
       0.5*conj(IT_0202)*IT_0273 + 0.5*conj(IT_0223)*IT_0274 + 0.5*IT_0202
      *IT_0277 + 0.5*IT_0118*IT_0292 + 0.5*IT_0225*IT_0297 + 0.5*conj(IT_0118)
      *IT_0298;
    const ccomplex_t IT_0300 = m_C_1*IT_0245;
    const ccomplex_t IT_0301 = 2*IT_0255;
    const ccomplex_t IT_0302 = IT_0044*IT_0250;
    const ccomplex_t IT_0303 = -IT_0302;
    const ccomplex_t IT_0304 = (-2)*IT_0253;
    const ccomplex_t IT_0305 = (-2)*conj(IT_0253);
    const ccomplex_t IT_0306 = s_12*s_24;
    const ccomplex_t IT_0307 = s_14*IT_0044;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*IT_0202*IT_0254 + conj(IT_0223)*IT_0255 
      + IT_0223*conj(IT_0255) + (-0.5)*IT_0118*IT_0256 + (-0.5)*conj(IT_0202)
      *IT_0271 + (-0.5)*conj(IT_0118)*IT_0301 + (-0.5)*IT_0224*IT_0304 + (-0.5)
      *IT_0225*IT_0305;
    const ccomplex_t IT_0311 = IT_0044*IT_0068;
    const ccomplex_t IT_0312 = (-3)*IT_0311;
    const ccomplex_t IT_0313 = (-2)*IT_0223;
    const ccomplex_t IT_0314 = (-2)*conj(IT_0223);
    const ccomplex_t IT_0315 = (-2)*IT_0202;
    const ccomplex_t IT_0316 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0317 = 3*IT_0311;
    const ccomplex_t IT_0318 = (-2)*IT_0118;
    const ccomplex_t IT_0319 = (-2)*conj(IT_0118);
    const ccomplex_t IT_0320 = IT_0264 + IT_0266;
    const ccomplex_t IT_0321 = IT_0074*IT_0320;
    const ccomplex_t IT_0322 = (-0.5)*IT_0321;
    const ccomplex_t IT_0323 = 0.166666666666667*IT_0252;
    const ccomplex_t IT_0324 = conj(IT_0073)*IT_0079*IT_0080 + IT_0086*(conj
      (IT_0118)*IT_0200 + IT_0118*IT_0201 + conj(IT_0202)*IT_0217 + IT_0202
      *IT_0218) + IT_0222*(IT_0201*IT_0223 + IT_0200*conj(IT_0223) + IT_0217
      *IT_0224 + IT_0218*IT_0225) + (IT_0242*IT_0246 + IT_0216*IT_0252)*IT_0254 
      + IT_0199*IT_0252*IT_0256 + (-8)*IT_0262*(conj(IT_0118)*IT_0223 + IT_0118
      *conj(IT_0223) + IT_0202*IT_0224 + (-0.125)*IT_0225*IT_0263) + (-4)
      *IT_0267*(IT_0118*conj(IT_0118) + IT_0202*conj(IT_0202) + IT_0223*conj
      (IT_0223) + IT_0224*(IT_0047 + (-0.25)*IT_0268)) + conj(IT_0253)*(IT_0270
      *IT_0271 + IT_0272*IT_0273) + conj(IT_0242)*(IT_0246*IT_0271 + IT_0079
      *IT_0274 + IT_0273*IT_0275) + (IT_0253*IT_0272 + IT_0242*IT_0275 + IT_0199
      *IT_0276)*IT_0277 + IT_0283*(IT_0224*IT_0284 + IT_0225*IT_0285 + conj
      (IT_0223)*IT_0286 + IT_0223*IT_0287) + (conj(IT_0202)*IT_0284 + IT_0202
      *IT_0285 + conj(IT_0118)*IT_0286 + IT_0118*IT_0287)*IT_0291 + IT_0216
      *IT_0276*IT_0292 + 2*IT_0296*IT_0299 + IT_0300*(IT_0216*(IT_0256 + IT_0277
      ) + IT_0199*(IT_0254 + IT_0292) + conj(IT_0199)*(IT_0271 + IT_0298) + conj
      (IT_0216)*(IT_0273 + IT_0301)) + conj(IT_0255)*(IT_0080*IT_0246 + IT_0272
      *IT_0298 + IT_0270*IT_0301 + IT_0303*IT_0304) + IT_0255*(IT_0272*IT_0292 +
       IT_0246*IT_0297 + IT_0303*IT_0305) + (-2)*IT_0309*IT_0310 + IT_0312*
      (IT_0224*IT_0313 + IT_0225*IT_0314 + conj(IT_0118)*IT_0315 + IT_0118
      *IT_0316) + IT_0317*(conj(IT_0223)*IT_0315 + IT_0223*IT_0316 + IT_0224
      *IT_0318 + IT_0225*IT_0319) + 6*conj(IT_0216)*(IT_0044*IT_0199 +
       0.166666666666667*IT_0276*IT_0298 + 0.666666666666667*IT_0216*(s_12 +
       IT_0322) + IT_0271*IT_0323) + 6*conj(IT_0199)*(IT_0044*IT_0216 +
       0.166666666666667*IT_0273*IT_0276 + 0.666666666666667*IT_0199*(s_12 +
       IT_0322) + IT_0301*IT_0323);
    return create_ccomplex_return(IT_0324);
}

