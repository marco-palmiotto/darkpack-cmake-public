#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_C_2_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_C_2_to_Z_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0000 = m_W*e_em;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -4);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = IT_0001*IT_0003*IT_0004*IT_0005;
    const ccomplex_t IT_0007 = cos(beta);
    const ccomplex_t IT_0008 = sin(alpha);
    const ccomplex_t IT_0009 = IT_0003*IT_0005*IT_0007*IT_0008;
    const ccomplex_t IT_0010 = cpow(IT_0005, -3);
    const ccomplex_t IT_0011 = IT_0001*IT_0004*IT_0010;
    const ccomplex_t IT_0012 = IT_0007*IT_0008*IT_0010;
    const ccomplex_t IT_0013 = cpow(IT_0002, -2);
    const ccomplex_t IT_0014 = cpow(IT_0005, -1);
    const ccomplex_t IT_0015 = IT_0001*IT_0004*IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0008*IT_0013*IT_0014;
    const ccomplex_t IT_0017 = cpow(IT_0005, -2);
    const ccomplex_t IT_0018 = cpow(IT_0013 + IT_0017, -1);
    const ccomplex_t IT_0019 = (IT_0006 + -IT_0009 + IT_0011 + -IT_0012 + 2
      *IT_0015 + (-2)*IT_0016)*IT_0018;
    const ccomplex_t IT_0020 = IT_0000*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = V_u2*U_Wm2;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = U_d2*V_Wp2;
    const ccomplex_t IT_0025 = IT_0008*IT_0024;
    const ccomplex_t IT_0026 = IT_0014*(IT_0023 + -IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = pow(m_C_2, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0031 + -reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0001*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0036 = IT_0003*IT_0004*IT_0005*IT_0008;
    const ccomplex_t IT_0037 = IT_0001*IT_0007*IT_0010;
    const ccomplex_t IT_0038 = IT_0004*IT_0008*IT_0010;
    const ccomplex_t IT_0039 = IT_0001*IT_0007*IT_0013*IT_0014;
    const ccomplex_t IT_0040 = IT_0004*IT_0008*IT_0013*IT_0014;
    const ccomplex_t IT_0041 = IT_0018*(IT_0035 + IT_0036 + IT_0037 + IT_0038 
      + 2*IT_0039 + 2*IT_0040);
    const ccomplex_t IT_0042 = IT_0000*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0001*IT_0024;
    const ccomplex_t IT_0045 = IT_0008*IT_0022;
    const ccomplex_t IT_0046 = IT_0014*(IT_0044 + IT_0045);
    const ccomplex_t IT_0047 = 1.4142135623731*e_em*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = IT_0043*IT_0049;
    const ccomplex_t IT_0051 = cpow((-2)*s_12 + (-2)*IT_0031 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = IT_0034 + IT_0053;
    const ccomplex_t IT_0055 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0056 = IT_0001*IT_0055;
    const ccomplex_t IT_0057 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0058 = IT_0008*IT_0057;
    const ccomplex_t IT_0059 = IT_0014*(IT_0056 + -IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0021*IT_0062;
    const ccomplex_t IT_0064 = IT_0032*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0001*IT_0057;
    const ccomplex_t IT_0067 = IT_0008*IT_0055;
    const ccomplex_t IT_0068 = IT_0014*(IT_0066 + IT_0067);
    const ccomplex_t IT_0069 = 1.4142135623731*e_em*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0070;
    const ccomplex_t IT_0072 = IT_0043*IT_0071;
    const ccomplex_t IT_0073 = IT_0051*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0065 + IT_0074;
    const ccomplex_t IT_0076 = pow(m_Z, -4);
    const ccomplex_t IT_0077 = pow(s_34, 2);
    const ccomplex_t IT_0078 = IT_0076*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = 1 + IT_0079;
    const ccomplex_t IT_0081 = s_12*IT_0080;
    const ccomplex_t IT_0082 = cpow(IT_0002, -1);
    const ccomplex_t IT_0083 = IT_0005*IT_0082;
    const ccomplex_t IT_0084 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = IT_0002*IT_0014;
    const ccomplex_t IT_0087 = IT_0084*IT_0086;
    const ccomplex_t IT_0088 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0089 = IT_0086*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0085 + -IT_0087 + (-2)
      *IT_0089);
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = cpow(IT_0091, 2);
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*IT_0092;
    const ccomplex_t IT_0094 = pow(m_Z, 2);
    const ccomplex_t IT_0095 = cpow((-2)*s_23 + IT_0094 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0093*IT_0095;
    const ccomplex_t IT_0097 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0098 = IT_0083*IT_0097;
    const ccomplex_t IT_0099 = IT_0086*IT_0097;
    const ccomplex_t IT_0100 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0101 = IT_0086*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + -IT_0099 + (-2)
      *IT_0101);
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0105 = IT_0083*IT_0104;
    const ccomplex_t IT_0106 = IT_0086*IT_0104;
    const ccomplex_t IT_0107 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0108 = IT_0086*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + -IT_0106 + (-2)
      *IT_0108);
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = IT_0103*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_23 + IT_0031 + IT_0094 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0111*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = 2*IT_0096 + 2*IT_0114;
    const ccomplex_t IT_0116 = s_14*s_23*s_34*IT_0076;
    const ccomplex_t IT_0117 = pow(m_Z, -2);
    const ccomplex_t IT_0118 = s_13*s_23*IT_0117;
    const ccomplex_t IT_0119 = s_14*s_24*IT_0117;
    const ccomplex_t IT_0120 = -IT_0118 + -IT_0119;
    const ccomplex_t IT_0121 = s_12 + IT_0116 + IT_0120;
    const ccomplex_t IT_0122 = m_C_2*IT_0121;
    const ccomplex_t IT_0123 = 2*conj(IT_0054);
    const ccomplex_t IT_0124 = s_13*s_14*s_34*IT_0076;
    const ccomplex_t IT_0125 = IT_0031 + IT_0124;
    const ccomplex_t IT_0126 = pow(s_13, 2);
    const ccomplex_t IT_0127 = IT_0117*IT_0126;
    const ccomplex_t IT_0128 = pow(s_14, 2);
    const ccomplex_t IT_0129 = IT_0117*IT_0128;
    const ccomplex_t IT_0130 = -IT_0127 + -IT_0129;
    const ccomplex_t IT_0131 = IT_0125 + IT_0130;
    const ccomplex_t IT_0132 = m_C_2*IT_0131;
    const ccomplex_t IT_0133 = -IT_0132;
    const ccomplex_t IT_0134 = 2*conj(IT_0075);
    const ccomplex_t IT_0135 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0136 = IT_0083*IT_0135;
    const ccomplex_t IT_0137 = IT_0086*IT_0135;
    const ccomplex_t IT_0138 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0139 = IT_0086*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + -IT_0137 + (-2)
      *IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0143 = IT_0083*IT_0142;
    const ccomplex_t IT_0144 = IT_0086*IT_0142;
    const ccomplex_t IT_0145 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0146 = IT_0086*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*(IT_0143 + -IT_0144 + (-2)
      *IT_0146);
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0141*IT_0148;
    const ccomplex_t IT_0150 = IT_0112*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0153 = IT_0083*IT_0152;
    const ccomplex_t IT_0154 = IT_0086*IT_0152;
    const ccomplex_t IT_0155 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0156 = IT_0086*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + -IT_0154 + (-2)
      *IT_0156);
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = cpow(IT_0158, 2);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0095*IT_0160;
    const ccomplex_t IT_0162 = 2*IT_0151 + 2*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0164 = IT_0110*IT_0148;
    const ccomplex_t IT_0165 = cpow((-2)*s_13 + IT_0031 + IT_0094 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0166 = IT_0163*IT_0164*IT_0165;
    const ccomplex_t IT_0167 = IT_0149*IT_0165;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = m_C_2*IT_0168;
    const ccomplex_t IT_0170 = cpow((-2)*s_13 + IT_0094 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0171 = IT_0091*IT_0158;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0173 = IT_0170*IT_0171*IT_0172;
    const ccomplex_t IT_0174 = IT_0160*IT_0170;
    const ccomplex_t IT_0175 = m_C_2*IT_0174;
    const ccomplex_t IT_0176 = -IT_0166 + -IT_0169 + -IT_0173 + -IT_0175;
    const ccomplex_t IT_0177 = s_34*IT_0076;
    const ccomplex_t IT_0178 = s_14*s_23;
    const ccomplex_t IT_0179 = s_13*s_24;
    const ccomplex_t IT_0180 = s_12*s_34;
    const ccomplex_t IT_0181 = -IT_0179 + -IT_0180;
    const ccomplex_t IT_0182 = IT_0178 + IT_0181;
    const ccomplex_t IT_0183 = IT_0177*IT_0182;
    const ccomplex_t IT_0184 = (-0.25)*IT_0183;
    const ccomplex_t IT_0185 = s_12*IT_0094;
    const ccomplex_t IT_0186 = IT_0117*IT_0185;
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = s_12 + IT_0184 + IT_0187;
    const ccomplex_t IT_0189 = 8*conj(IT_0054);
    const ccomplex_t IT_0190 = pow(m_Z, 4);
    const ccomplex_t IT_0191 = s_12*IT_0190;
    const ccomplex_t IT_0192 = s_13*s_23*IT_0094;
    const ccomplex_t IT_0193 = s_14*s_24*IT_0094;
    const ccomplex_t IT_0194 = s_14*s_23*s_34;
    const ccomplex_t IT_0195 = s_13*s_24*s_34;
    const ccomplex_t IT_0196 = s_12*IT_0077;
    const ccomplex_t IT_0197 = (-2)*IT_0192 + (-2)*IT_0193 + 2*IT_0194 + 2
      *IT_0195 + (-2)*IT_0196;
    const ccomplex_t IT_0198 = IT_0191 + IT_0197;
    const ccomplex_t IT_0199 = IT_0076*IT_0198;
    const ccomplex_t IT_0200 = (-0.25)*IT_0199;
    const ccomplex_t IT_0201 = -IT_0119;
    const ccomplex_t IT_0202 = -IT_0118;
    const ccomplex_t IT_0203 = s_12 + IT_0200 + IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = IT_0112*IT_0163*IT_0164;
    const ccomplex_t IT_0205 = IT_0095*IT_0171*IT_0172;
    const ccomplex_t IT_0206 = m_C_2*IT_0096;
    const ccomplex_t IT_0207 = -IT_0204 + -IT_0205 + -IT_0206;
    const ccomplex_t IT_0208 = 8*conj(IT_0207);
    const ccomplex_t IT_0209 = IT_0031*IT_0080;
    const ccomplex_t IT_0210 = IT_0103*IT_0141;
    const ccomplex_t IT_0211 = IT_0163*IT_0165*IT_0210;
    const ccomplex_t IT_0212 = IT_0111*IT_0165;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = m_C_2*IT_0213;
    const ccomplex_t IT_0215 = IT_0093*IT_0170;
    const ccomplex_t IT_0216 = m_C_2*IT_0215;
    const ccomplex_t IT_0217 = -IT_0173 + -IT_0211 + -IT_0214 + -IT_0216;
    const ccomplex_t IT_0218 = IT_0112*IT_0163*IT_0210;
    const ccomplex_t IT_0219 = m_C_2*IT_0151;
    const ccomplex_t IT_0220 = m_C_2*IT_0161;
    const ccomplex_t IT_0221 = -IT_0205 + -IT_0218 + -IT_0219 + -IT_0220;
    const ccomplex_t IT_0222 = (-4)*IT_0054;
    const ccomplex_t IT_0223 = (-4)*conj(IT_0054);
    const ccomplex_t IT_0224 = m_C_2*IT_0114;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = (-4)*IT_0075;
    const ccomplex_t IT_0227 = IT_0207 + IT_0225;
    const ccomplex_t IT_0228 = (-4)*conj(IT_0075);
    const ccomplex_t IT_0229 = (-2)*IT_0213 + (-2)*IT_0215;
    const ccomplex_t IT_0230 = s_13*s_23;
    const ccomplex_t IT_0231 = (-2)*IT_0230;
    const ccomplex_t IT_0232 = IT_0185 + IT_0231;
    const ccomplex_t IT_0233 = IT_0117*IT_0232;
    const ccomplex_t IT_0234 = (s_12 + -IT_0119)*(s_12 + (-0.5)*IT_0233);
    const ccomplex_t IT_0235 = (-2)*IT_0234;
    const ccomplex_t IT_0236 = IT_0031*IT_0094;
    const ccomplex_t IT_0237 = s_13*s_14*s_34*IT_0117;
    const ccomplex_t IT_0238 = 0.5*IT_0128 + (-0.5)*IT_0236 + -IT_0237;
    const ccomplex_t IT_0239 = IT_0126 + IT_0238;
    const ccomplex_t IT_0240 = IT_0117*IT_0239;
    const ccomplex_t IT_0241 = -IT_0129;
    const ccomplex_t IT_0242 = IT_0031 + IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = m_C_2*IT_0242;
    const ccomplex_t IT_0244 = IT_0031 + IT_0241;
    const ccomplex_t IT_0245 = IT_0031*IT_0244;
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = (-2)*IT_0168 + (-2)*IT_0174;
    const ccomplex_t IT_0248 = s_12 + IT_0201;
    const ccomplex_t IT_0249 = IT_0031*IT_0248;
    const ccomplex_t IT_0250 = s_14*s_24;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0185 + IT_0251;
    const ccomplex_t IT_0253 = IT_0117*IT_0252;
    const ccomplex_t IT_0254 = -IT_0119 + (-0.25)*IT_0253;
    const ccomplex_t IT_0255 = s_12 + IT_0254;
    const ccomplex_t IT_0256 = m_C_2*IT_0255;
    const ccomplex_t IT_0257 = (-4)*IT_0256;
    const ccomplex_t IT_0258 = s_14*s_23*s_34*IT_0117;
    const ccomplex_t IT_0259 = (-2)*IT_0230 + -IT_0250 + 2*IT_0258;
    const ccomplex_t IT_0260 = IT_0185 + IT_0259;
    const ccomplex_t IT_0261 = IT_0117*IT_0260;
    const ccomplex_t IT_0262 = -IT_0119 + (-0.5)*IT_0261;
    const ccomplex_t IT_0263 = s_12 + IT_0262;
    const ccomplex_t IT_0264 = m_C_2*IT_0263;
    const ccomplex_t IT_0265 = (-2)*IT_0264;
    const ccomplex_t IT_0266 = 2*IT_0054;
    const ccomplex_t IT_0267 = 2*IT_0075;
    const ccomplex_t IT_0268 = (IT_0031 + -IT_0129)*(s_12 + (-0.5)*IT_0233);
    const ccomplex_t IT_0269 = 2*IT_0268;
    const ccomplex_t IT_0270 = (-2)*IT_0115;
    const ccomplex_t IT_0271 = IT_0229*IT_0235 + 2*IT_0217*IT_0243 + 3*IT_0162
      *IT_0246 + (-3)*IT_0247*IT_0249 + -IT_0176*IT_0257 + IT_0227*IT_0265 + 0.5
      *IT_0122*IT_0266 + 0.5*IT_0133*IT_0267 + 0.5*IT_0269*IT_0270;
    const ccomplex_t IT_0272 = 2*conj(IT_0162);
    const ccomplex_t IT_0273 = IT_0178 + IT_0180;
    const ccomplex_t IT_0274 = -IT_0179;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = IT_0177*IT_0275;
    const ccomplex_t IT_0277 = 0.25*IT_0276;
    const ccomplex_t IT_0278 = s_12 + IT_0187 + IT_0277;
    const ccomplex_t IT_0279 = 8*conj(IT_0075);
    const ccomplex_t IT_0280 = (-2)*IT_0077;
    const ccomplex_t IT_0281 = IT_0190 + IT_0280;
    const ccomplex_t IT_0282 = IT_0076*IT_0281;
    const ccomplex_t IT_0283 = 0.25*IT_0282;
    const ccomplex_t IT_0284 = 1 + IT_0283;
    const ccomplex_t IT_0285 = IT_0031*IT_0284;
    const ccomplex_t IT_0286 = 8*conj(IT_0176);
    const ccomplex_t IT_0287 = -IT_0236;
    const ccomplex_t IT_0288 = IT_0128 + IT_0287;
    const ccomplex_t IT_0289 = IT_0117*IT_0288;
    const ccomplex_t IT_0290 = -IT_0031;
    const ccomplex_t IT_0291 = IT_0129 + IT_0290;
    const ccomplex_t IT_0292 = (-4)*IT_0291;
    const ccomplex_t IT_0293 = IT_0289 + IT_0292;
    const ccomplex_t IT_0294 = m_C_2*IT_0293;
    const ccomplex_t IT_0295 = (-2)*conj(IT_0115);
    const ccomplex_t IT_0296 = pow(s_24, 2);
    const ccomplex_t IT_0297 = IT_0287 + IT_0296;
    const ccomplex_t IT_0298 = IT_0117*IT_0297;
    const ccomplex_t IT_0299 = IT_0117*IT_0296;
    const ccomplex_t IT_0300 = IT_0290 + IT_0299;
    const ccomplex_t IT_0301 = (-4)*IT_0300;
    const ccomplex_t IT_0302 = IT_0298 + IT_0301;
    const ccomplex_t IT_0303 = m_C_2*IT_0302;
    const ccomplex_t IT_0304 = 2*conj(IT_0176);
    const ccomplex_t IT_0305 = s_23*s_24*s_34*IT_0076;
    const ccomplex_t IT_0306 = IT_0031 + IT_0305;
    const ccomplex_t IT_0307 = pow(s_23, 2);
    const ccomplex_t IT_0308 = IT_0117*IT_0307;
    const ccomplex_t IT_0309 = -IT_0299 + -IT_0308;
    const ccomplex_t IT_0310 = IT_0306 + IT_0309;
    const ccomplex_t IT_0311 = m_C_2*IT_0310;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = (-2)*conj(IT_0054);
    const ccomplex_t IT_0314 = s_13*s_24*s_34*IT_0076;
    const ccomplex_t IT_0315 = s_12 + IT_0120 + IT_0314;
    const ccomplex_t IT_0316 = m_C_2*IT_0315;
    const ccomplex_t IT_0317 = (-2)*conj(IT_0075);
    const ccomplex_t IT_0318 = -IT_0114;
    const ccomplex_t IT_0319 = -IT_0096;
    const ccomplex_t IT_0320 = conj(IT_0318) + conj(IT_0319);
    const ccomplex_t IT_0321 = s_24*IT_0094;
    const ccomplex_t IT_0322 = s_23*s_34;
    const ccomplex_t IT_0323 = (-2)*IT_0322;
    const ccomplex_t IT_0324 = IT_0321 + IT_0323;
    const ccomplex_t IT_0325 = IT_0117*IT_0324;
    const ccomplex_t IT_0326 = s_24*IT_0190;
    const ccomplex_t IT_0327 = s_23*s_34*IT_0094;
    const ccomplex_t IT_0328 = (-2)*IT_0327;
    const ccomplex_t IT_0329 = IT_0326 + IT_0328;
    const ccomplex_t IT_0330 = IT_0076*IT_0329;
    const ccomplex_t IT_0331 = (-0.666666666666667)*IT_0325 +
       0.166666666666667*IT_0330;
    const ccomplex_t IT_0332 = s_24 + IT_0331;
    const ccomplex_t IT_0333 = m_C_2*IT_0332;
    const ccomplex_t IT_0334 = 6*IT_0333;
    const ccomplex_t IT_0335 = (-2)*IT_0225;
    const ccomplex_t IT_0336 = (-2)*IT_0207;
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = s_23*IT_0094;
    const ccomplex_t IT_0339 = IT_0177*IT_0338;
    const ccomplex_t IT_0340 = IT_0117*IT_0321;
    const ccomplex_t IT_0341 = (-0.5)*IT_0325 + (-0.5)*IT_0339 + 0.5*IT_0340;
    const ccomplex_t IT_0342 = s_24 + IT_0341;
    const ccomplex_t IT_0343 = m_C_2*IT_0342;
    const ccomplex_t IT_0344 = (-2)*IT_0343;
    const ccomplex_t IT_0345 = m_C_2*IT_0117;
    const ccomplex_t IT_0346 = IT_0117*IT_0326;
    const ccomplex_t IT_0347 = (-2)*IT_0322 + (-0.5)*IT_0346;
    const ccomplex_t IT_0348 = IT_0321 + IT_0347;
    const ccomplex_t IT_0349 = IT_0345*IT_0348;
    const ccomplex_t IT_0350 = (-2)*IT_0349;
    const ccomplex_t IT_0351 = 2*IT_0176;
    const ccomplex_t IT_0352 = s_13*s_34*IT_0117;
    const ccomplex_t IT_0353 = s_14*IT_0190;
    const ccomplex_t IT_0354 = s_13*s_34*IT_0094;
    const ccomplex_t IT_0355 = (-2)*IT_0354;
    const ccomplex_t IT_0356 = IT_0353 + IT_0355;
    const ccomplex_t IT_0357 = IT_0076*IT_0356;
    const ccomplex_t IT_0358 = (-0.666666666666667)*IT_0352 + (
      -0.166666666666667)*IT_0357;
    const ccomplex_t IT_0359 = s_14 + IT_0358;
    const ccomplex_t IT_0360 = m_C_2*IT_0359;
    const ccomplex_t IT_0361 = (-6)*IT_0360;
    const ccomplex_t IT_0362 = 2*IT_0217;
    const ccomplex_t IT_0363 = s_14*IT_0031;
    const ccomplex_t IT_0364 = s_14*IT_0307;
    const ccomplex_t IT_0365 = s_12*s_24*IT_0094;
    const ccomplex_t IT_0366 = s_12*s_23*s_34;
    const ccomplex_t IT_0367 = s_14*IT_0031*IT_0094;
    const ccomplex_t IT_0368 = 0.5*IT_0365 + -IT_0366 + (-0.5)*IT_0367;
    const ccomplex_t IT_0369 = IT_0364 + IT_0368;
    const ccomplex_t IT_0370 = IT_0117*IT_0369;
    const ccomplex_t IT_0371 = s_12*s_24;
    const ccomplex_t IT_0372 = -IT_0371;
    const ccomplex_t IT_0373 = IT_0363 + IT_0370 + IT_0372;
    const ccomplex_t IT_0374 = 4*IT_0229;
    const ccomplex_t IT_0375 = pow(m_Z, 6);
    const ccomplex_t IT_0376 = s_12*IT_0375;
    const ccomplex_t IT_0377 = s_13*s_23*IT_0190;
    const ccomplex_t IT_0378 = (-2)*IT_0377;
    const ccomplex_t IT_0379 = IT_0376 + IT_0378;
    const ccomplex_t IT_0380 = IT_0076*IT_0379;
    const ccomplex_t IT_0381 = 0.5*IT_0380;
    const ccomplex_t IT_0382 = (-4)*IT_0250;
    const ccomplex_t IT_0383 = IT_0185 + IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = (-2)*IT_0192 + (-2)*IT_0193 + 4*IT_0194;
    const ccomplex_t IT_0385 = IT_0191 + IT_0384;
    const ccomplex_t IT_0386 = IT_0117*IT_0385;
    const ccomplex_t IT_0387 = -IT_0386;
    const ccomplex_t IT_0388 = IT_0383 + IT_0387;
    const ccomplex_t IT_0389 = (-4)*IT_0319;
    const ccomplex_t IT_0390 = s_12*s_14*IT_0094;
    const ccomplex_t IT_0391 = s_13*s_14*s_23;
    const ccomplex_t IT_0392 = s_24*IT_0031*IT_0094;
    const ccomplex_t IT_0393 = s_23*s_34*IT_0031;
    const ccomplex_t IT_0394 = (-2)*IT_0391 + -IT_0392 + 2*IT_0393;
    const ccomplex_t IT_0395 = IT_0390 + IT_0394;
    const ccomplex_t IT_0396 = IT_0117*IT_0395;
    const ccomplex_t IT_0397 = s_24*IT_0031;
    const ccomplex_t IT_0398 = 2*IT_0397;
    const ccomplex_t IT_0399 = s_12*s_14;
    const ccomplex_t IT_0400 = (-2)*IT_0399;
    const ccomplex_t IT_0401 = IT_0396 + IT_0398 + IT_0400;
    const ccomplex_t IT_0402 = 9*IT_0236;
    const ccomplex_t IT_0403 = -IT_0151 + -IT_0161;
    const ccomplex_t IT_0404 = (-2)*IT_0403;
    const ccomplex_t IT_0405 = IT_0031*IT_0117;
    const ccomplex_t IT_0406 = (-0.25)*IT_0190;
    const ccomplex_t IT_0407 = IT_0077 + IT_0406;
    const ccomplex_t IT_0408 = IT_0405*IT_0407;
    const ccomplex_t IT_0409 = (-4)*IT_0408;
    const ccomplex_t IT_0410 = IT_0168 + IT_0174;
    const ccomplex_t IT_0411 = (-2)*IT_0410;
    const ccomplex_t IT_0412 = s_13*IT_0094;
    const ccomplex_t IT_0413 = IT_0177*IT_0412;
    const ccomplex_t IT_0414 = s_14*IT_0094;
    const ccomplex_t IT_0415 = IT_0117*IT_0414;
    const ccomplex_t IT_0416 = s_13*s_34;
    const ccomplex_t IT_0417 = (-2)*IT_0416;
    const ccomplex_t IT_0418 = IT_0414 + IT_0417;
    const ccomplex_t IT_0419 = IT_0117*IT_0418;
    const ccomplex_t IT_0420 = (-0.5)*IT_0413 + 0.5*IT_0415 + (-0.5)*IT_0419;
    const ccomplex_t IT_0421 = s_14 + IT_0420;
    const ccomplex_t IT_0422 = m_C_2*IT_0421;
    const ccomplex_t IT_0423 = (-2)*IT_0422;
    const ccomplex_t IT_0424 = (-2)*IT_0075;
    const ccomplex_t IT_0425 = IT_0076*IT_0353;
    const ccomplex_t IT_0426 = 0.666666666666667*IT_0415 + (-0.166666666666667
      )*IT_0425;
    const ccomplex_t IT_0427 = s_14 + IT_0426;
    const ccomplex_t IT_0428 = m_C_2*IT_0427;
    const ccomplex_t IT_0429 = (-6)*IT_0428;
    const ccomplex_t IT_0430 = (-2)*IT_0221;
    const ccomplex_t IT_0431 = IT_0334*IT_0337 + IT_0266*IT_0344 + IT_0350
      *IT_0351 + IT_0361*IT_0362 + IT_0373*IT_0374 + IT_0388*IT_0389 + IT_0270
      *IT_0401 + IT_0402*IT_0404 + IT_0409*IT_0411 + IT_0423*IT_0424 + IT_0429
      *IT_0430;
    const ccomplex_t IT_0432 = IT_0213 + IT_0215;
    const ccomplex_t IT_0433 = 0.25*IT_0380;
    const ccomplex_t IT_0434 = IT_0185 + IT_0433;
    const ccomplex_t IT_0435 = (-2)*conj(IT_0225);
    const ccomplex_t IT_0436 = (-2)*conj(IT_0207);
    const ccomplex_t IT_0437 = IT_0435 + IT_0436;
    const ccomplex_t IT_0438 = 2*conj(IT_0217);
    const ccomplex_t IT_0439 = 4*conj(IT_0229);
    const ccomplex_t IT_0440 = (-2)*conj(IT_0403);
    const ccomplex_t IT_0441 = (-2)*conj(IT_0410);
    const ccomplex_t IT_0442 = (-2)*conj(IT_0221);
    const ccomplex_t IT_0443 = (-0.125)*IT_0123*IT_0344 + (-0.125)*IT_0304
      *IT_0350 + (-0.125)*IT_0295*IT_0401 + (-0.125)*IT_0317*IT_0423 + conj
      (IT_0432)*IT_0434 + (-0.125)*IT_0334*IT_0437 + (-0.125)*IT_0361*IT_0438 + 
      (-0.125)*IT_0373*IT_0439 + (-0.125)*IT_0402*IT_0440 + (-0.125)*IT_0409
      *IT_0441 + (-0.125)*IT_0429*IT_0442;
    const ccomplex_t IT_0444 = (-2)*IT_0162;
    const ccomplex_t IT_0445 = 4*IT_0247;
    const ccomplex_t IT_0446 = (-2)*IT_0054;
    const ccomplex_t IT_0447 = 2*conj(IT_0247);
    const ccomplex_t IT_0448 = (-2)*conj(IT_0162);
    const ccomplex_t IT_0449 = s_13*s_24*s_34*IT_0117;
    const ccomplex_t IT_0450 = (-2)*IT_0230 + -IT_0250 + 2*IT_0449;
    const ccomplex_t IT_0451 = IT_0185 + IT_0450;
    const ccomplex_t IT_0452 = IT_0117*IT_0451;
    const ccomplex_t IT_0453 = -IT_0119 + (-0.5)*IT_0452;
    const ccomplex_t IT_0454 = s_12 + IT_0453;
    const ccomplex_t IT_0455 = m_C_2*IT_0454;
    const ccomplex_t IT_0456 = (-2)*IT_0455;
    const ccomplex_t IT_0457 = -IT_0299;
    const ccomplex_t IT_0458 = IT_0031 + IT_0457;
    const ccomplex_t IT_0459 = IT_0031*IT_0458;
    const ccomplex_t IT_0460 = -IT_0459;
    const ccomplex_t IT_0461 = s_23*s_24*s_34*IT_0117;
    const ccomplex_t IT_0462 = (-0.5)*IT_0236 + 0.5*IT_0296 + -IT_0461;
    const ccomplex_t IT_0463 = IT_0307 + IT_0462;
    const ccomplex_t IT_0464 = IT_0117*IT_0463;
    const ccomplex_t IT_0465 = IT_0031 + IT_0457 + IT_0464;
    const ccomplex_t IT_0466 = m_C_2*IT_0465;
    const ccomplex_t IT_0467 = 2*IT_0162;
    const ccomplex_t IT_0468 = (s_12 + (-0.5)*IT_0233)*(IT_0031 + -IT_0299);
    const ccomplex_t IT_0469 = 2*IT_0468;
    const ccomplex_t IT_0470 = (-2)*IT_0247;
    const ccomplex_t IT_0471 = IT_0115*IT_0249 + (-0.166666666666667)*IT_0312
      *IT_0424 + (-0.166666666666667)*IT_0316*IT_0446 + 0.333333333333333
      *IT_0176*IT_0456 + -IT_0229*IT_0460 + 0.666666666666667*IT_0221*IT_0466 + 
      (-0.166666666666667)*IT_0235*IT_0467 + (-0.166666666666667)*IT_0469*IT_0470;
    const ccomplex_t IT_0472 = 1.125*IT_0031;
    const ccomplex_t IT_0473 = 8*IT_0054;
    const ccomplex_t IT_0474 = 8*IT_0207;
    const ccomplex_t IT_0475 = (-2)*IT_0250;
    const ccomplex_t IT_0476 = IT_0185 + IT_0475;
    const ccomplex_t IT_0477 = IT_0117*IT_0476;
    const ccomplex_t IT_0478 = (-0.5)*IT_0477;
    const ccomplex_t IT_0479 = s_12 + IT_0478;
    const ccomplex_t IT_0480 = (-0.5)*IT_0233;
    const ccomplex_t IT_0481 = (-2)*IT_0192 + (-2)*IT_0193 + 4*IT_0195;
    const ccomplex_t IT_0482 = IT_0191 + IT_0481;
    const ccomplex_t IT_0483 = IT_0076*IT_0482;
    const ccomplex_t IT_0484 = 0.25*IT_0483;
    const ccomplex_t IT_0485 = IT_0479 + IT_0480 + IT_0484;
    const ccomplex_t IT_0486 = 8*IT_0176;
    const ccomplex_t IT_0487 = (-2)*IT_0229;
    const ccomplex_t IT_0488 = IT_0115*IT_0235 + (-3)*IT_0162*IT_0249 +
       IT_0221*IT_0257 + 0.5*IT_0303*IT_0351 + 0.5*IT_0316*IT_0424 + 0.5*IT_0312
      *IT_0446 + -IT_0217*IT_0456 + 3*IT_0247*IT_0460 + (-2)*IT_0227*IT_0466 +
       0.5*IT_0469*IT_0487;
    const ccomplex_t IT_0489 = 0.125*IT_0474;
    const ccomplex_t IT_0490 = IT_0225 + IT_0489;
    const ccomplex_t IT_0491 = 2*IT_0247;
    const ccomplex_t IT_0492 = 8*IT_0075;
    const ccomplex_t IT_0493 = 8*IT_0217;
    const ccomplex_t IT_0494 = 8*IT_0221;
    const ccomplex_t IT_0495 = 4*IT_0115*IT_0243 + IT_0257*IT_0444 + (-16)
      *IT_0176*IT_0472 + IT_0456*IT_0487 + 8*IT_0285*IT_0490 + IT_0303*IT_0491 +
       IT_0188*IT_0492 + IT_0485*IT_0493 + IT_0203*IT_0494;
    const ccomplex_t IT_0496 = IT_0076*IT_0385;
    const ccomplex_t IT_0497 = 0.25*IT_0496;
    const ccomplex_t IT_0498 = IT_0479 + IT_0480 + IT_0497;
    const ccomplex_t IT_0499 = 2*IT_0229*IT_0257 + IT_0270*IT_0294 + (-4)
      *IT_0247*IT_0466 + IT_0265*IT_0467 + ((-16)*IT_0207 + (-16)*IT_0225)
      *IT_0472 + IT_0285*IT_0486 + IT_0278*IT_0492 + IT_0203*IT_0493 + IT_0494
      *IT_0498;
    const ccomplex_t IT_0500 = s_24*IT_0126;
    const ccomplex_t IT_0501 = s_12*s_13*s_34;
    const ccomplex_t IT_0502 = -IT_0392 + 2*IT_0500 + (-2)*IT_0501;
    const ccomplex_t IT_0503 = IT_0390 + IT_0502;
    const ccomplex_t IT_0504 = IT_0117*IT_0503;
    const ccomplex_t IT_0505 = (-0.5)*IT_0504;
    const ccomplex_t IT_0506 = -IT_0397;
    const ccomplex_t IT_0507 = IT_0399 + IT_0505 + IT_0506;
    const ccomplex_t IT_0508 = (-0.25)*IT_0423;
    const ccomplex_t IT_0509 = (-0.25)*IT_0344;
    const ccomplex_t IT_0510 = IT_0117*IT_0353;
    const ccomplex_t IT_0511 = (-2)*IT_0416 + (-0.5)*IT_0510;
    const ccomplex_t IT_0512 = IT_0414 + IT_0511;
    const ccomplex_t IT_0513 = IT_0345*IT_0512;
    const ccomplex_t IT_0514 = (-2)*IT_0513;
    const ccomplex_t IT_0515 = (-0.25)*IT_0514;
    const ccomplex_t IT_0516 = s_23*s_34*IT_0117;
    const ccomplex_t IT_0517 = (-0.166666666666667)*IT_0330 + (
      -0.666666666666667)*IT_0516;
    const ccomplex_t IT_0518 = s_24 + IT_0517;
    const ccomplex_t IT_0519 = m_C_2*IT_0518;
    const ccomplex_t IT_0520 = (-6)*IT_0519;
    const ccomplex_t IT_0521 = (-0.25)*IT_0520;
    const ccomplex_t IT_0522 = 0.166666666666667*IT_0357 + (-0.666666666666667
      )*IT_0419;
    const ccomplex_t IT_0523 = s_14 + IT_0522;
    const ccomplex_t IT_0524 = m_C_2*IT_0523;
    const ccomplex_t IT_0525 = 6*IT_0524;
    const ccomplex_t IT_0526 = (-0.25)*IT_0525;
    const ccomplex_t IT_0527 = IT_0076*IT_0326;
    const ccomplex_t IT_0528 = (-10)*s_24;
    const ccomplex_t IT_0529 = IT_0527 + IT_0528;
    const ccomplex_t IT_0530 = m_C_2*IT_0529;
    const ccomplex_t IT_0531 = (-0.25)*IT_0530;
    const ccomplex_t IT_0532 = s_13*s_23*s_24;
    const ccomplex_t IT_0533 = s_13*s_34*IT_0031;
    const ccomplex_t IT_0534 = -IT_0367 + (-2)*IT_0532 + 2*IT_0533;
    const ccomplex_t IT_0535 = IT_0365 + IT_0534;
    const ccomplex_t IT_0536 = IT_0117*IT_0535;
    const ccomplex_t IT_0537 = (-0.5)*IT_0536;
    const ccomplex_t IT_0538 = -IT_0363;
    const ccomplex_t IT_0539 = IT_0371 + IT_0537 + IT_0538;
    const ccomplex_t IT_0540 = (-0.25)*IT_0539;
    const ccomplex_t IT_0541 = IT_0117*IT_0482;
    const ccomplex_t IT_0542 = -IT_0541;
    const ccomplex_t IT_0543 = IT_0383 + IT_0542;
    const ccomplex_t IT_0544 = IT_0313*IT_0508;
    const ccomplex_t IT_0545 = IT_0134*IT_0509;
    const ccomplex_t IT_0546 = 4*conj(IT_0247);
    const ccomplex_t IT_0547 = (-8)*IT_0410*(conj(IT_0403)*IT_0434 + 1./2*conj
      (IT_0162)*IT_0507 + 1./2*IT_0437*IT_0515 + 1./2*IT_0442*IT_0521 + 1./2
      *IT_0304*IT_0526 + 1./2*IT_0438*IT_0531 + 1./2*conj(IT_0410)*IT_0543 + 1.
      /2*IT_0544 + 1./2*IT_0545 + 1./2*IT_0540*IT_0546);
    const ccomplex_t IT_0548 = (-0.25)*IT_0402;
    const ccomplex_t IT_0549 = (-0.25)*IT_0409;
    const ccomplex_t IT_0550 = (-4)*conj(IT_0432)*(IT_0115*IT_0507 + IT_0424
      *IT_0508 + IT_0266*IT_0509 + IT_0430*IT_0515 + IT_0337*IT_0521 + IT_0362
      *IT_0526 + IT_0351*IT_0531 + IT_0374*IT_0540 + IT_0411*IT_0548 + IT_0404
      *IT_0549);
    const ccomplex_t IT_0551 = IT_0317*IT_0508;
    const ccomplex_t IT_0552 = IT_0123*IT_0509;
    const ccomplex_t IT_0553 = 2*IT_0320*IT_0434 + conj(IT_0115)*IT_0507 +
       IT_0442*IT_0515 + IT_0437*IT_0521 + IT_0438*IT_0526 + IT_0304*IT_0531 +
       IT_0439*IT_0540 + conj(IT_0432)*IT_0543 + IT_0441*IT_0548 + IT_0440
      *IT_0549 + IT_0551 + IT_0552;
    const ccomplex_t IT_0554 = (-0.25)*IT_0429;
    const ccomplex_t IT_0555 = (-0.25)*IT_0334;
    const ccomplex_t IT_0556 = (-0.25)*IT_0361;
    const ccomplex_t IT_0557 = (-0.25)*IT_0350;
    const ccomplex_t IT_0558 = (-0.25)*IT_0401;
    const ccomplex_t IT_0559 = (-0.25)*IT_0373;
    const ccomplex_t IT_0560 = (-4)*IT_0403*(IT_0388*conj(IT_0403) + 2*conj
      (IT_0410)*IT_0434 + IT_0544 + IT_0545 + IT_0437*IT_0554 + IT_0442*IT_0555 
      + IT_0304*IT_0556 + IT_0438*IT_0557 + IT_0448*IT_0558 + IT_0546*IT_0559);
    const ccomplex_t IT_0561 = 0.5*IT_0320*IT_0388 + conj(IT_0432)*IT_0434 +
       0.5*IT_0440*IT_0548 + 0.5*IT_0441*IT_0549 + 0.5*IT_0551 + 0.5*IT_0552 +
       0.5*IT_0442*IT_0554 + 0.5*IT_0437*IT_0555 + 0.5*IT_0438*IT_0556 + 0.5
      *IT_0304*IT_0557 + 0.5*IT_0295*IT_0558 + 0.5*IT_0439*IT_0559;
    const ccomplex_t IT_0562 = IT_0189*IT_0278;
    const ccomplex_t IT_0563 = IT_0257*IT_0447;
    const ccomplex_t IT_0564 = IT_0562 + IT_0563;
    const ccomplex_t IT_0565 = IT_0294*IT_0448;
    const ccomplex_t IT_0566 = (-4)*IT_0229*IT_0466;
    const ccomplex_t IT_0567 = (-16)*IT_0221*IT_0472;
    const ccomplex_t IT_0568 = 2*IT_0115*IT_0265 + IT_0566 + IT_0567;
    const ccomplex_t IT_0569 = IT_0294*IT_0444;
    const ccomplex_t IT_0570 = IT_0257*IT_0491;
    const ccomplex_t IT_0571 = 4*(IT_0054*conj(IT_0054) + IT_0075*conj(IT_0075
      ))*IT_0081 + IT_0115*(IT_0122*IT_0123 + IT_0133*IT_0134) + (IT_0123
      *IT_0133 + IT_0122*IT_0134)*IT_0162 + IT_0176*(IT_0188*IT_0189 + IT_0203
      *IT_0208) + IT_0209*((conj(IT_0075) + conj(IT_0217) + conj(IT_0221))
      *IT_0222 + (IT_0075 + IT_0217 + IT_0221)*IT_0223 + (conj(IT_0176) + conj
      (IT_0207) + conj(IT_0225))*IT_0226 + (IT_0176 + IT_0227)*IT_0228) + 2*conj
      (IT_0115)*IT_0271 + IT_0221*(IT_0265*IT_0272 + IT_0278*IT_0279 + IT_0285
      *IT_0286 + IT_0294*IT_0295) + IT_0247*(IT_0235*IT_0272 + IT_0313*IT_0316 +
       IT_0312*IT_0317) + IT_0229*(IT_0303*IT_0304 + IT_0312*IT_0313 + IT_0316
      *IT_0317) + IT_0320*IT_0431 + (-8)*IT_0319*IT_0443 + 4*conj(IT_0162)*
      (IT_0176*IT_0243 + 1.5*IT_0115*IT_0246 + (-1.5)*IT_0229*IT_0249 + 0.25
      *IT_0133*IT_0266 + 0.25*IT_0122*IT_0267 + 0.25*IT_0269*IT_0444) + conj
      (IT_0403)*(IT_0267*IT_0344 + IT_0351*IT_0361 + IT_0350*IT_0362 + IT_0337
      *IT_0429 + IT_0334*IT_0430 + IT_0401*IT_0444 + IT_0373*IT_0445 + IT_0423
      *IT_0446) + IT_0217*(IT_0188*IT_0279 + IT_0208*IT_0285 + IT_0303*IT_0447 +
       IT_0257*IT_0448) + (-6)*conj(IT_0247)*IT_0471 + 4*conj(IT_0176)*(IT_0162
      *IT_0243 + -1./2*IT_0115*IT_0257 + 1./4*IT_0456*IT_0470 + (-4)*IT_0217
      *IT_0472 + 1./4*IT_0188*IT_0473 + 1./4*IT_0203*IT_0474 + 1./4*IT_0485
      *IT_0486) + 2*conj(IT_0229)*IT_0488 + conj(IT_0217)*IT_0495 + conj(IT_0221
      )*IT_0499 + (-4)*conj(IT_0410)*(IT_0162*IT_0507 + IT_0446*IT_0508 +
       IT_0267*IT_0509 + IT_0337*IT_0515 + IT_0430*IT_0521 + IT_0351*IT_0526 +
       IT_0362*IT_0531 + IT_0445*IT_0540) + IT_0547 + IT_0550 + (-4)*IT_0432
      *IT_0553 + IT_0560 + (-8)*IT_0318*IT_0561 + IT_0207*(IT_0564 + IT_0565) +
       IT_0225*(IT_0203*IT_0286 + IT_0208*IT_0498 + IT_0564 + IT_0565) + conj
      (IT_0207)*(IT_0278*IT_0473 + IT_0474*IT_0498 + IT_0568 + IT_0569 + IT_0570
      ) + 8*conj(IT_0225)*(IT_0217*IT_0285 + 0.125*IT_0278*IT_0473 + 0.125
      *IT_0203*IT_0486 + IT_0490*IT_0498 + 0.125*IT_0568 + 0.125*IT_0569 + 0.125
      *IT_0570);
    return create_ccomplex_return(IT_0571);
}

