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
    const ccomplex_t IT_0000 = pow(m_Z, -2);
    const ccomplex_t IT_0001 = pow(m_Z, 2);
    const ccomplex_t IT_0002 = s_24*IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = s_24 + IT_0004;
    const ccomplex_t IT_0006 = m_C_1*IT_0005;
    const ccomplex_t IT_0007 = 2*IT_0006;
    const ccomplex_t IT_0008 = cos(alpha);
    const ccomplex_t IT_0009 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (IT_0010 + -IT_0013)*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*e_em*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = cos(theta_W);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0014*IT_0021;
    const ccomplex_t IT_0023 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = IT_0015*IT_0020;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0024 + -IT_0026 + (-2)
      *IT_0028);
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = IT_0019*IT_0030;
    const ccomplex_t IT_0032 = pow(m_h, 2);
    const ccomplex_t IT_0033 = pow(m_C_1, 2);
    const ccomplex_t IT_0034 = cpow((-2)*s_13 + IT_0032 + IT_0033 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0031*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0032 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0038 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0039 = IT_0008*IT_0038;
    const ccomplex_t IT_0040 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0041 = IT_0011*IT_0040;
    const ccomplex_t IT_0042 = IT_0015*(IT_0039 + -IT_0041);
    const ccomplex_t IT_0043 = 1.4142135623731*e_em*IT_0042;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0047 = IT_0022*IT_0046;
    const ccomplex_t IT_0048 = IT_0025*IT_0046;
    const ccomplex_t IT_0049 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0050 = IT_0025*IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*(IT_0047 + -IT_0048 + (-2)
      *IT_0050);
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = IT_0045*IT_0052;
    const ccomplex_t IT_0054 = IT_0037*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*IT_0054;
    const ccomplex_t IT_0056 = -IT_0036 + -IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0058 = cpow((-2)*s_23 + IT_0032 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0059 = IT_0053*IT_0057*IT_0058;
    const ccomplex_t IT_0060 = V_u1*U_Wm1;
    const ccomplex_t IT_0061 = IT_0008*IT_0060;
    const ccomplex_t IT_0062 = U_d1*V_Wp1;
    const ccomplex_t IT_0063 = IT_0011*IT_0062;
    const ccomplex_t IT_0064 = IT_0015*(IT_0061 + -IT_0063);
    const ccomplex_t IT_0065 = 1.4142135623731*e_em*IT_0064;
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = IT_0052*IT_0067;
    const ccomplex_t IT_0069 = IT_0037*IT_0057*IT_0068;
    const ccomplex_t IT_0070 = m_W*e_em;
    const ccomplex_t IT_0071 = cpow(IT_0020, -4);
    const ccomplex_t IT_0072 = sin(beta);
    const ccomplex_t IT_0073 = IT_0008*IT_0014*IT_0071*IT_0072;
    const ccomplex_t IT_0074 = cos(beta);
    const ccomplex_t IT_0075 = IT_0011*IT_0014*IT_0071*IT_0074;
    const ccomplex_t IT_0076 = cpow(IT_0014, -3);
    const ccomplex_t IT_0077 = IT_0008*IT_0072*IT_0076;
    const ccomplex_t IT_0078 = IT_0011*IT_0074*IT_0076;
    const ccomplex_t IT_0079 = cpow(IT_0020, -2);
    const ccomplex_t IT_0080 = IT_0008*IT_0015*IT_0072*IT_0079;
    const ccomplex_t IT_0081 = IT_0011*IT_0015*IT_0074*IT_0079;
    const ccomplex_t IT_0082 = cpow(IT_0014, -2);
    const ccomplex_t IT_0083 = cpow(IT_0079 + IT_0082, -1);
    const ccomplex_t IT_0084 = (IT_0073 + -IT_0075 + IT_0077 + -IT_0078 + 2
      *IT_0080 + (-2)*IT_0081)*IT_0083;
    const ccomplex_t IT_0085 = IT_0070*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0052*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0033 + -reg_prop, -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0092 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0093 = IT_0008*IT_0092;
    const ccomplex_t IT_0094 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0095 = IT_0011*IT_0094;
    const ccomplex_t IT_0096 = IT_0015*(IT_0093 + -IT_0095);
    const ccomplex_t IT_0097 = 1.4142135623731*e_em*IT_0096;
    const ccomplex_t IT_0098 = 0.5*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0101 = IT_0022*IT_0100;
    const ccomplex_t IT_0102 = IT_0025*IT_0100;
    const ccomplex_t IT_0103 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0104 = IT_0025*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + -IT_0102 + (-2)
      *IT_0104);
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0099*IT_0106;
    const ccomplex_t IT_0108 = cpow((-2)*s_23 + IT_0032 + IT_0033 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0109 = IT_0091*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = V_u1*U_Wm2;
    const ccomplex_t IT_0111 = IT_0008*IT_0110;
    const ccomplex_t IT_0112 = U_d2*V_Wp1;
    const ccomplex_t IT_0113 = IT_0011*IT_0112;
    const ccomplex_t IT_0114 = IT_0015*(IT_0111 + -IT_0113);
    const ccomplex_t IT_0115 = 1.4142135623731*e_em*IT_0114;
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = IT_0030*IT_0117;
    const ccomplex_t IT_0119 = IT_0034*IT_0091*IT_0118;
    const ccomplex_t IT_0120 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0121 = IT_0022*IT_0120;
    const ccomplex_t IT_0122 = IT_0025*IT_0120;
    const ccomplex_t IT_0123 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0124 = IT_0025*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + -IT_0122 + (-2)
      *IT_0124);
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = IT_0067*IT_0126;
    const ccomplex_t IT_0128 = IT_0037*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = m_C_1*IT_0129;
    const ccomplex_t IT_0131 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0132 = IT_0022*IT_0131;
    const ccomplex_t IT_0133 = IT_0025*IT_0131;
    const ccomplex_t IT_0134 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0135 = IT_0025*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + -IT_0133 + (-2)
      *IT_0135);
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0117*IT_0137;
    const ccomplex_t IT_0139 = IT_0034*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_C_1*IT_0140;
    const ccomplex_t IT_0142 = IT_0045*IT_0126;
    const ccomplex_t IT_0143 = IT_0058*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = m_C_1*IT_0144;
    const ccomplex_t IT_0146 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0147 = IT_0022*IT_0146;
    const ccomplex_t IT_0148 = IT_0025*IT_0146;
    const ccomplex_t IT_0149 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0150 = IT_0025*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + -IT_0148 + (-2)
      *IT_0150);
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = IT_0099*IT_0152;
    const ccomplex_t IT_0154 = IT_0108*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_C_1*IT_0155;
    const ccomplex_t IT_0157 = -IT_0059 + -IT_0069 + IT_0090 + -IT_0109 + 
      -IT_0119 + -IT_0130 + -IT_0141 + -IT_0145 + -IT_0156;
    const ccomplex_t IT_0158 = 2*IT_0157;
    const ccomplex_t IT_0159 = 2*conj(IT_0157);
    const ccomplex_t IT_0160 = IT_0129 + IT_0140;
    const ccomplex_t IT_0161 = V_u2*U_Wm1;
    const ccomplex_t IT_0162 = IT_0008*IT_0161;
    const ccomplex_t IT_0163 = U_d1*V_Wp2;
    const ccomplex_t IT_0164 = IT_0011*IT_0163;
    const ccomplex_t IT_0165 = IT_0015*(IT_0162 + -IT_0164);
    const ccomplex_t IT_0166 = 1.4142135623731*e_em*IT_0165;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0152*IT_0168;
    const ccomplex_t IT_0170 = IT_0091*IT_0108*IT_0169;
    const ccomplex_t IT_0171 = IT_0019*IT_0137;
    const ccomplex_t IT_0172 = IT_0034*IT_0091*IT_0171;
    const ccomplex_t IT_0173 = IT_0037*IT_0057*IT_0142;
    const ccomplex_t IT_0174 = IT_0057*IT_0058*IT_0127;
    const ccomplex_t IT_0175 = IT_0058*IT_0068;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = m_C_1*IT_0176;
    const ccomplex_t IT_0178 = IT_0106*IT_0168;
    const ccomplex_t IT_0179 = IT_0108*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = m_C_1*IT_0180;
    const ccomplex_t IT_0182 = IT_0086*IT_0126;
    const ccomplex_t IT_0183 = IT_0088*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = m_C_1*IT_0055;
    const ccomplex_t IT_0186 = m_C_1*IT_0036;
    const ccomplex_t IT_0187 = IT_0170 + IT_0172 + IT_0173 + IT_0174 + IT_0177
       + IT_0181 + -IT_0184 + IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = 2*IT_0187;
    const ccomplex_t IT_0189 = 2*conj(IT_0187);
    const ccomplex_t IT_0190 = (-0.25)*IT_0003;
    const ccomplex_t IT_0191 = s_24 + IT_0190;
    const ccomplex_t IT_0192 = m_C_1*IT_0191;
    const ccomplex_t IT_0193 = (-4)*IT_0192;
    const ccomplex_t IT_0194 = -IT_0144 + -IT_0155;
    const ccomplex_t IT_0195 = conj(IT_0176) + conj(IT_0180);
    const ccomplex_t IT_0196 = IT_0176 + IT_0180;
    const ccomplex_t IT_0197 = e_em*IT_0015*IT_0021*(IT_0011*IT_0072 + IT_0008
      *IT_0074);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0060*IT_0074;
    const ccomplex_t IT_0200 = IT_0062*IT_0072;
    const ccomplex_t IT_0201 = IT_0015*(IT_0199 + IT_0200);
    const ccomplex_t IT_0202 = 1.4142135623731*e_em*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0202;
    const ccomplex_t IT_0204 = cpow((-2)*s_12 + (-2)*IT_0033 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0205 = IT_0203*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = IT_0198*IT_0206;
    const ccomplex_t IT_0208 = e_em*IT_0015*IT_0021*(IT_0008*IT_0072 + 
      -IT_0011*IT_0074);
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = IT_0062*IT_0074;
    const ccomplex_t IT_0211 = IT_0060*IT_0072;
    const ccomplex_t IT_0212 = IT_0015*(IT_0210 + -IT_0211);
    const ccomplex_t IT_0213 = 1.4142135623731*e_em*IT_0212;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = cpow((-2)*s_12 + IT_0001 + (-2)*IT_0033 + 
      -reg_prop, -1);
    const ccomplex_t IT_0216 = IT_0214*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0209*IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0176 + (-2)*IT_0180 + -IT_0207 + 
      -IT_0218;
    const ccomplex_t IT_0220 = pow(s_14, 2);
    const ccomplex_t IT_0221 = IT_0000*IT_0220;
    const ccomplex_t IT_0222 = -IT_0033;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = s_12*IT_0223;
    const ccomplex_t IT_0225 = 2*IT_0219;
    const ccomplex_t IT_0226 = s_14*s_24;
    const ccomplex_t IT_0227 = pow(m_Z, 4);
    const ccomplex_t IT_0228 = s_12*IT_0227;
    const ccomplex_t IT_0229 = IT_0000*IT_0228;
    const ccomplex_t IT_0230 = (-0.25)*IT_0229;
    const ccomplex_t IT_0231 = IT_0226 + IT_0230;
    const ccomplex_t IT_0232 = (-8)*conj(IT_0160);
    const ccomplex_t IT_0233 = s_12*IT_0001;
    const ccomplex_t IT_0234 = (-0.5)*IT_0229;
    const ccomplex_t IT_0235 = (-2)*IT_0226;
    const ccomplex_t IT_0236 = IT_0233 + IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = (-4)*IT_0176;
    const ccomplex_t IT_0238 = s_14*IT_0001;
    const ccomplex_t IT_0239 = IT_0000*IT_0238;
    const ccomplex_t IT_0240 = (-0.25)*IT_0239;
    const ccomplex_t IT_0241 = s_14 + IT_0240;
    const ccomplex_t IT_0242 = m_C_1*IT_0241;
    const ccomplex_t IT_0243 = (-4)*IT_0242;
    const ccomplex_t IT_0244 = (-2)*IT_0157;
    const ccomplex_t IT_0245 = (-2)*conj(IT_0157);
    const ccomplex_t IT_0246 = (-2)*IT_0187;
    const ccomplex_t IT_0247 = (-2)*conj(IT_0187);
    const ccomplex_t IT_0248 = 0.5*IT_0239;
    const ccomplex_t IT_0249 = s_14 + IT_0248;
    const ccomplex_t IT_0250 = m_C_1*IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0250;
    const ccomplex_t IT_0252 = m_C_1*IT_0223;
    const ccomplex_t IT_0253 = (-2)*IT_0219;
    const ccomplex_t IT_0254 = (-2)*IT_0129 + (-2)*IT_0140 + -IT_0207 + 
      -IT_0218;
    const ccomplex_t IT_0255 = pow(s_24, 2);
    const ccomplex_t IT_0256 = IT_0000*IT_0255;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0033 + IT_0257;
    const ccomplex_t IT_0259 = s_12*IT_0258;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = 2*IT_0254;
    const ccomplex_t IT_0262 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = s_12 + IT_0263;
    const ccomplex_t IT_0265 = s_12*IT_0264;
    const ccomplex_t IT_0266 = IT_0040*IT_0074;
    const ccomplex_t IT_0267 = IT_0038*IT_0072;
    const ccomplex_t IT_0268 = IT_0015*(IT_0266 + -IT_0267);
    const ccomplex_t IT_0269 = 1.4142135623731*e_em*IT_0268;
    const ccomplex_t IT_0270 = 0.5*IT_0269;
    const ccomplex_t IT_0271 = IT_0215*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = IT_0209*IT_0272;
    const ccomplex_t IT_0274 = IT_0038*IT_0074;
    const ccomplex_t IT_0275 = IT_0040*IT_0072;
    const ccomplex_t IT_0276 = IT_0015*(IT_0274 + IT_0275);
    const ccomplex_t IT_0277 = 1.4142135623731*e_em*IT_0276;
    const ccomplex_t IT_0278 = (-0.5)*IT_0277;
    const ccomplex_t IT_0279 = IT_0204*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = IT_0198*IT_0280;
    const ccomplex_t IT_0282 = 2*IT_0144 + 2*IT_0155 + -IT_0273 + -IT_0281;
    const ccomplex_t IT_0283 = 2*IT_0282;
    const ccomplex_t IT_0284 = IT_0033*IT_0264;
    const ccomplex_t IT_0285 = IT_0033*IT_0223;
    const ccomplex_t IT_0286 = (-2)*conj(IT_0219);
    const ccomplex_t IT_0287 = 2*conj(IT_0254);
    const ccomplex_t IT_0288 = m_C_1*IT_0264;
    const ccomplex_t IT_0289 = 2*IT_0036 + 2*IT_0055 + -IT_0273 + -IT_0281;
    const ccomplex_t IT_0290 = 2*IT_0289;
    const ccomplex_t IT_0291 = 2*conj(IT_0289);
    const ccomplex_t IT_0292 = m_C_1*IT_0258;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = (-2)*IT_0282;
    const ccomplex_t IT_0295 = s_12*s_14;
    const ccomplex_t IT_0296 = s_24*IT_0033;
    const ccomplex_t IT_0297 = -IT_0296;
    const ccomplex_t IT_0298 = IT_0295 + IT_0297;
    const ccomplex_t IT_0299 = 2*conj(IT_0219);
    const ccomplex_t IT_0300 = (-2)*conj(IT_0282);
    const ccomplex_t IT_0301 = 0.5*IT_0195*IT_0225 + 0.5*conj(IT_0160)*IT_0253
       + IT_0194*conj(IT_0282) + 0.5*conj(IT_0194)*IT_0283 + 0.5*IT_0160*IT_0286
       + 0.5*conj(IT_0056)*IT_0294 + 0.5*IT_0196*IT_0299 + 0.5*IT_0056*IT_0300;
    const ccomplex_t IT_0302 = IT_0033*IT_0258;
    const ccomplex_t IT_0303 = -IT_0302;
    const ccomplex_t IT_0304 = (-2)*IT_0254;
    const ccomplex_t IT_0305 = (-2)*conj(IT_0254);
    const ccomplex_t IT_0306 = s_12*s_24;
    const ccomplex_t IT_0307 = s_14*IT_0033;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = (-0.5)*conj(IT_0160)*IT_0261 + (-0.5)*IT_0160
      *IT_0287 + conj(IT_0194)*IT_0289 + IT_0194*conj(IT_0289) + (-0.5)*conj
      (IT_0056)*IT_0290 + (-0.5)*IT_0056*IT_0291 + (-0.5)*IT_0195*IT_0304 + (
      -0.5)*IT_0196*IT_0305;
    const ccomplex_t IT_0311 = IT_0001*IT_0033;
    const ccomplex_t IT_0312 = 3*IT_0311;
    const ccomplex_t IT_0313 = (-2)*IT_0160;
    const ccomplex_t IT_0314 = (-2)*conj(IT_0160);
    const ccomplex_t IT_0315 = (-2)*IT_0056;
    const ccomplex_t IT_0316 = (-2)*conj(IT_0056);
    const ccomplex_t IT_0317 = (-3)*IT_0311;
    const ccomplex_t IT_0318 = (-2)*IT_0194;
    const ccomplex_t IT_0319 = (-2)*conj(IT_0194);
    const ccomplex_t IT_0320 = IT_0233 + IT_0235;
    const ccomplex_t IT_0321 = IT_0000*IT_0320;
    const ccomplex_t IT_0322 = (-0.5)*IT_0321;
    const ccomplex_t IT_0323 = IT_0007*(conj(IT_0056)*IT_0158 + IT_0056
      *IT_0159 + conj(IT_0160)*IT_0188 + IT_0160*IT_0189) + IT_0193*(IT_0159
      *IT_0194 + IT_0158*conj(IT_0194) + IT_0188*IT_0195 + IT_0189*IT_0196) +
       conj(IT_0219)*IT_0224*IT_0225 + (-8)*IT_0231*(conj(IT_0056)*IT_0194 +
       IT_0056*conj(IT_0194) + IT_0160*IT_0195 + (-0.125)*IT_0196*IT_0232) + (-4
      )*IT_0236*(IT_0056*conj(IT_0056) + IT_0160*conj(IT_0160) + IT_0194*conj
      (IT_0194) + IT_0195*(IT_0180 + (-0.25)*IT_0237)) + IT_0243*(conj(IT_0160)
      *IT_0244 + IT_0160*IT_0245 + conj(IT_0056)*IT_0246 + IT_0056*IT_0247) + 
      (IT_0195*IT_0244 + IT_0196*IT_0245 + conj(IT_0194)*IT_0246 + IT_0194
      *IT_0247)*IT_0251 + conj(IT_0157)*IT_0252*IT_0253 + conj(IT_0254)*(IT_0260
      *IT_0261 + IT_0253*IT_0265) + conj(IT_0282)*(IT_0224*IT_0283 + IT_0261
      *IT_0284 + IT_0253*IT_0285) + IT_0254*IT_0265*IT_0286 + IT_0282*(IT_0285
      *IT_0286 + IT_0284*IT_0287) + IT_0288*(conj(IT_0157)*IT_0261 + IT_0245
      *IT_0282 + IT_0244*conj(IT_0282) + IT_0157*IT_0287 + conj(IT_0187)*
      (IT_0253 + IT_0290) + IT_0187*(IT_0286 + IT_0291)) + conj(IT_0187)*
      (IT_0261*IT_0293 + IT_0252*IT_0294) + 2*IT_0298*IT_0301 + conj(IT_0289)*
      (IT_0225*IT_0284 + IT_0260*IT_0290 + IT_0158*IT_0293 + IT_0265*IT_0294 +
       IT_0303*IT_0304) + IT_0289*(IT_0159*IT_0293 + IT_0284*IT_0299 + IT_0265
      *IT_0300 + IT_0303*IT_0305) + (-2)*IT_0309*IT_0310 + IT_0312*(conj(IT_0194
      )*IT_0313 + IT_0194*IT_0314 + IT_0195*IT_0315 + IT_0196*IT_0316) + IT_0317
      *(conj(IT_0056)*IT_0313 + IT_0056*IT_0314 + IT_0195*IT_0318 + IT_0196
      *IT_0319) + 6*IT_0157*(IT_0033*conj(IT_0187) + 0.166666666666667*IT_0252
      *IT_0286 + 0.666666666666667*conj(IT_0157)*(s_12 + IT_0322)) + 6*IT_0187*
      (IT_0033*conj(IT_0157) + 0.166666666666667*IT_0287*IT_0293 +
       0.166666666666667*IT_0252*IT_0300 + 0.666666666666667*conj(IT_0187)*(s_12
       + IT_0322));
    return create_ccomplex_return(IT_0323);
}

