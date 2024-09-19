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
    const ccomplex_t IT_0014 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = IT_0007*IT_0011;
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0019 = IT_0016*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0015 + -IT_0017 + (-2)
      *IT_0019);
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = IT_0010*IT_0021;
    const ccomplex_t IT_0023 = pow(m_A0, 2);
    const ccomplex_t IT_0024 = pow(m_C_2, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_13 + IT_0023 + IT_0024 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = cpow((-2)*s_13 + IT_0023 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0029 = V_u2*U_Wm2;
    const ccomplex_t IT_0030 = IT_0000*IT_0029;
    const ccomplex_t IT_0031 = U_d2*V_Wp2;
    const ccomplex_t IT_0032 = IT_0003*IT_0031;
    const ccomplex_t IT_0033 = IT_0007*(IT_0030 + IT_0032);
    const ccomplex_t IT_0034 = 1.4142135623731*e_em*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0037 = IT_0013*IT_0036;
    const ccomplex_t IT_0038 = IT_0016*IT_0036;
    const ccomplex_t IT_0039 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0040 = IT_0016*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0038 + (-2)
      *IT_0040);
    const ccomplex_t IT_0042 = 0.5*IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = IT_0028*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = pow(m_Z, 2);
    const ccomplex_t IT_0047 = s_12*IT_0046;
    const ccomplex_t IT_0048 = pow(m_Z, -2);
    const ccomplex_t IT_0049 = pow(m_Z, 4);
    const ccomplex_t IT_0050 = s_12*IT_0049;
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = s_14*s_24;
    const ccomplex_t IT_0054 = (-2)*IT_0053;
    const ccomplex_t IT_0055 = IT_0047 + IT_0052 + IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0023 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0057 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0058 = IT_0013*IT_0057;
    const ccomplex_t IT_0059 = IT_0016*IT_0057;
    const ccomplex_t IT_0060 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0061 = IT_0016*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + -IT_0059 + (-2)
      *IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0035*IT_0063;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = V_u1*U_Wm2;
    const ccomplex_t IT_0068 = IT_0000*IT_0067;
    const ccomplex_t IT_0069 = U_d2*V_Wp1;
    const ccomplex_t IT_0070 = IT_0003*IT_0069;
    const ccomplex_t IT_0071 = IT_0007*(IT_0068 + IT_0070);
    const ccomplex_t IT_0072 = 1.4142135623731*e_em*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0075 = IT_0013*IT_0074;
    const ccomplex_t IT_0076 = IT_0016*IT_0074;
    const ccomplex_t IT_0077 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0078 = IT_0016*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + -IT_0076 + (-2)
      *IT_0078);
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0073*IT_0080;
    const ccomplex_t IT_0082 = cpow((-2)*s_23 + IT_0023 + IT_0024 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = IT_0066 + IT_0084;
    const ccomplex_t IT_0086 = (-0.25)*IT_0051;
    const ccomplex_t IT_0087 = IT_0053 + IT_0086;
    const ccomplex_t IT_0088 = s_24*IT_0046;
    const ccomplex_t IT_0089 = IT_0048*IT_0088;
    const ccomplex_t IT_0090 = 0.5*IT_0089;
    const ccomplex_t IT_0091 = s_24 + IT_0090;
    const ccomplex_t IT_0092 = m_C_2*IT_0091;
    const ccomplex_t IT_0093 = 2*IT_0092;
    const ccomplex_t IT_0094 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0095 = IT_0000*IT_0094;
    const ccomplex_t IT_0096 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0097 = IT_0003*IT_0096;
    const ccomplex_t IT_0098 = IT_0007*(IT_0095 + IT_0097);
    const ccomplex_t IT_0099 = 1.4142135623731*e_em*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0042*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0103 = IT_0028*IT_0101*IT_0102;
    const ccomplex_t IT_0104 = IT_0043*IT_0056*IT_0102;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0106 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0107 = IT_0000*IT_0106;
    const ccomplex_t IT_0108 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0109 = IT_0003*IT_0108;
    const ccomplex_t IT_0110 = IT_0007*(IT_0107 + IT_0109);
    const ccomplex_t IT_0111 = 1.4142135623731*e_em*IT_0110;
    const ccomplex_t IT_0112 = (-0.5)*IT_0111;
    const ccomplex_t IT_0113 = IT_0021*IT_0112;
    const ccomplex_t IT_0114 = IT_0025*IT_0105*IT_0113;
    const ccomplex_t IT_0115 = IT_0063*IT_0100;
    const ccomplex_t IT_0116 = IT_0028*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = m_C_2*IT_0117;
    const ccomplex_t IT_0119 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0120 = IT_0013*IT_0119;
    const ccomplex_t IT_0121 = IT_0016*IT_0119;
    const ccomplex_t IT_0122 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0123 = IT_0016*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + -IT_0121 + (-2)
      *IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0112*IT_0125;
    const ccomplex_t IT_0127 = IT_0025*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = m_C_2*IT_0128;
    const ccomplex_t IT_0130 = m_C_2*IT_0084;
    const ccomplex_t IT_0131 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0132 = IT_0013*IT_0131;
    const ccomplex_t IT_0133 = IT_0016*IT_0131;
    const ccomplex_t IT_0134 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0135 = IT_0016*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + -IT_0133 + (-2)
      *IT_0135);
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0073*IT_0137;
    const ccomplex_t IT_0139 = IT_0082*IT_0105*IT_0138;
    const ccomplex_t IT_0140 = m_C_2*IT_0066;
    const ccomplex_t IT_0141 = IT_0103 + IT_0104 + IT_0114 + IT_0118 + IT_0129
       + IT_0130 + IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = cos(alpha);
    const ccomplex_t IT_0143 = sin(alpha);
    const ccomplex_t IT_0144 = e_em*IT_0007*IT_0012*(IT_0003*IT_0142 + 
      -IT_0000*IT_0143);
    const ccomplex_t IT_0145 = IT_0031*IT_0142;
    const ccomplex_t IT_0146 = IT_0029*IT_0143;
    const ccomplex_t IT_0147 = IT_0007*(IT_0145 + IT_0146);
    const ccomplex_t IT_0148 = 1.4142135623731*e_em*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = cpow((-2)*s_12 + (-2)*IT_0024 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0152 = IT_0150*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = IT_0144*IT_0153;
    const ccomplex_t IT_0155 = e_em*IT_0007*IT_0012*(IT_0000*IT_0142 + IT_0003
      *IT_0143);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0029*IT_0142;
    const ccomplex_t IT_0158 = IT_0031*IT_0143;
    const ccomplex_t IT_0159 = IT_0007*(IT_0157 + -IT_0158);
    const ccomplex_t IT_0160 = 1.4142135623731*e_em*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0164 = IT_0162*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = IT_0156*IT_0165;
    const ccomplex_t IT_0167 = (-2)*IT_0027 + (-2)*IT_0045 + IT_0154 + IT_0166;
    const ccomplex_t IT_0168 = s_12*s_24;
    const ccomplex_t IT_0169 = s_14*IT_0024;
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = 2*conj(IT_0027);
    const ccomplex_t IT_0173 = 2*conj(IT_0045);
    const ccomplex_t IT_0174 = s_14*IT_0046;
    const ccomplex_t IT_0175 = IT_0048*IT_0174;
    const ccomplex_t IT_0176 = (-0.25)*IT_0175;
    const ccomplex_t IT_0177 = s_14 + IT_0176;
    const ccomplex_t IT_0178 = m_C_2*IT_0177;
    const ccomplex_t IT_0179 = (-4)*IT_0178;
    const ccomplex_t IT_0180 = IT_0010*IT_0125;
    const ccomplex_t IT_0181 = IT_0025*IT_0105*IT_0180;
    const ccomplex_t IT_0182 = IT_0056*IT_0101;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = m_C_2*IT_0183;
    const ccomplex_t IT_0185 = m_C_2*IT_0045;
    const ccomplex_t IT_0186 = IT_0028*IT_0064*IT_0102;
    const ccomplex_t IT_0187 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0188 = IT_0000*IT_0187;
    const ccomplex_t IT_0189 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0190 = IT_0003*IT_0189;
    const ccomplex_t IT_0191 = IT_0007*(IT_0188 + IT_0190);
    const ccomplex_t IT_0192 = 1.4142135623731*e_em*IT_0191;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0080*IT_0193;
    const ccomplex_t IT_0195 = IT_0082*IT_0105*IT_0194;
    const ccomplex_t IT_0196 = IT_0137*IT_0193;
    const ccomplex_t IT_0197 = IT_0082*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = m_C_2*IT_0198;
    const ccomplex_t IT_0200 = m_C_2*IT_0027;
    const ccomplex_t IT_0201 = IT_0056*IT_0102*IT_0115;
    const ccomplex_t IT_0202 = -IT_0181 + -IT_0184 + -IT_0185 + -IT_0186 + 
      -IT_0195 + -IT_0199 + -IT_0200 + -IT_0201;
    const ccomplex_t IT_0203 = IT_0024*IT_0046;
    const ccomplex_t IT_0204 = (-3)*IT_0203;
    const ccomplex_t IT_0205 = -IT_0117 + -IT_0128;
    const ccomplex_t IT_0206 = 3*IT_0203;
    const ccomplex_t IT_0207 = -IT_0183 + -IT_0198;
    const ccomplex_t IT_0208 = (-2)*IT_0066 + (-2)*IT_0084 + IT_0154 + IT_0166;
    const ccomplex_t IT_0209 = s_12*s_14;
    const ccomplex_t IT_0210 = s_24*IT_0024;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + IT_0211;
    const ccomplex_t IT_0213 = (-2)*conj(IT_0027);
    const ccomplex_t IT_0214 = (-2)*conj(IT_0045);
    const ccomplex_t IT_0215 = (-2)*IT_0027;
    const ccomplex_t IT_0216 = (-2)*IT_0045;
    const ccomplex_t IT_0217 = IT_0215 + IT_0216;
    const ccomplex_t IT_0218 = IT_0094*IT_0142;
    const ccomplex_t IT_0219 = IT_0096*IT_0143;
    const ccomplex_t IT_0220 = IT_0007*(IT_0218 + -IT_0219);
    const ccomplex_t IT_0221 = 1.4142135623731*e_em*IT_0220;
    const ccomplex_t IT_0222 = 0.5*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = IT_0163*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = IT_0156*IT_0225;
    const ccomplex_t IT_0227 = IT_0096*IT_0142;
    const ccomplex_t IT_0228 = IT_0094*IT_0143;
    const ccomplex_t IT_0229 = IT_0007*(IT_0227 + IT_0228);
    const ccomplex_t IT_0230 = 1.4142135623731*e_em*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = 0.5*IT_0231;
    const ccomplex_t IT_0233 = IT_0151*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = IT_0144*IT_0234;
    const ccomplex_t IT_0236 = 2*IT_0117 + 2*IT_0128 + IT_0226 + IT_0235;
    const ccomplex_t IT_0237 = 2*IT_0171;
    const ccomplex_t IT_0238 = 2*IT_0093;
    const ccomplex_t IT_0239 = (-8)*IT_0087;
    const ccomplex_t IT_0240 = (-4)*IT_0055;
    const ccomplex_t IT_0241 = 2*IT_0183 + 2*IT_0198 + IT_0226 + IT_0235;
    const ccomplex_t IT_0242 = (-2)*IT_0212;
    const ccomplex_t IT_0243 = (-2)*IT_0206;
    const ccomplex_t IT_0244 = (-2)*IT_0179;
    const ccomplex_t IT_0245 = pow(s_24, 2);
    const ccomplex_t IT_0246 = IT_0048*IT_0245;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0024 + IT_0247;
    const ccomplex_t IT_0249 = s_12*IT_0248;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0250;
    const ccomplex_t IT_0252 = s_14*s_24*IT_0048;
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = s_12 + IT_0253;
    const ccomplex_t IT_0255 = m_C_2*IT_0254;
    const ccomplex_t IT_0256 = 2*IT_0255;
    const ccomplex_t IT_0257 = IT_0024*IT_0254;
    const ccomplex_t IT_0258 = 2*IT_0257;
    const ccomplex_t IT_0259 = m_C_2*IT_0248;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = 2*IT_0260;
    const ccomplex_t IT_0262 = (-2)*IT_0171;
    const ccomplex_t IT_0263 = s_12*IT_0254;
    const ccomplex_t IT_0264 = (-2)*IT_0263;
    const ccomplex_t IT_0265 = IT_0024*IT_0248;
    const ccomplex_t IT_0266 = -IT_0265;
    const ccomplex_t IT_0267 = (-2)*IT_0266;
    const ccomplex_t IT_0268 = (-0.25)*IT_0089;
    const ccomplex_t IT_0269 = s_24 + IT_0268;
    const ccomplex_t IT_0270 = m_C_2*IT_0269;
    const ccomplex_t IT_0271 = (-4)*IT_0270;
    const ccomplex_t IT_0272 = 2*IT_0271;
    const ccomplex_t IT_0273 = IT_0047 + IT_0054;
    const ccomplex_t IT_0274 = IT_0048*IT_0273;
    const ccomplex_t IT_0275 = 4*s_12 + (-2)*IT_0274;
    const ccomplex_t IT_0276 = 6*IT_0024;
    const ccomplex_t IT_0277 = pow(s_14, 2);
    const ccomplex_t IT_0278 = IT_0048*IT_0277;
    const ccomplex_t IT_0279 = -IT_0024;
    const ccomplex_t IT_0280 = IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = m_C_2*IT_0280;
    const ccomplex_t IT_0282 = (-2)*IT_0281;
    const ccomplex_t IT_0283 = (-2)*IT_0255;
    const ccomplex_t IT_0284 = 0.5*IT_0175;
    const ccomplex_t IT_0285 = s_14 + IT_0284;
    const ccomplex_t IT_0286 = m_C_2*IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = conj(IT_0202)*(IT_0179*IT_0217 + IT_0205
      *IT_0238 + IT_0167*IT_0256 + IT_0236*IT_0261 + IT_0207*IT_0272 + IT_0202
      *IT_0275 + IT_0141*IT_0276 + IT_0208*IT_0282 + IT_0241*IT_0283 + IT_0085
      *IT_0288);
    const ccomplex_t IT_0290 = 2*IT_0212;
    const ccomplex_t IT_0291 = (-2)*IT_0204;
    const ccomplex_t IT_0292 = (IT_0027 + IT_0045)*IT_0087 + (-0.125)*IT_0085
      *IT_0240 + (-0.125)*IT_0205*IT_0243 + (-0.125)*IT_0167*IT_0262 + (-0.125)
      *IT_0141*IT_0272 + (-0.125)*IT_0202*IT_0288 + (-0.125)*IT_0208*IT_0290 + (
      -0.125)*IT_0207*IT_0291;
    const ccomplex_t IT_0293 = s_12*IT_0280;
    const ccomplex_t IT_0294 = 0.5*IT_0167;
    const ccomplex_t IT_0295 = 0.5*IT_0202;
    const ccomplex_t IT_0296 = IT_0024*IT_0280;
    const ccomplex_t IT_0297 = (-2)*IT_0296;
    const ccomplex_t IT_0298 = 0.5*IT_0241;
    const ccomplex_t IT_0299 = 0.5*IT_0236;
    const ccomplex_t IT_0300 = 0.5*IT_0141;
    const ccomplex_t IT_0301 = 0.5*IT_0085;
    const ccomplex_t IT_0302 = 2*IT_0045;
    const ccomplex_t IT_0303 = 0.5*IT_0302;
    const ccomplex_t IT_0304 = IT_0027 + IT_0303;
    const ccomplex_t IT_0305 = 0.5*IT_0208;
    const ccomplex_t IT_0306 = 0.5*IT_0205;
    const ccomplex_t IT_0307 = 0.5*IT_0207;
    const ccomplex_t IT_0308 = 2*conj(IT_0141)*(IT_0261*IT_0294 + IT_0276
      *IT_0295 + IT_0282*IT_0298 + IT_0256*IT_0299 + IT_0275*IT_0300 + IT_0272
      *IT_0301 + IT_0093*IT_0304 + IT_0283*IT_0305 + IT_0244*IT_0306 + IT_0288
      *IT_0307);
    const ccomplex_t IT_0309 = (-4)*(conj(IT_0027) + conj(IT_0045))*((IT_0027 
      + IT_0045)*IT_0055 + 2*IT_0085*IT_0087) + (IT_0093*IT_0141 + IT_0167
      *IT_0171)*(IT_0172 + IT_0173) + (IT_0179*IT_0202 + IT_0204*IT_0205 +
       IT_0206*IT_0207 + IT_0208*IT_0212)*(IT_0213 + IT_0214) + conj(IT_0205)*
      (IT_0204*IT_0217 + IT_0236*IT_0237 + IT_0202*IT_0238 + IT_0207*IT_0239 +
       IT_0205*IT_0240 + IT_0241*IT_0242 + IT_0085*IT_0243 + IT_0141*IT_0244) +
       conj(IT_0236)*(IT_0205*IT_0237 + IT_0236*IT_0251 + IT_0141*IT_0256 +
       IT_0208*IT_0258 + IT_0202*IT_0261 + IT_0207*IT_0262 + IT_0241*IT_0264 +
       IT_0167*IT_0267) + IT_0289 + conj(IT_0207)*(IT_0206*IT_0217 + IT_0205
      *IT_0239 + IT_0207*IT_0240 + IT_0236*IT_0262 + IT_0202*IT_0272 + IT_0141
      *IT_0288 + IT_0241*IT_0290 + IT_0085*IT_0291) + (-8)*conj(IT_0085)*IT_0292
       + 2*conj(IT_0208)*(0.5*IT_0212*IT_0217 + IT_0208*IT_0293 + IT_0264
      *IT_0294 + IT_0282*IT_0295 + IT_0297*IT_0298 + IT_0258*IT_0299 + IT_0283
      *IT_0300 + IT_0290*IT_0301) + 2*conj(IT_0167)*(IT_0251*IT_0294 + IT_0256
      *IT_0295 + IT_0258*IT_0298 + IT_0267*IT_0299 + IT_0261*IT_0300 + IT_0262
      *IT_0301 + IT_0171*IT_0304 + IT_0264*IT_0305) + 2*conj(IT_0241)*(IT_0241
      *IT_0293 + IT_0258*IT_0294 + IT_0283*IT_0295 + IT_0264*IT_0299 + IT_0282
      *IT_0300 + IT_0297*IT_0305 + IT_0242*IT_0306 + IT_0290*IT_0307) + IT_0308;
    return create_ccomplex_return(IT_0309);
}

