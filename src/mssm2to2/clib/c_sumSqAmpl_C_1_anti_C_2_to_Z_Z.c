#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_Z_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_Z_Z(
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
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = s_12*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(m_Z, 4);
    const ccomplex_t IT_0004 = s_12*IT_0003;
    const ccomplex_t IT_0005 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0006 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0007 = s_13*s_24*s_34;
    const ccomplex_t IT_0008 = (-2)*IT_0005 + (-2)*IT_0006 + 4*IT_0007;
    const ccomplex_t IT_0009 = IT_0004 + IT_0008;
    const ccomplex_t IT_0010 = IT_0002*IT_0009;
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = pow(m_Z, -4);
    const ccomplex_t IT_0013 = pow(m_Z, 6);
    const ccomplex_t IT_0014 = s_12*IT_0013;
    const ccomplex_t IT_0015 = s_13*s_23*IT_0003;
    const ccomplex_t IT_0016 = (-2)*IT_0015;
    const ccomplex_t IT_0017 = IT_0014 + IT_0016;
    const ccomplex_t IT_0018 = IT_0012*IT_0017;
    const ccomplex_t IT_0019 = 0.5*IT_0018;
    const ccomplex_t IT_0020 = s_14*s_24;
    const ccomplex_t IT_0021 = (-4)*IT_0020;
    const ccomplex_t IT_0022 = IT_0001 + IT_0011 + IT_0019 + IT_0021;
    const ccomplex_t IT_0023 = cpow((-2)*s_13 + IT_0000 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0024 = cos(theta_W);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = sin(theta_W);
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = cpow(IT_0026, -1);
    const ccomplex_t IT_0031 = IT_0024*IT_0030;
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0029 + -IT_0032 + (-2)
      *IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0038 = IT_0027*IT_0037;
    const ccomplex_t IT_0039 = IT_0031*IT_0037;
    const ccomplex_t IT_0040 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0041 = IT_0031*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + -IT_0039 + (-2)
      *IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = IT_0036*IT_0043;
    const ccomplex_t IT_0045 = IT_0023*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0048 = IT_0027*IT_0047;
    const ccomplex_t IT_0049 = IT_0031*IT_0047;
    const ccomplex_t IT_0050 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0051 = IT_0031*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + -IT_0049 + (-2)
      *IT_0051);
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = IT_0043*IT_0053;
    const ccomplex_t IT_0055 = pow(m_C_1, 2);
    const ccomplex_t IT_0056 = cpow((-2)*s_13 + IT_0000 + IT_0055 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = conj(IT_0046) + conj(IT_0058);
    const ccomplex_t IT_0060 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0061 = IT_0027*IT_0060;
    const ccomplex_t IT_0062 = IT_0031*IT_0060;
    const ccomplex_t IT_0063 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0064 = IT_0031*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + -IT_0062 + (-2)
      *IT_0064);
    const ccomplex_t IT_0066 = 0.5*IT_0065;
    const ccomplex_t IT_0067 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0068 = IT_0027*IT_0067;
    const ccomplex_t IT_0069 = IT_0031*IT_0067;
    const ccomplex_t IT_0070 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0071 = IT_0031*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + -IT_0069 + (-2)
      *IT_0071);
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0066*IT_0073;
    const ccomplex_t IT_0075 = IT_0056*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0078 = IT_0027*IT_0077;
    const ccomplex_t IT_0079 = IT_0031*IT_0077;
    const ccomplex_t IT_0080 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0081 = IT_0031*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + -IT_0079 + (-2)
      *IT_0081);
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0073*IT_0083;
    const ccomplex_t IT_0085 = IT_0023*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0076 + IT_0086;
    const ccomplex_t IT_0088 = conj(IT_0076) + conj(IT_0086);
    const ccomplex_t IT_0089 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0090 = pow(s_13, 2);
    const ccomplex_t IT_0091 = s_24*IT_0090;
    const ccomplex_t IT_0092 = s_12*s_13*s_34;
    const ccomplex_t IT_0093 = s_24*IT_0000*IT_0055;
    const ccomplex_t IT_0094 = 2*IT_0091 + (-2)*IT_0092 + -IT_0093;
    const ccomplex_t IT_0095 = IT_0089 + IT_0094;
    const ccomplex_t IT_0096 = IT_0002*IT_0095;
    const ccomplex_t IT_0097 = s_24*IT_0055;
    const ccomplex_t IT_0098 = 2*IT_0097;
    const ccomplex_t IT_0099 = s_12*s_14;
    const ccomplex_t IT_0100 = (-2)*IT_0099;
    const ccomplex_t IT_0101 = IT_0096 + IT_0098 + IT_0100;
    const ccomplex_t IT_0102 = pow(m_C_2, 2);
    const ccomplex_t IT_0103 = cpow((-2)*s_23 + IT_0000 + IT_0102 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0044*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = cpow((-2)*s_23 + IT_0000 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0107 = IT_0054*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = 2*IT_0105 + 2*IT_0108;
    const ccomplex_t IT_0110 = 2*IT_0046;
    const ccomplex_t IT_0111 = 2*IT_0058;
    const ccomplex_t IT_0112 = IT_0110 + IT_0111;
    const ccomplex_t IT_0113 = IT_0084*IT_0103;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0074*IT_0106;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 2*IT_0114 + 2*IT_0116;
    const ccomplex_t IT_0118 = 2*IT_0076;
    const ccomplex_t IT_0119 = 2*IT_0086;
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = s_24*IT_0000;
    const ccomplex_t IT_0122 = IT_0002*IT_0121;
    const ccomplex_t IT_0123 = s_34*IT_0012;
    const ccomplex_t IT_0124 = s_23*IT_0000;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = s_23*s_34;
    const ccomplex_t IT_0127 = (-2)*IT_0126;
    const ccomplex_t IT_0128 = IT_0121 + IT_0127;
    const ccomplex_t IT_0129 = IT_0002*IT_0128;
    const ccomplex_t IT_0130 = 0.5*IT_0122 + (-0.5)*IT_0125 + (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = s_24 + IT_0130;
    const ccomplex_t IT_0132 = m_C_1*IT_0131;
    const ccomplex_t IT_0133 = (-2)*IT_0132;
    const ccomplex_t IT_0134 = cos(beta);
    const ccomplex_t IT_0135 = cos(alpha);
    const ccomplex_t IT_0136 = cpow(IT_0024, -4);
    const ccomplex_t IT_0137 = IT_0026*IT_0134*IT_0135*IT_0136;
    const ccomplex_t IT_0138 = sin(beta);
    const ccomplex_t IT_0139 = sin(alpha);
    const ccomplex_t IT_0140 = IT_0026*IT_0136*IT_0138*IT_0139;
    const ccomplex_t IT_0141 = cpow(IT_0026, -3);
    const ccomplex_t IT_0142 = IT_0134*IT_0135*IT_0141;
    const ccomplex_t IT_0143 = IT_0138*IT_0139*IT_0141;
    const ccomplex_t IT_0144 = cpow(IT_0024, -2);
    const ccomplex_t IT_0145 = IT_0030*IT_0134*IT_0135*IT_0144;
    const ccomplex_t IT_0146 = IT_0030*IT_0138*IT_0139*IT_0144;
    const ccomplex_t IT_0147 = cpow(IT_0026, -2);
    const ccomplex_t IT_0148 = cpow(IT_0144 + IT_0147, -1);
    const ccomplex_t IT_0149 = (IT_0137 + IT_0140 + IT_0142 + IT_0143 + 2
      *IT_0145 + 2*IT_0146)*IT_0148;
    const ccomplex_t IT_0150 = m_W*e_em;
    const ccomplex_t IT_0151 = IT_0149*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0154 = IT_0135*IT_0153;
    const ccomplex_t IT_0155 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0156 = IT_0139*IT_0155;
    const ccomplex_t IT_0157 = IT_0030*(IT_0154 + IT_0156);
    const ccomplex_t IT_0158 = 1.4142135623731*e_em*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0159;
    const ccomplex_t IT_0161 = IT_0152*IT_0160;
    const ccomplex_t IT_0162 = cpow((-2)*s_12 + -IT_0055 + -IT_0102 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0026*IT_0135*IT_0136*IT_0138;
    const ccomplex_t IT_0166 = IT_0026*IT_0134*IT_0136*IT_0139;
    const ccomplex_t IT_0167 = IT_0135*IT_0138*IT_0141;
    const ccomplex_t IT_0168 = IT_0134*IT_0139*IT_0141;
    const ccomplex_t IT_0169 = IT_0030*IT_0135*IT_0138*IT_0144;
    const ccomplex_t IT_0170 = IT_0030*IT_0134*IT_0139*IT_0144;
    const ccomplex_t IT_0171 = IT_0148*(IT_0165 + -IT_0166 + IT_0167 + 
      -IT_0168 + 2*IT_0169 + (-2)*IT_0170);
    const ccomplex_t IT_0172 = IT_0150*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = IT_0135*IT_0155;
    const ccomplex_t IT_0175 = IT_0139*IT_0153;
    const ccomplex_t IT_0176 = IT_0030*(IT_0174 + -IT_0175);
    const ccomplex_t IT_0177 = 1.4142135623731*e_em*IT_0176;
    const ccomplex_t IT_0178 = 0.5*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = IT_0173*IT_0179;
    const ccomplex_t IT_0181 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0055 + -IT_0102 + -reg_prop, -1);
    const ccomplex_t IT_0182 = IT_0180*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = IT_0164 + IT_0183;
    const ccomplex_t IT_0185 = V_u1*U_Wm2;
    const ccomplex_t IT_0186 = IT_0135*IT_0185;
    const ccomplex_t IT_0187 = U_d2*V_Wp1;
    const ccomplex_t IT_0188 = IT_0139*IT_0187;
    const ccomplex_t IT_0189 = IT_0030*(IT_0186 + -IT_0188);
    const ccomplex_t IT_0190 = 1.4142135623731*e_em*IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0173*IT_0192;
    const ccomplex_t IT_0194 = IT_0181*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = IT_0135*IT_0187;
    const ccomplex_t IT_0197 = IT_0139*IT_0185;
    const ccomplex_t IT_0198 = IT_0030*(IT_0196 + IT_0197);
    const ccomplex_t IT_0199 = 1.4142135623731*e_em*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = IT_0152*IT_0201;
    const ccomplex_t IT_0203 = IT_0162*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = IT_0195 + IT_0204;
    const ccomplex_t IT_0206 = m_C_1*m_C_2;
    const ccomplex_t IT_0207 = IT_0000*IT_0206;
    const ccomplex_t IT_0208 = 9*IT_0207;
    const ccomplex_t IT_0209 = (-2)*IT_0046;
    const ccomplex_t IT_0210 = (-2)*IT_0058;
    const ccomplex_t IT_0211 = IT_0209 + IT_0210;
    const ccomplex_t IT_0212 = (-2)*conj(IT_0046);
    const ccomplex_t IT_0213 = (-2)*conj(IT_0058);
    const ccomplex_t IT_0214 = IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = s_13*IT_0000;
    const ccomplex_t IT_0216 = IT_0123*IT_0215;
    const ccomplex_t IT_0217 = s_14*IT_0000;
    const ccomplex_t IT_0218 = IT_0002*IT_0217;
    const ccomplex_t IT_0219 = s_13*s_34;
    const ccomplex_t IT_0220 = (-2)*IT_0219;
    const ccomplex_t IT_0221 = IT_0217 + IT_0220;
    const ccomplex_t IT_0222 = IT_0002*IT_0221;
    const ccomplex_t IT_0223 = (-0.5)*IT_0216 + 0.5*IT_0218 + (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = s_14 + IT_0223;
    const ccomplex_t IT_0225 = m_C_2*IT_0224;
    const ccomplex_t IT_0226 = (-2)*IT_0225;
    const ccomplex_t IT_0227 = (-2)*IT_0076;
    const ccomplex_t IT_0228 = (-2)*IT_0086;
    const ccomplex_t IT_0229 = IT_0227 + IT_0228;
    const ccomplex_t IT_0230 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0231 = s_13*s_23*s_24;
    const ccomplex_t IT_0232 = s_14*IT_0000*IT_0102;
    const ccomplex_t IT_0233 = s_13*s_34*IT_0102;
    const ccomplex_t IT_0234 = (-2)*IT_0231 + -IT_0232 + 2*IT_0233;
    const ccomplex_t IT_0235 = IT_0230 + IT_0234;
    const ccomplex_t IT_0236 = IT_0002*IT_0235;
    const ccomplex_t IT_0237 = s_14*IT_0102;
    const ccomplex_t IT_0238 = 2*IT_0237;
    const ccomplex_t IT_0239 = s_12*s_24;
    const ccomplex_t IT_0240 = (-2)*IT_0239;
    const ccomplex_t IT_0241 = IT_0236 + IT_0238 + IT_0240;
    const ccomplex_t IT_0242 = (-2)*IT_0046 + (-2)*IT_0058;
    const ccomplex_t IT_0243 = (-2)*IT_0076 + (-2)*IT_0086;
    const ccomplex_t IT_0244 = 2*conj(IT_0076);
    const ccomplex_t IT_0245 = 2*conj(IT_0086);
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = s_14*s_23*s_34*IT_0012;
    const ccomplex_t IT_0248 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0249 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0250 = -IT_0248 + -IT_0249;
    const ccomplex_t IT_0251 = s_12 + IT_0247 + IT_0250;
    const ccomplex_t IT_0252 = m_C_1*IT_0251;
    const ccomplex_t IT_0253 = 2*IT_0252;
    const ccomplex_t IT_0254 = IT_0002*IT_0090;
    const ccomplex_t IT_0255 = pow(s_14, 2);
    const ccomplex_t IT_0256 = IT_0002*IT_0255;
    const ccomplex_t IT_0257 = IT_0254 + IT_0256;
    const ccomplex_t IT_0258 = s_13*s_14*s_34*IT_0012;
    const ccomplex_t IT_0259 = -IT_0055 + -IT_0258;
    const ccomplex_t IT_0260 = IT_0257 + IT_0259;
    const ccomplex_t IT_0261 = m_C_2*IT_0260;
    const ccomplex_t IT_0262 = 2*IT_0261;
    const ccomplex_t IT_0263 = s_13*s_23;
    const ccomplex_t IT_0264 = (-2)*IT_0263;
    const ccomplex_t IT_0265 = IT_0001 + IT_0264;
    const ccomplex_t IT_0266 = IT_0002*IT_0265;
    const ccomplex_t IT_0267 = (s_12 + -IT_0249)*(s_12 + (-0.5)*IT_0266);
    const ccomplex_t IT_0268 = (-2)*IT_0267;
    const ccomplex_t IT_0269 = 2*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0271 = IT_0036*IT_0073;
    const ccomplex_t IT_0272 = IT_0103*IT_0270*IT_0271;
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = IT_0043*IT_0066;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0276 = IT_0106*IT_0274*IT_0275;
    const ccomplex_t IT_0277 = m_C_1*IT_0114;
    const ccomplex_t IT_0278 = m_C_1*IT_0116;
    const ccomplex_t IT_0279 = -IT_0276 + -IT_0277 + -IT_0278;
    const ccomplex_t IT_0280 = conj(IT_0273) + conj(IT_0279);
    const ccomplex_t IT_0281 = s_14*s_23*s_34*IT_0002;
    const ccomplex_t IT_0282 = -IT_0020 + (-2)*IT_0263 + 2*IT_0281;
    const ccomplex_t IT_0283 = IT_0001 + IT_0282;
    const ccomplex_t IT_0284 = IT_0002*IT_0283;
    const ccomplex_t IT_0285 = -IT_0249 + (-0.5)*IT_0284;
    const ccomplex_t IT_0286 = s_12 + IT_0285;
    const ccomplex_t IT_0287 = m_C_1*IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = 2*IT_0288;
    const ccomplex_t IT_0290 = -IT_0114 + -IT_0116;
    const ccomplex_t IT_0291 = s_13*s_14*s_23;
    const ccomplex_t IT_0292 = s_23*s_34*IT_0055;
    const ccomplex_t IT_0293 = -IT_0093 + (-2)*IT_0291 + 2*IT_0292;
    const ccomplex_t IT_0294 = IT_0089 + IT_0293;
    const ccomplex_t IT_0295 = IT_0002*IT_0294;
    const ccomplex_t IT_0296 = (-0.5)*IT_0295;
    const ccomplex_t IT_0297 = -IT_0097;
    const ccomplex_t IT_0298 = IT_0099 + IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = 4*IT_0298;
    const ccomplex_t IT_0300 = -IT_0055;
    const ccomplex_t IT_0301 = IT_0256 + IT_0300;
    const ccomplex_t IT_0302 = IT_0206*IT_0301;
    const ccomplex_t IT_0303 = 6*IT_0302;
    const ccomplex_t IT_0304 = -IT_0249;
    const ccomplex_t IT_0305 = s_12 + IT_0304;
    const ccomplex_t IT_0306 = IT_0206*IT_0305;
    const ccomplex_t IT_0307 = (-6)*IT_0306;
    const ccomplex_t IT_0308 = -(IT_0055 + -IT_0256)*(s_12 + (-0.5)*IT_0266);
    const ccomplex_t IT_0309 = (-2)*IT_0308;
    const ccomplex_t IT_0310 = (-2)*IT_0309;
    const ccomplex_t IT_0311 = IT_0056*IT_0274*IT_0275;
    const ccomplex_t IT_0312 = m_C_2*IT_0058;
    const ccomplex_t IT_0313 = IT_0023*IT_0270*IT_0271;
    const ccomplex_t IT_0314 = -IT_0311 + -IT_0312 + -IT_0313;
    const ccomplex_t IT_0315 = m_C_2*IT_0046;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = conj(IT_0314) + conj(IT_0316);
    const ccomplex_t IT_0318 = -IT_0020;
    const ccomplex_t IT_0319 = IT_0001 + IT_0318;
    const ccomplex_t IT_0320 = IT_0002*IT_0319;
    const ccomplex_t IT_0321 = -IT_0249 + (-0.25)*IT_0320;
    const ccomplex_t IT_0322 = s_12 + IT_0321;
    const ccomplex_t IT_0323 = m_C_1*IT_0322;
    const ccomplex_t IT_0324 = (-4)*IT_0323;
    const ccomplex_t IT_0325 = (-2)*IT_0324;
    const ccomplex_t IT_0326 = IT_0117*(IT_0101*IT_0246 + conj(IT_0205)
      *IT_0253 + conj(IT_0184)*IT_0262 + conj(IT_0243)*IT_0269 + IT_0280*IT_0289
       + conj(IT_0290)*IT_0299 + conj(IT_0109)*IT_0303 + conj(IT_0242)*IT_0307 +
       conj(IT_0117)*IT_0310 + IT_0317*IT_0325);
    const ccomplex_t IT_0327 = 2*conj(IT_0046);
    const ccomplex_t IT_0328 = 2*conj(IT_0058);
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = m_C_1*IT_0108;
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = IT_0043*IT_0083;
    const ccomplex_t IT_0333 = IT_0103*IT_0270*IT_0332;
    const ccomplex_t IT_0334 = m_C_1*IT_0105;
    const ccomplex_t IT_0335 = IT_0053*IT_0073;
    const ccomplex_t IT_0336 = IT_0106*IT_0275*IT_0335;
    const ccomplex_t IT_0337 = -IT_0333 + -IT_0334 + -IT_0336;
    const ccomplex_t IT_0338 = conj(IT_0331) + conj(IT_0337);
    const ccomplex_t IT_0339 = s_13*s_14*s_34*IT_0002;
    const ccomplex_t IT_0340 = IT_0000*IT_0055;
    const ccomplex_t IT_0341 = 0.5*IT_0255 + -IT_0339 + (-0.5)*IT_0340;
    const ccomplex_t IT_0342 = IT_0090 + IT_0341;
    const ccomplex_t IT_0343 = IT_0002*IT_0342;
    const ccomplex_t IT_0344 = -IT_0256;
    const ccomplex_t IT_0345 = IT_0055 + IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = m_C_2*IT_0345;
    const ccomplex_t IT_0347 = 4*IT_0346;
    const ccomplex_t IT_0348 = -IT_0105 + -IT_0108;
    const ccomplex_t IT_0349 = IT_0056*IT_0275*IT_0335;
    const ccomplex_t IT_0350 = m_C_2*IT_0076;
    const ccomplex_t IT_0351 = -IT_0349 + -IT_0350;
    const ccomplex_t IT_0352 = IT_0023*IT_0270*IT_0332;
    const ccomplex_t IT_0353 = m_C_2*IT_0086;
    const ccomplex_t IT_0354 = -IT_0352 + -IT_0353;
    const ccomplex_t IT_0355 = conj(IT_0351) + conj(IT_0354);
    const ccomplex_t IT_0356 = -IT_0340;
    const ccomplex_t IT_0357 = IT_0255 + IT_0356;
    const ccomplex_t IT_0358 = IT_0002*IT_0357;
    const ccomplex_t IT_0359 = (-4)*IT_0301;
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = m_C_2*IT_0360;
    const ccomplex_t IT_0362 = (-2)*IT_0361;
    const ccomplex_t IT_0363 = conj(IT_0184)*IT_0253 + conj(IT_0205)*IT_0262 +
       conj(IT_0242)*IT_0269 + conj(IT_0117)*IT_0303 + conj(IT_0243)*IT_0307 +
       conj(IT_0109)*IT_0310 + IT_0101*IT_0329 + IT_0289*IT_0338 + IT_0317
      *IT_0347 + IT_0299*conj(IT_0348) + IT_0325*IT_0355 + IT_0280*IT_0362;
    const ccomplex_t IT_0364 = s_14*IT_0003;
    const ccomplex_t IT_0365 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0366 = (-2)*IT_0365;
    const ccomplex_t IT_0367 = IT_0364 + IT_0366;
    const ccomplex_t IT_0368 = IT_0012*IT_0367;
    const ccomplex_t IT_0369 = (-0.666666666666667)*IT_0222 +
       0.166666666666667*IT_0368;
    const ccomplex_t IT_0370 = s_14 + IT_0369;
    const ccomplex_t IT_0371 = m_C_2*IT_0370;
    const ccomplex_t IT_0372 = 6*IT_0371;
    const ccomplex_t IT_0373 = IT_0110*IT_0372;
    const ccomplex_t IT_0374 = IT_0111*IT_0372;
    const ccomplex_t IT_0375 = s_24*IT_0003;
    const ccomplex_t IT_0376 = IT_0012*IT_0375;
    const ccomplex_t IT_0377 = (-10)*s_24;
    const ccomplex_t IT_0378 = IT_0376 + IT_0377;
    const ccomplex_t IT_0379 = m_C_1*IT_0378;
    const ccomplex_t IT_0380 = IT_0118*IT_0379;
    const ccomplex_t IT_0381 = IT_0119*IT_0379;
    const ccomplex_t IT_0382 = pow(s_34, 2);
    const ccomplex_t IT_0383 = IT_0012*IT_0382;
    const ccomplex_t IT_0384 = 0.5*IT_0383;
    const ccomplex_t IT_0385 = 1 + IT_0384;
    const ccomplex_t IT_0386 = IT_0206*IT_0385;
    const ccomplex_t IT_0387 = (-4)*IT_0386;
    const ccomplex_t IT_0388 = IT_0184*IT_0387;
    const ccomplex_t IT_0389 = m_C_2*IT_0002;
    const ccomplex_t IT_0390 = IT_0002*IT_0364;
    const ccomplex_t IT_0391 = (-2)*IT_0219 + (-0.5)*IT_0390;
    const ccomplex_t IT_0392 = IT_0217 + IT_0391;
    const ccomplex_t IT_0393 = IT_0389*IT_0392;
    const ccomplex_t IT_0394 = (-2)*IT_0393;
    const ccomplex_t IT_0395 = IT_0209*IT_0394;
    const ccomplex_t IT_0396 = IT_0210*IT_0394;
    const ccomplex_t IT_0397 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0398 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0399 = (-2)*IT_0398;
    const ccomplex_t IT_0400 = IT_0375 + IT_0399;
    const ccomplex_t IT_0401 = IT_0012*IT_0400;
    const ccomplex_t IT_0402 = (-0.666666666666667)*IT_0397 + (
      -0.166666666666667)*IT_0401;
    const ccomplex_t IT_0403 = s_24 + IT_0402;
    const ccomplex_t IT_0404 = m_C_1*IT_0403;
    const ccomplex_t IT_0405 = (-6)*IT_0404;
    const ccomplex_t IT_0406 = IT_0227*IT_0405;
    const ccomplex_t IT_0407 = IT_0228*IT_0405;
    const ccomplex_t IT_0408 = IT_0110*IT_0379;
    const ccomplex_t IT_0409 = IT_0111*IT_0379;
    const ccomplex_t IT_0410 = IT_0118*IT_0372;
    const ccomplex_t IT_0411 = IT_0119*IT_0372;
    const ccomplex_t IT_0412 = IT_0117*IT_0347;
    const ccomplex_t IT_0413 = (-0.5)*IT_0266;
    const ccomplex_t IT_0414 = s_12 + IT_0413;
    const ccomplex_t IT_0415 = (-2)*IT_0020;
    const ccomplex_t IT_0416 = IT_0001 + IT_0415;
    const ccomplex_t IT_0417 = IT_0002*IT_0416;
    const ccomplex_t IT_0418 = (-0.5)*IT_0417;
    const ccomplex_t IT_0419 = IT_0009*IT_0012;
    const ccomplex_t IT_0420 = 0.25*IT_0419;
    const ccomplex_t IT_0421 = IT_0414 + IT_0418 + IT_0420;
    const ccomplex_t IT_0422 = 8*IT_0421;
    const ccomplex_t IT_0423 = IT_0351*IT_0422;
    const ccomplex_t IT_0424 = s_14*s_23*s_34;
    const ccomplex_t IT_0425 = s_12*IT_0382;
    const ccomplex_t IT_0426 = (-2)*IT_0005 + (-2)*IT_0006 + 2*IT_0007 + 2
      *IT_0424 + (-2)*IT_0425;
    const ccomplex_t IT_0427 = IT_0004 + IT_0426;
    const ccomplex_t IT_0428 = IT_0012*IT_0427;
    const ccomplex_t IT_0429 = (-0.25)*IT_0428;
    const ccomplex_t IT_0430 = -IT_0248;
    const ccomplex_t IT_0431 = s_12 + IT_0304 + IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 8*IT_0431;
    const ccomplex_t IT_0433 = IT_0337*IT_0432;
    const ccomplex_t IT_0434 = (-2)*IT_0382;
    const ccomplex_t IT_0435 = IT_0003 + IT_0434;
    const ccomplex_t IT_0436 = IT_0012*IT_0435;
    const ccomplex_t IT_0437 = 0.25*IT_0436;
    const ccomplex_t IT_0438 = 1 + IT_0437;
    const ccomplex_t IT_0439 = IT_0206*IT_0438;
    const ccomplex_t IT_0440 = 8*IT_0439;
    const ccomplex_t IT_0441 = IT_0273*IT_0440;
    const ccomplex_t IT_0442 = IT_0279*IT_0440;
    const ccomplex_t IT_0443 = 1.125*IT_0206;
    const ccomplex_t IT_0444 = (-16)*IT_0443;
    const ccomplex_t IT_0445 = IT_0316*IT_0444;
    const ccomplex_t IT_0446 = s_13*s_24*s_34*IT_0002;
    const ccomplex_t IT_0447 = -IT_0020 + (-2)*IT_0263 + 2*IT_0446;
    const ccomplex_t IT_0448 = IT_0001 + IT_0447;
    const ccomplex_t IT_0449 = IT_0002*IT_0448;
    const ccomplex_t IT_0450 = -IT_0249 + (-0.5)*IT_0449;
    const ccomplex_t IT_0451 = s_12 + IT_0450;
    const ccomplex_t IT_0452 = m_C_2*IT_0451;
    const ccomplex_t IT_0453 = (-2)*IT_0452;
    const ccomplex_t IT_0454 = (-2)*IT_0453;
    const ccomplex_t IT_0455 = IT_0243*IT_0454;
    const ccomplex_t IT_0456 = IT_0205*IT_0387;
    const ccomplex_t IT_0457 = m_C_2*IT_0322;
    const ccomplex_t IT_0458 = (-4)*IT_0457;
    const ccomplex_t IT_0459 = 2*IT_0458;
    const ccomplex_t IT_0460 = IT_0243*IT_0459;
    const ccomplex_t IT_0461 = (-2)*IT_0005 + (-2)*IT_0006 + 4*IT_0424;
    const ccomplex_t IT_0462 = IT_0004 + IT_0461;
    const ccomplex_t IT_0463 = IT_0012*IT_0462;
    const ccomplex_t IT_0464 = 0.25*IT_0463;
    const ccomplex_t IT_0465 = IT_0414 + IT_0418 + IT_0464;
    const ccomplex_t IT_0466 = 8*IT_0465;
    const ccomplex_t IT_0467 = IT_0337*IT_0466;
    const ccomplex_t IT_0468 = IT_0316*IT_0440;
    const ccomplex_t IT_0469 = IT_0273*IT_0444;
    const ccomplex_t IT_0470 = IT_0279*IT_0444;
    const ccomplex_t IT_0471 = IT_0209*IT_0405;
    const ccomplex_t IT_0472 = IT_0210*IT_0405;
    const ccomplex_t IT_0473 = IT_0227*IT_0394;
    const ccomplex_t IT_0474 = IT_0228*IT_0394;
    const ccomplex_t IT_0475 = IT_0117*IT_0362;
    const ccomplex_t IT_0476 = (-0.25)*IT_0133;
    const ccomplex_t IT_0477 = (-0.25)*conj(IT_0109);
    const ccomplex_t IT_0478 = s_13*s_24*s_34*IT_0012;
    const ccomplex_t IT_0479 = s_12 + IT_0250 + IT_0478;
    const ccomplex_t IT_0480 = m_C_2*IT_0479;
    const ccomplex_t IT_0481 = (-2)*IT_0480;
    const ccomplex_t IT_0482 = (-0.25)*conj(IT_0242);
    const ccomplex_t IT_0483 = (-0.25)*conj(IT_0117);
    const ccomplex_t IT_0484 = IT_0001*IT_0002;
    const ccomplex_t IT_0485 = (-0.5)*IT_0484;
    const ccomplex_t IT_0486 = s_14*s_23;
    const ccomplex_t IT_0487 = s_13*s_24;
    const ccomplex_t IT_0488 = s_12*s_34;
    const ccomplex_t IT_0489 = -IT_0487 + -IT_0488;
    const ccomplex_t IT_0490 = IT_0486 + IT_0489;
    const ccomplex_t IT_0491 = IT_0123*IT_0490;
    const ccomplex_t IT_0492 = (-0.25)*IT_0491;
    const ccomplex_t IT_0493 = s_12 + IT_0485 + IT_0492;
    const ccomplex_t IT_0494 = 8*IT_0493;
    const ccomplex_t IT_0495 = (-0.25)*IT_0317;
    const ccomplex_t IT_0496 = (-2)*IT_0226;
    const ccomplex_t IT_0497 = (-0.25)*conj(IT_0348);
    const ccomplex_t IT_0498 = s_12*IT_0385;
    const ccomplex_t IT_0499 = 4*IT_0498;
    const ccomplex_t IT_0500 = (-0.25)*conj(IT_0205);
    const ccomplex_t IT_0501 = pow(s_23, 2);
    const ccomplex_t IT_0502 = IT_0002*IT_0501;
    const ccomplex_t IT_0503 = pow(s_24, 2);
    const ccomplex_t IT_0504 = IT_0002*IT_0503;
    const ccomplex_t IT_0505 = IT_0502 + IT_0504;
    const ccomplex_t IT_0506 = s_23*s_24*s_34*IT_0012;
    const ccomplex_t IT_0507 = -IT_0102 + -IT_0506;
    const ccomplex_t IT_0508 = IT_0505 + IT_0507;
    const ccomplex_t IT_0509 = m_C_1*IT_0508;
    const ccomplex_t IT_0510 = (-2)*IT_0509;
    const ccomplex_t IT_0511 = (-0.25)*conj(IT_0243);
    const ccomplex_t IT_0512 = IT_0486 + IT_0488;
    const ccomplex_t IT_0513 = -IT_0487;
    const ccomplex_t IT_0514 = IT_0512 + IT_0513;
    const ccomplex_t IT_0515 = IT_0123*IT_0514;
    const ccomplex_t IT_0516 = 0.25*IT_0515;
    const ccomplex_t IT_0517 = s_12 + IT_0485 + IT_0516;
    const ccomplex_t IT_0518 = 8*IT_0517;
    const ccomplex_t IT_0519 = (-0.25)*IT_0280;
    const ccomplex_t IT_0520 = 2*IT_0133;
    const ccomplex_t IT_0521 = (-0.25)*conj(IT_0290);
    const ccomplex_t IT_0522 = (-0.25)*IT_0214*IT_0226 + conj(IT_0184)*IT_0386
       + IT_0246*IT_0476 + IT_0262*IT_0477 + IT_0481*IT_0482 + IT_0253*IT_0483 +
       IT_0494*IT_0495 + IT_0496*IT_0497 + IT_0499*IT_0500 + IT_0510*IT_0511 +
       IT_0518*IT_0519 + IT_0520*IT_0521;
    const ccomplex_t IT_0523 = (-2)*conj(IT_0076);
    const ccomplex_t IT_0524 = (-2)*conj(IT_0086);
    const ccomplex_t IT_0525 = IT_0523 + IT_0524;
    const ccomplex_t IT_0526 = (-0.25)*conj(IT_0184);
    const ccomplex_t IT_0527 = (-0.25)*IT_0338;
    const ccomplex_t IT_0528 = (-0.25)*IT_0355;
    const ccomplex_t IT_0529 = conj(IT_0205)*IT_0386 + IT_0329*IT_0476 +
       IT_0253*IT_0477 + IT_0262*IT_0483 + IT_0482*IT_0510 + IT_0481*IT_0511 +
       IT_0497*IT_0520 + IT_0496*IT_0521 + (-0.25)*IT_0226*IT_0525 + IT_0499
      *IT_0526 + IT_0518*IT_0527 + IT_0494*IT_0528;
    const ccomplex_t IT_0530 = s_23*s_24*s_34*IT_0002;
    const ccomplex_t IT_0531 = IT_0000*IT_0102;
    const ccomplex_t IT_0532 = 0.5*IT_0503 + -IT_0530 + (-0.5)*IT_0531;
    const ccomplex_t IT_0533 = IT_0501 + IT_0532;
    const ccomplex_t IT_0534 = IT_0002*IT_0533;
    const ccomplex_t IT_0535 = -IT_0504;
    const ccomplex_t IT_0536 = IT_0102 + IT_0534 + IT_0535;
    const ccomplex_t IT_0537 = m_C_1*IT_0536;
    const ccomplex_t IT_0538 = -(s_12 + (-0.5)*IT_0266)*(IT_0102 + -IT_0504);
    const ccomplex_t IT_0539 = (-2)*IT_0538;
    const ccomplex_t IT_0540 = (-2)*IT_0539;
    const ccomplex_t IT_0541 = s_14*IT_0501;
    const ccomplex_t IT_0542 = s_12*s_23*s_34;
    const ccomplex_t IT_0543 = -IT_0232 + 2*IT_0541 + (-2)*IT_0542;
    const ccomplex_t IT_0544 = IT_0230 + IT_0543;
    const ccomplex_t IT_0545 = IT_0002*IT_0544;
    const ccomplex_t IT_0546 = 0.5*IT_0545;
    const ccomplex_t IT_0547 = -IT_0546;
    const ccomplex_t IT_0548 = -IT_0237;
    const ccomplex_t IT_0549 = IT_0239 + IT_0547 + IT_0548;
    const ccomplex_t IT_0550 = (-4)*IT_0549;
    const ccomplex_t IT_0551 = -IT_0102;
    const ccomplex_t IT_0552 = IT_0504 + IT_0551;
    const ccomplex_t IT_0553 = IT_0206*IT_0552;
    const ccomplex_t IT_0554 = 6*IT_0553;
    const ccomplex_t IT_0555 = -IT_0531;
    const ccomplex_t IT_0556 = IT_0503 + IT_0555;
    const ccomplex_t IT_0557 = IT_0002*IT_0556;
    const ccomplex_t IT_0558 = (-4)*IT_0552;
    const ccomplex_t IT_0559 = IT_0557 + IT_0558;
    const ccomplex_t IT_0560 = m_C_1*IT_0559;
    const ccomplex_t IT_0561 = 2*IT_0560;
    const ccomplex_t IT_0562 = (-0.25)*IT_0241;
    const ccomplex_t IT_0563 = IT_0269*IT_0477 + IT_0307*IT_0483 + IT_0454
      *IT_0495 + IT_0481*IT_0500 + IT_0459*IT_0519 + IT_0510*IT_0526 + IT_0338
      *IT_0537 + IT_0482*IT_0540 + IT_0497*IT_0550 + IT_0511*IT_0554 + IT_0528
      *IT_0561 + IT_0214*IT_0562;
    const ccomplex_t IT_0564 = (-4)*IT_0243*(IT_0307*IT_0477 + IT_0269*IT_0483
       + IT_0500*IT_0510 + IT_0481*IT_0526 + IT_0280*IT_0537 + IT_0511*IT_0540 +
       IT_0521*IT_0550 + IT_0482*IT_0554 + IT_0495*IT_0561 + IT_0525*IT_0562);
    const ccomplex_t IT_0565 = conj(IT_0316)*IT_0432;
    const ccomplex_t IT_0566 = (-0.25)*IT_0565;
    const ccomplex_t IT_0567 = IT_0012*IT_0364;
    const ccomplex_t IT_0568 = 0.666666666666667*IT_0218 + (-0.166666666666667
      )*IT_0567;
    const ccomplex_t IT_0569 = s_14 + IT_0568;
    const ccomplex_t IT_0570 = m_C_2*IT_0569;
    const ccomplex_t IT_0571 = (-6)*IT_0570;
    const ccomplex_t IT_0572 = (-2)*IT_0571;
    const ccomplex_t IT_0573 = IT_0497*IT_0572;
    const ccomplex_t IT_0574 = IT_0466*IT_0519;
    const ccomplex_t IT_0575 = (-0.666666666666667)*IT_0129 +
       0.166666666666667*IT_0401;
    const ccomplex_t IT_0576 = s_24 + IT_0575;
    const ccomplex_t IT_0577 = m_C_1*IT_0576;
    const ccomplex_t IT_0578 = 6*IT_0577;
    const ccomplex_t IT_0579 = (-2)*IT_0578;
    const ccomplex_t IT_0580 = IT_0521*IT_0579;
    const ccomplex_t IT_0581 = conj(IT_0242)*IT_0459;
    const ccomplex_t IT_0582 = (-0.25)*IT_0581;
    const ccomplex_t IT_0583 = conj(IT_0117)*IT_0289;
    const ccomplex_t IT_0584 = (-0.25)*IT_0583;
    const ccomplex_t IT_0585 = conj(IT_0205)*IT_0518;
    const ccomplex_t IT_0586 = (-0.25)*IT_0585;
    const ccomplex_t IT_0587 = conj(IT_0314)*IT_0432;
    const ccomplex_t IT_0588 = (-0.25)*IT_0587;
    const ccomplex_t IT_0589 = IT_0212*IT_0394;
    const ccomplex_t IT_0590 = (-0.25)*IT_0589;
    const ccomplex_t IT_0591 = IT_0213*IT_0394;
    const ccomplex_t IT_0592 = (-0.25)*IT_0591;
    const ccomplex_t IT_0593 = IT_0405*IT_0523;
    const ccomplex_t IT_0594 = (-0.25)*IT_0593;
    const ccomplex_t IT_0595 = IT_0405*IT_0524;
    const ccomplex_t IT_0596 = (-0.25)*IT_0595;
    const ccomplex_t IT_0597 = conj(IT_0109)*IT_0362;
    const ccomplex_t IT_0598 = (-0.25)*IT_0597;
    const ccomplex_t IT_0599 = conj(IT_0184)*IT_0387;
    const ccomplex_t IT_0600 = (-0.25)*IT_0599;
    const ccomplex_t IT_0601 = IT_0497*IT_0579;
    const ccomplex_t IT_0602 = IT_0521*IT_0572;
    const ccomplex_t IT_0603 = conj(IT_0243)*IT_0459;
    const ccomplex_t IT_0604 = (-0.25)*IT_0603;
    const ccomplex_t IT_0605 = conj(IT_0109)*IT_0289;
    const ccomplex_t IT_0606 = (-0.25)*IT_0605;
    const ccomplex_t IT_0607 = conj(IT_0184)*IT_0518;
    const ccomplex_t IT_0608 = (-0.25)*IT_0607;
    const ccomplex_t IT_0609 = conj(IT_0314)*IT_0440;
    const ccomplex_t IT_0610 = (-0.25)*IT_0609;
    const ccomplex_t IT_0611 = conj(IT_0316)*IT_0440;
    const ccomplex_t IT_0612 = (-0.25)*IT_0611;
    const ccomplex_t IT_0613 = conj(IT_0273)*IT_0444;
    const ccomplex_t IT_0614 = (-0.25)*IT_0613;
    const ccomplex_t IT_0615 = conj(IT_0279)*IT_0444;
    const ccomplex_t IT_0616 = (-0.25)*IT_0615;
    const ccomplex_t IT_0617 = IT_0212*IT_0405;
    const ccomplex_t IT_0618 = (-0.25)*IT_0617;
    const ccomplex_t IT_0619 = IT_0213*IT_0405;
    const ccomplex_t IT_0620 = (-0.25)*IT_0619;
    const ccomplex_t IT_0621 = IT_0394*IT_0523;
    const ccomplex_t IT_0622 = (-0.25)*IT_0621;
    const ccomplex_t IT_0623 = IT_0394*IT_0524;
    const ccomplex_t IT_0624 = (-0.25)*IT_0623;
    const ccomplex_t IT_0625 = conj(IT_0117)*IT_0362;
    const ccomplex_t IT_0626 = (-0.25)*IT_0625;
    const ccomplex_t IT_0627 = conj(IT_0205)*IT_0387;
    const ccomplex_t IT_0628 = (-0.25)*IT_0627;
    const ccomplex_t IT_0629 = IT_0466*IT_0527 + IT_0432*IT_0528 + conj
      (IT_0242)*IT_0537 + IT_0601 + IT_0602 + IT_0604 + IT_0606 + IT_0608 +
       IT_0610 + IT_0612 + IT_0614 + IT_0616 + IT_0618 + IT_0620 + IT_0622 +
       IT_0624 + IT_0626 + IT_0628;
    const ccomplex_t IT_0630 = 0.25*IT_0018;
    const ccomplex_t IT_0631 = IT_0001 + IT_0630;
    const ccomplex_t IT_0632 = m_C_1*m_C_2*IT_0002;
    const ccomplex_t IT_0633 = (-0.25)*IT_0003;
    const ccomplex_t IT_0634 = IT_0382 + IT_0633;
    const ccomplex_t IT_0635 = IT_0632*IT_0634;
    const ccomplex_t IT_0636 = (-4)*IT_0635;
    const ccomplex_t IT_0637 = (-0.125)*IT_0636;
    const ccomplex_t IT_0638 = (-0.125)*conj(IT_0184);
    const ccomplex_t IT_0639 = m_C_1*IT_0002;
    const ccomplex_t IT_0640 = IT_0002*IT_0375;
    const ccomplex_t IT_0641 = (-2)*IT_0126 + (-0.5)*IT_0640;
    const ccomplex_t IT_0642 = IT_0121 + IT_0641;
    const ccomplex_t IT_0643 = IT_0639*IT_0642;
    const ccomplex_t IT_0644 = (-2)*IT_0643;
    const ccomplex_t IT_0645 = 2*IT_0644;
    const ccomplex_t IT_0646 = (-0.125)*IT_0317;
    const ccomplex_t IT_0647 = (-2)*IT_0208;
    const ccomplex_t IT_0648 = (-0.125)*conj(IT_0348);
    const ccomplex_t IT_0649 = (-0.125)*conj(IT_0205);
    const ccomplex_t IT_0650 = (-0.125)*IT_0280;
    const ccomplex_t IT_0651 = (-0.125)*IT_0338;
    const ccomplex_t IT_0652 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0653 = (-0.166666666666667)*IT_0368 + (
      -0.666666666666667)*IT_0652;
    const ccomplex_t IT_0654 = s_14 + IT_0653;
    const ccomplex_t IT_0655 = m_C_2*IT_0654;
    const ccomplex_t IT_0656 = (-6)*IT_0655;
    const ccomplex_t IT_0657 = 2*IT_0656;
    const ccomplex_t IT_0658 = (-0.125)*IT_0355;
    const ccomplex_t IT_0659 = IT_0002*IT_0462;
    const ccomplex_t IT_0660 = (-4)*IT_0001 + (-4)*IT_0019 + (-4)*IT_0021 + 4
      *IT_0659;
    const ccomplex_t IT_0661 = (-0.125)*conj(IT_0290);
    const ccomplex_t IT_0662 = (-0.125)*conj(IT_0117)*IT_0299 + (-0.125)*conj
      (IT_0243)*IT_0550 + IT_0088*IT_0631 + IT_0214*IT_0637 + IT_0496*IT_0638 +
       IT_0645*IT_0646 + IT_0647*IT_0648 + IT_0520*IT_0649 + IT_0579*IT_0650 +
       IT_0572*IT_0651 + IT_0657*IT_0658 + IT_0660*IT_0661;
    const ccomplex_t IT_0663 = (-0.125)*conj(IT_0109)*IT_0299 + (-0.125)*conj
      (IT_0242)*IT_0550 + IT_0059*IT_0631 + IT_0525*IT_0637 + IT_0520*IT_0638 +
       IT_0496*IT_0649 + IT_0572*IT_0650 + IT_0579*IT_0651 + IT_0646*IT_0657 +
       IT_0645*IT_0658 + IT_0648*IT_0660 + IT_0647*IT_0661;
    const ccomplex_t IT_0664 = IT_0327*IT_0379;
    const ccomplex_t IT_0665 = IT_0328*IT_0379;
    const ccomplex_t IT_0666 = IT_0244*IT_0372;
    const ccomplex_t IT_0667 = IT_0245*IT_0372;
    const ccomplex_t IT_0668 = conj(IT_0242)*IT_0561;
    const ccomplex_t IT_0669 = conj(IT_0117)*IT_0347;
    const ccomplex_t IT_0670 = conj(IT_0337)*IT_0432;
    const ccomplex_t IT_0671 = conj(IT_0273)*IT_0440;
    const ccomplex_t IT_0672 = conj(IT_0279)*IT_0440;
    const ccomplex_t IT_0673 = conj(IT_0314)*IT_0444;
    const ccomplex_t IT_0674 = conj(IT_0316)*IT_0444;
    const ccomplex_t IT_0675 = conj(IT_0243)*IT_0454;
    const ccomplex_t IT_0676 = conj(IT_0109)*IT_0325;
    const ccomplex_t IT_0677 = conj(IT_0184)*IT_0494;
    const ccomplex_t IT_0678 = conj(IT_0331)*IT_0432;
    const ccomplex_t IT_0679 = conj(IT_0348)*IT_0645;
    const ccomplex_t IT_0680 = IT_0355*IT_0422 + IT_0627 + conj(IT_0290)
      *IT_0657 + IT_0664 + IT_0665 + IT_0666 + IT_0667 + IT_0668 + IT_0669 +
       IT_0670 + IT_0671 + IT_0672 + IT_0673 + IT_0674 + IT_0675 + IT_0676 +
       IT_0677 + IT_0678 + IT_0679;
    const ccomplex_t IT_0681 = IT_0327*IT_0372;
    const ccomplex_t IT_0682 = IT_0328*IT_0372;
    const ccomplex_t IT_0683 = IT_0244*IT_0379;
    const ccomplex_t IT_0684 = IT_0245*IT_0379;
    const ccomplex_t IT_0685 = conj(IT_0243)*IT_0561;
    const ccomplex_t IT_0686 = conj(IT_0109)*IT_0347;
    const ccomplex_t IT_0687 = conj(IT_0242)*IT_0454;
    const ccomplex_t IT_0688 = conj(IT_0117)*IT_0325;
    const ccomplex_t IT_0689 = IT_0317*IT_0422;
    const ccomplex_t IT_0690 = conj(IT_0205)*IT_0494;
    const ccomplex_t IT_0691 = IT_0280*IT_0432;
    const ccomplex_t IT_0692 = IT_0314*(IT_0338*IT_0440 + IT_0355*IT_0444 +
       IT_0599 + conj(IT_0290)*IT_0645 + conj(IT_0348)*IT_0657 + IT_0681 +
       IT_0682 + IT_0683 + IT_0684 + IT_0685 + IT_0686 + IT_0687 + IT_0688 +
       IT_0689 + IT_0690 + IT_0691);
    const ccomplex_t IT_0693 = (-4)*IT_0022*((IT_0046 + IT_0058)*IT_0059 +
       IT_0087*IT_0088) + IT_0101*(conj(IT_0109)*IT_0112 + conj(IT_0117)*IT_0120
      ) + IT_0133*(IT_0112*conj(IT_0184) + IT_0120*conj(IT_0205)) + IT_0208*
      (IT_0088*IT_0211 + IT_0087*IT_0214) + IT_0226*(conj(IT_0205)*IT_0211 +
       conj(IT_0184)*IT_0229) + IT_0241*(IT_0211*conj(IT_0242) + IT_0229*conj
      (IT_0243)) + IT_0326 + IT_0109*IT_0363 + IT_0317*(IT_0373 + IT_0374 +
       IT_0380 + IT_0381 + IT_0388) + IT_0280*(IT_0388 + IT_0395 + IT_0396 +
       IT_0406 + IT_0407) + IT_0355*(IT_0408 + IT_0409 + IT_0410 + IT_0411 +
       IT_0412 + IT_0423 + IT_0433 + IT_0441 + IT_0442 + IT_0445 + IT_0455 +
       IT_0456) + IT_0338*(IT_0456 + IT_0460 + IT_0467 + IT_0468 + IT_0469 +
       IT_0470 + IT_0471 + IT_0472 + IT_0473 + IT_0474 + IT_0475) + (-4)*IT_0205
      *IT_0522 + (-4)*IT_0184*IT_0529 + (-4)*IT_0242*IT_0563 + IT_0564 + (-4)*
      (IT_0273 + IT_0279)*(conj(IT_0243)*IT_0537 + IT_0566 + IT_0573 + IT_0574 +
       IT_0580 + IT_0582 + IT_0584 + IT_0586 + IT_0588 + IT_0590 + IT_0592 +
       IT_0594 + IT_0596 + IT_0598 + IT_0600) + (-4)*IT_0337*(conj(IT_0242)
      *IT_0537 + IT_0601 + IT_0602 + IT_0604 + IT_0606 + IT_0608 + IT_0610 +
       IT_0612 + IT_0614 + IT_0616 + IT_0618 + IT_0620 + IT_0622 + IT_0624 +
       IT_0626 + IT_0628) + (-4)*IT_0331*IT_0629 + (-8)*conj(IT_0290)*(IT_0087
      *IT_0631 + IT_0211*IT_0637) + (-8)*conj(IT_0348)*((IT_0046 + IT_0058)
      *IT_0631 + IT_0229*IT_0637) + (-8)*IT_0290*IT_0662 + (-8)*IT_0348*IT_0663 
      + IT_0351*(IT_0627 + conj(IT_0290)*IT_0657 + IT_0664 + IT_0665 + IT_0666 +
       IT_0667 + IT_0668 + IT_0669 + IT_0670 + IT_0671 + IT_0672 + IT_0673 +
       IT_0674 + IT_0675 + IT_0676 + IT_0677 + IT_0678 + IT_0679) + IT_0354
      *IT_0680 + IT_0316*(IT_0599 + conj(IT_0290)*IT_0645 + conj(IT_0348)
      *IT_0657 + IT_0681 + IT_0682 + IT_0683 + IT_0684 + IT_0685 + IT_0686 +
       IT_0687 + IT_0688 + IT_0689 + IT_0690 + IT_0691) + IT_0692;
    return create_ccomplex_return(IT_0693);
}

