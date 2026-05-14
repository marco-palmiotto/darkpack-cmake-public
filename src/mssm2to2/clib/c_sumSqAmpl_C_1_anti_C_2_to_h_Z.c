#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_h_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_h_Z(
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
    const ccomplex_t IT_0000 = s_14*s_24;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_12*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = V_u1*U_Wm2;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = U_d2*V_Wp1;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (IT_0009 + -IT_0012)*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*e_em*IT_0015;
    const ccomplex_t IT_0017 = 0.5*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = IT_0013*IT_0020;
    const ccomplex_t IT_0022 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = IT_0014*IT_0019;
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0027 = IT_0024*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0023 + -IT_0025 + (-2)
      *IT_0027);
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0018*IT_0029;
    const ccomplex_t IT_0031 = pow(m_h, 2);
    const ccomplex_t IT_0032 = pow(m_C_1, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0031 + IT_0032 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0031 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0037 = V_u1*U_Wm1;
    const ccomplex_t IT_0038 = IT_0007*IT_0037;
    const ccomplex_t IT_0039 = U_d1*V_Wp1;
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = IT_0014*(IT_0038 + -IT_0040);
    const ccomplex_t IT_0042 = 1.4142135623731*e_em*IT_0041;
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0046 = IT_0021*IT_0045;
    const ccomplex_t IT_0047 = IT_0024*IT_0045;
    const ccomplex_t IT_0048 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0049 = IT_0024*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + -IT_0047 + (-2)
      *IT_0049);
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0044*IT_0051;
    const ccomplex_t IT_0053 = IT_0036*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0031 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0056 = V_u2*U_Wm2;
    const ccomplex_t IT_0057 = IT_0007*IT_0056;
    const ccomplex_t IT_0058 = U_d2*V_Wp2;
    const ccomplex_t IT_0059 = IT_0010*IT_0058;
    const ccomplex_t IT_0060 = IT_0014*(IT_0057 + -IT_0059);
    const ccomplex_t IT_0061 = 1.4142135623731*e_em*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0065 = IT_0021*IT_0064;
    const ccomplex_t IT_0066 = IT_0024*IT_0064;
    const ccomplex_t IT_0067 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0068 = IT_0024*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0065 + -IT_0066 + (-2)
      *IT_0068);
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0063*IT_0070;
    const ccomplex_t IT_0072 = IT_0055*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0075 = IT_0021*IT_0074;
    const ccomplex_t IT_0076 = IT_0024*IT_0074;
    const ccomplex_t IT_0077 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0078 = IT_0024*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + -IT_0076 + (-2)
      *IT_0078);
    const ccomplex_t IT_0080 = (-0.5)*IT_0079;
    const ccomplex_t IT_0081 = IT_0018*IT_0080;
    const ccomplex_t IT_0082 = pow(m_C_2, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_23 + IT_0031 + IT_0082 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = conj(IT_0073) + conj(IT_0085);
    const ccomplex_t IT_0087 = conj(IT_0035) + conj(IT_0054);
    const ccomplex_t IT_0088 = pow(m_Z, 2);
    const ccomplex_t IT_0089 = s_12*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0004;
    const ccomplex_t IT_0091 = (-2)*IT_0000;
    const ccomplex_t IT_0092 = IT_0089 + IT_0090 + IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0094 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0095 = IT_0007*IT_0094;
    const ccomplex_t IT_0096 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0097 = IT_0010*IT_0096;
    const ccomplex_t IT_0098 = IT_0014*(IT_0095 + -IT_0097);
    const ccomplex_t IT_0099 = 1.4142135623731*e_em*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = IT_0051*IT_0101;
    const ccomplex_t IT_0103 = IT_0036*IT_0093*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0105 = IT_0051*IT_0063;
    const ccomplex_t IT_0106 = IT_0055*IT_0104*IT_0105;
    const ccomplex_t IT_0107 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0108 = IT_0021*IT_0107;
    const ccomplex_t IT_0109 = IT_0024*IT_0107;
    const ccomplex_t IT_0110 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0111 = IT_0024*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*(IT_0108 + -IT_0109 + (-2)
      *IT_0111);
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0018*IT_0113;
    const ccomplex_t IT_0115 = IT_0083*IT_0093*IT_0114;
    const ccomplex_t IT_0116 = m_C_1*IT_0085;
    const ccomplex_t IT_0117 = m_W*e_em;
    const ccomplex_t IT_0118 = cpow(IT_0019, -4);
    const ccomplex_t IT_0119 = sin(beta);
    const ccomplex_t IT_0120 = IT_0007*IT_0013*IT_0118*IT_0119;
    const ccomplex_t IT_0121 = cos(beta);
    const ccomplex_t IT_0122 = IT_0010*IT_0013*IT_0118*IT_0121;
    const ccomplex_t IT_0123 = cpow(IT_0013, -3);
    const ccomplex_t IT_0124 = IT_0007*IT_0119*IT_0123;
    const ccomplex_t IT_0125 = IT_0010*IT_0121*IT_0123;
    const ccomplex_t IT_0126 = cpow(IT_0019, -2);
    const ccomplex_t IT_0127 = IT_0007*IT_0014*IT_0119*IT_0126;
    const ccomplex_t IT_0128 = IT_0010*IT_0014*IT_0121*IT_0126;
    const ccomplex_t IT_0129 = cpow(IT_0013, -2);
    const ccomplex_t IT_0130 = cpow(IT_0126 + IT_0129, -1);
    const ccomplex_t IT_0131 = (IT_0120 + -IT_0122 + IT_0124 + -IT_0125 + 2
      *IT_0127 + (-2)*IT_0128)*IT_0130;
    const ccomplex_t IT_0132 = IT_0117*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0051*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0032 + -IT_0082 + -reg_prop, -1);
    const ccomplex_t IT_0136 = IT_0134*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = IT_0070*IT_0101;
    const ccomplex_t IT_0139 = IT_0036*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_C_2*IT_0140;
    const ccomplex_t IT_0142 = m_C_1*IT_0073;
    const ccomplex_t IT_0143 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0144 = IT_0007*IT_0143;
    const ccomplex_t IT_0145 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0146 = IT_0010*IT_0145;
    const ccomplex_t IT_0147 = IT_0014*(IT_0144 + -IT_0146);
    const ccomplex_t IT_0148 = 1.4142135623731*e_em*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = IT_0029*IT_0150;
    const ccomplex_t IT_0152 = IT_0033*IT_0104*IT_0151;
    const ccomplex_t IT_0153 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0154 = IT_0021*IT_0153;
    const ccomplex_t IT_0155 = IT_0024*IT_0153;
    const ccomplex_t IT_0156 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0157 = IT_0024*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + -IT_0155 + (-2)
      *IT_0157);
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0150*IT_0159;
    const ccomplex_t IT_0161 = IT_0033*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = m_C_2*IT_0162;
    const ccomplex_t IT_0164 = IT_0103 + IT_0106 + IT_0115 + IT_0116 + 
      -IT_0137 + IT_0141 + IT_0142 + IT_0152 + IT_0163;
    const ccomplex_t IT_0165 = s_24*IT_0088;
    const ccomplex_t IT_0166 = IT_0001*IT_0165;
    const ccomplex_t IT_0167 = (-0.25)*IT_0166;
    const ccomplex_t IT_0168 = s_24 + IT_0167;
    const ccomplex_t IT_0169 = m_C_1*IT_0168;
    const ccomplex_t IT_0170 = (-4)*IT_0169;
    const ccomplex_t IT_0171 = 2*IT_0085;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = IT_0073 + IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0166;
    const ccomplex_t IT_0175 = s_24 + IT_0174;
    const ccomplex_t IT_0176 = m_C_1*IT_0175;
    const ccomplex_t IT_0177 = 2*IT_0176;
    const ccomplex_t IT_0178 = 2*IT_0035;
    const ccomplex_t IT_0179 = 2*IT_0054;
    const ccomplex_t IT_0180 = IT_0178 + IT_0179;
    const ccomplex_t IT_0181 = IT_0018*IT_0159;
    const ccomplex_t IT_0182 = IT_0033*IT_0104*IT_0181;
    const ccomplex_t IT_0183 = m_C_2*IT_0035;
    const ccomplex_t IT_0184 = IT_0044*IT_0070;
    const ccomplex_t IT_0185 = IT_0036*IT_0093*IT_0184;
    const ccomplex_t IT_0186 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0187 = IT_0007*IT_0186;
    const ccomplex_t IT_0188 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0189 = IT_0010*IT_0188;
    const ccomplex_t IT_0190 = IT_0014*(IT_0187 + -IT_0189);
    const ccomplex_t IT_0191 = 1.4142135623731*e_em*IT_0190;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = IT_0070*IT_0193;
    const ccomplex_t IT_0195 = IT_0055*IT_0104*IT_0194;
    const ccomplex_t IT_0196 = IT_0070*IT_0133;
    const ccomplex_t IT_0197 = IT_0135*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = IT_0080*IT_0150;
    const ccomplex_t IT_0200 = IT_0083*IT_0093*IT_0199;
    const ccomplex_t IT_0201 = IT_0113*IT_0150;
    const ccomplex_t IT_0202 = IT_0083*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = m_C_1*IT_0203;
    const ccomplex_t IT_0205 = m_C_2*IT_0054;
    const ccomplex_t IT_0206 = IT_0051*IT_0193;
    const ccomplex_t IT_0207 = IT_0055*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = m_C_1*IT_0208;
    const ccomplex_t IT_0210 = -IT_0182 + -IT_0183 + -IT_0185 + -IT_0195 +
       IT_0198 + -IT_0200 + -IT_0204 + -IT_0205 + -IT_0209;
    const ccomplex_t IT_0211 = s_14*IT_0088;
    const ccomplex_t IT_0212 = IT_0001*IT_0211;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = s_14 + IT_0213;
    const ccomplex_t IT_0215 = m_C_2*IT_0214;
    const ccomplex_t IT_0216 = 2*IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0073;
    const ccomplex_t IT_0218 = (-2)*IT_0085;
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = (-0.25)*IT_0212;
    const ccomplex_t IT_0221 = s_14 + IT_0220;
    const ccomplex_t IT_0222 = m_C_2*IT_0221;
    const ccomplex_t IT_0223 = (-4)*IT_0222;
    const ccomplex_t IT_0224 = (-2)*IT_0035;
    const ccomplex_t IT_0225 = (-2)*IT_0054;
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = m_C_1*m_C_2;
    const ccomplex_t IT_0228 = IT_0088*IT_0227;
    const ccomplex_t IT_0229 = (-3)*IT_0228;
    const ccomplex_t IT_0230 = -IT_0203 + -IT_0208;
    const ccomplex_t IT_0231 = -IT_0140 + -IT_0162;
    const ccomplex_t IT_0232 = 3*IT_0228;
    const ccomplex_t IT_0233 = e_em*IT_0014*IT_0020*(IT_0010*IT_0119 + IT_0007
      *IT_0121);
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0008*IT_0121;
    const ccomplex_t IT_0236 = IT_0011*IT_0119;
    const ccomplex_t IT_0237 = IT_0014*(IT_0235 + IT_0236);
    const ccomplex_t IT_0238 = 1.4142135623731*e_em*IT_0237;
    const ccomplex_t IT_0239 = 0.5*IT_0238;
    const ccomplex_t IT_0240 = cpow((-2)*s_12 + -IT_0032 + -IT_0082 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0241 = IT_0239*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = IT_0234*IT_0242;
    const ccomplex_t IT_0244 = e_em*IT_0014*IT_0020*(IT_0007*IT_0119 + 
      -IT_0010*IT_0121);
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = IT_0011*IT_0121;
    const ccomplex_t IT_0247 = IT_0008*IT_0119;
    const ccomplex_t IT_0248 = IT_0014*(IT_0246 + -IT_0247);
    const ccomplex_t IT_0249 = 1.4142135623731*e_em*IT_0248;
    const ccomplex_t IT_0250 = (-0.5)*IT_0249;
    const ccomplex_t IT_0251 = cpow((-2)*s_12 + -IT_0032 + -IT_0082 + IT_0088 
      + -reg_prop, -1);
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = IT_0245*IT_0253;
    const ccomplex_t IT_0255 = (-2)*IT_0035 + (-2)*IT_0054 + -IT_0243 + 
      -IT_0254;
    const ccomplex_t IT_0256 = s_12*s_24;
    const ccomplex_t IT_0257 = s_14*IT_0082;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = IT_0256 + IT_0258;
    const ccomplex_t IT_0260 = (-2)*IT_0073 + (-2)*IT_0085 + -IT_0243 + 
      -IT_0254;
    const ccomplex_t IT_0261 = s_12*s_14;
    const ccomplex_t IT_0262 = s_24*IT_0032;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = IT_0261 + IT_0263;
    const ccomplex_t IT_0265 = IT_0121*IT_0145;
    const ccomplex_t IT_0266 = IT_0119*IT_0143;
    const ccomplex_t IT_0267 = IT_0014*(IT_0265 + -IT_0266);
    const ccomplex_t IT_0268 = 1.4142135623731*e_em*IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = IT_0251*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0245*IT_0271;
    const ccomplex_t IT_0273 = IT_0121*IT_0143;
    const ccomplex_t IT_0274 = IT_0119*IT_0145;
    const ccomplex_t IT_0275 = IT_0014*(IT_0273 + IT_0274);
    const ccomplex_t IT_0276 = 1.4142135623731*e_em*IT_0275;
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = IT_0240*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0234*IT_0279;
    const ccomplex_t IT_0281 = 2*IT_0140 + 2*IT_0162 + -IT_0272 + -IT_0280;
    const ccomplex_t IT_0282 = 2*IT_0259;
    const ccomplex_t IT_0283 = pow(s_24, 2);
    const ccomplex_t IT_0284 = IT_0001*IT_0283;
    const ccomplex_t IT_0285 = -IT_0082;
    const ccomplex_t IT_0286 = IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = s_12*IT_0286;
    const ccomplex_t IT_0288 = 2*IT_0287;
    const ccomplex_t IT_0289 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = s_12 + IT_0290;
    const ccomplex_t IT_0292 = IT_0227*IT_0291;
    const ccomplex_t IT_0293 = 2*IT_0292;
    const ccomplex_t IT_0294 = m_C_2*IT_0291;
    const ccomplex_t IT_0295 = 2*IT_0294;
    const ccomplex_t IT_0296 = m_C_1*IT_0286;
    const ccomplex_t IT_0297 = 2*IT_0296;
    const ccomplex_t IT_0298 = (-2)*IT_0259;
    const ccomplex_t IT_0299 = 2*IT_0203 + 2*IT_0208 + -IT_0272 + -IT_0280;
    const ccomplex_t IT_0300 = s_12*IT_0291;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = IT_0227*IT_0286;
    const ccomplex_t IT_0303 = (-2)*IT_0302;
    const ccomplex_t IT_0304 = 2*conj(IT_0035);
    const ccomplex_t IT_0305 = 2*conj(IT_0054);
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = (-2)*conj(IT_0073);
    const ccomplex_t IT_0308 = (-2)*conj(IT_0085);
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = (-2)*conj(IT_0035);
    const ccomplex_t IT_0311 = (-2)*conj(IT_0054);
    const ccomplex_t IT_0312 = IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = 2*IT_0264;
    const ccomplex_t IT_0314 = 2*IT_0170;
    const ccomplex_t IT_0315 = (-8)*IT_0006;
    const ccomplex_t IT_0316 = (-4)*IT_0092;
    const ccomplex_t IT_0317 = (-2)*IT_0216;
    const ccomplex_t IT_0318 = 2*IT_0177;
    const ccomplex_t IT_0319 = (-2)*IT_0264;
    const ccomplex_t IT_0320 = (-2)*IT_0223;
    const ccomplex_t IT_0321 = IT_0089 + IT_0091;
    const ccomplex_t IT_0322 = IT_0001*IT_0321;
    const ccomplex_t IT_0323 = 4*s_12 + (-2)*IT_0322;
    const ccomplex_t IT_0324 = 6*IT_0227;
    const ccomplex_t IT_0325 = pow(s_14, 2);
    const ccomplex_t IT_0326 = IT_0001*IT_0325;
    const ccomplex_t IT_0327 = -IT_0032;
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = m_C_2*IT_0328;
    const ccomplex_t IT_0330 = (-2)*IT_0329;
    const ccomplex_t IT_0331 = m_C_1*IT_0291;
    const ccomplex_t IT_0332 = (-2)*IT_0331;
    const ccomplex_t IT_0333 = IT_0210*(conj(IT_0255)*IT_0295 + conj(IT_0281)
      *IT_0297 + IT_0216*IT_0309 + IT_0223*IT_0312 + conj(IT_0230)*IT_0314 +
       conj(IT_0231)*IT_0318 + conj(IT_0210)*IT_0323 + conj(IT_0164)*IT_0324 +
       conj(IT_0260)*IT_0330 + conj(IT_0299)*IT_0332);
    const ccomplex_t IT_0334 = s_12*IT_0328;
    const ccomplex_t IT_0335 = 2*IT_0334;
    const ccomplex_t IT_0336 = IT_0227*IT_0328;
    const ccomplex_t IT_0337 = (-2)*IT_0336;
    const ccomplex_t IT_0338 = 2*conj(IT_0085);
    const ccomplex_t IT_0339 = 0.5*IT_0338;
    const ccomplex_t IT_0340 = conj(IT_0073) + IT_0339;
    const ccomplex_t IT_0341 = 0.5*conj(IT_0260);
    const ccomplex_t IT_0342 = 0.5*conj(IT_0255);
    const ccomplex_t IT_0343 = 0.5*conj(IT_0299);
    const ccomplex_t IT_0344 = 0.5*conj(IT_0281);
    const ccomplex_t IT_0345 = 0.5*conj(IT_0210);
    const ccomplex_t IT_0346 = 0.5*conj(IT_0164);
    const ccomplex_t IT_0347 = 0.5*IT_0177*IT_0306 + 0.5*conj(IT_0230)*IT_0317
       + 0.5*conj(IT_0231)*IT_0320 + IT_0170*IT_0340 + IT_0332*IT_0341 + IT_0297
      *IT_0342 + IT_0330*IT_0343 + IT_0295*IT_0344 + IT_0324*IT_0345 + IT_0323
      *IT_0346;
    const ccomplex_t IT_0348 = (-8)*IT_0006*((IT_0035 + IT_0054)*IT_0086 + 
      (IT_0073 + IT_0085)*IT_0087) + (-4)*((IT_0073 + IT_0085)*IT_0086 + 
      (IT_0035 + IT_0054)*IT_0087)*IT_0092 + 2*conj(IT_0164)*(IT_0170*IT_0173 +
       0.5*IT_0177*IT_0180) + conj(IT_0210)*(IT_0216*IT_0219 + IT_0223*IT_0226) 
      + IT_0229*(IT_0219*conj(IT_0230) + IT_0226*conj(IT_0231)) + (IT_0226*conj
      (IT_0230) + IT_0219*conj(IT_0231))*IT_0232 + (IT_0180 + IT_0217 + IT_0218)
      *conj(IT_0255)*IT_0259 + 2*(IT_0173 + 0.5*IT_0224 + 0.5*IT_0225)*conj
      (IT_0260)*IT_0264 + IT_0281*(conj(IT_0231)*IT_0282 + conj(IT_0281)*IT_0288
       + conj(IT_0260)*IT_0293 + conj(IT_0164)*IT_0295 + conj(IT_0210)*IT_0297 +
       conj(IT_0230)*IT_0298 + conj(IT_0299)*IT_0301 + conj(IT_0255)*IT_0303) +
       IT_0255*(conj(IT_0255)*IT_0288 + conj(IT_0210)*IT_0295 + conj(IT_0164)
      *IT_0297 + IT_0293*conj(IT_0299) + conj(IT_0260)*IT_0301 + conj(IT_0281)
      *IT_0303 + IT_0259*(IT_0306 + IT_0307 + IT_0308)) + IT_0230*(conj(IT_0281)
      *IT_0298 + IT_0229*IT_0309 + IT_0232*IT_0312 + conj(IT_0299)*IT_0313 +
       conj(IT_0210)*IT_0314 + conj(IT_0231)*IT_0315 + conj(IT_0230)*IT_0316 +
       conj(IT_0164)*IT_0317) + IT_0231*(conj(IT_0281)*IT_0282 + IT_0232*IT_0309
       + IT_0229*IT_0312 + conj(IT_0230)*IT_0315 + conj(IT_0231)*IT_0316 + conj
      (IT_0210)*IT_0318 + conj(IT_0299)*IT_0319 + conj(IT_0164)*IT_0320) +
       IT_0333 + IT_0299*(conj(IT_0255)*IT_0293 + conj(IT_0281)*IT_0301 + conj
      (IT_0230)*IT_0313 + conj(IT_0231)*IT_0319 + conj(IT_0164)*IT_0330 + conj
      (IT_0210)*IT_0332 + conj(IT_0299)*IT_0335 + conj(IT_0260)*IT_0337) + 2
      *IT_0260*(IT_0264*(0.5*IT_0310 + 0.5*IT_0311 + IT_0340) + IT_0335*IT_0341 
      + IT_0301*IT_0342 + IT_0337*IT_0343 + IT_0293*IT_0344 + IT_0330*IT_0345 +
       IT_0332*IT_0346) + 2*IT_0164*IT_0347;
    return create_ccomplex_return(IT_0348);
}

