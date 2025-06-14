#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_H0_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_H0_Z(
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
    const ccomplex_t IT_0000 = s_14*s_24;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_12*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = U_d2*V_Wp1;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = V_u1*U_Wm2;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (IT_0009 + IT_0012)*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*e_em*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = 0.5*IT_0017;
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
    const ccomplex_t IT_0031 = pow(m_H0, 2);
    const ccomplex_t IT_0032 = pow(m_C_1, 2);
    const ccomplex_t IT_0033 = cpow((-2)*s_13 + IT_0031 + IT_0032 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0030*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = cpow((-2)*s_13 + IT_0031 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0037 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0038 = IT_0021*IT_0037;
    const ccomplex_t IT_0039 = IT_0024*IT_0037;
    const ccomplex_t IT_0040 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0041 = IT_0024*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + -IT_0039 + (-2)
      *IT_0041);
    const ccomplex_t IT_0043 = 0.5*IT_0042;
    const ccomplex_t IT_0044 = U_d1*V_Wp1;
    const ccomplex_t IT_0045 = IT_0007*IT_0044;
    const ccomplex_t IT_0046 = V_u1*U_Wm1;
    const ccomplex_t IT_0047 = IT_0010*IT_0046;
    const ccomplex_t IT_0048 = IT_0014*(IT_0045 + IT_0047);
    const ccomplex_t IT_0049 = 1.4142135623731*e_em*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = 0.5*IT_0050;
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = IT_0036*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0031 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0056 = U_d2*V_Wp2;
    const ccomplex_t IT_0057 = IT_0007*IT_0056;
    const ccomplex_t IT_0058 = V_u2*U_Wm2;
    const ccomplex_t IT_0059 = IT_0010*IT_0058;
    const ccomplex_t IT_0060 = IT_0014*(IT_0057 + IT_0059);
    const ccomplex_t IT_0061 = 1.4142135623731*e_em*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
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
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0094 = IT_0043*IT_0063;
    const ccomplex_t IT_0095 = IT_0055*IT_0093*IT_0094;
    const ccomplex_t IT_0096 = cos(beta);
    const ccomplex_t IT_0097 = cpow(IT_0019, -4);
    const ccomplex_t IT_0098 = IT_0007*IT_0013*IT_0096*IT_0097;
    const ccomplex_t IT_0099 = sin(beta);
    const ccomplex_t IT_0100 = IT_0010*IT_0013*IT_0097*IT_0099;
    const ccomplex_t IT_0101 = cpow(IT_0013, -3);
    const ccomplex_t IT_0102 = IT_0007*IT_0096*IT_0101;
    const ccomplex_t IT_0103 = IT_0010*IT_0099*IT_0101;
    const ccomplex_t IT_0104 = cpow(IT_0019, -2);
    const ccomplex_t IT_0105 = IT_0007*IT_0014*IT_0096*IT_0104;
    const ccomplex_t IT_0106 = IT_0010*IT_0014*IT_0099*IT_0104;
    const ccomplex_t IT_0107 = cpow(IT_0013, -2);
    const ccomplex_t IT_0108 = cpow(IT_0104 + IT_0107, -1);
    const ccomplex_t IT_0109 = (IT_0098 + IT_0100 + IT_0102 + IT_0103 + 2
      *IT_0105 + 2*IT_0106)*IT_0108;
    const ccomplex_t IT_0110 = m_W*e_em;
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0043*IT_0112;
    const ccomplex_t IT_0114 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0032 + -IT_0082 + -reg_prop, -1);
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0118 = IT_0021*IT_0117;
    const ccomplex_t IT_0119 = IT_0024*IT_0117;
    const ccomplex_t IT_0120 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0121 = IT_0024*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + -IT_0119 + (-2)
      *IT_0121);
    const ccomplex_t IT_0123 = 0.5*IT_0122;
    const ccomplex_t IT_0124 = IT_0018*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0126 = IT_0083*IT_0124*IT_0125;
    const ccomplex_t IT_0127 = m_C_1*IT_0085;
    const ccomplex_t IT_0128 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0129 = IT_0007*IT_0128;
    const ccomplex_t IT_0130 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0131 = IT_0010*IT_0130;
    const ccomplex_t IT_0132 = IT_0014*(IT_0129 + IT_0131);
    const ccomplex_t IT_0133 = 1.4142135623731*e_em*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = 0.5*IT_0134;
    const ccomplex_t IT_0136 = IT_0043*IT_0135;
    const ccomplex_t IT_0137 = IT_0036*IT_0125*IT_0136;
    const ccomplex_t IT_0138 = IT_0070*IT_0135;
    const ccomplex_t IT_0139 = IT_0036*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_C_2*IT_0140;
    const ccomplex_t IT_0142 = m_C_1*IT_0073;
    const ccomplex_t IT_0143 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0144 = IT_0007*IT_0143;
    const ccomplex_t IT_0145 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0146 = IT_0010*IT_0145;
    const ccomplex_t IT_0147 = IT_0014*(IT_0144 + IT_0146);
    const ccomplex_t IT_0148 = 1.4142135623731*e_em*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0029*IT_0150;
    const ccomplex_t IT_0152 = IT_0033*IT_0093*IT_0151;
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
    const ccomplex_t IT_0164 = IT_0095 + -IT_0116 + IT_0126 + IT_0127 +
       IT_0137 + IT_0141 + IT_0142 + IT_0152 + IT_0163;
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
    const ccomplex_t IT_0181 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0182 = IT_0007*IT_0181;
    const ccomplex_t IT_0183 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0184 = IT_0010*IT_0183;
    const ccomplex_t IT_0185 = IT_0014*(IT_0182 + IT_0184);
    const ccomplex_t IT_0186 = 1.4142135623731*e_em*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = IT_0070*IT_0188;
    const ccomplex_t IT_0190 = IT_0055*IT_0093*IT_0189;
    const ccomplex_t IT_0191 = IT_0080*IT_0150;
    const ccomplex_t IT_0192 = IT_0083*IT_0125*IT_0191;
    const ccomplex_t IT_0193 = IT_0051*IT_0070;
    const ccomplex_t IT_0194 = IT_0036*IT_0125*IT_0193;
    const ccomplex_t IT_0195 = IT_0070*IT_0112;
    const ccomplex_t IT_0196 = IT_0114*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = IT_0123*IT_0150;
    const ccomplex_t IT_0199 = IT_0083*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = m_C_1*IT_0200;
    const ccomplex_t IT_0202 = IT_0018*IT_0159;
    const ccomplex_t IT_0203 = IT_0033*IT_0093*IT_0202;
    const ccomplex_t IT_0204 = m_C_2*IT_0035;
    const ccomplex_t IT_0205 = m_C_2*IT_0054;
    const ccomplex_t IT_0206 = IT_0043*IT_0188;
    const ccomplex_t IT_0207 = IT_0055*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = m_C_1*IT_0208;
    const ccomplex_t IT_0210 = -IT_0190 + -IT_0192 + -IT_0194 + IT_0197 + 
      -IT_0201 + -IT_0203 + -IT_0204 + -IT_0205 + -IT_0209;
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
    const ccomplex_t IT_0230 = -IT_0200 + -IT_0208;
    const ccomplex_t IT_0231 = -IT_0140 + -IT_0162;
    const ccomplex_t IT_0232 = 3*IT_0228;
    const ccomplex_t IT_0233 = -e_em*IT_0014*IT_0020*(IT_0010*IT_0096 + 
      -IT_0007*IT_0099);
    const ccomplex_t IT_0234 = IT_0011*IT_0096;
    const ccomplex_t IT_0235 = IT_0008*IT_0099;
    const ccomplex_t IT_0236 = IT_0014*(IT_0234 + IT_0235);
    const ccomplex_t IT_0237 = 1.4142135623731*e_em*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0237;
    const ccomplex_t IT_0239 = cpow((-2)*s_12 + -IT_0032 + -IT_0082 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0240 = IT_0238*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = IT_0233*IT_0241;
    const ccomplex_t IT_0243 = e_em*IT_0014*IT_0020*(IT_0007*IT_0096 + IT_0010
      *IT_0099);
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = IT_0008*IT_0096;
    const ccomplex_t IT_0246 = IT_0011*IT_0099;
    const ccomplex_t IT_0247 = IT_0014*(IT_0245 + -IT_0246);
    const ccomplex_t IT_0248 = 1.4142135623731*e_em*IT_0247;
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = cpow((-2)*s_12 + -IT_0032 + -IT_0082 + IT_0088 
      + -reg_prop, -1);
    const ccomplex_t IT_0251 = IT_0249*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = IT_0244*IT_0252;
    const ccomplex_t IT_0254 = (-2)*IT_0035 + (-2)*IT_0054 + -IT_0242 + 
      -IT_0253;
    const ccomplex_t IT_0255 = s_12*s_24;
    const ccomplex_t IT_0256 = s_14*IT_0082;
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = (-2)*IT_0073 + (-2)*IT_0085 + -IT_0242 + 
      -IT_0253;
    const ccomplex_t IT_0260 = s_12*s_14;
    const ccomplex_t IT_0261 = s_24*IT_0032;
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = IT_0096*IT_0143;
    const ccomplex_t IT_0265 = IT_0099*IT_0145;
    const ccomplex_t IT_0266 = IT_0014*(IT_0264 + -IT_0265);
    const ccomplex_t IT_0267 = 1.4142135623731*e_em*IT_0266;
    const ccomplex_t IT_0268 = 0.5*IT_0267;
    const ccomplex_t IT_0269 = IT_0250*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = IT_0244*IT_0270;
    const ccomplex_t IT_0272 = IT_0096*IT_0145;
    const ccomplex_t IT_0273 = IT_0099*IT_0143;
    const ccomplex_t IT_0274 = IT_0014*(IT_0272 + IT_0273);
    const ccomplex_t IT_0275 = 1.4142135623731*e_em*IT_0274;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0239*IT_0276;
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*IT_0277;
    const ccomplex_t IT_0279 = IT_0233*IT_0278;
    const ccomplex_t IT_0280 = 2*IT_0140 + 2*IT_0162 + -IT_0271 + -IT_0279;
    const ccomplex_t IT_0281 = 2*IT_0258;
    const ccomplex_t IT_0282 = pow(s_24, 2);
    const ccomplex_t IT_0283 = IT_0001*IT_0282;
    const ccomplex_t IT_0284 = -IT_0082;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = s_12*IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0286;
    const ccomplex_t IT_0288 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = s_12 + IT_0289;
    const ccomplex_t IT_0291 = IT_0227*IT_0290;
    const ccomplex_t IT_0292 = 2*IT_0291;
    const ccomplex_t IT_0293 = m_C_2*IT_0290;
    const ccomplex_t IT_0294 = 2*IT_0293;
    const ccomplex_t IT_0295 = m_C_1*IT_0285;
    const ccomplex_t IT_0296 = 2*IT_0295;
    const ccomplex_t IT_0297 = (-2)*IT_0258;
    const ccomplex_t IT_0298 = 2*IT_0200 + 2*IT_0208 + -IT_0271 + -IT_0279;
    const ccomplex_t IT_0299 = s_12*IT_0290;
    const ccomplex_t IT_0300 = (-2)*IT_0299;
    const ccomplex_t IT_0301 = IT_0227*IT_0285;
    const ccomplex_t IT_0302 = (-2)*IT_0301;
    const ccomplex_t IT_0303 = 2*conj(IT_0035);
    const ccomplex_t IT_0304 = 2*conj(IT_0054);
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = (-2)*conj(IT_0073);
    const ccomplex_t IT_0307 = (-2)*conj(IT_0085);
    const ccomplex_t IT_0308 = IT_0306 + IT_0307;
    const ccomplex_t IT_0309 = (-2)*conj(IT_0035);
    const ccomplex_t IT_0310 = (-2)*conj(IT_0054);
    const ccomplex_t IT_0311 = IT_0309 + IT_0310;
    const ccomplex_t IT_0312 = 2*IT_0263;
    const ccomplex_t IT_0313 = 2*IT_0170;
    const ccomplex_t IT_0314 = (-8)*IT_0006;
    const ccomplex_t IT_0315 = (-4)*IT_0092;
    const ccomplex_t IT_0316 = (-2)*IT_0216;
    const ccomplex_t IT_0317 = 2*IT_0177;
    const ccomplex_t IT_0318 = (-2)*IT_0263;
    const ccomplex_t IT_0319 = (-2)*IT_0223;
    const ccomplex_t IT_0320 = IT_0089 + IT_0091;
    const ccomplex_t IT_0321 = IT_0001*IT_0320;
    const ccomplex_t IT_0322 = 4*s_12 + (-2)*IT_0321;
    const ccomplex_t IT_0323 = 6*IT_0227;
    const ccomplex_t IT_0324 = pow(s_14, 2);
    const ccomplex_t IT_0325 = IT_0001*IT_0324;
    const ccomplex_t IT_0326 = -IT_0032;
    const ccomplex_t IT_0327 = IT_0325 + IT_0326;
    const ccomplex_t IT_0328 = m_C_2*IT_0327;
    const ccomplex_t IT_0329 = (-2)*IT_0328;
    const ccomplex_t IT_0330 = m_C_1*IT_0290;
    const ccomplex_t IT_0331 = (-2)*IT_0330;
    const ccomplex_t IT_0332 = IT_0210*(conj(IT_0254)*IT_0294 + conj(IT_0280)
      *IT_0296 + IT_0216*IT_0308 + IT_0223*IT_0311 + conj(IT_0230)*IT_0313 +
       conj(IT_0231)*IT_0317 + conj(IT_0210)*IT_0322 + conj(IT_0164)*IT_0323 +
       conj(IT_0259)*IT_0329 + conj(IT_0298)*IT_0331);
    const ccomplex_t IT_0333 = s_12*IT_0327;
    const ccomplex_t IT_0334 = 2*IT_0333;
    const ccomplex_t IT_0335 = IT_0227*IT_0327;
    const ccomplex_t IT_0336 = (-2)*IT_0335;
    const ccomplex_t IT_0337 = 2*conj(IT_0085);
    const ccomplex_t IT_0338 = 0.5*IT_0337;
    const ccomplex_t IT_0339 = conj(IT_0073) + IT_0338;
    const ccomplex_t IT_0340 = 0.5*conj(IT_0259);
    const ccomplex_t IT_0341 = 0.5*conj(IT_0254);
    const ccomplex_t IT_0342 = 0.5*conj(IT_0298);
    const ccomplex_t IT_0343 = 0.5*conj(IT_0280);
    const ccomplex_t IT_0344 = 0.5*conj(IT_0210);
    const ccomplex_t IT_0345 = 0.5*conj(IT_0164);
    const ccomplex_t IT_0346 = 0.5*IT_0177*IT_0305 + 0.5*conj(IT_0230)*IT_0316
       + 0.5*conj(IT_0231)*IT_0319 + IT_0170*IT_0339 + IT_0331*IT_0340 + IT_0296
      *IT_0341 + IT_0329*IT_0342 + IT_0294*IT_0343 + IT_0323*IT_0344 + IT_0322
      *IT_0345;
    const ccomplex_t IT_0347 = (-8)*IT_0006*((IT_0035 + IT_0054)*IT_0086 + 
      (IT_0073 + IT_0085)*IT_0087) + (-4)*((IT_0073 + IT_0085)*IT_0086 + 
      (IT_0035 + IT_0054)*IT_0087)*IT_0092 + 2*conj(IT_0164)*(IT_0170*IT_0173 +
       0.5*IT_0177*IT_0180) + conj(IT_0210)*(IT_0216*IT_0219 + IT_0223*IT_0226) 
      + IT_0229*(IT_0219*conj(IT_0230) + IT_0226*conj(IT_0231)) + (IT_0226*conj
      (IT_0230) + IT_0219*conj(IT_0231))*IT_0232 + (IT_0180 + IT_0217 + IT_0218)
      *conj(IT_0254)*IT_0258 + 2*(IT_0173 + 0.5*IT_0224 + 0.5*IT_0225)*conj
      (IT_0259)*IT_0263 + IT_0280*(conj(IT_0231)*IT_0281 + conj(IT_0280)*IT_0287
       + conj(IT_0259)*IT_0292 + conj(IT_0164)*IT_0294 + conj(IT_0210)*IT_0296 +
       conj(IT_0230)*IT_0297 + conj(IT_0298)*IT_0300 + conj(IT_0254)*IT_0302) +
       IT_0254*(conj(IT_0254)*IT_0287 + conj(IT_0210)*IT_0294 + conj(IT_0164)
      *IT_0296 + IT_0292*conj(IT_0298) + conj(IT_0259)*IT_0300 + conj(IT_0280)
      *IT_0302 + IT_0258*(IT_0305 + IT_0306 + IT_0307)) + IT_0230*(conj(IT_0280)
      *IT_0297 + IT_0229*IT_0308 + IT_0232*IT_0311 + conj(IT_0298)*IT_0312 +
       conj(IT_0210)*IT_0313 + conj(IT_0231)*IT_0314 + conj(IT_0230)*IT_0315 +
       conj(IT_0164)*IT_0316) + IT_0231*(conj(IT_0280)*IT_0281 + IT_0232*IT_0308
       + IT_0229*IT_0311 + conj(IT_0230)*IT_0314 + conj(IT_0231)*IT_0315 + conj
      (IT_0210)*IT_0317 + conj(IT_0298)*IT_0318 + conj(IT_0164)*IT_0319) +
       IT_0332 + IT_0298*(conj(IT_0254)*IT_0292 + conj(IT_0280)*IT_0300 + conj
      (IT_0230)*IT_0312 + conj(IT_0231)*IT_0318 + conj(IT_0164)*IT_0329 + conj
      (IT_0210)*IT_0331 + conj(IT_0298)*IT_0334 + conj(IT_0259)*IT_0336) + 2
      *IT_0259*(IT_0263*(0.5*IT_0309 + 0.5*IT_0310 + IT_0339) + IT_0334*IT_0340 
      + IT_0300*IT_0341 + IT_0336*IT_0342 + IT_0292*IT_0343 + IT_0329*IT_0344 +
       IT_0331*IT_0345) + 2*IT_0164*IT_0346;
    return create_ccomplex_return(IT_0347);
}

