#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_1_to_A_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_1_to_A_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + 0.5*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = pow(m_N_1, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0020 = IT_0003*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0018 + (-0.5)*IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = IT_0013*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = s_13*s_23;
    const ccomplex_t IT_0028 = pow(m_W, -2);
    const ccomplex_t IT_0029 = pow(s_24, 2);
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = -IT_0012;
    const ccomplex_t IT_0032 = IT_0030 + IT_0031;
    const ccomplex_t IT_0033 = m_C_2*IT_0032;
    const ccomplex_t IT_0034 = (-2)*IT_0000;
    const ccomplex_t IT_0035 = IT_0015*IT_0034;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0000*IT_0015;
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = conj(IT_0036) + conj(IT_0038);
    const ccomplex_t IT_0040 = m_C_2*IT_0000;
    const ccomplex_t IT_0041 = IT_0015*IT_0040;
    const ccomplex_t IT_0042 = 2*IT_0041;
    const ccomplex_t IT_0043 = pow(m_W, 2);
    const ccomplex_t IT_0044 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0043 
      + -reg_prop, -1);
    const ccomplex_t IT_0045 = sin(beta);
    const ccomplex_t IT_0046 = IT_0003*IT_0045;
    const ccomplex_t IT_0047 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 2*IT_0048;
    const ccomplex_t IT_0050 = cos(theta_W);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0057 = IT_0046*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = IT_0055 + IT_0058;
    const ccomplex_t IT_0060 = IT_0049 + IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = IT_0044*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = m_N_1*IT_0000;
    const ccomplex_t IT_0068 = IT_0025*IT_0067;
    const ccomplex_t IT_0069 = -IT_0066 + IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0069;
    const ccomplex_t IT_0071 = IT_0042 + IT_0070;
    const ccomplex_t IT_0072 = IT_0036 + IT_0038;
    const ccomplex_t IT_0073 = 2*conj(IT_0041);
    const ccomplex_t IT_0074 = (-2)*conj(IT_0069);
    const ccomplex_t IT_0075 = IT_0073 + IT_0074;
    const ccomplex_t IT_0076 = IT_0025*IT_0034;
    const ccomplex_t IT_0077 = -IT_0076;
    const ccomplex_t IT_0078 = IT_0015*IT_0067;
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = 2*IT_0079;
    const ccomplex_t IT_0081 = cos(beta);
    const ccomplex_t IT_0082 = IT_0003*IT_0081;
    const ccomplex_t IT_0083 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = IT_0051*IT_0081;
    const ccomplex_t IT_0087 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0091 = IT_0082*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = IT_0089 + IT_0092;
    const ccomplex_t IT_0094 = IT_0085 + IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = IT_0063*IT_0096;
    const ccomplex_t IT_0098 = IT_0044*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0025*IT_0040;
    const ccomplex_t IT_0101 = IT_0099 + -IT_0100;
    const ccomplex_t IT_0102 = (-2)*IT_0101;
    const ccomplex_t IT_0103 = IT_0080 + IT_0102;
    const ccomplex_t IT_0104 = IT_0000*IT_0025;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0077 + IT_0105;
    const ccomplex_t IT_0107 = 2*conj(IT_0079);
    const ccomplex_t IT_0108 = (-2)*conj(IT_0101);
    const ccomplex_t IT_0109 = IT_0107 + IT_0108;
    const ccomplex_t IT_0110 = s_14*s_24*IT_0028;
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = s_12 + IT_0111;
    const ccomplex_t IT_0113 = m_N_1*IT_0112;
    const ccomplex_t IT_0114 = m_C_2*IT_0112;
    const ccomplex_t IT_0115 = IT_0000*IT_0010;
    const ccomplex_t IT_0116 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = IT_0039 + conj(IT_0118);
    const ccomplex_t IT_0120 = 2*IT_0069;
    const ccomplex_t IT_0121 = (-2)*IT_0041;
    const ccomplex_t IT_0122 = IT_0120 + IT_0121;
    const ccomplex_t IT_0123 = 2*conj(IT_0069);
    const ccomplex_t IT_0124 = (-2)*conj(IT_0041);
    const ccomplex_t IT_0125 = IT_0123 + IT_0124;
    const ccomplex_t IT_0126 = IT_0000*IT_0023;
    const ccomplex_t IT_0127 = IT_0116*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = -IT_0076 + IT_0128;
    const ccomplex_t IT_0130 = 2*IT_0101;
    const ccomplex_t IT_0131 = (-2)*IT_0079;
    const ccomplex_t IT_0132 = IT_0130 + IT_0131;
    const ccomplex_t IT_0133 = 2*conj(IT_0101);
    const ccomplex_t IT_0134 = (-2)*conj(IT_0079);
    const ccomplex_t IT_0135 = IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = pow(s_14, 2);
    const ccomplex_t IT_0137 = IT_0028*IT_0136;
    const ccomplex_t IT_0138 = -IT_0011;
    const ccomplex_t IT_0139 = IT_0137 + IT_0138;
    const ccomplex_t IT_0140 = m_N_1*IT_0139;
    const ccomplex_t IT_0141 = s_13*m_N_1;
    const ccomplex_t IT_0142 = -IT_0101;
    const ccomplex_t IT_0143 = IT_0079 + IT_0142;
    const ccomplex_t IT_0144 = -conj(IT_0101);
    const ccomplex_t IT_0145 = conj(IT_0079) + IT_0144;
    const ccomplex_t IT_0146 = 6*IT_0041;
    const ccomplex_t IT_0147 = (-6)*IT_0069;
    const ccomplex_t IT_0148 = IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = 6*conj(IT_0041);
    const ccomplex_t IT_0150 = (-6)*conj(IT_0069);
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = s_23*m_C_2;
    const ccomplex_t IT_0153 = 6*IT_0069;
    const ccomplex_t IT_0154 = (-6)*IT_0041;
    const ccomplex_t IT_0155 = IT_0153 + IT_0154;
    const ccomplex_t IT_0156 = 6*conj(IT_0069);
    const ccomplex_t IT_0157 = (-6)*conj(IT_0041);
    const ccomplex_t IT_0158 = IT_0156 + IT_0157;
    const ccomplex_t IT_0159 = m_C_2*m_N_1;
    const ccomplex_t IT_0160 = s_13*IT_0043;
    const ccomplex_t IT_0161 = IT_0028*IT_0160;
    const ccomplex_t IT_0162 = (-0.25)*IT_0161;
    const ccomplex_t IT_0163 = s_13 + IT_0162;
    const ccomplex_t IT_0164 = IT_0159*IT_0163;
    const ccomplex_t IT_0165 = (-4)*IT_0164;
    const ccomplex_t IT_0166 = s_23*IT_0043;
    const ccomplex_t IT_0167 = IT_0028*IT_0166;
    const ccomplex_t IT_0168 = (-0.25)*IT_0167;
    const ccomplex_t IT_0169 = s_23 + IT_0168;
    const ccomplex_t IT_0170 = IT_0159*IT_0169;
    const ccomplex_t IT_0171 = (-4)*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0128;
    const ccomplex_t IT_0173 = 2*IT_0036;
    const ccomplex_t IT_0174 = conj(IT_0172)*IT_0173;
    const ccomplex_t IT_0175 = 2*conj(IT_0036);
    const ccomplex_t IT_0176 = IT_0172*IT_0175;
    const ccomplex_t IT_0177 = IT_0012*IT_0112;
    const ccomplex_t IT_0178 = 2*conj(IT_0077);
    const ccomplex_t IT_0179 = conj(IT_0038)*IT_0173;
    const ccomplex_t IT_0180 = IT_0038*IT_0175;
    const ccomplex_t IT_0181 = IT_0011*IT_0112;
    const ccomplex_t IT_0182 = 2*conj(IT_0129);
    const ccomplex_t IT_0183 = 2*IT_0038;
    const ccomplex_t IT_0184 = 2*conj(IT_0038);
    const ccomplex_t IT_0185 = s_12*IT_0139;
    const ccomplex_t IT_0186 = s_12*IT_0032;
    const ccomplex_t IT_0187 = 4*IT_0105*conj(IT_0105);
    const ccomplex_t IT_0188 = 4*IT_0038*conj(IT_0038);
    const ccomplex_t IT_0189 = 6*s_12;
    const ccomplex_t IT_0190 = pow(s_34, 2);
    const ccomplex_t IT_0191 = s_12*IT_0028*IT_0190;
    const ccomplex_t IT_0192 = 0.5*IT_0191;
    const ccomplex_t IT_0193 = s_13*s_24*s_34*IT_0028;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = s_14*s_23*s_34*IT_0028;
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0027 + IT_0192 + IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = -IT_0000;
    const ccomplex_t IT_0199 = IT_0015*IT_0198;
    const ccomplex_t IT_0200 = 8*IT_0199;
    const ccomplex_t IT_0201 = 8*conj(IT_0199);
    const ccomplex_t IT_0202 = IT_0025*IT_0198;
    const ccomplex_t IT_0203 = 8*IT_0202;
    const ccomplex_t IT_0204 = 8*conj(IT_0202);
    const ccomplex_t IT_0205 = m_C_2*m_N_1*IT_0028*IT_0190;
    const ccomplex_t IT_0206 = 4*IT_0199;
    const ccomplex_t IT_0207 = 4*IT_0016;
    const ccomplex_t IT_0208 = 4*conj(IT_0016);
    const ccomplex_t IT_0209 = (-0.5)*IT_0118;
    const ccomplex_t IT_0210 = 4*conj(IT_0199)*(IT_0026 + 3*IT_0202) + conj
      (IT_0026)*IT_0206 + conj(IT_0202)*(12*IT_0199 + IT_0207) + conj(IT_0172)*
      (IT_0200 + IT_0207) + IT_0202*IT_0208 + IT_0172*(IT_0201 + IT_0208) + 4*
      (conj(IT_0026) + 0.25*IT_0204)*IT_0209 + 4*(IT_0026 + 0.25*IT_0203)*conj
      (IT_0209);
    const ccomplex_t IT_0211 = IT_0032*IT_0159;
    const ccomplex_t IT_0212 = conj(IT_0105)*IT_0173;
    const ccomplex_t IT_0213 = IT_0105*IT_0175;
    const ccomplex_t IT_0214 = 8*conj(IT_0038)*IT_0105;
    const ccomplex_t IT_0215 = IT_0139*IT_0159;
    const ccomplex_t IT_0216 = s_12*IT_0112;
    const ccomplex_t IT_0217 = (-8)*IT_0038;
    const ccomplex_t IT_0218 = (-6)*s_12;
    const ccomplex_t IT_0219 = (-6)*s_13;
    const ccomplex_t IT_0220 = (-6)*s_23;
    const ccomplex_t IT_0221 = (-6)*IT_0011;
    const ccomplex_t IT_0222 = (-6)*IT_0012;
    const ccomplex_t IT_0223 = IT_0218 + IT_0222;
    const ccomplex_t IT_0224 = IT_0219 + IT_0220;
    const ccomplex_t IT_0225 = conj(IT_0101)*IT_0148 + IT_0101*IT_0151 + conj
      (IT_0079)*IT_0155 + IT_0079*IT_0158 + (conj(IT_0129)*IT_0199 + IT_0129
      *conj(IT_0199) + conj(IT_0118)*IT_0202 + IT_0118*conj(IT_0202))*IT_0219 +
       conj(IT_0077)*(IT_0118*IT_0218 + IT_0199*IT_0220) + IT_0077*(conj(IT_0118
      )*IT_0218 + conj(IT_0199)*IT_0220) + (conj(IT_0118)*IT_0129 + IT_0118*conj
      (IT_0129))*IT_0221 + conj(IT_0036)*(IT_0129*(IT_0218 + IT_0221) + IT_0077
      *IT_0223 + IT_0202*IT_0224) + IT_0036*(conj(IT_0129)*(IT_0218 + IT_0221) +
       conj(IT_0077)*IT_0223 + conj(IT_0202)*IT_0224);
    const ccomplex_t IT_0226 = s_12*IT_0043;
    const ccomplex_t IT_0227 = s_14*s_24;
    const ccomplex_t IT_0228 = (-2)*IT_0227;
    const ccomplex_t IT_0229 = IT_0226 + IT_0228;
    const ccomplex_t IT_0230 = IT_0028*IT_0229;
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = s_12 + IT_0231;
    const ccomplex_t IT_0233 = (-4)*s_13;
    const ccomplex_t IT_0234 = (-4)*s_23;
    const ccomplex_t IT_0235 = (-4)*s_12;
    const ccomplex_t IT_0236 = (-4)*IT_0011;
    const ccomplex_t IT_0237 = IT_0036*IT_0236;
    const ccomplex_t IT_0238 = (-4)*IT_0012;
    const ccomplex_t IT_0239 = (conj(IT_0118)*IT_0199 + (IT_0036 + IT_0118)
      *conj(IT_0199) + conj(IT_0129)*IT_0202 + IT_0129*conj(IT_0202))*IT_0233 + 
      (IT_0036*conj(IT_0199) + conj(IT_0077)*IT_0202)*IT_0234 + IT_0129*(conj
      (IT_0077)*IT_0235 + conj(IT_0129)*IT_0236) + conj(IT_0118)*(IT_0036
      *IT_0235 + IT_0118*IT_0236 + IT_0237) + IT_0077*(conj(IT_0202)*IT_0234 +
       conj(IT_0129)*IT_0235 + conj(IT_0077)*IT_0238) + conj(IT_0036)*(IT_0199*
      (IT_0233 + IT_0234) + IT_0118*(IT_0235 + IT_0236) + IT_0237 + (-8)*IT_0036
      *(s_12 + (-0.125)*IT_0238));
    const ccomplex_t IT_0240 = s_14*s_34*IT_0028;
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = s_13 + IT_0241;
    const ccomplex_t IT_0243 = s_12*IT_0242;
    const ccomplex_t IT_0244 = s_24*s_34*IT_0028;
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = s_23 + IT_0245;
    const ccomplex_t IT_0247 = s_12*IT_0246;
    const ccomplex_t IT_0248 = (-4)*IT_0199;
    const ccomplex_t IT_0249 = conj(IT_0038)*IT_0248;
    const ccomplex_t IT_0250 = (-4)*conj(IT_0199);
    const ccomplex_t IT_0251 = IT_0038*IT_0250;
    const ccomplex_t IT_0252 = (-4)*IT_0202;
    const ccomplex_t IT_0253 = conj(IT_0105)*IT_0252;
    const ccomplex_t IT_0254 = (-4)*conj(IT_0202);
    const ccomplex_t IT_0255 = IT_0105*IT_0254;
    const ccomplex_t IT_0256 = s_23*IT_0028*IT_0136;
    const ccomplex_t IT_0257 = s_23*IT_0011;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = s_12*s_14*s_34*IT_0028;
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = s_12*s_13;
    const ccomplex_t IT_0262 = 2*IT_0261;
    const ccomplex_t IT_0263 = s_13*s_14*s_24*IT_0028;
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = IT_0256 + IT_0258 + IT_0260 + IT_0262 + IT_0264;
    const ccomplex_t IT_0266 = s_14*s_23*s_24*IT_0028;
    const ccomplex_t IT_0267 = s_13*IT_0012;
    const ccomplex_t IT_0268 = s_12*s_24*s_34*IT_0028;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = s_13*IT_0028*IT_0029;
    const ccomplex_t IT_0271 = -IT_0270;
    const ccomplex_t IT_0272 = IT_0266 + IT_0267 + IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = (-4)*IT_0105*conj(IT_0172);
    const ccomplex_t IT_0274 = (-4)*conj(IT_0105)*IT_0172;
    const ccomplex_t IT_0275 = (-4)*IT_0038*conj(IT_0209);
    const ccomplex_t IT_0276 = (-4)*conj(IT_0038)*IT_0209;
    const ccomplex_t IT_0277 = IT_0012*IT_0139;
    const ccomplex_t IT_0278 = (-2)*conj(IT_0077);
    const ccomplex_t IT_0279 = (-2)*IT_0036;
    const ccomplex_t IT_0280 = conj(IT_0038)*IT_0279;
    const ccomplex_t IT_0281 = (-2)*conj(IT_0036);
    const ccomplex_t IT_0282 = IT_0038*IT_0281;
    const ccomplex_t IT_0283 = IT_0112*IT_0159;
    const ccomplex_t IT_0284 = IT_0159*IT_0246;
    const ccomplex_t IT_0285 = 2*IT_0129;
    const ccomplex_t IT_0286 = 2*IT_0077;
    const ccomplex_t IT_0287 = (-8)*IT_0105;
    const ccomplex_t IT_0288 = (-4)*IT_0036;
    const ccomplex_t IT_0289 = IT_0011*IT_0032;
    const ccomplex_t IT_0290 = (-2)*IT_0129;
    const ccomplex_t IT_0291 = (-2)*IT_0077;
    const ccomplex_t IT_0292 = IT_0256 + IT_0258 + IT_0260 + IT_0263;
    const ccomplex_t IT_0293 = -IT_0267;
    const ccomplex_t IT_0294 = IT_0266 + IT_0269 + IT_0270 + IT_0293;
    const ccomplex_t IT_0295 = (-2)*IT_0202;
    const ccomplex_t IT_0296 = 0.5*IT_0033*IT_0103 + 0.5*IT_0071*IT_0113 + 0.5
      *IT_0114*IT_0132 + 0.5*IT_0122*IT_0140 + IT_0118*IT_0215 + 4*IT_0038*
      (IT_0211 + IT_0215 + (-2)*IT_0283) + (-5)*IT_0199*IT_0284 + 0.5*IT_0181
      *IT_0285 + 0.5*IT_0177*IT_0286 + 0.5*IT_0216*IT_0287 + 0.5*IT_0283*IT_0288
       + 0.5*IT_0289*IT_0290 + 0.5*IT_0277*IT_0291 + 0.5*(IT_0292 + IT_0294)
      *IT_0295;
    const ccomplex_t IT_0297 = (-4)*conj(IT_0036);
    const ccomplex_t IT_0298 = (-2)*IT_0105;
    const ccomplex_t IT_0299 = (-2)*conj(IT_0105);
    const ccomplex_t IT_0300 = (-2)*IT_0038;
    const ccomplex_t IT_0301 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0302 = (-2)*conj(IT_0129);
    const ccomplex_t IT_0303 = m_N_1*IT_0242;
    const ccomplex_t IT_0304 = IT_0011*IT_0246;
    const ccomplex_t IT_0305 = IT_0012*IT_0242;
    const ccomplex_t IT_0306 = m_C_2*IT_0246;
    const ccomplex_t IT_0307 = -IT_0143;
    const ccomplex_t IT_0308 = s_23*IT_0011*IT_0043;
    const ccomplex_t IT_0309 = s_24*s_34*IT_0011;
    const ccomplex_t IT_0310 = (-2)*IT_0309;
    const ccomplex_t IT_0311 = IT_0308 + IT_0310;
    const ccomplex_t IT_0312 = IT_0028*IT_0311;
    const ccomplex_t IT_0313 = (-0.5)*IT_0312;
    const ccomplex_t IT_0314 = IT_0257 + IT_0313;
    const ccomplex_t IT_0315 = s_12*s_23*IT_0043;
    const ccomplex_t IT_0316 = s_13*IT_0029;
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = s_14*s_23*s_24;
    const ccomplex_t IT_0319 = s_12*s_24*s_34;
    const ccomplex_t IT_0320 = s_13*IT_0012*IT_0043;
    const ccomplex_t IT_0321 = -IT_0318 + -IT_0319 + (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = IT_0317 + IT_0321;
    const ccomplex_t IT_0323 = IT_0028*IT_0322;
    const ccomplex_t IT_0324 = s_12*s_23;
    const ccomplex_t IT_0325 = (-2)*IT_0324;
    const ccomplex_t IT_0326 = IT_0267 + IT_0323 + IT_0325;
    const ccomplex_t IT_0327 = s_13*s_23*IT_0043;
    const ccomplex_t IT_0328 = s_12*IT_0190;
    const ccomplex_t IT_0329 = IT_0327 + IT_0328;
    const ccomplex_t IT_0330 = s_14*s_23*s_34;
    const ccomplex_t IT_0331 = s_13*s_24*s_34;
    const ccomplex_t IT_0332 = -IT_0330 + -IT_0331;
    const ccomplex_t IT_0333 = IT_0329 + IT_0332;
    const ccomplex_t IT_0334 = IT_0028*IT_0333;
    const ccomplex_t IT_0335 = (-0.5)*IT_0334;
    const ccomplex_t IT_0336 = IT_0027 + IT_0335;
    const ccomplex_t IT_0337 = (-2)*IT_0331;
    const ccomplex_t IT_0338 = IT_0327 + IT_0337;
    const ccomplex_t IT_0339 = IT_0028*IT_0338;
    const ccomplex_t IT_0340 = (-8)*IT_0339;
    const ccomplex_t IT_0341 = 16*IT_0027;
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = s_14*s_34;
    const ccomplex_t IT_0344 = (-2)*IT_0343;
    const ccomplex_t IT_0345 = IT_0160 + IT_0344;
    const ccomplex_t IT_0346 = IT_0028*IT_0345;
    const ccomplex_t IT_0347 = (-0.5)*IT_0346;
    const ccomplex_t IT_0348 = s_13 + IT_0347;
    const ccomplex_t IT_0349 = m_N_1*IT_0348;
    const ccomplex_t IT_0350 = (-2)*IT_0349;
    const ccomplex_t IT_0351 = 0.5*IT_0350;
    const ccomplex_t IT_0352 = s_24*s_34;
    const ccomplex_t IT_0353 = (-2)*IT_0352;
    const ccomplex_t IT_0354 = IT_0166 + IT_0353;
    const ccomplex_t IT_0355 = IT_0028*IT_0354;
    const ccomplex_t IT_0356 = (-0.5)*IT_0355;
    const ccomplex_t IT_0357 = s_23 + IT_0356;
    const ccomplex_t IT_0358 = m_C_2*IT_0357;
    const ccomplex_t IT_0359 = (-2)*IT_0358;
    const ccomplex_t IT_0360 = 0.5*IT_0359;
    const ccomplex_t IT_0361 = s_34*IT_0028;
    const ccomplex_t IT_0362 = s_14*s_23;
    const ccomplex_t IT_0363 = s_13*s_24;
    const ccomplex_t IT_0364 = s_12*s_34;
    const ccomplex_t IT_0365 = -IT_0363 + -IT_0364;
    const ccomplex_t IT_0366 = IT_0362 + IT_0365;
    const ccomplex_t IT_0367 = IT_0361*IT_0366;
    const ccomplex_t IT_0368 = 0.5*IT_0367;
    const ccomplex_t IT_0369 = 2*conj(IT_0172)*(IT_0118*IT_0165 + (-2)*IT_0129
      *IT_0314 + 2*IT_0077*IT_0326 + (-4)*IT_0026*IT_0336 + 0.5*IT_0172*IT_0342 
      + IT_0071*IT_0351 + IT_0132*IT_0360 + IT_0252*IT_0368);
    const ccomplex_t IT_0370 = (-0.25)*IT_0305;
    const ccomplex_t IT_0371 = (-0.5)*IT_0304;
    const ccomplex_t IT_0372 = (-0.25)*IT_0350;
    const ccomplex_t IT_0373 = (-0.25)*IT_0359;
    const ccomplex_t IT_0374 = (-0.25)*IT_0165;
    const ccomplex_t IT_0375 = (-0.25)*IT_0171;
    const ccomplex_t IT_0376 = (-0.25)*IT_0367;
    const ccomplex_t IT_0377 = -(IT_0118 + (-0.25)*IT_0288)*IT_0314 + IT_0036
      *IT_0326 + (-2)*IT_0016*IT_0336 + 0.25*IT_0209*IT_0342 + -IT_0103*IT_0372 
      + -IT_0122*IT_0373 + -IT_0285*IT_0374 + -IT_0286*IT_0375 + -IT_0248*IT_0376;
    const ccomplex_t IT_0378 = -(conj(IT_0118) + (-0.25)*IT_0297)*IT_0314 +
       conj(IT_0036)*IT_0326 + (-2)*conj(IT_0016)*IT_0336 + -IT_0109*IT_0372 + 
      -IT_0125*IT_0373 + -IT_0182*IT_0374 + -IT_0178*IT_0375 + -IT_0250*IT_0376;
    const ccomplex_t IT_0379 = (-2)*IT_0199;
    const ccomplex_t IT_0380 = conj(IT_0038)*IT_0379;
    const ccomplex_t IT_0381 = (-2)*conj(IT_0199);
    const ccomplex_t IT_0382 = IT_0038*IT_0381;
    const ccomplex_t IT_0383 = conj(IT_0016)*IT_0279;
    const ccomplex_t IT_0384 = IT_0016*IT_0281;
    const ccomplex_t IT_0385 = conj(IT_0016)*IT_0300;
    const ccomplex_t IT_0386 = IT_0016*IT_0301;
    const ccomplex_t IT_0387 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0388 = IT_0105*IT_0387;
    const ccomplex_t IT_0389 = (-2)*IT_0016;
    const ccomplex_t IT_0390 = (-2)*conj(IT_0016);
    const ccomplex_t IT_0391 = IT_0379 + IT_0389;
    const ccomplex_t IT_0392 = IT_0381 + IT_0390;
    const ccomplex_t IT_0393 = (-10)*IT_0105*conj(IT_0199);
    const ccomplex_t IT_0394 = (-10)*IT_0038*conj(IT_0202);
    const ccomplex_t IT_0395 = (-10)*conj(IT_0038)*IT_0202;
    const ccomplex_t IT_0396 = conj(IT_0209)*IT_0287;
    const ccomplex_t IT_0397 = (-8)*conj(IT_0105)*IT_0209;
    const ccomplex_t IT_0398 = conj(IT_0172)*IT_0217;
    const ccomplex_t IT_0399 = (-8)*conj(IT_0038)*IT_0172;
    const ccomplex_t IT_0400 = conj(IT_0026)*IT_0279;
    const ccomplex_t IT_0401 = IT_0026*IT_0281;
    const ccomplex_t IT_0402 = conj(IT_0026)*IT_0300;
    const ccomplex_t IT_0403 = IT_0026*IT_0301;
    const ccomplex_t IT_0404 = conj(IT_0105)*IT_0389;
    const ccomplex_t IT_0405 = IT_0105*IT_0390;
    const ccomplex_t IT_0406 = IT_0202*IT_0281;
    const ccomplex_t IT_0407 = conj(IT_0202)*IT_0279;
    const ccomplex_t IT_0408 = IT_0159*IT_0242;
    const ccomplex_t IT_0409 = (-10)*conj(IT_0105)*IT_0199 + (-2)*conj(IT_0118
      )*(IT_0026 + (-0.5)*IT_0295) + (-2)*IT_0118*(conj(IT_0026) + (-0.5)
      *IT_0387) + conj(IT_0129)*IT_0391 + IT_0129*IT_0392 + IT_0393 + IT_0394 +
       IT_0395 + IT_0396 + IT_0397 + IT_0398 + IT_0399 + IT_0400 + IT_0401 +
       IT_0402 + IT_0403 + IT_0404 + IT_0405 + IT_0406 + IT_0407;
    const ccomplex_t IT_0410 = (12*IT_0016*conj(IT_0016) + 12*IT_0026*conj
      (IT_0026))*IT_0027 + IT_0033*(IT_0039*IT_0071 + IT_0072*IT_0075 + conj
      (IT_0077)*IT_0103 + IT_0106*IT_0109) + (IT_0071*conj(IT_0077) + IT_0039
      *IT_0103 + IT_0075*IT_0106 + IT_0072*IT_0109)*IT_0113 + IT_0114*(IT_0119
      *IT_0122 + (IT_0072 + IT_0118)*IT_0125 + conj(IT_0129)*IT_0132 + (IT_0105 
      + IT_0129)*IT_0135) + (IT_0125*(IT_0105 + IT_0129) + IT_0122*conj(IT_0129)
       + IT_0119*IT_0132 + (IT_0072 + IT_0118)*IT_0135)*IT_0140 + 6*IT_0141*
      (conj(IT_0016)*IT_0143 + IT_0016*IT_0145 + 0.166666666666667*conj(IT_0026)
      *IT_0148 + 0.166666666666667*IT_0026*IT_0151) + (-6)*IT_0152*(conj(IT_0026
      )*IT_0143 + IT_0026*IT_0145 + (-0.166666666666667)*conj(IT_0016)*IT_0155 +
       (-0.166666666666667)*IT_0016*IT_0158) + (IT_0165 + IT_0171)*(IT_0174 +
       IT_0176) + IT_0177*(IT_0105*IT_0178 + IT_0179 + IT_0180) + IT_0181*
      (IT_0179 + IT_0180 + IT_0105*IT_0182 + conj(IT_0118)*IT_0183 + IT_0118
      *IT_0184) + (IT_0185 + IT_0186)*(IT_0187 + IT_0188) + (IT_0041*conj
      (IT_0041) + IT_0069*conj(IT_0069) + IT_0079*conj(IT_0079) + IT_0101*conj
      (IT_0101))*IT_0189 + IT_0197*(conj(IT_0016)*IT_0200 + IT_0016*IT_0201 +
       conj(IT_0026)*IT_0203 + IT_0026*IT_0204) + IT_0205*IT_0210 + IT_0211*
      (conj(IT_0077)*IT_0183 + IT_0077*IT_0184 + IT_0212 + IT_0213 + IT_0214) +
       2*(IT_0105*conj(IT_0118) + 0.5*conj(IT_0129)*IT_0183 + 0.5*IT_0129
      *IT_0184 + 0.5*IT_0212 + 0.5*IT_0213 + 0.5*IT_0214)*IT_0215 + conj(IT_0038
      )*IT_0216*IT_0217 + (conj(IT_0041)*IT_0069 + IT_0041*conj(IT_0069) + conj
      (IT_0079)*IT_0101 + IT_0079*conj(IT_0101))*IT_0218 + IT_0159*IT_0225 +
       IT_0232*IT_0239 + (IT_0243 + IT_0247)*(IT_0249 + IT_0251 + IT_0253 +
       IT_0255) + (IT_0265 + IT_0272)*(IT_0273 + IT_0274 + IT_0275 + IT_0276) +
       IT_0277*(IT_0105*IT_0278 + IT_0280 + IT_0282) + 2*conj(IT_0105)*IT_0296 +
       IT_0283*(IT_0105*((-16)*conj(IT_0038) + IT_0297) + conj(IT_0118)*IT_0298 
      + IT_0118*IT_0299 + (conj(IT_0077) + conj(IT_0129))*IT_0300 + (IT_0077 +
       IT_0129)*IT_0301) + IT_0289*(IT_0280 + IT_0282 + conj(IT_0118)*IT_0300 +
       IT_0118*IT_0301 + IT_0105*IT_0302) + 4*conj(IT_0202)*(IT_0191*IT_0202 +
       0.25*IT_0197*IT_0203 + (IT_0041 + -IT_0069)*IT_0303 + 0.25*IT_0285
      *IT_0304 + 0.25*IT_0286*IT_0305 + IT_0306*IT_0307) + 2*IT_0172*(conj
      (IT_0118)*IT_0165 + (-2)*conj(IT_0129)*IT_0314 + 2*conj(IT_0077)*IT_0326 +
       (-4)*conj(IT_0026)*IT_0336 + IT_0075*IT_0351 + IT_0135*IT_0360 + IT_0254
      *IT_0368) + IT_0369 + 4*IT_0202*((conj(IT_0041) + -conj(IT_0069))*IT_0303 
      + 0.25*IT_0182*IT_0304 + -IT_0145*IT_0306 + -IT_0178*IT_0370) + (-4)*conj
      (IT_0199)*((-0.25)*IT_0197*IT_0200 + (-0.25)*IT_0191*IT_0206 + (IT_0041 + 
      -IT_0069)*IT_0306 + IT_0303*IT_0307 + IT_0173*IT_0370 + (IT_0118 + 0.5
      *IT_0173)*IT_0371) + 4*IT_0199*(IT_0145*IT_0303 + -(conj(IT_0041) + -conj
      (IT_0069))*IT_0306 + -IT_0175*IT_0370 + -(conj(IT_0118) + 0.5*IT_0175)
      *IT_0371) + 4*conj(IT_0209)*IT_0377 + 4*IT_0209*IT_0378 + IT_0294*(conj
      (IT_0026)*(IT_0291 + IT_0298) + IT_0026*(IT_0278 + IT_0299) + IT_0380 +
       IT_0382 + IT_0383 + IT_0384 + IT_0385 + IT_0386 + IT_0388) + IT_0292*
      (conj(IT_0026)*(IT_0290 + IT_0298) + IT_0026*(IT_0299 + IT_0302) + IT_0380
       + IT_0382 + IT_0383 + IT_0384 + IT_0385 + IT_0386 + IT_0388 + conj
      (IT_0118)*IT_0389 + IT_0118*IT_0390) + IT_0284*(conj(IT_0077)*IT_0391 +
       IT_0077*IT_0392 + IT_0393 + IT_0394 + IT_0395 + IT_0396 + IT_0397 +
       IT_0398 + IT_0399 + IT_0400 + IT_0401 + IT_0402 + IT_0403 + IT_0404 +
       IT_0405 + IT_0406 + IT_0407) + IT_0408*IT_0409;
    return create_ccomplex_return(IT_0410);
}

