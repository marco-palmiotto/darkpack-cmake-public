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
    const ccomplex_t IT_0001 = m_C_2*IT_0000;
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
    const ccomplex_t IT_0017 = m_N_1*IT_0000;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, 2);
    const ccomplex_t IT_0021 = cpow((-2)*s_12 + -IT_0011 + -IT_0012 + IT_0020 
      + -reg_prop, -1);
    const ccomplex_t IT_0022 = sin(beta);
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0025 = IT_0023*IT_0024;
    const ccomplex_t IT_0026 = 2*IT_0025;
    const ccomplex_t IT_0027 = cos(theta_W);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = IT_0022*IT_0028;
    const ccomplex_t IT_0030 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0034 = IT_0023*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = IT_0032 + IT_0035;
    const ccomplex_t IT_0037 = IT_0026 + IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0041 = IT_0039*IT_0040;
    const ccomplex_t IT_0042 = IT_0021*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0045 = IT_0003*IT_0044;
    const ccomplex_t IT_0046 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0047 = IT_0003*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + (-0.5)*IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = IT_0013*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = IT_0017*IT_0052;
    const ccomplex_t IT_0054 = -IT_0043 + IT_0053;
    const ccomplex_t IT_0055 = cos(beta);
    const ccomplex_t IT_0056 = IT_0003*IT_0055;
    const ccomplex_t IT_0057 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (-2)*IT_0058;
    const ccomplex_t IT_0060 = IT_0028*IT_0055;
    const ccomplex_t IT_0061 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0065 = IT_0056*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = IT_0063 + IT_0066;
    const ccomplex_t IT_0068 = IT_0059 + IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = IT_0040*IT_0070;
    const ccomplex_t IT_0072 = IT_0021*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0001*IT_0052;
    const ccomplex_t IT_0075 = IT_0073 + -IT_0074;
    const ccomplex_t IT_0076 = pow(m_W, -2);
    const ccomplex_t IT_0077 = s_14*s_34*IT_0076;
    const ccomplex_t IT_0078 = -IT_0077;
    const ccomplex_t IT_0079 = s_13 + IT_0078;
    const ccomplex_t IT_0080 = IT_0012*IT_0079;
    const ccomplex_t IT_0081 = -IT_0000;
    const ccomplex_t IT_0082 = IT_0015*IT_0081;
    const ccomplex_t IT_0083 = (-2)*IT_0000;
    const ccomplex_t IT_0084 = IT_0015*IT_0083;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = 2*IT_0085;
    const ccomplex_t IT_0087 = conj(IT_0082)*IT_0086;
    const ccomplex_t IT_0088 = 2*conj(IT_0085);
    const ccomplex_t IT_0089 = IT_0082*IT_0088;
    const ccomplex_t IT_0090 = IT_0087 + IT_0089;
    const ccomplex_t IT_0091 = s_14*s_24*IT_0076;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = s_12 + IT_0092;
    const ccomplex_t IT_0094 = IT_0012*IT_0093;
    const ccomplex_t IT_0095 = IT_0000*IT_0015;
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0086*conj(IT_0096);
    const ccomplex_t IT_0098 = IT_0088*IT_0096;
    const ccomplex_t IT_0099 = IT_0097 + IT_0098;
    const ccomplex_t IT_0100 = pow(s_14, 2);
    const ccomplex_t IT_0101 = IT_0076*IT_0100;
    const ccomplex_t IT_0102 = -IT_0011;
    const ccomplex_t IT_0103 = IT_0101 + IT_0102;
    const ccomplex_t IT_0104 = IT_0012*IT_0103;
    const ccomplex_t IT_0105 = (-2)*IT_0085;
    const ccomplex_t IT_0106 = conj(IT_0096)*IT_0105;
    const ccomplex_t IT_0107 = (-2)*conj(IT_0085);
    const ccomplex_t IT_0108 = IT_0096*IT_0107;
    const ccomplex_t IT_0109 = IT_0106 + IT_0108;
    const ccomplex_t IT_0110 = s_14*s_23*s_24*IT_0076;
    const ccomplex_t IT_0111 = pow(s_24, 2);
    const ccomplex_t IT_0112 = s_13*IT_0076*IT_0111;
    const ccomplex_t IT_0113 = s_13*IT_0012;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = s_12*s_24*s_34*IT_0076;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0110 + IT_0112 + IT_0114 + IT_0116;
    const ccomplex_t IT_0118 = 2*IT_0000;
    const ccomplex_t IT_0119 = IT_0015*IT_0118;
    const ccomplex_t IT_0120 = IT_0105*conj(IT_0119);
    const ccomplex_t IT_0121 = IT_0107*IT_0119;
    const ccomplex_t IT_0122 = (-2)*IT_0082;
    const ccomplex_t IT_0123 = conj(IT_0096)*IT_0122;
    const ccomplex_t IT_0124 = (-2)*conj(IT_0082);
    const ccomplex_t IT_0125 = IT_0096*IT_0124;
    const ccomplex_t IT_0126 = (-2)*IT_0096;
    const ccomplex_t IT_0127 = conj(IT_0119)*IT_0126;
    const ccomplex_t IT_0128 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0129 = IT_0119*IT_0128;
    const ccomplex_t IT_0130 = IT_0120 + IT_0121 + IT_0123 + IT_0125 + IT_0127
       + IT_0129;
    const ccomplex_t IT_0131 = m_C_2*m_N_1;
    const ccomplex_t IT_0132 = s_13*IT_0020;
    const ccomplex_t IT_0133 = IT_0076*IT_0132;
    const ccomplex_t IT_0134 = (-0.25)*IT_0133;
    const ccomplex_t IT_0135 = s_13 + IT_0134;
    const ccomplex_t IT_0136 = IT_0131*IT_0135;
    const ccomplex_t IT_0137 = (-4)*IT_0136;
    const ccomplex_t IT_0138 = s_23*IT_0020;
    const ccomplex_t IT_0139 = IT_0076*IT_0138;
    const ccomplex_t IT_0140 = (-0.25)*IT_0139;
    const ccomplex_t IT_0141 = s_23 + IT_0140;
    const ccomplex_t IT_0142 = IT_0131*IT_0141;
    const ccomplex_t IT_0143 = (-4)*IT_0142;
    const ccomplex_t IT_0144 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0145 = IT_0000*IT_0050;
    const ccomplex_t IT_0146 = IT_0144*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = (-0.5)*IT_0147;
    const ccomplex_t IT_0149 = IT_0086*conj(IT_0148);
    const ccomplex_t IT_0150 = IT_0088*IT_0148;
    const ccomplex_t IT_0151 = s_24*s_34*IT_0076;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = s_23 + IT_0152;
    const ccomplex_t IT_0154 = IT_0011*IT_0153;
    const ccomplex_t IT_0155 = IT_0000*IT_0010;
    const ccomplex_t IT_0156 = IT_0144*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = 2*IT_0082;
    const ccomplex_t IT_0159 = 2*conj(IT_0082);
    const ccomplex_t IT_0160 = IT_0011*IT_0093;
    const ccomplex_t IT_0161 = 2*IT_0096;
    const ccomplex_t IT_0162 = 2*conj(IT_0096);
    const ccomplex_t IT_0163 = s_12*IT_0103;
    const ccomplex_t IT_0164 = IT_0076*IT_0111;
    const ccomplex_t IT_0165 = -IT_0012;
    const ccomplex_t IT_0166 = IT_0164 + IT_0165;
    const ccomplex_t IT_0167 = s_12*IT_0166;
    const ccomplex_t IT_0168 = IT_0163 + IT_0167;
    const ccomplex_t IT_0169 = 4*IT_0096*conj(IT_0096);
    const ccomplex_t IT_0170 = IT_0131*IT_0166;
    const ccomplex_t IT_0171 = IT_0052*IT_0083;
    const ccomplex_t IT_0172 = -IT_0171;
    const ccomplex_t IT_0173 = IT_0000*IT_0052;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = 8*conj(IT_0096);
    const ccomplex_t IT_0176 = IT_0174*IT_0175;
    const ccomplex_t IT_0177 = IT_0103*IT_0131;
    const ccomplex_t IT_0178 = IT_0147 + -IT_0171;
    const ccomplex_t IT_0179 = s_12*IT_0093;
    const ccomplex_t IT_0180 = (-8)*IT_0096;
    const ccomplex_t IT_0181 = m_C_2*IT_0093;
    const ccomplex_t IT_0182 = conj(IT_0085) + conj(IT_0096);
    const ccomplex_t IT_0183 = 2*IT_0054;
    const ccomplex_t IT_0184 = (-2)*IT_0016;
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = IT_0085 + IT_0096;
    const ccomplex_t IT_0187 = IT_0157 + IT_0186;
    const ccomplex_t IT_0188 = 2*conj(IT_0054);
    const ccomplex_t IT_0189 = (-2)*conj(IT_0016);
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0075;
    const ccomplex_t IT_0192 = (-2)*IT_0019;
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = 2*conj(IT_0075);
    const ccomplex_t IT_0195 = (-2)*conj(IT_0019);
    const ccomplex_t IT_0196 = IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = m_C_2*IT_0157;
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = 8*IT_0096;
    const ccomplex_t IT_0200 = (-8)*conj(IT_0096);
    const ccomplex_t IT_0201 = m_C_2*IT_0166;
    const ccomplex_t IT_0202 = 2*IT_0016;
    const ccomplex_t IT_0203 = (-2)*IT_0054;
    const ccomplex_t IT_0204 = IT_0202 + IT_0203;
    const ccomplex_t IT_0205 = 2*conj(IT_0016);
    const ccomplex_t IT_0206 = (-2)*conj(IT_0054);
    const ccomplex_t IT_0207 = IT_0205 + IT_0206;
    const ccomplex_t IT_0208 = 2*IT_0019;
    const ccomplex_t IT_0209 = (-2)*IT_0075;
    const ccomplex_t IT_0210 = IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = 2*conj(IT_0019);
    const ccomplex_t IT_0212 = (-2)*conj(IT_0075);
    const ccomplex_t IT_0213 = IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = (-6)*s_12;
    const ccomplex_t IT_0215 = s_13*m_N_1;
    const ccomplex_t IT_0216 = s_23*IT_0011;
    const ccomplex_t IT_0217 = s_23*IT_0011*IT_0020;
    const ccomplex_t IT_0218 = s_24*s_34*IT_0011;
    const ccomplex_t IT_0219 = (-2)*IT_0218;
    const ccomplex_t IT_0220 = IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = IT_0076*IT_0220;
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = IT_0216 + IT_0222;
    const ccomplex_t IT_0224 = s_12*s_23*IT_0020;
    const ccomplex_t IT_0225 = s_13*IT_0111;
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = s_14*s_23*s_24;
    const ccomplex_t IT_0228 = s_12*s_24*s_34;
    const ccomplex_t IT_0229 = s_13*IT_0012*IT_0020;
    const ccomplex_t IT_0230 = -IT_0227 + -IT_0228 + (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0226 + IT_0230;
    const ccomplex_t IT_0232 = IT_0076*IT_0231;
    const ccomplex_t IT_0233 = s_12*s_23;
    const ccomplex_t IT_0234 = (-2)*IT_0233;
    const ccomplex_t IT_0235 = IT_0113 + IT_0232 + IT_0234;
    const ccomplex_t IT_0236 = s_14*s_34;
    const ccomplex_t IT_0237 = (-2)*IT_0236;
    const ccomplex_t IT_0238 = IT_0132 + IT_0237;
    const ccomplex_t IT_0239 = IT_0076*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = s_13 + IT_0240;
    const ccomplex_t IT_0242 = m_N_1*IT_0241;
    const ccomplex_t IT_0243 = (-2)*IT_0242;
    const ccomplex_t IT_0244 = s_24*s_34;
    const ccomplex_t IT_0245 = (-2)*IT_0244;
    const ccomplex_t IT_0246 = IT_0138 + IT_0245;
    const ccomplex_t IT_0247 = IT_0076*IT_0246;
    const ccomplex_t IT_0248 = (-0.5)*IT_0247;
    const ccomplex_t IT_0249 = s_23 + IT_0248;
    const ccomplex_t IT_0250 = m_C_2*IT_0249;
    const ccomplex_t IT_0251 = (-2)*IT_0250;
    const ccomplex_t IT_0252 = pow(s_34, 2);
    const ccomplex_t IT_0253 = m_C_2*m_N_1*IT_0076*IT_0252;
    const ccomplex_t IT_0254 = 4*IT_0119;
    const ccomplex_t IT_0255 = 8*IT_0082;
    const ccomplex_t IT_0256 = s_34*IT_0076;
    const ccomplex_t IT_0257 = s_14*s_23;
    const ccomplex_t IT_0258 = s_13*s_24;
    const ccomplex_t IT_0259 = s_12*s_34;
    const ccomplex_t IT_0260 = -IT_0258 + -IT_0259;
    const ccomplex_t IT_0261 = IT_0257 + IT_0260;
    const ccomplex_t IT_0262 = IT_0256*IT_0261;
    const ccomplex_t IT_0263 = IT_0052*IT_0081;
    const ccomplex_t IT_0264 = (-4)*IT_0263;
    const ccomplex_t IT_0265 = 2*IT_0137*IT_0157 + (-12)*IT_0198*IT_0215 + (-4
      )*IT_0178*IT_0223 + 4*IT_0172*IT_0235 + IT_0204*IT_0243 + IT_0193*IT_0251 
      + IT_0253*(IT_0254 + IT_0255) + IT_0262*IT_0264;
    const ccomplex_t IT_0266 = IT_0052*IT_0118;
    const ccomplex_t IT_0267 = s_13*s_23;
    const ccomplex_t IT_0268 = s_13*s_23*IT_0020;
    const ccomplex_t IT_0269 = s_12*IT_0252;
    const ccomplex_t IT_0270 = IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = s_14*s_23*s_34;
    const ccomplex_t IT_0272 = s_13*s_24*s_34;
    const ccomplex_t IT_0273 = -IT_0271 + -IT_0272;
    const ccomplex_t IT_0274 = IT_0270 + IT_0273;
    const ccomplex_t IT_0275 = IT_0076*IT_0274;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = IT_0267 + IT_0276;
    const ccomplex_t IT_0278 = (-2)*IT_0272;
    const ccomplex_t IT_0279 = IT_0268 + IT_0278;
    const ccomplex_t IT_0280 = IT_0076*IT_0279;
    const ccomplex_t IT_0281 = (-8)*IT_0280;
    const ccomplex_t IT_0282 = 16*IT_0267;
    const ccomplex_t IT_0283 = IT_0281 + IT_0282;
    const ccomplex_t IT_0284 = 4*conj(IT_0119);
    const ccomplex_t IT_0285 = 8*conj(IT_0082);
    const ccomplex_t IT_0286 = (-4)*conj(IT_0263);
    const ccomplex_t IT_0287 = 2*IT_0137*conj(IT_0157) + (-4)*conj(IT_0178)
      *IT_0223 + 4*conj(IT_0172)*IT_0235 + IT_0207*IT_0243 + IT_0196*IT_0251 + (
      -8)*conj(IT_0266)*IT_0277 + conj(IT_0148)*IT_0283 + IT_0253*(IT_0284 +
       IT_0285) + IT_0262*IT_0286;
    const ccomplex_t IT_0288 = (-0.5)*IT_0157;
    const ccomplex_t IT_0289 = 2*IT_0172;
    const ccomplex_t IT_0290 = 2*IT_0178;
    const ccomplex_t IT_0291 = 4*IT_0085;
    const ccomplex_t IT_0292 = 8*IT_0263;
    const ccomplex_t IT_0293 = (-4)*IT_0085;
    const ccomplex_t IT_0294 = (-4)*IT_0082;
    const ccomplex_t IT_0295 = IT_0210*IT_0243 + 16*IT_0198*IT_0250 + IT_0185
      *IT_0251 + (-8)*IT_0119*IT_0277 + IT_0143*IT_0289 + IT_0137*IT_0290 +
       IT_0235*IT_0291 + IT_0253*IT_0292 + (-4)*IT_0223*(IT_0157 + (-0.25)
      *IT_0293) + IT_0262*IT_0294;
    const ccomplex_t IT_0296 = s_12*IT_0079;
    const ccomplex_t IT_0297 = s_12*IT_0153;
    const ccomplex_t IT_0298 = conj(IT_0174)*IT_0264;
    const ccomplex_t IT_0299 = IT_0174*IT_0286;
    const ccomplex_t IT_0300 = conj(IT_0096)*IT_0294;
    const ccomplex_t IT_0301 = (-4)*conj(IT_0082);
    const ccomplex_t IT_0302 = IT_0096*IT_0301;
    const ccomplex_t IT_0303 = s_23*IT_0076*IT_0100;
    const ccomplex_t IT_0304 = -IT_0216;
    const ccomplex_t IT_0305 = s_12*s_14*s_34*IT_0076;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = s_12*s_13;
    const ccomplex_t IT_0308 = 2*IT_0307;
    const ccomplex_t IT_0309 = s_13*s_14*s_24*IT_0076;
    const ccomplex_t IT_0310 = -IT_0309;
    const ccomplex_t IT_0311 = IT_0303 + IT_0304 + IT_0306 + IT_0308 + IT_0310;
    const ccomplex_t IT_0312 = -IT_0112;
    const ccomplex_t IT_0313 = IT_0110 + IT_0113 + IT_0116 + IT_0312;
    const ccomplex_t IT_0314 = (-4)*IT_0174;
    const ccomplex_t IT_0315 = conj(IT_0148)*IT_0314;
    const ccomplex_t IT_0316 = (-4)*conj(IT_0174);
    const ccomplex_t IT_0317 = IT_0148*IT_0316;
    const ccomplex_t IT_0318 = IT_0096*conj(IT_0288);
    const ccomplex_t IT_0319 = (-4)*IT_0318;
    const ccomplex_t IT_0320 = conj(IT_0096)*IT_0288;
    const ccomplex_t IT_0321 = (-4)*IT_0320;
    const ccomplex_t IT_0322 = m_N_1*IT_0103;
    const ccomplex_t IT_0323 = (conj(IT_0178)*IT_0185 + IT_0178*IT_0190 + 
      (conj(IT_0157) + IT_0182)*IT_0193 + IT_0187*IT_0196 + 4*IT_0174*(0.25
      *IT_0188 + 0.25*IT_0189 + conj(IT_0198)) + conj(IT_0198)*IT_0314 + IT_0198
      *IT_0316)*IT_0322;
    const ccomplex_t IT_0324 = s_12*IT_0020;
    const ccomplex_t IT_0325 = s_14*s_24;
    const ccomplex_t IT_0326 = (-2)*IT_0325;
    const ccomplex_t IT_0327 = IT_0324 + IT_0326;
    const ccomplex_t IT_0328 = IT_0076*IT_0327;
    const ccomplex_t IT_0329 = (-0.5)*IT_0328;
    const ccomplex_t IT_0330 = s_12 + IT_0329;
    const ccomplex_t IT_0331 = m_C_2*IT_0330;
    const ccomplex_t IT_0332 = (-2)*IT_0331;
    const ccomplex_t IT_0333 = IT_0011*IT_0020;
    const ccomplex_t IT_0334 = IT_0076*IT_0333;
    const ccomplex_t IT_0335 = (-0.25)*IT_0334;
    const ccomplex_t IT_0336 = IT_0011 + IT_0335;
    const ccomplex_t IT_0337 = m_N_1*IT_0336;
    const ccomplex_t IT_0338 = (-4)*IT_0337;
    const ccomplex_t IT_0339 = IT_0076*IT_0324;
    const ccomplex_t IT_0340 = (-0.25)*IT_0339;
    const ccomplex_t IT_0341 = s_12 + IT_0340;
    const ccomplex_t IT_0342 = m_N_1*IT_0341;
    const ccomplex_t IT_0343 = (-4)*IT_0342;
    const ccomplex_t IT_0344 = (-2)*IT_0198;
    const ccomplex_t IT_0345 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0346 = m_N_1*IT_0093;
    const ccomplex_t IT_0347 = 4*IT_0174;
    const ccomplex_t IT_0348 = (-4)*conj(IT_0198);
    const ccomplex_t IT_0349 = m_C_2*IT_0153;
    const ccomplex_t IT_0350 = -IT_0054;
    const ccomplex_t IT_0351 = IT_0016 + IT_0350;
    const ccomplex_t IT_0352 = IT_0012*IT_0020;
    const ccomplex_t IT_0353 = (-0.5)*IT_0352;
    const ccomplex_t IT_0354 = IT_0111 + IT_0353;
    const ccomplex_t IT_0355 = IT_0076*IT_0354;
    const ccomplex_t IT_0356 = IT_0012 + IT_0355;
    const ccomplex_t IT_0357 = m_C_2*IT_0356;
    const ccomplex_t IT_0358 = -conj(IT_0054);
    const ccomplex_t IT_0359 = conj(IT_0016) + IT_0358;
    const ccomplex_t IT_0360 = 4*conj(IT_0085);
    const ccomplex_t IT_0361 = IT_0093*IT_0131;
    const ccomplex_t IT_0362 = (-2)*IT_0174;
    const ccomplex_t IT_0363 = (-2)*conj(IT_0174);
    const ccomplex_t IT_0364 = IT_0011*IT_0166;
    const ccomplex_t IT_0365 = IT_0303 + IT_0304 + IT_0306 + IT_0309;
    const ccomplex_t IT_0366 = (-2)*IT_0119;
    const ccomplex_t IT_0367 = (-2)*conj(IT_0119);
    const ccomplex_t IT_0368 = IT_0131*IT_0153;
    const ccomplex_t IT_0369 = IT_0122 + IT_0366;
    const ccomplex_t IT_0370 = IT_0124 + IT_0367;
    const ccomplex_t IT_0371 = (-10)*conj(IT_0082)*IT_0174;
    const ccomplex_t IT_0372 = (-10)*IT_0096*conj(IT_0263);
    const ccomplex_t IT_0373 = (-10)*conj(IT_0096)*IT_0263;
    const ccomplex_t IT_0374 = (-8)*IT_0174;
    const ccomplex_t IT_0375 = conj(IT_0288)*IT_0374;
    const ccomplex_t IT_0376 = (-8)*conj(IT_0174)*IT_0288;
    const ccomplex_t IT_0377 = IT_0148*IT_0200;
    const ccomplex_t IT_0378 = IT_0105*conj(IT_0263);
    const ccomplex_t IT_0379 = IT_0105*conj(IT_0266);
    const ccomplex_t IT_0380 = IT_0107*IT_0263;
    const ccomplex_t IT_0381 = IT_0107*IT_0266;
    const ccomplex_t IT_0382 = IT_0126*conj(IT_0266);
    const ccomplex_t IT_0383 = IT_0128*IT_0266;
    const ccomplex_t IT_0384 = conj(IT_0174)*IT_0366;
    const ccomplex_t IT_0385 = IT_0174*IT_0367;
    const ccomplex_t IT_0386 = IT_0368*(conj(IT_0148)*IT_0180 + conj(IT_0172)
      *IT_0369 + IT_0172*IT_0370 + IT_0371 + IT_0372 + IT_0373 + IT_0375 +
       IT_0376 + IT_0377 + IT_0378 + IT_0379 + IT_0380 + IT_0381 + IT_0382 +
       IT_0383 + IT_0384 + IT_0385);
    const ccomplex_t IT_0387 = IT_0079*IT_0131;
    const ccomplex_t IT_0388 = (-2)*IT_0263;
    const ccomplex_t IT_0389 = (-2)*conj(IT_0263);
    const ccomplex_t IT_0390 = conj(IT_0148)*IT_0180 + conj(IT_0178)*IT_0369 +
       IT_0178*IT_0370 + IT_0371 + IT_0372 + IT_0373 + IT_0375 + IT_0376 +
       IT_0377 + IT_0378 + IT_0379 + IT_0380 + IT_0381 + IT_0382 + IT_0383 +
       IT_0384 + IT_0385 + (-2)*conj(IT_0157)*(IT_0266 + (-0.5)*IT_0388) + (-2)
      *IT_0157*(conj(IT_0266) + (-0.5)*IT_0389);
    const ccomplex_t IT_0391 = m_C_2*IT_0141;
    const ccomplex_t IT_0392 = (-4)*IT_0391;
    const ccomplex_t IT_0393 = s_12*IT_0076*IT_0252;
    const ccomplex_t IT_0394 = 0.5*IT_0393;
    const ccomplex_t IT_0395 = s_13*s_24*s_34*IT_0076;
    const ccomplex_t IT_0396 = (-0.5)*IT_0395;
    const ccomplex_t IT_0397 = s_14*s_23*s_34*IT_0076;
    const ccomplex_t IT_0398 = (-0.5)*IT_0397;
    const ccomplex_t IT_0399 = IT_0267 + IT_0394 + IT_0396 + IT_0398;
    const ccomplex_t IT_0400 = s_23*m_C_2;
    const ccomplex_t IT_0401 = (-6)*IT_0016;
    const ccomplex_t IT_0402 = (-6)*IT_0075;
    const ccomplex_t IT_0403 = 2*conj(IT_0119);
    const ccomplex_t IT_0404 = (-4)*s_13;
    const ccomplex_t IT_0405 = (-4)*s_12;
    const ccomplex_t IT_0406 = (-4)*s_23;
    const ccomplex_t IT_0407 = (-4)*IT_0011;
    const ccomplex_t IT_0408 = IT_0405 + IT_0407;
    const ccomplex_t IT_0409 = (-4)*IT_0012;
    const ccomplex_t IT_0410 = -conj(IT_0198);
    const ccomplex_t IT_0411 = (conj(IT_0082)*IT_0157 + IT_0082*(conj(IT_0085)
       + conj(IT_0157)) + conj(IT_0178)*IT_0263 + IT_0178*conj(IT_0263))*IT_0404
       + IT_0172*(conj(IT_0178)*IT_0405 + conj(IT_0263)*IT_0406) + (IT_0157*conj
      (IT_0157) + IT_0178*conj(IT_0178))*IT_0407 + conj(IT_0085)*(IT_0082
      *IT_0406 + IT_0157*IT_0408) + conj(IT_0172)*(IT_0178*IT_0405 + IT_0263
      *IT_0406 + IT_0172*IT_0409) + IT_0085*(conj(IT_0082)*(IT_0404 + IT_0406) +
       conj(IT_0157)*IT_0408 + (-8)*conj(IT_0085)*(s_12 + (-0.125)*IT_0407 + (
      -0.125)*IT_0409));
    const ccomplex_t IT_0412 = 2*conj(IT_0172);
    const ccomplex_t IT_0413 = 2*conj(IT_0178);
    const ccomplex_t IT_0414 = 8*conj(IT_0263);
    const ccomplex_t IT_0415 = (-4)*conj(IT_0085);
    const ccomplex_t IT_0416 = (-0.125)*IT_0213*IT_0243 + (-0.125)*IT_0190
      *IT_0251 + conj(IT_0119)*IT_0277 + (-0.125)*IT_0283*conj(IT_0288) + (
      -0.125)*IT_0262*IT_0301 + (-0.125)*IT_0235*IT_0360 + (-2)*IT_0250*(conj
      (IT_0198) + IT_0410) + (-0.125)*IT_0143*IT_0412 + (-0.125)*IT_0137*IT_0413
       + (-0.5)*IT_0253*(conj(IT_0266) + 0.25*IT_0414) + 0.5*IT_0223*(conj
      (IT_0157) + (-0.25)*IT_0415);
    const ccomplex_t IT_0417 = 4*IT_0082;
    const ccomplex_t IT_0418 = m_N_1*IT_0079;
    const ccomplex_t IT_0419 = -IT_0075;
    const ccomplex_t IT_0420 = m_N_1*IT_0135;
    const ccomplex_t IT_0421 = (-4)*IT_0420;
    const ccomplex_t IT_0422 = -conj(IT_0075);
    const ccomplex_t IT_0423 = (-6)*conj(IT_0075);
    const ccomplex_t IT_0424 = (-6)*conj(IT_0016);
    const ccomplex_t IT_0425 = 0.166666666666667*IT_0424;
    const ccomplex_t IT_0426 = (-2)*conj(IT_0172);
    const ccomplex_t IT_0427 = (-2)*conj(IT_0178);
    const ccomplex_t IT_0428 = (-6)*conj(IT_0019);
    const ccomplex_t IT_0429 = 6*IT_0400;
    const ccomplex_t IT_0430 = 6*IT_0215;
    const ccomplex_t IT_0431 = (-6)*conj(IT_0054);
    const ccomplex_t IT_0432 = 0.166666666666667*IT_0431;
    const ccomplex_t IT_0433 = 12*conj(IT_0266)*IT_0267 + (-8)*conj(IT_0148)
      *IT_0277 + 4*IT_0253*(conj(IT_0082) + conj(IT_0288)) + IT_0399*IT_0414 + 2
      *(conj(IT_0198) + 0.5*IT_0345)*IT_0421 + IT_0117*(IT_0363 + IT_0426) +
       IT_0365*(IT_0363 + IT_0427) + (conj(IT_0075) + 0.166666666666667*IT_0428)
      *IT_0429 + IT_0430*(conj(IT_0016) + IT_0432);
    const ccomplex_t IT_0434 = (-6)*s_13;
    const ccomplex_t IT_0435 = (-6)*s_23;
    const ccomplex_t IT_0436 = (-6)*IT_0011;
    const ccomplex_t IT_0437 = (-6)*IT_0012;
    const ccomplex_t IT_0438 = IT_0214 + IT_0437;
    const ccomplex_t IT_0439 = IT_0434 + IT_0435;
    const ccomplex_t IT_0440 = (-6)*IT_0019;
    const ccomplex_t IT_0441 = (-6)*IT_0054;
    const ccomplex_t IT_0442 = 0.166666666666667*IT_0441;
    const ccomplex_t IT_0443 = 6*conj(IT_0019)*(IT_0054 + IT_0198 +
       0.166666666666667*IT_0401) + conj(IT_0198)*IT_0402 + IT_0198*IT_0423 + 6
      *IT_0019*(conj(IT_0054) + conj(IT_0198) + IT_0425) + IT_0198*IT_0428 + 6
      *IT_0075*(conj(IT_0016) + conj(IT_0198) + IT_0432) + (conj(IT_0082)
      *IT_0178 + IT_0082*conj(IT_0178) + conj(IT_0157)*IT_0263 + IT_0157*conj
      (IT_0263))*IT_0434 + conj(IT_0172)*(IT_0157*IT_0214 + IT_0082*IT_0435) +
       IT_0172*(conj(IT_0157)*IT_0214 + conj(IT_0082)*IT_0435) + (conj(IT_0157)
      *IT_0178 + IT_0157*conj(IT_0178))*IT_0436 + conj(IT_0085)*(IT_0178*
      (IT_0214 + IT_0436) + IT_0172*IT_0438 + IT_0263*IT_0439) + IT_0085*(conj
      (IT_0178)*(IT_0214 + IT_0436) + conj(IT_0172)*IT_0438 + conj(IT_0263)
      *IT_0439) + conj(IT_0198)*IT_0440 + 6*conj(IT_0075)*(IT_0016 + IT_0198 +
       IT_0442);
    const ccomplex_t IT_0444 = 4*IT_0263;
    const ccomplex_t IT_0445 = 0.5*IT_0086;
    const ccomplex_t IT_0446 = 0.5*IT_0158*IT_0251 + 0.5*IT_0290*IT_0338 + 0.5
      *IT_0289*IT_0343 + 0.5*IT_0255*IT_0349 + (-4)*IT_0341*IT_0351 + 0.5
      *IT_0293*IT_0357 + 0.5*IT_0366*IT_0392 + IT_0263*IT_0421 + 0.5*IT_0418
      *IT_0444 + IT_0332*(IT_0157 + IT_0445);
    const ccomplex_t IT_0447 = (-2)*IT_0172;
    const ccomplex_t IT_0448 = IT_0117 + IT_0365;
    const ccomplex_t IT_0449 = (-2)*IT_0178;
    const ccomplex_t IT_0450 = (-4)*IT_0198;
    const ccomplex_t IT_0451 = 0.25*IT_0450;
    const ccomplex_t IT_0452 = 2*IT_0177;
    const ccomplex_t IT_0453 = IT_0181*IT_0193 + IT_0170*(IT_0086 + IT_0199) +
       IT_0201*IT_0210 + IT_0094*IT_0289 + IT_0160*IT_0290 + 4*(0.25*IT_0183 +
       0.25*IT_0184 + IT_0198)*IT_0322 + IT_0168*IT_0347 + ((-16)*IT_0096 +
       IT_0293)*IT_0361 + IT_0179*IT_0374 + IT_0082*((-10)*IT_0368 + (-10)
      *IT_0387) + IT_0104*IT_0447 + IT_0388*IT_0448 + IT_0364*IT_0449 + 4
      *IT_0346*(IT_0198 + 0.25*IT_0202 + 0.25*IT_0203 + IT_0451) + (IT_0157 +
       0.5*IT_0199 + IT_0445)*IT_0452;
    const ccomplex_t IT_0454 = 0.5*IT_0088;
    const ccomplex_t IT_0455 = 12*conj(IT_0148)*IT_0215 + IT_0159*IT_0251 + (
      -16)*IT_0250*conj(IT_0288) + IT_0285*IT_0349 + (-8)*IT_0341*IT_0359 +
       IT_0367*IT_0392 + IT_0343*IT_0412 + IT_0338*IT_0413 + IT_0357*IT_0415 + 2
      *conj(IT_0266)*IT_0421 + 4*conj(IT_0263)*(IT_0418 + 0.5*IT_0421) + 2
      *IT_0332*(conj(IT_0157) + IT_0454);
    const ccomplex_t IT_0456 = IT_0344*IT_0421;
    const ccomplex_t IT_0457 = IT_0292*IT_0399;
    const ccomplex_t IT_0458 = conj(IT_0263)*(IT_0253*(12*IT_0082 + IT_0254) +
       IT_0080*IT_0289 + IT_0154*IT_0290 + (-4)*IT_0349*(IT_0019 + IT_0419) +
       IT_0393*IT_0444 + 4*IT_0418*(IT_0351 + IT_0451) + IT_0456 + IT_0457);
    const ccomplex_t IT_0459 = 6*s_12*(IT_0016*conj(IT_0016) + IT_0019*conj
      (IT_0019) + IT_0054*conj(IT_0054) + IT_0075*conj(IT_0075)) + IT_0080
      *IT_0090 + IT_0094*IT_0099 + IT_0104*IT_0109 + IT_0117*IT_0130 + (IT_0137 
      + IT_0143)*(IT_0149 + IT_0150) + IT_0154*(IT_0090 + conj(IT_0157)*IT_0158 
      + IT_0157*IT_0159) + IT_0160*(IT_0099 + conj(IT_0157)*IT_0161 + IT_0157
      *IT_0162) + IT_0168*IT_0169 + IT_0170*(IT_0162*IT_0172 + IT_0161*conj
      (IT_0172) + IT_0176) + IT_0177*(IT_0176 + IT_0162*IT_0178 + IT_0161*conj
      (IT_0178)) + conj(IT_0096)*IT_0179*IT_0180 + IT_0181*((conj(IT_0157) +
       IT_0182)*IT_0185 + IT_0187*IT_0190 + conj(IT_0178)*IT_0193 + (IT_0174 +
       IT_0178)*IT_0196 + IT_0175*IT_0198 + IT_0180*conj(IT_0198) + conj(IT_0198
      )*IT_0199 + IT_0198*IT_0200) + IT_0201*(IT_0175*IT_0198 + IT_0180*conj
      (IT_0198) + conj(IT_0198)*IT_0199 + IT_0198*IT_0200 + IT_0182*IT_0204 +
       IT_0186*IT_0207 + conj(IT_0172)*IT_0210 + (IT_0172 + IT_0174)*IT_0213) + 
      (conj(IT_0016)*IT_0054 + IT_0016*conj(IT_0054) + conj(IT_0019)*IT_0075 +
       IT_0019*conj(IT_0075))*IT_0214 + conj(IT_0148)*IT_0265 + IT_0148*IT_0287 
      + conj(IT_0288)*IT_0295 + (IT_0296 + IT_0297)*(IT_0298 + IT_0299 + IT_0300
       + IT_0302) + (IT_0311 + IT_0313)*(IT_0315 + IT_0317 + IT_0319 + IT_0321) 
      + IT_0323 + (conj(IT_0157)*IT_0332 + conj(IT_0178)*IT_0338 + conj(IT_0172)
      *IT_0343)*IT_0344 + (IT_0157*IT_0332 + IT_0178*IT_0338 + IT_0172*IT_0343)
      *IT_0345 + IT_0346*(conj(IT_0172)*IT_0204 + IT_0172*IT_0207 + IT_0182
      *IT_0210 + IT_0186*IT_0213 + conj(IT_0198)*IT_0347 + IT_0174*(IT_0205 +
       IT_0206 + IT_0348)) + (-8)*conj(IT_0198)*((-0.125)*IT_0122*IT_0251 + (
      -0.125)*IT_0105*IT_0332 + IT_0082*IT_0349 + -IT_0341*IT_0351 + (-0.125)
      *IT_0291*IT_0357) + 8*IT_0198*(0.125*IT_0124*IT_0251 + 0.125*IT_0107
      *IT_0332 + IT_0341*IT_0359 + 0.125*IT_0357*IT_0360) + IT_0361*(IT_0128*
      (IT_0172 + IT_0178) + IT_0126*(conj(IT_0172) + conj(IT_0178)) + conj
      (IT_0157)*IT_0362 + IT_0157*IT_0363) + (IT_0109 + IT_0128*IT_0157 +
       IT_0126*conj(IT_0157))*IT_0364 + IT_0365*(IT_0130 + conj(IT_0157)*IT_0366
       + IT_0157*IT_0367) + IT_0386 + IT_0387*IT_0390 + (IT_0198*IT_0392 + 0.5
      *IT_0255*IT_0399 + 3*IT_0400*(IT_0054 + 0.166666666666667*IT_0401) + 3
      *IT_0215*(IT_0019 + 0.166666666666667*IT_0402))*IT_0403 + IT_0330*IT_0411 
      + (-8)*IT_0288*IT_0416 + (-4)*conj(IT_0082)*(IT_0349*(IT_0016 + 2*IT_0198 
      + IT_0350) + -1./4*IT_0255*IT_0399 + -1./4*IT_0393*IT_0417 + -IT_0418*
      (IT_0019 + IT_0419)) + IT_0263*(IT_0253*(12*conj(IT_0082) + IT_0284) +
       IT_0080*IT_0412 + IT_0154*IT_0413 + 4*(0.25*IT_0348 + IT_0359)*IT_0418 +
       IT_0345*IT_0421 + (-4)*IT_0349*(conj(IT_0019) + IT_0422)) + (-4)*IT_0082*
      (IT_0349*IT_0359 + -IT_0418*(conj(IT_0019) + IT_0422)) + 2*IT_0119*(conj
      (IT_0198)*IT_0392 + 0.5*IT_0285*IT_0399 + 3*IT_0267*IT_0403 + 3*IT_0215*
      (conj(IT_0019) + 0.166666666666667*IT_0423) + 3*IT_0400*(conj(IT_0054) +
       IT_0425)) + IT_0266*IT_0433 + IT_0131*IT_0443 + 2*conj(IT_0198)*IT_0446 +
       conj(IT_0174)*IT_0453 + IT_0174*(IT_0088*IT_0170 + IT_0094*IT_0412 +
       IT_0160*IT_0413 + IT_0361*((-16)*conj(IT_0096) + IT_0415) + IT_0104
      *IT_0426 + IT_0364*IT_0427 + IT_0389*IT_0448 + IT_0452*(conj(IT_0157) +
       IT_0454)) + IT_0198*IT_0455 + conj(IT_0266)*(IT_0253*IT_0417 + IT_0429*
      (IT_0075 + 0.166666666666667*IT_0440) + IT_0430*(IT_0016 + IT_0442) +
       IT_0117*(IT_0362 + IT_0447) + IT_0365*(IT_0362 + IT_0449) + IT_0456 +
       IT_0457) + IT_0458;
    return create_ccomplex_return(IT_0459);
}

