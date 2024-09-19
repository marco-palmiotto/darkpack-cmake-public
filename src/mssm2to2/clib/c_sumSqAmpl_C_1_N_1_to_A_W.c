#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_1_to_A_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_1_to_A_W(
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
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + (-0.5)*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = IT_0001*IT_0010;
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0013;
    const ccomplex_t IT_0015 = 2*IT_0001;
    const ccomplex_t IT_0016 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0017 = IT_0003*IT_0016;
    const ccomplex_t IT_0018 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = 1.4142135623731*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0017 + 0.5*IT_0020);
    const ccomplex_t IT_0022 = -IT_0021;
    const ccomplex_t IT_0023 = pow(m_C_1, 2);
    const ccomplex_t IT_0024 = pow(m_N_1, 2);
    const ccomplex_t IT_0025 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0023 + -IT_0024 + -reg_prop, -1);
    const ccomplex_t IT_0026 = IT_0022*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = IT_0015*IT_0027;
    const ccomplex_t IT_0029 = pow(m_W, -2);
    const ccomplex_t IT_0030 = pow(s_34, 2);
    const ccomplex_t IT_0031 = m_C_1*m_N_1*IT_0029*IT_0030;
    const ccomplex_t IT_0032 = 4*IT_0031;
    const ccomplex_t IT_0033 = IT_0028*IT_0032;
    const ccomplex_t IT_0034 = IT_0001*IT_0027;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = (-2)*IT_0039;
    const ccomplex_t IT_0041 = cos(theta_W);
    const ccomplex_t IT_0042 = cpow(IT_0041, -1);
    const ccomplex_t IT_0043 = IT_0036*IT_0042;
    const ccomplex_t IT_0044 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0048 = IT_0037*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = IT_0046 + IT_0049;
    const ccomplex_t IT_0051 = IT_0040 + IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = (-0.5)*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = pow(m_W, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + -IT_0023 + -IT_0024 + IT_0056 
      + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = m_C_1*IT_0001;
    const ccomplex_t IT_0061 = IT_0010*IT_0025;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = IT_0059 + -IT_0063;
    const ccomplex_t IT_0065 = pow(s_14, 2);
    const ccomplex_t IT_0066 = IT_0029*IT_0065;
    const ccomplex_t IT_0067 = -IT_0023;
    const ccomplex_t IT_0068 = IT_0066 + IT_0067;
    const ccomplex_t IT_0069 = m_N_1*IT_0068;
    const ccomplex_t IT_0070 = 2*IT_0069;
    const ccomplex_t IT_0071 = s_14*s_24*IT_0029;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = s_12 + IT_0072;
    const ccomplex_t IT_0074 = m_N_1*IT_0073;
    const ccomplex_t IT_0075 = (-2)*IT_0074;
    const ccomplex_t IT_0076 = IT_0070 + IT_0075;
    const ccomplex_t IT_0077 = IT_0027*IT_0060;
    const ccomplex_t IT_0078 = pow(s_24, 2);
    const ccomplex_t IT_0079 = IT_0029*IT_0078;
    const ccomplex_t IT_0080 = -IT_0024;
    const ccomplex_t IT_0081 = IT_0079 + IT_0080;
    const ccomplex_t IT_0082 = m_C_1*IT_0081;
    const ccomplex_t IT_0083 = 2*IT_0082;
    const ccomplex_t IT_0084 = m_C_1*IT_0073;
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = IT_0001*IT_0022;
    const ccomplex_t IT_0088 = IT_0000*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = (-2)*IT_0001;
    const ccomplex_t IT_0091 = IT_0027*IT_0090;
    const ccomplex_t IT_0092 = IT_0089 + -IT_0091;
    const ccomplex_t IT_0093 = IT_0023*IT_0073;
    const ccomplex_t IT_0094 = 2*IT_0093;
    const ccomplex_t IT_0095 = IT_0023*IT_0081;
    const ccomplex_t IT_0096 = (-2)*IT_0095;
    const ccomplex_t IT_0097 = IT_0094 + IT_0096;
    const ccomplex_t IT_0098 = -IT_0091;
    const ccomplex_t IT_0099 = IT_0024*IT_0073;
    const ccomplex_t IT_0100 = 2*IT_0099;
    const ccomplex_t IT_0101 = IT_0024*IT_0068;
    const ccomplex_t IT_0102 = (-2)*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = -IT_0001;
    const ccomplex_t IT_0105 = IT_0062*IT_0104;
    const ccomplex_t IT_0106 = m_C_1*m_N_1;
    const ccomplex_t IT_0107 = s_14*s_34*IT_0029;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = s_13 + IT_0108;
    const ccomplex_t IT_0110 = IT_0106*IT_0109;
    const ccomplex_t IT_0111 = (-10)*IT_0110;
    const ccomplex_t IT_0112 = s_24*s_34*IT_0029;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = s_23 + IT_0113;
    const ccomplex_t IT_0115 = IT_0106*IT_0114;
    const ccomplex_t IT_0116 = (-10)*IT_0115;
    const ccomplex_t IT_0117 = IT_0111 + IT_0116;
    const ccomplex_t IT_0118 = (-8)*IT_0110;
    const ccomplex_t IT_0119 = (-8)*IT_0115;
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = (-0.5)*IT_0089;
    const ccomplex_t IT_0122 = s_14*s_23*s_24*IT_0029;
    const ccomplex_t IT_0123 = s_13*IT_0024;
    const ccomplex_t IT_0124 = s_12*s_24*s_34*IT_0029;
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = s_13*IT_0029*IT_0078;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0122 + IT_0123 + IT_0125 + IT_0127;
    const ccomplex_t IT_0129 = (-4)*IT_0128;
    const ccomplex_t IT_0130 = s_23*IT_0029*IT_0065;
    const ccomplex_t IT_0131 = s_23*IT_0023;
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = s_12*s_14*s_34*IT_0029;
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = s_12*s_13;
    const ccomplex_t IT_0136 = 2*IT_0135;
    const ccomplex_t IT_0137 = s_13*s_14*s_24*IT_0029;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0130 + IT_0132 + IT_0134 + IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = (-4)*IT_0139;
    const ccomplex_t IT_0141 = IT_0129 + IT_0140;
    const ccomplex_t IT_0142 = IT_0027*IT_0104;
    const ccomplex_t IT_0143 = s_12*IT_0109;
    const ccomplex_t IT_0144 = (-4)*IT_0143;
    const ccomplex_t IT_0145 = s_12*IT_0114;
    const ccomplex_t IT_0146 = (-4)*IT_0145;
    const ccomplex_t IT_0147 = IT_0130 + IT_0132 + IT_0134 + IT_0137;
    const ccomplex_t IT_0148 = (-2)*IT_0147;
    const ccomplex_t IT_0149 = -IT_0123;
    const ccomplex_t IT_0150 = IT_0122 + IT_0125 + IT_0126 + IT_0149;
    const ccomplex_t IT_0151 = (-2)*IT_0150;
    const ccomplex_t IT_0152 = IT_0144 + IT_0146 + IT_0148 + IT_0151;
    const ccomplex_t IT_0153 = IT_0148 + IT_0151;
    const ccomplex_t IT_0154 = IT_0015*IT_0062;
    const ccomplex_t IT_0155 = (-2)*IT_0110;
    const ccomplex_t IT_0156 = (-2)*IT_0115;
    const ccomplex_t IT_0157 = IT_0155 + IT_0156;
    const ccomplex_t IT_0158 = IT_0062*IT_0090;
    const ccomplex_t IT_0159 = IT_0013 + -IT_0158;
    const ccomplex_t IT_0160 = IT_0068*IT_0106;
    const ccomplex_t IT_0161 = 2*IT_0160;
    const ccomplex_t IT_0162 = -IT_0158;
    const ccomplex_t IT_0163 = IT_0081*IT_0106;
    const ccomplex_t IT_0164 = 2*IT_0163;
    const ccomplex_t IT_0165 = m_N_1*IT_0001;
    const ccomplex_t IT_0166 = IT_0027*IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = 2*IT_0074;
    const ccomplex_t IT_0169 = (-2)*IT_0069;
    const ccomplex_t IT_0170 = sin(beta);
    const ccomplex_t IT_0171 = IT_0003*IT_0170;
    const ccomplex_t IT_0172 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0173 = IT_0171*IT_0172;
    const ccomplex_t IT_0174 = 2*IT_0173;
    const ccomplex_t IT_0175 = IT_0042*IT_0170;
    const ccomplex_t IT_0176 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0177 = IT_0175*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0180 = IT_0171*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = IT_0178 + IT_0181;
    const ccomplex_t IT_0183 = IT_0174 + IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = 0.5*IT_0184;
    const ccomplex_t IT_0186 = IT_0054*IT_0185;
    const ccomplex_t IT_0187 = IT_0057*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = IT_0062*IT_0165;
    const ccomplex_t IT_0190 = -IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0084;
    const ccomplex_t IT_0192 = (-2)*IT_0082;
    const ccomplex_t IT_0193 = IT_0064*IT_0076 + IT_0077*IT_0086 + IT_0092
      *IT_0097 + IT_0098*IT_0103 + IT_0105*IT_0117 + IT_0014*IT_0120 + IT_0121
      *IT_0141 + IT_0142*IT_0152 + IT_0028*IT_0153 + IT_0154*IT_0157 + IT_0159
      *IT_0161 + IT_0162*IT_0164 + IT_0167*(IT_0168 + IT_0169) + IT_0190*
      (IT_0191 + IT_0192);
    const ccomplex_t IT_0194 = IT_0064*IT_0191;
    const ccomplex_t IT_0195 = IT_0077*IT_0169;
    const ccomplex_t IT_0196 = IT_0028*IT_0155;
    const ccomplex_t IT_0197 = IT_0092*IT_0155;
    const ccomplex_t IT_0198 = IT_0098*IT_0156;
    const ccomplex_t IT_0199 = IT_0077*IT_0168;
    const ccomplex_t IT_0200 = IT_0064*IT_0192;
    const ccomplex_t IT_0201 = IT_0028*IT_0156;
    const ccomplex_t IT_0202 = IT_0001*IT_0062;
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = IT_0073*IT_0106;
    const ccomplex_t IT_0205 = (-2)*IT_0204;
    const ccomplex_t IT_0206 = IT_0161 + IT_0205;
    const ccomplex_t IT_0207 = IT_0164 + IT_0205;
    const ccomplex_t IT_0208 = IT_0120*IT_0121 + IT_0014*IT_0141 + IT_0117
      *IT_0142 + IT_0105*IT_0152 + IT_0153*IT_0154 + IT_0097*IT_0159 + IT_0103
      *IT_0162 + IT_0086*IT_0167 + IT_0076*IT_0190 + IT_0194 + IT_0195 + IT_0196
       + IT_0199 + IT_0200 + IT_0201 + IT_0092*IT_0206 + IT_0098*IT_0207;
    const ccomplex_t IT_0209 = IT_0024*IT_0109;
    const ccomplex_t IT_0210 = 2*IT_0209;
    const ccomplex_t IT_0211 = s_12*IT_0056;
    const ccomplex_t IT_0212 = s_14*s_24;
    const ccomplex_t IT_0213 = (-2)*IT_0212;
    const ccomplex_t IT_0214 = IT_0211 + IT_0213;
    const ccomplex_t IT_0215 = IT_0029*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = s_12 + IT_0216;
    const ccomplex_t IT_0218 = (-4)*s_23*IT_0217;
    const ccomplex_t IT_0219 = IT_0210 + IT_0218;
    const ccomplex_t IT_0220 = (-6)*s_23*IT_0106;
    const ccomplex_t IT_0221 = IT_0156 + IT_0220;
    const ccomplex_t IT_0222 = conj(IT_0077)*IT_0168;
    const ccomplex_t IT_0223 = IT_0083*conj(IT_0167);
    const ccomplex_t IT_0224 = s_23*IT_0056;
    const ccomplex_t IT_0225 = IT_0029*IT_0224;
    const ccomplex_t IT_0226 = (-0.25)*IT_0225;
    const ccomplex_t IT_0227 = s_23 + IT_0226;
    const ccomplex_t IT_0228 = IT_0106*IT_0227;
    const ccomplex_t IT_0229 = (-4)*IT_0228;
    const ccomplex_t IT_0230 = 2*IT_0229;
    const ccomplex_t IT_0231 = s_12*s_23*IT_0056;
    const ccomplex_t IT_0232 = s_13*IT_0078;
    const ccomplex_t IT_0233 = IT_0231 + IT_0232;
    const ccomplex_t IT_0234 = s_14*s_23*s_24;
    const ccomplex_t IT_0235 = s_12*s_24*s_34;
    const ccomplex_t IT_0236 = s_13*IT_0024*IT_0056;
    const ccomplex_t IT_0237 = -IT_0234 + -IT_0235 + (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0233 + IT_0237;
    const ccomplex_t IT_0239 = IT_0029*IT_0238;
    const ccomplex_t IT_0240 = s_12*s_23;
    const ccomplex_t IT_0241 = (-2)*IT_0240;
    const ccomplex_t IT_0242 = IT_0123 + IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = 4*IT_0242;
    const ccomplex_t IT_0244 = (-6)*s_12;
    const ccomplex_t IT_0245 = (-6)*IT_0024*IT_0106;
    const ccomplex_t IT_0246 = (-4)*s_12*IT_0217;
    const ccomplex_t IT_0247 = (-4)*IT_0024*IT_0217;
    const ccomplex_t IT_0248 = IT_0075*conj(IT_0190);
    const ccomplex_t IT_0249 = conj(IT_0064)*IT_0192;
    const ccomplex_t IT_0250 = conj(IT_0035)*IT_0205;
    const ccomplex_t IT_0251 = conj(IT_0028)*IT_0156;
    const ccomplex_t IT_0252 = IT_0151*conj(IT_0154) + conj(IT_0105)*IT_0219 +
       conj(IT_0142)*IT_0221 + IT_0222 + IT_0223 + conj(IT_0121)*IT_0230 + conj
      (IT_0014)*IT_0243 + conj(IT_0092)*IT_0106*IT_0244 + conj(IT_0098)*IT_0245 
      + conj(IT_0159)*IT_0246 + conj(IT_0162)*IT_0247 + IT_0248 + IT_0249 +
       IT_0250 + IT_0251;
    const ccomplex_t IT_0253 = s_13*s_23*IT_0056;
    const ccomplex_t IT_0254 = s_13*s_24*s_34;
    const ccomplex_t IT_0255 = (-2)*IT_0254;
    const ccomplex_t IT_0256 = IT_0253 + IT_0255;
    const ccomplex_t IT_0257 = IT_0029*IT_0256;
    const ccomplex_t IT_0258 = (-8)*IT_0257;
    const ccomplex_t IT_0259 = s_13*s_23;
    const ccomplex_t IT_0260 = 16*IT_0259;
    const ccomplex_t IT_0261 = IT_0258 + IT_0260;
    const ccomplex_t IT_0262 = s_24*s_34;
    const ccomplex_t IT_0263 = (-2)*IT_0262;
    const ccomplex_t IT_0264 = IT_0224 + IT_0263;
    const ccomplex_t IT_0265 = IT_0029*IT_0264;
    const ccomplex_t IT_0266 = (-0.5)*IT_0265;
    const ccomplex_t IT_0267 = s_23 + IT_0266;
    const ccomplex_t IT_0268 = m_C_1*IT_0267;
    const ccomplex_t IT_0269 = (-2)*IT_0268;
    const ccomplex_t IT_0270 = 2*IT_0269;
    const ccomplex_t IT_0271 = s_13*IT_0056;
    const ccomplex_t IT_0272 = s_14*s_34;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = IT_0271 + IT_0273;
    const ccomplex_t IT_0275 = IT_0029*IT_0274;
    const ccomplex_t IT_0276 = (-0.5)*IT_0275;
    const ccomplex_t IT_0277 = s_13 + IT_0276;
    const ccomplex_t IT_0278 = m_N_1*IT_0277;
    const ccomplex_t IT_0279 = (-2)*IT_0278;
    const ccomplex_t IT_0280 = 2*IT_0279;
    const ccomplex_t IT_0281 = IT_0029*IT_0271;
    const ccomplex_t IT_0282 = (-0.25)*IT_0281;
    const ccomplex_t IT_0283 = s_13 + IT_0282;
    const ccomplex_t IT_0284 = IT_0106*IT_0283;
    const ccomplex_t IT_0285 = (-4)*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0285;
    const ccomplex_t IT_0287 = conj(IT_0028)*IT_0032;
    const ccomplex_t IT_0288 = 8*IT_0031;
    const ccomplex_t IT_0289 = s_12*IT_0030;
    const ccomplex_t IT_0290 = IT_0253 + IT_0289;
    const ccomplex_t IT_0291 = s_14*s_23*s_34;
    const ccomplex_t IT_0292 = -IT_0254 + -IT_0291;
    const ccomplex_t IT_0293 = IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = IT_0029*IT_0293;
    const ccomplex_t IT_0295 = (-0.5)*IT_0294;
    const ccomplex_t IT_0296 = IT_0259 + IT_0295;
    const ccomplex_t IT_0297 = (-8)*IT_0296;
    const ccomplex_t IT_0298 = s_23*IT_0023*IT_0056;
    const ccomplex_t IT_0299 = s_24*s_34*IT_0023;
    const ccomplex_t IT_0300 = (-2)*IT_0299;
    const ccomplex_t IT_0301 = IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = IT_0029*IT_0301;
    const ccomplex_t IT_0303 = (-0.5)*IT_0302;
    const ccomplex_t IT_0304 = IT_0131 + IT_0303;
    const ccomplex_t IT_0305 = (-4)*IT_0304;
    const ccomplex_t IT_0306 = s_34*IT_0029;
    const ccomplex_t IT_0307 = s_14*s_23;
    const ccomplex_t IT_0308 = s_13*s_24;
    const ccomplex_t IT_0309 = s_12*s_34;
    const ccomplex_t IT_0310 = -IT_0308 + -IT_0309;
    const ccomplex_t IT_0311 = IT_0307 + IT_0310;
    const ccomplex_t IT_0312 = IT_0306*IT_0311;
    const ccomplex_t IT_0313 = (-4)*IT_0312;
    const ccomplex_t IT_0314 = (-2)*IT_0269;
    const ccomplex_t IT_0315 = (-2)*IT_0279;
    const ccomplex_t IT_0316 = conj(IT_0098)*IT_0230 + conj(IT_0162)*IT_0243 +
       conj(IT_0014)*IT_0261 + conj(IT_0064)*IT_0270 + conj(IT_0077)*IT_0280 +
       conj(IT_0092)*IT_0286 + IT_0287 + conj(IT_0142)*IT_0288 + conj(IT_0154)
      *IT_0297 + conj(IT_0159)*IT_0305 + conj(IT_0105)*IT_0313 + conj(IT_0167)
      *IT_0314 + conj(IT_0190)*IT_0315;
    const ccomplex_t IT_0317 = IT_0032*conj(IT_0154) + conj(IT_0162)*IT_0230 +
       conj(IT_0098)*IT_0243 + conj(IT_0121)*IT_0261 + conj(IT_0190)*IT_0270 +
       conj(IT_0167)*IT_0280 + conj(IT_0159)*IT_0286 + conj(IT_0105)*IT_0288 +
       conj(IT_0028)*IT_0297 + conj(IT_0092)*IT_0305 + conj(IT_0142)*IT_0313 +
       conj(IT_0077)*IT_0314 + conj(IT_0064)*IT_0315;
    const ccomplex_t IT_0318 = s_23*m_C_1;
    const ccomplex_t IT_0319 = 6*IT_0318;
    const ccomplex_t IT_0320 = s_13*m_N_1;
    const ccomplex_t IT_0321 = 6*IT_0320;
    const ccomplex_t IT_0322 = s_12*IT_0029*IT_0030;
    const ccomplex_t IT_0323 = 0.5*IT_0322;
    const ccomplex_t IT_0324 = s_13*s_24*s_34*IT_0029;
    const ccomplex_t IT_0325 = (-0.5)*IT_0324;
    const ccomplex_t IT_0326 = s_14*s_23*s_34*IT_0029;
    const ccomplex_t IT_0327 = (-0.5)*IT_0326;
    const ccomplex_t IT_0328 = IT_0259 + IT_0323 + IT_0325 + IT_0327;
    const ccomplex_t IT_0329 = 8*IT_0328;
    const ccomplex_t IT_0330 = (-6)*IT_0318;
    const ccomplex_t IT_0331 = (-6)*IT_0320;
    const ccomplex_t IT_0332 = conj(IT_0092)*IT_0155;
    const ccomplex_t IT_0333 = conj(IT_0098)*IT_0156;
    const ccomplex_t IT_0334 = 12*IT_0259;
    const ccomplex_t IT_0335 = IT_0032*(conj(IT_0121) + conj(IT_0142)) +
       IT_0148*conj(IT_0159) + IT_0151*conj(IT_0162) + conj(IT_0014)*IT_0297 +
       conj(IT_0064)*IT_0319 + conj(IT_0077)*IT_0321 + conj(IT_0105)*IT_0329 +
       conj(IT_0167)*IT_0330 + conj(IT_0190)*IT_0331 + IT_0332 + IT_0333 + conj
      (IT_0154)*IT_0334;
    const ccomplex_t IT_0336 = (-6)*IT_0106;
    const ccomplex_t IT_0337 = 0.166666666666667*conj(IT_0167);
    const ccomplex_t IT_0338 = m_C_1*IT_0114;
    const ccomplex_t IT_0339 = 4*IT_0338;
    const ccomplex_t IT_0340 = 0.166666666666667*conj(IT_0105);
    const ccomplex_t IT_0341 = 0.166666666666667*conj(IT_0154);
    const ccomplex_t IT_0342 = 0.166666666666667*conj(IT_0121);
    const ccomplex_t IT_0343 = 0.166666666666667*conj(IT_0092);
    const ccomplex_t IT_0344 = 0.166666666666667*conj(IT_0098);
    const ccomplex_t IT_0345 = m_N_1*IT_0109;
    const ccomplex_t IT_0346 = (-4)*IT_0345;
    const ccomplex_t IT_0347 = 0.166666666666667*conj(IT_0142);
    const ccomplex_t IT_0348 = 0.166666666666667*conj(IT_0028);
    const ccomplex_t IT_0349 = 0.166666666666667*conj(IT_0014);
    const ccomplex_t IT_0350 = s_12*conj(IT_0064) + conj(IT_0077)*IT_0106 +
       0.166666666666667*conj(IT_0190)*IT_0336 + IT_0244*IT_0337 + IT_0339
      *IT_0340 + IT_0319*IT_0341 + IT_0315*IT_0342 + IT_0070*IT_0343 + IT_0075
      *IT_0344 + IT_0346*IT_0347 + IT_0331*IT_0348 + IT_0270*IT_0349;
    const ccomplex_t IT_0351 = 4*IT_0345;
    const ccomplex_t IT_0352 = (-4)*IT_0338;
    const ccomplex_t IT_0353 = s_12*conj(IT_0077) + conj(IT_0064)*IT_0106 +
       0.166666666666667*conj(IT_0190)*IT_0244 + IT_0336*IT_0337 + IT_0321
      *IT_0341 + IT_0314*IT_0342 + IT_0085*IT_0343 + IT_0083*IT_0344 + IT_0330
      *IT_0348 + IT_0280*IT_0349 + IT_0340*IT_0351 + IT_0347*IT_0352;
    const ccomplex_t IT_0354 = 0.166666666666667*conj(IT_0064);
    const ccomplex_t IT_0355 = 0.166666666666667*conj(IT_0077);
    const ccomplex_t IT_0356 = 0.166666666666667*conj(IT_0159);
    const ccomplex_t IT_0357 = 0.166666666666667*conj(IT_0162);
    const ccomplex_t IT_0358 = s_12*conj(IT_0167) + IT_0106*conj(IT_0190) +
       IT_0330*IT_0341 + IT_0280*IT_0342 + IT_0169*IT_0343 + IT_0168*IT_0344 +
       IT_0321*IT_0348 + IT_0314*IT_0349 + IT_0347*IT_0351 + IT_0340*IT_0352 +
       IT_0244*IT_0354 + IT_0336*IT_0355 + IT_0085*IT_0356 + IT_0083*IT_0357;
    const ccomplex_t IT_0359 = IT_0106*conj(IT_0167) + s_12*conj(IT_0190) +
       IT_0331*IT_0341 + IT_0270*IT_0342 + IT_0191*IT_0343 + IT_0192*IT_0344 +
       IT_0340*IT_0346 + IT_0339*IT_0347 + IT_0319*IT_0348 + IT_0315*IT_0349 +
       IT_0336*IT_0354 + IT_0244*IT_0355 + IT_0070*IT_0356 + IT_0075*IT_0357;
    const ccomplex_t IT_0360 = 12*IT_0031;
    const ccomplex_t IT_0361 = 0.25*conj(IT_0167);
    const ccomplex_t IT_0362 = 0.25*conj(IT_0064);
    const ccomplex_t IT_0363 = 0.25*conj(IT_0154);
    const ccomplex_t IT_0364 = 0.25*conj(IT_0121);
    const ccomplex_t IT_0365 = 0.25*conj(IT_0077);
    const ccomplex_t IT_0366 = (-6)*s_13*IT_0106;
    const ccomplex_t IT_0367 = 0.25*conj(IT_0092);
    const ccomplex_t IT_0368 = 0.25*conj(IT_0098);
    const ccomplex_t IT_0369 = IT_0023*IT_0114;
    const ccomplex_t IT_0370 = 2*IT_0369;
    const ccomplex_t IT_0371 = (-4)*s_13*IT_0217;
    const ccomplex_t IT_0372 = IT_0370 + IT_0371;
    const ccomplex_t IT_0373 = 0.25*conj(IT_0159);
    const ccomplex_t IT_0374 = 0.25*conj(IT_0162);
    const ccomplex_t IT_0375 = 0.25*conj(IT_0014);
    const ccomplex_t IT_0376 = 0.25*IT_0329;
    const ccomplex_t IT_0377 = IT_0287 + IT_0332 + IT_0333 + conj(IT_0190)
      *IT_0346 + conj(IT_0142)*IT_0360 + 4*IT_0352*IT_0361 + 4*IT_0339*IT_0362 +
       4*IT_0329*IT_0363 + 4*IT_0288*IT_0364 + 4*IT_0351*IT_0365 + 4*IT_0366
      *IT_0367 + 4*IT_0220*IT_0368 + 4*IT_0372*IT_0373 + 4*IT_0219*IT_0374 + 4
      *IT_0313*IT_0375 + 4*conj(IT_0105)*(IT_0322 + IT_0376);
    const ccomplex_t IT_0378 = IT_0155 + IT_0366;
    const ccomplex_t IT_0379 = conj(IT_0028)*IT_0329 + conj(IT_0190)*IT_0339 +
       conj(IT_0105)*IT_0360 + 4*IT_0351*IT_0361 + 4*IT_0346*IT_0362 + 4*IT_0032
      *IT_0363 + 4*IT_0313*IT_0364 + 4*IT_0352*IT_0365 + 4*IT_0219*IT_0368 + 4
      *IT_0367*IT_0372 + 4*IT_0221*IT_0374 + 4*IT_0288*IT_0375 + 4*conj(IT_0142)
      *(IT_0322 + IT_0376) + 4*IT_0373*IT_0378;
    const ccomplex_t IT_0380 = conj(IT_0064)*IT_0191;
    const ccomplex_t IT_0381 = IT_0070*conj(IT_0190);
    const ccomplex_t IT_0382 = (-6)*IT_0023*IT_0106;
    const ccomplex_t IT_0383 = (-4)*IT_0023*IT_0217;
    const ccomplex_t IT_0384 = conj(IT_0077)*IT_0169;
    const ccomplex_t IT_0385 = conj(IT_0028)*IT_0155;
    const ccomplex_t IT_0386 = IT_0085*conj(IT_0167);
    const ccomplex_t IT_0387 = IT_0148*conj(IT_0154) + conj(IT_0098)*IT_0106
      *IT_0244 + conj(IT_0162)*IT_0246 + IT_0250 + conj(IT_0121)*IT_0286 + conj
      (IT_0014)*IT_0305 + conj(IT_0105)*IT_0372 + conj(IT_0142)*IT_0378 +
       IT_0380 + IT_0381 + conj(IT_0092)*IT_0382 + conj(IT_0159)*IT_0383 +
       IT_0384 + IT_0385 + IT_0386;
    const ccomplex_t IT_0388 = s_12*IT_0068;
    const ccomplex_t IT_0389 = 4*IT_0388;
    const ccomplex_t IT_0390 = s_12*IT_0081;
    const ccomplex_t IT_0391 = 4*IT_0390;
    const ccomplex_t IT_0392 = s_12*IT_0073;
    const ccomplex_t IT_0393 = (-8)*IT_0392;
    const ccomplex_t IT_0394 = IT_0389 + IT_0391 + IT_0393;
    const ccomplex_t IT_0395 = 8*IT_0160;
    const ccomplex_t IT_0396 = 8*IT_0163;
    const ccomplex_t IT_0397 = (-16)*IT_0204;
    const ccomplex_t IT_0398 = IT_0395 + IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = IT_0120*conj(IT_0121) + conj(IT_0014)*IT_0141 +
       IT_0117*conj(IT_0142) + conj(IT_0105)*IT_0152 + IT_0153*conj(IT_0154) +
       IT_0097*conj(IT_0159) + IT_0103*conj(IT_0162) + conj(IT_0092)*IT_0206 +
       conj(IT_0098)*IT_0207 + IT_0222 + IT_0223 + IT_0248 + IT_0249 + IT_0251 +
       IT_0380 + IT_0381 + IT_0384 + IT_0385 + IT_0386 + conj(IT_0203)*IT_0394 +
       conj(IT_0035)*IT_0398;
    const ccomplex_t IT_0400 = conj(IT_0077)*IT_0083;
    const ccomplex_t IT_0401 = conj(IT_0190)*IT_0192;
    const ccomplex_t IT_0402 = conj(IT_0064)*IT_0075;
    const ccomplex_t IT_0403 = conj(IT_0028)*IT_0151 + conj(IT_0167)*IT_0168 +
       conj(IT_0142)*IT_0219 + conj(IT_0105)*IT_0220 + conj(IT_0014)*IT_0230 +
       conj(IT_0121)*IT_0243 + IT_0106*conj(IT_0159)*IT_0244 + conj(IT_0162)
      *IT_0245 + conj(IT_0092)*IT_0246 + conj(IT_0098)*IT_0247 + IT_0400 +
       IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = conj(IT_0064)*IT_0070;
    const ccomplex_t IT_0405 = conj(IT_0167)*IT_0169;
    const ccomplex_t IT_0406 = conj(IT_0077)*IT_0085;
    const ccomplex_t IT_0407 = conj(IT_0028)*IT_0148 + conj(IT_0190)*IT_0191 +
       IT_0106*conj(IT_0162)*IT_0244 + conj(IT_0098)*IT_0246 + conj(IT_0014)
      *IT_0286 + conj(IT_0121)*IT_0305 + conj(IT_0105)*IT_0366 + conj(IT_0142)
      *IT_0372 + conj(IT_0159)*IT_0382 + conj(IT_0092)*IT_0383 + IT_0404 +
       IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = conj(IT_0092)*IT_0097 + conj(IT_0098)*IT_0103 +
       conj(IT_0105)*IT_0117 + conj(IT_0014)*IT_0120 + conj(IT_0121)*IT_0141 +
       conj(IT_0142)*IT_0152 + conj(IT_0028)*IT_0153 + conj(IT_0154)*IT_0157 +
       conj(IT_0167)*IT_0168 + conj(IT_0190)*IT_0191 + conj(IT_0159)*IT_0206 +
       conj(IT_0162)*IT_0207 + conj(IT_0035)*IT_0394 + conj(IT_0203)*IT_0398 +
       IT_0400 + IT_0401 + IT_0402 + IT_0404 + IT_0405 + IT_0406;
    const ccomplex_t IT_0409 = conj(IT_0014)*IT_0033 + conj(IT_0035)*IT_0193 +
       conj(IT_0159)*(IT_0194 + IT_0195 + IT_0196) + conj(IT_0154)*(IT_0197 +
       IT_0198) + conj(IT_0105)*(IT_0033 + IT_0197 + IT_0198) + conj(IT_0162)*
      (IT_0199 + IT_0200 + IT_0201) + conj(IT_0203)*IT_0208 + IT_0162*IT_0252 +
       IT_0014*IT_0316 + IT_0121*IT_0317 + IT_0028*(conj(IT_0092)*IT_0148 + conj
      (IT_0098)*IT_0151 + conj(IT_0121)*IT_0297 + conj(IT_0190)*IT_0319 + conj
      (IT_0167)*IT_0321 + conj(IT_0142)*IT_0329 + conj(IT_0077)*IT_0330 + conj
      (IT_0064)*IT_0331 + conj(IT_0028)*IT_0334) + IT_0154*IT_0335 + 6*IT_0064
      *IT_0350 + 6*IT_0077*IT_0353 + 6*IT_0167*IT_0358 + 6*IT_0190*IT_0359 +
       IT_0105*IT_0377 + IT_0142*IT_0379 + IT_0159*IT_0387 + IT_0203*IT_0399 +
       IT_0098*IT_0403 + IT_0092*IT_0407 + IT_0035*IT_0408;
    return create_ccomplex_return(IT_0409);
}

