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
    const ccomplex_t IT_0000 = pow(m_W, -2);
    const ccomplex_t IT_0001 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0002 = -IT_0001;
    const ccomplex_t IT_0003 = s_12 + IT_0002;
    const ccomplex_t IT_0004 = m_C_1*IT_0003;
    const ccomplex_t IT_0005 = pow(m_W, 2);
    const ccomplex_t IT_0006 = pow(m_C_1, 2);
    const ccomplex_t IT_0007 = pow(m_N_1, 2);
    const ccomplex_t IT_0008 = cpow((-2)*s_12 + IT_0005 + -IT_0006 + -IT_0007 
      + -reg_prop, -1);
    const ccomplex_t IT_0009 = sin(beta);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 2*IT_0014;
    const ccomplex_t IT_0016 = cos(theta_W);
    const ccomplex_t IT_0017 = cpow(IT_0016, -1);
    const ccomplex_t IT_0018 = IT_0009*IT_0017;
    const ccomplex_t IT_0019 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0020 = IT_0018*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = IT_0021 + IT_0024;
    const ccomplex_t IT_0026 = IT_0015 + IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*m_W*e_em;
    const ccomplex_t IT_0030 = IT_0028*IT_0029;
    const ccomplex_t IT_0031 = IT_0008*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0034 = m_N_1*IT_0033;
    const ccomplex_t IT_0035 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0036 = IT_0011*IT_0035;
    const ccomplex_t IT_0037 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0038 = IT_0011*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + (-0.5)*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0006 + -IT_0007 + -reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0034*IT_0044;
    const ccomplex_t IT_0046 = -IT_0032 + IT_0045;
    const ccomplex_t IT_0047 = pow(s_14, 2);
    const ccomplex_t IT_0048 = IT_0000*IT_0047;
    const ccomplex_t IT_0049 = -IT_0006;
    const ccomplex_t IT_0050 = IT_0048 + IT_0049;
    const ccomplex_t IT_0051 = m_N_1*IT_0050;
    const ccomplex_t IT_0052 = m_C_1*IT_0033;
    const ccomplex_t IT_0053 = IT_0044*IT_0052;
    const ccomplex_t IT_0054 = cos(beta);
    const ccomplex_t IT_0055 = IT_0011*IT_0054;
    const ccomplex_t IT_0056 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (-2)*IT_0057;
    const ccomplex_t IT_0059 = IT_0017*IT_0054;
    const ccomplex_t IT_0060 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0061 = IT_0059*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0064 = IT_0055*IT_0063;
    const ccomplex_t IT_0065 = 1.4142135623731*IT_0064;
    const ccomplex_t IT_0066 = IT_0062 + IT_0065;
    const ccomplex_t IT_0067 = IT_0058 + IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = IT_0029*IT_0069;
    const ccomplex_t IT_0071 = IT_0008*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = -IT_0053 + IT_0072;
    const ccomplex_t IT_0074 = s_24*s_34*IT_0000;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = s_23 + IT_0075;
    const ccomplex_t IT_0077 = IT_0006*IT_0076;
    const ccomplex_t IT_0078 = -IT_0033;
    const ccomplex_t IT_0079 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0080 = IT_0011*IT_0079;
    const ccomplex_t IT_0081 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0082 = IT_0011*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + 0.5*IT_0083);
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = IT_0042*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0078*IT_0087;
    const ccomplex_t IT_0089 = m_C_1*m_N_1;
    const ccomplex_t IT_0090 = s_13*IT_0005;
    const ccomplex_t IT_0091 = IT_0000*IT_0090;
    const ccomplex_t IT_0092 = (-0.25)*IT_0091;
    const ccomplex_t IT_0093 = s_13 + IT_0092;
    const ccomplex_t IT_0094 = IT_0089*IT_0093;
    const ccomplex_t IT_0095 = (-4)*IT_0094;
    const ccomplex_t IT_0096 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_1
      *Gamma_c1 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0033*IT_0041;
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0033*IT_0085;
    const ccomplex_t IT_0102 = IT_0096*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 2*IT_0103;
    const ccomplex_t IT_0105 = 2*IT_0033;
    const ccomplex_t IT_0106 = IT_0087*IT_0105;
    const ccomplex_t IT_0107 = pow(s_34, 2);
    const ccomplex_t IT_0108 = m_C_1*m_N_1*IT_0000*IT_0107;
    const ccomplex_t IT_0109 = 4*IT_0108;
    const ccomplex_t IT_0110 = IT_0106*IT_0109;
    const ccomplex_t IT_0111 = (-0.5)*IT_0103;
    const ccomplex_t IT_0112 = s_23*IT_0006;
    const ccomplex_t IT_0113 = s_23*IT_0005*IT_0006;
    const ccomplex_t IT_0114 = s_24*s_34*IT_0006;
    const ccomplex_t IT_0115 = (-2)*IT_0114;
    const ccomplex_t IT_0116 = IT_0113 + IT_0115;
    const ccomplex_t IT_0117 = IT_0000*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = IT_0112 + IT_0118;
    const ccomplex_t IT_0120 = (-2)*IT_0033;
    const ccomplex_t IT_0121 = IT_0044*IT_0120;
    const ccomplex_t IT_0122 = -IT_0121;
    const ccomplex_t IT_0123 = (-6)*s_12;
    const ccomplex_t IT_0124 = IT_0044*IT_0078;
    const ccomplex_t IT_0125 = (-6)*s_13;
    const ccomplex_t IT_0126 = IT_0099 + -IT_0121;
    const ccomplex_t IT_0127 = (-6)*IT_0006;
    const ccomplex_t IT_0128 = s_12*IT_0005;
    const ccomplex_t IT_0129 = s_14*s_24;
    const ccomplex_t IT_0130 = (-2)*IT_0129;
    const ccomplex_t IT_0131 = IT_0128 + IT_0130;
    const ccomplex_t IT_0132 = IT_0000*IT_0131;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = s_12 + IT_0133;
    const ccomplex_t IT_0135 = (-4)*s_13;
    const ccomplex_t IT_0136 = (-4)*IT_0006;
    const ccomplex_t IT_0137 = IT_0052*IT_0087;
    const ccomplex_t IT_0138 = IT_0034*IT_0087;
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = s_23*IT_0000*IT_0047;
    const ccomplex_t IT_0141 = s_13*s_14*s_24*IT_0000;
    const ccomplex_t IT_0142 = -IT_0112;
    const ccomplex_t IT_0143 = s_12*s_14*s_34*IT_0000;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0140 + IT_0141 + IT_0142 + IT_0144;
    const ccomplex_t IT_0146 = (-2)*IT_0103;
    const ccomplex_t IT_0147 = IT_0044*IT_0105;
    const ccomplex_t IT_0148 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = s_13 + IT_0149;
    const ccomplex_t IT_0151 = IT_0089*IT_0150;
    const ccomplex_t IT_0152 = IT_0146*IT_0151;
    const ccomplex_t IT_0153 = 2*IT_0004;
    const ccomplex_t IT_0154 = IT_0073*IT_0153;
    const ccomplex_t IT_0155 = (-2)*IT_0004;
    const ccomplex_t IT_0156 = IT_0139*IT_0155;
    const ccomplex_t IT_0157 = (-2)*IT_0151;
    const ccomplex_t IT_0158 = IT_0106*IT_0157;
    const ccomplex_t IT_0159 = IT_0033*IT_0087;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0087*IT_0120;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = IT_0003*IT_0006;
    const ccomplex_t IT_0164 = 2*IT_0163;
    const ccomplex_t IT_0165 = IT_0003*IT_0007;
    const ccomplex_t IT_0166 = 2*IT_0165;
    const ccomplex_t IT_0167 = pow(s_24, 2);
    const ccomplex_t IT_0168 = IT_0000*IT_0167;
    const ccomplex_t IT_0169 = -IT_0007;
    const ccomplex_t IT_0170 = IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = IT_0006*IT_0170;
    const ccomplex_t IT_0172 = (-2)*IT_0171;
    const ccomplex_t IT_0173 = IT_0007*IT_0050;
    const ccomplex_t IT_0174 = (-2)*IT_0173;
    const ccomplex_t IT_0175 = IT_0164 + IT_0166 + IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = s_12*IT_0050;
    const ccomplex_t IT_0177 = 4*IT_0176;
    const ccomplex_t IT_0178 = s_12*IT_0170;
    const ccomplex_t IT_0179 = 4*IT_0178;
    const ccomplex_t IT_0180 = s_12*IT_0003;
    const ccomplex_t IT_0181 = (-8)*IT_0180;
    const ccomplex_t IT_0182 = IT_0177 + IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = IT_0033*IT_0044;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = IT_0050*IT_0089;
    const ccomplex_t IT_0186 = 8*IT_0185;
    const ccomplex_t IT_0187 = IT_0089*IT_0170;
    const ccomplex_t IT_0188 = 8*IT_0187;
    const ccomplex_t IT_0189 = IT_0003*IT_0089;
    const ccomplex_t IT_0190 = (-16)*IT_0189;
    const ccomplex_t IT_0191 = IT_0186 + IT_0188 + IT_0190;
    const ccomplex_t IT_0192 = (-10)*IT_0151;
    const ccomplex_t IT_0193 = IT_0076*IT_0089;
    const ccomplex_t IT_0194 = (-10)*IT_0193;
    const ccomplex_t IT_0195 = IT_0192 + IT_0194;
    const ccomplex_t IT_0196 = (-8)*IT_0151;
    const ccomplex_t IT_0197 = (-8)*IT_0193;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = s_14*s_23*s_24*IT_0000;
    const ccomplex_t IT_0200 = s_13*IT_0007;
    const ccomplex_t IT_0201 = s_12*s_24*s_34*IT_0000;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = s_13*IT_0000*IT_0167;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = IT_0199 + IT_0200 + IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = (-4)*IT_0205;
    const ccomplex_t IT_0207 = s_12*s_13;
    const ccomplex_t IT_0208 = 2*IT_0207;
    const ccomplex_t IT_0209 = -IT_0141;
    const ccomplex_t IT_0210 = IT_0140 + IT_0142 + IT_0144 + IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = (-4)*IT_0210;
    const ccomplex_t IT_0212 = IT_0206 + IT_0211;
    const ccomplex_t IT_0213 = s_12*IT_0150;
    const ccomplex_t IT_0214 = (-4)*IT_0213;
    const ccomplex_t IT_0215 = s_12*IT_0076;
    const ccomplex_t IT_0216 = (-4)*IT_0215;
    const ccomplex_t IT_0217 = (-2)*IT_0145;
    const ccomplex_t IT_0218 = -IT_0200;
    const ccomplex_t IT_0219 = IT_0199 + IT_0202 + IT_0203 + IT_0218;
    const ccomplex_t IT_0220 = (-2)*IT_0219;
    const ccomplex_t IT_0221 = IT_0214 + IT_0216 + IT_0217 + IT_0220;
    const ccomplex_t IT_0222 = IT_0217 + IT_0220;
    const ccomplex_t IT_0223 = (-2)*IT_0193;
    const ccomplex_t IT_0224 = IT_0157 + IT_0223;
    const ccomplex_t IT_0225 = 2*conj(IT_0103);
    const ccomplex_t IT_0226 = conj(IT_0046)*IT_0153;
    const ccomplex_t IT_0227 = 2*IT_0051;
    const ccomplex_t IT_0228 = conj(IT_0073)*IT_0227;
    const ccomplex_t IT_0229 = m_C_1*IT_0170;
    const ccomplex_t IT_0230 = 2*IT_0229;
    const ccomplex_t IT_0231 = conj(IT_0137)*IT_0230;
    const ccomplex_t IT_0232 = m_N_1*IT_0003;
    const ccomplex_t IT_0233 = 2*IT_0232;
    const ccomplex_t IT_0234 = conj(IT_0139)*IT_0233;
    const ccomplex_t IT_0235 = (-2)*conj(IT_0103);
    const ccomplex_t IT_0236 = conj(IT_0137)*IT_0155;
    const ccomplex_t IT_0237 = (-2)*IT_0051;
    const ccomplex_t IT_0238 = conj(IT_0139)*IT_0237;
    const ccomplex_t IT_0239 = (-2)*IT_0229;
    const ccomplex_t IT_0240 = conj(IT_0046)*IT_0239;
    const ccomplex_t IT_0241 = (-2)*IT_0232;
    const ccomplex_t IT_0242 = conj(IT_0073)*IT_0241;
    const ccomplex_t IT_0243 = 2*IT_0185;
    const ccomplex_t IT_0244 = (-2)*IT_0189;
    const ccomplex_t IT_0245 = 2*IT_0187;
    const ccomplex_t IT_0246 = conj(IT_0162)*IT_0175 + conj(IT_0160)*IT_0182 +
       conj(IT_0184)*IT_0191 + conj(IT_0124)*IT_0195 + conj(IT_0100)*IT_0198 +
       conj(IT_0111)*IT_0212 + conj(IT_0088)*IT_0221 + conj(IT_0106)*IT_0222 +
       conj(IT_0147)*IT_0224 + IT_0163*IT_0225 + IT_0226 + IT_0228 + IT_0231 +
       IT_0234 + IT_0171*IT_0235 + IT_0236 + IT_0238 + IT_0240 + IT_0242 + conj
      (IT_0126)*(IT_0243 + IT_0244) + conj(IT_0122)*(IT_0244 + IT_0245);
    const ccomplex_t IT_0247 = 2*IT_0077;
    const ccomplex_t IT_0248 = IT_0134*IT_0135;
    const ccomplex_t IT_0249 = IT_0247 + IT_0248;
    const ccomplex_t IT_0250 = IT_0089*IT_0125;
    const ccomplex_t IT_0251 = IT_0157 + IT_0250;
    const ccomplex_t IT_0252 = conj(IT_0073)*IT_0153;
    const ccomplex_t IT_0253 = conj(IT_0046)*IT_0227;
    const ccomplex_t IT_0254 = 2*IT_0095;
    const ccomplex_t IT_0255 = (-4)*s_12;
    const ccomplex_t IT_0256 = IT_0134*IT_0255;
    const ccomplex_t IT_0257 = IT_0134*IT_0136;
    const ccomplex_t IT_0258 = (-4)*IT_0119;
    const ccomplex_t IT_0259 = conj(IT_0139)*IT_0155;
    const ccomplex_t IT_0260 = conj(IT_0137)*IT_0237;
    const ccomplex_t IT_0261 = conj(IT_0106)*IT_0157;
    const ccomplex_t IT_0262 = conj(IT_0160)*IT_0244;
    const ccomplex_t IT_0263 = IT_0089*conj(IT_0103)*IT_0127 + conj(IT_0147)
      *IT_0217 + conj(IT_0124)*IT_0249 + conj(IT_0088)*IT_0251 + IT_0252 +
       IT_0253 + conj(IT_0111)*IT_0254 + conj(IT_0122)*IT_0256 + conj(IT_0126)
      *IT_0257 + conj(IT_0100)*IT_0258 + IT_0259 + IT_0260 + IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = s_23*IT_0005;
    const ccomplex_t IT_0265 = IT_0000*IT_0264;
    const ccomplex_t IT_0266 = (-0.25)*IT_0265;
    const ccomplex_t IT_0267 = s_23 + IT_0266;
    const ccomplex_t IT_0268 = IT_0089*IT_0267;
    const ccomplex_t IT_0269 = (-4)*IT_0268;
    const ccomplex_t IT_0270 = 2*IT_0269;
    const ccomplex_t IT_0271 = IT_0254 + IT_0270;
    const ccomplex_t IT_0272 = s_12*s_23*IT_0005;
    const ccomplex_t IT_0273 = s_13*IT_0167;
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = s_14*s_23*s_24;
    const ccomplex_t IT_0276 = s_12*s_24*s_34;
    const ccomplex_t IT_0277 = s_13*IT_0005*IT_0007;
    const ccomplex_t IT_0278 = -IT_0275 + -IT_0276 + (-0.5)*IT_0277;
    const ccomplex_t IT_0279 = IT_0274 + IT_0278;
    const ccomplex_t IT_0280 = IT_0000*IT_0279;
    const ccomplex_t IT_0281 = s_12*s_23;
    const ccomplex_t IT_0282 = (-2)*IT_0281;
    const ccomplex_t IT_0283 = IT_0200 + IT_0280 + IT_0282;
    const ccomplex_t IT_0284 = 4*IT_0283;
    const ccomplex_t IT_0285 = IT_0258 + IT_0284;
    const ccomplex_t IT_0286 = IT_0136 + IT_0255;
    const ccomplex_t IT_0287 = IT_0089*IT_0127;
    const ccomplex_t IT_0288 = (-6)*IT_0007*IT_0089;
    const ccomplex_t IT_0289 = IT_0007*IT_0150;
    const ccomplex_t IT_0290 = 2*IT_0289;
    const ccomplex_t IT_0291 = (-4)*s_23*IT_0134;
    const ccomplex_t IT_0292 = (-4)*IT_0007*IT_0134;
    const ccomplex_t IT_0293 = (-4)*IT_0189;
    const ccomplex_t IT_0294 = (-6)*s_23*IT_0089;
    const ccomplex_t IT_0295 = IT_0089*IT_0123*(conj(IT_0122) + conj(IT_0126))
       + conj(IT_0160)*IT_0175 + conj(IT_0106)*IT_0222 + conj(IT_0147)*IT_0224 +
       IT_0226 + IT_0228 + IT_0231 + IT_0234 + IT_0236 + IT_0238 + IT_0240 +
       IT_0242 + conj(IT_0100)*IT_0271 + conj(IT_0111)*IT_0285 + (-8)*IT_0134*
      (s_12*conj(IT_0162) + (-0.125)*conj(IT_0103)*IT_0286) + conj(IT_0126)
      *IT_0287 + conj(IT_0122)*IT_0288 + conj(IT_0088)*(IT_0247 + IT_0248 +
       IT_0290 + IT_0291) + conj(IT_0162)*(IT_0257 + IT_0292) + conj(IT_0184)*
      (IT_0243 + IT_0245 + IT_0293) + conj(IT_0124)*(IT_0157 + IT_0223 + IT_0250
       + IT_0294);
    const ccomplex_t IT_0296 = IT_0139*IT_0230;
    const ccomplex_t IT_0297 = IT_0073*IT_0239;
    const ccomplex_t IT_0298 = IT_0106*IT_0223;
    const ccomplex_t IT_0299 = IT_0227 + IT_0241;
    const ccomplex_t IT_0300 = IT_0233 + IT_0237;
    const ccomplex_t IT_0301 = IT_0164 + IT_0172;
    const ccomplex_t IT_0302 = IT_0166 + IT_0174;
    const ccomplex_t IT_0303 = IT_0154 + IT_0156 + IT_0158 + IT_0104*IT_0185 +
       IT_0146*IT_0189 + IT_0088*IT_0195 + IT_0111*IT_0198 + IT_0100*IT_0212 +
       IT_0124*IT_0221 + IT_0147*IT_0222 + IT_0296 + IT_0297 + IT_0298 + IT_0046
      *IT_0299 + IT_0137*IT_0300 + IT_0126*IT_0301 + IT_0122*IT_0302;
    const ccomplex_t IT_0304 = conj(IT_0139)*IT_0230;
    const ccomplex_t IT_0305 = conj(IT_0137)*IT_0233;
    const ccomplex_t IT_0306 = conj(IT_0073)*IT_0239;
    const ccomplex_t IT_0307 = conj(IT_0046)*IT_0241;
    const ccomplex_t IT_0308 = conj(IT_0106)*IT_0223;
    const ccomplex_t IT_0309 = IT_0182*conj(IT_0184) + conj(IT_0160)*IT_0191 +
       conj(IT_0088)*IT_0195 + conj(IT_0111)*IT_0198 + conj(IT_0100)*IT_0212 +
       conj(IT_0124)*IT_0221 + conj(IT_0147)*IT_0222 + IT_0185*IT_0225 + IT_0189
      *IT_0235 + IT_0252 + IT_0253 + IT_0259 + IT_0260 + IT_0261 + conj(IT_0162)
      *(IT_0243 + IT_0245 + IT_0293) + conj(IT_0126)*IT_0301 + conj(IT_0122)
      *IT_0302 + IT_0304 + IT_0305 + IT_0306 + IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = IT_0290 + IT_0291;
    const ccomplex_t IT_0311 = IT_0223 + IT_0294;
    const ccomplex_t IT_0312 = IT_0089*conj(IT_0103)*IT_0123 + conj(IT_0147)
      *IT_0220 + conj(IT_0126)*IT_0256 + IT_0262 + conj(IT_0111)*IT_0270 + conj
      (IT_0100)*IT_0284 + conj(IT_0122)*IT_0292 + IT_0304 + IT_0305 + IT_0306 +
       IT_0307 + IT_0308 + conj(IT_0124)*IT_0310 + conj(IT_0088)*IT_0311;
    const ccomplex_t IT_0313 = s_24*s_34;
    const ccomplex_t IT_0314 = (-2)*IT_0313;
    const ccomplex_t IT_0315 = IT_0264 + IT_0314;
    const ccomplex_t IT_0316 = IT_0000*IT_0315;
    const ccomplex_t IT_0317 = (-0.5)*IT_0316;
    const ccomplex_t IT_0318 = s_23 + IT_0317;
    const ccomplex_t IT_0319 = m_C_1*IT_0318;
    const ccomplex_t IT_0320 = (-2)*IT_0319;
    const ccomplex_t IT_0321 = 2*IT_0320;
    const ccomplex_t IT_0322 = m_C_1*IT_0076;
    const ccomplex_t IT_0323 = 4*IT_0322;
    const ccomplex_t IT_0324 = 6*s_12;
    const ccomplex_t IT_0325 = 6*IT_0089;
    const ccomplex_t IT_0326 = s_23*m_C_1;
    const ccomplex_t IT_0327 = 6*IT_0326;
    const ccomplex_t IT_0328 = (-6)*IT_0089;
    const ccomplex_t IT_0329 = s_13*m_N_1;
    const ccomplex_t IT_0330 = (-6)*IT_0329;
    const ccomplex_t IT_0331 = m_N_1*IT_0150;
    const ccomplex_t IT_0332 = (-4)*IT_0331;
    const ccomplex_t IT_0333 = s_14*s_34;
    const ccomplex_t IT_0334 = (-2)*IT_0333;
    const ccomplex_t IT_0335 = IT_0090 + IT_0334;
    const ccomplex_t IT_0336 = IT_0000*IT_0335;
    const ccomplex_t IT_0337 = (-0.5)*IT_0336;
    const ccomplex_t IT_0338 = s_13 + IT_0337;
    const ccomplex_t IT_0339 = m_N_1*IT_0338;
    const ccomplex_t IT_0340 = (-2)*IT_0339;
    const ccomplex_t IT_0341 = (-2)*IT_0340;
    const ccomplex_t IT_0342 = IT_0123*conj(IT_0139) + IT_0051*IT_0225 + conj
      (IT_0100)*IT_0321 + conj(IT_0124)*IT_0323 + conj(IT_0073)*IT_0324 + conj
      (IT_0137)*IT_0325 + conj(IT_0147)*IT_0327 + conj(IT_0046)*IT_0328 + conj
      (IT_0106)*IT_0330 + conj(IT_0088)*IT_0332 + conj(IT_0111)*IT_0341;
    const ccomplex_t IT_0343 = IT_0123*conj(IT_0137) + IT_0004*IT_0225 + conj
      (IT_0126)*IT_0227 + conj(IT_0122)*IT_0241 + conj(IT_0111)*IT_0321 + conj
      (IT_0088)*IT_0323 + conj(IT_0046)*IT_0324 + conj(IT_0139)*IT_0325 + conj
      (IT_0106)*IT_0327 + conj(IT_0073)*IT_0328 + conj(IT_0147)*IT_0330 + conj
      (IT_0124)*IT_0332 + conj(IT_0100)*IT_0341;
    const ccomplex_t IT_0344 = s_13*s_23*IT_0005;
    const ccomplex_t IT_0345 = s_13*s_24*s_34;
    const ccomplex_t IT_0346 = (-2)*IT_0345;
    const ccomplex_t IT_0347 = IT_0344 + IT_0346;
    const ccomplex_t IT_0348 = IT_0000*IT_0347;
    const ccomplex_t IT_0349 = (-8)*IT_0348;
    const ccomplex_t IT_0350 = s_13*s_23;
    const ccomplex_t IT_0351 = 16*IT_0350;
    const ccomplex_t IT_0352 = IT_0349 + IT_0351;
    const ccomplex_t IT_0353 = 2*IT_0340;
    const ccomplex_t IT_0354 = 8*IT_0108;
    const ccomplex_t IT_0355 = s_12*IT_0107;
    const ccomplex_t IT_0356 = IT_0344 + IT_0355;
    const ccomplex_t IT_0357 = s_14*s_23*s_34;
    const ccomplex_t IT_0358 = -IT_0345 + -IT_0357;
    const ccomplex_t IT_0359 = IT_0356 + IT_0358;
    const ccomplex_t IT_0360 = IT_0000*IT_0359;
    const ccomplex_t IT_0361 = (-0.5)*IT_0360;
    const ccomplex_t IT_0362 = IT_0350 + IT_0361;
    const ccomplex_t IT_0363 = (-8)*IT_0362;
    const ccomplex_t IT_0364 = s_34*IT_0000;
    const ccomplex_t IT_0365 = s_14*s_23;
    const ccomplex_t IT_0366 = s_13*s_24;
    const ccomplex_t IT_0367 = s_12*s_34;
    const ccomplex_t IT_0368 = -IT_0366 + -IT_0367;
    const ccomplex_t IT_0369 = IT_0365 + IT_0368;
    const ccomplex_t IT_0370 = IT_0364*IT_0369;
    const ccomplex_t IT_0371 = (-4)*IT_0370;
    const ccomplex_t IT_0372 = (-2)*IT_0320;
    const ccomplex_t IT_0373 = conj(IT_0103)*IT_0119 + (-0.25)*IT_0109*conj
      (IT_0147) + (-0.25)*conj(IT_0126)*IT_0254 + (-0.25)*conj(IT_0122)*IT_0270 
      + (-0.25)*conj(IT_0046)*IT_0321 + (-0.25)*conj(IT_0073)*IT_0341 + (-0.25)
      *conj(IT_0111)*IT_0352 + (-0.25)*conj(IT_0139)*IT_0353 + (-0.25)*conj
      (IT_0124)*IT_0354 + (-0.25)*conj(IT_0106)*IT_0363 + (-0.25)*conj(IT_0088)
      *IT_0371 + (-0.25)*conj(IT_0137)*IT_0372;
    const ccomplex_t IT_0374 = conj(IT_0106)*IT_0109;
    const ccomplex_t IT_0375 = IT_0095*IT_0225 + conj(IT_0126)*IT_0258 + conj
      (IT_0122)*IT_0284 + conj(IT_0073)*IT_0321 + conj(IT_0046)*IT_0341 + conj
      (IT_0100)*IT_0352 + conj(IT_0137)*IT_0353 + conj(IT_0088)*IT_0354 + conj
      (IT_0147)*IT_0363 + conj(IT_0124)*IT_0371 + conj(IT_0139)*IT_0372 + IT_0374;
    const ccomplex_t IT_0376 = 4*IT_0331;
    const ccomplex_t IT_0377 = 6*IT_0329;
    const ccomplex_t IT_0378 = (-6)*IT_0326;
    const ccomplex_t IT_0379 = (-4)*IT_0322;
    const ccomplex_t IT_0380 = conj(IT_0046)*IT_0123 + conj(IT_0122)*IT_0233 +
       IT_0004*IT_0235 + conj(IT_0126)*IT_0237 + conj(IT_0137)*IT_0324 + conj
      (IT_0073)*IT_0325 + conj(IT_0139)*IT_0328 + conj(IT_0100)*IT_0353 + conj
      (IT_0111)*IT_0372 + conj(IT_0124)*IT_0376 + conj(IT_0147)*IT_0377 + conj
      (IT_0106)*IT_0378 + conj(IT_0088)*IT_0379;
    const ccomplex_t IT_0381 = conj(IT_0073)*IT_0123 + IT_0051*IT_0235 + conj
      (IT_0139)*IT_0324 + conj(IT_0046)*IT_0325 + conj(IT_0137)*IT_0328 + conj
      (IT_0111)*IT_0353 + conj(IT_0100)*IT_0372 + conj(IT_0088)*IT_0376 + conj
      (IT_0106)*IT_0377 + conj(IT_0147)*IT_0378 + conj(IT_0124)*IT_0379;
    const ccomplex_t IT_0382 = s_12*IT_0000*IT_0107;
    const ccomplex_t IT_0383 = 0.5*IT_0382;
    const ccomplex_t IT_0384 = s_13*s_24*s_34*IT_0000;
    const ccomplex_t IT_0385 = (-0.5)*IT_0384;
    const ccomplex_t IT_0386 = s_14*s_23*s_34*IT_0000;
    const ccomplex_t IT_0387 = (-0.5)*IT_0386;
    const ccomplex_t IT_0388 = IT_0350 + IT_0383 + IT_0385 + IT_0387;
    const ccomplex_t IT_0389 = 8*IT_0388;
    const ccomplex_t IT_0390 = 12*IT_0350;
    const ccomplex_t IT_0391 = IT_0151*IT_0235;
    const ccomplex_t IT_0392 = IT_0157*conj(IT_0162);
    const ccomplex_t IT_0393 = conj(IT_0162)*IT_0223;
    const ccomplex_t IT_0394 = IT_0109*(conj(IT_0088) + conj(IT_0111)) + conj
      (IT_0126)*IT_0217 + conj(IT_0122)*IT_0220 + conj(IT_0073)*IT_0327 + conj
      (IT_0046)*IT_0330 + conj(IT_0100)*IT_0363 + conj(IT_0137)*IT_0377 + conj
      (IT_0139)*IT_0378 + conj(IT_0124)*IT_0389 + conj(IT_0147)*IT_0390 +
       IT_0391 + IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = 12*IT_0108;
    const ccomplex_t IT_0396 = 0.25*conj(IT_0137);
    const ccomplex_t IT_0397 = 0.25*conj(IT_0139);
    const ccomplex_t IT_0398 = 0.25*conj(IT_0111);
    const ccomplex_t IT_0399 = 0.25*conj(IT_0046);
    const ccomplex_t IT_0400 = 0.25*conj(IT_0073);
    const ccomplex_t IT_0401 = 0.25*conj(IT_0126);
    const ccomplex_t IT_0402 = 0.25*conj(IT_0122);
    const ccomplex_t IT_0403 = 0.25*conj(IT_0147);
    const ccomplex_t IT_0404 = 0.25*conj(IT_0100);
    const ccomplex_t IT_0405 = 0.25*IT_0389;
    const ccomplex_t IT_0406 = IT_0089*conj(IT_0103)*IT_0125 + IT_0374 +
       IT_0391 + IT_0392 + IT_0393 + conj(IT_0088)*IT_0395 + 4*IT_0376*IT_0396 +
       4*IT_0379*IT_0397 + 4*IT_0354*IT_0398 + 4*IT_0332*IT_0399 + 4*IT_0323
      *IT_0400 + 4*IT_0249*IT_0401 + 4*IT_0310*IT_0402 + 4*IT_0389*IT_0403 + 4
      *IT_0371*IT_0404 + 4*conj(IT_0124)*(IT_0382 + IT_0405);
    const ccomplex_t IT_0407 = conj(IT_0103)*IT_0134*IT_0135 + IT_0077*IT_0225
       + conj(IT_0106)*IT_0389 + conj(IT_0124)*IT_0395 + 4*IT_0379*IT_0396 + 4
      *IT_0376*IT_0397 + 4*IT_0371*IT_0398 + 4*IT_0323*IT_0399 + 4*IT_0332
      *IT_0400 + 4*IT_0251*IT_0401 + 4*IT_0311*IT_0402 + 4*IT_0109*IT_0403 + 4
      *IT_0354*IT_0404 + 4*conj(IT_0088)*(IT_0382 + IT_0405);
    const ccomplex_t IT_0408 = IT_0155 + IT_0230;
    const ccomplex_t IT_0409 = IT_0137*IT_0408;
    const ccomplex_t IT_0410 = IT_0139*IT_0300;
    const ccomplex_t IT_0411 = IT_0153 + IT_0239;
    const ccomplex_t IT_0412 = IT_0046*IT_0411;
    const ccomplex_t IT_0413 = IT_0073*IT_0299;
    const ccomplex_t IT_0414 = IT_0106*IT_0222;
    const ccomplex_t IT_0415 = IT_0089*IT_0123*(IT_0122 + IT_0126) + IT_0100
      *IT_0271 + IT_0111*IT_0285 + IT_0103*IT_0134*IT_0286 + IT_0126*IT_0287 +
       IT_0122*IT_0288 + IT_0088*(IT_0247 + IT_0248 + IT_0290 + IT_0291) +
       IT_0124*(IT_0250 + IT_0294) + IT_0409 + IT_0410 + IT_0412 + IT_0413 +
       IT_0414;
    const ccomplex_t IT_0416 = IT_0104*IT_0163 + IT_0146*IT_0171 + IT_0124
      *IT_0195 + IT_0100*IT_0198 + IT_0111*IT_0212 + IT_0088*IT_0221 + IT_0147
      *IT_0224 + IT_0126*IT_0243 + IT_0122*IT_0245 + IT_0409 + IT_0410 + IT_0412
       + IT_0413 + IT_0414;
    const ccomplex_t IT_0417 = (IT_0004*conj(IT_0046) + IT_0051*conj(IT_0073) 
      + IT_0077*conj(IT_0088) + IT_0095*conj(IT_0100))*IT_0104 + conj(IT_0100)
      *IT_0110 + (-4)*IT_0103*(conj(IT_0111)*IT_0119 + (-0.25)*IT_0089*(conj
      (IT_0122)*IT_0123 + conj(IT_0124)*IT_0125 + conj(IT_0126)*IT_0127) + (
      -0.25)*IT_0134*(conj(IT_0088)*IT_0135 + conj(IT_0103)*IT_0136)) + (IT_0004
      *conj(IT_0137) + IT_0051*conj(IT_0139) + conj(IT_0106)*IT_0145)*IT_0146 +
       conj(IT_0147)*IT_0152 + conj(IT_0124)*(IT_0110 + IT_0152) + conj(IT_0126)
      *(IT_0154 + IT_0156 + IT_0158) + IT_0160*IT_0246 + IT_0126*IT_0263 +
       IT_0162*IT_0295 + conj(IT_0122)*(IT_0296 + IT_0297 + IT_0298) + conj
      (IT_0184)*IT_0303 + IT_0184*IT_0309 + IT_0122*IT_0312 + IT_0073*IT_0342 +
       IT_0046*IT_0343 + (-4)*IT_0111*IT_0373 + IT_0100*IT_0375 + IT_0137
      *IT_0380 + IT_0139*IT_0381 + IT_0106*(IT_0145*IT_0235 + conj(IT_0046)
      *IT_0327 + conj(IT_0073)*IT_0330 + conj(IT_0111)*IT_0363 + conj(IT_0139)
      *IT_0377 + conj(IT_0137)*IT_0378 + conj(IT_0088)*IT_0389 + conj(IT_0106)
      *IT_0390) + IT_0147*IT_0394 + IT_0124*IT_0406 + IT_0088*IT_0407 + conj
      (IT_0162)*IT_0415 + conj(IT_0160)*IT_0416;
    return create_ccomplex_return(IT_0417);
}

