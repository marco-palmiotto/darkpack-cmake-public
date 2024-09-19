#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_A_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_A_W(
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0001 = 2*IT_0000;
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0007 = IT_0003*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = (0 + _Complex_I*1)*(IT_0005 + 0.5*IT_0008);
    const ccomplex_t IT_0010 = -IT_0009;
    const ccomplex_t IT_0011 = pow(m_C_2, 2);
    const ccomplex_t IT_0012 = pow(m_N_3, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0011 + -IT_0012 + -reg_prop, -1);
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*IT_0014;
    const ccomplex_t IT_0016 = IT_0001*IT_0015;
    const ccomplex_t IT_0017 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0018 = IT_0003*IT_0017;
    const ccomplex_t IT_0019 = conj(N_u3)*V_u2*e_em;
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
    const ccomplex_t IT_0047 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = 2*IT_0048;
    const ccomplex_t IT_0050 = cos(theta_W);
    const ccomplex_t IT_0051 = cpow(IT_0050, -1);
    const ccomplex_t IT_0052 = IT_0045*IT_0051;
    const ccomplex_t IT_0053 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_W3*V_u2*e_em;
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
    const ccomplex_t IT_0067 = m_N_3*IT_0000;
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
    const ccomplex_t IT_0081 = IT_0025*IT_0040;
    const ccomplex_t IT_0082 = cos(beta);
    const ccomplex_t IT_0083 = IT_0003*IT_0082;
    const ccomplex_t IT_0084 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = (-2)*IT_0085;
    const ccomplex_t IT_0087 = IT_0051*IT_0082;
    const ccomplex_t IT_0088 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0092 = IT_0083*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = IT_0086 + IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = (-0.5)*IT_0096;
    const ccomplex_t IT_0098 = IT_0063*IT_0097;
    const ccomplex_t IT_0099 = IT_0044*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = -IT_0081 + IT_0100;
    const ccomplex_t IT_0102 = (-2)*IT_0101;
    const ccomplex_t IT_0103 = IT_0080 + IT_0102;
    const ccomplex_t IT_0104 = 2*conj(IT_0079);
    const ccomplex_t IT_0105 = (-2)*conj(IT_0101);
    const ccomplex_t IT_0106 = IT_0104 + IT_0105;
    const ccomplex_t IT_0107 = s_14*s_24*IT_0028;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = s_12 + IT_0108;
    const ccomplex_t IT_0110 = m_N_3*IT_0109;
    const ccomplex_t IT_0111 = m_C_2*IT_0109;
    const ccomplex_t IT_0112 = IT_0000*IT_0010;
    const ccomplex_t IT_0113 = cpow(s_13 + (0 + _Complex_I*(-0.5))*m_C_2
      *Gamma_c2 + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0112*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = 2*IT_0069;
    const ccomplex_t IT_0117 = (-2)*IT_0041;
    const ccomplex_t IT_0118 = IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = 2*conj(IT_0069);
    const ccomplex_t IT_0120 = (-2)*conj(IT_0041);
    const ccomplex_t IT_0121 = IT_0119 + IT_0120;
    const ccomplex_t IT_0122 = IT_0000*IT_0023;
    const ccomplex_t IT_0123 = IT_0113*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = -IT_0076 + IT_0124;
    const ccomplex_t IT_0126 = 2*IT_0101;
    const ccomplex_t IT_0127 = (-2)*IT_0079;
    const ccomplex_t IT_0128 = IT_0126 + IT_0127;
    const ccomplex_t IT_0129 = 2*conj(IT_0101);
    const ccomplex_t IT_0130 = (-2)*conj(IT_0079);
    const ccomplex_t IT_0131 = IT_0129 + IT_0130;
    const ccomplex_t IT_0132 = pow(s_14, 2);
    const ccomplex_t IT_0133 = IT_0028*IT_0132;
    const ccomplex_t IT_0134 = -IT_0011;
    const ccomplex_t IT_0135 = IT_0133 + IT_0134;
    const ccomplex_t IT_0136 = m_N_3*IT_0135;
    const ccomplex_t IT_0137 = s_13*m_N_3;
    const ccomplex_t IT_0138 = -IT_0101;
    const ccomplex_t IT_0139 = IT_0079 + IT_0138;
    const ccomplex_t IT_0140 = -conj(IT_0101);
    const ccomplex_t IT_0141 = conj(IT_0079) + IT_0140;
    const ccomplex_t IT_0142 = 6*IT_0041;
    const ccomplex_t IT_0143 = (-6)*IT_0069;
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = 6*conj(IT_0041);
    const ccomplex_t IT_0146 = (-6)*conj(IT_0069);
    const ccomplex_t IT_0147 = IT_0145 + IT_0146;
    const ccomplex_t IT_0148 = s_23*m_C_2;
    const ccomplex_t IT_0149 = 6*IT_0069;
    const ccomplex_t IT_0150 = (-6)*IT_0041;
    const ccomplex_t IT_0151 = IT_0149 + IT_0150;
    const ccomplex_t IT_0152 = 6*conj(IT_0069);
    const ccomplex_t IT_0153 = (-6)*conj(IT_0041);
    const ccomplex_t IT_0154 = IT_0152 + IT_0153;
    const ccomplex_t IT_0155 = IT_0011*IT_0109;
    const ccomplex_t IT_0156 = m_C_2*m_N_3;
    const ccomplex_t IT_0157 = IT_0135*IT_0156;
    const ccomplex_t IT_0158 = IT_0032*IT_0156;
    const ccomplex_t IT_0159 = 2*IT_0038;
    const ccomplex_t IT_0160 = 2*conj(IT_0038);
    const ccomplex_t IT_0161 = s_12*IT_0135;
    const ccomplex_t IT_0162 = s_12*IT_0032;
    const ccomplex_t IT_0163 = IT_0000*IT_0025;
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = 4*IT_0164*conj(IT_0164);
    const ccomplex_t IT_0166 = 4*IT_0038*conj(IT_0038);
    const ccomplex_t IT_0167 = 6*s_12;
    const ccomplex_t IT_0168 = pow(s_34, 2);
    const ccomplex_t IT_0169 = s_12*IT_0028*IT_0168;
    const ccomplex_t IT_0170 = 0.5*IT_0169;
    const ccomplex_t IT_0171 = s_13*s_24*s_34*IT_0028;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = s_14*s_23*s_34*IT_0028;
    const ccomplex_t IT_0174 = (-0.5)*IT_0173;
    const ccomplex_t IT_0175 = IT_0027 + IT_0170 + IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = -IT_0000;
    const ccomplex_t IT_0177 = IT_0015*IT_0176;
    const ccomplex_t IT_0178 = 8*IT_0177;
    const ccomplex_t IT_0179 = 8*conj(IT_0177);
    const ccomplex_t IT_0180 = IT_0025*IT_0176;
    const ccomplex_t IT_0181 = 8*IT_0180;
    const ccomplex_t IT_0182 = 8*conj(IT_0180);
    const ccomplex_t IT_0183 = m_C_2*m_N_3*IT_0028*IT_0168;
    const ccomplex_t IT_0184 = (-0.5)*IT_0124;
    const ccomplex_t IT_0185 = 4*IT_0016;
    const ccomplex_t IT_0186 = 4*conj(IT_0016);
    const ccomplex_t IT_0187 = (-0.5)*IT_0115;
    const ccomplex_t IT_0188 = 4*conj(IT_0177)*(IT_0026 + 3*IT_0180) + 4
      *IT_0177*(conj(IT_0026) + 3*conj(IT_0180)) + IT_0179*IT_0184 + IT_0178
      *conj(IT_0184) + (conj(IT_0180) + conj(IT_0184))*IT_0185 + (IT_0180 +
       IT_0184)*IT_0186 + 4*(conj(IT_0026) + 0.25*IT_0182)*IT_0187 + 4*(IT_0026 
      + 0.25*IT_0181)*conj(IT_0187);
    const ccomplex_t IT_0189 = 8*conj(IT_0038)*IT_0164;
    const ccomplex_t IT_0190 = (-6)*s_12;
    const ccomplex_t IT_0191 = IT_0036 + IT_0115;
    const ccomplex_t IT_0192 = conj(IT_0036) + conj(IT_0115);
    const ccomplex_t IT_0193 = (-6)*s_13;
    const ccomplex_t IT_0194 = (-6)*s_23;
    const ccomplex_t IT_0195 = (-6)*IT_0011;
    const ccomplex_t IT_0196 = (-6)*IT_0012;
    const ccomplex_t IT_0197 = conj(IT_0101)*IT_0144 + IT_0101*IT_0147 + conj
      (IT_0079)*IT_0151 + IT_0079*IT_0154 + IT_0190*(conj(IT_0036)*IT_0125 +
       IT_0036*conj(IT_0125) + conj(IT_0077)*IT_0191 + IT_0077*IT_0192) + (conj
      (IT_0125)*IT_0177 + IT_0125*conj(IT_0177) + conj(IT_0180)*IT_0191 +
       IT_0180*IT_0192)*IT_0193 + (conj(IT_0077)*IT_0177 + IT_0077*conj(IT_0177)
       + conj(IT_0036)*IT_0180 + IT_0036*conj(IT_0180))*IT_0194 + (conj(IT_0125)
      *IT_0191 + IT_0125*IT_0192)*IT_0195 + (conj(IT_0036)*IT_0077 + IT_0036
      *conj(IT_0077))*IT_0196;
    const ccomplex_t IT_0198 = s_12*IT_0043;
    const ccomplex_t IT_0199 = s_14*s_24;
    const ccomplex_t IT_0200 = (-2)*IT_0199;
    const ccomplex_t IT_0201 = IT_0198 + IT_0200;
    const ccomplex_t IT_0202 = IT_0028*IT_0201;
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = s_12 + IT_0203;
    const ccomplex_t IT_0205 = (-4)*s_12;
    const ccomplex_t IT_0206 = (-4)*s_13;
    const ccomplex_t IT_0207 = (-4)*s_23;
    const ccomplex_t IT_0208 = (-4)*IT_0011;
    const ccomplex_t IT_0209 = (-4)*IT_0012;
    const ccomplex_t IT_0210 = (conj(IT_0077)*IT_0125 + IT_0077*conj(IT_0125) 
      + conj(IT_0036)*IT_0191 + IT_0036*IT_0192)*IT_0205 + (conj(IT_0125)
      *IT_0180 + IT_0125*conj(IT_0180) + conj(IT_0177)*IT_0191 + IT_0177*IT_0192
      )*IT_0206 + (conj(IT_0036)*IT_0177 + IT_0036*conj(IT_0177) + conj(IT_0077)
      *IT_0180 + IT_0077*conj(IT_0180))*IT_0207 + (IT_0125*conj(IT_0125) +
       IT_0191*IT_0192)*IT_0208 + (IT_0036*conj(IT_0036) + IT_0077*conj(IT_0077)
      )*IT_0209;
    const ccomplex_t IT_0211 = s_23*IT_0011;
    const ccomplex_t IT_0212 = s_23*IT_0011*IT_0043;
    const ccomplex_t IT_0213 = s_24*s_34*IT_0011;
    const ccomplex_t IT_0214 = (-2)*IT_0213;
    const ccomplex_t IT_0215 = IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = IT_0028*IT_0215;
    const ccomplex_t IT_0217 = (-0.5)*IT_0216;
    const ccomplex_t IT_0218 = IT_0211 + IT_0217;
    const ccomplex_t IT_0219 = s_13*IT_0012;
    const ccomplex_t IT_0220 = s_12*s_23*IT_0043;
    const ccomplex_t IT_0221 = s_13*IT_0029;
    const ccomplex_t IT_0222 = IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = s_14*s_23*s_24;
    const ccomplex_t IT_0224 = s_12*s_24*s_34;
    const ccomplex_t IT_0225 = s_13*IT_0012*IT_0043;
    const ccomplex_t IT_0226 = -IT_0223 + -IT_0224 + (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = IT_0222 + IT_0226;
    const ccomplex_t IT_0228 = IT_0028*IT_0227;
    const ccomplex_t IT_0229 = s_12*s_23;
    const ccomplex_t IT_0230 = (-2)*IT_0229;
    const ccomplex_t IT_0231 = IT_0219 + IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = s_13*s_23*IT_0043;
    const ccomplex_t IT_0233 = s_12*IT_0168;
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = s_14*s_23*s_34;
    const ccomplex_t IT_0236 = s_13*s_24*s_34;
    const ccomplex_t IT_0237 = -IT_0235 + -IT_0236;
    const ccomplex_t IT_0238 = IT_0234 + IT_0237;
    const ccomplex_t IT_0239 = IT_0028*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = IT_0027 + IT_0240;
    const ccomplex_t IT_0242 = (-2)*IT_0236;
    const ccomplex_t IT_0243 = IT_0232 + IT_0242;
    const ccomplex_t IT_0244 = IT_0028*IT_0243;
    const ccomplex_t IT_0245 = (-8)*IT_0244;
    const ccomplex_t IT_0246 = 16*IT_0027;
    const ccomplex_t IT_0247 = IT_0245 + IT_0246;
    const ccomplex_t IT_0248 = s_13*IT_0043;
    const ccomplex_t IT_0249 = IT_0028*IT_0248;
    const ccomplex_t IT_0250 = (-0.25)*IT_0249;
    const ccomplex_t IT_0251 = s_13 + IT_0250;
    const ccomplex_t IT_0252 = IT_0156*IT_0251;
    const ccomplex_t IT_0253 = (-4)*IT_0252;
    const ccomplex_t IT_0254 = 2*IT_0036;
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = IT_0115 + IT_0255;
    const ccomplex_t IT_0257 = s_14*s_34;
    const ccomplex_t IT_0258 = (-2)*IT_0257;
    const ccomplex_t IT_0259 = IT_0248 + IT_0258;
    const ccomplex_t IT_0260 = IT_0028*IT_0259;
    const ccomplex_t IT_0261 = (-0.5)*IT_0260;
    const ccomplex_t IT_0262 = s_13 + IT_0261;
    const ccomplex_t IT_0263 = m_N_3*IT_0262;
    const ccomplex_t IT_0264 = (-2)*IT_0263;
    const ccomplex_t IT_0265 = s_23*IT_0043;
    const ccomplex_t IT_0266 = s_24*s_34;
    const ccomplex_t IT_0267 = (-2)*IT_0266;
    const ccomplex_t IT_0268 = IT_0265 + IT_0267;
    const ccomplex_t IT_0269 = IT_0028*IT_0268;
    const ccomplex_t IT_0270 = (-0.5)*IT_0269;
    const ccomplex_t IT_0271 = s_23 + IT_0270;
    const ccomplex_t IT_0272 = m_C_2*IT_0271;
    const ccomplex_t IT_0273 = (-2)*IT_0272;
    const ccomplex_t IT_0274 = IT_0028*IT_0265;
    const ccomplex_t IT_0275 = (-0.25)*IT_0274;
    const ccomplex_t IT_0276 = s_23 + IT_0275;
    const ccomplex_t IT_0277 = IT_0156*IT_0276;
    const ccomplex_t IT_0278 = (-4)*IT_0277;
    const ccomplex_t IT_0279 = s_34*IT_0028;
    const ccomplex_t IT_0280 = s_14*s_23;
    const ccomplex_t IT_0281 = s_13*s_24;
    const ccomplex_t IT_0282 = s_12*s_34;
    const ccomplex_t IT_0283 = -IT_0281 + -IT_0282;
    const ccomplex_t IT_0284 = IT_0280 + IT_0283;
    const ccomplex_t IT_0285 = IT_0279*IT_0284;
    const ccomplex_t IT_0286 = (-4)*IT_0180;
    const ccomplex_t IT_0287 = IT_0125*IT_0218 + -IT_0077*IT_0231 + 2*IT_0026
      *IT_0241 + (-0.25)*IT_0184*IT_0247 + (-0.5)*IT_0253*IT_0256 + (-0.25)
      *IT_0071*IT_0264 + (-0.25)*IT_0128*IT_0273 + (-0.25)*IT_0254*IT_0278 + (
      -0.25)*IT_0285*IT_0286;
    const ccomplex_t IT_0288 = s_14*s_34*IT_0028;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = s_13 + IT_0289;
    const ccomplex_t IT_0291 = s_12*IT_0290;
    const ccomplex_t IT_0292 = s_24*s_34*IT_0028;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = s_23 + IT_0293;
    const ccomplex_t IT_0295 = s_12*IT_0294;
    const ccomplex_t IT_0296 = (-4)*IT_0177;
    const ccomplex_t IT_0297 = conj(IT_0038)*IT_0296;
    const ccomplex_t IT_0298 = (-4)*conj(IT_0177);
    const ccomplex_t IT_0299 = IT_0038*IT_0298;
    const ccomplex_t IT_0300 = conj(IT_0164)*IT_0286;
    const ccomplex_t IT_0301 = (-4)*conj(IT_0180);
    const ccomplex_t IT_0302 = IT_0164*IT_0301;
    const ccomplex_t IT_0303 = s_23*IT_0028*IT_0132;
    const ccomplex_t IT_0304 = -IT_0211;
    const ccomplex_t IT_0305 = s_12*s_14*s_34*IT_0028;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = s_12*s_13;
    const ccomplex_t IT_0308 = 2*IT_0307;
    const ccomplex_t IT_0309 = s_13*s_14*s_24*IT_0028;
    const ccomplex_t IT_0310 = -IT_0309;
    const ccomplex_t IT_0311 = IT_0303 + IT_0304 + IT_0306 + IT_0308 + IT_0310;
    const ccomplex_t IT_0312 = s_14*s_23*s_24*IT_0028;
    const ccomplex_t IT_0313 = s_12*s_24*s_34*IT_0028;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = s_13*IT_0028*IT_0029;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = IT_0219 + IT_0312 + IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = (-4)*IT_0164*conj(IT_0184);
    const ccomplex_t IT_0319 = (-4)*conj(IT_0164)*IT_0184;
    const ccomplex_t IT_0320 = (-4)*IT_0038*conj(IT_0187);
    const ccomplex_t IT_0321 = (-4)*conj(IT_0038)*IT_0187;
    const ccomplex_t IT_0322 = IT_0012*IT_0109;
    const ccomplex_t IT_0323 = s_12*IT_0109;
    const ccomplex_t IT_0324 = (-8)*IT_0038;
    const ccomplex_t IT_0325 = IT_0011*IT_0032;
    const ccomplex_t IT_0326 = (-2)*IT_0036;
    const ccomplex_t IT_0327 = IT_0012*IT_0135;
    const ccomplex_t IT_0328 = IT_0109*IT_0156;
    const ccomplex_t IT_0329 = IT_0156*IT_0290;
    const ccomplex_t IT_0330 = IT_0156*IT_0294;
    const ccomplex_t IT_0331 = 2*IT_0077;
    const ccomplex_t IT_0332 = 2*IT_0125;
    const ccomplex_t IT_0333 = (-8)*IT_0164;
    const ccomplex_t IT_0334 = (-2)*IT_0077;
    const ccomplex_t IT_0335 = IT_0303 + IT_0304 + IT_0306 + IT_0309;
    const ccomplex_t IT_0336 = -IT_0219;
    const ccomplex_t IT_0337 = IT_0312 + IT_0314 + IT_0315 + IT_0336;
    const ccomplex_t IT_0338 = IT_0335 + IT_0337;
    const ccomplex_t IT_0339 = (-2)*IT_0180;
    const ccomplex_t IT_0340 = (-2)*IT_0125;
    const ccomplex_t IT_0341 = 0.125*IT_0033*IT_0103 + 0.125*IT_0071*IT_0110 +
       0.125*IT_0111*IT_0128 + 0.125*IT_0118*IT_0136 + 0.125*IT_0158*IT_0254 +
       0.25*IT_0157*IT_0256 + 0.25*IT_0326*IT_0328 + IT_0038*(IT_0157 + IT_0158 
      + (-2)*IT_0328) + 0.125*IT_0177*((-10)*IT_0329 + (-10)*IT_0330) + 0.125
      *IT_0322*IT_0331 + 0.125*IT_0155*IT_0332 + 0.125*IT_0323*IT_0333 + 0.125
      *IT_0327*IT_0334 + 0.125*IT_0338*IT_0339 + 0.125*IT_0325*IT_0340;
    const ccomplex_t IT_0342 = 2*conj(IT_0036);
    const ccomplex_t IT_0343 = 0.5*IT_0342;
    const ccomplex_t IT_0344 = conj(IT_0115) + IT_0343;
    const ccomplex_t IT_0345 = 2*conj(IT_0077);
    const ccomplex_t IT_0346 = 2*conj(IT_0125);
    const ccomplex_t IT_0347 = (-2)*conj(IT_0077);
    const ccomplex_t IT_0348 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0349 = (-2)*conj(IT_0125);
    const ccomplex_t IT_0350 = (-2)*conj(IT_0036);
    const ccomplex_t IT_0351 = 0.5*IT_0033*IT_0106 + 0.5*IT_0075*IT_0110 + 0.5
      *IT_0111*IT_0131 + 0.5*IT_0121*IT_0136 + 0.5*IT_0158*IT_0342 + IT_0157
      *IT_0344 + 0.5*IT_0322*IT_0345 + 0.5*IT_0155*IT_0346 + 0.5*IT_0327*IT_0347
       + 0.5*IT_0338*IT_0348 + 0.5*IT_0325*IT_0349 + 0.5*IT_0328*((-16)*conj
      (IT_0038) + 2*IT_0350);
    const ccomplex_t IT_0352 = (-2)*IT_0164;
    const ccomplex_t IT_0353 = (-2)*conj(IT_0164);
    const ccomplex_t IT_0354 = (-2)*IT_0038;
    const ccomplex_t IT_0355 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0356 = (-2)*IT_0177;
    const ccomplex_t IT_0357 = conj(IT_0038)*IT_0356;
    const ccomplex_t IT_0358 = (-2)*conj(IT_0177);
    const ccomplex_t IT_0359 = IT_0038*IT_0358;
    const ccomplex_t IT_0360 = conj(IT_0016)*IT_0354;
    const ccomplex_t IT_0361 = IT_0016*IT_0355;
    const ccomplex_t IT_0362 = (-2)*IT_0016;
    const ccomplex_t IT_0363 = (-2)*conj(IT_0016);
    const ccomplex_t IT_0364 = IT_0356 + IT_0362;
    const ccomplex_t IT_0365 = IT_0358 + IT_0363;
    const ccomplex_t IT_0366 = (-10)*IT_0164*conj(IT_0177);
    const ccomplex_t IT_0367 = (-10)*IT_0038*conj(IT_0180);
    const ccomplex_t IT_0368 = (-10)*conj(IT_0038)*IT_0180;
    const ccomplex_t IT_0369 = (-8)*conj(IT_0164)*IT_0187;
    const ccomplex_t IT_0370 = conj(IT_0184)*IT_0324;
    const ccomplex_t IT_0371 = (-8)*conj(IT_0038)*IT_0184;
    const ccomplex_t IT_0372 = conj(IT_0026)*IT_0354;
    const ccomplex_t IT_0373 = IT_0026*IT_0355;
    const ccomplex_t IT_0374 = conj(IT_0164)*IT_0362;
    const ccomplex_t IT_0375 = IT_0164*IT_0363;
    const ccomplex_t IT_0376 = conj(IT_0180)*IT_0326 + (-2)*conj(IT_0026)*
      (IT_0115 + (-0.5)*IT_0326) + conj(IT_0187)*IT_0333 + conj(IT_0115)*IT_0339
       + IT_0115*IT_0348 + IT_0180*IT_0350 + (-2)*IT_0026*(conj(IT_0115) + (-0.5
      )*IT_0350) + conj(IT_0125)*IT_0364 + IT_0125*IT_0365 + IT_0366 + IT_0367 +
       IT_0368 + IT_0369 + IT_0370 + IT_0371 + IT_0372 + IT_0373 + IT_0374 +
       IT_0375;
    const ccomplex_t IT_0377 = (conj(IT_0026) + conj(IT_0180))*IT_0326 + conj
      (IT_0187)*IT_0333 + (IT_0026 + IT_0180)*IT_0350 + conj(IT_0077)*IT_0364 +
       IT_0077*IT_0365 + IT_0366 + IT_0367 + IT_0368 + IT_0369 + IT_0370 +
       IT_0371 + IT_0372 + IT_0373 + IT_0374 + IT_0375;
    const ccomplex_t IT_0378 = m_C_2*IT_0294;
    const ccomplex_t IT_0379 = m_N_3*IT_0290;
    const ccomplex_t IT_0380 = IT_0011*IT_0294;
    const ccomplex_t IT_0381 = IT_0012*IT_0290;
    const ccomplex_t IT_0382 = 0.25*IT_0175;
    const ccomplex_t IT_0383 = (-0.25)*IT_0381;
    const ccomplex_t IT_0384 = (-0.25)*IT_0264;
    const ccomplex_t IT_0385 = (-0.25)*IT_0273;
    const ccomplex_t IT_0386 = (-0.25)*IT_0253;
    const ccomplex_t IT_0387 = (-0.25)*IT_0278;
    const ccomplex_t IT_0388 = (-0.25)*IT_0285;
    const ccomplex_t IT_0389 = -IT_0191*IT_0218 + IT_0036*IT_0231 + (-2)
      *IT_0016*IT_0241 + -IT_0103*IT_0384 + -IT_0118*IT_0385 + -IT_0332*IT_0386 
      + -IT_0331*IT_0387 + -IT_0296*IT_0388;
    const ccomplex_t IT_0390 = -IT_0192*IT_0218 + conj(IT_0036)*IT_0231 + (-2)
      *conj(IT_0016)*IT_0241 + 0.25*conj(IT_0187)*IT_0247 + -IT_0106*IT_0384 + 
      -IT_0121*IT_0385 + -IT_0346*IT_0386 + -IT_0345*IT_0387 + -IT_0298*IT_0388;
    const ccomplex_t IT_0391 = (-4)*IT_0184*(conj(IT_0125)*IT_0218 + -conj
      (IT_0077)*IT_0231 + 2*conj(IT_0026)*IT_0241 + (-0.5)*IT_0253*IT_0344 +
       IT_0075*IT_0384 + IT_0131*IT_0385 + IT_0342*IT_0387 + IT_0301*IT_0388);
    const ccomplex_t IT_0392 = (12*IT_0016*conj(IT_0016) + 12*IT_0026*conj
      (IT_0026))*IT_0027 + IT_0033*(IT_0039*IT_0071 + IT_0072*IT_0075 + conj
      (IT_0077)*IT_0103 + IT_0077*IT_0106) + (IT_0075*IT_0077 + IT_0071*conj
      (IT_0077) + IT_0039*IT_0103 + IT_0072*IT_0106)*IT_0110 + IT_0111*((conj
      (IT_0036) + conj(IT_0038) + conj(IT_0115))*IT_0118 + (IT_0036 + IT_0038 +
       IT_0115)*IT_0121 + conj(IT_0125)*IT_0128 + IT_0125*IT_0131) + (IT_0121
      *IT_0125 + IT_0118*conj(IT_0125) + (conj(IT_0036) + conj(IT_0038) + conj
      (IT_0115))*IT_0128 + (IT_0036 + IT_0038 + IT_0115)*IT_0131)*IT_0136 + 6
      *IT_0137*(conj(IT_0016)*IT_0139 + IT_0016*IT_0141 + 0.166666666666667*conj
      (IT_0026)*IT_0144 + 0.166666666666667*IT_0026*IT_0147) + (-6)*IT_0148*
      (conj(IT_0026)*IT_0139 + IT_0026*IT_0141 + (-0.166666666666667)*conj
      (IT_0016)*IT_0151 + (-0.166666666666667)*IT_0016*IT_0154) + (conj(IT_0115)
      *IT_0155 + conj(IT_0125)*IT_0157 + conj(IT_0077)*IT_0158)*IT_0159 + 
      (IT_0115*IT_0155 + IT_0125*IT_0157 + IT_0077*IT_0158)*IT_0160 + (IT_0161 +
       IT_0162)*(IT_0165 + IT_0166) + (IT_0041*conj(IT_0041) + IT_0069*conj
      (IT_0069) + IT_0079*conj(IT_0079) + IT_0101*conj(IT_0101))*IT_0167 +
       IT_0175*(conj(IT_0016)*IT_0178 + IT_0016*IT_0179 + conj(IT_0026)*IT_0181 
      + IT_0026*IT_0182) + IT_0183*IT_0188 + (IT_0157 + IT_0158)*IT_0189 + (conj
      (IT_0041)*IT_0069 + IT_0041*conj(IT_0069) + conj(IT_0079)*IT_0101 +
       IT_0079*conj(IT_0101))*IT_0190 + IT_0156*IT_0197 + IT_0204*IT_0210 + (-4)
      *conj(IT_0184)*IT_0287 + (IT_0291 + IT_0295)*(IT_0297 + IT_0299 + IT_0300 
      + IT_0302) + (IT_0311 + IT_0317)*(IT_0318 + IT_0319 + IT_0320 + IT_0321) +
       conj(IT_0038)*(IT_0155*IT_0254 + IT_0254*IT_0322 + IT_0323*IT_0324 +
       IT_0325*IT_0326 + IT_0326*IT_0327) + 8*conj(IT_0164)*IT_0341 + IT_0038*
      (IT_0155*IT_0342 + IT_0322*IT_0342 + IT_0325*IT_0350 + IT_0327*IT_0350) +
       2*IT_0164*IT_0351 + IT_0325*(conj(IT_0115)*IT_0354 + IT_0115*IT_0355) +
       IT_0328*(conj(IT_0115)*IT_0352 + IT_0115*IT_0353 + (conj(IT_0077) + conj
      (IT_0125))*IT_0354 + (IT_0077 + IT_0125)*IT_0355) + IT_0337*(conj(IT_0016)
      *IT_0326 + IT_0016*IT_0350 + conj(IT_0026)*(IT_0334 + IT_0352) + IT_0026*
      (IT_0347 + IT_0353) + IT_0357 + IT_0359 + IT_0360 + IT_0361) + IT_0335*
      (conj(IT_0016)*IT_0326 + IT_0016*IT_0350 + conj(IT_0026)*(IT_0340 +
       IT_0352) + IT_0026*(IT_0349 + IT_0353) + IT_0357 + IT_0359 + IT_0360 +
       IT_0361 + conj(IT_0115)*IT_0362 + IT_0115*IT_0363) + IT_0329*IT_0376 +
       IT_0330*IT_0377 + 4*conj(IT_0177)*(IT_0169*IT_0177 + -(IT_0041 + -IT_0069
      )*IT_0378 + IT_0139*IT_0379 + 0.5*IT_0256*IT_0380 + 0.25*IT_0254*IT_0381 +
       IT_0178*IT_0382) + 4*conj(IT_0180)*(IT_0169*IT_0180 + -IT_0139*IT_0378 + 
      (IT_0041 + -IT_0069)*IT_0379 + 0.25*IT_0332*IT_0380 + 0.25*IT_0331*IT_0381
       + IT_0181*IT_0382) + (-4)*IT_0177*((conj(IT_0041) + -conj(IT_0069))
      *IT_0378 + -IT_0141*IT_0379 + (-0.5)*IT_0344*IT_0380 + IT_0342*IT_0383) + 
      (-4)*IT_0180*(IT_0141*IT_0378 + -(conj(IT_0041) + -conj(IT_0069))*IT_0379 
      + (-0.25)*IT_0346*IT_0380 + IT_0345*IT_0383) + 4*conj(IT_0187)*IT_0389 + 4
      *IT_0187*IT_0390 + IT_0391;
    return create_ccomplex_return(IT_0392);
}

