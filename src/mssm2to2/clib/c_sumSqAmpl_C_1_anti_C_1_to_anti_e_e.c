#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_1_to_anti_e_e.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_1_to_anti_e_e(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_e = param->m_e;
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
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_nuel = param->Gamma_nuel;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_W, -1);
    const ccomplex_t IT_0001 = cos(beta);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = sin(alpha);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0003*IT_0005;
    const ccomplex_t IT_0007 = 0.5*IT_0006;
    const ccomplex_t IT_0008 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0009 = cos(alpha);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0012 = IT_0003*IT_0011;
    const ccomplex_t IT_0013 = IT_0005*(IT_0010 + -IT_0012);
    const ccomplex_t IT_0014 = 1.4142135623731*e_em*IT_0013;
    const ccomplex_t IT_0015 = 0.5*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*IT_0015;
    const ccomplex_t IT_0017 = IT_0007*IT_0016;
    const ccomplex_t IT_0018 = pow(m_C_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*m_e*e_em*IT_0000*IT_0002
      *IT_0005*IT_0009;
    const ccomplex_t IT_0023 = (-0.5)*IT_0022;
    const ccomplex_t IT_0024 = IT_0009*IT_0011;
    const ccomplex_t IT_0025 = IT_0003*IT_0008;
    const ccomplex_t IT_0026 = IT_0005*(IT_0024 + IT_0025);
    const ccomplex_t IT_0027 = 1.4142135623731*e_em*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = IT_0023*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + (-2)*IT_0018 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = -IT_0021 + -IT_0033;
    const ccomplex_t IT_0035 = -s_34;
    const ccomplex_t IT_0036 = pow(m_e, 2);
    const ccomplex_t IT_0037 = IT_0035 + IT_0036;
    const ccomplex_t IT_0038 = s_12*IT_0037;
    const ccomplex_t IT_0039 = V_u1*U_Wm1;
    const ccomplex_t IT_0040 = IT_0009*IT_0039;
    const ccomplex_t IT_0041 = U_d1*V_Wp1;
    const ccomplex_t IT_0042 = IT_0003*IT_0041;
    const ccomplex_t IT_0043 = IT_0005*(IT_0040 + -IT_0042);
    const ccomplex_t IT_0044 = 1.4142135623731*e_em*IT_0043;
    const ccomplex_t IT_0045 = 0.5*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = IT_0007*IT_0046;
    const ccomplex_t IT_0048 = IT_0019*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0009*IT_0041;
    const ccomplex_t IT_0051 = IT_0003*IT_0039;
    const ccomplex_t IT_0052 = IT_0005*(IT_0050 + IT_0051);
    const ccomplex_t IT_0053 = 1.4142135623731*e_em*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = IT_0023*IT_0055;
    const ccomplex_t IT_0057 = IT_0031*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0049 + -IT_0058;
    const ccomplex_t IT_0060 = IT_0018*IT_0037;
    const ccomplex_t IT_0061 = sin(beta);
    const ccomplex_t IT_0062 = m_e*e_em*IT_0000*IT_0002*IT_0005*IT_0061;
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = IT_0001*IT_0008;
    const ccomplex_t IT_0065 = IT_0011*IT_0061;
    const ccomplex_t IT_0066 = IT_0005*(IT_0064 + IT_0065);
    const ccomplex_t IT_0067 = 1.4142135623731*e_em*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = cpow((-2)*s_12 + (-2)*IT_0018 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = IT_0001*IT_0011;
    const ccomplex_t IT_0074 = IT_0008*IT_0061;
    const ccomplex_t IT_0075 = IT_0005*(IT_0073 + -IT_0074);
    const ccomplex_t IT_0076 = 1.4142135623731*e_em*IT_0075;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = m_e*e_em*IT_0000*IT_0005;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = IT_0077*IT_0079;
    const ccomplex_t IT_0081 = pow(m_Z, 2);
    const ccomplex_t IT_0082 = cpow((-2)*s_12 + (-2)*IT_0018 + IT_0081 + 
      -reg_prop, -1);
    const ccomplex_t IT_0083 = IT_0080*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*IT_0083;
    const ccomplex_t IT_0085 = -IT_0072 + -IT_0084;
    const ccomplex_t IT_0086 = s_34 + IT_0036;
    const ccomplex_t IT_0087 = s_12*IT_0086;
    const ccomplex_t IT_0088 = IT_0001*IT_0039;
    const ccomplex_t IT_0089 = IT_0041*IT_0061;
    const ccomplex_t IT_0090 = IT_0005*(IT_0088 + IT_0089);
    const ccomplex_t IT_0091 = 1.4142135623731*e_em*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0063*IT_0092;
    const ccomplex_t IT_0094 = IT_0070*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = IT_0001*IT_0041;
    const ccomplex_t IT_0097 = IT_0039*IT_0061;
    const ccomplex_t IT_0098 = IT_0005*(IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = 1.4142135623731*e_em*IT_0098;
    const ccomplex_t IT_0100 = (-0.5)*IT_0099;
    const ccomplex_t IT_0101 = IT_0079*IT_0100;
    const ccomplex_t IT_0102 = IT_0082*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = -IT_0095 + -IT_0103;
    const ccomplex_t IT_0105 = IT_0018*IT_0086;
    const ccomplex_t IT_0106 = cos(theta_W);
    const ccomplex_t IT_0107 = cpow(IT_0106, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0107;
    const ccomplex_t IT_0109 = IT_0004*IT_0107;
    const ccomplex_t IT_0110 = U_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = IT_0005*IT_0106;
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = e_em*U_Wm1*conj(U_Wm1);
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0111 + -IT_0113 + (-2)
      *IT_0115);
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = IT_0108*IT_0117;
    const ccomplex_t IT_0119 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0018 + -reg_prop, -1);
    const ccomplex_t IT_0120 = IT_0118*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = m_e*m_C_1;
    const ccomplex_t IT_0123 = -s_24;
    const ccomplex_t IT_0124 = s_23 + IT_0123;
    const ccomplex_t IT_0125 = IT_0122*IT_0124;
    const ccomplex_t IT_0126 = 4*IT_0125;
    const ccomplex_t IT_0127 = -s_14;
    const ccomplex_t IT_0128 = s_13 + IT_0127;
    const ccomplex_t IT_0129 = IT_0122*IT_0128;
    const ccomplex_t IT_0130 = (-4)*IT_0129;
    const ccomplex_t IT_0131 = IT_0034*IT_0130;
    const ccomplex_t IT_0132 = IT_0059*IT_0126 + IT_0131;
    const ccomplex_t IT_0133 = s_23 + s_24;
    const ccomplex_t IT_0134 = IT_0122*IT_0133;
    const ccomplex_t IT_0135 = 4*IT_0104;
    const ccomplex_t IT_0136 = IT_0018*IT_0036;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0005;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*e_em*conj(V_Wp1)*IT_0005;
    const ccomplex_t IT_0140 = IT_0138*IT_0139;
    const ccomplex_t IT_0141 = cpow((-2)*s_13 + IT_0018 + IT_0036 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_nuel), -1);
    const ccomplex_t IT_0142 = IT_0140*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*IT_0142;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = (-8)*IT_0144;
    const ccomplex_t IT_0146 = s_13 + s_14;
    const ccomplex_t IT_0147 = IT_0122*IT_0146;
    const ccomplex_t IT_0148 = (-4)*IT_0085;
    const ccomplex_t IT_0149 = e_em*IT_0109;
    const ccomplex_t IT_0150 = e_em*IT_0112;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0149 + -IT_0150);
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = V_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0154 = IT_0109*IT_0153;
    const ccomplex_t IT_0155 = IT_0112*IT_0153;
    const ccomplex_t IT_0156 = e_em*V_Wp1*conj(V_Wp1);
    const ccomplex_t IT_0157 = IT_0112*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + -IT_0155 + (-2)
      *IT_0157);
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0152*IT_0159;
    const ccomplex_t IT_0161 = IT_0119*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = IT_0059*IT_0130;
    const ccomplex_t IT_0165 = IT_0034*IT_0126 + IT_0164;
    const ccomplex_t IT_0166 = s_13*s_24;
    const ccomplex_t IT_0167 = IT_0108*IT_0159;
    const ccomplex_t IT_0168 = IT_0119*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = m_e*s_14*m_C_1;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1.4142135623731)*m_e*U_d1*e_em
      *IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0173 = 0.5*IT_0172;
    const ccomplex_t IT_0174 = IT_0138*IT_0173;
    const ccomplex_t IT_0175 = IT_0141*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = 4*IT_0177;
    const ccomplex_t IT_0179 = 4*IT_0085;
    const ccomplex_t IT_0180 = m_e*s_23*m_C_1;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1.4142135623731)*m_e*conj(U_d1)
      *e_em*IT_0000*IT_0002*IT_0005;
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0139*IT_0182;
    const ccomplex_t IT_0184 = IT_0141*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = 4*IT_0186;
    const ccomplex_t IT_0188 = s_14*s_23;
    const ccomplex_t IT_0189 = 16*IT_0170;
    const ccomplex_t IT_0190 = s_12*IT_0036;
    const ccomplex_t IT_0191 = (-4)*IT_0144;
    const ccomplex_t IT_0192 = (-4)*IT_0104;
    const ccomplex_t IT_0193 = IT_0117*IT_0152;
    const ccomplex_t IT_0194 = IT_0119*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = IT_0173*IT_0182;
    const ccomplex_t IT_0197 = IT_0141*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = 4*IT_0198;
    const ccomplex_t IT_0200 = s_34*IT_0018;
    const ccomplex_t IT_0201 = 2*IT_0104;
    const ccomplex_t IT_0202 = m_e*s_24*m_C_1;
    const ccomplex_t IT_0203 = m_e*s_13*m_C_1;
    const ccomplex_t IT_0204 = (-2)*IT_0085;
    const ccomplex_t IT_0205 = (-2)*IT_0125;
    const ccomplex_t IT_0206 = -IT_0059;
    const ccomplex_t IT_0207 = 2*IT_0129;
    const ccomplex_t IT_0208 = -IT_0034;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = IT_0209*IT_0210;
    const ccomplex_t IT_0212 = cpow(s_12 + IT_0018 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0213 = IT_0211*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = 0.5*s_14 + 0.5*s_23;
    const ccomplex_t IT_0217 = s_13 + s_24;
    const ccomplex_t IT_0218 = IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = IT_0122*IT_0218;
    const ccomplex_t IT_0220 = 2*IT_0136;
    const ccomplex_t IT_0221 = 2*IT_0188;
    const ccomplex_t IT_0222 = IT_0190 + IT_0200 + IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = -s_14 + -s_23;
    const ccomplex_t IT_0224 = IT_0217 + IT_0223;
    const ccomplex_t IT_0225 = IT_0122*IT_0224;
    const ccomplex_t IT_0226 = (-8)*IT_0225;
    const ccomplex_t IT_0227 = 2*IT_0166;
    const ccomplex_t IT_0228 = IT_0190 + IT_0200 + IT_0220 + IT_0227;
    const ccomplex_t IT_0229 = 4*IT_0144;
    const ccomplex_t IT_0230 = (-4)*conj(IT_0144);
    const ccomplex_t IT_0231 = 0.125*IT_0230;
    const ccomplex_t IT_0232 = (64*IT_0136 + 32*IT_0166 + 32*IT_0188 + 32
      *IT_0190 + 32*IT_0200)*IT_0215 + 8*(IT_0177 + -conj(IT_0177) + IT_0186 + 
      -conj(IT_0186))*IT_0219 + (-8)*(conj(IT_0170) + conj(IT_0195))*IT_0222 + -
      (conj(IT_0034) + conj(IT_0059) + IT_0206 + IT_0208)*IT_0226 + (-8)*IT_0228
      *(conj(IT_0121) + conj(IT_0163) + (-0.5)*IT_0198 + 0.5*conj(IT_0198) +
       0.125*IT_0229 + IT_0231);
    const ccomplex_t IT_0233 = -conj(IT_0085);
    const ccomplex_t IT_0234 = conj(IT_0163) + IT_0231;
    const ccomplex_t IT_0235 = 2*conj(IT_0104);
    const ccomplex_t IT_0236 = 4*conj(IT_0177);
    const ccomplex_t IT_0237 = 4*conj(IT_0186);
    const ccomplex_t IT_0238 = -conj(IT_0059);
    const ccomplex_t IT_0239 = -conj(IT_0034);
    const ccomplex_t IT_0240 = (-2)*conj(IT_0085);
    const ccomplex_t IT_0241 = 0.125*IT_0147;
    const ccomplex_t IT_0242 = conj(IT_0121)*IT_0166 + 0.5*conj(IT_0170)
      *IT_0200 + IT_0136*IT_0234 + 0.125*IT_0134*IT_0235 + 0.125*IT_0202*IT_0236
       + 0.125*IT_0203*IT_0237 + 0.125*IT_0205*IT_0238 + 0.125*IT_0207*IT_0239 +
       IT_0240*IT_0241;
    const ccomplex_t IT_0243 = 2*conj(IT_0121);
    const ccomplex_t IT_0244 = conj(IT_0198) + IT_0243;
    const ccomplex_t IT_0245 = (-8)*IT_0215;
    const ccomplex_t IT_0246 = 0.125*conj(IT_0059)*IT_0130;
    const ccomplex_t IT_0247 = 0.125*conj(IT_0034);
    const ccomplex_t IT_0248 = (-4)*conj(IT_0104);
    const ccomplex_t IT_0249 = 4*conj(IT_0085);
    const ccomplex_t IT_0250 = 0.125*IT_0134;
    const ccomplex_t IT_0251 = IT_0190*IT_0234 + 0.125*IT_0171*IT_0236 + 0.125
      *IT_0180*IT_0237 + 0.5*IT_0200*IT_0244 + 0.125*IT_0222*IT_0245 + IT_0246 +
       IT_0126*IT_0247 + IT_0241*IT_0248 + IT_0249*IT_0250;
    const ccomplex_t IT_0252 = (-8)*conj(IT_0144);
    const ccomplex_t IT_0253 = 0.0625*IT_0252;
    const ccomplex_t IT_0254 = conj(IT_0163) + IT_0253;
    const ccomplex_t IT_0255 = (-4)*conj(IT_0085);
    const ccomplex_t IT_0256 = 4*conj(IT_0104);
    const ccomplex_t IT_0257 = 0.125*IT_0228*IT_0245;
    const ccomplex_t IT_0258 = 0.125*conj(IT_0059)*IT_0126 + IT_0190*conj
      (IT_0195) + conj(IT_0170)*IT_0200 + conj(IT_0177)*IT_0202 + conj(IT_0186)
      *IT_0203 + IT_0166*IT_0244 + IT_0130*IT_0247 + 2*IT_0136*IT_0254 + IT_0241
      *IT_0255 + IT_0250*IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = conj(IT_0170)*IT_0190 + conj(IT_0195)*IT_0200 +
       conj(IT_0186)*IT_0202 + conj(IT_0177)*IT_0203 + IT_0136*IT_0244 + IT_0246
       + IT_0126*IT_0247 + 2*IT_0166*IT_0254 + IT_0250*IT_0255 + IT_0241*IT_0256
       + IT_0257;
    const ccomplex_t IT_0260 = 0.25*conj(IT_0034);
    const ccomplex_t IT_0261 = 0.0625*conj(IT_0059)*IT_0126 + IT_0136*conj
      (IT_0170) + IT_0188*conj(IT_0195) + 0.5*IT_0200*IT_0234 + 0.0625*IT_0180
      *IT_0236 + 0.0625*IT_0171*IT_0237 + 0.25*IT_0190*IT_0244 + 0.0625*IT_0222
      *IT_0245 + 0.0625*IT_0134*IT_0248 + 0.0625*IT_0147*IT_0249 + 0.25*IT_0130
      *IT_0260;
    const ccomplex_t IT_0262 = s_12*s_34;
    const ccomplex_t IT_0263 = -IT_0188;
    const ccomplex_t IT_0264 = IT_0166 + IT_0190 + IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = -IT_0262;
    const ccomplex_t IT_0266 = -IT_0166;
    const ccomplex_t IT_0267 = IT_0188 + IT_0190 + IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = (-0.5)*IT_0166;
    const ccomplex_t IT_0269 = (-2)*conj(IT_0186)*(IT_0059*IT_0060 + (-0.5)
      *IT_0136*IT_0178 + (-0.5)*IT_0105*IT_0201 + (-0.5)*IT_0191*IT_0202 + (-0.5
      )*IT_0204*IT_0264 + -IT_0034*IT_0267 + IT_0187*IT_0268);
    const ccomplex_t IT_0270 = (-2)*IT_0129;
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = 2*IT_0125;
    const ccomplex_t IT_0273 = (-0.5)*IT_0272;
    const ccomplex_t IT_0274 = (-0.25)*IT_0264;
    const ccomplex_t IT_0275 = (-8)*IT_0059*(IT_0038*conj(IT_0059) + -conj
      (IT_0034)*IT_0060) + (-8)*IT_0034*(conj(IT_0034)*IT_0038 + -conj(IT_0059)
      *IT_0060) + 8*IT_0085*(conj(IT_0085)*IT_0087 + -conj(IT_0104)*IT_0105) +
       conj(IT_0121)*(IT_0132 + IT_0134*IT_0135 + IT_0136*IT_0145 + IT_0147
      *IT_0148) + conj(IT_0163)*(IT_0135*IT_0147 + IT_0134*IT_0148 + IT_0165 +
       IT_0145*IT_0166) + conj(IT_0170)*(IT_0165 + IT_0171*IT_0178 + IT_0134
      *IT_0179 + IT_0180*IT_0187 + IT_0188*IT_0189 + IT_0190*IT_0191 + IT_0147
      *IT_0192) + conj(IT_0195)*(IT_0132 + IT_0147*IT_0179 + IT_0178*IT_0180 +
       IT_0171*IT_0187 + IT_0136*IT_0189 + IT_0134*IT_0192 + IT_0190*IT_0199 +
       IT_0191*IT_0200) + conj(IT_0198)*(IT_0136*IT_0191 + IT_0166*IT_0199 +
       IT_0134*IT_0201 + IT_0178*IT_0202 + IT_0187*IT_0203 + IT_0147*IT_0204 +
       IT_0205*IT_0206 + IT_0207*IT_0208) + -IT_0215*IT_0232 + 8*IT_0104*
      (IT_0087*conj(IT_0104) + IT_0105*IT_0233) + 8*IT_0198*IT_0242 + 8*IT_0170
      *IT_0251 + 8*IT_0121*IT_0258 + 8*IT_0163*IT_0259 + 16*IT_0195*IT_0261 + 8
      *IT_0186*((-0.25)*conj(IT_0059)*IT_0060 + conj(IT_0121)*IT_0203 + IT_0202
      *IT_0234 + 0.125*IT_0105*IT_0235 + 0.125*IT_0136*IT_0236 + 0.125*IT_0240
      *IT_0264 + IT_0260*IT_0267) + IT_0269 + 2*conj(IT_0144)*(IT_0085*IT_0134 +
       -IT_0104*IT_0147 + -IT_0229*IT_0268 + -IT_0206*IT_0271 + -IT_0208*IT_0273
      ) + (-2)*IT_0144*(conj(IT_0104)*IT_0147 + IT_0134*IT_0233 + IT_0238
      *IT_0271 + IT_0239*IT_0273) + (-4)*conj(IT_0177)*(0.5*IT_0034*IT_0060 +
       0.5*IT_0085*IT_0105 + (-0.25)*IT_0166*IT_0178 + IT_0144*IT_0203 + (-0.5)
      *IT_0059*IT_0267 + IT_0201*IT_0274) + (-4)*IT_0177*(0.5*conj(IT_0034)
      *IT_0060 + 0.5*conj(IT_0085)*IT_0105 + (-2)*conj(IT_0121)*IT_0202 + (conj
      (IT_0144) + (-2)*conj(IT_0163))*IT_0203 + (-0.5)*conj(IT_0059)*IT_0267 +
       IT_0235*IT_0274);
    return create_ccomplex_return(IT_0275);
}

