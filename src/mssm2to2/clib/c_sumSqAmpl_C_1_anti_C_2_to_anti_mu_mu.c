#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_mu_mu.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_mu_mu(
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
    const creal_t m_H0 = param->m_H0;
    const creal_t m_mu = param->m_mu;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_snu_mu = param->m_snu_mu;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_numl = param->Gamma_numl;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_mu, 2);
    const ccomplex_t IT_0001 = m_C_1*m_C_2*IT_0000;
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = e_em*IT_0005;
    const ccomplex_t IT_0007 = cpow(IT_0004, -1);
    const ccomplex_t IT_0008 = IT_0002*IT_0007;
    const ccomplex_t IT_0009 = e_em*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + -IT_0009);
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0013 = IT_0005*IT_0012;
    const ccomplex_t IT_0014 = IT_0008*IT_0012;
    const ccomplex_t IT_0015 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + -IT_0014 + (-2)
      *IT_0016);
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = pow(m_C_1, 2);
    const ccomplex_t IT_0021 = pow(m_C_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0020 + -IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0027 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0028 = IT_0005*IT_0027;
    const ccomplex_t IT_0029 = IT_0008*IT_0027;
    const ccomplex_t IT_0030 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0031 = IT_0008*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + -IT_0029 + (-2)
      *IT_0031);
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = IT_0026*IT_0033;
    const ccomplex_t IT_0035 = IT_0022*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = 16*IT_0036;
    const ccomplex_t IT_0038 = s_12*IT_0000;
    const ccomplex_t IT_0039 = IT_0011*IT_0033;
    const ccomplex_t IT_0040 = IT_0022*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = 8*IT_0041;
    const ccomplex_t IT_0043 = s_13*s_24;
    const ccomplex_t IT_0044 = m_mu*m_C_1;
    const ccomplex_t IT_0045 = s_23 + s_24;
    const ccomplex_t IT_0046 = IT_0044*IT_0045;
    const ccomplex_t IT_0047 = pow(m_W, -1);
    const ccomplex_t IT_0048 = cos(beta);
    const ccomplex_t IT_0049 = cpow(IT_0048, -1);
    const ccomplex_t IT_0050 = sin(beta);
    const ccomplex_t IT_0051 = e_em*m_mu*IT_0007*IT_0047*IT_0049*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = V_u1*U_Wm2;
    const ccomplex_t IT_0054 = IT_0048*IT_0053;
    const ccomplex_t IT_0055 = U_d2*V_Wp1;
    const ccomplex_t IT_0056 = IT_0050*IT_0055;
    const ccomplex_t IT_0057 = IT_0007*(IT_0054 + IT_0056);
    const ccomplex_t IT_0058 = 1.4142135623731*e_em*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = IT_0052*IT_0059;
    const ccomplex_t IT_0061 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0062 = IT_0060*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = pow(m_Z, 2);
    const ccomplex_t IT_0065 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + IT_0064 
      + -reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0048*IT_0055;
    const ccomplex_t IT_0067 = IT_0050*IT_0053;
    const ccomplex_t IT_0068 = IT_0007*(IT_0066 + -IT_0067);
    const ccomplex_t IT_0069 = 1.4142135623731*e_em*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = e_em*m_mu*IT_0007*IT_0047;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0070*IT_0072;
    const ccomplex_t IT_0074 = IT_0065*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = -IT_0063 + -IT_0075;
    const ccomplex_t IT_0077 = (-4)*IT_0041;
    const ccomplex_t IT_0078 = s_34*m_C_1*m_C_2;
    const ccomplex_t IT_0079 = IT_0018*IT_0026;
    const ccomplex_t IT_0080 = IT_0022*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = 8*IT_0036;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0007;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*e_em*conj(V_Wp2)*IT_0007;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = cpow((-2)*s_13 + IT_0000 + IT_0020 + reg_prop +
       -m_snu_mu*(m_snu_mu + (0 + _Complex_I*-1)*Gamma_numl), -1);
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*U_d2*e_em*m_mu
      *IT_0007*IT_0047*IT_0049;
    const ccomplex_t IT_0093 = 0.5*IT_0092;
    const ccomplex_t IT_0094 = IT_0085*IT_0093;
    const ccomplex_t IT_0095 = IT_0088*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = s_12*s_34;
    const ccomplex_t IT_0099 = s_14*s_23;
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0038 + IT_0043 + IT_0098 + IT_0100;
    const ccomplex_t IT_0102 = 2*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_mu*IT_0007*IT_0047*IT_0049;
    const ccomplex_t IT_0104 = (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = IT_0086*IT_0104;
    const ccomplex_t IT_0106 = IT_0088*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = -IT_0107;
    const ccomplex_t IT_0109 = m_C_1*m_C_2;
    const ccomplex_t IT_0110 = s_34 + IT_0000;
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = 2*IT_0111;
    const ccomplex_t IT_0113 = m_mu*m_C_2;
    const ccomplex_t IT_0114 = s_13 + s_14;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = 4*IT_0115;
    const ccomplex_t IT_0117 = s_12*IT_0110;
    const ccomplex_t IT_0118 = 8*IT_0117;
    const ccomplex_t IT_0119 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0120 = IT_0048*IT_0119;
    const ccomplex_t IT_0121 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0122 = IT_0050*IT_0121;
    const ccomplex_t IT_0123 = IT_0007*(IT_0120 + IT_0122);
    const ccomplex_t IT_0124 = 1.4142135623731*e_em*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0052*IT_0125;
    const ccomplex_t IT_0127 = IT_0061*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0048*IT_0121;
    const ccomplex_t IT_0130 = IT_0050*IT_0119;
    const ccomplex_t IT_0131 = IT_0007*(IT_0129 + -IT_0130);
    const ccomplex_t IT_0132 = 1.4142135623731*e_em*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0072*IT_0133;
    const ccomplex_t IT_0135 = IT_0065*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = -IT_0128 + -IT_0136;
    const ccomplex_t IT_0138 = (-8)*IT_0111;
    const ccomplex_t IT_0139 = IT_0093*IT_0104;
    const ccomplex_t IT_0140 = IT_0088*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = 2*conj(IT_0036);
    const ccomplex_t IT_0143 = conj(IT_0141) + IT_0142;
    const ccomplex_t IT_0144 = (-4)*conj(IT_0041);
    const ccomplex_t IT_0145 = (-4)*IT_0115;
    const ccomplex_t IT_0146 = (-2)*IT_0115;
    const ccomplex_t IT_0147 = cos(alpha);
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*e_em*m_mu*IT_0007*IT_0047
      *IT_0049*IT_0147;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = IT_0121*IT_0147;
    const ccomplex_t IT_0151 = sin(alpha);
    const ccomplex_t IT_0152 = IT_0119*IT_0151;
    const ccomplex_t IT_0153 = IT_0007*(IT_0150 + IT_0152);
    const ccomplex_t IT_0154 = 1.4142135623731*e_em*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = 0.5*IT_0155;
    const ccomplex_t IT_0157 = IT_0149*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*e_em*m_mu*IT_0007*IT_0047
      *IT_0049*IT_0151;
    const ccomplex_t IT_0162 = 0.5*IT_0161;
    const ccomplex_t IT_0163 = IT_0119*IT_0147;
    const ccomplex_t IT_0164 = IT_0121*IT_0151;
    const ccomplex_t IT_0165 = IT_0007*(IT_0163 + -IT_0164);
    const ccomplex_t IT_0166 = 1.4142135623731*e_em*IT_0165;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0162*IT_0168;
    const ccomplex_t IT_0170 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0020 + -IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0171 = IT_0169*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = -IT_0160 + -IT_0172;
    const ccomplex_t IT_0174 = -IT_0000;
    const ccomplex_t IT_0175 = s_34 + IT_0174;
    const ccomplex_t IT_0176 = IT_0109*IT_0175;
    const ccomplex_t IT_0177 = (-8)*IT_0176;
    const ccomplex_t IT_0178 = IT_0055*IT_0147;
    const ccomplex_t IT_0179 = IT_0053*IT_0151;
    const ccomplex_t IT_0180 = IT_0007*(IT_0178 + IT_0179);
    const ccomplex_t IT_0181 = 1.4142135623731*e_em*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = IT_0149*IT_0183;
    const ccomplex_t IT_0185 = IT_0158*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = IT_0053*IT_0147;
    const ccomplex_t IT_0188 = IT_0055*IT_0151;
    const ccomplex_t IT_0189 = IT_0007*(IT_0187 + -IT_0188);
    const ccomplex_t IT_0190 = 1.4142135623731*e_em*IT_0189;
    const ccomplex_t IT_0191 = 0.5*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0162*IT_0192;
    const ccomplex_t IT_0194 = IT_0170*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = -IT_0186 + -IT_0195;
    const ccomplex_t IT_0197 = -s_14;
    const ccomplex_t IT_0198 = s_13 + IT_0197;
    const ccomplex_t IT_0199 = IT_0113*IT_0198;
    const ccomplex_t IT_0200 = (-4)*IT_0199;
    const ccomplex_t IT_0201 = conj(IT_0036) + conj(IT_0041);
    const ccomplex_t IT_0202 = -s_24;
    const ccomplex_t IT_0203 = s_23 + IT_0202;
    const ccomplex_t IT_0204 = IT_0044*IT_0203;
    const ccomplex_t IT_0205 = 4*IT_0204;
    const ccomplex_t IT_0206 = conj(IT_0025) + conj(IT_0082);
    const ccomplex_t IT_0207 = -IT_0098;
    const ccomplex_t IT_0208 = -IT_0043;
    const ccomplex_t IT_0209 = IT_0038 + IT_0099 + IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = 2*IT_0209;
    const ccomplex_t IT_0211 = -s_34;
    const ccomplex_t IT_0212 = IT_0000 + IT_0211;
    const ccomplex_t IT_0213 = s_12*IT_0212;
    const ccomplex_t IT_0214 = (-8)*IT_0213;
    const ccomplex_t IT_0215 = 2*IT_0199;
    const ccomplex_t IT_0216 = -conj(IT_0141);
    const ccomplex_t IT_0217 = (-2)*IT_0176;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = 2*IT_0204;
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = 2*IT_0046;
    const ccomplex_t IT_0222 = (-8)*IT_0043;
    const ccomplex_t IT_0223 = (-4)*IT_0038;
    const ccomplex_t IT_0224 = m_mu*s_24*m_C_1;
    const ccomplex_t IT_0225 = (-4)*IT_0224;
    const ccomplex_t IT_0226 = m_mu*s_13*m_C_2;
    const ccomplex_t IT_0227 = (-4)*IT_0226;
    const ccomplex_t IT_0228 = (-2)*IT_0199;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = IT_0043*conj(IT_0091) + -IT_0001*IT_0143 + 0.25
      *IT_0078*IT_0144 + 0.25*conj(IT_0076)*IT_0146 + 0.25*conj(IT_0173)*IT_0220
       + 0.25*conj(IT_0137)*IT_0221 + 0.25*conj(IT_0025)*IT_0222 + 0.25*conj
      (IT_0082)*IT_0223 + 0.25*conj(IT_0108)*IT_0225 + 0.25*conj(IT_0097)
      *IT_0227 + 0.25*conj(IT_0196)*IT_0229;
    const ccomplex_t IT_0231 = (-2)*IT_0204;
    const ccomplex_t IT_0232 = 8*conj(IT_0141);
    const ccomplex_t IT_0233 = 8*conj(IT_0041);
    const ccomplex_t IT_0234 = 8*IT_0038;
    const ccomplex_t IT_0235 = 8*IT_0224;
    const ccomplex_t IT_0236 = 8*IT_0226;
    const ccomplex_t IT_0237 = (-4)*IT_0046;
    const ccomplex_t IT_0238 = conj(IT_0173)*IT_0205;
    const ccomplex_t IT_0239 = conj(IT_0196)*IT_0200;
    const ccomplex_t IT_0240 = 16*conj(IT_0025)*IT_0043 + conj(IT_0076)
      *IT_0116 + conj(IT_0091)*IT_0222 + IT_0001*(16*conj(IT_0036) + IT_0232) +
       IT_0078*IT_0233 + conj(IT_0082)*IT_0234 + conj(IT_0108)*IT_0235 + conj
      (IT_0097)*IT_0236 + conj(IT_0137)*IT_0237 + IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = 4*conj(IT_0141);
    const ccomplex_t IT_0242 = 8*conj(IT_0036);
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = 4*IT_0046;
    const ccomplex_t IT_0245 = m_mu*s_14*m_C_2;
    const ccomplex_t IT_0246 = 4*IT_0245;
    const ccomplex_t IT_0247 = m_mu*s_23*m_C_1;
    const ccomplex_t IT_0248 = 4*IT_0247;
    const ccomplex_t IT_0249 = 16*IT_0001*conj(IT_0041) + 16*conj(IT_0082)
      *IT_0099 + conj(IT_0076)*IT_0145 + conj(IT_0091)*IT_0223 + conj(IT_0025)
      *IT_0234 + IT_0238 + IT_0239 + IT_0078*IT_0243 + conj(IT_0137)*IT_0244 +
       conj(IT_0097)*IT_0246 + conj(IT_0108)*IT_0248;
    const ccomplex_t IT_0250 = (-2)*IT_0101;
    const ccomplex_t IT_0251 = (-2)*IT_0111;
    const ccomplex_t IT_0252 = IT_0137*(IT_0118*conj(IT_0137) + conj(IT_0076)
      *IT_0138 + 4*IT_0115*(conj(IT_0041) + -1./2*conj(IT_0141) + -1./2*IT_0142)
       + conj(IT_0091)*IT_0221 + conj(IT_0025)*IT_0237 + conj(IT_0082)*IT_0244 +
       conj(IT_0108)*IT_0250 + conj(IT_0097)*IT_0251);
    const ccomplex_t IT_0253 = 0.25*conj(IT_0196);
    const ccomplex_t IT_0254 = 0.25*conj(IT_0173);
    const ccomplex_t IT_0255 = 0.25*conj(IT_0076);
    const ccomplex_t IT_0256 = 0.25*conj(IT_0137);
    const ccomplex_t IT_0257 = 0.25*conj(IT_0091);
    const ccomplex_t IT_0258 = 0.25*conj(IT_0025);
    const ccomplex_t IT_0259 = 0.25*conj(IT_0082);
    const ccomplex_t IT_0260 = IT_0043*conj(IT_0097) + IT_0001*conj(IT_0108) +
       0.25*IT_0224*IT_0243 + conj(IT_0041)*IT_0247 + IT_0210*IT_0253 + IT_0218
      *IT_0254 + IT_0102*IT_0255 + IT_0251*IT_0256 + IT_0227*IT_0257 + IT_0236
      *IT_0258 + IT_0246*IT_0259;
    const ccomplex_t IT_0261 = IT_0001*conj(IT_0097) + IT_0043*conj(IT_0108) +
       0.25*IT_0226*IT_0243 + conj(IT_0041)*IT_0245 + IT_0218*IT_0253 + IT_0210
      *IT_0254 + IT_0112*IT_0255 + IT_0250*IT_0256 + IT_0225*IT_0257 + IT_0235
      *IT_0258 + IT_0248*IT_0259;
    const ccomplex_t IT_0262 = conj(IT_0097)*IT_0224;
    const ccomplex_t IT_0263 = conj(IT_0108)*IT_0226;
    const ccomplex_t IT_0264 = IT_0038*conj(IT_0041) + 0.5*IT_0046*conj
      (IT_0076) + IT_0078*conj(IT_0082) + 2*IT_0001*(conj(IT_0025) + -1./2*conj
      (IT_0091)) + (-0.5)*IT_0115*conj(IT_0137) + IT_0043*IT_0143 + (-0.25)*conj
      (IT_0173)*IT_0215 + (-0.25)*conj(IT_0196)*IT_0231 + IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = conj(IT_0196)*IT_0205;
    const ccomplex_t IT_0266 = conj(IT_0173)*IT_0200;
    const ccomplex_t IT_0267 = IT_0001*conj(IT_0025)*IT_0037 + conj(IT_0036)*
      (IT_0038*IT_0042 + IT_0037*IT_0043) + IT_0046*conj(IT_0076)*IT_0077 +
       IT_0078*(conj(IT_0025)*IT_0042 + conj(IT_0082)*IT_0083 + IT_0077*conj
      (IT_0091)) + IT_0076*(conj(IT_0097)*IT_0102 + conj(IT_0108)*IT_0112 + conj
      (IT_0025)*IT_0116 + conj(IT_0076)*IT_0118 + conj(IT_0137)*IT_0138 + 2
      *IT_0046*(IT_0143 + 0.5*IT_0144) + conj(IT_0082)*IT_0145 + conj(IT_0091)
      *IT_0146) + IT_0173*(IT_0177*conj(IT_0196) + IT_0200*IT_0201 + IT_0205
      *IT_0206 + conj(IT_0108)*IT_0210 + conj(IT_0173)*IT_0214 + IT_0215*IT_0216
       + conj(IT_0097)*IT_0218 + conj(IT_0091)*IT_0220) + 4*IT_0091*IT_0230 +
       IT_0196*(conj(IT_0173)*IT_0177 + IT_0201*IT_0205 + IT_0200*IT_0206 + conj
      (IT_0097)*IT_0210 + conj(IT_0196)*IT_0214 + conj(IT_0108)*IT_0218 + conj
      (IT_0091)*IT_0229 + IT_0216*IT_0231) + IT_0025*IT_0240 + IT_0082*IT_0249 +
       IT_0252 + 4*IT_0097*IT_0260 + 4*IT_0108*IT_0261 + IT_0083*(IT_0262 +
       IT_0263) + 4*IT_0141*IT_0264 + IT_0041*(16*IT_0001*conj(IT_0082) + 16
      *conj(IT_0041)*IT_0099 + 4*IT_0115*conj(IT_0137) + IT_0038*IT_0241 + 4
      *conj(IT_0108)*IT_0245 + 4*conj(IT_0097)*IT_0247 + IT_0265 + IT_0266) + 4
      *IT_0036*(IT_0046*conj(IT_0076) + (-2)*IT_0001*conj(IT_0091) + -IT_0115
      *conj(IT_0137) + 0.25*IT_0043*IT_0232 + 0.25*IT_0038*IT_0233 + 0.25
      *IT_0265 + 0.25*IT_0266);
    return create_ccomplex_return(IT_0267);
}

