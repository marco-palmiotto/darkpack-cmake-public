#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_anti_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_N_2, 2);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = -IT_0000;
    const ccomplex_t IT_0005 = IT_0003 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = sin(beta);
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = sin(alpha);
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*e_em*(IT_0007*IT_0008 + 
      -IT_0009*IT_0010)*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0007*IT_0015;
    const ccomplex_t IT_0017 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = IT_0010*IT_0015;
    const ccomplex_t IT_0020 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0007*IT_0012;
    const ccomplex_t IT_0023 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = IT_0010*IT_0012;
    const ccomplex_t IT_0026 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0018 + -IT_0021 + 
      -IT_0024 + IT_0027);
    const ccomplex_t IT_0029 = IT_0013*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = IT_0016*IT_0020;
    const ccomplex_t IT_0036 = IT_0017*IT_0019;
    const ccomplex_t IT_0037 = IT_0022*IT_0026;
    const ccomplex_t IT_0038 = IT_0023*IT_0025;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + IT_0036 + 
      -IT_0037 + -IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*e_em*(IT_0007*IT_0009 +
       IT_0008*IT_0010)*IT_0012;
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0040*IT_0042;
    const ccomplex_t IT_0044 = IT_0034*IT_0043;
    const ccomplex_t IT_0045 = IT_0009*IT_0012;
    const ccomplex_t IT_0046 = conj(N_u2)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = 2*IT_0047;
    const ccomplex_t IT_0049 = IT_0009*IT_0015;
    const ccomplex_t IT_0050 = conj(N_B2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = conj(N_W2)*conj(V_u1)*e_em;
    const ccomplex_t IT_0054 = IT_0045*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = IT_0052 + IT_0055;
    const ccomplex_t IT_0057 = IT_0048 + IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0061 = IT_0012*IT_0060;
    const ccomplex_t IT_0062 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0063 = IT_0012*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-0.5)*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0059*IT_0066;
    const ccomplex_t IT_0068 = pow(m_Hp, 2);
    const ccomplex_t IT_0069 = cpow((-2)*s_23 + IT_0000 + IT_0068 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0067*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0000 + IT_0068 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0073 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0074 = IT_0012*IT_0073;
    const ccomplex_t IT_0075 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0076 = IT_0012*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + (-0.5)*IT_0077);
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = conj(N_u2)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0081 = IT_0045*IT_0080;
    const ccomplex_t IT_0082 = 2*IT_0081;
    const ccomplex_t IT_0083 = conj(N_B2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0084 = IT_0049*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W2)*conj(V_u2)*e_em;
    const ccomplex_t IT_0087 = IT_0045*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = IT_0085 + IT_0088;
    const ccomplex_t IT_0090 = IT_0082 + IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = 0.5*IT_0091;
    const ccomplex_t IT_0093 = IT_0079*IT_0092;
    const ccomplex_t IT_0094 = IT_0072*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = IT_0020*IT_0049;
    const ccomplex_t IT_0099 = IT_0008*IT_0015;
    const ccomplex_t IT_0100 = IT_0017*IT_0099;
    const ccomplex_t IT_0101 = IT_0026*IT_0045;
    const ccomplex_t IT_0102 = IT_0008*IT_0012;
    const ccomplex_t IT_0103 = IT_0023*IT_0102;
    const ccomplex_t IT_0104 = IT_0098 + -IT_0100 + -IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = e_em*IT_0012;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = IT_0097*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0032 + (-0.5)*IT_0044 + 2*IT_0071 + 2
      *IT_0095 + (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0111 = IT_0012*IT_0110;
    const ccomplex_t IT_0112 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0113 = IT_0012*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + 0.5*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = N_d2*e_em*U_Wm1;
    const ccomplex_t IT_0118 = IT_0102*IT_0117;
    const ccomplex_t IT_0119 = (-2)*IT_0118;
    const ccomplex_t IT_0120 = N_B2*U_d1*e_em;
    const ccomplex_t IT_0121 = IT_0099*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = N_W2*U_d1*e_em;
    const ccomplex_t IT_0124 = IT_0102*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = IT_0122 + IT_0125;
    const ccomplex_t IT_0127 = IT_0119 + IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0128;
    const ccomplex_t IT_0130 = IT_0116*IT_0129;
    const ccomplex_t IT_0131 = IT_0069*IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0134 = IT_0012*IT_0133;
    const ccomplex_t IT_0135 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0136 = IT_0012*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + 0.5*IT_0137);
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = N_d2*e_em*U_Wm2;
    const ccomplex_t IT_0141 = IT_0102*IT_0140;
    const ccomplex_t IT_0142 = (-2)*IT_0141;
    const ccomplex_t IT_0143 = N_B2*U_d2*e_em;
    const ccomplex_t IT_0144 = IT_0099*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = N_W2*U_d2*e_em;
    const ccomplex_t IT_0147 = IT_0102*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = IT_0145 + IT_0148;
    const ccomplex_t IT_0150 = IT_0142 + IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = IT_0139*IT_0152;
    const ccomplex_t IT_0154 = IT_0072*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0157 = IT_0049*IT_0156;
    const ccomplex_t IT_0158 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0159 = IT_0099*IT_0158;
    const ccomplex_t IT_0160 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0161 = IT_0045*IT_0160;
    const ccomplex_t IT_0162 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0163 = IT_0102*IT_0162;
    const ccomplex_t IT_0164 = -IT_0157 + IT_0159 + IT_0161 + -IT_0163;
    const ccomplex_t IT_0165 = IT_0106*IT_0164;
    const ccomplex_t IT_0166 = IT_0097*IT_0165;
    const ccomplex_t IT_0167 = IT_0016*IT_0158;
    const ccomplex_t IT_0168 = IT_0019*IT_0156;
    const ccomplex_t IT_0169 = IT_0022*IT_0162;
    const ccomplex_t IT_0170 = IT_0025*IT_0160;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*(IT_0167 + -IT_0168 + 
      -IT_0169 + IT_0170);
    const ccomplex_t IT_0172 = IT_0013*IT_0171;
    const ccomplex_t IT_0173 = IT_0031*IT_0172;
    const ccomplex_t IT_0174 = IT_0016*IT_0156;
    const ccomplex_t IT_0175 = IT_0019*IT_0158;
    const ccomplex_t IT_0176 = IT_0022*IT_0160;
    const ccomplex_t IT_0177 = IT_0025*IT_0162;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + IT_0175 + 
      -IT_0176 + -IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0042*IT_0179;
    const ccomplex_t IT_0181 = IT_0034*IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0132 + 2*IT_0155 + 0.5*IT_0166 + 0.5
      *IT_0173 + 0.5*IT_0181;
    const ccomplex_t IT_0183 = 2*IT_0182;
    const ccomplex_t IT_0184 = (-0.5)*IT_0166 + (-0.5)*IT_0173 + (-0.5)*IT_0181;
    const ccomplex_t IT_0185 = (-2)*IT_0184;
    const ccomplex_t IT_0186 = IT_0183 + IT_0185;
    const ccomplex_t IT_0187 = 2*conj(IT_0182);
    const ccomplex_t IT_0188 = (-2)*conj(IT_0184);
    const ccomplex_t IT_0189 = IT_0187 + IT_0188;
    const ccomplex_t IT_0190 = 0.5*IT_0032 + 0.5*IT_0044 + 0.5*IT_0108;
    const ccomplex_t IT_0191 = 2*IT_0184;
    const ccomplex_t IT_0192 = (-2)*IT_0182;
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = 2*conj(IT_0184);
    const ccomplex_t IT_0195 = (-2)*conj(IT_0182);
    const ccomplex_t IT_0196 = IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = m_N_2*IT_0005;
    const ccomplex_t IT_0198 = IT_0059*IT_0116;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0200 = IT_0069*IT_0198*IT_0199;
    const ccomplex_t IT_0201 = IT_0066*IT_0129;
    const ccomplex_t IT_0202 = cpow((-2)*s_13 + IT_0000 + IT_0068 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0203 = IT_0199*IT_0201*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0205 = cpow((-2)*s_13 + IT_0000 + IT_0068 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0206 = IT_0079*IT_0152;
    const ccomplex_t IT_0207 = IT_0204*IT_0205*IT_0206;
    const ccomplex_t IT_0208 = IT_0092*IT_0139;
    const ccomplex_t IT_0209 = IT_0072*IT_0204*IT_0208;
    const ccomplex_t IT_0210 = -IT_0200 + -IT_0203 + -IT_0207 + -IT_0209;
    const ccomplex_t IT_0211 = IT_0130*IT_0202;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = IT_0153*IT_0205;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = IT_0071 + IT_0095 + IT_0212 + IT_0214;
    const ccomplex_t IT_0216 = m_N_2*IT_0215;
    const ccomplex_t IT_0217 = conj(IT_0210) + conj(IT_0216);
    const ccomplex_t IT_0218 = IT_0210 + IT_0216;
    const ccomplex_t IT_0219 = IT_0198*IT_0199*IT_0202;
    const ccomplex_t IT_0220 = IT_0069*IT_0199*IT_0201;
    const ccomplex_t IT_0221 = IT_0204*IT_0205*IT_0208;
    const ccomplex_t IT_0222 = IT_0072*IT_0204*IT_0206;
    const ccomplex_t IT_0223 = IT_0219 + IT_0220 + IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = IT_0067*IT_0202;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = IT_0093*IT_0205;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = IT_0132 + IT_0155 + IT_0225 + IT_0227;
    const ccomplex_t IT_0229 = m_N_2*IT_0228;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = conj(IT_0223) + conj(IT_0230);
    const ccomplex_t IT_0232 = 2*IT_0109;
    const ccomplex_t IT_0233 = (-2)*IT_0190;
    const ccomplex_t IT_0234 = IT_0232 + IT_0233;
    const ccomplex_t IT_0235 = 2*conj(IT_0109);
    const ccomplex_t IT_0236 = (-2)*conj(IT_0190);
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = s_12*s_14;
    const ccomplex_t IT_0239 = s_24*IT_0000;
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = IT_0238 + IT_0240;
    const ccomplex_t IT_0242 = -IT_0132 + -IT_0155;
    const ccomplex_t IT_0243 = IT_0212 + IT_0214;
    const ccomplex_t IT_0244 = conj(IT_0242) + conj(IT_0243);
    const ccomplex_t IT_0245 = IT_0242 + IT_0243;
    const ccomplex_t IT_0246 = IT_0225 + IT_0227;
    const ccomplex_t IT_0247 = -IT_0071 + -IT_0095;
    const ccomplex_t IT_0248 = conj(IT_0246) + conj(IT_0247);
    const ccomplex_t IT_0249 = IT_0246 + IT_0247;
    const ccomplex_t IT_0250 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = s_12 + IT_0251;
    const ccomplex_t IT_0253 = IT_0000*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0032 + 0.5*IT_0044 + 0.5*IT_0108 + (-2)
      *IT_0225 + (-2)*IT_0227;
    const ccomplex_t IT_0255 = (-0.5)*IT_0032 + (-0.5)*IT_0044 + (-0.5)*IT_0108;
    const ccomplex_t IT_0256 = 0.5*IT_0166 + 0.5*IT_0173 + 0.5*IT_0181;
    const ccomplex_t IT_0257 = 2*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0166 + (-0.5)*IT_0173 + (-0.5)
      *IT_0181 + (-2)*IT_0212 + (-2)*IT_0214;
    const ccomplex_t IT_0259 = (-2)*IT_0258;
    const ccomplex_t IT_0260 = IT_0257 + IT_0259;
    const ccomplex_t IT_0261 = 2*conj(IT_0256);
    const ccomplex_t IT_0262 = (-2)*conj(IT_0258);
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0258;
    const ccomplex_t IT_0265 = (-2)*IT_0256;
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = 2*conj(IT_0258);
    const ccomplex_t IT_0268 = (-2)*conj(IT_0256);
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = pow(s_24, 2);
    const ccomplex_t IT_0271 = IT_0001*IT_0270;
    const ccomplex_t IT_0272 = -IT_0271;
    const ccomplex_t IT_0273 = IT_0000 + IT_0272;
    const ccomplex_t IT_0274 = IT_0000*IT_0273;
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = s_12*s_24;
    const ccomplex_t IT_0277 = s_14*IT_0000;
    const ccomplex_t IT_0278 = -IT_0277;
    const ccomplex_t IT_0279 = IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = 2*IT_0254;
    const ccomplex_t IT_0281 = (-2)*IT_0255;
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = 2*conj(IT_0254);
    const ccomplex_t IT_0284 = (-2)*conj(IT_0255);
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = pow(m_W, 2);
    const ccomplex_t IT_0287 = IT_0000*IT_0286;
    const ccomplex_t IT_0288 = (-3)*IT_0287;
    const ccomplex_t IT_0289 = 2*conj(IT_0246);
    const ccomplex_t IT_0290 = 3*IT_0287;
    const ccomplex_t IT_0291 = s_12*IT_0286;
    const ccomplex_t IT_0292 = pow(m_W, 4);
    const ccomplex_t IT_0293 = s_12*IT_0292;
    const ccomplex_t IT_0294 = IT_0001*IT_0293;
    const ccomplex_t IT_0295 = (-0.5)*IT_0294;
    const ccomplex_t IT_0296 = s_14*s_24;
    const ccomplex_t IT_0297 = (-2)*IT_0296;
    const ccomplex_t IT_0298 = IT_0291 + IT_0295 + IT_0297;
    const ccomplex_t IT_0299 = (-0.25)*IT_0294;
    const ccomplex_t IT_0300 = IT_0296 + IT_0299;
    const ccomplex_t IT_0301 = s_24*IT_0286;
    const ccomplex_t IT_0302 = IT_0001*IT_0301;
    const ccomplex_t IT_0303 = 0.5*IT_0302;
    const ccomplex_t IT_0304 = s_24 + IT_0303;
    const ccomplex_t IT_0305 = m_N_2*IT_0304;
    const ccomplex_t IT_0306 = 2*IT_0305;
    const ccomplex_t IT_0307 = (-2)*IT_0223;
    const ccomplex_t IT_0308 = (-2)*IT_0230;
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0246;
    const ccomplex_t IT_0311 = s_14*IT_0286;
    const ccomplex_t IT_0312 = IT_0001*IT_0311;
    const ccomplex_t IT_0313 = (-0.25)*IT_0312;
    const ccomplex_t IT_0314 = s_14 + IT_0313;
    const ccomplex_t IT_0315 = m_N_2*IT_0314;
    const ccomplex_t IT_0316 = (-4)*IT_0315;
    const ccomplex_t IT_0317 = 2*IT_0210;
    const ccomplex_t IT_0318 = (-2)*conj(IT_0243)*(IT_0247*IT_0290 + 2*IT_0243
      *IT_0298 + (-4)*IT_0242*IT_0300 + (-0.5)*IT_0306*IT_0309 + (-0.5)*IT_0288
      *IT_0310 + -IT_0316*(IT_0216 + 0.5*IT_0317));
    const ccomplex_t IT_0319 = (-2)*conj(IT_0223);
    const ccomplex_t IT_0320 = (-2)*conj(IT_0230);
    const ccomplex_t IT_0321 = IT_0319 + IT_0320;
    const ccomplex_t IT_0322 = 2*conj(IT_0210);
    const ccomplex_t IT_0323 = s_12*IT_0252;
    const ccomplex_t IT_0324 = 2*IT_0190;
    const ccomplex_t IT_0325 = 0.5*IT_0237*IT_0254 + 0.5*IT_0234*conj(IT_0254)
       + -(conj(IT_0109) + -conj(IT_0190))*IT_0255 + 0.5*IT_0196*IT_0256 + 0.5
      *IT_0193*conj(IT_0256) + 0.5*IT_0189*IT_0258 + 0.5*IT_0186*conj(IT_0258) +
       -conj(IT_0255)*(IT_0109 + (-0.5)*IT_0324);
    const ccomplex_t IT_0326 = IT_0291 + IT_0297;
    const ccomplex_t IT_0327 = IT_0001*IT_0326;
    const ccomplex_t IT_0328 = (-0.5)*IT_0327;
    const ccomplex_t IT_0329 = s_12 + IT_0328;
    const ccomplex_t IT_0330 = 4*IT_0210;
    const ccomplex_t IT_0331 = 4*IT_0230;
    const ccomplex_t IT_0332 = 6*conj(IT_0210);
    const ccomplex_t IT_0333 = 0.5*IT_0312;
    const ccomplex_t IT_0334 = s_14 + IT_0333;
    const ccomplex_t IT_0335 = m_N_2*IT_0334;
    const ccomplex_t IT_0336 = 2*IT_0335;
    const ccomplex_t IT_0337 = (-2)*IT_0242;
    const ccomplex_t IT_0338 = (-2)*conj(IT_0242);
    const ccomplex_t IT_0339 = (-0.25)*IT_0302;
    const ccomplex_t IT_0340 = s_24 + IT_0339;
    const ccomplex_t IT_0341 = m_N_2*IT_0340;
    const ccomplex_t IT_0342 = (-4)*IT_0341;
    const ccomplex_t IT_0343 = 2*IT_0223;
    const ccomplex_t IT_0344 = 2*IT_0230;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = (-2)*IT_0210;
    const ccomplex_t IT_0347 = m_N_2*IT_0273;
    const ccomplex_t IT_0348 = -IT_0347;
    const ccomplex_t IT_0349 = (-2)*conj(IT_0210);
    const ccomplex_t IT_0350 = 2*conj(IT_0223);
    const ccomplex_t IT_0351 = 2*conj(IT_0230);
    const ccomplex_t IT_0352 = IT_0350 + IT_0351;
    const ccomplex_t IT_0353 = s_12*IT_0005;
    const ccomplex_t IT_0354 = s_12*IT_0273;
    const ccomplex_t IT_0355 = -IT_0354;
    const ccomplex_t IT_0356 = (-2)*(IT_0254*conj(IT_0255) + IT_0255*(conj
      (IT_0254) + -conj(IT_0255)) + (-0.5)*conj(IT_0256)*IT_0257 + (-0.5)*conj
      (IT_0258)*IT_0266 + (-0.5)*IT_0258*IT_0268 + (-0.5)*conj(IT_0254)*IT_0280)
      *IT_0355;
    const ccomplex_t IT_0357 = m_N_2*IT_0252;
    const ccomplex_t IT_0358 = conj(IT_0216)*(IT_0232 + IT_0233 + IT_0257 +
       IT_0259) + conj(IT_0210)*IT_0260 + IT_0216*(IT_0235 + IT_0236 + IT_0261 +
       IT_0262) + IT_0210*IT_0263 + IT_0231*(IT_0191 + IT_0192 + IT_0280 +
       IT_0281) + (IT_0223 + IT_0230)*(IT_0194 + IT_0195 + IT_0283 + IT_0284) +
       conj(IT_0109)*IT_0317 + IT_0109*IT_0322 + conj(IT_0190)*IT_0346 + IT_0190
      *IT_0349;
    const ccomplex_t IT_0359 = (-4)*IT_0246;
    const ccomplex_t IT_0360 = IT_0006*(conj(IT_0109)*IT_0186 + IT_0109
      *IT_0189 + conj(IT_0190)*IT_0193 + IT_0190*IT_0196) + IT_0197*(IT_0193
      *IT_0217 + IT_0196*IT_0218 + IT_0231*IT_0234 + (IT_0223 + IT_0230)*IT_0237
      ) + IT_0241*(IT_0186*IT_0244 + IT_0189*IT_0245 + IT_0234*IT_0248 + IT_0237
      *IT_0249) + IT_0253*(IT_0189*IT_0254 + IT_0186*conj(IT_0254) + IT_0196
      *IT_0255 + IT_0193*conj(IT_0255) + conj(IT_0190)*IT_0260 + IT_0190*IT_0263
       + conj(IT_0109)*IT_0266 + IT_0109*IT_0269) + (conj(IT_0255)*IT_0260 +
       IT_0255*IT_0263 + conj(IT_0254)*IT_0266 + IT_0254*IT_0269)*IT_0275 +
       IT_0279*(IT_0244*IT_0266 + IT_0245*IT_0269 + IT_0248*IT_0282 + IT_0249
      *IT_0285) + IT_0243*IT_0288*IT_0289 + IT_0318 + IT_0306*(IT_0216*IT_0289 +
       conj(IT_0216)*IT_0310 + conj(IT_0246)*IT_0317 + IT_0243*IT_0321 + IT_0246
      *IT_0322) + IT_0316*(conj(IT_0246)*IT_0309 + IT_0246*IT_0321 + 2*IT_0243*
      (conj(IT_0216) + 0.5*IT_0322)) + 2*IT_0323*IT_0325 + 4*IT_0329*(IT_0217*
      (IT_0216 + 0.25*IT_0330) + IT_0231*(IT_0223 + 0.25*IT_0331)) + 6*IT_0000*
      (IT_0218*IT_0231 + (IT_0223 + IT_0230)*(conj(IT_0216) + 0.166666666666667
      *IT_0332)) + (conj(IT_0246)*IT_0290 + conj(IT_0216)*IT_0336)*IT_0337 +
       IT_0216*IT_0336*IT_0338 + 2*conj(IT_0242)*(IT_0247*IT_0288 + (-2)*IT_0242
      *IT_0298 + 4*IT_0243*IT_0300 + 0.5*IT_0342*IT_0345 + 0.5*IT_0336*IT_0346) 
      + IT_0348*(IT_0231*IT_0260 + (IT_0223 + IT_0230)*IT_0263 + conj(IT_0216)
      *IT_0282 + IT_0216*IT_0285 + conj(IT_0254)*IT_0317 + IT_0254*IT_0322 +
       conj(IT_0255)*IT_0346 + IT_0255*IT_0349) + 8*IT_0247*(conj(IT_0246)
      *IT_0300 + (-0.25)*IT_0342*(conj(IT_0216) + (-0.5)*IT_0349) + 0.125
      *IT_0336*IT_0352) + IT_0242*(IT_0336*IT_0349 + IT_0342*IT_0352) + (conj
      (IT_0182)*IT_0183 + conj(IT_0184)*IT_0193 + IT_0184*IT_0195 + conj(IT_0109
      )*IT_0234 + IT_0109*IT_0236 + conj(IT_0190)*IT_0324)*IT_0353 + IT_0356 +
       IT_0357*IT_0358 + (conj(IT_0246)*IT_0298 + (-0.25)*IT_0290*IT_0338)
      *IT_0359 + 2*conj(IT_0247)*(IT_0242*IT_0288 + -IT_0243*IT_0290 + (-2)
      *IT_0247*IT_0298 + 0.5*IT_0336*IT_0345 + -IT_0342*(IT_0216 + (-0.5)
      *IT_0346) + -IT_0300*IT_0359);
    return create_ccomplex_return(IT_0360);
}

