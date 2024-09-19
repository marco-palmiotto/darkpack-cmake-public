#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_Hp_W(
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
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_N_4, 2);
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = -IT_0000;
    const ccomplex_t IT_0005 = IT_0003 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = cos(alpha);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = sin(alpha);
    const ccomplex_t IT_0014 = IT_0009*IT_0013;
    const ccomplex_t IT_0015 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = sin(theta_W);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0007*IT_0018;
    const ccomplex_t IT_0020 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0018;
    const ccomplex_t IT_0023 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0012 + IT_0016 + 
      -IT_0021 + -IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = cos(beta);
    const ccomplex_t IT_0028 = sin(beta);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*IT_0018*(IT_0007
      *IT_0027 + IT_0013*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0026*IT_0030;
    const ccomplex_t IT_0032 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*-1)*e_em*IT_0018*(IT_0013
      *IT_0027 + -IT_0007*IT_0028);
    const ccomplex_t IT_0036 = IT_0010*IT_0015;
    const ccomplex_t IT_0037 = IT_0011*IT_0014;
    const ccomplex_t IT_0038 = IT_0019*IT_0023;
    const ccomplex_t IT_0039 = IT_0020*IT_0022;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + -IT_0037 + 
      -IT_0038 + IT_0039);
    const ccomplex_t IT_0041 = IT_0035*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = IT_0009*IT_0027;
    const ccomplex_t IT_0046 = IT_0011*IT_0045;
    const ccomplex_t IT_0047 = IT_0009*IT_0028;
    const ccomplex_t IT_0048 = IT_0015*IT_0047;
    const ccomplex_t IT_0049 = IT_0018*IT_0027;
    const ccomplex_t IT_0050 = IT_0020*IT_0049;
    const ccomplex_t IT_0051 = IT_0018*IT_0028;
    const ccomplex_t IT_0052 = IT_0023*IT_0051;
    const ccomplex_t IT_0053 = IT_0046 + -IT_0048 + -IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = e_em*IT_0018;
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = IT_0056*IT_0058;
    const ccomplex_t IT_0060 = 0.5*IT_0034 + 0.5*IT_0044 + 0.5*IT_0059;
    const ccomplex_t IT_0061 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0062 = IT_0045*IT_0061;
    const ccomplex_t IT_0063 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0064 = IT_0047*IT_0063;
    const ccomplex_t IT_0065 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0066 = IT_0049*IT_0065;
    const ccomplex_t IT_0067 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0068 = IT_0051*IT_0067;
    const ccomplex_t IT_0069 = -IT_0062 + IT_0064 + IT_0066 + -IT_0068;
    const ccomplex_t IT_0070 = IT_0055*IT_0069;
    const ccomplex_t IT_0071 = IT_0058*IT_0070;
    const ccomplex_t IT_0072 = IT_0010*IT_0063;
    const ccomplex_t IT_0073 = IT_0014*IT_0061;
    const ccomplex_t IT_0074 = IT_0019*IT_0067;
    const ccomplex_t IT_0075 = IT_0022*IT_0065;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*(IT_0072 + -IT_0073 + 
      -IT_0074 + IT_0075);
    const ccomplex_t IT_0077 = IT_0035*IT_0076;
    const ccomplex_t IT_0078 = IT_0043*IT_0077;
    const ccomplex_t IT_0079 = IT_0010*IT_0061;
    const ccomplex_t IT_0080 = IT_0014*IT_0063;
    const ccomplex_t IT_0081 = IT_0019*IT_0065;
    const ccomplex_t IT_0082 = IT_0022*IT_0067;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + IT_0080 + 
      -IT_0081 + -IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = IT_0030*IT_0084;
    const ccomplex_t IT_0086 = IT_0033*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0071 + (-0.5)*IT_0078 + (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = 2*IT_0087;
    const ccomplex_t IT_0089 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0090 = IT_0018*IT_0089;
    const ccomplex_t IT_0091 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0092 = IT_0018*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + 0.5*IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = N_d4*e_em*U_Wm1;
    const ccomplex_t IT_0097 = IT_0051*IT_0096;
    const ccomplex_t IT_0098 = (-2)*IT_0097;
    const ccomplex_t IT_0099 = N_B4*U_d1*e_em;
    const ccomplex_t IT_0100 = IT_0047*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = N_W4*U_d1*e_em;
    const ccomplex_t IT_0103 = IT_0051*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = IT_0101 + IT_0104;
    const ccomplex_t IT_0106 = IT_0098 + IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0107;
    const ccomplex_t IT_0109 = IT_0095*IT_0108;
    const ccomplex_t IT_0110 = pow(m_Hp, 2);
    const ccomplex_t IT_0111 = cpow((-2)*s_23 + IT_0000 + IT_0110 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0112 = IT_0109*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0115 = IT_0018*IT_0114;
    const ccomplex_t IT_0116 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0117 = IT_0018*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + 0.5*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = N_d4*e_em*U_Wm2;
    const ccomplex_t IT_0122 = IT_0051*IT_0121;
    const ccomplex_t IT_0123 = (-2)*IT_0122;
    const ccomplex_t IT_0124 = N_B4*U_d2*e_em;
    const ccomplex_t IT_0125 = IT_0047*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = N_W4*U_d2*e_em;
    const ccomplex_t IT_0128 = IT_0051*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = IT_0126 + IT_0129;
    const ccomplex_t IT_0131 = IT_0123 + IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0120*IT_0133;
    const ccomplex_t IT_0135 = cpow((-2)*s_23 + IT_0000 + IT_0110 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0136 = IT_0134*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = 0.5*IT_0071 + 0.5*IT_0078 + 0.5*IT_0086 + 2
      *IT_0113 + 2*IT_0137;
    const ccomplex_t IT_0139 = (-2)*IT_0138;
    const ccomplex_t IT_0140 = IT_0088 + IT_0139;
    const ccomplex_t IT_0141 = 2*conj(IT_0087);
    const ccomplex_t IT_0142 = (-2)*conj(IT_0138);
    const ccomplex_t IT_0143 = IT_0141 + IT_0142;
    const ccomplex_t IT_0144 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0145 = IT_0018*IT_0144;
    const ccomplex_t IT_0146 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0147 = IT_0018*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + (-0.5)*IT_0148);
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = conj(N_B4)*conj(V_u2)*e_em;
    const ccomplex_t IT_0152 = IT_0045*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W4)*conj(V_u2)*e_em;
    const ccomplex_t IT_0155 = IT_0049*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = IT_0153 + IT_0156;
    const ccomplex_t IT_0158 = conj(N_u4)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0159 = IT_0049*IT_0158;
    const ccomplex_t IT_0160 = 2*IT_0159;
    const ccomplex_t IT_0161 = IT_0157 + IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = IT_0150*IT_0163;
    const ccomplex_t IT_0165 = IT_0135*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = conj(N_u4)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0168 = IT_0049*IT_0167;
    const ccomplex_t IT_0169 = 2*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B4)*conj(V_u1)*e_em;
    const ccomplex_t IT_0171 = IT_0045*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W4)*conj(V_u1)*e_em;
    const ccomplex_t IT_0174 = IT_0049*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = IT_0172 + IT_0175;
    const ccomplex_t IT_0177 = IT_0169 + IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0181 = IT_0018*IT_0180;
    const ccomplex_t IT_0182 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0183 = IT_0018*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0181 + (-0.5)*IT_0184);
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = IT_0179*IT_0186;
    const ccomplex_t IT_0188 = IT_0111*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = (-0.5)*IT_0034 + (-0.5)*IT_0044 + (-0.5)
      *IT_0059 + 2*IT_0166 + 2*IT_0189;
    const ccomplex_t IT_0191 = 2*IT_0138;
    const ccomplex_t IT_0192 = (-2)*IT_0087;
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = 2*conj(IT_0138);
    const ccomplex_t IT_0195 = (-2)*conj(IT_0087);
    const ccomplex_t IT_0196 = IT_0194 + IT_0195;
    const ccomplex_t IT_0197 = s_12*s_14;
    const ccomplex_t IT_0198 = s_24*IT_0000;
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = IT_0197 + IT_0199;
    const ccomplex_t IT_0201 = cpow((-2)*s_13 + IT_0000 + IT_0110 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0202 = IT_0134*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = cpow((-2)*s_13 + IT_0000 + IT_0110 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0205 = IT_0109*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = IT_0203 + IT_0206;
    const ccomplex_t IT_0208 = -IT_0113 + -IT_0137;
    const ccomplex_t IT_0209 = conj(IT_0207) + conj(IT_0208);
    const ccomplex_t IT_0210 = IT_0207 + IT_0208;
    const ccomplex_t IT_0211 = -IT_0166 + -IT_0189;
    const ccomplex_t IT_0212 = IT_0187*IT_0204;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = IT_0164*IT_0201;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0213 + IT_0215;
    const ccomplex_t IT_0217 = conj(IT_0211) + conj(IT_0216);
    const ccomplex_t IT_0218 = 2*IT_0190;
    const ccomplex_t IT_0219 = (-2)*IT_0060;
    const ccomplex_t IT_0220 = IT_0218 + IT_0219;
    const ccomplex_t IT_0221 = 2*conj(IT_0190);
    const ccomplex_t IT_0222 = (-2)*conj(IT_0060);
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = s_12 + IT_0225;
    const ccomplex_t IT_0227 = IT_0000*IT_0226;
    const ccomplex_t IT_0228 = (-0.5)*IT_0034 + (-0.5)*IT_0044 + (-0.5)*IT_0059;
    const ccomplex_t IT_0229 = 0.5*IT_0034 + 0.5*IT_0044 + 0.5*IT_0059 + (-2)
      *IT_0213 + (-2)*IT_0215;
    const ccomplex_t IT_0230 = (-0.5)*IT_0071 + (-0.5)*IT_0078 + (-0.5)
      *IT_0086 + (-2)*IT_0203 + (-2)*IT_0206;
    const ccomplex_t IT_0231 = 2*IT_0230;
    const ccomplex_t IT_0232 = 0.5*IT_0071 + 0.5*IT_0078 + 0.5*IT_0086;
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = IT_0231 + IT_0233;
    const ccomplex_t IT_0235 = 2*conj(IT_0230);
    const ccomplex_t IT_0236 = (-2)*conj(IT_0232);
    const ccomplex_t IT_0237 = IT_0235 + IT_0236;
    const ccomplex_t IT_0238 = 2*IT_0232;
    const ccomplex_t IT_0239 = (-2)*IT_0230;
    const ccomplex_t IT_0240 = IT_0238 + IT_0239;
    const ccomplex_t IT_0241 = 2*conj(IT_0232);
    const ccomplex_t IT_0242 = (-2)*conj(IT_0230);
    const ccomplex_t IT_0243 = IT_0241 + IT_0242;
    const ccomplex_t IT_0244 = pow(s_24, 2);
    const ccomplex_t IT_0245 = IT_0001*IT_0244;
    const ccomplex_t IT_0246 = -IT_0245;
    const ccomplex_t IT_0247 = IT_0000 + IT_0246;
    const ccomplex_t IT_0248 = IT_0000*IT_0247;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = s_12*s_24;
    const ccomplex_t IT_0251 = s_14*IT_0000;
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = IT_0250 + IT_0252;
    const ccomplex_t IT_0254 = 2*IT_0229;
    const ccomplex_t IT_0255 = (-2)*IT_0228;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 2*conj(IT_0229);
    const ccomplex_t IT_0258 = (-2)*conj(IT_0228);
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = m_N_4*IT_0247;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0095*IT_0179;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0264 = IT_0204*IT_0262*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0266 = IT_0133*IT_0150;
    const ccomplex_t IT_0267 = IT_0135*IT_0265*IT_0266;
    const ccomplex_t IT_0268 = IT_0120*IT_0163;
    const ccomplex_t IT_0269 = IT_0201*IT_0265*IT_0268;
    const ccomplex_t IT_0270 = IT_0108*IT_0186;
    const ccomplex_t IT_0271 = IT_0111*IT_0263*IT_0270;
    const ccomplex_t IT_0272 = IT_0264 + IT_0267 + IT_0269 + IT_0271;
    const ccomplex_t IT_0273 = IT_0113 + IT_0137 + IT_0213 + IT_0215;
    const ccomplex_t IT_0274 = m_N_4*IT_0273;
    const ccomplex_t IT_0275 = -IT_0274;
    const ccomplex_t IT_0276 = conj(IT_0272) + conj(IT_0275);
    const ccomplex_t IT_0277 = IT_0272 + IT_0275;
    const ccomplex_t IT_0278 = IT_0111*IT_0262*IT_0263;
    const ccomplex_t IT_0279 = IT_0204*IT_0263*IT_0270;
    const ccomplex_t IT_0280 = IT_0135*IT_0265*IT_0268;
    const ccomplex_t IT_0281 = IT_0201*IT_0265*IT_0266;
    const ccomplex_t IT_0282 = -IT_0278 + -IT_0279 + -IT_0280 + -IT_0281;
    const ccomplex_t IT_0283 = IT_0166 + IT_0189 + IT_0203 + IT_0206;
    const ccomplex_t IT_0284 = m_N_4*IT_0283;
    const ccomplex_t IT_0285 = conj(IT_0282) + conj(IT_0284);
    const ccomplex_t IT_0286 = IT_0282 + IT_0284;
    const ccomplex_t IT_0287 = pow(m_W, 2);
    const ccomplex_t IT_0288 = s_24*IT_0287;
    const ccomplex_t IT_0289 = IT_0001*IT_0288;
    const ccomplex_t IT_0290 = 0.5*IT_0289;
    const ccomplex_t IT_0291 = s_24 + IT_0290;
    const ccomplex_t IT_0292 = m_N_4*IT_0291;
    const ccomplex_t IT_0293 = 2*IT_0292;
    const ccomplex_t IT_0294 = (-2)*IT_0272;
    const ccomplex_t IT_0295 = (-2)*IT_0275;
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = s_14*IT_0287;
    const ccomplex_t IT_0298 = IT_0001*IT_0297;
    const ccomplex_t IT_0299 = (-0.25)*IT_0298;
    const ccomplex_t IT_0300 = s_14 + IT_0299;
    const ccomplex_t IT_0301 = m_N_4*IT_0300;
    const ccomplex_t IT_0302 = (-4)*IT_0301;
    const ccomplex_t IT_0303 = 2*IT_0282;
    const ccomplex_t IT_0304 = s_12*IT_0287;
    const ccomplex_t IT_0305 = pow(m_W, 4);
    const ccomplex_t IT_0306 = s_12*IT_0305;
    const ccomplex_t IT_0307 = IT_0001*IT_0306;
    const ccomplex_t IT_0308 = (-0.5)*IT_0307;
    const ccomplex_t IT_0309 = s_14*s_24;
    const ccomplex_t IT_0310 = (-2)*IT_0309;
    const ccomplex_t IT_0311 = IT_0304 + IT_0308 + IT_0310;
    const ccomplex_t IT_0312 = (-0.25)*IT_0307;
    const ccomplex_t IT_0313 = IT_0309 + IT_0312;
    const ccomplex_t IT_0314 = (-2)*conj(IT_0272);
    const ccomplex_t IT_0315 = (-2)*conj(IT_0275);
    const ccomplex_t IT_0316 = IT_0314 + IT_0315;
    const ccomplex_t IT_0317 = 2*conj(IT_0282);
    const ccomplex_t IT_0318 = IT_0304 + IT_0310;
    const ccomplex_t IT_0319 = IT_0001*IT_0318;
    const ccomplex_t IT_0320 = (-0.5)*IT_0319;
    const ccomplex_t IT_0321 = s_12 + IT_0320;
    const ccomplex_t IT_0322 = 4*IT_0282;
    const ccomplex_t IT_0323 = 6*conj(IT_0272);
    const ccomplex_t IT_0324 = m_N_4*IT_0005;
    const ccomplex_t IT_0325 = 2*IT_0272;
    const ccomplex_t IT_0326 = 2*conj(IT_0272);
    const ccomplex_t IT_0327 = (-0.25)*IT_0289;
    const ccomplex_t IT_0328 = s_24 + IT_0327;
    const ccomplex_t IT_0329 = m_N_4*IT_0328;
    const ccomplex_t IT_0330 = (-4)*IT_0329;
    const ccomplex_t IT_0331 = 2*IT_0275;
    const ccomplex_t IT_0332 = IT_0325 + IT_0331;
    const ccomplex_t IT_0333 = 0.5*IT_0298;
    const ccomplex_t IT_0334 = s_14 + IT_0333;
    const ccomplex_t IT_0335 = m_N_4*IT_0334;
    const ccomplex_t IT_0336 = 2*IT_0335;
    const ccomplex_t IT_0337 = (-2)*IT_0282;
    const ccomplex_t IT_0338 = IT_0000*IT_0287;
    const ccomplex_t IT_0339 = 3*IT_0338;
    const ccomplex_t IT_0340 = (-3)*IT_0338;
    const ccomplex_t IT_0341 = 2*conj(IT_0275);
    const ccomplex_t IT_0342 = IT_0326 + IT_0341;
    const ccomplex_t IT_0343 = (-2)*conj(IT_0282);
    const ccomplex_t IT_0344 = s_12*IT_0005;
    const ccomplex_t IT_0345 = 2*IT_0060;
    const ccomplex_t IT_0346 = (-2)*IT_0190;
    const ccomplex_t IT_0347 = IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = (-2)*conj(IT_0190);
    const ccomplex_t IT_0349 = s_12*IT_0226;
    const ccomplex_t IT_0350 = m_N_4*IT_0226;
    const ccomplex_t IT_0351 = IT_0143*IT_0272 + IT_0140*conj(IT_0272) + 
      (IT_0141 + IT_0142 + IT_0257 + IT_0258)*IT_0275 + (IT_0088 + IT_0139 +
       IT_0254 + IT_0255)*conj(IT_0275) + (IT_0218 + IT_0219 + IT_0238 + IT_0239
      )*IT_0285 + (IT_0221 + IT_0222 + IT_0241 + IT_0242)*IT_0286 + conj(IT_0228
      )*IT_0294 + IT_0228*IT_0314 + conj(IT_0229)*IT_0325 + IT_0229*IT_0326;
    const ccomplex_t IT_0352 = s_12*IT_0247;
    const ccomplex_t IT_0353 = -IT_0352;
    const ccomplex_t IT_0354 = 2*IT_0207;
    const ccomplex_t IT_0355 = (-2)*IT_0208;
    const ccomplex_t IT_0356 = 2*conj(IT_0207);
    const ccomplex_t IT_0357 = (-2)*conj(IT_0208);
    const ccomplex_t IT_0358 = 0.5*IT_0339;
    const ccomplex_t IT_0359 = 0.5*IT_0340;
    const ccomplex_t IT_0360 = 0.5*IT_0302;
    const ccomplex_t IT_0361 = (-4)*conj(IT_0216)*((-0.5)*IT_0293*(IT_0284 +
       0.5*IT_0303) + IT_0216*IT_0311 + (-2)*IT_0211*IT_0313 + (-0.5)*IT_0355
      *IT_0358 + (-0.5)*IT_0354*IT_0359 + (-0.5)*IT_0296*IT_0360);
    const ccomplex_t IT_0362 = IT_0006*(conj(IT_0060)*IT_0140 + IT_0060
      *IT_0143 + conj(IT_0190)*IT_0193 + IT_0190*IT_0196) + IT_0200*(IT_0193
      *IT_0209 + IT_0196*IT_0210 + IT_0217*IT_0220 + (IT_0211 + IT_0216)*IT_0223
      ) + IT_0227*(IT_0143*IT_0228 + IT_0140*conj(IT_0228) + IT_0196*IT_0229 +
       IT_0193*conj(IT_0229) + conj(IT_0190)*IT_0234 + IT_0190*IT_0237 + conj
      (IT_0060)*IT_0240 + IT_0060*IT_0243) + (conj(IT_0229)*IT_0234 + IT_0229
      *IT_0237 + conj(IT_0228)*IT_0240 + IT_0228*IT_0243)*IT_0249 + IT_0253*
      (IT_0209*IT_0234 + IT_0210*IT_0237 + IT_0217*IT_0256 + (IT_0211 + IT_0216)
      *IT_0259) + IT_0261*(IT_0240*IT_0276 + IT_0243*IT_0277 + IT_0256*IT_0285 +
       IT_0259*IT_0286) + conj(IT_0207)*(IT_0293*IT_0296 + IT_0302*IT_0303) + (
      -4)*IT_0207*(conj(IT_0207)*IT_0311 + (-2)*conj(IT_0208)*IT_0313 + (-0.25)
      *IT_0293*IT_0316 + (-0.25)*IT_0302*IT_0317) + 4*IT_0321*(IT_0276*IT_0277 +
       IT_0285*(IT_0284 + 0.25*IT_0322)) + 6*IT_0000*(IT_0277*IT_0285 + IT_0286*
      (conj(IT_0275) + 0.166666666666667*IT_0323)) + IT_0324*(IT_0223*IT_0275 +
       IT_0220*conj(IT_0275) + IT_0140*IT_0285 + IT_0143*IT_0286 + conj(IT_0060)
      *IT_0294 + IT_0060*IT_0314 + conj(IT_0190)*IT_0325 + IT_0190*IT_0326) +
       conj(IT_0208)*(IT_0330*IT_0332 + IT_0336*IT_0337) + (-4)*conj(IT_0211)*
      (IT_0211*IT_0311 + (-2)*IT_0216*IT_0313 + (-0.25)*IT_0332*IT_0336 + 1./2
      *IT_0330*(IT_0284 + (-0.5)*IT_0337) + 1./2*IT_0207*IT_0339 + -1./2*IT_0208
      *IT_0340) + (-4)*IT_0208*(conj(IT_0208)*IT_0311 + (-2)*conj(IT_0207)
      *IT_0313 + (-0.25)*IT_0330*IT_0342 + (-0.25)*IT_0336*IT_0343) + (-2)
      *IT_0211*(conj(IT_0207)*IT_0339 + -conj(IT_0208)*IT_0340 + (-0.5)*IT_0336
      *IT_0342 + IT_0330*(conj(IT_0284) + (-0.5)*IT_0343)) + IT_0344*(conj
      (IT_0087)*IT_0088 + conj(IT_0138)*IT_0193 + IT_0138*IT_0195 + conj(IT_0190
      )*IT_0218 + conj(IT_0060)*IT_0347 + IT_0060*IT_0348) + (IT_0223*IT_0229 +
       IT_0220*conj(IT_0229) + IT_0196*IT_0230 + IT_0193*conj(IT_0230) + IT_0143
      *IT_0232 + IT_0140*conj(IT_0232) + conj(IT_0228)*IT_0347 + 2*IT_0228*(conj
      (IT_0060) + 0.5*IT_0348))*IT_0349 + IT_0350*IT_0351 + 2*(IT_0228*conj
      (IT_0228) + 0.5*conj(IT_0230)*IT_0231 + 0.5*conj(IT_0232)*IT_0240 + 0.5
      *IT_0232*IT_0242 + 0.5*conj(IT_0229)*IT_0256 + 0.5*IT_0229*IT_0258)
      *IT_0353 + conj(IT_0284)*(IT_0302*IT_0354 + IT_0336*IT_0355) + IT_0284*
      (IT_0302*IT_0356 + IT_0336*IT_0357) + 2*IT_0216*(IT_0293*(conj(IT_0284) +
       0.5*IT_0317) + IT_0357*IT_0358 + IT_0356*IT_0359 + IT_0316*IT_0360) +
       IT_0361;
    return create_ccomplex_return(IT_0362);
}

