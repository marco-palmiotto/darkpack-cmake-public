#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_Hp_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_Hp_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0002 = pow(m_Hp, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = IT_0004*IT_0006;
    const ccomplex_t IT_0008 = conj(N_u1)*e_em*conj(V_Wp2);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 2*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = conj(N_B1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = conj(N_W1)*conj(V_u2)*e_em;
    const ccomplex_t IT_0018 = IT_0007*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = IT_0016 + IT_0019;
    const ccomplex_t IT_0021 = IT_0010 + IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = sin(beta);
    const ccomplex_t IT_0025 = IT_0006*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0027 = IT_0025*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = IT_0012*IT_0024;
    const ccomplex_t IT_0030 = conj(N_B1)*conj(U_d2)*e_em;
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = 1.4142135623731*IT_0031;
    const ccomplex_t IT_0033 = IT_0028 + IT_0032;
    const ccomplex_t IT_0034 = conj(N_d1)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0035 = IT_0025*IT_0034;
    const ccomplex_t IT_0036 = (-2)*IT_0035;
    const ccomplex_t IT_0037 = IT_0033 + IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = 0.5*IT_0038;
    const ccomplex_t IT_0040 = IT_0023*IT_0039;
    const ccomplex_t IT_0041 = IT_0001*IT_0003*IT_0040;
    const ccomplex_t IT_0042 = N_d1*e_em*U_Wm2;
    const ccomplex_t IT_0043 = IT_0025*IT_0042;
    const ccomplex_t IT_0044 = (-2)*IT_0043;
    const ccomplex_t IT_0045 = N_B1*U_d2*e_em;
    const ccomplex_t IT_0046 = IT_0029*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = N_W1*U_d2*e_em;
    const ccomplex_t IT_0049 = IT_0025*IT_0048;
    const ccomplex_t IT_0050 = 1.4142135623731*IT_0049;
    const ccomplex_t IT_0051 = IT_0047 + IT_0050;
    const ccomplex_t IT_0052 = IT_0044 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0039*IT_0054;
    const ccomplex_t IT_0056 = IT_0003*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_N_1*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0060 = cpow((-2)*s_13 + IT_0000 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0061 = conj(N_u1)*e_em*conj(V_Wp1);
    const ccomplex_t IT_0062 = IT_0007*IT_0061;
    const ccomplex_t IT_0063 = 2*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0065 = IT_0013*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = conj(N_W1)*conj(V_u1)*e_em;
    const ccomplex_t IT_0068 = IT_0007*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = IT_0066 + IT_0069;
    const ccomplex_t IT_0071 = IT_0063 + IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0075 = IT_0029*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0078 = IT_0025*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = IT_0076 + IT_0079;
    const ccomplex_t IT_0081 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0082 = IT_0025*IT_0081;
    const ccomplex_t IT_0083 = (-2)*IT_0082;
    const ccomplex_t IT_0084 = IT_0080 + IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0085;
    const ccomplex_t IT_0087 = IT_0073*IT_0086;
    const ccomplex_t IT_0088 = IT_0059*IT_0060*IT_0087;
    const ccomplex_t IT_0089 = N_d1*e_em*U_Wm1;
    const ccomplex_t IT_0090 = IT_0025*IT_0089;
    const ccomplex_t IT_0091 = (-2)*IT_0090;
    const ccomplex_t IT_0092 = N_B1*U_d1*e_em;
    const ccomplex_t IT_0093 = IT_0029*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = N_W1*U_d1*e_em;
    const ccomplex_t IT_0096 = IT_0025*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = IT_0094 + IT_0097;
    const ccomplex_t IT_0099 = IT_0091 + IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = 0.5*IT_0100;
    const ccomplex_t IT_0102 = IT_0086*IT_0101;
    const ccomplex_t IT_0103 = IT_0060*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = m_N_1*IT_0104;
    const ccomplex_t IT_0106 = IT_0041 + -IT_0058 + IT_0088 + -IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*-1)*e_em*(IT_0006*IT_0011 + 
      -IT_0005*IT_0012);
    const ccomplex_t IT_0110 = IT_0005*IT_0012;
    const ccomplex_t IT_0111 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = IT_0006*IT_0011;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0116 = IT_0110*IT_0115;
    const ccomplex_t IT_0117 = IT_0113*IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0112 + IT_0114 + 
      -IT_0116 + -IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = IT_0109*IT_0119;
    const ccomplex_t IT_0121 = m_N_1*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0109;
    const ccomplex_t IT_0123 = 0.5*IT_0118;
    const ccomplex_t IT_0124 = IT_0122*IT_0123;
    const ccomplex_t IT_0125 = m_N_1*IT_0124;
    const ccomplex_t IT_0126 = IT_0121 + IT_0125;
    const ccomplex_t IT_0127 = IT_0119*IT_0122;
    const ccomplex_t IT_0128 = m_N_1*IT_0127;
    const ccomplex_t IT_0129 = IT_0109*IT_0123;
    const ccomplex_t IT_0130 = m_N_1*IT_0129;
    const ccomplex_t IT_0131 = -IT_0128 + -IT_0130;
    const ccomplex_t IT_0132 = IT_0126 + IT_0131;
    const ccomplex_t IT_0133 = IT_0108*IT_0132;
    const ccomplex_t IT_0134 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0135 = cos(alpha);
    const ccomplex_t IT_0136 = IT_0012*IT_0135;
    const ccomplex_t IT_0137 = IT_0134*IT_0136;
    const ccomplex_t IT_0138 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0139 = sin(alpha);
    const ccomplex_t IT_0140 = IT_0012*IT_0139;
    const ccomplex_t IT_0141 = IT_0138*IT_0140;
    const ccomplex_t IT_0142 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0143 = IT_0006*IT_0135;
    const ccomplex_t IT_0144 = IT_0142*IT_0143;
    const ccomplex_t IT_0145 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0146 = IT_0006*IT_0139;
    const ccomplex_t IT_0147 = IT_0145*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0137 + IT_0141 + 
      -IT_0144 + -IT_0147);
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = cpow(IT_0011, -2);
    const ccomplex_t IT_0151 = IT_0005*IT_0150;
    const ccomplex_t IT_0152 = -IT_0006 + IT_0151;
    const ccomplex_t IT_0153 = cpow(IT_0024, 3);
    const ccomplex_t IT_0154 = 3*IT_0006 + IT_0151;
    const ccomplex_t IT_0155 = cpow(IT_0024, 2);
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*m_W*e_em*(IT_0135*IT_0152
      *IT_0153 + -IT_0004*(IT_0004*(IT_0004*IT_0139*IT_0152 + IT_0024*IT_0135
      *IT_0154) + -IT_0139*IT_0154*IT_0155));
    const ccomplex_t IT_0157 = 0.5*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0149*IT_0160;
    const ccomplex_t IT_0162 = N_u1*e_em*V_Wp2;
    const ccomplex_t IT_0163 = IT_0007*IT_0162;
    const ccomplex_t IT_0164 = 2*IT_0163;
    const ccomplex_t IT_0165 = N_B1*V_u2*e_em;
    const ccomplex_t IT_0166 = IT_0013*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = N_W1*V_u2*e_em;
    const ccomplex_t IT_0169 = IT_0007*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = IT_0167 + IT_0170;
    const ccomplex_t IT_0172 = IT_0164 + IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = 0.5*IT_0173;
    const ccomplex_t IT_0175 = IT_0023*IT_0174;
    const ccomplex_t IT_0176 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0177 = IT_0175*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = m_N_1*IT_0178;
    const ccomplex_t IT_0180 = IT_0136*IT_0138;
    const ccomplex_t IT_0181 = IT_0134*IT_0140;
    const ccomplex_t IT_0182 = IT_0143*IT_0145;
    const ccomplex_t IT_0183 = IT_0142*IT_0146;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*(IT_0180 + -IT_0181 + 
      -IT_0182 + IT_0183);
    const ccomplex_t IT_0185 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0186 = -IT_0006;
    const ccomplex_t IT_0187 = IT_0151 + IT_0186;
    const ccomplex_t IT_0188 = -IT_0004;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*m_W*e_em*(IT_0139*IT_0153
      *IT_0187 + IT_0188*(IT_0135*IT_0154*IT_0155 + -(IT_0024*IT_0139*IT_0154 + 
      -IT_0004*IT_0135*IT_0187)*IT_0188));
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = IT_0185*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = IT_0184*IT_0192;
    const ccomplex_t IT_0194 = cpow((-2)*s_23 + IT_0000 + IT_0002 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0195 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0196 = IT_0007*IT_0195;
    const ccomplex_t IT_0197 = 2*IT_0196;
    const ccomplex_t IT_0198 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0199 = IT_0013*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0202 = IT_0007*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = IT_0200 + IT_0203;
    const ccomplex_t IT_0205 = IT_0197 + IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = 0.5*IT_0206;
    const ccomplex_t IT_0208 = IT_0101*IT_0207;
    const ccomplex_t IT_0209 = IT_0059*IT_0194*IT_0208;
    const ccomplex_t IT_0210 = IT_0054*IT_0174;
    const ccomplex_t IT_0211 = IT_0001*IT_0176*IT_0210;
    const ccomplex_t IT_0212 = IT_0073*IT_0207;
    const ccomplex_t IT_0213 = IT_0194*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = m_N_1*IT_0214;
    const ccomplex_t IT_0216 = IT_0133 + -IT_0161 + IT_0179 + -IT_0193 + 
      -IT_0209 + -IT_0211 + IT_0215;
    const ccomplex_t IT_0217 = 2*IT_0216;
    const ccomplex_t IT_0218 = IT_0001*IT_0003*IT_0210;
    const ccomplex_t IT_0219 = IT_0060*IT_0212;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = m_N_1*IT_0220;
    const ccomplex_t IT_0222 = IT_0059*IT_0060*IT_0208;
    const ccomplex_t IT_0223 = IT_0003*IT_0175;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = m_N_1*IT_0224;
    const ccomplex_t IT_0226 = IT_0218 + -IT_0221 + IT_0222 + -IT_0225;
    const ccomplex_t IT_0227 = (-2)*IT_0226;
    const ccomplex_t IT_0228 = IT_0217 + IT_0227;
    const ccomplex_t IT_0229 = 2*conj(IT_0216);
    const ccomplex_t IT_0230 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0231 = IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0233 = IT_0136*IT_0232;
    const ccomplex_t IT_0234 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0235 = IT_0140*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0237 = IT_0143*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0239 = IT_0146*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0233 + IT_0235 + 
      -IT_0237 + -IT_0239);
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = IT_0160*IT_0241;
    const ccomplex_t IT_0243 = IT_0059*IT_0087*IT_0194;
    const ccomplex_t IT_0244 = IT_0136*IT_0234;
    const ccomplex_t IT_0245 = IT_0140*IT_0232;
    const ccomplex_t IT_0246 = IT_0143*IT_0238;
    const ccomplex_t IT_0247 = IT_0146*IT_0236;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0244 + -IT_0245 + 
      -IT_0246 + IT_0247);
    const ccomplex_t IT_0249 = IT_0192*IT_0248;
    const ccomplex_t IT_0250 = IT_0055*IT_0176;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = m_N_1*IT_0251;
    const ccomplex_t IT_0253 = IT_0001*IT_0040*IT_0176;
    const ccomplex_t IT_0254 = IT_0102*IT_0194;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = m_N_1*IT_0255;
    const ccomplex_t IT_0257 = -IT_0133 + -IT_0242 + -IT_0243 + -IT_0249 +
       IT_0252 + -IT_0253 + IT_0256;
    const ccomplex_t IT_0258 = 2*IT_0226;
    const ccomplex_t IT_0259 = (-2)*IT_0216;
    const ccomplex_t IT_0260 = IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = 2*conj(IT_0226);
    const ccomplex_t IT_0262 = (-2)*conj(IT_0216);
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 2*s_12;
    const ccomplex_t IT_0265 = (-2)*s_12;
    const ccomplex_t IT_0266 = s_23*m_N_1;
    const ccomplex_t IT_0267 = IT_0108*IT_0120;
    const ccomplex_t IT_0268 = -IT_0251 + -IT_0255 + IT_0267;
    const ccomplex_t IT_0269 = IT_0220 + IT_0224;
    const ccomplex_t IT_0270 = conj(IT_0268) + conj(IT_0269);
    const ccomplex_t IT_0271 = IT_0268 + IT_0269;
    const ccomplex_t IT_0272 = IT_0108*IT_0129;
    const ccomplex_t IT_0273 = -IT_0178 + -IT_0214 + IT_0272;
    const ccomplex_t IT_0274 = IT_0057 + IT_0104;
    const ccomplex_t IT_0275 = conj(IT_0273) + conj(IT_0274);
    const ccomplex_t IT_0276 = 2*IT_0106;
    const ccomplex_t IT_0277 = (-2)*IT_0257;
    const ccomplex_t IT_0278 = IT_0273 + IT_0274;
    const ccomplex_t IT_0279 = 2*conj(IT_0106);
    const ccomplex_t IT_0280 = (-2)*conj(IT_0257);
    const ccomplex_t IT_0281 = s_13*m_N_1;
    const ccomplex_t IT_0282 = 2*IT_0257;
    const ccomplex_t IT_0283 = (-2)*IT_0106;
    const ccomplex_t IT_0284 = 2*conj(IT_0257);
    const ccomplex_t IT_0285 = (-2)*conj(IT_0106);
    const ccomplex_t IT_0286 = s_12*IT_0002;
    const ccomplex_t IT_0287 = s_13*s_23;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = IT_0286 + IT_0288;
    const ccomplex_t IT_0290 = (-2)*IT_0268;
    const ccomplex_t IT_0291 = (-2)*conj(IT_0268);
    const ccomplex_t IT_0292 = (-2)*IT_0273;
    const ccomplex_t IT_0293 = (-2)*conj(IT_0273);
    const ccomplex_t IT_0294 = (-2)*IT_0274;
    const ccomplex_t IT_0295 = IT_0000*IT_0002;
    const ccomplex_t IT_0296 = IT_0000*(conj(IT_0106)*IT_0228 + IT_0106
      *IT_0231 + conj(IT_0257)*IT_0260 + IT_0257*IT_0263) + (IT_0106*conj
      (IT_0106) + IT_0216*conj(IT_0216) + IT_0226*conj(IT_0226) + IT_0257*conj
      (IT_0257))*IT_0264 + (conj(IT_0216)*IT_0226 + IT_0216*conj(IT_0226) + conj
      (IT_0106)*IT_0257 + IT_0106*conj(IT_0257))*IT_0265 + IT_0266*(IT_0260
      *IT_0270 + IT_0263*IT_0271 + IT_0275*(IT_0276 + IT_0277) + IT_0278*
      (IT_0279 + IT_0280)) + IT_0281*(IT_0228*IT_0275 + IT_0231*IT_0278 +
       IT_0270*(IT_0282 + IT_0283) + IT_0271*(IT_0284 + IT_0285)) + IT_0289*
      (IT_0270*IT_0290 + (-2)*IT_0269*(conj(IT_0269) + (-0.5)*IT_0291) + IT_0275
      *IT_0292 + IT_0274*IT_0293 + conj(IT_0274)*IT_0294) + (conj(IT_0273)
      *IT_0290 + (-2)*conj(IT_0274)*(IT_0269 + (-0.5)*IT_0290) + IT_0278*IT_0291
       + IT_0269*IT_0293 + conj(IT_0269)*(IT_0292 + IT_0294))*IT_0295;
    return create_ccomplex_return(IT_0296);
}

