#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_tau_tau(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = m_N_1*m_tau;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 4*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0007, -1);
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0014;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0010 + IT_0013 + 
      -IT_0015 + -IT_0016);
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = e_em*IT_0008;
    const ccomplex_t IT_0020 = e_em*IT_0012;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0019 + -IT_0020);
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = pow(m_N_1, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0018*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0007;
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0018*IT_0031;
    const ccomplex_t IT_0033 = 0.5*IT_0032;
    const ccomplex_t IT_0034 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0035 = cos(alpha);
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = IT_0034*IT_0036;
    const ccomplex_t IT_0038 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0039 = sin(alpha);
    const ccomplex_t IT_0040 = IT_0006*IT_0039;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0043 = IT_0011*IT_0035;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0046 = IT_0011*IT_0039;
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0037 + IT_0041 + 
      -IT_0044 + -IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = pow(m_W, -1);
    const ccomplex_t IT_0051 = cos(beta);
    const ccomplex_t IT_0052 = cpow(IT_0051, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*m_tau*IT_0011*IT_0039
      *IT_0050*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0049*IT_0057;
    const ccomplex_t IT_0059 = IT_0036*IT_0038;
    const ccomplex_t IT_0060 = IT_0034*IT_0040;
    const ccomplex_t IT_0061 = IT_0043*IT_0045;
    const ccomplex_t IT_0062 = IT_0042*IT_0046;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + -IT_0060 + 
      -IT_0061 + IT_0062);
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*e_em*m_tau*IT_0011*IT_0035
      *IT_0050*IT_0052;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0067 = IT_0065*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = IT_0063*IT_0068;
    const ccomplex_t IT_0070 = (-0.5)*IT_0058 + (-0.5)*IT_0069;
    const ccomplex_t IT_0071 = 0.5*IT_0058 + 0.5*IT_0069;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0070 + IT_0072;
    const ccomplex_t IT_0074 = -conj(IT_0071);
    const ccomplex_t IT_0075 = conj(IT_0070) + IT_0074;
    const ccomplex_t IT_0076 = IT_0028 + IT_0033;
    const ccomplex_t IT_0077 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0078 = IT_0036*IT_0077;
    const ccomplex_t IT_0079 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0080 = IT_0040*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0082 = IT_0043*IT_0081;
    const ccomplex_t IT_0083 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0084 = IT_0046*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0078 + IT_0080 + 
      -IT_0082 + -IT_0084);
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = IT_0057*IT_0086;
    const ccomplex_t IT_0088 = IT_0036*IT_0079;
    const ccomplex_t IT_0089 = IT_0040*IT_0077;
    const ccomplex_t IT_0090 = IT_0043*IT_0083;
    const ccomplex_t IT_0091 = IT_0046*IT_0081;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0088 + -IT_0089 + 
      -IT_0090 + IT_0091);
    const ccomplex_t IT_0093 = IT_0068*IT_0092;
    const ccomplex_t IT_0094 = (-0.5)*IT_0087 + (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = 0.5*IT_0087 + 0.5*IT_0093;
    const ccomplex_t IT_0096 = -conj(IT_0095);
    const ccomplex_t IT_0097 = conj(IT_0094) + IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0017;
    const ccomplex_t IT_0099 = IT_0031*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0026*IT_0098;
    const ccomplex_t IT_0102 = 0.5*IT_0101;
    const ccomplex_t IT_0103 = IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = -IT_0095;
    const ccomplex_t IT_0105 = IT_0094 + IT_0104;
    const ccomplex_t IT_0106 = conj(IT_0100) + conj(IT_0102);
    const ccomplex_t IT_0107 = -s_24;
    const ccomplex_t IT_0108 = s_23 + IT_0107;
    const ccomplex_t IT_0109 = IT_0000*IT_0108;
    const ccomplex_t IT_0110 = 4*IT_0109;
    const ccomplex_t IT_0111 = (-0.5)*IT_0099;
    const ccomplex_t IT_0112 = (-0.5)*IT_0101;
    const ccomplex_t IT_0113 = -conj(IT_0094);
    const ccomplex_t IT_0114 = conj(IT_0095) + IT_0113;
    const ccomplex_t IT_0115 = -IT_0094;
    const ccomplex_t IT_0116 = IT_0095 + IT_0115;
    const ccomplex_t IT_0117 = conj(IT_0111) + conj(IT_0112);
    const ccomplex_t IT_0118 = -conj(IT_0070);
    const ccomplex_t IT_0119 = conj(IT_0071) + IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0032;
    const ccomplex_t IT_0121 = (-0.5)*IT_0027;
    const ccomplex_t IT_0122 = IT_0120 + IT_0121;
    const ccomplex_t IT_0123 = -IT_0070;
    const ccomplex_t IT_0124 = IT_0071 + IT_0123;
    const ccomplex_t IT_0125 = conj(IT_0120) + conj(IT_0121);
    const ccomplex_t IT_0126 = 2*IT_0109;
    const ccomplex_t IT_0127 = pow(m_tau, 2);
    const ccomplex_t IT_0128 = cpow(s_23 + (-0.5)*IT_0023 + (-0.5)*IT_0127 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0129 = conj(N_B1)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0130 = IT_0006*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = IT_0011*IT_0052;
    const ccomplex_t IT_0133 = conj(N_d1)*e_em*m_tau*IT_0050*conj(U_stau_01);
    const ccomplex_t IT_0134 = IT_0132*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0131 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0139 = IT_0006*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = N_d1*e_em*m_tau*IT_0050*U_stau_01;
    const ccomplex_t IT_0142 = IT_0132*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0140 + 0.5*IT_0143);
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = IT_0137*IT_0145;
    const ccomplex_t IT_0147 = IT_0128*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = (-0.5)*IT_0148;
    const ccomplex_t IT_0150 = cpow(s_23 + (-0.5)*IT_0023 + (-0.5)*IT_0127 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0151 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0152 = IT_0006*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = N_d1*e_em*m_tau*IT_0050*U_stau_00;
    const ccomplex_t IT_0155 = IT_0132*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + 0.5*IT_0156);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0160 = IT_0006*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = conj(N_d1)*e_em*m_tau*IT_0050*conj(U_stau_00);
    const ccomplex_t IT_0163 = IT_0132*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0161 + 0.5*IT_0164);
    const ccomplex_t IT_0166 = -IT_0165;
    const ccomplex_t IT_0167 = IT_0158*IT_0166;
    const ccomplex_t IT_0168 = IT_0150*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = (-0.5)*IT_0169;
    const ccomplex_t IT_0171 = IT_0149 + IT_0170;
    const ccomplex_t IT_0172 = conj(IT_0149) + conj(IT_0170);
    const ccomplex_t IT_0173 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0174 = IT_0006*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0177 = IT_0011*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = conj(N_d1)*e_em*m_tau*IT_0050*U_stau_11;
    const ccomplex_t IT_0180 = IT_0132*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*(IT_0175 + IT_0178 + 
      -IT_0181);
    const ccomplex_t IT_0183 = (-0.5)*IT_0182;
    const ccomplex_t IT_0184 = N_B1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0185 = IT_0006*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = N_W1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0188 = IT_0011*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = N_d1*e_em*m_tau*IT_0050*conj(U_stau_11);
    const ccomplex_t IT_0191 = IT_0132*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*(IT_0186 + IT_0189 + 
      -IT_0192);
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = IT_0183*IT_0194;
    const ccomplex_t IT_0196 = IT_0128*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*IT_0196;
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0200 = IT_0006*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0203 = IT_0011*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = conj(N_d1)*e_em*m_tau*IT_0050*U_stau_10;
    const ccomplex_t IT_0206 = IT_0132*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0201 + IT_0204 + 
      -IT_0207);
    const ccomplex_t IT_0209 = (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0211 = IT_0006*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0214 = IT_0011*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = N_d1*e_em*m_tau*IT_0050*conj(U_stau_10);
    const ccomplex_t IT_0217 = IT_0132*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0212 + IT_0215 + 
      -IT_0218);
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = IT_0209*IT_0220;
    const ccomplex_t IT_0222 = IT_0150*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = IT_0198 + IT_0224;
    const ccomplex_t IT_0226 = conj(IT_0198) + conj(IT_0224);
    const ccomplex_t IT_0227 = (-2)*IT_0109;
    const ccomplex_t IT_0228 = cpow(s_13 + (-0.5)*IT_0023 + (-0.5)*IT_0127 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0229 = IT_0221*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = cpow(s_13 + (-0.5)*IT_0023 + (-0.5)*IT_0127 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0232 = IT_0195*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = 0.5*IT_0230 + 0.5*IT_0233;
    const ccomplex_t IT_0235 = IT_0146*IT_0231;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = IT_0167*IT_0228;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = 0.5*IT_0239;
    const ccomplex_t IT_0241 = IT_0237 + IT_0240;
    const ccomplex_t IT_0242 = conj(IT_0237) + conj(IT_0240);
    const ccomplex_t IT_0243 = 2*IT_0003;
    const ccomplex_t IT_0244 = s_12*IT_0127;
    const ccomplex_t IT_0245 = s_14*s_23;
    const ccomplex_t IT_0246 = s_13*s_24;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = IT_0245 + IT_0247;
    const ccomplex_t IT_0249 = s_12*s_34;
    const ccomplex_t IT_0250 = -IT_0249;
    const ccomplex_t IT_0251 = IT_0248 + IT_0250;
    const ccomplex_t IT_0252 = IT_0244 + IT_0251;
    const ccomplex_t IT_0253 = IT_0158*IT_0220;
    const ccomplex_t IT_0254 = IT_0228*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = IT_0145*IT_0194;
    const ccomplex_t IT_0257 = IT_0231*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = (-0.5)*IT_0255 + (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = 2*IT_0071;
    const ccomplex_t IT_0261 = (-2)*IT_0070;
    const ccomplex_t IT_0262 = IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = 2*conj(IT_0071);
    const ccomplex_t IT_0264 = (-2)*conj(IT_0070);
    const ccomplex_t IT_0265 = IT_0263 + IT_0264;
    const ccomplex_t IT_0266 = IT_0166*IT_0209;
    const ccomplex_t IT_0267 = IT_0228*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = IT_0137*IT_0183;
    const ccomplex_t IT_0270 = IT_0231*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = (-0.5)*IT_0268 + (-0.5)*IT_0271;
    const ccomplex_t IT_0273 = 2*IT_0095;
    const ccomplex_t IT_0274 = (-2)*IT_0094;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 2*conj(IT_0095);
    const ccomplex_t IT_0277 = (-2)*conj(IT_0094);
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = s_23 + s_24;
    const ccomplex_t IT_0280 = IT_0000*IT_0279;
    const ccomplex_t IT_0281 = IT_0006*IT_0051;
    const ccomplex_t IT_0282 = IT_0077*IT_0281;
    const ccomplex_t IT_0283 = sin(beta);
    const ccomplex_t IT_0284 = IT_0006*IT_0283;
    const ccomplex_t IT_0285 = IT_0079*IT_0284;
    const ccomplex_t IT_0286 = IT_0011*IT_0051;
    const ccomplex_t IT_0287 = IT_0081*IT_0286;
    const ccomplex_t IT_0288 = IT_0011*IT_0283;
    const ccomplex_t IT_0289 = IT_0083*IT_0288;
    const ccomplex_t IT_0290 = IT_0282 + -IT_0285 + -IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = e_em*m_tau*IT_0011*IT_0050*IT_0052*IT_0283;
    const ccomplex_t IT_0292 = (-0.5)*IT_0291;
    const ccomplex_t IT_0293 = cpow((-2)*s_12 + (-2)*IT_0023 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0294 = IT_0292*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*IT_0294;
    const ccomplex_t IT_0296 = IT_0290*IT_0295;
    const ccomplex_t IT_0297 = IT_0079*IT_0281;
    const ccomplex_t IT_0298 = IT_0077*IT_0284;
    const ccomplex_t IT_0299 = IT_0083*IT_0286;
    const ccomplex_t IT_0300 = IT_0081*IT_0288;
    const ccomplex_t IT_0301 = IT_0297 + IT_0298 + -IT_0299 + -IT_0300;
    const ccomplex_t IT_0302 = e_em*m_tau*IT_0011*IT_0050;
    const ccomplex_t IT_0303 = 0.5*IT_0302;
    const ccomplex_t IT_0304 = pow(m_Z, 2);
    const ccomplex_t IT_0305 = cpow((-2)*s_12 + (-2)*IT_0023 + IT_0304 + 
      -reg_prop, -1);
    const ccomplex_t IT_0306 = IT_0303*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*IT_0306;
    const ccomplex_t IT_0308 = IT_0301*IT_0307;
    const ccomplex_t IT_0309 = 0.5*IT_0296 + 0.5*IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0149;
    const ccomplex_t IT_0311 = 2*IT_0170;
    const ccomplex_t IT_0312 = IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = (-2)*IT_0234;
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = 2*conj(IT_0149);
    const ccomplex_t IT_0316 = 2*conj(IT_0170);
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = (-2)*conj(IT_0234);
    const ccomplex_t IT_0319 = IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = IT_0034*IT_0281;
    const ccomplex_t IT_0321 = IT_0038*IT_0284;
    const ccomplex_t IT_0322 = IT_0042*IT_0286;
    const ccomplex_t IT_0323 = IT_0045*IT_0288;
    const ccomplex_t IT_0324 = -IT_0320 + IT_0321 + IT_0322 + -IT_0323;
    const ccomplex_t IT_0325 = IT_0295*IT_0324;
    const ccomplex_t IT_0326 = IT_0038*IT_0281;
    const ccomplex_t IT_0327 = IT_0034*IT_0284;
    const ccomplex_t IT_0328 = IT_0045*IT_0286;
    const ccomplex_t IT_0329 = IT_0042*IT_0288;
    const ccomplex_t IT_0330 = -IT_0326 + -IT_0327 + IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = IT_0307*IT_0330;
    const ccomplex_t IT_0332 = (-0.5)*IT_0325 + (-0.5)*IT_0331;
    const ccomplex_t IT_0333 = 2*IT_0198;
    const ccomplex_t IT_0334 = 2*IT_0224;
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = (-2)*IT_0237;
    const ccomplex_t IT_0337 = (-2)*IT_0240;
    const ccomplex_t IT_0338 = IT_0335 + IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = 2*conj(IT_0198);
    const ccomplex_t IT_0340 = 2*conj(IT_0224);
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = (-2)*conj(IT_0237);
    const ccomplex_t IT_0343 = (-2)*conj(IT_0240);
    const ccomplex_t IT_0344 = IT_0341 + IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = 0.5*IT_0325 + 0.5*IT_0331;
    const ccomplex_t IT_0346 = 2*IT_0237;
    const ccomplex_t IT_0347 = 2*IT_0240;
    const ccomplex_t IT_0348 = IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = (-2)*IT_0198;
    const ccomplex_t IT_0350 = (-2)*IT_0224;
    const ccomplex_t IT_0351 = IT_0348 + IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 2*conj(IT_0237);
    const ccomplex_t IT_0353 = 2*conj(IT_0240);
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = (-2)*conj(IT_0198);
    const ccomplex_t IT_0356 = (-2)*conj(IT_0224);
    const ccomplex_t IT_0357 = IT_0354 + IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = (-0.5)*IT_0296 + (-0.5)*IT_0308;
    const ccomplex_t IT_0359 = 2*IT_0234;
    const ccomplex_t IT_0360 = (-2)*IT_0149;
    const ccomplex_t IT_0361 = (-2)*IT_0170;
    const ccomplex_t IT_0362 = IT_0359 + IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = 2*conj(IT_0234);
    const ccomplex_t IT_0364 = (-2)*conj(IT_0149);
    const ccomplex_t IT_0365 = (-2)*conj(IT_0170);
    const ccomplex_t IT_0366 = IT_0363 + IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = conj(IT_0100) + conj(IT_0120);
    const ccomplex_t IT_0368 = 4*IT_0345;
    const ccomplex_t IT_0369 = (-4)*IT_0332;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = IT_0100 + IT_0120;
    const ccomplex_t IT_0372 = 4*conj(IT_0345);
    const ccomplex_t IT_0373 = (-4)*conj(IT_0332);
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = conj(IT_0102) + conj(IT_0121);
    const ccomplex_t IT_0376 = 4*IT_0332;
    const ccomplex_t IT_0377 = (-4)*IT_0345;
    const ccomplex_t IT_0378 = IT_0376 + IT_0377;
    const ccomplex_t IT_0379 = IT_0102 + IT_0121;
    const ccomplex_t IT_0380 = 4*conj(IT_0332);
    const ccomplex_t IT_0381 = (-4)*conj(IT_0345);
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = conj(IT_0033) + conj(IT_0111);
    const ccomplex_t IT_0384 = 4*IT_0309;
    const ccomplex_t IT_0385 = (-4)*IT_0358;
    const ccomplex_t IT_0386 = IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = IT_0033 + IT_0111;
    const ccomplex_t IT_0388 = 4*conj(IT_0309);
    const ccomplex_t IT_0389 = (-4)*conj(IT_0358);
    const ccomplex_t IT_0390 = IT_0388 + IT_0389;
    const ccomplex_t IT_0391 = conj(IT_0028) + conj(IT_0112);
    const ccomplex_t IT_0392 = 4*IT_0358;
    const ccomplex_t IT_0393 = (-4)*IT_0309;
    const ccomplex_t IT_0394 = IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = IT_0028 + IT_0112;
    const ccomplex_t IT_0396 = 4*conj(IT_0358);
    const ccomplex_t IT_0397 = (-4)*conj(IT_0309);
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = conj(IT_0309)*IT_0314 + IT_0309*IT_0319 + conj
      (IT_0332)*IT_0338 + IT_0332*IT_0344 + conj(IT_0345)*IT_0351 + IT_0345
      *IT_0357 + conj(IT_0358)*IT_0362 + IT_0358*IT_0366 + IT_0367*IT_0370 +
       IT_0371*IT_0374 + IT_0375*IT_0378 + IT_0379*IT_0382 + IT_0383*IT_0386 +
       IT_0387*IT_0390 + IT_0391*IT_0394 + IT_0395*IT_0398;
    const ccomplex_t IT_0400 = s_13 + s_14;
    const ccomplex_t IT_0401 = IT_0000*IT_0400;
    const ccomplex_t IT_0402 = IT_0319*IT_0332 + IT_0314*conj(IT_0332) + conj
      (IT_0309)*IT_0338 + IT_0309*IT_0344 + IT_0357*IT_0358 + IT_0351*conj
      (IT_0358) + conj(IT_0345)*IT_0362 + IT_0345*IT_0366 + IT_0378*IT_0383 +
       IT_0375*IT_0386 + IT_0382*IT_0387 + IT_0379*IT_0390 + IT_0370*IT_0391 +
       IT_0367*IT_0394 + IT_0374*IT_0395 + IT_0371*IT_0398;
    const ccomplex_t IT_0403 = IT_0128*IT_0256;
    const ccomplex_t IT_0404 = (0 + _Complex_I*1)*IT_0403;
    const ccomplex_t IT_0405 = IT_0150*IT_0253;
    const ccomplex_t IT_0406 = (0 + _Complex_I*1)*IT_0405;
    const ccomplex_t IT_0407 = 0.5*IT_0404 + 0.5*IT_0406;
    const ccomplex_t IT_0408 = IT_0128*IT_0269;
    const ccomplex_t IT_0409 = (0 + _Complex_I*1)*IT_0408;
    const ccomplex_t IT_0410 = IT_0150*IT_0266;
    const ccomplex_t IT_0411 = (0 + _Complex_I*1)*IT_0410;
    const ccomplex_t IT_0412 = 0.5*IT_0409 + 0.5*IT_0411;
    const ccomplex_t IT_0413 = -IT_0244;
    const ccomplex_t IT_0414 = IT_0245 + IT_0246 + IT_0250;
    const ccomplex_t IT_0415 = 2*IT_0259;
    const ccomplex_t IT_0416 = 2*conj(IT_0259);
    const ccomplex_t IT_0417 = 2*IT_0272;
    const ccomplex_t IT_0418 = 2*conj(IT_0272);
    const ccomplex_t IT_0419 = s_24*m_N_1*m_tau;
    const ccomplex_t IT_0420 = 4*IT_0234;
    const ccomplex_t IT_0421 = 0.125*IT_0420;
    const ccomplex_t IT_0422 = IT_0395 + IT_0421;
    const ccomplex_t IT_0423 = 4*IT_0407;
    const ccomplex_t IT_0424 = 0.125*IT_0423;
    const ccomplex_t IT_0425 = IT_0259 + IT_0424;
    const ccomplex_t IT_0426 = 4*conj(IT_0407);
    const ccomplex_t IT_0427 = 0.125*IT_0426;
    const ccomplex_t IT_0428 = conj(IT_0259) + IT_0427;
    const ccomplex_t IT_0429 = 2*IT_0407;
    const ccomplex_t IT_0430 = 4*IT_0259;
    const ccomplex_t IT_0431 = IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 2*conj(IT_0407);
    const ccomplex_t IT_0433 = 4*conj(IT_0259);
    const ccomplex_t IT_0434 = IT_0432 + IT_0433;
    const ccomplex_t IT_0435 = 4*IT_0412;
    const ccomplex_t IT_0436 = 4*conj(IT_0412);
    const ccomplex_t IT_0437 = 4*conj(IT_0234);
    const ccomplex_t IT_0438 = 0.125*IT_0363*IT_0412 + 0.125*IT_0359*conj
      (IT_0412) + conj(IT_0272)*IT_0422 + IT_0367*IT_0425 + IT_0371*IT_0428 +
       0.125*IT_0242*IT_0431 + 0.125*IT_0241*IT_0434 + 0.125*IT_0391*IT_0435 +
       0.125*IT_0395*IT_0436 + IT_0272*(IT_0391 + 0.125*IT_0437);
    const ccomplex_t IT_0439 = -IT_0127;
    const ccomplex_t IT_0440 = s_34 + IT_0439;
    const ccomplex_t IT_0441 = s_12*IT_0440;
    const ccomplex_t IT_0442 = (-8)*conj(IT_0071);
    const ccomplex_t IT_0443 = IT_0251 + IT_0413;
    const ccomplex_t IT_0444 = (-2)*IT_0259;
    const ccomplex_t IT_0445 = (-2)*conj(IT_0259);
    const ccomplex_t IT_0446 = (-2)*IT_0272;
    const ccomplex_t IT_0447 = (-2)*conj(IT_0272);
    const ccomplex_t IT_0448 = IT_0244 + IT_0248 + IT_0249;
    const ccomplex_t IT_0449 = 2*IT_0412;
    const ccomplex_t IT_0450 = 2*conj(IT_0412);
    const ccomplex_t IT_0451 = (-2)*IT_0407;
    const ccomplex_t IT_0452 = (-2)*conj(IT_0407);
    const ccomplex_t IT_0453 = (-2)*IT_0412;
    const ccomplex_t IT_0454 = (-2)*conj(IT_0412);
    const ccomplex_t IT_0455 = (-4)*IT_0109;
    const ccomplex_t IT_0456 = (-4)*IT_0003;
    const ccomplex_t IT_0457 = (-2)*IT_0003;
    const ccomplex_t IT_0458 = s_23*m_N_1*m_tau;
    const ccomplex_t IT_0459 = 4*conj(IT_0170);
    const ccomplex_t IT_0460 = 0.125*IT_0459;
    const ccomplex_t IT_0461 = 4*conj(IT_0149);
    const ccomplex_t IT_0462 = 0.125*IT_0461;
    const ccomplex_t IT_0463 = IT_0383 + IT_0460 + IT_0462;
    const ccomplex_t IT_0464 = 0.125*IT_0433;
    const ccomplex_t IT_0465 = conj(IT_0407) + IT_0464;
    const ccomplex_t IT_0466 = 0.125*IT_0430;
    const ccomplex_t IT_0467 = IT_0407 + IT_0466;
    const ccomplex_t IT_0468 = IT_0415 + IT_0423;
    const ccomplex_t IT_0469 = IT_0416 + IT_0426;
    const ccomplex_t IT_0470 = 4*IT_0272;
    const ccomplex_t IT_0471 = 4*conj(IT_0272);
    const ccomplex_t IT_0472 = 4*IT_0170;
    const ccomplex_t IT_0473 = 0.125*IT_0472;
    const ccomplex_t IT_0474 = 4*IT_0149;
    const ccomplex_t IT_0475 = 0.125*IT_0474;
    const ccomplex_t IT_0476 = 0.125*conj(IT_0272)*IT_0312 + 0.125*IT_0272
      *IT_0317 + IT_0412*IT_0463 + IT_0379*IT_0465 + IT_0375*IT_0467 + 0.125
      *IT_0226*IT_0468 + 0.125*IT_0225*IT_0469 + 0.125*IT_0383*IT_0470 + 0.125
      *IT_0387*IT_0471 + conj(IT_0412)*(IT_0387 + IT_0473 + IT_0475);
    const ccomplex_t IT_0477 = s_34*IT_0023;
    const ccomplex_t IT_0478 = 4*conj(IT_0240);
    const ccomplex_t IT_0479 = 0.125*IT_0478;
    const ccomplex_t IT_0480 = 4*conj(IT_0237);
    const ccomplex_t IT_0481 = 0.125*IT_0480;
    const ccomplex_t IT_0482 = IT_0367 + IT_0479 + IT_0481;
    const ccomplex_t IT_0483 = 4*IT_0198;
    const ccomplex_t IT_0484 = 0.125*IT_0483;
    const ccomplex_t IT_0485 = IT_0379 + IT_0484;
    const ccomplex_t IT_0486 = 4*IT_0224;
    const ccomplex_t IT_0487 = 0.125*IT_0486;
    const ccomplex_t IT_0488 = IT_0485 + IT_0487;
    const ccomplex_t IT_0489 = 4*conj(IT_0224);
    const ccomplex_t IT_0490 = 0.125*IT_0489;
    const ccomplex_t IT_0491 = 4*conj(IT_0198);
    const ccomplex_t IT_0492 = 0.125*IT_0491;
    const ccomplex_t IT_0493 = IT_0375 + IT_0490 + IT_0492;
    const ccomplex_t IT_0494 = 4*IT_0237;
    const ccomplex_t IT_0495 = 0.125*IT_0494;
    const ccomplex_t IT_0496 = IT_0371 + IT_0495;
    const ccomplex_t IT_0497 = 4*IT_0240;
    const ccomplex_t IT_0498 = 0.125*IT_0497;
    const ccomplex_t IT_0499 = IT_0496 + IT_0498;
    const ccomplex_t IT_0500 = 8*IT_0383;
    const ccomplex_t IT_0501 = 0.125*conj(IT_0240)*IT_0310 + 0.125*IT_0240
      *IT_0315 + 0.125*conj(IT_0234)*IT_0333 + 0.125*IT_0234*IT_0339 + 0.125
      *IT_0172*IT_0346 + 0.125*conj(IT_0170)*IT_0347 + 0.125*IT_0171*IT_0352 +
       0.125*IT_0170*IT_0353 + 0.125*conj(IT_0224)*IT_0359 + 0.125*IT_0224
      *IT_0363 + 0.125*IT_0375*IT_0420 + 0.125*IT_0379*IT_0437 + 0.125*IT_0371*
      (IT_0459 + IT_0461) + 0.125*IT_0367*(IT_0472 + IT_0474) + IT_0387*IT_0482 
      + IT_0391*IT_0488 + IT_0395*IT_0493 + 0.125*IT_0499*IT_0500;
    const ccomplex_t IT_0502 = 16*IT_0102;
    const ccomplex_t IT_0503 = 0.125*IT_0502;
    const ccomplex_t IT_0504 = IT_0225 + IT_0503;
    const ccomplex_t IT_0505 = 16*conj(IT_0102);
    const ccomplex_t IT_0506 = 16*IT_0121;
    const ccomplex_t IT_0507 = 16*IT_0111;
    const ccomplex_t IT_0508 = 8*conj(IT_0149);
    const ccomplex_t IT_0509 = 16*conj(IT_0111);
    const ccomplex_t IT_0510 = 8*conj(IT_0170);
    const ccomplex_t IT_0511 = 0.125*conj(IT_0407)*IT_0423 + 0.125*conj
      (IT_0412)*IT_0435 + 0.125*IT_0170*IT_0461 + conj(IT_0198)*IT_0485 + conj
      (IT_0224)*IT_0488 + 0.125*IT_0224*IT_0491 + conj(IT_0102)*IT_0504 + 0.125
      *IT_0121*IT_0505 + conj(IT_0121)*(IT_0504 + 0.125*IT_0506) + 0.125*IT_0500
      *(IT_0171 + 0.125*IT_0507) + 0.125*(IT_0111 + IT_0475)*IT_0508 + 0.125
      *IT_0033*(16*conj(IT_0033) + IT_0508 + IT_0509 + IT_0510) + 0.125*(IT_0111
       + IT_0473 + IT_0475)*IT_0510;
    const ccomplex_t IT_0512 = s_34 + IT_0127;
    const ccomplex_t IT_0513 = s_12*IT_0512;
    const ccomplex_t IT_0514 = (-8)*IT_0345;
    const ccomplex_t IT_0515 = (-8)*conj(IT_0345);
    const ccomplex_t IT_0516 = -conj(IT_0358);
    const ccomplex_t IT_0517 = conj(IT_0309) + IT_0516;
    const ccomplex_t IT_0518 = -IT_0358;
    const ccomplex_t IT_0519 = IT_0023*IT_0512;
    const ccomplex_t IT_0520 = conj(IT_0358)*(IT_0415 + IT_0451) + IT_0358*
      (IT_0416 + IT_0452) + 8*conj(IT_0345)*(IT_0358 + 0.125*IT_0417 + 0.125
      *IT_0453) + 8*IT_0345*(conj(IT_0358) + 0.125*IT_0418 + 0.125*IT_0454) +
       conj(IT_0309)*(IT_0429 + IT_0444 + IT_0514) + IT_0309*(IT_0432 + IT_0445 
      + IT_0515) + 8*IT_0332*(0.125*IT_0447 + 0.125*IT_0450 + IT_0517) + 8*conj
      (IT_0332)*(IT_0309 + 0.125*IT_0446 + 0.125*IT_0449 + IT_0518);
    const ccomplex_t IT_0521 = IT_0023*IT_0440;
    const ccomplex_t IT_0522 = 0.125*IT_0272;
    const ccomplex_t IT_0523 = 0.125*conj(IT_0272);
    const ccomplex_t IT_0524 = 0.5*conj(IT_0259)*IT_0275 + 0.5*IT_0259*IT_0278
       + -IT_0114*IT_0407 + IT_0094*conj(IT_0407) + 4*conj(IT_0071)*(IT_0094 + (
      -0.25)*IT_0412) + 4*conj(IT_0070)*(IT_0116 + 1./4*IT_0412) + 4*(IT_0070 + 
      -IT_0071)*(IT_0114 + 1./4*conj(IT_0412)) + -IT_0095*(conj(IT_0407) + (-0.5
      )*IT_0442) + 4*IT_0265*IT_0522 + 4*IT_0262*IT_0523;
    const ccomplex_t IT_0525 = IT_0023*IT_0127;
    const ccomplex_t IT_0526 = 16*conj(IT_0028);
    const ccomplex_t IT_0527 = 16*conj(IT_0112);
    const ccomplex_t IT_0528 = 16*IT_0028;
    const ccomplex_t IT_0529 = 16*IT_0112;
    const ccomplex_t IT_0530 = 0.125*IT_0367;
    const ccomplex_t IT_0531 = 0.125*IT_0387;
    const ccomplex_t IT_0532 = 0.125*IT_0383;
    const ccomplex_t IT_0533 = 0.125*IT_0412;
    const ccomplex_t IT_0534 = 0.125*conj(IT_0412);
    const ccomplex_t IT_0535 = conj(IT_0102)*IT_0171 + IT_0102*IT_0172 +
       IT_0241*IT_0391 + IT_0242*IT_0395 + IT_0225*IT_0463 + IT_0226*(IT_0387 +
       IT_0473 + IT_0475) + IT_0234*IT_0482 + conj(IT_0234)*IT_0499 + 0.125*conj
      (IT_0033)*IT_0506 + 0.125*conj(IT_0121)*(16*IT_0033 + 8*IT_0149 + 8
      *IT_0170 + IT_0507) + IT_0121*(IT_0172 + 0.125*IT_0509) + IT_0434*IT_0522 
      + IT_0431*IT_0523 + 0.125*IT_0371*(IT_0526 + IT_0527) + (IT_0528 + IT_0529
      )*IT_0530 + IT_0505*IT_0531 + IT_0502*IT_0532 + IT_0469*IT_0533 + IT_0468
      *IT_0534;
    const ccomplex_t IT_0536 = s_14*m_N_1*m_tau;
    const ccomplex_t IT_0537 = 0.125*IT_0379;
    const ccomplex_t IT_0538 = 0.125*IT_0375;
    const ccomplex_t IT_0539 = IT_0387*IT_0465 + IT_0383*IT_0467 + 0.125
      *IT_0172*IT_0468 + 0.125*IT_0171*IT_0469 + conj(IT_0412)*IT_0488 + IT_0412
      *IT_0493 + IT_0341*IT_0522 + IT_0335*IT_0523 + IT_0471*IT_0537 + IT_0470
      *IT_0538;
    const ccomplex_t IT_0540 = s_13*m_N_1*m_tau;
    const ccomplex_t IT_0541 = 0.125*IT_0234;
    const ccomplex_t IT_0542 = 0.125*conj(IT_0234);
    const ccomplex_t IT_0543 = IT_0391*IT_0425 + IT_0395*IT_0428 + 0.125
      *IT_0371*IT_0436 + IT_0272*IT_0482 + conj(IT_0272)*IT_0499 + IT_0435
      *IT_0530 + IT_0354*IT_0533 + IT_0348*IT_0534 + IT_0434*IT_0541 + IT_0431
      *IT_0542;
    const ccomplex_t IT_0544 = 0.125*IT_0226*IT_0348 + 0.125*IT_0225*IT_0354 +
       IT_0395*IT_0463 + IT_0391*(IT_0387 + IT_0473 + IT_0475) + IT_0367*IT_0488
       + IT_0371*IT_0493 + IT_0437*IT_0531 + IT_0420*IT_0532 + (IT_0478 +
       IT_0480)*IT_0537 + (IT_0494 + IT_0497)*IT_0538 + IT_0317*IT_0541 +
       IT_0312*IT_0542;
    const ccomplex_t IT_0545 = 16*IT_0100;
    const ccomplex_t IT_0546 = 0.125*IT_0529;
    const ccomplex_t IT_0547 = conj(IT_0234)*IT_0422 + 0.125*conj(IT_0259)
      *IT_0430 + 0.125*conj(IT_0272)*IT_0470 + 0.125*IT_0240*IT_0480 + conj
      (IT_0237)*IT_0496 + conj(IT_0240)*IT_0499 + 0.125*IT_0028*IT_0527 + 0.125
      *IT_0367*(16*IT_0120 + 8*IT_0237 + 8*IT_0240 + IT_0545) + conj(IT_0112)*
      (IT_0234 + IT_0546) + conj(IT_0028)*(IT_0234 + 0.125*IT_0528 + IT_0546);
    const ccomplex_t IT_0548 = IT_0004*((conj(IT_0028) + conj(IT_0033))
      *IT_0073 + IT_0075*IT_0076 + IT_0097*IT_0103 + IT_0105*IT_0106) + IT_0110*
      ((IT_0111 + IT_0112)*IT_0114 + IT_0116*IT_0117 + IT_0119*IT_0122 + IT_0124
      *IT_0125) + IT_0126*(IT_0114*IT_0171 + IT_0116*IT_0172 + IT_0119*IT_0225 +
       IT_0124*IT_0226) + IT_0227*(IT_0097*IT_0234 + IT_0105*conj(IT_0234) +
       IT_0075*IT_0241 + IT_0073*IT_0242) + (IT_0075*IT_0234 + IT_0073*conj
      (IT_0234) + IT_0097*IT_0241 + IT_0105*IT_0242)*IT_0243 + IT_0252*(conj
      (IT_0259)*IT_0262 + IT_0259*IT_0265 + conj(IT_0272)*IT_0275 + IT_0272
      *IT_0278) + IT_0280*IT_0399 + IT_0401*IT_0402 + (IT_0265*IT_0407 + IT_0262
      *conj(IT_0407) + IT_0278*IT_0412 + IT_0275*conj(IT_0412))*(IT_0248 +
       IT_0249 + IT_0413) + IT_0414*(conj(IT_0407)*IT_0415 + IT_0407*IT_0416 +
       conj(IT_0412)*IT_0417 + IT_0412*IT_0418) + 8*IT_0419*IT_0438 + (-8)
      *IT_0441*(IT_0071*(conj(IT_0070) + -conj(IT_0071)) + -IT_0114*IT_0116 + 
      -IT_0070*(conj(IT_0070) + 0.125*IT_0442)) + IT_0443*(conj(IT_0332)*IT_0415
       + IT_0332*IT_0416 + conj(IT_0309)*IT_0417 + IT_0309*IT_0418 + conj
      (IT_0345)*IT_0444 + IT_0345*IT_0445 + conj(IT_0358)*IT_0446 + IT_0358
      *IT_0447) + IT_0448*(conj(IT_0332)*IT_0429 + IT_0332*IT_0432 + conj
      (IT_0309)*IT_0449 + IT_0309*IT_0450 + conj(IT_0345)*IT_0451 + IT_0345
      *IT_0452 + conj(IT_0358)*IT_0453 + IT_0358*IT_0454) + (IT_0076*(conj
      (IT_0094) + IT_0096) + IT_0075*IT_0103 + (conj(IT_0028) + conj(IT_0033))
      *IT_0105 + IT_0073*IT_0106)*IT_0455 + ((IT_0111 + IT_0112)*IT_0119 + (conj
      (IT_0095) + IT_0113)*IT_0122 + IT_0117*IT_0124 + IT_0116*IT_0125)*IT_0456 
      + (IT_0119*IT_0171 + IT_0124*IT_0172 + (conj(IT_0095) + IT_0113)*IT_0225 +
       IT_0116*IT_0226)*IT_0457 + 8*IT_0458*IT_0476 + 8*IT_0477*IT_0501 + 8
      *IT_0245*IT_0511 + 8*IT_0513*(IT_0345*conj(IT_0345) + 0.125*conj(IT_0332)
      *IT_0514 + IT_0332*(conj(IT_0332) + 0.125*IT_0515) + IT_0517*(IT_0309 +
       IT_0518)) + IT_0519*IT_0520 + 2*IT_0521*IT_0524 + 8*IT_0525*IT_0535 + 8
      *IT_0536*IT_0539 + 8*IT_0540*IT_0543 + 8*IT_0244*IT_0544 + 8*IT_0246
      *IT_0547;
    return create_ccomplex_return(IT_0548);
}

