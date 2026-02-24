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
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
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
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = (0 + _Complex_I*1)*IT_0012;
    const ccomplex_t IT_0014 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0015 = IT_0005*IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0014;
    const ccomplex_t IT_0017 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0018 = IT_0005*IT_0017;
    const ccomplex_t IT_0019 = IT_0008*IT_0017;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0015 + IT_0016 + 
      -IT_0018 + -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0021;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0025 = IT_0003*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = cos(beta);
    const ccomplex_t IT_0028 = cpow(IT_0027, -1);
    const ccomplex_t IT_0029 = IT_0007*IT_0028;
    const ccomplex_t IT_0030 = pow(m_W, -1);
    const ccomplex_t IT_0031 = conj(N_d1)*e_em*m_tau*IT_0030*conj(U_stau_00);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0026 + 0.5*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = N_d1*e_em*m_tau*IT_0030*U_stau_00;
    const ccomplex_t IT_0040 = IT_0029*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + 0.5*IT_0041);
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = IT_0035*IT_0043;
    const ccomplex_t IT_0045 = pow(m_tau, 2);
    const ccomplex_t IT_0046 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0045 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = 0.5*IT_0048;
    const ccomplex_t IT_0050 = s_12*IT_0045;
    const ccomplex_t IT_0051 = 4*IT_0050;
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = conj(N_B1)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0054 = IT_0003*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = conj(N_d1)*e_em*m_tau*IT_0030*conj(U_stau_01);
    const ccomplex_t IT_0057 = IT_0029*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0055 + 0.5*IT_0058);
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0062 = IT_0003*IT_0061;
    const ccomplex_t IT_0063 = 1.4142135623731*IT_0062;
    const ccomplex_t IT_0064 = N_d1*e_em*m_tau*IT_0030*U_stau_01;
    const ccomplex_t IT_0065 = IT_0029*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(IT_0063 + 0.5*IT_0066);
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = IT_0060*IT_0068;
    const ccomplex_t IT_0070 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0045 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0071 = IT_0069*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = IT_0051*IT_0073;
    const ccomplex_t IT_0075 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0076 = IT_0003*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0079 = IT_0007*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = conj(N_d1)*e_em*m_tau*IT_0030*U_stau_10;
    const ccomplex_t IT_0082 = IT_0029*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0077 + IT_0080 + 
      -IT_0083);
    const ccomplex_t IT_0085 = (-0.5)*IT_0084;
    const ccomplex_t IT_0086 = IT_0035*IT_0085;
    const ccomplex_t IT_0087 = IT_0046*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0090 = IT_0003*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0093 = IT_0007*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = conj(N_d1)*e_em*m_tau*IT_0030*U_stau_11;
    const ccomplex_t IT_0096 = IT_0029*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0091 + IT_0094 + 
      -IT_0097);
    const ccomplex_t IT_0099 = (-0.5)*IT_0098;
    const ccomplex_t IT_0100 = IT_0060*IT_0099;
    const ccomplex_t IT_0101 = IT_0070*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0088 + (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = s_14*m_N_1*m_tau;
    const ccomplex_t IT_0105 = 4*IT_0104;
    const ccomplex_t IT_0106 = IT_0103*IT_0105;
    const ccomplex_t IT_0107 = N_B1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0108 = IT_0003*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_W1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = N_d1*e_em*m_tau*IT_0030*conj(U_stau_11);
    const ccomplex_t IT_0114 = IT_0029*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0109 + IT_0112 + 
      -IT_0115);
    const ccomplex_t IT_0117 = (-0.5)*IT_0116;
    const ccomplex_t IT_0118 = IT_0068*IT_0117;
    const ccomplex_t IT_0119 = IT_0070*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0122 = IT_0003*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0125 = IT_0007*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = N_d1*e_em*m_tau*IT_0030*conj(U_stau_10);
    const ccomplex_t IT_0128 = IT_0029*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0123 + IT_0126 + 
      -IT_0129);
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0043*IT_0131;
    const ccomplex_t IT_0133 = IT_0046*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0120 + (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = s_23*m_N_1*m_tau;
    const ccomplex_t IT_0137 = 4*IT_0136;
    const ccomplex_t IT_0138 = IT_0135*IT_0137;
    const ccomplex_t IT_0139 = IT_0085*IT_0131;
    const ccomplex_t IT_0140 = IT_0046*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0099*IT_0117;
    const ccomplex_t IT_0143 = IT_0070*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = 0.5*IT_0141 + 0.5*IT_0144;
    const ccomplex_t IT_0146 = s_34*IT_0000;
    const ccomplex_t IT_0147 = 4*IT_0146;
    const ccomplex_t IT_0148 = IT_0145*IT_0147;
    const ccomplex_t IT_0149 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0150 = sin(beta);
    const ccomplex_t IT_0151 = e_em*m_tau*IT_0007*IT_0028*IT_0030*IT_0150;
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0149*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0003*IT_0027;
    const ccomplex_t IT_0156 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0157 = IT_0155*IT_0156;
    const ccomplex_t IT_0158 = IT_0003*IT_0150;
    const ccomplex_t IT_0159 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0160 = IT_0158*IT_0159;
    const ccomplex_t IT_0161 = IT_0007*IT_0027;
    const ccomplex_t IT_0162 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0163 = IT_0161*IT_0162;
    const ccomplex_t IT_0164 = IT_0007*IT_0150;
    const ccomplex_t IT_0165 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0166 = IT_0164*IT_0165;
    const ccomplex_t IT_0167 = IT_0157 + -IT_0160 + -IT_0163 + IT_0166;
    const ccomplex_t IT_0168 = IT_0154*IT_0167;
    const ccomplex_t IT_0169 = pow(m_Z, 2);
    const ccomplex_t IT_0170 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0169 + 
      -reg_prop, -1);
    const ccomplex_t IT_0171 = e_em*m_tau*IT_0007*IT_0030;
    const ccomplex_t IT_0172 = 0.5*IT_0171;
    const ccomplex_t IT_0173 = IT_0170*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = IT_0155*IT_0159;
    const ccomplex_t IT_0176 = IT_0156*IT_0158;
    const ccomplex_t IT_0177 = IT_0161*IT_0165;
    const ccomplex_t IT_0178 = IT_0162*IT_0164;
    const ccomplex_t IT_0179 = IT_0175 + IT_0176 + -IT_0177 + -IT_0178;
    const ccomplex_t IT_0180 = IT_0174*IT_0179;
    const ccomplex_t IT_0181 = 0.5*IT_0168 + 0.5*IT_0180;
    const ccomplex_t IT_0182 = m_N_1*m_tau;
    const ccomplex_t IT_0183 = s_13 + s_14;
    const ccomplex_t IT_0184 = IT_0182*IT_0183;
    const ccomplex_t IT_0185 = 4*IT_0184;
    const ccomplex_t IT_0186 = IT_0181*IT_0185;
    const ccomplex_t IT_0187 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0188 = IT_0155*IT_0187;
    const ccomplex_t IT_0189 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0190 = IT_0158*IT_0189;
    const ccomplex_t IT_0191 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0192 = IT_0161*IT_0191;
    const ccomplex_t IT_0193 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0194 = IT_0164*IT_0193;
    const ccomplex_t IT_0195 = -IT_0188 + IT_0190 + IT_0192 + -IT_0194;
    const ccomplex_t IT_0196 = IT_0154*IT_0195;
    const ccomplex_t IT_0197 = IT_0155*IT_0189;
    const ccomplex_t IT_0198 = IT_0158*IT_0187;
    const ccomplex_t IT_0199 = IT_0161*IT_0193;
    const ccomplex_t IT_0200 = IT_0164*IT_0191;
    const ccomplex_t IT_0201 = -IT_0197 + -IT_0198 + IT_0199 + IT_0200;
    const ccomplex_t IT_0202 = IT_0174*IT_0201;
    const ccomplex_t IT_0203 = 0.5*IT_0196 + 0.5*IT_0202;
    const ccomplex_t IT_0204 = s_23 + s_24;
    const ccomplex_t IT_0205 = IT_0182*IT_0204;
    const ccomplex_t IT_0206 = (-4)*IT_0205;
    const ccomplex_t IT_0207 = IT_0203*IT_0206;
    const ccomplex_t IT_0208 = (-0.5)*IT_0168 + (-0.5)*IT_0180;
    const ccomplex_t IT_0209 = (-4)*IT_0184;
    const ccomplex_t IT_0210 = IT_0208*IT_0209;
    const ccomplex_t IT_0211 = IT_0052 + IT_0074 + IT_0106 + IT_0138 + IT_0148
       + IT_0186 + IT_0207 + IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0213 = IT_0001*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = IT_0021*IT_0214;
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0045 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0218 = IT_0139*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = (-0.5)*IT_0219;
    const ccomplex_t IT_0221 = IT_0051*IT_0220;
    const ccomplex_t IT_0222 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0045 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0223 = IT_0142*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = IT_0051*IT_0225;
    const ccomplex_t IT_0227 = IT_0086*IT_0217;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = IT_0100*IT_0222;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = 0.5*IT_0228 + 0.5*IT_0230;
    const ccomplex_t IT_0232 = s_13*m_N_1*m_tau;
    const ccomplex_t IT_0233 = 4*IT_0232;
    const ccomplex_t IT_0234 = IT_0231*IT_0233;
    const ccomplex_t IT_0235 = IT_0044*IT_0217;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0147*IT_0237;
    const ccomplex_t IT_0239 = 4*IT_0205;
    const ccomplex_t IT_0240 = IT_0203*IT_0239;
    const ccomplex_t IT_0241 = IT_0185*IT_0208;
    const ccomplex_t IT_0242 = s_13*s_24;
    const ccomplex_t IT_0243 = 16*IT_0242;
    const ccomplex_t IT_0244 = IT_0216*IT_0243;
    const ccomplex_t IT_0245 = IT_0181*IT_0209;
    const ccomplex_t IT_0246 = IT_0221 + IT_0226 + IT_0234 + IT_0238 + IT_0240
       + IT_0241 + IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = IT_0069*IT_0222;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = conj(IT_0237) + conj(IT_0249);
    const ccomplex_t IT_0251 = 2*IT_0146;
    const ccomplex_t IT_0252 = IT_0073*IT_0251;
    const ccomplex_t IT_0253 = conj(IT_0049) + conj(IT_0073);
    const ccomplex_t IT_0254 = 4*IT_0242;
    const ccomplex_t IT_0255 = IT_0073*IT_0254;
    const ccomplex_t IT_0256 = 0.5*IT_0020;
    const ccomplex_t IT_0257 = IT_0214*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0013*IT_0256;
    const ccomplex_t IT_0260 = (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = (-0.5)*IT_0022;
    const ccomplex_t IT_0262 = 0.5*IT_0259;
    const ccomplex_t IT_0263 = conj(IT_0261) + conj(IT_0262);
    const ccomplex_t IT_0264 = conj(IT_0023) + conj(IT_0260);
    const ccomplex_t IT_0265 = conj(IT_0216) + conj(IT_0258);
    const ccomplex_t IT_0266 = (-0.5)*IT_0215;
    const ccomplex_t IT_0267 = 0.5*IT_0257;
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = conj(IT_0266) + conj(IT_0267);
    const ccomplex_t IT_0270 = conj(IT_0145)*(IT_0216 + IT_0258) + IT_0250*
      (IT_0023 + IT_0260) + IT_0253*(IT_0261 + IT_0262) + (IT_0049 + IT_0073)
      *IT_0263 + (IT_0237 + IT_0249)*IT_0264 + IT_0145*IT_0265 + (conj(IT_0220) 
      + conj(IT_0225))*IT_0268 + (IT_0220 + IT_0225)*IT_0269;
    const ccomplex_t IT_0271 = IT_0000*IT_0045;
    const ccomplex_t IT_0272 = 8*IT_0271;
    const ccomplex_t IT_0273 = 2*IT_0205;
    const ccomplex_t IT_0274 = 2*IT_0184;
    const ccomplex_t IT_0275 = s_34 + IT_0045;
    const ccomplex_t IT_0276 = IT_0000*IT_0275;
    const ccomplex_t IT_0277 = 2*IT_0276;
    const ccomplex_t IT_0278 = 8*IT_0276;
    const ccomplex_t IT_0279 = s_12*IT_0275;
    const ccomplex_t IT_0280 = 8*IT_0279;
    const ccomplex_t IT_0281 = (-0.5)*IT_0196 + (-0.5)*IT_0202;
    const ccomplex_t IT_0282 = (-8)*IT_0276;
    const ccomplex_t IT_0283 = (-8)*IT_0279;
    const ccomplex_t IT_0284 = s_12*s_34;
    const ccomplex_t IT_0285 = s_14*s_23;
    const ccomplex_t IT_0286 = -IT_0242;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = IT_0284 + IT_0287;
    const ccomplex_t IT_0289 = IT_0050 + IT_0288;
    const ccomplex_t IT_0290 = (-2)*IT_0289;
    const ccomplex_t IT_0291 = -IT_0284;
    const ccomplex_t IT_0292 = IT_0287 + IT_0291;
    const ccomplex_t IT_0293 = -IT_0050;
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = (-2)*IT_0294;
    const ccomplex_t IT_0296 = (-2)*IT_0205;
    const ccomplex_t IT_0297 = conj(IT_0220) + conj(IT_0225);
    const ccomplex_t IT_0298 = (-2)*IT_0184;
    const ccomplex_t IT_0299 = IT_0132*IT_0217;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0118*IT_0222;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = 0.5*IT_0300 + 0.5*IT_0302;
    const ccomplex_t IT_0304 = (-2)*IT_0276;
    const ccomplex_t IT_0305 = conj(IT_0145)*IT_0273 + IT_0253*IT_0274 + conj
      (IT_0135)*IT_0277 + conj(IT_0203)*IT_0278 + conj(IT_0208)*IT_0280 + conj
      (IT_0281)*IT_0282 + conj(IT_0181)*IT_0283 + conj(IT_0231)*IT_0290 + conj
      (IT_0103)*IT_0295 + IT_0250*IT_0296 + IT_0297*IT_0298 + conj(IT_0303)
      *IT_0304;
    const ccomplex_t IT_0306 = 2*IT_0289;
    const ccomplex_t IT_0307 = 2*IT_0294;
    const ccomplex_t IT_0308 = IT_0209*IT_0263 + IT_0239*IT_0264 + IT_0206
      *IT_0265 + IT_0185*IT_0269 + IT_0250*IT_0274 + conj(IT_0231)*IT_0277 +
       conj(IT_0181)*IT_0278 + IT_0280*conj(IT_0281) + conj(IT_0208)*IT_0282 +
       conj(IT_0203)*IT_0283 + IT_0253*IT_0296 + IT_0273*IT_0297 + conj(IT_0145)
      *IT_0298 + conj(IT_0103)*IT_0304 + conj(IT_0303)*IT_0306 + conj(IT_0135)
      *IT_0307;
    const ccomplex_t IT_0309 = IT_0185*IT_0263 + IT_0253*IT_0273 + conj
      (IT_0145)*IT_0274 + conj(IT_0103)*IT_0277 + conj(IT_0208)*IT_0278 + conj
      (IT_0203)*IT_0280 + conj(IT_0181)*IT_0282 + conj(IT_0281)*IT_0283 + conj
      (IT_0135)*IT_0295 + IT_0296*IT_0297 + IT_0250*IT_0298 + IT_0290*conj
      (IT_0303) + conj(IT_0231)*IT_0304;
    const ccomplex_t IT_0310 = IT_0242 + IT_0285 + IT_0291;
    const ccomplex_t IT_0311 = 2*IT_0310;
    const ccomplex_t IT_0312 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0313 = cos(alpha);
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*e_em*m_tau*IT_0007*IT_0028
      *IT_0030*IT_0313;
    const ccomplex_t IT_0315 = (-0.5)*IT_0314;
    const ccomplex_t IT_0316 = IT_0312*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = IT_0003*IT_0313;
    const ccomplex_t IT_0319 = IT_0159*IT_0318;
    const ccomplex_t IT_0320 = sin(alpha);
    const ccomplex_t IT_0321 = IT_0003*IT_0320;
    const ccomplex_t IT_0322 = IT_0156*IT_0321;
    const ccomplex_t IT_0323 = IT_0007*IT_0313;
    const ccomplex_t IT_0324 = IT_0165*IT_0323;
    const ccomplex_t IT_0325 = IT_0007*IT_0320;
    const ccomplex_t IT_0326 = IT_0162*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*(IT_0319 + -IT_0322 + 
      -IT_0324 + IT_0326);
    const ccomplex_t IT_0328 = IT_0317*IT_0327;
    const ccomplex_t IT_0329 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*e_em*m_tau*IT_0007*IT_0028
      *IT_0030*IT_0320;
    const ccomplex_t IT_0331 = 0.5*IT_0330;
    const ccomplex_t IT_0332 = IT_0329*IT_0331;
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*IT_0332;
    const ccomplex_t IT_0334 = IT_0156*IT_0318;
    const ccomplex_t IT_0335 = IT_0159*IT_0321;
    const ccomplex_t IT_0336 = IT_0162*IT_0323;
    const ccomplex_t IT_0337 = IT_0165*IT_0325;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*(IT_0334 + IT_0335 + 
      -IT_0336 + -IT_0337);
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = IT_0333*IT_0339;
    const ccomplex_t IT_0341 = 0.5*IT_0328 + 0.5*IT_0340;
    const ccomplex_t IT_0342 = IT_0050 + IT_0292;
    const ccomplex_t IT_0343 = 2*IT_0342;
    const ccomplex_t IT_0344 = 2*IT_0104;
    const ccomplex_t IT_0345 = 2*IT_0136;
    const ccomplex_t IT_0346 = 2*IT_0271;
    const ccomplex_t IT_0347 = IT_0187*IT_0318;
    const ccomplex_t IT_0348 = IT_0189*IT_0321;
    const ccomplex_t IT_0349 = IT_0191*IT_0323;
    const ccomplex_t IT_0350 = IT_0193*IT_0325;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*(IT_0347 + IT_0348 + 
      -IT_0349 + -IT_0350);
    const ccomplex_t IT_0352 = -IT_0351;
    const ccomplex_t IT_0353 = IT_0333*IT_0352;
    const ccomplex_t IT_0354 = IT_0189*IT_0318;
    const ccomplex_t IT_0355 = IT_0187*IT_0321;
    const ccomplex_t IT_0356 = IT_0193*IT_0323;
    const ccomplex_t IT_0357 = IT_0191*IT_0325;
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*(IT_0354 + -IT_0355 + 
      -IT_0356 + IT_0357);
    const ccomplex_t IT_0359 = IT_0317*IT_0358;
    const ccomplex_t IT_0360 = 0.5*IT_0353 + 0.5*IT_0359;
    const ccomplex_t IT_0361 = -IT_0045;
    const ccomplex_t IT_0362 = s_34 + IT_0361;
    const ccomplex_t IT_0363 = IT_0000*IT_0362;
    const ccomplex_t IT_0364 = 2*IT_0363;
    const ccomplex_t IT_0365 = s_24*m_N_1*m_tau;
    const ccomplex_t IT_0366 = 4*IT_0365;
    const ccomplex_t IT_0367 = 4*IT_0271;
    const ccomplex_t IT_0368 = 8*IT_0365;
    const ccomplex_t IT_0369 = 8*IT_0232;
    const ccomplex_t IT_0370 = (-0.5)*IT_0328 + (-0.5)*IT_0340;
    const ccomplex_t IT_0371 = (-2)*IT_0342;
    const ccomplex_t IT_0372 = (-0.5)*IT_0353 + (-0.5)*IT_0359;
    const ccomplex_t IT_0373 = (-2)*IT_0363;
    const ccomplex_t IT_0374 = IT_0233*IT_0253 + conj(IT_0103)*IT_0254 + conj
      (IT_0203)*IT_0277 + conj(IT_0208)*IT_0295 + conj(IT_0281)*IT_0304 + conj
      (IT_0181)*IT_0307 + conj(IT_0231)*IT_0311 + conj(IT_0341)*IT_0343 +
       IT_0297*IT_0344 + IT_0250*IT_0345 + conj(IT_0303)*IT_0346 + conj(IT_0360)
      *IT_0364 + conj(IT_0145)*IT_0366 + conj(IT_0135)*IT_0367 + IT_0263*IT_0368
       + IT_0265*IT_0369 + conj(IT_0370)*IT_0371 + conj(IT_0372)*IT_0373;
    const ccomplex_t IT_0375 = conj(IT_0145)*IT_0233 + conj(IT_0135)*IT_0254 +
       conj(IT_0208)*IT_0277 + conj(IT_0203)*IT_0295 + conj(IT_0181)*IT_0304 +
       conj(IT_0281)*IT_0307 + conj(IT_0303)*IT_0311 + IT_0250*IT_0344 + IT_0297
      *IT_0345 + conj(IT_0231)*IT_0346 + IT_0343*conj(IT_0360) + conj(IT_0341)
      *IT_0364 + IT_0253*IT_0366 + conj(IT_0103)*IT_0367 + IT_0265*IT_0368 +
       IT_0263*IT_0369 + IT_0371*conj(IT_0372) + conj(IT_0370)*IT_0373;
    const ccomplex_t IT_0376 = IT_0250*IT_0273 + conj(IT_0181)*IT_0280 +
       IT_0278*conj(IT_0281) + conj(IT_0203)*IT_0282 + conj(IT_0208)*IT_0283 +
       conj(IT_0145)*IT_0296 + IT_0274*IT_0297 + IT_0253*IT_0298 + IT_0277*conj
      (IT_0303) + conj(IT_0135)*IT_0304 + conj(IT_0231)*IT_0306 + conj(IT_0103)
      *IT_0307;
    const ccomplex_t IT_0377 = -s_24;
    const ccomplex_t IT_0378 = s_23 + IT_0377;
    const ccomplex_t IT_0379 = IT_0182*IT_0378;
    const ccomplex_t IT_0380 = (-2)*IT_0379;
    const ccomplex_t IT_0381 = -s_14;
    const ccomplex_t IT_0382 = s_13 + IT_0381;
    const ccomplex_t IT_0383 = IT_0182*IT_0382;
    const ccomplex_t IT_0384 = 2*IT_0383;
    const ccomplex_t IT_0385 = 2*IT_0050;
    const ccomplex_t IT_0386 = 2*IT_0365;
    const ccomplex_t IT_0387 = 2*IT_0232;
    const ccomplex_t IT_0388 = -IT_0380;
    const ccomplex_t IT_0389 = -IT_0384;
    const ccomplex_t IT_0390 = conj(IT_0135)*IT_0233 + conj(IT_0145)*IT_0254 +
       8*IT_0242*IT_0263 + conj(IT_0208)*IT_0273 + conj(IT_0203)*IT_0274 + conj
      (IT_0181)*IT_0296 + IT_0251*IT_0297 + conj(IT_0281)*IT_0298 + conj(IT_0103
      )*IT_0366 + IT_0253*IT_0367 + conj(IT_0370)*IT_0380 + conj(IT_0372)
      *IT_0384 + IT_0250*IT_0385 + conj(IT_0231)*IT_0386 + conj(IT_0303)*IT_0387
       + conj(IT_0341)*IT_0388 + conj(IT_0360)*IT_0389;
    const ccomplex_t IT_0391 = IT_0049*IT_0385;
    const ccomplex_t IT_0392 = IT_0073*IT_0385;
    const ccomplex_t IT_0393 = 4*IT_0285;
    const ccomplex_t IT_0394 = IT_0225*IT_0393;
    const ccomplex_t IT_0395 = IT_0237*IT_0367;
    const ccomplex_t IT_0396 = 2*IT_0379;
    const ccomplex_t IT_0397 = -IT_0396;
    const ccomplex_t IT_0398 = IT_0372*IT_0397;
    const ccomplex_t IT_0399 = (-2)*IT_0383;
    const ccomplex_t IT_0400 = IT_0341*IT_0399;
    const ccomplex_t IT_0401 = 4*IT_0379;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = IT_0372*IT_0402;
    const ccomplex_t IT_0404 = (-4)*IT_0383;
    const ccomplex_t IT_0405 = IT_0341*IT_0404;
    const ccomplex_t IT_0406 = IT_0051*IT_0237;
    const ccomplex_t IT_0407 = IT_0231*IT_0366;
    const ccomplex_t IT_0408 = IT_0147*IT_0220;
    const ccomplex_t IT_0409 = IT_0147*IT_0225;
    const ccomplex_t IT_0410 = IT_0208*IT_0239;
    const ccomplex_t IT_0411 = IT_0181*IT_0206;
    const ccomplex_t IT_0412 = IT_0406 + IT_0407 + IT_0408 + IT_0409 + IT_0410
       + IT_0411;
    const ccomplex_t IT_0413 = -IT_0404;
    const ccomplex_t IT_0414 = IT_0372*IT_0413;
    const ccomplex_t IT_0415 = IT_0341*IT_0401;
    const ccomplex_t IT_0416 = IT_0051*IT_0145;
    const ccomplex_t IT_0417 = IT_0105*IT_0135;
    const ccomplex_t IT_0418 = IT_0103*IT_0137;
    const ccomplex_t IT_0419 = IT_0049*IT_0147;
    const ccomplex_t IT_0420 = IT_0073*IT_0147;
    const ccomplex_t IT_0421 = IT_0181*IT_0239;
    const ccomplex_t IT_0422 = 16*IT_0285;
    const ccomplex_t IT_0423 = IT_0266*IT_0422;
    const ccomplex_t IT_0424 = 16*IT_0271;
    const ccomplex_t IT_0425 = IT_0023*IT_0424;
    const ccomplex_t IT_0426 = IT_0260*IT_0424;
    const ccomplex_t IT_0427 = IT_0206*IT_0208;
    const ccomplex_t IT_0428 = IT_0203*IT_0209;
    const ccomplex_t IT_0429 = (-4)*IT_0379;
    const ccomplex_t IT_0430 = -IT_0429;
    const ccomplex_t IT_0431 = IT_0341*IT_0430;
    const ccomplex_t IT_0432 = 4*IT_0383;
    const ccomplex_t IT_0433 = IT_0372*IT_0432;
    const ccomplex_t IT_0434 = s_12*IT_0362;
    const ccomplex_t IT_0435 = 0.125*IT_0384;
    const ccomplex_t IT_0436 = conj(IT_0023) + conj(IT_0216);
    const ccomplex_t IT_0437 = 0.125*IT_0429;
    const ccomplex_t IT_0438 = (-8)*IT_0434;
    const ccomplex_t IT_0439 = 0.125*conj(IT_0360);
    const ccomplex_t IT_0440 = (-8)*IT_0363;
    const ccomplex_t IT_0441 = 0.125*conj(IT_0370);
    const ccomplex_t IT_0442 = IT_0288 + IT_0293;
    const ccomplex_t IT_0443 = (-2)*IT_0442;
    const ccomplex_t IT_0444 = 0.125*conj(IT_0303);
    const ccomplex_t IT_0445 = -IT_0399;
    const ccomplex_t IT_0446 = 0.125*IT_0250;
    const ccomplex_t IT_0447 = 0.125*conj(IT_0103);
    const ccomplex_t IT_0448 = 0.125*conj(IT_0231);
    const ccomplex_t IT_0449 = 0.125*conj(IT_0135);
    const ccomplex_t IT_0450 = conj(IT_0341)*IT_0363 + 0.125*IT_0253*IT_0380 +
       conj(IT_0372)*IT_0434 + conj(IT_0145)*IT_0435 + IT_0436*IT_0437 + IT_0438
      *IT_0439 + IT_0440*IT_0441 + IT_0443*IT_0444 + IT_0445*IT_0446 + IT_0373
      *IT_0447 + IT_0364*IT_0448 + IT_0371*IT_0449;
    const ccomplex_t IT_0451 = conj(IT_0261) + conj(IT_0266);
    const ccomplex_t IT_0452 = conj(IT_0258) + conj(IT_0260);
    const ccomplex_t IT_0453 = conj(IT_0262) + conj(IT_0267);
    const ccomplex_t IT_0454 = -IT_0432;
    const ccomplex_t IT_0455 = 0.125*IT_0396;
    const ccomplex_t IT_0456 = 0.125*conj(IT_0341);
    const ccomplex_t IT_0457 = 2*IT_0442;
    const ccomplex_t IT_0458 = 0.125*conj(IT_0145);
    const ccomplex_t IT_0459 = 0.125*IT_0436;
    const ccomplex_t IT_0460 = 0.125*conj(IT_0372);
    const ccomplex_t IT_0461 = 0.125*IT_0253;
    const ccomplex_t IT_0462 = IT_0363*conj(IT_0370) + 0.125*IT_0250*IT_0399 +
       conj(IT_0360)*IT_0434 + IT_0364*IT_0447 + IT_0373*IT_0448 + IT_0343
      *IT_0449 + 0.125*IT_0404*IT_0451 + 0.125*IT_0401*IT_0452 + 0.125*IT_0453
      *IT_0454 + IT_0297*IT_0455 + IT_0440*IT_0456 + IT_0444*IT_0457 + IT_0389
      *IT_0458 + IT_0430*IT_0459 + IT_0438*IT_0460 + IT_0388*IT_0461;
    const ccomplex_t IT_0463 = IT_0363*conj(IT_0372) + conj(IT_0341)*IT_0434 +
       IT_0439*IT_0440 + IT_0438*IT_0441 + IT_0373*IT_0444 + IT_0343*IT_0447 +
       IT_0364*IT_0449 + IT_0250*IT_0455 + IT_0448*IT_0457 + IT_0388*IT_0458 +
       IT_0454*IT_0459 + IT_0389*IT_0461;
    const ccomplex_t IT_0464 = 0.125*conj(IT_0281);
    const ccomplex_t IT_0465 = 0.125*conj(IT_0203);
    const ccomplex_t IT_0466 = 0.125*conj(IT_0208);
    const ccomplex_t IT_0467 = 0.125*conj(IT_0181);
    const ccomplex_t IT_0468 = 0.125*IT_0297;
    const ccomplex_t IT_0469 = IT_0104*IT_0264 + IT_0136*IT_0269 + IT_0373
      *IT_0439 + IT_0441*IT_0443 + IT_0367*IT_0444 + IT_0137*IT_0446 + IT_0311
      *IT_0447 + IT_0393*IT_0448 + IT_0346*IT_0449 + IT_0456*IT_0457 + IT_0386
      *IT_0458 + IT_0364*IT_0460 + IT_0387*IT_0461 + IT_0277*IT_0464 + IT_0304
      *IT_0465 + IT_0290*IT_0466 + IT_0306*IT_0467 + IT_0105*IT_0468;
    const ccomplex_t IT_0470 = 0.125*IT_0233*IT_0263 + IT_0136*IT_0264 +
       IT_0104*IT_0269 + 0.125*IT_0265*IT_0366 + IT_0364*IT_0441 + IT_0393
      *IT_0444 + IT_0105*IT_0446 + IT_0346*IT_0447 + IT_0367*IT_0448 + IT_0311
      *IT_0449 + IT_0373*IT_0456 + IT_0439*IT_0457 + IT_0387*IT_0458 + IT_0443
      *IT_0460 + IT_0386*IT_0461 + IT_0306*IT_0464 + IT_0290*IT_0465 + IT_0304
      *IT_0466 + IT_0277*IT_0467 + IT_0137*IT_0468;
    const ccomplex_t IT_0471 = conj(IT_0360)*IT_0363 + 0.125*conj(IT_0145)
      *IT_0380 + conj(IT_0370)*IT_0434 + IT_0253*IT_0435 + 0.125*IT_0432*IT_0436
       + IT_0364*IT_0444 + IT_0397*IT_0446 + IT_0371*IT_0447 + IT_0443*IT_0448 +
       IT_0373*IT_0449 + 0.125*IT_0402*IT_0451 + 0.125*IT_0413*IT_0452 + IT_0437
      *IT_0453 + IT_0438*IT_0456 + IT_0440*IT_0460 + IT_0445*IT_0468;
    const ccomplex_t IT_0472 = IT_0344*IT_0447;
    const ccomplex_t IT_0473 = IT_0105*IT_0448;
    const ccomplex_t IT_0474 = IT_0296*IT_0465;
    const ccomplex_t IT_0475 = IT_0298*IT_0466;
    const ccomplex_t IT_0476 = conj(IT_0049)*IT_0385;
    const ccomplex_t IT_0477 = 0.125*IT_0476;
    const ccomplex_t IT_0478 = conj(IT_0073)*IT_0385;
    const ccomplex_t IT_0479 = 0.125*IT_0478;
    const ccomplex_t IT_0480 = conj(IT_0135)*IT_0345;
    const ccomplex_t IT_0481 = 0.125*IT_0480;
    const ccomplex_t IT_0482 = IT_0137*conj(IT_0303);
    const ccomplex_t IT_0483 = 0.125*IT_0482;
    const ccomplex_t IT_0484 = conj(IT_0249)*IT_0367;
    const ccomplex_t IT_0485 = 0.125*IT_0484;
    const ccomplex_t IT_0486 = conj(IT_0237)*IT_0367;
    const ccomplex_t IT_0487 = 0.125*IT_0486;
    const ccomplex_t IT_0488 = conj(IT_0370)*IT_0445;
    const ccomplex_t IT_0489 = 0.125*IT_0488;
    const ccomplex_t IT_0490 = conj(IT_0372)*IT_0397;
    const ccomplex_t IT_0491 = 0.125*IT_0490;
    const ccomplex_t IT_0492 = conj(IT_0341)*IT_0399;
    const ccomplex_t IT_0493 = 0.125*IT_0492;
    const ccomplex_t IT_0494 = conj(IT_0360)*IT_0396;
    const ccomplex_t IT_0495 = 0.125*IT_0494;
    const ccomplex_t IT_0496 = IT_0264*IT_0285;
    const ccomplex_t IT_0497 = IT_0251*IT_0458 + IT_0273*IT_0464 + IT_0274
      *IT_0467 + IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0477 + IT_0479 +
       IT_0481 + IT_0483 + IT_0485 + IT_0487 + IT_0489 + IT_0491 + IT_0493 +
       IT_0495 + IT_0496;
    const ccomplex_t IT_0498 = IT_0251*IT_0458 + IT_0273*IT_0464 + IT_0274
      *IT_0467 + IT_0393*IT_0468 + IT_0472 + IT_0473 + IT_0474 + IT_0475 +
       IT_0477 + IT_0479 + IT_0481 + IT_0483 + IT_0485 + IT_0487 + IT_0489 +
       IT_0491 + IT_0493 + IT_0495 + IT_0496;
    const ccomplex_t IT_0499 = conj(IT_0049)*IT_0251;
    const ccomplex_t IT_0500 = 0.125*IT_0499;
    const ccomplex_t IT_0501 = IT_0393*IT_0446;
    const ccomplex_t IT_0502 = IT_0345*IT_0447;
    const ccomplex_t IT_0503 = IT_0137*IT_0448;
    const ccomplex_t IT_0504 = IT_0274*IT_0464;
    const ccomplex_t IT_0505 = IT_0298*IT_0465;
    const ccomplex_t IT_0506 = IT_0296*IT_0466;
    const ccomplex_t IT_0507 = conj(IT_0135)*IT_0344;
    const ccomplex_t IT_0508 = 0.125*IT_0507;
    const ccomplex_t IT_0509 = conj(IT_0073)*IT_0251;
    const ccomplex_t IT_0510 = 0.125*IT_0509;
    const ccomplex_t IT_0511 = IT_0105*conj(IT_0303);
    const ccomplex_t IT_0512 = 0.125*IT_0511;
    const ccomplex_t IT_0513 = conj(IT_0372)*IT_0445;
    const ccomplex_t IT_0514 = 0.125*IT_0513;
    const ccomplex_t IT_0515 = conj(IT_0370)*IT_0397;
    const ccomplex_t IT_0516 = 0.125*IT_0515;
    const ccomplex_t IT_0517 = conj(IT_0360)*IT_0399;
    const ccomplex_t IT_0518 = 0.125*IT_0517;
    const ccomplex_t IT_0519 = conj(IT_0341)*IT_0396;
    const ccomplex_t IT_0520 = 0.125*IT_0519;
    const ccomplex_t IT_0521 = IT_0269*IT_0285;
    const ccomplex_t IT_0522 = 8*IT_0237*(IT_0385*IT_0458 + IT_0273*IT_0467 +
       IT_0500 + IT_0501 + IT_0502 + IT_0503 + IT_0504 + IT_0505 + IT_0506 +
       IT_0508 + IT_0510 + IT_0512 + IT_0514 + IT_0516 + IT_0518 + IT_0520 +
       IT_0521);
    const ccomplex_t IT_0523 = IT_0051*IT_0263 + IT_0147*IT_0265 + 8*IT_0385
      *IT_0458 + 8*IT_0273*IT_0467 + 8*IT_0367*IT_0468 + 8*IT_0500 + 8*IT_0501 +
       8*IT_0502 + 8*IT_0503 + 8*IT_0504 + 8*IT_0505 + 8*IT_0506 + 8*IT_0508 + 8
      *IT_0510 + 8*IT_0512 + 8*IT_0514 + 8*IT_0516 + 8*IT_0518 + 8*IT_0520 + 8
      *IT_0521;
    const ccomplex_t IT_0524 = IT_0233*IT_0447;
    const ccomplex_t IT_0525 = IT_0387*IT_0448;
    const ccomplex_t IT_0526 = IT_0296*IT_0464;
    const ccomplex_t IT_0527 = IT_0273*IT_0465;
    const ccomplex_t IT_0528 = IT_0274*IT_0466;
    const ccomplex_t IT_0529 = IT_0298*IT_0467;
    const ccomplex_t IT_0530 = conj(IT_0303)*IT_0386;
    const ccomplex_t IT_0531 = 0.125*IT_0530;
    const ccomplex_t IT_0532 = conj(IT_0135)*IT_0366;
    const ccomplex_t IT_0533 = 0.125*IT_0532;
    const ccomplex_t IT_0534 = conj(IT_0360)*IT_0388;
    const ccomplex_t IT_0535 = 0.125*IT_0534;
    const ccomplex_t IT_0536 = conj(IT_0341)*IT_0389;
    const ccomplex_t IT_0537 = 0.125*IT_0536;
    const ccomplex_t IT_0538 = conj(IT_0372)*IT_0380;
    const ccomplex_t IT_0539 = 0.125*IT_0538;
    const ccomplex_t IT_0540 = conj(IT_0370)*IT_0384;
    const ccomplex_t IT_0541 = 0.125*IT_0540;
    const ccomplex_t IT_0542 = IT_0242*IT_0265;
    const ccomplex_t IT_0543 = 8*IT_0049*(IT_0251*IT_0446 + IT_0367*IT_0458 +
       IT_0254*IT_0461 + IT_0524 + IT_0525 + IT_0526 + IT_0527 + IT_0528 +
       IT_0529 + IT_0531 + IT_0533 + IT_0535 + IT_0537 + IT_0539 + IT_0541 +
       IT_0542);
    const ccomplex_t IT_0544 = IT_0051*conj(IT_0249);
    const ccomplex_t IT_0545 = IT_0051*conj(IT_0237);
    const ccomplex_t IT_0546 = conj(IT_0231)*IT_0366;
    const ccomplex_t IT_0547 = IT_0233*conj(IT_0303);
    const ccomplex_t IT_0548 = IT_0147*conj(IT_0220);
    const ccomplex_t IT_0549 = IT_0147*conj(IT_0225);
    const ccomplex_t IT_0550 = conj(IT_0208)*IT_0239;
    const ccomplex_t IT_0551 = IT_0185*conj(IT_0203);
    const ccomplex_t IT_0552 = IT_0209*conj(IT_0281);
    const ccomplex_t IT_0553 = conj(IT_0372)*IT_0413;
    const ccomplex_t IT_0554 = conj(IT_0370)*IT_0402;
    const ccomplex_t IT_0555 = conj(IT_0360)*IT_0404;
    const ccomplex_t IT_0556 = conj(IT_0341)*IT_0401;
    const ccomplex_t IT_0557 = IT_0050*IT_0269;
    const ccomplex_t IT_0558 = conj(IT_0145)*IT_0242 + 0.125*IT_0243*IT_0263 +
       IT_0146*IT_0264 + 0.125*IT_0265*IT_0424 + IT_0368*IT_0447 + IT_0369
      *IT_0449 + IT_0206*IT_0467 + 0.125*IT_0544 + 0.125*IT_0545 + 0.125*IT_0546
       + 0.125*IT_0547 + 0.125*IT_0548 + 0.125*IT_0549 + 0.125*IT_0550 + 0.125
      *IT_0551 + 0.125*IT_0552 + 0.125*IT_0553 + 0.125*IT_0554 + 0.125*IT_0555 +
       0.125*IT_0556 + IT_0557;
    const ccomplex_t IT_0559 = conj(IT_0341)*IT_0430;
    const ccomplex_t IT_0560 = conj(IT_0360)*IT_0454;
    const ccomplex_t IT_0561 = conj(IT_0370)*IT_0429;
    const ccomplex_t IT_0562 = conj(IT_0372)*IT_0432;
    const ccomplex_t IT_0563 = conj(IT_0145)*IT_0242 + 0.125*IT_0243*IT_0263 +
       IT_0146*IT_0264 + 0.125*IT_0265*IT_0424 + IT_0368*IT_0447 + IT_0369
      *IT_0449 + IT_0206*IT_0467 + 0.125*IT_0544 + 0.125*IT_0545 + 0.125*IT_0546
       + 0.125*IT_0547 + 0.125*IT_0548 + 0.125*IT_0549 + 0.125*IT_0550 + 0.125
      *IT_0551 + 0.125*IT_0552 + IT_0557 + 0.125*IT_0559 + 0.125*IT_0560 + 0.125
      *IT_0561 + 0.125*IT_0562;
    const ccomplex_t IT_0564 = conj(IT_0049)*IT_0051;
    const ccomplex_t IT_0565 = IT_0051*conj(IT_0073);
    const ccomplex_t IT_0566 = conj(IT_0103)*IT_0105;
    const ccomplex_t IT_0567 = conj(IT_0135)*IT_0137;
    const ccomplex_t IT_0568 = conj(IT_0145)*IT_0147;
    const ccomplex_t IT_0569 = IT_0239*conj(IT_0281);
    const ccomplex_t IT_0570 = conj(IT_0181)*IT_0185;
    const ccomplex_t IT_0571 = conj(IT_0203)*IT_0206;
    const ccomplex_t IT_0572 = conj(IT_0208)*IT_0209;
    const ccomplex_t IT_0573 = conj(IT_0370)*IT_0413;
    const ccomplex_t IT_0574 = conj(IT_0372)*IT_0402;
    const ccomplex_t IT_0575 = conj(IT_0341)*IT_0404;
    const ccomplex_t IT_0576 = conj(IT_0360)*IT_0401;
    const ccomplex_t IT_0577 = IT_0285*IT_0297;
    const ccomplex_t IT_0578 = IT_0050*IT_0265;
    const ccomplex_t IT_0579 = IT_0104*conj(IT_0231);
    const ccomplex_t IT_0580 = IT_0146*IT_0263;
    const ccomplex_t IT_0581 = IT_0136*conj(IT_0303) + 0.125*IT_0264*IT_0422 +
       0.125*IT_0564 + 0.125*IT_0565 + 0.125*IT_0566 + 0.125*IT_0567 + 0.125
      *IT_0568 + 0.125*IT_0569 + 0.125*IT_0570 + 0.125*IT_0571 + 0.125*IT_0572 +
       0.125*IT_0573 + 0.125*IT_0574 + 0.125*IT_0575 + 0.125*IT_0576 + IT_0577 +
       IT_0578 + IT_0579 + IT_0580;
    const ccomplex_t IT_0582 = conj(IT_0360)*IT_0430;
    const ccomplex_t IT_0583 = conj(IT_0341)*IT_0454;
    const ccomplex_t IT_0584 = conj(IT_0372)*IT_0429;
    const ccomplex_t IT_0585 = conj(IT_0370)*IT_0432;
    const ccomplex_t IT_0586 = IT_0136*conj(IT_0303) + 0.125*IT_0264*IT_0422 +
       0.125*IT_0564 + 0.125*IT_0565 + 0.125*IT_0566 + 0.125*IT_0567 + 0.125
      *IT_0568 + 0.125*IT_0569 + 0.125*IT_0570 + 0.125*IT_0571 + 0.125*IT_0572 +
       IT_0577 + IT_0578 + IT_0579 + IT_0580 + 0.125*IT_0582 + 0.125*IT_0583 +
       0.125*IT_0584 + 0.125*IT_0585;
    const ccomplex_t IT_0587 = IT_0051*conj(IT_0145);
    const ccomplex_t IT_0588 = IT_0105*conj(IT_0135);
    const ccomplex_t IT_0589 = conj(IT_0103)*IT_0137;
    const ccomplex_t IT_0590 = conj(IT_0049)*IT_0147;
    const ccomplex_t IT_0591 = conj(IT_0073)*IT_0147;
    const ccomplex_t IT_0592 = conj(IT_0181)*IT_0239;
    const ccomplex_t IT_0593 = IT_0185*conj(IT_0281);
    const ccomplex_t IT_0594 = conj(IT_0023)*IT_0424;
    const ccomplex_t IT_0595 = conj(IT_0260)*IT_0424;
    const ccomplex_t IT_0596 = IT_0206*conj(IT_0208);
    const ccomplex_t IT_0597 = conj(IT_0203)*IT_0209;
    const ccomplex_t IT_0598 = IT_0050*IT_0263;
    const ccomplex_t IT_0599 = IT_0136*conj(IT_0231);
    const ccomplex_t IT_0600 = IT_0146*IT_0265;
    const ccomplex_t IT_0601 = IT_0250*IT_0285 + IT_0104*conj(IT_0303) + 0.125
      *IT_0553 + 0.125*IT_0554 + 0.125*IT_0555 + 0.125*IT_0556 + 0.125*IT_0587 +
       0.125*IT_0588 + 0.125*IT_0589 + 0.125*IT_0590 + 0.125*IT_0591 + 0.125
      *IT_0592 + 0.125*IT_0593 + 0.125*IT_0594 + 0.125*IT_0595 + 0.125*IT_0596 +
       0.125*IT_0597 + IT_0598 + IT_0599 + IT_0600;
    const ccomplex_t IT_0602 = IT_0250*IT_0285 + IT_0104*conj(IT_0303) + 0.125
      *IT_0269*IT_0422 + 0.125*IT_0559 + 0.125*IT_0560 + 0.125*IT_0561 + 0.125
      *IT_0562 + 0.125*IT_0587 + 0.125*IT_0588 + 0.125*IT_0589 + 0.125*IT_0590 +
       0.125*IT_0591 + 0.125*IT_0592 + 0.125*IT_0593 + 0.125*IT_0594 + 0.125
      *IT_0595 + 0.125*IT_0596 + 0.125*IT_0597 + IT_0598 + IT_0599 + IT_0600;
    const ccomplex_t IT_0603 = IT_0051*conj(IT_0220);
    const ccomplex_t IT_0604 = IT_0051*conj(IT_0225);
    const ccomplex_t IT_0605 = conj(IT_0303)*IT_0366;
    const ccomplex_t IT_0606 = conj(IT_0231)*IT_0233;
    const ccomplex_t IT_0607 = IT_0147*conj(IT_0249);
    const ccomplex_t IT_0608 = IT_0147*conj(IT_0237);
    const ccomplex_t IT_0609 = conj(IT_0203)*IT_0239;
    const ccomplex_t IT_0610 = IT_0185*conj(IT_0208);
    const ccomplex_t IT_0611 = conj(IT_0261)*IT_0424;
    const ccomplex_t IT_0612 = IT_0206*conj(IT_0281);
    const ccomplex_t IT_0613 = conj(IT_0181)*IT_0209;
    const ccomplex_t IT_0614 = IT_0242*IT_0253;
    const ccomplex_t IT_0615 = IT_0146*IT_0269;
    const ccomplex_t IT_0616 = IT_0050*IT_0264 + 0.125*IT_0243*IT_0265 + 0.125
      *conj(IT_0262)*IT_0424 + IT_0369*IT_0447 + IT_0368*IT_0449 + 0.125*IT_0573
       + 0.125*IT_0574 + 0.125*IT_0575 + 0.125*IT_0576 + 0.125*IT_0603 + 0.125
      *IT_0604 + 0.125*IT_0605 + 0.125*IT_0606 + 0.125*IT_0607 + 0.125*IT_0608 +
       0.125*IT_0609 + 0.125*IT_0610 + 0.125*IT_0611 + 0.125*IT_0612 + 0.125
      *IT_0613 + IT_0614 + IT_0615;
    const ccomplex_t IT_0617 = IT_0050*IT_0264 + 0.125*conj(IT_0262)*IT_0424 +
       IT_0369*IT_0447 + IT_0368*IT_0449 + 0.125*IT_0582 + 0.125*IT_0583 + 0.125
      *IT_0584 + 0.125*IT_0585 + 0.125*IT_0603 + 0.125*IT_0604 + 0.125*IT_0605 +
       0.125*IT_0606 + 0.125*IT_0607 + 0.125*IT_0608 + 0.125*IT_0609 + 0.125
      *IT_0610 + 0.125*IT_0611 + 0.125*IT_0612 + 0.125*IT_0613 + IT_0614 +
       IT_0615;
    const ccomplex_t IT_0618 = conj(IT_0023)*IT_0211 + conj(IT_0216)*IT_0246 +
       IT_0250*IT_0252 + IT_0253*IT_0255 + IT_0270*IT_0272 + IT_0208*IT_0305 +
       IT_0281*IT_0308 + IT_0203*IT_0309 + IT_0103*IT_0374 + IT_0135*IT_0375 +
       IT_0181*IT_0376 + IT_0145*IT_0390 + IT_0297*(IT_0391 + IT_0392 + IT_0394 
      + IT_0395 + IT_0398 + IT_0400) + conj(IT_0260)*(IT_0211 + IT_0403 +
       IT_0405) + conj(IT_0258)*(IT_0246 + IT_0403 + IT_0405) + conj(IT_0261)*
      (IT_0412 + IT_0414 + IT_0415) + conj(IT_0266)*(IT_0414 + IT_0415 + IT_0416
       + IT_0417 + IT_0418 + IT_0419 + IT_0420 + IT_0421 + IT_0423 + IT_0425 +
       IT_0426 + IT_0427 + IT_0428) + conj(IT_0262)*(IT_0412 + IT_0431 + IT_0433
      ) + conj(IT_0267)*(IT_0416 + IT_0417 + IT_0418 + IT_0419 + IT_0420 +
       IT_0421 + IT_0423 + IT_0425 + IT_0426 + IT_0427 + IT_0428 + IT_0431 +
       IT_0433) + 8*IT_0372*IT_0450 + 8*IT_0360*IT_0462 + 8*IT_0341*IT_0463 + 8
      *IT_0231*IT_0469 + 8*IT_0303*IT_0470 + 8*IT_0370*IT_0471 + 8*IT_0225
      *IT_0497 + 8*IT_0220*IT_0498 + IT_0522 + IT_0249*IT_0523 + 8*IT_0073*
      (IT_0367*IT_0458 + IT_0524 + IT_0525 + IT_0526 + IT_0527 + IT_0528 +
       IT_0529 + IT_0531 + IT_0533 + IT_0535 + IT_0537 + IT_0539 + IT_0541 +
       IT_0542) + IT_0543 + 8*IT_0261*IT_0558 + 8*IT_0262*IT_0563 + 8*IT_0260
      *IT_0581 + 8*IT_0023*IT_0586 + 8*IT_0266*IT_0601 + 8*IT_0267*IT_0602 + 8
      *IT_0258*IT_0616 + 8*IT_0216*IT_0617;
    return create_ccomplex_return(IT_0618);
}

