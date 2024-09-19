#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_tau_tau(
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
    const creal_t m_N_4 = param->m_N_4;
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
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = m_N_4*m_tau;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 4*IT_0003;
    const ccomplex_t IT_0005 = cos(alpha);
    const ccomplex_t IT_0006 = cos(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = IT_0005*IT_0007;
    const ccomplex_t IT_0009 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = sin(alpha);
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0005*IT_0016;
    const ccomplex_t IT_0018 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0016;
    const ccomplex_t IT_0021 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*(IT_0010 + IT_0014 + 
      -IT_0019 + -IT_0022);
    const ccomplex_t IT_0024 = -IT_0023;
    const ccomplex_t IT_0025 = pow(m_W, -1);
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = cpow(IT_0026, -1);
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*e_em*m_tau*IT_0011*IT_0016
      *IT_0025*IT_0027;
    const ccomplex_t IT_0029 = 0.5*IT_0028;
    const ccomplex_t IT_0030 = pow(m_N_4, 2);
    const ccomplex_t IT_0031 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0030 + -reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = IT_0024*IT_0033;
    const ccomplex_t IT_0035 = IT_0008*IT_0013;
    const ccomplex_t IT_0036 = IT_0009*IT_0012;
    const ccomplex_t IT_0037 = IT_0017*IT_0021;
    const ccomplex_t IT_0038 = IT_0018*IT_0020;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + -IT_0036 + 
      -IT_0037 + IT_0038);
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*e_em*m_tau*IT_0005*IT_0016
      *IT_0025*IT_0027;
    const ccomplex_t IT_0041 = (-0.5)*IT_0040;
    const ccomplex_t IT_0042 = cpow((-2)*s_12 + (-2)*IT_0030 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0041*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0039*IT_0044;
    const ccomplex_t IT_0046 = (-0.5)*IT_0034 + (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = 0.5*IT_0034 + 0.5*IT_0045;
    const ccomplex_t IT_0048 = -conj(IT_0047);
    const ccomplex_t IT_0049 = conj(IT_0046) + IT_0048;
    const ccomplex_t IT_0050 = IT_0007*IT_0015;
    const ccomplex_t IT_0051 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = IT_0006*IT_0016;
    const ccomplex_t IT_0054 = IT_0051*IT_0053;
    const ccomplex_t IT_0055 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0056 = IT_0053*IT_0055;
    const ccomplex_t IT_0057 = IT_0050*IT_0055;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0052 + IT_0054 + 
      -IT_0056 + -IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = -IT_0059;
    const ccomplex_t IT_0061 = -IT_0060;
    const ccomplex_t IT_0062 = 0.5*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*e_em*IT_0007*IT_0015;
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0030 + -reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = IT_0062*IT_0066;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = e_em*IT_0050;
    const ccomplex_t IT_0070 = e_em*IT_0053;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = IT_0064*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0062*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0068 + IT_0076;
    const ccomplex_t IT_0078 = -IT_0047;
    const ccomplex_t IT_0079 = IT_0046 + IT_0078;
    const ccomplex_t IT_0080 = conj(IT_0068) + conj(IT_0076);
    const ccomplex_t IT_0081 = (-0.5)*IT_0061;
    const ccomplex_t IT_0082 = IT_0074*IT_0081;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0066*IT_0081;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = conj(IT_0083) + conj(IT_0085);
    const ccomplex_t IT_0087 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0088 = IT_0008*IT_0087;
    const ccomplex_t IT_0089 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0090 = IT_0012*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0092 = IT_0017*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0094 = IT_0020*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0088 + IT_0090 + 
      -IT_0092 + -IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0033*IT_0096;
    const ccomplex_t IT_0098 = IT_0008*IT_0089;
    const ccomplex_t IT_0099 = IT_0012*IT_0087;
    const ccomplex_t IT_0100 = IT_0017*IT_0093;
    const ccomplex_t IT_0101 = IT_0020*IT_0091;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + -IT_0099 + 
      -IT_0100 + IT_0101);
    const ccomplex_t IT_0103 = IT_0044*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0097 + (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = 0.5*IT_0097 + 0.5*IT_0103;
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0104 + IT_0106;
    const ccomplex_t IT_0108 = IT_0083 + IT_0085;
    const ccomplex_t IT_0109 = -conj(IT_0105);
    const ccomplex_t IT_0110 = conj(IT_0104) + IT_0109;
    const ccomplex_t IT_0111 = -s_24;
    const ccomplex_t IT_0112 = s_23 + IT_0111;
    const ccomplex_t IT_0113 = IT_0000*IT_0112;
    const ccomplex_t IT_0114 = (-4)*IT_0113;
    const ccomplex_t IT_0115 = 2*IT_0113;
    const ccomplex_t IT_0116 = -conj(IT_0104);
    const ccomplex_t IT_0117 = conj(IT_0105) + IT_0116;
    const ccomplex_t IT_0118 = pow(m_tau, 2);
    const ccomplex_t IT_0119 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0118 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0120 = conj(N_B4)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0121 = IT_0007*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = IT_0016*IT_0027;
    const ccomplex_t IT_0124 = conj(N_d4)*e_em*m_tau*IT_0025*conj(U_stau_01);
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0122 + 0.5*IT_0126);
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = N_d4*e_em*m_tau*IT_0025*U_stau_01;
    const ccomplex_t IT_0130 = IT_0123*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0133 = IT_0007*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + 2*IT_0134);
    const ccomplex_t IT_0136 = 0.5*IT_0135;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0128*IT_0137;
    const ccomplex_t IT_0139 = IT_0119*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0118 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0143 = conj(N_B4)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0144 = IT_0007*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_d4)*e_em*m_tau*IT_0025*conj(U_stau_00);
    const ccomplex_t IT_0147 = IT_0123*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + 0.5*IT_0148);
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0152 = IT_0007*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = N_d4*e_em*m_tau*IT_0025*U_stau_00;
    const ccomplex_t IT_0155 = IT_0123*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + 0.5*IT_0156);
    const ccomplex_t IT_0158 = 2*IT_0157;
    const ccomplex_t IT_0159 = 0.5*IT_0158;
    const ccomplex_t IT_0160 = 2*IT_0159;
    const ccomplex_t IT_0161 = 0.5*IT_0160;
    const ccomplex_t IT_0162 = 2*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0162;
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0150*IT_0164;
    const ccomplex_t IT_0166 = IT_0142*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = IT_0141 + IT_0168;
    const ccomplex_t IT_0170 = -IT_0104;
    const ccomplex_t IT_0171 = IT_0105 + IT_0170;
    const ccomplex_t IT_0172 = conj(IT_0141) + conj(IT_0168);
    const ccomplex_t IT_0173 = -conj(IT_0046);
    const ccomplex_t IT_0174 = conj(IT_0047) + IT_0173;
    const ccomplex_t IT_0175 = N_B4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0176 = IT_0007*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = N_W4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0179 = IT_0016*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = N_d4*e_em*m_tau*IT_0025*conj(U_stau_11);
    const ccomplex_t IT_0182 = IT_0123*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*(IT_0177 + IT_0180 + 
      -IT_0183);
    const ccomplex_t IT_0185 = (-0.5)*IT_0184;
    const ccomplex_t IT_0186 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0187 = IT_0007*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0190 = IT_0016*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_d4)*e_em*m_tau*IT_0025*U_stau_11;
    const ccomplex_t IT_0193 = IT_0123*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*(IT_0188 + IT_0191 + 
      -IT_0194);
    const ccomplex_t IT_0196 = (-0.5)*IT_0195;
    const ccomplex_t IT_0197 = IT_0185*IT_0196;
    const ccomplex_t IT_0198 = IT_0119*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = (-0.5)*IT_0199;
    const ccomplex_t IT_0201 = N_B4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0202 = IT_0007*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = N_W4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0205 = IT_0016*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = N_d4*e_em*m_tau*IT_0025*conj(U_stau_10);
    const ccomplex_t IT_0208 = IT_0123*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(IT_0203 + IT_0206 + 
      -IT_0209);
    const ccomplex_t IT_0211 = (-0.5)*IT_0210;
    const ccomplex_t IT_0212 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0213 = IT_0007*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0216 = IT_0016*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = conj(N_d4)*e_em*m_tau*IT_0025*U_stau_10;
    const ccomplex_t IT_0219 = IT_0123*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0214 + IT_0217 + 
      -IT_0220);
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = IT_0211*IT_0222;
    const ccomplex_t IT_0224 = IT_0142*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = IT_0200 + IT_0226;
    const ccomplex_t IT_0228 = -IT_0046;
    const ccomplex_t IT_0229 = IT_0047 + IT_0228;
    const ccomplex_t IT_0230 = conj(IT_0200) + conj(IT_0226);
    const ccomplex_t IT_0231 = 4*IT_0113;
    const ccomplex_t IT_0232 = (-0.5)*IT_0084;
    const ccomplex_t IT_0233 = (-0.5)*IT_0082;
    const ccomplex_t IT_0234 = (-0.5)*IT_0075;
    const ccomplex_t IT_0235 = (-0.5)*IT_0067;
    const ccomplex_t IT_0236 = IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = conj(IT_0234) + conj(IT_0235);
    const ccomplex_t IT_0238 = (-2)*IT_0113;
    const ccomplex_t IT_0239 = cpow(s_13 + (-0.5)*IT_0030 + (-0.5)*IT_0118 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0240 = IT_0223*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*IT_0240;
    const ccomplex_t IT_0242 = cpow(s_13 + (-0.5)*IT_0030 + (-0.5)*IT_0118 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0243 = IT_0197*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0241 + 0.5*IT_0244;
    const ccomplex_t IT_0246 = IT_0138*IT_0242;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = 0.5*IT_0247;
    const ccomplex_t IT_0249 = IT_0165*IT_0239;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = 0.5*IT_0250;
    const ccomplex_t IT_0252 = IT_0248 + IT_0251;
    const ccomplex_t IT_0253 = conj(IT_0248) + conj(IT_0251);
    const ccomplex_t IT_0254 = 2*IT_0003;
    const ccomplex_t IT_0255 = s_12*IT_0118;
    const ccomplex_t IT_0256 = s_14*s_23;
    const ccomplex_t IT_0257 = s_13*s_24;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = IT_0256 + IT_0258;
    const ccomplex_t IT_0260 = s_12*s_34;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0259 + IT_0261;
    const ccomplex_t IT_0263 = IT_0255 + IT_0262;
    const ccomplex_t IT_0264 = IT_0164*IT_0211;
    const ccomplex_t IT_0265 = IT_0239*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*IT_0265;
    const ccomplex_t IT_0267 = IT_0137*IT_0185;
    const ccomplex_t IT_0268 = IT_0242*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = (-0.5)*IT_0266 + (-0.5)*IT_0269;
    const ccomplex_t IT_0271 = 2*IT_0047;
    const ccomplex_t IT_0272 = (-2)*IT_0046;
    const ccomplex_t IT_0273 = IT_0271 + IT_0272;
    const ccomplex_t IT_0274 = 2*conj(IT_0047);
    const ccomplex_t IT_0275 = (-2)*conj(IT_0046);
    const ccomplex_t IT_0276 = IT_0274 + IT_0275;
    const ccomplex_t IT_0277 = IT_0150*IT_0222;
    const ccomplex_t IT_0278 = IT_0239*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0128*IT_0196;
    const ccomplex_t IT_0281 = IT_0242*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = (-0.5)*IT_0279 + (-0.5)*IT_0282;
    const ccomplex_t IT_0284 = 2*IT_0105;
    const ccomplex_t IT_0285 = (-2)*IT_0104;
    const ccomplex_t IT_0286 = IT_0284 + IT_0285;
    const ccomplex_t IT_0287 = 2*conj(IT_0105);
    const ccomplex_t IT_0288 = (-2)*conj(IT_0104);
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = s_13 + s_14;
    const ccomplex_t IT_0291 = IT_0000*IT_0290;
    const ccomplex_t IT_0292 = IT_0007*IT_0026;
    const ccomplex_t IT_0293 = IT_0013*IT_0292;
    const ccomplex_t IT_0294 = sin(beta);
    const ccomplex_t IT_0295 = IT_0007*IT_0294;
    const ccomplex_t IT_0296 = IT_0009*IT_0295;
    const ccomplex_t IT_0297 = IT_0016*IT_0026;
    const ccomplex_t IT_0298 = IT_0021*IT_0297;
    const ccomplex_t IT_0299 = IT_0016*IT_0294;
    const ccomplex_t IT_0300 = IT_0018*IT_0299;
    const ccomplex_t IT_0301 = -IT_0293 + -IT_0296 + IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = e_em*m_tau*IT_0016*IT_0025;
    const ccomplex_t IT_0303 = 0.5*IT_0302;
    const ccomplex_t IT_0304 = pow(m_Z, 2);
    const ccomplex_t IT_0305 = cpow((-2)*s_12 + (-2)*IT_0030 + IT_0304 + 
      -reg_prop, -1);
    const ccomplex_t IT_0306 = IT_0303*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*IT_0306;
    const ccomplex_t IT_0308 = IT_0301*IT_0307;
    const ccomplex_t IT_0309 = IT_0009*IT_0292;
    const ccomplex_t IT_0310 = IT_0013*IT_0295;
    const ccomplex_t IT_0311 = IT_0018*IT_0297;
    const ccomplex_t IT_0312 = IT_0021*IT_0299;
    const ccomplex_t IT_0313 = -IT_0309 + IT_0310 + IT_0311 + -IT_0312;
    const ccomplex_t IT_0314 = e_em*m_tau*IT_0016*IT_0025*IT_0027*IT_0294;
    const ccomplex_t IT_0315 = (-0.5)*IT_0314;
    const ccomplex_t IT_0316 = cpow((-2)*s_12 + (-2)*IT_0030 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0317 = IT_0315*IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*IT_0317;
    const ccomplex_t IT_0319 = IT_0313*IT_0318;
    const ccomplex_t IT_0320 = (-0.5)*IT_0308 + (-0.5)*IT_0319;
    const ccomplex_t IT_0321 = 2*IT_0141;
    const ccomplex_t IT_0322 = 2*IT_0168;
    const ccomplex_t IT_0323 = IT_0321 + IT_0322;
    const ccomplex_t IT_0324 = (-2)*IT_0245;
    const ccomplex_t IT_0325 = IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = 2*conj(IT_0141);
    const ccomplex_t IT_0327 = 2*conj(IT_0168);
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = (-2)*conj(IT_0245);
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = IT_0089*IT_0292;
    const ccomplex_t IT_0332 = IT_0087*IT_0295;
    const ccomplex_t IT_0333 = IT_0093*IT_0297;
    const ccomplex_t IT_0334 = IT_0091*IT_0299;
    const ccomplex_t IT_0335 = IT_0331 + IT_0332 + -IT_0333 + -IT_0334;
    const ccomplex_t IT_0336 = IT_0307*IT_0335;
    const ccomplex_t IT_0337 = IT_0087*IT_0292;
    const ccomplex_t IT_0338 = IT_0089*IT_0295;
    const ccomplex_t IT_0339 = IT_0091*IT_0297;
    const ccomplex_t IT_0340 = IT_0093*IT_0299;
    const ccomplex_t IT_0341 = IT_0337 + -IT_0338 + -IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = IT_0318*IT_0341;
    const ccomplex_t IT_0343 = 0.5*IT_0336 + 0.5*IT_0342;
    const ccomplex_t IT_0344 = 2*IT_0200;
    const ccomplex_t IT_0345 = 2*IT_0226;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = (-2)*IT_0248;
    const ccomplex_t IT_0348 = (-2)*IT_0251;
    const ccomplex_t IT_0349 = IT_0346 + IT_0347 + IT_0348;
    const ccomplex_t IT_0350 = 2*conj(IT_0200);
    const ccomplex_t IT_0351 = 2*conj(IT_0226);
    const ccomplex_t IT_0352 = IT_0350 + IT_0351;
    const ccomplex_t IT_0353 = (-2)*conj(IT_0248);
    const ccomplex_t IT_0354 = (-2)*conj(IT_0251);
    const ccomplex_t IT_0355 = IT_0352 + IT_0353 + IT_0354;
    const ccomplex_t IT_0356 = (-0.5)*IT_0336 + (-0.5)*IT_0342;
    const ccomplex_t IT_0357 = 2*IT_0248;
    const ccomplex_t IT_0358 = 2*IT_0251;
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = (-2)*IT_0200;
    const ccomplex_t IT_0361 = (-2)*IT_0226;
    const ccomplex_t IT_0362 = IT_0359 + IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = 2*conj(IT_0248);
    const ccomplex_t IT_0364 = 2*conj(IT_0251);
    const ccomplex_t IT_0365 = IT_0363 + IT_0364;
    const ccomplex_t IT_0366 = (-2)*conj(IT_0200);
    const ccomplex_t IT_0367 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0368 = IT_0365 + IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = 0.5*IT_0308 + 0.5*IT_0319;
    const ccomplex_t IT_0370 = 2*IT_0245;
    const ccomplex_t IT_0371 = (-2)*IT_0141;
    const ccomplex_t IT_0372 = (-2)*IT_0168;
    const ccomplex_t IT_0373 = IT_0370 + IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = 2*conj(IT_0245);
    const ccomplex_t IT_0375 = (-2)*conj(IT_0141);
    const ccomplex_t IT_0376 = (-2)*conj(IT_0168);
    const ccomplex_t IT_0377 = IT_0374 + IT_0375 + IT_0376;
    const ccomplex_t IT_0378 = conj(IT_0068) + conj(IT_0232);
    const ccomplex_t IT_0379 = 4*IT_0320;
    const ccomplex_t IT_0380 = (-4)*IT_0369;
    const ccomplex_t IT_0381 = IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = IT_0068 + IT_0232;
    const ccomplex_t IT_0383 = 4*conj(IT_0320);
    const ccomplex_t IT_0384 = (-4)*conj(IT_0369);
    const ccomplex_t IT_0385 = IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = conj(IT_0076) + conj(IT_0233);
    const ccomplex_t IT_0387 = 4*IT_0369;
    const ccomplex_t IT_0388 = (-4)*IT_0320;
    const ccomplex_t IT_0389 = IT_0387 + IT_0388;
    const ccomplex_t IT_0390 = IT_0076 + IT_0233;
    const ccomplex_t IT_0391 = 4*conj(IT_0369);
    const ccomplex_t IT_0392 = (-4)*conj(IT_0320);
    const ccomplex_t IT_0393 = IT_0391 + IT_0392;
    const ccomplex_t IT_0394 = conj(IT_0085) + conj(IT_0235);
    const ccomplex_t IT_0395 = 4*IT_0356;
    const ccomplex_t IT_0396 = (-4)*IT_0343;
    const ccomplex_t IT_0397 = IT_0395 + IT_0396;
    const ccomplex_t IT_0398 = IT_0085 + IT_0235;
    const ccomplex_t IT_0399 = 4*conj(IT_0356);
    const ccomplex_t IT_0400 = (-4)*conj(IT_0343);
    const ccomplex_t IT_0401 = IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = conj(IT_0083) + conj(IT_0234);
    const ccomplex_t IT_0403 = 4*IT_0343;
    const ccomplex_t IT_0404 = (-4)*IT_0356;
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = IT_0083 + IT_0234;
    const ccomplex_t IT_0407 = 4*conj(IT_0343);
    const ccomplex_t IT_0408 = (-4)*conj(IT_0356);
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = conj(IT_0320)*IT_0325 + IT_0320*IT_0330 + conj
      (IT_0343)*IT_0349 + IT_0343*IT_0355 + conj(IT_0356)*IT_0362 + IT_0356
      *IT_0368 + conj(IT_0369)*IT_0373 + IT_0369*IT_0377 + IT_0378*IT_0381 +
       IT_0382*IT_0385 + IT_0386*IT_0389 + IT_0390*IT_0393 + IT_0394*IT_0397 +
       IT_0398*IT_0401 + IT_0402*IT_0405 + IT_0406*IT_0409;
    const ccomplex_t IT_0411 = s_23 + s_24;
    const ccomplex_t IT_0412 = IT_0000*IT_0411;
    const ccomplex_t IT_0413 = IT_0330*IT_0343 + IT_0325*conj(IT_0343) + conj
      (IT_0320)*IT_0349 + IT_0320*IT_0355 + IT_0368*IT_0369 + IT_0362*conj
      (IT_0369) + conj(IT_0356)*IT_0373 + IT_0356*IT_0377 + IT_0389*IT_0394 +
       IT_0386*IT_0397 + IT_0393*IT_0398 + IT_0390*IT_0401 + IT_0381*IT_0402 +
       IT_0378*IT_0405 + IT_0385*IT_0406 + IT_0382*IT_0409;
    const ccomplex_t IT_0414 = IT_0119*IT_0267;
    const ccomplex_t IT_0415 = (0 + _Complex_I*1)*IT_0414;
    const ccomplex_t IT_0416 = IT_0142*IT_0264;
    const ccomplex_t IT_0417 = (0 + _Complex_I*1)*IT_0416;
    const ccomplex_t IT_0418 = 0.5*IT_0415 + 0.5*IT_0417;
    const ccomplex_t IT_0419 = IT_0142*IT_0277;
    const ccomplex_t IT_0420 = (0 + _Complex_I*1)*IT_0419;
    const ccomplex_t IT_0421 = IT_0119*IT_0280;
    const ccomplex_t IT_0422 = (0 + _Complex_I*1)*IT_0421;
    const ccomplex_t IT_0423 = 0.5*IT_0420 + 0.5*IT_0422;
    const ccomplex_t IT_0424 = -IT_0255;
    const ccomplex_t IT_0425 = IT_0256 + IT_0257 + IT_0261;
    const ccomplex_t IT_0426 = 2*IT_0270;
    const ccomplex_t IT_0427 = 2*conj(IT_0270);
    const ccomplex_t IT_0428 = 2*IT_0283;
    const ccomplex_t IT_0429 = 2*conj(IT_0283);
    const ccomplex_t IT_0430 = s_23*m_N_4*m_tau;
    const ccomplex_t IT_0431 = 4*IT_0168;
    const ccomplex_t IT_0432 = 0.125*IT_0431;
    const ccomplex_t IT_0433 = 4*IT_0141;
    const ccomplex_t IT_0434 = 0.125*IT_0433;
    const ccomplex_t IT_0435 = IT_0382 + IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = 4*conj(IT_0168);
    const ccomplex_t IT_0437 = 0.125*IT_0436;
    const ccomplex_t IT_0438 = 4*conj(IT_0141);
    const ccomplex_t IT_0439 = 0.125*IT_0438;
    const ccomplex_t IT_0440 = IT_0378 + IT_0437 + IT_0439;
    const ccomplex_t IT_0441 = 4*IT_0270;
    const ccomplex_t IT_0442 = 0.125*IT_0441;
    const ccomplex_t IT_0443 = IT_0418 + IT_0442;
    const ccomplex_t IT_0444 = 4*conj(IT_0270);
    const ccomplex_t IT_0445 = 0.125*IT_0444;
    const ccomplex_t IT_0446 = conj(IT_0418) + IT_0445;
    const ccomplex_t IT_0447 = 4*IT_0418;
    const ccomplex_t IT_0448 = IT_0426 + IT_0447;
    const ccomplex_t IT_0449 = 4*conj(IT_0418);
    const ccomplex_t IT_0450 = IT_0427 + IT_0449;
    const ccomplex_t IT_0451 = 4*IT_0283;
    const ccomplex_t IT_0452 = 4*conj(IT_0283);
    const ccomplex_t IT_0453 = 0.125*conj(IT_0283)*IT_0323 + 0.125*IT_0283
      *IT_0328 + conj(IT_0423)*IT_0435 + IT_0423*IT_0440 + IT_0402*IT_0443 +
       IT_0406*IT_0446 + 0.125*IT_0230*IT_0448 + 0.125*IT_0227*IT_0450 + 0.125
      *IT_0378*IT_0451 + 0.125*IT_0382*IT_0452;
    const ccomplex_t IT_0454 = s_24*m_N_4*m_tau;
    const ccomplex_t IT_0455 = 4*IT_0245;
    const ccomplex_t IT_0456 = 0.125*IT_0455;
    const ccomplex_t IT_0457 = IT_0390 + IT_0456;
    const ccomplex_t IT_0458 = 0.125*IT_0449;
    const ccomplex_t IT_0459 = conj(IT_0270) + IT_0458;
    const ccomplex_t IT_0460 = 0.125*IT_0447;
    const ccomplex_t IT_0461 = IT_0270 + IT_0460;
    const ccomplex_t IT_0462 = 2*IT_0418;
    const ccomplex_t IT_0463 = IT_0441 + IT_0462;
    const ccomplex_t IT_0464 = 2*conj(IT_0418);
    const ccomplex_t IT_0465 = IT_0444 + IT_0464;
    const ccomplex_t IT_0466 = 4*IT_0423;
    const ccomplex_t IT_0467 = 4*conj(IT_0423);
    const ccomplex_t IT_0468 = 4*conj(IT_0245);
    const ccomplex_t IT_0469 = 0.125*IT_0374*IT_0423 + 0.125*IT_0370*conj
      (IT_0423) + conj(IT_0283)*IT_0457 + IT_0398*IT_0459 + IT_0394*IT_0461 +
       0.125*IT_0253*IT_0463 + 0.125*IT_0252*IT_0465 + 0.125*IT_0386*IT_0466 +
       0.125*IT_0390*IT_0467 + IT_0283*(IT_0386 + 0.125*IT_0468);
    const ccomplex_t IT_0470 = -IT_0118;
    const ccomplex_t IT_0471 = s_34 + IT_0470;
    const ccomplex_t IT_0472 = s_12*IT_0471;
    const ccomplex_t IT_0473 = (-8)*conj(IT_0046);
    const ccomplex_t IT_0474 = IT_0262 + IT_0424;
    const ccomplex_t IT_0475 = (-2)*IT_0270;
    const ccomplex_t IT_0476 = (-2)*conj(IT_0270);
    const ccomplex_t IT_0477 = (-2)*IT_0283;
    const ccomplex_t IT_0478 = (-2)*conj(IT_0283);
    const ccomplex_t IT_0479 = IT_0255 + IT_0259 + IT_0260;
    const ccomplex_t IT_0480 = 2*IT_0423;
    const ccomplex_t IT_0481 = 2*conj(IT_0423);
    const ccomplex_t IT_0482 = (-2)*IT_0418;
    const ccomplex_t IT_0483 = (-2)*conj(IT_0418);
    const ccomplex_t IT_0484 = (-2)*IT_0423;
    const ccomplex_t IT_0485 = (-2)*conj(IT_0423);
    const ccomplex_t IT_0486 = (-2)*IT_0003;
    const ccomplex_t IT_0487 = (-4)*IT_0003;
    const ccomplex_t IT_0488 = 16*IT_0232;
    const ccomplex_t IT_0489 = 0.125*IT_0488;
    const ccomplex_t IT_0490 = IT_0169 + IT_0489;
    const ccomplex_t IT_0491 = 16*IT_0083;
    const ccomplex_t IT_0492 = 0.125*IT_0491;
    const ccomplex_t IT_0493 = IT_0227 + IT_0492;
    const ccomplex_t IT_0494 = 16*IT_0234;
    const ccomplex_t IT_0495 = 0.125*IT_0494;
    const ccomplex_t IT_0496 = IT_0493 + IT_0495;
    const ccomplex_t IT_0497 = IT_0232 + IT_0434;
    const ccomplex_t IT_0498 = IT_0432 + IT_0497;
    const ccomplex_t IT_0499 = 4*conj(IT_0200);
    const ccomplex_t IT_0500 = 16*conj(IT_0083);
    const ccomplex_t IT_0501 = 16*conj(IT_0232);
    const ccomplex_t IT_0502 = 4*IT_0200;
    const ccomplex_t IT_0503 = 0.125*IT_0502;
    const ccomplex_t IT_0504 = 4*IT_0226;
    const ccomplex_t IT_0505 = 0.125*IT_0504;
    const ccomplex_t IT_0506 = 0.125*IT_0168*IT_0438 + 0.125*conj(IT_0418)
      *IT_0447 + 0.125*conj(IT_0423)*IT_0466 + IT_0378*IT_0490 + conj(IT_0083)
      *IT_0493 + conj(IT_0234)*IT_0496 + conj(IT_0141)*IT_0497 + conj(IT_0168)
      *IT_0498 + 0.125*IT_0226*IT_0499 + 0.125*IT_0234*IT_0500 + 0.125*IT_0068*
      (16*conj(IT_0068) + 8*conj(IT_0141) + 8*conj(IT_0168) + IT_0501) + conj
      (IT_0200)*(IT_0406 + IT_0503) + conj(IT_0226)*(IT_0406 + IT_0503 + IT_0505);
    const ccomplex_t IT_0507 = s_34*IT_0030;
    const ccomplex_t IT_0508 = 4*conj(IT_0226);
    const ccomplex_t IT_0509 = 0.125*IT_0508;
    const ccomplex_t IT_0510 = 0.125*IT_0499;
    const ccomplex_t IT_0511 = IT_0402 + IT_0509 + IT_0510;
    const ccomplex_t IT_0512 = 4*IT_0248;
    const ccomplex_t IT_0513 = 0.125*IT_0512;
    const ccomplex_t IT_0514 = IT_0398 + IT_0513;
    const ccomplex_t IT_0515 = 4*IT_0251;
    const ccomplex_t IT_0516 = 0.125*IT_0515;
    const ccomplex_t IT_0517 = IT_0514 + IT_0516;
    const ccomplex_t IT_0518 = 4*conj(IT_0251);
    const ccomplex_t IT_0519 = 0.125*IT_0518;
    const ccomplex_t IT_0520 = 4*conj(IT_0248);
    const ccomplex_t IT_0521 = 0.125*IT_0520;
    const ccomplex_t IT_0522 = IT_0394 + IT_0519 + IT_0521;
    const ccomplex_t IT_0523 = 0.125*conj(IT_0251)*IT_0321 + 0.125*IT_0251
      *IT_0326 + 0.125*conj(IT_0245)*IT_0344 + 0.125*IT_0245*IT_0350 + 0.125
      *IT_0172*IT_0357 + 0.125*conj(IT_0168)*IT_0358 + 0.125*IT_0169*IT_0363 +
       0.125*IT_0168*IT_0364 + 0.125*conj(IT_0226)*IT_0370 + 0.125*IT_0226
      *IT_0374 + 0.125*IT_0394*(IT_0431 + IT_0433) + 0.125*IT_0398*(IT_0436 +
       IT_0438) + 0.125*IT_0402*IT_0455 + 0.125*IT_0406*IT_0468 + IT_0386*
      (IT_0406 + IT_0503 + IT_0505) + IT_0390*IT_0511 + IT_0378*IT_0517 +
       IT_0382*IT_0522;
    const ccomplex_t IT_0524 = s_34 + IT_0118;
    const ccomplex_t IT_0525 = s_12*IT_0524;
    const ccomplex_t IT_0526 = (-8)*IT_0320;
    const ccomplex_t IT_0527 = (-8)*conj(IT_0320);
    const ccomplex_t IT_0528 = -conj(IT_0343);
    const ccomplex_t IT_0529 = conj(IT_0356) + IT_0528;
    const ccomplex_t IT_0530 = -IT_0343;
    const ccomplex_t IT_0531 = IT_0030*IT_0524;
    const ccomplex_t IT_0532 = conj(IT_0343)*(IT_0462 + IT_0475) + IT_0343*
      (IT_0464 + IT_0476) + 8*conj(IT_0320)*(IT_0343 + 0.125*IT_0477 + 0.125
      *IT_0480) + 8*IT_0320*(conj(IT_0343) + 0.125*IT_0478 + 0.125*IT_0481) +
       conj(IT_0356)*(IT_0426 + IT_0482 + IT_0526) + IT_0356*(IT_0427 + IT_0483 
      + IT_0527) + 8*IT_0369*(0.125*IT_0429 + 0.125*IT_0485 + IT_0529) + 8*conj
      (IT_0369)*(IT_0356 + 0.125*IT_0428 + 0.125*IT_0484 + IT_0530);
    const ccomplex_t IT_0533 = IT_0030*IT_0471;
    const ccomplex_t IT_0534 = 0.125*IT_0283;
    const ccomplex_t IT_0535 = 0.125*conj(IT_0283);
    const ccomplex_t IT_0536 = (-0.5)*conj(IT_0270)*IT_0286 + (-0.5)*IT_0270
      *IT_0289 + -IT_0110*IT_0418 + IT_0105*conj(IT_0418) + (-4)*conj(IT_0046)*
      (IT_0105 + 0.25*IT_0423) + (-4)*conj(IT_0047)*(IT_0107 + -1./4*IT_0423) +
       4*(IT_0046 + -IT_0047)*(IT_0110 + -1./4*conj(IT_0423)) + -IT_0104*(conj
      (IT_0418) + 0.5*IT_0473) + (-4)*IT_0276*IT_0534 + (-4)*IT_0273*IT_0535;
    const ccomplex_t IT_0537 = s_14*m_N_4*m_tau;
    const ccomplex_t IT_0538 = IT_0378*IT_0443 + IT_0382*IT_0446 + 0.125
      *IT_0172*IT_0448 + 0.125*IT_0169*IT_0450 + 0.125*IT_0402*IT_0451 + 0.125
      *IT_0406*IT_0452 + conj(IT_0423)*(IT_0406 + IT_0503 + IT_0505) + IT_0423
      *IT_0511 + IT_0352*IT_0534 + IT_0346*IT_0535;
    const ccomplex_t IT_0539 = s_13*m_N_4*m_tau;
    const ccomplex_t IT_0540 = 0.125*IT_0423;
    const ccomplex_t IT_0541 = 0.125*conj(IT_0423);
    const ccomplex_t IT_0542 = 0.125*IT_0245;
    const ccomplex_t IT_0543 = 0.125*conj(IT_0245);
    const ccomplex_t IT_0544 = 0.125*IT_0398;
    const ccomplex_t IT_0545 = 0.125*IT_0394;
    const ccomplex_t IT_0546 = 8*IT_0539*(IT_0390*IT_0459 + IT_0386*IT_0461 +
       conj(IT_0283)*IT_0517 + IT_0283*IT_0522 + IT_0365*IT_0540 + IT_0359
      *IT_0541 + IT_0465*IT_0542 + IT_0463*IT_0543 + IT_0467*IT_0544 + IT_0466
      *IT_0545);
    const ccomplex_t IT_0547 = 0.125*IT_0230*IT_0359 + 0.125*IT_0227*IT_0365 +
       IT_0386*IT_0435 + IT_0390*IT_0440 + 0.125*IT_0378*IT_0455 + 0.125*IT_0382
      *IT_0468 + IT_0402*IT_0517 + IT_0406*IT_0522 + IT_0328*IT_0542 + IT_0323
      *IT_0543 + (IT_0499 + IT_0508)*IT_0544 + (IT_0502 + IT_0504)*IT_0545;
    const ccomplex_t IT_0548 = IT_0030*IT_0118;
    const ccomplex_t IT_0549 = 16*conj(IT_0085);
    const ccomplex_t IT_0550 = 16*IT_0076;
    const ccomplex_t IT_0551 = 16*conj(IT_0233);
    const ccomplex_t IT_0552 = 16*IT_0233;
    const ccomplex_t IT_0553 = 2*IT_0235;
    const ccomplex_t IT_0554 = 16*IT_0085;
    const ccomplex_t IT_0555 = 0.125*IT_0554;
    const ccomplex_t IT_0556 = conj(IT_0083)*IT_0169 + IT_0083*IT_0172 + conj
      (IT_0233)*IT_0252 + IT_0233*IT_0253 + IT_0227*(conj(IT_0232) + IT_0437 +
       IT_0439) + conj(IT_0234)*IT_0490 + 0.125*conj(IT_0232)*IT_0491 + conj
      (IT_0068)*IT_0496 + IT_0230*IT_0498 + 0.125*IT_0232*IT_0500 + 2*IT_0068*
      (0.5*IT_0230 + conj(IT_0234) + 0.0625*IT_0500) + IT_0234*(IT_0172 + 0.125
      *IT_0501) + conj(IT_0245)*IT_0517 + IT_0245*IT_0522 + IT_0465*IT_0534 +
       IT_0463*IT_0535 + IT_0450*IT_0540 + IT_0448*IT_0541 + IT_0076*(IT_0253 +
       0.125*IT_0549) + 0.125*conj(IT_0235)*IT_0550 + IT_0544*IT_0551 + IT_0545
      *IT_0552 + conj(IT_0076)*(IT_0252 + IT_0553 + IT_0555);
    const ccomplex_t IT_0557 = 0.125*IT_0552;
    const ccomplex_t IT_0558 = IT_0245 + IT_0557;
    const ccomplex_t IT_0559 = 0.125*conj(IT_0270)*IT_0441 + 0.125*conj
      (IT_0283)*IT_0451 + conj(IT_0245)*IT_0457 + conj(IT_0248)*IT_0514 + conj
      (IT_0251)*IT_0517 + 0.125*IT_0251*IT_0520 + 0.125*IT_0235*IT_0549 + 0.125
      *IT_0076*IT_0551 + conj(IT_0085)*(IT_0252 + IT_0555) + conj(IT_0235)*
      (IT_0252 + IT_0553 + IT_0555) + conj(IT_0233)*IT_0558 + conj(IT_0076)*
      (0.125*IT_0550 + IT_0558);
    const ccomplex_t IT_0560 = IT_0004*(IT_0049*IT_0077 + IT_0079*IT_0080 +
       IT_0086*IT_0107 + IT_0108*IT_0110) + (IT_0079*IT_0086 + IT_0080*IT_0107 +
       IT_0049*IT_0108 + IT_0077*IT_0110)*IT_0114 + IT_0115*(IT_0117*IT_0169 +
       IT_0171*IT_0172 + IT_0174*IT_0227 + IT_0229*IT_0230) + IT_0231*(IT_0117*
      (IT_0232 + IT_0233) + IT_0171*(conj(IT_0232) + conj(IT_0233)) + IT_0174
      *IT_0236 + IT_0229*IT_0237) + IT_0238*(IT_0110*IT_0245 + IT_0107*conj
      (IT_0245) + IT_0049*IT_0252 + IT_0079*IT_0253) + (IT_0049*IT_0245 +
       IT_0079*conj(IT_0245) + IT_0110*IT_0252 + IT_0107*IT_0253)*IT_0254 +
       IT_0263*(conj(IT_0270)*IT_0273 + IT_0270*IT_0276 + conj(IT_0283)*IT_0286 
      + IT_0283*IT_0289) + IT_0291*IT_0410 + IT_0412*IT_0413 + (IT_0276*IT_0418 
      + IT_0273*conj(IT_0418) + IT_0289*IT_0423 + IT_0286*conj(IT_0423))*
      (IT_0259 + IT_0260 + IT_0424) + IT_0425*(conj(IT_0418)*IT_0426 + IT_0418
      *IT_0427 + conj(IT_0423)*IT_0428 + IT_0423*IT_0429) + 8*IT_0430*IT_0453 +
       8*IT_0454*IT_0469 + 8*IT_0472*(IT_0046*(conj(IT_0046) + -conj(IT_0047)) +
       IT_0107*IT_0110 + IT_0047*(conj(IT_0047) + 0.125*IT_0473)) + IT_0474*
      (conj(IT_0320)*IT_0426 + IT_0320*IT_0427 + conj(IT_0343)*IT_0428 + IT_0343
      *IT_0429 + conj(IT_0369)*IT_0475 + IT_0369*IT_0476 + conj(IT_0356)*IT_0477
       + IT_0356*IT_0478) + IT_0479*(conj(IT_0320)*IT_0462 + IT_0320*IT_0464 +
       conj(IT_0343)*IT_0480 + IT_0343*IT_0481 + conj(IT_0369)*IT_0482 + IT_0369
      *IT_0483 + conj(IT_0356)*IT_0484 + IT_0356*IT_0485) + (IT_0169*IT_0174 + 
      (conj(IT_0105) + IT_0116)*IT_0227 + IT_0172*IT_0229 + IT_0171*IT_0230)
      *IT_0486 + (IT_0174*(IT_0232 + IT_0233) + IT_0229*(conj(IT_0232) + conj
      (IT_0233)) + (conj(IT_0105) + IT_0116)*IT_0236 + IT_0171*IT_0237)*IT_0487 
      + 8*IT_0256*IT_0506 + 8*IT_0507*IT_0523 + 8*IT_0525*(IT_0320*conj(IT_0320)
       + 0.125*conj(IT_0369)*IT_0526 + IT_0369*(conj(IT_0369) + 0.125*IT_0527) +
       IT_0529*(IT_0356 + IT_0530)) + IT_0531*IT_0532 + (-2)*IT_0533*IT_0536 + 8
      *IT_0537*IT_0538 + IT_0546 + 8*IT_0255*IT_0547 + 8*IT_0548*IT_0556 + 8
      *IT_0257*IT_0559;
    return create_ccomplex_return(IT_0560);
}

