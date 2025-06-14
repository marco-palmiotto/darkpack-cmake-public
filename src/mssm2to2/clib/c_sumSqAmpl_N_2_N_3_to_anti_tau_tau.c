#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_3_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_3_to_anti_tau_tau(
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
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
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
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0005 + IT_0009 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0016;
    const ccomplex_t IT_0021 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0014 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = cos(beta);
    const ccomplex_t IT_0033 = cpow(IT_0032, -1);
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0016
      *IT_0031*IT_0033;
    const ccomplex_t IT_0035 = 0.5*IT_0034;
    const ccomplex_t IT_0036 = pow(m_N_2, 2);
    const ccomplex_t IT_0037 = pow(m_N_3, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0036 + -IT_0037 + -reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0035*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0030*IT_0040;
    const ccomplex_t IT_0042 = IT_0003*IT_0008;
    const ccomplex_t IT_0043 = IT_0018*IT_0020;
    const ccomplex_t IT_0044 = IT_0003*IT_0012;
    const ccomplex_t IT_0045 = IT_0020*IT_0023;
    const ccomplex_t IT_0046 = IT_0042 + IT_0043 + IT_0044 + IT_0045;
    const ccomplex_t IT_0047 = IT_0004*IT_0007;
    const ccomplex_t IT_0048 = IT_0017*IT_0021;
    const ccomplex_t IT_0049 = IT_0007*IT_0010;
    const ccomplex_t IT_0050 = IT_0017*IT_0025;
    const ccomplex_t IT_0051 = -IT_0047 + -IT_0048 + -IT_0049 + -IT_0050;
    const ccomplex_t IT_0052 = IT_0046 + IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*e_em*m_tau*IT_0000*IT_0016
      *IT_0031*IT_0033;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = cpow((-2)*s_12 + -IT_0036 + -IT_0037 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0056*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0054*IT_0059;
    const ccomplex_t IT_0061 = -IT_0041 + -IT_0060;
    const ccomplex_t IT_0062 = pow(m_tau, 2);
    const ccomplex_t IT_0063 = cpow(s_13 + (-0.5)*IT_0036 + (-0.5)*IT_0062 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0064 = conj(N_B2)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0065 = IT_0002*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = IT_0016*IT_0033;
    const ccomplex_t IT_0068 = conj(N_d2)*e_em*m_tau*IT_0031*conj(U_stau_01);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + 0.5*IT_0070);
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0074 = IT_0002*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_d3*e_em*m_tau*IT_0031*U_stau_01;
    const ccomplex_t IT_0077 = IT_0067*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + 0.5*IT_0078);
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = IT_0072*IT_0080;
    const ccomplex_t IT_0082 = IT_0063*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = cpow(s_13 + (-0.5)*IT_0036 + (-0.5)*IT_0062 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0085 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0086 = IT_0002*IT_0085;
    const ccomplex_t IT_0087 = 1.4142135623731*IT_0086;
    const ccomplex_t IT_0088 = N_d3*e_em*m_tau*IT_0031*U_stau_00;
    const ccomplex_t IT_0089 = IT_0067*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + 0.5*IT_0090);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0094 = IT_0002*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_d2)*e_em*m_tau*IT_0031*conj(U_stau_00);
    const ccomplex_t IT_0097 = IT_0067*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*(IT_0095 + 0.5*IT_0098);
    const ccomplex_t IT_0100 = -IT_0099;
    const ccomplex_t IT_0101 = IT_0092*IT_0100;
    const ccomplex_t IT_0102 = IT_0084*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = 0.5*IT_0083 + 0.5*IT_0103;
    const ccomplex_t IT_0105 = m_N_3*m_tau;
    const ccomplex_t IT_0106 = -s_14;
    const ccomplex_t IT_0107 = s_13 + IT_0106;
    const ccomplex_t IT_0108 = IT_0105*IT_0107;
    const ccomplex_t IT_0109 = 2*IT_0108;
    const ccomplex_t IT_0110 = N_B2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0111 = IT_0002*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = N_W2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0114 = IT_0016*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = N_d2*e_em*m_tau*IT_0031*conj(U_stau_11);
    const ccomplex_t IT_0117 = IT_0067*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + IT_0115 + 
      -IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0122 = IT_0002*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0125 = IT_0016*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = conj(N_d3)*e_em*m_tau*IT_0031*U_stau_11;
    const ccomplex_t IT_0128 = IT_0067*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0123 + IT_0126 + 
      -IT_0129);
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0120*IT_0131;
    const ccomplex_t IT_0133 = IT_0063*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0136 = IT_0002*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0139 = IT_0016*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = conj(N_d3)*e_em*m_tau*IT_0031*U_stau_10;
    const ccomplex_t IT_0142 = IT_0067*IT_0141;
    const ccomplex_t IT_0143 = 1.4142135623731*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*(IT_0137 + IT_0140 + 
      -IT_0143);
    const ccomplex_t IT_0145 = (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0147 = IT_0002*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0150 = IT_0016*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = N_d2*e_em*m_tau*IT_0031*conj(U_stau_10);
    const ccomplex_t IT_0153 = IT_0067*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0148 + IT_0151 + 
      -IT_0154);
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = IT_0145*IT_0156;
    const ccomplex_t IT_0158 = IT_0084*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = 0.5*IT_0134 + 0.5*IT_0159;
    const ccomplex_t IT_0161 = m_N_2*m_tau;
    const ccomplex_t IT_0162 = -s_24;
    const ccomplex_t IT_0163 = s_23 + IT_0162;
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = (-2)*IT_0164;
    const ccomplex_t IT_0166 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0167 = IT_0003*IT_0166;
    const ccomplex_t IT_0168 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0169 = IT_0007*IT_0168;
    const ccomplex_t IT_0170 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0171 = IT_0003*IT_0170;
    const ccomplex_t IT_0172 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0173 = IT_0007*IT_0172;
    const ccomplex_t IT_0174 = IT_0167 + IT_0169 + IT_0171 + IT_0173;
    const ccomplex_t IT_0175 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0176 = IT_0017*IT_0175;
    const ccomplex_t IT_0177 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0178 = IT_0020*IT_0177;
    const ccomplex_t IT_0179 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0180 = IT_0017*IT_0179;
    const ccomplex_t IT_0181 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0182 = IT_0020*IT_0181;
    const ccomplex_t IT_0183 = -IT_0176 + -IT_0178 + -IT_0180 + -IT_0182;
    const ccomplex_t IT_0184 = IT_0174 + IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = IT_0040*IT_0186;
    const ccomplex_t IT_0188 = IT_0003*IT_0168;
    const ccomplex_t IT_0189 = IT_0020*IT_0175;
    const ccomplex_t IT_0190 = IT_0003*IT_0172;
    const ccomplex_t IT_0191 = IT_0020*IT_0179;
    const ccomplex_t IT_0192 = IT_0188 + IT_0189 + IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = IT_0007*IT_0166;
    const ccomplex_t IT_0194 = IT_0017*IT_0177;
    const ccomplex_t IT_0195 = IT_0007*IT_0170;
    const ccomplex_t IT_0196 = IT_0017*IT_0181;
    const ccomplex_t IT_0197 = -IT_0193 + -IT_0194 + -IT_0195 + -IT_0196;
    const ccomplex_t IT_0198 = IT_0192 + IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = 0.5*IT_0199;
    const ccomplex_t IT_0201 = IT_0059*IT_0200;
    const ccomplex_t IT_0202 = -IT_0187 + -IT_0201;
    const ccomplex_t IT_0203 = s_14*s_23;
    const ccomplex_t IT_0204 = s_13*s_24;
    const ccomplex_t IT_0205 = s_12*s_34;
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = IT_0203 + IT_0204 + IT_0206;
    const ccomplex_t IT_0208 = cpow(s_23 + (-0.5)*IT_0037 + (-0.5)*IT_0062 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0209 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0210 = IT_0002*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = N_d2*e_em*m_tau*IT_0031*U_stau_00;
    const ccomplex_t IT_0213 = IT_0067*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0211 + 0.5*IT_0214);
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0218 = IT_0002*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0221 = IT_0016*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = N_d3*e_em*m_tau*IT_0031*conj(U_stau_10);
    const ccomplex_t IT_0224 = IT_0067*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0219 + IT_0222 + 
      -IT_0225);
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = IT_0216*IT_0227;
    const ccomplex_t IT_0229 = IT_0208*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = 0.5*IT_0230;
    const ccomplex_t IT_0232 = IT_0080*IT_0120;
    const ccomplex_t IT_0233 = IT_0063*IT_0232;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1)*IT_0233;
    const ccomplex_t IT_0235 = (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = 2*IT_0235;
    const ccomplex_t IT_0237 = 2*conj(IT_0235);
    const ccomplex_t IT_0238 = IT_0100*IT_0145;
    const ccomplex_t IT_0239 = IT_0084*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = cpow(s_23 + (-0.5)*IT_0037 + (-0.5)*IT_0062 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0243 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0244 = IT_0002*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0247 = IT_0016*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = conj(N_d2)*e_em*m_tau*IT_0031*U_stau_11;
    const ccomplex_t IT_0250 = IT_0067*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*(IT_0245 + IT_0248 + 
      -IT_0251);
    const ccomplex_t IT_0253 = (-0.5)*IT_0252;
    const ccomplex_t IT_0254 = conj(N_B3)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0255 = IT_0002*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = conj(N_d3)*e_em*m_tau*IT_0031*conj(U_stau_01);
    const ccomplex_t IT_0258 = IT_0067*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*(IT_0256 + 0.5*IT_0259);
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0253*IT_0261;
    const ccomplex_t IT_0263 = IT_0242*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*IT_0263;
    const ccomplex_t IT_0265 = 0.5*IT_0264;
    const ccomplex_t IT_0266 = 2*IT_0265;
    const ccomplex_t IT_0267 = 2*conj(IT_0265);
    const ccomplex_t IT_0268 = IT_0092*IT_0156;
    const ccomplex_t IT_0269 = IT_0084*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*IT_0269;
    const ccomplex_t IT_0271 = (-0.5)*IT_0270;
    const ccomplex_t IT_0272 = conj(IT_0235) + conj(IT_0271);
    const ccomplex_t IT_0273 = N_B3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0274 = IT_0002*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = N_W3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0277 = IT_0016*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = N_d3*e_em*m_tau*IT_0031*conj(U_stau_11);
    const ccomplex_t IT_0280 = IT_0067*IT_0279;
    const ccomplex_t IT_0281 = 1.4142135623731*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*(IT_0275 + IT_0278 + 
      -IT_0281);
    const ccomplex_t IT_0283 = (-0.5)*IT_0282;
    const ccomplex_t IT_0284 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0285 = IT_0002*IT_0284;
    const ccomplex_t IT_0286 = 1.4142135623731*IT_0285;
    const ccomplex_t IT_0287 = N_d2*e_em*m_tau*IT_0031*U_stau_01;
    const ccomplex_t IT_0288 = IT_0067*IT_0287;
    const ccomplex_t IT_0289 = 1.4142135623731*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*(IT_0286 + 0.5*IT_0289);
    const ccomplex_t IT_0291 = -IT_0290;
    const ccomplex_t IT_0292 = IT_0283*IT_0291;
    const ccomplex_t IT_0293 = IT_0242*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = 0.5*IT_0294;
    const ccomplex_t IT_0296 = 2*IT_0295;
    const ccomplex_t IT_0297 = IT_0235 + IT_0271;
    const ccomplex_t IT_0298 = 2*conj(IT_0295);
    const ccomplex_t IT_0299 = 2*IT_0231;
    const ccomplex_t IT_0300 = 2*conj(IT_0231);
    const ccomplex_t IT_0301 = IT_0072*IT_0131;
    const ccomplex_t IT_0302 = IT_0063*IT_0301;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*IT_0302;
    const ccomplex_t IT_0304 = (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0304;
    const ccomplex_t IT_0306 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0307 = IT_0002*IT_0306;
    const ccomplex_t IT_0308 = 1.4142135623731*IT_0307;
    const ccomplex_t IT_0309 = conj(N_d3)*e_em*m_tau*IT_0031*conj(U_stau_00);
    const ccomplex_t IT_0310 = IT_0067*IT_0309;
    const ccomplex_t IT_0311 = 1.4142135623731*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*(IT_0308 + 0.5*IT_0311);
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0315 = IT_0002*IT_0314;
    const ccomplex_t IT_0316 = 1.4142135623731*IT_0315;
    const ccomplex_t IT_0317 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0318 = IT_0016*IT_0317;
    const ccomplex_t IT_0319 = 1.4142135623731*IT_0318;
    const ccomplex_t IT_0320 = conj(N_d2)*e_em*m_tau*IT_0031*U_stau_10;
    const ccomplex_t IT_0321 = IT_0067*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*(IT_0316 + IT_0319 + 
      -IT_0322);
    const ccomplex_t IT_0324 = (-0.5)*IT_0323;
    const ccomplex_t IT_0325 = IT_0313*IT_0324;
    const ccomplex_t IT_0326 = IT_0208*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*IT_0326;
    const ccomplex_t IT_0328 = 0.5*IT_0327;
    const ccomplex_t IT_0329 = 2*conj(IT_0304);
    const ccomplex_t IT_0330 = conj(IT_0231)*IT_0236 + IT_0231*IT_0237 + conj
      (IT_0241)*IT_0266 + IT_0241*IT_0267 + IT_0272*IT_0296 + IT_0297*IT_0298 +
       conj(IT_0271)*IT_0299 + IT_0271*IT_0300 + conj(IT_0265)*IT_0305 + 2*
      (IT_0241 + 0.5*IT_0305)*conj(IT_0328) + IT_0265*IT_0329 + 2*IT_0328*(conj
      (IT_0241) + 0.5*IT_0329);
    const ccomplex_t IT_0331 = IT_0227*IT_0324;
    const ccomplex_t IT_0332 = IT_0208*IT_0331;
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*IT_0332;
    const ccomplex_t IT_0334 = IT_0253*IT_0283;
    const ccomplex_t IT_0335 = IT_0242*IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1)*IT_0335;
    const ccomplex_t IT_0337 = (-0.5)*IT_0333 + (-0.5)*IT_0336;
    const ccomplex_t IT_0338 = IT_0265 + IT_0328;
    const ccomplex_t IT_0339 = conj(IT_0265) + conj(IT_0328);
    const ccomplex_t IT_0340 = IT_0231 + IT_0295;
    const ccomplex_t IT_0341 = conj(IT_0231) + conj(IT_0295);
    const ccomplex_t IT_0342 = IT_0261*IT_0291;
    const ccomplex_t IT_0343 = IT_0242*IT_0342;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*IT_0343;
    const ccomplex_t IT_0345 = (-0.5)*IT_0344;
    const ccomplex_t IT_0346 = IT_0216*IT_0313;
    const ccomplex_t IT_0347 = IT_0208*IT_0346;
    const ccomplex_t IT_0348 = (0 + _Complex_I*1)*IT_0347;
    const ccomplex_t IT_0349 = (-0.5)*IT_0348;
    const ccomplex_t IT_0350 = IT_0345 + IT_0349;
    const ccomplex_t IT_0351 = conj(IT_0345) + conj(IT_0349);
    const ccomplex_t IT_0352 = IT_0002*IT_0015;
    const ccomplex_t IT_0353 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0354 = IT_0352*IT_0353;
    const ccomplex_t IT_0355 = IT_0001*IT_0016;
    const ccomplex_t IT_0356 = IT_0353*IT_0355;
    const ccomplex_t IT_0357 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0358 = IT_0352*IT_0357;
    const ccomplex_t IT_0359 = IT_0355*IT_0357;
    const ccomplex_t IT_0360 = (0 + _Complex_I*1)*(IT_0354 + IT_0356 + 
      -IT_0358 + -IT_0359);
    const ccomplex_t IT_0361 = (-0.5)*IT_0360;
    const ccomplex_t IT_0362 = e_em*IT_0352;
    const ccomplex_t IT_0363 = e_em*IT_0355;
    const ccomplex_t IT_0364 = (0 + _Complex_I*1)*(IT_0362 + -IT_0363);
    const ccomplex_t IT_0365 = 0.5*IT_0364;
    const ccomplex_t IT_0366 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0036 + -IT_0037 + -reg_prop, -1);
    const ccomplex_t IT_0367 = IT_0365*IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*IT_0367;
    const ccomplex_t IT_0369 = IT_0361*IT_0368;
    const ccomplex_t IT_0370 = -IT_0369;
    const ccomplex_t IT_0371 = 8*IT_0337;
    const ccomplex_t IT_0372 = 8*conj(IT_0337);
    const ccomplex_t IT_0373 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0374 = IT_0352*IT_0373;
    const ccomplex_t IT_0375 = IT_0355*IT_0373;
    const ccomplex_t IT_0376 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0377 = IT_0352*IT_0376;
    const ccomplex_t IT_0378 = IT_0355*IT_0376;
    const ccomplex_t IT_0379 = (0 + _Complex_I*1)*(IT_0374 + IT_0375 + 
      -IT_0377 + -IT_0378);
    const ccomplex_t IT_0380 = (-0.5)*IT_0379;
    const ccomplex_t IT_0381 = (0 + _Complex_I*1)*e_em*IT_0002*IT_0015;
    const ccomplex_t IT_0382 = IT_0366*IT_0381;
    const ccomplex_t IT_0383 = (0 + _Complex_I*1)*IT_0382;
    const ccomplex_t IT_0384 = IT_0380*IT_0383;
    const ccomplex_t IT_0385 = -IT_0384;
    const ccomplex_t IT_0386 = 8*IT_0345;
    const ccomplex_t IT_0387 = 8*IT_0349;
    const ccomplex_t IT_0388 = 8*conj(IT_0345);
    const ccomplex_t IT_0389 = 8*conj(IT_0349);
    const ccomplex_t IT_0390 = IT_0337*conj(IT_0337) + IT_0338*IT_0339 +
       IT_0340*IT_0341 + IT_0350*IT_0351 + 0.25*conj(IT_0370)*IT_0371 + 4
      *IT_0370*(conj(IT_0370) + 0.0625*IT_0372) + 0.25*conj(IT_0385)*(IT_0386 +
       IT_0387) + 4*IT_0385*(conj(IT_0385) + 0.0625*IT_0388 + 0.0625*IT_0389);
    const ccomplex_t IT_0391 = m_N_2*m_N_3*IT_0062;
    const ccomplex_t IT_0392 = 2*IT_0271;
    const ccomplex_t IT_0393 = 0.25*IT_0392;
    const ccomplex_t IT_0394 = 0.25*IT_0236;
    const ccomplex_t IT_0395 = IT_0340 + IT_0393 + IT_0394;
    const ccomplex_t IT_0396 = 2*conj(IT_0271);
    const ccomplex_t IT_0397 = 0.25*IT_0396;
    const ccomplex_t IT_0398 = 0.25*IT_0237;
    const ccomplex_t IT_0399 = IT_0341 + IT_0397 + IT_0398;
    const ccomplex_t IT_0400 = 0.25*IT_0300;
    const ccomplex_t IT_0401 = 0.25*IT_0298;
    const ccomplex_t IT_0402 = IT_0272 + IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = IT_0241 + IT_0304;
    const ccomplex_t IT_0404 = 0.25*IT_0299;
    const ccomplex_t IT_0405 = 0.25*IT_0296;
    const ccomplex_t IT_0406 = IT_0297 + IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = conj(IT_0241) + conj(IT_0304);
    const ccomplex_t IT_0408 = IT_0361*IT_0383;
    const ccomplex_t IT_0409 = -IT_0408;
    const ccomplex_t IT_0410 = 8*IT_0160;
    const ccomplex_t IT_0411 = IT_0368*IT_0380;
    const ccomplex_t IT_0412 = -IT_0411;
    const ccomplex_t IT_0413 = 16*IT_0412;
    const ccomplex_t IT_0414 = IT_0410 + IT_0413;
    const ccomplex_t IT_0415 = 8*conj(IT_0160);
    const ccomplex_t IT_0416 = 8*IT_0104;
    const ccomplex_t IT_0417 = 8*conj(IT_0104);
    const ccomplex_t IT_0418 = conj(IT_0104)*IT_0160 + IT_0104*conj(IT_0160) +
       conj(IT_0337)*IT_0350 + IT_0337*IT_0351 + 0.25*IT_0372*IT_0385 + 0.25
      *IT_0371*conj(IT_0385) + 0.25*conj(IT_0370)*(16*IT_0385 + IT_0386 +
       IT_0387) + 0.25*IT_0370*(16*conj(IT_0385) + IT_0388 + IT_0389) + IT_0339
      *IT_0395 + IT_0338*IT_0399 + IT_0402*IT_0403 + IT_0406*IT_0407 + 0.25*conj
      (IT_0409)*IT_0414 + 0.25*IT_0409*(16*conj(IT_0412) + IT_0415) + 0.25*conj
      (IT_0412)*IT_0416 + 0.25*IT_0412*IT_0417;
    const ccomplex_t IT_0419 = IT_0104*conj(IT_0104) + IT_0160*conj(IT_0160) +
       IT_0272*IT_0297 + IT_0403*IT_0407 + 0.25*conj(IT_0412)*IT_0414 + 0.25
      *IT_0412*IT_0415 + 0.25*conj(IT_0409)*IT_0416 + 4*IT_0409*(conj(IT_0409) +
       0.0625*IT_0417);
    const ccomplex_t IT_0420 = s_13 + s_14;
    const ccomplex_t IT_0421 = IT_0105*IT_0420;
    const ccomplex_t IT_0422 = IT_0002*IT_0032;
    const ccomplex_t IT_0423 = IT_0166*IT_0422;
    const ccomplex_t IT_0424 = sin(beta);
    const ccomplex_t IT_0425 = IT_0002*IT_0424;
    const ccomplex_t IT_0426 = IT_0168*IT_0425;
    const ccomplex_t IT_0427 = IT_0016*IT_0032;
    const ccomplex_t IT_0428 = IT_0175*IT_0427;
    const ccomplex_t IT_0429 = IT_0016*IT_0424;
    const ccomplex_t IT_0430 = IT_0177*IT_0429;
    const ccomplex_t IT_0431 = IT_0170*IT_0422;
    const ccomplex_t IT_0432 = IT_0172*IT_0425;
    const ccomplex_t IT_0433 = IT_0179*IT_0427;
    const ccomplex_t IT_0434 = IT_0181*IT_0429;
    const ccomplex_t IT_0435 = (-0.5)*IT_0423 + 0.5*IT_0426 + 0.5*IT_0428 + (
      -0.5)*IT_0430 + (-0.5)*IT_0431 + 0.5*IT_0432 + 0.5*IT_0433 + (-0.5)*IT_0434;
    const ccomplex_t IT_0436 = e_em*m_tau*IT_0016*IT_0031*IT_0033*IT_0424;
    const ccomplex_t IT_0437 = (-0.5)*IT_0436;
    const ccomplex_t IT_0438 = cpow((-2)*s_12 + -IT_0036 + -IT_0037 + m_A0*
      (m_A0 + (0 + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0439 = IT_0437*IT_0438;
    const ccomplex_t IT_0440 = (0 + _Complex_I*1)*IT_0439;
    const ccomplex_t IT_0441 = IT_0435*IT_0440;
    const ccomplex_t IT_0442 = IT_0168*IT_0422;
    const ccomplex_t IT_0443 = IT_0166*IT_0425;
    const ccomplex_t IT_0444 = IT_0177*IT_0427;
    const ccomplex_t IT_0445 = IT_0175*IT_0429;
    const ccomplex_t IT_0446 = IT_0172*IT_0422;
    const ccomplex_t IT_0447 = IT_0170*IT_0425;
    const ccomplex_t IT_0448 = IT_0181*IT_0427;
    const ccomplex_t IT_0449 = IT_0179*IT_0429;
    const ccomplex_t IT_0450 = (-0.5)*IT_0442 + (-0.5)*IT_0443 + 0.5*IT_0444 +
       0.5*IT_0445 + (-0.5)*IT_0446 + (-0.5)*IT_0447 + 0.5*IT_0448 + 0.5*IT_0449;
    const ccomplex_t IT_0451 = e_em*m_tau*IT_0016*IT_0031;
    const ccomplex_t IT_0452 = 0.5*IT_0451;
    const ccomplex_t IT_0453 = pow(m_Z, 2);
    const ccomplex_t IT_0454 = cpow((-2)*s_12 + -IT_0036 + -IT_0037 + IT_0453 
      + -reg_prop, -1);
    const ccomplex_t IT_0455 = IT_0452*IT_0454;
    const ccomplex_t IT_0456 = (0 + _Complex_I*1)*IT_0455;
    const ccomplex_t IT_0457 = IT_0450*IT_0456;
    const ccomplex_t IT_0458 = -IT_0441 + -IT_0457;
    const ccomplex_t IT_0459 = (-4)*conj(IT_0458);
    const ccomplex_t IT_0460 = IT_0004*IT_0422;
    const ccomplex_t IT_0461 = IT_0008*IT_0425;
    const ccomplex_t IT_0462 = IT_0018*IT_0427;
    const ccomplex_t IT_0463 = IT_0021*IT_0429;
    const ccomplex_t IT_0464 = IT_0010*IT_0422;
    const ccomplex_t IT_0465 = IT_0012*IT_0425;
    const ccomplex_t IT_0466 = IT_0023*IT_0427;
    const ccomplex_t IT_0467 = IT_0025*IT_0429;
    const ccomplex_t IT_0468 = 0.5*IT_0460 + (-0.5)*IT_0461 + (-0.5)*IT_0462 +
       0.5*IT_0463 + 0.5*IT_0464 + (-0.5)*IT_0465 + (-0.5)*IT_0466 + 0.5*IT_0467;
    const ccomplex_t IT_0469 = IT_0440*IT_0468;
    const ccomplex_t IT_0470 = IT_0008*IT_0422;
    const ccomplex_t IT_0471 = IT_0004*IT_0425;
    const ccomplex_t IT_0472 = IT_0021*IT_0427;
    const ccomplex_t IT_0473 = IT_0018*IT_0429;
    const ccomplex_t IT_0474 = IT_0012*IT_0422;
    const ccomplex_t IT_0475 = IT_0010*IT_0425;
    const ccomplex_t IT_0476 = IT_0025*IT_0427;
    const ccomplex_t IT_0477 = IT_0023*IT_0429;
    const ccomplex_t IT_0478 = 0.5*IT_0470 + 0.5*IT_0471 + (-0.5)*IT_0472 + (
      -0.5)*IT_0473 + 0.5*IT_0474 + 0.5*IT_0475 + (-0.5)*IT_0476 + (-0.5)*IT_0477;
    const ccomplex_t IT_0479 = IT_0456*IT_0478;
    const ccomplex_t IT_0480 = -IT_0469 + -IT_0479;
    const ccomplex_t IT_0481 = (-4)*IT_0480;
    const ccomplex_t IT_0482 = (-4)*conj(IT_0480);
    const ccomplex_t IT_0483 = s_23 + s_24;
    const ccomplex_t IT_0484 = IT_0161*IT_0483;
    const ccomplex_t IT_0485 = -IT_0204;
    const ccomplex_t IT_0486 = IT_0203 + IT_0485;
    const ccomplex_t IT_0487 = s_12*IT_0062;
    const ccomplex_t IT_0488 = -IT_0487;
    const ccomplex_t IT_0489 = (-2)*IT_0202;
    const ccomplex_t IT_0490 = (-2)*conj(IT_0202);
    const ccomplex_t IT_0491 = (-2)*IT_0061;
    const ccomplex_t IT_0492 = (-2)*conj(IT_0061);
    const ccomplex_t IT_0493 = IT_0206 + IT_0486;
    const ccomplex_t IT_0494 = IT_0487 + IT_0493;
    const ccomplex_t IT_0495 = (-4)*IT_0108;
    const ccomplex_t IT_0496 = conj(IT_0385) + conj(IT_0412);
    const ccomplex_t IT_0497 = -IT_0202;
    const ccomplex_t IT_0498 = IT_0385 + IT_0412;
    const ccomplex_t IT_0499 = -conj(IT_0202);
    const ccomplex_t IT_0500 = -IT_0061;
    const ccomplex_t IT_0501 = -conj(IT_0061);
    const ccomplex_t IT_0502 = (-2)*IT_0108;
    const ccomplex_t IT_0503 = 4*IT_0164;
    const ccomplex_t IT_0504 = m_N_2*m_N_3;
    const ccomplex_t IT_0505 = -IT_0062;
    const ccomplex_t IT_0506 = s_34 + IT_0505;
    const ccomplex_t IT_0507 = IT_0504*IT_0506;
    const ccomplex_t IT_0508 = (-2)*IT_0507;
    const ccomplex_t IT_0509 = 2*IT_0164;
    const ccomplex_t IT_0510 = s_12*IT_0506;
    const ccomplex_t IT_0511 = (-8)*IT_0507;
    const ccomplex_t IT_0512 = 0.125*conj(IT_0160);
    const ccomplex_t IT_0513 = 0.125*conj(IT_0104);
    const ccomplex_t IT_0514 = s_24*m_N_2*m_tau;
    const ccomplex_t IT_0515 = 8*IT_0241;
    const ccomplex_t IT_0516 = 0.25*IT_0515;
    const ccomplex_t IT_0517 = 8*IT_0304;
    const ccomplex_t IT_0518 = 0.25*IT_0517;
    const ccomplex_t IT_0519 = IT_0338 + IT_0516 + IT_0518;
    const ccomplex_t IT_0520 = 8*conj(IT_0241);
    const ccomplex_t IT_0521 = 0.25*IT_0520;
    const ccomplex_t IT_0522 = 8*conj(IT_0304);
    const ccomplex_t IT_0523 = 0.25*IT_0522;
    const ccomplex_t IT_0524 = IT_0339 + IT_0521 + IT_0523;
    const ccomplex_t IT_0525 = 8*IT_0271;
    const ccomplex_t IT_0526 = 0.25*IT_0525;
    const ccomplex_t IT_0527 = 8*IT_0235;
    const ccomplex_t IT_0528 = 0.25*IT_0527;
    const ccomplex_t IT_0529 = IT_0340 + IT_0526 + IT_0528;
    const ccomplex_t IT_0530 = 8*conj(IT_0271);
    const ccomplex_t IT_0531 = 0.25*IT_0530;
    const ccomplex_t IT_0532 = 8*conj(IT_0235);
    const ccomplex_t IT_0533 = 0.25*IT_0532;
    const ccomplex_t IT_0534 = IT_0341 + IT_0531 + IT_0533;
    const ccomplex_t IT_0535 = 2*conj(IT_0160);
    const ccomplex_t IT_0536 = 0.25*IT_0338;
    const ccomplex_t IT_0537 = 2*IT_0160;
    const ccomplex_t IT_0538 = 0.25*IT_0339;
    const ccomplex_t IT_0539 = 4*IT_0514*(IT_0104*IT_0402 + conj(IT_0160)
      *IT_0403 + conj(IT_0104)*IT_0406 + IT_0160*IT_0407 + conj(IT_0412)*IT_0519
       + IT_0412*IT_0524 + conj(IT_0409)*IT_0529 + IT_0409*IT_0534 + IT_0535
      *IT_0536 + IT_0537*IT_0538);
    const ccomplex_t IT_0540 = s_13*m_N_3*m_tau;
    const ccomplex_t IT_0541 = 2*conj(IT_0104);
    const ccomplex_t IT_0542 = 2*IT_0104;
    const ccomplex_t IT_0543 = 4*IT_0540*(IT_0160*IT_0402 + conj(IT_0104)
      *IT_0403 + conj(IT_0160)*IT_0406 + IT_0104*IT_0407 + conj(IT_0409)*IT_0519
       + IT_0409*IT_0524 + conj(IT_0412)*IT_0529 + IT_0412*IT_0534 + IT_0536
      *IT_0541 + IT_0538*IT_0542);
    const ccomplex_t IT_0544 = s_14*m_N_3*m_tau;
    const ccomplex_t IT_0545 = 8*IT_0231;
    const ccomplex_t IT_0546 = 0.25*IT_0545;
    const ccomplex_t IT_0547 = 8*IT_0295;
    const ccomplex_t IT_0548 = 0.25*IT_0547;
    const ccomplex_t IT_0549 = IT_0297 + IT_0546 + IT_0548;
    const ccomplex_t IT_0550 = 8*conj(IT_0231);
    const ccomplex_t IT_0551 = 0.25*IT_0550;
    const ccomplex_t IT_0552 = 8*conj(IT_0295);
    const ccomplex_t IT_0553 = 0.25*IT_0552;
    const ccomplex_t IT_0554 = IT_0272 + IT_0551 + IT_0553;
    const ccomplex_t IT_0555 = 8*IT_0328;
    const ccomplex_t IT_0556 = 0.25*IT_0555;
    const ccomplex_t IT_0557 = 8*IT_0265;
    const ccomplex_t IT_0558 = 0.25*IT_0557;
    const ccomplex_t IT_0559 = IT_0403 + IT_0556 + IT_0558;
    const ccomplex_t IT_0560 = 8*conj(IT_0328);
    const ccomplex_t IT_0561 = 0.25*IT_0560;
    const ccomplex_t IT_0562 = 8*conj(IT_0265);
    const ccomplex_t IT_0563 = 0.25*IT_0562;
    const ccomplex_t IT_0564 = IT_0407 + IT_0561 + IT_0563;
    const ccomplex_t IT_0565 = 2*conj(IT_0337);
    const ccomplex_t IT_0566 = 0.25*IT_0403;
    const ccomplex_t IT_0567 = 2*IT_0337;
    const ccomplex_t IT_0568 = 0.25*IT_0407;
    const ccomplex_t IT_0569 = 4*IT_0544*(conj(IT_0337)*IT_0338 + IT_0337
      *IT_0339 + IT_0351*IT_0395 + IT_0350*IT_0399 + conj(IT_0385)*IT_0549 +
       IT_0385*IT_0554 + conj(IT_0370)*IT_0559 + IT_0370*IT_0564 + IT_0565
      *IT_0566 + IT_0567*IT_0568);
    const ccomplex_t IT_0570 = s_23*m_N_2*m_tau;
    const ccomplex_t IT_0571 = 2*conj(IT_0345);
    const ccomplex_t IT_0572 = 2*conj(IT_0349);
    const ccomplex_t IT_0573 = 2*IT_0345;
    const ccomplex_t IT_0574 = 2*IT_0349;
    const ccomplex_t IT_0575 = IT_0339*IT_0350 + IT_0338*IT_0351 + conj
      (IT_0337)*IT_0395 + IT_0337*IT_0399 + conj(IT_0370)*IT_0549 + IT_0370
      *IT_0554 + conj(IT_0385)*IT_0559 + IT_0385*IT_0564 + IT_0566*(IT_0571 +
       IT_0572) + IT_0568*(IT_0573 + IT_0574);
    const ccomplex_t IT_0576 = s_34*m_N_2*m_N_3;
    const ccomplex_t IT_0577 = 8*IT_0412;
    const ccomplex_t IT_0578 = 0.25*IT_0577;
    const ccomplex_t IT_0579 = IT_0160 + IT_0578;
    const ccomplex_t IT_0580 = 8*conj(IT_0412);
    const ccomplex_t IT_0581 = 0.25*IT_0580;
    const ccomplex_t IT_0582 = conj(IT_0160) + IT_0581;
    const ccomplex_t IT_0583 = 0.25*IT_0571;
    const ccomplex_t IT_0584 = 0.5*IT_0104;
    const ccomplex_t IT_0585 = 0.25*IT_0573;
    const ccomplex_t IT_0586 = 0.5*conj(IT_0104);
    const ccomplex_t IT_0587 = (0.5*conj(IT_0345) + 0.5*conj(IT_0349) + conj
      (IT_0385))*IT_0409 + (0.5*IT_0345 + 0.5*IT_0349 + IT_0385)*conj(IT_0409) +
       0.5*conj(IT_0337)*IT_0412 + 0.5*IT_0337*conj(IT_0412) + 0.125*IT_0349
      *IT_0541 + 0.125*conj(IT_0349)*IT_0542 + 0.125*IT_0160*IT_0565 + IT_0512
      *IT_0567 + 0.5*conj(IT_0370)*IT_0579 + 0.5*IT_0370*IT_0582 + (conj(IT_0385
      ) + IT_0583)*IT_0584 + (IT_0385 + IT_0585)*IT_0586;
    const ccomplex_t IT_0588 = conj(IT_0337)*IT_0409 + IT_0337*conj(IT_0409) +
       IT_0350*(conj(IT_0412) + 0.25*IT_0535) + IT_0351*(IT_0412 + 0.25*IT_0537)
       + 0.25*IT_0104*IT_0565 + 0.25*conj(IT_0104)*IT_0567 + conj(IT_0385)
      *IT_0579 + IT_0385*IT_0582 + 2*conj(IT_0370)*(IT_0409 + IT_0584) + 2
      *IT_0370*(conj(IT_0409) + IT_0586);
    const ccomplex_t IT_0589 = s_34 + IT_0062;
    const ccomplex_t IT_0590 = IT_0504*IT_0589;
    const ccomplex_t IT_0591 = (-2)*IT_0304;
    const ccomplex_t IT_0592 = (-2)*IT_0241;
    const ccomplex_t IT_0593 = 2*IT_0484;
    const ccomplex_t IT_0594 = 2*IT_0421;
    const ccomplex_t IT_0595 = IT_0205 + IT_0486 + IT_0487;
    const ccomplex_t IT_0596 = 0.5*IT_0595;
    const ccomplex_t IT_0597 = IT_0488 + IT_0493;
    const ccomplex_t IT_0598 = 0.5*IT_0597;
    const ccomplex_t IT_0599 = IT_0590*(0.5*IT_0266 + IT_0328 + 0.5*IT_0591 +
       0.5*IT_0592) + ((-0.5)*IT_0104 + IT_0370 + 0.25*IT_0567)*IT_0593 + ((-0.5
      )*IT_0160 + IT_0385 + 0.25*IT_0574 + IT_0585)*IT_0594 + (IT_0296 + IT_0299
      )*IT_0596 + (IT_0236 + IT_0392)*IT_0598;
    const ccomplex_t IT_0600 = (-2)*conj(IT_0304);
    const ccomplex_t IT_0601 = (-2)*conj(IT_0241);
    const ccomplex_t IT_0602 = s_12*IT_0589;
    const ccomplex_t IT_0603 = 4*conj(IT_0458);
    const ccomplex_t IT_0604 = ((-0.5)*conj(IT_0104) + conj(IT_0370) + -conj
      (IT_0409) + 0.25*IT_0565)*IT_0593 + ((-0.5)*conj(IT_0160) + conj(IT_0385) 
      + -conj(IT_0412) + 0.25*IT_0572 + IT_0583)*IT_0594 + (IT_0298 + IT_0300)
      *IT_0596 + (IT_0237 + IT_0396)*IT_0598 + IT_0590*(0.5*IT_0267 + conj
      (IT_0328) + IT_0482 + 0.5*IT_0600 + 0.5*IT_0601) + IT_0602*IT_0603;
    const ccomplex_t IT_0605 = (-2)*IT_0421;
    const ccomplex_t IT_0606 = (-0.5)*IT_0597;
    const ccomplex_t IT_0607 = (-2)*IT_0484*((-0.5)*IT_0350 + IT_0412 + 0.25
      *IT_0537) + ((-0.5)*IT_0236 + IT_0340 + (-0.5)*IT_0392)*IT_0590 + IT_0338
      *IT_0595 + ((-0.5)*IT_0337 + IT_0409 + 0.25*IT_0542)*IT_0605 + (IT_0591 +
       IT_0592)*IT_0606;
    const ccomplex_t IT_0608 = (-2)*IT_0484*((-0.5)*IT_0351 + conj(IT_0412) +
       0.25*IT_0535) + IT_0339*IT_0595 + IT_0482*IT_0602 + IT_0590*((-0.5)
      *IT_0237 + IT_0341 + (-0.5)*IT_0396 + IT_0603) + ((-0.5)*conj(IT_0337) +
       conj(IT_0409) + 0.25*IT_0541)*IT_0605 + (IT_0600 + IT_0601)*IT_0606;
    const ccomplex_t IT_0609 = conj(IT_0061)*(IT_0104*IT_0109 + IT_0160
      *IT_0165) + (IT_0109*IT_0160 + IT_0104*IT_0165)*conj(IT_0202) + IT_0207
      *IT_0330 + 4*IT_0203*IT_0390 + 4*IT_0391*IT_0418 + 4*IT_0204*IT_0419 +
       IT_0421*(IT_0412*IT_0459 + conj(IT_0370)*IT_0481 + IT_0370*IT_0482) + 
      (IT_0409*IT_0459 + conj(IT_0385)*IT_0481 + IT_0385*IT_0482)*IT_0484 + 
      (IT_0205 + IT_0486 + IT_0488)*(IT_0341*IT_0489 + IT_0340*IT_0490 + IT_0339
      *IT_0491 + IT_0338*IT_0492) + (IT_0272*IT_0489 + IT_0297*IT_0490 + IT_0407
      *IT_0491 + IT_0403*IT_0492)*IT_0494 + IT_0495*(IT_0496*IT_0497 + IT_0498
      *IT_0499 + (conj(IT_0370) + conj(IT_0409))*IT_0500 + (IT_0370 + IT_0409)
      *IT_0501) + (IT_0351*IT_0497 + IT_0350*IT_0499 + conj(IT_0337)*IT_0500 +
       IT_0337*IT_0501)*IT_0502 + ((conj(IT_0370) + conj(IT_0409))*IT_0497 + 
      (IT_0370 + IT_0409)*IT_0499 + IT_0496*IT_0500 + IT_0498*IT_0501)*IT_0503 +
       (IT_0061*IT_0272 + conj(IT_0061)*IT_0297 + conj(IT_0202)*IT_0403 +
       IT_0202*IT_0407 + IT_0339*IT_0497 + IT_0338*IT_0499 + IT_0341*IT_0500 +
       IT_0340*IT_0501)*IT_0508 + (conj(IT_0337)*IT_0497 + IT_0337*IT_0499 +
       IT_0351*IT_0500 + IT_0350*IT_0501)*IT_0509 + 8*IT_0061*(conj(IT_0061)
      *IT_0510 + 0.125*conj(IT_0202)*IT_0511 + IT_0165*IT_0512 + IT_0109*IT_0513
      ) + 8*IT_0202*(conj(IT_0202)*IT_0510 + 0.125*conj(IT_0061)*IT_0511 +
       IT_0109*IT_0512 + IT_0165*IT_0513) + IT_0539 + IT_0543 + IT_0569 + 4
      *IT_0570*IT_0575 + 8*IT_0576*IT_0587 + 4*IT_0487*IT_0588 + 2*conj(IT_0458)
      *IT_0599 + 2*IT_0458*IT_0604 + (-2)*conj(IT_0480)*IT_0607 + (-2)*IT_0480
      *IT_0608;
    return create_ccomplex_return(IT_0609);
}

