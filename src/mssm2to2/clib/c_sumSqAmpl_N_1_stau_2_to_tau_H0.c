#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_2_to_tau_H0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_2_to_tau_H0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
    const creal_t alpha = param->alpha;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_tau = param->m_tau;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = pow(m_tau, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = cos(beta);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = N_d1*e_em*m_tau*IT_0013*U_stau_01;
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0007 + 0.5*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = cos(alpha);
    const ccomplex_t IT_0020 = mu_h*m_tau*IT_0013;
    const ccomplex_t IT_0021 = sin(alpha);
    const ccomplex_t IT_0022 = sin(beta);
    const ccomplex_t IT_0023 = -IT_0021*IT_0022;
    const ccomplex_t IT_0024 = IT_0008*IT_0019;
    const ccomplex_t IT_0025 = IT_0023 + IT_0024;
    const ccomplex_t IT_0026 = cpow(IT_0003, -2);
    const ccomplex_t IT_0027 = (-2)*U_stau_01*((IT_0001*IT_0009*IT_0011
      *IT_0013*IT_0019 + (-0.5)*m_W*IT_0025*(IT_0011 + -IT_0010*IT_0026))*conj
      (U_stau_01) + -1./2*IT_0009*IT_0011*IT_0020*IT_0021*conj(U_stau_11)) +
       U_stau_11*(IT_0009*IT_0011*IT_0020*IT_0021*conj(U_stau_01) + (-2)*
      (IT_0001*IT_0009*IT_0011*IT_0013*IT_0019 + -m_W*IT_0010*IT_0025*IT_0026)
      *conj(U_stau_11));
    const ccomplex_t IT_0028 = (0 + _Complex_I*-1)*(e_em*IT_0027 + (
      -1.4142135623731)*A_tau*IT_0019*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0018*IT_0029;
    const ccomplex_t IT_0031 = IT_0002*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0034 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0035 = IT_0004*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_d1*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0038 = IT_0012*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = (-2)*conj(U_stau_00)*((IT_0001*IT_0009*IT_0011
      *IT_0013*IT_0019 + (-0.5)*m_W*IT_0025*(IT_0011 + -IT_0010*IT_0026))
      *U_stau_01 + -1./2*IT_0009*IT_0011*IT_0020*IT_0021*U_stau_11) + conj
      (U_stau_10)*(IT_0009*IT_0011*IT_0020*IT_0021*U_stau_01 + (-2)*(IT_0001
      *IT_0009*IT_0011*IT_0013*IT_0019 + -m_W*IT_0010*IT_0025*IT_0026)*U_stau_11);
    const ccomplex_t IT_0043 = (0 + _Complex_I*-1)*(e_em*IT_0042 + (
      -1.4142135623731)*A_tau*IT_0019*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = IT_0033*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0049 = IT_0004*IT_0019;
    const ccomplex_t IT_0050 = IT_0048*IT_0049;
    const ccomplex_t IT_0051 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0052 = IT_0004*IT_0021;
    const ccomplex_t IT_0053 = IT_0051*IT_0052;
    const ccomplex_t IT_0054 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0055 = IT_0011*IT_0019;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0058 = IT_0011*IT_0021;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = (0 + _Complex_I*1)*(IT_0050 + -IT_0053 + 
      -IT_0056 + IT_0059);
    const ccomplex_t IT_0061 = m_N_1*IT_0060;
    const ccomplex_t IT_0062 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0063 = IT_0049*IT_0062;
    const ccomplex_t IT_0064 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0065 = IT_0052*IT_0064;
    const ccomplex_t IT_0066 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0067 = IT_0055*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0069 = IT_0058*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0063 + -IT_0065 + 
      -IT_0067 + IT_0069);
    const ccomplex_t IT_0071 = m_N_1*IT_0070;
    const ccomplex_t IT_0072 = IT_0061 + IT_0071;
    const ccomplex_t IT_0073 = pow(m_stau_2, 2);
    const ccomplex_t IT_0074 = cpow(s_23 + (-0.5)*IT_0001 + (-0.5)*IT_0073 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0018*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0072*IT_0076;
    const ccomplex_t IT_0078 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0079 = IT_0004*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = N_d2*e_em*m_tau*IT_0013*U_stau_01;
    const ccomplex_t IT_0082 = IT_0012*IT_0081;
    const ccomplex_t IT_0083 = 1.4142135623731*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*(IT_0080 + 0.5*IT_0083);
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = cpow(s_23 + (-0.5)*IT_0001 + (-0.5)*IT_0073 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0090 = IT_0049*IT_0089;
    const ccomplex_t IT_0091 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0092 = IT_0058*IT_0091;
    const ccomplex_t IT_0093 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0094 = IT_0049*IT_0093;
    const ccomplex_t IT_0095 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0096 = IT_0058*IT_0095;
    const ccomplex_t IT_0097 = IT_0090 + IT_0092 + IT_0094 + IT_0096;
    const ccomplex_t IT_0098 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0099 = IT_0052*IT_0098;
    const ccomplex_t IT_0100 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0101 = IT_0055*IT_0100;
    const ccomplex_t IT_0102 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0103 = IT_0052*IT_0102;
    const ccomplex_t IT_0104 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0105 = IT_0055*IT_0104;
    const ccomplex_t IT_0106 = -IT_0099 + -IT_0101 + -IT_0103 + -IT_0105;
    const ccomplex_t IT_0107 = IT_0097 + IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = m_N_2*IT_0109;
    const ccomplex_t IT_0111 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0112 = IT_0049*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0114 = IT_0058*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0116 = IT_0049*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0118 = IT_0058*IT_0117;
    const ccomplex_t IT_0119 = IT_0112 + IT_0114 + IT_0116 + IT_0118;
    const ccomplex_t IT_0120 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0121 = IT_0052*IT_0120;
    const ccomplex_t IT_0122 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0123 = IT_0055*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0125 = IT_0052*IT_0124;
    const ccomplex_t IT_0126 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0127 = IT_0055*IT_0126;
    const ccomplex_t IT_0128 = -IT_0121 + -IT_0123 + -IT_0125 + -IT_0127;
    const ccomplex_t IT_0129 = IT_0119 + IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = 0.5*IT_0130;
    const ccomplex_t IT_0132 = m_N_1*IT_0131;
    const ccomplex_t IT_0133 = IT_0110 + IT_0132;
    const ccomplex_t IT_0134 = IT_0088*IT_0133;
    const ccomplex_t IT_0135 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0136 = IT_0049*IT_0135;
    const ccomplex_t IT_0137 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0138 = IT_0058*IT_0137;
    const ccomplex_t IT_0139 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0140 = IT_0049*IT_0139;
    const ccomplex_t IT_0141 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0142 = IT_0058*IT_0141;
    const ccomplex_t IT_0143 = IT_0136 + IT_0138 + IT_0140 + IT_0142;
    const ccomplex_t IT_0144 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0145 = IT_0052*IT_0144;
    const ccomplex_t IT_0146 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0147 = IT_0055*IT_0146;
    const ccomplex_t IT_0148 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0149 = IT_0052*IT_0148;
    const ccomplex_t IT_0150 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0151 = IT_0055*IT_0150;
    const ccomplex_t IT_0152 = -IT_0145 + -IT_0147 + -IT_0149 + -IT_0151;
    const ccomplex_t IT_0153 = IT_0143 + IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = 0.5*IT_0154;
    const ccomplex_t IT_0156 = m_N_3*IT_0155;
    const ccomplex_t IT_0157 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0158 = IT_0049*IT_0157;
    const ccomplex_t IT_0159 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0160 = IT_0058*IT_0159;
    const ccomplex_t IT_0161 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0162 = IT_0049*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0164 = IT_0058*IT_0163;
    const ccomplex_t IT_0165 = IT_0158 + IT_0160 + IT_0162 + IT_0164;
    const ccomplex_t IT_0166 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0167 = IT_0052*IT_0166;
    const ccomplex_t IT_0168 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0169 = IT_0055*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0171 = IT_0052*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0173 = IT_0055*IT_0172;
    const ccomplex_t IT_0174 = -IT_0167 + -IT_0169 + -IT_0171 + -IT_0173;
    const ccomplex_t IT_0175 = IT_0165 + IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = 0.5*IT_0176;
    const ccomplex_t IT_0178 = m_N_1*IT_0177;
    const ccomplex_t IT_0179 = IT_0156 + IT_0178;
    const ccomplex_t IT_0180 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0181 = IT_0004*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = N_d3*e_em*m_tau*IT_0013*U_stau_01;
    const ccomplex_t IT_0184 = IT_0012*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + 0.5*IT_0185);
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = cpow(s_23 + (-0.5)*IT_0001 + (-0.5)*IT_0073 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0189 = IT_0187*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = IT_0179*IT_0190;
    const ccomplex_t IT_0192 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0193 = IT_0004*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = N_d4*e_em*m_tau*IT_0013*U_stau_01;
    const ccomplex_t IT_0196 = IT_0012*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0194 + 0.5*IT_0197);
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = cpow(s_23 + (-0.5)*IT_0001 + (-0.5)*IT_0073 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0201 = IT_0199*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0204 = IT_0049*IT_0203;
    const ccomplex_t IT_0205 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0206 = IT_0058*IT_0205;
    const ccomplex_t IT_0207 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0208 = IT_0049*IT_0207;
    const ccomplex_t IT_0209 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0210 = IT_0058*IT_0209;
    const ccomplex_t IT_0211 = IT_0204 + IT_0206 + IT_0208 + IT_0210;
    const ccomplex_t IT_0212 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0213 = IT_0052*IT_0212;
    const ccomplex_t IT_0214 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0215 = IT_0055*IT_0214;
    const ccomplex_t IT_0216 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0217 = IT_0052*IT_0216;
    const ccomplex_t IT_0218 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0219 = IT_0055*IT_0218;
    const ccomplex_t IT_0220 = -IT_0213 + -IT_0215 + -IT_0217 + -IT_0219;
    const ccomplex_t IT_0221 = IT_0211 + IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0222;
    const ccomplex_t IT_0224 = m_N_4*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0226 = IT_0049*IT_0225;
    const ccomplex_t IT_0227 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0228 = IT_0058*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0230 = IT_0049*IT_0229;
    const ccomplex_t IT_0231 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0232 = IT_0058*IT_0231;
    const ccomplex_t IT_0233 = IT_0226 + IT_0228 + IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0235 = IT_0052*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0237 = IT_0055*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0239 = IT_0052*IT_0238;
    const ccomplex_t IT_0240 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0241 = IT_0055*IT_0240;
    const ccomplex_t IT_0242 = -IT_0235 + -IT_0237 + -IT_0239 + -IT_0241;
    const ccomplex_t IT_0243 = IT_0233 + IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = m_N_1*IT_0245;
    const ccomplex_t IT_0247 = IT_0224 + IT_0246;
    const ccomplex_t IT_0248 = IT_0202*IT_0247;
    const ccomplex_t IT_0249 = cpow(s_12 + 0.5*IT_0000 + (-0.5)*IT_0001 + 0.5
      *IT_0073 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*e_em*m_tau*IT_0009*IT_0011
      *IT_0013*IT_0019;
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = IT_0018*IT_0251;
    const ccomplex_t IT_0253 = IT_0249*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = m_tau*IT_0254;
    const ccomplex_t IT_0256 = (-0.5)*IT_0032 + (-0.5)*IT_0047 + (-0.5)
      *IT_0077 + (-0.5)*IT_0134 + (-0.5)*IT_0191 + (-0.5)*IT_0248 + IT_0255;
    const ccomplex_t IT_0257 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0258 = IT_0004*IT_0257;
    const ccomplex_t IT_0259 = 1.4142135623731*IT_0258;
    const ccomplex_t IT_0260 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0261 = IT_0011*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = conj(N_d1)*e_em*m_tau*IT_0013*U_stau_11;
    const ccomplex_t IT_0264 = IT_0012*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = (0 + _Complex_I*1)*(IT_0259 + IT_0262 + 
      -IT_0265);
    const ccomplex_t IT_0267 = (-0.5)*IT_0266;
    const ccomplex_t IT_0268 = IT_0074*IT_0267;
    const ccomplex_t IT_0269 = (0 + _Complex_I*1)*IT_0268;
    const ccomplex_t IT_0270 = IT_0072*IT_0269;
    const ccomplex_t IT_0271 = IT_0251*IT_0267;
    const ccomplex_t IT_0272 = IT_0249*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = m_tau*IT_0273;
    const ccomplex_t IT_0275 = (-2)*U_stau_01*((IT_0001*IT_0009*IT_0011
      *IT_0013*IT_0019 + (-0.5)*m_W*IT_0025*(IT_0011 + -IT_0010*IT_0026))*conj
      (U_stau_01) + -1./2*IT_0009*IT_0011*IT_0020*IT_0021*conj(U_stau_11)) +
       U_stau_11*(IT_0009*IT_0011*IT_0020*IT_0021*conj(U_stau_01) + (-2)*
      (IT_0001*IT_0009*IT_0011*IT_0013*IT_0019 + -m_W*IT_0010*IT_0025*IT_0026)
      *conj(U_stau_11));
    const ccomplex_t IT_0276 = (0 + _Complex_I*-1)*(e_em*IT_0275 + (
      -1.4142135623731)*A_tau*IT_0019*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = IT_0267*IT_0277;
    const ccomplex_t IT_0279 = IT_0002*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*IT_0279;
    const ccomplex_t IT_0281 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0282 = IT_0004*IT_0281;
    const ccomplex_t IT_0283 = 1.4142135623731*IT_0282;
    const ccomplex_t IT_0284 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0285 = IT_0011*IT_0284;
    const ccomplex_t IT_0286 = 1.4142135623731*IT_0285;
    const ccomplex_t IT_0287 = conj(N_d1)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0288 = IT_0012*IT_0287;
    const ccomplex_t IT_0289 = 1.4142135623731*IT_0288;
    const ccomplex_t IT_0290 = (0 + _Complex_I*1)*(IT_0283 + IT_0286 + 
      -IT_0289);
    const ccomplex_t IT_0291 = (-0.5)*IT_0290;
    const ccomplex_t IT_0292 = (-2)*conj(U_stau_00)*((IT_0001*IT_0009*IT_0011
      *IT_0013*IT_0019 + (-0.5)*m_W*IT_0025*(IT_0011 + -IT_0010*IT_0026))
      *U_stau_01 + -1./2*IT_0009*IT_0011*IT_0020*IT_0021*U_stau_11) + conj
      (U_stau_10)*(IT_0009*IT_0011*IT_0020*IT_0021*U_stau_01 + (-2)*(IT_0001
      *IT_0009*IT_0011*IT_0013*IT_0019 + -m_W*IT_0010*IT_0025*IT_0026)*U_stau_11);
    const ccomplex_t IT_0293 = (0 + _Complex_I*-1)*(e_em*IT_0292 + (
      -1.4142135623731)*A_tau*IT_0019*(U_stau_01*conj(U_stau_10) + conj
      (U_stau_00)*U_stau_11));
    const ccomplex_t IT_0294 = (-0.5)*IT_0293;
    const ccomplex_t IT_0295 = IT_0291*IT_0294;
    const ccomplex_t IT_0296 = IT_0033*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*IT_0296;
    const ccomplex_t IT_0298 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0299 = IT_0004*IT_0298;
    const ccomplex_t IT_0300 = 1.4142135623731*IT_0299;
    const ccomplex_t IT_0301 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0302 = IT_0011*IT_0301;
    const ccomplex_t IT_0303 = 1.4142135623731*IT_0302;
    const ccomplex_t IT_0304 = conj(N_d2)*e_em*m_tau*IT_0013*U_stau_11;
    const ccomplex_t IT_0305 = IT_0012*IT_0304;
    const ccomplex_t IT_0306 = 1.4142135623731*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*(IT_0300 + IT_0303 + 
      -IT_0306);
    const ccomplex_t IT_0308 = (-0.5)*IT_0307;
    const ccomplex_t IT_0309 = IT_0086*IT_0308;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*IT_0309;
    const ccomplex_t IT_0311 = m_N_2*IT_0131;
    const ccomplex_t IT_0312 = m_N_1*IT_0109;
    const ccomplex_t IT_0313 = IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = IT_0310*IT_0313;
    const ccomplex_t IT_0315 = m_N_1*IT_0155;
    const ccomplex_t IT_0316 = m_N_3*IT_0177;
    const ccomplex_t IT_0317 = IT_0315 + IT_0316;
    const ccomplex_t IT_0318 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0319 = IT_0004*IT_0318;
    const ccomplex_t IT_0320 = 1.4142135623731*IT_0319;
    const ccomplex_t IT_0321 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0322 = IT_0011*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = conj(N_d3)*e_em*m_tau*IT_0013*U_stau_11;
    const ccomplex_t IT_0325 = IT_0012*IT_0324;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*(IT_0320 + IT_0323 + 
      -IT_0326);
    const ccomplex_t IT_0328 = (-0.5)*IT_0327;
    const ccomplex_t IT_0329 = IT_0188*IT_0328;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*IT_0329;
    const ccomplex_t IT_0331 = IT_0317*IT_0330;
    const ccomplex_t IT_0332 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0333 = IT_0004*IT_0332;
    const ccomplex_t IT_0334 = 1.4142135623731*IT_0333;
    const ccomplex_t IT_0335 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0336 = IT_0011*IT_0335;
    const ccomplex_t IT_0337 = 1.4142135623731*IT_0336;
    const ccomplex_t IT_0338 = conj(N_d4)*e_em*m_tau*IT_0013*U_stau_11;
    const ccomplex_t IT_0339 = IT_0012*IT_0338;
    const ccomplex_t IT_0340 = 1.4142135623731*IT_0339;
    const ccomplex_t IT_0341 = (0 + _Complex_I*1)*(IT_0334 + IT_0337 + 
      -IT_0340);
    const ccomplex_t IT_0342 = (-0.5)*IT_0341;
    const ccomplex_t IT_0343 = IT_0200*IT_0342;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*IT_0343;
    const ccomplex_t IT_0345 = m_N_1*IT_0223;
    const ccomplex_t IT_0346 = m_N_4*IT_0245;
    const ccomplex_t IT_0347 = IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = IT_0344*IT_0347;
    const ccomplex_t IT_0349 = 0.5*IT_0270 + -IT_0274 + 0.5*IT_0280 + 0.5
      *IT_0297 + 0.5*IT_0314 + 0.5*IT_0331 + 0.5*IT_0348;
    const ccomplex_t IT_0350 = 2*s_13;
    const ccomplex_t IT_0351 = m_N_1*m_tau;
    const ccomplex_t IT_0352 = 2*IT_0256;
    const ccomplex_t IT_0353 = IT_0060*IT_0269;
    const ccomplex_t IT_0354 = IT_0109*IT_0310;
    const ccomplex_t IT_0355 = IT_0155*IT_0330;
    const ccomplex_t IT_0356 = IT_0223*IT_0344;
    const ccomplex_t IT_0357 = 0.5*IT_0254 + (-0.5)*IT_0353 + (-0.5)*IT_0354 +
       (-0.5)*IT_0355 + (-0.5)*IT_0356;
    const ccomplex_t IT_0358 = s_14*m_tau;
    const ccomplex_t IT_0359 = IT_0070*IT_0076;
    const ccomplex_t IT_0360 = IT_0088*IT_0131;
    const ccomplex_t IT_0361 = IT_0177*IT_0190;
    const ccomplex_t IT_0362 = IT_0202*IT_0245;
    const ccomplex_t IT_0363 = (-0.5)*IT_0273 + 0.5*IT_0359 + 0.5*IT_0360 +
       0.5*IT_0361 + 0.5*IT_0362;
    const ccomplex_t IT_0364 = s_34*m_N_1;
    const ccomplex_t IT_0365 = 2*conj(IT_0256);
    const ccomplex_t IT_0366 = 2*conj(IT_0349);
    const ccomplex_t IT_0367 = pow(m_H0, 2);
    const ccomplex_t IT_0368 = m_N_1*m_tau*IT_0367;
    const ccomplex_t IT_0369 = s_13*IT_0367;
    const ccomplex_t IT_0370 = s_14*s_34;
    const ccomplex_t IT_0371 = (-2)*IT_0370;
    const ccomplex_t IT_0372 = IT_0369 + IT_0371;
    const ccomplex_t IT_0373 = 2*IT_0349;
    const ccomplex_t IT_0374 = 0.5*IT_0364;
    const ccomplex_t IT_0375 = (IT_0256*conj(IT_0256) + IT_0349*conj(IT_0349))
      *IT_0350 + conj(IT_0349)*IT_0351*IT_0352 + (IT_0349*IT_0351 + IT_0357
      *IT_0358 + IT_0363*IT_0364)*IT_0365 + (IT_0358*IT_0363 + IT_0357*IT_0364)
      *IT_0366 + 2*conj(IT_0363)*(IT_0357*IT_0368 + -IT_0363*IT_0372 + 0.5
      *IT_0358*IT_0373 + IT_0352*IT_0374) + 2*conj(IT_0357)*(0.5*IT_0352*IT_0358
       + IT_0363*IT_0368 + -IT_0357*IT_0372 + IT_0373*IT_0374);
    return create_ccomplex_return(IT_0375);
}

