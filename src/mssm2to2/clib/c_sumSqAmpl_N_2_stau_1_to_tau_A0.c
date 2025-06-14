#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_stau_1_to_tau_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_stau_1_to_tau_A0(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t mu_h = param->mu_h;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t A_tau = param->A_tau;
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = pow(m_W, -1);
    const ccomplex_t IT_0003 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = e_em*mu_h*m_tau*IT_0002*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (-0.5)*IT_0004 + 0.5*IT_0006 + (-0.5)*IT_0010 +
       0.5*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = N_d2*e_em*m_tau*IT_0002*U_stau_00;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0019 + 0.5*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = pow(m_N_2, 2);
    const ccomplex_t IT_0030 = pow(m_tau, 2);
    const ccomplex_t IT_0031 = cpow(s_13 + (-0.5)*IT_0029 + (-0.5)*IT_0030 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0035 = IT_0016*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_d1*e_em*m_tau*IT_0002*U_stau_00;
    const ccomplex_t IT_0038 = IT_0022*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = pow(m_stau_1, 2);
    const ccomplex_t IT_0043 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0016*IT_0020;
    const ccomplex_t IT_0047 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = IT_0007*IT_0016;
    const ccomplex_t IT_0050 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0051 = IT_0049*IT_0050;
    const ccomplex_t IT_0052 = IT_0001*IT_0020;
    const ccomplex_t IT_0053 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0054 = IT_0052*IT_0053;
    const ccomplex_t IT_0055 = IT_0001*IT_0007;
    const ccomplex_t IT_0056 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0059 = IT_0046*IT_0058;
    const ccomplex_t IT_0060 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0061 = IT_0049*IT_0060;
    const ccomplex_t IT_0062 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0063 = IT_0052*IT_0062;
    const ccomplex_t IT_0064 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0065 = IT_0055*IT_0064;
    const ccomplex_t IT_0066 = (-0.5)*IT_0048 + 0.5*IT_0051 + 0.5*IT_0054 + (
      -0.5)*IT_0057 + (-0.5)*IT_0059 + 0.5*IT_0061 + 0.5*IT_0063 + (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = m_N_1*IT_0066;
    const ccomplex_t IT_0068 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0069 = IT_0046*IT_0068;
    const ccomplex_t IT_0070 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0071 = IT_0049*IT_0070;
    const ccomplex_t IT_0072 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0073 = IT_0052*IT_0072;
    const ccomplex_t IT_0074 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0075 = IT_0055*IT_0074;
    const ccomplex_t IT_0076 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0077 = IT_0046*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0079 = IT_0049*IT_0078;
    const ccomplex_t IT_0080 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0081 = IT_0052*IT_0080;
    const ccomplex_t IT_0082 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0083 = IT_0055*IT_0082;
    const ccomplex_t IT_0084 = 0.5*IT_0069 + (-0.5)*IT_0071 + (-0.5)*IT_0073 +
       0.5*IT_0075 + 0.5*IT_0077 + (-0.5)*IT_0079 + (-0.5)*IT_0081 + 0.5*IT_0083;
    const ccomplex_t IT_0085 = m_N_2*IT_0084;
    const ccomplex_t IT_0086 = IT_0067 + IT_0085;
    const ccomplex_t IT_0087 = IT_0045*IT_0086;
    const ccomplex_t IT_0088 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0089 = IT_0016*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = N_d2*e_em*m_tau*IT_0002*U_stau_01;
    const ccomplex_t IT_0092 = IT_0022*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + 0.5*IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_01)
      *U_stau_10;
    const ccomplex_t IT_0097 = IT_0001*IT_0096;
    const ccomplex_t IT_0098 = e_em*mu_h*m_tau*IT_0002*U_stau_00*conj
      (U_stau_11);
    const ccomplex_t IT_0099 = IT_0001*IT_0098;
    const ccomplex_t IT_0100 = A_tau*conj(U_stau_01)*U_stau_10;
    const ccomplex_t IT_0101 = IT_0007*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = A_tau*U_stau_00*conj(U_stau_11);
    const ccomplex_t IT_0104 = IT_0007*IT_0103;
    const ccomplex_t IT_0105 = 1.4142135623731*IT_0104;
    const ccomplex_t IT_0106 = (-0.5)*IT_0097 + 0.5*IT_0099 + (-0.5)*IT_0102 +
       0.5*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = cpow(s_13 + (-0.5)*IT_0029 + (-0.5)*IT_0030 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0109 = IT_0107*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*IT_0109;
    const ccomplex_t IT_0111 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0112 = IT_0046*IT_0111;
    const ccomplex_t IT_0113 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0114 = IT_0049*IT_0113;
    const ccomplex_t IT_0115 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0116 = IT_0052*IT_0115;
    const ccomplex_t IT_0117 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0118 = IT_0055*IT_0117;
    const ccomplex_t IT_0119 = -IT_0112 + IT_0114 + IT_0116 + -IT_0118;
    const ccomplex_t IT_0120 = m_N_2*IT_0119;
    const ccomplex_t IT_0121 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0122 = IT_0046*IT_0121;
    const ccomplex_t IT_0123 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0124 = IT_0049*IT_0123;
    const ccomplex_t IT_0125 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0126 = IT_0052*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0128 = IT_0055*IT_0127;
    const ccomplex_t IT_0129 = IT_0122 + -IT_0124 + -IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = m_N_2*IT_0129;
    const ccomplex_t IT_0131 = IT_0120 + IT_0130;
    const ccomplex_t IT_0132 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0027*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0131*IT_0134;
    const ccomplex_t IT_0136 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0137 = IT_0016*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = N_d3*e_em*m_tau*IT_0002*U_stau_00;
    const ccomplex_t IT_0140 = IT_0022*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 0.5*IT_0141);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0145 = IT_0143*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0148 = IT_0046*IT_0147;
    const ccomplex_t IT_0149 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0150 = IT_0049*IT_0149;
    const ccomplex_t IT_0151 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0152 = IT_0052*IT_0151;
    const ccomplex_t IT_0153 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0154 = IT_0055*IT_0153;
    const ccomplex_t IT_0155 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0156 = IT_0046*IT_0155;
    const ccomplex_t IT_0157 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0158 = IT_0049*IT_0157;
    const ccomplex_t IT_0159 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0160 = IT_0052*IT_0159;
    const ccomplex_t IT_0161 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0162 = IT_0055*IT_0161;
    const ccomplex_t IT_0163 = (-0.5)*IT_0148 + 0.5*IT_0150 + 0.5*IT_0152 + (
      -0.5)*IT_0154 + (-0.5)*IT_0156 + 0.5*IT_0158 + 0.5*IT_0160 + (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = m_N_3*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0166 = IT_0046*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0168 = IT_0049*IT_0167;
    const ccomplex_t IT_0169 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0170 = IT_0052*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0172 = IT_0055*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0174 = IT_0046*IT_0173;
    const ccomplex_t IT_0175 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0176 = IT_0049*IT_0175;
    const ccomplex_t IT_0177 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0178 = IT_0052*IT_0177;
    const ccomplex_t IT_0179 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0180 = IT_0055*IT_0179;
    const ccomplex_t IT_0181 = 0.5*IT_0166 + (-0.5)*IT_0168 + (-0.5)*IT_0170 +
       0.5*IT_0172 + 0.5*IT_0174 + (-0.5)*IT_0176 + (-0.5)*IT_0178 + 0.5*IT_0180;
    const ccomplex_t IT_0182 = m_N_2*IT_0181;
    const ccomplex_t IT_0183 = IT_0164 + IT_0182;
    const ccomplex_t IT_0184 = IT_0146*IT_0183;
    const ccomplex_t IT_0185 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0186 = IT_0016*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = N_d4*e_em*m_tau*IT_0002*U_stau_00;
    const ccomplex_t IT_0189 = IT_0022*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*(IT_0187 + 0.5*IT_0190);
    const ccomplex_t IT_0192 = -IT_0191;
    const ccomplex_t IT_0193 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0194 = IT_0192*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*IT_0194;
    const ccomplex_t IT_0196 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0197 = IT_0046*IT_0196;
    const ccomplex_t IT_0198 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0199 = IT_0049*IT_0198;
    const ccomplex_t IT_0200 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0201 = IT_0052*IT_0200;
    const ccomplex_t IT_0202 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0203 = IT_0055*IT_0202;
    const ccomplex_t IT_0204 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0205 = IT_0046*IT_0204;
    const ccomplex_t IT_0206 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0207 = IT_0049*IT_0206;
    const ccomplex_t IT_0208 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0209 = IT_0052*IT_0208;
    const ccomplex_t IT_0210 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0211 = IT_0055*IT_0210;
    const ccomplex_t IT_0212 = (-0.5)*IT_0197 + 0.5*IT_0199 + 0.5*IT_0201 + (
      -0.5)*IT_0203 + (-0.5)*IT_0205 + 0.5*IT_0207 + 0.5*IT_0209 + (-0.5)*IT_0211;
    const ccomplex_t IT_0213 = m_N_4*IT_0212;
    const ccomplex_t IT_0214 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0215 = IT_0046*IT_0214;
    const ccomplex_t IT_0216 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0217 = IT_0049*IT_0216;
    const ccomplex_t IT_0218 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0219 = IT_0052*IT_0218;
    const ccomplex_t IT_0220 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0221 = IT_0055*IT_0220;
    const ccomplex_t IT_0222 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0223 = IT_0046*IT_0222;
    const ccomplex_t IT_0224 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0225 = IT_0049*IT_0224;
    const ccomplex_t IT_0226 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0227 = IT_0052*IT_0226;
    const ccomplex_t IT_0228 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0229 = IT_0055*IT_0228;
    const ccomplex_t IT_0230 = 0.5*IT_0215 + (-0.5)*IT_0217 + (-0.5)*IT_0219 +
       0.5*IT_0221 + 0.5*IT_0223 + (-0.5)*IT_0225 + (-0.5)*IT_0227 + 0.5*IT_0229;
    const ccomplex_t IT_0231 = m_N_2*IT_0230;
    const ccomplex_t IT_0232 = IT_0213 + IT_0231;
    const ccomplex_t IT_0233 = IT_0195*IT_0232;
    const ccomplex_t IT_0234 = (-0.5)*IT_0033 + (-0.5)*IT_0087 + (-0.5)
      *IT_0110 + (-0.5)*IT_0135 + (-0.5)*IT_0184 + (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = 2*s_13;
    const ccomplex_t IT_0236 = IT_0146*IT_0181;
    const ccomplex_t IT_0237 = IT_0195*IT_0230;
    const ccomplex_t IT_0238 = cpow(s_12 + 0.5*IT_0029 + (-0.5)*IT_0030 + 0.5
      *IT_0042 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0239 = e_em*m_tau*IT_0001*IT_0002*IT_0007*IT_0021;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0242 = IT_0016*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0245 = IT_0001*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = conj(N_d2)*e_em*m_tau*IT_0002*U_stau_10;
    const ccomplex_t IT_0248 = IT_0022*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0243 + IT_0246 + 
      -IT_0249);
    const ccomplex_t IT_0251 = (-0.5)*IT_0250;
    const ccomplex_t IT_0252 = IT_0240*IT_0251;
    const ccomplex_t IT_0253 = IT_0238*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0129*IT_0134;
    const ccomplex_t IT_0256 = IT_0045*IT_0084;
    const ccomplex_t IT_0257 = 0.5*IT_0236 + 0.5*IT_0237 + 0.5*IT_0254 + 0.5
      *IT_0255 + 0.5*IT_0256;
    const ccomplex_t IT_0258 = s_14*m_tau;
    const ccomplex_t IT_0259 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0260 = IT_0016*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0263 = IT_0001*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = conj(N_d3)*e_em*m_tau*IT_0002*U_stau_10;
    const ccomplex_t IT_0266 = IT_0022*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0261 + IT_0264 + 
      -IT_0267);
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = IT_0144*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0163*IT_0271;
    const ccomplex_t IT_0273 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0274 = IT_0016*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0277 = IT_0001*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = conj(N_d4)*e_em*m_tau*IT_0002*U_stau_10;
    const ccomplex_t IT_0280 = IT_0022*IT_0279;
    const ccomplex_t IT_0281 = 1.4142135623731*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*(IT_0275 + IT_0278 + 
      -IT_0281);
    const ccomplex_t IT_0283 = (-0.5)*IT_0282;
    const ccomplex_t IT_0284 = IT_0193*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = IT_0212*IT_0285;
    const ccomplex_t IT_0287 = IT_0132*IT_0251;
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*IT_0287;
    const ccomplex_t IT_0289 = IT_0119*IT_0288;
    const ccomplex_t IT_0290 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0291 = IT_0016*IT_0290;
    const ccomplex_t IT_0292 = 1.4142135623731*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0294 = IT_0001*IT_0293;
    const ccomplex_t IT_0295 = 1.4142135623731*IT_0294;
    const ccomplex_t IT_0296 = conj(N_d1)*e_em*m_tau*IT_0002*U_stau_10;
    const ccomplex_t IT_0297 = IT_0022*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*(IT_0292 + IT_0295 + 
      -IT_0298);
    const ccomplex_t IT_0300 = (-0.5)*IT_0299;
    const ccomplex_t IT_0301 = IT_0043*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*IT_0301;
    const ccomplex_t IT_0303 = IT_0066*IT_0302;
    const ccomplex_t IT_0304 = IT_0027*IT_0240;
    const ccomplex_t IT_0305 = IT_0238*IT_0304;
    const ccomplex_t IT_0306 = (0 + _Complex_I*1)*IT_0305;
    const ccomplex_t IT_0307 = (-0.5)*IT_0272 + (-0.5)*IT_0286 + (-0.5)
      *IT_0289 + (-0.5)*IT_0303 + 0.5*IT_0306;
    const ccomplex_t IT_0308 = s_34*m_N_2;
    const ccomplex_t IT_0309 = IT_0131*IT_0288;
    const ccomplex_t IT_0310 = m_N_2*IT_0066;
    const ccomplex_t IT_0311 = m_N_1*IT_0084;
    const ccomplex_t IT_0312 = IT_0310 + IT_0311;
    const ccomplex_t IT_0313 = IT_0302*IT_0312;
    const ccomplex_t IT_0314 = m_N_2*IT_0163;
    const ccomplex_t IT_0315 = m_N_3*IT_0181;
    const ccomplex_t IT_0316 = IT_0314 + IT_0315;
    const ccomplex_t IT_0317 = IT_0271*IT_0316;
    const ccomplex_t IT_0318 = IT_0014*IT_0251;
    const ccomplex_t IT_0319 = IT_0031*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0322 = IT_0016*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0325 = IT_0001*IT_0324;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = conj(N_d2)*e_em*m_tau*IT_0002*U_stau_11;
    const ccomplex_t IT_0328 = IT_0022*IT_0327;
    const ccomplex_t IT_0329 = 1.4142135623731*IT_0328;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*(IT_0323 + IT_0326 + 
      -IT_0329);
    const ccomplex_t IT_0331 = (-0.5)*IT_0330;
    const ccomplex_t IT_0332 = IT_0106*IT_0331;
    const ccomplex_t IT_0333 = IT_0108*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*IT_0333;
    const ccomplex_t IT_0335 = m_N_2*IT_0212;
    const ccomplex_t IT_0336 = m_N_4*IT_0230;
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = IT_0285*IT_0337;
    const ccomplex_t IT_0339 = 0.5*IT_0309 + 0.5*IT_0313 + 0.5*IT_0317 + 0.5
      *IT_0320 + 0.5*IT_0334 + 0.5*IT_0338;
    const ccomplex_t IT_0340 = 2*conj(IT_0339);
    const ccomplex_t IT_0341 = m_N_2*m_tau;
    const ccomplex_t IT_0342 = 2*IT_0234;
    const ccomplex_t IT_0343 = 2*conj(IT_0234);
    const ccomplex_t IT_0344 = pow(m_A0, 2);
    const ccomplex_t IT_0345 = m_N_2*m_tau*IT_0344;
    const ccomplex_t IT_0346 = s_13*IT_0344;
    const ccomplex_t IT_0347 = s_14*s_34;
    const ccomplex_t IT_0348 = (-2)*IT_0347;
    const ccomplex_t IT_0349 = IT_0346 + IT_0348;
    const ccomplex_t IT_0350 = 2*IT_0339;
    const ccomplex_t IT_0351 = 0.5*IT_0308;
    const ccomplex_t IT_0352 = IT_0234*conj(IT_0234)*IT_0235 + (IT_0257
      *IT_0258 + IT_0307*IT_0308)*IT_0340 + conj(IT_0339)*(IT_0235*IT_0339 +
       IT_0341*IT_0342) + (IT_0258*IT_0307 + IT_0257*IT_0308 + IT_0339*IT_0341)
      *IT_0343 + 2*conj(IT_0257)*(IT_0307*IT_0345 + -IT_0257*IT_0349 + 0.5
      *IT_0258*IT_0350 + IT_0342*IT_0351) + 2*conj(IT_0307)*(0.5*IT_0258*IT_0342
       + IT_0257*IT_0345 + -IT_0307*IT_0349 + IT_0350*IT_0351);
    return create_ccomplex_return(IT_0352);
}

