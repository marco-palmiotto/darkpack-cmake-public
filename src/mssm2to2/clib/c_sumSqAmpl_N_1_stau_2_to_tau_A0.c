#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_2_to_tau_A0.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_2_to_tau_A0(
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
    const ccomplex_t IT_0003 = e_em*mu_h*m_tau*IT_0002*U_stau_01*conj
      (U_stau_10);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_00)
      *U_stau_11;
    const ccomplex_t IT_0006 = IT_0001*IT_0005;
    const ccomplex_t IT_0007 = sin(beta);
    const ccomplex_t IT_0008 = A_tau*U_stau_01*conj(U_stau_10);
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = 1.4142135623731*IT_0009;
    const ccomplex_t IT_0011 = A_tau*conj(U_stau_00)*U_stau_11;
    const ccomplex_t IT_0012 = IT_0007*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = 0.5*IT_0004 + (-0.5)*IT_0006 + 0.5*IT_0010 + (
      -0.5)*IT_0013;
    const ccomplex_t IT_0015 = cos(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = cos(beta);
    const ccomplex_t IT_0021 = cpow(IT_0020, -1);
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = N_d1*e_em*m_tau*IT_0002*U_stau_00;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0019 + 0.5*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = pow(m_N_1, 2);
    const ccomplex_t IT_0030 = pow(m_tau, 2);
    const ccomplex_t IT_0031 = cpow(s_13 + (-0.5)*IT_0029 + (-0.5)*IT_0030 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0032 = IT_0028*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0035 = IT_0016*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_d2*e_em*m_tau*IT_0002*U_stau_01;
    const ccomplex_t IT_0038 = IT_0022*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = pow(m_stau_2, 2);
    const ccomplex_t IT_0043 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
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
    const ccomplex_t IT_0067 = m_N_2*IT_0066;
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
    const ccomplex_t IT_0085 = m_N_1*IT_0084;
    const ccomplex_t IT_0086 = IT_0067 + IT_0085;
    const ccomplex_t IT_0087 = IT_0045*IT_0086;
    const ccomplex_t IT_0088 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0089 = IT_0016*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = N_d3*e_em*m_tau*IT_0002*U_stau_01;
    const ccomplex_t IT_0092 = IT_0022*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0090 + 0.5*IT_0093);
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0100 = IT_0046*IT_0099;
    const ccomplex_t IT_0101 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0102 = IT_0049*IT_0101;
    const ccomplex_t IT_0103 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0104 = IT_0052*IT_0103;
    const ccomplex_t IT_0105 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0106 = IT_0055*IT_0105;
    const ccomplex_t IT_0107 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0108 = IT_0046*IT_0107;
    const ccomplex_t IT_0109 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0110 = IT_0049*IT_0109;
    const ccomplex_t IT_0111 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0112 = IT_0052*IT_0111;
    const ccomplex_t IT_0113 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0114 = IT_0055*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0100 + 0.5*IT_0102 + 0.5*IT_0104 + (
      -0.5)*IT_0106 + (-0.5)*IT_0108 + 0.5*IT_0110 + 0.5*IT_0112 + (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = m_N_3*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0118 = IT_0046*IT_0117;
    const ccomplex_t IT_0119 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0120 = IT_0049*IT_0119;
    const ccomplex_t IT_0121 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0122 = IT_0052*IT_0121;
    const ccomplex_t IT_0123 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0124 = IT_0055*IT_0123;
    const ccomplex_t IT_0125 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0126 = IT_0046*IT_0125;
    const ccomplex_t IT_0127 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0128 = IT_0049*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0130 = IT_0052*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0132 = IT_0055*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0118 + (-0.5)*IT_0120 + (-0.5)*IT_0122 +
       0.5*IT_0124 + 0.5*IT_0126 + (-0.5)*IT_0128 + (-0.5)*IT_0130 + 0.5*IT_0132;
    const ccomplex_t IT_0134 = m_N_1*IT_0133;
    const ccomplex_t IT_0135 = IT_0116 + IT_0134;
    const ccomplex_t IT_0136 = IT_0098*IT_0135;
    const ccomplex_t IT_0137 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0138 = IT_0016*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = N_d4*e_em*m_tau*IT_0002*U_stau_01;
    const ccomplex_t IT_0141 = IT_0022*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0139 + 0.5*IT_0142);
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0146 = IT_0144*IT_0145;
    const ccomplex_t IT_0147 = (0 + _Complex_I*1)*IT_0146;
    const ccomplex_t IT_0148 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0149 = IT_0046*IT_0148;
    const ccomplex_t IT_0150 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0151 = IT_0049*IT_0150;
    const ccomplex_t IT_0152 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0153 = IT_0052*IT_0152;
    const ccomplex_t IT_0154 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0155 = IT_0055*IT_0154;
    const ccomplex_t IT_0156 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0157 = IT_0046*IT_0156;
    const ccomplex_t IT_0158 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0159 = IT_0049*IT_0158;
    const ccomplex_t IT_0160 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0161 = IT_0052*IT_0160;
    const ccomplex_t IT_0162 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0163 = IT_0055*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0149 + 0.5*IT_0151 + 0.5*IT_0153 + (
      -0.5)*IT_0155 + (-0.5)*IT_0157 + 0.5*IT_0159 + 0.5*IT_0161 + (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = m_N_4*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0167 = IT_0046*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0169 = IT_0049*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0171 = IT_0052*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0173 = IT_0055*IT_0172;
    const ccomplex_t IT_0174 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0175 = IT_0046*IT_0174;
    const ccomplex_t IT_0176 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0177 = IT_0049*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0179 = IT_0052*IT_0178;
    const ccomplex_t IT_0180 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0181 = IT_0055*IT_0180;
    const ccomplex_t IT_0182 = 0.5*IT_0167 + (-0.5)*IT_0169 + (-0.5)*IT_0171 +
       0.5*IT_0173 + 0.5*IT_0175 + (-0.5)*IT_0177 + (-0.5)*IT_0179 + 0.5*IT_0181;
    const ccomplex_t IT_0183 = m_N_1*IT_0182;
    const ccomplex_t IT_0184 = IT_0165 + IT_0183;
    const ccomplex_t IT_0185 = IT_0147*IT_0184;
    const ccomplex_t IT_0186 = cpow(s_13 + (-0.5)*IT_0029 + (-0.5)*IT_0030 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0187 = e_em*mu_h*m_tau*IT_0002*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0188 = IT_0001*IT_0187;
    const ccomplex_t IT_0189 = e_em*mu_h*m_tau*IT_0002*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0190 = IT_0001*IT_0189;
    const ccomplex_t IT_0191 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0192 = IT_0007*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0195 = IT_0007*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (-0.5)*IT_0188 + 0.5*IT_0190 + (-0.5)*IT_0193 +
       0.5*IT_0196;
    const ccomplex_t IT_0198 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0199 = IT_0016*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = N_d1*e_em*m_tau*IT_0002*U_stau_01;
    const ccomplex_t IT_0202 = IT_0022*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0200 + 0.5*IT_0203);
    const ccomplex_t IT_0205 = -IT_0204;
    const ccomplex_t IT_0206 = IT_0197*IT_0205;
    const ccomplex_t IT_0207 = IT_0186*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0210 = IT_0046*IT_0209;
    const ccomplex_t IT_0211 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0212 = IT_0049*IT_0211;
    const ccomplex_t IT_0213 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0214 = IT_0052*IT_0213;
    const ccomplex_t IT_0215 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0216 = IT_0055*IT_0215;
    const ccomplex_t IT_0217 = -IT_0210 + IT_0212 + IT_0214 + -IT_0216;
    const ccomplex_t IT_0218 = m_N_1*IT_0217;
    const ccomplex_t IT_0219 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0220 = IT_0046*IT_0219;
    const ccomplex_t IT_0221 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0222 = IT_0049*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0224 = IT_0052*IT_0223;
    const ccomplex_t IT_0225 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0226 = IT_0055*IT_0225;
    const ccomplex_t IT_0227 = IT_0220 + -IT_0222 + -IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = m_N_1*IT_0227;
    const ccomplex_t IT_0229 = IT_0218 + IT_0228;
    const ccomplex_t IT_0230 = cpow(s_23 + (-0.5)*IT_0030 + (-0.5)*IT_0042 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0231 = IT_0205*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = IT_0229*IT_0232;
    const ccomplex_t IT_0234 = (-0.5)*IT_0033 + (-0.5)*IT_0087 + (-0.5)
      *IT_0136 + (-0.5)*IT_0185 + (-0.5)*IT_0208 + (-0.5)*IT_0233;
    const ccomplex_t IT_0235 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0236 = IT_0016*IT_0235;
    const ccomplex_t IT_0237 = 1.4142135623731*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0239 = IT_0001*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = conj(N_d1)*e_em*m_tau*IT_0002*U_stau_10;
    const ccomplex_t IT_0242 = IT_0022*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0237 + IT_0240 + 
      -IT_0243);
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = IT_0014*IT_0245;
    const ccomplex_t IT_0247 = IT_0031*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0250 = IT_0016*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0253 = IT_0001*IT_0252;
    const ccomplex_t IT_0254 = 1.4142135623731*IT_0253;
    const ccomplex_t IT_0255 = conj(N_d1)*e_em*m_tau*IT_0002*U_stau_11;
    const ccomplex_t IT_0256 = IT_0022*IT_0255;
    const ccomplex_t IT_0257 = 1.4142135623731*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*(IT_0251 + IT_0254 + 
      -IT_0257);
    const ccomplex_t IT_0259 = (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = IT_0230*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = IT_0229*IT_0261;
    const ccomplex_t IT_0263 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0264 = IT_0016*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0267 = IT_0001*IT_0266;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = conj(N_d4)*e_em*m_tau*IT_0002*U_stau_11;
    const ccomplex_t IT_0270 = IT_0022*IT_0269;
    const ccomplex_t IT_0271 = 1.4142135623731*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*(IT_0265 + IT_0268 + 
      -IT_0271);
    const ccomplex_t IT_0273 = (-0.5)*IT_0272;
    const ccomplex_t IT_0274 = IT_0145*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = m_N_1*IT_0164;
    const ccomplex_t IT_0277 = m_N_4*IT_0182;
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = IT_0275*IT_0278;
    const ccomplex_t IT_0280 = IT_0197*IT_0259;
    const ccomplex_t IT_0281 = IT_0186*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = m_N_1*IT_0066;
    const ccomplex_t IT_0284 = m_N_2*IT_0084;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0287 = IT_0016*IT_0286;
    const ccomplex_t IT_0288 = 1.4142135623731*IT_0287;
    const ccomplex_t IT_0289 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0290 = IT_0001*IT_0289;
    const ccomplex_t IT_0291 = 1.4142135623731*IT_0290;
    const ccomplex_t IT_0292 = conj(N_d2)*e_em*m_tau*IT_0002*U_stau_11;
    const ccomplex_t IT_0293 = IT_0022*IT_0292;
    const ccomplex_t IT_0294 = 1.4142135623731*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*(IT_0288 + IT_0291 + 
      -IT_0294);
    const ccomplex_t IT_0296 = (-0.5)*IT_0295;
    const ccomplex_t IT_0297 = IT_0043*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = IT_0285*IT_0298;
    const ccomplex_t IT_0300 = m_N_1*IT_0115;
    const ccomplex_t IT_0301 = m_N_3*IT_0133;
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0304 = IT_0016*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0307 = IT_0001*IT_0306;
    const ccomplex_t IT_0308 = 1.4142135623731*IT_0307;
    const ccomplex_t IT_0309 = conj(N_d3)*e_em*m_tau*IT_0002*U_stau_11;
    const ccomplex_t IT_0310 = IT_0022*IT_0309;
    const ccomplex_t IT_0311 = 1.4142135623731*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*(IT_0305 + IT_0308 + 
      -IT_0311);
    const ccomplex_t IT_0313 = (-0.5)*IT_0312;
    const ccomplex_t IT_0314 = IT_0096*IT_0313;
    const ccomplex_t IT_0315 = (0 + _Complex_I*1)*IT_0314;
    const ccomplex_t IT_0316 = IT_0302*IT_0315;
    const ccomplex_t IT_0317 = 0.5*IT_0248 + 0.5*IT_0262 + 0.5*IT_0279 + 0.5
      *IT_0282 + 0.5*IT_0299 + 0.5*IT_0316;
    const ccomplex_t IT_0318 = 2*s_13;
    const ccomplex_t IT_0319 = m_N_1*m_tau;
    const ccomplex_t IT_0320 = 2*IT_0234;
    const ccomplex_t IT_0321 = IT_0066*IT_0298;
    const ccomplex_t IT_0322 = IT_0115*IT_0315;
    const ccomplex_t IT_0323 = IT_0164*IT_0275;
    const ccomplex_t IT_0324 = cpow(s_12 + 0.5*IT_0029 + (-0.5)*IT_0030 + 0.5
      *IT_0042 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0325 = e_em*m_tau*IT_0001*IT_0002*IT_0007*IT_0021;
    const ccomplex_t IT_0326 = (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = IT_0205*IT_0326;
    const ccomplex_t IT_0328 = IT_0324*IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*IT_0328;
    const ccomplex_t IT_0330 = IT_0217*IT_0261;
    const ccomplex_t IT_0331 = (-0.5)*IT_0321 + (-0.5)*IT_0322 + (-0.5)
      *IT_0323 + 0.5*IT_0329 + (-0.5)*IT_0330;
    const ccomplex_t IT_0332 = s_14*m_tau;
    const ccomplex_t IT_0333 = IT_0045*IT_0084;
    const ccomplex_t IT_0334 = IT_0098*IT_0133;
    const ccomplex_t IT_0335 = IT_0147*IT_0182;
    const ccomplex_t IT_0336 = IT_0227*IT_0232;
    const ccomplex_t IT_0337 = IT_0259*IT_0326;
    const ccomplex_t IT_0338 = IT_0324*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*IT_0338;
    const ccomplex_t IT_0340 = 0.5*IT_0333 + 0.5*IT_0334 + 0.5*IT_0335 + 0.5
      *IT_0336 + 0.5*IT_0339;
    const ccomplex_t IT_0341 = s_34*m_N_1;
    const ccomplex_t IT_0342 = 2*conj(IT_0234);
    const ccomplex_t IT_0343 = 2*conj(IT_0317);
    const ccomplex_t IT_0344 = pow(m_A0, 2);
    const ccomplex_t IT_0345 = m_N_1*m_tau*IT_0344;
    const ccomplex_t IT_0346 = s_13*IT_0344;
    const ccomplex_t IT_0347 = s_14*s_34;
    const ccomplex_t IT_0348 = (-2)*IT_0347;
    const ccomplex_t IT_0349 = IT_0346 + IT_0348;
    const ccomplex_t IT_0350 = 2*IT_0317;
    const ccomplex_t IT_0351 = 0.5*IT_0341;
    const ccomplex_t IT_0352 = (IT_0234*conj(IT_0234) + IT_0317*conj(IT_0317))
      *IT_0318 + conj(IT_0317)*IT_0319*IT_0320 + (IT_0317*IT_0319 + IT_0331
      *IT_0332 + IT_0340*IT_0341)*IT_0342 + (IT_0332*IT_0340 + IT_0331*IT_0341)
      *IT_0343 + 2*conj(IT_0340)*(IT_0331*IT_0345 + -IT_0340*IT_0349 + 0.5
      *IT_0332*IT_0350 + IT_0320*IT_0351) + 2*conj(IT_0331)*(0.5*IT_0320*IT_0332
       + IT_0340*IT_0345 + -IT_0331*IT_0349 + IT_0350*IT_0351);
    return create_ccomplex_return(IT_0352);
}

