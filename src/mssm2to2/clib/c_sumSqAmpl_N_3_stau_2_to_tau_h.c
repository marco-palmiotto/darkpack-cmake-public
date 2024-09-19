#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_stau_2_to_tau_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_stau_2_to_tau_h(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = cos(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = pow(m_W, -1);
    const ccomplex_t IT_0011 = N_d2*e_em*m_tau*IT_0010*U_stau_01;
    const ccomplex_t IT_0012 = IT_0009*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0004 + 0.5*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = pow(m_tau, 2);
    const ccomplex_t IT_0017 = pow(m_stau_2, 2);
    const ccomplex_t IT_0018 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0017 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0019 = IT_0015*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*IT_0019;
    const ccomplex_t IT_0021 = cos(alpha);
    const ccomplex_t IT_0022 = IT_0001*IT_0021;
    const ccomplex_t IT_0023 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0024 = IT_0022*IT_0023;
    const ccomplex_t IT_0025 = sin(alpha);
    const ccomplex_t IT_0026 = IT_0001*IT_0025;
    const ccomplex_t IT_0027 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0028 = IT_0026*IT_0027;
    const ccomplex_t IT_0029 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0024 + IT_0028 + IT_0030 + IT_0032;
    const ccomplex_t IT_0034 = IT_0008*IT_0021;
    const ccomplex_t IT_0035 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0008*IT_0025;
    const ccomplex_t IT_0038 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0041 = IT_0034*IT_0040;
    const ccomplex_t IT_0042 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0043 = IT_0037*IT_0042;
    const ccomplex_t IT_0044 = -IT_0036 + -IT_0039 + -IT_0041 + -IT_0043;
    const ccomplex_t IT_0045 = IT_0033 + IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*IT_0045;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = m_N_2*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0050 = IT_0022*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0052 = IT_0026*IT_0051;
    const ccomplex_t IT_0053 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0054 = IT_0022*IT_0053;
    const ccomplex_t IT_0055 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0056 = IT_0026*IT_0055;
    const ccomplex_t IT_0057 = IT_0050 + IT_0052 + IT_0054 + IT_0056;
    const ccomplex_t IT_0058 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0059 = IT_0034*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0061 = IT_0037*IT_0060;
    const ccomplex_t IT_0062 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0063 = IT_0034*IT_0062;
    const ccomplex_t IT_0064 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0065 = IT_0037*IT_0064;
    const ccomplex_t IT_0066 = -IT_0059 + -IT_0061 + -IT_0063 + -IT_0065;
    const ccomplex_t IT_0067 = IT_0057 + IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*IT_0067;
    const ccomplex_t IT_0069 = (-0.5)*IT_0068;
    const ccomplex_t IT_0070 = m_N_3*IT_0069;
    const ccomplex_t IT_0071 = IT_0048 + IT_0070;
    const ccomplex_t IT_0072 = IT_0020*IT_0071;
    const ccomplex_t IT_0073 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0074 = IT_0022*IT_0073;
    const ccomplex_t IT_0075 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0076 = IT_0026*IT_0075;
    const ccomplex_t IT_0077 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0078 = IT_0034*IT_0077;
    const ccomplex_t IT_0079 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0080 = IT_0037*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0074 + IT_0076 + 
      -IT_0078 + -IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = m_N_3*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0085 = IT_0022*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0087 = IT_0026*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0089 = IT_0034*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0091 = IT_0037*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + IT_0087 + 
      -IT_0089 + -IT_0091);
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = m_N_3*IT_0093;
    const ccomplex_t IT_0095 = IT_0083 + IT_0094;
    const ccomplex_t IT_0096 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0097 = IT_0001*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = N_d3*e_em*m_tau*IT_0010*U_stau_01;
    const ccomplex_t IT_0100 = IT_0009*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + 0.5*IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0017 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0109 = IT_0001*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = N_d1*e_em*m_tau*IT_0010*U_stau_01;
    const ccomplex_t IT_0112 = IT_0009*IT_0111;
    const ccomplex_t IT_0113 = 1.4142135623731*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0110 + 0.5*IT_0113);
    const ccomplex_t IT_0115 = -IT_0114;
    const ccomplex_t IT_0116 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0017 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0117 = IT_0115*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0120 = IT_0022*IT_0119;
    const ccomplex_t IT_0121 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0122 = IT_0026*IT_0121;
    const ccomplex_t IT_0123 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0124 = IT_0022*IT_0123;
    const ccomplex_t IT_0125 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0126 = IT_0026*IT_0125;
    const ccomplex_t IT_0127 = IT_0120 + IT_0122 + IT_0124 + IT_0126;
    const ccomplex_t IT_0128 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0129 = IT_0034*IT_0128;
    const ccomplex_t IT_0130 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0131 = IT_0037*IT_0130;
    const ccomplex_t IT_0132 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0133 = IT_0034*IT_0132;
    const ccomplex_t IT_0134 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0135 = IT_0037*IT_0134;
    const ccomplex_t IT_0136 = -IT_0129 + -IT_0131 + -IT_0133 + -IT_0135;
    const ccomplex_t IT_0137 = IT_0127 + IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = (-0.5)*IT_0138;
    const ccomplex_t IT_0140 = m_N_1*IT_0139;
    const ccomplex_t IT_0141 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0142 = IT_0022*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0144 = IT_0026*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0146 = IT_0022*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0148 = IT_0026*IT_0147;
    const ccomplex_t IT_0149 = IT_0142 + IT_0144 + IT_0146 + IT_0148;
    const ccomplex_t IT_0150 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0151 = IT_0034*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0153 = IT_0037*IT_0152;
    const ccomplex_t IT_0154 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0155 = IT_0034*IT_0154;
    const ccomplex_t IT_0156 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0157 = IT_0037*IT_0156;
    const ccomplex_t IT_0158 = -IT_0151 + -IT_0153 + -IT_0155 + -IT_0157;
    const ccomplex_t IT_0159 = IT_0149 + IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = (-0.5)*IT_0160;
    const ccomplex_t IT_0162 = m_N_3*IT_0161;
    const ccomplex_t IT_0163 = IT_0140 + IT_0162;
    const ccomplex_t IT_0164 = IT_0118*IT_0163;
    const ccomplex_t IT_0165 = pow(m_N_3, 2);
    const ccomplex_t IT_0166 = cpow(s_13 + (-0.5)*IT_0016 + (-0.5)*IT_0165 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0167 = mu_h*m_tau*IT_0010;
    const ccomplex_t IT_0168 = IT_0005*IT_0025;
    const ccomplex_t IT_0169 = sin(beta);
    const ccomplex_t IT_0170 = IT_0021*IT_0169;
    const ccomplex_t IT_0171 = IT_0168 + IT_0170;
    const ccomplex_t IT_0172 = cpow(IT_0000, -2);
    const ccomplex_t IT_0173 = 2*U_stau_01*((IT_0006*IT_0008*IT_0010*IT_0016
      *IT_0025 + (-0.5)*m_W*IT_0171*(IT_0008 + -IT_0007*IT_0172))*conj(U_stau_01
      ) + 1./2*IT_0006*IT_0008*IT_0021*IT_0167*conj(U_stau_11)) + U_stau_11*
      (IT_0006*IT_0008*IT_0021*IT_0167*conj(U_stau_01) + 2*(IT_0006*IT_0008
      *IT_0010*IT_0016*IT_0025 + -m_W*IT_0007*IT_0171*IT_0172)*conj(U_stau_11));
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(e_em*IT_0173 +
       1.4142135623731*A_tau*IT_0025*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0175 = 0.5*IT_0174;
    const ccomplex_t IT_0176 = IT_0103*IT_0175;
    const ccomplex_t IT_0177 = IT_0166*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = cpow(s_13 + (-0.5)*IT_0016 + (-0.5)*IT_0165 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0180 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0181 = IT_0001*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = N_d3*e_em*m_tau*IT_0010*U_stau_00;
    const ccomplex_t IT_0184 = IT_0009*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + 0.5*IT_0185);
    const ccomplex_t IT_0187 = -IT_0186;
    const ccomplex_t IT_0188 = 2*conj(U_stau_00)*((IT_0006*IT_0008*IT_0010
      *IT_0016*IT_0025 + (-0.5)*m_W*IT_0171*(IT_0008 + -IT_0007*IT_0172))
      *U_stau_01 + 1./2*IT_0006*IT_0008*IT_0021*IT_0167*U_stau_11) + conj
      (U_stau_10)*(IT_0006*IT_0008*IT_0021*IT_0167*U_stau_01 + 2*(IT_0006
      *IT_0008*IT_0010*IT_0016*IT_0025 + -m_W*IT_0007*IT_0171*IT_0172)*U_stau_11);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(e_em*IT_0188 +
       1.4142135623731*A_tau*IT_0025*(U_stau_01*conj(U_stau_10) + conj(U_stau_00
      )*U_stau_11));
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = IT_0187*IT_0190;
    const ccomplex_t IT_0192 = IT_0179*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0195 = IT_0001*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = N_d4*e_em*m_tau*IT_0010*U_stau_01;
    const ccomplex_t IT_0198 = IT_0009*IT_0197;
    const ccomplex_t IT_0199 = 1.4142135623731*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*(IT_0196 + 0.5*IT_0199);
    const ccomplex_t IT_0201 = -IT_0200;
    const ccomplex_t IT_0202 = cpow(s_23 + (-0.5)*IT_0016 + (-0.5)*IT_0017 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0203 = IT_0201*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0206 = IT_0022*IT_0205;
    const ccomplex_t IT_0207 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0208 = IT_0026*IT_0207;
    const ccomplex_t IT_0209 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0210 = IT_0022*IT_0209;
    const ccomplex_t IT_0211 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0212 = IT_0026*IT_0211;
    const ccomplex_t IT_0213 = IT_0206 + IT_0208 + IT_0210 + IT_0212;
    const ccomplex_t IT_0214 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0215 = IT_0034*IT_0214;
    const ccomplex_t IT_0216 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0217 = IT_0037*IT_0216;
    const ccomplex_t IT_0218 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0219 = IT_0034*IT_0218;
    const ccomplex_t IT_0220 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0221 = IT_0037*IT_0220;
    const ccomplex_t IT_0222 = -IT_0215 + -IT_0217 + -IT_0219 + -IT_0221;
    const ccomplex_t IT_0223 = IT_0213 + IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = (-0.5)*IT_0224;
    const ccomplex_t IT_0226 = m_N_4*IT_0225;
    const ccomplex_t IT_0227 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0228 = IT_0022*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0230 = IT_0026*IT_0229;
    const ccomplex_t IT_0231 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0232 = IT_0022*IT_0231;
    const ccomplex_t IT_0233 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0234 = IT_0026*IT_0233;
    const ccomplex_t IT_0235 = IT_0228 + IT_0230 + IT_0232 + IT_0234;
    const ccomplex_t IT_0236 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0237 = IT_0034*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0239 = IT_0037*IT_0238;
    const ccomplex_t IT_0240 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0241 = IT_0034*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0243 = IT_0037*IT_0242;
    const ccomplex_t IT_0244 = -IT_0237 + -IT_0239 + -IT_0241 + -IT_0243;
    const ccomplex_t IT_0245 = IT_0235 + IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = m_N_3*IT_0247;
    const ccomplex_t IT_0249 = IT_0226 + IT_0248;
    const ccomplex_t IT_0250 = IT_0204*IT_0249;
    const ccomplex_t IT_0251 = cpow(s_12 + (-0.5)*IT_0016 + 0.5*IT_0017 + 0.5
      *IT_0165 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0008
      *IT_0010*IT_0025;
    const ccomplex_t IT_0253 = 0.5*IT_0252;
    const ccomplex_t IT_0254 = IT_0103*IT_0253;
    const ccomplex_t IT_0255 = IT_0251*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = m_tau*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0072 + (-0.5)*IT_0107 + (-0.5)
      *IT_0164 + (-0.5)*IT_0178 + (-0.5)*IT_0193 + (-0.5)*IT_0250 + IT_0257;
    const ccomplex_t IT_0259 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0260 = IT_0001*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0263 = IT_0008*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = conj(N_d3)*e_em*m_tau*IT_0010*U_stau_11;
    const ccomplex_t IT_0266 = IT_0009*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0261 + IT_0264 + 
      -IT_0267);
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = 2*U_stau_01*((IT_0006*IT_0008*IT_0010*IT_0016
      *IT_0025 + (-0.5)*m_W*IT_0171*(IT_0008 + -IT_0007*IT_0172))*conj(U_stau_01
      ) + 1./2*IT_0006*IT_0008*IT_0021*IT_0167*conj(U_stau_11)) + U_stau_11*
      (IT_0006*IT_0008*IT_0021*IT_0167*conj(U_stau_01) + 2*(IT_0006*IT_0008
      *IT_0010*IT_0016*IT_0025 + -m_W*IT_0007*IT_0171*IT_0172)*conj(U_stau_11));
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*(e_em*IT_0270 +
       1.4142135623731*A_tau*IT_0025*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0272 = 0.5*IT_0271;
    const ccomplex_t IT_0273 = IT_0269*IT_0272;
    const ccomplex_t IT_0274 = IT_0166*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0277 = IT_0001*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0280 = IT_0008*IT_0279;
    const ccomplex_t IT_0281 = 1.4142135623731*IT_0280;
    const ccomplex_t IT_0282 = conj(N_d3)*e_em*m_tau*IT_0010*U_stau_10;
    const ccomplex_t IT_0283 = IT_0009*IT_0282;
    const ccomplex_t IT_0284 = 1.4142135623731*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*(IT_0278 + IT_0281 + 
      -IT_0284);
    const ccomplex_t IT_0286 = (-0.5)*IT_0285;
    const ccomplex_t IT_0287 = 2*conj(U_stau_00)*((IT_0006*IT_0008*IT_0010
      *IT_0016*IT_0025 + (-0.5)*m_W*IT_0171*(IT_0008 + -IT_0007*IT_0172))
      *U_stau_01 + 1./2*IT_0006*IT_0008*IT_0021*IT_0167*U_stau_11) + conj
      (U_stau_10)*(IT_0006*IT_0008*IT_0021*IT_0167*U_stau_01 + 2*(IT_0006
      *IT_0008*IT_0010*IT_0016*IT_0025 + -m_W*IT_0007*IT_0171*IT_0172)*U_stau_11);
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*(e_em*IT_0287 +
       1.4142135623731*A_tau*IT_0025*(U_stau_01*conj(U_stau_10) + conj(U_stau_00
      )*U_stau_11));
    const ccomplex_t IT_0289 = 0.5*IT_0288;
    const ccomplex_t IT_0290 = IT_0286*IT_0289;
    const ccomplex_t IT_0291 = IT_0179*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = IT_0104*IT_0269;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = IT_0095*IT_0294;
    const ccomplex_t IT_0296 = m_N_3*IT_0139;
    const ccomplex_t IT_0297 = m_N_1*IT_0161;
    const ccomplex_t IT_0298 = IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0300 = IT_0001*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0303 = IT_0008*IT_0302;
    const ccomplex_t IT_0304 = 1.4142135623731*IT_0303;
    const ccomplex_t IT_0305 = conj(N_d1)*e_em*m_tau*IT_0010*U_stau_11;
    const ccomplex_t IT_0306 = IT_0009*IT_0305;
    const ccomplex_t IT_0307 = 1.4142135623731*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*(IT_0301 + IT_0304 + 
      -IT_0307);
    const ccomplex_t IT_0309 = (-0.5)*IT_0308;
    const ccomplex_t IT_0310 = IT_0116*IT_0309;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*IT_0310;
    const ccomplex_t IT_0312 = IT_0298*IT_0311;
    const ccomplex_t IT_0313 = m_N_3*IT_0047;
    const ccomplex_t IT_0314 = m_N_2*IT_0069;
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0317 = IT_0001*IT_0316;
    const ccomplex_t IT_0318 = 1.4142135623731*IT_0317;
    const ccomplex_t IT_0319 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0320 = IT_0008*IT_0319;
    const ccomplex_t IT_0321 = 1.4142135623731*IT_0320;
    const ccomplex_t IT_0322 = conj(N_d2)*e_em*m_tau*IT_0010*U_stau_11;
    const ccomplex_t IT_0323 = IT_0009*IT_0322;
    const ccomplex_t IT_0324 = 1.4142135623731*IT_0323;
    const ccomplex_t IT_0325 = (0 + _Complex_I*1)*(IT_0318 + IT_0321 + 
      -IT_0324);
    const ccomplex_t IT_0326 = (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = IT_0018*IT_0326;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*IT_0327;
    const ccomplex_t IT_0329 = IT_0315*IT_0328;
    const ccomplex_t IT_0330 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0331 = IT_0001*IT_0330;
    const ccomplex_t IT_0332 = 1.4142135623731*IT_0331;
    const ccomplex_t IT_0333 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0334 = IT_0008*IT_0333;
    const ccomplex_t IT_0335 = 1.4142135623731*IT_0334;
    const ccomplex_t IT_0336 = conj(N_d4)*e_em*m_tau*IT_0010*U_stau_11;
    const ccomplex_t IT_0337 = IT_0009*IT_0336;
    const ccomplex_t IT_0338 = 1.4142135623731*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*(IT_0332 + IT_0335 + 
      -IT_0338);
    const ccomplex_t IT_0340 = (-0.5)*IT_0339;
    const ccomplex_t IT_0341 = IT_0202*IT_0340;
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*IT_0341;
    const ccomplex_t IT_0343 = m_N_3*IT_0225;
    const ccomplex_t IT_0344 = m_N_4*IT_0247;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = IT_0342*IT_0345;
    const ccomplex_t IT_0347 = IT_0253*IT_0269;
    const ccomplex_t IT_0348 = IT_0251*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = m_tau*IT_0349;
    const ccomplex_t IT_0351 = 0.5*IT_0275 + 0.5*IT_0292 + 0.5*IT_0295 + 0.5
      *IT_0312 + 0.5*IT_0329 + 0.5*IT_0346 + -IT_0350;
    const ccomplex_t IT_0352 = 2*s_13;
    const ccomplex_t IT_0353 = m_N_3*m_tau;
    const ccomplex_t IT_0354 = 2*IT_0258;
    const ccomplex_t IT_0355 = IT_0082*IT_0294;
    const ccomplex_t IT_0356 = IT_0139*IT_0311;
    const ccomplex_t IT_0357 = IT_0047*IT_0328;
    const ccomplex_t IT_0358 = IT_0225*IT_0342;
    const ccomplex_t IT_0359 = 0.5*IT_0256 + (-0.5)*IT_0355 + (-0.5)*IT_0356 +
       (-0.5)*IT_0357 + (-0.5)*IT_0358;
    const ccomplex_t IT_0360 = s_14*m_tau;
    const ccomplex_t IT_0361 = IT_0020*IT_0069;
    const ccomplex_t IT_0362 = IT_0093*IT_0106;
    const ccomplex_t IT_0363 = IT_0118*IT_0161;
    const ccomplex_t IT_0364 = IT_0204*IT_0247;
    const ccomplex_t IT_0365 = (-0.5)*IT_0349 + 0.5*IT_0361 + 0.5*IT_0362 +
       0.5*IT_0363 + 0.5*IT_0364;
    const ccomplex_t IT_0366 = s_34*m_N_3;
    const ccomplex_t IT_0367 = 2*conj(IT_0258);
    const ccomplex_t IT_0368 = 2*conj(IT_0351);
    const ccomplex_t IT_0369 = pow(m_h, 2);
    const ccomplex_t IT_0370 = m_N_3*m_tau*IT_0369;
    const ccomplex_t IT_0371 = s_13*IT_0369;
    const ccomplex_t IT_0372 = s_14*s_34;
    const ccomplex_t IT_0373 = (-2)*IT_0372;
    const ccomplex_t IT_0374 = IT_0371 + IT_0373;
    const ccomplex_t IT_0375 = 2*IT_0351;
    const ccomplex_t IT_0376 = 0.5*IT_0366;
    const ccomplex_t IT_0377 = (IT_0258*conj(IT_0258) + IT_0351*conj(IT_0351))
      *IT_0352 + conj(IT_0351)*IT_0353*IT_0354 + (IT_0351*IT_0353 + IT_0359
      *IT_0360 + IT_0365*IT_0366)*IT_0367 + (IT_0360*IT_0365 + IT_0359*IT_0366)
      *IT_0368 + 2*conj(IT_0365)*(IT_0359*IT_0370 + -IT_0365*IT_0374 + 0.5
      *IT_0360*IT_0375 + IT_0354*IT_0376) + 2*conj(IT_0359)*(0.5*IT_0354*IT_0360
       + IT_0365*IT_0370 + -IT_0359*IT_0374 + IT_0375*IT_0376);
    return create_ccomplex_return(IT_0377);
}

