#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_1_to_tau_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_1_to_tau_h(
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
    const ccomplex_t IT_0000 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0001 = cos(alpha);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = IT_0000*IT_0004;
    const ccomplex_t IT_0006 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0007 = sin(alpha);
    const ccomplex_t IT_0008 = IT_0003*IT_0007;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0011 = sin(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = IT_0001*IT_0012;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0016 = IT_0007*IT_0012;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0005 + IT_0009 + 
      -IT_0014 + -IT_0017);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = m_N_1*IT_0019;
    const ccomplex_t IT_0021 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0022 = IT_0004*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0028 = IT_0016*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0022 + IT_0024 + 
      -IT_0026 + -IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = m_N_1*IT_0030;
    const ccomplex_t IT_0032 = IT_0020 + IT_0031;
    const ccomplex_t IT_0033 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0034 = IT_0003*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = cos(beta);
    const ccomplex_t IT_0037 = cpow(IT_0036, -1);
    const ccomplex_t IT_0038 = IT_0012*IT_0037;
    const ccomplex_t IT_0039 = pow(m_W, -1);
    const ccomplex_t IT_0040 = N_d1*e_em*m_tau*IT_0039*U_stau_00;
    const ccomplex_t IT_0041 = IT_0038*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0035 + 0.5*IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = pow(m_tau, 2);
    const ccomplex_t IT_0046 = pow(m_stau_1, 2);
    const ccomplex_t IT_0047 = cpow(s_23 + (-0.5)*IT_0045 + (-0.5)*IT_0046 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0044*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = IT_0032*IT_0049;
    const ccomplex_t IT_0051 = pow(m_N_1, 2);
    const ccomplex_t IT_0052 = cpow(s_12 + (-0.5)*IT_0045 + 0.5*IT_0046 + 0.5
      *IT_0051 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*e_em*m_tau*IT_0007*IT_0012
      *IT_0037*IT_0039;
    const ccomplex_t IT_0054 = 0.5*IT_0053;
    const ccomplex_t IT_0055 = IT_0044*IT_0054;
    const ccomplex_t IT_0056 = IT_0052*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_tau*IT_0057;
    const ccomplex_t IT_0059 = cpow(s_13 + (-0.5)*IT_0045 + (-0.5)*IT_0051 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0060 = mu_h*m_tau*IT_0039;
    const ccomplex_t IT_0061 = IT_0007*IT_0036;
    const ccomplex_t IT_0062 = sin(beta);
    const ccomplex_t IT_0063 = IT_0001*IT_0062;
    const ccomplex_t IT_0064 = IT_0061 + IT_0063;
    const ccomplex_t IT_0065 = cpow(IT_0002, -2);
    const ccomplex_t IT_0066 = 2*U_stau_00*((IT_0007*IT_0012*IT_0037*IT_0039
      *IT_0045 + (-0.5)*m_W*IT_0064*(IT_0012 + -IT_0011*IT_0065))*conj(U_stau_00
      ) + 1./2*IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_10)) + U_stau_10*
      (IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_00) + 2*(IT_0007*IT_0012
      *IT_0037*IT_0039*IT_0045 + -m_W*IT_0011*IT_0064*IT_0065)*conj(U_stau_10));
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*(e_em*IT_0066 +
       1.4142135623731*A_tau*IT_0007*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0044*IT_0068;
    const ccomplex_t IT_0070 = IT_0059*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = cpow(s_13 + (-0.5)*IT_0045 + (-0.5)*IT_0051 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0073 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0074 = IT_0003*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_d1*e_em*m_tau*IT_0039*U_stau_01;
    const ccomplex_t IT_0077 = IT_0038*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*(IT_0075 + 0.5*IT_0078);
    const ccomplex_t IT_0080 = -IT_0079;
    const ccomplex_t IT_0081 = 2*U_stau_00*((IT_0007*IT_0012*IT_0037*IT_0039
      *IT_0045 + (-0.5)*m_W*IT_0064*(IT_0012 + -IT_0011*IT_0065))*conj(U_stau_01
      ) + 1./2*IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_11)) + U_stau_10*
      (IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_01) + 2*(IT_0007*IT_0012
      *IT_0037*IT_0039*IT_0045 + -m_W*IT_0011*IT_0064*IT_0065)*conj(U_stau_11));
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(e_em*IT_0081 +
       1.4142135623731*A_tau*IT_0007*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0080*IT_0083;
    const ccomplex_t IT_0085 = IT_0072*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0088 = IT_0003*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = N_d2*e_em*m_tau*IT_0039*U_stau_00;
    const ccomplex_t IT_0091 = IT_0038*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + 0.5*IT_0092);
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = cpow(s_23 + (-0.5)*IT_0045 + (-0.5)*IT_0046 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0094*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0099 = IT_0004*IT_0098;
    const ccomplex_t IT_0100 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0101 = IT_0008*IT_0100;
    const ccomplex_t IT_0102 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0103 = IT_0004*IT_0102;
    const ccomplex_t IT_0104 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0105 = IT_0008*IT_0104;
    const ccomplex_t IT_0106 = IT_0099 + IT_0101 + IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0108 = IT_0013*IT_0107;
    const ccomplex_t IT_0109 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0110 = IT_0016*IT_0109;
    const ccomplex_t IT_0111 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0112 = IT_0013*IT_0111;
    const ccomplex_t IT_0113 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0114 = IT_0016*IT_0113;
    const ccomplex_t IT_0115 = -IT_0108 + -IT_0110 + -IT_0112 + -IT_0114;
    const ccomplex_t IT_0116 = IT_0106 + IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = m_N_2*IT_0118;
    const ccomplex_t IT_0120 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0121 = IT_0004*IT_0120;
    const ccomplex_t IT_0122 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0123 = IT_0008*IT_0122;
    const ccomplex_t IT_0124 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0125 = IT_0004*IT_0124;
    const ccomplex_t IT_0126 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0127 = IT_0008*IT_0126;
    const ccomplex_t IT_0128 = IT_0121 + IT_0123 + IT_0125 + IT_0127;
    const ccomplex_t IT_0129 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0130 = IT_0013*IT_0129;
    const ccomplex_t IT_0131 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0132 = IT_0016*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0134 = IT_0013*IT_0133;
    const ccomplex_t IT_0135 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0136 = IT_0016*IT_0135;
    const ccomplex_t IT_0137 = -IT_0130 + -IT_0132 + -IT_0134 + -IT_0136;
    const ccomplex_t IT_0138 = IT_0128 + IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = m_N_1*IT_0140;
    const ccomplex_t IT_0142 = IT_0119 + IT_0141;
    const ccomplex_t IT_0143 = IT_0097*IT_0142;
    const ccomplex_t IT_0144 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0145 = IT_0003*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = N_d3*e_em*m_tau*IT_0039*U_stau_00;
    const ccomplex_t IT_0148 = IT_0038*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + 0.5*IT_0149);
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = cpow(s_23 + (-0.5)*IT_0045 + (-0.5)*IT_0046 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0156 = IT_0004*IT_0155;
    const ccomplex_t IT_0157 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0158 = IT_0008*IT_0157;
    const ccomplex_t IT_0159 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0160 = IT_0004*IT_0159;
    const ccomplex_t IT_0161 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0162 = IT_0008*IT_0161;
    const ccomplex_t IT_0163 = IT_0156 + IT_0158 + IT_0160 + IT_0162;
    const ccomplex_t IT_0164 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0165 = IT_0013*IT_0164;
    const ccomplex_t IT_0166 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0167 = IT_0016*IT_0166;
    const ccomplex_t IT_0168 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0169 = IT_0013*IT_0168;
    const ccomplex_t IT_0170 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0171 = IT_0016*IT_0170;
    const ccomplex_t IT_0172 = -IT_0165 + -IT_0167 + -IT_0169 + -IT_0171;
    const ccomplex_t IT_0173 = IT_0163 + IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = m_N_3*IT_0175;
    const ccomplex_t IT_0177 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0178 = IT_0004*IT_0177;
    const ccomplex_t IT_0179 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0180 = IT_0008*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0182 = IT_0004*IT_0181;
    const ccomplex_t IT_0183 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0184 = IT_0008*IT_0183;
    const ccomplex_t IT_0185 = IT_0178 + IT_0180 + IT_0182 + IT_0184;
    const ccomplex_t IT_0186 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0187 = IT_0013*IT_0186;
    const ccomplex_t IT_0188 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0189 = IT_0016*IT_0188;
    const ccomplex_t IT_0190 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0191 = IT_0013*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0193 = IT_0016*IT_0192;
    const ccomplex_t IT_0194 = -IT_0187 + -IT_0189 + -IT_0191 + -IT_0193;
    const ccomplex_t IT_0195 = IT_0185 + IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*IT_0195;
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = m_N_1*IT_0197;
    const ccomplex_t IT_0199 = IT_0176 + IT_0198;
    const ccomplex_t IT_0200 = IT_0154*IT_0199;
    const ccomplex_t IT_0201 = N_B4*N_u1*e_em;
    const ccomplex_t IT_0202 = IT_0004*IT_0201;
    const ccomplex_t IT_0203 = N_B4*N_d1*e_em;
    const ccomplex_t IT_0204 = IT_0008*IT_0203;
    const ccomplex_t IT_0205 = N_B1*N_u4*e_em;
    const ccomplex_t IT_0206 = IT_0004*IT_0205;
    const ccomplex_t IT_0207 = N_B1*N_d4*e_em;
    const ccomplex_t IT_0208 = IT_0008*IT_0207;
    const ccomplex_t IT_0209 = IT_0202 + IT_0204 + IT_0206 + IT_0208;
    const ccomplex_t IT_0210 = N_W4*N_u1*e_em;
    const ccomplex_t IT_0211 = IT_0013*IT_0210;
    const ccomplex_t IT_0212 = N_W4*N_d1*e_em;
    const ccomplex_t IT_0213 = IT_0016*IT_0212;
    const ccomplex_t IT_0214 = N_W1*N_u4*e_em;
    const ccomplex_t IT_0215 = IT_0013*IT_0214;
    const ccomplex_t IT_0216 = N_W1*N_d4*e_em;
    const ccomplex_t IT_0217 = IT_0016*IT_0216;
    const ccomplex_t IT_0218 = -IT_0211 + -IT_0213 + -IT_0215 + -IT_0217;
    const ccomplex_t IT_0219 = IT_0209 + IT_0218;
    const ccomplex_t IT_0220 = (0 + _Complex_I*1)*IT_0219;
    const ccomplex_t IT_0221 = (-0.5)*IT_0220;
    const ccomplex_t IT_0222 = m_N_4*IT_0221;
    const ccomplex_t IT_0223 = conj(N_B1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0224 = IT_0004*IT_0223;
    const ccomplex_t IT_0225 = conj(N_B1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0226 = IT_0008*IT_0225;
    const ccomplex_t IT_0227 = conj(N_B4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0228 = IT_0004*IT_0227;
    const ccomplex_t IT_0229 = conj(N_B4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0230 = IT_0008*IT_0229;
    const ccomplex_t IT_0231 = IT_0224 + IT_0226 + IT_0228 + IT_0230;
    const ccomplex_t IT_0232 = conj(N_W1)*conj(N_u4)*e_em;
    const ccomplex_t IT_0233 = IT_0013*IT_0232;
    const ccomplex_t IT_0234 = conj(N_W1)*conj(N_d4)*e_em;
    const ccomplex_t IT_0235 = IT_0016*IT_0234;
    const ccomplex_t IT_0236 = conj(N_W4)*conj(N_u1)*e_em;
    const ccomplex_t IT_0237 = IT_0013*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W4)*conj(N_d1)*e_em;
    const ccomplex_t IT_0239 = IT_0016*IT_0238;
    const ccomplex_t IT_0240 = -IT_0233 + -IT_0235 + -IT_0237 + -IT_0239;
    const ccomplex_t IT_0241 = IT_0231 + IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = m_N_1*IT_0243;
    const ccomplex_t IT_0245 = IT_0222 + IT_0244;
    const ccomplex_t IT_0246 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0247 = IT_0003*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = N_d4*e_em*m_tau*IT_0039*U_stau_00;
    const ccomplex_t IT_0250 = IT_0038*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*(IT_0248 + 0.5*IT_0251);
    const ccomplex_t IT_0253 = -IT_0252;
    const ccomplex_t IT_0254 = cpow(s_23 + (-0.5)*IT_0045 + (-0.5)*IT_0046 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0255 = IT_0253*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = IT_0245*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0050 + IT_0058 + (-0.5)*IT_0071 + (
      -0.5)*IT_0086 + (-0.5)*IT_0143 + (-0.5)*IT_0200 + (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0260 = IT_0003*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0263 = IT_0012*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = conj(N_d1)*e_em*m_tau*IT_0039*U_stau_10;
    const ccomplex_t IT_0266 = IT_0038*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0261 + IT_0264 + 
      -IT_0267);
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = IT_0054*IT_0269;
    const ccomplex_t IT_0271 = IT_0052*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = m_tau*IT_0272;
    const ccomplex_t IT_0274 = 2*U_stau_00*((IT_0007*IT_0012*IT_0037*IT_0039
      *IT_0045 + (-0.5)*m_W*IT_0064*(IT_0012 + -IT_0011*IT_0065))*conj(U_stau_00
      ) + 1./2*IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_10)) + U_stau_10*
      (IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_00) + 2*(IT_0007*IT_0012
      *IT_0037*IT_0039*IT_0045 + -m_W*IT_0011*IT_0064*IT_0065)*conj(U_stau_10));
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*(e_em*IT_0274 +
       1.4142135623731*A_tau*IT_0007*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0276 = 0.5*IT_0275;
    const ccomplex_t IT_0277 = IT_0269*IT_0276;
    const ccomplex_t IT_0278 = IT_0059*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0281 = IT_0003*IT_0280;
    const ccomplex_t IT_0282 = 1.4142135623731*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0284 = IT_0012*IT_0283;
    const ccomplex_t IT_0285 = 1.4142135623731*IT_0284;
    const ccomplex_t IT_0286 = conj(N_d1)*e_em*m_tau*IT_0039*U_stau_11;
    const ccomplex_t IT_0287 = IT_0038*IT_0286;
    const ccomplex_t IT_0288 = 1.4142135623731*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*(IT_0282 + IT_0285 + 
      -IT_0288);
    const ccomplex_t IT_0290 = (-0.5)*IT_0289;
    const ccomplex_t IT_0291 = 2*U_stau_00*((IT_0007*IT_0012*IT_0037*IT_0039
      *IT_0045 + (-0.5)*m_W*IT_0064*(IT_0012 + -IT_0011*IT_0065))*conj(U_stau_01
      ) + 1./2*IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_11)) + U_stau_10*
      (IT_0001*IT_0012*IT_0037*IT_0060*conj(U_stau_01) + 2*(IT_0007*IT_0012
      *IT_0037*IT_0039*IT_0045 + -m_W*IT_0011*IT_0064*IT_0065)*conj(U_stau_11));
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*(e_em*IT_0291 +
       1.4142135623731*A_tau*IT_0007*(conj(U_stau_01)*U_stau_10 + U_stau_00*conj
      (U_stau_11)));
    const ccomplex_t IT_0293 = 0.5*IT_0292;
    const ccomplex_t IT_0294 = IT_0290*IT_0293;
    const ccomplex_t IT_0295 = IT_0072*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = IT_0047*IT_0269;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = IT_0032*IT_0298;
    const ccomplex_t IT_0300 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0301 = IT_0003*IT_0300;
    const ccomplex_t IT_0302 = 1.4142135623731*IT_0301;
    const ccomplex_t IT_0303 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0304 = IT_0012*IT_0303;
    const ccomplex_t IT_0305 = 1.4142135623731*IT_0304;
    const ccomplex_t IT_0306 = conj(N_d2)*e_em*m_tau*IT_0039*U_stau_10;
    const ccomplex_t IT_0307 = IT_0038*IT_0306;
    const ccomplex_t IT_0308 = 1.4142135623731*IT_0307;
    const ccomplex_t IT_0309 = (0 + _Complex_I*1)*(IT_0302 + IT_0305 + 
      -IT_0308);
    const ccomplex_t IT_0310 = (-0.5)*IT_0309;
    const ccomplex_t IT_0311 = IT_0095*IT_0310;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*IT_0311;
    const ccomplex_t IT_0313 = m_N_1*IT_0118;
    const ccomplex_t IT_0314 = m_N_2*IT_0140;
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0312*IT_0315;
    const ccomplex_t IT_0317 = m_N_1*IT_0175;
    const ccomplex_t IT_0318 = m_N_3*IT_0197;
    const ccomplex_t IT_0319 = IT_0317 + IT_0318;
    const ccomplex_t IT_0320 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0321 = IT_0003*IT_0320;
    const ccomplex_t IT_0322 = 1.4142135623731*IT_0321;
    const ccomplex_t IT_0323 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0324 = IT_0012*IT_0323;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = conj(N_d3)*e_em*m_tau*IT_0039*U_stau_10;
    const ccomplex_t IT_0327 = IT_0038*IT_0326;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*(IT_0322 + IT_0325 + 
      -IT_0328);
    const ccomplex_t IT_0330 = (-0.5)*IT_0329;
    const ccomplex_t IT_0331 = IT_0152*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = IT_0319*IT_0332;
    const ccomplex_t IT_0334 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0335 = IT_0003*IT_0334;
    const ccomplex_t IT_0336 = 1.4142135623731*IT_0335;
    const ccomplex_t IT_0337 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0338 = IT_0012*IT_0337;
    const ccomplex_t IT_0339 = 1.4142135623731*IT_0338;
    const ccomplex_t IT_0340 = conj(N_d4)*e_em*m_tau*IT_0039*U_stau_10;
    const ccomplex_t IT_0341 = IT_0038*IT_0340;
    const ccomplex_t IT_0342 = 1.4142135623731*IT_0341;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*(IT_0336 + IT_0339 + 
      -IT_0342);
    const ccomplex_t IT_0344 = (-0.5)*IT_0343;
    const ccomplex_t IT_0345 = IT_0254*IT_0344;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*IT_0345;
    const ccomplex_t IT_0347 = m_N_4*IT_0243;
    const ccomplex_t IT_0348 = m_N_1*IT_0221;
    const ccomplex_t IT_0349 = IT_0347 + IT_0348;
    const ccomplex_t IT_0350 = IT_0346*IT_0349;
    const ccomplex_t IT_0351 = -IT_0273 + 0.5*IT_0279 + 0.5*IT_0296 + 0.5
      *IT_0299 + 0.5*IT_0316 + 0.5*IT_0333 + 0.5*IT_0350;
    const ccomplex_t IT_0352 = 2*s_13;
    const ccomplex_t IT_0353 = m_N_1*m_tau;
    const ccomplex_t IT_0354 = 2*IT_0258;
    const ccomplex_t IT_0355 = IT_0175*IT_0332;
    const ccomplex_t IT_0356 = IT_0221*IT_0346;
    const ccomplex_t IT_0357 = IT_0019*IT_0298;
    const ccomplex_t IT_0358 = IT_0118*IT_0312;
    const ccomplex_t IT_0359 = 0.5*IT_0057 + (-0.5)*IT_0355 + (-0.5)*IT_0356 +
       (-0.5)*IT_0357 + (-0.5)*IT_0358;
    const ccomplex_t IT_0360 = s_14*m_tau;
    const ccomplex_t IT_0361 = IT_0097*IT_0140;
    const ccomplex_t IT_0362 = IT_0154*IT_0197;
    const ccomplex_t IT_0363 = IT_0030*IT_0049;
    const ccomplex_t IT_0364 = IT_0243*IT_0256;
    const ccomplex_t IT_0365 = (-0.5)*IT_0272 + 0.5*IT_0361 + 0.5*IT_0362 +
       0.5*IT_0363 + 0.5*IT_0364;
    const ccomplex_t IT_0366 = s_34*m_N_1;
    const ccomplex_t IT_0367 = 2*conj(IT_0258);
    const ccomplex_t IT_0368 = 2*conj(IT_0351);
    const ccomplex_t IT_0369 = pow(m_h, 2);
    const ccomplex_t IT_0370 = m_N_1*m_tau*IT_0369;
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

