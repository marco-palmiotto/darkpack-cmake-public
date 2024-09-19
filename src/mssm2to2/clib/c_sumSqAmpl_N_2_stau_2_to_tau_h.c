#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_stau_2_to_tau_h.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_stau_2_to_tau_h(
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0005 + IT_0009 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0016;
    const ccomplex_t IT_0021 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0014 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = m_N_3*IT_0030;
    const ccomplex_t IT_0032 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0033 = IT_0003*IT_0032;
    const ccomplex_t IT_0034 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0035 = IT_0007*IT_0034;
    const ccomplex_t IT_0036 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0037 = IT_0003*IT_0036;
    const ccomplex_t IT_0038 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0039 = IT_0007*IT_0038;
    const ccomplex_t IT_0040 = IT_0033 + IT_0035 + IT_0037 + IT_0039;
    const ccomplex_t IT_0041 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0042 = IT_0017*IT_0041;
    const ccomplex_t IT_0043 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0044 = IT_0020*IT_0043;
    const ccomplex_t IT_0045 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0046 = IT_0017*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = -IT_0042 + -IT_0044 + -IT_0046 + -IT_0048;
    const ccomplex_t IT_0050 = IT_0040 + IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*IT_0050;
    const ccomplex_t IT_0052 = (-0.5)*IT_0051;
    const ccomplex_t IT_0053 = m_N_2*IT_0052;
    const ccomplex_t IT_0054 = IT_0031 + IT_0053;
    const ccomplex_t IT_0055 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0056 = IT_0002*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = cos(beta);
    const ccomplex_t IT_0059 = cpow(IT_0058, -1);
    const ccomplex_t IT_0060 = IT_0016*IT_0059;
    const ccomplex_t IT_0061 = pow(m_W, -1);
    const ccomplex_t IT_0062 = N_d3*e_em*m_tau*IT_0061*U_stau_01;
    const ccomplex_t IT_0063 = IT_0060*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0057 + 0.5*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = pow(m_tau, 2);
    const ccomplex_t IT_0068 = pow(m_stau_2, 2);
    const ccomplex_t IT_0069 = cpow(s_23 + (-0.5)*IT_0067 + (-0.5)*IT_0068 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0070 = IT_0066*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*IT_0070;
    const ccomplex_t IT_0072 = IT_0054*IT_0071;
    const ccomplex_t IT_0073 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0074 = IT_0003*IT_0073;
    const ccomplex_t IT_0075 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0076 = IT_0007*IT_0075;
    const ccomplex_t IT_0077 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0078 = IT_0003*IT_0077;
    const ccomplex_t IT_0079 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0080 = IT_0007*IT_0079;
    const ccomplex_t IT_0081 = IT_0074 + IT_0076 + IT_0078 + IT_0080;
    const ccomplex_t IT_0082 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0083 = IT_0017*IT_0082;
    const ccomplex_t IT_0084 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0085 = IT_0020*IT_0084;
    const ccomplex_t IT_0086 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0087 = IT_0017*IT_0086;
    const ccomplex_t IT_0088 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0089 = IT_0020*IT_0088;
    const ccomplex_t IT_0090 = -IT_0083 + -IT_0085 + -IT_0087 + -IT_0089;
    const ccomplex_t IT_0091 = IT_0081 + IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = m_N_1*IT_0093;
    const ccomplex_t IT_0095 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0096 = IT_0003*IT_0095;
    const ccomplex_t IT_0097 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0098 = IT_0007*IT_0097;
    const ccomplex_t IT_0099 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0100 = IT_0003*IT_0099;
    const ccomplex_t IT_0101 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0102 = IT_0007*IT_0101;
    const ccomplex_t IT_0103 = IT_0096 + IT_0098 + IT_0100 + IT_0102;
    const ccomplex_t IT_0104 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0105 = IT_0017*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0107 = IT_0020*IT_0106;
    const ccomplex_t IT_0108 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0109 = IT_0017*IT_0108;
    const ccomplex_t IT_0110 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0111 = IT_0020*IT_0110;
    const ccomplex_t IT_0112 = -IT_0105 + -IT_0107 + -IT_0109 + -IT_0111;
    const ccomplex_t IT_0113 = IT_0103 + IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0114;
    const ccomplex_t IT_0116 = m_N_2*IT_0115;
    const ccomplex_t IT_0117 = IT_0094 + IT_0116;
    const ccomplex_t IT_0118 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0119 = IT_0002*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = N_d1*e_em*m_tau*IT_0061*U_stau_01;
    const ccomplex_t IT_0122 = IT_0060*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0120 + 0.5*IT_0123);
    const ccomplex_t IT_0125 = -IT_0124;
    const ccomplex_t IT_0126 = cpow(s_23 + (-0.5)*IT_0067 + (-0.5)*IT_0068 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0127 = IT_0125*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = IT_0117*IT_0128;
    const ccomplex_t IT_0130 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0131 = IT_0003*IT_0130;
    const ccomplex_t IT_0132 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0133 = IT_0007*IT_0132;
    const ccomplex_t IT_0134 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0135 = IT_0003*IT_0134;
    const ccomplex_t IT_0136 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0137 = IT_0007*IT_0136;
    const ccomplex_t IT_0138 = IT_0131 + IT_0133 + IT_0135 + IT_0137;
    const ccomplex_t IT_0139 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0140 = IT_0017*IT_0139;
    const ccomplex_t IT_0141 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0142 = IT_0020*IT_0141;
    const ccomplex_t IT_0143 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0144 = IT_0017*IT_0143;
    const ccomplex_t IT_0145 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0146 = IT_0020*IT_0145;
    const ccomplex_t IT_0147 = -IT_0140 + -IT_0142 + -IT_0144 + -IT_0146;
    const ccomplex_t IT_0148 = IT_0138 + IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*IT_0148;
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = m_N_4*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0153 = IT_0003*IT_0152;
    const ccomplex_t IT_0154 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0155 = IT_0007*IT_0154;
    const ccomplex_t IT_0156 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0157 = IT_0003*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0159 = IT_0007*IT_0158;
    const ccomplex_t IT_0160 = IT_0153 + IT_0155 + IT_0157 + IT_0159;
    const ccomplex_t IT_0161 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0162 = IT_0017*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0164 = IT_0020*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0166 = IT_0017*IT_0165;
    const ccomplex_t IT_0167 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0168 = IT_0020*IT_0167;
    const ccomplex_t IT_0169 = -IT_0162 + -IT_0164 + -IT_0166 + -IT_0168;
    const ccomplex_t IT_0170 = IT_0160 + IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = m_N_2*IT_0172;
    const ccomplex_t IT_0174 = IT_0151 + IT_0173;
    const ccomplex_t IT_0175 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0176 = IT_0002*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = N_d4*e_em*m_tau*IT_0061*U_stau_01;
    const ccomplex_t IT_0179 = IT_0060*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*(IT_0177 + 0.5*IT_0180);
    const ccomplex_t IT_0182 = -IT_0181;
    const ccomplex_t IT_0183 = cpow(s_23 + (-0.5)*IT_0067 + (-0.5)*IT_0068 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0184 = IT_0182*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*IT_0184;
    const ccomplex_t IT_0186 = IT_0174*IT_0185;
    const ccomplex_t IT_0187 = pow(m_N_2, 2);
    const ccomplex_t IT_0188 = cpow(s_12 + (-0.5)*IT_0067 + 0.5*IT_0068 + 0.5
      *IT_0187 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0016
      *IT_0059*IT_0061;
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0192 = IT_0002*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = N_d2*e_em*m_tau*IT_0061*U_stau_01;
    const ccomplex_t IT_0195 = IT_0060*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + 0.5*IT_0196);
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0190*IT_0198;
    const ccomplex_t IT_0200 = IT_0188*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = m_tau*IT_0201;
    const ccomplex_t IT_0203 = cpow(s_13 + (-0.5)*IT_0067 + (-0.5)*IT_0187 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0204 = mu_h*m_tau*IT_0061;
    const ccomplex_t IT_0205 = IT_0006*IT_0058;
    const ccomplex_t IT_0206 = sin(beta);
    const ccomplex_t IT_0207 = IT_0000*IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = cpow(IT_0001, -2);
    const ccomplex_t IT_0210 = 2*U_stau_01*((IT_0006*IT_0016*IT_0059*IT_0061
      *IT_0067 + (-0.5)*m_W*IT_0208*(IT_0016 + -IT_0015*IT_0209))*conj(U_stau_01
      ) + 1./2*IT_0000*IT_0016*IT_0059*IT_0204*conj(U_stau_11)) + U_stau_11*
      (IT_0000*IT_0016*IT_0059*IT_0204*conj(U_stau_01) + 2*(IT_0006*IT_0016
      *IT_0059*IT_0061*IT_0067 + -m_W*IT_0015*IT_0208*IT_0209)*conj(U_stau_11));
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(e_em*IT_0210 +
       1.4142135623731*A_tau*IT_0006*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0212 = 0.5*IT_0211;
    const ccomplex_t IT_0213 = IT_0198*IT_0212;
    const ccomplex_t IT_0214 = IT_0203*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = cpow(s_13 + (-0.5)*IT_0067 + (-0.5)*IT_0187 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0217 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0218 = IT_0002*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = N_d2*e_em*m_tau*IT_0061*U_stau_00;
    const ccomplex_t IT_0221 = IT_0060*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*(IT_0219 + 0.5*IT_0222);
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = 2*conj(U_stau_00)*((IT_0006*IT_0016*IT_0059
      *IT_0061*IT_0067 + (-0.5)*m_W*IT_0208*(IT_0016 + -IT_0015*IT_0209))
      *U_stau_01 + 1./2*IT_0000*IT_0016*IT_0059*IT_0204*U_stau_11) + conj
      (U_stau_10)*(IT_0000*IT_0016*IT_0059*IT_0204*U_stau_01 + 2*(IT_0006
      *IT_0016*IT_0059*IT_0061*IT_0067 + -m_W*IT_0015*IT_0208*IT_0209)*U_stau_11);
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(e_em*IT_0225 +
       1.4142135623731*A_tau*IT_0006*(U_stau_01*conj(U_stau_10) + conj(U_stau_00
      )*U_stau_11));
    const ccomplex_t IT_0227 = 0.5*IT_0226;
    const ccomplex_t IT_0228 = IT_0224*IT_0227;
    const ccomplex_t IT_0229 = IT_0216*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0232 = IT_0003*IT_0231;
    const ccomplex_t IT_0233 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0234 = IT_0007*IT_0233;
    const ccomplex_t IT_0235 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0236 = IT_0017*IT_0235;
    const ccomplex_t IT_0237 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0238 = IT_0020*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0232 + IT_0234 + 
      -IT_0236 + -IT_0238);
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = m_N_2*IT_0240;
    const ccomplex_t IT_0242 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0243 = IT_0003*IT_0242;
    const ccomplex_t IT_0244 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0245 = IT_0007*IT_0244;
    const ccomplex_t IT_0246 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0247 = IT_0017*IT_0246;
    const ccomplex_t IT_0248 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0249 = IT_0020*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0243 + IT_0245 + 
      -IT_0247 + -IT_0249);
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = m_N_2*IT_0251;
    const ccomplex_t IT_0253 = IT_0241 + IT_0252;
    const ccomplex_t IT_0254 = cpow(s_23 + (-0.5)*IT_0067 + (-0.5)*IT_0068 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0255 = IT_0198*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = IT_0253*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0072 + (-0.5)*IT_0129 + (-0.5)
      *IT_0186 + IT_0202 + (-0.5)*IT_0215 + (-0.5)*IT_0230 + (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0260 = IT_0002*IT_0259;
    const ccomplex_t IT_0261 = 1.4142135623731*IT_0260;
    const ccomplex_t IT_0262 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0263 = IT_0016*IT_0262;
    const ccomplex_t IT_0264 = 1.4142135623731*IT_0263;
    const ccomplex_t IT_0265 = conj(N_d2)*e_em*m_tau*IT_0061*U_stau_11;
    const ccomplex_t IT_0266 = IT_0060*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0261 + IT_0264 + 
      -IT_0267);
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = 2*U_stau_01*((IT_0006*IT_0016*IT_0059*IT_0061
      *IT_0067 + (-0.5)*m_W*IT_0208*(IT_0016 + -IT_0015*IT_0209))*conj(U_stau_01
      ) + 1./2*IT_0000*IT_0016*IT_0059*IT_0204*conj(U_stau_11)) + U_stau_11*
      (IT_0000*IT_0016*IT_0059*IT_0204*conj(U_stau_01) + 2*(IT_0006*IT_0016
      *IT_0059*IT_0061*IT_0067 + -m_W*IT_0015*IT_0208*IT_0209)*conj(U_stau_11));
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*(e_em*IT_0270 +
       1.4142135623731*A_tau*IT_0006*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0272 = 0.5*IT_0271;
    const ccomplex_t IT_0273 = IT_0269*IT_0272;
    const ccomplex_t IT_0274 = IT_0203*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*IT_0274;
    const ccomplex_t IT_0276 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0277 = IT_0002*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0280 = IT_0016*IT_0279;
    const ccomplex_t IT_0281 = 1.4142135623731*IT_0280;
    const ccomplex_t IT_0282 = conj(N_d2)*e_em*m_tau*IT_0061*U_stau_10;
    const ccomplex_t IT_0283 = IT_0060*IT_0282;
    const ccomplex_t IT_0284 = 1.4142135623731*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*(IT_0278 + IT_0281 + 
      -IT_0284);
    const ccomplex_t IT_0286 = (-0.5)*IT_0285;
    const ccomplex_t IT_0287 = 2*conj(U_stau_00)*((IT_0006*IT_0016*IT_0059
      *IT_0061*IT_0067 + (-0.5)*m_W*IT_0208*(IT_0016 + -IT_0015*IT_0209))
      *U_stau_01 + 1./2*IT_0000*IT_0016*IT_0059*IT_0204*U_stau_11) + conj
      (U_stau_10)*(IT_0000*IT_0016*IT_0059*IT_0204*U_stau_01 + 2*(IT_0006
      *IT_0016*IT_0059*IT_0061*IT_0067 + -m_W*IT_0015*IT_0208*IT_0209)*U_stau_11);
    const ccomplex_t IT_0288 = (0 + _Complex_I*1)*(e_em*IT_0287 +
       1.4142135623731*A_tau*IT_0006*(U_stau_01*conj(U_stau_10) + conj(U_stau_00
      )*U_stau_11));
    const ccomplex_t IT_0289 = 0.5*IT_0288;
    const ccomplex_t IT_0290 = IT_0286*IT_0289;
    const ccomplex_t IT_0291 = IT_0216*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = m_N_2*IT_0030;
    const ccomplex_t IT_0294 = m_N_3*IT_0052;
    const ccomplex_t IT_0295 = IT_0293 + IT_0294;
    const ccomplex_t IT_0296 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0297 = IT_0002*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0300 = IT_0016*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = conj(N_d3)*e_em*m_tau*IT_0061*U_stau_11;
    const ccomplex_t IT_0303 = IT_0060*IT_0302;
    const ccomplex_t IT_0304 = 1.4142135623731*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*(IT_0298 + IT_0301 + 
      -IT_0304);
    const ccomplex_t IT_0306 = (-0.5)*IT_0305;
    const ccomplex_t IT_0307 = IT_0069*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = IT_0295*IT_0308;
    const ccomplex_t IT_0310 = IT_0254*IT_0269;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*IT_0310;
    const ccomplex_t IT_0312 = IT_0253*IT_0311;
    const ccomplex_t IT_0313 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0314 = IT_0002*IT_0313;
    const ccomplex_t IT_0315 = 1.4142135623731*IT_0314;
    const ccomplex_t IT_0316 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0317 = IT_0016*IT_0316;
    const ccomplex_t IT_0318 = 1.4142135623731*IT_0317;
    const ccomplex_t IT_0319 = conj(N_d1)*e_em*m_tau*IT_0061*U_stau_11;
    const ccomplex_t IT_0320 = IT_0060*IT_0319;
    const ccomplex_t IT_0321 = 1.4142135623731*IT_0320;
    const ccomplex_t IT_0322 = (0 + _Complex_I*1)*(IT_0315 + IT_0318 + 
      -IT_0321);
    const ccomplex_t IT_0323 = (-0.5)*IT_0322;
    const ccomplex_t IT_0324 = IT_0126*IT_0323;
    const ccomplex_t IT_0325 = (0 + _Complex_I*1)*IT_0324;
    const ccomplex_t IT_0326 = m_N_1*IT_0115;
    const ccomplex_t IT_0327 = m_N_2*IT_0093;
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = IT_0325*IT_0328;
    const ccomplex_t IT_0330 = m_N_2*IT_0150;
    const ccomplex_t IT_0331 = m_N_4*IT_0172;
    const ccomplex_t IT_0332 = IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0334 = IT_0002*IT_0333;
    const ccomplex_t IT_0335 = 1.4142135623731*IT_0334;
    const ccomplex_t IT_0336 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0337 = IT_0016*IT_0336;
    const ccomplex_t IT_0338 = 1.4142135623731*IT_0337;
    const ccomplex_t IT_0339 = conj(N_d4)*e_em*m_tau*IT_0061*U_stau_11;
    const ccomplex_t IT_0340 = IT_0060*IT_0339;
    const ccomplex_t IT_0341 = 1.4142135623731*IT_0340;
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*(IT_0335 + IT_0338 + 
      -IT_0341);
    const ccomplex_t IT_0343 = (-0.5)*IT_0342;
    const ccomplex_t IT_0344 = IT_0183*IT_0343;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*IT_0344;
    const ccomplex_t IT_0346 = IT_0332*IT_0345;
    const ccomplex_t IT_0347 = IT_0190*IT_0269;
    const ccomplex_t IT_0348 = IT_0188*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = m_tau*IT_0349;
    const ccomplex_t IT_0351 = 0.5*IT_0275 + 0.5*IT_0292 + 0.5*IT_0309 + 0.5
      *IT_0312 + 0.5*IT_0329 + 0.5*IT_0346 + -IT_0350;
    const ccomplex_t IT_0352 = 2*s_13;
    const ccomplex_t IT_0353 = m_N_2*m_tau;
    const ccomplex_t IT_0354 = 2*IT_0258;
    const ccomplex_t IT_0355 = IT_0150*IT_0345;
    const ccomplex_t IT_0356 = IT_0093*IT_0325;
    const ccomplex_t IT_0357 = IT_0030*IT_0308;
    const ccomplex_t IT_0358 = IT_0240*IT_0311;
    const ccomplex_t IT_0359 = 0.5*IT_0201 + (-0.5)*IT_0355 + (-0.5)*IT_0356 +
       (-0.5)*IT_0357 + (-0.5)*IT_0358;
    const ccomplex_t IT_0360 = s_14*m_tau;
    const ccomplex_t IT_0361 = IT_0251*IT_0256;
    const ccomplex_t IT_0362 = IT_0115*IT_0128;
    const ccomplex_t IT_0363 = IT_0052*IT_0071;
    const ccomplex_t IT_0364 = IT_0172*IT_0185;
    const ccomplex_t IT_0365 = (-0.5)*IT_0349 + 0.5*IT_0361 + 0.5*IT_0362 +
       0.5*IT_0363 + 0.5*IT_0364;
    const ccomplex_t IT_0366 = s_34*m_N_2;
    const ccomplex_t IT_0367 = 2*conj(IT_0258);
    const ccomplex_t IT_0368 = 2*conj(IT_0351);
    const ccomplex_t IT_0369 = pow(m_h, 2);
    const ccomplex_t IT_0370 = m_N_2*m_tau*IT_0369;
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

