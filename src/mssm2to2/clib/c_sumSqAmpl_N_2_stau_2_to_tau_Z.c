#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_stau_2_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_stau_2_to_tau_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
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
    const ccomplex_t IT_0000 = s_14*s_34;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(m_Z, 4);
    const ccomplex_t IT_0003 = s_13*IT_0002;
    const ccomplex_t IT_0004 = IT_0001*IT_0003;
    const ccomplex_t IT_0005 = (-0.25)*IT_0004;
    const ccomplex_t IT_0006 = IT_0000 + IT_0005;
    const ccomplex_t IT_0007 = cos(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0009, -1);
    const ccomplex_t IT_0014 = IT_0007*IT_0013;
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0017 = IT_0010*IT_0016;
    const ccomplex_t IT_0018 = IT_0014*IT_0016;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0012 + IT_0015 + 
      -IT_0017 + -IT_0018);
    const ccomplex_t IT_0020 = (-0.5)*IT_0019;
    const ccomplex_t IT_0021 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0022 = IT_0008*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = cos(beta);
    const ccomplex_t IT_0025 = cpow(IT_0024, -1);
    const ccomplex_t IT_0026 = IT_0013*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = N_d2*e_em*m_tau*IT_0027*U_stau_01;
    const ccomplex_t IT_0029 = IT_0026*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0030);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = pow(m_tau, 2);
    const ccomplex_t IT_0034 = pow(m_stau_2, 2);
    const ccomplex_t IT_0035 = cpow(s_23 + (-0.5)*IT_0033 + (-0.5)*IT_0034 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0036 = IT_0032*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = IT_0020*IT_0037;
    const ccomplex_t IT_0039 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0040 = IT_0010*IT_0039;
    const ccomplex_t IT_0041 = IT_0014*IT_0039;
    const ccomplex_t IT_0042 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0043 = IT_0010*IT_0042;
    const ccomplex_t IT_0044 = IT_0014*IT_0042;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*(IT_0040 + IT_0041 + 
      -IT_0043 + -IT_0044);
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0048 = IT_0008*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = N_d3*e_em*m_tau*IT_0027*U_stau_01;
    const ccomplex_t IT_0051 = IT_0026*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + 0.5*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = cpow(s_23 + (-0.5)*IT_0033 + (-0.5)*IT_0034 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0046*IT_0057;
    const ccomplex_t IT_0059 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0060 = IT_0010*IT_0059;
    const ccomplex_t IT_0061 = IT_0014*IT_0059;
    const ccomplex_t IT_0062 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0063 = IT_0010*IT_0062;
    const ccomplex_t IT_0064 = IT_0014*IT_0062;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0060 + IT_0061 + 
      -IT_0063 + -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0068 = IT_0008*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = N_d4*e_em*m_tau*IT_0027*U_stau_01;
    const ccomplex_t IT_0071 = IT_0026*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 0.5*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = cpow(s_23 + (-0.5)*IT_0033 + (-0.5)*IT_0034 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0066*IT_0077;
    const ccomplex_t IT_0079 = 0.25*IT_0038 + 0.5*IT_0058 + 0.5*IT_0078;
    const ccomplex_t IT_0080 = 0.5*IT_0019;
    const ccomplex_t IT_0081 = IT_0037*IT_0080;
    const ccomplex_t IT_0082 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0083 = IT_0010*IT_0082;
    const ccomplex_t IT_0084 = IT_0014*IT_0082;
    const ccomplex_t IT_0085 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0086 = IT_0010*IT_0085;
    const ccomplex_t IT_0087 = IT_0014*IT_0085;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*(IT_0083 + IT_0084 + 
      -IT_0086 + -IT_0087);
    const ccomplex_t IT_0089 = (-0.5)*IT_0088;
    const ccomplex_t IT_0090 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0091 = IT_0008*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = N_d1*e_em*m_tau*IT_0027*U_stau_01;
    const ccomplex_t IT_0094 = IT_0026*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0092 + 0.5*IT_0095);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = cpow(s_23 + (-0.5)*IT_0033 + (-0.5)*IT_0034 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0099 = IT_0097*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*IT_0099;
    const ccomplex_t IT_0101 = IT_0089*IT_0100;
    const ccomplex_t IT_0102 = 0.25*IT_0081 + 0.5*IT_0101;
    const ccomplex_t IT_0103 = pow(m_N_2, 2);
    const ccomplex_t IT_0104 = cpow(s_12 + (-0.5)*IT_0033 + 0.5*IT_0034 + 0.5
      *IT_0103 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0009;
    const ccomplex_t IT_0106 = IT_0032*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0108;
    const ccomplex_t IT_0110 = -conj(IT_0102);
    const ccomplex_t IT_0111 = conj(IT_0079) + IT_0110;
    const ccomplex_t IT_0112 = e_em*IT_0010;
    const ccomplex_t IT_0113 = e_em*IT_0014;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*(IT_0112 + -IT_0113);
    const ccomplex_t IT_0115 = 0.5*IT_0114;
    const ccomplex_t IT_0116 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0117 = IT_0008*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0120 = IT_0013*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = conj(N_d2)*e_em*m_tau*IT_0027*U_stau_11;
    const ccomplex_t IT_0123 = IT_0026*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0118 + IT_0121 + 
      -IT_0124);
    const ccomplex_t IT_0126 = (-0.5)*IT_0125;
    const ccomplex_t IT_0127 = IT_0115*IT_0126;
    const ccomplex_t IT_0128 = IT_0104*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = IT_0035*IT_0126;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = IT_0020*IT_0132;
    const ccomplex_t IT_0134 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0135 = IT_0010*IT_0134;
    const ccomplex_t IT_0136 = IT_0014*IT_0134;
    const ccomplex_t IT_0137 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0138 = IT_0010*IT_0137;
    const ccomplex_t IT_0139 = IT_0014*IT_0137;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0135 + IT_0136 + 
      -IT_0138 + -IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0143 = IT_0008*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0146 = IT_0013*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = conj(N_d1)*e_em*m_tau*IT_0027*U_stau_11;
    const ccomplex_t IT_0149 = IT_0026*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0144 + IT_0147 + 
      -IT_0150);
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0098*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0141*IT_0154;
    const ccomplex_t IT_0156 = (-0.25)*IT_0133 + (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0158 = IT_0010*IT_0157;
    const ccomplex_t IT_0159 = IT_0014*IT_0157;
    const ccomplex_t IT_0160 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0161 = IT_0010*IT_0160;
    const ccomplex_t IT_0162 = IT_0014*IT_0160;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0158 + IT_0159 + 
      -IT_0161 + -IT_0162);
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0166 = IT_0008*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0169 = IT_0013*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_d4)*e_em*m_tau*IT_0027*U_stau_11;
    const ccomplex_t IT_0172 = IT_0026*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + IT_0170 + 
      -IT_0173);
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = IT_0075*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0164*IT_0177;
    const ccomplex_t IT_0179 = IT_0080*IT_0132;
    const ccomplex_t IT_0180 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0181 = IT_0010*IT_0180;
    const ccomplex_t IT_0182 = IT_0014*IT_0180;
    const ccomplex_t IT_0183 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0184 = IT_0010*IT_0183;
    const ccomplex_t IT_0185 = IT_0014*IT_0183;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0181 + IT_0182 + 
      -IT_0184 + -IT_0185);
    const ccomplex_t IT_0187 = (-0.5)*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0189 = IT_0008*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0192 = IT_0013*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = conj(N_d3)*e_em*m_tau*IT_0027*U_stau_11;
    const ccomplex_t IT_0195 = IT_0026*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0190 + IT_0193 + 
      -IT_0196);
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = IT_0055*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = IT_0187*IT_0200;
    const ccomplex_t IT_0202 = (-0.5)*IT_0178 + (-0.25)*IT_0179 + (-0.5)
      *IT_0201;
    const ccomplex_t IT_0203 = IT_0156 + -IT_0202;
    const ccomplex_t IT_0204 = conj(IT_0156) + -conj(IT_0202);
    const ccomplex_t IT_0205 = pow(m_Z, 2);
    const ccomplex_t IT_0206 = s_13*IT_0205;
    const ccomplex_t IT_0207 = (-0.5)*IT_0004;
    const ccomplex_t IT_0208 = (-2)*IT_0000;
    const ccomplex_t IT_0209 = IT_0206 + IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = m_N_2*m_tau;
    const ccomplex_t IT_0211 = pow(s_14, 2);
    const ccomplex_t IT_0212 = IT_0001*IT_0211;
    const ccomplex_t IT_0213 = -IT_0103;
    const ccomplex_t IT_0214 = IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = IT_0210*IT_0214;
    const ccomplex_t IT_0216 = 0.5*IT_0133 + IT_0155;
    const ccomplex_t IT_0217 = (-0.5)*IT_0081 + -IT_0101;
    const ccomplex_t IT_0218 = 2*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0008*IT_0009 +
       -IT_0007*IT_0013)*U_stau_01*conj(U_stau_01) + IT_0008*IT_0009*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0220 = (-2)*IT_0219;
    const ccomplex_t IT_0221 = cpow(s_13 + (-0.5)*IT_0033 + (-0.5)*IT_0103 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0222 = IT_0032*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = IT_0220*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0008*IT_0009 +
       -IT_0007*IT_0013)*conj(U_stau_00)*U_stau_01 + IT_0008*IT_0009*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0226 = (-2)*IT_0225;
    const ccomplex_t IT_0227 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0228 = IT_0008*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = N_d2*e_em*m_tau*IT_0027*U_stau_00;
    const ccomplex_t IT_0231 = IT_0026*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0229 + 0.5*IT_0232);
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = cpow(s_13 + (-0.5)*IT_0033 + (-0.5)*IT_0103 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0236 = IT_0234*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = IT_0226*IT_0237;
    const ccomplex_t IT_0239 = (-0.5)*IT_0038 + -IT_0058 + -IT_0078 + (-0.5)
      *IT_0224 + (-0.5)*IT_0238;
    const ccomplex_t IT_0240 = (-2)*IT_0239;
    const ccomplex_t IT_0241 = IT_0218 + IT_0240;
    const ccomplex_t IT_0242 = 2*conj(IT_0217);
    const ccomplex_t IT_0243 = (-2)*conj(IT_0239);
    const ccomplex_t IT_0244 = IT_0242 + IT_0243;
    const ccomplex_t IT_0245 = IT_0126*IT_0221;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = IT_0220*IT_0246;
    const ccomplex_t IT_0248 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0249 = IT_0008*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0252 = IT_0013*IT_0251;
    const ccomplex_t IT_0253 = 1.4142135623731*IT_0252;
    const ccomplex_t IT_0254 = conj(N_d2)*e_em*m_tau*IT_0027*U_stau_10;
    const ccomplex_t IT_0255 = IT_0026*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*(IT_0250 + IT_0253 + 
      -IT_0256);
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = IT_0235*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*IT_0259;
    const ccomplex_t IT_0261 = IT_0226*IT_0260;
    const ccomplex_t IT_0262 = IT_0178 + 0.5*IT_0179 + IT_0201 + 0.5*IT_0247 +
       0.5*IT_0261;
    const ccomplex_t IT_0263 = 2*IT_0239;
    const ccomplex_t IT_0264 = (-2)*IT_0217;
    const ccomplex_t IT_0265 = IT_0263 + IT_0264;
    const ccomplex_t IT_0266 = 2*conj(IT_0239);
    const ccomplex_t IT_0267 = (-2)*conj(IT_0217);
    const ccomplex_t IT_0268 = IT_0266 + IT_0267;
    const ccomplex_t IT_0269 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = s_13 + IT_0270;
    const ccomplex_t IT_0272 = s_13*IT_0271;
    const ccomplex_t IT_0273 = IT_0108 + 0.5*IT_0224 + 0.5*IT_0238;
    const ccomplex_t IT_0274 = -IT_0129 + (-0.5)*IT_0247 + (-0.5)*IT_0261;
    const ccomplex_t IT_0275 = 2*IT_0216;
    const ccomplex_t IT_0276 = (-2)*IT_0262;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = 2*conj(IT_0216);
    const ccomplex_t IT_0279 = (-2)*conj(IT_0262);
    const ccomplex_t IT_0280 = IT_0278 + IT_0279;
    const ccomplex_t IT_0281 = m_N_2*IT_0271;
    const ccomplex_t IT_0282 = m_N_2*IT_0020;
    const ccomplex_t IT_0283 = m_N_2*IT_0080;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = IT_0282 + IT_0284;
    const ccomplex_t IT_0286 = IT_0037*IT_0285;
    const ccomplex_t IT_0287 = m_N_1*IT_0141;
    const ccomplex_t IT_0288 = m_N_2*IT_0089;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = IT_0100*IT_0290;
    const ccomplex_t IT_0292 = 0.25*IT_0286 + 0.5*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0294 = IT_0105*IT_0126;
    const ccomplex_t IT_0295 = IT_0104*IT_0293*IT_0294;
    const ccomplex_t IT_0296 = m_tau*IT_0129;
    const ccomplex_t IT_0297 = m_N_3*IT_0187;
    const ccomplex_t IT_0298 = m_N_2*IT_0046;
    const ccomplex_t IT_0299 = -IT_0298;
    const ccomplex_t IT_0300 = IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = IT_0057*IT_0300;
    const ccomplex_t IT_0302 = m_N_4*IT_0164;
    const ccomplex_t IT_0303 = m_N_2*IT_0066;
    const ccomplex_t IT_0304 = -IT_0303;
    const ccomplex_t IT_0305 = IT_0302 + IT_0304;
    const ccomplex_t IT_0306 = IT_0077*IT_0305;
    const ccomplex_t IT_0307 = 0.25*IT_0286 + (-0.5)*IT_0295 + 0.5*IT_0296 + (
      -0.5)*IT_0301 + (-0.5)*IT_0306;
    const ccomplex_t IT_0308 = conj(IT_0292) + conj(IT_0307);
    const ccomplex_t IT_0309 = IT_0292 + IT_0307;
    const ccomplex_t IT_0310 = IT_0132*IT_0285;
    const ccomplex_t IT_0311 = m_N_1*IT_0089;
    const ccomplex_t IT_0312 = m_N_2*IT_0141;
    const ccomplex_t IT_0313 = -IT_0312;
    const ccomplex_t IT_0314 = IT_0311 + IT_0313;
    const ccomplex_t IT_0315 = IT_0154*IT_0314;
    const ccomplex_t IT_0316 = 0.25*IT_0310 + (-0.5)*IT_0315;
    const ccomplex_t IT_0317 = IT_0032*IT_0115;
    const ccomplex_t IT_0318 = IT_0104*IT_0293*IT_0317;
    const ccomplex_t IT_0319 = m_N_3*IT_0046;
    const ccomplex_t IT_0320 = m_N_2*IT_0187;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0319 + IT_0321;
    const ccomplex_t IT_0323 = IT_0200*IT_0322;
    const ccomplex_t IT_0324 = m_N_4*IT_0066;
    const ccomplex_t IT_0325 = m_N_2*IT_0164;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = IT_0324 + IT_0326;
    const ccomplex_t IT_0328 = IT_0177*IT_0327;
    const ccomplex_t IT_0329 = m_tau*IT_0108;
    const ccomplex_t IT_0330 = 0.25*IT_0310 + 0.5*IT_0318 + 0.5*IT_0323 + 0.5
      *IT_0328 + (-0.5)*IT_0329;
    const ccomplex_t IT_0331 = conj(IT_0316) + conj(IT_0330);
    const ccomplex_t IT_0332 = IT_0316 + IT_0330;
    const ccomplex_t IT_0333 = s_13*s_14;
    const ccomplex_t IT_0334 = s_34*IT_0103;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = conj(IT_0102) + conj(IT_0109);
    const ccomplex_t IT_0338 = IT_0102 + IT_0109;
    const ccomplex_t IT_0339 = conj(IT_0130) + conj(IT_0156);
    const ccomplex_t IT_0340 = IT_0130 + IT_0156;
    const ccomplex_t IT_0341 = 2*IT_0262;
    const ccomplex_t IT_0342 = (-2)*IT_0216;
    const ccomplex_t IT_0343 = IT_0341 + IT_0342;
    const ccomplex_t IT_0344 = 2*conj(IT_0262);
    const ccomplex_t IT_0345 = (-2)*conj(IT_0216);
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = s_14*IT_0205;
    const ccomplex_t IT_0348 = IT_0001*IT_0347;
    const ccomplex_t IT_0349 = (-0.25)*IT_0348;
    const ccomplex_t IT_0350 = s_14 + IT_0349;
    const ccomplex_t IT_0351 = m_tau*IT_0350;
    const ccomplex_t IT_0352 = (-4)*IT_0351;
    const ccomplex_t IT_0353 = 2*IT_0292;
    const ccomplex_t IT_0354 = 2*IT_0307;
    const ccomplex_t IT_0355 = IT_0353 + IT_0354;
    const ccomplex_t IT_0356 = 2*conj(IT_0292);
    const ccomplex_t IT_0357 = 2*conj(IT_0307);
    const ccomplex_t IT_0358 = IT_0356 + IT_0357;
    const ccomplex_t IT_0359 = 2*IT_0316;
    const ccomplex_t IT_0360 = 2*IT_0330;
    const ccomplex_t IT_0361 = IT_0359 + IT_0360;
    const ccomplex_t IT_0362 = 2*conj(IT_0316);
    const ccomplex_t IT_0363 = 2*conj(IT_0330);
    const ccomplex_t IT_0364 = IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = s_34*IT_0205;
    const ccomplex_t IT_0366 = IT_0001*IT_0365;
    const ccomplex_t IT_0367 = 0.5*IT_0366;
    const ccomplex_t IT_0368 = s_34 + IT_0367;
    const ccomplex_t IT_0369 = m_N_2*IT_0368;
    const ccomplex_t IT_0370 = 2*IT_0369;
    const ccomplex_t IT_0371 = 0.5*IT_0348;
    const ccomplex_t IT_0372 = s_14 + IT_0371;
    const ccomplex_t IT_0373 = m_tau*IT_0372;
    const ccomplex_t IT_0374 = 2*IT_0373;
    const ccomplex_t IT_0375 = (-2)*IT_0292;
    const ccomplex_t IT_0376 = (-2)*IT_0307;
    const ccomplex_t IT_0377 = IT_0375 + IT_0376;
    const ccomplex_t IT_0378 = (-2)*conj(IT_0292);
    const ccomplex_t IT_0379 = (-2)*conj(IT_0307);
    const ccomplex_t IT_0380 = IT_0378 + IT_0379;
    const ccomplex_t IT_0381 = (-2)*IT_0316;
    const ccomplex_t IT_0382 = (-2)*IT_0330;
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = (-2)*conj(IT_0316);
    const ccomplex_t IT_0385 = (-2)*conj(IT_0330);
    const ccomplex_t IT_0386 = IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = (-0.25)*IT_0366;
    const ccomplex_t IT_0388 = s_34 + IT_0387;
    const ccomplex_t IT_0389 = m_N_2*IT_0388;
    const ccomplex_t IT_0390 = (-4)*IT_0389;
    const ccomplex_t IT_0391 = pow(s_34, 2);
    const ccomplex_t IT_0392 = IT_0001*IT_0391;
    const ccomplex_t IT_0393 = -IT_0033;
    const ccomplex_t IT_0394 = IT_0392 + IT_0393;
    const ccomplex_t IT_0395 = IT_0210*IT_0394;
    const ccomplex_t IT_0396 = 2*IT_0273;
    const ccomplex_t IT_0397 = 2*conj(IT_0273);
    const ccomplex_t IT_0398 = s_13*IT_0394;
    const ccomplex_t IT_0399 = 2*IT_0274;
    const ccomplex_t IT_0400 = m_N_2*IT_0394;
    const ccomplex_t IT_0401 = 2*conj(IT_0274);
    const ccomplex_t IT_0402 = IT_0206 + IT_0208;
    const ccomplex_t IT_0403 = IT_0001*IT_0402;
    const ccomplex_t IT_0404 = (-0.5)*IT_0403;
    const ccomplex_t IT_0405 = s_13 + IT_0404;
    const ccomplex_t IT_0406 = 4*IT_0316;
    const ccomplex_t IT_0407 = (-6)*conj(IT_0292);
    const ccomplex_t IT_0408 = s_13*IT_0214;
    const ccomplex_t IT_0409 = IT_0210*IT_0271;
    const ccomplex_t IT_0410 = (-2)*IT_0273;
    const ccomplex_t IT_0411 = (-2)*conj(IT_0273);
    const ccomplex_t IT_0412 = s_13*s_34;
    const ccomplex_t IT_0413 = s_14*IT_0033;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0412 + IT_0414;
    const ccomplex_t IT_0416 = (-2)*IT_0274;
    const ccomplex_t IT_0417 = (-2)*conj(IT_0274);
    const ccomplex_t IT_0418 = m_tau*IT_0214;
    const ccomplex_t IT_0419 = m_tau*IT_0271;
    const ccomplex_t IT_0420 = IT_0205*IT_0210;
    const ccomplex_t IT_0421 = (-3)*IT_0420;
    const ccomplex_t IT_0422 = 2*IT_0079;
    const ccomplex_t IT_0423 = (-2)*IT_0102;
    const ccomplex_t IT_0424 = IT_0422 + IT_0423;
    const ccomplex_t IT_0425 = 2*conj(IT_0079);
    const ccomplex_t IT_0426 = (-2)*conj(IT_0102);
    const ccomplex_t IT_0427 = IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = 2*IT_0109;
    const ccomplex_t IT_0429 = 2*conj(IT_0109);
    const ccomplex_t IT_0430 = (-2)*IT_0156;
    const ccomplex_t IT_0431 = (-2)*conj(IT_0156);
    const ccomplex_t IT_0432 = conj(IT_0102)*IT_0156 + IT_0102*conj(IT_0156) +
       0.5*conj(IT_0202)*IT_0424 + 0.5*IT_0202*IT_0427 + 0.5*conj(IT_0130)
      *IT_0428 + 0.5*IT_0130*IT_0429 + 0.5*conj(IT_0079)*IT_0430 + 0.5*IT_0079
      *IT_0431;
    const ccomplex_t IT_0433 = 3*IT_0420;
    const ccomplex_t IT_0434 = (-8)*IT_0006*((IT_0079 + -IT_0102)*conj(IT_0109
      ) + IT_0109*IT_0111 + -conj(IT_0130)*IT_0203 + -IT_0130*IT_0204) + (-4)*
      (IT_0109*conj(IT_0109) + (IT_0079 + -IT_0102)*IT_0111 + IT_0130*conj
      (IT_0130) + IT_0203*IT_0204)*IT_0209 + IT_0215*(conj(IT_0216)*IT_0241 +
       IT_0216*IT_0244 + conj(IT_0262)*IT_0265 + IT_0262*IT_0268) + IT_0272*
      (IT_0244*IT_0273 + IT_0241*conj(IT_0273) + conj(IT_0274)*IT_0277 + IT_0274
      *IT_0280) + IT_0281*(IT_0241*IT_0308 + IT_0244*IT_0309 + IT_0277*IT_0331 +
       IT_0280*IT_0332) + IT_0336*(conj(IT_0079)*IT_0265 + IT_0079*IT_0268 +
       IT_0241*IT_0337 + IT_0244*IT_0338 + IT_0277*IT_0339 + IT_0280*IT_0340 +
       conj(IT_0202)*IT_0343 + IT_0202*IT_0346) + IT_0352*(conj(IT_0130)*IT_0355
       + IT_0130*IT_0358 + conj(IT_0109)*IT_0361 + IT_0109*IT_0364) + (conj
      (IT_0109)*IT_0355 + IT_0109*IT_0358 + conj(IT_0130)*IT_0361 + IT_0130
      *IT_0364)*IT_0370 + IT_0374*(conj(IT_0202)*IT_0355 + IT_0202*IT_0358 +
       conj(IT_0079)*IT_0361 + IT_0079*IT_0364 + conj(IT_0156)*IT_0377 + IT_0156
      *IT_0380 + conj(IT_0102)*IT_0383 + IT_0102*IT_0386) + (conj(IT_0079)
      *IT_0355 + IT_0079*IT_0358 + conj(IT_0202)*IT_0361 + IT_0202*IT_0364 +
       conj(IT_0102)*IT_0377 + IT_0102*IT_0380 + conj(IT_0156)*IT_0383 + IT_0156
      *IT_0386)*IT_0390 + IT_0395*(conj(IT_0274)*IT_0396 + IT_0274*IT_0397) +
       IT_0398*(conj(IT_0273)*IT_0396 + conj(IT_0274)*IT_0399) + IT_0400*
      (IT_0308*IT_0396 + IT_0309*IT_0397 + IT_0331*IT_0399 + IT_0332*IT_0401) +
       4*IT_0405*(IT_0308*IT_0309 + IT_0331*(IT_0330 + 0.25*IT_0406)) + (-6)
      *IT_0210*(IT_0309*IT_0331 + IT_0332*(conj(IT_0307) + (-0.166666666666667)
      *IT_0407)) + (conj(IT_0217)*IT_0218 + conj(IT_0239)*IT_0265 + IT_0239
      *IT_0267 + conj(IT_0216)*IT_0275 + conj(IT_0262)*IT_0343 + IT_0262*IT_0345
      )*IT_0408 + IT_0409*(IT_0244*IT_0274 + IT_0241*conj(IT_0274) + conj
      (IT_0216)*IT_0396 + IT_0216*IT_0397 + conj(IT_0262)*IT_0410 + IT_0262
      *IT_0411) + IT_0415*(IT_0337*IT_0396 + IT_0338*IT_0397 + IT_0339*IT_0399 +
       IT_0340*IT_0401 + conj(IT_0079)*IT_0410 + IT_0079*IT_0411 + conj(IT_0202)
      *IT_0416 + IT_0202*IT_0417) + (IT_0265*IT_0331 + IT_0268*IT_0332 + conj
      (IT_0307)*IT_0343 + IT_0307*IT_0346 + conj(IT_0262)*IT_0353 + IT_0262
      *IT_0356 + conj(IT_0216)*IT_0375 + IT_0216*IT_0378)*IT_0418 + (conj
      (IT_0274)*IT_0375 + IT_0274*IT_0378 + IT_0331*IT_0410 + IT_0332*IT_0411 +
       conj(IT_0307)*IT_0416 + IT_0307*IT_0417)*IT_0419 + 2*IT_0421*IT_0432 + 
      (conj(IT_0130)*IT_0424 + IT_0130*IT_0427 + conj(IT_0202)*IT_0428 + IT_0202
      *IT_0429 + conj(IT_0109)*IT_0430 + IT_0109*IT_0431)*IT_0433;
    return create_ccomplex_return(IT_0434);
}

