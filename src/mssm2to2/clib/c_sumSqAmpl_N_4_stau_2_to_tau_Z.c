#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_stau_2_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_stau_2_to_tau_Z(
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
    const ccomplex_t IT_0000 = m_N_4*m_tau;
    const ccomplex_t IT_0001 = pow(m_Z, -2);
    const ccomplex_t IT_0002 = pow(s_14, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = pow(m_N_4, 2);
    const ccomplex_t IT_0005 = -IT_0004;
    const ccomplex_t IT_0006 = IT_0003 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0013 = IT_0011*IT_0012;
    const ccomplex_t IT_0014 = cpow(IT_0010, -1);
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0017;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0013 + IT_0016 + 
      -IT_0018 + -IT_0019);
    const ccomplex_t IT_0021 = -IT_0020;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0024 = IT_0009*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0027 = IT_0014*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = cos(beta);
    const ccomplex_t IT_0030 = cpow(IT_0029, -1);
    const ccomplex_t IT_0031 = IT_0014*IT_0030;
    const ccomplex_t IT_0032 = pow(m_W, -1);
    const ccomplex_t IT_0033 = conj(N_d2)*e_em*m_tau*IT_0032*U_stau_11;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0025 + IT_0028 + 
      -IT_0035);
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = pow(m_tau, 2);
    const ccomplex_t IT_0039 = pow(m_stau_2, 2);
    const ccomplex_t IT_0040 = cpow(s_23 + (-0.5)*IT_0038 + (-0.5)*IT_0039 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0022*IT_0042;
    const ccomplex_t IT_0044 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0045 = IT_0011*IT_0044;
    const ccomplex_t IT_0046 = IT_0015*IT_0044;
    const ccomplex_t IT_0047 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0048 = IT_0011*IT_0047;
    const ccomplex_t IT_0049 = IT_0015*IT_0047;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0045 + IT_0046 + 
      -IT_0048 + -IT_0049);
    const ccomplex_t IT_0051 = (-0.5)*IT_0050;
    const ccomplex_t IT_0052 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0053 = IT_0009*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0056 = IT_0014*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = conj(N_d3)*e_em*m_tau*IT_0032*U_stau_11;
    const ccomplex_t IT_0059 = IT_0031*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0054 + IT_0057 + 
      -IT_0060);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = cpow(s_23 + (-0.5)*IT_0038 + (-0.5)*IT_0039 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0051*IT_0065;
    const ccomplex_t IT_0067 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0068 = IT_0015*IT_0067;
    const ccomplex_t IT_0069 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0070 = IT_0011*IT_0069;
    const ccomplex_t IT_0071 = IT_0015*IT_0069;
    const ccomplex_t IT_0072 = IT_0011*IT_0067;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + -IT_0070 + 
      -IT_0071 + IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = (-0.5)*IT_0074;
    const ccomplex_t IT_0076 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0077 = IT_0009*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0080 = IT_0014*IT_0079;
    const ccomplex_t IT_0081 = 1.4142135623731*IT_0080;
    const ccomplex_t IT_0082 = conj(N_d4)*e_em*m_tau*IT_0032*U_stau_11;
    const ccomplex_t IT_0083 = IT_0031*IT_0082;
    const ccomplex_t IT_0084 = 1.4142135623731*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*(IT_0078 + IT_0081 + 
      -IT_0084);
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = cpow(s_23 + (-0.5)*IT_0038 + (-0.5)*IT_0039 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0088 = IT_0086*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = IT_0075*IT_0089;
    const ccomplex_t IT_0091 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0092 = IT_0011*IT_0091;
    const ccomplex_t IT_0093 = IT_0015*IT_0091;
    const ccomplex_t IT_0094 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0095 = IT_0011*IT_0094;
    const ccomplex_t IT_0096 = IT_0015*IT_0094;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0092 + IT_0093 + 
      -IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0100 = IT_0009*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0103 = IT_0014*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = conj(N_d1)*e_em*m_tau*IT_0032*U_stau_11;
    const ccomplex_t IT_0106 = IT_0031*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0101 + IT_0104 + 
      -IT_0107);
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = cpow(s_23 + (-0.5)*IT_0038 + (-0.5)*IT_0039 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0098*IT_0112;
    const ccomplex_t IT_0114 = IT_0043 + IT_0066 + 0.5*IT_0090 + IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0074;
    const ccomplex_t IT_0116 = N_d4*e_em*m_tau*IT_0032*U_stau_01;
    const ccomplex_t IT_0117 = IT_0031*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0120 = IT_0009*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(0.5*IT_0118 + IT_0121);
    const ccomplex_t IT_0123 = -IT_0122;
    const ccomplex_t IT_0124 = IT_0087*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = IT_0115*IT_0125;
    const ccomplex_t IT_0127 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0128 = IT_0011*IT_0127;
    const ccomplex_t IT_0129 = IT_0015*IT_0127;
    const ccomplex_t IT_0130 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0131 = IT_0015*IT_0130;
    const ccomplex_t IT_0132 = IT_0011*IT_0130;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0128 + IT_0129 + 
      -IT_0131 + -IT_0132);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = (-0.5)*IT_0134;
    const ccomplex_t IT_0136 = N_d1*e_em*m_tau*IT_0032*U_stau_01;
    const ccomplex_t IT_0137 = IT_0031*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0140 = IT_0009*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 2*IT_0141);
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = -IT_0143;
    const ccomplex_t IT_0145 = IT_0110*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = IT_0135*IT_0146;
    const ccomplex_t IT_0148 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0149 = IT_0011*IT_0148;
    const ccomplex_t IT_0150 = IT_0015*IT_0148;
    const ccomplex_t IT_0151 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0152 = IT_0015*IT_0151;
    const ccomplex_t IT_0153 = IT_0011*IT_0151;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0149 + IT_0150 + 
      -IT_0152 + -IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = (-0.5)*IT_0157;
    const ccomplex_t IT_0159 = N_d2*e_em*m_tau*IT_0032*U_stau_01;
    const ccomplex_t IT_0160 = IT_0031*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0163 = IT_0009*IT_0162;
    const ccomplex_t IT_0164 = 1.4142135623731*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0161 + 2*IT_0164);
    const ccomplex_t IT_0166 = 0.5*IT_0165;
    const ccomplex_t IT_0167 = -IT_0166;
    const ccomplex_t IT_0168 = IT_0040*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*IT_0168;
    const ccomplex_t IT_0170 = IT_0158*IT_0169;
    const ccomplex_t IT_0171 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0172 = IT_0015*IT_0171;
    const ccomplex_t IT_0173 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0174 = IT_0011*IT_0173;
    const ccomplex_t IT_0175 = IT_0015*IT_0173;
    const ccomplex_t IT_0176 = IT_0011*IT_0171;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0172 + -IT_0174 + 
      -IT_0175 + IT_0176);
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = N_d3*e_em*m_tau*IT_0032*U_stau_01;
    const ccomplex_t IT_0180 = IT_0031*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0183 = IT_0009*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0181 + 2*IT_0184);
    const ccomplex_t IT_0186 = 0.5*IT_0185;
    const ccomplex_t IT_0187 = 2*IT_0186;
    const ccomplex_t IT_0188 = 0.5*IT_0187;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0063*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = IT_0178*IT_0191;
    const ccomplex_t IT_0193 = (-0.5)*IT_0126 + -IT_0147 + -IT_0170 + -IT_0192;
    const ccomplex_t IT_0194 = 2*IT_0193;
    const ccomplex_t IT_0195 = IT_0075*IT_0125;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0009*IT_0010 +
       -IT_0008*IT_0014)*U_stau_01*conj(U_stau_01) + IT_0009*IT_0010*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0197 = (-2)*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0118 + 2*IT_0121);
    const ccomplex_t IT_0199 = 0.5*IT_0198;
    const ccomplex_t IT_0200 = 2*IT_0199;
    const ccomplex_t IT_0201 = 0.5*IT_0200;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = cpow(s_13 + (-0.5)*IT_0004 + (-0.5)*IT_0038 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0204 = IT_0202*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*IT_0204;
    const ccomplex_t IT_0206 = IT_0197*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0009*IT_0010 +
       -IT_0008*IT_0014)*conj(U_stau_00)*U_stau_01 + IT_0009*IT_0010*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0208 = (-2)*IT_0207;
    const ccomplex_t IT_0209 = N_d4*e_em*m_tau*IT_0032*U_stau_00;
    const ccomplex_t IT_0210 = IT_0031*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0213 = IT_0009*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0211 + 2*IT_0214);
    const ccomplex_t IT_0216 = 0.5*IT_0215;
    const ccomplex_t IT_0217 = 2*IT_0216;
    const ccomplex_t IT_0218 = 0.5*IT_0217;
    const ccomplex_t IT_0219 = 2*IT_0218;
    const ccomplex_t IT_0220 = 0.5*IT_0219;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = cpow(s_13 + (-0.5)*IT_0004 + (-0.5)*IT_0038 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0223 = IT_0221*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0208*IT_0224;
    const ccomplex_t IT_0226 = (-0.5)*IT_0195 + (-0.5)*IT_0206 + (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = (-2)*IT_0226;
    const ccomplex_t IT_0228 = IT_0194 + IT_0227;
    const ccomplex_t IT_0229 = 2*conj(IT_0193);
    const ccomplex_t IT_0230 = (-2)*conj(IT_0226);
    const ccomplex_t IT_0231 = IT_0229 + IT_0230;
    const ccomplex_t IT_0232 = IT_0086*IT_0203;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = IT_0197*IT_0233;
    const ccomplex_t IT_0235 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0236 = IT_0009*IT_0235;
    const ccomplex_t IT_0237 = 1.4142135623731*IT_0236;
    const ccomplex_t IT_0238 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0239 = IT_0014*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = conj(N_d4)*e_em*m_tau*IT_0032*U_stau_10;
    const ccomplex_t IT_0242 = IT_0031*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0237 + IT_0240 + 
      -IT_0243);
    const ccomplex_t IT_0245 = (-0.5)*IT_0244;
    const ccomplex_t IT_0246 = IT_0222*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*IT_0246;
    const ccomplex_t IT_0248 = IT_0208*IT_0247;
    const ccomplex_t IT_0249 = IT_0089*IT_0115;
    const ccomplex_t IT_0250 = 0.5*IT_0234 + 0.5*IT_0248 + 0.5*IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0226;
    const ccomplex_t IT_0252 = (-2)*IT_0193;
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 2*conj(IT_0226);
    const ccomplex_t IT_0255 = (-2)*conj(IT_0193);
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = s_13 + IT_0258;
    const ccomplex_t IT_0260 = s_13*IT_0259;
    const ccomplex_t IT_0261 = cpow(s_12 + 0.5*IT_0004 + (-0.5)*IT_0038 + 0.5
      *IT_0039 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0010;
    const ccomplex_t IT_0263 = IT_0202*IT_0262;
    const ccomplex_t IT_0264 = IT_0261*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = 2*IT_0122;
    const ccomplex_t IT_0267 = 0.5*IT_0266;
    const ccomplex_t IT_0268 = 2*IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = IT_0203*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = IT_0197*IT_0272;
    const ccomplex_t IT_0274 = 0.5*IT_0225 + IT_0265 + 0.5*IT_0273;
    const ccomplex_t IT_0275 = e_em*IT_0011;
    const ccomplex_t IT_0276 = e_em*IT_0015;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*(IT_0275 + -IT_0276);
    const ccomplex_t IT_0278 = 0.5*IT_0277;
    const ccomplex_t IT_0279 = IT_0086*IT_0278;
    const ccomplex_t IT_0280 = IT_0261*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = (-0.5)*IT_0234 + (-0.5)*IT_0248 + -IT_0281;
    const ccomplex_t IT_0283 = 2*IT_0114;
    const ccomplex_t IT_0284 = (-2)*IT_0250;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = 2*conj(IT_0114);
    const ccomplex_t IT_0287 = (-2)*conj(IT_0250);
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = m_N_4*IT_0259;
    const ccomplex_t IT_0290 = -IT_0267;
    const ccomplex_t IT_0291 = IT_0087*IT_0290;
    const ccomplex_t IT_0292 = (0 + _Complex_I*1)*IT_0291;
    const ccomplex_t IT_0293 = m_N_4*IT_0115;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = m_N_4*IT_0075;
    const ccomplex_t IT_0296 = IT_0294 + IT_0295;
    const ccomplex_t IT_0297 = IT_0292*IT_0296;
    const ccomplex_t IT_0298 = m_N_1*IT_0098;
    const ccomplex_t IT_0299 = m_N_4*IT_0135;
    const ccomplex_t IT_0300 = -IT_0299;
    const ccomplex_t IT_0301 = IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = IT_0146*IT_0301;
    const ccomplex_t IT_0303 = m_N_2*IT_0022;
    const ccomplex_t IT_0304 = m_N_4*IT_0158;
    const ccomplex_t IT_0305 = -IT_0304;
    const ccomplex_t IT_0306 = IT_0303 + IT_0305;
    const ccomplex_t IT_0307 = IT_0169*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*(0.5*IT_0181 + IT_0184);
    const ccomplex_t IT_0309 = -IT_0308;
    const ccomplex_t IT_0310 = IT_0063*IT_0309;
    const ccomplex_t IT_0311 = (0 + _Complex_I*1)*IT_0310;
    const ccomplex_t IT_0312 = m_N_3*IT_0051;
    const ccomplex_t IT_0313 = m_N_4*IT_0178;
    const ccomplex_t IT_0314 = -IT_0313;
    const ccomplex_t IT_0315 = IT_0312 + IT_0314;
    const ccomplex_t IT_0316 = IT_0311*IT_0315;
    const ccomplex_t IT_0317 = 0.25*IT_0297 + 0.5*IT_0302 + 0.5*IT_0307 + 0.5
      *IT_0316;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0319 = IT_0086*IT_0262;
    const ccomplex_t IT_0320 = IT_0261*IT_0318*IT_0319;
    const ccomplex_t IT_0321 = m_tau*IT_0281;
    const ccomplex_t IT_0322 = IT_0087*IT_0270;
    const ccomplex_t IT_0323 = (0 + _Complex_I*1)*IT_0322;
    const ccomplex_t IT_0324 = IT_0296*IT_0323;
    const ccomplex_t IT_0325 = (-0.5)*IT_0320 + 0.5*IT_0321 + 0.25*IT_0324;
    const ccomplex_t IT_0326 = conj(IT_0317) + conj(IT_0325);
    const ccomplex_t IT_0327 = IT_0317 + IT_0325;
    const ccomplex_t IT_0328 = IT_0089*IT_0296;
    const ccomplex_t IT_0329 = m_N_1*IT_0135;
    const ccomplex_t IT_0330 = m_N_4*IT_0098;
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = IT_0329 + IT_0331;
    const ccomplex_t IT_0333 = IT_0112*IT_0332;
    const ccomplex_t IT_0334 = m_N_4*IT_0022;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = m_N_2*IT_0158;
    const ccomplex_t IT_0337 = IT_0335 + IT_0336;
    const ccomplex_t IT_0338 = IT_0042*IT_0337;
    const ccomplex_t IT_0339 = m_N_4*IT_0051;
    const ccomplex_t IT_0340 = -IT_0339;
    const ccomplex_t IT_0341 = m_N_3*IT_0178;
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = IT_0065*IT_0342;
    const ccomplex_t IT_0344 = 0.25*IT_0328 + (-0.5)*IT_0333 + (-0.5)*IT_0338 
      + (-0.5)*IT_0343;
    const ccomplex_t IT_0345 = IT_0123*IT_0278;
    const ccomplex_t IT_0346 = IT_0261*IT_0318*IT_0345;
    const ccomplex_t IT_0347 = IT_0262*IT_0290;
    const ccomplex_t IT_0348 = IT_0261*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = m_tau*IT_0349;
    const ccomplex_t IT_0351 = 0.25*IT_0328 + 0.5*IT_0346 + (-0.5)*IT_0350;
    const ccomplex_t IT_0352 = conj(IT_0344) + conj(IT_0351);
    const ccomplex_t IT_0353 = IT_0344 + IT_0351;
    const ccomplex_t IT_0354 = s_13*s_14;
    const ccomplex_t IT_0355 = s_34*IT_0004;
    const ccomplex_t IT_0356 = -IT_0355;
    const ccomplex_t IT_0357 = IT_0354 + IT_0356;
    const ccomplex_t IT_0358 = IT_0115*IT_0323;
    const ccomplex_t IT_0359 = 0.5*IT_0147 + 0.5*IT_0170 + 0.5*IT_0192 + 0.25
      *IT_0358;
    const ccomplex_t IT_0360 = 0.5*IT_0349;
    const ccomplex_t IT_0361 = IT_0075*IT_0292;
    const ccomplex_t IT_0362 = 0.25*IT_0361;
    const ccomplex_t IT_0363 = (-0.5)*IT_0043 + (-0.5)*IT_0066 + (-0.25)
      *IT_0090 + (-0.5)*IT_0113;
    const ccomplex_t IT_0364 = (-0.5)*IT_0281;
    const ccomplex_t IT_0365 = conj(IT_0363) + conj(IT_0364);
    const ccomplex_t IT_0366 = IT_0363 + IT_0364;
    const ccomplex_t IT_0367 = (-0.25)*IT_0249;
    const ccomplex_t IT_0368 = 2*IT_0250;
    const ccomplex_t IT_0369 = (-2)*IT_0114;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = 2*conj(IT_0250);
    const ccomplex_t IT_0372 = (-2)*conj(IT_0114);
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = pow(m_Z, 2);
    const ccomplex_t IT_0375 = s_14*IT_0374;
    const ccomplex_t IT_0376 = IT_0001*IT_0375;
    const ccomplex_t IT_0377 = (-0.25)*IT_0376;
    const ccomplex_t IT_0378 = s_14 + IT_0377;
    const ccomplex_t IT_0379 = m_tau*IT_0378;
    const ccomplex_t IT_0380 = (-4)*IT_0379;
    const ccomplex_t IT_0381 = 2*IT_0317;
    const ccomplex_t IT_0382 = 2*IT_0325;
    const ccomplex_t IT_0383 = IT_0381 + IT_0382;
    const ccomplex_t IT_0384 = 2*conj(IT_0317);
    const ccomplex_t IT_0385 = 2*conj(IT_0325);
    const ccomplex_t IT_0386 = IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = 2*IT_0344;
    const ccomplex_t IT_0388 = 2*IT_0351;
    const ccomplex_t IT_0389 = IT_0387 + IT_0388;
    const ccomplex_t IT_0390 = 2*conj(IT_0344);
    const ccomplex_t IT_0391 = 2*conj(IT_0351);
    const ccomplex_t IT_0392 = IT_0390 + IT_0391;
    const ccomplex_t IT_0393 = s_34*IT_0374;
    const ccomplex_t IT_0394 = IT_0001*IT_0393;
    const ccomplex_t IT_0395 = 0.5*IT_0394;
    const ccomplex_t IT_0396 = s_34 + IT_0395;
    const ccomplex_t IT_0397 = m_N_4*IT_0396;
    const ccomplex_t IT_0398 = 2*IT_0397;
    const ccomplex_t IT_0399 = 0.5*IT_0376;
    const ccomplex_t IT_0400 = s_14 + IT_0399;
    const ccomplex_t IT_0401 = m_tau*IT_0400;
    const ccomplex_t IT_0402 = 2*IT_0401;
    const ccomplex_t IT_0403 = (-2)*IT_0317;
    const ccomplex_t IT_0404 = (-2)*IT_0325;
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = (-2)*conj(IT_0317);
    const ccomplex_t IT_0407 = (-2)*conj(IT_0325);
    const ccomplex_t IT_0408 = IT_0406 + IT_0407;
    const ccomplex_t IT_0409 = (-2)*IT_0344;
    const ccomplex_t IT_0410 = (-2)*IT_0351;
    const ccomplex_t IT_0411 = IT_0409 + IT_0410;
    const ccomplex_t IT_0412 = (-2)*conj(IT_0344);
    const ccomplex_t IT_0413 = (-2)*conj(IT_0351);
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = (-0.25)*IT_0394;
    const ccomplex_t IT_0416 = s_34 + IT_0415;
    const ccomplex_t IT_0417 = m_N_4*IT_0416;
    const ccomplex_t IT_0418 = (-4)*IT_0417;
    const ccomplex_t IT_0419 = pow(s_34, 2);
    const ccomplex_t IT_0420 = IT_0001*IT_0419;
    const ccomplex_t IT_0421 = -IT_0038;
    const ccomplex_t IT_0422 = IT_0420 + IT_0421;
    const ccomplex_t IT_0423 = IT_0000*IT_0422;
    const ccomplex_t IT_0424 = 2*IT_0274;
    const ccomplex_t IT_0425 = 2*conj(IT_0274);
    const ccomplex_t IT_0426 = s_13*IT_0422;
    const ccomplex_t IT_0427 = 2*IT_0282;
    const ccomplex_t IT_0428 = m_N_4*IT_0422;
    const ccomplex_t IT_0429 = 2*conj(IT_0282);
    const ccomplex_t IT_0430 = s_13*IT_0374;
    const ccomplex_t IT_0431 = s_14*s_34;
    const ccomplex_t IT_0432 = (-2)*IT_0431;
    const ccomplex_t IT_0433 = IT_0430 + IT_0432;
    const ccomplex_t IT_0434 = IT_0001*IT_0433;
    const ccomplex_t IT_0435 = (-0.5)*IT_0434;
    const ccomplex_t IT_0436 = s_13 + IT_0435;
    const ccomplex_t IT_0437 = 4*IT_0344;
    const ccomplex_t IT_0438 = (-6)*conj(IT_0317);
    const ccomplex_t IT_0439 = s_13*IT_0006;
    const ccomplex_t IT_0440 = IT_0000*IT_0259;
    const ccomplex_t IT_0441 = (-2)*IT_0274;
    const ccomplex_t IT_0442 = (-2)*conj(IT_0274);
    const ccomplex_t IT_0443 = s_13*s_34;
    const ccomplex_t IT_0444 = s_14*IT_0038;
    const ccomplex_t IT_0445 = -IT_0444;
    const ccomplex_t IT_0446 = IT_0443 + IT_0445;
    const ccomplex_t IT_0447 = (-2)*IT_0282;
    const ccomplex_t IT_0448 = (-2)*conj(IT_0282);
    const ccomplex_t IT_0449 = m_tau*IT_0006;
    const ccomplex_t IT_0450 = m_tau*IT_0259;
    const ccomplex_t IT_0451 = IT_0000*IT_0374;
    const ccomplex_t IT_0452 = 3*IT_0451;
    const ccomplex_t IT_0453 = 2*IT_0362;
    const ccomplex_t IT_0454 = (-2)*IT_0359;
    const ccomplex_t IT_0455 = IT_0453 + IT_0454;
    const ccomplex_t IT_0456 = 2*conj(IT_0362);
    const ccomplex_t IT_0457 = (-2)*conj(IT_0359);
    const ccomplex_t IT_0458 = IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = 2*IT_0360;
    const ccomplex_t IT_0460 = 2*conj(IT_0360);
    const ccomplex_t IT_0461 = (-2)*IT_0363;
    const ccomplex_t IT_0462 = (-2)*conj(IT_0363);
    const ccomplex_t IT_0463 = (-3)*IT_0451;
    const ccomplex_t IT_0464 = conj(IT_0359)*IT_0363 + IT_0359*conj(IT_0363) +
       0.5*conj(IT_0367)*IT_0455 + 0.5*IT_0367*IT_0458 + 0.5*conj(IT_0364)
      *IT_0459 + 0.5*IT_0364*IT_0460 + 0.5*conj(IT_0362)*IT_0461 + 0.5*IT_0362
      *IT_0462;
    const ccomplex_t IT_0465 = pow(m_Z, 4);
    const ccomplex_t IT_0466 = s_13*IT_0465;
    const ccomplex_t IT_0467 = IT_0001*IT_0466;
    const ccomplex_t IT_0468 = (-0.25)*IT_0467;
    const ccomplex_t IT_0469 = IT_0431 + IT_0468;
    const ccomplex_t IT_0470 = conj(IT_0359) + -conj(IT_0362);
    const ccomplex_t IT_0471 = IT_0363 + -IT_0367;
    const ccomplex_t IT_0472 = conj(IT_0363) + -conj(IT_0367);
    const ccomplex_t IT_0473 = -IT_0362;
    const ccomplex_t IT_0474 = (-0.5)*IT_0467;
    const ccomplex_t IT_0475 = IT_0430 + IT_0432 + IT_0474;
    const ccomplex_t IT_0476 = IT_0007*(conj(IT_0114)*IT_0228 + IT_0114
      *IT_0231 + conj(IT_0250)*IT_0253 + IT_0250*IT_0256) + IT_0260*(IT_0231
      *IT_0274 + IT_0228*conj(IT_0274) + conj(IT_0282)*IT_0285 + IT_0282*IT_0288
      ) + IT_0289*(IT_0228*IT_0326 + IT_0231*IT_0327 + IT_0285*IT_0352 + IT_0288
      *IT_0353) + IT_0357*(IT_0231*(IT_0359 + IT_0360) + IT_0228*(conj(IT_0359) 
      + conj(IT_0360)) + IT_0256*IT_0362 + IT_0253*conj(IT_0362) + IT_0285
      *IT_0365 + IT_0288*IT_0366 + conj(IT_0367)*IT_0370 + IT_0367*IT_0373) +
       IT_0380*(conj(IT_0364)*IT_0383 + IT_0364*IT_0386 + conj(IT_0360)*IT_0389 
      + IT_0360*IT_0392) + (conj(IT_0360)*IT_0383 + IT_0360*IT_0386 + conj
      (IT_0364)*IT_0389 + IT_0364*IT_0392)*IT_0398 + IT_0402*(conj(IT_0367)
      *IT_0383 + IT_0367*IT_0386 + conj(IT_0362)*IT_0389 + IT_0362*IT_0392 +
       conj(IT_0363)*IT_0405 + IT_0363*IT_0408 + conj(IT_0359)*IT_0411 + IT_0359
      *IT_0414) + (conj(IT_0362)*IT_0383 + IT_0362*IT_0386 + conj(IT_0367)
      *IT_0389 + IT_0367*IT_0392 + conj(IT_0359)*IT_0405 + IT_0359*IT_0408 +
       conj(IT_0363)*IT_0411 + IT_0363*IT_0414)*IT_0418 + IT_0423*(conj(IT_0282)
      *IT_0424 + IT_0282*IT_0425) + IT_0426*(conj(IT_0274)*IT_0424 + conj
      (IT_0282)*IT_0427) + IT_0428*(IT_0326*IT_0424 + IT_0327*IT_0425 + IT_0352
      *IT_0427 + IT_0353*IT_0429) + 4*IT_0436*(IT_0326*IT_0327 + IT_0352*
      (IT_0351 + 0.25*IT_0437)) + (-6)*IT_0000*(IT_0327*IT_0352 + IT_0353*(conj
      (IT_0325) + (-0.166666666666667)*IT_0438)) + (conj(IT_0193)*IT_0194 + conj
      (IT_0226)*IT_0253 + IT_0226*IT_0255 + conj(IT_0114)*IT_0283 + conj(IT_0250
      )*IT_0370 + IT_0250*IT_0372)*IT_0439 + IT_0440*(IT_0231*IT_0282 + IT_0228
      *conj(IT_0282) + conj(IT_0114)*IT_0424 + IT_0114*IT_0425 + conj(IT_0250)
      *IT_0441 + IT_0250*IT_0442) + IT_0446*((conj(IT_0359) + conj(IT_0360))
      *IT_0424 + (IT_0359 + IT_0360)*IT_0425 + IT_0365*IT_0427 + IT_0366*IT_0429
       + conj(IT_0362)*IT_0441 + IT_0362*IT_0442 + conj(IT_0367)*IT_0447 +
       IT_0367*IT_0448) + (IT_0253*IT_0352 + IT_0256*IT_0353 + conj(IT_0325)
      *IT_0370 + IT_0325*IT_0373 + conj(IT_0250)*IT_0381 + IT_0250*IT_0384 +
       conj(IT_0114)*IT_0403 + IT_0114*IT_0406)*IT_0449 + (conj(IT_0282)*IT_0403
       + IT_0282*IT_0406 + IT_0352*IT_0441 + IT_0353*IT_0442 + conj(IT_0325)
      *IT_0447 + IT_0325*IT_0448)*IT_0450 + IT_0452*(conj(IT_0364)*IT_0455 +
       IT_0364*IT_0458 + conj(IT_0367)*IT_0459 + IT_0367*IT_0460 + conj(IT_0360)
      *IT_0461 + IT_0360*IT_0462) + 2*IT_0463*IT_0464 + 8*IT_0469*(IT_0360
      *IT_0470 + conj(IT_0364)*IT_0471 + IT_0364*IT_0472 + conj(IT_0360)*
      (IT_0359 + IT_0473)) + (-4)*(IT_0360*conj(IT_0360) + IT_0364*conj(IT_0364)
       + IT_0471*IT_0472 + IT_0470*(IT_0359 + IT_0473))*IT_0475;
    return create_ccomplex_return(IT_0476);
}

