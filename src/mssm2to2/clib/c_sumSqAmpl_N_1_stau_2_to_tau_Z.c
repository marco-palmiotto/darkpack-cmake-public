#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_stau_2_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_stau_2_to_tau_Z(
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
    const ccomplex_t IT_0007 = pow(m_N_1, 2);
    const ccomplex_t IT_0008 = pow(m_tau, 2);
    const ccomplex_t IT_0009 = pow(m_stau_2, 2);
    const ccomplex_t IT_0010 = cpow(s_12 + 0.5*IT_0007 + (-0.5)*IT_0008 + 0.5
      *IT_0009 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = cpow(IT_0011, -1);
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*IT_0012*IT_0013;
    const ccomplex_t IT_0015 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0016 = IT_0012*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = cos(beta);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = cpow(IT_0013, -1);
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = pow(m_W, -1);
    const ccomplex_t IT_0023 = N_d1*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*(IT_0017 + 0.5*IT_0025);
    const ccomplex_t IT_0027 = -IT_0026;
    const ccomplex_t IT_0028 = IT_0014*IT_0027;
    const ccomplex_t IT_0029 = IT_0010*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = IT_0012*IT_0013;
    const ccomplex_t IT_0033 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = IT_0011*IT_0020;
    const ccomplex_t IT_0036 = IT_0033*IT_0035;
    const ccomplex_t IT_0037 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0038 = IT_0032*IT_0037;
    const ccomplex_t IT_0039 = IT_0035*IT_0037;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0034 + IT_0036 + 
      -IT_0038 + -IT_0039);
    const ccomplex_t IT_0041 = 0.5*IT_0040;
    const ccomplex_t IT_0042 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0043 = IT_0027*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*IT_0043;
    const ccomplex_t IT_0045 = IT_0041*IT_0044;
    const ccomplex_t IT_0046 = 0.25*IT_0045;
    const ccomplex_t IT_0047 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0048 = IT_0032*IT_0047;
    const ccomplex_t IT_0049 = IT_0035*IT_0047;
    const ccomplex_t IT_0050 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0051 = IT_0032*IT_0050;
    const ccomplex_t IT_0052 = IT_0035*IT_0050;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0048 + IT_0049 + 
      -IT_0051 + -IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0056 = IT_0012*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = N_d2*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0059 = IT_0021*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 0.5*IT_0060);
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0068 = IT_0032*IT_0067;
    const ccomplex_t IT_0069 = IT_0035*IT_0067;
    const ccomplex_t IT_0070 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0071 = IT_0032*IT_0070;
    const ccomplex_t IT_0072 = IT_0035*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0076 = IT_0012*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_d3*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 0.5*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0088 = IT_0032*IT_0087;
    const ccomplex_t IT_0089 = IT_0035*IT_0087;
    const ccomplex_t IT_0090 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0091 = IT_0032*IT_0090;
    const ccomplex_t IT_0092 = IT_0035*IT_0090;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0088 + IT_0089 + 
      -IT_0091 + -IT_0092);
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0096 = IT_0012*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = N_d4*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0099 = IT_0021*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + 0.5*IT_0100);
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0094*IT_0105;
    const ccomplex_t IT_0107 = 0.5*IT_0066 + 0.5*IT_0086 + 0.5*IT_0106;
    const ccomplex_t IT_0108 = (-0.5)*IT_0040;
    const ccomplex_t IT_0109 = IT_0044*IT_0108;
    const ccomplex_t IT_0110 = 0.25*IT_0109;
    const ccomplex_t IT_0111 = e_em*IT_0032;
    const ccomplex_t IT_0112 = e_em*IT_0035;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0111 + -IT_0112);
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0116 = IT_0012*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0119 = IT_0020*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = conj(N_d1)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0122 = IT_0021*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0117 + IT_0120 + 
      -IT_0123);
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = IT_0114*IT_0125;
    const ccomplex_t IT_0127 = IT_0010*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = IT_0042*IT_0125;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0108*IT_0131;
    const ccomplex_t IT_0133 = (-0.25)*IT_0132;
    const ccomplex_t IT_0134 = IT_0041*IT_0131;
    const ccomplex_t IT_0135 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0136 = IT_0032*IT_0135;
    const ccomplex_t IT_0137 = IT_0035*IT_0135;
    const ccomplex_t IT_0138 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0139 = IT_0032*IT_0138;
    const ccomplex_t IT_0140 = IT_0035*IT_0138;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*(IT_0136 + IT_0137 + 
      -IT_0139 + -IT_0140);
    const ccomplex_t IT_0142 = (-0.5)*IT_0141;
    const ccomplex_t IT_0143 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0144 = IT_0012*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0147 = IT_0020*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = conj(N_d2)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0150 = IT_0021*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0145 + IT_0148 + 
      -IT_0151);
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = IT_0063*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = IT_0142*IT_0155;
    const ccomplex_t IT_0157 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0158 = IT_0032*IT_0157;
    const ccomplex_t IT_0159 = IT_0035*IT_0157;
    const ccomplex_t IT_0160 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0161 = IT_0032*IT_0160;
    const ccomplex_t IT_0162 = IT_0035*IT_0160;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0158 + IT_0159 + 
      -IT_0161 + -IT_0162);
    const ccomplex_t IT_0164 = (-0.5)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0166 = IT_0012*IT_0165;
    const ccomplex_t IT_0167 = 1.4142135623731*IT_0166;
    const ccomplex_t IT_0168 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0169 = IT_0020*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_d3)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0172 = IT_0021*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + IT_0170 + 
      -IT_0173);
    const ccomplex_t IT_0175 = (-0.5)*IT_0174;
    const ccomplex_t IT_0176 = IT_0083*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0164*IT_0177;
    const ccomplex_t IT_0179 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0180 = IT_0032*IT_0179;
    const ccomplex_t IT_0181 = IT_0035*IT_0179;
    const ccomplex_t IT_0182 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0183 = IT_0032*IT_0182;
    const ccomplex_t IT_0184 = IT_0035*IT_0182;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0180 + IT_0181 + 
      -IT_0183 + -IT_0184);
    const ccomplex_t IT_0186 = (-0.5)*IT_0185;
    const ccomplex_t IT_0187 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0188 = IT_0012*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0191 = IT_0020*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = conj(N_d4)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0194 = IT_0021*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0189 + IT_0192 + 
      -IT_0195);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0103*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*IT_0198;
    const ccomplex_t IT_0200 = IT_0186*IT_0199;
    const ccomplex_t IT_0201 = (-0.25)*IT_0134 + (-0.5)*IT_0156 + (-0.5)
      *IT_0178 + (-0.5)*IT_0200;
    const ccomplex_t IT_0202 = IT_0133 + -IT_0201;
    const ccomplex_t IT_0203 = conj(IT_0133) + -conj(IT_0201);
    const ccomplex_t IT_0204 = m_N_1*m_tau;
    const ccomplex_t IT_0205 = pow(s_14, 2);
    const ccomplex_t IT_0206 = IT_0001*IT_0205;
    const ccomplex_t IT_0207 = -IT_0007;
    const ccomplex_t IT_0208 = IT_0206 + IT_0207;
    const ccomplex_t IT_0209 = IT_0204*IT_0208;
    const ccomplex_t IT_0210 = 0.5*IT_0132;
    const ccomplex_t IT_0211 = (-0.5)*IT_0045;
    const ccomplex_t IT_0212 = 2*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0012*IT_0013 +
       -IT_0011*IT_0020)*U_stau_01*conj(U_stau_01) + IT_0012*IT_0013*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0214 = (-2)*IT_0213;
    const ccomplex_t IT_0215 = cpow(s_13 + (-0.5)*IT_0007 + (-0.5)*IT_0008 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0216 = IT_0027*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0214*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0012*IT_0013 +
       -IT_0011*IT_0020)*conj(U_stau_00)*U_stau_01 + IT_0012*IT_0013*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0220 = (-2)*IT_0219;
    const ccomplex_t IT_0221 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0222 = IT_0012*IT_0221;
    const ccomplex_t IT_0223 = 1.4142135623731*IT_0222;
    const ccomplex_t IT_0224 = N_d1*e_em*m_tau*IT_0022*U_stau_00;
    const ccomplex_t IT_0225 = IT_0021*IT_0224;
    const ccomplex_t IT_0226 = 1.4142135623731*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*(IT_0223 + 0.5*IT_0226);
    const ccomplex_t IT_0228 = -IT_0227;
    const ccomplex_t IT_0229 = cpow(s_13 + (-0.5)*IT_0007 + (-0.5)*IT_0008 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0230 = IT_0228*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = IT_0220*IT_0231;
    const ccomplex_t IT_0233 = -IT_0066 + -IT_0086 + -IT_0106 + (-0.5)*IT_0109
       + (-0.5)*IT_0218 + (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = (-2)*IT_0233;
    const ccomplex_t IT_0235 = IT_0212 + IT_0234;
    const ccomplex_t IT_0236 = 2*conj(IT_0211);
    const ccomplex_t IT_0237 = (-2)*conj(IT_0233);
    const ccomplex_t IT_0238 = IT_0236 + IT_0237;
    const ccomplex_t IT_0239 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0240 = IT_0012*IT_0239;
    const ccomplex_t IT_0241 = 1.4142135623731*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0243 = IT_0020*IT_0242;
    const ccomplex_t IT_0244 = 1.4142135623731*IT_0243;
    const ccomplex_t IT_0245 = conj(N_d1)*e_em*m_tau*IT_0022*U_stau_10;
    const ccomplex_t IT_0246 = IT_0021*IT_0245;
    const ccomplex_t IT_0247 = 1.4142135623731*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0241 + IT_0244 + 
      -IT_0247);
    const ccomplex_t IT_0249 = (-0.5)*IT_0248;
    const ccomplex_t IT_0250 = IT_0229*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0220*IT_0251;
    const ccomplex_t IT_0253 = IT_0125*IT_0215;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0214*IT_0254;
    const ccomplex_t IT_0256 = 0.5*IT_0134 + IT_0156 + IT_0178 + IT_0200 + 0.5
      *IT_0252 + 0.5*IT_0255;
    const ccomplex_t IT_0257 = 2*IT_0233;
    const ccomplex_t IT_0258 = (-2)*IT_0211;
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = 2*conj(IT_0233);
    const ccomplex_t IT_0261 = (-2)*conj(IT_0211);
    const ccomplex_t IT_0262 = IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = s_13 + IT_0264;
    const ccomplex_t IT_0266 = s_13*IT_0265;
    const ccomplex_t IT_0267 = IT_0030 + 0.5*IT_0218 + 0.5*IT_0232;
    const ccomplex_t IT_0268 = -IT_0128 + (-0.5)*IT_0252 + (-0.5)*IT_0255;
    const ccomplex_t IT_0269 = 2*IT_0210;
    const ccomplex_t IT_0270 = (-2)*IT_0256;
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = 2*conj(IT_0210);
    const ccomplex_t IT_0273 = (-2)*conj(IT_0256);
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = m_N_1*IT_0265;
    const ccomplex_t IT_0276 = m_N_1*IT_0108;
    const ccomplex_t IT_0277 = m_N_1*IT_0041;
    const ccomplex_t IT_0278 = IT_0044*(IT_0276 + -IT_0277);
    const ccomplex_t IT_0279 = 0.25*IT_0278;
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0281 = IT_0014*IT_0125;
    const ccomplex_t IT_0282 = IT_0010*IT_0280*IT_0281;
    const ccomplex_t IT_0283 = m_tau*IT_0128;
    const ccomplex_t IT_0284 = -IT_0277;
    const ccomplex_t IT_0285 = IT_0276 + IT_0284;
    const ccomplex_t IT_0286 = IT_0044*IT_0285;
    const ccomplex_t IT_0287 = m_N_2*IT_0142;
    const ccomplex_t IT_0288 = m_N_1*IT_0054;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = IT_0065*IT_0290;
    const ccomplex_t IT_0292 = m_N_3*IT_0164;
    const ccomplex_t IT_0293 = m_N_1*IT_0074;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = IT_0085*IT_0295;
    const ccomplex_t IT_0297 = m_N_4*IT_0186;
    const ccomplex_t IT_0298 = m_N_1*IT_0094;
    const ccomplex_t IT_0299 = -IT_0298;
    const ccomplex_t IT_0300 = IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = IT_0105*IT_0300;
    const ccomplex_t IT_0302 = (-0.5)*IT_0282 + 0.5*IT_0283 + 0.25*IT_0286 + (
      -0.5)*IT_0291 + (-0.5)*IT_0296 + (-0.5)*IT_0301;
    const ccomplex_t IT_0303 = conj(IT_0279) + conj(IT_0302);
    const ccomplex_t IT_0304 = IT_0279 + IT_0302;
    const ccomplex_t IT_0305 = IT_0131*(IT_0276 + -IT_0277);
    const ccomplex_t IT_0306 = 0.25*IT_0305;
    const ccomplex_t IT_0307 = IT_0027*IT_0114;
    const ccomplex_t IT_0308 = IT_0010*IT_0280*IT_0307;
    const ccomplex_t IT_0309 = IT_0131*IT_0285;
    const ccomplex_t IT_0310 = m_N_2*IT_0054;
    const ccomplex_t IT_0311 = m_N_1*IT_0142;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = IT_0310 + IT_0312;
    const ccomplex_t IT_0314 = IT_0155*IT_0313;
    const ccomplex_t IT_0315 = m_N_3*IT_0074;
    const ccomplex_t IT_0316 = m_N_1*IT_0164;
    const ccomplex_t IT_0317 = -IT_0316;
    const ccomplex_t IT_0318 = IT_0315 + IT_0317;
    const ccomplex_t IT_0319 = IT_0177*IT_0318;
    const ccomplex_t IT_0320 = m_N_4*IT_0094;
    const ccomplex_t IT_0321 = m_N_1*IT_0186;
    const ccomplex_t IT_0322 = -IT_0321;
    const ccomplex_t IT_0323 = IT_0320 + IT_0322;
    const ccomplex_t IT_0324 = IT_0199*IT_0323;
    const ccomplex_t IT_0325 = m_tau*IT_0030;
    const ccomplex_t IT_0326 = 0.5*IT_0308 + 0.25*IT_0309 + 0.5*IT_0314 + 0.5
      *IT_0319 + 0.5*IT_0324 + (-0.5)*IT_0325;
    const ccomplex_t IT_0327 = conj(IT_0306) + conj(IT_0326);
    const ccomplex_t IT_0328 = IT_0306 + IT_0326;
    const ccomplex_t IT_0329 = s_13*s_14;
    const ccomplex_t IT_0330 = s_34*IT_0007;
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = IT_0329 + IT_0331;
    const ccomplex_t IT_0333 = conj(IT_0031) + conj(IT_0046);
    const ccomplex_t IT_0334 = IT_0031 + IT_0046;
    const ccomplex_t IT_0335 = conj(IT_0107) + conj(IT_0110);
    const ccomplex_t IT_0336 = IT_0107 + IT_0110;
    const ccomplex_t IT_0337 = conj(IT_0129) + conj(IT_0133);
    const ccomplex_t IT_0338 = IT_0129 + IT_0133;
    const ccomplex_t IT_0339 = 2*IT_0256;
    const ccomplex_t IT_0340 = (-2)*IT_0210;
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = 2*conj(IT_0256);
    const ccomplex_t IT_0343 = (-2)*conj(IT_0210);
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = pow(m_Z, 2);
    const ccomplex_t IT_0346 = s_34*IT_0345;
    const ccomplex_t IT_0347 = IT_0001*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0347;
    const ccomplex_t IT_0349 = s_34 + IT_0348;
    const ccomplex_t IT_0350 = m_N_1*IT_0349;
    const ccomplex_t IT_0351 = 2*IT_0350;
    const ccomplex_t IT_0352 = 2*IT_0279;
    const ccomplex_t IT_0353 = 2*IT_0302;
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = 2*conj(IT_0279);
    const ccomplex_t IT_0356 = 2*conj(IT_0302);
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = 2*IT_0306;
    const ccomplex_t IT_0359 = 2*IT_0326;
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = 2*conj(IT_0306);
    const ccomplex_t IT_0362 = 2*conj(IT_0326);
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = s_14*IT_0345;
    const ccomplex_t IT_0365 = IT_0001*IT_0364;
    const ccomplex_t IT_0366 = (-0.25)*IT_0365;
    const ccomplex_t IT_0367 = s_14 + IT_0366;
    const ccomplex_t IT_0368 = m_tau*IT_0367;
    const ccomplex_t IT_0369 = (-4)*IT_0368;
    const ccomplex_t IT_0370 = 0.5*IT_0365;
    const ccomplex_t IT_0371 = s_14 + IT_0370;
    const ccomplex_t IT_0372 = m_tau*IT_0371;
    const ccomplex_t IT_0373 = 2*IT_0372;
    const ccomplex_t IT_0374 = (-2)*IT_0279;
    const ccomplex_t IT_0375 = (-2)*IT_0302;
    const ccomplex_t IT_0376 = IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = (-2)*conj(IT_0279);
    const ccomplex_t IT_0378 = (-2)*conj(IT_0302);
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = (-2)*IT_0306;
    const ccomplex_t IT_0381 = (-2)*IT_0326;
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = (-2)*conj(IT_0306);
    const ccomplex_t IT_0384 = (-2)*conj(IT_0326);
    const ccomplex_t IT_0385 = IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = (-0.25)*IT_0347;
    const ccomplex_t IT_0387 = s_34 + IT_0386;
    const ccomplex_t IT_0388 = m_N_1*IT_0387;
    const ccomplex_t IT_0389 = (-4)*IT_0388;
    const ccomplex_t IT_0390 = pow(s_34, 2);
    const ccomplex_t IT_0391 = IT_0001*IT_0390;
    const ccomplex_t IT_0392 = -IT_0008;
    const ccomplex_t IT_0393 = IT_0391 + IT_0392;
    const ccomplex_t IT_0394 = IT_0204*IT_0393;
    const ccomplex_t IT_0395 = 2*IT_0267;
    const ccomplex_t IT_0396 = 2*conj(IT_0267);
    const ccomplex_t IT_0397 = s_13*IT_0393;
    const ccomplex_t IT_0398 = 2*IT_0268;
    const ccomplex_t IT_0399 = m_N_1*IT_0393;
    const ccomplex_t IT_0400 = 2*conj(IT_0268);
    const ccomplex_t IT_0401 = s_13*IT_0345;
    const ccomplex_t IT_0402 = (-2)*IT_0000;
    const ccomplex_t IT_0403 = IT_0401 + IT_0402;
    const ccomplex_t IT_0404 = IT_0001*IT_0403;
    const ccomplex_t IT_0405 = (-0.5)*IT_0404;
    const ccomplex_t IT_0406 = s_13 + IT_0405;
    const ccomplex_t IT_0407 = 4*IT_0306;
    const ccomplex_t IT_0408 = (-6)*conj(IT_0279);
    const ccomplex_t IT_0409 = (-0.5)*IT_0004;
    const ccomplex_t IT_0410 = IT_0401 + IT_0402 + IT_0409;
    const ccomplex_t IT_0411 = 4*conj(IT_0046);
    const ccomplex_t IT_0412 = 4*IT_0046;
    const ccomplex_t IT_0413 = (-4)*IT_0107;
    const ccomplex_t IT_0414 = m_tau*IT_0208;
    const ccomplex_t IT_0415 = s_13*IT_0208;
    const ccomplex_t IT_0416 = IT_0204*IT_0265;
    const ccomplex_t IT_0417 = (-2)*IT_0267;
    const ccomplex_t IT_0418 = (-2)*conj(IT_0267);
    const ccomplex_t IT_0419 = m_tau*IT_0265;
    const ccomplex_t IT_0420 = (-2)*IT_0268;
    const ccomplex_t IT_0421 = (-2)*conj(IT_0268);
    const ccomplex_t IT_0422 = s_13*s_34;
    const ccomplex_t IT_0423 = s_14*IT_0008;
    const ccomplex_t IT_0424 = -IT_0423;
    const ccomplex_t IT_0425 = IT_0422 + IT_0424;
    const ccomplex_t IT_0426 = IT_0204*IT_0345;
    const ccomplex_t IT_0427 = 3*IT_0426;
    const ccomplex_t IT_0428 = 2*IT_0107;
    const ccomplex_t IT_0429 = 2*IT_0110;
    const ccomplex_t IT_0430 = (-2)*IT_0046;
    const ccomplex_t IT_0431 = IT_0428 + IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 2*conj(IT_0107);
    const ccomplex_t IT_0433 = 2*conj(IT_0110);
    const ccomplex_t IT_0434 = (-2)*conj(IT_0046);
    const ccomplex_t IT_0435 = IT_0432 + IT_0433 + IT_0434;
    const ccomplex_t IT_0436 = 2*IT_0031;
    const ccomplex_t IT_0437 = 2*conj(IT_0031);
    const ccomplex_t IT_0438 = (-2)*IT_0133;
    const ccomplex_t IT_0439 = (-2)*conj(IT_0133);
    const ccomplex_t IT_0440 = (-3)*IT_0426;
    const ccomplex_t IT_0441 = conj(IT_0046)*IT_0133 + IT_0046*conj(IT_0133) +
       0.5*conj(IT_0201)*IT_0431 + 0.5*IT_0201*IT_0435 + 0.5*conj(IT_0129)
      *IT_0436 + 0.5*IT_0129*IT_0437 + 0.5*IT_0335*IT_0438 + 0.5*IT_0336*IT_0439;
    const ccomplex_t IT_0442 = 8*IT_0006*(conj(IT_0031)*(IT_0046 + -IT_0107 + 
      -IT_0110) + IT_0031*(conj(IT_0046) + -conj(IT_0107) + -conj(IT_0110)) +
       conj(IT_0129)*IT_0202 + IT_0129*IT_0203) + IT_0209*(conj(IT_0210)*IT_0235
       + IT_0210*IT_0238 + conj(IT_0256)*IT_0259 + IT_0256*IT_0262) + IT_0266*
      (IT_0238*IT_0267 + IT_0235*conj(IT_0267) + conj(IT_0268)*IT_0271 + IT_0268
      *IT_0274) + IT_0275*(IT_0235*IT_0303 + IT_0238*IT_0304 + IT_0271*IT_0327 +
       IT_0274*IT_0328) + IT_0332*(IT_0235*IT_0333 + IT_0238*IT_0334 + IT_0259
      *IT_0335 + IT_0262*IT_0336 + IT_0271*IT_0337 + IT_0274*IT_0338 + conj
      (IT_0201)*IT_0341 + IT_0201*IT_0344) + IT_0351*(conj(IT_0031)*IT_0354 +
       IT_0031*IT_0357 + conj(IT_0129)*IT_0360 + IT_0129*IT_0363) + (conj
      (IT_0129)*IT_0354 + IT_0129*IT_0357 + conj(IT_0031)*IT_0360 + IT_0031
      *IT_0363)*IT_0369 + IT_0373*(conj(IT_0201)*IT_0354 + IT_0201*IT_0357 +
       IT_0335*IT_0360 + IT_0336*IT_0363 + conj(IT_0133)*IT_0376 + IT_0133
      *IT_0379 + conj(IT_0046)*IT_0382 + IT_0046*IT_0385) + (IT_0335*IT_0354 +
       IT_0336*IT_0357 + conj(IT_0201)*IT_0360 + IT_0201*IT_0363 + conj(IT_0046)
      *IT_0376 + IT_0046*IT_0379 + conj(IT_0133)*IT_0382 + IT_0133*IT_0385)
      *IT_0389 + IT_0394*(conj(IT_0268)*IT_0395 + IT_0268*IT_0396) + IT_0397*
      (conj(IT_0267)*IT_0395 + conj(IT_0268)*IT_0398) + IT_0399*(IT_0303*IT_0395
       + IT_0304*IT_0396 + IT_0327*IT_0398 + IT_0328*IT_0400) + 4*IT_0406*
      (IT_0303*IT_0304 + IT_0327*(IT_0326 + 0.25*IT_0407)) + (-6)*IT_0204*
      (IT_0304*IT_0327 + IT_0328*(conj(IT_0302) + (-0.166666666666667)*IT_0408))
       + (-4)*IT_0410*(IT_0031*conj(IT_0031) + IT_0046*conj(IT_0046) + IT_0129
      *conj(IT_0129) + IT_0202*IT_0203 + (-0.25)*IT_0336*IT_0411 + IT_0335*
      (IT_0110 + (-0.25)*IT_0412 + (-0.25)*IT_0413)) + (IT_0259*IT_0327 +
       IT_0262*IT_0328 + conj(IT_0302)*IT_0341 + IT_0302*IT_0344 + conj(IT_0256)
      *IT_0352 + IT_0256*IT_0355 + conj(IT_0210)*IT_0374 + IT_0210*IT_0377)
      *IT_0414 + (conj(IT_0211)*IT_0212 + conj(IT_0233)*IT_0259 + IT_0233
      *IT_0261 + conj(IT_0210)*IT_0269 + conj(IT_0256)*IT_0341 + IT_0256*IT_0343
      )*IT_0415 + IT_0416*(IT_0238*IT_0268 + IT_0235*conj(IT_0268) + conj
      (IT_0210)*IT_0395 + IT_0210*IT_0396 + conj(IT_0256)*IT_0417 + IT_0256
      *IT_0418) + IT_0419*(conj(IT_0268)*IT_0374 + IT_0268*IT_0377 + IT_0327
      *IT_0417 + IT_0328*IT_0418 + conj(IT_0302)*IT_0420 + IT_0302*IT_0421) + 
      (IT_0333*IT_0395 + IT_0334*IT_0396 + IT_0337*IT_0398 + IT_0338*IT_0400 +
       IT_0335*IT_0417 + IT_0336*IT_0418 + conj(IT_0201)*IT_0420 + IT_0201
      *IT_0421)*IT_0425 + IT_0427*(conj(IT_0129)*IT_0431 + IT_0129*IT_0435 +
       conj(IT_0201)*IT_0436 + IT_0201*IT_0437 + conj(IT_0031)*IT_0438 + IT_0031
      *IT_0439) + 2*IT_0440*IT_0441;
    return create_ccomplex_return(IT_0442);
}

