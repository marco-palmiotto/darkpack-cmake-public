#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_stau_1_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_stau_1_to_tau_Z(
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
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = IT_0007*IT_0009;
    const ccomplex_t IT_0011 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = cpow(IT_0007, -1);
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0016 = IT_0014*IT_0015;
    const ccomplex_t IT_0017 = IT_0010*IT_0015;
    const ccomplex_t IT_0018 = IT_0011*IT_0014;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*(IT_0012 + -IT_0016 + 
      -IT_0017 + IT_0018);
    const ccomplex_t IT_0020 = -IT_0019;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0023 = IT_0013*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = cos(beta);
    const ccomplex_t IT_0026 = cpow(IT_0025, -1);
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = pow(m_W, -1);
    const ccomplex_t IT_0029 = N_d4*e_em*m_tau*IT_0028*U_stau_00;
    const ccomplex_t IT_0030 = IT_0027*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0024 + 0.5*IT_0031);
    const ccomplex_t IT_0033 = 2*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = pow(m_tau, 2);
    const ccomplex_t IT_0037 = pow(m_stau_1, 2);
    const ccomplex_t IT_0038 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0035*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = IT_0021*IT_0040;
    const ccomplex_t IT_0042 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0043 = IT_0010*IT_0042;
    const ccomplex_t IT_0044 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0045 = IT_0014*IT_0044;
    const ccomplex_t IT_0046 = IT_0010*IT_0044;
    const ccomplex_t IT_0047 = IT_0014*IT_0042;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0043 + -IT_0045 + 
      -IT_0046 + IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = N_d1*e_em*m_tau*IT_0028*U_stau_00;
    const ccomplex_t IT_0051 = IT_0027*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0054 = IT_0013*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + 2*IT_0055);
    const ccomplex_t IT_0057 = 0.5*IT_0056;
    const ccomplex_t IT_0058 = -IT_0057;
    const ccomplex_t IT_0059 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = IT_0049*IT_0061;
    const ccomplex_t IT_0063 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0064 = IT_0014*IT_0063;
    const ccomplex_t IT_0065 = IT_0010*IT_0063;
    const ccomplex_t IT_0066 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0067 = IT_0010*IT_0066;
    const ccomplex_t IT_0068 = IT_0014*IT_0066;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*(IT_0064 + IT_0065 + 
      -IT_0067 + -IT_0068);
    const ccomplex_t IT_0070 = -IT_0069;
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = (-0.5)*IT_0072;
    const ccomplex_t IT_0074 = N_d2*e_em*m_tau*IT_0028*U_stau_00;
    const ccomplex_t IT_0075 = IT_0027*IT_0074;
    const ccomplex_t IT_0076 = 1.4142135623731*IT_0075;
    const ccomplex_t IT_0077 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0078 = IT_0013*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*(IT_0076 + 2*IT_0079);
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0073*IT_0085;
    const ccomplex_t IT_0087 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0088 = IT_0010*IT_0087;
    const ccomplex_t IT_0089 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0090 = IT_0014*IT_0089;
    const ccomplex_t IT_0091 = IT_0010*IT_0089;
    const ccomplex_t IT_0092 = IT_0014*IT_0087;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0088 + -IT_0090 + 
      -IT_0091 + IT_0092);
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0096 = IT_0013*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = N_d3*e_em*m_tau*IT_0028*U_stau_00;
    const ccomplex_t IT_0099 = IT_0027*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*(IT_0097 + 0.5*IT_0100);
    const ccomplex_t IT_0102 = -IT_0101;
    const ccomplex_t IT_0103 = cpow(s_23 + (-0.5)*IT_0036 + (-0.5)*IT_0037 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = IT_0094*IT_0105;
    const ccomplex_t IT_0107 = 0.25*IT_0041 + 0.5*IT_0062 + 0.5*IT_0086 + 0.5
      *IT_0106;
    const ccomplex_t IT_0108 = (-0.5)*IT_0020;
    const ccomplex_t IT_0109 = IT_0040*IT_0108;
    const ccomplex_t IT_0110 = 0.25*IT_0109;
    const ccomplex_t IT_0111 = pow(m_N_4, 2);
    const ccomplex_t IT_0112 = cpow(s_12 + (-0.5)*IT_0036 + 0.5*IT_0037 + 0.5
      *IT_0111 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0013;
    const ccomplex_t IT_0114 = IT_0035*IT_0113;
    const ccomplex_t IT_0115 = IT_0112*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = 0.5*IT_0116;
    const ccomplex_t IT_0118 = -conj(IT_0110);
    const ccomplex_t IT_0119 = conj(IT_0107) + IT_0118;
    const ccomplex_t IT_0120 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0121 = IT_0013*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0124 = IT_0009*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = conj(N_d4)*e_em*m_tau*IT_0028*U_stau_10;
    const ccomplex_t IT_0127 = IT_0027*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0122 + IT_0125 + 
      -IT_0128);
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = e_em*IT_0014;
    const ccomplex_t IT_0132 = e_em*IT_0010;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0131 + -IT_0132);
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = IT_0130*IT_0134;
    const ccomplex_t IT_0136 = IT_0112*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = (-0.5)*IT_0137;
    const ccomplex_t IT_0139 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0140 = IT_0014*IT_0139;
    const ccomplex_t IT_0141 = IT_0010*IT_0139;
    const ccomplex_t IT_0142 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0143 = IT_0014*IT_0142;
    const ccomplex_t IT_0144 = IT_0010*IT_0142;
    const ccomplex_t IT_0145 = (0 + _Complex_I*1)*(IT_0140 + IT_0141 + 
      -IT_0143 + -IT_0144);
    const ccomplex_t IT_0146 = (-0.5)*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0148 = IT_0013*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0151 = IT_0009*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = conj(N_d3)*e_em*m_tau*IT_0028*U_stau_10;
    const ccomplex_t IT_0154 = IT_0027*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0149 + IT_0152 + 
      -IT_0155);
    const ccomplex_t IT_0157 = (-0.5)*IT_0156;
    const ccomplex_t IT_0158 = IT_0103*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = IT_0146*IT_0159;
    const ccomplex_t IT_0161 = IT_0038*IT_0130;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*IT_0161;
    const ccomplex_t IT_0163 = IT_0108*IT_0162;
    const ccomplex_t IT_0164 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0165 = IT_0014*IT_0164;
    const ccomplex_t IT_0166 = IT_0010*IT_0164;
    const ccomplex_t IT_0167 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0168 = IT_0014*IT_0167;
    const ccomplex_t IT_0169 = IT_0010*IT_0167;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0165 + IT_0166 + 
      -IT_0168 + -IT_0169);
    const ccomplex_t IT_0171 = (-0.5)*IT_0170;
    const ccomplex_t IT_0172 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0173 = IT_0013*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0176 = IT_0009*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = conj(N_d1)*e_em*m_tau*IT_0028*U_stau_10;
    const ccomplex_t IT_0179 = IT_0027*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*(IT_0174 + IT_0177 + 
      -IT_0180);
    const ccomplex_t IT_0182 = (-0.5)*IT_0181;
    const ccomplex_t IT_0183 = IT_0059*IT_0182;
    const ccomplex_t IT_0184 = (0 + _Complex_I*1)*IT_0183;
    const ccomplex_t IT_0185 = IT_0171*IT_0184;
    const ccomplex_t IT_0186 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0187 = IT_0014*IT_0186;
    const ccomplex_t IT_0188 = IT_0010*IT_0186;
    const ccomplex_t IT_0189 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0190 = IT_0010*IT_0189;
    const ccomplex_t IT_0191 = IT_0014*IT_0189;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0187 + IT_0188 + 
      -IT_0190 + -IT_0191);
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = (-0.5)*IT_0193;
    const ccomplex_t IT_0195 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0196 = IT_0013*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0199 = IT_0009*IT_0198;
    const ccomplex_t IT_0200 = 1.4142135623731*IT_0199;
    const ccomplex_t IT_0201 = conj(N_d2)*e_em*m_tau*IT_0028*U_stau_10;
    const ccomplex_t IT_0202 = IT_0027*IT_0201;
    const ccomplex_t IT_0203 = 1.4142135623731*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0197 + IT_0200 + 
      -IT_0203);
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0083*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = IT_0194*IT_0207;
    const ccomplex_t IT_0209 = (-0.5)*IT_0160 + (-0.25)*IT_0163 + (-0.5)
      *IT_0185 + (-0.5)*IT_0208;
    const ccomplex_t IT_0210 = IT_0021*IT_0162;
    const ccomplex_t IT_0211 = (-0.25)*IT_0210;
    const ccomplex_t IT_0212 = IT_0209 + -IT_0211;
    const ccomplex_t IT_0213 = conj(IT_0209) + -conj(IT_0211);
    const ccomplex_t IT_0214 = pow(m_Z, 2);
    const ccomplex_t IT_0215 = s_13*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0004;
    const ccomplex_t IT_0217 = (-2)*IT_0000;
    const ccomplex_t IT_0218 = IT_0215 + IT_0216 + IT_0217;
    const ccomplex_t IT_0219 = m_N_4*m_tau;
    const ccomplex_t IT_0220 = pow(s_14, 2);
    const ccomplex_t IT_0221 = IT_0001*IT_0220;
    const ccomplex_t IT_0222 = -IT_0111;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = IT_0219*IT_0223;
    const ccomplex_t IT_0225 = IT_0160 + 0.5*IT_0163 + IT_0185 + IT_0208;
    const ccomplex_t IT_0226 = -IT_0032;
    const ccomplex_t IT_0227 = IT_0038*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = IT_0021*IT_0228;
    const ccomplex_t IT_0230 = -IT_0062 + -IT_0086 + -IT_0106 + (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = 2*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0007
      *IT_0009 + -IT_0008*IT_0013)*U_stau_00*conj(U_stau_01) + IT_0008*IT_0013
      *U_stau_10*conj(U_stau_11));
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0235 = IT_0013*IT_0234;
    const ccomplex_t IT_0236 = 1.4142135623731*IT_0235;
    const ccomplex_t IT_0237 = N_d4*e_em*m_tau*IT_0028*U_stau_01;
    const ccomplex_t IT_0238 = IT_0027*IT_0237;
    const ccomplex_t IT_0239 = 1.4142135623731*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0236 + 0.5*IT_0239);
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = cpow(s_13 + (-0.5)*IT_0036 + (-0.5)*IT_0111 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0243 = IT_0241*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = IT_0233*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*2)*(IT_0024 + 1./2*IT_0031);
    const ccomplex_t IT_0247 = 0.5*IT_0246;
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = IT_0038*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = IT_0108*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*e_em*((-0.5)*(IT_0007
      *IT_0009 + -IT_0008*IT_0013)*U_stau_00*conj(U_stau_00) + IT_0008*IT_0013
      *U_stau_10*conj(U_stau_10));
    const ccomplex_t IT_0253 = (-2)*IT_0252;
    const ccomplex_t IT_0254 = cpow(s_13 + (-0.5)*IT_0036 + (-0.5)*IT_0111 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0255 = IT_0248*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = IT_0253*IT_0256;
    const ccomplex_t IT_0258 = (-0.5)*IT_0245 + (-0.5)*IT_0251 + (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = (-2)*IT_0258;
    const ccomplex_t IT_0260 = IT_0231 + IT_0259;
    const ccomplex_t IT_0261 = 2*conj(IT_0230);
    const ccomplex_t IT_0262 = (-2)*conj(IT_0258);
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = IT_0130*IT_0254;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = IT_0253*IT_0265;
    const ccomplex_t IT_0267 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0268 = IT_0013*IT_0267;
    const ccomplex_t IT_0269 = 1.4142135623731*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0271 = IT_0009*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = conj(N_d4)*e_em*m_tau*IT_0028*U_stau_11;
    const ccomplex_t IT_0274 = IT_0027*IT_0273;
    const ccomplex_t IT_0275 = 1.4142135623731*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*(IT_0269 + IT_0272 + 
      -IT_0275);
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = IT_0242*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = IT_0233*IT_0279;
    const ccomplex_t IT_0281 = 0.5*IT_0210 + 0.5*IT_0266 + 0.5*IT_0280;
    const ccomplex_t IT_0282 = 2*IT_0258;
    const ccomplex_t IT_0283 = (-2)*IT_0230;
    const ccomplex_t IT_0284 = IT_0282 + IT_0283;
    const ccomplex_t IT_0285 = 2*conj(IT_0258);
    const ccomplex_t IT_0286 = (-2)*conj(IT_0230);
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = s_13 + IT_0289;
    const ccomplex_t IT_0291 = s_13*IT_0290;
    const ccomplex_t IT_0292 = IT_0113*IT_0248;
    const ccomplex_t IT_0293 = IT_0112*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = 0.5*IT_0245 + 0.5*IT_0257 + IT_0294;
    const ccomplex_t IT_0296 = -IT_0137 + (-0.5)*IT_0266 + (-0.5)*IT_0280;
    const ccomplex_t IT_0297 = 2*IT_0225;
    const ccomplex_t IT_0298 = (-2)*IT_0281;
    const ccomplex_t IT_0299 = IT_0297 + IT_0298;
    const ccomplex_t IT_0300 = 2*conj(IT_0225);
    const ccomplex_t IT_0301 = (-2)*conj(IT_0281);
    const ccomplex_t IT_0302 = IT_0300 + IT_0301;
    const ccomplex_t IT_0303 = m_N_4*IT_0290;
    const ccomplex_t IT_0304 = m_N_2*IT_0194;
    const ccomplex_t IT_0305 = m_N_4*IT_0073;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0304 + IT_0306;
    const ccomplex_t IT_0308 = IT_0085*IT_0307;
    const ccomplex_t IT_0309 = m_N_4*IT_0108;
    const ccomplex_t IT_0310 = m_N_4*IT_0021;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = IT_0309 + IT_0311;
    const ccomplex_t IT_0313 = IT_0040*IT_0312;
    const ccomplex_t IT_0314 = m_N_1*IT_0171;
    const ccomplex_t IT_0315 = m_N_4*IT_0049;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = IT_0061*IT_0317;
    const ccomplex_t IT_0319 = m_N_3*IT_0146;
    const ccomplex_t IT_0320 = m_N_4*IT_0094;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0319 + IT_0321;
    const ccomplex_t IT_0323 = IT_0105*IT_0322;
    const ccomplex_t IT_0324 = 0.5*IT_0308 + 0.25*IT_0313 + 0.5*IT_0318 + 0.5
      *IT_0323;
    const ccomplex_t IT_0325 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0326 = IT_0113*IT_0130;
    const ccomplex_t IT_0327 = IT_0112*IT_0325*IT_0326;
    const ccomplex_t IT_0328 = m_tau*IT_0137;
    const ccomplex_t IT_0329 = 0.25*IT_0313 + (-0.5)*IT_0327 + 0.5*IT_0328;
    const ccomplex_t IT_0330 = conj(IT_0324) + conj(IT_0329);
    const ccomplex_t IT_0331 = IT_0324 + IT_0329;
    const ccomplex_t IT_0332 = IT_0162*IT_0312;
    const ccomplex_t IT_0333 = m_N_1*IT_0049;
    const ccomplex_t IT_0334 = m_N_4*IT_0171;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0333 + IT_0335;
    const ccomplex_t IT_0337 = IT_0184*IT_0336;
    const ccomplex_t IT_0338 = m_N_4*IT_0194;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = m_N_2*IT_0073;
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = IT_0207*IT_0341;
    const ccomplex_t IT_0343 = m_N_4*IT_0146;
    const ccomplex_t IT_0344 = -IT_0343;
    const ccomplex_t IT_0345 = m_N_3*IT_0094;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = IT_0159*IT_0346;
    const ccomplex_t IT_0348 = 0.25*IT_0332 + (-0.5)*IT_0337 + (-0.5)*IT_0342 
      + (-0.5)*IT_0347;
    const ccomplex_t IT_0349 = IT_0035*IT_0134;
    const ccomplex_t IT_0350 = IT_0112*IT_0325*IT_0349;
    const ccomplex_t IT_0351 = m_tau*IT_0116;
    const ccomplex_t IT_0352 = 0.25*IT_0332 + 0.5*IT_0350 + (-0.5)*IT_0351;
    const ccomplex_t IT_0353 = conj(IT_0348) + conj(IT_0352);
    const ccomplex_t IT_0354 = IT_0348 + IT_0352;
    const ccomplex_t IT_0355 = s_13*s_14;
    const ccomplex_t IT_0356 = s_34*IT_0111;
    const ccomplex_t IT_0357 = -IT_0356;
    const ccomplex_t IT_0358 = IT_0355 + IT_0357;
    const ccomplex_t IT_0359 = conj(IT_0107) + conj(IT_0117);
    const ccomplex_t IT_0360 = conj(IT_0138) + conj(IT_0209);
    const ccomplex_t IT_0361 = IT_0138 + IT_0209;
    const ccomplex_t IT_0362 = 2*IT_0281;
    const ccomplex_t IT_0363 = (-2)*IT_0225;
    const ccomplex_t IT_0364 = IT_0362 + IT_0363;
    const ccomplex_t IT_0365 = 2*conj(IT_0281);
    const ccomplex_t IT_0366 = (-2)*conj(IT_0225);
    const ccomplex_t IT_0367 = IT_0365 + IT_0366;
    const ccomplex_t IT_0368 = s_14*IT_0214;
    const ccomplex_t IT_0369 = IT_0001*IT_0368;
    const ccomplex_t IT_0370 = (-0.25)*IT_0369;
    const ccomplex_t IT_0371 = s_14 + IT_0370;
    const ccomplex_t IT_0372 = m_tau*IT_0371;
    const ccomplex_t IT_0373 = (-4)*IT_0372;
    const ccomplex_t IT_0374 = 2*IT_0324;
    const ccomplex_t IT_0375 = 2*IT_0329;
    const ccomplex_t IT_0376 = IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = 2*conj(IT_0324);
    const ccomplex_t IT_0378 = 2*conj(IT_0329);
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = 2*IT_0348;
    const ccomplex_t IT_0381 = 2*IT_0352;
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = 2*conj(IT_0348);
    const ccomplex_t IT_0384 = 2*conj(IT_0352);
    const ccomplex_t IT_0385 = IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = s_34*IT_0214;
    const ccomplex_t IT_0387 = IT_0001*IT_0386;
    const ccomplex_t IT_0388 = 0.5*IT_0387;
    const ccomplex_t IT_0389 = s_34 + IT_0388;
    const ccomplex_t IT_0390 = m_N_4*IT_0389;
    const ccomplex_t IT_0391 = 2*IT_0390;
    const ccomplex_t IT_0392 = 0.5*IT_0369;
    const ccomplex_t IT_0393 = s_14 + IT_0392;
    const ccomplex_t IT_0394 = m_tau*IT_0393;
    const ccomplex_t IT_0395 = 2*IT_0394;
    const ccomplex_t IT_0396 = (-2)*IT_0324;
    const ccomplex_t IT_0397 = (-2)*IT_0329;
    const ccomplex_t IT_0398 = IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = (-2)*conj(IT_0324);
    const ccomplex_t IT_0400 = (-2)*conj(IT_0329);
    const ccomplex_t IT_0401 = IT_0399 + IT_0400;
    const ccomplex_t IT_0402 = (-2)*IT_0348;
    const ccomplex_t IT_0403 = (-2)*IT_0352;
    const ccomplex_t IT_0404 = IT_0402 + IT_0403;
    const ccomplex_t IT_0405 = (-2)*conj(IT_0348);
    const ccomplex_t IT_0406 = (-2)*conj(IT_0352);
    const ccomplex_t IT_0407 = IT_0405 + IT_0406;
    const ccomplex_t IT_0408 = (-0.25)*IT_0387;
    const ccomplex_t IT_0409 = s_34 + IT_0408;
    const ccomplex_t IT_0410 = m_N_4*IT_0409;
    const ccomplex_t IT_0411 = (-4)*IT_0410;
    const ccomplex_t IT_0412 = pow(s_34, 2);
    const ccomplex_t IT_0413 = IT_0001*IT_0412;
    const ccomplex_t IT_0414 = -IT_0036;
    const ccomplex_t IT_0415 = IT_0413 + IT_0414;
    const ccomplex_t IT_0416 = IT_0219*IT_0415;
    const ccomplex_t IT_0417 = 2*IT_0295;
    const ccomplex_t IT_0418 = 2*conj(IT_0295);
    const ccomplex_t IT_0419 = s_13*IT_0415;
    const ccomplex_t IT_0420 = 2*IT_0296;
    const ccomplex_t IT_0421 = m_N_4*IT_0415;
    const ccomplex_t IT_0422 = 2*conj(IT_0296);
    const ccomplex_t IT_0423 = IT_0215 + IT_0217;
    const ccomplex_t IT_0424 = IT_0001*IT_0423;
    const ccomplex_t IT_0425 = (-0.5)*IT_0424;
    const ccomplex_t IT_0426 = s_13 + IT_0425;
    const ccomplex_t IT_0427 = 4*IT_0348;
    const ccomplex_t IT_0428 = (-6)*conj(IT_0324);
    const ccomplex_t IT_0429 = s_13*IT_0223;
    const ccomplex_t IT_0430 = IT_0219*IT_0290;
    const ccomplex_t IT_0431 = (-2)*IT_0295;
    const ccomplex_t IT_0432 = (-2)*conj(IT_0295);
    const ccomplex_t IT_0433 = s_13*s_34;
    const ccomplex_t IT_0434 = s_14*IT_0036;
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = IT_0433 + IT_0435;
    const ccomplex_t IT_0437 = (-2)*IT_0296;
    const ccomplex_t IT_0438 = (-2)*conj(IT_0296);
    const ccomplex_t IT_0439 = m_tau*IT_0223;
    const ccomplex_t IT_0440 = m_tau*IT_0290;
    const ccomplex_t IT_0441 = IT_0214*IT_0219;
    const ccomplex_t IT_0442 = (-3)*IT_0441;
    const ccomplex_t IT_0443 = 2*IT_0110;
    const ccomplex_t IT_0444 = (-2)*IT_0107;
    const ccomplex_t IT_0445 = IT_0443 + IT_0444;
    const ccomplex_t IT_0446 = 2*conj(IT_0110);
    const ccomplex_t IT_0447 = (-2)*conj(IT_0107);
    const ccomplex_t IT_0448 = IT_0446 + IT_0447;
    const ccomplex_t IT_0449 = 2*IT_0117;
    const ccomplex_t IT_0450 = 2*conj(IT_0117);
    const ccomplex_t IT_0451 = (-2)*IT_0209;
    const ccomplex_t IT_0452 = (-2)*conj(IT_0209);
    const ccomplex_t IT_0453 = conj(IT_0107)*IT_0209 + IT_0107*conj(IT_0209) +
       0.5*conj(IT_0211)*IT_0445 + 0.5*IT_0211*IT_0448 + 0.5*conj(IT_0138)
      *IT_0449 + 0.5*IT_0138*IT_0450 + 0.5*conj(IT_0110)*IT_0451 + 0.5*IT_0110
      *IT_0452;
    const ccomplex_t IT_0454 = 3*IT_0441;
    const ccomplex_t IT_0455 = 8*IT_0006*((IT_0107 + -IT_0110)*conj(IT_0117) +
       IT_0117*IT_0119 + conj(IT_0138)*IT_0212 + IT_0138*IT_0213) + (-4)*
      (IT_0117*conj(IT_0117) + (IT_0107 + -IT_0110)*IT_0119 + IT_0138*conj
      (IT_0138) + IT_0212*IT_0213)*IT_0218 + IT_0224*(conj(IT_0225)*IT_0260 +
       IT_0225*IT_0263 + conj(IT_0281)*IT_0284 + IT_0281*IT_0287) + IT_0291*
      (IT_0263*IT_0295 + IT_0260*conj(IT_0295) + conj(IT_0296)*IT_0299 + IT_0296
      *IT_0302) + IT_0303*(IT_0260*IT_0330 + IT_0263*IT_0331 + IT_0299*IT_0353 +
       IT_0302*IT_0354) + IT_0358*((IT_0107 + IT_0117)*IT_0263 + conj(IT_0110)
      *IT_0284 + IT_0110*IT_0287 + IT_0260*IT_0359 + IT_0299*IT_0360 + IT_0302
      *IT_0361 + conj(IT_0211)*IT_0364 + IT_0211*IT_0367) + IT_0373*(conj
      (IT_0138)*IT_0376 + IT_0138*IT_0379 + conj(IT_0117)*IT_0382 + IT_0117
      *IT_0385) + (conj(IT_0117)*IT_0376 + IT_0117*IT_0379 + conj(IT_0138)
      *IT_0382 + IT_0138*IT_0385)*IT_0391 + IT_0395*(conj(IT_0211)*IT_0376 +
       IT_0211*IT_0379 + conj(IT_0110)*IT_0382 + IT_0110*IT_0385 + conj(IT_0209)
      *IT_0398 + IT_0209*IT_0401 + conj(IT_0107)*IT_0404 + IT_0107*IT_0407) + 
      (conj(IT_0110)*IT_0376 + IT_0110*IT_0379 + conj(IT_0211)*IT_0382 + IT_0211
      *IT_0385 + conj(IT_0107)*IT_0398 + IT_0107*IT_0401 + conj(IT_0209)*IT_0404
       + IT_0209*IT_0407)*IT_0411 + IT_0416*(conj(IT_0296)*IT_0417 + IT_0296
      *IT_0418) + IT_0419*(conj(IT_0295)*IT_0417 + conj(IT_0296)*IT_0420) +
       IT_0421*(IT_0330*IT_0417 + IT_0331*IT_0418 + IT_0353*IT_0420 + IT_0354
      *IT_0422) + 4*IT_0426*(IT_0330*IT_0331 + IT_0353*(IT_0352 + 0.25*IT_0427))
       + (-6)*IT_0219*(IT_0331*IT_0353 + IT_0354*(conj(IT_0329) + (
      -0.166666666666667)*IT_0428)) + (conj(IT_0230)*IT_0231 + conj(IT_0258)
      *IT_0284 + IT_0258*IT_0286 + conj(IT_0225)*IT_0297 + conj(IT_0281)*IT_0364
       + IT_0281*IT_0366)*IT_0429 + IT_0430*(IT_0263*IT_0296 + IT_0260*conj
      (IT_0296) + conj(IT_0225)*IT_0417 + IT_0225*IT_0418 + conj(IT_0281)
      *IT_0431 + IT_0281*IT_0432) + IT_0436*(IT_0359*IT_0417 + (IT_0107 +
       IT_0117)*IT_0418 + IT_0360*IT_0420 + IT_0361*IT_0422 + conj(IT_0110)
      *IT_0431 + IT_0110*IT_0432 + conj(IT_0211)*IT_0437 + IT_0211*IT_0438) + 
      (IT_0284*IT_0353 + IT_0287*IT_0354 + conj(IT_0329)*IT_0364 + IT_0329
      *IT_0367 + conj(IT_0281)*IT_0374 + IT_0281*IT_0377 + conj(IT_0225)*IT_0396
       + IT_0225*IT_0399)*IT_0439 + (conj(IT_0296)*IT_0396 + IT_0296*IT_0399 +
       IT_0353*IT_0431 + IT_0354*IT_0432 + conj(IT_0329)*IT_0437 + IT_0329
      *IT_0438)*IT_0440 + 2*IT_0442*IT_0453 + (conj(IT_0138)*IT_0445 + IT_0138
      *IT_0448 + conj(IT_0211)*IT_0449 + IT_0211*IT_0450 + conj(IT_0117)*IT_0451
       + IT_0117*IT_0452)*IT_0454;
    return create_ccomplex_return(IT_0455);
}

