#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_3_stau_2_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_3_stau_2_to_tau_Z(
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
    const ccomplex_t IT_0000 = pow(m_Z, 2);
    const ccomplex_t IT_0001 = s_13*IT_0000;
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(m_Z, 4);
    const ccomplex_t IT_0004 = s_13*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = s_14*s_34;
    const ccomplex_t IT_0008 = (-2)*IT_0007;
    const ccomplex_t IT_0009 = IT_0001 + IT_0006 + IT_0008;
    const ccomplex_t IT_0010 = pow(m_N_3, 2);
    const ccomplex_t IT_0011 = pow(m_tau, 2);
    const ccomplex_t IT_0012 = pow(m_stau_2, 2);
    const ccomplex_t IT_0013 = cpow(s_12 + 0.5*IT_0010 + (-0.5)*IT_0011 + 0.5
      *IT_0012 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = e_em*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0016, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = e_em*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0018 + -IT_0021);
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0028 = IT_0019*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = cos(beta);
    const ccomplex_t IT_0031 = cpow(IT_0030, -1);
    const ccomplex_t IT_0032 = IT_0019*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -1);
    const ccomplex_t IT_0034 = conj(N_d3)*e_em*m_tau*IT_0033*U_stau_11;
    const ccomplex_t IT_0035 = IT_0032*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0026 + IT_0029 + 
      -IT_0036);
    const ccomplex_t IT_0038 = (-0.5)*IT_0037;
    const ccomplex_t IT_0039 = IT_0023*IT_0038;
    const ccomplex_t IT_0040 = IT_0013*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0044 = IT_0015*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = N_d3*e_em*m_tau*IT_0033*U_stau_01;
    const ccomplex_t IT_0047 = IT_0032*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0045 + 0.5*IT_0048);
    const ccomplex_t IT_0050 = -IT_0049;
    const ccomplex_t IT_0051 = (0 + _Complex_I*1)*e_em*IT_0015*IT_0016;
    const ccomplex_t IT_0052 = IT_0050*IT_0051;
    const ccomplex_t IT_0053 = IT_0013*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = 0.5*IT_0054;
    const ccomplex_t IT_0056 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0057 = IT_0017*IT_0056;
    const ccomplex_t IT_0058 = IT_0020*IT_0056;
    const ccomplex_t IT_0059 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0060 = IT_0017*IT_0059;
    const ccomplex_t IT_0061 = IT_0020*IT_0059;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0057 + IT_0058 + 
      -IT_0060 + -IT_0061);
    const ccomplex_t IT_0063 = (-0.5)*IT_0062;
    const ccomplex_t IT_0064 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0065 = IT_0015*IT_0064;
    const ccomplex_t IT_0066 = 1.4142135623731*IT_0065;
    const ccomplex_t IT_0067 = N_d1*e_em*m_tau*IT_0033*U_stau_01;
    const ccomplex_t IT_0068 = IT_0032*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + 0.5*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = cpow(s_23 + (-0.5)*IT_0011 + (-0.5)*IT_0012 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0063*IT_0074;
    const ccomplex_t IT_0076 = N_u3*conj(N_u3)*e_em;
    const ccomplex_t IT_0077 = IT_0017*IT_0076;
    const ccomplex_t IT_0078 = IT_0020*IT_0076;
    const ccomplex_t IT_0079 = N_d3*conj(N_d3)*e_em;
    const ccomplex_t IT_0080 = IT_0017*IT_0079;
    const ccomplex_t IT_0081 = IT_0020*IT_0079;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0077 + IT_0078 + 
      -IT_0080 + -IT_0081);
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = cpow(s_23 + (-0.5)*IT_0011 + (-0.5)*IT_0012 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0050*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0083*IT_0086;
    const ccomplex_t IT_0088 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0089 = IT_0017*IT_0088;
    const ccomplex_t IT_0090 = IT_0020*IT_0088;
    const ccomplex_t IT_0091 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0092 = IT_0017*IT_0091;
    const ccomplex_t IT_0093 = IT_0020*IT_0091;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*(IT_0089 + IT_0090 + 
      -IT_0092 + -IT_0093);
    const ccomplex_t IT_0095 = (-0.5)*IT_0094;
    const ccomplex_t IT_0096 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0097 = IT_0015*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = N_d2*e_em*m_tau*IT_0033*U_stau_01;
    const ccomplex_t IT_0100 = IT_0032*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + 0.5*IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = cpow(s_23 + (-0.5)*IT_0011 + (-0.5)*IT_0012 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0095*IT_0106;
    const ccomplex_t IT_0108 = 0.5*IT_0075 + 0.25*IT_0087 + 0.5*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0082;
    const ccomplex_t IT_0110 = IT_0086*IT_0109;
    const ccomplex_t IT_0111 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0112 = IT_0017*IT_0111;
    const ccomplex_t IT_0113 = IT_0020*IT_0111;
    const ccomplex_t IT_0114 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0115 = IT_0017*IT_0114;
    const ccomplex_t IT_0116 = IT_0020*IT_0114;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*(IT_0112 + IT_0113 + 
      -IT_0115 + -IT_0116);
    const ccomplex_t IT_0118 = (-0.5)*IT_0117;
    const ccomplex_t IT_0119 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0120 = IT_0015*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = N_d4*e_em*m_tau*IT_0033*U_stau_01;
    const ccomplex_t IT_0123 = IT_0032*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = cpow(s_23 + (-0.5)*IT_0011 + (-0.5)*IT_0012 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0128 = IT_0126*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*IT_0128;
    const ccomplex_t IT_0130 = IT_0118*IT_0129;
    const ccomplex_t IT_0131 = 0.25*IT_0110 + 0.5*IT_0130;
    const ccomplex_t IT_0132 = IT_0108 + -IT_0131;
    const ccomplex_t IT_0133 = conj(IT_0108) + -conj(IT_0131);
    const ccomplex_t IT_0134 = IT_0038*IT_0084;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = IT_0109*IT_0135;
    const ccomplex_t IT_0137 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0138 = IT_0017*IT_0137;
    const ccomplex_t IT_0139 = IT_0020*IT_0137;
    const ccomplex_t IT_0140 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0141 = IT_0017*IT_0140;
    const ccomplex_t IT_0142 = IT_0020*IT_0140;
    const ccomplex_t IT_0143 = (0 + _Complex_I*1)*(IT_0138 + IT_0139 + 
      -IT_0141 + -IT_0142);
    const ccomplex_t IT_0144 = (-0.5)*IT_0143;
    const ccomplex_t IT_0145 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0146 = IT_0015*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0149 = IT_0019*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = conj(N_d1)*e_em*m_tau*IT_0033*U_stau_11;
    const ccomplex_t IT_0152 = IT_0032*IT_0151;
    const ccomplex_t IT_0153 = 1.4142135623731*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0147 + IT_0150 + 
      -IT_0153);
    const ccomplex_t IT_0155 = (-0.5)*IT_0154;
    const ccomplex_t IT_0156 = IT_0072*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = IT_0144*IT_0157;
    const ccomplex_t IT_0159 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0160 = IT_0017*IT_0159;
    const ccomplex_t IT_0161 = IT_0020*IT_0159;
    const ccomplex_t IT_0162 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0163 = IT_0017*IT_0162;
    const ccomplex_t IT_0164 = IT_0020*IT_0162;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0160 + IT_0161 + 
      -IT_0163 + -IT_0164);
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0168 = IT_0015*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0171 = IT_0019*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_d2)*e_em*m_tau*IT_0033*U_stau_11;
    const ccomplex_t IT_0174 = IT_0032*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0169 + IT_0172 + 
      -IT_0175);
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0104*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = IT_0166*IT_0179;
    const ccomplex_t IT_0181 = (-0.25)*IT_0136 + (-0.5)*IT_0158 + (-0.5)
      *IT_0180;
    const ccomplex_t IT_0182 = IT_0083*IT_0135;
    const ccomplex_t IT_0183 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0184 = IT_0017*IT_0183;
    const ccomplex_t IT_0185 = IT_0020*IT_0183;
    const ccomplex_t IT_0186 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0187 = IT_0017*IT_0186;
    const ccomplex_t IT_0188 = IT_0020*IT_0186;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0184 + IT_0185 + 
      -IT_0187 + -IT_0188);
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0192 = IT_0015*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0195 = IT_0019*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = conj(N_d4)*e_em*m_tau*IT_0033*U_stau_11;
    const ccomplex_t IT_0198 = IT_0032*IT_0197;
    const ccomplex_t IT_0199 = 1.4142135623731*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*(IT_0193 + IT_0196 + 
      -IT_0199);
    const ccomplex_t IT_0201 = (-0.5)*IT_0200;
    const ccomplex_t IT_0202 = IT_0127*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = IT_0190*IT_0203;
    const ccomplex_t IT_0205 = (-0.25)*IT_0182 + (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = IT_0181 + -IT_0205;
    const ccomplex_t IT_0207 = conj(IT_0181) + -conj(IT_0205);
    const ccomplex_t IT_0208 = m_N_3*m_tau;
    const ccomplex_t IT_0209 = pow(s_14, 2);
    const ccomplex_t IT_0210 = IT_0002*IT_0209;
    const ccomplex_t IT_0211 = -IT_0010;
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = IT_0208*IT_0212;
    const ccomplex_t IT_0214 = 0.5*IT_0136 + IT_0158 + IT_0180;
    const ccomplex_t IT_0215 = -IT_0075 + (-0.5)*IT_0087 + -IT_0107;
    const ccomplex_t IT_0216 = 2*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0015*IT_0016 +
       -IT_0014*IT_0019)*U_stau_01*conj(U_stau_01) + IT_0015*IT_0016*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0218 = (-2)*IT_0217;
    const ccomplex_t IT_0219 = cpow(s_13 + (-0.5)*IT_0010 + (-0.5)*IT_0011 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0220 = IT_0050*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = IT_0218*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0015*IT_0016 +
       -IT_0014*IT_0019)*conj(U_stau_00)*U_stau_01 + IT_0015*IT_0016*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0224 = (-2)*IT_0223;
    const ccomplex_t IT_0225 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0226 = IT_0015*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = N_d3*e_em*m_tau*IT_0033*U_stau_00;
    const ccomplex_t IT_0229 = IT_0032*IT_0228;
    const ccomplex_t IT_0230 = 1.4142135623731*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*(IT_0227 + 0.5*IT_0230);
    const ccomplex_t IT_0232 = -IT_0231;
    const ccomplex_t IT_0233 = cpow(s_13 + (-0.5)*IT_0010 + (-0.5)*IT_0011 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0234 = IT_0232*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = IT_0224*IT_0235;
    const ccomplex_t IT_0237 = (-0.5)*IT_0110 + -IT_0130 + (-0.5)*IT_0222 + (
      -0.5)*IT_0236;
    const ccomplex_t IT_0238 = (-2)*IT_0237;
    const ccomplex_t IT_0239 = IT_0216 + IT_0238;
    const ccomplex_t IT_0240 = 2*conj(IT_0215);
    const ccomplex_t IT_0241 = (-2)*conj(IT_0237);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = IT_0038*IT_0219;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*IT_0243;
    const ccomplex_t IT_0245 = IT_0218*IT_0244;
    const ccomplex_t IT_0246 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0247 = IT_0015*IT_0246;
    const ccomplex_t IT_0248 = 1.4142135623731*IT_0247;
    const ccomplex_t IT_0249 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0250 = IT_0019*IT_0249;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = conj(N_d3)*e_em*m_tau*IT_0033*U_stau_10;
    const ccomplex_t IT_0253 = IT_0032*IT_0252;
    const ccomplex_t IT_0254 = 1.4142135623731*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*(IT_0248 + IT_0251 + 
      -IT_0254);
    const ccomplex_t IT_0256 = (-0.5)*IT_0255;
    const ccomplex_t IT_0257 = IT_0233*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = IT_0224*IT_0258;
    const ccomplex_t IT_0260 = 0.5*IT_0182 + IT_0204 + 0.5*IT_0245 + 0.5
      *IT_0259;
    const ccomplex_t IT_0261 = 2*IT_0237;
    const ccomplex_t IT_0262 = (-2)*IT_0215;
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 2*conj(IT_0237);
    const ccomplex_t IT_0265 = (-2)*conj(IT_0215);
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = s_14*s_34*IT_0002;
    const ccomplex_t IT_0268 = -IT_0267;
    const ccomplex_t IT_0269 = s_13 + IT_0268;
    const ccomplex_t IT_0270 = s_13*IT_0269;
    const ccomplex_t IT_0271 = IT_0054 + 0.5*IT_0222 + 0.5*IT_0236;
    const ccomplex_t IT_0272 = -IT_0041 + (-0.5)*IT_0245 + (-0.5)*IT_0259;
    const ccomplex_t IT_0273 = 2*IT_0214;
    const ccomplex_t IT_0274 = (-2)*IT_0260;
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 2*conj(IT_0214);
    const ccomplex_t IT_0277 = (-2)*conj(IT_0260);
    const ccomplex_t IT_0278 = IT_0276 + IT_0277;
    const ccomplex_t IT_0279 = m_N_3*IT_0269;
    const ccomplex_t IT_0280 = m_N_3*IT_0109;
    const ccomplex_t IT_0281 = m_N_3*IT_0083;
    const ccomplex_t IT_0282 = -IT_0281;
    const ccomplex_t IT_0283 = IT_0280 + IT_0282;
    const ccomplex_t IT_0284 = IT_0086*IT_0283;
    const ccomplex_t IT_0285 = m_N_1*IT_0144;
    const ccomplex_t IT_0286 = m_N_3*IT_0063;
    const ccomplex_t IT_0287 = -IT_0286;
    const ccomplex_t IT_0288 = IT_0285 + IT_0287;
    const ccomplex_t IT_0289 = IT_0074*IT_0288;
    const ccomplex_t IT_0290 = m_N_2*IT_0166;
    const ccomplex_t IT_0291 = m_N_3*IT_0095;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = IT_0290 + IT_0292;
    const ccomplex_t IT_0294 = IT_0106*IT_0293;
    const ccomplex_t IT_0295 = 0.25*IT_0284 + 0.5*IT_0289 + 0.5*IT_0294;
    const ccomplex_t IT_0296 = m_N_4*IT_0190;
    const ccomplex_t IT_0297 = m_N_3*IT_0118;
    const ccomplex_t IT_0298 = -IT_0297;
    const ccomplex_t IT_0299 = IT_0296 + IT_0298;
    const ccomplex_t IT_0300 = IT_0129*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0302 = IT_0038*IT_0051;
    const ccomplex_t IT_0303 = IT_0013*IT_0301*IT_0302;
    const ccomplex_t IT_0304 = m_tau*IT_0041;
    const ccomplex_t IT_0305 = 0.25*IT_0284 + (-0.5)*IT_0300 + (-0.5)*IT_0303 
      + 0.5*IT_0304;
    const ccomplex_t IT_0306 = conj(IT_0295) + conj(IT_0305);
    const ccomplex_t IT_0307 = IT_0295 + IT_0305;
    const ccomplex_t IT_0308 = IT_0135*IT_0283;
    const ccomplex_t IT_0309 = m_N_1*IT_0063;
    const ccomplex_t IT_0310 = m_N_3*IT_0144;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = IT_0309 + IT_0311;
    const ccomplex_t IT_0313 = IT_0157*IT_0312;
    const ccomplex_t IT_0314 = m_N_2*IT_0095;
    const ccomplex_t IT_0315 = m_N_3*IT_0166;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = IT_0179*IT_0317;
    const ccomplex_t IT_0319 = 0.25*IT_0308 + (-0.5)*IT_0313 + (-0.5)*IT_0318;
    const ccomplex_t IT_0320 = IT_0023*IT_0050;
    const ccomplex_t IT_0321 = IT_0013*IT_0301*IT_0320;
    const ccomplex_t IT_0322 = m_N_4*IT_0118;
    const ccomplex_t IT_0323 = m_N_3*IT_0190;
    const ccomplex_t IT_0324 = -IT_0323;
    const ccomplex_t IT_0325 = IT_0322 + IT_0324;
    const ccomplex_t IT_0326 = IT_0203*IT_0325;
    const ccomplex_t IT_0327 = m_tau*IT_0054;
    const ccomplex_t IT_0328 = 0.25*IT_0308 + 0.5*IT_0321 + 0.5*IT_0326 + (
      -0.5)*IT_0327;
    const ccomplex_t IT_0329 = conj(IT_0319) + conj(IT_0328);
    const ccomplex_t IT_0330 = IT_0319 + IT_0328;
    const ccomplex_t IT_0331 = s_13*s_14;
    const ccomplex_t IT_0332 = s_34*IT_0010;
    const ccomplex_t IT_0333 = -IT_0332;
    const ccomplex_t IT_0334 = IT_0331 + IT_0333;
    const ccomplex_t IT_0335 = conj(IT_0042) + conj(IT_0181);
    const ccomplex_t IT_0336 = 2*IT_0260;
    const ccomplex_t IT_0337 = (-2)*IT_0214;
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = 2*conj(IT_0260);
    const ccomplex_t IT_0340 = (-2)*conj(IT_0214);
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = s_14*IT_0000;
    const ccomplex_t IT_0343 = IT_0002*IT_0342;
    const ccomplex_t IT_0344 = 0.5*IT_0343;
    const ccomplex_t IT_0345 = s_14 + IT_0344;
    const ccomplex_t IT_0346 = m_tau*IT_0345;
    const ccomplex_t IT_0347 = 2*IT_0346;
    const ccomplex_t IT_0348 = 2*IT_0295;
    const ccomplex_t IT_0349 = 2*IT_0305;
    const ccomplex_t IT_0350 = IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = 2*conj(IT_0295);
    const ccomplex_t IT_0352 = 2*conj(IT_0305);
    const ccomplex_t IT_0353 = IT_0351 + IT_0352;
    const ccomplex_t IT_0354 = 2*IT_0319;
    const ccomplex_t IT_0355 = 2*IT_0328;
    const ccomplex_t IT_0356 = IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = 2*conj(IT_0319);
    const ccomplex_t IT_0358 = 2*conj(IT_0328);
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = (-2)*IT_0295;
    const ccomplex_t IT_0361 = (-2)*IT_0305;
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = (-2)*conj(IT_0295);
    const ccomplex_t IT_0364 = (-2)*conj(IT_0305);
    const ccomplex_t IT_0365 = IT_0363 + IT_0364;
    const ccomplex_t IT_0366 = (-2)*IT_0319;
    const ccomplex_t IT_0367 = (-2)*IT_0328;
    const ccomplex_t IT_0368 = IT_0366 + IT_0367;
    const ccomplex_t IT_0369 = (-2)*conj(IT_0319);
    const ccomplex_t IT_0370 = (-2)*conj(IT_0328);
    const ccomplex_t IT_0371 = IT_0369 + IT_0370;
    const ccomplex_t IT_0372 = s_34*IT_0000;
    const ccomplex_t IT_0373 = IT_0002*IT_0372;
    const ccomplex_t IT_0374 = (-0.25)*IT_0373;
    const ccomplex_t IT_0375 = s_34 + IT_0374;
    const ccomplex_t IT_0376 = m_N_3*IT_0375;
    const ccomplex_t IT_0377 = (-4)*IT_0376;
    const ccomplex_t IT_0378 = pow(s_34, 2);
    const ccomplex_t IT_0379 = IT_0002*IT_0378;
    const ccomplex_t IT_0380 = -IT_0011;
    const ccomplex_t IT_0381 = IT_0379 + IT_0380;
    const ccomplex_t IT_0382 = IT_0208*IT_0381;
    const ccomplex_t IT_0383 = 2*IT_0271;
    const ccomplex_t IT_0384 = 2*conj(IT_0271);
    const ccomplex_t IT_0385 = s_13*IT_0381;
    const ccomplex_t IT_0386 = 2*IT_0272;
    const ccomplex_t IT_0387 = m_N_3*IT_0381;
    const ccomplex_t IT_0388 = 2*conj(IT_0272);
    const ccomplex_t IT_0389 = IT_0000*IT_0208;
    const ccomplex_t IT_0390 = 3*IT_0389;
    const ccomplex_t IT_0391 = 2*IT_0055;
    const ccomplex_t IT_0392 = 2*conj(IT_0055);
    const ccomplex_t IT_0393 = IT_0001 + IT_0008;
    const ccomplex_t IT_0394 = IT_0002*IT_0393;
    const ccomplex_t IT_0395 = (-0.5)*IT_0394;
    const ccomplex_t IT_0396 = s_13 + IT_0395;
    const ccomplex_t IT_0397 = 4*IT_0319;
    const ccomplex_t IT_0398 = (-6)*conj(IT_0295);
    const ccomplex_t IT_0399 = s_13*IT_0212;
    const ccomplex_t IT_0400 = IT_0208*IT_0269;
    const ccomplex_t IT_0401 = (-2)*IT_0271;
    const ccomplex_t IT_0402 = (-2)*conj(IT_0271);
    const ccomplex_t IT_0403 = s_13*s_34;
    const ccomplex_t IT_0404 = s_14*IT_0011;
    const ccomplex_t IT_0405 = -IT_0404;
    const ccomplex_t IT_0406 = IT_0403 + IT_0405;
    const ccomplex_t IT_0407 = (-2)*IT_0272;
    const ccomplex_t IT_0408 = (-2)*conj(IT_0272);
    const ccomplex_t IT_0409 = m_tau*IT_0212;
    const ccomplex_t IT_0410 = m_tau*IT_0269;
    const ccomplex_t IT_0411 = (-3)*IT_0389;
    const ccomplex_t IT_0412 = 2*IT_0131;
    const ccomplex_t IT_0413 = (-2)*IT_0108;
    const ccomplex_t IT_0414 = IT_0412 + IT_0413;
    const ccomplex_t IT_0415 = 2*conj(IT_0131);
    const ccomplex_t IT_0416 = (-2)*conj(IT_0108);
    const ccomplex_t IT_0417 = IT_0415 + IT_0416;
    const ccomplex_t IT_0418 = (-2)*IT_0181;
    const ccomplex_t IT_0419 = (-2)*conj(IT_0181);
    const ccomplex_t IT_0420 = conj(IT_0108)*IT_0181 + IT_0108*conj(IT_0181) +
       0.5*conj(IT_0042)*IT_0391 + 0.5*IT_0042*IT_0392 + 0.5*conj(IT_0205)
      *IT_0414 + 0.5*IT_0205*IT_0417 + 0.5*conj(IT_0131)*IT_0418 + 0.5*IT_0131
      *IT_0419;
    const ccomplex_t IT_0421 = (-0.25)*IT_0005;
    const ccomplex_t IT_0422 = IT_0007 + IT_0421;
    const ccomplex_t IT_0423 = 0.125*IT_0390;
    const ccomplex_t IT_0424 = 0.5*IT_0373;
    const ccomplex_t IT_0425 = s_34 + IT_0424;
    const ccomplex_t IT_0426 = m_N_3*IT_0425;
    const ccomplex_t IT_0427 = 2*IT_0426;
    const ccomplex_t IT_0428 = 0.125*IT_0427;
    const ccomplex_t IT_0429 = (-0.25)*IT_0343;
    const ccomplex_t IT_0430 = s_14 + IT_0429;
    const ccomplex_t IT_0431 = m_tau*IT_0430;
    const ccomplex_t IT_0432 = (-4)*IT_0431;
    const ccomplex_t IT_0433 = 0.125*IT_0432;
    const ccomplex_t IT_0434 = (-4)*IT_0009*(IT_0042*conj(IT_0042) + IT_0055
      *conj(IT_0055) + IT_0132*IT_0133 + IT_0206*IT_0207) + IT_0213*(conj
      (IT_0214)*IT_0239 + IT_0214*IT_0242 + conj(IT_0260)*IT_0263 + IT_0260
      *IT_0266) + IT_0270*(IT_0242*IT_0271 + IT_0239*conj(IT_0271) + conj
      (IT_0272)*IT_0275 + IT_0272*IT_0278) + IT_0279*(IT_0239*IT_0306 + IT_0242
      *IT_0307 + IT_0275*IT_0329 + IT_0278*IT_0330) + IT_0334*((conj(IT_0055) +
       conj(IT_0108))*IT_0239 + (IT_0055 + IT_0108)*IT_0242 + conj(IT_0131)
      *IT_0263 + IT_0131*IT_0266 + (IT_0042 + IT_0181)*IT_0278 + IT_0275*IT_0335
       + conj(IT_0205)*IT_0338 + IT_0205*IT_0341) + IT_0347*(conj(IT_0205)
      *IT_0350 + IT_0205*IT_0353 + conj(IT_0131)*IT_0356 + IT_0131*IT_0359 +
       conj(IT_0181)*IT_0362 + IT_0181*IT_0365 + conj(IT_0108)*IT_0368 + IT_0108
      *IT_0371) + (conj(IT_0131)*IT_0350 + IT_0131*IT_0353 + conj(IT_0205)
      *IT_0356 + IT_0205*IT_0359 + conj(IT_0108)*IT_0362 + IT_0108*IT_0365 +
       conj(IT_0181)*IT_0368 + IT_0181*IT_0371)*IT_0377 + IT_0382*(conj(IT_0272)
      *IT_0383 + IT_0272*IT_0384) + IT_0385*(conj(IT_0271)*IT_0383 + conj
      (IT_0272)*IT_0386) + IT_0387*(IT_0306*IT_0383 + IT_0307*IT_0384 + IT_0329
      *IT_0386 + IT_0330*IT_0388) + IT_0390*(conj(IT_0205)*IT_0391 + IT_0205
      *IT_0392) + 4*IT_0396*(IT_0306*IT_0307 + IT_0329*(IT_0328 + 0.25*IT_0397))
       + (-6)*IT_0208*(IT_0307*IT_0329 + IT_0330*(conj(IT_0305) + (
      -0.166666666666667)*IT_0398)) + (conj(IT_0215)*IT_0216 + conj(IT_0237)
      *IT_0263 + IT_0237*IT_0265 + conj(IT_0214)*IT_0273 + conj(IT_0260)*IT_0338
       + IT_0260*IT_0340)*IT_0399 + IT_0400*(IT_0242*IT_0272 + IT_0239*conj
      (IT_0272) + conj(IT_0214)*IT_0383 + IT_0214*IT_0384 + conj(IT_0260)
      *IT_0401 + IT_0260*IT_0402) + IT_0406*((conj(IT_0055) + conj(IT_0108))
      *IT_0383 + (IT_0055 + IT_0108)*IT_0384 + IT_0335*IT_0386 + (IT_0042 +
       IT_0181)*IT_0388 + conj(IT_0131)*IT_0401 + IT_0131*IT_0402 + conj(IT_0205
      )*IT_0407 + IT_0205*IT_0408) + (IT_0263*IT_0329 + IT_0266*IT_0330 + conj
      (IT_0305)*IT_0338 + IT_0305*IT_0341 + conj(IT_0260)*IT_0348 + IT_0260
      *IT_0351 + conj(IT_0214)*IT_0360 + IT_0214*IT_0363)*IT_0409 + (conj
      (IT_0272)*IT_0360 + IT_0272*IT_0363 + IT_0329*IT_0401 + IT_0330*IT_0402 +
       conj(IT_0305)*IT_0407 + IT_0305*IT_0408)*IT_0410 + 2*IT_0411*IT_0420 + 8
      *conj(IT_0042)*(IT_0206*IT_0422 + IT_0414*IT_0423 + IT_0356*IT_0428 +
       IT_0350*IT_0433) + 8*IT_0042*(IT_0207*IT_0422 + IT_0417*IT_0423 + IT_0359
      *IT_0428 + IT_0353*IT_0433) + 8*conj(IT_0055)*(IT_0132*IT_0422 + IT_0418
      *IT_0423 + IT_0350*IT_0428 + IT_0356*IT_0433) + 8*IT_0055*(IT_0133*IT_0422
       + IT_0419*IT_0423 + IT_0353*IT_0428 + IT_0359*IT_0433);
    return create_ccomplex_return(IT_0434);
}

