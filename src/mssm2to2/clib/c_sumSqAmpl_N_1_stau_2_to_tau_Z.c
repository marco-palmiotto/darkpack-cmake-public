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
    const ccomplex_t IT_0047 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0048 = IT_0032*IT_0047;
    const ccomplex_t IT_0049 = IT_0035*IT_0047;
    const ccomplex_t IT_0050 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0051 = IT_0032*IT_0050;
    const ccomplex_t IT_0052 = IT_0035*IT_0050;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0048 + IT_0049 + 
      -IT_0051 + -IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0056 = IT_0012*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = N_d3*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0059 = IT_0021*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 0.5*IT_0060);
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0068 = IT_0032*IT_0067;
    const ccomplex_t IT_0069 = IT_0035*IT_0067;
    const ccomplex_t IT_0070 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0071 = IT_0032*IT_0070;
    const ccomplex_t IT_0072 = IT_0035*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0076 = IT_0012*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_d4*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0079 = IT_0021*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 0.5*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0040;
    const ccomplex_t IT_0088 = IT_0044*IT_0087;
    const ccomplex_t IT_0089 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0090 = IT_0032*IT_0089;
    const ccomplex_t IT_0091 = IT_0035*IT_0089;
    const ccomplex_t IT_0092 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0093 = IT_0032*IT_0092;
    const ccomplex_t IT_0094 = IT_0035*IT_0092;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0090 + IT_0091 + 
      -IT_0093 + -IT_0094);
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0098 = IT_0012*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = N_d2*e_em*m_tau*IT_0022*U_stau_01;
    const ccomplex_t IT_0101 = IT_0021*IT_0100;
    const ccomplex_t IT_0102 = 1.4142135623731*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*(IT_0099 + 0.5*IT_0102);
    const ccomplex_t IT_0104 = -IT_0103;
    const ccomplex_t IT_0105 = cpow(s_23 + (-0.5)*IT_0008 + (-0.5)*IT_0009 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = IT_0096*IT_0107;
    const ccomplex_t IT_0109 = 0.5*IT_0066 + 0.5*IT_0086 + 0.25*IT_0088 + 0.5
      *IT_0108;
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0046 + IT_0110;
    const ccomplex_t IT_0112 = conj(IT_0046) + -conj(IT_0109);
    const ccomplex_t IT_0113 = e_em*IT_0032;
    const ccomplex_t IT_0114 = e_em*IT_0035;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0113 + -IT_0114);
    const ccomplex_t IT_0116 = 0.5*IT_0115;
    const ccomplex_t IT_0117 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0118 = IT_0012*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0121 = IT_0020*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = conj(N_d1)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0124 = IT_0021*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0119 + IT_0122 + 
      -IT_0125);
    const ccomplex_t IT_0127 = (-0.5)*IT_0126;
    const ccomplex_t IT_0128 = IT_0116*IT_0127;
    const ccomplex_t IT_0129 = IT_0010*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = (-0.5)*IT_0130;
    const ccomplex_t IT_0132 = IT_0042*IT_0127;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = IT_0087*IT_0133;
    const ccomplex_t IT_0135 = (-0.25)*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0137 = IT_0032*IT_0136;
    const ccomplex_t IT_0138 = IT_0035*IT_0136;
    const ccomplex_t IT_0139 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0140 = IT_0032*IT_0139;
    const ccomplex_t IT_0141 = IT_0035*IT_0139;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0137 + IT_0138 + 
      -IT_0140 + -IT_0141);
    const ccomplex_t IT_0143 = (-0.5)*IT_0142;
    const ccomplex_t IT_0144 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0145 = IT_0012*IT_0144;
    const ccomplex_t IT_0146 = 1.4142135623731*IT_0145;
    const ccomplex_t IT_0147 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0148 = IT_0020*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = conj(N_d2)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0151 = IT_0021*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0146 + IT_0149 + 
      -IT_0152);
    const ccomplex_t IT_0154 = (-0.5)*IT_0153;
    const ccomplex_t IT_0155 = IT_0105*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = IT_0143*IT_0156;
    const ccomplex_t IT_0158 = IT_0041*IT_0133;
    const ccomplex_t IT_0159 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0160 = IT_0032*IT_0159;
    const ccomplex_t IT_0161 = IT_0035*IT_0159;
    const ccomplex_t IT_0162 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0163 = IT_0032*IT_0162;
    const ccomplex_t IT_0164 = IT_0035*IT_0162;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*(IT_0160 + IT_0161 + 
      -IT_0163 + -IT_0164);
    const ccomplex_t IT_0166 = (-0.5)*IT_0165;
    const ccomplex_t IT_0167 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0168 = IT_0012*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0171 = IT_0020*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_d3)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0174 = IT_0021*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0169 + IT_0172 + 
      -IT_0175);
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0063*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = IT_0166*IT_0179;
    const ccomplex_t IT_0181 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0182 = IT_0032*IT_0181;
    const ccomplex_t IT_0183 = IT_0035*IT_0181;
    const ccomplex_t IT_0184 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0185 = IT_0032*IT_0184;
    const ccomplex_t IT_0186 = IT_0035*IT_0184;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0182 + IT_0183 + 
      -IT_0185 + -IT_0186);
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0190 = IT_0012*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0193 = IT_0020*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = conj(N_d4)*e_em*m_tau*IT_0022*U_stau_11;
    const ccomplex_t IT_0196 = IT_0021*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0191 + IT_0194 + 
      -IT_0197);
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0083*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = IT_0188*IT_0201;
    const ccomplex_t IT_0203 = (-0.5)*IT_0157 + (-0.25)*IT_0158 + (-0.5)
      *IT_0180 + (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = IT_0135 + -IT_0203;
    const ccomplex_t IT_0205 = conj(IT_0135) + -conj(IT_0203);
    const ccomplex_t IT_0206 = pow(m_Z, 2);
    const ccomplex_t IT_0207 = s_13*IT_0206;
    const ccomplex_t IT_0208 = (-0.5)*IT_0004;
    const ccomplex_t IT_0209 = (-2)*IT_0000;
    const ccomplex_t IT_0210 = IT_0207 + IT_0208 + IT_0209;
    const ccomplex_t IT_0211 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = s_13 + IT_0212;
    const ccomplex_t IT_0214 = s_13*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0012*IT_0013 +
       -IT_0011*IT_0020)*conj(U_stau_00)*U_stau_01 + IT_0012*IT_0013*conj
      (U_stau_10)*U_stau_11);
    const ccomplex_t IT_0216 = (-2)*IT_0215;
    const ccomplex_t IT_0217 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0218 = IT_0012*IT_0217;
    const ccomplex_t IT_0219 = 1.4142135623731*IT_0218;
    const ccomplex_t IT_0220 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0221 = IT_0020*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = conj(N_d1)*e_em*m_tau*IT_0022*U_stau_10;
    const ccomplex_t IT_0224 = IT_0021*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0219 + IT_0222 + 
      -IT_0225);
    const ccomplex_t IT_0227 = (-0.5)*IT_0226;
    const ccomplex_t IT_0228 = cpow(s_13 + (-0.5)*IT_0007 + (-0.5)*IT_0008 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0229 = IT_0227*IT_0228;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = IT_0216*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0012*IT_0013 +
       -IT_0011*IT_0020)*U_stau_01*conj(U_stau_01) + IT_0012*IT_0013*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0233 = (-2)*IT_0232;
    const ccomplex_t IT_0234 = cpow(s_13 + (-0.5)*IT_0007 + (-0.5)*IT_0008 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0235 = IT_0127*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0233*IT_0236;
    const ccomplex_t IT_0238 = -IT_0130 + (-0.5)*IT_0231 + (-0.5)*IT_0237;
    const ccomplex_t IT_0239 = 0.5*IT_0134;
    const ccomplex_t IT_0240 = 2*IT_0239;
    const ccomplex_t IT_0241 = IT_0157 + 0.5*IT_0158 + IT_0180 + IT_0202 + 0.5
      *IT_0231 + 0.5*IT_0237;
    const ccomplex_t IT_0242 = (-2)*IT_0241;
    const ccomplex_t IT_0243 = IT_0240 + IT_0242;
    const ccomplex_t IT_0244 = 2*conj(IT_0239);
    const ccomplex_t IT_0245 = (-2)*conj(IT_0241);
    const ccomplex_t IT_0246 = IT_0244 + IT_0245;
    const ccomplex_t IT_0247 = IT_0027*IT_0234;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0233*IT_0248;
    const ccomplex_t IT_0250 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0251 = IT_0012*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = N_d1*e_em*m_tau*IT_0022*U_stau_00;
    const ccomplex_t IT_0254 = IT_0021*IT_0253;
    const ccomplex_t IT_0255 = 1.4142135623731*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*(IT_0252 + 0.5*IT_0255);
    const ccomplex_t IT_0257 = -IT_0256;
    const ccomplex_t IT_0258 = IT_0228*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = IT_0216*IT_0259;
    const ccomplex_t IT_0261 = IT_0030 + 0.5*IT_0249 + 0.5*IT_0260;
    const ccomplex_t IT_0262 = (-0.5)*IT_0045;
    const ccomplex_t IT_0263 = 2*IT_0262;
    const ccomplex_t IT_0264 = -IT_0066 + -IT_0086 + (-0.5)*IT_0088 + -IT_0108
       + (-0.5)*IT_0249 + (-0.5)*IT_0260;
    const ccomplex_t IT_0265 = (-2)*IT_0264;
    const ccomplex_t IT_0266 = IT_0263 + IT_0265;
    const ccomplex_t IT_0267 = 2*conj(IT_0262);
    const ccomplex_t IT_0268 = (-2)*conj(IT_0264);
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = m_N_1*IT_0213;
    const ccomplex_t IT_0271 = m_N_1*IT_0087;
    const ccomplex_t IT_0272 = m_N_1*IT_0041;
    const ccomplex_t IT_0273 = IT_0133*(IT_0271 + -IT_0272);
    const ccomplex_t IT_0274 = 0.25*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0276 = IT_0027*IT_0116;
    const ccomplex_t IT_0277 = IT_0010*IT_0275*IT_0276;
    const ccomplex_t IT_0278 = m_tau*IT_0030;
    const ccomplex_t IT_0279 = m_N_3*IT_0054;
    const ccomplex_t IT_0280 = m_N_1*IT_0166;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = IT_0279 + IT_0281;
    const ccomplex_t IT_0283 = IT_0179*IT_0282;
    const ccomplex_t IT_0284 = -IT_0272;
    const ccomplex_t IT_0285 = IT_0271 + IT_0284;
    const ccomplex_t IT_0286 = IT_0133*IT_0285;
    const ccomplex_t IT_0287 = m_N_2*IT_0096;
    const ccomplex_t IT_0288 = m_N_1*IT_0143;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = IT_0156*IT_0290;
    const ccomplex_t IT_0292 = m_N_4*IT_0074;
    const ccomplex_t IT_0293 = m_N_1*IT_0188;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = IT_0292 + IT_0294;
    const ccomplex_t IT_0296 = IT_0201*IT_0295;
    const ccomplex_t IT_0297 = 0.5*IT_0277 + (-0.5)*IT_0278 + 0.5*IT_0283 +
       0.25*IT_0286 + 0.5*IT_0291 + 0.5*IT_0296;
    const ccomplex_t IT_0298 = conj(IT_0274) + conj(IT_0297);
    const ccomplex_t IT_0299 = IT_0274 + IT_0297;
    const ccomplex_t IT_0300 = IT_0044*(IT_0271 + -IT_0272);
    const ccomplex_t IT_0301 = 0.25*IT_0300;
    const ccomplex_t IT_0302 = IT_0014*IT_0127;
    const ccomplex_t IT_0303 = IT_0010*IT_0275*IT_0302;
    const ccomplex_t IT_0304 = m_N_2*IT_0143;
    const ccomplex_t IT_0305 = m_N_1*IT_0096;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0304 + IT_0306;
    const ccomplex_t IT_0308 = IT_0107*IT_0307;
    const ccomplex_t IT_0309 = m_N_3*IT_0166;
    const ccomplex_t IT_0310 = m_N_1*IT_0054;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = IT_0309 + IT_0311;
    const ccomplex_t IT_0313 = IT_0065*IT_0312;
    const ccomplex_t IT_0314 = m_tau*IT_0130;
    const ccomplex_t IT_0315 = IT_0044*IT_0285;
    const ccomplex_t IT_0316 = m_N_4*IT_0188;
    const ccomplex_t IT_0317 = m_N_1*IT_0074;
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = IT_0316 + IT_0318;
    const ccomplex_t IT_0320 = IT_0085*IT_0319;
    const ccomplex_t IT_0321 = (-0.5)*IT_0303 + (-0.5)*IT_0308 + (-0.5)
      *IT_0313 + 0.5*IT_0314 + 0.25*IT_0315 + (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = conj(IT_0301) + conj(IT_0321);
    const ccomplex_t IT_0323 = IT_0301 + IT_0321;
    const ccomplex_t IT_0324 = m_N_1*m_tau;
    const ccomplex_t IT_0325 = pow(s_14, 2);
    const ccomplex_t IT_0326 = IT_0001*IT_0325;
    const ccomplex_t IT_0327 = -IT_0007;
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = IT_0324*IT_0328;
    const ccomplex_t IT_0330 = 2*IT_0264;
    const ccomplex_t IT_0331 = (-2)*IT_0262;
    const ccomplex_t IT_0332 = IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = 2*conj(IT_0264);
    const ccomplex_t IT_0334 = (-2)*conj(IT_0262);
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = s_13*s_14;
    const ccomplex_t IT_0337 = s_34*IT_0007;
    const ccomplex_t IT_0338 = -IT_0337;
    const ccomplex_t IT_0339 = IT_0336 + IT_0338;
    const ccomplex_t IT_0340 = conj(IT_0131) + conj(IT_0135);
    const ccomplex_t IT_0341 = IT_0131 + IT_0135;
    const ccomplex_t IT_0342 = conj(IT_0031) + conj(IT_0046);
    const ccomplex_t IT_0343 = 2*IT_0241;
    const ccomplex_t IT_0344 = (-2)*IT_0239;
    const ccomplex_t IT_0345 = IT_0343 + IT_0344;
    const ccomplex_t IT_0346 = 2*conj(IT_0241);
    const ccomplex_t IT_0347 = (-2)*conj(IT_0239);
    const ccomplex_t IT_0348 = IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = s_14*IT_0206;
    const ccomplex_t IT_0350 = IT_0001*IT_0349;
    const ccomplex_t IT_0351 = (-0.25)*IT_0350;
    const ccomplex_t IT_0352 = s_14 + IT_0351;
    const ccomplex_t IT_0353 = m_tau*IT_0352;
    const ccomplex_t IT_0354 = (-4)*IT_0353;
    const ccomplex_t IT_0355 = 2*IT_0301;
    const ccomplex_t IT_0356 = 2*IT_0321;
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = 2*conj(IT_0301);
    const ccomplex_t IT_0359 = 2*conj(IT_0321);
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = 2*IT_0274;
    const ccomplex_t IT_0362 = 2*IT_0297;
    const ccomplex_t IT_0363 = IT_0361 + IT_0362;
    const ccomplex_t IT_0364 = 2*conj(IT_0274);
    const ccomplex_t IT_0365 = 2*conj(IT_0297);
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = s_34*IT_0206;
    const ccomplex_t IT_0368 = IT_0001*IT_0367;
    const ccomplex_t IT_0369 = 0.5*IT_0368;
    const ccomplex_t IT_0370 = s_34 + IT_0369;
    const ccomplex_t IT_0371 = m_N_1*IT_0370;
    const ccomplex_t IT_0372 = 2*IT_0371;
    const ccomplex_t IT_0373 = 0.5*IT_0350;
    const ccomplex_t IT_0374 = s_14 + IT_0373;
    const ccomplex_t IT_0375 = m_tau*IT_0374;
    const ccomplex_t IT_0376 = 2*IT_0375;
    const ccomplex_t IT_0377 = (-2)*IT_0301;
    const ccomplex_t IT_0378 = (-2)*IT_0321;
    const ccomplex_t IT_0379 = IT_0377 + IT_0378;
    const ccomplex_t IT_0380 = (-2)*conj(IT_0301);
    const ccomplex_t IT_0381 = (-2)*conj(IT_0321);
    const ccomplex_t IT_0382 = IT_0380 + IT_0381;
    const ccomplex_t IT_0383 = (-2)*IT_0274;
    const ccomplex_t IT_0384 = (-2)*IT_0297;
    const ccomplex_t IT_0385 = IT_0383 + IT_0384;
    const ccomplex_t IT_0386 = (-2)*conj(IT_0274);
    const ccomplex_t IT_0387 = (-2)*conj(IT_0297);
    const ccomplex_t IT_0388 = IT_0386 + IT_0387;
    const ccomplex_t IT_0389 = (-0.25)*IT_0368;
    const ccomplex_t IT_0390 = s_34 + IT_0389;
    const ccomplex_t IT_0391 = m_N_1*IT_0390;
    const ccomplex_t IT_0392 = (-4)*IT_0391;
    const ccomplex_t IT_0393 = pow(s_34, 2);
    const ccomplex_t IT_0394 = IT_0001*IT_0393;
    const ccomplex_t IT_0395 = -IT_0008;
    const ccomplex_t IT_0396 = IT_0394 + IT_0395;
    const ccomplex_t IT_0397 = IT_0324*IT_0396;
    const ccomplex_t IT_0398 = 2*IT_0261;
    const ccomplex_t IT_0399 = 2*conj(IT_0261);
    const ccomplex_t IT_0400 = s_13*IT_0396;
    const ccomplex_t IT_0401 = 2*IT_0238;
    const ccomplex_t IT_0402 = m_N_1*IT_0396;
    const ccomplex_t IT_0403 = 2*conj(IT_0238);
    const ccomplex_t IT_0404 = IT_0207 + IT_0209;
    const ccomplex_t IT_0405 = IT_0001*IT_0404;
    const ccomplex_t IT_0406 = (-0.5)*IT_0405;
    const ccomplex_t IT_0407 = s_13 + IT_0406;
    const ccomplex_t IT_0408 = 4*IT_0274;
    const ccomplex_t IT_0409 = (-6)*conj(IT_0301);
    const ccomplex_t IT_0410 = s_13*IT_0328;
    const ccomplex_t IT_0411 = IT_0213*IT_0324;
    const ccomplex_t IT_0412 = (-2)*IT_0261;
    const ccomplex_t IT_0413 = (-2)*conj(IT_0261);
    const ccomplex_t IT_0414 = s_13*s_34;
    const ccomplex_t IT_0415 = s_14*IT_0008;
    const ccomplex_t IT_0416 = -IT_0415;
    const ccomplex_t IT_0417 = IT_0414 + IT_0416;
    const ccomplex_t IT_0418 = (-2)*IT_0238;
    const ccomplex_t IT_0419 = (-2)*conj(IT_0238);
    const ccomplex_t IT_0420 = m_tau*IT_0328;
    const ccomplex_t IT_0421 = m_tau*IT_0213;
    const ccomplex_t IT_0422 = IT_0206*IT_0324;
    const ccomplex_t IT_0423 = 3*IT_0422;
    const ccomplex_t IT_0424 = 2*IT_0109;
    const ccomplex_t IT_0425 = (-2)*IT_0046;
    const ccomplex_t IT_0426 = IT_0424 + IT_0425;
    const ccomplex_t IT_0427 = 2*conj(IT_0109);
    const ccomplex_t IT_0428 = (-2)*conj(IT_0046);
    const ccomplex_t IT_0429 = IT_0427 + IT_0428;
    const ccomplex_t IT_0430 = 2*IT_0031;
    const ccomplex_t IT_0431 = 2*conj(IT_0031);
    const ccomplex_t IT_0432 = (-2)*IT_0135;
    const ccomplex_t IT_0433 = (-2)*conj(IT_0135);
    const ccomplex_t IT_0434 = (-3)*IT_0422;
    const ccomplex_t IT_0435 = conj(IT_0046)*IT_0135 + IT_0046*conj(IT_0135) +
       0.5*conj(IT_0203)*IT_0426 + 0.5*IT_0203*IT_0429 + 0.5*conj(IT_0131)
      *IT_0430 + 0.5*IT_0131*IT_0431 + 0.5*conj(IT_0109)*IT_0432 + 0.5*IT_0109
      *IT_0433;
    const ccomplex_t IT_0436 = 8*IT_0006*(conj(IT_0031)*IT_0111 + IT_0031
      *IT_0112 + conj(IT_0131)*IT_0204 + IT_0131*IT_0205) + (-4)*(IT_0031*conj
      (IT_0031) + IT_0111*IT_0112 + IT_0131*conj(IT_0131) + IT_0204*IT_0205)
      *IT_0210 + IT_0214*(conj(IT_0238)*IT_0243 + IT_0238*IT_0246 + conj(IT_0261
      )*IT_0266 + IT_0261*IT_0269) + IT_0270*(IT_0243*IT_0298 + IT_0246*IT_0299 
      + IT_0266*IT_0322 + IT_0269*IT_0323) + IT_0329*(conj(IT_0239)*IT_0266 +
       IT_0239*IT_0269 + conj(IT_0241)*IT_0332 + IT_0241*IT_0335) + IT_0339*(
      (IT_0031 + IT_0046)*IT_0269 + conj(IT_0109)*IT_0332 + IT_0109*IT_0335 +
       IT_0243*IT_0340 + IT_0246*IT_0341 + IT_0266*IT_0342 + conj(IT_0203)
      *IT_0345 + IT_0203*IT_0348) + IT_0354*(conj(IT_0131)*IT_0357 + IT_0131
      *IT_0360 + conj(IT_0031)*IT_0363 + IT_0031*IT_0366) + (conj(IT_0031)
      *IT_0357 + IT_0031*IT_0360 + conj(IT_0131)*IT_0363 + IT_0131*IT_0366)
      *IT_0372 + IT_0376*(conj(IT_0203)*IT_0357 + IT_0203*IT_0360 + conj(IT_0109
      )*IT_0363 + IT_0109*IT_0366 + conj(IT_0135)*IT_0379 + IT_0135*IT_0382 +
       conj(IT_0046)*IT_0385 + IT_0046*IT_0388) + (conj(IT_0109)*IT_0357 +
       IT_0109*IT_0360 + conj(IT_0203)*IT_0363 + IT_0203*IT_0366 + conj(IT_0046)
      *IT_0379 + IT_0046*IT_0382 + conj(IT_0135)*IT_0385 + IT_0135*IT_0388)
      *IT_0392 + IT_0397*(conj(IT_0238)*IT_0398 + IT_0238*IT_0399) + IT_0400*
      (conj(IT_0261)*IT_0398 + conj(IT_0238)*IT_0401) + IT_0402*(IT_0322*IT_0398
       + IT_0323*IT_0399 + IT_0298*IT_0401 + IT_0299*IT_0403) + 4*IT_0407*
      (IT_0322*IT_0323 + IT_0298*(IT_0297 + 0.25*IT_0408)) + (-6)*IT_0324*
      (IT_0298*IT_0323 + IT_0299*(conj(IT_0321) + (-0.166666666666667)*IT_0409))
       + (conj(IT_0239)*IT_0240 + conj(IT_0262)*IT_0263 + conj(IT_0264)*IT_0332 
      + IT_0264*IT_0334 + conj(IT_0241)*IT_0345 + IT_0241*IT_0347)*IT_0410 +
       IT_0411*(conj(IT_0238)*IT_0266 + IT_0238*IT_0269 + conj(IT_0239)*IT_0398 
      + IT_0239*IT_0399 + conj(IT_0241)*IT_0412 + IT_0241*IT_0413) + IT_0417*
      (IT_0342*IT_0398 + (IT_0031 + IT_0046)*IT_0399 + IT_0340*IT_0401 + IT_0341
      *IT_0403 + conj(IT_0109)*IT_0412 + IT_0109*IT_0413 + conj(IT_0203)*IT_0418
       + IT_0203*IT_0419) + (IT_0298*IT_0332 + IT_0299*IT_0335 + conj(IT_0321)
      *IT_0345 + IT_0321*IT_0348 + conj(IT_0241)*IT_0355 + IT_0241*IT_0358 +
       conj(IT_0239)*IT_0377 + IT_0239*IT_0380)*IT_0420 + (conj(IT_0238)*IT_0377
       + IT_0238*IT_0380 + IT_0298*IT_0412 + IT_0299*IT_0413 + conj(IT_0321)
      *IT_0418 + IT_0321*IT_0419)*IT_0421 + IT_0423*(conj(IT_0131)*IT_0426 +
       IT_0131*IT_0429 + conj(IT_0203)*IT_0430 + IT_0203*IT_0431 + conj(IT_0031)
      *IT_0432 + IT_0031*IT_0433) + 2*IT_0434*IT_0435;
    return create_ccomplex_return(IT_0436);
}

