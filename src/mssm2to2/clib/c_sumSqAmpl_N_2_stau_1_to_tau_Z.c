#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_stau_1_to_tau_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_stau_1_to_tau_Z(
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = sin(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_01) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_11));
    const ccomplex_t IT_0005 = (-2)*IT_0004;
    const ccomplex_t IT_0006 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = 1.4142135623731*IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = N_d2*e_em*m_tau*IT_0012*U_stau_01;
    const ccomplex_t IT_0014 = IT_0011*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0008 + 0.5*IT_0015);
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = pow(m_N_2, 2);
    const ccomplex_t IT_0019 = pow(m_tau, 2);
    const ccomplex_t IT_0020 = cpow(s_13 + (-0.5)*IT_0018 + (-0.5)*IT_0019 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = IT_0005*IT_0022;
    const ccomplex_t IT_0024 = IT_0001*IT_0002;
    const ccomplex_t IT_0025 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0026 = IT_0024*IT_0025;
    const ccomplex_t IT_0027 = IT_0000*IT_0003;
    const ccomplex_t IT_0028 = IT_0025*IT_0027;
    const ccomplex_t IT_0029 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0030 = IT_0024*IT_0029;
    const ccomplex_t IT_0031 = IT_0027*IT_0029;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0026 + IT_0028 + 
      -IT_0030 + -IT_0031);
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = N_d3*e_em*m_tau*IT_0012*U_stau_00;
    const ccomplex_t IT_0038 = IT_0011*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = pow(m_stau_1, 2);
    const ccomplex_t IT_0043 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0042 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0041*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = IT_0033*IT_0045;
    const ccomplex_t IT_0047 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0048 = IT_0024*IT_0047;
    const ccomplex_t IT_0049 = IT_0027*IT_0047;
    const ccomplex_t IT_0050 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0051 = IT_0024*IT_0050;
    const ccomplex_t IT_0052 = IT_0027*IT_0050;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0048 + IT_0049 + 
      -IT_0051 + -IT_0052);
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0056 = IT_0001*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = N_d4*e_em*m_tau*IT_0012*U_stau_00;
    const ccomplex_t IT_0059 = IT_0011*IT_0058;
    const ccomplex_t IT_0060 = 1.4142135623731*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0057 + 0.5*IT_0060);
    const ccomplex_t IT_0062 = -IT_0061;
    const ccomplex_t IT_0063 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0042 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*IT_0064;
    const ccomplex_t IT_0066 = IT_0054*IT_0065;
    const ccomplex_t IT_0067 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0068 = IT_0024*IT_0067;
    const ccomplex_t IT_0069 = IT_0027*IT_0067;
    const ccomplex_t IT_0070 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0071 = IT_0024*IT_0070;
    const ccomplex_t IT_0072 = IT_0027*IT_0070;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0068 + IT_0069 + 
      -IT_0071 + -IT_0072);
    const ccomplex_t IT_0074 = (-0.5)*IT_0073;
    const ccomplex_t IT_0075 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0076 = IT_0001*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_d1*e_em*m_tau*IT_0012*U_stau_00;
    const ccomplex_t IT_0079 = IT_0011*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 0.5*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0042 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0084 = IT_0082*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0074*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0002 +
       -IT_0000*IT_0003)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0002*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0088 = (-2)*IT_0087;
    const ccomplex_t IT_0089 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0090 = IT_0001*IT_0089;
    const ccomplex_t IT_0091 = 1.4142135623731*IT_0090;
    const ccomplex_t IT_0092 = N_d2*e_em*m_tau*IT_0012*U_stau_00;
    const ccomplex_t IT_0093 = IT_0011*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + 0.5*IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = cpow(s_13 + (-0.5)*IT_0018 + (-0.5)*IT_0019 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = IT_0088*IT_0099;
    const ccomplex_t IT_0101 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0102 = IT_0024*IT_0101;
    const ccomplex_t IT_0103 = IT_0027*IT_0101;
    const ccomplex_t IT_0104 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0105 = IT_0024*IT_0104;
    const ccomplex_t IT_0106 = IT_0027*IT_0104;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*(IT_0102 + IT_0103 + 
      -IT_0105 + -IT_0106);
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = cpow(s_23 + (-0.5)*IT_0019 + (-0.5)*IT_0042 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0110 = IT_0096*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = IT_0108*IT_0111;
    const ccomplex_t IT_0113 = 0.5*IT_0107;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = (-0.5)*IT_0023 + -IT_0046 + -IT_0066 + IT_0086 
      + (-0.5)*IT_0100 + (-0.5)*IT_0112 + 0.5*IT_0114;
    const ccomplex_t IT_0116 = pow(m_Z, -2);
    const ccomplex_t IT_0117 = pow(s_14, 2);
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = -IT_0018;
    const ccomplex_t IT_0120 = IT_0118 + IT_0119;
    const ccomplex_t IT_0121 = s_13*IT_0120;
    const ccomplex_t IT_0122 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0123 = IT_0024*IT_0122;
    const ccomplex_t IT_0124 = IT_0027*IT_0122;
    const ccomplex_t IT_0125 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0126 = IT_0024*IT_0125;
    const ccomplex_t IT_0127 = IT_0027*IT_0125;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0123 + IT_0124 + 
      -IT_0126 + -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0131 = IT_0001*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0134 = IT_0003*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = conj(N_d3)*e_em*m_tau*IT_0012*U_stau_10;
    const ccomplex_t IT_0137 = IT_0011*IT_0136;
    const ccomplex_t IT_0138 = 1.4142135623731*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*(IT_0132 + IT_0135 + 
      -IT_0138);
    const ccomplex_t IT_0140 = (-0.5)*IT_0139;
    const ccomplex_t IT_0141 = IT_0043*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = IT_0129*IT_0142;
    const ccomplex_t IT_0144 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0145 = IT_0024*IT_0144;
    const ccomplex_t IT_0146 = IT_0027*IT_0144;
    const ccomplex_t IT_0147 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0148 = IT_0024*IT_0147;
    const ccomplex_t IT_0149 = IT_0027*IT_0147;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0145 + IT_0146 + 
      -IT_0148 + -IT_0149);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0153 = IT_0001*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0156 = IT_0003*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = conj(N_d1)*e_em*m_tau*IT_0012*U_stau_10;
    const ccomplex_t IT_0159 = IT_0011*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0154 + IT_0157 + 
      -IT_0160);
    const ccomplex_t IT_0162 = (-0.5)*IT_0161;
    const ccomplex_t IT_0163 = IT_0083*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = IT_0151*IT_0164;
    const ccomplex_t IT_0166 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0167 = IT_0001*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0170 = IT_0003*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = conj(N_d2)*e_em*m_tau*IT_0012*U_stau_10;
    const ccomplex_t IT_0173 = IT_0011*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0168 + IT_0171 + 
      -IT_0174);
    const ccomplex_t IT_0176 = (-0.5)*IT_0175;
    const ccomplex_t IT_0177 = IT_0109*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = IT_0108*IT_0178;
    const ccomplex_t IT_0180 = IT_0113*IT_0178;
    const ccomplex_t IT_0181 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0182 = IT_0024*IT_0181;
    const ccomplex_t IT_0183 = IT_0027*IT_0181;
    const ccomplex_t IT_0184 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0185 = IT_0024*IT_0184;
    const ccomplex_t IT_0186 = IT_0027*IT_0184;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0182 + IT_0183 + 
      -IT_0185 + -IT_0186);
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0190 = IT_0001*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0193 = IT_0003*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = conj(N_d4)*e_em*m_tau*IT_0012*U_stau_10;
    const ccomplex_t IT_0196 = IT_0011*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*(IT_0191 + IT_0194 + 
      -IT_0197);
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0063*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = IT_0188*IT_0201;
    const ccomplex_t IT_0203 = IT_0097*IT_0176;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*IT_0203;
    const ccomplex_t IT_0205 = IT_0088*IT_0204;
    const ccomplex_t IT_0206 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0207 = IT_0001*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0210 = IT_0003*IT_0209;
    const ccomplex_t IT_0211 = 1.4142135623731*IT_0210;
    const ccomplex_t IT_0212 = conj(N_d2)*e_em*m_tau*IT_0012*U_stau_11;
    const ccomplex_t IT_0213 = IT_0011*IT_0212;
    const ccomplex_t IT_0214 = 1.4142135623731*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*(IT_0208 + IT_0211 + 
      -IT_0214);
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = IT_0020*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = IT_0005*IT_0218;
    const ccomplex_t IT_0220 = IT_0143 + -IT_0165 + (-0.5)*IT_0179 + 0.5
      *IT_0180 + IT_0202 + 0.5*IT_0205 + 0.5*IT_0219;
    const ccomplex_t IT_0221 = m_N_2*m_tau;
    const ccomplex_t IT_0222 = IT_0120*IT_0221;
    const ccomplex_t IT_0223 = 0.5*IT_0046 + 0.5*IT_0066 + (-0.5)*IT_0086 +
       0.25*IT_0112 + (-0.25)*IT_0114;
    const ccomplex_t IT_0224 = s_13*s_14;
    const ccomplex_t IT_0225 = s_34*IT_0018;
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = m_tau*IT_0120;
    const ccomplex_t IT_0229 = m_N_2*IT_0113;
    const ccomplex_t IT_0230 = IT_0178*IT_0229;
    const ccomplex_t IT_0231 = m_N_1*IT_0074;
    const ccomplex_t IT_0232 = IT_0164*IT_0231;
    const ccomplex_t IT_0233 = (-0.25)*IT_0230 + (-0.5)*IT_0232;
    const ccomplex_t IT_0234 = m_N_3*IT_0033;
    const ccomplex_t IT_0235 = m_N_2*IT_0129;
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = IT_0234 + IT_0236;
    const ccomplex_t IT_0238 = IT_0142*IT_0237;
    const ccomplex_t IT_0239 = cpow(s_12 + 0.5*IT_0018 + (-0.5)*IT_0019 + 0.5
      *IT_0042 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0002;
    const ccomplex_t IT_0241 = IT_0096*IT_0240;
    const ccomplex_t IT_0242 = IT_0239*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = m_tau*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*m_tau;
    const ccomplex_t IT_0246 = e_em*IT_0024;
    const ccomplex_t IT_0247 = e_em*IT_0027;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*(IT_0246 + -IT_0247);
    const ccomplex_t IT_0249 = 0.5*IT_0248;
    const ccomplex_t IT_0250 = IT_0096*IT_0249;
    const ccomplex_t IT_0251 = IT_0239*IT_0245*IT_0250;
    const ccomplex_t IT_0252 = m_N_2*IT_0108;
    const ccomplex_t IT_0253 = IT_0178*IT_0252;
    const ccomplex_t IT_0254 = -IT_0252;
    const ccomplex_t IT_0255 = IT_0229 + IT_0254;
    const ccomplex_t IT_0256 = IT_0178*IT_0255;
    const ccomplex_t IT_0257 = m_N_2*IT_0151;
    const ccomplex_t IT_0258 = IT_0164*IT_0257;
    const ccomplex_t IT_0259 = m_N_4*IT_0054;
    const ccomplex_t IT_0260 = m_N_2*IT_0188;
    const ccomplex_t IT_0261 = -IT_0260;
    const ccomplex_t IT_0262 = IT_0259 + IT_0261;
    const ccomplex_t IT_0263 = IT_0201*IT_0262;
    const ccomplex_t IT_0264 = 0.5*IT_0238 + (-0.5)*IT_0244 + 0.5*IT_0251 +
       0.25*IT_0253 + (-0.25)*IT_0256 + 0.5*IT_0258 + 0.5*IT_0263;
    const ccomplex_t IT_0265 = conj(IT_0233) + conj(IT_0264);
    const ccomplex_t IT_0266 = pow(m_Z, 2);
    const ccomplex_t IT_0267 = s_13*IT_0266;
    const ccomplex_t IT_0268 = s_14*s_34;
    const ccomplex_t IT_0269 = (-2)*IT_0268;
    const ccomplex_t IT_0270 = IT_0267 + IT_0269;
    const ccomplex_t IT_0271 = IT_0116*IT_0270;
    const ccomplex_t IT_0272 = (-0.5)*IT_0271;
    const ccomplex_t IT_0273 = s_13 + IT_0272;
    const ccomplex_t IT_0274 = 4*IT_0233;
    const ccomplex_t IT_0275 = 0.5*IT_0023 + 0.5*IT_0100 + IT_0243;
    const ccomplex_t IT_0276 = s_14*s_34*IT_0116;
    const ccomplex_t IT_0277 = -IT_0276;
    const ccomplex_t IT_0278 = s_13 + IT_0277;
    const ccomplex_t IT_0279 = s_13*IT_0278;
    const ccomplex_t IT_0280 = IT_0176*IT_0249;
    const ccomplex_t IT_0281 = IT_0239*IT_0280;
    const ccomplex_t IT_0282 = (0 + _Complex_I*1)*IT_0281;
    const ccomplex_t IT_0283 = (-0.5)*IT_0205 + (-0.5)*IT_0219 + -IT_0282;
    const ccomplex_t IT_0284 = IT_0221*IT_0278;
    const ccomplex_t IT_0285 = 0.5*IT_0243;
    const ccomplex_t IT_0286 = m_N_2*IT_0278;
    const ccomplex_t IT_0287 = IT_0111*IT_0252;
    const ccomplex_t IT_0288 = m_N_1*IT_0151;
    const ccomplex_t IT_0289 = IT_0085*IT_0288;
    const ccomplex_t IT_0290 = 0.25*IT_0287 + 0.5*IT_0289;
    const ccomplex_t IT_0291 = m_N_2*IT_0074;
    const ccomplex_t IT_0292 = IT_0085*IT_0291;
    const ccomplex_t IT_0293 = m_N_3*IT_0129;
    const ccomplex_t IT_0294 = m_N_2*IT_0033;
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = IT_0293 + IT_0295;
    const ccomplex_t IT_0297 = IT_0045*IT_0296;
    const ccomplex_t IT_0298 = m_N_4*IT_0188;
    const ccomplex_t IT_0299 = m_N_2*IT_0054;
    const ccomplex_t IT_0300 = -IT_0299;
    const ccomplex_t IT_0301 = IT_0298 + IT_0300;
    const ccomplex_t IT_0302 = IT_0065*IT_0301;
    const ccomplex_t IT_0303 = IT_0176*IT_0240;
    const ccomplex_t IT_0304 = IT_0239*IT_0245*IT_0303;
    const ccomplex_t IT_0305 = m_tau*IT_0282;
    const ccomplex_t IT_0306 = IT_0111*IT_0229;
    const ccomplex_t IT_0307 = IT_0111*IT_0255;
    const ccomplex_t IT_0308 = (-0.5)*IT_0292 + (-0.5)*IT_0297 + (-0.5)
      *IT_0302 + (-0.5)*IT_0304 + 0.5*IT_0305 + (-0.25)*IT_0306 + (-0.25)*IT_0307;
    const ccomplex_t IT_0309 = conj(IT_0290) + conj(IT_0308);
    const ccomplex_t IT_0310 = (-2)*IT_0115;
    const ccomplex_t IT_0311 = pow(s_34, 2);
    const ccomplex_t IT_0312 = IT_0116*IT_0311;
    const ccomplex_t IT_0313 = -IT_0019;
    const ccomplex_t IT_0314 = IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = s_13*IT_0314;
    const ccomplex_t IT_0316 = IT_0221*IT_0314;
    const ccomplex_t IT_0317 = s_13*s_34;
    const ccomplex_t IT_0318 = s_14*IT_0019;
    const ccomplex_t IT_0319 = -IT_0318;
    const ccomplex_t IT_0320 = IT_0317 + IT_0319;
    const ccomplex_t IT_0321 = m_N_2*IT_0314;
    const ccomplex_t IT_0322 = (-2)*conj(IT_0115);
    const ccomplex_t IT_0323 = (-0.5)*IT_0282;
    const ccomplex_t IT_0324 = (-2)*IT_0220;
    const ccomplex_t IT_0325 = m_tau*IT_0278;
    const ccomplex_t IT_0326 = (-2)*IT_0275;
    const ccomplex_t IT_0327 = (-0.5)*IT_0143 + 0.5*IT_0165 + 0.25*IT_0179 + (
      -0.25)*IT_0180 + (-0.5)*IT_0202;
    const ccomplex_t IT_0328 = s_34*IT_0266;
    const ccomplex_t IT_0329 = IT_0116*IT_0328;
    const ccomplex_t IT_0330 = (-0.25)*IT_0329;
    const ccomplex_t IT_0331 = s_34 + IT_0330;
    const ccomplex_t IT_0332 = m_N_2*IT_0331;
    const ccomplex_t IT_0333 = (-4)*IT_0332;
    const ccomplex_t IT_0334 = 0.5*IT_0329;
    const ccomplex_t IT_0335 = s_34 + IT_0334;
    const ccomplex_t IT_0336 = m_N_2*IT_0335;
    const ccomplex_t IT_0337 = 2*IT_0336;
    const ccomplex_t IT_0338 = s_14*IT_0266;
    const ccomplex_t IT_0339 = IT_0116*IT_0338;
    const ccomplex_t IT_0340 = 0.5*IT_0339;
    const ccomplex_t IT_0341 = s_14 + IT_0340;
    const ccomplex_t IT_0342 = m_tau*IT_0341;
    const ccomplex_t IT_0343 = 2*IT_0342;
    const ccomplex_t IT_0344 = (-0.25)*IT_0339;
    const ccomplex_t IT_0345 = s_14 + IT_0344;
    const ccomplex_t IT_0346 = m_tau*IT_0345;
    const ccomplex_t IT_0347 = (-4)*IT_0346;
    const ccomplex_t IT_0348 = (-6)*conj(IT_0290);
    const ccomplex_t IT_0349 = (-0.166666666666667)*IT_0348;
    const ccomplex_t IT_0350 = conj(IT_0308) + IT_0349;
    const ccomplex_t IT_0351 = (-2)*conj(IT_0220);
    const ccomplex_t IT_0352 = (-2)*conj(IT_0275);
    const ccomplex_t IT_0353 = conj(IT_0115)*IT_0228 + conj(IT_0283)*IT_0321 +
       conj(IT_0327)*IT_0333 + conj(IT_0323)*IT_0337 + conj(IT_0223)*IT_0343 +
       conj(IT_0285)*IT_0347 + (-3)*IT_0221*IT_0350 + 0.5*IT_0286*IT_0351 + 0.5
      *IT_0325*IT_0352;
    const ccomplex_t IT_0354 = (-2)*IT_0283;
    const ccomplex_t IT_0355 = 0.5*IT_0284;
    const ccomplex_t IT_0356 = IT_0221*IT_0266;
    const ccomplex_t IT_0357 = (-3)*IT_0356;
    const ccomplex_t IT_0358 = 3*IT_0356;
    const ccomplex_t IT_0359 = pow(m_Z, 4);
    const ccomplex_t IT_0360 = s_13*IT_0359;
    const ccomplex_t IT_0361 = IT_0116*IT_0360;
    const ccomplex_t IT_0362 = (-0.5)*IT_0361;
    const ccomplex_t IT_0363 = IT_0267 + IT_0269 + IT_0362;
    const ccomplex_t IT_0364 = (-0.25)*IT_0361;
    const ccomplex_t IT_0365 = IT_0268 + IT_0364;
    const ccomplex_t IT_0366 = 0.5*IT_0227;
    const ccomplex_t IT_0367 = 0.5*IT_0320;
    const ccomplex_t IT_0368 = (-2)*conj(IT_0283);
    const ccomplex_t IT_0369 = (-0.333333333333333)*conj(IT_0275)*IT_0321;
    const ccomplex_t IT_0370 = (-0.333333333333333)*conj(IT_0220)*IT_0228;
    const ccomplex_t IT_0371 = IT_0221*(conj(IT_0233) + conj(IT_0264)) +
       IT_0369 + IT_0370;
    const ccomplex_t IT_0372 = (-0.666666666666667)*IT_0273*IT_0309;
    const ccomplex_t IT_0373 = (-0.333333333333333)*conj(IT_0223);
    const ccomplex_t IT_0374 = (-0.333333333333333)*conj(IT_0285);
    const ccomplex_t IT_0375 = (-0.333333333333333)*conj(IT_0327);
    const ccomplex_t IT_0376 = IT_0343*IT_0375;
    const ccomplex_t IT_0377 = (-0.333333333333333)*conj(IT_0323);
    const ccomplex_t IT_0378 = IT_0347*IT_0377;
    const ccomplex_t IT_0379 = (-0.166666666666667)*IT_0286;
    const ccomplex_t IT_0380 = (-2)*IT_0290*(conj(IT_0283)*IT_0325 + 3*IT_0371
       + 3*IT_0372 + 3*IT_0333*IT_0373 + 3*IT_0337*IT_0374 + 3*IT_0376 + 3
      *IT_0378 + 3*IT_0322*IT_0379);
    const ccomplex_t IT_0381 = (-0.166666666666667)*IT_0325;
    const ccomplex_t IT_0382 = conj(IT_0115)*IT_0228 + 2*IT_0265*IT_0273 +
       conj(IT_0283)*IT_0321 + (-3)*IT_0221*IT_0350 + (-3)*IT_0343*IT_0373 + (-3
      )*IT_0347*IT_0374 + (-3)*IT_0333*IT_0375 + (-3)*IT_0337*IT_0377 + (-3)
      *IT_0351*IT_0379 + (-3)*IT_0352*IT_0381;
    const ccomplex_t IT_0383 = 2*IT_0115*(conj(IT_0115)*IT_0121 + conj(IT_0220
      )*IT_0222 + conj(IT_0223)*IT_0227 + IT_0228*IT_0265) + IT_0265*IT_0273
      *IT_0274 + (conj(IT_0275)*IT_0279 + conj(IT_0283)*IT_0284 + IT_0227*conj
      (IT_0285) + IT_0286*IT_0309)*IT_0310 + 2*IT_0275*(conj(IT_0275)*IT_0315 +
       conj(IT_0283)*IT_0316 + conj(IT_0285)*IT_0320 + IT_0309*IT_0321 + 0.5
      *IT_0279*IT_0322) + (IT_0279*conj(IT_0283) + IT_0265*IT_0286 + IT_0227
      *conj(IT_0323))*IT_0324 + (conj(IT_0220)*IT_0284 + conj(IT_0223)*IT_0320 +
       IT_0265*IT_0325)*IT_0326 + 2*IT_0233*IT_0353 + (conj(IT_0308)*IT_0325 +
       IT_0320*conj(IT_0327))*IT_0354 + 2*IT_0283*(conj(IT_0283)*IT_0315 + conj
      (IT_0275)*IT_0316 + IT_0265*IT_0321 + IT_0320*conj(IT_0323) + -conj
      (IT_0290)*IT_0325 + 0.5*IT_0279*IT_0351 + IT_0322*IT_0355) + 2*IT_0220*
      (IT_0121*conj(IT_0220) + conj(IT_0115)*IT_0222 + IT_0228*IT_0309 + IT_0227
      *conj(IT_0327) + IT_0352*IT_0355) + 2*IT_0285*(conj(IT_0275)*IT_0320 +
       IT_0309*IT_0337 + IT_0265*IT_0347 + conj(IT_0323)*IT_0357 + conj(IT_0327)
      *IT_0358 + (-2)*conj(IT_0285)*IT_0363 + (-4)*conj(IT_0223)*IT_0365 +
       IT_0322*IT_0366) + 2*IT_0323*(conj(IT_0283)*IT_0320 + IT_0265*IT_0337 +
       IT_0309*IT_0347 + conj(IT_0285)*IT_0357 + conj(IT_0223)*IT_0358 + (-2)
      *conj(IT_0323)*IT_0363 + (-4)*conj(IT_0327)*IT_0365 + IT_0351*IT_0366) + 2
      *IT_0223*(conj(IT_0115)*IT_0227 + IT_0309*IT_0333 + IT_0265*IT_0343 + conj
      (IT_0327)*IT_0357 + conj(IT_0323)*IT_0358 + (-2)*conj(IT_0223)*IT_0363 + (
      -4)*conj(IT_0285)*IT_0365 + IT_0352*IT_0367) + 2*IT_0327*(conj(IT_0220)
      *IT_0227 + IT_0265*IT_0333 + IT_0309*IT_0343 + conj(IT_0223)*IT_0357 +
       conj(IT_0285)*IT_0358 + (-2)*conj(IT_0327)*IT_0363 + (-4)*conj(IT_0323)
      *IT_0365 + IT_0367*IT_0368) + IT_0380 + (-6)*IT_0308*(IT_0371 + IT_0372 +
       IT_0333*IT_0373 + IT_0337*IT_0374 + IT_0376 + IT_0378 + IT_0322*IT_0379 +
       IT_0368*IT_0381) + 2*IT_0264*IT_0382;
    return create_ccomplex_return(IT_0383);
}

