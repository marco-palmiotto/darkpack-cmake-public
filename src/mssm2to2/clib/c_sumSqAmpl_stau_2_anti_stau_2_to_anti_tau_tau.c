#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_stau_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_stau_2_to_anti_tau_tau(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = -IT_0000;
    const ccomplex_t IT_0002 = s_34 + IT_0001;
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_B3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = N_W3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = N_d3*e_em*m_tau*IT_0016*conj(U_stau_11);
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0007 + IT_0012 + 
      -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_d3*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0026 = IT_0015*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + 0.5*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0021*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0032 = pow(m_stau_2, 2);
    const ccomplex_t IT_0033 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0032 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0034 = IT_0030*IT_0031*IT_0033;
    const ccomplex_t IT_0035 = N_B4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0036 = IT_0004*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = N_W4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0039 = IT_0009*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = N_d4*e_em*m_tau*IT_0016*conj(U_stau_11);
    const ccomplex_t IT_0042 = IT_0015*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0037 + IT_0040 + 
      -IT_0043);
    const ccomplex_t IT_0045 = (-0.5)*IT_0044;
    const ccomplex_t IT_0046 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0047 = IT_0004*IT_0046;
    const ccomplex_t IT_0048 = 1.4142135623731*IT_0047;
    const ccomplex_t IT_0049 = N_d4*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0050 = IT_0015*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + 0.5*IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0045*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0056 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0032 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0054*IT_0055*IT_0056;
    const ccomplex_t IT_0058 = IT_0004*IT_0008;
    const ccomplex_t IT_0059 = e_em*IT_0058;
    const ccomplex_t IT_0060 = IT_0003*IT_0009;
    const ccomplex_t IT_0061 = e_em*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0059 + -IT_0061);
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0004*IT_0008 +
       -IT_0003*IT_0009)*U_stau_01*conj(U_stau_01) + IT_0004*IT_0008*U_stau_11
      *conj(U_stau_11));
    const ccomplex_t IT_0068 = -IT_0067;
    const ccomplex_t IT_0069 = m_tau*IT_0068;
    const ccomplex_t IT_0070 = (-2)*IT_0067;
    const ccomplex_t IT_0071 = m_tau*IT_0070;
    const ccomplex_t IT_0072 = -IT_0071;
    const ccomplex_t IT_0073 = IT_0069 + IT_0072;
    const ccomplex_t IT_0074 = IT_0066*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0008;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = IT_0073*IT_0077;
    const ccomplex_t IT_0079 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0032 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0080 = N_B1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_W1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0084 = IT_0009*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = N_d1*e_em*m_tau*IT_0016*conj(U_stau_11);
    const ccomplex_t IT_0087 = IT_0015*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*(IT_0082 + IT_0085 + 
      -IT_0088);
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0092 = IT_0004*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0095 = IT_0009*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = conj(N_d1)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0098 = IT_0015*IT_0097;
    const ccomplex_t IT_0099 = 1.4142135623731*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*(IT_0093 + IT_0096 + 
      -IT_0099);
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0090*IT_0101;
    const ccomplex_t IT_0103 = IT_0079*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = m_tau*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0107 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0108 = IT_0004*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_d1*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0109 + 0.5*IT_0112);
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0090*IT_0114;
    const ccomplex_t IT_0116 = IT_0079*IT_0106*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0118 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0032 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0119 = N_B2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0120 = IT_0004*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = N_W2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0123 = IT_0009*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = N_d2*e_em*m_tau*IT_0016*conj(U_stau_11);
    const ccomplex_t IT_0126 = IT_0015*IT_0125;
    const ccomplex_t IT_0127 = 1.4142135623731*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*(IT_0121 + IT_0124 + 
      -IT_0127);
    const ccomplex_t IT_0129 = (-0.5)*IT_0128;
    const ccomplex_t IT_0130 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0131 = IT_0004*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = N_d2*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0134 = IT_0015*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0132 + 0.5*IT_0135);
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = IT_0129*IT_0137;
    const ccomplex_t IT_0139 = IT_0117*IT_0118*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0141 = IT_0004*IT_0140;
    const ccomplex_t IT_0142 = 1.4142135623731*IT_0141;
    const ccomplex_t IT_0143 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0144 = IT_0009*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = conj(N_d3)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0147 = IT_0015*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0142 + IT_0145 + 
      -IT_0148);
    const ccomplex_t IT_0150 = (-0.5)*IT_0149;
    const ccomplex_t IT_0151 = IT_0021*IT_0150;
    const ccomplex_t IT_0152 = IT_0033*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = m_tau*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0156 = IT_0004*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0159 = IT_0009*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_d4)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0162 = IT_0015*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0157 + IT_0160 + 
      -IT_0163);
    const ccomplex_t IT_0165 = (-0.5)*IT_0164;
    const ccomplex_t IT_0166 = IT_0045*IT_0165;
    const ccomplex_t IT_0167 = IT_0056*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = m_tau*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0171 = IT_0004*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0174 = IT_0009*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = conj(N_d2)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0177 = IT_0015*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0172 + IT_0175 + 
      -IT_0178);
    const ccomplex_t IT_0180 = (-0.5)*IT_0179;
    const ccomplex_t IT_0181 = IT_0129*IT_0180;
    const ccomplex_t IT_0182 = IT_0118*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = m_tau*IT_0183;
    const ccomplex_t IT_0185 = 0.5*IT_0034 + 0.5*IT_0057 + -IT_0074 + IT_0078 
      + (-0.5)*IT_0105 + 0.5*IT_0116 + 0.5*IT_0139 + (-0.5)*IT_0154 + (-0.5)
      *IT_0169 + (-0.5)*IT_0184;
    const ccomplex_t IT_0186 = conj(N_B4)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0187 = IT_0004*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = conj(N_d4)*e_em*m_tau*IT_0016*conj(U_stau_01);
    const ccomplex_t IT_0190 = IT_0015*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*(IT_0188 + 0.5*IT_0191);
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = IT_0165*IT_0193;
    const ccomplex_t IT_0195 = IT_0055*IT_0056*IT_0194;
    const ccomplex_t IT_0196 = conj(N_B1)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0197 = IT_0004*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = conj(N_d1)*e_em*m_tau*IT_0016*conj(U_stau_01);
    const ccomplex_t IT_0200 = IT_0015*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*(IT_0198 + 0.5*IT_0201);
    const ccomplex_t IT_0203 = -IT_0202;
    const ccomplex_t IT_0204 = IT_0101*IT_0203;
    const ccomplex_t IT_0205 = IT_0079*IT_0106*IT_0204;
    const ccomplex_t IT_0206 = IT_0053*IT_0193;
    const ccomplex_t IT_0207 = IT_0056*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = m_tau*IT_0208;
    const ccomplex_t IT_0210 = conj(N_B2)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0211 = IT_0004*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = conj(N_d2)*e_em*m_tau*IT_0016*conj(U_stau_01);
    const ccomplex_t IT_0214 = IT_0015*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0212 + 0.5*IT_0215);
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0180*IT_0217;
    const ccomplex_t IT_0219 = IT_0117*IT_0118*IT_0218;
    const ccomplex_t IT_0220 = conj(N_B3)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0221 = IT_0004*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = conj(N_d3)*e_em*m_tau*IT_0016*conj(U_stau_01);
    const ccomplex_t IT_0224 = IT_0015*IT_0223;
    const ccomplex_t IT_0225 = 1.4142135623731*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*(IT_0222 + 0.5*IT_0225);
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0150*IT_0227;
    const ccomplex_t IT_0229 = IT_0031*IT_0033*IT_0228;
    const ccomplex_t IT_0230 = IT_0114*IT_0203;
    const ccomplex_t IT_0231 = IT_0079*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*IT_0231;
    const ccomplex_t IT_0233 = m_tau*IT_0232;
    const ccomplex_t IT_0234 = IT_0137*IT_0217;
    const ccomplex_t IT_0235 = IT_0118*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = m_tau*IT_0236;
    const ccomplex_t IT_0238 = IT_0029*IT_0227;
    const ccomplex_t IT_0239 = IT_0033*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = m_tau*IT_0240;
    const ccomplex_t IT_0242 = IT_0074 + -IT_0078 + 0.5*IT_0195 + 0.5*IT_0205 
      + (-0.5)*IT_0209 + 0.5*IT_0219 + 0.5*IT_0229 + (-0.5)*IT_0233 + (-0.5)
      *IT_0237 + (-0.5)*IT_0241;
    const ccomplex_t IT_0243 = sin(alpha);
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*e_em*m_tau*IT_0009*IT_0014
      *IT_0016*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0244;
    const ccomplex_t IT_0246 = mu_h*m_tau*IT_0016;
    const ccomplex_t IT_0247 = cos(alpha);
    const ccomplex_t IT_0248 = IT_0013*IT_0243;
    const ccomplex_t IT_0249 = sin(beta);
    const ccomplex_t IT_0250 = IT_0247*IT_0249;
    const ccomplex_t IT_0251 = IT_0248 + IT_0250;
    const ccomplex_t IT_0252 = cpow(IT_0003, -2);
    const ccomplex_t IT_0253 = 2*U_stau_01*((IT_0000*IT_0009*IT_0014*IT_0016
      *IT_0243 + (-0.5)*m_W*IT_0251*(IT_0009 + -IT_0008*IT_0252))*conj(U_stau_01
      ) + 1./2*IT_0009*IT_0014*IT_0246*IT_0247*conj(U_stau_11)) + U_stau_11*
      (IT_0009*IT_0014*IT_0246*IT_0247*conj(U_stau_01) + 2*(IT_0000*IT_0009
      *IT_0014*IT_0016*IT_0243 + -m_W*IT_0008*IT_0251*IT_0252)*conj(U_stau_11));
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*(e_em*IT_0253 +
       1.4142135623731*A_tau*IT_0243*(conj(U_stau_01)*U_stau_11 + U_stau_01*conj
      (U_stau_11)));
    const ccomplex_t IT_0255 = 0.5*IT_0254;
    const ccomplex_t IT_0256 = IT_0245*IT_0255;
    const ccomplex_t IT_0257 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0032 + -reg_prop, -1);
    const ccomplex_t IT_0258 = IT_0256*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = (0 + _Complex_I*1)*e_em*m_tau*IT_0009*IT_0014
      *IT_0016*IT_0247;
    const ccomplex_t IT_0261 = (-0.5)*IT_0260;
    const ccomplex_t IT_0262 = -IT_0243*IT_0249;
    const ccomplex_t IT_0263 = IT_0013*IT_0247;
    const ccomplex_t IT_0264 = IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = (-2)*U_stau_01*((IT_0000*IT_0009*IT_0014
      *IT_0016*IT_0247 + (-0.5)*m_W*(IT_0009 + -IT_0008*IT_0252)*IT_0264)*conj
      (U_stau_01) + -1./2*IT_0009*IT_0014*IT_0243*IT_0246*conj(U_stau_11)) +
       U_stau_11*(IT_0009*IT_0014*IT_0243*IT_0246*conj(U_stau_01) + (-2)*
      (IT_0000*IT_0009*IT_0014*IT_0016*IT_0247 + -m_W*IT_0008*IT_0252*IT_0264)
      *conj(U_stau_11));
    const ccomplex_t IT_0266 = (0 + _Complex_I*-1)*(e_em*IT_0265 + (
      -1.4142135623731)*A_tau*IT_0247*(conj(U_stau_01)*U_stau_11 + U_stau_01
      *conj(U_stau_11)));
    const ccomplex_t IT_0267 = (-0.5)*IT_0266;
    const ccomplex_t IT_0268 = IT_0261*IT_0267;
    const ccomplex_t IT_0269 = cpow((-2)*s_12 + (-2)*IT_0032 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0270 = IT_0268*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = -IT_0259 + -IT_0271;
    const ccomplex_t IT_0273 = conj(IT_0185) + conj(IT_0242);
    const ccomplex_t IT_0274 = IT_0066*IT_0070;
    const ccomplex_t IT_0275 = 0.5*IT_0104 + 0.5*IT_0153 + 0.5*IT_0168 + 0.5
      *IT_0183 + -IT_0274;
    const ccomplex_t IT_0276 = s_13*m_tau;
    const ccomplex_t IT_0277 = 2*s_34;
    const ccomplex_t IT_0278 = s_13*s_14;
    const ccomplex_t IT_0279 = s_34*IT_0032;
    const ccomplex_t IT_0280 = (-0.5)*IT_0279;
    const ccomplex_t IT_0281 = IT_0278 + IT_0280;
    const ccomplex_t IT_0282 = IT_0000*IT_0032;
    const ccomplex_t IT_0283 = (-0.5)*IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0286 = -IT_0285;
    const ccomplex_t IT_0287 = 2*IT_0285;
    const ccomplex_t IT_0288 = IT_0286*IT_0287;
    const ccomplex_t IT_0289 = cpow(s_12 + IT_0032 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0290 = IT_0288*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*IT_0290;
    const ccomplex_t IT_0292 = 0.5*IT_0291;
    const ccomplex_t IT_0293 = cpow(IT_0292, 2);
    const ccomplex_t IT_0294 = IT_0070*IT_0077;
    const ccomplex_t IT_0295 = 0.5*IT_0208 + 0.5*IT_0232 + 0.5*IT_0236 + 0.5
      *IT_0240 + -IT_0294;
    const ccomplex_t IT_0296 = 4*IT_0292;
    const ccomplex_t IT_0297 = s_14*m_tau;
    const ccomplex_t IT_0298 = (-2)*IT_0242;
    const ccomplex_t IT_0299 = -s_14;
    const ccomplex_t IT_0300 = s_13 + IT_0299;
    const ccomplex_t IT_0301 = m_tau*IT_0300;
    const ccomplex_t IT_0302 = (-2)*IT_0185;
    const ccomplex_t IT_0303 = e_em*m_tau*IT_0009*IT_0014*IT_0016*IT_0249;
    const ccomplex_t IT_0304 = (-0.5)*IT_0303;
    const ccomplex_t IT_0305 = e_em*mu_h*m_tau*IT_0016*conj(U_stau_01)
      *U_stau_11;
    const ccomplex_t IT_0306 = IT_0009*IT_0305;
    const ccomplex_t IT_0307 = e_em*mu_h*m_tau*IT_0016*U_stau_01*conj
      (U_stau_11);
    const ccomplex_t IT_0308 = IT_0009*IT_0307;
    const ccomplex_t IT_0309 = A_tau*conj(U_stau_01)*U_stau_11;
    const ccomplex_t IT_0310 = IT_0249*IT_0309;
    const ccomplex_t IT_0311 = 1.4142135623731*IT_0310;
    const ccomplex_t IT_0312 = A_tau*U_stau_01*conj(U_stau_11);
    const ccomplex_t IT_0313 = IT_0249*IT_0312;
    const ccomplex_t IT_0314 = 1.4142135623731*IT_0313;
    const ccomplex_t IT_0315 = (-0.5)*IT_0306 + 0.5*IT_0308 + (-0.5)*IT_0311 +
       0.5*IT_0314;
    const ccomplex_t IT_0316 = IT_0304*IT_0315;
    const ccomplex_t IT_0317 = cpow((-2)*s_12 + (-2)*IT_0032 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0318 = IT_0316*IT_0317;
    const ccomplex_t IT_0319 = (0 + _Complex_I*1)*IT_0318;
    const ccomplex_t IT_0320 = e_em*m_tau*IT_0009*IT_0016;
    const ccomplex_t IT_0321 = 0.5*IT_0320;
    const ccomplex_t IT_0322 = IT_0009*IT_0014*IT_0249;
    const ccomplex_t IT_0323 = IT_0305*IT_0322;
    const ccomplex_t IT_0324 = IT_0307*IT_0322;
    const ccomplex_t IT_0325 = IT_0013*IT_0309;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = IT_0013*IT_0312;
    const ccomplex_t IT_0328 = 1.4142135623731*IT_0327;
    const ccomplex_t IT_0329 = (-0.5)*IT_0323 + 0.5*IT_0324 + 0.5*IT_0326 + (
      -0.5)*IT_0328;
    const ccomplex_t IT_0330 = IT_0321*IT_0329;
    const ccomplex_t IT_0331 = pow(m_Z, 2);
    const ccomplex_t IT_0332 = cpow((-2)*s_12 + (-2)*IT_0032 + IT_0331 + 
      -reg_prop, -1);
    const ccomplex_t IT_0333 = IT_0330*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*IT_0333;
    const ccomplex_t IT_0335 = -IT_0319 + -IT_0334;
    const ccomplex_t IT_0336 = s_13 + s_14;
    const ccomplex_t IT_0337 = m_tau*IT_0336;
    const ccomplex_t IT_0338 = (-2)*conj(IT_0185);
    const ccomplex_t IT_0339 = (-2)*IT_0284*IT_0292;
    const ccomplex_t IT_0340 = -conj(IT_0295);
    const ccomplex_t IT_0341 = s_34 + IT_0000;
    const ccomplex_t IT_0342 = 0.5*conj(IT_0242);
    const ccomplex_t IT_0343 = conj(IT_0335)*IT_0341;
    const ccomplex_t IT_0344 = 2*IT_0002*((IT_0185 + IT_0242)*conj(IT_0272) +
       2*IT_0272*(conj(IT_0272) + 0.5*IT_0273)) + 2*IT_0185*(conj(IT_0275)
      *IT_0276 + 0.5*conj(IT_0185)*IT_0277) + (-8)*IT_0284*(IT_0293 + (-0.125)*
      (conj(IT_0275) + conj(IT_0295))*IT_0296) + conj(IT_0275)*IT_0297*IT_0298 +
       2*IT_0301*(conj(IT_0272)*(IT_0275 + IT_0295) + IT_0272*(conj(IT_0275) +
       conj(IT_0295)) + (-2)*IT_0292*(IT_0272 + -conj(IT_0272) + (-0.5)*IT_0273 
      + (-0.25)*IT_0298 + (-0.25)*IT_0302)) + 2*IT_0295*(conj(IT_0242)*IT_0276 +
       -conj(IT_0275)*IT_0282 + 2*IT_0281*conj(IT_0295) + conj(IT_0335)*IT_0337 
      + 0.5*IT_0297*IT_0338 + IT_0339) + 2*IT_0275*(conj(IT_0185)*IT_0276 + 2
      *conj(IT_0275)*IT_0281 + -conj(IT_0242)*IT_0297 + -conj(IT_0335)*IT_0337 +
       IT_0339 + IT_0282*IT_0340) + (-2)*IT_0335*(IT_0337*(conj(IT_0275) +
       IT_0340) + (-2)*IT_0341*(conj(IT_0335) + 0.25*IT_0338 + IT_0342)) +
       IT_0302*(IT_0000*conj(IT_0242) + conj(IT_0295)*IT_0297 + IT_0343) + 2
      *IT_0242*(IT_0276*conj(IT_0295) + 0.5*IT_0000*IT_0338 + IT_0277*IT_0342 +
       IT_0343);
    return create_ccomplex_return(IT_0344);
}

