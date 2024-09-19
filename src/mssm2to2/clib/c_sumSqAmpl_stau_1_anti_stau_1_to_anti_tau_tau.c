#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_stau_1_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_stau_1_to_anti_tau_tau(
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
    const creal_t m_stau_1 = param->m_stau_1;
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
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = cos(beta);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0006*IT_0011;
    const ccomplex_t IT_0013 = pow(m_W, -1);
    const ccomplex_t IT_0014 = N_d2*e_em*m_tau*IT_0013*conj(U_stau_10);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0004 + IT_0009 + 
      -IT_0016);
    const ccomplex_t IT_0018 = (-0.5)*IT_0017;
    const ccomplex_t IT_0019 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0020 = IT_0001*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = N_d2*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0021 + 0.5*IT_0024);
    const ccomplex_t IT_0026 = -IT_0025;
    const ccomplex_t IT_0027 = IT_0018*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0029 = pow(m_tau, 2);
    const ccomplex_t IT_0030 = pow(m_stau_1, 2);
    const ccomplex_t IT_0031 = cpow(s_23 + (-0.5)*IT_0029 + (-0.5)*IT_0030 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0027*IT_0028*IT_0031;
    const ccomplex_t IT_0033 = N_B3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0034 = IT_0001*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = N_W3*e_em*conj(U_stau_00);
    const ccomplex_t IT_0037 = IT_0006*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = N_d3*e_em*m_tau*IT_0013*conj(U_stau_10);
    const ccomplex_t IT_0040 = IT_0012*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0035 + IT_0038 + 
      -IT_0041);
    const ccomplex_t IT_0043 = (-0.5)*IT_0042;
    const ccomplex_t IT_0044 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0045 = IT_0001*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_d3*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0048 = IT_0012*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + 0.5*IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0054 = cpow(s_23 + (-0.5)*IT_0029 + (-0.5)*IT_0030 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0055 = IT_0052*IT_0053*IT_0054;
    const ccomplex_t IT_0056 = N_B4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0057 = IT_0001*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = N_W4*e_em*conj(U_stau_00);
    const ccomplex_t IT_0060 = IT_0006*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = N_d4*e_em*m_tau*IT_0013*conj(U_stau_10);
    const ccomplex_t IT_0063 = IT_0012*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0058 + IT_0061 + 
      -IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0068 = IT_0001*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = N_d4*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0071 = IT_0012*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 0.5*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = IT_0066*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0077 = cpow(s_23 + (-0.5)*IT_0029 + (-0.5)*IT_0030 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0078 = IT_0075*IT_0076*IT_0077;
    const ccomplex_t IT_0079 = IT_0001*IT_0005;
    const ccomplex_t IT_0080 = e_em*IT_0079;
    const ccomplex_t IT_0081 = IT_0000*IT_0006;
    const ccomplex_t IT_0082 = e_em*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0080 + -IT_0082);
    const ccomplex_t IT_0084 = 0.5*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0030 + -reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*e_em*(0.5*(IT_0001*IT_0005 +
       -IT_0000*IT_0006)*U_stau_00*conj(U_stau_00) + IT_0001*IT_0005*U_stau_10
      *conj(U_stau_10));
    const ccomplex_t IT_0089 = -IT_0088;
    const ccomplex_t IT_0090 = m_tau*IT_0089;
    const ccomplex_t IT_0091 = (-2)*IT_0088;
    const ccomplex_t IT_0092 = m_tau*IT_0091;
    const ccomplex_t IT_0093 = -IT_0092;
    const ccomplex_t IT_0094 = IT_0090 + IT_0093;
    const ccomplex_t IT_0095 = IT_0087*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0005;
    const ccomplex_t IT_0097 = IT_0085*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*IT_0097;
    const ccomplex_t IT_0099 = IT_0094*IT_0098;
    const ccomplex_t IT_0100 = cpow(s_23 + (-0.5)*IT_0029 + (-0.5)*IT_0030 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0101 = N_B1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0102 = IT_0001*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_W1*e_em*conj(U_stau_00);
    const ccomplex_t IT_0105 = IT_0006*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = N_d1*e_em*m_tau*IT_0013*conj(U_stau_10);
    const ccomplex_t IT_0108 = IT_0012*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0103 + IT_0106 + 
      -IT_0109);
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0113 = IT_0001*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0116 = IT_0006*IT_0115;
    const ccomplex_t IT_0117 = 1.4142135623731*IT_0116;
    const ccomplex_t IT_0118 = conj(N_d1)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0119 = IT_0012*IT_0118;
    const ccomplex_t IT_0120 = 1.4142135623731*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*(IT_0114 + IT_0117 + 
      -IT_0120);
    const ccomplex_t IT_0122 = (-0.5)*IT_0121;
    const ccomplex_t IT_0123 = IT_0111*IT_0122;
    const ccomplex_t IT_0124 = IT_0100*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = m_tau*IT_0125;
    const ccomplex_t IT_0127 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0128 = IT_0001*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0131 = IT_0006*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = conj(N_d3)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0134 = IT_0012*IT_0133;
    const ccomplex_t IT_0135 = 1.4142135623731*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*(IT_0129 + IT_0132 + 
      -IT_0135);
    const ccomplex_t IT_0137 = (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = IT_0043*IT_0137;
    const ccomplex_t IT_0139 = IT_0054*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*IT_0139;
    const ccomplex_t IT_0141 = m_tau*IT_0140;
    const ccomplex_t IT_0142 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0143 = IT_0001*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0146 = IT_0006*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = conj(N_d4)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0149 = IT_0012*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0144 + IT_0147 + 
      -IT_0150);
    const ccomplex_t IT_0152 = (-0.5)*IT_0151;
    const ccomplex_t IT_0153 = IT_0066*IT_0152;
    const ccomplex_t IT_0154 = IT_0077*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = m_tau*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0158 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0159 = IT_0001*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = N_d1*e_em*m_tau*IT_0013*U_stau_00;
    const ccomplex_t IT_0162 = IT_0012*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + 0.5*IT_0163);
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = IT_0111*IT_0165;
    const ccomplex_t IT_0167 = IT_0100*IT_0157*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0169 = IT_0001*IT_0168;
    const ccomplex_t IT_0170 = 1.4142135623731*IT_0169;
    const ccomplex_t IT_0171 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0172 = IT_0006*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = conj(N_d2)*e_em*m_tau*IT_0013*U_stau_10;
    const ccomplex_t IT_0175 = IT_0012*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0170 + IT_0173 + 
      -IT_0176);
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = IT_0018*IT_0178;
    const ccomplex_t IT_0180 = IT_0031*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = m_tau*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0032 + 0.5*IT_0055 + 0.5*IT_0078 + 
      -IT_0095 + IT_0099 + (-0.5)*IT_0126 + (-0.5)*IT_0141 + (-0.5)*IT_0156 +
       0.5*IT_0167 + (-0.5)*IT_0182;
    const ccomplex_t IT_0184 = IT_0087*IT_0091;
    const ccomplex_t IT_0185 = 0.5*IT_0125 + 0.5*IT_0140 + 0.5*IT_0155 + 0.5
      *IT_0181 + -IT_0184;
    const ccomplex_t IT_0186 = s_13*m_tau;
    const ccomplex_t IT_0187 = 2*s_34;
    const ccomplex_t IT_0188 = s_13*s_14;
    const ccomplex_t IT_0189 = s_34*IT_0030;
    const ccomplex_t IT_0190 = (-0.5)*IT_0189;
    const ccomplex_t IT_0191 = IT_0188 + IT_0190;
    const ccomplex_t IT_0192 = IT_0029*IT_0030;
    const ccomplex_t IT_0193 = (-0.5)*IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*e_em;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = 2*IT_0195;
    const ccomplex_t IT_0198 = IT_0196*IT_0197;
    const ccomplex_t IT_0199 = cpow(s_12 + IT_0030 + 0.5*reg_prop, -1);
    const ccomplex_t IT_0200 = IT_0198*IT_0199;
    const ccomplex_t IT_0201 = (0 + _Complex_I*1)*IT_0200;
    const ccomplex_t IT_0202 = 0.5*IT_0201;
    const ccomplex_t IT_0203 = cpow(IT_0202, 2);
    const ccomplex_t IT_0204 = conj(N_B3)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0205 = IT_0001*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = conj(N_d3)*e_em*m_tau*IT_0013*conj(U_stau_00);
    const ccomplex_t IT_0208 = IT_0012*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(IT_0206 + 0.5*IT_0209);
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0051*IT_0211;
    const ccomplex_t IT_0213 = IT_0054*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = conj(N_B4)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0216 = IT_0001*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = conj(N_d4)*e_em*m_tau*IT_0013*conj(U_stau_00);
    const ccomplex_t IT_0219 = IT_0012*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0217 + 0.5*IT_0220);
    const ccomplex_t IT_0222 = -IT_0221;
    const ccomplex_t IT_0223 = IT_0074*IT_0222;
    const ccomplex_t IT_0224 = IT_0077*IT_0223;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = IT_0091*IT_0098;
    const ccomplex_t IT_0227 = conj(N_B1)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0228 = IT_0001*IT_0227;
    const ccomplex_t IT_0229 = 1.4142135623731*IT_0228;
    const ccomplex_t IT_0230 = conj(N_d1)*e_em*m_tau*IT_0013*conj(U_stau_00);
    const ccomplex_t IT_0231 = IT_0012*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0229 + 0.5*IT_0232);
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0165*IT_0234;
    const ccomplex_t IT_0236 = IT_0100*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0239 = IT_0001*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = conj(N_d2)*e_em*m_tau*IT_0013*conj(U_stau_00);
    const ccomplex_t IT_0242 = IT_0012*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*(IT_0240 + 0.5*IT_0243);
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = IT_0026*IT_0245;
    const ccomplex_t IT_0247 = IT_0031*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = 0.5*IT_0214 + 0.5*IT_0225 + -IT_0226 + 0.5
      *IT_0237 + 0.5*IT_0248;
    const ccomplex_t IT_0250 = 4*IT_0202;
    const ccomplex_t IT_0251 = s_14*m_tau;
    const ccomplex_t IT_0252 = IT_0178*IT_0245;
    const ccomplex_t IT_0253 = IT_0028*IT_0031*IT_0252;
    const ccomplex_t IT_0254 = IT_0137*IT_0211;
    const ccomplex_t IT_0255 = IT_0053*IT_0054*IT_0254;
    const ccomplex_t IT_0256 = IT_0152*IT_0222;
    const ccomplex_t IT_0257 = IT_0076*IT_0077*IT_0256;
    const ccomplex_t IT_0258 = m_tau*IT_0248;
    const ccomplex_t IT_0259 = m_tau*IT_0214;
    const ccomplex_t IT_0260 = m_tau*IT_0225;
    const ccomplex_t IT_0261 = IT_0122*IT_0234;
    const ccomplex_t IT_0262 = IT_0100*IT_0157*IT_0261;
    const ccomplex_t IT_0263 = m_tau*IT_0237;
    const ccomplex_t IT_0264 = IT_0095 + -IT_0099 + 0.5*IT_0253 + 0.5*IT_0255 
      + 0.5*IT_0257 + (-0.5)*IT_0258 + (-0.5)*IT_0259 + (-0.5)*IT_0260 + 0.5
      *IT_0262 + (-0.5)*IT_0263;
    const ccomplex_t IT_0265 = (-2)*IT_0264;
    const ccomplex_t IT_0266 = -s_14;
    const ccomplex_t IT_0267 = s_13 + IT_0266;
    const ccomplex_t IT_0268 = m_tau*IT_0267;
    const ccomplex_t IT_0269 = sin(alpha);
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0011
      *IT_0013*IT_0269;
    const ccomplex_t IT_0271 = 0.5*IT_0270;
    const ccomplex_t IT_0272 = mu_h*m_tau*IT_0013;
    const ccomplex_t IT_0273 = cos(alpha);
    const ccomplex_t IT_0274 = IT_0010*IT_0269;
    const ccomplex_t IT_0275 = sin(beta);
    const ccomplex_t IT_0276 = IT_0273*IT_0275;
    const ccomplex_t IT_0277 = IT_0274 + IT_0276;
    const ccomplex_t IT_0278 = cpow(IT_0000, -2);
    const ccomplex_t IT_0279 = 2*U_stau_00*((IT_0006*IT_0011*IT_0013*IT_0029
      *IT_0269 + (-0.5)*m_W*IT_0277*(IT_0006 + -IT_0005*IT_0278))*conj(U_stau_00
      ) + 1./2*IT_0006*IT_0011*IT_0272*IT_0273*conj(U_stau_10)) + U_stau_10*
      (IT_0006*IT_0011*IT_0272*IT_0273*conj(U_stau_00) + 2*(IT_0006*IT_0011
      *IT_0013*IT_0029*IT_0269 + -m_W*IT_0005*IT_0277*IT_0278)*conj(U_stau_10));
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*(e_em*IT_0279 +
       1.4142135623731*A_tau*IT_0269*(conj(U_stau_00)*U_stau_10 + U_stau_00*conj
      (U_stau_10)));
    const ccomplex_t IT_0281 = 0.5*IT_0280;
    const ccomplex_t IT_0282 = IT_0271*IT_0281;
    const ccomplex_t IT_0283 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0030 + -reg_prop, -1);
    const ccomplex_t IT_0284 = IT_0282*IT_0283;
    const ccomplex_t IT_0285 = (0 + _Complex_I*1)*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*e_em*m_tau*IT_0006*IT_0011
      *IT_0013*IT_0273;
    const ccomplex_t IT_0287 = (-0.5)*IT_0286;
    const ccomplex_t IT_0288 = -IT_0269*IT_0275;
    const ccomplex_t IT_0289 = IT_0010*IT_0273;
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = (-2)*U_stau_00*((IT_0006*IT_0011*IT_0013
      *IT_0029*IT_0273 + (-0.5)*m_W*(IT_0006 + -IT_0005*IT_0278)*IT_0290)*conj
      (U_stau_00) + -1./2*IT_0006*IT_0011*IT_0269*IT_0272*conj(U_stau_10)) +
       U_stau_10*(IT_0006*IT_0011*IT_0269*IT_0272*conj(U_stau_00) + (-2)*
      (IT_0006*IT_0011*IT_0013*IT_0029*IT_0273 + -m_W*IT_0005*IT_0278*IT_0290)
      *conj(U_stau_10));
    const ccomplex_t IT_0292 = (0 + _Complex_I*-1)*(e_em*IT_0291 + (
      -1.4142135623731)*A_tau*IT_0273*(conj(U_stau_00)*U_stau_10 + U_stau_00
      *conj(U_stau_10)));
    const ccomplex_t IT_0293 = (-0.5)*IT_0292;
    const ccomplex_t IT_0294 = IT_0287*IT_0293;
    const ccomplex_t IT_0295 = cpow((-2)*s_12 + (-2)*IT_0030 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0296 = IT_0294*IT_0295;
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*IT_0296;
    const ccomplex_t IT_0298 = -IT_0285 + -IT_0297;
    const ccomplex_t IT_0299 = (-2)*IT_0183;
    const ccomplex_t IT_0300 = 2*IT_0268*((conj(IT_0185) + conj(IT_0249))
      *IT_0298 + (IT_0185 + IT_0249)*conj(IT_0298) + (-2)*IT_0202*((-0.5)*conj
      (IT_0183) + (-0.5)*conj(IT_0264) + (-0.25)*IT_0265 + IT_0298 + -conj
      (IT_0298) + (-0.25)*IT_0299));
    const ccomplex_t IT_0301 = e_em*m_tau*IT_0006*IT_0013;
    const ccomplex_t IT_0302 = 0.5*IT_0301;
    const ccomplex_t IT_0303 = e_em*mu_h*m_tau*IT_0013*conj(U_stau_00)
      *U_stau_10;
    const ccomplex_t IT_0304 = IT_0006*IT_0011*IT_0275;
    const ccomplex_t IT_0305 = IT_0303*IT_0304;
    const ccomplex_t IT_0306 = e_em*mu_h*m_tau*IT_0013*U_stau_00*conj
      (U_stau_10);
    const ccomplex_t IT_0307 = IT_0304*IT_0306;
    const ccomplex_t IT_0308 = A_tau*conj(U_stau_00)*U_stau_10;
    const ccomplex_t IT_0309 = IT_0010*IT_0308;
    const ccomplex_t IT_0310 = 1.4142135623731*IT_0309;
    const ccomplex_t IT_0311 = A_tau*U_stau_00*conj(U_stau_10);
    const ccomplex_t IT_0312 = IT_0010*IT_0311;
    const ccomplex_t IT_0313 = 1.4142135623731*IT_0312;
    const ccomplex_t IT_0314 = (-0.5)*IT_0305 + 0.5*IT_0307 + 0.5*IT_0310 + (
      -0.5)*IT_0313;
    const ccomplex_t IT_0315 = IT_0302*IT_0314;
    const ccomplex_t IT_0316 = pow(m_Z, 2);
    const ccomplex_t IT_0317 = cpow((-2)*s_12 + (-2)*IT_0030 + IT_0316 + 
      -reg_prop, -1);
    const ccomplex_t IT_0318 = IT_0315*IT_0317;
    const ccomplex_t IT_0319 = (0 + _Complex_I*1)*IT_0318;
    const ccomplex_t IT_0320 = e_em*m_tau*IT_0006*IT_0011*IT_0013*IT_0275;
    const ccomplex_t IT_0321 = (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = IT_0006*IT_0303;
    const ccomplex_t IT_0323 = IT_0006*IT_0306;
    const ccomplex_t IT_0324 = IT_0275*IT_0308;
    const ccomplex_t IT_0325 = 1.4142135623731*IT_0324;
    const ccomplex_t IT_0326 = IT_0275*IT_0311;
    const ccomplex_t IT_0327 = 1.4142135623731*IT_0326;
    const ccomplex_t IT_0328 = (-0.5)*IT_0322 + 0.5*IT_0323 + (-0.5)*IT_0325 +
       0.5*IT_0327;
    const ccomplex_t IT_0329 = IT_0321*IT_0328;
    const ccomplex_t IT_0330 = cpow((-2)*s_12 + (-2)*IT_0030 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0331 = IT_0329*IT_0330;
    const ccomplex_t IT_0332 = (0 + _Complex_I*1)*IT_0331;
    const ccomplex_t IT_0333 = -IT_0319 + -IT_0332;
    const ccomplex_t IT_0334 = s_34 + IT_0029;
    const ccomplex_t IT_0335 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0336 = s_13 + s_14;
    const ccomplex_t IT_0337 = m_tau*IT_0336;
    const ccomplex_t IT_0338 = -conj(IT_0249);
    const ccomplex_t IT_0339 = -IT_0029;
    const ccomplex_t IT_0340 = s_34 + IT_0339;
    const ccomplex_t IT_0341 = 0.5*conj(IT_0264);
    const ccomplex_t IT_0342 = (-4)*IT_0202;
    const ccomplex_t IT_0343 = 0.5*IT_0194*IT_0342;
    const ccomplex_t IT_0344 = conj(IT_0333)*IT_0334;
    const ccomplex_t IT_0345 = 2*IT_0183*(conj(IT_0185)*IT_0186 + 0.5*conj
      (IT_0183)*IT_0187) + (-8)*IT_0194*(IT_0203 + (-0.125)*(conj(IT_0185) +
       conj(IT_0249))*IT_0250) + conj(IT_0185)*IT_0251*IT_0265 + IT_0300 + 2
      *IT_0333*(IT_0334*(conj(IT_0264) + 2*conj(IT_0333) + 0.5*IT_0335) + 
      -IT_0337*(conj(IT_0185) + IT_0338)) + 2*IT_0340*((IT_0183 + IT_0264)*conj
      (IT_0298) + 2*IT_0298*(0.5*conj(IT_0183) + conj(IT_0298) + IT_0341)) + 2
      *IT_0185*(conj(IT_0183)*IT_0186 + 2*conj(IT_0185)*IT_0191 + -IT_0251*conj
      (IT_0264) + -conj(IT_0333)*IT_0337 + IT_0192*IT_0338 + IT_0343) + (-2)
      *IT_0249*(conj(IT_0185)*IT_0192 + (-2)*IT_0191*conj(IT_0249) + -IT_0186
      *conj(IT_0264) + (-0.5)*IT_0251*IT_0335 + -conj(IT_0333)*IT_0337 + 
      -IT_0343) + IT_0299*(conj(IT_0249)*IT_0251 + IT_0029*conj(IT_0264) +
       IT_0344) + 2*IT_0264*(IT_0186*conj(IT_0249) + 0.5*IT_0029*IT_0335 +
       IT_0187*IT_0341 + IT_0344);
    return create_ccomplex_return(IT_0345);
}

