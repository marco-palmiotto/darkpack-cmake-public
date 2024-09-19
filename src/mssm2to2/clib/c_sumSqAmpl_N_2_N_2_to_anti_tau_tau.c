#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_2_to_anti_tau_tau.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_2_to_anti_tau_tau(
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
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_tau = param->m_tau;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t m_stau_1 = param->m_stau_1;
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_tau1 = param->Gamma_tau1;
    const creal_t Gamma_tau2 = param->Gamma_tau2;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = m_N_2*m_tau;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = (-2)*IT_0003;
    const ccomplex_t IT_0005 = pow(m_N_2, 2);
    const ccomplex_t IT_0006 = pow(m_tau, 2);
    const ccomplex_t IT_0007 = cpow(s_23 + (-0.5)*IT_0005 + (-0.5)*IT_0006 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_B2)*e_em*conj(U_stau_10);
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0014*IT_0016;
    const ccomplex_t IT_0018 = pow(m_W, -1);
    const ccomplex_t IT_0019 = conj(N_d2)*e_em*m_tau*IT_0018*conj(U_stau_00);
    const ccomplex_t IT_0020 = IT_0017*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*(IT_0012 + 0.5*IT_0021);
    const ccomplex_t IT_0023 = -IT_0022;
    const ccomplex_t IT_0024 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0025 = IT_0009*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = N_d2*e_em*m_tau*IT_0018*U_stau_00;
    const ccomplex_t IT_0028 = IT_0017*IT_0027;
    const ccomplex_t IT_0029 = 1.4142135623731*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*(IT_0026 + 0.5*IT_0029);
    const ccomplex_t IT_0031 = -IT_0030;
    const ccomplex_t IT_0032 = IT_0023*IT_0031;
    const ccomplex_t IT_0033 = IT_0007*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*IT_0033;
    const ccomplex_t IT_0035 = cpow(s_23 + (-0.5)*IT_0005 + (-0.5)*IT_0006 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0036 = conj(N_B2)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0037 = IT_0009*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = conj(N_d2)*e_em*m_tau*IT_0018*conj(U_stau_01);
    const ccomplex_t IT_0040 = IT_0017*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*(IT_0038 + 0.5*IT_0041);
    const ccomplex_t IT_0043 = -IT_0042;
    const ccomplex_t IT_0044 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0045 = IT_0009*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = N_d2*e_em*m_tau*IT_0018*U_stau_01;
    const ccomplex_t IT_0048 = IT_0017*IT_0047;
    const ccomplex_t IT_0049 = 1.4142135623731*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*(IT_0046 + 0.5*IT_0049);
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = IT_0043*IT_0051;
    const ccomplex_t IT_0053 = IT_0035*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*IT_0053;
    const ccomplex_t IT_0055 = (-0.5)*IT_0034 + (-0.5)*IT_0054;
    const ccomplex_t IT_0056 = cos(alpha);
    const ccomplex_t IT_0057 = IT_0009*IT_0056;
    const ccomplex_t IT_0058 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = sin(alpha);
    const ccomplex_t IT_0061 = IT_0009*IT_0060;
    const ccomplex_t IT_0062 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = IT_0016*IT_0056;
    const ccomplex_t IT_0065 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = IT_0016*IT_0060;
    const ccomplex_t IT_0068 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0059 + IT_0063 + 
      -IT_0066 + -IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*e_em*m_tau*IT_0014*IT_0016
      *IT_0018*IT_0060;
    const ccomplex_t IT_0073 = 0.5*IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0071*IT_0076;
    const ccomplex_t IT_0078 = IT_0057*IT_0062;
    const ccomplex_t IT_0079 = IT_0058*IT_0061;
    const ccomplex_t IT_0080 = IT_0064*IT_0068;
    const ccomplex_t IT_0081 = IT_0065*IT_0067;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0078 + -IT_0079 + 
      -IT_0080 + IT_0081);
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*e_em*m_tau*IT_0014*IT_0016
      *IT_0018*IT_0056;
    const ccomplex_t IT_0084 = (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = cpow((-2)*s_12 + (-2)*IT_0005 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0086 = IT_0084*IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = IT_0082*IT_0087;
    const ccomplex_t IT_0089 = 0.5*IT_0077 + 0.5*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0077 + (-0.5)*IT_0088;
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = -conj(IT_0090);
    const ccomplex_t IT_0094 = conj(IT_0089) + IT_0093;
    const ccomplex_t IT_0095 = N_B2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0096 = IT_0009*IT_0095;
    const ccomplex_t IT_0097 = 1.4142135623731*IT_0096;
    const ccomplex_t IT_0098 = N_W2*e_em*conj(U_stau_00);
    const ccomplex_t IT_0099 = IT_0016*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = N_d2*e_em*m_tau*IT_0018*conj(U_stau_10);
    const ccomplex_t IT_0102 = IT_0017*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0097 + IT_0100 + 
      -IT_0103);
    const ccomplex_t IT_0105 = (-0.5)*IT_0104;
    const ccomplex_t IT_0106 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0107 = IT_0009*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0110 = IT_0016*IT_0109;
    const ccomplex_t IT_0111 = 1.4142135623731*IT_0110;
    const ccomplex_t IT_0112 = conj(N_d2)*e_em*m_tau*IT_0018*U_stau_10;
    const ccomplex_t IT_0113 = IT_0017*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0108 + IT_0111 + 
      -IT_0114);
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = IT_0105*IT_0116;
    const ccomplex_t IT_0118 = IT_0007*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = N_B2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0121 = IT_0009*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = N_W2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0124 = IT_0016*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = N_d2*e_em*m_tau*IT_0018*conj(U_stau_11);
    const ccomplex_t IT_0127 = IT_0017*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0122 + IT_0125 + 
      -IT_0128);
    const ccomplex_t IT_0130 = (-0.5)*IT_0129;
    const ccomplex_t IT_0131 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0132 = IT_0009*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0135 = IT_0016*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = conj(N_d2)*e_em*m_tau*IT_0018*U_stau_11;
    const ccomplex_t IT_0138 = IT_0017*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0133 + IT_0136 + 
      -IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = IT_0130*IT_0141;
    const ccomplex_t IT_0143 = IT_0035*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = (-0.5)*IT_0119 + (-0.5)*IT_0144;
    const ccomplex_t IT_0146 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0147 = IT_0057*IT_0146;
    const ccomplex_t IT_0148 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0149 = IT_0061*IT_0148;
    const ccomplex_t IT_0150 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0151 = IT_0064*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0153 = IT_0067*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*(IT_0147 + IT_0149 + 
      -IT_0151 + -IT_0153);
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = IT_0076*IT_0155;
    const ccomplex_t IT_0157 = IT_0057*IT_0148;
    const ccomplex_t IT_0158 = IT_0061*IT_0146;
    const ccomplex_t IT_0159 = IT_0064*IT_0152;
    const ccomplex_t IT_0160 = IT_0067*IT_0150;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*(IT_0157 + -IT_0158 + 
      -IT_0159 + IT_0160);
    const ccomplex_t IT_0162 = IT_0087*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0156 + 0.5*IT_0162;
    const ccomplex_t IT_0164 = (-0.5)*IT_0156 + (-0.5)*IT_0162;
    const ccomplex_t IT_0165 = -IT_0164;
    const ccomplex_t IT_0166 = IT_0163 + IT_0165;
    const ccomplex_t IT_0167 = -conj(IT_0164);
    const ccomplex_t IT_0168 = conj(IT_0163) + IT_0167;
    const ccomplex_t IT_0169 = -s_24;
    const ccomplex_t IT_0170 = s_23 + IT_0169;
    const ccomplex_t IT_0171 = IT_0000*IT_0170;
    const ccomplex_t IT_0172 = 2*IT_0171;
    const ccomplex_t IT_0173 = (-2)*IT_0171;
    const ccomplex_t IT_0174 = cpow(s_13 + (-0.5)*IT_0005 + (-0.5)*IT_0006 + (
      -0.5)*reg_prop + 0.5*m_stau_1*(m_stau_1 + (0 + _Complex_I*-1)*Gamma_tau1),
       -1);
    const ccomplex_t IT_0175 = IT_0032*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = cpow(s_13 + (-0.5)*IT_0005 + (-0.5)*IT_0006 + (
      -0.5)*reg_prop + 0.5*m_stau_2*(m_stau_2 + (0 + _Complex_I*-1)*Gamma_tau2),
       -1);
    const ccomplex_t IT_0178 = IT_0052*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*IT_0178;
    const ccomplex_t IT_0180 = 0.5*IT_0176 + 0.5*IT_0179;
    const ccomplex_t IT_0181 = -IT_0089;
    const ccomplex_t IT_0182 = IT_0090 + IT_0181;
    const ccomplex_t IT_0183 = -conj(IT_0089);
    const ccomplex_t IT_0184 = conj(IT_0090) + IT_0183;
    const ccomplex_t IT_0185 = IT_0117*IT_0174;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = IT_0142*IT_0177;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*IT_0187;
    const ccomplex_t IT_0189 = 0.5*IT_0186 + 0.5*IT_0188;
    const ccomplex_t IT_0190 = -IT_0163;
    const ccomplex_t IT_0191 = IT_0164 + IT_0190;
    const ccomplex_t IT_0192 = -conj(IT_0163);
    const ccomplex_t IT_0193 = conj(IT_0164) + IT_0192;
    const ccomplex_t IT_0194 = 2*IT_0003;
    const ccomplex_t IT_0195 = 4*IT_0171;
    const ccomplex_t IT_0196 = IT_0008*IT_0016;
    const ccomplex_t IT_0197 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0198 = IT_0196*IT_0197;
    const ccomplex_t IT_0199 = IT_0009*IT_0015;
    const ccomplex_t IT_0200 = IT_0197*IT_0199;
    const ccomplex_t IT_0201 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0202 = IT_0199*IT_0201;
    const ccomplex_t IT_0203 = IT_0196*IT_0201;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*(IT_0198 + IT_0200 + 
      -IT_0202 + -IT_0203);
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = e_em*IT_0199;
    const ccomplex_t IT_0207 = e_em*IT_0196;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*(IT_0206 + -IT_0207);
    const ccomplex_t IT_0209 = 0.5*IT_0208;
    const ccomplex_t IT_0210 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0005 + -reg_prop, -1);
    const ccomplex_t IT_0211 = IT_0209*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = IT_0205*IT_0212;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*e_em*IT_0009*IT_0015;
    const ccomplex_t IT_0216 = IT_0210*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = IT_0205*IT_0217;
    const ccomplex_t IT_0219 = (-0.5)*IT_0218;
    const ccomplex_t IT_0220 = IT_0214 + IT_0219;
    const ccomplex_t IT_0221 = conj(IT_0214) + conj(IT_0219);
    const ccomplex_t IT_0222 = 0.5*IT_0204;
    const ccomplex_t IT_0223 = IT_0212*IT_0222;
    const ccomplex_t IT_0224 = (-0.5)*IT_0223;
    const ccomplex_t IT_0225 = IT_0217*IT_0222;
    const ccomplex_t IT_0226 = (-0.5)*IT_0225;
    const ccomplex_t IT_0227 = IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = conj(IT_0224) + conj(IT_0226);
    const ccomplex_t IT_0229 = 4*IT_0003;
    const ccomplex_t IT_0230 = 0.5*IT_0213;
    const ccomplex_t IT_0231 = 0.5*IT_0218;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = conj(IT_0230) + conj(IT_0231);
    const ccomplex_t IT_0234 = 0.5*IT_0223;
    const ccomplex_t IT_0235 = 0.5*IT_0225;
    const ccomplex_t IT_0236 = IT_0234 + IT_0235;
    const ccomplex_t IT_0237 = conj(IT_0234) + conj(IT_0235);
    const ccomplex_t IT_0238 = s_12*IT_0006;
    const ccomplex_t IT_0239 = s_14*s_23;
    const ccomplex_t IT_0240 = s_13*s_24;
    const ccomplex_t IT_0241 = -IT_0240;
    const ccomplex_t IT_0242 = IT_0239 + IT_0241;
    const ccomplex_t IT_0243 = s_12*s_34;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = IT_0242 + IT_0244;
    const ccomplex_t IT_0246 = IT_0238 + IT_0245;
    const ccomplex_t IT_0247 = IT_0051*IT_0130;
    const ccomplex_t IT_0248 = IT_0177*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = IT_0031*IT_0105;
    const ccomplex_t IT_0251 = IT_0174*IT_0250;
    const ccomplex_t IT_0252 = (0 + _Complex_I*1)*IT_0251;
    const ccomplex_t IT_0253 = (-0.5)*IT_0249 + (-0.5)*IT_0252;
    const ccomplex_t IT_0254 = 2*IT_0089;
    const ccomplex_t IT_0255 = (-2)*IT_0090;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 2*conj(IT_0089);
    const ccomplex_t IT_0258 = (-2)*conj(IT_0090);
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = IT_0023*IT_0116;
    const ccomplex_t IT_0261 = IT_0174*IT_0260;
    const ccomplex_t IT_0262 = (0 + _Complex_I*1)*IT_0261;
    const ccomplex_t IT_0263 = IT_0043*IT_0141;
    const ccomplex_t IT_0264 = IT_0177*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = (-0.5)*IT_0262 + (-0.5)*IT_0265;
    const ccomplex_t IT_0267 = 2*IT_0163;
    const ccomplex_t IT_0268 = (-2)*IT_0164;
    const ccomplex_t IT_0269 = IT_0267 + IT_0268;
    const ccomplex_t IT_0270 = 2*conj(IT_0163);
    const ccomplex_t IT_0271 = (-2)*conj(IT_0164);
    const ccomplex_t IT_0272 = IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = s_13 + s_14;
    const ccomplex_t IT_0274 = IT_0000*IT_0273;
    const ccomplex_t IT_0275 = IT_0009*IT_0013;
    const ccomplex_t IT_0276 = IT_0146*IT_0275;
    const ccomplex_t IT_0277 = sin(beta);
    const ccomplex_t IT_0278 = IT_0009*IT_0277;
    const ccomplex_t IT_0279 = IT_0148*IT_0278;
    const ccomplex_t IT_0280 = IT_0013*IT_0016;
    const ccomplex_t IT_0281 = IT_0150*IT_0280;
    const ccomplex_t IT_0282 = IT_0016*IT_0277;
    const ccomplex_t IT_0283 = IT_0152*IT_0282;
    const ccomplex_t IT_0284 = IT_0276 + -IT_0279 + -IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = e_em*m_tau*IT_0014*IT_0016*IT_0018*IT_0277;
    const ccomplex_t IT_0286 = (-0.5)*IT_0285;
    const ccomplex_t IT_0287 = cpow((-2)*s_12 + (-2)*IT_0005 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0288 = IT_0286*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0284*IT_0289;
    const ccomplex_t IT_0291 = IT_0148*IT_0275;
    const ccomplex_t IT_0292 = IT_0146*IT_0278;
    const ccomplex_t IT_0293 = IT_0152*IT_0280;
    const ccomplex_t IT_0294 = IT_0150*IT_0282;
    const ccomplex_t IT_0295 = IT_0291 + IT_0292 + -IT_0293 + -IT_0294;
    const ccomplex_t IT_0296 = e_em*m_tau*IT_0016*IT_0018;
    const ccomplex_t IT_0297 = 0.5*IT_0296;
    const ccomplex_t IT_0298 = pow(m_Z, 2);
    const ccomplex_t IT_0299 = cpow((-2)*s_12 + (-2)*IT_0005 + IT_0298 + 
      -reg_prop, -1);
    const ccomplex_t IT_0300 = IT_0297*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*IT_0300;
    const ccomplex_t IT_0302 = IT_0295*IT_0301;
    const ccomplex_t IT_0303 = 0.5*IT_0290 + 0.5*IT_0302;
    const ccomplex_t IT_0304 = 2*IT_0145;
    const ccomplex_t IT_0305 = (-2)*IT_0180;
    const ccomplex_t IT_0306 = IT_0304 + IT_0305;
    const ccomplex_t IT_0307 = 2*conj(IT_0145);
    const ccomplex_t IT_0308 = (-2)*conj(IT_0180);
    const ccomplex_t IT_0309 = IT_0307 + IT_0308;
    const ccomplex_t IT_0310 = IT_0058*IT_0275;
    const ccomplex_t IT_0311 = IT_0062*IT_0278;
    const ccomplex_t IT_0312 = IT_0065*IT_0280;
    const ccomplex_t IT_0313 = IT_0068*IT_0282;
    const ccomplex_t IT_0314 = -IT_0310 + IT_0311 + IT_0312 + -IT_0313;
    const ccomplex_t IT_0315 = IT_0289*IT_0314;
    const ccomplex_t IT_0316 = IT_0062*IT_0275;
    const ccomplex_t IT_0317 = IT_0058*IT_0278;
    const ccomplex_t IT_0318 = IT_0068*IT_0280;
    const ccomplex_t IT_0319 = IT_0065*IT_0282;
    const ccomplex_t IT_0320 = -IT_0316 + -IT_0317 + IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = IT_0301*IT_0320;
    const ccomplex_t IT_0322 = (-0.5)*IT_0315 + (-0.5)*IT_0321;
    const ccomplex_t IT_0323 = 2*IT_0055;
    const ccomplex_t IT_0324 = (-2)*IT_0189;
    const ccomplex_t IT_0325 = IT_0323 + IT_0324;
    const ccomplex_t IT_0326 = 2*conj(IT_0055);
    const ccomplex_t IT_0327 = (-2)*conj(IT_0189);
    const ccomplex_t IT_0328 = IT_0326 + IT_0327;
    const ccomplex_t IT_0329 = 0.5*IT_0315 + 0.5*IT_0321;
    const ccomplex_t IT_0330 = 2*IT_0189;
    const ccomplex_t IT_0331 = (-2)*IT_0055;
    const ccomplex_t IT_0332 = IT_0330 + IT_0331;
    const ccomplex_t IT_0333 = 2*conj(IT_0189);
    const ccomplex_t IT_0334 = (-2)*conj(IT_0055);
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = (-0.5)*IT_0290 + (-0.5)*IT_0302;
    const ccomplex_t IT_0337 = 2*IT_0180;
    const ccomplex_t IT_0338 = (-2)*IT_0145;
    const ccomplex_t IT_0339 = IT_0337 + IT_0338;
    const ccomplex_t IT_0340 = 2*conj(IT_0180);
    const ccomplex_t IT_0341 = (-2)*conj(IT_0145);
    const ccomplex_t IT_0342 = IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = conj(IT_0219) + conj(IT_0235);
    const ccomplex_t IT_0344 = 4*IT_0322;
    const ccomplex_t IT_0345 = (-4)*IT_0329;
    const ccomplex_t IT_0346 = IT_0344 + IT_0345;
    const ccomplex_t IT_0347 = IT_0219 + IT_0235;
    const ccomplex_t IT_0348 = 4*conj(IT_0322);
    const ccomplex_t IT_0349 = (-4)*conj(IT_0329);
    const ccomplex_t IT_0350 = IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = conj(IT_0224) + conj(IT_0230);
    const ccomplex_t IT_0352 = 4*IT_0303;
    const ccomplex_t IT_0353 = (-4)*IT_0336;
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = IT_0224 + IT_0230;
    const ccomplex_t IT_0356 = 4*conj(IT_0303);
    const ccomplex_t IT_0357 = (-4)*conj(IT_0336);
    const ccomplex_t IT_0358 = IT_0356 + IT_0357;
    const ccomplex_t IT_0359 = conj(IT_0226) + conj(IT_0231);
    const ccomplex_t IT_0360 = 4*IT_0336;
    const ccomplex_t IT_0361 = (-4)*IT_0303;
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = IT_0226 + IT_0231;
    const ccomplex_t IT_0364 = 4*conj(IT_0336);
    const ccomplex_t IT_0365 = (-4)*conj(IT_0303);
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = conj(IT_0214) + conj(IT_0234);
    const ccomplex_t IT_0368 = 4*IT_0329;
    const ccomplex_t IT_0369 = (-4)*IT_0322;
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = IT_0214 + IT_0234;
    const ccomplex_t IT_0372 = 4*conj(IT_0329);
    const ccomplex_t IT_0373 = (-4)*conj(IT_0322);
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = conj(IT_0303)*IT_0306 + IT_0303*IT_0309 + conj
      (IT_0322)*IT_0325 + IT_0322*IT_0328 + conj(IT_0329)*IT_0332 + IT_0329
      *IT_0335 + conj(IT_0336)*IT_0339 + IT_0336*IT_0342 + IT_0343*IT_0346 +
       IT_0347*IT_0350 + IT_0351*IT_0354 + IT_0355*IT_0358 + IT_0359*IT_0362 +
       IT_0363*IT_0366 + IT_0367*IT_0370 + IT_0371*IT_0374;
    const ccomplex_t IT_0376 = s_23 + s_24;
    const ccomplex_t IT_0377 = IT_0000*IT_0376;
    const ccomplex_t IT_0378 = IT_0309*IT_0322 + IT_0306*conj(IT_0322) + conj
      (IT_0303)*IT_0325 + IT_0303*IT_0328 + IT_0335*IT_0336 + IT_0332*conj
      (IT_0336) + conj(IT_0329)*IT_0339 + IT_0329*IT_0342 + IT_0346*IT_0351 +
       IT_0343*IT_0354 + IT_0350*IT_0355 + IT_0347*IT_0358 + IT_0362*IT_0367 +
       IT_0359*IT_0370 + IT_0366*IT_0371 + IT_0363*IT_0374;
    const ccomplex_t IT_0379 = IT_0035*IT_0247;
    const ccomplex_t IT_0380 = (0 + _Complex_I*1)*IT_0379;
    const ccomplex_t IT_0381 = IT_0007*IT_0250;
    const ccomplex_t IT_0382 = (0 + _Complex_I*1)*IT_0381;
    const ccomplex_t IT_0383 = 0.5*IT_0380 + 0.5*IT_0382;
    const ccomplex_t IT_0384 = IT_0007*IT_0260;
    const ccomplex_t IT_0385 = (0 + _Complex_I*1)*IT_0384;
    const ccomplex_t IT_0386 = IT_0035*IT_0263;
    const ccomplex_t IT_0387 = (0 + _Complex_I*1)*IT_0386;
    const ccomplex_t IT_0388 = 0.5*IT_0385 + 0.5*IT_0387;
    const ccomplex_t IT_0389 = -IT_0238;
    const ccomplex_t IT_0390 = IT_0239 + IT_0240 + IT_0244;
    const ccomplex_t IT_0391 = 2*IT_0253;
    const ccomplex_t IT_0392 = 2*conj(IT_0253);
    const ccomplex_t IT_0393 = 2*IT_0266;
    const ccomplex_t IT_0394 = 2*conj(IT_0266);
    const ccomplex_t IT_0395 = IT_0238 + IT_0242 + IT_0243;
    const ccomplex_t IT_0396 = 2*IT_0388;
    const ccomplex_t IT_0397 = 2*conj(IT_0388);
    const ccomplex_t IT_0398 = 2*IT_0383;
    const ccomplex_t IT_0399 = 2*conj(IT_0383);
    const ccomplex_t IT_0400 = (-2)*IT_0388;
    const ccomplex_t IT_0401 = (-2)*conj(IT_0388);
    const ccomplex_t IT_0402 = (-2)*IT_0383;
    const ccomplex_t IT_0403 = (-2)*conj(IT_0383);
    const ccomplex_t IT_0404 = IT_0245 + IT_0389;
    const ccomplex_t IT_0405 = (-2)*IT_0253;
    const ccomplex_t IT_0406 = (-2)*conj(IT_0253);
    const ccomplex_t IT_0407 = (-2)*IT_0266;
    const ccomplex_t IT_0408 = (-2)*conj(IT_0266);
    const ccomplex_t IT_0409 = (-4)*IT_0171;
    const ccomplex_t IT_0410 = -IT_0006;
    const ccomplex_t IT_0411 = s_34 + IT_0410;
    const ccomplex_t IT_0412 = s_12*IT_0411;
    const ccomplex_t IT_0413 = (-8)*conj(IT_0089);
    const ccomplex_t IT_0414 = (-4)*IT_0003;
    const ccomplex_t IT_0415 = s_34*IT_0005;
    const ccomplex_t IT_0416 = 4*IT_0145;
    const ccomplex_t IT_0417 = 0.125*IT_0416;
    const ccomplex_t IT_0418 = IT_0355 + IT_0417;
    const ccomplex_t IT_0419 = 4*conj(IT_0145);
    const ccomplex_t IT_0420 = 0.125*IT_0419;
    const ccomplex_t IT_0421 = IT_0351 + IT_0420;
    const ccomplex_t IT_0422 = 4*IT_0055;
    const ccomplex_t IT_0423 = 4*conj(IT_0055);
    const ccomplex_t IT_0424 = 4*IT_0189;
    const ccomplex_t IT_0425 = 4*conj(IT_0189);
    const ccomplex_t IT_0426 = 4*conj(IT_0180);
    const ccomplex_t IT_0427 = 0.125*IT_0426;
    const ccomplex_t IT_0428 = IT_0359 + IT_0427;
    const ccomplex_t IT_0429 = 8*IT_0347;
    const ccomplex_t IT_0430 = 4*IT_0180;
    const ccomplex_t IT_0431 = 0.125*IT_0430;
    const ccomplex_t IT_0432 = IT_0363 + IT_0431;
    const ccomplex_t IT_0433 = 8*IT_0343;
    const ccomplex_t IT_0434 = 0.125*conj(IT_0145)*IT_0330 + 0.125*IT_0145
      *IT_0333 + 0.125*conj(IT_0055)*IT_0337 + 0.125*IT_0055*IT_0340 + IT_0367
      *IT_0418 + IT_0371*IT_0421 + 0.125*IT_0359*IT_0422 + 0.125*IT_0363*IT_0423
       + 0.125*IT_0351*IT_0424 + 0.125*IT_0355*IT_0425 + 0.125*IT_0428*IT_0429 +
       0.125*IT_0432*IT_0433;
    const ccomplex_t IT_0435 = s_14*m_N_2*m_tau;
    const ccomplex_t IT_0436 = 4*IT_0383;
    const ccomplex_t IT_0437 = IT_0391 + IT_0436;
    const ccomplex_t IT_0438 = 4*conj(IT_0383);
    const ccomplex_t IT_0439 = IT_0392 + IT_0438;
    const ccomplex_t IT_0440 = 4*IT_0266;
    const ccomplex_t IT_0441 = 4*conj(IT_0266);
    const ccomplex_t IT_0442 = 4*conj(IT_0253);
    const ccomplex_t IT_0443 = 0.125*IT_0442;
    const ccomplex_t IT_0444 = conj(IT_0383) + IT_0443;
    const ccomplex_t IT_0445 = 4*IT_0253;
    const ccomplex_t IT_0446 = 0.125*IT_0445;
    const ccomplex_t IT_0447 = IT_0383 + IT_0446;
    const ccomplex_t IT_0448 = 0.125*conj(IT_0266)*IT_0304 + 0.125*IT_0266
      *IT_0307 + conj(IT_0388)*IT_0418 + IT_0388*IT_0421 + 0.125*conj(IT_0055)
      *IT_0437 + 0.125*IT_0055*IT_0439 + 0.125*IT_0351*IT_0440 + 0.125*IT_0355
      *IT_0441 + 0.125*IT_0429*IT_0444 + 0.125*IT_0433*IT_0447;
    const ccomplex_t IT_0449 = 16*IT_0230;
    const ccomplex_t IT_0450 = 0.125*IT_0449;
    const ccomplex_t IT_0451 = IT_0145 + IT_0450;
    const ccomplex_t IT_0452 = 16*IT_0224;
    const ccomplex_t IT_0453 = 0.125*IT_0452;
    const ccomplex_t IT_0454 = IT_0451 + IT_0453;
    const ccomplex_t IT_0455 = 4*IT_0388;
    const ccomplex_t IT_0456 = 16*conj(IT_0230);
    const ccomplex_t IT_0457 = 16*IT_0219;
    const ccomplex_t IT_0458 = 8*IT_0055;
    const ccomplex_t IT_0459 = 0.125*IT_0422;
    const ccomplex_t IT_0460 = IT_0347 + IT_0459;
    const ccomplex_t IT_0461 = 8*conj(IT_0055);
    const ccomplex_t IT_0462 = conj(IT_0145)*IT_0418 + 0.125*conj(IT_0383)
      *IT_0436 + conj(IT_0230)*IT_0451 + conj(IT_0224)*IT_0454 + 0.125*conj
      (IT_0388)*IT_0455 + 0.125*IT_0224*IT_0456 + 0.125*IT_0343*(16*IT_0235 +
       IT_0457 + IT_0458) + 0.125*IT_0460*IT_0461;
    const ccomplex_t IT_0463 = IT_0005*IT_0006;
    const ccomplex_t IT_0464 = 16*IT_0214;
    const ccomplex_t IT_0465 = 0.125*IT_0464;
    const ccomplex_t IT_0466 = IT_0189 + IT_0465;
    const ccomplex_t IT_0467 = 16*IT_0234;
    const ccomplex_t IT_0468 = 0.125*IT_0467;
    const ccomplex_t IT_0469 = IT_0466 + IT_0468;
    const ccomplex_t IT_0470 = IT_0398 + IT_0445;
    const ccomplex_t IT_0471 = IT_0399 + IT_0442;
    const ccomplex_t IT_0472 = 16*conj(IT_0214);
    const ccomplex_t IT_0473 = 16*conj(IT_0234);
    const ccomplex_t IT_0474 = 16*conj(IT_0224);
    const ccomplex_t IT_0475 = 0.125*IT_0425;
    const ccomplex_t IT_0476 = IT_0367 + IT_0475;
    const ccomplex_t IT_0477 = 8*IT_0180;
    const ccomplex_t IT_0478 = 0.125*IT_0424;
    const ccomplex_t IT_0479 = IT_0371 + IT_0478;
    const ccomplex_t IT_0480 = 8*conj(IT_0180);
    const ccomplex_t IT_0481 = 0.125*conj(IT_0388)*IT_0437 + 0.125*IT_0388
      *IT_0439 + 0.125*IT_0433*IT_0454 + 0.125*IT_0421*IT_0458 + 0.125*IT_0418
      *IT_0461 + IT_0359*IT_0469 + 0.125*conj(IT_0266)*IT_0470 + 0.125*IT_0266
      *IT_0471 + IT_0363*(conj(IT_0189) + 0.125*IT_0472 + 0.125*IT_0473) + 0.125
      *IT_0429*(conj(IT_0145) + 0.125*IT_0456 + 0.125*IT_0474) + 0.125*IT_0476
      *IT_0477 + 0.125*IT_0479*IT_0480;
    const ccomplex_t IT_0482 = 16*IT_0231;
    const ccomplex_t IT_0483 = 0.125*conj(IT_0266)*IT_0440 + 0.125*conj
      (IT_0253)*IT_0445 + conj(IT_0214)*IT_0466 + conj(IT_0234)*IT_0469 + 0.125
      *IT_0234*IT_0472 + conj(IT_0189)*IT_0479 + 0.125*IT_0432*IT_0480 + 0.125
      *IT_0359*(16*IT_0226 + IT_0477 + IT_0482);
    const ccomplex_t IT_0484 = s_34 + IT_0006;
    const ccomplex_t IT_0485 = s_12*IT_0484;
    const ccomplex_t IT_0486 = (-8)*IT_0329;
    const ccomplex_t IT_0487 = (-8)*conj(IT_0329);
    const ccomplex_t IT_0488 = -conj(IT_0336);
    const ccomplex_t IT_0489 = conj(IT_0303) + IT_0488;
    const ccomplex_t IT_0490 = -IT_0336;
    const ccomplex_t IT_0491 = IT_0005*IT_0484;
    const ccomplex_t IT_0492 = conj(IT_0329)*(IT_0393 + IT_0400) + IT_0329*
      (IT_0394 + IT_0401) + 8*conj(IT_0336)*(IT_0329 + 0.125*IT_0391 + 0.125
      *IT_0402) + 8*IT_0336*(conj(IT_0329) + 0.125*IT_0392 + 0.125*IT_0403) +
       conj(IT_0303)*(IT_0398 + IT_0405 + IT_0486) + IT_0303*(IT_0399 + IT_0406 
      + IT_0487) + 8*IT_0322*(0.125*IT_0397 + 0.125*IT_0408 + IT_0489) + 8*conj
      (IT_0322)*(IT_0303 + 0.125*IT_0396 + 0.125*IT_0407 + IT_0490);
    const ccomplex_t IT_0493 = IT_0005*IT_0411;
    const ccomplex_t IT_0494 = 0.125*IT_0266;
    const ccomplex_t IT_0495 = 0.125*conj(IT_0266);
    const ccomplex_t IT_0496 = 0.125*conj(IT_0253)*IT_0269 + 0.125*IT_0253
      *IT_0272 + (-0.25)*IT_0168*IT_0383 + 0.25*IT_0164*conj(IT_0383) + conj
      (IT_0090)*(IT_0163 + IT_0165 + 0.25*IT_0388) + conj(IT_0089)*(IT_0164 + (
      -0.25)*IT_0388) + -(IT_0089 + -IT_0090)*(IT_0168 + 1./4*conj(IT_0388)) + (
      -0.25)*IT_0163*(conj(IT_0383) + (-0.5)*IT_0413) + IT_0259*IT_0494 +
       IT_0256*IT_0495;
    const ccomplex_t IT_0497 = s_23*m_N_2*m_tau;
    const ccomplex_t IT_0498 = 0.125*IT_0423;
    const ccomplex_t IT_0499 = IT_0343 + IT_0498;
    const ccomplex_t IT_0500 = 0.125*IT_0347;
    const ccomplex_t IT_0501 = 0.125*IT_0343;
    const ccomplex_t IT_0502 = 0.125*conj(IT_0145)*IT_0437 + 0.125*IT_0145
      *IT_0439 + IT_0355*IT_0444 + IT_0351*IT_0447 + conj(IT_0388)*IT_0460 +
       IT_0326*IT_0494 + IT_0323*IT_0495 + IT_0388*IT_0499 + IT_0441*IT_0500 +
       IT_0440*IT_0501;
    const ccomplex_t IT_0503 = 0.125*IT_0363;
    const ccomplex_t IT_0504 = 0.125*IT_0359;
    const ccomplex_t IT_0505 = 0.125*IT_0180;
    const ccomplex_t IT_0506 = 0.125*conj(IT_0180);
    const ccomplex_t IT_0507 = 0.125*conj(IT_0055)*IT_0330 + 0.125*IT_0055
      *IT_0333 + IT_0355*IT_0428 + IT_0351*IT_0432 + IT_0367*IT_0460 + IT_0371
      *IT_0499 + IT_0425*IT_0500 + IT_0424*IT_0501 + IT_0419*IT_0503 + IT_0416
      *IT_0504 + IT_0307*IT_0505 + IT_0304*IT_0506;
    const ccomplex_t IT_0508 = s_13*m_N_2*m_tau;
    const ccomplex_t IT_0509 = 0.125*conj(IT_0388);
    const ccomplex_t IT_0510 = 0.125*IT_0371;
    const ccomplex_t IT_0511 = 4*conj(IT_0388);
    const ccomplex_t IT_0512 = 0.125*IT_0436;
    const ccomplex_t IT_0513 = conj(IT_0253)*IT_0371 + 0.125*IT_0340*IT_0388 +
       IT_0266*IT_0428 + conj(IT_0266)*IT_0432 + 0.125*conj(IT_0189)*IT_0470 +
       0.125*IT_0189*IT_0471 + IT_0455*IT_0504 + IT_0337*IT_0509 + IT_0438
      *IT_0510 + IT_0503*IT_0511 + IT_0367*(IT_0253 + IT_0512);
    const ccomplex_t IT_0514 = s_24*m_N_2*m_tau;
    const ccomplex_t IT_0515 = 0.125*IT_0333*IT_0388 + IT_0363*(conj(IT_0253) 
      + 0.125*IT_0438) + 0.125*IT_0367*IT_0455 + IT_0266*IT_0476 + conj(IT_0266)
      *IT_0479 + IT_0471*IT_0505 + IT_0470*IT_0506 + IT_0330*IT_0509 + IT_0510
      *IT_0511 + IT_0359*(IT_0253 + IT_0512);
    const ccomplex_t IT_0516 = IT_0004*(conj(IT_0055)*IT_0092 + IT_0055
      *IT_0094 + conj(IT_0145)*IT_0166 + IT_0145*IT_0168) + (IT_0094*IT_0145 +
       IT_0092*conj(IT_0145) + conj(IT_0055)*IT_0166 + IT_0055*IT_0168)*IT_0172 
      + IT_0173*(conj(IT_0180)*IT_0182 + IT_0180*IT_0184 + conj(IT_0189)*IT_0191
       + IT_0189*IT_0193) + (IT_0184*IT_0189 + IT_0182*conj(IT_0189) + conj
      (IT_0180)*IT_0191 + IT_0180*IT_0193)*IT_0194 + IT_0195*(IT_0168*IT_0220 +
       IT_0166*IT_0221 + IT_0094*IT_0227 + IT_0092*IT_0228) + IT_0229*(IT_0193
      *IT_0232 + IT_0191*IT_0233 + IT_0184*IT_0236 + IT_0182*IT_0237) + IT_0246*
      (conj(IT_0253)*IT_0256 + IT_0253*IT_0259 + conj(IT_0266)*IT_0269 + IT_0266
      *IT_0272) + IT_0274*IT_0375 + IT_0377*IT_0378 + (IT_0259*IT_0383 + IT_0256
      *conj(IT_0383) + IT_0272*IT_0388 + IT_0269*conj(IT_0388))*(IT_0242 +
       IT_0243 + IT_0389) + IT_0390*(conj(IT_0383)*IT_0391 + IT_0383*IT_0392 +
       conj(IT_0388)*IT_0393 + IT_0388*IT_0394) + IT_0395*(conj(IT_0303)*IT_0396
       + IT_0303*IT_0397 + conj(IT_0322)*IT_0398 + IT_0322*IT_0399 + conj
      (IT_0336)*IT_0400 + IT_0336*IT_0401 + conj(IT_0329)*IT_0402 + IT_0329
      *IT_0403) + IT_0404*(conj(IT_0322)*IT_0391 + IT_0322*IT_0392 + conj
      (IT_0303)*IT_0393 + IT_0303*IT_0394 + conj(IT_0329)*IT_0405 + IT_0329
      *IT_0406 + conj(IT_0336)*IT_0407 + IT_0336*IT_0408) + (IT_0184*IT_0232 +
       IT_0182*IT_0233 + (conj(IT_0164) + IT_0192)*IT_0236 + IT_0191*IT_0237)
      *IT_0409 + 8*IT_0412*(IT_0089*(conj(IT_0089) + -conj(IT_0090)) + (IT_0163 
      + IT_0165)*IT_0168 + IT_0090*(conj(IT_0090) + 0.125*IT_0413)) + (IT_0094
      *IT_0220 + IT_0092*IT_0221 + (conj(IT_0163) + IT_0167)*IT_0227 + IT_0166
      *IT_0228)*IT_0414 + 8*IT_0415*IT_0434 + 8*IT_0435*IT_0448 + 8*IT_0239
      *IT_0462 + 8*IT_0463*IT_0481 + 8*IT_0240*IT_0483 + 8*IT_0485*(IT_0329*conj
      (IT_0329) + 0.125*conj(IT_0322)*IT_0486 + IT_0322*(conj(IT_0322) + 0.125
      *IT_0487) + IT_0489*(IT_0303 + IT_0490)) + IT_0491*IT_0492 + 8*IT_0493
      *IT_0496 + 8*IT_0497*IT_0502 + 8*IT_0238*IT_0507 + 8*IT_0508*IT_0513 + 8
      *IT_0514*IT_0515;
    return create_ccomplex_return(IT_0516);
}

