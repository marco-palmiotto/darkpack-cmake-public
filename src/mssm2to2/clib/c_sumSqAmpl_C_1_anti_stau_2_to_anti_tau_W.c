#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_stau_2_to_anti_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_stau_2_to_anti_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
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
    const creal_t m_stau_2 = param->m_stau_2;
    const creal_t reg_prop = param->reg_prop;
    const creal_t m_snu_tau = param->m_snu_tau;
    const creal_t Gamma_nutaul = param->Gamma_nutaul;
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
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_stau_2, 2);
    const ccomplex_t IT_0002 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 0.5*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = conj(N_B1)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0004*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = conj(N_d1)*e_em*m_tau*IT_0020*conj(U_stau_01);
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + 0.5*IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0030 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0031 = IT_0004*IT_0030;
    const ccomplex_t IT_0032 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0033 = IT_0004*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + 0.5*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = conj(N_B2)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0038 = IT_0013*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_d2)*e_em*m_tau*IT_0020*conj(U_stau_01);
    const ccomplex_t IT_0041 = IT_0019*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 0.5*IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = IT_0029*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0049 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0050 = IT_0004*IT_0049;
    const ccomplex_t IT_0051 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + 0.5*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = conj(N_B3)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0057 = IT_0013*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = conj(N_d3)*e_em*m_tau*IT_0020*conj(U_stau_01);
    const ccomplex_t IT_0060 = IT_0019*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 0.5*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = IT_0055*IT_0063;
    const ccomplex_t IT_0065 = IT_0048*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0068 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0071 = IT_0004*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + 0.5*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = conj(N_B4)*e_em*conj(U_stau_11);
    const ccomplex_t IT_0076 = IT_0013*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = conj(N_d4)*e_em*m_tau*IT_0020*conj(U_stau_01);
    const ccomplex_t IT_0079 = IT_0019*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 0.5*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = IT_0067*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = IT_0028 + IT_0047 + IT_0066 + IT_0085;
    const ccomplex_t IT_0087 = m_C_1*IT_0086;
    const ccomplex_t IT_0088 = (-0.5)*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0090 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0091 = IT_0004*IT_0090;
    const ccomplex_t IT_0092 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0093 = IT_0004*IT_0092;
    const ccomplex_t IT_0094 = 1.4142135623731*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*(IT_0091 + (-0.5)*IT_0094);
    const ccomplex_t IT_0096 = -IT_0095;
    const ccomplex_t IT_0097 = IT_0025*IT_0096;
    const ccomplex_t IT_0098 = IT_0002*IT_0089*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0100 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0101 = IT_0004*IT_0100;
    const ccomplex_t IT_0102 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0103 = IT_0004*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + (-0.5)*IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0044*IT_0106;
    const ccomplex_t IT_0108 = IT_0029*IT_0099*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0110 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0111 = IT_0004*IT_0110;
    const ccomplex_t IT_0112 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0113 = IT_0004*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + (-0.5)*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0063*IT_0116;
    const ccomplex_t IT_0118 = IT_0048*IT_0109*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0120 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0121 = IT_0004*IT_0120;
    const ccomplex_t IT_0122 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0123 = IT_0004*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + (-0.5)*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0082*IT_0126;
    const ccomplex_t IT_0128 = IT_0067*IT_0119*IT_0127;
    const ccomplex_t IT_0129 = 0.5*IT_0098 + 0.5*IT_0108 + 0.5*IT_0118 + 0.5
      *IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004
      *conj(U_stau_01);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1.4142135623731)*conj(U_d1)
      *e_em*m_tau*IT_0004*IT_0018*IT_0020;
    const ccomplex_t IT_0133 = (-0.5)*IT_0132;
    const ccomplex_t IT_0134 = pow(m_C_1, 2);
    const ccomplex_t IT_0135 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0134 + (
      -0.5)*reg_prop + 0.5*m_snu_tau*(m_snu_tau + (0 + _Complex_I*-1)
      *Gamma_nutaul), -1);
    const ccomplex_t IT_0136 = IT_0133*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*IT_0136;
    const ccomplex_t IT_0138 = IT_0131*IT_0137;
    const ccomplex_t IT_0139 = -IT_0028 + -IT_0047 + -IT_0066 + -IT_0085 + (
      -0.5)*IT_0138;
    const ccomplex_t IT_0140 = pow(m_W, -2);
    const ccomplex_t IT_0141 = s_14*s_34*IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = s_13 + IT_0142;
    const ccomplex_t IT_0144 = m_C_1*IT_0143;
    const ccomplex_t IT_0145 = 2*IT_0144;
    const ccomplex_t IT_0146 = IT_0139*IT_0145;
    const ccomplex_t IT_0147 = pow(m_W, 2);
    const ccomplex_t IT_0148 = s_13*IT_0147;
    const ccomplex_t IT_0149 = s_14*s_34;
    const ccomplex_t IT_0150 = (-2)*IT_0149;
    const ccomplex_t IT_0151 = IT_0148 + IT_0150;
    const ccomplex_t IT_0152 = IT_0140*IT_0151;
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = s_13 + IT_0153;
    const ccomplex_t IT_0155 = 4*IT_0154;
    const ccomplex_t IT_0156 = IT_0088*IT_0155;
    const ccomplex_t IT_0157 = 0.5*IT_0138;
    const ccomplex_t IT_0158 = pow(s_34, 2);
    const ccomplex_t IT_0159 = IT_0140*IT_0158;
    const ccomplex_t IT_0160 = -IT_0000;
    const ccomplex_t IT_0161 = IT_0159 + IT_0160;
    const ccomplex_t IT_0162 = m_C_1*IT_0161;
    const ccomplex_t IT_0163 = (-2)*IT_0162;
    const ccomplex_t IT_0164 = IT_0157*IT_0163;
    const ccomplex_t IT_0165 = 0.5*IT_0028 + 0.5*IT_0047 + 0.5*IT_0066 + 0.5
      *IT_0085;
    const ccomplex_t IT_0166 = s_34*IT_0147;
    const ccomplex_t IT_0167 = IT_0140*IT_0166;
    const ccomplex_t IT_0168 = (-0.25)*IT_0167;
    const ccomplex_t IT_0169 = s_34 + IT_0168;
    const ccomplex_t IT_0170 = m_C_1*IT_0169;
    const ccomplex_t IT_0171 = (-4)*IT_0170;
    const ccomplex_t IT_0172 = (-2)*IT_0171;
    const ccomplex_t IT_0173 = IT_0165*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*e_em*V_Wp1*IT_0004;
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = IT_0135*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = IT_0131*IT_0177;
    const ccomplex_t IT_0179 = 0.5*IT_0178;
    const ccomplex_t IT_0180 = pow(s_14, 2);
    const ccomplex_t IT_0181 = IT_0140*IT_0180;
    const ccomplex_t IT_0182 = -IT_0134;
    const ccomplex_t IT_0183 = IT_0181 + IT_0182;
    const ccomplex_t IT_0184 = m_tau*IT_0183;
    const ccomplex_t IT_0185 = (-2)*IT_0184;
    const ccomplex_t IT_0186 = IT_0179*IT_0185;
    const ccomplex_t IT_0187 = N_B1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0188 = IT_0013*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = N_W1*e_em*conj(U_stau_01);
    const ccomplex_t IT_0191 = IT_0004*IT_0190;
    const ccomplex_t IT_0192 = 1.4142135623731*IT_0191;
    const ccomplex_t IT_0193 = N_d1*e_em*m_tau*IT_0020*conj(U_stau_11);
    const ccomplex_t IT_0194 = IT_0019*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0189 + IT_0192 + 
      -IT_0195);
    const ccomplex_t IT_0197 = (-0.5)*IT_0196;
    const ccomplex_t IT_0198 = IT_0011*IT_0197;
    const ccomplex_t IT_0199 = IT_0002*IT_0089*IT_0198;
    const ccomplex_t IT_0200 = N_B2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0201 = IT_0013*IT_0200;
    const ccomplex_t IT_0202 = 1.4142135623731*IT_0201;
    const ccomplex_t IT_0203 = N_W2*e_em*conj(U_stau_01);
    const ccomplex_t IT_0204 = IT_0004*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = N_d2*e_em*m_tau*IT_0020*conj(U_stau_11);
    const ccomplex_t IT_0207 = IT_0019*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0202 + IT_0205 + 
      -IT_0208);
    const ccomplex_t IT_0210 = (-0.5)*IT_0209;
    const ccomplex_t IT_0211 = IT_0036*IT_0210;
    const ccomplex_t IT_0212 = IT_0029*IT_0099*IT_0211;
    const ccomplex_t IT_0213 = N_B3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0214 = IT_0013*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = N_W3*e_em*conj(U_stau_01);
    const ccomplex_t IT_0217 = IT_0004*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = N_d3*e_em*m_tau*IT_0020*conj(U_stau_11);
    const ccomplex_t IT_0220 = IT_0019*IT_0219;
    const ccomplex_t IT_0221 = 1.4142135623731*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*(IT_0215 + IT_0218 + 
      -IT_0221);
    const ccomplex_t IT_0223 = (-0.5)*IT_0222;
    const ccomplex_t IT_0224 = IT_0055*IT_0223;
    const ccomplex_t IT_0225 = IT_0048*IT_0109*IT_0224;
    const ccomplex_t IT_0226 = N_B4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0227 = IT_0013*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = N_W4*e_em*conj(U_stau_01);
    const ccomplex_t IT_0230 = IT_0004*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = N_d4*e_em*m_tau*IT_0020*conj(U_stau_11);
    const ccomplex_t IT_0233 = IT_0019*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*(IT_0228 + IT_0231 + 
      -IT_0234);
    const ccomplex_t IT_0236 = (-0.5)*IT_0235;
    const ccomplex_t IT_0237 = IT_0074*IT_0236;
    const ccomplex_t IT_0238 = IT_0067*IT_0119*IT_0237;
    const ccomplex_t IT_0239 = (-0.5)*IT_0199 + (-0.5)*IT_0212 + (-0.5)
      *IT_0225 + (-0.5)*IT_0238;
    const ccomplex_t IT_0240 = IT_0096*IT_0197;
    const ccomplex_t IT_0241 = IT_0002*IT_0240;
    const ccomplex_t IT_0242 = (0 + _Complex_I*1)*IT_0241;
    const ccomplex_t IT_0243 = IT_0106*IT_0210;
    const ccomplex_t IT_0244 = IT_0029*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = IT_0116*IT_0223;
    const ccomplex_t IT_0247 = IT_0048*IT_0246;
    const ccomplex_t IT_0248 = (0 + _Complex_I*1)*IT_0247;
    const ccomplex_t IT_0249 = IT_0126*IT_0236;
    const ccomplex_t IT_0250 = IT_0067*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*IT_0250;
    const ccomplex_t IT_0252 = IT_0242 + IT_0245 + IT_0248 + IT_0251;
    const ccomplex_t IT_0253 = m_C_1*IT_0252;
    const ccomplex_t IT_0254 = 0.5*IT_0253;
    const ccomplex_t IT_0255 = IT_0239 + IT_0254;
    const ccomplex_t IT_0256 = m_C_1*m_tau;
    const ccomplex_t IT_0257 = (-6)*IT_0256;
    const ccomplex_t IT_0258 = IT_0255*IT_0257;
    const ccomplex_t IT_0259 = IT_0129*IT_0155;
    const ccomplex_t IT_0260 = 0.5*IT_0242 + 0.5*IT_0245 + 0.5*IT_0248 + 0.5
      *IT_0251;
    const ccomplex_t IT_0261 = s_14*IT_0147;
    const ccomplex_t IT_0262 = IT_0140*IT_0261;
    const ccomplex_t IT_0263 = 0.5*IT_0262;
    const ccomplex_t IT_0264 = s_14 + IT_0263;
    const ccomplex_t IT_0265 = m_tau*IT_0264;
    const ccomplex_t IT_0266 = 2*IT_0265;
    const ccomplex_t IT_0267 = 2*IT_0266;
    const ccomplex_t IT_0268 = IT_0260*IT_0267;
    const ccomplex_t IT_0269 = (-0.5)*IT_0178;
    const ccomplex_t IT_0270 = m_tau*IT_0143;
    const ccomplex_t IT_0271 = 2*IT_0270;
    const ccomplex_t IT_0272 = IT_0269*IT_0271;
    const ccomplex_t IT_0273 = -IT_0242 + -IT_0245 + -IT_0248 + -IT_0251;
    const ccomplex_t IT_0274 = 2*IT_0184;
    const ccomplex_t IT_0275 = IT_0273*IT_0274;
    const ccomplex_t IT_0276 = IT_0146 + IT_0156 + IT_0164 + IT_0173 + IT_0186
       + IT_0258 + IT_0259 + IT_0268 + IT_0272 + IT_0275;
    const ccomplex_t IT_0277 = IT_0145*conj(IT_0179);
    const ccomplex_t IT_0278 = (-2)*IT_0144;
    const ccomplex_t IT_0279 = conj(IT_0273)*IT_0278;
    const ccomplex_t IT_0280 = IT_0163*conj(IT_0269);
    const ccomplex_t IT_0281 = s_13*s_34;
    const ccomplex_t IT_0282 = s_14*IT_0000;
    const ccomplex_t IT_0283 = -IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = 2*IT_0284;
    const ccomplex_t IT_0286 = s_13*IT_0143;
    const ccomplex_t IT_0287 = 2*IT_0286;
    const ccomplex_t IT_0288 = s_13*IT_0161;
    const ccomplex_t IT_0289 = 2*IT_0288;
    const ccomplex_t IT_0290 = IT_0161*IT_0256;
    const ccomplex_t IT_0291 = 2*IT_0290;
    const ccomplex_t IT_0292 = (-2)*IT_0286;
    const ccomplex_t IT_0293 = IT_0143*IT_0256;
    const ccomplex_t IT_0294 = (-2)*IT_0293;
    const ccomplex_t IT_0295 = IT_0088 + IT_0129;
    const ccomplex_t IT_0296 = conj(IT_0139)*IT_0145;
    const ccomplex_t IT_0297 = conj(IT_0273)*IT_0274;
    const ccomplex_t IT_0298 = conj(IT_0269)*IT_0271;
    const ccomplex_t IT_0299 = conj(IT_0239)*IT_0257;
    const ccomplex_t IT_0300 = conj(IT_0157)*IT_0163;
    const ccomplex_t IT_0301 = conj(IT_0179)*IT_0185;
    const ccomplex_t IT_0302 = s_13*IT_0183;
    const ccomplex_t IT_0303 = s_13*s_14;
    const ccomplex_t IT_0304 = s_34*IT_0134;
    const ccomplex_t IT_0305 = -IT_0304;
    const ccomplex_t IT_0306 = IT_0303 + IT_0305;
    const ccomplex_t IT_0307 = 2*IT_0306;
    const ccomplex_t IT_0308 = 0.5*IT_0269;
    const ccomplex_t IT_0309 = 0.5*IT_0157;
    const ccomplex_t IT_0310 = IT_0183*IT_0256;
    const ccomplex_t IT_0311 = 2*IT_0310;
    const ccomplex_t IT_0312 = 0.5*IT_0179;
    const ccomplex_t IT_0313 = (-2)*IT_0310;
    const ccomplex_t IT_0314 = 0.5*IT_0273;
    const ccomplex_t IT_0315 = (-2)*IT_0302;
    const ccomplex_t IT_0316 = 0.5*IT_0139;
    const ccomplex_t IT_0317 = 2*IT_0293;
    const ccomplex_t IT_0318 = (-0.5)*IT_0273;
    const ccomplex_t IT_0319 = (-0.5)*IT_0139;
    const ccomplex_t IT_0320 = (-0.5)*IT_0255;
    const ccomplex_t IT_0321 = pow(m_W, 4);
    const ccomplex_t IT_0322 = s_13*IT_0321;
    const ccomplex_t IT_0323 = IT_0140*IT_0322;
    const ccomplex_t IT_0324 = (-0.5)*IT_0323;
    const ccomplex_t IT_0325 = IT_0148 + IT_0150 + IT_0324;
    const ccomplex_t IT_0326 = (-4)*IT_0325;
    const ccomplex_t IT_0327 = IT_0147*IT_0256;
    const ccomplex_t IT_0328 = (-3)*IT_0327;
    const ccomplex_t IT_0329 = (-2)*IT_0328;
    const ccomplex_t IT_0330 = 2*IT_0171;
    const ccomplex_t IT_0331 = (-2)*IT_0266;
    const ccomplex_t IT_0332 = IT_0260*IT_0330;
    const ccomplex_t IT_0333 = IT_0157*IT_0271;
    const ccomplex_t IT_0334 = IT_0145*IT_0179;
    const ccomplex_t IT_0335 = IT_0155*(IT_0239 + IT_0254) + IT_0332 + IT_0333
       + IT_0334;
    const ccomplex_t IT_0336 = IT_0163*IT_0269;
    const ccomplex_t IT_0337 = IT_0273*IT_0278;
    const ccomplex_t IT_0338 = IT_0139*IT_0185;
    const ccomplex_t IT_0339 = IT_0165*IT_0331;
    const ccomplex_t IT_0340 = (conj(IT_0088) + conj(IT_0129))*IT_0276 +
       IT_0255*(IT_0277 + IT_0279 + IT_0280) + conj(IT_0269)*(IT_0260*IT_0285 +
       IT_0273*IT_0287 + IT_0269*IT_0289 + IT_0157*IT_0291 + IT_0179*IT_0292 +
       IT_0139*IT_0294) + IT_0295*(IT_0296 + IT_0297 + IT_0298 + IT_0299 +
       IT_0300 + IT_0301) + 2*conj(IT_0139)*(0.5*IT_0185*IT_0255 + IT_0139
      *IT_0302 + 0.5*IT_0165*IT_0307 + IT_0294*IT_0308 + IT_0292*IT_0309 +
       IT_0311*IT_0312 + IT_0313*IT_0314) + 2*conj(IT_0179)*(IT_0179*IT_0302 + 
      -IT_0260*IT_0306 + IT_0292*IT_0308 + IT_0294*IT_0309 + IT_0314*IT_0315 +
       IT_0311*IT_0316) + 2*conj(IT_0273)*(IT_0273*IT_0302 + 0.5*IT_0260*IT_0307
       + IT_0287*IT_0308 + IT_0312*IT_0315 + IT_0313*IT_0316 + IT_0309*IT_0317) 
      + (-2)*conj(IT_0157)*(IT_0165*IT_0284 + (-0.5)*IT_0157*IT_0289 + (-0.5)
      *IT_0269*IT_0291 + (-0.5)*IT_0179*IT_0294 + IT_0317*IT_0318 + IT_0292
      *IT_0319 + IT_0271*IT_0320) + (-2)*conj(IT_0260)*((-0.5)*IT_0269*IT_0285 +
       (-0.5)*IT_0267*IT_0295 + IT_0179*IT_0306 + IT_0307*IT_0318 + (-0.5)
      *IT_0260*IT_0326 + (-0.5)*IT_0165*IT_0329 + IT_0320*IT_0330) + (-2)*conj
      (IT_0165)*(IT_0157*IT_0284 + (-0.5)*IT_0172*IT_0295 + IT_0307*IT_0319 + (
      -0.5)*IT_0165*IT_0326 + (-0.5)*IT_0260*IT_0329 + IT_0320*IT_0331) + conj
      (IT_0239)*(IT_0335 + IT_0336 + IT_0337 + IT_0338 + IT_0339) + conj(IT_0254
      )*(IT_0257*IT_0295 + IT_0335 + IT_0336 + IT_0337 + IT_0338 + IT_0339);
    return create_ccomplex_return(IT_0340);
}

