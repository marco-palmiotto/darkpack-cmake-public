#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_stau_1_to_tau_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_stau_1_to_tau_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
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
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_stau_1, 2);
    const ccomplex_t IT_0002 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + (-0.5)*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = cos(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0015 = IT_0013*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = cos(beta);
    const ccomplex_t IT_0018 = cpow(IT_0017, -1);
    const ccomplex_t IT_0019 = IT_0004*IT_0018;
    const ccomplex_t IT_0020 = pow(m_W, -1);
    const ccomplex_t IT_0021 = N_d1*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0022 = IT_0019*IT_0021;
    const ccomplex_t IT_0023 = 1.4142135623731*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*(IT_0016 + 0.5*IT_0023);
    const ccomplex_t IT_0025 = -IT_0024;
    const ccomplex_t IT_0026 = IT_0011*IT_0025;
    const ccomplex_t IT_0027 = IT_0002*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*IT_0027;
    const ccomplex_t IT_0029 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0030 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0031 = IT_0004*IT_0030;
    const ccomplex_t IT_0032 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0033 = IT_0004*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + (-0.5)*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0038 = IT_0013*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = N_d4*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0041 = IT_0019*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 0.5*IT_0042);
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0036*IT_0044;
    const ccomplex_t IT_0046 = IT_0029*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*IT_0046;
    const ccomplex_t IT_0048 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0049 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0050 = IT_0004*IT_0049;
    const ccomplex_t IT_0051 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0052 = IT_0004*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = (0 + _Complex_I*1)*(IT_0050 + (-0.5)*IT_0053);
    const ccomplex_t IT_0055 = -IT_0054;
    const ccomplex_t IT_0056 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0057 = IT_0013*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = N_d2*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0060 = IT_0019*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + 0.5*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = IT_0055*IT_0063;
    const ccomplex_t IT_0065 = IT_0048*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = cpow(s_23 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0068 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0069 = IT_0004*IT_0068;
    const ccomplex_t IT_0070 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0071 = IT_0004*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*(IT_0069 + (-0.5)*IT_0072);
    const ccomplex_t IT_0074 = -IT_0073;
    const ccomplex_t IT_0075 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0076 = IT_0013*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = N_d3*e_em*m_tau*IT_0020*U_stau_00;
    const ccomplex_t IT_0079 = IT_0019*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*(IT_0077 + 0.5*IT_0080);
    const ccomplex_t IT_0082 = -IT_0081;
    const ccomplex_t IT_0083 = IT_0074*IT_0082;
    const ccomplex_t IT_0084 = IT_0067*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = 0.5*IT_0028 + 0.5*IT_0047 + 0.5*IT_0066 + 0.5
      *IT_0085;
    const ccomplex_t IT_0087 = pow(m_W, -2);
    const ccomplex_t IT_0088 = pow(m_W, 2);
    const ccomplex_t IT_0089 = s_14*IT_0088;
    const ccomplex_t IT_0090 = IT_0087*IT_0089;
    const ccomplex_t IT_0091 = 0.5*IT_0090;
    const ccomplex_t IT_0092 = s_14 + IT_0091;
    const ccomplex_t IT_0093 = m_tau*IT_0092;
    const ccomplex_t IT_0094 = 2*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0096 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0097 = IT_0013*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0100 = IT_0004*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = conj(N_d1)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0103 = IT_0019*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0098 + IT_0101 + 
      -IT_0104);
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0011*IT_0106;
    const ccomplex_t IT_0108 = IT_0002*IT_0095*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0110 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0111 = IT_0013*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0114 = IT_0004*IT_0113;
    const ccomplex_t IT_0115 = 1.4142135623731*IT_0114;
    const ccomplex_t IT_0116 = conj(N_d2)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0117 = IT_0019*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + IT_0115 + 
      -IT_0118);
    const ccomplex_t IT_0120 = (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = IT_0055*IT_0120;
    const ccomplex_t IT_0122 = IT_0048*IT_0109*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0124 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0125 = IT_0013*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0128 = IT_0004*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = conj(N_d3)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0131 = IT_0019*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0126 + IT_0129 + 
      -IT_0132);
    const ccomplex_t IT_0134 = (-0.5)*IT_0133;
    const ccomplex_t IT_0135 = IT_0074*IT_0134;
    const ccomplex_t IT_0136 = IT_0067*IT_0123*IT_0135;
    const ccomplex_t IT_0137 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0138 = IT_0004*IT_0137;
    const ccomplex_t IT_0139 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0140 = IT_0004*IT_0139;
    const ccomplex_t IT_0141 = 1.4142135623731*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*(IT_0138 + 0.5*IT_0141);
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = IT_0134*IT_0143;
    const ccomplex_t IT_0145 = IT_0067*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*IT_0145;
    const ccomplex_t IT_0147 = m_C_1*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0149 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0150 = IT_0013*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0153 = IT_0004*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = conj(N_d4)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0156 = IT_0019*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0151 + IT_0154 + 
      -IT_0157);
    const ccomplex_t IT_0159 = (-0.5)*IT_0158;
    const ccomplex_t IT_0160 = IT_0036*IT_0159;
    const ccomplex_t IT_0161 = IT_0029*IT_0148*IT_0160;
    const ccomplex_t IT_0162 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0163 = IT_0004*IT_0162;
    const ccomplex_t IT_0164 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0165 = IT_0004*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + 0.5*IT_0166);
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = IT_0106*IT_0168;
    const ccomplex_t IT_0170 = IT_0002*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_C_1*IT_0171;
    const ccomplex_t IT_0173 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0174 = IT_0004*IT_0173;
    const ccomplex_t IT_0175 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0176 = IT_0004*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*(IT_0174 + 0.5*IT_0177);
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0120*IT_0179;
    const ccomplex_t IT_0181 = IT_0048*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = m_C_1*IT_0182;
    const ccomplex_t IT_0184 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0185 = IT_0004*IT_0184;
    const ccomplex_t IT_0186 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0187 = IT_0004*IT_0186;
    const ccomplex_t IT_0188 = 1.4142135623731*IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*(IT_0185 + 0.5*IT_0188);
    const ccomplex_t IT_0190 = -IT_0189;
    const ccomplex_t IT_0191 = IT_0159*IT_0190;
    const ccomplex_t IT_0192 = IT_0029*IT_0191;
    const ccomplex_t IT_0193 = (0 + _Complex_I*1)*IT_0192;
    const ccomplex_t IT_0194 = m_C_1*IT_0193;
    const ccomplex_t IT_0195 = 0.5*IT_0108 + 0.5*IT_0122 + 0.5*IT_0136 + (-0.5
      )*IT_0147 + 0.5*IT_0161 + (-0.5)*IT_0172 + (-0.5)*IT_0183 + (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = s_34*IT_0088;
    const ccomplex_t IT_0197 = IT_0087*IT_0196;
    const ccomplex_t IT_0198 = (-0.25)*IT_0197;
    const ccomplex_t IT_0199 = s_34 + IT_0198;
    const ccomplex_t IT_0200 = m_C_1*IT_0199;
    const ccomplex_t IT_0201 = (-4)*IT_0200;
    const ccomplex_t IT_0202 = m_C_1*IT_0085;
    const ccomplex_t IT_0203 = IT_0025*IT_0168;
    const ccomplex_t IT_0204 = IT_0002*IT_0095*IT_0203;
    const ccomplex_t IT_0205 = IT_0063*IT_0179;
    const ccomplex_t IT_0206 = IT_0048*IT_0109*IT_0205;
    const ccomplex_t IT_0207 = IT_0082*IT_0143;
    const ccomplex_t IT_0208 = IT_0067*IT_0123*IT_0207;
    const ccomplex_t IT_0209 = IT_0044*IT_0190;
    const ccomplex_t IT_0210 = IT_0029*IT_0148*IT_0209;
    const ccomplex_t IT_0211 = m_C_1*IT_0028;
    const ccomplex_t IT_0212 = m_C_1*IT_0066;
    const ccomplex_t IT_0213 = m_C_1*IT_0047;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0004;
    const ccomplex_t IT_0215 = 0.5*IT_0214;
    const ccomplex_t IT_0216 = e_em*conj(U_Wm1)*U_stau_00;
    const ccomplex_t IT_0217 = IT_0004*IT_0216;
    const ccomplex_t IT_0218 = conj(U_d1)*e_em*m_tau*IT_0020*U_stau_10;
    const ccomplex_t IT_0219 = IT_0019*IT_0218;
    const ccomplex_t IT_0220 = 1.4142135623731*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0217 + (-0.5)*IT_0220);
    const ccomplex_t IT_0222 = IT_0215*IT_0221;
    const ccomplex_t IT_0223 = pow(m_C_1, 2);
    const ccomplex_t IT_0224 = cpow(s_12 + 0.5*IT_0001 + 0.5*IT_0223 + 0.5
      *reg_prop, -1);
    const ccomplex_t IT_0225 = IT_0222*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = m_tau*IT_0226;
    const ccomplex_t IT_0228 = 0.5*IT_0202 + (-0.5)*IT_0204 + (-0.5)*IT_0206 +
       (-0.5)*IT_0208 + (-0.5)*IT_0210 + 0.5*IT_0211 + 0.5*IT_0212 + 0.5*IT_0213
       + 0.5*IT_0227;
    const ccomplex_t IT_0229 = m_C_1*m_tau;
    const ccomplex_t IT_0230 = IT_0088*IT_0229;
    const ccomplex_t IT_0231 = 3*IT_0230;
    const ccomplex_t IT_0232 = (-0.5)*IT_0226;
    const ccomplex_t IT_0233 = (-3)*IT_0230;
    const ccomplex_t IT_0234 = (-0.5)*IT_0146 + (-0.5)*IT_0171 + (-0.5)
      *IT_0182 + (-0.5)*IT_0193;
    const ccomplex_t IT_0235 = -IT_0028 + -IT_0047 + -IT_0066 + -IT_0085;
    const ccomplex_t IT_0236 = s_13*s_14;
    const ccomplex_t IT_0237 = s_34*IT_0223;
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = IT_0236 + IT_0238;
    const ccomplex_t IT_0240 = s_13*IT_0088;
    const ccomplex_t IT_0241 = pow(m_W, 4);
    const ccomplex_t IT_0242 = s_13*IT_0241;
    const ccomplex_t IT_0243 = IT_0087*IT_0242;
    const ccomplex_t IT_0244 = (-0.5)*IT_0243;
    const ccomplex_t IT_0245 = s_14*s_34;
    const ccomplex_t IT_0246 = (-2)*IT_0245;
    const ccomplex_t IT_0247 = IT_0240 + IT_0244 + IT_0246;
    const ccomplex_t IT_0248 = (-4)*IT_0247;
    const ccomplex_t IT_0249 = -IT_0226;
    const ccomplex_t IT_0250 = s_14*s_34*IT_0087;
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = s_13 + IT_0251;
    const ccomplex_t IT_0253 = IT_0229*IT_0252;
    const ccomplex_t IT_0254 = pow(s_34, 2);
    const ccomplex_t IT_0255 = IT_0087*IT_0254;
    const ccomplex_t IT_0256 = -IT_0000;
    const ccomplex_t IT_0257 = IT_0255 + IT_0256;
    const ccomplex_t IT_0258 = s_13*IT_0257;
    const ccomplex_t IT_0259 = m_C_1*IT_0257;
    const ccomplex_t IT_0260 = s_13*s_34;
    const ccomplex_t IT_0261 = s_14*IT_0000;
    const ccomplex_t IT_0262 = -IT_0261;
    const ccomplex_t IT_0263 = IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = (-2)*IT_0263;
    const ccomplex_t IT_0265 = IT_0146 + IT_0171 + IT_0182 + IT_0193;
    const ccomplex_t IT_0266 = s_13*IT_0252;
    const ccomplex_t IT_0267 = (-2)*IT_0266;
    const ccomplex_t IT_0268 = m_tau*IT_0252;
    const ccomplex_t IT_0269 = (-2)*IT_0268;
    const ccomplex_t IT_0270 = (-2)*IT_0249*(conj(IT_0235)*IT_0253 + -conj
      (IT_0249)*IT_0258 + -conj(IT_0195)*IT_0259 + -conj(IT_0232)*IT_0263 + (
      -0.5)*conj(IT_0234)*IT_0264 + (-0.5)*conj(IT_0265)*IT_0267 + (-0.5)*conj
      (IT_0228)*IT_0269);
    const ccomplex_t IT_0271 = pow(s_14, 2);
    const ccomplex_t IT_0272 = IT_0087*IT_0271;
    const ccomplex_t IT_0273 = -IT_0223;
    const ccomplex_t IT_0274 = IT_0272 + IT_0273;
    const ccomplex_t IT_0275 = m_tau*IT_0274;
    const ccomplex_t IT_0276 = IT_0240 + IT_0246;
    const ccomplex_t IT_0277 = IT_0087*IT_0276;
    const ccomplex_t IT_0278 = (-0.25)*IT_0090;
    const ccomplex_t IT_0279 = s_14 + IT_0278;
    const ccomplex_t IT_0280 = m_tau*IT_0279;
    const ccomplex_t IT_0281 = (-4)*IT_0280;
    const ccomplex_t IT_0282 = m_C_1*IT_0252;
    const ccomplex_t IT_0283 = (-2)*IT_0282;
    const ccomplex_t IT_0284 = 0.5*conj(IT_0249);
    const ccomplex_t IT_0285 = (-6)*IT_0229;
    const ccomplex_t IT_0286 = 0.5*conj(IT_0195);
    const ccomplex_t IT_0287 = 2*IT_0228*(conj(IT_0086)*IT_0201 + IT_0094*conj
      (IT_0234) + conj(IT_0265)*IT_0275 + 2*conj(IT_0228)*(s_13 + (-0.5)*IT_0277
      ) + conj(IT_0232)*IT_0281 + 0.5*conj(IT_0235)*IT_0283 + IT_0269*IT_0284 +
       IT_0285*IT_0286);
    const ccomplex_t IT_0288 = s_13*IT_0274;
    const ccomplex_t IT_0289 = IT_0229*IT_0274;
    const ccomplex_t IT_0290 = 2*IT_0289;
    const ccomplex_t IT_0291 = 0.5*IT_0197;
    const ccomplex_t IT_0292 = s_34 + IT_0291;
    const ccomplex_t IT_0293 = m_C_1*IT_0292;
    const ccomplex_t IT_0294 = 2*IT_0293;
    const ccomplex_t IT_0295 = 0.5*conj(IT_0265);
    const ccomplex_t IT_0296 = 0.5*conj(IT_0228);
    const ccomplex_t IT_0297 = 2*IT_0195*(conj(IT_0086)*IT_0094 + IT_0201*conj
      (IT_0234) + conj(IT_0249)*IT_0259 + conj(IT_0235)*IT_0275 + 2*conj(IT_0195
      )*(s_13 + (-0.5)*IT_0277) + conj(IT_0232)*IT_0294 + IT_0283*IT_0295 +
       IT_0285*IT_0296);
    const ccomplex_t IT_0298 = 0.5*conj(IT_0232);
    const ccomplex_t IT_0299 = (-0.25)*IT_0243;
    const ccomplex_t IT_0300 = IT_0245 + IT_0299;
    const ccomplex_t IT_0301 = (-8)*IT_0300;
    const ccomplex_t IT_0302 = 0.5*conj(IT_0234);
    const ccomplex_t IT_0303 = 2*IT_0086*(IT_0094*conj(IT_0195) + IT_0201*conj
      (IT_0228) + IT_0231*conj(IT_0232) + IT_0233*conj(IT_0234) + conj(IT_0235)
      *IT_0239 + 0.5*conj(IT_0086)*IT_0248) + IT_0270 + IT_0287 + (-2)*IT_0265*(
      (conj(IT_0232) + -conj(IT_0234))*IT_0239 + -conj(IT_0228)*IT_0275 + 
      -IT_0267*IT_0284 + -IT_0283*IT_0286 + -conj(IT_0265)*IT_0288 + (-0.5)*conj
      (IT_0235)*IT_0290) + 2*IT_0235*(conj(IT_0086)*IT_0239 + -conj(IT_0249)
      *IT_0253 + conj(IT_0195)*IT_0275 + conj(IT_0235)*IT_0288 + IT_0290*IT_0295
       + IT_0283*IT_0296) + IT_0297 + 2*IT_0234*(conj(IT_0195)*IT_0201 + IT_0094
      *conj(IT_0228) + conj(IT_0086)*IT_0233 + IT_0239*conj(IT_0265) + IT_0264
      *IT_0284 + IT_0298*IT_0301 + IT_0248*IT_0302) + 2*IT_0232*(conj(IT_0086)
      *IT_0231 + conj(IT_0249)*IT_0263 + -IT_0239*conj(IT_0265) + conj(IT_0228)
      *IT_0281 + conj(IT_0195)*IT_0294 + IT_0248*IT_0298 + IT_0301*IT_0302);
    return create_ccomplex_return(IT_0303);
}

