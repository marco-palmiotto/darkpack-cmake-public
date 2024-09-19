#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_1_anti_sb_1_to_tau_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_1_anti_sb_1_to_tau_anti_b(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_stau_00 = param->U_stau_00;
    const ccomplex_t U_stau_10 = param->U_stau_10;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_stau_1, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B1)*e_em*U_stau_00;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W1)*e_em*U_stau_00;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = conj(N_d1)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0007 + IT_0012 + 
      -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = m_b*N_d1*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0023 = IT_0015*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0029 = IT_0009*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0024 + 0.333333333333333
      *IT_0027 + -IT_0030);
    const ccomplex_t IT_0032 = -IT_0031;
    const ccomplex_t IT_0033 = (-3)*IT_0032;
    const ccomplex_t IT_0034 = 0.166666666666667*IT_0033;
    const ccomplex_t IT_0035 = IT_0021*IT_0034;
    const ccomplex_t IT_0036 = IT_0002*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*IT_0036;
    const ccomplex_t IT_0038 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0039 = conj(N_B2)*e_em*U_stau_00;
    const ccomplex_t IT_0040 = IT_0004*IT_0039;
    const ccomplex_t IT_0041 = 1.4142135623731*IT_0040;
    const ccomplex_t IT_0042 = conj(N_W2)*e_em*U_stau_00;
    const ccomplex_t IT_0043 = IT_0009*IT_0042;
    const ccomplex_t IT_0044 = 1.4142135623731*IT_0043;
    const ccomplex_t IT_0045 = conj(N_d2)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0046 = IT_0015*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0041 + IT_0044 + 
      -IT_0047);
    const ccomplex_t IT_0049 = (-0.5)*IT_0048;
    const ccomplex_t IT_0050 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0051 = IT_0009*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = m_b*N_d2*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0054 = IT_0015*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0057 = IT_0004*IT_0056;
    const ccomplex_t IT_0058 = 1.4142135623731*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*(IT_0052 + -IT_0055 + (
      -0.333333333333333)*IT_0058);
    const ccomplex_t IT_0060 = (-3)*IT_0059;
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0049*IT_0061;
    const ccomplex_t IT_0063 = IT_0038*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0066 = conj(N_B3)*e_em*U_stau_00;
    const ccomplex_t IT_0067 = IT_0004*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W3)*e_em*U_stau_00;
    const ccomplex_t IT_0070 = IT_0009*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_d3)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0073 = IT_0015*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + IT_0071 + 
      -IT_0074);
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0078 = IT_0009*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = m_b*N_d3*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0081 = IT_0015*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0084 = IT_0004*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0079 + -IT_0082 + (
      -0.333333333333333)*IT_0085);
    const ccomplex_t IT_0087 = (-3)*IT_0086;
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = IT_0076*IT_0088;
    const ccomplex_t IT_0090 = IT_0065*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0093 = conj(N_B4)*e_em*U_stau_00;
    const ccomplex_t IT_0094 = IT_0004*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W4)*e_em*U_stau_00;
    const ccomplex_t IT_0097 = IT_0009*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_d4)*e_em*m_tau*IT_0016*U_stau_10;
    const ccomplex_t IT_0100 = IT_0015*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0095 + IT_0098 + 
      -IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0105 = IT_0009*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = m_b*N_d4*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0108 = IT_0015*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0111 = IT_0004*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0106 + -IT_0109 + (
      -0.333333333333333)*IT_0112);
    const ccomplex_t IT_0114 = (-3)*IT_0113;
    const ccomplex_t IT_0115 = 0.166666666666667*IT_0114;
    const ccomplex_t IT_0116 = IT_0103*IT_0115;
    const ccomplex_t IT_0117 = IT_0092*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = 0.5*IT_0037 + 0.5*IT_0064 + 0.5*IT_0091 + 0.5
      *IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0121 = N_B2*e_em*U_stau_10;
    const ccomplex_t IT_0122 = IT_0004*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = N_d2*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0125 = IT_0015*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + 0.5*IT_0126);
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = IT_0061*IT_0128;
    const ccomplex_t IT_0130 = IT_0038*IT_0120*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0132 = N_B1*e_em*U_stau_10;
    const ccomplex_t IT_0133 = IT_0004*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = N_d1*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0136 = IT_0015*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + 0.5*IT_0137);
    const ccomplex_t IT_0139 = -IT_0138;
    const ccomplex_t IT_0140 = IT_0034*IT_0139;
    const ccomplex_t IT_0141 = IT_0002*IT_0131*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0143 = N_B4*e_em*U_stau_10;
    const ccomplex_t IT_0144 = IT_0004*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = N_d4*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0147 = IT_0015*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + 0.5*IT_0148);
    const ccomplex_t IT_0150 = -IT_0149;
    const ccomplex_t IT_0151 = IT_0115*IT_0150;
    const ccomplex_t IT_0152 = IT_0092*IT_0142*IT_0151;
    const ccomplex_t IT_0153 = m_tau*IT_0037;
    const ccomplex_t IT_0154 = m_tau*IT_0064;
    const ccomplex_t IT_0155 = m_tau*IT_0091;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0157 = N_B3*e_em*U_stau_10;
    const ccomplex_t IT_0158 = IT_0004*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = N_d3*e_em*m_tau*IT_0016*U_stau_00;
    const ccomplex_t IT_0161 = IT_0015*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0159 + 0.5*IT_0162);
    const ccomplex_t IT_0164 = -IT_0163;
    const ccomplex_t IT_0165 = IT_0088*IT_0164;
    const ccomplex_t IT_0166 = IT_0065*IT_0156*IT_0165;
    const ccomplex_t IT_0167 = m_tau*IT_0118;
    const ccomplex_t IT_0168 = 0.5*IT_0130 + 0.5*IT_0141 + 0.5*IT_0152 + (-0.5
      )*IT_0153 + (-0.5)*IT_0154 + (-0.5)*IT_0155 + 0.5*IT_0166 + (-0.5)*IT_0167;
    const ccomplex_t IT_0169 = s_14*m_tau;
    const ccomplex_t IT_0170 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0171 = IT_0004*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = m_b*conj(N_d2)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0174 = IT_0015*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0172 + 1.5*IT_0175);
    const ccomplex_t IT_0177 = 0.666666666666667*IT_0176;
    const ccomplex_t IT_0178 = 1.5*IT_0177;
    const ccomplex_t IT_0179 = (-0.333333333333333)*IT_0178;
    const ccomplex_t IT_0180 = IT_0049*IT_0179;
    const ccomplex_t IT_0181 = IT_0038*IT_0120*IT_0180;
    const ccomplex_t IT_0182 = m_b*conj(N_d4)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0183 = IT_0015*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0186 = IT_0004*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0184 + 0.666666666666667
      *IT_0187);
    const ccomplex_t IT_0189 = 1.5*IT_0188;
    const ccomplex_t IT_0190 = (-0.333333333333333)*IT_0189;
    const ccomplex_t IT_0191 = IT_0103*IT_0190;
    const ccomplex_t IT_0192 = IT_0092*IT_0142*IT_0191;
    const ccomplex_t IT_0193 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0194 = IT_0004*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = m_b*conj(N_d1)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0197 = IT_0015*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + 1.5*IT_0198);
    const ccomplex_t IT_0200 = (-0.333333333333333)*IT_0199;
    const ccomplex_t IT_0201 = IT_0139*IT_0200;
    const ccomplex_t IT_0202 = IT_0002*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = m_tau*IT_0203;
    const ccomplex_t IT_0205 = m_b*conj(N_d3)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0206 = IT_0015*IT_0205;
    const ccomplex_t IT_0207 = 1.4142135623731*IT_0206;
    const ccomplex_t IT_0208 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0209 = IT_0004*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0207 + 0.666666666666667
      *IT_0210);
    const ccomplex_t IT_0212 = 1.5*IT_0211;
    const ccomplex_t IT_0213 = (-0.333333333333333)*IT_0212;
    const ccomplex_t IT_0214 = IT_0076*IT_0213;
    const ccomplex_t IT_0215 = IT_0065*IT_0156*IT_0214;
    const ccomplex_t IT_0216 = IT_0164*IT_0213;
    const ccomplex_t IT_0217 = IT_0065*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = m_tau*IT_0218;
    const ccomplex_t IT_0220 = IT_0021*IT_0200;
    const ccomplex_t IT_0221 = IT_0002*IT_0131*IT_0220;
    const ccomplex_t IT_0222 = IT_0128*IT_0179;
    const ccomplex_t IT_0223 = IT_0038*IT_0222;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = m_tau*IT_0224;
    const ccomplex_t IT_0226 = IT_0150*IT_0190;
    const ccomplex_t IT_0227 = IT_0092*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = m_tau*IT_0228;
    const ccomplex_t IT_0230 = 0.5*IT_0181 + 0.5*IT_0192 + (-0.5)*IT_0204 +
       0.5*IT_0215 + (-0.5)*IT_0219 + 0.5*IT_0221 + (-0.5)*IT_0225 + (-0.5)
      *IT_0229;
    const ccomplex_t IT_0231 = m_b*s_13;
    const ccomplex_t IT_0232 = 0.5*IT_0203 + 0.5*IT_0218 + 0.5*IT_0224 + 0.5
      *IT_0228;
    const ccomplex_t IT_0233 = m_b*m_tau*IT_0001;
    const ccomplex_t IT_0234 = s_34*IT_0001;
    const ccomplex_t IT_0235 = (-6)*IT_0234;
    const ccomplex_t IT_0236 = s_13*s_14;
    const ccomplex_t IT_0237 = 12*IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = m_b*m_tau;
    const ccomplex_t IT_0240 = 6*s_34;
    const ccomplex_t IT_0241 = -conj(IT_0119);
    const ccomplex_t IT_0242 = (-6)*conj(IT_0168);
    const ccomplex_t IT_0243 = 6*IT_0119*(conj(IT_0168)*IT_0169 + -conj
      (IT_0230)*IT_0231 + -conj(IT_0232)*IT_0233 + 0.166666666666667*conj
      (IT_0119)*IT_0238) + (-6)*IT_0168*(IT_0231*conj(IT_0232) + conj(IT_0230)
      *IT_0239 + (-0.166666666666667)*conj(IT_0168)*IT_0240 + IT_0169*IT_0241) +
       6*IT_0232*(IT_0169*conj(IT_0230) + 0.166666666666667*conj(IT_0232)
      *IT_0238 + IT_0233*IT_0241 + 0.166666666666667*IT_0231*IT_0242) + 6
      *IT_0230*(IT_0169*conj(IT_0232) + 0.166666666666667*conj(IT_0230)*IT_0240 
      + IT_0231*IT_0241 + 0.166666666666667*IT_0239*IT_0242);
    return create_ccomplex_return(IT_0243);
}

