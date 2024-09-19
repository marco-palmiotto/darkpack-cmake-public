#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_stau_2_anti_sb_1_to_tau_anti_b.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_stau_2_anti_sb_1_to_tau_anti_b(
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
    const ccomplex_t U_sb_00 = param->U_sb_00;
    const ccomplex_t U_sb_10 = param->U_sb_10;
    const ccomplex_t U_stau_01 = param->U_stau_01;
    const ccomplex_t U_stau_11 = param->U_stau_11;
    const ccomplex_t IT_0000 = pow(m_tau, 2);
    const ccomplex_t IT_0001 = pow(m_stau_2, 2);
    const ccomplex_t IT_0002 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_1*(m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0003 = cos(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_B1)*e_em*U_stau_01;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = 1.4142135623731*IT_0006;
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = conj(N_W1)*e_em*U_stau_01;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = 1.4142135623731*IT_0011;
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0009*IT_0014;
    const ccomplex_t IT_0016 = pow(m_W, -1);
    const ccomplex_t IT_0017 = conj(N_d1)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0018 = IT_0015*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*(IT_0007 + IT_0012 + 
      -IT_0019);
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = m_b*N_d1*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0023 = IT_0015*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_W1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0026 = IT_0009*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = N_B1*e_em*conj(U_sb_00);
    const ccomplex_t IT_0029 = IT_0004*IT_0028;
    const ccomplex_t IT_0030 = 1.4142135623731*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*(IT_0024 + -IT_0027 +
       0.333333333333333*IT_0030);
    const ccomplex_t IT_0032 = 3*IT_0031;
    const ccomplex_t IT_0033 = 0.166666666666667*IT_0032;
    const ccomplex_t IT_0034 = IT_0021*IT_0033;
    const ccomplex_t IT_0035 = IT_0002*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_2*(m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0038 = conj(N_B2)*e_em*U_stau_01;
    const ccomplex_t IT_0039 = IT_0004*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W2)*e_em*U_stau_01;
    const ccomplex_t IT_0042 = IT_0009*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = conj(N_d2)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0045 = IT_0015*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0040 + IT_0043 + 
      -IT_0046);
    const ccomplex_t IT_0048 = (-0.5)*IT_0047;
    const ccomplex_t IT_0049 = m_b*N_d2*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0050 = IT_0015*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = N_W2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0053 = IT_0009*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = N_B2*e_em*conj(U_sb_00);
    const ccomplex_t IT_0056 = IT_0004*IT_0055;
    const ccomplex_t IT_0057 = 1.4142135623731*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*(IT_0051 + -IT_0054 +
       0.333333333333333*IT_0057);
    const ccomplex_t IT_0059 = -IT_0058;
    const ccomplex_t IT_0060 = (-3)*IT_0059;
    const ccomplex_t IT_0061 = 0.166666666666667*IT_0060;
    const ccomplex_t IT_0062 = IT_0048*IT_0061;
    const ccomplex_t IT_0063 = IT_0037*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_3*(m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0066 = conj(N_B3)*e_em*U_stau_01;
    const ccomplex_t IT_0067 = IT_0004*IT_0066;
    const ccomplex_t IT_0068 = 1.4142135623731*IT_0067;
    const ccomplex_t IT_0069 = conj(N_W3)*e_em*U_stau_01;
    const ccomplex_t IT_0070 = IT_0009*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = conj(N_d3)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0073 = IT_0015*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0068 + IT_0071 + 
      -IT_0074);
    const ccomplex_t IT_0076 = (-0.5)*IT_0075;
    const ccomplex_t IT_0077 = m_b*N_d3*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0078 = IT_0015*IT_0077;
    const ccomplex_t IT_0079 = 1.4142135623731*IT_0078;
    const ccomplex_t IT_0080 = N_B3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0081 = IT_0004*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_W3*e_em*conj(U_sb_00);
    const ccomplex_t IT_0084 = IT_0009*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0079 + 0.333333333333333
      *IT_0082 + -IT_0085);
    const ccomplex_t IT_0087 = 3*IT_0086;
    const ccomplex_t IT_0088 = 0.166666666666667*IT_0087;
    const ccomplex_t IT_0089 = IT_0076*IT_0088;
    const ccomplex_t IT_0090 = IT_0065*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow(s_13 + (-0.5)*IT_0000 + (-0.5)*IT_0001 +
       0.5*m_N_4*(m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + (-0.5)*reg_prop, -1);
    const ccomplex_t IT_0093 = conj(N_B4)*e_em*U_stau_01;
    const ccomplex_t IT_0094 = IT_0004*IT_0093;
    const ccomplex_t IT_0095 = 1.4142135623731*IT_0094;
    const ccomplex_t IT_0096 = conj(N_W4)*e_em*U_stau_01;
    const ccomplex_t IT_0097 = IT_0009*IT_0096;
    const ccomplex_t IT_0098 = 1.4142135623731*IT_0097;
    const ccomplex_t IT_0099 = conj(N_d4)*e_em*m_tau*IT_0016*U_stau_11;
    const ccomplex_t IT_0100 = IT_0015*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0095 + IT_0098 + 
      -IT_0101);
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = N_B4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0105 = IT_0004*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = N_W4*e_em*conj(U_sb_00);
    const ccomplex_t IT_0108 = IT_0009*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = m_b*N_d4*e_em*IT_0016*conj(U_sb_10);
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = 1.4142135623731*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*(IT_0106 + (-3)*IT_0109 + 3
      *IT_0112);
    const ccomplex_t IT_0114 = 0.166666666666667*IT_0113;
    const ccomplex_t IT_0115 = IT_0103*IT_0114;
    const ccomplex_t IT_0116 = IT_0092*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = 0.5*IT_0036 + 0.5*IT_0064 + 0.5*IT_0091 + 0.5
      *IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0120 = N_B1*e_em*U_stau_11;
    const ccomplex_t IT_0121 = IT_0004*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = N_d1*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0124 = IT_0015*IT_0123;
    const ccomplex_t IT_0125 = 1.4142135623731*IT_0124;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*(IT_0122 + 0.5*IT_0125);
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0033*IT_0127;
    const ccomplex_t IT_0129 = IT_0002*IT_0119*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0131 = N_B3*e_em*U_stau_11;
    const ccomplex_t IT_0132 = IT_0004*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = N_d3*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0135 = IT_0015*IT_0134;
    const ccomplex_t IT_0136 = 1.4142135623731*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + 0.5*IT_0136);
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0088*IT_0138;
    const ccomplex_t IT_0140 = IT_0065*IT_0130*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0142 = N_B4*e_em*U_stau_11;
    const ccomplex_t IT_0143 = IT_0004*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = N_d4*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0146 = IT_0015*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0144 + 0.5*IT_0147);
    const ccomplex_t IT_0149 = -IT_0148;
    const ccomplex_t IT_0150 = IT_0114*IT_0149;
    const ccomplex_t IT_0151 = IT_0092*IT_0141*IT_0150;
    const ccomplex_t IT_0152 = m_tau*IT_0036;
    const ccomplex_t IT_0153 = m_tau*IT_0064;
    const ccomplex_t IT_0154 = m_tau*IT_0091;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0156 = N_B2*e_em*U_stau_11;
    const ccomplex_t IT_0157 = IT_0004*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = N_d2*e_em*m_tau*IT_0016*U_stau_01;
    const ccomplex_t IT_0160 = IT_0015*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0158 + 0.5*IT_0161);
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = IT_0061*IT_0163;
    const ccomplex_t IT_0165 = IT_0037*IT_0155*IT_0164;
    const ccomplex_t IT_0166 = m_tau*IT_0117;
    const ccomplex_t IT_0167 = 0.5*IT_0129 + 0.5*IT_0140 + 0.5*IT_0151 + (-0.5
      )*IT_0152 + (-0.5)*IT_0153 + (-0.5)*IT_0154 + 0.5*IT_0165 + (-0.5)*IT_0166;
    const ccomplex_t IT_0168 = s_14*m_tau;
    const ccomplex_t IT_0169 = conj(N_B3)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0170 = IT_0004*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = m_b*conj(N_d3)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0173 = IT_0015*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + 1.5*IT_0174);
    const ccomplex_t IT_0176 = (-0.333333333333333)*IT_0175;
    const ccomplex_t IT_0177 = IT_0076*IT_0176;
    const ccomplex_t IT_0178 = IT_0065*IT_0130*IT_0177;
    const ccomplex_t IT_0179 = conj(N_B2)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0180 = IT_0004*IT_0179;
    const ccomplex_t IT_0181 = 1.4142135623731*IT_0180;
    const ccomplex_t IT_0182 = m_b*conj(N_d2)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0183 = IT_0015*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = (0 + _Complex_I*1)*(IT_0181 + 1.5*IT_0184);
    const ccomplex_t IT_0186 = (-0.333333333333333)*IT_0185;
    const ccomplex_t IT_0187 = IT_0048*IT_0186;
    const ccomplex_t IT_0188 = IT_0037*IT_0155*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B4)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0190 = IT_0004*IT_0189;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = m_b*conj(N_d4)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0193 = IT_0015*IT_0192;
    const ccomplex_t IT_0194 = 1.4142135623731*IT_0193;
    const ccomplex_t IT_0195 = (0 + _Complex_I*1)*(IT_0191 + 1.5*IT_0194);
    const ccomplex_t IT_0196 = (-0.333333333333333)*IT_0195;
    const ccomplex_t IT_0197 = IT_0103*IT_0196;
    const ccomplex_t IT_0198 = IT_0092*IT_0141*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B1)*e_em*conj(U_sb_10);
    const ccomplex_t IT_0200 = IT_0004*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = m_b*conj(N_d1)*e_em*IT_0016*conj(U_sb_00);
    const ccomplex_t IT_0203 = IT_0015*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0201 + 1.5*IT_0204);
    const ccomplex_t IT_0206 = (-0.333333333333333)*IT_0205;
    const ccomplex_t IT_0207 = IT_0127*IT_0206;
    const ccomplex_t IT_0208 = IT_0002*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = m_tau*IT_0209;
    const ccomplex_t IT_0211 = IT_0163*IT_0186;
    const ccomplex_t IT_0212 = IT_0037*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = m_tau*IT_0213;
    const ccomplex_t IT_0215 = IT_0138*IT_0176;
    const ccomplex_t IT_0216 = IT_0065*IT_0215;
    const ccomplex_t IT_0217 = (0 + _Complex_I*1)*IT_0216;
    const ccomplex_t IT_0218 = m_tau*IT_0217;
    const ccomplex_t IT_0219 = IT_0149*IT_0196;
    const ccomplex_t IT_0220 = IT_0092*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = m_tau*IT_0221;
    const ccomplex_t IT_0223 = IT_0021*IT_0206;
    const ccomplex_t IT_0224 = IT_0002*IT_0119*IT_0223;
    const ccomplex_t IT_0225 = 0.5*IT_0178 + 0.5*IT_0188 + 0.5*IT_0198 + (-0.5
      )*IT_0210 + (-0.5)*IT_0214 + (-0.5)*IT_0218 + (-0.5)*IT_0222 + 0.5*IT_0224;
    const ccomplex_t IT_0226 = m_b*s_13;
    const ccomplex_t IT_0227 = 0.5*IT_0209 + 0.5*IT_0213 + 0.5*IT_0217 + 0.5
      *IT_0221;
    const ccomplex_t IT_0228 = m_b*m_tau*IT_0001;
    const ccomplex_t IT_0229 = s_34*IT_0001;
    const ccomplex_t IT_0230 = (-6)*IT_0229;
    const ccomplex_t IT_0231 = s_13*s_14;
    const ccomplex_t IT_0232 = 12*IT_0231;
    const ccomplex_t IT_0233 = IT_0230 + IT_0232;
    const ccomplex_t IT_0234 = m_b*m_tau;
    const ccomplex_t IT_0235 = 6*s_34;
    const ccomplex_t IT_0236 = -conj(IT_0118);
    const ccomplex_t IT_0237 = (-6)*conj(IT_0167);
    const ccomplex_t IT_0238 = 6*IT_0118*(conj(IT_0167)*IT_0168 + -conj
      (IT_0225)*IT_0226 + -conj(IT_0227)*IT_0228 + 0.166666666666667*conj
      (IT_0118)*IT_0233) + (-6)*IT_0167*(IT_0226*conj(IT_0227) + conj(IT_0225)
      *IT_0234 + (-0.166666666666667)*conj(IT_0167)*IT_0235 + IT_0168*IT_0236) +
       6*IT_0227*(IT_0168*conj(IT_0225) + 0.166666666666667*conj(IT_0227)
      *IT_0233 + IT_0228*IT_0236 + 0.166666666666667*IT_0226*IT_0237) + 6
      *IT_0225*(IT_0168*conj(IT_0227) + 0.166666666666667*conj(IT_0225)*IT_0235 
      + IT_0226*IT_0236 + 0.166666666666667*IT_0234*IT_0237);
    return create_ccomplex_return(IT_0238);
}

