#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_C_1_to_Hp_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_C_1_to_Hp_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
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
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + (-0.5)*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = cos(beta);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0012 = IT_0010*IT_0011;
    const ccomplex_t IT_0013 = 2*IT_0012;
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0021 = IT_0010*IT_0020;
    const ccomplex_t IT_0022 = 1.4142135623731*IT_0021;
    const ccomplex_t IT_0023 = IT_0019 + IT_0022;
    const ccomplex_t IT_0024 = IT_0013 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = pow(m_Hp, 2);
    const ccomplex_t IT_0029 = pow(m_C_1, 2);
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0027*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0034 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0035 = IT_0001*IT_0034;
    const ccomplex_t IT_0036 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0037 = IT_0001*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + (-0.5)*IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0042 = IT_0016*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0045 = IT_0010*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = IT_0043 + IT_0046;
    const ccomplex_t IT_0048 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0049 = IT_0010*IT_0048;
    const ccomplex_t IT_0050 = 2*IT_0049;
    const ccomplex_t IT_0051 = IT_0047 + IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0040*IT_0053;
    const ccomplex_t IT_0055 = IT_0033*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0058 = IT_0001*IT_0057;
    const ccomplex_t IT_0059 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0060 = IT_0001*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-0.5)*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0065 = IT_0010*IT_0064;
    const ccomplex_t IT_0066 = 2*IT_0065;
    const ccomplex_t IT_0067 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0068 = IT_0016*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0071 = IT_0010*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = IT_0069 + IT_0072;
    const ccomplex_t IT_0074 = IT_0066 + IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0063*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0082 = IT_0001*IT_0081;
    const ccomplex_t IT_0083 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0084 = IT_0001*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-0.5)*IT_0085);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0089 = IT_0010*IT_0088;
    const ccomplex_t IT_0090 = 2*IT_0089;
    const ccomplex_t IT_0091 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0092 = IT_0016*IT_0091;
    const ccomplex_t IT_0093 = 1.4142135623731*IT_0092;
    const ccomplex_t IT_0094 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0095 = IT_0010*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = IT_0093 + IT_0096;
    const ccomplex_t IT_0098 = IT_0090 + IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = 0.5*IT_0099;
    const ccomplex_t IT_0101 = IT_0087*IT_0100;
    const ccomplex_t IT_0102 = cpow((-2)*s_23 + IT_0028 + IT_0029 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0103 = IT_0101*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = 2*IT_0032 + 2*IT_0056 + 2*IT_0080 + 2*IT_0104;
    const ccomplex_t IT_0106 = pow(m_W, -2);
    const ccomplex_t IT_0107 = pow(s_14, 2);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = -IT_0029;
    const ccomplex_t IT_0110 = IT_0108 + IT_0109;
    const ccomplex_t IT_0111 = s_12*IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0054*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0116 = IT_0101*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1)*IT_0116;
    const ccomplex_t IT_0118 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0119 = IT_0027*IT_0118;
    const ccomplex_t IT_0120 = (0 + _Complex_I*1)*IT_0119;
    const ccomplex_t IT_0121 = cpow((-2)*s_13 + IT_0028 + IT_0029 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0122 = IT_0077*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = (-2)*IT_0114 + (-2)*IT_0117 + (-2)*IT_0120 + (
      -2)*IT_0123;
    const ccomplex_t IT_0125 = s_14*s_24*IT_0106;
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = s_12 + IT_0126;
    const ccomplex_t IT_0128 = s_12*IT_0127;
    const ccomplex_t IT_0129 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0130 = IT_0001*IT_0129;
    const ccomplex_t IT_0131 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0132 = IT_0001*IT_0131;
    const ccomplex_t IT_0133 = 1.4142135623731*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*(IT_0130 + 0.5*IT_0133);
    const ccomplex_t IT_0135 = -IT_0134;
    const ccomplex_t IT_0136 = sin(beta);
    const ccomplex_t IT_0137 = IT_0015*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = 1.4142135623731*IT_0139;
    const ccomplex_t IT_0141 = IT_0001*IT_0136;
    const ccomplex_t IT_0142 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0143 = IT_0141*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = IT_0140 + IT_0144;
    const ccomplex_t IT_0146 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0147 = IT_0141*IT_0146;
    const ccomplex_t IT_0148 = (-2)*IT_0147;
    const ccomplex_t IT_0149 = IT_0145 + IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0135*IT_0151;
    const ccomplex_t IT_0153 = IT_0030*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0156 = IT_0141*IT_0155;
    const ccomplex_t IT_0157 = (-2)*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0159 = IT_0137*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0162 = IT_0141*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = IT_0160 + IT_0163;
    const ccomplex_t IT_0165 = IT_0157 + IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = 0.5*IT_0166;
    const ccomplex_t IT_0168 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0169 = IT_0001*IT_0168;
    const ccomplex_t IT_0170 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0171 = IT_0001*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0169 + 0.5*IT_0172);
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = IT_0167*IT_0174;
    const ccomplex_t IT_0176 = IT_0033*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*IT_0176;
    const ccomplex_t IT_0178 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0179 = IT_0141*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0182 = IT_0137*IT_0181;
    const ccomplex_t IT_0183 = 1.4142135623731*IT_0182;
    const ccomplex_t IT_0184 = IT_0180 + IT_0183;
    const ccomplex_t IT_0185 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0186 = IT_0141*IT_0185;
    const ccomplex_t IT_0187 = (-2)*IT_0186;
    const ccomplex_t IT_0188 = IT_0184 + IT_0187;
    const ccomplex_t IT_0189 = (0 + _Complex_I*1)*IT_0188;
    const ccomplex_t IT_0190 = 0.5*IT_0189;
    const ccomplex_t IT_0191 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0192 = IT_0001*IT_0191;
    const ccomplex_t IT_0193 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0194 = IT_0001*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0192 + 0.5*IT_0195);
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = IT_0190*IT_0197;
    const ccomplex_t IT_0199 = IT_0078*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0202 = IT_0141*IT_0201;
    const ccomplex_t IT_0203 = (-2)*IT_0202;
    const ccomplex_t IT_0204 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0205 = IT_0137*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0208 = IT_0141*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = IT_0206 + IT_0209;
    const ccomplex_t IT_0211 = IT_0203 + IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = 0.5*IT_0212;
    const ccomplex_t IT_0214 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0215 = IT_0001*IT_0214;
    const ccomplex_t IT_0216 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0217 = IT_0001*IT_0216;
    const ccomplex_t IT_0218 = 1.4142135623731*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*(IT_0215 + 0.5*IT_0218);
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = IT_0213*IT_0220;
    const ccomplex_t IT_0222 = IT_0102*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = 2*IT_0154 + 2*IT_0177 + 2*IT_0200 + 2*IT_0223;
    const ccomplex_t IT_0225 = IT_0029*IT_0110;
    const ccomplex_t IT_0226 = IT_0121*IT_0198;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = IT_0118*IT_0152;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = IT_0112*IT_0175;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = IT_0115*IT_0221;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = (-2)*IT_0227 + (-2)*IT_0229 + (-2)*IT_0231 + (
      -2)*IT_0233;
    const ccomplex_t IT_0235 = IT_0029*IT_0127;
    const ccomplex_t IT_0236 = 2*IT_0105;
    const ccomplex_t IT_0237 = 2*conj(IT_0105);
    const ccomplex_t IT_0238 = 2*IT_0224;
    const ccomplex_t IT_0239 = 2*conj(IT_0224);
    const ccomplex_t IT_0240 = 2*IT_0124;
    const ccomplex_t IT_0241 = 2*conj(IT_0124);
    const ccomplex_t IT_0242 = pow(s_24, 2);
    const ccomplex_t IT_0243 = IT_0106*IT_0242;
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = IT_0029 + IT_0244;
    const ccomplex_t IT_0246 = IT_0029*IT_0245;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = s_12*IT_0245;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = 2*IT_0234;
    const ccomplex_t IT_0251 = IT_0114 + IT_0117 + IT_0120 + IT_0123 + IT_0154
       + IT_0177 + IT_0200 + IT_0223;
    const ccomplex_t IT_0252 = m_C_1*IT_0251;
    const ccomplex_t IT_0253 = pow(m_W, 2);
    const ccomplex_t IT_0254 = s_24*IT_0253;
    const ccomplex_t IT_0255 = IT_0106*IT_0254;
    const ccomplex_t IT_0256 = 0.5*IT_0255;
    const ccomplex_t IT_0257 = s_24 + IT_0256;
    const ccomplex_t IT_0258 = m_C_1*IT_0257;
    const ccomplex_t IT_0259 = 2*IT_0258;
    const ccomplex_t IT_0260 = IT_0227 + IT_0229 + IT_0231 + IT_0233;
    const ccomplex_t IT_0261 = 2*IT_0260;
    const ccomplex_t IT_0262 = IT_0114 + IT_0117 + IT_0120 + IT_0123;
    const ccomplex_t IT_0263 = IT_0032 + IT_0056 + IT_0080 + IT_0104 + IT_0227
       + IT_0229 + IT_0231 + IT_0233;
    const ccomplex_t IT_0264 = m_C_1*IT_0263;
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0268 = IT_0008*IT_0151;
    const ccomplex_t IT_0269 = IT_0118*IT_0267*IT_0268;
    const ccomplex_t IT_0270 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0271 = IT_0053*IT_0174;
    const ccomplex_t IT_0272 = IT_0033*IT_0270*IT_0271;
    const ccomplex_t IT_0273 = IT_0026*IT_0135;
    const ccomplex_t IT_0274 = IT_0030*IT_0267*IT_0273;
    const ccomplex_t IT_0275 = IT_0040*IT_0167;
    const ccomplex_t IT_0276 = IT_0112*IT_0270*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0278 = IT_0076*IT_0197;
    const ccomplex_t IT_0279 = IT_0078*IT_0277*IT_0278;
    const ccomplex_t IT_0280 = IT_0063*IT_0190;
    const ccomplex_t IT_0281 = IT_0121*IT_0277*IT_0280;
    const ccomplex_t IT_0282 = IT_0100*IT_0220;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0284 = IT_0102*IT_0282*IT_0283;
    const ccomplex_t IT_0285 = IT_0087*IT_0213;
    const ccomplex_t IT_0286 = IT_0115*IT_0283*IT_0285;
    const ccomplex_t IT_0287 = IT_0269 + IT_0272 + IT_0274 + IT_0276 + IT_0279
       + IT_0281 + IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = IT_0266 + IT_0288;
    const ccomplex_t IT_0290 = s_14*IT_0253;
    const ccomplex_t IT_0291 = IT_0106*IT_0290;
    const ccomplex_t IT_0292 = (-0.25)*IT_0291;
    const ccomplex_t IT_0293 = s_14 + IT_0292;
    const ccomplex_t IT_0294 = m_C_1*IT_0293;
    const ccomplex_t IT_0295 = (-4)*IT_0294;
    const ccomplex_t IT_0296 = IT_0030*IT_0267*IT_0268;
    const ccomplex_t IT_0297 = IT_0112*IT_0270*IT_0271;
    const ccomplex_t IT_0298 = IT_0115*IT_0282*IT_0283;
    const ccomplex_t IT_0299 = IT_0102*IT_0283*IT_0285;
    const ccomplex_t IT_0300 = IT_0118*IT_0267*IT_0273;
    const ccomplex_t IT_0301 = IT_0033*IT_0270*IT_0275;
    const ccomplex_t IT_0302 = IT_0121*IT_0277*IT_0278;
    const ccomplex_t IT_0303 = IT_0078*IT_0277*IT_0280;
    const ccomplex_t IT_0304 = -IT_0296 + -IT_0297 + -IT_0298 + -IT_0299 + 
      -IT_0300 + -IT_0301 + -IT_0302 + -IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0304;
    const ccomplex_t IT_0306 = s_12*s_14;
    const ccomplex_t IT_0307 = s_24*IT_0029;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = s_12*s_24;
    const ccomplex_t IT_0311 = s_14*IT_0029;
    const ccomplex_t IT_0312 = -IT_0311;
    const ccomplex_t IT_0313 = IT_0310 + IT_0312;
    const ccomplex_t IT_0314 = IT_0029*IT_0253;
    const ccomplex_t IT_0315 = (-3)*IT_0314;
    const ccomplex_t IT_0316 = s_12*IT_0253;
    const ccomplex_t IT_0317 = s_14*s_24;
    const ccomplex_t IT_0318 = (-2)*IT_0317;
    const ccomplex_t IT_0319 = IT_0316 + IT_0318;
    const ccomplex_t IT_0320 = IT_0106*IT_0319;
    const ccomplex_t IT_0321 = (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = s_12 + IT_0321;
    const ccomplex_t IT_0323 = conj(IT_0252) + conj(IT_0304);
    const ccomplex_t IT_0324 = conj(IT_0265) + conj(IT_0287);
    const ccomplex_t IT_0325 = 4*IT_0265;
    const ccomplex_t IT_0326 = IT_0265 + IT_0287;
    const ccomplex_t IT_0327 = 6*conj(IT_0304);
    const ccomplex_t IT_0328 = pow(m_W, 4);
    const ccomplex_t IT_0329 = s_12*IT_0328;
    const ccomplex_t IT_0330 = IT_0106*IT_0329;
    const ccomplex_t IT_0331 = (-0.5)*IT_0330;
    const ccomplex_t IT_0332 = IT_0316 + IT_0318 + IT_0331;
    const ccomplex_t IT_0333 = -IT_0154 + -IT_0177 + -IT_0200 + -IT_0223;
    const ccomplex_t IT_0334 = (-0.25)*IT_0330;
    const ccomplex_t IT_0335 = IT_0317 + IT_0334;
    const ccomplex_t IT_0336 = (-2)*conj(IT_0265);
    const ccomplex_t IT_0337 = (-2)*conj(IT_0287);
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = 2*conj(IT_0304);
    const ccomplex_t IT_0340 = 2*conj(IT_0234);
    const ccomplex_t IT_0341 = 3*IT_0314;
    const ccomplex_t IT_0342 = -IT_0032 + -IT_0056 + -IT_0080 + -IT_0104;
    const ccomplex_t IT_0343 = (-2)*conj(IT_0342);
    const ccomplex_t IT_0344 = (-4)*IT_0260*((-0.25)*IT_0239*IT_0309 + conj
      (IT_0260)*IT_0332 + (-2)*conj(IT_0333)*IT_0335 + (-0.25)*IT_0295*IT_0338 +
       (-0.25)*IT_0259*IT_0339 + (-0.25)*IT_0313*IT_0340 + (-0.25)*IT_0341
      *IT_0343);
    const ccomplex_t IT_0345 = (-0.25)*IT_0255;
    const ccomplex_t IT_0346 = s_24 + IT_0345;
    const ccomplex_t IT_0347 = m_C_1*IT_0346;
    const ccomplex_t IT_0348 = (-4)*IT_0347;
    const ccomplex_t IT_0349 = 2*IT_0265;
    const ccomplex_t IT_0350 = 2*IT_0287;
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = 0.5*IT_0291;
    const ccomplex_t IT_0353 = s_14 + IT_0352;
    const ccomplex_t IT_0354 = m_C_1*IT_0353;
    const ccomplex_t IT_0355 = 2*IT_0354;
    const ccomplex_t IT_0356 = (-2)*IT_0304;
    const ccomplex_t IT_0357 = m_C_1*IT_0110;
    const ccomplex_t IT_0358 = (-2)*IT_0105;
    const ccomplex_t IT_0359 = (-2)*conj(IT_0105);
    const ccomplex_t IT_0360 = m_C_1*IT_0245;
    const ccomplex_t IT_0361 = -IT_0360;
    const ccomplex_t IT_0362 = 2*conj(IT_0260);
    const ccomplex_t IT_0363 = conj(IT_0262)*IT_0295 + 2*IT_0322*IT_0323 + 3
      *IT_0029*IT_0324 + -conj(IT_0333)*IT_0348 + 0.5*IT_0343*IT_0355 + 0.5
      *IT_0357*IT_0359 + 0.5*IT_0340*IT_0361 + 0.5*IT_0259*IT_0362;
    const ccomplex_t IT_0364 = (-2)*IT_0124;
    const ccomplex_t IT_0365 = (-2)*conj(IT_0124);
    const ccomplex_t IT_0366 = m_C_1*IT_0127;
    const ccomplex_t IT_0367 = (conj(IT_0224)*IT_0305 + IT_0224*IT_0339 +
       IT_0324*(IT_0250 + IT_0358) + IT_0326*(IT_0340 + IT_0359) + conj(IT_0304)
      *IT_0364 + conj(IT_0252)*(IT_0238 + IT_0364) + IT_0304*IT_0365 + IT_0252*
      (IT_0239 + IT_0365))*IT_0366;
    const ccomplex_t IT_0368 = (-2)*conj(IT_0304);
    const ccomplex_t IT_0369 = (-0.5)*IT_0368;
    const ccomplex_t IT_0370 = conj(IT_0252) + IT_0369;
    const ccomplex_t IT_0371 = 2*conj(IT_0265);
    const ccomplex_t IT_0372 = 2*conj(IT_0287);
    const ccomplex_t IT_0373 = IT_0371 + IT_0372;
    const ccomplex_t IT_0374 = (-0.5)*IT_0309;
    const ccomplex_t IT_0375 = (-0.5)*IT_0313;
    const ccomplex_t IT_0376 = 2*IT_0342*(IT_0315*conj(IT_0333) + 4*conj
      (IT_0262)*IT_0335 + -conj(IT_0260)*IT_0341 + (-2)*IT_0332*conj(IT_0342) + 
      -IT_0355*IT_0370 + 0.5*IT_0348*IT_0373 + -IT_0237*IT_0374 + -IT_0241
      *IT_0375);
    const ccomplex_t IT_0377 = 2*conj(IT_0262)*IT_0332 + (-0.5)*IT_0259
      *IT_0338 + -IT_0295*(conj(IT_0252) + 0.5*IT_0339) + conj(IT_0333)*IT_0341 
      + (-4)*IT_0335*conj(IT_0342) + (-0.5)*IT_0315*IT_0362 + IT_0237*IT_0374 +
       IT_0241*IT_0375;
    const ccomplex_t IT_0378 = (-4)*IT_0333*(IT_0332*conj(IT_0333) + (-2)*conj
      (IT_0260)*IT_0335 + 1./2*conj(IT_0262)*IT_0341 + -1./2*IT_0315*conj
      (IT_0342) + 1./2*IT_0348*IT_0370 + (-0.25)*IT_0355*IT_0373 + 1./2*IT_0239
      *IT_0374 + 1./2*IT_0340*IT_0375);
    const ccomplex_t IT_0379 = IT_0238*IT_0309;
    const ccomplex_t IT_0380 = IT_0250*IT_0313;
    const ccomplex_t IT_0381 = (conj(IT_0105)*IT_0111 + conj(IT_0124)*IT_0128 
      + conj(IT_0224)*IT_0225 + conj(IT_0234)*IT_0235)*IT_0236 + IT_0224*IT_0225
      *IT_0237 + IT_0111*conj(IT_0224)*IT_0238 + IT_0128*(IT_0124*IT_0237 + conj
      (IT_0234)*IT_0238 + IT_0234*IT_0239) + IT_0235*(IT_0234*IT_0237 + conj
      (IT_0224)*IT_0240 + IT_0224*IT_0241) + (conj(IT_0234)*IT_0240 + IT_0234
      *IT_0241)*IT_0247 + IT_0249*(conj(IT_0124)*IT_0240 + conj(IT_0234)*IT_0250
      ) + conj(IT_0252)*IT_0259*IT_0261 + conj(IT_0262)*(IT_0259*IT_0289 +
       IT_0295*IT_0305 + IT_0236*IT_0309 + IT_0240*IT_0313 + IT_0261*IT_0315) +
       4*IT_0322*(IT_0304*IT_0323 + IT_0324*(IT_0287 + 0.25*IT_0325)) + 6
      *IT_0029*(IT_0304*IT_0324 + IT_0326*(conj(IT_0252) + 0.166666666666667
      *IT_0327)) + IT_0344 + conj(IT_0342)*(IT_0236*IT_0309 + IT_0240*IT_0313 +
       IT_0348*IT_0351 + IT_0355*IT_0356) + IT_0357*(IT_0238*IT_0324 + IT_0239
      *IT_0326 + IT_0323*IT_0358 + IT_0304*IT_0359) + 2*IT_0252*IT_0363 +
       IT_0361*(IT_0250*conj(IT_0252) + conj(IT_0234)*IT_0305 + IT_0234*IT_0339 
      + IT_0324*IT_0364 + IT_0326*IT_0365) + IT_0367 + IT_0376 + (-2)*IT_0262
      *IT_0377 + IT_0378 + conj(IT_0260)*(IT_0289*IT_0295 + IT_0259*IT_0305 +
       IT_0379 + IT_0380) + conj(IT_0333)*(IT_0351*IT_0355 + IT_0348*IT_0356 +
       IT_0379 + IT_0380);
    return create_ccomplex_return(IT_0381);
}

