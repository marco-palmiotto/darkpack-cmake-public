#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_A0_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_Hp = param->Gamma_Hp;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = cos(beta);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0025 = IT_0015*IT_0024;
    const ccomplex_t IT_0026 = (-0.5)*IT_0005 + 0.5*IT_0009 + 0.5*IT_0014 + (
      -0.5)*IT_0017 + (-0.5)*IT_0019 + 0.5*IT_0021 + 0.5*IT_0023 + (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0028 = IT_0011*IT_0027;
    const ccomplex_t IT_0029 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0030 = IT_0011*IT_0029;
    const ccomplex_t IT_0031 = 1.4142135623731*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*(IT_0028 + (-0.5)*IT_0031);
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = pow(m_A0, 2);
    const ccomplex_t IT_0035 = pow(m_N_3, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0034 + IT_0035 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0026*IT_0038;
    const ccomplex_t IT_0040 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0041 = IT_0003*IT_0040;
    const ccomplex_t IT_0042 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0045 = IT_0012*IT_0044;
    const ccomplex_t IT_0046 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0047 = IT_0015*IT_0046;
    const ccomplex_t IT_0048 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0053 = IT_0012*IT_0052;
    const ccomplex_t IT_0054 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0055 = IT_0015*IT_0054;
    const ccomplex_t IT_0056 = (-0.5)*IT_0041 + 0.5*IT_0043 + 0.5*IT_0045 + (
      -0.5)*IT_0047 + (-0.5)*IT_0049 + 0.5*IT_0051 + 0.5*IT_0053 + (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0058 = IT_0011*IT_0057;
    const ccomplex_t IT_0059 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0060 = IT_0011*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-0.5)*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0034 + IT_0035 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = IT_0056*IT_0066;
    const ccomplex_t IT_0068 = e_em*IT_0011;
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0071 = IT_0012*IT_0070;
    const ccomplex_t IT_0072 = 2*IT_0071;
    const ccomplex_t IT_0073 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0074 = IT_0003*IT_0073;
    const ccomplex_t IT_0075 = 1.4142135623731*IT_0074;
    const ccomplex_t IT_0076 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0077 = IT_0012*IT_0076;
    const ccomplex_t IT_0078 = 1.4142135623731*IT_0077;
    const ccomplex_t IT_0079 = IT_0075 + IT_0078;
    const ccomplex_t IT_0080 = IT_0072 + IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = 0.5*IT_0081;
    const ccomplex_t IT_0083 = pow(m_C_2, 2);
    const ccomplex_t IT_0084 = cpow((-2)*s_12 + -IT_0035 + -IT_0083 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0085 = IT_0082*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = IT_0069*IT_0086;
    const ccomplex_t IT_0088 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0089 = IT_0003*IT_0088;
    const ccomplex_t IT_0090 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0091 = IT_0007*IT_0090;
    const ccomplex_t IT_0092 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0093 = IT_0012*IT_0092;
    const ccomplex_t IT_0094 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0095 = IT_0015*IT_0094;
    const ccomplex_t IT_0096 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0097 = IT_0003*IT_0096;
    const ccomplex_t IT_0098 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0099 = IT_0007*IT_0098;
    const ccomplex_t IT_0100 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0101 = IT_0012*IT_0100;
    const ccomplex_t IT_0102 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0103 = IT_0015*IT_0102;
    const ccomplex_t IT_0104 = (-0.5)*IT_0089 + 0.5*IT_0091 + 0.5*IT_0093 + (
      -0.5)*IT_0095 + (-0.5)*IT_0097 + 0.5*IT_0099 + 0.5*IT_0101 + (-0.5)*IT_0103;
    const ccomplex_t IT_0105 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0106 = IT_0011*IT_0105;
    const ccomplex_t IT_0107 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0108 = IT_0011*IT_0107;
    const ccomplex_t IT_0109 = 1.4142135623731*IT_0108;
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*(IT_0106 + (-0.5)*IT_0109);
    const ccomplex_t IT_0111 = -IT_0110;
    const ccomplex_t IT_0112 = cpow((-2)*s_23 + IT_0034 + IT_0035 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0111*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0104*IT_0114;
    const ccomplex_t IT_0116 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0117 = IT_0003*IT_0116;
    const ccomplex_t IT_0118 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0119 = IT_0007*IT_0118;
    const ccomplex_t IT_0120 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0121 = IT_0012*IT_0120;
    const ccomplex_t IT_0122 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0123 = IT_0015*IT_0122;
    const ccomplex_t IT_0124 = -IT_0117 + IT_0119 + IT_0121 + -IT_0123;
    const ccomplex_t IT_0125 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0126 = IT_0011*IT_0125;
    const ccomplex_t IT_0127 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0128 = IT_0011*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + (-0.5)*IT_0129);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = cpow((-2)*s_23 + IT_0034 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0124*IT_0134;
    const ccomplex_t IT_0136 = (-2)*IT_0039 + (-2)*IT_0067 + IT_0087 + (-2)
      *IT_0115 + (-2)*IT_0135;
    const ccomplex_t IT_0137 = pow(m_W, -2);
    const ccomplex_t IT_0138 = pow(s_14, 2);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = -IT_0083;
    const ccomplex_t IT_0141 = IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = s_12*IT_0141;
    const ccomplex_t IT_0143 = cpow((-2)*s_13 + IT_0034 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0144 = V_u2*U_Wm2;
    const ccomplex_t IT_0145 = IT_0000*IT_0144;
    const ccomplex_t IT_0146 = U_d2*V_Wp2;
    const ccomplex_t IT_0147 = IT_0006*IT_0146;
    const ccomplex_t IT_0148 = IT_0011*(IT_0145 + IT_0147);
    const ccomplex_t IT_0149 = 1.4142135623731*e_em*IT_0148;
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0152 = IT_0011*IT_0151;
    const ccomplex_t IT_0153 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0154 = IT_0011*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*(IT_0152 + 0.5*IT_0155);
    const ccomplex_t IT_0157 = -IT_0156;
    const ccomplex_t IT_0158 = IT_0150*IT_0157;
    const ccomplex_t IT_0159 = IT_0143*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = 2*IT_0160;
    const ccomplex_t IT_0162 = s_14*s_24*IT_0137;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = s_12 + IT_0163;
    const ccomplex_t IT_0165 = s_12*IT_0164;
    const ccomplex_t IT_0166 = pow(s_24, 2);
    const ccomplex_t IT_0167 = IT_0137*IT_0166;
    const ccomplex_t IT_0168 = -IT_0035;
    const ccomplex_t IT_0169 = IT_0167 + IT_0168;
    const ccomplex_t IT_0170 = s_12*IT_0169;
    const ccomplex_t IT_0171 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0172 = IT_0003*IT_0171;
    const ccomplex_t IT_0173 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0174 = IT_0007*IT_0173;
    const ccomplex_t IT_0175 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0176 = IT_0012*IT_0175;
    const ccomplex_t IT_0177 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0178 = IT_0015*IT_0177;
    const ccomplex_t IT_0179 = IT_0172 + -IT_0174 + -IT_0176 + IT_0178;
    const ccomplex_t IT_0180 = IT_0132*IT_0157;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = IT_0179*IT_0181;
    const ccomplex_t IT_0183 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0184 = IT_0003*IT_0183;
    const ccomplex_t IT_0185 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0186 = IT_0007*IT_0185;
    const ccomplex_t IT_0187 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0188 = IT_0012*IT_0187;
    const ccomplex_t IT_0189 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0190 = IT_0015*IT_0189;
    const ccomplex_t IT_0191 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0192 = IT_0003*IT_0191;
    const ccomplex_t IT_0193 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0194 = IT_0007*IT_0193;
    const ccomplex_t IT_0195 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0196 = IT_0012*IT_0195;
    const ccomplex_t IT_0197 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0198 = IT_0015*IT_0197;
    const ccomplex_t IT_0199 = 0.5*IT_0184 + (-0.5)*IT_0186 + (-0.5)*IT_0188 +
       0.5*IT_0190 + 0.5*IT_0192 + (-0.5)*IT_0194 + (-0.5)*IT_0196 + 0.5*IT_0198;
    const ccomplex_t IT_0200 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0201 = IT_0011*IT_0200;
    const ccomplex_t IT_0202 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0203 = IT_0011*IT_0202;
    const ccomplex_t IT_0204 = 1.4142135623731*IT_0203;
    const ccomplex_t IT_0205 = (0 + _Complex_I*1)*(IT_0201 + 0.5*IT_0204);
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = IT_0112*IT_0206;
    const ccomplex_t IT_0208 = (0 + _Complex_I*1)*IT_0207;
    const ccomplex_t IT_0209 = IT_0199*IT_0208;
    const ccomplex_t IT_0210 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0211 = IT_0003*IT_0210;
    const ccomplex_t IT_0212 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0213 = IT_0007*IT_0212;
    const ccomplex_t IT_0214 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0215 = IT_0012*IT_0214;
    const ccomplex_t IT_0216 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0217 = IT_0015*IT_0216;
    const ccomplex_t IT_0218 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0219 = IT_0003*IT_0218;
    const ccomplex_t IT_0220 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0221 = IT_0007*IT_0220;
    const ccomplex_t IT_0222 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0223 = IT_0012*IT_0222;
    const ccomplex_t IT_0224 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0225 = IT_0015*IT_0224;
    const ccomplex_t IT_0226 = 0.5*IT_0211 + (-0.5)*IT_0213 + (-0.5)*IT_0215 +
       0.5*IT_0217 + 0.5*IT_0219 + (-0.5)*IT_0221 + (-0.5)*IT_0223 + 0.5*IT_0225;
    const ccomplex_t IT_0227 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0228 = IT_0011*IT_0227;
    const ccomplex_t IT_0229 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0230 = IT_0011*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + 0.5*IT_0231);
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0036*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*IT_0234;
    const ccomplex_t IT_0236 = IT_0226*IT_0235;
    const ccomplex_t IT_0237 = 2*IT_0182 + 2*IT_0209 + 2*IT_0236;
    const ccomplex_t IT_0238 = m_C_2*m_N_3;
    const ccomplex_t IT_0239 = IT_0141*IT_0238;
    const ccomplex_t IT_0240 = IT_0164*IT_0238;
    const ccomplex_t IT_0241 = V_u2*U_Wm1;
    const ccomplex_t IT_0242 = IT_0000*IT_0241;
    const ccomplex_t IT_0243 = U_d1*V_Wp2;
    const ccomplex_t IT_0244 = IT_0006*IT_0243;
    const ccomplex_t IT_0245 = IT_0011*(IT_0242 + IT_0244);
    const ccomplex_t IT_0246 = 1.4142135623731*e_em*IT_0245;
    const ccomplex_t IT_0247 = 0.5*IT_0246;
    const ccomplex_t IT_0248 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0249 = IT_0011*IT_0248;
    const ccomplex_t IT_0250 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0251 = IT_0011*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0249 + 0.5*IT_0252);
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = IT_0247*IT_0254;
    const ccomplex_t IT_0256 = cpow((-2)*s_13 + IT_0034 + IT_0083 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0257 = IT_0255*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = IT_0087 + (-2)*IT_0258;
    const ccomplex_t IT_0260 = (-2)*IT_0259;
    const ccomplex_t IT_0261 = 0.5*IT_0260;
    const ccomplex_t IT_0262 = IT_0161 + IT_0261;
    const ccomplex_t IT_0263 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0264 = IT_0000*IT_0263;
    const ccomplex_t IT_0265 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0266 = IT_0006*IT_0265;
    const ccomplex_t IT_0267 = IT_0011*(IT_0264 + IT_0266);
    const ccomplex_t IT_0268 = 1.4142135623731*e_em*IT_0267;
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0271 = IT_0011*IT_0270;
    const ccomplex_t IT_0272 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0273 = IT_0011*IT_0272;
    const ccomplex_t IT_0274 = 1.4142135623731*IT_0273;
    const ccomplex_t IT_0275 = (0 + _Complex_I*1)*(IT_0271 + (-0.5)*IT_0274);
    const ccomplex_t IT_0276 = -IT_0275;
    const ccomplex_t IT_0277 = IT_0269*IT_0276;
    const ccomplex_t IT_0278 = IT_0256*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*IT_0278;
    const ccomplex_t IT_0280 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0281 = IT_0015*IT_0280;
    const ccomplex_t IT_0282 = (-2)*IT_0281;
    const ccomplex_t IT_0283 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0284 = IT_0007*IT_0283;
    const ccomplex_t IT_0285 = 1.4142135623731*IT_0284;
    const ccomplex_t IT_0286 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0287 = IT_0015*IT_0286;
    const ccomplex_t IT_0288 = 1.4142135623731*IT_0287;
    const ccomplex_t IT_0289 = IT_0285 + IT_0288;
    const ccomplex_t IT_0290 = IT_0282 + IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*IT_0290;
    const ccomplex_t IT_0292 = 0.5*IT_0291;
    const ccomplex_t IT_0293 = IT_0084*IT_0292;
    const ccomplex_t IT_0294 = (0 + _Complex_I*1)*IT_0293;
    const ccomplex_t IT_0295 = IT_0069*IT_0294;
    const ccomplex_t IT_0296 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0297 = IT_0000*IT_0296;
    const ccomplex_t IT_0298 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0299 = IT_0006*IT_0298;
    const ccomplex_t IT_0300 = IT_0011*(IT_0297 + IT_0299);
    const ccomplex_t IT_0301 = 1.4142135623731*e_em*IT_0300;
    const ccomplex_t IT_0302 = (-0.5)*IT_0301;
    const ccomplex_t IT_0303 = IT_0131*IT_0302;
    const ccomplex_t IT_0304 = IT_0143*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*IT_0304;
    const ccomplex_t IT_0306 = (-2)*IT_0279 + -IT_0295 + (-2)*IT_0305;
    const ccomplex_t IT_0307 = IT_0169*IT_0238;
    const ccomplex_t IT_0308 = (-2)*conj(IT_0259);
    const ccomplex_t IT_0309 = 0.5*IT_0308;
    const ccomplex_t IT_0310 = conj(IT_0161) + IT_0309;
    const ccomplex_t IT_0311 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0312 = IT_0003*IT_0311;
    const ccomplex_t IT_0313 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0314 = IT_0007*IT_0313;
    const ccomplex_t IT_0315 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0316 = IT_0012*IT_0315;
    const ccomplex_t IT_0317 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0318 = IT_0015*IT_0317;
    const ccomplex_t IT_0319 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0320 = IT_0003*IT_0319;
    const ccomplex_t IT_0321 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0322 = IT_0007*IT_0321;
    const ccomplex_t IT_0323 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0324 = IT_0012*IT_0323;
    const ccomplex_t IT_0325 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0326 = IT_0015*IT_0325;
    const ccomplex_t IT_0327 = 0.5*IT_0312 + (-0.5)*IT_0314 + (-0.5)*IT_0316 +
       0.5*IT_0318 + 0.5*IT_0320 + (-0.5)*IT_0322 + (-0.5)*IT_0324 + 0.5*IT_0326;
    const ccomplex_t IT_0328 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0329 = IT_0011*IT_0328;
    const ccomplex_t IT_0330 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0331 = IT_0011*IT_0330;
    const ccomplex_t IT_0332 = 1.4142135623731*IT_0331;
    const ccomplex_t IT_0333 = (0 + _Complex_I*1)*(IT_0329 + 0.5*IT_0332);
    const ccomplex_t IT_0334 = -IT_0333;
    const ccomplex_t IT_0335 = IT_0064*IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1)*IT_0335;
    const ccomplex_t IT_0337 = IT_0327*IT_0336;
    const ccomplex_t IT_0338 = -IT_0295 + (-2)*IT_0337;
    const ccomplex_t IT_0339 = (-2)*conj(IT_0338);
    const ccomplex_t IT_0340 = 0.5*IT_0339;
    const ccomplex_t IT_0341 = conj(IT_0237) + IT_0340;
    const ccomplex_t IT_0342 = m_N_3*IT_0164;
    const ccomplex_t IT_0343 = IT_0157*IT_0302;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0345 = IT_0143*IT_0343*IT_0344;
    const ccomplex_t IT_0346 = m_N_1*IT_0026;
    const ccomplex_t IT_0347 = IT_0235*IT_0346;
    const ccomplex_t IT_0348 = m_N_3*IT_0124;
    const ccomplex_t IT_0349 = IT_0181*IT_0348;
    const ccomplex_t IT_0350 = IT_0254*IT_0269;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0352 = IT_0256*IT_0350*IT_0351;
    const ccomplex_t IT_0353 = m_N_2*IT_0104;
    const ccomplex_t IT_0354 = IT_0208*IT_0353;
    const ccomplex_t IT_0355 = IT_0345 + -IT_0347 + -IT_0349 + IT_0352 + 
      -IT_0354;
    const ccomplex_t IT_0356 = m_C_2*IT_0026;
    const ccomplex_t IT_0357 = IT_0038*IT_0356;
    const ccomplex_t IT_0358 = m_C_2*IT_0104;
    const ccomplex_t IT_0359 = IT_0114*IT_0358;
    const ccomplex_t IT_0360 = m_N_3*IT_0279;
    const ccomplex_t IT_0361 = m_C_2*IT_0124;
    const ccomplex_t IT_0362 = IT_0134*IT_0361;
    const ccomplex_t IT_0363 = m_N_3*IT_0305;
    const ccomplex_t IT_0364 = m_N_4*IT_0056;
    const ccomplex_t IT_0365 = IT_0336*IT_0364;
    const ccomplex_t IT_0366 = m_C_2*IT_0056;
    const ccomplex_t IT_0367 = IT_0066*IT_0366;
    const ccomplex_t IT_0368 = IT_0357 + IT_0359 + -IT_0360 + IT_0362 + 
      -IT_0363 + -IT_0365 + IT_0367;
    const ccomplex_t IT_0369 = conj(IT_0355) + conj(IT_0368);
    const ccomplex_t IT_0370 = m_C_2*IT_0169;
    const ccomplex_t IT_0371 = m_C_2*IT_0226;
    const ccomplex_t IT_0372 = IT_0235*IT_0371;
    const ccomplex_t IT_0373 = m_C_2*IT_0199;
    const ccomplex_t IT_0374 = IT_0208*IT_0373;
    const ccomplex_t IT_0375 = m_C_2*IT_0179;
    const ccomplex_t IT_0376 = IT_0181*IT_0375;
    const ccomplex_t IT_0377 = m_N_3*IT_0160;
    const ccomplex_t IT_0378 = IT_0372 + IT_0374 + IT_0376 + -IT_0377;
    const ccomplex_t IT_0379 = IT_0247*IT_0276;
    const ccomplex_t IT_0380 = IT_0256*IT_0351*IT_0379;
    const ccomplex_t IT_0381 = IT_0131*IT_0150;
    const ccomplex_t IT_0382 = IT_0143*IT_0344*IT_0381;
    const ccomplex_t IT_0383 = m_N_1*IT_0226;
    const ccomplex_t IT_0384 = IT_0038*IT_0383;
    const ccomplex_t IT_0385 = m_N_2*IT_0199;
    const ccomplex_t IT_0386 = IT_0114*IT_0385;
    const ccomplex_t IT_0387 = m_N_3*IT_0258;
    const ccomplex_t IT_0388 = m_N_4*IT_0327;
    const ccomplex_t IT_0389 = IT_0066*IT_0388;
    const ccomplex_t IT_0390 = m_N_3*IT_0179;
    const ccomplex_t IT_0391 = IT_0134*IT_0390;
    const ccomplex_t IT_0392 = m_C_2*IT_0327;
    const ccomplex_t IT_0393 = IT_0336*IT_0392;
    const ccomplex_t IT_0394 = IT_0380 + IT_0382 + -IT_0384 + -IT_0386 + 
      -IT_0387 + -IT_0389 + -IT_0391 + IT_0393;
    const ccomplex_t IT_0395 = conj(IT_0378) + conj(IT_0394);
    const ccomplex_t IT_0396 = pow(m_W, 2);
    const ccomplex_t IT_0397 = s_24*IT_0396;
    const ccomplex_t IT_0398 = IT_0137*IT_0397;
    const ccomplex_t IT_0399 = 0.5*IT_0398;
    const ccomplex_t IT_0400 = s_24 + IT_0399;
    const ccomplex_t IT_0401 = m_C_2*IT_0400;
    const ccomplex_t IT_0402 = 2*IT_0401;
    const ccomplex_t IT_0403 = IT_0279 + IT_0305;
    const ccomplex_t IT_0404 = IT_0378 + IT_0394;
    const ccomplex_t IT_0405 = conj(IT_0279) + conj(IT_0305);
    const ccomplex_t IT_0406 = -IT_0160;
    const ccomplex_t IT_0407 = (-2)*IT_0355;
    const ccomplex_t IT_0408 = (-2)*IT_0368;
    const ccomplex_t IT_0409 = IT_0407 + IT_0408;
    const ccomplex_t IT_0410 = (-2)*conj(IT_0355);
    const ccomplex_t IT_0411 = (-2)*conj(IT_0368);
    const ccomplex_t IT_0412 = IT_0410 + IT_0411;
    const ccomplex_t IT_0413 = (-0.25)*IT_0398;
    const ccomplex_t IT_0414 = s_24 + IT_0413;
    const ccomplex_t IT_0415 = m_C_2*IT_0414;
    const ccomplex_t IT_0416 = (-4)*IT_0415;
    const ccomplex_t IT_0417 = IT_0039 + IT_0067 + IT_0115 + IT_0135;
    const ccomplex_t IT_0418 = -IT_0182 + -IT_0209 + -IT_0236;
    const ccomplex_t IT_0419 = (-2)*IT_0378;
    const ccomplex_t IT_0420 = (-2)*IT_0394;
    const ccomplex_t IT_0421 = IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = (-2)*conj(IT_0378);
    const ccomplex_t IT_0423 = (-2)*conj(IT_0394);
    const ccomplex_t IT_0424 = IT_0422 + IT_0423;
    const ccomplex_t IT_0425 = s_12*IT_0396;
    const ccomplex_t IT_0426 = s_14*s_24;
    const ccomplex_t IT_0427 = (-2)*IT_0426;
    const ccomplex_t IT_0428 = IT_0425 + IT_0427;
    const ccomplex_t IT_0429 = IT_0137*IT_0428;
    const ccomplex_t IT_0430 = (-0.5)*IT_0429;
    const ccomplex_t IT_0431 = s_12 + IT_0430;
    const ccomplex_t IT_0432 = 4*IT_0378;
    const ccomplex_t IT_0433 = pow(m_W, 4);
    const ccomplex_t IT_0434 = s_12*IT_0433;
    const ccomplex_t IT_0435 = IT_0137*IT_0434;
    const ccomplex_t IT_0436 = (-0.25)*IT_0435;
    const ccomplex_t IT_0437 = IT_0426 + IT_0436;
    const ccomplex_t IT_0438 = -IT_0258;
    const ccomplex_t IT_0439 = IT_0406 + IT_0438;
    const ccomplex_t IT_0440 = -conj(IT_0258);
    const ccomplex_t IT_0441 = conj(IT_0406) + IT_0440;
    const ccomplex_t IT_0442 = 8*IT_0418;
    const ccomplex_t IT_0443 = (-8)*IT_0337;
    const ccomplex_t IT_0444 = 8*conj(IT_0418);
    const ccomplex_t IT_0445 = (-8)*conj(IT_0337);
    const ccomplex_t IT_0446 = (-0.5)*IT_0435;
    const ccomplex_t IT_0447 = IT_0425 + IT_0427 + IT_0446;
    const ccomplex_t IT_0448 = -IT_0337;
    const ccomplex_t IT_0449 = IT_0418 + IT_0448;
    const ccomplex_t IT_0450 = -conj(IT_0337);
    const ccomplex_t IT_0451 = conj(IT_0418) + IT_0450;
    const ccomplex_t IT_0452 = (-4)*IT_0305;
    const ccomplex_t IT_0453 = m_C_2*IT_0164;
    const ccomplex_t IT_0454 = m_N_3*IT_0141;
    const ccomplex_t IT_0455 = (-2)*IT_0136;
    const ccomplex_t IT_0456 = (-2)*IT_0338;
    const ccomplex_t IT_0457 = s_12*s_24;
    const ccomplex_t IT_0458 = s_14*IT_0035;
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = IT_0457 + IT_0459;
    const ccomplex_t IT_0461 = (-2)*conj(IT_0161);
    const ccomplex_t IT_0462 = 0.5*IT_0461;
    const ccomplex_t IT_0463 = conj(IT_0259) + IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0161;
    const ccomplex_t IT_0465 = IT_0406 + IT_0417;
    const ccomplex_t IT_0466 = IT_0262*(conj(IT_0406) + conj(IT_0417)) + conj
      (IT_0306)*(IT_0279 + IT_0305 + IT_0449) + IT_0306*(conj(IT_0279) + conj
      (IT_0305) + IT_0451) + IT_0258*IT_0463 + conj(IT_0258)*(IT_0259 + 0.5
      *IT_0464) + (conj(IT_0161) + 0.5*IT_0308)*IT_0465;
    const ccomplex_t IT_0467 = (-2)*IT_0237;
    const ccomplex_t IT_0468 = IT_0238*IT_0396;
    const ccomplex_t IT_0469 = 3*IT_0468;
    const ccomplex_t IT_0470 = (-2)*IT_0417;
    const ccomplex_t IT_0471 = (-2)*conj(IT_0417);
    const ccomplex_t IT_0472 = 6*conj(IT_0355);
    const ccomplex_t IT_0473 = 0.166666666666667*IT_0472;
    const ccomplex_t IT_0474 = conj(IT_0368) + IT_0473;
    const ccomplex_t IT_0475 = (-2)*conj(IT_0136);
    const ccomplex_t IT_0476 = 0.333333333333333*conj(IT_0306);
    const ccomplex_t IT_0477 = s_14*IT_0396;
    const ccomplex_t IT_0478 = IT_0137*IT_0477;
    const ccomplex_t IT_0479 = 0.5*IT_0478;
    const ccomplex_t IT_0480 = s_14 + IT_0479;
    const ccomplex_t IT_0481 = m_N_3*IT_0480;
    const ccomplex_t IT_0482 = 2*IT_0481;
    const ccomplex_t IT_0483 = 0.5*IT_0337;
    const ccomplex_t IT_0484 = 0.5*conj(IT_0337);
    const ccomplex_t IT_0485 = (-3)*IT_0468;
    const ccomplex_t IT_0486 = (-2)*IT_0258;
    const ccomplex_t IT_0487 = (-2)*conj(IT_0258);
    const ccomplex_t IT_0488 = (-0.25)*IT_0478;
    const ccomplex_t IT_0489 = s_14 + IT_0488;
    const ccomplex_t IT_0490 = m_N_3*IT_0489;
    const ccomplex_t IT_0491 = (-4)*IT_0490;
    const ccomplex_t IT_0492 = 0.5*IT_0258;
    const ccomplex_t IT_0493 = 0.5*conj(IT_0258);
    const ccomplex_t IT_0494 = s_12*s_14;
    const ccomplex_t IT_0495 = s_24*IT_0083;
    const ccomplex_t IT_0496 = -IT_0495;
    const ccomplex_t IT_0497 = IT_0494 + IT_0496;
    const ccomplex_t IT_0498 = (-2)*conj(IT_0237);
    const ccomplex_t IT_0499 = 0.5*IT_0498;
    const ccomplex_t IT_0500 = conj(IT_0338) + IT_0499;
    const ccomplex_t IT_0501 = IT_0136*(conj(IT_0406) + conj(IT_0417)) +
       IT_0341*(IT_0279 + IT_0305 + IT_0418) + (conj(IT_0279) + conj(IT_0305) +
       conj(IT_0418))*(IT_0237 + 0.5*IT_0456) + conj(IT_0136)*IT_0465 + conj
      (IT_0337)*(IT_0338 + 0.5*IT_0467) + IT_0475*IT_0492 + IT_0455*IT_0493 +
       IT_0337*IT_0500;
    const ccomplex_t IT_0502 = 0.5*IT_0454;
    const ccomplex_t IT_0503 = IT_0342*IT_0476;
    const ccomplex_t IT_0504 = 0.333333333333333*conj(IT_0418)*IT_0482;
    const ccomplex_t IT_0505 = 0.333333333333333*IT_0416*conj(IT_0417);
    const ccomplex_t IT_0506 = 0.333333333333333*IT_0370*IT_0463;
    const ccomplex_t IT_0507 = 0.666666666666667*IT_0369*IT_0431;
    const ccomplex_t IT_0508 = IT_0395*IT_0453;
    const ccomplex_t IT_0509 = IT_0240*conj(IT_0338);
    const ccomplex_t IT_0510 = IT_0342*IT_0395;
    const ccomplex_t IT_0511 = IT_0369*IT_0370;
    const ccomplex_t IT_0512 = 2*IT_0136*(conj(IT_0136)*IT_0142 + conj(IT_0161
      )*IT_0165) + 2*IT_0161*(conj(IT_0136)*IT_0165 + conj(IT_0161)*IT_0170) + 2
      *conj(IT_0237)*(IT_0136*IT_0239 + IT_0240*IT_0262) + 2*IT_0306*(conj
      (IT_0136)*IT_0240 + IT_0170*conj(IT_0306) + IT_0307*IT_0310 + IT_0165
      *IT_0341 + IT_0342*IT_0369 + IT_0370*IT_0395) + 2*IT_0402*(conj(IT_0258)*
      (IT_0355 + IT_0368) + IT_0258*IT_0369 + IT_0395*IT_0403 + IT_0404*IT_0405 
      + 0.5*conj(IT_0406)*IT_0409 + 0.5*IT_0406*IT_0412) + 2*IT_0416*(IT_0337
      *IT_0395 + IT_0369*IT_0417 + 0.5*conj(IT_0418)*IT_0421 + 0.5*IT_0418
      *IT_0424) + IT_0395*IT_0431*IT_0432 + 8*IT_0437*(conj(IT_0417)*IT_0439 +
       IT_0417*IT_0441 + 0.125*IT_0405*(IT_0442 + IT_0443) + 0.125*IT_0403*
      (IT_0444 + IT_0445)) + (-4)*IT_0447*(IT_0417*conj(IT_0417) + IT_0439
      *IT_0441 + IT_0449*IT_0451 + IT_0405*(IT_0279 + (-0.25)*IT_0452)) + 
      (IT_0165*conj(IT_0259) + IT_0239*conj(IT_0338) + IT_0369*IT_0453 + IT_0395
      *IT_0454)*IT_0455 + 2*conj(IT_0306)*(IT_0136*IT_0240 + IT_0262*IT_0307 +
       0.5*IT_0165*IT_0456) + 2*IT_0460*IT_0466 + conj(IT_0338)*(IT_0407*IT_0454
       + IT_0142*IT_0467) + (-2)*IT_0469*(IT_0441*IT_0449 + IT_0439*IT_0451 + (
      -0.5)*IT_0405*IT_0470 + (-0.5)*IT_0403*IT_0471) + 2*IT_0394*(conj(IT_0337)
      *IT_0416 + 2*IT_0395*IT_0431 + IT_0341*IT_0453 + IT_0342*IT_0463 + 3
      *IT_0238*IT_0474 + 0.5*IT_0454*IT_0475 + 3*IT_0370*IT_0476) + 2*IT_0482*
      (IT_0369*IT_0418 + 0.5*conj(IT_0417)*IT_0421 + 0.5*IT_0417*IT_0424 +
       IT_0412*IT_0483 + IT_0409*IT_0484) + 2*IT_0485*(conj(IT_0417)*IT_0418 +
       IT_0417*conj(IT_0418) + IT_0471*IT_0483 + IT_0470*IT_0484 + IT_0405*
      (IT_0406 + 0.5*IT_0486) + IT_0403*(conj(IT_0406) + 0.5*IT_0487)) + 2
      *IT_0491*(IT_0395*IT_0406 + IT_0404*conj(IT_0406) + 0.5*IT_0405*IT_0409 +
       0.5*IT_0403*IT_0412 + IT_0424*IT_0492 + IT_0421*IT_0493) + 2*IT_0497
      *IT_0501 + 2*IT_0338*(IT_0240*IT_0463 + 0.5*IT_0239*IT_0475 + IT_0142
      *IT_0500 + IT_0410*IT_0502) + 2*IT_0378*(conj(IT_0306)*IT_0370 + conj
      (IT_0337)*IT_0416 + IT_0341*IT_0453 + IT_0342*IT_0463 + 3*IT_0238*IT_0474 
      + IT_0475*IT_0502) + 6*IT_0355*(IT_0238*IT_0395 + 1./3*conj(IT_0237)
      *IT_0454 + 0.166666666666667*IT_0453*IT_0475 + IT_0503 + IT_0504 + IT_0505
       + IT_0506 + IT_0507) + 6*IT_0368*(IT_0238*IT_0395 + 1./3*IT_0341*IT_0454 
      + 0.166666666666667*IT_0453*IT_0475 + IT_0503 + IT_0504 + IT_0505 +
       IT_0506 + IT_0507) + IT_0456*(conj(IT_0368)*IT_0454 + IT_0508) + 2
      *IT_0237*(IT_0142*conj(IT_0237) + conj(IT_0136)*IT_0239 + IT_0165*conj
      (IT_0306) + IT_0240*IT_0310 + IT_0369*IT_0454 + IT_0508) + IT_0464*
      (IT_0170*conj(IT_0259) + IT_0509 + IT_0510 + IT_0511) + 2*IT_0259*(IT_0170
      *IT_0463 + 0.5*IT_0165*IT_0475 + IT_0509 + IT_0510 + IT_0511);
    return create_ccomplex_return(IT_0512);
}

