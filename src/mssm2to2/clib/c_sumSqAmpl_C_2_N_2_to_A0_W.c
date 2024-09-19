#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_2_to_A0_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_2_to_A0_W(
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
    const ccomplex_t IT_0004 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = sin(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0000*IT_0011;
    const ccomplex_t IT_0013 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = IT_0006*IT_0011;
    const ccomplex_t IT_0016 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0021 = IT_0007*IT_0020;
    const ccomplex_t IT_0022 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0023 = IT_0012*IT_0022;
    const ccomplex_t IT_0024 = N_W1*N_d2*e_em;
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
    const ccomplex_t IT_0035 = pow(m_N_2, 2);
    const ccomplex_t IT_0036 = cpow((-2)*s_23 + IT_0034 + IT_0035 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0037 = IT_0033*IT_0036;
    const ccomplex_t IT_0038 = (0 + _Complex_I*1)*IT_0037;
    const ccomplex_t IT_0039 = IT_0026*IT_0038;
    const ccomplex_t IT_0040 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0041 = IT_0003*IT_0040;
    const ccomplex_t IT_0042 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0043 = IT_0007*IT_0042;
    const ccomplex_t IT_0044 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0045 = IT_0012*IT_0044;
    const ccomplex_t IT_0046 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0047 = IT_0015*IT_0046;
    const ccomplex_t IT_0048 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0053 = IT_0012*IT_0052;
    const ccomplex_t IT_0054 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0055 = IT_0015*IT_0054;
    const ccomplex_t IT_0056 = (-0.5)*IT_0041 + 0.5*IT_0043 + 0.5*IT_0045 + (
      -0.5)*IT_0047 + (-0.5)*IT_0049 + 0.5*IT_0051 + 0.5*IT_0053 + (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0058 = IT_0011*IT_0057;
    const ccomplex_t IT_0059 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0060 = IT_0011*IT_0059;
    const ccomplex_t IT_0061 = 1.4142135623731*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*(IT_0058 + (-0.5)*IT_0061);
    const ccomplex_t IT_0063 = -IT_0062;
    const ccomplex_t IT_0064 = cpow((-2)*s_23 + IT_0034 + IT_0035 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0065 = IT_0063*IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = IT_0056*IT_0066;
    const ccomplex_t IT_0068 = N_B4*N_u2*e_em;
    const ccomplex_t IT_0069 = IT_0003*IT_0068;
    const ccomplex_t IT_0070 = N_B4*N_d2*e_em;
    const ccomplex_t IT_0071 = IT_0007*IT_0070;
    const ccomplex_t IT_0072 = N_W4*N_u2*e_em;
    const ccomplex_t IT_0073 = IT_0012*IT_0072;
    const ccomplex_t IT_0074 = N_W4*N_d2*e_em;
    const ccomplex_t IT_0075 = IT_0015*IT_0074;
    const ccomplex_t IT_0076 = N_B2*N_u4*e_em;
    const ccomplex_t IT_0077 = IT_0003*IT_0076;
    const ccomplex_t IT_0078 = N_B2*N_d4*e_em;
    const ccomplex_t IT_0079 = IT_0007*IT_0078;
    const ccomplex_t IT_0080 = N_W2*N_u4*e_em;
    const ccomplex_t IT_0081 = IT_0012*IT_0080;
    const ccomplex_t IT_0082 = N_W2*N_d4*e_em;
    const ccomplex_t IT_0083 = IT_0015*IT_0082;
    const ccomplex_t IT_0084 = (-0.5)*IT_0069 + 0.5*IT_0071 + 0.5*IT_0073 + (
      -0.5)*IT_0075 + (-0.5)*IT_0077 + 0.5*IT_0079 + 0.5*IT_0081 + (-0.5)*IT_0083;
    const ccomplex_t IT_0085 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0086 = IT_0011*IT_0085;
    const ccomplex_t IT_0087 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0088 = IT_0011*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + (-0.5)*IT_0089);
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0034 + IT_0035 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0091*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = IT_0084*IT_0094;
    const ccomplex_t IT_0096 = e_em*IT_0011;
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0099 = IT_0012*IT_0098;
    const ccomplex_t IT_0100 = 2*IT_0099;
    const ccomplex_t IT_0101 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0102 = IT_0003*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0105 = IT_0012*IT_0104;
    const ccomplex_t IT_0106 = 1.4142135623731*IT_0105;
    const ccomplex_t IT_0107 = IT_0103 + IT_0106;
    const ccomplex_t IT_0108 = IT_0100 + IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = 0.5*IT_0109;
    const ccomplex_t IT_0111 = pow(m_C_2, 2);
    const ccomplex_t IT_0112 = cpow((-2)*s_12 + -IT_0035 + -IT_0111 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0113 = IT_0110*IT_0112;
    const ccomplex_t IT_0114 = (0 + _Complex_I*1)*IT_0113;
    const ccomplex_t IT_0115 = IT_0097*IT_0114;
    const ccomplex_t IT_0116 = N_B2*N_u2*e_em;
    const ccomplex_t IT_0117 = IT_0003*IT_0116;
    const ccomplex_t IT_0118 = N_B2*N_d2*e_em;
    const ccomplex_t IT_0119 = IT_0007*IT_0118;
    const ccomplex_t IT_0120 = N_W2*N_u2*e_em;
    const ccomplex_t IT_0121 = IT_0012*IT_0120;
    const ccomplex_t IT_0122 = N_W2*N_d2*e_em;
    const ccomplex_t IT_0123 = IT_0015*IT_0122;
    const ccomplex_t IT_0124 = -IT_0117 + IT_0119 + IT_0121 + -IT_0123;
    const ccomplex_t IT_0125 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0126 = IT_0011*IT_0125;
    const ccomplex_t IT_0127 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0128 = IT_0011*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + (-0.5)*IT_0129);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = cpow((-2)*s_23 + IT_0034 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0133 = IT_0131*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0124*IT_0134;
    const ccomplex_t IT_0136 = (-2)*IT_0039 + (-2)*IT_0067 + (-2)*IT_0095 +
       IT_0115 + (-2)*IT_0135;
    const ccomplex_t IT_0137 = pow(m_W, -2);
    const ccomplex_t IT_0138 = pow(s_14, 2);
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = -IT_0111;
    const ccomplex_t IT_0141 = IT_0139 + IT_0140;
    const ccomplex_t IT_0142 = s_12*IT_0141;
    const ccomplex_t IT_0143 = V_u2*U_Wm1;
    const ccomplex_t IT_0144 = IT_0000*IT_0143;
    const ccomplex_t IT_0145 = U_d1*V_Wp2;
    const ccomplex_t IT_0146 = IT_0006*IT_0145;
    const ccomplex_t IT_0147 = IT_0011*(IT_0144 + IT_0146);
    const ccomplex_t IT_0148 = 1.4142135623731*e_em*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0151 = IT_0011*IT_0150;
    const ccomplex_t IT_0152 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0153 = IT_0011*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + 0.5*IT_0154);
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0149*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_13 + IT_0034 + IT_0111 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0157*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_13 + IT_0034 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0162 = V_u2*U_Wm2;
    const ccomplex_t IT_0163 = IT_0000*IT_0162;
    const ccomplex_t IT_0164 = U_d2*V_Wp2;
    const ccomplex_t IT_0165 = IT_0006*IT_0164;
    const ccomplex_t IT_0166 = IT_0011*(IT_0163 + IT_0165);
    const ccomplex_t IT_0167 = 1.4142135623731*e_em*IT_0166;
    const ccomplex_t IT_0168 = 0.5*IT_0167;
    const ccomplex_t IT_0169 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0170 = IT_0011*IT_0169;
    const ccomplex_t IT_0171 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0172 = IT_0011*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0170 + 0.5*IT_0173);
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = IT_0168*IT_0175;
    const ccomplex_t IT_0177 = IT_0161*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = 2*IT_0160 + 2*IT_0178;
    const ccomplex_t IT_0180 = s_14*s_24*IT_0137;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = s_12 + IT_0181;
    const ccomplex_t IT_0183 = s_12*IT_0182;
    const ccomplex_t IT_0184 = pow(s_24, 2);
    const ccomplex_t IT_0185 = IT_0137*IT_0184;
    const ccomplex_t IT_0186 = -IT_0035;
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = s_12*IT_0187;
    const ccomplex_t IT_0189 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0190 = IT_0003*IT_0189;
    const ccomplex_t IT_0191 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0192 = IT_0007*IT_0191;
    const ccomplex_t IT_0193 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0194 = IT_0012*IT_0193;
    const ccomplex_t IT_0195 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0196 = IT_0015*IT_0195;
    const ccomplex_t IT_0197 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0198 = IT_0003*IT_0197;
    const ccomplex_t IT_0199 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0200 = IT_0007*IT_0199;
    const ccomplex_t IT_0201 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0202 = IT_0012*IT_0201;
    const ccomplex_t IT_0203 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0204 = IT_0015*IT_0203;
    const ccomplex_t IT_0205 = 0.5*IT_0190 + (-0.5)*IT_0192 + (-0.5)*IT_0194 +
       0.5*IT_0196 + 0.5*IT_0198 + (-0.5)*IT_0200 + (-0.5)*IT_0202 + 0.5*IT_0204;
    const ccomplex_t IT_0206 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0207 = IT_0011*IT_0206;
    const ccomplex_t IT_0208 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0209 = IT_0011*IT_0208;
    const ccomplex_t IT_0210 = 1.4142135623731*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*(IT_0207 + 0.5*IT_0210);
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0036*IT_0212;
    const ccomplex_t IT_0214 = (0 + _Complex_I*1)*IT_0213;
    const ccomplex_t IT_0215 = IT_0205*IT_0214;
    const ccomplex_t IT_0216 = conj(N_B2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0217 = IT_0003*IT_0216;
    const ccomplex_t IT_0218 = conj(N_B2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0219 = IT_0007*IT_0218;
    const ccomplex_t IT_0220 = conj(N_W2)*conj(N_u2)*e_em;
    const ccomplex_t IT_0221 = IT_0012*IT_0220;
    const ccomplex_t IT_0222 = conj(N_W2)*conj(N_d2)*e_em;
    const ccomplex_t IT_0223 = IT_0015*IT_0222;
    const ccomplex_t IT_0224 = IT_0217 + -IT_0219 + -IT_0221 + IT_0223;
    const ccomplex_t IT_0225 = IT_0132*IT_0175;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0224*IT_0226;
    const ccomplex_t IT_0228 = 2*IT_0215 + 2*IT_0227;
    const ccomplex_t IT_0229 = m_C_2*m_N_2;
    const ccomplex_t IT_0230 = IT_0141*IT_0229;
    const ccomplex_t IT_0231 = IT_0182*IT_0229;
    const ccomplex_t IT_0232 = (-2)*conj(IT_0115);
    const ccomplex_t IT_0233 = 0.5*IT_0232;
    const ccomplex_t IT_0234 = conj(IT_0179) + IT_0233;
    const ccomplex_t IT_0235 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0236 = IT_0015*IT_0235;
    const ccomplex_t IT_0237 = (-2)*IT_0236;
    const ccomplex_t IT_0238 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0239 = IT_0007*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0242 = IT_0015*IT_0241;
    const ccomplex_t IT_0243 = 1.4142135623731*IT_0242;
    const ccomplex_t IT_0244 = IT_0240 + IT_0243;
    const ccomplex_t IT_0245 = IT_0237 + IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*IT_0245;
    const ccomplex_t IT_0247 = 0.5*IT_0246;
    const ccomplex_t IT_0248 = IT_0112*IT_0247;
    const ccomplex_t IT_0249 = (0 + _Complex_I*1)*IT_0248;
    const ccomplex_t IT_0250 = IT_0097*IT_0249;
    const ccomplex_t IT_0251 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0252 = IT_0000*IT_0251;
    const ccomplex_t IT_0253 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0254 = IT_0006*IT_0253;
    const ccomplex_t IT_0255 = IT_0011*(IT_0252 + IT_0254);
    const ccomplex_t IT_0256 = 1.4142135623731*e_em*IT_0255;
    const ccomplex_t IT_0257 = (-0.5)*IT_0256;
    const ccomplex_t IT_0258 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0259 = IT_0011*IT_0258;
    const ccomplex_t IT_0260 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0261 = IT_0011*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*(IT_0259 + (-0.5)*IT_0262);
    const ccomplex_t IT_0264 = -IT_0263;
    const ccomplex_t IT_0265 = IT_0257*IT_0264;
    const ccomplex_t IT_0266 = IT_0158*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0269 = IT_0000*IT_0268;
    const ccomplex_t IT_0270 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0271 = IT_0006*IT_0270;
    const ccomplex_t IT_0272 = IT_0011*(IT_0269 + IT_0271);
    const ccomplex_t IT_0273 = 1.4142135623731*e_em*IT_0272;
    const ccomplex_t IT_0274 = (-0.5)*IT_0273;
    const ccomplex_t IT_0275 = IT_0131*IT_0274;
    const ccomplex_t IT_0276 = IT_0161*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = -IT_0250 + (-2)*IT_0267 + (-2)*IT_0277;
    const ccomplex_t IT_0279 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0280 = IT_0003*IT_0279;
    const ccomplex_t IT_0281 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0282 = IT_0007*IT_0281;
    const ccomplex_t IT_0283 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0284 = IT_0012*IT_0283;
    const ccomplex_t IT_0285 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0286 = IT_0015*IT_0285;
    const ccomplex_t IT_0287 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0288 = IT_0003*IT_0287;
    const ccomplex_t IT_0289 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0290 = IT_0007*IT_0289;
    const ccomplex_t IT_0291 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0292 = IT_0012*IT_0291;
    const ccomplex_t IT_0293 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0294 = IT_0015*IT_0293;
    const ccomplex_t IT_0295 = 0.5*IT_0280 + (-0.5)*IT_0282 + (-0.5)*IT_0284 +
       0.5*IT_0286 + 0.5*IT_0288 + (-0.5)*IT_0290 + (-0.5)*IT_0292 + 0.5*IT_0294;
    const ccomplex_t IT_0296 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0297 = IT_0011*IT_0296;
    const ccomplex_t IT_0298 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0299 = IT_0011*IT_0298;
    const ccomplex_t IT_0300 = 1.4142135623731*IT_0299;
    const ccomplex_t IT_0301 = (0 + _Complex_I*1)*(IT_0297 + 0.5*IT_0300);
    const ccomplex_t IT_0302 = -IT_0301;
    const ccomplex_t IT_0303 = IT_0064*IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*IT_0303;
    const ccomplex_t IT_0305 = IT_0295*IT_0304;
    const ccomplex_t IT_0306 = conj(N_B2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0307 = IT_0003*IT_0306;
    const ccomplex_t IT_0308 = conj(N_B2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0309 = IT_0007*IT_0308;
    const ccomplex_t IT_0310 = conj(N_W2)*conj(N_u4)*e_em;
    const ccomplex_t IT_0311 = IT_0012*IT_0310;
    const ccomplex_t IT_0312 = conj(N_W2)*conj(N_d4)*e_em;
    const ccomplex_t IT_0313 = IT_0015*IT_0312;
    const ccomplex_t IT_0314 = conj(N_B4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0315 = IT_0003*IT_0314;
    const ccomplex_t IT_0316 = conj(N_B4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0317 = IT_0007*IT_0316;
    const ccomplex_t IT_0318 = conj(N_W4)*conj(N_u2)*e_em;
    const ccomplex_t IT_0319 = IT_0012*IT_0318;
    const ccomplex_t IT_0320 = conj(N_W4)*conj(N_d2)*e_em;
    const ccomplex_t IT_0321 = IT_0015*IT_0320;
    const ccomplex_t IT_0322 = 0.5*IT_0307 + (-0.5)*IT_0309 + (-0.5)*IT_0311 +
       0.5*IT_0313 + 0.5*IT_0315 + (-0.5)*IT_0317 + (-0.5)*IT_0319 + 0.5*IT_0321;
    const ccomplex_t IT_0323 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0324 = IT_0011*IT_0323;
    const ccomplex_t IT_0325 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0326 = IT_0011*IT_0325;
    const ccomplex_t IT_0327 = 1.4142135623731*IT_0326;
    const ccomplex_t IT_0328 = (0 + _Complex_I*1)*(IT_0324 + 0.5*IT_0327);
    const ccomplex_t IT_0329 = -IT_0328;
    const ccomplex_t IT_0330 = IT_0092*IT_0329;
    const ccomplex_t IT_0331 = (0 + _Complex_I*1)*IT_0330;
    const ccomplex_t IT_0332 = IT_0322*IT_0331;
    const ccomplex_t IT_0333 = -IT_0250 + (-2)*IT_0305 + (-2)*IT_0332;
    const ccomplex_t IT_0334 = (-2)*IT_0333;
    const ccomplex_t IT_0335 = 0.5*IT_0334;
    const ccomplex_t IT_0336 = IT_0228 + IT_0335;
    const ccomplex_t IT_0337 = IT_0187*IT_0229;
    const ccomplex_t IT_0338 = (-2)*IT_0115;
    const ccomplex_t IT_0339 = 0.5*IT_0338;
    const ccomplex_t IT_0340 = IT_0179 + IT_0339;
    const ccomplex_t IT_0341 = m_N_2*IT_0182;
    const ccomplex_t IT_0342 = IT_0175*IT_0274;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0344 = IT_0161*IT_0342*IT_0343;
    const ccomplex_t IT_0345 = IT_0156*IT_0257;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0347 = IT_0158*IT_0345*IT_0346;
    const ccomplex_t IT_0348 = m_N_1*IT_0026;
    const ccomplex_t IT_0349 = IT_0214*IT_0348;
    const ccomplex_t IT_0350 = m_N_2*IT_0124;
    const ccomplex_t IT_0351 = IT_0226*IT_0350;
    const ccomplex_t IT_0352 = IT_0344 + IT_0347 + -IT_0349 + -IT_0351;
    const ccomplex_t IT_0353 = m_N_2*IT_0267;
    const ccomplex_t IT_0354 = m_N_2*IT_0277;
    const ccomplex_t IT_0355 = m_C_2*IT_0026;
    const ccomplex_t IT_0356 = IT_0038*IT_0355;
    const ccomplex_t IT_0357 = m_C_2*IT_0124;
    const ccomplex_t IT_0358 = IT_0134*IT_0357;
    const ccomplex_t IT_0359 = m_N_3*IT_0056;
    const ccomplex_t IT_0360 = IT_0304*IT_0359;
    const ccomplex_t IT_0361 = m_C_2*IT_0056;
    const ccomplex_t IT_0362 = IT_0066*IT_0361;
    const ccomplex_t IT_0363 = m_N_4*IT_0084;
    const ccomplex_t IT_0364 = IT_0331*IT_0363;
    const ccomplex_t IT_0365 = m_C_2*IT_0084;
    const ccomplex_t IT_0366 = IT_0094*IT_0365;
    const ccomplex_t IT_0367 = -IT_0353 + -IT_0354 + IT_0356 + IT_0358 + 
      -IT_0360 + IT_0362 + -IT_0364 + IT_0366;
    const ccomplex_t IT_0368 = IT_0352 + IT_0367;
    const ccomplex_t IT_0369 = m_C_2*IT_0187;
    const ccomplex_t IT_0370 = m_C_2*IT_0205;
    const ccomplex_t IT_0371 = IT_0214*IT_0370;
    const ccomplex_t IT_0372 = m_N_2*IT_0160;
    const ccomplex_t IT_0373 = m_C_2*IT_0224;
    const ccomplex_t IT_0374 = IT_0226*IT_0373;
    const ccomplex_t IT_0375 = m_N_2*IT_0178;
    const ccomplex_t IT_0376 = IT_0371 + -IT_0372 + IT_0374 + -IT_0375;
    const ccomplex_t IT_0377 = m_N_2*IT_0224;
    const ccomplex_t IT_0378 = IT_0134*IT_0377;
    const ccomplex_t IT_0379 = m_C_2*IT_0295;
    const ccomplex_t IT_0380 = IT_0304*IT_0379;
    const ccomplex_t IT_0381 = m_N_3*IT_0295;
    const ccomplex_t IT_0382 = IT_0066*IT_0381;
    const ccomplex_t IT_0383 = m_C_2*IT_0322;
    const ccomplex_t IT_0384 = IT_0331*IT_0383;
    const ccomplex_t IT_0385 = m_N_4*IT_0322;
    const ccomplex_t IT_0386 = IT_0094*IT_0385;
    const ccomplex_t IT_0387 = IT_0131*IT_0168;
    const ccomplex_t IT_0388 = IT_0161*IT_0343*IT_0387;
    const ccomplex_t IT_0389 = IT_0149*IT_0264;
    const ccomplex_t IT_0390 = IT_0158*IT_0346*IT_0389;
    const ccomplex_t IT_0391 = m_N_1*IT_0205;
    const ccomplex_t IT_0392 = IT_0038*IT_0391;
    const ccomplex_t IT_0393 = -IT_0378 + IT_0380 + -IT_0382 + IT_0384 + 
      -IT_0386 + IT_0388 + IT_0390 + -IT_0392;
    const ccomplex_t IT_0394 = IT_0376 + IT_0393;
    const ccomplex_t IT_0395 = pow(m_W, 2);
    const ccomplex_t IT_0396 = s_12*IT_0395;
    const ccomplex_t IT_0397 = pow(m_W, 4);
    const ccomplex_t IT_0398 = s_12*IT_0397;
    const ccomplex_t IT_0399 = IT_0137*IT_0398;
    const ccomplex_t IT_0400 = (-0.5)*IT_0399;
    const ccomplex_t IT_0401 = s_14*s_24;
    const ccomplex_t IT_0402 = (-2)*IT_0401;
    const ccomplex_t IT_0403 = IT_0396 + IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = IT_0039 + IT_0067 + IT_0095 + IT_0135;
    const ccomplex_t IT_0405 = -IT_0215 + -IT_0227;
    const ccomplex_t IT_0406 = IT_0305 + IT_0332;
    const ccomplex_t IT_0407 = -IT_0406;
    const ccomplex_t IT_0408 = IT_0405 + IT_0407;
    const ccomplex_t IT_0409 = -conj(IT_0406);
    const ccomplex_t IT_0410 = conj(IT_0405) + IT_0409;
    const ccomplex_t IT_0411 = s_24*IT_0395;
    const ccomplex_t IT_0412 = IT_0137*IT_0411;
    const ccomplex_t IT_0413 = (-0.25)*IT_0412;
    const ccomplex_t IT_0414 = s_24 + IT_0413;
    const ccomplex_t IT_0415 = m_C_2*IT_0414;
    const ccomplex_t IT_0416 = (-4)*IT_0415;
    const ccomplex_t IT_0417 = (-2)*IT_0376;
    const ccomplex_t IT_0418 = (-2)*IT_0393;
    const ccomplex_t IT_0419 = IT_0417 + IT_0418;
    const ccomplex_t IT_0420 = (-2)*conj(IT_0376);
    const ccomplex_t IT_0421 = (-2)*conj(IT_0393);
    const ccomplex_t IT_0422 = IT_0420 + IT_0421;
    const ccomplex_t IT_0423 = 6*conj(IT_0352);
    const ccomplex_t IT_0424 = s_12*s_24;
    const ccomplex_t IT_0425 = s_14*IT_0035;
    const ccomplex_t IT_0426 = -IT_0425;
    const ccomplex_t IT_0427 = IT_0424 + IT_0426;
    const ccomplex_t IT_0428 = -IT_0160 + -IT_0178;
    const ccomplex_t IT_0429 = conj(IT_0404) + conj(IT_0428);
    const ccomplex_t IT_0430 = IT_0267 + IT_0277;
    const ccomplex_t IT_0431 = IT_0404 + IT_0428;
    const ccomplex_t IT_0432 = m_N_2*IT_0141;
    const ccomplex_t IT_0433 = (-2)*IT_0179;
    const ccomplex_t IT_0434 = 0.5*IT_0433;
    const ccomplex_t IT_0435 = IT_0115 + IT_0434;
    const ccomplex_t IT_0436 = IT_0396 + IT_0402;
    const ccomplex_t IT_0437 = IT_0137*IT_0436;
    const ccomplex_t IT_0438 = (-0.5)*IT_0437;
    const ccomplex_t IT_0439 = s_12 + IT_0438;
    const ccomplex_t IT_0440 = m_C_2*IT_0182;
    const ccomplex_t IT_0441 = (-2)*IT_0136;
    const ccomplex_t IT_0442 = s_14*IT_0395;
    const ccomplex_t IT_0443 = IT_0137*IT_0442;
    const ccomplex_t IT_0444 = 0.5*IT_0443;
    const ccomplex_t IT_0445 = s_14 + IT_0444;
    const ccomplex_t IT_0446 = m_N_2*IT_0445;
    const ccomplex_t IT_0447 = 2*IT_0446;
    const ccomplex_t IT_0448 = 2*conj(IT_0367)*(IT_0278*IT_0341 + 3*IT_0229
      *IT_0394 + IT_0404*IT_0416 + IT_0336*IT_0432 + IT_0369*IT_0435 + 2*IT_0368
      *IT_0439 + 0.5*IT_0440*IT_0441 + IT_0405*IT_0447);
    const ccomplex_t IT_0449 = 0.5*IT_0412;
    const ccomplex_t IT_0450 = s_24 + IT_0449;
    const ccomplex_t IT_0451 = m_C_2*IT_0450;
    const ccomplex_t IT_0452 = 2*IT_0451;
    const ccomplex_t IT_0453 = (-0.25)*IT_0443;
    const ccomplex_t IT_0454 = s_14 + IT_0453;
    const ccomplex_t IT_0455 = m_N_2*IT_0454;
    const ccomplex_t IT_0456 = (-4)*IT_0455;
    const ccomplex_t IT_0457 = 4*IT_0376;
    const ccomplex_t IT_0458 = 3*IT_0229*IT_0368 + IT_0278*IT_0369 + IT_0406
      *IT_0416 + IT_0341*IT_0435 + IT_0336*IT_0440 + 0.5*IT_0432*IT_0441 +
       IT_0430*IT_0452 + IT_0428*IT_0456 + 2*IT_0439*(IT_0393 + 0.25*IT_0457);
    const ccomplex_t IT_0459 = (-2)*conj(IT_0136);
    const ccomplex_t IT_0460 = s_12*s_14;
    const ccomplex_t IT_0461 = s_24*IT_0111;
    const ccomplex_t IT_0462 = -IT_0461;
    const ccomplex_t IT_0463 = IT_0460 + IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0228;
    const ccomplex_t IT_0465 = 0.5*IT_0464;
    const ccomplex_t IT_0466 = IT_0333 + IT_0465;
    const ccomplex_t IT_0467 = (-2)*conj(IT_0228);
    const ccomplex_t IT_0468 = IT_0405 + IT_0430;
    const ccomplex_t IT_0469 = (-2)*conj(IT_0333);
    const ccomplex_t IT_0470 = (-2)*IT_0352;
    const ccomplex_t IT_0471 = (-2)*conj(IT_0352);
    const ccomplex_t IT_0472 = IT_0229*IT_0395;
    const ccomplex_t IT_0473 = (-3)*IT_0472;
    const ccomplex_t IT_0474 = (-0.25)*IT_0399;
    const ccomplex_t IT_0475 = IT_0401 + IT_0474;
    const ccomplex_t IT_0476 = (-2)*IT_0367;
    const ccomplex_t IT_0477 = IT_0470 + IT_0476;
    const ccomplex_t IT_0478 = 3*IT_0472;
    const ccomplex_t IT_0479 = (-2)*IT_0404;
    const ccomplex_t IT_0480 = (-2)*conj(IT_0367);
    const ccomplex_t IT_0481 = IT_0471 + IT_0480;
    const ccomplex_t IT_0482 = (-2)*conj(IT_0404);
    const ccomplex_t IT_0483 = 0.5*IT_0406;
    const ccomplex_t IT_0484 = 0.5*conj(IT_0406);
    const ccomplex_t IT_0485 = 0.5*IT_0183;
    const ccomplex_t IT_0486 = (-0.5)*IT_0452;
    const ccomplex_t IT_0487 = IT_0394*IT_0440;
    const ccomplex_t IT_0488 = (-2)*conj(IT_0179);
    const ccomplex_t IT_0489 = IT_0231*IT_0333;
    const ccomplex_t IT_0490 = IT_0341*IT_0394;
    const ccomplex_t IT_0491 = IT_0368*IT_0369;
    const ccomplex_t IT_0492 = 2*IT_0136*(conj(IT_0136)*IT_0142 + conj(IT_0179
      )*IT_0183) + 2*IT_0179*(conj(IT_0136)*IT_0183 + conj(IT_0179)*IT_0188) + 2
      *IT_0228*(conj(IT_0136)*IT_0230 + IT_0231*IT_0234) + 2*conj(IT_0278)*
      (IT_0136*IT_0231 + IT_0188*IT_0278 + IT_0183*IT_0336 + IT_0337*IT_0340 +
       IT_0341*IT_0368 + IT_0369*IT_0394) + (-4)*IT_0403*(IT_0404*conj(IT_0404) 
      + IT_0408*IT_0410) + 2*IT_0416*(conj(IT_0352)*IT_0404 + IT_0368*conj
      (IT_0404) + IT_0394*conj(IT_0406) + 0.5*conj(IT_0405)*IT_0419 + 0.5
      *IT_0405*IT_0422) + IT_0229*IT_0394*IT_0423 + 2*IT_0427*(IT_0340*IT_0429 +
       conj(IT_0278)*(IT_0408 + IT_0430) + IT_0278*(IT_0410 + conj(IT_0430)) + 
      (conj(IT_0179) + 0.5*IT_0232)*IT_0431) + 2*conj(IT_0352)*(IT_0278*IT_0341 
      + IT_0228*IT_0432 + IT_0369*IT_0435 + 2*IT_0368*IT_0439 + 0.5*IT_0440
      *IT_0441) + IT_0448 + 2*(conj(IT_0376) + conj(IT_0393))*IT_0458 + (IT_0115
      *IT_0183 + IT_0230*IT_0333 + IT_0394*IT_0432 + IT_0368*IT_0440)*IT_0459 +
       2*IT_0463*(IT_0136*IT_0429 + IT_0336*(conj(IT_0405) + conj(IT_0430)) +
       conj(IT_0136)*IT_0431 + conj(IT_0406)*IT_0466 + IT_0406*(conj(IT_0333) +
       0.5*IT_0467) + IT_0468*(conj(IT_0228) + 0.5*IT_0469)) + 2*conj(IT_0333)*
      (IT_0231*IT_0435 + 0.5*IT_0230*IT_0441 + IT_0142*IT_0466 + 0.5*IT_0432
      *IT_0470) + IT_0333*(IT_0142*IT_0467 + IT_0432*IT_0471) + (-4)*conj
      (IT_0430)*(IT_0403*IT_0430 + -1./2*IT_0394*IT_0452 + -1./2*IT_0428*IT_0473
       + (-2)*IT_0408*IT_0475 + (-0.25)*IT_0456*IT_0477 + (-0.25)*IT_0478
      *IT_0479) + 8*IT_0430*(IT_0410*IT_0475 + 0.125*IT_0456*IT_0481 + 0.125
      *IT_0478*IT_0482) + 2*IT_0447*(conj(IT_0352)*IT_0405 + IT_0368*conj
      (IT_0405) + 0.5*conj(IT_0404)*IT_0419 + 0.5*IT_0404*IT_0422 + IT_0481
      *IT_0483 + IT_0477*IT_0484) + 2*IT_0473*(conj(IT_0404)*IT_0405 + IT_0404
      *conj(IT_0405) + IT_0482*IT_0483 + IT_0479*IT_0484) + 2*IT_0278*(conj
      (IT_0136)*IT_0231 + IT_0234*IT_0337 + IT_0469*IT_0485) + (-4)*conj(IT_0428
      )*(IT_0403*IT_0428 + -1./2*IT_0394*IT_0456 + -1./2*IT_0430*IT_0473 + (-2)
      *IT_0404*IT_0475 + 1./2*IT_0408*IT_0478 + 1./2*IT_0477*IT_0486) + 8
      *IT_0428*(conj(IT_0404)*IT_0475 + (-0.25)*IT_0410*IT_0478 + (-0.25)
      *IT_0481*IT_0486) + IT_0469*(IT_0367*IT_0432 + IT_0487) + 2*conj(IT_0228)*
      (IT_0142*IT_0228 + IT_0136*IT_0230 + IT_0183*IT_0278 + IT_0231*IT_0340 +
       IT_0368*IT_0432 + IT_0487) + IT_0488*(IT_0115*IT_0188 + IT_0489 + IT_0490
       + IT_0491) + 2*conj(IT_0115)*(IT_0188*IT_0435 + IT_0441*IT_0485 + IT_0489
       + IT_0490 + IT_0491);
    return create_ccomplex_return(IT_0492);
}

