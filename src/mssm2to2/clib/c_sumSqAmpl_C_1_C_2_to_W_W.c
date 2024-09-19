#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_C_2_to_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_C_2_to_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
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
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 0.5*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + (-0.5)*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = IT_0002*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + (-0.5)*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0028*IT_0035;
    const ccomplex_t IT_0037 = pow(m_C_2, 2);
    const ccomplex_t IT_0038 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0039 = IT_0036*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*IT_0039;
    const ccomplex_t IT_0041 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 0.5*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.5)*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0060 = IT_0019*IT_0059;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0000 + IT_0037 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0063 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0064 = IT_0004*IT_0063;
    const ccomplex_t IT_0065 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0066 = IT_0004*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + 0.5*IT_0067);
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0071 = IT_0004*IT_0070;
    const ccomplex_t IT_0072 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0073 = IT_0004*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + (-0.5)*IT_0074);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0069*IT_0076;
    const ccomplex_t IT_0078 = IT_0062*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = -IT_0040 + -IT_0058 + -IT_0061 + -IT_0079;
    const ccomplex_t IT_0081 = s_12*IT_0000;
    const ccomplex_t IT_0082 = pow(m_W, -4);
    const ccomplex_t IT_0083 = pow(m_W, 6);
    const ccomplex_t IT_0084 = s_12*IT_0083;
    const ccomplex_t IT_0085 = pow(m_W, 4);
    const ccomplex_t IT_0086 = s_13*s_23*IT_0085;
    const ccomplex_t IT_0087 = (-2)*IT_0086;
    const ccomplex_t IT_0088 = IT_0084 + IT_0087;
    const ccomplex_t IT_0089 = IT_0082*IT_0088;
    const ccomplex_t IT_0090 = 0.25*IT_0089;
    const ccomplex_t IT_0091 = IT_0081 + IT_0090;
    const ccomplex_t IT_0092 = 2*IT_0040 + 2*IT_0058 + 2*IT_0061 + 2*IT_0079;
    const ccomplex_t IT_0093 = s_12*s_14;
    const ccomplex_t IT_0094 = pow(m_W, -2);
    const ccomplex_t IT_0095 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0096 = s_24*IT_0000*IT_0001;
    const ccomplex_t IT_0097 = pow(s_13, 2);
    const ccomplex_t IT_0098 = s_24*IT_0097;
    const ccomplex_t IT_0099 = s_12*s_13*s_34;
    const ccomplex_t IT_0100 = -IT_0096 + 2*IT_0098 + (-2)*IT_0099;
    const ccomplex_t IT_0101 = IT_0095 + IT_0100;
    const ccomplex_t IT_0102 = IT_0094*IT_0101;
    const ccomplex_t IT_0103 = (-0.5)*IT_0102;
    const ccomplex_t IT_0104 = s_24*IT_0001;
    const ccomplex_t IT_0105 = -IT_0104;
    const ccomplex_t IT_0106 = IT_0093 + IT_0103 + IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0055*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0077*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0114 = IT_0036*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = (-2)*IT_0021 + (-2)*IT_0109 + (-2)*IT_0112 + (
      -2)*IT_0115;
    const ccomplex_t IT_0117 = s_12*s_24;
    const ccomplex_t IT_0118 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0119 = s_14*IT_0000*IT_0037;
    const ccomplex_t IT_0120 = s_13*s_23*s_24;
    const ccomplex_t IT_0121 = s_13*s_34*IT_0037;
    const ccomplex_t IT_0122 = -IT_0119 + (-2)*IT_0120 + 2*IT_0121;
    const ccomplex_t IT_0123 = IT_0118 + IT_0122;
    const ccomplex_t IT_0124 = IT_0094*IT_0123;
    const ccomplex_t IT_0125 = (-0.5)*IT_0124;
    const ccomplex_t IT_0126 = s_14*IT_0037;
    const ccomplex_t IT_0127 = -IT_0126;
    const ccomplex_t IT_0128 = IT_0117 + IT_0125 + IT_0127;
    const ccomplex_t IT_0129 = IT_0109 + IT_0112 + IT_0115;
    const ccomplex_t IT_0130 = conj(IT_0021) + conj(IT_0129);
    const ccomplex_t IT_0131 = s_14*s_24;
    const ccomplex_t IT_0132 = 0.5*IT_0089;
    const ccomplex_t IT_0133 = IT_0081 + (-4)*IT_0131 + IT_0132;
    const ccomplex_t IT_0134 = s_12*IT_0085;
    const ccomplex_t IT_0135 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0136 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0137 = s_13*s_24*s_34;
    const ccomplex_t IT_0138 = (-2)*IT_0135 + (-2)*IT_0136 + 4*IT_0137;
    const ccomplex_t IT_0139 = IT_0134 + IT_0138;
    const ccomplex_t IT_0140 = IT_0094*IT_0139;
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0133 + IT_0141;
    const ccomplex_t IT_0143 = m_C_2*IT_0094;
    const ccomplex_t IT_0144 = s_14*IT_0085;
    const ccomplex_t IT_0145 = IT_0094*IT_0144;
    const ccomplex_t IT_0146 = s_14*IT_0000;
    const ccomplex_t IT_0147 = s_13*s_34;
    const ccomplex_t IT_0148 = (-2)*IT_0146 + 4*IT_0147;
    const ccomplex_t IT_0149 = IT_0145 + IT_0148;
    const ccomplex_t IT_0150 = IT_0143*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0152 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0153 = IT_0004*IT_0152;
    const ccomplex_t IT_0154 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0155 = IT_0004*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + (-0.5)*IT_0156);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0018*IT_0158;
    const ccomplex_t IT_0160 = IT_0059*IT_0151*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0162 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0163 = IT_0004*IT_0162;
    const ccomplex_t IT_0164 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0165 = IT_0004*IT_0164;
    const ccomplex_t IT_0166 = 1.4142135623731*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*(IT_0163 + (-0.5)*IT_0166);
    const ccomplex_t IT_0168 = -IT_0167;
    const ccomplex_t IT_0169 = IT_0076*IT_0168;
    const ccomplex_t IT_0170 = IT_0062*IT_0161*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0172 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0173 = IT_0004*IT_0172;
    const ccomplex_t IT_0174 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0175 = IT_0004*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0173 + (-0.5)*IT_0176);
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = IT_0035*IT_0178;
    const ccomplex_t IT_0180 = IT_0038*IT_0171*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0182 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0183 = IT_0004*IT_0182;
    const ccomplex_t IT_0184 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0185 = IT_0004*IT_0184;
    const ccomplex_t IT_0186 = 1.4142135623731*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*(IT_0183 + (-0.5)*IT_0186);
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = IT_0054*IT_0188;
    const ccomplex_t IT_0190 = IT_0056*IT_0181*IT_0189;
    const ccomplex_t IT_0191 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0192 = IT_0004*IT_0191;
    const ccomplex_t IT_0193 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0194 = IT_0004*IT_0193;
    const ccomplex_t IT_0195 = 1.4142135623731*IT_0194;
    const ccomplex_t IT_0196 = (0 + _Complex_I*1)*(IT_0192 + 0.5*IT_0195);
    const ccomplex_t IT_0197 = -IT_0196;
    const ccomplex_t IT_0198 = IT_0158*IT_0197;
    const ccomplex_t IT_0199 = IT_0059*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*IT_0199;
    const ccomplex_t IT_0201 = m_C_1*IT_0200;
    const ccomplex_t IT_0202 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0203 = IT_0004*IT_0202;
    const ccomplex_t IT_0204 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0205 = IT_0004*IT_0204;
    const ccomplex_t IT_0206 = 1.4142135623731*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*(IT_0203 + 0.5*IT_0206);
    const ccomplex_t IT_0208 = -IT_0207;
    const ccomplex_t IT_0209 = IT_0188*IT_0208;
    const ccomplex_t IT_0210 = IT_0056*IT_0209;
    const ccomplex_t IT_0211 = (0 + _Complex_I*1)*IT_0210;
    const ccomplex_t IT_0212 = m_C_1*IT_0211;
    const ccomplex_t IT_0213 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0214 = IT_0004*IT_0213;
    const ccomplex_t IT_0215 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0216 = IT_0004*IT_0215;
    const ccomplex_t IT_0217 = 1.4142135623731*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*(IT_0214 + 0.5*IT_0217);
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0168*IT_0219;
    const ccomplex_t IT_0221 = IT_0062*IT_0220;
    const ccomplex_t IT_0222 = (0 + _Complex_I*1)*IT_0221;
    const ccomplex_t IT_0223 = m_C_1*IT_0222;
    const ccomplex_t IT_0224 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0225 = IT_0004*IT_0224;
    const ccomplex_t IT_0226 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0227 = IT_0004*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0225 + 0.5*IT_0228);
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = IT_0178*IT_0230;
    const ccomplex_t IT_0232 = IT_0038*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*IT_0232;
    const ccomplex_t IT_0234 = m_C_1*IT_0233;
    const ccomplex_t IT_0235 = -IT_0160 + -IT_0170 + -IT_0180 + -IT_0190 +
       IT_0201 + IT_0212 + IT_0223 + IT_0234;
    const ccomplex_t IT_0236 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0237 = (-2)*IT_0236;
    const ccomplex_t IT_0238 = IT_0144 + IT_0237;
    const ccomplex_t IT_0239 = IT_0082*IT_0238;
    const ccomplex_t IT_0240 = (-2)*IT_0147;
    const ccomplex_t IT_0241 = IT_0146 + IT_0240;
    const ccomplex_t IT_0242 = IT_0094*IT_0241;
    const ccomplex_t IT_0243 = 0.166666666666667*IT_0239 + (-0.666666666666667
      )*IT_0242;
    const ccomplex_t IT_0244 = s_14 + IT_0243;
    const ccomplex_t IT_0245 = m_C_2*IT_0244;
    const ccomplex_t IT_0246 = 6*IT_0245;
    const ccomplex_t IT_0247 = IT_0113*IT_0171*IT_0179;
    const ccomplex_t IT_0248 = IT_0002*IT_0151*IT_0159;
    const ccomplex_t IT_0249 = IT_0107*IT_0181*IT_0189;
    const ccomplex_t IT_0250 = IT_0110*IT_0161*IT_0169;
    const ccomplex_t IT_0251 = m_C_2*IT_0021;
    const ccomplex_t IT_0252 = m_C_2*IT_0109;
    const ccomplex_t IT_0253 = m_C_2*IT_0112;
    const ccomplex_t IT_0254 = m_C_2*IT_0115;
    const ccomplex_t IT_0255 = IT_0247 + IT_0248 + IT_0249 + IT_0250 + 
      -IT_0251 + -IT_0252 + -IT_0253 + -IT_0254;
    const ccomplex_t IT_0256 = 2*IT_0021;
    const ccomplex_t IT_0257 = m_C_1*m_C_2;
    const ccomplex_t IT_0258 = IT_0000*IT_0257;
    const ccomplex_t IT_0259 = 9*IT_0258;
    const ccomplex_t IT_0260 = IT_0107*IT_0209;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = IT_0110*IT_0220;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = IT_0113*IT_0231;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*IT_0264;
    const ccomplex_t IT_0266 = IT_0002*IT_0198;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = IT_0261 + IT_0263 + IT_0265 + IT_0267;
    const ccomplex_t IT_0269 = s_23*s_34*IT_0094;
    const ccomplex_t IT_0270 = s_24*IT_0085;
    const ccomplex_t IT_0271 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0272 = (-2)*IT_0271;
    const ccomplex_t IT_0273 = IT_0270 + IT_0272;
    const ccomplex_t IT_0274 = IT_0082*IT_0273;
    const ccomplex_t IT_0275 = (-0.666666666666667)*IT_0269 + (
      -0.166666666666667)*IT_0274;
    const ccomplex_t IT_0276 = s_24 + IT_0275;
    const ccomplex_t IT_0277 = m_C_1*IT_0276;
    const ccomplex_t IT_0278 = (-6)*IT_0277;
    const ccomplex_t IT_0279 = IT_0040 + IT_0058 + IT_0061 + IT_0079;
    const ccomplex_t IT_0280 = m_C_1*IT_0279;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = m_C_1*m_C_2*IT_0094;
    const ccomplex_t IT_0283 = pow(s_34, 2);
    const ccomplex_t IT_0284 = (-0.25)*IT_0085;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = IT_0282*IT_0285;
    const ccomplex_t IT_0287 = (-4)*IT_0286;
    const ccomplex_t IT_0288 = -IT_0200 + -IT_0211 + -IT_0222 + -IT_0233;
    const ccomplex_t IT_0289 = (-2)*IT_0021;
    const ccomplex_t IT_0290 = IT_0028*IT_0230;
    const ccomplex_t IT_0291 = IT_0113*IT_0171*IT_0290;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = (-2)*IT_0261 + (-2)*IT_0263 + (-2)*IT_0265 + (
      -2)*IT_0267;
    const ccomplex_t IT_0294 = s_14*s_24*IT_0094;
    const ccomplex_t IT_0295 = s_13*s_23;
    const ccomplex_t IT_0296 = s_13*s_24*s_34*IT_0094;
    const ccomplex_t IT_0297 = -IT_0131 + (-2)*IT_0295 + 2*IT_0296;
    const ccomplex_t IT_0298 = IT_0081 + IT_0297;
    const ccomplex_t IT_0299 = IT_0094*IT_0298;
    const ccomplex_t IT_0300 = -IT_0294 + (-0.5)*IT_0299;
    const ccomplex_t IT_0301 = s_12 + IT_0300;
    const ccomplex_t IT_0302 = m_C_2*IT_0301;
    const ccomplex_t IT_0303 = (-2)*IT_0302;
    const ccomplex_t IT_0304 = 2*IT_0303;
    const ccomplex_t IT_0305 = IT_0293*IT_0304;
    const ccomplex_t IT_0306 = (-2)*IT_0283;
    const ccomplex_t IT_0307 = IT_0085 + IT_0306;
    const ccomplex_t IT_0308 = IT_0082*IT_0307;
    const ccomplex_t IT_0309 = 0.25*IT_0308;
    const ccomplex_t IT_0310 = 1 + IT_0309;
    const ccomplex_t IT_0311 = IT_0257*IT_0310;
    const ccomplex_t IT_0312 = 8*IT_0311;
    const ccomplex_t IT_0313 = IT_0235*IT_0312;
    const ccomplex_t IT_0314 = (-2)*IT_0295;
    const ccomplex_t IT_0315 = IT_0081 + IT_0314;
    const ccomplex_t IT_0316 = IT_0094*IT_0315;
    const ccomplex_t IT_0317 = (-0.5)*IT_0316;
    const ccomplex_t IT_0318 = s_12 + IT_0317;
    const ccomplex_t IT_0319 = (-2)*IT_0131;
    const ccomplex_t IT_0320 = IT_0081 + IT_0319;
    const ccomplex_t IT_0321 = IT_0094*IT_0320;
    const ccomplex_t IT_0322 = (-0.5)*IT_0321;
    const ccomplex_t IT_0323 = IT_0082*IT_0139;
    const ccomplex_t IT_0324 = 0.25*IT_0323;
    const ccomplex_t IT_0325 = IT_0318 + IT_0322 + IT_0324;
    const ccomplex_t IT_0326 = 8*IT_0325;
    const ccomplex_t IT_0327 = IT_0292*IT_0326;
    const ccomplex_t IT_0328 = IT_0082*IT_0270;
    const ccomplex_t IT_0329 = (-10)*s_24;
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = m_C_1*IT_0330;
    const ccomplex_t IT_0332 = IT_0289*IT_0331;
    const ccomplex_t IT_0333 = IT_0069*IT_0219;
    const ccomplex_t IT_0334 = IT_0110*IT_0161*IT_0333;
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0038*IT_0171*IT_0290;
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = s_13*s_23*IT_0094;
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = s_14*s_23*s_34;
    const ccomplex_t IT_0341 = s_12*IT_0283;
    const ccomplex_t IT_0342 = (-2)*IT_0135 + (-2)*IT_0136 + 2*IT_0137 + 2
      *IT_0340 + (-2)*IT_0341;
    const ccomplex_t IT_0343 = IT_0134 + IT_0342;
    const ccomplex_t IT_0344 = IT_0082*IT_0343;
    const ccomplex_t IT_0345 = (-0.25)*IT_0344;
    const ccomplex_t IT_0346 = -IT_0294;
    const ccomplex_t IT_0347 = s_12 + IT_0339 + IT_0345 + IT_0346;
    const ccomplex_t IT_0348 = 8*IT_0347;
    const ccomplex_t IT_0349 = IT_0337*IT_0348;
    const ccomplex_t IT_0350 = IT_0326*IT_0335;
    const ccomplex_t IT_0351 = IT_0047*IT_0208;
    const ccomplex_t IT_0352 = IT_0107*IT_0181*IT_0351;
    const ccomplex_t IT_0353 = m_C_2*IT_0267;
    const ccomplex_t IT_0354 = m_C_2*IT_0263;
    const ccomplex_t IT_0355 = m_C_2*IT_0265;
    const ccomplex_t IT_0356 = -IT_0352 + IT_0353 + IT_0354 + IT_0355;
    const ccomplex_t IT_0357 = IT_0011*IT_0197;
    const ccomplex_t IT_0358 = IT_0002*IT_0151*IT_0357;
    const ccomplex_t IT_0359 = m_C_2*IT_0261;
    const ccomplex_t IT_0360 = -IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = conj(IT_0356) + conj(IT_0360);
    const ccomplex_t IT_0362 = IT_0056*IT_0181*IT_0351;
    const ccomplex_t IT_0363 = IT_0062*IT_0161*IT_0333;
    const ccomplex_t IT_0364 = -IT_0362 + -IT_0363;
    const ccomplex_t IT_0365 = IT_0348*IT_0364;
    const ccomplex_t IT_0366 = IT_0326*IT_0356;
    const ccomplex_t IT_0367 = 2*IT_0200 + 2*IT_0211 + 2*IT_0222 + 2*IT_0233;
    const ccomplex_t IT_0368 = (s_12 + -IT_0294)*(s_12 + (-0.5)*IT_0316);
    const ccomplex_t IT_0369 = (-2)*IT_0368;
    const ccomplex_t IT_0370 = 2*IT_0369;
    const ccomplex_t IT_0371 = -IT_0131;
    const ccomplex_t IT_0372 = IT_0081 + IT_0371;
    const ccomplex_t IT_0373 = IT_0094*IT_0372;
    const ccomplex_t IT_0374 = -IT_0294 + (-0.25)*IT_0373;
    const ccomplex_t IT_0375 = s_12 + IT_0374;
    const ccomplex_t IT_0376 = m_C_2*IT_0375;
    const ccomplex_t IT_0377 = (-4)*IT_0376;
    const ccomplex_t IT_0378 = 2*IT_0377;
    const ccomplex_t IT_0379 = pow(s_24, 2);
    const ccomplex_t IT_0380 = IT_0000*IT_0037;
    const ccomplex_t IT_0381 = -IT_0380;
    const ccomplex_t IT_0382 = IT_0379 + IT_0381;
    const ccomplex_t IT_0383 = IT_0094*IT_0382;
    const ccomplex_t IT_0384 = IT_0094*IT_0379;
    const ccomplex_t IT_0385 = -IT_0037;
    const ccomplex_t IT_0386 = IT_0384 + IT_0385;
    const ccomplex_t IT_0387 = (-4)*IT_0386;
    const ccomplex_t IT_0388 = IT_0383 + IT_0387;
    const ccomplex_t IT_0389 = m_C_1*IT_0388;
    const ccomplex_t IT_0390 = 2*IT_0389;
    const ccomplex_t IT_0391 = pow(s_23, 2);
    const ccomplex_t IT_0392 = s_14*IT_0391;
    const ccomplex_t IT_0393 = s_12*s_23*s_34;
    const ccomplex_t IT_0394 = 0.5*IT_0118 + (-0.5)*IT_0119 + -IT_0393;
    const ccomplex_t IT_0395 = IT_0392 + IT_0394;
    const ccomplex_t IT_0396 = IT_0094*IT_0395;
    const ccomplex_t IT_0397 = -IT_0117;
    const ccomplex_t IT_0398 = IT_0126 + IT_0396 + IT_0397;
    const ccomplex_t IT_0399 = 4*IT_0398;
    const ccomplex_t IT_0400 = s_23*s_24*s_34*IT_0094;
    const ccomplex_t IT_0401 = 0.5*IT_0379 + (-0.5)*IT_0380 + -IT_0400;
    const ccomplex_t IT_0402 = IT_0391 + IT_0401;
    const ccomplex_t IT_0403 = IT_0094*IT_0402;
    const ccomplex_t IT_0404 = -IT_0384;
    const ccomplex_t IT_0405 = IT_0037 + IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = m_C_1*IT_0405;
    const ccomplex_t IT_0407 = 4*IT_0406;
    const ccomplex_t IT_0408 = 4*IT_0128;
    const ccomplex_t IT_0409 = IT_0257*IT_0386;
    const ccomplex_t IT_0410 = 6*IT_0409;
    const ccomplex_t IT_0411 = s_12 + IT_0346;
    const ccomplex_t IT_0412 = IT_0257*IT_0411;
    const ccomplex_t IT_0413 = (-6)*IT_0412;
    const ccomplex_t IT_0414 = -(s_12 + (-0.5)*IT_0316)*(IT_0037 + -IT_0384);
    const ccomplex_t IT_0415 = (-2)*IT_0414;
    const ccomplex_t IT_0416 = (-2)*IT_0415;
    const ccomplex_t IT_0417 = 2*conj(IT_0021);
    const ccomplex_t IT_0418 = 2*IT_0150;
    const ccomplex_t IT_0419 = 2*IT_0278;
    const ccomplex_t IT_0420 = pow(s_14, 2);
    const ccomplex_t IT_0421 = IT_0000*IT_0001;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0420 + IT_0422;
    const ccomplex_t IT_0424 = IT_0094*IT_0423;
    const ccomplex_t IT_0425 = IT_0094*IT_0420;
    const ccomplex_t IT_0426 = -IT_0001;
    const ccomplex_t IT_0427 = IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = (-4)*IT_0427;
    const ccomplex_t IT_0429 = IT_0424 + IT_0428;
    const ccomplex_t IT_0430 = m_C_2*IT_0429;
    const ccomplex_t IT_0431 = 2*IT_0430;
    const ccomplex_t IT_0432 = s_24*IT_0000;
    const ccomplex_t IT_0433 = s_23*s_34;
    const ccomplex_t IT_0434 = (-2)*IT_0433;
    const ccomplex_t IT_0435 = IT_0432 + IT_0434;
    const ccomplex_t IT_0436 = IT_0094*IT_0435;
    const ccomplex_t IT_0437 = 0.166666666666667*IT_0274 + (-0.666666666666667
      )*IT_0436;
    const ccomplex_t IT_0438 = s_24 + IT_0437;
    const ccomplex_t IT_0439 = m_C_1*IT_0438;
    const ccomplex_t IT_0440 = 6*IT_0439;
    const ccomplex_t IT_0441 = 2*IT_0440;
    const ccomplex_t IT_0442 = IT_0082*IT_0144;
    const ccomplex_t IT_0443 = (-10)*s_14;
    const ccomplex_t IT_0444 = IT_0442 + IT_0443;
    const ccomplex_t IT_0445 = m_C_2*IT_0444;
    const ccomplex_t IT_0446 = 2*IT_0445;
    const ccomplex_t IT_0447 = (-2)*IT_0135 + (-2)*IT_0136 + 4*IT_0340;
    const ccomplex_t IT_0448 = IT_0134 + IT_0447;
    const ccomplex_t IT_0449 = IT_0082*IT_0448;
    const ccomplex_t IT_0450 = 0.25*IT_0449;
    const ccomplex_t IT_0451 = IT_0318 + IT_0322 + IT_0450;
    const ccomplex_t IT_0452 = 8*IT_0451;
    const ccomplex_t IT_0453 = 1.125*IT_0257;
    const ccomplex_t IT_0454 = 16*IT_0453;
    const ccomplex_t IT_0455 = (-8)*IT_0347;
    const ccomplex_t IT_0456 = s_14*s_23*s_34*IT_0094;
    const ccomplex_t IT_0457 = -IT_0131 + (-2)*IT_0295 + 2*IT_0456;
    const ccomplex_t IT_0458 = IT_0081 + IT_0457;
    const ccomplex_t IT_0459 = IT_0094*IT_0458;
    const ccomplex_t IT_0460 = -IT_0294 + (-0.5)*IT_0459;
    const ccomplex_t IT_0461 = s_12 + IT_0460;
    const ccomplex_t IT_0462 = m_C_1*IT_0461;
    const ccomplex_t IT_0463 = (-2)*IT_0462;
    const ccomplex_t IT_0464 = (-2)*IT_0463;
    const ccomplex_t IT_0465 = (-2)*IT_0377;
    const ccomplex_t IT_0466 = conj(IT_0293)*IT_0407 + IT_0150*IT_0417 + conj
      (IT_0129)*IT_0418 + conj(IT_0268)*IT_0419 + conj(IT_0092)*IT_0431 + conj
      (IT_0288)*IT_0441 + conj(IT_0080)*IT_0446 + conj(IT_0235)*IT_0452 + conj
      (IT_0281)*IT_0454 + conj(IT_0255)*IT_0455 + conj(IT_0367)*IT_0464 + conj
      (IT_0116)*IT_0465;
    const ccomplex_t IT_0467 = IT_0256*IT_0278;
    const ccomplex_t IT_0468 = IT_0367*IT_0431;
    const ccomplex_t IT_0469 = IT_0292*IT_0348;
    const ccomplex_t IT_0470 = IT_0467 + IT_0468 + IT_0469;
    const ccomplex_t IT_0471 = IT_0337*IT_0452;
    const ccomplex_t IT_0472 = (-16)*IT_0453;
    const ccomplex_t IT_0473 = IT_0235*IT_0472;
    const ccomplex_t IT_0474 = IT_0293*IT_0465;
    const ccomplex_t IT_0475 = IT_0059*IT_0151*IT_0357;
    const ccomplex_t IT_0476 = -IT_0475;
    const ccomplex_t IT_0477 = conj(IT_0364) + conj(IT_0476);
    const ccomplex_t IT_0478 = IT_0335*IT_0348;
    const ccomplex_t IT_0479 = IT_0364*IT_0452;
    const ccomplex_t IT_0480 = 2*IT_0331;
    const ccomplex_t IT_0481 = 2*IT_0246;
    const ccomplex_t IT_0482 = s_13*s_34*IT_0094;
    const ccomplex_t IT_0483 = (-0.166666666666667)*IT_0239 + (
      -0.666666666666667)*IT_0482;
    const ccomplex_t IT_0484 = s_14 + IT_0483;
    const ccomplex_t IT_0485 = m_C_2*IT_0484;
    const ccomplex_t IT_0486 = (-6)*IT_0485;
    const ccomplex_t IT_0487 = 2*IT_0486;
    const ccomplex_t IT_0488 = m_C_1*IT_0094;
    const ccomplex_t IT_0489 = IT_0094*IT_0270;
    const ccomplex_t IT_0490 = (-2)*IT_0433 + (-0.5)*IT_0489;
    const ccomplex_t IT_0491 = IT_0432 + IT_0490;
    const ccomplex_t IT_0492 = IT_0488*IT_0491;
    const ccomplex_t IT_0493 = (-2)*IT_0492;
    const ccomplex_t IT_0494 = 2*IT_0493;
    const ccomplex_t IT_0495 = s_13*s_14*s_34*IT_0094;
    const ccomplex_t IT_0496 = 0.5*IT_0420 + (-0.5)*IT_0421 + -IT_0495;
    const ccomplex_t IT_0497 = IT_0097 + IT_0496;
    const ccomplex_t IT_0498 = IT_0094*IT_0497;
    const ccomplex_t IT_0499 = -IT_0425;
    const ccomplex_t IT_0500 = IT_0001 + IT_0498 + IT_0499;
    const ccomplex_t IT_0501 = m_C_2*IT_0500;
    const ccomplex_t IT_0502 = 4*IT_0501;
    const ccomplex_t IT_0503 = conj(IT_0292) + conj(IT_0335) + conj(IT_0356) +
       conj(IT_0360);
    const ccomplex_t IT_0504 = conj(IT_0337) + conj(IT_0364) + conj(IT_0476);
    const ccomplex_t IT_0505 = (-8)*IT_0311;
    const ccomplex_t IT_0506 = (-2)*IT_0303;
    const ccomplex_t IT_0507 = m_C_1*IT_0375;
    const ccomplex_t IT_0508 = (-4)*IT_0507;
    const ccomplex_t IT_0509 = (-2)*IT_0508;
    const ccomplex_t IT_0510 = conj(IT_0281)*IT_0312 + conj(IT_0255)*IT_0326 +
       conj(IT_0293)*IT_0390 + IT_0246*IT_0417 + conj(IT_0235)*IT_0455 + conj
      (IT_0268)*IT_0480 + conj(IT_0129)*IT_0481 + conj(IT_0080)*IT_0487 + conj
      (IT_0288)*IT_0494 + conj(IT_0092)*IT_0502 + IT_0454*IT_0503 + IT_0504
      *IT_0505 + conj(IT_0116)*IT_0506 + conj(IT_0367)*IT_0509;
    const ccomplex_t IT_0511 = IT_0094*IT_0448;
    const ccomplex_t IT_0512 = -IT_0511;
    const ccomplex_t IT_0513 = IT_0133 + IT_0512;
    const ccomplex_t IT_0514 = (-8)*IT_0091;
    const ccomplex_t IT_0515 = s_13*s_14*s_23;
    const ccomplex_t IT_0516 = s_23*s_34*IT_0001;
    const ccomplex_t IT_0517 = -IT_0096 + (-2)*IT_0515 + 2*IT_0516;
    const ccomplex_t IT_0518 = IT_0095 + IT_0517;
    const ccomplex_t IT_0519 = IT_0094*IT_0518;
    const ccomplex_t IT_0520 = 2*IT_0104;
    const ccomplex_t IT_0521 = (-2)*IT_0093;
    const ccomplex_t IT_0522 = IT_0519 + IT_0520 + IT_0521;
    const ccomplex_t IT_0523 = (-2)*IT_0522;
    const ccomplex_t IT_0524 = (-2)*IT_0287;
    const ccomplex_t IT_0525 = (-2)*IT_0440;
    const ccomplex_t IT_0526 = (-2)*IT_0493;
    const ccomplex_t IT_0527 = (-2)*IT_0259;
    const ccomplex_t IT_0528 = conj(IT_0021)*IT_0091 + (-0.125)*conj(IT_0116)
      *IT_0399 + (-0.125)*conj(IT_0235)*IT_0446 + (-0.125)*conj(IT_0255)*IT_0487
       + (-0.125)*IT_0441*IT_0504 + 0.5*conj(IT_0080)*IT_0513 + (-0.125)*conj
      (IT_0129)*IT_0514 + (-0.125)*conj(IT_0092)*IT_0523 + (-0.125)*conj(IT_0268
      )*IT_0524 + (-0.125)*conj(IT_0281)*IT_0525 + (-0.125)*IT_0503*IT_0526 + (
      -0.125)*conj(IT_0288)*IT_0527;
    const ccomplex_t IT_0529 = (-2)*conj(IT_0021);
    const ccomplex_t IT_0530 = (-2)*IT_0486;
    const ccomplex_t IT_0531 = (-2)*IT_0445;
    const ccomplex_t IT_0532 = (-0.25)*conj(IT_0293)*IT_0399 + (-0.25)*conj
      (IT_0235)*IT_0441 + (-0.25)*conj(IT_0255)*IT_0494 + (-0.25)*IT_0446
      *IT_0504 + conj(IT_0288)*IT_0513 + (-0.25)*conj(IT_0268)*IT_0514 + (-0.25)
      *conj(IT_0367)*IT_0523 + (-0.25)*conj(IT_0129)*IT_0524 + (-0.25)*conj
      (IT_0080)*IT_0527 + (-0.25)*IT_0287*IT_0529 + (-0.25)*IT_0503*IT_0530 + (
      -0.25)*conj(IT_0281)*IT_0531;
    const ccomplex_t IT_0533 = -(s_12 + (-0.5)*IT_0316)*(IT_0001 + -IT_0425);
    const ccomplex_t IT_0534 = (-2)*IT_0533;
    const ccomplex_t IT_0535 = (-2)*IT_0534;
    const ccomplex_t IT_0536 = (-0.25)*conj(IT_0367);
    const ccomplex_t IT_0537 = (-0.25)*conj(IT_0293);
    const ccomplex_t IT_0538 = (-0.25)*conj(IT_0235);
    const ccomplex_t IT_0539 = IT_0257*IT_0427;
    const ccomplex_t IT_0540 = 6*IT_0539;
    const ccomplex_t IT_0541 = (-0.25)*conj(IT_0092);
    const ccomplex_t IT_0542 = (-0.25)*conj(IT_0116);
    const ccomplex_t IT_0543 = (-0.25)*conj(IT_0255);
    const ccomplex_t IT_0544 = (-2)*IT_0430;
    const ccomplex_t IT_0545 = (-0.25)*conj(IT_0281);
    const ccomplex_t IT_0546 = (-0.25)*conj(IT_0288);
    const ccomplex_t IT_0547 = (-4)*IT_0367*(IT_0106*conj(IT_0268) + IT_0501
      *IT_0503 + IT_0535*IT_0536 + IT_0370*IT_0537 + IT_0464*IT_0538 + IT_0540
      *IT_0541 + IT_0413*IT_0542 + IT_0509*IT_0543 + IT_0544*IT_0545 + IT_0523
      *IT_0546);
    const ccomplex_t IT_0548 = (-0.25)*IT_0504;
    const ccomplex_t IT_0549 = (-2)*IT_0389;
    const ccomplex_t IT_0550 = (-0.25)*IT_0503;
    const ccomplex_t IT_0551 = (-0.25)*conj(IT_0129);
    const ccomplex_t IT_0552 = (-0.25)*conj(IT_0080);
    const ccomplex_t IT_0553 = -conj(IT_0021)*IT_0128 + conj(IT_0281)*IT_0406 
      + IT_0413*IT_0536 + IT_0410*IT_0537 + IT_0465*IT_0538 + IT_0370*IT_0541 +
       IT_0416*IT_0542 + IT_0506*IT_0543 + IT_0407*IT_0548 + IT_0549*IT_0550 +
       IT_0408*IT_0551 + IT_0399*IT_0552;
    const ccomplex_t IT_0554 = (-0.25)*conj(IT_0268);
    const ccomplex_t IT_0555 = (-2)*IT_0331;
    const ccomplex_t IT_0556 = (-2)*IT_0278;
    const ccomplex_t IT_0557 = conj(IT_0092)*IT_0106 + IT_0130*IT_0142 +
       IT_0418*IT_0538 + IT_0408*IT_0542 + IT_0481*IT_0543 + IT_0524*IT_0546 +
       IT_0419*IT_0548 + IT_0514*IT_0552 + IT_0527*IT_0554 + IT_0550*IT_0555 +
       IT_0545*IT_0556;
    const ccomplex_t IT_0558 = 2*IT_0508;
    const ccomplex_t IT_0559 = 2*IT_0463;
    const ccomplex_t IT_0560 = IT_0106*IT_0130 + IT_0413*IT_0537 + IT_0431
      *IT_0538 + IT_0536*IT_0540 + IT_0535*IT_0541 + IT_0370*IT_0542 + IT_0502
      *IT_0543 + IT_0464*IT_0548 + IT_0523*IT_0552 + IT_0550*IT_0558 + IT_0545
      *IT_0559;
    const ccomplex_t IT_0561 = (-2)*IT_0246;
    const ccomplex_t IT_0562 = (-2)*IT_0150;
    const ccomplex_t IT_0563 = IT_0142*conj(IT_0268) + IT_0106*conj(IT_0367) +
       (-0.25)*IT_0259*IT_0529 + IT_0408*IT_0537 + IT_0419*IT_0538 + IT_0480
      *IT_0543 + IT_0514*IT_0546 + IT_0418*IT_0548 + IT_0527*IT_0551 + IT_0524
      *IT_0552 + IT_0550*IT_0561 + IT_0545*IT_0562;
    const ccomplex_t IT_0564 = (-8)*IT_0451;
    const ccomplex_t IT_0565 = conj(IT_0116)*IT_0406 + (-0.25)*IT_0278*IT_0529
       + IT_0378*IT_0537 + IT_0454*IT_0538 + IT_0312*IT_0543 + IT_0536*IT_0544 +
       IT_0452*IT_0545 + IT_0531*IT_0546 + IT_0455*IT_0550 + IT_0525*IT_0552 +
       IT_0551*IT_0556 + IT_0541*IT_0559 + IT_0554*IT_0562 + IT_0548*IT_0564;
    const ccomplex_t IT_0566 = IT_0554*IT_0561;
    const ccomplex_t IT_0567 = IT_0541*IT_0558;
    const ccomplex_t IT_0568 = IT_0542*IT_0549;
    const ccomplex_t IT_0569 = IT_0454*IT_0543;
    const ccomplex_t IT_0570 = IT_0455*IT_0545;
    const ccomplex_t IT_0571 = IT_0530*IT_0546;
    const ccomplex_t IT_0572 = IT_0551*IT_0555;
    const ccomplex_t IT_0573 = IT_0526*IT_0552;
    const ccomplex_t IT_0574 = conj(IT_0293)*IT_0304;
    const ccomplex_t IT_0575 = (-0.25)*IT_0574;
    const ccomplex_t IT_0576 = conj(IT_0235)*IT_0312;
    const ccomplex_t IT_0577 = (-0.25)*IT_0576;
    const ccomplex_t IT_0578 = IT_0331*IT_0529;
    const ccomplex_t IT_0579 = (-0.25)*IT_0578;
    const ccomplex_t IT_0580 = conj(IT_0367)*IT_0501;
    const ccomplex_t IT_0581 = conj(IT_0337)*IT_0348;
    const ccomplex_t IT_0582 = (-0.25)*IT_0581;
    const ccomplex_t IT_0583 = conj(IT_0292)*IT_0326;
    const ccomplex_t IT_0584 = (-0.25)*IT_0583;
    const ccomplex_t IT_0585 = IT_0348*conj(IT_0476);
    const ccomplex_t IT_0586 = (-0.25)*IT_0585;
    const ccomplex_t IT_0587 = IT_0348*conj(IT_0364);
    const ccomplex_t IT_0588 = (-0.25)*IT_0587;
    const ccomplex_t IT_0589 = IT_0326*conj(IT_0335);
    const ccomplex_t IT_0590 = (-0.25)*IT_0589;
    const ccomplex_t IT_0591 = IT_0326*IT_0361 + (-4)*IT_0566 + (-4)*IT_0567 +
       (-4)*IT_0568 + (-4)*IT_0569 + (-4)*IT_0570 + (-4)*IT_0571 + (-4)*IT_0572 
      + (-4)*IT_0573 + (-4)*IT_0575 + (-4)*IT_0577 + (-4)*IT_0579 + (-4)*IT_0580
       + (-4)*IT_0582 + (-4)*IT_0584 + (-4)*IT_0586 + (-4)*IT_0588 + (-4)*IT_0590;
    const ccomplex_t IT_0592 = IT_0278*IT_0417;
    const ccomplex_t IT_0593 = conj(IT_0367)*IT_0431;
    const ccomplex_t IT_0594 = conj(IT_0292)*IT_0348;
    const ccomplex_t IT_0595 = conj(IT_0235)*IT_0472;
    const ccomplex_t IT_0596 = conj(IT_0293)*IT_0465;
    const ccomplex_t IT_0597 = conj(IT_0268)*IT_0418;
    const ccomplex_t IT_0598 = conj(IT_0092)*IT_0464;
    const ccomplex_t IT_0599 = conj(IT_0116)*IT_0407;
    const ccomplex_t IT_0600 = conj(IT_0255)*IT_0505;
    const ccomplex_t IT_0601 = conj(IT_0281)*IT_0564;
    const ccomplex_t IT_0602 = conj(IT_0288)*IT_0446;
    const ccomplex_t IT_0603 = conj(IT_0129)*IT_0419;
    const ccomplex_t IT_0604 = conj(IT_0080)*IT_0441;
    const ccomplex_t IT_0605 = conj(IT_0335)*IT_0348;
    const ccomplex_t IT_0606 = conj(IT_0337)*IT_0452;
    const ccomplex_t IT_0607 = (-8)*IT_0021*(conj(IT_0080)*IT_0091 + 0.5*conj
      (IT_0092)*IT_0106 + (-0.5)*conj(IT_0116)*IT_0128 + 0.5*IT_0130*IT_0142) + 
      (IT_0150*conj(IT_0235) + IT_0246*conj(IT_0255))*IT_0256 + (IT_0259*conj
      (IT_0268) + IT_0278*conj(IT_0281) + IT_0287*conj(IT_0288))*IT_0289 + conj
      (IT_0292)*(IT_0305 + IT_0313 + IT_0327 + IT_0332) + conj(IT_0335)*(IT_0305
       + IT_0313 + IT_0327 + IT_0332 + IT_0349 + IT_0350) + IT_0361*(IT_0305 +
       IT_0313 + IT_0327 + IT_0332 + IT_0349 + IT_0350 + IT_0365 + IT_0366) +
       IT_0293*(conj(IT_0367)*IT_0370 + conj(IT_0281)*IT_0378 + conj(IT_0255)
      *IT_0390 + conj(IT_0288)*IT_0399 + conj(IT_0235)*IT_0407 + conj(IT_0268)
      *IT_0408 + conj(IT_0116)*IT_0410 + conj(IT_0092)*IT_0413 + conj(IT_0293)
      *IT_0416) + IT_0235*IT_0466 + conj(IT_0337)*(IT_0470 + IT_0471 + IT_0473 +
       IT_0474) + IT_0477*(IT_0470 + IT_0471 + IT_0473 + IT_0474 + IT_0478 +
       IT_0479) + IT_0255*IT_0510 + (-8)*IT_0080*IT_0528 + (-4)*IT_0288*IT_0532 
      + IT_0547 + (-4)*IT_0116*IT_0553 + (-4)*IT_0129*IT_0557 + (-4)*IT_0092
      *IT_0560 + (-4)*IT_0268*IT_0563 + (-4)*IT_0281*IT_0565 + (-4)*IT_0292*
      (IT_0566 + IT_0567 + IT_0568 + IT_0569 + IT_0570 + IT_0571 + IT_0572 +
       IT_0573 + IT_0575 + IT_0577 + IT_0579 + IT_0580) + (-4)*IT_0335*(IT_0566 
      + IT_0567 + IT_0568 + IT_0569 + IT_0570 + IT_0571 + IT_0572 + IT_0573 +
       IT_0575 + IT_0577 + IT_0579 + IT_0580 + IT_0582 + IT_0584) + (-4)*IT_0356
      *(IT_0566 + IT_0567 + IT_0568 + IT_0569 + IT_0570 + IT_0571 + IT_0572 +
       IT_0573 + IT_0575 + IT_0577 + IT_0579 + IT_0580 + IT_0582 + IT_0584 +
       IT_0586 + IT_0588 + IT_0590) + IT_0360*IT_0591 + IT_0337*(IT_0592 +
       IT_0593 + IT_0594 + IT_0595 + IT_0596 + IT_0597 + IT_0598 + IT_0599 +
       IT_0600 + IT_0601 + IT_0602 + IT_0603 + IT_0604) + IT_0364*(IT_0592 +
       IT_0593 + IT_0594 + IT_0595 + IT_0596 + IT_0597 + IT_0598 + IT_0599 +
       IT_0600 + IT_0601 + IT_0602 + IT_0603 + IT_0604 + IT_0605 + IT_0606) +
       IT_0476*(IT_0348*IT_0361 + IT_0452*IT_0477 + IT_0592 + IT_0593 + IT_0594 
      + IT_0595 + IT_0596 + IT_0597 + IT_0598 + IT_0599 + IT_0600 + IT_0601 +
       IT_0602 + IT_0603 + IT_0604 + IT_0605 + IT_0606);
    return create_ccomplex_return(IT_0607);
}

