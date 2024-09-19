#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_C_2_to_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_C_2_to_W_W(
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_C_2, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 0.5*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + (-0.5)*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = IT_0002*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0023 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0024 = IT_0004*IT_0023;
    const ccomplex_t IT_0025 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0026 = IT_0004*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + 0.5*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0031 = IT_0004*IT_0030;
    const ccomplex_t IT_0032 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0033 = IT_0004*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0031 + (-0.5)*IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = IT_0029*IT_0036;
    const ccomplex_t IT_0038 = IT_0022*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0041 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0042 = IT_0004*IT_0041;
    const ccomplex_t IT_0043 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0044 = IT_0004*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + 0.5*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0049 = IT_0004*IT_0048;
    const ccomplex_t IT_0050 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0051 = IT_0004*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.5)*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = IT_0040*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0059 = IT_0004*IT_0058;
    const ccomplex_t IT_0060 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0061 = IT_0004*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + 0.5*IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0066 = IT_0004*IT_0065;
    const ccomplex_t IT_0067 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0068 = IT_0004*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-0.5)*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = IT_0064*IT_0071;
    const ccomplex_t IT_0073 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0074 = IT_0072*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = (-2)*IT_0021 + (-2)*IT_0039 + (-2)*IT_0057 + (
      -2)*IT_0075;
    const ccomplex_t IT_0077 = pow(m_W, -2);
    const ccomplex_t IT_0078 = s_14*s_24*IT_0077;
    const ccomplex_t IT_0079 = s_12*IT_0000;
    const ccomplex_t IT_0080 = s_13*s_23;
    const ccomplex_t IT_0081 = (-2)*IT_0080;
    const ccomplex_t IT_0082 = IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = IT_0077*IT_0082;
    const ccomplex_t IT_0084 = (s_12 + -IT_0078)*(s_12 + (-0.5)*IT_0083);
    const ccomplex_t IT_0085 = (-2)*IT_0084;
    const ccomplex_t IT_0086 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0087 = IT_0072*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0019*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0093 = IT_0055*IT_0092;
    const ccomplex_t IT_0094 = (0 + _Complex_I*1)*IT_0093;
    const ccomplex_t IT_0095 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0096 = IT_0037*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = 2*IT_0088 + 2*IT_0091 + 2*IT_0094 + 2*IT_0097;
    const ccomplex_t IT_0099 = 2*conj(IT_0098);
    const ccomplex_t IT_0100 = 1.125*IT_0001;
    const ccomplex_t IT_0101 = cpow(IT_0071, 2);
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*m_N_1*IT_0101;
    const ccomplex_t IT_0103 = IT_0073*IT_0102;
    const ccomplex_t IT_0104 = cpow(IT_0054, 2);
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*m_N_3*IT_0104;
    const ccomplex_t IT_0106 = IT_0040*IT_0105;
    const ccomplex_t IT_0107 = cpow(IT_0018, 2);
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_2*IT_0107;
    const ccomplex_t IT_0109 = IT_0002*IT_0108;
    const ccomplex_t IT_0110 = cpow(IT_0036, 2);
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*m_N_4*IT_0110;
    const ccomplex_t IT_0112 = IT_0022*IT_0111;
    const ccomplex_t IT_0113 = m_C_2*IT_0075;
    const ccomplex_t IT_0114 = m_C_2*IT_0021;
    const ccomplex_t IT_0115 = m_C_2*IT_0057;
    const ccomplex_t IT_0116 = m_C_2*IT_0039;
    const ccomplex_t IT_0117 = IT_0103 + IT_0106 + IT_0109 + IT_0112 + 
      -IT_0113 + -IT_0114 + -IT_0115 + -IT_0116;
    const ccomplex_t IT_0118 = cpow(IT_0029, 2);
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*m_N_4*IT_0118;
    const ccomplex_t IT_0120 = IT_0022*IT_0119;
    const ccomplex_t IT_0121 = cpow(IT_0064, 2);
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*m_N_1*IT_0121;
    const ccomplex_t IT_0123 = IT_0073*IT_0122;
    const ccomplex_t IT_0124 = cpow(IT_0011, 2);
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*m_N_2*IT_0124;
    const ccomplex_t IT_0126 = IT_0002*IT_0125;
    const ccomplex_t IT_0127 = cpow(IT_0047, 2);
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_N_3*IT_0127;
    const ccomplex_t IT_0129 = IT_0040*IT_0128;
    const ccomplex_t IT_0130 = IT_0113 + IT_0114 + IT_0115 + IT_0116 + 
      -IT_0120 + -IT_0123 + -IT_0126 + -IT_0129;
    const ccomplex_t IT_0131 = IT_0086*IT_0102;
    const ccomplex_t IT_0132 = IT_0089*IT_0108;
    const ccomplex_t IT_0133 = IT_0092*IT_0105;
    const ccomplex_t IT_0134 = IT_0095*IT_0111;
    const ccomplex_t IT_0135 = m_C_2*IT_0088;
    const ccomplex_t IT_0136 = m_C_2*IT_0091;
    const ccomplex_t IT_0137 = m_C_2*IT_0094;
    const ccomplex_t IT_0138 = m_C_2*IT_0097;
    const ccomplex_t IT_0139 = -IT_0131 + -IT_0132 + -IT_0133 + -IT_0134 +
       IT_0135 + IT_0136 + IT_0137 + IT_0138;
    const ccomplex_t IT_0140 = IT_0088 + IT_0091 + IT_0094 + IT_0097;
    const ccomplex_t IT_0141 = m_C_2*IT_0140;
    const ccomplex_t IT_0142 = -IT_0141;
    const ccomplex_t IT_0143 = IT_0092*IT_0128;
    const ccomplex_t IT_0144 = IT_0095*IT_0119;
    const ccomplex_t IT_0145 = -IT_0143 + -IT_0144;
    const ccomplex_t IT_0146 = IT_0089*IT_0125;
    const ccomplex_t IT_0147 = IT_0086*IT_0122;
    const ccomplex_t IT_0148 = -IT_0146 + -IT_0147;
    const ccomplex_t IT_0149 = IT_0145 + IT_0148;
    const ccomplex_t IT_0150 = (-16)*conj(IT_0139);
    const ccomplex_t IT_0151 = s_13*s_23*IT_0077;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = pow(m_W, -4);
    const ccomplex_t IT_0154 = pow(m_W, 4);
    const ccomplex_t IT_0155 = s_12*IT_0154;
    const ccomplex_t IT_0156 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0157 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0158 = s_13*s_24*s_34;
    const ccomplex_t IT_0159 = s_14*s_23*s_34;
    const ccomplex_t IT_0160 = pow(s_34, 2);
    const ccomplex_t IT_0161 = s_12*IT_0160;
    const ccomplex_t IT_0162 = (-2)*IT_0156 + (-2)*IT_0157 + 2*IT_0158 + 2
      *IT_0159 + (-2)*IT_0161;
    const ccomplex_t IT_0163 = IT_0155 + IT_0162;
    const ccomplex_t IT_0164 = IT_0153*IT_0163;
    const ccomplex_t IT_0165 = (-0.25)*IT_0164;
    const ccomplex_t IT_0166 = -IT_0078;
    const ccomplex_t IT_0167 = s_12 + IT_0152 + IT_0165 + IT_0166;
    const ccomplex_t IT_0168 = 8*IT_0130;
    const ccomplex_t IT_0169 = (-8)*IT_0117;
    const ccomplex_t IT_0170 = (-8)*conj(IT_0117);
    const ccomplex_t IT_0171 = 8*conj(IT_0145);
    const ccomplex_t IT_0172 = (-8)*conj(IT_0142);
    const ccomplex_t IT_0173 = IT_0000*IT_0001;
    const ccomplex_t IT_0174 = pow(s_24, 2);
    const ccomplex_t IT_0175 = s_23*s_24*s_34*IT_0077;
    const ccomplex_t IT_0176 = (-0.5)*IT_0173 + 0.5*IT_0174 + -IT_0175;
    const ccomplex_t IT_0177 = pow(s_23, 2);
    const ccomplex_t IT_0178 = IT_0176 + IT_0177;
    const ccomplex_t IT_0179 = IT_0077*IT_0178;
    const ccomplex_t IT_0180 = IT_0077*IT_0174;
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = IT_0001 + IT_0179 + IT_0181;
    const ccomplex_t IT_0183 = m_C_2*IT_0182;
    const ccomplex_t IT_0184 = conj(IT_0145) + conj(IT_0148);
    const ccomplex_t IT_0185 = 4*IT_0076;
    const ccomplex_t IT_0186 = 4*conj(IT_0076);
    const ccomplex_t IT_0187 = (-4)*IT_0076;
    const ccomplex_t IT_0188 = (-4)*conj(IT_0076);
    const ccomplex_t IT_0189 = pow(s_14, 2);
    const ccomplex_t IT_0190 = s_13*s_14*s_34*IT_0077;
    const ccomplex_t IT_0191 = (-0.5)*IT_0173 + 0.5*IT_0189 + -IT_0190;
    const ccomplex_t IT_0192 = pow(s_13, 2);
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = IT_0077*IT_0193;
    const ccomplex_t IT_0195 = IT_0077*IT_0189;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = IT_0001 + IT_0194 + IT_0196;
    const ccomplex_t IT_0198 = m_C_2*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0083;
    const ccomplex_t IT_0200 = s_12 + IT_0199;
    const ccomplex_t IT_0201 = s_14*s_24;
    const ccomplex_t IT_0202 = (-2)*IT_0201;
    const ccomplex_t IT_0203 = IT_0079 + IT_0202;
    const ccomplex_t IT_0204 = IT_0077*IT_0203;
    const ccomplex_t IT_0205 = (-0.5)*IT_0204;
    const ccomplex_t IT_0206 = (-2)*IT_0156 + (-2)*IT_0157 + 4*IT_0158;
    const ccomplex_t IT_0207 = IT_0155 + IT_0206;
    const ccomplex_t IT_0208 = IT_0153*IT_0207;
    const ccomplex_t IT_0209 = 0.25*IT_0208;
    const ccomplex_t IT_0210 = IT_0200 + IT_0205 + IT_0209;
    const ccomplex_t IT_0211 = s_13*s_24*s_34*IT_0077;
    const ccomplex_t IT_0212 = (-2)*IT_0080 + -IT_0201 + 2*IT_0211;
    const ccomplex_t IT_0213 = IT_0079 + IT_0212;
    const ccomplex_t IT_0214 = IT_0077*IT_0213;
    const ccomplex_t IT_0215 = -IT_0078 + (-0.5)*IT_0214;
    const ccomplex_t IT_0216 = s_12 + IT_0215;
    const ccomplex_t IT_0217 = m_C_2*IT_0216;
    const ccomplex_t IT_0218 = (-2)*IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0076;
    const ccomplex_t IT_0220 = (-2)*IT_0160;
    const ccomplex_t IT_0221 = IT_0154 + IT_0220;
    const ccomplex_t IT_0222 = IT_0153*IT_0221;
    const ccomplex_t IT_0223 = 0.25*IT_0222;
    const ccomplex_t IT_0224 = 1 + IT_0223;
    const ccomplex_t IT_0225 = IT_0001*IT_0224;
    const ccomplex_t IT_0226 = 8*IT_0145;
    const ccomplex_t IT_0227 = (-8)*IT_0142;
    const ccomplex_t IT_0228 = -IT_0173;
    const ccomplex_t IT_0229 = IT_0174 + IT_0228;
    const ccomplex_t IT_0230 = IT_0077*IT_0229;
    const ccomplex_t IT_0231 = -IT_0001;
    const ccomplex_t IT_0232 = IT_0180 + IT_0231;
    const ccomplex_t IT_0233 = (-4)*IT_0232;
    const ccomplex_t IT_0234 = IT_0230 + IT_0233;
    const ccomplex_t IT_0235 = m_C_2*IT_0234;
    const ccomplex_t IT_0236 = (-4)*IT_0100*IT_0117 + IT_0098*IT_0198 + (-0.25
      )*IT_0168*IT_0210 + (-0.5)*IT_0076*IT_0218 + (-2)*IT_0139*IT_0225 + (-2)
      *IT_0167*(IT_0148 + 0.125*IT_0226 + 0.125*IT_0227) + (-0.25)*IT_0219
      *IT_0235;
    const ccomplex_t IT_0237 = (-2)*conj(IT_0076);
    const ccomplex_t IT_0238 = 4*conj(IT_0098);
    const ccomplex_t IT_0239 = IT_0021 + IT_0039 + IT_0057 + IT_0075;
    const ccomplex_t IT_0240 = s_23*s_34*IT_0077;
    const ccomplex_t IT_0241 = s_24*IT_0154;
    const ccomplex_t IT_0242 = s_23*s_34*IT_0000;
    const ccomplex_t IT_0243 = (-2)*IT_0242;
    const ccomplex_t IT_0244 = IT_0241 + IT_0243;
    const ccomplex_t IT_0245 = IT_0153*IT_0244;
    const ccomplex_t IT_0246 = (-0.666666666666667)*IT_0240 + (
      -0.166666666666667)*IT_0245;
    const ccomplex_t IT_0247 = s_24 + IT_0246;
    const ccomplex_t IT_0248 = m_C_2*IT_0247;
    const ccomplex_t IT_0249 = (-6)*IT_0248;
    const ccomplex_t IT_0250 = 2*IT_0145;
    const ccomplex_t IT_0251 = 2*IT_0148;
    const ccomplex_t IT_0252 = (-2)*IT_0142;
    const ccomplex_t IT_0253 = IT_0250 + IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = m_C_2*IT_0077;
    const ccomplex_t IT_0255 = s_14*IT_0000;
    const ccomplex_t IT_0256 = s_13*s_34;
    const ccomplex_t IT_0257 = s_14*IT_0154;
    const ccomplex_t IT_0258 = IT_0077*IT_0257;
    const ccomplex_t IT_0259 = (-2)*IT_0256 + (-0.5)*IT_0258;
    const ccomplex_t IT_0260 = IT_0255 + IT_0259;
    const ccomplex_t IT_0261 = IT_0254*IT_0260;
    const ccomplex_t IT_0262 = (-2)*IT_0261;
    const ccomplex_t IT_0263 = 2*IT_0139;
    const ccomplex_t IT_0264 = s_12*s_14*IT_0000;
    const ccomplex_t IT_0265 = s_24*IT_0192;
    const ccomplex_t IT_0266 = s_12*s_13*s_34;
    const ccomplex_t IT_0267 = s_24*IT_0000*IT_0001;
    const ccomplex_t IT_0268 = 2*IT_0265 + (-2)*IT_0266 + -IT_0267;
    const ccomplex_t IT_0269 = IT_0264 + IT_0268;
    const ccomplex_t IT_0270 = IT_0077*IT_0269;
    const ccomplex_t IT_0271 = s_24*IT_0001;
    const ccomplex_t IT_0272 = 2*IT_0271;
    const ccomplex_t IT_0273 = s_12*s_14;
    const ccomplex_t IT_0274 = (-2)*IT_0273;
    const ccomplex_t IT_0275 = IT_0270 + IT_0272 + IT_0274;
    const ccomplex_t IT_0276 = 2*IT_0098;
    const ccomplex_t IT_0277 = s_13*s_34*IT_0000;
    const ccomplex_t IT_0278 = (-2)*IT_0277;
    const ccomplex_t IT_0279 = IT_0257 + IT_0278;
    const ccomplex_t IT_0280 = IT_0153*IT_0279;
    const ccomplex_t IT_0281 = (-2)*IT_0256;
    const ccomplex_t IT_0282 = IT_0255 + IT_0281;
    const ccomplex_t IT_0283 = IT_0077*IT_0282;
    const ccomplex_t IT_0284 = 0.166666666666667*IT_0280 + (-0.666666666666667
      )*IT_0283;
    const ccomplex_t IT_0285 = s_14 + IT_0284;
    const ccomplex_t IT_0286 = m_C_2*IT_0285;
    const ccomplex_t IT_0287 = 6*IT_0286;
    const ccomplex_t IT_0288 = 2*IT_0117;
    const ccomplex_t IT_0289 = s_12*s_24*IT_0000;
    const ccomplex_t IT_0290 = s_13*s_23*s_24;
    const ccomplex_t IT_0291 = s_14*IT_0000*IT_0001;
    const ccomplex_t IT_0292 = s_13*s_34*IT_0001;
    const ccomplex_t IT_0293 = (-2)*IT_0290 + -IT_0291 + 2*IT_0292;
    const ccomplex_t IT_0294 = IT_0289 + IT_0293;
    const ccomplex_t IT_0295 = IT_0077*IT_0294;
    const ccomplex_t IT_0296 = s_14*IT_0001;
    const ccomplex_t IT_0297 = 2*IT_0296;
    const ccomplex_t IT_0298 = s_12*s_24;
    const ccomplex_t IT_0299 = (-2)*IT_0298;
    const ccomplex_t IT_0300 = IT_0295 + IT_0297 + IT_0299;
    const ccomplex_t IT_0301 = IT_0153*IT_0241;
    const ccomplex_t IT_0302 = (-10)*s_24;
    const ccomplex_t IT_0303 = IT_0301 + IT_0302;
    const ccomplex_t IT_0304 = m_C_2*IT_0303;
    const ccomplex_t IT_0305 = (-2)*IT_0130;
    const ccomplex_t IT_0306 = IT_0189 + IT_0228;
    const ccomplex_t IT_0307 = IT_0077*IT_0306;
    const ccomplex_t IT_0308 = IT_0195 + IT_0231;
    const ccomplex_t IT_0309 = (-4)*IT_0308;
    const ccomplex_t IT_0310 = IT_0307 + IT_0309;
    const ccomplex_t IT_0311 = m_C_2*IT_0310;
    const ccomplex_t IT_0312 = (-2)*conj(IT_0142);
    const ccomplex_t IT_0313 = 9*IT_0173;
    const ccomplex_t IT_0314 = (-2)*IT_0239;
    const ccomplex_t IT_0315 = IT_0001*IT_0077;
    const ccomplex_t IT_0316 = (-0.25)*IT_0154;
    const ccomplex_t IT_0317 = IT_0160 + IT_0316;
    const ccomplex_t IT_0318 = IT_0315*IT_0317;
    const ccomplex_t IT_0319 = (-4)*IT_0318;
    const ccomplex_t IT_0320 = -IT_0088 + -IT_0091 + -IT_0094 + -IT_0097;
    const ccomplex_t IT_0321 = (-2)*IT_0320;
    const ccomplex_t IT_0322 = 0.125*IT_0184;
    const ccomplex_t IT_0323 = (-2)*IT_0156 + (-2)*IT_0157 + 4*IT_0159;
    const ccomplex_t IT_0324 = IT_0155 + IT_0323;
    const ccomplex_t IT_0325 = IT_0153*IT_0324;
    const ccomplex_t IT_0326 = 0.25*IT_0325;
    const ccomplex_t IT_0327 = IT_0200 + IT_0205 + IT_0326;
    const ccomplex_t IT_0328 = IT_0001*IT_0308;
    const ccomplex_t IT_0329 = (-2)*IT_0098;
    const ccomplex_t IT_0330 = -(s_12 + (-0.5)*IT_0083)*(IT_0001 + -IT_0195);
    const ccomplex_t IT_0331 = (-2)*IT_0330;
    const ccomplex_t IT_0332 = 0.166666666666667*IT_0331;
    const ccomplex_t IT_0333 = s_12 + IT_0166;
    const ccomplex_t IT_0334 = IT_0001*IT_0333;
    const ccomplex_t IT_0335 = s_13*s_14*s_23;
    const ccomplex_t IT_0336 = s_23*s_34*IT_0001;
    const ccomplex_t IT_0337 = -IT_0267 + (-2)*IT_0335 + 2*IT_0336;
    const ccomplex_t IT_0338 = IT_0264 + IT_0337;
    const ccomplex_t IT_0339 = IT_0077*IT_0338;
    const ccomplex_t IT_0340 = (-0.5)*IT_0339;
    const ccomplex_t IT_0341 = -IT_0271;
    const ccomplex_t IT_0342 = IT_0273 + IT_0340 + IT_0341;
    const ccomplex_t IT_0343 = s_14*IT_0177;
    const ccomplex_t IT_0344 = s_12*s_23*s_34;
    const ccomplex_t IT_0345 = -IT_0291 + 2*IT_0343 + (-2)*IT_0344;
    const ccomplex_t IT_0346 = IT_0289 + IT_0345;
    const ccomplex_t IT_0347 = IT_0077*IT_0346;
    const ccomplex_t IT_0348 = 0.5*IT_0347;
    const ccomplex_t IT_0349 = -IT_0348;
    const ccomplex_t IT_0350 = -IT_0296;
    const ccomplex_t IT_0351 = IT_0298 + IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = s_13*s_34*IT_0077;
    const ccomplex_t IT_0353 = (-0.166666666666667)*IT_0280 + (
      -0.666666666666667)*IT_0352;
    const ccomplex_t IT_0354 = s_14 + IT_0353;
    const ccomplex_t IT_0355 = m_C_2*IT_0354;
    const ccomplex_t IT_0356 = (-6)*IT_0355;
    const ccomplex_t IT_0357 = 0.25*IT_0356;
    const ccomplex_t IT_0358 = IT_0077*IT_0241;
    const ccomplex_t IT_0359 = s_24*IT_0000;
    const ccomplex_t IT_0360 = s_23*s_34;
    const ccomplex_t IT_0361 = (-2)*IT_0359 + 4*IT_0360;
    const ccomplex_t IT_0362 = IT_0358 + IT_0361;
    const ccomplex_t IT_0363 = IT_0254*IT_0362;
    const ccomplex_t IT_0364 = 0.25*IT_0363;
    const ccomplex_t IT_0365 = (-2)*IT_0360;
    const ccomplex_t IT_0366 = IT_0359 + IT_0365;
    const ccomplex_t IT_0367 = IT_0077*IT_0366;
    const ccomplex_t IT_0368 = 0.166666666666667*IT_0245 + (-0.666666666666667
      )*IT_0367;
    const ccomplex_t IT_0369 = s_24 + IT_0368;
    const ccomplex_t IT_0370 = m_C_2*IT_0369;
    const ccomplex_t IT_0371 = 6*IT_0370;
    const ccomplex_t IT_0372 = 0.25*IT_0371;
    const ccomplex_t IT_0373 = IT_0153*IT_0257;
    const ccomplex_t IT_0374 = (-10)*s_14;
    const ccomplex_t IT_0375 = IT_0373 + IT_0374;
    const ccomplex_t IT_0376 = m_C_2*IT_0375;
    const ccomplex_t IT_0377 = 0.25*IT_0376;
    const ccomplex_t IT_0378 = 4*conj(IT_0320)*(0.25*IT_0314*IT_0319 + 0.25
      *IT_0313*IT_0321 + IT_0098*IT_0342 + -IT_0076*IT_0351 + IT_0305*IT_0357 +
       IT_0288*IT_0364 + IT_0263*IT_0372 + IT_0253*IT_0377);
    const ccomplex_t IT_0379 = -IT_0201;
    const ccomplex_t IT_0380 = IT_0079 + IT_0379;
    const ccomplex_t IT_0381 = IT_0077*IT_0380;
    const ccomplex_t IT_0382 = -IT_0078 + (-0.25)*IT_0381;
    const ccomplex_t IT_0383 = s_12 + IT_0382;
    const ccomplex_t IT_0384 = m_C_2*IT_0383;
    const ccomplex_t IT_0385 = (-4)*IT_0384;
    const ccomplex_t IT_0386 = (-0.5)*IT_0149;
    const ccomplex_t IT_0387 = (-0.5)*IT_0184;
    const ccomplex_t IT_0388 = conj(IT_0098)*IT_0117 + IT_0098*conj(IT_0117) +
       (-0.5)*IT_0099*IT_0130 + -conj(IT_0076)*IT_0142 + -IT_0076*conj(IT_0142) 
      + (-0.5)*conj(IT_0139)*IT_0219 + (-0.5)*IT_0139*IT_0237 + (-0.5)*conj
      (IT_0130)*IT_0276 + IT_0237*IT_0386 + IT_0219*IT_0387;
    const ccomplex_t IT_0389 = s_14*s_23*s_34*IT_0077;
    const ccomplex_t IT_0390 = (-2)*IT_0080 + -IT_0201 + 2*IT_0389;
    const ccomplex_t IT_0391 = IT_0079 + IT_0390;
    const ccomplex_t IT_0392 = IT_0077*IT_0391;
    const ccomplex_t IT_0393 = -IT_0078 + (-0.5)*IT_0392;
    const ccomplex_t IT_0394 = s_12 + IT_0393;
    const ccomplex_t IT_0395 = m_C_2*IT_0394;
    const ccomplex_t IT_0396 = (-2)*IT_0395;
    const ccomplex_t IT_0397 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0398 = pow(m_W, 6);
    const ccomplex_t IT_0399 = s_12*IT_0398;
    const ccomplex_t IT_0400 = s_13*s_23*IT_0154;
    const ccomplex_t IT_0401 = (-2)*IT_0400;
    const ccomplex_t IT_0402 = IT_0399 + IT_0401;
    const ccomplex_t IT_0403 = IT_0153*IT_0402;
    const ccomplex_t IT_0404 = 0.5*IT_0403;
    const ccomplex_t IT_0405 = (-4)*IT_0201;
    const ccomplex_t IT_0406 = IT_0079 + IT_0404 + IT_0405;
    const ccomplex_t IT_0407 = IT_0077*IT_0324;
    const ccomplex_t IT_0408 = -IT_0407;
    const ccomplex_t IT_0409 = IT_0406 + IT_0408;
    const ccomplex_t IT_0410 = 0.25*IT_0403;
    const ccomplex_t IT_0411 = IT_0079 + IT_0410;
    const ccomplex_t IT_0412 = (-2)*conj(IT_0130);
    const ccomplex_t IT_0413 = (-2)*conj(IT_0239);
    const ccomplex_t IT_0414 = (-2)*conj(IT_0320);
    const ccomplex_t IT_0415 = 2*conj(IT_0117);
    const ccomplex_t IT_0416 = (-0.25)*IT_0363;
    const ccomplex_t IT_0417 = 2*conj(IT_0139);
    const ccomplex_t IT_0418 = (-0.25)*IT_0371;
    const ccomplex_t IT_0419 = 2*conj(IT_0145);
    const ccomplex_t IT_0420 = 2*conj(IT_0148);
    const ccomplex_t IT_0421 = IT_0312 + IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = (-0.25)*IT_0376;
    const ccomplex_t IT_0423 = -conj(IT_0098)*IT_0342 + conj(IT_0076)*IT_0351 
      + conj(IT_0320)*IT_0409 + 2*conj(IT_0239)*IT_0411 + (-0.25)*IT_0356
      *IT_0412 + (-0.25)*IT_0319*IT_0413 + (-0.25)*IT_0313*IT_0414 + IT_0415
      *IT_0416 + IT_0417*IT_0418 + IT_0421*IT_0422;
    const ccomplex_t IT_0424 = (-4)*IT_0320*((-0.25)*IT_0238*IT_0342 + (-0.25)
      *IT_0188*IT_0351 + conj(IT_0320)*IT_0409 + 2*conj(IT_0239)*IT_0411 + (
      -0.25)*IT_0356*IT_0415 + IT_0412*IT_0416 + IT_0418*IT_0421 + IT_0417
      *IT_0422);
    const ccomplex_t IT_0425 = (-0.166666666666667)*IT_0085;
    const ccomplex_t IT_0426 = (s_12 + (-0.5)*IT_0083)*(IT_0001 + -IT_0180);
    const ccomplex_t IT_0427 = 2*IT_0426;
    const ccomplex_t IT_0428 = (-0.166666666666667)*IT_0427;
    const ccomplex_t IT_0429 = IT_0001 + IT_0181;
    const ccomplex_t IT_0430 = IT_0001*IT_0429;
    const ccomplex_t IT_0431 = -IT_0430;
    const ccomplex_t IT_0432 = (-0.166666666666667)*IT_0235;
    const ccomplex_t IT_0433 = IT_0077*IT_0207;
    const ccomplex_t IT_0434 = -IT_0433;
    const ccomplex_t IT_0435 = IT_0406 + IT_0434;
    const ccomplex_t IT_0436 = (-0.125)*IT_0304;
    const ccomplex_t IT_0437 = (-0.125)*IT_0249;
    const ccomplex_t IT_0438 = (-0.125)*IT_0262;
    const ccomplex_t IT_0439 = (-0.125)*IT_0275;
    const ccomplex_t IT_0440 = (-0.125)*IT_0300;
    const ccomplex_t IT_0441 = conj(IT_0320)*IT_0411 + (-0.125)*IT_0287
      *IT_0412 + (-0.125)*IT_0313*IT_0413 + (-0.125)*IT_0319*IT_0414 + 0.5*conj
      (IT_0239)*IT_0435 + IT_0415*IT_0436 + IT_0417*IT_0437 + IT_0421*IT_0438 +
       IT_0099*IT_0439 + IT_0237*IT_0440;
    const ccomplex_t IT_0442 = IT_0076*IT_0085*IT_0099 + IT_0100*(16*conj
      (IT_0117)*IT_0130 + 16*conj(IT_0139)*IT_0142 + IT_0139*(16*conj(IT_0142) +
       (-16)*conj(IT_0145) + (-16)*conj(IT_0148)) + IT_0149*IT_0150) + IT_0167*
      (conj(IT_0148)*IT_0168 + conj(IT_0139)*IT_0169 + IT_0139*IT_0170 + IT_0130
      *(IT_0171 + IT_0172)) + 4*IT_0183*(conj(IT_0076)*IT_0139 + IT_0076*conj
      (IT_0139) + 0.25*IT_0184*IT_0185 + 0.25*IT_0149*IT_0186 + 0.25*conj
      (IT_0142)*IT_0187 + 0.25*IT_0142*IT_0188) + 4*conj(IT_0117)*(IT_0098
      *IT_0198 + 2*IT_0117*IT_0210 + 0.25*IT_0218*IT_0219) + (-4)*conj(IT_0130)
      *IT_0236 + IT_0130*IT_0235*IT_0237 + IT_0117*(IT_0218*IT_0237 + IT_0198
      *IT_0238) + conj(IT_0239)*(IT_0249*IT_0253 + IT_0262*IT_0263 + IT_0275
      *IT_0276 + IT_0287*IT_0288 + IT_0219*IT_0300 + IT_0304*IT_0305) + IT_0311*
      (IT_0099*IT_0139 + IT_0099*IT_0149 + conj(IT_0098)*IT_0252 + conj(IT_0139)
      *IT_0276 + IT_0184*IT_0276 + IT_0098*IT_0312) + conj(IT_0239)*(IT_0253
      *IT_0262 + IT_0249*IT_0263 + IT_0275*IT_0276 + IT_0219*IT_0300 + IT_0288
      *IT_0304 + IT_0287*IT_0305 + IT_0313*IT_0314 + IT_0319*IT_0321) + 8
      *IT_0225*(IT_0130*conj(IT_0139) + conj(IT_0117)*IT_0142 + IT_0117*conj
      (IT_0142) + 0.125*IT_0149*IT_0170 + IT_0169*IT_0322) + 8*(IT_0139*conj
      (IT_0139) + IT_0142*conj(IT_0142) + 0.125*IT_0145*IT_0172 + IT_0148*(conj
      (IT_0148) + 0.125*IT_0171 + 0.125*IT_0172) + (IT_0226 + IT_0227)*IT_0322)
      *IT_0327 + 6*conj(IT_0098)*(IT_0098*IT_0328 + IT_0329*IT_0332) + 6*conj
      (IT_0098)*((-0.666666666666667)*IT_0130*IT_0198 + IT_0098*IT_0328 +
       IT_0329*IT_0332 + -IT_0076*IT_0334) + 4*conj(IT_0320)*(IT_0098*IT_0342 +
       0.25*IT_0187*IT_0351 + IT_0288*IT_0357 + IT_0305*IT_0364 + IT_0253
      *IT_0372 + IT_0263*IT_0377) + IT_0378 + (-2)*IT_0385*IT_0388 + (-2)
      *IT_0396*(conj(IT_0098)*IT_0139 + (-0.5)*IT_0099*IT_0142 + (-0.5)*conj
      (IT_0142)*IT_0276 + (-0.5)*conj(IT_0139)*IT_0329 + IT_0329*IT_0387 +
       IT_0386*IT_0397) + (-4)*IT_0320*IT_0423 + IT_0424 + (-6)*conj(IT_0076)*
      (IT_0098*IT_0334 + IT_0276*IT_0425 + IT_0219*IT_0428) + (-6)*conj(IT_0076)
      *((-0.333333333333333)*IT_0130*IT_0218 + IT_0098*IT_0334 + IT_0276*IT_0425
       + IT_0219*IT_0428 + -IT_0076*IT_0431 + IT_0288*IT_0432) + (-6)*IT_0076*
      (conj(IT_0098)*IT_0334 + IT_0099*IT_0425 + -conj(IT_0076)*IT_0431 +
       IT_0415*IT_0432) + (-8)*IT_0239*(conj(IT_0320)*IT_0411 + (-0.125)*IT_0287
      *IT_0415 + 0.5*conj(IT_0239)*IT_0435 + IT_0412*IT_0436 + IT_0421*IT_0437 +
       IT_0417*IT_0438 + IT_0099*IT_0439 + IT_0237*IT_0440) + (-8)*IT_0239
      *IT_0441;
    return create_ccomplex_return(IT_0442);
}

