#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_2_to_Z_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_2_to_Z_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t beta = param->beta;
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
    const creal_t Gamma_W = param->Gamma_W;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
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
    const ccomplex_t IT_0000 = m_C_2*m_N_2;
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(s_34, 2);
    const ccomplex_t IT_0004 = IT_0001*IT_0002*IT_0003;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = 1 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = cos(theta_W);
    const ccomplex_t IT_0009 = sin(theta_W);
    const ccomplex_t IT_0010 = cpow(IT_0009, -1);
    const ccomplex_t IT_0011 = (0 + _Complex_I*1)*e_em*IT_0008*IT_0010;
    const ccomplex_t IT_0012 = m_C_2*IT_0011;
    const ccomplex_t IT_0013 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0016 = IT_0010*IT_0015;
    const ccomplex_t IT_0017 = 1.4142135623731*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0014 + (-0.5)*IT_0017);
    const ccomplex_t IT_0019 = -IT_0018;
    const ccomplex_t IT_0020 = pow(m_C_2, 2);
    const ccomplex_t IT_0021 = pow(m_N_2, 2);
    const ccomplex_t IT_0022 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0020 + -IT_0021 + -reg_prop, -1);
    const ccomplex_t IT_0023 = IT_0019*IT_0022;
    const ccomplex_t IT_0024 = (0 + _Complex_I*1)*IT_0023;
    const ccomplex_t IT_0025 = IT_0012*IT_0024;
    const ccomplex_t IT_0026 = cos(beta);
    const ccomplex_t IT_0027 = IT_0010*IT_0026;
    const ccomplex_t IT_0028 = conj(N_d2)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0029 = IT_0027*IT_0028;
    const ccomplex_t IT_0030 = (-2)*IT_0029;
    const ccomplex_t IT_0031 = cpow(IT_0008, -1);
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = conj(N_B2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0034 = IT_0032*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W2)*conj(U_d2)*e_em;
    const ccomplex_t IT_0037 = IT_0027*IT_0036;
    const ccomplex_t IT_0038 = 1.4142135623731*IT_0037;
    const ccomplex_t IT_0039 = IT_0035 + IT_0038;
    const ccomplex_t IT_0040 = IT_0030 + IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*m_W*e_em*IT_0009*IT_0031;
    const ccomplex_t IT_0044 = -IT_0043;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = pow(m_W, 2);
    const ccomplex_t IT_0047 = cpow((-2)*s_12 + -IT_0020 + -IT_0021 + IT_0046 
      + -reg_prop, -1);
    const ccomplex_t IT_0048 = IT_0045*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*IT_0048;
    const ccomplex_t IT_0050 = -IT_0025 + IT_0049;
    const ccomplex_t IT_0051 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0052 = IT_0010*IT_0051;
    const ccomplex_t IT_0053 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0054 = IT_0010*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + 0.5*IT_0055);
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0022*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0012*IT_0059;
    const ccomplex_t IT_0061 = m_N_2*IT_0011;
    const ccomplex_t IT_0062 = IT_0024*IT_0061;
    const ccomplex_t IT_0063 = sin(beta);
    const ccomplex_t IT_0064 = IT_0010*IT_0063;
    const ccomplex_t IT_0065 = N_u2*e_em*V_Wp2;
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = 2*IT_0066;
    const ccomplex_t IT_0068 = IT_0031*IT_0063;
    const ccomplex_t IT_0069 = N_B2*V_u2*e_em;
    const ccomplex_t IT_0070 = IT_0068*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = N_W2*V_u2*e_em;
    const ccomplex_t IT_0073 = IT_0064*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = IT_0071 + IT_0074;
    const ccomplex_t IT_0076 = IT_0067 + IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0077;
    const ccomplex_t IT_0079 = IT_0044*IT_0078;
    const ccomplex_t IT_0080 = IT_0047*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = IT_0062 + -IT_0081;
    const ccomplex_t IT_0083 = (-4)*IT_0082;
    const ccomplex_t IT_0084 = 0.25*IT_0083;
    const ccomplex_t IT_0085 = IT_0060 + IT_0084;
    const ccomplex_t IT_0086 = (-4)*conj(IT_0082);
    const ccomplex_t IT_0087 = 0.25*IT_0086;
    const ccomplex_t IT_0088 = conj(IT_0060) + IT_0087;
    const ccomplex_t IT_0089 = IT_0059*IT_0061;
    const ccomplex_t IT_0090 = -IT_0089;
    const ccomplex_t IT_0091 = (-4)*IT_0060;
    const ccomplex_t IT_0092 = 0.25*IT_0091;
    const ccomplex_t IT_0093 = IT_0082 + IT_0092;
    const ccomplex_t IT_0094 = (-4)*conj(IT_0060);
    const ccomplex_t IT_0095 = 0.25*IT_0094;
    const ccomplex_t IT_0096 = conj(IT_0082) + IT_0095;
    const ccomplex_t IT_0097 = s_12*s_23;
    const ccomplex_t IT_0098 = s_13*IT_0021;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0097 + IT_0099;
    const ccomplex_t IT_0101 = s_12*s_23*IT_0046;
    const ccomplex_t IT_0102 = s_14*s_23*s_24;
    const ccomplex_t IT_0103 = s_13*IT_0021*IT_0046;
    const ccomplex_t IT_0104 = s_14*s_34*IT_0021;
    const ccomplex_t IT_0105 = (-2)*IT_0102 + -IT_0103 + 2*IT_0104;
    const ccomplex_t IT_0106 = IT_0101 + IT_0105;
    const ccomplex_t IT_0107 = IT_0001*IT_0106;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0100 + IT_0108;
    const ccomplex_t IT_0110 = IT_0009*IT_0031;
    const ccomplex_t IT_0111 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = IT_0008*IT_0010;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0116 = IT_0110*IT_0115;
    const ccomplex_t IT_0117 = IT_0113*IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0112 + IT_0114 + 
      -IT_0116 + -IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0121 = IT_0010*IT_0120;
    const ccomplex_t IT_0122 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0123 = IT_0010*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = pow(m_Z, 2);
    const ccomplex_t IT_0128 = cpow((-2)*s_23 + IT_0021 + IT_0127 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0129 = IT_0126*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0119*IT_0130;
    const ccomplex_t IT_0132 = (-2)*IT_0011;
    const ccomplex_t IT_0133 = IT_0059*IT_0132;
    const ccomplex_t IT_0134 = N_u2*conj(N_u2)*e_em;
    const ccomplex_t IT_0135 = IT_0110*IT_0134;
    const ccomplex_t IT_0136 = IT_0113*IT_0134;
    const ccomplex_t IT_0137 = N_d2*conj(N_d2)*e_em;
    const ccomplex_t IT_0138 = IT_0110*IT_0137;
    const ccomplex_t IT_0139 = IT_0113*IT_0137;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0135 + IT_0136 + 
      -IT_0138 + -IT_0139);
    const ccomplex_t IT_0141 = (-0.5)*IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_23 + IT_0127 + (0 + _Complex_I*1)
      *m_N_2*Gamma_o2 + reg_prop, -1);
    const ccomplex_t IT_0143 = IT_0057*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = IT_0141*IT_0144;
    const ccomplex_t IT_0146 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0147 = IT_0110*IT_0146;
    const ccomplex_t IT_0148 = IT_0113*IT_0146;
    const ccomplex_t IT_0149 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0150 = IT_0110*IT_0149;
    const ccomplex_t IT_0151 = IT_0113*IT_0149;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0147 + IT_0148 + 
      -IT_0150 + -IT_0151);
    const ccomplex_t IT_0153 = (-0.5)*IT_0152;
    const ccomplex_t IT_0154 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0155 = IT_0010*IT_0154;
    const ccomplex_t IT_0156 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0157 = IT_0010*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + 0.5*IT_0158);
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = cpow((-2)*s_23 + IT_0021 + IT_0127 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0162 = IT_0160*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = IT_0153*IT_0163;
    const ccomplex_t IT_0165 = (-2)*IT_0131 + -IT_0133 + -IT_0145 + (-2)
      *IT_0164;
    const ccomplex_t IT_0166 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0167 = IT_0110*IT_0166;
    const ccomplex_t IT_0168 = IT_0113*IT_0166;
    const ccomplex_t IT_0169 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0170 = IT_0110*IT_0169;
    const ccomplex_t IT_0171 = IT_0113*IT_0169;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*(IT_0167 + IT_0168 + 
      -IT_0170 + -IT_0171);
    const ccomplex_t IT_0173 = (-0.5)*IT_0172;
    const ccomplex_t IT_0174 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0175 = IT_0010*IT_0174;
    const ccomplex_t IT_0176 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0177 = IT_0010*IT_0176;
    const ccomplex_t IT_0178 = 1.4142135623731*IT_0177;
    const ccomplex_t IT_0179 = (0 + _Complex_I*1)*(IT_0175 + 0.5*IT_0178);
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = cpow((-2)*s_23 + IT_0021 + IT_0127 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0182 = IT_0180*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*IT_0182;
    const ccomplex_t IT_0184 = IT_0173*IT_0183;
    const ccomplex_t IT_0185 = 0.5*IT_0140;
    const ccomplex_t IT_0186 = IT_0144*IT_0185;
    const ccomplex_t IT_0187 = 2*IT_0184 + IT_0186;
    const ccomplex_t IT_0188 = conj(IT_0165) + conj(IT_0187);
    const ccomplex_t IT_0189 = -IT_0184 + (-0.5)*IT_0186;
    const ccomplex_t IT_0190 = 0.5*IT_0145;
    const ccomplex_t IT_0191 = IT_0131 + IT_0164;
    const ccomplex_t IT_0192 = IT_0190 + IT_0191;
    const ccomplex_t IT_0193 = IT_0189 + IT_0192;
    const ccomplex_t IT_0194 = IT_0165 + IT_0187;
    const ccomplex_t IT_0195 = conj(IT_0190) + conj(IT_0191);
    const ccomplex_t IT_0196 = conj(IT_0189) + IT_0195;
    const ccomplex_t IT_0197 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0198 = IT_0110*IT_0197;
    const ccomplex_t IT_0199 = IT_0113*IT_0197;
    const ccomplex_t IT_0200 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0201 = IT_0110*IT_0200;
    const ccomplex_t IT_0202 = IT_0113*IT_0200;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*(IT_0198 + IT_0199 + 
      -IT_0201 + -IT_0202);
    const ccomplex_t IT_0204 = (-0.5)*IT_0203;
    const ccomplex_t IT_0205 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0206 = IT_0010*IT_0205;
    const ccomplex_t IT_0207 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0208 = IT_0010*IT_0207;
    const ccomplex_t IT_0209 = 1.4142135623731*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*(IT_0206 + (-0.5)*IT_0209);
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0128*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = IT_0204*IT_0213;
    const ccomplex_t IT_0215 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0216 = IT_0110*IT_0215;
    const ccomplex_t IT_0217 = IT_0113*IT_0215;
    const ccomplex_t IT_0218 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0219 = IT_0110*IT_0218;
    const ccomplex_t IT_0220 = IT_0113*IT_0218;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*(IT_0216 + IT_0217 + 
      -IT_0219 + -IT_0220);
    const ccomplex_t IT_0222 = (-0.5)*IT_0221;
    const ccomplex_t IT_0223 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0224 = IT_0010*IT_0223;
    const ccomplex_t IT_0225 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0226 = IT_0010*IT_0225;
    const ccomplex_t IT_0227 = 1.4142135623731*IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*(IT_0224 + (-0.5)*IT_0227);
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = IT_0161*IT_0229;
    const ccomplex_t IT_0231 = (0 + _Complex_I*1)*IT_0230;
    const ccomplex_t IT_0232 = IT_0222*IT_0231;
    const ccomplex_t IT_0233 = (-2)*IT_0214 + (-2)*IT_0232;
    const ccomplex_t IT_0234 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0235 = IT_0110*IT_0234;
    const ccomplex_t IT_0236 = IT_0113*IT_0234;
    const ccomplex_t IT_0237 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0238 = IT_0110*IT_0237;
    const ccomplex_t IT_0239 = IT_0113*IT_0237;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*(IT_0235 + IT_0236 + 
      -IT_0238 + -IT_0239);
    const ccomplex_t IT_0241 = (-0.5)*IT_0240;
    const ccomplex_t IT_0242 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0243 = IT_0010*IT_0242;
    const ccomplex_t IT_0244 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0245 = IT_0010*IT_0244;
    const ccomplex_t IT_0246 = 1.4142135623731*IT_0245;
    const ccomplex_t IT_0247 = (0 + _Complex_I*1)*(IT_0243 + (-0.5)*IT_0246);
    const ccomplex_t IT_0248 = -IT_0247;
    const ccomplex_t IT_0249 = IT_0181*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = IT_0241*IT_0250;
    const ccomplex_t IT_0252 = IT_0024*IT_0132;
    const ccomplex_t IT_0253 = IT_0019*IT_0142;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = IT_0185*IT_0254;
    const ccomplex_t IT_0256 = 2*IT_0251 + -IT_0252 + (-2)*IT_0255;
    const ccomplex_t IT_0257 = conj(IT_0233) + conj(IT_0256);
    const ccomplex_t IT_0258 = IT_0214 + IT_0232;
    const ccomplex_t IT_0259 = -IT_0251;
    const ccomplex_t IT_0260 = IT_0255 + IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = IT_0233 + IT_0256;
    const ccomplex_t IT_0262 = conj(IT_0255) + conj(IT_0258) + conj(IT_0259);
    const ccomplex_t IT_0263 = IT_0001*IT_0002;
    const ccomplex_t IT_0264 = pow(m_Z, 4);
    const ccomplex_t IT_0265 = IT_0046*IT_0264;
    const ccomplex_t IT_0266 = IT_0263*IT_0265;
    const ccomplex_t IT_0267 = IT_0046*IT_0127;
    const ccomplex_t IT_0268 = IT_0001*IT_0267;
    const ccomplex_t IT_0269 = 12*IT_0127 + (-4)*IT_0268;
    const ccomplex_t IT_0270 = IT_0266 + IT_0269;
    const ccomplex_t IT_0271 = IT_0000*IT_0270;
    const ccomplex_t IT_0272 = IT_0255 + IT_0259;
    const ccomplex_t IT_0273 = conj(IT_0255) + conj(IT_0259);
    const ccomplex_t IT_0274 = pow(s_23, 2);
    const ccomplex_t IT_0275 = IT_0002*IT_0274;
    const ccomplex_t IT_0276 = pow(s_24, 2);
    const ccomplex_t IT_0277 = IT_0001*IT_0276;
    const ccomplex_t IT_0278 = IT_0275 + IT_0277;
    const ccomplex_t IT_0279 = s_23*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0280 = -IT_0021 + -IT_0279;
    const ccomplex_t IT_0281 = IT_0278 + IT_0280;
    const ccomplex_t IT_0282 = m_C_2*IT_0281;
    const ccomplex_t IT_0283 = IT_0011*IT_0059;
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = 2*IT_0060;
    const ccomplex_t IT_0286 = (-2)*IT_0082;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = 2*conj(IT_0060);
    const ccomplex_t IT_0289 = (-2)*conj(IT_0082);
    const ccomplex_t IT_0290 = IT_0288 + IT_0289;
    const ccomplex_t IT_0291 = 2*IT_0090;
    const ccomplex_t IT_0292 = (-2)*IT_0050;
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = IT_0011*IT_0024;
    const ccomplex_t IT_0295 = -IT_0294;
    const ccomplex_t IT_0296 = 2*conj(IT_0090);
    const ccomplex_t IT_0297 = (-2)*conj(IT_0050);
    const ccomplex_t IT_0298 = IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = s_14*s_23*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0300 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0301 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0302 = -IT_0300 + -IT_0301;
    const ccomplex_t IT_0303 = s_12 + IT_0299 + IT_0302;
    const ccomplex_t IT_0304 = m_C_2*IT_0303;
    const ccomplex_t IT_0305 = 2*conj(IT_0082);
    const ccomplex_t IT_0306 = (-2)*conj(IT_0060);
    const ccomplex_t IT_0307 = IT_0305 + IT_0306;
    const ccomplex_t IT_0308 = s_13*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0309 = s_12 + IT_0302 + IT_0308;
    const ccomplex_t IT_0310 = m_N_2*IT_0309;
    const ccomplex_t IT_0311 = s_23*IT_0127;
    const ccomplex_t IT_0312 = IT_0002*IT_0311;
    const ccomplex_t IT_0313 = s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0314 = s_24*IT_0127;
    const ccomplex_t IT_0315 = IT_0313*IT_0314;
    const ccomplex_t IT_0316 = s_23*IT_0046;
    const ccomplex_t IT_0317 = s_24*s_34;
    const ccomplex_t IT_0318 = (-2)*IT_0317;
    const ccomplex_t IT_0319 = IT_0316 + IT_0318;
    const ccomplex_t IT_0320 = IT_0001*IT_0319;
    const ccomplex_t IT_0321 = 0.5*IT_0312 + (-0.5)*IT_0315 + (-0.5)*IT_0320;
    const ccomplex_t IT_0322 = s_23 + IT_0321;
    const ccomplex_t IT_0323 = m_C_2*IT_0322;
    const ccomplex_t IT_0324 = 2*IT_0323;
    const ccomplex_t IT_0325 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0326 = IT_0010*IT_0325;
    const ccomplex_t IT_0327 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0328 = IT_0010*IT_0327;
    const ccomplex_t IT_0329 = 1.4142135623731*IT_0328;
    const ccomplex_t IT_0330 = (0 + _Complex_I*1)*(IT_0326 + 0.5*IT_0329);
    const ccomplex_t IT_0331 = -IT_0330;
    const ccomplex_t IT_0332 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0333 = IT_0110*IT_0332;
    const ccomplex_t IT_0334 = IT_0113*IT_0332;
    const ccomplex_t IT_0335 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0336 = IT_0113*IT_0335;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*(IT_0333 + -IT_0334 + (-2)
      *IT_0336);
    const ccomplex_t IT_0338 = 0.5*IT_0337;
    const ccomplex_t IT_0339 = IT_0331*IT_0338;
    const ccomplex_t IT_0340 = cpow((-2)*s_13 + IT_0020 + IT_0127 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0341 = IT_0339*IT_0340;
    const ccomplex_t IT_0342 = (0 + _Complex_I*1)*IT_0341;
    const ccomplex_t IT_0343 = -IT_0342;
    const ccomplex_t IT_0344 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0345 = IT_0110*IT_0344;
    const ccomplex_t IT_0346 = IT_0113*IT_0344;
    const ccomplex_t IT_0347 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0348 = IT_0113*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*(IT_0345 + -IT_0346 + (-2)
      *IT_0348);
    const ccomplex_t IT_0350 = 0.5*IT_0349;
    const ccomplex_t IT_0351 = IT_0057*IT_0350;
    const ccomplex_t IT_0352 = cpow((-2)*s_13 + IT_0127 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0353 = IT_0351*IT_0352;
    const ccomplex_t IT_0354 = (0 + _Complex_I*1)*IT_0353;
    const ccomplex_t IT_0355 = -IT_0354;
    const ccomplex_t IT_0356 = conj(IT_0189) + conj(IT_0190) + conj(IT_0191) +
       conj(IT_0343) + conj(IT_0355);
    const ccomplex_t IT_0357 = IT_0189 + IT_0190 + IT_0191 + IT_0343 + IT_0355;
    const ccomplex_t IT_0358 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0359 = IT_0010*IT_0358;
    const ccomplex_t IT_0360 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0361 = IT_0010*IT_0360;
    const ccomplex_t IT_0362 = 1.4142135623731*IT_0361;
    const ccomplex_t IT_0363 = (0 + _Complex_I*1)*(IT_0359 + (-0.5)*IT_0362);
    const ccomplex_t IT_0364 = -IT_0363;
    const ccomplex_t IT_0365 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0366 = IT_0110*IT_0365;
    const ccomplex_t IT_0367 = IT_0113*IT_0365;
    const ccomplex_t IT_0368 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0369 = IT_0113*IT_0368;
    const ccomplex_t IT_0370 = (0 + _Complex_I*1)*(IT_0366 + -IT_0367 + (-2)
      *IT_0369);
    const ccomplex_t IT_0371 = (-0.5)*IT_0370;
    const ccomplex_t IT_0372 = IT_0364*IT_0371;
    const ccomplex_t IT_0373 = IT_0340*IT_0372;
    const ccomplex_t IT_0374 = (0 + _Complex_I*1)*IT_0373;
    const ccomplex_t IT_0375 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0376 = IT_0110*IT_0375;
    const ccomplex_t IT_0377 = IT_0113*IT_0375;
    const ccomplex_t IT_0378 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0379 = IT_0113*IT_0378;
    const ccomplex_t IT_0380 = (0 + _Complex_I*1)*(IT_0376 + -IT_0377 + (-2)
      *IT_0379);
    const ccomplex_t IT_0381 = (-0.5)*IT_0380;
    const ccomplex_t IT_0382 = IT_0019*IT_0381;
    const ccomplex_t IT_0383 = IT_0352*IT_0382;
    const ccomplex_t IT_0384 = (0 + _Complex_I*1)*IT_0383;
    const ccomplex_t IT_0385 = IT_0262 + conj(IT_0374) + conj(IT_0384);
    const ccomplex_t IT_0386 = IT_0260 + IT_0374 + IT_0384;
    const ccomplex_t IT_0387 = pow(s_13, 2);
    const ccomplex_t IT_0388 = IT_0002*IT_0387;
    const ccomplex_t IT_0389 = pow(s_14, 2);
    const ccomplex_t IT_0390 = IT_0001*IT_0389;
    const ccomplex_t IT_0391 = IT_0388 + IT_0390;
    const ccomplex_t IT_0392 = s_13*s_14*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0393 = -IT_0020 + -IT_0392;
    const ccomplex_t IT_0394 = IT_0391 + IT_0393;
    const ccomplex_t IT_0395 = m_N_2*IT_0394;
    const ccomplex_t IT_0396 = -IT_0133 + 2*IT_0342 + 2*IT_0354;
    const ccomplex_t IT_0397 = 2*IT_0050;
    const ccomplex_t IT_0398 = (-2)*IT_0090;
    const ccomplex_t IT_0399 = IT_0397 + IT_0398;
    const ccomplex_t IT_0400 = 2*conj(IT_0050);
    const ccomplex_t IT_0401 = (-2)*conj(IT_0090);
    const ccomplex_t IT_0402 = IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = s_13*IT_0127;
    const ccomplex_t IT_0404 = IT_0002*IT_0403;
    const ccomplex_t IT_0405 = s_14*IT_0127;
    const ccomplex_t IT_0406 = IT_0313*IT_0405;
    const ccomplex_t IT_0407 = s_13*IT_0046;
    const ccomplex_t IT_0408 = s_14*s_34;
    const ccomplex_t IT_0409 = (-2)*IT_0408;
    const ccomplex_t IT_0410 = IT_0407 + IT_0409;
    const ccomplex_t IT_0411 = IT_0001*IT_0410;
    const ccomplex_t IT_0412 = 0.5*IT_0404 + (-0.5)*IT_0406 + (-0.5)*IT_0411;
    const ccomplex_t IT_0413 = s_13 + IT_0412;
    const ccomplex_t IT_0414 = m_N_2*IT_0413;
    const ccomplex_t IT_0415 = 2*IT_0414;
    const ccomplex_t IT_0416 = 2*IT_0082;
    const ccomplex_t IT_0417 = (-2)*IT_0060;
    const ccomplex_t IT_0418 = IT_0416 + IT_0417;
    const ccomplex_t IT_0419 = s_13*IT_0046*IT_0127;
    const ccomplex_t IT_0420 = s_14*s_34*IT_0127;
    const ccomplex_t IT_0421 = (-2)*IT_0420;
    const ccomplex_t IT_0422 = IT_0419 + IT_0421;
    const ccomplex_t IT_0423 = IT_0263*IT_0422;
    const ccomplex_t IT_0424 = (-0.666666666666667)*IT_0411 +
       0.166666666666667*IT_0423;
    const ccomplex_t IT_0425 = s_13 + IT_0424;
    const ccomplex_t IT_0426 = m_N_2*IT_0425;
    const ccomplex_t IT_0427 = 6*IT_0426;
    const ccomplex_t IT_0428 = m_N_2*IT_0141;
    const ccomplex_t IT_0429 = m_N_2*IT_0185;
    const ccomplex_t IT_0430 = IT_0254*(IT_0428 + -IT_0429);
    const ccomplex_t IT_0431 = 2*IT_0430;
    const ccomplex_t IT_0432 = m_N_1*IT_0173;
    const ccomplex_t IT_0433 = m_N_2*IT_0241;
    const ccomplex_t IT_0434 = IT_0250*(IT_0432 + -IT_0433);
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = 2*IT_0435;
    const ccomplex_t IT_0437 = IT_0431 + IT_0436;
    const ccomplex_t IT_0438 = m_N_3*IT_0119;
    const ccomplex_t IT_0439 = m_N_2*IT_0204;
    const ccomplex_t IT_0440 = -IT_0439;
    const ccomplex_t IT_0441 = IT_0438 + IT_0440;
    const ccomplex_t IT_0442 = IT_0213*IT_0441;
    const ccomplex_t IT_0443 = m_N_4*IT_0153;
    const ccomplex_t IT_0444 = m_N_2*IT_0222;
    const ccomplex_t IT_0445 = -IT_0444;
    const ccomplex_t IT_0446 = IT_0443 + IT_0445;
    const ccomplex_t IT_0447 = IT_0231*IT_0446;
    const ccomplex_t IT_0448 = -IT_0442 + -IT_0447;
    const ccomplex_t IT_0449 = (-2)*IT_0448;
    const ccomplex_t IT_0450 = IT_0437 + IT_0449;
    const ccomplex_t IT_0451 = 2*conj(IT_0430);
    const ccomplex_t IT_0452 = 2*conj(IT_0435);
    const ccomplex_t IT_0453 = IT_0451 + IT_0452;
    const ccomplex_t IT_0454 = (-2)*conj(IT_0448);
    const ccomplex_t IT_0455 = IT_0453 + IT_0454;
    const ccomplex_t IT_0456 = IT_0144*(IT_0428 + -IT_0429);
    const ccomplex_t IT_0457 = (-0.5)*IT_0456;
    const ccomplex_t IT_0458 = 2*IT_0457;
    const ccomplex_t IT_0459 = m_N_1*IT_0241;
    const ccomplex_t IT_0460 = m_N_2*IT_0173;
    const ccomplex_t IT_0461 = IT_0183*(IT_0459 + -IT_0460);
    const ccomplex_t IT_0462 = -IT_0461;
    const ccomplex_t IT_0463 = 2*IT_0462;
    const ccomplex_t IT_0464 = IT_0458 + IT_0463;
    const ccomplex_t IT_0465 = -IT_0429;
    const ccomplex_t IT_0466 = IT_0428 + IT_0465;
    const ccomplex_t IT_0467 = IT_0144*IT_0466;
    const ccomplex_t IT_0468 = m_N_3*IT_0204;
    const ccomplex_t IT_0469 = m_N_2*IT_0119;
    const ccomplex_t IT_0470 = -IT_0469;
    const ccomplex_t IT_0471 = IT_0468 + IT_0470;
    const ccomplex_t IT_0472 = IT_0130*IT_0471;
    const ccomplex_t IT_0473 = m_N_4*IT_0222;
    const ccomplex_t IT_0474 = m_N_2*IT_0153;
    const ccomplex_t IT_0475 = -IT_0474;
    const ccomplex_t IT_0476 = IT_0473 + IT_0475;
    const ccomplex_t IT_0477 = IT_0163*IT_0476;
    const ccomplex_t IT_0478 = 0.5*IT_0467 + -IT_0472 + -IT_0477;
    const ccomplex_t IT_0479 = (-2)*IT_0478;
    const ccomplex_t IT_0480 = IT_0464 + IT_0479;
    const ccomplex_t IT_0481 = 2*conj(IT_0457);
    const ccomplex_t IT_0482 = 2*conj(IT_0462);
    const ccomplex_t IT_0483 = IT_0481 + IT_0482;
    const ccomplex_t IT_0484 = (-2)*conj(IT_0478);
    const ccomplex_t IT_0485 = IT_0483 + IT_0484;
    const ccomplex_t IT_0486 = m_C_2*m_N_2*IT_0001;
    const ccomplex_t IT_0487 = IT_0002*IT_0265;
    const ccomplex_t IT_0488 = (-0.25)*IT_0487;
    const ccomplex_t IT_0489 = IT_0003 + IT_0488;
    const ccomplex_t IT_0490 = IT_0486*IT_0489;
    const ccomplex_t IT_0491 = (-4)*IT_0490;
    const ccomplex_t IT_0492 = (-2)*IT_0343;
    const ccomplex_t IT_0493 = (-2)*IT_0355;
    const ccomplex_t IT_0494 = IT_0492 + IT_0493;
    const ccomplex_t IT_0495 = (-2)*conj(IT_0343);
    const ccomplex_t IT_0496 = (-2)*conj(IT_0355);
    const ccomplex_t IT_0497 = IT_0495 + IT_0496;
    const ccomplex_t IT_0498 = m_C_2*IT_0309;
    const ccomplex_t IT_0499 = s_24*s_34*IT_0001*IT_0002*IT_0387;
    const ccomplex_t IT_0500 = s_12*s_13;
    const ccomplex_t IT_0501 = s_13*s_14*s_23*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0502 = s_23*IT_0020;
    const ccomplex_t IT_0503 = IT_0500 + IT_0501 + IT_0502;
    const ccomplex_t IT_0504 = s_23*IT_0002*IT_0387;
    const ccomplex_t IT_0505 = (-2)*IT_0504;
    const ccomplex_t IT_0506 = s_12*s_13*IT_0001*IT_0002*IT_0003;
    const ccomplex_t IT_0507 = -IT_0506;
    const ccomplex_t IT_0508 = s_24*s_34*IT_0001*IT_0020;
    const ccomplex_t IT_0509 = -IT_0508;
    const ccomplex_t IT_0510 = IT_0499 + IT_0503 + IT_0505 + IT_0507 + IT_0509;
    const ccomplex_t IT_0511 = -IT_0011;
    const ccomplex_t IT_0512 = IT_0059*IT_0511;
    const ccomplex_t IT_0513 = 2*conj(IT_0512);
    const ccomplex_t IT_0514 = s_34*m_C_2*IT_0001;
    const ccomplex_t IT_0515 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0516 = -IT_0515;
    const ccomplex_t IT_0517 = s_24 + IT_0516;
    const ccomplex_t IT_0518 = IT_0514*IT_0517;
    const ccomplex_t IT_0519 = s_34*m_N_2*IT_0001;
    const ccomplex_t IT_0520 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0521 = -IT_0520;
    const ccomplex_t IT_0522 = s_14 + IT_0521;
    const ccomplex_t IT_0523 = IT_0519*IT_0522;
    const ccomplex_t IT_0524 = IT_0001*IT_0003;
    const ccomplex_t IT_0525 = s_12*IT_0127;
    const ccomplex_t IT_0526 = s_13*s_23;
    const ccomplex_t IT_0527 = (-2)*IT_0526;
    const ccomplex_t IT_0528 = IT_0525 + IT_0527;
    const ccomplex_t IT_0529 = IT_0002*IT_0528;
    const ccomplex_t IT_0530 = (IT_0127 + -IT_0524)*(s_12 + (-0.5)*IT_0529);
    const ccomplex_t IT_0531 = (-2)*IT_0530;
    const ccomplex_t IT_0532 = 2*IT_0512;
    const ccomplex_t IT_0533 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0534 = (s_12 + (-0.5)*IT_0529)*(s_13 + -IT_0533);
    const ccomplex_t IT_0535 = (-2)*IT_0534;
    const ccomplex_t IT_0536 = s_24*s_34*IT_0001;
    const ccomplex_t IT_0537 = (s_12 + (-0.5)*IT_0529)*(s_23 + -IT_0536);
    const ccomplex_t IT_0538 = (-2)*IT_0537;
    const ccomplex_t IT_0539 = IT_0535 + IT_0538;
    const ccomplex_t IT_0540 = 2*IT_0284;
    const ccomplex_t IT_0541 = 2*conj(IT_0284);
    const ccomplex_t IT_0542 = 2*IT_0501;
    const ccomplex_t IT_0543 = s_13*s_14*s_24*IT_0001;
    const ccomplex_t IT_0544 = -IT_0543;
    const ccomplex_t IT_0545 = s_23*IT_0001*IT_0389;
    const ccomplex_t IT_0546 = -IT_0545;
    const ccomplex_t IT_0547 = IT_0500 + IT_0502 + IT_0505 + IT_0542 + IT_0544
       + IT_0546;
    const ccomplex_t IT_0548 = s_13*s_23*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0549 = 2*IT_0548;
    const ccomplex_t IT_0550 = s_13*IT_0001*IT_0276;
    const ccomplex_t IT_0551 = -IT_0550;
    const ccomplex_t IT_0552 = s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0553 = -IT_0552;
    const ccomplex_t IT_0554 = s_13*IT_0002*IT_0274;
    const ccomplex_t IT_0555 = (-2)*IT_0554;
    const ccomplex_t IT_0556 = IT_0097 + IT_0098 + IT_0549 + IT_0551 + IT_0553
       + IT_0555;
    const ccomplex_t IT_0557 = IT_0547 + IT_0556;
    const ccomplex_t IT_0558 = 2*IT_0011;
    const ccomplex_t IT_0559 = IT_0059*IT_0558;
    const ccomplex_t IT_0560 = IT_0540*conj(IT_0559);
    const ccomplex_t IT_0561 = IT_0541*IT_0559;
    const ccomplex_t IT_0562 = IT_0001*IT_0002*IT_0274*IT_0389;
    const ccomplex_t IT_0563 = s_13*s_14*s_23*s_24*IT_0001*IT_0002;
    const ccomplex_t IT_0564 = IT_0020*IT_0021;
    const ccomplex_t IT_0565 = IT_0001*IT_0021*IT_0389;
    const ccomplex_t IT_0566 = -IT_0565;
    const ccomplex_t IT_0567 = IT_0002*IT_0020*IT_0274;
    const ccomplex_t IT_0568 = -IT_0567;
    const ccomplex_t IT_0569 = s_12*s_14*s_23*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0570 = -IT_0569;
    const ccomplex_t IT_0571 = IT_0562 + IT_0563 + IT_0564 + IT_0566 + IT_0568
       + IT_0570;
    const ccomplex_t IT_0572 = s_14*s_24*IT_0001*IT_0002*IT_0274;
    const ccomplex_t IT_0573 = s_13*s_23*IT_0001*IT_0002*IT_0276;
    const ccomplex_t IT_0574 = s_12*IT_0021;
    const ccomplex_t IT_0575 = s_14*s_24*IT_0001*IT_0021;
    const ccomplex_t IT_0576 = -IT_0575;
    const ccomplex_t IT_0577 = s_13*s_23*IT_0002*IT_0021;
    const ccomplex_t IT_0578 = -IT_0577;
    const ccomplex_t IT_0579 = s_12*s_23*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0580 = -IT_0579;
    const ccomplex_t IT_0581 = IT_0572 + IT_0573 + IT_0574 + IT_0576 + IT_0578
       + IT_0580;
    const ccomplex_t IT_0582 = 2*IT_0233;
    const ccomplex_t IT_0583 = conj(IT_0295)*IT_0582;
    const ccomplex_t IT_0584 = 2*conj(IT_0233);
    const ccomplex_t IT_0585 = IT_0295*IT_0584;
    const ccomplex_t IT_0586 = conj(IT_0165)*IT_0540;
    const ccomplex_t IT_0587 = conj(IT_0187)*IT_0540;
    const ccomplex_t IT_0588 = IT_0165*IT_0541;
    const ccomplex_t IT_0589 = IT_0187*IT_0541;
    const ccomplex_t IT_0590 = 2*IT_0295;
    const ccomplex_t IT_0591 = conj(IT_0256)*IT_0590;
    const ccomplex_t IT_0592 = 2*conj(IT_0295);
    const ccomplex_t IT_0593 = IT_0256*IT_0592;
    const ccomplex_t IT_0594 = IT_0024*IT_0511;
    const ccomplex_t IT_0595 = IT_0590*conj(IT_0594);
    const ccomplex_t IT_0596 = IT_0592*IT_0594;
    const ccomplex_t IT_0597 = s_13*s_23*IT_0001*IT_0002*IT_0389;
    const ccomplex_t IT_0598 = s_14*s_24*IT_0001*IT_0002*IT_0387;
    const ccomplex_t IT_0599 = s_12*IT_0020;
    const ccomplex_t IT_0600 = s_14*s_24*IT_0001*IT_0020;
    const ccomplex_t IT_0601 = -IT_0600;
    const ccomplex_t IT_0602 = s_13*s_23*IT_0002*IT_0020;
    const ccomplex_t IT_0603 = -IT_0602;
    const ccomplex_t IT_0604 = s_12*s_13*s_14*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0605 = -IT_0604;
    const ccomplex_t IT_0606 = IT_0597 + IT_0598 + IT_0599 + IT_0601 + IT_0603
       + IT_0605;
    const ccomplex_t IT_0607 = IT_0001*IT_0002*IT_0276*IT_0387;
    const ccomplex_t IT_0608 = IT_0002*IT_0021*IT_0387;
    const ccomplex_t IT_0609 = -IT_0608;
    const ccomplex_t IT_0610 = IT_0001*IT_0020*IT_0276;
    const ccomplex_t IT_0611 = -IT_0610;
    const ccomplex_t IT_0612 = s_12*s_13*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0613 = -IT_0612;
    const ccomplex_t IT_0614 = IT_0563 + IT_0564 + IT_0607 + IT_0609 + IT_0611
       + IT_0613;
    const ccomplex_t IT_0615 = conj(IT_0396)*IT_0540;
    const ccomplex_t IT_0616 = IT_0396*IT_0541;
    const ccomplex_t IT_0617 = -IT_0252 + (-2)*IT_0374 + (-2)*IT_0384;
    const ccomplex_t IT_0618 = IT_0590*conj(IT_0617);
    const ccomplex_t IT_0619 = IT_0592*IT_0617;
    const ccomplex_t IT_0620 = s_14*s_34*IT_0001*IT_0127;
    const ccomplex_t IT_0621 = -IT_0620;
    const ccomplex_t IT_0622 = IT_0403 + IT_0621;
    const ccomplex_t IT_0623 = IT_0002*IT_0622;
    const ccomplex_t IT_0624 = -IT_0533 + (-0.25)*IT_0623;
    const ccomplex_t IT_0625 = s_13 + IT_0624;
    const ccomplex_t IT_0626 = IT_0000*IT_0625;
    const ccomplex_t IT_0627 = (-4)*IT_0626;
    const ccomplex_t IT_0628 = s_24*s_34*IT_0001*IT_0127;
    const ccomplex_t IT_0629 = -IT_0628;
    const ccomplex_t IT_0630 = IT_0311 + IT_0629;
    const ccomplex_t IT_0631 = IT_0002*IT_0630;
    const ccomplex_t IT_0632 = -IT_0536 + (-0.25)*IT_0631;
    const ccomplex_t IT_0633 = s_23 + IT_0632;
    const ccomplex_t IT_0634 = IT_0000*IT_0633;
    const ccomplex_t IT_0635 = (-4)*IT_0634;
    const ccomplex_t IT_0636 = conj(IT_0258)*IT_0540;
    const ccomplex_t IT_0637 = conj(IT_0255)*IT_0540;
    const ccomplex_t IT_0638 = conj(IT_0259)*IT_0540;
    const ccomplex_t IT_0639 = conj(IT_0374)*IT_0540;
    const ccomplex_t IT_0640 = conj(IT_0384)*IT_0540;
    const ccomplex_t IT_0641 = IT_0258*IT_0541;
    const ccomplex_t IT_0642 = IT_0255*IT_0541;
    const ccomplex_t IT_0643 = IT_0259*IT_0541;
    const ccomplex_t IT_0644 = IT_0374*IT_0541;
    const ccomplex_t IT_0645 = IT_0384*IT_0541;
    const ccomplex_t IT_0646 = conj(IT_0343)*IT_0590;
    const ccomplex_t IT_0647 = conj(IT_0190)*IT_0590;
    const ccomplex_t IT_0648 = conj(IT_0355)*IT_0590;
    const ccomplex_t IT_0649 = conj(IT_0191)*IT_0590;
    const ccomplex_t IT_0650 = conj(IT_0189)*IT_0590;
    const ccomplex_t IT_0651 = IT_0343*IT_0592;
    const ccomplex_t IT_0652 = IT_0190*IT_0592;
    const ccomplex_t IT_0653 = IT_0355*IT_0592;
    const ccomplex_t IT_0654 = IT_0191*IT_0592;
    const ccomplex_t IT_0655 = IT_0189*IT_0592;
    const ccomplex_t IT_0656 = (IT_0627 + IT_0635)*(IT_0636 + IT_0637 +
       IT_0638 + IT_0639 + IT_0640 + IT_0641 + IT_0642 + IT_0643 + IT_0644 +
       IT_0645 + IT_0646 + IT_0647 + IT_0648 + IT_0649 + IT_0650 + IT_0651 +
       IT_0652 + IT_0653 + IT_0654 + IT_0655);
    const ccomplex_t IT_0657 = IT_0001*IT_0003*IT_0127;
    const ccomplex_t IT_0658 = -IT_0657;
    const ccomplex_t IT_0659 = IT_0264 + IT_0658;
    const ccomplex_t IT_0660 = IT_0002*IT_0659;
    const ccomplex_t IT_0661 = -IT_0127 + 0.25*IT_0660;
    const ccomplex_t IT_0662 = IT_0524 + IT_0661;
    const ccomplex_t IT_0663 = IT_0000*IT_0662;
    const ccomplex_t IT_0664 = 4*IT_0663;
    const ccomplex_t IT_0665 = 2*IT_0594;
    const ccomplex_t IT_0666 = 2*conj(IT_0594);
    const ccomplex_t IT_0667 = m_N_2*IT_0303;
    const ccomplex_t IT_0668 = s_14*s_34*IT_0001*IT_0002*IT_0274;
    const ccomplex_t IT_0669 = IT_0097 + IT_0098 + IT_0548;
    const ccomplex_t IT_0670 = s_12*s_23*IT_0001*IT_0002*IT_0003;
    const ccomplex_t IT_0671 = -IT_0670;
    const ccomplex_t IT_0672 = s_14*s_34*IT_0001*IT_0021;
    const ccomplex_t IT_0673 = -IT_0672;
    const ccomplex_t IT_0674 = IT_0555 + IT_0668 + IT_0669 + IT_0671 + IT_0673;
    const ccomplex_t IT_0675 = s_12*IT_0046;
    const ccomplex_t IT_0676 = s_14*s_24;
    const ccomplex_t IT_0677 = (-2)*IT_0676;
    const ccomplex_t IT_0678 = IT_0675 + IT_0677;
    const ccomplex_t IT_0679 = IT_0001*IT_0678;
    const ccomplex_t IT_0680 = (s_12 + -IT_0300)*(s_12 + (-0.5)*IT_0679);
    const ccomplex_t IT_0681 = (-2)*IT_0680;
    const ccomplex_t IT_0682 = 2*conj(IT_0396);
    const ccomplex_t IT_0683 = s_34*IT_0127;
    const ccomplex_t IT_0684 = IT_0313*IT_0683;
    const ccomplex_t IT_0685 = (-2)*IT_0003;
    const ccomplex_t IT_0686 = IT_0267 + IT_0685;
    const ccomplex_t IT_0687 = IT_0001*IT_0686;
    const ccomplex_t IT_0688 = (-0.333333333333333)*IT_0684 + (
      -0.333333333333333)*IT_0687;
    const ccomplex_t IT_0689 = IT_0127 + IT_0688;
    const ccomplex_t IT_0690 = IT_0000*IT_0689;
    const ccomplex_t IT_0691 = 3*IT_0690;
    const ccomplex_t IT_0692 = IT_0024*IT_0558;
    const ccomplex_t IT_0693 = 2*IT_0692;
    const ccomplex_t IT_0694 = 2*conj(IT_0692);
    const ccomplex_t IT_0695 = 2*IT_0559;
    const ccomplex_t IT_0696 = 2*conj(IT_0559);
    const ccomplex_t IT_0697 = (s_12 + -IT_0301)*(s_12 + (-0.5)*IT_0529);
    const ccomplex_t IT_0698 = (-2)*IT_0697;
    const ccomplex_t IT_0699 = 4*IT_0284;
    const ccomplex_t IT_0700 = s_12*s_14*s_34*IT_0001;
    const ccomplex_t IT_0701 = -IT_0502;
    const ccomplex_t IT_0702 = s_23*IT_0001*IT_0127*IT_0389;
    const ccomplex_t IT_0703 = s_13*s_14*s_24*IT_0001*IT_0127;
    const ccomplex_t IT_0704 = IT_0702 + IT_0703;
    const ccomplex_t IT_0705 = s_12*s_14*s_34*IT_0001*IT_0127;
    const ccomplex_t IT_0706 = s_23*IT_0020*IT_0127;
    const ccomplex_t IT_0707 = -IT_0705 + -IT_0706;
    const ccomplex_t IT_0708 = IT_0704 + IT_0707;
    const ccomplex_t IT_0709 = IT_0002*IT_0708;
    const ccomplex_t IT_0710 = 0.5*IT_0709;
    const ccomplex_t IT_0711 = IT_0544 + IT_0545 + IT_0700 + IT_0701 + IT_0710;
    const ccomplex_t IT_0712 = 4*conj(IT_0284);
    const ccomplex_t IT_0713 = 4*IT_0295;
    const ccomplex_t IT_0714 = 4*conj(IT_0295);
    const ccomplex_t IT_0715 = -IT_0020;
    const ccomplex_t IT_0716 = IT_0390 + IT_0715;
    const ccomplex_t IT_0717 = IT_0000*IT_0716;
    const ccomplex_t IT_0718 = -IT_0021;
    const ccomplex_t IT_0719 = IT_0277 + IT_0718;
    const ccomplex_t IT_0720 = IT_0000*IT_0719;
    const ccomplex_t IT_0721 = conj(IT_0284)*IT_0295;
    const ccomplex_t IT_0722 = 6*IT_0721;
    const ccomplex_t IT_0723 = s_12*s_24*s_34*IT_0001;
    const ccomplex_t IT_0724 = (-0.5)*IT_0723;
    const ccomplex_t IT_0725 = s_14*s_23*s_24*IT_0001*IT_0127;
    const ccomplex_t IT_0726 = s_13*IT_0001*IT_0127*IT_0276;
    const ccomplex_t IT_0727 = IT_0725 + IT_0726;
    const ccomplex_t IT_0728 = s_12*s_24*s_34*IT_0001*IT_0127;
    const ccomplex_t IT_0729 = s_13*IT_0021*IT_0127;
    const ccomplex_t IT_0730 = -IT_0728 + -IT_0729;
    const ccomplex_t IT_0731 = IT_0727 + IT_0730;
    const ccomplex_t IT_0732 = IT_0002*IT_0731;
    const ccomplex_t IT_0733 = (-0.25)*IT_0732;
    const ccomplex_t IT_0734 = (-0.5)*IT_0098;
    const ccomplex_t IT_0735 = 0.5*IT_0550;
    const ccomplex_t IT_0736 = (-0.5)*IT_0552;
    const ccomplex_t IT_0737 = IT_0097 + IT_0724 + IT_0733 + IT_0734 + IT_0735
       + IT_0736;
    const ccomplex_t IT_0738 = (-8)*IT_0284;
    const ccomplex_t IT_0739 = (-8)*conj(IT_0284);
    const ccomplex_t IT_0740 = (-8)*IT_0295;
    const ccomplex_t IT_0741 = (-8)*conj(IT_0295);
    const ccomplex_t IT_0742 = s_12*IT_0046*IT_0264;
    const ccomplex_t IT_0743 = s_14*s_24*IT_0264;
    const ccomplex_t IT_0744 = (-2)*IT_0743;
    const ccomplex_t IT_0745 = IT_0742 + IT_0744;
    const ccomplex_t IT_0746 = IT_0263*IT_0745;
    const ccomplex_t IT_0747 = 0.25*IT_0746;
    const ccomplex_t IT_0748 = IT_0525 + IT_0747;
    const ccomplex_t IT_0749 = conj(IT_0189) + conj(IT_0191);
    const ccomplex_t IT_0750 = (-8)*IT_0355;
    const ccomplex_t IT_0751 = IT_0189 + IT_0191;
    const ccomplex_t IT_0752 = (-8)*conj(IT_0355);
    const ccomplex_t IT_0753 = -IT_0533;
    const ccomplex_t IT_0754 = s_13 + IT_0753;
    const ccomplex_t IT_0755 = IT_0000*IT_0754;
    const ccomplex_t IT_0756 = -IT_0536;
    const ccomplex_t IT_0757 = s_23 + IT_0756;
    const ccomplex_t IT_0758 = IT_0000*IT_0757;
    const ccomplex_t IT_0759 = IT_0755 + IT_0758;
    const ccomplex_t IT_0760 = (-6)*conj(IT_0512);
    const ccomplex_t IT_0761 = IT_0295*IT_0760;
    const ccomplex_t IT_0762 = IT_0284*conj(IT_0594);
    const ccomplex_t IT_0763 = (-6)*IT_0762;
    const ccomplex_t IT_0764 = conj(IT_0284)*IT_0594;
    const ccomplex_t IT_0765 = (-6)*IT_0764;
    const ccomplex_t IT_0766 = s_12*IT_0006;
    const ccomplex_t IT_0767 = (-4)*conj(IT_0090);
    const ccomplex_t IT_0768 = 0.25*IT_0767;
    const ccomplex_t IT_0769 = conj(IT_0050) + IT_0768;
    const ccomplex_t IT_0770 = (-4)*IT_0090;
    const ccomplex_t IT_0771 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0772 = s_12*IT_0001*IT_0003;
    const ccomplex_t IT_0773 = -IT_0772;
    const ccomplex_t IT_0774 = s_12*IT_0264;
    const ccomplex_t IT_0775 = s_14*s_23*s_34*IT_0001*IT_0127;
    const ccomplex_t IT_0776 = s_13*s_24*s_34*IT_0001*IT_0127;
    const ccomplex_t IT_0777 = IT_0774 + IT_0775 + IT_0776;
    const ccomplex_t IT_0778 = s_13*s_23*IT_0127;
    const ccomplex_t IT_0779 = s_12*IT_0001*IT_0003*IT_0127;
    const ccomplex_t IT_0780 = (-2)*IT_0778 + -IT_0779;
    const ccomplex_t IT_0781 = IT_0777 + IT_0780;
    const ccomplex_t IT_0782 = IT_0002*IT_0781;
    const ccomplex_t IT_0783 = (-0.5)*IT_0782;
    const ccomplex_t IT_0784 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0785 = -IT_0784;
    const ccomplex_t IT_0786 = IT_0525 + IT_0771 + IT_0773 + IT_0783 + IT_0785;
    const ccomplex_t IT_0787 = (-4)*IT_0512;
    const ccomplex_t IT_0788 = (-4)*conj(IT_0512);
    const ccomplex_t IT_0789 = (-4)*IT_0594;
    const ccomplex_t IT_0790 = (-4)*conj(IT_0594);
    const ccomplex_t IT_0791 = IT_0500 + IT_0701;
    const ccomplex_t IT_0792 = s_12*s_13*IT_0046;
    const ccomplex_t IT_0793 = s_23*IT_0389;
    const ccomplex_t IT_0794 = s_12*s_14*s_34;
    const ccomplex_t IT_0795 = s_23*IT_0020*IT_0046;
    const ccomplex_t IT_0796 = 2*IT_0793 + (-2)*IT_0794 + -IT_0795;
    const ccomplex_t IT_0797 = IT_0792 + IT_0796;
    const ccomplex_t IT_0798 = IT_0001*IT_0797;
    const ccomplex_t IT_0799 = (-0.5)*IT_0798;
    const ccomplex_t IT_0800 = IT_0791 + IT_0799;
    const ccomplex_t IT_0801 = (-4)*IT_0396;
    const ccomplex_t IT_0802 = (-4)*conj(IT_0396);
    const ccomplex_t IT_0803 = (-4)*IT_0617;
    const ccomplex_t IT_0804 = (-4)*conj(IT_0617);
    const ccomplex_t IT_0805 = 0.5*IT_0746;
    const ccomplex_t IT_0806 = (-4)*IT_0526;
    const ccomplex_t IT_0807 = IT_0525 + IT_0805 + IT_0806;
    const ccomplex_t IT_0808 = s_12*IT_0046*IT_0127;
    const ccomplex_t IT_0809 = s_13*s_23*IT_0046;
    const ccomplex_t IT_0810 = s_14*s_23*s_34;
    const ccomplex_t IT_0811 = s_14*s_24*IT_0127;
    const ccomplex_t IT_0812 = (-2)*IT_0809 + 4*IT_0810 + (-2)*IT_0811;
    const ccomplex_t IT_0813 = IT_0808 + IT_0812;
    const ccomplex_t IT_0814 = IT_0001*IT_0813;
    const ccomplex_t IT_0815 = -IT_0814;
    const ccomplex_t IT_0816 = IT_0807 + IT_0815;
    const ccomplex_t IT_0817 = (-4)*IT_0190;
    const ccomplex_t IT_0818 = (-4)*conj(IT_0190);
    const ccomplex_t IT_0819 = (-4)*IT_0191;
    const ccomplex_t IT_0820 = (-4)*IT_0258;
    const ccomplex_t IT_0821 = (-4)*conj(IT_0258);
    const ccomplex_t IT_0822 = (-4)*IT_0255;
    const ccomplex_t IT_0823 = conj(IT_0190)*IT_0817 + IT_0751*IT_0818 + (-4)
      *IT_0749*(IT_0189 + (-0.25)*IT_0817 + (-0.25)*IT_0819) + conj(IT_0258)
      *IT_0820 + IT_0272*IT_0821 + (-4)*IT_0273*(IT_0259 + (-0.25)*IT_0820 + (
      -0.25)*IT_0822);
    const ccomplex_t IT_0824 = s_34*m_C_2*m_N_2*IT_0001;
    const ccomplex_t IT_0825 = IT_0517*IT_0824;
    const ccomplex_t IT_0826 = IT_0532 + IT_0695;
    const ccomplex_t IT_0827 = s_13*s_23*IT_0002*IT_0046;
    const ccomplex_t IT_0828 = -IT_0827;
    const ccomplex_t IT_0829 = IT_0675 + IT_0828;
    const ccomplex_t IT_0830 = IT_0001*IT_0829;
    const ccomplex_t IT_0831 = -IT_0300 + (-0.25)*IT_0830;
    const ccomplex_t IT_0832 = s_12 + IT_0831;
    const ccomplex_t IT_0833 = m_N_2*IT_0832;
    const ccomplex_t IT_0834 = (-4)*IT_0833;
    const ccomplex_t IT_0835 = m_C_2*IT_0354;
    const ccomplex_t IT_0836 = (-2)*IT_0835;
    const ccomplex_t IT_0837 = m_C_2*IT_0342;
    const ccomplex_t IT_0838 = (-2)*IT_0837;
    const ccomplex_t IT_0839 = IT_0836 + IT_0838;
    const ccomplex_t IT_0840 = IT_0000*IT_0281;
    const ccomplex_t IT_0841 = s_14*s_23*s_34*IT_0002;
    const ccomplex_t IT_0842 = (-2)*IT_0676 + -IT_0827 + 2*IT_0841;
    const ccomplex_t IT_0843 = IT_0675 + IT_0842;
    const ccomplex_t IT_0844 = IT_0001*IT_0843;
    const ccomplex_t IT_0845 = -IT_0300 + (-0.5)*IT_0844;
    const ccomplex_t IT_0846 = s_12 + IT_0845;
    const ccomplex_t IT_0847 = m_N_2*IT_0846;
    const ccomplex_t IT_0848 = (-2)*IT_0847;
    const ccomplex_t IT_0849 = 2*IT_0448;
    const ccomplex_t IT_0850 = IT_0275 + IT_0718;
    const ccomplex_t IT_0851 = IT_0000*IT_0850;
    const ccomplex_t IT_0852 = 6*IT_0165;
    const ccomplex_t IT_0853 = -IT_0300;
    const ccomplex_t IT_0854 = s_12 + IT_0853;
    const ccomplex_t IT_0855 = IT_0000*IT_0854;
    const ccomplex_t IT_0856 = (-6)*IT_0396;
    const ccomplex_t IT_0857 = -(IT_0021 + -IT_0275)*(s_12 + (-0.5)*IT_0679);
    const ccomplex_t IT_0858 = (-2)*IT_0857;
    const ccomplex_t IT_0859 = (-2)*IT_0233;
    const ccomplex_t IT_0860 = IT_0000*IT_0303;
    const ccomplex_t IT_0861 = (-2)*IT_0284;
    const ccomplex_t IT_0862 = -(IT_0020 + -IT_0390)*(s_12 + (-0.5)*IT_0529);
    const ccomplex_t IT_0863 = (-2)*IT_0862;
    const ccomplex_t IT_0864 = -(IT_0021 + -IT_0277)*(s_12 + (-0.5)*IT_0529);
    const ccomplex_t IT_0865 = (-2)*IT_0864;
    const ccomplex_t IT_0866 = IT_0863 + IT_0865;
    const ccomplex_t IT_0867 = conj(IT_0284)*IT_0861;
    const ccomplex_t IT_0868 = 2*IT_0396;
    const ccomplex_t IT_0869 = (-2)*IT_0165;
    const ccomplex_t IT_0870 = IT_0002*IT_0525;
    const ccomplex_t IT_0871 = (-0.25)*IT_0870;
    const ccomplex_t IT_0872 = IT_0001*IT_0675;
    const ccomplex_t IT_0873 = (-0.25)*IT_0872;
    const ccomplex_t IT_0874 = s_12 + IT_0871 + IT_0873;
    const ccomplex_t IT_0875 = s_14*s_23;
    const ccomplex_t IT_0876 = s_12*s_34;
    const ccomplex_t IT_0877 = IT_0875 + IT_0876;
    const ccomplex_t IT_0878 = s_13*s_24;
    const ccomplex_t IT_0879 = -IT_0878;
    const ccomplex_t IT_0880 = IT_0877 + IT_0879;
    const ccomplex_t IT_0881 = IT_0313*IT_0880;
    const ccomplex_t IT_0882 = 0.25*IT_0881;
    const ccomplex_t IT_0883 = IT_0874 + IT_0882;
    const ccomplex_t IT_0884 = 8*conj(IT_0060);
    const ccomplex_t IT_0885 = (-8)*conj(IT_0082);
    const ccomplex_t IT_0886 = IT_0884 + IT_0885;
    const ccomplex_t IT_0887 = IT_0311 + IT_0628;
    const ccomplex_t IT_0888 = s_23*IT_0001*IT_0003;
    const ccomplex_t IT_0889 = (-2)*IT_0888;
    const ccomplex_t IT_0890 = IT_0887 + IT_0889;
    const ccomplex_t IT_0891 = IT_0002*IT_0890;
    const ccomplex_t IT_0892 = -IT_0536 + 0.5*IT_0891;
    const ccomplex_t IT_0893 = s_23 + IT_0892;
    const ccomplex_t IT_0894 = m_C_2*IT_0893;
    const ccomplex_t IT_0895 = 2*IT_0894;
    const ccomplex_t IT_0896 = m_C_2*IT_0832;
    const ccomplex_t IT_0897 = (-4)*IT_0896;
    const ccomplex_t IT_0898 = IT_0001*IT_0316;
    const ccomplex_t IT_0899 = s_23*s_34;
    const ccomplex_t IT_0900 = (-2)*IT_0899;
    const ccomplex_t IT_0901 = IT_0314 + IT_0900;
    const ccomplex_t IT_0902 = IT_0313*IT_0901;
    const ccomplex_t IT_0903 = (-0.25)*IT_0312 + (-0.25)*IT_0898 + (-0.25)
      *IT_0902;
    const ccomplex_t IT_0904 = s_23 + IT_0903;
    const ccomplex_t IT_0905 = m_C_2*IT_0904;
    const ccomplex_t IT_0906 = (-4)*IT_0905;
    const ccomplex_t IT_0907 = IT_0001*IT_0127*IT_0276;
    const ccomplex_t IT_0908 = IT_0021*IT_0127;
    const ccomplex_t IT_0909 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0910 = 0.5*IT_0907 + (-0.5)*IT_0908 + -IT_0909;
    const ccomplex_t IT_0911 = IT_0274 + IT_0910;
    const ccomplex_t IT_0912 = IT_0002*IT_0911;
    const ccomplex_t IT_0913 = -IT_0277;
    const ccomplex_t IT_0914 = IT_0021 + IT_0912 + IT_0913;
    const ccomplex_t IT_0915 = m_C_2*IT_0914;
    const ccomplex_t IT_0916 = s_13*s_24*s_34;
    const ccomplex_t IT_0917 = s_12*IT_0003;
    const ccomplex_t IT_0918 = (-2)*IT_0809 + 2*IT_0810 + (-2)*IT_0811 + 2
      *IT_0916 + (-2)*IT_0917;
    const ccomplex_t IT_0919 = IT_0808 + IT_0918;
    const ccomplex_t IT_0920 = IT_0263*IT_0919;
    const ccomplex_t IT_0921 = (-0.25)*IT_0920;
    const ccomplex_t IT_0922 = -IT_0301;
    const ccomplex_t IT_0923 = s_12 + IT_0853 + IT_0921 + IT_0922;
    const ccomplex_t IT_0924 = 8*conj(IT_0835);
    const ccomplex_t IT_0925 = s_14*s_24*IT_0001*IT_0127;
    const ccomplex_t IT_0926 = (-2)*IT_0526 + 2*IT_0784 + -IT_0925;
    const ccomplex_t IT_0927 = IT_0525 + IT_0926;
    const ccomplex_t IT_0928 = IT_0002*IT_0927;
    const ccomplex_t IT_0929 = -IT_0301 + (-0.5)*IT_0928;
    const ccomplex_t IT_0930 = s_12 + IT_0929;
    const ccomplex_t IT_0931 = m_C_2*IT_0930;
    const ccomplex_t IT_0932 = (-2)*IT_0931;
    const ccomplex_t IT_0933 = (-2)*conj(IT_0284);
    const ccomplex_t IT_0934 = IT_0002*IT_0046*IT_0274;
    const ccomplex_t IT_0935 = IT_0021*IT_0046;
    const ccomplex_t IT_0936 = -IT_0935;
    const ccomplex_t IT_0937 = IT_0934 + IT_0936;
    const ccomplex_t IT_0938 = IT_0001*IT_0937;
    const ccomplex_t IT_0939 = (-4)*IT_0850;
    const ccomplex_t IT_0940 = IT_0938 + IT_0939;
    const ccomplex_t IT_0941 = m_C_2*IT_0940;
    const ccomplex_t IT_0942 = (-2)*conj(IT_0165);
    const ccomplex_t IT_0943 = -(IT_0020 + -IT_0388)*(s_12 + (-0.5)*IT_0679);
    const ccomplex_t IT_0944 = (-2)*IT_0943;
    const ccomplex_t IT_0945 = (-2)*IT_0396;
    const ccomplex_t IT_0946 = m_C_2*IT_0374;
    const ccomplex_t IT_0947 = m_C_2*IT_0384;
    const ccomplex_t IT_0948 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0949 = IT_0338*IT_0364;
    const ccomplex_t IT_0950 = IT_0340*IT_0948*IT_0949;
    const ccomplex_t IT_0951 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0952 = IT_0019*IT_0350;
    const ccomplex_t IT_0953 = IT_0352*IT_0951*IT_0952;
    const ccomplex_t IT_0954 = IT_0946 + IT_0947 + IT_0950 + IT_0953;
    const ccomplex_t IT_0955 = IT_0057*IT_0381;
    const ccomplex_t IT_0956 = IT_0352*IT_0951*IT_0955;
    const ccomplex_t IT_0957 = IT_0331*IT_0371;
    const ccomplex_t IT_0958 = IT_0340*IT_0948*IT_0957;
    const ccomplex_t IT_0959 = -IT_0956 + -IT_0958;
    const ccomplex_t IT_0960 = IT_0263*IT_0267;
    const ccomplex_t IT_0961 = 0.125*IT_0960;
    const ccomplex_t IT_0962 = 1 + IT_0961;
    const ccomplex_t IT_0963 = IT_0000*IT_0962;
    const ccomplex_t IT_0964 = (-8)*IT_0478;
    const ccomplex_t IT_0965 = 0.125*IT_0964;
    const ccomplex_t IT_0966 = 8*IT_0457;
    const ccomplex_t IT_0967 = 0.125*IT_0966;
    const ccomplex_t IT_0968 = IT_0462 + IT_0965 + IT_0967;
    const ccomplex_t IT_0969 = 2*IT_0165;
    const ccomplex_t IT_0970 = 2*IT_0187;
    const ccomplex_t IT_0971 = IT_0969 + IT_0970;
    const ccomplex_t IT_0972 = -IT_0876 + -IT_0878;
    const ccomplex_t IT_0973 = IT_0875 + IT_0972;
    const ccomplex_t IT_0974 = IT_0313*IT_0973;
    const ccomplex_t IT_0975 = (-0.25)*IT_0974;
    const ccomplex_t IT_0976 = IT_0874 + IT_0975;
    const ccomplex_t IT_0977 = 8*IT_0050;
    const ccomplex_t IT_0978 = (-8)*IT_0090;
    const ccomplex_t IT_0979 = IT_0977 + IT_0978;
    const ccomplex_t IT_0980 = s_23*s_24*s_34*IT_0002;
    const ccomplex_t IT_0981 = 2*IT_0276 + -IT_0935 + (-2)*IT_0980;
    const ccomplex_t IT_0982 = IT_0934 + IT_0981;
    const ccomplex_t IT_0983 = IT_0001*IT_0982;
    const ccomplex_t IT_0984 = -IT_0021 + (-0.5)*IT_0983;
    const ccomplex_t IT_0985 = IT_0275 + IT_0984;
    const ccomplex_t IT_0986 = m_C_2*IT_0985;
    const ccomplex_t IT_0987 = (-2)*IT_0986;
    const ccomplex_t IT_0988 = (-2)*IT_0256;
    const ccomplex_t IT_0989 = IT_0859 + IT_0988;
    const ccomplex_t IT_0990 = IT_0403 + IT_0620;
    const ccomplex_t IT_0991 = s_13*IT_0001*IT_0003;
    const ccomplex_t IT_0992 = (-2)*IT_0991;
    const ccomplex_t IT_0993 = IT_0990 + IT_0992;
    const ccomplex_t IT_0994 = IT_0002*IT_0993;
    const ccomplex_t IT_0995 = -IT_0533 + 0.5*IT_0994;
    const ccomplex_t IT_0996 = s_13 + IT_0995;
    const ccomplex_t IT_0997 = m_N_2*IT_0996;
    const ccomplex_t IT_0998 = 2*IT_0997;
    const ccomplex_t IT_0999 = -IT_0908;
    const ccomplex_t IT_1000 = IT_0907 + IT_0999;
    const ccomplex_t IT_1001 = IT_0002*IT_1000;
    const ccomplex_t IT_1002 = (-4)*IT_0719;
    const ccomplex_t IT_1003 = IT_1001 + IT_1002;
    const ccomplex_t IT_1004 = m_C_2*IT_1003;
    const ccomplex_t IT_1005 = s_13*s_24*s_34*IT_0002;
    const ccomplex_t IT_1006 = (-2)*IT_0676 + -IT_0827 + 2*IT_1005;
    const ccomplex_t IT_1007 = IT_0675 + IT_1006;
    const ccomplex_t IT_1008 = IT_0001*IT_1007;
    const ccomplex_t IT_1009 = -IT_0300 + (-0.5)*IT_1008;
    const ccomplex_t IT_1010 = s_12 + IT_1009;
    const ccomplex_t IT_1011 = m_C_2*IT_1010;
    const ccomplex_t IT_1012 = (-2)*IT_1011;
    const ccomplex_t IT_1013 = 2*IT_0617;
    const ccomplex_t IT_1014 = IT_0001*IT_0407;
    const ccomplex_t IT_1015 = s_13*s_34;
    const ccomplex_t IT_1016 = (-2)*IT_1015;
    const ccomplex_t IT_1017 = IT_0405 + IT_1016;
    const ccomplex_t IT_1018 = IT_0313*IT_1017;
    const ccomplex_t IT_1019 = (-0.25)*IT_0404 + (-0.25)*IT_1014 + (-0.25)
      *IT_1018;
    const ccomplex_t IT_1020 = s_13 + IT_1019;
    const ccomplex_t IT_1021 = m_N_2*IT_1020;
    const ccomplex_t IT_1022 = (-4)*IT_1021;
    const ccomplex_t IT_1023 = IT_0001*IT_0127*IT_0389;
    const ccomplex_t IT_1024 = IT_0020*IT_0127;
    const ccomplex_t IT_1025 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_1026 = 0.5*IT_1023 + (-0.5)*IT_1024 + -IT_1025;
    const ccomplex_t IT_1027 = IT_0387 + IT_1026;
    const ccomplex_t IT_1028 = IT_0002*IT_1027;
    const ccomplex_t IT_1029 = -IT_0390;
    const ccomplex_t IT_1030 = IT_0020 + IT_1028 + IT_1029;
    const ccomplex_t IT_1031 = m_N_2*IT_1030;
    const ccomplex_t IT_1032 = (-2)*IT_0809 + (-2)*IT_0811 + 4*IT_0916;
    const ccomplex_t IT_1033 = IT_0808 + IT_1032;
    const ccomplex_t IT_1034 = IT_0263*IT_1033;
    const ccomplex_t IT_1035 = 0.25*IT_1034;
    const ccomplex_t IT_1036 = (-0.5)*IT_0679;
    const ccomplex_t IT_1037 = (-0.5)*IT_0529;
    const ccomplex_t IT_1038 = s_12 + IT_1035 + IT_1036 + IT_1037;
    const ccomplex_t IT_1039 = 8*IT_0954;
    const ccomplex_t IT_1040 = IT_0263*IT_0686;
    const ccomplex_t IT_1041 = 0.25*IT_1040;
    const ccomplex_t IT_1042 = 1 + IT_1041;
    const ccomplex_t IT_1043 = IT_0000*IT_1042;
    const ccomplex_t IT_1044 = 8*IT_0430;
    const ccomplex_t IT_1045 = 8*IT_0435;
    const ccomplex_t IT_1046 = (-8)*IT_0448;
    const ccomplex_t IT_1047 = (-0.25)*IT_0312 + 0.25*IT_0315 + (-0.25)*IT_0898;
    const ccomplex_t IT_1048 = s_23 + IT_1047;
    const ccomplex_t IT_1049 = m_C_2*IT_1048;
    const ccomplex_t IT_1050 = (-4)*IT_1049;
    const ccomplex_t IT_1051 = (-2)*IT_0692;
    const ccomplex_t IT_1052 = (-2)*IT_0526 + 2*IT_0771 + -IT_0925;
    const ccomplex_t IT_1053 = IT_0525 + IT_1052;
    const ccomplex_t IT_1054 = IT_0002*IT_1053;
    const ccomplex_t IT_1055 = -IT_0301 + (-0.5)*IT_1054;
    const ccomplex_t IT_1056 = s_12 + IT_1055;
    const ccomplex_t IT_1057 = m_N_2*IT_1056;
    const ccomplex_t IT_1058 = (-2)*IT_1057;
    const ccomplex_t IT_1059 = -IT_0925;
    const ccomplex_t IT_1060 = IT_0525 + IT_1059;
    const ccomplex_t IT_1061 = IT_0002*IT_1060;
    const ccomplex_t IT_1062 = -IT_0301 + (-0.25)*IT_1061;
    const ccomplex_t IT_1063 = s_12 + IT_1062;
    const ccomplex_t IT_1064 = m_C_2*IT_1063;
    const ccomplex_t IT_1065 = (-4)*IT_1064;
    const ccomplex_t IT_1066 = (-2)*IT_0295;
    const ccomplex_t IT_1067 = m_C_2*IT_0633;
    const ccomplex_t IT_1068 = (-4)*IT_1067;
    const ccomplex_t IT_1069 = (-2)*IT_0594;
    const ccomplex_t IT_1070 = IT_0002*IT_0046*IT_0387;
    const ccomplex_t IT_1071 = IT_0020*IT_0046;
    const ccomplex_t IT_1072 = -IT_1071;
    const ccomplex_t IT_1073 = IT_1070 + IT_1072;
    const ccomplex_t IT_1074 = IT_0001*IT_1073;
    const ccomplex_t IT_1075 = IT_0388 + IT_0715;
    const ccomplex_t IT_1076 = (-4)*IT_1075;
    const ccomplex_t IT_1077 = IT_1074 + IT_1076;
    const ccomplex_t IT_1078 = m_N_2*IT_1077;
    const ccomplex_t IT_1079 = (16*IT_0835 + 16*IT_0837 + (-16)*IT_0959)
      *IT_0963 + 8*IT_0923*IT_0968 + IT_0834*IT_0971 + IT_0976*IT_0979 + IT_0987
      *IT_0989 + IT_0532*IT_0998 + IT_0590*IT_1004 + IT_1012*IT_1013 + IT_0695
      *IT_1022 + IT_0699*IT_1031 + IT_1038*IT_1039 + IT_1043*(IT_1044 + IT_1045 
      + IT_1046) + IT_1050*IT_1051 + IT_0861*IT_1058 + IT_1065*IT_1066 + IT_1068
      *IT_1069 + IT_0945*IT_1078;
    const ccomplex_t IT_1080 = 2*IT_0256;
    const ccomplex_t IT_1081 = IT_0582 + IT_1080;
    const ccomplex_t IT_1082 = 8*IT_0959;
    const ccomplex_t IT_1083 = (-8)*IT_0835;
    const ccomplex_t IT_1084 = (-8)*IT_0837;
    const ccomplex_t IT_1085 = IT_1082 + IT_1083 + IT_1084;
    const ccomplex_t IT_1086 = 16*IT_0448;
    const ccomplex_t IT_1087 = (-16)*IT_0430;
    const ccomplex_t IT_1088 = (-16)*IT_0435;
    const ccomplex_t IT_1089 = IT_1086 + IT_1087 + IT_1088;
    const ccomplex_t IT_1090 = (-2)*IT_0187;
    const ccomplex_t IT_1091 = IT_0869 + IT_1090;
    const ccomplex_t IT_1092 = m_N_2*IT_0625;
    const ccomplex_t IT_1093 = (-4)*IT_1092;
    const ccomplex_t IT_1094 = m_N_2*IT_1063;
    const ccomplex_t IT_1095 = (-4)*IT_1094;
    const ccomplex_t IT_1096 = s_13*s_14*s_34*IT_0002;
    const ccomplex_t IT_1097 = 2*IT_0389 + -IT_1071 + (-2)*IT_1096;
    const ccomplex_t IT_1098 = IT_1070 + IT_1097;
    const ccomplex_t IT_1099 = IT_0001*IT_1098;
    const ccomplex_t IT_1100 = -IT_0020 + (-0.5)*IT_1099;
    const ccomplex_t IT_1101 = IT_0388 + IT_1100;
    const ccomplex_t IT_1102 = m_N_2*IT_1101;
    const ccomplex_t IT_1103 = (-2)*IT_1102;
    const ccomplex_t IT_1104 = (-0.25)*IT_0404 + 0.25*IT_0406 + (-0.25)*IT_1014;
    const ccomplex_t IT_1105 = s_13 + IT_1104;
    const ccomplex_t IT_1106 = m_N_2*IT_1105;
    const ccomplex_t IT_1107 = (-4)*IT_1106;
    const ccomplex_t IT_1108 = IT_0263*IT_0813;
    const ccomplex_t IT_1109 = 0.25*IT_1108;
    const ccomplex_t IT_1110 = s_12 + IT_1036 + IT_1037 + IT_1109;
    const ccomplex_t IT_1111 = (-4)*IT_0295;
    const ccomplex_t IT_1112 = -IT_1024;
    const ccomplex_t IT_1113 = IT_1023 + IT_1112;
    const ccomplex_t IT_1114 = IT_0002*IT_1113;
    const ccomplex_t IT_1115 = (-4)*IT_0716;
    const ccomplex_t IT_1116 = IT_1114 + IT_1115;
    const ccomplex_t IT_1117 = m_N_2*IT_1116;
    const ccomplex_t IT_1118 = (-2)*IT_0617;
    const ccomplex_t IT_1119 = IT_0590*IT_0932 + IT_0883*IT_0979 + IT_0906
      *IT_1051 + IT_0895*IT_1069 + IT_0941*IT_1081 + IT_1043*IT_1085 + IT_0963
      *IT_1089 + IT_0848*IT_1091 + IT_0532*IT_1093 + IT_0540*IT_1095 + IT_0868
      *IT_1103 + IT_0695*IT_1107 + (IT_0964 + IT_0966)*IT_1110 + IT_0915*IT_1111
       + IT_0861*IT_1117 + IT_0897*IT_1118;
    const ccomplex_t IT_1120 = 8*IT_0082;
    const ccomplex_t IT_1121 = (-8)*IT_0060;
    const ccomplex_t IT_1122 = IT_1120 + IT_1121;
    const ccomplex_t IT_1123 = IT_1044 + IT_1046;
    const ccomplex_t IT_1124 = (-2)*IT_0512;
    const ccomplex_t IT_1125 = (-2)*IT_0559;
    const ccomplex_t IT_1126 = IT_0665*IT_0998 + IT_0540*IT_1004 + IT_0868
      *IT_1012 + IT_0693*IT_1022 + IT_0713*IT_1031 + IT_0964*IT_1043 + IT_0861
      *IT_1065 + IT_1058*IT_1066 + IT_0834*IT_1081 + IT_1038*IT_1085 + IT_0987
      *IT_1091 + IT_1078*IT_1118 + IT_0976*IT_1122 + IT_0923*IT_1123 + IT_1068
      *IT_1124 + IT_1050*IT_1125;
    const ccomplex_t IT_1127 = -IT_0050 + IT_0090;
    const ccomplex_t IT_1128 = 8*IT_0835;
    const ccomplex_t IT_1129 = 8*IT_0837;
    const ccomplex_t IT_1130 = IT_1128 + IT_1129;
    const ccomplex_t IT_1131 = IT_0859*IT_0941 + (-16)*IT_0448*IT_0963 +
       IT_0848*IT_0969 + IT_0861*IT_1095 + IT_0945*IT_1103 + IT_0540*IT_1117 +
       IT_1093*IT_1124 + IT_1107*IT_1125 + 8*IT_0883*IT_1127 + IT_1043*IT_1130;
    const ccomplex_t IT_1132 = 2*conj(IT_0256);
    const ccomplex_t IT_1133 = IT_0584 + IT_1132;
    const ccomplex_t IT_1134 = 8*conj(IT_0082);
    const ccomplex_t IT_1135 = (-8)*conj(IT_0060);
    const ccomplex_t IT_1136 = IT_1134 + IT_1135;
    const ccomplex_t IT_1137 = 8*conj(IT_0430);
    const ccomplex_t IT_1138 = (-8)*conj(IT_0448);
    const ccomplex_t IT_1139 = IT_1137 + IT_1138;
    const ccomplex_t IT_1140 = (-8)*conj(IT_0835);
    const ccomplex_t IT_1141 = (-8)*conj(IT_0837);
    const ccomplex_t IT_1142 = IT_1140 + IT_1141;
    const ccomplex_t IT_1143 = (-2)*conj(IT_0187);
    const ccomplex_t IT_1144 = IT_0942 + IT_1143;
    const ccomplex_t IT_1145 = (-8)*conj(IT_0478);
    const ccomplex_t IT_1146 = (-2)*conj(IT_0512);
    const ccomplex_t IT_1147 = (-2)*conj(IT_0295);
    const ccomplex_t IT_1148 = (-2)*conj(IT_0617);
    const ccomplex_t IT_1149 = (-2)*conj(IT_0559);
    const ccomplex_t IT_1150 = IT_0666*IT_0998 + IT_0541*IT_1004 + IT_0682
      *IT_1012 + IT_0694*IT_1022 + IT_0714*IT_1031 + IT_0933*IT_1065 + IT_0834
      *IT_1133 + IT_0976*IT_1136 + IT_0923*IT_1139 + IT_1038*IT_1142 + IT_0987
      *IT_1144 + IT_1043*IT_1145 + IT_1068*IT_1146 + IT_1058*IT_1147 + IT_1078
      *IT_1148 + IT_1050*IT_1149;
    const ccomplex_t IT_1151 = m_N_2*IT_0001;
    const ccomplex_t IT_1152 = IT_0002*IT_0419;
    const ccomplex_t IT_1153 = (-2)*IT_0408 + (-0.5)*IT_1152;
    const ccomplex_t IT_1154 = IT_0407 + IT_1153;
    const ccomplex_t IT_1155 = IT_1151*IT_1154;
    const ccomplex_t IT_1156 = (-2)*IT_1155;
    const ccomplex_t IT_1157 = 2*IT_0835;
    const ccomplex_t IT_1158 = 2*IT_0837;
    const ccomplex_t IT_1159 = IT_1157 + IT_1158;
    const ccomplex_t IT_1160 = (-2)*IT_0959;
    const ccomplex_t IT_1161 = IT_1159 + IT_1160;
    const ccomplex_t IT_1162 = 2*conj(IT_0835);
    const ccomplex_t IT_1163 = 2*conj(IT_0837);
    const ccomplex_t IT_1164 = IT_1162 + IT_1163;
    const ccomplex_t IT_1165 = (-2)*conj(IT_0959);
    const ccomplex_t IT_1166 = IT_1164 + IT_1165;
    const ccomplex_t IT_1167 = (-2)*IT_0954;
    const ccomplex_t IT_1168 = (-2)*conj(IT_0954);
    const ccomplex_t IT_1169 = s_23*IT_0046*IT_0127;
    const ccomplex_t IT_1170 = s_24*s_34*IT_0127;
    const ccomplex_t IT_1171 = (-2)*IT_1170;
    const ccomplex_t IT_1172 = IT_1169 + IT_1171;
    const ccomplex_t IT_1173 = IT_0263*IT_1172;
    const ccomplex_t IT_1174 = (-0.666666666666667)*IT_0536 + (
      -0.166666666666667)*IT_1173;
    const ccomplex_t IT_1175 = s_23 + IT_1174;
    const ccomplex_t IT_1176 = m_C_2*IT_1175;
    const ccomplex_t IT_1177 = (-6)*IT_1176;
    const ccomplex_t IT_1178 = IT_0522*IT_0824;
    const ccomplex_t IT_1179 = 2*IT_0478;
    const ccomplex_t IT_1180 = -IT_0524;
    const ccomplex_t IT_1181 = IT_0127 + IT_1180;
    const ccomplex_t IT_1182 = IT_0000*IT_1181;
    const ccomplex_t IT_1183 = (-6)*IT_0512;
    const ccomplex_t IT_1184 = conj(IT_0594)*(IT_0293*IT_0518 + IT_0418
      *IT_0523 + IT_0531*IT_0665 + IT_0582*IT_0674 + IT_0825*IT_0969 + IT_0839
      *IT_0998 + IT_0449*IT_1093 + IT_0868*IT_1178 + IT_0895*IT_1179 + IT_1182
      *IT_1183);
    const ccomplex_t IT_1185 = (-2)*conj(IT_0835);
    const ccomplex_t IT_1186 = (-2)*conj(IT_0837);
    const ccomplex_t IT_1187 = IT_1185 + IT_1186;
    const ccomplex_t IT_1188 = 2*conj(IT_0165);
    const ccomplex_t IT_1189 = 2*conj(IT_0478);
    const ccomplex_t IT_1190 = m_C_2*m_N_2*IT_0127;
    const ccomplex_t IT_1191 = IT_0006*IT_1190;
    const ccomplex_t IT_1192 = s_13*m_N_2;
    const ccomplex_t IT_1193 = IT_0006*IT_1192;
    const ccomplex_t IT_1194 = s_23*m_C_2;
    const ccomplex_t IT_1195 = IT_0006*IT_1194;
    const ccomplex_t IT_1196 = (-2)*IT_0778;
    const ccomplex_t IT_1197 = IT_0774 + IT_1196;
    const ccomplex_t IT_1198 = IT_0002*IT_1197;
    const ccomplex_t IT_1199 = 0.5*IT_1198;
    const ccomplex_t IT_1200 = s_14*s_23*IT_0127;
    const ccomplex_t IT_1201 = s_13*s_24*IT_0127;
    const ccomplex_t IT_1202 = IT_1200 + IT_1201;
    const ccomplex_t IT_1203 = s_12*s_34*IT_0127;
    const ccomplex_t IT_1204 = -IT_1203;
    const ccomplex_t IT_1205 = IT_1202 + IT_1204;
    const ccomplex_t IT_1206 = IT_0313*IT_1205;
    const ccomplex_t IT_1207 = 0.5*IT_1206;
    const ccomplex_t IT_1208 = (-2)*IT_0809 + 2*IT_0810 + 2*IT_0916 + (-2)
      *IT_0917;
    const ccomplex_t IT_1209 = IT_0808 + IT_1208;
    const ccomplex_t IT_1210 = IT_0001*IT_1209;
    const ccomplex_t IT_1211 = (-0.5)*IT_1210;
    const ccomplex_t IT_1212 = IT_0525 + IT_0527 + IT_1199 + IT_1207 + IT_1211;
    const ccomplex_t IT_1213 = IT_0551 + IT_0553 + IT_0555 + IT_0668 + IT_0669
       + IT_0671 + IT_0723;
    const ccomplex_t IT_1214 = 2*conj(IT_0187);
    const ccomplex_t IT_1215 = IT_1188 + IT_1214;
    const ccomplex_t IT_1216 = IT_0499 + IT_0503 + IT_0505 + IT_0507 + IT_0544
       + IT_0546 + IT_0700;
    const ccomplex_t IT_1217 = 2*conj(IT_0617);
    const ccomplex_t IT_1218 = s_34*IT_0001;
    const ccomplex_t IT_1219 = IT_0875 + IT_0878;
    const ccomplex_t IT_1220 = s_13*s_23*s_34*IT_0002;
    const ccomplex_t IT_1221 = (-2)*IT_1220;
    const ccomplex_t IT_1222 = IT_1219 + IT_1221;
    const ccomplex_t IT_1223 = IT_1218*IT_1222;
    const ccomplex_t IT_1224 = (-2)*conj(IT_0594);
    const ccomplex_t IT_1225 = (-0.25)*IT_0557*IT_0592 + (-0.25)*IT_0454
      *IT_1107 + (-0.25)*IT_0682*IT_1178 + (-0.25)*IT_1022*IT_1187 + (-0.25)
      *IT_0906*IT_1189 + conj(IT_0559)*IT_1191 + -IT_0088*IT_1193 + -IT_0769
      *IT_1195 + -IT_0262*IT_1212 + (-0.25)*IT_1133*IT_1213 + (-0.25)*IT_0825
      *IT_1215 + (-0.25)*IT_1216*IT_1217 + (-0.25)*IT_1223*IT_1224;
    const ccomplex_t IT_1226 = IT_0263*IT_1169;
    const ccomplex_t IT_1227 = 0.666666666666667*IT_0312 + (-0.166666666666667
      )*IT_1226;
    const ccomplex_t IT_1228 = s_23 + IT_1227;
    const ccomplex_t IT_1229 = m_C_2*IT_1228;
    const ccomplex_t IT_1230 = (-6)*IT_1229;
    const ccomplex_t IT_1231 = s_12 + IT_0922;
    const ccomplex_t IT_1232 = IT_0000*IT_1231;
    const ccomplex_t IT_1233 = conj(IT_0835) + conj(IT_0837);
    const ccomplex_t IT_1234 = IT_0000*IT_0394;
    const ccomplex_t IT_1235 = 2*conj(IT_0448);
    const ccomplex_t IT_1236 = IT_0000*IT_0309;
    const ccomplex_t IT_1237 = (-2)*conj(IT_0396);
    const ccomplex_t IT_1238 = IT_0290*IT_0310 + IT_0307*IT_0395 + IT_0304
      *IT_0402 + 4*conj(IT_0478)*IT_0915 + IT_0484*IT_0932 + IT_0860*IT_0942 +
       IT_0454*IT_1095 + IT_1058*IT_1164 + IT_0840*IT_1188 + (-12)*conj(IT_0284)
      *IT_1232 + (-4)*IT_1031*IT_1233 + IT_0682*IT_1234 + IT_1117*IT_1235 +
       IT_1236*IT_1237;
    const ccomplex_t IT_1239 = (-10)*s_23;
    const ccomplex_t IT_1240 = IT_1226 + IT_1239;
    const ccomplex_t IT_1241 = m_C_2*IT_1240;
    const ccomplex_t IT_1242 = 0.125*IT_0266;
    const ccomplex_t IT_1243 = IT_0127 + IT_1242;
    const ccomplex_t IT_1244 = IT_0000*IT_1243;
    const ccomplex_t IT_1245 = 8*IT_1244;
    const ccomplex_t IT_1246 = (-2)*IT_0190;
    const ccomplex_t IT_1247 = (-2)*IT_0191;
    const ccomplex_t IT_1248 = IT_1246 + IT_1247;
    const ccomplex_t IT_1249 = (-2)*conj(IT_0190);
    const ccomplex_t IT_1250 = (-2)*conj(IT_0191);
    const ccomplex_t IT_1251 = IT_1249 + IT_1250;
    const ccomplex_t IT_1252 = (-4)*IT_0007;
    const ccomplex_t IT_1253 = -IT_0060;
    const ccomplex_t IT_1254 = -conj(IT_0060);
    const ccomplex_t IT_1255 = -IT_0090;
    const ccomplex_t IT_1256 = -conj(IT_0090);
    const ccomplex_t IT_1257 = -IT_0050;
    const ccomplex_t IT_1258 = IT_0835 + IT_0837;
    const ccomplex_t IT_1259 = -conj(IT_0050);
    const ccomplex_t IT_1260 = -IT_0478;
    const ccomplex_t IT_1261 = conj(IT_0060)*IT_0478 + IT_0060*conj(IT_0478) +
       IT_0082*(conj(IT_0457) + conj(IT_0462) + -conj(IT_0478) + conj(IT_0954)) 
      + (conj(IT_0457) + conj(IT_0462) + conj(IT_0954))*IT_1253 + (IT_0457 +
       IT_0462 + IT_0954)*IT_1254 + (conj(IT_0430) + conj(IT_0435) + conj
      (IT_0959))*(IT_0050 + IT_1255) + (IT_0430 + IT_0435 + IT_0959)*(conj
      (IT_0050) + IT_1256) + (conj(IT_0448) + IT_1233)*(IT_0090 + IT_1257) + 
      (IT_0448 + IT_1258)*(conj(IT_0090) + IT_1259) + conj(IT_0082)*(IT_0457 +
       IT_0462 + IT_0954 + IT_1260);
    const ccomplex_t IT_1262 = 0.166666666666667*IT_0282*IT_0293 +
       0.166666666666667*IT_0287*IT_0310 + 0.166666666666667*IT_0304*IT_0399 +
       0.166666666666667*IT_0395*IT_0418 + 0.166666666666667*IT_0698*IT_0713 +
       0.166666666666667*IT_0860*IT_0869 + 0.166666666666667*IT_0479*IT_0932 +
       0.166666666666667*IT_0840*IT_0969 + 0.166666666666667*IT_0866*IT_1066 +
       0.166666666666667*IT_0449*IT_1095 + 0.166666666666667*IT_0849*IT_1117 +
       0.166666666666667*IT_1058*IT_1159 + 0.166666666666667*IT_0759*IT_1183 +
       IT_0284*(IT_0717 + IT_0720 + (-2)*IT_1232) + 0.166666666666667*IT_0868
      *IT_1234 + 0.166666666666667*IT_0945*IT_1236 + (-0.666666666666667)
      *IT_1031*IT_1258 + (-0.666666666666667)*IT_0915*IT_1260;
    const ccomplex_t IT_1263 = conj(IT_0090) + IT_1259;
    const ccomplex_t IT_1264 = 0.125*IT_1117;
    const ccomplex_t IT_1265 = 0.125*IT_1095;
    const ccomplex_t IT_1266 = 0.125*IT_1093;
    const ccomplex_t IT_1267 = 0.125*IT_1103;
    const ccomplex_t IT_1268 = (-2)*conj(IT_0233);
    const ccomplex_t IT_1269 = 0.125*IT_0941;
    const ccomplex_t IT_1270 = 0.125*IT_0848;
    const ccomplex_t IT_1271 = 0.125*IT_1107;
    const ccomplex_t IT_1272 = conj(IT_0448)*IT_0963 + (-0.5)*(conj(IT_0837) +
       0.125*IT_0924)*IT_1043 + (-0.5)*conj(IT_0478)*IT_1110 + (-0.5)*IT_0883
      *IT_1263 + (-0.5)*IT_0541*IT_1264 + (-0.5)*IT_0933*IT_1265 + (-0.5)
      *IT_1146*IT_1266 + (-0.5)*IT_1237*IT_1267 + (-0.5)*IT_1268*IT_1269 + (-0.5
      )*IT_1188*IT_1270 + (-0.5)*IT_1149*IT_1271;
    const ccomplex_t IT_1273 = 8*IT_0448;
    const ccomplex_t IT_1274 = 0.125*IT_0932;
    const ccomplex_t IT_1275 = 0.125*IT_0915;
    const ccomplex_t IT_1276 = 0.125*IT_0895;
    const ccomplex_t IT_1277 = 0.125*IT_0897;
    const ccomplex_t IT_1278 = 0.125*IT_0906;
    const ccomplex_t IT_1279 = 8*IT_0060;
    const ccomplex_t IT_1280 = (-8)*IT_0082;
    const ccomplex_t IT_1281 = IT_1279 + IT_1280;
    const ccomplex_t IT_1282 = 0.125*IT_0883;
    const ccomplex_t IT_1283 = IT_0923*(IT_0837 + 0.125*IT_1128) + IT_0869
      *IT_1269 + 0.125*IT_1110*IT_1273 + IT_0861*IT_1274 + IT_0699*IT_1275 +
       IT_0532*IT_1276 + IT_0868*IT_1277 + IT_0695*IT_1278 + IT_1281*IT_1282;
    const ccomplex_t IT_1284 = 8*conj(IT_0050);
    const ccomplex_t IT_1285 = (-8)*conj(IT_0090);
    const ccomplex_t IT_1286 = IT_1284 + IT_1285;
    const ccomplex_t IT_1287 = (-2)*conj(IT_0256);
    const ccomplex_t IT_1288 = IT_1268 + IT_1287;
    const ccomplex_t IT_1289 = (-2)*conj(IT_0692);
    const ccomplex_t IT_1290 = 8*conj(IT_0457);
    const ccomplex_t IT_1291 = IT_1145 + IT_1290;
    const ccomplex_t IT_1292 = 0.125*IT_0923;
    const ccomplex_t IT_1293 = (-0.0625)*IT_0513*IT_0998 + (-0.0625)*IT_0592
      *IT_1004 + (-0.0625)*IT_0696*IT_1022 + (-0.0625)*IT_0712*IT_1031 + (
      -0.0625)*IT_0933*IT_1058 + (-0.5)*IT_1043*(conj(IT_0435) + 0.125*IT_1137 +
       0.125*IT_1138) + (-0.0625)*IT_1065*IT_1147 + (-0.0625)*IT_0834*IT_1215 + 
      (-0.0625)*IT_1012*IT_1217 + (-0.0625)*IT_1068*IT_1224 + IT_0963*(conj
      (IT_0959) + -IT_1233) + (-0.0625)*IT_1078*IT_1237 + (-0.0625)*IT_0976
      *IT_1286 + (-0.0625)*IT_0987*IT_1288 + (-0.0625)*IT_1050*IT_1289 + (-0.5)
      *IT_1291*IT_1292;
    const ccomplex_t IT_1294 = 0.125*IT_1043*IT_1085 + 0.125*IT_0963*IT_1089 +
       IT_0968*IT_1110 + IT_0861*IT_1264 + IT_0540*IT_1265 + IT_0532*IT_1266 +
       IT_0868*IT_1267 + IT_1081*IT_1269 + IT_1091*IT_1270 + IT_0695*IT_1271 +
       IT_0590*IT_1274 + IT_1111*IT_1275 + IT_1069*IT_1276 + IT_1118*IT_1277 +
       IT_1051*IT_1278 + IT_0979*IT_1282 + IT_1039*IT_1292;
    const ccomplex_t IT_1295 = 0.166666666666667*IT_0667;
    const ccomplex_t IT_1296 = 0.166666666666667*IT_0860;
    const ccomplex_t IT_1297 = 0.166666666666667*IT_0840;
    const ccomplex_t IT_1298 = 0.166666666666667*IT_0825;
    const ccomplex_t IT_1299 = 0.166666666666667*IT_0681;
    const ccomplex_t IT_1300 = 0.166666666666667*IT_0858;
    const ccomplex_t IT_1301 = 0.166666666666667*IT_0987;
    const ccomplex_t IT_1302 = 0.166666666666667*IT_0941;
    const ccomplex_t IT_1303 = 0.166666666666667*IT_0848;
    const ccomplex_t IT_1304 = 0.166666666666667*IT_1213;
    const ccomplex_t IT_1305 = 0.166666666666667*IT_0674;
    const ccomplex_t IT_1306 = IT_0261*IT_0851 + IT_0293*IT_1295 + IT_1066
      *IT_1296 + IT_0590*IT_1297 + IT_0665*IT_1298 + IT_0868*IT_1299 + IT_1091
      *IT_1300 + IT_1159*IT_1301 + IT_0449*IT_1302 + IT_1179*IT_1303 + IT_0695
      *IT_1304 + IT_0532*IT_1305;
    const ccomplex_t IT_1307 = conj(IT_0233)*IT_0851 + IT_0298*IT_1295 +
       IT_1147*IT_1296 + IT_0592*IT_1297 + IT_0666*IT_1298 + IT_0682*IT_1299 +
       IT_0942*IT_1300 + IT_1164*IT_1301 + IT_0454*IT_1302 + IT_1189*IT_1303 +
       IT_0696*IT_1304 + IT_0513*IT_1305;
    const ccomplex_t IT_1308 = IT_0834*IT_0839 + 6*IT_0851*(IT_0187 +
       0.166666666666667*IT_0852) + IT_0855*IT_0856 + 6*IT_0287*IT_1295 + 6
      *IT_0861*IT_1296 + 6*IT_0540*IT_1297 + 6*IT_0826*IT_1298 + 6*IT_1013
      *IT_1299 + 6*IT_0989*IT_1300 + 6*IT_0479*IT_1302 + 6*IT_0849*IT_1303 + 6
      *IT_0665*IT_1305;
    const ccomplex_t IT_1309 = s_13*s_14*s_24;
    const ccomplex_t IT_1310 = s_24*s_34*IT_0020;
    const ccomplex_t IT_1311 = -IT_0795 + (-2)*IT_1309 + 2*IT_1310;
    const ccomplex_t IT_1312 = IT_0792 + IT_1311;
    const ccomplex_t IT_1313 = IT_0001*IT_1312;
    const ccomplex_t IT_1314 = (-0.5)*IT_1313;
    const ccomplex_t IT_1315 = IT_0791 + IT_1314;
    const ccomplex_t IT_1316 = 0.25*IT_0664;
    const ccomplex_t IT_1317 = 0.25*IT_0691;
    const ccomplex_t IT_1318 = (-0.666666666666667)*IT_0320 +
       0.166666666666667*IT_1173;
    const ccomplex_t IT_1319 = s_23 + IT_1318;
    const ccomplex_t IT_1320 = m_C_2*IT_1319;
    const ccomplex_t IT_1321 = 6*IT_1320;
    const ccomplex_t IT_1322 = 0.25*IT_1321;
    const ccomplex_t IT_1323 = IT_0263*IT_0419;
    const ccomplex_t IT_1324 = (-10)*s_13;
    const ccomplex_t IT_1325 = IT_1323 + IT_1324;
    const ccomplex_t IT_1326 = m_N_2*IT_1325;
    const ccomplex_t IT_1327 = 0.25*IT_1326;
    const ccomplex_t IT_1328 = m_C_2*IT_0001;
    const ccomplex_t IT_1329 = IT_0002*IT_1169;
    const ccomplex_t IT_1330 = (-2)*IT_0317 + (-0.5)*IT_1329;
    const ccomplex_t IT_1331 = IT_0316 + IT_1330;
    const ccomplex_t IT_1332 = IT_1328*IT_1331;
    const ccomplex_t IT_1333 = (-2)*IT_1332;
    const ccomplex_t IT_1334 = 0.25*IT_1333;
    const ccomplex_t IT_1335 = (-0.166666666666667)*IT_0423 + (
      -0.666666666666667)*IT_0533;
    const ccomplex_t IT_1336 = s_13 + IT_1335;
    const ccomplex_t IT_1337 = m_N_2*IT_1336;
    const ccomplex_t IT_1338 = (-6)*IT_1337;
    const ccomplex_t IT_1339 = 0.25*IT_1338;
    const ccomplex_t IT_1340 = -IT_0771;
    const ccomplex_t IT_1341 = IT_0525 + IT_0773 + IT_0783 + IT_0784 + IT_1340;
    const ccomplex_t IT_1342 = 0.25*IT_1341;
    const ccomplex_t IT_1343 = (-8)*conj(IT_0190);
    const ccomplex_t IT_1344 = 0.25*IT_0748;
    const ccomplex_t IT_1345 = (-0.5)*IT_0700;
    const ccomplex_t IT_1346 = (-0.25)*IT_0709;
    const ccomplex_t IT_1347 = (-0.5)*IT_0502;
    const ccomplex_t IT_1348 = (-0.5)*IT_0543;
    const ccomplex_t IT_1349 = 0.5*IT_0545;
    const ccomplex_t IT_1350 = IT_0500 + IT_1345 + IT_1346 + IT_1347 + IT_1348
       + IT_1349;
    const ccomplex_t IT_1351 = 0.25*IT_1350;
    const ccomplex_t IT_1352 = (-4)*conj(IT_0165);
    const ccomplex_t IT_1353 = (-4)*conj(IT_0187);
    const ccomplex_t IT_1354 = IT_1352 + IT_1353;
    const ccomplex_t IT_1355 = s_13*IT_0276;
    const ccomplex_t IT_1356 = s_12*s_24*s_34;
    const ccomplex_t IT_1357 = -IT_0103 + 2*IT_1355 + (-2)*IT_1356;
    const ccomplex_t IT_1358 = IT_0101 + IT_1357;
    const ccomplex_t IT_1359 = IT_0001*IT_1358;
    const ccomplex_t IT_1360 = (-0.5)*IT_1359;
    const ccomplex_t IT_1361 = IT_0100 + IT_1360;
    const ccomplex_t IT_1362 = 0.25*IT_1361;
    const ccomplex_t IT_1363 = (-4)*conj(IT_0284);
    const ccomplex_t IT_1364 = -IT_0723;
    const ccomplex_t IT_1365 = (-0.5)*IT_0732;
    const ccomplex_t IT_1366 = IT_0098 + IT_0551 + IT_0552 + IT_1364 + IT_1365;
    const ccomplex_t IT_1367 = 0.25*IT_1366;
    const ccomplex_t IT_1368 = conj(IT_0559)*IT_1212 + conj(IT_0396)*IT_1315 +
       IT_0666*IT_1316 + IT_0694*IT_1317 + IT_1166*IT_1322 + IT_1168*IT_1327 +
       IT_0455*IT_1334 + IT_0485*IT_1339 + IT_0788*IT_1342 + IT_1343*IT_1344 +
       IT_0739*IT_1351 + IT_1354*IT_1362 + IT_1363*IT_1367;
    const ccomplex_t IT_1369 = 0.666666666666667*IT_0404 + (-0.166666666666667
      )*IT_1323;
    const ccomplex_t IT_1370 = s_13 + IT_1369;
    const ccomplex_t IT_1371 = m_N_2*IT_1370;
    const ccomplex_t IT_1372 = (-6)*IT_1371;
    const ccomplex_t IT_1373 = (-2)*conj(IT_0189);
    const ccomplex_t IT_1374 = (-8)*conj(IT_0258);
    const ccomplex_t IT_1375 = (-8)*conj(IT_0255);
    const ccomplex_t IT_1376 = (-8)*conj(IT_0259);
    const ccomplex_t IT_1377 = (-4)*conj(IT_0233);
    const ccomplex_t IT_1378 = (-4)*conj(IT_0256);
    const ccomplex_t IT_1379 = (-4)*conj(IT_0295);
    const ccomplex_t IT_1380 = conj(IT_0692)*IT_1212 + 0.25*IT_0497*IT_1245 +
       conj(IT_0617)*IT_1315 + IT_0513*IT_1316 + IT_0696*IT_1317 + IT_1168
      *IT_1322 + IT_1164*IT_1327 + IT_0485*IT_1334 + IT_0455*IT_1339 + IT_0790
      *IT_1342 + IT_0741*IT_1351 + 0.25*IT_1165*IT_1372 + 0.25*IT_0491*(IT_1251 
      + IT_1373) + IT_1344*(IT_1374 + IT_1375 + IT_1376) + IT_1362*(IT_1377 +
       IT_1378) + IT_1367*IT_1379;
    const ccomplex_t IT_1381 = (1 + 0.5*IT_0004)*(IT_0525 + (-2)*IT_0526);
    const ccomplex_t IT_1382 = (-4)*IT_0559;
    const ccomplex_t IT_1383 = 0.25*IT_1223;
    const ccomplex_t IT_1384 = 0.25*IT_1216;
    const ccomplex_t IT_1385 = 0.25*IT_1213;
    const ccomplex_t IT_1386 = (-0.25)*IT_0557*IT_0590 + (-0.25)*IT_0825
      *IT_0971 + (-0.25)*IT_0839*IT_1022 + (-0.25)*IT_0449*IT_1107 + (-0.25)
      *IT_0868*IT_1178 + (-0.25)*IT_0906*IT_1179 + -IT_0085*IT_1193 + -(IT_0050 
      + 0.25*IT_0770)*IT_1195 + -IT_0260*IT_1212 + IT_0692*IT_1381 + (-0.25)
      *IT_1191*IT_1382 + -IT_1069*IT_1383 + -IT_1013*IT_1384 + -IT_1081*IT_1385;
    const ccomplex_t IT_1387 = IT_0001*IT_1033;
    const ccomplex_t IT_1388 = -IT_1387;
    const ccomplex_t IT_1389 = IT_0807 + IT_1388;
    const ccomplex_t IT_1390 = -conj(IT_0396);
    const ccomplex_t IT_1391 = (-0.25)*IT_0664;
    const ccomplex_t IT_1392 = (-0.25)*IT_0691;
    const ccomplex_t IT_1393 = (-0.25)*IT_1321;
    const ccomplex_t IT_1394 = (-0.25)*IT_1326;
    const ccomplex_t IT_1395 = (-0.25)*IT_1333;
    const ccomplex_t IT_1396 = (-0.25)*IT_1338;
    const ccomplex_t IT_1397 = (-0.25)*IT_1341;
    const ccomplex_t IT_1398 = (-0.25)*IT_1350;
    const ccomplex_t IT_1399 = (-0.25)*IT_1361;
    const ccomplex_t IT_1400 = (-0.25)*IT_1366;
    const ccomplex_t IT_1401 = -conj(IT_0559)*IT_1212 + 2*IT_0748*(IT_0749 + (
      -0.125)*IT_1343) + (conj(IT_0343) + conj(IT_0355))*IT_1389 + IT_1315
      *IT_1390 + IT_0666*IT_1391 + IT_0694*IT_1392 + IT_1166*IT_1393 + IT_1168
      *IT_1394 + IT_0455*IT_1395 + IT_0485*IT_1396 + IT_0788*IT_1397 + IT_0739
      *IT_1398 + IT_1354*IT_1399 + IT_1363*IT_1400;
    const ccomplex_t IT_1402 = (-8)*IT_0258;
    const ccomplex_t IT_1403 = (-8)*IT_0255;
    const ccomplex_t IT_1404 = (-8)*IT_0259;
    const ccomplex_t IT_1405 = (-4)*IT_0374;
    const ccomplex_t IT_1406 = (-2)*IT_0189;
    const ccomplex_t IT_1407 = (-4)*IT_0233;
    const ccomplex_t IT_1408 = (-4)*IT_0256;
    const ccomplex_t IT_1409 = IT_0692*IT_1212 + 0.25*IT_0494*IT_1245 +
       IT_0617*IT_1315 + 0.25*IT_1160*IT_1372 + -IT_0532*IT_1391 + -IT_0695
      *IT_1392 + -IT_1167*IT_1393 + -IT_1159*IT_1394 + -IT_0480*IT_1395 + 
      -IT_0450*IT_1396 + -IT_0789*IT_1397 + -IT_0740*IT_1398 + -IT_1111*IT_1400 
      + 0.25*IT_0748*(IT_1402 + IT_1403 + IT_1404) + -IT_1389*(IT_0384 + (-0.25)
      *IT_1405) + 0.25*IT_0491*(IT_1248 + IT_1406) + -IT_1399*(IT_1407 + IT_1408);
    const ccomplex_t IT_1410 = IT_0000*IT_1075;
    const ccomplex_t IT_1411 = (-0.166666666666667)*IT_0395;
    const ccomplex_t IT_1412 = (-0.166666666666667)*IT_1234;
    const ccomplex_t IT_1413 = (-0.166666666666667)*IT_1236;
    const ccomplex_t IT_1414 = (-0.166666666666667)*IT_1178;
    const ccomplex_t IT_1415 = (-0.166666666666667)*IT_0681;
    const ccomplex_t IT_1416 = (-0.166666666666667)*IT_0897;
    const ccomplex_t IT_1417 = (-0.166666666666667)*IT_1103;
    const ccomplex_t IT_1418 = (-0.166666666666667)*IT_1078;
    const ccomplex_t IT_1419 = (-0.166666666666667)*IT_0510;
    const ccomplex_t IT_1420 = 0.166666666666667*IT_0399*IT_0498 + -IT_0194
      *IT_0855 + 0.166666666666667*IT_0944*IT_1118 + IT_0396*IT_1410 + -IT_0418
      *IT_1411 + -IT_0540*IT_1412 + -IT_0861*IT_1413 + -IT_0826*IT_1414 + 
      -IT_0582*IT_1415 + -IT_1179*IT_1416 + -IT_0449*IT_1417 + -IT_1159*IT_1418 
      + -IT_0665*IT_1419;
    const ccomplex_t IT_1421 = IT_0513 + IT_0696;
    const ccomplex_t IT_1422 = (-0.166666666666667)*IT_0402*IT_0498 + IT_0188
      *IT_0855 + IT_1390*IT_1410 + IT_0307*IT_1411 + IT_0541*IT_1412 + IT_0933
      *IT_1413 + IT_0584*IT_1415 + IT_1189*IT_1416 + IT_0454*IT_1417 + IT_1164
      *IT_1418 + IT_0666*IT_1419 + IT_1414*IT_1421;
    const ccomplex_t IT_1423 = IT_0665*IT_1316;
    const ccomplex_t IT_1424 = IT_0693*IT_1317;
    const ccomplex_t IT_1425 = IT_1161*IT_1322;
    const ccomplex_t IT_1426 = IT_1167*IT_1327;
    const ccomplex_t IT_1427 = IT_0450*IT_1334;
    const ccomplex_t IT_1428 = IT_0480*IT_1339;
    const ccomplex_t IT_1429 = (-4)*IT_0355;
    const ccomplex_t IT_1430 = IT_1389*IT_1429;
    const ccomplex_t IT_1431 = 0.25*IT_1430;
    const ccomplex_t IT_1432 = IT_0787*IT_1342;
    const ccomplex_t IT_1433 = (-8)*IT_0190;
    const ccomplex_t IT_1434 = IT_0738*IT_1351;
    const ccomplex_t IT_1435 = (-4)*IT_0165;
    const ccomplex_t IT_1436 = (-4)*IT_0187;
    const ccomplex_t IT_1437 = IT_1435 + IT_1436;
    const ccomplex_t IT_1438 = IT_1362*IT_1437;
    const ccomplex_t IT_1439 = (-4)*IT_0284;
    const ccomplex_t IT_1440 = IT_1367*IT_1439;
    const ccomplex_t IT_1441 = IT_0396*IT_1315;
    const ccomplex_t IT_1442 = IT_0559*IT_1212;
    const ccomplex_t IT_1443 = (-0.5)*IT_1423 + (-0.5)*IT_1424 + (-0.5)
      *IT_1425 + (-0.5)*IT_1426 + (-0.5)*IT_1427 + (-0.5)*IT_1428 + (-0.5)
      *IT_1431 + (-0.5)*IT_1432 + IT_0748*(IT_0751 + (-0.125)*IT_1433) + (-0.5)
      *IT_1434 + (-0.5)*IT_1438 + (-0.5)*IT_1440 + (-0.5)*IT_1441 + (-0.5)
      *IT_1442;
    const ccomplex_t IT_1444 = IT_0290*IT_0667;
    const ccomplex_t IT_1445 = IT_0860*IT_0933;
    const ccomplex_t IT_1446 = IT_0541*IT_0840;
    const ccomplex_t IT_1447 = IT_0825*IT_1421;
    const ccomplex_t IT_1448 = (-6)*conj(IT_0396);
    const ccomplex_t IT_1449 = IT_0855*IT_1448;
    const ccomplex_t IT_1450 = 6*conj(IT_0165);
    const ccomplex_t IT_1451 = IT_0851*IT_1450;
    const ccomplex_t IT_1452 = IT_0848*IT_1235;
    const ccomplex_t IT_1453 = IT_0834*IT_1187;
    const ccomplex_t IT_1454 = IT_0478*IT_0963;
    const ccomplex_t IT_1455 = 16*IT_1454;
    const ccomplex_t IT_1456 = IT_1066*IT_1117;
    const ccomplex_t IT_1457 = IT_0915*IT_1439;
    const ccomplex_t IT_1458 = IT_0590*IT_1095;
    const ccomplex_t IT_1459 = IT_0895*IT_1124;
    const ccomplex_t IT_1460 = IT_0665*IT_1093;
    const ccomplex_t IT_1461 = IT_0897*IT_0945;
    const ccomplex_t IT_1462 = IT_1013*IT_1103;
    const ccomplex_t IT_1463 = IT_0941*IT_0971;
    const ccomplex_t IT_1464 = IT_0906*IT_1125;
    const ccomplex_t IT_1465 = IT_0540*IT_0932 + IT_0848*IT_0989 + IT_0923*
      (IT_1083 + IT_1084) + IT_0693*IT_1107 + IT_0883*IT_1122 + IT_1110*IT_1123 
      + IT_1455 + IT_1456 + IT_1457 + IT_1458 + IT_1459 + IT_1460 + IT_1461 +
       IT_1462 + IT_1463 + IT_1464;
    const ccomplex_t IT_1466 = IT_0540*IT_0932 + IT_0848*IT_0989 + IT_0923
      *IT_1085 + IT_0693*IT_1107 + (IT_1044 + IT_1045 + IT_1046)*IT_1110 +
       IT_0883*IT_1122 + IT_1455 + IT_1456 + IT_1457 + IT_1458 + IT_1459 +
       IT_1460 + IT_1461 + IT_1462 + IT_1463 + IT_1464;
    const ccomplex_t IT_1467 = IT_0695*IT_1050;
    const ccomplex_t IT_1468 = IT_0969*IT_0987;
    const ccomplex_t IT_1469 = IT_0540*IT_1065;
    const ccomplex_t IT_1470 = IT_0532*IT_1068;
    const ccomplex_t IT_1471 = IT_1467 + IT_1468 + IT_1469 + IT_1470;
    const ccomplex_t IT_1472 = IT_0861*IT_1004;
    const ccomplex_t IT_1473 = IT_0945*IT_1012;
    const ccomplex_t IT_1474 = IT_0976*IT_1281;
    const ccomplex_t IT_1475 = IT_0696*IT_1050;
    const ccomplex_t IT_1476 = IT_0987*IT_1188;
    const ccomplex_t IT_1477 = IT_0541*IT_1065;
    const ccomplex_t IT_1478 = IT_0513*IT_1068;
    const ccomplex_t IT_1479 = IT_1475 + IT_1476 + IT_1477 + IT_1478;
    const ccomplex_t IT_1480 = IT_0933*IT_1004;
    const ccomplex_t IT_1481 = IT_1012*IT_1237;
    const ccomplex_t IT_1482 = IT_0886*IT_0976;
    const ccomplex_t IT_1483 = 8*conj(IT_0959);
    const ccomplex_t IT_1484 = IT_1140 + IT_1141 + IT_1483;
    const ccomplex_t IT_1485 = 16*conj(IT_0478);
    const ccomplex_t IT_1486 = IT_0541*IT_0932;
    const ccomplex_t IT_1487 = IT_1117*IT_1147;
    const ccomplex_t IT_1488 = IT_0915*IT_1363;
    const ccomplex_t IT_1489 = IT_0592*IT_1095;
    const ccomplex_t IT_1490 = IT_0895*IT_1146;
    const ccomplex_t IT_1491 = IT_0666*IT_1093;
    const ccomplex_t IT_1492 = IT_0897*IT_1237;
    const ccomplex_t IT_1493 = IT_1103*IT_1217;
    const ccomplex_t IT_1494 = IT_0941*IT_1215;
    const ccomplex_t IT_1495 = IT_0848*IT_1288;
    const ccomplex_t IT_1496 = IT_0906*IT_1149;
    const ccomplex_t IT_1497 = IT_0694*IT_1107;
    const ccomplex_t IT_1498 = IT_0883*IT_1136;
    const ccomplex_t IT_1499 = IT_0592*IT_0932;
    const ccomplex_t IT_1500 = IT_0933*IT_1117;
    const ccomplex_t IT_1501 = IT_0915*IT_1379;
    const ccomplex_t IT_1502 = IT_0541*IT_1095;
    const ccomplex_t IT_1503 = IT_0895*IT_1224;
    const ccomplex_t IT_1504 = IT_0513*IT_1093;
    const ccomplex_t IT_1505 = IT_0897*IT_1148;
    const ccomplex_t IT_1506 = IT_0682*IT_1103;
    const ccomplex_t IT_1507 = IT_0941*IT_1133;
    const ccomplex_t IT_1508 = IT_0848*IT_1144;
    const ccomplex_t IT_1509 = IT_0906*IT_1289;
    const ccomplex_t IT_1510 = IT_0696*IT_1107;
    const ccomplex_t IT_1511 = IT_1043*IT_1484;
    const ccomplex_t IT_1512 = 16*conj(IT_0448);
    const ccomplex_t IT_1513 = (-16)*conj(IT_0430);
    const ccomplex_t IT_1514 = (-16)*conj(IT_0435);
    const ccomplex_t IT_1515 = IT_1512 + IT_1513 + IT_1514;
    const ccomplex_t IT_1516 = IT_0963*IT_1515;
    const ccomplex_t IT_1517 = IT_0883*IT_1286;
    const ccomplex_t IT_1518 = 4*IT_0007*(conj(IT_0050)*IT_0085 + IT_0050
      *IT_0088 + conj(IT_0090)*IT_0093 + IT_0090*IT_0096) + 4*IT_0109*(IT_0188
      *IT_0193 + IT_0194*IT_0196 + IT_0257*IT_0260 + IT_0261*IT_0262) + (-2)
      *IT_0271*(conj(IT_0189)*IT_0272 + IT_0189*IT_0273) + IT_0282*((conj
      (IT_0165) + conj(IT_0187) + conj(IT_0284))*IT_0287 + (IT_0165 + IT_0187 +
       IT_0284)*IT_0290 + IT_0257*IT_0293 + (IT_0233 + IT_0256 + IT_0295)
      *IT_0298) + IT_0284*(IT_0304*IT_0307 + IT_0298*IT_0310) + IT_0324*(IT_0287
      *IT_0356 + IT_0290*IT_0357 + IT_0293*IT_0385 + IT_0298*IT_0386) + IT_0395*
      ((conj(IT_0284) + conj(IT_0396))*IT_0399 + (IT_0284 + IT_0396)*IT_0402) +
       IT_0415*(IT_0307*IT_0386 + IT_0356*IT_0399 + IT_0357*IT_0402 + IT_0385
      *IT_0418) + IT_0427*(IT_0262*IT_0450 + IT_0260*IT_0455 + IT_0196*IT_0480 +
       IT_0193*IT_0485) + IT_0491*(IT_0262*IT_0494 + IT_0260*IT_0497) + IT_0396*
      (IT_0307*IT_0498 + IT_0510*IT_0513) + conj(IT_0512)*(IT_0287*IT_0518 +
       IT_0399*IT_0523 + IT_0531*IT_0532 + IT_0539*IT_0540) + IT_0512*(IT_0290
      *IT_0518 + IT_0402*IT_0523 + IT_0539*IT_0541) + IT_0557*(IT_0560 + IT_0561
      ) + (IT_0571 + IT_0581)*(IT_0583 + IT_0585 + IT_0586 + IT_0587 + IT_0588 +
       IT_0589 + IT_0591 + IT_0593) + IT_0539*(IT_0595 + IT_0596) + (IT_0606 +
       IT_0614)*(IT_0615 + IT_0616 + IT_0618 + IT_0619) + IT_0656 + IT_0664*
      (IT_0260*IT_0513 + IT_0262*IT_0532 + IT_0196*IT_0665 + IT_0193*IT_0666) +
       IT_0165*(IT_0298*IT_0667 + IT_0513*IT_0674 + IT_0681*IT_0682) + IT_0691*
      (IT_0196*IT_0693 + IT_0193*IT_0694 + IT_0262*IT_0695 + IT_0260*IT_0696) +
       conj(IT_0284)*(IT_0293*IT_0310 + IT_0304*IT_0418 + IT_0698*IT_0699) +
       IT_0711*(IT_0196*IT_0699 + IT_0193*IT_0712 + IT_0262*IT_0713 + IT_0260
      *IT_0714) + (IT_0717 + IT_0720)*IT_0722 + IT_0737*(IT_0196*IT_0738 +
       IT_0193*IT_0739 + IT_0262*IT_0740 + IT_0260*IT_0741) + IT_0748*(IT_0749
      *IT_0750 + IT_0751*IT_0752) + IT_0759*(IT_0761 + IT_0763 + IT_0765) + 4
      *IT_0766*(IT_0060*conj(IT_0060) + IT_0090*conj(IT_0090) + 0.25*conj
      (IT_0082)*IT_0091 + IT_0082*IT_0096 + IT_0050*IT_0769 + 0.25*conj(IT_0050)
      *IT_0770) + IT_0786*(IT_0196*IT_0787 + IT_0193*IT_0788 + IT_0262*IT_0789 +
       IT_0260*IT_0790) + IT_0800*(IT_0196*IT_0801 + IT_0193*IT_0802 + IT_0262
      *IT_0803 + IT_0260*IT_0804) + IT_0816*IT_0823 + conj(IT_0233)*(IT_0287
      *IT_0667 + IT_0825*IT_0826 + IT_0834*IT_0839 + IT_0540*IT_0840 + IT_0848
      *IT_0849 + IT_0851*IT_0852 + IT_0855*IT_0856 + IT_0858*IT_0859 + IT_0860
      *IT_0861) + IT_0866*IT_0867 + conj(IT_0165)*(IT_0293*IT_0667 + IT_0532
      *IT_0674 + IT_0681*IT_0868 + IT_0858*IT_0869) + IT_0448*(IT_0883*IT_0886 +
       IT_0513*IT_0895 + IT_0682*IT_0897 + IT_0696*IT_0906 + IT_0712*IT_0915 +
       IT_0923*IT_0924 + IT_0932*IT_0933 + IT_0941*IT_0942) + conj(IT_0396)*
      (IT_0418*IT_0498 + IT_0510*IT_0532 + IT_0944*IT_0945) + conj(IT_0954)
      *IT_1079 + conj(IT_0457)*IT_1119 + conj(IT_0959)*IT_1126 + conj(IT_0478)
      *IT_1131 + IT_0959*IT_1150 + IT_1156*(IT_0262*IT_1161 + IT_0260*IT_1166 +
       IT_0196*IT_1167 + IT_0193*IT_1168) + (IT_0196*IT_1161 + IT_0193*IT_1166 +
       IT_0262*IT_1167 + IT_0260*IT_1168)*IT_1177 + IT_1184 + IT_0594*(IT_0298
      *IT_0518 + IT_0307*IT_0523 + IT_0584*IT_0674 + IT_0454*IT_1093 + IT_0682
      *IT_1178 + IT_0760*IT_1182 + IT_0998*IT_1187 + IT_0825*IT_1188 + IT_0895
      *IT_1189) + (-4)*IT_0692*IT_1225 + (IT_0195*IT_0449 + IT_0192*IT_0454 +
       conj(IT_0258)*IT_0479 + IT_0258*IT_0484)*IT_1230 + IT_0295*IT_1238 + 
      (IT_0196*IT_0437 + conj(IT_0189)*IT_0449 + IT_0193*IT_0453 + IT_0189
      *IT_0454 + IT_0262*IT_0464 + IT_0273*IT_0479 + IT_0260*IT_0483 + IT_0272
      *IT_0484)*IT_1241 + (-2)*IT_1245*(conj(IT_0258)*(IT_0189 + (-0.5)*IT_1246 
      + (-0.5)*IT_1247) + -1./2*IT_0273*IT_1248 + IT_0258*(conj(IT_0189) + (-0.5
      )*IT_1249 + (-0.5)*IT_1250) + -1./2*IT_0272*IT_1251) + IT_1252*IT_1261 + 6
      *conj(IT_0295)*IT_1262 + (-16)*IT_0478*IT_1272 + 8*conj(IT_0448)*IT_1283 +
       (-16)*IT_0954*IT_1293 + 8*conj(IT_0462)*IT_1294 + 6*conj(IT_0187)*IT_1306
       + 6*IT_0187*IT_1307 + conj(IT_0256)*IT_1308 + 4*IT_0355*IT_1368 + 4*
      (IT_0374 + IT_0384)*IT_1380 + 4*conj(IT_0559)*(IT_1127*IT_1193 + IT_0093
      *IT_1195 + IT_0193*IT_1212 + 0.25*IT_1381*IT_1382 + IT_1124*IT_1383 +
       IT_0868*IT_1384 + IT_0969*IT_1385) + 4*IT_0559*(IT_0096*IT_1195 + IT_0196
      *IT_1212 + IT_1193*IT_1263 + IT_1146*IT_1383 + IT_0682*IT_1384 + IT_1188
      *IT_1385) + (-4)*conj(IT_0692)*IT_1386 + (-4)*IT_0343*IT_1401 + 4*(conj
      (IT_0374) + conj(IT_0384))*IT_1409 + 6*conj(IT_0617)*IT_1420 + (-6)
      *IT_0617*IT_1422 + 4*conj(IT_0355)*(IT_1423 + IT_1424 + IT_1425 + IT_1426 
      + IT_1427 + IT_1428 + IT_1431 + IT_1432 + IT_1344*IT_1433 + IT_1434 +
       IT_1438 + IT_1440 + IT_1441 + IT_1442) + (-8)*conj(IT_0343)*IT_1443 +
       IT_0233*(IT_1444 + IT_1445 + IT_1446 + IT_1447 + IT_1449 + IT_1451 +
       IT_1452 + IT_1453) + IT_0256*(IT_0666*IT_0674 + IT_0484*IT_0941 + IT_0681
      *IT_1217 + IT_0858*IT_1268 + IT_1444 + IT_1445 + IT_1446 + IT_1447 +
       IT_1449 + IT_1451 + IT_1452 + IT_1453) + conj(IT_0430)*IT_1465 + conj
      (IT_0435)*IT_1466 + conj(IT_0835)*(IT_1038*IT_1128 + IT_1471 + IT_1472 +
       IT_1473 + IT_1474) + conj(IT_0837)*(IT_1038*IT_1130 + IT_0923*IT_1273 +
       IT_1471 + IT_1472 + IT_1473 + IT_1474) + IT_0835*(IT_1479 + IT_1480 +
       IT_1481 + IT_1482) + IT_0837*(IT_0924*IT_1038 + IT_1479 + IT_1480 +
       IT_1481 + IT_1482) + IT_0430*(IT_1110*IT_1138 + IT_0923*IT_1142 + IT_0963
      *IT_1485 + IT_1486 + IT_1487 + IT_1488 + IT_1489 + IT_1490 + IT_1491 +
       IT_1492 + IT_1493 + IT_1494 + IT_1495 + IT_1496 + IT_1497 + IT_1498) +
       IT_0435*(IT_1110*IT_1139 + IT_0923*IT_1484 + IT_0963*IT_1485 + IT_1486 +
       IT_1487 + IT_1488 + IT_1489 + IT_1490 + IT_1491 + IT_1492 + IT_1493 +
       IT_1494 + IT_1495 + IT_1496 + IT_1497 + IT_1498) + IT_0457*(IT_1110
      *IT_1145 + IT_1499 + IT_1500 + IT_1501 + IT_1502 + IT_1503 + IT_1504 +
       IT_1505 + IT_1506 + IT_1507 + IT_1508 + IT_1509 + IT_1510 + IT_1511 +
       IT_1516 + IT_1517) + IT_0462*(IT_1110*IT_1291 + IT_1499 + IT_1500 +
       IT_1501 + IT_1502 + IT_1503 + IT_1504 + IT_1505 + IT_1506 + IT_1507 +
       IT_1508 + IT_1509 + IT_1510 + IT_1511 + IT_1516 + IT_1517);
    return create_ccomplex_return(IT_1518);
}

