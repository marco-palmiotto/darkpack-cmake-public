#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_4_to_Z_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_4_to_Z_W(
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
    const ccomplex_t IT_0000 = m_C_2*m_N_4;
    const ccomplex_t IT_0001 = pow(m_W, -2);
    const ccomplex_t IT_0002 = pow(m_Z, -2);
    const ccomplex_t IT_0003 = pow(s_34, 2);
    const ccomplex_t IT_0004 = IT_0001*IT_0002*IT_0003;
    const ccomplex_t IT_0005 = 0.5*IT_0004;
    const ccomplex_t IT_0006 = 1 + IT_0005;
    const ccomplex_t IT_0007 = IT_0000*IT_0006;
    const ccomplex_t IT_0008 = pow(m_W, 2);
    const ccomplex_t IT_0009 = pow(m_C_2, 2);
    const ccomplex_t IT_0010 = pow(m_N_4, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_12 + IT_0008 + -IT_0009 + -IT_0010 
      + -reg_prop, -1);
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = sin(theta_W);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = IT_0012*IT_0014;
    const ccomplex_t IT_0016 = conj(N_d4)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = (-2)*IT_0017;
    const ccomplex_t IT_0019 = cos(theta_W);
    const ccomplex_t IT_0020 = cpow(IT_0019, -1);
    const ccomplex_t IT_0021 = IT_0012*IT_0020;
    const ccomplex_t IT_0022 = conj(N_B4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0023 = IT_0021*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W4)*conj(U_d2)*e_em;
    const ccomplex_t IT_0026 = IT_0015*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = IT_0024 + IT_0027;
    const ccomplex_t IT_0029 = IT_0018 + IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1)*IT_0029;
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*m_W*e_em*IT_0013*IT_0020;
    const ccomplex_t IT_0033 = -IT_0032;
    const ccomplex_t IT_0034 = IT_0031*IT_0033;
    const ccomplex_t IT_0035 = IT_0011*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*e_em*IT_0014*IT_0019;
    const ccomplex_t IT_0038 = m_C_2*IT_0037;
    const ccomplex_t IT_0039 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0040 = IT_0014*IT_0039;
    const ccomplex_t IT_0041 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0042 = IT_0014*IT_0041;
    const ccomplex_t IT_0043 = 1.4142135623731*IT_0042;
    const ccomplex_t IT_0044 = (0 + _Complex_I*1)*(IT_0040 + (-0.5)*IT_0043);
    const ccomplex_t IT_0045 = -IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0009 + -IT_0010 + -reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = IT_0038*IT_0048;
    const ccomplex_t IT_0050 = IT_0036 + -IT_0049;
    const ccomplex_t IT_0051 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0052 = IT_0014*IT_0051;
    const ccomplex_t IT_0053 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0054 = IT_0014*IT_0053;
    const ccomplex_t IT_0055 = 1.4142135623731*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*(IT_0052 + 0.5*IT_0055);
    const ccomplex_t IT_0057 = -IT_0056;
    const ccomplex_t IT_0058 = IT_0046*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = IT_0038*IT_0059;
    const ccomplex_t IT_0061 = sin(beta);
    const ccomplex_t IT_0062 = IT_0014*IT_0061;
    const ccomplex_t IT_0063 = N_u4*e_em*V_Wp2;
    const ccomplex_t IT_0064 = IT_0062*IT_0063;
    const ccomplex_t IT_0065 = 2*IT_0064;
    const ccomplex_t IT_0066 = IT_0020*IT_0061;
    const ccomplex_t IT_0067 = N_B4*V_u2*e_em;
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = N_W4*V_u2*e_em;
    const ccomplex_t IT_0071 = IT_0062*IT_0070;
    const ccomplex_t IT_0072 = 1.4142135623731*IT_0071;
    const ccomplex_t IT_0073 = IT_0069 + IT_0072;
    const ccomplex_t IT_0074 = IT_0065 + IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*IT_0074;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0033*IT_0076;
    const ccomplex_t IT_0078 = IT_0011*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = m_N_4*IT_0037;
    const ccomplex_t IT_0081 = IT_0048*IT_0080;
    const ccomplex_t IT_0082 = -IT_0079 + IT_0081;
    const ccomplex_t IT_0083 = (-4)*IT_0082;
    const ccomplex_t IT_0084 = IT_0060 + 0.25*IT_0083;
    const ccomplex_t IT_0085 = (-4)*conj(IT_0082);
    const ccomplex_t IT_0086 = conj(IT_0060) + 0.25*IT_0085;
    const ccomplex_t IT_0087 = IT_0059*IT_0080;
    const ccomplex_t IT_0088 = -IT_0087;
    const ccomplex_t IT_0089 = (-4)*IT_0060;
    const ccomplex_t IT_0090 = 0.25*IT_0089;
    const ccomplex_t IT_0091 = IT_0082 + IT_0090;
    const ccomplex_t IT_0092 = (-4)*conj(IT_0060);
    const ccomplex_t IT_0093 = 0.25*IT_0092;
    const ccomplex_t IT_0094 = conj(IT_0082) + IT_0093;
    const ccomplex_t IT_0095 = pow(m_Z, 2);
    const ccomplex_t IT_0096 = s_23*IT_0095;
    const ccomplex_t IT_0097 = IT_0002*IT_0096;
    const ccomplex_t IT_0098 = s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0099 = s_24*IT_0095;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = s_23*IT_0008;
    const ccomplex_t IT_0102 = s_24*s_34;
    const ccomplex_t IT_0103 = (-2)*IT_0102;
    const ccomplex_t IT_0104 = IT_0101 + IT_0103;
    const ccomplex_t IT_0105 = IT_0001*IT_0104;
    const ccomplex_t IT_0106 = 0.5*IT_0097 + (-0.5)*IT_0100 + (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = s_23 + IT_0106;
    const ccomplex_t IT_0108 = m_C_2*IT_0107;
    const ccomplex_t IT_0109 = 2*IT_0108;
    const ccomplex_t IT_0110 = IT_0013*IT_0020;
    const ccomplex_t IT_0111 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = IT_0014*IT_0019;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0116 = IT_0110*IT_0115;
    const ccomplex_t IT_0117 = IT_0113*IT_0115;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*(IT_0112 + IT_0114 + 
      -IT_0116 + -IT_0117);
    const ccomplex_t IT_0119 = (-0.5)*IT_0118;
    const ccomplex_t IT_0120 = cpow((-2)*s_23 + IT_0095 + (0 + _Complex_I*1)
      *m_N_4*Gamma_o4 + reg_prop, -1);
    const ccomplex_t IT_0121 = IT_0057*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = IT_0119*IT_0122;
    const ccomplex_t IT_0124 = 0.5*IT_0123;
    const ccomplex_t IT_0125 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0126 = IT_0014*IT_0125;
    const ccomplex_t IT_0127 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0128 = IT_0014*IT_0127;
    const ccomplex_t IT_0129 = 1.4142135623731*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*(IT_0126 + 0.5*IT_0129);
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0133 = IT_0110*IT_0132;
    const ccomplex_t IT_0134 = IT_0113*IT_0132;
    const ccomplex_t IT_0135 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0136 = IT_0113*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1)*(IT_0133 + -IT_0134 + (-2)
      *IT_0136);
    const ccomplex_t IT_0138 = 0.5*IT_0137;
    const ccomplex_t IT_0139 = IT_0131*IT_0138;
    const ccomplex_t IT_0140 = cpow((-2)*s_13 + IT_0009 + IT_0095 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0141 = IT_0139*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0145 = IT_0110*IT_0144;
    const ccomplex_t IT_0146 = IT_0113*IT_0144;
    const ccomplex_t IT_0147 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0148 = IT_0113*IT_0147;
    const ccomplex_t IT_0149 = (0 + _Complex_I*1)*(IT_0145 + -IT_0146 + (-2)
      *IT_0148);
    const ccomplex_t IT_0150 = 0.5*IT_0149;
    const ccomplex_t IT_0151 = IT_0057*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_13 + IT_0095 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = -IT_0154;
    const ccomplex_t IT_0156 = conj(N_d3)*N_d4*e_em;
    const ccomplex_t IT_0157 = IT_0110*IT_0156;
    const ccomplex_t IT_0158 = IT_0113*IT_0156;
    const ccomplex_t IT_0159 = conj(N_u3)*N_u4*e_em;
    const ccomplex_t IT_0160 = IT_0110*IT_0159;
    const ccomplex_t IT_0161 = IT_0113*IT_0159;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0157 + IT_0158 + 
      -IT_0160 + -IT_0161);
    const ccomplex_t IT_0163 = (-0.5)*IT_0162;
    const ccomplex_t IT_0164 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0165 = IT_0014*IT_0164;
    const ccomplex_t IT_0166 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0167 = IT_0014*IT_0166;
    const ccomplex_t IT_0168 = 1.4142135623731*IT_0167;
    const ccomplex_t IT_0169 = (0 + _Complex_I*1)*(IT_0165 + 0.5*IT_0168);
    const ccomplex_t IT_0170 = -IT_0169;
    const ccomplex_t IT_0171 = cpow((-2)*s_23 + IT_0010 + IT_0095 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0172 = IT_0170*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*IT_0172;
    const ccomplex_t IT_0174 = IT_0163*IT_0173;
    const ccomplex_t IT_0175 = 0.5*IT_0118;
    const ccomplex_t IT_0176 = IT_0122*IT_0175;
    const ccomplex_t IT_0177 = conj(N_d1)*N_d4*e_em;
    const ccomplex_t IT_0178 = IT_0110*IT_0177;
    const ccomplex_t IT_0179 = IT_0113*IT_0177;
    const ccomplex_t IT_0180 = conj(N_u1)*N_u4*e_em;
    const ccomplex_t IT_0181 = IT_0110*IT_0180;
    const ccomplex_t IT_0182 = IT_0113*IT_0180;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0178 + IT_0179 + 
      -IT_0181 + -IT_0182);
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0186 = IT_0014*IT_0185;
    const ccomplex_t IT_0187 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0188 = IT_0014*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*(IT_0186 + 0.5*IT_0189);
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = cpow((-2)*s_23 + IT_0010 + IT_0095 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0193 = IT_0191*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0184*IT_0194;
    const ccomplex_t IT_0196 = conj(N_d2)*N_d4*e_em;
    const ccomplex_t IT_0197 = IT_0110*IT_0196;
    const ccomplex_t IT_0198 = IT_0113*IT_0196;
    const ccomplex_t IT_0199 = conj(N_u2)*N_u4*e_em;
    const ccomplex_t IT_0200 = IT_0110*IT_0199;
    const ccomplex_t IT_0201 = IT_0113*IT_0199;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*(IT_0197 + IT_0198 + 
      -IT_0200 + -IT_0201);
    const ccomplex_t IT_0203 = (-0.5)*IT_0202;
    const ccomplex_t IT_0204 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0205 = IT_0014*IT_0204;
    const ccomplex_t IT_0206 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0207 = IT_0014*IT_0206;
    const ccomplex_t IT_0208 = 1.4142135623731*IT_0207;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0205 + 0.5*IT_0208);
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = cpow((-2)*s_23 + IT_0010 + IT_0095 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0212 = IT_0210*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*IT_0212;
    const ccomplex_t IT_0214 = IT_0203*IT_0213;
    const ccomplex_t IT_0215 = -IT_0174 + (-0.5)*IT_0176 + -IT_0195 + -IT_0214;
    const ccomplex_t IT_0216 = conj(IT_0124) + conj(IT_0143) + conj(IT_0155) +
       conj(IT_0215);
    const ccomplex_t IT_0217 = 2*IT_0060;
    const ccomplex_t IT_0218 = (-2)*IT_0082;
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = IT_0124 + IT_0143 + IT_0155 + IT_0215;
    const ccomplex_t IT_0221 = 2*conj(IT_0060);
    const ccomplex_t IT_0222 = (-2)*conj(IT_0082);
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0225 = IT_0014*IT_0224;
    const ccomplex_t IT_0226 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0227 = IT_0014*IT_0226;
    const ccomplex_t IT_0228 = 1.4142135623731*IT_0227;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*(IT_0225 + (-0.5)*IT_0228);
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0232 = IT_0110*IT_0231;
    const ccomplex_t IT_0233 = IT_0113*IT_0231;
    const ccomplex_t IT_0234 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0235 = IT_0113*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*(IT_0232 + -IT_0233 + (-2)
      *IT_0235);
    const ccomplex_t IT_0237 = (-0.5)*IT_0236;
    const ccomplex_t IT_0238 = IT_0230*IT_0237;
    const ccomplex_t IT_0239 = IT_0140*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*IT_0239;
    const ccomplex_t IT_0241 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0242 = IT_0110*IT_0241;
    const ccomplex_t IT_0243 = IT_0113*IT_0241;
    const ccomplex_t IT_0244 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0245 = IT_0113*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + -IT_0243 + (-2)
      *IT_0245);
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0045*IT_0247;
    const ccomplex_t IT_0249 = IT_0152*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*IT_0249;
    const ccomplex_t IT_0251 = N_u3*conj(N_u4)*e_em;
    const ccomplex_t IT_0252 = IT_0110*IT_0251;
    const ccomplex_t IT_0253 = IT_0113*IT_0251;
    const ccomplex_t IT_0254 = N_d3*conj(N_d4)*e_em;
    const ccomplex_t IT_0255 = IT_0110*IT_0254;
    const ccomplex_t IT_0256 = IT_0113*IT_0254;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*(IT_0252 + IT_0253 + 
      -IT_0255 + -IT_0256);
    const ccomplex_t IT_0258 = (-0.5)*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0260 = IT_0014*IT_0259;
    const ccomplex_t IT_0261 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0262 = IT_0014*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0260 + (-0.5)*IT_0263);
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = IT_0171*IT_0265;
    const ccomplex_t IT_0267 = (0 + _Complex_I*1)*IT_0266;
    const ccomplex_t IT_0268 = IT_0258*IT_0267;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = N_u2*conj(N_u4)*e_em;
    const ccomplex_t IT_0271 = IT_0110*IT_0270;
    const ccomplex_t IT_0272 = IT_0113*IT_0270;
    const ccomplex_t IT_0273 = N_d2*conj(N_d4)*e_em;
    const ccomplex_t IT_0274 = IT_0110*IT_0273;
    const ccomplex_t IT_0275 = IT_0113*IT_0273;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*(IT_0271 + IT_0272 + 
      -IT_0274 + -IT_0275);
    const ccomplex_t IT_0277 = (-0.5)*IT_0276;
    const ccomplex_t IT_0278 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0279 = IT_0014*IT_0278;
    const ccomplex_t IT_0280 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0281 = IT_0014*IT_0280;
    const ccomplex_t IT_0282 = 1.4142135623731*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*(IT_0279 + (-0.5)*IT_0282);
    const ccomplex_t IT_0284 = -IT_0283;
    const ccomplex_t IT_0285 = IT_0211*IT_0284;
    const ccomplex_t IT_0286 = (0 + _Complex_I*1)*IT_0285;
    const ccomplex_t IT_0287 = IT_0277*IT_0286;
    const ccomplex_t IT_0288 = IT_0045*IT_0120;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*IT_0288;
    const ccomplex_t IT_0290 = IT_0175*IT_0289;
    const ccomplex_t IT_0291 = -IT_0287 + IT_0290;
    const ccomplex_t IT_0292 = N_u1*conj(N_u4)*e_em;
    const ccomplex_t IT_0293 = IT_0110*IT_0292;
    const ccomplex_t IT_0294 = IT_0113*IT_0292;
    const ccomplex_t IT_0295 = N_d1*conj(N_d4)*e_em;
    const ccomplex_t IT_0296 = IT_0110*IT_0295;
    const ccomplex_t IT_0297 = IT_0113*IT_0295;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*(IT_0293 + IT_0294 + 
      -IT_0296 + -IT_0297);
    const ccomplex_t IT_0299 = (-0.5)*IT_0298;
    const ccomplex_t IT_0300 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0301 = IT_0014*IT_0300;
    const ccomplex_t IT_0302 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0303 = IT_0014*IT_0302;
    const ccomplex_t IT_0304 = 1.4142135623731*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*(IT_0301 + (-0.5)*IT_0304);
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0192*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = IT_0299*IT_0308;
    const ccomplex_t IT_0310 = -IT_0309;
    const ccomplex_t IT_0311 = conj(IT_0269) + conj(IT_0291) + conj(IT_0310);
    const ccomplex_t IT_0312 = conj(IT_0240) + conj(IT_0250) + IT_0311;
    const ccomplex_t IT_0313 = 2*IT_0088;
    const ccomplex_t IT_0314 = (-2)*IT_0050;
    const ccomplex_t IT_0315 = IT_0313 + IT_0314;
    const ccomplex_t IT_0316 = IT_0269 + IT_0291 + IT_0310;
    const ccomplex_t IT_0317 = IT_0240 + IT_0250 + IT_0316;
    const ccomplex_t IT_0318 = 2*conj(IT_0088);
    const ccomplex_t IT_0319 = (-2)*conj(IT_0050);
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = pow(s_13, 2);
    const ccomplex_t IT_0322 = IT_0002*IT_0321;
    const ccomplex_t IT_0323 = pow(s_14, 2);
    const ccomplex_t IT_0324 = IT_0001*IT_0323;
    const ccomplex_t IT_0325 = IT_0322 + IT_0324;
    const ccomplex_t IT_0326 = s_13*s_14*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0327 = -IT_0009 + -IT_0326;
    const ccomplex_t IT_0328 = IT_0325 + IT_0327;
    const ccomplex_t IT_0329 = m_N_4*IT_0328;
    const ccomplex_t IT_0330 = (-2)*IT_0037;
    const ccomplex_t IT_0331 = IT_0048*IT_0330;
    const ccomplex_t IT_0332 = (-2)*IT_0240 + (-2)*IT_0250 + -IT_0331;
    const ccomplex_t IT_0333 = 2*IT_0082;
    const ccomplex_t IT_0334 = (-2)*IT_0060;
    const ccomplex_t IT_0335 = IT_0333 + IT_0334;
    const ccomplex_t IT_0336 = IT_0037*IT_0048;
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = 2*conj(IT_0082);
    const ccomplex_t IT_0339 = (-2)*conj(IT_0060);
    const ccomplex_t IT_0340 = IT_0338 + IT_0339;
    const ccomplex_t IT_0341 = IT_0037*IT_0059;
    const ccomplex_t IT_0342 = -IT_0341;
    const ccomplex_t IT_0343 = IT_0059*IT_0330;
    const ccomplex_t IT_0344 = 2*IT_0142 + -IT_0343;
    const ccomplex_t IT_0345 = 2*IT_0154;
    const ccomplex_t IT_0346 = 2*IT_0050;
    const ccomplex_t IT_0347 = (-2)*IT_0088;
    const ccomplex_t IT_0348 = IT_0346 + IT_0347;
    const ccomplex_t IT_0349 = 2*conj(IT_0050);
    const ccomplex_t IT_0350 = (-2)*conj(IT_0088);
    const ccomplex_t IT_0351 = IT_0349 + IT_0350;
    const ccomplex_t IT_0352 = s_13*IT_0095;
    const ccomplex_t IT_0353 = IT_0002*IT_0352;
    const ccomplex_t IT_0354 = s_14*IT_0095;
    const ccomplex_t IT_0355 = IT_0098*IT_0354;
    const ccomplex_t IT_0356 = s_13*IT_0008;
    const ccomplex_t IT_0357 = s_14*s_34;
    const ccomplex_t IT_0358 = (-2)*IT_0357;
    const ccomplex_t IT_0359 = IT_0356 + IT_0358;
    const ccomplex_t IT_0360 = IT_0001*IT_0359;
    const ccomplex_t IT_0361 = 0.5*IT_0353 + (-0.5)*IT_0355 + (-0.5)*IT_0360;
    const ccomplex_t IT_0362 = s_13 + IT_0361;
    const ccomplex_t IT_0363 = m_N_4*IT_0362;
    const ccomplex_t IT_0364 = 2*IT_0363;
    const ccomplex_t IT_0365 = IT_0001*IT_0002;
    const ccomplex_t IT_0366 = s_13*IT_0008*IT_0095;
    const ccomplex_t IT_0367 = s_14*s_34*IT_0095;
    const ccomplex_t IT_0368 = (-2)*IT_0367;
    const ccomplex_t IT_0369 = IT_0366 + IT_0368;
    const ccomplex_t IT_0370 = IT_0365*IT_0369;
    const ccomplex_t IT_0371 = (-0.666666666666667)*IT_0360 +
       0.166666666666667*IT_0370;
    const ccomplex_t IT_0372 = s_13 + IT_0371;
    const ccomplex_t IT_0373 = m_N_4*IT_0372;
    const ccomplex_t IT_0374 = 6*IT_0373;
    const ccomplex_t IT_0375 = m_N_2*IT_0203;
    const ccomplex_t IT_0376 = m_N_4*IT_0277;
    const ccomplex_t IT_0377 = -IT_0376;
    const ccomplex_t IT_0378 = IT_0375 + IT_0377;
    const ccomplex_t IT_0379 = IT_0286*IT_0378;
    const ccomplex_t IT_0380 = m_N_3*IT_0163;
    const ccomplex_t IT_0381 = m_N_4*IT_0258;
    const ccomplex_t IT_0382 = -IT_0381;
    const ccomplex_t IT_0383 = IT_0380 + IT_0382;
    const ccomplex_t IT_0384 = IT_0267*IT_0383;
    const ccomplex_t IT_0385 = m_N_4*IT_0119;
    const ccomplex_t IT_0386 = m_N_4*IT_0175;
    const ccomplex_t IT_0387 = -IT_0386;
    const ccomplex_t IT_0388 = IT_0385 + IT_0387;
    const ccomplex_t IT_0389 = IT_0289*IT_0388;
    const ccomplex_t IT_0390 = -IT_0379 + -IT_0384 + IT_0389;
    const ccomplex_t IT_0391 = 2*IT_0390;
    const ccomplex_t IT_0392 = m_N_1*IT_0184;
    const ccomplex_t IT_0393 = m_N_4*IT_0299;
    const ccomplex_t IT_0394 = IT_0308*(IT_0392 + -IT_0393);
    const ccomplex_t IT_0395 = -IT_0394;
    const ccomplex_t IT_0396 = 2*IT_0395;
    const ccomplex_t IT_0397 = IT_0391 + IT_0396;
    const ccomplex_t IT_0398 = 2*conj(IT_0390);
    const ccomplex_t IT_0399 = 2*conj(IT_0395);
    const ccomplex_t IT_0400 = IT_0398 + IT_0399;
    const ccomplex_t IT_0401 = conj(IT_0124) + conj(IT_0215);
    const ccomplex_t IT_0402 = m_N_2*IT_0277;
    const ccomplex_t IT_0403 = m_N_4*IT_0203;
    const ccomplex_t IT_0404 = -IT_0403;
    const ccomplex_t IT_0405 = IT_0402 + IT_0404;
    const ccomplex_t IT_0406 = IT_0213*IT_0405;
    const ccomplex_t IT_0407 = m_N_3*IT_0258;
    const ccomplex_t IT_0408 = m_N_4*IT_0163;
    const ccomplex_t IT_0409 = -IT_0408;
    const ccomplex_t IT_0410 = IT_0407 + IT_0409;
    const ccomplex_t IT_0411 = IT_0173*IT_0410;
    const ccomplex_t IT_0412 = IT_0122*IT_0388;
    const ccomplex_t IT_0413 = -IT_0406 + -IT_0411 + (-0.5)*IT_0412;
    const ccomplex_t IT_0414 = 2*IT_0413;
    const ccomplex_t IT_0415 = m_N_1*IT_0299;
    const ccomplex_t IT_0416 = m_N_4*IT_0184;
    const ccomplex_t IT_0417 = IT_0194*(IT_0415 + -IT_0416);
    const ccomplex_t IT_0418 = -IT_0417;
    const ccomplex_t IT_0419 = 2*IT_0418;
    const ccomplex_t IT_0420 = IT_0122*(IT_0385 + -IT_0386);
    const ccomplex_t IT_0421 = 0.5*IT_0420;
    const ccomplex_t IT_0422 = (-2)*IT_0421;
    const ccomplex_t IT_0423 = IT_0414 + IT_0419 + IT_0422;
    const ccomplex_t IT_0424 = IT_0124 + IT_0215;
    const ccomplex_t IT_0425 = 2*conj(IT_0413);
    const ccomplex_t IT_0426 = 2*conj(IT_0418);
    const ccomplex_t IT_0427 = (-2)*conj(IT_0421);
    const ccomplex_t IT_0428 = IT_0425 + IT_0426 + IT_0427;
    const ccomplex_t IT_0429 = s_23*IT_0008*IT_0095;
    const ccomplex_t IT_0430 = IT_0365*IT_0429;
    const ccomplex_t IT_0431 = (-10)*s_23;
    const ccomplex_t IT_0432 = IT_0430 + IT_0431;
    const ccomplex_t IT_0433 = m_C_2*IT_0432;
    const ccomplex_t IT_0434 = -IT_0037;
    const ccomplex_t IT_0435 = IT_0059*IT_0434;
    const ccomplex_t IT_0436 = s_34*m_C_2*IT_0001;
    const ccomplex_t IT_0437 = s_23*s_34*IT_0002;
    const ccomplex_t IT_0438 = -IT_0437;
    const ccomplex_t IT_0439 = s_24 + IT_0438;
    const ccomplex_t IT_0440 = IT_0436*IT_0439;
    const ccomplex_t IT_0441 = s_34*m_N_4*IT_0001;
    const ccomplex_t IT_0442 = s_13*s_34*IT_0002;
    const ccomplex_t IT_0443 = -IT_0442;
    const ccomplex_t IT_0444 = s_14 + IT_0443;
    const ccomplex_t IT_0445 = IT_0441*IT_0444;
    const ccomplex_t IT_0446 = s_14*s_34*IT_0001;
    const ccomplex_t IT_0447 = s_12*IT_0095;
    const ccomplex_t IT_0448 = s_13*s_23;
    const ccomplex_t IT_0449 = (-2)*IT_0448;
    const ccomplex_t IT_0450 = IT_0447 + IT_0449;
    const ccomplex_t IT_0451 = IT_0002*IT_0450;
    const ccomplex_t IT_0452 = (s_13 + -IT_0446)*(s_12 + (-0.5)*IT_0451);
    const ccomplex_t IT_0453 = (-2)*IT_0452;
    const ccomplex_t IT_0454 = s_24*s_34*IT_0001;
    const ccomplex_t IT_0455 = (s_12 + (-0.5)*IT_0451)*(s_23 + -IT_0454);
    const ccomplex_t IT_0456 = (-2)*IT_0455;
    const ccomplex_t IT_0457 = IT_0453 + IT_0456;
    const ccomplex_t IT_0458 = 2*conj(IT_0342);
    const ccomplex_t IT_0459 = pow(s_23, 2);
    const ccomplex_t IT_0460 = IT_0001*IT_0002*IT_0323*IT_0459;
    const ccomplex_t IT_0461 = s_13*s_14*s_23*s_24*IT_0001*IT_0002;
    const ccomplex_t IT_0462 = IT_0009*IT_0010;
    const ccomplex_t IT_0463 = IT_0001*IT_0010*IT_0323;
    const ccomplex_t IT_0464 = -IT_0463;
    const ccomplex_t IT_0465 = IT_0002*IT_0009*IT_0459;
    const ccomplex_t IT_0466 = -IT_0465;
    const ccomplex_t IT_0467 = s_12*s_14*s_23*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0468 = -IT_0467;
    const ccomplex_t IT_0469 = IT_0460 + IT_0461 + IT_0462 + IT_0464 + IT_0466
       + IT_0468;
    const ccomplex_t IT_0470 = s_14*s_24*IT_0001*IT_0002*IT_0459;
    const ccomplex_t IT_0471 = pow(s_24, 2);
    const ccomplex_t IT_0472 = s_13*s_23*IT_0001*IT_0002*IT_0471;
    const ccomplex_t IT_0473 = s_12*IT_0010;
    const ccomplex_t IT_0474 = s_14*s_24*IT_0001*IT_0010;
    const ccomplex_t IT_0475 = -IT_0474;
    const ccomplex_t IT_0476 = s_13*s_23*IT_0002*IT_0010;
    const ccomplex_t IT_0477 = -IT_0476;
    const ccomplex_t IT_0478 = s_12*s_23*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0479 = -IT_0478;
    const ccomplex_t IT_0480 = IT_0470 + IT_0472 + IT_0473 + IT_0475 + IT_0477
       + IT_0479;
    const ccomplex_t IT_0481 = IT_0469 + IT_0480;
    const ccomplex_t IT_0482 = -IT_0123 + -IT_0343;
    const ccomplex_t IT_0483 = 2*IT_0342;
    const ccomplex_t IT_0484 = conj(IT_0482)*IT_0483;
    const ccomplex_t IT_0485 = 2*IT_0174 + IT_0176 + 2*IT_0195 + 2*IT_0214;
    const ccomplex_t IT_0486 = IT_0483*conj(IT_0485);
    const ccomplex_t IT_0487 = IT_0458*IT_0482;
    const ccomplex_t IT_0488 = IT_0458*IT_0485;
    const ccomplex_t IT_0489 = IT_0001*IT_0003;
    const ccomplex_t IT_0490 = (s_12 + (-0.5)*IT_0451)*(IT_0095 + -IT_0489);
    const ccomplex_t IT_0491 = (-2)*IT_0490;
    const ccomplex_t IT_0492 = 2*IT_0435;
    const ccomplex_t IT_0493 = IT_0048*IT_0434;
    const ccomplex_t IT_0494 = 2*IT_0337;
    const ccomplex_t IT_0495 = conj(IT_0493)*IT_0494;
    const ccomplex_t IT_0496 = 2*conj(IT_0337);
    const ccomplex_t IT_0497 = IT_0493*IT_0496;
    const ccomplex_t IT_0498 = s_13*s_23*IT_0001*IT_0002*IT_0323;
    const ccomplex_t IT_0499 = s_14*s_24*IT_0001*IT_0002*IT_0321;
    const ccomplex_t IT_0500 = s_12*IT_0009;
    const ccomplex_t IT_0501 = s_14*s_24*IT_0001*IT_0009;
    const ccomplex_t IT_0502 = -IT_0501;
    const ccomplex_t IT_0503 = s_13*s_23*IT_0002*IT_0009;
    const ccomplex_t IT_0504 = -IT_0503;
    const ccomplex_t IT_0505 = s_12*s_13*s_14*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0506 = -IT_0505;
    const ccomplex_t IT_0507 = IT_0498 + IT_0499 + IT_0500 + IT_0502 + IT_0504
       + IT_0506;
    const ccomplex_t IT_0508 = IT_0001*IT_0002*IT_0321*IT_0471;
    const ccomplex_t IT_0509 = IT_0002*IT_0010*IT_0321;
    const ccomplex_t IT_0510 = -IT_0509;
    const ccomplex_t IT_0511 = IT_0001*IT_0009*IT_0471;
    const ccomplex_t IT_0512 = -IT_0511;
    const ccomplex_t IT_0513 = s_12*s_13*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0514 = -IT_0513;
    const ccomplex_t IT_0515 = IT_0461 + IT_0462 + IT_0508 + IT_0510 + IT_0512
       + IT_0514;
    const ccomplex_t IT_0516 = conj(IT_0344)*IT_0483;
    const ccomplex_t IT_0517 = conj(IT_0345)*IT_0483;
    const ccomplex_t IT_0518 = IT_0344*IT_0458;
    const ccomplex_t IT_0519 = IT_0345*IT_0458;
    const ccomplex_t IT_0520 = conj(IT_0332)*IT_0494;
    const ccomplex_t IT_0521 = IT_0332*IT_0496;
    const ccomplex_t IT_0522 = s_12*s_13;
    const ccomplex_t IT_0523 = s_23*IT_0009;
    const ccomplex_t IT_0524 = s_13*s_14*s_23*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0525 = 2*IT_0524;
    const ccomplex_t IT_0526 = s_13*s_14*s_24*IT_0001;
    const ccomplex_t IT_0527 = -IT_0526;
    const ccomplex_t IT_0528 = s_23*IT_0001*IT_0323;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = s_23*IT_0002*IT_0321;
    const ccomplex_t IT_0531 = (-2)*IT_0530;
    const ccomplex_t IT_0532 = IT_0522 + IT_0523 + IT_0525 + IT_0527 + IT_0529
       + IT_0531;
    const ccomplex_t IT_0533 = s_12*s_23;
    const ccomplex_t IT_0534 = s_13*IT_0010;
    const ccomplex_t IT_0535 = s_13*s_23*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0536 = 2*IT_0535;
    const ccomplex_t IT_0537 = s_13*IT_0001*IT_0471;
    const ccomplex_t IT_0538 = -IT_0537;
    const ccomplex_t IT_0539 = s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0540 = -IT_0539;
    const ccomplex_t IT_0541 = s_13*IT_0002*IT_0459;
    const ccomplex_t IT_0542 = (-2)*IT_0541;
    const ccomplex_t IT_0543 = IT_0533 + IT_0534 + IT_0536 + IT_0538 + IT_0540
       + IT_0542;
    const ccomplex_t IT_0544 = 2*IT_0037;
    const ccomplex_t IT_0545 = IT_0059*IT_0544;
    const ccomplex_t IT_0546 = IT_0483*conj(IT_0545);
    const ccomplex_t IT_0547 = IT_0458*IT_0545;
    const ccomplex_t IT_0548 = IT_0048*IT_0544;
    const ccomplex_t IT_0549 = IT_0494*conj(IT_0548);
    const ccomplex_t IT_0550 = IT_0496*IT_0548;
    const ccomplex_t IT_0551 = s_14*s_34*IT_0001*IT_0095;
    const ccomplex_t IT_0552 = -IT_0551;
    const ccomplex_t IT_0553 = IT_0352 + IT_0552;
    const ccomplex_t IT_0554 = IT_0002*IT_0553;
    const ccomplex_t IT_0555 = -IT_0446 + (-0.25)*IT_0554;
    const ccomplex_t IT_0556 = s_13 + IT_0555;
    const ccomplex_t IT_0557 = IT_0000*IT_0556;
    const ccomplex_t IT_0558 = (-4)*IT_0557;
    const ccomplex_t IT_0559 = s_24*s_34*IT_0001*IT_0095;
    const ccomplex_t IT_0560 = -IT_0559;
    const ccomplex_t IT_0561 = IT_0096 + IT_0560;
    const ccomplex_t IT_0562 = IT_0002*IT_0561;
    const ccomplex_t IT_0563 = -IT_0454 + (-0.25)*IT_0562;
    const ccomplex_t IT_0564 = s_23 + IT_0563;
    const ccomplex_t IT_0565 = IT_0000*IT_0564;
    const ccomplex_t IT_0566 = (-4)*IT_0565;
    const ccomplex_t IT_0567 = conj(IT_0269)*IT_0483;
    const ccomplex_t IT_0568 = conj(IT_0291)*IT_0483;
    const ccomplex_t IT_0569 = conj(IT_0310)*IT_0483;
    const ccomplex_t IT_0570 = conj(IT_0240)*IT_0483;
    const ccomplex_t IT_0571 = conj(IT_0250)*IT_0483;
    const ccomplex_t IT_0572 = IT_0269*IT_0458;
    const ccomplex_t IT_0573 = IT_0291*IT_0458;
    const ccomplex_t IT_0574 = IT_0310*IT_0458;
    const ccomplex_t IT_0575 = IT_0240*IT_0458;
    const ccomplex_t IT_0576 = IT_0250*IT_0458;
    const ccomplex_t IT_0577 = conj(IT_0124)*IT_0494;
    const ccomplex_t IT_0578 = conj(IT_0143)*IT_0494;
    const ccomplex_t IT_0579 = conj(IT_0155)*IT_0494;
    const ccomplex_t IT_0580 = conj(IT_0215)*IT_0494;
    const ccomplex_t IT_0581 = IT_0124*IT_0496;
    const ccomplex_t IT_0582 = IT_0143*IT_0496;
    const ccomplex_t IT_0583 = IT_0155*IT_0496;
    const ccomplex_t IT_0584 = IT_0215*IT_0496;
    const ccomplex_t IT_0585 = IT_0001*IT_0471;
    const ccomplex_t IT_0586 = -IT_0010;
    const ccomplex_t IT_0587 = IT_0585 + IT_0586;
    const ccomplex_t IT_0588 = IT_0000*IT_0587;
    const ccomplex_t IT_0589 = -IT_0009;
    const ccomplex_t IT_0590 = IT_0324 + IT_0589;
    const ccomplex_t IT_0591 = IT_0000*IT_0590;
    const ccomplex_t IT_0592 = IT_0588 + IT_0591;
    const ccomplex_t IT_0593 = IT_0337*conj(IT_0342);
    const ccomplex_t IT_0594 = 6*IT_0593;
    const ccomplex_t IT_0595 = -IT_0446;
    const ccomplex_t IT_0596 = s_13 + IT_0595;
    const ccomplex_t IT_0597 = IT_0000*IT_0596;
    const ccomplex_t IT_0598 = -IT_0454;
    const ccomplex_t IT_0599 = s_23 + IT_0598;
    const ccomplex_t IT_0600 = IT_0000*IT_0599;
    const ccomplex_t IT_0601 = IT_0597 + IT_0600;
    const ccomplex_t IT_0602 = conj(IT_0342)*IT_0493;
    const ccomplex_t IT_0603 = (-6)*IT_0602;
    const ccomplex_t IT_0604 = (-6)*IT_0435;
    const ccomplex_t IT_0605 = conj(IT_0337)*IT_0604;
    const ccomplex_t IT_0606 = (-6)*conj(IT_0435);
    const ccomplex_t IT_0607 = IT_0337*IT_0606;
    const ccomplex_t IT_0608 = s_12*IT_0006;
    const ccomplex_t IT_0609 = (-4)*conj(IT_0088);
    const ccomplex_t IT_0610 = conj(IT_0050) + 0.25*IT_0609;
    const ccomplex_t IT_0611 = (-4)*IT_0088;
    const ccomplex_t IT_0612 = IT_0002*IT_0459;
    const ccomplex_t IT_0613 = IT_0585 + IT_0612;
    const ccomplex_t IT_0614 = s_23*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0615 = -IT_0010 + -IT_0614;
    const ccomplex_t IT_0616 = IT_0613 + IT_0615;
    const ccomplex_t IT_0617 = m_C_2*IT_0616;
    const ccomplex_t IT_0618 = s_14*s_23*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0619 = s_13*s_23*IT_0002;
    const ccomplex_t IT_0620 = s_14*s_24*IT_0001;
    const ccomplex_t IT_0621 = -IT_0619 + -IT_0620;
    const ccomplex_t IT_0622 = s_12 + IT_0618 + IT_0621;
    const ccomplex_t IT_0623 = m_C_2*IT_0622;
    const ccomplex_t IT_0624 = s_13*s_24*s_34*IT_0001*IT_0002;
    const ccomplex_t IT_0625 = s_12 + IT_0621 + IT_0624;
    const ccomplex_t IT_0626 = m_N_4*IT_0625;
    const ccomplex_t IT_0627 = (s_12 + (-0.5)*IT_0451)*(s_12 + -IT_0620);
    const ccomplex_t IT_0628 = (-2)*IT_0627;
    const ccomplex_t IT_0629 = 4*IT_0342;
    const ccomplex_t IT_0630 = -(IT_0009 + -IT_0324)*(s_12 + (-0.5)*IT_0451);
    const ccomplex_t IT_0631 = (-2)*IT_0630;
    const ccomplex_t IT_0632 = -(s_12 + (-0.5)*IT_0451)*(IT_0010 + -IT_0585);
    const ccomplex_t IT_0633 = (-2)*IT_0632;
    const ccomplex_t IT_0634 = IT_0631 + IT_0633;
    const ccomplex_t IT_0635 = (-2)*IT_0342;
    const ccomplex_t IT_0636 = m_N_4*IT_0622;
    const ccomplex_t IT_0637 = s_14*s_34*IT_0001*IT_0002*IT_0459;
    const ccomplex_t IT_0638 = IT_0533 + IT_0534 + IT_0535;
    const ccomplex_t IT_0639 = s_12*s_23*IT_0001*IT_0002*IT_0003;
    const ccomplex_t IT_0640 = -IT_0639;
    const ccomplex_t IT_0641 = s_14*s_34*IT_0001*IT_0010;
    const ccomplex_t IT_0642 = -IT_0641;
    const ccomplex_t IT_0643 = IT_0542 + IT_0637 + IT_0638 + IT_0640 + IT_0642;
    const ccomplex_t IT_0644 = s_12*IT_0008;
    const ccomplex_t IT_0645 = s_14*s_24;
    const ccomplex_t IT_0646 = (-2)*IT_0645;
    const ccomplex_t IT_0647 = IT_0644 + IT_0646;
    const ccomplex_t IT_0648 = IT_0001*IT_0647;
    const ccomplex_t IT_0649 = (s_12 + -IT_0619)*(s_12 + (-0.5)*IT_0648);
    const ccomplex_t IT_0650 = (-2)*IT_0649;
    const ccomplex_t IT_0651 = 2*IT_0344;
    const ccomplex_t IT_0652 = -(IT_0010 + -IT_0612)*(s_12 + (-0.5)*IT_0648);
    const ccomplex_t IT_0653 = (-2)*IT_0652;
    const ccomplex_t IT_0654 = (-2)*IT_0482;
    const ccomplex_t IT_0655 = m_C_2*IT_0625;
    const ccomplex_t IT_0656 = s_24*s_34*IT_0001*IT_0002*IT_0321;
    const ccomplex_t IT_0657 = IT_0522 + IT_0523 + IT_0524;
    const ccomplex_t IT_0658 = s_12*s_13*IT_0001*IT_0002*IT_0003;
    const ccomplex_t IT_0659 = -IT_0658;
    const ccomplex_t IT_0660 = s_24*s_34*IT_0001*IT_0009;
    const ccomplex_t IT_0661 = -IT_0660;
    const ccomplex_t IT_0662 = IT_0531 + IT_0656 + IT_0657 + IT_0659 + IT_0661;
    const ccomplex_t IT_0663 = -(IT_0009 + -IT_0322)*(s_12 + (-0.5)*IT_0648);
    const ccomplex_t IT_0664 = (-2)*IT_0663;
    const ccomplex_t IT_0665 = (-2)*IT_0344;
    const ccomplex_t IT_0666 = IT_0057*IT_0247;
    const ccomplex_t IT_0667 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0668 = IT_0152*IT_0666*IT_0667;
    const ccomplex_t IT_0669 = m_C_2*IT_0154;
    const ccomplex_t IT_0670 = IT_0131*IT_0237;
    const ccomplex_t IT_0671 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0672 = IT_0140*IT_0670*IT_0671;
    const ccomplex_t IT_0673 = -IT_0668 + -IT_0669 + -IT_0672;
    const ccomplex_t IT_0674 = s_13*s_23*IT_0002*IT_0008;
    const ccomplex_t IT_0675 = s_13*s_24*s_34*IT_0002;
    const ccomplex_t IT_0676 = (-2)*IT_0645 + -IT_0674 + 2*IT_0675;
    const ccomplex_t IT_0677 = IT_0644 + IT_0676;
    const ccomplex_t IT_0678 = IT_0001*IT_0677;
    const ccomplex_t IT_0679 = -IT_0619 + (-0.5)*IT_0678;
    const ccomplex_t IT_0680 = s_12 + IT_0679;
    const ccomplex_t IT_0681 = m_C_2*IT_0680;
    const ccomplex_t IT_0682 = (-2)*IT_0681;
    const ccomplex_t IT_0683 = 2*IT_0345;
    const ccomplex_t IT_0684 = IT_0651 + IT_0683;
    const ccomplex_t IT_0685 = IT_0002*IT_0447;
    const ccomplex_t IT_0686 = (-0.25)*IT_0685;
    const ccomplex_t IT_0687 = IT_0001*IT_0644;
    const ccomplex_t IT_0688 = (-0.25)*IT_0687;
    const ccomplex_t IT_0689 = s_12 + IT_0686 + IT_0688;
    const ccomplex_t IT_0690 = s_14*s_23;
    const ccomplex_t IT_0691 = s_13*s_24;
    const ccomplex_t IT_0692 = s_12*s_34;
    const ccomplex_t IT_0693 = -IT_0691 + -IT_0692;
    const ccomplex_t IT_0694 = IT_0690 + IT_0693;
    const ccomplex_t IT_0695 = IT_0098*IT_0694;
    const ccomplex_t IT_0696 = (-0.25)*IT_0695;
    const ccomplex_t IT_0697 = IT_0689 + IT_0696;
    const ccomplex_t IT_0698 = 8*IT_0082;
    const ccomplex_t IT_0699 = (-8)*IT_0060;
    const ccomplex_t IT_0700 = IT_0698 + IT_0699;
    const ccomplex_t IT_0701 = s_12*IT_0008*IT_0095;
    const ccomplex_t IT_0702 = s_13*s_23*IT_0008;
    const ccomplex_t IT_0703 = s_13*s_24*s_34;
    const ccomplex_t IT_0704 = s_14*s_24*IT_0095;
    const ccomplex_t IT_0705 = (-2)*IT_0702 + 4*IT_0703 + (-2)*IT_0704;
    const ccomplex_t IT_0706 = IT_0701 + IT_0705;
    const ccomplex_t IT_0707 = IT_0365*IT_0706;
    const ccomplex_t IT_0708 = 0.25*IT_0707;
    const ccomplex_t IT_0709 = (-0.5)*IT_0648;
    const ccomplex_t IT_0710 = (-0.5)*IT_0451;
    const ccomplex_t IT_0711 = s_12 + IT_0708 + IT_0709 + IT_0710;
    const ccomplex_t IT_0712 = 8*IT_0673;
    const ccomplex_t IT_0713 = m_C_2*IT_0142;
    const ccomplex_t IT_0714 = (-8)*IT_0713;
    const ccomplex_t IT_0715 = IT_0712 + IT_0714;
    const ccomplex_t IT_0716 = IT_0002*IT_0008*IT_0459;
    const ccomplex_t IT_0717 = s_23*s_24*s_34*IT_0002;
    const ccomplex_t IT_0718 = IT_0008*IT_0010;
    const ccomplex_t IT_0719 = 2*IT_0471 + (-2)*IT_0717 + -IT_0718;
    const ccomplex_t IT_0720 = IT_0716 + IT_0719;
    const ccomplex_t IT_0721 = IT_0001*IT_0720;
    const ccomplex_t IT_0722 = -IT_0010 + (-0.5)*IT_0721;
    const ccomplex_t IT_0723 = IT_0612 + IT_0722;
    const ccomplex_t IT_0724 = m_C_2*IT_0723;
    const ccomplex_t IT_0725 = (-2)*IT_0724;
    const ccomplex_t IT_0726 = (-2)*IT_0485;
    const ccomplex_t IT_0727 = IT_0654 + IT_0726;
    const ccomplex_t IT_0728 = IT_0001*IT_0095*IT_0471;
    const ccomplex_t IT_0729 = IT_0010*IT_0095;
    const ccomplex_t IT_0730 = -IT_0729;
    const ccomplex_t IT_0731 = IT_0728 + IT_0730;
    const ccomplex_t IT_0732 = IT_0002*IT_0731;
    const ccomplex_t IT_0733 = (-4)*IT_0587;
    const ccomplex_t IT_0734 = IT_0732 + IT_0733;
    const ccomplex_t IT_0735 = m_C_2*IT_0734;
    const ccomplex_t IT_0736 = -IT_0674;
    const ccomplex_t IT_0737 = IT_0644 + IT_0736;
    const ccomplex_t IT_0738 = IT_0001*IT_0737;
    const ccomplex_t IT_0739 = -IT_0619 + (-0.25)*IT_0738;
    const ccomplex_t IT_0740 = s_12 + IT_0739;
    const ccomplex_t IT_0741 = m_N_4*IT_0740;
    const ccomplex_t IT_0742 = (-4)*IT_0741;
    const ccomplex_t IT_0743 = 2*IT_0268 + 2*IT_0287 + (-2)*IT_0290 + 2
      *IT_0309 + -IT_0331;
    const ccomplex_t IT_0744 = 2*IT_0743;
    const ccomplex_t IT_0745 = IT_0352 + IT_0551;
    const ccomplex_t IT_0746 = s_13*IT_0001*IT_0003;
    const ccomplex_t IT_0747 = (-2)*IT_0746;
    const ccomplex_t IT_0748 = IT_0745 + IT_0747;
    const ccomplex_t IT_0749 = IT_0002*IT_0748;
    const ccomplex_t IT_0750 = -IT_0446 + 0.5*IT_0749;
    const ccomplex_t IT_0751 = s_13 + IT_0750;
    const ccomplex_t IT_0752 = m_N_4*IT_0751;
    const ccomplex_t IT_0753 = 2*IT_0752;
    const ccomplex_t IT_0754 = 2*IT_0493;
    const ccomplex_t IT_0755 = IT_0001*IT_0356;
    const ccomplex_t IT_0756 = s_13*s_34;
    const ccomplex_t IT_0757 = (-2)*IT_0756;
    const ccomplex_t IT_0758 = IT_0354 + IT_0757;
    const ccomplex_t IT_0759 = IT_0098*IT_0758;
    const ccomplex_t IT_0760 = (-0.25)*IT_0353 + (-0.25)*IT_0755 + (-0.25)
      *IT_0759;
    const ccomplex_t IT_0761 = s_13 + IT_0760;
    const ccomplex_t IT_0762 = m_N_4*IT_0761;
    const ccomplex_t IT_0763 = (-4)*IT_0762;
    const ccomplex_t IT_0764 = 2*IT_0548;
    const ccomplex_t IT_0765 = IT_0001*IT_0095*IT_0323;
    const ccomplex_t IT_0766 = IT_0009*IT_0095;
    const ccomplex_t IT_0767 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0768 = 0.5*IT_0765 + (-0.5)*IT_0766 + -IT_0767;
    const ccomplex_t IT_0769 = IT_0321 + IT_0768;
    const ccomplex_t IT_0770 = IT_0002*IT_0769;
    const ccomplex_t IT_0771 = -IT_0324;
    const ccomplex_t IT_0772 = IT_0009 + IT_0770 + IT_0771;
    const ccomplex_t IT_0773 = m_N_4*IT_0772;
    const ccomplex_t IT_0774 = 4*IT_0337;
    const ccomplex_t IT_0775 = s_14*s_23*s_34;
    const ccomplex_t IT_0776 = s_12*IT_0003;
    const ccomplex_t IT_0777 = (-2)*IT_0702 + 2*IT_0703 + (-2)*IT_0704 + 2
      *IT_0775 + (-2)*IT_0776;
    const ccomplex_t IT_0778 = IT_0701 + IT_0777;
    const ccomplex_t IT_0779 = IT_0365*IT_0778;
    const ccomplex_t IT_0780 = (-0.25)*IT_0779;
    const ccomplex_t IT_0781 = -IT_0620;
    const ccomplex_t IT_0782 = -IT_0619;
    const ccomplex_t IT_0783 = s_12 + IT_0780 + IT_0781 + IT_0782;
    const ccomplex_t IT_0784 = 8*IT_0390;
    const ccomplex_t IT_0785 = IT_0008*IT_0095;
    const ccomplex_t IT_0786 = (-2)*IT_0003;
    const ccomplex_t IT_0787 = IT_0785 + IT_0786;
    const ccomplex_t IT_0788 = IT_0365*IT_0787;
    const ccomplex_t IT_0789 = 0.25*IT_0788;
    const ccomplex_t IT_0790 = 1 + IT_0789;
    const ccomplex_t IT_0791 = IT_0000*IT_0790;
    const ccomplex_t IT_0792 = (-8)*IT_0421;
    const ccomplex_t IT_0793 = s_14*s_24*IT_0001*IT_0095;
    const ccomplex_t IT_0794 = -IT_0793;
    const ccomplex_t IT_0795 = IT_0447 + IT_0794;
    const ccomplex_t IT_0796 = IT_0002*IT_0795;
    const ccomplex_t IT_0797 = -IT_0620 + (-0.25)*IT_0796;
    const ccomplex_t IT_0798 = s_12 + IT_0797;
    const ccomplex_t IT_0799 = m_C_2*IT_0798;
    const ccomplex_t IT_0800 = (-4)*IT_0799;
    const ccomplex_t IT_0801 = m_C_2*IT_0564;
    const ccomplex_t IT_0802 = (-4)*IT_0801;
    const ccomplex_t IT_0803 = (-2)*IT_0435;
    const ccomplex_t IT_0804 = IT_0001*IT_0101;
    const ccomplex_t IT_0805 = (-0.25)*IT_0097 + 0.25*IT_0100 + (-0.25)*IT_0804;
    const ccomplex_t IT_0806 = s_23 + IT_0805;
    const ccomplex_t IT_0807 = m_C_2*IT_0806;
    const ccomplex_t IT_0808 = (-4)*IT_0807;
    const ccomplex_t IT_0809 = (-2)*IT_0545;
    const ccomplex_t IT_0810 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0811 = (-2)*IT_0448 + -IT_0793 + 2*IT_0810;
    const ccomplex_t IT_0812 = IT_0447 + IT_0811;
    const ccomplex_t IT_0813 = IT_0002*IT_0812;
    const ccomplex_t IT_0814 = -IT_0620 + (-0.5)*IT_0813;
    const ccomplex_t IT_0815 = s_12 + IT_0814;
    const ccomplex_t IT_0816 = m_N_4*IT_0815;
    const ccomplex_t IT_0817 = (-2)*IT_0816;
    const ccomplex_t IT_0818 = (-2)*IT_0337;
    const ccomplex_t IT_0819 = IT_0002*IT_0008*IT_0321;
    const ccomplex_t IT_0820 = IT_0008*IT_0009;
    const ccomplex_t IT_0821 = -IT_0820;
    const ccomplex_t IT_0822 = IT_0819 + IT_0821;
    const ccomplex_t IT_0823 = IT_0001*IT_0822;
    const ccomplex_t IT_0824 = IT_0322 + IT_0589;
    const ccomplex_t IT_0825 = (-4)*IT_0824;
    const ccomplex_t IT_0826 = IT_0823 + IT_0825;
    const ccomplex_t IT_0827 = m_N_4*IT_0826;
    const ccomplex_t IT_0828 = (-2)*IT_0332;
    const ccomplex_t IT_0829 = IT_0682*IT_0684 + IT_0697*IT_0700 + IT_0711
      *IT_0715 + IT_0725*IT_0727 + IT_0483*IT_0735 + IT_0742*IT_0744 + IT_0753
      *IT_0754 + IT_0763*IT_0764 + IT_0773*IT_0774 + IT_0783*IT_0784 + IT_0791
      *IT_0792 + IT_0635*IT_0800 + IT_0802*IT_0803 + IT_0808*IT_0809 + IT_0817
      *IT_0818 + IT_0827*IT_0828;
    const ccomplex_t IT_0830 = 2*conj(IT_0344);
    const ccomplex_t IT_0831 = 2*conj(IT_0345);
    const ccomplex_t IT_0832 = IT_0830 + IT_0831;
    const ccomplex_t IT_0833 = 8*conj(IT_0082);
    const ccomplex_t IT_0834 = (-8)*conj(IT_0060);
    const ccomplex_t IT_0835 = IT_0833 + IT_0834;
    const ccomplex_t IT_0836 = (-2)*conj(IT_0482);
    const ccomplex_t IT_0837 = (-2)*conj(IT_0485);
    const ccomplex_t IT_0838 = IT_0836 + IT_0837;
    const ccomplex_t IT_0839 = 2*conj(IT_0743);
    const ccomplex_t IT_0840 = 2*conj(IT_0493);
    const ccomplex_t IT_0841 = 2*conj(IT_0548);
    const ccomplex_t IT_0842 = 4*conj(IT_0337);
    const ccomplex_t IT_0843 = 8*conj(IT_0390);
    const ccomplex_t IT_0844 = (-8)*conj(IT_0713);
    const ccomplex_t IT_0845 = (-8)*conj(IT_0421);
    const ccomplex_t IT_0846 = (-2)*conj(IT_0342);
    const ccomplex_t IT_0847 = (-2)*conj(IT_0435);
    const ccomplex_t IT_0848 = (-2)*conj(IT_0545);
    const ccomplex_t IT_0849 = (-2)*conj(IT_0337);
    const ccomplex_t IT_0850 = (-2)*conj(IT_0332);
    const ccomplex_t IT_0851 = IT_0458*IT_0735 + IT_0682*IT_0832 + IT_0697
      *IT_0835 + IT_0725*IT_0838 + IT_0742*IT_0839 + IT_0753*IT_0840 + IT_0763
      *IT_0841 + IT_0773*IT_0842 + IT_0783*IT_0843 + IT_0711*IT_0844 + IT_0791
      *IT_0845 + IT_0800*IT_0846 + IT_0802*IT_0847 + IT_0808*IT_0848 + IT_0817
      *IT_0849 + IT_0827*IT_0850;
    const ccomplex_t IT_0852 = IT_0138*IT_0230;
    const ccomplex_t IT_0853 = IT_0140*IT_0671*IT_0852;
    const ccomplex_t IT_0854 = IT_0045*IT_0150;
    const ccomplex_t IT_0855 = IT_0152*IT_0667*IT_0854;
    const ccomplex_t IT_0856 = m_C_2*IT_0240;
    const ccomplex_t IT_0857 = m_C_2*IT_0250;
    const ccomplex_t IT_0858 = IT_0853 + IT_0855 + IT_0856 + IT_0857;
    const ccomplex_t IT_0859 = IT_0365*IT_0785;
    const ccomplex_t IT_0860 = 0.125*IT_0859;
    const ccomplex_t IT_0861 = 1 + IT_0860;
    const ccomplex_t IT_0862 = IT_0000*IT_0861;
    const ccomplex_t IT_0863 = 0.125*IT_0792;
    const ccomplex_t IT_0864 = 8*IT_0413;
    const ccomplex_t IT_0865 = 0.125*IT_0864;
    const ccomplex_t IT_0866 = IT_0418 + IT_0863 + IT_0865;
    const ccomplex_t IT_0867 = 2*IT_0482;
    const ccomplex_t IT_0868 = 2*IT_0485;
    const ccomplex_t IT_0869 = IT_0867 + IT_0868;
    const ccomplex_t IT_0870 = 8*IT_0050;
    const ccomplex_t IT_0871 = (-8)*IT_0088;
    const ccomplex_t IT_0872 = IT_0870 + IT_0871;
    const ccomplex_t IT_0873 = 8*IT_0395;
    const ccomplex_t IT_0874 = IT_0784 + IT_0873;
    const ccomplex_t IT_0875 = (-2)*IT_0345;
    const ccomplex_t IT_0876 = IT_0665 + IT_0875;
    const ccomplex_t IT_0877 = 2*IT_0545;
    const ccomplex_t IT_0878 = 2*IT_0332;
    const ccomplex_t IT_0879 = 8*IT_0858;
    const ccomplex_t IT_0880 = (-2)*IT_0743;
    const ccomplex_t IT_0881 = (-2)*IT_0493;
    const ccomplex_t IT_0882 = (-2)*IT_0548;
    const ccomplex_t IT_0883 = IT_0494*IT_0735 + IT_0492*IT_0753 + IT_0629
      *IT_0773 + IT_0635*IT_0817 + IT_0800*IT_0818 + ((-16)*IT_0673 + 16*IT_0713
      )*IT_0862 + 8*IT_0783*IT_0866 + IT_0742*IT_0869 + IT_0697*IT_0872 +
       IT_0791*IT_0874 + IT_0827*IT_0876 + IT_0763*IT_0877 + IT_0682*IT_0878 +
       IT_0711*IT_0879 + IT_0725*IT_0880 + IT_0802*IT_0881 + IT_0808*IT_0882;
    const ccomplex_t IT_0884 = s_13*s_14*s_34*IT_0002;
    const ccomplex_t IT_0885 = 2*IT_0323 + -IT_0820 + (-2)*IT_0884;
    const ccomplex_t IT_0886 = IT_0819 + IT_0885;
    const ccomplex_t IT_0887 = IT_0001*IT_0886;
    const ccomplex_t IT_0888 = -IT_0009 + (-0.5)*IT_0887;
    const ccomplex_t IT_0889 = IT_0322 + IT_0888;
    const ccomplex_t IT_0890 = m_N_4*IT_0889;
    const ccomplex_t IT_0891 = (-2)*IT_0890;
    const ccomplex_t IT_0892 = IT_0690 + IT_0692;
    const ccomplex_t IT_0893 = -IT_0691;
    const ccomplex_t IT_0894 = IT_0892 + IT_0893;
    const ccomplex_t IT_0895 = IT_0098*IT_0894;
    const ccomplex_t IT_0896 = 0.25*IT_0895;
    const ccomplex_t IT_0897 = IT_0689 + IT_0896;
    const ccomplex_t IT_0898 = (-16)*IT_0390;
    const ccomplex_t IT_0899 = (-16)*IT_0395;
    const ccomplex_t IT_0900 = IT_0898 + IT_0899;
    const ccomplex_t IT_0901 = s_14*s_23*s_34*IT_0002;
    const ccomplex_t IT_0902 = (-2)*IT_0645 + -IT_0674 + 2*IT_0901;
    const ccomplex_t IT_0903 = IT_0644 + IT_0902;
    const ccomplex_t IT_0904 = IT_0001*IT_0903;
    const ccomplex_t IT_0905 = -IT_0619 + (-0.5)*IT_0904;
    const ccomplex_t IT_0906 = s_12 + IT_0905;
    const ccomplex_t IT_0907 = m_N_4*IT_0906;
    const ccomplex_t IT_0908 = (-2)*IT_0907;
    const ccomplex_t IT_0909 = m_N_4*IT_0798;
    const ccomplex_t IT_0910 = (-4)*IT_0909;
    const ccomplex_t IT_0911 = m_N_4*IT_0556;
    const ccomplex_t IT_0912 = (-4)*IT_0911;
    const ccomplex_t IT_0913 = (-0.25)*IT_0353 + 0.25*IT_0355 + (-0.25)*IT_0755;
    const ccomplex_t IT_0914 = s_13 + IT_0913;
    const ccomplex_t IT_0915 = m_N_4*IT_0914;
    const ccomplex_t IT_0916 = (-4)*IT_0915;
    const ccomplex_t IT_0917 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0918 = (-2)*IT_0448 + -IT_0793 + 2*IT_0917;
    const ccomplex_t IT_0919 = IT_0447 + IT_0918;
    const ccomplex_t IT_0920 = IT_0002*IT_0919;
    const ccomplex_t IT_0921 = -IT_0620 + (-0.5)*IT_0920;
    const ccomplex_t IT_0922 = s_12 + IT_0921;
    const ccomplex_t IT_0923 = m_C_2*IT_0922;
    const ccomplex_t IT_0924 = (-2)*IT_0923;
    const ccomplex_t IT_0925 = -IT_0718;
    const ccomplex_t IT_0926 = IT_0716 + IT_0925;
    const ccomplex_t IT_0927 = IT_0001*IT_0926;
    const ccomplex_t IT_0928 = IT_0586 + IT_0612;
    const ccomplex_t IT_0929 = (-4)*IT_0928;
    const ccomplex_t IT_0930 = IT_0927 + IT_0929;
    const ccomplex_t IT_0931 = m_C_2*IT_0930;
    const ccomplex_t IT_0932 = (-2)*IT_0702 + (-2)*IT_0704 + 4*IT_0775;
    const ccomplex_t IT_0933 = IT_0701 + IT_0932;
    const ccomplex_t IT_0934 = IT_0365*IT_0933;
    const ccomplex_t IT_0935 = 0.25*IT_0934;
    const ccomplex_t IT_0936 = s_12 + IT_0709 + IT_0710 + IT_0935;
    const ccomplex_t IT_0937 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0938 = 0.5*IT_0728 + (-0.5)*IT_0729 + -IT_0937;
    const ccomplex_t IT_0939 = IT_0459 + IT_0938;
    const ccomplex_t IT_0940 = IT_0002*IT_0939;
    const ccomplex_t IT_0941 = -IT_0585;
    const ccomplex_t IT_0942 = IT_0010 + IT_0940 + IT_0941;
    const ccomplex_t IT_0943 = m_C_2*IT_0942;
    const ccomplex_t IT_0944 = (-4)*IT_0337;
    const ccomplex_t IT_0945 = -IT_0766;
    const ccomplex_t IT_0946 = IT_0765 + IT_0945;
    const ccomplex_t IT_0947 = IT_0002*IT_0946;
    const ccomplex_t IT_0948 = (-4)*IT_0590;
    const ccomplex_t IT_0949 = IT_0947 + IT_0948;
    const ccomplex_t IT_0950 = m_N_4*IT_0949;
    const ccomplex_t IT_0951 = IT_0096 + IT_0559;
    const ccomplex_t IT_0952 = s_23*IT_0001*IT_0003;
    const ccomplex_t IT_0953 = (-2)*IT_0952;
    const ccomplex_t IT_0954 = IT_0951 + IT_0953;
    const ccomplex_t IT_0955 = IT_0002*IT_0954;
    const ccomplex_t IT_0956 = -IT_0454 + 0.5*IT_0955;
    const ccomplex_t IT_0957 = s_23 + IT_0956;
    const ccomplex_t IT_0958 = m_C_2*IT_0957;
    const ccomplex_t IT_0959 = 2*IT_0958;
    const ccomplex_t IT_0960 = m_C_2*IT_0740;
    const ccomplex_t IT_0961 = (-4)*IT_0960;
    const ccomplex_t IT_0962 = s_23*s_34;
    const ccomplex_t IT_0963 = (-2)*IT_0962;
    const ccomplex_t IT_0964 = IT_0099 + IT_0963;
    const ccomplex_t IT_0965 = IT_0098*IT_0964;
    const ccomplex_t IT_0966 = (-0.25)*IT_0097 + (-0.25)*IT_0804 + (-0.25)
      *IT_0965;
    const ccomplex_t IT_0967 = s_23 + IT_0966;
    const ccomplex_t IT_0968 = m_C_2*IT_0967;
    const ccomplex_t IT_0969 = (-4)*IT_0968;
    const ccomplex_t IT_0970 = IT_0715*IT_0791 + IT_0684*IT_0891 + IT_0872
      *IT_0897 + IT_0862*IT_0900 + IT_0727*IT_0908 + IT_0483*IT_0910 + IT_0492
      *IT_0912 + IT_0877*IT_0916 + IT_0494*IT_0924 + IT_0744*IT_0931 + (IT_0792 
      + IT_0864)*IT_0936 + IT_0943*IT_0944 + IT_0635*IT_0950 + IT_0881*IT_0959 +
       IT_0828*IT_0961 + IT_0882*IT_0969;
    const ccomplex_t IT_0971 = s_34*m_C_2*m_N_4*IT_0001;
    const ccomplex_t IT_0972 = IT_0439*IT_0971;
    const ccomplex_t IT_0973 = 2*conj(IT_0482);
    const ccomplex_t IT_0974 = IT_0444*IT_0971;
    const ccomplex_t IT_0975 = 2*conj(IT_0421);
    const ccomplex_t IT_0976 = -IT_0489;
    const ccomplex_t IT_0977 = IT_0095 + IT_0976;
    const ccomplex_t IT_0978 = IT_0000*IT_0977;
    const ccomplex_t IT_0979 = (-2)*conj(IT_0713);
    const ccomplex_t IT_0980 = s_12 + IT_0781;
    const ccomplex_t IT_0981 = IT_0000*IT_0980;
    const ccomplex_t IT_0982 = IT_0000*IT_0616;
    const ccomplex_t IT_0983 = IT_0000*IT_0328;
    const ccomplex_t IT_0984 = 2*conj(IT_0713);
    const ccomplex_t IT_0985 = IT_0000*IT_0622;
    const ccomplex_t IT_0986 = IT_0000*IT_0625;
    const ccomplex_t IT_0987 = (-2)*conj(IT_0344);
    const ccomplex_t IT_0988 = IT_0320*IT_0617 + IT_0351*IT_0623 + IT_0223
      *IT_0626 + (-4)*conj(IT_0713)*IT_0773 + IT_0427*IT_0924 + 4*conj(IT_0421)
      *IT_0943 + (-12)*conj(IT_0342)*IT_0981 + IT_0973*IT_0982 + IT_0830*IT_0983
       + IT_0817*IT_0984 + IT_0836*IT_0985 + IT_0986*IT_0987;
    const ccomplex_t IT_0989 = m_N_4*IT_0001;
    const ccomplex_t IT_0990 = IT_0002*IT_0366;
    const ccomplex_t IT_0991 = (-2)*IT_0357 + (-0.5)*IT_0990;
    const ccomplex_t IT_0992 = IT_0356 + IT_0991;
    const ccomplex_t IT_0993 = IT_0989*IT_0992;
    const ccomplex_t IT_0994 = (-2)*IT_0993;
    const ccomplex_t IT_0995 = 2*IT_0713;
    const ccomplex_t IT_0996 = (-2)*IT_0673;
    const ccomplex_t IT_0997 = IT_0995 + IT_0996;
    const ccomplex_t IT_0998 = IT_0001*IT_0003*IT_0095;
    const ccomplex_t IT_0999 = -IT_0998;
    const ccomplex_t IT_1000 = pow(m_Z, 4);
    const ccomplex_t IT_1001 = IT_0999 + IT_1000;
    const ccomplex_t IT_1002 = IT_0002*IT_1001;
    const ccomplex_t IT_1003 = -IT_0095 + 0.25*IT_1002;
    const ccomplex_t IT_1004 = IT_0489 + IT_1003;
    const ccomplex_t IT_1005 = IT_0000*IT_1004;
    const ccomplex_t IT_1006 = 4*IT_1005;
    const ccomplex_t IT_1007 = s_34*IT_0095;
    const ccomplex_t IT_1008 = IT_0098*IT_1007;
    const ccomplex_t IT_1009 = IT_0001*IT_0787;
    const ccomplex_t IT_1010 = (-0.333333333333333)*IT_1008 + (
      -0.333333333333333)*IT_1009;
    const ccomplex_t IT_1011 = IT_0095 + IT_1010;
    const ccomplex_t IT_1012 = IT_0000*IT_1011;
    const ccomplex_t IT_1013 = 3*IT_1012;
    const ccomplex_t IT_1014 = s_12*s_14*s_34*IT_0001;
    const ccomplex_t IT_1015 = -IT_0523;
    const ccomplex_t IT_1016 = s_23*IT_0001*IT_0095*IT_0323;
    const ccomplex_t IT_1017 = s_13*s_14*s_24*IT_0001*IT_0095;
    const ccomplex_t IT_1018 = IT_1016 + IT_1017;
    const ccomplex_t IT_1019 = s_12*s_14*s_34*IT_0001*IT_0095;
    const ccomplex_t IT_1020 = s_23*IT_0009*IT_0095;
    const ccomplex_t IT_1021 = -IT_1019 + -IT_1020;
    const ccomplex_t IT_1022 = IT_1018 + IT_1021;
    const ccomplex_t IT_1023 = IT_0002*IT_1022;
    const ccomplex_t IT_1024 = 0.5*IT_1023;
    const ccomplex_t IT_1025 = IT_0527 + IT_0528 + IT_1014 + IT_1015 + IT_1024;
    const ccomplex_t IT_1026 = s_12*s_24*s_34*IT_0001;
    const ccomplex_t IT_1027 = (-0.5)*IT_1026;
    const ccomplex_t IT_1028 = s_14*s_23*s_24*IT_0001*IT_0095;
    const ccomplex_t IT_1029 = s_13*IT_0001*IT_0095*IT_0471;
    const ccomplex_t IT_1030 = IT_1028 + IT_1029;
    const ccomplex_t IT_1031 = s_12*s_24*s_34*IT_0001*IT_0095;
    const ccomplex_t IT_1032 = s_13*IT_0010*IT_0095;
    const ccomplex_t IT_1033 = -IT_1031 + -IT_1032;
    const ccomplex_t IT_1034 = IT_1030 + IT_1033;
    const ccomplex_t IT_1035 = IT_0002*IT_1034;
    const ccomplex_t IT_1036 = (-0.25)*IT_1035;
    const ccomplex_t IT_1037 = (-0.5)*IT_0534;
    const ccomplex_t IT_1038 = 0.5*IT_0537;
    const ccomplex_t IT_1039 = (-0.5)*IT_0539;
    const ccomplex_t IT_1040 = IT_0533 + IT_1027 + IT_1036 + IT_1037 + IT_1038
       + IT_1039;
    const ccomplex_t IT_1041 = (-8)*IT_0337;
    const ccomplex_t IT_1042 = s_12*IT_0001*IT_0003;
    const ccomplex_t IT_1043 = -IT_1042;
    const ccomplex_t IT_1044 = s_12*IT_1000;
    const ccomplex_t IT_1045 = s_14*s_23*s_34*IT_0001*IT_0095;
    const ccomplex_t IT_1046 = s_13*s_24*s_34*IT_0001*IT_0095;
    const ccomplex_t IT_1047 = IT_1044 + IT_1045 + IT_1046;
    const ccomplex_t IT_1048 = s_13*s_23*IT_0095;
    const ccomplex_t IT_1049 = s_12*IT_0001*IT_0003*IT_0095;
    const ccomplex_t IT_1050 = (-2)*IT_1048 + -IT_1049;
    const ccomplex_t IT_1051 = IT_1047 + IT_1050;
    const ccomplex_t IT_1052 = IT_0002*IT_1051;
    const ccomplex_t IT_1053 = (-0.5)*IT_1052;
    const ccomplex_t IT_1054 = -IT_0917;
    const ccomplex_t IT_1055 = IT_0447 + IT_0810 + IT_1043 + IT_1053 + IT_1054;
    const ccomplex_t IT_1056 = (-4)*IT_0493;
    const ccomplex_t IT_1057 = IT_0522 + IT_1015;
    const ccomplex_t IT_1058 = s_12*s_13*IT_0008;
    const ccomplex_t IT_1059 = s_23*IT_0323;
    const ccomplex_t IT_1060 = s_12*s_14*s_34;
    const ccomplex_t IT_1061 = s_23*IT_0008*IT_0009;
    const ccomplex_t IT_1062 = 2*IT_1059 + (-2)*IT_1060 + -IT_1061;
    const ccomplex_t IT_1063 = IT_1058 + IT_1062;
    const ccomplex_t IT_1064 = IT_0001*IT_1063;
    const ccomplex_t IT_1065 = (-0.5)*IT_1064;
    const ccomplex_t IT_1066 = IT_1057 + IT_1065;
    const ccomplex_t IT_1067 = (-4)*IT_0332;
    const ccomplex_t IT_1068 = s_12*IT_0008*IT_1000;
    const ccomplex_t IT_1069 = s_14*s_24*IT_1000;
    const ccomplex_t IT_1070 = (-2)*IT_1069;
    const ccomplex_t IT_1071 = IT_1068 + IT_1070;
    const ccomplex_t IT_1072 = IT_0365*IT_1071;
    const ccomplex_t IT_1073 = 0.5*IT_1072;
    const ccomplex_t IT_1074 = (-4)*IT_0448;
    const ccomplex_t IT_1075 = IT_0447 + IT_1073 + IT_1074;
    const ccomplex_t IT_1076 = IT_0001*IT_0933;
    const ccomplex_t IT_1077 = -IT_1076;
    const ccomplex_t IT_1078 = IT_1075 + IT_1077;
    const ccomplex_t IT_1079 = (-4)*IT_0269;
    const ccomplex_t IT_1080 = s_24*s_34*IT_0095;
    const ccomplex_t IT_1081 = (-2)*IT_1080;
    const ccomplex_t IT_1082 = IT_0429 + IT_1081;
    const ccomplex_t IT_1083 = IT_0365*IT_1082;
    const ccomplex_t IT_1084 = (-0.666666666666667)*IT_0454 + (
      -0.166666666666667)*IT_1083;
    const ccomplex_t IT_1085 = s_23 + IT_1084;
    const ccomplex_t IT_1086 = m_C_2*IT_1085;
    const ccomplex_t IT_1087 = (-6)*IT_1086;
    const ccomplex_t IT_1088 = (-2)*IT_0858;
    const ccomplex_t IT_1089 = IT_0008*IT_1000;
    const ccomplex_t IT_1090 = IT_0365*IT_1089;
    const ccomplex_t IT_1091 = 0.125*IT_1090;
    const ccomplex_t IT_1092 = IT_0095 + IT_1091;
    const ccomplex_t IT_1093 = IT_0000*IT_1092;
    const ccomplex_t IT_1094 = 8*IT_1093;
    const ccomplex_t IT_1095 = (-2)*IT_0124;
    const ccomplex_t IT_1096 = m_C_2*m_N_4*IT_0001;
    const ccomplex_t IT_1097 = IT_0002*IT_1089;
    const ccomplex_t IT_1098 = (-0.25)*IT_1097;
    const ccomplex_t IT_1099 = IT_0003 + IT_1098;
    const ccomplex_t IT_1100 = IT_1096*IT_1099;
    const ccomplex_t IT_1101 = (-4)*IT_1100;
    const ccomplex_t IT_1102 = (-2)*IT_0143;
    const ccomplex_t IT_1103 = IT_0994*IT_0997 + IT_0492*IT_1006 + IT_0877
      *IT_1013 + IT_0774*IT_1025 + IT_1040*IT_1041 + IT_1055*IT_1056 + IT_1066
      *IT_1067 + IT_1078*IT_1079 + IT_1087*IT_1088 + IT_1094*IT_1095 + IT_1101
      *IT_1102;
    const ccomplex_t IT_1104 = (-4)*IT_0007;
    const ccomplex_t IT_1105 = -IT_0060;
    const ccomplex_t IT_1106 = -conj(IT_0060);
    const ccomplex_t IT_1107 = -IT_0088;
    const ccomplex_t IT_1108 = -conj(IT_0088);
    const ccomplex_t IT_1109 = -conj(IT_0673);
    const ccomplex_t IT_1110 = conj(IT_0713) + IT_1109;
    const ccomplex_t IT_1111 = -IT_0050;
    const ccomplex_t IT_1112 = -conj(IT_0050);
    const ccomplex_t IT_1113 = conj(IT_0060)*IT_0421 + IT_0060*conj(IT_0421) +
       conj(IT_0088)*IT_0713 + IT_0088*conj(IT_0713) + conj(IT_0082)*(IT_0413 +
       IT_0418 + -IT_0421 + IT_0858) + IT_0082*(conj(IT_0413) + conj(IT_0418) + 
      -conj(IT_0421) + conj(IT_0858)) + (conj(IT_0413) + conj(IT_0418) + conj
      (IT_0858))*IT_1105 + (IT_0413 + IT_0418 + IT_0858)*IT_1106 + (conj(IT_0390
      ) + conj(IT_0395) + conj(IT_0673))*IT_1107 + (IT_0390 + IT_0395 + IT_0673)
      *IT_1108 + -(conj(IT_0390) + conj(IT_0395) + -IT_1110)*IT_1111 + -(IT_0390
       + IT_0395 + IT_0673 + -IT_0713)*IT_1112;
    const ccomplex_t IT_1114 = 8*IT_0713;
    const ccomplex_t IT_1115 = 0.125*IT_0800;
    const ccomplex_t IT_1116 = 0.125*IT_0735;
    const ccomplex_t IT_1117 = 0.125*IT_0802;
    const ccomplex_t IT_1118 = 0.125*IT_0682;
    const ccomplex_t IT_1119 = 0.125*IT_0725;
    const ccomplex_t IT_1120 = 0.125*IT_0808;
    const ccomplex_t IT_1121 = (IT_0060 + -IT_0082)*IT_0697 + IT_0421*IT_0791 
      + 0.125*IT_0711*IT_1114 + IT_0483*IT_1115 + IT_0635*IT_1116 + IT_0492
      *IT_1117 + IT_0665*IT_1118 + IT_0867*IT_1119 + IT_0877*IT_1120;
    const ccomplex_t IT_1122 = 2*conj(IT_0435);
    const ccomplex_t IT_1123 = 2*conj(IT_0545);
    const ccomplex_t IT_1124 = conj(IT_0088) + IT_1112;
    const ccomplex_t IT_1125 = 0.125*IT_0950;
    const ccomplex_t IT_1126 = 0.125*IT_0910;
    const ccomplex_t IT_1127 = 0.125*IT_0912;
    const ccomplex_t IT_1128 = 0.125*IT_0891;
    const ccomplex_t IT_1129 = 0.125*IT_0908;
    const ccomplex_t IT_1130 = 0.125*IT_0916;
    const ccomplex_t IT_1131 = IT_0088 + IT_1111;
    const ccomplex_t IT_1132 = 0.125*IT_0791;
    const ccomplex_t IT_1133 = 2*conj(IT_0485);
    const ccomplex_t IT_1134 = IT_0973 + IT_1133;
    const ccomplex_t IT_1135 = 8*conj(IT_0050);
    const ccomplex_t IT_1136 = (-8)*conj(IT_0088);
    const ccomplex_t IT_1137 = IT_1135 + IT_1136;
    const ccomplex_t IT_1138 = (-2)*conj(IT_0345);
    const ccomplex_t IT_1139 = IT_0987 + IT_1138;
    const ccomplex_t IT_1140 = 4*conj(IT_0342);
    const ccomplex_t IT_1141 = (-2)*conj(IT_0493);
    const ccomplex_t IT_1142 = 2*conj(IT_0332);
    const ccomplex_t IT_1143 = (-2)*conj(IT_0743);
    const ccomplex_t IT_1144 = (-2)*conj(IT_0548);
    const ccomplex_t IT_1145 = 8*conj(IT_0413);
    const ccomplex_t IT_1146 = IT_0845 + IT_1145;
    const ccomplex_t IT_1147 = 0.125*IT_0783;
    const ccomplex_t IT_1148 = 0.5*IT_0791*(conj(IT_0395) + 0.125*IT_0843) +
       0.0625*IT_0817*IT_0846 + IT_0862*IT_1110 + 0.5*IT_0849*IT_1115 + 0.5
      *IT_0496*IT_1116 + 0.0625*IT_0753*IT_1122 + 0.0625*IT_0763*IT_1123 +
       0.0625*IT_0742*IT_1134 + 0.0625*IT_0697*IT_1137 + 0.0625*IT_0827*IT_1139 
      + 0.0625*IT_0773*IT_1140 + 0.5*IT_1117*IT_1141 + 0.5*IT_1118*IT_1142 + 0.5
      *IT_1119*IT_1143 + 0.5*IT_1120*IT_1144 + 0.5*IT_1146*IT_1147;
    const ccomplex_t IT_1149 = 0.125*IT_0872*IT_0897 + 0.125*IT_0862*IT_0900 +
       0.125*IT_0494*IT_0924 + 0.125*IT_0744*IT_0931 + IT_0866*IT_0936 + 0.125
      *IT_0943*IT_0944 + 0.125*IT_0881*IT_0959 + 0.125*IT_0828*IT_0961 + 0.125
      *IT_0882*IT_0969 + IT_0635*IT_1125 + IT_0483*IT_1126 + IT_0492*IT_1127 +
       IT_0684*IT_1128 + IT_0727*IT_1129 + IT_0877*IT_1130 + IT_0715*IT_1132 +
       IT_0879*IT_1147;
    const ccomplex_t IT_1150 = 0.166666666666667*IT_0617;
    const ccomplex_t IT_1151 = IT_0315*IT_1150;
    const ccomplex_t IT_1152 = 0.166666666666667*IT_0983;
    const ccomplex_t IT_1153 = 0.166666666666667*IT_0985;
    const ccomplex_t IT_1154 = 0.166666666666667*IT_0986;
    const ccomplex_t IT_1155 = 0.166666666666667*IT_0982;
    const ccomplex_t IT_1156 = (-0.25)*IT_0329*IT_0335 + (-0.25)*IT_0348
      *IT_0623 + (-0.25)*IT_0219*IT_0626 + IT_0713*IT_0773 + (-0.25)*IT_0628
      *IT_0774 + (-0.25)*IT_0634*IT_0818 + (-0.25)*IT_0422*IT_0924 + -IT_0421
      *IT_0943 + 3*IT_0342*((-0.5)*IT_0592 + IT_0981) + (-0.25)*IT_0817*IT_0995 
      + (-1.5)*IT_1151 + (-1.5)*IT_0651*IT_1152 + (-1.5)*IT_0654*IT_1153 + (-1.5
      )*IT_0665*IT_1154 + (-1.5)*IT_0867*IT_1155;
    const ccomplex_t IT_1157 = IT_0000*IT_0824;
    const ccomplex_t IT_1158 = IT_0527 + IT_0529 + IT_0531 + IT_0656 + IT_0657
       + IT_0659 + IT_1014;
    const ccomplex_t IT_1159 = (-2)*IT_0713;
    const ccomplex_t IT_1160 = 0.166666666666667*IT_0655;
    const ccomplex_t IT_1161 = 0.166666666666667*IT_0974;
    const ccomplex_t IT_1162 = 0.166666666666667*IT_0664;
    const ccomplex_t IT_1163 = 0.166666666666667*IT_0650;
    const ccomplex_t IT_1164 = 0.166666666666667*IT_0662;
    const ccomplex_t IT_1165 = 0.166666666666667*IT_0422*IT_0891 + IT_0494
      *IT_1152 + IT_0818*IT_1154 + IT_0332*IT_1157 + 0.166666666666667*IT_0877
      *IT_1158 + 0.166666666666667*IT_0682*IT_1159 + IT_0335*IT_1160 + IT_0754
      *IT_1161 + IT_0876*IT_1162 + IT_0867*IT_1163 + IT_0492*IT_1164;
    const ccomplex_t IT_1166 = s_12 + IT_0782;
    const ccomplex_t IT_1167 = IT_0000*IT_1166;
    const ccomplex_t IT_1168 = IT_0344 + IT_0345;
    const ccomplex_t IT_1169 = IT_0492 + IT_0877;
    const ccomplex_t IT_1170 = 2*IT_0421;
    const ccomplex_t IT_1171 = 0.166666666666667*IT_0961;
    const ccomplex_t IT_1172 = 0.166666666666667*IT_0827;
    const ccomplex_t IT_1173 = -IT_0483*IT_1152 + -IT_0635*IT_1154 + -IT_0348
      *IT_1160 + -IT_0828*IT_1162 + -IT_0754*IT_1164 + IT_0482*IT_1167 + 
      -IT_1157*IT_1168 + -IT_1161*IT_1169 + -IT_1170*IT_1171 + -IT_0995*IT_1172;
    const ccomplex_t IT_1174 = IT_1122 + IT_1123;
    const ccomplex_t IT_1175 = 6*IT_0332*(IT_0458*IT_1152 + IT_0846*IT_1154 +
       conj(IT_0344)*IT_1157 + IT_0351*IT_1160 + IT_0840*IT_1164 + -conj(IT_0482
      )*IT_1167 + IT_0975*IT_1171 + IT_0984*IT_1172 + IT_1161*IT_1174);
    const ccomplex_t IT_1176 = IT_0000*IT_0928;
    const ccomplex_t IT_1177 = IT_0482 + IT_0485;
    const ccomplex_t IT_1178 = -IT_0534;
    const ccomplex_t IT_1179 = IT_0533 + IT_1178;
    const ccomplex_t IT_1180 = s_12*s_23*IT_0008;
    const ccomplex_t IT_1181 = s_14*s_23*s_24;
    const ccomplex_t IT_1182 = s_13*IT_0008*IT_0010;
    const ccomplex_t IT_1183 = s_14*s_34*IT_0010;
    const ccomplex_t IT_1184 = (-2)*IT_1181 + -IT_1182 + 2*IT_1183;
    const ccomplex_t IT_1185 = IT_1180 + IT_1184;
    const ccomplex_t IT_1186 = IT_0001*IT_1185;
    const ccomplex_t IT_1187 = (-0.5)*IT_1186;
    const ccomplex_t IT_1188 = IT_1179 + IT_1187;
    const ccomplex_t IT_1189 = 0.166666666666667*IT_0636;
    const ccomplex_t IT_1190 = 0.166666666666667*IT_0972;
    const ccomplex_t IT_1191 = 0.166666666666667*IT_0931;
    const ccomplex_t IT_1192 = 0.166666666666667*IT_0742;
    const ccomplex_t IT_1193 = 0.166666666666667*IT_0481;
    const ccomplex_t IT_1194 = 0.166666666666667*IT_0643;
    const ccomplex_t IT_1195 = 0.166666666666667*IT_0653*IT_0880 + IT_1151 +
       IT_0635*IT_1153 + IT_0483*IT_1155 + IT_0878*IT_1163 + -IT_1167*IT_1168 +
       IT_1176*IT_1177 + 0.666666666666667*IT_0316*IT_1188 + IT_0219*IT_1189 +
       IT_1169*IT_1190 + IT_0422*IT_1191 + IT_1159*IT_1192 + IT_0494*IT_1193 +
       IT_0754*IT_1194;
    const ccomplex_t IT_1196 = conj(IT_0482) + conj(IT_0485);
    const ccomplex_t IT_1197 = conj(IT_0344) + conj(IT_0345);
    const ccomplex_t IT_1198 = -IT_1197;
    const ccomplex_t IT_1199 = IT_0320*IT_1150 + IT_0846*IT_1153 + IT_0458
      *IT_1155 + IT_1142*IT_1163 + 0.666666666666667*IT_0311*IT_1188 + IT_0223
      *IT_1189 + IT_1174*IT_1190 + IT_0427*IT_1191 + IT_0979*IT_1192 + IT_0496
      *IT_1193 + IT_0840*IT_1194 + IT_1176*IT_1196 + IT_1167*IT_1198;
    const ccomplex_t IT_1200 = (-2)*IT_1048;
    const ccomplex_t IT_1201 = IT_1044 + IT_1200;
    const ccomplex_t IT_1202 = IT_0002*IT_1201;
    const ccomplex_t IT_1203 = 0.5*IT_1202;
    const ccomplex_t IT_1204 = s_14*s_23*IT_0095;
    const ccomplex_t IT_1205 = s_13*s_24*IT_0095;
    const ccomplex_t IT_1206 = IT_1204 + IT_1205;
    const ccomplex_t IT_1207 = s_12*s_34*IT_0095;
    const ccomplex_t IT_1208 = -IT_1207;
    const ccomplex_t IT_1209 = IT_1206 + IT_1208;
    const ccomplex_t IT_1210 = IT_0098*IT_1209;
    const ccomplex_t IT_1211 = 0.5*IT_1210;
    const ccomplex_t IT_1212 = (-2)*IT_0702 + 2*IT_0703 + 2*IT_0775 + (-2)
      *IT_0776;
    const ccomplex_t IT_1213 = IT_0701 + IT_1212;
    const ccomplex_t IT_1214 = IT_0001*IT_1213;
    const ccomplex_t IT_1215 = (-0.5)*IT_1214;
    const ccomplex_t IT_1216 = IT_0447 + IT_0449 + IT_1203 + IT_1211 + IT_1215;
    const ccomplex_t IT_1217 = conj(IT_0545)*IT_1216;
    const ccomplex_t IT_1218 = 0.25*IT_1006;
    const ccomplex_t IT_1219 = IT_0840*IT_1218;
    const ccomplex_t IT_1220 = 0.25*IT_1013;
    const ccomplex_t IT_1221 = IT_0841*IT_1220;
    const ccomplex_t IT_1222 = (-2)*conj(IT_0673);
    const ccomplex_t IT_1223 = IT_0984 + IT_1222;
    const ccomplex_t IT_1224 = 0.25*IT_1087;
    const ccomplex_t IT_1225 = (-2)*conj(IT_0858);
    const ccomplex_t IT_1226 = 0.25*IT_0994;
    const ccomplex_t IT_1227 = (-4)*conj(IT_0435);
    const ccomplex_t IT_1228 = 0.25*IT_1055;
    const ccomplex_t IT_1229 = (-4)*conj(IT_0344);
    const ccomplex_t IT_1230 = (-4)*conj(IT_0345);
    const ccomplex_t IT_1231 = IT_1229 + IT_1230;
    const ccomplex_t IT_1232 = 0.25*IT_1066;
    const ccomplex_t IT_1233 = 0.25*IT_1025;
    const ccomplex_t IT_1234 = (-8)*conj(IT_0342);
    const ccomplex_t IT_1235 = 0.25*IT_1040;
    const ccomplex_t IT_1236 = s_13*s_14*s_24;
    const ccomplex_t IT_1237 = s_24*s_34*IT_0009;
    const ccomplex_t IT_1238 = -IT_1061 + (-2)*IT_1236 + 2*IT_1237;
    const ccomplex_t IT_1239 = IT_1058 + IT_1238;
    const ccomplex_t IT_1240 = IT_0001*IT_1239;
    const ccomplex_t IT_1241 = (-0.5)*IT_1240;
    const ccomplex_t IT_1242 = IT_1057 + IT_1241;
    const ccomplex_t IT_1243 = IT_0365*IT_0366;
    const ccomplex_t IT_1244 = (-10)*s_13;
    const ccomplex_t IT_1245 = IT_1243 + IT_1244;
    const ccomplex_t IT_1246 = m_N_4*IT_1245;
    const ccomplex_t IT_1247 = (-8)*conj(IT_0124);
    const ccomplex_t IT_1248 = (-0.125)*IT_1247;
    const ccomplex_t IT_1249 = conj(IT_0215) + IT_1248;
    const ccomplex_t IT_1250 = 0.25*IT_1072;
    const ccomplex_t IT_1251 = IT_0447 + IT_1250;
    const ccomplex_t IT_1252 = (-2)*IT_1251;
    const ccomplex_t IT_1253 = (-0.666666666666667)*IT_0105 +
       0.166666666666667*IT_1083;
    const ccomplex_t IT_1254 = s_23 + IT_1253;
    const ccomplex_t IT_1255 = m_C_2*IT_1254;
    const ccomplex_t IT_1256 = 6*IT_1255;
    const ccomplex_t IT_1257 = 0.25*IT_1256;
    const ccomplex_t IT_1258 = (-0.166666666666667)*IT_0370 + (
      -0.666666666666667)*IT_0446;
    const ccomplex_t IT_1259 = s_13 + IT_1258;
    const ccomplex_t IT_1260 = m_N_4*IT_1259;
    const ccomplex_t IT_1261 = (-6)*IT_1260;
    const ccomplex_t IT_1262 = 0.25*IT_1261;
    const ccomplex_t IT_1263 = m_C_2*IT_0001;
    const ccomplex_t IT_1264 = IT_0002*IT_0429;
    const ccomplex_t IT_1265 = (-2)*IT_0102 + (-0.5)*IT_1264;
    const ccomplex_t IT_1266 = IT_0101 + IT_1265;
    const ccomplex_t IT_1267 = IT_1263*IT_1266;
    const ccomplex_t IT_1268 = (-2)*IT_1267;
    const ccomplex_t IT_1269 = 0.25*IT_1268;
    const ccomplex_t IT_1270 = -IT_0810;
    const ccomplex_t IT_1271 = IT_0447 + IT_0917 + IT_1043 + IT_1053 + IT_1270;
    const ccomplex_t IT_1272 = 0.25*IT_1271;
    const ccomplex_t IT_1273 = (-0.5)*IT_1014;
    const ccomplex_t IT_1274 = (-0.25)*IT_1023;
    const ccomplex_t IT_1275 = (-0.5)*IT_0523;
    const ccomplex_t IT_1276 = (-0.5)*IT_0526;
    const ccomplex_t IT_1277 = 0.5*IT_0528;
    const ccomplex_t IT_1278 = IT_0522 + IT_1273 + IT_1274 + IT_1275 + IT_1276
       + IT_1277;
    const ccomplex_t IT_1279 = 0.25*IT_1278;
    const ccomplex_t IT_1280 = (-4)*conj(IT_0482);
    const ccomplex_t IT_1281 = (-4)*conj(IT_0485);
    const ccomplex_t IT_1282 = IT_1280 + IT_1281;
    const ccomplex_t IT_1283 = s_13*IT_0471;
    const ccomplex_t IT_1284 = s_12*s_24*s_34;
    const ccomplex_t IT_1285 = -IT_1182 + 2*IT_1283 + (-2)*IT_1284;
    const ccomplex_t IT_1286 = IT_1180 + IT_1285;
    const ccomplex_t IT_1287 = IT_0001*IT_1286;
    const ccomplex_t IT_1288 = (-0.5)*IT_1287;
    const ccomplex_t IT_1289 = IT_1179 + IT_1288;
    const ccomplex_t IT_1290 = 0.25*IT_1289;
    const ccomplex_t IT_1291 = (-4)*conj(IT_0342);
    const ccomplex_t IT_1292 = -IT_1026;
    const ccomplex_t IT_1293 = (-0.5)*IT_1035;
    const ccomplex_t IT_1294 = IT_0534 + IT_0538 + IT_0539 + IT_1292 + IT_1293;
    const ccomplex_t IT_1295 = 0.25*IT_1294;
    const ccomplex_t IT_1296 = IT_1217 + IT_1219 + IT_1221 + IT_1197*IT_1242 +
       0.25*IT_1225*IT_1246 + IT_1249*IT_1252 + IT_1223*IT_1257 + IT_0428
      *IT_1262 + IT_0400*IT_1269 + IT_1227*IT_1272 + IT_1234*IT_1279 + IT_1282
      *IT_1290 + IT_1291*IT_1295;
    const ccomplex_t IT_1297 = (-8)*conj(IT_0269);
    const ccomplex_t IT_1298 = (-8)*conj(IT_0291);
    const ccomplex_t IT_1299 = (-8)*conj(IT_0310);
    const ccomplex_t IT_1300 = 0.666666666666667*IT_0353 + (-0.166666666666667
      )*IT_1243;
    const ccomplex_t IT_1301 = s_13 + IT_1300;
    const ccomplex_t IT_1302 = m_N_4*IT_1301;
    const ccomplex_t IT_1303 = (-6)*IT_1302;
    const ccomplex_t IT_1304 = (-2)*conj(IT_0143);
    const ccomplex_t IT_1305 = (-2)*conj(IT_0155);
    const ccomplex_t IT_1306 = (-2)*conj(IT_0124);
    const ccomplex_t IT_1307 = (-2)*conj(IT_0215);
    const ccomplex_t IT_1308 = (-4)*conj(IT_0493);
    const ccomplex_t IT_1309 = (-8)*conj(IT_0337);
    const ccomplex_t IT_1310 = (-4)*conj(IT_0743);
    const ccomplex_t IT_1311 = (-4)*conj(IT_0337);
    const ccomplex_t IT_1312 = conj(IT_0548)*IT_1216 + IT_1122*IT_1218 +
       IT_1123*IT_1220 + conj(IT_0332)*IT_1242 + 0.25*IT_0984*IT_1246 + IT_1225
      *IT_1257 + IT_0400*IT_1262 + IT_0428*IT_1269 + 0.25*IT_1251*(IT_1297 +
       IT_1298 + IT_1299) + 0.25*IT_1222*IT_1303 + 0.25*IT_1094*(IT_1304 +
       IT_1305) + 0.25*IT_1101*(IT_1306 + IT_1307) + IT_1272*IT_1308 + IT_1279
      *IT_1309 + IT_1290*IT_1310 + IT_1295*IT_1311;
    const ccomplex_t IT_1313 = s_23*m_C_2;
    const ccomplex_t IT_1314 = IT_0006*IT_1313;
    const ccomplex_t IT_1315 = s_13*m_N_4;
    const ccomplex_t IT_1316 = IT_0006*IT_1315;
    const ccomplex_t IT_1317 = (1 + 0.5*IT_0004)*(IT_0447 + (-2)*IT_0448);
    const ccomplex_t IT_1318 = (-4)*IT_0545;
    const ccomplex_t IT_1319 = s_34*IT_0001;
    const ccomplex_t IT_1320 = IT_0690 + IT_0691;
    const ccomplex_t IT_1321 = s_13*s_23*s_34*IT_0002;
    const ccomplex_t IT_1322 = (-2)*IT_1321;
    const ccomplex_t IT_1323 = IT_1320 + IT_1322;
    const ccomplex_t IT_1324 = IT_1319*IT_1323;
    const ccomplex_t IT_1325 = 0.25*IT_1324;
    const ccomplex_t IT_1326 = 0.25*IT_1158;
    const ccomplex_t IT_1327 = IT_0538 + IT_0540 + IT_0542 + IT_0637 + IT_0638
       + IT_0640 + IT_1026;
    const ccomplex_t IT_1328 = 0.25*IT_1327;
    const ccomplex_t IT_1329 = m_C_2*m_N_4*IT_0095;
    const ccomplex_t IT_1330 = IT_0006*IT_1329;
    const ccomplex_t IT_1331 = 0.25*IT_0869*IT_0972 + 0.25*IT_0684*IT_0974 +
       0.25*IT_0763*IT_1159 + 0.25*IT_0969*IT_1170 + IT_0316*IT_1216 + (IT_0050 
      + 0.25*IT_0611)*IT_1314 + IT_0084*IT_1316 + IT_0881*IT_1325 + IT_0878
      *IT_1326 + IT_0744*IT_1328 + 0.25*IT_1318*IT_1330;
    const ccomplex_t IT_1332 = (-0.25)*IT_0832*IT_0974 + (-0.25)*IT_0969
      *IT_0975 + (-0.25)*IT_0763*IT_0979 + (-0.25)*IT_0972*IT_1134 + -IT_0311
      *IT_1216 + -IT_0610*IT_1314 + -IT_0086*IT_1316 + conj(IT_0548)*IT_1317 + 
      -IT_1141*IT_1325 + -IT_1142*IT_1326 + -IT_0839*IT_1328 + conj(IT_0545)
      *IT_1330;
    const ccomplex_t IT_1333 = IT_0001*IT_0785;
    const ccomplex_t IT_1334 = 12*IT_0095 + (-4)*IT_1333;
    const ccomplex_t IT_1335 = IT_1090 + IT_1334;
    const ccomplex_t IT_1336 = IT_0000*IT_1335;
    const ccomplex_t IT_1337 = (-2)*conj(IT_0269);
    const ccomplex_t IT_1338 = (-2)*conj(IT_0291);
    const ccomplex_t IT_1339 = (-2)*conj(IT_0310);
    const ccomplex_t IT_1340 = IT_1337 + IT_1338 + IT_1339;
    const ccomplex_t IT_1341 = -IT_1217;
    const ccomplex_t IT_1342 = (-0.25)*IT_1006;
    const ccomplex_t IT_1343 = IT_0840*IT_1342;
    const ccomplex_t IT_1344 = (-0.25)*IT_1013;
    const ccomplex_t IT_1345 = IT_0841*IT_1344;
    const ccomplex_t IT_1346 = (-0.25)*IT_1087;
    const ccomplex_t IT_1347 = (-0.25)*IT_0994;
    const ccomplex_t IT_1348 = (-0.25)*IT_1055;
    const ccomplex_t IT_1349 = (-0.25)*IT_1066;
    const ccomplex_t IT_1350 = (-0.25)*IT_1025;
    const ccomplex_t IT_1351 = (-0.25)*IT_1040;
    const ccomplex_t IT_1352 = -IT_0401*IT_1078 + IT_1188*IT_1196 + 0.25
      *IT_1336*IT_1340 + -IT_1341 + -IT_1343 + -IT_1345 + -IT_1223*IT_1346 + 
      -IT_1225*IT_1347 + -IT_1227*IT_1348 + -IT_1231*IT_1349 + -IT_1140*IT_1350 
      + -IT_1234*IT_1351;
    const ccomplex_t IT_1353 = (-4)*IT_0291;
    const ccomplex_t IT_1354 = (-0.25)*IT_1094;
    const ccomplex_t IT_1355 = (-0.25)*IT_1101;
    const ccomplex_t IT_1356 = IT_0492*IT_1342 + IT_0877*IT_1344 + IT_1088
      *IT_1346 + IT_0997*IT_1347 + IT_1056*IT_1348 + IT_1067*IT_1349 + IT_0774
      *IT_1350 + IT_1041*IT_1351 + IT_1078*(IT_0310 + (-0.25)*IT_1079 + (-0.25)
      *IT_1353) + IT_1095*IT_1354 + IT_1102*IT_1355;
    const ccomplex_t IT_1357 = IT_0001*IT_0706;
    const ccomplex_t IT_1358 = -IT_1357;
    const ccomplex_t IT_1359 = IT_1075 + IT_1358;
    const ccomplex_t IT_1360 = (-0.25)*IT_1256;
    const ccomplex_t IT_1361 = (-0.25)*IT_1261;
    const ccomplex_t IT_1362 = (-0.25)*IT_1268;
    const ccomplex_t IT_1363 = (-0.25)*IT_1303;
    const ccomplex_t IT_1364 = (-0.25)*IT_1271;
    const ccomplex_t IT_1365 = (-0.25)*IT_1278;
    const ccomplex_t IT_1366 = (-0.25)*IT_1289;
    const ccomplex_t IT_1367 = (-0.25)*IT_1294;
    const ccomplex_t IT_1368 = IT_1198*IT_1242 + 2*IT_1249*IT_1251 + IT_1341 +
       IT_1343 + IT_1345 + IT_1340*IT_1355 + (conj(IT_0143) + conj(IT_0155))
      *IT_1359 + IT_1223*IT_1360 + IT_0428*IT_1361 + IT_0400*IT_1362 + IT_1225
      *IT_1363 + IT_1227*IT_1364 + IT_1234*IT_1365 + IT_1282*IT_1366 + IT_1291
      *IT_1367;
    const ccomplex_t IT_1369 = (-8)*IT_0269;
    const ccomplex_t IT_1370 = (-8)*IT_0291;
    const ccomplex_t IT_1371 = (-8)*IT_0310;
    const ccomplex_t IT_1372 = (-4)*IT_0240;
    const ccomplex_t IT_1373 = (-2)*IT_0155;
    const ccomplex_t IT_1374 = (-2)*IT_0215;
    const ccomplex_t IT_1375 = (-4)*IT_0743;
    const ccomplex_t IT_1376 = IT_0548*IT_1216 + IT_0332*IT_1242 + 0.25
      *IT_0995*IT_1246 + -IT_0492*IT_1342 + -IT_0877*IT_1344 + -IT_1088*IT_1360 
      + -IT_0397*IT_1361 + -IT_0423*IT_1362 + -IT_0996*IT_1363 + -IT_1056
      *IT_1364 + -IT_1041*IT_1365 + -IT_0944*IT_1367 + 0.25*IT_1251*(IT_1369 +
       IT_1370 + IT_1371) + -IT_1359*(IT_0250 + (-0.25)*IT_1372) + -IT_1354*
      (IT_1102 + IT_1373) + -IT_1355*(IT_1095 + IT_1374) + -IT_1366*IT_1375;
    const ccomplex_t IT_1377 = (-0.166666666666667)*IT_0972;
    const ccomplex_t IT_1378 = (-0.166666666666667)*IT_0315*IT_0440 + (
      -0.166666666666667)*IT_0335*IT_0445 + IT_0342*IT_0601 + (
      -0.166666666666667)*IT_0491*IT_0754 + (-0.166666666666667)*IT_0651*IT_0974
       + (-0.166666666666667)*IT_0604*IT_0978 + (-0.166666666666667)*IT_0753
      *IT_1159 + (-0.166666666666667)*IT_0959*IT_1170 + IT_0867*IT_1377;
    const ccomplex_t IT_1379 = (-0.166666666666667)*IT_0636;
    const ccomplex_t IT_1380 = (-0.166666666666667)*IT_0985;
    const ccomplex_t IT_1381 = (-0.166666666666667)*IT_0982;
    const ccomplex_t IT_1382 = (-0.166666666666667)*IT_0650;
    const ccomplex_t IT_1383 = (-0.166666666666667)*IT_0653;
    const ccomplex_t IT_1384 = (-0.166666666666667)*IT_0725;
    const ccomplex_t IT_1385 = (-0.166666666666667)*IT_0908;
    const ccomplex_t IT_1386 = (-0.166666666666667)*IT_1327;
    const ccomplex_t IT_1387 = (-0.166666666666667)*IT_0643;
    const ccomplex_t IT_1388 = (-0.166666666666667)*IT_0219*IT_0617 + IT_0332
      *IT_1167 + (-0.666666666666667)*IT_0124*IT_1188 + IT_0754*IT_1377 +
       IT_0315*IT_1379 + IT_0818*IT_1380 + IT_0494*IT_1381 + IT_0684*IT_1382 +
       IT_0727*IT_1383 + IT_0995*IT_1384 + IT_1170*IT_1385 + IT_0877*IT_1386 +
       IT_0492*IT_1387;
    const ccomplex_t IT_1389 = (-0.166666666666667)*IT_0223*IT_0617 + conj
      (IT_0332)*IT_1167 + (-0.666666666666667)*conj(IT_0124)*IT_1188 + IT_0840
      *IT_1377 + IT_0320*IT_1379 + IT_0849*IT_1380 + IT_0496*IT_1381 + IT_0832
      *IT_1382 + IT_0836*IT_1383 + IT_0984*IT_1384 + IT_0975*IT_1385 + IT_1123
      *IT_1386 + IT_1122*IT_1387;
    const ccomplex_t IT_1390 = IT_0754*IT_1218;
    const ccomplex_t IT_1391 = IT_0764*IT_1220;
    const ccomplex_t IT_1392 = (-4)*IT_0124;
    const ccomplex_t IT_1393 = IT_1078*IT_1392;
    const ccomplex_t IT_1394 = 0.25*IT_1393;
    const ccomplex_t IT_1395 = (-4)*IT_0435;
    const ccomplex_t IT_1396 = (-4)*IT_0344;
    const ccomplex_t IT_1397 = (-4)*IT_0345;
    const ccomplex_t IT_1398 = IT_1396 + IT_1397;
    const ccomplex_t IT_1399 = (-8)*IT_0342;
    const ccomplex_t IT_1400 = IT_0545*IT_1216;
    const ccomplex_t IT_1401 = (-2)*IT_0269;
    const ccomplex_t IT_1402 = (-2)*IT_0291;
    const ccomplex_t IT_1403 = (-2)*IT_0310;
    const ccomplex_t IT_1404 = IT_1401 + IT_1402 + IT_1403;
    const ccomplex_t IT_1405 = IT_1177*IT_1188 + IT_0997*IT_1224 + IT_1088
      *IT_1226 + IT_0629*IT_1233 + IT_1390 + IT_1391 + IT_1394 + IT_1228*IT_1395
       + IT_1232*IT_1398 + IT_1235*IT_1399 + IT_1400 + 0.25*IT_1336*IT_1404;
    const ccomplex_t IT_1406 = (-4)*IT_0143;
    const ccomplex_t IT_1407 = (-8)*IT_0124;
    const ccomplex_t IT_1408 = IT_0215 + (-0.125)*IT_1407;
    const ccomplex_t IT_1409 = (-4)*IT_0482;
    const ccomplex_t IT_1410 = (-4)*IT_0485;
    const ccomplex_t IT_1411 = IT_1409 + IT_1410;
    const ccomplex_t IT_1412 = (-4)*IT_0342;
    const ccomplex_t IT_1413 = IT_1168*IT_1242;
    const ccomplex_t IT_1414 = IT_1088*IT_1246 + 4*IT_0997*IT_1257 + 4*IT_0423
      *IT_1262 + 4*IT_0397*IT_1269 + 4*IT_1390 + 4*IT_1391 + 4*IT_1272*IT_1395 +
       4*IT_1279*IT_1399 + 4*IT_1400 + IT_1359*IT_1406 + 4*IT_1252*IT_1408 + 4
      *IT_1290*IT_1411 + 4*IT_1295*IT_1412 + 4*IT_1413;
    const ccomplex_t IT_1415 = 4*IT_0997*IT_1257 + 4*IT_0423*IT_1262 + 4
      *IT_0397*IT_1269 + IT_1088*IT_1303 + 4*IT_1390 + 4*IT_1391 + 4*IT_1272
      *IT_1395 + 4*IT_1279*IT_1399 + 4*IT_1400 + IT_1101*IT_1404 + IT_1359
      *IT_1406 + 4*IT_1252*IT_1408 + 4*IT_1290*IT_1411 + 4*IT_1295*IT_1412 + 4
      *IT_1413;
    const ccomplex_t IT_1416 = IT_0223*IT_0617;
    const ccomplex_t IT_1417 = IT_0662*IT_1122;
    const ccomplex_t IT_1418 = IT_0340*IT_0655;
    const ccomplex_t IT_1419 = IT_0421*IT_0862;
    const ccomplex_t IT_1420 = 16*IT_1419;
    const ccomplex_t IT_1421 = IT_0818*IT_0950;
    const ccomplex_t IT_1422 = IT_0483*IT_0924;
    const ccomplex_t IT_1423 = IT_0494*IT_0910;
    const ccomplex_t IT_1424 = IT_0943*IT_1412;
    const ccomplex_t IT_1425 = IT_0754*IT_0912;
    const ccomplex_t IT_1426 = IT_0803*IT_0959;
    const ccomplex_t IT_1427 = IT_0878*IT_0891;
    const ccomplex_t IT_1428 = IT_0876*IT_0961;
    const ccomplex_t IT_1429 = IT_0880*IT_0908;
    const ccomplex_t IT_1430 = IT_0869*IT_0931;
    const ccomplex_t IT_1431 = IT_0764*IT_0916;
    const ccomplex_t IT_1432 = IT_0809*IT_0969;
    const ccomplex_t IT_1433 = IT_0700*IT_0897;
    const ccomplex_t IT_1434 = 8*conj(IT_0673);
    const ccomplex_t IT_1435 = IT_0844 + IT_1434;
    const ccomplex_t IT_1436 = IT_0849*IT_0950;
    const ccomplex_t IT_1437 = IT_0458*IT_0924;
    const ccomplex_t IT_1438 = IT_0496*IT_0910;
    const ccomplex_t IT_1439 = IT_0943*IT_1291;
    const ccomplex_t IT_1440 = IT_0840*IT_0912;
    const ccomplex_t IT_1441 = IT_0847*IT_0959;
    const ccomplex_t IT_1442 = IT_0891*IT_1142;
    const ccomplex_t IT_1443 = IT_0961*IT_1139;
    const ccomplex_t IT_1444 = IT_0908*IT_1143;
    const ccomplex_t IT_1445 = IT_0931*IT_1134;
    const ccomplex_t IT_1446 = IT_0841*IT_0916;
    const ccomplex_t IT_1447 = IT_0848*IT_0969;
    const ccomplex_t IT_1448 = 16*conj(IT_0421);
    const ccomplex_t IT_1449 = IT_0862*IT_1448;
    const ccomplex_t IT_1450 = IT_0835*IT_0897;
    const ccomplex_t IT_1451 = IT_0846*IT_0950;
    const ccomplex_t IT_1452 = IT_0496*IT_0924;
    const ccomplex_t IT_1453 = IT_0458*IT_0910;
    const ccomplex_t IT_1454 = IT_0943*IT_1311;
    const ccomplex_t IT_1455 = IT_0912*IT_1122;
    const ccomplex_t IT_1456 = IT_0959*IT_1141;
    const ccomplex_t IT_1457 = IT_0832*IT_0891;
    const ccomplex_t IT_1458 = IT_0850*IT_0961;
    const ccomplex_t IT_1459 = IT_0838*IT_0908;
    const ccomplex_t IT_1460 = IT_0839*IT_0931;
    const ccomplex_t IT_1461 = IT_0916*IT_1123;
    const ccomplex_t IT_1462 = IT_0969*IT_1144;
    const ccomplex_t IT_1463 = IT_0791*IT_1435;
    const ccomplex_t IT_1464 = (-16)*conj(IT_0390);
    const ccomplex_t IT_1465 = (-16)*conj(IT_0395);
    const ccomplex_t IT_1466 = IT_1464 + IT_1465;
    const ccomplex_t IT_1467 = IT_0862*IT_1466;
    const ccomplex_t IT_1468 = IT_0897*IT_1137;
    const ccomplex_t IT_1469 = IT_1006*IT_1122;
    const ccomplex_t IT_1470 = IT_1013*IT_1123;
    const ccomplex_t IT_1471 = IT_1087*IT_1225;
    const ccomplex_t IT_1472 = IT_0994*IT_1223;
    const ccomplex_t IT_1473 = IT_1094*IT_1306;
    const ccomplex_t IT_1474 = IT_1101*IT_1304;
    const ccomplex_t IT_1475 = IT_1055*IT_1308;
    const ccomplex_t IT_1476 = (-4)*conj(IT_0332);
    const ccomplex_t IT_1477 = IT_1066*IT_1476;
    const ccomplex_t IT_1478 = IT_0842*IT_1025;
    const ccomplex_t IT_1479 = IT_1040*IT_1309;
    const ccomplex_t IT_1480 = (-4)*conj(IT_0269);
    const ccomplex_t IT_1481 = 4*IT_0007*(conj(IT_0050)*IT_0084 + IT_0050
      *IT_0086 + conj(IT_0088)*IT_0091 + IT_0088*IT_0094) + IT_0109*(IT_0216
      *IT_0219 + IT_0220*IT_0223 + IT_0312*IT_0315 + IT_0317*IT_0320) + IT_0329*
      (conj(IT_0332)*IT_0335 + (IT_0332 + IT_0337)*IT_0340 + (conj(IT_0342) +
       conj(IT_0344) + conj(IT_0345))*IT_0348 + (IT_0342 + IT_0344 + IT_0345)
      *IT_0351) + (IT_0312*IT_0335 + IT_0317*IT_0340 + IT_0216*IT_0348 + IT_0220
      *IT_0351)*IT_0364 + IT_0374*(IT_0311*IT_0397 + IT_0316*IT_0400 + IT_0401
      *IT_0423 + IT_0424*IT_0428) + (IT_0397*IT_0401 + IT_0311*IT_0423 + IT_0400
      *IT_0424 + IT_0316*IT_0428)*IT_0433 + IT_0435*(IT_0223*IT_0440 + IT_0351
      *IT_0445 + IT_0457*IT_0458) + IT_0481*(IT_0484 + IT_0486 + IT_0487 +
       IT_0488) + conj(IT_0435)*(IT_0219*IT_0440 + IT_0348*IT_0445 + IT_0457
      *IT_0483 + IT_0491*IT_0492) + IT_0457*(IT_0495 + IT_0497) + (IT_0507 +
       IT_0515)*(IT_0516 + IT_0517 + IT_0518 + IT_0519 + IT_0520 + IT_0521) + 
      (IT_0532 + IT_0543)*(IT_0546 + IT_0547 + IT_0549 + IT_0550) + (IT_0558 +
       IT_0566)*(IT_0567 + IT_0568 + IT_0569 + IT_0570 + IT_0571 + IT_0572 +
       IT_0573 + IT_0574 + IT_0575 + IT_0576 + IT_0577 + IT_0578 + IT_0579 +
       IT_0580 + IT_0581 + IT_0582 + IT_0583 + IT_0584) + IT_0592*IT_0594 +
       IT_0601*(IT_0603 + IT_0605 + IT_0607) + 4*IT_0608*(IT_0060*conj(IT_0060) 
      + IT_0088*conj(IT_0088) + 0.25*conj(IT_0082)*IT_0089 + IT_0082*IT_0094 +
       IT_0050*IT_0610 + 0.25*conj(IT_0050)*IT_0611) + conj(IT_0342)*(IT_0219
      *IT_0617 + IT_0335*IT_0623 + IT_0315*IT_0626 + IT_0628*IT_0629 + IT_0634
      *IT_0635) + conj(IT_0482)*(IT_0219*IT_0617 + IT_0315*IT_0636 + IT_0492
      *IT_0643 + IT_0650*IT_0651 + IT_0653*IT_0654) + conj(IT_0344)*(IT_0335
      *IT_0655 + IT_0492*IT_0662 + IT_0664*IT_0665) + conj(IT_0673)*IT_0829 +
       IT_0673*IT_0851 + conj(IT_0858)*IT_0883 + conj(IT_0413)*IT_0970 + IT_0493
      *(IT_0320*IT_0440 + IT_0340*IT_0445 + IT_0972*IT_0973 + IT_0830*IT_0974 +
       IT_0959*IT_0975 + IT_0606*IT_0978 + IT_0753*IT_0979) + IT_0337*IT_0988 +
       conj(IT_0269)*IT_1103 + IT_1104*IT_1113 + 8*conj(IT_0713)*IT_1121 + 8
      *IT_0713*((conj(IT_0060) + -conj(IT_0082))*IT_0697 + IT_0458*IT_1115 +
       IT_0846*IT_1116 + IT_0987*IT_1118 + IT_0973*IT_1119 + IT_1117*IT_1122 +
       IT_1120*IT_1123) + 8*IT_0421*(IT_0897*IT_1124 + IT_0458*IT_1125 + IT_0846
      *IT_1126 + IT_0847*IT_1127 + IT_0987*IT_1128 + IT_0973*IT_1129 + IT_0848
      *IT_1130) + 8*conj(IT_0421)*(IT_0421*IT_0936 + IT_0483*IT_1125 + IT_0635
      *IT_1126 + IT_0803*IT_1127 + IT_0665*IT_1128 + IT_0867*IT_1129 + IT_0809
      *IT_1130 + IT_0897*IT_1131 + IT_1114*IT_1132) + 16*IT_0858*IT_1148 + 8
      *conj(IT_0418)*IT_1149 + (-4)*conj(IT_0337)*IT_1156 + 6*conj(IT_0345)
      *IT_1165 + (-6)*conj(IT_0332)*IT_1173 + IT_1175 + 6*conj(IT_0743)*IT_1195 
      + 6*IT_0743*IT_1199 + 4*IT_0124*(conj(IT_0482)*IT_1188 + IT_1217 + IT_1219
       + IT_1221 + IT_1223*IT_1224 + IT_1225*IT_1226 + IT_1227*IT_1228 + IT_1231
      *IT_1232 + IT_1140*IT_1233 + IT_1234*IT_1235) + 4*IT_0143*IT_1296 + 4*
      (IT_0240 + IT_0250)*IT_1312 + 4*conj(IT_0545)*(IT_0091*IT_1314 + IT_1131
      *IT_1316 + 0.25*IT_1317*IT_1318 + IT_0803*IT_1325 + IT_0651*IT_1326 +
       IT_0867*IT_1328) + 4*IT_0545*(IT_0094*IT_1314 + IT_1124*IT_1316 + IT_0847
      *IT_1325 + IT_0830*IT_1326 + IT_0973*IT_1328) + 4*conj(IT_0548)*IT_1331 + 
      (-4)*IT_0548*IT_1332 + 4*IT_0215*IT_1352 + (-4)*(conj(IT_0291) + conj
      (IT_0310))*IT_1356 + (-4)*IT_0155*IT_1368 + 4*(conj(IT_0240) + conj
      (IT_0250))*IT_1376 + (-6)*conj(IT_0493)*IT_1378 + (-6)*conj(IT_0485)
      *IT_1388 + (-6)*IT_0485*IT_1389 + 4*conj(IT_0124)*(IT_0482*IT_1188 +
       IT_0997*IT_1224 + IT_1088*IT_1226 + IT_0629*IT_1233 + IT_1390 + IT_1391 +
       IT_1394 + IT_1228*IT_1395 + IT_1232*IT_1398 + IT_1235*IT_1399 + IT_1400) 
      + 4*conj(IT_0215)*IT_1405 + conj(IT_0143)*IT_1414 + conj(IT_0155)*IT_1415 
      + IT_0342*(IT_0340*IT_0623 + IT_0320*IT_0626 + IT_1416) + IT_0482*(IT_0320
      *IT_0636 + IT_0650*IT_0830 + IT_0643*IT_1122 + IT_1416) + IT_0344*(IT_1417
       + IT_1418) + IT_0345*(IT_0427*IT_0891 + IT_0650*IT_0973 + IT_0840*IT_0974
       + IT_0682*IT_0979 + IT_0496*IT_0983 + IT_0849*IT_0986 + IT_0664*IT_0987 +
       IT_1123*IT_1158 + IT_1417 + IT_1418) + conj(IT_0390)*(IT_0714*IT_0783 +
       IT_0784*IT_0936 + IT_1420 + IT_1421 + IT_1422 + IT_1423 + IT_1424 +
       IT_1425 + IT_1426 + IT_1427 + IT_1428 + IT_1429 + IT_1430 + IT_1431 +
       IT_1432 + IT_1433) + conj(IT_0395)*(IT_0715*IT_0783 + IT_0874*IT_0936 +
       IT_1420 + IT_1421 + IT_1422 + IT_1423 + IT_1424 + IT_1425 + IT_1426 +
       IT_1427 + IT_1428 + IT_1429 + IT_1430 + IT_1431 + IT_1432 + IT_1433) +
       IT_0390*(IT_0783*IT_0844 + IT_1436 + IT_1437 + IT_1438 + IT_1439 +
       IT_1440 + IT_1441 + IT_1442 + IT_1443 + IT_1444 + IT_1445 + IT_1446 +
       IT_1447 + IT_1449 + IT_1450) + IT_0395*(IT_0843*IT_0936 + IT_0783*IT_1435
       + IT_1436 + IT_1437 + IT_1438 + IT_1439 + IT_1440 + IT_1441 + IT_1442 +
       IT_1443 + IT_1444 + IT_1445 + IT_1446 + IT_1447 + IT_1449 + IT_1450) +
       IT_0413*(IT_0845*IT_0936 + IT_1451 + IT_1452 + IT_1453 + IT_1454 +
       IT_1455 + IT_1456 + IT_1457 + IT_1458 + IT_1459 + IT_1460 + IT_1461 +
       IT_1462 + IT_1463 + IT_1467 + IT_1468) + IT_0418*(IT_0936*IT_1146 +
       IT_1451 + IT_1452 + IT_1453 + IT_1454 + IT_1455 + IT_1456 + IT_1457 +
       IT_1458 + IT_1459 + IT_1460 + IT_1461 + IT_1462 + IT_1463 + IT_1467 +
       IT_1468) + IT_0269*(IT_1469 + IT_1470 + IT_1471 + IT_1472 + IT_1473 +
       IT_1474 + IT_1475 + IT_1477 + IT_1478 + IT_1479) + (IT_0291 + IT_0310)*
      (IT_1469 + IT_1470 + IT_1471 + IT_1472 + IT_1473 + IT_1474 + IT_1475 +
       IT_1477 + IT_1478 + IT_1479 + IT_1078*IT_1480);
    return create_ccomplex_return(IT_1481);
}

