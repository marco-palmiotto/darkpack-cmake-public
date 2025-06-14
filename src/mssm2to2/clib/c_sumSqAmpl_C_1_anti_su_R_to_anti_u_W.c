#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_su_R_to_anti_u_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_su_R_to_anti_u_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_ud = param->V_ud;
    const creal_t V_us = param->V_us;
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
    const creal_t m_su_R = param->m_su_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t V_ub_mod = param->V_ub_mod;
    const creal_t reg_prop = param->reg_prop;
    const creal_t delta_wolf = param->delta_wolf;
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
    const ccomplex_t IT_0000 = pow(m_b, 2);
    const ccomplex_t IT_0001 = pow(m_C_1, 2);
    const ccomplex_t IT_0002 = pow(m_su_R, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_12 + IT_0000 + -IT_0001 + -IT_0002 
      + -reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0006 = sin(beta);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = sin(theta_W);
    const ccomplex_t IT_0009 = cpow(IT_0008, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*e_em
      *IT_0004*IT_0005*IT_0007*IT_0009*V_ub_mod;
    const ccomplex_t IT_0011 = 0.5*IT_0010;
    const ccomplex_t IT_0012 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0013 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0009
      *IT_0012*V_ub_mod;
    const ccomplex_t IT_0014 = 0.5*IT_0013;
    const ccomplex_t IT_0015 = IT_0011*IT_0014;
    const ccomplex_t IT_0016 = IT_0003*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*IT_0016;
    const ccomplex_t IT_0018 = pow(m_d, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + IT_0018 
      + -reg_prop, -1);
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_ud
      *e_em*IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0021 = 0.5*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0009;
    const ccomplex_t IT_0023 = 0.5*IT_0022;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_s, 2);
    const ccomplex_t IT_0028 = cpow((-2)*s_12 + -IT_0001 + -IT_0002 + IT_0027 
      + -reg_prop, -1);
    const ccomplex_t IT_0029 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_us
      *e_em*IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0030 = 0.5*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0009;
    const ccomplex_t IT_0032 = 0.5*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = IT_0028*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*IT_0034;
    const ccomplex_t IT_0036 = 2*IT_0017 + 2*IT_0026 + 2*IT_0035;
    const ccomplex_t IT_0037 = pow(m_W, -2);
    const ccomplex_t IT_0038 = pow(s_34, 2);
    const ccomplex_t IT_0039 = IT_0037*IT_0038;
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = pow(m_u, 2);
    const ccomplex_t IT_0042 = IT_0040 + IT_0041;
    const ccomplex_t IT_0043 = s_13*IT_0042;
    const ccomplex_t IT_0044 = m_C_1*IT_0042;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0046 = 0.5*IT_0045;
    const ccomplex_t IT_0047 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0048 = IT_0009*IT_0047;
    const ccomplex_t IT_0049 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0050 = IT_0009*IT_0049;
    const ccomplex_t IT_0051 = 1.4142135623731*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + (-0.5)*IT_0051);
    const ccomplex_t IT_0053 = -IT_0052;
    const ccomplex_t IT_0054 = IT_0046*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0002 + IT_0041 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0059 = 0.5*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0061 = IT_0009*IT_0060;
    const ccomplex_t IT_0062 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0063 = IT_0009*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-0.5)*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0059*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0002 + IT_0041 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0069 = IT_0067*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*IT_0069;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0072 = 0.5*IT_0071;
    const ccomplex_t IT_0073 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0074 = IT_0009*IT_0073;
    const ccomplex_t IT_0075 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0076 = IT_0009*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*(IT_0074 + (-0.5)*IT_0077);
    const ccomplex_t IT_0079 = -IT_0078;
    const ccomplex_t IT_0080 = IT_0072*IT_0079;
    const ccomplex_t IT_0081 = cpow((-2)*s_23 + IT_0002 + IT_0041 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0082 = IT_0080*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*IT_0082;
    const ccomplex_t IT_0084 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0004*IT_0007*IT_0009;
    const ccomplex_t IT_0085 = 0.5*IT_0084;
    const ccomplex_t IT_0086 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0087 = IT_0009*IT_0086;
    const ccomplex_t IT_0088 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0089 = IT_0009*IT_0088;
    const ccomplex_t IT_0090 = 1.4142135623731*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0087 + (-0.5)*IT_0090);
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0085*IT_0092;
    const ccomplex_t IT_0094 = cpow((-2)*s_23 + IT_0002 + IT_0041 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = IT_0057 + IT_0070 + IT_0083 + IT_0096;
    const ccomplex_t IT_0098 = m_C_1*IT_0097;
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0101 = IT_0009*IT_0100;
    const ccomplex_t IT_0102 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0103 = IT_0009*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + 0.5*IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = IT_0046*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0109 = IT_0055*IT_0107*IT_0108;
    const ccomplex_t IT_0110 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0111 = IT_0009*IT_0110;
    const ccomplex_t IT_0112 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0113 = IT_0009*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*(IT_0111 + 0.5*IT_0114);
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = IT_0059*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0119 = IT_0068*IT_0117*IT_0118;
    const ccomplex_t IT_0120 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0121 = IT_0009*IT_0120;
    const ccomplex_t IT_0122 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0123 = IT_0009*IT_0122;
    const ccomplex_t IT_0124 = 1.4142135623731*IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*(IT_0121 + 0.5*IT_0124);
    const ccomplex_t IT_0126 = -IT_0125;
    const ccomplex_t IT_0127 = IT_0072*IT_0126;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0129 = IT_0081*IT_0127*IT_0128;
    const ccomplex_t IT_0130 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0131 = IT_0009*IT_0130;
    const ccomplex_t IT_0132 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0133 = IT_0009*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0131 + 0.5*IT_0134);
    const ccomplex_t IT_0136 = -IT_0135;
    const ccomplex_t IT_0137 = IT_0085*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0139 = IT_0094*IT_0137*IT_0138;
    const ccomplex_t IT_0140 = IT_0109 + IT_0119 + IT_0129 + IT_0139;
    const ccomplex_t IT_0141 = conj(IT_0099) + conj(IT_0140);
    const ccomplex_t IT_0142 = s_14*s_34*IT_0037;
    const ccomplex_t IT_0143 = -IT_0142;
    const ccomplex_t IT_0144 = s_13 + IT_0143;
    const ccomplex_t IT_0145 = m_u*IT_0144;
    const ccomplex_t IT_0146 = cos(theta_W);
    const ccomplex_t IT_0147 = cpow(IT_0146, -1);
    const ccomplex_t IT_0148 = (0 + _Complex_I*1.4142135623731)*conj(N_B2)
      *e_em*IT_0147;
    const ccomplex_t IT_0149 = 0.666666666666667*IT_0148;
    const ccomplex_t IT_0150 = IT_0066*IT_0149;
    const ccomplex_t IT_0151 = IT_0068*IT_0118*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*conj(N_B3)
      *e_em*IT_0147;
    const ccomplex_t IT_0153 = 0.666666666666667*IT_0152;
    const ccomplex_t IT_0154 = IT_0079*IT_0153;
    const ccomplex_t IT_0155 = IT_0081*IT_0128*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1.4142135623731)*conj(N_B4)
      *e_em*IT_0147;
    const ccomplex_t IT_0157 = 0.666666666666667*IT_0156;
    const ccomplex_t IT_0158 = IT_0092*IT_0157;
    const ccomplex_t IT_0159 = IT_0094*IT_0138*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1.4142135623731)*conj(N_B1)
      *e_em*IT_0147;
    const ccomplex_t IT_0161 = 0.666666666666667*IT_0160;
    const ccomplex_t IT_0162 = IT_0053*IT_0161;
    const ccomplex_t IT_0163 = IT_0055*IT_0108*IT_0162;
    const ccomplex_t IT_0164 = -IT_0151 + -IT_0155 + -IT_0159 + -IT_0163;
    const ccomplex_t IT_0165 = IT_0116*IT_0149;
    const ccomplex_t IT_0166 = IT_0068*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = m_C_1*IT_0167;
    const ccomplex_t IT_0169 = IT_0126*IT_0153;
    const ccomplex_t IT_0170 = IT_0081*IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0170;
    const ccomplex_t IT_0172 = m_C_1*IT_0171;
    const ccomplex_t IT_0173 = IT_0136*IT_0157;
    const ccomplex_t IT_0174 = IT_0094*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = m_C_1*IT_0175;
    const ccomplex_t IT_0177 = m_u*IT_0017;
    const ccomplex_t IT_0178 = m_u*IT_0026;
    const ccomplex_t IT_0179 = m_u*IT_0035;
    const ccomplex_t IT_0180 = IT_0106*IT_0161;
    const ccomplex_t IT_0181 = IT_0055*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = m_C_1*IT_0182;
    const ccomplex_t IT_0184 = IT_0168 + IT_0172 + IT_0176 + IT_0177 + IT_0178
       + IT_0179 + IT_0183;
    const ccomplex_t IT_0185 = conj(IT_0164) + conj(IT_0184);
    const ccomplex_t IT_0186 = 2*IT_0057 + 2*IT_0070 + 2*IT_0083 + 2*IT_0096;
    const ccomplex_t IT_0187 = s_13*IT_0144;
    const ccomplex_t IT_0188 = 6*IT_0187;
    const ccomplex_t IT_0189 = -IT_0057 + -IT_0070 + -IT_0083 + -IT_0096;
    const ccomplex_t IT_0190 = IT_0017 + IT_0026 + IT_0035;
    const ccomplex_t IT_0191 = conj(IT_0189) + conj(IT_0190);
    const ccomplex_t IT_0192 = s_14*IT_0041;
    const ccomplex_t IT_0193 = s_13*s_34;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = IT_0192 + IT_0194;
    const ccomplex_t IT_0196 = (-6)*IT_0195;
    const ccomplex_t IT_0197 = 2*IT_0167 + 2*IT_0171 + 2*IT_0175 + 2*IT_0182;
    const ccomplex_t IT_0198 = m_u*m_C_1;
    const ccomplex_t IT_0199 = IT_0144*IT_0198;
    const ccomplex_t IT_0200 = (-6)*IT_0199;
    const ccomplex_t IT_0201 = pow(m_W, 2);
    const ccomplex_t IT_0202 = s_13*IT_0201;
    const ccomplex_t IT_0203 = s_14*s_34;
    const ccomplex_t IT_0204 = (-2)*IT_0203;
    const ccomplex_t IT_0205 = IT_0202 + IT_0204;
    const ccomplex_t IT_0206 = IT_0037*IT_0205;
    const ccomplex_t IT_0207 = (-6)*IT_0206;
    const ccomplex_t IT_0208 = 12*s_13;
    const ccomplex_t IT_0209 = IT_0207 + IT_0208;
    const ccomplex_t IT_0210 = IT_0164*IT_0209;
    const ccomplex_t IT_0211 = pow(s_14, 2);
    const ccomplex_t IT_0212 = IT_0037*IT_0211;
    const ccomplex_t IT_0213 = -IT_0001;
    const ccomplex_t IT_0214 = IT_0212 + IT_0213;
    const ccomplex_t IT_0215 = s_13*IT_0214;
    const ccomplex_t IT_0216 = -IT_0167 + -IT_0171 + -IT_0175 + -IT_0182;
    const ccomplex_t IT_0217 = s_13*s_14;
    const ccomplex_t IT_0218 = s_34*IT_0001;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0217 + IT_0219;
    const ccomplex_t IT_0221 = m_C_1*IT_0144;
    const ccomplex_t IT_0222 = IT_0198*IT_0214;
    const ccomplex_t IT_0223 = (-6)*IT_0222;
    const ccomplex_t IT_0224 = m_u*IT_0214;
    const ccomplex_t IT_0225 = (-6)*IT_0224;
    const ccomplex_t IT_0226 = 0.166666666666667*IT_0141;
    const ccomplex_t IT_0227 = (-6)*IT_0221;
    const ccomplex_t IT_0228 = IT_0198*IT_0201;
    const ccomplex_t IT_0229 = 3*IT_0228;
    const ccomplex_t IT_0230 = pow(m_W, 4);
    const ccomplex_t IT_0231 = s_13*IT_0230;
    const ccomplex_t IT_0232 = IT_0037*IT_0231;
    const ccomplex_t IT_0233 = 24*IT_0203;
    const ccomplex_t IT_0234 = (-12)*IT_0202;
    const ccomplex_t IT_0235 = 6*IT_0232 + IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = s_34*IT_0201;
    const ccomplex_t IT_0237 = IT_0037*IT_0236;
    const ccomplex_t IT_0238 = (-0.25)*IT_0237;
    const ccomplex_t IT_0239 = s_34 + IT_0238;
    const ccomplex_t IT_0240 = m_C_1*IT_0239;
    const ccomplex_t IT_0241 = (-4)*IT_0240;
    const ccomplex_t IT_0242 = (-6)*IT_0241;
    const ccomplex_t IT_0243 = (-3)*IT_0228;
    const ccomplex_t IT_0244 = (-6)*IT_0243;
    const ccomplex_t IT_0245 = s_14*IT_0201;
    const ccomplex_t IT_0246 = IT_0037*IT_0245;
    const ccomplex_t IT_0247 = 0.5*IT_0246;
    const ccomplex_t IT_0248 = s_14 + IT_0247;
    const ccomplex_t IT_0249 = m_u*IT_0248;
    const ccomplex_t IT_0250 = 2*IT_0249;
    const ccomplex_t IT_0251 = (-6)*IT_0250;
    const ccomplex_t IT_0252 = conj(IT_0189)*IT_0250;
    const ccomplex_t IT_0253 = conj(IT_0036)*IT_0145;
    const ccomplex_t IT_0254 = conj(IT_0197)*IT_0221 + conj(IT_0186)*IT_0224 +
       IT_0252 + IT_0253;
    const ccomplex_t IT_0255 = (-18)*IT_0198;
    const ccomplex_t IT_0256 = (-0.25)*IT_0246;
    const ccomplex_t IT_0257 = s_14 + IT_0256;
    const ccomplex_t IT_0258 = m_u*IT_0257;
    const ccomplex_t IT_0259 = (-4)*IT_0258;
    const ccomplex_t IT_0260 = (-6)*IT_0259;
    const ccomplex_t IT_0261 = 0.166666666666667*conj(IT_0190);
    const ccomplex_t IT_0262 = conj(IT_0164)*IT_0255;
    const ccomplex_t IT_0263 = conj(IT_0197)*IT_0225;
    const ccomplex_t IT_0264 = conj(IT_0186)*IT_0227;
    const ccomplex_t IT_0265 = 0.5*IT_0237;
    const ccomplex_t IT_0266 = s_34 + IT_0265;
    const ccomplex_t IT_0267 = m_C_1*IT_0266;
    const ccomplex_t IT_0268 = 2*IT_0267;
    const ccomplex_t IT_0269 = (-6)*IT_0268;
    const ccomplex_t IT_0270 = 6*(IT_0099 + IT_0140)*(conj(IT_0036)*IT_0044 +
       IT_0209*IT_0226 + conj(IT_0189)*IT_0241 + 0.166666666666667*conj(IT_0216)
      *IT_0251 + 0.166666666666667*conj(IT_0184)*IT_0255 + 0.166666666666667
      *IT_0262 + 0.166666666666667*IT_0263 + 0.166666666666667*IT_0264 + IT_0261
      *IT_0269);
    const ccomplex_t IT_0271 = (-6)*IT_0232;
    const ccomplex_t IT_0272 = IT_0233 + IT_0271;
    const ccomplex_t IT_0273 = conj(IT_0186)*IT_0220;
    const ccomplex_t IT_0274 = (-6)*IT_0036*(conj(IT_0036)*IT_0043 + -IT_0044
      *IT_0141 + -IT_0145*IT_0185 + (-0.166666666666667)*conj(IT_0186)*IT_0188 +
       (-0.166666666666667)*IT_0191*IT_0196 + (-0.166666666666667)*conj(IT_0197)
      *IT_0200) + IT_0185*IT_0210 + 6*IT_0197*(0.166666666666667*conj(IT_0036)
      *IT_0200 + conj(IT_0197)*IT_0215 + conj(IT_0216)*IT_0220 + IT_0185*IT_0221
       + 0.166666666666667*conj(IT_0186)*IT_0223 + IT_0225*IT_0226) + 6*IT_0186*
      (0.166666666666667*conj(IT_0036)*IT_0188 + conj(IT_0186)*IT_0215 + IT_0191
      *IT_0220 + 0.166666666666667*conj(IT_0197)*IT_0223 + IT_0185*IT_0224 +
       IT_0226*IT_0227) + 6*IT_0216*(conj(IT_0197)*IT_0220 + conj(IT_0190)
      *IT_0229 + 0.166666666666667*conj(IT_0216)*IT_0235 + 0.166666666666667
      *IT_0185*IT_0242 + 0.166666666666667*conj(IT_0189)*IT_0244 + IT_0226
      *IT_0251) + 6*IT_0164*(0.166666666666667*conj(IT_0216)*IT_0242 + IT_0254 +
       IT_0226*IT_0255 + IT_0260*IT_0261) + 6*IT_0184*(0.166666666666667*IT_0185
      *IT_0209 + 0.166666666666667*conj(IT_0216)*IT_0242 + IT_0254 + IT_0226
      *IT_0255 + IT_0260*IT_0261) + IT_0270 + 6*IT_0190*(0.166666666666667*conj
      (IT_0036)*IT_0196 + conj(IT_0216)*IT_0229 + 0.166666666666667*IT_0185
      *IT_0260 + IT_0235*IT_0261 + IT_0226*IT_0269 + 0.166666666666667*conj
      (IT_0189)*IT_0272 + IT_0273) + IT_0189*(conj(IT_0036)*IT_0196 + conj
      (IT_0189)*IT_0235 + 6*IT_0141*IT_0241 + conj(IT_0216)*IT_0244 + 6*IT_0185
      *IT_0250 + 6*IT_0261*IT_0272 + 6*IT_0273);
    return create_ccomplex_return(IT_0274);
}

