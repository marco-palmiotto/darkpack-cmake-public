#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_su_L_to_anti_u_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_su_L_to_anti_u_W(
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
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_su_L = param->m_su_L;
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
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_u, 2);
    const ccomplex_t IT_0001 = pow(m_su_L, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = sin(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = (-0.5)*IT_0008;
    const ccomplex_t IT_0010 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0011 = IT_0007*IT_0010;
    const ccomplex_t IT_0012 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = 1.4142135623731*IT_0013;
    const ccomplex_t IT_0015 = (0 + _Complex_I*1)*(IT_0011 + 0.5*IT_0014);
    const ccomplex_t IT_0016 = -IT_0015;
    const ccomplex_t IT_0017 = IT_0009*IT_0016;
    const ccomplex_t IT_0018 = IT_0002*IT_0017;
    const ccomplex_t IT_0019 = (0 + _Complex_I*1)*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0023 = IT_0007*IT_0022;
    const ccomplex_t IT_0024 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0025 = IT_0007*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = IT_0021*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0031 = IT_0029*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0034 = (-0.5)*IT_0033;
    const ccomplex_t IT_0035 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0036 = IT_0007*IT_0035;
    const ccomplex_t IT_0037 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0038 = IT_0007*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = (0 + _Complex_I*1)*(IT_0036 + 0.5*IT_0039);
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = IT_0034*IT_0041;
    const ccomplex_t IT_0043 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1)*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1.4142135623731)*m_u*conj(N_u4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0047 = (-0.5)*IT_0046;
    const ccomplex_t IT_0048 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0049 = IT_0007*IT_0048;
    const ccomplex_t IT_0050 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + 0.5*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = (0 + _Complex_I*1)*IT_0057;
    const ccomplex_t IT_0059 = -IT_0019 + -IT_0032 + -IT_0045 + -IT_0058;
    const ccomplex_t IT_0060 = cos(theta_W);
    const ccomplex_t IT_0061 = cpow(IT_0060, -1);
    const ccomplex_t IT_0062 = N_B1*e_em;
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = N_W1*e_em;
    const ccomplex_t IT_0066 = IT_0007*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + 3*IT_0067);
    const ccomplex_t IT_0069 = 0.166666666666667*IT_0068;
    const ccomplex_t IT_0070 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0071 = IT_0007*IT_0070;
    const ccomplex_t IT_0072 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0073 = IT_0007*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1)*(IT_0071 + (-0.5)*IT_0074);
    const ccomplex_t IT_0076 = -IT_0075;
    const ccomplex_t IT_0077 = IT_0069*IT_0076;
    const ccomplex_t IT_0078 = IT_0002*IT_0077;
    const ccomplex_t IT_0079 = (0 + _Complex_I*1)*IT_0078;
    const ccomplex_t IT_0080 = N_B3*e_em;
    const ccomplex_t IT_0081 = IT_0061*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = N_W3*e_em;
    const ccomplex_t IT_0084 = IT_0007*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + 3*IT_0085);
    const ccomplex_t IT_0087 = 0.166666666666667*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0089 = IT_0007*IT_0088;
    const ccomplex_t IT_0090 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0091 = IT_0007*IT_0090;
    const ccomplex_t IT_0092 = 1.4142135623731*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0089 + (-0.5)*IT_0092);
    const ccomplex_t IT_0094 = -IT_0093;
    const ccomplex_t IT_0095 = IT_0087*IT_0094;
    const ccomplex_t IT_0096 = IT_0043*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = N_B4*e_em;
    const ccomplex_t IT_0099 = IT_0061*IT_0098;
    const ccomplex_t IT_0100 = 1.4142135623731*IT_0099;
    const ccomplex_t IT_0101 = N_W4*e_em;
    const ccomplex_t IT_0102 = IT_0007*IT_0101;
    const ccomplex_t IT_0103 = 1.4142135623731*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*(IT_0100 + 3*IT_0103);
    const ccomplex_t IT_0105 = 0.166666666666667*IT_0104;
    const ccomplex_t IT_0106 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0107 = IT_0007*IT_0106;
    const ccomplex_t IT_0108 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0109 = IT_0007*IT_0108;
    const ccomplex_t IT_0110 = 1.4142135623731*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*(IT_0107 + (-0.5)*IT_0110);
    const ccomplex_t IT_0112 = -IT_0111;
    const ccomplex_t IT_0113 = IT_0105*IT_0112;
    const ccomplex_t IT_0114 = IT_0056*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = N_B2*e_em;
    const ccomplex_t IT_0117 = IT_0061*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = N_W2*e_em;
    const ccomplex_t IT_0120 = IT_0007*IT_0119;
    const ccomplex_t IT_0121 = 1.4142135623731*IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*(IT_0118 + 3*IT_0121);
    const ccomplex_t IT_0123 = 0.166666666666667*IT_0122;
    const ccomplex_t IT_0124 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0125 = IT_0007*IT_0124;
    const ccomplex_t IT_0126 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0127 = IT_0007*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1)*(IT_0125 + (-0.5)*IT_0128);
    const ccomplex_t IT_0130 = -IT_0129;
    const ccomplex_t IT_0131 = IT_0123*IT_0130;
    const ccomplex_t IT_0132 = IT_0030*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = -IT_0079 + -IT_0097 + -IT_0115 + -IT_0133;
    const ccomplex_t IT_0135 = m_u*m_C_2;
    const ccomplex_t IT_0136 = pow(m_W, 2);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (-3)*IT_0137;
    const ccomplex_t IT_0139 = pow(m_d, 2);
    const ccomplex_t IT_0140 = pow(m_C_2, 2);
    const ccomplex_t IT_0141 = cpow((-2)*s_12 + -IT_0001 + IT_0139 + -IT_0140 
      + -reg_prop, -1);
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*V_ud*e_em
      *IT_0007;
    const ccomplex_t IT_0143 = 0.5*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*V_ud*e_em*V_Wp2*IT_0007;
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = IT_0141*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*IT_0147;
    const ccomplex_t IT_0149 = pow(m_b, 2);
    const ccomplex_t IT_0150 = cpow((-2)*s_12 + -IT_0001 + -IT_0140 + IT_0149 
      + -reg_prop, -1);
    const ccomplex_t IT_0151 = cexp((0 + _Complex_I*1)*delta_wolf);
    const ccomplex_t IT_0152 = (0 + _Complex_I*1.4142135623731)*e_em*IT_0007
      *IT_0151*V_ub_mod;
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = cexp((0 + _Complex_I*-1)*delta_wolf);
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*e_em*V_Wp2*IT_0007*IT_0154
      *V_ub_mod;
    const ccomplex_t IT_0156 = -IT_0155;
    const ccomplex_t IT_0157 = IT_0153*IT_0156;
    const ccomplex_t IT_0158 = IT_0150*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = pow(m_s, 2);
    const ccomplex_t IT_0161 = cpow((-2)*s_12 + -IT_0001 + -IT_0140 + IT_0160 
      + -reg_prop, -1);
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*V_us*e_em*V_Wp2*IT_0007;
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1.4142135623731)*V_us*e_em
      *IT_0007;
    const ccomplex_t IT_0165 = 0.5*IT_0164;
    const ccomplex_t IT_0166 = IT_0163*IT_0165;
    const ccomplex_t IT_0167 = IT_0161*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0148 + IT_0159 + IT_0168;
    const ccomplex_t IT_0170 = 3*IT_0137;
    const ccomplex_t IT_0171 = 2*IT_0019 + 2*IT_0032 + 2*IT_0045 + 2*IT_0058;
    const ccomplex_t IT_0172 = s_13*s_14;
    const ccomplex_t IT_0173 = s_34*IT_0140;
    const ccomplex_t IT_0174 = -IT_0173;
    const ccomplex_t IT_0175 = IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = pow(m_W, -2);
    const ccomplex_t IT_0177 = pow(m_W, 4);
    const ccomplex_t IT_0178 = s_13*IT_0177;
    const ccomplex_t IT_0179 = IT_0176*IT_0178;
    const ccomplex_t IT_0180 = s_14*s_34;
    const ccomplex_t IT_0181 = 24*IT_0180;
    const ccomplex_t IT_0182 = s_13*IT_0136;
    const ccomplex_t IT_0183 = (-12)*IT_0182;
    const ccomplex_t IT_0184 = 6*IT_0179 + IT_0181 + IT_0183;
    const ccomplex_t IT_0185 = s_34*IT_0136;
    const ccomplex_t IT_0186 = IT_0176*IT_0185;
    const ccomplex_t IT_0187 = (-0.25)*IT_0186;
    const ccomplex_t IT_0188 = s_34 + IT_0187;
    const ccomplex_t IT_0189 = m_C_2*IT_0188;
    const ccomplex_t IT_0190 = (-4)*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0192 = IT_0009*IT_0076;
    const ccomplex_t IT_0193 = IT_0002*IT_0191*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0195 = IT_0021*IT_0130;
    const ccomplex_t IT_0196 = IT_0030*IT_0194*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0198 = IT_0034*IT_0094;
    const ccomplex_t IT_0199 = IT_0043*IT_0197*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0201 = IT_0047*IT_0112;
    const ccomplex_t IT_0202 = IT_0056*IT_0200*IT_0201;
    const ccomplex_t IT_0203 = -IT_0193 + -IT_0196 + -IT_0199 + -IT_0202;
    const ccomplex_t IT_0204 = (-6)*IT_0203;
    const ccomplex_t IT_0205 = m_C_2*IT_0032;
    const ccomplex_t IT_0206 = m_C_2*IT_0045;
    const ccomplex_t IT_0207 = m_u*IT_0159;
    const ccomplex_t IT_0208 = m_u*IT_0148;
    const ccomplex_t IT_0209 = m_C_2*IT_0058;
    const ccomplex_t IT_0210 = m_C_2*IT_0019;
    const ccomplex_t IT_0211 = m_u*IT_0168;
    const ccomplex_t IT_0212 = IT_0205 + IT_0206 + IT_0207 + IT_0208 + IT_0209
       + IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = (-6)*IT_0212;
    const ccomplex_t IT_0214 = IT_0204 + IT_0213;
    const ccomplex_t IT_0215 = s_14*IT_0136;
    const ccomplex_t IT_0216 = IT_0176*IT_0215;
    const ccomplex_t IT_0217 = 0.5*IT_0216;
    const ccomplex_t IT_0218 = s_14 + IT_0217;
    const ccomplex_t IT_0219 = m_u*IT_0218;
    const ccomplex_t IT_0220 = 2*IT_0219;
    const ccomplex_t IT_0221 = IT_0079 + IT_0097 + IT_0115 + IT_0133;
    const ccomplex_t IT_0222 = m_C_2*IT_0221;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = (-6)*IT_0223;
    const ccomplex_t IT_0225 = IT_0028*IT_0123;
    const ccomplex_t IT_0226 = IT_0030*IT_0194*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0228 = cos(beta);
    const ccomplex_t IT_0229 = cpow(IT_0228, -1);
    const ccomplex_t IT_0230 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d2)
      *V_ud*e_em*IT_0003*IT_0007*IT_0229;
    const ccomplex_t IT_0231 = (-0.5)*IT_0230;
    const ccomplex_t IT_0232 = IT_0143*IT_0231;
    const ccomplex_t IT_0233 = IT_0141*IT_0227*IT_0232;
    const ccomplex_t IT_0234 = IT_0041*IT_0087;
    const ccomplex_t IT_0235 = IT_0043*IT_0197*IT_0234;
    const ccomplex_t IT_0236 = IT_0054*IT_0105;
    const ccomplex_t IT_0237 = IT_0056*IT_0200*IT_0236;
    const ccomplex_t IT_0238 = IT_0016*IT_0069;
    const ccomplex_t IT_0239 = IT_0002*IT_0191*IT_0238;
    const ccomplex_t IT_0240 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d2)
      *V_us*e_em*IT_0003*IT_0007*IT_0229;
    const ccomplex_t IT_0242 = (-0.5)*IT_0241;
    const ccomplex_t IT_0243 = IT_0165*IT_0242;
    const ccomplex_t IT_0244 = IT_0161*IT_0240*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *e_em*IT_0003*IT_0007*IT_0154*IT_0229*V_ub_mod;
    const ccomplex_t IT_0247 = (-0.5)*IT_0246;
    const ccomplex_t IT_0248 = IT_0153*IT_0247;
    const ccomplex_t IT_0249 = IT_0150*IT_0245*IT_0248;
    const ccomplex_t IT_0250 = IT_0226 + IT_0233 + IT_0235 + IT_0237 + IT_0239
       + IT_0244 + IT_0249;
    const ccomplex_t IT_0251 = (-6)*IT_0250;
    const ccomplex_t IT_0252 = IT_0224 + IT_0251;
    const ccomplex_t IT_0253 = (-6)*conj(IT_0203);
    const ccomplex_t IT_0254 = (-6)*conj(IT_0212);
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = (-6)*conj(IT_0223);
    const ccomplex_t IT_0257 = (-6)*conj(IT_0250);
    const ccomplex_t IT_0258 = IT_0256 + IT_0257;
    const ccomplex_t IT_0259 = 2*IT_0079 + 2*IT_0097 + 2*IT_0115 + 2*IT_0133;
    const ccomplex_t IT_0260 = pow(s_14, 2);
    const ccomplex_t IT_0261 = IT_0176*IT_0260;
    const ccomplex_t IT_0262 = -IT_0140;
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = s_13*IT_0263;
    const ccomplex_t IT_0265 = 6*IT_0259;
    const ccomplex_t IT_0266 = s_14*s_34*IT_0176;
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = s_13 + IT_0267;
    const ccomplex_t IT_0269 = m_C_2*IT_0268;
    const ccomplex_t IT_0270 = 2*IT_0148 + 2*IT_0159 + 2*IT_0168;
    const ccomplex_t IT_0271 = m_u*IT_0268;
    const ccomplex_t IT_0272 = (-2)*IT_0180;
    const ccomplex_t IT_0273 = IT_0182 + IT_0272;
    const ccomplex_t IT_0274 = IT_0176*IT_0273;
    const ccomplex_t IT_0275 = (-6)*IT_0274;
    const ccomplex_t IT_0276 = 12*s_13;
    const ccomplex_t IT_0277 = IT_0275 + IT_0276;
    const ccomplex_t IT_0278 = IT_0203 + IT_0212;
    const ccomplex_t IT_0279 = m_u*IT_0263;
    const ccomplex_t IT_0280 = s_13*IT_0268;
    const ccomplex_t IT_0281 = 6*conj(IT_0259);
    const ccomplex_t IT_0282 = IT_0135*IT_0263;
    const ccomplex_t IT_0283 = IT_0223 + IT_0250;
    const ccomplex_t IT_0284 = (-6)*conj(IT_0259);
    const ccomplex_t IT_0285 = (-18)*conj(IT_0203);
    const ccomplex_t IT_0286 = (-6)*conj(IT_0171);
    const ccomplex_t IT_0287 = pow(s_34, 2);
    const ccomplex_t IT_0288 = IT_0176*IT_0287;
    const ccomplex_t IT_0289 = -IT_0288;
    const ccomplex_t IT_0290 = IT_0000 + IT_0289;
    const ccomplex_t IT_0291 = m_C_2*IT_0290;
    const ccomplex_t IT_0292 = s_13*IT_0290;
    const ccomplex_t IT_0293 = (-6)*IT_0270;
    const ccomplex_t IT_0294 = (-6)*IT_0259;
    const ccomplex_t IT_0295 = IT_0135*IT_0268;
    const ccomplex_t IT_0296 = (-6)*IT_0179;
    const ccomplex_t IT_0297 = IT_0181 + IT_0296;
    const ccomplex_t IT_0298 = (-0.25)*IT_0216;
    const ccomplex_t IT_0299 = s_14 + IT_0298;
    const ccomplex_t IT_0300 = m_u*IT_0299;
    const ccomplex_t IT_0301 = (-4)*IT_0300;
    const ccomplex_t IT_0302 = 0.5*IT_0186;
    const ccomplex_t IT_0303 = s_34 + IT_0302;
    const ccomplex_t IT_0304 = m_C_2*IT_0303;
    const ccomplex_t IT_0305 = 2*IT_0304;
    const ccomplex_t IT_0306 = s_14*IT_0000;
    const ccomplex_t IT_0307 = s_13*s_34;
    const ccomplex_t IT_0308 = -IT_0307;
    const ccomplex_t IT_0309 = IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = 6*conj(IT_0169)*(IT_0059*IT_0170 +
       0.166666666666667*IT_0169*IT_0184 + 0.166666666666667*IT_0175*IT_0265 +
       0.166666666666667*IT_0134*IT_0297 + 0.166666666666667*IT_0214*IT_0301 +
       0.166666666666667*IT_0252*IT_0305 + 0.166666666666667*IT_0293*IT_0309);
    const ccomplex_t IT_0311 = (-6)*conj(IT_0134)*(IT_0059*IT_0138 + (
      -0.166666666666667)*IT_0134*IT_0184 + (-0.166666666666667)*IT_0175*IT_0265
       + -IT_0220*IT_0278 + -IT_0190*IT_0283 + (-0.166666666666667)*IT_0169
      *IT_0297 + (-0.166666666666667)*IT_0293*IT_0309);
    const ccomplex_t IT_0312 = (-6)*conj(IT_0270);
    const ccomplex_t IT_0313 = 6*IT_0134;
    const ccomplex_t IT_0314 = (-6)*conj(IT_0059)*(IT_0134*IT_0138 + -IT_0169
      *IT_0170 + -IT_0171*IT_0175 + (-0.166666666666667)*IT_0059*IT_0184 + (
      -0.166666666666667)*IT_0190*IT_0214 + (-0.166666666666667)*IT_0220*IT_0252
      ) + IT_0059*(IT_0190*IT_0255 + IT_0220*IT_0258) + conj(IT_0259)*IT_0264
      *IT_0265 + 6*conj(IT_0203)*(IT_0134*IT_0220 + IT_0171*IT_0269 + IT_0270
      *IT_0271 + 0.166666666666667*IT_0277*IT_0278 + 0.166666666666667*IT_0265
      *IT_0279) + (IT_0134*IT_0175 + IT_0278*IT_0279 + IT_0270*IT_0280)*IT_0281 
      + (IT_0171*IT_0282 + IT_0269*IT_0283)*IT_0284 + IT_0283*(IT_0135*IT_0285 +
       IT_0279*IT_0286) + 6*conj(IT_0270)*(IT_0271*IT_0278 + 0.166666666666667
      *IT_0265*IT_0280 + IT_0283*IT_0291 + 0.166666666666667*IT_0292*IT_0293) +
       6*conj(IT_0171)*(IT_0059*IT_0175 + IT_0171*IT_0264 + IT_0269*IT_0278 +
       0.166666666666667*IT_0282*IT_0294 + 0.166666666666667*IT_0293*IT_0295) +
       IT_0310 + IT_0311 + (IT_0171*IT_0295 + IT_0134*IT_0309)*IT_0312 + IT_0169
      *(IT_0175*IT_0281 + IT_0255*IT_0301 + IT_0258*IT_0305 + IT_0309*IT_0312) +
       (conj(IT_0223) + conj(IT_0250))*(IT_0135*((-18)*IT_0203 + (-18)*IT_0212) 
      + (-6)*IT_0171*IT_0279 + IT_0277*IT_0283 + 6*IT_0270*IT_0291 + IT_0269
      *IT_0294 + IT_0190*IT_0313) + conj(IT_0212)*(IT_0135*((-18)*IT_0223 + (-18
      )*IT_0250) + 6*IT_0171*IT_0269 + 6*IT_0270*IT_0271 + IT_0277*IT_0278 +
       IT_0265*IT_0279 + IT_0220*IT_0313);
    return create_ccomplex_return(IT_0314);
}

