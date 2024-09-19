#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_N_3_to_h_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_N_3_to_h_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_Hp = param->m_Hp;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_W = param->Gamma_W;
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
    const ccomplex_t IT_0000 = cos(alpha);
    const ccomplex_t IT_0001 = cos(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0005 + IT_0009 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0016;
    const ccomplex_t IT_0021 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0014 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = m_C_1*IT_0030;
    const ccomplex_t IT_0032 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0033 = IT_0016*IT_0032;
    const ccomplex_t IT_0034 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0035 = IT_0016*IT_0034;
    const ccomplex_t IT_0036 = 1.4142135623731*IT_0035;
    const ccomplex_t IT_0037 = (0 + _Complex_I*1)*(IT_0033 + 0.5*IT_0036);
    const ccomplex_t IT_0038 = -IT_0037;
    const ccomplex_t IT_0039 = pow(m_h, 2);
    const ccomplex_t IT_0040 = pow(m_N_3, 2);
    const ccomplex_t IT_0041 = cpow((-2)*s_23 + IT_0039 + IT_0040 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0042 = IT_0038*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = IT_0031*IT_0043;
    const ccomplex_t IT_0045 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0046 = IT_0003*IT_0045;
    const ccomplex_t IT_0047 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0048 = IT_0007*IT_0047;
    const ccomplex_t IT_0049 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0050 = IT_0003*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0052 = IT_0007*IT_0051;
    const ccomplex_t IT_0053 = IT_0046 + IT_0048 + IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0055 = IT_0017*IT_0054;
    const ccomplex_t IT_0056 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0057 = IT_0020*IT_0056;
    const ccomplex_t IT_0058 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0059 = IT_0017*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0061 = IT_0020*IT_0060;
    const ccomplex_t IT_0062 = -IT_0055 + -IT_0057 + -IT_0059 + -IT_0061;
    const ccomplex_t IT_0063 = IT_0053 + IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = (-0.5)*IT_0064;
    const ccomplex_t IT_0066 = m_C_1*IT_0065;
    const ccomplex_t IT_0067 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0068 = IT_0016*IT_0067;
    const ccomplex_t IT_0069 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0070 = IT_0016*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + 0.5*IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = cpow((-2)*s_23 + IT_0039 + IT_0040 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0075 = IT_0073*IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1)*IT_0075;
    const ccomplex_t IT_0077 = IT_0066*IT_0076;
    const ccomplex_t IT_0078 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0079 = IT_0003*IT_0078;
    const ccomplex_t IT_0080 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0081 = IT_0007*IT_0080;
    const ccomplex_t IT_0082 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0083 = IT_0003*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0085 = IT_0007*IT_0084;
    const ccomplex_t IT_0086 = IT_0079 + IT_0081 + IT_0083 + IT_0085;
    const ccomplex_t IT_0087 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0088 = IT_0017*IT_0087;
    const ccomplex_t IT_0089 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0090 = IT_0020*IT_0089;
    const ccomplex_t IT_0091 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0092 = IT_0017*IT_0091;
    const ccomplex_t IT_0093 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0094 = IT_0020*IT_0093;
    const ccomplex_t IT_0095 = -IT_0088 + -IT_0090 + -IT_0092 + -IT_0094;
    const ccomplex_t IT_0096 = IT_0086 + IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (-0.5)*IT_0097;
    const ccomplex_t IT_0099 = m_C_1*IT_0098;
    const ccomplex_t IT_0100 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0101 = IT_0016*IT_0100;
    const ccomplex_t IT_0102 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0103 = IT_0016*IT_0102;
    const ccomplex_t IT_0104 = 1.4142135623731*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0101 + 0.5*IT_0104);
    const ccomplex_t IT_0106 = -IT_0105;
    const ccomplex_t IT_0107 = cpow((-2)*s_23 + IT_0039 + IT_0040 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0108 = IT_0106*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = IT_0099*IT_0109;
    const ccomplex_t IT_0111 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0112 = IT_0003*IT_0111;
    const ccomplex_t IT_0113 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0114 = IT_0007*IT_0113;
    const ccomplex_t IT_0115 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0116 = IT_0017*IT_0115;
    const ccomplex_t IT_0117 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0118 = IT_0020*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0112 + IT_0114 + 
      -IT_0116 + -IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = m_C_1*IT_0120;
    const ccomplex_t IT_0122 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0123 = IT_0016*IT_0122;
    const ccomplex_t IT_0124 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0125 = IT_0016*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + 0.5*IT_0126);
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = cpow((-2)*s_23 + IT_0039 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0130 = IT_0128*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = IT_0121*IT_0131;
    const ccomplex_t IT_0133 = cpow((-2)*s_13 + IT_0039 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0134 = V_u1*U_Wm1;
    const ccomplex_t IT_0135 = IT_0000*IT_0134;
    const ccomplex_t IT_0136 = U_d1*V_Wp1;
    const ccomplex_t IT_0137 = IT_0006*IT_0136;
    const ccomplex_t IT_0138 = IT_0016*(IT_0135 + -IT_0137);
    const ccomplex_t IT_0139 = 1.4142135623731*e_em*IT_0138;
    const ccomplex_t IT_0140 = 0.5*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = IT_0128*IT_0141;
    const ccomplex_t IT_0143 = IT_0133*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = m_N_3*IT_0144;
    const ccomplex_t IT_0146 = V_u1*U_Wm2;
    const ccomplex_t IT_0147 = IT_0000*IT_0146;
    const ccomplex_t IT_0148 = U_d2*V_Wp1;
    const ccomplex_t IT_0149 = IT_0006*IT_0148;
    const ccomplex_t IT_0150 = IT_0016*(IT_0147 + -IT_0149);
    const ccomplex_t IT_0151 = 1.4142135623731*e_em*IT_0150;
    const ccomplex_t IT_0152 = 0.5*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0154 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0155 = IT_0016*IT_0154;
    const ccomplex_t IT_0156 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0157 = IT_0016*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + 0.5*IT_0158);
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = IT_0153*IT_0160;
    const ccomplex_t IT_0162 = pow(m_C_1, 2);
    const ccomplex_t IT_0163 = cpow((-2)*s_13 + IT_0039 + IT_0162 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0164 = IT_0161*IT_0163;
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = m_N_3*IT_0165;
    const ccomplex_t IT_0167 = IT_0044 + IT_0077 + IT_0110 + IT_0132 + 
      -IT_0145 + -IT_0166;
    const ccomplex_t IT_0168 = sin(beta);
    const ccomplex_t IT_0169 = cos(beta);
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*m_W*e_em*IT_0016*(IT_0000
      *IT_0168 + -IT_0006*IT_0169);
    const ccomplex_t IT_0171 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0172 = IT_0016*IT_0171;
    const ccomplex_t IT_0173 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0174 = IT_0016*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0172 + (-0.5)*IT_0175);
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = IT_0170*IT_0177;
    const ccomplex_t IT_0179 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0040 + -IT_0162 + -reg_prop, -1);
    const ccomplex_t IT_0180 = IT_0178*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = m_N_1*IT_0065;
    const ccomplex_t IT_0183 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0184 = IT_0016*IT_0183;
    const ccomplex_t IT_0185 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0186 = IT_0016*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0184 + (-0.5)*IT_0187);
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0074*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = IT_0182*IT_0191;
    const ccomplex_t IT_0193 = m_N_2*IT_0030;
    const ccomplex_t IT_0194 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0195 = IT_0016*IT_0194;
    const ccomplex_t IT_0196 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0197 = IT_0016*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = (0 + _Complex_I*1)*(IT_0195 + (-0.5)*IT_0198);
    const ccomplex_t IT_0200 = -IT_0199;
    const ccomplex_t IT_0201 = IT_0041*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0193*IT_0202;
    const ccomplex_t IT_0204 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0205 = IT_0141*IT_0177;
    const ccomplex_t IT_0206 = IT_0133*IT_0204*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0208 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0209 = IT_0016*IT_0208;
    const ccomplex_t IT_0210 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0211 = IT_0016*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*(IT_0209 + (-0.5)*IT_0212);
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0153*IT_0214;
    const ccomplex_t IT_0216 = IT_0163*IT_0207*IT_0215;
    const ccomplex_t IT_0217 = m_N_4*IT_0098;
    const ccomplex_t IT_0218 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0219 = IT_0016*IT_0218;
    const ccomplex_t IT_0220 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0221 = IT_0016*IT_0220;
    const ccomplex_t IT_0222 = 1.4142135623731*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*(IT_0219 + (-0.5)*IT_0222);
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = IT_0107*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0217*IT_0226;
    const ccomplex_t IT_0228 = m_N_3*IT_0120;
    const ccomplex_t IT_0229 = IT_0129*IT_0177;
    const ccomplex_t IT_0230 = (0 + _Complex_I*1)*IT_0229;
    const ccomplex_t IT_0231 = IT_0228*IT_0230;
    const ccomplex_t IT_0232 = -IT_0181 + -IT_0192 + -IT_0203 + IT_0206 +
       IT_0216 + -IT_0227 + -IT_0231;
    const ccomplex_t IT_0233 = pow(m_W, -2);
    const ccomplex_t IT_0234 = pow(m_W, 2);
    const ccomplex_t IT_0235 = s_12*IT_0234;
    const ccomplex_t IT_0236 = s_14*s_24;
    const ccomplex_t IT_0237 = (-2)*IT_0236;
    const ccomplex_t IT_0238 = IT_0235 + IT_0237;
    const ccomplex_t IT_0239 = IT_0233*IT_0238;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = s_12 + IT_0240;
    const ccomplex_t IT_0242 = 4*IT_0241;
    const ccomplex_t IT_0243 = IT_0167*IT_0242;
    const ccomplex_t IT_0244 = (0 + _Complex_I*1)*e_em*IT_0016*(IT_0000
      *IT_0168 + -IT_0006*IT_0169);
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = IT_0016*IT_0169;
    const ccomplex_t IT_0247 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0248 = IT_0246*IT_0247;
    const ccomplex_t IT_0249 = (-2)*IT_0248;
    const ccomplex_t IT_0250 = IT_0002*IT_0169;
    const ccomplex_t IT_0251 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0252 = IT_0250*IT_0251;
    const ccomplex_t IT_0253 = 1.4142135623731*IT_0252;
    const ccomplex_t IT_0254 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0255 = IT_0246*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = IT_0253 + IT_0256;
    const ccomplex_t IT_0258 = IT_0249 + IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = (-0.5)*IT_0259;
    const ccomplex_t IT_0261 = cpow((-2)*s_12 + -IT_0040 + -IT_0162 + IT_0234 
      + -reg_prop, -1);
    const ccomplex_t IT_0262 = IT_0260*IT_0261;
    const ccomplex_t IT_0263 = (0 + _Complex_I*1)*IT_0262;
    const ccomplex_t IT_0264 = IT_0245*IT_0263;
    const ccomplex_t IT_0265 = (0 + _Complex_I*1)*e_em*IT_0016*(IT_0006
      *IT_0168 + IT_0000*IT_0169);
    const ccomplex_t IT_0266 = -IT_0265;
    const ccomplex_t IT_0267 = IT_0016*IT_0168;
    const ccomplex_t IT_0268 = IT_0254*IT_0267;
    const ccomplex_t IT_0269 = 1.4142135623731*IT_0268;
    const ccomplex_t IT_0270 = IT_0002*IT_0168;
    const ccomplex_t IT_0271 = IT_0251*IT_0270;
    const ccomplex_t IT_0272 = 1.4142135623731*IT_0271;
    const ccomplex_t IT_0273 = IT_0269 + IT_0272;
    const ccomplex_t IT_0274 = IT_0247*IT_0267;
    const ccomplex_t IT_0275 = (-2)*IT_0274;
    const ccomplex_t IT_0276 = IT_0273 + IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = 0.5*IT_0277;
    const ccomplex_t IT_0279 = cpow((-2)*s_12 + -IT_0040 + -IT_0162 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0280 = IT_0278*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = IT_0266*IT_0281;
    const ccomplex_t IT_0283 = -IT_0264 + -IT_0282;
    const ccomplex_t IT_0284 = s_14*s_24*IT_0233;
    const ccomplex_t IT_0285 = -IT_0284;
    const ccomplex_t IT_0286 = s_12 + IT_0285;
    const ccomplex_t IT_0287 = m_C_1*IT_0286;
    const ccomplex_t IT_0288 = (-2)*IT_0287;
    const ccomplex_t IT_0289 = IT_0283*IT_0288;
    const ccomplex_t IT_0290 = IT_0243 + IT_0289;
    const ccomplex_t IT_0291 = m_C_1*m_N_3;
    const ccomplex_t IT_0292 = IT_0286*IT_0291;
    const ccomplex_t IT_0293 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0294 = IT_0267*IT_0293;
    const ccomplex_t IT_0295 = 2*IT_0294;
    const ccomplex_t IT_0296 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0297 = IT_0270*IT_0296;
    const ccomplex_t IT_0298 = 1.4142135623731*IT_0297;
    const ccomplex_t IT_0299 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0300 = IT_0267*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = IT_0298 + IT_0301;
    const ccomplex_t IT_0303 = IT_0295 + IT_0302;
    const ccomplex_t IT_0304 = (0 + _Complex_I*1)*IT_0303;
    const ccomplex_t IT_0305 = 0.5*IT_0304;
    const ccomplex_t IT_0306 = IT_0261*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*IT_0306;
    const ccomplex_t IT_0308 = IT_0245*IT_0307;
    const ccomplex_t IT_0309 = IT_0246*IT_0293;
    const ccomplex_t IT_0310 = 2*IT_0309;
    const ccomplex_t IT_0311 = IT_0250*IT_0296;
    const ccomplex_t IT_0312 = 1.4142135623731*IT_0311;
    const ccomplex_t IT_0313 = IT_0246*IT_0299;
    const ccomplex_t IT_0314 = 1.4142135623731*IT_0313;
    const ccomplex_t IT_0315 = IT_0312 + IT_0314;
    const ccomplex_t IT_0316 = IT_0310 + IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = 0.5*IT_0317;
    const ccomplex_t IT_0319 = IT_0279*IT_0318;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = IT_0266*IT_0320;
    const ccomplex_t IT_0322 = IT_0308 + IT_0321;
    const ccomplex_t IT_0323 = pow(s_14, 2);
    const ccomplex_t IT_0324 = IT_0233*IT_0323;
    const ccomplex_t IT_0325 = -IT_0162;
    const ccomplex_t IT_0326 = IT_0324 + IT_0325;
    const ccomplex_t IT_0327 = s_12*IT_0326;
    const ccomplex_t IT_0328 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0329 = IT_0000*IT_0328;
    const ccomplex_t IT_0330 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0331 = IT_0006*IT_0330;
    const ccomplex_t IT_0332 = IT_0016*(IT_0329 + -IT_0331);
    const ccomplex_t IT_0333 = 1.4142135623731*e_em*IT_0332;
    const ccomplex_t IT_0334 = 0.5*IT_0333;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0177*IT_0335;
    const ccomplex_t IT_0337 = IT_0133*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0340 = IT_0000*IT_0339;
    const ccomplex_t IT_0341 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0342 = IT_0006*IT_0341;
    const ccomplex_t IT_0343 = IT_0016*(IT_0340 + -IT_0342);
    const ccomplex_t IT_0344 = 1.4142135623731*e_em*IT_0343;
    const ccomplex_t IT_0345 = 0.5*IT_0344;
    const ccomplex_t IT_0346 = (0 + _Complex_I*1)*IT_0345;
    const ccomplex_t IT_0347 = IT_0214*IT_0346;
    const ccomplex_t IT_0348 = IT_0163*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = IT_0338 + IT_0349;
    const ccomplex_t IT_0351 = IT_0030*IT_0043;
    const ccomplex_t IT_0352 = IT_0098*IT_0109;
    const ccomplex_t IT_0353 = IT_0120*IT_0131;
    const ccomplex_t IT_0354 = IT_0065*IT_0076;
    const ccomplex_t IT_0355 = -IT_0351 + -IT_0352 + -IT_0353 + -IT_0354;
    const ccomplex_t IT_0356 = conj(IT_0350) + conj(IT_0355);
    const ccomplex_t IT_0357 = s_12*s_14;
    const ccomplex_t IT_0358 = s_24*IT_0162;
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = IT_0357 + IT_0359;
    const ccomplex_t IT_0361 = (-2)*IT_0360;
    const ccomplex_t IT_0362 = 2*IT_0351 + 2*IT_0352 + 2*IT_0353 + 2*IT_0354;
    const ccomplex_t IT_0363 = (-2)*IT_0327;
    const ccomplex_t IT_0364 = -IT_0264 + -IT_0282 + (-2)*IT_0338 + (-2)
      *IT_0349;
    const ccomplex_t IT_0365 = s_12*IT_0286;
    const ccomplex_t IT_0366 = (-2)*IT_0365;
    const ccomplex_t IT_0367 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0368 = IT_0003*IT_0367;
    const ccomplex_t IT_0369 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0370 = IT_0007*IT_0369;
    const ccomplex_t IT_0371 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0372 = IT_0003*IT_0371;
    const ccomplex_t IT_0373 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0374 = IT_0007*IT_0373;
    const ccomplex_t IT_0375 = IT_0368 + IT_0370 + IT_0372 + IT_0374;
    const ccomplex_t IT_0376 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0377 = IT_0017*IT_0376;
    const ccomplex_t IT_0378 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0379 = IT_0020*IT_0378;
    const ccomplex_t IT_0380 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0381 = IT_0017*IT_0380;
    const ccomplex_t IT_0382 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0383 = IT_0020*IT_0382;
    const ccomplex_t IT_0384 = -IT_0377 + -IT_0379 + -IT_0381 + -IT_0383;
    const ccomplex_t IT_0385 = IT_0375 + IT_0384;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*IT_0385;
    const ccomplex_t IT_0387 = (-0.5)*IT_0386;
    const ccomplex_t IT_0388 = IT_0191*IT_0387;
    const ccomplex_t IT_0389 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0390 = IT_0003*IT_0389;
    const ccomplex_t IT_0391 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0392 = IT_0007*IT_0391;
    const ccomplex_t IT_0393 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0394 = IT_0003*IT_0393;
    const ccomplex_t IT_0395 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0396 = IT_0007*IT_0395;
    const ccomplex_t IT_0397 = IT_0390 + IT_0392 + IT_0394 + IT_0396;
    const ccomplex_t IT_0398 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0399 = IT_0017*IT_0398;
    const ccomplex_t IT_0400 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0401 = IT_0020*IT_0400;
    const ccomplex_t IT_0402 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0403 = IT_0017*IT_0402;
    const ccomplex_t IT_0404 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0405 = IT_0020*IT_0404;
    const ccomplex_t IT_0406 = -IT_0399 + -IT_0401 + -IT_0403 + -IT_0405;
    const ccomplex_t IT_0407 = IT_0397 + IT_0406;
    const ccomplex_t IT_0408 = (0 + _Complex_I*1)*IT_0407;
    const ccomplex_t IT_0409 = (-0.5)*IT_0408;
    const ccomplex_t IT_0410 = IT_0202*IT_0409;
    const ccomplex_t IT_0411 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0412 = IT_0003*IT_0411;
    const ccomplex_t IT_0413 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0414 = IT_0007*IT_0413;
    const ccomplex_t IT_0415 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0416 = IT_0017*IT_0415;
    const ccomplex_t IT_0417 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0418 = IT_0020*IT_0417;
    const ccomplex_t IT_0419 = (0 + _Complex_I*1)*(IT_0412 + IT_0414 + 
      -IT_0416 + -IT_0418);
    const ccomplex_t IT_0420 = -IT_0419;
    const ccomplex_t IT_0421 = IT_0230*IT_0420;
    const ccomplex_t IT_0422 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0423 = IT_0003*IT_0422;
    const ccomplex_t IT_0424 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0425 = IT_0007*IT_0424;
    const ccomplex_t IT_0426 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0427 = IT_0003*IT_0426;
    const ccomplex_t IT_0428 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0429 = IT_0007*IT_0428;
    const ccomplex_t IT_0430 = IT_0423 + IT_0425 + IT_0427 + IT_0429;
    const ccomplex_t IT_0431 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0432 = IT_0017*IT_0431;
    const ccomplex_t IT_0433 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0434 = IT_0020*IT_0433;
    const ccomplex_t IT_0435 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0436 = IT_0017*IT_0435;
    const ccomplex_t IT_0437 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0438 = IT_0020*IT_0437;
    const ccomplex_t IT_0439 = -IT_0432 + -IT_0434 + -IT_0436 + -IT_0438;
    const ccomplex_t IT_0440 = IT_0430 + IT_0439;
    const ccomplex_t IT_0441 = (0 + _Complex_I*1)*IT_0440;
    const ccomplex_t IT_0442 = (-0.5)*IT_0441;
    const ccomplex_t IT_0443 = IT_0226*IT_0442;
    const ccomplex_t IT_0444 = IT_0308 + IT_0321 + (-2)*IT_0388 + (-2)*IT_0410
       + (-2)*IT_0421 + (-2)*IT_0443;
    const ccomplex_t IT_0445 = IT_0291*IT_0326;
    const ccomplex_t IT_0446 = (-2)*IT_0445;
    const ccomplex_t IT_0447 = 2*IT_0144 + 2*IT_0165;
    const ccomplex_t IT_0448 = (-2)*IT_0292;
    const ccomplex_t IT_0449 = IT_0128*IT_0170;
    const ccomplex_t IT_0450 = IT_0179*IT_0449;
    const ccomplex_t IT_0451 = (0 + _Complex_I*1)*IT_0450;
    const ccomplex_t IT_0452 = IT_0160*IT_0346;
    const ccomplex_t IT_0453 = IT_0163*IT_0207*IT_0452;
    const ccomplex_t IT_0454 = IT_0128*IT_0335;
    const ccomplex_t IT_0455 = IT_0133*IT_0204*IT_0454;
    const ccomplex_t IT_0456 = m_N_1*IT_0387;
    const ccomplex_t IT_0457 = IT_0076*IT_0456;
    const ccomplex_t IT_0458 = m_N_2*IT_0409;
    const ccomplex_t IT_0459 = IT_0043*IT_0458;
    const ccomplex_t IT_0460 = m_N_3*IT_0420;
    const ccomplex_t IT_0461 = IT_0131*IT_0460;
    const ccomplex_t IT_0462 = -IT_0451 + IT_0453 + IT_0455 + -IT_0457 + 
      -IT_0459 + (-0.5)*IT_0461;
    const ccomplex_t IT_0463 = m_C_1*IT_0409;
    const ccomplex_t IT_0464 = IT_0202*IT_0463;
    const ccomplex_t IT_0465 = m_C_1*IT_0420;
    const ccomplex_t IT_0466 = IT_0230*IT_0465;
    const ccomplex_t IT_0467 = m_N_3*IT_0338;
    const ccomplex_t IT_0468 = m_N_4*IT_0442;
    const ccomplex_t IT_0469 = IT_0109*IT_0468;
    const ccomplex_t IT_0470 = m_C_1*IT_0442;
    const ccomplex_t IT_0471 = IT_0226*IT_0470;
    const ccomplex_t IT_0472 = m_C_1*IT_0387;
    const ccomplex_t IT_0473 = IT_0191*IT_0472;
    const ccomplex_t IT_0474 = m_N_3*IT_0349;
    const ccomplex_t IT_0475 = (-0.5)*IT_0461 + IT_0464 + IT_0466 + -IT_0467 +
       -IT_0469 + IT_0471 + IT_0473 + -IT_0474;
    const ccomplex_t IT_0476 = conj(IT_0462) + conj(IT_0475);
    const ccomplex_t IT_0477 = m_N_3*IT_0326;
    const ccomplex_t IT_0478 = (-2)*IT_0477;
    const ccomplex_t IT_0479 = IT_0292*conj(IT_0322) + conj(IT_0283)*IT_0327 +
       0.5*IT_0356*IT_0361 + 0.5*conj(IT_0362)*IT_0363 + 0.5*conj(IT_0364)
      *IT_0366 + 0.5*conj(IT_0444)*IT_0446 + 0.5*conj(IT_0447)*IT_0448 + 0.5
      *IT_0476*IT_0478;
    const ccomplex_t IT_0480 = pow(s_24, 2);
    const ccomplex_t IT_0481 = IT_0233*IT_0480;
    const ccomplex_t IT_0482 = -IT_0040;
    const ccomplex_t IT_0483 = IT_0481 + IT_0482;
    const ccomplex_t IT_0484 = m_C_1*IT_0483;
    const ccomplex_t IT_0485 = m_N_3*IT_0286;
    const ccomplex_t IT_0486 = conj(IT_0167) + conj(IT_0232);
    const ccomplex_t IT_0487 = -IT_0144 + -IT_0165;
    const ccomplex_t IT_0488 = IT_0388 + IT_0410 + IT_0421 + IT_0443;
    const ccomplex_t IT_0489 = conj(IT_0487) + conj(IT_0488);
    const ccomplex_t IT_0490 = s_12*s_24;
    const ccomplex_t IT_0491 = s_14*IT_0040;
    const ccomplex_t IT_0492 = -IT_0491;
    const ccomplex_t IT_0493 = IT_0490 + IT_0492;
    const ccomplex_t IT_0494 = (-2)*IT_0493;
    const ccomplex_t IT_0495 = IT_0291*IT_0483;
    const ccomplex_t IT_0496 = (-2)*IT_0495;
    const ccomplex_t IT_0497 = s_12*IT_0483;
    const ccomplex_t IT_0498 = -conj(IT_0322);
    const ccomplex_t IT_0499 = conj(IT_0283)*IT_0292 + 0.5*IT_0366*conj
      (IT_0444) + 0.5*conj(IT_0362)*IT_0448 + IT_0476*IT_0484 + IT_0485*IT_0486 
      + 0.5*IT_0489*IT_0494 + 0.5*conj(IT_0364)*IT_0496 + -IT_0497*(conj(IT_0447
      ) + IT_0498);
    const ccomplex_t IT_0500 = 2*IT_0365;
    const ccomplex_t IT_0501 = 0.5*conj(IT_0322);
    const ccomplex_t IT_0502 = 0.5*conj(IT_0283);
    const ccomplex_t IT_0503 = 2*IT_0362*(IT_0356*IT_0360 + IT_0327*conj
      (IT_0362) + conj(IT_0444)*IT_0445 + IT_0292*conj(IT_0447) + IT_0476
      *IT_0477 + IT_0287*IT_0486 + 0.5*conj(IT_0364)*IT_0500 + IT_0448*IT_0501 +
       IT_0363*IT_0502);
    const ccomplex_t IT_0504 = 2*IT_0364*(IT_0292*conj(IT_0444) + IT_0476
      *IT_0485 + IT_0484*IT_0486 + IT_0356*IT_0493 + conj(IT_0447)*IT_0495 +
       conj(IT_0364)*IT_0497 + 0.5*conj(IT_0362)*IT_0500 + IT_0496*IT_0501 +
       IT_0366*IT_0502);
    const ccomplex_t IT_0505 = 0.5*IT_0486;
    const ccomplex_t IT_0506 = 2*IT_0444*(IT_0292*conj(IT_0364) + IT_0327*conj
      (IT_0444) + conj(IT_0362)*IT_0445 + 0.5*IT_0288*IT_0476 + IT_0360*IT_0489 
      + 0.5*conj(IT_0447)*IT_0500 + IT_0366*IT_0501 + IT_0446*IT_0502 + IT_0478
      *IT_0505);
    const ccomplex_t IT_0507 = (-2)*IT_0484;
    const ccomplex_t IT_0508 = (-2)*IT_0485;
    const ccomplex_t IT_0509 = 2*IT_0447*(IT_0292*conj(IT_0362) + IT_0489
      *IT_0493 + conj(IT_0364)*IT_0495 + IT_0497*(conj(IT_0447) + IT_0498) + 0.5
      *conj(IT_0444)*IT_0500 + IT_0448*IT_0502 + 0.5*IT_0476*IT_0507 + IT_0505
      *IT_0508);
    const ccomplex_t IT_0510 = 6*IT_0291;
    const ccomplex_t IT_0511 = 0.5*conj(IT_0475)*IT_0510;
    const ccomplex_t IT_0512 = s_24*IT_0234;
    const ccomplex_t IT_0513 = IT_0233*IT_0512;
    const ccomplex_t IT_0514 = (-0.25)*IT_0513;
    const ccomplex_t IT_0515 = s_24 + IT_0514;
    const ccomplex_t IT_0516 = m_C_1*IT_0515;
    const ccomplex_t IT_0517 = (-4)*IT_0516;
    const ccomplex_t IT_0518 = (-2)*IT_0517;
    const ccomplex_t IT_0519 = 0.5*conj(IT_0355);
    const ccomplex_t IT_0520 = IT_0518*IT_0519;
    const ccomplex_t IT_0521 = s_14*IT_0234;
    const ccomplex_t IT_0522 = IT_0233*IT_0521;
    const ccomplex_t IT_0523 = 0.5*IT_0522;
    const ccomplex_t IT_0524 = s_14 + IT_0523;
    const ccomplex_t IT_0525 = m_N_3*IT_0524;
    const ccomplex_t IT_0526 = 2*IT_0525;
    const ccomplex_t IT_0527 = (-2)*IT_0526;
    const ccomplex_t IT_0528 = 0.5*conj(IT_0488)*IT_0527;
    const ccomplex_t IT_0529 = conj(IT_0462)*IT_0510;
    const ccomplex_t IT_0530 = 0.5*IT_0529;
    const ccomplex_t IT_0531 = conj(IT_0283)*IT_0288;
    const ccomplex_t IT_0532 = 0.5*IT_0531;
    const ccomplex_t IT_0533 = conj(IT_0444)*IT_0478;
    const ccomplex_t IT_0534 = 0.5*IT_0533;
    const ccomplex_t IT_0535 = conj(IT_0447)*IT_0508;
    const ccomplex_t IT_0536 = 0.5*IT_0535;
    const ccomplex_t IT_0537 = conj(IT_0364)*IT_0484;
    const ccomplex_t IT_0538 = 0.5*IT_0513;
    const ccomplex_t IT_0539 = s_24 + IT_0538;
    const ccomplex_t IT_0540 = m_C_1*IT_0539;
    const ccomplex_t IT_0541 = 2*IT_0540;
    const ccomplex_t IT_0542 = conj(IT_0350)*IT_0541;
    const ccomplex_t IT_0543 = (-0.25)*IT_0522;
    const ccomplex_t IT_0544 = s_14 + IT_0543;
    const ccomplex_t IT_0545 = m_N_3*IT_0544;
    const ccomplex_t IT_0546 = (-4)*IT_0545;
    const ccomplex_t IT_0547 = conj(IT_0487)*IT_0546;
    const ccomplex_t IT_0548 = IT_0288*conj(IT_0444);
    const ccomplex_t IT_0549 = conj(IT_0447)*IT_0507;
    const ccomplex_t IT_0550 = (-2)*IT_0541;
    const ccomplex_t IT_0551 = (-2)*IT_0546;
    const ccomplex_t IT_0552 = conj(IT_0355)*IT_0526;
    const ccomplex_t IT_0553 = conj(IT_0488)*IT_0517;
    const ccomplex_t IT_0554 = 0.5*IT_0242*IT_0476 + conj(IT_0362)*IT_0477 +
       conj(IT_0322)*IT_0484 + conj(IT_0364)*IT_0485 + IT_0478*IT_0502 + IT_0505
      *IT_0510 + 0.5*IT_0548 + 0.5*IT_0549 + 0.5*conj(IT_0487)*IT_0550 + 0.5
      *conj(IT_0350)*IT_0551 + IT_0552 + IT_0553;
    const ccomplex_t IT_0555 = IT_0234*IT_0291;
    const ccomplex_t IT_0556 = (-3)*IT_0555;
    const ccomplex_t IT_0557 = pow(m_W, 4);
    const ccomplex_t IT_0558 = s_12*IT_0557;
    const ccomplex_t IT_0559 = IT_0233*IT_0558;
    const ccomplex_t IT_0560 = (-0.25)*IT_0559;
    const ccomplex_t IT_0561 = IT_0236 + IT_0560;
    const ccomplex_t IT_0562 = 8*IT_0561;
    const ccomplex_t IT_0563 = (-0.5)*IT_0559;
    const ccomplex_t IT_0564 = IT_0235 + IT_0237 + IT_0563;
    const ccomplex_t IT_0565 = (-4)*IT_0564;
    const ccomplex_t IT_0566 = conj(IT_0322)*IT_0494;
    const ccomplex_t IT_0567 = 3*IT_0555;
    const ccomplex_t IT_0568 = (-2)*IT_0567;
    const ccomplex_t IT_0569 = conj(IT_0447)*IT_0493;
    const ccomplex_t IT_0570 = IT_0360*conj(IT_0444);
    const ccomplex_t IT_0571 = 2*IT_0488*(IT_0476*IT_0517 + IT_0505*IT_0527 +
       conj(IT_0355)*IT_0556 + 0.5*conj(IT_0487)*IT_0562 + 0.5*conj(IT_0488)
      *IT_0565 + 0.5*IT_0566 + 0.5*conj(IT_0350)*IT_0568 + IT_0569 + IT_0570);
    const ccomplex_t IT_0572 = 2*IT_0487*(IT_0486*IT_0546 + 0.5*IT_0476
      *IT_0550 + conj(IT_0350)*IT_0556 + 0.5*conj(IT_0488)*IT_0562 + 0.5*conj
      (IT_0487)*IT_0565 + 0.5*IT_0566 + IT_0519*IT_0568 + IT_0569 + IT_0570);
    const ccomplex_t IT_0573 = IT_0360*conj(IT_0362);
    const ccomplex_t IT_0574 = conj(IT_0364)*IT_0493;
    const ccomplex_t IT_0575 = 2*IT_0350*(IT_0361*IT_0502 + IT_0486*IT_0541 +
       0.5*IT_0476*IT_0551 + conj(IT_0487)*IT_0556 + IT_0519*IT_0562 + 0.5*conj
      (IT_0350)*IT_0565 + 0.5*conj(IT_0488)*IT_0568 + IT_0573 + IT_0574);
    const ccomplex_t IT_0576 = (conj(IT_0167) + conj(IT_0232))*IT_0290 + 2
      *IT_0283*IT_0479 + 2*IT_0322*IT_0499 + IT_0503 + IT_0504 + IT_0506 +
       IT_0509 + 2*IT_0167*(IT_0287*conj(IT_0362) + conj(IT_0322)*IT_0485 +
       IT_0511 + IT_0520 + IT_0528 + IT_0530 + IT_0532 + IT_0534 + IT_0536 +
       IT_0537 + IT_0542 + IT_0547) + 2*IT_0232*(IT_0287*conj(IT_0362) + conj
      (IT_0322)*IT_0485 + IT_0242*IT_0505 + IT_0511 + IT_0520 + IT_0528 +
       IT_0530 + IT_0532 + IT_0534 + IT_0536 + IT_0537 + IT_0542 + IT_0547) + 2*
      (IT_0462 + IT_0475)*IT_0554 + IT_0571 + IT_0572 + 2*IT_0355*(IT_0361
      *IT_0502 + IT_0505*IT_0518 + IT_0476*IT_0526 + conj(IT_0488)*IT_0556 + 0.5
      *conj(IT_0350)*IT_0562 + IT_0519*IT_0565 + 0.5*conj(IT_0487)*IT_0568 +
       IT_0573 + IT_0574) + IT_0575;
    return create_ccomplex_return(IT_0576);
}

