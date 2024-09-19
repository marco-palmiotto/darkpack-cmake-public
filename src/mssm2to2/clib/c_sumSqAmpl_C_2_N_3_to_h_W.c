#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_N_3_to_h_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_N_3_to_h_W(
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
    const ccomplex_t IT_0004 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0005 = IT_0003*IT_0004;
    const ccomplex_t IT_0006 = sin(alpha);
    const ccomplex_t IT_0007 = IT_0002*IT_0006;
    const ccomplex_t IT_0008 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0009 = IT_0007*IT_0008;
    const ccomplex_t IT_0010 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0011 = IT_0003*IT_0010;
    const ccomplex_t IT_0012 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0013 = IT_0007*IT_0012;
    const ccomplex_t IT_0014 = IT_0005 + IT_0009 + IT_0011 + IT_0013;
    const ccomplex_t IT_0015 = sin(theta_W);
    const ccomplex_t IT_0016 = cpow(IT_0015, -1);
    const ccomplex_t IT_0017 = IT_0000*IT_0016;
    const ccomplex_t IT_0018 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0006*IT_0016;
    const ccomplex_t IT_0021 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0024 = IT_0017*IT_0023;
    const ccomplex_t IT_0025 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0026 = IT_0020*IT_0025;
    const ccomplex_t IT_0027 = -IT_0019 + -IT_0022 + -IT_0024 + -IT_0026;
    const ccomplex_t IT_0028 = IT_0014 + IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (-0.5)*IT_0029;
    const ccomplex_t IT_0031 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0032 = IT_0016*IT_0031;
    const ccomplex_t IT_0033 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0034 = IT_0016*IT_0033;
    const ccomplex_t IT_0035 = 1.4142135623731*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*(IT_0032 + (-0.5)*IT_0035);
    const ccomplex_t IT_0037 = -IT_0036;
    const ccomplex_t IT_0038 = pow(m_h, 2);
    const ccomplex_t IT_0039 = pow(m_N_3, 2);
    const ccomplex_t IT_0040 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0041 = IT_0037*IT_0040;
    const ccomplex_t IT_0042 = (0 + _Complex_I*1)*IT_0041;
    const ccomplex_t IT_0043 = IT_0030*IT_0042;
    const ccomplex_t IT_0044 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0045 = IT_0003*IT_0044;
    const ccomplex_t IT_0046 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0047 = IT_0007*IT_0046;
    const ccomplex_t IT_0048 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0049 = IT_0003*IT_0048;
    const ccomplex_t IT_0050 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0051 = IT_0007*IT_0050;
    const ccomplex_t IT_0052 = IT_0045 + IT_0047 + IT_0049 + IT_0051;
    const ccomplex_t IT_0053 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0054 = IT_0017*IT_0053;
    const ccomplex_t IT_0055 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0056 = IT_0020*IT_0055;
    const ccomplex_t IT_0057 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0058 = IT_0017*IT_0057;
    const ccomplex_t IT_0059 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0060 = IT_0020*IT_0059;
    const ccomplex_t IT_0061 = -IT_0054 + -IT_0056 + -IT_0058 + -IT_0060;
    const ccomplex_t IT_0062 = IT_0052 + IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*IT_0062;
    const ccomplex_t IT_0064 = (-0.5)*IT_0063;
    const ccomplex_t IT_0065 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0066 = IT_0016*IT_0065;
    const ccomplex_t IT_0067 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0068 = IT_0016*IT_0067;
    const ccomplex_t IT_0069 = 1.4142135623731*IT_0068;
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*(IT_0066 + (-0.5)*IT_0069);
    const ccomplex_t IT_0071 = -IT_0070;
    const ccomplex_t IT_0072 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0064*IT_0074;
    const ccomplex_t IT_0076 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0077 = IT_0003*IT_0076;
    const ccomplex_t IT_0078 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0079 = IT_0007*IT_0078;
    const ccomplex_t IT_0080 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0081 = IT_0003*IT_0080;
    const ccomplex_t IT_0082 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0083 = IT_0007*IT_0082;
    const ccomplex_t IT_0084 = IT_0077 + IT_0079 + IT_0081 + IT_0083;
    const ccomplex_t IT_0085 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0086 = IT_0017*IT_0085;
    const ccomplex_t IT_0087 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0088 = IT_0020*IT_0087;
    const ccomplex_t IT_0089 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0090 = IT_0017*IT_0089;
    const ccomplex_t IT_0091 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0092 = IT_0020*IT_0091;
    const ccomplex_t IT_0093 = -IT_0086 + -IT_0088 + -IT_0090 + -IT_0092;
    const ccomplex_t IT_0094 = IT_0084 + IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = (-0.5)*IT_0095;
    const ccomplex_t IT_0097 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0098 = IT_0016*IT_0097;
    const ccomplex_t IT_0099 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0100 = IT_0016*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + (-0.5)*IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = cpow((-2)*s_23 + IT_0038 + IT_0039 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0105 = IT_0103*IT_0104;
    const ccomplex_t IT_0106 = (0 + _Complex_I*1)*IT_0105;
    const ccomplex_t IT_0107 = IT_0096*IT_0106;
    const ccomplex_t IT_0108 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0109 = IT_0003*IT_0108;
    const ccomplex_t IT_0110 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0111 = IT_0007*IT_0110;
    const ccomplex_t IT_0112 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0113 = IT_0017*IT_0112;
    const ccomplex_t IT_0114 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0115 = IT_0020*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*(IT_0109 + IT_0111 + 
      -IT_0113 + -IT_0115);
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0119 = IT_0016*IT_0118;
    const ccomplex_t IT_0120 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0121 = IT_0016*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*(IT_0119 + (-0.5)*IT_0122);
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = cpow((-2)*s_23 + IT_0038 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0126 = IT_0124*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*IT_0126;
    const ccomplex_t IT_0128 = IT_0117*IT_0127;
    const ccomplex_t IT_0129 = IT_0043 + IT_0075 + IT_0107 + IT_0128;
    const ccomplex_t IT_0130 = pow(m_W, 2);
    const ccomplex_t IT_0131 = s_12*IT_0130;
    const ccomplex_t IT_0132 = pow(m_W, -2);
    const ccomplex_t IT_0133 = pow(m_W, 4);
    const ccomplex_t IT_0134 = s_12*IT_0133;
    const ccomplex_t IT_0135 = IT_0132*IT_0134;
    const ccomplex_t IT_0136 = (-0.5)*IT_0135;
    const ccomplex_t IT_0137 = s_14*s_24;
    const ccomplex_t IT_0138 = (-2)*IT_0137;
    const ccomplex_t IT_0139 = IT_0131 + IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = m_C_2*m_N_3;
    const ccomplex_t IT_0141 = pow(s_24, 2);
    const ccomplex_t IT_0142 = IT_0132*IT_0141;
    const ccomplex_t IT_0143 = -IT_0039;
    const ccomplex_t IT_0144 = IT_0142 + IT_0143;
    const ccomplex_t IT_0145 = IT_0140*IT_0144;
    const ccomplex_t IT_0146 = cpow((-2)*s_13 + IT_0038 + (0 + _Complex_I*1)
      *m_C_2*Gamma_c2 + reg_prop, -1);
    const ccomplex_t IT_0147 = conj(V_u2)*conj(U_Wm2);
    const ccomplex_t IT_0148 = IT_0000*IT_0147;
    const ccomplex_t IT_0149 = conj(U_d2)*conj(V_Wp2);
    const ccomplex_t IT_0150 = IT_0006*IT_0149;
    const ccomplex_t IT_0151 = IT_0016*(IT_0148 + -IT_0150);
    const ccomplex_t IT_0152 = 1.4142135623731*e_em*IT_0151;
    const ccomplex_t IT_0153 = 0.5*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0124*IT_0154;
    const ccomplex_t IT_0156 = IT_0146*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = cos(beta);
    const ccomplex_t IT_0159 = sin(beta);
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*e_em*IT_0016*(IT_0000
      *IT_0158 + IT_0006*IT_0159);
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0016*IT_0159;
    const ccomplex_t IT_0163 = conj(N_d3)*e_em*conj(U_Wm2);
    const ccomplex_t IT_0164 = IT_0162*IT_0163;
    const ccomplex_t IT_0165 = (-2)*IT_0164;
    const ccomplex_t IT_0166 = IT_0002*IT_0159;
    const ccomplex_t IT_0167 = conj(N_B3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0168 = IT_0166*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W3)*conj(U_d2)*e_em;
    const ccomplex_t IT_0171 = IT_0162*IT_0170;
    const ccomplex_t IT_0172 = 1.4142135623731*IT_0171;
    const ccomplex_t IT_0173 = IT_0169 + IT_0172;
    const ccomplex_t IT_0174 = IT_0165 + IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = 0.5*IT_0175;
    const ccomplex_t IT_0177 = pow(m_C_2, 2);
    const ccomplex_t IT_0178 = cpow((-2)*s_12 + -IT_0039 + -IT_0177 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0179 = IT_0176*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0161*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*-1)*e_em*IT_0016*(IT_0006
      *IT_0158 + -IT_0000*IT_0159);
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = IT_0016*IT_0158;
    const ccomplex_t IT_0185 = IT_0163*IT_0184;
    const ccomplex_t IT_0186 = (-2)*IT_0185;
    const ccomplex_t IT_0187 = IT_0002*IT_0158;
    const ccomplex_t IT_0188 = IT_0167*IT_0187;
    const ccomplex_t IT_0189 = 1.4142135623731*IT_0188;
    const ccomplex_t IT_0190 = IT_0170*IT_0184;
    const ccomplex_t IT_0191 = 1.4142135623731*IT_0190;
    const ccomplex_t IT_0192 = IT_0189 + IT_0191;
    const ccomplex_t IT_0193 = IT_0186 + IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = (-0.5)*IT_0194;
    const ccomplex_t IT_0196 = cpow((-2)*s_12 + -IT_0039 + IT_0130 + -IT_0177 
      + -reg_prop, -1);
    const ccomplex_t IT_0197 = IT_0195*IT_0196;
    const ccomplex_t IT_0198 = (0 + _Complex_I*1)*IT_0197;
    const ccomplex_t IT_0199 = IT_0183*IT_0198;
    const ccomplex_t IT_0200 = conj(V_u1)*conj(U_Wm2);
    const ccomplex_t IT_0201 = IT_0000*IT_0200;
    const ccomplex_t IT_0202 = conj(U_d2)*conj(V_Wp1);
    const ccomplex_t IT_0203 = IT_0006*IT_0202;
    const ccomplex_t IT_0204 = IT_0016*(IT_0201 + -IT_0203);
    const ccomplex_t IT_0205 = 1.4142135623731*e_em*IT_0204;
    const ccomplex_t IT_0206 = 0.5*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0209 = IT_0016*IT_0208;
    const ccomplex_t IT_0210 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0211 = IT_0016*IT_0210;
    const ccomplex_t IT_0212 = 1.4142135623731*IT_0211;
    const ccomplex_t IT_0213 = (0 + _Complex_I*1)*(IT_0209 + (-0.5)*IT_0212);
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = IT_0207*IT_0214;
    const ccomplex_t IT_0216 = cpow((-2)*s_13 + IT_0038 + IT_0177 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0217 = IT_0215*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = (-2)*IT_0157 + -IT_0181 + -IT_0199 + (-2)
      *IT_0218;
    const ccomplex_t IT_0220 = V_u2*U_Wm2;
    const ccomplex_t IT_0221 = IT_0000*IT_0220;
    const ccomplex_t IT_0222 = U_d2*V_Wp2;
    const ccomplex_t IT_0223 = IT_0006*IT_0222;
    const ccomplex_t IT_0224 = IT_0016*(IT_0221 + -IT_0223);
    const ccomplex_t IT_0225 = 1.4142135623731*e_em*IT_0224;
    const ccomplex_t IT_0226 = 0.5*IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0229 = IT_0016*IT_0228;
    const ccomplex_t IT_0230 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0231 = IT_0016*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0229 + 0.5*IT_0232);
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0227*IT_0234;
    const ccomplex_t IT_0236 = IT_0146*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*IT_0236;
    const ccomplex_t IT_0238 = V_u2*U_Wm1;
    const ccomplex_t IT_0239 = IT_0000*IT_0238;
    const ccomplex_t IT_0240 = U_d1*V_Wp2;
    const ccomplex_t IT_0241 = IT_0006*IT_0240;
    const ccomplex_t IT_0242 = IT_0016*(IT_0239 + -IT_0241);
    const ccomplex_t IT_0243 = 1.4142135623731*e_em*IT_0242;
    const ccomplex_t IT_0244 = 0.5*IT_0243;
    const ccomplex_t IT_0245 = (0 + _Complex_I*1)*IT_0244;
    const ccomplex_t IT_0246 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0247 = IT_0016*IT_0246;
    const ccomplex_t IT_0248 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0249 = IT_0016*IT_0248;
    const ccomplex_t IT_0250 = 1.4142135623731*IT_0249;
    const ccomplex_t IT_0251 = (0 + _Complex_I*1)*(IT_0247 + 0.5*IT_0250);
    const ccomplex_t IT_0252 = -IT_0251;
    const ccomplex_t IT_0253 = IT_0245*IT_0252;
    const ccomplex_t IT_0254 = IT_0216*IT_0253;
    const ccomplex_t IT_0255 = (0 + _Complex_I*1)*IT_0254;
    const ccomplex_t IT_0256 = 2*IT_0237 + 2*IT_0255;
    const ccomplex_t IT_0257 = N_u3*e_em*V_Wp2;
    const ccomplex_t IT_0258 = IT_0162*IT_0257;
    const ccomplex_t IT_0259 = 2*IT_0258;
    const ccomplex_t IT_0260 = N_B3*V_u2*e_em;
    const ccomplex_t IT_0261 = IT_0166*IT_0260;
    const ccomplex_t IT_0262 = 1.4142135623731*IT_0261;
    const ccomplex_t IT_0263 = N_W3*V_u2*e_em;
    const ccomplex_t IT_0264 = IT_0162*IT_0263;
    const ccomplex_t IT_0265 = 1.4142135623731*IT_0264;
    const ccomplex_t IT_0266 = IT_0262 + IT_0265;
    const ccomplex_t IT_0267 = IT_0259 + IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*IT_0267;
    const ccomplex_t IT_0269 = 0.5*IT_0268;
    const ccomplex_t IT_0270 = IT_0196*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*IT_0270;
    const ccomplex_t IT_0272 = IT_0183*IT_0271;
    const ccomplex_t IT_0273 = IT_0184*IT_0257;
    const ccomplex_t IT_0274 = 2*IT_0273;
    const ccomplex_t IT_0275 = IT_0187*IT_0260;
    const ccomplex_t IT_0276 = 1.4142135623731*IT_0275;
    const ccomplex_t IT_0277 = IT_0184*IT_0263;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = IT_0276 + IT_0278;
    const ccomplex_t IT_0280 = IT_0274 + IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = 0.5*IT_0281;
    const ccomplex_t IT_0283 = IT_0178*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0161*IT_0284;
    const ccomplex_t IT_0286 = IT_0272 + IT_0285;
    const ccomplex_t IT_0287 = (-2)*IT_0286;
    const ccomplex_t IT_0288 = 0.5*IT_0287;
    const ccomplex_t IT_0289 = IT_0256 + IT_0288;
    const ccomplex_t IT_0290 = (-2)*conj(IT_0286);
    const ccomplex_t IT_0291 = 0.5*IT_0290;
    const ccomplex_t IT_0292 = conj(IT_0256) + IT_0291;
    const ccomplex_t IT_0293 = s_14*s_24*IT_0132;
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = s_12 + IT_0294;
    const ccomplex_t IT_0296 = IT_0140*IT_0295;
    const ccomplex_t IT_0297 = (-2)*IT_0043 + (-2)*IT_0075 + (-2)*IT_0107 + (
      -2)*IT_0128 + IT_0272 + IT_0285;
    const ccomplex_t IT_0298 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0299 = IT_0003*IT_0298;
    const ccomplex_t IT_0300 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0301 = IT_0007*IT_0300;
    const ccomplex_t IT_0302 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0303 = IT_0003*IT_0302;
    const ccomplex_t IT_0304 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0305 = IT_0007*IT_0304;
    const ccomplex_t IT_0306 = IT_0299 + IT_0301 + IT_0303 + IT_0305;
    const ccomplex_t IT_0307 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0308 = IT_0017*IT_0307;
    const ccomplex_t IT_0309 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0310 = IT_0020*IT_0309;
    const ccomplex_t IT_0311 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0312 = IT_0017*IT_0311;
    const ccomplex_t IT_0313 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0314 = IT_0020*IT_0313;
    const ccomplex_t IT_0315 = -IT_0308 + -IT_0310 + -IT_0312 + -IT_0314;
    const ccomplex_t IT_0316 = IT_0306 + IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*IT_0316;
    const ccomplex_t IT_0318 = (-0.5)*IT_0317;
    const ccomplex_t IT_0319 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0320 = IT_0016*IT_0319;
    const ccomplex_t IT_0321 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0322 = IT_0016*IT_0321;
    const ccomplex_t IT_0323 = 1.4142135623731*IT_0322;
    const ccomplex_t IT_0324 = (0 + _Complex_I*1)*(IT_0320 + 0.5*IT_0323);
    const ccomplex_t IT_0325 = -IT_0324;
    const ccomplex_t IT_0326 = IT_0040*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*IT_0326;
    const ccomplex_t IT_0328 = IT_0318*IT_0327;
    const ccomplex_t IT_0329 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0330 = IT_0003*IT_0329;
    const ccomplex_t IT_0331 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0332 = IT_0007*IT_0331;
    const ccomplex_t IT_0333 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0334 = IT_0017*IT_0333;
    const ccomplex_t IT_0335 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0336 = IT_0020*IT_0335;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*(IT_0330 + IT_0332 + 
      -IT_0334 + -IT_0336);
    const ccomplex_t IT_0338 = -IT_0337;
    const ccomplex_t IT_0339 = IT_0125*IT_0234;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*IT_0339;
    const ccomplex_t IT_0341 = IT_0338*IT_0340;
    const ccomplex_t IT_0342 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0343 = IT_0003*IT_0342;
    const ccomplex_t IT_0344 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0345 = IT_0007*IT_0344;
    const ccomplex_t IT_0346 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0347 = IT_0003*IT_0346;
    const ccomplex_t IT_0348 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0349 = IT_0007*IT_0348;
    const ccomplex_t IT_0350 = IT_0343 + IT_0345 + IT_0347 + IT_0349;
    const ccomplex_t IT_0351 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0352 = IT_0017*IT_0351;
    const ccomplex_t IT_0353 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0354 = IT_0020*IT_0353;
    const ccomplex_t IT_0355 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0356 = IT_0017*IT_0355;
    const ccomplex_t IT_0357 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0358 = IT_0020*IT_0357;
    const ccomplex_t IT_0359 = -IT_0352 + -IT_0354 + -IT_0356 + -IT_0358;
    const ccomplex_t IT_0360 = IT_0350 + IT_0359;
    const ccomplex_t IT_0361 = (0 + _Complex_I*1)*IT_0360;
    const ccomplex_t IT_0362 = (-0.5)*IT_0361;
    const ccomplex_t IT_0363 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0364 = IT_0016*IT_0363;
    const ccomplex_t IT_0365 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0366 = IT_0016*IT_0365;
    const ccomplex_t IT_0367 = 1.4142135623731*IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*(IT_0364 + 0.5*IT_0367);
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = IT_0104*IT_0369;
    const ccomplex_t IT_0371 = (0 + _Complex_I*1)*IT_0370;
    const ccomplex_t IT_0372 = IT_0362*IT_0371;
    const ccomplex_t IT_0373 = 2*IT_0328 + 2*IT_0341 + 2*IT_0372;
    const ccomplex_t IT_0374 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0375 = IT_0003*IT_0374;
    const ccomplex_t IT_0376 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0377 = IT_0007*IT_0376;
    const ccomplex_t IT_0378 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0379 = IT_0003*IT_0378;
    const ccomplex_t IT_0380 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0381 = IT_0007*IT_0380;
    const ccomplex_t IT_0382 = IT_0375 + IT_0377 + IT_0379 + IT_0381;
    const ccomplex_t IT_0383 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0384 = IT_0017*IT_0383;
    const ccomplex_t IT_0385 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0386 = IT_0020*IT_0385;
    const ccomplex_t IT_0387 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0388 = IT_0017*IT_0387;
    const ccomplex_t IT_0389 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0390 = IT_0020*IT_0389;
    const ccomplex_t IT_0391 = -IT_0384 + -IT_0386 + -IT_0388 + -IT_0390;
    const ccomplex_t IT_0392 = IT_0382 + IT_0391;
    const ccomplex_t IT_0393 = (0 + _Complex_I*1)*IT_0392;
    const ccomplex_t IT_0394 = (-0.5)*IT_0393;
    const ccomplex_t IT_0395 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0396 = IT_0016*IT_0395;
    const ccomplex_t IT_0397 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0398 = IT_0016*IT_0397;
    const ccomplex_t IT_0399 = 1.4142135623731*IT_0398;
    const ccomplex_t IT_0400 = (0 + _Complex_I*1)*(IT_0396 + 0.5*IT_0399);
    const ccomplex_t IT_0401 = -IT_0400;
    const ccomplex_t IT_0402 = IT_0072*IT_0401;
    const ccomplex_t IT_0403 = (0 + _Complex_I*1)*IT_0402;
    const ccomplex_t IT_0404 = IT_0394*IT_0403;
    const ccomplex_t IT_0405 = -IT_0181 + -IT_0199 + (-2)*IT_0404;
    const ccomplex_t IT_0406 = (-2)*IT_0256;
    const ccomplex_t IT_0407 = 0.5*IT_0406;
    const ccomplex_t IT_0408 = IT_0286 + IT_0407;
    const ccomplex_t IT_0409 = (-2)*conj(IT_0256);
    const ccomplex_t IT_0410 = 0.5*IT_0409;
    const ccomplex_t IT_0411 = conj(IT_0286) + IT_0410;
    const ccomplex_t IT_0412 = -IT_0237 + -IT_0255;
    const ccomplex_t IT_0413 = (-0.25)*IT_0135;
    const ccomplex_t IT_0414 = IT_0137 + IT_0413;
    const ccomplex_t IT_0415 = s_24*IT_0130;
    const ccomplex_t IT_0416 = IT_0132*IT_0415;
    const ccomplex_t IT_0417 = 0.5*IT_0416;
    const ccomplex_t IT_0418 = s_24 + IT_0417;
    const ccomplex_t IT_0419 = m_C_2*IT_0418;
    const ccomplex_t IT_0420 = 2*IT_0419;
    const ccomplex_t IT_0421 = IT_0154*IT_0234;
    const ccomplex_t IT_0422 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0423 = IT_0146*IT_0421*IT_0422;
    const ccomplex_t IT_0424 = m_N_3*IT_0117;
    const ccomplex_t IT_0425 = IT_0340*IT_0424;
    const ccomplex_t IT_0426 = (0 + _Complex_I*-1)*m_W*e_em*IT_0016*(IT_0006
      *IT_0158 + -IT_0000*IT_0159);
    const ccomplex_t IT_0427 = IT_0234*IT_0426;
    const ccomplex_t IT_0428 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0039 + -IT_0177 + -reg_prop, -1);
    const ccomplex_t IT_0429 = IT_0427*IT_0428;
    const ccomplex_t IT_0430 = (0 + _Complex_I*1)*IT_0429;
    const ccomplex_t IT_0431 = IT_0207*IT_0252;
    const ccomplex_t IT_0432 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0433 = IT_0216*IT_0431*IT_0432;
    const ccomplex_t IT_0434 = m_N_1*IT_0096;
    const ccomplex_t IT_0435 = IT_0371*IT_0434;
    const ccomplex_t IT_0436 = m_N_2*IT_0030;
    const ccomplex_t IT_0437 = IT_0327*IT_0436;
    const ccomplex_t IT_0438 = IT_0423 + -IT_0425 + -IT_0430 + IT_0433 + 
      -IT_0435 + -IT_0437;
    const ccomplex_t IT_0439 = (-2)*IT_0438;
    const ccomplex_t IT_0440 = m_C_2*IT_0096;
    const ccomplex_t IT_0441 = IT_0106*IT_0440;
    const ccomplex_t IT_0442 = m_C_2*IT_0030;
    const ccomplex_t IT_0443 = IT_0042*IT_0442;
    const ccomplex_t IT_0444 = m_N_3*IT_0218;
    const ccomplex_t IT_0445 = m_N_4*IT_0064;
    const ccomplex_t IT_0446 = IT_0403*IT_0445;
    const ccomplex_t IT_0447 = m_C_2*IT_0117;
    const ccomplex_t IT_0448 = IT_0127*IT_0447;
    const ccomplex_t IT_0449 = m_N_3*IT_0157;
    const ccomplex_t IT_0450 = m_C_2*IT_0064;
    const ccomplex_t IT_0451 = IT_0074*IT_0450;
    const ccomplex_t IT_0452 = IT_0441 + IT_0443 + -IT_0444 + -IT_0446 +
       IT_0448 + -IT_0449 + IT_0451;
    const ccomplex_t IT_0453 = (-2)*IT_0452;
    const ccomplex_t IT_0454 = IT_0439 + IT_0453;
    const ccomplex_t IT_0455 = IT_0157 + IT_0218;
    const ccomplex_t IT_0456 = IT_0130*IT_0140;
    const ccomplex_t IT_0457 = (-3)*IT_0456;
    const ccomplex_t IT_0458 = 3*IT_0456;
    const ccomplex_t IT_0459 = -IT_0328 + -IT_0341 + -IT_0372;
    const ccomplex_t IT_0460 = -conj(IT_0459);
    const ccomplex_t IT_0461 = conj(IT_0404) + IT_0460;
    const ccomplex_t IT_0462 = (-2)*conj(IT_0438);
    const ccomplex_t IT_0463 = (-2)*conj(IT_0452);
    const ccomplex_t IT_0464 = IT_0462 + IT_0463;
    const ccomplex_t IT_0465 = s_12*s_14;
    const ccomplex_t IT_0466 = s_24*IT_0177;
    const ccomplex_t IT_0467 = -IT_0466;
    const ccomplex_t IT_0468 = IT_0465 + IT_0467;
    const ccomplex_t IT_0469 = (-2)*IT_0405;
    const ccomplex_t IT_0470 = 0.5*IT_0469;
    const ccomplex_t IT_0471 = IT_0373 + IT_0470;
    const ccomplex_t IT_0472 = (-2)*conj(IT_0373);
    const ccomplex_t IT_0473 = 0.5*IT_0472;
    const ccomplex_t IT_0474 = conj(IT_0405) + IT_0473;
    const ccomplex_t IT_0475 = IT_0129 + IT_0412;
    const ccomplex_t IT_0476 = conj(IT_0129) + conj(IT_0412);
    const ccomplex_t IT_0477 = (-2)*IT_0373;
    const ccomplex_t IT_0478 = IT_0455 + IT_0459;
    const ccomplex_t IT_0479 = (-2)*conj(IT_0405);
    const ccomplex_t IT_0480 = s_12*s_24;
    const ccomplex_t IT_0481 = s_14*IT_0039;
    const ccomplex_t IT_0482 = -IT_0481;
    const ccomplex_t IT_0483 = IT_0480 + IT_0482;
    const ccomplex_t IT_0484 = (-2)*(conj(IT_0219)*(IT_0404 + -IT_0455 + 
      -IT_0459) + -IT_0219*(conj(IT_0455) + -IT_0461) + -(conj(IT_0256) + 0.5
      *IT_0290)*IT_0475 + -IT_0289*IT_0476)*IT_0483;
    const ccomplex_t IT_0485 = pow(s_14, 2);
    const ccomplex_t IT_0486 = IT_0132*IT_0485;
    const ccomplex_t IT_0487 = -IT_0177;
    const ccomplex_t IT_0488 = IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = m_N_3*IT_0488;
    const ccomplex_t IT_0490 = s_14*IT_0130;
    const ccomplex_t IT_0491 = IT_0132*IT_0490;
    const ccomplex_t IT_0492 = 0.5*IT_0491;
    const ccomplex_t IT_0493 = s_14 + IT_0492;
    const ccomplex_t IT_0494 = m_N_3*IT_0493;
    const ccomplex_t IT_0495 = 2*IT_0494;
    const ccomplex_t IT_0496 = conj(IT_0438) + conj(IT_0452);
    const ccomplex_t IT_0497 = m_C_2*IT_0318;
    const ccomplex_t IT_0498 = IT_0327*IT_0497;
    const ccomplex_t IT_0499 = m_C_2*IT_0362;
    const ccomplex_t IT_0500 = IT_0371*IT_0499;
    const ccomplex_t IT_0501 = m_N_3*IT_0255;
    const ccomplex_t IT_0502 = m_C_2*IT_0338;
    const ccomplex_t IT_0503 = IT_0340*IT_0502;
    const ccomplex_t IT_0504 = m_N_3*IT_0237;
    const ccomplex_t IT_0505 = IT_0498 + IT_0500 + -IT_0501 + IT_0503 + 
      -IT_0504;
    const ccomplex_t IT_0506 = (-2)*conj(IT_0505);
    const ccomplex_t IT_0507 = IT_0124*IT_0227;
    const ccomplex_t IT_0508 = IT_0146*IT_0422*IT_0507;
    const ccomplex_t IT_0509 = IT_0124*IT_0426;
    const ccomplex_t IT_0510 = IT_0428*IT_0509;
    const ccomplex_t IT_0511 = (0 + _Complex_I*1)*IT_0510;
    const ccomplex_t IT_0512 = m_N_3*IT_0338;
    const ccomplex_t IT_0513 = IT_0127*IT_0512;
    const ccomplex_t IT_0514 = IT_0214*IT_0245;
    const ccomplex_t IT_0515 = IT_0216*IT_0432*IT_0514;
    const ccomplex_t IT_0516 = m_N_1*IT_0362;
    const ccomplex_t IT_0517 = IT_0106*IT_0516;
    const ccomplex_t IT_0518 = m_N_2*IT_0318;
    const ccomplex_t IT_0519 = IT_0042*IT_0518;
    const ccomplex_t IT_0520 = m_C_2*IT_0394;
    const ccomplex_t IT_0521 = IT_0403*IT_0520;
    const ccomplex_t IT_0522 = m_N_4*IT_0394;
    const ccomplex_t IT_0523 = IT_0074*IT_0522;
    const ccomplex_t IT_0524 = IT_0508 + -IT_0511 + -IT_0513 + IT_0515 + 
      -IT_0517 + -IT_0519 + IT_0521 + -IT_0523;
    const ccomplex_t IT_0525 = (-2)*conj(IT_0524);
    const ccomplex_t IT_0526 = IT_0506 + IT_0525;
    const ccomplex_t IT_0527 = (-2)*IT_0524;
    const ccomplex_t IT_0528 = (-0.25)*IT_0416;
    const ccomplex_t IT_0529 = s_24 + IT_0528;
    const ccomplex_t IT_0530 = m_C_2*IT_0529;
    const ccomplex_t IT_0531 = (-4)*IT_0530;
    const ccomplex_t IT_0532 = (-2)*IT_0129;
    const ccomplex_t IT_0533 = (-0.25)*IT_0491;
    const ccomplex_t IT_0534 = s_14 + IT_0533;
    const ccomplex_t IT_0535 = m_N_3*IT_0534;
    const ccomplex_t IT_0536 = (-4)*IT_0535;
    const ccomplex_t IT_0537 = (-2)*conj(IT_0129);
    const ccomplex_t IT_0538 = conj(IT_0505) + conj(IT_0524);
    const ccomplex_t IT_0539 = 0.5*IT_0479;
    const ccomplex_t IT_0540 = conj(IT_0373) + IT_0539;
    const ccomplex_t IT_0541 = m_C_2*IT_0144;
    const ccomplex_t IT_0542 = IT_0131 + IT_0138;
    const ccomplex_t IT_0543 = IT_0132*IT_0542;
    const ccomplex_t IT_0544 = (-0.5)*IT_0543;
    const ccomplex_t IT_0545 = s_12 + IT_0544;
    const ccomplex_t IT_0546 = m_C_2*IT_0295;
    const ccomplex_t IT_0547 = (-2)*conj(IT_0297);
    const ccomplex_t IT_0548 = m_N_3*IT_0295;
    const ccomplex_t IT_0549 = 0.333333333333333*conj(IT_0219);
    const ccomplex_t IT_0550 = 2*IT_0452*(conj(IT_0459)*IT_0495 + conj(IT_0129
      )*IT_0531 + 3*IT_0140*IT_0538 + IT_0489*IT_0540 + IT_0411*IT_0541 + 2
      *IT_0496*IT_0545 + 0.5*IT_0546*IT_0547 + 3*IT_0548*IT_0549);
    const ccomplex_t IT_0551 = (-4)*IT_0414*conj(IT_0455);
    const ccomplex_t IT_0552 = (-2)*IT_0139*IT_0461;
    const ccomplex_t IT_0553 = (-2)*IT_0297;
    const ccomplex_t IT_0554 = 0.166666666666667*IT_0489;
    const ccomplex_t IT_0555 = 0.333333333333333*IT_0546;
    const ccomplex_t IT_0556 = 0.333333333333333*IT_0548;
    const ccomplex_t IT_0557 = 6*IT_0538*(IT_0140*IT_0438 + 0.333333333333333
      *IT_0420*IT_0455 + 0.333333333333333*IT_0404*IT_0531 + 0.333333333333333
      *IT_0412*IT_0536 + 0.333333333333333*IT_0219*IT_0541 + IT_0553*IT_0554 +
       IT_0471*IT_0555 + IT_0408*IT_0556);
    const ccomplex_t IT_0558 = IT_0140*IT_0488;
    const ccomplex_t IT_0559 = 0.5*conj(IT_0405);
    const ccomplex_t IT_0560 = s_12*IT_0488;
    const ccomplex_t IT_0561 = s_12*IT_0295;
    const ccomplex_t IT_0562 = 0.5*conj(IT_0286);
    const ccomplex_t IT_0563 = s_12*IT_0144;
    const ccomplex_t IT_0564 = 0.5*IT_0546;
    const ccomplex_t IT_0565 = 0.333333333333333*IT_0420*conj(IT_0455);
    const ccomplex_t IT_0566 = IT_0541*IT_0549 + IT_0565;
    const ccomplex_t IT_0567 = IT_0547*IT_0554;
    const ccomplex_t IT_0568 = 0.333333333333333*conj(IT_0412)*IT_0536;
    const ccomplex_t IT_0569 = IT_0540*IT_0555;
    const ccomplex_t IT_0570 = IT_0411*IT_0556;
    const ccomplex_t IT_0571 = 6*conj(IT_0438);
    const ccomplex_t IT_0572 = 0.166666666666667*IT_0571;
    const ccomplex_t IT_0573 = conj(IT_0452) + IT_0572;
    const ccomplex_t IT_0574 = IT_0140*IT_0573;
    const ccomplex_t IT_0575 = (-2)*IT_0505*(conj(IT_0129)*IT_0495 + -IT_0461
      *IT_0531 + (-2)*IT_0538*IT_0545 + (-3)*IT_0566 + (-3)*IT_0567 + (-3)
      *IT_0568 + (-3)*IT_0569 + (-3)*IT_0570 + (-3)*IT_0574);
    const ccomplex_t IT_0576 = IT_0553*IT_0564;
    const ccomplex_t IT_0577 = 4*IT_0438;
    const ccomplex_t IT_0578 = 0.5*IT_0545*IT_0577;
    const ccomplex_t IT_0579 = IT_0129*IT_0531;
    const ccomplex_t IT_0580 = IT_0408*IT_0541;
    const ccomplex_t IT_0581 = (-4)*IT_0129*conj(IT_0129)*IT_0139 + 2*IT_0145*
      (conj(IT_0219)*IT_0289 + IT_0219*IT_0292) + 2*IT_0296*(conj(IT_0219)
      *IT_0297 + IT_0219*conj(IT_0297) + IT_0292*IT_0373 + IT_0289*conj(IT_0373)
       + conj(IT_0405)*IT_0408 + IT_0405*IT_0411) + 8*conj(IT_0412)*(IT_0129
      *IT_0414 + 0.125*IT_0420*IT_0454) + (-4)*IT_0412*(IT_0139*conj(IT_0412) + 
      (-2)*conj(IT_0129)*IT_0414 + -1./2*conj(IT_0455)*IT_0457 + -1./2*IT_0458
      *IT_0461 + (-0.25)*IT_0420*IT_0464) + 2*IT_0468*((conj(IT_0455) + conj
      (IT_0459))*IT_0471 + IT_0404*IT_0474 + conj(IT_0297)*IT_0475 + IT_0297
      *IT_0476 + conj(IT_0404)*(IT_0405 + 0.5*IT_0477) + IT_0478*(conj(IT_0373) 
      + 0.5*IT_0479)) + IT_0484 + (conj(IT_0405)*IT_0439 + IT_0405*IT_0462)
      *IT_0489 + 2*IT_0495*(0.5*conj(IT_0404)*IT_0454 + IT_0438*conj(IT_0459) +
       0.5*IT_0404*IT_0464 + IT_0459*IT_0496 + 0.5*IT_0129*IT_0526 + 0.5*conj
      (IT_0129)*IT_0527) + 2*conj(IT_0459)*(IT_0129*IT_0457 + 0.5*IT_0527
      *IT_0531) + conj(IT_0404)*IT_0457*IT_0532 + conj(IT_0455)*(IT_0458*IT_0532
       + IT_0454*IT_0536) + (-4)*IT_0455*(IT_0139*conj(IT_0455) + -1./2*conj
      (IT_0412)*IT_0457 + 2*IT_0414*IT_0461 + (-0.25)*IT_0464*IT_0536 + (-0.25)
      *IT_0458*IT_0537) + IT_0550 + 2*IT_0404*(conj(IT_0412)*IT_0458 + 0.5
      *IT_0457*IT_0537 + IT_0551 + IT_0552) + 2*IT_0459*(conj(IT_0129)*IT_0457 +
       -conj(IT_0412)*IT_0458 + 0.5*IT_0526*IT_0531 + -IT_0551 + -IT_0552) +
       IT_0557 + 2*IT_0558*(conj(IT_0297)*IT_0373 + IT_0297*conj(IT_0373) + 0.5
      *IT_0405*IT_0547 + IT_0553*IT_0559) + 2*(IT_0297*conj(IT_0297) + IT_0373
      *conj(IT_0373) + IT_0405*IT_0474 + IT_0477*IT_0559)*IT_0560 + 2*IT_0561*
      (conj(IT_0256)*IT_0297 + IT_0256*conj(IT_0297) + conj(IT_0219)*IT_0471 +
       IT_0219*IT_0540 + 0.5*IT_0286*IT_0547 + IT_0553*IT_0562) + 2*(IT_0219
      *conj(IT_0219) + IT_0256*conj(IT_0256) + IT_0286*IT_0411 + IT_0406*IT_0562
      )*IT_0563 + 2*IT_0438*(conj(IT_0373)*IT_0489 + conj(IT_0129)*IT_0531 +
       IT_0411*IT_0541 + conj(IT_0219)*IT_0548 + IT_0547*IT_0564) + 2*IT_0524*
      (conj(IT_0404)*IT_0531 + 2*IT_0538*IT_0545 + 3*IT_0566 + 3*IT_0567 + 3
      *IT_0568 + 3*IT_0569 + 3*IT_0570 + 3*IT_0574) + IT_0575 + 2*conj(IT_0438)*
      (IT_0373*IT_0489 + IT_0219*IT_0548 + IT_0576 + IT_0578 + IT_0579 + IT_0580
      ) + 2*conj(IT_0452)*(IT_0471*IT_0489 + IT_0219*IT_0548 + IT_0576 + IT_0578
       + IT_0579 + IT_0580);
    return create_ccomplex_return(IT_0581);
}

