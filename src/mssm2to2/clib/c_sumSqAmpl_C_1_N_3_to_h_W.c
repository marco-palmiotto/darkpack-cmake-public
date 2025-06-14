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
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = s_12*IT_0000;
    const ccomplex_t IT_0002 = pow(m_W, -2);
    const ccomplex_t IT_0003 = pow(m_W, 4);
    const ccomplex_t IT_0004 = s_12*IT_0003;
    const ccomplex_t IT_0005 = IT_0002*IT_0004;
    const ccomplex_t IT_0006 = (-0.5)*IT_0005;
    const ccomplex_t IT_0007 = s_14*s_24;
    const ccomplex_t IT_0008 = (-2)*IT_0007;
    const ccomplex_t IT_0009 = IT_0001 + IT_0006 + IT_0008;
    const ccomplex_t IT_0010 = pow(m_h, 2);
    const ccomplex_t IT_0011 = cpow((-2)*s_13 + IT_0010 + (0 + _Complex_I*1)
      *m_C_1*Gamma_c1 + reg_prop, -1);
    const ccomplex_t IT_0012 = cos(alpha);
    const ccomplex_t IT_0013 = V_u1*U_Wm1;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = sin(alpha);
    const ccomplex_t IT_0016 = U_d1*V_Wp1;
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = sin(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = (IT_0014 + -IT_0017)*IT_0019;
    const ccomplex_t IT_0021 = 1.4142135623731*e_em*IT_0020;
    const ccomplex_t IT_0022 = 0.5*IT_0021;
    const ccomplex_t IT_0023 = (0 + _Complex_I*1)*IT_0022;
    const ccomplex_t IT_0024 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0027 = IT_0019*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + 0.5*IT_0028);
    const ccomplex_t IT_0030 = -IT_0029;
    const ccomplex_t IT_0031 = IT_0023*IT_0030;
    const ccomplex_t IT_0032 = IT_0011*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*IT_0032;
    const ccomplex_t IT_0034 = V_u1*U_Wm2;
    const ccomplex_t IT_0035 = IT_0012*IT_0034;
    const ccomplex_t IT_0036 = U_d2*V_Wp1;
    const ccomplex_t IT_0037 = IT_0015*IT_0036;
    const ccomplex_t IT_0038 = IT_0019*(IT_0035 + -IT_0037);
    const ccomplex_t IT_0039 = 1.4142135623731*e_em*IT_0038;
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*IT_0040;
    const ccomplex_t IT_0042 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0043 = IT_0019*IT_0042;
    const ccomplex_t IT_0044 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0045 = IT_0019*IT_0044;
    const ccomplex_t IT_0046 = 1.4142135623731*IT_0045;
    const ccomplex_t IT_0047 = (0 + _Complex_I*1)*(IT_0043 + 0.5*IT_0046);
    const ccomplex_t IT_0048 = -IT_0047;
    const ccomplex_t IT_0049 = IT_0041*IT_0048;
    const ccomplex_t IT_0050 = pow(m_C_1, 2);
    const ccomplex_t IT_0051 = cpow((-2)*s_13 + IT_0010 + IT_0050 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0052 = IT_0049*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*IT_0052;
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = conj(V_u1)*conj(U_Wm1);
    const ccomplex_t IT_0056 = IT_0012*IT_0055;
    const ccomplex_t IT_0057 = conj(U_d1)*conj(V_Wp1);
    const ccomplex_t IT_0058 = IT_0015*IT_0057;
    const ccomplex_t IT_0059 = IT_0019*(IT_0056 + -IT_0058);
    const ccomplex_t IT_0060 = 1.4142135623731*e_em*IT_0059;
    const ccomplex_t IT_0061 = 0.5*IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0064 = IT_0019*IT_0063;
    const ccomplex_t IT_0065 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0066 = IT_0019*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = (0 + _Complex_I*1)*(IT_0064 + (-0.5)*IT_0067);
    const ccomplex_t IT_0069 = -IT_0068;
    const ccomplex_t IT_0070 = IT_0062*IT_0069;
    const ccomplex_t IT_0071 = IT_0011*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*IT_0071;
    const ccomplex_t IT_0073 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0074 = IT_0012*IT_0073;
    const ccomplex_t IT_0075 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0076 = IT_0015*IT_0075;
    const ccomplex_t IT_0077 = IT_0019*(IT_0074 + -IT_0076);
    const ccomplex_t IT_0078 = 1.4142135623731*e_em*IT_0077;
    const ccomplex_t IT_0079 = 0.5*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0082 = IT_0019*IT_0081;
    const ccomplex_t IT_0083 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0084 = IT_0019*IT_0083;
    const ccomplex_t IT_0085 = 1.4142135623731*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*(IT_0082 + (-0.5)*IT_0085);
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0080*IT_0087;
    const ccomplex_t IT_0089 = IT_0051*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = IT_0072 + IT_0090;
    const ccomplex_t IT_0092 = conj(IT_0072) + conj(IT_0090);
    const ccomplex_t IT_0093 = s_12*s_14;
    const ccomplex_t IT_0094 = s_24*IT_0050;
    const ccomplex_t IT_0095 = -IT_0094;
    const ccomplex_t IT_0096 = IT_0093 + IT_0095;
    const ccomplex_t IT_0097 = cos(theta_W);
    const ccomplex_t IT_0098 = cpow(IT_0097, -1);
    const ccomplex_t IT_0099 = IT_0012*IT_0098;
    const ccomplex_t IT_0100 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0101 = IT_0099*IT_0100;
    const ccomplex_t IT_0102 = IT_0015*IT_0098;
    const ccomplex_t IT_0103 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0104 = IT_0102*IT_0103;
    const ccomplex_t IT_0105 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0106 = IT_0099*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0108 = IT_0102*IT_0107;
    const ccomplex_t IT_0109 = IT_0101 + IT_0104 + IT_0106 + IT_0108;
    const ccomplex_t IT_0110 = IT_0012*IT_0019;
    const ccomplex_t IT_0111 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = IT_0015*IT_0019;
    const ccomplex_t IT_0114 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0115 = IT_0113*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0117 = IT_0110*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0119 = IT_0113*IT_0118;
    const ccomplex_t IT_0120 = -IT_0112 + -IT_0115 + -IT_0117 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0109 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = pow(m_N_3, 2);
    const ccomplex_t IT_0125 = cpow((-2)*s_23 + IT_0010 + IT_0124 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0126 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0127 = IT_0019*IT_0126;
    const ccomplex_t IT_0128 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0129 = IT_0019*IT_0128;
    const ccomplex_t IT_0130 = 1.4142135623731*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*(IT_0127 + 0.5*IT_0130);
    const ccomplex_t IT_0132 = -IT_0131;
    const ccomplex_t IT_0133 = IT_0125*IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*IT_0133;
    const ccomplex_t IT_0135 = IT_0123*IT_0134;
    const ccomplex_t IT_0136 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0137 = IT_0099*IT_0136;
    const ccomplex_t IT_0138 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0139 = IT_0102*IT_0138;
    const ccomplex_t IT_0140 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0141 = IT_0099*IT_0140;
    const ccomplex_t IT_0142 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0143 = IT_0102*IT_0142;
    const ccomplex_t IT_0144 = IT_0137 + IT_0139 + IT_0141 + IT_0143;
    const ccomplex_t IT_0145 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0146 = IT_0110*IT_0145;
    const ccomplex_t IT_0147 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0148 = IT_0113*IT_0147;
    const ccomplex_t IT_0149 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0150 = IT_0110*IT_0149;
    const ccomplex_t IT_0151 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0152 = IT_0113*IT_0151;
    const ccomplex_t IT_0153 = -IT_0146 + -IT_0148 + -IT_0150 + -IT_0152;
    const ccomplex_t IT_0154 = IT_0144 + IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*IT_0154;
    const ccomplex_t IT_0156 = (-0.5)*IT_0155;
    const ccomplex_t IT_0157 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0158 = IT_0019*IT_0157;
    const ccomplex_t IT_0159 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0160 = IT_0019*IT_0159;
    const ccomplex_t IT_0161 = 1.4142135623731*IT_0160;
    const ccomplex_t IT_0162 = (0 + _Complex_I*1)*(IT_0158 + 0.5*IT_0161);
    const ccomplex_t IT_0163 = -IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_23 + IT_0010 + IT_0124 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0165 = IT_0163*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*IT_0165;
    const ccomplex_t IT_0167 = IT_0156*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0169 = IT_0099*IT_0168;
    const ccomplex_t IT_0170 = conj(N_B3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0171 = IT_0102*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W3)*conj(N_u3)*e_em;
    const ccomplex_t IT_0173 = IT_0110*IT_0172;
    const ccomplex_t IT_0174 = conj(N_W3)*conj(N_d3)*e_em;
    const ccomplex_t IT_0175 = IT_0113*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*(IT_0169 + IT_0171 + 
      -IT_0173 + -IT_0175);
    const ccomplex_t IT_0177 = -IT_0176;
    const ccomplex_t IT_0178 = cpow((-2)*s_23 + IT_0010 + (0 + _Complex_I*1)
      *m_N_3*Gamma_o3 + reg_prop, -1);
    const ccomplex_t IT_0179 = IT_0030*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = IT_0177*IT_0180;
    const ccomplex_t IT_0182 = 2*IT_0135 + 2*IT_0167 + IT_0181;
    const ccomplex_t IT_0183 = 2*IT_0072;
    const ccomplex_t IT_0184 = 2*IT_0090;
    const ccomplex_t IT_0185 = IT_0183 + IT_0184;
    const ccomplex_t IT_0186 = cos(beta);
    const ccomplex_t IT_0187 = sin(beta);
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*e_em*IT_0019*(IT_0012
      *IT_0186 + IT_0015*IT_0187);
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0019*IT_0187;
    const ccomplex_t IT_0191 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0192 = IT_0190*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = IT_0098*IT_0187;
    const ccomplex_t IT_0195 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0196 = IT_0194*IT_0195;
    const ccomplex_t IT_0197 = 1.4142135623731*IT_0196;
    const ccomplex_t IT_0198 = IT_0193 + IT_0197;
    const ccomplex_t IT_0199 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0200 = IT_0190*IT_0199;
    const ccomplex_t IT_0201 = (-2)*IT_0200;
    const ccomplex_t IT_0202 = IT_0198 + IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = cpow((-2)*s_12 + -IT_0050 + -IT_0124 + m_Hp*
      (m_Hp + (0 + _Complex_I*-1)*Gamma_Hp) + -reg_prop, -1);
    const ccomplex_t IT_0206 = IT_0204*IT_0205;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = IT_0189*IT_0207;
    const ccomplex_t IT_0209 = conj(N_B3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0210 = IT_0099*IT_0209;
    const ccomplex_t IT_0211 = conj(N_B3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0212 = IT_0102*IT_0211;
    const ccomplex_t IT_0213 = conj(N_B4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0214 = IT_0099*IT_0213;
    const ccomplex_t IT_0215 = conj(N_B4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0216 = IT_0102*IT_0215;
    const ccomplex_t IT_0217 = IT_0210 + IT_0212 + IT_0214 + IT_0216;
    const ccomplex_t IT_0218 = conj(N_W3)*conj(N_u4)*e_em;
    const ccomplex_t IT_0219 = IT_0110*IT_0218;
    const ccomplex_t IT_0220 = conj(N_W3)*conj(N_d4)*e_em;
    const ccomplex_t IT_0221 = IT_0113*IT_0220;
    const ccomplex_t IT_0222 = conj(N_W4)*conj(N_u3)*e_em;
    const ccomplex_t IT_0223 = IT_0110*IT_0222;
    const ccomplex_t IT_0224 = conj(N_W4)*conj(N_d3)*e_em;
    const ccomplex_t IT_0225 = IT_0113*IT_0224;
    const ccomplex_t IT_0226 = -IT_0219 + -IT_0221 + -IT_0223 + -IT_0225;
    const ccomplex_t IT_0227 = IT_0217 + IT_0226;
    const ccomplex_t IT_0228 = (0 + _Complex_I*1)*IT_0227;
    const ccomplex_t IT_0229 = (-0.5)*IT_0228;
    const ccomplex_t IT_0230 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0231 = IT_0019*IT_0230;
    const ccomplex_t IT_0232 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0233 = IT_0019*IT_0232;
    const ccomplex_t IT_0234 = 1.4142135623731*IT_0233;
    const ccomplex_t IT_0235 = (0 + _Complex_I*1)*(IT_0231 + 0.5*IT_0234);
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = cpow((-2)*s_23 + IT_0010 + IT_0124 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0238 = IT_0236*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*IT_0238;
    const ccomplex_t IT_0240 = IT_0229*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*-1)*e_em*IT_0019*(IT_0015
      *IT_0186 + -IT_0012*IT_0187);
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = cpow((-2)*s_12 + IT_0000 + -IT_0050 + -IT_0124 
      + -reg_prop, -1);
    const ccomplex_t IT_0244 = IT_0019*IT_0186;
    const ccomplex_t IT_0245 = IT_0199*IT_0244;
    const ccomplex_t IT_0246 = (-2)*IT_0245;
    const ccomplex_t IT_0247 = IT_0098*IT_0186;
    const ccomplex_t IT_0248 = IT_0195*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = IT_0191*IT_0244;
    const ccomplex_t IT_0251 = 1.4142135623731*IT_0250;
    const ccomplex_t IT_0252 = IT_0249 + IT_0251;
    const ccomplex_t IT_0253 = IT_0246 + IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = (-0.5)*IT_0254;
    const ccomplex_t IT_0256 = IT_0243*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*IT_0256;
    const ccomplex_t IT_0258 = IT_0242*IT_0257;
    const ccomplex_t IT_0259 = -IT_0181 + -IT_0208 + (-2)*IT_0240 + -IT_0258;
    const ccomplex_t IT_0260 = (-2)*IT_0072;
    const ccomplex_t IT_0261 = (-2)*IT_0090;
    const ccomplex_t IT_0262 = IT_0260 + IT_0261;
    const ccomplex_t IT_0263 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0264 = IT_0244*IT_0263;
    const ccomplex_t IT_0265 = 2*IT_0264;
    const ccomplex_t IT_0266 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0267 = IT_0247*IT_0266;
    const ccomplex_t IT_0268 = 1.4142135623731*IT_0267;
    const ccomplex_t IT_0269 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0270 = IT_0244*IT_0269;
    const ccomplex_t IT_0271 = 1.4142135623731*IT_0270;
    const ccomplex_t IT_0272 = IT_0268 + IT_0271;
    const ccomplex_t IT_0273 = IT_0265 + IT_0272;
    const ccomplex_t IT_0274 = (0 + _Complex_I*1)*IT_0273;
    const ccomplex_t IT_0275 = 0.5*IT_0274;
    const ccomplex_t IT_0276 = IT_0205*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = IT_0189*IT_0277;
    const ccomplex_t IT_0279 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0280 = IT_0099*IT_0279;
    const ccomplex_t IT_0281 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0282 = IT_0102*IT_0281;
    const ccomplex_t IT_0283 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0284 = IT_0099*IT_0283;
    const ccomplex_t IT_0285 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0286 = IT_0102*IT_0285;
    const ccomplex_t IT_0287 = IT_0280 + IT_0282 + IT_0284 + IT_0286;
    const ccomplex_t IT_0288 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0289 = IT_0110*IT_0288;
    const ccomplex_t IT_0290 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0291 = IT_0113*IT_0290;
    const ccomplex_t IT_0292 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0293 = IT_0110*IT_0292;
    const ccomplex_t IT_0294 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0295 = IT_0113*IT_0294;
    const ccomplex_t IT_0296 = -IT_0289 + -IT_0291 + -IT_0293 + -IT_0295;
    const ccomplex_t IT_0297 = IT_0287 + IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = (-0.5)*IT_0298;
    const ccomplex_t IT_0300 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0301 = IT_0019*IT_0300;
    const ccomplex_t IT_0302 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0303 = IT_0019*IT_0302;
    const ccomplex_t IT_0304 = 1.4142135623731*IT_0303;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*(IT_0301 + (-0.5)*IT_0304);
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = IT_0125*IT_0306;
    const ccomplex_t IT_0308 = (0 + _Complex_I*1)*IT_0307;
    const ccomplex_t IT_0309 = IT_0299*IT_0308;
    const ccomplex_t IT_0310 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0311 = IT_0099*IT_0310;
    const ccomplex_t IT_0312 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0313 = IT_0102*IT_0312;
    const ccomplex_t IT_0314 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0315 = IT_0099*IT_0314;
    const ccomplex_t IT_0316 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0317 = IT_0102*IT_0316;
    const ccomplex_t IT_0318 = IT_0311 + IT_0313 + IT_0315 + IT_0317;
    const ccomplex_t IT_0319 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0320 = IT_0110*IT_0319;
    const ccomplex_t IT_0321 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0322 = IT_0113*IT_0321;
    const ccomplex_t IT_0323 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0324 = IT_0110*IT_0323;
    const ccomplex_t IT_0325 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0326 = IT_0113*IT_0325;
    const ccomplex_t IT_0327 = -IT_0320 + -IT_0322 + -IT_0324 + -IT_0326;
    const ccomplex_t IT_0328 = IT_0318 + IT_0327;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*IT_0328;
    const ccomplex_t IT_0330 = (-0.5)*IT_0329;
    const ccomplex_t IT_0331 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0332 = IT_0019*IT_0331;
    const ccomplex_t IT_0333 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0334 = IT_0019*IT_0333;
    const ccomplex_t IT_0335 = 1.4142135623731*IT_0334;
    const ccomplex_t IT_0336 = (0 + _Complex_I*1)*(IT_0332 + (-0.5)*IT_0335);
    const ccomplex_t IT_0337 = -IT_0336;
    const ccomplex_t IT_0338 = IT_0164*IT_0337;
    const ccomplex_t IT_0339 = (0 + _Complex_I*1)*IT_0338;
    const ccomplex_t IT_0340 = IT_0330*IT_0339;
    const ccomplex_t IT_0341 = IT_0190*IT_0263;
    const ccomplex_t IT_0342 = 2*IT_0341;
    const ccomplex_t IT_0343 = IT_0194*IT_0266;
    const ccomplex_t IT_0344 = 1.4142135623731*IT_0343;
    const ccomplex_t IT_0345 = IT_0190*IT_0269;
    const ccomplex_t IT_0346 = 1.4142135623731*IT_0345;
    const ccomplex_t IT_0347 = IT_0344 + IT_0346;
    const ccomplex_t IT_0348 = IT_0342 + IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = 0.5*IT_0349;
    const ccomplex_t IT_0351 = IT_0243*IT_0350;
    const ccomplex_t IT_0352 = (0 + _Complex_I*1)*IT_0351;
    const ccomplex_t IT_0353 = IT_0242*IT_0352;
    const ccomplex_t IT_0354 = N_B3*N_u3*e_em;
    const ccomplex_t IT_0355 = IT_0099*IT_0354;
    const ccomplex_t IT_0356 = N_B3*N_d3*e_em;
    const ccomplex_t IT_0357 = IT_0102*IT_0356;
    const ccomplex_t IT_0358 = N_W3*N_u3*e_em;
    const ccomplex_t IT_0359 = IT_0110*IT_0358;
    const ccomplex_t IT_0360 = N_W3*N_d3*e_em;
    const ccomplex_t IT_0361 = IT_0113*IT_0360;
    const ccomplex_t IT_0362 = (0 + _Complex_I*1)*(IT_0355 + IT_0357 + 
      -IT_0359 + -IT_0361);
    const ccomplex_t IT_0363 = -IT_0362;
    const ccomplex_t IT_0364 = IT_0069*IT_0178;
    const ccomplex_t IT_0365 = (0 + _Complex_I*1)*IT_0364;
    const ccomplex_t IT_0366 = IT_0363*IT_0365;
    const ccomplex_t IT_0367 = N_B4*N_u3*e_em;
    const ccomplex_t IT_0368 = IT_0099*IT_0367;
    const ccomplex_t IT_0369 = N_B4*N_d3*e_em;
    const ccomplex_t IT_0370 = IT_0102*IT_0369;
    const ccomplex_t IT_0371 = N_B3*N_u4*e_em;
    const ccomplex_t IT_0372 = IT_0099*IT_0371;
    const ccomplex_t IT_0373 = N_B3*N_d4*e_em;
    const ccomplex_t IT_0374 = IT_0102*IT_0373;
    const ccomplex_t IT_0375 = IT_0368 + IT_0370 + IT_0372 + IT_0374;
    const ccomplex_t IT_0376 = N_W4*N_u3*e_em;
    const ccomplex_t IT_0377 = IT_0110*IT_0376;
    const ccomplex_t IT_0378 = N_W4*N_d3*e_em;
    const ccomplex_t IT_0379 = IT_0113*IT_0378;
    const ccomplex_t IT_0380 = N_W3*N_u4*e_em;
    const ccomplex_t IT_0381 = IT_0110*IT_0380;
    const ccomplex_t IT_0382 = N_W3*N_d4*e_em;
    const ccomplex_t IT_0383 = IT_0113*IT_0382;
    const ccomplex_t IT_0384 = -IT_0377 + -IT_0379 + -IT_0381 + -IT_0383;
    const ccomplex_t IT_0385 = IT_0375 + IT_0384;
    const ccomplex_t IT_0386 = (0 + _Complex_I*1)*IT_0385;
    const ccomplex_t IT_0387 = (-0.5)*IT_0386;
    const ccomplex_t IT_0388 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0389 = IT_0019*IT_0388;
    const ccomplex_t IT_0390 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0391 = IT_0019*IT_0390;
    const ccomplex_t IT_0392 = 1.4142135623731*IT_0391;
    const ccomplex_t IT_0393 = (0 + _Complex_I*1)*(IT_0389 + (-0.5)*IT_0392);
    const ccomplex_t IT_0394 = -IT_0393;
    const ccomplex_t IT_0395 = IT_0237*IT_0394;
    const ccomplex_t IT_0396 = (0 + _Complex_I*1)*IT_0395;
    const ccomplex_t IT_0397 = IT_0387*IT_0396;
    const ccomplex_t IT_0398 = IT_0278 + (-2)*IT_0309 + (-2)*IT_0340 + IT_0353
       + (-2)*IT_0366 + (-2)*IT_0397;
    const ccomplex_t IT_0399 = (-2)*IT_0033;
    const ccomplex_t IT_0400 = s_12*s_24;
    const ccomplex_t IT_0401 = s_14*IT_0124;
    const ccomplex_t IT_0402 = -IT_0401;
    const ccomplex_t IT_0403 = IT_0400 + IT_0402;
    const ccomplex_t IT_0404 = (-2)*IT_0072 + (-2)*IT_0090 + -IT_0208 + 
      -IT_0258;
    const ccomplex_t IT_0405 = (-2)*IT_0033 + IT_0278 + IT_0353;
    const ccomplex_t IT_0406 = 2*IT_0033;
    const ccomplex_t IT_0407 = 2*IT_0053;
    const ccomplex_t IT_0408 = m_C_1*m_N_3;
    const ccomplex_t IT_0409 = IT_0000*IT_0408;
    const ccomplex_t IT_0410 = (-3)*IT_0409;
    const ccomplex_t IT_0411 = (-2)*conj(IT_0033);
    const ccomplex_t IT_0412 = (0 + _Complex_I*-1)*m_W*e_em*IT_0019*(IT_0015
      *IT_0186 + -IT_0012*IT_0187);
    const ccomplex_t IT_0413 = IT_0030*IT_0412;
    const ccomplex_t IT_0414 = cpow((-2)*s_12 + m_W*(m_W + (0 + _Complex_I*-1)
      *Gamma_W) + -IT_0050 + -IT_0124 + -reg_prop, -1);
    const ccomplex_t IT_0415 = IT_0413*IT_0414;
    const ccomplex_t IT_0416 = (0 + _Complex_I*1)*IT_0415;
    const ccomplex_t IT_0417 = m_N_1*IT_0299;
    const ccomplex_t IT_0418 = IT_0134*IT_0417;
    const ccomplex_t IT_0419 = m_N_2*IT_0330;
    const ccomplex_t IT_0420 = IT_0166*IT_0419;
    const ccomplex_t IT_0421 = IT_0030*IT_0062;
    const ccomplex_t IT_0422 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0423 = IT_0011*IT_0421*IT_0422;
    const ccomplex_t IT_0424 = IT_0048*IT_0080;
    const ccomplex_t IT_0425 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0426 = IT_0051*IT_0424*IT_0425;
    const ccomplex_t IT_0427 = m_N_3*IT_0363;
    const ccomplex_t IT_0428 = IT_0180*IT_0427;
    const ccomplex_t IT_0429 = -IT_0416 + -IT_0418 + -IT_0420 + IT_0423 +
       IT_0426 + (-0.5)*IT_0428;
    const ccomplex_t IT_0430 = m_C_1*IT_0387;
    const ccomplex_t IT_0431 = IT_0396*IT_0430;
    const ccomplex_t IT_0432 = m_C_1*IT_0299;
    const ccomplex_t IT_0433 = IT_0308*IT_0432;
    const ccomplex_t IT_0434 = m_C_1*IT_0330;
    const ccomplex_t IT_0435 = IT_0339*IT_0434;
    const ccomplex_t IT_0436 = m_C_1*IT_0363;
    const ccomplex_t IT_0437 = IT_0365*IT_0436;
    const ccomplex_t IT_0438 = m_N_3*IT_0072;
    const ccomplex_t IT_0439 = m_N_3*IT_0090;
    const ccomplex_t IT_0440 = m_N_4*IT_0387;
    const ccomplex_t IT_0441 = IT_0239*IT_0440;
    const ccomplex_t IT_0442 = (-0.5)*IT_0428 + IT_0431 + IT_0433 + IT_0435 +
       IT_0437 + -IT_0438 + -IT_0439 + -IT_0441;
    const ccomplex_t IT_0443 = conj(IT_0429) + conj(IT_0442);
    const ccomplex_t IT_0444 = s_24*IT_0000;
    const ccomplex_t IT_0445 = IT_0002*IT_0444;
    const ccomplex_t IT_0446 = 0.5*IT_0445;
    const ccomplex_t IT_0447 = s_24 + IT_0446;
    const ccomplex_t IT_0448 = m_C_1*IT_0447;
    const ccomplex_t IT_0449 = 2*IT_0448;
    const ccomplex_t IT_0450 = s_14*IT_0000;
    const ccomplex_t IT_0451 = IT_0002*IT_0450;
    const ccomplex_t IT_0452 = (-0.25)*IT_0451;
    const ccomplex_t IT_0453 = s_14 + IT_0452;
    const ccomplex_t IT_0454 = m_N_3*IT_0453;
    const ccomplex_t IT_0455 = (-4)*IT_0454;
    const ccomplex_t IT_0456 = IT_0261*IT_0455;
    const ccomplex_t IT_0457 = 2*conj(IT_0072);
    const ccomplex_t IT_0458 = 2*conj(IT_0090);
    const ccomplex_t IT_0459 = IT_0457 + IT_0458;
    const ccomplex_t IT_0460 = 2*IT_0096;
    const ccomplex_t IT_0461 = conj(IT_0398)*IT_0460;
    const ccomplex_t IT_0462 = 2*IT_0403;
    const ccomplex_t IT_0463 = conj(IT_0407)*IT_0462;
    const ccomplex_t IT_0464 = 0.5*IT_0181;
    const ccomplex_t IT_0465 = conj(IT_0240) + conj(IT_0464);
    const ccomplex_t IT_0466 = 3*IT_0409;
    const ccomplex_t IT_0467 = 2*IT_0466;
    const ccomplex_t IT_0468 = m_C_1*IT_0177;
    const ccomplex_t IT_0469 = IT_0180*IT_0468;
    const ccomplex_t IT_0470 = m_C_1*IT_0123;
    const ccomplex_t IT_0471 = IT_0134*IT_0470;
    const ccomplex_t IT_0472 = m_N_3*IT_0053;
    const ccomplex_t IT_0473 = m_C_1*IT_0156;
    const ccomplex_t IT_0474 = IT_0166*IT_0473;
    const ccomplex_t IT_0475 = 0.5*IT_0469 + IT_0471 + -IT_0472 + IT_0474;
    const ccomplex_t IT_0476 = IT_0023*IT_0069;
    const ccomplex_t IT_0477 = IT_0011*IT_0422*IT_0476;
    const ccomplex_t IT_0478 = IT_0041*IT_0087;
    const ccomplex_t IT_0479 = IT_0051*IT_0425*IT_0478;
    const ccomplex_t IT_0480 = m_N_3*IT_0033;
    const ccomplex_t IT_0481 = IT_0069*IT_0412;
    const ccomplex_t IT_0482 = IT_0414*IT_0481;
    const ccomplex_t IT_0483 = (0 + _Complex_I*1)*IT_0482;
    const ccomplex_t IT_0484 = m_N_1*IT_0123;
    const ccomplex_t IT_0485 = IT_0308*IT_0484;
    const ccomplex_t IT_0486 = m_N_2*IT_0156;
    const ccomplex_t IT_0487 = IT_0339*IT_0486;
    const ccomplex_t IT_0488 = m_N_3*IT_0177;
    const ccomplex_t IT_0489 = IT_0365*IT_0488;
    const ccomplex_t IT_0490 = m_N_4*IT_0229;
    const ccomplex_t IT_0491 = IT_0396*IT_0490;
    const ccomplex_t IT_0492 = m_C_1*IT_0229;
    const ccomplex_t IT_0493 = IT_0239*IT_0492;
    const ccomplex_t IT_0494 = 0.5*IT_0469 + IT_0477 + IT_0479 + -IT_0480 + 
      -IT_0483 + -IT_0485 + -IT_0487 + -IT_0489 + -IT_0491 + IT_0493;
    const ccomplex_t IT_0495 = conj(IT_0475) + conj(IT_0494);
    const ccomplex_t IT_0496 = 2*IT_0455;
    const ccomplex_t IT_0497 = IT_0309 + IT_0340 + IT_0366 + IT_0397;
    const ccomplex_t IT_0498 = (-0.25)*IT_0005;
    const ccomplex_t IT_0499 = IT_0007 + IT_0498;
    const ccomplex_t IT_0500 = 8*IT_0499;
    const ccomplex_t IT_0501 = (-4)*IT_0009;
    const ccomplex_t IT_0502 = (-2)*IT_0403;
    const ccomplex_t IT_0503 = conj(IT_0405)*IT_0502;
    const ccomplex_t IT_0504 = -IT_0135 + -IT_0167 + (-0.5)*IT_0181;
    const ccomplex_t IT_0505 = (-2)*IT_0466;
    const ccomplex_t IT_0506 = (-2)*IT_0449;
    const ccomplex_t IT_0507 = IT_0009*conj(IT_0033) + 0.25*IT_0410*IT_0459 +
       0.25*IT_0461 + 0.25*IT_0463 + 0.25*IT_0465*IT_0467 + 0.25*IT_0495*IT_0496
       + 0.25*conj(IT_0497)*IT_0500 + 0.25*conj(IT_0054)*IT_0501 + 0.25*IT_0503 
      + 0.25*conj(IT_0504)*IT_0505 + 0.25*IT_0443*IT_0506;
    const ccomplex_t IT_0508 = IT_0183*IT_0449;
    const ccomplex_t IT_0509 = IT_0184*IT_0449;
    const ccomplex_t IT_0510 = pow(s_24, 2);
    const ccomplex_t IT_0511 = IT_0002*IT_0510;
    const ccomplex_t IT_0512 = -IT_0124;
    const ccomplex_t IT_0513 = IT_0511 + IT_0512;
    const ccomplex_t IT_0514 = m_C_1*IT_0513;
    const ccomplex_t IT_0515 = 2*IT_0514;
    const ccomplex_t IT_0516 = IT_0404*IT_0515;
    const ccomplex_t IT_0517 = IT_0001 + IT_0008;
    const ccomplex_t IT_0518 = IT_0002*IT_0517;
    const ccomplex_t IT_0519 = (-0.5)*IT_0518;
    const ccomplex_t IT_0520 = s_12 + IT_0519;
    const ccomplex_t IT_0521 = 4*IT_0520;
    const ccomplex_t IT_0522 = IT_0475*IT_0521;
    const ccomplex_t IT_0523 = s_14*s_24*IT_0002;
    const ccomplex_t IT_0524 = -IT_0523;
    const ccomplex_t IT_0525 = s_12 + IT_0524;
    const ccomplex_t IT_0526 = m_C_1*IT_0525;
    const ccomplex_t IT_0527 = (-2)*IT_0526;
    const ccomplex_t IT_0528 = IT_0259*IT_0527;
    const ccomplex_t IT_0529 = pow(s_14, 2);
    const ccomplex_t IT_0530 = IT_0002*IT_0529;
    const ccomplex_t IT_0531 = -IT_0050;
    const ccomplex_t IT_0532 = IT_0530 + IT_0531;
    const ccomplex_t IT_0533 = s_12*IT_0532;
    const ccomplex_t IT_0534 = 2*IT_0533;
    const ccomplex_t IT_0535 = s_12*IT_0525;
    const ccomplex_t IT_0536 = 2*IT_0535;
    const ccomplex_t IT_0537 = 2*IT_0526;
    const ccomplex_t IT_0538 = IT_0408*IT_0532;
    const ccomplex_t IT_0539 = 2*IT_0538;
    const ccomplex_t IT_0540 = IT_0408*IT_0525;
    const ccomplex_t IT_0541 = 2*IT_0540;
    const ccomplex_t IT_0542 = m_N_3*IT_0532;
    const ccomplex_t IT_0543 = 2*IT_0542;
    const ccomplex_t IT_0544 = (-2)*IT_0096;
    const ccomplex_t IT_0545 = (-2)*IT_0533;
    const ccomplex_t IT_0546 = (-2)*IT_0540;
    const ccomplex_t IT_0547 = IT_0096*IT_0459 + IT_0460*conj(IT_0504) + conj
      (IT_0182)*IT_0534 + conj(IT_0404)*IT_0536 + IT_0495*IT_0537 + conj(IT_0398
      )*IT_0539 + conj(IT_0407)*IT_0541 + IT_0443*IT_0543 + IT_0465*IT_0544 +
       conj(IT_0259)*IT_0545 + conj(IT_0405)*IT_0546;
    const ccomplex_t IT_0548 = (-2)*conj(IT_0072);
    const ccomplex_t IT_0549 = (-2)*conj(IT_0090);
    const ccomplex_t IT_0550 = IT_0548 + IT_0549;
    const ccomplex_t IT_0551 = (-2)*IT_0535;
    const ccomplex_t IT_0552 = (-2)*IT_0538;
    const ccomplex_t IT_0553 = (-2)*IT_0542;
    const ccomplex_t IT_0554 = conj(IT_0054) + conj(IT_0497);
    const ccomplex_t IT_0555 = IT_0398*(IT_0096*IT_0411 + IT_0443*IT_0527 +
       conj(IT_0398)*IT_0534 + conj(IT_0407)*IT_0536 + conj(IT_0182)*IT_0539 +
       conj(IT_0404)*IT_0541 + conj(IT_0405)*IT_0551 + conj(IT_0259)*IT_0552 +
       IT_0495*IT_0553 + IT_0460*IT_0554);
    const ccomplex_t IT_0556 = s_12*IT_0513;
    const ccomplex_t IT_0557 = 2*IT_0556;
    const ccomplex_t IT_0558 = IT_0408*IT_0513;
    const ccomplex_t IT_0559 = 2*IT_0558;
    const ccomplex_t IT_0560 = m_N_3*IT_0525;
    const ccomplex_t IT_0561 = 2*IT_0560;
    const ccomplex_t IT_0562 = (-2)*IT_0558;
    const ccomplex_t IT_0563 = 2*conj(IT_0033);
    const ccomplex_t IT_0564 = (-2)*IT_0556;
    const ccomplex_t IT_0565 = IT_0405*(IT_0443*IT_0515 + conj(IT_0259)
      *IT_0541 + conj(IT_0182)*IT_0546 + conj(IT_0398)*IT_0551 + IT_0502*IT_0554
       + conj(IT_0405)*IT_0557 + IT_0495*IT_0561 + conj(IT_0404)*IT_0562 +
       IT_0403*IT_0563 + conj(IT_0407)*IT_0564);
    const ccomplex_t IT_0566 = (-2)*IT_0514;
    const ccomplex_t IT_0567 = (-2)*IT_0560;
    const ccomplex_t IT_0568 = IT_0407*(IT_0403*IT_0411 + conj(IT_0398)
      *IT_0536 + conj(IT_0182)*IT_0541 + conj(IT_0259)*IT_0546 + IT_0462*IT_0554
       + conj(IT_0407)*IT_0557 + conj(IT_0404)*IT_0559 + conj(IT_0405)*IT_0564 +
       IT_0443*IT_0566 + IT_0495*IT_0567);
    const ccomplex_t IT_0569 = conj(IT_0259)*IT_0460;
    const ccomplex_t IT_0570 = conj(IT_0054)*IT_0467;
    const ccomplex_t IT_0571 = (-0.25)*IT_0445;
    const ccomplex_t IT_0572 = s_24 + IT_0571;
    const ccomplex_t IT_0573 = m_C_1*IT_0572;
    const ccomplex_t IT_0574 = (-4)*IT_0573;
    const ccomplex_t IT_0575 = 2*IT_0574;
    const ccomplex_t IT_0576 = conj(IT_0475)*IT_0575;
    const ccomplex_t IT_0577 = conj(IT_0494)*IT_0575;
    const ccomplex_t IT_0578 = 4*IT_0009;
    const ccomplex_t IT_0579 = conj(IT_0504)*IT_0578;
    const ccomplex_t IT_0580 = (-8)*conj(IT_0072)*IT_0499;
    const ccomplex_t IT_0581 = (-8)*conj(IT_0090)*IT_0499;
    const ccomplex_t IT_0582 = IT_0411*IT_0466;
    const ccomplex_t IT_0583 = conj(IT_0182)*IT_0544;
    const ccomplex_t IT_0584 = conj(IT_0404)*IT_0502;
    const ccomplex_t IT_0585 = (-2)*IT_0410;
    const ccomplex_t IT_0586 = conj(IT_0497)*IT_0585;
    const ccomplex_t IT_0587 = 0.5*IT_0451;
    const ccomplex_t IT_0588 = s_14 + IT_0587;
    const ccomplex_t IT_0589 = m_N_3*IT_0588;
    const ccomplex_t IT_0590 = 2*IT_0589;
    const ccomplex_t IT_0591 = (-2)*IT_0590;
    const ccomplex_t IT_0592 = conj(IT_0429)*IT_0591;
    const ccomplex_t IT_0593 = conj(IT_0442)*IT_0591;
    const ccomplex_t IT_0594 = 0.125*IT_0466;
    const ccomplex_t IT_0595 = 2*IT_0410;
    const ccomplex_t IT_0596 = 2*IT_0590;
    const ccomplex_t IT_0597 = (-2)*IT_0574;
    const ccomplex_t IT_0598 = 0.125*conj(IT_0182)*IT_0460 + 0.125*conj
      (IT_0404)*IT_0462 + IT_0092*IT_0499 + 0.125*IT_0501*conj(IT_0504) + 0.125
      *conj(IT_0054)*IT_0505 + 0.125*conj(IT_0259)*IT_0544 + 0.125*IT_0465
      *IT_0578 + IT_0563*IT_0594 + 0.125*conj(IT_0497)*IT_0595 + 0.125*IT_0443
      *IT_0596 + 0.125*IT_0495*IT_0597;
    const ccomplex_t IT_0599 = (-0.125)*IT_0466;
    const ccomplex_t IT_0600 = (-0.125)*IT_0461 + (-0.125)*IT_0463 + conj
      (IT_0033)*IT_0499 + (-0.125)*conj(IT_0054)*IT_0500 + (-0.125)*conj(IT_0497
      )*IT_0501 + (-0.125)*IT_0503 + (-0.125)*IT_0443*IT_0575 + (-0.125)*IT_0465
      *IT_0585 + (-0.125)*IT_0495*IT_0591 + (-0.125)*conj(IT_0504)*IT_0595 +
       IT_0550*IT_0599;
    const ccomplex_t IT_0601 = IT_0259*IT_0460;
    const ccomplex_t IT_0602 = IT_0475*IT_0575;
    const ccomplex_t IT_0603 = IT_0494*IT_0575;
    const ccomplex_t IT_0604 = (-8)*IT_0090*IT_0499;
    const ccomplex_t IT_0605 = IT_0240*IT_0501;
    const ccomplex_t IT_0606 = IT_0404*IT_0502;
    const ccomplex_t IT_0607 = IT_0442*IT_0591;
    const ccomplex_t IT_0608 = IT_0449*IT_0563;
    const ccomplex_t IT_0609 = conj(IT_0405)*IT_0515;
    const ccomplex_t IT_0610 = IT_0455*IT_0548;
    const ccomplex_t IT_0611 = IT_0455*IT_0549;
    const ccomplex_t IT_0612 = conj(IT_0398)*IT_0527;
    const ccomplex_t IT_0613 = conj(IT_0407)*IT_0566;
    const ccomplex_t IT_0614 = conj(IT_0182)*IT_0543;
    const ccomplex_t IT_0615 = conj(IT_0259)*IT_0553;
    const ccomplex_t IT_0616 = conj(IT_0404)*IT_0561;
    const ccomplex_t IT_0617 = 6*IT_0408;
    const ccomplex_t IT_0618 = IT_0495*IT_0617;
    const ccomplex_t IT_0619 = IT_0443*IT_0521;
    const ccomplex_t IT_0620 = conj(IT_0054)*IT_0506;
    const ccomplex_t IT_0621 = conj(IT_0504)*IT_0596;
    const ccomplex_t IT_0622 = conj(IT_0497)*IT_0575;
    const ccomplex_t IT_0623 = IT_0449*IT_0457;
    const ccomplex_t IT_0624 = IT_0449*IT_0458;
    const ccomplex_t IT_0625 = conj(IT_0182)*IT_0537;
    const ccomplex_t IT_0626 = conj(IT_0404)*IT_0515;
    const ccomplex_t IT_0627 = conj(IT_0405)*IT_0561;
    const ccomplex_t IT_0628 = conj(IT_0429)*IT_0617;
    const ccomplex_t IT_0629 = IT_0411*IT_0455;
    const ccomplex_t IT_0630 = conj(IT_0259)*IT_0527;
    const ccomplex_t IT_0631 = conj(IT_0398)*IT_0553;
    const ccomplex_t IT_0632 = conj(IT_0407)*IT_0567;
    const ccomplex_t IT_0633 = conj(IT_0442)*IT_0617;
    const ccomplex_t IT_0634 = conj(IT_0054)*IT_0496;
    const ccomplex_t IT_0635 = conj(IT_0504)*IT_0597;
    const ccomplex_t IT_0636 = conj(IT_0497)*IT_0591;
    const ccomplex_t IT_0637 = (-4)*IT_0009*(IT_0033*(conj(IT_0033) + -conj
      (IT_0054)) + IT_0091*IT_0092) + IT_0096*(conj(IT_0182)*IT_0185 + conj
      (IT_0259)*IT_0262 + conj(IT_0398)*IT_0399) + IT_0403*(IT_0185*conj(IT_0404
      ) + conj(IT_0405)*IT_0406 + IT_0399*conj(IT_0407)) + IT_0410*(conj(IT_0054
      )*IT_0185 + IT_0092*IT_0399 + IT_0091*IT_0411) + IT_0443*(IT_0406*IT_0449 
      + IT_0260*IT_0455 + IT_0456) + 4*IT_0054*IT_0507 + IT_0495*(IT_0399
      *IT_0455 + IT_0508 + IT_0509 + IT_0516 + IT_0522 + IT_0528) + IT_0182
      *IT_0547 + IT_0259*(conj(IT_0259)*IT_0534 + conj(IT_0405)*IT_0541 + conj
      (IT_0504)*IT_0544 + conj(IT_0182)*IT_0545 + conj(IT_0407)*IT_0546 +
       IT_0096*IT_0550 + conj(IT_0404)*IT_0551 + conj(IT_0398)*IT_0552 + IT_0443
      *IT_0553) + IT_0555 + IT_0404*(IT_0403*IT_0459 + IT_0462*conj(IT_0504) +
       conj(IT_0182)*IT_0536 + conj(IT_0398)*IT_0541 + conj(IT_0259)*IT_0551 +
       conj(IT_0404)*IT_0557 + conj(IT_0407)*IT_0559 + IT_0443*IT_0561 + conj
      (IT_0405)*IT_0562) + IT_0565 + IT_0568 + IT_0240*(IT_0569 + IT_0570 +
       IT_0576 + IT_0577 + IT_0579 + IT_0580 + IT_0581 + IT_0582 + IT_0583 +
       IT_0584 + IT_0586 + IT_0592 + IT_0593) + IT_0464*(IT_0465*IT_0501 +
       IT_0569 + IT_0570 + IT_0576 + IT_0577 + IT_0579 + IT_0580 + IT_0581 +
       IT_0582 + IT_0583 + IT_0584 + IT_0586 + IT_0592 + IT_0593) + 8*conj
      (IT_0504)*(IT_0091*IT_0499 + IT_0406*IT_0594) + 8*IT_0504*IT_0598 + (-8)
      *conj(IT_0497)*(IT_0033*IT_0499 + IT_0262*IT_0599) + (-8)*IT_0497*IT_0600 
      + (-8)*IT_0465*(IT_0072*IT_0499 + IT_0399*IT_0599 + (-0.125)*IT_0601 + (
      -0.125)*IT_0602 + (-0.125)*IT_0603 + (-0.125)*IT_0604 + (-0.125)*IT_0605 +
       (-0.125)*IT_0606 + (-0.125)*IT_0607) + IT_0442*(IT_0608 + IT_0609 +
       IT_0610 + IT_0611 + IT_0612 + IT_0613 + IT_0614 + IT_0615 + IT_0616 +
       IT_0618 + IT_0619 + IT_0620 + IT_0621 + IT_0622) + IT_0429*(IT_0465
      *IT_0591 + IT_0608 + IT_0609 + IT_0610 + IT_0611 + IT_0612 + IT_0613 +
       IT_0614 + IT_0615 + IT_0616 + IT_0618 + IT_0619 + IT_0620 + IT_0621 +
       IT_0622) + IT_0475*(IT_0623 + IT_0624 + IT_0625 + IT_0626 + IT_0627 +
       IT_0628 + IT_0629 + IT_0630 + IT_0631 + IT_0632 + IT_0633 + IT_0634 +
       IT_0635 + IT_0636) + IT_0494*(IT_0495*IT_0521 + IT_0623 + IT_0624 +
       IT_0625 + IT_0626 + IT_0627 + IT_0628 + IT_0629 + IT_0630 + IT_0631 +
       IT_0632 + IT_0633 + IT_0634 + IT_0635 + IT_0636);
    return create_ccomplex_return(IT_0637);
}

