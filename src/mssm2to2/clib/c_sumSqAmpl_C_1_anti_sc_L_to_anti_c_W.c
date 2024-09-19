#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_sc_L_to_anti_c_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_sc_L_to_anti_c_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_b = param->m_b;
    const creal_t m_c = param->m_c;
    const creal_t m_d = param->m_d;
    const creal_t m_s = param->m_s;
    const creal_t V_cb = param->V_cb;
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
    const creal_t m_sc_L = param->m_sc_L;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t V_cd = param->V_cd;
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0001 = pow(m_c, 2);
    const ccomplex_t IT_0002 = pow(m_sc_L, 2);
    const ccomplex_t IT_0003 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0004 = pow(m_W, -1);
    const ccomplex_t IT_0005 = sin(beta);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = cpow(IT_0007, -1);
    const ccomplex_t IT_0009 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u3)
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0010 = (-0.5)*IT_0009;
    const ccomplex_t IT_0011 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0012 = IT_0008*IT_0011;
    const ccomplex_t IT_0013 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0014 = IT_0008*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*(IT_0012 + (-0.5)*IT_0015);
    const ccomplex_t IT_0017 = -IT_0016;
    const ccomplex_t IT_0018 = IT_0010*IT_0017;
    const ccomplex_t IT_0019 = IT_0000*IT_0003*IT_0018;
    const ccomplex_t IT_0020 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u1)
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0022 = (-0.5)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0024 = IT_0008*IT_0023;
    const ccomplex_t IT_0025 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0026 = IT_0008*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + (-0.5)*IT_0027);
    const ccomplex_t IT_0029 = -IT_0028;
    const ccomplex_t IT_0030 = IT_0022*IT_0029;
    const ccomplex_t IT_0031 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0032 = IT_0020*IT_0030*IT_0031;
    const ccomplex_t IT_0033 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u2)
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0035 = (-0.5)*IT_0034;
    const ccomplex_t IT_0036 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0037 = IT_0008*IT_0036;
    const ccomplex_t IT_0038 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0039 = IT_0008*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + (-0.5)*IT_0040);
    const ccomplex_t IT_0042 = -IT_0041;
    const ccomplex_t IT_0043 = IT_0035*IT_0042;
    const ccomplex_t IT_0044 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0045 = IT_0033*IT_0043*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0047 = cpow((-2)*s_23 + IT_0001 + IT_0002 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0048 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0049 = IT_0008*IT_0048;
    const ccomplex_t IT_0050 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0051 = IT_0008*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.5)*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1.4142135623731)*m_c*conj(N_u4)
      *e_em*IT_0004*IT_0006*IT_0008;
    const ccomplex_t IT_0056 = (-0.5)*IT_0055;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = IT_0046*IT_0047*IT_0057;
    const ccomplex_t IT_0059 = -IT_0019 + -IT_0032 + -IT_0045 + -IT_0058;
    const ccomplex_t IT_0060 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0061 = IT_0008*IT_0060;
    const ccomplex_t IT_0062 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0063 = IT_0008*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + 0.5*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = IT_0035*IT_0066;
    const ccomplex_t IT_0068 = IT_0044*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = m_C_1*IT_0069;
    const ccomplex_t IT_0071 = pow(m_b, 2);
    const ccomplex_t IT_0072 = pow(m_C_1, 2);
    const ccomplex_t IT_0073 = cpow((-2)*s_12 + -IT_0002 + IT_0071 + -IT_0072 
      + -reg_prop, -1);
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*V_cb*e_em*V_Wp1*IT_0008;
    const ccomplex_t IT_0075 = -IT_0074;
    const ccomplex_t IT_0076 = (0 + _Complex_I*1.4142135623731)*V_cb*e_em
      *IT_0008;
    const ccomplex_t IT_0077 = 0.5*IT_0076;
    const ccomplex_t IT_0078 = IT_0075*IT_0077;
    const ccomplex_t IT_0079 = IT_0073*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = m_c*IT_0080;
    const ccomplex_t IT_0082 = pow(m_d, 2);
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + -IT_0002 + -IT_0072 + IT_0082 
      + -reg_prop, -1);
    const ccomplex_t IT_0084 = (0 + _Complex_I*1)*V_cd*e_em*V_Wp1*IT_0008;
    const ccomplex_t IT_0085 = -IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*conj(V_cd)
      *e_em*IT_0008;
    const ccomplex_t IT_0087 = 0.5*IT_0086;
    const ccomplex_t IT_0088 = IT_0085*IT_0087;
    const ccomplex_t IT_0089 = IT_0083*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = m_c*IT_0090;
    const ccomplex_t IT_0092 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0093 = IT_0008*IT_0092;
    const ccomplex_t IT_0094 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0095 = IT_0008*IT_0094;
    const ccomplex_t IT_0096 = 1.4142135623731*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*(IT_0093 + 0.5*IT_0096);
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = IT_0022*IT_0098;
    const ccomplex_t IT_0100 = IT_0031*IT_0099;
    const ccomplex_t IT_0101 = (0 + _Complex_I*1)*IT_0100;
    const ccomplex_t IT_0102 = m_C_1*IT_0101;
    const ccomplex_t IT_0103 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0104 = IT_0008*IT_0103;
    const ccomplex_t IT_0105 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0106 = IT_0008*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*(IT_0104 + 0.5*IT_0107);
    const ccomplex_t IT_0109 = -IT_0108;
    const ccomplex_t IT_0110 = IT_0010*IT_0109;
    const ccomplex_t IT_0111 = IT_0003*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = m_C_1*IT_0112;
    const ccomplex_t IT_0114 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0115 = IT_0008*IT_0114;
    const ccomplex_t IT_0116 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0117 = IT_0008*IT_0116;
    const ccomplex_t IT_0118 = 1.4142135623731*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*(IT_0115 + 0.5*IT_0118);
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = IT_0056*IT_0120;
    const ccomplex_t IT_0122 = IT_0047*IT_0121;
    const ccomplex_t IT_0123 = (0 + _Complex_I*1)*IT_0122;
    const ccomplex_t IT_0124 = m_C_1*IT_0123;
    const ccomplex_t IT_0125 = pow(m_s, 2);
    const ccomplex_t IT_0126 = cpow((-2)*s_12 + -IT_0002 + -IT_0072 + IT_0125 
      + -reg_prop, -1);
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*V_cs*e_em*V_Wp1*IT_0008;
    const ccomplex_t IT_0128 = -IT_0127;
    const ccomplex_t IT_0129 = (0 + _Complex_I*1.4142135623731)*conj(V_cs)
      *e_em*IT_0008;
    const ccomplex_t IT_0130 = 0.5*IT_0129;
    const ccomplex_t IT_0131 = IT_0128*IT_0130;
    const ccomplex_t IT_0132 = IT_0126*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = m_c*IT_0133;
    const ccomplex_t IT_0135 = IT_0070 + IT_0081 + IT_0091 + IT_0102 + IT_0113
       + IT_0124 + IT_0134;
    const ccomplex_t IT_0136 = pow(m_W, -2);
    const ccomplex_t IT_0137 = s_14*s_34*IT_0136;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = s_13 + IT_0138;
    const ccomplex_t IT_0140 = m_C_1*IT_0139;
    const ccomplex_t IT_0141 = 2*IT_0069 + 2*IT_0101 + 2*IT_0112 + 2*IT_0123;
    const ccomplex_t IT_0142 = pow(s_14, 2);
    const ccomplex_t IT_0143 = IT_0136*IT_0142;
    const ccomplex_t IT_0144 = -IT_0072;
    const ccomplex_t IT_0145 = IT_0143 + IT_0144;
    const ccomplex_t IT_0146 = m_c*IT_0145;
    const ccomplex_t IT_0147 = cos(theta_W);
    const ccomplex_t IT_0148 = cpow(IT_0147, -1);
    const ccomplex_t IT_0149 = N_B2*e_em;
    const ccomplex_t IT_0150 = IT_0148*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = N_W2*e_em;
    const ccomplex_t IT_0153 = IT_0008*IT_0152;
    const ccomplex_t IT_0154 = 1.4142135623731*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1)*(IT_0151 + 3*IT_0154);
    const ccomplex_t IT_0156 = 0.166666666666667*IT_0155;
    const ccomplex_t IT_0157 = IT_0042*IT_0156;
    const ccomplex_t IT_0158 = IT_0044*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0160 = N_B3*e_em;
    const ccomplex_t IT_0161 = IT_0148*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = N_W3*e_em;
    const ccomplex_t IT_0164 = IT_0008*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = (0 + _Complex_I*1)*(IT_0162 + 3*IT_0165);
    const ccomplex_t IT_0167 = 0.166666666666667*IT_0166;
    const ccomplex_t IT_0168 = IT_0017*IT_0167;
    const ccomplex_t IT_0169 = IT_0003*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*IT_0169;
    const ccomplex_t IT_0171 = N_B1*e_em;
    const ccomplex_t IT_0172 = IT_0148*IT_0171;
    const ccomplex_t IT_0173 = 1.4142135623731*IT_0172;
    const ccomplex_t IT_0174 = N_W1*e_em;
    const ccomplex_t IT_0175 = IT_0008*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0173 + 3*IT_0176);
    const ccomplex_t IT_0178 = 0.166666666666667*IT_0177;
    const ccomplex_t IT_0179 = IT_0029*IT_0178;
    const ccomplex_t IT_0180 = IT_0031*IT_0179;
    const ccomplex_t IT_0181 = (0 + _Complex_I*1)*IT_0180;
    const ccomplex_t IT_0182 = N_B4*e_em;
    const ccomplex_t IT_0183 = IT_0148*IT_0182;
    const ccomplex_t IT_0184 = 1.4142135623731*IT_0183;
    const ccomplex_t IT_0185 = N_W4*e_em;
    const ccomplex_t IT_0186 = IT_0008*IT_0185;
    const ccomplex_t IT_0187 = 1.4142135623731*IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1)*(IT_0184 + 3*IT_0187);
    const ccomplex_t IT_0189 = 0.166666666666667*IT_0188;
    const ccomplex_t IT_0190 = IT_0054*IT_0189;
    const ccomplex_t IT_0191 = IT_0047*IT_0190;
    const ccomplex_t IT_0192 = (0 + _Complex_I*1)*IT_0191;
    const ccomplex_t IT_0193 = 2*IT_0159 + 2*IT_0170 + 2*IT_0181 + 2*IT_0192;
    const ccomplex_t IT_0194 = pow(m_W, 2);
    const ccomplex_t IT_0195 = s_14*IT_0194;
    const ccomplex_t IT_0196 = IT_0136*IT_0195;
    const ccomplex_t IT_0197 = (-0.25)*IT_0196;
    const ccomplex_t IT_0198 = s_14 + IT_0197;
    const ccomplex_t IT_0199 = m_c*IT_0198;
    const ccomplex_t IT_0200 = (-4)*IT_0199;
    const ccomplex_t IT_0201 = (-6)*IT_0200;
    const ccomplex_t IT_0202 = IT_0080 + IT_0090 + IT_0133;
    const ccomplex_t IT_0203 = 0.166666666666667*conj(IT_0202);
    const ccomplex_t IT_0204 = IT_0201*IT_0203;
    const ccomplex_t IT_0205 = s_34*IT_0194;
    const ccomplex_t IT_0206 = IT_0136*IT_0205;
    const ccomplex_t IT_0207 = (-0.25)*IT_0206;
    const ccomplex_t IT_0208 = s_34 + IT_0207;
    const ccomplex_t IT_0209 = m_C_1*IT_0208;
    const ccomplex_t IT_0210 = (-4)*IT_0209;
    const ccomplex_t IT_0211 = (-6)*IT_0210;
    const ccomplex_t IT_0212 = -IT_0069 + -IT_0101 + -IT_0112 + -IT_0123;
    const ccomplex_t IT_0213 = 0.166666666666667*conj(IT_0212);
    const ccomplex_t IT_0214 = IT_0211*IT_0213;
    const ccomplex_t IT_0215 = m_c*m_C_1;
    const ccomplex_t IT_0216 = (-18)*IT_0215;
    const ccomplex_t IT_0217 = IT_0159 + IT_0170 + IT_0181 + IT_0192;
    const ccomplex_t IT_0218 = m_C_1*IT_0217;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = IT_0098*IT_0178;
    const ccomplex_t IT_0221 = IT_0020*IT_0031*IT_0220;
    const ccomplex_t IT_0222 = IT_0066*IT_0156;
    const ccomplex_t IT_0223 = IT_0033*IT_0044*IT_0222;
    const ccomplex_t IT_0224 = IT_0109*IT_0167;
    const ccomplex_t IT_0225 = IT_0000*IT_0003*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0227 = cos(beta);
    const ccomplex_t IT_0228 = cpow(IT_0227, -1);
    const ccomplex_t IT_0229 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_cb*e_em*IT_0004*IT_0008*IT_0228;
    const ccomplex_t IT_0230 = (-0.5)*IT_0229;
    const ccomplex_t IT_0231 = IT_0077*IT_0230;
    const ccomplex_t IT_0232 = IT_0073*IT_0226*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*m_d;
    const ccomplex_t IT_0234 = (0 + _Complex_I*1.4142135623731)*m_d*conj(U_d1)
      *V_cd*e_em*IT_0004*IT_0008*IT_0228;
    const ccomplex_t IT_0235 = (-0.5)*IT_0234;
    const ccomplex_t IT_0236 = IT_0087*IT_0235;
    const ccomplex_t IT_0237 = IT_0083*IT_0233*IT_0236;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*m_s;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0004*IT_0008*IT_0228;
    const ccomplex_t IT_0240 = (-0.5)*IT_0239;
    const ccomplex_t IT_0241 = IT_0130*IT_0240;
    const ccomplex_t IT_0242 = IT_0126*IT_0238*IT_0241;
    const ccomplex_t IT_0243 = IT_0120*IT_0189;
    const ccomplex_t IT_0244 = IT_0046*IT_0047*IT_0243;
    const ccomplex_t IT_0245 = IT_0221 + IT_0223 + IT_0225 + IT_0232 + IT_0237
       + IT_0242 + IT_0244;
    const ccomplex_t IT_0246 = conj(IT_0219) + conj(IT_0245);
    const ccomplex_t IT_0247 = 0.166666666666667*IT_0246;
    const ccomplex_t IT_0248 = s_13*IT_0194;
    const ccomplex_t IT_0249 = s_14*s_34;
    const ccomplex_t IT_0250 = (-2)*IT_0249;
    const ccomplex_t IT_0251 = IT_0248 + IT_0250;
    const ccomplex_t IT_0252 = IT_0136*IT_0251;
    const ccomplex_t IT_0253 = (-6)*IT_0252;
    const ccomplex_t IT_0254 = 12*s_13;
    const ccomplex_t IT_0255 = IT_0253 + IT_0254;
    const ccomplex_t IT_0256 = conj(IT_0059) + conj(IT_0135);
    const ccomplex_t IT_0257 = 0.166666666666667*IT_0256;
    const ccomplex_t IT_0258 = 0.5*IT_0196;
    const ccomplex_t IT_0259 = s_14 + IT_0258;
    const ccomplex_t IT_0260 = m_c*IT_0259;
    const ccomplex_t IT_0261 = 2*IT_0260;
    const ccomplex_t IT_0262 = -IT_0159 + -IT_0170 + -IT_0181 + -IT_0192;
    const ccomplex_t IT_0263 = IT_0261*conj(IT_0262);
    const ccomplex_t IT_0264 = m_c*IT_0139;
    const ccomplex_t IT_0265 = 2*IT_0080 + 2*IT_0090 + 2*IT_0133;
    const ccomplex_t IT_0266 = IT_0264*conj(IT_0265);
    const ccomplex_t IT_0267 = IT_0140*conj(IT_0141) + IT_0146*conj(IT_0193) +
       IT_0204 + IT_0214 + IT_0216*IT_0247 + IT_0255*IT_0257 + IT_0263 + IT_0266;
    const ccomplex_t IT_0268 = pow(s_34, 2);
    const ccomplex_t IT_0269 = IT_0136*IT_0268;
    const ccomplex_t IT_0270 = -IT_0269;
    const ccomplex_t IT_0271 = IT_0001 + IT_0270;
    const ccomplex_t IT_0272 = s_13*IT_0271;
    const ccomplex_t IT_0273 = m_C_1*IT_0271;
    const ccomplex_t IT_0274 = s_13*IT_0139;
    const ccomplex_t IT_0275 = 6*IT_0274;
    const ccomplex_t IT_0276 = conj(IT_0202) + conj(IT_0262);
    const ccomplex_t IT_0277 = s_14*IT_0001;
    const ccomplex_t IT_0278 = s_13*s_34;
    const ccomplex_t IT_0279 = -IT_0278;
    const ccomplex_t IT_0280 = IT_0277 + IT_0279;
    const ccomplex_t IT_0281 = (-6)*IT_0280;
    const ccomplex_t IT_0282 = IT_0139*IT_0215;
    const ccomplex_t IT_0283 = (-6)*IT_0282;
    const ccomplex_t IT_0284 = IT_0219*IT_0255;
    const ccomplex_t IT_0285 = IT_0210*conj(IT_0262);
    const ccomplex_t IT_0286 = conj(IT_0265)*IT_0273;
    const ccomplex_t IT_0287 = IT_0285 + IT_0286;
    const ccomplex_t IT_0288 = conj(IT_0059)*IT_0216;
    const ccomplex_t IT_0289 = (-6)*IT_0146;
    const ccomplex_t IT_0290 = conj(IT_0141)*IT_0289;
    const ccomplex_t IT_0291 = (-6)*IT_0140;
    const ccomplex_t IT_0292 = conj(IT_0193)*IT_0291;
    const ccomplex_t IT_0293 = 0.5*IT_0206;
    const ccomplex_t IT_0294 = s_34 + IT_0293;
    const ccomplex_t IT_0295 = m_C_1*IT_0294;
    const ccomplex_t IT_0296 = 2*IT_0295;
    const ccomplex_t IT_0297 = (-6)*IT_0296;
    const ccomplex_t IT_0298 = (-6)*IT_0261;
    const ccomplex_t IT_0299 = s_13*IT_0145;
    const ccomplex_t IT_0300 = s_13*s_14;
    const ccomplex_t IT_0301 = s_34*IT_0072;
    const ccomplex_t IT_0302 = -IT_0301;
    const ccomplex_t IT_0303 = IT_0300 + IT_0302;
    const ccomplex_t IT_0304 = IT_0145*IT_0215;
    const ccomplex_t IT_0305 = (-6)*IT_0304;
    const ccomplex_t IT_0306 = 0.166666666666667*conj(IT_0265);
    const ccomplex_t IT_0307 = IT_0194*IT_0215;
    const ccomplex_t IT_0308 = 3*IT_0307;
    const ccomplex_t IT_0309 = pow(m_W, 4);
    const ccomplex_t IT_0310 = s_13*IT_0309;
    const ccomplex_t IT_0311 = IT_0136*IT_0310;
    const ccomplex_t IT_0312 = 24*IT_0249;
    const ccomplex_t IT_0313 = (-12)*IT_0248;
    const ccomplex_t IT_0314 = 6*IT_0311 + IT_0312 + IT_0313;
    const ccomplex_t IT_0315 = (-3)*IT_0307;
    const ccomplex_t IT_0316 = (-6)*IT_0315;
    const ccomplex_t IT_0317 = 0.166666666666667*conj(IT_0262);
    const ccomplex_t IT_0318 = (-6)*IT_0311;
    const ccomplex_t IT_0319 = IT_0312 + IT_0318;
    const ccomplex_t IT_0320 = conj(IT_0193)*IT_0303;
    const ccomplex_t IT_0321 = 6*(IT_0059 + IT_0135)*IT_0267 + 6*IT_0265*
      (IT_0256*IT_0264 + -conj(IT_0265)*IT_0272 + IT_0246*IT_0273 +
       0.166666666666667*conj(IT_0193)*IT_0275 + 0.166666666666667*IT_0276
      *IT_0281 + 0.166666666666667*conj(IT_0141)*IT_0283) + IT_0246*IT_0284 + 6
      *IT_0219*(0.166666666666667*conj(IT_0135)*IT_0216 + IT_0287 +
       0.166666666666667*IT_0288 + 0.166666666666667*IT_0290 + 0.166666666666667
      *IT_0292 + IT_0203*IT_0297 + IT_0213*IT_0298) + 6*IT_0245*
      (0.166666666666667*conj(IT_0135)*IT_0216 + IT_0247*IT_0255 + IT_0287 +
       0.166666666666667*IT_0288 + 0.166666666666667*IT_0290 + 0.166666666666667
      *IT_0292 + IT_0203*IT_0297 + IT_0213*IT_0298) + 6*IT_0193*(IT_0146*IT_0256
       + IT_0247*IT_0291 + conj(IT_0193)*IT_0299 + IT_0276*IT_0303 +
       0.166666666666667*conj(IT_0141)*IT_0305 + IT_0275*IT_0306) + 6*IT_0141*
      (IT_0140*IT_0256 + IT_0247*IT_0289 + conj(IT_0141)*IT_0299 + conj(IT_0212)
      *IT_0303 + 0.166666666666667*conj(IT_0193)*IT_0305 + IT_0283*IT_0306) + 6
      *IT_0212*(IT_0211*IT_0257 + IT_0247*IT_0298 + conj(IT_0141)*IT_0303 + conj
      (IT_0202)*IT_0308 + IT_0213*IT_0314 + IT_0316*IT_0317) + 6*IT_0262*
      (IT_0210*IT_0246 + IT_0256*IT_0261 + IT_0281*IT_0306 + IT_0213*IT_0316 +
       IT_0314*IT_0317 + IT_0203*IT_0319 + IT_0320) + 6*IT_0202*(IT_0201*IT_0257
       + IT_0247*IT_0297 + IT_0281*IT_0306 + conj(IT_0212)*IT_0308 + IT_0203
      *IT_0314 + IT_0317*IT_0319 + IT_0320);
    return create_ccomplex_return(IT_0321);
}

