#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_su_L_ss_L_to_s_u.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_su_L_ss_L_to_s_u(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t g_s = param->g_s;
    const creal_t m_W = param->m_W;
    const creal_t m_s = param->m_s;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_sG = param->m_sG;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_L = param->m_ss_L;
    const creal_t m_su_L = param->m_su_L;
    const creal_t theta_W = param->theta_W;
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
    const ccomplex_t IT_0000 = cos(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_B1)*e_em;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 1.4142135623731*IT_0003;
    const ccomplex_t IT_0005 = sin(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = conj(N_W1)*e_em;
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0004 + 3*IT_0009);
    const ccomplex_t IT_0011 = 0.166666666666667*IT_0010;
    const ccomplex_t IT_0012 = pow(m_W, -1);
    const ccomplex_t IT_0013 = cos(beta);
    const ccomplex_t IT_0014 = cpow(IT_0013, -1);
    const ccomplex_t IT_0015 = (0 + _Complex_I*1.4142135623731)*m_s*N_d1*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0016 = (-0.5)*IT_0015;
    const ccomplex_t IT_0017 = IT_0011*IT_0016;
    const ccomplex_t IT_0018 = pow(m_s, 2);
    const ccomplex_t IT_0019 = pow(m_ss_L, 2);
    const ccomplex_t IT_0020 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0021 = IT_0017*IT_0020;
    const ccomplex_t IT_0022 = (0 + _Complex_I*1)*IT_0021;
    const ccomplex_t IT_0023 = conj(N_B2)*e_em;
    const ccomplex_t IT_0024 = IT_0001*IT_0023;
    const ccomplex_t IT_0025 = 1.4142135623731*IT_0024;
    const ccomplex_t IT_0026 = conj(N_W2)*e_em;
    const ccomplex_t IT_0027 = IT_0006*IT_0026;
    const ccomplex_t IT_0028 = 1.4142135623731*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*(IT_0025 + 3*IT_0028);
    const ccomplex_t IT_0030 = 0.166666666666667*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1.4142135623731)*m_s*N_d2*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0032 = (-0.5)*IT_0031;
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0035 = IT_0033*IT_0034;
    const ccomplex_t IT_0036 = (0 + _Complex_I*1)*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B3)*e_em;
    const ccomplex_t IT_0038 = IT_0001*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W3)*e_em;
    const ccomplex_t IT_0041 = IT_0006*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*(IT_0039 + 3*IT_0042);
    const ccomplex_t IT_0044 = 0.166666666666667*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_s*N_d3*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = IT_0044*IT_0046;
    const ccomplex_t IT_0048 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*IT_0049;
    const ccomplex_t IT_0051 = conj(N_B4)*e_em;
    const ccomplex_t IT_0052 = IT_0001*IT_0051;
    const ccomplex_t IT_0053 = 1.4142135623731*IT_0052;
    const ccomplex_t IT_0054 = conj(N_W4)*e_em;
    const ccomplex_t IT_0055 = IT_0006*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*(IT_0053 + 3*IT_0056);
    const ccomplex_t IT_0058 = 0.166666666666667*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1.4142135623731)*m_s*N_d4*e_em
      *IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0060 = (-0.5)*IT_0059;
    const ccomplex_t IT_0061 = IT_0058*IT_0060;
    const ccomplex_t IT_0062 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0063 = IT_0061*IT_0062;
    const ccomplex_t IT_0064 = (0 + _Complex_I*1)*IT_0063;
    const ccomplex_t IT_0065 = -IT_0022 + -IT_0036 + -IT_0050 + -IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm1)*IT_0006;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1.4142135623731)*m_s*U_d1*V_us
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0068 = 0.5*IT_0067;
    const ccomplex_t IT_0069 = IT_0066*IT_0068;
    const ccomplex_t IT_0070 = pow(m_su_L, 2);
    const ccomplex_t IT_0071 = cpow((-2)*s_13 + IT_0018 + IT_0070 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0069*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*V_us*e_em*conj(U_Wm2)*IT_0006;
    const ccomplex_t IT_0075 = (0 + _Complex_I*1.4142135623731)*m_s*U_d2*V_us
      *e_em*IT_0006*IT_0012*IT_0014;
    const ccomplex_t IT_0076 = 0.5*IT_0075;
    const ccomplex_t IT_0077 = IT_0074*IT_0076;
    const ccomplex_t IT_0078 = cpow((-2)*s_13 + IT_0018 + IT_0070 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = IT_0073 + IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1)*(IT_0025 + (-3)*IT_0028);
    const ccomplex_t IT_0083 = 0.166666666666667*IT_0082;
    const ccomplex_t IT_0084 = sin(beta);
    const ccomplex_t IT_0085 = cpow(IT_0084, -1);
    const ccomplex_t IT_0086 = (0 + _Complex_I*1.4142135623731)*m_u*N_u2*e_em
      *IT_0006*IT_0012*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = IT_0083*IT_0087;
    const ccomplex_t IT_0089 = IT_0034*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*(IT_0004 + (-3)*IT_0009);
    const ccomplex_t IT_0092 = 0.166666666666667*IT_0091;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1.4142135623731)*m_u*N_u1*e_em
      *IT_0006*IT_0012*IT_0085;
    const ccomplex_t IT_0094 = (-0.5)*IT_0093;
    const ccomplex_t IT_0095 = IT_0092*IT_0094;
    const ccomplex_t IT_0096 = IT_0020*IT_0095;
    const ccomplex_t IT_0097 = (0 + _Complex_I*1)*IT_0096;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0039 + (-3)*IT_0042);
    const ccomplex_t IT_0099 = 0.166666666666667*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*m_u*N_u3*e_em
      *IT_0006*IT_0012*IT_0085;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = IT_0099*IT_0101;
    const ccomplex_t IT_0103 = IT_0048*IT_0102;
    const ccomplex_t IT_0104 = (0 + _Complex_I*1)*IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*(IT_0053 + (-3)*IT_0056);
    const ccomplex_t IT_0106 = 0.166666666666667*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1.4142135623731)*m_u*N_u4*e_em
      *IT_0006*IT_0012*IT_0085;
    const ccomplex_t IT_0108 = (-0.5)*IT_0107;
    const ccomplex_t IT_0109 = IT_0106*IT_0108;
    const ccomplex_t IT_0110 = IT_0062*IT_0109;
    const ccomplex_t IT_0111 = (0 + _Complex_I*1)*IT_0110;
    const ccomplex_t IT_0112 = -IT_0090 + -IT_0097 + -IT_0104 + -IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1.4142135623731)*m_u*V_u1*V_us
      *e_em*IT_0006*IT_0012*IT_0085;
    const ccomplex_t IT_0114 = 0.5*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*V_us*e_em*conj(V_Wp1)*IT_0006;
    const ccomplex_t IT_0116 = IT_0114*IT_0115;
    const ccomplex_t IT_0117 = IT_0071*IT_0116;
    const ccomplex_t IT_0118 = (0 + _Complex_I*1)*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1.4142135623731)*m_u*V_u2*V_us
      *e_em*IT_0006*IT_0012*IT_0085;
    const ccomplex_t IT_0120 = 0.5*IT_0119;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*V_us*e_em*conj(V_Wp2)*IT_0006;
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = IT_0078*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*IT_0123;
    const ccomplex_t IT_0125 = IT_0118 + IT_0124;
    const ccomplex_t IT_0126 = s_34*IT_0070;
    const ccomplex_t IT_0127 = (-6)*IT_0126;
    const ccomplex_t IT_0128 = s_13*s_14;
    const ccomplex_t IT_0129 = 12*IT_0128;
    const ccomplex_t IT_0130 = IT_0127 + IT_0129;
    const ccomplex_t IT_0131 = (-18)*IT_0126;
    const ccomplex_t IT_0132 = 36*IT_0128;
    const ccomplex_t IT_0133 = IT_0131 + IT_0132;
    const ccomplex_t IT_0134 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0135 = IT_0044*IT_0099;
    const ccomplex_t IT_0136 = IT_0048*IT_0134*IT_0135;
    const ccomplex_t IT_0137 = -IT_0136;
    const ccomplex_t IT_0138 = 18*s_34;
    const ccomplex_t IT_0139 = IT_0137*IT_0138;
    const ccomplex_t IT_0140 = pow(m_sG, 2);
    const ccomplex_t IT_0141 = cpow((-2)*s_23 + IT_0018 + IT_0019 + -IT_0140 +
       reg_prop, -1);
    const ccomplex_t IT_0142 = (0 + _Complex_I*1.4142135623731)*g_s;
    const ccomplex_t IT_0143 = cpow(IT_0142, 2);
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_sG*IT_0143;
    const ccomplex_t IT_0145 = IT_0141*IT_0144;
    const ccomplex_t IT_0146 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0147 = IT_0058*IT_0106;
    const ccomplex_t IT_0148 = IT_0062*IT_0146*IT_0147;
    const ccomplex_t IT_0149 = m_u*IT_0104;
    const ccomplex_t IT_0150 = m_u*IT_0111;
    const ccomplex_t IT_0151 = (-0.166666666666667)*IT_0145 + IT_0148 + 
      -IT_0149 + -IT_0150;
    const ccomplex_t IT_0152 = (-18)*s_34;
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = IT_0139 + IT_0153;
    const ccomplex_t IT_0155 = IT_0066*IT_0115;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0157 = IT_0071*IT_0155*IT_0156;
    const ccomplex_t IT_0158 = IT_0090 + IT_0097;
    const ccomplex_t IT_0159 = m_u*IT_0158;
    const ccomplex_t IT_0160 = -IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0162 = IT_0030*IT_0083;
    const ccomplex_t IT_0163 = IT_0034*IT_0161*IT_0162;
    const ccomplex_t IT_0164 = IT_0068*IT_0114;
    const ccomplex_t IT_0165 = IT_0071*IT_0156*IT_0164;
    const ccomplex_t IT_0166 = IT_0046*IT_0101;
    const ccomplex_t IT_0167 = IT_0048*IT_0134*IT_0166;
    const ccomplex_t IT_0168 = m_u*IT_0022;
    const ccomplex_t IT_0169 = m_u*IT_0036;
    const ccomplex_t IT_0170 = m_u*IT_0050;
    const ccomplex_t IT_0171 = m_u*IT_0064;
    const ccomplex_t IT_0172 = IT_0167 + -IT_0168 + -IT_0169 + -IT_0170 + 
      -IT_0171;
    const ccomplex_t IT_0173 = IT_0032*IT_0087;
    const ccomplex_t IT_0174 = IT_0034*IT_0161*IT_0173;
    const ccomplex_t IT_0175 = IT_0172 + IT_0174;
    const ccomplex_t IT_0176 = conj(IT_0172) + conj(IT_0174);
    const ccomplex_t IT_0177 = conj(IT_0160) + conj(IT_0163);
    const ccomplex_t IT_0178 = IT_0060*IT_0108;
    const ccomplex_t IT_0179 = IT_0062*IT_0146*IT_0178;
    const ccomplex_t IT_0180 = -IT_0179;
    const ccomplex_t IT_0181 = IT_0016*IT_0094;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0183 = IT_0020*IT_0181*IT_0182;
    const ccomplex_t IT_0184 = -IT_0183;
    const ccomplex_t IT_0185 = conj(IT_0180) + conj(IT_0184);
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0187 = IT_0076*IT_0120;
    const ccomplex_t IT_0188 = IT_0078*IT_0186*IT_0187;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = m_s*IT_0125;
    const ccomplex_t IT_0191 = IT_0189 + IT_0190;
    const ccomplex_t IT_0192 = IT_0180 + IT_0184;
    const ccomplex_t IT_0193 = conj(IT_0189) + conj(IT_0190);
    const ccomplex_t IT_0194 = 6*s_34;
    const ccomplex_t IT_0195 = m_s*IT_0081;
    const ccomplex_t IT_0196 = IT_0074*IT_0121;
    const ccomplex_t IT_0197 = IT_0078*IT_0186*IT_0196;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = conj(IT_0195) + conj(IT_0198);
    const ccomplex_t IT_0200 = IT_0011*IT_0092;
    const ccomplex_t IT_0201 = IT_0020*IT_0182*IT_0200;
    const ccomplex_t IT_0202 = -IT_0201;
    const ccomplex_t IT_0203 = conj(IT_0137) + conj(IT_0202);
    const ccomplex_t IT_0204 = 0.5*IT_0145;
    const ccomplex_t IT_0205 = IT_0138*IT_0204;
    const ccomplex_t IT_0206 = IT_0138*IT_0189;
    const ccomplex_t IT_0207 = IT_0138*IT_0151;
    const ccomplex_t IT_0208 = cpow(IT_0204, 2);
    const ccomplex_t IT_0209 = IT_0195 + IT_0198;
    const ccomplex_t IT_0210 = IT_0152*conj(IT_0172);
    const ccomplex_t IT_0211 = conj(IT_0151)*IT_0152;
    const ccomplex_t IT_0212 = IT_0194*IT_0204;
    const ccomplex_t IT_0213 = IT_0160 + IT_0163;
    const ccomplex_t IT_0214 = IT_0137 + IT_0202;
    const ccomplex_t IT_0215 = (-6)*s_34;
    const ccomplex_t IT_0216 = conj(IT_0172)*IT_0215;
    const ccomplex_t IT_0217 = IT_0204*IT_0215;
    const ccomplex_t IT_0218 = IT_0138*IT_0160;
    const ccomplex_t IT_0219 = IT_0137*IT_0152;
    const ccomplex_t IT_0220 = IT_0152*IT_0202;
    const ccomplex_t IT_0221 = IT_0198*IT_0215;
    const ccomplex_t IT_0222 = IT_0138*conj(IT_0151);
    const ccomplex_t IT_0223 = conj(IT_0137)*IT_0152;
    const ccomplex_t IT_0224 = IT_0152*conj(IT_0202);
    const ccomplex_t IT_0225 = conj(IT_0198)*IT_0215;
    const ccomplex_t IT_0226 = IT_0137*IT_0194;
    const ccomplex_t IT_0227 = IT_0194*IT_0202;
    const ccomplex_t IT_0228 = IT_0138*IT_0198;
    const ccomplex_t IT_0229 = IT_0152*IT_0204;
    const ccomplex_t IT_0230 = IT_0151*IT_0215;
    const ccomplex_t IT_0231 = conj(IT_0137)*IT_0194;
    const ccomplex_t IT_0232 = IT_0194*conj(IT_0202);
    const ccomplex_t IT_0233 = conj(IT_0151)*IT_0215;
    const ccomplex_t IT_0234 = conj(IT_0151)*IT_0194;
    const ccomplex_t IT_0235 = conj(IT_0137)*IT_0215;
    const ccomplex_t IT_0236 = m_s*m_u*IT_0070;
    const ccomplex_t IT_0237 = (-18)*IT_0081;
    const ccomplex_t IT_0238 = (-6)*IT_0065;
    const ccomplex_t IT_0239 = IT_0237 + IT_0238;
    const ccomplex_t IT_0240 = (-18)*conj(IT_0081);
    const ccomplex_t IT_0241 = (-6)*conj(IT_0065);
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = (-18)*IT_0112;
    const ccomplex_t IT_0244 = (-18)*conj(IT_0112);
    const ccomplex_t IT_0245 = (-6)*IT_0112;
    const ccomplex_t IT_0246 = (-6)*conj(IT_0112);
    const ccomplex_t IT_0247 = m_u*s_13;
    const ccomplex_t IT_0248 = 6*IT_0112;
    const ccomplex_t IT_0249 = 18*IT_0125;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = 6*conj(IT_0112);
    const ccomplex_t IT_0252 = 18*conj(IT_0125);
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = 6*IT_0065;
    const ccomplex_t IT_0255 = 18*IT_0081;
    const ccomplex_t IT_0256 = IT_0254 + IT_0255;
    const ccomplex_t IT_0257 = 6*conj(IT_0065);
    const ccomplex_t IT_0258 = 18*conj(IT_0081);
    const ccomplex_t IT_0259 = IT_0257 + IT_0258;
    const ccomplex_t IT_0260 = (-18)*IT_0125;
    const ccomplex_t IT_0261 = IT_0245 + IT_0260;
    const ccomplex_t IT_0262 = (-18)*conj(IT_0125);
    const ccomplex_t IT_0263 = IT_0246 + IT_0262;
    const ccomplex_t IT_0264 = 6*IT_0189;
    const ccomplex_t IT_0265 = 18*IT_0180;
    const ccomplex_t IT_0266 = 18*IT_0184;
    const ccomplex_t IT_0267 = (-18)*IT_0172;
    const ccomplex_t IT_0268 = 6*conj(IT_0189);
    const ccomplex_t IT_0269 = 18*conj(IT_0180);
    const ccomplex_t IT_0270 = 18*conj(IT_0184);
    const ccomplex_t IT_0271 = (-18)*conj(IT_0172);
    const ccomplex_t IT_0272 = 6*IT_0204;
    const ccomplex_t IT_0273 = 6*conj(IT_0198);
    const ccomplex_t IT_0274 = 18*conj(IT_0137);
    const ccomplex_t IT_0275 = 18*conj(IT_0202);
    const ccomplex_t IT_0276 = (-18)*conj(IT_0160);
    const ccomplex_t IT_0277 = (-18)*conj(IT_0151);
    const ccomplex_t IT_0278 = 6*IT_0180;
    const ccomplex_t IT_0279 = 6*IT_0184;
    const ccomplex_t IT_0280 = 18*IT_0189;
    const ccomplex_t IT_0281 = (-6)*IT_0172;
    const ccomplex_t IT_0282 = 6*conj(IT_0180);
    const ccomplex_t IT_0283 = 6*conj(IT_0184);
    const ccomplex_t IT_0284 = 18*conj(IT_0189);
    const ccomplex_t IT_0285 = (-6)*conj(IT_0172);
    const ccomplex_t IT_0286 = 6*IT_0198;
    const ccomplex_t IT_0287 = 18*IT_0137;
    const ccomplex_t IT_0288 = 18*IT_0202;
    const ccomplex_t IT_0289 = (-18)*IT_0160;
    const ccomplex_t IT_0290 = (-18)*IT_0151;
    const ccomplex_t IT_0291 = (-6)*IT_0204;
    const ccomplex_t IT_0292 = 6*IT_0137;
    const ccomplex_t IT_0293 = 6*IT_0202;
    const ccomplex_t IT_0294 = 18*IT_0198;
    const ccomplex_t IT_0295 = (-18)*IT_0204;
    const ccomplex_t IT_0296 = (-6)*IT_0160;
    const ccomplex_t IT_0297 = (-6)*IT_0151;
    const ccomplex_t IT_0298 = 6*conj(IT_0137);
    const ccomplex_t IT_0299 = 6*conj(IT_0202);
    const ccomplex_t IT_0300 = 18*IT_0204;
    const ccomplex_t IT_0301 = 18*conj(IT_0198);
    const ccomplex_t IT_0302 = (-6)*conj(IT_0160);
    const ccomplex_t IT_0303 = (-6)*conj(IT_0151);
    const ccomplex_t IT_0304 = conj(IT_0165)*IT_0239 + IT_0165*IT_0242 + conj
      (IT_0163)*IT_0243 + IT_0163*IT_0244 + conj(IT_0195)*IT_0250 + IT_0195
      *IT_0253 + conj(IT_0190)*IT_0256 + IT_0190*IT_0259 + conj(IT_0157)*IT_0261
       + IT_0157*IT_0263 + conj(IT_0065)*((-18)*IT_0174 + IT_0264 + IT_0265 +
       IT_0266 + IT_0267) + IT_0065*((-18)*conj(IT_0174) + IT_0268 + IT_0269 +
       IT_0270 + IT_0271) + IT_0112*(IT_0272 + IT_0273 + IT_0274 + IT_0275 +
       IT_0276 + IT_0277) + (-6)*conj(IT_0081)*(IT_0174 + (-0.166666666666667)
      *IT_0278 + (-0.166666666666667)*IT_0279 + (-0.166666666666667)*IT_0280 + (
      -0.166666666666667)*IT_0281) + (-6)*IT_0081*(conj(IT_0174) + (
      -0.166666666666667)*IT_0282 + (-0.166666666666667)*IT_0283 + (
      -0.166666666666667)*IT_0284 + (-0.166666666666667)*IT_0285) + conj(IT_0112
      )*(IT_0286 + IT_0287 + IT_0288 + IT_0289 + IT_0290 + IT_0291) + (-6)*conj
      (IT_0125)*(IT_0163 + (-0.166666666666667)*IT_0292 + (-0.166666666666667)
      *IT_0293 + (-0.166666666666667)*IT_0294 + (-0.166666666666667)*IT_0295 + (
      -0.166666666666667)*IT_0296 + (-0.166666666666667)*IT_0297) + (-6)*IT_0125
      *(conj(IT_0163) + (-0.166666666666667)*IT_0298 + (-0.166666666666667)
      *IT_0299 + (-0.166666666666667)*IT_0300 + (-0.166666666666667)*IT_0301 + (
      -0.166666666666667)*IT_0302 + (-0.166666666666667)*IT_0303);
    const ccomplex_t IT_0305 = m_s*s_14;
    const ccomplex_t IT_0306 = 6*IT_0172;
    const ccomplex_t IT_0307 = (-18)*IT_0189;
    const ccomplex_t IT_0308 = (-6)*IT_0180;
    const ccomplex_t IT_0309 = (-6)*IT_0184;
    const ccomplex_t IT_0310 = 6*conj(IT_0172);
    const ccomplex_t IT_0311 = (-18)*conj(IT_0189);
    const ccomplex_t IT_0312 = (-6)*conj(IT_0180);
    const ccomplex_t IT_0313 = (-6)*conj(IT_0184);
    const ccomplex_t IT_0314 = 18*IT_0172;
    const ccomplex_t IT_0315 = (-18)*IT_0180;
    const ccomplex_t IT_0316 = (-18)*IT_0184;
    const ccomplex_t IT_0317 = (-6)*IT_0189;
    const ccomplex_t IT_0318 = 18*conj(IT_0172);
    const ccomplex_t IT_0319 = (-18)*conj(IT_0180);
    const ccomplex_t IT_0320 = (-18)*conj(IT_0184);
    const ccomplex_t IT_0321 = (-6)*conj(IT_0189);
    const ccomplex_t IT_0322 = 6*IT_0160;
    const ccomplex_t IT_0323 = 6*IT_0151;
    const ccomplex_t IT_0324 = (-18)*IT_0198;
    const ccomplex_t IT_0325 = (-6)*IT_0137;
    const ccomplex_t IT_0326 = (-6)*IT_0202;
    const ccomplex_t IT_0327 = 6*conj(IT_0160);
    const ccomplex_t IT_0328 = 6*conj(IT_0151);
    const ccomplex_t IT_0329 = (-18)*conj(IT_0198);
    const ccomplex_t IT_0330 = (-6)*conj(IT_0137);
    const ccomplex_t IT_0331 = (-6)*conj(IT_0202);
    const ccomplex_t IT_0332 = (-0.166666666666667)*IT_0272;
    const ccomplex_t IT_0333 = (-0.166666666666667)*IT_0291;
    const ccomplex_t IT_0334 = conj(IT_0195)*IT_0239 + IT_0195*IT_0242 + conj
      (IT_0165)*IT_0250 + IT_0165*IT_0253 + conj(IT_0157)*IT_0256 + IT_0157
      *IT_0259 + conj(IT_0190)*IT_0261 + IT_0190*IT_0263 + 6*conj(IT_0125)*
      (IT_0174 + 0.166666666666667*IT_0306 + 0.166666666666667*IT_0307 +
       0.166666666666667*IT_0308 + 0.166666666666667*IT_0309) + 6*IT_0125*(conj
      (IT_0174) + 0.166666666666667*IT_0310 + 0.166666666666667*IT_0311 +
       0.166666666666667*IT_0312 + 0.166666666666667*IT_0313) + 18*conj(IT_0112)
      *(IT_0174 + 0.0555555555555556*IT_0314 + 0.0555555555555556*IT_0315 +
       0.0555555555555556*IT_0316 + 0.0555555555555556*IT_0317) + 18*IT_0112*
      (conj(IT_0174) + 0.0555555555555556*IT_0318 + 0.0555555555555556*IT_0319 +
       0.0555555555555556*IT_0320 + 0.0555555555555556*IT_0321) + 6*conj(IT_0081
      )*(IT_0163 + 0.166666666666667*IT_0300 + 0.166666666666667*IT_0322 +
       0.166666666666667*IT_0323 + 0.166666666666667*IT_0324 + 0.166666666666667
      *IT_0325 + 0.166666666666667*IT_0326) + 6*IT_0081*(conj(IT_0163) +
       0.166666666666667*IT_0295 + 0.166666666666667*IT_0327 + 0.166666666666667
      *IT_0328 + 0.166666666666667*IT_0329 + 0.166666666666667*IT_0330 +
       0.166666666666667*IT_0331) + (-18)*conj(IT_0065)*(IT_0137 + -IT_0151 + 
      -IT_0160 + -IT_0163 + 1./3*IT_0198 + IT_0202 + 1./3*IT_0332) + (-18)
      *IT_0065*(conj(IT_0137) + -conj(IT_0151) + -conj(IT_0160) + -conj(IT_0163)
       + 1./3*conj(IT_0198) + conj(IT_0202) + 1./3*IT_0333);
    const ccomplex_t IT_0335 = m_s*m_u;
    const ccomplex_t IT_0336 = (conj(IT_0151) + conj(IT_0160) + conj(IT_0163))
      *(IT_0264 + IT_0265 + IT_0266 + IT_0267) + (IT_0151 + IT_0160 + IT_0163)*
      (IT_0268 + IT_0269 + IT_0270 + IT_0271) + conj(IT_0157)*(IT_0278 + IT_0279
       + IT_0280 + IT_0281) + IT_0157*(IT_0282 + IT_0283 + IT_0284 + IT_0285) +
       conj(IT_0165)*((-18)*IT_0157 + (-6)*IT_0163 + 18*IT_0195 + IT_0292 +
       IT_0293 + IT_0294 + IT_0295 + IT_0296 + IT_0297) + IT_0165*((-18)*conj
      (IT_0157) + (-6)*conj(IT_0163) + 18*conj(IT_0195) + IT_0298 + IT_0299 +
       IT_0300 + IT_0301 + IT_0302 + IT_0303) + IT_0199*(IT_0306 + IT_0307 +
       IT_0308 + IT_0309) + IT_0204*(IT_0282 + IT_0283 + IT_0284 + IT_0285 +
       IT_0306 + IT_0307 + IT_0308 + IT_0309) + IT_0209*(IT_0310 + IT_0311 +
       IT_0312 + IT_0313) + IT_0203*(IT_0314 + IT_0315 + IT_0316 + IT_0317) +
       IT_0214*(IT_0318 + IT_0319 + IT_0320 + IT_0321) + conj(IT_0190)*(18
      *IT_0157 + 6*IT_0163 + (-18)*IT_0195 + IT_0300 + IT_0322 + IT_0323 +
       IT_0324 + IT_0325 + IT_0326) + IT_0190*(18*conj(IT_0157) + 6*conj(IT_0163
      ) + (-18)*conj(IT_0195) + IT_0295 + IT_0327 + IT_0328 + IT_0329 + IT_0330 
      + IT_0331) + (-6)*IT_0174*(conj(IT_0157) + 3*conj(IT_0163) + -conj(IT_0195
      ) + (-0.166666666666667)*IT_0273 + (-0.166666666666667)*IT_0274 + (
      -0.166666666666667)*IT_0275 + (-0.166666666666667)*IT_0276 + (
      -0.166666666666667)*IT_0277 + IT_0332) + (-6)*conj(IT_0174)*(IT_0157 + 3
      *IT_0163 + -IT_0195 + (-0.166666666666667)*IT_0286 + (-0.166666666666667)
      *IT_0287 + (-0.166666666666667)*IT_0288 + (-0.166666666666667)*IT_0289 + (
      -0.166666666666667)*IT_0290 + IT_0333);
    const ccomplex_t IT_0337 = (conj(IT_0065)*IT_0081 + IT_0065*conj(IT_0081) 
      + conj(IT_0112)*IT_0125 + IT_0112*conj(IT_0125))*IT_0130 + (IT_0065*conj
      (IT_0065) + IT_0081*conj(IT_0081) + IT_0112*conj(IT_0112) + IT_0125*conj
      (IT_0125))*IT_0133 + conj(IT_0137)*IT_0154 + (conj(IT_0157)*(IT_0151 +
       IT_0160 + IT_0163) + conj(IT_0165)*IT_0175 + IT_0165*IT_0176 + IT_0157
      *IT_0177 + IT_0185*IT_0191 + IT_0192*IT_0193)*IT_0194 + IT_0138*(IT_0157
      *conj(IT_0157) + IT_0165*conj(IT_0165) + IT_0175*IT_0176 + IT_0163*IT_0177
       + IT_0185*IT_0192 + IT_0190*IT_0193 + IT_0195*IT_0199 + IT_0202*IT_0203) 
      + conj(IT_0157)*IT_0205 + IT_0193*IT_0206 + conj(IT_0151)*IT_0207 +
       IT_0152*(IT_0175*IT_0185 + conj(IT_0165)*IT_0191 + conj(IT_0174)*IT_0192 
      + IT_0165*IT_0193 + IT_0157*IT_0199 + IT_0208 + conj(IT_0157)*IT_0209) +
       IT_0192*IT_0210 + IT_0137*IT_0211 + IT_0202*(IT_0211 + IT_0212) + 
      (IT_0165*IT_0185 + conj(IT_0174)*IT_0191 + conj(IT_0165)*IT_0192 + IT_0175
      *IT_0193 + IT_0177*IT_0195 + IT_0157*conj(IT_0202) + conj(IT_0195)*IT_0213
       + conj(IT_0157)*IT_0214)*IT_0215 + IT_0191*IT_0216 + conj(IT_0202)*
      (IT_0154 + IT_0217) + IT_0177*(IT_0207 + IT_0212 + IT_0218 + IT_0219 +
       IT_0220 + IT_0221) + IT_0213*(IT_0217 + IT_0222 + IT_0223 + IT_0224 +
       IT_0225) + IT_0199*(IT_0226 + IT_0227 + IT_0228 + IT_0229 + IT_0230) +
       IT_0209*(IT_0205 + IT_0231 + IT_0232 + IT_0233) + IT_0157*(IT_0229 +
       IT_0234 + IT_0235) + IT_0204*(IT_0226 + IT_0230 + IT_0234 + IT_0235) +
       IT_0236*(conj(IT_0125)*IT_0239 + IT_0125*IT_0242 + conj(IT_0065)*IT_0243 
      + IT_0065*IT_0244 + conj(IT_0081)*IT_0245 + IT_0081*IT_0246) + IT_0247
      *IT_0304 + IT_0305*IT_0334 + IT_0335*IT_0336;
    return create_ccomplex_return(IT_0337);
}

