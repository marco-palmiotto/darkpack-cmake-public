#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_ss_R_to_s_Hp.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_ss_R_to_s_Hp(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_c = param->m_c;
    const creal_t m_s = param->m_s;
    const creal_t m_t = param->m_t;
    const creal_t m_u = param->m_u;
    const creal_t V_us = param->V_us;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_23 = param->s_23;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_ss_R = param->m_ss_R;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
    const creal_t Gamma_top = param->Gamma_top;
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
    const ccomplex_t V_cs = param->V_cs;
    const ccomplex_t V_ts = param->V_ts;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t IT_0000 = pow(m_s, 2);
    const ccomplex_t IT_0001 = pow(m_ss_R, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0003 = pow(m_W, -1);
    const ccomplex_t IT_0004 = cos(beta);
    const ccomplex_t IT_0005 = cpow(IT_0004, -1);
    const ccomplex_t IT_0006 = sin(theta_W);
    const ccomplex_t IT_0007 = cpow(IT_0006, -1);
    const ccomplex_t IT_0008 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d2)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0009 = 0.5*IT_0008;
    const ccomplex_t IT_0010 = cos(theta_W);
    const ccomplex_t IT_0011 = cpow(IT_0010, -1);
    const ccomplex_t IT_0012 = IT_0004*IT_0011;
    const ccomplex_t IT_0013 = N_B2*V_u1*e_em;
    const ccomplex_t IT_0014 = IT_0012*IT_0013;
    const ccomplex_t IT_0015 = 1.4142135623731*IT_0014;
    const ccomplex_t IT_0016 = IT_0004*IT_0007;
    const ccomplex_t IT_0017 = N_W2*V_u1*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = 1.4142135623731*IT_0018;
    const ccomplex_t IT_0020 = IT_0015 + IT_0019;
    const ccomplex_t IT_0021 = N_u2*e_em*V_Wp1;
    const ccomplex_t IT_0022 = IT_0016*IT_0021;
    const ccomplex_t IT_0023 = 2*IT_0022;
    const ccomplex_t IT_0024 = IT_0020 + IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*IT_0024;
    const ccomplex_t IT_0026 = 0.5*IT_0025;
    const ccomplex_t IT_0027 = IT_0009*IT_0026;
    const ccomplex_t IT_0028 = IT_0002*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*IT_0028;
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d3)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0031 = 0.5*IT_0030;
    const ccomplex_t IT_0032 = N_u3*e_em*V_Wp1;
    const ccomplex_t IT_0033 = IT_0016*IT_0032;
    const ccomplex_t IT_0034 = 2*IT_0033;
    const ccomplex_t IT_0035 = N_B3*V_u1*e_em;
    const ccomplex_t IT_0036 = IT_0012*IT_0035;
    const ccomplex_t IT_0037 = 1.4142135623731*IT_0036;
    const ccomplex_t IT_0038 = N_W3*V_u1*e_em;
    const ccomplex_t IT_0039 = IT_0016*IT_0038;
    const ccomplex_t IT_0040 = 1.4142135623731*IT_0039;
    const ccomplex_t IT_0041 = IT_0037 + IT_0040;
    const ccomplex_t IT_0042 = IT_0034 + IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = 0.5*IT_0043;
    const ccomplex_t IT_0045 = IT_0031*IT_0044;
    const ccomplex_t IT_0046 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0047 = IT_0045*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d4)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0050 = 0.5*IT_0049;
    const ccomplex_t IT_0051 = N_u4*e_em*V_Wp1;
    const ccomplex_t IT_0052 = IT_0016*IT_0051;
    const ccomplex_t IT_0053 = 2*IT_0052;
    const ccomplex_t IT_0054 = N_B4*V_u1*e_em;
    const ccomplex_t IT_0055 = IT_0012*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = N_W4*V_u1*e_em;
    const ccomplex_t IT_0058 = IT_0016*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = IT_0056 + IT_0059;
    const ccomplex_t IT_0061 = IT_0053 + IT_0060;
    const ccomplex_t IT_0062 = (0 + _Complex_I*1)*IT_0061;
    const ccomplex_t IT_0063 = 0.5*IT_0062;
    const ccomplex_t IT_0064 = IT_0050*IT_0063;
    const ccomplex_t IT_0065 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0066 = IT_0064*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = cpow((-2)*s_23 + IT_0000 + IT_0001 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0069 = N_u1*e_em*V_Wp1;
    const ccomplex_t IT_0070 = IT_0016*IT_0069;
    const ccomplex_t IT_0071 = 2*IT_0070;
    const ccomplex_t IT_0072 = N_B1*V_u1*e_em;
    const ccomplex_t IT_0073 = IT_0012*IT_0072;
    const ccomplex_t IT_0074 = 1.4142135623731*IT_0073;
    const ccomplex_t IT_0075 = N_W1*V_u1*e_em;
    const ccomplex_t IT_0076 = IT_0016*IT_0075;
    const ccomplex_t IT_0077 = 1.4142135623731*IT_0076;
    const ccomplex_t IT_0078 = IT_0074 + IT_0077;
    const ccomplex_t IT_0079 = IT_0071 + IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = (0 + _Complex_I*1.4142135623731)*m_s*conj(N_d1)
      *e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = IT_0081*IT_0083;
    const ccomplex_t IT_0085 = IT_0068*IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = -IT_0029 + -IT_0048 + -IT_0067 + -IT_0086;
    const ccomplex_t IT_0088 = m_s*m_C_1*IT_0001;
    const ccomplex_t IT_0089 = pow(m_c, 2);
    const ccomplex_t IT_0090 = pow(m_C_1, 2);
    const ccomplex_t IT_0091 = cpow((-2)*s_12 + -IT_0001 + IT_0089 + -IT_0090 
      + -reg_prop, -1);
    const ccomplex_t IT_0092 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_cs*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0093 = (-0.5)*IT_0092;
    const ccomplex_t IT_0094 = sin(beta);
    const ccomplex_t IT_0095 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_cs)
      *e_em*IT_0003*IT_0005*IT_0007*IT_0094;
    const ccomplex_t IT_0096 = 0.5*IT_0095;
    const ccomplex_t IT_0097 = IT_0093*IT_0096;
    const ccomplex_t IT_0098 = IT_0091*IT_0097;
    const ccomplex_t IT_0099 = (0 + _Complex_I*1)*IT_0098;
    const ccomplex_t IT_0100 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_ts*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0101 = (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1.4142135623731)*m_s*conj(V_ts)
      *e_em*IT_0003*IT_0005*IT_0007*IT_0094;
    const ccomplex_t IT_0103 = 0.5*IT_0102;
    const ccomplex_t IT_0104 = IT_0101*IT_0103;
    const ccomplex_t IT_0105 = cpow((-2)*s_12 + -IT_0001 + -IT_0090 + 
      -reg_prop + m_t*(m_t + (0 + _Complex_I*-1)*Gamma_top), -1);
    const ccomplex_t IT_0106 = IT_0104*IT_0105;
    const ccomplex_t IT_0107 = (0 + _Complex_I*1)*IT_0106;
    const ccomplex_t IT_0108 = pow(m_u, 2);
    const ccomplex_t IT_0109 = cpow((-2)*s_12 + -IT_0001 + -IT_0090 + IT_0108 
      + -reg_prop, -1);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1.4142135623731)*m_s*conj(U_d1)
      *V_us*e_em*IT_0003*IT_0005*IT_0007;
    const ccomplex_t IT_0111 = (-0.5)*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1.4142135623731)*m_s*V_us*e_em
      *IT_0003*IT_0005*IT_0007*IT_0094;
    const ccomplex_t IT_0113 = 0.5*IT_0112;
    const ccomplex_t IT_0114 = IT_0111*IT_0113;
    const ccomplex_t IT_0115 = IT_0109*IT_0114;
    const ccomplex_t IT_0116 = (0 + _Complex_I*1)*IT_0115;
    const ccomplex_t IT_0117 = (0 + _Complex_I*1.4142135623731)*N_B1*e_em
      *IT_0011;
    const ccomplex_t IT_0118 = (-0.333333333333333)*IT_0117;
    const ccomplex_t IT_0119 = IT_0011*IT_0094;
    const ccomplex_t IT_0120 = conj(N_B1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0121 = IT_0119*IT_0120;
    const ccomplex_t IT_0122 = 1.4142135623731*IT_0121;
    const ccomplex_t IT_0123 = IT_0007*IT_0094;
    const ccomplex_t IT_0124 = conj(N_W1)*conj(U_d1)*e_em;
    const ccomplex_t IT_0125 = IT_0123*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = IT_0122 + IT_0126;
    const ccomplex_t IT_0128 = conj(N_d1)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0129 = IT_0123*IT_0128;
    const ccomplex_t IT_0130 = (-2)*IT_0129;
    const ccomplex_t IT_0131 = IT_0127 + IT_0130;
    const ccomplex_t IT_0132 = (0 + _Complex_I*1)*IT_0131;
    const ccomplex_t IT_0133 = 0.5*IT_0132;
    const ccomplex_t IT_0134 = IT_0118*IT_0133;
    const ccomplex_t IT_0135 = IT_0068*IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = (0 + _Complex_I*1.4142135623731)*N_B4*e_em
      *IT_0011;
    const ccomplex_t IT_0138 = (-0.333333333333333)*IT_0137;
    const ccomplex_t IT_0139 = conj(N_d4)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0140 = IT_0123*IT_0139;
    const ccomplex_t IT_0141 = (-2)*IT_0140;
    const ccomplex_t IT_0142 = conj(N_B4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0143 = IT_0119*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = conj(N_W4)*conj(U_d1)*e_em;
    const ccomplex_t IT_0146 = IT_0123*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = IT_0144 + IT_0147;
    const ccomplex_t IT_0149 = IT_0141 + IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = 0.5*IT_0150;
    const ccomplex_t IT_0152 = IT_0138*IT_0151;
    const ccomplex_t IT_0153 = IT_0065*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = (0 + _Complex_I*1.4142135623731)*N_B2*e_em
      *IT_0011;
    const ccomplex_t IT_0156 = (-0.333333333333333)*IT_0155;
    const ccomplex_t IT_0157 = conj(N_d2)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0158 = IT_0123*IT_0157;
    const ccomplex_t IT_0159 = (-2)*IT_0158;
    const ccomplex_t IT_0160 = conj(N_B2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0161 = IT_0119*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = conj(N_W2)*conj(U_d1)*e_em;
    const ccomplex_t IT_0164 = IT_0123*IT_0163;
    const ccomplex_t IT_0165 = 1.4142135623731*IT_0164;
    const ccomplex_t IT_0166 = IT_0162 + IT_0165;
    const ccomplex_t IT_0167 = IT_0159 + IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = 0.5*IT_0168;
    const ccomplex_t IT_0170 = IT_0156*IT_0169;
    const ccomplex_t IT_0171 = IT_0002*IT_0170;
    const ccomplex_t IT_0172 = (0 + _Complex_I*1)*IT_0171;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1.4142135623731)*N_B3*e_em
      *IT_0011;
    const ccomplex_t IT_0174 = (-0.333333333333333)*IT_0173;
    const ccomplex_t IT_0175 = conj(N_W3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0176 = IT_0123*IT_0175;
    const ccomplex_t IT_0177 = 1.4142135623731*IT_0176;
    const ccomplex_t IT_0178 = conj(N_B3)*conj(U_d1)*e_em;
    const ccomplex_t IT_0179 = IT_0119*IT_0178;
    const ccomplex_t IT_0180 = 1.4142135623731*IT_0179;
    const ccomplex_t IT_0181 = IT_0177 + IT_0180;
    const ccomplex_t IT_0182 = conj(N_d3)*e_em*conj(U_Wm1);
    const ccomplex_t IT_0183 = IT_0123*IT_0182;
    const ccomplex_t IT_0184 = (-2)*IT_0183;
    const ccomplex_t IT_0185 = IT_0181 + IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*IT_0185;
    const ccomplex_t IT_0187 = 0.5*IT_0186;
    const ccomplex_t IT_0188 = IT_0174*IT_0187;
    const ccomplex_t IT_0189 = IT_0046*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = IT_0099 + IT_0107 + IT_0116 + -IT_0136 + 
      -IT_0154 + -IT_0172 + -IT_0190;
    const ccomplex_t IT_0192 = m_s*s_12;
    const ccomplex_t IT_0193 = 6*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0195 = IT_0081*IT_0118;
    const ccomplex_t IT_0196 = IT_0068*IT_0194*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0198 = IT_0026*IT_0156;
    const ccomplex_t IT_0199 = IT_0002*IT_0197*IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0201 = IT_0044*IT_0174;
    const ccomplex_t IT_0202 = IT_0046*IT_0200*IT_0201;
    const ccomplex_t IT_0203 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0204 = IT_0063*IT_0138;
    const ccomplex_t IT_0205 = IT_0065*IT_0203*IT_0204;
    const ccomplex_t IT_0206 = m_s*IT_0029;
    const ccomplex_t IT_0207 = m_s*IT_0048;
    const ccomplex_t IT_0208 = m_s*IT_0067;
    const ccomplex_t IT_0209 = m_C_1*IT_0099;
    const ccomplex_t IT_0210 = m_C_1*IT_0107;
    const ccomplex_t IT_0211 = m_C_1*IT_0116;
    const ccomplex_t IT_0212 = m_s*IT_0086;
    const ccomplex_t IT_0213 = -IT_0196 + -IT_0199 + -IT_0202 + -IT_0205 +
       IT_0206 + IT_0207 + IT_0208 + IT_0209 + IT_0210 + IT_0211 + IT_0212;
    const ccomplex_t IT_0214 = 0.166666666666667*conj(IT_0213);
    const ccomplex_t IT_0215 = s_23*m_C_1;
    const ccomplex_t IT_0216 = 6*IT_0215;
    const ccomplex_t IT_0217 = cpow(IT_0094, -1);
    const ccomplex_t IT_0218 = (0 + _Complex_I*1.4142135623731)*m_c*conj(V_cs)
      *e_em*IT_0003*IT_0004*IT_0007*IT_0217;
    const ccomplex_t IT_0219 = 0.5*IT_0218;
    const ccomplex_t IT_0220 = IT_0093*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*m_c;
    const ccomplex_t IT_0222 = IT_0091*IT_0220*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1.4142135623731)*m_u*V_us*e_em
      *IT_0003*IT_0004*IT_0007*IT_0217;
    const ccomplex_t IT_0224 = 0.5*IT_0223;
    const ccomplex_t IT_0225 = IT_0111*IT_0224;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*m_u;
    const ccomplex_t IT_0227 = IT_0109*IT_0225*IT_0226;
    const ccomplex_t IT_0228 = IT_0083*IT_0133;
    const ccomplex_t IT_0229 = IT_0068*IT_0194*IT_0228;
    const ccomplex_t IT_0230 = IT_0009*IT_0169;
    const ccomplex_t IT_0231 = IT_0002*IT_0197*IT_0230;
    const ccomplex_t IT_0232 = IT_0031*IT_0187;
    const ccomplex_t IT_0233 = IT_0046*IT_0200*IT_0232;
    const ccomplex_t IT_0234 = IT_0050*IT_0151;
    const ccomplex_t IT_0235 = IT_0065*IT_0203*IT_0234;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1.4142135623731)*m_t*conj(V_ts)
      *e_em*IT_0003*IT_0004*IT_0007*IT_0217;
    const ccomplex_t IT_0237 = 0.5*IT_0236;
    const ccomplex_t IT_0238 = IT_0101*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*m_t;
    const ccomplex_t IT_0240 = IT_0105*IT_0238*IT_0239;
    const ccomplex_t IT_0241 = m_s*IT_0136;
    const ccomplex_t IT_0242 = m_s*IT_0154;
    const ccomplex_t IT_0243 = m_s*IT_0172;
    const ccomplex_t IT_0244 = m_s*IT_0190;
    const ccomplex_t IT_0245 = IT_0222 + IT_0227 + -IT_0229 + -IT_0231 + 
      -IT_0233 + -IT_0235 + IT_0240 + IT_0241 + IT_0242 + IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = 0.166666666666667*conj(IT_0245);
    const ccomplex_t IT_0247 = s_13*IT_0001;
    const ccomplex_t IT_0248 = (-6)*IT_0247;
    const ccomplex_t IT_0249 = s_12*s_23;
    const ccomplex_t IT_0250 = 12*IT_0249;
    const ccomplex_t IT_0251 = IT_0248 + IT_0250;
    const ccomplex_t IT_0252 = 0.166666666666667*conj(IT_0087);
    const ccomplex_t IT_0253 = 0.166666666666667*conj(IT_0191);
    const ccomplex_t IT_0254 = m_s*m_C_1;
    const ccomplex_t IT_0255 = 6*IT_0254;
    const ccomplex_t IT_0256 = 6*IT_0087*(IT_0088*conj(IT_0191) + IT_0193
      *IT_0214 + IT_0216*IT_0246 + IT_0251*IT_0252) + 6*IT_0191*(conj(IT_0087)
      *IT_0088 + IT_0214*IT_0216 + IT_0193*IT_0246 + IT_0251*IT_0253) + 6
      *IT_0245*(s_13*conj(IT_0245) + IT_0216*IT_0252 + IT_0193*IT_0253 + IT_0214
      *IT_0255) + 6*IT_0213*(s_13*conj(IT_0213) + IT_0193*IT_0252 + IT_0216
      *IT_0253 + IT_0246*IT_0255);
    return create_ccomplex_return(IT_0256);
}

