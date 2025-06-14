#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_1_anti_C_2_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_1_anti_C_2_to_anti_W_W(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_N_1 = param->m_N_1;
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_o1 = param->Gamma_o1;
    const creal_t Gamma_o2 = param->Gamma_o2;
    const creal_t Gamma_o3 = param->Gamma_o3;
    const creal_t Gamma_o4 = param->Gamma_o4;
    const creal_t reg_prop = param->reg_prop;
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
    const ccomplex_t IT_0001 = sin(theta_W);
    const ccomplex_t IT_0002 = cpow(IT_0001, -1);
    const ccomplex_t IT_0003 = (0 + _Complex_I*1)*e_em*IT_0000*IT_0002;
    const ccomplex_t IT_0004 = cpow(IT_0000, -1);
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = conj(U_d1)*U_d2*e_em;
    const ccomplex_t IT_0007 = IT_0005*IT_0006;
    const ccomplex_t IT_0008 = IT_0000*IT_0002;
    const ccomplex_t IT_0009 = IT_0006*IT_0008;
    const ccomplex_t IT_0010 = e_em*conj(U_Wm1)*U_Wm2;
    const ccomplex_t IT_0011 = IT_0008*IT_0010;
    const ccomplex_t IT_0012 = (0 + _Complex_I*1)*(IT_0007 + -IT_0009 + (-2)
      *IT_0011);
    const ccomplex_t IT_0013 = 0.5*IT_0012;
    const ccomplex_t IT_0014 = pow(m_C_1, 2);
    const ccomplex_t IT_0015 = pow(m_C_2, 2);
    const ccomplex_t IT_0016 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0014 + -IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0017 = IT_0013*IT_0016;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*IT_0017;
    const ccomplex_t IT_0019 = IT_0003*IT_0018;
    const ccomplex_t IT_0020 = (-2)*IT_0003;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = 2*IT_0003;
    const ccomplex_t IT_0023 = V_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0024 = IT_0005*IT_0023;
    const ccomplex_t IT_0025 = IT_0008*IT_0023;
    const ccomplex_t IT_0026 = e_em*V_Wp1*conj(V_Wp2);
    const ccomplex_t IT_0027 = IT_0008*IT_0026;
    const ccomplex_t IT_0028 = (0 + _Complex_I*1)*(IT_0024 + -IT_0025 + (-2)
      *IT_0027);
    const ccomplex_t IT_0029 = (-0.5)*IT_0028;
    const ccomplex_t IT_0030 = IT_0016*IT_0029;
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0022*IT_0031;
    const ccomplex_t IT_0033 = pow(m_W, -2);
    const ccomplex_t IT_0034 = s_34*m_C_1*m_C_2*IT_0033;
    const ccomplex_t IT_0035 = s_13*s_34*IT_0033;
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = s_14 + IT_0036;
    const ccomplex_t IT_0038 = IT_0034*IT_0037;
    const ccomplex_t IT_0039 = 2*IT_0038;
    const ccomplex_t IT_0040 = IT_0032*IT_0039;
    const ccomplex_t IT_0041 = pow(m_W, 2);
    const ccomplex_t IT_0042 = cpow((-2)*s_23 + IT_0015 + IT_0041 + -m_N_1*
      (m_N_1 + (0 + _Complex_I*-1)*Gamma_o1) + reg_prop, -1);
    const ccomplex_t IT_0043 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0044 = IT_0002*IT_0043;
    const ccomplex_t IT_0045 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0046 = IT_0002*IT_0045;
    const ccomplex_t IT_0047 = 1.4142135623731*IT_0046;
    const ccomplex_t IT_0048 = (0 + _Complex_I*1)*(IT_0044 + (-0.5)*IT_0047);
    const ccomplex_t IT_0049 = -IT_0048;
    const ccomplex_t IT_0050 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0051 = IT_0002*IT_0050;
    const ccomplex_t IT_0052 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0053 = IT_0002*IT_0052;
    const ccomplex_t IT_0054 = 1.4142135623731*IT_0053;
    const ccomplex_t IT_0055 = (0 + _Complex_I*1)*(IT_0051 + (-0.5)*IT_0054);
    const ccomplex_t IT_0056 = -IT_0055;
    const ccomplex_t IT_0057 = IT_0049*IT_0056;
    const ccomplex_t IT_0058 = IT_0042*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0061 = IT_0002*IT_0060;
    const ccomplex_t IT_0062 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0063 = IT_0002*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-0.5)*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0068 = IT_0002*IT_0067;
    const ccomplex_t IT_0069 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0070 = IT_0002*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + (-0.5)*IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = IT_0066*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0015 + IT_0041 + -m_N_2*
      (m_N_2 + (0 + _Complex_I*-1)*Gamma_o2) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0079 = IT_0002*IT_0078;
    const ccomplex_t IT_0080 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0081 = IT_0002*IT_0080;
    const ccomplex_t IT_0082 = 1.4142135623731*IT_0081;
    const ccomplex_t IT_0083 = (0 + _Complex_I*1)*(IT_0079 + (-0.5)*IT_0082);
    const ccomplex_t IT_0084 = -IT_0083;
    const ccomplex_t IT_0085 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0086 = IT_0002*IT_0085;
    const ccomplex_t IT_0087 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0088 = IT_0002*IT_0087;
    const ccomplex_t IT_0089 = 1.4142135623731*IT_0088;
    const ccomplex_t IT_0090 = (0 + _Complex_I*1)*(IT_0086 + (-0.5)*IT_0089);
    const ccomplex_t IT_0091 = -IT_0090;
    const ccomplex_t IT_0092 = IT_0084*IT_0091;
    const ccomplex_t IT_0093 = cpow((-2)*s_23 + IT_0015 + IT_0041 + -m_N_3*
      (m_N_3 + (0 + _Complex_I*-1)*Gamma_o3) + reg_prop, -1);
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*IT_0094;
    const ccomplex_t IT_0096 = cpow((-2)*s_23 + IT_0015 + IT_0041 + -m_N_4*
      (m_N_4 + (0 + _Complex_I*-1)*Gamma_o4) + reg_prop, -1);
    const ccomplex_t IT_0097 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0098 = IT_0002*IT_0097;
    const ccomplex_t IT_0099 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0100 = IT_0002*IT_0099;
    const ccomplex_t IT_0101 = 1.4142135623731*IT_0100;
    const ccomplex_t IT_0102 = (0 + _Complex_I*1)*(IT_0098 + (-0.5)*IT_0101);
    const ccomplex_t IT_0103 = -IT_0102;
    const ccomplex_t IT_0104 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0105 = IT_0002*IT_0104;
    const ccomplex_t IT_0106 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0107 = IT_0002*IT_0106;
    const ccomplex_t IT_0108 = 1.4142135623731*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + (-0.5)*IT_0108);
    const ccomplex_t IT_0110 = -IT_0109;
    const ccomplex_t IT_0111 = IT_0103*IT_0110;
    const ccomplex_t IT_0112 = IT_0096*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0032 + (-2)*IT_0059 + (-2)*IT_0077 + (-2)
      *IT_0095 + (-2)*IT_0113;
    const ccomplex_t IT_0115 = s_23*s_34*IT_0033;
    const ccomplex_t IT_0116 = -IT_0115;
    const ccomplex_t IT_0117 = s_24 + IT_0116;
    const ccomplex_t IT_0118 = IT_0034*IT_0117;
    const ccomplex_t IT_0119 = 2*IT_0118;
    const ccomplex_t IT_0120 = IT_0114*IT_0119;
    const ccomplex_t IT_0121 = IT_0040 + IT_0120;
    const ccomplex_t IT_0122 = cos(beta);
    const ccomplex_t IT_0123 = cos(alpha);
    const ccomplex_t IT_0124 = sin(beta);
    const ccomplex_t IT_0125 = sin(alpha);
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*(IT_0122
      *IT_0123 + IT_0124*IT_0125);
    const ccomplex_t IT_0127 = U_d2*V_Wp1;
    const ccomplex_t IT_0128 = IT_0123*IT_0127;
    const ccomplex_t IT_0129 = V_u1*U_Wm2;
    const ccomplex_t IT_0130 = IT_0125*IT_0129;
    const ccomplex_t IT_0131 = IT_0002*(IT_0128 + IT_0130);
    const ccomplex_t IT_0132 = 1.4142135623731*e_em*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*IT_0132;
    const ccomplex_t IT_0134 = 0.5*IT_0133;
    const ccomplex_t IT_0135 = IT_0126*IT_0134;
    const ccomplex_t IT_0136 = cpow((-2)*s_12 + -IT_0014 + -IT_0015 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0137 = IT_0135*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*IT_0137;
    const ccomplex_t IT_0139 = -IT_0003;
    const ccomplex_t IT_0140 = m_C_1*IT_0139;
    const ccomplex_t IT_0141 = IT_0018*IT_0140;
    const ccomplex_t IT_0142 = m_C_2*IT_0139;
    const ccomplex_t IT_0143 = IT_0031*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*m_W*e_em*IT_0002*(IT_0123
      *IT_0124 + -IT_0122*IT_0125);
    const ccomplex_t IT_0145 = IT_0123*IT_0129;
    const ccomplex_t IT_0146 = IT_0125*IT_0127;
    const ccomplex_t IT_0147 = IT_0002*(IT_0145 + -IT_0146);
    const ccomplex_t IT_0148 = 1.4142135623731*e_em*IT_0147;
    const ccomplex_t IT_0149 = 0.5*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*IT_0149;
    const ccomplex_t IT_0151 = IT_0144*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0014 + -IT_0015 + -reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0138 + -IT_0141 + -IT_0143 + IT_0154;
    const ccomplex_t IT_0156 = pow(s_13, 2);
    const ccomplex_t IT_0157 = IT_0033*IT_0156;
    const ccomplex_t IT_0158 = pow(s_14, 2);
    const ccomplex_t IT_0159 = IT_0033*IT_0158;
    const ccomplex_t IT_0160 = IT_0157 + IT_0159;
    const ccomplex_t IT_0161 = pow(m_W, -4);
    const ccomplex_t IT_0162 = s_13*s_14*s_34*IT_0161;
    const ccomplex_t IT_0163 = -IT_0014 + -IT_0162;
    const ccomplex_t IT_0164 = IT_0160 + IT_0163;
    const ccomplex_t IT_0165 = m_C_2*IT_0164;
    const ccomplex_t IT_0166 = 2*IT_0165;
    const ccomplex_t IT_0167 = IT_0155*IT_0166;
    const ccomplex_t IT_0168 = IT_0031*IT_0139;
    const ccomplex_t IT_0169 = IT_0059 + IT_0077 + IT_0095 + IT_0113;
    const ccomplex_t IT_0170 = s_12*s_23;
    const ccomplex_t IT_0171 = s_12*s_24*s_34*IT_0033;
    const ccomplex_t IT_0172 = (-0.5)*IT_0171;
    const ccomplex_t IT_0173 = s_14*s_23*s_24;
    const ccomplex_t IT_0174 = pow(s_24, 2);
    const ccomplex_t IT_0175 = s_13*IT_0174;
    const ccomplex_t IT_0176 = IT_0173 + IT_0175;
    const ccomplex_t IT_0177 = s_12*s_24*s_34;
    const ccomplex_t IT_0178 = s_13*IT_0015*IT_0041;
    const ccomplex_t IT_0179 = -IT_0177 + -IT_0178;
    const ccomplex_t IT_0180 = IT_0176 + IT_0179;
    const ccomplex_t IT_0181 = IT_0033*IT_0180;
    const ccomplex_t IT_0182 = (-0.25)*IT_0181;
    const ccomplex_t IT_0183 = s_13*IT_0015;
    const ccomplex_t IT_0184 = (-0.5)*IT_0183;
    const ccomplex_t IT_0185 = s_13*IT_0033*IT_0174;
    const ccomplex_t IT_0186 = 0.5*IT_0185;
    const ccomplex_t IT_0187 = s_14*s_23*s_24*IT_0033;
    const ccomplex_t IT_0188 = (-0.5)*IT_0187;
    const ccomplex_t IT_0189 = IT_0170 + IT_0172 + IT_0182 + IT_0184 + IT_0186
       + IT_0188;
    const ccomplex_t IT_0190 = s_23*IT_0033*IT_0158;
    const ccomplex_t IT_0191 = s_12*s_14*s_34*IT_0033;
    const ccomplex_t IT_0192 = s_13*s_14*s_24*IT_0033;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = s_23*IT_0014;
    const ccomplex_t IT_0195 = -IT_0194;
    const ccomplex_t IT_0196 = s_23*IT_0158;
    const ccomplex_t IT_0197 = s_13*s_14*s_24;
    const ccomplex_t IT_0198 = IT_0196 + IT_0197;
    const ccomplex_t IT_0199 = s_12*s_14*s_34;
    const ccomplex_t IT_0200 = s_23*IT_0014*IT_0041;
    const ccomplex_t IT_0201 = -IT_0199 + -IT_0200;
    const ccomplex_t IT_0202 = IT_0198 + IT_0201;
    const ccomplex_t IT_0203 = IT_0033*IT_0202;
    const ccomplex_t IT_0204 = 0.5*IT_0203;
    const ccomplex_t IT_0205 = IT_0190 + IT_0191 + IT_0193 + IT_0195 + IT_0204;
    const ccomplex_t IT_0206 = 4*IT_0205;
    const ccomplex_t IT_0207 = IT_0189 + (-0.125)*IT_0206;
    const ccomplex_t IT_0208 = s_13*s_23*IT_0158*IT_0161;
    const ccomplex_t IT_0209 = s_14*s_24*IT_0156*IT_0161;
    const ccomplex_t IT_0210 = s_12*IT_0014;
    const ccomplex_t IT_0211 = s_14*s_24*IT_0014*IT_0033;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = s_13*s_23*IT_0014*IT_0033;
    const ccomplex_t IT_0214 = -IT_0213;
    const ccomplex_t IT_0215 = s_12*s_13*s_14*s_34*IT_0161;
    const ccomplex_t IT_0216 = -IT_0215;
    const ccomplex_t IT_0217 = IT_0208 + IT_0209 + IT_0210 + IT_0212 + IT_0214
       + IT_0216;
    const ccomplex_t IT_0218 = 2*IT_0217;
    const ccomplex_t IT_0219 = s_13*s_14*s_23*s_24*IT_0161;
    const ccomplex_t IT_0220 = IT_0014*IT_0015;
    const ccomplex_t IT_0221 = IT_0156*IT_0161*IT_0174;
    const ccomplex_t IT_0222 = IT_0015*IT_0033*IT_0156;
    const ccomplex_t IT_0223 = -IT_0222;
    const ccomplex_t IT_0224 = IT_0014*IT_0033*IT_0174;
    const ccomplex_t IT_0225 = -IT_0224;
    const ccomplex_t IT_0226 = s_12*s_13*s_24*s_34*IT_0161;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0219 + IT_0220 + IT_0221 + IT_0223 + IT_0225
       + IT_0227;
    const ccomplex_t IT_0229 = 2*IT_0228;
    const ccomplex_t IT_0230 = IT_0218 + IT_0229;
    const ccomplex_t IT_0231 = pow(s_23, 2);
    const ccomplex_t IT_0232 = IT_0158*IT_0161*IT_0231;
    const ccomplex_t IT_0233 = IT_0015*IT_0033*IT_0158;
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0014*IT_0033*IT_0231;
    const ccomplex_t IT_0236 = -IT_0235;
    const ccomplex_t IT_0237 = s_12*s_14*s_23*s_34*IT_0161;
    const ccomplex_t IT_0238 = -IT_0237;
    const ccomplex_t IT_0239 = IT_0219 + IT_0220 + IT_0232 + IT_0234 + IT_0236
       + IT_0238;
    const ccomplex_t IT_0240 = 2*IT_0239;
    const ccomplex_t IT_0241 = s_14*s_24*IT_0161*IT_0231;
    const ccomplex_t IT_0242 = s_13*s_23*IT_0161*IT_0174;
    const ccomplex_t IT_0243 = s_12*IT_0015;
    const ccomplex_t IT_0244 = s_14*s_24*IT_0015*IT_0033;
    const ccomplex_t IT_0245 = -IT_0244;
    const ccomplex_t IT_0246 = s_13*s_23*IT_0015*IT_0033;
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = s_12*s_23*s_24*s_34*IT_0161;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0241 + IT_0242 + IT_0243 + IT_0245 + IT_0247
       + IT_0249;
    const ccomplex_t IT_0251 = 2*IT_0250;
    const ccomplex_t IT_0252 = IT_0240 + IT_0251;
    const ccomplex_t IT_0253 = IT_0020*IT_0031;
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = s_12*s_13;
    const ccomplex_t IT_0256 = s_23*IT_0033*IT_0156;
    const ccomplex_t IT_0257 = (-2)*IT_0256;
    const ccomplex_t IT_0258 = IT_0255 + IT_0257;
    const ccomplex_t IT_0259 = -IT_0190;
    const ccomplex_t IT_0260 = IT_0193 + IT_0258 + IT_0259;
    const ccomplex_t IT_0261 = s_13*s_14*s_23*s_34*IT_0161;
    const ccomplex_t IT_0262 = 2*IT_0261;
    const ccomplex_t IT_0263 = IT_0194 + IT_0260 + IT_0262;
    const ccomplex_t IT_0264 = 2*IT_0263;
    const ccomplex_t IT_0265 = s_13*IT_0033*IT_0231;
    const ccomplex_t IT_0266 = (-2)*IT_0265;
    const ccomplex_t IT_0267 = IT_0170 + IT_0266;
    const ccomplex_t IT_0268 = -IT_0185;
    const ccomplex_t IT_0269 = -IT_0187;
    const ccomplex_t IT_0270 = IT_0267 + IT_0268 + IT_0269;
    const ccomplex_t IT_0271 = s_13*s_23*s_24*s_34*IT_0161;
    const ccomplex_t IT_0272 = 2*IT_0271;
    const ccomplex_t IT_0273 = IT_0183 + IT_0270 + IT_0272;
    const ccomplex_t IT_0274 = 2*IT_0273;
    const ccomplex_t IT_0275 = IT_0264 + IT_0274;
    const ccomplex_t IT_0276 = IT_0003*IT_0031;
    const ccomplex_t IT_0277 = -IT_0276;
    const ccomplex_t IT_0278 = s_14*s_34*IT_0033;
    const ccomplex_t IT_0279 = s_12*IT_0041;
    const ccomplex_t IT_0280 = s_13*s_23;
    const ccomplex_t IT_0281 = (-2)*IT_0280;
    const ccomplex_t IT_0282 = IT_0279 + IT_0281;
    const ccomplex_t IT_0283 = IT_0033*IT_0282;
    const ccomplex_t IT_0284 = (s_13 + -IT_0278)*(s_12 + (-0.5)*IT_0283);
    const ccomplex_t IT_0285 = (-2)*IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0285;
    const ccomplex_t IT_0287 = s_24*s_34*IT_0033;
    const ccomplex_t IT_0288 = (s_12 + (-0.5)*IT_0283)*(s_23 + -IT_0287);
    const ccomplex_t IT_0289 = (-2)*IT_0288;
    const ccomplex_t IT_0290 = 2*IT_0289;
    const ccomplex_t IT_0291 = IT_0286 + IT_0290;
    const ccomplex_t IT_0292 = IT_0018*IT_0022;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = m_C_1*m_C_2;
    const ccomplex_t IT_0295 = IT_0164*IT_0294;
    const ccomplex_t IT_0296 = 2*IT_0295;
    const ccomplex_t IT_0297 = s_13*s_24*s_34*IT_0161;
    const ccomplex_t IT_0298 = s_13*s_23*IT_0033;
    const ccomplex_t IT_0299 = s_14*s_24*IT_0033;
    const ccomplex_t IT_0300 = -IT_0298 + -IT_0299;
    const ccomplex_t IT_0301 = s_12 + IT_0297 + IT_0300;
    const ccomplex_t IT_0302 = IT_0294*IT_0301;
    const ccomplex_t IT_0303 = (-2)*IT_0302;
    const ccomplex_t IT_0304 = IT_0296 + IT_0303;
    const ccomplex_t IT_0305 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0306 = IT_0002*IT_0305;
    const ccomplex_t IT_0307 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0308 = IT_0002*IT_0307;
    const ccomplex_t IT_0309 = 1.4142135623731*IT_0308;
    const ccomplex_t IT_0310 = (0 + _Complex_I*1)*(IT_0306 + 0.5*IT_0309);
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0313 = IT_0002*IT_0312;
    const ccomplex_t IT_0314 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0315 = IT_0002*IT_0314;
    const ccomplex_t IT_0316 = 1.4142135623731*IT_0315;
    const ccomplex_t IT_0317 = (0 + _Complex_I*1)*(IT_0313 + 0.5*IT_0316);
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = IT_0311*IT_0318;
    const ccomplex_t IT_0320 = IT_0093*IT_0319;
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*IT_0320;
    const ccomplex_t IT_0322 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0323 = IT_0002*IT_0322;
    const ccomplex_t IT_0324 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0325 = IT_0002*IT_0324;
    const ccomplex_t IT_0326 = 1.4142135623731*IT_0325;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*(IT_0323 + 0.5*IT_0326);
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0330 = IT_0002*IT_0329;
    const ccomplex_t IT_0331 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0332 = IT_0002*IT_0331;
    const ccomplex_t IT_0333 = 1.4142135623731*IT_0332;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*(IT_0330 + 0.5*IT_0333);
    const ccomplex_t IT_0335 = -IT_0334;
    const ccomplex_t IT_0336 = IT_0328*IT_0335;
    const ccomplex_t IT_0337 = IT_0096*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0340 = IT_0002*IT_0339;
    const ccomplex_t IT_0341 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0342 = IT_0002*IT_0341;
    const ccomplex_t IT_0343 = 1.4142135623731*IT_0342;
    const ccomplex_t IT_0344 = (0 + _Complex_I*1)*(IT_0340 + 0.5*IT_0343);
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0347 = IT_0002*IT_0346;
    const ccomplex_t IT_0348 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0349 = IT_0002*IT_0348;
    const ccomplex_t IT_0350 = 1.4142135623731*IT_0349;
    const ccomplex_t IT_0351 = (0 + _Complex_I*1)*(IT_0347 + 0.5*IT_0350);
    const ccomplex_t IT_0352 = -IT_0351;
    const ccomplex_t IT_0353 = IT_0345*IT_0352;
    const ccomplex_t IT_0354 = IT_0042*IT_0353;
    const ccomplex_t IT_0355 = (0 + _Complex_I*1)*IT_0354;
    const ccomplex_t IT_0356 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0357 = IT_0002*IT_0356;
    const ccomplex_t IT_0358 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0359 = IT_0002*IT_0358;
    const ccomplex_t IT_0360 = 1.4142135623731*IT_0359;
    const ccomplex_t IT_0361 = (0 + _Complex_I*1)*(IT_0357 + 0.5*IT_0360);
    const ccomplex_t IT_0362 = -IT_0361;
    const ccomplex_t IT_0363 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0364 = IT_0002*IT_0363;
    const ccomplex_t IT_0365 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0366 = IT_0002*IT_0365;
    const ccomplex_t IT_0367 = 1.4142135623731*IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*(IT_0364 + 0.5*IT_0367);
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = IT_0362*IT_0369;
    const ccomplex_t IT_0371 = IT_0075*IT_0370;
    const ccomplex_t IT_0372 = (0 + _Complex_I*1)*IT_0371;
    const ccomplex_t IT_0373 = -IT_0292 + (-2)*IT_0321 + (-2)*IT_0338 + (-2)
      *IT_0355 + (-2)*IT_0372;
    const ccomplex_t IT_0374 = IT_0033*IT_0231;
    const ccomplex_t IT_0375 = IT_0033*IT_0174;
    const ccomplex_t IT_0376 = IT_0374 + IT_0375;
    const ccomplex_t IT_0377 = s_23*s_24*s_34*IT_0161;
    const ccomplex_t IT_0378 = -IT_0015 + -IT_0377;
    const ccomplex_t IT_0379 = IT_0376 + IT_0378;
    const ccomplex_t IT_0380 = IT_0294*IT_0379;
    const ccomplex_t IT_0381 = 2*IT_0380;
    const ccomplex_t IT_0382 = s_14*s_23*s_34*IT_0161;
    const ccomplex_t IT_0383 = s_12 + IT_0300 + IT_0382;
    const ccomplex_t IT_0384 = IT_0294*IT_0383;
    const ccomplex_t IT_0385 = (-2)*IT_0384;
    const ccomplex_t IT_0386 = IT_0381 + IT_0385;
    const ccomplex_t IT_0387 = (0 + _Complex_I*1)*m_N_3;
    const ccomplex_t IT_0388 = IT_0084*IT_0318;
    const ccomplex_t IT_0389 = IT_0093*IT_0387*IT_0388;
    const ccomplex_t IT_0390 = (0 + _Complex_I*1)*m_N_4;
    const ccomplex_t IT_0391 = IT_0103*IT_0335;
    const ccomplex_t IT_0392 = IT_0096*IT_0390*IT_0391;
    const ccomplex_t IT_0393 = m_C_2*IT_0059;
    const ccomplex_t IT_0394 = (0 + _Complex_I*1)*m_N_1;
    const ccomplex_t IT_0395 = IT_0049*IT_0345;
    const ccomplex_t IT_0396 = IT_0042*IT_0394*IT_0395;
    const ccomplex_t IT_0397 = (0 + _Complex_I*1)*m_N_2;
    const ccomplex_t IT_0398 = IT_0066*IT_0369;
    const ccomplex_t IT_0399 = IT_0075*IT_0397*IT_0398;
    const ccomplex_t IT_0400 = m_C_2*IT_0077;
    const ccomplex_t IT_0401 = m_C_2*IT_0095;
    const ccomplex_t IT_0402 = m_C_2*IT_0113;
    const ccomplex_t IT_0403 = IT_0389 + IT_0392 + -IT_0393 + IT_0396 +
       IT_0399 + -IT_0400 + -IT_0401 + -IT_0402;
    const ccomplex_t IT_0404 = s_14*s_24;
    const ccomplex_t IT_0405 = -IT_0404;
    const ccomplex_t IT_0406 = IT_0279 + IT_0405;
    const ccomplex_t IT_0407 = IT_0033*IT_0406;
    const ccomplex_t IT_0408 = -IT_0299 + (-0.25)*IT_0407;
    const ccomplex_t IT_0409 = s_12 + IT_0408;
    const ccomplex_t IT_0410 = m_C_2*IT_0409;
    const ccomplex_t IT_0411 = (-4)*IT_0410;
    const ccomplex_t IT_0412 = 2*IT_0411;
    const ccomplex_t IT_0413 = IT_0014*IT_0041;
    const ccomplex_t IT_0414 = -IT_0413;
    const ccomplex_t IT_0415 = IT_0158 + IT_0414;
    const ccomplex_t IT_0416 = IT_0033*IT_0415;
    const ccomplex_t IT_0417 = -IT_0014;
    const ccomplex_t IT_0418 = IT_0159 + IT_0417;
    const ccomplex_t IT_0419 = (-4)*IT_0418;
    const ccomplex_t IT_0420 = IT_0416 + IT_0419;
    const ccomplex_t IT_0421 = m_C_2*IT_0420;
    const ccomplex_t IT_0422 = (-2)*IT_0421;
    const ccomplex_t IT_0423 = IT_0412 + IT_0422;
    const ccomplex_t IT_0424 = IT_0110*IT_0328;
    const ccomplex_t IT_0425 = IT_0096*IT_0390*IT_0424;
    const ccomplex_t IT_0426 = m_C_2*IT_0338;
    const ccomplex_t IT_0427 = IT_0056*IT_0352;
    const ccomplex_t IT_0428 = IT_0042*IT_0394*IT_0427;
    const ccomplex_t IT_0429 = IT_0073*IT_0362;
    const ccomplex_t IT_0430 = IT_0075*IT_0397*IT_0429;
    const ccomplex_t IT_0431 = IT_0091*IT_0311;
    const ccomplex_t IT_0432 = IT_0093*IT_0387*IT_0431;
    const ccomplex_t IT_0433 = m_C_2*IT_0355;
    const ccomplex_t IT_0434 = m_C_2*IT_0372;
    const ccomplex_t IT_0435 = m_C_2*IT_0321;
    const ccomplex_t IT_0436 = IT_0425 + -IT_0426 + IT_0428 + IT_0430 +
       IT_0432 + -IT_0433 + -IT_0434 + -IT_0435;
    const ccomplex_t IT_0437 = s_14*s_23*s_34*IT_0033;
    const ccomplex_t IT_0438 = (-2)*IT_0280 + -IT_0404 + 2*IT_0437;
    const ccomplex_t IT_0439 = IT_0279 + IT_0438;
    const ccomplex_t IT_0440 = IT_0033*IT_0439;
    const ccomplex_t IT_0441 = -IT_0299 + (-0.5)*IT_0440;
    const ccomplex_t IT_0442 = s_12 + IT_0441;
    const ccomplex_t IT_0443 = m_C_1*IT_0442;
    const ccomplex_t IT_0444 = (-2)*IT_0443;
    const ccomplex_t IT_0445 = 2*IT_0444;
    const ccomplex_t IT_0446 = s_23*s_24*s_34*IT_0033;
    const ccomplex_t IT_0447 = IT_0015*IT_0041;
    const ccomplex_t IT_0448 = 0.5*IT_0174 + -IT_0446 + (-0.5)*IT_0447;
    const ccomplex_t IT_0449 = IT_0231 + IT_0448;
    const ccomplex_t IT_0450 = IT_0033*IT_0449;
    const ccomplex_t IT_0451 = -IT_0375;
    const ccomplex_t IT_0452 = IT_0015 + IT_0450 + IT_0451;
    const ccomplex_t IT_0453 = m_C_1*IT_0452;
    const ccomplex_t IT_0454 = (-4)*IT_0453;
    const ccomplex_t IT_0455 = IT_0445 + IT_0454;
    const ccomplex_t IT_0456 = IT_0321 + IT_0338 + IT_0355 + IT_0372;
    const ccomplex_t IT_0457 = s_13*IT_0041;
    const ccomplex_t IT_0458 = s_14*s_34;
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = IT_0457 + IT_0459;
    const ccomplex_t IT_0461 = IT_0033*IT_0460;
    const ccomplex_t IT_0462 = -IT_0278 + (-0.25)*IT_0461;
    const ccomplex_t IT_0463 = s_13 + IT_0462;
    const ccomplex_t IT_0464 = IT_0294*IT_0463;
    const ccomplex_t IT_0465 = (-4)*IT_0464;
    const ccomplex_t IT_0466 = 2*IT_0465;
    const ccomplex_t IT_0467 = s_23*IT_0041;
    const ccomplex_t IT_0468 = s_24*s_34;
    const ccomplex_t IT_0469 = -IT_0468;
    const ccomplex_t IT_0470 = IT_0467 + IT_0469;
    const ccomplex_t IT_0471 = IT_0033*IT_0470;
    const ccomplex_t IT_0472 = -IT_0287 + (-0.25)*IT_0471;
    const ccomplex_t IT_0473 = s_23 + IT_0472;
    const ccomplex_t IT_0474 = IT_0294*IT_0473;
    const ccomplex_t IT_0475 = (-4)*IT_0474;
    const ccomplex_t IT_0476 = 2*IT_0475;
    const ccomplex_t IT_0477 = IT_0466 + IT_0476;
    const ccomplex_t IT_0478 = IT_0018*IT_0139;
    const ccomplex_t IT_0479 = -IT_0478;
    const ccomplex_t IT_0480 = IT_0294*IT_0418;
    const ccomplex_t IT_0481 = 6*IT_0480;
    const ccomplex_t IT_0482 = -IT_0299;
    const ccomplex_t IT_0483 = s_12 + IT_0482;
    const ccomplex_t IT_0484 = IT_0294*IT_0483;
    const ccomplex_t IT_0485 = (-6)*IT_0484;
    const ccomplex_t IT_0486 = IT_0481 + IT_0485;
    const ccomplex_t IT_0487 = -IT_0015;
    const ccomplex_t IT_0488 = IT_0375 + IT_0487;
    const ccomplex_t IT_0489 = IT_0294*IT_0488;
    const ccomplex_t IT_0490 = 6*IT_0489;
    const ccomplex_t IT_0491 = IT_0485 + IT_0490;
    const ccomplex_t IT_0492 = -IT_0278;
    const ccomplex_t IT_0493 = s_13 + IT_0492;
    const ccomplex_t IT_0494 = IT_0294*IT_0493;
    const ccomplex_t IT_0495 = (-6)*IT_0494;
    const ccomplex_t IT_0496 = -IT_0287;
    const ccomplex_t IT_0497 = s_23 + IT_0496;
    const ccomplex_t IT_0498 = IT_0294*IT_0497;
    const ccomplex_t IT_0499 = (-6)*IT_0498;
    const ccomplex_t IT_0500 = IT_0495 + IT_0499;
    const ccomplex_t IT_0501 = m_C_2*IT_0301;
    const ccomplex_t IT_0502 = (-2)*IT_0501;
    const ccomplex_t IT_0503 = conj(V_u2)*conj(U_Wm1);
    const ccomplex_t IT_0504 = IT_0123*IT_0503;
    const ccomplex_t IT_0505 = conj(U_d1)*conj(V_Wp2);
    const ccomplex_t IT_0506 = IT_0125*IT_0505;
    const ccomplex_t IT_0507 = IT_0002*(IT_0504 + -IT_0506);
    const ccomplex_t IT_0508 = 1.4142135623731*e_em*IT_0507;
    const ccomplex_t IT_0509 = 0.5*IT_0508;
    const ccomplex_t IT_0510 = (0 + _Complex_I*1)*IT_0509;
    const ccomplex_t IT_0511 = IT_0144*IT_0510;
    const ccomplex_t IT_0512 = IT_0152*IT_0511;
    const ccomplex_t IT_0513 = (0 + _Complex_I*1)*IT_0512;
    const ccomplex_t IT_0514 = IT_0123*IT_0505;
    const ccomplex_t IT_0515 = IT_0125*IT_0503;
    const ccomplex_t IT_0516 = IT_0002*(IT_0514 + IT_0515);
    const ccomplex_t IT_0517 = 1.4142135623731*e_em*IT_0516;
    const ccomplex_t IT_0518 = (0 + _Complex_I*1)*IT_0517;
    const ccomplex_t IT_0519 = 0.5*IT_0518;
    const ccomplex_t IT_0520 = IT_0126*IT_0519;
    const ccomplex_t IT_0521 = IT_0136*IT_0520;
    const ccomplex_t IT_0522 = (0 + _Complex_I*1)*IT_0521;
    const ccomplex_t IT_0523 = IT_0018*IT_0142;
    const ccomplex_t IT_0524 = IT_0031*IT_0140;
    const ccomplex_t IT_0525 = IT_0513 + IT_0522 + IT_0523 + IT_0524;
    const ccomplex_t IT_0526 = m_C_1*IT_0383;
    const ccomplex_t IT_0527 = 2*IT_0526;
    const ccomplex_t IT_0528 = m_C_1*IT_0379;
    const ccomplex_t IT_0529 = (-2)*IT_0528;
    const ccomplex_t IT_0530 = (-0.125)*IT_0167 + IT_0169*IT_0207 + (-0.125)
      *IT_0032*IT_0230 + (-0.125)*IT_0114*IT_0252 + (-0.125)*IT_0254*IT_0275 + (
      -0.125)*IT_0277*IT_0291 + (-0.125)*IT_0293*IT_0304 + (-0.125)*IT_0373
      *IT_0386 + (-0.125)*IT_0403*IT_0423 + (-0.125)*IT_0436*IT_0455 + (-0.125)
      *IT_0456*IT_0477 + (-0.125)*IT_0479*IT_0486 + (-0.125)*IT_0479*IT_0491 + (
      -0.125)*IT_0019*IT_0500 + (-0.125)*IT_0155*IT_0502 + (-0.125)*IT_0525*
      (IT_0527 + IT_0529);
    const ccomplex_t IT_0531 = (s_12 + (-0.5)*IT_0283)*(s_12 + -IT_0299);
    const ccomplex_t IT_0532 = (-2)*IT_0531;
    const ccomplex_t IT_0533 = 2*IT_0532;
    const ccomplex_t IT_0534 = conj(IT_0155)*IT_0529;
    const ccomplex_t IT_0535 = -(s_12 + (-0.5)*IT_0283)*(IT_0015 + -IT_0375);
    const ccomplex_t IT_0536 = (-2)*IT_0535;
    const ccomplex_t IT_0537 = (-2)*IT_0536;
    const ccomplex_t IT_0538 = (-0.125)*conj(IT_0021)*IT_0274 + (-0.125)*conj
      (IT_0019)*IT_0290 + (-0.125)*IT_0229*conj(IT_0293) + (-0.125)*conj(IT_0032
      )*IT_0303 + (-0.125)*IT_0251*conj(IT_0373) + (-0.125)*conj(IT_0114)
      *IT_0381 + (-0.125)*IT_0412*conj(IT_0436) + (-0.125)*conj(IT_0403)*IT_0454
       + IT_0189*conj(IT_0456) + (-0.125)*conj(IT_0169)*IT_0476 + (-0.125)*conj
      (IT_0277)*IT_0499 + (-0.125)*IT_0502*conj(IT_0525) + (-0.125)*conj(IT_0479
      )*IT_0533 + (-0.125)*IT_0534 + (-0.125)*conj(IT_0479)*IT_0537;
    const ccomplex_t IT_0539 = -(IT_0014 + -IT_0159)*(s_12 + (-0.5)*IT_0283);
    const ccomplex_t IT_0540 = (-2)*IT_0539;
    const ccomplex_t IT_0541 = (-2)*IT_0540;
    const ccomplex_t IT_0542 = conj(IT_0021)*IT_0264 + conj(IT_0019)*IT_0286 +
       IT_0218*conj(IT_0293) + conj(IT_0032)*IT_0296 + IT_0240*conj(IT_0373) +
       conj(IT_0114)*IT_0385 + IT_0422*conj(IT_0436) + conj(IT_0403)*IT_0445 +
       IT_0206*conj(IT_0456) + conj(IT_0169)*IT_0466 + conj(IT_0277)*IT_0495 +
       IT_0166*conj(IT_0525) + conj(IT_0155)*IT_0527 + conj(IT_0479)*IT_0533 +
       conj(IT_0479)*IT_0541;
    const ccomplex_t IT_0543 = s_24*s_34*IT_0156*IT_0161;
    const ccomplex_t IT_0544 = pow(s_34, 2);
    const ccomplex_t IT_0545 = s_12*s_13*IT_0161*IT_0544;
    const ccomplex_t IT_0546 = -IT_0545;
    const ccomplex_t IT_0547 = s_24*s_34*IT_0014*IT_0033;
    const ccomplex_t IT_0548 = -IT_0547;
    const ccomplex_t IT_0549 = IT_0194 + IT_0258 + IT_0261 + IT_0543 + IT_0546
       + IT_0548;
    const ccomplex_t IT_0550 = 2*IT_0549;
    const ccomplex_t IT_0551 = s_14*s_34*IT_0161*IT_0231;
    const ccomplex_t IT_0552 = s_12*s_23*IT_0161*IT_0544;
    const ccomplex_t IT_0553 = -IT_0552;
    const ccomplex_t IT_0554 = s_14*s_34*IT_0015*IT_0033;
    const ccomplex_t IT_0555 = -IT_0554;
    const ccomplex_t IT_0556 = IT_0183 + IT_0267 + IT_0271 + IT_0551 + IT_0553
       + IT_0555;
    const ccomplex_t IT_0557 = 2*IT_0556;
    const ccomplex_t IT_0558 = s_34*m_C_2*IT_0033;
    const ccomplex_t IT_0559 = IT_0037*IT_0558;
    const ccomplex_t IT_0560 = 2*IT_0559;
    const ccomplex_t IT_0561 = IT_0033*IT_0544;
    const ccomplex_t IT_0562 = (s_12 + (-0.5)*IT_0283)*(IT_0041 + -IT_0561);
    const ccomplex_t IT_0563 = (-2)*IT_0562;
    const ccomplex_t IT_0564 = 2*IT_0563;
    const ccomplex_t IT_0565 = conj(IT_0032)*IT_0039;
    const ccomplex_t IT_0566 = conj(IT_0114)*IT_0119;
    const ccomplex_t IT_0567 = m_C_2*IT_0463;
    const ccomplex_t IT_0568 = (-4)*IT_0567;
    const ccomplex_t IT_0569 = 2*IT_0568;
    const ccomplex_t IT_0570 = -IT_0544;
    const ccomplex_t IT_0571 = pow(m_W, 4);
    const ccomplex_t IT_0572 = IT_0570 + IT_0571;
    const ccomplex_t IT_0573 = IT_0033*IT_0572;
    const ccomplex_t IT_0574 = -IT_0041 + 0.25*IT_0573;
    const ccomplex_t IT_0575 = IT_0561 + IT_0574;
    const ccomplex_t IT_0576 = IT_0294*IT_0575;
    const ccomplex_t IT_0577 = 4*IT_0576;
    const ccomplex_t IT_0578 = 2*IT_0577;
    const ccomplex_t IT_0579 = -IT_0561;
    const ccomplex_t IT_0580 = IT_0041 + IT_0579;
    const ccomplex_t IT_0581 = IT_0294*IT_0580;
    const ccomplex_t IT_0582 = (-6)*IT_0581;
    const ccomplex_t IT_0583 = s_13*s_24*s_34*IT_0033;
    const ccomplex_t IT_0584 = s_12*IT_0033*IT_0544;
    const ccomplex_t IT_0585 = -IT_0584;
    const ccomplex_t IT_0586 = s_12*IT_0571;
    const ccomplex_t IT_0587 = s_14*s_23*s_34;
    const ccomplex_t IT_0588 = s_13*s_24*s_34;
    const ccomplex_t IT_0589 = IT_0586 + IT_0587 + IT_0588;
    const ccomplex_t IT_0590 = s_13*s_23*IT_0041;
    const ccomplex_t IT_0591 = s_12*IT_0544;
    const ccomplex_t IT_0592 = (-2)*IT_0590 + -IT_0591;
    const ccomplex_t IT_0593 = IT_0589 + IT_0592;
    const ccomplex_t IT_0594 = IT_0033*IT_0593;
    const ccomplex_t IT_0595 = (-0.5)*IT_0594;
    const ccomplex_t IT_0596 = -IT_0437;
    const ccomplex_t IT_0597 = IT_0279 + IT_0583 + IT_0585 + IT_0595 + IT_0596;
    const ccomplex_t IT_0598 = (-4)*IT_0597;
    const ccomplex_t IT_0599 = s_34*m_C_1*IT_0033;
    const ccomplex_t IT_0600 = IT_0117*IT_0599;
    const ccomplex_t IT_0601 = (-2)*IT_0600;
    const ccomplex_t IT_0602 = s_34*IT_0033;
    const ccomplex_t IT_0603 = s_14*s_23;
    const ccomplex_t IT_0604 = s_13*s_24;
    const ccomplex_t IT_0605 = IT_0603 + IT_0604;
    const ccomplex_t IT_0606 = s_13*s_23*s_34*IT_0033;
    const ccomplex_t IT_0607 = (-2)*IT_0606;
    const ccomplex_t IT_0608 = IT_0605 + IT_0607;
    const ccomplex_t IT_0609 = IT_0602*IT_0608;
    const ccomplex_t IT_0610 = (-2)*IT_0609;
    const ccomplex_t IT_0611 = IT_0467 + IT_0468;
    const ccomplex_t IT_0612 = s_23*IT_0033*IT_0544;
    const ccomplex_t IT_0613 = (-2)*IT_0612;
    const ccomplex_t IT_0614 = IT_0611 + IT_0613;
    const ccomplex_t IT_0615 = IT_0033*IT_0614;
    const ccomplex_t IT_0616 = -IT_0287 + 0.5*IT_0615;
    const ccomplex_t IT_0617 = s_23 + IT_0616;
    const ccomplex_t IT_0618 = m_C_1*IT_0617;
    const ccomplex_t IT_0619 = 2*IT_0618;
    const ccomplex_t IT_0620 = (-2)*IT_0619;
    const ccomplex_t IT_0621 = IT_0286*conj(IT_0479) + IT_0290*conj(IT_0479) +
       conj(IT_0293)*IT_0550 + conj(IT_0373)*IT_0557 + conj(IT_0525)*IT_0560 +
       conj(IT_0019)*IT_0564 + IT_0565 + IT_0566 + conj(IT_0436)*IT_0569 + conj
      (IT_0169)*IT_0578 + conj(IT_0277)*IT_0582 + conj(IT_0456)*IT_0598 + conj
      (IT_0155)*IT_0601 + conj(IT_0021)*IT_0610 + conj(IT_0403)*IT_0620;
    const ccomplex_t IT_0622 = IT_0039*conj(IT_0293) + IT_0119*conj(IT_0373) +
       conj(IT_0479)*IT_0495 + conj(IT_0479)*IT_0499 + conj(IT_0032)*IT_0550 +
       conj(IT_0114)*IT_0557 + conj(IT_0155)*IT_0560 + conj(IT_0277)*IT_0564 +
       conj(IT_0403)*IT_0569 + conj(IT_0456)*IT_0578 + conj(IT_0019)*IT_0582 +
       conj(IT_0169)*IT_0598 + conj(IT_0525)*IT_0601 + conj(IT_0254)*IT_0610 +
       conj(IT_0436)*IT_0620;
    const ccomplex_t IT_0623 = IT_0161*IT_0544;
    const ccomplex_t IT_0624 = 0.5*IT_0623;
    const ccomplex_t IT_0625 = 1 + IT_0624;
    const ccomplex_t IT_0626 = IT_0294*IT_0625;
    const ccomplex_t IT_0627 = (-4)*IT_0626;
    const ccomplex_t IT_0628 = s_13*s_14*s_34*IT_0033;
    const ccomplex_t IT_0629 = 2*IT_0158 + -IT_0413 + (-2)*IT_0628;
    const ccomplex_t IT_0630 = IT_0156 + IT_0629;
    const ccomplex_t IT_0631 = IT_0033*IT_0630;
    const ccomplex_t IT_0632 = -IT_0014 + (-0.5)*IT_0631;
    const ccomplex_t IT_0633 = IT_0157 + IT_0632;
    const ccomplex_t IT_0634 = m_C_2*IT_0633;
    const ccomplex_t IT_0635 = (-2)*IT_0634;
    const ccomplex_t IT_0636 = 2*IT_0635;
    const ccomplex_t IT_0637 = -IT_0447;
    const ccomplex_t IT_0638 = IT_0231 + IT_0637;
    const ccomplex_t IT_0639 = IT_0033*IT_0638;
    const ccomplex_t IT_0640 = IT_0374 + IT_0487;
    const ccomplex_t IT_0641 = (-4)*IT_0640;
    const ccomplex_t IT_0642 = IT_0639 + IT_0641;
    const ccomplex_t IT_0643 = m_C_1*IT_0642;
    const ccomplex_t IT_0644 = 2*IT_0643;
    const ccomplex_t IT_0645 = (-2)*IT_0458;
    const ccomplex_t IT_0646 = IT_0457 + IT_0645;
    const ccomplex_t IT_0647 = IT_0033*IT_0646;
    const ccomplex_t IT_0648 = s_13*IT_0571;
    const ccomplex_t IT_0649 = s_14*s_34*IT_0041;
    const ccomplex_t IT_0650 = (-2)*IT_0649;
    const ccomplex_t IT_0651 = IT_0648 + IT_0650;
    const ccomplex_t IT_0652 = IT_0161*IT_0651;
    const ccomplex_t IT_0653 = (-0.666666666666667)*IT_0647 +
       0.166666666666667*IT_0652;
    const ccomplex_t IT_0654 = s_13 + IT_0653;
    const ccomplex_t IT_0655 = m_C_2*IT_0654;
    const ccomplex_t IT_0656 = 6*IT_0655;
    const ccomplex_t IT_0657 = 2*IT_0656;
    const ccomplex_t IT_0658 = s_23*IT_0571;
    const ccomplex_t IT_0659 = IT_0161*IT_0658;
    const ccomplex_t IT_0660 = (-10)*s_23;
    const ccomplex_t IT_0661 = IT_0659 + IT_0660;
    const ccomplex_t IT_0662 = m_C_1*IT_0661;
    const ccomplex_t IT_0663 = 2*IT_0662;
    const ccomplex_t IT_0664 = IT_0033*IT_0279;
    const ccomplex_t IT_0665 = (-0.5)*IT_0664;
    const ccomplex_t IT_0666 = s_34*IT_0161;
    const ccomplex_t IT_0667 = s_12*s_34;
    const ccomplex_t IT_0668 = IT_0603 + IT_0667;
    const ccomplex_t IT_0669 = -IT_0604;
    const ccomplex_t IT_0670 = IT_0668 + IT_0669;
    const ccomplex_t IT_0671 = IT_0666*IT_0670;
    const ccomplex_t IT_0672 = 0.25*IT_0671;
    const ccomplex_t IT_0673 = s_12 + IT_0665 + IT_0672;
    const ccomplex_t IT_0674 = 8*IT_0673;
    const ccomplex_t IT_0675 = (-2)*IT_0404;
    const ccomplex_t IT_0676 = IT_0279 + IT_0675;
    const ccomplex_t IT_0677 = IT_0033*IT_0676;
    const ccomplex_t IT_0678 = s_14*s_24*IT_0041;
    const ccomplex_t IT_0679 = 4*IT_0587 + (-2)*IT_0590 + (-2)*IT_0678;
    const ccomplex_t IT_0680 = IT_0586 + IT_0679;
    const ccomplex_t IT_0681 = IT_0161*IT_0680;
    const ccomplex_t IT_0682 = 8*s_12 + (-4)*IT_0283 + (-4)*IT_0677 + 2*IT_0681;
    const ccomplex_t IT_0683 = 1.125*IT_0294;
    const ccomplex_t IT_0684 = (-16)*IT_0683;
    const ccomplex_t IT_0685 = -IT_0280 + (-2)*IT_0404 + 2*IT_0437;
    const ccomplex_t IT_0686 = IT_0279 + IT_0685;
    const ccomplex_t IT_0687 = IT_0033*IT_0686;
    const ccomplex_t IT_0688 = -IT_0298 + (-0.5)*IT_0687;
    const ccomplex_t IT_0689 = s_12 + IT_0688;
    const ccomplex_t IT_0690 = m_C_2*IT_0689;
    const ccomplex_t IT_0691 = (-2)*IT_0690;
    const ccomplex_t IT_0692 = (-2)*IT_0691;
    const ccomplex_t IT_0693 = -IT_0280;
    const ccomplex_t IT_0694 = IT_0279 + IT_0693;
    const ccomplex_t IT_0695 = IT_0033*IT_0694;
    const ccomplex_t IT_0696 = -IT_0298 + (-0.25)*IT_0695;
    const ccomplex_t IT_0697 = s_12 + IT_0696;
    const ccomplex_t IT_0698 = m_C_1*IT_0697;
    const ccomplex_t IT_0699 = (-4)*IT_0698;
    const ccomplex_t IT_0700 = (-2)*IT_0699;
    const ccomplex_t IT_0701 = IT_0033*IT_0457;
    const ccomplex_t IT_0702 = s_14*IT_0041;
    const ccomplex_t IT_0703 = IT_0666*IT_0702;
    const ccomplex_t IT_0704 = (-0.5)*IT_0701 + 0.25*IT_0703;
    const ccomplex_t IT_0705 = s_13 + IT_0704;
    const ccomplex_t IT_0706 = m_C_2*IT_0705;
    const ccomplex_t IT_0707 = 4*IT_0706;
    const ccomplex_t IT_0708 = (-2)*IT_0707;
    const ccomplex_t IT_0709 = IT_0033*IT_0467;
    const ccomplex_t IT_0710 = s_24*IT_0041;
    const ccomplex_t IT_0711 = s_23*s_34;
    const ccomplex_t IT_0712 = (-2)*IT_0711;
    const ccomplex_t IT_0713 = IT_0710 + IT_0712;
    const ccomplex_t IT_0714 = IT_0666*IT_0713;
    const ccomplex_t IT_0715 = (-0.5)*IT_0709 + (-0.25)*IT_0714;
    const ccomplex_t IT_0716 = s_23 + IT_0715;
    const ccomplex_t IT_0717 = m_C_1*IT_0716;
    const ccomplex_t IT_0718 = (-4)*IT_0717;
    const ccomplex_t IT_0719 = (-2)*IT_0718;
    const ccomplex_t IT_0720 = IT_0412*conj(IT_0479) + IT_0422*conj(IT_0479) +
       conj(IT_0019)*IT_0569 + conj(IT_0277)*IT_0620 + conj(IT_0155)*IT_0627 +
       conj(IT_0293)*IT_0636 + conj(IT_0114)*IT_0644 + conj(IT_0456)*IT_0657 +
       conj(IT_0169)*IT_0663 + conj(IT_0525)*IT_0674 + conj(IT_0436)*IT_0682 +
       conj(IT_0403)*IT_0684 + conj(IT_0373)*IT_0692 + conj(IT_0032)*IT_0700 +
       conj(IT_0021)*IT_0708 + conj(IT_0254)*IT_0719;
    const ccomplex_t IT_0721 = IT_0445*conj(IT_0479) + IT_0454*conj(IT_0479) +
       conj(IT_0277)*IT_0569 + conj(IT_0019)*IT_0620 + conj(IT_0525)*IT_0627 +
       conj(IT_0032)*IT_0636 + conj(IT_0373)*IT_0644 + conj(IT_0169)*IT_0657 +
       conj(IT_0456)*IT_0663 + conj(IT_0155)*IT_0674 + conj(IT_0403)*IT_0682 +
       conj(IT_0436)*IT_0684 + conj(IT_0114)*IT_0692 + conj(IT_0293)*IT_0700 +
       conj(IT_0254)*IT_0708 + conj(IT_0021)*IT_0719;
    const ccomplex_t IT_0722 = (-0.5)*IT_0647 + 0.5*IT_0701 + (-0.5)*IT_0703;
    const ccomplex_t IT_0723 = s_13 + IT_0722;
    const ccomplex_t IT_0724 = m_C_2*IT_0723;
    const ccomplex_t IT_0725 = 2*IT_0724;
    const ccomplex_t IT_0726 = 2*IT_0725;
    const ccomplex_t IT_0727 = s_34*IT_0041;
    const ccomplex_t IT_0728 = IT_0666*IT_0727;
    const ccomplex_t IT_0729 = (-2)*IT_0544;
    const ccomplex_t IT_0730 = IT_0571 + IT_0729;
    const ccomplex_t IT_0731 = IT_0033*IT_0730;
    const ccomplex_t IT_0732 = (-0.333333333333333)*IT_0728 + (
      -0.333333333333333)*IT_0731;
    const ccomplex_t IT_0733 = IT_0041 + IT_0732;
    const ccomplex_t IT_0734 = IT_0294*IT_0733;
    const ccomplex_t IT_0735 = 3*IT_0734;
    const ccomplex_t IT_0736 = 2*IT_0735;
    const ccomplex_t IT_0737 = s_12*s_23*IT_0041;
    const ccomplex_t IT_0738 = s_14*s_34*IT_0015;
    const ccomplex_t IT_0739 = 0.5*IT_0178 + (-0.5)*IT_0737 + -IT_0738;
    const ccomplex_t IT_0740 = IT_0173 + IT_0739;
    const ccomplex_t IT_0741 = IT_0033*IT_0740;
    const ccomplex_t IT_0742 = (-2)*IT_0741;
    const ccomplex_t IT_0743 = (-0.5)*IT_0742;
    const ccomplex_t IT_0744 = -IT_0183;
    const ccomplex_t IT_0745 = IT_0170 + IT_0743 + IT_0744;
    const ccomplex_t IT_0746 = 4*IT_0745;
    const ccomplex_t IT_0747 = (-2)*IT_0590;
    const ccomplex_t IT_0748 = IT_0586 + IT_0747;
    const ccomplex_t IT_0749 = IT_0033*IT_0748;
    const ccomplex_t IT_0750 = 0.5*IT_0749;
    const ccomplex_t IT_0751 = 2*IT_0587 + 2*IT_0588 + (-2)*IT_0590 + (-2)
      *IT_0591;
    const ccomplex_t IT_0752 = IT_0586 + IT_0751;
    const ccomplex_t IT_0753 = IT_0033*IT_0752;
    const ccomplex_t IT_0754 = (-0.5)*IT_0753;
    const ccomplex_t IT_0755 = s_14*s_23*IT_0041;
    const ccomplex_t IT_0756 = s_13*s_24*IT_0041;
    const ccomplex_t IT_0757 = IT_0755 + IT_0756;
    const ccomplex_t IT_0758 = s_12*s_34*IT_0041;
    const ccomplex_t IT_0759 = -IT_0758;
    const ccomplex_t IT_0760 = IT_0757 + IT_0759;
    const ccomplex_t IT_0761 = IT_0666*IT_0760;
    const ccomplex_t IT_0762 = 0.5*IT_0761;
    const ccomplex_t IT_0763 = IT_0279 + IT_0281 + IT_0750 + IT_0754 + IT_0762;
    const ccomplex_t IT_0764 = 4*IT_0763;
    const ccomplex_t IT_0765 = s_12*s_13*IT_0041;
    const ccomplex_t IT_0766 = -IT_0199 + (-0.5)*IT_0200 + 0.5*IT_0765;
    const ccomplex_t IT_0767 = IT_0196 + IT_0766;
    const ccomplex_t IT_0768 = IT_0033*IT_0767;
    const ccomplex_t IT_0769 = 2*IT_0768;
    const ccomplex_t IT_0770 = (-0.5)*IT_0769;
    const ccomplex_t IT_0771 = IT_0195 + IT_0255 + IT_0770;
    const ccomplex_t IT_0772 = (-4)*IT_0771;
    const ccomplex_t IT_0773 = pow(m_W, 6);
    const ccomplex_t IT_0774 = s_12*IT_0773;
    const ccomplex_t IT_0775 = s_14*s_24*IT_0571;
    const ccomplex_t IT_0776 = (-2)*IT_0775;
    const ccomplex_t IT_0777 = IT_0774 + IT_0776;
    const ccomplex_t IT_0778 = IT_0161*IT_0777;
    const ccomplex_t IT_0779 = IT_0033*IT_0680;
    const ccomplex_t IT_0780 = (-4)*IT_0279 + 16*IT_0280 + (-2)*IT_0778 + 4
      *IT_0779;
    const ccomplex_t IT_0781 = IT_0666*IT_0710;
    const ccomplex_t IT_0782 = (-2)*IT_0468;
    const ccomplex_t IT_0783 = IT_0467 + IT_0782;
    const ccomplex_t IT_0784 = IT_0033*IT_0783;
    const ccomplex_t IT_0785 = 0.5*IT_0709 + (-0.5)*IT_0781 + (-0.5)*IT_0784;
    const ccomplex_t IT_0786 = s_23 + IT_0785;
    const ccomplex_t IT_0787 = m_C_1*IT_0786;
    const ccomplex_t IT_0788 = 2*IT_0787;
    const ccomplex_t IT_0789 = (-2)*IT_0788;
    const ccomplex_t IT_0790 = IT_0041*IT_0294;
    const ccomplex_t IT_0791 = 9*IT_0790;
    const ccomplex_t IT_0792 = (-2)*IT_0791;
    const ccomplex_t IT_0793 = IT_0189*conj(IT_0479) + (-0.125)*IT_0206*conj
      (IT_0479) + (-0.125)*conj(IT_0277)*IT_0578 + (-0.125)*conj(IT_0019)
      *IT_0598 + (-0.125)*conj(IT_0436)*IT_0657 + (-0.125)*conj(IT_0403)*IT_0663
       + (-0.125)*conj(IT_0525)*IT_0726 + (-0.125)*conj(IT_0254)*IT_0736 + (
      -0.125)*conj(IT_0373)*IT_0746 + (-0.125)*conj(IT_0021)*IT_0764 + (-0.125)
      *conj(IT_0293)*IT_0772 + (-0.125)*conj(IT_0456)*IT_0780 + (-0.125)*conj
      (IT_0155)*IT_0789 + (-0.125)*conj(IT_0169)*IT_0792;
    const ccomplex_t IT_0794 = IT_0466*conj(IT_0479) + IT_0476*conj(IT_0479) +
       conj(IT_0019)*IT_0578 + conj(IT_0277)*IT_0598 + conj(IT_0403)*IT_0657 +
       conj(IT_0436)*IT_0663 + conj(IT_0155)*IT_0726 + conj(IT_0021)*IT_0736 +
       conj(IT_0114)*IT_0746 + conj(IT_0254)*IT_0764 + conj(IT_0032)*IT_0772 +
       conj(IT_0169)*IT_0780 + conj(IT_0525)*IT_0789 + conj(IT_0456)*IT_0792;
    const ccomplex_t IT_0795 = conj(IT_0525)*IT_0529;
    const ccomplex_t IT_0796 = 0.25*conj(IT_0155);
    const ccomplex_t IT_0797 = 0.25*conj(IT_0525);
    const ccomplex_t IT_0798 = 0.25*conj(IT_0277);
    const ccomplex_t IT_0799 = 0.25*conj(IT_0032);
    const ccomplex_t IT_0800 = 0.25*conj(IT_0114);
    const ccomplex_t IT_0801 = 0.25*conj(IT_0479);
    const ccomplex_t IT_0802 = 0.25*conj(IT_0019);
    const ccomplex_t IT_0803 = 0.25*conj(IT_0293);
    const ccomplex_t IT_0804 = 0.25*conj(IT_0373);
    const ccomplex_t IT_0805 = 0.25*conj(IT_0254);
    const ccomplex_t IT_0806 = 0.25*conj(IT_0403);
    const ccomplex_t IT_0807 = 0.25*conj(IT_0436);
    const ccomplex_t IT_0808 = 0.25*conj(IT_0456);
    const ccomplex_t IT_0809 = conj(IT_0169)*IT_0207 + (-0.5)*conj(IT_0168)*
      (IT_0532 + 0.25*IT_0537 + 0.25*IT_0541) + (-0.125)*IT_0795 + (-0.5)*
      (IT_0166 + IT_0502)*IT_0796 + (-0.5)*IT_0527*IT_0797 + (-0.5)*IT_0291
      *IT_0798 + (-0.5)*IT_0230*IT_0799 + (-0.5)*IT_0252*IT_0800 + (-0.5)
      *IT_0486*IT_0801 + (-0.5)*IT_0491*IT_0801 + (-0.5)*IT_0500*IT_0802 + (-0.5
      )*IT_0304*IT_0803 + (-0.5)*IT_0386*IT_0804 + (-0.5)*IT_0275*IT_0805 + (
      -0.5)*IT_0423*IT_0806 + (-0.5)*IT_0455*IT_0807 + (-0.5)*IT_0477*IT_0808;
    const ccomplex_t IT_0810 = s_12*IT_0625;
    const ccomplex_t IT_0811 = m_C_1*IT_0301;
    const ccomplex_t IT_0812 = 2*IT_0811;
    const ccomplex_t IT_0813 = m_C_2*IT_0383;
    const ccomplex_t IT_0814 = (-2)*IT_0813;
    const ccomplex_t IT_0815 = s_23*m_C_1;
    const ccomplex_t IT_0816 = IT_0625*IT_0815;
    const ccomplex_t IT_0817 = 4*IT_0816;
    const ccomplex_t IT_0818 = s_13*m_C_2;
    const ccomplex_t IT_0819 = IT_0625*IT_0818;
    const ccomplex_t IT_0820 = (-4)*IT_0819;
    const ccomplex_t IT_0821 = 0.25*conj(IT_0021);
    const ccomplex_t IT_0822 = 0.25*conj(IT_0169);
    const ccomplex_t IT_0823 = 0.25*IT_0166*(conj(IT_0293) + conj(IT_0479)) +
       0.25*conj(IT_0403)*IT_0627 + IT_0627*IT_0796 + IT_0601*IT_0798 + IT_0529
      *IT_0800 + IT_0502*IT_0801 + IT_0560*IT_0802 + IT_0674*IT_0807 + IT_0726
      *IT_0808 + conj(IT_0525)*IT_0810 + IT_0799*IT_0812 + IT_0804*IT_0814 +
       IT_0805*IT_0817 + IT_0820*IT_0821 + IT_0789*IT_0822;
    const ccomplex_t IT_0824 = 0.25*(conj(IT_0373) + conj(IT_0479))*IT_0529 +
       0.25*conj(IT_0436)*IT_0627 + IT_0627*IT_0797 + IT_0560*IT_0798 + IT_0527
      *IT_0801 + IT_0601*IT_0802 + IT_0674*IT_0806 + IT_0789*IT_0808 + conj
      (IT_0155)*IT_0810 + IT_0803*IT_0812 + IT_0800*IT_0814 + IT_0805*IT_0820 +
       IT_0817*IT_0821 + IT_0726*IT_0822;
    const ccomplex_t IT_0825 = IT_0171 + IT_0183 + IT_0270 + IT_0271 + IT_0551
       + IT_0553;
    const ccomplex_t IT_0826 = 0.5*conj(IT_0155);
    const ccomplex_t IT_0827 = 0.5*conj(IT_0277);
    const ccomplex_t IT_0828 = (s_12 + -IT_0298)*(s_12 + (-0.5)*IT_0677);
    const ccomplex_t IT_0829 = (-2)*IT_0828;
    const ccomplex_t IT_0830 = 2*IT_0829;
    const ccomplex_t IT_0831 = 0.5*conj(IT_0032);
    const ccomplex_t IT_0832 = -(IT_0015 + -IT_0374)*(s_12 + (-0.5)*IT_0677);
    const ccomplex_t IT_0833 = (-2)*IT_0832;
    const ccomplex_t IT_0834 = (-2)*IT_0833;
    const ccomplex_t IT_0835 = 0.5*conj(IT_0114);
    const ccomplex_t IT_0836 = 0.5*conj(IT_0479);
    const ccomplex_t IT_0837 = -IT_0298;
    const ccomplex_t IT_0838 = s_12 + IT_0837;
    const ccomplex_t IT_0839 = IT_0294*IT_0838;
    const ccomplex_t IT_0840 = (-6)*IT_0839;
    const ccomplex_t IT_0841 = 0.5*conj(IT_0293);
    const ccomplex_t IT_0842 = IT_0294*IT_0640;
    const ccomplex_t IT_0843 = 6*IT_0842;
    const ccomplex_t IT_0844 = 0.5*conj(IT_0373);
    const ccomplex_t IT_0845 = 0.5*conj(IT_0403);
    const ccomplex_t IT_0846 = 0.5*conj(IT_0436);
    const ccomplex_t IT_0847 = 0.5*conj(IT_0169);
    const ccomplex_t IT_0848 = 0.5*IT_0795 + conj(IT_0254)*IT_0825 + IT_0814
      *IT_0826 + IT_0557*IT_0827 + IT_0830*IT_0831 + IT_0834*IT_0835 + IT_0381
      *IT_0836 + IT_0385*IT_0836 + IT_0840*IT_0841 + IT_0843*IT_0844 + IT_0692
      *IT_0845 + IT_0644*IT_0846 + IT_0746*IT_0847;
    const ccomplex_t IT_0849 = IT_0191 + IT_0194 + IT_0260 + IT_0261 + IT_0543
       + IT_0546;
    const ccomplex_t IT_0850 = 0.5*conj(IT_0525);
    const ccomplex_t IT_0851 = -(IT_0014 + -IT_0157)*(s_12 + (-0.5)*IT_0677);
    const ccomplex_t IT_0852 = (-2)*IT_0851;
    const ccomplex_t IT_0853 = (-2)*IT_0852;
    const ccomplex_t IT_0854 = IT_0157 + IT_0417;
    const ccomplex_t IT_0855 = IT_0294*IT_0854;
    const ccomplex_t IT_0856 = 6*IT_0855;
    const ccomplex_t IT_0857 = IT_0166*IT_0826 + IT_0550*IT_0827 + IT_0830
      *IT_0835 + IT_0296*IT_0836 + IT_0303*IT_0836 + IT_0840*IT_0844 + IT_0636
      *IT_0845 + IT_0700*IT_0846 + IT_0772*IT_0847 + conj(IT_0254)*IT_0849 +
       IT_0812*IT_0850 + IT_0831*IT_0853 + IT_0841*IT_0856;
    const ccomplex_t IT_0858 = 0.5*conj(IT_0019);
    const ccomplex_t IT_0859 = conj(IT_0254) + conj(IT_0277);
    const ccomplex_t IT_0860 = 0.5*IT_0859;
    const ccomplex_t IT_0861 = 0.5*conj(IT_0456);
    const ccomplex_t IT_0862 = 0.5*IT_0534 + conj(IT_0021)*IT_0825 + IT_0240
      *IT_0836 + IT_0251*IT_0836 + IT_0831*IT_0840 + IT_0830*IT_0841 + IT_0835
      *IT_0843 + IT_0834*IT_0844 + IT_0644*IT_0845 + IT_0692*IT_0846 + IT_0814
      *IT_0850 + IT_0557*IT_0858 + IT_0119*IT_0860 + IT_0746*IT_0861;
    const ccomplex_t IT_0863 = IT_0812*IT_0826 + IT_0218*IT_0836 + IT_0229
      *IT_0836 + IT_0835*IT_0840 + IT_0830*IT_0844 + IT_0700*IT_0845 + IT_0636
      *IT_0846 + conj(IT_0021)*IT_0849 + IT_0166*IT_0850 + IT_0841*IT_0853 +
       IT_0831*IT_0856 + IT_0550*IT_0858 + IT_0039*IT_0860 + IT_0772*IT_0861;
    const ccomplex_t IT_0864 = (IT_0279 + (-2)*IT_0280)*(1 + 0.5*IT_0623);
    const ccomplex_t IT_0865 = m_C_1*m_C_2*IT_0041;
    const ccomplex_t IT_0866 = IT_0625*IT_0865;
    const ccomplex_t IT_0867 = (-4)*IT_0866;
    const ccomplex_t IT_0868 = (-0.25)*conj(IT_0525);
    const ccomplex_t IT_0869 = (-0.25)*conj(IT_0403);
    const ccomplex_t IT_0870 = (-0.25)*conj(IT_0436);
    const ccomplex_t IT_0871 = (-0.25)*conj(IT_0169);
    const ccomplex_t IT_0872 = (-0.25)*conj(IT_0456);
    const ccomplex_t IT_0873 = (-0.25)*IT_0039*conj(IT_0293) + (-0.25)*IT_0119
      *conj(IT_0373) + (-0.25)*conj(IT_0277)*IT_0610 + (-0.25)*conj(IT_0155)
      *IT_0820 + (-0.5)*conj(IT_0114)*IT_0825 + (-0.5)*conj(IT_0032)*IT_0849 +
       conj(IT_0254)*IT_0864 + (-0.25)*conj(IT_0021)*IT_0867 + IT_0817*IT_0868 +
       IT_0708*IT_0869 + IT_0719*IT_0870 + IT_0764*IT_0871 + IT_0736*IT_0872;
    const ccomplex_t IT_0874 = (-0.25)*IT_0264*conj(IT_0479) + (-0.25)*IT_0274
      *conj(IT_0479) + (-0.25)*IT_0565 + (-0.25)*IT_0566 + (-0.25)*conj(IT_0019)
      *IT_0610 + (-0.25)*conj(IT_0155)*IT_0817 + (-0.5)*conj(IT_0373)*IT_0825 + 
      (-0.5)*conj(IT_0293)*IT_0849 + conj(IT_0021)*IT_0864 + (-0.25)*conj
      (IT_0254)*IT_0867 + IT_0820*IT_0868 + IT_0719*IT_0869 + IT_0708*IT_0870 +
       IT_0736*IT_0871 + IT_0764*IT_0872;
    const ccomplex_t IT_0875 = (conj(IT_0019) + conj(IT_0021))*IT_0121 + conj
      (IT_0032)*IT_0167 + (-8)*conj(IT_0168)*IT_0530 + (-8)*IT_0479*IT_0538 +
       IT_0479*IT_0542 + IT_0019*IT_0621 + IT_0277*IT_0622 + IT_0436*IT_0720 +
       IT_0403*IT_0721 + (-8)*IT_0456*IT_0793 + IT_0169*IT_0794 + (-8)*IT_0168
      *IT_0809 + 4*IT_0525*IT_0823 + 4*IT_0155*IT_0824 + 2*IT_0114*IT_0848 + 2
      *IT_0032*IT_0857 + 2*IT_0373*IT_0862 + 2*IT_0293*IT_0863 + (-4)*IT_0254
      *IT_0873 + (-4)*IT_0021*IT_0874;
    return create_ccomplex_return(IT_0875);
}

