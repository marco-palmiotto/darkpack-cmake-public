#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_W_W(
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
    const creal_t m_N_4 = param->m_N_4;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_N_4, 2);
    const ccomplex_t IT_0001 = pow(m_W, -4);
    const ccomplex_t IT_0002 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0003 = IT_0000 + IT_0002;
    const ccomplex_t IT_0004 = pow(m_W, -2);
    const ccomplex_t IT_0005 = pow(s_23, 2);
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = pow(s_24, 2);
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = -IT_0006 + -IT_0008;
    const ccomplex_t IT_0010 = IT_0003 + IT_0009;
    const ccomplex_t IT_0011 = m_N_4*IT_0010;
    const ccomplex_t IT_0012 = -IT_0011;
    const ccomplex_t IT_0013 = cos(theta_W);
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1)*e_em*IT_0013*IT_0015;
    const ccomplex_t IT_0017 = 2*IT_0016;
    const ccomplex_t IT_0018 = cpow(IT_0013, -1);
    const ccomplex_t IT_0019 = IT_0014*IT_0018;
    const ccomplex_t IT_0020 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0021 = IT_0019*IT_0020;
    const ccomplex_t IT_0022 = IT_0013*IT_0015;
    const ccomplex_t IT_0023 = IT_0020*IT_0022;
    const ccomplex_t IT_0024 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0025 = IT_0019*IT_0024;
    const ccomplex_t IT_0026 = IT_0022*IT_0024;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0021 + IT_0023 + 
      -IT_0025 + -IT_0026);
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = IT_0017*IT_0028;
    const ccomplex_t IT_0030 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0031 = (0 + _Complex_I*1)*IT_0030;
    const ccomplex_t IT_0032 = IT_0029*IT_0031;
    const ccomplex_t IT_0033 = (-0.5)*IT_0032;
    const ccomplex_t IT_0034 = -IT_0016;
    const ccomplex_t IT_0035 = (-0.5)*IT_0027;
    const ccomplex_t IT_0036 = IT_0034*IT_0035;
    const ccomplex_t IT_0037 = IT_0031*IT_0036;
    const ccomplex_t IT_0038 = 0.5*IT_0037;
    const ccomplex_t IT_0039 = IT_0017*IT_0035;
    const ccomplex_t IT_0040 = IT_0031*IT_0039;
    const ccomplex_t IT_0041 = conj(N_W4)*e_em*V_Wp1;
    const ccomplex_t IT_0042 = IT_0015*IT_0041;
    const ccomplex_t IT_0043 = conj(N_u4)*V_u1*e_em;
    const ccomplex_t IT_0044 = IT_0015*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0042 + (-0.5)*IT_0045);
    const ccomplex_t IT_0047 = -IT_0046;
    const ccomplex_t IT_0048 = N_W4*e_em*conj(V_Wp1);
    const ccomplex_t IT_0049 = IT_0015*IT_0048;
    const ccomplex_t IT_0050 = N_u4*conj(V_u1)*e_em;
    const ccomplex_t IT_0051 = IT_0015*IT_0050;
    const ccomplex_t IT_0052 = 1.4142135623731*IT_0051;
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*(IT_0049 + (-0.5)*IT_0052);
    const ccomplex_t IT_0054 = -IT_0053;
    const ccomplex_t IT_0055 = IT_0047*IT_0054;
    const ccomplex_t IT_0056 = pow(m_W, 2);
    const ccomplex_t IT_0057 = cpow((-2)*s_23 + IT_0000 + IT_0056 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0058 = IT_0055*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*IT_0058;
    const ccomplex_t IT_0060 = conj(N_W4)*e_em*V_Wp2;
    const ccomplex_t IT_0061 = IT_0015*IT_0060;
    const ccomplex_t IT_0062 = conj(N_u4)*V_u2*e_em;
    const ccomplex_t IT_0063 = IT_0015*IT_0062;
    const ccomplex_t IT_0064 = 1.4142135623731*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + (-0.5)*IT_0064);
    const ccomplex_t IT_0066 = -IT_0065;
    const ccomplex_t IT_0067 = N_W4*e_em*conj(V_Wp2);
    const ccomplex_t IT_0068 = IT_0015*IT_0067;
    const ccomplex_t IT_0069 = N_u4*conj(V_u2)*e_em;
    const ccomplex_t IT_0070 = IT_0015*IT_0069;
    const ccomplex_t IT_0071 = 1.4142135623731*IT_0070;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0068 + (-0.5)*IT_0071);
    const ccomplex_t IT_0073 = -IT_0072;
    const ccomplex_t IT_0074 = IT_0066*IT_0073;
    const ccomplex_t IT_0075 = cpow((-2)*s_23 + IT_0000 + IT_0056 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0074*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = 0.5*IT_0040 + 2*IT_0059 + 2*IT_0077;
    const ccomplex_t IT_0079 = IT_0028*IT_0034;
    const ccomplex_t IT_0080 = IT_0031*IT_0079;
    const ccomplex_t IT_0081 = (-0.5)*IT_0080;
    const ccomplex_t IT_0082 = cos(alpha);
    const ccomplex_t IT_0083 = IT_0018*IT_0082;
    const ccomplex_t IT_0084 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0085 = IT_0083*IT_0084;
    const ccomplex_t IT_0086 = sin(alpha);
    const ccomplex_t IT_0087 = IT_0018*IT_0086;
    const ccomplex_t IT_0088 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0089 = IT_0087*IT_0088;
    const ccomplex_t IT_0090 = IT_0015*IT_0082;
    const ccomplex_t IT_0091 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0092 = IT_0090*IT_0091;
    const ccomplex_t IT_0093 = IT_0015*IT_0086;
    const ccomplex_t IT_0094 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*(IT_0085 + IT_0089 + 
      -IT_0092 + -IT_0095);
    const ccomplex_t IT_0097 = -IT_0096;
    const ccomplex_t IT_0098 = sin(beta);
    const ccomplex_t IT_0099 = cos(beta);
    const ccomplex_t IT_0100 = (0 + _Complex_I*1)*m_W*e_em*IT_0015*(IT_0082
      *IT_0098 + -IT_0086*IT_0099);
    const ccomplex_t IT_0101 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0102 = IT_0100*IT_0101;
    const ccomplex_t IT_0103 = (0 + _Complex_I*1)*IT_0102;
    const ccomplex_t IT_0104 = IT_0097*IT_0103;
    const ccomplex_t IT_0105 = IT_0083*IT_0088;
    const ccomplex_t IT_0106 = IT_0084*IT_0087;
    const ccomplex_t IT_0107 = IT_0090*IT_0094;
    const ccomplex_t IT_0108 = IT_0091*IT_0093;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*(IT_0105 + -IT_0106 + 
      -IT_0107 + IT_0108);
    const ccomplex_t IT_0110 = (0 + _Complex_I*1)*m_W*e_em*IT_0015*(IT_0086
      *IT_0098 + IT_0082*IT_0099);
    const ccomplex_t IT_0111 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0112 = IT_0110*IT_0111;
    const ccomplex_t IT_0113 = (0 + _Complex_I*1)*IT_0112;
    const ccomplex_t IT_0114 = IT_0109*IT_0113;
    const ccomplex_t IT_0115 = m_N_4*IT_0036;
    const ccomplex_t IT_0116 = m_N_4*IT_0079;
    const ccomplex_t IT_0117 = -IT_0116;
    const ccomplex_t IT_0118 = IT_0115 + IT_0117;
    const ccomplex_t IT_0119 = IT_0031*IT_0118;
    const ccomplex_t IT_0120 = (-0.5)*IT_0104 + (-0.5)*IT_0114 + (-0.5)*IT_0119;
    const ccomplex_t IT_0121 = 2*IT_0120;
    const ccomplex_t IT_0122 = 0.5*IT_0104 + 0.5*IT_0114 + 0.5*IT_0119;
    const ccomplex_t IT_0123 = (-2)*IT_0122;
    const ccomplex_t IT_0124 = IT_0121 + IT_0123;
    const ccomplex_t IT_0125 = 2*conj(IT_0120);
    const ccomplex_t IT_0126 = (-2)*conj(IT_0122);
    const ccomplex_t IT_0127 = IT_0125 + IT_0126;
    const ccomplex_t IT_0128 = N_W4*e_em*conj(U_Wm1);
    const ccomplex_t IT_0129 = IT_0015*IT_0128;
    const ccomplex_t IT_0130 = N_d4*conj(U_d1)*e_em;
    const ccomplex_t IT_0131 = IT_0015*IT_0130;
    const ccomplex_t IT_0132 = 1.4142135623731*IT_0131;
    const ccomplex_t IT_0133 = (0 + _Complex_I*1)*(IT_0129 + 0.5*IT_0132);
    const ccomplex_t IT_0134 = -IT_0133;
    const ccomplex_t IT_0135 = conj(N_W4)*e_em*U_Wm1;
    const ccomplex_t IT_0136 = IT_0015*IT_0135;
    const ccomplex_t IT_0137 = conj(N_d4)*U_d1*e_em;
    const ccomplex_t IT_0138 = IT_0015*IT_0137;
    const ccomplex_t IT_0139 = 1.4142135623731*IT_0138;
    const ccomplex_t IT_0140 = (0 + _Complex_I*1)*(IT_0136 + 0.5*IT_0139);
    const ccomplex_t IT_0141 = -IT_0140;
    const ccomplex_t IT_0142 = IT_0134*IT_0141;
    const ccomplex_t IT_0143 = IT_0057*IT_0142;
    const ccomplex_t IT_0144 = (0 + _Complex_I*1)*IT_0143;
    const ccomplex_t IT_0145 = N_W4*e_em*conj(U_Wm2);
    const ccomplex_t IT_0146 = IT_0015*IT_0145;
    const ccomplex_t IT_0147 = N_d4*conj(U_d2)*e_em;
    const ccomplex_t IT_0148 = IT_0015*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1)*(IT_0146 + 0.5*IT_0149);
    const ccomplex_t IT_0151 = -IT_0150;
    const ccomplex_t IT_0152 = conj(N_W4)*e_em*U_Wm2;
    const ccomplex_t IT_0153 = IT_0015*IT_0152;
    const ccomplex_t IT_0154 = conj(N_d4)*U_d2*e_em;
    const ccomplex_t IT_0155 = IT_0015*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*(IT_0153 + 0.5*IT_0156);
    const ccomplex_t IT_0158 = -IT_0157;
    const ccomplex_t IT_0159 = IT_0151*IT_0158;
    const ccomplex_t IT_0160 = IT_0075*IT_0159;
    const ccomplex_t IT_0161 = (0 + _Complex_I*1)*IT_0160;
    const ccomplex_t IT_0162 = (-0.5)*IT_0040 + 2*IT_0144 + 2*IT_0161;
    const ccomplex_t IT_0163 = 0.5*IT_0080;
    const ccomplex_t IT_0164 = 0.5*IT_0032;
    const ccomplex_t IT_0165 = (-0.5)*IT_0037;
    const ccomplex_t IT_0166 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0167 = IT_0083*IT_0166;
    const ccomplex_t IT_0168 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0169 = IT_0087*IT_0168;
    const ccomplex_t IT_0170 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0171 = IT_0090*IT_0170;
    const ccomplex_t IT_0172 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0173 = IT_0093*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*(IT_0167 + -IT_0169 + 
      -IT_0171 + IT_0173);
    const ccomplex_t IT_0175 = IT_0113*IT_0174;
    const ccomplex_t IT_0176 = IT_0083*IT_0168;
    const ccomplex_t IT_0177 = IT_0087*IT_0166;
    const ccomplex_t IT_0178 = IT_0090*IT_0172;
    const ccomplex_t IT_0179 = IT_0093*IT_0170;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*(IT_0176 + IT_0177 + 
      -IT_0178 + -IT_0179);
    const ccomplex_t IT_0181 = -IT_0180;
    const ccomplex_t IT_0182 = IT_0103*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0119 + (-0.5)*IT_0175 + (-0.5)*IT_0182;
    const ccomplex_t IT_0184 = 2*IT_0183;
    const ccomplex_t IT_0185 = (-0.5)*IT_0119 + 0.5*IT_0175 + 0.5*IT_0182;
    const ccomplex_t IT_0186 = (-2)*IT_0185;
    const ccomplex_t IT_0187 = IT_0184 + IT_0186;
    const ccomplex_t IT_0188 = 2*conj(IT_0183);
    const ccomplex_t IT_0189 = (-2)*conj(IT_0185);
    const ccomplex_t IT_0190 = IT_0188 + IT_0189;
    const ccomplex_t IT_0191 = IT_0012*((conj(IT_0033) + conj(IT_0038) + conj
      (IT_0078) + conj(IT_0081))*IT_0124 + (IT_0033 + IT_0038 + IT_0078 +
       IT_0081)*IT_0127 + (conj(IT_0162) + conj(IT_0163) + conj(IT_0164) + conj
      (IT_0165))*IT_0187 + (IT_0162 + IT_0163 + IT_0164 + IT_0165)*IT_0190);
    const ccomplex_t IT_0192 = s_34*m_N_4*IT_0004;
    const ccomplex_t IT_0193 = s_23*s_34*IT_0004;
    const ccomplex_t IT_0194 = -IT_0193;
    const ccomplex_t IT_0195 = s_24 + IT_0194;
    const ccomplex_t IT_0196 = IT_0192*IT_0195;
    const ccomplex_t IT_0197 = IT_0016*IT_0035;
    const ccomplex_t IT_0198 = IT_0031*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = IT_0016*IT_0028;
    const ccomplex_t IT_0201 = IT_0031*IT_0200;
    const ccomplex_t IT_0202 = 0.5*IT_0201;
    const ccomplex_t IT_0203 = conj(IT_0199) + conj(IT_0202);
    const ccomplex_t IT_0204 = IT_0199 + IT_0202;
    const ccomplex_t IT_0205 = (-0.5)*IT_0201;
    const ccomplex_t IT_0206 = 0.5*IT_0198;
    const ccomplex_t IT_0207 = conj(IT_0205) + conj(IT_0206);
    const ccomplex_t IT_0208 = IT_0205 + IT_0206;
    const ccomplex_t IT_0209 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0210 = s_13*s_23*IT_0004;
    const ccomplex_t IT_0211 = s_14*s_24*IT_0004;
    const ccomplex_t IT_0212 = -IT_0210 + -IT_0211;
    const ccomplex_t IT_0213 = s_12 + IT_0209 + IT_0212;
    const ccomplex_t IT_0214 = m_N_4*IT_0213;
    const ccomplex_t IT_0215 = conj(IT_0163) + conj(IT_0165);
    const ccomplex_t IT_0216 = IT_0163 + IT_0165;
    const ccomplex_t IT_0217 = cpow((-2)*s_13 + IT_0000 + IT_0056 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0218 = IT_0159*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = (-0.5)*IT_0032 + (-2)*IT_0219;
    const ccomplex_t IT_0221 = cpow((-2)*s_13 + IT_0000 + IT_0056 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0222 = IT_0142*IT_0221;
    const ccomplex_t IT_0223 = (0 + _Complex_I*1)*IT_0222;
    const ccomplex_t IT_0224 = 0.5*IT_0040 + (-2)*IT_0223;
    const ccomplex_t IT_0225 = conj(IT_0220) + conj(IT_0224);
    const ccomplex_t IT_0226 = 2*IT_0122;
    const ccomplex_t IT_0227 = (-2)*IT_0120;
    const ccomplex_t IT_0228 = IT_0226 + IT_0227;
    const ccomplex_t IT_0229 = IT_0220 + IT_0224;
    const ccomplex_t IT_0230 = 2*conj(IT_0122);
    const ccomplex_t IT_0231 = (-2)*conj(IT_0120);
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = conj(IT_0038) + conj(IT_0081);
    const ccomplex_t IT_0234 = IT_0038 + IT_0081;
    const ccomplex_t IT_0235 = IT_0055*IT_0221;
    const ccomplex_t IT_0236 = (0 + _Complex_I*1)*IT_0235;
    const ccomplex_t IT_0237 = IT_0074*IT_0217;
    const ccomplex_t IT_0238 = (0 + _Complex_I*1)*IT_0237;
    const ccomplex_t IT_0239 = (-0.5)*IT_0040 + (-2)*IT_0236 + (-2)*IT_0238;
    const ccomplex_t IT_0240 = 2*IT_0185;
    const ccomplex_t IT_0241 = (-2)*IT_0183;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = 2*conj(IT_0185);
    const ccomplex_t IT_0244 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0245 = IT_0243 + IT_0244;
    const ccomplex_t IT_0246 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0247 = IT_0000 + IT_0246;
    const ccomplex_t IT_0248 = pow(s_13, 2);
    const ccomplex_t IT_0249 = IT_0004*IT_0248;
    const ccomplex_t IT_0250 = pow(s_14, 2);
    const ccomplex_t IT_0251 = IT_0004*IT_0250;
    const ccomplex_t IT_0252 = -IT_0249 + -IT_0251;
    const ccomplex_t IT_0253 = IT_0247 + IT_0252;
    const ccomplex_t IT_0254 = m_N_4*IT_0253;
    const ccomplex_t IT_0255 = -IT_0254;
    const ccomplex_t IT_0256 = (IT_0232*(IT_0163 + IT_0164 + IT_0165 + IT_0239
      ) + IT_0228*(conj(IT_0163) + conj(IT_0164) + conj(IT_0165) + conj(IT_0239)
      ) + (conj(IT_0038) + conj(IT_0081) + conj(IT_0220) + conj(IT_0224))
      *IT_0242 + (IT_0038 + IT_0081 + IT_0220 + IT_0224)*IT_0245)*IT_0255;
    const ccomplex_t IT_0257 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0258 = s_12 + IT_0212 + IT_0257;
    const ccomplex_t IT_0259 = m_N_4*IT_0258;
    const ccomplex_t IT_0260 = conj(IT_0033) + conj(IT_0078);
    const ccomplex_t IT_0261 = IT_0033 + IT_0078;
    const ccomplex_t IT_0262 = s_13*s_34*IT_0004;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = s_14 + IT_0263;
    const ccomplex_t IT_0265 = IT_0192*IT_0264;
    const ccomplex_t IT_0266 = IT_0000*IT_0253;
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = IT_0000*IT_0010;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = 2*IT_0038;
    const ccomplex_t IT_0271 = conj(IT_0164)*IT_0270;
    const ccomplex_t IT_0272 = 2*conj(IT_0038);
    const ccomplex_t IT_0273 = IT_0164*IT_0272;
    const ccomplex_t IT_0274 = 2*IT_0164;
    const ccomplex_t IT_0275 = conj(IT_0081)*IT_0274;
    const ccomplex_t IT_0276 = 2*conj(IT_0164);
    const ccomplex_t IT_0277 = IT_0081*IT_0276;
    const ccomplex_t IT_0278 = s_12*s_14*s_34*IT_0004;
    const ccomplex_t IT_0279 = s_12*s_23;
    const ccomplex_t IT_0280 = s_13*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0281 = s_13*IT_0000;
    const ccomplex_t IT_0282 = s_12*s_24*s_34*IT_0004;
    const ccomplex_t IT_0283 = s_24*s_34*IT_0001*IT_0248;
    const ccomplex_t IT_0284 = s_14*s_34*IT_0001*IT_0005;
    const ccomplex_t IT_0285 = s_12*s_13;
    const ccomplex_t IT_0286 = s_13*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0287 = s_23*IT_0000;
    const ccomplex_t IT_0288 = IT_0285 + IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = s_13*s_14*s_24*IT_0004;
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = s_23*IT_0004*IT_0250;
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = s_23*IT_0004*IT_0248;
    const ccomplex_t IT_0294 = (-2)*IT_0293;
    const ccomplex_t IT_0295 = pow(s_34, 2);
    const ccomplex_t IT_0296 = s_12*s_13*IT_0001*IT_0295;
    const ccomplex_t IT_0297 = -IT_0296;
    const ccomplex_t IT_0298 = s_13*IT_0004*IT_0007;
    const ccomplex_t IT_0299 = -IT_0298;
    const ccomplex_t IT_0300 = s_14*s_23*s_24*IT_0004;
    const ccomplex_t IT_0301 = -IT_0300;
    const ccomplex_t IT_0302 = s_13*IT_0004*IT_0005;
    const ccomplex_t IT_0303 = (-2)*IT_0302;
    const ccomplex_t IT_0304 = s_12*s_23*IT_0001*IT_0295;
    const ccomplex_t IT_0305 = -IT_0304;
    const ccomplex_t IT_0306 = (-2)*IT_0016;
    const ccomplex_t IT_0307 = IT_0028*IT_0306;
    const ccomplex_t IT_0308 = IT_0031*IT_0307;
    const ccomplex_t IT_0309 = (-0.5)*IT_0308;
    const ccomplex_t IT_0310 = IT_0274*conj(IT_0309);
    const ccomplex_t IT_0311 = IT_0276*IT_0309;
    const ccomplex_t IT_0312 = s_12*IT_0000;
    const ccomplex_t IT_0313 = s_14*s_24*IT_0001*IT_0005;
    const ccomplex_t IT_0314 = s_13*s_23*IT_0001*IT_0007;
    const ccomplex_t IT_0315 = s_14*s_24*IT_0000*IT_0004;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = s_13*s_23*IT_0000*IT_0004;
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = s_12*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = IT_0312 + IT_0313 + IT_0314 + IT_0316 + IT_0318
       + IT_0320;
    const ccomplex_t IT_0322 = IT_0001*IT_0005*IT_0250;
    const ccomplex_t IT_0323 = s_13*s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0324 = pow(m_N_4, 4);
    const ccomplex_t IT_0325 = IT_0000*IT_0004*IT_0005;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = IT_0000*IT_0004*IT_0250;
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = s_12*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0330 = -IT_0329;
    const ccomplex_t IT_0331 = IT_0322 + IT_0323 + IT_0324 + IT_0326 + IT_0328
       + IT_0330;
    const ccomplex_t IT_0332 = conj(IT_0033)*IT_0270;
    const ccomplex_t IT_0333 = conj(IT_0078)*IT_0270;
    const ccomplex_t IT_0334 = IT_0033*IT_0272;
    const ccomplex_t IT_0335 = IT_0078*IT_0272;
    const ccomplex_t IT_0336 = 2*IT_0078;
    const ccomplex_t IT_0337 = conj(IT_0081)*IT_0336;
    const ccomplex_t IT_0338 = 2*conj(IT_0078);
    const ccomplex_t IT_0339 = IT_0081*IT_0338;
    const ccomplex_t IT_0340 = 2*IT_0163;
    const ccomplex_t IT_0341 = conj(IT_0162)*IT_0340;
    const ccomplex_t IT_0342 = conj(IT_0164)*IT_0340;
    const ccomplex_t IT_0343 = 2*conj(IT_0163);
    const ccomplex_t IT_0344 = IT_0162*IT_0343;
    const ccomplex_t IT_0345 = IT_0164*IT_0343;
    const ccomplex_t IT_0346 = conj(IT_0165)*IT_0274;
    const ccomplex_t IT_0347 = IT_0165*IT_0276;
    const ccomplex_t IT_0348 = 2*IT_0165;
    const ccomplex_t IT_0349 = conj(IT_0162)*IT_0348;
    const ccomplex_t IT_0350 = 2*conj(IT_0165);
    const ccomplex_t IT_0351 = IT_0162*IT_0350;
    const ccomplex_t IT_0352 = 2*IT_0081;
    const ccomplex_t IT_0353 = conj(IT_0033)*IT_0352;
    const ccomplex_t IT_0354 = 2*conj(IT_0081);
    const ccomplex_t IT_0355 = IT_0033*IT_0354;
    const ccomplex_t IT_0356 = s_14*s_34*IT_0004;
    const ccomplex_t IT_0357 = s_13*IT_0056;
    const ccomplex_t IT_0358 = s_14*s_34;
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = IT_0357 + IT_0359;
    const ccomplex_t IT_0361 = IT_0004*IT_0360;
    const ccomplex_t IT_0362 = -IT_0356 + (-0.25)*IT_0361;
    const ccomplex_t IT_0363 = s_13 + IT_0362;
    const ccomplex_t IT_0364 = IT_0000*IT_0363;
    const ccomplex_t IT_0365 = (-4)*IT_0364;
    const ccomplex_t IT_0366 = s_24*s_34*IT_0004;
    const ccomplex_t IT_0367 = s_23*IT_0056;
    const ccomplex_t IT_0368 = s_24*s_34;
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = IT_0367 + IT_0369;
    const ccomplex_t IT_0371 = IT_0004*IT_0370;
    const ccomplex_t IT_0372 = -IT_0366 + (-0.25)*IT_0371;
    const ccomplex_t IT_0373 = s_23 + IT_0372;
    const ccomplex_t IT_0374 = IT_0000*IT_0373;
    const ccomplex_t IT_0375 = (-4)*IT_0374;
    const ccomplex_t IT_0376 = IT_0236 + IT_0238;
    const ccomplex_t IT_0377 = IT_0270*conj(IT_0376);
    const ccomplex_t IT_0378 = -IT_0144 + -IT_0161;
    const ccomplex_t IT_0379 = IT_0270*conj(IT_0378);
    const ccomplex_t IT_0380 = IT_0272*IT_0376;
    const ccomplex_t IT_0381 = IT_0272*IT_0378;
    const ccomplex_t IT_0382 = -IT_0059 + -IT_0077;
    const ccomplex_t IT_0383 = IT_0340*conj(IT_0382);
    const ccomplex_t IT_0384 = conj(IT_0223)*IT_0340;
    const ccomplex_t IT_0385 = conj(IT_0219)*IT_0340;
    const ccomplex_t IT_0386 = IT_0343*IT_0382;
    const ccomplex_t IT_0387 = IT_0223*IT_0343;
    const ccomplex_t IT_0388 = IT_0219*IT_0343;
    const ccomplex_t IT_0389 = IT_0348*conj(IT_0382);
    const ccomplex_t IT_0390 = conj(IT_0223)*IT_0348;
    const ccomplex_t IT_0391 = conj(IT_0219)*IT_0348;
    const ccomplex_t IT_0392 = IT_0350*IT_0382;
    const ccomplex_t IT_0393 = IT_0223*IT_0350;
    const ccomplex_t IT_0394 = IT_0219*IT_0350;
    const ccomplex_t IT_0395 = IT_0352*conj(IT_0376);
    const ccomplex_t IT_0396 = IT_0352*conj(IT_0378);
    const ccomplex_t IT_0397 = IT_0354*IT_0376;
    const ccomplex_t IT_0398 = IT_0354*IT_0378;
    const ccomplex_t IT_0399 = (IT_0365 + IT_0375)*(IT_0377 + IT_0379 +
       IT_0380 + IT_0381 + IT_0383 + IT_0384 + IT_0385 + IT_0386 + IT_0387 +
       IT_0388 + IT_0389 + IT_0390 + IT_0391 + IT_0392 + IT_0393 + IT_0394 +
       IT_0395 + IT_0396 + IT_0397 + IT_0398);
    const ccomplex_t IT_0400 = IT_0004*IT_0295;
    const ccomplex_t IT_0401 = s_12*IT_0056;
    const ccomplex_t IT_0402 = s_13*s_23;
    const ccomplex_t IT_0403 = (-2)*IT_0402;
    const ccomplex_t IT_0404 = IT_0401 + IT_0403;
    const ccomplex_t IT_0405 = IT_0004*IT_0404;
    const ccomplex_t IT_0406 = (IT_0056 + -IT_0400)*(s_12 + (-0.5)*IT_0405);
    const ccomplex_t IT_0407 = (-2)*IT_0406;
    const ccomplex_t IT_0408 = 2*IT_0199;
    const ccomplex_t IT_0409 = IT_0279 + IT_0280 + IT_0281;
    const ccomplex_t IT_0410 = s_14*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0411 = -IT_0410;
    const ccomplex_t IT_0412 = IT_0284 + IT_0303 + IT_0305 + IT_0409 + IT_0411;
    const ccomplex_t IT_0413 = conj(IT_0206)*IT_0274;
    const ccomplex_t IT_0414 = IT_0206*IT_0276;
    const ccomplex_t IT_0415 = IT_0413 + IT_0414;
    const ccomplex_t IT_0416 = 2*conj(IT_0199);
    const ccomplex_t IT_0417 = 2*IT_0205;
    const ccomplex_t IT_0418 = conj(IT_0164)*IT_0417;
    const ccomplex_t IT_0419 = 2*conj(IT_0205);
    const ccomplex_t IT_0420 = IT_0164*IT_0419;
    const ccomplex_t IT_0421 = s_24*s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0422 = -IT_0421;
    const ccomplex_t IT_0423 = IT_0283 + IT_0288 + IT_0294 + IT_0297 + IT_0422;
    const ccomplex_t IT_0424 = (s_13 + -IT_0356)*(s_12 + (-0.5)*IT_0405);
    const ccomplex_t IT_0425 = (-2)*IT_0424;
    const ccomplex_t IT_0426 = (s_23 + -IT_0366)*(s_12 + (-0.5)*IT_0405);
    const ccomplex_t IT_0427 = (-2)*IT_0426;
    const ccomplex_t IT_0428 = conj(IT_0199)*IT_0270;
    const ccomplex_t IT_0429 = conj(IT_0202)*IT_0270;
    const ccomplex_t IT_0430 = IT_0199*IT_0272;
    const ccomplex_t IT_0431 = IT_0202*IT_0272;
    const ccomplex_t IT_0432 = conj(IT_0205)*IT_0340;
    const ccomplex_t IT_0433 = conj(IT_0206)*IT_0340;
    const ccomplex_t IT_0434 = IT_0205*IT_0343;
    const ccomplex_t IT_0435 = IT_0206*IT_0343;
    const ccomplex_t IT_0436 = conj(IT_0206)*IT_0348;
    const ccomplex_t IT_0437 = IT_0206*IT_0350;
    const ccomplex_t IT_0438 = conj(IT_0202)*IT_0352;
    const ccomplex_t IT_0439 = IT_0202*IT_0354;
    const ccomplex_t IT_0440 = conj(IT_0081)*IT_0408;
    const ccomplex_t IT_0441 = IT_0081*IT_0416;
    const ccomplex_t IT_0442 = conj(IT_0165)*IT_0417;
    const ccomplex_t IT_0443 = IT_0165*IT_0419;
    const ccomplex_t IT_0444 = s_14*s_24;
    const ccomplex_t IT_0445 = (-2)*IT_0444;
    const ccomplex_t IT_0446 = IT_0401 + IT_0445;
    const ccomplex_t IT_0447 = IT_0004*IT_0446;
    const ccomplex_t IT_0448 = (s_12 + -IT_0210)*(s_12 + (-0.5)*IT_0447);
    const ccomplex_t IT_0449 = (-2)*IT_0448;
    const ccomplex_t IT_0450 = 2*conj(IT_0224);
    const ccomplex_t IT_0451 = s_13*s_23*IT_0001*IT_0250;
    const ccomplex_t IT_0452 = s_14*s_24*IT_0001*IT_0248;
    const ccomplex_t IT_0453 = s_12*s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0454 = -IT_0453;
    const ccomplex_t IT_0455 = IT_0312 + IT_0316 + IT_0318 + IT_0451 + IT_0452
       + IT_0454;
    const ccomplex_t IT_0456 = IT_0001*IT_0007*IT_0248;
    const ccomplex_t IT_0457 = IT_0000*IT_0004*IT_0007;
    const ccomplex_t IT_0458 = -IT_0457;
    const ccomplex_t IT_0459 = IT_0000*IT_0004*IT_0248;
    const ccomplex_t IT_0460 = -IT_0459;
    const ccomplex_t IT_0461 = s_12*s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0462 = -IT_0461;
    const ccomplex_t IT_0463 = IT_0323 + IT_0324 + IT_0456 + IT_0458 + IT_0460
       + IT_0462;
    const ccomplex_t IT_0464 = IT_0455 + IT_0463;
    const ccomplex_t IT_0465 = conj(IT_0220)*IT_0270;
    const ccomplex_t IT_0466 = conj(IT_0224)*IT_0270;
    const ccomplex_t IT_0467 = IT_0220*IT_0272;
    const ccomplex_t IT_0468 = IT_0224*IT_0272;
    const ccomplex_t IT_0469 = 2*IT_0224;
    const ccomplex_t IT_0470 = conj(IT_0081)*IT_0469;
    const ccomplex_t IT_0471 = IT_0081*IT_0450;
    const ccomplex_t IT_0472 = s_34*IT_0000*IT_0004;
    const ccomplex_t IT_0473 = IT_0195*IT_0472;
    const ccomplex_t IT_0474 = IT_0035*IT_0306;
    const ccomplex_t IT_0475 = IT_0031*IT_0474;
    const ccomplex_t IT_0476 = (-0.5)*IT_0475;
    const ccomplex_t IT_0477 = 2*IT_0476;
    const ccomplex_t IT_0478 = IT_0408 + IT_0477;
    const ccomplex_t IT_0479 = 2*conj(IT_0476);
    const ccomplex_t IT_0480 = IT_0416 + IT_0479;
    const ccomplex_t IT_0481 = 2*IT_0206;
    const ccomplex_t IT_0482 = IT_0417 + IT_0481;
    const ccomplex_t IT_0483 = 2*IT_0309;
    const ccomplex_t IT_0484 = IT_0482 + IT_0483;
    const ccomplex_t IT_0485 = 2*conj(IT_0206);
    const ccomplex_t IT_0486 = IT_0419 + IT_0485;
    const ccomplex_t IT_0487 = 2*conj(IT_0309);
    const ccomplex_t IT_0488 = IT_0486 + IT_0487;
    const ccomplex_t IT_0489 = 0.5*IT_0308;
    const ccomplex_t IT_0490 = conj(IT_0202) + conj(IT_0489);
    const ccomplex_t IT_0491 = 2*IT_0162;
    const ccomplex_t IT_0492 = IT_0202 + IT_0489;
    const ccomplex_t IT_0493 = 2*conj(IT_0162);
    const ccomplex_t IT_0494 = 0.5*IT_0475;
    const ccomplex_t IT_0495 = 2*IT_0033;
    const ccomplex_t IT_0496 = 2*conj(IT_0033);
    const ccomplex_t IT_0497 = conj(IT_0205) + conj(IT_0206) + conj(IT_0309) +
       conj(IT_0494);
    const ccomplex_t IT_0498 = conj(IT_0202)*IT_0274;
    const ccomplex_t IT_0499 = IT_0274*conj(IT_0489);
    const ccomplex_t IT_0500 = IT_0202*IT_0276;
    const ccomplex_t IT_0501 = IT_0276*IT_0489;
    const ccomplex_t IT_0502 = conj(IT_0164)*IT_0408;
    const ccomplex_t IT_0503 = IT_0164*IT_0416;
    const ccomplex_t IT_0504 = conj(IT_0164)*IT_0477;
    const ccomplex_t IT_0505 = IT_0164*IT_0479;
    const ccomplex_t IT_0506 = conj(IT_0162)*IT_0478 + IT_0162*IT_0480 + conj
      (IT_0033)*IT_0484 + IT_0033*IT_0488 + IT_0490*IT_0491 + IT_0492*IT_0493 +
       IT_0338*(IT_0205 + IT_0206 + IT_0309 + IT_0494) + conj(IT_0494)*IT_0495 +
       IT_0494*IT_0496 + IT_0336*IT_0497 + IT_0498 + IT_0499 + IT_0500 + IT_0501
       + IT_0502 + IT_0503 + IT_0504 + IT_0505;
    const ccomplex_t IT_0507 = IT_0264*IT_0472;
    const ccomplex_t IT_0508 = 2*IT_0239;
    const ccomplex_t IT_0509 = 2*conj(IT_0239);
    const ccomplex_t IT_0510 = 2*IT_0220;
    const ccomplex_t IT_0511 = 2*conj(IT_0220);
    const ccomplex_t IT_0512 = conj(IT_0239)*IT_0478 + IT_0239*IT_0480 + conj
      (IT_0220)*IT_0484 + IT_0220*IT_0488 + IT_0450*(IT_0205 + IT_0206 + IT_0309
       + IT_0494) + IT_0469*IT_0497 + IT_0498 + IT_0499 + IT_0500 + IT_0501 +
       IT_0502 + IT_0503 + IT_0504 + IT_0505 + IT_0490*IT_0508 + IT_0492*IT_0509
       + conj(IT_0494)*IT_0510 + IT_0494*IT_0511;
    const ccomplex_t IT_0513 = 2*IT_0286;
    const ccomplex_t IT_0514 = IT_0285 + IT_0287 + IT_0290 + IT_0292 + IT_0294
       + IT_0513;
    const ccomplex_t IT_0515 = 2*IT_0280;
    const ccomplex_t IT_0516 = IT_0279 + IT_0281 + IT_0299 + IT_0301 + IT_0303
       + IT_0515;
    const ccomplex_t IT_0517 = IT_0270*conj(IT_0476);
    const ccomplex_t IT_0518 = IT_0270*conj(IT_0489);
    const ccomplex_t IT_0519 = IT_0272*IT_0476;
    const ccomplex_t IT_0520 = IT_0272*IT_0489;
    const ccomplex_t IT_0521 = IT_0340*conj(IT_0494);
    const ccomplex_t IT_0522 = conj(IT_0309)*IT_0340;
    const ccomplex_t IT_0523 = IT_0343*IT_0494;
    const ccomplex_t IT_0524 = IT_0309*IT_0343;
    const ccomplex_t IT_0525 = IT_0348*conj(IT_0494);
    const ccomplex_t IT_0526 = IT_0350*IT_0494;
    const ccomplex_t IT_0527 = IT_0352*conj(IT_0489);
    const ccomplex_t IT_0528 = IT_0354*IT_0489;
    const ccomplex_t IT_0529 = conj(IT_0081)*IT_0477;
    const ccomplex_t IT_0530 = IT_0081*IT_0479;
    const ccomplex_t IT_0531 = conj(IT_0165)*IT_0483;
    const ccomplex_t IT_0532 = IT_0165*IT_0487;
    const ccomplex_t IT_0533 = (s_12 + -IT_0211)*(s_12 + (-0.5)*IT_0405);
    const ccomplex_t IT_0534 = (-2)*IT_0533;
    const ccomplex_t IT_0535 = 4*IT_0038;
    const ccomplex_t IT_0536 = -IT_0251;
    const ccomplex_t IT_0537 = IT_0000 + IT_0536;
    const ccomplex_t IT_0538 = IT_0000*IT_0537;
    const ccomplex_t IT_0539 = -IT_0538;
    const ccomplex_t IT_0540 = -IT_0008;
    const ccomplex_t IT_0541 = IT_0000 + IT_0540;
    const ccomplex_t IT_0542 = IT_0000*IT_0541;
    const ccomplex_t IT_0543 = -IT_0542;
    const ccomplex_t IT_0544 = IT_0539 + IT_0543;
    const ccomplex_t IT_0545 = 6*IT_0038;
    const ccomplex_t IT_0546 = conj(IT_0165)*IT_0545;
    const ccomplex_t IT_0547 = 6*conj(IT_0038);
    const ccomplex_t IT_0548 = IT_0165*IT_0547;
    const ccomplex_t IT_0549 = conj(IT_0081)*IT_0163;
    const ccomplex_t IT_0550 = 6*IT_0549;
    const ccomplex_t IT_0551 = IT_0081*conj(IT_0163);
    const ccomplex_t IT_0552 = 6*IT_0551;
    const ccomplex_t IT_0553 = conj(IT_0081)*IT_0165;
    const ccomplex_t IT_0554 = 6*IT_0553;
    const ccomplex_t IT_0555 = IT_0081*conj(IT_0165);
    const ccomplex_t IT_0556 = 6*IT_0555;
    const ccomplex_t IT_0557 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0558 = IT_0066*IT_0158;
    const ccomplex_t IT_0559 = IT_0217*IT_0557*IT_0558;
    const ccomplex_t IT_0560 = m_N_4*IT_0236;
    const ccomplex_t IT_0561 = m_N_4*IT_0238;
    const ccomplex_t IT_0562 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0563 = IT_0047*IT_0141;
    const ccomplex_t IT_0564 = IT_0221*IT_0562*IT_0563;
    const ccomplex_t IT_0565 = -IT_0559 + IT_0560 + IT_0561 + -IT_0564;
    const ccomplex_t IT_0566 = m_N_4*IT_0223;
    const ccomplex_t IT_0567 = -IT_0566;
    const ccomplex_t IT_0568 = 1.125*IT_0000;
    const ccomplex_t IT_0569 = s_34*IT_0001;
    const ccomplex_t IT_0570 = s_14*s_23;
    const ccomplex_t IT_0571 = s_13*s_24;
    const ccomplex_t IT_0572 = s_12*s_34;
    const ccomplex_t IT_0573 = -IT_0571 + -IT_0572;
    const ccomplex_t IT_0574 = IT_0570 + IT_0573;
    const ccomplex_t IT_0575 = IT_0569*IT_0574;
    const ccomplex_t IT_0576 = (-0.25)*IT_0575;
    const ccomplex_t IT_0577 = IT_0004*IT_0401;
    const ccomplex_t IT_0578 = (-0.5)*IT_0577;
    const ccomplex_t IT_0579 = s_12 + IT_0576 + IT_0578;
    const ccomplex_t IT_0580 = -conj(IT_0185);
    const ccomplex_t IT_0581 = conj(IT_0183) + IT_0580;
    const ccomplex_t IT_0582 = IT_0000*IT_0056;
    const ccomplex_t IT_0583 = s_13*s_14*s_34*IT_0004;
    const ccomplex_t IT_0584 = 0.5*IT_0250 + (-0.5)*IT_0582 + -IT_0583;
    const ccomplex_t IT_0585 = IT_0248 + IT_0584;
    const ccomplex_t IT_0586 = IT_0004*IT_0585;
    const ccomplex_t IT_0587 = IT_0000 + IT_0536 + IT_0586;
    const ccomplex_t IT_0588 = m_N_4*IT_0587;
    const ccomplex_t IT_0589 = 4*conj(IT_0038);
    const ccomplex_t IT_0590 = 0.25*IT_0589;
    const ccomplex_t IT_0591 = conj(IT_0081) + IT_0590;
    const ccomplex_t IT_0592 = s_13*s_24*s_34*IT_0004;
    const ccomplex_t IT_0593 = -IT_0402 + (-2)*IT_0444 + 2*IT_0592;
    const ccomplex_t IT_0594 = IT_0401 + IT_0593;
    const ccomplex_t IT_0595 = IT_0004*IT_0594;
    const ccomplex_t IT_0596 = -IT_0210 + (-0.5)*IT_0595;
    const ccomplex_t IT_0597 = s_12 + IT_0596;
    const ccomplex_t IT_0598 = m_N_4*IT_0597;
    const ccomplex_t IT_0599 = (-2)*IT_0598;
    const ccomplex_t IT_0600 = IT_0276 + IT_0509;
    const ccomplex_t IT_0601 = -IT_0402;
    const ccomplex_t IT_0602 = IT_0401 + IT_0601;
    const ccomplex_t IT_0603 = IT_0004*IT_0602;
    const ccomplex_t IT_0604 = -IT_0210 + (-0.25)*IT_0603;
    const ccomplex_t IT_0605 = s_12 + IT_0604;
    const ccomplex_t IT_0606 = m_N_4*IT_0605;
    const ccomplex_t IT_0607 = (-4)*IT_0606;
    const ccomplex_t IT_0608 = IT_0338 + IT_0496;
    const ccomplex_t IT_0609 = IT_0004*IT_0367;
    const ccomplex_t IT_0610 = s_24*IT_0056;
    const ccomplex_t IT_0611 = IT_0569*IT_0610;
    const ccomplex_t IT_0612 = (-0.5)*IT_0609 + 0.25*IT_0611;
    const ccomplex_t IT_0613 = s_23 + IT_0612;
    const ccomplex_t IT_0614 = m_N_4*IT_0613;
    const ccomplex_t IT_0615 = 4*IT_0614;
    const ccomplex_t IT_0616 = 2*conj(IT_0494);
    const ccomplex_t IT_0617 = IT_0487 + IT_0616;
    const ccomplex_t IT_0618 = -IT_0582;
    const ccomplex_t IT_0619 = IT_0007 + IT_0618;
    const ccomplex_t IT_0620 = IT_0004*IT_0619;
    const ccomplex_t IT_0621 = -IT_0000;
    const ccomplex_t IT_0622 = IT_0008 + IT_0621;
    const ccomplex_t IT_0623 = (-4)*IT_0622;
    const ccomplex_t IT_0624 = IT_0620 + IT_0623;
    const ccomplex_t IT_0625 = m_N_4*IT_0624;
    const ccomplex_t IT_0626 = IT_0343 + IT_0350;
    const ccomplex_t IT_0627 = IT_0357 + IT_0358;
    const ccomplex_t IT_0628 = s_13*IT_0004*IT_0295;
    const ccomplex_t IT_0629 = (-2)*IT_0628;
    const ccomplex_t IT_0630 = IT_0627 + IT_0629;
    const ccomplex_t IT_0631 = IT_0004*IT_0630;
    const ccomplex_t IT_0632 = -IT_0356 + 0.5*IT_0631;
    const ccomplex_t IT_0633 = s_13 + IT_0632;
    const ccomplex_t IT_0634 = m_N_4*IT_0633;
    const ccomplex_t IT_0635 = 2*IT_0634;
    const ccomplex_t IT_0636 = 2*conj(IT_0202);
    const ccomplex_t IT_0637 = IT_0416 + IT_0636;
    const ccomplex_t IT_0638 = IT_0004*IT_0357;
    const ccomplex_t IT_0639 = s_14*IT_0056;
    const ccomplex_t IT_0640 = s_13*s_34;
    const ccomplex_t IT_0641 = (-2)*IT_0640;
    const ccomplex_t IT_0642 = IT_0639 + IT_0641;
    const ccomplex_t IT_0643 = IT_0569*IT_0642;
    const ccomplex_t IT_0644 = (-0.5)*IT_0638 + (-0.25)*IT_0643;
    const ccomplex_t IT_0645 = s_13 + IT_0644;
    const ccomplex_t IT_0646 = m_N_4*IT_0645;
    const ccomplex_t IT_0647 = (-4)*IT_0646;
    const ccomplex_t IT_0648 = 2*conj(IT_0489);
    const ccomplex_t IT_0649 = IT_0479 + IT_0648;
    const ccomplex_t IT_0650 = pow(m_W, 4);
    const ccomplex_t IT_0651 = (-2)*IT_0295;
    const ccomplex_t IT_0652 = IT_0650 + IT_0651;
    const ccomplex_t IT_0653 = IT_0001*IT_0652;
    const ccomplex_t IT_0654 = 0.25*IT_0653;
    const ccomplex_t IT_0655 = 1 + IT_0654;
    const ccomplex_t IT_0656 = IT_0000*IT_0655;
    const ccomplex_t IT_0657 = m_N_4*IT_0161;
    const ccomplex_t IT_0658 = 8*conj(IT_0657);
    const ccomplex_t IT_0659 = IT_0054*IT_0134;
    const ccomplex_t IT_0660 = IT_0057*IT_0562*IT_0659;
    const ccomplex_t IT_0661 = IT_0073*IT_0151;
    const ccomplex_t IT_0662 = IT_0075*IT_0557*IT_0661;
    const ccomplex_t IT_0663 = IT_0660 + IT_0662;
    const ccomplex_t IT_0664 = (-8)*conj(IT_0663);
    const ccomplex_t IT_0665 = IT_0658 + IT_0664;
    const ccomplex_t IT_0666 = s_23*s_24*s_34*IT_0004;
    const ccomplex_t IT_0667 = 2*IT_0007 + -IT_0582 + (-2)*IT_0666;
    const ccomplex_t IT_0668 = IT_0005 + IT_0667;
    const ccomplex_t IT_0669 = IT_0004*IT_0668;
    const ccomplex_t IT_0670 = -IT_0006 + 0.5*IT_0669;
    const ccomplex_t IT_0671 = IT_0000 + IT_0670;
    const ccomplex_t IT_0672 = m_N_4*IT_0671;
    const ccomplex_t IT_0673 = 2*IT_0672;
    const ccomplex_t IT_0674 = (-2)*conj(IT_0162);
    const ccomplex_t IT_0675 = (-2)*conj(IT_0164);
    const ccomplex_t IT_0676 = IT_0674 + IT_0675;
    const ccomplex_t IT_0677 = IT_0248 + IT_0618;
    const ccomplex_t IT_0678 = IT_0004*IT_0677;
    const ccomplex_t IT_0679 = IT_0249 + IT_0621;
    const ccomplex_t IT_0680 = (-4)*IT_0679;
    const ccomplex_t IT_0681 = IT_0678 + IT_0680;
    const ccomplex_t IT_0682 = m_N_4*IT_0681;
    const ccomplex_t IT_0683 = (-2)*conj(IT_0220);
    const ccomplex_t IT_0684 = (-2)*conj(IT_0224);
    const ccomplex_t IT_0685 = IT_0683 + IT_0684;
    const ccomplex_t IT_0686 = (-2)*IT_0402 + -IT_0444 + 2*IT_0592;
    const ccomplex_t IT_0687 = IT_0401 + IT_0686;
    const ccomplex_t IT_0688 = IT_0004*IT_0687;
    const ccomplex_t IT_0689 = -IT_0211 + (-0.5)*IT_0688;
    const ccomplex_t IT_0690 = s_12 + IT_0689;
    const ccomplex_t IT_0691 = m_N_4*IT_0690;
    const ccomplex_t IT_0692 = (-2)*IT_0691;
    const ccomplex_t IT_0693 = (-2)*conj(IT_0038);
    const ccomplex_t IT_0694 = (-2)*conj(IT_0081);
    const ccomplex_t IT_0695 = IT_0693 + IT_0694;
    const ccomplex_t IT_0696 = -IT_0444;
    const ccomplex_t IT_0697 = IT_0401 + IT_0696;
    const ccomplex_t IT_0698 = IT_0004*IT_0697;
    const ccomplex_t IT_0699 = -IT_0211 + (-0.25)*IT_0698;
    const ccomplex_t IT_0700 = s_12 + IT_0699;
    const ccomplex_t IT_0701 = m_N_4*IT_0700;
    const ccomplex_t IT_0702 = (-4)*IT_0701;
    const ccomplex_t IT_0703 = (-2)*conj(IT_0163);
    const ccomplex_t IT_0704 = (-2)*conj(IT_0165);
    const ccomplex_t IT_0705 = IT_0703 + IT_0704;
    const ccomplex_t IT_0706 = m_N_4*IT_0373;
    const ccomplex_t IT_0707 = (-4)*IT_0706;
    const ccomplex_t IT_0708 = (-2)*conj(IT_0205);
    const ccomplex_t IT_0709 = (-2)*conj(IT_0206);
    const ccomplex_t IT_0710 = IT_0708 + IT_0709;
    const ccomplex_t IT_0711 = s_12*IT_0650;
    const ccomplex_t IT_0712 = s_13*s_23*IT_0056;
    const ccomplex_t IT_0713 = s_14*s_23*s_34;
    const ccomplex_t IT_0714 = s_13*s_24*s_34;
    const ccomplex_t IT_0715 = s_12*IT_0295;
    const ccomplex_t IT_0716 = s_14*s_24*IT_0056;
    const ccomplex_t IT_0717 = (-2)*IT_0712 + 2*IT_0713 + 2*IT_0714 + (-2)
      *IT_0715 + (-2)*IT_0716;
    const ccomplex_t IT_0718 = IT_0711 + IT_0717;
    const ccomplex_t IT_0719 = IT_0001*IT_0718;
    const ccomplex_t IT_0720 = (-0.25)*IT_0719;
    const ccomplex_t IT_0721 = -IT_0211;
    const ccomplex_t IT_0722 = -IT_0210;
    const ccomplex_t IT_0723 = s_12 + IT_0720 + IT_0721 + IT_0722;
    const ccomplex_t IT_0724 = IT_0057*IT_0562*IT_0563;
    const ccomplex_t IT_0725 = IT_0075*IT_0557*IT_0558;
    const ccomplex_t IT_0726 = m_N_4*IT_0059;
    const ccomplex_t IT_0727 = m_N_4*IT_0077;
    const ccomplex_t IT_0728 = IT_0724 + IT_0725 + -IT_0726 + -IT_0727;
    const ccomplex_t IT_0729 = (-8)*conj(IT_0728);
    const ccomplex_t IT_0730 = 16*conj(IT_0567)*IT_0568 + (-8)*IT_0579*IT_0581
       + 4*IT_0588*IT_0591 + IT_0599*IT_0600 + IT_0607*IT_0608 + IT_0615*IT_0617
       + IT_0625*IT_0626 + IT_0635*IT_0637 + IT_0647*IT_0649 + IT_0656*IT_0665 +
       IT_0673*IT_0676 + IT_0682*IT_0685 + IT_0692*IT_0695 + IT_0702*IT_0705 +
       IT_0707*IT_0710 + IT_0723*IT_0729;
    const ccomplex_t IT_0731 = -IT_0356;
    const ccomplex_t IT_0732 = s_13 + IT_0731;
    const ccomplex_t IT_0733 = IT_0000*IT_0732;
    const ccomplex_t IT_0734 = -IT_0366;
    const ccomplex_t IT_0735 = s_23 + IT_0734;
    const ccomplex_t IT_0736 = IT_0000*IT_0735;
    const ccomplex_t IT_0737 = (-6)*IT_0038;
    const ccomplex_t IT_0738 = conj(IT_0206)*IT_0737;
    const ccomplex_t IT_0739 = (-6)*conj(IT_0038);
    const ccomplex_t IT_0740 = IT_0206*IT_0739;
    const ccomplex_t IT_0741 = IT_0163*conj(IT_0202);
    const ccomplex_t IT_0742 = (-6)*IT_0741;
    const ccomplex_t IT_0743 = conj(IT_0163)*IT_0202;
    const ccomplex_t IT_0744 = (-6)*IT_0743;
    const ccomplex_t IT_0745 = IT_0165*conj(IT_0202);
    const ccomplex_t IT_0746 = (-6)*IT_0745;
    const ccomplex_t IT_0747 = conj(IT_0165)*IT_0202;
    const ccomplex_t IT_0748 = (-6)*IT_0747;
    const ccomplex_t IT_0749 = (-6)*IT_0199;
    const ccomplex_t IT_0750 = conj(IT_0165)*IT_0749;
    const ccomplex_t IT_0751 = (-6)*conj(IT_0199);
    const ccomplex_t IT_0752 = IT_0165*IT_0751;
    const ccomplex_t IT_0753 = (-6)*IT_0205;
    const ccomplex_t IT_0754 = conj(IT_0081)*IT_0753;
    const ccomplex_t IT_0755 = (-6)*conj(IT_0205);
    const ccomplex_t IT_0756 = IT_0081*IT_0755;
    const ccomplex_t IT_0757 = (-6)*IT_0206;
    const ccomplex_t IT_0758 = conj(IT_0081)*IT_0757;
    const ccomplex_t IT_0759 = (-6)*conj(IT_0206);
    const ccomplex_t IT_0760 = IT_0081*IT_0759;
    const ccomplex_t IT_0761 = (IT_0733 + IT_0736)*(conj(IT_0205)*IT_0737 +
       IT_0738 + IT_0205*IT_0739 + IT_0740 + IT_0742 + IT_0744 + IT_0746 +
       IT_0748 + conj(IT_0163)*IT_0749 + IT_0750 + IT_0163*IT_0751 + IT_0752 +
       IT_0754 + IT_0756 + IT_0758 + IT_0760);
    const ccomplex_t IT_0762 = IT_0001*IT_0295;
    const ccomplex_t IT_0763 = 0.5*IT_0762;
    const ccomplex_t IT_0764 = 1 + IT_0763;
    const ccomplex_t IT_0765 = s_12*IT_0764;
    const ccomplex_t IT_0766 = (-4)*conj(IT_0120);
    const ccomplex_t IT_0767 = 0.25*IT_0766;
    const ccomplex_t IT_0768 = conj(IT_0122) + IT_0767;
    const ccomplex_t IT_0769 = (-4)*conj(IT_0183);
    const ccomplex_t IT_0770 = 0.25*IT_0769;
    const ccomplex_t IT_0771 = conj(IT_0185) + IT_0770;
    const ccomplex_t IT_0772 = 4*IT_0120;
    const ccomplex_t IT_0773 = 4*IT_0183;
    const ccomplex_t IT_0774 = (-4)*IT_0120;
    const ccomplex_t IT_0775 = (-4)*IT_0183;
    const ccomplex_t IT_0776 = IT_0000*IT_0764;
    const ccomplex_t IT_0777 = 0.25*IT_0774;
    const ccomplex_t IT_0778 = IT_0122 + IT_0777;
    const ccomplex_t IT_0779 = (-4)*IT_0122;
    const ccomplex_t IT_0780 = IT_0772 + IT_0779;
    const ccomplex_t IT_0781 = 4*conj(IT_0120);
    const ccomplex_t IT_0782 = (-4)*conj(IT_0122);
    const ccomplex_t IT_0783 = IT_0781 + IT_0782;
    const ccomplex_t IT_0784 = m_N_4*IT_0219;
    const ccomplex_t IT_0785 = IT_0221*IT_0562*IT_0659;
    const ccomplex_t IT_0786 = IT_0217*IT_0557*IT_0661;
    const ccomplex_t IT_0787 = -IT_0785 + -IT_0786;
    const ccomplex_t IT_0788 = m_N_4*IT_0144;
    const ccomplex_t IT_0789 = (-4)*IT_0185;
    const ccomplex_t IT_0790 = IT_0773 + IT_0789;
    const ccomplex_t IT_0791 = 4*conj(IT_0183);
    const ccomplex_t IT_0792 = (-4)*conj(IT_0185);
    const ccomplex_t IT_0793 = IT_0791 + IT_0792;
    const ccomplex_t IT_0794 = conj(IT_0567) + conj(IT_0663);
    const ccomplex_t IT_0795 = (IT_0183 + IT_0728)*IT_0768 + 0.25*(IT_0567 +
       IT_0663)*IT_0769 + (conj(IT_0183) + conj(IT_0728))*IT_0778 + conj(IT_0185
      )*(IT_0567 + IT_0663 + 0.25*IT_0772 + 0.25*IT_0779) + 0.25*conj(IT_0565)
      *IT_0780 + 0.25*IT_0565*IT_0783 + 0.25*(conj(IT_0657) + conj(IT_0784) +
       conj(IT_0787) + conj(IT_0788))*IT_0790 + 0.25*(IT_0657 + IT_0784 +
       IT_0787 + IT_0788)*IT_0793 + IT_0185*(0.25*IT_0781 + 0.25*IT_0782 +
       IT_0794) + 0.25*IT_0775*IT_0794;
    const ccomplex_t IT_0796 = pow(m_W, 6);
    const ccomplex_t IT_0797 = s_12*IT_0796;
    const ccomplex_t IT_0798 = s_14*s_24*IT_0650;
    const ccomplex_t IT_0799 = (-2)*IT_0798;
    const ccomplex_t IT_0800 = IT_0797 + IT_0799;
    const ccomplex_t IT_0801 = IT_0001*IT_0800;
    const ccomplex_t IT_0802 = 0.5*IT_0801;
    const ccomplex_t IT_0803 = (-4)*IT_0402;
    const ccomplex_t IT_0804 = IT_0401 + IT_0802 + IT_0803;
    const ccomplex_t IT_0805 = (-2)*IT_0712 + 4*IT_0714 + (-2)*IT_0716;
    const ccomplex_t IT_0806 = IT_0711 + IT_0805;
    const ccomplex_t IT_0807 = IT_0004*IT_0806;
    const ccomplex_t IT_0808 = -IT_0807;
    const ccomplex_t IT_0809 = IT_0804 + IT_0808;
    const ccomplex_t IT_0810 = (-4)*IT_0223;
    const ccomplex_t IT_0811 = (-0.25)*IT_0810;
    const ccomplex_t IT_0812 = IT_0219 + IT_0811;
    const ccomplex_t IT_0813 = s_34*IT_0056;
    const ccomplex_t IT_0814 = IT_0569*IT_0813;
    const ccomplex_t IT_0815 = IT_0004*IT_0652;
    const ccomplex_t IT_0816 = (-0.333333333333333)*IT_0814 + (
      -0.333333333333333)*IT_0815;
    const ccomplex_t IT_0817 = IT_0056 + IT_0816;
    const ccomplex_t IT_0818 = IT_0000*IT_0817;
    const ccomplex_t IT_0819 = 3*IT_0818;
    const ccomplex_t IT_0820 = 2*IT_0494;
    const ccomplex_t IT_0821 = IT_0483 + IT_0820;
    const ccomplex_t IT_0822 = m_N_4*IT_0004;
    const ccomplex_t IT_0823 = s_23*IT_0650;
    const ccomplex_t IT_0824 = IT_0004*IT_0823;
    const ccomplex_t IT_0825 = (-2)*IT_0368 + (-0.5)*IT_0824;
    const ccomplex_t IT_0826 = IT_0367 + IT_0825;
    const ccomplex_t IT_0827 = IT_0822*IT_0826;
    const ccomplex_t IT_0828 = (-2)*IT_0827;
    const ccomplex_t IT_0829 = 2*IT_0657;
    const ccomplex_t IT_0830 = 2*IT_0788;
    const ccomplex_t IT_0831 = IT_0829 + IT_0830;
    const ccomplex_t IT_0832 = (-2)*IT_0663;
    const ccomplex_t IT_0833 = IT_0831 + IT_0832;
    const ccomplex_t IT_0834 = (-2)*IT_0368;
    const ccomplex_t IT_0835 = IT_0367 + IT_0834;
    const ccomplex_t IT_0836 = IT_0004*IT_0835;
    const ccomplex_t IT_0837 = 0.5*IT_0609 + (-0.5)*IT_0611 + (-0.5)*IT_0836;
    const ccomplex_t IT_0838 = s_23 + IT_0837;
    const ccomplex_t IT_0839 = m_N_4*IT_0838;
    const ccomplex_t IT_0840 = 2*IT_0839;
    const ccomplex_t IT_0841 = IT_0569*IT_0639;
    const ccomplex_t IT_0842 = (-2)*IT_0358;
    const ccomplex_t IT_0843 = IT_0357 + IT_0842;
    const ccomplex_t IT_0844 = IT_0004*IT_0843;
    const ccomplex_t IT_0845 = 0.5*IT_0638 + (-0.5)*IT_0841 + (-0.5)*IT_0844;
    const ccomplex_t IT_0846 = s_13 + IT_0845;
    const ccomplex_t IT_0847 = m_N_4*IT_0846;
    const ccomplex_t IT_0848 = 2*IT_0847;
    const ccomplex_t IT_0849 = s_24*s_34*IT_0056;
    const ccomplex_t IT_0850 = (-2)*IT_0849;
    const ccomplex_t IT_0851 = IT_0823 + IT_0850;
    const ccomplex_t IT_0852 = IT_0001*IT_0851;
    const ccomplex_t IT_0853 = (-0.666666666666667)*IT_0836 +
       0.166666666666667*IT_0852;
    const ccomplex_t IT_0854 = s_23 + IT_0853;
    const ccomplex_t IT_0855 = m_N_4*IT_0854;
    const ccomplex_t IT_0856 = 6*IT_0855;
    const ccomplex_t IT_0857 = 2*IT_0567;
    const ccomplex_t IT_0858 = (-2)*IT_0784;
    const ccomplex_t IT_0859 = (-2)*IT_0787;
    const ccomplex_t IT_0860 = IT_0857 + IT_0858 + IT_0859;
    const ccomplex_t IT_0861 = -IT_0295;
    const ccomplex_t IT_0862 = IT_0650 + IT_0861;
    const ccomplex_t IT_0863 = IT_0004*IT_0862;
    const ccomplex_t IT_0864 = -IT_0056 + 0.25*IT_0863;
    const ccomplex_t IT_0865 = IT_0400 + IT_0864;
    const ccomplex_t IT_0866 = IT_0000*IT_0865;
    const ccomplex_t IT_0867 = 4*IT_0866;
    const ccomplex_t IT_0868 = -IT_0287;
    const ccomplex_t IT_0869 = IT_0285 + IT_0868;
    const ccomplex_t IT_0870 = s_13*s_14*s_24;
    const ccomplex_t IT_0871 = s_23*IT_0000*IT_0056;
    const ccomplex_t IT_0872 = s_12*s_13*IT_0056;
    const ccomplex_t IT_0873 = s_24*s_34*IT_0000;
    const ccomplex_t IT_0874 = 0.5*IT_0871 + (-0.5)*IT_0872 + -IT_0873;
    const ccomplex_t IT_0875 = IT_0870 + IT_0874;
    const ccomplex_t IT_0876 = IT_0004*IT_0875;
    const ccomplex_t IT_0877 = (-2)*IT_0876;
    const ccomplex_t IT_0878 = (-0.5)*IT_0877;
    const ccomplex_t IT_0879 = IT_0869 + IT_0878;
    const ccomplex_t IT_0880 = 4*IT_0220;
    const ccomplex_t IT_0881 = 4*IT_0224;
    const ccomplex_t IT_0882 = IT_0880 + IT_0881;
    const ccomplex_t IT_0883 = (-2)*IT_0712;
    const ccomplex_t IT_0884 = IT_0711 + IT_0883;
    const ccomplex_t IT_0885 = IT_0004*IT_0884;
    const ccomplex_t IT_0886 = 0.5*IT_0885;
    const ccomplex_t IT_0887 = (-2)*IT_0712 + 2*IT_0713 + 2*IT_0714 + (-2)
      *IT_0715;
    const ccomplex_t IT_0888 = IT_0711 + IT_0887;
    const ccomplex_t IT_0889 = IT_0004*IT_0888;
    const ccomplex_t IT_0890 = (-0.5)*IT_0889;
    const ccomplex_t IT_0891 = s_14*s_23*IT_0056;
    const ccomplex_t IT_0892 = s_13*s_24*IT_0056;
    const ccomplex_t IT_0893 = IT_0891 + IT_0892;
    const ccomplex_t IT_0894 = s_12*s_34*IT_0056;
    const ccomplex_t IT_0895 = -IT_0894;
    const ccomplex_t IT_0896 = IT_0893 + IT_0895;
    const ccomplex_t IT_0897 = IT_0569*IT_0896;
    const ccomplex_t IT_0898 = 0.5*IT_0897;
    const ccomplex_t IT_0899 = IT_0401 + IT_0403 + IT_0886 + IT_0890 + IT_0898;
    const ccomplex_t IT_0900 = 4*IT_0476;
    const ccomplex_t IT_0901 = 4*IT_0489;
    const ccomplex_t IT_0902 = IT_0900 + IT_0901;
    const ccomplex_t IT_0903 = (-0.5)*IT_0278;
    const ccomplex_t IT_0904 = s_23*IT_0250;
    const ccomplex_t IT_0905 = IT_0870 + IT_0904;
    const ccomplex_t IT_0906 = s_12*s_14*s_34;
    const ccomplex_t IT_0907 = -IT_0871 + -IT_0906;
    const ccomplex_t IT_0908 = IT_0905 + IT_0907;
    const ccomplex_t IT_0909 = IT_0004*IT_0908;
    const ccomplex_t IT_0910 = (-0.25)*IT_0909;
    const ccomplex_t IT_0911 = (-0.5)*IT_0287;
    const ccomplex_t IT_0912 = (-0.5)*IT_0289;
    const ccomplex_t IT_0913 = 0.5*IT_0291;
    const ccomplex_t IT_0914 = IT_0285 + IT_0903 + IT_0910 + IT_0911 + IT_0912
       + IT_0913;
    const ccomplex_t IT_0915 = (-8)*IT_0038;
    const ccomplex_t IT_0916 = (-8)*IT_0081;
    const ccomplex_t IT_0917 = IT_0915 + IT_0916;
    const ccomplex_t IT_0918 = -IT_0281;
    const ccomplex_t IT_0919 = IT_0279 + IT_0918;
    const ccomplex_t IT_0920 = s_13*IT_0007;
    const ccomplex_t IT_0921 = s_12*s_24*s_34;
    const ccomplex_t IT_0922 = s_13*IT_0000*IT_0056;
    const ccomplex_t IT_0923 = s_12*s_23*IT_0056;
    const ccomplex_t IT_0924 = -IT_0921 + (-0.5)*IT_0922 + 0.5*IT_0923;
    const ccomplex_t IT_0925 = IT_0920 + IT_0924;
    const ccomplex_t IT_0926 = IT_0004*IT_0925;
    const ccomplex_t IT_0927 = 2*IT_0926;
    const ccomplex_t IT_0928 = (-0.5)*IT_0927;
    const ccomplex_t IT_0929 = IT_0919 + IT_0928;
    const ccomplex_t IT_0930 = (-4)*IT_0033;
    const ccomplex_t IT_0931 = (-4)*IT_0078;
    const ccomplex_t IT_0932 = IT_0930 + IT_0931;
    const ccomplex_t IT_0933 = -IT_0282;
    const ccomplex_t IT_0934 = s_14*s_23*s_24;
    const ccomplex_t IT_0935 = IT_0920 + IT_0934;
    const ccomplex_t IT_0936 = -IT_0921 + -IT_0922;
    const ccomplex_t IT_0937 = IT_0935 + IT_0936;
    const ccomplex_t IT_0938 = IT_0004*IT_0937;
    const ccomplex_t IT_0939 = (-0.5)*IT_0938;
    const ccomplex_t IT_0940 = IT_0281 + IT_0299 + IT_0300 + IT_0933 + IT_0939;
    const ccomplex_t IT_0941 = (-4)*IT_0038;
    const ccomplex_t IT_0942 = (-4)*IT_0081;
    const ccomplex_t IT_0943 = IT_0941 + IT_0942;
    const ccomplex_t IT_0944 = s_14*s_23*s_34*IT_0004;
    const ccomplex_t IT_0945 = s_12*IT_0004*IT_0295;
    const ccomplex_t IT_0946 = -IT_0945;
    const ccomplex_t IT_0947 = IT_0711 + IT_0713 + IT_0714;
    const ccomplex_t IT_0948 = (-2)*IT_0712 + -IT_0715;
    const ccomplex_t IT_0949 = IT_0947 + IT_0948;
    const ccomplex_t IT_0950 = IT_0004*IT_0949;
    const ccomplex_t IT_0951 = (-0.5)*IT_0950;
    const ccomplex_t IT_0952 = -IT_0592;
    const ccomplex_t IT_0953 = IT_0401 + IT_0944 + IT_0946 + IT_0951 + IT_0952;
    const ccomplex_t IT_0954 = (-4)*IT_0199;
    const ccomplex_t IT_0955 = (-4)*IT_0202;
    const ccomplex_t IT_0956 = IT_0954 + IT_0955;
    const ccomplex_t IT_0957 = 0.25*IT_0801;
    const ccomplex_t IT_0958 = IT_0401 + IT_0957;
    const ccomplex_t IT_0959 = (-8)*IT_0382;
    const ccomplex_t IT_0960 = s_13*IT_0650;
    const ccomplex_t IT_0961 = s_14*s_34*IT_0056;
    const ccomplex_t IT_0962 = (-2)*IT_0961;
    const ccomplex_t IT_0963 = IT_0960 + IT_0962;
    const ccomplex_t IT_0964 = IT_0001*IT_0963;
    const ccomplex_t IT_0965 = (-0.666666666666667)*IT_0356 + (
      -0.166666666666667)*IT_0964;
    const ccomplex_t IT_0966 = s_13 + IT_0965;
    const ccomplex_t IT_0967 = m_N_4*IT_0966;
    const ccomplex_t IT_0968 = (-6)*IT_0967;
    const ccomplex_t IT_0969 = (-2)*IT_0728;
    const ccomplex_t IT_0970 = IT_0001*IT_0960;
    const ccomplex_t IT_0971 = 0.666666666666667*IT_0638 + (-0.166666666666667
      )*IT_0970;
    const ccomplex_t IT_0972 = s_13 + IT_0971;
    const ccomplex_t IT_0973 = m_N_4*IT_0972;
    const ccomplex_t IT_0974 = (-6)*IT_0973;
    const ccomplex_t IT_0975 = (-2)*IT_0565;
    const ccomplex_t IT_0976 = 9*IT_0582;
    const ccomplex_t IT_0977 = (-2)*IT_0376;
    const ccomplex_t IT_0978 = IT_0000*IT_0004;
    const ccomplex_t IT_0979 = (-0.25)*IT_0650;
    const ccomplex_t IT_0980 = IT_0295 + IT_0979;
    const ccomplex_t IT_0981 = IT_0978*IT_0980;
    const ccomplex_t IT_0982 = (-4)*IT_0981;
    const ccomplex_t IT_0983 = (-2)*IT_0378;
    const ccomplex_t IT_0984 = IT_0809*IT_0812 + (-0.25)*IT_0819*IT_0821 + (
      -0.25)*IT_0828*IT_0833 + (-0.25)*IT_0124*IT_0840 + (-0.25)*IT_0242*IT_0848
       + (-0.25)*IT_0856*IT_0860 + (-0.25)*IT_0482*IT_0867 + (-0.25)*IT_0879
      *IT_0882 + (-0.25)*IT_0899*IT_0902 + (-0.25)*IT_0914*IT_0917 + (-0.25)
      *IT_0929*IT_0932 + (-0.25)*IT_0940*IT_0943 + (-0.25)*IT_0953*IT_0956 + (
      -0.25)*IT_0958*IT_0959 + (-0.25)*IT_0968*IT_0969 + (-0.25)*IT_0974*IT_0975
       + (-0.25)*IT_0976*IT_0977 + (-0.25)*IT_0982*IT_0983;
    const ccomplex_t IT_0985 = s_12 + IT_0721;
    const ccomplex_t IT_0986 = IT_0000*IT_0985;
    const ccomplex_t IT_0987 = 0.5*IT_0007 + (-0.5)*IT_0582 + -IT_0666;
    const ccomplex_t IT_0988 = IT_0005 + IT_0987;
    const ccomplex_t IT_0989 = IT_0004*IT_0988;
    const ccomplex_t IT_0990 = IT_0000 + IT_0540 + IT_0989;
    const ccomplex_t IT_0991 = m_N_4*IT_0990;
    const ccomplex_t IT_0992 = 4*IT_0081;
    const ccomplex_t IT_0993 = IT_0535 + IT_0992;
    const ccomplex_t IT_0994 = IT_0250 + IT_0618;
    const ccomplex_t IT_0995 = IT_0004*IT_0994;
    const ccomplex_t IT_0996 = IT_0251 + IT_0621;
    const ccomplex_t IT_0997 = (-4)*IT_0996;
    const ccomplex_t IT_0998 = IT_0995 + IT_0997;
    const ccomplex_t IT_0999 = m_N_4*IT_0998;
    const ccomplex_t IT_1000 = 2*IT_0728;
    const ccomplex_t IT_1001 = (-2)*IT_0402 + -IT_0444 + 2*IT_0944;
    const ccomplex_t IT_1002 = IT_0401 + IT_1001;
    const ccomplex_t IT_1003 = IT_0004*IT_1002;
    const ccomplex_t IT_1004 = -IT_0211 + (-0.5)*IT_1003;
    const ccomplex_t IT_1005 = s_12 + IT_1004;
    const ccomplex_t IT_1006 = m_N_4*IT_1005;
    const ccomplex_t IT_1007 = (-2)*IT_1006;
    const ccomplex_t IT_1008 = IT_0000*IT_0213;
    const ccomplex_t IT_1009 = (-2)*IT_0239;
    const ccomplex_t IT_1010 = IT_0000*IT_0258;
    const ccomplex_t IT_1011 = (-2)*IT_0162;
    const ccomplex_t IT_1012 = (-2)*IT_0567;
    const ccomplex_t IT_1013 = IT_0269*IT_0491 + IT_0267*IT_0508 + IT_0702*
      (IT_0857 + IT_0969) + ((-12)*IT_0163 + (-12)*IT_0165)*IT_0986 + 4*IT_0663
      *IT_0991 + IT_0534*IT_0993 + IT_0999*IT_1000 + IT_0832*IT_1007 + IT_1008
      *IT_1009 + IT_1010*IT_1011 + IT_0625*IT_1012;
    const ccomplex_t IT_1014 = 2*conj(IT_0728);
    const ccomplex_t IT_1015 = (-2)*conj(IT_0663);
    const ccomplex_t IT_1016 = 2*conj(IT_0567);
    const ccomplex_t IT_1017 = (-2)*conj(IT_0728);
    const ccomplex_t IT_1018 = (-2)*conj(IT_0239);
    const ccomplex_t IT_1019 = (-2)*conj(IT_0567);
    const ccomplex_t IT_1020 = (IT_0000 + -IT_0006)*(s_12 + (-0.5)*IT_0447);
    const ccomplex_t IT_1021 = 2*IT_1020;
    const ccomplex_t IT_1022 = (-2)*IT_0078;
    const ccomplex_t IT_1023 = (-2)*IT_0712 + 4*IT_0713 + (-2)*IT_0716;
    const ccomplex_t IT_1024 = IT_0711 + IT_1023;
    const ccomplex_t IT_1025 = IT_0001*IT_1024;
    const ccomplex_t IT_1026 = 0.25*IT_1025;
    const ccomplex_t IT_1027 = (-0.5)*IT_0447;
    const ccomplex_t IT_1028 = (-0.5)*IT_0405;
    const ccomplex_t IT_1029 = s_12 + IT_1026 + IT_1027 + IT_1028;
    const ccomplex_t IT_1030 = IT_0570 + IT_0572;
    const ccomplex_t IT_1031 = -IT_0571;
    const ccomplex_t IT_1032 = IT_1030 + IT_1031;
    const ccomplex_t IT_1033 = IT_0569*IT_1032;
    const ccomplex_t IT_1034 = 0.25*IT_1033;
    const ccomplex_t IT_1035 = s_12 + IT_0578 + IT_1034;
    const ccomplex_t IT_1036 = (-8)*IT_0122;
    const ccomplex_t IT_1037 = 0.125*IT_1036;
    const ccomplex_t IT_1038 = IT_0120 + IT_1037;
    const ccomplex_t IT_1039 = IT_0367 + IT_0368;
    const ccomplex_t IT_1040 = s_23*IT_0004*IT_0295;
    const ccomplex_t IT_1041 = (-2)*IT_1040;
    const ccomplex_t IT_1042 = IT_1039 + IT_1041;
    const ccomplex_t IT_1043 = IT_0004*IT_1042;
    const ccomplex_t IT_1044 = -IT_0366 + 0.5*IT_1043;
    const ccomplex_t IT_1045 = s_23 + IT_1044;
    const ccomplex_t IT_1046 = m_N_4*IT_1045;
    const ccomplex_t IT_1047 = 2*IT_1046;
    const ccomplex_t IT_1048 = s_23*s_34;
    const ccomplex_t IT_1049 = (-2)*IT_1048;
    const ccomplex_t IT_1050 = IT_0610 + IT_1049;
    const ccomplex_t IT_1051 = IT_0569*IT_1050;
    const ccomplex_t IT_1052 = (-0.5)*IT_0609 + (-0.25)*IT_1051;
    const ccomplex_t IT_1053 = s_23 + IT_1052;
    const ccomplex_t IT_1054 = m_N_4*IT_1053;
    const ccomplex_t IT_1055 = (-4)*IT_1054;
    const ccomplex_t IT_1056 = (-2)*IT_0038;
    const ccomplex_t IT_1057 = IT_0005 + IT_0618;
    const ccomplex_t IT_1058 = IT_0004*IT_1057;
    const ccomplex_t IT_1059 = IT_0006 + IT_0621;
    const ccomplex_t IT_1060 = (-4)*IT_1059;
    const ccomplex_t IT_1061 = IT_1058 + IT_1060;
    const ccomplex_t IT_1062 = m_N_4*IT_1061;
    const ccomplex_t IT_1063 = 0.125*IT_0469*IT_0607 + IT_0567*IT_0723 + 0.5*
      (IT_0081 + 0.25*IT_0535)*IT_0991 + IT_0663*IT_1029 + IT_1035*IT_1038 +
       0.125*IT_0408*IT_1047 + 0.125*IT_0477*IT_1055 + 0.125*IT_1007*IT_1056 +
       0.125*IT_1022*IT_1062;
    const ccomplex_t IT_1064 = 8*conj(IT_0120);
    const ccomplex_t IT_1065 = (-8)*conj(IT_0122);
    const ccomplex_t IT_1066 = IT_1064 + IT_1065;
    const ccomplex_t IT_1067 = (-2)*conj(IT_0078);
    const ccomplex_t IT_1068 = (-2)*IT_0164;
    const ccomplex_t IT_1069 = conj(IT_0164)*IT_1068;
    const ccomplex_t IT_1070 = conj(IT_0164)*IT_1056;
    const ccomplex_t IT_1071 = IT_0164*IT_0693;
    const ccomplex_t IT_1072 = conj(IT_0081)*IT_1068;
    const ccomplex_t IT_1073 = IT_0081*IT_0675;
    const ccomplex_t IT_1074 = (IT_0000 + -IT_0251)*(s_12 + (-0.5)*IT_0405);
    const ccomplex_t IT_1075 = 2*IT_1074;
    const ccomplex_t IT_1076 = (IT_0000 + -IT_0008)*(s_12 + (-0.5)*IT_0405);
    const ccomplex_t IT_1077 = 2*IT_1076;
    const ccomplex_t IT_1078 = conj(IT_0038)*IT_1056;
    const ccomplex_t IT_1079 = conj(IT_0081)*IT_1056;
    const ccomplex_t IT_1080 = IT_0081*IT_0693;
    const ccomplex_t IT_1081 = (-2)*IT_0163;
    const ccomplex_t IT_1082 = conj(IT_0163)*IT_1081;
    const ccomplex_t IT_1083 = conj(IT_0165)*IT_1081;
    const ccomplex_t IT_1084 = IT_0165*IT_0703;
    const ccomplex_t IT_1085 = (-2)*IT_0165;
    const ccomplex_t IT_1086 = conj(IT_0165)*IT_1085;
    const ccomplex_t IT_1087 = (-2)*IT_0081;
    const ccomplex_t IT_1088 = conj(IT_0081)*IT_1087;
    const ccomplex_t IT_1089 = (IT_0000 + -IT_0249)*(s_12 + (-0.5)*IT_0447);
    const ccomplex_t IT_1090 = 2*IT_1089;
    const ccomplex_t IT_1091 = (-2)*IT_0224;
    const ccomplex_t IT_1092 = 2*IT_0202;
    const ccomplex_t IT_1093 = IT_0408 + IT_1092;
    const ccomplex_t IT_1094 = 2*IT_0663;
    const ccomplex_t IT_1095 = -IT_0400;
    const ccomplex_t IT_1096 = IT_0056 + IT_1095;
    const ccomplex_t IT_1097 = IT_0000*IT_1096;
    const ccomplex_t IT_1098 = m_N_4*IT_0363;
    const ccomplex_t IT_1099 = (-4)*IT_1098;
    const ccomplex_t IT_1100 = s_34*IT_0004;
    const ccomplex_t IT_1101 = IT_0570 + IT_0571;
    const ccomplex_t IT_1102 = s_13*s_23*s_34*IT_0004;
    const ccomplex_t IT_1103 = (-2)*IT_1102;
    const ccomplex_t IT_1104 = IT_1101 + IT_1103;
    const ccomplex_t IT_1105 = IT_1100*IT_1104;
    const ccomplex_t IT_1106 = (-2)*IT_0476;
    const ccomplex_t IT_1107 = 2*conj(IT_0663);
    const ccomplex_t IT_1108 = (-2)*conj(IT_0476);
    const ccomplex_t IT_1109 = s_13*m_N_4;
    const ccomplex_t IT_1110 = IT_0764*IT_1109;
    const ccomplex_t IT_1111 = IT_0282 + IT_0284 + IT_0299 + IT_0301 + IT_0303
       + IT_0305 + IT_0409;
    const ccomplex_t IT_1112 = IT_0278 + IT_0283 + IT_0288 + IT_0290 + IT_0292
       + IT_0294 + IT_0297;
    const ccomplex_t IT_1113 = (-2)*conj(IT_0199);
    const ccomplex_t IT_1114 = s_23*m_N_4;
    const ccomplex_t IT_1115 = IT_0764*IT_1114;
    const ccomplex_t IT_1116 = 4*IT_1115;
    const ccomplex_t IT_1117 = IT_0582*IT_0764;
    const ccomplex_t IT_1118 = (-4)*IT_0476;
    const ccomplex_t IT_1119 = (-0.25)*IT_1118;
    const ccomplex_t IT_1120 = IT_0489 + IT_1119;
    const ccomplex_t IT_1121 = IT_0274 + IT_0508;
    const ccomplex_t IT_1122 = IT_0274 + IT_0491;
    const ccomplex_t IT_1123 = (IT_0401 + (-2)*IT_0402)*(1 + 0.5*IT_0762);
    const ccomplex_t IT_1124 = (-4)*IT_0494;
    const ccomplex_t IT_1125 = (-4)*IT_0309;
    const ccomplex_t IT_1126 = IT_1124 + IT_1125;
    const ccomplex_t IT_1127 = (-2)*IT_0205;
    const ccomplex_t IT_1128 = (-2)*IT_0206;
    const ccomplex_t IT_1129 = IT_1127 + IT_1128;
    const ccomplex_t IT_1130 = (-0.5)*IT_0638 + 0.25*IT_0841;
    const ccomplex_t IT_1131 = s_13 + IT_1130;
    const ccomplex_t IT_1132 = m_N_4*IT_1131;
    const ccomplex_t IT_1133 = 4*IT_1132;
    const ccomplex_t IT_1134 = 0.25*IT_0775;
    const ccomplex_t IT_1135 = IT_0185 + IT_1134;
    const ccomplex_t IT_1136 = (-0.25)*IT_0647*IT_1012 + (-0.25)*IT_1000
      *IT_1055 + (-0.25)*IT_0780*IT_1110 + IT_1117*IT_1120 + (-0.25)*IT_1112
      *IT_1121 + (-0.25)*IT_1111*IT_1122 + (-0.25)*IT_1123*IT_1126 + (-0.25)
      *IT_1105*IT_1129 + (-0.25)*IT_1094*IT_1133 + (-0.25)*IT_1116*IT_1135;
    const ccomplex_t IT_1137 = IT_0276 + IT_0493;
    const ccomplex_t IT_1138 = (-4)*conj(IT_0476);
    const ccomplex_t IT_1139 = (-4)*conj(IT_0309);
    const ccomplex_t IT_1140 = IT_0494*(IT_0647*IT_1019 + IT_1014*IT_1055 +
       IT_0710*IT_1105 + IT_0783*IT_1110 + IT_0600*IT_1112 + IT_0771*IT_1116 +
       IT_1107*IT_1133 + IT_1111*IT_1137 + IT_1117*IT_1138 + IT_1123*IT_1139);
    const ccomplex_t IT_1141 = IT_0183 + -IT_0185;
    const ccomplex_t IT_1142 = -IT_0402 + (-2)*IT_0444 + 2*IT_0944;
    const ccomplex_t IT_1143 = IT_0401 + IT_1142;
    const ccomplex_t IT_1144 = IT_0004*IT_1143;
    const ccomplex_t IT_1145 = -IT_0210 + (-0.5)*IT_1144;
    const ccomplex_t IT_1146 = s_12 + IT_1145;
    const ccomplex_t IT_1147 = m_N_4*IT_1146;
    const ccomplex_t IT_1148 = (-2)*IT_1147;
    const ccomplex_t IT_1149 = 8*IT_0567;
    const ccomplex_t IT_1150 = (-2)*IT_0199;
    const ccomplex_t IT_1151 = 2*IT_0250 + -IT_0582 + (-2)*IT_0583;
    const ccomplex_t IT_1152 = IT_0248 + IT_1151;
    const ccomplex_t IT_1153 = IT_0004*IT_1152;
    const ccomplex_t IT_1154 = -IT_0249 + 0.5*IT_1153;
    const ccomplex_t IT_1155 = IT_0000 + IT_1154;
    const ccomplex_t IT_1156 = m_N_4*IT_1155;
    const ccomplex_t IT_1157 = 2*IT_1156;
    const ccomplex_t IT_1158 = 4*IT_0163;
    const ccomplex_t IT_1159 = 4*IT_0991;
    const ccomplex_t IT_1160 = IT_0274*IT_0607 + (-16)*IT_0568*IT_0663 +
       IT_0270*IT_0999 + IT_0417*IT_1047 + IT_0483*IT_1055 + IT_0702*IT_1056 +
       IT_1062*IT_1068 + IT_1007*IT_1081 + IT_0477*IT_1133 + 8*IT_1035*IT_1141 +
       IT_0336*IT_1148 + IT_0656*IT_1149 + IT_1099*IT_1150 + IT_1091*IT_1157 + 
      (IT_0165 + 0.25*IT_1158)*IT_1159;
    const ccomplex_t IT_1161 = IT_0336 + IT_0495;
    const ccomplex_t IT_1162 = IT_0340 + IT_0348;
    const ccomplex_t IT_1163 = 2*IT_0489;
    const ccomplex_t IT_1164 = IT_0477 + IT_1163;
    const ccomplex_t IT_1165 = IT_1011 + IT_1068;
    const ccomplex_t IT_1166 = (-2)*IT_0220;
    const ccomplex_t IT_1167 = IT_1091 + IT_1166;
    const ccomplex_t IT_1168 = IT_1056 + IT_1087;
    const ccomplex_t IT_1169 = IT_1081 + IT_1085;
    const ccomplex_t IT_1170 = 8*IT_0657;
    const ccomplex_t IT_1171 = 8*IT_0788;
    const ccomplex_t IT_1172 = (-8)*IT_0663;
    const ccomplex_t IT_1173 = (-8)*IT_0728;
    const ccomplex_t IT_1174 = IT_0001*IT_0806;
    const ccomplex_t IT_1175 = 0.25*IT_1174;
    const ccomplex_t IT_1176 = s_12 + IT_1027 + IT_1028 + IT_1175;
    const ccomplex_t IT_1177 = 8*IT_0565;
    const ccomplex_t IT_1178 = IT_0568*(16*IT_0567 + (-16)*IT_0784 + (-16)
      *IT_0787) + IT_0615*IT_0821 + IT_0588*IT_0993 + IT_0635*IT_1093 + IT_0599
      *IT_1121 + IT_0707*IT_1129 + (-8)*IT_0579*IT_1141 + IT_0607*IT_1161 +
       IT_0625*IT_1162 + IT_0647*IT_1164 + IT_0673*IT_1165 + IT_0682*IT_1167 +
       IT_0692*IT_1168 + IT_0702*IT_1169 + IT_0656*(IT_1170 + IT_1171 + IT_1172)
       + IT_0723*IT_1173 + IT_1176*IT_1177;
    const ccomplex_t IT_1179 = 0.125*IT_0336*IT_0673 + 0.125*IT_0270*IT_0702 +
       0.125*IT_0408*IT_0707 + IT_0663*IT_0723 + IT_0579*IT_1038 + 0.125*IT_0625
      *IT_1056 + 0.125*IT_0599*IT_1091 + 0.125*IT_0615*IT_1106 + 0.125*IT_1149
      *IT_1176;
    const ccomplex_t IT_1180 = 8*conj(IT_0567);
    const ccomplex_t IT_1181 = 4*conj(IT_0163);
    const ccomplex_t IT_1182 = 0.25*IT_1181;
    const ccomplex_t IT_1183 = conj(IT_0165) + IT_1182;
    const ccomplex_t IT_1184 = IT_0276*IT_0607 + (-16)*IT_0568*conj(IT_0663) +
       IT_0693*IT_0702 + IT_0272*IT_0999 + IT_0703*IT_1007 + 8*conj(IT_0728)
      *IT_1029 + 8*IT_0581*IT_1035 + IT_0419*IT_1047 + IT_0487*IT_1055 + IT_0675
      *IT_1062 + IT_1099*IT_1113 + IT_0479*IT_1133 + IT_0338*IT_1148 + IT_0684
      *IT_1157 + IT_0656*IT_1180 + IT_1159*IT_1183;
    const ccomplex_t IT_1185 = -IT_0249;
    const ccomplex_t IT_1186 = IT_0000 + IT_1185;
    const ccomplex_t IT_1187 = IT_0000*IT_1186;
    const ccomplex_t IT_1188 = -IT_1187;
    const ccomplex_t IT_1189 = 6*conj(IT_0164);
    const ccomplex_t IT_1190 = (-6)*conj(IT_0164);
    const ccomplex_t IT_1191 = s_12 + IT_0722;
    const ccomplex_t IT_1192 = IT_0000*IT_1191;
    const ccomplex_t IT_1193 = (-6)*IT_1192;
    const ccomplex_t IT_1194 = IT_0338*IT_0449 + IT_0354*IT_0464 + IT_0267
      *IT_0626 + IT_0423*IT_0637 + IT_0705*IT_1008 + IT_0599*IT_1019 + IT_0684
      *IT_1090 + IT_0607*IT_1107 + IT_0479*IT_1112 + IT_1017*IT_1157 + IT_1188
      *IT_1189 + (conj(IT_0162) + (-0.166666666666667)*IT_1190)*IT_1193;
    const ccomplex_t IT_1195 = -IT_0006;
    const ccomplex_t IT_1196 = IT_0000 + IT_1195;
    const ccomplex_t IT_1197 = IT_0000*IT_1196;
    const ccomplex_t IT_1198 = -IT_1197;
    const ccomplex_t IT_1199 = IT_0450 + IT_0511;
    const ccomplex_t IT_1200 = 0.166666666666667*IT_0449;
    const ccomplex_t IT_1201 = IT_0269*IT_0626 + IT_0412*IT_0637 + IT_0705
      *IT_1010 + IT_0673*IT_1016 + IT_1015*IT_1062 + IT_1021*IT_1067 + IT_0479
      *IT_1111 + IT_1014*IT_1148 + IT_1190*IT_1192 + 6*(conj(IT_0162) +
       0.166666666666667*IT_1189)*IT_1198 + 6*IT_1199*IT_1200;
    const ccomplex_t IT_1202 = IT_1009 + IT_1068;
    const ccomplex_t IT_1203 = 6*IT_0224;
    const ccomplex_t IT_1204 = (-6)*IT_0078;
    const ccomplex_t IT_1205 = IT_0267*IT_0270 + IT_0423*IT_0482 + IT_0682
      *IT_0857 + IT_0607*IT_1000 + IT_1008*IT_1056 + IT_0483*IT_1112 + IT_0832
      *IT_1157 + IT_0464*IT_1162 + 6*IT_0274*IT_1200 + IT_1090*IT_1202 + 6
      *IT_1188*(IT_0220 + 0.166666666666667*IT_1203) + (-6)*IT_1192*(IT_0033 + (
      -0.166666666666667)*IT_1204);
    const ccomplex_t IT_1206 = 0.25*IT_1110;
    const ccomplex_t IT_1207 = 0.25*IT_1105;
    const ccomplex_t IT_1208 = 0.25*IT_1123;
    const ccomplex_t IT_1209 = 0.25*IT_1117;
    const ccomplex_t IT_1210 = 0.25*IT_1133;
    const ccomplex_t IT_1211 = 0.25*IT_0647;
    const ccomplex_t IT_1212 = IT_0695*IT_0702;
    const ccomplex_t IT_1213 = 0.25*IT_1212;
    const ccomplex_t IT_1214 = IT_0272 + IT_0354;
    const ccomplex_t IT_1215 = IT_0625*IT_1214;
    const ccomplex_t IT_1216 = 0.25*IT_1215;
    const ccomplex_t IT_1217 = IT_0692*IT_0705;
    const ccomplex_t IT_1218 = 0.25*IT_1217;
    const ccomplex_t IT_1219 = (-2)*conj(IT_0202);
    const ccomplex_t IT_1220 = IT_1113 + IT_1219;
    const ccomplex_t IT_1221 = IT_0707*IT_1220;
    const ccomplex_t IT_1222 = 0.25*IT_1221;
    const ccomplex_t IT_1223 = IT_0486*IT_0635;
    const ccomplex_t IT_1224 = 0.25*IT_1223;
    const ccomplex_t IT_1225 = IT_0607*IT_1137;
    const ccomplex_t IT_1226 = 0.25*IT_1225;
    const ccomplex_t IT_1227 = IT_0599*IT_1199;
    const ccomplex_t IT_1228 = 0.25*IT_1227;
    const ccomplex_t IT_1229 = IT_0675 + IT_1018;
    const ccomplex_t IT_1230 = IT_0682*IT_1229;
    const ccomplex_t IT_1231 = 0.25*IT_1230;
    const ccomplex_t IT_1232 = (-2)*conj(IT_0033);
    const ccomplex_t IT_1233 = IT_1067 + IT_1232;
    const ccomplex_t IT_1234 = IT_0673*IT_1233;
    const ccomplex_t IT_1235 = 0.25*IT_1234;
    const ccomplex_t IT_1236 = 0.25*IT_0615;
    const ccomplex_t IT_1237 = IT_0656*IT_0729;
    const ccomplex_t IT_1238 = 0.25*IT_1237;
    const ccomplex_t IT_1239 = 8*conj(IT_0122);
    const ccomplex_t IT_1240 = (-8)*conj(IT_0120);
    const ccomplex_t IT_1241 = IT_1239 + IT_1240;
    const ccomplex_t IT_1242 = IT_0579*IT_1241;
    const ccomplex_t IT_1243 = 0.25*IT_1242;
    const ccomplex_t IT_1244 = 0.25*IT_0723;
    const ccomplex_t IT_1245 = 8*conj(IT_0784);
    const ccomplex_t IT_1246 = (-8)*conj(IT_0567);
    const ccomplex_t IT_1247 = IT_1245 + IT_1246;
    const ccomplex_t IT_1248 = 0.25*IT_1176;
    const ccomplex_t IT_1249 = IT_0588*IT_1183 + IT_0617*IT_1211 + IT_1213 +
       IT_1216 + IT_1218 + IT_1222 + IT_1224 + IT_1226 + IT_1228 + IT_1231 +
       IT_1235 + IT_0649*IT_1236 + IT_1238 + IT_1243 + IT_0665*IT_1244 + IT_1247
      *IT_1248;
    const ccomplex_t IT_1250 = IT_0588*IT_1183 + IT_0617*IT_1211 + IT_1213 +
       IT_1216 + IT_1218 + IT_1222 + IT_1224 + IT_1226 + IT_1228 + IT_1231 +
       IT_1235 + IT_0649*IT_1236 + IT_1238 + IT_1243 + IT_0664*IT_1244 + IT_1246
      *IT_1248;
    const ccomplex_t IT_1251 = IT_0790*IT_1206;
    const ccomplex_t IT_1252 = (-2)*IT_0202;
    const ccomplex_t IT_1253 = IT_1150 + IT_1252;
    const ccomplex_t IT_1254 = 0.25*IT_1055;
    const ccomplex_t IT_1255 = IT_0469 + IT_0510;
    const ccomplex_t IT_1256 = 0.25*IT_1112;
    const ccomplex_t IT_1257 = 0.25*IT_1111;
    const ccomplex_t IT_1258 = 4*conj(IT_0489)*(IT_0778*IT_1115 + -IT_1120
      *IT_1123 + IT_1126*IT_1209 + IT_1000*IT_1210 + IT_1012*IT_1236 + IT_1251 +
       IT_1207*IT_1253 + IT_1094*IT_1254 + IT_1255*IT_1256 + IT_1161*IT_1257);
    const ccomplex_t IT_1259 = 4*IT_0489*(IT_0768*IT_1115 + IT_0793*IT_1206 +
       IT_1138*IT_1208 + IT_1139*IT_1209 + IT_1014*IT_1210 + IT_1207*IT_1220 +
       IT_1019*IT_1236 + IT_1107*IT_1254 + IT_1199*IT_1256 + IT_0608*IT_1257);
    const ccomplex_t IT_1260 = (-0.25)*IT_0840;
    const ccomplex_t IT_1261 = IT_0124*IT_1260;
    const ccomplex_t IT_1262 = (-0.25)*IT_0848;
    const ccomplex_t IT_1263 = IT_0242*IT_1262;
    const ccomplex_t IT_1264 = (-0.25)*IT_0867;
    const ccomplex_t IT_1265 = IT_0482*IT_1264;
    const ccomplex_t IT_1266 = (-0.25)*IT_0819;
    const ccomplex_t IT_1267 = IT_0821*IT_1266;
    const ccomplex_t IT_1268 = (-0.25)*IT_0828;
    const ccomplex_t IT_1269 = (-0.25)*IT_0968;
    const ccomplex_t IT_1270 = (-0.25)*IT_0856;
    const ccomplex_t IT_1271 = (-10)*s_13;
    const ccomplex_t IT_1272 = IT_0970 + IT_1271;
    const ccomplex_t IT_1273 = m_N_4*IT_1272;
    const ccomplex_t IT_1274 = (-0.25)*IT_1273;
    const ccomplex_t IT_1275 = (-0.25)*IT_0982;
    const ccomplex_t IT_1276 = (-0.25)*IT_0899;
    const ccomplex_t IT_1277 = IT_0902*IT_1276;
    const ccomplex_t IT_1278 = (-0.25)*IT_0953;
    const ccomplex_t IT_1279 = (-0.25)*IT_0914;
    const ccomplex_t IT_1280 = IT_0809*IT_0812 + (-0.25)*IT_0879*IT_0882 + (
      -0.25)*IT_0929*IT_0932 + (-0.25)*IT_0940*IT_0943 + (-0.25)*IT_0958*IT_0959
       + (-0.25)*IT_0976*IT_0977 + IT_1261 + IT_1263 + IT_1265 + IT_1267 +
       IT_0833*IT_1268 + IT_0969*IT_1269 + IT_0860*IT_1270 + IT_0975*IT_1274 +
       IT_0983*IT_1275 + IT_1277 + IT_0956*IT_1278 + IT_0917*IT_1279;
    const ccomplex_t IT_1281 = (-0.5)*IT_0871 + 0.5*IT_0872 + -IT_0906;
    const ccomplex_t IT_1282 = IT_0904 + IT_1281;
    const ccomplex_t IT_1283 = IT_0004*IT_1282;
    const ccomplex_t IT_1284 = 2*IT_1283;
    const ccomplex_t IT_1285 = (-0.5)*IT_1284;
    const ccomplex_t IT_1286 = IT_0869 + IT_1285;
    const ccomplex_t IT_1287 = s_14*s_34*IT_0000;
    const ccomplex_t IT_1288 = 0.5*IT_0922 + (-0.5)*IT_0923 + -IT_1287;
    const ccomplex_t IT_1289 = IT_0934 + IT_1288;
    const ccomplex_t IT_1290 = IT_0004*IT_1289;
    const ccomplex_t IT_1291 = (-2)*IT_1290;
    const ccomplex_t IT_1292 = (-0.5)*IT_1291;
    const ccomplex_t IT_1293 = IT_0919 + IT_1292;
    const ccomplex_t IT_1294 = IT_0001*IT_0823;
    const ccomplex_t IT_1295 = 0.666666666666667*IT_0609 + (-0.166666666666667
      )*IT_1294;
    const ccomplex_t IT_1296 = s_23 + IT_1295;
    const ccomplex_t IT_1297 = m_N_4*IT_1296;
    const ccomplex_t IT_1298 = (-6)*IT_1297;
    const ccomplex_t IT_1299 = (-0.25)*IT_1298;
    const ccomplex_t IT_1300 = (-0.666666666666667)*IT_0844 +
       0.166666666666667*IT_0964;
    const ccomplex_t IT_1301 = s_13 + IT_1300;
    const ccomplex_t IT_1302 = m_N_4*IT_1301;
    const ccomplex_t IT_1303 = 6*IT_1302;
    const ccomplex_t IT_1304 = (-0.25)*IT_1303;
    const ccomplex_t IT_1305 = (-0.666666666666667)*IT_0366 + (
      -0.166666666666667)*IT_0852;
    const ccomplex_t IT_1306 = s_23 + IT_1305;
    const ccomplex_t IT_1307 = m_N_4*IT_1306;
    const ccomplex_t IT_1308 = (-6)*IT_1307;
    const ccomplex_t IT_1309 = (-0.25)*IT_1308;
    const ccomplex_t IT_1310 = IT_0004*IT_0960;
    const ccomplex_t IT_1311 = (-2)*IT_0358 + (-0.5)*IT_1310;
    const ccomplex_t IT_1312 = IT_0357 + IT_1311;
    const ccomplex_t IT_1313 = IT_0822*IT_1312;
    const ccomplex_t IT_1314 = (-2)*IT_1313;
    const ccomplex_t IT_1315 = (-0.25)*IT_1314;
    const ccomplex_t IT_1316 = (-10)*s_23;
    const ccomplex_t IT_1317 = IT_1294 + IT_1316;
    const ccomplex_t IT_1318 = m_N_4*IT_1317;
    const ccomplex_t IT_1319 = (-0.25)*IT_1318;
    const ccomplex_t IT_1320 = -IT_0944;
    const ccomplex_t IT_1321 = IT_0401 + IT_0592 + IT_0946 + IT_0951 + IT_1320;
    const ccomplex_t IT_1322 = (-0.25)*IT_1321;
    const ccomplex_t IT_1323 = 0.5*IT_0909;
    const ccomplex_t IT_1324 = IT_0278 + IT_0290 + IT_0291 + IT_0868 + IT_1323;
    const ccomplex_t IT_1325 = (-0.25)*IT_1324;
    const ccomplex_t IT_1326 = (-0.5)*IT_0282;
    const ccomplex_t IT_1327 = (-0.25)*IT_0938;
    const ccomplex_t IT_1328 = (-0.5)*IT_0281;
    const ccomplex_t IT_1329 = 0.5*IT_0298;
    const ccomplex_t IT_1330 = (-0.5)*IT_0300;
    const ccomplex_t IT_1331 = IT_0279 + IT_1326 + IT_1327 + IT_1328 + IT_1329
       + IT_1330;
    const ccomplex_t IT_1332 = (-0.25)*IT_1331;
    const ccomplex_t IT_1333 = IT_1261 + IT_1263 + IT_1265 + IT_1267 + IT_1277
       + IT_0229*IT_1286 + -IT_0261*IT_1293 + IT_0832*IT_1299 + IT_0969*IT_1304 
      + IT_0860*IT_1309 + IT_0975*IT_1315 + IT_0831*IT_1319 + IT_0956*IT_1322 +
       IT_0993*IT_1325 + IT_0917*IT_1332;
    const ccomplex_t IT_1334 = IT_0004*IT_1024;
    const ccomplex_t IT_1335 = -IT_1334;
    const ccomplex_t IT_1336 = IT_0804 + IT_1335;
    const ccomplex_t IT_1337 = (-2)*conj(IT_0784);
    const ccomplex_t IT_1338 = (-2)*conj(IT_0787);
    const ccomplex_t IT_1339 = IT_1016 + IT_1337 + IT_1338;
    const ccomplex_t IT_1340 = (-2)*conj(IT_0565);
    const ccomplex_t IT_1341 = 2*conj(IT_0657);
    const ccomplex_t IT_1342 = 2*conj(IT_0788);
    const ccomplex_t IT_1343 = IT_1341 + IT_1342;
    const ccomplex_t IT_1344 = 4*conj(IT_0476);
    const ccomplex_t IT_1345 = 4*conj(IT_0489);
    const ccomplex_t IT_1346 = IT_1344 + IT_1345;
    const ccomplex_t IT_1347 = (-4)*conj(IT_0199);
    const ccomplex_t IT_1348 = (-4)*conj(IT_0202);
    const ccomplex_t IT_1349 = IT_1347 + IT_1348;
    const ccomplex_t IT_1350 = (-8)*conj(IT_0038);
    const ccomplex_t IT_1351 = (-8)*conj(IT_0081);
    const ccomplex_t IT_1352 = IT_1350 + IT_1351;
    const ccomplex_t IT_1353 = IT_0127*IT_1260 + IT_0245*IT_1262 + IT_0486
      *IT_1264 + IT_0617*IT_1266 + IT_0225*IT_1286 + -IT_0260*IT_1293 + IT_1015
      *IT_1299 + IT_1017*IT_1304 + -IT_0591*IT_1324 + conj(IT_0382)*IT_1336 +
       IT_1309*IT_1339 + IT_1315*IT_1340 + IT_1319*IT_1343 + IT_1276*IT_1346 +
       IT_1322*IT_1349 + IT_1332*IT_1352;
    const ccomplex_t IT_1354 = (-2)*IT_0382;
    const ccomplex_t IT_1355 = 4*IT_0164;
    const ccomplex_t IT_1356 = 0.25*IT_1355;
    const ccomplex_t IT_1357 = IT_0187*IT_1260;
    const ccomplex_t IT_1358 = IT_1093*IT_1264;
    const ccomplex_t IT_1359 = IT_1164*IT_1266;
    const ccomplex_t IT_1360 = 4*IT_0494;
    const ccomplex_t IT_1361 = 4*IT_0309;
    const ccomplex_t IT_1362 = IT_1360 + IT_1361;
    const ccomplex_t IT_1363 = (-4)*IT_0205;
    const ccomplex_t IT_1364 = (-4)*IT_0206;
    const ccomplex_t IT_1365 = IT_1363 + IT_1364;
    const ccomplex_t IT_1366 = 4*IT_0165;
    const ccomplex_t IT_1367 = IT_1158 + IT_1366;
    const ccomplex_t IT_1368 = (-8)*IT_0163;
    const ccomplex_t IT_1369 = (-8)*IT_0165;
    const ccomplex_t IT_1370 = IT_1368 + IT_1369;
    const ccomplex_t IT_1371 = (-4)*IT_0164;
    const ccomplex_t IT_1372 = (-0.25)*IT_1371;
    const ccomplex_t IT_1373 = (-4)*IT_0228*IT_1262 + (-4)*IT_0833*IT_1304 + (
      -4)*IT_0975*IT_1309 + (-4)*IT_0860*IT_1315 + (-4)*IT_0969*IT_1319 +
       IT_0976*IT_1354 + 4*IT_1293*(IT_0162 + IT_1356) + (-4)*IT_1357 + (-4)
      *IT_1358 + (-4)*IT_1359 + (-4)*IT_1276*IT_1362 + (-4)*IT_1322*IT_1365 + (
      -4)*IT_1325*IT_1367 + (-4)*IT_1332*IT_1370 + (-4)*IT_1286*(IT_0239 +
       IT_1372);
    const ccomplex_t IT_1374 = -IT_0879;
    const ccomplex_t IT_1375 = (-0.25)*IT_0974;
    const ccomplex_t IT_1376 = IT_0216*IT_0940 + IT_0228*IT_1262 + IT_0969
      *IT_1268 + IT_0833*IT_1269 + IT_0975*IT_1270 + IT_0857*IT_1274 + IT_1275
      *IT_1354 + IT_1357 + IT_1358 + IT_1359 + IT_1276*IT_1362 + IT_1278*IT_1365
       + IT_1279*IT_1370 + IT_0929*(IT_0162 + IT_1372) + (IT_0239 + IT_1356)
      *IT_1374 + (IT_0858 + IT_0859)*IT_1375;
    const ccomplex_t IT_1377 = IT_1015 + IT_1343;
    const ccomplex_t IT_1378 = 4*conj(IT_0494);
    const ccomplex_t IT_1379 = 4*conj(IT_0309);
    const ccomplex_t IT_1380 = IT_1378 + IT_1379;
    const ccomplex_t IT_1381 = (-8)*conj(IT_0163);
    const ccomplex_t IT_1382 = (-8)*conj(IT_0165);
    const ccomplex_t IT_1383 = IT_1381 + IT_1382;
    const ccomplex_t IT_1384 = (-4)*conj(IT_0205);
    const ccomplex_t IT_1385 = (-4)*conj(IT_0206);
    const ccomplex_t IT_1386 = IT_1384 + IT_1385;
    const ccomplex_t IT_1387 = (-2)*conj(IT_0382);
    const ccomplex_t IT_1388 = 4*conj(IT_0164);
    const ccomplex_t IT_1389 = 0.25*IT_1388;
    const ccomplex_t IT_1390 = (-4)*conj(IT_0164);
    const ccomplex_t IT_1391 = (-0.25)*IT_1390;
    const ccomplex_t IT_1392 = (-0.125)*IT_0649*IT_0819 + (-0.125)*IT_0190
      *IT_0840 + (-0.125)*IT_0232*IT_0848 + (-0.125)*IT_0637*IT_0867 + conj
      (IT_0376)*IT_0958 + (-0.125)*IT_1017*IT_1318 + (-0.5)*IT_1183*IT_1324 +
       0.5*conj(IT_0378)*IT_1336 + (-0.125)*IT_1314*IT_1339 + (-0.125)*IT_1308
      *IT_1340 + (-0.125)*IT_1303*IT_1377 + (-0.125)*IT_0899*IT_1380 + (-0.125)
      *IT_1331*IT_1383 + (-0.125)*IT_1321*IT_1386 + (-0.125)*IT_0976*IT_1387 + (
      -0.5)*IT_1293*(conj(IT_0162) + IT_1389) + 0.5*IT_1286*(conj(IT_0239) +
       IT_1391);
    const ccomplex_t IT_1393 = conj(IT_0376)*IT_0809 + IT_0215*IT_0940 + 2
      *conj(IT_0378)*IT_0958 + IT_0190*IT_1260 + IT_0232*IT_1262 + IT_0637
      *IT_1264 + IT_0649*IT_1266 + IT_1017*IT_1268 + IT_1016*IT_1274 + IT_1270
      *IT_1340 + (IT_1337 + IT_1338)*IT_1375 + IT_1269*IT_1377 + IT_1276*IT_1380
       + IT_1279*IT_1383 + IT_1278*IT_1386 + IT_1275*IT_1387 + IT_1374*(conj
      (IT_0239) + IT_1389) + IT_0929*(conj(IT_0162) + IT_1391);
    const ccomplex_t IT_1394 = (-2)*IT_0033;
    const ccomplex_t IT_1395 = IT_1022 + IT_1394;
    const ccomplex_t IT_1396 = 6*IT_0164;
    const ccomplex_t IT_1397 = 0.166666666666667*IT_1396;
    const ccomplex_t IT_1398 = (-6)*IT_0164;
    const ccomplex_t IT_1399 = (-0.166666666666667)*IT_1398;
    const ccomplex_t IT_1400 = IT_0673*IT_0857 + IT_0832*IT_1062 + IT_0412
      *IT_1093 + IT_0477*IT_1111 + IT_1000*IT_1148 + IT_0269*IT_1162 + IT_1010
      *IT_1169 + IT_0449*IT_1255 + IT_1021*IT_1395 + 6*IT_1198*(IT_0162 +
       IT_1397) + (-6)*IT_1192*(IT_0239 + IT_1399);
    const ccomplex_t IT_1401 = IT_0336*IT_0449 + IT_0352*IT_0464 + IT_0599
      *IT_1012 + IT_0423*IT_1093 + IT_0607*IT_1094 + IT_0477*IT_1112 + IT_0969
      *IT_1157 + IT_0267*IT_1162 + IT_1090*IT_1167 + IT_1008*IT_1169 + 6*IT_1188
      *(IT_0239 + IT_1397) + IT_1193*(IT_0162 + IT_1399);
    const ccomplex_t IT_1402 = (-12)*IT_0081;
    const ccomplex_t IT_1403 = (-12)*IT_0038 + IT_1402;
    const ccomplex_t IT_1404 = IT_0986*IT_1403;
    const ccomplex_t IT_1405 = IT_0267*IT_0469;
    const ccomplex_t IT_1406 = IT_1010*IT_1022;
    const ccomplex_t IT_1407 = (-4)*IT_0567;
    const ccomplex_t IT_1408 = IT_0588*IT_1407;
    const ccomplex_t IT_1409 = IT_1008*IT_1091;
    const ccomplex_t IT_1410 = IT_0269*IT_0336;
    const ccomplex_t IT_1411 = IT_0702*IT_1015;
    const ccomplex_t IT_1412 = (-12)*conj(IT_0038);
    const ccomplex_t IT_1413 = IT_0986*IT_1412;
    const ccomplex_t IT_1414 = IT_0267*IT_0450;
    const ccomplex_t IT_1415 = IT_1010*IT_1067;
    const ccomplex_t IT_1416 = IT_0999*IT_1107;
    const ccomplex_t IT_1417 = (-4)*conj(IT_0567);
    const ccomplex_t IT_1418 = IT_0588*IT_1417;
    const ccomplex_t IT_1419 = IT_0684*IT_1008;
    const ccomplex_t IT_1420 = IT_0269*IT_0338;
    const ccomplex_t IT_1421 = IT_0692*IT_1016;
    const ccomplex_t IT_1422 = IT_0635*IT_1012;
    const ccomplex_t IT_1423 = (-2)*IT_0309;
    const ccomplex_t IT_1424 = IT_0751*IT_1097;
    const ccomplex_t IT_1425 = IT_1015*IT_1099;
    const ccomplex_t IT_1426 = IT_0635*IT_1019;
    const ccomplex_t IT_1427 = (-2)*conj(IT_0309);
    const ccomplex_t IT_1428 = 6*IT_0078;
    const ccomplex_t IT_1429 = (-6)*IT_0224;
    const ccomplex_t IT_1430 = IT_0607*IT_1012;
    const ccomplex_t IT_1431 = IT_0269*IT_0270 + IT_0412*IT_0482 + IT_1010
      *IT_1056 + IT_0969*IT_1062 + IT_0483*IT_1111 + IT_0449*IT_1121 + IT_1094
      *IT_1148 + IT_1021*IT_1165 + IT_1198*IT_1428 + IT_1192*IT_1429 + IT_1430;
    const ccomplex_t IT_1432 = (-6)*conj(IT_0078);
    const ccomplex_t IT_1433 = 6*conj(IT_0224);
    const ccomplex_t IT_1434 = IT_1188*IT_1433;
    const ccomplex_t IT_1435 = IT_1015*IT_1157;
    const ccomplex_t IT_1436 = IT_0682*IT_1016;
    const ccomplex_t IT_1437 = IT_0267*IT_0272 + IT_0276*IT_0449 + IT_0423
      *IT_0486 + IT_0464*IT_0626 + IT_0693*IT_1008 + IT_0607*IT_1014 + IT_0675
      *IT_1090 + IT_0487*IT_1112 + IT_1192*IT_1432 + IT_1434 + IT_1435 + IT_1436;
    const ccomplex_t IT_1438 = (-6)*conj(IT_0224);
    const ccomplex_t IT_1439 = IT_0607*IT_1019;
    const ccomplex_t IT_1440 = 6*conj(IT_0078);
    const ccomplex_t IT_1441 = IT_1198*IT_1440;
    const ccomplex_t IT_1442 = IT_1107*IT_1148;
    const ccomplex_t IT_1443 = IT_0270 + IT_0352;
    const ccomplex_t IT_1444 = 8*IT_0784;
    const ccomplex_t IT_1445 = (-8)*IT_0567;
    const ccomplex_t IT_1446 = IT_1444 + IT_1445;
    const ccomplex_t IT_1447 = 8*IT_0122;
    const ccomplex_t IT_1448 = (-8)*IT_0120;
    const ccomplex_t IT_1449 = IT_1447 + IT_1448;
    const ccomplex_t IT_1450 = IT_0702*IT_1162;
    const ccomplex_t IT_1451 = IT_0999*IT_1169;
    const ccomplex_t IT_1452 = IT_0943*IT_0991;
    const ccomplex_t IT_1453 = IT_1047*IT_1253;
    const ccomplex_t IT_1454 = IT_0607*IT_1167;
    const ccomplex_t IT_1455 = IT_1148*IT_1165;
    const ccomplex_t IT_1456 = (-2)*IT_0489;
    const ccomplex_t IT_1457 = IT_1106 + IT_1456;
    const ccomplex_t IT_1458 = IT_1055*IT_1457;
    const ccomplex_t IT_1459 = (-2)*IT_0494;
    const ccomplex_t IT_1460 = IT_1423 + IT_1459;
    const ccomplex_t IT_1461 = IT_1133*IT_1460;
    const ccomplex_t IT_1462 = 16*IT_0568*IT_0728 + IT_0482*IT_1099 + IT_1121
      *IT_1157 + IT_1062*IT_1161 + IT_1029*(IT_1170 + IT_1172) + IT_1007*IT_1443
       + IT_0723*IT_1446 + IT_1035*IT_1449 + IT_1450 + IT_1451 + IT_1452 +
       IT_1453 + IT_1454 + IT_1455 + IT_1458 + IT_1461;
    const ccomplex_t IT_1463 = 8*IT_0787;
    const ccomplex_t IT_1464 = 16*IT_0568*IT_0728 + IT_0482*IT_1099 + IT_1121
      *IT_1157 + IT_1062*IT_1161 + IT_1029*(IT_1170 + IT_1171 + IT_1172) +
       IT_0656*IT_1177 + IT_1007*IT_1443 + IT_1035*IT_1449 + IT_1450 + IT_1451 +
       IT_1452 + IT_1453 + IT_1454 + IT_1455 + IT_1458 + IT_1461 + IT_0723*
      (IT_1444 + IT_1445 + IT_1463);
    const ccomplex_t IT_1465 = IT_0565*IT_0568;
    const ccomplex_t IT_1466 = (-16)*IT_1465;
    const ccomplex_t IT_1467 = IT_0702*IT_1168;
    const ccomplex_t IT_1468 = IT_0588*IT_1367;
    const ccomplex_t IT_1469 = IT_0692*IT_1169;
    const ccomplex_t IT_1470 = IT_0707*IT_1253;
    const ccomplex_t IT_1471 = IT_0607*IT_1122;
    const ccomplex_t IT_1472 = IT_0599*IT_1255;
    const ccomplex_t IT_1473 = IT_0682*IT_1202;
    const ccomplex_t IT_1474 = IT_0673*IT_1395;
    const ccomplex_t IT_1475 = IT_0615*IT_1164;
    const ccomplex_t IT_1476 = IT_0656*IT_1173;
    const ccomplex_t IT_1477 = IT_0579*IT_1449;
    const ccomplex_t IT_1478 = IT_0482*IT_0635 + IT_0647*IT_0821 + IT_0723
      *IT_1172 + IT_0625*IT_1443 + IT_1176*IT_1446 + IT_1466 + IT_1467 + IT_1468
       + IT_1469 + IT_1470 + IT_1471 + IT_1472 + IT_1473 + IT_1474 + IT_1475 +
       IT_1476 + IT_1477;
    const ccomplex_t IT_1479 = IT_0482*IT_0635 + IT_0647*IT_0821 + 8*IT_0723*
      (IT_0788 + 0.125*IT_1170 + 0.125*IT_1172) + IT_0625*IT_1443 + IT_1176*
      (IT_1444 + IT_1445 + IT_1463) + IT_1466 + IT_1467 + IT_1468 + IT_1469 +
       IT_1470 + IT_1471 + IT_1472 + IT_1473 + IT_1474 + IT_1475 + IT_1476 +
       IT_1477;
    const ccomplex_t IT_1480 = (-2)*conj(IT_0376);
    const ccomplex_t IT_1481 = IT_0976*IT_1480;
    const ccomplex_t IT_1482 = IT_0127*IT_0840;
    const ccomplex_t IT_1483 = IT_0245*IT_0848;
    const ccomplex_t IT_1484 = IT_0486*IT_0867;
    const ccomplex_t IT_1485 = IT_0617*IT_0819;
    const ccomplex_t IT_1486 = IT_0828*IT_1377;
    const ccomplex_t IT_1487 = IT_0968*IT_1017;
    const ccomplex_t IT_1488 = IT_0856*IT_1339;
    const ccomplex_t IT_1489 = (-2)*conj(IT_0378);
    const ccomplex_t IT_1490 = IT_0982*IT_1489;
    const ccomplex_t IT_1491 = IT_0899*IT_1346;
    const ccomplex_t IT_1492 = IT_0953*IT_1349;
    const ccomplex_t IT_1493 = (-8)*conj(IT_0382);
    const ccomplex_t IT_1494 = IT_0958*IT_1493;
    const ccomplex_t IT_1495 = IT_0914*IT_1352;
    const ccomplex_t IT_1496 = 4*conj(IT_0220);
    const ccomplex_t IT_1497 = 4*conj(IT_0224);
    const ccomplex_t IT_1498 = IT_1496 + IT_1497;
    const ccomplex_t IT_1499 = IT_0879*IT_1498;
    const ccomplex_t IT_1500 = (-4)*conj(IT_0033);
    const ccomplex_t IT_1501 = (-4)*conj(IT_0078);
    const ccomplex_t IT_1502 = IT_1500 + IT_1501;
    const ccomplex_t IT_1503 = IT_0929*IT_1502;
    const ccomplex_t IT_1504 = (-4)*conj(IT_0038);
    const ccomplex_t IT_1505 = (-4)*conj(IT_0081);
    const ccomplex_t IT_1506 = IT_1504 + IT_1505;
    const ccomplex_t IT_1507 = IT_0940*IT_1506;
    const ccomplex_t IT_1508 = 16*conj(IT_0728);
    const ccomplex_t IT_1509 = IT_0568*IT_1508;
    const ccomplex_t IT_1510 = IT_1007*IT_1214;
    const ccomplex_t IT_1511 = IT_0626*IT_0702;
    const ccomplex_t IT_1512 = IT_0705*IT_0999;
    const ccomplex_t IT_1513 = IT_0991*IT_1506;
    const ccomplex_t IT_1514 = IT_1047*IT_1220;
    const ccomplex_t IT_1515 = IT_0486*IT_1099;
    const ccomplex_t IT_1516 = IT_0607*IT_0685;
    const ccomplex_t IT_1517 = IT_0600*IT_1157;
    const ccomplex_t IT_1518 = IT_0608*IT_1062;
    const ccomplex_t IT_1519 = IT_0676*IT_1148;
    const ccomplex_t IT_1520 = (-2)*conj(IT_0489);
    const ccomplex_t IT_1521 = IT_1108 + IT_1520;
    const ccomplex_t IT_1522 = IT_1055*IT_1521;
    const ccomplex_t IT_1523 = (-2)*conj(IT_0494);
    const ccomplex_t IT_1524 = IT_1427 + IT_1523;
    const ccomplex_t IT_1525 = IT_1133*IT_1524;
    const ccomplex_t IT_1526 = IT_1035*IT_1241;
    const ccomplex_t IT_1527 = IT_0723*IT_1247;
    const ccomplex_t IT_1528 = IT_0191 + IT_0196*(IT_0124*IT_0203 + IT_0127
      *IT_0204 + IT_0187*IT_0207 + IT_0190*IT_0208) + IT_0214*(IT_0124*IT_0215 +
       IT_0127*IT_0216 + IT_0225*IT_0228 + IT_0229*IT_0232 + IT_0187*IT_0233 +
       IT_0190*IT_0234 + (conj(IT_0164) + conj(IT_0239))*IT_0242 + (IT_0164 +
       IT_0239)*IT_0245) + IT_0256 + IT_0259*(IT_0127*(IT_0162 + IT_0164) +
       IT_0124*(conj(IT_0162) + conj(IT_0164)) + IT_0228*IT_0233 + IT_0232
      *IT_0234 + IT_0215*IT_0242 + IT_0216*IT_0245 + IT_0187*IT_0260 + IT_0190
      *IT_0261) + (IT_0207*IT_0228 + IT_0208*IT_0232 + IT_0203*IT_0242 + IT_0204
      *IT_0245)*IT_0265 + (IT_0267 + IT_0269)*(IT_0271 + IT_0273 + IT_0275 +
       IT_0277) + (IT_0278 + IT_0279 + IT_0280 + IT_0281 + IT_0282 + IT_0283 +
       IT_0284 + IT_0288 + IT_0290 + IT_0292 + IT_0294 + IT_0297 + IT_0299 +
       IT_0301 + IT_0303 + IT_0305)*(IT_0310 + IT_0311) + (IT_0321 + IT_0331)*
      (IT_0332 + IT_0333 + IT_0334 + IT_0335 + IT_0337 + IT_0339 + IT_0341 +
       IT_0342 + IT_0344 + IT_0345 + IT_0346 + IT_0347 + IT_0349 + IT_0351 +
       IT_0353 + IT_0355) + IT_0399 + conj(IT_0199)*IT_0407*IT_0408 + IT_0412*
      (conj(IT_0078)*IT_0408 + IT_0415 + IT_0078*IT_0416 + IT_0418 + IT_0420) + 
      (conj(IT_0224)*IT_0408 + IT_0415 + IT_0224*IT_0416 + IT_0418 + IT_0420)
      *IT_0423 + (IT_0425 + IT_0427)*(IT_0428 + IT_0429 + IT_0430 + IT_0431 +
       IT_0432 + IT_0433 + IT_0434 + IT_0435 + IT_0436 + IT_0437 + IT_0438 +
       IT_0439 + IT_0440 + IT_0441 + IT_0442 + IT_0443) + IT_0078*IT_0449
      *IT_0450 + IT_0464*(IT_0342 + IT_0345 + IT_0346 + IT_0347 + IT_0465 +
       IT_0466 + IT_0467 + IT_0468 + IT_0470 + IT_0471) + IT_0473*IT_0506 +
       IT_0507*IT_0512 + (IT_0514 + IT_0516)*(IT_0517 + IT_0518 + IT_0519 +
       IT_0520 + IT_0521 + IT_0522 + IT_0523 + IT_0524 + IT_0525 + IT_0526 +
       IT_0527 + IT_0528 + IT_0529 + IT_0530 + IT_0531 + IT_0532) + conj(IT_0038
      )*IT_0534*IT_0535 + IT_0544*(IT_0546 + IT_0548 + IT_0550 + IT_0552 +
       IT_0554 + IT_0556) + IT_0565*IT_0730 + IT_0761 + 4*IT_0765*(IT_0122
      *IT_0768 + IT_0185*IT_0771 + 0.25*conj(IT_0120)*IT_0772 + 0.25*conj
      (IT_0183)*IT_0773 + 0.25*conj(IT_0122)*IT_0774 + 0.25*conj(IT_0185)
      *IT_0775) + 4*IT_0776*IT_0795 + (-4)*conj(IT_0219)*IT_0984 + conj(IT_0081)
      *IT_1013 + IT_0081*(IT_0269*IT_0493 + IT_0267*IT_0509 + IT_0534*IT_0589 +
       IT_0674*IT_1010 + IT_0999*IT_1014 + IT_1007*IT_1015 + IT_0702*(IT_1016 +
       IT_1017) + IT_1008*IT_1018 + IT_0625*IT_1019) + conj(IT_0078)*(IT_0449
      *IT_0469 + IT_1021*IT_1022) + 8*conj(IT_0663)*IT_1063 + IT_0663*(IT_0450
      *IT_0607 + IT_0589*IT_0991 + IT_0693*IT_1007 + IT_0416*IT_1047 + IT_0479
      *IT_1055 + IT_1035*IT_1066 + IT_1062*IT_1067) + IT_1021*IT_1069 + (IT_1008
       + IT_1010)*(IT_1070 + IT_1071 + IT_1072 + IT_1073) + (IT_1075 + IT_1077)*
      (IT_1078 + IT_1079 + IT_1080 + IT_1082 + IT_1083 + IT_1084 + IT_1086 +
       IT_1088) + IT_1090*(IT_1069 + conj(IT_0224)*IT_1091) + conj(IT_0202)*
      (IT_0336*IT_0412 + IT_0423*IT_0469 + IT_0707*IT_0857 + IT_0407*IT_1093 +
       IT_1047*IT_1094 + (IT_0753 + IT_0757)*IT_1097 + IT_0969*IT_1099 + IT_1105
      *IT_1106) + IT_0567*(IT_0338*IT_0673 + IT_0599*IT_0684 + IT_0625*IT_0693 +
       IT_0272*IT_0702 + IT_0416*IT_0707 + IT_0579*IT_1066 + IT_0615*IT_1108) +
       IT_0202*(IT_0338*IT_0412 + IT_0407*IT_0416 + IT_0423*IT_0450 + IT_0707
      *IT_1016 + (IT_0755 + IT_0759)*IT_1097 + IT_1017*IT_1099 + IT_1047*IT_1107
       + IT_1105*IT_1108) + IT_0476*(IT_0793*IT_1110 + IT_0338*IT_1111 + IT_0450
      *IT_1112 + IT_1105*IT_1113 + IT_0768*IT_1116) + (-4)*conj(IT_0494)*IT_1136
       + IT_1140 + conj(IT_0728)*IT_1160 + conj(IT_0565)*IT_1178 + IT_1179
      *IT_1180 + IT_0728*IT_1184 + IT_0220*IT_1194 + IT_0033*IT_1201 + conj
      (IT_0239)*IT_1205 + 4*conj(IT_0309)*(IT_1115*IT_1135 + IT_0780*IT_1206 +
       IT_1127*IT_1207 + IT_1125*IT_1208 + IT_1118*IT_1209 + IT_1094*IT_1210 +
       IT_1012*IT_1211) + 4*IT_0309*(IT_0771*IT_1115 + IT_0783*IT_1206 + IT_0708
      *IT_1207 + IT_1138*IT_1209 + IT_1107*IT_1210 + IT_1019*IT_1211) + 4
      *IT_0787*IT_1249 + 4*IT_0784*IT_1250 + 4*conj(IT_0476)*(IT_0778*IT_1115 +
       IT_1150*IT_1207 + IT_1118*IT_1208 + IT_1251 + IT_0469*IT_1256 + IT_0336
      *IT_1257) + IT_1258 + IT_1259 + (-4)*conj(IT_0223)*IT_1280 + (-4)*conj
      (IT_0382)*IT_1333 + (-4)*IT_0382*IT_1353 + conj(IT_0378)*IT_1373 + (-4)
      *conj(IT_0376)*IT_1376 + (-8)*IT_0378*IT_1392 + (-4)*IT_0376*IT_1393 +
       conj(IT_0033)*IT_1400 + conj(IT_0220)*IT_1401 + conj(IT_0163)*(IT_0544
      *IT_0545 + IT_0702*IT_0832 + IT_0692*IT_0857 + IT_0999*IT_1094 + IT_0534
      *IT_1158 + IT_1404 + IT_1405 + IT_1406 + IT_1408 + IT_1409 + IT_1410) +
       conj(IT_0165)*(IT_0702*IT_0832 + IT_0692*IT_0857 + IT_0969*IT_1007 +
       IT_0999*IT_1094 + IT_0534*IT_1367 + IT_1404 + IT_1405 + IT_1406 + IT_1408
       + IT_1409 + IT_1410) + IT_0163*(IT_0544*IT_0547 + IT_1411 + IT_1413 +
       IT_1414 + IT_1415 + IT_1416 + IT_1418 + IT_1419 + IT_1420 + IT_1421) +
       IT_0165*(IT_1007*IT_1017 + IT_0534*IT_1181 + IT_1411 + IT_1413 + IT_1414 
      + IT_1415 + IT_1416 + IT_1418 + IT_1419 + IT_1420 + IT_1421) + conj
      (IT_0205)*(IT_0407*IT_0417 + IT_0749*IT_1097 + IT_0832*IT_1099 + IT_1422) 
      + conj(IT_0206)*(IT_0407*IT_0482 + IT_1000*IT_1047 + IT_0749*IT_1097 +
       IT_0832*IT_1099 + IT_1422 + IT_1105*IT_1423) + IT_0205*(IT_1424 + IT_1425
       + IT_1426) + IT_0206*(IT_0407*IT_0419 + IT_1014*IT_1047 + IT_1424 +
       IT_1425 + IT_1426 + IT_1105*IT_1427) + conj(IT_0164)*(IT_0682*IT_0857 +
       IT_1094*IT_1148 + IT_0832*IT_1157 + IT_1188*IT_1203 + IT_0449*IT_1355 +
       IT_1198*IT_1428 + IT_1192*(IT_1204 + IT_1429) + IT_1430) + conj(IT_0162)
      *IT_1431 + IT_0239*IT_1437 + IT_0162*(IT_0269*IT_0272 + IT_0412*IT_0486 +
       IT_0449*IT_0600 + IT_0693*IT_1010 + IT_0675*IT_1021 + IT_1017*IT_1062 +
       IT_0487*IT_1111 + IT_1192*IT_1438 + IT_1439 + IT_1441 + IT_1442) +
       IT_0164*(IT_1434 + IT_1435 + IT_1436 + IT_1192*(IT_1432 + IT_1438) +
       IT_1439 + IT_1441 + IT_1442) + conj(IT_0657)*IT_1462 + conj(IT_0788)
      *IT_1464 + conj(IT_0784)*IT_1478 + conj(IT_0787)*IT_1479 + IT_0219*
      (IT_0974*IT_1340 + IT_1481 + IT_1482 + IT_1483 + IT_1484 + IT_1485 +
       IT_1486 + IT_1487 + IT_1488 + IT_1490 + IT_1491 + IT_1492 + IT_1494 +
       IT_1495 + IT_1499 + IT_1503 + IT_1507) + IT_0223*(IT_1273*IT_1340 +
       IT_1481 + IT_1482 + IT_1483 + IT_1484 + IT_1485 + IT_1486 + IT_1487 +
       IT_1488 + IT_1490 + IT_1491 + IT_1492 + IT_1494 + IT_1495 + IT_1499 +
       IT_1503 + IT_1507) + IT_0657*(IT_0664*IT_1029 + IT_1509 + IT_1510 +
       IT_1511 + IT_1512 + IT_1513 + IT_1514 + IT_1515 + IT_1516 + IT_1517 +
       IT_1518 + IT_1519 + IT_1522 + IT_1525 + IT_1526 + IT_1527) + IT_0788*
      (IT_0665*IT_1029 + IT_1509 + IT_1510 + IT_1511 + IT_1512 + IT_1513 +
       IT_1514 + IT_1515 + IT_1516 + IT_1517 + IT_1518 + IT_1519 + IT_1522 +
       IT_1525 + IT_1526 + IT_1527);
    return create_ccomplex_return(IT_1528);
}

