#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_2_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_2_to_anti_W_W(
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
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = sin(theta_W);
    const ccomplex_t IT_0001 = cpow(IT_0000, -1);
    const ccomplex_t IT_0002 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0005 = IT_0001*IT_0004;
    const ccomplex_t IT_0006 = 1.4142135623731*IT_0005;
    const ccomplex_t IT_0007 = (0 + _Complex_I*1)*(IT_0003 + 0.5*IT_0006);
    const ccomplex_t IT_0008 = -IT_0007;
    const ccomplex_t IT_0009 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0010 = IT_0001*IT_0009;
    const ccomplex_t IT_0011 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0012 = IT_0001*IT_0011;
    const ccomplex_t IT_0013 = 1.4142135623731*IT_0012;
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*(IT_0010 + 0.5*IT_0013);
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0008*IT_0015;
    const ccomplex_t IT_0017 = pow(m_W, 2);
    const ccomplex_t IT_0018 = pow(m_N_1, 2);
    const ccomplex_t IT_0019 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0020 = IT_0016*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0023 = IT_0001*IT_0022;
    const ccomplex_t IT_0024 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0025 = IT_0001*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0030 = IT_0001*IT_0029;
    const ccomplex_t IT_0031 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0032 = IT_0001*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + 0.5*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0028*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0017 + IT_0018 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*IT_0038;
    const ccomplex_t IT_0040 = (-2)*IT_0021 + (-2)*IT_0039;
    const ccomplex_t IT_0041 = pow(m_W, -4);
    const ccomplex_t IT_0042 = s_13*s_24*s_34*IT_0041;
    const ccomplex_t IT_0043 = pow(m_W, -2);
    const ccomplex_t IT_0044 = s_13*s_23*IT_0043;
    const ccomplex_t IT_0045 = s_14*s_24*IT_0043;
    const ccomplex_t IT_0046 = -IT_0044 + -IT_0045;
    const ccomplex_t IT_0047 = s_12 + IT_0042 + IT_0046;
    const ccomplex_t IT_0048 = m_N_1*IT_0047;
    const ccomplex_t IT_0049 = cos(theta_W);
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*IT_0001*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = cpow(IT_0049, -1);
    const ccomplex_t IT_0053 = IT_0000*IT_0052;
    const ccomplex_t IT_0054 = N_u1*conj(N_u2)*e_em;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = IT_0001*IT_0049;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = N_d1*conj(N_d2)*e_em;
    const ccomplex_t IT_0059 = IT_0053*IT_0058;
    const ccomplex_t IT_0060 = IT_0056*IT_0058;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0055 + IT_0057 + 
      -IT_0059 + -IT_0060);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0051*IT_0062;
    const ccomplex_t IT_0064 = m_N_2*IT_0063;
    const ccomplex_t IT_0065 = conj(N_d1)*N_d2*e_em;
    const ccomplex_t IT_0066 = IT_0053*IT_0065;
    const ccomplex_t IT_0067 = IT_0056*IT_0065;
    const ccomplex_t IT_0068 = conj(N_u1)*N_u2*e_em;
    const ccomplex_t IT_0069 = IT_0053*IT_0068;
    const ccomplex_t IT_0070 = IT_0056*IT_0068;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + IT_0067 + 
      -IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0051*IT_0072;
    const ccomplex_t IT_0074 = m_N_1*IT_0073;
    const ccomplex_t IT_0075 = pow(m_N_2, 2);
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0018 + -IT_0075 + -reg_prop, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (IT_0064 + -IT_0074)*IT_0077;
    const ccomplex_t IT_0079 = 2*conj(IT_0078);
    const ccomplex_t IT_0080 = cos(alpha);
    const ccomplex_t IT_0081 = IT_0052*IT_0080;
    const ccomplex_t IT_0082 = N_B2*N_u1*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = sin(alpha);
    const ccomplex_t IT_0085 = IT_0052*IT_0084;
    const ccomplex_t IT_0086 = N_B2*N_d1*e_em;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = N_B1*N_u2*e_em;
    const ccomplex_t IT_0089 = IT_0081*IT_0088;
    const ccomplex_t IT_0090 = N_B1*N_d2*e_em;
    const ccomplex_t IT_0091 = IT_0085*IT_0090;
    const ccomplex_t IT_0092 = IT_0083 + IT_0087 + IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = IT_0001*IT_0080;
    const ccomplex_t IT_0094 = N_W2*N_u1*e_em;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = IT_0001*IT_0084;
    const ccomplex_t IT_0097 = N_W2*N_d1*e_em;
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = N_W1*N_u2*e_em;
    const ccomplex_t IT_0100 = IT_0093*IT_0099;
    const ccomplex_t IT_0101 = N_W1*N_d2*e_em;
    const ccomplex_t IT_0102 = IT_0096*IT_0101;
    const ccomplex_t IT_0103 = -IT_0095 + -IT_0098 + -IT_0100 + -IT_0102;
    const ccomplex_t IT_0104 = IT_0092 + IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = sin(beta);
    const ccomplex_t IT_0108 = cos(beta);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*(IT_0080
      *IT_0107 + -IT_0084*IT_0108);
    const ccomplex_t IT_0110 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0018 + -IT_0075 + -reg_prop, -1);
    const ccomplex_t IT_0111 = IT_0109*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = IT_0106*IT_0112;
    const ccomplex_t IT_0114 = IT_0081*IT_0086;
    const ccomplex_t IT_0115 = IT_0094*IT_0096;
    const ccomplex_t IT_0116 = IT_0081*IT_0090;
    const ccomplex_t IT_0117 = IT_0096*IT_0099;
    const ccomplex_t IT_0118 = IT_0114 + IT_0115 + IT_0116 + IT_0117;
    const ccomplex_t IT_0119 = IT_0082*IT_0085;
    const ccomplex_t IT_0120 = IT_0093*IT_0097;
    const ccomplex_t IT_0121 = IT_0085*IT_0088;
    const ccomplex_t IT_0122 = IT_0093*IT_0101;
    const ccomplex_t IT_0123 = -IT_0119 + -IT_0120 + -IT_0121 + -IT_0122;
    const ccomplex_t IT_0124 = IT_0118 + IT_0123;
    const ccomplex_t IT_0125 = (0 + _Complex_I*1)*IT_0124;
    const ccomplex_t IT_0126 = 0.5*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_W*e_em*IT_0001*(IT_0084
      *IT_0107 + IT_0080*IT_0108);
    const ccomplex_t IT_0128 = cpow((-2)*s_12 + -IT_0018 + -IT_0075 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0129 = IT_0127*IT_0128;
    const ccomplex_t IT_0130 = (0 + _Complex_I*1)*IT_0129;
    const ccomplex_t IT_0131 = IT_0126*IT_0130;
    const ccomplex_t IT_0132 = -IT_0113 + -IT_0131;
    const ccomplex_t IT_0133 = (-2)*conj(IT_0132);
    const ccomplex_t IT_0134 = IT_0079 + IT_0133;
    const ccomplex_t IT_0135 = pow(s_13, 2);
    const ccomplex_t IT_0136 = IT_0043*IT_0135;
    const ccomplex_t IT_0137 = pow(s_14, 2);
    const ccomplex_t IT_0138 = IT_0043*IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = s_13*s_14*s_34*IT_0041;
    const ccomplex_t IT_0141 = -IT_0018 + -IT_0140;
    const ccomplex_t IT_0142 = IT_0139 + IT_0141;
    const ccomplex_t IT_0143 = m_N_2*IT_0142;
    const ccomplex_t IT_0144 = m_N_1*IT_0063;
    const ccomplex_t IT_0145 = m_N_2*IT_0073;
    const ccomplex_t IT_0146 = IT_0077*(IT_0144 + -IT_0145);
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = 2*conj(IT_0147);
    const ccomplex_t IT_0149 = conj(N_B1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0150 = IT_0081*IT_0149;
    const ccomplex_t IT_0151 = conj(N_B1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0152 = IT_0085*IT_0151;
    const ccomplex_t IT_0153 = conj(N_B2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0154 = IT_0081*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0156 = IT_0085*IT_0155;
    const ccomplex_t IT_0157 = IT_0150 + IT_0152 + IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = conj(N_W1)*conj(N_u2)*e_em;
    const ccomplex_t IT_0159 = IT_0093*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W1)*conj(N_d2)*e_em;
    const ccomplex_t IT_0161 = IT_0096*IT_0160;
    const ccomplex_t IT_0162 = conj(N_W2)*conj(N_u1)*e_em;
    const ccomplex_t IT_0163 = IT_0093*IT_0162;
    const ccomplex_t IT_0164 = conj(N_W2)*conj(N_d1)*e_em;
    const ccomplex_t IT_0165 = IT_0096*IT_0164;
    const ccomplex_t IT_0166 = -IT_0159 + -IT_0161 + -IT_0163 + -IT_0165;
    const ccomplex_t IT_0167 = IT_0157 + IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = (-0.5)*IT_0168;
    const ccomplex_t IT_0170 = IT_0112*IT_0169;
    const ccomplex_t IT_0171 = IT_0081*IT_0151;
    const ccomplex_t IT_0172 = IT_0096*IT_0158;
    const ccomplex_t IT_0173 = IT_0081*IT_0155;
    const ccomplex_t IT_0174 = IT_0096*IT_0162;
    const ccomplex_t IT_0175 = IT_0171 + IT_0172 + IT_0173 + IT_0174;
    const ccomplex_t IT_0176 = IT_0085*IT_0149;
    const ccomplex_t IT_0177 = IT_0093*IT_0160;
    const ccomplex_t IT_0178 = IT_0085*IT_0153;
    const ccomplex_t IT_0179 = IT_0093*IT_0164;
    const ccomplex_t IT_0180 = -IT_0176 + -IT_0177 + -IT_0178 + -IT_0179;
    const ccomplex_t IT_0181 = IT_0175 + IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = 0.5*IT_0182;
    const ccomplex_t IT_0184 = IT_0130*IT_0183;
    const ccomplex_t IT_0185 = -IT_0170 + -IT_0184;
    const ccomplex_t IT_0186 = (-2)*conj(IT_0185);
    const ccomplex_t IT_0187 = IT_0148 + IT_0186;
    const ccomplex_t IT_0188 = s_34*IT_0041;
    const ccomplex_t IT_0189 = s_14*s_23;
    const ccomplex_t IT_0190 = s_13*s_24;
    const ccomplex_t IT_0191 = s_12*s_34;
    const ccomplex_t IT_0192 = -IT_0190 + -IT_0191;
    const ccomplex_t IT_0193 = IT_0189 + IT_0192;
    const ccomplex_t IT_0194 = IT_0188*IT_0193;
    const ccomplex_t IT_0195 = (-0.25)*IT_0194;
    const ccomplex_t IT_0196 = s_12*IT_0017;
    const ccomplex_t IT_0197 = IT_0043*IT_0196;
    const ccomplex_t IT_0198 = (-0.5)*IT_0197;
    const ccomplex_t IT_0199 = s_12 + IT_0195 + IT_0198;
    const ccomplex_t IT_0200 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0201 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0202 = IT_0001*IT_0201;
    const ccomplex_t IT_0203 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0204 = IT_0001*IT_0203;
    const ccomplex_t IT_0205 = 1.4142135623731*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*(IT_0202 + (-0.5)*IT_0205);
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = IT_0035*IT_0207;
    const ccomplex_t IT_0209 = IT_0037*IT_0200*IT_0208;
    const ccomplex_t IT_0210 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0211 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0212 = IT_0001*IT_0211;
    const ccomplex_t IT_0213 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0214 = IT_0001*IT_0213;
    const ccomplex_t IT_0215 = 1.4142135623731*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*(IT_0212 + (-0.5)*IT_0215);
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0015*IT_0217;
    const ccomplex_t IT_0219 = IT_0019*IT_0210*IT_0218;
    const ccomplex_t IT_0220 = -IT_0209 + -IT_0219;
    const ccomplex_t IT_0221 = 8*IT_0078;
    const ccomplex_t IT_0222 = (-8)*IT_0132;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = 8*conj(IT_0078);
    const ccomplex_t IT_0225 = (-8)*conj(IT_0132);
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = IT_0021 + IT_0039;
    const ccomplex_t IT_0228 = m_N_1*IT_0227;
    const ccomplex_t IT_0229 = -IT_0228;
    const ccomplex_t IT_0230 = 8*IT_0132;
    const ccomplex_t IT_0231 = (-8)*IT_0078;
    const ccomplex_t IT_0232 = IT_0230 + IT_0231;
    const ccomplex_t IT_0233 = 8*conj(IT_0132);
    const ccomplex_t IT_0234 = (-8)*conj(IT_0078);
    const ccomplex_t IT_0235 = IT_0233 + IT_0234;
    const ccomplex_t IT_0236 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0237 = IT_0001*IT_0236;
    const ccomplex_t IT_0238 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0239 = IT_0001*IT_0238;
    const ccomplex_t IT_0240 = 1.4142135623731*IT_0239;
    const ccomplex_t IT_0241 = (0 + _Complex_I*1)*(IT_0237 + (-0.5)*IT_0240);
    const ccomplex_t IT_0242 = -IT_0241;
    const ccomplex_t IT_0243 = IT_0008*IT_0242;
    const ccomplex_t IT_0244 = IT_0019*IT_0210*IT_0243;
    const ccomplex_t IT_0245 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0246 = IT_0001*IT_0245;
    const ccomplex_t IT_0247 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0248 = IT_0001*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0246 + (-0.5)*IT_0249);
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = IT_0207*IT_0251;
    const ccomplex_t IT_0253 = IT_0037*IT_0252;
    const ccomplex_t IT_0254 = (0 + _Complex_I*1)*IT_0253;
    const ccomplex_t IT_0255 = m_N_1*IT_0254;
    const ccomplex_t IT_0256 = -IT_0244 + IT_0255;
    const ccomplex_t IT_0257 = IT_0217*IT_0242;
    const ccomplex_t IT_0258 = IT_0019*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*IT_0258;
    const ccomplex_t IT_0260 = m_N_1*IT_0259;
    const ccomplex_t IT_0261 = IT_0028*IT_0251;
    const ccomplex_t IT_0262 = IT_0037*IT_0200*IT_0261;
    const ccomplex_t IT_0263 = -IT_0262;
    const ccomplex_t IT_0264 = conj(IT_0256) + conj(IT_0260) + conj(IT_0263);
    const ccomplex_t IT_0265 = 8*IT_0147;
    const ccomplex_t IT_0266 = (-8)*IT_0185;
    const ccomplex_t IT_0267 = IT_0265 + IT_0266;
    const ccomplex_t IT_0268 = IT_0256 + IT_0260 + IT_0263;
    const ccomplex_t IT_0269 = 8*conj(IT_0147);
    const ccomplex_t IT_0270 = (-8)*conj(IT_0185);
    const ccomplex_t IT_0271 = IT_0269 + IT_0270;
    const ccomplex_t IT_0272 = s_14*s_34*IT_0043;
    const ccomplex_t IT_0273 = s_13*s_23;
    const ccomplex_t IT_0274 = (-2)*IT_0273;
    const ccomplex_t IT_0275 = IT_0196 + IT_0274;
    const ccomplex_t IT_0276 = IT_0043*IT_0275;
    const ccomplex_t IT_0277 = (s_13 + -IT_0272)*(s_12 + (-0.5)*IT_0276);
    const ccomplex_t IT_0278 = (-2)*IT_0277;
    const ccomplex_t IT_0279 = s_24*s_34*IT_0043;
    const ccomplex_t IT_0280 = (s_12 + (-0.5)*IT_0276)*(s_23 + -IT_0279);
    const ccomplex_t IT_0281 = (-2)*IT_0280;
    const ccomplex_t IT_0282 = IT_0050*IT_0062;
    const ccomplex_t IT_0283 = IT_0077*IT_0282;
    const ccomplex_t IT_0284 = IT_0063*IT_0077;
    const ccomplex_t IT_0285 = -IT_0284;
    const ccomplex_t IT_0286 = 2*IT_0285;
    const ccomplex_t IT_0287 = conj(IT_0283)*IT_0286;
    const ccomplex_t IT_0288 = 2*conj(IT_0285);
    const ccomplex_t IT_0289 = IT_0283*IT_0288;
    const ccomplex_t IT_0290 = IT_0050*IT_0072;
    const ccomplex_t IT_0291 = IT_0077*IT_0290;
    const ccomplex_t IT_0292 = IT_0073*IT_0077;
    const ccomplex_t IT_0293 = -IT_0292;
    const ccomplex_t IT_0294 = 2*IT_0293;
    const ccomplex_t IT_0295 = conj(IT_0291)*IT_0294;
    const ccomplex_t IT_0296 = 2*conj(IT_0293);
    const ccomplex_t IT_0297 = IT_0291*IT_0296;
    const ccomplex_t IT_0298 = s_23*IT_0018;
    const ccomplex_t IT_0299 = s_12*s_13;
    const ccomplex_t IT_0300 = s_23*IT_0043*IT_0135;
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = IT_0299 + IT_0301;
    const ccomplex_t IT_0303 = s_13*s_14*s_24*IT_0043;
    const ccomplex_t IT_0304 = -IT_0303;
    const ccomplex_t IT_0305 = s_23*IT_0043*IT_0137;
    const ccomplex_t IT_0306 = -IT_0305;
    const ccomplex_t IT_0307 = s_13*s_14*s_23*s_34*IT_0041;
    const ccomplex_t IT_0308 = 2*IT_0307;
    const ccomplex_t IT_0309 = IT_0298 + IT_0302 + IT_0304 + IT_0306 + IT_0308;
    const ccomplex_t IT_0310 = s_13*IT_0075;
    const ccomplex_t IT_0311 = s_12*s_23;
    const ccomplex_t IT_0312 = pow(s_23, 2);
    const ccomplex_t IT_0313 = s_13*IT_0043*IT_0312;
    const ccomplex_t IT_0314 = (-2)*IT_0313;
    const ccomplex_t IT_0315 = IT_0311 + IT_0314;
    const ccomplex_t IT_0316 = pow(s_24, 2);
    const ccomplex_t IT_0317 = s_13*IT_0043*IT_0316;
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = s_14*s_23*s_24*IT_0043;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = s_13*s_23*s_24*s_34*IT_0041;
    const ccomplex_t IT_0322 = 2*IT_0321;
    const ccomplex_t IT_0323 = IT_0310 + IT_0315 + IT_0318 + IT_0320 + IT_0322;
    const ccomplex_t IT_0324 = (-2)*IT_0050;
    const ccomplex_t IT_0325 = IT_0062*IT_0324;
    const ccomplex_t IT_0326 = IT_0077*IT_0325;
    const ccomplex_t IT_0327 = IT_0286*conj(IT_0326);
    const ccomplex_t IT_0328 = IT_0288*IT_0326;
    const ccomplex_t IT_0329 = IT_0072*IT_0324;
    const ccomplex_t IT_0330 = IT_0077*IT_0329;
    const ccomplex_t IT_0331 = IT_0294*conj(IT_0330);
    const ccomplex_t IT_0332 = IT_0296*IT_0330;
    const ccomplex_t IT_0333 = m_N_1*m_N_2;
    const ccomplex_t IT_0334 = s_13*IT_0017;
    const ccomplex_t IT_0335 = s_14*s_34;
    const ccomplex_t IT_0336 = -IT_0335;
    const ccomplex_t IT_0337 = IT_0334 + IT_0336;
    const ccomplex_t IT_0338 = IT_0043*IT_0337;
    const ccomplex_t IT_0339 = -IT_0272 + (-0.25)*IT_0338;
    const ccomplex_t IT_0340 = s_13 + IT_0339;
    const ccomplex_t IT_0341 = IT_0333*IT_0340;
    const ccomplex_t IT_0342 = (-4)*IT_0341;
    const ccomplex_t IT_0343 = s_23*IT_0017;
    const ccomplex_t IT_0344 = s_24*s_34;
    const ccomplex_t IT_0345 = -IT_0344;
    const ccomplex_t IT_0346 = IT_0343 + IT_0345;
    const ccomplex_t IT_0347 = IT_0043*IT_0346;
    const ccomplex_t IT_0348 = -IT_0279 + (-0.25)*IT_0347;
    const ccomplex_t IT_0349 = s_23 + IT_0348;
    const ccomplex_t IT_0350 = IT_0333*IT_0349;
    const ccomplex_t IT_0351 = (-4)*IT_0350;
    const ccomplex_t IT_0352 = cpow((-2)*s_23 + IT_0017 + IT_0075 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0353 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0354 = IT_0001*IT_0353;
    const ccomplex_t IT_0355 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0356 = IT_0001*IT_0355;
    const ccomplex_t IT_0357 = 1.4142135623731*IT_0356;
    const ccomplex_t IT_0358 = (0 + _Complex_I*1)*(IT_0354 + (-0.5)*IT_0357);
    const ccomplex_t IT_0359 = -IT_0358;
    const ccomplex_t IT_0360 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0361 = IT_0001*IT_0360;
    const ccomplex_t IT_0362 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0363 = IT_0001*IT_0362;
    const ccomplex_t IT_0364 = 1.4142135623731*IT_0363;
    const ccomplex_t IT_0365 = (0 + _Complex_I*1)*(IT_0361 + (-0.5)*IT_0364);
    const ccomplex_t IT_0366 = -IT_0365;
    const ccomplex_t IT_0367 = IT_0359*IT_0366;
    const ccomplex_t IT_0368 = IT_0352*IT_0367;
    const ccomplex_t IT_0369 = (0 + _Complex_I*1)*IT_0368;
    const ccomplex_t IT_0370 = cpow((-2)*s_23 + IT_0017 + IT_0075 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0371 = N_W2*e_em*conj(V_Wp1);
    const ccomplex_t IT_0372 = IT_0001*IT_0371;
    const ccomplex_t IT_0373 = N_u2*conj(V_u1)*e_em;
    const ccomplex_t IT_0374 = IT_0001*IT_0373;
    const ccomplex_t IT_0375 = 1.4142135623731*IT_0374;
    const ccomplex_t IT_0376 = (0 + _Complex_I*1)*(IT_0372 + (-0.5)*IT_0375);
    const ccomplex_t IT_0377 = -IT_0376;
    const ccomplex_t IT_0378 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0379 = IT_0001*IT_0378;
    const ccomplex_t IT_0380 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0381 = IT_0001*IT_0380;
    const ccomplex_t IT_0382 = 1.4142135623731*IT_0381;
    const ccomplex_t IT_0383 = (0 + _Complex_I*1)*(IT_0379 + (-0.5)*IT_0382);
    const ccomplex_t IT_0384 = -IT_0383;
    const ccomplex_t IT_0385 = IT_0377*IT_0384;
    const ccomplex_t IT_0386 = IT_0370*IT_0385;
    const ccomplex_t IT_0387 = (0 + _Complex_I*1)*IT_0386;
    const ccomplex_t IT_0388 = -IT_0369 + -IT_0387;
    const ccomplex_t IT_0389 = IT_0286*conj(IT_0388);
    const ccomplex_t IT_0390 = conj(IT_0039)*IT_0286;
    const ccomplex_t IT_0391 = conj(IT_0021)*IT_0286;
    const ccomplex_t IT_0392 = IT_0288*IT_0388;
    const ccomplex_t IT_0393 = IT_0039*IT_0288;
    const ccomplex_t IT_0394 = IT_0021*IT_0288;
    const ccomplex_t IT_0395 = IT_0254 + IT_0259;
    const ccomplex_t IT_0396 = IT_0294*conj(IT_0395);
    const ccomplex_t IT_0397 = conj(N_W2)*e_em*U_Wm1;
    const ccomplex_t IT_0398 = IT_0001*IT_0397;
    const ccomplex_t IT_0399 = conj(N_d2)*U_d1*e_em;
    const ccomplex_t IT_0400 = IT_0001*IT_0399;
    const ccomplex_t IT_0401 = 1.4142135623731*IT_0400;
    const ccomplex_t IT_0402 = (0 + _Complex_I*1)*(IT_0398 + 0.5*IT_0401);
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0405 = IT_0001*IT_0404;
    const ccomplex_t IT_0406 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0407 = IT_0001*IT_0406;
    const ccomplex_t IT_0408 = 1.4142135623731*IT_0407;
    const ccomplex_t IT_0409 = (0 + _Complex_I*1)*(IT_0405 + 0.5*IT_0408);
    const ccomplex_t IT_0410 = -IT_0409;
    const ccomplex_t IT_0411 = IT_0403*IT_0410;
    const ccomplex_t IT_0412 = IT_0370*IT_0411;
    const ccomplex_t IT_0413 = (0 + _Complex_I*1)*IT_0412;
    const ccomplex_t IT_0414 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0415 = IT_0001*IT_0414;
    const ccomplex_t IT_0416 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0417 = IT_0001*IT_0416;
    const ccomplex_t IT_0418 = 1.4142135623731*IT_0417;
    const ccomplex_t IT_0419 = (0 + _Complex_I*1)*(IT_0415 + 0.5*IT_0418);
    const ccomplex_t IT_0420 = -IT_0419;
    const ccomplex_t IT_0421 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0422 = IT_0001*IT_0421;
    const ccomplex_t IT_0423 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0424 = IT_0001*IT_0423;
    const ccomplex_t IT_0425 = 1.4142135623731*IT_0424;
    const ccomplex_t IT_0426 = (0 + _Complex_I*1)*(IT_0422 + 0.5*IT_0425);
    const ccomplex_t IT_0427 = -IT_0426;
    const ccomplex_t IT_0428 = IT_0420*IT_0427;
    const ccomplex_t IT_0429 = IT_0352*IT_0428;
    const ccomplex_t IT_0430 = (0 + _Complex_I*1)*IT_0429;
    const ccomplex_t IT_0431 = -IT_0413 + -IT_0430;
    const ccomplex_t IT_0432 = IT_0294*conj(IT_0431);
    const ccomplex_t IT_0433 = IT_0296*IT_0395;
    const ccomplex_t IT_0434 = IT_0296*IT_0431;
    const ccomplex_t IT_0435 = s_34*m_N_1*m_N_2*IT_0043;
    const ccomplex_t IT_0436 = s_13*s_34*IT_0043;
    const ccomplex_t IT_0437 = -IT_0436;
    const ccomplex_t IT_0438 = s_14 + IT_0437;
    const ccomplex_t IT_0439 = IT_0435*IT_0438;
    const ccomplex_t IT_0440 = s_23*s_34*IT_0043;
    const ccomplex_t IT_0441 = -IT_0440;
    const ccomplex_t IT_0442 = s_24 + IT_0441;
    const ccomplex_t IT_0443 = IT_0435*IT_0442;
    const ccomplex_t IT_0444 = IT_0439 + IT_0443;
    const ccomplex_t IT_0445 = 2*IT_0050;
    const ccomplex_t IT_0446 = IT_0072*IT_0445;
    const ccomplex_t IT_0447 = IT_0077*IT_0446;
    const ccomplex_t IT_0448 = -IT_0447;
    const ccomplex_t IT_0449 = 2*IT_0448;
    const ccomplex_t IT_0450 = conj(IT_0326)*IT_0449;
    const ccomplex_t IT_0451 = 2*conj(IT_0448);
    const ccomplex_t IT_0452 = IT_0326*IT_0451;
    const ccomplex_t IT_0453 = s_24*s_34*IT_0041*IT_0135;
    const ccomplex_t IT_0454 = s_12*s_14*s_34*IT_0043;
    const ccomplex_t IT_0455 = s_14*s_34*IT_0041*IT_0312;
    const ccomplex_t IT_0456 = s_12*s_24*s_34*IT_0043;
    const ccomplex_t IT_0457 = pow(s_34, 2);
    const ccomplex_t IT_0458 = s_12*s_13*IT_0041*IT_0457;
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = s_12*s_23*IT_0041*IT_0457;
    const ccomplex_t IT_0461 = -IT_0460;
    const ccomplex_t IT_0462 = conj(IT_0330)*IT_0449;
    const ccomplex_t IT_0463 = IT_0330*IT_0451;
    const ccomplex_t IT_0464 = 2*IT_0369 + 2*IT_0387;
    const ccomplex_t IT_0465 = IT_0043*IT_0312;
    const ccomplex_t IT_0466 = IT_0043*IT_0316;
    const ccomplex_t IT_0467 = IT_0465 + IT_0466;
    const ccomplex_t IT_0468 = s_23*s_24*s_34*IT_0041;
    const ccomplex_t IT_0469 = -IT_0075 + -IT_0468;
    const ccomplex_t IT_0470 = IT_0467 + IT_0469;
    const ccomplex_t IT_0471 = m_N_1*IT_0470;
    const ccomplex_t IT_0472 = 2*conj(IT_0132);
    const ccomplex_t IT_0473 = (-2)*conj(IT_0078);
    const ccomplex_t IT_0474 = IT_0472 + IT_0473;
    const ccomplex_t IT_0475 = s_14*s_23*s_34*IT_0041;
    const ccomplex_t IT_0476 = s_12 + IT_0046 + IT_0475;
    const ccomplex_t IT_0477 = m_N_2*IT_0476;
    const ccomplex_t IT_0478 = 2*conj(IT_0185);
    const ccomplex_t IT_0479 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0480 = IT_0478 + IT_0479;
    const ccomplex_t IT_0481 = s_14*s_24;
    const ccomplex_t IT_0482 = (-2)*IT_0481;
    const ccomplex_t IT_0483 = IT_0196 + IT_0482;
    const ccomplex_t IT_0484 = IT_0043*IT_0483;
    const ccomplex_t IT_0485 = (s_12 + -IT_0044)*(s_12 + (-0.5)*IT_0484);
    const ccomplex_t IT_0486 = (-2)*IT_0485;
    const ccomplex_t IT_0487 = 2*conj(IT_0040);
    const ccomplex_t IT_0488 = s_13*s_23*IT_0041*IT_0137;
    const ccomplex_t IT_0489 = s_14*s_24*IT_0041*IT_0135;
    const ccomplex_t IT_0490 = s_12*IT_0018;
    const ccomplex_t IT_0491 = s_14*s_24*IT_0018*IT_0043;
    const ccomplex_t IT_0492 = -IT_0491;
    const ccomplex_t IT_0493 = s_13*s_23*IT_0018*IT_0043;
    const ccomplex_t IT_0494 = -IT_0493;
    const ccomplex_t IT_0495 = s_12*s_13*s_14*s_34*IT_0041;
    const ccomplex_t IT_0496 = -IT_0495;
    const ccomplex_t IT_0497 = IT_0488 + IT_0489 + IT_0490 + IT_0492 + IT_0494
       + IT_0496;
    const ccomplex_t IT_0498 = s_13*s_14*s_23*s_24*IT_0041;
    const ccomplex_t IT_0499 = IT_0041*IT_0135*IT_0316;
    const ccomplex_t IT_0500 = IT_0018*IT_0075;
    const ccomplex_t IT_0501 = IT_0043*IT_0075*IT_0135;
    const ccomplex_t IT_0502 = -IT_0501;
    const ccomplex_t IT_0503 = IT_0018*IT_0043*IT_0316;
    const ccomplex_t IT_0504 = -IT_0503;
    const ccomplex_t IT_0505 = s_12*s_13*s_24*s_34*IT_0041;
    const ccomplex_t IT_0506 = -IT_0505;
    const ccomplex_t IT_0507 = IT_0498 + IT_0499 + IT_0500 + IT_0502 + IT_0504
       + IT_0506;
    const ccomplex_t IT_0508 = IT_0062*IT_0445;
    const ccomplex_t IT_0509 = IT_0077*IT_0508;
    const ccomplex_t IT_0510 = (-2)*IT_0254 + (-2)*IT_0259 + -IT_0509;
    const ccomplex_t IT_0511 = IT_0286*conj(IT_0510);
    const ccomplex_t IT_0512 = IT_0288*IT_0510;
    const ccomplex_t IT_0513 = IT_0294*conj(IT_0448);
    const ccomplex_t IT_0514 = IT_0296*IT_0448;
    const ccomplex_t IT_0515 = 2*IT_0040;
    const ccomplex_t IT_0516 = conj(IT_0293)*IT_0515;
    const ccomplex_t IT_0517 = IT_0293*IT_0487;
    const ccomplex_t IT_0518 = IT_0041*IT_0137*IT_0312;
    const ccomplex_t IT_0519 = IT_0043*IT_0075*IT_0137;
    const ccomplex_t IT_0520 = -IT_0519;
    const ccomplex_t IT_0521 = IT_0018*IT_0043*IT_0312;
    const ccomplex_t IT_0522 = -IT_0521;
    const ccomplex_t IT_0523 = s_12*s_14*s_23*s_34*IT_0041;
    const ccomplex_t IT_0524 = -IT_0523;
    const ccomplex_t IT_0525 = IT_0498 + IT_0500 + IT_0518 + IT_0520 + IT_0522
       + IT_0524;
    const ccomplex_t IT_0526 = s_14*s_24*IT_0041*IT_0312;
    const ccomplex_t IT_0527 = s_13*s_23*IT_0041*IT_0316;
    const ccomplex_t IT_0528 = s_12*IT_0075;
    const ccomplex_t IT_0529 = s_14*s_24*IT_0043*IT_0075;
    const ccomplex_t IT_0530 = -IT_0529;
    const ccomplex_t IT_0531 = s_13*s_23*IT_0043*IT_0075;
    const ccomplex_t IT_0532 = -IT_0531;
    const ccomplex_t IT_0533 = s_12*s_23*s_24*s_34*IT_0041;
    const ccomplex_t IT_0534 = -IT_0533;
    const ccomplex_t IT_0535 = IT_0526 + IT_0527 + IT_0528 + IT_0530 + IT_0532
       + IT_0534;
    const ccomplex_t IT_0536 = 2*IT_0413 + 2*IT_0430 + -IT_0509;
    const ccomplex_t IT_0537 = IT_0286*conj(IT_0536);
    const ccomplex_t IT_0538 = IT_0288*IT_0536;
    const ccomplex_t IT_0539 = 2*IT_0464;
    const ccomplex_t IT_0540 = conj(IT_0293)*IT_0539;
    const ccomplex_t IT_0541 = 2*conj(IT_0464);
    const ccomplex_t IT_0542 = IT_0293*IT_0541;
    const ccomplex_t IT_0543 = s_24*s_34*IT_0018*IT_0043;
    const ccomplex_t IT_0544 = -IT_0543;
    const ccomplex_t IT_0545 = IT_0298 + IT_0302 + IT_0307 + IT_0453 + IT_0459
       + IT_0544;
    const ccomplex_t IT_0546 = s_14*s_34*IT_0043*IT_0075;
    const ccomplex_t IT_0547 = -IT_0546;
    const ccomplex_t IT_0548 = IT_0310 + IT_0315 + IT_0321 + IT_0455 + IT_0461
       + IT_0547;
    const ccomplex_t IT_0549 = 2*IT_0291;
    const ccomplex_t IT_0550 = conj(IT_0448)*IT_0549;
    const ccomplex_t IT_0551 = 2*conj(IT_0291);
    const ccomplex_t IT_0552 = IT_0448*IT_0551;
    const ccomplex_t IT_0553 = -IT_0018;
    const ccomplex_t IT_0554 = IT_0138 + IT_0553;
    const ccomplex_t IT_0555 = IT_0333*IT_0554;
    const ccomplex_t IT_0556 = -IT_0075;
    const ccomplex_t IT_0557 = IT_0466 + IT_0556;
    const ccomplex_t IT_0558 = IT_0333*IT_0557;
    const ccomplex_t IT_0559 = conj(IT_0285)*IT_0293;
    const ccomplex_t IT_0560 = 6*IT_0559;
    const ccomplex_t IT_0561 = 1.125*IT_0333;
    const ccomplex_t IT_0562 = IT_0369 + IT_0387;
    const ccomplex_t IT_0563 = m_N_2*IT_0562;
    const ccomplex_t IT_0564 = -IT_0563;
    const ccomplex_t IT_0565 = m_N_2*IT_0430;
    const ccomplex_t IT_0566 = m_N_2*IT_0413;
    const ccomplex_t IT_0567 = IT_0366*IT_0427;
    const ccomplex_t IT_0568 = IT_0210*IT_0352*IT_0567;
    const ccomplex_t IT_0569 = IT_0377*IT_0410;
    const ccomplex_t IT_0570 = IT_0200*IT_0370*IT_0569;
    const ccomplex_t IT_0571 = IT_0384*IT_0403;
    const ccomplex_t IT_0572 = IT_0200*IT_0370*IT_0571;
    const ccomplex_t IT_0573 = -IT_0572;
    const ccomplex_t IT_0574 = IT_0359*IT_0420;
    const ccomplex_t IT_0575 = IT_0210*IT_0352*IT_0574;
    const ccomplex_t IT_0576 = -IT_0575;
    const ccomplex_t IT_0577 = 16*conj(IT_0564);
    const ccomplex_t IT_0578 = conj(IT_0568) + conj(IT_0570);
    const ccomplex_t IT_0579 = 16*IT_0573;
    const ccomplex_t IT_0580 = 16*IT_0576;
    const ccomplex_t IT_0581 = IT_0568 + IT_0570;
    const ccomplex_t IT_0582 = 16*conj(IT_0573);
    const ccomplex_t IT_0583 = 16*conj(IT_0576);
    const ccomplex_t IT_0584 = (-16)*conj(IT_0564);
    const ccomplex_t IT_0585 = conj(IT_0573) + conj(IT_0576);
    const ccomplex_t IT_0586 = (-16)*IT_0566;
    const ccomplex_t IT_0587 = IT_0573 + IT_0576;
    const ccomplex_t IT_0588 = (-16)*conj(IT_0565);
    const ccomplex_t IT_0589 = (-16)*conj(IT_0566);
    const ccomplex_t IT_0590 = 16*IT_0229;
    const ccomplex_t IT_0591 = (-16)*IT_0220;
    const ccomplex_t IT_0592 = 16*conj(IT_0229);
    const ccomplex_t IT_0593 = (-16)*conj(IT_0220);
    const ccomplex_t IT_0594 = IT_0564*(16*conj(IT_0565) + 16*conj(IT_0566) + 
      (-16)*conj(IT_0568) + (-16)*conj(IT_0570)) + IT_0565*((-16)*conj(IT_0573) 
      + (-16)*conj(IT_0576) + IT_0577) + IT_0566*IT_0577 + IT_0578*(IT_0579 +
       IT_0580) + IT_0581*(IT_0582 + IT_0583 + IT_0584) + IT_0585*IT_0586 +
       IT_0587*(IT_0588 + IT_0589) + IT_0264*(IT_0590 + IT_0591) + IT_0268*
      (IT_0592 + IT_0593);
    const ccomplex_t IT_0595 = pow(m_W, 4);
    const ccomplex_t IT_0596 = s_12*IT_0595;
    const ccomplex_t IT_0597 = s_13*s_23*IT_0017;
    const ccomplex_t IT_0598 = s_14*s_24*IT_0017;
    const ccomplex_t IT_0599 = s_14*s_23*s_34;
    const ccomplex_t IT_0600 = s_13*s_24*s_34;
    const ccomplex_t IT_0601 = s_12*IT_0457;
    const ccomplex_t IT_0602 = (-2)*IT_0597 + (-2)*IT_0598 + 2*IT_0599 + 2
      *IT_0600 + (-2)*IT_0601;
    const ccomplex_t IT_0603 = IT_0596 + IT_0602;
    const ccomplex_t IT_0604 = IT_0041*IT_0603;
    const ccomplex_t IT_0605 = (-0.25)*IT_0604;
    const ccomplex_t IT_0606 = -IT_0045;
    const ccomplex_t IT_0607 = -IT_0044;
    const ccomplex_t IT_0608 = 8*IT_0573;
    const ccomplex_t IT_0609 = (-8)*IT_0564;
    const ccomplex_t IT_0610 = IT_0608 + IT_0609;
    const ccomplex_t IT_0611 = 8*IT_0220;
    const ccomplex_t IT_0612 = (-8)*IT_0229;
    const ccomplex_t IT_0613 = IT_0611 + IT_0612;
    const ccomplex_t IT_0614 = 8*conj(IT_0220);
    const ccomplex_t IT_0615 = (-8)*conj(IT_0229);
    const ccomplex_t IT_0616 = IT_0614 + IT_0615;
    const ccomplex_t IT_0617 = 8*IT_0565;
    const ccomplex_t IT_0618 = 8*conj(IT_0565);
    const ccomplex_t IT_0619 = 8*conj(IT_0573);
    const ccomplex_t IT_0620 = 8*IT_0256;
    const ccomplex_t IT_0621 = 8*conj(IT_0256);
    const ccomplex_t IT_0622 = 8*IT_0260;
    const ccomplex_t IT_0623 = 8*conj(IT_0260);
    const ccomplex_t IT_0624 = 8*IT_0576;
    const ccomplex_t IT_0625 = (-8)*conj(IT_0564);
    const ccomplex_t IT_0626 = 8*IT_0229;
    const ccomplex_t IT_0627 = (-8)*IT_0220;
    const ccomplex_t IT_0628 = 8*conj(IT_0229);
    const ccomplex_t IT_0629 = (-8)*conj(IT_0220);
    const ccomplex_t IT_0630 = conj(IT_0256)*IT_0609 + conj(IT_0260)*IT_0610 +
       conj(IT_0565)*IT_0612 + conj(IT_0566)*IT_0613 + IT_0565*IT_0615 + IT_0566
      *IT_0616 + conj(IT_0220)*IT_0617 + IT_0220*IT_0618 + (IT_0260 + IT_0263)
      *IT_0619 + conj(IT_0573)*IT_0620 + IT_0587*IT_0621 + 8*conj(IT_0576)*
      (IT_0263 + 0.125*IT_0620 + 0.125*IT_0622) + IT_0576*IT_0623 + conj(IT_0263
      )*(IT_0608 + IT_0609 + IT_0624) + IT_0268*IT_0625 + IT_0578*(IT_0626 +
       IT_0627) + IT_0581*(IT_0628 + IT_0629);
    const ccomplex_t IT_0631 = (-2)*IT_0457;
    const ccomplex_t IT_0632 = IT_0595 + IT_0631;
    const ccomplex_t IT_0633 = IT_0041*IT_0632;
    const ccomplex_t IT_0634 = 0.25*IT_0633;
    const ccomplex_t IT_0635 = 1 + IT_0634;
    const ccomplex_t IT_0636 = IT_0333*IT_0635;
    const ccomplex_t IT_0637 = 8*IT_0566;
    const ccomplex_t IT_0638 = 8*conj(IT_0566);
    const ccomplex_t IT_0639 = (-8)*IT_0256;
    const ccomplex_t IT_0640 = (-8)*IT_0260;
    const ccomplex_t IT_0641 = (-8)*IT_0263;
    const ccomplex_t IT_0642 = (-8)*conj(IT_0256);
    const ccomplex_t IT_0643 = (-8)*conj(IT_0260);
    const ccomplex_t IT_0644 = (-8)*conj(IT_0263);
    const ccomplex_t IT_0645 = conj(IT_0220)*IT_0609 + IT_0585*IT_0613 +
       IT_0587*IT_0616 + IT_0220*IT_0625 + conj(IT_0564)*IT_0626 + IT_0564
      *IT_0628 + IT_0264*(IT_0617 + IT_0637) + IT_0268*(IT_0618 + IT_0638) +
       IT_0578*(IT_0639 + IT_0640 + IT_0641) + IT_0581*(IT_0642 + IT_0643 +
       IT_0644);
    const ccomplex_t IT_0646 = s_34*m_N_1*IT_0043;
    const ccomplex_t IT_0647 = IT_0442*IT_0646;
    const ccomplex_t IT_0648 = 2*IT_0132;
    const ccomplex_t IT_0649 = (-2)*IT_0078;
    const ccomplex_t IT_0650 = IT_0648 + IT_0649;
    const ccomplex_t IT_0651 = s_34*m_N_2*IT_0043;
    const ccomplex_t IT_0652 = IT_0438*IT_0651;
    const ccomplex_t IT_0653 = 2*IT_0147;
    const ccomplex_t IT_0654 = (-2)*IT_0185;
    const ccomplex_t IT_0655 = IT_0653 + IT_0654;
    const ccomplex_t IT_0656 = 2*IT_0510;
    const ccomplex_t IT_0657 = 2*IT_0536;
    const ccomplex_t IT_0658 = m_N_1*IT_0349;
    const ccomplex_t IT_0659 = (-4)*IT_0658;
    const ccomplex_t IT_0660 = 2*IT_0229;
    const ccomplex_t IT_0661 = IT_0043*IT_0457;
    const ccomplex_t IT_0662 = (s_12 + (-0.5)*IT_0276)*(IT_0017 + -IT_0661);
    const ccomplex_t IT_0663 = (-2)*IT_0662;
    const ccomplex_t IT_0664 = -IT_0661;
    const ccomplex_t IT_0665 = IT_0017 + IT_0664;
    const ccomplex_t IT_0666 = IT_0333*IT_0665;
    const ccomplex_t IT_0667 = (-6)*IT_0283;
    const ccomplex_t IT_0668 = 2*conj(IT_0510);
    const ccomplex_t IT_0669 = 2*conj(IT_0536);
    const ccomplex_t IT_0670 = 2*conj(IT_0229);
    const ccomplex_t IT_0671 = (-6)*conj(IT_0283);
    const ccomplex_t IT_0672 = -IT_0272;
    const ccomplex_t IT_0673 = s_13 + IT_0672;
    const ccomplex_t IT_0674 = IT_0333*IT_0673;
    const ccomplex_t IT_0675 = -IT_0279;
    const ccomplex_t IT_0676 = s_23 + IT_0675;
    const ccomplex_t IT_0677 = IT_0333*IT_0676;
    const ccomplex_t IT_0678 = IT_0674 + IT_0677;
    const ccomplex_t IT_0679 = conj(IT_0285)*IT_0291;
    const ccomplex_t IT_0680 = (-6)*IT_0679;
    const ccomplex_t IT_0681 = conj(IT_0293)*IT_0667;
    const ccomplex_t IT_0682 = IT_0293*IT_0671;
    const ccomplex_t IT_0683 = IT_0017*IT_0075;
    const ccomplex_t IT_0684 = s_23*s_24*s_34*IT_0043;
    const ccomplex_t IT_0685 = 0.5*IT_0316 + (-0.5)*IT_0683 + -IT_0684;
    const ccomplex_t IT_0686 = IT_0312 + IT_0685;
    const ccomplex_t IT_0687 = IT_0043*IT_0686;
    const ccomplex_t IT_0688 = -IT_0466;
    const ccomplex_t IT_0689 = IT_0075 + IT_0687 + IT_0688;
    const ccomplex_t IT_0690 = m_N_1*IT_0689;
    const ccomplex_t IT_0691 = (-4)*IT_0285;
    const ccomplex_t IT_0692 = (-4)*conj(IT_0285);
    const ccomplex_t IT_0693 = conj(IT_0565) + conj(IT_0566);
    const ccomplex_t IT_0694 = (-4)*IT_0293;
    const ccomplex_t IT_0695 = IT_0565 + IT_0566;
    const ccomplex_t IT_0696 = (-4)*conj(IT_0293);
    const ccomplex_t IT_0697 = IT_0041*IT_0457;
    const ccomplex_t IT_0698 = 0.5*IT_0697;
    const ccomplex_t IT_0699 = 1 + IT_0698;
    const ccomplex_t IT_0700 = IT_0333*IT_0699;
    const ccomplex_t IT_0701 = (-4)*conj(IT_0132);
    const ccomplex_t IT_0702 = 0.25*IT_0701;
    const ccomplex_t IT_0703 = conj(IT_0078) + IT_0702;
    const ccomplex_t IT_0704 = (-4)*conj(IT_0078);
    const ccomplex_t IT_0705 = 0.25*IT_0704;
    const ccomplex_t IT_0706 = conj(IT_0132) + IT_0705;
    const ccomplex_t IT_0707 = -conj(IT_0147);
    const ccomplex_t IT_0708 = conj(IT_0185) + IT_0707;
    const ccomplex_t IT_0709 = (-4)*IT_0132;
    const ccomplex_t IT_0710 = s_12*IT_0699;
    const ccomplex_t IT_0711 = (-4)*conj(IT_0185);
    const ccomplex_t IT_0712 = 0.25*IT_0711;
    const ccomplex_t IT_0713 = conj(IT_0147) + IT_0712;
    const ccomplex_t IT_0714 = (-4)*IT_0185;
    const ccomplex_t IT_0715 = -IT_0481;
    const ccomplex_t IT_0716 = IT_0196 + IT_0715;
    const ccomplex_t IT_0717 = IT_0043*IT_0716;
    const ccomplex_t IT_0718 = -IT_0045 + (-0.25)*IT_0717;
    const ccomplex_t IT_0719 = s_12 + IT_0718;
    const ccomplex_t IT_0720 = m_N_2*IT_0719;
    const ccomplex_t IT_0721 = (-4)*IT_0720;
    const ccomplex_t IT_0722 = (-2)*IT_0285;
    const ccomplex_t IT_0723 = (-2)*conj(IT_0285);
    const ccomplex_t IT_0724 = IT_0047*IT_0333;
    const ccomplex_t IT_0725 = IT_0333*IT_0476;
    const ccomplex_t IT_0726 = conj(IT_0448)*IT_0722;
    const ccomplex_t IT_0727 = IT_0448*IT_0723;
    const ccomplex_t IT_0728 = IT_0017*IT_0018;
    const ccomplex_t IT_0729 = -IT_0728;
    const ccomplex_t IT_0730 = IT_0137 + IT_0729;
    const ccomplex_t IT_0731 = IT_0043*IT_0730;
    const ccomplex_t IT_0732 = (-4)*IT_0554;
    const ccomplex_t IT_0733 = IT_0731 + IT_0732;
    const ccomplex_t IT_0734 = m_N_2*IT_0733;
    const ccomplex_t IT_0735 = (-2)*IT_0293;
    const ccomplex_t IT_0736 = (-2)*conj(IT_0293);
    const ccomplex_t IT_0737 = s_14*s_23*s_34*IT_0043;
    const ccomplex_t IT_0738 = (-2)*IT_0273 + -IT_0481 + 2*IT_0737;
    const ccomplex_t IT_0739 = IT_0196 + IT_0738;
    const ccomplex_t IT_0740 = IT_0043*IT_0739;
    const ccomplex_t IT_0741 = -IT_0045 + (-0.5)*IT_0740;
    const ccomplex_t IT_0742 = s_12 + IT_0741;
    const ccomplex_t IT_0743 = m_N_1*IT_0742;
    const ccomplex_t IT_0744 = (-2)*IT_0743;
    const ccomplex_t IT_0745 = 2*IT_0078;
    const ccomplex_t IT_0746 = (-2)*IT_0132;
    const ccomplex_t IT_0747 = IT_0745 + IT_0746;
    const ccomplex_t IT_0748 = -(IT_0018 + -IT_0136)*(s_12 + (-0.5)*IT_0484);
    const ccomplex_t IT_0749 = (-2)*IT_0748;
    const ccomplex_t IT_0750 = (-2)*IT_0040;
    const ccomplex_t IT_0751 = 2*IT_0316 + -IT_0683 + (-2)*IT_0684;
    const ccomplex_t IT_0752 = IT_0312 + IT_0751;
    const ccomplex_t IT_0753 = IT_0043*IT_0752;
    const ccomplex_t IT_0754 = -IT_0075 + (-0.5)*IT_0753;
    const ccomplex_t IT_0755 = IT_0465 + IT_0754;
    const ccomplex_t IT_0756 = m_N_1*IT_0755;
    const ccomplex_t IT_0757 = (-2)*IT_0756;
    const ccomplex_t IT_0758 = (-0.5)*IT_0276;
    const ccomplex_t IT_0759 = (-0.5)*IT_0484;
    const ccomplex_t IT_0760 = (-2)*IT_0597 + (-2)*IT_0598 + 4*IT_0600;
    const ccomplex_t IT_0761 = IT_0596 + IT_0760;
    const ccomplex_t IT_0762 = IT_0041*IT_0761;
    const ccomplex_t IT_0763 = 0.25*IT_0762;
    const ccomplex_t IT_0764 = s_12 + IT_0758 + IT_0759 + IT_0763;
    const ccomplex_t IT_0765 = s_13*s_24*s_34*IT_0043;
    const ccomplex_t IT_0766 = -IT_0273 + (-2)*IT_0481 + 2*IT_0765;
    const ccomplex_t IT_0767 = IT_0196 + IT_0766;
    const ccomplex_t IT_0768 = IT_0043*IT_0767;
    const ccomplex_t IT_0769 = -IT_0044 + (-0.5)*IT_0768;
    const ccomplex_t IT_0770 = s_12 + IT_0769;
    const ccomplex_t IT_0771 = m_N_1*IT_0770;
    const ccomplex_t IT_0772 = (-2)*IT_0771;
    const ccomplex_t IT_0773 = (-2)*conj(IT_0040);
    const ccomplex_t IT_0774 = -IT_0273;
    const ccomplex_t IT_0775 = IT_0196 + IT_0774;
    const ccomplex_t IT_0776 = IT_0043*IT_0775;
    const ccomplex_t IT_0777 = -IT_0044 + (-0.25)*IT_0776;
    const ccomplex_t IT_0778 = s_12 + IT_0777;
    const ccomplex_t IT_0779 = m_N_1*IT_0778;
    const ccomplex_t IT_0780 = (-4)*IT_0779;
    const ccomplex_t IT_0781 = IT_0449 + IT_0515;
    const ccomplex_t IT_0782 = IT_0451 + IT_0487;
    const ccomplex_t IT_0783 = (-2)*IT_0448;
    const ccomplex_t IT_0784 = IT_0750 + IT_0783;
    const ccomplex_t IT_0785 = (-2)*conj(IT_0448);
    const ccomplex_t IT_0786 = IT_0773 + IT_0785;
    const ccomplex_t IT_0787 = 2*IT_0564;
    const ccomplex_t IT_0788 = 2*conj(IT_0564);
    const ccomplex_t IT_0789 = (-2)*IT_0510;
    const ccomplex_t IT_0790 = (-2)*conj(IT_0510);
    const ccomplex_t IT_0791 = s_13*s_14*s_34*IT_0043;
    const ccomplex_t IT_0792 = 2*IT_0137 + -IT_0728 + (-2)*IT_0791;
    const ccomplex_t IT_0793 = IT_0135 + IT_0792;
    const ccomplex_t IT_0794 = IT_0043*IT_0793;
    const ccomplex_t IT_0795 = -IT_0018 + (-0.5)*IT_0794;
    const ccomplex_t IT_0796 = IT_0136 + IT_0795;
    const ccomplex_t IT_0797 = m_N_2*IT_0796;
    const ccomplex_t IT_0798 = (-2)*IT_0797;
    const ccomplex_t IT_0799 = 2*IT_0185;
    const ccomplex_t IT_0800 = (-2)*IT_0147;
    const ccomplex_t IT_0801 = IT_0799 + IT_0800;
    const ccomplex_t IT_0802 = -(IT_0075 + -IT_0465)*(s_12 + (-0.5)*IT_0484);
    const ccomplex_t IT_0803 = (-2)*IT_0802;
    const ccomplex_t IT_0804 = (-2)*IT_0464;
    const ccomplex_t IT_0805 = -IT_0273 + (-2)*IT_0481 + 2*IT_0737;
    const ccomplex_t IT_0806 = IT_0196 + IT_0805;
    const ccomplex_t IT_0807 = IT_0043*IT_0806;
    const ccomplex_t IT_0808 = -IT_0044 + (-0.5)*IT_0807;
    const ccomplex_t IT_0809 = s_12 + IT_0808;
    const ccomplex_t IT_0810 = m_N_2*IT_0809;
    const ccomplex_t IT_0811 = (-2)*IT_0810;
    const ccomplex_t IT_0812 = IT_0783 + IT_0804;
    const ccomplex_t IT_0813 = (-2)*conj(IT_0464);
    const ccomplex_t IT_0814 = IT_0785 + IT_0813;
    const ccomplex_t IT_0815 = (-2)*IT_0536;
    const ccomplex_t IT_0816 = (-2)*conj(IT_0536);
    const ccomplex_t IT_0817 = -IT_0683;
    const ccomplex_t IT_0818 = IT_0312 + IT_0817;
    const ccomplex_t IT_0819 = IT_0043*IT_0818;
    const ccomplex_t IT_0820 = IT_0465 + IT_0556;
    const ccomplex_t IT_0821 = (-4)*IT_0820;
    const ccomplex_t IT_0822 = IT_0819 + IT_0821;
    const ccomplex_t IT_0823 = m_N_1*IT_0822;
    const ccomplex_t IT_0824 = IT_0449 + IT_0539;
    const ccomplex_t IT_0825 = IT_0451 + IT_0541;
    const ccomplex_t IT_0826 = (-2)*IT_0564;
    const ccomplex_t IT_0827 = (-2)*conj(IT_0564);
    const ccomplex_t IT_0828 = 2*IT_0283;
    const ccomplex_t IT_0829 = IT_0334 + IT_0335;
    const ccomplex_t IT_0830 = s_13*IT_0043*IT_0457;
    const ccomplex_t IT_0831 = (-2)*IT_0830;
    const ccomplex_t IT_0832 = IT_0829 + IT_0831;
    const ccomplex_t IT_0833 = IT_0043*IT_0832;
    const ccomplex_t IT_0834 = -IT_0272 + 0.5*IT_0833;
    const ccomplex_t IT_0835 = s_13 + IT_0834;
    const ccomplex_t IT_0836 = m_N_2*IT_0835;
    const ccomplex_t IT_0837 = 2*IT_0836;
    const ccomplex_t IT_0838 = (-2)*IT_0229;
    const ccomplex_t IT_0839 = (-2)*conj(IT_0229);
    const ccomplex_t IT_0840 = s_13*m_N_2;
    const ccomplex_t IT_0841 = IT_0699*IT_0840;
    const ccomplex_t IT_0842 = m_N_1*m_N_2*IT_0017;
    const ccomplex_t IT_0843 = IT_0699*IT_0842;
    const ccomplex_t IT_0844 = IT_0298 + IT_0302 + IT_0304 + IT_0306 + IT_0307
       + IT_0453 + IT_0454 + IT_0459;
    const ccomplex_t IT_0845 = IT_0310 + IT_0315 + IT_0318 + IT_0320 + IT_0321
       + IT_0455 + IT_0456 + IT_0461;
    const ccomplex_t IT_0846 = s_23*m_N_1;
    const ccomplex_t IT_0847 = IT_0699*IT_0846;
    const ccomplex_t IT_0848 = (IT_0196 + (-2)*IT_0273)*(1 + 0.5*IT_0697);
    const ccomplex_t IT_0849 = (-4)*IT_0326;
    const ccomplex_t IT_0850 = IT_0043*IT_0334;
    const ccomplex_t IT_0851 = s_14*IT_0017;
    const ccomplex_t IT_0852 = s_13*s_34;
    const ccomplex_t IT_0853 = (-2)*IT_0852;
    const ccomplex_t IT_0854 = IT_0851 + IT_0853;
    const ccomplex_t IT_0855 = IT_0188*IT_0854;
    const ccomplex_t IT_0856 = (-0.5)*IT_0850 + (-0.25)*IT_0855;
    const ccomplex_t IT_0857 = s_13 + IT_0856;
    const ccomplex_t IT_0858 = m_N_2*IT_0857;
    const ccomplex_t IT_0859 = (-4)*IT_0858;
    const ccomplex_t IT_0860 = s_34*IT_0043;
    const ccomplex_t IT_0861 = IT_0189 + IT_0190;
    const ccomplex_t IT_0862 = s_13*s_23*s_34*IT_0043;
    const ccomplex_t IT_0863 = (-2)*IT_0862;
    const ccomplex_t IT_0864 = IT_0861 + IT_0863;
    const ccomplex_t IT_0865 = IT_0860*IT_0864;
    const ccomplex_t IT_0866 = (-2)*IT_0283;
    const ccomplex_t IT_0867 = (-0.25)*IT_0439*IT_0515 + (-0.25)*IT_0443
      *IT_0539 + (IT_0078 + -IT_0132)*IT_0841 + IT_0330*IT_0843 + (-0.25)
      *IT_0656*IT_0844 + (-0.25)*IT_0657*IT_0845 + -(IT_0147 + 0.25*IT_0714)
      *IT_0847 + (-0.25)*IT_0848*IT_0849 + (-0.25)*IT_0838*IT_0859 + (-0.25)
      *IT_0865*IT_0866;
    const ccomplex_t IT_0868 = m_N_2*IT_0778;
    const ccomplex_t IT_0869 = (-4)*IT_0868;
    const ccomplex_t IT_0870 = IT_0135 + IT_0729;
    const ccomplex_t IT_0871 = IT_0043*IT_0870;
    const ccomplex_t IT_0872 = IT_0136 + IT_0553;
    const ccomplex_t IT_0873 = (-4)*IT_0872;
    const ccomplex_t IT_0874 = IT_0871 + IT_0873;
    const ccomplex_t IT_0875 = m_N_2*IT_0874;
    const ccomplex_t IT_0876 = IT_0316 + IT_0817;
    const ccomplex_t IT_0877 = IT_0043*IT_0876;
    const ccomplex_t IT_0878 = (-4)*IT_0557;
    const ccomplex_t IT_0879 = IT_0877 + IT_0878;
    const ccomplex_t IT_0880 = m_N_1*IT_0879;
    const ccomplex_t IT_0881 = IT_0043*IT_0343;
    const ccomplex_t IT_0882 = s_24*IT_0017;
    const ccomplex_t IT_0883 = IT_0188*IT_0882;
    const ccomplex_t IT_0884 = (-0.5)*IT_0881 + 0.25*IT_0883;
    const ccomplex_t IT_0885 = s_23 + IT_0884;
    const ccomplex_t IT_0886 = m_N_1*IT_0885;
    const ccomplex_t IT_0887 = 4*IT_0886;
    const ccomplex_t IT_0888 = 2*conj(IT_0326);
    const ccomplex_t IT_0889 = 2*conj(IT_0330);
    const ccomplex_t IT_0890 = m_N_1*IT_0719;
    const ccomplex_t IT_0891 = (-4)*IT_0890;
    const ccomplex_t IT_0892 = (-2)*IT_0273 + -IT_0481 + 2*IT_0765;
    const ccomplex_t IT_0893 = IT_0196 + IT_0892;
    const ccomplex_t IT_0894 = IT_0043*IT_0893;
    const ccomplex_t IT_0895 = -IT_0045 + (-0.5)*IT_0894;
    const ccomplex_t IT_0896 = s_12 + IT_0895;
    const ccomplex_t IT_0897 = m_N_2*IT_0896;
    const ccomplex_t IT_0898 = (-2)*IT_0897;
    const ccomplex_t IT_0899 = (-2)*conj(IT_0283);
    const ccomplex_t IT_0900 = 8*(conj(IT_0263) + 0.125*IT_0621 + 0.125
      *IT_0623)*IT_0764 + IT_0668*IT_0772 + IT_0757*IT_0816 + IT_0551*IT_0837 +
       IT_0825*IT_0869 + IT_0786*IT_0875 + IT_0288*IT_0880 + IT_0887*IT_0888 +
       IT_0859*IT_0889 + IT_0723*IT_0891 + IT_0736*IT_0898 + IT_0659*IT_0899;
    const ccomplex_t IT_0901 = m_N_2*IT_0340;
    const ccomplex_t IT_0902 = (-4)*IT_0901;
    const ccomplex_t IT_0903 = 2*conj(IT_0283);
    const ccomplex_t IT_0904 = 2*IT_0326;
    const ccomplex_t IT_0905 = 2*IT_0330;
    const ccomplex_t IT_0906 = (-2)*IT_0291;
    const ccomplex_t IT_0907 = IT_0613*IT_0764 + IT_0772*IT_0781 + IT_0757
      *IT_0812 + IT_0828*IT_0837 + IT_0657*IT_0869 + IT_0789*IT_0875 + IT_0294
      *IT_0880 + IT_0735*IT_0891 + IT_0722*IT_0898 + IT_0859*IT_0904 + IT_0887
      *IT_0905 + IT_0659*IT_0906;
    const ccomplex_t IT_0908 = (-2)*conj(IT_0291);
    const ccomplex_t IT_0909 = IT_0615*IT_0764 + IT_0772*IT_0782 + IT_0757
      *IT_0814 + IT_0669*IT_0869 + IT_0790*IT_0875 + IT_0296*IT_0880 + IT_0859
      *IT_0888 + IT_0887*IT_0889 + IT_0736*IT_0891 + IT_0723*IT_0898 + IT_0837
      *IT_0903 + IT_0659*IT_0908;
    const ccomplex_t IT_0910 = IT_0343 + IT_0344;
    const ccomplex_t IT_0911 = s_23*IT_0043*IT_0457;
    const ccomplex_t IT_0912 = (-2)*IT_0911;
    const ccomplex_t IT_0913 = IT_0910 + IT_0912;
    const ccomplex_t IT_0914 = IT_0043*IT_0913;
    const ccomplex_t IT_0915 = -IT_0279 + 0.5*IT_0914;
    const ccomplex_t IT_0916 = s_23 + IT_0915;
    const ccomplex_t IT_0917 = m_N_1*IT_0916;
    const ccomplex_t IT_0918 = 2*IT_0917;
    const ccomplex_t IT_0919 = s_23*s_34;
    const ccomplex_t IT_0920 = (-2)*IT_0919;
    const ccomplex_t IT_0921 = IT_0882 + IT_0920;
    const ccomplex_t IT_0922 = IT_0188*IT_0921;
    const ccomplex_t IT_0923 = (-0.5)*IT_0881 + (-0.25)*IT_0922;
    const ccomplex_t IT_0924 = s_23 + IT_0923;
    const ccomplex_t IT_0925 = m_N_1*IT_0924;
    const ccomplex_t IT_0926 = (-4)*IT_0925;
    const ccomplex_t IT_0927 = (-2)*IT_0326;
    const ccomplex_t IT_0928 = (-2)*conj(IT_0326);
    const ccomplex_t IT_0929 = (-2)*IT_0330;
    const ccomplex_t IT_0930 = (-2)*conj(IT_0330);
    const ccomplex_t IT_0931 = IT_0188*IT_0851;
    const ccomplex_t IT_0932 = (-0.5)*IT_0850 + 0.25*IT_0931;
    const ccomplex_t IT_0933 = s_13 + IT_0932;
    const ccomplex_t IT_0934 = m_N_2*IT_0933;
    const ccomplex_t IT_0935 = 4*IT_0934;
    const ccomplex_t IT_0936 = s_23*IT_0595;
    const ccomplex_t IT_0937 = IT_0041*IT_0936;
    const ccomplex_t IT_0938 = 0.666666666666667*IT_0881 + (-0.166666666666667
      )*IT_0937;
    const ccomplex_t IT_0939 = s_23 + IT_0938;
    const ccomplex_t IT_0940 = m_N_1*IT_0939;
    const ccomplex_t IT_0941 = (-6)*IT_0940;
    const ccomplex_t IT_0942 = (-2)*IT_0388;
    const ccomplex_t IT_0943 = (-2)*conj(IT_0388);
    const ccomplex_t IT_0944 = s_13*IT_0595;
    const ccomplex_t IT_0945 = s_14*s_34*IT_0017;
    const ccomplex_t IT_0946 = (-2)*IT_0945;
    const ccomplex_t IT_0947 = IT_0944 + IT_0946;
    const ccomplex_t IT_0948 = IT_0041*IT_0947;
    const ccomplex_t IT_0949 = (-0.666666666666667)*IT_0272 + (
      -0.166666666666667)*IT_0948;
    const ccomplex_t IT_0950 = s_13 + IT_0949;
    const ccomplex_t IT_0951 = m_N_2*IT_0950;
    const ccomplex_t IT_0952 = (-6)*IT_0951;
    const ccomplex_t IT_0953 = (-2)*IT_0395;
    const ccomplex_t IT_0954 = (-2)*conj(IT_0395);
    const ccomplex_t IT_0955 = (-2)*IT_0335;
    const ccomplex_t IT_0956 = IT_0334 + IT_0955;
    const ccomplex_t IT_0957 = IT_0043*IT_0956;
    const ccomplex_t IT_0958 = 0.166666666666667*IT_0948 + (-0.666666666666667
      )*IT_0957;
    const ccomplex_t IT_0959 = s_13 + IT_0958;
    const ccomplex_t IT_0960 = m_N_2*IT_0959;
    const ccomplex_t IT_0961 = 6*IT_0960;
    const ccomplex_t IT_0962 = (-2)*IT_0431;
    const ccomplex_t IT_0963 = (-2)*conj(IT_0431);
    const ccomplex_t IT_0964 = (-4)*IT_0700;
    const ccomplex_t IT_0965 = -IT_0132;
    const ccomplex_t IT_0966 = -conj(IT_0132);
    const ccomplex_t IT_0967 = -IT_0185;
    const ccomplex_t IT_0968 = -conj(IT_0185);
    const ccomplex_t IT_0969 = -IT_0147;
    const ccomplex_t IT_0970 = conj(IT_0132)*IT_0564 + IT_0132*conj(IT_0564) +
       conj(IT_0078)*(IT_0256 + IT_0260 + IT_0263 + -IT_0564 + IT_0573 + IT_0576
      ) + IT_0078*(conj(IT_0256) + conj(IT_0260) + conj(IT_0263) + -conj(IT_0564
      ) + conj(IT_0573) + conj(IT_0576)) + (IT_0229 + IT_0568 + IT_0570)*(conj
      (IT_0185) + IT_0707) + (conj(IT_0256) + conj(IT_0260) + conj(IT_0263) +
       IT_0585)*IT_0965 + (IT_0256 + IT_0260 + IT_0263 + IT_0587)*IT_0966 + 
      (conj(IT_0220) + IT_0693)*(IT_0147 + IT_0967) + (IT_0220 + IT_0695)*(conj
      (IT_0147) + IT_0968) + (conj(IT_0229) + conj(IT_0568) + conj(IT_0570))*
      (IT_0185 + IT_0969);
    const ccomplex_t IT_0971 = IT_0189 + IT_0191;
    const ccomplex_t IT_0972 = -IT_0190;
    const ccomplex_t IT_0973 = IT_0971 + IT_0972;
    const ccomplex_t IT_0974 = IT_0188*IT_0973;
    const ccomplex_t IT_0975 = 0.25*IT_0974;
    const ccomplex_t IT_0976 = s_12 + IT_0198 + IT_0975;
    const ccomplex_t IT_0977 = IT_0185 + IT_0969;
    const ccomplex_t IT_0978 = 0.125*IT_0693;
    const ccomplex_t IT_0979 = 0.125*IT_0581;
    const ccomplex_t IT_0980 = 8*IT_0976*(0.125*IT_0232*IT_0578 + 0.125
      *IT_0267*IT_0585 + 0.125*IT_0271*IT_0587 + 0.125*IT_0226*IT_0695 + IT_0564
      *IT_0708 + conj(IT_0564)*IT_0977 + IT_0223*IT_0978 + IT_0235*IT_0979);
    const ccomplex_t IT_0981 = (-2)*IT_0597 + (-2)*IT_0598 + 4*IT_0599;
    const ccomplex_t IT_0982 = IT_0596 + IT_0981;
    const ccomplex_t IT_0983 = IT_0041*IT_0982;
    const ccomplex_t IT_0984 = 0.25*IT_0983;
    const ccomplex_t IT_0985 = 8*IT_0568;
    const ccomplex_t IT_0986 = (-8)*IT_0565;
    const ccomplex_t IT_0987 = (-8)*IT_0566;
    const ccomplex_t IT_0988 = (-8)*conj(IT_0565);
    const ccomplex_t IT_0989 = (-8)*conj(IT_0566);
    const ccomplex_t IT_0990 = IT_0564*conj(IT_0564) + 0.125*conj(IT_0573)
      *IT_0610 + 0.125*IT_0566*IT_0618 + 0.125*conj(IT_0576)*(IT_0608 + IT_0609 
      + IT_0624) + 0.125*IT_0573*IT_0625 + 0.125*IT_0576*(IT_0619 + IT_0625) +
       0.125*conj(IT_0566)*IT_0637 + IT_0617*IT_0978 + IT_0578*(IT_0570 + 0.125
      *IT_0985 + 0.125*IT_0986 + 0.125*IT_0987) + IT_0979*(IT_0988 + IT_0989);
    const ccomplex_t IT_0991 = IT_0333*IT_0872;
    const ccomplex_t IT_0992 = s_12 + IT_0607;
    const ccomplex_t IT_0993 = IT_0333*IT_0992;
    const ccomplex_t IT_0994 = -IT_0993;
    const ccomplex_t IT_0995 = IT_0991 + IT_0994;
    const ccomplex_t IT_0996 = IT_0333*IT_0820;
    const ccomplex_t IT_0997 = -IT_0996;
    const ccomplex_t IT_0998 = IT_0993 + IT_0997;
    const ccomplex_t IT_0999 = -IT_0536;
    const ccomplex_t IT_1000 = 0.166666666666667*IT_0048;
    const ccomplex_t IT_1001 = 0.166666666666667*IT_0471;
    const ccomplex_t IT_1002 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_1003 = 0.166666666666667*IT_0477;
    const ccomplex_t IT_1004 = 0.166666666666667*IT_0749;
    const ccomplex_t IT_1005 = 4*IT_0448;
    const ccomplex_t IT_1006 = 0.166666666666667*IT_0486;
    const ccomplex_t IT_1007 = 0.166666666666667*IT_0772;
    const ccomplex_t IT_1008 = 0.166666666666667*IT_0798;
    const ccomplex_t IT_1009 = 0.166666666666667*IT_0803;
    const ccomplex_t IT_1010 = 0.166666666666667*IT_0757;
    const ccomplex_t IT_1011 = 0.166666666666667*IT_0811;
    const ccomplex_t IT_1012 = IT_0142*IT_0333;
    const ccomplex_t IT_1013 = IT_0333*IT_0470;
    const ccomplex_t IT_1014 = IT_1012 + IT_1013;
    const ccomplex_t IT_1015 = 0.166666666666667*IT_1014;
    const ccomplex_t IT_1016 = 0.166666666666667*IT_0444;
    const ccomplex_t IT_1017 = IT_0510*IT_0995 + IT_0998*IT_0999 + IT_0747
      *IT_1000 + IT_0650*IT_1001 + IT_0655*IT_1002 + IT_0801*IT_1003 + IT_0784
      *IT_1004 + (IT_0515 + IT_0539 + IT_1005)*IT_1006 + IT_0838*IT_1007 +
       IT_0826*IT_1008 + IT_0812*IT_1009 + IT_0660*IT_1010 + IT_0787*IT_1011 +
       IT_0286*IT_1015 + IT_0828*IT_1016;
    const ccomplex_t IT_1018 = -conj(IT_0536);
    const ccomplex_t IT_1019 = conj(IT_0510)*IT_0995 + IT_0134*IT_1000 +
       IT_0474*IT_1001 + IT_0187*IT_1002 + IT_0480*IT_1003 + IT_0773*IT_1004 + 
      (IT_0487 + IT_0541)*IT_1006 + IT_0839*IT_1007 + IT_0827*IT_1008 + IT_0813
      *IT_1009 + IT_0670*IT_1010 + IT_0788*IT_1011 + IT_0288*IT_1015 + IT_0903
      *IT_1016 + IT_0998*IT_1018;
    const ccomplex_t IT_1020 = -IT_0464;
    const ccomplex_t IT_1021 = 0.166666666666667*IT_0875;
    const ccomplex_t IT_1022 = 0.166666666666667*IT_0545;
    const ccomplex_t IT_1023 = -conj(IT_0464);
    const ccomplex_t IT_1024 = 0.5*IT_0137 + (-0.5)*IT_0728 + -IT_0791;
    const ccomplex_t IT_1025 = IT_0135 + IT_1024;
    const ccomplex_t IT_1026 = IT_0043*IT_1025;
    const ccomplex_t IT_1027 = -IT_0138;
    const ccomplex_t IT_1028 = IT_0018 + IT_1026 + IT_1027;
    const ccomplex_t IT_1029 = m_N_2*IT_1028;
    const ccomplex_t IT_1030 = 0.25*IT_0471;
    const ccomplex_t IT_1031 = 0.25*IT_0143;
    const ccomplex_t IT_1032 = m_N_2*IT_0047;
    const ccomplex_t IT_1033 = 0.25*IT_1032;
    const ccomplex_t IT_1034 = m_N_1*IT_0476;
    const ccomplex_t IT_1035 = 0.25*IT_1034;
    const ccomplex_t IT_1036 = 0.25*IT_1012;
    const ccomplex_t IT_1037 = 0.25*IT_0724;
    const ccomplex_t IT_1038 = 0.25*IT_0725;
    const ccomplex_t IT_1039 = 0.25*IT_1013;
    const ccomplex_t IT_1040 = 0.25*IT_0744;
    const ccomplex_t IT_1041 = 0.25*IT_0898;
    const ccomplex_t IT_1042 = (-1.5)*conj(IT_0291)*IT_0678 + conj(IT_0564)
      *IT_0690 + (conj(IT_0220) + -conj(IT_0229))*IT_1029 + IT_0480*IT_1030 +
       IT_0134*IT_1031 + IT_0474*IT_1033 + IT_0187*IT_1035 + IT_0487*IT_1036 +
       IT_0773*IT_1037 + IT_0813*IT_1038 + IT_0541*IT_1039 + IT_0827*IT_1040 +
       IT_0670*IT_1041;
    const ccomplex_t IT_1043 = s_12 + IT_0606;
    const ccomplex_t IT_1044 = IT_0333*IT_1043;
    const ccomplex_t IT_1045 = 0.25*IT_0734;
    const ccomplex_t IT_1046 = 0.25*IT_0721;
    const ccomplex_t IT_1047 = 0.25*IT_0891;
    const ccomplex_t IT_1048 = 0.25*IT_0880;
    const ccomplex_t IT_1049 = (-0.333333333333333)*IT_0578*IT_0690 + (
      -0.333333333333333)*IT_0264*IT_1029 + (-0.333333333333333)*IT_0474*IT_1030
       + (-0.333333333333333)*IT_0187*IT_1031 + (-0.333333333333333)*IT_0480
      *IT_1033 + (-0.333333333333333)*IT_0134*IT_1035 + (-0.333333333333333)
      *IT_0668*IT_1036 + (-0.333333333333333)*IT_0790*IT_1037 + (
      -0.333333333333333)*IT_0816*IT_1038 + (-0.333333333333333)*IT_0669*IT_1039
       + conj(IT_0285)*IT_1044 + (-0.333333333333333)*IT_0788*IT_1045 + (
      -0.333333333333333)*IT_0827*IT_1046 + (-0.333333333333333)*IT_0670*IT_1047
       + (-0.333333333333333)*IT_0839*IT_1048;
    const ccomplex_t IT_1050 = (s_12 + -IT_0045)*(s_12 + (-0.5)*IT_0276);
    const ccomplex_t IT_1051 = (-2)*IT_1050;
    const ccomplex_t IT_1052 = -(IT_0018 + -IT_0138)*(s_12 + (-0.5)*IT_0276);
    const ccomplex_t IT_1053 = (-2)*IT_1052;
    const ccomplex_t IT_1054 = -(s_12 + (-0.5)*IT_0276)*(IT_0075 + -IT_0466);
    const ccomplex_t IT_1055 = (-2)*IT_1054;
    const ccomplex_t IT_1056 = IT_1053 + IT_1055;
    const ccomplex_t IT_1057 = 0.25*IT_1056;
    const ccomplex_t IT_1058 = IT_0564*IT_0690 + (IT_0220 + -IT_0229)*IT_1029 
      + IT_0801*IT_1030 + IT_0747*IT_1031 + IT_0650*IT_1033 + IT_0655*IT_1035 +
       IT_0515*IT_1036 + IT_0750*IT_1037 + IT_0804*IT_1038 + IT_0539*IT_1039 +
       IT_0826*IT_1040 + IT_0660*IT_1041 + IT_0285*IT_1051 + IT_0722*IT_1057;
    const ccomplex_t IT_1059 = IT_0581*IT_0690 + IT_0268*IT_1029 + IT_0650
      *IT_1030 + IT_0655*IT_1031 + IT_0801*IT_1033 + IT_0747*IT_1035 + IT_0656
      *IT_1036 + IT_0789*IT_1037 + IT_0815*IT_1038 + IT_0657*IT_1039 + 1.5
      *IT_0285*(IT_0555 + IT_0558 + (-2)*IT_1044) + IT_0787*IT_1045 + IT_0826
      *IT_1046 + IT_0660*IT_1047 + IT_0838*IT_1048 + IT_0293*IT_1051 + IT_0735
      *IT_1057;
    const ccomplex_t IT_1060 = (-2)*IT_0597;
    const ccomplex_t IT_1061 = IT_0596 + IT_1060;
    const ccomplex_t IT_1062 = IT_0043*IT_1061;
    const ccomplex_t IT_1063 = 0.5*IT_1062;
    const ccomplex_t IT_1064 = (-2)*IT_0597 + 2*IT_0599 + 2*IT_0600 + (-2)
      *IT_0601;
    const ccomplex_t IT_1065 = IT_0596 + IT_1064;
    const ccomplex_t IT_1066 = IT_0043*IT_1065;
    const ccomplex_t IT_1067 = (-0.5)*IT_1066;
    const ccomplex_t IT_1068 = s_14*s_23*IT_0017;
    const ccomplex_t IT_1069 = s_13*s_24*IT_0017;
    const ccomplex_t IT_1070 = IT_1068 + IT_1069;
    const ccomplex_t IT_1071 = s_12*s_34*IT_0017;
    const ccomplex_t IT_1072 = -IT_1071;
    const ccomplex_t IT_1073 = IT_1070 + IT_1072;
    const ccomplex_t IT_1074 = IT_0188*IT_1073;
    const ccomplex_t IT_1075 = 0.5*IT_1074;
    const ccomplex_t IT_1076 = IT_0196 + IT_0274 + IT_1063 + IT_1067 + IT_1075;
    const ccomplex_t IT_1077 = -IT_0298;
    const ccomplex_t IT_1078 = IT_0299 + IT_1077;
    const ccomplex_t IT_1079 = s_12*s_13*IT_0017;
    const ccomplex_t IT_1080 = s_23*IT_0017*IT_0018;
    const ccomplex_t IT_1081 = s_23*IT_0137;
    const ccomplex_t IT_1082 = s_12*s_14*s_34;
    const ccomplex_t IT_1083 = -IT_1080 + 2*IT_1081 + (-2)*IT_1082;
    const ccomplex_t IT_1084 = IT_1079 + IT_1083;
    const ccomplex_t IT_1085 = IT_0043*IT_1084;
    const ccomplex_t IT_1086 = (-0.5)*IT_1085;
    const ccomplex_t IT_1087 = IT_1078 + IT_1086;
    const ccomplex_t IT_1088 = -IT_0310;
    const ccomplex_t IT_1089 = IT_0311 + IT_1088;
    const ccomplex_t IT_1090 = s_12*s_23*IT_0017;
    const ccomplex_t IT_1091 = s_13*IT_0017*IT_0075;
    const ccomplex_t IT_1092 = s_14*s_23*s_24;
    const ccomplex_t IT_1093 = s_14*s_34*IT_0075;
    const ccomplex_t IT_1094 = -IT_1091 + (-2)*IT_1092 + 2*IT_1093;
    const ccomplex_t IT_1095 = IT_1090 + IT_1094;
    const ccomplex_t IT_1096 = IT_0043*IT_1095;
    const ccomplex_t IT_1097 = (-0.5)*IT_1096;
    const ccomplex_t IT_1098 = IT_1089 + IT_1097;
    const ccomplex_t IT_1099 = IT_0017*IT_0333;
    const ccomplex_t IT_1100 = 9*IT_1099;
    const ccomplex_t IT_1101 = 0.25*IT_1100;
    const ccomplex_t IT_1102 = (-2)*IT_0344;
    const ccomplex_t IT_1103 = IT_0343 + IT_1102;
    const ccomplex_t IT_1104 = IT_0043*IT_1103;
    const ccomplex_t IT_1105 = 0.5*IT_0881 + (-0.5)*IT_0883 + (-0.5)*IT_1104;
    const ccomplex_t IT_1106 = s_23 + IT_1105;
    const ccomplex_t IT_1107 = m_N_1*IT_1106;
    const ccomplex_t IT_1108 = 2*IT_1107;
    const ccomplex_t IT_1109 = 0.25*IT_1108;
    const ccomplex_t IT_1110 = IT_0801*IT_1109;
    const ccomplex_t IT_1111 = 0.5*IT_0850 + (-0.5)*IT_0931 + (-0.5)*IT_0957;
    const ccomplex_t IT_1112 = s_13 + IT_1111;
    const ccomplex_t IT_1113 = m_N_2*IT_1112;
    const ccomplex_t IT_1114 = 2*IT_1113;
    const ccomplex_t IT_1115 = 0.25*IT_1114;
    const ccomplex_t IT_1116 = 2*IT_0565;
    const ccomplex_t IT_1117 = 2*IT_0566;
    const ccomplex_t IT_1118 = IT_1116 + IT_1117;
    const ccomplex_t IT_1119 = 0.25*IT_0961;
    const ccomplex_t IT_1120 = (-2)*IT_0256;
    const ccomplex_t IT_1121 = (-2)*IT_0260;
    const ccomplex_t IT_1122 = (-2)*IT_0263;
    const ccomplex_t IT_1123 = IT_1120 + IT_1121 + IT_1122;
    const ccomplex_t IT_1124 = s_24*s_34*IT_0017;
    const ccomplex_t IT_1125 = (-2)*IT_1124;
    const ccomplex_t IT_1126 = IT_0936 + IT_1125;
    const ccomplex_t IT_1127 = IT_0041*IT_1126;
    const ccomplex_t IT_1128 = (-0.666666666666667)*IT_0279 + (
      -0.166666666666667)*IT_1127;
    const ccomplex_t IT_1129 = s_23 + IT_1128;
    const ccomplex_t IT_1130 = m_N_1*IT_1129;
    const ccomplex_t IT_1131 = (-6)*IT_1130;
    const ccomplex_t IT_1132 = 0.25*IT_1131;
    const ccomplex_t IT_1133 = (-2)*IT_0220;
    const ccomplex_t IT_1134 = IT_0660 + IT_1133;
    const ccomplex_t IT_1135 = m_N_2*IT_0043;
    const ccomplex_t IT_1136 = IT_0043*IT_0944;
    const ccomplex_t IT_1137 = (-2)*IT_0335 + (-0.5)*IT_1136;
    const ccomplex_t IT_1138 = IT_0334 + IT_1137;
    const ccomplex_t IT_1139 = IT_1135*IT_1138;
    const ccomplex_t IT_1140 = (-2)*IT_1139;
    const ccomplex_t IT_1141 = 0.25*IT_1140;
    const ccomplex_t IT_1142 = 2*IT_0573;
    const ccomplex_t IT_1143 = 2*IT_0576;
    const ccomplex_t IT_1144 = IT_0826 + IT_1142 + IT_1143;
    const ccomplex_t IT_1145 = (-10)*s_23;
    const ccomplex_t IT_1146 = IT_0937 + IT_1145;
    const ccomplex_t IT_1147 = m_N_1*IT_1146;
    const ccomplex_t IT_1148 = 0.25*IT_1147;
    const ccomplex_t IT_1149 = -IT_0457;
    const ccomplex_t IT_1150 = IT_0595 + IT_1149;
    const ccomplex_t IT_1151 = IT_0043*IT_1150;
    const ccomplex_t IT_1152 = -IT_0017 + 0.25*IT_1151;
    const ccomplex_t IT_1153 = IT_0661 + IT_1152;
    const ccomplex_t IT_1154 = IT_0333*IT_1153;
    const ccomplex_t IT_1155 = 4*IT_1154;
    const ccomplex_t IT_1156 = 0.25*IT_1155;
    const ccomplex_t IT_1157 = IT_0549*IT_1156;
    const ccomplex_t IT_1158 = s_34*IT_0017;
    const ccomplex_t IT_1159 = IT_0188*IT_1158;
    const ccomplex_t IT_1160 = IT_0043*IT_0632;
    const ccomplex_t IT_1161 = (-0.333333333333333)*IT_1159 + (
      -0.333333333333333)*IT_1160;
    const ccomplex_t IT_1162 = IT_0017 + IT_1161;
    const ccomplex_t IT_1163 = IT_0333*IT_1162;
    const ccomplex_t IT_1164 = 3*IT_1163;
    const ccomplex_t IT_1165 = 0.25*IT_1164;
    const ccomplex_t IT_1166 = IT_0905*IT_1165;
    const ccomplex_t IT_1167 = (-4)*IT_0283;
    const ccomplex_t IT_1168 = s_12*IT_0043*IT_0457;
    const ccomplex_t IT_1169 = -IT_1168;
    const ccomplex_t IT_1170 = IT_0596 + IT_0599 + IT_0600;
    const ccomplex_t IT_1171 = (-2)*IT_0597 + -IT_0601;
    const ccomplex_t IT_1172 = IT_1170 + IT_1171;
    const ccomplex_t IT_1173 = IT_0043*IT_1172;
    const ccomplex_t IT_1174 = (-0.5)*IT_1173;
    const ccomplex_t IT_1175 = -IT_0737;
    const ccomplex_t IT_1176 = IT_0196 + IT_0765 + IT_1169 + IT_1174 + IT_1175;
    const ccomplex_t IT_1177 = 0.25*IT_1176;
    const ccomplex_t IT_1178 = s_13*s_14*s_24;
    const ccomplex_t IT_1179 = IT_1081 + IT_1178;
    const ccomplex_t IT_1180 = -IT_1080 + -IT_1082;
    const ccomplex_t IT_1181 = IT_1179 + IT_1180;
    const ccomplex_t IT_1182 = IT_0043*IT_1181;
    const ccomplex_t IT_1183 = (-0.5)*IT_1182;
    const ccomplex_t IT_1184 = -IT_0454;
    const ccomplex_t IT_1185 = IT_0298 + IT_0303 + IT_0306 + IT_1183 + IT_1184;
    const ccomplex_t IT_1186 = 0.25*IT_1185;
    const ccomplex_t IT_1187 = (-8)*IT_0285;
    const ccomplex_t IT_1188 = s_13*IT_0316;
    const ccomplex_t IT_1189 = IT_1092 + IT_1188;
    const ccomplex_t IT_1190 = s_12*s_24*s_34;
    const ccomplex_t IT_1191 = -IT_1091 + -IT_1190;
    const ccomplex_t IT_1192 = IT_1189 + IT_1191;
    const ccomplex_t IT_1193 = IT_0043*IT_1192;
    const ccomplex_t IT_1194 = (-0.25)*IT_1193;
    const ccomplex_t IT_1195 = (-0.5)*IT_0456;
    const ccomplex_t IT_1196 = 0.5*IT_0317;
    const ccomplex_t IT_1197 = (-0.5)*IT_0319;
    const ccomplex_t IT_1198 = (-0.5)*IT_0310;
    const ccomplex_t IT_1199 = IT_0311 + IT_1194 + IT_1195 + IT_1196 + IT_1197
       + IT_1198;
    const ccomplex_t IT_1200 = 0.25*IT_1199;
    const ccomplex_t IT_1201 = IT_0326*IT_1076 + -IT_0510*IT_1087 + IT_0536
      *IT_1098 + IT_0942*IT_1101 + IT_1110 + IT_0747*IT_1115 + IT_1118*IT_1119 +
       IT_1123*IT_1132 + IT_1134*IT_1141 + IT_1144*IT_1148 + IT_1157 + IT_1166 +
       IT_1167*IT_1177 + IT_0691*IT_1186 + IT_1187*IT_1200;
    const ccomplex_t IT_1202 = s_24*s_34*IT_0018;
    const ccomplex_t IT_1203 = -IT_1080 + (-2)*IT_1178 + 2*IT_1202;
    const ccomplex_t IT_1204 = IT_1079 + IT_1203;
    const ccomplex_t IT_1205 = IT_0043*IT_1204;
    const ccomplex_t IT_1206 = (-0.5)*IT_1205;
    const ccomplex_t IT_1207 = IT_1078 + IT_1206;
    const ccomplex_t IT_1208 = m_N_1*IT_0043;
    const ccomplex_t IT_1209 = IT_0043*IT_0936;
    const ccomplex_t IT_1210 = (-2)*IT_0344 + (-0.5)*IT_1209;
    const ccomplex_t IT_1211 = IT_0343 + IT_1210;
    const ccomplex_t IT_1212 = IT_1208*IT_1211;
    const ccomplex_t IT_1213 = (-2)*IT_1212;
    const ccomplex_t IT_1214 = IT_0041*IT_0944;
    const ccomplex_t IT_1215 = 0.666666666666667*IT_0850 + (-0.166666666666667
      )*IT_1214;
    const ccomplex_t IT_1216 = s_13 + IT_1215;
    const ccomplex_t IT_1217 = m_N_2*IT_1216;
    const ccomplex_t IT_1218 = (-6)*IT_1217;
    const ccomplex_t IT_1219 = -IT_1091 + 2*IT_1188 + (-2)*IT_1190;
    const ccomplex_t IT_1220 = IT_1090 + IT_1219;
    const ccomplex_t IT_1221 = IT_0043*IT_1220;
    const ccomplex_t IT_1222 = (-0.5)*IT_1221;
    const ccomplex_t IT_1223 = IT_1089 + IT_1222;
    const ccomplex_t IT_1224 = 0.25*IT_0952;
    const ccomplex_t IT_1225 = (-0.666666666666667)*IT_1104 +
       0.166666666666667*IT_1127;
    const ccomplex_t IT_1226 = s_23 + IT_1225;
    const ccomplex_t IT_1227 = m_N_1*IT_1226;
    const ccomplex_t IT_1228 = 6*IT_1227;
    const ccomplex_t IT_1229 = 0.25*IT_1228;
    const ccomplex_t IT_1230 = (-10)*s_13;
    const ccomplex_t IT_1231 = IT_1214 + IT_1230;
    const ccomplex_t IT_1232 = m_N_2*IT_1231;
    const ccomplex_t IT_1233 = 0.25*IT_1232;
    const ccomplex_t IT_1234 = m_N_1*m_N_2*IT_0043;
    const ccomplex_t IT_1235 = (-0.25)*IT_0595;
    const ccomplex_t IT_1236 = IT_0457 + IT_1235;
    const ccomplex_t IT_1237 = IT_1234*IT_1236;
    const ccomplex_t IT_1238 = (-4)*IT_1237;
    const ccomplex_t IT_1239 = 0.25*IT_1238;
    const ccomplex_t IT_1240 = -IT_0765;
    const ccomplex_t IT_1241 = IT_0196 + IT_0737 + IT_1169 + IT_1174 + IT_1240;
    const ccomplex_t IT_1242 = 0.25*IT_1241;
    const ccomplex_t IT_1243 = (-0.25)*IT_1182;
    const ccomplex_t IT_1244 = (-0.5)*IT_0454;
    const ccomplex_t IT_1245 = (-0.5)*IT_0303;
    const ccomplex_t IT_1246 = 0.5*IT_0305;
    const ccomplex_t IT_1247 = (-0.5)*IT_0298;
    const ccomplex_t IT_1248 = IT_0299 + IT_1243 + IT_1244 + IT_1245 + IT_1246
       + IT_1247;
    const ccomplex_t IT_1249 = 0.25*IT_1248;
    const ccomplex_t IT_1250 = (-0.5)*IT_1193;
    const ccomplex_t IT_1251 = -IT_0456;
    const ccomplex_t IT_1252 = IT_0310 + IT_0318 + IT_0319 + IT_1250 + IT_1251;
    const ccomplex_t IT_1253 = 0.25*IT_1252;
    const ccomplex_t IT_1254 = IT_0326*IT_1076 + IT_1110 + IT_0747*IT_1115 +
       IT_1157 + IT_1166 + IT_0510*IT_1207 + 0.25*IT_1144*IT_1213 + 0.25*IT_1133
      *IT_1218 + IT_0999*IT_1223 + IT_1118*IT_1224 + IT_1123*IT_1229 + IT_0660
      *IT_1233 + IT_0942*IT_1239 + IT_1167*IT_1242 + IT_1187*IT_1249 + IT_0691
      *IT_1253;
    const ccomplex_t IT_1255 = 0.25*IT_0439;
    const ccomplex_t IT_1256 = 0.25*IT_0443;
    const ccomplex_t IT_1257 = 0.25*IT_0887;
    const ccomplex_t IT_1258 = 0.25*IT_0865;
    const ccomplex_t IT_1259 = 0.25*IT_0844;
    const ccomplex_t IT_1260 = 0.25*IT_0845;
    const ccomplex_t IT_1261 = -(IT_0078 + 0.25*IT_0709)*IT_0847 + IT_0330
      *IT_0848 + (-0.25)*IT_0843*IT_0849 + -IT_0841*IT_0977 + -IT_0656*IT_1255 +
       -IT_0657*IT_1256 + -IT_0838*IT_1257 + -IT_0906*IT_1258 + -IT_0515*IT_1259
       + -IT_0539*IT_1260;
    const ccomplex_t IT_1262 = 2*conj(IT_0565);
    const ccomplex_t IT_1263 = 2*conj(IT_0566);
    const ccomplex_t IT_1264 = (-2)*conj(IT_0568);
    const ccomplex_t IT_1265 = (-0.5)*IT_1108;
    const ccomplex_t IT_1266 = (-0.5)*IT_1114;
    const ccomplex_t IT_1267 = 2*conj(IT_0573);
    const ccomplex_t IT_1268 = 2*conj(IT_0576);
    const ccomplex_t IT_1269 = IT_0827 + IT_1267 + IT_1268;
    const ccomplex_t IT_1270 = (-0.5)*IT_0952;
    const ccomplex_t IT_1271 = (-2)*conj(IT_0220);
    const ccomplex_t IT_1272 = IT_0670 + IT_1271;
    const ccomplex_t IT_1273 = (-0.5)*IT_1228;
    const ccomplex_t IT_1274 = (-2)*conj(IT_0256);
    const ccomplex_t IT_1275 = (-2)*conj(IT_0260);
    const ccomplex_t IT_1276 = (-2)*conj(IT_0263);
    const ccomplex_t IT_1277 = IT_1274 + IT_1275 + IT_1276;
    const ccomplex_t IT_1278 = (-0.5)*IT_1232;
    const ccomplex_t IT_1279 = (-0.5)*IT_1155;
    const ccomplex_t IT_1280 = (-0.5)*IT_1164;
    const ccomplex_t IT_1281 = (-0.5)*IT_1238;
    const ccomplex_t IT_1282 = (-4)*conj(IT_0291);
    const ccomplex_t IT_1283 = (-0.5)*IT_1241;
    const ccomplex_t IT_1284 = (-8)*conj(IT_0388);
    const ccomplex_t IT_1285 = pow(m_W, 6);
    const ccomplex_t IT_1286 = s_12*IT_1285;
    const ccomplex_t IT_1287 = s_14*s_24*IT_0595;
    const ccomplex_t IT_1288 = (-2)*IT_1287;
    const ccomplex_t IT_1289 = IT_1286 + IT_1288;
    const ccomplex_t IT_1290 = IT_0041*IT_1289;
    const ccomplex_t IT_1291 = 0.25*IT_1290;
    const ccomplex_t IT_1292 = IT_0196 + IT_1291;
    const ccomplex_t IT_1293 = (-0.5)*IT_1292;
    const ccomplex_t IT_1294 = 4*conj(IT_0448);
    const ccomplex_t IT_1295 = 4*conj(IT_0040);
    const ccomplex_t IT_1296 = (-0.5)*IT_1207;
    const ccomplex_t IT_1297 = (-8)*conj(IT_0293);
    const ccomplex_t IT_1298 = (-0.5)*IT_1248;
    const ccomplex_t IT_1299 = (-4)*conj(IT_0448);
    const ccomplex_t IT_1300 = (-4)*conj(IT_0464);
    const ccomplex_t IT_1301 = (-0.5)*IT_1223;
    const ccomplex_t IT_1302 = (-0.5)*IT_1252;
    const ccomplex_t IT_1303 = conj(IT_0330)*IT_1076 + 0.25*IT_0954*IT_1100 + 
      (-0.5)*IT_1213*(conj(IT_0570) + (-0.5)*IT_1262 + (-0.5)*IT_1263 + (-0.5)
      *IT_1264) + (-0.5)*IT_0474*IT_1265 + (-0.5)*IT_0187*IT_1266 + (-0.5)
      *IT_1269*IT_1270 + (-0.5)*IT_1272*IT_1273 + (-0.5)*IT_1277*IT_1278 + (-0.5
      )*IT_0903*IT_1279 + (-0.5)*IT_0888*IT_1280 + (-0.5)*IT_0963*IT_1281 + (
      -0.5)*IT_1282*IT_1283 + (-0.5)*IT_1284*IT_1293 + (-0.5)*(IT_1294 + IT_1295
      )*IT_1296 + (-0.5)*IT_1297*IT_1298 + (-0.5)*(IT_1299 + IT_1300)*IT_1301 + 
      (-0.5)*IT_0696*IT_1302;
    const ccomplex_t IT_1304 = 0.5*IT_1290;
    const ccomplex_t IT_1305 = (-4)*IT_0273;
    const ccomplex_t IT_1306 = IT_0196 + IT_1304 + IT_1305;
    const ccomplex_t IT_1307 = IT_0043*IT_0761;
    const ccomplex_t IT_1308 = -IT_1307;
    const ccomplex_t IT_1309 = IT_1306 + IT_1308;
    const ccomplex_t IT_1310 = (-4)*IT_0039;
    const ccomplex_t IT_1311 = (-0.25)*IT_1310;
    const ccomplex_t IT_1312 = IT_0021 + IT_1311;
    const ccomplex_t IT_1313 = (-4)*IT_0291;
    const ccomplex_t IT_1314 = (-8)*IT_0388;
    const ccomplex_t IT_1315 = 4*IT_0040;
    const ccomplex_t IT_1316 = IT_1005 + IT_1315;
    const ccomplex_t IT_1317 = (-8)*IT_0293;
    const ccomplex_t IT_1318 = (-4)*IT_0448;
    const ccomplex_t IT_1319 = (-4)*IT_0464;
    const ccomplex_t IT_1320 = IT_1318 + IT_1319;
    const ccomplex_t IT_1321 = (-0.5)*IT_1116;
    const ccomplex_t IT_1322 = (-0.5)*IT_1117;
    const ccomplex_t IT_1323 = (-2)*IT_0568;
    const ccomplex_t IT_1324 = (-0.5)*IT_1323;
    const ccomplex_t IT_1325 = IT_0330*IT_1076 + 0.25*IT_0953*IT_1100 + (-0.5)
      *IT_0650*IT_1265 + (-0.5)*IT_0655*IT_1266 + (-0.5)*IT_1144*IT_1270 + (-0.5
      )*IT_1134*IT_1273 + (-0.5)*IT_1123*IT_1278 + (-0.5)*IT_0828*IT_1279 + (
      -0.5)*IT_0904*IT_1280 + (-0.5)*IT_0962*IT_1281 + (-0.5)*IT_0694*IT_1302 + 
      -IT_1309*IT_1312 + (-0.5)*IT_1283*IT_1313 + (-0.5)*IT_1293*IT_1314 + (-0.5
      )*IT_1296*IT_1316 + (-0.5)*IT_1298*IT_1317 + (-0.5)*IT_1301*IT_1320 + (
      -0.5)*IT_1213*(IT_0570 + IT_1321 + IT_1322 + IT_1324);
    const ccomplex_t IT_1326 = IT_1262 + IT_1263;
    const ccomplex_t IT_1327 = (-8)*conj(IT_0285);
    const ccomplex_t IT_1328 = (-4)*conj(IT_0283);
    const ccomplex_t IT_1329 = conj(IT_0326)*IT_1076;
    const ccomplex_t IT_1330 = -IT_1329;
    const ccomplex_t IT_1331 = (-0.25)*IT_1108;
    const ccomplex_t IT_1332 = (-0.25)*IT_1114;
    const ccomplex_t IT_1333 = (-0.25)*IT_1155;
    const ccomplex_t IT_1334 = IT_0551*IT_1333;
    const ccomplex_t IT_1335 = (-0.25)*IT_1164;
    const ccomplex_t IT_1336 = IT_0889*IT_1335;
    const ccomplex_t IT_1337 = -conj(IT_0510)*IT_1207 + conj(IT_0536)*IT_1223 
      + (-0.25)*IT_0670*IT_1232 + (-0.25)*IT_0943*IT_1238 + (-0.25)*IT_0692
      *IT_1252 + (-0.25)*IT_1213*IT_1269 + (-0.25)*IT_1218*IT_1271 + (-0.25)
      *IT_1228*IT_1277 + 2*conj(IT_0431)*IT_1292 + conj(IT_0395)*IT_1309 + (
      -0.25)*IT_0952*IT_1326 + (-0.25)*IT_1248*IT_1327 + (-0.25)*IT_1241*IT_1328
       + IT_1330 + IT_0480*IT_1331 + IT_0134*IT_1332 + IT_1334 + IT_1336;
    const ccomplex_t IT_1338 = IT_0043*IT_0982;
    const ccomplex_t IT_1339 = -IT_1338;
    const ccomplex_t IT_1340 = IT_1306 + IT_1339;
    const ccomplex_t IT_1341 = (-0.25)*IT_0961;
    const ccomplex_t IT_1342 = (-0.25)*IT_1131;
    const ccomplex_t IT_1343 = (-0.25)*IT_1140;
    const ccomplex_t IT_1344 = (-0.25)*IT_1147;
    const ccomplex_t IT_1345 = (-0.25)*IT_1176;
    const ccomplex_t IT_1346 = (-0.25)*IT_1185;
    const ccomplex_t IT_1347 = (-0.25)*IT_1199;
    const ccomplex_t IT_1348 = conj(IT_0510)*IT_1087 + IT_1018*IT_1098 + (
      -0.25)*IT_0943*IT_1100 + 2*conj(IT_0395)*IT_1292 + IT_1330 + IT_0480
      *IT_1331 + IT_0134*IT_1332 + IT_1334 + IT_1336 + conj(IT_0431)*IT_1340 +
       IT_1326*IT_1341 + IT_1277*IT_1342 + IT_1272*IT_1343 + IT_1269*IT_1344 +
       IT_1328*IT_1345 + IT_0692*IT_1346 + IT_1327*IT_1347;
    const ccomplex_t IT_1349 = -conj(IT_0330)*IT_1076 + -IT_1098*(conj(IT_0464
      ) + 0.25*IT_1294) + IT_1087*(conj(IT_0040) + (-0.25)*IT_1299) + IT_0474
      *IT_1331 + IT_0187*IT_1332 + IT_0903*IT_1333 + IT_0888*IT_1335 + conj
      (IT_0388)*IT_1340 + IT_1269*IT_1341 + IT_1272*IT_1342 + IT_1277*IT_1343 +
       IT_1326*IT_1344 + IT_1282*IT_1345 + IT_0696*IT_1346 + IT_1297*IT_1347;
    const ccomplex_t IT_1350 = (-0.166666666666667)*IT_0477;
    const ccomplex_t IT_1351 = (-0.166666666666667)*IT_0486;
    const ccomplex_t IT_1352 = (-0.166666666666667)*IT_0869;
    const ccomplex_t IT_1353 = (-0.166666666666667)*IT_0548;
    const ccomplex_t IT_1354 = IT_0655*IT_1115;
    const ccomplex_t IT_1355 = IT_0330*IT_1076;
    const ccomplex_t IT_1356 = (IT_0464 + 0.25*IT_1005)*IT_1098 + IT_0650
      *IT_1109 + IT_1132*IT_1134 + IT_1123*IT_1141 + IT_1119*IT_1144 + IT_1118
      *IT_1148 + IT_0828*IT_1156 + IT_0904*IT_1165 + IT_0694*IT_1186 + IT_1177
      *IT_1313 + IT_1200*IT_1317 + -IT_1087*(IT_0040 + (-0.25)*IT_1318) +
       IT_1354 + IT_1355;
    const ccomplex_t IT_1357 = -IT_0953*IT_1101 + -IT_0650*IT_1109 + -IT_0828
      *IT_1156 + -IT_0904*IT_1165 + -IT_1144*IT_1224 + -IT_1134*IT_1229 + 
      -IT_1123*IT_1233 + -IT_0962*IT_1239 + -IT_0694*IT_1253 + IT_1309*IT_1312 +
       -IT_1242*IT_1313 + (-0.25)*IT_1292*IT_1314 + (-0.25)*IT_1207*IT_1316 + 
      -IT_1249*IT_1317 + (-0.25)*IT_1223*IT_1320 + 0.5*IT_1213*(IT_0570 +
       IT_1321 + IT_1322 + IT_1324) + -IT_1354 + -IT_1355;
    const ccomplex_t IT_1358 = IT_0784*IT_0875;
    const ccomplex_t IT_1359 = IT_0757*IT_0815;
    const ccomplex_t IT_1360 = IT_0824*IT_0869;
    const ccomplex_t IT_1361 = IT_0735*IT_0898;
    const ccomplex_t IT_1362 = IT_0722*IT_0891;
    const ccomplex_t IT_1363 = IT_0286*IT_0880;
    const ccomplex_t IT_1364 = IT_0659*IT_0866;
    const ccomplex_t IT_1365 = IT_0859*IT_0905;
    const ccomplex_t IT_1366 = IT_0887*IT_0904;
    const ccomplex_t IT_1367 = IT_0668*IT_0772;
    const ccomplex_t IT_1368 = IT_0786*IT_0875;
    const ccomplex_t IT_1369 = IT_0757*IT_0816;
    const ccomplex_t IT_1370 = IT_0825*IT_0869;
    const ccomplex_t IT_1371 = IT_0736*IT_0898;
    const ccomplex_t IT_1372 = IT_0551*IT_0837;
    const ccomplex_t IT_1373 = IT_0723*IT_0891;
    const ccomplex_t IT_1374 = IT_0288*IT_0880;
    const ccomplex_t IT_1375 = IT_0659*IT_0899;
    const ccomplex_t IT_1376 = IT_0859*IT_0889;
    const ccomplex_t IT_1377 = IT_0887*IT_0888;
    const ccomplex_t IT_1378 = IT_0040*(IT_0048*IT_0134 + IT_0143*IT_0187) +
       IT_0199*(conj(IT_0220)*IT_0223 + IT_0220*IT_0226 + conj(IT_0229)*IT_0232 
      + IT_0229*IT_0235 + IT_0264*IT_0267 + IT_0268*IT_0271) + (IT_0278 +
       IT_0281)*(IT_0287 + IT_0289 + IT_0295 + IT_0297) + (IT_0309 + IT_0323)*
      (IT_0327 + IT_0328 + IT_0331 + IT_0332) + (IT_0342 + IT_0351)*(IT_0389 +
       IT_0390 + IT_0391 + IT_0392 + IT_0393 + IT_0394 + IT_0396 + IT_0432 +
       IT_0433 + IT_0434) + IT_0444*(IT_0450 + IT_0452) + (IT_0298 + IT_0302 +
       IT_0304 + IT_0306 + IT_0307 + IT_0310 + IT_0315 + IT_0318 + IT_0320 +
       IT_0321 + IT_0453 + IT_0454 + IT_0455 + IT_0456 + IT_0459 + IT_0461)*
      (IT_0462 + IT_0463) + IT_0464*(IT_0471*IT_0474 + IT_0477*IT_0480 + IT_0486
      *IT_0487) + (IT_0497 + IT_0507)*(IT_0511 + IT_0512 + IT_0513 + IT_0514 +
       IT_0516 + IT_0517) + (IT_0525 + IT_0535)*(IT_0513 + IT_0514 + IT_0537 +
       IT_0538 + IT_0540 + IT_0542) + (IT_0545 + IT_0548)*(IT_0550 + IT_0552) + 
      (IT_0555 + IT_0558)*IT_0560 + IT_0561*IT_0594 + (s_12 + IT_0605 + IT_0606 
      + IT_0607)*IT_0630 + IT_0636*IT_0645 + conj(IT_0291)*(IT_0515*IT_0545 +
       IT_0539*IT_0548 + IT_0647*IT_0650 + IT_0652*IT_0655 + IT_0439*IT_0656 +
       IT_0443*IT_0657 + IT_0659*IT_0660 + IT_0549*IT_0663 + IT_0666*IT_0667) +
       IT_0291*(IT_0487*IT_0545 + IT_0541*IT_0548 + IT_0474*IT_0647 + IT_0187
      *IT_0652 + IT_0439*IT_0668 + IT_0443*IT_0669 + IT_0659*IT_0670 + IT_0666
      *IT_0671) + IT_0678*(IT_0680 + IT_0681 + IT_0682) + IT_0690*(IT_0585
      *IT_0691 + IT_0587*IT_0692 + IT_0693*IT_0694 + IT_0695*IT_0696) + 4
      *IT_0700*(IT_0132*conj(IT_0147) + IT_0185*IT_0703 + IT_0147*IT_0706 +
       IT_0078*IT_0708 + 0.25*conj(IT_0185)*IT_0709) + 4*IT_0710*(IT_0132*conj
      (IT_0132) + IT_0185*conj(IT_0185) + IT_0078*IT_0703 + 0.25*conj(IT_0078)
      *IT_0709 + IT_0147*IT_0713 + 0.25*conj(IT_0147)*IT_0714) + IT_0721*
      (IT_0294*IT_0585 + IT_0296*IT_0587 + IT_0286*IT_0693 + IT_0288*IT_0695 +
       IT_0578*IT_0722 + IT_0581*IT_0723) + (IT_0724 + IT_0725)*(IT_0726 +
       IT_0727) + IT_0734*(IT_0286*IT_0578 + IT_0288*IT_0581 + IT_0693*IT_0722 +
       IT_0695*IT_0723 + IT_0585*IT_0735 + IT_0587*IT_0736) + (IT_0286*IT_0585 +
       IT_0288*IT_0587 + IT_0294*IT_0693 + IT_0296*IT_0695 + IT_0578*IT_0735 +
       IT_0581*IT_0736)*IT_0744 + conj(IT_0040)*(IT_0143*IT_0655 + IT_0048
      *IT_0747 + IT_0749*IT_0750) + conj(IT_0229)*(IT_0539*IT_0757 + IT_0626
      *IT_0764 + IT_0750*IT_0772) + IT_0229*(IT_0541*IT_0757 + IT_0772*IT_0773) 
      + IT_0780*(IT_0578*IT_0781 + IT_0581*IT_0782 + IT_0693*IT_0784 + IT_0695
      *IT_0786 + conj(IT_0510)*IT_0787 + IT_0510*IT_0788 + IT_0585*IT_0789 +
       IT_0587*IT_0790) + (IT_0656*IT_0693 + IT_0668*IT_0695 + conj(IT_0564)
      *IT_0750 + IT_0564*IT_0773 + IT_0585*IT_0781 + IT_0587*IT_0782 + IT_0578
      *IT_0789 + IT_0581*IT_0790)*IT_0798 + conj(IT_0464)*(IT_0486*IT_0515 +
       IT_0471*IT_0650 + IT_0477*IT_0801 + IT_0803*IT_0804) + IT_0811*(IT_0541
      *IT_0564 + IT_0539*conj(IT_0564) + IT_0578*IT_0657 + IT_0581*IT_0669 +
       IT_0585*IT_0812 + IT_0587*IT_0814 + IT_0693*IT_0815 + IT_0695*IT_0816) +
       IT_0823*(IT_0585*IT_0657 + IT_0587*IT_0669 + IT_0578*IT_0812 + IT_0581
      *IT_0814 + IT_0693*IT_0824 + IT_0695*IT_0825 + conj(IT_0536)*IT_0826 +
       IT_0536*IT_0827) + conj(IT_0283)*(IT_0439*IT_0515 + IT_0443*IT_0539 +
       IT_0652*IT_0747 + IT_0647*IT_0801 + IT_0663*IT_0828 + IT_0837*IT_0838) +
       IT_0283*(IT_0439*IT_0487 + IT_0443*IT_0541 + IT_0480*IT_0647 + IT_0134
      *IT_0652 + IT_0837*IT_0839) + (-4)*conj(IT_0326)*IT_0867 + IT_0263*IT_0900
       + IT_0902*(IT_0549*IT_0585 + IT_0551*IT_0587 + conj(IT_0291)*IT_0826 +
       IT_0291*IT_0827 + IT_0693*IT_0828 + IT_0578*IT_0866 + IT_0581*IT_0899 +
       IT_0695*IT_0903) + conj(IT_0220)*IT_0907 + IT_0220*IT_0909 + (IT_0549
      *IT_0578 + IT_0551*IT_0581 + conj(IT_0283)*IT_0787 + IT_0283*IT_0788 +
       IT_0585*IT_0866 + IT_0587*IT_0899 + IT_0693*IT_0906 + IT_0695*IT_0908)
      *IT_0918 + IT_0926*(conj(IT_0326)*IT_0787 + IT_0326*IT_0788 + IT_0581
      *IT_0889 + IT_0578*IT_0905 + IT_0585*IT_0927 + IT_0587*IT_0928 + IT_0693
      *IT_0929 + IT_0695*IT_0930) + (conj(IT_0330)*IT_0787 + IT_0330*IT_0788 +
       IT_0581*IT_0888 + IT_0578*IT_0904 + IT_0693*IT_0927 + IT_0695*IT_0928 +
       IT_0585*IT_0929 + IT_0587*IT_0930)*IT_0935 + IT_0941*(IT_0578*IT_0942 +
       IT_0581*IT_0943) + IT_0952*(IT_0578*IT_0953 + IT_0581*IT_0954) + IT_0961*
      (IT_0578*IT_0962 + IT_0581*IT_0963) + IT_0964*IT_0970 + IT_0980 + 8*(s_12 
      + IT_0758 + IT_0759 + IT_0984)*IT_0990 + 6*conj(IT_0448)*IT_1017 + 6
      *IT_0448*IT_1019 + 6*conj(IT_0510)*(IT_0040*IT_0991 + IT_0655*IT_1000 +
       IT_0747*IT_1002 + IT_0789*IT_1004 + IT_0993*IT_1020 + IT_0660*IT_1021 +
       IT_0828*IT_1022) + 6*IT_0510*(conj(IT_0040)*IT_0991 + IT_0187*IT_1000 +
       IT_0134*IT_1002 + IT_0670*IT_1021 + IT_0903*IT_1022 + IT_0993*IT_1023) +
       4*IT_0285*IT_1042 + (-12)*IT_0293*IT_1049 + 4*conj(IT_0285)*IT_1058 + 4
      *conj(IT_0293)*IT_1059 + 4*conj(IT_0431)*IT_1201 + 4*conj(IT_0395)*IT_1254
       + 4*IT_0330*(IT_0708*IT_0841 + IT_0703*IT_0847 + IT_0668*IT_1255 +
       IT_0669*IT_1256 + IT_0839*IT_1257 + IT_0908*IT_1258 + IT_0487*IT_1259 +
       IT_0541*IT_1260) + 4*IT_0326*(IT_0706*IT_0841 + IT_0713*IT_0847 + 0.25
      *IT_0839*IT_0859 + IT_0487*IT_1255 + IT_0541*IT_1256 + IT_0899*IT_1258 +
       IT_0668*IT_1259 + IT_0669*IT_1260) + (-4)*conj(IT_0330)*IT_1261 + 4*
      (IT_0021 + IT_0039)*IT_1303 + 4*conj(IT_0039)*IT_1325 + (-4)*IT_0395
      *IT_1337 + (-4)*IT_0431*IT_1348 + (-4)*IT_0388*IT_1349 + (-6)*conj(IT_0536
      )*((-0.166666666666667)*IT_0471*IT_0801 + (-0.166666666666667)*IT_0803
      *IT_0815 + IT_0040*IT_0993 + IT_0996*IT_1020 + IT_0650*IT_1350 + IT_0656
      *IT_1351 + IT_0838*IT_1352 + IT_0828*IT_1353) + (-6)*IT_0536*((
      -0.166666666666667)*IT_0471*IT_0480 + conj(IT_0040)*IT_0993 + IT_0996
      *IT_1023 + IT_0474*IT_1350 + IT_0668*IT_1351 + IT_0839*IT_1352 + IT_0903
      *IT_1353) + 4*conj(IT_0388)*IT_1356 + (-4)*conj(IT_0021)*IT_1357 + conj
      (IT_0256)*(IT_0620*IT_0764 + IT_0656*IT_0772 + IT_0549*IT_0837 + IT_1358 +
       IT_1359 + IT_1360 + IT_1361 + IT_1362 + IT_1363 + IT_1364 + IT_1365 +
       IT_1366) + (conj(IT_0260) + conj(IT_0263))*((IT_0620 + IT_0622)*IT_0764 +
       IT_0656*IT_0772 + IT_0549*IT_0837 + IT_1358 + IT_1359 + IT_1360 + IT_1361
       + IT_1362 + IT_1363 + IT_1364 + IT_1365 + IT_1366) + IT_0256*(IT_1367 +
       IT_1368 + IT_1369 + IT_1370 + IT_1371 + IT_1372 + IT_1373 + IT_1374 +
       IT_1375 + IT_1376 + IT_1377) + IT_0260*(IT_0621*IT_0764 + IT_1367 +
       IT_1368 + IT_1369 + IT_1370 + IT_1371 + IT_1372 + IT_1373 + IT_1374 +
       IT_1375 + IT_1376 + IT_1377);
    return create_ccomplex_return(IT_1378);
}

