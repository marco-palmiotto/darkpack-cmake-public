#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_2_N_3_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_2_N_3_to_anti_W_W(
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
    const creal_t m_N_2 = param->m_N_2;
    const creal_t m_N_3 = param->m_N_3;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B2 = param->N_B2;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W2 = param->N_W2;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d2 = param->N_d2;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u2 = param->N_u2;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = m_N_2*m_N_3;
    const ccomplex_t IT_0001 = pow(m_W, -4);
    const ccomplex_t IT_0002 = pow(s_34, 2);
    const ccomplex_t IT_0003 = IT_0001*IT_0002;
    const ccomplex_t IT_0004 = 0.5*IT_0003;
    const ccomplex_t IT_0005 = 1 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = pow(m_N_2, 2);
    const ccomplex_t IT_0008 = pow(m_N_3, 2);
    const ccomplex_t IT_0009 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*IT_0009;
    const ccomplex_t IT_0011 = cos(theta_W);
    const ccomplex_t IT_0012 = sin(theta_W);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = (0 + _Complex_I*1)*e_em*IT_0011*IT_0013;
    const ccomplex_t IT_0015 = -IT_0014;
    const ccomplex_t IT_0016 = IT_0011*IT_0013;
    const ccomplex_t IT_0017 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0018 = IT_0016*IT_0017;
    const ccomplex_t IT_0019 = cpow(IT_0011, -1);
    const ccomplex_t IT_0020 = IT_0012*IT_0019;
    const ccomplex_t IT_0021 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0016*IT_0021;
    const ccomplex_t IT_0024 = IT_0017*IT_0020;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0018 + -IT_0022 + 
      -IT_0023 + IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = m_N_3*IT_0027;
    const ccomplex_t IT_0029 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0030 = IT_0020*IT_0029;
    const ccomplex_t IT_0031 = IT_0016*IT_0029;
    const ccomplex_t IT_0032 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0033 = IT_0020*IT_0032;
    const ccomplex_t IT_0034 = IT_0016*IT_0032;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0030 + IT_0031 + 
      -IT_0033 + -IT_0034);
    const ccomplex_t IT_0036 = -IT_0035;
    const ccomplex_t IT_0037 = (-0.5)*IT_0036;
    const ccomplex_t IT_0038 = IT_0015*IT_0037;
    const ccomplex_t IT_0039 = m_N_2*IT_0038;
    const ccomplex_t IT_0040 = IT_0010*(IT_0028 + -IT_0039);
    const ccomplex_t IT_0041 = cos(alpha);
    const ccomplex_t IT_0042 = IT_0019*IT_0041;
    const ccomplex_t IT_0043 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0044 = IT_0042*IT_0043;
    const ccomplex_t IT_0045 = sin(alpha);
    const ccomplex_t IT_0046 = IT_0013*IT_0045;
    const ccomplex_t IT_0047 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0050 = IT_0042*IT_0049;
    const ccomplex_t IT_0051 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0052 = IT_0046*IT_0051;
    const ccomplex_t IT_0053 = IT_0044 + IT_0048 + IT_0050 + IT_0052;
    const ccomplex_t IT_0054 = IT_0019*IT_0045;
    const ccomplex_t IT_0055 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = IT_0013*IT_0041;
    const ccomplex_t IT_0058 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0059 = IT_0057*IT_0058;
    const ccomplex_t IT_0060 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0061 = IT_0054*IT_0060;
    const ccomplex_t IT_0062 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0063 = IT_0057*IT_0062;
    const ccomplex_t IT_0064 = -IT_0056 + -IT_0059 + -IT_0061 + -IT_0063;
    const ccomplex_t IT_0065 = IT_0053 + IT_0064;
    const ccomplex_t IT_0066 = (0 + _Complex_I*1)*IT_0065;
    const ccomplex_t IT_0067 = 0.5*IT_0066;
    const ccomplex_t IT_0068 = cos(beta);
    const ccomplex_t IT_0069 = sin(beta);
    const ccomplex_t IT_0070 = (0 + _Complex_I*1)*m_W*e_em*IT_0013*(IT_0041
      *IT_0068 + IT_0045*IT_0069);
    const ccomplex_t IT_0071 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0072 = IT_0070*IT_0071;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*IT_0072;
    const ccomplex_t IT_0074 = IT_0067*IT_0073;
    const ccomplex_t IT_0075 = IT_0042*IT_0055;
    const ccomplex_t IT_0076 = IT_0043*IT_0054;
    const ccomplex_t IT_0077 = IT_0042*IT_0060;
    const ccomplex_t IT_0078 = IT_0049*IT_0054;
    const ccomplex_t IT_0079 = IT_0075 + IT_0076 + IT_0077 + IT_0078;
    const ccomplex_t IT_0080 = IT_0047*IT_0057;
    const ccomplex_t IT_0081 = IT_0046*IT_0058;
    const ccomplex_t IT_0082 = IT_0051*IT_0057;
    const ccomplex_t IT_0083 = IT_0046*IT_0062;
    const ccomplex_t IT_0084 = -IT_0080 + -IT_0081 + -IT_0082 + -IT_0083;
    const ccomplex_t IT_0085 = IT_0079 + IT_0084;
    const ccomplex_t IT_0086 = (0 + _Complex_I*1)*IT_0085;
    const ccomplex_t IT_0087 = (-0.5)*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*-1)*m_W*e_em*IT_0013*(IT_0045
      *IT_0068 + -IT_0041*IT_0069);
    const ccomplex_t IT_0089 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0090 = IT_0088*IT_0089;
    const ccomplex_t IT_0091 = (0 + _Complex_I*1)*IT_0090;
    const ccomplex_t IT_0092 = IT_0087*IT_0091;
    const ccomplex_t IT_0093 = -IT_0074 + -IT_0092;
    const ccomplex_t IT_0094 = m_N_3*IT_0038;
    const ccomplex_t IT_0095 = m_N_2*IT_0027;
    const ccomplex_t IT_0096 = IT_0010*(IT_0094 + -IT_0095);
    const ccomplex_t IT_0097 = (-4)*IT_0096;
    const ccomplex_t IT_0098 = 0.25*IT_0097;
    const ccomplex_t IT_0099 = IT_0093 + IT_0098;
    const ccomplex_t IT_0100 = (-4)*conj(IT_0096);
    const ccomplex_t IT_0101 = 0.25*IT_0100;
    const ccomplex_t IT_0102 = conj(IT_0093) + IT_0101;
    const ccomplex_t IT_0103 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0104 = IT_0042*IT_0103;
    const ccomplex_t IT_0105 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0106 = IT_0054*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0108 = IT_0042*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0110 = IT_0054*IT_0109;
    const ccomplex_t IT_0111 = IT_0104 + IT_0106 + IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0113 = IT_0057*IT_0112;
    const ccomplex_t IT_0114 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0115 = IT_0046*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0117 = IT_0057*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0119 = IT_0046*IT_0118;
    const ccomplex_t IT_0120 = -IT_0113 + -IT_0115 + -IT_0117 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0111 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0091*IT_0123;
    const ccomplex_t IT_0125 = IT_0042*IT_0105;
    const ccomplex_t IT_0126 = IT_0046*IT_0112;
    const ccomplex_t IT_0127 = IT_0042*IT_0109;
    const ccomplex_t IT_0128 = IT_0046*IT_0116;
    const ccomplex_t IT_0129 = IT_0125 + IT_0126 + IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = IT_0054*IT_0103;
    const ccomplex_t IT_0131 = IT_0057*IT_0114;
    const ccomplex_t IT_0132 = IT_0054*IT_0107;
    const ccomplex_t IT_0133 = IT_0057*IT_0118;
    const ccomplex_t IT_0134 = -IT_0130 + -IT_0131 + -IT_0132 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0129 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0073*IT_0137;
    const ccomplex_t IT_0139 = -IT_0124 + -IT_0138;
    const ccomplex_t IT_0140 = (-4)*IT_0093;
    const ccomplex_t IT_0141 = 0.25*IT_0140;
    const ccomplex_t IT_0142 = IT_0096 + IT_0141;
    const ccomplex_t IT_0143 = (-4)*conj(IT_0093);
    const ccomplex_t IT_0144 = 0.25*IT_0143;
    const ccomplex_t IT_0145 = conj(IT_0096) + IT_0144;
    const ccomplex_t IT_0146 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0147 = IT_0013*IT_0146;
    const ccomplex_t IT_0148 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0149 = IT_0013*IT_0148;
    const ccomplex_t IT_0150 = 1.4142135623731*IT_0149;
    const ccomplex_t IT_0151 = (0 + _Complex_I*1)*(IT_0147 + 0.5*IT_0150);
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0154 = IT_0013*IT_0153;
    const ccomplex_t IT_0155 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0156 = IT_0013*IT_0155;
    const ccomplex_t IT_0157 = 1.4142135623731*IT_0156;
    const ccomplex_t IT_0158 = (0 + _Complex_I*1)*(IT_0154 + 0.5*IT_0157);
    const ccomplex_t IT_0159 = -IT_0158;
    const ccomplex_t IT_0160 = IT_0152*IT_0159;
    const ccomplex_t IT_0161 = pow(m_W, 2);
    const ccomplex_t IT_0162 = cpow((-2)*s_13 + IT_0007 + IT_0161 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0163 = IT_0160*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*IT_0163;
    const ccomplex_t IT_0165 = conj(N_W2)*e_em*U_Wm1;
    const ccomplex_t IT_0166 = IT_0013*IT_0165;
    const ccomplex_t IT_0167 = conj(N_d2)*U_d1*e_em;
    const ccomplex_t IT_0168 = IT_0013*IT_0167;
    const ccomplex_t IT_0169 = 1.4142135623731*IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*(IT_0166 + 0.5*IT_0169);
    const ccomplex_t IT_0171 = -IT_0170;
    const ccomplex_t IT_0172 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0173 = IT_0013*IT_0172;
    const ccomplex_t IT_0174 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0175 = IT_0013*IT_0174;
    const ccomplex_t IT_0176 = 1.4142135623731*IT_0175;
    const ccomplex_t IT_0177 = (0 + _Complex_I*1)*(IT_0173 + 0.5*IT_0176);
    const ccomplex_t IT_0178 = -IT_0177;
    const ccomplex_t IT_0179 = IT_0171*IT_0178;
    const ccomplex_t IT_0180 = cpow((-2)*s_13 + IT_0007 + IT_0161 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0181 = IT_0179*IT_0180;
    const ccomplex_t IT_0182 = (0 + _Complex_I*1)*IT_0181;
    const ccomplex_t IT_0183 = (-2)*IT_0164 + (-2)*IT_0182;
    const ccomplex_t IT_0184 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0185 = pow(m_W, -2);
    const ccomplex_t IT_0186 = s_13*s_23*IT_0185;
    const ccomplex_t IT_0187 = s_14*s_24*IT_0185;
    const ccomplex_t IT_0188 = -IT_0186 + -IT_0187;
    const ccomplex_t IT_0189 = s_12 + IT_0184 + IT_0188;
    const ccomplex_t IT_0190 = m_N_2*IT_0189;
    const ccomplex_t IT_0191 = 2*conj(IT_0096);
    const ccomplex_t IT_0192 = (-2)*conj(IT_0093);
    const ccomplex_t IT_0193 = IT_0191 + IT_0192;
    const ccomplex_t IT_0194 = pow(s_13, 2);
    const ccomplex_t IT_0195 = IT_0185*IT_0194;
    const ccomplex_t IT_0196 = pow(s_14, 2);
    const ccomplex_t IT_0197 = IT_0185*IT_0196;
    const ccomplex_t IT_0198 = IT_0195 + IT_0197;
    const ccomplex_t IT_0199 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0200 = -IT_0007 + -IT_0199;
    const ccomplex_t IT_0201 = IT_0198 + IT_0200;
    const ccomplex_t IT_0202 = m_N_3*IT_0201;
    const ccomplex_t IT_0203 = 2*conj(IT_0040);
    const ccomplex_t IT_0204 = (-2)*conj(IT_0139);
    const ccomplex_t IT_0205 = IT_0203 + IT_0204;
    const ccomplex_t IT_0206 = s_34*IT_0001;
    const ccomplex_t IT_0207 = s_14*s_23;
    const ccomplex_t IT_0208 = s_13*s_24;
    const ccomplex_t IT_0209 = s_12*s_34;
    const ccomplex_t IT_0210 = -IT_0208 + -IT_0209;
    const ccomplex_t IT_0211 = IT_0207 + IT_0210;
    const ccomplex_t IT_0212 = IT_0206*IT_0211;
    const ccomplex_t IT_0213 = (-0.25)*IT_0212;
    const ccomplex_t IT_0214 = s_12*IT_0161;
    const ccomplex_t IT_0215 = IT_0185*IT_0214;
    const ccomplex_t IT_0216 = (-0.5)*IT_0215;
    const ccomplex_t IT_0217 = s_12 + IT_0213 + IT_0216;
    const ccomplex_t IT_0218 = IT_0164 + IT_0182;
    const ccomplex_t IT_0219 = m_N_2*IT_0218;
    const ccomplex_t IT_0220 = -IT_0219;
    const ccomplex_t IT_0221 = 8*IT_0093;
    const ccomplex_t IT_0222 = (-8)*IT_0096;
    const ccomplex_t IT_0223 = IT_0221 + IT_0222;
    const ccomplex_t IT_0224 = 8*conj(IT_0093);
    const ccomplex_t IT_0225 = (-8)*conj(IT_0096);
    const ccomplex_t IT_0226 = IT_0224 + IT_0225;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0228 = N_W2*e_em*conj(V_Wp1);
    const ccomplex_t IT_0229 = IT_0013*IT_0228;
    const ccomplex_t IT_0230 = N_u2*conj(V_u1)*e_em;
    const ccomplex_t IT_0231 = IT_0013*IT_0230;
    const ccomplex_t IT_0232 = 1.4142135623731*IT_0231;
    const ccomplex_t IT_0233 = (0 + _Complex_I*1)*(IT_0229 + (-0.5)*IT_0232);
    const ccomplex_t IT_0234 = -IT_0233;
    const ccomplex_t IT_0235 = IT_0178*IT_0234;
    const ccomplex_t IT_0236 = IT_0180*IT_0227*IT_0235;
    const ccomplex_t IT_0237 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0238 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0239 = IT_0013*IT_0238;
    const ccomplex_t IT_0240 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0241 = IT_0013*IT_0240;
    const ccomplex_t IT_0242 = 1.4142135623731*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*(IT_0239 + (-0.5)*IT_0242);
    const ccomplex_t IT_0244 = -IT_0243;
    const ccomplex_t IT_0245 = IT_0159*IT_0244;
    const ccomplex_t IT_0246 = IT_0162*IT_0237*IT_0245;
    const ccomplex_t IT_0247 = -IT_0236 + -IT_0246;
    const ccomplex_t IT_0248 = 8*IT_0096;
    const ccomplex_t IT_0249 = (-8)*IT_0093;
    const ccomplex_t IT_0250 = IT_0248 + IT_0249;
    const ccomplex_t IT_0251 = 8*conj(IT_0096);
    const ccomplex_t IT_0252 = (-8)*conj(IT_0093);
    const ccomplex_t IT_0253 = IT_0251 + IT_0252;
    const ccomplex_t IT_0254 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0255 = IT_0013*IT_0254;
    const ccomplex_t IT_0256 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0257 = IT_0013*IT_0256;
    const ccomplex_t IT_0258 = 1.4142135623731*IT_0257;
    const ccomplex_t IT_0259 = (0 + _Complex_I*1)*(IT_0255 + (-0.5)*IT_0258);
    const ccomplex_t IT_0260 = -IT_0259;
    const ccomplex_t IT_0261 = IT_0152*IT_0260;
    const ccomplex_t IT_0262 = IT_0162*IT_0237*IT_0261;
    const ccomplex_t IT_0263 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0264 = IT_0013*IT_0263;
    const ccomplex_t IT_0265 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0266 = IT_0013*IT_0265;
    const ccomplex_t IT_0267 = 1.4142135623731*IT_0266;
    const ccomplex_t IT_0268 = (0 + _Complex_I*1)*(IT_0264 + (-0.5)*IT_0267);
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = IT_0234*IT_0269;
    const ccomplex_t IT_0271 = IT_0180*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*IT_0271;
    const ccomplex_t IT_0273 = m_N_2*IT_0272;
    const ccomplex_t IT_0274 = -IT_0262 + IT_0273;
    const ccomplex_t IT_0275 = IT_0244*IT_0260;
    const ccomplex_t IT_0276 = IT_0162*IT_0275;
    const ccomplex_t IT_0277 = (0 + _Complex_I*1)*IT_0276;
    const ccomplex_t IT_0278 = m_N_2*IT_0277;
    const ccomplex_t IT_0279 = IT_0171*IT_0269;
    const ccomplex_t IT_0280 = IT_0180*IT_0227*IT_0279;
    const ccomplex_t IT_0281 = -IT_0280;
    const ccomplex_t IT_0282 = conj(IT_0274) + conj(IT_0278) + conj(IT_0281);
    const ccomplex_t IT_0283 = 8*IT_0040;
    const ccomplex_t IT_0284 = (-8)*IT_0139;
    const ccomplex_t IT_0285 = IT_0283 + IT_0284;
    const ccomplex_t IT_0286 = IT_0274 + IT_0278 + IT_0281;
    const ccomplex_t IT_0287 = 8*conj(IT_0040);
    const ccomplex_t IT_0288 = (-8)*conj(IT_0139);
    const ccomplex_t IT_0289 = IT_0287 + IT_0288;
    const ccomplex_t IT_0290 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0291 = IT_0013*IT_0290;
    const ccomplex_t IT_0292 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0293 = IT_0013*IT_0292;
    const ccomplex_t IT_0294 = 1.4142135623731*IT_0293;
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*(IT_0291 + (-0.5)*IT_0294);
    const ccomplex_t IT_0296 = -IT_0295;
    const ccomplex_t IT_0297 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0298 = IT_0013*IT_0297;
    const ccomplex_t IT_0299 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0300 = IT_0013*IT_0299;
    const ccomplex_t IT_0301 = 1.4142135623731*IT_0300;
    const ccomplex_t IT_0302 = (0 + _Complex_I*1)*(IT_0298 + (-0.5)*IT_0301);
    const ccomplex_t IT_0303 = -IT_0302;
    const ccomplex_t IT_0304 = IT_0296*IT_0303;
    const ccomplex_t IT_0305 = cpow((-2)*s_23 + IT_0008 + IT_0161 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0306 = IT_0304*IT_0305;
    const ccomplex_t IT_0307 = (0 + _Complex_I*1)*IT_0306;
    const ccomplex_t IT_0308 = cpow((-2)*s_23 + IT_0008 + IT_0161 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0309 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0310 = IT_0013*IT_0309;
    const ccomplex_t IT_0311 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0312 = IT_0013*IT_0311;
    const ccomplex_t IT_0313 = 1.4142135623731*IT_0312;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*(IT_0310 + (-0.5)*IT_0313);
    const ccomplex_t IT_0315 = -IT_0314;
    const ccomplex_t IT_0316 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0317 = IT_0013*IT_0316;
    const ccomplex_t IT_0318 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0319 = IT_0013*IT_0318;
    const ccomplex_t IT_0320 = 1.4142135623731*IT_0319;
    const ccomplex_t IT_0321 = (0 + _Complex_I*1)*(IT_0317 + (-0.5)*IT_0320);
    const ccomplex_t IT_0322 = -IT_0321;
    const ccomplex_t IT_0323 = IT_0315*IT_0322;
    const ccomplex_t IT_0324 = IT_0308*IT_0323;
    const ccomplex_t IT_0325 = (0 + _Complex_I*1)*IT_0324;
    const ccomplex_t IT_0326 = 2*IT_0307 + 2*IT_0325;
    const ccomplex_t IT_0327 = pow(s_23, 2);
    const ccomplex_t IT_0328 = IT_0185*IT_0327;
    const ccomplex_t IT_0329 = pow(s_24, 2);
    const ccomplex_t IT_0330 = IT_0185*IT_0329;
    const ccomplex_t IT_0331 = IT_0328 + IT_0330;
    const ccomplex_t IT_0332 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0333 = -IT_0008 + -IT_0332;
    const ccomplex_t IT_0334 = IT_0331 + IT_0333;
    const ccomplex_t IT_0335 = m_N_2*IT_0334;
    const ccomplex_t IT_0336 = 2*conj(IT_0093);
    const ccomplex_t IT_0337 = (-2)*conj(IT_0096);
    const ccomplex_t IT_0338 = IT_0336 + IT_0337;
    const ccomplex_t IT_0339 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0340 = s_12 + IT_0188 + IT_0339;
    const ccomplex_t IT_0341 = m_N_3*IT_0340;
    const ccomplex_t IT_0342 = 2*conj(IT_0139);
    const ccomplex_t IT_0343 = (-2)*conj(IT_0040);
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = s_14*s_24;
    const ccomplex_t IT_0346 = (-2)*IT_0345;
    const ccomplex_t IT_0347 = IT_0214 + IT_0346;
    const ccomplex_t IT_0348 = IT_0185*IT_0347;
    const ccomplex_t IT_0349 = (s_12 + -IT_0186)*(s_12 + (-0.5)*IT_0348);
    const ccomplex_t IT_0350 = (-2)*IT_0349;
    const ccomplex_t IT_0351 = 2*conj(IT_0183);
    const ccomplex_t IT_0352 = s_13*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0353 = s_24*s_34*IT_0001*IT_0194;
    const ccomplex_t IT_0354 = s_23*IT_0007;
    const ccomplex_t IT_0355 = s_12*s_14*s_34*IT_0185;
    const ccomplex_t IT_0356 = s_14*s_34*IT_0001*IT_0327;
    const ccomplex_t IT_0357 = s_13*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0358 = s_13*IT_0008;
    const ccomplex_t IT_0359 = s_12*s_24*s_34*IT_0185;
    const ccomplex_t IT_0360 = s_12*s_13;
    const ccomplex_t IT_0361 = s_23*IT_0185*IT_0194;
    const ccomplex_t IT_0362 = (-2)*IT_0361;
    const ccomplex_t IT_0363 = IT_0360 + IT_0362;
    const ccomplex_t IT_0364 = s_12*s_23;
    const ccomplex_t IT_0365 = s_13*IT_0185*IT_0327;
    const ccomplex_t IT_0366 = (-2)*IT_0365;
    const ccomplex_t IT_0367 = IT_0364 + IT_0366;
    const ccomplex_t IT_0368 = s_12*s_13*IT_0001*IT_0002;
    const ccomplex_t IT_0369 = -IT_0368;
    const ccomplex_t IT_0370 = s_13*s_14*s_24*IT_0185;
    const ccomplex_t IT_0371 = -IT_0370;
    const ccomplex_t IT_0372 = s_23*IT_0185*IT_0196;
    const ccomplex_t IT_0373 = -IT_0372;
    const ccomplex_t IT_0374 = s_12*s_23*IT_0001*IT_0002;
    const ccomplex_t IT_0375 = -IT_0374;
    const ccomplex_t IT_0376 = s_13*IT_0185*IT_0329;
    const ccomplex_t IT_0377 = -IT_0376;
    const ccomplex_t IT_0378 = s_14*s_23*s_24*IT_0185;
    const ccomplex_t IT_0379 = -IT_0378;
    const ccomplex_t IT_0380 = (-2)*IT_0014;
    const ccomplex_t IT_0381 = IT_0026*IT_0380;
    const ccomplex_t IT_0382 = IT_0010*IT_0381;
    const ccomplex_t IT_0383 = 2*IT_0014;
    const ccomplex_t IT_0384 = IT_0026*IT_0383;
    const ccomplex_t IT_0385 = IT_0010*IT_0384;
    const ccomplex_t IT_0386 = -IT_0385;
    const ccomplex_t IT_0387 = 2*IT_0386;
    const ccomplex_t IT_0388 = conj(IT_0382)*IT_0387;
    const ccomplex_t IT_0389 = 2*conj(IT_0386);
    const ccomplex_t IT_0390 = IT_0382*IT_0389;
    const ccomplex_t IT_0391 = IT_0000*IT_0201;
    const ccomplex_t IT_0392 = IT_0000*IT_0334;
    const ccomplex_t IT_0393 = IT_0010*IT_0038;
    const ccomplex_t IT_0394 = -IT_0393;
    const ccomplex_t IT_0395 = 2*IT_0394;
    const ccomplex_t IT_0396 = conj(IT_0386)*IT_0395;
    const ccomplex_t IT_0397 = 2*conj(IT_0394);
    const ccomplex_t IT_0398 = IT_0386*IT_0397;
    const ccomplex_t IT_0399 = s_13*s_23*IT_0001*IT_0196;
    const ccomplex_t IT_0400 = s_14*s_24*IT_0001*IT_0194;
    const ccomplex_t IT_0401 = s_12*IT_0007;
    const ccomplex_t IT_0402 = s_14*s_24*IT_0007*IT_0185;
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = s_13*s_23*IT_0007*IT_0185;
    const ccomplex_t IT_0405 = -IT_0404;
    const ccomplex_t IT_0406 = s_12*s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0407 = -IT_0406;
    const ccomplex_t IT_0408 = IT_0399 + IT_0400 + IT_0401 + IT_0403 + IT_0405
       + IT_0407;
    const ccomplex_t IT_0409 = s_13*s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0410 = IT_0001*IT_0194*IT_0329;
    const ccomplex_t IT_0411 = IT_0007*IT_0008;
    const ccomplex_t IT_0412 = IT_0008*IT_0185*IT_0194;
    const ccomplex_t IT_0413 = -IT_0412;
    const ccomplex_t IT_0414 = IT_0007*IT_0185*IT_0329;
    const ccomplex_t IT_0415 = -IT_0414;
    const ccomplex_t IT_0416 = s_12*s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0417 = -IT_0416;
    const ccomplex_t IT_0418 = IT_0409 + IT_0410 + IT_0411 + IT_0413 + IT_0415
       + IT_0417;
    const ccomplex_t IT_0419 = IT_0408 + IT_0418;
    const ccomplex_t IT_0420 = IT_0010*IT_0027;
    const ccomplex_t IT_0421 = -IT_0420;
    const ccomplex_t IT_0422 = 2*IT_0183;
    const ccomplex_t IT_0423 = conj(IT_0421)*IT_0422;
    const ccomplex_t IT_0424 = IT_0351*IT_0421;
    const ccomplex_t IT_0425 = 2*IT_0421;
    const ccomplex_t IT_0426 = conj(IT_0386)*IT_0425;
    const ccomplex_t IT_0427 = 2*conj(IT_0421);
    const ccomplex_t IT_0428 = IT_0386*IT_0427;
    const ccomplex_t IT_0429 = IT_0001*IT_0196*IT_0327;
    const ccomplex_t IT_0430 = IT_0008*IT_0185*IT_0196;
    const ccomplex_t IT_0431 = -IT_0430;
    const ccomplex_t IT_0432 = IT_0007*IT_0185*IT_0327;
    const ccomplex_t IT_0433 = -IT_0432;
    const ccomplex_t IT_0434 = s_12*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0435 = -IT_0434;
    const ccomplex_t IT_0436 = IT_0409 + IT_0411 + IT_0429 + IT_0431 + IT_0433
       + IT_0435;
    const ccomplex_t IT_0437 = s_14*s_24*IT_0001*IT_0327;
    const ccomplex_t IT_0438 = s_13*s_23*IT_0001*IT_0329;
    const ccomplex_t IT_0439 = s_12*IT_0008;
    const ccomplex_t IT_0440 = s_14*s_24*IT_0008*IT_0185;
    const ccomplex_t IT_0441 = -IT_0440;
    const ccomplex_t IT_0442 = s_13*s_23*IT_0008*IT_0185;
    const ccomplex_t IT_0443 = -IT_0442;
    const ccomplex_t IT_0444 = s_12*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0445 = -IT_0444;
    const ccomplex_t IT_0446 = IT_0437 + IT_0438 + IT_0439 + IT_0441 + IT_0443
       + IT_0445;
    const ccomplex_t IT_0447 = 2*IT_0326;
    const ccomplex_t IT_0448 = conj(IT_0421)*IT_0447;
    const ccomplex_t IT_0449 = 2*conj(IT_0326);
    const ccomplex_t IT_0450 = IT_0421*IT_0449;
    const ccomplex_t IT_0451 = IT_0037*IT_0383;
    const ccomplex_t IT_0452 = IT_0010*IT_0451;
    const ccomplex_t IT_0453 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0454 = IT_0013*IT_0453;
    const ccomplex_t IT_0455 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0456 = IT_0013*IT_0455;
    const ccomplex_t IT_0457 = 1.4142135623731*IT_0456;
    const ccomplex_t IT_0458 = (0 + _Complex_I*1)*(IT_0454 + 0.5*IT_0457);
    const ccomplex_t IT_0459 = -IT_0458;
    const ccomplex_t IT_0460 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0461 = IT_0013*IT_0460;
    const ccomplex_t IT_0462 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0463 = IT_0013*IT_0462;
    const ccomplex_t IT_0464 = 1.4142135623731*IT_0463;
    const ccomplex_t IT_0465 = (0 + _Complex_I*1)*(IT_0461 + 0.5*IT_0464);
    const ccomplex_t IT_0466 = -IT_0465;
    const ccomplex_t IT_0467 = IT_0459*IT_0466;
    const ccomplex_t IT_0468 = IT_0305*IT_0467;
    const ccomplex_t IT_0469 = (0 + _Complex_I*1)*IT_0468;
    const ccomplex_t IT_0470 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0471 = IT_0013*IT_0470;
    const ccomplex_t IT_0472 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0473 = IT_0013*IT_0472;
    const ccomplex_t IT_0474 = 1.4142135623731*IT_0473;
    const ccomplex_t IT_0475 = (0 + _Complex_I*1)*(IT_0471 + 0.5*IT_0474);
    const ccomplex_t IT_0476 = -IT_0475;
    const ccomplex_t IT_0477 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0478 = IT_0013*IT_0477;
    const ccomplex_t IT_0479 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0480 = IT_0013*IT_0479;
    const ccomplex_t IT_0481 = 1.4142135623731*IT_0480;
    const ccomplex_t IT_0482 = (0 + _Complex_I*1)*(IT_0478 + 0.5*IT_0481);
    const ccomplex_t IT_0483 = -IT_0482;
    const ccomplex_t IT_0484 = IT_0476*IT_0483;
    const ccomplex_t IT_0485 = IT_0308*IT_0484;
    const ccomplex_t IT_0486 = (0 + _Complex_I*1)*IT_0485;
    const ccomplex_t IT_0487 = -IT_0452 + 2*IT_0469 + 2*IT_0486;
    const ccomplex_t IT_0488 = IT_0395*conj(IT_0487);
    const ccomplex_t IT_0489 = IT_0397*IT_0487;
    const ccomplex_t IT_0490 = s_14*s_34*IT_0185;
    const ccomplex_t IT_0491 = s_13*s_23;
    const ccomplex_t IT_0492 = (-2)*IT_0491;
    const ccomplex_t IT_0493 = IT_0214 + IT_0492;
    const ccomplex_t IT_0494 = IT_0185*IT_0493;
    const ccomplex_t IT_0495 = (s_13 + -IT_0490)*(s_12 + (-0.5)*IT_0494);
    const ccomplex_t IT_0496 = (-2)*IT_0495;
    const ccomplex_t IT_0497 = s_24*s_34*IT_0185;
    const ccomplex_t IT_0498 = (s_12 + (-0.5)*IT_0494)*(s_23 + -IT_0497);
    const ccomplex_t IT_0499 = (-2)*IT_0498;
    const ccomplex_t IT_0500 = IT_0014*IT_0037;
    const ccomplex_t IT_0501 = IT_0010*IT_0500;
    const ccomplex_t IT_0502 = IT_0395*conj(IT_0501);
    const ccomplex_t IT_0503 = IT_0397*IT_0501;
    const ccomplex_t IT_0504 = IT_0014*IT_0026;
    const ccomplex_t IT_0505 = IT_0010*IT_0504;
    const ccomplex_t IT_0506 = IT_0425*conj(IT_0505);
    const ccomplex_t IT_0507 = IT_0427*IT_0505;
    const ccomplex_t IT_0508 = 2*IT_0352;
    const ccomplex_t IT_0509 = IT_0354 + IT_0363 + IT_0371 + IT_0373 + IT_0508;
    const ccomplex_t IT_0510 = 2*IT_0357;
    const ccomplex_t IT_0511 = IT_0358 + IT_0367 + IT_0377 + IT_0379 + IT_0510;
    const ccomplex_t IT_0512 = IT_0037*IT_0380;
    const ccomplex_t IT_0513 = IT_0010*IT_0512;
    const ccomplex_t IT_0514 = IT_0395*conj(IT_0513);
    const ccomplex_t IT_0515 = IT_0397*IT_0513;
    const ccomplex_t IT_0516 = conj(IT_0382)*IT_0425;
    const ccomplex_t IT_0517 = IT_0382*IT_0427;
    const ccomplex_t IT_0518 = s_13*IT_0161;
    const ccomplex_t IT_0519 = s_14*s_34;
    const ccomplex_t IT_0520 = -IT_0519;
    const ccomplex_t IT_0521 = IT_0518 + IT_0520;
    const ccomplex_t IT_0522 = IT_0185*IT_0521;
    const ccomplex_t IT_0523 = -IT_0490 + (-0.25)*IT_0522;
    const ccomplex_t IT_0524 = s_13 + IT_0523;
    const ccomplex_t IT_0525 = IT_0000*IT_0524;
    const ccomplex_t IT_0526 = (-4)*IT_0525;
    const ccomplex_t IT_0527 = s_23*IT_0161;
    const ccomplex_t IT_0528 = s_24*s_34;
    const ccomplex_t IT_0529 = -IT_0528;
    const ccomplex_t IT_0530 = IT_0527 + IT_0529;
    const ccomplex_t IT_0531 = IT_0185*IT_0530;
    const ccomplex_t IT_0532 = -IT_0497 + (-0.25)*IT_0531;
    const ccomplex_t IT_0533 = s_23 + IT_0532;
    const ccomplex_t IT_0534 = IT_0000*IT_0533;
    const ccomplex_t IT_0535 = (-4)*IT_0534;
    const ccomplex_t IT_0536 = -IT_0307 + -IT_0325;
    const ccomplex_t IT_0537 = IT_0395*conj(IT_0536);
    const ccomplex_t IT_0538 = conj(IT_0164)*IT_0395;
    const ccomplex_t IT_0539 = conj(IT_0182)*IT_0395;
    const ccomplex_t IT_0540 = IT_0397*IT_0536;
    const ccomplex_t IT_0541 = IT_0164*IT_0397;
    const ccomplex_t IT_0542 = IT_0182*IT_0397;
    const ccomplex_t IT_0543 = IT_0272 + IT_0277;
    const ccomplex_t IT_0544 = IT_0425*conj(IT_0543);
    const ccomplex_t IT_0545 = -IT_0469 + -IT_0486;
    const ccomplex_t IT_0546 = IT_0425*conj(IT_0545);
    const ccomplex_t IT_0547 = IT_0427*IT_0543;
    const ccomplex_t IT_0548 = IT_0427*IT_0545;
    const ccomplex_t IT_0549 = s_34*m_N_2*m_N_3*IT_0185;
    const ccomplex_t IT_0550 = s_13*s_34*IT_0185;
    const ccomplex_t IT_0551 = -IT_0550;
    const ccomplex_t IT_0552 = s_14 + IT_0551;
    const ccomplex_t IT_0553 = IT_0549*IT_0552;
    const ccomplex_t IT_0554 = s_23*s_34*IT_0185;
    const ccomplex_t IT_0555 = -IT_0554;
    const ccomplex_t IT_0556 = s_24 + IT_0555;
    const ccomplex_t IT_0557 = IT_0549*IT_0556;
    const ccomplex_t IT_0558 = IT_0387*conj(IT_0513);
    const ccomplex_t IT_0559 = IT_0389*IT_0513;
    const ccomplex_t IT_0560 = 2*IT_0501;
    const ccomplex_t IT_0561 = conj(IT_0386)*IT_0560;
    const ccomplex_t IT_0562 = 2*conj(IT_0501);
    const ccomplex_t IT_0563 = IT_0386*IT_0562;
    const ccomplex_t IT_0564 = s_24*s_34*IT_0007*IT_0185;
    const ccomplex_t IT_0565 = -IT_0564;
    const ccomplex_t IT_0566 = IT_0352 + IT_0353 + IT_0354 + IT_0363 + IT_0369
       + IT_0565;
    const ccomplex_t IT_0567 = s_14*s_34*IT_0008*IT_0185;
    const ccomplex_t IT_0568 = -IT_0567;
    const ccomplex_t IT_0569 = IT_0356 + IT_0357 + IT_0358 + IT_0367 + IT_0375
       + IT_0568;
    const ccomplex_t IT_0570 = 2*IT_0505;
    const ccomplex_t IT_0571 = conj(IT_0386)*IT_0570;
    const ccomplex_t IT_0572 = 2*conj(IT_0505);
    const ccomplex_t IT_0573 = IT_0386*IT_0572;
    const ccomplex_t IT_0574 = -IT_0007;
    const ccomplex_t IT_0575 = IT_0197 + IT_0574;
    const ccomplex_t IT_0576 = IT_0000*IT_0575;
    const ccomplex_t IT_0577 = -IT_0008;
    const ccomplex_t IT_0578 = IT_0330 + IT_0577;
    const ccomplex_t IT_0579 = IT_0000*IT_0578;
    const ccomplex_t IT_0580 = IT_0394*conj(IT_0421);
    const ccomplex_t IT_0581 = 6*IT_0580;
    const ccomplex_t IT_0582 = conj(IT_0394)*IT_0421;
    const ccomplex_t IT_0583 = 6*IT_0582;
    const ccomplex_t IT_0584 = 1.125*IT_0000;
    const ccomplex_t IT_0585 = IT_0307 + IT_0325;
    const ccomplex_t IT_0586 = m_N_3*IT_0585;
    const ccomplex_t IT_0587 = -IT_0586;
    const ccomplex_t IT_0588 = m_N_3*IT_0486;
    const ccomplex_t IT_0589 = m_N_3*IT_0469;
    const ccomplex_t IT_0590 = IT_0303*IT_0459;
    const ccomplex_t IT_0591 = IT_0227*IT_0305*IT_0590;
    const ccomplex_t IT_0592 = IT_0315*IT_0476;
    const ccomplex_t IT_0593 = IT_0237*IT_0308*IT_0592;
    const ccomplex_t IT_0594 = IT_0322*IT_0483;
    const ccomplex_t IT_0595 = IT_0237*IT_0308*IT_0594;
    const ccomplex_t IT_0596 = -IT_0595;
    const ccomplex_t IT_0597 = IT_0296*IT_0466;
    const ccomplex_t IT_0598 = IT_0227*IT_0305*IT_0597;
    const ccomplex_t IT_0599 = -IT_0598;
    const ccomplex_t IT_0600 = 16*conj(IT_0587);
    const ccomplex_t IT_0601 = conj(IT_0591) + conj(IT_0593);
    const ccomplex_t IT_0602 = 16*IT_0596;
    const ccomplex_t IT_0603 = 16*IT_0599;
    const ccomplex_t IT_0604 = IT_0591 + IT_0593;
    const ccomplex_t IT_0605 = 16*conj(IT_0596);
    const ccomplex_t IT_0606 = 16*conj(IT_0599);
    const ccomplex_t IT_0607 = (-16)*conj(IT_0587);
    const ccomplex_t IT_0608 = conj(IT_0596) + conj(IT_0599);
    const ccomplex_t IT_0609 = (-16)*IT_0589;
    const ccomplex_t IT_0610 = IT_0596 + IT_0599;
    const ccomplex_t IT_0611 = (-16)*conj(IT_0588);
    const ccomplex_t IT_0612 = (-16)*conj(IT_0589);
    const ccomplex_t IT_0613 = 16*IT_0220;
    const ccomplex_t IT_0614 = (-16)*IT_0247;
    const ccomplex_t IT_0615 = 16*conj(IT_0220);
    const ccomplex_t IT_0616 = (-16)*conj(IT_0247);
    const ccomplex_t IT_0617 = IT_0587*(16*conj(IT_0588) + 16*conj(IT_0589) + 
      (-16)*conj(IT_0591) + (-16)*conj(IT_0593)) + IT_0588*((-16)*conj(IT_0596) 
      + (-16)*conj(IT_0599) + IT_0600) + IT_0589*IT_0600 + IT_0601*(IT_0602 +
       IT_0603) + IT_0604*(IT_0605 + IT_0606 + IT_0607) + IT_0608*IT_0609 +
       IT_0610*(IT_0611 + IT_0612) + IT_0282*(IT_0613 + IT_0614) + IT_0286*
      (IT_0615 + IT_0616);
    const ccomplex_t IT_0618 = pow(m_W, 4);
    const ccomplex_t IT_0619 = s_12*IT_0618;
    const ccomplex_t IT_0620 = s_13*s_23*IT_0161;
    const ccomplex_t IT_0621 = s_14*s_24*IT_0161;
    const ccomplex_t IT_0622 = s_14*s_23*s_34;
    const ccomplex_t IT_0623 = s_13*s_24*s_34;
    const ccomplex_t IT_0624 = s_12*IT_0002;
    const ccomplex_t IT_0625 = (-2)*IT_0620 + (-2)*IT_0621 + 2*IT_0622 + 2
      *IT_0623 + (-2)*IT_0624;
    const ccomplex_t IT_0626 = IT_0619 + IT_0625;
    const ccomplex_t IT_0627 = IT_0001*IT_0626;
    const ccomplex_t IT_0628 = (-0.25)*IT_0627;
    const ccomplex_t IT_0629 = -IT_0187;
    const ccomplex_t IT_0630 = -IT_0186;
    const ccomplex_t IT_0631 = 8*conj(IT_0278);
    const ccomplex_t IT_0632 = 0.125*IT_0631;
    const ccomplex_t IT_0633 = 8*conj(IT_0274);
    const ccomplex_t IT_0634 = 0.125*IT_0633;
    const ccomplex_t IT_0635 = conj(IT_0281) + IT_0632 + IT_0634;
    const ccomplex_t IT_0636 = 8*IT_0247;
    const ccomplex_t IT_0637 = (-8)*IT_0220;
    const ccomplex_t IT_0638 = IT_0636 + IT_0637;
    const ccomplex_t IT_0639 = 8*conj(IT_0247);
    const ccomplex_t IT_0640 = (-8)*conj(IT_0220);
    const ccomplex_t IT_0641 = IT_0639 + IT_0640;
    const ccomplex_t IT_0642 = 8*IT_0588;
    const ccomplex_t IT_0643 = 8*conj(IT_0588);
    const ccomplex_t IT_0644 = conj(IT_0278) + conj(IT_0281);
    const ccomplex_t IT_0645 = 8*IT_0599;
    const ccomplex_t IT_0646 = 8*conj(IT_0599);
    const ccomplex_t IT_0647 = 8*IT_0274;
    const ccomplex_t IT_0648 = 8*IT_0278;
    const ccomplex_t IT_0649 = (-8)*IT_0587;
    const ccomplex_t IT_0650 = (-8)*conj(IT_0587);
    const ccomplex_t IT_0651 = 8*IT_0220;
    const ccomplex_t IT_0652 = (-8)*IT_0247;
    const ccomplex_t IT_0653 = 8*conj(IT_0220);
    const ccomplex_t IT_0654 = (-8)*conj(IT_0247);
    const ccomplex_t IT_0655 = 0.125*IT_0599*IT_0633 + IT_0596*IT_0635 + 0.125
      *conj(IT_0588)*IT_0637 + 0.125*conj(IT_0589)*IT_0638 + 0.125*IT_0588
      *IT_0640 + 0.125*IT_0589*IT_0641 + 0.125*conj(IT_0247)*IT_0642 + 0.125
      *IT_0247*IT_0643 + 0.125*IT_0644*IT_0645 + 0.125*(IT_0278 + IT_0281)
      *IT_0646 + 0.125*conj(IT_0599)*IT_0647 + conj(IT_0596)*(IT_0281 + 0.125
      *IT_0647 + 0.125*IT_0648) + 0.125*IT_0282*IT_0649 + 0.125*IT_0286*IT_0650 
      + 0.125*IT_0601*(IT_0651 + IT_0652) + 0.125*IT_0604*(IT_0653 + IT_0654);
    const ccomplex_t IT_0656 = (-2)*IT_0002;
    const ccomplex_t IT_0657 = IT_0618 + IT_0656;
    const ccomplex_t IT_0658 = IT_0001*IT_0657;
    const ccomplex_t IT_0659 = 0.25*IT_0658;
    const ccomplex_t IT_0660 = 1 + IT_0659;
    const ccomplex_t IT_0661 = IT_0000*IT_0660;
    const ccomplex_t IT_0662 = 8*IT_0589;
    const ccomplex_t IT_0663 = 8*conj(IT_0589);
    const ccomplex_t IT_0664 = (-8)*IT_0274;
    const ccomplex_t IT_0665 = (-8)*IT_0278;
    const ccomplex_t IT_0666 = (-8)*IT_0281;
    const ccomplex_t IT_0667 = (-8)*conj(IT_0274);
    const ccomplex_t IT_0668 = (-8)*conj(IT_0278);
    const ccomplex_t IT_0669 = (-8)*conj(IT_0281);
    const ccomplex_t IT_0670 = IT_0608*IT_0638 + IT_0610*IT_0641 + conj
      (IT_0247)*IT_0649 + IT_0247*IT_0650 + conj(IT_0587)*IT_0651 + IT_0587
      *IT_0653 + IT_0282*(IT_0642 + IT_0662) + IT_0286*(IT_0643 + IT_0663) +
       IT_0601*(IT_0664 + IT_0665 + IT_0666) + IT_0604*(IT_0667 + IT_0668 +
       IT_0669);
    const ccomplex_t IT_0671 = s_34*m_N_2*IT_0185;
    const ccomplex_t IT_0672 = IT_0556*IT_0671;
    const ccomplex_t IT_0673 = 2*IT_0093;
    const ccomplex_t IT_0674 = (-2)*IT_0096;
    const ccomplex_t IT_0675 = IT_0673 + IT_0674;
    const ccomplex_t IT_0676 = s_34*m_N_3*IT_0185;
    const ccomplex_t IT_0677 = IT_0552*IT_0676;
    const ccomplex_t IT_0678 = 2*IT_0040;
    const ccomplex_t IT_0679 = (-2)*IT_0139;
    const ccomplex_t IT_0680 = IT_0678 + IT_0679;
    const ccomplex_t IT_0681 = m_N_2*IT_0533;
    const ccomplex_t IT_0682 = (-4)*IT_0681;
    const ccomplex_t IT_0683 = 2*IT_0220;
    const ccomplex_t IT_0684 = (-2)*IT_0272 + (-2)*IT_0277 + -IT_0452;
    const ccomplex_t IT_0685 = 2*IT_0684;
    const ccomplex_t IT_0686 = 2*IT_0487;
    const ccomplex_t IT_0687 = IT_0002*IT_0185;
    const ccomplex_t IT_0688 = (s_12 + (-0.5)*IT_0494)*(IT_0161 + -IT_0687);
    const ccomplex_t IT_0689 = (-2)*IT_0688;
    const ccomplex_t IT_0690 = -IT_0687;
    const ccomplex_t IT_0691 = IT_0161 + IT_0690;
    const ccomplex_t IT_0692 = IT_0000*IT_0691;
    const ccomplex_t IT_0693 = (-6)*IT_0501;
    const ccomplex_t IT_0694 = 2*conj(IT_0220);
    const ccomplex_t IT_0695 = 2*conj(IT_0684);
    const ccomplex_t IT_0696 = 2*conj(IT_0487);
    const ccomplex_t IT_0697 = (-6)*conj(IT_0501);
    const ccomplex_t IT_0698 = -IT_0490;
    const ccomplex_t IT_0699 = s_13 + IT_0698;
    const ccomplex_t IT_0700 = IT_0000*IT_0699;
    const ccomplex_t IT_0701 = -IT_0497;
    const ccomplex_t IT_0702 = s_23 + IT_0701;
    const ccomplex_t IT_0703 = IT_0000*IT_0702;
    const ccomplex_t IT_0704 = IT_0394*conj(IT_0505);
    const ccomplex_t IT_0705 = (-6)*IT_0704;
    const ccomplex_t IT_0706 = conj(IT_0394)*IT_0505;
    const ccomplex_t IT_0707 = (-6)*IT_0706;
    const ccomplex_t IT_0708 = conj(IT_0421)*IT_0693;
    const ccomplex_t IT_0709 = IT_0421*IT_0697;
    const ccomplex_t IT_0710 = s_12*IT_0005;
    const ccomplex_t IT_0711 = (-4)*conj(IT_0139);
    const ccomplex_t IT_0712 = 0.25*IT_0711;
    const ccomplex_t IT_0713 = conj(IT_0040) + IT_0712;
    const ccomplex_t IT_0714 = (-4)*IT_0139;
    const ccomplex_t IT_0715 = IT_0008*IT_0161;
    const ccomplex_t IT_0716 = s_23*s_24*s_34*IT_0185;
    const ccomplex_t IT_0717 = 0.5*IT_0329 + (-0.5)*IT_0715 + -IT_0716;
    const ccomplex_t IT_0718 = IT_0327 + IT_0717;
    const ccomplex_t IT_0719 = IT_0185*IT_0718;
    const ccomplex_t IT_0720 = -IT_0330;
    const ccomplex_t IT_0721 = IT_0008 + IT_0719 + IT_0720;
    const ccomplex_t IT_0722 = m_N_2*IT_0721;
    const ccomplex_t IT_0723 = (-4)*IT_0394;
    const ccomplex_t IT_0724 = (-4)*conj(IT_0394);
    const ccomplex_t IT_0725 = conj(IT_0588) + conj(IT_0589);
    const ccomplex_t IT_0726 = (-4)*IT_0421;
    const ccomplex_t IT_0727 = IT_0588 + IT_0589;
    const ccomplex_t IT_0728 = (-4)*conj(IT_0421);
    const ccomplex_t IT_0729 = 2*IT_0096;
    const ccomplex_t IT_0730 = (-2)*IT_0093;
    const ccomplex_t IT_0731 = IT_0729 + IT_0730;
    const ccomplex_t IT_0732 = 2*IT_0139;
    const ccomplex_t IT_0733 = (-2)*IT_0040;
    const ccomplex_t IT_0734 = IT_0732 + IT_0733;
    const ccomplex_t IT_0735 = IT_0518 + IT_0519;
    const ccomplex_t IT_0736 = s_13*IT_0002*IT_0185;
    const ccomplex_t IT_0737 = (-2)*IT_0736;
    const ccomplex_t IT_0738 = IT_0735 + IT_0737;
    const ccomplex_t IT_0739 = IT_0185*IT_0738;
    const ccomplex_t IT_0740 = -IT_0490 + 0.5*IT_0739;
    const ccomplex_t IT_0741 = s_13 + IT_0740;
    const ccomplex_t IT_0742 = m_N_3*IT_0741;
    const ccomplex_t IT_0743 = 2*IT_0742;
    const ccomplex_t IT_0744 = (-2)*IT_0220;
    const ccomplex_t IT_0745 = (-2)*conj(IT_0220);
    const ccomplex_t IT_0746 = -(IT_0007 + -IT_0195)*(s_12 + (-0.5)*IT_0348);
    const ccomplex_t IT_0747 = (-2)*IT_0746;
    const ccomplex_t IT_0748 = (-2)*IT_0183;
    const ccomplex_t IT_0749 = 2*IT_0329 + -IT_0715 + (-2)*IT_0716;
    const ccomplex_t IT_0750 = IT_0327 + IT_0749;
    const ccomplex_t IT_0751 = IT_0185*IT_0750;
    const ccomplex_t IT_0752 = -IT_0008 + (-0.5)*IT_0751;
    const ccomplex_t IT_0753 = IT_0328 + IT_0752;
    const ccomplex_t IT_0754 = m_N_2*IT_0753;
    const ccomplex_t IT_0755 = (-2)*IT_0754;
    const ccomplex_t IT_0756 = (-0.5)*IT_0494;
    const ccomplex_t IT_0757 = (-0.5)*IT_0348;
    const ccomplex_t IT_0758 = (-2)*IT_0620 + (-2)*IT_0621 + 4*IT_0623;
    const ccomplex_t IT_0759 = IT_0619 + IT_0758;
    const ccomplex_t IT_0760 = IT_0001*IT_0759;
    const ccomplex_t IT_0761 = 0.25*IT_0760;
    const ccomplex_t IT_0762 = s_12 + IT_0756 + IT_0757 + IT_0761;
    const ccomplex_t IT_0763 = s_13*s_24*s_34*IT_0185;
    const ccomplex_t IT_0764 = (-2)*IT_0345 + -IT_0491 + 2*IT_0763;
    const ccomplex_t IT_0765 = IT_0214 + IT_0764;
    const ccomplex_t IT_0766 = IT_0185*IT_0765;
    const ccomplex_t IT_0767 = -IT_0186 + (-0.5)*IT_0766;
    const ccomplex_t IT_0768 = s_12 + IT_0767;
    const ccomplex_t IT_0769 = m_N_2*IT_0768;
    const ccomplex_t IT_0770 = (-2)*IT_0769;
    const ccomplex_t IT_0771 = (-2)*conj(IT_0183);
    const ccomplex_t IT_0772 = -IT_0491;
    const ccomplex_t IT_0773 = IT_0214 + IT_0772;
    const ccomplex_t IT_0774 = IT_0185*IT_0773;
    const ccomplex_t IT_0775 = -IT_0186 + (-0.25)*IT_0774;
    const ccomplex_t IT_0776 = s_12 + IT_0775;
    const ccomplex_t IT_0777 = m_N_2*IT_0776;
    const ccomplex_t IT_0778 = (-4)*IT_0777;
    const ccomplex_t IT_0779 = IT_0387 + IT_0422;
    const ccomplex_t IT_0780 = IT_0351 + IT_0389;
    const ccomplex_t IT_0781 = (-2)*IT_0386;
    const ccomplex_t IT_0782 = IT_0748 + IT_0781;
    const ccomplex_t IT_0783 = (-2)*conj(IT_0386);
    const ccomplex_t IT_0784 = IT_0771 + IT_0783;
    const ccomplex_t IT_0785 = 2*IT_0587;
    const ccomplex_t IT_0786 = 2*conj(IT_0587);
    const ccomplex_t IT_0787 = (-2)*IT_0684;
    const ccomplex_t IT_0788 = (-2)*conj(IT_0684);
    const ccomplex_t IT_0789 = s_13*s_14*s_34*IT_0185;
    const ccomplex_t IT_0790 = IT_0007*IT_0161;
    const ccomplex_t IT_0791 = 2*IT_0196 + (-2)*IT_0789 + -IT_0790;
    const ccomplex_t IT_0792 = IT_0194 + IT_0791;
    const ccomplex_t IT_0793 = IT_0185*IT_0792;
    const ccomplex_t IT_0794 = -IT_0007 + (-0.5)*IT_0793;
    const ccomplex_t IT_0795 = IT_0195 + IT_0794;
    const ccomplex_t IT_0796 = m_N_3*IT_0795;
    const ccomplex_t IT_0797 = (-2)*IT_0796;
    const ccomplex_t IT_0798 = -(IT_0008 + -IT_0328)*(s_12 + (-0.5)*IT_0348);
    const ccomplex_t IT_0799 = (-2)*IT_0798;
    const ccomplex_t IT_0800 = (-2)*IT_0326;
    const ccomplex_t IT_0801 = -IT_0345;
    const ccomplex_t IT_0802 = IT_0214 + IT_0801;
    const ccomplex_t IT_0803 = IT_0185*IT_0802;
    const ccomplex_t IT_0804 = -IT_0187 + (-0.25)*IT_0803;
    const ccomplex_t IT_0805 = s_12 + IT_0804;
    const ccomplex_t IT_0806 = m_N_3*IT_0805;
    const ccomplex_t IT_0807 = (-4)*IT_0806;
    const ccomplex_t IT_0808 = (-2)*IT_0394;
    const ccomplex_t IT_0809 = (-2)*conj(IT_0394);
    const ccomplex_t IT_0810 = IT_0000*IT_0189;
    const ccomplex_t IT_0811 = IT_0000*IT_0340;
    const ccomplex_t IT_0812 = conj(IT_0386)*IT_0808;
    const ccomplex_t IT_0813 = IT_0386*IT_0809;
    const ccomplex_t IT_0814 = s_14*s_23*s_34*IT_0185;
    const ccomplex_t IT_0815 = (-2)*IT_0345 + -IT_0491 + 2*IT_0814;
    const ccomplex_t IT_0816 = IT_0214 + IT_0815;
    const ccomplex_t IT_0817 = IT_0185*IT_0816;
    const ccomplex_t IT_0818 = -IT_0186 + (-0.5)*IT_0817;
    const ccomplex_t IT_0819 = s_12 + IT_0818;
    const ccomplex_t IT_0820 = m_N_3*IT_0819;
    const ccomplex_t IT_0821 = (-2)*IT_0820;
    const ccomplex_t IT_0822 = IT_0781 + IT_0800;
    const ccomplex_t IT_0823 = (-2)*conj(IT_0326);
    const ccomplex_t IT_0824 = IT_0783 + IT_0823;
    const ccomplex_t IT_0825 = (-2)*IT_0487;
    const ccomplex_t IT_0826 = (-2)*conj(IT_0487);
    const ccomplex_t IT_0827 = -(IT_0007 + -IT_0197)*(s_12 + (-0.5)*IT_0494);
    const ccomplex_t IT_0828 = (-2)*IT_0827;
    const ccomplex_t IT_0829 = -(IT_0008 + -IT_0330)*(s_12 + (-0.5)*IT_0494);
    const ccomplex_t IT_0830 = (-2)*IT_0829;
    const ccomplex_t IT_0831 = conj(IT_0394)*IT_0808;
    const ccomplex_t IT_0832 = (-2)*IT_0421;
    const ccomplex_t IT_0833 = conj(IT_0421)*IT_0832;
    const ccomplex_t IT_0834 = -IT_0790;
    const ccomplex_t IT_0835 = IT_0196 + IT_0834;
    const ccomplex_t IT_0836 = IT_0185*IT_0835;
    const ccomplex_t IT_0837 = (-4)*IT_0575;
    const ccomplex_t IT_0838 = IT_0836 + IT_0837;
    const ccomplex_t IT_0839 = m_N_3*IT_0838;
    const ccomplex_t IT_0840 = (-2)*conj(IT_0421);
    const ccomplex_t IT_0841 = -IT_0345 + (-2)*IT_0491 + 2*IT_0814;
    const ccomplex_t IT_0842 = IT_0214 + IT_0841;
    const ccomplex_t IT_0843 = IT_0185*IT_0842;
    const ccomplex_t IT_0844 = -IT_0187 + (-0.5)*IT_0843;
    const ccomplex_t IT_0845 = s_12 + IT_0844;
    const ccomplex_t IT_0846 = m_N_2*IT_0845;
    const ccomplex_t IT_0847 = (-2)*IT_0846;
    const ccomplex_t IT_0848 = -IT_0715;
    const ccomplex_t IT_0849 = IT_0327 + IT_0848;
    const ccomplex_t IT_0850 = IT_0185*IT_0849;
    const ccomplex_t IT_0851 = IT_0328 + IT_0577;
    const ccomplex_t IT_0852 = (-4)*IT_0851;
    const ccomplex_t IT_0853 = IT_0850 + IT_0852;
    const ccomplex_t IT_0854 = m_N_2*IT_0853;
    const ccomplex_t IT_0855 = IT_0387 + IT_0447;
    const ccomplex_t IT_0856 = IT_0389 + IT_0449;
    const ccomplex_t IT_0857 = (-2)*IT_0587;
    const ccomplex_t IT_0858 = (-2)*conj(IT_0587);
    const ccomplex_t IT_0859 = m_N_3*IT_0776;
    const ccomplex_t IT_0860 = (-4)*IT_0859;
    const ccomplex_t IT_0861 = IT_0194 + IT_0834;
    const ccomplex_t IT_0862 = IT_0185*IT_0861;
    const ccomplex_t IT_0863 = IT_0195 + IT_0574;
    const ccomplex_t IT_0864 = (-4)*IT_0863;
    const ccomplex_t IT_0865 = IT_0862 + IT_0864;
    const ccomplex_t IT_0866 = m_N_3*IT_0865;
    const ccomplex_t IT_0867 = IT_0329 + IT_0848;
    const ccomplex_t IT_0868 = IT_0185*IT_0867;
    const ccomplex_t IT_0869 = (-4)*IT_0578;
    const ccomplex_t IT_0870 = IT_0868 + IT_0869;
    const ccomplex_t IT_0871 = m_N_2*IT_0870;
    const ccomplex_t IT_0872 = IT_0185*IT_0527;
    const ccomplex_t IT_0873 = s_24*IT_0161;
    const ccomplex_t IT_0874 = IT_0206*IT_0873;
    const ccomplex_t IT_0875 = (-0.5)*IT_0872 + 0.25*IT_0874;
    const ccomplex_t IT_0876 = s_23 + IT_0875;
    const ccomplex_t IT_0877 = m_N_2*IT_0876;
    const ccomplex_t IT_0878 = 4*IT_0877;
    const ccomplex_t IT_0879 = 2*conj(IT_0513);
    const ccomplex_t IT_0880 = IT_0185*IT_0518;
    const ccomplex_t IT_0881 = s_14*IT_0161;
    const ccomplex_t IT_0882 = s_13*s_34;
    const ccomplex_t IT_0883 = (-2)*IT_0882;
    const ccomplex_t IT_0884 = IT_0881 + IT_0883;
    const ccomplex_t IT_0885 = IT_0206*IT_0884;
    const ccomplex_t IT_0886 = (-0.5)*IT_0880 + (-0.25)*IT_0885;
    const ccomplex_t IT_0887 = s_13 + IT_0886;
    const ccomplex_t IT_0888 = m_N_3*IT_0887;
    const ccomplex_t IT_0889 = (-4)*IT_0888;
    const ccomplex_t IT_0890 = 2*conj(IT_0382);
    const ccomplex_t IT_0891 = m_N_2*IT_0805;
    const ccomplex_t IT_0892 = (-4)*IT_0891;
    const ccomplex_t IT_0893 = -IT_0345 + (-2)*IT_0491 + 2*IT_0763;
    const ccomplex_t IT_0894 = IT_0214 + IT_0893;
    const ccomplex_t IT_0895 = IT_0185*IT_0894;
    const ccomplex_t IT_0896 = -IT_0187 + (-0.5)*IT_0895;
    const ccomplex_t IT_0897 = s_12 + IT_0896;
    const ccomplex_t IT_0898 = m_N_3*IT_0897;
    const ccomplex_t IT_0899 = (-2)*IT_0898;
    const ccomplex_t IT_0900 = (-2)*conj(IT_0501);
    const ccomplex_t IT_0901 = 0.125*IT_0572*IT_0743 + IT_0635*IT_0762 + 0.125
      *IT_0695*IT_0770 + 0.125*IT_0755*IT_0826 + 0.125*IT_0856*IT_0860 + 0.125
      *IT_0784*IT_0866 + 0.125*IT_0397*IT_0871 + 0.125*IT_0878*IT_0879 + 0.125
      *IT_0889*IT_0890 + 0.125*IT_0809*IT_0892 + 0.125*IT_0840*IT_0899 + 0.125
      *IT_0682*IT_0900;
    const ccomplex_t IT_0902 = m_N_3*IT_0524;
    const ccomplex_t IT_0903 = (-4)*IT_0902;
    const ccomplex_t IT_0904 = (-2)*IT_0501;
    const ccomplex_t IT_0905 = 2*IT_0513;
    const ccomplex_t IT_0906 = 2*IT_0382;
    const ccomplex_t IT_0907 = (-2)*IT_0505;
    const ccomplex_t IT_0908 = IT_0560*IT_0743 + IT_0638*IT_0762 + IT_0770
      *IT_0779 + IT_0755*IT_0822 + IT_0686*IT_0860 + IT_0787*IT_0866 + IT_0425
      *IT_0871 + IT_0832*IT_0892 + IT_0808*IT_0899 + IT_0889*IT_0905 + IT_0878
      *IT_0906 + IT_0682*IT_0907;
    const ccomplex_t IT_0909 = (-2)*conj(IT_0505);
    const ccomplex_t IT_0910 = IT_0562*IT_0743 + IT_0640*IT_0762 + IT_0770
      *IT_0780 + IT_0755*IT_0824 + IT_0696*IT_0860 + IT_0788*IT_0866 + IT_0427
      *IT_0871 + IT_0879*IT_0889 + IT_0878*IT_0890 + IT_0840*IT_0892 + IT_0809
      *IT_0899 + IT_0682*IT_0909;
    const ccomplex_t IT_0911 = IT_0527 + IT_0528;
    const ccomplex_t IT_0912 = s_23*IT_0002*IT_0185;
    const ccomplex_t IT_0913 = (-2)*IT_0912;
    const ccomplex_t IT_0914 = IT_0911 + IT_0913;
    const ccomplex_t IT_0915 = IT_0185*IT_0914;
    const ccomplex_t IT_0916 = -IT_0497 + 0.5*IT_0915;
    const ccomplex_t IT_0917 = s_23 + IT_0916;
    const ccomplex_t IT_0918 = m_N_2*IT_0917;
    const ccomplex_t IT_0919 = 2*IT_0918;
    const ccomplex_t IT_0920 = IT_0206*IT_0881;
    const ccomplex_t IT_0921 = (-0.5)*IT_0880 + 0.25*IT_0920;
    const ccomplex_t IT_0922 = s_13 + IT_0921;
    const ccomplex_t IT_0923 = m_N_3*IT_0922;
    const ccomplex_t IT_0924 = 4*IT_0923;
    const ccomplex_t IT_0925 = (-2)*IT_0513;
    const ccomplex_t IT_0926 = (-2)*conj(IT_0513);
    const ccomplex_t IT_0927 = (-2)*IT_0382;
    const ccomplex_t IT_0928 = (-2)*conj(IT_0382);
    const ccomplex_t IT_0929 = s_23*s_34;
    const ccomplex_t IT_0930 = (-2)*IT_0929;
    const ccomplex_t IT_0931 = IT_0873 + IT_0930;
    const ccomplex_t IT_0932 = IT_0206*IT_0931;
    const ccomplex_t IT_0933 = (-0.5)*IT_0872 + (-0.25)*IT_0932;
    const ccomplex_t IT_0934 = s_23 + IT_0933;
    const ccomplex_t IT_0935 = m_N_2*IT_0934;
    const ccomplex_t IT_0936 = (-4)*IT_0935;
    const ccomplex_t IT_0937 = s_23*IT_0618;
    const ccomplex_t IT_0938 = IT_0001*IT_0937;
    const ccomplex_t IT_0939 = 0.666666666666667*IT_0872 + (-0.166666666666667
      )*IT_0938;
    const ccomplex_t IT_0940 = s_23 + IT_0939;
    const ccomplex_t IT_0941 = m_N_2*IT_0940;
    const ccomplex_t IT_0942 = (-6)*IT_0941;
    const ccomplex_t IT_0943 = (-2)*IT_0536;
    const ccomplex_t IT_0944 = (-2)*conj(IT_0536);
    const ccomplex_t IT_0945 = s_13*IT_0618;
    const ccomplex_t IT_0946 = s_14*s_34*IT_0161;
    const ccomplex_t IT_0947 = (-2)*IT_0946;
    const ccomplex_t IT_0948 = IT_0945 + IT_0947;
    const ccomplex_t IT_0949 = IT_0001*IT_0948;
    const ccomplex_t IT_0950 = (-0.666666666666667)*IT_0490 + (
      -0.166666666666667)*IT_0949;
    const ccomplex_t IT_0951 = s_13 + IT_0950;
    const ccomplex_t IT_0952 = m_N_3*IT_0951;
    const ccomplex_t IT_0953 = (-6)*IT_0952;
    const ccomplex_t IT_0954 = (-2)*IT_0543;
    const ccomplex_t IT_0955 = (-2)*conj(IT_0543);
    const ccomplex_t IT_0956 = (-2)*IT_0519;
    const ccomplex_t IT_0957 = IT_0518 + IT_0956;
    const ccomplex_t IT_0958 = IT_0185*IT_0957;
    const ccomplex_t IT_0959 = 0.166666666666667*IT_0949 + (-0.666666666666667
      )*IT_0958;
    const ccomplex_t IT_0960 = s_13 + IT_0959;
    const ccomplex_t IT_0961 = m_N_3*IT_0960;
    const ccomplex_t IT_0962 = 6*IT_0961;
    const ccomplex_t IT_0963 = (-2)*IT_0545;
    const ccomplex_t IT_0964 = (-2)*conj(IT_0545);
    const ccomplex_t IT_0965 = (-4)*IT_0006;
    const ccomplex_t IT_0966 = -IT_0093;
    const ccomplex_t IT_0967 = -conj(IT_0093);
    const ccomplex_t IT_0968 = -IT_0139;
    const ccomplex_t IT_0969 = -conj(IT_0139);
    const ccomplex_t IT_0970 = -IT_0040;
    const ccomplex_t IT_0971 = -conj(IT_0040);
    const ccomplex_t IT_0972 = conj(IT_0093)*IT_0587 + IT_0093*conj(IT_0587) +
       conj(IT_0096)*(IT_0274 + IT_0278 + IT_0281 + -IT_0587 + IT_0596 + IT_0599
      ) + IT_0096*(conj(IT_0274) + conj(IT_0278) + conj(IT_0281) + -conj(IT_0587
      ) + conj(IT_0596) + conj(IT_0599)) + (conj(IT_0274) + conj(IT_0278) + conj
      (IT_0281) + IT_0608)*IT_0966 + (IT_0274 + IT_0278 + IT_0281 + IT_0610)
      *IT_0967 + (conj(IT_0247) + IT_0725)*(IT_0040 + IT_0968) + (IT_0247 +
       IT_0727)*(conj(IT_0040) + IT_0969) + (conj(IT_0220) + conj(IT_0591) +
       conj(IT_0593))*(IT_0139 + IT_0970) + (IT_0220 + IT_0591 + IT_0593)*(conj
      (IT_0139) + IT_0971);
    const ccomplex_t IT_0973 = IT_0207 + IT_0209;
    const ccomplex_t IT_0974 = -IT_0208;
    const ccomplex_t IT_0975 = IT_0973 + IT_0974;
    const ccomplex_t IT_0976 = IT_0206*IT_0975;
    const ccomplex_t IT_0977 = 0.25*IT_0976;
    const ccomplex_t IT_0978 = s_12 + IT_0216 + IT_0977;
    const ccomplex_t IT_0979 = -IT_0040 + IT_0139;
    const ccomplex_t IT_0980 = conj(IT_0139) + IT_0971;
    const ccomplex_t IT_0981 = 0.125*IT_0725;
    const ccomplex_t IT_0982 = 0.125*IT_0608;
    const ccomplex_t IT_0983 = 0.125*IT_0601;
    const ccomplex_t IT_0984 = 8*IT_0978*(0.125*IT_0226*IT_0604 + 0.125
      *IT_0289*IT_0610 + 0.125*IT_0253*IT_0727 + conj(IT_0587)*IT_0979 + IT_0587
      *IT_0980 + IT_0250*IT_0981 + IT_0285*IT_0982 + IT_0223*IT_0983);
    const ccomplex_t IT_0985 = (-2)*IT_0620 + (-2)*IT_0621 + 4*IT_0622;
    const ccomplex_t IT_0986 = IT_0619 + IT_0985;
    const ccomplex_t IT_0987 = IT_0001*IT_0986;
    const ccomplex_t IT_0988 = 0.25*IT_0987;
    const ccomplex_t IT_0989 = (-8)*conj(IT_0588);
    const ccomplex_t IT_0990 = (-8)*conj(IT_0589);
    const ccomplex_t IT_0991 = 8*IT_0593;
    const ccomplex_t IT_0992 = (-8)*IT_0588;
    const ccomplex_t IT_0993 = (-8)*IT_0589;
    const ccomplex_t IT_0994 = IT_0587*conj(IT_0587) + 0.125*IT_0589*IT_0643 +
       0.125*IT_0599*IT_0650 + IT_0596*(conj(IT_0596) + 0.125*IT_0646 + 0.125
      *IT_0650) + 0.125*conj(IT_0589)*IT_0662 + IT_0642*IT_0981 + (IT_0645 +
       IT_0649)*IT_0982 + 0.125*IT_0593*(IT_0989 + IT_0990) + IT_0591*(IT_0601 +
       0.125*IT_0989 + 0.125*IT_0990) + IT_0983*(IT_0991 + IT_0992 + IT_0993);
    const ccomplex_t IT_0995 = IT_0000*IT_0863;
    const ccomplex_t IT_0996 = s_12 + IT_0630;
    const ccomplex_t IT_0997 = IT_0000*IT_0996;
    const ccomplex_t IT_0998 = -IT_0997;
    const ccomplex_t IT_0999 = IT_0995 + IT_0998;
    const ccomplex_t IT_1000 = IT_0000*IT_0851;
    const ccomplex_t IT_1001 = -IT_1000;
    const ccomplex_t IT_1002 = IT_0997 + IT_1001;
    const ccomplex_t IT_1003 = -IT_0487;
    const ccomplex_t IT_1004 = 0.166666666666667*IT_0190;
    const ccomplex_t IT_1005 = 0.166666666666667*IT_0335;
    const ccomplex_t IT_1006 = 0.166666666666667*IT_0202;
    const ccomplex_t IT_1007 = 0.166666666666667*IT_0341;
    const ccomplex_t IT_1008 = 0.166666666666667*IT_0747;
    const ccomplex_t IT_1009 = 4*IT_0386;
    const ccomplex_t IT_1010 = 0.166666666666667*IT_0350;
    const ccomplex_t IT_1011 = 0.166666666666667*IT_0770;
    const ccomplex_t IT_1012 = 0.166666666666667*IT_0797;
    const ccomplex_t IT_1013 = 0.166666666666667*IT_0799;
    const ccomplex_t IT_1014 = 0.166666666666667*IT_0755;
    const ccomplex_t IT_1015 = 0.166666666666667*IT_0821;
    const ccomplex_t IT_1016 = IT_0684*IT_0999 + IT_1002*IT_1003 + IT_0731
      *IT_1004 + IT_0675*IT_1005 + IT_0680*IT_1006 + IT_0734*IT_1007 + IT_0782
      *IT_1008 + (IT_0422 + IT_0447 + IT_1009)*IT_1010 + IT_0744*IT_1011 +
       IT_0857*IT_1012 + IT_0822*IT_1013 + IT_0683*IT_1014 + IT_0785*IT_1015;
    const ccomplex_t IT_1017 = -conj(IT_0487);
    const ccomplex_t IT_1018 = conj(IT_0684)*IT_0999 + IT_0193*IT_1004 +
       IT_0338*IT_1005 + IT_0205*IT_1006 + IT_0344*IT_1007 + IT_0771*IT_1008 + 
      (IT_0351 + IT_0449)*IT_1010 + IT_0745*IT_1011 + IT_0858*IT_1012 + IT_0823
      *IT_1013 + IT_0694*IT_1014 + IT_0786*IT_1015 + IT_1002*IT_1017;
    const ccomplex_t IT_1019 = -IT_0326;
    const ccomplex_t IT_1020 = 0.166666666666667*IT_0866;
    const ccomplex_t IT_1021 = 0.166666666666667*IT_0419;
    const ccomplex_t IT_1022 = 0.166666666666667*IT_0566;
    const ccomplex_t IT_1023 = -conj(IT_0326);
    const ccomplex_t IT_1024 = 0.5*IT_0196 + -IT_0789 + (-0.5)*IT_0790;
    const ccomplex_t IT_1025 = IT_0194 + IT_1024;
    const ccomplex_t IT_1026 = IT_0185*IT_1025;
    const ccomplex_t IT_1027 = -IT_0197;
    const ccomplex_t IT_1028 = IT_0007 + IT_1026 + IT_1027;
    const ccomplex_t IT_1029 = m_N_3*IT_1028;
    const ccomplex_t IT_1030 = (s_12 + -IT_0187)*(s_12 + (-0.5)*IT_0494);
    const ccomplex_t IT_1031 = (-2)*IT_1030;
    const ccomplex_t IT_1032 = 0.25*IT_0335;
    const ccomplex_t IT_1033 = 0.25*IT_0202;
    const ccomplex_t IT_1034 = m_N_3*IT_0189;
    const ccomplex_t IT_1035 = 0.25*IT_1034;
    const ccomplex_t IT_1036 = m_N_2*IT_0340;
    const ccomplex_t IT_1037 = 0.25*IT_1036;
    const ccomplex_t IT_1038 = 0.25*IT_0391;
    const ccomplex_t IT_1039 = 0.25*IT_0810;
    const ccomplex_t IT_1040 = 0.25*IT_0811;
    const ccomplex_t IT_1041 = 0.25*IT_0392;
    const ccomplex_t IT_1042 = 0.25*IT_0847;
    const ccomplex_t IT_1043 = 0.25*IT_0899;
    const ccomplex_t IT_1044 = IT_0587*IT_0722 + -(IT_0220 + -IT_0247)*IT_1029
       + IT_0394*IT_1031 + IT_0734*IT_1032 + IT_0731*IT_1033 + IT_0675*IT_1035 +
       IT_0680*IT_1037 + IT_0422*IT_1038 + IT_0748*IT_1039 + IT_0800*IT_1040 +
       IT_0447*IT_1041 + IT_0857*IT_1042 + IT_0683*IT_1043;
    const ccomplex_t IT_1045 = conj(IT_0587)*IT_0722 + -(conj(IT_0220) + -conj
      (IT_0247))*IT_1029 + IT_0344*IT_1032 + IT_0193*IT_1033 + IT_0338*IT_1035 +
       IT_0205*IT_1037 + IT_0351*IT_1038 + IT_0771*IT_1039 + IT_0823*IT_1040 +
       IT_0449*IT_1041 + IT_0858*IT_1042 + IT_0694*IT_1043;
    const ccomplex_t IT_1046 = s_12 + IT_0629;
    const ccomplex_t IT_1047 = IT_0000*IT_1046;
    const ccomplex_t IT_1048 = 0.25*IT_0839;
    const ccomplex_t IT_1049 = 0.25*IT_0807;
    const ccomplex_t IT_1050 = 0.25*IT_0892;
    const ccomplex_t IT_1051 = 0.25*IT_0871;
    const ccomplex_t IT_1052 = IT_0604*IT_0722 + IT_0286*IT_1029 + IT_0421
      *IT_1031 + IT_0675*IT_1032 + IT_0680*IT_1033 + IT_0734*IT_1035 + IT_0731
      *IT_1037 + IT_0685*IT_1038 + IT_0787*IT_1039 + IT_0825*IT_1040 + IT_0686
      *IT_1041 + (-3)*IT_0394*IT_1047 + IT_0785*IT_1048 + IT_0857*IT_1049 +
       IT_0683*IT_1050 + IT_0744*IT_1051;
    const ccomplex_t IT_1053 = (-0.333333333333333)*IT_0601*IT_0722 + (
      -0.333333333333333)*IT_0282*IT_1029 + (-0.333333333333333)*IT_0338*IT_1032
       + (-0.333333333333333)*IT_0205*IT_1033 + (-0.333333333333333)*IT_0344
      *IT_1035 + (-0.333333333333333)*IT_0193*IT_1037 + (-0.333333333333333)
      *IT_0695*IT_1038 + (-0.333333333333333)*IT_0788*IT_1039 + (
      -0.333333333333333)*IT_0826*IT_1040 + (-0.333333333333333)*IT_0696*IT_1041
       + conj(IT_0394)*IT_1047 + (-0.333333333333333)*IT_0786*IT_1048 + (
      -0.333333333333333)*IT_0858*IT_1049 + (-0.333333333333333)*IT_0694*IT_1050
       + (-0.333333333333333)*IT_0745*IT_1051;
    const ccomplex_t IT_1054 = (-2)*IT_0620;
    const ccomplex_t IT_1055 = IT_0619 + IT_1054;
    const ccomplex_t IT_1056 = IT_0185*IT_1055;
    const ccomplex_t IT_1057 = 0.5*IT_1056;
    const ccomplex_t IT_1058 = (-2)*IT_0620 + 2*IT_0622 + 2*IT_0623 + (-2)
      *IT_0624;
    const ccomplex_t IT_1059 = IT_0619 + IT_1058;
    const ccomplex_t IT_1060 = IT_0185*IT_1059;
    const ccomplex_t IT_1061 = (-0.5)*IT_1060;
    const ccomplex_t IT_1062 = s_14*s_23*IT_0161;
    const ccomplex_t IT_1063 = s_13*s_24*IT_0161;
    const ccomplex_t IT_1064 = IT_1062 + IT_1063;
    const ccomplex_t IT_1065 = s_12*s_34*IT_0161;
    const ccomplex_t IT_1066 = -IT_1065;
    const ccomplex_t IT_1067 = IT_1064 + IT_1066;
    const ccomplex_t IT_1068 = IT_0206*IT_1067;
    const ccomplex_t IT_1069 = 0.5*IT_1068;
    const ccomplex_t IT_1070 = IT_0214 + IT_0492 + IT_1057 + IT_1061 + IT_1069;
    const ccomplex_t IT_1071 = -IT_0354;
    const ccomplex_t IT_1072 = IT_0360 + IT_1071;
    const ccomplex_t IT_1073 = s_12*s_13*IT_0161;
    const ccomplex_t IT_1074 = s_23*IT_0007*IT_0161;
    const ccomplex_t IT_1075 = s_23*IT_0196;
    const ccomplex_t IT_1076 = s_12*s_14*s_34;
    const ccomplex_t IT_1077 = -IT_1074 + 2*IT_1075 + (-2)*IT_1076;
    const ccomplex_t IT_1078 = IT_1073 + IT_1077;
    const ccomplex_t IT_1079 = IT_0185*IT_1078;
    const ccomplex_t IT_1080 = (-0.5)*IT_1079;
    const ccomplex_t IT_1081 = IT_1072 + IT_1080;
    const ccomplex_t IT_1082 = -IT_0358;
    const ccomplex_t IT_1083 = IT_0364 + IT_1082;
    const ccomplex_t IT_1084 = s_12*s_23*IT_0161;
    const ccomplex_t IT_1085 = s_13*IT_0008*IT_0161;
    const ccomplex_t IT_1086 = s_14*s_23*s_24;
    const ccomplex_t IT_1087 = s_14*s_34*IT_0008;
    const ccomplex_t IT_1088 = -IT_1085 + (-2)*IT_1086 + 2*IT_1087;
    const ccomplex_t IT_1089 = IT_1084 + IT_1088;
    const ccomplex_t IT_1090 = IT_0185*IT_1089;
    const ccomplex_t IT_1091 = (-0.5)*IT_1090;
    const ccomplex_t IT_1092 = IT_1083 + IT_1091;
    const ccomplex_t IT_1093 = IT_0000*IT_0161;
    const ccomplex_t IT_1094 = 9*IT_1093;
    const ccomplex_t IT_1095 = 0.25*IT_1094;
    const ccomplex_t IT_1096 = (-2)*IT_0528;
    const ccomplex_t IT_1097 = IT_0527 + IT_1096;
    const ccomplex_t IT_1098 = IT_0185*IT_1097;
    const ccomplex_t IT_1099 = 0.5*IT_0872 + (-0.5)*IT_0874 + (-0.5)*IT_1098;
    const ccomplex_t IT_1100 = s_23 + IT_1099;
    const ccomplex_t IT_1101 = m_N_2*IT_1100;
    const ccomplex_t IT_1102 = 2*IT_1101;
    const ccomplex_t IT_1103 = 0.25*IT_1102;
    const ccomplex_t IT_1104 = IT_0734*IT_1103;
    const ccomplex_t IT_1105 = 0.5*IT_0880 + (-0.5)*IT_0920 + (-0.5)*IT_0958;
    const ccomplex_t IT_1106 = s_13 + IT_1105;
    const ccomplex_t IT_1107 = m_N_3*IT_1106;
    const ccomplex_t IT_1108 = 2*IT_1107;
    const ccomplex_t IT_1109 = 0.25*IT_1108;
    const ccomplex_t IT_1110 = 2*IT_0588;
    const ccomplex_t IT_1111 = 2*IT_0589;
    const ccomplex_t IT_1112 = IT_1110 + IT_1111;
    const ccomplex_t IT_1113 = 0.25*IT_0962;
    const ccomplex_t IT_1114 = (-2)*IT_0274;
    const ccomplex_t IT_1115 = (-2)*IT_0278;
    const ccomplex_t IT_1116 = (-2)*IT_0281;
    const ccomplex_t IT_1117 = IT_1114 + IT_1115 + IT_1116;
    const ccomplex_t IT_1118 = s_24*s_34*IT_0161;
    const ccomplex_t IT_1119 = (-2)*IT_1118;
    const ccomplex_t IT_1120 = IT_0937 + IT_1119;
    const ccomplex_t IT_1121 = IT_0001*IT_1120;
    const ccomplex_t IT_1122 = (-0.666666666666667)*IT_0497 + (
      -0.166666666666667)*IT_1121;
    const ccomplex_t IT_1123 = s_23 + IT_1122;
    const ccomplex_t IT_1124 = m_N_2*IT_1123;
    const ccomplex_t IT_1125 = (-6)*IT_1124;
    const ccomplex_t IT_1126 = 0.25*IT_1125;
    const ccomplex_t IT_1127 = (-2)*IT_0247;
    const ccomplex_t IT_1128 = IT_0683 + IT_1127;
    const ccomplex_t IT_1129 = m_N_3*IT_0185;
    const ccomplex_t IT_1130 = IT_0185*IT_0945;
    const ccomplex_t IT_1131 = (-2)*IT_0519 + (-0.5)*IT_1130;
    const ccomplex_t IT_1132 = IT_0518 + IT_1131;
    const ccomplex_t IT_1133 = IT_1129*IT_1132;
    const ccomplex_t IT_1134 = (-2)*IT_1133;
    const ccomplex_t IT_1135 = 0.25*IT_1134;
    const ccomplex_t IT_1136 = 2*IT_0596;
    const ccomplex_t IT_1137 = 2*IT_0599;
    const ccomplex_t IT_1138 = IT_0857 + IT_1136 + IT_1137;
    const ccomplex_t IT_1139 = (-10)*s_23;
    const ccomplex_t IT_1140 = IT_0938 + IT_1139;
    const ccomplex_t IT_1141 = m_N_2*IT_1140;
    const ccomplex_t IT_1142 = 0.25*IT_1141;
    const ccomplex_t IT_1143 = -IT_0002;
    const ccomplex_t IT_1144 = IT_0618 + IT_1143;
    const ccomplex_t IT_1145 = IT_0185*IT_1144;
    const ccomplex_t IT_1146 = -IT_0161 + 0.25*IT_1145;
    const ccomplex_t IT_1147 = IT_0687 + IT_1146;
    const ccomplex_t IT_1148 = IT_0000*IT_1147;
    const ccomplex_t IT_1149 = 4*IT_1148;
    const ccomplex_t IT_1150 = 0.25*IT_1149;
    const ccomplex_t IT_1151 = IT_0570*IT_1150;
    const ccomplex_t IT_1152 = s_34*IT_0161;
    const ccomplex_t IT_1153 = IT_0206*IT_1152;
    const ccomplex_t IT_1154 = IT_0185*IT_0657;
    const ccomplex_t IT_1155 = (-0.333333333333333)*IT_1153 + (
      -0.333333333333333)*IT_1154;
    const ccomplex_t IT_1156 = IT_0161 + IT_1155;
    const ccomplex_t IT_1157 = IT_0000*IT_1156;
    const ccomplex_t IT_1158 = 3*IT_1157;
    const ccomplex_t IT_1159 = 0.25*IT_1158;
    const ccomplex_t IT_1160 = IT_0906*IT_1159;
    const ccomplex_t IT_1161 = (-4)*IT_0501;
    const ccomplex_t IT_1162 = s_12*IT_0002*IT_0185;
    const ccomplex_t IT_1163 = -IT_1162;
    const ccomplex_t IT_1164 = IT_0619 + IT_0622 + IT_0623;
    const ccomplex_t IT_1165 = (-2)*IT_0620 + -IT_0624;
    const ccomplex_t IT_1166 = IT_1164 + IT_1165;
    const ccomplex_t IT_1167 = IT_0185*IT_1166;
    const ccomplex_t IT_1168 = (-0.5)*IT_1167;
    const ccomplex_t IT_1169 = -IT_0814;
    const ccomplex_t IT_1170 = IT_0214 + IT_0763 + IT_1163 + IT_1168 + IT_1169;
    const ccomplex_t IT_1171 = 0.25*IT_1170;
    const ccomplex_t IT_1172 = s_13*s_14*s_24;
    const ccomplex_t IT_1173 = IT_1075 + IT_1172;
    const ccomplex_t IT_1174 = -IT_1074 + -IT_1076;
    const ccomplex_t IT_1175 = IT_1173 + IT_1174;
    const ccomplex_t IT_1176 = IT_0185*IT_1175;
    const ccomplex_t IT_1177 = (-0.5)*IT_1176;
    const ccomplex_t IT_1178 = -IT_0355;
    const ccomplex_t IT_1179 = IT_0354 + IT_0370 + IT_0373 + IT_1177 + IT_1178;
    const ccomplex_t IT_1180 = 0.25*IT_1179;
    const ccomplex_t IT_1181 = (-8)*IT_0394;
    const ccomplex_t IT_1182 = s_13*IT_0329;
    const ccomplex_t IT_1183 = IT_1086 + IT_1182;
    const ccomplex_t IT_1184 = s_12*s_24*s_34;
    const ccomplex_t IT_1185 = -IT_1085 + -IT_1184;
    const ccomplex_t IT_1186 = IT_1183 + IT_1185;
    const ccomplex_t IT_1187 = IT_0185*IT_1186;
    const ccomplex_t IT_1188 = (-0.25)*IT_1187;
    const ccomplex_t IT_1189 = (-0.5)*IT_0359;
    const ccomplex_t IT_1190 = 0.5*IT_0376;
    const ccomplex_t IT_1191 = (-0.5)*IT_0378;
    const ccomplex_t IT_1192 = (-0.5)*IT_0358;
    const ccomplex_t IT_1193 = IT_0364 + IT_1188 + IT_1189 + IT_1190 + IT_1191
       + IT_1192;
    const ccomplex_t IT_1194 = 0.25*IT_1193;
    const ccomplex_t IT_1195 = IT_0513*IT_1070 + -IT_0684*IT_1081 + IT_0487
      *IT_1092 + IT_0943*IT_1095 + IT_1104 + IT_0731*IT_1109 + IT_1112*IT_1113 +
       IT_1117*IT_1126 + IT_1128*IT_1135 + IT_1138*IT_1142 + IT_1151 + IT_1160 +
       IT_1161*IT_1171 + IT_0723*IT_1180 + IT_1181*IT_1194;
    const ccomplex_t IT_1196 = (-4)*IT_0386;
    const ccomplex_t IT_1197 = IT_0680*IT_1109;
    const ccomplex_t IT_1198 = (-4)*IT_0505;
    const ccomplex_t IT_1199 = (-8)*IT_0421;
    const ccomplex_t IT_1200 = IT_0382*IT_1070 + (IT_0326 + 0.25*IT_1009)
      *IT_1092 + IT_0675*IT_1103 + IT_1126*IT_1128 + IT_1117*IT_1135 + IT_1113
      *IT_1138 + IT_1112*IT_1142 + IT_0560*IT_1150 + IT_0905*IT_1159 + IT_0726
      *IT_1180 + -IT_1081*(IT_0183 + (-0.25)*IT_1196) + IT_1197 + IT_1171
      *IT_1198 + IT_1194*IT_1199;
    const ccomplex_t IT_1201 = s_24*s_34*IT_0007;
    const ccomplex_t IT_1202 = -IT_1074 + (-2)*IT_1172 + 2*IT_1201;
    const ccomplex_t IT_1203 = IT_1073 + IT_1202;
    const ccomplex_t IT_1204 = IT_0185*IT_1203;
    const ccomplex_t IT_1205 = (-0.5)*IT_1204;
    const ccomplex_t IT_1206 = IT_1072 + IT_1205;
    const ccomplex_t IT_1207 = m_N_2*IT_0185;
    const ccomplex_t IT_1208 = IT_0185*IT_0937;
    const ccomplex_t IT_1209 = (-2)*IT_0528 + (-0.5)*IT_1208;
    const ccomplex_t IT_1210 = IT_0527 + IT_1209;
    const ccomplex_t IT_1211 = IT_1207*IT_1210;
    const ccomplex_t IT_1212 = (-2)*IT_1211;
    const ccomplex_t IT_1213 = IT_0001*IT_0945;
    const ccomplex_t IT_1214 = 0.666666666666667*IT_0880 + (-0.166666666666667
      )*IT_1213;
    const ccomplex_t IT_1215 = s_13 + IT_1214;
    const ccomplex_t IT_1216 = m_N_3*IT_1215;
    const ccomplex_t IT_1217 = (-6)*IT_1216;
    const ccomplex_t IT_1218 = -IT_1085 + 2*IT_1182 + (-2)*IT_1184;
    const ccomplex_t IT_1219 = IT_1084 + IT_1218;
    const ccomplex_t IT_1220 = IT_0185*IT_1219;
    const ccomplex_t IT_1221 = (-0.5)*IT_1220;
    const ccomplex_t IT_1222 = IT_1083 + IT_1221;
    const ccomplex_t IT_1223 = 0.25*IT_0953;
    const ccomplex_t IT_1224 = (-0.666666666666667)*IT_1098 +
       0.166666666666667*IT_1121;
    const ccomplex_t IT_1225 = s_23 + IT_1224;
    const ccomplex_t IT_1226 = m_N_2*IT_1225;
    const ccomplex_t IT_1227 = 6*IT_1226;
    const ccomplex_t IT_1228 = 0.25*IT_1227;
    const ccomplex_t IT_1229 = (-10)*s_13;
    const ccomplex_t IT_1230 = IT_1213 + IT_1229;
    const ccomplex_t IT_1231 = m_N_3*IT_1230;
    const ccomplex_t IT_1232 = 0.25*IT_1231;
    const ccomplex_t IT_1233 = m_N_2*m_N_3*IT_0185;
    const ccomplex_t IT_1234 = (-0.25)*IT_0618;
    const ccomplex_t IT_1235 = IT_0002 + IT_1234;
    const ccomplex_t IT_1236 = IT_1233*IT_1235;
    const ccomplex_t IT_1237 = (-4)*IT_1236;
    const ccomplex_t IT_1238 = 0.25*IT_1237;
    const ccomplex_t IT_1239 = -IT_0763;
    const ccomplex_t IT_1240 = IT_0214 + IT_0814 + IT_1163 + IT_1168 + IT_1239;
    const ccomplex_t IT_1241 = 0.25*IT_1240;
    const ccomplex_t IT_1242 = (-0.25)*IT_1176;
    const ccomplex_t IT_1243 = (-0.5)*IT_0355;
    const ccomplex_t IT_1244 = (-0.5)*IT_0370;
    const ccomplex_t IT_1245 = 0.5*IT_0372;
    const ccomplex_t IT_1246 = (-0.5)*IT_0354;
    const ccomplex_t IT_1247 = IT_0360 + IT_1242 + IT_1243 + IT_1244 + IT_1245
       + IT_1246;
    const ccomplex_t IT_1248 = 0.25*IT_1247;
    const ccomplex_t IT_1249 = (-0.5)*IT_1187;
    const ccomplex_t IT_1250 = -IT_0359;
    const ccomplex_t IT_1251 = IT_0358 + IT_0377 + IT_0378 + IT_1249 + IT_1250;
    const ccomplex_t IT_1252 = 0.25*IT_1251;
    const ccomplex_t IT_1253 = IT_0513*IT_1070 + IT_1104 + IT_0731*IT_1109 +
       IT_1151 + IT_1160 + IT_0684*IT_1206 + 0.25*IT_1138*IT_1212 + 0.25*IT_1127
      *IT_1217 + IT_1003*IT_1222 + IT_1112*IT_1223 + IT_1117*IT_1228 + IT_0683
      *IT_1232 + IT_0943*IT_1238 + IT_1161*IT_1241 + IT_1181*IT_1248 + IT_0723
      *IT_1252;
    const ccomplex_t IT_1254 = pow(m_W, 6);
    const ccomplex_t IT_1255 = s_12*IT_1254;
    const ccomplex_t IT_1256 = s_14*s_24*IT_0618;
    const ccomplex_t IT_1257 = (-2)*IT_1256;
    const ccomplex_t IT_1258 = IT_1255 + IT_1257;
    const ccomplex_t IT_1259 = IT_0001*IT_1258;
    const ccomplex_t IT_1260 = 0.5*IT_1259;
    const ccomplex_t IT_1261 = (-4)*IT_0491;
    const ccomplex_t IT_1262 = IT_0214 + IT_1260 + IT_1261;
    const ccomplex_t IT_1263 = IT_0185*IT_0759;
    const ccomplex_t IT_1264 = -IT_1263;
    const ccomplex_t IT_1265 = IT_1262 + IT_1264;
    const ccomplex_t IT_1266 = (-4)*IT_0182;
    const ccomplex_t IT_1267 = (-0.25)*IT_1266;
    const ccomplex_t IT_1268 = IT_0164 + IT_1267;
    const ccomplex_t IT_1269 = (-2)*IT_0593;
    const ccomplex_t IT_1270 = (-0.5)*IT_1269;
    const ccomplex_t IT_1271 = (-0.5)*IT_1111;
    const ccomplex_t IT_1272 = (-0.5)*IT_1110;
    const ccomplex_t IT_1273 = IT_0591 + IT_1270 + IT_1271 + IT_1272;
    const ccomplex_t IT_1274 = 4*IT_0183;
    const ccomplex_t IT_1275 = IT_1009 + IT_1274;
    const ccomplex_t IT_1276 = (-4)*IT_0326;
    const ccomplex_t IT_1277 = IT_1196 + IT_1276;
    const ccomplex_t IT_1278 = 0.25*IT_1259;
    const ccomplex_t IT_1279 = IT_0214 + IT_1278;
    const ccomplex_t IT_1280 = (-8)*IT_0536;
    const ccomplex_t IT_1281 = IT_0382*IT_1070 + IT_0954*IT_1095 + IT_0675
      *IT_1103 + IT_0560*IT_1150 + IT_0905*IT_1159 + IT_1197 + IT_1138*IT_1223 +
       IT_1128*IT_1228 + IT_1117*IT_1232 + IT_0963*IT_1238 + IT_1198*IT_1241 +
       IT_1199*IT_1248 + IT_0726*IT_1252 + -IT_1265*IT_1268 + (-0.5)*IT_1212
      *IT_1273 + 0.25*IT_1206*IT_1275 + 0.25*IT_1222*IT_1277 + 0.25*IT_1279
      *IT_1280;
    const ccomplex_t IT_1282 = s_23*m_N_2;
    const ccomplex_t IT_1283 = IT_0005*IT_1282;
    const ccomplex_t IT_1284 = s_13*m_N_3;
    const ccomplex_t IT_1285 = IT_0005*IT_1284;
    const ccomplex_t IT_1286 = m_N_2*m_N_3*IT_0161;
    const ccomplex_t IT_1287 = IT_0005*IT_1286;
    const ccomplex_t IT_1288 = (-4)*IT_0513;
    const ccomplex_t IT_1289 = (1 + 0.5*IT_0003)*(IT_0214 + (-2)*IT_0491);
    const ccomplex_t IT_1290 = -IT_0382;
    const ccomplex_t IT_1291 = 0.25*IT_0553;
    const ccomplex_t IT_1292 = 0.25*IT_0557;
    const ccomplex_t IT_1293 = 0.25*IT_0878;
    const ccomplex_t IT_1294 = s_34*IT_0185;
    const ccomplex_t IT_1295 = IT_0207 + IT_0208;
    const ccomplex_t IT_1296 = s_13*s_23*s_34*IT_0185;
    const ccomplex_t IT_1297 = (-2)*IT_1296;
    const ccomplex_t IT_1298 = IT_1295 + IT_1297;
    const ccomplex_t IT_1299 = IT_1294*IT_1298;
    const ccomplex_t IT_1300 = 0.25*IT_1299;
    const ccomplex_t IT_1301 = IT_0352 + IT_0353 + IT_0354 + IT_0355 + IT_0363
       + IT_0369 + IT_0371 + IT_0373;
    const ccomplex_t IT_1302 = 0.25*IT_1301;
    const ccomplex_t IT_1303 = IT_0356 + IT_0357 + IT_0358 + IT_0359 + IT_0367
       + IT_0375 + IT_0377 + IT_0379;
    const ccomplex_t IT_1304 = 0.25*IT_1303;
    const ccomplex_t IT_1305 = IT_0142*IT_1283 + IT_0979*IT_1285 + 0.25
      *IT_1287*IT_1288 + IT_1289*IT_1290 + IT_0685*IT_1291 + IT_0686*IT_1292 +
       IT_0744*IT_1293 + IT_0907*IT_1300 + IT_0422*IT_1302 + IT_0447*IT_1304;
    const ccomplex_t IT_1306 = 0.25*IT_0889;
    const ccomplex_t IT_1307 = (IT_0040 + 0.25*IT_0714)*IT_1283 + IT_0099
      *IT_1285 + 0.25*IT_1288*IT_1289 + IT_1287*IT_1290 + IT_0422*IT_1291 +
       IT_0447*IT_1292 + IT_0904*IT_1300 + IT_0685*IT_1302 + IT_0686*IT_1304 +
       IT_0744*IT_1306;
    const ccomplex_t IT_1308 = (-0.5)*IT_1102;
    const ccomplex_t IT_1309 = (-0.5)*IT_1108;
    const ccomplex_t IT_1310 = (-0.5)*IT_0953;
    const ccomplex_t IT_1311 = (-0.5)*IT_1227;
    const ccomplex_t IT_1312 = (-0.5)*IT_1231;
    const ccomplex_t IT_1313 = (-0.5)*IT_1149;
    const ccomplex_t IT_1314 = (-0.5)*IT_1158;
    const ccomplex_t IT_1315 = (-0.5)*IT_1237;
    const ccomplex_t IT_1316 = (-0.5)*IT_1240;
    const ccomplex_t IT_1317 = (-0.5)*IT_1279;
    const ccomplex_t IT_1318 = (-0.5)*IT_1206;
    const ccomplex_t IT_1319 = (-0.5)*IT_1247;
    const ccomplex_t IT_1320 = (-0.5)*IT_1222;
    const ccomplex_t IT_1321 = (-0.5)*IT_1251;
    const ccomplex_t IT_1322 = IT_0382*IT_1070 + 0.25*IT_0954*IT_1094 + 
      -IT_1265*IT_1268 + (-0.5)*IT_1212*IT_1273 + (-0.5)*IT_0675*IT_1308 + (-0.5
      )*IT_0680*IT_1309 + (-0.5)*IT_1138*IT_1310 + (-0.5)*IT_1128*IT_1311 + (
      -0.5)*IT_1117*IT_1312 + (-0.5)*IT_0560*IT_1313 + (-0.5)*IT_0905*IT_1314 + 
      (-0.5)*IT_0963*IT_1315 + (-0.5)*IT_1198*IT_1316 + (-0.5)*IT_1280*IT_1317 +
       (-0.5)*IT_1275*IT_1318 + (-0.5)*IT_1199*IT_1319 + (-0.5)*IT_1277*IT_1320 
      + (-0.5)*IT_0726*IT_1321;
    const ccomplex_t IT_1323 = 2*conj(IT_0588);
    const ccomplex_t IT_1324 = 2*conj(IT_0589);
    const ccomplex_t IT_1325 = (-2)*conj(IT_0593);
    const ccomplex_t IT_1326 = 2*conj(IT_0596);
    const ccomplex_t IT_1327 = 2*conj(IT_0599);
    const ccomplex_t IT_1328 = IT_0858 + IT_1326 + IT_1327;
    const ccomplex_t IT_1329 = (-2)*conj(IT_0247);
    const ccomplex_t IT_1330 = IT_0694 + IT_1329;
    const ccomplex_t IT_1331 = (-2)*conj(IT_0274);
    const ccomplex_t IT_1332 = (-2)*conj(IT_0278);
    const ccomplex_t IT_1333 = (-2)*conj(IT_0281);
    const ccomplex_t IT_1334 = IT_1331 + IT_1332 + IT_1333;
    const ccomplex_t IT_1335 = (-4)*conj(IT_0505);
    const ccomplex_t IT_1336 = (-8)*conj(IT_0536);
    const ccomplex_t IT_1337 = 4*conj(IT_0183);
    const ccomplex_t IT_1338 = 4*conj(IT_0386);
    const ccomplex_t IT_1339 = (-8)*conj(IT_0421);
    const ccomplex_t IT_1340 = (-4)*conj(IT_0386);
    const ccomplex_t IT_1341 = (-4)*conj(IT_0326);
    const ccomplex_t IT_1342 = conj(IT_0382)*IT_1070 + 0.25*IT_0955*IT_1094 + 
      (-0.5)*IT_0338*IT_1308 + (-0.5)*IT_0205*IT_1309 + (-0.5)*IT_0562*IT_1313 +
       (-0.5)*IT_0879*IT_1314 + (-0.5)*IT_0964*IT_1315 + (-0.5)*IT_0728*IT_1321 
      + (-0.5)*IT_1212*(conj(IT_0591) + (-0.5)*IT_1323 + (-0.5)*IT_1324 + (-0.5)
      *IT_1325) + (-0.5)*IT_1310*IT_1328 + (-0.5)*IT_1311*IT_1330 + (-0.5)
      *IT_1312*IT_1334 + (-0.5)*IT_1316*IT_1335 + (-0.5)*IT_1317*IT_1336 + (-0.5
      )*IT_1318*(IT_1337 + IT_1338) + (-0.5)*IT_1319*IT_1339 + (-0.5)*IT_1320*
      (IT_1340 + IT_1341);
    const ccomplex_t IT_1343 = IT_1323 + IT_1324;
    const ccomplex_t IT_1344 = (-8)*conj(IT_0394);
    const ccomplex_t IT_1345 = (-4)*conj(IT_0501);
    const ccomplex_t IT_1346 = conj(IT_0513)*IT_1070;
    const ccomplex_t IT_1347 = -IT_1346;
    const ccomplex_t IT_1348 = IT_0344*IT_1102;
    const ccomplex_t IT_1349 = (-0.25)*IT_1348;
    const ccomplex_t IT_1350 = (-0.25)*IT_1108;
    const ccomplex_t IT_1351 = (-0.25)*IT_1149;
    const ccomplex_t IT_1352 = IT_0572*IT_1351;
    const ccomplex_t IT_1353 = (-0.25)*IT_1158;
    const ccomplex_t IT_1354 = IT_0890*IT_1353;
    const ccomplex_t IT_1355 = -conj(IT_0684)*IT_1206 + conj(IT_0487)*IT_1222 
      + (-0.25)*IT_0694*IT_1231 + (-0.25)*IT_0944*IT_1237 + (-0.25)*IT_0724
      *IT_1251 + conj(IT_0543)*IT_1265 + 2*conj(IT_0545)*IT_1279 + (-0.25)
      *IT_1212*IT_1328 + (-0.25)*IT_1217*IT_1329 + (-0.25)*IT_1227*IT_1334 + (
      -0.25)*IT_0953*IT_1343 + (-0.25)*IT_1247*IT_1344 + (-0.25)*IT_1240*IT_1345
       + IT_1347 + IT_1349 + IT_0193*IT_1350 + IT_1352 + IT_1354;
    const ccomplex_t IT_1356 = IT_0185*IT_0986;
    const ccomplex_t IT_1357 = -IT_1356;
    const ccomplex_t IT_1358 = IT_1262 + IT_1357;
    const ccomplex_t IT_1359 = (-0.25)*IT_0962;
    const ccomplex_t IT_1360 = (-0.25)*IT_1125;
    const ccomplex_t IT_1361 = (-0.25)*IT_1134;
    const ccomplex_t IT_1362 = (-0.25)*IT_1141;
    const ccomplex_t IT_1363 = (-0.25)*IT_1170;
    const ccomplex_t IT_1364 = (-0.25)*IT_1179;
    const ccomplex_t IT_1365 = (-0.25)*IT_1193;
    const ccomplex_t IT_1366 = conj(IT_0684)*IT_1081 + IT_1017*IT_1092 + (
      -0.25)*IT_0944*IT_1094 + 2*conj(IT_0543)*IT_1279 + IT_1347 + IT_1349 +
       IT_0193*IT_1350 + IT_1352 + IT_1354 + conj(IT_0545)*IT_1358 + IT_1343
      *IT_1359 + IT_1334*IT_1360 + IT_1330*IT_1361 + IT_1328*IT_1362 + IT_1345
      *IT_1363 + IT_0724*IT_1364 + IT_1344*IT_1365;
    const ccomplex_t IT_1367 = -conj(IT_0382)*IT_1070 + (-0.25)*IT_0338
      *IT_1102 + -IT_1092*(conj(IT_0326) + 0.25*IT_1338) + IT_1081*(conj(IT_0183
      ) + (-0.25)*IT_1340) + IT_0205*IT_1350 + IT_0562*IT_1351 + IT_0879*IT_1353
       + conj(IT_0536)*IT_1358 + IT_1328*IT_1359 + IT_1330*IT_1360 + IT_1334
      *IT_1361 + IT_1343*IT_1362 + IT_1335*IT_1363 + IT_0728*IT_1364 + IT_1339
      *IT_1365;
    const ccomplex_t IT_1368 = (-0.166666666666667)*IT_0341;
    const ccomplex_t IT_1369 = (-0.166666666666667)*IT_0350;
    const ccomplex_t IT_1370 = (-0.166666666666667)*IT_0860;
    const ccomplex_t IT_1371 = (-0.166666666666667)*IT_0569;
    const ccomplex_t IT_1372 = IT_0782*IT_0866;
    const ccomplex_t IT_1373 = IT_0755*IT_0825;
    const ccomplex_t IT_1374 = IT_0855*IT_0860;
    const ccomplex_t IT_1375 = IT_0832*IT_0899;
    const ccomplex_t IT_1376 = IT_0808*IT_0892;
    const ccomplex_t IT_1377 = IT_0395*IT_0871;
    const ccomplex_t IT_1378 = IT_0682*IT_0904;
    const ccomplex_t IT_1379 = IT_0889*IT_0906;
    const ccomplex_t IT_1380 = IT_0878*IT_0905;
    const ccomplex_t IT_1381 = IT_0695*IT_0770;
    const ccomplex_t IT_1382 = IT_0784*IT_0866;
    const ccomplex_t IT_1383 = IT_0755*IT_0826;
    const ccomplex_t IT_1384 = IT_0856*IT_0860;
    const ccomplex_t IT_1385 = IT_0840*IT_0899;
    const ccomplex_t IT_1386 = IT_0572*IT_0743;
    const ccomplex_t IT_1387 = IT_0809*IT_0892;
    const ccomplex_t IT_1388 = IT_0397*IT_0871;
    const ccomplex_t IT_1389 = IT_0682*IT_0900;
    const ccomplex_t IT_1390 = IT_0889*IT_0890;
    const ccomplex_t IT_1391 = IT_0878*IT_0879;
    const ccomplex_t IT_1392 = 4*IT_0006*(conj(IT_0040)*IT_0099 + IT_0040
      *IT_0102 + conj(IT_0139)*IT_0142 + IT_0139*IT_0145) + IT_0183*(IT_0190
      *IT_0193 + IT_0202*IT_0205) + IT_0217*(conj(IT_0220)*IT_0223 + IT_0220
      *IT_0226 + conj(IT_0247)*IT_0250 + IT_0247*IT_0253 + IT_0282*IT_0285 +
       IT_0286*IT_0289) + IT_0326*(IT_0335*IT_0338 + IT_0341*IT_0344 + IT_0350
      *IT_0351) + (IT_0352 + IT_0353 + IT_0354 + IT_0355 + IT_0356 + IT_0357 +
       IT_0358 + IT_0359 + IT_0363 + IT_0367 + IT_0369 + IT_0371 + IT_0373 +
       IT_0375 + IT_0377 + IT_0379)*(IT_0388 + IT_0390) + (IT_0391 + IT_0392)*
      (IT_0396 + IT_0398) + IT_0419*(IT_0423 + IT_0424 + IT_0426 + IT_0428) + 
      (IT_0436 + IT_0446)*(IT_0426 + IT_0428 + IT_0448 + IT_0450 + IT_0488 +
       IT_0489) + (IT_0496 + IT_0499)*(IT_0502 + IT_0503 + IT_0506 + IT_0507) + 
      (IT_0509 + IT_0511)*(IT_0514 + IT_0515 + IT_0516 + IT_0517) + (IT_0526 +
       IT_0535)*(IT_0537 + IT_0538 + IT_0539 + IT_0540 + IT_0541 + IT_0542 +
       IT_0544 + IT_0546 + IT_0547 + IT_0548) + (IT_0553 + IT_0557)*(IT_0558 +
       IT_0559 + IT_0561 + IT_0563) + (IT_0566 + IT_0569)*(IT_0571 + IT_0573) + 
      (IT_0576 + IT_0579)*(IT_0581 + IT_0583) + IT_0584*IT_0617 + 8*(s_12 +
       IT_0628 + IT_0629 + IT_0630)*IT_0655 + IT_0661*IT_0670 + conj(IT_0505)*
      (IT_0422*IT_0566 + IT_0447*IT_0569 + IT_0672*IT_0675 + IT_0677*IT_0680 +
       IT_0682*IT_0683 + IT_0553*IT_0685 + IT_0557*IT_0686 + IT_0570*IT_0689 +
       IT_0692*IT_0693) + IT_0505*(IT_0351*IT_0566 + IT_0449*IT_0569 + IT_0338
      *IT_0672 + IT_0205*IT_0677 + IT_0682*IT_0694 + IT_0553*IT_0695 + IT_0557
      *IT_0696 + IT_0692*IT_0697) + (IT_0700 + IT_0703)*(IT_0705 + IT_0707 +
       IT_0708 + IT_0709) + 4*IT_0710*(IT_0093*conj(IT_0093) + IT_0139*conj
      (IT_0139) + 0.25*conj(IT_0096)*IT_0140 + IT_0096*IT_0145 + IT_0040*IT_0713
       + 0.25*conj(IT_0040)*IT_0714) + IT_0722*(IT_0608*IT_0723 + IT_0610
      *IT_0724 + IT_0725*IT_0726 + IT_0727*IT_0728) + conj(IT_0501)*(IT_0422
      *IT_0553 + IT_0447*IT_0557 + IT_0560*IT_0689 + IT_0677*IT_0731 + IT_0672
      *IT_0734 + IT_0743*IT_0744) + IT_0501*(IT_0351*IT_0553 + IT_0449*IT_0557 +
       IT_0344*IT_0672 + IT_0193*IT_0677 + IT_0743*IT_0745) + conj(IT_0183)*
      (IT_0202*IT_0680 + IT_0190*IT_0731 + IT_0747*IT_0748) + conj(IT_0220)*
      (IT_0447*IT_0755 + IT_0651*IT_0762 + IT_0748*IT_0770) + IT_0220*(IT_0449
      *IT_0755 + IT_0770*IT_0771) + IT_0778*(IT_0601*IT_0779 + IT_0604*IT_0780 +
       IT_0725*IT_0782 + IT_0727*IT_0784 + conj(IT_0684)*IT_0785 + IT_0684
      *IT_0786 + IT_0608*IT_0787 + IT_0610*IT_0788) + (IT_0685*IT_0725 + IT_0695
      *IT_0727 + conj(IT_0587)*IT_0748 + IT_0587*IT_0771 + IT_0608*IT_0779 +
       IT_0610*IT_0780 + IT_0601*IT_0787 + IT_0604*IT_0788)*IT_0797 + conj
      (IT_0326)*(IT_0350*IT_0422 + IT_0335*IT_0675 + IT_0341*IT_0734 + IT_0799
      *IT_0800) + IT_0807*(IT_0425*IT_0608 + IT_0427*IT_0610 + IT_0395*IT_0725 +
       IT_0397*IT_0727 + IT_0601*IT_0808 + IT_0604*IT_0809) + (IT_0810 + IT_0811
      )*(IT_0812 + IT_0813) + IT_0821*(IT_0449*IT_0587 + IT_0447*conj(IT_0587) +
       IT_0601*IT_0686 + IT_0604*IT_0696 + IT_0608*IT_0822 + IT_0610*IT_0824 +
       IT_0725*IT_0825 + IT_0727*IT_0826) + (IT_0828 + IT_0830)*(IT_0831 +
       IT_0833) + IT_0839*(IT_0395*IT_0601 + IT_0397*IT_0604 + IT_0725*IT_0808 +
       IT_0727*IT_0809 + IT_0608*IT_0832 + IT_0610*IT_0840) + (IT_0395*IT_0608 +
       IT_0397*IT_0610 + IT_0425*IT_0725 + IT_0427*IT_0727 + IT_0601*IT_0832 +
       IT_0604*IT_0840)*IT_0847 + IT_0854*(IT_0608*IT_0686 + IT_0610*IT_0696 +
       IT_0601*IT_0822 + IT_0604*IT_0824 + IT_0725*IT_0855 + IT_0727*IT_0856 +
       conj(IT_0487)*IT_0857 + IT_0487*IT_0858) + 8*IT_0281*IT_0901 + IT_0903*
      (IT_0570*IT_0608 + IT_0572*IT_0610 + IT_0560*IT_0725 + IT_0562*IT_0727 +
       conj(IT_0505)*IT_0857 + IT_0505*IT_0858 + IT_0604*IT_0900 + IT_0601
      *IT_0904) + conj(IT_0247)*IT_0908 + IT_0247*IT_0910 + (IT_0570*IT_0601 +
       IT_0572*IT_0604 + conj(IT_0501)*IT_0785 + IT_0501*IT_0786 + IT_0610
      *IT_0900 + IT_0608*IT_0904 + IT_0725*IT_0907 + IT_0727*IT_0909)*IT_0919 +
       IT_0924*(conj(IT_0382)*IT_0785 + IT_0382*IT_0786 + IT_0604*IT_0879 +
       IT_0601*IT_0905 + IT_0725*IT_0925 + IT_0727*IT_0926 + IT_0608*IT_0927 +
       IT_0610*IT_0928) + (conj(IT_0513)*IT_0785 + IT_0513*IT_0786 + IT_0604
      *IT_0890 + IT_0601*IT_0906 + IT_0608*IT_0925 + IT_0610*IT_0926 + IT_0725
      *IT_0927 + IT_0727*IT_0928)*IT_0936 + IT_0942*(IT_0601*IT_0943 + IT_0604
      *IT_0944) + IT_0953*(IT_0601*IT_0954 + IT_0604*IT_0955) + IT_0962*(IT_0601
      *IT_0963 + IT_0604*IT_0964) + IT_0965*IT_0972 + IT_0984 + 8*(s_12 +
       IT_0756 + IT_0757 + IT_0988)*IT_0994 + 6*conj(IT_0386)*IT_1016 + 6
      *IT_0386*IT_1018 + 6*conj(IT_0684)*(IT_0183*IT_0995 + IT_0680*IT_1004 +
       IT_0731*IT_1006 + IT_0787*IT_1008 + IT_0997*IT_1019 + IT_0683*IT_1020 +
       IT_0395*IT_1021 + IT_0560*IT_1022) + 6*IT_0684*(conj(IT_0183)*IT_0995 +
       IT_0205*IT_1004 + IT_0193*IT_1006 + IT_0694*IT_1020 + IT_0397*IT_1021 +
       IT_0562*IT_1022 + IT_0997*IT_1023) + 4*conj(IT_0394)*IT_1044 + 4*IT_0394
      *IT_1045 + 4*conj(IT_0421)*IT_1052 + (-12)*IT_0421*IT_1053 + 4*conj
      (IT_0545)*IT_1195 + 4*conj(IT_0536)*IT_1200 + 4*conj(IT_0543)*IT_1253 + 4
      *conj(IT_0182)*IT_1281 + 4*IT_0382*(IT_0145*IT_1283 + IT_0980*IT_1285 +
       IT_0695*IT_1291 + IT_0696*IT_1292 + IT_0745*IT_1293 + IT_0909*IT_1300 +
       IT_0351*IT_1302 + IT_0449*IT_1304) + 4*conj(IT_0382)*IT_1305 + 4*IT_0513*
      (IT_0713*IT_1283 + IT_0102*IT_1285 + IT_0351*IT_1291 + IT_0449*IT_1292 +
       IT_0900*IT_1300 + IT_0695*IT_1302 + IT_0696*IT_1304 + IT_0745*IT_1306) +
       4*conj(IT_0513)*IT_1307 + 4*conj(IT_0164)*IT_1322 + 4*(IT_0164 + IT_0182)
      *IT_1342 + (-4)*IT_0543*IT_1355 + (-4)*IT_0545*IT_1366 + (-4)*IT_0536
      *IT_1367 + (-6)*conj(IT_0487)*((-0.166666666666667)*IT_0335*IT_0734 + (
      -0.166666666666667)*IT_0799*IT_0825 + IT_0183*IT_0997 + IT_1000*IT_1019 +
       IT_0675*IT_1368 + IT_0685*IT_1369 + IT_0744*IT_1370 + IT_0560*IT_1371) + 
      (-6)*IT_0487*((-0.166666666666667)*IT_0335*IT_0344 + conj(IT_0183)*IT_0997
       + IT_1000*IT_1023 + IT_0338*IT_1368 + IT_0695*IT_1369 + IT_0745*IT_1370 +
       IT_0562*IT_1371) + conj(IT_0274)*(IT_0570*IT_0743 + IT_0647*IT_0762 +
       IT_0685*IT_0770 + IT_1372 + IT_1373 + IT_1374 + IT_1375 + IT_1376 +
       IT_1377 + IT_1378 + IT_1379 + IT_1380) + IT_0644*(IT_0570*IT_0743 + 
      (IT_0647 + IT_0648)*IT_0762 + IT_0685*IT_0770 + IT_1372 + IT_1373 +
       IT_1374 + IT_1375 + IT_1376 + IT_1377 + IT_1378 + IT_1379 + IT_1380) +
       IT_0274*(IT_1381 + IT_1382 + IT_1383 + IT_1384 + IT_1385 + IT_1386 +
       IT_1387 + IT_1388 + IT_1389 + IT_1390 + IT_1391) + IT_0278*(IT_0633
      *IT_0762 + IT_1381 + IT_1382 + IT_1383 + IT_1384 + IT_1385 + IT_1386 +
       IT_1387 + IT_1388 + IT_1389 + IT_1390 + IT_1391);
    return create_ccomplex_return(IT_1392);
}

