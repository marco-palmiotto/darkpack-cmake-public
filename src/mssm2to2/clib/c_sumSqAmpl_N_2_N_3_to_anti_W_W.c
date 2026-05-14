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
    const ccomplex_t IT_0016 = cpow(IT_0011, -1);
    const ccomplex_t IT_0017 = IT_0012*IT_0016;
    const ccomplex_t IT_0018 = conj(N_d2)*N_d3*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = IT_0011*IT_0013;
    const ccomplex_t IT_0021 = IT_0018*IT_0020;
    const ccomplex_t IT_0022 = conj(N_u2)*N_u3*e_em;
    const ccomplex_t IT_0023 = IT_0017*IT_0022;
    const ccomplex_t IT_0024 = IT_0020*IT_0022;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0019 + IT_0021 + 
      -IT_0023 + -IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = IT_0015*IT_0026;
    const ccomplex_t IT_0028 = m_N_3*IT_0027;
    const ccomplex_t IT_0029 = N_u2*conj(N_u3)*e_em;
    const ccomplex_t IT_0030 = IT_0017*IT_0029;
    const ccomplex_t IT_0031 = IT_0020*IT_0029;
    const ccomplex_t IT_0032 = N_d2*conj(N_d3)*e_em;
    const ccomplex_t IT_0033 = IT_0017*IT_0032;
    const ccomplex_t IT_0034 = IT_0020*IT_0032;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0030 + IT_0031 + 
      -IT_0033 + -IT_0034);
    const ccomplex_t IT_0036 = (-0.5)*IT_0035;
    const ccomplex_t IT_0037 = IT_0015*IT_0036;
    const ccomplex_t IT_0038 = m_N_2*IT_0037;
    const ccomplex_t IT_0039 = IT_0010*(IT_0028 + -IT_0038);
    const ccomplex_t IT_0040 = -IT_0039;
    const ccomplex_t IT_0041 = -IT_0040;
    const ccomplex_t IT_0042 = cos(alpha);
    const ccomplex_t IT_0043 = IT_0016*IT_0042;
    const ccomplex_t IT_0044 = N_B3*N_u2*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = sin(alpha);
    const ccomplex_t IT_0047 = IT_0016*IT_0046;
    const ccomplex_t IT_0048 = N_B3*N_d2*e_em;
    const ccomplex_t IT_0049 = IT_0047*IT_0048;
    const ccomplex_t IT_0050 = N_B2*N_u3*e_em;
    const ccomplex_t IT_0051 = IT_0043*IT_0050;
    const ccomplex_t IT_0052 = N_B2*N_d3*e_em;
    const ccomplex_t IT_0053 = IT_0047*IT_0052;
    const ccomplex_t IT_0054 = IT_0045 + IT_0049 + IT_0051 + IT_0053;
    const ccomplex_t IT_0055 = IT_0013*IT_0042;
    const ccomplex_t IT_0056 = N_W3*N_u2*e_em;
    const ccomplex_t IT_0057 = IT_0055*IT_0056;
    const ccomplex_t IT_0058 = IT_0013*IT_0046;
    const ccomplex_t IT_0059 = N_W3*N_d2*e_em;
    const ccomplex_t IT_0060 = IT_0058*IT_0059;
    const ccomplex_t IT_0061 = N_W2*N_u3*e_em;
    const ccomplex_t IT_0062 = IT_0055*IT_0061;
    const ccomplex_t IT_0063 = N_W2*N_d3*e_em;
    const ccomplex_t IT_0064 = IT_0058*IT_0063;
    const ccomplex_t IT_0065 = -IT_0057 + -IT_0060 + -IT_0062 + -IT_0064;
    const ccomplex_t IT_0066 = IT_0054 + IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = (-0.5)*IT_0067;
    const ccomplex_t IT_0069 = sin(beta);
    const ccomplex_t IT_0070 = cos(beta);
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*m_W*e_em*IT_0013*(IT_0042
      *IT_0069 + -IT_0046*IT_0070);
    const ccomplex_t IT_0072 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0007 + -IT_0008 + -reg_prop, -1);
    const ccomplex_t IT_0073 = IT_0071*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*IT_0073;
    const ccomplex_t IT_0075 = IT_0068*IT_0074;
    const ccomplex_t IT_0076 = IT_0043*IT_0048;
    const ccomplex_t IT_0077 = IT_0056*IT_0058;
    const ccomplex_t IT_0078 = IT_0043*IT_0052;
    const ccomplex_t IT_0079 = IT_0058*IT_0061;
    const ccomplex_t IT_0080 = IT_0076 + IT_0077 + IT_0078 + IT_0079;
    const ccomplex_t IT_0081 = IT_0044*IT_0047;
    const ccomplex_t IT_0082 = IT_0055*IT_0059;
    const ccomplex_t IT_0083 = IT_0047*IT_0050;
    const ccomplex_t IT_0084 = IT_0055*IT_0063;
    const ccomplex_t IT_0085 = -IT_0081 + -IT_0082 + -IT_0083 + -IT_0084;
    const ccomplex_t IT_0086 = IT_0080 + IT_0085;
    const ccomplex_t IT_0087 = (0 + _Complex_I*1)*IT_0086;
    const ccomplex_t IT_0088 = 0.5*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*m_W*e_em*IT_0013*(IT_0046
      *IT_0069 + IT_0042*IT_0070);
    const ccomplex_t IT_0090 = cpow((-2)*s_12 + -IT_0007 + -IT_0008 + m_H0*
      (m_H0 + (0 + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*IT_0091;
    const ccomplex_t IT_0093 = IT_0088*IT_0092;
    const ccomplex_t IT_0094 = -IT_0075 + -IT_0093;
    const ccomplex_t IT_0095 = m_N_2*IT_0027;
    const ccomplex_t IT_0096 = m_N_3*IT_0037;
    const ccomplex_t IT_0097 = IT_0010*(IT_0095 + -IT_0096);
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = (-4)*IT_0098;
    const ccomplex_t IT_0100 = IT_0094 + 0.25*IT_0099;
    const ccomplex_t IT_0101 = (-4)*conj(IT_0098);
    const ccomplex_t IT_0102 = conj(IT_0094) + 0.25*IT_0101;
    const ccomplex_t IT_0103 = conj(N_B2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0104 = IT_0043*IT_0103;
    const ccomplex_t IT_0105 = conj(N_B2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0106 = IT_0047*IT_0105;
    const ccomplex_t IT_0107 = conj(N_B3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0108 = IT_0043*IT_0107;
    const ccomplex_t IT_0109 = conj(N_B3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0110 = IT_0047*IT_0109;
    const ccomplex_t IT_0111 = IT_0104 + IT_0106 + IT_0108 + IT_0110;
    const ccomplex_t IT_0112 = conj(N_W2)*conj(N_u3)*e_em;
    const ccomplex_t IT_0113 = IT_0055*IT_0112;
    const ccomplex_t IT_0114 = conj(N_W2)*conj(N_d3)*e_em;
    const ccomplex_t IT_0115 = IT_0058*IT_0114;
    const ccomplex_t IT_0116 = conj(N_W3)*conj(N_u2)*e_em;
    const ccomplex_t IT_0117 = IT_0055*IT_0116;
    const ccomplex_t IT_0118 = conj(N_W3)*conj(N_d2)*e_em;
    const ccomplex_t IT_0119 = IT_0058*IT_0118;
    const ccomplex_t IT_0120 = -IT_0113 + -IT_0115 + -IT_0117 + -IT_0119;
    const ccomplex_t IT_0121 = IT_0111 + IT_0120;
    const ccomplex_t IT_0122 = (0 + _Complex_I*1)*IT_0121;
    const ccomplex_t IT_0123 = (-0.5)*IT_0122;
    const ccomplex_t IT_0124 = IT_0074*IT_0123;
    const ccomplex_t IT_0125 = IT_0043*IT_0105;
    const ccomplex_t IT_0126 = IT_0058*IT_0112;
    const ccomplex_t IT_0127 = IT_0043*IT_0109;
    const ccomplex_t IT_0128 = IT_0058*IT_0116;
    const ccomplex_t IT_0129 = IT_0125 + IT_0126 + IT_0127 + IT_0128;
    const ccomplex_t IT_0130 = IT_0047*IT_0103;
    const ccomplex_t IT_0131 = IT_0055*IT_0114;
    const ccomplex_t IT_0132 = IT_0047*IT_0107;
    const ccomplex_t IT_0133 = IT_0055*IT_0118;
    const ccomplex_t IT_0134 = -IT_0130 + -IT_0131 + -IT_0132 + -IT_0133;
    const ccomplex_t IT_0135 = IT_0129 + IT_0134;
    const ccomplex_t IT_0136 = (0 + _Complex_I*1)*IT_0135;
    const ccomplex_t IT_0137 = 0.5*IT_0136;
    const ccomplex_t IT_0138 = IT_0092*IT_0137;
    const ccomplex_t IT_0139 = -IT_0124 + -IT_0138;
    const ccomplex_t IT_0140 = (-4)*IT_0094;
    const ccomplex_t IT_0141 = 0.25*IT_0140;
    const ccomplex_t IT_0142 = IT_0098 + IT_0141;
    const ccomplex_t IT_0143 = (-4)*conj(IT_0094);
    const ccomplex_t IT_0144 = conj(IT_0098) + 0.25*IT_0143;
    const ccomplex_t IT_0145 = pow(m_W, -2);
    const ccomplex_t IT_0146 = s_34*m_N_2*m_N_3*IT_0145;
    const ccomplex_t IT_0147 = s_13*s_34*IT_0145;
    const ccomplex_t IT_0148 = -IT_0147;
    const ccomplex_t IT_0149 = s_14 + IT_0148;
    const ccomplex_t IT_0150 = IT_0146*IT_0149;
    const ccomplex_t IT_0151 = s_23*s_34*IT_0145;
    const ccomplex_t IT_0152 = -IT_0151;
    const ccomplex_t IT_0153 = s_24 + IT_0152;
    const ccomplex_t IT_0154 = IT_0146*IT_0153;
    const ccomplex_t IT_0155 = IT_0150 + IT_0154;
    const ccomplex_t IT_0156 = (-2)*IT_0014;
    const ccomplex_t IT_0157 = IT_0036*IT_0156;
    const ccomplex_t IT_0158 = IT_0010*IT_0157;
    const ccomplex_t IT_0159 = 2*IT_0014;
    const ccomplex_t IT_0160 = IT_0026*IT_0159;
    const ccomplex_t IT_0161 = IT_0010*IT_0160;
    const ccomplex_t IT_0162 = -IT_0161;
    const ccomplex_t IT_0163 = 2*IT_0162;
    const ccomplex_t IT_0164 = conj(IT_0158)*IT_0163;
    const ccomplex_t IT_0165 = 2*conj(IT_0162);
    const ccomplex_t IT_0166 = IT_0158*IT_0165;
    const ccomplex_t IT_0167 = s_13*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0168 = pow(s_13, 2);
    const ccomplex_t IT_0169 = s_24*s_34*IT_0001*IT_0168;
    const ccomplex_t IT_0170 = s_23*IT_0007;
    const ccomplex_t IT_0171 = s_12*s_14*s_34*IT_0145;
    const ccomplex_t IT_0172 = pow(s_23, 2);
    const ccomplex_t IT_0173 = s_14*s_34*IT_0001*IT_0172;
    const ccomplex_t IT_0174 = s_13*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0175 = s_13*IT_0008;
    const ccomplex_t IT_0176 = s_12*s_24*s_34*IT_0145;
    const ccomplex_t IT_0177 = s_12*s_13;
    const ccomplex_t IT_0178 = s_23*IT_0145*IT_0168;
    const ccomplex_t IT_0179 = (-2)*IT_0178;
    const ccomplex_t IT_0180 = IT_0177 + IT_0179;
    const ccomplex_t IT_0181 = s_12*s_23;
    const ccomplex_t IT_0182 = s_13*IT_0145*IT_0172;
    const ccomplex_t IT_0183 = (-2)*IT_0182;
    const ccomplex_t IT_0184 = IT_0181 + IT_0183;
    const ccomplex_t IT_0185 = s_12*s_13*IT_0001*IT_0002;
    const ccomplex_t IT_0186 = -IT_0185;
    const ccomplex_t IT_0187 = s_13*s_14*s_24*IT_0145;
    const ccomplex_t IT_0188 = -IT_0187;
    const ccomplex_t IT_0189 = pow(s_14, 2);
    const ccomplex_t IT_0190 = s_23*IT_0145*IT_0189;
    const ccomplex_t IT_0191 = -IT_0190;
    const ccomplex_t IT_0192 = s_12*s_23*IT_0001*IT_0002;
    const ccomplex_t IT_0193 = -IT_0192;
    const ccomplex_t IT_0194 = pow(s_24, 2);
    const ccomplex_t IT_0195 = s_13*IT_0145*IT_0194;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = s_14*s_23*s_24*IT_0145;
    const ccomplex_t IT_0198 = -IT_0197;
    const ccomplex_t IT_0199 = IT_0026*IT_0156;
    const ccomplex_t IT_0200 = IT_0010*IT_0199;
    const ccomplex_t IT_0201 = IT_0163*conj(IT_0200);
    const ccomplex_t IT_0202 = IT_0165*IT_0200;
    const ccomplex_t IT_0203 = s_13*s_23*IT_0001*IT_0189;
    const ccomplex_t IT_0204 = s_14*s_24*IT_0001*IT_0168;
    const ccomplex_t IT_0205 = s_12*IT_0007;
    const ccomplex_t IT_0206 = s_14*s_24*IT_0007*IT_0145;
    const ccomplex_t IT_0207 = -IT_0206;
    const ccomplex_t IT_0208 = s_13*s_23*IT_0007*IT_0145;
    const ccomplex_t IT_0209 = -IT_0208;
    const ccomplex_t IT_0210 = s_12*s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0211 = -IT_0210;
    const ccomplex_t IT_0212 = IT_0203 + IT_0204 + IT_0205 + IT_0207 + IT_0209
       + IT_0211;
    const ccomplex_t IT_0213 = s_13*s_14*s_23*s_24*IT_0001;
    const ccomplex_t IT_0214 = IT_0001*IT_0168*IT_0194;
    const ccomplex_t IT_0215 = IT_0007*IT_0008;
    const ccomplex_t IT_0216 = IT_0008*IT_0145*IT_0168;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0007*IT_0145*IT_0194;
    const ccomplex_t IT_0219 = -IT_0218;
    const ccomplex_t IT_0220 = s_12*s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = IT_0213 + IT_0214 + IT_0215 + IT_0217 + IT_0219
       + IT_0221;
    const ccomplex_t IT_0223 = IT_0010*IT_0027;
    const ccomplex_t IT_0224 = -IT_0223;
    const ccomplex_t IT_0225 = pow(m_W, 2);
    const ccomplex_t IT_0226 = cpow((-2)*s_13 + IT_0007 + IT_0225 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0227 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0228 = IT_0013*IT_0227;
    const ccomplex_t IT_0229 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0230 = IT_0013*IT_0229;
    const ccomplex_t IT_0231 = 1.4142135623731*IT_0230;
    const ccomplex_t IT_0232 = (0 + _Complex_I*1)*(IT_0228 + 0.5*IT_0231);
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = conj(N_W2)*e_em*U_Wm1;
    const ccomplex_t IT_0235 = IT_0013*IT_0234;
    const ccomplex_t IT_0236 = conj(N_d2)*U_d1*e_em;
    const ccomplex_t IT_0237 = IT_0013*IT_0236;
    const ccomplex_t IT_0238 = 1.4142135623731*IT_0237;
    const ccomplex_t IT_0239 = (0 + _Complex_I*1)*(IT_0235 + 0.5*IT_0238);
    const ccomplex_t IT_0240 = -IT_0239;
    const ccomplex_t IT_0241 = IT_0233*IT_0240;
    const ccomplex_t IT_0242 = IT_0226*IT_0241;
    const ccomplex_t IT_0243 = (0 + _Complex_I*1)*IT_0242;
    const ccomplex_t IT_0244 = cpow((-2)*s_13 + IT_0007 + IT_0225 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0245 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0246 = IT_0013*IT_0245;
    const ccomplex_t IT_0247 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0248 = IT_0013*IT_0247;
    const ccomplex_t IT_0249 = 1.4142135623731*IT_0248;
    const ccomplex_t IT_0250 = (0 + _Complex_I*1)*(IT_0246 + 0.5*IT_0249);
    const ccomplex_t IT_0251 = -IT_0250;
    const ccomplex_t IT_0252 = conj(N_W2)*e_em*U_Wm2;
    const ccomplex_t IT_0253 = IT_0013*IT_0252;
    const ccomplex_t IT_0254 = conj(N_d2)*U_d2*e_em;
    const ccomplex_t IT_0255 = IT_0013*IT_0254;
    const ccomplex_t IT_0256 = 1.4142135623731*IT_0255;
    const ccomplex_t IT_0257 = (0 + _Complex_I*1)*(IT_0253 + 0.5*IT_0256);
    const ccomplex_t IT_0258 = -IT_0257;
    const ccomplex_t IT_0259 = IT_0251*IT_0258;
    const ccomplex_t IT_0260 = IT_0244*IT_0259;
    const ccomplex_t IT_0261 = (0 + _Complex_I*1)*IT_0260;
    const ccomplex_t IT_0262 = (-2)*IT_0243 + (-2)*IT_0261;
    const ccomplex_t IT_0263 = 2*IT_0262;
    const ccomplex_t IT_0264 = conj(IT_0224)*IT_0263;
    const ccomplex_t IT_0265 = 2*conj(IT_0262);
    const ccomplex_t IT_0266 = IT_0224*IT_0265;
    const ccomplex_t IT_0267 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0268 = IT_0013*IT_0267;
    const ccomplex_t IT_0269 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0270 = IT_0013*IT_0269;
    const ccomplex_t IT_0271 = 1.4142135623731*IT_0270;
    const ccomplex_t IT_0272 = (0 + _Complex_I*1)*(IT_0268 + (-0.5)*IT_0271);
    const ccomplex_t IT_0273 = -IT_0272;
    const ccomplex_t IT_0274 = N_W2*e_em*conj(V_Wp1);
    const ccomplex_t IT_0275 = IT_0013*IT_0274;
    const ccomplex_t IT_0276 = N_u2*conj(V_u1)*e_em;
    const ccomplex_t IT_0277 = IT_0013*IT_0276;
    const ccomplex_t IT_0278 = 1.4142135623731*IT_0277;
    const ccomplex_t IT_0279 = (0 + _Complex_I*1)*(IT_0275 + (-0.5)*IT_0278);
    const ccomplex_t IT_0280 = -IT_0279;
    const ccomplex_t IT_0281 = IT_0273*IT_0280;
    const ccomplex_t IT_0282 = IT_0226*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = N_W2*e_em*conj(V_Wp2);
    const ccomplex_t IT_0285 = IT_0013*IT_0284;
    const ccomplex_t IT_0286 = N_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0287 = IT_0013*IT_0286;
    const ccomplex_t IT_0288 = 1.4142135623731*IT_0287;
    const ccomplex_t IT_0289 = (0 + _Complex_I*1)*(IT_0285 + (-0.5)*IT_0288);
    const ccomplex_t IT_0290 = -IT_0289;
    const ccomplex_t IT_0291 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0292 = IT_0013*IT_0291;
    const ccomplex_t IT_0293 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0294 = IT_0013*IT_0293;
    const ccomplex_t IT_0295 = 1.4142135623731*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*(IT_0292 + (-0.5)*IT_0295);
    const ccomplex_t IT_0297 = -IT_0296;
    const ccomplex_t IT_0298 = IT_0290*IT_0297;
    const ccomplex_t IT_0299 = IT_0244*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0036*IT_0159;
    const ccomplex_t IT_0302 = IT_0010*IT_0301;
    const ccomplex_t IT_0303 = (-2)*IT_0283 + (-2)*IT_0300 + -IT_0302;
    const ccomplex_t IT_0304 = IT_0010*IT_0037;
    const ccomplex_t IT_0305 = -IT_0304;
    const ccomplex_t IT_0306 = 2*IT_0305;
    const ccomplex_t IT_0307 = conj(IT_0303)*IT_0306;
    const ccomplex_t IT_0308 = 2*conj(IT_0305);
    const ccomplex_t IT_0309 = IT_0303*IT_0308;
    const ccomplex_t IT_0310 = 2*IT_0224;
    const ccomplex_t IT_0311 = conj(IT_0162)*IT_0310;
    const ccomplex_t IT_0312 = 2*conj(IT_0224);
    const ccomplex_t IT_0313 = IT_0162*IT_0312;
    const ccomplex_t IT_0314 = IT_0001*IT_0172*IT_0189;
    const ccomplex_t IT_0315 = IT_0008*IT_0145*IT_0189;
    const ccomplex_t IT_0316 = -IT_0315;
    const ccomplex_t IT_0317 = IT_0007*IT_0145*IT_0172;
    const ccomplex_t IT_0318 = -IT_0317;
    const ccomplex_t IT_0319 = s_12*s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0320 = -IT_0319;
    const ccomplex_t IT_0321 = IT_0213 + IT_0215 + IT_0314 + IT_0316 + IT_0318
       + IT_0320;
    const ccomplex_t IT_0322 = s_14*s_24*IT_0001*IT_0172;
    const ccomplex_t IT_0323 = s_13*s_23*IT_0001*IT_0194;
    const ccomplex_t IT_0324 = s_12*IT_0008;
    const ccomplex_t IT_0325 = s_14*s_24*IT_0008*IT_0145;
    const ccomplex_t IT_0326 = -IT_0325;
    const ccomplex_t IT_0327 = s_13*s_23*IT_0008*IT_0145;
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = s_12*s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0330 = -IT_0329;
    const ccomplex_t IT_0331 = IT_0322 + IT_0323 + IT_0324 + IT_0326 + IT_0328
       + IT_0330;
    const ccomplex_t IT_0332 = cpow((-2)*s_23 + IT_0008 + IT_0225 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0333 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0334 = IT_0013*IT_0333;
    const ccomplex_t IT_0335 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0336 = IT_0013*IT_0335;
    const ccomplex_t IT_0337 = 1.4142135623731*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*(IT_0334 + (-0.5)*IT_0337);
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = conj(N_W2)*e_em*V_Wp1;
    const ccomplex_t IT_0341 = IT_0013*IT_0340;
    const ccomplex_t IT_0342 = conj(N_u2)*V_u1*e_em;
    const ccomplex_t IT_0343 = IT_0013*IT_0342;
    const ccomplex_t IT_0344 = 1.4142135623731*IT_0343;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*(IT_0341 + (-0.5)*IT_0344);
    const ccomplex_t IT_0346 = -IT_0345;
    const ccomplex_t IT_0347 = IT_0339*IT_0346;
    const ccomplex_t IT_0348 = IT_0332*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = cpow((-2)*s_23 + IT_0008 + IT_0225 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0351 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0352 = IT_0013*IT_0351;
    const ccomplex_t IT_0353 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0354 = IT_0013*IT_0353;
    const ccomplex_t IT_0355 = 1.4142135623731*IT_0354;
    const ccomplex_t IT_0356 = (0 + _Complex_I*1)*(IT_0352 + (-0.5)*IT_0355);
    const ccomplex_t IT_0357 = -IT_0356;
    const ccomplex_t IT_0358 = conj(N_W2)*e_em*V_Wp2;
    const ccomplex_t IT_0359 = IT_0013*IT_0358;
    const ccomplex_t IT_0360 = conj(N_u2)*V_u2*e_em;
    const ccomplex_t IT_0361 = IT_0013*IT_0360;
    const ccomplex_t IT_0362 = 1.4142135623731*IT_0361;
    const ccomplex_t IT_0363 = (0 + _Complex_I*1)*(IT_0359 + (-0.5)*IT_0362);
    const ccomplex_t IT_0364 = -IT_0363;
    const ccomplex_t IT_0365 = IT_0357*IT_0364;
    const ccomplex_t IT_0366 = IT_0350*IT_0365;
    const ccomplex_t IT_0367 = (0 + _Complex_I*1)*IT_0366;
    const ccomplex_t IT_0368 = 2*IT_0349 + 2*IT_0367;
    const ccomplex_t IT_0369 = 2*IT_0368;
    const ccomplex_t IT_0370 = conj(IT_0224)*IT_0369;
    const ccomplex_t IT_0371 = 2*conj(IT_0368);
    const ccomplex_t IT_0372 = IT_0224*IT_0371;
    const ccomplex_t IT_0373 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0374 = IT_0013*IT_0373;
    const ccomplex_t IT_0375 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0376 = IT_0013*IT_0375;
    const ccomplex_t IT_0377 = 1.4142135623731*IT_0376;
    const ccomplex_t IT_0378 = (0 + _Complex_I*1)*(IT_0374 + 0.5*IT_0377);
    const ccomplex_t IT_0379 = -IT_0378;
    const ccomplex_t IT_0380 = N_W2*e_em*conj(U_Wm1);
    const ccomplex_t IT_0381 = IT_0013*IT_0380;
    const ccomplex_t IT_0382 = N_d2*conj(U_d1)*e_em;
    const ccomplex_t IT_0383 = IT_0013*IT_0382;
    const ccomplex_t IT_0384 = 1.4142135623731*IT_0383;
    const ccomplex_t IT_0385 = (0 + _Complex_I*1)*(IT_0381 + 0.5*IT_0384);
    const ccomplex_t IT_0386 = -IT_0385;
    const ccomplex_t IT_0387 = IT_0379*IT_0386;
    const ccomplex_t IT_0388 = IT_0332*IT_0387;
    const ccomplex_t IT_0389 = (0 + _Complex_I*1)*IT_0388;
    const ccomplex_t IT_0390 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0391 = IT_0013*IT_0390;
    const ccomplex_t IT_0392 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0393 = IT_0013*IT_0392;
    const ccomplex_t IT_0394 = 1.4142135623731*IT_0393;
    const ccomplex_t IT_0395 = (0 + _Complex_I*1)*(IT_0391 + 0.5*IT_0394);
    const ccomplex_t IT_0396 = -IT_0395;
    const ccomplex_t IT_0397 = N_W2*e_em*conj(U_Wm2);
    const ccomplex_t IT_0398 = IT_0013*IT_0397;
    const ccomplex_t IT_0399 = N_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0400 = IT_0013*IT_0399;
    const ccomplex_t IT_0401 = 1.4142135623731*IT_0400;
    const ccomplex_t IT_0402 = (0 + _Complex_I*1)*(IT_0398 + 0.5*IT_0401);
    const ccomplex_t IT_0403 = -IT_0402;
    const ccomplex_t IT_0404 = IT_0396*IT_0403;
    const ccomplex_t IT_0405 = IT_0350*IT_0404;
    const ccomplex_t IT_0406 = (0 + _Complex_I*1)*IT_0405;
    const ccomplex_t IT_0407 = -IT_0302 + 2*IT_0389 + 2*IT_0406;
    const ccomplex_t IT_0408 = IT_0306*conj(IT_0407);
    const ccomplex_t IT_0409 = IT_0308*IT_0407;
    const ccomplex_t IT_0410 = s_14*s_34*IT_0145;
    const ccomplex_t IT_0411 = s_12*IT_0225;
    const ccomplex_t IT_0412 = s_13*s_23;
    const ccomplex_t IT_0413 = (-2)*IT_0412;
    const ccomplex_t IT_0414 = IT_0411 + IT_0413;
    const ccomplex_t IT_0415 = IT_0145*IT_0414;
    const ccomplex_t IT_0416 = (s_13 + -IT_0410)*(s_12 + (-0.5)*IT_0415);
    const ccomplex_t IT_0417 = (-2)*IT_0416;
    const ccomplex_t IT_0418 = s_24*s_34*IT_0145;
    const ccomplex_t IT_0419 = (s_12 + (-0.5)*IT_0415)*(s_23 + -IT_0418);
    const ccomplex_t IT_0420 = (-2)*IT_0419;
    const ccomplex_t IT_0421 = IT_0417 + IT_0420;
    const ccomplex_t IT_0422 = IT_0014*IT_0026;
    const ccomplex_t IT_0423 = IT_0010*IT_0422;
    const ccomplex_t IT_0424 = IT_0310*conj(IT_0423);
    const ccomplex_t IT_0425 = IT_0312*IT_0423;
    const ccomplex_t IT_0426 = 2*IT_0167;
    const ccomplex_t IT_0427 = IT_0170 + IT_0180 + IT_0188 + IT_0191 + IT_0426;
    const ccomplex_t IT_0428 = 2*IT_0174;
    const ccomplex_t IT_0429 = IT_0175 + IT_0184 + IT_0196 + IT_0198 + IT_0428;
    const ccomplex_t IT_0430 = conj(IT_0158)*IT_0306;
    const ccomplex_t IT_0431 = IT_0158*IT_0308;
    const ccomplex_t IT_0432 = conj(IT_0200)*IT_0310;
    const ccomplex_t IT_0433 = IT_0200*IT_0312;
    const ccomplex_t IT_0434 = s_13*IT_0225;
    const ccomplex_t IT_0435 = s_14*s_34;
    const ccomplex_t IT_0436 = -IT_0435;
    const ccomplex_t IT_0437 = IT_0434 + IT_0436;
    const ccomplex_t IT_0438 = IT_0145*IT_0437;
    const ccomplex_t IT_0439 = -IT_0410 + (-0.25)*IT_0438;
    const ccomplex_t IT_0440 = s_13 + IT_0439;
    const ccomplex_t IT_0441 = IT_0000*IT_0440;
    const ccomplex_t IT_0442 = (-4)*IT_0441;
    const ccomplex_t IT_0443 = s_23*IT_0225;
    const ccomplex_t IT_0444 = s_24*s_34;
    const ccomplex_t IT_0445 = -IT_0444;
    const ccomplex_t IT_0446 = IT_0443 + IT_0445;
    const ccomplex_t IT_0447 = IT_0145*IT_0446;
    const ccomplex_t IT_0448 = -IT_0418 + (-0.25)*IT_0447;
    const ccomplex_t IT_0449 = s_23 + IT_0448;
    const ccomplex_t IT_0450 = IT_0000*IT_0449;
    const ccomplex_t IT_0451 = (-4)*IT_0450;
    const ccomplex_t IT_0452 = -IT_0349 + -IT_0367;
    const ccomplex_t IT_0453 = IT_0306*conj(IT_0452);
    const ccomplex_t IT_0454 = IT_0243 + IT_0261;
    const ccomplex_t IT_0455 = IT_0306*conj(IT_0454);
    const ccomplex_t IT_0456 = IT_0308*IT_0452;
    const ccomplex_t IT_0457 = IT_0308*IT_0454;
    const ccomplex_t IT_0458 = IT_0283 + IT_0300;
    const ccomplex_t IT_0459 = IT_0310*conj(IT_0458);
    const ccomplex_t IT_0460 = -IT_0389 + -IT_0406;
    const ccomplex_t IT_0461 = IT_0310*conj(IT_0460);
    const ccomplex_t IT_0462 = IT_0312*IT_0458;
    const ccomplex_t IT_0463 = IT_0312*IT_0460;
    const ccomplex_t IT_0464 = s_24*s_34*IT_0007*IT_0145;
    const ccomplex_t IT_0465 = -IT_0464;
    const ccomplex_t IT_0466 = IT_0167 + IT_0169 + IT_0170 + IT_0180 + IT_0186
       + IT_0465;
    const ccomplex_t IT_0467 = s_14*s_34*IT_0008*IT_0145;
    const ccomplex_t IT_0468 = -IT_0467;
    const ccomplex_t IT_0469 = IT_0173 + IT_0174 + IT_0175 + IT_0184 + IT_0193
       + IT_0468;
    const ccomplex_t IT_0470 = 2*IT_0423;
    const ccomplex_t IT_0471 = conj(IT_0162)*IT_0470;
    const ccomplex_t IT_0472 = 2*conj(IT_0423);
    const ccomplex_t IT_0473 = IT_0162*IT_0472;
    const ccomplex_t IT_0474 = IT_0145*IT_0189;
    const ccomplex_t IT_0475 = -IT_0007;
    const ccomplex_t IT_0476 = IT_0474 + IT_0475;
    const ccomplex_t IT_0477 = IT_0000*IT_0476;
    const ccomplex_t IT_0478 = IT_0145*IT_0194;
    const ccomplex_t IT_0479 = -IT_0008;
    const ccomplex_t IT_0480 = IT_0478 + IT_0479;
    const ccomplex_t IT_0481 = IT_0000*IT_0480;
    const ccomplex_t IT_0482 = conj(IT_0224)*IT_0305;
    const ccomplex_t IT_0483 = 6*IT_0482;
    const ccomplex_t IT_0484 = IT_0224*conj(IT_0305);
    const ccomplex_t IT_0485 = 6*IT_0484;
    const ccomplex_t IT_0486 = -IT_0410;
    const ccomplex_t IT_0487 = s_13 + IT_0486;
    const ccomplex_t IT_0488 = IT_0000*IT_0487;
    const ccomplex_t IT_0489 = -IT_0418;
    const ccomplex_t IT_0490 = s_23 + IT_0489;
    const ccomplex_t IT_0491 = IT_0000*IT_0490;
    const ccomplex_t IT_0492 = IT_0014*IT_0036;
    const ccomplex_t IT_0493 = IT_0010*IT_0492;
    const ccomplex_t IT_0494 = (-6)*IT_0493;
    const ccomplex_t IT_0495 = conj(IT_0224)*IT_0494;
    const ccomplex_t IT_0496 = (-6)*conj(IT_0493);
    const ccomplex_t IT_0497 = IT_0224*IT_0496;
    const ccomplex_t IT_0498 = IT_0305*conj(IT_0423);
    const ccomplex_t IT_0499 = (-6)*IT_0498;
    const ccomplex_t IT_0500 = conj(IT_0305)*IT_0423;
    const ccomplex_t IT_0501 = (-6)*IT_0500;
    const ccomplex_t IT_0502 = s_12*IT_0005;
    const ccomplex_t IT_0503 = (-4)*conj(IT_0139);
    const ccomplex_t IT_0504 = conj(IT_0041) + 0.25*IT_0503;
    const ccomplex_t IT_0505 = (-4)*IT_0139;
    const ccomplex_t IT_0506 = s_13*s_24*s_34*IT_0001;
    const ccomplex_t IT_0507 = s_13*s_23*IT_0145;
    const ccomplex_t IT_0508 = s_14*s_24*IT_0145;
    const ccomplex_t IT_0509 = -IT_0507 + -IT_0508;
    const ccomplex_t IT_0510 = s_12 + IT_0506 + IT_0509;
    const ccomplex_t IT_0511 = m_N_2*IT_0510;
    const ccomplex_t IT_0512 = 2*IT_0098;
    const ccomplex_t IT_0513 = (-2)*IT_0094;
    const ccomplex_t IT_0514 = IT_0512 + IT_0513;
    const ccomplex_t IT_0515 = IT_0145*IT_0168;
    const ccomplex_t IT_0516 = IT_0474 + IT_0515;
    const ccomplex_t IT_0517 = s_13*s_14*s_34*IT_0001;
    const ccomplex_t IT_0518 = -IT_0007 + -IT_0517;
    const ccomplex_t IT_0519 = IT_0516 + IT_0518;
    const ccomplex_t IT_0520 = m_N_3*IT_0519;
    const ccomplex_t IT_0521 = 2*IT_0041;
    const ccomplex_t IT_0522 = (-2)*IT_0139;
    const ccomplex_t IT_0523 = IT_0521 + IT_0522;
    const ccomplex_t IT_0524 = s_14*s_24;
    const ccomplex_t IT_0525 = (-2)*IT_0524;
    const ccomplex_t IT_0526 = IT_0411 + IT_0525;
    const ccomplex_t IT_0527 = IT_0145*IT_0526;
    const ccomplex_t IT_0528 = -(IT_0007 + -IT_0515)*(s_12 + (-0.5)*IT_0527);
    const ccomplex_t IT_0529 = (-2)*IT_0528;
    const ccomplex_t IT_0530 = (-2)*IT_0262;
    const ccomplex_t IT_0531 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0532 = IT_0233*IT_0280;
    const ccomplex_t IT_0533 = IT_0226*IT_0531*IT_0532;
    const ccomplex_t IT_0534 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0535 = IT_0251*IT_0290;
    const ccomplex_t IT_0536 = IT_0244*IT_0534*IT_0535;
    const ccomplex_t IT_0537 = m_N_2*IT_0243;
    const ccomplex_t IT_0538 = m_N_2*IT_0261;
    const ccomplex_t IT_0539 = IT_0533 + IT_0536 + -IT_0537 + -IT_0538;
    const ccomplex_t IT_0540 = IT_0145*IT_0411;
    const ccomplex_t IT_0541 = (-0.5)*IT_0540;
    const ccomplex_t IT_0542 = s_34*IT_0001;
    const ccomplex_t IT_0543 = s_14*s_23;
    const ccomplex_t IT_0544 = s_13*s_24;
    const ccomplex_t IT_0545 = s_12*s_34;
    const ccomplex_t IT_0546 = -IT_0544 + -IT_0545;
    const ccomplex_t IT_0547 = IT_0543 + IT_0546;
    const ccomplex_t IT_0548 = IT_0542*IT_0547;
    const ccomplex_t IT_0549 = (-0.25)*IT_0548;
    const ccomplex_t IT_0550 = s_12 + IT_0541 + IT_0549;
    const ccomplex_t IT_0551 = -conj(IT_0098);
    const ccomplex_t IT_0552 = conj(IT_0094) + IT_0551;
    const ccomplex_t IT_0553 = (-4)*IT_0006;
    const ccomplex_t IT_0554 = -conj(IT_0041);
    const ccomplex_t IT_0555 = conj(IT_0139) + IT_0554;
    const ccomplex_t IT_0556 = IT_0145*IT_0172;
    const ccomplex_t IT_0557 = s_23*s_24*s_34*IT_0145;
    const ccomplex_t IT_0558 = IT_0008*IT_0225;
    const ccomplex_t IT_0559 = 2*IT_0194 + (-2)*IT_0557 + -IT_0558;
    const ccomplex_t IT_0560 = IT_0172 + IT_0559;
    const ccomplex_t IT_0561 = IT_0145*IT_0560;
    const ccomplex_t IT_0562 = -IT_0008 + (-0.5)*IT_0561;
    const ccomplex_t IT_0563 = IT_0556 + IT_0562;
    const ccomplex_t IT_0564 = m_N_2*IT_0563;
    const ccomplex_t IT_0565 = (-2)*IT_0564;
    const ccomplex_t IT_0566 = s_13*s_24*s_34*IT_0145;
    const ccomplex_t IT_0567 = -IT_0412 + (-2)*IT_0524 + 2*IT_0566;
    const ccomplex_t IT_0568 = IT_0411 + IT_0567;
    const ccomplex_t IT_0569 = IT_0145*IT_0568;
    const ccomplex_t IT_0570 = -IT_0507 + (-0.5)*IT_0569;
    const ccomplex_t IT_0571 = s_12 + IT_0570;
    const ccomplex_t IT_0572 = m_N_2*IT_0571;
    const ccomplex_t IT_0573 = (-2)*IT_0572;
    const ccomplex_t IT_0574 = (-2)*conj(IT_0262);
    const ccomplex_t IT_0575 = IT_0349 + IT_0367;
    const ccomplex_t IT_0576 = m_N_3*IT_0575;
    const ccomplex_t IT_0577 = -IT_0576;
    const ccomplex_t IT_0578 = IT_0543 + IT_0545;
    const ccomplex_t IT_0579 = -IT_0544;
    const ccomplex_t IT_0580 = IT_0578 + IT_0579;
    const ccomplex_t IT_0581 = IT_0542*IT_0580;
    const ccomplex_t IT_0582 = 0.25*IT_0581;
    const ccomplex_t IT_0583 = s_12 + IT_0541 + IT_0582;
    const ccomplex_t IT_0584 = s_14*s_23*s_34*IT_0145;
    const ccomplex_t IT_0585 = -IT_0412 + (-2)*IT_0524 + 2*IT_0584;
    const ccomplex_t IT_0586 = IT_0411 + IT_0585;
    const ccomplex_t IT_0587 = IT_0145*IT_0586;
    const ccomplex_t IT_0588 = -IT_0507 + (-0.5)*IT_0587;
    const ccomplex_t IT_0589 = s_12 + IT_0588;
    const ccomplex_t IT_0590 = m_N_3*IT_0589;
    const ccomplex_t IT_0591 = (-2)*IT_0590;
    const ccomplex_t IT_0592 = s_13*s_14*s_34*IT_0145;
    const ccomplex_t IT_0593 = IT_0007*IT_0225;
    const ccomplex_t IT_0594 = 2*IT_0189 + (-2)*IT_0592 + -IT_0593;
    const ccomplex_t IT_0595 = IT_0168 + IT_0594;
    const ccomplex_t IT_0596 = IT_0145*IT_0595;
    const ccomplex_t IT_0597 = -IT_0007 + (-0.5)*IT_0596;
    const ccomplex_t IT_0598 = IT_0515 + IT_0597;
    const ccomplex_t IT_0599 = m_N_3*IT_0598;
    const ccomplex_t IT_0600 = (-2)*IT_0599;
    const ccomplex_t IT_0601 = s_34*m_N_3*IT_0145;
    const ccomplex_t IT_0602 = IT_0149*IT_0601;
    const ccomplex_t IT_0603 = s_34*m_N_2*IT_0145;
    const ccomplex_t IT_0604 = IT_0153*IT_0603;
    const ccomplex_t IT_0605 = 2*IT_0139;
    const ccomplex_t IT_0606 = (-2)*IT_0041;
    const ccomplex_t IT_0607 = IT_0605 + IT_0606;
    const ccomplex_t IT_0608 = IT_0002*IT_0145;
    const ccomplex_t IT_0609 = (s_12 + (-0.5)*IT_0415)*(IT_0225 + -IT_0608);
    const ccomplex_t IT_0610 = (-2)*IT_0609;
    const ccomplex_t IT_0611 = 2*IT_0493;
    const ccomplex_t IT_0612 = IT_0443 + IT_0444;
    const ccomplex_t IT_0613 = s_23*IT_0002*IT_0145;
    const ccomplex_t IT_0614 = (-2)*IT_0613;
    const ccomplex_t IT_0615 = IT_0612 + IT_0614;
    const ccomplex_t IT_0616 = IT_0145*IT_0615;
    const ccomplex_t IT_0617 = -IT_0418 + 0.5*IT_0616;
    const ccomplex_t IT_0618 = s_23 + IT_0617;
    const ccomplex_t IT_0619 = m_N_2*IT_0618;
    const ccomplex_t IT_0620 = 2*IT_0619;
    const ccomplex_t IT_0621 = 2*IT_0577;
    const ccomplex_t IT_0622 = IT_0434 + IT_0435;
    const ccomplex_t IT_0623 = s_13*IT_0002*IT_0145;
    const ccomplex_t IT_0624 = (-2)*IT_0623;
    const ccomplex_t IT_0625 = IT_0622 + IT_0624;
    const ccomplex_t IT_0626 = IT_0145*IT_0625;
    const ccomplex_t IT_0627 = -IT_0410 + 0.5*IT_0626;
    const ccomplex_t IT_0628 = s_13 + IT_0627;
    const ccomplex_t IT_0629 = m_N_3*IT_0628;
    const ccomplex_t IT_0630 = 2*IT_0629;
    const ccomplex_t IT_0631 = (-2)*IT_0539;
    const ccomplex_t IT_0632 = 2*conj(IT_0098);
    const ccomplex_t IT_0633 = (-2)*conj(IT_0094);
    const ccomplex_t IT_0634 = IT_0632 + IT_0633;
    const ccomplex_t IT_0635 = 2*conj(IT_0139);
    const ccomplex_t IT_0636 = (-2)*conj(IT_0041);
    const ccomplex_t IT_0637 = IT_0635 + IT_0636;
    const ccomplex_t IT_0638 = 2*conj(IT_0577);
    const ccomplex_t IT_0639 = (-2)*conj(IT_0539);
    const ccomplex_t IT_0640 = IT_0478 + IT_0556;
    const ccomplex_t IT_0641 = s_23*s_24*s_34*IT_0001;
    const ccomplex_t IT_0642 = -IT_0008 + -IT_0641;
    const ccomplex_t IT_0643 = IT_0640 + IT_0642;
    const ccomplex_t IT_0644 = m_N_2*IT_0643;
    const ccomplex_t IT_0645 = 2*IT_0094;
    const ccomplex_t IT_0646 = (-2)*IT_0098;
    const ccomplex_t IT_0647 = IT_0645 + IT_0646;
    const ccomplex_t IT_0648 = s_14*s_23*s_34*IT_0001;
    const ccomplex_t IT_0649 = s_12 + IT_0509 + IT_0648;
    const ccomplex_t IT_0650 = m_N_3*IT_0649;
    const ccomplex_t IT_0651 = (s_12 + -IT_0507)*(s_12 + (-0.5)*IT_0527);
    const ccomplex_t IT_0652 = (-2)*IT_0651;
    const ccomplex_t IT_0653 = -(s_12 + (-0.5)*IT_0527)*(IT_0008 + -IT_0556);
    const ccomplex_t IT_0654 = (-2)*IT_0653;
    const ccomplex_t IT_0655 = (-2)*IT_0368;
    const ccomplex_t IT_0656 = -(s_12 + (-0.5)*IT_0415)*(IT_0007 + -IT_0474);
    const ccomplex_t IT_0657 = (-2)*IT_0656;
    const ccomplex_t IT_0658 = -(s_12 + (-0.5)*IT_0415)*(IT_0008 + -IT_0478);
    const ccomplex_t IT_0659 = (-2)*IT_0658;
    const ccomplex_t IT_0660 = IT_0657 + IT_0659;
    const ccomplex_t IT_0661 = (-2)*IT_0305;
    const ccomplex_t IT_0662 = conj(IT_0305)*IT_0661;
    const ccomplex_t IT_0663 = IT_0000*IT_0510;
    const ccomplex_t IT_0664 = IT_0000*IT_0649;
    const ccomplex_t IT_0665 = conj(IT_0162)*IT_0661;
    const ccomplex_t IT_0666 = (-2)*conj(IT_0305);
    const ccomplex_t IT_0667 = IT_0162*IT_0666;
    const ccomplex_t IT_0668 = (s_12 + (-0.5)*IT_0415)*(s_12 + -IT_0508);
    const ccomplex_t IT_0669 = (-2)*IT_0668;
    const ccomplex_t IT_0670 = -IT_0508;
    const ccomplex_t IT_0671 = s_12 + IT_0670;
    const ccomplex_t IT_0672 = IT_0000*IT_0671;
    const ccomplex_t IT_0673 = m_N_2*IT_0649;
    const ccomplex_t IT_0674 = m_N_3*IT_0510;
    const ccomplex_t IT_0675 = -IT_0524;
    const ccomplex_t IT_0676 = IT_0411 + IT_0675;
    const ccomplex_t IT_0677 = IT_0145*IT_0676;
    const ccomplex_t IT_0678 = -IT_0508 + (-0.25)*IT_0677;
    const ccomplex_t IT_0679 = s_12 + IT_0678;
    const ccomplex_t IT_0680 = m_N_2*IT_0679;
    const ccomplex_t IT_0681 = (-4)*IT_0680;
    const ccomplex_t IT_0682 = 2*IT_0539;
    const ccomplex_t IT_0683 = IT_0000*IT_0519;
    const ccomplex_t IT_0684 = 2*IT_0303;
    const ccomplex_t IT_0685 = IT_0000*IT_0643;
    const ccomplex_t IT_0686 = 2*IT_0407;
    const ccomplex_t IT_0687 = -IT_0593;
    const ccomplex_t IT_0688 = IT_0189 + IT_0687;
    const ccomplex_t IT_0689 = IT_0145*IT_0688;
    const ccomplex_t IT_0690 = (-4)*IT_0476;
    const ccomplex_t IT_0691 = IT_0689 + IT_0690;
    const ccomplex_t IT_0692 = m_N_3*IT_0691;
    const ccomplex_t IT_0693 = -IT_0558;
    const ccomplex_t IT_0694 = IT_0194 + IT_0693;
    const ccomplex_t IT_0695 = IT_0145*IT_0694;
    const ccomplex_t IT_0696 = (-4)*IT_0480;
    const ccomplex_t IT_0697 = IT_0695 + IT_0696;
    const ccomplex_t IT_0698 = m_N_2*IT_0697;
    const ccomplex_t IT_0699 = (-2)*IT_0303;
    const ccomplex_t IT_0700 = (-2)*IT_0407;
    const ccomplex_t IT_0701 = (-2)*IT_0224;
    const ccomplex_t IT_0702 = m_N_3*IT_0679;
    const ccomplex_t IT_0703 = (-4)*IT_0702;
    const ccomplex_t IT_0704 = (-2)*IT_0577;
    const ccomplex_t IT_0705 = 0.25*IT_0520*IT_0523 + 0.25*IT_0644*IT_0647 +
       IT_0224*IT_0669 + (-3)*IT_0305*IT_0672 + 0.25*IT_0514*IT_0673 + 0.25
      *IT_0607*IT_0674 + 0.25*IT_0681*IT_0682 + 0.25*IT_0683*IT_0684 + 0.25
      *IT_0685*IT_0686 + 0.25*IT_0621*IT_0692 + 0.25*IT_0631*IT_0698 + 0.25
      *IT_0663*IT_0699 + 0.25*IT_0664*IT_0700 + 0.25*IT_0660*IT_0701 + 0.25
      *IT_0703*IT_0704;
    const ccomplex_t IT_0706 = m_N_2*IT_0449;
    const ccomplex_t IT_0707 = (-4)*IT_0706;
    const ccomplex_t IT_0708 = -IT_0608;
    const ccomplex_t IT_0709 = IT_0225 + IT_0708;
    const ccomplex_t IT_0710 = IT_0000*IT_0709;
    const ccomplex_t IT_0711 = m_N_3*IT_0440;
    const ccomplex_t IT_0712 = (-4)*IT_0711;
    const ccomplex_t IT_0713 = conj(IT_0423)*(IT_0263*IT_0466 + IT_0369
      *IT_0469 + IT_0523*IT_0602 + IT_0470*IT_0610 + IT_0604*IT_0647 + IT_0150
      *IT_0684 + IT_0154*IT_0686 + IT_0682*IT_0707 + IT_0494*IT_0710 + IT_0704
      *IT_0712);
    const ccomplex_t IT_0714 = 2*conj(IT_0094);
    const ccomplex_t IT_0715 = (-2)*conj(IT_0098);
    const ccomplex_t IT_0716 = IT_0714 + IT_0715;
    const ccomplex_t IT_0717 = 2*conj(IT_0041);
    const ccomplex_t IT_0718 = (-2)*conj(IT_0139);
    const ccomplex_t IT_0719 = IT_0717 + IT_0718;
    const ccomplex_t IT_0720 = 2*conj(IT_0539);
    const ccomplex_t IT_0721 = 2*conj(IT_0303);
    const ccomplex_t IT_0722 = 2*conj(IT_0407);
    const ccomplex_t IT_0723 = (-2)*conj(IT_0577);
    const ccomplex_t IT_0724 = IT_0339*IT_0386;
    const ccomplex_t IT_0725 = IT_0332*IT_0531*IT_0724;
    const ccomplex_t IT_0726 = IT_0357*IT_0403;
    const ccomplex_t IT_0727 = IT_0350*IT_0534*IT_0726;
    const ccomplex_t IT_0728 = m_N_3*IT_0389;
    const ccomplex_t IT_0729 = m_N_3*IT_0406;
    const ccomplex_t IT_0730 = -IT_0725 + -IT_0727 + IT_0728 + IT_0729;
    const ccomplex_t IT_0731 = 1.125*IT_0000;
    const ccomplex_t IT_0732 = IT_0172 + IT_0693;
    const ccomplex_t IT_0733 = IT_0145*IT_0732;
    const ccomplex_t IT_0734 = IT_0479 + IT_0556;
    const ccomplex_t IT_0735 = (-4)*IT_0734;
    const ccomplex_t IT_0736 = IT_0733 + IT_0735;
    const ccomplex_t IT_0737 = m_N_2*IT_0736;
    const ccomplex_t IT_0738 = IT_0165 + IT_0371;
    const ccomplex_t IT_0739 = -IT_0412;
    const ccomplex_t IT_0740 = IT_0411 + IT_0739;
    const ccomplex_t IT_0741 = IT_0145*IT_0740;
    const ccomplex_t IT_0742 = -IT_0507 + (-0.25)*IT_0741;
    const ccomplex_t IT_0743 = s_12 + IT_0742;
    const ccomplex_t IT_0744 = m_N_2*IT_0743;
    const ccomplex_t IT_0745 = (-4)*IT_0744;
    const ccomplex_t IT_0746 = (-2)*conj(IT_0162);
    const ccomplex_t IT_0747 = IT_0574 + IT_0746;
    const ccomplex_t IT_0748 = 2*conj(IT_0493);
    const ccomplex_t IT_0749 = (-2)*IT_0412 + -IT_0524 + 2*IT_0584;
    const ccomplex_t IT_0750 = IT_0411 + IT_0749;
    const ccomplex_t IT_0751 = IT_0145*IT_0750;
    const ccomplex_t IT_0752 = -IT_0508 + (-0.5)*IT_0751;
    const ccomplex_t IT_0753 = s_12 + IT_0752;
    const ccomplex_t IT_0754 = m_N_2*IT_0753;
    const ccomplex_t IT_0755 = (-2)*IT_0754;
    const ccomplex_t IT_0756 = pow(m_W, 4);
    const ccomplex_t IT_0757 = s_12*IT_0756;
    const ccomplex_t IT_0758 = s_13*s_23*IT_0225;
    const ccomplex_t IT_0759 = s_14*s_24*IT_0225;
    const ccomplex_t IT_0760 = s_13*s_24*s_34;
    const ccomplex_t IT_0761 = s_14*s_23*s_34;
    const ccomplex_t IT_0762 = s_12*IT_0002;
    const ccomplex_t IT_0763 = (-2)*IT_0758 + (-2)*IT_0759 + 2*IT_0760 + 2
      *IT_0761 + (-2)*IT_0762;
    const ccomplex_t IT_0764 = IT_0757 + IT_0763;
    const ccomplex_t IT_0765 = IT_0001*IT_0764;
    const ccomplex_t IT_0766 = (-0.25)*IT_0765;
    const ccomplex_t IT_0767 = -IT_0507;
    const ccomplex_t IT_0768 = s_12 + IT_0670 + IT_0766 + IT_0767;
    const ccomplex_t IT_0769 = (-8)*conj(IT_0539);
    const ccomplex_t IT_0770 = 0.5*IT_0194 + -IT_0557 + (-0.5)*IT_0558;
    const ccomplex_t IT_0771 = IT_0172 + IT_0770;
    const ccomplex_t IT_0772 = IT_0145*IT_0771;
    const ccomplex_t IT_0773 = -IT_0478;
    const ccomplex_t IT_0774 = IT_0008 + IT_0772 + IT_0773;
    const ccomplex_t IT_0775 = m_N_2*IT_0774;
    const ccomplex_t IT_0776 = (-4)*conj(IT_0224);
    const ccomplex_t IT_0777 = (-2)*conj(IT_0407);
    const ccomplex_t IT_0778 = (-2)*conj(IT_0423);
    const ccomplex_t IT_0779 = s_14*IT_0225;
    const ccomplex_t IT_0780 = IT_0542*IT_0779;
    const ccomplex_t IT_0781 = IT_0145*IT_0434;
    const ccomplex_t IT_0782 = 0.25*IT_0780 + (-0.5)*IT_0781;
    const ccomplex_t IT_0783 = s_13 + IT_0782;
    const ccomplex_t IT_0784 = m_N_3*IT_0783;
    const ccomplex_t IT_0785 = 4*IT_0784;
    const ccomplex_t IT_0786 = (-2)*conj(IT_0158);
    const ccomplex_t IT_0787 = IT_0145*IT_0443;
    const ccomplex_t IT_0788 = s_24*IT_0225;
    const ccomplex_t IT_0789 = s_23*s_34;
    const ccomplex_t IT_0790 = (-2)*IT_0789;
    const ccomplex_t IT_0791 = IT_0788 + IT_0790;
    const ccomplex_t IT_0792 = IT_0542*IT_0791;
    const ccomplex_t IT_0793 = (-0.5)*IT_0787 + (-0.25)*IT_0792;
    const ccomplex_t IT_0794 = s_23 + IT_0793;
    const ccomplex_t IT_0795 = m_N_2*IT_0794;
    const ccomplex_t IT_0796 = (-4)*IT_0795;
    const ccomplex_t IT_0797 = (-2)*conj(IT_0200);
    const ccomplex_t IT_0798 = -IT_0553*IT_0555 + (-8)*IT_0552*IT_0583 +
       IT_0666*IT_0692 + IT_0308*IT_0703 + IT_0600*IT_0721 + 16*conj(IT_0577)
      *IT_0731 + IT_0737*IT_0738 + IT_0745*IT_0747 + IT_0712*IT_0748 + IT_0312
      *IT_0755 + IT_0768*IT_0769 + IT_0775*IT_0776 + IT_0591*IT_0777 + IT_0620
      *IT_0778 + IT_0785*IT_0786 + IT_0796*IT_0797;
    const ccomplex_t IT_0799 = IT_0346*IT_0379;
    const ccomplex_t IT_0800 = IT_0332*IT_0531*IT_0799;
    const ccomplex_t IT_0801 = IT_0364*IT_0396;
    const ccomplex_t IT_0802 = IT_0350*IT_0534*IT_0801;
    const ccomplex_t IT_0803 = -IT_0800 + -IT_0802;
    const ccomplex_t IT_0804 = -IT_0041;
    const ccomplex_t IT_0805 = IT_0139 + IT_0804;
    const ccomplex_t IT_0806 = IT_0163 + IT_0369;
    const ccomplex_t IT_0807 = (-2)*IT_0162;
    const ccomplex_t IT_0808 = IT_0530 + IT_0807;
    const ccomplex_t IT_0809 = (-0.5)*IT_0415;
    const ccomplex_t IT_0810 = s_12 + IT_0809;
    const ccomplex_t IT_0811 = (-0.5)*IT_0527;
    const ccomplex_t IT_0812 = (-2)*IT_0758 + (-2)*IT_0759 + 4*IT_0761;
    const ccomplex_t IT_0813 = IT_0757 + IT_0812;
    const ccomplex_t IT_0814 = IT_0001*IT_0813;
    const ccomplex_t IT_0815 = 0.25*IT_0814;
    const ccomplex_t IT_0816 = IT_0810 + IT_0811 + IT_0815;
    const ccomplex_t IT_0817 = 8*IT_0730;
    const ccomplex_t IT_0818 = (-4)*IT_0224;
    const ccomplex_t IT_0819 = (-2)*IT_0423;
    const ccomplex_t IT_0820 = (-2)*IT_0158;
    const ccomplex_t IT_0821 = (-2)*IT_0200;
    const ccomplex_t IT_0822 = -IT_0098;
    const ccomplex_t IT_0823 = IT_0600*IT_0684 + IT_0661*IT_0692 + IT_0591
      *IT_0700 + IT_0306*IT_0703 + IT_0611*IT_0712 + IT_0310*IT_0755 + (-8)
      *IT_0539*IT_0768 + IT_0731*(16*IT_0577 + (-16)*IT_0803) + -IT_0553*IT_0805
       + IT_0737*IT_0806 + IT_0745*IT_0808 + IT_0816*IT_0817 + IT_0775*IT_0818 +
       IT_0620*IT_0819 + IT_0785*IT_0820 + IT_0796*IT_0821 + (-8)*IT_0583*
      (IT_0094 + IT_0822);
    const ccomplex_t IT_0824 = (-2)*IT_0758 + (-2)*IT_0759 + 4*IT_0760;
    const ccomplex_t IT_0825 = IT_0757 + IT_0824;
    const ccomplex_t IT_0826 = IT_0001*IT_0825;
    const ccomplex_t IT_0827 = 0.25*IT_0826;
    const ccomplex_t IT_0828 = (-2)*IT_0002;
    const ccomplex_t IT_0829 = IT_0756 + IT_0828;
    const ccomplex_t IT_0830 = IT_0001*IT_0829;
    const ccomplex_t IT_0831 = 0.25*IT_0830;
    const ccomplex_t IT_0832 = 1 + IT_0831;
    const ccomplex_t IT_0833 = IT_0000*IT_0832;
    const ccomplex_t IT_0834 = 8*IT_0577;
    const ccomplex_t IT_0835 = IT_0258*IT_0297;
    const ccomplex_t IT_0836 = IT_0244*IT_0534*IT_0835;
    const ccomplex_t IT_0837 = IT_0240*IT_0273;
    const ccomplex_t IT_0838 = IT_0226*IT_0531*IT_0837;
    const ccomplex_t IT_0839 = m_N_2*IT_0283;
    const ccomplex_t IT_0840 = m_N_2*IT_0300;
    const ccomplex_t IT_0841 = -IT_0836 + -IT_0838 + IT_0839 + IT_0840;
    const ccomplex_t IT_0842 = 0.5*IT_0189 + -IT_0592 + (-0.5)*IT_0593;
    const ccomplex_t IT_0843 = IT_0168 + IT_0842;
    const ccomplex_t IT_0844 = IT_0145*IT_0843;
    const ccomplex_t IT_0845 = -IT_0474;
    const ccomplex_t IT_0846 = IT_0007 + IT_0844 + IT_0845;
    const ccomplex_t IT_0847 = m_N_3*IT_0846;
    const ccomplex_t IT_0848 = m_N_3*IT_0743;
    const ccomplex_t IT_0849 = (-4)*IT_0848;
    const ccomplex_t IT_0850 = 8*IT_0041;
    const ccomplex_t IT_0851 = (-8)*IT_0139;
    const ccomplex_t IT_0852 = IT_0850 + IT_0851;
    const ccomplex_t IT_0853 = 8*IT_0803;
    const ccomplex_t IT_0854 = (-8)*IT_0577;
    const ccomplex_t IT_0855 = IT_0853 + IT_0854;
    const ccomplex_t IT_0856 = IT_0168 + IT_0687;
    const ccomplex_t IT_0857 = IT_0145*IT_0856;
    const ccomplex_t IT_0858 = IT_0475 + IT_0515;
    const ccomplex_t IT_0859 = (-4)*IT_0858;
    const ccomplex_t IT_0860 = IT_0857 + IT_0859;
    const ccomplex_t IT_0861 = m_N_3*IT_0860;
    const ccomplex_t IT_0862 = IT_0542*IT_0788;
    const ccomplex_t IT_0863 = (-0.5)*IT_0787 + 0.25*IT_0862;
    const ccomplex_t IT_0864 = s_23 + IT_0863;
    const ccomplex_t IT_0865 = m_N_2*IT_0864;
    const ccomplex_t IT_0866 = 4*IT_0865;
    const ccomplex_t IT_0867 = 2*IT_0158;
    const ccomplex_t IT_0868 = s_13*s_34;
    const ccomplex_t IT_0869 = (-2)*IT_0868;
    const ccomplex_t IT_0870 = IT_0779 + IT_0869;
    const ccomplex_t IT_0871 = IT_0542*IT_0870;
    const ccomplex_t IT_0872 = (-0.5)*IT_0781 + (-0.25)*IT_0871;
    const ccomplex_t IT_0873 = s_13 + IT_0872;
    const ccomplex_t IT_0874 = m_N_3*IT_0873;
    const ccomplex_t IT_0875 = (-4)*IT_0874;
    const ccomplex_t IT_0876 = 2*IT_0200;
    const ccomplex_t IT_0877 = (-2)*IT_0493;
    const ccomplex_t IT_0878 = (-2)*IT_0412 + -IT_0524 + 2*IT_0566;
    const ccomplex_t IT_0879 = IT_0411 + IT_0878;
    const ccomplex_t IT_0880 = IT_0145*IT_0879;
    const ccomplex_t IT_0881 = -IT_0508 + (-0.5)*IT_0880;
    const ccomplex_t IT_0882 = s_12 + IT_0881;
    const ccomplex_t IT_0883 = m_N_3*IT_0882;
    const ccomplex_t IT_0884 = (-2)*IT_0883;
    const ccomplex_t IT_0885 = (IT_0094 + -IT_0098)*IT_0553;
    const ccomplex_t IT_0886 = -IT_0885;
    const ccomplex_t IT_0887 = IT_0470*IT_0630 + IT_0661*IT_0681 + IT_0573
      *IT_0684 + IT_0306*IT_0698 + IT_0565*IT_0700 + 16*IT_0539*IT_0731 +
       IT_0817*IT_0833 + 4*IT_0224*IT_0847 + IT_0806*IT_0849 + IT_0550*IT_0852 +
       IT_0768*IT_0855 + IT_0808*IT_0861 + IT_0866*IT_0867 + IT_0875*IT_0876 +
       IT_0707*IT_0877 + IT_0701*IT_0884 + IT_0886;
    const ccomplex_t IT_0888 = IT_0655 + IT_0807;
    const ccomplex_t IT_0889 = (-4)*IT_0305;
    const ccomplex_t IT_0890 = (IT_0163 + IT_0263)*IT_0600 + IT_0692*IT_0701 +
       IT_0310*IT_0703 + IT_0470*IT_0712 + (-16)*IT_0730*IT_0731 + IT_0686
      *IT_0737 + IT_0699*IT_0745 + IT_0306*IT_0755 + IT_0796*IT_0820 + IT_0785
      *IT_0821 + (-8)*IT_0539*IT_0833 + IT_0583*IT_0852 + IT_0816*IT_0855 +
       IT_0620*IT_0877 + IT_0886 + IT_0591*IT_0888 + IT_0775*IT_0889;
    const ccomplex_t IT_0891 = IT_0000*IT_0858;
    const ccomplex_t IT_0892 = s_12 + IT_0767;
    const ccomplex_t IT_0893 = IT_0000*IT_0892;
    const ccomplex_t IT_0894 = -IT_0893;
    const ccomplex_t IT_0895 = IT_0891 + IT_0894;
    const ccomplex_t IT_0896 = IT_0000*IT_0734;
    const ccomplex_t IT_0897 = -IT_0896;
    const ccomplex_t IT_0898 = IT_0893 + IT_0897;
    const ccomplex_t IT_0899 = -IT_0407;
    const ccomplex_t IT_0900 = 0.166666666666667*IT_0511;
    const ccomplex_t IT_0901 = 0.166666666666667*IT_0644;
    const ccomplex_t IT_0902 = 0.166666666666667*IT_0520;
    const ccomplex_t IT_0903 = 0.166666666666667*IT_0650;
    const ccomplex_t IT_0904 = 0.166666666666667*IT_0529;
    const ccomplex_t IT_0905 = 4*IT_0162;
    const ccomplex_t IT_0906 = 0.166666666666667*IT_0652;
    const ccomplex_t IT_0907 = 0.166666666666667*IT_0573;
    const ccomplex_t IT_0908 = 0.166666666666667*IT_0600;
    const ccomplex_t IT_0909 = 0.166666666666667*IT_0654;
    const ccomplex_t IT_0910 = 0.166666666666667*IT_0565;
    const ccomplex_t IT_0911 = 0.166666666666667*IT_0591;
    const ccomplex_t IT_0912 = IT_0683 + IT_0685;
    const ccomplex_t IT_0913 = 0.166666666666667*IT_0912;
    const ccomplex_t IT_0914 = 0.166666666666667*IT_0155;
    const ccomplex_t IT_0915 = IT_0303*IT_0895 + IT_0898*IT_0899 + IT_0514
      *IT_0900 + IT_0647*IT_0901 + IT_0523*IT_0902 + IT_0607*IT_0903 + IT_0808
      *IT_0904 + (IT_0263 + IT_0369 + IT_0905)*IT_0906 + IT_0631*IT_0907 +
       IT_0704*IT_0908 + IT_0888*IT_0909 + IT_0682*IT_0910 + IT_0621*IT_0911 +
       IT_0306*IT_0913 + IT_0611*IT_0914;
    const ccomplex_t IT_0916 = (-2)*conj(IT_0368);
    const ccomplex_t IT_0917 = conj(IT_0303)*IT_0895 + -conj(IT_0407)*IT_0898 
      + IT_0634*IT_0900 + IT_0716*IT_0901 + IT_0719*IT_0902 + IT_0637*IT_0903 +
       IT_0574*IT_0904 + (IT_0265 + IT_0371)*IT_0906 + IT_0639*IT_0907 + IT_0723
      *IT_0908 + IT_0720*IT_0910 + IT_0638*IT_0911 + IT_0308*IT_0913 + IT_0748
      *IT_0914 + IT_0909*IT_0916;
    const ccomplex_t IT_0918 = -IT_0368;
    const ccomplex_t IT_0919 = 0.166666666666667*IT_0745;
    const ccomplex_t IT_0920 = 0.166666666666667*IT_0861;
    const ccomplex_t IT_0921 = 0.166666666666667*IT_0466;
    const ccomplex_t IT_0922 = -conj(IT_0368);
    const ccomplex_t IT_0923 = (-2)*IT_0758;
    const ccomplex_t IT_0924 = IT_0757 + IT_0923;
    const ccomplex_t IT_0925 = IT_0145*IT_0924;
    const ccomplex_t IT_0926 = 0.5*IT_0925;
    const ccomplex_t IT_0927 = s_14*s_23*IT_0225;
    const ccomplex_t IT_0928 = s_13*s_24*IT_0225;
    const ccomplex_t IT_0929 = IT_0927 + IT_0928;
    const ccomplex_t IT_0930 = s_12*s_34*IT_0225;
    const ccomplex_t IT_0931 = -IT_0930;
    const ccomplex_t IT_0932 = IT_0929 + IT_0931;
    const ccomplex_t IT_0933 = IT_0542*IT_0932;
    const ccomplex_t IT_0934 = 0.5*IT_0933;
    const ccomplex_t IT_0935 = (-2)*IT_0758 + 2*IT_0760 + 2*IT_0761 + (-2)
      *IT_0762;
    const ccomplex_t IT_0936 = IT_0757 + IT_0935;
    const ccomplex_t IT_0937 = IT_0145*IT_0936;
    const ccomplex_t IT_0938 = (-0.5)*IT_0937;
    const ccomplex_t IT_0939 = IT_0411 + IT_0413 + IT_0926 + IT_0934 + IT_0938;
    const ccomplex_t IT_0940 = pow(m_W, 6);
    const ccomplex_t IT_0941 = s_12*IT_0940;
    const ccomplex_t IT_0942 = s_14*s_24*IT_0756;
    const ccomplex_t IT_0943 = (-2)*IT_0942;
    const ccomplex_t IT_0944 = IT_0941 + IT_0943;
    const ccomplex_t IT_0945 = IT_0001*IT_0944;
    const ccomplex_t IT_0946 = 0.25*IT_0945;
    const ccomplex_t IT_0947 = IT_0411 + IT_0946;
    const ccomplex_t IT_0948 = -IT_0175;
    const ccomplex_t IT_0949 = IT_0181 + IT_0948;
    const ccomplex_t IT_0950 = s_12*s_23*IT_0225;
    const ccomplex_t IT_0951 = s_14*s_23*s_24;
    const ccomplex_t IT_0952 = s_13*IT_0008*IT_0225;
    const ccomplex_t IT_0953 = s_14*s_34*IT_0008;
    const ccomplex_t IT_0954 = (-2)*IT_0951 + -IT_0952 + 2*IT_0953;
    const ccomplex_t IT_0955 = IT_0950 + IT_0954;
    const ccomplex_t IT_0956 = IT_0145*IT_0955;
    const ccomplex_t IT_0957 = (-0.5)*IT_0956;
    const ccomplex_t IT_0958 = IT_0949 + IT_0957;
    const ccomplex_t IT_0959 = 0.5*IT_0945;
    const ccomplex_t IT_0960 = (-4)*IT_0412;
    const ccomplex_t IT_0961 = IT_0411 + IT_0959 + IT_0960;
    const ccomplex_t IT_0962 = IT_0145*IT_0813;
    const ccomplex_t IT_0963 = -IT_0962;
    const ccomplex_t IT_0964 = IT_0961 + IT_0963;
    const ccomplex_t IT_0965 = -IT_0170;
    const ccomplex_t IT_0966 = IT_0177 + IT_0965;
    const ccomplex_t IT_0967 = s_12*s_13*IT_0225;
    const ccomplex_t IT_0968 = s_23*IT_0189;
    const ccomplex_t IT_0969 = s_12*s_14*s_34;
    const ccomplex_t IT_0970 = s_23*IT_0007*IT_0225;
    const ccomplex_t IT_0971 = 2*IT_0968 + (-2)*IT_0969 + -IT_0970;
    const ccomplex_t IT_0972 = IT_0967 + IT_0971;
    const ccomplex_t IT_0973 = IT_0145*IT_0972;
    const ccomplex_t IT_0974 = (-0.5)*IT_0973;
    const ccomplex_t IT_0975 = IT_0966 + IT_0974;
    const ccomplex_t IT_0976 = -conj(IT_0303);
    const ccomplex_t IT_0977 = (-2)*IT_0444;
    const ccomplex_t IT_0978 = IT_0443 + IT_0977;
    const ccomplex_t IT_0979 = IT_0145*IT_0978;
    const ccomplex_t IT_0980 = 0.5*IT_0787 + (-0.5)*IT_0862 + (-0.5)*IT_0979;
    const ccomplex_t IT_0981 = s_23 + IT_0980;
    const ccomplex_t IT_0982 = m_N_2*IT_0981;
    const ccomplex_t IT_0983 = 2*IT_0982;
    const ccomplex_t IT_0984 = 0.25*IT_0983;
    const ccomplex_t IT_0985 = (-2)*IT_0435;
    const ccomplex_t IT_0986 = IT_0434 + IT_0985;
    const ccomplex_t IT_0987 = IT_0145*IT_0986;
    const ccomplex_t IT_0988 = (-0.5)*IT_0780 + 0.5*IT_0781 + (-0.5)*IT_0987;
    const ccomplex_t IT_0989 = s_13 + IT_0988;
    const ccomplex_t IT_0990 = m_N_3*IT_0989;
    const ccomplex_t IT_0991 = 2*IT_0990;
    const ccomplex_t IT_0992 = 0.25*IT_0991;
    const ccomplex_t IT_0993 = (-2)*conj(IT_0841);
    const ccomplex_t IT_0994 = s_23*IT_0756;
    const ccomplex_t IT_0995 = s_24*s_34*IT_0225;
    const ccomplex_t IT_0996 = (-2)*IT_0995;
    const ccomplex_t IT_0997 = IT_0994 + IT_0996;
    const ccomplex_t IT_0998 = IT_0001*IT_0997;
    const ccomplex_t IT_0999 = (-0.666666666666667)*IT_0418 + (
      -0.166666666666667)*IT_0998;
    const ccomplex_t IT_1000 = s_23 + IT_0999;
    const ccomplex_t IT_1001 = m_N_2*IT_1000;
    const ccomplex_t IT_1002 = (-6)*IT_1001;
    const ccomplex_t IT_1003 = 0.25*IT_1002;
    const ccomplex_t IT_1004 = m_N_3*IT_0145;
    const ccomplex_t IT_1005 = s_13*IT_0756;
    const ccomplex_t IT_1006 = IT_0145*IT_1005;
    const ccomplex_t IT_1007 = (-2)*IT_0435 + (-0.5)*IT_1006;
    const ccomplex_t IT_1008 = IT_0434 + IT_1007;
    const ccomplex_t IT_1009 = IT_1004*IT_1008;
    const ccomplex_t IT_1010 = (-2)*IT_1009;
    const ccomplex_t IT_1011 = 0.25*IT_1010;
    const ccomplex_t IT_1012 = 2*conj(IT_0730);
    const ccomplex_t IT_1013 = s_14*s_34*IT_0225;
    const ccomplex_t IT_1014 = (-2)*IT_1013;
    const ccomplex_t IT_1015 = IT_1005 + IT_1014;
    const ccomplex_t IT_1016 = IT_0001*IT_1015;
    const ccomplex_t IT_1017 = (-0.666666666666667)*IT_0987 +
       0.166666666666667*IT_1016;
    const ccomplex_t IT_1018 = s_13 + IT_1017;
    const ccomplex_t IT_1019 = m_N_3*IT_1018;
    const ccomplex_t IT_1020 = 6*IT_1019;
    const ccomplex_t IT_1021 = 0.25*IT_1020;
    const ccomplex_t IT_1022 = 2*conj(IT_0803);
    const ccomplex_t IT_1023 = IT_0723 + IT_1022;
    const ccomplex_t IT_1024 = IT_0001*IT_0994;
    const ccomplex_t IT_1025 = (-10)*s_23;
    const ccomplex_t IT_1026 = IT_1024 + IT_1025;
    const ccomplex_t IT_1027 = m_N_2*IT_1026;
    const ccomplex_t IT_1028 = 0.25*IT_1027;
    const ccomplex_t IT_1029 = -IT_0002;
    const ccomplex_t IT_1030 = IT_0756 + IT_1029;
    const ccomplex_t IT_1031 = IT_0145*IT_1030;
    const ccomplex_t IT_1032 = -IT_0225 + 0.25*IT_1031;
    const ccomplex_t IT_1033 = IT_0608 + IT_1032;
    const ccomplex_t IT_1034 = IT_0000*IT_1033;
    const ccomplex_t IT_1035 = 4*IT_1034;
    const ccomplex_t IT_1036 = 0.25*IT_1035;
    const ccomplex_t IT_1037 = 2*conj(IT_0200);
    const ccomplex_t IT_1038 = s_34*IT_0225;
    const ccomplex_t IT_1039 = IT_0542*IT_1038;
    const ccomplex_t IT_1040 = IT_0145*IT_0829;
    const ccomplex_t IT_1041 = (-0.333333333333333)*IT_1039 + (
      -0.333333333333333)*IT_1040;
    const ccomplex_t IT_1042 = IT_0225 + IT_1041;
    const ccomplex_t IT_1043 = IT_0000*IT_1042;
    const ccomplex_t IT_1044 = 3*IT_1043;
    const ccomplex_t IT_1045 = 0.25*IT_1044;
    const ccomplex_t IT_1046 = (-2)*conj(IT_0454);
    const ccomplex_t IT_1047 = m_N_2*m_N_3*IT_0145;
    const ccomplex_t IT_1048 = (-0.25)*IT_0756;
    const ccomplex_t IT_1049 = IT_0002 + IT_1048;
    const ccomplex_t IT_1050 = IT_1047*IT_1049;
    const ccomplex_t IT_1051 = (-4)*IT_1050;
    const ccomplex_t IT_1052 = 0.25*IT_1051;
    const ccomplex_t IT_1053 = (-2)*conj(IT_0452);
    const ccomplex_t IT_1054 = IT_0000*IT_0225;
    const ccomplex_t IT_1055 = 9*IT_1054;
    const ccomplex_t IT_1056 = 0.25*IT_1055;
    const ccomplex_t IT_1057 = (-4)*conj(IT_0493);
    const ccomplex_t IT_1058 = s_12*IT_0002*IT_0145;
    const ccomplex_t IT_1059 = -IT_1058;
    const ccomplex_t IT_1060 = IT_0757 + IT_0760 + IT_0761;
    const ccomplex_t IT_1061 = (-2)*IT_0758 + -IT_0762;
    const ccomplex_t IT_1062 = IT_1060 + IT_1061;
    const ccomplex_t IT_1063 = IT_0145*IT_1062;
    const ccomplex_t IT_1064 = (-0.5)*IT_1063;
    const ccomplex_t IT_1065 = -IT_0584;
    const ccomplex_t IT_1066 = IT_0411 + IT_0566 + IT_1059 + IT_1064 + IT_1065;
    const ccomplex_t IT_1067 = 0.25*IT_1066;
    const ccomplex_t IT_1068 = (-4)*conj(IT_0305);
    const ccomplex_t IT_1069 = s_13*s_14*s_24;
    const ccomplex_t IT_1070 = IT_0968 + IT_1069;
    const ccomplex_t IT_1071 = -IT_0969 + -IT_0970;
    const ccomplex_t IT_1072 = IT_1070 + IT_1071;
    const ccomplex_t IT_1073 = IT_0145*IT_1072;
    const ccomplex_t IT_1074 = (-0.5)*IT_1073;
    const ccomplex_t IT_1075 = -IT_0171;
    const ccomplex_t IT_1076 = IT_0170 + IT_0187 + IT_0191 + IT_1074 + IT_1075;
    const ccomplex_t IT_1077 = 0.25*IT_1076;
    const ccomplex_t IT_1078 = (-8)*conj(IT_0305);
    const ccomplex_t IT_1079 = s_13*IT_0194;
    const ccomplex_t IT_1080 = IT_0951 + IT_1079;
    const ccomplex_t IT_1081 = s_12*s_24*s_34;
    const ccomplex_t IT_1082 = -IT_0952 + -IT_1081;
    const ccomplex_t IT_1083 = IT_1080 + IT_1082;
    const ccomplex_t IT_1084 = IT_0145*IT_1083;
    const ccomplex_t IT_1085 = (-0.25)*IT_1084;
    const ccomplex_t IT_1086 = (-0.5)*IT_0176;
    const ccomplex_t IT_1087 = 0.5*IT_0195;
    const ccomplex_t IT_1088 = (-0.5)*IT_0197;
    const ccomplex_t IT_1089 = (-0.5)*IT_0175;
    const ccomplex_t IT_1090 = IT_0181 + IT_1085 + IT_1086 + IT_1087 + IT_1088
       + IT_1089;
    const ccomplex_t IT_1091 = 0.25*IT_1090;
    const ccomplex_t IT_1092 = conj(IT_0158)*IT_0939 + (-2)*conj(IT_0458)
      *IT_0947 + conj(IT_0407)*IT_0958 + -conj(IT_0460)*IT_0964 + IT_0975
      *IT_0976 + IT_0637*IT_0984 + IT_0634*IT_0992 + IT_0993*IT_1003 + IT_0720
      *IT_1011 + IT_1012*IT_1021 + IT_1023*IT_1028 + IT_0472*IT_1036 + IT_1037
      *IT_1045 + IT_1046*IT_1052 + IT_1053*IT_1056 + IT_1057*IT_1067 + IT_1068
      *IT_1077 + IT_1078*IT_1091;
    const ccomplex_t IT_1093 = s_23*m_N_2;
    const ccomplex_t IT_1094 = IT_0005*IT_1093;
    const ccomplex_t IT_1095 = s_13*m_N_3;
    const ccomplex_t IT_1096 = IT_0005*IT_1095;
    const ccomplex_t IT_1097 = m_N_2*m_N_3*IT_0225;
    const ccomplex_t IT_1098 = IT_0005*IT_1097;
    const ccomplex_t IT_1099 = (-4)*IT_0158;
    const ccomplex_t IT_1100 = (1 + 0.5*IT_0003)*(IT_0411 + (-2)*IT_0412);
    const ccomplex_t IT_1101 = -IT_0200;
    const ccomplex_t IT_1102 = 0.25*IT_0150;
    const ccomplex_t IT_1103 = 0.25*IT_0154;
    const ccomplex_t IT_1104 = 0.25*IT_0866;
    const ccomplex_t IT_1105 = 0.25*IT_0785;
    const ccomplex_t IT_1106 = s_34*IT_0145;
    const ccomplex_t IT_1107 = IT_0543 + IT_0544;
    const ccomplex_t IT_1108 = s_13*s_23*s_34*IT_0145;
    const ccomplex_t IT_1109 = (-2)*IT_1108;
    const ccomplex_t IT_1110 = IT_1107 + IT_1109;
    const ccomplex_t IT_1111 = IT_1106*IT_1110;
    const ccomplex_t IT_1112 = 0.25*IT_1111;
    const ccomplex_t IT_1113 = IT_0167 + IT_0169 + IT_0170 + IT_0171 + IT_0180
       + IT_0186 + IT_0188 + IT_0191;
    const ccomplex_t IT_1114 = 0.25*IT_1113;
    const ccomplex_t IT_1115 = IT_0173 + IT_0174 + IT_0175 + IT_0176 + IT_0184
       + IT_0193 + IT_0196 + IT_0198;
    const ccomplex_t IT_1116 = 0.25*IT_1115;
    const ccomplex_t IT_1117 = IT_0142*IT_1094 + IT_0805*IT_1096 + 0.25
      *IT_1098*IT_1099 + IT_1100*IT_1101 + IT_0684*IT_1102 + IT_0686*IT_1103 +
       IT_0631*IT_1104 + IT_0621*IT_1105 + IT_0819*IT_1112 + IT_0263*IT_1114 +
       IT_0369*IT_1116;
    const ccomplex_t IT_1118 = 0.25*IT_0875;
    const ccomplex_t IT_1119 = 0.25*IT_0796;
    const ccomplex_t IT_1120 = (IT_0041 + 0.25*IT_0505)*IT_1094 + IT_0100
      *IT_1096 + 0.25*IT_1099*IT_1100 + IT_1098*IT_1101 + IT_0263*IT_1102 +
       IT_0369*IT_1103 + IT_0877*IT_1112 + IT_0684*IT_1114 + IT_0686*IT_1116 +
       IT_0631*IT_1118 + IT_0621*IT_1119;
    const ccomplex_t IT_1121 = (-2)*conj(IT_0493);
    const ccomplex_t IT_1122 = (-0.25)*IT_0520;
    const ccomplex_t IT_1123 = (-0.25)*IT_0674;
    const ccomplex_t IT_1124 = (-0.25)*IT_0673;
    const ccomplex_t IT_1125 = (-0.25)*IT_0683;
    const ccomplex_t IT_1126 = (-0.25)*IT_0663;
    const ccomplex_t IT_1127 = (-0.25)*IT_0664;
    const ccomplex_t IT_1128 = (-0.25)*IT_0685;
    const ccomplex_t IT_1129 = (-0.25)*IT_0884;
    const ccomplex_t IT_1130 = (-0.25)*IT_0755;
    const ccomplex_t IT_1131 = (-0.25)*IT_0607*IT_0644 + -IT_0305*IT_0669 + 3
      *IT_0224*IT_0672 + -IT_0577*IT_0775 + IT_0539*IT_0847 + IT_0514*IT_1122 +
       IT_0647*IT_1123 + IT_0523*IT_1124 + IT_0263*IT_1125 + IT_0530*IT_1126 +
       IT_0655*IT_1127 + IT_0369*IT_1128 + IT_0682*IT_1129 + IT_0704*IT_1130;
    const ccomplex_t IT_1132 = (-0.25)*IT_0637*IT_0644 + -conj(IT_0577)
      *IT_0775 + conj(IT_0539)*IT_0847 + IT_0634*IT_1122 + IT_0716*IT_1123 +
       IT_0719*IT_1124 + IT_0265*IT_1125 + IT_0574*IT_1126 + IT_0916*IT_1127 +
       IT_0371*IT_1128 + IT_0720*IT_1129 + IT_0723*IT_1130;
    const ccomplex_t IT_1133 = IT_0145*IT_0825;
    const ccomplex_t IT_1134 = -IT_1133;
    const ccomplex_t IT_1135 = IT_0961 + IT_1134;
    const ccomplex_t IT_1136 = -IT_0952 + 2*IT_1079 + (-2)*IT_1081;
    const ccomplex_t IT_1137 = IT_0950 + IT_1136;
    const ccomplex_t IT_1138 = IT_0145*IT_1137;
    const ccomplex_t IT_1139 = (-0.5)*IT_1138;
    const ccomplex_t IT_1140 = IT_0949 + IT_1139;
    const ccomplex_t IT_1141 = s_24*s_34*IT_0007;
    const ccomplex_t IT_1142 = -IT_0970 + (-2)*IT_1069 + 2*IT_1141;
    const ccomplex_t IT_1143 = IT_0967 + IT_1142;
    const ccomplex_t IT_1144 = IT_0145*IT_1143;
    const ccomplex_t IT_1145 = (-0.5)*IT_1144;
    const ccomplex_t IT_1146 = IT_0966 + IT_1145;
    const ccomplex_t IT_1147 = (-0.25)*IT_0983;
    const ccomplex_t IT_1148 = (-0.25)*IT_0991;
    const ccomplex_t IT_1149 = (-0.666666666666667)*IT_0979 +
       0.166666666666667*IT_0998;
    const ccomplex_t IT_1150 = s_23 + IT_1149;
    const ccomplex_t IT_1151 = m_N_2*IT_1150;
    const ccomplex_t IT_1152 = 6*IT_1151;
    const ccomplex_t IT_1153 = (-0.25)*IT_1152;
    const ccomplex_t IT_1154 = IT_0001*IT_1005;
    const ccomplex_t IT_1155 = (-10)*s_13;
    const ccomplex_t IT_1156 = IT_1154 + IT_1155;
    const ccomplex_t IT_1157 = m_N_3*IT_1156;
    const ccomplex_t IT_1158 = (-0.25)*IT_1157;
    const ccomplex_t IT_1159 = (-0.666666666666667)*IT_0410 + (
      -0.166666666666667)*IT_1016;
    const ccomplex_t IT_1160 = s_13 + IT_1159;
    const ccomplex_t IT_1161 = m_N_3*IT_1160;
    const ccomplex_t IT_1162 = (-6)*IT_1161;
    const ccomplex_t IT_1163 = (-0.25)*IT_1162;
    const ccomplex_t IT_1164 = m_N_2*IT_0145;
    const ccomplex_t IT_1165 = IT_0145*IT_0994;
    const ccomplex_t IT_1166 = (-2)*IT_0444 + (-0.5)*IT_1165;
    const ccomplex_t IT_1167 = IT_0443 + IT_1166;
    const ccomplex_t IT_1168 = IT_1164*IT_1167;
    const ccomplex_t IT_1169 = (-2)*IT_1168;
    const ccomplex_t IT_1170 = (-0.25)*IT_1169;
    const ccomplex_t IT_1171 = (-0.25)*IT_1035;
    const ccomplex_t IT_1172 = (-0.25)*IT_1044;
    const ccomplex_t IT_1173 = -IT_0566;
    const ccomplex_t IT_1174 = IT_0411 + IT_0584 + IT_1059 + IT_1064 + IT_1173;
    const ccomplex_t IT_1175 = (-0.25)*IT_1174;
    const ccomplex_t IT_1176 = (-0.25)*IT_1073;
    const ccomplex_t IT_1177 = (-0.5)*IT_0171;
    const ccomplex_t IT_1178 = (-0.5)*IT_0187;
    const ccomplex_t IT_1179 = 0.5*IT_0190;
    const ccomplex_t IT_1180 = (-0.5)*IT_0170;
    const ccomplex_t IT_1181 = IT_0177 + IT_1176 + IT_1177 + IT_1178 + IT_1179
       + IT_1180;
    const ccomplex_t IT_1182 = (-0.25)*IT_1181;
    const ccomplex_t IT_1183 = (-0.5)*IT_1084;
    const ccomplex_t IT_1184 = -IT_0176;
    const ccomplex_t IT_1185 = IT_0175 + IT_0196 + IT_0197 + IT_1183 + IT_1184;
    const ccomplex_t IT_1186 = (-0.25)*IT_1185;
    const ccomplex_t IT_1187 = -conj(IT_0158)*IT_0939 + 2*conj(IT_0460)
      *IT_0947 + (-0.25)*IT_1051*IT_1053 + (-0.25)*IT_1046*IT_1055 + conj
      (IT_0458)*IT_1135 + conj(IT_0407)*IT_1140 + IT_0976*IT_1146 + IT_0637
      *IT_1147 + IT_0634*IT_1148 + IT_0993*IT_1153 + IT_0720*IT_1158 + IT_1012
      *IT_1163 + IT_1023*IT_1170 + IT_0472*IT_1171 + IT_1037*IT_1172 + IT_1057
      *IT_1175 + IT_1078*IT_1182 + IT_1068*IT_1186;
    const ccomplex_t IT_1188 = 2*conj(IT_0158);
    const ccomplex_t IT_1189 = (-4)*conj(IT_0423);
    const ccomplex_t IT_1190 = (-8)*conj(IT_0224);
    const ccomplex_t IT_1191 = 4*conj(IT_0162);
    const ccomplex_t IT_1192 = 0.25*IT_1191;
    const ccomplex_t IT_1193 = -IT_0975;
    const ccomplex_t IT_1194 = (-4)*conj(IT_0162);
    const ccomplex_t IT_1195 = (-0.25)*IT_1194;
    const ccomplex_t IT_1196 = conj(IT_0200)*IT_0939 + (-2)*conj(IT_0454)
      *IT_0947 + -conj(IT_0452)*IT_0964 + IT_0716*IT_0984 + IT_0719*IT_0992 +
       IT_0720*IT_1003 + IT_0993*IT_1011 + IT_1021*IT_1023 + IT_1012*IT_1028 +
       IT_0748*IT_1036 + IT_0776*IT_1077 + IT_1045*IT_1188 + IT_1067*IT_1189 +
       IT_1091*IT_1190 + IT_0958*(conj(IT_0368) + IT_1192) + IT_1193*(conj
      (IT_0262) + IT_1195);
    const ccomplex_t IT_1197 = -conj(IT_0200)*IT_0939 + 2*conj(IT_0452)
      *IT_0947 + conj(IT_0454)*IT_1135 + IT_0716*IT_1147 + IT_0719*IT_1148 +
       IT_0720*IT_1153 + IT_0993*IT_1158 + IT_1023*IT_1163 + IT_1012*IT_1170 +
       IT_0748*IT_1171 + IT_0776*IT_1186 + IT_1172*IT_1188 + IT_1175*IT_1189 +
       IT_1182*IT_1190 + -IT_1146*(conj(IT_0262) + IT_1192) + IT_1140*(conj
      (IT_0368) + IT_1195);
    const ccomplex_t IT_1198 = (-0.166666666666667)*IT_0650;
    const ccomplex_t IT_1199 = (-0.166666666666667)*IT_0652;
    const ccomplex_t IT_1200 = (-0.166666666666667)*IT_0737;
    const ccomplex_t IT_1201 = (-0.166666666666667)*IT_0849;
    const ccomplex_t IT_1202 = (-0.166666666666667)*IT_0469;
    const ccomplex_t IT_1203 = (-6)*conj(IT_0407)*((-0.166666666666667)
      *IT_0607*IT_0644 + (-0.166666666666667)*IT_0654*IT_0700 + IT_0262*IT_0893 
      + IT_0896*IT_0918 + IT_0647*IT_1198 + IT_0684*IT_1199 + IT_0704*IT_1200 +
       IT_0631*IT_1201 + IT_0611*IT_1202);
    const ccomplex_t IT_1204 = 2*IT_0803;
    const ccomplex_t IT_1205 = IT_0704 + IT_1204;
    const ccomplex_t IT_1206 = (-2)*IT_0841;
    const ccomplex_t IT_1207 = IT_0607*IT_0984;
    const ccomplex_t IT_1208 = 0.25*IT_1157;
    const ccomplex_t IT_1209 = 2*IT_0730;
    const ccomplex_t IT_1210 = 0.25*IT_1162;
    const ccomplex_t IT_1211 = IT_0470*IT_1036;
    const ccomplex_t IT_1212 = IT_0876*IT_1045;
    const ccomplex_t IT_1213 = (-2)*IT_0452;
    const ccomplex_t IT_1214 = (-2)*IT_0454;
    const ccomplex_t IT_1215 = (-4)*IT_0493;
    const ccomplex_t IT_1216 = 0.25*IT_1174;
    const ccomplex_t IT_1217 = (-8)*IT_0305;
    const ccomplex_t IT_1218 = 0.25*IT_1181;
    const ccomplex_t IT_1219 = 0.25*IT_1185;
    const ccomplex_t IT_1220 = IT_0158*IT_0939;
    const ccomplex_t IT_1221 = IT_0514*IT_0992 + IT_0899*IT_1140 + IT_0303
      *IT_1146 + 0.25*IT_1169*IT_1205 + 0.25*IT_1152*IT_1206 + IT_1207 + IT_0682
      *IT_1208 + IT_1209*IT_1210 + IT_1211 + IT_1212 + IT_1052*IT_1213 + IT_1056
      *IT_1214 + IT_1215*IT_1216 + IT_1217*IT_1218 + IT_0889*IT_1219 + IT_1220;
    const ccomplex_t IT_1222 = -IT_0407*IT_0958 + IT_0303*IT_0975 + -IT_0514
      *IT_0992 + -IT_0682*IT_1011 + -IT_0889*IT_1077 + -IT_1028*IT_1205 + 
      -IT_1003*IT_1206 + -IT_1207 + -IT_1021*IT_1209 + -IT_1211 + -IT_1212 + 
      -IT_1056*IT_1213 + -IT_1052*IT_1214 + -IT_1067*IT_1215 + -IT_1091*IT_1217 
      + -IT_1220;
    const ccomplex_t IT_1223 = (-4)*IT_0423;
    const ccomplex_t IT_1224 = (-8)*IT_0224;
    const ccomplex_t IT_1225 = (-4)*IT_0162;
    const ccomplex_t IT_1226 = (-0.25)*IT_1225;
    const ccomplex_t IT_1227 = IT_0200*IT_0939;
    const ccomplex_t IT_1228 = (IT_0368 + 0.25*IT_0905)*IT_0958 + IT_0647
      *IT_0984 + IT_0523*IT_0992 + IT_0682*IT_1003 + IT_0611*IT_1036 + IT_0867
      *IT_1045 + IT_0818*IT_1077 + IT_1021*IT_1205 + IT_1011*IT_1206 + IT_1028
      *IT_1209 + IT_1067*IT_1223 + IT_1091*IT_1224 + IT_1193*(IT_0262 + IT_1226)
       + IT_1227;
    const ccomplex_t IT_1229 = 4*IT_0647*IT_0984 + 4*IT_0523*IT_0992 + 4
      *IT_0611*IT_1036 + 4*IT_0867*IT_1045 + 4*(IT_0262 + 0.25*IT_0905)*IT_1146 
      + IT_0682*IT_1152 + 4*IT_1206*IT_1208 + IT_1169*IT_1209 + 4*IT_1205
      *IT_1210 + 4*IT_0818*IT_1219 + 4*IT_1216*IT_1223 + 4*IT_1218*IT_1224 + (-4
      )*IT_1140*(IT_0368 + IT_1226) + 4*IT_1227;
    const ccomplex_t IT_1230 = IT_0644*IT_0716;
    const ccomplex_t IT_1231 = 8*conj(IT_0041);
    const ccomplex_t IT_1232 = (-8)*conj(IT_0139);
    const ccomplex_t IT_1233 = IT_1231 + IT_1232;
    const ccomplex_t IT_1234 = (-8)*conj(IT_0577);
    const ccomplex_t IT_1235 = (-2)*conj(IT_0224);
    const ccomplex_t IT_1236 = -conj(IT_0094);
    const ccomplex_t IT_1237 = conj(IT_0098) + IT_1236;
    const ccomplex_t IT_1238 = IT_0553*IT_1237;
    const ccomplex_t IT_1239 = IT_0472*IT_0630 + IT_0666*IT_0681 + IT_0308
      *IT_0698 + IT_0573*IT_0721 + 16*conj(IT_0539)*IT_0731 + IT_0565*IT_0777 +
       8*conj(IT_0730)*IT_0833 + 8*(IT_0810 + IT_0811 + IT_0827)*conj(IT_0841) +
       4*conj(IT_0224)*IT_0847 + IT_0738*IT_0849 + IT_0747*IT_0861 + IT_0875
      *IT_1037 + IT_0707*IT_1121 + IT_0866*IT_1188 + IT_0550*IT_1233 + 8*IT_0768
      *(conj(IT_0803) + 0.125*IT_1234) + IT_0884*IT_1235 + IT_1238;
    const ccomplex_t IT_1240 = (-2)*conj(IT_0303);
    const ccomplex_t IT_1241 = (IT_0165 + IT_0265)*IT_0600 + IT_0312*IT_0703 +
       IT_0472*IT_0712 + IT_0722*IT_0737 + IT_0308*IT_0755 + IT_0786*IT_0796 +
       IT_0785*IT_0797 + IT_0769*IT_0833 + IT_0591*(IT_0746 + IT_0916) + IT_0775
      *IT_1068 + IT_0620*IT_1121 + IT_0583*IT_1233 + IT_0816*IT_1234 + IT_0692
      *IT_1235 + IT_1238 + IT_0745*IT_1240;
    const ccomplex_t IT_1242 = IT_0520*IT_0719;
    const ccomplex_t IT_1243 = IT_0634*IT_0673 + IT_0637*IT_0674 + IT_0638
      *IT_0692 + IT_0639*IT_0698 + IT_0681*IT_0720 + IT_0683*IT_0721 + IT_0685
      *IT_0722 + IT_0703*IT_0723 + IT_0664*IT_0777 + IT_1230 + IT_0663*IT_1240 +
       IT_1242;
    const ccomplex_t IT_1244 = 4*IT_0006*(conj(IT_0041)*IT_0100 + IT_0041
      *IT_0102 + conj(IT_0139)*IT_0142 + IT_0139*IT_0144) + IT_0155*(IT_0164 +
       IT_0166) + (IT_0167 + IT_0169 + IT_0170 + IT_0171 + IT_0173 + IT_0174 +
       IT_0175 + IT_0176 + IT_0180 + IT_0184 + IT_0186 + IT_0188 + IT_0191 +
       IT_0193 + IT_0196 + IT_0198)*(IT_0201 + IT_0202) + (IT_0212 + IT_0222)*
      (IT_0264 + IT_0266 + IT_0307 + IT_0309 + IT_0311 + IT_0313) + (IT_0321 +
       IT_0331)*(IT_0311 + IT_0313 + IT_0370 + IT_0372 + IT_0408 + IT_0409) +
       IT_0421*(IT_0424 + IT_0425) + (IT_0427 + IT_0429)*(IT_0430 + IT_0431 +
       IT_0432 + IT_0433) + (IT_0442 + IT_0451)*(IT_0453 + IT_0455 + IT_0456 +
       IT_0457 + IT_0459 + IT_0461 + IT_0462 + IT_0463) + (IT_0466 + IT_0469)*
      (IT_0471 + IT_0473) + (IT_0477 + IT_0481)*(IT_0483 + IT_0485) + (IT_0488 +
       IT_0491)*(IT_0495 + IT_0497 + IT_0499 + IT_0501) + 4*IT_0502*(IT_0094
      *conj(IT_0094) + IT_0139*conj(IT_0139) + 0.25*conj(IT_0098)*IT_0140 +
       IT_0098*IT_0144 + IT_0041*IT_0504 + 0.25*conj(IT_0041)*IT_0505) + conj
      (IT_0262)*(IT_0511*IT_0514 + IT_0520*IT_0523 + IT_0529*IT_0530) + 8
      *IT_0539*(IT_0550*IT_0552 + 0.125*IT_0553*IT_0555 + 0.125*IT_0371*IT_0565 
      + 0.125*IT_0573*IT_0574) + IT_0577*(IT_0552*IT_0553 + 8*IT_0555*IT_0583 +
       IT_0371*IT_0591 + IT_0574*IT_0600) + conj(IT_0493)*(IT_0150*IT_0263 +
       IT_0154*IT_0369 + IT_0306*IT_0421 + IT_0514*IT_0602 + IT_0604*IT_0607 +
       IT_0610*IT_0611 + IT_0620*IT_0621 + IT_0630*IT_0631) + IT_0493*(IT_0150
      *IT_0265 + IT_0154*IT_0371 + IT_0308*IT_0421 + IT_0602*IT_0634 + IT_0604
      *IT_0637 + IT_0620*IT_0638 + IT_0630*IT_0639) + conj(IT_0368)*(IT_0644
      *IT_0647 + IT_0607*IT_0650 + IT_0263*IT_0652 + IT_0654*IT_0655) + IT_0660
      *IT_0662 + (IT_0663 + IT_0664)*(IT_0665 + IT_0667) + 4*conj(IT_0224)
      *IT_0705 + IT_0713 + IT_0423*(IT_0265*IT_0466 + IT_0371*IT_0469 + IT_0496
      *IT_0710 + IT_0604*IT_0716 + IT_0602*IT_0719 + IT_0707*IT_0720 + IT_0150
      *IT_0721 + IT_0154*IT_0722 + IT_0712*IT_0723) + IT_0730*IT_0798 + conj
      (IT_0730)*IT_0823 + 8*conj(IT_0577)*(0.125*IT_0369*IT_0591 + 0.125*IT_0530
      *IT_0600 + IT_0583*IT_0805 + 0.125*IT_0553*(IT_0094 + IT_0822) + IT_0539
      *IT_0833 + 0.125*IT_0816*IT_0834) + conj(IT_0539)*(IT_0369*IT_0565 +
       IT_0530*IT_0573 + IT_0553*IT_0805 + 8*IT_0550*(IT_0094 + IT_0822) + 8
      *IT_0539*(IT_0810 + IT_0811 + IT_0827) + IT_0833*IT_0834) + conj(IT_0841)
      *IT_0887 + conj(IT_0803)*IT_0890 + 6*conj(IT_0162)*IT_0915 + 6*IT_0162
      *IT_0917 + 6*conj(IT_0303)*(IT_0262*IT_0891 + IT_0523*IT_0900 + IT_0514
      *IT_0902 + IT_0699*IT_0904 + IT_0893*IT_0918 + IT_0621*IT_0919 + IT_0682
      *IT_0920 + IT_0611*IT_0921) + 6*IT_0303*(conj(IT_0262)*IT_0891 + IT_0719
      *IT_0900 + IT_0634*IT_0902 + IT_0638*IT_0919 + IT_0720*IT_0920 + IT_0748
      *IT_0921 + IT_0893*IT_0922) + 4*IT_0460*IT_1092 + 4*IT_0200*(IT_0144
      *IT_1094 + IT_0555*IT_1096 + IT_0721*IT_1102 + IT_0722*IT_1103 + IT_0639
      *IT_1104 + IT_0638*IT_1105 + IT_0778*IT_1112 + IT_0265*IT_1114 + IT_0371
      *IT_1116) + 4*conj(IT_0200)*IT_1117 + 4*conj(IT_0158)*IT_1120 + 4*IT_0158*
      (IT_0504*IT_1094 + IT_0102*IT_1096 + IT_0265*IT_1102 + IT_0371*IT_1103 +
       IT_0721*IT_1114 + IT_0722*IT_1116 + IT_0639*IT_1118 + IT_0638*IT_1119 +
       IT_1112*IT_1121) + (-4)*conj(IT_0305)*IT_1131 + (-4)*IT_0305*IT_1132 + (
      -4)*IT_0458*IT_1187 + 4*IT_0452*IT_1196 + (-4)*IT_0454*IT_1197 + (-6)
      *IT_0407*((-0.166666666666667)*IT_0637*IT_0644 + conj(IT_0262)*IT_0893 +
       IT_0896*IT_0922 + IT_0716*IT_1198 + IT_0721*IT_1199 + IT_0723*IT_1200 +
       IT_0639*IT_1201 + IT_0748*IT_1202) + IT_1203 + 4*conj(IT_0458)*IT_1221 + 
      (-4)*conj(IT_0460)*IT_1222 + 4*conj(IT_0452)*IT_1228 + conj(IT_0454)
      *IT_1229 + IT_0368*(IT_0637*IT_0650 + IT_0265*IT_0652 + IT_1230) + IT_0841
      *IT_1239 + IT_0803*IT_1241 + IT_0262*(IT_0511*IT_0634 + IT_1242) + IT_0224
      *IT_1243;
    return create_ccomplex_return(IT_1244);
}

