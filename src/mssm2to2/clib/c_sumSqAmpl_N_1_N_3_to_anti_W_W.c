#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_3_to_anti_W_W.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_3_to_anti_W_W(
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
    const creal_t m_N_3 = param->m_N_3;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_B3 = param->N_B3;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_W3 = param->N_W3;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_d3 = param->N_d3;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t N_u3 = param->N_u3;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t IT_0000 = pow(m_W, 2);
    const ccomplex_t IT_0001 = pow(m_N_1, 2);
    const ccomplex_t IT_0002 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0003 = sin(theta_W);
    const ccomplex_t IT_0004 = cpow(IT_0003, -1);
    const ccomplex_t IT_0005 = conj(N_W1)*e_em*U_Wm1;
    const ccomplex_t IT_0006 = IT_0004*IT_0005;
    const ccomplex_t IT_0007 = conj(N_d1)*U_d1*e_em;
    const ccomplex_t IT_0008 = IT_0004*IT_0007;
    const ccomplex_t IT_0009 = 1.4142135623731*IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*(IT_0006 + 0.5*IT_0009);
    const ccomplex_t IT_0011 = -IT_0010;
    const ccomplex_t IT_0012 = N_W3*e_em*conj(U_Wm1);
    const ccomplex_t IT_0013 = IT_0004*IT_0012;
    const ccomplex_t IT_0014 = N_d3*conj(U_d1)*e_em;
    const ccomplex_t IT_0015 = IT_0004*IT_0014;
    const ccomplex_t IT_0016 = 1.4142135623731*IT_0015;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0013 + 0.5*IT_0016);
    const ccomplex_t IT_0018 = -IT_0017;
    const ccomplex_t IT_0019 = IT_0011*IT_0018;
    const ccomplex_t IT_0020 = IT_0002*IT_0019;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*IT_0020;
    const ccomplex_t IT_0022 = conj(N_W1)*e_em*U_Wm2;
    const ccomplex_t IT_0023 = IT_0004*IT_0022;
    const ccomplex_t IT_0024 = conj(N_d1)*U_d2*e_em;
    const ccomplex_t IT_0025 = IT_0004*IT_0024;
    const ccomplex_t IT_0026 = 1.4142135623731*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0023 + 0.5*IT_0026);
    const ccomplex_t IT_0028 = -IT_0027;
    const ccomplex_t IT_0029 = N_W3*e_em*conj(U_Wm2);
    const ccomplex_t IT_0030 = IT_0004*IT_0029;
    const ccomplex_t IT_0031 = N_d3*conj(U_d2)*e_em;
    const ccomplex_t IT_0032 = IT_0004*IT_0031;
    const ccomplex_t IT_0033 = 1.4142135623731*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*(IT_0030 + 0.5*IT_0033);
    const ccomplex_t IT_0035 = -IT_0034;
    const ccomplex_t IT_0036 = IT_0028*IT_0035;
    const ccomplex_t IT_0037 = cpow((-2)*s_13 + IT_0000 + IT_0001 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
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
    const ccomplex_t IT_0050 = (0 + _Complex_I*1)*e_em*IT_0004*IT_0049;
    const ccomplex_t IT_0051 = -IT_0050;
    const ccomplex_t IT_0052 = cpow(IT_0049, -1);
    const ccomplex_t IT_0053 = IT_0003*IT_0052;
    const ccomplex_t IT_0054 = N_u1*conj(N_u3)*e_em;
    const ccomplex_t IT_0055 = IT_0053*IT_0054;
    const ccomplex_t IT_0056 = IT_0004*IT_0049;
    const ccomplex_t IT_0057 = IT_0054*IT_0056;
    const ccomplex_t IT_0058 = N_d1*conj(N_d3)*e_em;
    const ccomplex_t IT_0059 = IT_0053*IT_0058;
    const ccomplex_t IT_0060 = IT_0056*IT_0058;
    const ccomplex_t IT_0061 = (0 + _Complex_I*1)*(IT_0055 + IT_0057 + 
      -IT_0059 + -IT_0060);
    const ccomplex_t IT_0062 = (-0.5)*IT_0061;
    const ccomplex_t IT_0063 = IT_0051*IT_0062;
    const ccomplex_t IT_0064 = m_N_3*IT_0063;
    const ccomplex_t IT_0065 = conj(N_d1)*N_d3*e_em;
    const ccomplex_t IT_0066 = IT_0053*IT_0065;
    const ccomplex_t IT_0067 = IT_0056*IT_0065;
    const ccomplex_t IT_0068 = conj(N_u1)*N_u3*e_em;
    const ccomplex_t IT_0069 = IT_0053*IT_0068;
    const ccomplex_t IT_0070 = IT_0056*IT_0068;
    const ccomplex_t IT_0071 = (0 + _Complex_I*1)*(IT_0066 + IT_0067 + 
      -IT_0069 + -IT_0070);
    const ccomplex_t IT_0072 = (-0.5)*IT_0071;
    const ccomplex_t IT_0073 = IT_0051*IT_0072;
    const ccomplex_t IT_0074 = m_N_1*IT_0073;
    const ccomplex_t IT_0075 = pow(m_N_3, 2);
    const ccomplex_t IT_0076 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + -IT_0001 + -IT_0075 + -reg_prop, -1);
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = (IT_0064 + -IT_0074)*IT_0077;
    const ccomplex_t IT_0079 = 2*conj(IT_0078);
    const ccomplex_t IT_0080 = cos(alpha);
    const ccomplex_t IT_0081 = IT_0052*IT_0080;
    const ccomplex_t IT_0082 = N_B3*N_u1*e_em;
    const ccomplex_t IT_0083 = IT_0081*IT_0082;
    const ccomplex_t IT_0084 = sin(alpha);
    const ccomplex_t IT_0085 = IT_0052*IT_0084;
    const ccomplex_t IT_0086 = N_B3*N_d1*e_em;
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = N_B1*N_u3*e_em;
    const ccomplex_t IT_0089 = IT_0081*IT_0088;
    const ccomplex_t IT_0090 = N_B1*N_d3*e_em;
    const ccomplex_t IT_0091 = IT_0085*IT_0090;
    const ccomplex_t IT_0092 = IT_0083 + IT_0087 + IT_0089 + IT_0091;
    const ccomplex_t IT_0093 = IT_0004*IT_0080;
    const ccomplex_t IT_0094 = N_W3*N_u1*e_em;
    const ccomplex_t IT_0095 = IT_0093*IT_0094;
    const ccomplex_t IT_0096 = IT_0004*IT_0084;
    const ccomplex_t IT_0097 = N_W3*N_d1*e_em;
    const ccomplex_t IT_0098 = IT_0096*IT_0097;
    const ccomplex_t IT_0099 = N_W1*N_u3*e_em;
    const ccomplex_t IT_0100 = IT_0093*IT_0099;
    const ccomplex_t IT_0101 = N_W1*N_d3*e_em;
    const ccomplex_t IT_0102 = IT_0096*IT_0101;
    const ccomplex_t IT_0103 = -IT_0095 + -IT_0098 + -IT_0100 + -IT_0102;
    const ccomplex_t IT_0104 = IT_0092 + IT_0103;
    const ccomplex_t IT_0105 = (0 + _Complex_I*1)*IT_0104;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = sin(beta);
    const ccomplex_t IT_0108 = cos(beta);
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*(IT_0080
      *IT_0107 + -IT_0084*IT_0108);
    const ccomplex_t IT_0110 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + -IT_0001 + -IT_0075 + -reg_prop, -1);
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
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*m_W*e_em*IT_0004*(IT_0084
      *IT_0107 + IT_0080*IT_0108);
    const ccomplex_t IT_0128 = cpow((-2)*s_12 + -IT_0001 + -IT_0075 + m_H0*
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
    const ccomplex_t IT_0141 = -IT_0001 + -IT_0140;
    const ccomplex_t IT_0142 = IT_0139 + IT_0141;
    const ccomplex_t IT_0143 = m_N_3*IT_0142;
    const ccomplex_t IT_0144 = m_N_1*IT_0063;
    const ccomplex_t IT_0145 = m_N_3*IT_0073;
    const ccomplex_t IT_0146 = IT_0077*(IT_0144 + -IT_0145);
    const ccomplex_t IT_0147 = -IT_0146;
    const ccomplex_t IT_0148 = 2*conj(IT_0147);
    const ccomplex_t IT_0149 = conj(N_B1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0150 = IT_0081*IT_0149;
    const ccomplex_t IT_0151 = conj(N_B1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0152 = IT_0085*IT_0151;
    const ccomplex_t IT_0153 = conj(N_B3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0154 = IT_0081*IT_0153;
    const ccomplex_t IT_0155 = conj(N_B3)*conj(N_d1)*e_em;
    const ccomplex_t IT_0156 = IT_0085*IT_0155;
    const ccomplex_t IT_0157 = IT_0150 + IT_0152 + IT_0154 + IT_0156;
    const ccomplex_t IT_0158 = conj(N_W1)*conj(N_u3)*e_em;
    const ccomplex_t IT_0159 = IT_0093*IT_0158;
    const ccomplex_t IT_0160 = conj(N_W1)*conj(N_d3)*e_em;
    const ccomplex_t IT_0161 = IT_0096*IT_0160;
    const ccomplex_t IT_0162 = conj(N_W3)*conj(N_u1)*e_em;
    const ccomplex_t IT_0163 = IT_0093*IT_0162;
    const ccomplex_t IT_0164 = conj(N_W3)*conj(N_d1)*e_em;
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
    const ccomplex_t IT_0188 = s_13*s_14*s_23*s_34*IT_0041;
    const ccomplex_t IT_0189 = s_24*s_34*IT_0041*IT_0135;
    const ccomplex_t IT_0190 = s_23*IT_0001;
    const ccomplex_t IT_0191 = s_12*s_13;
    const ccomplex_t IT_0192 = s_23*IT_0043*IT_0135;
    const ccomplex_t IT_0193 = (-2)*IT_0192;
    const ccomplex_t IT_0194 = IT_0191 + IT_0193;
    const ccomplex_t IT_0195 = s_24*s_34*IT_0001*IT_0043;
    const ccomplex_t IT_0196 = -IT_0195;
    const ccomplex_t IT_0197 = pow(s_34, 2);
    const ccomplex_t IT_0198 = s_12*s_13*IT_0041*IT_0197;
    const ccomplex_t IT_0199 = -IT_0198;
    const ccomplex_t IT_0200 = IT_0188 + IT_0189 + IT_0190 + IT_0194 + IT_0196
       + IT_0199;
    const ccomplex_t IT_0201 = pow(s_23, 2);
    const ccomplex_t IT_0202 = s_14*s_34*IT_0041*IT_0201;
    const ccomplex_t IT_0203 = s_13*s_23*s_24*s_34*IT_0041;
    const ccomplex_t IT_0204 = s_13*IT_0075;
    const ccomplex_t IT_0205 = s_12*s_23;
    const ccomplex_t IT_0206 = s_13*IT_0043*IT_0201;
    const ccomplex_t IT_0207 = (-2)*IT_0206;
    const ccomplex_t IT_0208 = IT_0205 + IT_0207;
    const ccomplex_t IT_0209 = s_14*s_34*IT_0043*IT_0075;
    const ccomplex_t IT_0210 = -IT_0209;
    const ccomplex_t IT_0211 = s_12*s_23*IT_0041*IT_0197;
    const ccomplex_t IT_0212 = -IT_0211;
    const ccomplex_t IT_0213 = IT_0202 + IT_0203 + IT_0204 + IT_0208 + IT_0210
       + IT_0212;
    const ccomplex_t IT_0214 = 2*IT_0050;
    const ccomplex_t IT_0215 = IT_0072*IT_0214;
    const ccomplex_t IT_0216 = IT_0077*IT_0215;
    const ccomplex_t IT_0217 = -IT_0216;
    const ccomplex_t IT_0218 = IT_0050*IT_0072;
    const ccomplex_t IT_0219 = IT_0077*IT_0218;
    const ccomplex_t IT_0220 = 2*IT_0219;
    const ccomplex_t IT_0221 = conj(IT_0217)*IT_0220;
    const ccomplex_t IT_0222 = 2*conj(IT_0219);
    const ccomplex_t IT_0223 = IT_0217*IT_0222;
    const ccomplex_t IT_0224 = s_34*m_N_1*m_N_3*IT_0043;
    const ccomplex_t IT_0225 = s_13*s_34*IT_0043;
    const ccomplex_t IT_0226 = -IT_0225;
    const ccomplex_t IT_0227 = s_14 + IT_0226;
    const ccomplex_t IT_0228 = IT_0224*IT_0227;
    const ccomplex_t IT_0229 = s_23*s_34*IT_0043;
    const ccomplex_t IT_0230 = -IT_0229;
    const ccomplex_t IT_0231 = s_24 + IT_0230;
    const ccomplex_t IT_0232 = IT_0224*IT_0231;
    const ccomplex_t IT_0233 = IT_0228 + IT_0232;
    const ccomplex_t IT_0234 = (-2)*IT_0050;
    const ccomplex_t IT_0235 = IT_0062*IT_0234;
    const ccomplex_t IT_0236 = IT_0077*IT_0235;
    const ccomplex_t IT_0237 = 2*IT_0217;
    const ccomplex_t IT_0238 = conj(IT_0236)*IT_0237;
    const ccomplex_t IT_0239 = 2*conj(IT_0217);
    const ccomplex_t IT_0240 = IT_0236*IT_0239;
    const ccomplex_t IT_0241 = conj(N_W1)*e_em*V_Wp1;
    const ccomplex_t IT_0242 = IT_0004*IT_0241;
    const ccomplex_t IT_0243 = conj(N_u1)*V_u1*e_em;
    const ccomplex_t IT_0244 = IT_0004*IT_0243;
    const ccomplex_t IT_0245 = 1.4142135623731*IT_0244;
    const ccomplex_t IT_0246 = (0 + _Complex_I*1)*(IT_0242 + (-0.5)*IT_0245);
    const ccomplex_t IT_0247 = -IT_0246;
    const ccomplex_t IT_0248 = N_W3*e_em*conj(V_Wp1);
    const ccomplex_t IT_0249 = IT_0004*IT_0248;
    const ccomplex_t IT_0250 = N_u3*conj(V_u1)*e_em;
    const ccomplex_t IT_0251 = IT_0004*IT_0250;
    const ccomplex_t IT_0252 = 1.4142135623731*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*(IT_0249 + (-0.5)*IT_0252);
    const ccomplex_t IT_0254 = -IT_0253;
    const ccomplex_t IT_0255 = IT_0247*IT_0254;
    const ccomplex_t IT_0256 = cpow((-2)*s_23 + IT_0000 + IT_0075 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0257 = IT_0255*IT_0256;
    const ccomplex_t IT_0258 = (0 + _Complex_I*1)*IT_0257;
    const ccomplex_t IT_0259 = conj(N_W1)*e_em*V_Wp2;
    const ccomplex_t IT_0260 = IT_0004*IT_0259;
    const ccomplex_t IT_0261 = conj(N_u1)*V_u2*e_em;
    const ccomplex_t IT_0262 = IT_0004*IT_0261;
    const ccomplex_t IT_0263 = 1.4142135623731*IT_0262;
    const ccomplex_t IT_0264 = (0 + _Complex_I*1)*(IT_0260 + (-0.5)*IT_0263);
    const ccomplex_t IT_0265 = -IT_0264;
    const ccomplex_t IT_0266 = N_W3*e_em*conj(V_Wp2);
    const ccomplex_t IT_0267 = IT_0004*IT_0266;
    const ccomplex_t IT_0268 = N_u3*conj(V_u2)*e_em;
    const ccomplex_t IT_0269 = IT_0004*IT_0268;
    const ccomplex_t IT_0270 = 1.4142135623731*IT_0269;
    const ccomplex_t IT_0271 = (0 + _Complex_I*1)*(IT_0267 + (-0.5)*IT_0270);
    const ccomplex_t IT_0272 = -IT_0271;
    const ccomplex_t IT_0273 = IT_0265*IT_0272;
    const ccomplex_t IT_0274 = cpow((-2)*s_23 + IT_0000 + IT_0075 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0275 = IT_0273*IT_0274;
    const ccomplex_t IT_0276 = (0 + _Complex_I*1)*IT_0275;
    const ccomplex_t IT_0277 = 2*IT_0258 + 2*IT_0276;
    const ccomplex_t IT_0278 = IT_0043*IT_0201;
    const ccomplex_t IT_0279 = pow(s_24, 2);
    const ccomplex_t IT_0280 = IT_0043*IT_0279;
    const ccomplex_t IT_0281 = IT_0278 + IT_0280;
    const ccomplex_t IT_0282 = s_23*s_24*s_34*IT_0041;
    const ccomplex_t IT_0283 = -IT_0075 + -IT_0282;
    const ccomplex_t IT_0284 = IT_0281 + IT_0283;
    const ccomplex_t IT_0285 = m_N_1*IT_0284;
    const ccomplex_t IT_0286 = 2*conj(IT_0132);
    const ccomplex_t IT_0287 = (-2)*conj(IT_0078);
    const ccomplex_t IT_0288 = IT_0286 + IT_0287;
    const ccomplex_t IT_0289 = s_14*s_23*s_34*IT_0041;
    const ccomplex_t IT_0290 = s_12 + IT_0046 + IT_0289;
    const ccomplex_t IT_0291 = m_N_3*IT_0290;
    const ccomplex_t IT_0292 = 2*conj(IT_0185);
    const ccomplex_t IT_0293 = (-2)*conj(IT_0147);
    const ccomplex_t IT_0294 = IT_0292 + IT_0293;
    const ccomplex_t IT_0295 = s_12*IT_0000;
    const ccomplex_t IT_0296 = s_14*s_24;
    const ccomplex_t IT_0297 = (-2)*IT_0296;
    const ccomplex_t IT_0298 = IT_0295 + IT_0297;
    const ccomplex_t IT_0299 = IT_0043*IT_0298;
    const ccomplex_t IT_0300 = (s_12 + -IT_0044)*(s_12 + (-0.5)*IT_0299);
    const ccomplex_t IT_0301 = (-2)*IT_0300;
    const ccomplex_t IT_0302 = 2*conj(IT_0040);
    const ccomplex_t IT_0303 = s_13*s_23*IT_0041*IT_0137;
    const ccomplex_t IT_0304 = s_14*s_24*IT_0041*IT_0135;
    const ccomplex_t IT_0305 = s_12*IT_0001;
    const ccomplex_t IT_0306 = s_14*s_24*IT_0001*IT_0043;
    const ccomplex_t IT_0307 = -IT_0306;
    const ccomplex_t IT_0308 = s_13*s_23*IT_0001*IT_0043;
    const ccomplex_t IT_0309 = -IT_0308;
    const ccomplex_t IT_0310 = s_12*s_13*s_14*s_34*IT_0041;
    const ccomplex_t IT_0311 = -IT_0310;
    const ccomplex_t IT_0312 = IT_0303 + IT_0304 + IT_0305 + IT_0307 + IT_0309
       + IT_0311;
    const ccomplex_t IT_0313 = s_13*s_14*s_23*s_24*IT_0041;
    const ccomplex_t IT_0314 = IT_0041*IT_0135*IT_0279;
    const ccomplex_t IT_0315 = IT_0001*IT_0075;
    const ccomplex_t IT_0316 = IT_0043*IT_0075*IT_0135;
    const ccomplex_t IT_0317 = -IT_0316;
    const ccomplex_t IT_0318 = IT_0001*IT_0043*IT_0279;
    const ccomplex_t IT_0319 = -IT_0318;
    const ccomplex_t IT_0320 = s_12*s_13*s_24*s_34*IT_0041;
    const ccomplex_t IT_0321 = -IT_0320;
    const ccomplex_t IT_0322 = IT_0313 + IT_0314 + IT_0315 + IT_0317 + IT_0319
       + IT_0321;
    const ccomplex_t IT_0323 = IT_0073*IT_0077;
    const ccomplex_t IT_0324 = -IT_0323;
    const ccomplex_t IT_0325 = 2*IT_0324;
    const ccomplex_t IT_0326 = conj(IT_0217)*IT_0325;
    const ccomplex_t IT_0327 = 2*conj(IT_0324);
    const ccomplex_t IT_0328 = IT_0217*IT_0327;
    const ccomplex_t IT_0329 = IT_0326 + IT_0328;
    const ccomplex_t IT_0330 = 2*IT_0040;
    const ccomplex_t IT_0331 = conj(IT_0324)*IT_0330;
    const ccomplex_t IT_0332 = IT_0302*IT_0324;
    const ccomplex_t IT_0333 = N_W1*e_em*conj(V_Wp1);
    const ccomplex_t IT_0334 = IT_0004*IT_0333;
    const ccomplex_t IT_0335 = N_u1*conj(V_u1)*e_em;
    const ccomplex_t IT_0336 = IT_0004*IT_0335;
    const ccomplex_t IT_0337 = 1.4142135623731*IT_0336;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*(IT_0334 + (-0.5)*IT_0337);
    const ccomplex_t IT_0339 = -IT_0338;
    const ccomplex_t IT_0340 = conj(N_W3)*e_em*V_Wp1;
    const ccomplex_t IT_0341 = IT_0004*IT_0340;
    const ccomplex_t IT_0342 = conj(N_u3)*V_u1*e_em;
    const ccomplex_t IT_0343 = IT_0004*IT_0342;
    const ccomplex_t IT_0344 = 1.4142135623731*IT_0343;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*(IT_0341 + (-0.5)*IT_0344);
    const ccomplex_t IT_0346 = -IT_0345;
    const ccomplex_t IT_0347 = IT_0339*IT_0346;
    const ccomplex_t IT_0348 = IT_0002*IT_0347;
    const ccomplex_t IT_0349 = (0 + _Complex_I*1)*IT_0348;
    const ccomplex_t IT_0350 = IT_0062*IT_0214;
    const ccomplex_t IT_0351 = IT_0077*IT_0350;
    const ccomplex_t IT_0352 = N_W1*e_em*conj(V_Wp2);
    const ccomplex_t IT_0353 = IT_0004*IT_0352;
    const ccomplex_t IT_0354 = N_u1*conj(V_u2)*e_em;
    const ccomplex_t IT_0355 = IT_0004*IT_0354;
    const ccomplex_t IT_0356 = 1.4142135623731*IT_0355;
    const ccomplex_t IT_0357 = (0 + _Complex_I*1)*(IT_0353 + (-0.5)*IT_0356);
    const ccomplex_t IT_0358 = -IT_0357;
    const ccomplex_t IT_0359 = conj(N_W3)*e_em*V_Wp2;
    const ccomplex_t IT_0360 = IT_0004*IT_0359;
    const ccomplex_t IT_0361 = conj(N_u3)*V_u2*e_em;
    const ccomplex_t IT_0362 = IT_0004*IT_0361;
    const ccomplex_t IT_0363 = 1.4142135623731*IT_0362;
    const ccomplex_t IT_0364 = (0 + _Complex_I*1)*(IT_0360 + (-0.5)*IT_0363);
    const ccomplex_t IT_0365 = -IT_0364;
    const ccomplex_t IT_0366 = IT_0358*IT_0365;
    const ccomplex_t IT_0367 = IT_0037*IT_0366;
    const ccomplex_t IT_0368 = (0 + _Complex_I*1)*IT_0367;
    const ccomplex_t IT_0369 = (-2)*IT_0349 + -IT_0351 + (-2)*IT_0368;
    const ccomplex_t IT_0370 = IT_0063*IT_0077;
    const ccomplex_t IT_0371 = -IT_0370;
    const ccomplex_t IT_0372 = 2*IT_0371;
    const ccomplex_t IT_0373 = conj(IT_0369)*IT_0372;
    const ccomplex_t IT_0374 = 2*conj(IT_0371);
    const ccomplex_t IT_0375 = IT_0369*IT_0374;
    const ccomplex_t IT_0376 = IT_0041*IT_0137*IT_0201;
    const ccomplex_t IT_0377 = IT_0043*IT_0075*IT_0137;
    const ccomplex_t IT_0378 = -IT_0377;
    const ccomplex_t IT_0379 = IT_0001*IT_0043*IT_0201;
    const ccomplex_t IT_0380 = -IT_0379;
    const ccomplex_t IT_0381 = s_12*s_14*s_23*s_34*IT_0041;
    const ccomplex_t IT_0382 = -IT_0381;
    const ccomplex_t IT_0383 = IT_0313 + IT_0315 + IT_0376 + IT_0378 + IT_0380
       + IT_0382;
    const ccomplex_t IT_0384 = s_14*s_24*IT_0041*IT_0201;
    const ccomplex_t IT_0385 = s_13*s_23*IT_0041*IT_0279;
    const ccomplex_t IT_0386 = s_12*IT_0075;
    const ccomplex_t IT_0387 = s_14*s_24*IT_0043*IT_0075;
    const ccomplex_t IT_0388 = -IT_0387;
    const ccomplex_t IT_0389 = s_13*s_23*IT_0043*IT_0075;
    const ccomplex_t IT_0390 = -IT_0389;
    const ccomplex_t IT_0391 = s_12*s_23*s_24*s_34*IT_0041;
    const ccomplex_t IT_0392 = -IT_0391;
    const ccomplex_t IT_0393 = IT_0384 + IT_0385 + IT_0386 + IT_0388 + IT_0390
       + IT_0392;
    const ccomplex_t IT_0394 = IT_0383 + IT_0393;
    const ccomplex_t IT_0395 = conj(N_W3)*e_em*U_Wm1;
    const ccomplex_t IT_0396 = IT_0004*IT_0395;
    const ccomplex_t IT_0397 = conj(N_d3)*U_d1*e_em;
    const ccomplex_t IT_0398 = IT_0004*IT_0397;
    const ccomplex_t IT_0399 = 1.4142135623731*IT_0398;
    const ccomplex_t IT_0400 = (0 + _Complex_I*1)*(IT_0396 + 0.5*IT_0399);
    const ccomplex_t IT_0401 = -IT_0400;
    const ccomplex_t IT_0402 = N_W1*e_em*conj(U_Wm1);
    const ccomplex_t IT_0403 = IT_0004*IT_0402;
    const ccomplex_t IT_0404 = N_d1*conj(U_d1)*e_em;
    const ccomplex_t IT_0405 = IT_0004*IT_0404;
    const ccomplex_t IT_0406 = 1.4142135623731*IT_0405;
    const ccomplex_t IT_0407 = (0 + _Complex_I*1)*(IT_0403 + 0.5*IT_0406);
    const ccomplex_t IT_0408 = -IT_0407;
    const ccomplex_t IT_0409 = IT_0401*IT_0408;
    const ccomplex_t IT_0410 = IT_0256*IT_0409;
    const ccomplex_t IT_0411 = (0 + _Complex_I*1)*IT_0410;
    const ccomplex_t IT_0412 = N_W1*e_em*conj(U_Wm2);
    const ccomplex_t IT_0413 = IT_0004*IT_0412;
    const ccomplex_t IT_0414 = N_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0415 = IT_0004*IT_0414;
    const ccomplex_t IT_0416 = 1.4142135623731*IT_0415;
    const ccomplex_t IT_0417 = (0 + _Complex_I*1)*(IT_0413 + 0.5*IT_0416);
    const ccomplex_t IT_0418 = -IT_0417;
    const ccomplex_t IT_0419 = conj(N_W3)*e_em*U_Wm2;
    const ccomplex_t IT_0420 = IT_0004*IT_0419;
    const ccomplex_t IT_0421 = conj(N_d3)*U_d2*e_em;
    const ccomplex_t IT_0422 = IT_0004*IT_0421;
    const ccomplex_t IT_0423 = 1.4142135623731*IT_0422;
    const ccomplex_t IT_0424 = (0 + _Complex_I*1)*(IT_0420 + 0.5*IT_0423);
    const ccomplex_t IT_0425 = -IT_0424;
    const ccomplex_t IT_0426 = IT_0418*IT_0425;
    const ccomplex_t IT_0427 = IT_0274*IT_0426;
    const ccomplex_t IT_0428 = (0 + _Complex_I*1)*IT_0427;
    const ccomplex_t IT_0429 = -IT_0351 + 2*IT_0411 + 2*IT_0428;
    const ccomplex_t IT_0430 = IT_0372*conj(IT_0429);
    const ccomplex_t IT_0431 = IT_0374*IT_0429;
    const ccomplex_t IT_0432 = s_14*s_34*IT_0043;
    const ccomplex_t IT_0433 = s_13*s_23;
    const ccomplex_t IT_0434 = (-2)*IT_0433;
    const ccomplex_t IT_0435 = IT_0295 + IT_0434;
    const ccomplex_t IT_0436 = IT_0043*IT_0435;
    const ccomplex_t IT_0437 = (s_13 + -IT_0432)*(s_12 + (-0.5)*IT_0436);
    const ccomplex_t IT_0438 = (-2)*IT_0437;
    const ccomplex_t IT_0439 = s_24*s_34*IT_0043;
    const ccomplex_t IT_0440 = (s_12 + (-0.5)*IT_0436)*(s_23 + -IT_0439);
    const ccomplex_t IT_0441 = (-2)*IT_0440;
    const ccomplex_t IT_0442 = IT_0438 + IT_0441;
    const ccomplex_t IT_0443 = IT_0050*IT_0062;
    const ccomplex_t IT_0444 = IT_0077*IT_0443;
    const ccomplex_t IT_0445 = IT_0372*conj(IT_0444);
    const ccomplex_t IT_0446 = IT_0374*IT_0444;
    const ccomplex_t IT_0447 = s_13*s_14*s_24*IT_0043;
    const ccomplex_t IT_0448 = -IT_0447;
    const ccomplex_t IT_0449 = s_23*IT_0043*IT_0137;
    const ccomplex_t IT_0450 = -IT_0449;
    const ccomplex_t IT_0451 = 2*IT_0188;
    const ccomplex_t IT_0452 = IT_0190 + IT_0194 + IT_0448 + IT_0450 + IT_0451;
    const ccomplex_t IT_0453 = s_13*IT_0043*IT_0279;
    const ccomplex_t IT_0454 = -IT_0453;
    const ccomplex_t IT_0455 = s_14*s_23*s_24*IT_0043;
    const ccomplex_t IT_0456 = -IT_0455;
    const ccomplex_t IT_0457 = 2*IT_0203;
    const ccomplex_t IT_0458 = IT_0204 + IT_0208 + IT_0454 + IT_0456 + IT_0457;
    const ccomplex_t IT_0459 = IT_0072*IT_0234;
    const ccomplex_t IT_0460 = IT_0077*IT_0459;
    const ccomplex_t IT_0461 = IT_0325*conj(IT_0460);
    const ccomplex_t IT_0462 = IT_0327*IT_0460;
    const ccomplex_t IT_0463 = conj(IT_0236)*IT_0372;
    const ccomplex_t IT_0464 = IT_0236*IT_0374;
    const ccomplex_t IT_0465 = m_N_1*m_N_3;
    const ccomplex_t IT_0466 = s_13*IT_0000;
    const ccomplex_t IT_0467 = s_14*s_34;
    const ccomplex_t IT_0468 = -IT_0467;
    const ccomplex_t IT_0469 = IT_0466 + IT_0468;
    const ccomplex_t IT_0470 = IT_0043*IT_0469;
    const ccomplex_t IT_0471 = -IT_0432 + (-0.25)*IT_0470;
    const ccomplex_t IT_0472 = s_13 + IT_0471;
    const ccomplex_t IT_0473 = IT_0465*IT_0472;
    const ccomplex_t IT_0474 = (-4)*IT_0473;
    const ccomplex_t IT_0475 = s_23*IT_0000;
    const ccomplex_t IT_0476 = s_24*s_34;
    const ccomplex_t IT_0477 = -IT_0476;
    const ccomplex_t IT_0478 = IT_0475 + IT_0477;
    const ccomplex_t IT_0479 = IT_0043*IT_0478;
    const ccomplex_t IT_0480 = -IT_0439 + (-0.25)*IT_0479;
    const ccomplex_t IT_0481 = s_23 + IT_0480;
    const ccomplex_t IT_0482 = IT_0465*IT_0481;
    const ccomplex_t IT_0483 = (-4)*IT_0482;
    const ccomplex_t IT_0484 = IT_0349 + IT_0368;
    const ccomplex_t IT_0485 = IT_0325*conj(IT_0484);
    const ccomplex_t IT_0486 = -IT_0411 + -IT_0428;
    const ccomplex_t IT_0487 = IT_0325*conj(IT_0486);
    const ccomplex_t IT_0488 = IT_0327*IT_0484;
    const ccomplex_t IT_0489 = IT_0327*IT_0486;
    const ccomplex_t IT_0490 = IT_0021 + IT_0039;
    const ccomplex_t IT_0491 = IT_0372*conj(IT_0490);
    const ccomplex_t IT_0492 = -IT_0258 + -IT_0276;
    const ccomplex_t IT_0493 = IT_0372*conj(IT_0492);
    const ccomplex_t IT_0494 = IT_0374*IT_0490;
    const ccomplex_t IT_0495 = IT_0374*IT_0492;
    const ccomplex_t IT_0496 = -IT_0001;
    const ccomplex_t IT_0497 = IT_0138 + IT_0496;
    const ccomplex_t IT_0498 = IT_0465*IT_0497;
    const ccomplex_t IT_0499 = -IT_0075;
    const ccomplex_t IT_0500 = IT_0280 + IT_0499;
    const ccomplex_t IT_0501 = IT_0465*IT_0500;
    const ccomplex_t IT_0502 = conj(IT_0324)*IT_0371;
    const ccomplex_t IT_0503 = 6*IT_0502;
    const ccomplex_t IT_0504 = IT_0324*conj(IT_0371);
    const ccomplex_t IT_0505 = 6*IT_0504;
    const ccomplex_t IT_0506 = -IT_0432;
    const ccomplex_t IT_0507 = s_13 + IT_0506;
    const ccomplex_t IT_0508 = IT_0465*IT_0507;
    const ccomplex_t IT_0509 = -IT_0439;
    const ccomplex_t IT_0510 = s_23 + IT_0509;
    const ccomplex_t IT_0511 = IT_0465*IT_0510;
    const ccomplex_t IT_0512 = conj(IT_0219)*IT_0371;
    const ccomplex_t IT_0513 = (-6)*IT_0512;
    const ccomplex_t IT_0514 = IT_0219*conj(IT_0371);
    const ccomplex_t IT_0515 = (-6)*IT_0514;
    const ccomplex_t IT_0516 = (-6)*IT_0444;
    const ccomplex_t IT_0517 = conj(IT_0324)*IT_0516;
    const ccomplex_t IT_0518 = (-6)*conj(IT_0444);
    const ccomplex_t IT_0519 = IT_0324*IT_0518;
    const ccomplex_t IT_0520 = IT_0041*IT_0197;
    const ccomplex_t IT_0521 = 0.5*IT_0520;
    const ccomplex_t IT_0522 = 1 + IT_0521;
    const ccomplex_t IT_0523 = s_12*IT_0522;
    const ccomplex_t IT_0524 = (-4)*conj(IT_0132);
    const ccomplex_t IT_0525 = 0.25*IT_0524;
    const ccomplex_t IT_0526 = conj(IT_0078) + IT_0525;
    const ccomplex_t IT_0527 = (-4)*conj(IT_0185);
    const ccomplex_t IT_0528 = 0.25*IT_0527;
    const ccomplex_t IT_0529 = conj(IT_0147) + IT_0528;
    const ccomplex_t IT_0530 = (-4)*IT_0132;
    const ccomplex_t IT_0531 = (-4)*IT_0185;
    const ccomplex_t IT_0532 = 2*IT_0132;
    const ccomplex_t IT_0533 = (-2)*IT_0078;
    const ccomplex_t IT_0534 = IT_0532 + IT_0533;
    const ccomplex_t IT_0535 = 2*IT_0185;
    const ccomplex_t IT_0536 = (-2)*IT_0147;
    const ccomplex_t IT_0537 = IT_0535 + IT_0536;
    const ccomplex_t IT_0538 = -(IT_0075 + -IT_0278)*(s_12 + (-0.5)*IT_0299);
    const ccomplex_t IT_0539 = (-2)*IT_0538;
    const ccomplex_t IT_0540 = (-2)*IT_0277;
    const ccomplex_t IT_0541 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0542 = IT_0254*IT_0408;
    const ccomplex_t IT_0543 = IT_0256*IT_0541*IT_0542;
    const ccomplex_t IT_0544 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0545 = IT_0272*IT_0418;
    const ccomplex_t IT_0546 = IT_0274*IT_0544*IT_0545;
    const ccomplex_t IT_0547 = IT_0543 + IT_0546;
    const ccomplex_t IT_0548 = -IT_0433;
    const ccomplex_t IT_0549 = IT_0295 + IT_0548;
    const ccomplex_t IT_0550 = IT_0043*IT_0549;
    const ccomplex_t IT_0551 = -IT_0044 + (-0.25)*IT_0550;
    const ccomplex_t IT_0552 = s_12 + IT_0551;
    const ccomplex_t IT_0553 = m_N_1*IT_0552;
    const ccomplex_t IT_0554 = (-4)*IT_0553;
    const ccomplex_t IT_0555 = pow(m_W, 4);
    const ccomplex_t IT_0556 = s_12*IT_0555;
    const ccomplex_t IT_0557 = s_13*s_23*IT_0000;
    const ccomplex_t IT_0558 = s_14*s_24*IT_0000;
    const ccomplex_t IT_0559 = s_13*s_24*s_34;
    const ccomplex_t IT_0560 = s_14*s_23*s_34;
    const ccomplex_t IT_0561 = s_12*IT_0197;
    const ccomplex_t IT_0562 = (-2)*IT_0557 + (-2)*IT_0558 + 2*IT_0559 + 2
      *IT_0560 + (-2)*IT_0561;
    const ccomplex_t IT_0563 = IT_0556 + IT_0562;
    const ccomplex_t IT_0564 = IT_0041*IT_0563;
    const ccomplex_t IT_0565 = (-0.25)*IT_0564;
    const ccomplex_t IT_0566 = -IT_0045;
    const ccomplex_t IT_0567 = -IT_0044;
    const ccomplex_t IT_0568 = s_12 + IT_0565 + IT_0566 + IT_0567;
    const ccomplex_t IT_0569 = IT_0018*IT_0339;
    const ccomplex_t IT_0570 = IT_0002*IT_0541*IT_0569;
    const ccomplex_t IT_0571 = IT_0035*IT_0358;
    const ccomplex_t IT_0572 = IT_0037*IT_0544*IT_0571;
    const ccomplex_t IT_0573 = m_N_1*IT_0021;
    const ccomplex_t IT_0574 = m_N_1*IT_0039;
    const ccomplex_t IT_0575 = IT_0570 + IT_0572 + -IT_0573 + -IT_0574;
    const ccomplex_t IT_0576 = 8*conj(IT_0575);
    const ccomplex_t IT_0577 = IT_0000*IT_0075;
    const ccomplex_t IT_0578 = -IT_0577;
    const ccomplex_t IT_0579 = IT_0201 + IT_0578;
    const ccomplex_t IT_0580 = IT_0043*IT_0579;
    const ccomplex_t IT_0581 = IT_0278 + IT_0499;
    const ccomplex_t IT_0582 = (-4)*IT_0581;
    const ccomplex_t IT_0583 = IT_0580 + IT_0582;
    const ccomplex_t IT_0584 = m_N_1*IT_0583;
    const ccomplex_t IT_0585 = (-2)*conj(IT_0277);
    const ccomplex_t IT_0586 = 2*IT_0078;
    const ccomplex_t IT_0587 = (-2)*IT_0132;
    const ccomplex_t IT_0588 = IT_0586 + IT_0587;
    const ccomplex_t IT_0589 = 2*IT_0147;
    const ccomplex_t IT_0590 = (-2)*IT_0185;
    const ccomplex_t IT_0591 = IT_0589 + IT_0590;
    const ccomplex_t IT_0592 = -(IT_0001 + -IT_0136)*(s_12 + (-0.5)*IT_0299);
    const ccomplex_t IT_0593 = (-2)*IT_0592;
    const ccomplex_t IT_0594 = (-2)*IT_0040;
    const ccomplex_t IT_0595 = IT_0043*IT_0295;
    const ccomplex_t IT_0596 = (-0.5)*IT_0595;
    const ccomplex_t IT_0597 = s_34*IT_0041;
    const ccomplex_t IT_0598 = s_14*s_23;
    const ccomplex_t IT_0599 = s_13*s_24;
    const ccomplex_t IT_0600 = s_12*s_34;
    const ccomplex_t IT_0601 = -IT_0599 + -IT_0600;
    const ccomplex_t IT_0602 = IT_0598 + IT_0601;
    const ccomplex_t IT_0603 = IT_0597*IT_0602;
    const ccomplex_t IT_0604 = (-0.25)*IT_0603;
    const ccomplex_t IT_0605 = s_12 + IT_0596 + IT_0604;
    const ccomplex_t IT_0606 = (-8)*IT_0078;
    const ccomplex_t IT_0607 = 0.125*IT_0606;
    const ccomplex_t IT_0608 = IT_0132 + IT_0607;
    const ccomplex_t IT_0609 = s_23*s_24*s_34*IT_0043;
    const ccomplex_t IT_0610 = 2*IT_0279 + -IT_0577 + (-2)*IT_0609;
    const ccomplex_t IT_0611 = IT_0201 + IT_0610;
    const ccomplex_t IT_0612 = IT_0043*IT_0611;
    const ccomplex_t IT_0613 = -IT_0075 + (-0.5)*IT_0612;
    const ccomplex_t IT_0614 = IT_0278 + IT_0613;
    const ccomplex_t IT_0615 = m_N_1*IT_0614;
    const ccomplex_t IT_0616 = (-2)*IT_0615;
    const ccomplex_t IT_0617 = 2*IT_0277;
    const ccomplex_t IT_0618 = (-0.5)*IT_0436;
    const ccomplex_t IT_0619 = s_12 + IT_0618;
    const ccomplex_t IT_0620 = (-2)*IT_0557 + (-2)*IT_0558 + 4*IT_0559;
    const ccomplex_t IT_0621 = IT_0556 + IT_0620;
    const ccomplex_t IT_0622 = IT_0041*IT_0621;
    const ccomplex_t IT_0623 = 0.25*IT_0622;
    const ccomplex_t IT_0624 = (-0.5)*IT_0299;
    const ccomplex_t IT_0625 = 8*IT_0575;
    const ccomplex_t IT_0626 = s_13*s_24*s_34*IT_0043;
    const ccomplex_t IT_0627 = (-2)*IT_0296 + -IT_0433 + 2*IT_0626;
    const ccomplex_t IT_0628 = IT_0295 + IT_0627;
    const ccomplex_t IT_0629 = IT_0043*IT_0628;
    const ccomplex_t IT_0630 = -IT_0044 + (-0.5)*IT_0629;
    const ccomplex_t IT_0631 = s_12 + IT_0630;
    const ccomplex_t IT_0632 = m_N_1*IT_0631;
    const ccomplex_t IT_0633 = (-2)*IT_0632;
    const ccomplex_t IT_0634 = 8*conj(IT_0132);
    const ccomplex_t IT_0635 = (-8)*conj(IT_0078);
    const ccomplex_t IT_0636 = IT_0634 + IT_0635;
    const ccomplex_t IT_0637 = 2*conj(IT_0277);
    const ccomplex_t IT_0638 = (-2)*conj(IT_0040);
    const ccomplex_t IT_0639 = IT_0258 + IT_0276;
    const ccomplex_t IT_0640 = m_N_3*IT_0639;
    const ccomplex_t IT_0641 = -IT_0640;
    const ccomplex_t IT_0642 = s_14*s_23*s_34*IT_0043;
    const ccomplex_t IT_0643 = (-2)*IT_0296 + -IT_0433 + 2*IT_0642;
    const ccomplex_t IT_0644 = IT_0295 + IT_0643;
    const ccomplex_t IT_0645 = IT_0043*IT_0644;
    const ccomplex_t IT_0646 = -IT_0044 + (-0.5)*IT_0645;
    const ccomplex_t IT_0647 = s_12 + IT_0646;
    const ccomplex_t IT_0648 = m_N_3*IT_0647;
    const ccomplex_t IT_0649 = (-2)*IT_0648;
    const ccomplex_t IT_0650 = (-2)*IT_0197;
    const ccomplex_t IT_0651 = IT_0555 + IT_0650;
    const ccomplex_t IT_0652 = IT_0041*IT_0651;
    const ccomplex_t IT_0653 = 0.25*IT_0652;
    const ccomplex_t IT_0654 = 1 + IT_0653;
    const ccomplex_t IT_0655 = IT_0465*IT_0654;
    const ccomplex_t IT_0656 = s_13*s_14*s_34*IT_0043;
    const ccomplex_t IT_0657 = IT_0000*IT_0001;
    const ccomplex_t IT_0658 = 2*IT_0137 + (-2)*IT_0656 + -IT_0657;
    const ccomplex_t IT_0659 = IT_0135 + IT_0658;
    const ccomplex_t IT_0660 = IT_0043*IT_0659;
    const ccomplex_t IT_0661 = -IT_0001 + (-0.5)*IT_0660;
    const ccomplex_t IT_0662 = IT_0136 + IT_0661;
    const ccomplex_t IT_0663 = m_N_3*IT_0662;
    const ccomplex_t IT_0664 = (-2)*IT_0663;
    const ccomplex_t IT_0665 = -(IT_0001 + -IT_0138)*(s_12 + (-0.5)*IT_0436);
    const ccomplex_t IT_0666 = (-2)*IT_0665;
    const ccomplex_t IT_0667 = -(IT_0075 + -IT_0280)*(s_12 + (-0.5)*IT_0436);
    const ccomplex_t IT_0668 = (-2)*IT_0667;
    const ccomplex_t IT_0669 = IT_0666 + IT_0668;
    const ccomplex_t IT_0670 = (-2)*IT_0371;
    const ccomplex_t IT_0671 = conj(IT_0371)*IT_0670;
    const ccomplex_t IT_0672 = IT_0047*IT_0465;
    const ccomplex_t IT_0673 = IT_0290*IT_0465;
    const ccomplex_t IT_0674 = conj(IT_0217)*IT_0670;
    const ccomplex_t IT_0675 = (-2)*conj(IT_0371);
    const ccomplex_t IT_0676 = IT_0217*IT_0675;
    const ccomplex_t IT_0677 = s_34*m_N_1*IT_0043;
    const ccomplex_t IT_0678 = IT_0231*IT_0677;
    const ccomplex_t IT_0679 = s_34*m_N_3*IT_0043;
    const ccomplex_t IT_0680 = IT_0227*IT_0679;
    const ccomplex_t IT_0681 = IT_0475 + IT_0476;
    const ccomplex_t IT_0682 = s_23*IT_0043*IT_0197;
    const ccomplex_t IT_0683 = (-2)*IT_0682;
    const ccomplex_t IT_0684 = IT_0681 + IT_0683;
    const ccomplex_t IT_0685 = IT_0043*IT_0684;
    const ccomplex_t IT_0686 = -IT_0439 + 0.5*IT_0685;
    const ccomplex_t IT_0687 = s_23 + IT_0686;
    const ccomplex_t IT_0688 = m_N_1*IT_0687;
    const ccomplex_t IT_0689 = 2*IT_0688;
    const ccomplex_t IT_0690 = 2*IT_0641;
    const ccomplex_t IT_0691 = IT_0043*IT_0197;
    const ccomplex_t IT_0692 = (s_12 + (-0.5)*IT_0436)*(IT_0000 + -IT_0691);
    const ccomplex_t IT_0693 = (-2)*IT_0692;
    const ccomplex_t IT_0694 = 2*IT_0444;
    const ccomplex_t IT_0695 = IT_0466 + IT_0467;
    const ccomplex_t IT_0696 = s_13*IT_0043*IT_0197;
    const ccomplex_t IT_0697 = (-2)*IT_0696;
    const ccomplex_t IT_0698 = IT_0695 + IT_0697;
    const ccomplex_t IT_0699 = IT_0043*IT_0698;
    const ccomplex_t IT_0700 = -IT_0432 + 0.5*IT_0699;
    const ccomplex_t IT_0701 = s_13 + IT_0700;
    const ccomplex_t IT_0702 = m_N_3*IT_0701;
    const ccomplex_t IT_0703 = 2*IT_0702;
    const ccomplex_t IT_0704 = (-2)*IT_0575;
    const ccomplex_t IT_0705 = m_N_3*IT_0472;
    const ccomplex_t IT_0706 = (-4)*IT_0705;
    const ccomplex_t IT_0707 = (-2)*IT_0547;
    const ccomplex_t IT_0708 = IT_0228*IT_0302;
    const ccomplex_t IT_0709 = IT_0232*IT_0637;
    const ccomplex_t IT_0710 = IT_0708 + IT_0709;
    const ccomplex_t IT_0711 = 2*conj(IT_0641);
    const ccomplex_t IT_0712 = (-2)*conj(IT_0575);
    const ccomplex_t IT_0713 = (-2)*conj(IT_0547);
    const ccomplex_t IT_0714 = 2*IT_0369;
    const ccomplex_t IT_0715 = 2*IT_0429;
    const ccomplex_t IT_0716 = m_N_1*IT_0481;
    const ccomplex_t IT_0717 = (-4)*IT_0716;
    const ccomplex_t IT_0718 = 2*IT_0575;
    const ccomplex_t IT_0719 = 2*IT_0547;
    const ccomplex_t IT_0720 = -IT_0691;
    const ccomplex_t IT_0721 = IT_0000 + IT_0720;
    const ccomplex_t IT_0722 = IT_0465*IT_0721;
    const ccomplex_t IT_0723 = (-2)*IT_0641;
    const ccomplex_t IT_0724 = IT_0200*IT_0330 + IT_0325*IT_0442 + IT_0213
      *IT_0617 + IT_0534*IT_0678 + IT_0591*IT_0680 + IT_0220*IT_0693 + IT_0228
      *IT_0714 + IT_0232*IT_0715 + IT_0717*IT_0718 + IT_0689*IT_0719 + IT_0516
      *IT_0722 + IT_0706*IT_0723;
    const ccomplex_t IT_0725 = 2*conj(IT_0369);
    const ccomplex_t IT_0726 = 2*conj(IT_0429);
    const ccomplex_t IT_0727 = 2*conj(IT_0575);
    const ccomplex_t IT_0728 = 2*conj(IT_0547);
    const ccomplex_t IT_0729 = (-2)*conj(IT_0641);
    const ccomplex_t IT_0730 = IT_0200*IT_0302 + IT_0327*IT_0442 + IT_0213
      *IT_0637 + IT_0288*IT_0678 + IT_0187*IT_0680 + IT_0518*IT_0722 + IT_0228
      *IT_0725 + IT_0232*IT_0726 + IT_0717*IT_0727 + IT_0689*IT_0728 + IT_0706
      *IT_0729;
    const ccomplex_t IT_0731 = s_13*m_N_3;
    const ccomplex_t IT_0732 = IT_0522*IT_0731;
    const ccomplex_t IT_0733 = (-4)*conj(IT_0078);
    const ccomplex_t IT_0734 = 0.25*IT_0733;
    const ccomplex_t IT_0735 = conj(IT_0132) + IT_0734;
    const ccomplex_t IT_0736 = s_23*m_N_1;
    const ccomplex_t IT_0737 = IT_0522*IT_0736;
    const ccomplex_t IT_0738 = s_12*s_14*s_34*IT_0043;
    const ccomplex_t IT_0739 = IT_0188 + IT_0189 + IT_0190 + IT_0194 + IT_0199
       + IT_0448 + IT_0450 + IT_0738;
    const ccomplex_t IT_0740 = s_12*s_24*s_34*IT_0043;
    const ccomplex_t IT_0741 = IT_0202 + IT_0203 + IT_0204 + IT_0208 + IT_0212
       + IT_0454 + IT_0456 + IT_0740;
    const ccomplex_t IT_0742 = s_14*IT_0000;
    const ccomplex_t IT_0743 = IT_0597*IT_0742;
    const ccomplex_t IT_0744 = IT_0043*IT_0466;
    const ccomplex_t IT_0745 = 0.25*IT_0743 + (-0.5)*IT_0744;
    const ccomplex_t IT_0746 = s_13 + IT_0745;
    const ccomplex_t IT_0747 = m_N_3*IT_0746;
    const ccomplex_t IT_0748 = 4*IT_0747;
    const ccomplex_t IT_0749 = IT_0043*IT_0475;
    const ccomplex_t IT_0750 = s_24*IT_0000;
    const ccomplex_t IT_0751 = s_23*s_34;
    const ccomplex_t IT_0752 = (-2)*IT_0751;
    const ccomplex_t IT_0753 = IT_0750 + IT_0752;
    const ccomplex_t IT_0754 = IT_0597*IT_0753;
    const ccomplex_t IT_0755 = (-0.5)*IT_0749 + (-0.25)*IT_0754;
    const ccomplex_t IT_0756 = s_23 + IT_0755;
    const ccomplex_t IT_0757 = m_N_1*IT_0756;
    const ccomplex_t IT_0758 = (-4)*IT_0757;
    const ccomplex_t IT_0759 = s_13*s_34;
    const ccomplex_t IT_0760 = (-2)*IT_0759;
    const ccomplex_t IT_0761 = IT_0742 + IT_0760;
    const ccomplex_t IT_0762 = IT_0597*IT_0761;
    const ccomplex_t IT_0763 = (-0.5)*IT_0744 + (-0.25)*IT_0762;
    const ccomplex_t IT_0764 = s_13 + IT_0763;
    const ccomplex_t IT_0765 = m_N_3*IT_0764;
    const ccomplex_t IT_0766 = (-4)*IT_0765;
    const ccomplex_t IT_0767 = s_34*IT_0043;
    const ccomplex_t IT_0768 = IT_0598 + IT_0599;
    const ccomplex_t IT_0769 = s_13*s_23*s_34*IT_0043;
    const ccomplex_t IT_0770 = (-2)*IT_0769;
    const ccomplex_t IT_0771 = IT_0768 + IT_0770;
    const ccomplex_t IT_0772 = IT_0767*IT_0771;
    const ccomplex_t IT_0773 = (-2)*conj(IT_0444);
    const ccomplex_t IT_0774 = 0.25*IT_0710 + IT_0732*IT_0735 + IT_0529
      *IT_0737 + 0.25*IT_0725*IT_0739 + 0.25*IT_0726*IT_0741 + 0.25*IT_0728
      *IT_0748 + 0.25*IT_0711*IT_0758 + 0.25*IT_0712*IT_0766 + 0.25*IT_0772
      *IT_0773;
    const ccomplex_t IT_0775 = IT_0465*IT_0522;
    const ccomplex_t IT_0776 = (-4)*IT_0775;
    const ccomplex_t IT_0777 = -IT_0078 + IT_0132;
    const ccomplex_t IT_0778 = -IT_0147;
    const ccomplex_t IT_0779 = IT_0185 + IT_0778;
    const ccomplex_t IT_0780 = -conj(IT_0147);
    const ccomplex_t IT_0781 = conj(IT_0185) + IT_0780;
    const ccomplex_t IT_0782 = IT_0598 + IT_0600;
    const ccomplex_t IT_0783 = -IT_0599;
    const ccomplex_t IT_0784 = IT_0782 + IT_0783;
    const ccomplex_t IT_0785 = IT_0597*IT_0784;
    const ccomplex_t IT_0786 = 0.25*IT_0785;
    const ccomplex_t IT_0787 = s_12 + IT_0596 + IT_0786;
    const ccomplex_t IT_0788 = 8*IT_0787;
    const ccomplex_t IT_0789 = (-2)*IT_0557 + (-2)*IT_0558 + 4*IT_0560;
    const ccomplex_t IT_0790 = IT_0556 + IT_0789;
    const ccomplex_t IT_0791 = IT_0041*IT_0790;
    const ccomplex_t IT_0792 = 0.25*IT_0791;
    const ccomplex_t IT_0793 = IT_0619 + IT_0624 + IT_0792;
    const ccomplex_t IT_0794 = 1.125*IT_0465;
    const ccomplex_t IT_0795 = IT_0547*IT_0794;
    const ccomplex_t IT_0796 = (-16)*IT_0795;
    const ccomplex_t IT_0797 = (-16)*IT_0641;
    const ccomplex_t IT_0798 = IT_0411 + IT_0428;
    const ccomplex_t IT_0799 = m_N_3*IT_0798;
    const ccomplex_t IT_0800 = m_N_1*IT_0043;
    const ccomplex_t IT_0801 = s_23*IT_0555;
    const ccomplex_t IT_0802 = IT_0043*IT_0801;
    const ccomplex_t IT_0803 = (-2)*IT_0476 + (-0.5)*IT_0802;
    const ccomplex_t IT_0804 = IT_0475 + IT_0803;
    const ccomplex_t IT_0805 = IT_0800*IT_0804;
    const ccomplex_t IT_0806 = (-2)*IT_0805;
    const ccomplex_t IT_0807 = s_13*IT_0555;
    const ccomplex_t IT_0808 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0809 = (-2)*IT_0808;
    const ccomplex_t IT_0810 = IT_0807 + IT_0809;
    const ccomplex_t IT_0811 = IT_0041*IT_0810;
    const ccomplex_t IT_0812 = (-0.666666666666667)*IT_0432 + (
      -0.166666666666667)*IT_0811;
    const ccomplex_t IT_0813 = s_13 + IT_0812;
    const ccomplex_t IT_0814 = m_N_3*IT_0813;
    const ccomplex_t IT_0815 = (-6)*IT_0814;
    const ccomplex_t IT_0816 = (-2)*IT_0467;
    const ccomplex_t IT_0817 = IT_0466 + IT_0816;
    const ccomplex_t IT_0818 = IT_0043*IT_0817;
    const ccomplex_t IT_0819 = 0.166666666666667*IT_0811 + (-0.666666666666667
      )*IT_0818;
    const ccomplex_t IT_0820 = s_13 + IT_0819;
    const ccomplex_t IT_0821 = m_N_3*IT_0820;
    const ccomplex_t IT_0822 = 6*IT_0821;
    const ccomplex_t IT_0823 = IT_0041*IT_0801;
    const ccomplex_t IT_0824 = (-10)*s_23;
    const ccomplex_t IT_0825 = IT_0823 + IT_0824;
    const ccomplex_t IT_0826 = m_N_1*IT_0825;
    const ccomplex_t IT_0827 = IT_0237 + IT_0617;
    const ccomplex_t IT_0828 = (-2)*IT_0217;
    const ccomplex_t IT_0829 = IT_0594 + IT_0828;
    const ccomplex_t IT_0830 = -IT_0296 + (-2)*IT_0433 + 2*IT_0642;
    const ccomplex_t IT_0831 = IT_0295 + IT_0830;
    const ccomplex_t IT_0832 = IT_0043*IT_0831;
    const ccomplex_t IT_0833 = -IT_0045 + (-0.5)*IT_0832;
    const ccomplex_t IT_0834 = s_12 + IT_0833;
    const ccomplex_t IT_0835 = m_N_1*IT_0834;
    const ccomplex_t IT_0836 = (-2)*IT_0835;
    const ccomplex_t IT_0837 = -IT_0296;
    const ccomplex_t IT_0838 = IT_0295 + IT_0837;
    const ccomplex_t IT_0839 = IT_0043*IT_0838;
    const ccomplex_t IT_0840 = -IT_0045 + (-0.25)*IT_0839;
    const ccomplex_t IT_0841 = s_12 + IT_0840;
    const ccomplex_t IT_0842 = m_N_3*IT_0841;
    const ccomplex_t IT_0843 = (-4)*IT_0842;
    const ccomplex_t IT_0844 = m_N_1*IT_0349;
    const ccomplex_t IT_0845 = IT_0011*IT_0346;
    const ccomplex_t IT_0846 = IT_0002*IT_0541*IT_0845;
    const ccomplex_t IT_0847 = m_N_1*IT_0368;
    const ccomplex_t IT_0848 = IT_0028*IT_0365;
    const ccomplex_t IT_0849 = IT_0037*IT_0544*IT_0848;
    const ccomplex_t IT_0850 = IT_0844 + -IT_0846 + IT_0847 + -IT_0849;
    const ccomplex_t IT_0851 = 8*IT_0850;
    const ccomplex_t IT_0852 = (-8)*IT_0575;
    const ccomplex_t IT_0853 = (-8)*IT_0547;
    const ccomplex_t IT_0854 = 0.5*IT_0279 + (-0.5)*IT_0577 + -IT_0609;
    const ccomplex_t IT_0855 = IT_0201 + IT_0854;
    const ccomplex_t IT_0856 = IT_0043*IT_0855;
    const ccomplex_t IT_0857 = -IT_0280;
    const ccomplex_t IT_0858 = IT_0075 + IT_0856 + IT_0857;
    const ccomplex_t IT_0859 = m_N_1*IT_0858;
    const ccomplex_t IT_0860 = (-4)*IT_0324;
    const ccomplex_t IT_0861 = (-2)*IT_0219;
    const ccomplex_t IT_0862 = (-2)*IT_0460;
    const ccomplex_t IT_0863 = -IT_0657;
    const ccomplex_t IT_0864 = IT_0137 + IT_0863;
    const ccomplex_t IT_0865 = IT_0043*IT_0864;
    const ccomplex_t IT_0866 = (-4)*IT_0497;
    const ccomplex_t IT_0867 = IT_0865 + IT_0866;
    const ccomplex_t IT_0868 = m_N_3*IT_0867;
    const ccomplex_t IT_0869 = (-2)*IT_0429;
    const ccomplex_t IT_0870 = (-2)*IT_0236;
    const ccomplex_t IT_0871 = IT_0247*IT_0401;
    const ccomplex_t IT_0872 = IT_0256*IT_0541*IT_0871;
    const ccomplex_t IT_0873 = IT_0265*IT_0425;
    const ccomplex_t IT_0874 = IT_0274*IT_0544*IT_0873;
    const ccomplex_t IT_0875 = -IT_0872 + -IT_0874;
    const ccomplex_t IT_0876 = 0.5*IT_0694*IT_0706 + 0.5*IT_0664*IT_0714 + (
      -0.5)*IT_0776*IT_0779 + (-4)*IT_0777*IT_0787 + IT_0490*IT_0806 + IT_0484
      *IT_0815 + IT_0486*IT_0822 + IT_0492*IT_0826 + 0.5*IT_0584*IT_0827 + 0.5
      *IT_0554*IT_0829 + 0.5*IT_0325*IT_0836 + 0.5*IT_0372*IT_0843 + 0.5*IT_0655
      *IT_0851 + 0.5*IT_0568*IT_0852 + 0.5*IT_0793*IT_0853 + 0.5*IT_0859*IT_0860
       + 0.5*IT_0689*IT_0861 + 0.5*IT_0758*IT_0862 + 0.5*IT_0670*IT_0868 + 0.5
      *IT_0649*IT_0869 + 0.5*IT_0748*IT_0870 + (-0.5)*IT_0794*(IT_0797 + 16
      *IT_0875);
    const ccomplex_t IT_0877 = IT_0239 + IT_0637;
    const ccomplex_t IT_0878 = (-2)*conj(IT_0217);
    const ccomplex_t IT_0879 = IT_0638 + IT_0878;
    const ccomplex_t IT_0880 = 2*conj(IT_0444);
    const ccomplex_t IT_0881 = 8*conj(IT_0850);
    const ccomplex_t IT_0882 = (-8)*conj(IT_0575);
    const ccomplex_t IT_0883 = (-8)*conj(IT_0547);
    const ccomplex_t IT_0884 = (-4)*conj(IT_0324);
    const ccomplex_t IT_0885 = (-2)*conj(IT_0219);
    const ccomplex_t IT_0886 = (-2)*conj(IT_0460);
    const ccomplex_t IT_0887 = (-2)*conj(IT_0429);
    const ccomplex_t IT_0888 = (-2)*conj(IT_0236);
    const ccomplex_t IT_0889 = (-2)*conj(IT_0492);
    const ccomplex_t IT_0890 = (-2)*conj(IT_0486);
    const ccomplex_t IT_0891 = -IT_0664*IT_0725 + IT_0776*IT_0781 + (-8)*(conj
      (IT_0078) + -conj(IT_0132))*IT_0787 + (-2)*conj(IT_0490)*IT_0806 + (-2)
      *conj(IT_0484)*IT_0815 + -IT_0327*IT_0836 + -IT_0374*IT_0843 + -IT_0675
      *IT_0868 + IT_0794*((-16)*conj(IT_0641) + 16*conj(IT_0875)) + -IT_0584
      *IT_0877 + -IT_0554*IT_0879 + -IT_0706*IT_0880 + -IT_0655*IT_0881 + 
      -IT_0568*IT_0882 + (-8)*IT_0793*(conj(IT_0799) + 0.125*IT_0883) + -IT_0859
      *IT_0884 + -IT_0689*IT_0885 + -IT_0758*IT_0886 + -IT_0649*IT_0887 + 
      -IT_0748*IT_0888 + IT_0826*IT_0889 + IT_0822*IT_0890;
    const ccomplex_t IT_0892 = IT_0136 + IT_0496;
    const ccomplex_t IT_0893 = IT_0465*IT_0892;
    const ccomplex_t IT_0894 = s_12 + IT_0567;
    const ccomplex_t IT_0895 = IT_0465*IT_0894;
    const ccomplex_t IT_0896 = IT_0893 + -IT_0895;
    const ccomplex_t IT_0897 = IT_0465*IT_0581;
    const ccomplex_t IT_0898 = IT_0895 + -IT_0897;
    const ccomplex_t IT_0899 = 0.166666666666667*IT_0048;
    const ccomplex_t IT_0900 = 0.166666666666667*IT_0285;
    const ccomplex_t IT_0901 = 0.166666666666667*IT_0143;
    const ccomplex_t IT_0902 = 0.166666666666667*IT_0291;
    const ccomplex_t IT_0903 = 0.166666666666667*IT_0593;
    const ccomplex_t IT_0904 = 4*IT_0217;
    const ccomplex_t IT_0905 = 0.166666666666667*IT_0301;
    const ccomplex_t IT_0906 = 0.166666666666667*IT_0633;
    const ccomplex_t IT_0907 = 0.166666666666667*IT_0554;
    const ccomplex_t IT_0908 = 0.166666666666667*IT_0664;
    const ccomplex_t IT_0909 = IT_0540 + IT_0828;
    const ccomplex_t IT_0910 = 0.166666666666667*IT_0539;
    const ccomplex_t IT_0911 = 0.166666666666667*IT_0616;
    const ccomplex_t IT_0912 = 0.166666666666667*IT_0584;
    const ccomplex_t IT_0913 = 0.166666666666667*IT_0649;
    const ccomplex_t IT_0914 = IT_0142*IT_0465;
    const ccomplex_t IT_0915 = IT_0284*IT_0465;
    const ccomplex_t IT_0916 = IT_0914 + IT_0915;
    const ccomplex_t IT_0917 = 0.166666666666667*IT_0916;
    const ccomplex_t IT_0918 = 0.166666666666667*IT_0233;
    const ccomplex_t IT_0919 = IT_0369*IT_0896 + -IT_0429*IT_0898 + IT_0588
      *IT_0899 + IT_0534*IT_0900 + IT_0591*IT_0901 + IT_0537*IT_0902 + IT_0829
      *IT_0903 + (IT_0330 + IT_0617 + IT_0904)*IT_0905 + IT_0704*IT_0906 +
       IT_0719*IT_0907 + IT_0723*IT_0908 + IT_0909*IT_0910 + IT_0718*IT_0911 +
       IT_0707*IT_0912 + IT_0690*IT_0913 + IT_0372*IT_0917 + IT_0694*IT_0918;
    const ccomplex_t IT_0920 = -conj(IT_0429);
    const ccomplex_t IT_0921 = conj(IT_0369)*IT_0896 + IT_0134*IT_0899 +
       IT_0288*IT_0900 + IT_0187*IT_0901 + IT_0294*IT_0902 + IT_0638*IT_0903 + 
      (IT_0302 + IT_0637)*IT_0905 + IT_0712*IT_0906 + IT_0728*IT_0907 + IT_0729
      *IT_0908 + IT_0585*IT_0910 + IT_0727*IT_0911 + IT_0713*IT_0912 + IT_0711
      *IT_0913 + IT_0374*IT_0917 + IT_0880*IT_0918 + IT_0898*IT_0920;
    const ccomplex_t IT_0922 = (-2)*IT_0369;
    const ccomplex_t IT_0923 = IT_0135 + IT_0863;
    const ccomplex_t IT_0924 = IT_0043*IT_0923;
    const ccomplex_t IT_0925 = (-4)*IT_0892;
    const ccomplex_t IT_0926 = IT_0924 + IT_0925;
    const ccomplex_t IT_0927 = m_N_3*IT_0926;
    const ccomplex_t IT_0928 = 0.166666666666667*IT_0927;
    const ccomplex_t IT_0929 = 0.166666666666667*IT_0200;
    const ccomplex_t IT_0930 = 6*conj(IT_0369)*(IT_0040*IT_0893 + -IT_0277
      *IT_0895 + IT_0591*IT_0899 + IT_0588*IT_0901 + IT_0690*IT_0907 + IT_0707
      *IT_0908 + IT_0903*IT_0922 + IT_0718*IT_0928 + IT_0694*IT_0929);
    const ccomplex_t IT_0931 = 0.5*IT_0137 + -IT_0656 + (-0.5)*IT_0657;
    const ccomplex_t IT_0932 = IT_0135 + IT_0931;
    const ccomplex_t IT_0933 = IT_0043*IT_0932;
    const ccomplex_t IT_0934 = -IT_0138;
    const ccomplex_t IT_0935 = IT_0001 + IT_0933 + IT_0934;
    const ccomplex_t IT_0936 = m_N_3*IT_0935;
    const ccomplex_t IT_0937 = (s_12 + -IT_0045)*(s_12 + (-0.5)*IT_0436);
    const ccomplex_t IT_0938 = (-2)*IT_0937;
    const ccomplex_t IT_0939 = s_12 + IT_0566;
    const ccomplex_t IT_0940 = IT_0465*IT_0939;
    const ccomplex_t IT_0941 = (-2)*IT_0324;
    const ccomplex_t IT_0942 = 0.25*IT_0285;
    const ccomplex_t IT_0943 = 0.25*IT_0143;
    const ccomplex_t IT_0944 = m_N_3*IT_0047;
    const ccomplex_t IT_0945 = 0.25*IT_0944;
    const ccomplex_t IT_0946 = m_N_1*IT_0290;
    const ccomplex_t IT_0947 = 0.25*IT_0946;
    const ccomplex_t IT_0948 = 0.25*IT_0914;
    const ccomplex_t IT_0949 = 0.25*IT_0672;
    const ccomplex_t IT_0950 = 0.25*IT_0673;
    const ccomplex_t IT_0951 = 0.25*IT_0915;
    const ccomplex_t IT_0952 = m_N_1*IT_0841;
    const ccomplex_t IT_0953 = (-4)*IT_0952;
    const ccomplex_t IT_0954 = 0.25*IT_0953;
    const ccomplex_t IT_0955 = IT_0279 + IT_0578;
    const ccomplex_t IT_0956 = IT_0043*IT_0955;
    const ccomplex_t IT_0957 = (-4)*IT_0500;
    const ccomplex_t IT_0958 = IT_0956 + IT_0957;
    const ccomplex_t IT_0959 = m_N_1*IT_0958;
    const ccomplex_t IT_0960 = 0.25*IT_0959;
    const ccomplex_t IT_0961 = 0.25*IT_0868;
    const ccomplex_t IT_0962 = 0.25*IT_0843;
    const ccomplex_t IT_0963 = 0.25*IT_0836;
    const ccomplex_t IT_0964 = 0.25*IT_0394;
    const ccomplex_t IT_0965 = IT_0547*IT_0859 + IT_0850*IT_0936 + IT_0324
      *IT_0938 + (-3)*IT_0371*IT_0940 + 0.25*IT_0669*IT_0941 + IT_0534*IT_0942 +
       IT_0591*IT_0943 + IT_0537*IT_0945 + IT_0588*IT_0947 + IT_0714*IT_0948 +
       IT_0922*IT_0949 + IT_0869*IT_0950 + IT_0715*IT_0951 + IT_0718*IT_0954 +
       IT_0704*IT_0960 + IT_0690*IT_0961 + IT_0723*IT_0962 + IT_0707*IT_0963 +
       IT_0617*IT_0964;
    const ccomplex_t IT_0966 = (-2)*conj(IT_0369);
    const ccomplex_t IT_0967 = conj(IT_0547)*IT_0859 + (-3)*conj(IT_0371)
      *IT_0940 + IT_0288*IT_0942 + IT_0187*IT_0943 + IT_0294*IT_0945 + IT_0134
      *IT_0947 + IT_0725*IT_0948 + IT_0887*IT_0950 + IT_0726*IT_0951 + IT_0727
      *IT_0954 + IT_0712*IT_0960 + IT_0711*IT_0961 + IT_0729*IT_0962 + IT_0713
      *IT_0963 + IT_0637*IT_0964 + IT_0949*IT_0966;
    const ccomplex_t IT_0968 = (-2)*IT_0557;
    const ccomplex_t IT_0969 = IT_0556 + IT_0968;
    const ccomplex_t IT_0970 = IT_0043*IT_0969;
    const ccomplex_t IT_0971 = 0.5*IT_0970;
    const ccomplex_t IT_0972 = (-2)*IT_0557 + 2*IT_0559 + 2*IT_0560 + (-2)
      *IT_0561;
    const ccomplex_t IT_0973 = IT_0556 + IT_0972;
    const ccomplex_t IT_0974 = IT_0043*IT_0973;
    const ccomplex_t IT_0975 = (-0.5)*IT_0974;
    const ccomplex_t IT_0976 = s_14*s_23*IT_0000;
    const ccomplex_t IT_0977 = s_13*s_24*IT_0000;
    const ccomplex_t IT_0978 = IT_0976 + IT_0977;
    const ccomplex_t IT_0979 = s_12*s_34*IT_0000;
    const ccomplex_t IT_0980 = -IT_0979;
    const ccomplex_t IT_0981 = IT_0978 + IT_0980;
    const ccomplex_t IT_0982 = IT_0597*IT_0981;
    const ccomplex_t IT_0983 = 0.5*IT_0982;
    const ccomplex_t IT_0984 = IT_0295 + IT_0434 + IT_0971 + IT_0975 + IT_0983;
    const ccomplex_t IT_0985 = -IT_0190;
    const ccomplex_t IT_0986 = IT_0191 + IT_0985;
    const ccomplex_t IT_0987 = s_12*s_13*IT_0000;
    const ccomplex_t IT_0988 = s_13*s_14*s_24;
    const ccomplex_t IT_0989 = s_23*IT_0000*IT_0001;
    const ccomplex_t IT_0990 = s_24*s_34*IT_0001;
    const ccomplex_t IT_0991 = (-2)*IT_0988 + -IT_0989 + 2*IT_0990;
    const ccomplex_t IT_0992 = IT_0987 + IT_0991;
    const ccomplex_t IT_0993 = IT_0043*IT_0992;
    const ccomplex_t IT_0994 = (-0.5)*IT_0993;
    const ccomplex_t IT_0995 = IT_0986 + IT_0994;
    const ccomplex_t IT_0996 = pow(m_W, 6);
    const ccomplex_t IT_0997 = s_12*IT_0996;
    const ccomplex_t IT_0998 = s_14*s_24*IT_0555;
    const ccomplex_t IT_0999 = (-2)*IT_0998;
    const ccomplex_t IT_1000 = IT_0997 + IT_0999;
    const ccomplex_t IT_1001 = IT_0041*IT_1000;
    const ccomplex_t IT_1002 = 0.5*IT_1001;
    const ccomplex_t IT_1003 = (-4)*IT_0433;
    const ccomplex_t IT_1004 = IT_0295 + IT_1002 + IT_1003;
    const ccomplex_t IT_1005 = IT_0043*IT_0621;
    const ccomplex_t IT_1006 = -IT_1005;
    const ccomplex_t IT_1007 = IT_1004 + IT_1006;
    const ccomplex_t IT_1008 = 0.25*IT_1001;
    const ccomplex_t IT_1009 = IT_0295 + IT_1008;
    const ccomplex_t IT_1010 = -IT_0204;
    const ccomplex_t IT_1011 = IT_0205 + IT_1010;
    const ccomplex_t IT_1012 = s_12*s_23*IT_0000;
    const ccomplex_t IT_1013 = s_13*IT_0279;
    const ccomplex_t IT_1014 = s_12*s_24*s_34;
    const ccomplex_t IT_1015 = s_13*IT_0000*IT_0075;
    const ccomplex_t IT_1016 = 2*IT_1013 + (-2)*IT_1014 + -IT_1015;
    const ccomplex_t IT_1017 = IT_1012 + IT_1016;
    const ccomplex_t IT_1018 = IT_0043*IT_1017;
    const ccomplex_t IT_1019 = (-0.5)*IT_1018;
    const ccomplex_t IT_1020 = IT_1011 + IT_1019;
    const ccomplex_t IT_1021 = IT_0597*IT_0750;
    const ccomplex_t IT_1022 = (-2)*IT_0476;
    const ccomplex_t IT_1023 = IT_0475 + IT_1022;
    const ccomplex_t IT_1024 = IT_0043*IT_1023;
    const ccomplex_t IT_1025 = 0.5*IT_0749 + (-0.5)*IT_1021 + (-0.5)*IT_1024;
    const ccomplex_t IT_1026 = s_23 + IT_1025;
    const ccomplex_t IT_1027 = m_N_1*IT_1026;
    const ccomplex_t IT_1028 = 2*IT_1027;
    const ccomplex_t IT_1029 = 0.25*IT_1028;
    const ccomplex_t IT_1030 = (-0.5)*IT_0743 + 0.5*IT_0744 + (-0.5)*IT_0818;
    const ccomplex_t IT_1031 = s_13 + IT_1030;
    const ccomplex_t IT_1032 = m_N_3*IT_1031;
    const ccomplex_t IT_1033 = 2*IT_1032;
    const ccomplex_t IT_1034 = 0.25*IT_1033;
    const ccomplex_t IT_1035 = (-2)*conj(IT_0850);
    const ccomplex_t IT_1036 = s_24*s_34*IT_0000;
    const ccomplex_t IT_1037 = (-2)*IT_1036;
    const ccomplex_t IT_1038 = IT_0801 + IT_1037;
    const ccomplex_t IT_1039 = IT_0041*IT_1038;
    const ccomplex_t IT_1040 = (-0.666666666666667)*IT_1024 +
       0.166666666666667*IT_1039;
    const ccomplex_t IT_1041 = s_23 + IT_1040;
    const ccomplex_t IT_1042 = m_N_1*IT_1041;
    const ccomplex_t IT_1043 = 6*IT_1042;
    const ccomplex_t IT_1044 = 0.25*IT_1043;
    const ccomplex_t IT_1045 = IT_0041*IT_0807;
    const ccomplex_t IT_1046 = (-10)*s_13;
    const ccomplex_t IT_1047 = IT_1045 + IT_1046;
    const ccomplex_t IT_1048 = m_N_3*IT_1047;
    const ccomplex_t IT_1049 = 0.25*IT_1048;
    const ccomplex_t IT_1050 = 2*conj(IT_0875);
    const ccomplex_t IT_1051 = IT_0729 + IT_1050;
    const ccomplex_t IT_1052 = 0.25*IT_0806;
    const ccomplex_t IT_1053 = 0.25*IT_0815;
    const ccomplex_t IT_1054 = -IT_0197;
    const ccomplex_t IT_1055 = IT_0555 + IT_1054;
    const ccomplex_t IT_1056 = IT_0043*IT_1055;
    const ccomplex_t IT_1057 = -IT_0000 + 0.25*IT_1056;
    const ccomplex_t IT_1058 = IT_0691 + IT_1057;
    const ccomplex_t IT_1059 = IT_0465*IT_1058;
    const ccomplex_t IT_1060 = 4*IT_1059;
    const ccomplex_t IT_1061 = 0.25*IT_1060;
    const ccomplex_t IT_1062 = 2*conj(IT_0460);
    const ccomplex_t IT_1063 = s_34*IT_0000;
    const ccomplex_t IT_1064 = IT_0597*IT_1063;
    const ccomplex_t IT_1065 = IT_0043*IT_0651;
    const ccomplex_t IT_1066 = (-0.333333333333333)*IT_1064 + (
      -0.333333333333333)*IT_1065;
    const ccomplex_t IT_1067 = IT_0000 + IT_1066;
    const ccomplex_t IT_1068 = IT_0465*IT_1067;
    const ccomplex_t IT_1069 = 3*IT_1068;
    const ccomplex_t IT_1070 = 0.25*IT_1069;
    const ccomplex_t IT_1071 = (-2)*conj(IT_0490);
    const ccomplex_t IT_1072 = IT_0000*IT_0465;
    const ccomplex_t IT_1073 = 9*IT_1072;
    const ccomplex_t IT_1074 = 0.25*IT_1073;
    const ccomplex_t IT_1075 = m_N_1*m_N_3*IT_0043;
    const ccomplex_t IT_1076 = (-0.25)*IT_0555;
    const ccomplex_t IT_1077 = IT_0197 + IT_1076;
    const ccomplex_t IT_1078 = IT_1075*IT_1077;
    const ccomplex_t IT_1079 = (-4)*IT_1078;
    const ccomplex_t IT_1080 = 0.25*IT_1079;
    const ccomplex_t IT_1081 = (-4)*conj(IT_0444);
    const ccomplex_t IT_1082 = s_12*IT_0043*IT_0197;
    const ccomplex_t IT_1083 = -IT_1082;
    const ccomplex_t IT_1084 = IT_0556 + IT_0559 + IT_0560;
    const ccomplex_t IT_1085 = (-2)*IT_0557 + -IT_0561;
    const ccomplex_t IT_1086 = IT_1084 + IT_1085;
    const ccomplex_t IT_1087 = IT_0043*IT_1086;
    const ccomplex_t IT_1088 = (-0.5)*IT_1087;
    const ccomplex_t IT_1089 = -IT_0626;
    const ccomplex_t IT_1090 = IT_0295 + IT_0642 + IT_1083 + IT_1088 + IT_1089;
    const ccomplex_t IT_1091 = 0.25*IT_1090;
    const ccomplex_t IT_1092 = (-8)*conj(IT_0371);
    const ccomplex_t IT_1093 = s_23*IT_0137;
    const ccomplex_t IT_1094 = IT_0988 + IT_1093;
    const ccomplex_t IT_1095 = s_12*s_14*s_34;
    const ccomplex_t IT_1096 = -IT_0989 + -IT_1095;
    const ccomplex_t IT_1097 = IT_1094 + IT_1096;
    const ccomplex_t IT_1098 = IT_0043*IT_1097;
    const ccomplex_t IT_1099 = (-0.25)*IT_1098;
    const ccomplex_t IT_1100 = (-0.5)*IT_0738;
    const ccomplex_t IT_1101 = (-0.5)*IT_0447;
    const ccomplex_t IT_1102 = 0.5*IT_0449;
    const ccomplex_t IT_1103 = (-0.5)*IT_0190;
    const ccomplex_t IT_1104 = IT_0191 + IT_1099 + IT_1100 + IT_1101 + IT_1102
       + IT_1103;
    const ccomplex_t IT_1105 = 0.25*IT_1104;
    const ccomplex_t IT_1106 = (-4)*conj(IT_0371);
    const ccomplex_t IT_1107 = s_14*s_23*s_24;
    const ccomplex_t IT_1108 = IT_1013 + IT_1107;
    const ccomplex_t IT_1109 = -IT_1014 + -IT_1015;
    const ccomplex_t IT_1110 = IT_1108 + IT_1109;
    const ccomplex_t IT_1111 = IT_0043*IT_1110;
    const ccomplex_t IT_1112 = (-0.5)*IT_1111;
    const ccomplex_t IT_1113 = -IT_0740;
    const ccomplex_t IT_1114 = IT_0204 + IT_0454 + IT_0455 + IT_1112 + IT_1113;
    const ccomplex_t IT_1115 = 0.25*IT_1114;
    const ccomplex_t IT_1116 = conj(IT_0236)*IT_0984 + conj(IT_0369)*IT_0995 +
       -conj(IT_0484)*IT_1007 + IT_0890*IT_1009 + IT_0920*IT_1020 + IT_0294
      *IT_1029 + IT_0134*IT_1034 + IT_1035*IT_1044 + IT_0727*IT_1049 + IT_1051
      *IT_1052 + IT_0713*IT_1053 + IT_0222*IT_1061 + IT_1062*IT_1070 + IT_1071
      *IT_1074 + IT_0889*IT_1080 + IT_1081*IT_1091 + IT_1092*IT_1105 + IT_1106
      *IT_1115;
    const ccomplex_t IT_1117 = 0.25*IT_0530;
    const ccomplex_t IT_1118 = IT_0078 + IT_1117;
    const ccomplex_t IT_1119 = m_N_1*m_N_3*IT_0000;
    const ccomplex_t IT_1120 = IT_0522*IT_1119;
    const ccomplex_t IT_1121 = (-4)*IT_0236;
    const ccomplex_t IT_1122 = 0.25*IT_0228;
    const ccomplex_t IT_1123 = 0.25*IT_0232;
    const ccomplex_t IT_1124 = (-0.5)*IT_0749 + 0.25*IT_1021;
    const ccomplex_t IT_1125 = s_23 + IT_1124;
    const ccomplex_t IT_1126 = m_N_1*IT_1125;
    const ccomplex_t IT_1127 = 4*IT_1126;
    const ccomplex_t IT_1128 = 0.25*IT_1127;
    const ccomplex_t IT_1129 = 0.25*IT_0748;
    const ccomplex_t IT_1130 = 0.25*IT_0758;
    const ccomplex_t IT_1131 = 0.25*IT_0772;
    const ccomplex_t IT_1132 = IT_0237 + IT_0330;
    const ccomplex_t IT_1133 = 0.25*IT_0739;
    const ccomplex_t IT_1134 = 0.25*IT_0741;
    const ccomplex_t IT_1135 = IT_0732*IT_0779 + IT_0737*IT_1118 + 0.25
      *IT_1120*IT_1121 + IT_0714*IT_1122 + IT_0715*IT_1123 + IT_0704*IT_1128 +
       IT_0690*IT_1129 + IT_0719*IT_1130 + IT_0861*IT_1131 + IT_1132*IT_1133 +
       IT_0827*IT_1134;
    const ccomplex_t IT_1136 = (IT_0295 + (-2)*IT_0433)*(1 + 0.5*IT_0520);
    const ccomplex_t IT_1137 = -conj(IT_0460);
    const ccomplex_t IT_1138 = IT_0239 + IT_0302;
    const ccomplex_t IT_1139 = -IT_0526*IT_0737 + -IT_0732*IT_0781 + conj
      (IT_0236)*IT_1120 + -IT_0725*IT_1122 + -IT_0726*IT_1123 + -IT_0712*IT_1128
       + -IT_0711*IT_1129 + -IT_0728*IT_1130 + -IT_0885*IT_1131 + -IT_0877
      *IT_1134 + -IT_1136*IT_1137 + -IT_1133*IT_1138;
    const ccomplex_t IT_1140 = (-4)*IT_0078;
    const ccomplex_t IT_1141 = 0.25*IT_1140;
    const ccomplex_t IT_1142 = (-2)*IT_0444;
    const ccomplex_t IT_1143 = (IT_0147 + 0.25*IT_0531)*IT_0737 + 0.25*IT_0704
      *IT_0766 + IT_0330*IT_1122 + IT_0617*IT_1123 + IT_0719*IT_1129 + IT_0690
      *IT_1130 + IT_0714*IT_1133 + IT_0715*IT_1134 + 0.25*IT_1121*IT_1136 +
       IT_0732*(IT_0132 + IT_1141) + IT_1131*IT_1142;
    const ccomplex_t IT_1144 = (-0.25)*IT_0143;
    const ccomplex_t IT_1145 = (-0.25)*IT_0944;
    const ccomplex_t IT_1146 = (-0.25)*IT_0946;
    const ccomplex_t IT_1147 = (-0.25)*IT_0914;
    const ccomplex_t IT_1148 = (-0.25)*IT_0672;
    const ccomplex_t IT_1149 = (-0.25)*IT_0673;
    const ccomplex_t IT_1150 = (-0.25)*IT_0915;
    const ccomplex_t IT_1151 = -IT_0296 + (-2)*IT_0433 + 2*IT_0626;
    const ccomplex_t IT_1152 = IT_0295 + IT_1151;
    const ccomplex_t IT_1153 = IT_0043*IT_1152;
    const ccomplex_t IT_1154 = -IT_0045 + (-0.5)*IT_1153;
    const ccomplex_t IT_1155 = s_12 + IT_1154;
    const ccomplex_t IT_1156 = m_N_3*IT_1155;
    const ccomplex_t IT_1157 = (-2)*IT_1156;
    const ccomplex_t IT_1158 = (-0.25)*IT_1157;
    const ccomplex_t IT_1159 = (-0.25)*IT_0868;
    const ccomplex_t IT_1160 = (-0.25)*IT_0843;
    const ccomplex_t IT_1161 = (-0.25)*IT_0836;
    const ccomplex_t IT_1162 = (-0.25)*IT_0285*IT_0537 + -IT_0641*IT_0859 +
       IT_0575*IT_0936 + -IT_0371*IT_0938 + IT_0588*IT_1144 + IT_0534*IT_1145 +
       IT_0591*IT_1146 + IT_0330*IT_1147 + IT_0594*IT_1148 + IT_0540*IT_1149 +
       IT_0617*IT_1150 + IT_0718*IT_1158 + IT_0719*IT_1159 + IT_0707*IT_1160 +
       IT_0723*IT_1161;
    const ccomplex_t IT_1163 = (-0.25)*IT_0285*IT_0294 + -conj(IT_0641)
      *IT_0859 + conj(IT_0575)*IT_0936 + IT_0134*IT_1144 + IT_0288*IT_1145 +
       IT_0187*IT_1146 + IT_0302*IT_1147 + IT_0638*IT_1148 + IT_0585*IT_1149 +
       IT_0637*IT_1150 + IT_0727*IT_1158 + IT_0728*IT_1159 + IT_0713*IT_1160 +
       IT_0729*IT_1161;
    const ccomplex_t IT_1164 = IT_0043*IT_0790;
    const ccomplex_t IT_1165 = -IT_1164;
    const ccomplex_t IT_1166 = IT_1004 + IT_1165;
    const ccomplex_t IT_1167 = -IT_0989 + 2*IT_1093 + (-2)*IT_1095;
    const ccomplex_t IT_1168 = IT_0987 + IT_1167;
    const ccomplex_t IT_1169 = IT_0043*IT_1168;
    const ccomplex_t IT_1170 = (-0.5)*IT_1169;
    const ccomplex_t IT_1171 = IT_0986 + IT_1170;
    const ccomplex_t IT_1172 = s_14*s_34*IT_0075;
    const ccomplex_t IT_1173 = -IT_1015 + (-2)*IT_1107 + 2*IT_1172;
    const ccomplex_t IT_1174 = IT_1012 + IT_1173;
    const ccomplex_t IT_1175 = IT_0043*IT_1174;
    const ccomplex_t IT_1176 = (-0.5)*IT_1175;
    const ccomplex_t IT_1177 = IT_1011 + IT_1176;
    const ccomplex_t IT_1178 = (-0.25)*IT_1028;
    const ccomplex_t IT_1179 = (-0.25)*IT_1033;
    const ccomplex_t IT_1180 = (-0.666666666666667)*IT_0439 + (
      -0.166666666666667)*IT_1039;
    const ccomplex_t IT_1181 = s_23 + IT_1180;
    const ccomplex_t IT_1182 = m_N_1*IT_1181;
    const ccomplex_t IT_1183 = (-6)*IT_1182;
    const ccomplex_t IT_1184 = (-0.25)*IT_1183;
    const ccomplex_t IT_1185 = m_N_3*IT_0043;
    const ccomplex_t IT_1186 = IT_0043*IT_0807;
    const ccomplex_t IT_1187 = (-2)*IT_0467 + (-0.5)*IT_1186;
    const ccomplex_t IT_1188 = IT_0466 + IT_1187;
    const ccomplex_t IT_1189 = IT_1185*IT_1188;
    const ccomplex_t IT_1190 = (-2)*IT_1189;
    const ccomplex_t IT_1191 = (-0.25)*IT_1190;
    const ccomplex_t IT_1192 = (-0.25)*IT_0822;
    const ccomplex_t IT_1193 = (-0.25)*IT_1060;
    const ccomplex_t IT_1194 = (-0.25)*IT_1069;
    const ccomplex_t IT_1195 = -IT_0642;
    const ccomplex_t IT_1196 = IT_0295 + IT_0626 + IT_1083 + IT_1088 + IT_1195;
    const ccomplex_t IT_1197 = (-0.25)*IT_1196;
    const ccomplex_t IT_1198 = (-0.5)*IT_1098;
    const ccomplex_t IT_1199 = -IT_0738;
    const ccomplex_t IT_1200 = IT_0190 + IT_0447 + IT_0450 + IT_1198 + IT_1199;
    const ccomplex_t IT_1201 = (-0.25)*IT_1200;
    const ccomplex_t IT_1202 = (-0.25)*IT_1111;
    const ccomplex_t IT_1203 = (-0.5)*IT_0740;
    const ccomplex_t IT_1204 = 0.5*IT_0453;
    const ccomplex_t IT_1205 = (-0.5)*IT_0455;
    const ccomplex_t IT_1206 = (-0.5)*IT_0204;
    const ccomplex_t IT_1207 = IT_0205 + IT_1202 + IT_1203 + IT_1204 + IT_1205
       + IT_1206;
    const ccomplex_t IT_1208 = (-0.25)*IT_1207;
    const ccomplex_t IT_1209 = -conj(IT_0236)*IT_0984 + 2*conj(IT_0484)
      *IT_1009 + (-0.25)*IT_0826*IT_1051 + (-0.25)*IT_0889*IT_1073 + (-0.25)
      *IT_1071*IT_1079 + conj(IT_0486)*IT_1166 + conj(IT_0369)*IT_1171 + IT_0920
      *IT_1177 + IT_0294*IT_1178 + IT_0134*IT_1179 + IT_1035*IT_1184 + IT_0727
      *IT_1191 + IT_0713*IT_1192 + IT_0222*IT_1193 + IT_1062*IT_1194 + IT_1081
      *IT_1197 + IT_1106*IT_1201 + IT_1092*IT_1208;
    const ccomplex_t IT_1210 = 2*conj(IT_0236);
    const ccomplex_t IT_1211 = (-4)*conj(IT_0219);
    const ccomplex_t IT_1212 = (-8)*conj(IT_0324);
    const ccomplex_t IT_1213 = 4*conj(IT_0217);
    const ccomplex_t IT_1214 = 0.25*IT_1213;
    const ccomplex_t IT_1215 = -IT_1020;
    const ccomplex_t IT_1216 = (-4)*conj(IT_0217);
    const ccomplex_t IT_1217 = (-0.25)*IT_1216;
    const ccomplex_t IT_1218 = conj(IT_0460)*IT_0984 + -conj(IT_0490)*IT_1007 
      + IT_0889*IT_1009 + IT_0288*IT_1029 + IT_0187*IT_1034 + IT_0727*IT_1044 +
       IT_1035*IT_1049 + IT_0713*IT_1052 + IT_1051*IT_1053 + IT_0880*IT_1061 +
       IT_0884*IT_1115 + IT_1070*IT_1210 + IT_1091*IT_1211 + IT_1105*IT_1212 +
       IT_0995*(conj(IT_0040) + IT_1214) + IT_1215*(conj(IT_0277) + IT_1217);
    const ccomplex_t IT_1219 = 0.666666666666667*IT_0749 + (-0.166666666666667
      )*IT_0823;
    const ccomplex_t IT_1220 = s_23 + IT_1219;
    const ccomplex_t IT_1221 = m_N_1*IT_1220;
    const ccomplex_t IT_1222 = (-6)*IT_1221;
    const ccomplex_t IT_1223 = 2*conj(IT_0490)*IT_1009 + IT_0984*IT_1137 +
       conj(IT_0492)*IT_1166 + IT_0288*IT_1178 + IT_0187*IT_1179 + IT_0727
      *IT_1184 + IT_1035*IT_1191 + IT_1051*IT_1192 + IT_0880*IT_1193 + IT_0884
      *IT_1201 + IT_1194*IT_1210 + IT_1197*IT_1211 + IT_1208*IT_1212 + -IT_1177*
      (conj(IT_0277) + IT_1214) + IT_1171*(conj(IT_0040) + IT_1217) + (-0.25)
      *IT_0713*IT_1222;
    const ccomplex_t IT_1224 = (-0.166666666666667)*IT_0291;
    const ccomplex_t IT_1225 = (-0.166666666666667)*IT_0301;
    const ccomplex_t IT_1226 = (-0.166666666666667)*IT_0584;
    const ccomplex_t IT_1227 = (-0.166666666666667)*IT_0649;
    const ccomplex_t IT_1228 = m_N_3*IT_0552;
    const ccomplex_t IT_1229 = (-4)*IT_1228;
    const ccomplex_t IT_1230 = (-0.166666666666667)*IT_1229;
    const ccomplex_t IT_1231 = (-0.166666666666667)*IT_0213;
    const ccomplex_t IT_1232 = (-0.166666666666667)*IT_0285*IT_0537 + (
      -0.166666666666667)*IT_0539*IT_0869 + IT_0040*IT_0895 + -IT_0277*IT_0897 +
       IT_0534*IT_1224 + IT_0714*IT_1225 + IT_0723*IT_1226 + IT_0719*IT_1227 +
       IT_0704*IT_1230 + IT_0694*IT_1231;
    const ccomplex_t IT_1233 = (-6)*IT_0429*((-0.166666666666667)*IT_0285
      *IT_0294 + conj(IT_0040)*IT_0895 + -conj(IT_0277)*IT_0897 + IT_0288
      *IT_1224 + IT_0725*IT_1225 + IT_0729*IT_1226 + IT_0728*IT_1227 + IT_0712
      *IT_1230 + IT_0880*IT_1231);
    const ccomplex_t IT_1234 = IT_0534*IT_1029;
    const ccomplex_t IT_1235 = IT_0591*IT_1034;
    const ccomplex_t IT_1236 = (-2)*IT_0850;
    const ccomplex_t IT_1237 = 2*IT_0875;
    const ccomplex_t IT_1238 = IT_0723 + IT_1237;
    const ccomplex_t IT_1239 = IT_0694*IT_1061;
    const ccomplex_t IT_1240 = 2*IT_0236;
    const ccomplex_t IT_1241 = IT_1070*IT_1240;
    const ccomplex_t IT_1242 = (-4)*IT_0219;
    const ccomplex_t IT_1243 = (-8)*IT_0324;
    const ccomplex_t IT_1244 = 0.25*IT_0904;
    const ccomplex_t IT_1245 = (-4)*IT_0217;
    const ccomplex_t IT_1246 = (-0.25)*IT_1245;
    const ccomplex_t IT_1247 = IT_0460*IT_0984;
    const ccomplex_t IT_1248 = IT_0718*IT_1044 + IT_0707*IT_1052 + IT_0860
      *IT_1115 + IT_1234 + IT_1235 + IT_1049*IT_1236 + IT_1053*IT_1238 + IT_1239
       + IT_1241 + IT_1091*IT_1242 + IT_1105*IT_1243 + IT_0995*(IT_0040 +
       IT_1244) + IT_1215*(IT_0277 + IT_1246) + IT_1247;
    const ccomplex_t IT_1249 = 0.25*IT_1183;
    const ccomplex_t IT_1250 = 0.25*IT_1190;
    const ccomplex_t IT_1251 = 0.25*IT_0822;
    const ccomplex_t IT_1252 = 0.25*IT_1196;
    const ccomplex_t IT_1253 = 0.25*IT_1200;
    const ccomplex_t IT_1254 = 0.25*IT_1207;
    const ccomplex_t IT_1255 = IT_0707*IT_1222 + 4*IT_1234 + 4*IT_1235 + 4
      *IT_1239 + 4*IT_1241 + 4*IT_1177*(IT_0277 + IT_1244) + (-4)*IT_1171*
      (IT_0040 + IT_1246) + 4*IT_1247 + 4*IT_0718*IT_1249 + 4*IT_1236*IT_1250 +
       4*IT_1238*IT_1251 + 4*IT_1242*IT_1252 + 4*IT_0860*IT_1253 + 4*IT_1243
      *IT_1254;
    const ccomplex_t IT_1256 = (-8)*IT_0641;
    const ccomplex_t IT_1257 = (-4)*IT_0371;
    const ccomplex_t IT_1258 = 8*IT_0147;
    const ccomplex_t IT_1259 = (-0.125)*IT_1258;
    const ccomplex_t IT_1260 = IT_0185 + IT_1259;
    const ccomplex_t IT_1261 = IT_0776*IT_0777;
    const ccomplex_t IT_1262 = IT_0220*IT_0706 + IT_0584*IT_0715 + -IT_0796 +
       IT_0372*IT_0836 + IT_0325*IT_0843 + IT_0568*IT_0851 + IT_0655*IT_0852 +
       IT_0748*IT_0862 + IT_0758*IT_0870 + IT_0649*IT_0909 + IT_0554*IT_0922 +
       IT_0868*IT_0941 + IT_0664*IT_1132 + IT_0689*IT_1142 + IT_0793*IT_1256 +
       IT_0859*IT_1257 + -IT_0788*IT_1260 + -IT_1261;
    const ccomplex_t IT_1263 = 2*IT_0460;
    const ccomplex_t IT_1264 = 0.25*IT_0220*IT_0703 + 0.25*IT_0633*IT_0714 + 4
      *IT_0575*IT_0794 + 0.25*(IT_0619 + IT_0623 + IT_0624)*IT_0851 + 0.25
      *IT_0655*IT_0853 + 0.25*IT_0616*IT_0869 + 0.25*IT_0829*IT_0927 + IT_0324
      *IT_0936 + 0.25*IT_0670*IT_0953 + 0.25*IT_0372*IT_0959 + 0.25*IT_0717
      *IT_1142 + 0.25*IT_0941*IT_1157 + 0.25*IT_0827*IT_1229 + 0.25*IT_1127
      *IT_1240 + 0.25*IT_0568*IT_1256 + (-2)*IT_0605*IT_1260 + (-0.25)*IT_1261 +
       0.25*IT_0766*IT_1263;
    const ccomplex_t IT_1265 = IT_0537*IT_1029;
    const ccomplex_t IT_1266 = (-2)*IT_0490;
    const ccomplex_t IT_1267 = (-2)*IT_0492;
    const ccomplex_t IT_1268 = (-4)*IT_0444;
    const ccomplex_t IT_1269 = (-8)*IT_0371;
    const ccomplex_t IT_1270 = IT_0236*IT_0984;
    const ccomplex_t IT_1271 = IT_0369*IT_0995 + -IT_0429*IT_1020 + IT_0588
      *IT_1034 + IT_0718*IT_1049 + IT_0707*IT_1053 + IT_0220*IT_1061 + IT_1044
      *IT_1236 + IT_1052*IT_1238 + IT_1115*IT_1257 + IT_1070*IT_1263 + IT_1265 +
       IT_1074*IT_1266 + IT_1080*IT_1267 + IT_1091*IT_1268 + IT_1105*IT_1269 +
       IT_1270;
    const ccomplex_t IT_1272 = -IT_0588*IT_1034 + -IT_0220*IT_1061 + IT_0369
      *IT_1171 + -IT_0429*IT_1177 + (-0.25)*IT_0826*IT_1238 + -IT_1236*IT_1249 +
       -IT_0718*IT_1250 + -IT_0707*IT_1251 + -IT_1253*IT_1257 + -IT_1070*IT_1263
       + -IT_1265 + -IT_1080*IT_1266 + -IT_1074*IT_1267 + -IT_1252*IT_1268 + 
      -IT_1254*IT_1269 + -IT_1270;
    const ccomplex_t IT_1273 = 8*conj(IT_0147);
    const ccomplex_t IT_1274 = (-8)*conj(IT_0185);
    const ccomplex_t IT_1275 = IT_1273 + IT_1274;
    const ccomplex_t IT_1276 = (-8)*conj(IT_0641);
    const ccomplex_t IT_1277 = (-2)*conj(IT_0324);
    const ccomplex_t IT_1278 = -conj(IT_0132);
    const ccomplex_t IT_1279 = conj(IT_0078) + IT_1278;
    const ccomplex_t IT_1280 = IT_0776*IT_1279;
    const ccomplex_t IT_1281 = IT_0222*IT_0703 + IT_0633*IT_0725 + IT_0717
      *IT_0773 + 16*conj(IT_0575)*IT_0794 + IT_0655*IT_0883 + IT_0616*IT_0887 +
       IT_0879*IT_0927 + IT_0675*IT_0953 + IT_0374*IT_0959 + IT_0766*IT_1062 +
       IT_1127*IT_1210 + IT_0877*IT_1229 + IT_0605*IT_1275 + IT_0568*IT_1276 +
       IT_1157*IT_1277 + IT_1280;
    const ccomplex_t IT_1282 = IT_0222*IT_0706 + IT_0584*IT_0726 + IT_0689
      *IT_0773 + 16*conj(IT_0547)*IT_0794 + IT_0374*IT_0836 + IT_0327*IT_0843 +
       IT_0649*(IT_0585 + IT_0878) + IT_0568*IT_0881 + IT_0655*IT_0882 + IT_0748
      *IT_0886 + IT_0758*IT_0888 + IT_0554*IT_0966 + IT_0859*IT_1106 + IT_0664
      *IT_1138 + IT_0787*IT_1275 + 8*IT_0793*(conj(IT_0875) + 0.125*IT_1276) +
       IT_0868*IT_1277 + IT_1280;
    const ccomplex_t IT_1283 = IT_0040*(IT_0048*IT_0134 + IT_0143*IT_0187) + 
      (IT_0200 + IT_0213)*(IT_0221 + IT_0223) + IT_0233*(IT_0238 + IT_0240) +
       IT_0277*(IT_0285*IT_0288 + IT_0291*IT_0294 + IT_0301*IT_0302) + (IT_0312 
      + IT_0322)*(IT_0329 + IT_0331 + IT_0332 + IT_0373 + IT_0375) + IT_0394*
      (IT_0329 + IT_0430 + IT_0431) + IT_0442*(IT_0445 + IT_0446) + (IT_0452 +
       IT_0458)*(IT_0461 + IT_0462 + IT_0463 + IT_0464) + (IT_0474 + IT_0483)*
      (IT_0485 + IT_0487 + IT_0488 + IT_0489 + IT_0491 + IT_0493 + IT_0494 +
       IT_0495) + (IT_0498 + IT_0501)*(IT_0503 + IT_0505) + (IT_0508 + IT_0511)*
      (IT_0513 + IT_0515 + IT_0517 + IT_0519) + 4*IT_0523*(IT_0132*conj(IT_0132)
       + IT_0185*conj(IT_0185) + IT_0078*IT_0526 + IT_0147*IT_0529 + 0.25*conj
      (IT_0078)*IT_0530 + 0.25*conj(IT_0147)*IT_0531) + conj(IT_0277)*(IT_0301
      *IT_0330 + IT_0285*IT_0534 + IT_0291*IT_0537 + IT_0539*IT_0540) + IT_0547*
      (IT_0302*IT_0554 + IT_0568*IT_0576 + IT_0584*IT_0585) + conj(IT_0040)*
      (IT_0048*IT_0588 + IT_0143*IT_0591 + IT_0593*IT_0594) + 8*conj(IT_0575)*
      (IT_0605*IT_0608 + 0.125*IT_0616*IT_0617 + 0.125*(IT_0619 + IT_0623 +
       IT_0624)*IT_0625 + 0.125*IT_0594*IT_0633) + IT_0575*(IT_0605*IT_0636 +
       IT_0616*IT_0637 + IT_0633*IT_0638) + IT_0641*(IT_0637*IT_0649 + IT_0576
      *IT_0655 + IT_0638*IT_0664) + IT_0669*IT_0671 + (IT_0672 + IT_0673)*
      (IT_0674 + IT_0676) + conj(IT_0444)*(IT_0228*IT_0330 + IT_0232*IT_0617 +
       IT_0537*IT_0678 + IT_0588*IT_0680 + IT_0689*IT_0690 + IT_0693*IT_0694 +
       IT_0703*IT_0704 + IT_0706*IT_0707) + IT_0444*(IT_0294*IT_0678 + IT_0134
      *IT_0680 + IT_0710 + IT_0689*IT_0711 + IT_0703*IT_0712 + IT_0706*IT_0713) 
      + conj(IT_0219)*IT_0724 + IT_0219*IT_0730 + 4*IT_0236*IT_0774 + -IT_0776*(
      (conj(IT_0078) + -conj(IT_0132))*IT_0641 + -conj(IT_0641)*IT_0777 + -(conj
      (IT_0547) + conj(IT_0575))*IT_0779 + -(IT_0547 + IT_0575)*(conj(IT_0185) +
       IT_0780)) + (conj(IT_0547)*IT_0608 + 0.125*IT_0547*IT_0636 + conj(IT_0641
      )*IT_0779 + IT_0641*IT_0781)*IT_0788 + 8*conj(IT_0641)*(0.125*IT_0617
      *IT_0649 + 0.125*IT_0625*IT_0655 + 0.125*IT_0594*IT_0664 + IT_0641*IT_0793
       + 0.125*IT_0796) + 8*conj(IT_0547)*(0.125*IT_0330*IT_0554 + 0.125*IT_0540
      *IT_0584 + 0.125*IT_0568*IT_0625 + IT_0547*IT_0793 + 0.125*IT_0794*IT_0797
      ) + 2*conj(IT_0799)*IT_0876 + -IT_0799*IT_0891 + 6*conj(IT_0217)*IT_0919 +
       6*IT_0217*IT_0921 + 6*IT_0369*(conj(IT_0040)*IT_0893 + -conj(IT_0277)
      *IT_0895 + IT_0187*IT_0899 + IT_0134*IT_0901 + IT_0711*IT_0907 + IT_0713
      *IT_0908 + IT_0727*IT_0928 + IT_0880*IT_0929) + IT_0930 + 4*conj(IT_0324)
      *IT_0965 + 4*IT_0324*IT_0967 + 4*IT_0484*IT_1116 + 4*conj(IT_0460)*IT_1135
       + (-4)*IT_0460*IT_1139 + 4*IT_0775*(IT_0185*IT_0526 + IT_0147*IT_0735 +
       conj(IT_0185)*IT_1118 + conj(IT_0147)*(IT_0132 + IT_1141)) + 4*conj
      (IT_0236)*IT_1143 + (-4)*conj(IT_0371)*IT_1162 + (-4)*IT_0371*IT_1163 + (
      -4)*IT_0486*IT_1209 + 4*IT_0490*IT_1218 + (-4)*IT_0492*IT_1223 + (-6)*conj
      (IT_0429)*IT_1232 + IT_1233 + 4*conj(IT_0490)*IT_1248 + conj(IT_0492)
      *IT_1255 + conj(IT_0875)*IT_1262 + 4*conj(IT_0850)*IT_1264 + 4*conj
      (IT_0484)*IT_1271 + (-4)*conj(IT_0486)*IT_1272 + IT_0850*IT_1281 + IT_0875
      *IT_1282;
    return create_ccomplex_return(IT_1283);
}

