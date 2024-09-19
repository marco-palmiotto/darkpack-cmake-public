#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_4_N_4_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_4_N_4_to_anti_t_t(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_h = param->m_h;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t m_A0 = param->m_A0;
    const creal_t m_H0 = param->m_H0;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_24 = param->s_24;
    const creal_t s_34 = param->s_34;
    const creal_t alpha = param->alpha;
    const creal_t m_N_4 = param->m_N_4;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t Gamma_Z = param->Gamma_Z;
    const creal_t Gamma_h = param->Gamma_h;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_A0 = param->Gamma_A0;
    const creal_t Gamma_H0 = param->Gamma_H0;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t N_B4 = param->N_B4;
    const ccomplex_t N_W4 = param->N_W4;
    const ccomplex_t N_d4 = param->N_d4;
    const ccomplex_t N_u4 = param->N_u4;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = m_t*m_N_4;
    const ccomplex_t IT_0001 = -s_14;
    const ccomplex_t IT_0002 = s_13 + IT_0001;
    const ccomplex_t IT_0003 = IT_0000*IT_0002;
    const ccomplex_t IT_0004 = 4*IT_0003;
    const ccomplex_t IT_0005 = cos(theta_W);
    const ccomplex_t IT_0006 = cpow(IT_0005, -1);
    const ccomplex_t IT_0007 = sin(theta_W);
    const ccomplex_t IT_0008 = IT_0006*IT_0007;
    const ccomplex_t IT_0009 = N_u4*conj(N_u4)*e_em;
    const ccomplex_t IT_0010 = IT_0008*IT_0009;
    const ccomplex_t IT_0011 = cpow(IT_0007, -1);
    const ccomplex_t IT_0012 = IT_0005*IT_0011;
    const ccomplex_t IT_0013 = IT_0009*IT_0012;
    const ccomplex_t IT_0014 = N_d4*conj(N_d4)*e_em;
    const ccomplex_t IT_0015 = IT_0008*IT_0014;
    const ccomplex_t IT_0016 = IT_0012*IT_0014;
    const ccomplex_t IT_0017 = (0 + _Complex_I*1)*(IT_0010 + IT_0013 + 
      -IT_0015 + -IT_0016);
    const ccomplex_t IT_0018 = 0.5*IT_0017;
    const ccomplex_t IT_0019 = e_em*IT_0008;
    const ccomplex_t IT_0020 = e_em*IT_0012;
    const ccomplex_t IT_0021 = (0 + _Complex_I*1)*(IT_0019 + (-3)*IT_0020);
    const ccomplex_t IT_0022 = (-0.166666666666667)*IT_0021;
    const ccomplex_t IT_0023 = pow(m_N_4, 2);
    const ccomplex_t IT_0024 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0025 = IT_0022*IT_0024;
    const ccomplex_t IT_0026 = (0 + _Complex_I*1)*IT_0025;
    const ccomplex_t IT_0027 = IT_0018*IT_0026;
    const ccomplex_t IT_0028 = 0.5*IT_0027;
    const ccomplex_t IT_0029 = (0 + _Complex_I*1)*e_em*IT_0006*IT_0007;
    const ccomplex_t IT_0030 = (-0.666666666666667)*IT_0029;
    const ccomplex_t IT_0031 = IT_0024*IT_0030;
    const ccomplex_t IT_0032 = (0 + _Complex_I*1)*IT_0031;
    const ccomplex_t IT_0033 = IT_0018*IT_0032;
    const ccomplex_t IT_0034 = 0.5*IT_0033;
    const ccomplex_t IT_0035 = cos(alpha);
    const ccomplex_t IT_0036 = IT_0006*IT_0035;
    const ccomplex_t IT_0037 = N_B4*N_d4*e_em;
    const ccomplex_t IT_0038 = IT_0036*IT_0037;
    const ccomplex_t IT_0039 = sin(alpha);
    const ccomplex_t IT_0040 = IT_0006*IT_0039;
    const ccomplex_t IT_0041 = N_B4*N_u4*e_em;
    const ccomplex_t IT_0042 = IT_0040*IT_0041;
    const ccomplex_t IT_0043 = IT_0011*IT_0035;
    const ccomplex_t IT_0044 = N_W4*N_d4*e_em;
    const ccomplex_t IT_0045 = IT_0043*IT_0044;
    const ccomplex_t IT_0046 = IT_0011*IT_0039;
    const ccomplex_t IT_0047 = N_W4*N_u4*e_em;
    const ccomplex_t IT_0048 = IT_0046*IT_0047;
    const ccomplex_t IT_0049 = (0 + _Complex_I*1)*(IT_0038 + -IT_0042 + 
      -IT_0045 + IT_0048);
    const ccomplex_t IT_0050 = pow(m_W, -1);
    const ccomplex_t IT_0051 = sin(beta);
    const ccomplex_t IT_0052 = cpow(IT_0051, -1);
    const ccomplex_t IT_0053 = (0 + _Complex_I*1)*m_t*e_em*IT_0011*IT_0039
      *IT_0050*IT_0052;
    const ccomplex_t IT_0054 = (-0.5)*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_12 + (-2)*IT_0023 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = IT_0049*IT_0057;
    const ccomplex_t IT_0059 = IT_0036*IT_0041;
    const ccomplex_t IT_0060 = IT_0037*IT_0040;
    const ccomplex_t IT_0061 = IT_0043*IT_0047;
    const ccomplex_t IT_0062 = IT_0044*IT_0046;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + IT_0060 + 
      -IT_0061 + -IT_0062);
    const ccomplex_t IT_0064 = -IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*m_t*e_em*IT_0011*IT_0035
      *IT_0050*IT_0052;
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0023 + -reg_prop, -1);
    const ccomplex_t IT_0068 = IT_0066*IT_0067;
    const ccomplex_t IT_0069 = (0 + _Complex_I*1)*IT_0068;
    const ccomplex_t IT_0070 = IT_0064*IT_0069;
    const ccomplex_t IT_0071 = (-0.5)*IT_0058 + (-0.5)*IT_0070;
    const ccomplex_t IT_0072 = 3*IT_0071;
    const ccomplex_t IT_0073 = 0.5*IT_0058 + 0.5*IT_0070;
    const ccomplex_t IT_0074 = (-3)*IT_0073;
    const ccomplex_t IT_0075 = IT_0072 + IT_0074;
    const ccomplex_t IT_0076 = 3*conj(IT_0071);
    const ccomplex_t IT_0077 = (-3)*conj(IT_0073);
    const ccomplex_t IT_0078 = IT_0076 + IT_0077;
    const ccomplex_t IT_0079 = (-0.5)*IT_0017;
    const ccomplex_t IT_0080 = IT_0026*IT_0079;
    const ccomplex_t IT_0081 = 0.5*IT_0080;
    const ccomplex_t IT_0082 = IT_0032*IT_0079;
    const ccomplex_t IT_0083 = 0.5*IT_0082;
    const ccomplex_t IT_0084 = conj(N_B4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0085 = IT_0036*IT_0084;
    const ccomplex_t IT_0086 = conj(N_B4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0087 = IT_0040*IT_0086;
    const ccomplex_t IT_0088 = conj(N_W4)*conj(N_d4)*e_em;
    const ccomplex_t IT_0089 = IT_0043*IT_0088;
    const ccomplex_t IT_0090 = conj(N_W4)*conj(N_u4)*e_em;
    const ccomplex_t IT_0091 = IT_0046*IT_0090;
    const ccomplex_t IT_0092 = (0 + _Complex_I*1)*(IT_0085 + -IT_0087 + 
      -IT_0089 + IT_0091);
    const ccomplex_t IT_0093 = IT_0057*IT_0092;
    const ccomplex_t IT_0094 = IT_0036*IT_0086;
    const ccomplex_t IT_0095 = IT_0040*IT_0084;
    const ccomplex_t IT_0096 = IT_0043*IT_0090;
    const ccomplex_t IT_0097 = IT_0046*IT_0088;
    const ccomplex_t IT_0098 = (0 + _Complex_I*1)*(IT_0094 + IT_0095 + 
      -IT_0096 + -IT_0097);
    const ccomplex_t IT_0099 = -IT_0098;
    const ccomplex_t IT_0100 = IT_0069*IT_0099;
    const ccomplex_t IT_0101 = (-0.5)*IT_0093 + (-0.5)*IT_0100;
    const ccomplex_t IT_0102 = 3*IT_0101;
    const ccomplex_t IT_0103 = 0.5*IT_0093 + 0.5*IT_0100;
    const ccomplex_t IT_0104 = (-3)*IT_0103;
    const ccomplex_t IT_0105 = IT_0102 + IT_0104;
    const ccomplex_t IT_0106 = 3*conj(IT_0101);
    const ccomplex_t IT_0107 = (-3)*conj(IT_0103);
    const ccomplex_t IT_0108 = IT_0106 + IT_0107;
    const ccomplex_t IT_0109 = -s_24;
    const ccomplex_t IT_0110 = s_23 + IT_0109;
    const ccomplex_t IT_0111 = IT_0000*IT_0110;
    const ccomplex_t IT_0112 = (-4)*IT_0111;
    const ccomplex_t IT_0113 = 2*IT_0003;
    const ccomplex_t IT_0114 = N_B4*e_em*conj(U_st_01);
    const ccomplex_t IT_0115 = IT_0006*IT_0114;
    const ccomplex_t IT_0116 = 1.4142135623731*IT_0115;
    const ccomplex_t IT_0117 = N_W4*e_em*conj(U_st_01);
    const ccomplex_t IT_0118 = IT_0011*IT_0117;
    const ccomplex_t IT_0119 = 1.4142135623731*IT_0118;
    const ccomplex_t IT_0120 = IT_0011*IT_0052;
    const ccomplex_t IT_0121 = m_t*N_u4*e_em*IT_0050*conj(U_st_11);
    const ccomplex_t IT_0122 = IT_0120*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = (0 + _Complex_I*1)*(IT_0116 + 3*IT_0119 + 3
      *IT_0123);
    const ccomplex_t IT_0125 = 0.166666666666667*IT_0124;
    const ccomplex_t IT_0126 = conj(N_B4)*e_em*U_st_01;
    const ccomplex_t IT_0127 = IT_0006*IT_0126;
    const ccomplex_t IT_0128 = 1.4142135623731*IT_0127;
    const ccomplex_t IT_0129 = conj(N_W4)*e_em*U_st_01;
    const ccomplex_t IT_0130 = IT_0011*IT_0129;
    const ccomplex_t IT_0131 = 1.4142135623731*IT_0130;
    const ccomplex_t IT_0132 = m_t*conj(N_u4)*e_em*IT_0050*U_st_11;
    const ccomplex_t IT_0133 = IT_0120*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*(IT_0128 + 3*IT_0131 + 3
      *IT_0134);
    const ccomplex_t IT_0136 = 0.166666666666667*IT_0135;
    const ccomplex_t IT_0137 = IT_0125*IT_0136;
    const ccomplex_t IT_0138 = pow(m_t, 2);
    const ccomplex_t IT_0139 = cpow((-2)*s_13 + IT_0023 + IT_0138 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0140 = IT_0137*IT_0139;
    const ccomplex_t IT_0141 = (0 + _Complex_I*1)*IT_0140;
    const ccomplex_t IT_0142 = cpow((-2)*s_13 + IT_0023 + IT_0138 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0143 = N_B4*e_em*conj(U_st_00);
    const ccomplex_t IT_0144 = IT_0006*IT_0143;
    const ccomplex_t IT_0145 = 1.4142135623731*IT_0144;
    const ccomplex_t IT_0146 = N_W4*e_em*conj(U_st_00);
    const ccomplex_t IT_0147 = IT_0011*IT_0146;
    const ccomplex_t IT_0148 = 1.4142135623731*IT_0147;
    const ccomplex_t IT_0149 = m_t*N_u4*e_em*IT_0050*conj(U_st_10);
    const ccomplex_t IT_0150 = IT_0120*IT_0149;
    const ccomplex_t IT_0151 = 1.4142135623731*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*(IT_0145 + 3*IT_0148 + 3
      *IT_0151);
    const ccomplex_t IT_0153 = 0.166666666666667*IT_0152;
    const ccomplex_t IT_0154 = conj(N_B4)*e_em*U_st_00;
    const ccomplex_t IT_0155 = IT_0006*IT_0154;
    const ccomplex_t IT_0156 = 1.4142135623731*IT_0155;
    const ccomplex_t IT_0157 = conj(N_W4)*e_em*U_st_00;
    const ccomplex_t IT_0158 = IT_0011*IT_0157;
    const ccomplex_t IT_0159 = 1.4142135623731*IT_0158;
    const ccomplex_t IT_0160 = m_t*conj(N_u4)*e_em*IT_0050*U_st_10;
    const ccomplex_t IT_0161 = IT_0120*IT_0160;
    const ccomplex_t IT_0162 = 1.4142135623731*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*(IT_0156 + 3*IT_0159 + 3
      *IT_0162);
    const ccomplex_t IT_0164 = 0.166666666666667*IT_0163;
    const ccomplex_t IT_0165 = IT_0153*IT_0164;
    const ccomplex_t IT_0166 = IT_0142*IT_0165;
    const ccomplex_t IT_0167 = (0 + _Complex_I*1)*IT_0166;
    const ccomplex_t IT_0168 = -IT_0141 + -IT_0167;
    const ccomplex_t IT_0169 = conj(N_B4)*e_em*conj(U_st_10);
    const ccomplex_t IT_0170 = IT_0006*IT_0169;
    const ccomplex_t IT_0171 = 1.4142135623731*IT_0170;
    const ccomplex_t IT_0172 = m_t*conj(N_u4)*e_em*IT_0050*conj(U_st_00);
    const ccomplex_t IT_0173 = IT_0120*IT_0172;
    const ccomplex_t IT_0174 = 1.4142135623731*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*(IT_0171 + (-0.75)*IT_0174);
    const ccomplex_t IT_0176 = 0.666666666666667*IT_0175;
    const ccomplex_t IT_0177 = N_B4*e_em*U_st_10;
    const ccomplex_t IT_0178 = IT_0006*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = m_t*N_u4*e_em*IT_0050*U_st_00;
    const ccomplex_t IT_0181 = IT_0120*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = (0 + _Complex_I*1)*(IT_0179 + (-0.75)*IT_0182);
    const ccomplex_t IT_0184 = 0.666666666666667*IT_0183;
    const ccomplex_t IT_0185 = IT_0176*IT_0184;
    const ccomplex_t IT_0186 = IT_0142*IT_0185;
    const ccomplex_t IT_0187 = (0 + _Complex_I*1)*IT_0186;
    const ccomplex_t IT_0188 = conj(N_B4)*e_em*conj(U_st_11);
    const ccomplex_t IT_0189 = IT_0006*IT_0188;
    const ccomplex_t IT_0190 = 1.4142135623731*IT_0189;
    const ccomplex_t IT_0191 = m_t*conj(N_u4)*e_em*IT_0050*conj(U_st_01);
    const ccomplex_t IT_0192 = IT_0120*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*(IT_0190 + (-0.75)*IT_0193);
    const ccomplex_t IT_0195 = 0.666666666666667*IT_0194;
    const ccomplex_t IT_0196 = N_B4*e_em*U_st_11;
    const ccomplex_t IT_0197 = IT_0006*IT_0196;
    const ccomplex_t IT_0198 = 1.4142135623731*IT_0197;
    const ccomplex_t IT_0199 = m_t*N_u4*e_em*IT_0050*U_st_01;
    const ccomplex_t IT_0200 = IT_0120*IT_0199;
    const ccomplex_t IT_0201 = 1.4142135623731*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*(IT_0198 + (-0.75)*IT_0201);
    const ccomplex_t IT_0203 = 0.666666666666667*IT_0202;
    const ccomplex_t IT_0204 = IT_0195*IT_0203;
    const ccomplex_t IT_0205 = IT_0139*IT_0204;
    const ccomplex_t IT_0206 = (0 + _Complex_I*1)*IT_0205;
    const ccomplex_t IT_0207 = -IT_0187 + -IT_0206;
    const ccomplex_t IT_0208 = (-2)*IT_0111;
    const ccomplex_t IT_0209 = (-2)*IT_0003;
    const ccomplex_t IT_0210 = cpow((-2)*s_23 + IT_0023 + IT_0138 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0211 = IT_0204*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = cpow((-2)*s_23 + IT_0023 + IT_0138 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0214 = IT_0185*IT_0213;
    const ccomplex_t IT_0215 = (0 + _Complex_I*1)*IT_0214;
    const ccomplex_t IT_0216 = IT_0212 + IT_0215;
    const ccomplex_t IT_0217 = 3*IT_0073;
    const ccomplex_t IT_0218 = (-3)*IT_0071;
    const ccomplex_t IT_0219 = IT_0217 + IT_0218;
    const ccomplex_t IT_0220 = 3*conj(IT_0073);
    const ccomplex_t IT_0221 = (-3)*conj(IT_0071);
    const ccomplex_t IT_0222 = IT_0220 + IT_0221;
    const ccomplex_t IT_0223 = IT_0137*IT_0210;
    const ccomplex_t IT_0224 = (0 + _Complex_I*1)*IT_0223;
    const ccomplex_t IT_0225 = IT_0165*IT_0213;
    const ccomplex_t IT_0226 = (0 + _Complex_I*1)*IT_0225;
    const ccomplex_t IT_0227 = IT_0224 + IT_0226;
    const ccomplex_t IT_0228 = 3*IT_0103;
    const ccomplex_t IT_0229 = (-3)*IT_0101;
    const ccomplex_t IT_0230 = IT_0228 + IT_0229;
    const ccomplex_t IT_0231 = 3*conj(IT_0103);
    const ccomplex_t IT_0232 = (-3)*conj(IT_0101);
    const ccomplex_t IT_0233 = IT_0231 + IT_0232;
    const ccomplex_t IT_0234 = 2*IT_0111;
    const ccomplex_t IT_0235 = 4*IT_0111;
    const ccomplex_t IT_0236 = (-0.5)*IT_0027;
    const ccomplex_t IT_0237 = (-0.5)*IT_0033;
    const ccomplex_t IT_0238 = conj(IT_0236) + conj(IT_0237);
    const ccomplex_t IT_0239 = IT_0236 + IT_0237;
    const ccomplex_t IT_0240 = (-0.5)*IT_0080;
    const ccomplex_t IT_0241 = (-0.5)*IT_0082;
    const ccomplex_t IT_0242 = conj(IT_0240) + conj(IT_0241);
    const ccomplex_t IT_0243 = IT_0240 + IT_0241;
    const ccomplex_t IT_0244 = (-4)*IT_0003;
    const ccomplex_t IT_0245 = s_12*IT_0138;
    const ccomplex_t IT_0246 = s_13*s_24;
    const ccomplex_t IT_0247 = s_12*s_34;
    const ccomplex_t IT_0248 = s_14*s_23;
    const ccomplex_t IT_0249 = -IT_0248;
    const ccomplex_t IT_0250 = IT_0245 + IT_0246 + IT_0247 + IT_0249;
    const ccomplex_t IT_0251 = IT_0125*IT_0203;
    const ccomplex_t IT_0252 = IT_0139*IT_0251;
    const ccomplex_t IT_0253 = (0 + _Complex_I*1)*IT_0252;
    const ccomplex_t IT_0254 = IT_0153*IT_0184;
    const ccomplex_t IT_0255 = IT_0142*IT_0254;
    const ccomplex_t IT_0256 = (0 + _Complex_I*1)*IT_0255;
    const ccomplex_t IT_0257 = IT_0253 + IT_0256;
    const ccomplex_t IT_0258 = cos(beta);
    const ccomplex_t IT_0259 = IT_0006*IT_0258;
    const ccomplex_t IT_0260 = IT_0037*IT_0259;
    const ccomplex_t IT_0261 = IT_0006*IT_0051;
    const ccomplex_t IT_0262 = IT_0041*IT_0261;
    const ccomplex_t IT_0263 = IT_0011*IT_0258;
    const ccomplex_t IT_0264 = IT_0044*IT_0263;
    const ccomplex_t IT_0265 = IT_0011*IT_0051;
    const ccomplex_t IT_0266 = IT_0047*IT_0265;
    const ccomplex_t IT_0267 = -IT_0260 + -IT_0262 + IT_0264 + IT_0266;
    const ccomplex_t IT_0268 = m_t*e_em*IT_0011*IT_0050;
    const ccomplex_t IT_0269 = (-0.5)*IT_0268;
    const ccomplex_t IT_0270 = pow(m_Z, 2);
    const ccomplex_t IT_0271 = cpow((-2)*s_12 + (-2)*IT_0023 + IT_0270 + 
      -reg_prop, -1);
    const ccomplex_t IT_0272 = IT_0269*IT_0271;
    const ccomplex_t IT_0273 = (0 + _Complex_I*1)*IT_0272;
    const ccomplex_t IT_0274 = IT_0267*IT_0273;
    const ccomplex_t IT_0275 = IT_0041*IT_0259;
    const ccomplex_t IT_0276 = IT_0037*IT_0261;
    const ccomplex_t IT_0277 = IT_0047*IT_0263;
    const ccomplex_t IT_0278 = IT_0044*IT_0265;
    const ccomplex_t IT_0279 = -IT_0275 + IT_0276 + IT_0277 + -IT_0278;
    const ccomplex_t IT_0280 = m_t*e_em*IT_0011*IT_0050*IT_0052*IT_0258;
    const ccomplex_t IT_0281 = (-0.5)*IT_0280;
    const ccomplex_t IT_0282 = cpow((-2)*s_12 + (-2)*IT_0023 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0283 = IT_0281*IT_0282;
    const ccomplex_t IT_0284 = (0 + _Complex_I*1)*IT_0283;
    const ccomplex_t IT_0285 = IT_0279*IT_0284;
    const ccomplex_t IT_0286 = 0.5*IT_0274 + 0.5*IT_0285;
    const ccomplex_t IT_0287 = 6*IT_0286;
    const ccomplex_t IT_0288 = (-0.5)*IT_0274 + (-0.5)*IT_0285;
    const ccomplex_t IT_0289 = (-6)*IT_0288;
    const ccomplex_t IT_0290 = IT_0287 + IT_0289;
    const ccomplex_t IT_0291 = 6*conj(IT_0286);
    const ccomplex_t IT_0292 = (-6)*conj(IT_0288);
    const ccomplex_t IT_0293 = IT_0291 + IT_0292;
    const ccomplex_t IT_0294 = IT_0136*IT_0195;
    const ccomplex_t IT_0295 = IT_0139*IT_0294;
    const ccomplex_t IT_0296 = (0 + _Complex_I*1)*IT_0295;
    const ccomplex_t IT_0297 = IT_0164*IT_0176;
    const ccomplex_t IT_0298 = IT_0142*IT_0297;
    const ccomplex_t IT_0299 = (0 + _Complex_I*1)*IT_0298;
    const ccomplex_t IT_0300 = IT_0296 + IT_0299;
    const ccomplex_t IT_0301 = IT_0086*IT_0259;
    const ccomplex_t IT_0302 = IT_0084*IT_0261;
    const ccomplex_t IT_0303 = IT_0090*IT_0263;
    const ccomplex_t IT_0304 = IT_0088*IT_0265;
    const ccomplex_t IT_0305 = IT_0301 + -IT_0302 + -IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = IT_0284*IT_0305;
    const ccomplex_t IT_0307 = IT_0084*IT_0259;
    const ccomplex_t IT_0308 = IT_0086*IT_0261;
    const ccomplex_t IT_0309 = IT_0088*IT_0263;
    const ccomplex_t IT_0310 = IT_0090*IT_0265;
    const ccomplex_t IT_0311 = IT_0307 + IT_0308 + -IT_0309 + -IT_0310;
    const ccomplex_t IT_0312 = IT_0273*IT_0311;
    const ccomplex_t IT_0313 = (-0.5)*IT_0306 + (-0.5)*IT_0312;
    const ccomplex_t IT_0314 = 6*IT_0313;
    const ccomplex_t IT_0315 = 0.5*IT_0306 + 0.5*IT_0312;
    const ccomplex_t IT_0316 = (-6)*IT_0315;
    const ccomplex_t IT_0317 = IT_0314 + IT_0316;
    const ccomplex_t IT_0318 = 6*conj(IT_0313);
    const ccomplex_t IT_0319 = (-6)*conj(IT_0315);
    const ccomplex_t IT_0320 = IT_0318 + IT_0319;
    const ccomplex_t IT_0321 = -IT_0246;
    const ccomplex_t IT_0322 = IT_0245 + IT_0247 + IT_0248 + IT_0321;
    const ccomplex_t IT_0323 = IT_0210*IT_0251;
    const ccomplex_t IT_0324 = (0 + _Complex_I*1)*IT_0323;
    const ccomplex_t IT_0325 = IT_0213*IT_0254;
    const ccomplex_t IT_0326 = (0 + _Complex_I*1)*IT_0325;
    const ccomplex_t IT_0327 = -IT_0324 + -IT_0326;
    const ccomplex_t IT_0328 = 6*IT_0288;
    const ccomplex_t IT_0329 = (-6)*IT_0286;
    const ccomplex_t IT_0330 = IT_0328 + IT_0329;
    const ccomplex_t IT_0331 = 6*conj(IT_0288);
    const ccomplex_t IT_0332 = (-6)*conj(IT_0286);
    const ccomplex_t IT_0333 = IT_0331 + IT_0332;
    const ccomplex_t IT_0334 = IT_0213*IT_0297;
    const ccomplex_t IT_0335 = (0 + _Complex_I*1)*IT_0334;
    const ccomplex_t IT_0336 = IT_0210*IT_0294;
    const ccomplex_t IT_0337 = (0 + _Complex_I*1)*IT_0336;
    const ccomplex_t IT_0338 = -IT_0335 + -IT_0337;
    const ccomplex_t IT_0339 = 6*IT_0315;
    const ccomplex_t IT_0340 = (-6)*IT_0313;
    const ccomplex_t IT_0341 = IT_0339 + IT_0340;
    const ccomplex_t IT_0342 = 6*conj(IT_0315);
    const ccomplex_t IT_0343 = (-6)*conj(IT_0313);
    const ccomplex_t IT_0344 = IT_0342 + IT_0343;
    const ccomplex_t IT_0345 = s_23 + s_24;
    const ccomplex_t IT_0346 = IT_0000*IT_0345;
    const ccomplex_t IT_0347 = conj(IT_0081) + conj(IT_0236);
    const ccomplex_t IT_0348 = 12*IT_0288;
    const ccomplex_t IT_0349 = (-12)*IT_0286;
    const ccomplex_t IT_0350 = IT_0348 + IT_0349;
    const ccomplex_t IT_0351 = IT_0081 + IT_0236;
    const ccomplex_t IT_0352 = 12*conj(IT_0288);
    const ccomplex_t IT_0353 = (-12)*conj(IT_0286);
    const ccomplex_t IT_0354 = IT_0352 + IT_0353;
    const ccomplex_t IT_0355 = conj(IT_0083) + conj(IT_0237);
    const ccomplex_t IT_0356 = 12*IT_0286;
    const ccomplex_t IT_0357 = (-12)*IT_0288;
    const ccomplex_t IT_0358 = IT_0356 + IT_0357;
    const ccomplex_t IT_0359 = IT_0083 + IT_0237;
    const ccomplex_t IT_0360 = 12*conj(IT_0286);
    const ccomplex_t IT_0361 = (-12)*conj(IT_0288);
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = conj(IT_0028) + conj(IT_0240);
    const ccomplex_t IT_0364 = 12*IT_0313;
    const ccomplex_t IT_0365 = (-12)*IT_0315;
    const ccomplex_t IT_0366 = IT_0364 + IT_0365;
    const ccomplex_t IT_0367 = IT_0028 + IT_0240;
    const ccomplex_t IT_0368 = 12*conj(IT_0313);
    const ccomplex_t IT_0369 = (-12)*conj(IT_0315);
    const ccomplex_t IT_0370 = IT_0368 + IT_0369;
    const ccomplex_t IT_0371 = conj(IT_0034) + conj(IT_0241);
    const ccomplex_t IT_0372 = 12*IT_0315;
    const ccomplex_t IT_0373 = (-12)*IT_0313;
    const ccomplex_t IT_0374 = IT_0372 + IT_0373;
    const ccomplex_t IT_0375 = IT_0034 + IT_0241;
    const ccomplex_t IT_0376 = 12*conj(IT_0315);
    const ccomplex_t IT_0377 = (-12)*conj(IT_0313);
    const ccomplex_t IT_0378 = IT_0376 + IT_0377;
    const ccomplex_t IT_0379 = conj(IT_0207)*IT_0290 + IT_0207*IT_0293 + conj
      (IT_0168)*IT_0317 + IT_0168*IT_0320 + conj(IT_0227)*IT_0330 + IT_0227
      *IT_0333 + conj(IT_0216)*IT_0341 + IT_0216*IT_0344 + IT_0347*IT_0350 +
       IT_0351*IT_0354 + IT_0355*IT_0358 + IT_0359*IT_0362 + IT_0363*IT_0366 +
       IT_0367*IT_0370 + IT_0371*IT_0374 + IT_0375*IT_0378;
    const ccomplex_t IT_0380 = -IT_0247;
    const ccomplex_t IT_0381 = IT_0246 + IT_0248 + IT_0380;
    const ccomplex_t IT_0382 = 6*IT_0257;
    const ccomplex_t IT_0383 = 6*conj(IT_0257);
    const ccomplex_t IT_0384 = s_13 + s_14;
    const ccomplex_t IT_0385 = IT_0000*IT_0384;
    const ccomplex_t IT_0386 = 6*IT_0207;
    const ccomplex_t IT_0387 = 6*conj(IT_0207);
    const ccomplex_t IT_0388 = 6*IT_0227;
    const ccomplex_t IT_0389 = 6*conj(IT_0227);
    const ccomplex_t IT_0390 = conj(IT_0168)*IT_0290 + IT_0168*IT_0293 + conj
      (IT_0216)*IT_0330 + IT_0216*IT_0333 + IT_0358*IT_0363 + IT_0355*IT_0366 +
       IT_0362*IT_0367 + IT_0359*IT_0370 + IT_0350*IT_0371 + IT_0347*IT_0374 +
       IT_0354*IT_0375 + IT_0351*IT_0378 + (-6)*conj(IT_0313)*(IT_0227 + (
      -0.166666666666667)*IT_0386) + (-6)*IT_0313*(conj(IT_0227) + (
      -0.166666666666667)*IT_0387) + (-6)*conj(IT_0315)*(IT_0207 + (
      -0.166666666666667)*IT_0388) + (-6)*IT_0315*(conj(IT_0207) + (
      -0.166666666666667)*IT_0389);
    const ccomplex_t IT_0391 = m_t*s_14*m_N_4;
    const ccomplex_t IT_0392 = 12*IT_0227;
    const ccomplex_t IT_0393 = 24*IT_0236;
    const ccomplex_t IT_0394 = 24*IT_0081;
    const ccomplex_t IT_0395 = IT_0392 + IT_0393 + IT_0394;
    const ccomplex_t IT_0396 = 12*conj(IT_0227);
    const ccomplex_t IT_0397 = 24*conj(IT_0081);
    const ccomplex_t IT_0398 = 24*conj(IT_0236);
    const ccomplex_t IT_0399 = IT_0396 + IT_0397 + IT_0398;
    const ccomplex_t IT_0400 = 12*IT_0257;
    const ccomplex_t IT_0401 = 24*IT_0327;
    const ccomplex_t IT_0402 = IT_0400 + IT_0401;
    const ccomplex_t IT_0403 = 12*conj(IT_0257);
    const ccomplex_t IT_0404 = 24*conj(IT_0327);
    const ccomplex_t IT_0405 = IT_0403 + IT_0404;
    const ccomplex_t IT_0406 = 12*IT_0327;
    const ccomplex_t IT_0407 = IT_0382 + IT_0406;
    const ccomplex_t IT_0408 = 12*conj(IT_0327);
    const ccomplex_t IT_0409 = IT_0383 + IT_0408;
    const ccomplex_t IT_0410 = 12*IT_0300;
    const ccomplex_t IT_0411 = 12*conj(IT_0300);
    const ccomplex_t IT_0412 = IT_0391*(conj(IT_0300)*IT_0388 + IT_0300
      *IT_0389 + conj(IT_0338)*IT_0395 + IT_0338*IT_0399 + IT_0371*IT_0402 +
       IT_0375*IT_0405 + conj(IT_0216)*IT_0407 + IT_0216*IT_0409 + IT_0347
      *IT_0410 + IT_0351*IT_0411);
    const ccomplex_t IT_0413 = m_t*s_23*m_N_4;
    const ccomplex_t IT_0414 = 12*IT_0216;
    const ccomplex_t IT_0415 = 24*IT_0241;
    const ccomplex_t IT_0416 = 24*IT_0034;
    const ccomplex_t IT_0417 = IT_0414 + IT_0415 + IT_0416;
    const ccomplex_t IT_0418 = 12*conj(IT_0216);
    const ccomplex_t IT_0419 = 24*conj(IT_0241);
    const ccomplex_t IT_0420 = 24*conj(IT_0034);
    const ccomplex_t IT_0421 = IT_0418 + IT_0419 + IT_0420;
    const ccomplex_t IT_0422 = conj(IT_0216)*IT_0300 + IT_0216*conj(IT_0300) +
       0.166666666666667*IT_0347*IT_0402 + 0.166666666666667*IT_0351*IT_0405 +
       0.166666666666667*conj(IT_0227)*IT_0407 + 0.166666666666667*IT_0227
      *IT_0409 + 0.166666666666667*IT_0371*IT_0410 + 0.166666666666667*IT_0375
      *IT_0411 + 0.166666666666667*conj(IT_0338)*IT_0417 + 0.166666666666667
      *IT_0338*IT_0421;
    const ccomplex_t IT_0423 = m_t*s_24*m_N_4;
    const ccomplex_t IT_0424 = 12*IT_0168;
    const ccomplex_t IT_0425 = 24*IT_0028;
    const ccomplex_t IT_0426 = 24*IT_0240;
    const ccomplex_t IT_0427 = IT_0424 + IT_0425 + IT_0426;
    const ccomplex_t IT_0428 = 12*conj(IT_0168);
    const ccomplex_t IT_0429 = 24*conj(IT_0028);
    const ccomplex_t IT_0430 = 24*conj(IT_0240);
    const ccomplex_t IT_0431 = IT_0428 + IT_0429 + IT_0430;
    const ccomplex_t IT_0432 = 24*IT_0257;
    const ccomplex_t IT_0433 = IT_0406 + IT_0432;
    const ccomplex_t IT_0434 = 24*conj(IT_0257);
    const ccomplex_t IT_0435 = IT_0408 + IT_0434;
    const ccomplex_t IT_0436 = 6*IT_0327;
    const ccomplex_t IT_0437 = IT_0400 + IT_0436;
    const ccomplex_t IT_0438 = 6*conj(IT_0327);
    const ccomplex_t IT_0439 = IT_0403 + IT_0438;
    const ccomplex_t IT_0440 = 6*IT_0168;
    const ccomplex_t IT_0441 = 6*conj(IT_0168);
    const ccomplex_t IT_0442 = 12*IT_0338;
    const ccomplex_t IT_0443 = 12*conj(IT_0338);
    const ccomplex_t IT_0444 = IT_0423*(conj(IT_0300)*IT_0427 + IT_0300
      *IT_0431 + IT_0355*IT_0433 + IT_0359*IT_0435 + conj(IT_0207)*IT_0437 +
       IT_0207*IT_0439 + conj(IT_0338)*IT_0440 + IT_0338*IT_0441 + IT_0363
      *IT_0442 + IT_0367*IT_0443);
    const ccomplex_t IT_0445 = s_34*IT_0023;
    const ccomplex_t IT_0446 = 12*IT_0207;
    const ccomplex_t IT_0447 = 12*conj(IT_0207);
    const ccomplex_t IT_0448 = conj(IT_0216)*IT_0386 + IT_0216*IT_0387 + conj
      (IT_0168)*IT_0388 + IT_0168*IT_0389 + IT_0363*IT_0392 + IT_0028*(IT_0396 +
       IT_0397) + IT_0240*IT_0399 + IT_0355*IT_0417 + IT_0359*IT_0421 + IT_0347
      *IT_0424 + conj(IT_0236)*IT_0425 + (IT_0081 + IT_0236)*IT_0431 + IT_0371
      *IT_0446 + IT_0375*IT_0447;
    const ccomplex_t IT_0449 = -s_34;
    const ccomplex_t IT_0450 = IT_0138 + IT_0449;
    const ccomplex_t IT_0451 = s_12*IT_0450;
    const ccomplex_t IT_0452 = 24*conj(IT_0103);
    const ccomplex_t IT_0453 = (-24)*conj(IT_0101);
    const ccomplex_t IT_0454 = IT_0452 + IT_0453;
    const ccomplex_t IT_0455 = (-24)*IT_0071;
    const ccomplex_t IT_0456 = s_34 + IT_0138;
    const ccomplex_t IT_0457 = s_12*IT_0456;
    const ccomplex_t IT_0458 = (-24)*IT_0288;
    const ccomplex_t IT_0459 = (-24)*conj(IT_0288);
    const ccomplex_t IT_0460 = 6*conj(IT_0073);
    const ccomplex_t IT_0461 = (-6)*conj(IT_0071);
    const ccomplex_t IT_0462 = IT_0460 + IT_0461;
    const ccomplex_t IT_0463 = 6*IT_0103;
    const ccomplex_t IT_0464 = (-6)*IT_0101;
    const ccomplex_t IT_0465 = IT_0463 + IT_0464;
    const ccomplex_t IT_0466 = 6*conj(IT_0103);
    const ccomplex_t IT_0467 = (-6)*conj(IT_0101);
    const ccomplex_t IT_0468 = IT_0466 + IT_0467;
    const ccomplex_t IT_0469 = (-6)*IT_0071;
    const ccomplex_t IT_0470 = IT_0023*IT_0450;
    const ccomplex_t IT_0471 = 6*conj(IT_0071);
    const ccomplex_t IT_0472 = (-6)*conj(IT_0073);
    const ccomplex_t IT_0473 = IT_0471 + IT_0472;
    const ccomplex_t IT_0474 = 6*IT_0101;
    const ccomplex_t IT_0475 = (-6)*IT_0103;
    const ccomplex_t IT_0476 = IT_0474 + IT_0475;
    const ccomplex_t IT_0477 = 6*IT_0071;
    const ccomplex_t IT_0478 = (-6)*conj(IT_0103);
    const ccomplex_t IT_0479 = conj(IT_0073)*((-24)*IT_0101 + 24*IT_0103) +
       conj(IT_0071)*(24*IT_0101 + (-24)*IT_0103) + conj(IT_0101)*(24*IT_0071 +
       6*IT_0257) + (-6)*IT_0073*(conj(IT_0300) + -conj(IT_0338) + (
      -0.166666666666667)*IT_0454) + conj(IT_0103)*IT_0455 + IT_0338*IT_0462 +
       conj(IT_0327)*IT_0465 + IT_0327*IT_0468 + conj(IT_0338)*IT_0469 + IT_0300
      *IT_0473 + conj(IT_0257)*IT_0476 + conj(IT_0300)*IT_0477 + IT_0257*IT_0478;
    const ccomplex_t IT_0480 = 24*conj(IT_0237);
    const ccomplex_t IT_0481 = 24*conj(IT_0083);
    const ccomplex_t IT_0482 = IT_0447 + IT_0480 + IT_0481;
    const ccomplex_t IT_0483 = conj(IT_0227)*IT_0386 + IT_0227*IT_0387 +
       IT_0355*IT_0392 + IT_0359*IT_0399 + IT_0363*IT_0414 + IT_0367*IT_0418 +
       IT_0371*(IT_0424 + IT_0425 + IT_0426) + IT_0375*IT_0431 + conj(IT_0216)
      *IT_0440 + IT_0216*IT_0441 + IT_0347*IT_0446 + IT_0351*IT_0482;
    const ccomplex_t IT_0484 = m_t*s_13*m_N_4;
    const ccomplex_t IT_0485 = 24*IT_0083;
    const ccomplex_t IT_0486 = 24*IT_0237;
    const ccomplex_t IT_0487 = conj(IT_0338)*IT_0386 + IT_0338*IT_0387 +
       IT_0363*IT_0433 + IT_0367*IT_0435 + conj(IT_0168)*IT_0437 + IT_0168
      *IT_0439 + IT_0355*IT_0442 + IT_0359*IT_0443 + IT_0300*IT_0482 + conj
      (IT_0300)*(IT_0446 + IT_0485 + IT_0486);
    const ccomplex_t IT_0488 = 48*IT_0083;
    const ccomplex_t IT_0489 = 24*IT_0207;
    const ccomplex_t IT_0490 = IT_0488 + IT_0489;
    const ccomplex_t IT_0491 = 48*IT_0028;
    const ccomplex_t IT_0492 = 24*IT_0168;
    const ccomplex_t IT_0493 = IT_0491 + IT_0492;
    const ccomplex_t IT_0494 = 48*conj(IT_0028);
    const ccomplex_t IT_0495 = 48*IT_0240;
    const ccomplex_t IT_0496 = 48*conj(IT_0083);
    const ccomplex_t IT_0497 = 24*conj(IT_0207);
    const ccomplex_t IT_0498 = conj(IT_0257)*IT_0400 + conj(IT_0300)*IT_0410 +
       conj(IT_0168)*IT_0427 + conj(IT_0207)*(IT_0446 + IT_0485) + IT_0355
      *IT_0490 + conj(IT_0028)*IT_0493 + IT_0240*IT_0494 + conj(IT_0240)*
      (IT_0493 + IT_0495) + IT_0237*(48*conj(IT_0237) + IT_0496 + IT_0497);
    const ccomplex_t IT_0499 = IT_0023*IT_0138;
    const ccomplex_t IT_0500 = 48*IT_0081;
    const ccomplex_t IT_0501 = 24*IT_0227;
    const ccomplex_t IT_0502 = IT_0500 + IT_0501;
    const ccomplex_t IT_0503 = 48*conj(IT_0034);
    const ccomplex_t IT_0504 = 48*IT_0034;
    const ccomplex_t IT_0505 = 48*IT_0241;
    const ccomplex_t IT_0506 = 48*conj(IT_0081);
    const ccomplex_t IT_0507 = 24*conj(IT_0227);
    const ccomplex_t IT_0508 = conj(IT_0216)*IT_0395 + IT_0216*IT_0399 + conj
      (IT_0338)*IT_0407 + IT_0338*IT_0409 + conj(IT_0300)*IT_0437 + IT_0300
      *IT_0439 + IT_0168*IT_0482 + conj(IT_0168)*(IT_0446 + IT_0485 + IT_0486) +
       conj(IT_0240)*(48*IT_0237 + IT_0488 + IT_0489) + conj(IT_0028)*IT_0490 +
       IT_0237*IT_0494 + conj(IT_0237)*(IT_0491 + IT_0495) + IT_0367*(IT_0496 +
       IT_0497) + conj(IT_0241)*(48*IT_0236 + IT_0500 + IT_0501) + conj(IT_0034)
      *IT_0502 + IT_0236*IT_0503 + conj(IT_0236)*(IT_0504 + IT_0505) + IT_0375*
      (IT_0506 + IT_0507);
    const ccomplex_t IT_0509 = 24*IT_0216;
    const ccomplex_t IT_0510 = IT_0504 + IT_0509;
    const ccomplex_t IT_0511 = conj(IT_0227)*(IT_0392 + IT_0394) + conj
      (IT_0327)*IT_0406 + conj(IT_0216)*IT_0417 + conj(IT_0338)*IT_0442 +
       IT_0347*IT_0502 + IT_0241*IT_0503 + IT_0236*(48*conj(IT_0236) + IT_0506 +
       IT_0507) + conj(IT_0034)*IT_0510 + conj(IT_0241)*(IT_0505 + IT_0510);
    const ccomplex_t IT_0512 = (-0.166666666666667)*conj(IT_0338);
    const ccomplex_t IT_0513 = IT_0023*IT_0456;
    const ccomplex_t IT_0514 = IT_0293*IT_0300 + IT_0290*conj(IT_0300) +
       IT_0333*IT_0338 + (-6)*conj(IT_0315)*(IT_0257 + 4*IT_0286 + (-4)*IT_0288 
      + (-0.166666666666667)*IT_0436) + (-6)*IT_0315*(conj(IT_0257) + 4*conj
      (IT_0286) + (-4)*conj(IT_0288) + (-0.166666666666667)*IT_0438) + 24*conj
      (IT_0313)*(IT_0286 + -1./4*IT_0327 + 0.0416666666666667*IT_0382 +
       0.0416666666666667*IT_0458) + 24*IT_0313*(conj(IT_0286) + -1./4*conj
      (IT_0327) + 0.0416666666666667*IT_0383 + 0.0416666666666667*IT_0459) + (-6
      )*IT_0330*IT_0512;
    const ccomplex_t IT_0515 = IT_0004*((conj(IT_0028) + conj(IT_0034))
      *IT_0075 + (IT_0028 + IT_0034)*IT_0078 + (conj(IT_0081) + conj(IT_0083))
      *IT_0105 + (IT_0081 + IT_0083)*IT_0108) + (IT_0078*(IT_0081 + IT_0083) +
       IT_0075*(conj(IT_0081) + conj(IT_0083)) + (conj(IT_0028) + conj(IT_0034))
      *IT_0105 + (IT_0028 + IT_0034)*IT_0108)*IT_0112 + IT_0113*(IT_0078*IT_0168
       + IT_0075*conj(IT_0168) + IT_0108*IT_0207 + IT_0105*conj(IT_0207)) + 
      (IT_0108*IT_0168 + IT_0105*conj(IT_0168) + IT_0078*IT_0207 + IT_0075*conj
      (IT_0207))*IT_0208 + IT_0209*(conj(IT_0216)*IT_0219 + IT_0216*IT_0222 +
       conj(IT_0227)*IT_0230 + IT_0227*IT_0233) + (IT_0222*IT_0227 + IT_0219
      *conj(IT_0227) + conj(IT_0216)*IT_0230 + IT_0216*IT_0233)*IT_0234 +
       IT_0235*(IT_0219*IT_0238 + IT_0222*IT_0239 + IT_0230*IT_0242 + IT_0233
      *IT_0243) + (IT_0230*IT_0238 + IT_0233*IT_0239 + IT_0219*IT_0242 + IT_0222
      *IT_0243)*IT_0244 + IT_0250*(conj(IT_0257)*IT_0290 + IT_0257*IT_0293 +
       conj(IT_0300)*IT_0317 + IT_0300*IT_0320) + IT_0322*(conj(IT_0327)*IT_0330
       + IT_0327*IT_0333 + conj(IT_0338)*IT_0341 + IT_0338*IT_0344) + IT_0346
      *IT_0379 + 6*IT_0381*(conj(IT_0300)*IT_0338 + IT_0300*conj(IT_0338) +
       0.166666666666667*conj(IT_0327)*IT_0382 + 0.166666666666667*IT_0327
      *IT_0383) + IT_0385*IT_0390 + IT_0412 + 6*IT_0413*IT_0422 + IT_0444 +
       IT_0445*IT_0448 + IT_0451*(24*IT_0071*conj(IT_0073) + IT_0073*(24*conj
      (IT_0071) + (-24)*conj(IT_0073)) + IT_0103*(24*conj(IT_0101) + (-24)*conj
      (IT_0103)) + IT_0101*IT_0454 + conj(IT_0071)*IT_0455) + IT_0457*(24
      *IT_0288*conj(IT_0288) + IT_0315*((-24)*conj(IT_0313) + 24*conj(IT_0315)) 
      + IT_0313*(24*conj(IT_0313) + (-24)*conj(IT_0315)) + conj(IT_0286)*IT_0458
       + IT_0286*(24*conj(IT_0286) + IT_0459)) + (IT_0245 + IT_0248 + IT_0321 +
       IT_0380)*(IT_0257*IT_0462 + conj(IT_0300)*IT_0465 + IT_0300*IT_0468 + 6
      *conj(IT_0257)*(IT_0073 + 0.166666666666667*IT_0469)) + IT_0470*IT_0479 +
       IT_0245*IT_0483 + IT_0484*IT_0487 + IT_0246*IT_0498 + IT_0499*IT_0508 +
       IT_0248*IT_0511 + (IT_0245 + IT_0246 + IT_0249 + IT_0380)*(IT_0327
      *IT_0473 + (-6)*conj(IT_0327)*(IT_0073 + (-0.166666666666667)*IT_0477) + 6
      *IT_0338*(conj(IT_0101) + 0.166666666666667*IT_0478) + (-6)*IT_0476
      *IT_0512) + IT_0513*IT_0514;
    return create_ccomplex_return(IT_0515);
}

