#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_N_1_N_1_to_anti_t_t.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_N_1_N_1_to_anti_t_t(
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
    const creal_t m_N_1 = param->m_N_1;
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
    const ccomplex_t N_B1 = param->N_B1;
    const ccomplex_t N_W1 = param->N_W1;
    const ccomplex_t N_d1 = param->N_d1;
    const ccomplex_t N_u1 = param->N_u1;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_N_1, 2);
    const ccomplex_t IT_0001 = cpow((-2)*s_12 + m_Z*(m_Z + (0 + _Complex_I*-1)
      *Gamma_Z) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0002 = cos(theta_W);
    const ccomplex_t IT_0003 = cpow(IT_0002, -1);
    const ccomplex_t IT_0004 = sin(theta_W);
    const ccomplex_t IT_0005 = (0 + _Complex_I*1)*e_em*IT_0003*IT_0004;
    const ccomplex_t IT_0006 = (-0.666666666666667)*IT_0005;
    const ccomplex_t IT_0007 = IT_0001*IT_0006;
    const ccomplex_t IT_0008 = (0 + _Complex_I*1)*IT_0007;
    const ccomplex_t IT_0009 = IT_0003*IT_0004;
    const ccomplex_t IT_0010 = N_u1*conj(N_u1)*e_em;
    const ccomplex_t IT_0011 = IT_0009*IT_0010;
    const ccomplex_t IT_0012 = cpow(IT_0004, -1);
    const ccomplex_t IT_0013 = IT_0002*IT_0012;
    const ccomplex_t IT_0014 = IT_0010*IT_0013;
    const ccomplex_t IT_0015 = N_d1*conj(N_d1)*e_em;
    const ccomplex_t IT_0016 = IT_0009*IT_0015;
    const ccomplex_t IT_0017 = IT_0013*IT_0015;
    const ccomplex_t IT_0018 = (0 + _Complex_I*1)*(IT_0011 + IT_0014 + 
      -IT_0016 + -IT_0017);
    const ccomplex_t IT_0019 = (-0.5)*IT_0018;
    const ccomplex_t IT_0020 = IT_0008*IT_0019;
    const ccomplex_t IT_0021 = (-0.5)*IT_0020;
    const ccomplex_t IT_0022 = N_B1*e_em*conj(U_st_01);
    const ccomplex_t IT_0023 = IT_0003*IT_0022;
    const ccomplex_t IT_0024 = 1.4142135623731*IT_0023;
    const ccomplex_t IT_0025 = N_W1*e_em*conj(U_st_01);
    const ccomplex_t IT_0026 = IT_0012*IT_0025;
    const ccomplex_t IT_0027 = 1.4142135623731*IT_0026;
    const ccomplex_t IT_0028 = sin(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = IT_0012*IT_0029;
    const ccomplex_t IT_0031 = pow(m_W, -1);
    const ccomplex_t IT_0032 = m_t*N_u1*e_em*IT_0031*conj(U_st_11);
    const ccomplex_t IT_0033 = IT_0030*IT_0032;
    const ccomplex_t IT_0034 = 1.4142135623731*IT_0033;
    const ccomplex_t IT_0035 = (0 + _Complex_I*1)*(IT_0024 + 3*IT_0027 + 3
      *IT_0034);
    const ccomplex_t IT_0036 = 0.166666666666667*IT_0035;
    const ccomplex_t IT_0037 = conj(N_B1)*e_em*U_st_01;
    const ccomplex_t IT_0038 = IT_0003*IT_0037;
    const ccomplex_t IT_0039 = 1.4142135623731*IT_0038;
    const ccomplex_t IT_0040 = conj(N_W1)*e_em*U_st_01;
    const ccomplex_t IT_0041 = IT_0012*IT_0040;
    const ccomplex_t IT_0042 = 1.4142135623731*IT_0041;
    const ccomplex_t IT_0043 = m_t*conj(N_u1)*e_em*IT_0031*U_st_11;
    const ccomplex_t IT_0044 = IT_0030*IT_0043;
    const ccomplex_t IT_0045 = 1.4142135623731*IT_0044;
    const ccomplex_t IT_0046 = (0 + _Complex_I*1)*(IT_0039 + 3*IT_0042 + 3
      *IT_0045);
    const ccomplex_t IT_0047 = 0.166666666666667*IT_0046;
    const ccomplex_t IT_0048 = IT_0036*IT_0047;
    const ccomplex_t IT_0049 = pow(m_t, 2);
    const ccomplex_t IT_0050 = cpow((-2)*s_13 + IT_0000 + IT_0049 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0051 = IT_0048*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*IT_0051;
    const ccomplex_t IT_0053 = cpow((-2)*s_13 + IT_0000 + IT_0049 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0054 = N_B1*e_em*conj(U_st_00);
    const ccomplex_t IT_0055 = IT_0003*IT_0054;
    const ccomplex_t IT_0056 = 1.4142135623731*IT_0055;
    const ccomplex_t IT_0057 = N_W1*e_em*conj(U_st_00);
    const ccomplex_t IT_0058 = IT_0012*IT_0057;
    const ccomplex_t IT_0059 = 1.4142135623731*IT_0058;
    const ccomplex_t IT_0060 = m_t*N_u1*e_em*IT_0031*conj(U_st_10);
    const ccomplex_t IT_0061 = IT_0030*IT_0060;
    const ccomplex_t IT_0062 = 1.4142135623731*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0056 + 3*IT_0059 + 3
      *IT_0062);
    const ccomplex_t IT_0064 = 0.166666666666667*IT_0063;
    const ccomplex_t IT_0065 = conj(N_B1)*e_em*U_st_00;
    const ccomplex_t IT_0066 = IT_0003*IT_0065;
    const ccomplex_t IT_0067 = 1.4142135623731*IT_0066;
    const ccomplex_t IT_0068 = conj(N_W1)*e_em*U_st_00;
    const ccomplex_t IT_0069 = IT_0012*IT_0068;
    const ccomplex_t IT_0070 = 1.4142135623731*IT_0069;
    const ccomplex_t IT_0071 = m_t*conj(N_u1)*e_em*IT_0031*U_st_10;
    const ccomplex_t IT_0072 = IT_0030*IT_0071;
    const ccomplex_t IT_0073 = 1.4142135623731*IT_0072;
    const ccomplex_t IT_0074 = (0 + _Complex_I*1)*(IT_0067 + 3*IT_0070 + 3
      *IT_0073);
    const ccomplex_t IT_0075 = 0.166666666666667*IT_0074;
    const ccomplex_t IT_0076 = IT_0064*IT_0075;
    const ccomplex_t IT_0077 = IT_0053*IT_0076;
    const ccomplex_t IT_0078 = (0 + _Complex_I*1)*IT_0077;
    const ccomplex_t IT_0079 = -IT_0052 + -IT_0078;
    const ccomplex_t IT_0080 = s_12*IT_0049;
    const ccomplex_t IT_0081 = 12*IT_0080;
    const ccomplex_t IT_0082 = IT_0079*IT_0081;
    const ccomplex_t IT_0083 = cpow((-2)*s_12 + (-2)*IT_0000 + m_A0*(m_A0 + (0
       + _Complex_I*-1)*Gamma_A0) + -reg_prop, -1);
    const ccomplex_t IT_0084 = cos(beta);
    const ccomplex_t IT_0085 = m_t*e_em*IT_0012*IT_0029*IT_0031*IT_0084;
    const ccomplex_t IT_0086 = (-0.5)*IT_0085;
    const ccomplex_t IT_0087 = IT_0083*IT_0086;
    const ccomplex_t IT_0088 = (0 + _Complex_I*1)*IT_0087;
    const ccomplex_t IT_0089 = IT_0003*IT_0084;
    const ccomplex_t IT_0090 = conj(N_B1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0091 = IT_0089*IT_0090;
    const ccomplex_t IT_0092 = IT_0003*IT_0028;
    const ccomplex_t IT_0093 = conj(N_B1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0094 = IT_0092*IT_0093;
    const ccomplex_t IT_0095 = IT_0012*IT_0084;
    const ccomplex_t IT_0096 = conj(N_W1)*conj(N_u1)*e_em;
    const ccomplex_t IT_0097 = IT_0095*IT_0096;
    const ccomplex_t IT_0098 = IT_0012*IT_0028;
    const ccomplex_t IT_0099 = conj(N_W1)*conj(N_d1)*e_em;
    const ccomplex_t IT_0100 = IT_0098*IT_0099;
    const ccomplex_t IT_0101 = IT_0091 + -IT_0094 + -IT_0097 + IT_0100;
    const ccomplex_t IT_0102 = IT_0088*IT_0101;
    const ccomplex_t IT_0103 = pow(m_Z, 2);
    const ccomplex_t IT_0104 = cpow((-2)*s_12 + (-2)*IT_0000 + IT_0103 + 
      -reg_prop, -1);
    const ccomplex_t IT_0105 = m_t*e_em*IT_0012*IT_0031;
    const ccomplex_t IT_0106 = (-0.5)*IT_0105;
    const ccomplex_t IT_0107 = IT_0104*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = IT_0089*IT_0093;
    const ccomplex_t IT_0110 = IT_0090*IT_0092;
    const ccomplex_t IT_0111 = IT_0095*IT_0099;
    const ccomplex_t IT_0112 = IT_0096*IT_0098;
    const ccomplex_t IT_0113 = IT_0109 + IT_0110 + -IT_0111 + -IT_0112;
    const ccomplex_t IT_0114 = IT_0108*IT_0113;
    const ccomplex_t IT_0115 = 0.5*IT_0102 + 0.5*IT_0114;
    const ccomplex_t IT_0116 = m_t*m_N_1;
    const ccomplex_t IT_0117 = s_23 + s_24;
    const ccomplex_t IT_0118 = IT_0116*IT_0117;
    const ccomplex_t IT_0119 = 12*IT_0118;
    const ccomplex_t IT_0120 = IT_0115*IT_0119;
    const ccomplex_t IT_0121 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0122 = IT_0003*IT_0121;
    const ccomplex_t IT_0123 = 1.4142135623731*IT_0122;
    const ccomplex_t IT_0124 = m_t*N_u1*e_em*IT_0031*U_st_01;
    const ccomplex_t IT_0125 = IT_0030*IT_0124;
    const ccomplex_t IT_0126 = 1.4142135623731*IT_0125;
    const ccomplex_t IT_0127 = (0 + _Complex_I*1)*(IT_0123 + (-0.75)*IT_0126);
    const ccomplex_t IT_0128 = 0.666666666666667*IT_0127;
    const ccomplex_t IT_0129 = IT_0036*IT_0128;
    const ccomplex_t IT_0130 = IT_0050*IT_0129;
    const ccomplex_t IT_0131 = (0 + _Complex_I*1)*IT_0130;
    const ccomplex_t IT_0132 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0133 = IT_0003*IT_0132;
    const ccomplex_t IT_0134 = 1.4142135623731*IT_0133;
    const ccomplex_t IT_0135 = m_t*N_u1*e_em*IT_0031*U_st_00;
    const ccomplex_t IT_0136 = IT_0030*IT_0135;
    const ccomplex_t IT_0137 = 1.4142135623731*IT_0136;
    const ccomplex_t IT_0138 = (0 + _Complex_I*1)*(IT_0134 + (-0.75)*IT_0137);
    const ccomplex_t IT_0139 = 0.666666666666667*IT_0138;
    const ccomplex_t IT_0140 = IT_0064*IT_0139;
    const ccomplex_t IT_0141 = IT_0053*IT_0140;
    const ccomplex_t IT_0142 = (0 + _Complex_I*1)*IT_0141;
    const ccomplex_t IT_0143 = IT_0131 + IT_0142;
    const ccomplex_t IT_0144 = m_t*s_14*m_N_1;
    const ccomplex_t IT_0145 = 12*IT_0144;
    const ccomplex_t IT_0146 = IT_0143*IT_0145;
    const ccomplex_t IT_0147 = conj(N_B1)*e_em*conj(U_st_11);
    const ccomplex_t IT_0148 = IT_0003*IT_0147;
    const ccomplex_t IT_0149 = 1.4142135623731*IT_0148;
    const ccomplex_t IT_0150 = m_t*conj(N_u1)*e_em*IT_0031*conj(U_st_01);
    const ccomplex_t IT_0151 = IT_0030*IT_0150;
    const ccomplex_t IT_0152 = 1.4142135623731*IT_0151;
    const ccomplex_t IT_0153 = (0 + _Complex_I*1)*(IT_0149 + (-0.75)*IT_0152);
    const ccomplex_t IT_0154 = 0.666666666666667*IT_0153;
    const ccomplex_t IT_0155 = IT_0047*IT_0154;
    const ccomplex_t IT_0156 = IT_0050*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*IT_0156;
    const ccomplex_t IT_0158 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0159 = IT_0003*IT_0158;
    const ccomplex_t IT_0160 = 1.4142135623731*IT_0159;
    const ccomplex_t IT_0161 = m_t*conj(N_u1)*e_em*IT_0031*conj(U_st_00);
    const ccomplex_t IT_0162 = IT_0030*IT_0161;
    const ccomplex_t IT_0163 = 1.4142135623731*IT_0162;
    const ccomplex_t IT_0164 = (0 + _Complex_I*1)*(IT_0160 + (-0.75)*IT_0163);
    const ccomplex_t IT_0165 = 0.666666666666667*IT_0164;
    const ccomplex_t IT_0166 = IT_0075*IT_0165;
    const ccomplex_t IT_0167 = IT_0053*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = IT_0157 + IT_0168;
    const ccomplex_t IT_0170 = m_t*s_23*m_N_1;
    const ccomplex_t IT_0171 = 12*IT_0170;
    const ccomplex_t IT_0172 = IT_0169*IT_0171;
    const ccomplex_t IT_0173 = IT_0128*IT_0154;
    const ccomplex_t IT_0174 = IT_0050*IT_0173;
    const ccomplex_t IT_0175 = (0 + _Complex_I*1)*IT_0174;
    const ccomplex_t IT_0176 = IT_0139*IT_0165;
    const ccomplex_t IT_0177 = IT_0053*IT_0176;
    const ccomplex_t IT_0178 = (0 + _Complex_I*1)*IT_0177;
    const ccomplex_t IT_0179 = -IT_0175 + -IT_0178;
    const ccomplex_t IT_0180 = s_34*IT_0000;
    const ccomplex_t IT_0181 = 12*IT_0180;
    const ccomplex_t IT_0182 = IT_0179*IT_0181;
    const ccomplex_t IT_0183 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0184 = IT_0089*IT_0183;
    const ccomplex_t IT_0185 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0186 = IT_0092*IT_0185;
    const ccomplex_t IT_0187 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0188 = IT_0095*IT_0187;
    const ccomplex_t IT_0189 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0190 = IT_0098*IT_0189;
    const ccomplex_t IT_0191 = -IT_0184 + IT_0186 + IT_0188 + -IT_0190;
    const ccomplex_t IT_0192 = IT_0088*IT_0191;
    const ccomplex_t IT_0193 = IT_0089*IT_0185;
    const ccomplex_t IT_0194 = IT_0092*IT_0183;
    const ccomplex_t IT_0195 = IT_0095*IT_0189;
    const ccomplex_t IT_0196 = IT_0098*IT_0187;
    const ccomplex_t IT_0197 = -IT_0193 + -IT_0194 + IT_0195 + IT_0196;
    const ccomplex_t IT_0198 = IT_0108*IT_0197;
    const ccomplex_t IT_0199 = (-0.5)*IT_0192 + (-0.5)*IT_0198;
    const ccomplex_t IT_0200 = s_13 + s_14;
    const ccomplex_t IT_0201 = IT_0116*IT_0200;
    const ccomplex_t IT_0202 = 12*IT_0201;
    const ccomplex_t IT_0203 = IT_0199*IT_0202;
    const ccomplex_t IT_0204 = (-0.5)*IT_0102 + (-0.5)*IT_0114;
    const ccomplex_t IT_0205 = (-12)*IT_0118;
    const ccomplex_t IT_0206 = IT_0204*IT_0205;
    const ccomplex_t IT_0207 = e_em*IT_0009;
    const ccomplex_t IT_0208 = e_em*IT_0013;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*(IT_0207 + (-3)*IT_0208);
    const ccomplex_t IT_0210 = (-0.166666666666667)*IT_0209;
    const ccomplex_t IT_0211 = IT_0001*IT_0210;
    const ccomplex_t IT_0212 = (0 + _Complex_I*1)*IT_0211;
    const ccomplex_t IT_0213 = IT_0019*IT_0212;
    const ccomplex_t IT_0214 = (-0.5)*IT_0213;
    const ccomplex_t IT_0215 = IT_0119*IT_0204;
    const ccomplex_t IT_0216 = cpow((-2)*s_23 + IT_0000 + IT_0049 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0217 = IT_0166*IT_0216;
    const ccomplex_t IT_0218 = (0 + _Complex_I*1)*IT_0217;
    const ccomplex_t IT_0219 = cpow((-2)*s_23 + IT_0000 + IT_0049 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0220 = IT_0155*IT_0219;
    const ccomplex_t IT_0221 = (0 + _Complex_I*1)*IT_0220;
    const ccomplex_t IT_0222 = -IT_0218 + -IT_0221;
    const ccomplex_t IT_0223 = m_t*s_24*m_N_1;
    const ccomplex_t IT_0224 = 12*IT_0223;
    const ccomplex_t IT_0225 = IT_0222*IT_0224;
    const ccomplex_t IT_0226 = IT_0048*IT_0219;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = IT_0076*IT_0216;
    const ccomplex_t IT_0229 = (0 + _Complex_I*1)*IT_0228;
    const ccomplex_t IT_0230 = IT_0227 + IT_0229;
    const ccomplex_t IT_0231 = IT_0181*IT_0230;
    const ccomplex_t IT_0232 = IT_0115*IT_0205;
    const ccomplex_t IT_0233 = (-12)*IT_0201;
    const ccomplex_t IT_0234 = IT_0199*IT_0233;
    const ccomplex_t IT_0235 = 0.5*IT_0020;
    const ccomplex_t IT_0236 = IT_0081*IT_0230;
    const ccomplex_t IT_0237 = m_t*s_13*m_N_1;
    const ccomplex_t IT_0238 = 12*IT_0237;
    const ccomplex_t IT_0239 = IT_0222*IT_0238;
    const ccomplex_t IT_0240 = IT_0202*IT_0204;
    const ccomplex_t IT_0241 = 0.5*IT_0018;
    const ccomplex_t IT_0242 = IT_0008*IT_0241;
    const ccomplex_t IT_0243 = (-0.5)*IT_0242;
    const ccomplex_t IT_0244 = s_13*s_24;
    const ccomplex_t IT_0245 = 48*IT_0244;
    const ccomplex_t IT_0246 = IT_0243*IT_0245;
    const ccomplex_t IT_0247 = IT_0199*IT_0205;
    const ccomplex_t IT_0248 = IT_0115*IT_0233;
    const ccomplex_t IT_0249 = 0.5*IT_0213;
    const ccomplex_t IT_0250 = IT_0081*IT_0179;
    const ccomplex_t IT_0251 = IT_0119*IT_0199;
    const ccomplex_t IT_0252 = IT_0145*IT_0169;
    const ccomplex_t IT_0253 = IT_0143*IT_0171;
    const ccomplex_t IT_0254 = IT_0079*IT_0181;
    const ccomplex_t IT_0255 = IT_0115*IT_0202;
    const ccomplex_t IT_0256 = IT_0000*IT_0049;
    const ccomplex_t IT_0257 = 48*IT_0256;
    const ccomplex_t IT_0258 = IT_0021*IT_0257;
    const ccomplex_t IT_0259 = IT_0204*IT_0233;
    const ccomplex_t IT_0260 = 0.5*IT_0242;
    const ccomplex_t IT_0261 = 24*conj(IT_0260);
    const ccomplex_t IT_0262 = 24*conj(IT_0021);
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 24*conj(IT_0249);
    const ccomplex_t IT_0265 = IT_0212*IT_0241;
    const ccomplex_t IT_0266 = (-0.5)*IT_0265;
    const ccomplex_t IT_0267 = 24*conj(IT_0266);
    const ccomplex_t IT_0268 = IT_0264 + IT_0267;
    const ccomplex_t IT_0269 = s_12*s_34;
    const ccomplex_t IT_0270 = s_14*s_23;
    const ccomplex_t IT_0271 = -IT_0244;
    const ccomplex_t IT_0272 = IT_0080 + IT_0270 + IT_0271;
    const ccomplex_t IT_0273 = IT_0269 + IT_0272;
    const ccomplex_t IT_0274 = 6*IT_0273;
    const ccomplex_t IT_0275 = -IT_0269;
    const ccomplex_t IT_0276 = IT_0244 + IT_0270 + IT_0275;
    const ccomplex_t IT_0277 = 6*IT_0276;
    const ccomplex_t IT_0278 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0279 = cos(alpha);
    const ccomplex_t IT_0280 = (0 + _Complex_I*1)*m_t*e_em*IT_0012*IT_0029
      *IT_0031*IT_0279;
    const ccomplex_t IT_0281 = (-0.5)*IT_0280;
    const ccomplex_t IT_0282 = IT_0278*IT_0281;
    const ccomplex_t IT_0283 = (0 + _Complex_I*1)*IT_0282;
    const ccomplex_t IT_0284 = IT_0003*IT_0279;
    const ccomplex_t IT_0285 = IT_0090*IT_0284;
    const ccomplex_t IT_0286 = sin(alpha);
    const ccomplex_t IT_0287 = IT_0003*IT_0286;
    const ccomplex_t IT_0288 = IT_0093*IT_0287;
    const ccomplex_t IT_0289 = IT_0012*IT_0279;
    const ccomplex_t IT_0290 = IT_0096*IT_0289;
    const ccomplex_t IT_0291 = IT_0012*IT_0286;
    const ccomplex_t IT_0292 = IT_0099*IT_0291;
    const ccomplex_t IT_0293 = (0 + _Complex_I*1)*(IT_0285 + IT_0288 + 
      -IT_0290 + -IT_0292);
    const ccomplex_t IT_0294 = -IT_0293;
    const ccomplex_t IT_0295 = IT_0283*IT_0294;
    const ccomplex_t IT_0296 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0297 = (0 + _Complex_I*1)*m_t*e_em*IT_0012*IT_0029
      *IT_0031*IT_0286;
    const ccomplex_t IT_0298 = (-0.5)*IT_0297;
    const ccomplex_t IT_0299 = IT_0296*IT_0298;
    const ccomplex_t IT_0300 = (0 + _Complex_I*1)*IT_0299;
    const ccomplex_t IT_0301 = IT_0093*IT_0284;
    const ccomplex_t IT_0302 = IT_0090*IT_0287;
    const ccomplex_t IT_0303 = IT_0099*IT_0289;
    const ccomplex_t IT_0304 = IT_0096*IT_0291;
    const ccomplex_t IT_0305 = (0 + _Complex_I*1)*(IT_0301 + -IT_0302 + 
      -IT_0303 + IT_0304);
    const ccomplex_t IT_0306 = IT_0300*IT_0305;
    const ccomplex_t IT_0307 = (-0.5)*IT_0295 + (-0.5)*IT_0306;
    const ccomplex_t IT_0308 = -IT_0270;
    const ccomplex_t IT_0309 = IT_0080 + IT_0244 + IT_0308;
    const ccomplex_t IT_0310 = IT_0275 + IT_0309;
    const ccomplex_t IT_0311 = 6*IT_0310;
    const ccomplex_t IT_0312 = s_34 + IT_0049;
    const ccomplex_t IT_0313 = IT_0000*IT_0312;
    const ccomplex_t IT_0314 = 6*IT_0313;
    const ccomplex_t IT_0315 = 6*IT_0237;
    const ccomplex_t IT_0316 = 6*IT_0223;
    const ccomplex_t IT_0317 = 6*IT_0256;
    const ccomplex_t IT_0318 = IT_0183*IT_0284;
    const ccomplex_t IT_0319 = IT_0185*IT_0287;
    const ccomplex_t IT_0320 = IT_0187*IT_0289;
    const ccomplex_t IT_0321 = IT_0189*IT_0291;
    const ccomplex_t IT_0322 = (0 + _Complex_I*1)*(IT_0318 + IT_0319 + 
      -IT_0320 + -IT_0321);
    const ccomplex_t IT_0323 = -IT_0322;
    const ccomplex_t IT_0324 = IT_0283*IT_0323;
    const ccomplex_t IT_0325 = IT_0185*IT_0284;
    const ccomplex_t IT_0326 = IT_0183*IT_0287;
    const ccomplex_t IT_0327 = IT_0189*IT_0289;
    const ccomplex_t IT_0328 = IT_0187*IT_0291;
    const ccomplex_t IT_0329 = (0 + _Complex_I*1)*(IT_0325 + -IT_0326 + 
      -IT_0327 + IT_0328);
    const ccomplex_t IT_0330 = IT_0300*IT_0329;
    const ccomplex_t IT_0331 = 0.5*IT_0324 + 0.5*IT_0330;
    const ccomplex_t IT_0332 = -s_34;
    const ccomplex_t IT_0333 = IT_0049 + IT_0332;
    const ccomplex_t IT_0334 = IT_0000*IT_0333;
    const ccomplex_t IT_0335 = 6*IT_0334;
    const ccomplex_t IT_0336 = 12*IT_0270;
    const ccomplex_t IT_0337 = IT_0173*IT_0219;
    const ccomplex_t IT_0338 = (0 + _Complex_I*1)*IT_0337;
    const ccomplex_t IT_0339 = IT_0176*IT_0216;
    const ccomplex_t IT_0340 = (0 + _Complex_I*1)*IT_0339;
    const ccomplex_t IT_0341 = IT_0338 + IT_0340;
    const ccomplex_t IT_0342 = IT_0129*IT_0219;
    const ccomplex_t IT_0343 = (0 + _Complex_I*1)*IT_0342;
    const ccomplex_t IT_0344 = IT_0140*IT_0216;
    const ccomplex_t IT_0345 = (0 + _Complex_I*1)*IT_0344;
    const ccomplex_t IT_0346 = -IT_0343 + -IT_0345;
    const ccomplex_t IT_0347 = 12*IT_0256;
    const ccomplex_t IT_0348 = (-6)*IT_0273;
    const ccomplex_t IT_0349 = 0.5*IT_0295 + 0.5*IT_0306;
    const ccomplex_t IT_0350 = (-6)*IT_0310;
    const ccomplex_t IT_0351 = 0.5*IT_0192 + 0.5*IT_0198;
    const ccomplex_t IT_0352 = (-6)*IT_0313;
    const ccomplex_t IT_0353 = (-0.5)*IT_0324 + (-0.5)*IT_0330;
    const ccomplex_t IT_0354 = (-6)*IT_0334;
    const ccomplex_t IT_0355 = IT_0145*conj(IT_0230) + IT_0170*IT_0263 +
       IT_0144*IT_0268 + conj(IT_0115)*IT_0274 + conj(IT_0169)*IT_0277 + conj
      (IT_0307)*IT_0311 + conj(IT_0199)*IT_0314 + conj(IT_0179)*IT_0315 + conj
      (IT_0079)*IT_0316 + conj(IT_0143)*IT_0317 + conj(IT_0331)*IT_0335 + conj
      (IT_0222)*IT_0336 + IT_0171*conj(IT_0341) + conj(IT_0346)*IT_0347 + conj
      (IT_0204)*IT_0348 + conj(IT_0349)*IT_0350 + conj(IT_0351)*IT_0352 + conj
      (IT_0353)*IT_0354;
    const ccomplex_t IT_0356 = 0.5*IT_0265;
    const ccomplex_t IT_0357 = 24*conj(IT_0356);
    const ccomplex_t IT_0358 = 24*conj(IT_0214);
    const ccomplex_t IT_0359 = IT_0357 + IT_0358;
    const ccomplex_t IT_0360 = 24*conj(IT_0235);
    const ccomplex_t IT_0361 = 24*conj(IT_0243);
    const ccomplex_t IT_0362 = IT_0360 + IT_0361;
    const ccomplex_t IT_0363 = 6*IT_0269 + 6*IT_0309;
    const ccomplex_t IT_0364 = IT_0272 + IT_0275;
    const ccomplex_t IT_0365 = 6*IT_0364;
    const ccomplex_t IT_0366 = 6*IT_0144;
    const ccomplex_t IT_0367 = 6*IT_0170;
    const ccomplex_t IT_0368 = 12*IT_0244;
    const ccomplex_t IT_0369 = (-6)*IT_0269 + (-6)*IT_0309;
    const ccomplex_t IT_0370 = (-6)*IT_0364;
    const ccomplex_t IT_0371 = conj(IT_0079)*IT_0224 + conj(IT_0179)*IT_0238 +
       conj(IT_0222)*IT_0277 + IT_0317*conj(IT_0346) + conj(IT_0143)*IT_0347 +
       IT_0314*conj(IT_0351) + conj(IT_0199)*IT_0352 + IT_0335*conj(IT_0353) +
       conj(IT_0331)*IT_0354 + IT_0223*IT_0359 + IT_0237*IT_0362 + conj(IT_0204)
      *IT_0363 + conj(IT_0349)*IT_0365 + conj(IT_0230)*IT_0366 + conj(IT_0341)
      *IT_0367 + conj(IT_0169)*IT_0368 + conj(IT_0115)*IT_0369 + conj(IT_0307)
      *IT_0370;
    const ccomplex_t IT_0372 = conj(IT_0179)*IT_0224 + conj(IT_0079)*IT_0238 +
       conj(IT_0204)*IT_0314 + conj(IT_0222)*IT_0317 + conj(IT_0307)*IT_0335 +
       IT_0277*conj(IT_0346) + conj(IT_0169)*IT_0347 + conj(IT_0115)*IT_0352 +
       conj(IT_0349)*IT_0354 + IT_0237*IT_0359 + IT_0223*IT_0362 + conj(IT_0351)
      *IT_0363 + conj(IT_0331)*IT_0365 + conj(IT_0341)*IT_0366 + conj(IT_0230)
      *IT_0367 + conj(IT_0143)*IT_0368 + conj(IT_0199)*IT_0369 + conj(IT_0353)
      *IT_0370;
    const ccomplex_t IT_0373 = conj(IT_0214) + conj(IT_0356);
    const ccomplex_t IT_0374 = conj(IT_0235) + conj(IT_0243);
    const ccomplex_t IT_0375 = IT_0171*conj(IT_0230) + IT_0144*IT_0263 +
       IT_0170*IT_0268 + conj(IT_0199)*IT_0274 + conj(IT_0143)*IT_0277 + conj
      (IT_0115)*IT_0314 + conj(IT_0079)*IT_0315 + conj(IT_0179)*IT_0316 + conj
      (IT_0169)*IT_0317 + IT_0145*conj(IT_0341) + IT_0336*conj(IT_0346) + conj
      (IT_0222)*IT_0347 + IT_0335*conj(IT_0349) + conj(IT_0331)*IT_0350 +
       IT_0348*conj(IT_0351) + conj(IT_0204)*IT_0352 + IT_0311*conj(IT_0353) +
       conj(IT_0307)*IT_0354 + IT_0238*IT_0373 + IT_0224*IT_0374;
    const ccomplex_t IT_0376 = -s_14;
    const ccomplex_t IT_0377 = s_13 + IT_0376;
    const ccomplex_t IT_0378 = IT_0116*IT_0377;
    const ccomplex_t IT_0379 = (-2)*IT_0378;
    const ccomplex_t IT_0380 = 3*IT_0379;
    const ccomplex_t IT_0381 = -s_24;
    const ccomplex_t IT_0382 = s_23 + IT_0381;
    const ccomplex_t IT_0383 = IT_0116*IT_0382;
    const ccomplex_t IT_0384 = 2*IT_0383;
    const ccomplex_t IT_0385 = 3*IT_0384;
    const ccomplex_t IT_0386 = 6*IT_0080;
    const ccomplex_t IT_0387 = 6*IT_0118;
    const ccomplex_t IT_0388 = 6*IT_0180;
    const ccomplex_t IT_0389 = 6*IT_0201;
    const ccomplex_t IT_0390 = (-6)*IT_0118;
    const ccomplex_t IT_0391 = (-6)*IT_0201;
    const ccomplex_t IT_0392 = (-3)*IT_0379;
    const ccomplex_t IT_0393 = (-3)*IT_0384;
    const ccomplex_t IT_0394 = IT_0145*conj(IT_0222) + IT_0256*IT_0263 +
       IT_0268*IT_0270 + conj(IT_0230)*IT_0336 + IT_0171*conj(IT_0346) + conj
      (IT_0341)*IT_0347 + conj(IT_0169)*IT_0366 + conj(IT_0143)*IT_0367 + conj
      (IT_0349)*IT_0380 + conj(IT_0331)*IT_0385 + conj(IT_0179)*IT_0386 + conj
      (IT_0199)*IT_0387 + conj(IT_0079)*IT_0388 + conj(IT_0115)*IT_0389 + conj
      (IT_0351)*IT_0390 + conj(IT_0204)*IT_0391 + conj(IT_0307)*IT_0392 + conj
      (IT_0353)*IT_0393;
    const ccomplex_t IT_0395 = IT_0171*conj(IT_0222) + IT_0256*IT_0268 +
       IT_0263*IT_0270 + IT_0336*conj(IT_0341) + IT_0145*conj(IT_0346) + conj
      (IT_0230)*IT_0347 + conj(IT_0143)*IT_0366 + conj(IT_0169)*IT_0367 +
       IT_0081*IT_0373 + IT_0181*IT_0374 + conj(IT_0331)*IT_0380 + conj(IT_0349)
      *IT_0385 + conj(IT_0079)*IT_0386 + conj(IT_0115)*IT_0387 + conj(IT_0179)
      *IT_0388 + conj(IT_0199)*IT_0389 + conj(IT_0204)*IT_0390 + conj(IT_0351)
      *IT_0391 + conj(IT_0353)*IT_0392 + conj(IT_0307)*IT_0393;
    const ccomplex_t IT_0396 = 2*IT_0378;
    const ccomplex_t IT_0397 = 3*IT_0396;
    const ccomplex_t IT_0398 = (-2)*IT_0383;
    const ccomplex_t IT_0399 = 3*IT_0398;
    const ccomplex_t IT_0400 = (-3)*IT_0396;
    const ccomplex_t IT_0401 = (-3)*IT_0398;
    const ccomplex_t IT_0402 = conj(IT_0143)*IT_0224 + conj(IT_0169)*IT_0238 +
       conj(IT_0222)*IT_0315 + IT_0316*conj(IT_0346) + conj(IT_0079)*IT_0347 +
       IT_0256*IT_0359 + IT_0244*IT_0362 + conj(IT_0179)*IT_0368 + conj(IT_0230)
      *IT_0386 + conj(IT_0351)*IT_0387 + conj(IT_0341)*IT_0388 + conj(IT_0204)
      *IT_0389 + conj(IT_0199)*IT_0390 + conj(IT_0115)*IT_0391 + conj(IT_0307)
      *IT_0397 + conj(IT_0353)*IT_0399 + conj(IT_0349)*IT_0400 + conj(IT_0331)
      *IT_0401;
    const ccomplex_t IT_0403 = conj(IT_0169)*IT_0224 + conj(IT_0143)*IT_0238 +
       conj(IT_0222)*IT_0316 + IT_0315*conj(IT_0346) + conj(IT_0179)*IT_0347 +
       IT_0244*IT_0359 + IT_0256*IT_0362 + conj(IT_0079)*IT_0368 + conj(IT_0341)
      *IT_0386 + conj(IT_0204)*IT_0387 + conj(IT_0230)*IT_0388 + conj(IT_0351)
      *IT_0389 + conj(IT_0115)*IT_0390 + conj(IT_0199)*IT_0391 + conj(IT_0353)
      *IT_0397 + conj(IT_0307)*IT_0399 + conj(IT_0331)*IT_0400 + conj(IT_0349)
      *IT_0401;
    const ccomplex_t IT_0404 = (-4)*IT_0378;
    const ccomplex_t IT_0405 = 3*IT_0404;
    const ccomplex_t IT_0406 = IT_0349*IT_0405;
    const ccomplex_t IT_0407 = (-3)*IT_0404;
    const ccomplex_t IT_0408 = IT_0307*IT_0407;
    const ccomplex_t IT_0409 = (-4)*IT_0383;
    const ccomplex_t IT_0410 = 3*IT_0409;
    const ccomplex_t IT_0411 = IT_0307*IT_0410;
    const ccomplex_t IT_0412 = (-3)*IT_0409;
    const ccomplex_t IT_0413 = IT_0349*IT_0412;
    const ccomplex_t IT_0414 = 24*IT_0313;
    const ccomplex_t IT_0415 = s_12*IT_0312;
    const ccomplex_t IT_0416 = (-24)*IT_0415;
    const ccomplex_t IT_0417 = (-24)*IT_0313;
    const ccomplex_t IT_0418 = conj(IT_0249) + conj(IT_0266);
    const ccomplex_t IT_0419 = conj(IT_0021) + conj(IT_0260);
    const ccomplex_t IT_0420 = 24*IT_0415;
    const ccomplex_t IT_0421 = conj(IT_0169)*IT_0314 + conj(IT_0346)*IT_0348 +
       conj(IT_0222)*IT_0352 + conj(IT_0143)*IT_0363 + IT_0202*IT_0373 + IT_0119
      *IT_0374 + conj(IT_0179)*IT_0387 + conj(IT_0079)*IT_0389 + conj(IT_0230)
      *IT_0390 + conj(IT_0341)*IT_0391 + conj(IT_0204)*IT_0414 + conj(IT_0199)
      *IT_0416 + conj(IT_0115)*IT_0417 + IT_0205*IT_0418 + IT_0233*IT_0419 +
       conj(IT_0351)*IT_0420;
    const ccomplex_t IT_0422 = conj(IT_0222)*IT_0314 + IT_0274*conj(IT_0346) +
       conj(IT_0169)*IT_0352 + conj(IT_0143)*IT_0369 + conj(IT_0230)*IT_0387 +
       conj(IT_0341)*IT_0389 + conj(IT_0179)*IT_0390 + conj(IT_0079)*IT_0391 +
       conj(IT_0115)*IT_0414 + conj(IT_0351)*IT_0416 + conj(IT_0204)*IT_0417 +
       conj(IT_0199)*IT_0420;
    const ccomplex_t IT_0423 = conj(IT_0222)*IT_0274 + IT_0314*conj(IT_0346) +
       conj(IT_0143)*IT_0352 + conj(IT_0169)*IT_0369 + conj(IT_0341)*IT_0387 +
       conj(IT_0230)*IT_0389 + conj(IT_0079)*IT_0390 + conj(IT_0179)*IT_0391 +
       conj(IT_0199)*IT_0414 + conj(IT_0204)*IT_0416 + conj(IT_0351)*IT_0417 +
       conj(IT_0115)*IT_0420;
    const ccomplex_t IT_0424 = conj(IT_0143)*IT_0314 + conj(IT_0222)*IT_0348 +
       conj(IT_0346)*IT_0352 + conj(IT_0169)*IT_0363 + conj(IT_0079)*IT_0387 +
       conj(IT_0179)*IT_0389 + conj(IT_0341)*IT_0390 + conj(IT_0230)*IT_0391 +
       conj(IT_0351)*IT_0414 + conj(IT_0115)*IT_0416 + conj(IT_0199)*IT_0417 +
       conj(IT_0204)*IT_0420;
    const ccomplex_t IT_0425 = conj(IT_0235) + conj(IT_0249);
    const ccomplex_t IT_0426 = 4*IT_0378;
    const ccomplex_t IT_0427 = 3*IT_0426;
    const ccomplex_t IT_0428 = s_12*IT_0333;
    const ccomplex_t IT_0429 = (-24)*IT_0428;
    const ccomplex_t IT_0430 = (-24)*IT_0334;
    const ccomplex_t IT_0431 = conj(IT_0021) + conj(IT_0214);
    const ccomplex_t IT_0432 = 4*IT_0383;
    const ccomplex_t IT_0433 = (-3)*IT_0432;
    const ccomplex_t IT_0434 = 24*IT_0428;
    const ccomplex_t IT_0435 = 24*IT_0334;
    const ccomplex_t IT_0436 = conj(IT_0222)*IT_0311 + conj(IT_0143)*IT_0335 +
       conj(IT_0346)*IT_0354 + conj(IT_0169)*IT_0370 + conj(IT_0230)*IT_0392 +
       conj(IT_0341)*IT_0393 + conj(IT_0179)*IT_0397 + conj(IT_0079)*IT_0399 +
       IT_0425*IT_0427 + conj(IT_0307)*IT_0429 + conj(IT_0331)*IT_0430 + IT_0431
      *IT_0433 + conj(IT_0349)*IT_0434 + conj(IT_0353)*IT_0435;
    const ccomplex_t IT_0437 = 3*IT_0432;
    const ccomplex_t IT_0438 = (-3)*IT_0426;
    const ccomplex_t IT_0439 = IT_0335*conj(IT_0346) + conj(IT_0222)*IT_0350 +
       conj(IT_0143)*IT_0354 + conj(IT_0169)*IT_0365 + conj(IT_0230)*IT_0380 +
       conj(IT_0341)*IT_0385 + conj(IT_0179)*IT_0400 + conj(IT_0079)*IT_0401 +
       conj(IT_0349)*IT_0429 + conj(IT_0353)*IT_0430 + conj(IT_0307)*IT_0434 +
       conj(IT_0331)*IT_0435 + IT_0431*IT_0437 + IT_0425*IT_0438;
    const ccomplex_t IT_0440 = conj(IT_0260) + conj(IT_0356);
    const ccomplex_t IT_0441 = conj(IT_0222)*IT_0335 + conj(IT_0346)*IT_0350 +
       conj(IT_0169)*IT_0354 + conj(IT_0143)*IT_0365 + conj(IT_0341)*IT_0380 +
       conj(IT_0230)*IT_0385 + conj(IT_0079)*IT_0400 + conj(IT_0179)*IT_0401 +
       IT_0412*IT_0425 + conj(IT_0331)*IT_0429 + conj(IT_0307)*IT_0430 + IT_0405
      *IT_0431 + conj(IT_0353)*IT_0434 + conj(IT_0349)*IT_0435 + (conj(IT_0243) 
      + conj(IT_0266))*IT_0437 + IT_0438*IT_0440;
    const ccomplex_t IT_0442 = conj(IT_0169)*IT_0335 + IT_0311*conj(IT_0346) +
       conj(IT_0222)*IT_0354 + conj(IT_0143)*IT_0370 + conj(IT_0341)*IT_0392 +
       conj(IT_0230)*IT_0393 + conj(IT_0079)*IT_0397 + conj(IT_0179)*IT_0399 +
       IT_0410*IT_0425 + conj(IT_0353)*IT_0429 + conj(IT_0349)*IT_0430 + IT_0407
      *IT_0431 + (conj(IT_0243) + conj(IT_0266))*IT_0433 + conj(IT_0331)*IT_0434
       + conj(IT_0307)*IT_0435 + IT_0427*IT_0440;
    const ccomplex_t IT_0443 = conj(IT_0331)*IT_0405;
    const ccomplex_t IT_0444 = conj(IT_0349)*IT_0437;
    const ccomplex_t IT_0445 = conj(IT_0079)*IT_0081;
    const ccomplex_t IT_0446 = conj(IT_0115)*IT_0119;
    const ccomplex_t IT_0447 = conj(IT_0143)*IT_0145;
    const ccomplex_t IT_0448 = conj(IT_0169)*IT_0171;
    const ccomplex_t IT_0449 = conj(IT_0179)*IT_0181;
    const ccomplex_t IT_0450 = conj(IT_0199)*IT_0202;
    const ccomplex_t IT_0451 = 48*IT_0270;
    const ccomplex_t IT_0452 = conj(IT_0204)*IT_0205;
    const ccomplex_t IT_0453 = IT_0233*conj(IT_0351);
    const ccomplex_t IT_0454 = conj(IT_0353)*IT_0407;
    const ccomplex_t IT_0455 = conj(IT_0307)*IT_0433;
    const ccomplex_t IT_0456 = 24*IT_0270;
    const ccomplex_t IT_0457 = conj(IT_0341)*IT_0456;
    const ccomplex_t IT_0458 = 24*IT_0144;
    const ccomplex_t IT_0459 = 24*IT_0170;
    const ccomplex_t IT_0460 = conj(IT_0222)*IT_0459;
    const ccomplex_t IT_0461 = 24*IT_0256;
    const ccomplex_t IT_0462 = conj(IT_0230)*IT_0461;
    const ccomplex_t IT_0463 = IT_0080*IT_0359;
    const ccomplex_t IT_0464 = IT_0180*IT_0362;
    const ccomplex_t IT_0465 = IT_0443 + IT_0444 + IT_0445 + IT_0446 + IT_0447
       + IT_0448 + IT_0449 + IT_0450 + IT_0419*IT_0451 + IT_0452 + IT_0453 +
       IT_0454 + IT_0455 + IT_0457 + conj(IT_0346)*IT_0458 + IT_0460 + IT_0462 +
       IT_0463 + IT_0464;
    const ccomplex_t IT_0466 = conj(IT_0307)*IT_0410;
    const ccomplex_t IT_0467 = conj(IT_0353)*IT_0427;
    const ccomplex_t IT_0468 = conj(IT_0349)*IT_0412;
    const ccomplex_t IT_0469 = conj(IT_0331)*IT_0438;
    const ccomplex_t IT_0470 = IT_0257*IT_0418 + IT_0445 + IT_0446 + IT_0447 +
       IT_0448 + IT_0449 + IT_0450 + IT_0419*IT_0451 + IT_0452 + IT_0453 +
       IT_0457 + conj(IT_0346)*IT_0458 + IT_0460 + IT_0462 + IT_0463 + IT_0464 +
       IT_0466 + IT_0467 + IT_0468 + IT_0469;
    const ccomplex_t IT_0471 = conj(IT_0349)*IT_0405;
    const ccomplex_t IT_0472 = conj(IT_0331)*IT_0437;
    const ccomplex_t IT_0473 = IT_0081*conj(IT_0179);
    const ccomplex_t IT_0474 = IT_0119*conj(IT_0199);
    const ccomplex_t IT_0475 = IT_0145*conj(IT_0169);
    const ccomplex_t IT_0476 = conj(IT_0143)*IT_0171;
    const ccomplex_t IT_0477 = conj(IT_0079)*IT_0181;
    const ccomplex_t IT_0478 = conj(IT_0115)*IT_0202;
    const ccomplex_t IT_0479 = conj(IT_0021)*IT_0257;
    const ccomplex_t IT_0480 = IT_0205*conj(IT_0351);
    const ccomplex_t IT_0481 = conj(IT_0204)*IT_0233;
    const ccomplex_t IT_0482 = conj(IT_0307)*IT_0407;
    const ccomplex_t IT_0483 = conj(IT_0353)*IT_0433;
    const ccomplex_t IT_0484 = conj(IT_0222)*IT_0458;
    const ccomplex_t IT_0485 = IT_0080*IT_0362;
    const ccomplex_t IT_0486 = IT_0180*IT_0359;
    const ccomplex_t IT_0487 = IT_0257*conj(IT_0260) + IT_0418*IT_0451 + conj
      (IT_0230)*IT_0456 + conj(IT_0346)*IT_0459 + conj(IT_0341)*IT_0461 +
       IT_0471 + IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 +
       IT_0478 + IT_0479 + IT_0480 + IT_0481 + IT_0482 + IT_0483 + IT_0484 +
       IT_0485 + IT_0486;
    const ccomplex_t IT_0488 = conj(IT_0353)*IT_0410;
    const ccomplex_t IT_0489 = conj(IT_0307)*IT_0427;
    const ccomplex_t IT_0490 = conj(IT_0331)*IT_0412;
    const ccomplex_t IT_0491 = conj(IT_0349)*IT_0438;
    const ccomplex_t IT_0492 = IT_0257*conj(IT_0260) + IT_0418*IT_0451 + conj
      (IT_0230)*IT_0456 + conj(IT_0346)*IT_0459 + conj(IT_0341)*IT_0461 +
       IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 + IT_0478 + IT_0479 +
       IT_0480 + IT_0481 + IT_0484 + IT_0485 + IT_0486 + IT_0488 + IT_0489 +
       IT_0490 + IT_0491;
    const ccomplex_t IT_0493 = IT_0081*conj(IT_0230);
    const ccomplex_t IT_0494 = IT_0119*conj(IT_0351);
    const ccomplex_t IT_0495 = conj(IT_0222)*IT_0238;
    const ccomplex_t IT_0496 = IT_0224*conj(IT_0346);
    const ccomplex_t IT_0497 = IT_0181*conj(IT_0341);
    const ccomplex_t IT_0498 = IT_0202*conj(IT_0204);
    const ccomplex_t IT_0499 = conj(IT_0214)*IT_0257;
    const ccomplex_t IT_0500 = conj(IT_0199)*IT_0205;
    const ccomplex_t IT_0501 = conj(IT_0115)*IT_0233;
    const ccomplex_t IT_0502 = 24*IT_0244;
    const ccomplex_t IT_0503 = conj(IT_0179)*IT_0502;
    const ccomplex_t IT_0504 = 24*IT_0237;
    const ccomplex_t IT_0505 = conj(IT_0169)*IT_0504;
    const ccomplex_t IT_0506 = 24*IT_0223;
    const ccomplex_t IT_0507 = conj(IT_0079)*IT_0461;
    const ccomplex_t IT_0508 = IT_0180*IT_0263;
    const ccomplex_t IT_0509 = IT_0080*IT_0268 + IT_0257*conj(IT_0356) +
       IT_0471 + IT_0472 + IT_0482 + IT_0483 + IT_0493 + IT_0494 + IT_0495 +
       IT_0496 + IT_0497 + IT_0498 + IT_0499 + IT_0500 + IT_0501 + IT_0503 +
       IT_0505 + conj(IT_0143)*IT_0506 + IT_0507 + IT_0508;
    const ccomplex_t IT_0510 = IT_0080*IT_0268 + IT_0257*conj(IT_0356) +
       IT_0245*IT_0374 + IT_0488 + IT_0489 + IT_0490 + IT_0491 + IT_0493 +
       IT_0494 + IT_0495 + IT_0496 + IT_0497 + IT_0498 + IT_0499 + IT_0500 +
       IT_0501 + IT_0503 + IT_0505 + conj(IT_0143)*IT_0506 + IT_0507 + IT_0508;
    const ccomplex_t IT_0511 = IT_0119*conj(IT_0204);
    const ccomplex_t IT_0512 = IT_0238*conj(IT_0346);
    const ccomplex_t IT_0513 = conj(IT_0222)*IT_0224;
    const ccomplex_t IT_0514 = IT_0181*conj(IT_0230);
    const ccomplex_t IT_0515 = IT_0202*conj(IT_0351);
    const ccomplex_t IT_0516 = conj(IT_0199)*IT_0233;
    const ccomplex_t IT_0517 = conj(IT_0169)*IT_0506;
    const ccomplex_t IT_0518 = conj(IT_0179)*IT_0461;
    const ccomplex_t IT_0519 = IT_0080*IT_0263;
    const ccomplex_t IT_0520 = conj(IT_0115)*IT_0205;
    const ccomplex_t IT_0521 = IT_0257*IT_0374;
    const ccomplex_t IT_0522 = IT_0245*IT_0373;
    const ccomplex_t IT_0523 = IT_0081*conj(IT_0341);
    const ccomplex_t IT_0524 = IT_0180*IT_0268 + IT_0443 + IT_0444 + IT_0454 +
       IT_0455 + conj(IT_0079)*IT_0502 + conj(IT_0143)*IT_0504 + IT_0511 +
       IT_0512 + IT_0513 + IT_0514 + IT_0515 + IT_0516 + IT_0517 + IT_0518 +
       IT_0519 + IT_0520 + IT_0521 + IT_0522 + IT_0523;
    const ccomplex_t IT_0525 = IT_0180*IT_0268 + IT_0466 + IT_0467 + IT_0468 +
       IT_0469 + conj(IT_0079)*IT_0502 + conj(IT_0143)*IT_0504 + IT_0511 +
       IT_0512 + IT_0513 + IT_0514 + IT_0515 + IT_0516 + IT_0517 + IT_0518 +
       IT_0519 + IT_0520 + IT_0521 + IT_0522 + IT_0523;
    const ccomplex_t IT_0526 = conj(IT_0021)*(IT_0082 + IT_0120 + IT_0146 +
       IT_0172 + IT_0182 + IT_0203 + IT_0206) + conj(IT_0214)*(IT_0215 + IT_0225
       + IT_0231 + IT_0232 + IT_0234) + conj(IT_0235)*(IT_0236 + IT_0239 +
       IT_0240 + IT_0246 + IT_0247 + IT_0248) + conj(IT_0249)*(IT_0250 + IT_0251
       + IT_0252 + IT_0253 + IT_0254 + IT_0255 + IT_0258 + IT_0259) + IT_0222
      *IT_0355 + IT_0169*IT_0371 + IT_0143*IT_0372 + IT_0346*IT_0375 + IT_0230
      *IT_0394 + IT_0341*IT_0395 + IT_0179*IT_0402 + IT_0079*IT_0403 + conj
      (IT_0243)*(IT_0236 + IT_0239 + IT_0240 + IT_0246 + IT_0247 + IT_0248 +
       IT_0406 + IT_0408) + conj(IT_0266)*(IT_0250 + IT_0251 + IT_0252 + IT_0253
       + IT_0254 + IT_0255 + IT_0258 + IT_0259 + IT_0406 + IT_0408) + conj
      (IT_0260)*(IT_0082 + IT_0120 + IT_0146 + IT_0172 + IT_0182 + IT_0203 +
       IT_0206 + IT_0411 + IT_0413) + conj(IT_0356)*(IT_0215 + IT_0225 + IT_0231
       + IT_0232 + IT_0234 + IT_0411 + IT_0413) + IT_0351*IT_0421 + IT_0199
      *IT_0422 + IT_0115*IT_0423 + IT_0204*IT_0424 + IT_0307*IT_0436 + IT_0349
      *IT_0439 + IT_0331*IT_0441 + IT_0353*IT_0442 + IT_0021*IT_0465 + IT_0260
      *IT_0470 + IT_0266*IT_0487 + IT_0249*IT_0492 + IT_0243*IT_0509 + IT_0235
      *IT_0510 + IT_0214*IT_0524 + IT_0356*IT_0525;
    return create_ccomplex_return(IT_0526);
}

