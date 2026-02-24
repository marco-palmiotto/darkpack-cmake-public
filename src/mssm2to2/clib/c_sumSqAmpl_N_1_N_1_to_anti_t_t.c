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
    const ccomplex_t IT_0121 = N_B1*N_d1*e_em;
    const ccomplex_t IT_0122 = IT_0089*IT_0121;
    const ccomplex_t IT_0123 = N_B1*N_u1*e_em;
    const ccomplex_t IT_0124 = IT_0092*IT_0123;
    const ccomplex_t IT_0125 = N_W1*N_d1*e_em;
    const ccomplex_t IT_0126 = IT_0095*IT_0125;
    const ccomplex_t IT_0127 = N_W1*N_u1*e_em;
    const ccomplex_t IT_0128 = IT_0098*IT_0127;
    const ccomplex_t IT_0129 = -IT_0122 + -IT_0124 + IT_0126 + IT_0128;
    const ccomplex_t IT_0130 = IT_0108*IT_0129;
    const ccomplex_t IT_0131 = IT_0089*IT_0123;
    const ccomplex_t IT_0132 = IT_0092*IT_0121;
    const ccomplex_t IT_0133 = IT_0095*IT_0127;
    const ccomplex_t IT_0134 = IT_0098*IT_0125;
    const ccomplex_t IT_0135 = -IT_0131 + IT_0132 + IT_0133 + -IT_0134;
    const ccomplex_t IT_0136 = IT_0088*IT_0135;
    const ccomplex_t IT_0137 = (-0.5)*IT_0130 + (-0.5)*IT_0136;
    const ccomplex_t IT_0138 = s_13 + s_14;
    const ccomplex_t IT_0139 = IT_0116*IT_0138;
    const ccomplex_t IT_0140 = 12*IT_0139;
    const ccomplex_t IT_0141 = IT_0137*IT_0140;
    const ccomplex_t IT_0142 = conj(N_B1)*e_em*conj(U_st_10);
    const ccomplex_t IT_0143 = IT_0003*IT_0142;
    const ccomplex_t IT_0144 = 1.4142135623731*IT_0143;
    const ccomplex_t IT_0145 = m_t*conj(N_u1)*e_em*IT_0031*conj(U_st_00);
    const ccomplex_t IT_0146 = IT_0030*IT_0145;
    const ccomplex_t IT_0147 = 1.4142135623731*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*(IT_0144 + (-0.75)*IT_0147);
    const ccomplex_t IT_0149 = 0.666666666666667*IT_0148;
    const ccomplex_t IT_0150 = IT_0075*IT_0149;
    const ccomplex_t IT_0151 = IT_0053*IT_0150;
    const ccomplex_t IT_0152 = (0 + _Complex_I*1)*IT_0151;
    const ccomplex_t IT_0153 = conj(N_B1)*e_em*conj(U_st_11);
    const ccomplex_t IT_0154 = IT_0003*IT_0153;
    const ccomplex_t IT_0155 = 1.4142135623731*IT_0154;
    const ccomplex_t IT_0156 = m_t*conj(N_u1)*e_em*IT_0031*conj(U_st_01);
    const ccomplex_t IT_0157 = IT_0030*IT_0156;
    const ccomplex_t IT_0158 = 1.4142135623731*IT_0157;
    const ccomplex_t IT_0159 = (0 + _Complex_I*1)*(IT_0155 + (-0.75)*IT_0158);
    const ccomplex_t IT_0160 = 0.666666666666667*IT_0159;
    const ccomplex_t IT_0161 = IT_0047*IT_0160;
    const ccomplex_t IT_0162 = IT_0050*IT_0161;
    const ccomplex_t IT_0163 = (0 + _Complex_I*1)*IT_0162;
    const ccomplex_t IT_0164 = IT_0152 + IT_0163;
    const ccomplex_t IT_0165 = m_t*s_23*m_N_1;
    const ccomplex_t IT_0166 = 12*IT_0165;
    const ccomplex_t IT_0167 = IT_0164*IT_0166;
    const ccomplex_t IT_0168 = (-0.5)*IT_0102 + (-0.5)*IT_0114;
    const ccomplex_t IT_0169 = (-12)*IT_0118;
    const ccomplex_t IT_0170 = IT_0168*IT_0169;
    const ccomplex_t IT_0171 = e_em*IT_0009;
    const ccomplex_t IT_0172 = e_em*IT_0013;
    const ccomplex_t IT_0173 = (0 + _Complex_I*1)*(IT_0171 + (-3)*IT_0172);
    const ccomplex_t IT_0174 = (-0.166666666666667)*IT_0173;
    const ccomplex_t IT_0175 = IT_0001*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = IT_0019*IT_0176;
    const ccomplex_t IT_0178 = (-0.5)*IT_0177;
    const ccomplex_t IT_0179 = cpow((-2)*s_23 + IT_0000 + IT_0049 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0180 = N_B1*e_em*U_st_10;
    const ccomplex_t IT_0181 = IT_0003*IT_0180;
    const ccomplex_t IT_0182 = 1.4142135623731*IT_0181;
    const ccomplex_t IT_0183 = m_t*N_u1*e_em*IT_0031*U_st_00;
    const ccomplex_t IT_0184 = IT_0030*IT_0183;
    const ccomplex_t IT_0185 = 1.4142135623731*IT_0184;
    const ccomplex_t IT_0186 = (0 + _Complex_I*1)*(IT_0182 + (-0.75)*IT_0185);
    const ccomplex_t IT_0187 = 0.666666666666667*IT_0186;
    const ccomplex_t IT_0188 = IT_0149*IT_0187;
    const ccomplex_t IT_0189 = IT_0179*IT_0188;
    const ccomplex_t IT_0190 = (0 + _Complex_I*1)*IT_0189;
    const ccomplex_t IT_0191 = N_B1*e_em*U_st_11;
    const ccomplex_t IT_0192 = IT_0003*IT_0191;
    const ccomplex_t IT_0193 = 1.4142135623731*IT_0192;
    const ccomplex_t IT_0194 = m_t*N_u1*e_em*IT_0031*U_st_01;
    const ccomplex_t IT_0195 = IT_0030*IT_0194;
    const ccomplex_t IT_0196 = 1.4142135623731*IT_0195;
    const ccomplex_t IT_0197 = (0 + _Complex_I*1)*(IT_0193 + (-0.75)*IT_0196);
    const ccomplex_t IT_0198 = 0.666666666666667*IT_0197;
    const ccomplex_t IT_0199 = IT_0160*IT_0198;
    const ccomplex_t IT_0200 = cpow((-2)*s_23 + IT_0000 + IT_0049 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0201 = IT_0199*IT_0200;
    const ccomplex_t IT_0202 = (0 + _Complex_I*1)*IT_0201;
    const ccomplex_t IT_0203 = IT_0190 + IT_0202;
    const ccomplex_t IT_0204 = IT_0081*IT_0203;
    const ccomplex_t IT_0205 = IT_0119*IT_0168;
    const ccomplex_t IT_0206 = IT_0161*IT_0200;
    const ccomplex_t IT_0207 = (0 + _Complex_I*1)*IT_0206;
    const ccomplex_t IT_0208 = IT_0150*IT_0179;
    const ccomplex_t IT_0209 = (0 + _Complex_I*1)*IT_0208;
    const ccomplex_t IT_0210 = -IT_0207 + -IT_0209;
    const ccomplex_t IT_0211 = m_t*s_24*m_N_1;
    const ccomplex_t IT_0212 = 12*IT_0211;
    const ccomplex_t IT_0213 = IT_0210*IT_0212;
    const ccomplex_t IT_0214 = IT_0064*IT_0187;
    const ccomplex_t IT_0215 = IT_0179*IT_0214;
    const ccomplex_t IT_0216 = (0 + _Complex_I*1)*IT_0215;
    const ccomplex_t IT_0217 = IT_0036*IT_0198;
    const ccomplex_t IT_0218 = IT_0200*IT_0217;
    const ccomplex_t IT_0219 = (0 + _Complex_I*1)*IT_0218;
    const ccomplex_t IT_0220 = -IT_0216 + -IT_0219;
    const ccomplex_t IT_0221 = m_t*s_13*m_N_1;
    const ccomplex_t IT_0222 = 12*IT_0221;
    const ccomplex_t IT_0223 = IT_0220*IT_0222;
    const ccomplex_t IT_0224 = IT_0048*IT_0200;
    const ccomplex_t IT_0225 = (0 + _Complex_I*1)*IT_0224;
    const ccomplex_t IT_0226 = IT_0076*IT_0179;
    const ccomplex_t IT_0227 = (0 + _Complex_I*1)*IT_0226;
    const ccomplex_t IT_0228 = IT_0225 + IT_0227;
    const ccomplex_t IT_0229 = s_34*IT_0000;
    const ccomplex_t IT_0230 = 12*IT_0229;
    const ccomplex_t IT_0231 = IT_0228*IT_0230;
    const ccomplex_t IT_0232 = 0.5*IT_0018;
    const ccomplex_t IT_0233 = IT_0176*IT_0232;
    const ccomplex_t IT_0234 = 0.5*IT_0233;
    const ccomplex_t IT_0235 = s_13*s_24;
    const ccomplex_t IT_0236 = 48*IT_0235;
    const ccomplex_t IT_0237 = IT_0234*IT_0236;
    const ccomplex_t IT_0238 = 0.5*IT_0020;
    const ccomplex_t IT_0239 = IT_0000*IT_0049;
    const ccomplex_t IT_0240 = 48*IT_0239;
    const ccomplex_t IT_0241 = IT_0238*IT_0240;
    const ccomplex_t IT_0242 = IT_0115*IT_0169;
    const ccomplex_t IT_0243 = (-12)*IT_0139;
    const ccomplex_t IT_0244 = IT_0137*IT_0243;
    const ccomplex_t IT_0245 = 0.5*IT_0177;
    const ccomplex_t IT_0246 = IT_0119*IT_0137;
    const ccomplex_t IT_0247 = IT_0115*IT_0140;
    const ccomplex_t IT_0248 = m_t*s_14*m_N_1;
    const ccomplex_t IT_0249 = 12*IT_0248;
    const ccomplex_t IT_0250 = IT_0164*IT_0249;
    const ccomplex_t IT_0251 = IT_0079*IT_0230;
    const ccomplex_t IT_0252 = IT_0168*IT_0243;
    const ccomplex_t IT_0253 = IT_0081*IT_0228;
    const ccomplex_t IT_0254 = IT_0140*IT_0168;
    const ccomplex_t IT_0255 = IT_0212*IT_0220;
    const ccomplex_t IT_0256 = IT_0210*IT_0222;
    const ccomplex_t IT_0257 = IT_0203*IT_0230;
    const ccomplex_t IT_0258 = IT_0137*IT_0169;
    const ccomplex_t IT_0259 = IT_0115*IT_0243;
    const ccomplex_t IT_0260 = (-0.5)*IT_0233;
    const ccomplex_t IT_0261 = 24*conj(IT_0260);
    const ccomplex_t IT_0262 = 24*conj(IT_0245);
    const ccomplex_t IT_0263 = IT_0261 + IT_0262;
    const ccomplex_t IT_0264 = 24*conj(IT_0021);
    const ccomplex_t IT_0265 = IT_0008*IT_0232;
    const ccomplex_t IT_0266 = 0.5*IT_0265;
    const ccomplex_t IT_0267 = 24*conj(IT_0266);
    const ccomplex_t IT_0268 = IT_0264 + IT_0267;
    const ccomplex_t IT_0269 = s_14*s_23;
    const ccomplex_t IT_0270 = s_12*s_34;
    const ccomplex_t IT_0271 = -IT_0235;
    const ccomplex_t IT_0272 = 6*IT_0080 + 6*IT_0269 + 6*IT_0270 + 6*IT_0271;
    const ccomplex_t IT_0273 = -IT_0270;
    const ccomplex_t IT_0274 = IT_0235 + IT_0269 + IT_0273;
    const ccomplex_t IT_0275 = 6*IT_0274;
    const ccomplex_t IT_0276 = cpow((-2)*s_12 + m_h*(m_h + (0 + _Complex_I*-1)
      *Gamma_h) + (-2)*IT_0000 + -reg_prop, -1);
    const ccomplex_t IT_0277 = cos(alpha);
    const ccomplex_t IT_0278 = (0 + _Complex_I*1)*m_t*e_em*IT_0012*IT_0029
      *IT_0031*IT_0277;
    const ccomplex_t IT_0279 = (-0.5)*IT_0278;
    const ccomplex_t IT_0280 = IT_0276*IT_0279;
    const ccomplex_t IT_0281 = (0 + _Complex_I*1)*IT_0280;
    const ccomplex_t IT_0282 = IT_0003*IT_0277;
    const ccomplex_t IT_0283 = IT_0090*IT_0282;
    const ccomplex_t IT_0284 = sin(alpha);
    const ccomplex_t IT_0285 = IT_0003*IT_0284;
    const ccomplex_t IT_0286 = IT_0093*IT_0285;
    const ccomplex_t IT_0287 = IT_0012*IT_0277;
    const ccomplex_t IT_0288 = IT_0096*IT_0287;
    const ccomplex_t IT_0289 = IT_0012*IT_0284;
    const ccomplex_t IT_0290 = IT_0099*IT_0289;
    const ccomplex_t IT_0291 = (0 + _Complex_I*1)*(IT_0283 + IT_0286 + 
      -IT_0288 + -IT_0290);
    const ccomplex_t IT_0292 = -IT_0291;
    const ccomplex_t IT_0293 = IT_0281*IT_0292;
    const ccomplex_t IT_0294 = cpow((-2)*s_12 + (-2)*IT_0000 + m_H0*(m_H0 + (0
       + _Complex_I*-1)*Gamma_H0) + -reg_prop, -1);
    const ccomplex_t IT_0295 = (0 + _Complex_I*1)*m_t*e_em*IT_0012*IT_0029
      *IT_0031*IT_0284;
    const ccomplex_t IT_0296 = (-0.5)*IT_0295;
    const ccomplex_t IT_0297 = IT_0294*IT_0296;
    const ccomplex_t IT_0298 = (0 + _Complex_I*1)*IT_0297;
    const ccomplex_t IT_0299 = IT_0093*IT_0282;
    const ccomplex_t IT_0300 = IT_0090*IT_0285;
    const ccomplex_t IT_0301 = IT_0099*IT_0287;
    const ccomplex_t IT_0302 = IT_0096*IT_0289;
    const ccomplex_t IT_0303 = (0 + _Complex_I*1)*(IT_0299 + -IT_0300 + 
      -IT_0301 + IT_0302);
    const ccomplex_t IT_0304 = IT_0298*IT_0303;
    const ccomplex_t IT_0305 = (-0.5)*IT_0293 + (-0.5)*IT_0304;
    const ccomplex_t IT_0306 = -IT_0269;
    const ccomplex_t IT_0307 = IT_0080 + IT_0235 + IT_0306;
    const ccomplex_t IT_0308 = IT_0273 + IT_0307;
    const ccomplex_t IT_0309 = 6*IT_0308;
    const ccomplex_t IT_0310 = 6*IT_0211;
    const ccomplex_t IT_0311 = IT_0053*IT_0188;
    const ccomplex_t IT_0312 = (0 + _Complex_I*1)*IT_0311;
    const ccomplex_t IT_0313 = IT_0050*IT_0199;
    const ccomplex_t IT_0314 = (0 + _Complex_I*1)*IT_0313;
    const ccomplex_t IT_0315 = -IT_0312 + -IT_0314;
    const ccomplex_t IT_0316 = 6*IT_0221;
    const ccomplex_t IT_0317 = IT_0053*IT_0214;
    const ccomplex_t IT_0318 = (0 + _Complex_I*1)*IT_0317;
    const ccomplex_t IT_0319 = IT_0050*IT_0217;
    const ccomplex_t IT_0320 = (0 + _Complex_I*1)*IT_0319;
    const ccomplex_t IT_0321 = IT_0318 + IT_0320;
    const ccomplex_t IT_0322 = 6*IT_0239;
    const ccomplex_t IT_0323 = IT_0123*IT_0282;
    const ccomplex_t IT_0324 = IT_0121*IT_0285;
    const ccomplex_t IT_0325 = IT_0127*IT_0287;
    const ccomplex_t IT_0326 = IT_0125*IT_0289;
    const ccomplex_t IT_0327 = (0 + _Complex_I*1)*(IT_0323 + IT_0324 + 
      -IT_0325 + -IT_0326);
    const ccomplex_t IT_0328 = -IT_0327;
    const ccomplex_t IT_0329 = IT_0281*IT_0328;
    const ccomplex_t IT_0330 = IT_0121*IT_0282;
    const ccomplex_t IT_0331 = IT_0123*IT_0285;
    const ccomplex_t IT_0332 = IT_0125*IT_0287;
    const ccomplex_t IT_0333 = IT_0127*IT_0289;
    const ccomplex_t IT_0334 = (0 + _Complex_I*1)*(IT_0330 + -IT_0331 + 
      -IT_0332 + IT_0333);
    const ccomplex_t IT_0335 = IT_0298*IT_0334;
    const ccomplex_t IT_0336 = 0.5*IT_0329 + 0.5*IT_0335;
    const ccomplex_t IT_0337 = -s_34;
    const ccomplex_t IT_0338 = IT_0049 + IT_0337;
    const ccomplex_t IT_0339 = IT_0000*IT_0338;
    const ccomplex_t IT_0340 = 6*IT_0339;
    const ccomplex_t IT_0341 = s_34 + IT_0049;
    const ccomplex_t IT_0342 = IT_0000*IT_0341;
    const ccomplex_t IT_0343 = 6*IT_0342;
    const ccomplex_t IT_0344 = 12*IT_0269;
    const ccomplex_t IT_0345 = 12*IT_0239;
    const ccomplex_t IT_0346 = (-6)*IT_0080 + (-6)*IT_0269 + (-6)*IT_0270 + (
      -6)*IT_0271;
    const ccomplex_t IT_0347 = 0.5*IT_0293 + 0.5*IT_0304;
    const ccomplex_t IT_0348 = (-6)*IT_0308;
    const ccomplex_t IT_0349 = (-0.5)*IT_0329 + (-0.5)*IT_0335;
    const ccomplex_t IT_0350 = (-6)*IT_0339;
    const ccomplex_t IT_0351 = 0.5*IT_0130 + 0.5*IT_0136;
    const ccomplex_t IT_0352 = (-6)*IT_0342;
    const ccomplex_t IT_0353 = IT_0166*conj(IT_0203) + conj(IT_0228)*IT_0249 +
       IT_0248*IT_0263 + IT_0165*IT_0268 + conj(IT_0115)*IT_0272 + conj(IT_0164)
      *IT_0275 + conj(IT_0305)*IT_0309 + conj(IT_0079)*IT_0310 + conj(IT_0315)
      *IT_0316 + conj(IT_0321)*IT_0322 + conj(IT_0336)*IT_0340 + conj(IT_0137)
      *IT_0343 + conj(IT_0210)*IT_0344 + conj(IT_0220)*IT_0345 + conj(IT_0168)
      *IT_0346 + conj(IT_0347)*IT_0348 + conj(IT_0349)*IT_0350 + conj(IT_0351)
      *IT_0352;
    const ccomplex_t IT_0354 = (-0.5)*IT_0265;
    const ccomplex_t IT_0355 = 24*conj(IT_0354);
    const ccomplex_t IT_0356 = 24*conj(IT_0238);
    const ccomplex_t IT_0357 = IT_0355 + IT_0356;
    const ccomplex_t IT_0358 = 24*conj(IT_0178);
    const ccomplex_t IT_0359 = 24*conj(IT_0234);
    const ccomplex_t IT_0360 = IT_0358 + IT_0359;
    const ccomplex_t IT_0361 = IT_0270 + IT_0307;
    const ccomplex_t IT_0362 = 6*IT_0361;
    const ccomplex_t IT_0363 = IT_0080 + IT_0269 + IT_0271 + IT_0273;
    const ccomplex_t IT_0364 = 6*IT_0363;
    const ccomplex_t IT_0365 = 6*IT_0165;
    const ccomplex_t IT_0366 = 6*IT_0248;
    const ccomplex_t IT_0367 = 12*IT_0235;
    const ccomplex_t IT_0368 = (-6)*IT_0361;
    const ccomplex_t IT_0369 = (-6)*IT_0363;
    const ccomplex_t IT_0370 = conj(IT_0079)*IT_0212 + conj(IT_0210)*IT_0275 +
       IT_0222*conj(IT_0315) + conj(IT_0220)*IT_0322 + conj(IT_0321)*IT_0345 +
       IT_0340*conj(IT_0349) + conj(IT_0336)*IT_0350 + IT_0343*conj(IT_0351) +
       conj(IT_0137)*IT_0352 + IT_0221*IT_0357 + IT_0211*IT_0360 + conj(IT_0168)
      *IT_0362 + conj(IT_0347)*IT_0364 + conj(IT_0203)*IT_0365 + conj(IT_0228)
      *IT_0366 + conj(IT_0164)*IT_0367 + conj(IT_0115)*IT_0368 + conj(IT_0305)
      *IT_0369;
    const ccomplex_t IT_0371 = IT_0166*conj(IT_0228) + conj(IT_0203)*IT_0249 +
       IT_0165*IT_0263 + IT_0248*IT_0268 + conj(IT_0137)*IT_0272 + IT_0310*conj
      (IT_0315) + conj(IT_0079)*IT_0316 + IT_0275*conj(IT_0321) + conj(IT_0164)
      *IT_0322 + conj(IT_0115)*IT_0343 + conj(IT_0220)*IT_0344 + conj(IT_0210)
      *IT_0345 + IT_0340*conj(IT_0347) + conj(IT_0336)*IT_0348 + IT_0309*conj
      (IT_0349) + conj(IT_0305)*IT_0350 + IT_0346*conj(IT_0351) + conj(IT_0168)
      *IT_0352;
    const ccomplex_t IT_0372 = conj(IT_0245) + conj(IT_0260);
    const ccomplex_t IT_0373 = conj(IT_0021) + conj(IT_0266);
    const ccomplex_t IT_0374 = conj(IT_0079)*IT_0222 + conj(IT_0220)*IT_0275 +
       IT_0212*conj(IT_0315) + conj(IT_0210)*IT_0322 + conj(IT_0305)*IT_0340 +
       conj(IT_0168)*IT_0343 + conj(IT_0164)*IT_0345 + conj(IT_0347)*IT_0350 +
       conj(IT_0115)*IT_0352 + IT_0211*IT_0357 + IT_0221*IT_0360 + conj(IT_0351)
      *IT_0362 + conj(IT_0336)*IT_0364 + conj(IT_0228)*IT_0365 + conj(IT_0203)
      *IT_0366 + conj(IT_0321)*IT_0367 + conj(IT_0137)*IT_0368 + conj(IT_0349)
      *IT_0369 + IT_0166*IT_0372 + IT_0249*IT_0373;
    const ccomplex_t IT_0375 = -s_24;
    const ccomplex_t IT_0376 = s_23 + IT_0375;
    const ccomplex_t IT_0377 = IT_0116*IT_0376;
    const ccomplex_t IT_0378 = (-2)*IT_0377;
    const ccomplex_t IT_0379 = 3*IT_0378;
    const ccomplex_t IT_0380 = -s_14;
    const ccomplex_t IT_0381 = s_13 + IT_0380;
    const ccomplex_t IT_0382 = IT_0116*IT_0381;
    const ccomplex_t IT_0383 = 2*IT_0382;
    const ccomplex_t IT_0384 = 3*IT_0383;
    const ccomplex_t IT_0385 = 6*IT_0080;
    const ccomplex_t IT_0386 = 6*IT_0118;
    const ccomplex_t IT_0387 = 6*IT_0139;
    const ccomplex_t IT_0388 = 6*IT_0229;
    const ccomplex_t IT_0389 = (-6)*IT_0118;
    const ccomplex_t IT_0390 = (-6)*IT_0139;
    const ccomplex_t IT_0391 = (-3)*IT_0378;
    const ccomplex_t IT_0392 = (-3)*IT_0383;
    const ccomplex_t IT_0393 = conj(IT_0164)*IT_0212 + conj(IT_0210)*IT_0310 +
       conj(IT_0220)*IT_0316 + IT_0222*conj(IT_0321) + conj(IT_0315)*IT_0345 +
       IT_0239*IT_0357 + IT_0235*IT_0360 + conj(IT_0079)*IT_0367 + conj(IT_0305)
      *IT_0379 + conj(IT_0349)*IT_0384 + conj(IT_0203)*IT_0385 + conj(IT_0168)
      *IT_0386 + conj(IT_0351)*IT_0387 + conj(IT_0228)*IT_0388 + conj(IT_0115)
      *IT_0389 + conj(IT_0137)*IT_0390 + conj(IT_0347)*IT_0391 + conj(IT_0336)
      *IT_0392;
    const ccomplex_t IT_0394 = conj(IT_0164)*IT_0222 + conj(IT_0220)*IT_0310 +
       conj(IT_0210)*IT_0316 + IT_0212*conj(IT_0321) + conj(IT_0079)*IT_0345 +
       IT_0235*IT_0357 + IT_0239*IT_0360 + conj(IT_0315)*IT_0367 + IT_0081
      *IT_0372 + IT_0230*IT_0373 + conj(IT_0349)*IT_0379 + conj(IT_0305)*IT_0384
       + conj(IT_0228)*IT_0385 + conj(IT_0351)*IT_0386 + conj(IT_0168)*IT_0387 +
       conj(IT_0203)*IT_0388 + conj(IT_0137)*IT_0389 + conj(IT_0115)*IT_0390 +
       conj(IT_0336)*IT_0391 + conj(IT_0347)*IT_0392;
    const ccomplex_t IT_0395 = 2*IT_0377;
    const ccomplex_t IT_0396 = 3*IT_0395;
    const ccomplex_t IT_0397 = (-2)*IT_0382;
    const ccomplex_t IT_0398 = 3*IT_0397;
    const ccomplex_t IT_0399 = (-3)*IT_0395;
    const ccomplex_t IT_0400 = (-3)*IT_0397;
    const ccomplex_t IT_0401 = IT_0166*conj(IT_0210) + conj(IT_0220)*IT_0249 +
       IT_0239*IT_0263 + IT_0268*IT_0269 + conj(IT_0203)*IT_0344 + conj(IT_0228)
      *IT_0345 + conj(IT_0164)*IT_0365 + conj(IT_0321)*IT_0366 + conj(IT_0079)
      *IT_0385 + conj(IT_0115)*IT_0386 + conj(IT_0137)*IT_0387 + conj(IT_0315)
      *IT_0388 + conj(IT_0168)*IT_0389 + conj(IT_0351)*IT_0390 + conj(IT_0347)
      *IT_0396 + conj(IT_0336)*IT_0398 + conj(IT_0305)*IT_0399 + conj(IT_0349)
      *IT_0400;
    const ccomplex_t IT_0402 = IT_0166*conj(IT_0220) + conj(IT_0210)*IT_0249 +
       IT_0239*IT_0268 + IT_0263*IT_0269 + conj(IT_0228)*IT_0344 + conj(IT_0203)
      *IT_0345 + conj(IT_0321)*IT_0365 + conj(IT_0164)*IT_0366 + conj(IT_0315)
      *IT_0385 + conj(IT_0137)*IT_0386 + conj(IT_0115)*IT_0387 + conj(IT_0079)
      *IT_0388 + conj(IT_0351)*IT_0389 + conj(IT_0168)*IT_0390 + conj(IT_0336)
      *IT_0396 + conj(IT_0347)*IT_0398 + conj(IT_0349)*IT_0399 + conj(IT_0305)
      *IT_0400;
    const ccomplex_t IT_0403 = (-4)*IT_0377;
    const ccomplex_t IT_0404 = 3*IT_0403;
    const ccomplex_t IT_0405 = IT_0305*IT_0404;
    const ccomplex_t IT_0406 = (-3)*IT_0403;
    const ccomplex_t IT_0407 = IT_0347*IT_0406;
    const ccomplex_t IT_0408 = (-4)*IT_0382;
    const ccomplex_t IT_0409 = 3*IT_0408;
    const ccomplex_t IT_0410 = IT_0347*IT_0409;
    const ccomplex_t IT_0411 = (-3)*IT_0408;
    const ccomplex_t IT_0412 = IT_0305*IT_0411;
    const ccomplex_t IT_0413 = conj(IT_0238) + conj(IT_0354);
    const ccomplex_t IT_0414 = conj(IT_0178) + conj(IT_0234);
    const ccomplex_t IT_0415 = 24*IT_0342;
    const ccomplex_t IT_0416 = s_12*IT_0341;
    const ccomplex_t IT_0417 = (-24)*IT_0416;
    const ccomplex_t IT_0418 = (-24)*IT_0342;
    const ccomplex_t IT_0419 = 24*IT_0416;
    const ccomplex_t IT_0420 = conj(IT_0164)*IT_0343 + conj(IT_0220)*IT_0346 +
       conj(IT_0210)*IT_0352 + conj(IT_0321)*IT_0362 + IT_0169*IT_0372 + IT_0243
      *IT_0373 + conj(IT_0315)*IT_0386 + conj(IT_0079)*IT_0387 + conj(IT_0228)
      *IT_0389 + conj(IT_0203)*IT_0390 + IT_0119*IT_0413 + IT_0140*IT_0414 +
       conj(IT_0168)*IT_0415 + conj(IT_0137)*IT_0417 + conj(IT_0115)*IT_0418 +
       conj(IT_0351)*IT_0419;
    const ccomplex_t IT_0421 = conj(IT_0220)*IT_0272 + conj(IT_0210)*IT_0343 +
       conj(IT_0164)*IT_0352 + conj(IT_0321)*IT_0368 + conj(IT_0228)*IT_0386 +
       conj(IT_0203)*IT_0387 + conj(IT_0315)*IT_0389 + conj(IT_0079)*IT_0390 +
       conj(IT_0115)*IT_0415 + conj(IT_0351)*IT_0417 + conj(IT_0168)*IT_0418 +
       conj(IT_0137)*IT_0419;
    const ccomplex_t IT_0422 = conj(IT_0321)*IT_0343 + conj(IT_0210)*IT_0346 +
       conj(IT_0220)*IT_0352 + conj(IT_0164)*IT_0362 + conj(IT_0079)*IT_0386 +
       conj(IT_0315)*IT_0387 + conj(IT_0203)*IT_0389 + conj(IT_0228)*IT_0390 +
       conj(IT_0351)*IT_0415 + conj(IT_0115)*IT_0417 + conj(IT_0137)*IT_0418 +
       conj(IT_0168)*IT_0419;
    const ccomplex_t IT_0423 = conj(IT_0210)*IT_0272 + conj(IT_0220)*IT_0343 +
       conj(IT_0321)*IT_0352 + conj(IT_0164)*IT_0368 + conj(IT_0203)*IT_0386 +
       conj(IT_0228)*IT_0387 + conj(IT_0079)*IT_0389 + conj(IT_0315)*IT_0390 +
       conj(IT_0137)*IT_0415 + conj(IT_0168)*IT_0417 + conj(IT_0351)*IT_0418 +
       conj(IT_0115)*IT_0419;
    const ccomplex_t IT_0424 = conj(IT_0238) + conj(IT_0245);
    const ccomplex_t IT_0425 = 4*IT_0382;
    const ccomplex_t IT_0426 = 3*IT_0425;
    const ccomplex_t IT_0427 = s_12*IT_0338;
    const ccomplex_t IT_0428 = (-24)*IT_0427;
    const ccomplex_t IT_0429 = (-24)*IT_0339;
    const ccomplex_t IT_0430 = conj(IT_0021) + conj(IT_0178);
    const ccomplex_t IT_0431 = 4*IT_0377;
    const ccomplex_t IT_0432 = (-3)*IT_0431;
    const ccomplex_t IT_0433 = 24*IT_0427;
    const ccomplex_t IT_0434 = 24*IT_0339;
    const ccomplex_t IT_0435 = conj(IT_0210)*IT_0309 + conj(IT_0321)*IT_0340 +
       conj(IT_0220)*IT_0350 + conj(IT_0164)*IT_0369 + conj(IT_0079)*IT_0379 +
       conj(IT_0315)*IT_0384 + conj(IT_0203)*IT_0399 + conj(IT_0228)*IT_0400 +
       IT_0424*IT_0426 + conj(IT_0305)*IT_0428 + conj(IT_0336)*IT_0429 + IT_0430
      *IT_0432 + conj(IT_0347)*IT_0433 + conj(IT_0349)*IT_0434;
    const ccomplex_t IT_0436 = 3*IT_0431;
    const ccomplex_t IT_0437 = (-3)*IT_0425;
    const ccomplex_t IT_0438 = conj(IT_0220)*IT_0340 + conj(IT_0210)*IT_0348 +
       conj(IT_0321)*IT_0350 + conj(IT_0164)*IT_0364 + conj(IT_0079)*IT_0391 +
       conj(IT_0315)*IT_0392 + conj(IT_0203)*IT_0396 + conj(IT_0228)*IT_0398 +
       conj(IT_0347)*IT_0428 + conj(IT_0349)*IT_0429 + conj(IT_0305)*IT_0433 +
       conj(IT_0336)*IT_0434 + IT_0430*IT_0436 + IT_0424*IT_0437;
    const ccomplex_t IT_0439 = conj(IT_0260) + conj(IT_0354);
    const ccomplex_t IT_0440 = conj(IT_0220)*IT_0309 + conj(IT_0164)*IT_0340 +
       conj(IT_0210)*IT_0350 + conj(IT_0321)*IT_0369 + conj(IT_0315)*IT_0379 +
       conj(IT_0079)*IT_0384 + conj(IT_0228)*IT_0399 + conj(IT_0203)*IT_0400 +
       IT_0404*IT_0424 + (conj(IT_0234) + conj(IT_0266))*IT_0426 + conj(IT_0349)
      *IT_0428 + conj(IT_0347)*IT_0429 + IT_0411*IT_0430 + conj(IT_0336)*IT_0433
       + conj(IT_0305)*IT_0434 + IT_0432*IT_0439;
    const ccomplex_t IT_0441 = conj(IT_0210)*IT_0340 + conj(IT_0220)*IT_0348 +
       conj(IT_0164)*IT_0350 + conj(IT_0321)*IT_0364 + conj(IT_0315)*IT_0391 +
       conj(IT_0079)*IT_0392 + conj(IT_0228)*IT_0396 + conj(IT_0203)*IT_0398 +
       IT_0406*IT_0424 + conj(IT_0336)*IT_0428 + conj(IT_0305)*IT_0429 + IT_0409
      *IT_0430 + conj(IT_0349)*IT_0433 + conj(IT_0347)*IT_0434 + (conj(IT_0234) 
      + conj(IT_0266))*IT_0437 + IT_0436*IT_0439;
    const ccomplex_t IT_0442 = conj(IT_0305)*IT_0404;
    const ccomplex_t IT_0443 = conj(IT_0349)*IT_0426;
    const ccomplex_t IT_0444 = IT_0081*conj(IT_0203);
    const ccomplex_t IT_0445 = IT_0119*conj(IT_0168);
    const ccomplex_t IT_0446 = IT_0140*conj(IT_0351);
    const ccomplex_t IT_0447 = conj(IT_0210)*IT_0212;
    const ccomplex_t IT_0448 = conj(IT_0220)*IT_0222;
    const ccomplex_t IT_0449 = conj(IT_0228)*IT_0230;
    const ccomplex_t IT_0450 = conj(IT_0238)*IT_0240;
    const ccomplex_t IT_0451 = conj(IT_0115)*IT_0169;
    const ccomplex_t IT_0452 = conj(IT_0137)*IT_0243;
    const ccomplex_t IT_0453 = conj(IT_0347)*IT_0406;
    const ccomplex_t IT_0454 = conj(IT_0336)*IT_0437;
    const ccomplex_t IT_0455 = 24*IT_0235;
    const ccomplex_t IT_0456 = conj(IT_0079)*IT_0455;
    const ccomplex_t IT_0457 = 24*IT_0211;
    const ccomplex_t IT_0458 = conj(IT_0164)*IT_0457;
    const ccomplex_t IT_0459 = 24*IT_0221;
    const ccomplex_t IT_0460 = 24*IT_0239;
    const ccomplex_t IT_0461 = IT_0080*IT_0268;
    const ccomplex_t IT_0462 = IT_0229*IT_0263;
    const ccomplex_t IT_0463 = IT_0240*conj(IT_0354) + IT_0442 + IT_0443 +
       IT_0444 + IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 + IT_0450 +
       IT_0451 + IT_0452 + IT_0453 + IT_0454 + IT_0456 + IT_0458 + conj(IT_0321)
      *IT_0459 + conj(IT_0315)*IT_0460 + IT_0461 + IT_0462;
    const ccomplex_t IT_0464 = conj(IT_0336)*IT_0409;
    const ccomplex_t IT_0465 = conj(IT_0347)*IT_0436;
    const ccomplex_t IT_0466 = conj(IT_0349)*IT_0411;
    const ccomplex_t IT_0467 = conj(IT_0305)*IT_0432;
    const ccomplex_t IT_0468 = IT_0240*conj(IT_0354) + IT_0236*IT_0414 +
       IT_0444 + IT_0445 + IT_0446 + IT_0447 + IT_0448 + IT_0449 + IT_0450 +
       IT_0451 + IT_0452 + IT_0456 + IT_0458 + conj(IT_0321)*IT_0459 + conj
      (IT_0315)*IT_0460 + IT_0461 + IT_0462 + IT_0464 + IT_0465 + IT_0466 +
       IT_0467;
    const ccomplex_t IT_0469 = conj(IT_0079)*IT_0081;
    const ccomplex_t IT_0470 = conj(IT_0115)*IT_0119;
    const ccomplex_t IT_0471 = conj(IT_0137)*IT_0140;
    const ccomplex_t IT_0472 = conj(IT_0164)*IT_0166;
    const ccomplex_t IT_0473 = IT_0249*conj(IT_0321);
    const ccomplex_t IT_0474 = IT_0230*conj(IT_0315);
    const ccomplex_t IT_0475 = IT_0240*conj(IT_0245);
    const ccomplex_t IT_0476 = conj(IT_0168)*IT_0169;
    const ccomplex_t IT_0477 = IT_0243*conj(IT_0351);
    const ccomplex_t IT_0478 = 24*IT_0269;
    const ccomplex_t IT_0479 = 24*IT_0165;
    const ccomplex_t IT_0480 = conj(IT_0210)*IT_0479;
    const ccomplex_t IT_0481 = 24*IT_0248;
    const ccomplex_t IT_0482 = conj(IT_0228)*IT_0460;
    const ccomplex_t IT_0483 = IT_0229*IT_0357;
    const ccomplex_t IT_0484 = 48*IT_0269;
    const ccomplex_t IT_0485 = IT_0373*IT_0484;
    const ccomplex_t IT_0486 = IT_0240*conj(IT_0260) + IT_0080*IT_0360 +
       IT_0442 + IT_0443 + IT_0453 + IT_0454 + IT_0469 + IT_0470 + IT_0471 +
       IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 + conj(IT_0203)
      *IT_0478 + IT_0480 + conj(IT_0220)*IT_0481 + IT_0482 + IT_0483 + IT_0485;
    const ccomplex_t IT_0487 = IT_0240*conj(IT_0260) + IT_0080*IT_0360 +
       IT_0464 + IT_0465 + IT_0466 + IT_0467 + IT_0469 + IT_0470 + IT_0471 +
       IT_0472 + IT_0473 + IT_0474 + IT_0475 + IT_0476 + IT_0477 + conj(IT_0203)
      *IT_0478 + IT_0480 + conj(IT_0220)*IT_0481 + IT_0482 + IT_0483 + IT_0485;
    const ccomplex_t IT_0488 = conj(IT_0349)*IT_0404;
    const ccomplex_t IT_0489 = conj(IT_0305)*IT_0426;
    const ccomplex_t IT_0490 = IT_0081*conj(IT_0228);
    const ccomplex_t IT_0491 = IT_0119*conj(IT_0351);
    const ccomplex_t IT_0492 = IT_0140*conj(IT_0168);
    const ccomplex_t IT_0493 = IT_0212*conj(IT_0220);
    const ccomplex_t IT_0494 = conj(IT_0210)*IT_0222;
    const ccomplex_t IT_0495 = conj(IT_0203)*IT_0230;
    const ccomplex_t IT_0496 = conj(IT_0137)*IT_0169;
    const ccomplex_t IT_0497 = conj(IT_0115)*IT_0243;
    const ccomplex_t IT_0498 = conj(IT_0336)*IT_0406;
    const ccomplex_t IT_0499 = conj(IT_0347)*IT_0437;
    const ccomplex_t IT_0500 = conj(IT_0164)*IT_0459;
    const ccomplex_t IT_0501 = conj(IT_0079)*IT_0460;
    const ccomplex_t IT_0502 = IT_0080*IT_0263;
    const ccomplex_t IT_0503 = IT_0229*IT_0268;
    const ccomplex_t IT_0504 = IT_0236*IT_0413;
    const ccomplex_t IT_0505 = conj(IT_0315)*IT_0455 + conj(IT_0321)*IT_0457 +
       IT_0488 + IT_0489 + IT_0490 + IT_0491 + IT_0492 + IT_0493 + IT_0494 +
       IT_0495 + IT_0496 + IT_0497 + IT_0498 + IT_0499 + IT_0500 + IT_0501 +
       IT_0502 + IT_0503 + IT_0504;
    const ccomplex_t IT_0506 = conj(IT_0347)*IT_0409;
    const ccomplex_t IT_0507 = conj(IT_0336)*IT_0436;
    const ccomplex_t IT_0508 = conj(IT_0305)*IT_0411;
    const ccomplex_t IT_0509 = conj(IT_0349)*IT_0432;
    const ccomplex_t IT_0510 = IT_0240*IT_0414 + conj(IT_0315)*IT_0455 + conj
      (IT_0321)*IT_0457 + IT_0490 + IT_0491 + IT_0492 + IT_0493 + IT_0494 +
       IT_0495 + IT_0496 + IT_0497 + IT_0500 + IT_0501 + IT_0502 + IT_0503 +
       IT_0504 + IT_0506 + IT_0507 + IT_0508 + IT_0509;
    const ccomplex_t IT_0511 = IT_0119*conj(IT_0137);
    const ccomplex_t IT_0512 = conj(IT_0115)*IT_0140;
    const ccomplex_t IT_0513 = IT_0166*conj(IT_0321);
    const ccomplex_t IT_0514 = conj(IT_0164)*IT_0249;
    const ccomplex_t IT_0515 = conj(IT_0079)*IT_0230;
    const ccomplex_t IT_0516 = conj(IT_0168)*IT_0243;
    const ccomplex_t IT_0517 = conj(IT_0228)*IT_0478;
    const ccomplex_t IT_0518 = conj(IT_0210)*IT_0481;
    const ccomplex_t IT_0519 = conj(IT_0203)*IT_0460;
    const ccomplex_t IT_0520 = IT_0080*IT_0357;
    const ccomplex_t IT_0521 = IT_0169*conj(IT_0351);
    const ccomplex_t IT_0522 = IT_0240*IT_0373;
    const ccomplex_t IT_0523 = IT_0372*IT_0484;
    const ccomplex_t IT_0524 = IT_0081*conj(IT_0315);
    const ccomplex_t IT_0525 = IT_0229*IT_0360 + conj(IT_0220)*IT_0479 +
       IT_0488 + IT_0489 + IT_0498 + IT_0499 + IT_0511 + IT_0512 + IT_0513 +
       IT_0514 + IT_0515 + IT_0516 + IT_0517 + IT_0518 + IT_0519 + IT_0520 +
       IT_0521 + IT_0522 + IT_0523 + IT_0524;
    const ccomplex_t IT_0526 = IT_0229*IT_0360 + conj(IT_0220)*IT_0479 +
       IT_0506 + IT_0507 + IT_0508 + IT_0509 + IT_0511 + IT_0512 + IT_0513 +
       IT_0514 + IT_0515 + IT_0516 + IT_0517 + IT_0518 + IT_0519 + IT_0520 +
       IT_0521 + IT_0522 + IT_0523 + IT_0524;
    const ccomplex_t IT_0527 = conj(IT_0021)*(IT_0082 + IT_0120 + IT_0141 +
       IT_0167 + IT_0170) + conj(IT_0178)*(IT_0204 + IT_0205 + IT_0213 + IT_0223
       + IT_0231 + IT_0237 + IT_0241 + IT_0242 + IT_0244) + conj(IT_0245)*
      (IT_0246 + IT_0247 + IT_0250 + IT_0251 + IT_0252) + conj(IT_0238)*(IT_0253
       + IT_0254 + IT_0255 + IT_0256 + IT_0257 + IT_0258 + IT_0259) + IT_0210
      *IT_0353 + IT_0164*IT_0370 + IT_0220*IT_0371 + IT_0321*IT_0374 + IT_0079
      *IT_0393 + IT_0315*IT_0394 + IT_0203*IT_0401 + IT_0228*IT_0402 + conj
      (IT_0266)*(IT_0082 + IT_0120 + IT_0141 + IT_0167 + IT_0170 + IT_0405 +
       IT_0407) + conj(IT_0234)*(IT_0204 + IT_0205 + IT_0213 + IT_0223 + IT_0231
       + IT_0237 + IT_0241 + IT_0242 + IT_0244 + IT_0405 + IT_0407) + conj
      (IT_0260)*(IT_0246 + IT_0247 + IT_0250 + IT_0251 + IT_0252 + IT_0410 +
       IT_0412) + conj(IT_0354)*(IT_0253 + IT_0254 + IT_0255 + IT_0256 + IT_0257
       + IT_0258 + IT_0259 + IT_0410 + IT_0412) + IT_0351*IT_0420 + IT_0137
      *IT_0421 + IT_0168*IT_0422 + IT_0115*IT_0423 + IT_0305*IT_0435 + IT_0347
      *IT_0438 + IT_0349*IT_0440 + IT_0336*IT_0441 + IT_0234*IT_0463 + IT_0178
      *IT_0468 + IT_0266*IT_0486 + IT_0021*IT_0487 + IT_0238*IT_0505 + IT_0354
      *IT_0510 + IT_0245*IT_0525 + IT_0260*IT_0526;
    return create_ccomplex_return(IT_0527);
}

