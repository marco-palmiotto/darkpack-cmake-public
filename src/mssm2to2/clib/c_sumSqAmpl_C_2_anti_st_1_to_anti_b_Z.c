#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_st_1_to_anti_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_st_1_to_anti_b_Z(
        cparam_t const *param
        )
{
    clearcache();
    const creal_t m_W = param->m_W;
    const creal_t m_Z = param->m_Z;
    const creal_t m_b = param->m_b;
    const creal_t m_t = param->m_t;
    const creal_t beta = param->beta;
    const creal_t e_em = param->e_em;
    const creal_t s_12 = param->s_12;
    const creal_t s_13 = param->s_13;
    const creal_t s_14 = param->s_14;
    const creal_t s_23 = param->s_23;
    const creal_t s_34 = param->s_34;
    const creal_t m_C_1 = param->m_C_1;
    const creal_t m_C_2 = param->m_C_2;
    const creal_t m_st_1 = param->m_st_1;
    const creal_t m_st_2 = param->m_st_2;
    const creal_t theta_W = param->theta_W;
    const creal_t Gamma_c1 = param->Gamma_c1;
    const creal_t Gamma_c2 = param->Gamma_c2;
    const creal_t Gamma_t1 = param->Gamma_t1;
    const creal_t Gamma_t2 = param->Gamma_t2;
    const creal_t reg_prop = param->reg_prop;
    const ccomplex_t U_d1 = param->U_d1;
    const ccomplex_t U_d2 = param->U_d2;
    const ccomplex_t V_tb = param->V_tb;
    const ccomplex_t V_u1 = param->V_u1;
    const ccomplex_t V_u2 = param->V_u2;
    const ccomplex_t U_Wm1 = param->U_Wm1;
    const ccomplex_t U_Wm2 = param->U_Wm2;
    const ccomplex_t V_Wp1 = param->V_Wp1;
    const ccomplex_t V_Wp2 = param->V_Wp2;
    const ccomplex_t U_st_00 = param->U_st_00;
    const ccomplex_t U_st_01 = param->U_st_01;
    const ccomplex_t U_st_10 = param->U_st_10;
    const ccomplex_t U_st_11 = param->U_st_11;
    const ccomplex_t IT_0000 = pow(m_Z, -2);
    const ccomplex_t IT_0001 = pow(m_Z, 2);
    const ccomplex_t IT_0002 = s_13*IT_0001;
    const ccomplex_t IT_0003 = s_14*s_34;
    const ccomplex_t IT_0004 = (-2)*IT_0003;
    const ccomplex_t IT_0005 = IT_0002 + IT_0004;
    const ccomplex_t IT_0006 = IT_0000*IT_0005;
    const ccomplex_t IT_0007 = (-6)*IT_0006;
    const ccomplex_t IT_0008 = 12*s_13;
    const ccomplex_t IT_0009 = IT_0007 + IT_0008;
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0011 = pow(m_W, -1);
    const ccomplex_t IT_0012 = cos(beta);
    const ccomplex_t IT_0013 = cpow(IT_0012, -1);
    const ccomplex_t IT_0014 = sin(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0011*IT_0013*IT_0015*conj(U_st_00);
    const ccomplex_t IT_0017 = (-0.5)*IT_0016;
    const ccomplex_t IT_0018 = cos(theta_W);
    const ccomplex_t IT_0019 = cpow(IT_0018, -1);
    const ccomplex_t IT_0020 = IT_0014*IT_0019;
    const ccomplex_t IT_0021 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0022 = IT_0020*IT_0021;
    const ccomplex_t IT_0023 = IT_0015*IT_0018;
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0026 = IT_0023*IT_0025;
    const ccomplex_t IT_0027 = (0 + _Complex_I*1)*(IT_0022 + -IT_0024 + (-2)
      *IT_0026);
    const ccomplex_t IT_0028 = (-0.5)*IT_0027;
    const ccomplex_t IT_0029 = IT_0017*IT_0028;
    const ccomplex_t IT_0030 = pow(m_b, 2);
    const ccomplex_t IT_0031 = pow(m_st_1, 2);
    const ccomplex_t IT_0032 = cpow((-2)*s_23 + IT_0030 + IT_0031 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0033 = IT_0010*IT_0029*IT_0032;
    const ccomplex_t IT_0034 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0035 = IT_0020*IT_0034;
    const ccomplex_t IT_0036 = IT_0023*IT_0034;
    const ccomplex_t IT_0037 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0038 = IT_0023*IT_0037;
    const ccomplex_t IT_0039 = (0 + _Complex_I*1)*(IT_0035 + -IT_0036 + (-2)
      *IT_0038);
    const ccomplex_t IT_0040 = 0.5*IT_0039;
    const ccomplex_t IT_0041 = IT_0017*IT_0040;
    const ccomplex_t IT_0042 = IT_0032*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1)*IT_0042;
    const ccomplex_t IT_0044 = m_C_2*IT_0043;
    const ccomplex_t IT_0045 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0011*IT_0013*IT_0015*conj(U_st_00);
    const ccomplex_t IT_0046 = (-0.5)*IT_0045;
    const ccomplex_t IT_0047 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0048 = IT_0020*IT_0047;
    const ccomplex_t IT_0049 = IT_0023*IT_0047;
    const ccomplex_t IT_0050 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0051 = IT_0023*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + -IT_0049 + (-2)
      *IT_0051);
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0046*IT_0053;
    const ccomplex_t IT_0055 = cpow((-2)*s_23 + IT_0030 + IT_0031 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0056 = IT_0054*IT_0055;
    const ccomplex_t IT_0057 = (0 + _Complex_I*1)*IT_0056;
    const ccomplex_t IT_0058 = m_C_2*IT_0057;
    const ccomplex_t IT_0059 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0060 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0061 = IT_0020*IT_0060;
    const ccomplex_t IT_0062 = IT_0023*IT_0060;
    const ccomplex_t IT_0063 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0064 = IT_0023*IT_0063;
    const ccomplex_t IT_0065 = (0 + _Complex_I*1)*(IT_0061 + -IT_0062 + (-2)
      *IT_0064);
    const ccomplex_t IT_0066 = (-0.5)*IT_0065;
    const ccomplex_t IT_0067 = IT_0046*IT_0066;
    const ccomplex_t IT_0068 = IT_0055*IT_0059*IT_0067;
    const ccomplex_t IT_0069 = IT_0033 + IT_0044 + IT_0058 + IT_0068;
    const ccomplex_t IT_0070 = e_em*IT_0020;
    const ccomplex_t IT_0071 = e_em*IT_0023;
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*(IT_0070 + 3*IT_0071);
    const ccomplex_t IT_0073 = (-0.166666666666667)*IT_0072;
    const ccomplex_t IT_0074 = pow(m_C_2, 2);
    const ccomplex_t IT_0075 = cpow((-2)*s_12 + IT_0030 + -IT_0031 + -IT_0074 
      + -reg_prop, -1);
    const ccomplex_t IT_0076 = IT_0073*IT_0075;
    const ccomplex_t IT_0077 = (0 + _Complex_I*1)*IT_0076;
    const ccomplex_t IT_0078 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0079 = IT_0015*IT_0078;
    const ccomplex_t IT_0080 = (0 + _Complex_I*1)*IT_0079;
    const ccomplex_t IT_0081 = -IT_0080;
    const ccomplex_t IT_0082 = m_b*IT_0081;
    const ccomplex_t IT_0083 = sin(beta);
    const ccomplex_t IT_0084 = cpow(IT_0083, -1);
    const ccomplex_t IT_0085 = IT_0015*IT_0084;
    const ccomplex_t IT_0086 = m_t*V_tb*V_u2*e_em*IT_0011*conj(U_st_10);
    const ccomplex_t IT_0087 = IT_0085*IT_0086;
    const ccomplex_t IT_0088 = 1.4142135623731*IT_0087;
    const ccomplex_t IT_0089 = (0 + _Complex_I*1)*IT_0088;
    const ccomplex_t IT_0090 = (-0.5)*IT_0089;
    const ccomplex_t IT_0091 = m_b*IT_0090;
    const ccomplex_t IT_0092 = -IT_0091;
    const ccomplex_t IT_0093 = IT_0082 + IT_0092;
    const ccomplex_t IT_0094 = IT_0077*IT_0093;
    const ccomplex_t IT_0095 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*e_em*IT_0014*IT_0019;
    const ccomplex_t IT_0097 = 0.333333333333333*IT_0096;
    const ccomplex_t IT_0098 = -IT_0081;
    const ccomplex_t IT_0099 = IT_0090 + IT_0098;
    const ccomplex_t IT_0100 = IT_0097*IT_0099;
    const ccomplex_t IT_0101 = IT_0075*IT_0095*IT_0100;
    const ccomplex_t IT_0102 = IT_0094 + IT_0101;
    const ccomplex_t IT_0103 = IT_0069 + IT_0102;
    const ccomplex_t IT_0104 = conj(IT_0069) + conj(IT_0102);
    const ccomplex_t IT_0105 = m_t*V_tb*V_u1*e_em*IT_0011*conj(U_st_10);
    const ccomplex_t IT_0106 = IT_0085*IT_0105;
    const ccomplex_t IT_0107 = 1.4142135623731*IT_0106;
    const ccomplex_t IT_0108 = (0 + _Complex_I*1)*IT_0107;
    const ccomplex_t IT_0109 = (-0.5)*IT_0108;
    const ccomplex_t IT_0110 = V_tb*e_em*V_Wp1*conj(U_st_00);
    const ccomplex_t IT_0111 = IT_0015*IT_0110;
    const ccomplex_t IT_0112 = (0 + _Complex_I*1)*IT_0111;
    const ccomplex_t IT_0113 = -IT_0112;
    const ccomplex_t IT_0114 = -IT_0113;
    const ccomplex_t IT_0115 = IT_0109 + IT_0114;
    const ccomplex_t IT_0116 = IT_0040*IT_0115;
    const ccomplex_t IT_0117 = IT_0010*IT_0032*IT_0116;
    const ccomplex_t IT_0118 = IT_0053*IT_0099;
    const ccomplex_t IT_0119 = IT_0055*IT_0059*IT_0118;
    const ccomplex_t IT_0120 = IT_0055*IT_0066;
    const ccomplex_t IT_0121 = (0 + _Complex_I*1)*IT_0120;
    const ccomplex_t IT_0122 = m_C_2*IT_0081;
    const ccomplex_t IT_0123 = m_C_2*IT_0090;
    const ccomplex_t IT_0124 = -IT_0123;
    const ccomplex_t IT_0125 = IT_0122 + IT_0124;
    const ccomplex_t IT_0126 = IT_0121*IT_0125;
    const ccomplex_t IT_0127 = IT_0028*IT_0032;
    const ccomplex_t IT_0128 = (0 + _Complex_I*1)*IT_0127;
    const ccomplex_t IT_0129 = m_C_2*IT_0113;
    const ccomplex_t IT_0130 = m_C_2*IT_0109;
    const ccomplex_t IT_0131 = -IT_0130;
    const ccomplex_t IT_0132 = IT_0129 + IT_0131;
    const ccomplex_t IT_0133 = IT_0128*IT_0132;
    const ccomplex_t IT_0134 = -IT_0117 + -IT_0119 + IT_0126 + IT_0133;
    const ccomplex_t IT_0135 = IT_0046*IT_0073;
    const ccomplex_t IT_0136 = IT_0075*IT_0095*IT_0135;
    const ccomplex_t IT_0137 = IT_0046*IT_0097;
    const ccomplex_t IT_0138 = IT_0075*IT_0137;
    const ccomplex_t IT_0139 = (0 + _Complex_I*1)*IT_0138;
    const ccomplex_t IT_0140 = m_b*IT_0139;
    const ccomplex_t IT_0141 = IT_0136 + -IT_0140;
    const ccomplex_t IT_0142 = IT_0134 + IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0134) + conj(IT_0141);
    const ccomplex_t IT_0144 = m_b*m_C_2;
    const ccomplex_t IT_0145 = pow(s_14, 2);
    const ccomplex_t IT_0146 = IT_0000*IT_0145;
    const ccomplex_t IT_0147 = -IT_0074;
    const ccomplex_t IT_0148 = IT_0146 + IT_0147;
    const ccomplex_t IT_0149 = IT_0144*IT_0148;
    const ccomplex_t IT_0150 = (-2)*IT_0043 + (-2)*IT_0057;
    const ccomplex_t IT_0151 = IT_0099*IT_0121;
    const ccomplex_t IT_0152 = IT_0115*IT_0128;
    const ccomplex_t IT_0153 = 2*IT_0151 + 2*IT_0152;
    const ccomplex_t IT_0154 = 6*IT_0153;
    const ccomplex_t IT_0155 = cpow((-2)*s_13 + IT_0030 + IT_0074 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*IT_0155;
    const ccomplex_t IT_0157 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0015
      *IT_0018 + -1./3*IT_0014*IT_0019)*U_st_00*conj(U_st_00) + IT_0014*IT_0019
      *U_st_10*conj(U_st_10));
    const ccomplex_t IT_0158 = 1.33333333333333*IT_0157;
    const ccomplex_t IT_0159 = IT_0081*IT_0158;
    const ccomplex_t IT_0160 = IT_0090*IT_0158;
    const ccomplex_t IT_0161 = -IT_0160;
    const ccomplex_t IT_0162 = IT_0159 + IT_0161;
    const ccomplex_t IT_0163 = IT_0156*IT_0162;
    const ccomplex_t IT_0164 = cpow((-2)*s_13 + IT_0030 + IT_0074 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0165 = (0 + _Complex_I*1)*IT_0164;
    const ccomplex_t IT_0166 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0167 = IT_0015*IT_0166;
    const ccomplex_t IT_0168 = (0 + _Complex_I*1)*IT_0167;
    const ccomplex_t IT_0169 = -IT_0168;
    const ccomplex_t IT_0170 = (0 + _Complex_I*1)*e_em*((-0.75)*(IT_0015
      *IT_0018 + -1./3*IT_0014*IT_0019)*conj(U_st_00)*U_st_01 + IT_0014*IT_0019
      *conj(U_st_10)*U_st_11);
    const ccomplex_t IT_0171 = 1.33333333333333*IT_0170;
    const ccomplex_t IT_0172 = IT_0169*IT_0171;
    const ccomplex_t IT_0173 = m_t*V_tb*V_u2*e_em*IT_0011*conj(U_st_11);
    const ccomplex_t IT_0174 = IT_0085*IT_0173;
    const ccomplex_t IT_0175 = 1.4142135623731*IT_0174;
    const ccomplex_t IT_0176 = (0 + _Complex_I*1)*IT_0175;
    const ccomplex_t IT_0177 = (-0.5)*IT_0176;
    const ccomplex_t IT_0178 = IT_0171*IT_0177;
    const ccomplex_t IT_0179 = -IT_0178;
    const ccomplex_t IT_0180 = IT_0172 + IT_0179;
    const ccomplex_t IT_0181 = IT_0165*IT_0180;
    const ccomplex_t IT_0182 = -IT_0163 + -IT_0181;
    const ccomplex_t IT_0183 = (-6)*IT_0182;
    const ccomplex_t IT_0184 = IT_0154 + IT_0183;
    const ccomplex_t IT_0185 = 6*conj(IT_0153);
    const ccomplex_t IT_0186 = (-6)*conj(IT_0182);
    const ccomplex_t IT_0187 = IT_0185 + IT_0186;
    const ccomplex_t IT_0188 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0011*IT_0013*IT_0015*conj(U_st_01);
    const ccomplex_t IT_0189 = (-0.5)*IT_0188;
    const ccomplex_t IT_0190 = IT_0164*IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0190;
    const ccomplex_t IT_0192 = IT_0171*IT_0191;
    const ccomplex_t IT_0193 = IT_0046*IT_0155;
    const ccomplex_t IT_0194 = (0 + _Complex_I*1)*IT_0193;
    const ccomplex_t IT_0195 = IT_0158*IT_0194;
    const ccomplex_t IT_0196 = IT_0192 + IT_0195;
    const ccomplex_t IT_0197 = 6*IT_0182;
    const ccomplex_t IT_0198 = (-6)*IT_0153;
    const ccomplex_t IT_0199 = IT_0197 + IT_0198;
    const ccomplex_t IT_0200 = 6*conj(IT_0182);
    const ccomplex_t IT_0201 = (-6)*conj(IT_0153);
    const ccomplex_t IT_0202 = IT_0200 + IT_0201;
    const ccomplex_t IT_0203 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0204 = -IT_0203;
    const ccomplex_t IT_0205 = s_13 + IT_0204;
    const ccomplex_t IT_0206 = s_13*IT_0205;
    const ccomplex_t IT_0207 = IT_0077*IT_0099;
    const ccomplex_t IT_0208 = IT_0163 + IT_0181 + (-2)*IT_0207;
    const ccomplex_t IT_0209 = (-2)*IT_0139 + -IT_0192 + -IT_0195;
    const ccomplex_t IT_0210 = 6*IT_0150;
    const ccomplex_t IT_0211 = (-6)*IT_0196;
    const ccomplex_t IT_0212 = IT_0210 + IT_0211;
    const ccomplex_t IT_0213 = 6*conj(IT_0150);
    const ccomplex_t IT_0214 = (-6)*conj(IT_0196);
    const ccomplex_t IT_0215 = IT_0213 + IT_0214;
    const ccomplex_t IT_0216 = s_13*s_14;
    const ccomplex_t IT_0217 = s_34*IT_0074;
    const ccomplex_t IT_0218 = -IT_0217;
    const ccomplex_t IT_0219 = IT_0216 + IT_0218;
    const ccomplex_t IT_0220 = -IT_0077*(IT_0081 + -IT_0090);
    const ccomplex_t IT_0221 = -IT_0220;
    const ccomplex_t IT_0222 = -IT_0151 + -IT_0152;
    const ccomplex_t IT_0223 = IT_0043 + IT_0057;
    const ccomplex_t IT_0224 = -IT_0139;
    const ccomplex_t IT_0225 = m_b*IT_0148;
    const ccomplex_t IT_0226 = m_b*IT_0205;
    const ccomplex_t IT_0227 = 6*IT_0208;
    const ccomplex_t IT_0228 = 6*conj(IT_0208);
    const ccomplex_t IT_0229 = 6*IT_0209;
    const ccomplex_t IT_0230 = 6*conj(IT_0209);
    const ccomplex_t IT_0231 = pow(s_34, 2);
    const ccomplex_t IT_0232 = IT_0000*IT_0231;
    const ccomplex_t IT_0233 = -IT_0232;
    const ccomplex_t IT_0234 = IT_0030 + IT_0233;
    const ccomplex_t IT_0235 = m_C_2*IT_0234;
    const ccomplex_t IT_0236 = 6*IT_0069;
    const ccomplex_t IT_0237 = 6*conj(IT_0069);
    const ccomplex_t IT_0238 = (-18)*conj(IT_0069);
    const ccomplex_t IT_0239 = s_13*IT_0148;
    const ccomplex_t IT_0240 = 6*IT_0196;
    const ccomplex_t IT_0241 = (-6)*IT_0150;
    const ccomplex_t IT_0242 = IT_0240 + IT_0241;
    const ccomplex_t IT_0243 = (-6)*conj(IT_0150);
    const ccomplex_t IT_0244 = IT_0144*IT_0205;
    const ccomplex_t IT_0245 = (-6)*IT_0208;
    const ccomplex_t IT_0246 = (-6)*conj(IT_0208);
    const ccomplex_t IT_0247 = IT_0144*IT_0234;
    const ccomplex_t IT_0248 = s_14*IT_0001;
    const ccomplex_t IT_0249 = IT_0000*IT_0248;
    const ccomplex_t IT_0250 = 0.5*IT_0249;
    const ccomplex_t IT_0251 = s_14 + IT_0250;
    const ccomplex_t IT_0252 = m_b*IT_0251;
    const ccomplex_t IT_0253 = 2*IT_0252;
    const ccomplex_t IT_0254 = 6*conj(IT_0102);
    const ccomplex_t IT_0255 = IT_0237 + IT_0254;
    const ccomplex_t IT_0256 = s_34*IT_0001;
    const ccomplex_t IT_0257 = IT_0000*IT_0256;
    const ccomplex_t IT_0258 = (-0.25)*IT_0257;
    const ccomplex_t IT_0259 = s_34 + IT_0258;
    const ccomplex_t IT_0260 = m_C_2*IT_0259;
    const ccomplex_t IT_0261 = (-4)*IT_0260;
    const ccomplex_t IT_0262 = 6*conj(IT_0134);
    const ccomplex_t IT_0263 = 6*conj(IT_0141);
    const ccomplex_t IT_0264 = IT_0262 + IT_0263;
    const ccomplex_t IT_0265 = s_14*IT_0030;
    const ccomplex_t IT_0266 = s_13*s_34;
    const ccomplex_t IT_0267 = -IT_0266;
    const ccomplex_t IT_0268 = IT_0265 + IT_0267;
    const ccomplex_t IT_0269 = (-0.25)*IT_0249;
    const ccomplex_t IT_0270 = s_14 + IT_0269;
    const ccomplex_t IT_0271 = m_b*IT_0270;
    const ccomplex_t IT_0272 = (-4)*IT_0271;
    const ccomplex_t IT_0273 = (-6)*conj(IT_0069);
    const ccomplex_t IT_0274 = (-6)*conj(IT_0102);
    const ccomplex_t IT_0275 = IT_0273 + IT_0274;
    const ccomplex_t IT_0276 = 0.5*IT_0257;
    const ccomplex_t IT_0277 = s_34 + IT_0276;
    const ccomplex_t IT_0278 = m_C_2*IT_0277;
    const ccomplex_t IT_0279 = 2*IT_0278;
    const ccomplex_t IT_0280 = (-6)*conj(IT_0134);
    const ccomplex_t IT_0281 = (-6)*conj(IT_0141);
    const ccomplex_t IT_0282 = IT_0280 + IT_0281;
    const ccomplex_t IT_0283 = s_13*IT_0234;
    const ccomplex_t IT_0284 = (-6)*IT_0209;
    const ccomplex_t IT_0285 = m_C_2*IT_0205;
    const ccomplex_t IT_0286 = (-6)*IT_0069;
    const ccomplex_t IT_0287 = IT_0001*IT_0144;
    const ccomplex_t IT_0288 = (-3)*IT_0287;
    const ccomplex_t IT_0289 = 3*IT_0287;
    const ccomplex_t IT_0290 = 24*IT_0003;
    const ccomplex_t IT_0291 = pow(m_Z, 4);
    const ccomplex_t IT_0292 = s_13*IT_0291;
    const ccomplex_t IT_0293 = IT_0000*IT_0292;
    const ccomplex_t IT_0294 = (-6)*IT_0293;
    const ccomplex_t IT_0295 = IT_0290 + IT_0294;
    const ccomplex_t IT_0296 = 6*IT_0293;
    const ccomplex_t IT_0297 = (-12)*IT_0002;
    const ccomplex_t IT_0298 = IT_0290 + IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = 0.166666666666667*IT_0222;
    const ccomplex_t IT_0300 = 6*IT_0102;
    const ccomplex_t IT_0301 = IT_0236 + IT_0300;
    const ccomplex_t IT_0302 = 0.166666666666667*IT_0253;
    const ccomplex_t IT_0303 = 6*IT_0134;
    const ccomplex_t IT_0304 = 6*IT_0141;
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0261;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0268;
    const ccomplex_t IT_0308 = (-6)*IT_0102;
    const ccomplex_t IT_0309 = IT_0286 + IT_0308;
    const ccomplex_t IT_0310 = 0.166666666666667*IT_0272;
    const ccomplex_t IT_0311 = (-6)*IT_0134;
    const ccomplex_t IT_0312 = (-6)*IT_0141;
    const ccomplex_t IT_0313 = IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = 0.166666666666667*IT_0279;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0223;
    const ccomplex_t IT_0316 = 0.166666666666667*IT_0224;
    const ccomplex_t IT_0317 = IT_0284*IT_0307;
    const ccomplex_t IT_0318 = (-6)*conj(IT_0209);
    const ccomplex_t IT_0319 = IT_0268*IT_0318;
    const ccomplex_t IT_0320 = IT_0009*(IT_0103*IT_0104 + IT_0142*IT_0143) +
       IT_0149*(conj(IT_0150)*IT_0184 + IT_0150*IT_0187 + conj(IT_0196)*IT_0199 
      + IT_0196*IT_0202) + IT_0206*(IT_0187*IT_0208 + IT_0184*conj(IT_0208) +
       conj(IT_0209)*IT_0212 + IT_0209*IT_0215) + IT_0219*(IT_0187*(IT_0221 +
       IT_0222) + IT_0184*(conj(IT_0221) + conj(IT_0222)) + IT_0215*(IT_0223 +
       IT_0224) + IT_0212*(conj(IT_0223) + conj(IT_0224))) + (IT_0104*IT_0184 +
       IT_0103*IT_0187 + IT_0143*IT_0212 + IT_0142*IT_0215)*IT_0225 + IT_0226*
      (IT_0104*IT_0227 + IT_0103*IT_0228 + IT_0143*IT_0229 + IT_0142*IT_0230) +
       IT_0235*(IT_0143*IT_0227 + IT_0142*IT_0228 + conj(IT_0102)*IT_0229 +
       IT_0102*IT_0230 + conj(IT_0209)*IT_0236 + IT_0209*IT_0237) + IT_0144*
      (IT_0103*((-18)*conj(IT_0134) + (-18)*conj(IT_0141)) + IT_0142*((-18)*conj
      (IT_0102) + IT_0238)) + IT_0239*(conj(IT_0153)*IT_0154 + conj(IT_0182)
      *IT_0199 + IT_0182*IT_0201 + conj(IT_0150)*IT_0210 + conj(IT_0196)*IT_0242
       + IT_0196*IT_0243) + IT_0244*(IT_0187*IT_0209 + IT_0184*conj(IT_0209) +
       conj(IT_0150)*IT_0227 + IT_0150*IT_0228 + conj(IT_0196)*IT_0245 + IT_0196
      *IT_0246) + (conj(IT_0209)*IT_0245 + IT_0209*IT_0246)*IT_0247 + IT_0222*
      (IT_0253*IT_0255 + IT_0261*IT_0264 + IT_0246*IT_0268) + IT_0221*(IT_0246
      *IT_0268 + IT_0272*IT_0275 + IT_0279*IT_0282) + IT_0283*(conj(IT_0208)
      *IT_0245 + conj(IT_0209)*IT_0284) + IT_0285*(IT_0143*IT_0199 + IT_0142
      *IT_0202 + 6*conj(IT_0196)*(IT_0102 + 0.166666666666667*IT_0236) + IT_0196
      *IT_0237 + conj(IT_0102)*IT_0242 + IT_0102*IT_0243 + IT_0150*IT_0273 +
       conj(IT_0150)*IT_0286) + 6*conj(IT_0222)*(IT_0223*IT_0288 + -IT_0224
      *IT_0289 + 0.166666666666667*IT_0221*IT_0295 + IT_0298*IT_0299 + IT_0301
      *IT_0302 + IT_0305*IT_0306 + IT_0245*IT_0307) + 6*conj(IT_0221)*(IT_0224
      *IT_0288 + -IT_0223*IT_0289 + 0.166666666666667*IT_0221*IT_0298 + IT_0295
      *IT_0299 + IT_0245*IT_0307 + IT_0309*IT_0310 + IT_0313*IT_0314) + 6*conj
      (IT_0223)*(IT_0222*IT_0288 + -IT_0221*IT_0289 + IT_0302*IT_0305 + IT_0301
      *IT_0306 + IT_0298*IT_0315 + IT_0295*IT_0316 + IT_0317) + 6*conj(IT_0224)*
      (IT_0221*IT_0288 + -IT_0222*IT_0289 + IT_0310*IT_0313 + IT_0309*IT_0314 +
       IT_0295*IT_0315 + IT_0298*IT_0316 + IT_0317) + IT_0223*(IT_0255*IT_0261 +
       IT_0253*IT_0264 + IT_0319) + IT_0224*(IT_0275*IT_0279 + IT_0272*IT_0282 +
       IT_0319);
    return create_ccomplex_return(IT_0320);
}

