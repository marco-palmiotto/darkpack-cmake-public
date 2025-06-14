#include "clooptools.h"
#include "marty/core/looptools_init.h"
#include <math.h>
#include "stdbool.h"
#include "c_sumSqAmpl_C_2_anti_st_2_to_anti_b_Z.h"
#include "ccommon.h"

#include "cparams.h"

ccomplex_return_t c_sumSqAmpl_C_2_anti_st_2_to_anti_b_Z(
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
    const ccomplex_t IT_0010 = (0 + _Complex_I*1)*m_C_2;
    const ccomplex_t IT_0011 = pow(m_b, 2);
    const ccomplex_t IT_0012 = pow(m_st_2, 2);
    const ccomplex_t IT_0013 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_C_2*
      (m_C_2 + (0 + _Complex_I*-1)*Gamma_c2) + reg_prop, -1);
    const ccomplex_t IT_0014 = cos(theta_W);
    const ccomplex_t IT_0015 = cpow(IT_0014, -1);
    const ccomplex_t IT_0016 = sin(theta_W);
    const ccomplex_t IT_0017 = IT_0015*IT_0016;
    const ccomplex_t IT_0018 = V_u2*conj(V_u2)*e_em;
    const ccomplex_t IT_0019 = IT_0017*IT_0018;
    const ccomplex_t IT_0020 = cpow(IT_0016, -1);
    const ccomplex_t IT_0021 = IT_0014*IT_0020;
    const ccomplex_t IT_0022 = IT_0018*IT_0021;
    const ccomplex_t IT_0023 = e_em*V_Wp2*conj(V_Wp2);
    const ccomplex_t IT_0024 = IT_0021*IT_0023;
    const ccomplex_t IT_0025 = (0 + _Complex_I*1)*(IT_0019 + -IT_0022 + (-2)
      *IT_0024);
    const ccomplex_t IT_0026 = (-0.5)*IT_0025;
    const ccomplex_t IT_0027 = pow(m_W, -1);
    const ccomplex_t IT_0028 = cos(beta);
    const ccomplex_t IT_0029 = cpow(IT_0028, -1);
    const ccomplex_t IT_0030 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0020*IT_0027*IT_0029*conj(U_st_01);
    const ccomplex_t IT_0031 = (-0.5)*IT_0030;
    const ccomplex_t IT_0032 = IT_0026*IT_0031;
    const ccomplex_t IT_0033 = IT_0010*IT_0013*IT_0032;
    const ccomplex_t IT_0034 = (0 + _Complex_I*1)*m_C_1;
    const ccomplex_t IT_0035 = cpow((-2)*s_23 + IT_0011 + IT_0012 + -m_C_1*
      (m_C_1 + (0 + _Complex_I*-1)*Gamma_c1) + reg_prop, -1);
    const ccomplex_t IT_0036 = conj(V_u1)*V_u2*e_em;
    const ccomplex_t IT_0037 = IT_0017*IT_0036;
    const ccomplex_t IT_0038 = IT_0021*IT_0036;
    const ccomplex_t IT_0039 = e_em*conj(V_Wp1)*V_Wp2;
    const ccomplex_t IT_0040 = IT_0021*IT_0039;
    const ccomplex_t IT_0041 = (0 + _Complex_I*1)*(IT_0037 + -IT_0038 + (-2)
      *IT_0040);
    const ccomplex_t IT_0042 = (-0.5)*IT_0041;
    const ccomplex_t IT_0043 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d1)
      *V_tb*e_em*IT_0020*IT_0027*IT_0029*conj(U_st_01);
    const ccomplex_t IT_0044 = (-0.5)*IT_0043;
    const ccomplex_t IT_0045 = IT_0042*IT_0044;
    const ccomplex_t IT_0046 = IT_0034*IT_0035*IT_0045;
    const ccomplex_t IT_0047 = U_d1*conj(U_d2)*e_em;
    const ccomplex_t IT_0048 = IT_0017*IT_0047;
    const ccomplex_t IT_0049 = IT_0021*IT_0047;
    const ccomplex_t IT_0050 = e_em*U_Wm1*conj(U_Wm2);
    const ccomplex_t IT_0051 = IT_0021*IT_0050;
    const ccomplex_t IT_0052 = (0 + _Complex_I*1)*(IT_0048 + -IT_0049 + (-2)
      *IT_0051);
    const ccomplex_t IT_0053 = 0.5*IT_0052;
    const ccomplex_t IT_0054 = IT_0044*IT_0053;
    const ccomplex_t IT_0055 = IT_0035*IT_0054;
    const ccomplex_t IT_0056 = (0 + _Complex_I*1)*IT_0055;
    const ccomplex_t IT_0057 = m_C_2*IT_0056;
    const ccomplex_t IT_0058 = U_d2*conj(U_d2)*e_em;
    const ccomplex_t IT_0059 = IT_0017*IT_0058;
    const ccomplex_t IT_0060 = IT_0021*IT_0058;
    const ccomplex_t IT_0061 = e_em*U_Wm2*conj(U_Wm2);
    const ccomplex_t IT_0062 = IT_0021*IT_0061;
    const ccomplex_t IT_0063 = (0 + _Complex_I*1)*(IT_0059 + -IT_0060 + (-2)
      *IT_0062);
    const ccomplex_t IT_0064 = 0.5*IT_0063;
    const ccomplex_t IT_0065 = IT_0031*IT_0064;
    const ccomplex_t IT_0066 = IT_0013*IT_0065;
    const ccomplex_t IT_0067 = (0 + _Complex_I*1)*IT_0066;
    const ccomplex_t IT_0068 = m_C_2*IT_0067;
    const ccomplex_t IT_0069 = IT_0033 + IT_0046 + IT_0057 + IT_0068;
    const ccomplex_t IT_0070 = pow(m_C_2, 2);
    const ccomplex_t IT_0071 = cpow((-2)*s_12 + IT_0011 + -IT_0012 + -IT_0070 
      + -reg_prop, -1);
    const ccomplex_t IT_0072 = (0 + _Complex_I*1)*m_b;
    const ccomplex_t IT_0073 = (0 + _Complex_I*1)*e_em*IT_0015*IT_0016;
    const ccomplex_t IT_0074 = 0.333333333333333*IT_0073;
    const ccomplex_t IT_0075 = sin(beta);
    const ccomplex_t IT_0076 = cpow(IT_0075, -1);
    const ccomplex_t IT_0077 = IT_0020*IT_0076;
    const ccomplex_t IT_0078 = m_t*V_tb*V_u2*e_em*IT_0027*conj(U_st_11);
    const ccomplex_t IT_0079 = IT_0077*IT_0078;
    const ccomplex_t IT_0080 = 1.4142135623731*IT_0079;
    const ccomplex_t IT_0081 = (0 + _Complex_I*1)*IT_0080;
    const ccomplex_t IT_0082 = (-0.5)*IT_0081;
    const ccomplex_t IT_0083 = V_tb*e_em*V_Wp2*conj(U_st_01);
    const ccomplex_t IT_0084 = IT_0020*IT_0083;
    const ccomplex_t IT_0085 = (0 + _Complex_I*1)*IT_0084;
    const ccomplex_t IT_0086 = -IT_0085;
    const ccomplex_t IT_0087 = -IT_0086;
    const ccomplex_t IT_0088 = IT_0082 + IT_0087;
    const ccomplex_t IT_0089 = IT_0074*IT_0088;
    const ccomplex_t IT_0090 = IT_0071*IT_0072*IT_0089;
    const ccomplex_t IT_0091 = e_em*IT_0017;
    const ccomplex_t IT_0092 = e_em*IT_0021;
    const ccomplex_t IT_0093 = (0 + _Complex_I*1)*(IT_0091 + 3*IT_0092);
    const ccomplex_t IT_0094 = (-0.166666666666667)*IT_0093;
    const ccomplex_t IT_0095 = IT_0071*IT_0094;
    const ccomplex_t IT_0096 = (0 + _Complex_I*1)*IT_0095;
    const ccomplex_t IT_0097 = m_b*IT_0082;
    const ccomplex_t IT_0098 = -IT_0097;
    const ccomplex_t IT_0099 = m_b*IT_0086;
    const ccomplex_t IT_0100 = IT_0098 + IT_0099;
    const ccomplex_t IT_0101 = IT_0096*IT_0100;
    const ccomplex_t IT_0102 = IT_0090 + IT_0101;
    const ccomplex_t IT_0103 = IT_0069 + IT_0102;
    const ccomplex_t IT_0104 = conj(IT_0069) + conj(IT_0102);
    const ccomplex_t IT_0105 = IT_0031*IT_0094;
    const ccomplex_t IT_0106 = IT_0071*IT_0072*IT_0105;
    const ccomplex_t IT_0107 = IT_0031*IT_0074;
    const ccomplex_t IT_0108 = IT_0071*IT_0107;
    const ccomplex_t IT_0109 = (0 + _Complex_I*1)*IT_0108;
    const ccomplex_t IT_0110 = m_b*IT_0109;
    const ccomplex_t IT_0111 = IT_0106 + -IT_0110;
    const ccomplex_t IT_0112 = m_t*V_tb*V_u1*e_em*IT_0027*conj(U_st_11);
    const ccomplex_t IT_0113 = IT_0077*IT_0112;
    const ccomplex_t IT_0114 = 1.4142135623731*IT_0113;
    const ccomplex_t IT_0115 = (0 + _Complex_I*1)*IT_0114;
    const ccomplex_t IT_0116 = (-0.5)*IT_0115;
    const ccomplex_t IT_0117 = V_tb*e_em*V_Wp1*conj(U_st_01);
    const ccomplex_t IT_0118 = IT_0020*IT_0117;
    const ccomplex_t IT_0119 = (0 + _Complex_I*1)*IT_0118;
    const ccomplex_t IT_0120 = -IT_0119;
    const ccomplex_t IT_0121 = -IT_0120;
    const ccomplex_t IT_0122 = IT_0116 + IT_0121;
    const ccomplex_t IT_0123 = IT_0053*IT_0122;
    const ccomplex_t IT_0124 = IT_0034*IT_0035*IT_0123;
    const ccomplex_t IT_0125 = IT_0035*IT_0042;
    const ccomplex_t IT_0126 = (0 + _Complex_I*1)*IT_0125;
    const ccomplex_t IT_0127 = m_C_2*IT_0120;
    const ccomplex_t IT_0128 = m_C_2*IT_0116;
    const ccomplex_t IT_0129 = -IT_0128;
    const ccomplex_t IT_0130 = IT_0127 + IT_0129;
    const ccomplex_t IT_0131 = IT_0126*IT_0130;
    const ccomplex_t IT_0132 = IT_0064*IT_0088;
    const ccomplex_t IT_0133 = IT_0010*IT_0013*IT_0132;
    const ccomplex_t IT_0134 = IT_0013*IT_0026;
    const ccomplex_t IT_0135 = (0 + _Complex_I*1)*IT_0134;
    const ccomplex_t IT_0136 = m_C_2*IT_0086;
    const ccomplex_t IT_0137 = m_C_2*IT_0082;
    const ccomplex_t IT_0138 = -IT_0137;
    const ccomplex_t IT_0139 = IT_0136 + IT_0138;
    const ccomplex_t IT_0140 = IT_0135*IT_0139;
    const ccomplex_t IT_0141 = -IT_0124 + IT_0131 + -IT_0133 + IT_0140;
    const ccomplex_t IT_0142 = IT_0111 + IT_0141;
    const ccomplex_t IT_0143 = conj(IT_0111) + conj(IT_0141);
    const ccomplex_t IT_0144 = s_14*s_34*IT_0000;
    const ccomplex_t IT_0145 = -IT_0144;
    const ccomplex_t IT_0146 = s_13 + IT_0145;
    const ccomplex_t IT_0147 = s_13*IT_0146;
    const ccomplex_t IT_0148 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0015*IT_0016 
      + (-3)*IT_0014*IT_0020)*U_st_00*conj(U_st_01) + IT_0015*IT_0016*U_st_10
      *conj(U_st_11));
    const ccomplex_t IT_0149 = 1.33333333333333*IT_0148;
    const ccomplex_t IT_0150 = (0 + _Complex_I*1.4142135623731)*m_b*conj(U_d2)
      *V_tb*e_em*IT_0020*IT_0027*IT_0029*conj(U_st_00);
    const ccomplex_t IT_0151 = (-0.5)*IT_0150;
    const ccomplex_t IT_0152 = cpow((-2)*s_13 + IT_0011 + IT_0070 + -m_st_1*
      (m_st_1 + (0 + _Complex_I*-1)*Gamma_t1) + reg_prop, -1);
    const ccomplex_t IT_0153 = IT_0151*IT_0152;
    const ccomplex_t IT_0154 = (0 + _Complex_I*1)*IT_0153;
    const ccomplex_t IT_0155 = IT_0149*IT_0154;
    const ccomplex_t IT_0156 = (0 + _Complex_I*1)*e_em*(0.25*(IT_0015*IT_0016 
      + (-3)*IT_0014*IT_0020)*U_st_01*conj(U_st_01) + IT_0015*IT_0016*U_st_11
      *conj(U_st_11));
    const ccomplex_t IT_0157 = 1.33333333333333*IT_0156;
    const ccomplex_t IT_0158 = cpow((-2)*s_13 + IT_0011 + IT_0070 + -m_st_2*
      (m_st_2 + (0 + _Complex_I*-1)*Gamma_t2) + reg_prop, -1);
    const ccomplex_t IT_0159 = IT_0031*IT_0158;
    const ccomplex_t IT_0160 = (0 + _Complex_I*1)*IT_0159;
    const ccomplex_t IT_0161 = IT_0157*IT_0160;
    const ccomplex_t IT_0162 = (-2)*IT_0109 + -IT_0155 + -IT_0161;
    const ccomplex_t IT_0163 = (-2)*IT_0056 + (-2)*IT_0067;
    const ccomplex_t IT_0164 = 6*IT_0163;
    const ccomplex_t IT_0165 = IT_0155 + IT_0161;
    const ccomplex_t IT_0166 = (-6)*IT_0165;
    const ccomplex_t IT_0167 = IT_0164 + IT_0166;
    const ccomplex_t IT_0168 = 6*conj(IT_0163);
    const ccomplex_t IT_0169 = (-6)*conj(IT_0165);
    const ccomplex_t IT_0170 = IT_0168 + IT_0169;
    const ccomplex_t IT_0171 = (0 + _Complex_I*1)*IT_0152;
    const ccomplex_t IT_0172 = V_tb*e_em*V_Wp2*conj(U_st_00);
    const ccomplex_t IT_0173 = IT_0020*IT_0172;
    const ccomplex_t IT_0174 = (0 + _Complex_I*1)*IT_0173;
    const ccomplex_t IT_0175 = -IT_0174;
    const ccomplex_t IT_0176 = IT_0149*IT_0175;
    const ccomplex_t IT_0177 = m_t*V_tb*V_u2*e_em*IT_0027*conj(U_st_10);
    const ccomplex_t IT_0178 = IT_0077*IT_0177;
    const ccomplex_t IT_0179 = 1.4142135623731*IT_0178;
    const ccomplex_t IT_0180 = (0 + _Complex_I*1)*IT_0179;
    const ccomplex_t IT_0181 = (-0.5)*IT_0180;
    const ccomplex_t IT_0182 = IT_0149*IT_0181;
    const ccomplex_t IT_0183 = -IT_0182;
    const ccomplex_t IT_0184 = IT_0176 + IT_0183;
    const ccomplex_t IT_0185 = IT_0171*IT_0184;
    const ccomplex_t IT_0186 = IT_0088*IT_0096;
    const ccomplex_t IT_0187 = IT_0086*IT_0157;
    const ccomplex_t IT_0188 = IT_0082*IT_0157;
    const ccomplex_t IT_0189 = -IT_0188;
    const ccomplex_t IT_0190 = IT_0187 + IT_0189;
    const ccomplex_t IT_0191 = (0 + _Complex_I*1)*IT_0158;
    const ccomplex_t IT_0192 = IT_0190*IT_0191;
    const ccomplex_t IT_0193 = IT_0185 + (-2)*IT_0186 + IT_0192;
    const ccomplex_t IT_0194 = IT_0088*IT_0135;
    const ccomplex_t IT_0195 = IT_0122*IT_0126;
    const ccomplex_t IT_0196 = 2*IT_0194 + 2*IT_0195;
    const ccomplex_t IT_0197 = 6*IT_0196;
    const ccomplex_t IT_0198 = -IT_0185 + -IT_0192;
    const ccomplex_t IT_0199 = (-6)*IT_0198;
    const ccomplex_t IT_0200 = IT_0197 + IT_0199;
    const ccomplex_t IT_0201 = 6*conj(IT_0196);
    const ccomplex_t IT_0202 = (-6)*conj(IT_0198);
    const ccomplex_t IT_0203 = IT_0201 + IT_0202;
    const ccomplex_t IT_0204 = s_13*s_14;
    const ccomplex_t IT_0205 = s_34*IT_0070;
    const ccomplex_t IT_0206 = -IT_0205;
    const ccomplex_t IT_0207 = IT_0204 + IT_0206;
    const ccomplex_t IT_0208 = -IT_0109;
    const ccomplex_t IT_0209 = IT_0056 + IT_0067;
    const ccomplex_t IT_0210 = -(IT_0082 + -IT_0086)*IT_0096;
    const ccomplex_t IT_0211 = -IT_0194 + -IT_0195;
    const ccomplex_t IT_0212 = pow(s_14, 2);
    const ccomplex_t IT_0213 = IT_0000*IT_0212;
    const ccomplex_t IT_0214 = -IT_0070;
    const ccomplex_t IT_0215 = IT_0213 + IT_0214;
    const ccomplex_t IT_0216 = m_b*IT_0215;
    const ccomplex_t IT_0217 = m_b*m_C_2;
    const ccomplex_t IT_0218 = IT_0215*IT_0217;
    const ccomplex_t IT_0219 = 6*IT_0198;
    const ccomplex_t IT_0220 = (-6)*IT_0196;
    const ccomplex_t IT_0221 = IT_0219 + IT_0220;
    const ccomplex_t IT_0222 = 6*conj(IT_0198);
    const ccomplex_t IT_0223 = (-6)*conj(IT_0196);
    const ccomplex_t IT_0224 = IT_0222 + IT_0223;
    const ccomplex_t IT_0225 = pow(s_34, 2);
    const ccomplex_t IT_0226 = IT_0000*IT_0225;
    const ccomplex_t IT_0227 = -IT_0226;
    const ccomplex_t IT_0228 = IT_0011 + IT_0227;
    const ccomplex_t IT_0229 = m_C_2*IT_0228;
    const ccomplex_t IT_0230 = 6*IT_0069;
    const ccomplex_t IT_0231 = 6*conj(IT_0069);
    const ccomplex_t IT_0232 = 6*IT_0193;
    const ccomplex_t IT_0233 = 6*conj(IT_0193);
    const ccomplex_t IT_0234 = 6*IT_0162;
    const ccomplex_t IT_0235 = 6*conj(IT_0162);
    const ccomplex_t IT_0236 = m_b*IT_0146;
    const ccomplex_t IT_0237 = (-18)*conj(IT_0069);
    const ccomplex_t IT_0238 = m_C_2*IT_0146;
    const ccomplex_t IT_0239 = 6*IT_0165;
    const ccomplex_t IT_0240 = (-6)*IT_0163;
    const ccomplex_t IT_0241 = IT_0239 + IT_0240;
    const ccomplex_t IT_0242 = (-6)*IT_0069;
    const ccomplex_t IT_0243 = (-6)*conj(IT_0069);
    const ccomplex_t IT_0244 = (-6)*conj(IT_0163);
    const ccomplex_t IT_0245 = s_13*IT_0215;
    const ccomplex_t IT_0246 = IT_0146*IT_0217;
    const ccomplex_t IT_0247 = (-6)*IT_0193;
    const ccomplex_t IT_0248 = (-6)*conj(IT_0193);
    const ccomplex_t IT_0249 = IT_0217*IT_0228;
    const ccomplex_t IT_0250 = s_14*IT_0001;
    const ccomplex_t IT_0251 = IT_0000*IT_0250;
    const ccomplex_t IT_0252 = 0.5*IT_0251;
    const ccomplex_t IT_0253 = s_14 + IT_0252;
    const ccomplex_t IT_0254 = m_b*IT_0253;
    const ccomplex_t IT_0255 = 2*IT_0254;
    const ccomplex_t IT_0256 = 6*conj(IT_0102);
    const ccomplex_t IT_0257 = IT_0231 + IT_0256;
    const ccomplex_t IT_0258 = s_34*IT_0001;
    const ccomplex_t IT_0259 = IT_0000*IT_0258;
    const ccomplex_t IT_0260 = (-0.25)*IT_0259;
    const ccomplex_t IT_0261 = s_34 + IT_0260;
    const ccomplex_t IT_0262 = m_C_2*IT_0261;
    const ccomplex_t IT_0263 = (-4)*IT_0262;
    const ccomplex_t IT_0264 = 6*conj(IT_0111);
    const ccomplex_t IT_0265 = 6*conj(IT_0141);
    const ccomplex_t IT_0266 = IT_0264 + IT_0265;
    const ccomplex_t IT_0267 = s_14*IT_0011;
    const ccomplex_t IT_0268 = s_13*s_34;
    const ccomplex_t IT_0269 = -IT_0268;
    const ccomplex_t IT_0270 = IT_0267 + IT_0269;
    const ccomplex_t IT_0271 = (-0.25)*IT_0251;
    const ccomplex_t IT_0272 = s_14 + IT_0271;
    const ccomplex_t IT_0273 = m_b*IT_0272;
    const ccomplex_t IT_0274 = (-4)*IT_0273;
    const ccomplex_t IT_0275 = (-6)*conj(IT_0102);
    const ccomplex_t IT_0276 = IT_0243 + IT_0275;
    const ccomplex_t IT_0277 = 0.5*IT_0259;
    const ccomplex_t IT_0278 = s_34 + IT_0277;
    const ccomplex_t IT_0279 = m_C_2*IT_0278;
    const ccomplex_t IT_0280 = 2*IT_0279;
    const ccomplex_t IT_0281 = (-6)*conj(IT_0111);
    const ccomplex_t IT_0282 = (-6)*conj(IT_0141);
    const ccomplex_t IT_0283 = IT_0281 + IT_0282;
    const ccomplex_t IT_0284 = s_13*IT_0228;
    const ccomplex_t IT_0285 = (-6)*IT_0162;
    const ccomplex_t IT_0286 = IT_0001*IT_0217;
    const ccomplex_t IT_0287 = (-3)*IT_0286;
    const ccomplex_t IT_0288 = 3*IT_0286;
    const ccomplex_t IT_0289 = 24*IT_0003;
    const ccomplex_t IT_0290 = pow(m_Z, 4);
    const ccomplex_t IT_0291 = s_13*IT_0290;
    const ccomplex_t IT_0292 = IT_0000*IT_0291;
    const ccomplex_t IT_0293 = (-6)*IT_0292;
    const ccomplex_t IT_0294 = IT_0289 + IT_0293;
    const ccomplex_t IT_0295 = 0.166666666666667*IT_0210;
    const ccomplex_t IT_0296 = 6*IT_0292;
    const ccomplex_t IT_0297 = (-12)*IT_0002;
    const ccomplex_t IT_0298 = IT_0289 + IT_0296 + IT_0297;
    const ccomplex_t IT_0299 = 0.166666666666667*IT_0211;
    const ccomplex_t IT_0300 = 6*IT_0102;
    const ccomplex_t IT_0301 = IT_0230 + IT_0300;
    const ccomplex_t IT_0302 = 0.166666666666667*IT_0255;
    const ccomplex_t IT_0303 = 6*IT_0111;
    const ccomplex_t IT_0304 = 6*IT_0141;
    const ccomplex_t IT_0305 = IT_0303 + IT_0304;
    const ccomplex_t IT_0306 = 0.166666666666667*IT_0263;
    const ccomplex_t IT_0307 = 0.166666666666667*IT_0270;
    const ccomplex_t IT_0308 = (-6)*IT_0102;
    const ccomplex_t IT_0309 = IT_0242 + IT_0308;
    const ccomplex_t IT_0310 = 0.166666666666667*IT_0274;
    const ccomplex_t IT_0311 = (-6)*IT_0111;
    const ccomplex_t IT_0312 = (-6)*IT_0141;
    const ccomplex_t IT_0313 = IT_0311 + IT_0312;
    const ccomplex_t IT_0314 = 0.166666666666667*IT_0280;
    const ccomplex_t IT_0315 = 0.166666666666667*IT_0208;
    const ccomplex_t IT_0316 = 0.166666666666667*IT_0209;
    const ccomplex_t IT_0317 = IT_0285*IT_0307;
    const ccomplex_t IT_0318 = (-6)*conj(IT_0162);
    const ccomplex_t IT_0319 = IT_0270*IT_0318;
    const ccomplex_t IT_0320 = IT_0009*(IT_0103*IT_0104 + IT_0142*IT_0143) +
       IT_0147*(conj(IT_0162)*IT_0167 + IT_0162*IT_0170 + conj(IT_0193)*IT_0200 
      + IT_0193*IT_0203) + IT_0207*(IT_0170*(IT_0208 + IT_0209) + IT_0167*(conj
      (IT_0208) + conj(IT_0209)) + IT_0203*(IT_0210 + IT_0211) + IT_0200*(conj
      (IT_0210) + conj(IT_0211))) + (IT_0143*IT_0167 + IT_0142*IT_0170 + IT_0104
      *IT_0200 + IT_0103*IT_0203)*IT_0216 + IT_0218*(conj(IT_0163)*IT_0200 +
       IT_0163*IT_0203 + conj(IT_0165)*IT_0221 + IT_0165*IT_0224) + IT_0229*
      (conj(IT_0162)*IT_0230 + IT_0162*IT_0231 + IT_0143*IT_0232 + IT_0142
      *IT_0233 + conj(IT_0102)*IT_0234 + IT_0102*IT_0235) + (IT_0104*IT_0232 +
       IT_0103*IT_0233 + IT_0143*IT_0234 + IT_0142*IT_0235)*IT_0236 + IT_0217*
      (IT_0103*((-18)*conj(IT_0111) + (-18)*conj(IT_0141)) + IT_0142*((-18)*conj
      (IT_0102) + IT_0237)) + IT_0238*(IT_0143*IT_0221 + IT_0142*IT_0224 + 6
      *conj(IT_0165)*(IT_0102 + 0.166666666666667*IT_0230) + IT_0165*IT_0231 +
       conj(IT_0102)*IT_0241 + conj(IT_0163)*IT_0242 + IT_0163*IT_0243 + IT_0102
      *IT_0244) + (conj(IT_0163)*IT_0164 + conj(IT_0196)*IT_0197 + conj(IT_0198)
      *IT_0221 + IT_0198*IT_0223 + conj(IT_0165)*IT_0241 + IT_0165*IT_0244)
      *IT_0245 + IT_0246*(conj(IT_0162)*IT_0200 + IT_0162*IT_0203 + conj(IT_0163
      )*IT_0232 + IT_0163*IT_0233 + conj(IT_0165)*IT_0247 + IT_0165*IT_0248) + 
      (conj(IT_0162)*IT_0247 + IT_0162*IT_0248)*IT_0249 + IT_0211*(IT_0255
      *IT_0257 + IT_0263*IT_0266 + IT_0248*IT_0270) + IT_0210*(IT_0248*IT_0270 +
       IT_0274*IT_0276 + IT_0280*IT_0283) + IT_0284*(conj(IT_0193)*IT_0247 +
       conj(IT_0162)*IT_0285) + 6*conj(IT_0211)*(IT_0209*IT_0287 + -IT_0208
      *IT_0288 + IT_0294*IT_0295 + IT_0298*IT_0299 + IT_0301*IT_0302 + IT_0305
      *IT_0306 + IT_0247*IT_0307) + 6*conj(IT_0210)*(IT_0208*IT_0287 + -IT_0209
      *IT_0288 + IT_0295*IT_0298 + IT_0294*IT_0299 + IT_0247*IT_0307 + IT_0309
      *IT_0310 + IT_0313*IT_0314) + 6*conj(IT_0208)*(IT_0210*IT_0287 + -IT_0211
      *IT_0288 + IT_0310*IT_0313 + IT_0309*IT_0314 + IT_0298*IT_0315 + IT_0294
      *IT_0316 + IT_0317) + 6*conj(IT_0209)*(IT_0211*IT_0287 + -IT_0210*IT_0288 
      + IT_0302*IT_0305 + IT_0301*IT_0306 + IT_0294*IT_0315 + IT_0298*IT_0316 +
       IT_0317) + IT_0209*(IT_0257*IT_0263 + IT_0255*IT_0266 + IT_0319) +
       IT_0208*(IT_0276*IT_0280 + IT_0274*IT_0283 + IT_0319);
    return create_ccomplex_return(IT_0320);
}

